#!/usr/bin/env python3
"""Run buflow_run until SIMPLEC iteration 350 and verify physical ranges.

The car mesh is large enough that running all 3000 SIMPLEC iterations can be slow.
This helper streams solver output, waits for the requested SIMPLEC iteration and
its diagnostic block, and checks that the reported max speed, pressure range, and
boundary mass-flux balance remain in physically plausible low-Mach ranges.
"""

from __future__ import annotations

import argparse
import os
import re
import signal
import subprocess
import sys
import time
from pathlib import Path

SIMPLEC_RE = re.compile(
    r"SIMPLEC\s+(?P<iter>\d+)\s+res\(u,v,p\)=.*?"
    r"Umax=\s*(?P<umax>[*+\-0-9.Ee]+)\s+dp=\s*(?P<dp>[*+\-0-9.Ee]+)"
)
DIAG_RE = re.compile(r"diag\s+(?P<iter>\d+)\b.*rawDp=\s*(?P<raw_dp>[*+\-0-9.Ee]+)")
BC_RE = re.compile(r"\bbc\s+(?P<index>\d+)\s+.*?flux=\s*(?P<flux>[*+\-0-9.Ee]+)")
PHYS_RE = re.compile(
    r"phys\s+(?P<iter>\d+)\s+vMax=\s*(?P<vmax>[*+\-0-9.Ee]+)\s+"
    r"vRms=\s*(?P<vrms>[*+\-0-9.Ee]+)\s+frontP=\s*(?P<frontp>[*+\-0-9.Ee]+)\s+"
    r"roofP=\s*(?P<roofp>[*+\-0-9.Ee]+)\s+rearP=\s*(?P<rearp>[*+\-0-9.Ee]+)\s+"
    r"wakeUx=\s*(?P<wakeux>[*+\-0-9.Ee]+)\s+wakeDef=\s*(?P<wakedef>[*+\-0-9.Ee]+)"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", default="./buflow_run", help="solver executable to run")
    parser.add_argument("--target-iter", type=int, default=350, help="iteration to verify")
    parser.add_argument("--timeout", type=float, default=420.0, help="wall-clock timeout in seconds")
    parser.add_argument("--expected-umax", type=float, default=None, help="optional exact expected Umax at target")
    parser.add_argument("--expected-dp", type=float, default=None, help="optional exact expected pressure range at target")
    parser.add_argument("--umax-tol", type=float, default=0.05, help="absolute Umax tolerance for --expected-umax")
    parser.add_argument("--dp-tol", type=float, default=0.05, help="absolute pressure-range tolerance for --expected-dp")
    parser.add_argument("--min-umax", type=float, default=8.0, help="minimum physically plausible Umax at target")
    parser.add_argument("--max-umax", type=float, default=15.0, help="maximum physically plausible Umax at target")
    parser.add_argument("--min-dp", type=float, default=0.0, help="minimum physically plausible pressure range at target")
    parser.add_argument("--max-dp", type=float, default=150.0, help="maximum physically plausible pressure range at target")
    parser.add_argument(
        "--max-mass-imbalance",
        type=float,
        default=1.0e-6,
        help="maximum absolute sum of target diagnostic boundary fluxes",
    )
    parser.add_argument(
        "--expected-boundaries",
        type=int,
        default=5,
        help="number of boundary diagnostic flux lines expected after target iteration",
    )
    parser.add_argument(
        "--skip-diagnostics-check",
        action="store_true",
        help="stop at the target SIMPLEC line without waiting for the diagnostic block",
    )
    parser.add_argument(
        "--min-vmax",
        type=float,
        default=5.0e-2,
        help="minimum target cross-stream velocity magnitude; catches frozen 1-D flow fields",
    )
    parser.add_argument(
        "--max-wake-accel",
        type=float,
        default=0.5,
        help="maximum allowed negative wake deficit, i.e. wake acceleration over upstream Ux",
    )
    return parser.parse_args()


def stop_process(proc: subprocess.Popen[str]) -> None:
    if proc.poll() is not None:
        return
    try:
        proc.terminate()
        proc.wait(timeout=5.0)
    except subprocess.TimeoutExpired:
        proc.kill()
        proc.wait(timeout=5.0)


def finite_float(text: str, label: str, line: str) -> float:
    if "*" in text:
        raise ValueError(f"overflow in {label}: {line.rstrip()}")
    value = float(text)
    if value != value or value in (float("inf"), float("-inf")):
        raise ValueError(f"non-finite {label}: {line.rstrip()}")
    return value


def validate_ranges(args: argparse.Namespace, last_iter: int, last_umax: float, last_dp: float) -> tuple[bool, list[str]]:
    checks_ok = True
    messages: list[str] = []

    if args.expected_umax is not None:
        umax_ok = abs(last_umax - args.expected_umax) <= args.umax_tol
        checks_ok = checks_ok and umax_ok
        messages.append(f"Umax={last_umax:.3f} expected {args.expected_umax}±{args.umax_tol}")
    else:
        umax_ok = args.min_umax <= last_umax <= args.max_umax
        checks_ok = checks_ok and umax_ok
        messages.append(f"Umax={last_umax:.3f} in [{args.min_umax}, {args.max_umax}]")

    if args.expected_dp is not None:
        dp_ok = abs(last_dp - args.expected_dp) <= args.dp_tol
        checks_ok = checks_ok and dp_ok
        messages.append(f"dp={last_dp:.3f} expected {args.expected_dp}±{args.dp_tol}")
    else:
        dp_ok = args.min_dp <= last_dp <= args.max_dp
        checks_ok = checks_ok and dp_ok
        messages.append(f"dp={last_dp:.3f} in [{args.min_dp}, {args.max_dp}]")

    messages.insert(0, f"SIMPLEC {last_iter}")
    return checks_ok, messages


def main() -> int:
    args = parse_args()
    exe = Path(args.exe)
    if not exe.exists():
        print(f"ERROR: executable not found: {exe}", file=sys.stderr)
        return 2
    exe_cmd = str(exe.resolve())

    start = time.monotonic()
    proc = subprocess.Popen(
        [exe_cmd],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
        preexec_fn=os.setsid,
    )

    last_iter: int | None = None
    last_umax: float | None = None
    last_dp: float | None = None
    target_reached = False
    target_diag_seen = False
    mass_check_done = False
    mass_imbalance: float | None = None
    target_fluxes: list[float] = []

    assert proc.stdout is not None
    try:
        while True:
            if time.monotonic() - start > args.timeout:
                print(
                    f"ERROR: timeout after {args.timeout:.1f}s before completing SIMPLEC {args.target_iter}; "
                    f"last iter={last_iter}, Umax={last_umax}, dp={last_dp}, "
                    f"target_diag_seen={target_diag_seen}, flux_lines={len(target_fluxes)}, "
                    f"mass_check_done={mass_check_done}",
                    file=sys.stderr,
                )
                stop_process(proc)
                return 124

            line = proc.stdout.readline()
            if line == "":
                rc = proc.poll()
                if rc is not None:
                    print(
                        f"ERROR: solver exited with code {rc} before completing SIMPLEC {args.target_iter}; "
                        f"last iter={last_iter}, Umax={last_umax}, dp={last_dp}",
                        file=sys.stderr,
                    )
                    return rc if rc != 0 else 1
                time.sleep(0.05)
                continue

            print(line, end="")
            if "ERROR:" in line or "NaN" in line:
                print("ERROR: solver reported failure/NaN", file=sys.stderr)
                stop_process(proc)
                return 1
            if target_reached and "diag warning:" in line:
                print(f"ERROR: target diagnostic warning: {line.rstrip()}", file=sys.stderr)
                stop_process(proc)
                return 1

            simplec_match = SIMPLEC_RE.search(line)
            if simplec_match is not None:
                last_iter = int(simplec_match.group("iter"))
                try:
                    last_umax = finite_float(simplec_match.group("umax"), "Umax", line)
                    last_dp = finite_float(simplec_match.group("dp"), "dp", line)
                except ValueError as exc:
                    print(f"ERROR: {exc}", file=sys.stderr)
                    stop_process(proc)
                    return 1

                if last_iter >= args.target_iter:
                    target_reached = True
                    checks_ok, messages = validate_ranges(args, last_iter, last_umax, last_dp)
                    if not checks_ok:
                        stop_process(proc)
                        print("ERROR: physical range check failed: " + ", ".join(messages), file=sys.stderr)
                        return 1
                    if args.skip_diagnostics_check:
                        stop_process(proc)
                        print("PASS: " + ", ".join(messages))
                        return 0
                continue

            if not target_reached:
                continue

            diag_match = DIAG_RE.search(line)
            if diag_match is not None and int(diag_match.group("iter")) >= args.target_iter:
                target_diag_seen = True
                target_fluxes = []
                continue

            bc_match = BC_RE.search(line)
            if target_diag_seen and bc_match is not None:
                try:
                    target_fluxes.append(finite_float(bc_match.group("flux"), "boundary flux", line))
                except ValueError as exc:
                    print(f"ERROR: {exc}", file=sys.stderr)
                    stop_process(proc)
                    return 1

                if len(target_fluxes) >= args.expected_boundaries and not mass_check_done:
                    mass_imbalance = abs(sum(target_fluxes))
                    mass_check_done = True
                    if mass_imbalance > args.max_mass_imbalance:
                        stop_process(proc)
                        _, messages = validate_ranges(args, last_iter or args.target_iter, last_umax or 0.0, last_dp or 0.0)
                        messages.append(f"mass imbalance={mass_imbalance:.3e} <= {args.max_mass_imbalance:.3e}")
                        print("ERROR: mass balance check failed: " + ", ".join(messages), file=sys.stderr)
                        return 1
                    continue

            phys_match = PHYS_RE.search(line)
            if target_reached and phys_match is not None and int(phys_match.group("iter")) >= args.target_iter:
                try:
                    vmax = finite_float(phys_match.group("vmax"), "vMax", line)
                    frontp = finite_float(phys_match.group("frontp"), "frontP", line)
                    roofp = finite_float(phys_match.group("roofp"), "roofP", line)
                    rearp = finite_float(phys_match.group("rearp"), "rearP", line)
                    wakedef = finite_float(phys_match.group("wakedef"), "wakeDef", line)
                except ValueError as exc:
                    print(f"ERROR: {exc}", file=sys.stderr)
                    stop_process(proc)
                    return 1

                _, messages = validate_ranges(args, last_iter or args.target_iter, last_umax or 0.0, last_dp or 0.0)
                messages.append(f"mass imbalance={mass_imbalance:.3e} <= {args.max_mass_imbalance:.3e}")
                messages.append(f"vMax={vmax:.3e} >= {args.min_vmax:.3e}")
                messages.append(f"frontP={frontp:.3e} > roofP={roofp:.3e} > rearP={rearp:.3e}")
                messages.append(f"wakeDef={wakedef:.3e} >= -{args.max_wake_accel:.3e}")

                physics_ok = mass_check_done and mass_imbalance is not None
                physics_ok = physics_ok and vmax >= args.min_vmax
                physics_ok = physics_ok and frontp > roofp and roofp > rearp
                physics_ok = physics_ok and wakedef >= -args.max_wake_accel
                stop_process(proc)
                if physics_ok:
                    print("PASS: " + ", ".join(messages))
                    return 0
                print("ERROR: contour-physics check failed: " + ", ".join(messages), file=sys.stderr)
                return 1
    except KeyboardInterrupt:
        try:
            os.killpg(proc.pid, signal.SIGTERM)
        finally:
            stop_process(proc)
        raise


if __name__ == "__main__":
    raise SystemExit(main())
