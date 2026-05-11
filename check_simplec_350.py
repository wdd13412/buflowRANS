#!/usr/bin/env python3
"""Run buflow_run until SIMPLEC iteration 350 and verify bounded car targets.

The car mesh is large enough that running all 3000 SIMPLEC iterations can be slow.
This helper streams solver output, stops as soon as the requested iteration line is
seen, and can optionally check reported max speed and pressure range against
user-provided bounds. By default it only checks that the solver reaches the
target iteration without NaN/ERROR output.
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", default="./buflow_run", help="solver executable to run")
    parser.add_argument("--target-iter", type=int, default=350, help="iteration to verify")
    parser.add_argument("--timeout", type=float, default=420.0, help="wall-clock timeout in seconds")
    parser.add_argument("--expected-umax", type=float, default=None, help="optional expected Umax at target")
    parser.add_argument("--expected-dp", type=float, default=None, help="optional expected pressure range at target")
    parser.add_argument("--umax-tol", type=float, default=0.05, help="absolute Umax tolerance")
    parser.add_argument("--dp-tol", type=float, default=0.05, help="absolute pressure-range tolerance")
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

    assert proc.stdout is not None
    try:
        while True:
            if time.monotonic() - start > args.timeout:
                print(
                    f"ERROR: timeout after {args.timeout:.1f}s before SIMPLEC {args.target_iter}; "
                    f"last iter={last_iter}, Umax={last_umax}, dp={last_dp}",
                    file=sys.stderr,
                )
                stop_process(proc)
                return 124

            line = proc.stdout.readline()
            if line == "":
                rc = proc.poll()
                if rc is not None:
                    print(
                        f"ERROR: solver exited with code {rc} before SIMPLEC {args.target_iter}; "
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

            match = SIMPLEC_RE.search(line)
            if match is None:
                continue

            last_iter = int(match.group("iter"))
            try:
                last_umax = float(match.group("umax"))
                last_dp = float(match.group("dp"))
            except ValueError:
                print(f"ERROR: non-finite/overflow SIMPLEC line: {line.rstrip()}", file=sys.stderr)
                stop_process(proc)
                return 1

            if last_iter >= args.target_iter:
                stop_process(proc)
                checks_ok = True
                messages = []
                if args.expected_umax is not None:
                    umax_ok = abs(last_umax - args.expected_umax) <= args.umax_tol
                    checks_ok = checks_ok and umax_ok
                    messages.append(f"Umax={last_umax:.3f} expected {args.expected_umax}±{args.umax_tol}")
                else:
                    messages.append(f"Umax={last_umax:.3f}")
                if args.expected_dp is not None:
                    dp_ok = abs(last_dp - args.expected_dp) <= args.dp_tol
                    checks_ok = checks_ok and dp_ok
                    messages.append(f"dp={last_dp:.3f} expected {args.expected_dp}±{args.dp_tol}")
                else:
                    messages.append(f"dp={last_dp:.3f}")
                if checks_ok:
                    print(f"PASS: SIMPLEC {last_iter} reached target iteration: " + ", ".join(messages))
                    return 0
                print(
                    f"ERROR: SIMPLEC {last_iter} outside requested bounds: " + ", ".join(messages),
                    file=sys.stderr,
                )
                return 1
    except KeyboardInterrupt:
        try:
            os.killpg(proc.pid, signal.SIGTERM)
        finally:
            stop_process(proc)
        raise


if __name__ == "__main__":
    raise SystemExit(main())
