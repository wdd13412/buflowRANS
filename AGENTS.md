# AGENTS.md

## Cursor Cloud specific instructions

### Overview

buflowRANS is a RANS k-omega CFD solver written in Fortran 90. There are two source files:
- `buflowRANS.f90` — main module (~4300 lines)
- `run_parameter.f90` — entry point

There is no build system (no Makefile/CMake/fpm). No automated test suite exists.

### System dependency

`gfortran` must be installed (`apt-get install -y gfortran`). The update script handles this.

### Build

```bash
gfortran -O2 -o buflowRANS buflowRANS.f90 run_parameter.f90
```

### Run

The solver expects an OpenFOAM polyMesh directory at `mesh/OFairfoilMesh` (hardcoded path relative to CWD). That directory must contain `points`, `faces`, `owner`, `neighbour`, and `boundary` files. The repository does **not** ship mesh data. Use `python3 generate_mesh.py` (committed in this repo) to generate a small test mesh.

**Boundary file caveat**: The Fortran parser cannot handle semicolons (`;`) in the boundary file. The `generate_mesh.py` script produces files without semicolons. If you use a different mesh source, strip semicolons from the `boundary` file.

Once mesh data exists:
```bash
./buflowRANS
```

The solver writes:
- `FvCFDRestart.txt` — restart state
- `solution_RANS.*.vtk` — VTK visualization output

### Lint / Static checks

No linter is configured. Compilation itself (`gfortran`) is the primary correctness check. Add `-Wall -Wextra` for warnings:
```bash
gfortran -Wall -Wextra -O2 -o buflowRANS buflowRANS.f90 run_parameter.f90
```

### Tests

No automated test suite exists. Verification is done by compiling and running the solver with valid mesh data and checking for convergence/valid output.
