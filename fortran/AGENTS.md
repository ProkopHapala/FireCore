# Fortran Fireball/DFTB Solver

## Purpose

Fortran reference implementation of the Fireball DFTB solver: Hamiltonian assembly, SCF, forces, and sparse matrix exports. This is the **source of truth** for all GPU/OpenCL parity work.

## Ownership

- `MAIN/` — Main programs and `libFireCore.f90` C/Fortran bindings
- `ASSEMBLERS/` — 2-center and 3-center Hamiltonian/overlap integrals
- `DASSEMBLERS/` — Force derivatives (dH/dR, dS/dR)
- `GRID/` — Grid-based orbital evaluation and density projection
- `INITIALIZERS/` — System initialization and parameter loading
- `INTERACTIONS/` — Non-bonded and exchange-correlation interactions
- `INTERPOLATERS/` — Spline interpolation for radial integrals
- `MODULES/` — Global modules, data structures, and debug utilities
- `NEIGHBORS/` — Neighbor list construction
- `READFILES/` — Input file parsers (XYZ, Fdata, etc.)
- `ROTATIONS/` — Angular momentum coupling and rotation matrices
- `MATH/` — Mathematical utilities
- `doc/` — Fortran-specific documentation

## Local Contracts

- **Fortran reference is truth:** Any OpenCL implementation must be based on careful analysis of Fortran code, NOT random guessing/flipping/brute-force mapping.
- **Debug exports only in `fortran/MODULES/debug.f90`:** Gated by `idebugWrite`, marked with `! DEBUG : TO EXPORT` blocks.
- **No permanent allocations in `allocate_h.f90`/`reallocate_h.f90`:** Use temporary exports only.
- **Minimal modifications:** Only tiny, fully commented, verbosity-gated changes allowed in production Fortran files.
- **Column-major, 1-based indexing:** Fortran arrays are column-major with 1-based indices. When reshaping in Python/C, reverse axes then transpose.

## Work Guidance

### Key Entry Points
- `libFireCore.f90` — C bindings: `firecore_get_HS_sparse()`, `firecore_get_wfcoef()`, `firecore_getGridDens()`
- `MAIN/` — Standalone Fireball programs

### Hamiltonian Assembly
- `ASSEMBLERS/` — `assemble_ca_2c.f90`, `assemble_ca_3c.f90` and variants
- Use `assemble_ca_2c_new.f90` as reference for debug instrumentation style

### Sparse Matrix Exports
- `h_mat_out(mu, nu, ineigh, iatom)` where `mu` runs orbitals on iatom, `nu` on jatom
- In Python: `H_blocks[iatom, ineigh, :nnu, :nmu]`; reconstruct dense via `blk.T`

### Grid Density
- `GRID/` — `orb2points()` samples at voxel centers (not corners)
- OpenCL grid kernels sample at corners; shift origin by `0.5*(dA+dB+dC)` for parity

## Verification

- Build via `make.sh` or project build scripts from repo root; never invoke `make` directly in this directory
- Run Fortran-only tests (e.g., `tests/pyFireball/test_h2o_fortran_only.py`) to verify before comparing OpenCL

## Child DOX Index

- `MODULES/` — Global data and debug export utilities
- `ASSEMBLERS/` — Hamiltonian integral assembly
- `GRID/` — Grid-based orbital/density evaluation
