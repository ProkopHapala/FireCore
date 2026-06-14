# pyFireball Tests: DFTB, STM, Orbital Projection

## Purpose

Fireball/DFTB reference calculations, STM orbital projection, Green's function transport, NEB hydrogen transfer, and GPU density projection parity tests.

## Ownership

- Fireball SCF tests and Fortran reference comparisons
- STM orbital projection and Green's function solvers
- MO coefficient and real-space density exports (`get_wfcoef`, `getGridDens`)
- GPU OpenCL Hamiltonian assembly parity vs Fortran
- NEB H-transfer (molecular and ribbon/k-point)
- Grid density projection and MO vs LDOS validation

## Local Contracts

- **Run from this directory** — scripts use relative paths to data and kernels.
- **Use `run.sh` / `make.sh` scripts** when available; never invoke `make` directly.
- **Fortran reference is truth:** Any OpenCL implementation must be based on careful analysis of Fortran code, NOT random guessing/flipping/brute-force mapping.
- **Fortran debug exports:** Only in `fortran/MODULES/debug.f90`, gated by `idebugWrite`, marked with `! DEBUG : TO EXPORT` blocks.
- **Orbital ordering conventions:**
  - Fortran/Fireball: `[s, py, pz, px]`
  - OpenCL Hamiltonian: `[s, px, py, pz]` (use `_PERM_FORT_TO_HAM = [0,3,1,2]`)
  - OpenCL Grid projection: `[px, py, pz, s]` (use `_PERM_FORT_TO_GRID = [3,1,2,0]`)
- **Sparse block layout:** `H_blocks[iatom, ineigh, :nnu, :nmu]` where `mu` runs on iatom and `nu` on jatom. Reconstruct dense: `blk = H_blocks[i,ineigh,:nj,:ni]; M[i0:i0+ni, j0:j0+nj] += blk.T`.
- **Grid sampling parity:** Fortran `orb2points()` samples at voxel centers; OpenCL grid kernels sample at voxel corners. Shift origin by `0.5*(dA+dB+dC)` for parity.

## Work Guidance

### Fireball Reference & Parity
- `test_h2o_fortran_only.py` — pure Fortran reference for H2O
- `test_h2o_opencl_vs_fortran.py` — OpenCL Hamiltonian assembly parity
- `test_h2o_opencl_debug.py` — debug instrumentation for OpenCL
- `export_HS_sparse.py` — export sparse Hamiltonian/overlap blocks

### STM / Orbital Projection
- `test_stm_orbital_projection.py` — orbital projection STM
- `test_stm_orbital_rotated.py` / `test_stm_orbital_rotated_2mol.py` — rotated orbital projection
- `test_stm_gf_dyson_2mol.py` / `test_stm_gf_dyson_2mol_ocl.py` — Green's function Dyson solver (CPU vs OpenCL)
- `test_mo_vs_ldos.py` — MO coefficient vs LDOS comparison
- `test_grid_projection.py` — GPU grid density projection

### Density / Response
- `test_dens2points.py` — density to grid points
- `test_response_function.py` / `test_response_function_rotated.py` — response function tests
- `run_density_on_samples.py` / `sample_density_points.py` — density sampling

### NEB / H-transfer
- `neb_h_transfer_molecules.py` — NEB for molecular H-transfer (gamma-point)
- `neb_h_transfer.py` — NEB for ribbon H-transfer (k-point, PBC)
- `build_two_ribbons.py` — ribbon unit cell builder
- `scan_constrained.py` / `scan_ribbon.py` / `scan_LHb.py` — constrained scans

## Verification

- Run relevant test script from this directory
- OpenCL parity should match Fortran to near-machine-precision (~1e-5 or better)
- STM orbital projection: verify against Fireball reference data
- NEB: check energy profiles are smooth and barriers are physical

## Child DOX Index

- No child subtrees
