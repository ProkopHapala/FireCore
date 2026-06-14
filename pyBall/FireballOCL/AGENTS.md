# FireballOCL: OpenCL Hamiltonian, STM, Grid Projection

## Purpose

OpenCL-accelerated Fireball Hamiltonian assembly, STM current orbital projection, real-space density grid projection, and sparse/dense orbital data utilities.

## Ownership

- OpenCL Hamiltonian kernel wrappers (`OCL_Hamiltonian.py`)
- STM orbital projection and Green's function solvers (`STM.py`, `STM_utils.py`)
- Grid density projection (`Grid.py`)
- Fdata parsing and Fireball reference comparisons (`FdataParser.py`, `Check_Fireball_wrt_Fotran.py`)
- CheFSI spectral solver (`CheFSI.py`, `OMM_ocl.py`)

## Local Contracts

- **Fortran reference is truth:** Any OpenCL implementation must be based on careful analysis of Fortran code, NOT random guessing/flipping/brute-force mapping.
- **Never implement per-point Python loops calling many small GPU kernels:** Production paths must use massively parallel GPU-side evaluation in one kernel or a few large batched kernels. Keep any per-point Python path only as explicit debug/reference/deprecated code.
- **Orbital ordering conventions (critical):**
  - Fortran/Fireball: `[s, py, pz, px]`
  - OpenCL Hamiltonian: `[s, px, py, pz]` (use `_PERM_FORT_TO_HAM = [0,3,1,2]`)
  - OpenCL Grid projection: `[px, py, pz, s]` (use `_PERM_FORT_TO_GRID = [3,1,2,0]`)
- **Sparse block layout:** `H_blocks[iatom, ineigh, :nnu, :nmu]` where `mu` runs on iatom, `nu` on jatom. Dense reconstruction: `blk = H_blocks[i,ineigh,:nj,:ni]; M[i0:i0+ni, j0:j0+nj] += blk.T`.
- **Dense→sparse for Grid.cl:** `blocks[i,ineigh,:nj,:ni] = blk.T` because Grid kernels expect C-order `rho[iatom][ineigh][inu][imu]` (inu major, imu minor).
- **Grid sampling parity:** Fortran `orb2points()` samples at voxel centers; OpenCL grid kernels sample at voxel corners. Shift origin by `0.5*(dA+dB+dC)` for parity.
- **Per-atom orbital counts:** Use `get_orbital_layout(sparse_data, natoms)` from `STM_utils.py`. `n_orb_atom[ia]` from mapping `iatyp[ia]` → `nzx` → `num_orb`. H=1 orb, O/C/N typically 4 (sp).
- **Qneutral_shell:** Per-species, shape `(nsh_max, nspecies_fdata)`. Slice to `nspecies` and zero-fill unused slots.

## Work Guidance

### Hamiltonian Assembly
- `OCL_Hamiltonian.py` — OpenCL kernel wrappers for 2-center and 3-center Hamiltonian integrals
- `Check_Fireball_wrt_Fotran.py` — parity checker against Fortran reference
- `FdataParser.py` — Slater-Koster table parsing

### STM / Transport
- `STM.py` — STM current and orbital projection
- `STM_utils.py` — sparse/dense packing, grid spec building, orbital layout utilities
- `test_CheFSI.py` — spectral solver test

### Grid Projection
- `Grid.py` — real-space density grid projection kernels
- `OMM_ocl.py` / `CheFSI.py` — orbital minimization and spectral solvers
- `FireballPlots.py` — plotting utilities for Fireball data

## Verification

- Run `python3 test_CheFSI.py` or relevant parity script from this directory
- Hamiltonian parity should match Fortran to near-machine-precision
- Grid projection parity: verify against Fortran `orb2points()` reference

## Child DOX Index

- No child subtrees
