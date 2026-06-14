# DFTB+ Tests and Integration

## Purpose

DFTB+ standalone program, C API, and Python wrapper tests. Library interfaces, parsers, eigenvector export, and OpenCL grid projection integration.

## Ownership

- DFTB+ Python library interface (`dftb.py`)
- C-API tests and Python wrapper validation
- Waveplot orbital visualization and cube file I/O
- Hessian/vibration calculations via ASE
- 3D grid density projection and dense projection tests
- DFTB+ eigenvector export for OpenCL orbital projection

## Local Contracts

- **Run from this directory** — scripts use relative paths to data and kernels.
- **Use `run.sh` / `make.sh` scripts** when available; never invoke `make` directly.
- **DFTB+ executable:** External dependency; ensure `dftb+` is in PATH or specified explicitly.
- **Slater-Koster files:** Required (`3ob-3-1` or similar); paths must be set in HSD input.

## Work Guidance

### Library & API
- `example_dftb_lib.py` — basic DFTB+ library usage
- `test_python_api.py` — Python API validation
- `dftb.py` — main DFTB+ interface module

### Waveplot & Orbitals
- `example_orbitals.py` — orbital visualization
- `test_waveplot_dftb.py` / `test_waveplot_dftbcore.py` — waveplot integration
- `compare_waveplot_lib.py` — compare waveplot library outputs

### Density & Projection
- `test_3d_grid_density.py` — 3D grid density evaluation
- `test_dense_projection.py` — dense matrix projection
- `compare_density_multizeta.py` — multi-zeta density comparison

### Scans & Jobs
- `dftb_scan.py` / `dftb_scan_2.py` / `dftb_scan_getE.py` — potential energy scans
- `dftb_jobs_frags.py` — fragment-based job generation
- `dftb_post_proc.py` / `dftb_postproc.py` — post-processing

## Verification

- Run `python3 test_python_api.py` or `python3 example_dftb_lib.py` from this directory
- Waveplot tests require DFTB+ with waveplot support compiled in
- Grid density tests should produce smooth, physical densities

## Child DOX Index

- No child subtrees
