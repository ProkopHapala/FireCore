# Indranil Branch Summary

## Overview
This branch contains all essential functionality from `debug/grid_geration_and_test` merged into a clean branch based on `master`. Debug prints and profiler files have been removed while preserving all functional changes.

## Commits (9 total, from master)

| Commit | Description |
|--------|-------------|
| `0992f3a0` | Add scan test scripts and utilities |
| `68a8dc6a` | Add test scripts and utilities |
| `9a666cc5` | Add remaining C++ files |
| `c3c0cd77` | Add supporting C++ files |
| `6b168381` | Add C++ core files |
| `6e39b668` | Add complete tGridFF folder |
| `52a2a04c` | Clean Python core files |
| `79567ded` | Clean ocl_GridFF_new.py |
| `3b0e7395` | Clean GridFF.py |

## Files Changed Summary

### C++ Core Files (from debug branch)
- `cpp/common/molecular/GridFF.h` - Grid force field with Morse/Coulomb
- `cpp/common/molecular/MolWorld_sp3.h` - Molecular world with scan_rigid_uff, scan_constr
- `cpp/common/molecular/UFF.h` - UFF implementation (tested vs LAMMPS)
- `cpp/common/molecular/MMFFparams.h` - Parameter handling
- `cpp/libs/Molecular/MMFF_lib.cpp` - Library interface
- `cpp/common_resources/cl/GridFF.cl` - OpenCL kernels (poissonW, slabPotential_zyx)

### Python Core Files (cleaned)
- `pyBall/OCL/GridFF.py` - OpenCL GridFF (debug prints removed)
- `pyBall/MMFF.py` - MMFF interface (debug prints removed)
- `pyBall/atomicUtils.py` - getVdWparams uses bDict=True, loadElementTypes key fix
- `pyBall/tests/ocl_GridFF_new.py` - OCL test (debug prints removed)
- `pyBall/tests/utils.py` - plot1Dcut with pos0 parameter

### Test Scripts (from debug branch)
- `tests/tGridFF/` - Complete folder with generate_grid.py, generate_scans.py, all_scan.py
- `tests/tMMFF/` - run.py, batch_compare.py, generate_scans.py, scan_*.py
- `tests/tUFF/run.py` - UFF vs LAMMPS test
- `tests/tPTCDA/`, `tests/tPTCDA2/` - PTCDA tests

### Documentation
- `doc/py/` - Visualization and plotting utilities

## Key Functional Changes from Debug Branch

1. **GPU Ewald Calculation** - Grid shifted by 1 pixel for CPU/LAMMPS agreement (commit cdf029af)
2. **UFF Implementation** - Tested vs LAMMPS, scan_rigid_uff function (commits 3657e5e1, 31f4d5d9)
3. **Grid Origin Fix** - Consistent z-coordinate calculation (commit caf4870b)
4. **VdW Parameter Lookup** - Uses atomic number instead of element name (commit 4f4b50f6)
5. **Poisson Kernel** - Corrected k-kernel normalization (commit 76ca776a)

## Excluded Files (per MERGE_STRATEGY.md)
- All profiler files (`*profil*`, `debug_*`)
- OpenCL cache files (`source.cl` in cache directories)
- Temporary plot files in cpp/common_resources subdirectories

## Differences from Debug Branch
The only differences are intentionally removed debug prints (~88 print statements in GridFF.py, ~20 in MMFF.py, etc.) while preserving all functional code.

## Next Steps
1. Test the branch to ensure functionality matches debug branch
2. Merge to master when ready
