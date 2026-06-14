# Shared Resources

## Purpose

Central repository for shared assets used across C++, Python, JavaScript, and OpenCL: example molecules, force field parameters, GPU kernels, GLSL shaders, and precalculated GridFF grids.

## Ownership

- Canonical force field parameter files (loaded by C++, Python, and JS parsers)
- Canonical OpenCL compute kernels (compiled by C++ and PyOpenCL harnesses)
- Example molecules and crystal structures (.mol2, .xyz, .cif, .poscar)
- Precalculated GridFF .npy files for various substrates
- GLSL shaders for WebGL visualization

## Structure

### Molecules & Crystals
- `mol/` — ~40 example molecules in .mol and .mol2 (benzene, water, pentacene, polymers, etc.)
- `xyz/` — ~120 example systems in .xyz (molecules, surfaces, test geometries)
- `crystals/` — CIF, POSCAR, and XYZ files for Si, diamond, CaF2, CaCO3; conversion scripts
- `Substrates/` — Substrate structures (CaF2_6L_Ni3.cif/.xyz) and generated_rect/

### Force Field Parameters (canonical .dat files)
- `ElementTypes.dat` — Element symbols, Z, radii, colors, UFF/QEq params
- `AtomTypes.dat` — Hybridized types (C_sp2, O_hydroxyl, etc.), valence, pi orbitals, MMFF params
- `BondTypes.dat` — Bond length `l0` and stiffness `k` by atom type pair + order
- `AngleTypes.dat` — Equilibrium angle `ang0` and stiffness `k` by atom type triple
- `DihedralTypes.dat` — Torsion barrier `k`, phase `ang0`, periodicity by atom type quadruple + order
- `_FitREQ` / `_PH` variants — Alternative parameter sets for specific fitting methods

**Loaded by**: `MMFFparams.h` (C++), `MMFF.py`/`MMparams` (Python), `MMParams.js` (JS)

### OpenCL Kernels (`cl/`)
Canonical compute kernels compiled by both C++ OpenCL harness and PyOpenCL:
- `relax_multi.cl` — Unified force field kernel (bonds, angles, pi-align, LJ, Morse, H-bond, multi-system)
- `UFF.cl` — UFF force evaluation
- `GridFF.cl` — B-spline grid interpolation and sampling
- `Surface.cl` — Surface interaction (Morse/LJ/Coulomb), FoldedAtomicFunctions, Ewald2D
- `eFF.cl` / `eFF_new.cl` — Electron force field kernels
- `FMM.cl` — Fast Multipole Method
- `Rigid.cl` — Rigid body dynamics
- `LFF.cl` — Local force field

**Note**: Other `cl/` folders in `pyBall/*/cl/` contain domain-specific kernels (DFTB/Grid, PME, hubbard, MQCA, Assembly) and are NOT duplicates of these canonical kernels.

### GLSL Shaders (`shaders/`)
WebGL rendering shaders: atom/bond/label vertex and fragment shaders, MD predictor/corrector, ProjectiveDynamics.

### Precalculated Grids
Large `.npy` B-spline PLQd grids for common substrates:
- `NaCl_1x1_L3/` — Small NaCl(001) reference grid (~120 MB)
- `NaCl_8x8_L3_ClHole/` — Large NaCl(001) with Cl vacancy (~480 MB)
- `CaF2_6L_Ni3_rect_nx2_nz1_L2_top/` — CaF2(111) grid (~400 MB Bspline_PLQd.npy)

**Caution**: These are large binary files. Do not add new grids to git without explicit approval.

## Local Contracts

- **This is the canonical source** for `.dat` parameters and `.cl` kernels. Updates here affect all language bindings.
- **Never duplicate parameter files** into language-specific folders; all parsers load from here.
- **Kernel debug macros** (e.g., `DBG_UFF`) must be enabled by default in the source to guarantee debug prints on the C++ path.
- **Grid files** are runtime artifacts; the generation scripts live in `tests/tMMFF/` and `pyBall/OCL/`.

## Verification

- Parameter parsers: verify all three languages (C++, Python, JS) load the same `.dat` files and produce identical type assignments.
- Kernels: parity tests in `tests/tUFF/`, `tests/tMMFF/`, `tests/tEFF/` compile these kernels.

## Child DOX Index

- No child subtrees
