

Below is a concise report (pre-FDBM), based on [doc/Topics/AFM/AFM_migration_discusion.md](cci:7://file:///home/prokophapala/git/FireCore/doc/Topics/AFM/AFM_migration_discusion.md:0:0-0:0) and the work done so far.

## Goal
Build a simple, reusable PyOpenCL AFM tool:
1) Morse/LJ + fixed charges AFM imaging (probe-particle relaxation, df maps) on planar molecules (e.g., PTCDA).
2) Add QEq-derived charges.
3) Later: Full Density (Fireball DFT) AFM via FFT-based Pauli/Coulomb convolutions.

## What we implemented (pre-FDBM)
- **AFMulator class** ([pyBall/OCL/AFM.py](cci:7://file:///home/prokophapala/git/FireCore/pyBall/OCL/AFM.py:0:0-0:0)): wraps existing kernels, context, and FFT from `GridFF.py`; loads `relax.cl` for force-field generation and probe relaxation.
- **Test harness** (`tests/tAFM/test_ptcda.py`, `tests/tAFM/run.sh`): CLI to run LJ/Morse, optional QEq, plots and npy outputs.
- **Kernels**: Added `evalMorseC_QZs_toImg` to `cpp/common_resources/cl/relax.cl` (Morse variant of LJ kernel).
- **Validation (GTX 1650)**:
  - LJ + fixed Q: Fz ~0–2 eV/Å; PTCDA ring pattern visible.
  - Morse + fixed Q: Correct attractive minimum then repulsion (sign bug in `getMorse` fixed).
  - LJ + QEq: Charge sum preserved; large charges noted for PTCDA parameters.
- **Coordinate convention**: Molecule shifted into [0,L] kernel space; normalized sampler uses `dinvA=(1/Lx,0,0)`, per-voxel step `dA=(Lx/nx,0,0)` in `evalLJC_QZs_toImg`.

## Key problems and fixes
- **Morse force sign was wrong** (`relax.cl/getMorse`): caused global repulsion and huge forces; fixed sign to restore attraction → repulsion.
- **Duplication risk**: Multiple host/kernels (GridFF Bspline vs relax trilinear). Strategy: keep `relax.cl` as canonical AFM/PP kernel set; keep `GridFF.cl` only for Bspline/FFT paths; avoid forking kernels.
- **Charge handling**: QEq charges can be large; keep ElementTypes.dat as single source of truth; QEq solver to be refined/ported carefully.
- **Interpolation**: For “fast” AFM, use trilinear sampling in `relax.cl` (avoid Bspline fitting overhead in `GridFF.cl`).

## Remaining work
- **FDBM/Fireball (Phase 3)**: Not yet integrated in [AFMulator](cci:2://file:///home/prokophapala/git/FireCore/pyBall/OCL/AFM.py:6:0-344:38). Planned path:
  - Project density via [pyBall/FireballOCL/Grid.py](cci:7://file:///home/prokophapala/git/FireCore/pyBall/FireballOCL/Grid.py:0:0-0:0) + [cl/Grid.cl](cci:7://file:///home/prokophapala/git/FireCore/pyBall/FireballOCL/cl/Grid.cl:0:0-0:0) (or C++ `OCL_DFT.h`).
  - Use FFT Poisson/convolution from `GridFF.cl` (`poissonW`, `convolution`) to build Pauli/Coulomb grids.
  - Feed summed force grid into `relaxStrokesTilted` → df maps.
- **QEq robustness**: Revisit parameterization and solver stability; consider GPU port after CPU correctness.
- **Height alignment**: When comparing to Morse/LJ, be mindful of shifts due to density tails.

## Insights / takeaways for future debugging
- **Kernel reuse over forking**: Keep `relax.cl` as AFM canonical; `GridFF.cl` only for FFT/Bspline/Poisson. Add variants via flags/macros, not copies.
- **Sign and sampler checks**: Small sign errors (as in `getMorse`) dominate forces; verify sampler normalization (`dinvA`) whenever grids are shifted.
- **Single source for parameters**: ElementTypes.dat drives LJ/Morse/Q/QEq in both C++ and Python; avoid duplicating constants.
- **Trilinear is enough for AFM imaging**: Bspline path adds fitting overhead; use it only when high interpolation accuracy is needed.
- **QEq can explode**: Validate charge magnitudes early (sum=0, ranges reasonable) before coupling to force-field kernels.