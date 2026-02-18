

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

## Phase 3: Full Density-Based Model (FDBM) AFM (Fireball DFT)

### What we achieved
- Fixed density projection normalization:
  - `.wf` cutoffs are in **Bohr**; applied Bohr→Å conversion in `Grid.py::load_basis`.
  - Added real spherical-harmonic prefactors (PREF_S/PREF_P) in `Grid.cl`.
  - Verified single-atom integrals: H=1.000 e, C=4.000 e; PTCDA SCF/NA integrals = 140.01 e, Δρ ≈ 0.
- Implemented full-density AFM forces in `tests/tAFM/test_fdbm.py`:
  - Pauli = FFT convolution of sample density with Gaussian tip density (σ=0.7 Å), force = -∇E.
  - Electrostatics = tip density convolved with Hartree potential from Δρ (Poisson FFT), force = -∇E.
  - London C6/r^6 retained; tuned to C6_CO=30 for balanced well.
  - PP relaxation uses these forces; raw and relaxed Fz/df grids with diagnostics are saved.
- Diagnostics: energy-field slices (E_Pauli/E_ES), component traces, raw/relax Fz, df, and FDBM vs Morse comparisons show clear PTCDA ring/bond contrast at h≈3.0–3.4 Å.

### Problems we hit
- Huge density overcount (∫ρ ≈ 6.748×): root cause was treating `.wf` rcutoff as Å instead of Bohr; also missing Ylm prefactors.
- PP forces too repulsive when Pauli dominated; balanced by increasing C6 and using Gaussian tip convolution (σ=0.7 Å) to extend reach.

### Takeaways
- Always confirm units of radial tables (.wf): Bohr cutoffs require explicit Bohr→Å conversion before resampling.
- Apply Ylm prefactors explicitly; radial tables store R(r), not u(r).
- Convolution with a realistic tip density (Gaussian) yields correct height scaling (~3–3.5 Å) and smoother fields than point gradients.
- Keep A_pauli and C6 tunable; small adjustments shift the attractive well significantly.

### Files created/touched
- `pyBall/FireballOCL/cl/Grid.cl` — add PREF_S/PREF_P.
- `pyBall/FireballOCL/Grid.py` — Bohr→Å in load_basis, resampling, normalization.
- `tests/tAFM/test_single_atom.py` — single-atom projection harness.
- `tests/pyFireball/test_grid_projection.py` — RCUT fix.
- `tests/tAFM/test_fdbm.py` — full FDBM force model (Pauli/ES FFT conv., C6 tune), diagnostics/plots.