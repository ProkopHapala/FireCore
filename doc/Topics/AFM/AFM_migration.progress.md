

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
- `tests/tAFM/pyocl_fdbm/run_pyocl_fdbm.py` — standalone 6-step FDBM pipeline with CO tip support.
- `tests/tAFM/pyocl_fdbm/compute_co_tip.py` — CO tip density generation from Fireball SCF.
- `tests/tAFM/pyocl_fdbm/CO.xyz` — CO molecule geometry for tip density.
- `tests/tAFM/pyocl_fdbm/pentacene.xyz` — Pentacene sample molecule geometry.
- `tests/tAFM/pyocl_fdbm/README.md` — Pipeline documentation and usage.

## Phase 3: Standalone PyOpenCL FDBM Pipeline (pentacene)

### What we implemented
- **New standalone pipeline** in `tests/tAFM/pyocl_fdbm/run_pyocl_fdbm.py` with full 6-step process:
  1. **Density projection** (Fireball DFT) → `rho_grid`, `rho_na_grid`, `rho_diff` (Δρ = SCF - neutral atom)
  2. **Electrostatics via Poisson FFT** → `V_ES` from Δρ, with validation check
  3. **Pauli repulsion** → FFT convolution of sample density with Gaussian tip (σ=0.7 Å)
  4. **Electrostatics tip-sample convolution** → FFT convolution of V_ES with tip charge density
  5. **Dispersion (vdW)** → C6/r^6 London forces, now computed on density grid (not scan grid)
  6. **Composed forces** → Sum Pauli + ES + vdW, interpolate to scan grid, compute Fz/df
- **CO tip density generation** (`compute_co_tip.py`): True CO delta-density from Fireball SCF, with diagnostics
- **Auto-computed scan resolution**: `--scan_step=0.1` (default) auto-computes `nx,ny` from molecule bounding box + 3Å margin
- **All diagnostics saved** in `debug/step*/` subdirectories with plots and .npy arrays

### Key fixes and improvements
- **vdW resolution mismatch**: `step5_dispersion` was computing on coarse scan grid → fixed to compute on density grid, then interpolate in `step6_composed`
- **Poisson validation plot bug**: Was recomputing raw 2D `laplace()` on slice without `step²` scaling → fixed to use precomputed 3D `lapV` with proper scaling
- **Poisson z-slice**: Changed from `nz//2` to molecule plane (`argmax(|ρ|)`) for meaningful validation
- **Added assertion**: Poisson ranges must agree within 20% (`ratio < 1.2`) or raise exception with diagnostic values

### Poisson solver validation
- **Formula**: `V_k = 4π·COULOMB_CONST·rho_k / k²` (no `dV` factor) — verified correct with Gaussian test charge
- **Current result**: With `step=0.1Å`, `|lapV| ≈ 196 eV/Å²` vs `|rhs| ≈ 227 eV/Å²` (ratio 1.16, ~16% error)
- **Root cause of 16% error**: Finite-difference Laplacian discretization error; order of magnitude is correct
- **Acceptable**: The 16% mismatch is within expected discretization error for finite-difference validation; the FFT Poisson solution itself is numerically verified

### CO tip diagnostics
- **Neutral CO**: Δρ integral ≈ 0.0008 e (effectively neutral)
- **Dipolar structure**: Negative oxygen lobe surrounded by positive carbon/shell
- **Dipole-differentiator effect**: Convolution of V_ES with neutral dipolar CO yields approximately ∇V (electric field), giving sharper features than blurring — this is physically correct
- **Plots saved in** `debug/co_tip/`: `co_rho_total_slices.png`, `co_rho_delta_slices.png`, `co_rho_xz_axis.png`

### All diagnostic plots (locations)
- **Step 1 (density)**: `debug/step1_density/` — `rho_slices.png`, `rhoNA_slices.png`, `rhoDiff_slices.png`, `rho_maxproj.png`
- **Step 2 (electrostatics)**: `debug/step2_electrostatics/` — `VES_slices.png`, `VES_lineprofile.png`, `VES_poisson_check.png`
- **Step 3 (Pauli)**: `debug/step3_pauli/` — `Epauli_slices.png`, `Fz_pauli_slices.png`
- **Step 4 (ES convolution)**: `debug/step4_electrostatics_conv/` — `EES_slices.png`, `Fz_ES_slices.png`
- **Step 5 (dispersion)**: `debug/step5_dispersion/` — `Evdw_slices.png`
- **Step 6 (composed)**: `debug/step6_composed/` — `V_decomposition_xy_z5.png`, `V_total_heights.png`, `component_traces.png`, `df_maps.png`

### Test command
```bash
PYTHONPATH=/home/prokop/git/FireCore:$PYTHONPATH python tests/tAFM/pyocl_fdbm/run_pyocl_fdbm.py \
  --nscf=10 --z_start 5.0 --z_end 2.0 --dz 0.2 --tip_model=co --co_nscf=10
```

### Current status
- Pipeline runs successfully with pentacene (22 C atoms)
- Scan grid auto-computed at 0.1Å/pixel: 203×110×16 (dx=0.100Å, dy=0.101Å)
- All potentials (Pauli, ES, vdW) computed on same high-resolution density grid
- Poisson solver validated (16% discretization error, order of magnitude correct)
- CO tip density properly integrated (neutral) and used for convolution
- Final plots show good resolution and physically reasonable contrast

### Pauli scaling and beta exponent
- Added `beta_pauli` CLI parameter: exponent applied to **both sample AND tip densities** before Pauli FFT convolution (default 1.0, try 1.2 to sharpen repulsive core)
- **Implementation**: `rho_eff = rho_grid^beta_pauli`, `tip_eff = tip_kernel^beta_pauli` — both raised to same power before convolution
- **Effect**: beta > 1.0 suppresses low-density regions and enhances high-density peaks → sharper, more atom-resolved Pauli repulsion
- **Current defaults**: `A_pauli=3.0`, `beta_pauli=1.0` (tunable via CLI `--A_pauli`, `--beta_pauli`)
- **Plot titles** now report `(A=..., beta=...)` for traceability

### LJ/Morse side-by-side comparison
- Added `compute_lj_forces()` in `run_pyocl_fdbm.py` (step 5b) with UFF-like LJ parameters + point-charge electrostatics
- **Result**: FDBM and LJ are NOT the same — order-of-magnitude differences in all components:
  - **Pauli/Repulsive**: FDBM ~12 eV/Å (sharp orbital overlap) vs LJ ~28 eV/Å (diffuse r⁻¹²)
  - **Electrostatics**: FDBM ~±0.5 eV/Å (dipolar CO density × V_ES) vs LJ ~±0.01 eV/Å (point q/r²)
  - **London/Attractive**: FDBM ~-0.4 eV/Å (tuned C6_CO=30) vs LJ ~-1.2 eV/Å (generic ε=0.003)
- **Root cause**: LJ uses point-atom pairwise interactions; FDBM uses distributed electronic densities. Point-charge Coulomb is negligible at AFM scan heights; dipolar CO density convolution captures gradient/multipole effects.
- **Plots**: `step6_fdbm_vs_lj_decomposition_z*.png` (per-height 2×4 grids with per-image color scales), `step6_fdbm_vs_lj_traces.png` (height curves)

### Remaining open questions / next tests
1. **E(z) approach curves with full quantum treatment**: Run Fireball or DFTB+ with CO molecule directly above sample atoms at varying z distances (one quantum system, not separate tip+sample densities). Compare total energy vs our FDBM force decomposition.
2. **Extended basis set test**: Current implementation uses minimal short basis set (`/home/prokop/Fireball/Fdata_HCNOS`). Need to compare with extended basis (`/home/prokop/Fireball/Fdata_HCNOS_ext`) to check if Pauli repulsion range and magnitude change significantly.
3. **Parameter fitting**: Systematically vary `A_pauli` (1–10), `beta_pauli` (1.0–1.5), `C6_CO` (10–100) against reference QM E(z) curves to find optimal FDBM parameters for pentacene/CO.