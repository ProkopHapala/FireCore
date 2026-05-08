# PyOpenCL FDBM AFM Implementation

Self-contained PyOpenCL reimplementation of Full Density-Based Model (FDBM) AFM.

## Overview

Calls Fireball only for SCF density matrix extraction; all subsequent operations
(density projection, Poisson solver, Pauli/ES convolution, vdW, PP relaxation)
are implemented in Python with PyOpenCL/NumPy.

## Usage

```bash
PYTHONPATH=/path/to/FireCore:$PYTHONPATH python run_pyocl_fdbm.py \
    --nscf 50 --step 0.15 --nxy 40 40 \
    --z_start 5.0 --z_end 2.0 --dz 0.2
```

### Parameters

- `--nscf`: Fireball SCF iterations (default: 200)
- `--step`: Density grid step [Å] (default: 0.15)
- `--nxy`: Scan grid resolution (default: 60 60)
- `--z_start`, `--z_end`, `--dz`: Scan height range [Å]
- `--A_pauli`: Pauli amplitude (default: 16.0)
- `--C6_CO`: CO-tip C6 parameter (default: 30.0)
- `--q_CO`: CO-tip charge (default: -0.05)
- `--sigma_tip`: Tip Gaussian width [Å] (default: 0.7)
- `--tip_model`: `'gaussian'` (default) or `'co'` for true CO density from Fireball SCF
- `--co_nscf`: SCF iterations for CO tip density computation (default: 50)
- `--stop_after`: Stop after step N (1-5) for debugging

## Implementation Steps

### Step 1: Density Projection
- Run Fireball SCF on pentacene
- Project SCF density to 3D grid using PyOpenCL (`GridProjector`)
- Project neutral-atom density for delta_rho = rho_SCF - rho_NA
- **Validation**: Charge conservation check (∫ delta_rho dV ≈ 0)

### Step 2: Electrostatics (Poisson)
- Solve ∇²V = -4πρ via FFT on uniform grid
- **Validation**: Compare ∇²V with -4πρ via finite difference

### Step 3: Pauli Repulsion
- **Gaussian model**: Build normalized Gaussian tip (σ=0.7 Å)
- **CO model** (`--tip_model=co`): Run Fireball SCF on CO molecule, project total electron density to grid. O is the probe center (origin). FFT convolution requires the tip to be centered at array index 0.
- FFT convolution: E_pauli = A * ∫ ρ_tip * ρ_sample dV
- Compute forces via FFT gradient: F = -∇E

### Step 4: Electrostatics Tip-Sample
- **Gaussian model**: Gaussian tip scaled by q_CO
- **CO model**: Use CO delta-density (ρ_SCF - ρ_NA) as tip charge distribution. For neutral CO, integral ≈ 0, but spatial charge redistribution (dipole) creates ES contrast.
- FFT convolution: E_es = q_CO * ∫ V_ES * ρ_tip dV (Gaussian) or E_es = ∫ V_ES * Δρ_CO dV (CO)
- Compute forces via FFT gradient

### Step 5: Dispersion (vdW)
- Pairwise C6/r^6 with regularization (RA2 = 1.5² Å²)
- Prevents divergence at close distances

### Step 6: Composition + PP Relaxation
- Interpolate Pauli/ES forces from density grid to scan grid
- Sum with vdW forces (already on scan grid)
- PP relaxation using FIRE-like algorithm (`pp_relax_2d`)
- Compute frequency shift: df = -dFz/dz

## Critical Issues Found & Fixed

1. **vdW Divergence**: Pure C6/r^6 blows up at z < 2 Å. Fixed with RA2=2.25 Å² regularization.
2. **Sign Convention**: Mixed energy-gradient vs force storage caused sign flip. Fixed by storing ∇E consistently everywhere, computing F = -∇E only in composition step.
3. **np.gradient Sign**: `probe_heights` is decreasing, so `np.gradient` needs actual (negative) spacing.
4. **Charge Conservation**: Grid projection with 0.2 Å step gives ~0.2e residual. Use ≤0.15 Å for production.
5. **Plot Range Problem**: Using global vmin/vmax across entire 3D array made XY slices at distant z invisible. Fixed with per-slice adaptive ranges (`per_slice_range=True`).
6. **FFT Tip Centering**: CO tip projected to grid center, but FFT convolution expects tip at array index 0. Fixed by circularly shifting CO density to index 0 to match `build_gaussian_tip()` convention.

## Validation Results

### Gaussian Tip Model (Pentacene, 40×40, z=5.0-2.0 Å)

| Component | Fz Range [eV/Å] | Character |
|---|---|---|
| Pauli | [0.0, 2.24] | Repulsive, z < 3.2 Å |
| vdW | [-0.36, -0.005] | Attractive, z > 3.5 Å |
| ES | [-0.008, 0.007] | Negligible (q_CO=-0.05) |
| Total | [-0.037, 1.88] | Balance ~3.2-3.5 Å |

### True CO Tip Model (Pentacene, 40×40, z=5.0-2.0 Å)

| Component | Fz Range [eV/Å] | Character |
|---|---|---|
| Pauli | [0.0, 59.0] | Much stronger than Gaussian; short-ranged |
| vdW | [-0.35, -0.006] | Similar to Gaussian |
| ES | [-0.48, 0.40] | Bipolar, from CO dipole |
| Total | [-0.03, 58.7] | Pauli dominates at z < 3.5 Å |

**Key Difference**: True CO density is compact (atomic-like) so Pauli falls off faster than the broad Gaussian. At z=5.0Å Pauli is negligible (~0.02 eV); at z=2.0Å it reaches ~120 eV. The Gaussian model overestimates Pauli range because σ=0.7Å is artificially broad.

## CO Tip Diagnostics

CO density is computed in a separate subprocess (`compute_co_tip.py`) because Fireball cannot reallocate for different molecules.

- `debug/co_tip/co_rho_total_slices.png` - Total electron density of CO
- `debug/co_tip/co_rho_delta_slices.png` - Delta density (SCF - NA) showing charge transfer C→O
- `debug/co_tip/co_rho_xz_axis.png` - xz cut through molecular axis

CO delta-density shows strong charge polarization: negative (electron accumulation) around O, positive (depletion) around C.

## Debug Outputs

All steps save `.npy` arrays and `.png` plots to `debug/step{1-6}/`:
- `step1_density/`: rho_grid, rho_na_grid, rho_diff, slices, line profiles
- `step2_electrostatics/`: V_ES, Poisson validation plots
- `step3_pauli/`: E_Pauli, F_Pauli, slices
- `step4_electrostatics_conv/`: E_ES, F_ES, slices
- `step5_dispersion/`: E_vdw, F_vdw, slices
- `step6_composed/`: V_total, V_pauli, V_es, V_vdw (potential energy maps), Fz_raw, Fz_relax, df, component traces, df maps
  - `step6_V_decomposition_xy_z5.png` - Potential energy decomposition at z=5.0Å
  - `step6_V_total_heights.png` - V_total at several scan heights
  - `step6_component_traces.png` - Force components vs height at scan center

## Next Steps

- [ ] Implement proper LJ repulsive wall instead of simple RA2 regularization
- [ ] Test with polar molecules where ES is significant
- [ ] Compare with classical Morse/LJ reference on same system
- [ ] Optimize density grid size for faster computation
- [ ] Add tip tilt effects (not just vertical probe)
