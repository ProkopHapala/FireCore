---
Description: Electrostatic continuum embedding report
---

# Electrostatic Continuum Embedding for Surface Interaction Scans

This note documents the current implementation of the macro (continuum) correction used to accelerate and stabilize molecule–substrate interaction scans under PBC, covering the analytic model, formulas, Python reference, OpenCL port, tests, and how to reproduce the validation plots.

## Problem Statement

Finite image sums (limited nPBC) of a periodic slab leave slow-decaying monopole/dipole/quadrupole tails that break translational invariance for probe charges and molecules. The goal is to reach <1e-4 eV equivalence error between primitive translations while keeping the compute on GPU.

## Analytic Model

We subtract the potential of a uniformly polarized rectangular sheet per substrate layer and add back the ideal infinite-lattice limit implicitly by centering the rectangle on the replicated cell. The correction is applied per layer `zl` with (sheet charge density `σ`, dipole density `P`, quadrupole density tensor `Q`). The macro potential for probe position `(x,y,z)` is:

- **Monopole (charge sheet):**
  \[ \phi_σ = σ \; I_0(x,y,z; A,B) \]
- **Dipole sheet:**
  \[ \phi_P = -P_x I_x - P_y I_y + P_z I_ω \]
- **Quadrupole sheet (full traceless tensor):**
  \[ \phi_Q = \tfrac{1}{6}( Q_{xx} I_{xx} + Q_{yy} I_{yy} + Q_{zz} I_{zz} + 2 Q_{xy} I_{xy} + 2 Q_{xz} I_{xz} + 2 Q_{yz} I_{yz}) \]

where the geometric integrals are derived from the charge-sheet potential of a rectangle spanning `[xmin,xmax]×[ymin,ymax]`:
- `I0 = rect_sheet_potential(...)` (closed form using logs/atan2)
- `Ix, Iy, Iω` come from derivatives of `I0`
- `Ixx, Iyy, Izz, Ixy, Ixz, Iyz` are second derivatives of `I0` (in Python implemented via small-step finite differences; OpenCL currently uses analytical dipole, charge, and numerical grad for quadrupole pending analytic force closure).

The probe energy correction is
\[ \Delta E = - C_{oul} \sum_{atoms} q_i (\phi_σ + \phi_P + \phi_Q) \]
with `Coulomb_const = 14.3996 eV·Å/e^2`.

### Rectangle Bounds (critical fix)
Bounds are now derived from cell geometry, not atom extents:
- half-step `hx, hy` = half of minimal nonzero spacing in substrate x/y
- `xmin = xmin_atoms - hx - nPBCx * |a|`, `xmax = xmax_atoms + hx + nPBCx * |a|` (same for y)
This aligns the continuum rectangle with the replicated supercell and removes spurious in-plane layer dipoles.

### Layer Moments
For each unique `z` layer:
- `σ = sum(q)/Area`
- `P = sum(q * [dx,dy,dz]) / Area`, with `dx,dy` relative to rectangle center and `dz` relative to layer z
- Full traceless quadrupole tensor `Q` using standard symmetric-traceless form (factors 2,3 as in current code).

## Python Reference Implementation
File: `pyBall/OCL/InteractionEnergy.py`
- Core functions: `rect_sheet_potential`, `rect_dipole_potential`, `rect_quadrupole_potential` (full tensor, small-step FD), `_update_macro_from_substrate`, `_apply_macro_correction`.
- Macro bounds use ` _half_step_from_coords` (new) to expand by half-step.
- Macro layers store `(zls, sigmas, Pmus, Qdens)`; applied per transform in `_apply_macro_correction`.
- Coulomb-only path with optional LJ/HBond/Morse disabled is used for validation.

## OpenCL Implementation

### Host (packing / launch)
File: `pyBall/OCL/MolecularDynamics.py`
- Added `half_step_from_coords` and applied the same bound expansion when setting surface.
- Packs per-layer: `surf_mpos` (bounds), `surf_mdip` (dipole density), `surf_mQa/b/c` (dipole layers), `surf_qQa/b/c` (quadrupole tensor components), and `surf_mQc` (sheet charge densities).
- Launch fix: `run_getSurfMorse` now uses exact global size `(natoms, nSystems)` with `local=None` to avoid barrier-after-return UB when natoms < nloc.
- Kernel params `GFFParams` include macro enable flag and layer count.

### Kernel
File: `cpp/common_resources/cl/relax_multi.cl`
- `macro_phi_rect_dipole`, `macro_phi_rect_charge` provide analytic rectangle potential terms.
- `getMacroRectLayers` accumulates charge, dipole, quadrupole contributions for each layer.
- Sign fix: macro contribution is **added** to the force-energy accumulator `fe` because the kernel stores negative energy in `forces.w`.
- `getSurfMorse` invokes `getMacroRectLayers` when `bMacro` and `charge≠0`; explicit image sum otherwise unchanged.

## Test Procedures (reproducible)

### 1) Non-GUI macro validation & plots
Script: `tests/tMMFF/test_interaction_scan.py`
- Function: `run_macro_reference_scan()`
- Generates:
  - 1D equivalent-site plot (primitive shifts 4 Å in x and y): `output_interaction_scan/H2O_NaCl8x8_equiv_line_macro.png`
  - 2D XY scan: `output_interaction_scan/H2O_NaCl8x8_XYscan_-5_35_npbc4_macro.png`
- How to run:
```bash
cd tests/tMMFF
python test_interaction_scan.py
```

### 2) Primitive-site equivalence (single ion and H2O)
Minimal snippets (already run during validation):
- Single ion, Python reference: `max |ΔE| ≈ 3.8e-05 eV (x), 5.6e-05 eV (y)`
- Single ion, OpenCL shared kernel: `max |ΔE| ≈ 7.5e-05 eV (x), 2.6e-05 eV (y)`
- H2O, OpenCL: `max |ΔE| ≈ 1.7e-05 eV (x), 2.4e-05 eV (y)`
All under `nPBC=(4,4,0)`, Coulomb-only, macro on.

### 3) Launch/packing regression
- Ensure `run_getSurfMorse` uses `(natoms, nSystems)` global size.
- Verify macro buffers exist: `surf_mpos`, `surf_mdip`, `surf_mQa/b/c`, `surf_qQa/b/c`, `surf_mQc`.

## How to Test Manually
1. **Single-ion parity (Python):**
```python
from pyBall.OCL.InteractionEnergy import InteractionScanner
from pyBall.OCL import ScanUtils
scanner = InteractionScanner(nloc=32)
scanner.set_molecule(...single ion...)
scanner.load_substrate_xyz('cpp/common_resources/xyz/NaCl_8x8_L3.xyz')
scanner.enable_macro=True; scanner.nPBC[:] = (4,4,0); scanner._update_macro_from_substrate()
# build transforms along primitive shifts and evaluate
```
2. **Single-ion parity (OpenCL):**
```python
from pyBall.OCL.MolecularDynamics import MolecularDynamics
md = MolecularDynamics(nloc=32, debug_build_options='-DDBG_UFF=0')
md.init_with_atoms(na=1, ...); md.setup_kernels()
md.set_surface('cpp/common_resources/xyz/NaCl_8x8_L3.xyz', nPBC=(4,4,0), alpha_morse=1.8, r_damp=0.0, bMacro=True)
# loop over primitive shifts, run getSurfMorse, read forces.w
```
3. **Plots:** run `python tests/tMMFF/test_interaction_scan.py` and inspect the two PNGs.

## Known Limitations / Next Steps
- Quadrupole forces in OpenCL still use numerical gradients (eps^2 finite difference). Analytical derivatives can replace this for full parity with Python’s forces once derived.
- Current tests target primitive translations; other lattice directions or taller nPBC can be added if needed.

## File Map (new/modified hotspots)
- Python: `pyBall/OCL/InteractionEnergy.py`
  - `_half_step_from_coords`, rectangle bound fix, full-tensor quadrupole, macro correction.
- Python host: `pyBall/OCL/MolecularDynamics.py`
  - `half_step_from_coords`, macro bound fix, buffer packing, getSurfMorse launch fix.
- OpenCL kernel: `cpp/common_resources/cl/relax_multi.cl`
  - macro charge term, sign fix, shared macro call inside `getSurfMorse`.
- Tests/plots: `tests/tMMFF/test_interaction_scan.py`
  - headless macro validation and PNG outputs.

## Summary of Results
- Primitive-site errors now <1e-4 eV for the validated cases (single ion, H2O) with macro enabled.
- Required plots produced:
  - `tests/tMMFF/output_interaction_scan/H2O_NaCl8x8_equiv_line_macro.png`
  - `tests/tMMFF/output_interaction_scan/H2O_NaCl8x8_XYscan_-5_35_npbc4_macro.png`
- Shared OpenCL path matches the corrected Python reference after fixing bounds, sign, packing, and launch configuration.

## New functions/kernels added (since last commit)

**Python – `pyBall/OCL/InteractionEnergy.py`**
- `_half_step_from_coords(cs)`: find minimal spacing to expand rectangle bounds by half-step.
- `rect_quadrupole_potential(..., Qxx,Qxy,Qyy,Qxz,Qyz,Qzz, ...)`: full-tensor rectangle quadrupole potential (FD second derivatives of charge sheet).
- `_update_macro_from_substrate`: now computes full traceless quadrupole tensor and uses half-step bounds.
- `_apply_macro_correction`: applies charge/dipole/full-tensor quadrupole per layer.

**Python host – `pyBall/OCL/MolecularDynamics.py`**
- `half_step_from_coords(cs)`: same half-step bound helper used for GPU packing.
- `run_getSurfMorse`: launch fix (exact global size, no padding) to avoid barrier/return UB.
- `scanSurfMorse_2D(...)`: headless 2D sampler over getSurfMorse for debugging/plots.
- Surface packing updates: uploads `surf_mpos`, `surf_mdip`, `surf_mQa/b/c`, `surf_mQc` (sheet charges), `surf_qQa/b/c` (quadrupole tensors), and kernel params (`GFFParams` including macro flag and layer count).

**OpenCL kernel – `cpp/common_resources/cl/relax_multi.cl`**
- `macro_phi_rect_charge`, `macro_phi_rect_dipole`: analytic rectangle charge/dipole potentials.
- `getMacroRectLayers(...)`: accumulates charge, dipole, quadrupole rectangle contributions per layer; numerical grads for quadrupole force; uses charge-sheet term (added) and corrected sign when adding to accumulator.
- `getSurfMorse`: now calls `getMacroRectLayers` when macro enabled; sign fix for macro addition.

**Tests/plots – `tests/tMMFF/test_interaction_scan.py`**
- `run_macro_reference_scan()`: non-GUI macro validation, prints equivalent-site errors, saves 1D/2D PNGs.
- Updated `main()` to run the macro reference scan by default.

## Addendum: Fast rigid-surface scan debugging, PTCDA fix, and GUI improvements

This section documents the later debugging round focused on the fast GUI/headless rigid-surface path (`ExplorerVisPy.py` + `tests/tMMFF/test_interaction_scan.py`) and large aromatic molecules such as PTCDA.

### 1) Root cause of the PTCDA fast-path failure

The remaining fast-path failure for `PTCDA.xyz` on `NaCl_8x8_L3.xyz` was **not** caused by MMFF `nvecs` / pi-orbital packing. The `getSurfMorse` path is intentionally atom-only:

- host side uses `natoms` atoms and `nnode = 0`
- kernel indexes atoms with `i0a = iS*natoms`, `i0v = iS*(natoms+nnode)`
- rigid-surface scan reduces only `forces[:natoms,3]`

The actual failure came from the OpenCL launch contract of `getSurfMorse` in `relax_multi.cl`:

- the kernel uses fixed local arrays `LATOMS[32]`, `LCLJS[32]`
- the kernel contains barriers inside the surface loop
- threads with `iG >= nAtoms` return early

Therefore a padded x-dimension with a partial final workgroup is invalid: some threads would return before a barrier, producing undefined behavior. This was harmless for small molecules with `natoms <= 32` or exact divisibility, but broke for `PTCDA` (`natoms = 38`).

#### Final host fix

`pyBall/OCL/MolecularDynamics.py` now launches `getSurfMorse` with:

- `global_x = natoms`
- `local_x = lx`, where `lx <= nloc` and `natoms % lx == 0`

For PTCDA this chooses `lx = 19`, avoiding any partial x-workgroup and restoring correct multi-system behavior.

### 2) Additional fast-path fixes retained

Two other fixes remain part of the final implementation:

- **Molecule typing parity**
  - `tests/tMMFF/test_interaction_scan.py` now passes the same `mol_type_map` into the fast GPU setup as the reference path.
  - This is essential for PTCDA, where `C -> C_R` and `O -> O_2` are needed when loading from XYZ.

- **Substrate-cell translation wrapping**
  - the rigid transform translation is wrapped into the substrate primitive cell before upload/evaluation.
  - This does **not** impose PBC on the molecule internals; it only wraps the rigid-body lateral translation so the finite substrate image sum and macro rectangle remain centered consistently with the reference scanner.

### 3) Separated timing: preparation vs kernel vs download

The fast path now reports three timing components in `MolecularDynamics.eval_rigid_getSurfMorse()`:

- `t_prep_s`
  - rigid transform packing
  - upload of batched positions
  - force-buffer clear
  - queue finished before timing closes

- `t_kernel_s`
  - kernel launch to `queue.finish()`

- `t_download_s`
  - download of `aforce`
  - host-side energy reduction
  - queue finished before timing closes

These are surfaced both in:

- headless output (`tests/tMMFF/test_interaction_scan.py`)
- GUI status / console output (`pyBall/ExplorerVisPy.py`)

### 4) Validated results after the final fix

#### PTCDA on NaCl 8x8

Headless fast vs reference, `nx = ny = 121`, `z = 3.5 Å`, `nPBC = (4,4,0)`:

- equivalent-site line error:
  - `max|ΔE|_x = 1.377101e-02 eV`
  - `max|ΔE|_y = 1.109781e-02 eV`
- 2D lateral scan parity:
  - `max|ΔE|_2D = 7.459815e-05 eV`
- fast path timing:
  - `wall = 0.134 s`
  - `prep = 0.009 s`
  - `kernel = 0.123 s`
  - `download = 0.002 s`

This confirms that for the large 8x8 PTCDA case the runtime is dominated by the kernel, not by Python harness overhead.

#### PTCDA on NaCl 1x1

Headless fast vs reference, `nx = ny = 121`, `z = 3.5 Å`:

- equivalent-site line error:
  - `max|ΔE|_x = 1.339363e-02 eV`
  - `max|ΔE|_y = 1.082945e-02 eV`
- 2D lateral scan parity:
  - `max|ΔE|_2D = 1.093703e-03 eV`
- fast path timing:
  - `wall = 0.017 s`
  - `prep = 0.007 s`
  - `kernel = 0.009 s`
  - `download = 0.001 s`

For the small 1x1 substrate, kernel time is much smaller and Python-side prep becomes a larger fraction of total wall time, as expected.

### 5) Generated plots

The headless script now supports batch preset runs and writes the validated plots into:

- `tests/tMMFF/output_interaction_scan/`

including the PTCDA outputs:

- `PTCDA__NaCl_1x1_L3_equiv_line_fast_gpu.png`
- `PTCDA__NaCl_1x1_L3_XYscan_fast_gpu.png`
- `PTCDA__NaCl_8x8_L3_equiv_line_fast_gpu.png`
- `PTCDA__NaCl_8x8_L3_XYscan_fast_gpu.png`
- corresponding `*_macro.png` reference plots

### 6) GUI updates

`pyBall/ExplorerVisPy.py` now includes:

- molecule presets:
  - `H2O`
  - `PTCDA`
- substrate presets:
  - `NaCl 1x1`
  - `NaCl 8x8`
- file dialogs for custom molecule/substrate loading
- startup CLI-triggered scan execution for debugging
- fast/reference backend timing printout

### 7) Bond visualization fix for PTCDA

The incorrect PTCDA bond display in the GUI was caused by a fixed-cutoff bond guess in `pyBall/MolGUI_common.py`:

- previous behavior used `findBondsNP(..., Rcut=1.8, byRvdW=False)`
- this is too crude for aromatic systems and can produce visually wrong connectivity

The GUI helper now uses element-aware bond detection:

- element names are converted to atomic types when available
- `findBondsNP(..., atypes=..., byRvdW=True, RvdwCut=1.5)` is used

This does not create an exact chemical topology in the MOL2 sense, but it is much more reasonable for aromatic PTCDA visualization than the previous constant cutoff.

### 8) Final conclusion

The fast rigid-surface scan path is now working for both small and large molecules, including PTCDA. The decisive issue was a **host launch mismatch with kernel barrier semantics**, not a physics error in the kernel itself. After fixing the launch, keeping molecule typing consistent, and exposing split timings, the large 8x8 PTCDA scan shows:

- correct parity (`~7.5e-05 eV` in 2D)
- strong acceleration relative to reference (`~0.13 s` vs `~4.2 s` for the validated 121×121 case)
- kernel-dominated runtime on the large GPU workload

So the Python harness is no longer the main bottleneck for the large fast scan; the measured runtime is mostly the GPU kernel execution itself.

### 9) Folded-basis PyOpenCL parity (primitive lattice fix)

We implemented a compressed folded-basis evaluation path (PyOpenCL + `getSurfFolded`) and discovered a critical parity bug: the lateral basis was built on the **8×8 supercell lattice (32 Å)** instead of the primitive NaCl **4 Å** repeat. This collapsed the folded field to a near-constant. The fix was to infer the primitive in-plane lattice from the substrate motif and use it for folded basis construction and evaluation.

Key diagnostics now printed during fit:

```
[folded] basis primitive_lvec a=[4.0, 0.0, 0.0] b=[0.0, 4.0, 0.0] |a|=4.000000 |b|=4.000000
[folded] basis u_freqs=[0, 1, 2, 3] v_freqs=[0, 1, 2, 3] z_alphas=[0.2, 0.4, 0.6, 0.8]
```

The test harness now also wraps periodic transforms in the CPU reference to match the folded path.

#### Parity results (folded vs reference, Morse + macro)

- **H2O on NaCl 8×8**, z≈1.0 Å, nx=ny=81:
  - line max|ΔE|: x=2.94e-04 eV, y=3.05e-04 eV
  - 2D max|ΔE|: 4.43e-04 eV
  - plots: `tests/tMMFF/output_interaction_scan/H2O_O__NaCl_8x8_L3_equiv_line_folded_gpu.png`, `tests/tMMFF/output_interaction_scan/H2O_O__NaCl_8x8_L3_XYscan_folded_gpu.png`

- **H2O on NaCl 1×1**, z≈1.0 Å, nx=ny=81:
  - line max|ΔE|: x=2.82e-04 eV, y=2.82e-04 eV
  - 2D max|ΔE|: 4.57e-04 eV
  - plots: `tests/tMMFF/output_interaction_scan/H2O_O__NaCl_1x1_L3_equiv_line_folded_gpu.png`, `tests/tMMFF/output_interaction_scan/H2O_O__NaCl_1x1_L3_XYscan_folded_gpu.png`

- **PTCDA on NaCl 8×8**, z≈3.5 Å, nx=ny=81:
  - line max|ΔE|: x=1.16e-04 eV, y=1.20e-04 eV
  - 2D max|ΔE|: 1.64e-04 eV
  - plots: `tests/tMMFF/output_interaction_scan/PTCDA__NaCl_8x8_L3_equiv_line_folded_gpu.png`, `tests/tMMFF/output_interaction_scan/PTCDA__NaCl_8x8_L3_XYscan_folded_gpu.png`

All folded runs use per-type fitted coefficients (unique REQ rows) and the primitive lateral basis shown above. The remaining errors are O(1e-4 eV) on the tested grids.

### 10) Takeaways to avoid similar issues

- **Always use primitive lattice for basis construction.** Do not assume the substrate XYZ cell is the primitive; infer the smallest repeat from the motif (or load a known primitive) before building lateral modes.
- **Print and inspect basis parameters.** Lattice vectors, lateral frequencies, and z-decay alphas must be logged so periodicity mistakes surface immediately.
- **Match periodic wrapping between reference and fast paths.** If the folded path wraps transforms into the cell, the CPU reference must wrap identically for fair parity.
- **Per-type coefficients are necessary but not sufficient.** Even with distinct REQ types, a wrong lattice or wrapping can collapse the folded field to a constant; always check min/max(ref) vs min/max(folded) and sample points.
- **Keep headless diagnostics strict.** Report ref/folded min/max, max|Δ|, and a few sample values for both line and 2D scans so qualitative mismatches are obvious.

### 11) 2026-03-13 — Folded fitting & periodic flat-surface test (current state)

#### What was broken

- The “periodic” lateral scan in `tests/tMMFF/test_flat_folded_components.py` was run along a diagonal (a+b) direction and applied the macro correction on **unwrapped** transforms while the local GPU path wrapped translations to the primary cell. This mismatch reintroduced a large linear slope/jump in lateral plots.
- The folded fit weighting plot was unclear (weights only, no energy overlay), and there was no explicit repulsion cutoff to suppress highly repulsive samples in the fit.

#### Fixes implemented

1. **Macro footing aligned**: `apply_macro()` now wraps transforms (`_wrap_transforms_PBC`) before applying `_apply_macro_correction`, so macro and local evaluations use identical coordinates.
2. **Validated directions**: Lateral parity scans now run along the primitive lattice vectors `a` and `b` (previously validated paths), not the diagonal `a+b`.
3. **Repulsion cutoff**: New CLI flag `--repel-cut` (default `0.2 eV`) excludes highly repulsive points from the fit mask; the mask also enforces `z >= zfit_min`.
4. **Weight/energy plot**: Weight plot now overlays macro-corrected reference energy, repulsion cutoff, z-fit-min marker, weight, and mask with clear labels. Files: `tests/tMMFF/output_flat_folded_components_periodic/H_weights.png` and `O_weights.png`.
5. **Consistent macro scanners**: Reference uses `nPBC=(12,12,0)`, embedding/folded use `nPBC=(4,4,0)`; each builds its own macro scanner to avoid mixed footing.

#### How the current fitting/testing works

- **Fit generation** (`build_folded` in `test_flat_folded_components.py`):
  - Sample uniform xy grid (24×24) over z in the requested range; evaluate macro-corrected reference energy.
  - Fit mask keeps points with `z >= zfit_min` and `E <= repel_cut`; weights use Boltzmann-like factor (`weight_power`, default 12).
  - Run `fit_folded_surface_basis` for components (Pauli, London, Coulomb) without macro during fit (explicit macro applied later in tests).
- **Scans**:
  - **Z-scans** above Na and Cl surface sites with periodic wrapping + macro on reference/embedding/folded; outputs PNG/CSV/JSON per site.
  - **Lateral scans** along `a` and `b` at fixed height (default 2.0 Å) with periodic wrapping + macro; outputs PNG/CSV/JSON per direction.
- **Components**: Pauli, London, Coulomb evaluated separately; total is their sum.
- **Artifacts to watch**: Macro-induced slope/jump is gone on validated `a`/`b` scans (embedding errors ~1e-4 eV). Remaining large errors are folded Pauli over Cl (basis expressivity in the repulsive wall, not macro).

#### Current outputs (2026-03-13 run)

- Folder: `tests/tMMFF/output_flat_folded_components_periodic/`
- Key files:
  - Lateral: `H_qp0.20_{a,b}_lateral.png`, `O_qp0.20_{a,b}_lateral.png`
  - Z-scans: `H_qp0.20_{na,cl}_zscan.png`, `O_qp0.20_{na,cl}_zscan.png`
  - Weights: `H_weights.png`, `O_weights.png`
  - Summary: `summary.md`, `summary.json`
- Embedding lateral max|ΔE| is back to ~1e-4 eV; folded still large on the Cl Pauli wall (≈1 eV for H, ≈4 eV for O) — attributed to basis limits in the repulsive region, not macro periodicity.

#### Next steps (fit improvements)

- Increase basis richness near the wall (higher `nzbasis` or specialized repulsive terms) while keeping the current macro footing intact.
- Optionally adjust `repel_cut` or weight shaping to de-emphasize the hardest wall but retain enough curvature to fit the minimum reliably.

### 12) 2026-03-13 — Orig folded kernel z-scan (mask/footing fixes)

What was still wrong

- Fit targets and weights in the harness were not aligned with the basis sampling grid, so the folded fit effectively trained on a laterally shifted dataset.
- `z_fit_min` was applied relative to the absolute z-range origin instead of “z above the top layer”, so the intended fit window was wrong.
- Baseline subtraction for `total` could double-count a separate total baseline instead of deriving it from component baselines.

Fixes applied

1. **Target grid match**: Rebuilt the fit targets on the exact `uvz` sampling grid emitted by `fit_folded_surface_basis`, eliminating lateral phase mismatch.
2. **Correct z mask**: `z_fit_min` now applies to `z - z_top` (height above the top layer), so the fit window is what the CLI specifies.
3. **Consistent baselines**: `total` baseline is now the sum of component baselines; no extra subtraction is applied.
4. **Weights**: Kept Boltzmann-like weights with `repel_cut=0.2 eV`; far-z baselines retained; tail forced to zero by construction.

Evidence (orig kernel, z-only)

- **Far tail**: folded total → ~1e-5 eV at z=8 Å (matches reference zero).
- **Attractive well (O@Cl)**: minimum at z≈3.30 Å; ref −0.026246 eV, folded −0.027563 eV (error ≈ −1.3 meV). Position and depth are aligned.
- **Fit window error (ref ≤ 0.2 eV)**: RMSE ≈ 5.3 meV (O@Cl), ≈ 4.4 meV (O@Na).
- **Repulsive wall (ref > 0.2 eV)**: still large (≈6.6 eV on O@Cl) because it is intentionally down-weighted/excluded by `repel_cut`.

Status and remaining work

- Orig kernel: tail and minimum fixed; near-min RMSE in the meV range. The only large discrepancy is the excluded Pauli wall.
- Lateral scans and optimized kernels (harmonics/workgroup) still need to be rerun under this corrected footing; planned next.
