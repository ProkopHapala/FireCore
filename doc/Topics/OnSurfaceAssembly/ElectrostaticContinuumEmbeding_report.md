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
