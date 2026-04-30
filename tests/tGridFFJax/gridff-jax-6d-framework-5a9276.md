# GridFF-JAX 6D Framework: Systematic Sampling, Fitting & MLIP Diagnostics   
**Comment:  this is not only limited to Ag or anu molecules, mad-surf is for the coign metal molecules @can see its paper if required, present here, and also itshould be extendable to any other MLIP. SO MAKE THAT ROBUST AND GENERIC FRAMEWORK. Also the sample generations shoulbe able to take molecules name as input from a file and then use it to generate samples for those molecules, you know what i mean right. Also we should have proper documentations and how to use this framework. as it is there**

Build a generic 6D (x,y,z,alpha,beta,gamma) sampling + fitting + diagnostic framework in JAX that replaces the current z-scan-only fitting pipeline, produces transferable per-element-pair parameters, and provides rigorous MLIP quality assessment tools.

---

## Diagnosis: Why Lateral Scans Fail

After reading the full codebase, the **root cause** is clear:

1. **Training data has NO lateral variation.** The current `benchmark_ag_zscan.py` builds poses at 3 fixed (u,v) sites (top/bridge/hollow) varying only z. The `build_pose_batch` in `rigid.py` adds tiny `jitter_uv=0.03` but this is noise, not systematic coverage.
2. **The fitted parameters encode only z-dependent physics.** REQ parameters (radius_offset, energy_scale) are optimized to reproduce z-profiles. They have no reason to capture lateral corrugation correctly.
3. **Lateral corrugation comes from PLQ grid shape, not parameters.** The PLQ grids themselves carry lateral variation, but the per-element coupling (REQ) was never trained to weight that variation correctly across the (x,y) plane.
4. **Force deviation from dE/dr**: If MAD-SURF forces don't exactly match numerical gradients of MAD-SURF energy, this inconsistency propagates into the fit. This hasn't been systematically checked.

---

## Implementation Plan

### Step 1: Run Baseline Diagnostics (Quick, ~30 min)

Produce reference numbers on this machine before changing anything.

- **1a.** Run `benchmark_ag_zscan.py` for CHONH2 with current best settings (damped PLQ + raw-r6) → z-scan RMSE baseline
- **1b.** Run `lateral_scan.py` with the z-scan-fitted params → lateral 1D/2D RMSE baseline
- **1c.** Compute MAD-SURF force-energy consistency: finite-difference dE/dx,dE/dy,dE/dz vs reported F at ~50 poses
- Save all baseline metrics to `tests/tGridFFJax/baseline_diagnostics/`

### Step 2: 6D Pose Sampler (`pyBall/gridff_jax/pose_sampling/sampler_6d.py`) — NEW FILE

A generic, configurable 6D rigid-body configuration generator.

**Dimensions:**
| Dim | Coordinate | Current | New |
|-----|-----------|---------|-----|
| 1,2 | x,y (lateral, fractional u,v) | Fixed at 3 sites + jitter | Systematic grid over unit cell |
| 3 | z (height) | Biased linspace at sites | Same, but at EVERY (u,v) point |
| 4,5,6 | alpha,beta,gamma (Euler / quaternion) | 4-8 systematic + random | Systematic Fibonacci sphere + yaw grid |

**Key design decisions:**
- **Lateral**: uniform (u,v) grid over one surface unit cell, configurable `n_u x n_v` (e.g., 10x10 = 100 lateral points)
- **Height**: biased linspace per (u,v) point, configurable `n_z` (e.g., 15-30 points)
- **Angular**: For monoatomic probes → identity only. For linear molecules → yaw grid + tilt grid. For polyatomics → Fibonacci sphere quaternions + yaw, configurable `n_orient`
- **Total poses**: `n_u * n_v * n_z * n_orient` — e.g., 10*10*20*8 = 16,000 poses
- **Extensibility hooks**: `SamplerConfig` dataclass with `extra_dimensions: list[ExtraDimension]` for future molecular deformation DOFs (bond stretch, angle bend, dihedral). Each `ExtraDimension` has `name`, `range`, `n_points`, `sampling_mode`

**Output**: `PoseBatch6D` — extends current `PoseBatch` with full 6D parametrization stored in `pose_params` and metadata tracking which dimension each index corresponds to.

### Step 3: Multi-Dimensional Training Pipeline (`pyBall/gridff_jax/fit/fit_6d.py`) — NEW FILE

Refactored fitting that uses 6D data.

- **3a.** Separate reference data generation from fitting: `generate_reference_data()` → saves `dataset_6d.npz` with energies + forces from teacher (MAD-SURF or DFT)
- **3b.** Train/validation/test split that is **stratified by (u,v) region and z-window**, not random — ensures every region of configuration space is represented in validation
- **3c.** Loss function includes all dimensions:
  ```
  L = w_E * MSE(E_pred - E_ref) + w_F * MSE(F_pred - F_ref) 
      + w_reg * regularization
  ```
  with optional per-dimension weighting (weight close-to-surface z more)
- **3d.** JAX vmap over poses for GPU parallelism
- **3e.** Support both single-molecule and multi-molecule joint fitting for transferability

### Step 4: Atom Type System for Transferable Parameters

Extend the current per-element parameters toward UFF-style per-atom-type parameters.

**Current**: `req_radius_offset["H"] = 0.5` — one value for ALL hydrogen atoms.

**New (Phase 1 — element-pair)**:
- Per-element parameters remain, but fitted jointly across multiple molecules
- `H` in CH4, H2O, NH3, CHONH2 all contribute to the same `H-Ag` parameters
- The fitting sees diverse chemical environments → parameters become transferable
- This is stored as `element_pair_params["H-Ag"]`, `element_pair_params["C-Ag"]`, etc.

**New (Phase 2 — atom-type, FUTURE scope annotation)**:
- Code structure allows `atom_type_params["H_alkyl-Ag"]`, `atom_type_params["H_hydroxyl-Ag"]`
- Atom typing from connectivity analysis (like UFF's `C_R`, `C_3`, etc.)
- This requires: (a) bond graph for the adsorbate, (b) typing rules
- The `atomtypes.ini` in `tests/tDFT_pentacene/` provides the per-element REQ data
- NOT implemented in first version — but the data structures will support it

**Implementation**: Add `AtomTypeConfig` dataclass and `atom_type_registry` module alongside existing `molecules.py`. The 6D sampler and fitter will accept either `element` or `atom_type` granularity.

### Step 5: MLIP Diagnostic & Benchmarking Tool (`pyBall/gridff_jax/diagnostics/`) — NEW PACKAGE

A standalone tool to rigorously assess any MLIP (MAD-SURF now, others later) and compare with DFT/GridFF.

**5a. `smoothness.py`** — Energy surface smoothness analysis
- Compute dE/dx, dE/dy, dE/dz via finite differences at configurable step sizes
- Compare with MLIP-reported forces
- Report: max |F_analytic - F_numerical|, RMS inconsistency, worst-case locations
- Heatmaps of inconsistency over (x,y) at fixed z, and along z-lines

**5b. `force_consistency.py`** — Force-energy derivative consistency
- For each pose: compute E(r), E(r+dx), E(r-dx) → numerical force
- Compare with MLIP analytical force
- Report per-atom, per-component statistics
- This is the "smoking gun" for MLIP smoothness problems

**5c. `corrugation.py`** — Lateral corrugation analysis
- At each z-height: compute E(u,v) over unit cell
- Report corrugation amplitude, symmetry, comparison between MLIP and GridFF
- Identify sites where MLIP shows unphysical features

**5d. `comparison.py`** — Multi-method comparison framework
- Generic interface: `MethodResult(energies, forces, metadata)`
- Compare any pair: MAD-SURF vs DFT, GridFF vs MAD-SURF, GridFF vs DFT
- Produce parity plots, error maps, z-resolved statistics, site-resolved statistics

**5e. `report.py`** — Automated report generation
- Aggregate all diagnostics into a single JSON + PDF summary
- Standard metrics: RMSE by dimension, by z-window, by site, by orientation

All tools use JAX where beneficial (vmap for batch finite differences).

### Step 6: Runner Scripts in `tests/tGridFFJax/`

- **`run_6d_sampling.py`** — Generate 6D reference dataset for a given molecule+substrate+teacher
- **`run_6d_fit.py`** — Fit transferable parameters from 6D dataset  
- **`run_6d_validate.py`** — Validate fitted parameters with held-out 6D data + lateral scans
- **`run_mlip_diagnostics.py`** — Run the full MLIP diagnostic suite
- **`run_transferability_test.py`** — Fit on Tier-1 molecules, test on held-out molecules

### Step 7: Integration with Existing GridFF Export

- Fitted parameters → `Bspline_PLQd.npy` via existing `export/firecore.py`
- The 6D-fitted parameters should produce grids that work with existing OpenCL/C++ runtime
- No changes to the C++ side — the improvement is entirely in how parameters are determined

---

## File Changes Summary

| Action | Path | Description |
|--------|------|-------------|
| **NEW** | `pyBall/gridff_jax/pose_sampling/sampler_6d.py` | 6D systematic pose generator |
| **NEW** | `pyBall/gridff_jax/fit/fit_6d.py` | Multi-dimensional fitting pipeline |
| **NEW** | `pyBall/gridff_jax/diagnostics/__init__.py` | Diagnostics package |
| **NEW** | `pyBall/gridff_jax/diagnostics/smoothness.py` | Energy surface smoothness |
| **NEW** | `pyBall/gridff_jax/diagnostics/force_consistency.py` | F vs dE/dr check |
| **NEW** | `pyBall/gridff_jax/diagnostics/corrugation.py` | Lateral corrugation analysis |
| **NEW** | `pyBall/gridff_jax/diagnostics/comparison.py` | Multi-method comparison |
| **NEW** | `pyBall/gridff_jax/diagnostics/report.py` | Automated report |
| **NEW** | `pyBall/gridff_jax/atom_types.py` | Atom type registry (future-ready) |
| **EDIT** | `pyBall/gridff_jax/config.py` | Add `Sampler6DConfig`, `DiagnosticsConfig` |
| **EDIT** | `pyBall/gridff_jax/interfaces.py` | Extend `PoseBatch` with 6D metadata |
| **NEW** | `tests/tGridFFJax/run_6d_sampling.py` | 6D dataset builder script |
| **NEW** | `tests/tGridFFJax/run_6d_fit.py` | 6D fitting script |
| **NEW** | `tests/tGridFFJax/run_6d_validate.py` | 6D validation script |
| **NEW** | `tests/tGridFFJax/run_mlip_diagnostics.py` | MLIP diagnostic runner |

---

## Execution Order

1. **Step 1**: Run baselines (~30 min compute, immediate insight)
2. **Step 2**: Build 6D sampler (core new code, ~2-3 hours)
3. **Step 3**: Build 6D fitting pipeline (~2 hours)
4. **Step 5**: Build diagnostics package (~2-3 hours)
5. **Step 4**: Atom type system (lighter, ~1 hour, mostly data structures)
6. **Step 6**: Runner scripts (~1 hour)
7. **Step 7**: Integration test with export pipeline (~30 min)

Total estimated: ~10-12 hours of implementation.

---

## Performance Notes

- **JAX on GPU/TPU**: All new code uses `jax.vmap` for batch evaluation. The 6D sampler generates poses as NumPy arrays (CPU), but evaluation + fitting runs on GPU via JAX.
- **Multi-CPU**: JAX's XLA backend automatically parallelizes across CPU cores when GPU is unavailable.
- **vs C++/OpenCL**: JAX JIT-compiled code is typically 2-10x slower than hand-tuned CUDA/OpenCL for simple kernels, BUT the fitting loop (gradient computation + parameter update) is where JAX excels — autograd + GPU makes it 50-100x faster than scipy on CPU. The runtime GridFF evaluation will still use the existing C++/OpenCL path.
- **Memory**: 16K poses × 6 atoms × 3 coords = ~300 KB per batch — trivial. The grids (~500 MB for PLQ) are the memory bottleneck, unchanged.

---

## Scope Boundaries

**IN scope (this plan):**
- 6D rigid-body sampling (x,y,z,alpha,beta,gamma)
- MAD-SURF as teacher backend
- Per-element parameters (H-Ag, C-Ag, N-Ag, O-Ag)
- MLIP diagnostic suite
- Tier-1 organic molecules (CH4, C2H4, C6H6, NH3, H2O, CO, CHONH2, CH2O, pyridine, furan, methanol, acetonitrile)
- Ag(111) substrate (with config hooks for Cu/Au/Pt)

**OUT of scope (future extension points annotated in code):**
- Molecular deformation dimensions (7D+) — data structure ready, sampling not implemented
- DFT teacher backend — interface defined, VASP single-point wrapper not implemented  
- Atom-type-aware parameters (C_R vs C_3) — registry defined, typing rules not implemented
- Multi-substrate joint fitting — config supports it, not tested
- Non-rigid substrate corrections
