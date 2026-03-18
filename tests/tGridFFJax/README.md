# tGridFFJax

JAX/GPU distillation workflow for building FireCore-style `GridFF` grids for rigid metal slabs from either DFT volumetric data or a teacher potential such as MAD-SURF.

This directory is for the new workflow only. The original runtime in FireCore remains unchanged.

## What This Is

The target workflow is:

1. read slab volumetric data from `CHGCAR` and `LOCPOT` or an ML density backend
2. build substrate `P/L/Q` fields on a 3D grid
3. sample rigid adsorbate poses above the slab
4. label those poses with a teacher potential such as MAD-SURF
5. fit a fast student model in JAX
6. export a static `Bspline_PLQd.npy` that FireCore can load like the original GridFF

For now, the recommended path is:

- focus on `PLQ` only
- keep the slab rigid
- use fixed adsorbate charges
- keep `sample_shift_z` and `coulomb_shift_z` at zero for Ag strict `PLQ`
- use cubic B-spline coefficients and cubic B-spline sampling

The extra `CT/image/reactive` terms exist in the codebase, but they are not yet production-ready for chemisorbed Ag systems.

## Architecture

```text
         DFT / ML density side                           teacher side
  CHGCAR + LOCPOT or ML rho/V(r)                   MAD-SURF / future MLIP
               |                                           |
               v                                           |
   density_backends/vasp_volumetric.py                     |
               |                                           |
               v                                           v
     substrate_fields.py builds zyx grids     pose_sampling/ rigid SE(3) poses
      parity_core or metal_dft_plq                     (u, v, z, quat)
        Pauli / London / Coulomb
               |                                           |
               +--------------------+----------------------+
                                    |
                                    v
                     hybrid_energy/model.py
          JAX student with cubic B-spline grid sampling
             and optional CT/image/reactive extensions
                                    |
                                    v
                         fit/optimize.py
                    Optax + JAX on GPU or CPU
                                    |
                +-------------------+------------------+
                |                                      |
                v                                      v
   export/firecore.py                         validation/ + plots
   Bspline_PLQd.npy                           parity / error / z-scans
   ReactiveChannels.npy                       grid slices / timings
   MetalResponse.npz
   GridMeta.json
```

## Training vs Runtime

This is the key point.

MAD-SURF is needed during grid generation and validation, not during later GridFF usage.

The intended sequence is:

1. generate teacher labels with MAD-SURF for many rigid poses
2. fit the GridFF student once
3. export `Bspline_PLQd.npy`
4. use that exported grid later without MAD-SURF

For polyatomic rigid adsorbates, the teacher force used for fitting must be the interaction force:

`F_interaction = F_slab+mol(on adsorbate atoms) - F_isolated_molecule(rotated into the same pose)`

That subtraction is now implemented in `teacher_backends/madsurf.py`. Without it, `CO` and `H2O` force labels include internal gas-phase restoring forces that the surface GridFF cannot reproduce.

After step 3, the static `PLQ` grid can be used like the original FireCore GridFF for:

- dragging molecules on the surface
- translations and rotations
- rigid-surface manipulations
- single-molecule or multi-molecule runs, if you also provide molecule-molecule and intramolecular force terms elsewhere

What the exported static grid does not yet replace:

- molecule-molecule interactions between adsorbates
- bond flexibility inside the adsorbate
- dynamic charge transfer
- image-charge response at runtime

So for multiple molecules on the surface:

- the substrate part can come from the generated GridFF
- intramolecular and intermolecular terms still need MMFF or another force field
- the surface is assumed rigid

## Current Code Layout

Implementation package:

- `pyBall/gridff_jax/density_backends/`
- `pyBall/gridff_jax/teacher_backends/`
- `pyBall/gridff_jax/pose_sampling/`
- `pyBall/gridff_jax/hybrid_energy/`
- `pyBall/gridff_jax/fit/`
- `pyBall/gridff_jax/export/`
- `pyBall/gridff_jax/validation/`
- `pyBall/gridff_jax/benchmarking.py`

Runner scripts:

- `tests/tGridFFJax/run_nacl_parity.py`
- `tests/tGridFFJax/build_ag_dataset.py`
- `tests/tGridFFJax/fit_hybrid_gridff.py`
- `tests/tGridFFJax/validate_hybrid_gridff.py`
- `tests/tGridFFJax/run_ag_probe_ladder.py`
- `tests/tGridFFJax/run_ablations.py`
- `tests/tGridFFJax/benchmark_pose_paths.py`
- `tests/tGridFFJax/benchmark_ag_zscan.py`
- `tests/tGridFFJax/run_ag_transfer_suite.py`
- `tests/tGridFFJax/run_ag_anchor_diagnostics.py`
- `tests/tGridFFJax/run_madsurf_teacher_sanity.py`
- `tests/tGridFFJax/run_full_suite.py`

## Local Assets Used Here

- Ag `CHGCAR`:
  `/home/niel/git/ORR_HER_Ag_Colab/results/Ag_ORR_HER/slab_clean/final_scf_12x12x1/CHGCAR`
- Ag `LOCPOT`:
  `/home/niel/git/ORR_HER_Ag_Colab/results/Ag_ORR_HER/slab_clean/workfunc_12x12x1/LOCPOT`
- primary MAD-SURF model:
  `tests/tGridFFJax/mad-surf_data/models/full_dataset_config_weights/MACE_model.model`
- full extxyz:
  `tests/tGridFFJax/mad-surf_data/full_train_test_std_config_types.extxyz`
- small extxyz:
  `tests/tGridFFJax/mad-surf_data/dataset/test_small_dataset_std.extxyz`

## What Data Comes From Where

For this workflow there are three distinct data sources:

1. substrate fields
   from the slab only, for example Ag `CHGCAR` and `LOCPOT`
2. teacher labels
   from MAD-SURF, evaluated on the exact rigid poses we generate
3. optional reference labels
   from the MAD-SURF extxyz bundle where `REF_energy` and `REF_forces` are stored

That means:

- MAD-SURF is not used to create the substrate `P/L/Q` grid directly
- MAD-SURF is used to label rigid adsorbate-on-surface poses
- the local extxyz files can be used to sanity-check MAD-SURF against DFT reference data when the chemistry overlaps

Concrete local evidence:

- `tests/tGridFFJax/mad-surf_data/dataset/test_small_dataset_std.extxyz`
  contains `REF_energy` and `REF_forces`
- `tests/tGridFFJax/mad-surf_data/full_train_test_std_config_types.extxyz`
  contains `REF_energy`, `REF_forces`, and `config_type`

Not every file in `mad-surf_data/` is a uniform reference-label benchmark.
Some figure folders contain structures, trajectories, or LAMMPS inputs used for the paper figures.

## Environment

For this machine, the recommended install path is:

```bash
python3 -m pip install --user --upgrade "jax[cuda12]" optax
python3 -m pip install --user --upgrade numpy scipy matplotlib ase mace-torch
```

Why `jax[cuda12]` and not `jax[cuda12-local]` here:

- local `nvcc` is CUDA `12.0`
- the wheel path is more robust on this setup
- JAX successfully sees `CudaDevice(id=0)` outside the sandbox with this install

Check JAX GPU visibility:

```bash
python3 -c "import jax; print(jax.devices()); print(jax.default_backend())"
```

If matplotlib complains about config paths:

```bash
export MPLCONFIGDIR=/tmp/mpl_gridff
```

## Important Implementation Notes

### B-spline path

The current JAX student now uses cubic B-spline sampling, not the earlier trilinear prototype.

- raw slab fields are built on the `zyx` grid
- cubic B-spline coefficients are prefiltered in Python
- JAX evaluates those coefficient grids with the same cubic basis used in the OpenCL GridFF kernels
- `Bspline_PLQd.npy` now exports coefficient-space grids, not raw field values

### Builder modes

Phase-1 strict `PLQ` now uses explicit substrate builder modes from `GridConfig.builder_mode`:

- `metal_density_plq` **(default, recommended for metals)**
  builds `P/L` from the DFT electron density `rho(r)^power` from `CHGCAR`;
  uses `LOCPOT` for `Q`. Physically correct for delocalized metal electrons.
- `parity_core`
  uses GridFF-style pairwise `REQ -> PLQ` substrate fields and direct periodic sums
- `metal_dft_plq`
  uses pairwise `P/L` Morse construction on atomic positions, prefers `LOCPOT` for `Q`;
  kept for backward compatibility
- `surrogate_density`
  older `rho^power` prototype normalized by `rho.max()`; fallback only

For NaCl parity work, use `surface_xyz + parity_core`.
For Ag/Au/Cu/Pt from real DFT volumetrics, use `vasp_volumetric + metal_density_plq`.

### Fixed charges in PLQ

The strict `PLQ` path now carries fixed adsorbate charges through the dataset and the student model.

Current built-in adsorbates:

- `H`: `[0.0]`
- `CO`: `[-0.0207, 0.0207]`
- `H2O`: `[-0.834, 0.417, 0.417]`
- `CHONH2`: `[0.616, -0.571, -0.862, 0.065, 0.376, 0.376]` (formamide, AMBER GAFF2 RESP)

`fit_static_charge` is off by default in the Ag strict-`PLQ` path. The goal in phase 1 is to keep the molecule-side `Q` term fixed and improve the substrate `PLQ` baseline first.

### Systematic Ag sampling

The Ag workflow no longer uses random orientation tied to the systematic `z` ladder.

- systematic site scans use deterministic orientation libraries
- `z` sampling is biased toward short range with `z_bias_power > 1`
- `H/C/O` probe fits can be used to initialize the `H/CO/H2O` molecule fit

This is the cleanest phase-1 route for diagnosing whether missing error comes from data coverage or from the `PLQ` model itself.

## Phase-2 Ag Protocol

The current recommended Ag protocol is:

- `builder_mode = metal_density_plq`
- `use_req_plq = True`
- fixed adsorbate charges
- rigid slab
- no `CT`, no image charge, no reactive residual
- fit on the `primary` molecular window `2.0-5.6 A`
- report but do not optimize directly on the `stress` window `1.8-2.0 A`

Each benchmark now reports three explicit windows:

- `primary`: `2.0-5.6 A`
- `stress`: `1.8-2.0 A`
- `full`: `1.8-5.6 A`

Training points are chosen deterministically per site:

- `18` points in the dense chemisorption window `2.0-2.8 A`
- `12` points in the outer window `2.8-5.6 A`

The optimizer now:

- fits `req_radius_offset` directly
- fits `req_energy_scale` in log space
- applies explicit quadratic regularization
- reports whether the final fit is `constraint_limited`

`constraint_limited = true` means at least one fitted `REQ` parameter finished within `5%` of its allowed range. That is treated as a useful warning, not as success.

## Focused CHON z-Scan

For step-by-step debugging of Ag short-range physics, use:

```bash
MPLCONFIGDIR=/tmp/mpl_gridff python3 tests/tGridFFJax/benchmark_ag_zscan.py \
  --device cpu \
  --out-dir /tmp/ag_transfer_phase2_real/CHONH2 \
  --adsorbate-name CHONH2 \
  --stress-z-min 1.8 \
  --z-min 2.0 \
  --z-max 5.6 \
  --eval-z-points 121 \
  --train-low-count 18 \
  --train-high-count 12 \
  --heldout-tilt-x-deg 20 \
  --heldout-tilt-y-deg 10 \
  --max-steps 600 \
  --force-weight 10.0 \
  --learning-rate 1.0e-2 \
  --no-fit-static-charge \
  --prefer-jax
```

What it does:

- loads the Ag slab from `CHGCAR` and `LOCPOT`
- uses the built-in `CHONH2` formamide-like adsorbate
- fits only on the primary `z` window
- evaluates on the full window plus one held-out tilted orientation
- writes train/holdout, primary/stress/full, and bound-distance diagnostics

Current real CHON result:

- output: `/tmp/ag_transfer_phase2_real/CHONH2`
- teacher: MAD-SURF on CPU
- student: JAX/Optax on `gpu:0`
- dense scan: `121` points per site, `363` poses total
- primary energy RMSE: `0.0867 eV`
- primary force RMSE: `0.2857 eV/A`
- stress energy RMSE: `0.4935 eV`
- stress force RMSE: `1.1661 eV/A`
- held-out tilted primary energy RMSE: `0.0787 eV`
- held-out tilted primary force RMSE: `0.3762 eV/A`

Important caveat:

- this CHON fit is still `constraint_limited`
- `req_radius_offset` for `C` and `H` finished at the upper bound `0.8 A`
- `req_energy_scale` for `H` finished at the upper bound `8.0`

Interpretation:

- the new windowed protocol is real and stable
- the primary-window CHON target is now strong
- transfer to a tilted orientation is also good
- but the molecule-side REQ fit is still compensating hard enough that the fit cannot yet be called physically relaxed

Important teacher-backend note:

- the focused CHON `z` scan still hangs on `--device cuda`
- the CPU teacher path completes cleanly
- the student JAX path still runs on `gpu:0` outside the sandbox

## Transfer Suite

To run the same protocol on `CHONH2`, `CO`, and `H2O` with one command:

```bash
MPLCONFIGDIR=/tmp/mpl_gridff python3 tests/tGridFFJax/run_ag_transfer_suite.py \
  --device cpu \
  --out-dir /tmp/ag_transfer_phase2_real \
  --eval-z-points 121 \
  --train-low-count 18 \
  --train-high-count 12 \
  --heldout-tilt-x-deg 20 \
  --heldout-tilt-y-deg 10 \
  --max-steps 600 \
  --force-weight 10.0 \
  --learning-rate 1.0e-2 \
  --prefer-jax
```

Current real transfer results:

- `CHONH2`
  primary energy `0.0867 eV`, primary force `0.2857 eV/A`
  held-out tilted primary energy `0.0787 eV`, held-out tilted primary force `0.3762 eV/A`
- `CO`
  primary energy `0.4867 eV`, primary force `0.2559 eV/A`
  held-out tilted primary energy `0.4695 eV`, held-out tilted primary force `0.2651 eV/A`
- `H2O`
  primary energy `0.0475 eV`, primary force `0.0925 eV/A`
  held-out tilted primary energy `0.0459 eV`, held-out tilted primary force `0.0923 eV/A`

What this means:

- `CHONH2` and `H2O` transfer well under the current strict-`PLQ` setup
- `CO` still has a clear energy residual even though its force is already quite good
- all three fits are still `constraint_limited`, so this is not yet the final physically relaxed parameterization

## Atomic Anchor Diagnostics

To localize remaining error by element channel:

```bash
MPLCONFIGDIR=/tmp/mpl_gridff python3 tests/tGridFFJax/run_ag_anchor_diagnostics.py \
  --device cpu \
  --out-dir /tmp/ag_anchor_phase2_real \
  --eval-z-points 121 \
  --train-low-count 18 \
  --train-high-count 12 \
  --max-steps 400 \
  --force-weight 10.0 \
  --learning-rate 1.0e-2 \
  --prefer-jax
```

Current real atomic-anchor results:

- `H`
  primary energy `0.4515 eV`, primary force `0.3047 eV/A`
- `C`
  primary energy `0.6434 eV`, primary force `0.5999 eV/A`
- `N`
  primary energy `0.4131 eV`, primary force `0.3657 eV/A`
- `O`
  primary energy `0.8656 eV`, primary force `0.6304 eV/A`

Interpretation:

- all four anchor fits are still `constraint_limited`
- `req_energy_scale` is pinned at `8.0` for `H`, `C`, `N`, and `O`
- `C` and `O` are the worst primary-window energy channels
- that is consistent with the current `CO` transfer problem being a genuine short-range substrate/REQ mismatch, not just a random CO-only fitting artifact

## Current Best Real Ag Strict-PLQ Result

Output directory:

- `/tmp/ag_strict_plq_forcefixed_fit`

Dataset:

- real Ag `CHGCAR` + `LOCPOT`
- MAD-SURF on `cuda`
- `654` poses total:
  `H=158`, `CO=248`, `H2O=248`
- systematic settings:
  `samples_per_site=10`,
  `systematic_orientations=4`,
  `random_samples=128`,
  `z_bias_power=2.5`

Validation against MAD-SURF interaction energies and corrected interaction forces:

- aggregate energy RMSE: `0.819 eV`
- aggregate force RMSE: `0.710 eV/A`
- `H`: energy `1.219 eV`, force `0.416 eV/A`
- `CO`: energy `0.828 eV`, force `0.697 eV/A`
- `H2O`: energy `0.372 eV`, force `0.766 eV/A`

Main plots:

- `/tmp/ag_strict_plq_forcefixed_fit/validation/CO_energy_parity.png`
- `/tmp/ag_strict_plq_forcefixed_fit/validation/CO_force_parity.png`
- `/tmp/ag_strict_plq_forcefixed_fit/validation/H2O_energy_parity.png`
- `/tmp/ag_strict_plq_forcefixed_fit/validation/H2O_force_parity.png`
- `/tmp/ag_strict_plq_forcefixed_fit/validation/training_loss_log.png`

## Deterministic Fixed-Path Benchmark

The mixed-pose dataset validation above is necessary, but it is not sufficient for judging smoothness.
The correct diagnostic is to compare MAD-SURF and GridFF on exactly the same fixed paths:

- fixed-site `z` lines
- fixed-height `x`, `y`, diagonal, and skew `xy` lines
- `xyz` skew lines
- yaw and tilt rotations for polyatomics
- fixed-height `xy` planes
- `xz` planes

Real output directories:

- `H`:
  `/tmp/ag_strict_plq_forcefixed_scans/H/path_benchmark.json`
- `CO`:
  `/tmp/ag_strict_plq_forcefixed_scans_chunked_v2/CO/path_benchmark.json`
- `H2O`:
  `/tmp/ag_strict_plq_forcefixed_scans_chunked_v2/H2O/path_benchmark.json`

Representative results:

- `H`
  average path RMSE:
  energy `0.597 eV`, force `0.387 eV/A`
  worst plane:
  `top_xy_plane_z2p2` with energy `1.030 eV`, force `0.727 eV/A`
- `CO`
  average path RMSE:
  energy `0.383 eV`, force `0.385 eV/A`
  worst plane:
  `top_xy_plane_z2p2` with energy `0.899 eV`, force `0.575 eV/A`
- `H2O`
  average path RMSE:
  energy `0.094 eV`, force `0.194 eV/A`
  worst path:
  `top_z_line` with energy `0.271 eV`, force `0.668 eV/A`

Interpretation:

- `H` and `CO` still show a mostly positive energy bias
  the current Ag strict-`PLQ` grid is too repulsive or not attractive enough
- the largest errors are concentrated in near-surface lateral corrugation
  especially the `xy` plane around `z = 2.2 A`
- `H2O` is already much closer to the teacher than `H` or `CO`
- GridFF is often smoother than MAD-SURF on fixed paths for `H` and `CO`
  which is acceptable if it removes MLIP short-range roughness without destroying the physically relevant corrugation
- for `H2O`, GridFF is not always smoother than MAD-SURF
  so smoothness alone is not the optimization target; the correct target is low error with physically sensible corrugation

The chunked path benchmark now also stores raw traces:

- `*_trace.npz`

These contain the exact path coordinate, MAD-SURF energies/forces, and GridFF energies/forces, so later analysis can be done without rerunning the teacher.

## How to Read the `z` Profiles

Be careful with apparent `wiggles` in `*_z_profile.png`.

In the Ag workflow, those plots are generated from mixed rigid poses:

- multiple deterministic orientations
- random orientations
- top / bridge / hollow registries
- small lateral jitter

So a raw line drawn through all points sorted only by `z` is not a true 1D scan. It mixes several different physical curves. The validation plotter now shows:

- sample scatter
- binned teacher mean
- binned student mean

If you want to judge whether the teacher itself is smooth, the correct diagnostic is a fixed-site, fixed-orientation path scan, not the mixed-pose dataset profile.

### No `try/except`

The new `pyBall/gridff_jax` and `tests/tGridFFJax` code was kept without `try/except` blocks, per project rule.

## Main Commands

### 1. Build a real Ag teacher dataset

Example: focused `H2O/Ag(111)` dataset for strict `PLQ`

```bash
env MPLCONFIGDIR=/tmp/mpl_gridff \
python3 tests/tGridFFJax/build_ag_dataset.py \
  --teacher madsurf \
  --device cuda \
  --model-path tests/tGridFFJax/mad-surf_data/models/full_dataset_config_weights/MACE_model.model \
  --adsorbates H2O \
  --samples-per-site 16 \
  --random-samples 96 \
  --representative-sites-per-label 1 \
  --out-dir /tmp/tGridFFJax_plq_h2o_dataset
```

For reduced-grid probe tests on the real Ag slab, add for example:

```bash
  --grid-shape 24,24,60
```

### 2. Fit the student in strict `PLQ` mode

```bash
env MPLCONFIGDIR=/tmp/mpl_gridff \
python3 tests/tGridFFJax/fit_hybrid_gridff.py \
  --config /tmp/tGridFFJax_plq_h2o_dataset/config.json \
  --dataset-dir /tmp/tGridFFJax_plq_h2o_dataset/datasets \
  --out-dir /tmp/tGridFFJax_plq_h2o_fit \
  --mode plq \
  --max-steps 600 \
  --force-weight 5.0 \
  --prefer-jax
```

To initialize from an `H/C/O` probe fit:

```bash
  --init-fit-dir /tmp/ag_probe_fit
```

### 3. Validate the fit

```bash
env MPLCONFIGDIR=/tmp/mpl_gridff \
python3 tests/tGridFFJax/validate_hybrid_gridff.py \
  --config /tmp/tGridFFJax_plq_h2o_fit/config_used.json \
  --dataset-dir /tmp/tGridFFJax_plq_h2o_dataset/datasets \
  --fit-dir /tmp/tGridFFJax_plq_h2o_fit \
  --out-dir /tmp/tGridFFJax_plq_h2o_fit/validation \
  --mode plq \
  --prefer-jax
```

### 4. Benchmark smooth dragging and rotation

```bash
env MPLCONFIGDIR=/tmp/mpl_gridff \
python3 tests/tGridFFJax/benchmark_pose_paths.py \
  --config /tmp/tGridFFJax_plq_h2o_fit/config_used.json \
  --fit-params /tmp/tGridFFJax_plq_h2o_fit/fit_params.json \
  --adsorbates H2O \
  --out-dir /tmp/tGridFFJax_plq_h2o_fit/path_benchmarks_warm \
  --prefer-jax
```

### 5. Run the strict-PLQ Ag curriculum

This is the recommended real benchmark path on this machine:

```bash
env MPLCONFIGDIR=/tmp/mpl_gridff \
python3 tests/tGridFFJax/run_ag_strict_plq_curriculum.py \
  --teacher-device cuda \
  --work-dir /tmp/ag_strict_plq_curriculum_real \
  --prefer-jax
```

It does:

1. `H/C/O` probe dataset build
2. `H/C/O` strict-`PLQ` fit
3. `H/CO/H2O` dataset build
4. `H/CO/H2O` fit initialized from the probe fit
5. validation plots and metrics

### 5. Build a NaCl parity-core grid

```bash
env MPLCONFIGDIR=/tmp/mpl_gridff \
python3 tests/tGridFFJax/run_nacl_parity.py \
  --xyz cpp/common_resources/xyz/NaCl_1x1_L3.xyz \
  --desired-voxel 0.1 \
  --out-dir /tmp/tGridFFJax_nacl_parity
```

If you already have a reference `Bspline_PLQd.npy` from the current OpenCL generator, add:

```bash
  --reference-grid /path/to/Bspline_PLQd.npy
```

### 6. Run the Ag probe ladder

This is the smallest real-DFT Ag entry point for phase-1 strict `PLQ`.

```bash
env MPLCONFIGDIR=/tmp/mpl_gridff \
python3 tests/tGridFFJax/run_ag_probe_ladder.py \
  --teacher-device cuda \
  --adsorbates H,C,O,CO,H2O \
  --grid-shape 24,24,60 \
  --prefer-jax
```

## What Was Actually Measured

### Earlier full hybrid Ag run

The broad `H/CO/H2O` hybrid run with `CT/image/reactive` enabled was not acceptable.

Aggregate validation from that run:

- energy RMSE: `2.821 eV`
- force RMSE: `344.823 eV/A`

The ablation result showed the present `CT/QEq` branch is destabilizing:

- `PLQ`: energy RMSE `0.711 eV`, force RMSE `1.936 eV/A`
- `full`: energy RMSE `2.821 eV`, force RMSE `344.925 eV/A`

Conclusion: do not trust the extra dynamic terms yet on Ag.

### Focused strict PLQ benchmark: `H2O/Ag(111)`

Real run using:

- teacher: MAD-SURF on `cuda`
- student: JAX + Optax on `gpu:0`
- slab: Ag from real `CHGCAR` and `LOCPOT`
- points: `120` total poses

Dataset generation:

- teacher total: `16.09 s`
- teacher per pose: `0.134 s`

Fit:

- backend: `optax_adam_jax`
- device: `gpu:0`
- total fit time: `22.93 s`
- best step count: `442`

Validation:

- all-sample energy RMSE: `1.111 eV`
- all-sample force RMSE: `2.656 eV/A`
- test-split energy RMSE: `0.083 eV`
- test-split force RMSE: `2.407 eV/A`

The important z-resolved result:

- `z = 4.0-5.5 A`: energy RMSE `0.002 eV`
- `z = 3.0-4.0 A`: energy RMSE `0.030 eV`
- `z = 2.5-3.0 A`: energy RMSE `0.073 eV`
- `z = 2.0-2.5 A`: energy RMSE `0.367 eV`
- `z = 1.5-2.0 A`: energy RMSE `3.022 eV`

This shows the current strict `PLQ` path is already reasonable in the far and mid range, but fails at short range near the surface.

### Interpolation self-consistency

The cubic B-spline representation itself is not the main failure.

Sparse voxel-center reproduction errors for the exported/sampled fields:

- Pauli RMSE: `3.52e-4`
- London RMSE: `1.09e-3`
- Coulomb RMSE: `1.85e-1 eV`

For the Coulomb field, the grid standard deviation is about `10.19 eV`, so the interpolation error is small compared with the field variation.

Conclusion: the B-spline interpolation/export path is now acceptable. The remaining mismatch is mostly in the physical field construction and fitting model, not in the interpolator.

### Steady-state path timing on `H2O`

After excluding the first-call JAX compile cost, the student is much faster than the teacher:

- `z_drag`: about `3521x` faster
- `xy_drag`: about `5307x` faster
- `yaw_rotation`: about `4318x` faster
- `tilt_rotation`: about `2144x` faster

These are steady-state timings. If you include the first JIT compile call, the first path can look slower than MAD-SURF. That is expected and should not be used as the runtime comparison.

## What Is Still Wrong

The current main error is not “too few points”.

The current main error is that the substrate `P/L` fields for metals are still too crude.

Right now the metal fields are built as:

- `Pauli ~ rho^pauli_power`
- `London ~ -rho^london_power`
- `Coulomb ~ LOCPOT or Poisson(rho)`

That is a workable prototype, but it is not yet the same physical decomposition used in the original NaCl GridFF workflow.

Evidence:

- in the focused `H2O/Ag` fit, the optimizer pushed `Pauli` to the upper bound `5.0`
- it pushed `London` to the lower bound `0.05`
- loss decreased only weakly over hundreds of steps

That means the student is trying to compensate for the wrong field shapes with coefficient extremes.

So the next accuracy step is not “just add more points”.

The next accuracy step is:

1. keep strict `PLQ`
2. keep the cubic B-spline path
3. improve how `P` and `L` are constructed from slab electronic structure
4. only after that revisit `CT/image/reactive`

## Practical Meaning for FireCore

If you stop at strict `PLQ`, then yes:

- you can generate a `Bspline_PLQd.npy`
- FireCore can load that static grid
- you can do rigid-surface dragging, scans, and MD-like manipulations using the grid

What you still need alongside it:

- the molecule’s own force field
- molecule-molecule interactions if more than one adsorbate is present
- any special reactive or charge-transfer physics not encoded in the static grid

## Recommended Next Step

Do not expand to more chemistry yet.

First make `PLQ` physically correct on one clean Ag benchmark, for example:

- `H2O/Ag(111)` for electrostatics and orientation
- optionally `H/Ag(111)` as a pure short-range check

Then improve the metal `P/L` construction until:

- far-range energy is already good
- mid-range force errors drop
- short-range z-scan no longer shows the large low-z collapse

Only then turn `CT/image/reactive` back on.

---

## What Was Changed and How to Use It (2026-03-17)

### Summary of Changes

Four improvements were made to fix the Ag GridFF accuracy and extend to new molecules:

1. **New `metal_density_plq` substrate builder** in `pyBall/gridff_jax/substrate_fields.py`
   Builds Pauli and London grids from the DFT electron density `rho(r)` from `CHGCAR`
   (previously: pairwise Morse on atomic positions, identical to the NaCl ionic model).
   Coulomb stays on LOCPOT as before.
   This is now the default `builder_mode`.

2. **Config extensions** in `pyBall/gridff_jax/config.py`
   New `GridConfig` parameters:
   - `metal_density_pauli_power = 1.0`  power for `P = rho_norm^power`
   - `metal_density_london_power = 0.5`  power for `L = -rho_norm^power`
   - `metal_density_rho_smoothing_sigma = 0.0`  Gaussian pre-smooth in voxel units
   - `metal_bulk_electron_density = 0.0530`  e/Å³, normalizer (Ag default; change for Au/Cu/Pt)

   New `FeatureToggles` parameter:
   - `use_image_charge_fixed = False`  standalone image charge, no QEq required

3. **Image charge bug fix** in `pyBall/gridff_jax/hybrid_energy/model.py`
   `e_image` was computed but never added to the total energy (`total = ... + e_image` was missing).
   Fixed in both the NumPy and JAX paths.
   Also added a standalone image charge path (`use_image_charge_fixed`) that works
   without enabling the broken QEq solver.

4. **CHONH2 (formamide) built-in adsorbate** in `pyBall/gridff_jax/pose_sampling/molecules.py`
   Available as `"CHONH2"` in `benchmark_adsorbates()` with standard B3LYP geometry
   and AMBER GAFF2 RESP partial charges (sum = 0.000 e, 6 atoms: C, O, N, H, H, H).

### Test Results: CHONH2/Ag(111), metal_density_plq

Run on 2026-03-17 using MAD-SURF teacher (CPU), 183 poses total (top + bridge + hollow z-scan,
z from 1.8 to 5.5 Å, 61 eval points per site):

```
builder_mode = metal_density_plq (NEW, default)
  all-sample energy RMSE : 0.288 eV
  all-sample force  RMSE : 0.462 eV/Å (components)
  holdout energy RMSE    : 0.246 eV
  far-z (z > 3.0 Å) energy RMSE : 0.049 eV   <- good
  mid-z (2.0 - 3.0 Å) energy RMSE : 0.203 eV
  low-z (z < 2.0 Å) energy RMSE   : 1.039 eV  <- short-range still poor
  pauli coefficients     : NOT saturated at bounds (all near 1.0)
  london coefficients    : NOT saturated at bounds (all near 1.0)
```

The key validation: with `use_req_plq = True` (default in all scripts), pauli and london
scale factors are no longer pushed to the optimizer bounds (5.0 / 0.05), confirming the
new substrate field shapes are physically compatible with the optimizer.

Short-range errors at `z < 2.0 Å` remain because the molecular exchange-repulsion profile
is still not fully matched by the current REQ parametrization (see `req_energy_scale`
hitting its upper bound). This is the next target, separate from the substrate field fix.

### How to Use the New Builder

The new `metal_density_plq` mode is the default. No changes to existing command lines are needed.
The only requirement is that `CHGCAR` is provided (not empty):

```bash
env MPLCONFIGDIR=/tmp/mpl_gridff \
python3 tests/tGridFFJax/calibrate_ag_zscan_plq.py \
  --device cpu \
  --adsorbate-name CHONH2 \
  --xyz-path cpp/common_resources/xyz/CHONH2.xyz \
  --use-input-charges \
  --sites top,bridge,hollow \
  --stress-z-min 1.8 --z-min 2.0 --z-max 5.6 \
  --eval-z-points 121 \
  --train-low-count 18 \
  --train-high-count 12 \
  --max-steps 300 \
  --force-weight 10.0 \
  --prefer-jax \
  --out-dir /tmp/ag_chonh2_metal_density
```

To use a different molecule from the built-in library (`H`, `CO`, `H2O`, `CHONH2`),
pass it with `--adsorbate-name` (and omit `--xyz-path` to use the built-in geometry):

```bash
  --adsorbate-name H2O
```

To use `CHONH2` with DFT-optimized geometry from the xyz file (recommended):

```bash
  --adsorbate-name CHONH2 \
  --xyz-path cpp/common_resources/xyz/CHONH2.xyz \
  --use-input-charges
```

### How to Tune the Builder for Other Metals

For Au, Cu, or Pt, provide the correct `CHGCAR`/`LOCPOT` paths and change the
bulk electron density normalizer in the config JSON or by overriding in the script:

```python
config.grid.metal_bulk_electron_density = 0.0590  # Au
config.grid.metal_bulk_electron_density = 0.0845  # Cu
config.grid.metal_bulk_electron_density = 0.0660  # Pt
```

To tune the density-to-field mapping, sweep `metal_density_pauli_power` and
`metal_density_london_power` (values between 0.5 and 1.5):

```python
config.grid.metal_density_pauli_power  = 1.0   # default; lower = softer repulsion
config.grid.metal_density_london_power = 0.5   # default; lower = broader attraction
```

To apply Gaussian pre-smoothing to the density before raising to the power (useful
if CHGCAR has grid noise near the vacuum boundary):

```python
config.grid.metal_density_rho_smoothing_sigma = 1.0  # in voxel units
```

### How to Enable Standalone Image Charge

Without QEq (stable), add to the config before fitting:

```python
config.toggles.use_image_charge_fixed = True
config.training.fit_image_plane       = True   # optional: let optimizer find z_image
```

This computes `E_image = -image_scale * sum(Q_i^2 / (4*(z_i - z_image)))` analytically.
Useful for polar molecules like H2O and CHONH2 where the metal image charge is ~0.05-0.1 eV.
Does nothing for neutral non-polar adsorbates (net charge contribution near zero).
