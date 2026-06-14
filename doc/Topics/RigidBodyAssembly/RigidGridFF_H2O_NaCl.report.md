---
description: Rigid GridFF H2O-on-NaCl and PTCDA-on-NaCl multi-replica relaxation report and tutorial
---

# Rigid GridFF — H2O on NaCl and PTCDA on NaCl (multi-replica) — Report and Tutorial

This document consolidates two related efforts:

1) **Single-molecule H2O on NaCl fast relaxation** (baseline rigid runner, convergence metrics, restart).
2) **Production multi-replica PTCDA on NaCl** with GridFF, GPU batching, stability diagnostics, and clustering of distinct minima.

It serves both as a progress report (what changed, what worked, what failed, and how we fixed it) and as a tutorial for running and interpreting the workflows.

## 1) Objectives
- Achieve **fast, stable rigid-body relaxation** on GridFF surfaces (H2O baseline, PTCDA production).
- Provide **explicit, physics-based convergence metrics** (body force/torque norms) and **lightweight stability diagnostics** (jump/height checks).
- Support **restartable two-stage workflows** (fast drop → settle/refine) without CPU/GPU thrash.
- Enable **massively parallel multi-replica runs** (PTCDA ~5k poses feasible) with minimal syncs and high throughput.
- Export **full XYZ stacks** (initial/final/converged/cluster reps) with metadata, and **cluster minima** by full 6D pose (translation with PBC + quaternion).

## 2) Implemented changes (code)

### Core rigid wrapper
- `pyBall/OCL/RigidBodyDynamics.py`
  - Added host-side storage of atom-body data, inertia, total mass; enabled **`reset_pose()`** for multi-stage restarts.
  - Added **`download_selected`** to cut GPU→CPU sync cost when only a subset of buffers is needed mid-loop.
  - Added **separate effective masses**: `mass_trans` (pos.w) and `mass_rot` (scales inverse inertia as `Iinv_relax = Iinv * (mtot/mass_rot)`).

### H2O runner (single-body baseline)
- `tests/tMMFF/test_rigid_gridff_surface.py`
  - Convergence metrics and plotting (`convergence.png`, semilogy |F|, |T|; configurable thresholds; early stop).
  - `--init-from-npz` to restart for settle passes; trajectory/plots remain.
  - Type-map preserved (default `O:O_3,H:H_OH`); z-scan option remains.

### PTCDA multi-replica runner (production)
- `tests/tMMFF/test_rigid_gridff_ptcda_batch.py`
  - Batch tiling of translations/rotations; restart from NPZ (`--init-from-npz`).
  - **Sampled stability diagnostics**: sampled z-min check and jump detection (`sample_check_stride`, `sample_every_steps`, `zmin_limit`, `jump_limit`). Flags `unstable` without heavy sync.
  - **Convergence criteria**: body force and torque norms compared to thresholds (`fconv`, `tconv`), using physical units (no hidden scaling).
  - **Clustering bug fixed**: clustering now uses **final relaxed fractional translation** + **final quaternion** (previously mixed initial translation with final rotation). PBC-aware in-plane distance, quaternion distance threshold.
  - **XYZ exports**: initial, final, converged, cluster representatives. Comments carry energy, |F|, |T|, converged flag, unstable flag, rot_id, initial frac (`frac0`), final frac (`fracf`).
  - **Outputs logged**: converged indices, unstable indices, NPZ with forces/torques/energies, cluster summary.

## 3) Workflows and commands

### A) H2O on NaCl (single body, fast drop + settle)
1) Fast-drop scan (example):
```bash
python3 test_rigid_gridff_surface.py --nsteps 1500 --stop-on-conv \
  --report-every 250 --dump-every 250 --plot-every 50 --conv-every 10 \
  --outdir output_rigid_gridff_fastscan_e \
  --type-map O:O_3,H:H_OH --x0 0.0 --y0 0.0 --z0 2.0 \
  --rx 0 --ry 35 --rz 30 --dt 0.08 --lin-damp 0.88 --ang-damp 0.80 \
  --force-scale 1.8 --torque-scale 1.8
```

2) Settle from saved pose:
```bash
python3 test_rigid_gridff_surface.py \
  --init-from-npz output_rigid_gridff_fastscan_e/rigid_scan_final.npz \
  --nsteps 1200 --stop-on-conv --report-every 100 --dump-every 100 \
  --plot-every 25 --conv-every 10 --outdir output_rigid_gridff_settle_e \
  --dt 0.03 --lin-damp 0.75 --ang-damp 0.70
```

### B) PTCDA on NaCl (multi-replica production, two-stage)
Parameters that proved stable and physical:
- `dt=0.02`, `lin_damp=0.99`, `ang_damp=0.97`, `mass_trans=1.0`, `mass_rot=4.0`, `force_scale=torque_scale=1.0`, `z0=5.2`.

1) Stage 1: systematic grid x,y and rotations (example 20×20×8 = 3200 replicas):
```bash
python3 test_rigid_gridff_ptcda_batch.py --nx 20 --ny 20 --nrot 8 --max-bodies 3200 \
  --nsteps 16000 --chunk-steps 400 --report-every-steps 4000 \
  --sample-check-stride 64 --sample-every-steps 400 \
  --z0 5.2 --dt 0.02 --lin-damp 0.99 --ang-damp 0.97 \
  --mass-trans 1.0 --mass-rot 4.0 --force-scale 1.0 --torque-scale 1.0 \
  --outdir output_rigid_gridff_ptcda_prod_stage1
```

2) Stage 2: refinement from stage-1 NPZ (same settings):
```bash
python3 test_rigid_gridff_ptcda_batch.py \
  --init-from-npz output_rigid_gridff_ptcda_prod_stage1/rigid_ptcda_batch_final.npz \
  --nsteps 20000 --chunk-steps 400 --report-every-steps 4000 \
  --sample-check-stride 64 --sample-every-steps 400 \
  --dt 0.02 --lin-damp 0.99 --ang-damp 0.97 \
  --mass-trans 1.0 --mass-rot 4.0 --force-scale 1.0 --torque-scale 1.0 \
  --outdir output_rigid_gridff_ptcda_prod_stage2_refine
```

## 4) Stability, convergence, and “physically meaningful” criteria

- **Convergence test (per replica):** `|F_body| < fconv` AND `|T_body| < tconv` (defaults 1e-3). Forces/torques are physical (no hidden scaling) — we divide by `force_scale/torque_scale` before testing.
- **Stability diagnostics (sampled):** every `sample_every_steps`, we download atoms (stride `sample_check_stride`) and flag:
  - `zmin < zmin_limit` (default 2.0 Å) → unstable
  - `jump > jump_limit` (default 1.0 Å between samples) → unstable
- **Physically meaningful geometry cross-check:** after transient, stable runs show z-height ≈ 3.1265 Å for PTCDA; sampled jumps drop to zero; converged poses lie flat on the surface.

## 5) Clustering of distinct minima (PTCDA)
- Input set: only **converged** replicas (and not flagged unstable).
- Pose distance:
  - Translation: final fractional coords with 2D PBC; in-plane metric via lattice vectors; threshold `cluster_pos_eps` (default 0.1 Å).
  - Rotation: quaternion distance with threshold `cluster_quat_eps` (default 0.02).
- Algorithm: energy-sorted greedy assignment — lowest-energy replica seeds a cluster; subsequent replicas join if within both thresholds; otherwise start a new cluster. Counts and representatives are written to `cluster_summary.txt`; representative atoms to `cluster_representatives.xyz`.
- Important fix: clustering now uses **final relaxed translation** (fracf), not the initial translation.

## 6) Results summary

### H2O baseline (rigid, current typing)
- Fast-drop then settle plateaus at residuals O(1e-2); reaching 1e-3 likely needs local pose search or flexible model.

### PTCDA production (20×20×8 = 3200 replicas)
- Stage 1 (stable regime): 980 / 3200 converged; 8 distinct minima.
- Stage 2 (refine from stage 1): 1154 / 3200 converged; 6 distinct minima.
- Stable diagnostics: no post-transient jumps; zmin tight at 3.1265 Å.
- Dominant minima:
  - Deep adsorbed flat state: E ≈ -0.885 eV, z ≈ 3.1265 Å (largest cluster, hundreds of members).
  - Higher adsorbed rotated states: E ≈ -0.699 eV, z ≈ 3.264–3.265 Å.

## 7) Outputs and how to read them (PTCDA)
- `initial_all.xyz` — all replicas at start; comment has `rot_id`, `frac0`.
- `final_all.xyz` — all replicas at end; comment fields:
  - `E`, `F`, `T`, `converged` (0/1), `unstable` (0/1), `rot_id`, `frac0` (initial fractional xy), `fracf` (final fractional xy).
- `converged_all.xyz` — only converged replicas; same metadata minus the flags.
- `cluster_representatives.xyz` — representative pose for each cluster; comment `cluster <id> count <n> E F T`.
- `cluster_summary.txt` — text summary: number of bodies, number converged, number of clusters, thresholds used, and for each cluster: representative replica id, count, E/F/T, position, quaternion.
- `converged_indices.txt` — per converged replica: index, |F|, |T|, energy, frac0, rot_id.
- `unstable_indices.txt` — replicas flagged by jump/zmin diagnostics.
- `rigid_ptcda_batch_final.npz` — restartable state (pos, quats, atom_positions, forces/torques, energies, flags, frac0, rot_id, reached step for force/torque thresholds).

## 8) Challenges and fixes
- **Rotational relaxation instability**: initial inverse inertia scaling was inverted. Fixed by scaling `Iinv` by `(mtot / mass_rot)`, and exposing `mass_rot` separately from translation.
- **Force/torque overscaling artifacts**: removed implicit scaling; use physical forces/torques with explicit masses/damping.
- **Clustering bug**: previously mixed initial translation with final rotation → unphysical clustering. Fixed to use final fractional translation + final quaternion with PBC-aware in-plane metric.
- **Detecting hidden blow-ups without heavy sync**: added sampled `zmin`/`jump` checks that touch atoms only at chosen intervals, marking `unstable` replicas.

## 9) Takeaways for future development
- Keep translational vs rotational effective masses separable; tune `mass_rot` and `ang_damp` for torque decay.
- Use sampled diagnostics for large batches; full downloads only when needed.
- Cluster in full 6D pose with PBC-aware translation; never mix initial and final pose components.
- For rigid systems plateauing at ~1e-2, consider local pose refinement (6D grid) or limited flexibility to reach 1e-3.

## 10) File map
- This report: `doc/Topics/RigidBodyAssembly/RigidGridFF_H2O_NaCl_report.md`
- Rigid wrapper: `pyBall/OCL/RigidBodyDynamics.py`
- H2O runner: `tests/tMMFF/test_rigid_gridff_surface.py`
- PTCDA multi-replica runner: `tests/tMMFF/test_rigid_gridff_ptcda_batch.py`
- Example outputs:
  - H2O: `tests/tMMFF/output_rigid_gridff_fastscan_*`, `output_rigid_gridff_settle_*`
  - PTCDA stage 1: `tests/tMMFF/output_rigid_gridff_ptcda_prod_stage1/`
  - PTCDA stage 2 refine: `tests/tMMFF/output_rigid_gridff_ptcda_prod_stage2_refine/`

---

# Rigid GridFF — PTCDA on CaF2 (multi-replica) — Report and Tutorial

This mirrors the NaCl workflow but uses the **validated CaF2 grid** built by the CaF2 tutorial scripts. Focus: validate the potential, locate the adsorption window, run two-stage batches, and filter true adsorbed minima.

## A) Goals (CaF2)
- Validate CaF2 GridFF field (signs/minima) with ion z-scans before any PTCDA batches.
- Run rigid PTCDA adsorption z-scan to locate the adsorption height window.
- Execute two-stage multi-replica relaxation on the validated CaF2 grid (separate outputs from NaCl).
- Separate physically adsorbed minima from shallow high-z stationary states.

## B) Sources / inspiration
- `tests/tMMFF/GridFF_CaF2_doc_tutorial.md`
- `tests/tMMFF/run_test_GridFF_CaF2.py`
- `tests/tMMFF/test_rigid_gridff_ptcda_batch.py`

## C) Challenges and fixes
- **Ca typing mismatch**: slab uses `Ca`, table uses `Ca+2`. Fixed by parsing/propagating `--type-map Ca:Ca+2` in batch runner and `RigidBodyDynamics.from_xyz_and_grid`.
- **Zero forces / huge z**: caused by wrong grid/typing. Fixed by using validated CaF2 grid `Bspline_PLQd_sigma_0p500.npy` and correct type-map.
- **High-z converged states**: CaF2 yields shallow stationary points far from surface; apply energy/z filter to isolate adsorbed minima.

## D) What we ran (and outputs)
- Validated grid (sigma=0.5): `tests/tMMFF/data/CaF2_6L_Ni3_rect_nx2_nz1_L2_top/Bspline_PLQd_sigma_0p500.npy`
- Substrate xyz: `cpp/common_resources/Substrates/generated_rect/CaF2_6L_Ni3_rect_nx2_nz1_L2_top.xyz`
- PTCDA xyz: `cpp/common_resources/xyz/PTCDA.xyz`
- Medium batch (6x6x8=288)
  - Stage 1: `tests/tMMFF/output_rigid_gridff_ptcda_caf2_midcheck/`
  - Stage 2 refine: `tests/tMMFF/output_rigid_gridff_ptcda_caf2_midcheck_refine/`
  - Converged 140/288; deep adsorbed subset 127 with E < -2.5 eV, z ≈ 6.97–7.00 Å
- Production (12x12x8=1152)
  - Stage 1: `tests/tMMFF/output_rigid_gridff_ptcda_caf2_prod_stage1/`
  - Stage 2 refine: `tests/tMMFF/output_rigid_gridff_ptcda_caf2_prod_stage2_refine/`
  - Converged 448/1152; distinct minima 222
  - Adsorbed subset (recommended filter): 344 with **E < -2.5 eV** and **z < 8.0 Å** (all z ≈ 6.965–6.999 Å)

## E) Run tutorial (CaF2)

### 1) Build/validate CaF2 grid + ion z-scans (sigma=0.5)
From `tests/tMMFF/`:
```bash
bash run_test_GridFF_CaF2.sh \
  --sigma 0.5 --sigma_scan 0.5 \
  --probe_specs H:0.2,O:-0.4 \
  --a 2.0 --auto_total_center 0.0 --auto_comp_center 0.0
```
Outputs: `tests/tMMFF/data/CaF2_6L_Ni3_rect_nx2_nz1_L2_top/` (plots, z-scan NPZs, `Bspline_PLQd_sigma_0p500.npy`).

### 2) Rigid PTCDA z-scan on validated CaF2 grid
Expectation: adsorption window COM z ~7.5 Å (surface top z ~4.88 Å). Example driver (simplified):
```python
# inputs: mol, sub, grid as above; type_map = {'Ca':'Ca+2'}
# scan_z = linspace(5.0, 20.0, 121)
# for each top Ca/F site: reset pose at (x,y,z), run_gridff once, sum atom_force[...,3] as energy
```

### 3) PTCDA CaF2 batch — Stage 1 (production tiling)
```bash
python3 test_rigid_gridff_ptcda_batch.py \
  --mol /home/prokop/git/FireCore/cpp/common_resources/xyz/PTCDA.xyz \
  --sub /home/prokop/git/FireCore/cpp/common_resources/Substrates/generated_rect/CaF2_6L_Ni3_rect_nx2_nz1_L2_top.xyz \
  --grid /home/prokop/git/FireCore/tests/tMMFF/data/CaF2_6L_Ni3_rect_nx2_nz1_L2_top/Bspline_PLQd_sigma_0p500.npy \
  --type-map Ca:Ca+2 \
  --nx 12 --ny 12 --nrot 8 --max-bodies 1152 \
  --nsteps 16000 --chunk-steps 400 --report-every-steps 4000 \
  --sample-check-stride 64 --sample-every-steps 400 \
  --z0 7.8 --dt 0.02 --lin-damp 0.99 --ang-damp 0.97 \
  --mass-trans 1.0 --mass-rot 4.0 --force-scale 1.0 --torque-scale 1.0 \
  --outdir output_rigid_gridff_ptcda_caf2_prod_stage1
```

### 4) PTCDA CaF2 batch — Stage 2 refinement (restart)
```bash
python3 test_rigid_gridff_ptcda_batch.py \
  --mol /home/prokop/git/FireCore/cpp/common_resources/xyz/PTCDA.xyz \
  --sub /home/prokop/git/FireCore/cpp/common_resources/Substrates/generated_rect/CaF2_6L_Ni3_rect_nx2_nz1_L2_top.xyz \
  --grid /home/prokop/git/FireCore/tests/tMMFF/data/CaF2_6L_Ni3_rect_nx2_nz1_L2_top/Bspline_PLQd_sigma_0p500.npy \
  --type-map Ca:Ca+2 \
  --init-from-npz output_rigid_gridff_ptcda_caf2_prod_stage1/rigid_ptcda_batch_final.npz \
  --nsteps 20000 --chunk-steps 400 --report-every-steps 4000 \
  --sample-check-stride 64 --sample-every-steps 400 \
  --dt 0.02 --lin-damp 0.99 --ang-damp 0.97 \
  --mass-trans 1.0 --mass-rot 4.0 --force-scale 1.0 --torque-scale 1.0 \
  --outdir output_rigid_gridff_ptcda_caf2_prod_stage2_refine
```

## F) Results (CaF2)
- Stage 1 (1152 reps): converged 295 / 1152; zmin stabilized at 6.933–6.934 Å; no post-transient jumps.
- Stage 2 refine: converged 448 / 1152; distinct minima 222.
- **Adsorbed subset (physically meaningful)**: 344 with E < -2.5 eV, z < 8.0 Å (all at z ≈ 6.965–6.999 Å).
- Remaining converged states are shallow/high-z stationary points; keep them separate from the adsorbed basin.

## G) Practical filters for CaF2 outputs
- `converged == 1`
- `unstable == 0`
- `energy < -2.5` eV
- `pos_z < 8.0` Å

## H) Takeaways (CaF2)
- Always validate GridFF on the surface (ion z-scans) before molecular batches.
- Type-map is essential: `Ca:Ca+2` for this slab/grid pair.
- CaF2 produces shallow high-z stationary states; apply an energy/z filter to isolate true adsorbed minima.
- Stable relaxation uses the same physical parameters as NaCl (dt 0.02, lin_damp 0.99, ang_damp 0.97, mass_rot 4.0, no force/torque overscaling) with higher starting height (z0 ≈ 7.8) matching the CaF2 adsorption window.

