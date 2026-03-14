---
description: Rigid GridFF H2O-on-NaCl fast relaxation session report
---

# Rigid GridFF H2O on NaCl — Fast Convergence Session Report

## 1) Objectives
- Achieve **fast, stable rigid-body relaxation** of an H2O molecule on the NaCl(001) GridFF surface.
- Add **explicit convergence metrics** (net body force/torque norms, decay history, threshold hits, early stop) to the rigid runner.
- Enable **reusable restart** from a saved pose for a two-stage workflow (fast drop → settle).
- Search integration parameters (dt, damping, force/torque scales) for rapid convergence well below 2000 steps.

## 2) Implemented changes (code)
- `pyBall/OCL/RigidBodyDynamics.py`
  - Stored host-side state (atom body positions, inertia, total mass) and added **`reset_pose()`** to restart a rigid pose with zeroed (or provided) momenta.
- `tests/tMMFF/test_rigid_gridff_surface.py`
  - Added **convergence metrics** collection and `convergence.png` plotting (|F_body|, |T_body|; semilogy).
  - Added **`--stop-on-conv`**, force/torque thresholds, and final summary of steps-to-threshold.
  - Added **`--init-from-npz`** to restart from a saved rigid pose (enables two-stage drop+settle).
  - Existing diagnostics retained: trajectory plot, per-step prints, z-scan option, type-map (`O:O_3,H:H_OH`).

## 3) Workflow tested
1) **Fast-drop scans** (vary dt, damping, scales, initial pose) to reach significant displacement/rotation quickly.
   - Example commands (from `tests/tMMFF/`):
     ```bash
     python3 test_rigid_gridff_surface.py --nsteps 1500 --stop-on-conv \
       --report-every 250 --dump-every 250 --plot-every 50 --conv-every 10 \
       --outdir output_rigid_gridff_fastscan_e \
       --type-map O:O_3,H:H_OH --x0 0.0 --y0 0.0 --z0 2.0 \
       --rx 0 --ry 35 --rz 30 --dt 0.08 --lin-damp 0.88 --ang-damp 0.80 \
       --force-scale 1.8 --torque-scale 1.8
     ```

2) **Settle phase** from best pose (uses `--init-from-npz` saved from step 1):
   ```bash
   python3 test_rigid_gridff_surface.py \
     --init-from-npz output_rigid_gridff_fastscan_e/rigid_scan_final.npz \
     --nsteps 1200 --stop-on-conv --report-every 100 --dump-every 100 \
     --plot-every 25 --conv-every 10 --outdir output_rigid_gridff_settle_e \
     --dt 0.03 --lin-damp 0.75 --ang-damp 0.70
   ```

## 4) Results observed
- **Fast-drop success (significant motion)**
  - Run `fastscan_e`: displacement ≈ **1.34 Å**, rotation ≈ **122°**, residuals **|F| ≈ 2.1e-2**, **|T| ≈ 9.1e-3** after 1500 steps.
  - Run `fastscan_d`: displacement ≈ **1.67 Å**, rotation ≈ **67°**, residuals **|F| ≈ 3.1e-2**, **|T| ≈ 2.8e-2**.
- **Settle-phase restart near adsorbed pose**
  - Run `settle_e` (restart from `fastscan_e` final pose): minimal further motion (≤0.014 Å), residuals plateau at **|F| ≈ 1.34e-2**, **|T| ≈ 7.4e-3**.
- **Key takeaway:** With current rigid GridFF + water typing, the apparent adsorbed pose plateaus around **1e-2** residual; simply reducing timestep/damping does **not** push to 1e-3. Likely needs a local pose search (x/y/z + small rotations) or relaxed flexibility to find a true force minimum.

## 5) How to use the updated runner
- **Location:** `tests/tMMFF/test_rigid_gridff_surface.py`
- **Common flags:**
  - `--nsteps`, `--dt`, `--lin-damp`, `--ang-damp`, `--force-scale`, `--torque-scale`
  - `--fconv`, `--tconv`, `--stop-on-conv`, `--conv-every`
  - `--init-from-npz` (restart from saved pose), `--type-map` (default `O:O_3,H:H_OH`)
  - `--report-every`, `--dump-every`, `--plot-every` (trajectory/convergence plots)
  - `--do-zscan` for height diagnostics (Na/Cl sites)
- **Outputs:**
  - `rigid_h2o_relax.xyz` (trajectory), `rigid_scan_final.npz` (restartable state), `trajectory.png`, `convergence.png`, final pose XYZ/NPZ.

## 6) Next recommended step
- Run a **local 6D pose refinement** around the adsorbed state (grid over x,y,z and small rx,ry,rz) to locate a lower-force basin; then a short settle pass. This will determine if `1e-3` residual is physically reachable with the rigid model.

## 7) File map
- Report (this file): `doc/Topics/RigidBodyAssembly/RigidGridFF_H2O_NaCl_report.md`
- Runner: `tests/tMMFF/test_rigid_gridff_surface.py`
- Rigid wrapper: `pyBall/OCL/RigidBodyDynamics.py`
- Example outputs: `tests/tMMFF/output_rigid_gridff_fastscan_*`, `output_rigid_gridff_settle_*`
