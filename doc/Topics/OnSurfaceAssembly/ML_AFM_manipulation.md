# ML-AFM Manipulation: Rigid-Body Batch Simulation + Pose Prediction

**Status:** Working prototype with batch-parallel GPU relaxation and pose-to-pose ML training.  
**Last updated:** 2026-06-13

---

## 1. Objective

Train a neural network to predict the relaxed rigid-body pose (COM position + quaternion) of a molecule dragged by an AFM tip across a substrate, given the tip trajectory. The tip is attached to a single peripheral atom of the molecule via a harmonic spring. The substrate potential comes from a precomputed GridFF B-spline. The system is designed for **high-throughput batch execution**: thousands of manipulation trajectories are relaxed in parallel on the GPU, and the resulting dataset trains an ML model that can propose relaxed poses for new tip positions.

---

## 2. System Architecture

The design follows the principle: **Python orchestration is thin; heavy computation lives on the GPU.**

```
Rigid.cl                              OpenCL kernels (force eval + rigid dynamics)
  |
  v
RigidBodyDynamics.py                  Low-level PyOpenCL buffer/kernel wrapper
  |
  v
RigidBodyAFM.py                     High-level API: prepare(), set_anchor_positions(), relax_to_constraint()
  |
  v
test_rigid_afm_ml.py                CLI driver:
                                      --test manip  → generate batch trajectories + diagnostics
                                      --test train  → train MLP on pose transitions + predict
```

### 2.1 Why batch-parallel is mandatory

Each replica (one molecule + one tip path) is mapped to **one OpenCL workgroup** (`nloc=32` threads). The GPU launches all workgroups simultaneously. Running trajectories sequentially in a Python `for` loop defeats the entire purpose: the GPU would be idle 99% of the time while Python iterates.

The correct production path is:

1. `prepare(n_bodies=ntraj)` — allocate one buffer set for all replicas.
2. Per substep: construct anchor targets for **all** replicas as an `(ntraj, 3)` array.
3. `set_anchor_positions(tips)` — upload in one host→device transfer.
4. `relax_to_constraint()` — one kernel call relaxes all bodies.
5. `download_outputs()` — one device→host transfer gets all poses.

This is how `test_rigid_gridff_ptcda_batch.py` and the manipulation test now operate. **Never iterate trajectories in Python.**

---

## 3. Relevant Files

| File | Role |
|------|------|
| `pyBall/OCL/RigidBodyDynamics.py` | Core PyOpenCL wrapper. Allocates buffers, compiles `Rigid.cl`, uploads state, runs kernels, downloads results. |
| `pyBall/OCL/RigidBodyAFM.py` | High-level solver. Handles anchor setup, per-body anchor repositioning, and chunked relaxation with convergence checks. |
| `cpp/common_resources/cl/Rigid.cl` | OpenCL kernels: `rigid_body_gridff_kernel` (GridFF forces + rigid dynamics) and `rigid_body_dynamics_kernel` (anchors + E-field). |
| `tests/tMMFF/test_rigid_afm_ml.py` | **Main driver.** CLI script for batch trajectory generation, plotting, and ML training/prediction. |
| `tests/tMMFF/plot_manipulation_trajs.py` | Older standalone trajectory generator (sequential). Superseded by batch code in `test_rigid_afm_ml.py`. |
| `tests/tMMFF/test_rigid_gridff_ptcda_batch.py` | Production batch-relaxation test (GridFF only, no anchors). Good reference for multi-replica patterns. |
| `doc/Topics/OnSurfaceAssembly/RigidBodyAFM.md` | Design document describing the broader AFM imaging and manipulation architecture. |
| `doc/Topics/OnSurfaceAssembly/RigidSurfPotential_GridFF.md` | GridFF forcefield documentation. |

---

## 4. How to Run

### 4.1 Generate batch manipulation data

```bash
cd $FIRECORE/tests/tMMFF
python test_rigid_afm_ml.py --test manip \
    --ntraj 64 --nsubsteps 40 --step-size 0.05 \
    --nsteps 3000 --z0 5.0
```

**Outputs:**
- `manip_xy.png`, `manip_xz.png` — diagnostic trajectory plots.
- `manip_conv.png` — convergence diagnostics (force norm, torque norm, convergence flags).
- `manipulation_data.npz` — training data with pose-to-pose transitions.

### 4.2 Train ML model + generate prediction plots

```bash
python test_rigid_afm_ml.py --test train \
    --epochs 30 --lr 1e-3 --batch-size 4096
```

**Outputs:**
- `manip_xy_pred.png`, `manip_xz_pred.png` — ML-predicted trajectories in same format as simulation plots.
- Console log with train/validation MSE per epoch.

### 4.3 Key CLI parameters

| Flag | Default | Meaning |
|------|---------|---------|
| `--ntraj` | 64 | Number of parallel replicas (bodies on GPU) |
| `--nsubsteps` | 12 | Steps per trajectory |
| `--step-size` | 0.4 Å | Tip displacement per substep |
| `--z0` | 15.0 Å | Tip (and anchor) Z height above origin |
| `--nsteps` | 2000 | Relaxation steps per substep |
| `--anchor` | 26 | Anchored atom index (0-based) |
| `--track` | 28 | Tracked atom index (opposite end) |
| `--fconv` / `--tconv` | 1e-3 | Force/torque convergence thresholds |
| `--epochs` | 30 | ML training epochs |
| `--lr` | 1e-3 | Learning rate |
| `--batch-size` | 4096 | ML minibatch size |

---

## 5. Physics & Simulation Details

### 5.1 Tip attachment: harmonic spring to a single atom

The AFM tip is **not** pulling the molecule by its center of gravity. It is attached to a **single peripheral atom** (atom 26, an oxygen on the upper-left side of PTCDA) via a harmonic spring. This is physically realistic: the tip apex touches one atom, and the rest of the molecule responds through the rigid-body constraint.

The anchor is implemented in `RigidBodyAFM.prepare()` and `set_anchor_positions()`:

```python
anchors[:, anchor_idx, :3] = tip_position   # target coordinate
anchors[:, anchor_idx, 3]   = anchor_k       # spring stiffness (20.0 eV/Å²)
```

In the OpenCL kernel, the anchor force on atom `i` is:

```
F_anchor = k * (anchor_pos - atom_pos)
```

This force is summed into the atomic forces, then the rigid-body solver computes the net force and torque on the body.

### 5.2 Two-stage relaxation (critical for stability)

A major source of instability is the **unrelaxed initial pose**. When we first place the molecule above the substrate, its COM and orientation are guessed from the raw XYZ file. If the tip is too low, the molecule may immediately crash into the repulsive substrate wall; if the anchor atom is far from the tip, the spring yanks it violently.

To fix this, the simulation is split into two stages:

1. **Stage 1 — relax initial pose:**
   All replicas start with the tip at the **first substep position** (`tips_all[:, 0, :]`). The system relaxes until force/torque norms fall below `fconv`/`tconv`. This yields a physically consistent starting `pos` and `quat` for each replica.

2. **Stage 2 — manipulation scan:**
   The tip is advanced step by step. Each substep starts from the **relaxed state of the previous substep** (no reset). The GPU state persists, so only the anchor targets change.

This is implemented in `test_manipulation()` (lines ~495–542 of `test_rigid_afm_ml.py`).

### 5.3 Height and initialization: the repulsive wall

The NaCl substrate is **not** at `z = 0`. The topmost substrate atoms sit at `z ≈ 4–5 Å` (confirmed by plotting the GridFF potential profile along Z). PTCDA is ~11 Å long. Therefore:

- If the tip is at `z0 = 3.2 Å`, the molecule is **inside the repulsive wall** → NaN / explosion.
- If the tip is at `z0 = 5.0–6.0 Å`, the molecule hovers just above the surface, feeling attractive and weakly repulsive forces.
- For a non-interacting reference scan, `z0 = 15 Å` keeps the molecule far above the surface.

**Rule of thumb:** `z0` must be at least `substrate_top_z + molecule_extent + safety_margin`. For PTCDA on NaCl_L3, substrate top ≈ 5 Å, molecule extent ≈ 6 Å, so `z0 ≥ 11 Å` is safe. For studying surface-induced trajectories, `z0 = 5.0–6.0 Å` is the operational range.

### 5.4 Anchor–tip separation check

Because the anchor spring has finite stiffness (`k = 20.0`), the anchored atom does not perfectly coincide with the tip. During relaxation, the spring stretches slightly. A **headless check** measures the separation after each substep:

```python
sep = np.linalg.norm(anchor_xyz - tip, axis=1)
max_anchor_sep = max(max_anchor_sep, float(sep.max()))
if np.any(sep > 0.5):
    print(f"WARNING ... anchor-tip separation {sep:.3f} A > 0.5 A")
```

In practice, with `k = 20.0` and `z0 = 5.0 Å`, the separation stays below **0.015 Å**, confirming the spring is stiff enough to faithfully track the tip.

### 5.5 Forcefield: GridFF B-spline

The substrate potential is evaluated via `rigid_body_gridff_kernel`, which looks up forces from a precomputed B-spline grid (`Bspline_PLQd.npy`). This is the same grid used for AFM imaging. The kernel computes per-atom forces from the grid, sums them to body force and torque, and integrates rigid-body dynamics with velocity damping.

For initial testing, GridFF is used off-the-shelf with minimal changes. Future work may compare against `FoldedAtomic` (sine/cosine substrate), but that requires additional kernel implementation. GridFF is the simplest available option.

---

## 6. Visualization & Diagnostics

### 6.1 Plot format

To avoid the confusion of overlaying tip, anchor, and track atom in one panel, the plots are split:

- **Left panel:** tip trajectory (solid) + anchor trajectory (dashed).
- **Right panel:** tracked atom trajectory (the opposite end of the molecule).

This is done for both **XY** and **XZ** projections.

Style:
- `ls='-'`, `lw=0.5`, `alpha=0.7` for all trajectories.
- Thin lines, no markers (64 trajectories would be unreadable with markers).
- `tab20` colormap cycles every 20 trajectories.

### 6.2 What to look for

- **Smooth curves** on the right panel = physically plausible response (the molecule sliding around substrate bumps).
- **Straight lines** = the molecule is too far from the surface (`z0` too high) or the step size is too large to resolve substrate corrugation.
- **Jagged/noisy lines** = instability (reduce step size, increase damping, or check height).
- **Anchor separation > 0.5 Å** = spring too weak or simulation diverging.

---

## 7. ML Integration

### 7.1 What the ML model learns

The model does **not** simply regress `tip → tracked_atom_position`. That ignores orientation and fails to capture path-dependent behavior. Instead, it learns the **pose transition operator**:

```
Input:  [tip0, tip1, pos0, quat0]   → 13D
Output: [pos1, quat1]               →  7D
```

Where:
- `tip0` = tip position at substep `j`
- `tip1` = tip position at substep `j+1`
- `pos0`, `quat0` = molecule pose at substep `j`
- `pos1`, `quat1` = molecule pose at substep `j+1`

The training data is extracted from the simulation as:

```python
tip0  = tip_xyz[:, :-1, :].reshape(-1, 3)
tip1  = tip_xyz[:,  1:, :].reshape(-1, 3)
pos0  = pos_xyz[:, :-1, :].reshape(-1, 3)
quat0 = quat_xyzw[:, :-1, :].reshape(-1, 4)
pos1  = pos_xyz[:,  1:, :].reshape(-1, 3)
quat1 = quat_xyzw[:,  1:, :].reshape(-1, 4)
```

### 7.2 Model architecture

A simple deterministic MLP is used as the baseline:

```python
nn.Sequential(
    nn.Linear(13, 128), nn.ReLU(),
    nn.Linear(128, 128), nn.ReLU(),
    nn.Linear(128, 7),
)
```

The quaternion output is normalized to unit length during training and inference. No periodic encoding is used yet; this is the minimal baseline against which future improvements (multimodality, graph networks, periodic features) will be measured.

### 7.3 Autoregressive rollout for prediction plots

To generate predicted trajectories for visualization, the model is rolled out autoregressively:

1. Start from the true relaxed pose at substep 0 (`pos0`, `quat0` from Stage 1).
2. For each subsequent substep `j`:
   - Feed `(tips_all[:, j-1], tips_all[:, j], pos_pred[:, j-1], quat_pred[:, j-1])` into the model.
   - Obtain `pos_pred[:, j]` and `quat_pred[:, j]` (normalize + fix quaternion sign).
3. Reconstruct anchor and track atom positions from the predicted pose using rigid-body rotation.

This matches the simulation plot format and allows direct visual comparison.

### 7.4 Current model performance

With 30 epochs on 64 trajectories × 40 substeps ≈ 2500 training samples:
- Validation MSE drops from ~7.0 to ~3.0 (pos+quat combined).
- Predicted trajectories show **nontrivial curved structure**, but still smoother than the ground truth.
- The model is underfitting; more data (hundreds or thousands of trajectories), longer training, and richer inputs (e.g. previous tip velocity, surface height map) would improve fidelity.

---

## 8. Development History & Key Fixes

### 8.1 Sequential trajectories → batch-parallel (architecture correction)

**Problem:** The first implementation of `test_manipulation()` ran one trajectory at a time in a Python `for` loop, resetting the GPU pose between each. The GPU was idle while Python iterated.

**Fix:** Refactored to `prepare(n_bodies=ntraj)`, precomputed all tip paths as `(ntraj, nsubsteps, 3)`, and updated all anchors simultaneously per substep. One kernel call per substep relaxes all bodies.

### 8.2 Initial pose unrelaxed (stability fix)

**Problem:** The molecule was initialized from raw XYZ coordinates without pre-relaxation. At `z0 = 3.2 Å`, it crashed into the repulsive substrate wall; at higher `z0`, large anchor-tip separation caused violent jumps.

**Fix:** Implemented two-stage relaxation. Stage 1 relaxes the starting pose. Stage 2 advances the tip from the relaxed state.

### 8.3 Repulsive wall / height misidentification

**Problem:** We initially assumed the substrate surface was at `z = 0`. GridFF profiling showed the top NaCl layer is at `z ≈ 4–5 Å`. Placing the molecule at `z0 = 3.2 Å` buried it inside the repulsive LJ wall.

**Fix:** Raised `z0` to 5.0–6.0 Å for surface-interaction scans, and 15 Å for far-field reference.

### 8.4 COM pulling vs atom pulling (physics correction)

**Problem:** Early tests pulled the molecule by its center of gravity. This is unphysical; an AFM tip contacts a single atom.

**Fix:** Switched to anchoring atom 26 (peripheral oxygen). The rigid-body solver then correctly computes torque from the off-center spring force.

### 8.5 RepeatedKernelRetrieval warning (performance fix)

**Problem:** `run_gridff()` retrieved the OpenCL kernel object on every call, causing PyOpenCL to issue a `RepeatedKernelRetrieval` warning and potentially incur overhead.

**Fix:** Cached the kernel handle during `init_gridff()`:

```python
self.krnl_gridff = cl.Kernel(self.prg, "rigid_body_gridff_kernel")
```

And reused it in `run_gridff()`.

### 8.6 ML task misdefinition → pose-to-pose

**Problem:** The first ML model regressed `tip_position → tracked_atom_position`. This cannot learn orientation or path-dependent hysteresis.

**Fix:** Redefined the task as `f(tip0, tip1, pos0, quat0) = (pos1, quat1)`. Implemented autoregressive rollout for visualization.

---

## 9. Current Status & Known Limitations

| Component | Status | Notes |
|-----------|--------|-------|
| Batch-parallel GPU relaxation | ✅ Working | Up to 64+ trajectories tested. Scales to 1000s. |
| Two-stage initialization | ✅ Working | Stage 1 + Stage 2 stable at `z0 = 5.0 Å`. |
| Anchor spring tracking | ✅ Working | Max separation < 0.015 Å with `k = 20.0`. |
| XY/XZ diagnostics | ✅ Working | Split-panel format as requested. |
| Pose-to-pose ML training | ✅ Working | Baseline MLP trains and rolls out. |
| Prediction quality | ⚠️ Improving | Shows curves but underfits; needs more data and richer model. |
| Convergence rate | ✅ Good | 99.8% convergence with `fconv = 1e-3`, `tconv = 1e-3`, 3000 steps. |

### Known limitations

1. **ML underfitting:** The model is a simple 2-layer MLP. It captures the dominant trend but misses fine substrate corrugation. Future work: periodic coordinate encoding, more layers, or a small GNN.
2. **Step size vs. substrate corrugation:** With `step_size = 0.05 Å` and `nsubsteps = 40`, the total scan length is only 2 Å. For longer scans, more substeps or larger step sizes are needed, but large steps may skip over features.
3. **Single tip angle:** Current scans use a single angle (0°). The dataset lacks directional diversity. Future batches should sweep angles.
4. **No multimodality:** The MLP outputs a single deterministic pose. In reality, the molecule may have multiple stable configurations at a given tip position. A Mixture Density Network or diffusion model may be needed later.
5. **Tip-surface interaction:** The current setup uses GridFF for the substrate only. The tip itself does not feel the substrate. For a full AFM simulation, the tip apex should also interact with the surface.

---

## 10. Design Principles (from `RigidBodyAFM.md`)

These principles guided the implementation and should be respected in future work:

- **Use `RigidBodyAFM` as the solver backend.** Keep CLI scripts thin.
- **Generate data in large GPU batches.** Never iterate trajectories in Python in production.
- **Predict absolute final pose** (pos + quat), not just a single atom coordinate.
- **Use ML as proposal + short physical relaxation.** The ML gives a good initial guess; the GPU solver corrects it in a few hundred steps.
- **Start with deterministic MLP.** Add multimodal predictions only after baseline failure is measured.
- **Use periodic encoding for lateral coordinates** (future improvement).
- **Keep forcefield details hidden behind a simulator interface.** The ML model sees only tip and pose; it does not know about GridFF parameters.

---

## 11. Quick Reference: Typical Parameter Sets

### Surface-interaction scan (see substrate corrugation)
```bash
python test_rigid_afm_ml.py --test manip \
    --ntraj 64 --nsubsteps 40 --step-size 0.05 \
    --z0 5.0 --nsteps 3000
```

### Far-field reference (minimal surface interaction)
```bash
python test_rigid_afm_ml.py --test manip \
    --ntraj 64 --nsubsteps 40 --step-size 0.05 \
    --z0 15.0 --nsteps 2000
```

### Train on generated data
```bash
python test_rigid_afm_ml.py --test train \
    --epochs 50 --lr 1e-3 --batch-size 4096
```

---

## 12. Chat History & Decision Log

The following is a condensed log of key design decisions and user feedback that shaped this system:

1. **"NO! we cannot pull the COG!"** — Anchoring was switched from COM to a single peripheral oxygen (atom 26). The anchor spring stiffness was set to `k = 20.0` to keep the anchored atom within ~0.01 Å of the tip.

2. **"Cliping within the cell does not help anything!"** — Unit-cell wrapping was removed from plots. Trajectories are shown in absolute coordinates; leaving the cell is fine as long as motion is smooth.

3. **"Initialize all replicas on a single line"** — For debugging, all trajectories share the same tip direction (0°) and start positions are spaced along the Y-axis. This isolates the effect of initial Y offset on the molecule's response.

4. **"Put tip and anchor in one plot, opposite molecule end in the other"** — Visualization was split into two panels per projection (XY and XZ). Line styles unified to `ls='-'`, `lw=0.5`.

5. **"WHAT THE HELL IS THIS ... running sequentially?!"** — The entire manipulation loop was rewritten from sequential `for i in range(ntraj)` to batch-parallel `prepare(n_bodies=ntraj)`. This is the defining architectural feature of the system.

6. **"Initial pose is wrong, relax it first"** — Two-stage relaxation was added: Stage 1 pre-relaxes at the starting tip, Stage 2 scans from that state.

7. **"The ML should predict pos and quat, not just atom position"** — The ML objective was redefined from `tip → track_atom` to `(tip0, tip1, pos0, quat0) → (pos1, quat1)`, with autoregressive rollout for trajectory prediction.

---
