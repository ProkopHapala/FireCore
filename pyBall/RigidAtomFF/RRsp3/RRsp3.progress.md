
# RRsp3 Progress Report (Session Summary)

## Background, goals, and physical motivation (from the consolidation discussions)

- **Algorithmic context**: Cluster-sorted, Jacobi-style Position-Based Dynamics (PBD) with “ports” (ARAP-style constraints) and explicit “recoil” buffers to maintain momentum conservation. Each workgroup handles one cluster (e.g., 64 atoms). Nodes appear first, caps later; padding is allowed but treated as inactive (mass=0, NaN positions).
- **Pipeline (cluster PBD)**:
  1. **update_bboxes_rigid**: per-cluster AABB for broad-phase.
  2. **build_local_topology_rigid**: finds ghosts in overlapping clusters; remaps global neighbor/exclusion indices to local (0..63 cluster, 64.. ghosts).
  3. **compute_collision_cluster_rigid**: uses local memory (cluster + ghosts) and local exclusions; writes collision deltas.
  4. **compute_ports_cluster_rigid**: node-only ARAP/ports; writes node deltas, angular deltas, and recoil `dpos_neigh` for neighbors.
  5. **apply_corrections_rigid_ports**: gathers collision + node + recoils (bkSlots) and applies to pos/quat with relaxation.
- **Motivation for ghost/local mapping**: Avoid divergent global memory access; collisions/ports operate purely on local indices in LDS; ghost mapping includes exclusions to keep collision loop O(local).
- **Exclusion strategy**: Bonded and 2nd neighbors must be skipped in collisions; local exclusion mapping preferred to avoid global lookups. Register-based exclusions are best in fused variants; in split path, topology maps global exclusions to local indices.
- **Sorting strategy**: Strong recommendation: cluster-sorted global layout; within each cluster: nodes first, caps next, then padding. This minimizes divergence and improves coalescing. Padding is allowed but must be inert (mass=0, invM=0).
- **Recoil rationale**: Recoil `dpos_neigh` is kept to preserve linear/angular momentum in a Jacobi scheme; apply kernel gathers these via bkSlots.
- **Momentum conservation test**: Check Σ m·dx = 0 and Σ (r×m·dx) + Σ I·dθ = 0 after a step; padding (m=0) excluded. Recoil is critical for conservation.

## What we implemented/fixed in this session

### Kernel/Harness
- Added **fixmask** buffer with bitwise control:
  - bit1/2/4: fix X/Y/Z; bit8: clamp Z to 0.
  - Padding atoms are pinned by default; pos set to NaN when invM=0 to expose accidental use.
- **Collision kernel**: barrier-safe (no early return before barrier); skips padding/fixed atoms; ghost building skips invalid atoms; collisions skip invalid masses/radii.
- **Apply kernel**: respects fixmask (axes pinned, clamp z); skips padding atoms.
- **Host** ([RRsp3.py](cci:7://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/RRsp3.py:0:0-0:0)):
  - [upload_fixmask](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/RRsp3.py:218:4-222:68), passed to apply; fixmask initialized to zeros.
  - Exports NaN padding when invM=0.
  - XYZ export filtered to real atoms only in momentum test; padding masked out of invariants.

### Tests
- `test_RRsp3_debug.py`: PASS after fixes; logs expected COLLIDE/EXCL.
- [test_RRsp3_momentum.py](cci:7://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/test_RRsp3_momentum.py:0:0-0:0): PASS after fixes; multi-step run with 6-atom XYZ only; invariants ~1e-8.

### GUI (Vispy)
- **Manual camera**: LMB pick-only; RMB rotate; wheel zoom (manual); debug prints confirm zoom.
- **Run/Stop loop**: QTimer calls solver steps; while dragging, the dragged atom is re-applied before/after each step so it stays under the cursor during relaxation.
- **Padding not pickable**: padding fixed and masked; picker ignores NaN/fixed atoms.
- Debug prints show wheel and drag-hold events firing.

## Remaining issue
- User reports wheel “not responding” previously; added robust delta decoding and terminal logs; current run shows `[WHEEL] ...` and `[CAM] zoom ...` with changing distances, so zoom events are now confirmed. If still unresponsive visually, we can add a UI zoom spinbox and reset-view button, but logs show the handler fires.

## Equations and invariants
- Momentum checks (Jacobi correction invariants):
  - Linear: ΔP = Σ m·dx = 0
  - Angular: ΔL = Σ r×(m·dx) + Σ I·dθ = 0
  - I_iso ≈ 0.4 m R² (used in kernel and test)
- Collision resolution (per pair i,j):
  - n = d/|d|, w_tot = w_i + w_j (+ eps)
  - dl = (r_sum - dist) / w_tot; collision correction for i: dx_i += n * dl * w_i (with the 0.5 factor in kernel for symmetry)
- Ports (ARAP/constraint impulse):
  - Impulse magnitude ~ dist / (w_i + w_j + w_ang + α), with α = 1/(K dt²); recoil stored in `dpos_neigh` and gathered by bkSlots.

## Files touched (for reference)
- Kernels: [RRsp3.cl](cci:7://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/RRsp3.cl:0:0-0:0)
- Harness: [RRsp3.py](cci:7://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/RRsp3.py:0:0-0:0)
- Tests: `test_RRsp3_debug.py`, [test_RRsp3_momentum.py](cci:7://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/test_RRsp3_momentum.py:0:0-0:0)
- GUI: [VispyUtils.py](cci:7://file:///home/prokop/git/FireCore/pyBall/VispyUtils.py:0:0-0:0), [test_RRsp3_vispy.py](cci:7://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/test_RRsp3_vispy.py:0:0-0:0)
- Utilities: [XPTB_utils.py](cci:7://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:0:0-0:0) (trajectory/pick helpers reused)

### Purpose of key files added/updated this session
- `RRsp3_progress.md` — this progress log summarizing goals, fixes, tests, and remaining issues.
- `VispyUtils.py` — reusable Vispy/Qt viewer with manual camera (LMB pick-only, RMB rotate, wheel zoom), 2D/3D pick modes, fixmask awareness, drag-state signals.
- `test_RRsp3_debug.py` — synthetic interaction test (two H2O clusters) parsing kernel debug logs for COLLIDE vs SKIP_EXCL; ensures ghost/exclusion mapping works.
- `test_RRsp3_momentum.py` — momentum conservation test with multi-step option and XYZ export (real atoms only); checks Σ m·dx and Σ r×m·dx+Σ I·dθ invariants.
- `test_RRsp3_vispy.py` — interactive GUI: pick/drag, run/stop loop, fixmask pinning, clamp toggle, manual camera; wires to `VispyUtils` and RRsp3 kernels.
- `RRsp3_XPBD_verification_strategy.md` — strategy notes for synthetic interaction/momentum tests and kernel debug gates.


## 2026-02-15 (today) (stretch + vispy + scan)

### Vispy improvements (today)
- Kept GUI intact while adding convergence changes; re-verified `VISPY_LOG_LEVEL=debug VISPY_DEBUG_GL=1 python3 test_RRsp3_vispy.py` runs clean.
- Interaction remains: LMB pick/drag, RMB rotate, wheel zoom with debug logs; dragged atom re-applied before/after solver step to stay under cursor; padding/fixed atoms ignored by picker.

### Headless convergence & scans (today)
- Added robust scan mode to `test_RRsp3_convergence.py`: sweeps bmix/dt/relax, early-stop on `max_err` thresholds (1e-2 / 1e-3), detects blow-ups, saves per-run `errors.csv` + `convergence.png` + full `traj.xyz` (configurable `--save_every`), and aggregates `summary.csv` under `scan_outputs/<timestamp>/`.
- Distortion workflow for polyacetylene: freeze ports from reference, then stretch along X by `--fstretch` (default 1.5), collisions off by default (`k_coll=0`).
- Fixed padded-atom NaN poisoning: `RRsp3.upload_state(..., nan_padding=False)` now used by headless scripts so padded slots stay finite; optional `--nan_padding` keeps old debug behavior. This enabled errors to decay to ~1e-3 for strong stretches when scanning bmix/dt/relax.
- Example scan (stretch 1.5, k_coll=0) produced converged runs to ~1e-3: e.g., bmix=0.9, dt=1.0, relax=1.0 reached 1e-3 at ~222 iterations (see `scan_outputs/scan_20260215_140055/summary.csv`).

### Visualization details (expanded)

- Fixed absolute radius scaling (collision glyphs) in `VispyUtils.AtomScene._px_per_world_ortho`: pixels-per-world now derived from the inverse scene transform of a 1-pixel screen delta (no try/except, no fallbacks, orientation-independent). set_zoom triggers `_redraw()` so glyphs rescale immediately. Result: R_coll=1.0 now visually matches 1–1.5 Å bond lengths (disks overlap as expected).
- Planned/added overlays (via `test_RRsp3_vispy.py` checkboxes):
  - Neighbor bonds (from global neighs, dedup j>i).
  - Node → port-tip rays (rotate `port_local` by node quaternion).
  - Port-tip → target atom bonds (bkSlots mapping).
  - Debug vectors: `dpos` (total) and `dpos_neigh` (slot recoil) per atom.
- Halo/inbox connectors recolor to cluster palette when “color by groups” is enabled; otherwise magenta. Padding/fixed atoms are skipped by picker/overlays; drag-hold re-applies grabbed atom around solver steps.
- Hardcoded H2O bonds removed; bonds now come from packed molecule data. GUI verified with `VISPY_LOG_LEVEL=debug VISPY_DEBUG_GL=1 python3 test_RRsp3_vispy.py` in the current worktree.

### Momentum/physics design docs (new)
- `RRsp3_momentum_design.md` — comprehensive design for XPBD heavy-ball momentum on pos/quat, compliance/weights, testing strategy (XYZ + convergence plots).
- `RRsp3_momentum_design_old.md` — prior draft kept for reference.
- `RRsp3_discussion_2.md` — historical discussion on port PBD vs XPBD and cluster handling (retained as context).

## Physics & Convergence Analysis (2026-02-16)

### Problem Identification
- Initial scan showed slow convergence (~73-140 iterations to reach 1e-2 error) even with reasonably stiff constraints.
- **Root Cause 1 (Fixed)**: Missing neighbor angular mass in `w_total`. 
  - The PBD constraint impulse denominator `w_total` originally included only `w_i` (linear), `w_j` (linear), and `w_ang_i` (angular for self).
  - For Node-Node bonds (backbone), the neighbor `j` also rotates. Neglecting `w_ang_j` made `w_total` too small, causing `impulse_mag` to be too large (overshooting/instability).
  - *Fix Applied*: Added `if (j_isnode) w_total += w_ang;` (assuming symmetric inertia) to `RRsp3.cl`. This improves stability allowing stiffer `alpha`.

- **Root Cause 2 (Identified - CRITICAL)**: **Torque Under-driving via Double-Counting Logic**.
  - To handle the fact that both thread `i` and thread `j` process the same bond, the linear impulse `P` is scaled by 0.5 (`impulse_mag *= 0.5`).
  - Linear updates are correctly reconstructed: `i` adds `0.5*P` (direct) + `0.5*P` (recoil from `j`).
  - **Angular updates are NOT shared**: Thread `i` applies torque to `i`. Thread `j` applies torque to `j`. Thread `j` does *not* send angular recoil to `i`.
  - **Result**: Thread `i` applies torque derived from `0.5*P`. Total torque on `i` is 50% of the required value. This under-damps rotation significantly, leading to slow geometric convergence.
  - *Solution*: Decouple scaling. Use `0.5 * P` for linear updates (`dpos_node`, `dpos_neigh`), but `1.0 * P` for angular updates (`drot_node`).

### Convergence Scan Results (Preliminary)
- Run: `scan_outputs/scan_20260216_113142/`
- Settings: `stiffness=5000` (stiff), `dt=2.0`, `invM=1` (uniform).
- With Cause 1 fixed (but Cause 2 present):
  - `bmix=0.9, dt=2.0`: Reaches 1e-2 in ~143 iterations.
  - `bmix=0.9, dt=5.0`: Reaches 1e-2 in ~136 iterations.
  - High `bmix` (0.99) is unstable with high stiffness.
- Created D3.js visualization: `scan_outputs/scan_20260216_113142/convergence_viz.html` (open in browser to explore Bmix vs Iterations).

### Scan 2 Analysis (Torque Decoupling Applied)
- Run: `scan_outputs/scan_20260216_114304/`
- **Result**: Convergence speed for the **Stretch** test remained largely unchanged (~140-150 iterations for best stable cases) despite the torque fix.
- **Physics Interpretation**: 
  - The test case is a pure linear stretch aligned with the bond axis. In this configuration, $\mathbf{r} \times \mathbf{n} \approx 0$, so torque is negligible. Thus, improving torque scaling (angular response) does not accelerate relaxation of a pure linear mode.
  - The bottleneck is confirmed to be **Jacobi Diffusion** on a linear chain. Information propagates diffisively ($O(N^2)$ steps) because the solver is parallel (Jacobi).
  - `bmix` (momentum) helps (reducing iters from >500 to ~150), but cannot beat the fundamental diffusion limit for long chains without a hierarchical or serial (Gauss-Seidel) component.

### Recommendations for "Few-Iteration" Convergence
To achieve game-physics-like speeds (3-5 iterations) for this system:
1. **Inner Local Iterations (Block Gauss-Seidel)**: 
   - Since the cluster (64 atoms) fits in Local Memory (LDS), we can perform multiple *inner* iterations of the solver loop *inside the kernel* using LDS.
   - This effectively solves the "local" system (the whole C16 molecule) rigidly before writing back.
   - **Plan**: Modify `compute_ports_cluster_rigid` to iterate updates on `l_pos` locally before outputting `dpos`. This should collapse the iteration count to ~1-5 per frame.

### Visualization
- Interactive plot generated: `scan_outputs/scan_20260216_114304/convergence_viz.html`
  - Open in browser to correlate `bmix`, `dt`, and `relaxation` with iteration count.


## 2026-02-16 (today) (rotation solver variants + scans + XYZ fix)

### Kernel refactors
- Added shared helpers in `RRsp3.cl`: `quat_delta_rotvec`, `quat_conj`, `add_hessian_rr`, `mat3_outer_add`; rewired `substep_optimized`, `shapematch`, and `eigen` kernels to remove duplicate math. Restored correct Frobenius scaling in shapematch polar decomposition.

### Rotation kernels overview
- `compute_ports_cluster_rigid_orig`: one-pass XPBD; halved linear impulse for node–node double-count; torque implicitly halved; neighbor angular inertia ignored. Lightest but under-drives rotation.
- `compute_ports_cluster_rigid_substep_optimized`: “Born-Oppenheimer-ish” — register-only Newton substeps tighten rotation ignoring linear masses/recoil, then a normal XPBD linear pass. Fast/cheap, but not angular-momentum conserving.
- `compute_ports_cluster_rigid_shapematch`: builds covariance, polar decomposition (Newton–Schulz), best-fit rotation, then XPBD linear pass. Stable/accurate rotation, heaviest register/ALU, massless rotation.
- `compute_ports_cluster_rigid_eigen`: Davenport q-method (implicit K) via quaternion power iteration; lighter than shapematch, best-fit rotation, massless rotation, still only self inertia in linear denom.

### Convergence scans run (headless)
- Grid: dt {0.1, 1.0, 10.0}, bmix {0.3, 0.5, 0.8, 0.9}, relax=1.0, K=5000, masses=1, k_coll=0, kernels {eigen, current/orig}, rms {1.0, 0.0 (fixed rotation)}, save_every=5.
- ICs:
  - Stretch (1.5×): outputs in `scan_outputs/scan_20260216_180710/` (per-run `errors.csv`, `traj.xyz`, `convergence.png`, aggregate `summary.csv`, `convergence_viz.html`).
  - Rand (0.2 Å, seed=0): outputs in `scan_outputs/scan_20260216_180829/` (same artifacts as above).

### XYZ atom-name fix
- `test_RRsp3_convergence.py` now strips element suffixes when writing XYZ (e.g., `C_2 -> C`) for viewer compatibility (Jmol, etc.).

### Run notes
- From repo root: `PYTHONPATH=/home/prokop/git/FireCore:$PYTHONPATH python3 pyBall/XPDB_AVBD/test_RRsp3_convergence.py --scan ...`
- From `pyBall/XPDB_AVBD`: `PYTHONPATH=..:$PYTHONPATH python3 test_RRsp3_convergence.py --scan --mol ../cpp/common_resources/mol/polyacetylene_C16.mol2 ...` (no extra `pyBall/` prefix).

Here’s a focused, detailed recap of the three “new” rotation-solvers vs. `compute_ports_cluster_rigid_orig`, including the Born-Oppenheimer–like massless rotation behavior and the Newton substep method.

### Baseline: `compute_ports_cluster_rigid_orig`
- **Rotation/translation** solved together in one XPBD pass.
- Denominator: `w_total = w_i + w_j + w_ang + alpha`; `w_ang` uses only the current node’s inertia.
- Double-count guard: linear impulse halved if neighbor is a node; torque implicitly halved too (so rotation under-driven).
- Pros: simple, light registers. Cons: rotation slow/under-driven; ignores neighbor rotational inertia; linear and angular coupled in the same (often large) denominator.

### `compute_ports_cluster_rigid_substep_optimized` (Newton substeps; “Born-Oppenheimer-ish”)
- **Two-phase approach**:
  1) **Rotation-only Newton substeps** in registers, using neighbors from LDS, ignoring linear masses and recoil (denom ~ `w_ang + alpha_rot`). This is the “massless rotation / Born-Oppenheimer-like” pre-tightening: it instantaneously finds a better orientation without sharing that torque with neighbors.
  2) Then a **full XPBD linear pass** using the updated `qi` with normal denom (`w_i + w_j + w_ang + alpha`) and linear double-count halving; torque derived from that linear impulse (neighbor angular inertia still not included).
- Delta rotation output is the total `qi` vs `qi_orig`.
- Pros: cheap, register-only inner loop; aggressively tightens orientation before linear solve; low global traffic in the substep.  
- Cons: angular momentum not conserved (rotation adjusted without recoil); quality depends on substep count/eps; still only own inertia in denom; not rigorously XPBD for rotation.

### `compute_ports_cluster_rigid_shapematch` (Polar decomposition; “infinite rotational stiffness”)
- Builds weighted covariance `A = Σ w * (xj−xi) * r_localᵀ` (stiffness as weight), runs Newton–Schulz polar decomposition (4 iters) → optimal rotation matrix, converts to quaternion.
- Then a standard XPBD linear pass using that optimal rotation; double-count halving on linear.
- Pros: near-exact best-fit rotation for current geometry; very stable per step.  
- Cons: heaviest register/ALU footprint (mat3 ops, 4 iterations); rotation is effectively massless/instant (angular momentum not conserved); still uses only own inertia in the linear denom.

### `compute_ports_cluster_rigid_eigen` (Davenport q-method / implicit K-matrix; quaternion power iteration)
- Centers targets by weighted centroid; accumulates covariance (6 unique terms) and applies the implicit K-matrix via **power iteration in quaternion space** (4 iters).
- Directly returns a unit quaternion (no mat→quat conversion); then standard XPBD linear pass with double-count halving.
- Pros: lighter than shapematch (no 3×3 polar loop), direct quaternion eigensolve, stable, enforces unit-norm.  
- Cons: still “massless” rotation (no recoil), still only own inertia in linear denom; heavier than substep.

### Massless rotation vs. coupled rotation
- **Massless/BO-like**: Shapematch and Eigen rotate to the best-fit orientation without sharing torque; substep_optimized does a pre-rotation phase ignoring linear masses/recoil. All three are not angular-momentum conserving.
- **Coupled rotation** (only in orig/current linear pass): rotation and translation share the same impulse; torque is derived from the linear impulse. Under-driven in orig because impulse is halved for nodes.

### Rotational weighting differences
- Orig/current: own inertia only; torque implicitly halved if neighbor is a node.
- Substep_optimized: rotational substeps ignore masses entirely; final linear uses own inertia only.
- Shapematch/Eigen: rotation from geometry only; final linear uses own inertia only.

### Practical pros/cons summary
- **Fastest/cheapest**: substep_optimized (register-only inner loop). Good when you want quick rotational tightening with minimal cost, but not momentum-consistent.
- **Most accurate geometric rotation**: shapematch (exact polar) and eigen (Davenport). Eigen is lighter than shapematch; both are massless rotations.
- **Most conservative/minimal registers**: orig/current (but rotation slow/under-driven, ignores neighbor inertia).

### Testing order (as discussed)
1) Simple small system to ensure they run (no crashes): collision off, known ports, check constraint error drops.
2) Compare fixed-rotation vs full rotation on pure stretch (translation-only should converge; good check of linear path).
3) Random small quaternion perturbations with dpos fixed vs free to isolate rotational tightening.
4) Only then performance/convergence on C16: measure iterations to 1e-2/1e-3, and consider separate residual metrics for translational vs rotational DOFs (dpos vs drot).

If you want next steps: I can script a small harness to A/B these kernels on a tiny molecule (e.g., ethylene or a 4-node toy) to confirm correctness and report per-kernel iteration counts and residuals split into translational vs rotational parts.

---

## 2026-06-17: Momentum Conservation Regression Fixed + Rigorous Testing Framework

### Problem Found
The `current` kernel (default) failed angular momentum conservation on **asymmetric molecules** (CH3OH = L2), while `orig` passed. On symmetric H2O (L0), both passed — the bug was hidden by symmetry.

### Root Cause
In `compute_ports_cluster_rigid` (current kernel), the **node-node double-counting guard was asymmetric**:
- **Linear impulse** was correctly halved (`impulse_lin *= 0.5`) since both node threads evaluate the same bond.
- **Angular impulse / torque** was **NOT halved** (kept full magnitude).

This meant torque on node-node constraints was applied **twice**, breaking angular momentum conservation for any molecule where torque matters (i.e., asymmetric inertia tensors).

### Fix Applied
In `RRsp3.cl` `compute_ports_cluster_rigid`:
- Removed the decoupled `impulse_lin` vs `impulse_ang` split.
- Both linear and angular contributions now use the **same halved impulse** for node-node bonds:

```c
if (j_isnode) { impulse_mag *= 0.5f; }
float3 P = n * impulse_mag;
sum_dpos   += P * w_i;
sum_dtheta += cross(r_arm, P) * invI;
dpos_neigh[idx] = (float4)(-P * w_j, 0.0f);
```

### Testing Infrastructure Added

#### 1. Shared utility module: `RRsp3_utils.py`
- Molecule builders: `make_ch3oh_geometry()`, `make_ch2nh_geometry()`, `make_tri3_geometry()`
- Momentum analysis: `compute_momentum_deltas()` with `massless_rot` flag
- Topology truth generation: `build_interaction_truth()`
- Debug log parsing: `parse_kernel_debug_logs()`, `validate_interactions()`
- 2D planar constraint helpers: `perturb_state_planar()`, `make_planar_fixmask()`

#### 2. Consolidated headless test: `test_RRsp3_headless.py`
Subcommands:
- `momentum` — linear + angular momentum conservation per step
- `topology` — collision/exclusion mapping via kernel debug prints
- `smoke` — basic NaN + range check
- `suite` — full matrix across systems and kernels

#### 3. Kernel debug print enhancement
Verbosity levels + component bitmask + workgroup targeting in `RRsp3.cl`:
```c
-DENABLE_DEBUG_PRINTS -DDEBUG_VERBOSITY=3 -DDEBUG_COMPONENTS=7 -DDEBUG_TARGET_WG=0
```

#### 4. Two-tier momentum check
- **Massless alignment kernels** (`substep_optimized`, `shapematch`, `eigen`): rotation is geometric optimization, not physical dynamics → check only `dL_trans` (linear displacement contributes to angular momentum), exclude `dL_rot`.
- **Full rotational dynamics** (`current`, `orig`): include both `dL_trans + dL_rot`.

#### 5. Nodes-only test system: `tri3` (L5)
3-node triangle ring, no capping atoms, no `bkSlots` recoil complexity. Used to isolate core kernel momentum behavior from cap-atom assembly issues.

#### 6. Kernel caching
Fixed `RepeatedKernelRetrieval` warnings by caching `cl.Kernel` instances in `RRsp3._kernels` dict.

### Current Test Results (5 steps, epsP=1e-4, epsL=1e-3)

| Kernel | L0 H2O (sym) | L2 CH3OH (asym) | L5 tri3 (nodes-only) | Interpretation |
|--------|-------------|------------------|----------------------|----------------|
| `current` | **PASS** | **PASS** | **PASS** | Fixed; now conserved |
| `orig` | **PASS** | **PASS** | **PASS** | Baseline correct |
| `substep_optimized` | **PASS** | **FAIL** dL=0.005 | **PASS** | Linear part non-conserving on asymmetry |
| `shapematch` | **FAIL** dL=0.067 | **FAIL** dL=0.102 | **FAIL** dL=0.240 | Even linear part violates conservation |
| `eigen` | **FAIL** dL=0.032 | **FAIL** dL=0.080 | **FAIL** dL=0.126 | Even linear part violates conservation |

### What This Means
- **`current` and `orig` are now the only momentum-conserving kernels** across symmetric, asymmetric, and nodes-only systems.
- **`substep_optimized`**: fails on asymmetric molecules because its linear displacement does not conserve angular momentum (massless rotation + XPBD linear pass mismatch).
- **`shapematch` and `eigen`**: fail on **all** systems, including symmetric H2O and nodes-only tri3. This indicates a deeper issue: their linear recoil/writeback does not satisfy momentum conservation even for pure linear constraints. The massless rotation is not the only problem.

### Outstanding Issues
1. **`shapematch` / `eigen` linear momentum violation**: Even on symmetric/node-only systems, `|dL| ~ 0.03–0.24` with `epsL=1e-3`. This is a **real bug**, not just "massless rotation by design". Their recoil writeback or linear impulse scaling is wrong.
2. **`substep_optimized` asymmetric failure**: On L2 CH3OH, `|dL|=0.005`. This may be acceptable for some use cases but should be understood.
3. **Topology test subcommand**: Needs debug print capture fix (kernel printf via subprocess not reliable).
4. **Interactive GUI consolidation**: `test_RRsp3_gui.py` (Vispy) not yet built.

### Files Modified
- `RRsp3.cl` — torque double-counting fix, enhanced debug macros
- `RRsp3.py` — kernel caching, `download_drot_node()`
- `RRsp3_utils.py` — new shared module (molecule builders, momentum, topology truth)
- `test_RRsp3_headless.py` — new consolidated test script

---

## 2026-06-18: Convergence Speed Test Added

### New Subcommand: `convergence`
Added `convergence` subcommand to `test_RRsp3_headless.py` that loads an XYZ file (e.g. pentacene), computes port targets from the **original unstretched** geometry, then stretches the initial positions in-plane and measures how many iterations each kernel needs to relax back.

### Metrics
- **mean_step**: mean(|dpos_node|) + mean(|dpos_neigh|) — solver step size, universal across kernels
- **geo_err**: mean ||pos[j] - tip_i|| — geometric constraint violation (tip-to-atom distance)

### Key Findings

**Massfull kernels (`current`, `orig`)**: Converge well even with 2x stretch. After 500 iterations with stretch=2.0, dt=0.1, beta=0.5:
- `current`: step=4.0e-04, geo_err=1.48e-02
- `orig`: step=2.7e-04, geo_err=9.19e-03

**Massless kernels (`substep_optimized`, `shapematch`, `eigen`)**: Fail to converge for large stretches. After 500 iterations with same params, geo_err stays at ~0.3–0.5 Å. Root cause: central-force projection (impulses only along atom centerlines) cannot provide the **lateral** displacement component needed for large geometric corrections.

**Momentum acceleration (beta) effect** on `current`, stretch=2.0, dt=0.1:
- beta=0.0: step=7.4e-03 (slowest)
- beta=0.5: step=3.5e-03
- beta=0.9: step=1.3e-05 (near-converged in 200 iters)

**dt effect** on `current`, beta=0.5, stretch=2.0:
- dt=0.05: step=3.0e-03, geo_err=3.3e-01
- dt=0.1: step=3.5e-03, geo_err=1.3e-01
- dt=0.5: step=1.7e-03, geo_err=2.3e-02 (best)

Larger dt + higher beta gives faster convergence for massfull kernels because the physical inertia term stabilizes the momentum acceleration.

### Plotting
The `--plot` flag generates a dual-panel figure (solver step size + geometric error) with all parameter sweeps overlaid. Labels include beta and dt values; title shows sweep ranges.

Example usage:
```bash
python test_RRsp3_headless.py convergence --kernel current --beta_list "0.0,0.5,0.9" --dt_list "0.05,0.1,0.5" --plot out.png
```

### Files Modified
- `test_RRsp3_headless.py` — `cmd_convergence`, `--beta_list`, `--dt_list`, `--plot`
- `RRsp3_utils.py` — `build_packed_system_from_xyz`, `quat_rotate_vec`, `compute_mean_constraint_error`
