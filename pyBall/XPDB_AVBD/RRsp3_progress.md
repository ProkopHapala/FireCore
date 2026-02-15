
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
- `test_RRsp3_debug_runner.py` — helper to run the debug test in a clean subprocess with build options and NPY temp payloads.
- `test_RRsp3_momentum.py` — momentum conservation test with multi-step option and XYZ export (real atoms only); checks Σ m·dx and Σ r×m·dx+Σ I·dθ invariants.
- `test_RRsp3_vispy.py` — interactive GUI: pick/drag, run/stop loop, fixmask pinning, clamp toggle, manual camera; wires to `VispyUtils` and RRsp3 kernels.
- `RRsp3_XPBD_verification_strategy.md` — strategy notes for synthetic interaction/momentum tests and kernel debug gates.


## **2026-02-15 (today)** (stretch + vispy + scan)

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
