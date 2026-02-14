# USER

look in dicsusion we have here @RRsp3_consolidation_and_optimization_discussion.md 

and look in our implementation in @XPDB_new.cl and the pyOpenCL harness in @XPDB_new.py 

it started to be extremely long, therefore I think it is time to extract just the importaint parts and consolidate them in the modules 

RRsp3.cl and RRsp3.py

we will move there only the cluster-sorted Posotion-based dynamics parts (not the fore/impulse based)

name relevant krenels from @XPDB_new.cl which we should migrate, and eventually what new kernels we should create to implement the two versions we discuss in @RRsp3_consolidation_and_optimization_discussion.md , read the whole dicusion there carefully you understand the tradeoffs and the motivation. 

---

 # GPT-5.1-Codex

Here’s the extraction and a proposal for RRsp3.{cl,py} focusing only on the cluster-sorted PBD path (no force/impulse path).

**Kernels to migrate from XPDB_new.cl (cluster PBD path)**
- `update_bboxes_rigid` – per-cluster AABB reduction to drive ghost search @pyBall/XPDB_AVBD/XPDB_new.cl#122-159.
- `build_local_topology_rigid` – critical: finds ghost atoms and remaps global neighbor indices to local (cluster+ghost) indices @pyBall/XPDB_AVBD/XPDB_new.cl#162-253.
- `compute_collision_cluster_rigid` – cluster-local collisions in LDS, bonded-neighbor exclusion via `neighs_local` @pyBall/XPDB_AVBD/XPDB_new.cl#256-335.
- `compute_ports_cluster_rigid` – node-only ARAP/ports in clusters using local indices, writes node deltas + recoil slots @pyBall/XPDB_AVBD/XPDB_new.cl#411-518.
- `apply_corrections_rigid_ports` – gathers collision + node + recoil deltas (bkSlots) and integrates pos/quat, node-vs-cap split @pyBall/XPDB_AVBD/XPDB_new.cl#521-569.
- (Optional carry-over) `compute_corrections_rigid_ports` – the non-cluster, global PBD version if you still want a non-cluster fallback @pyBall/XPDB_AVBD/XPDB_new.cl#338-408.

These are the only kernels in XPDB_new.cl tied to the clustered position-based dynamics path; the rest are math helpers.

**Harness pieces to mirror in RRsp3.py**
- Cluster buffer setup and launch order in `rigid_cluster_relax_ports_step` (AABB → topology → collisions → ports → apply) @pyBall/XPDB_AVBD/XPDB_new.py#814-911.
- Cluster buffer allocators: `_init_rigid_cluster_buffers`, `_ensure_rigid_cluster_node_buffers`, and the upload helpers for clustered port/stiffness data @pyBall/XPDB_AVBD/XPDB_new.py#320-431.
- bkSlots builder for clustered layout: `make_bk_slots_clustered` @pyBall/XPDB_AVBD/XPDB_new.py#56-88.

**New kernels to add (per the two versions discussed in the markdown)**
1) **Split path (recommended for recoil + global sync)**
- Reuse/migrate the existing split set above, but add:
  - `build_local_topology_rigid_excl` (or extend the current topology kernel) to also map 1st/2nd-neighbor exclusions into local indices for collisions, so the inner collision loop can skip via registers without global indirection.
  - (Optional) A light `clear_cluster_deltas` kernel to zero `dpos_coll`/node delta buffers per inner iteration (currently done via host fills).

2) **Fused experimental path (single cluster kernel)**
- New kernel `solve_cluster_fused` that:
  - Loads cluster+ghosts once, uses the same local index map for both collisions and ports.
  - Accumulates collision and port corrections in one pass (or two blocks inside one kernel) and writes deltas/updates once.
  - Still needs a pre-pass to build local indices (can reuse `build_local_topology_rigid` or the extended version with exclusions).
- If recoil must be preserved without atomics, keep `dpos_neigh` slots inside the fused kernel and apply locally; if cross-cluster recoil is required, the fused version will still need a global apply kernel (lighter than the split path).

**Summary of what moves into RRsp3.cl**
- Core cluster PBD kernels: bbox, topology (with exclusions), collisions, ports, apply.
- Optional fused kernel for experimentation.
- Keep math helpers needed by these kernels.

**Summary of what moves into RRsp3.py**
- Minimal harness to allocate cluster buffers, upload clustered ports/stiffness, bkSlots, and run the split (default) or fused (experimental) sequence.
- Retain the clustered bkSlots constructor and clustered upload helpers; drop unrelated force/impulse paths.

---

# USER

1) we will focus now purely on cluster-sorted variatns
2) look into Gemini discussion in @RRsp3_consolidation_and_optimization_discussion.md and try to map/compare it with existing kernels, which kernels we can basically coppy without changes or with minimal changes and which we need to make from scrath. Respectivally, I would prefer if you copy the kernels from RRsp3.cl (the one which is closest to the intend) and then modify them to intended form discussed. Make detailed plan for each kernel about coresponnede of the original kernel nd the one dicsussed and what are the changes we should make.
is there anything we need to wrote from scrath?
Notes: the kernels in @XPDB_new.cl are already tested to large degree (although not fully), we should keep the recoils dpos, it is very importaint to think deeply about layour and sorting, we encountered many problems when we implemented it last time with this sorting. 
 
 ---

 # GPT-5.1-Codex

Plan: map Gemini discussion → existing cluster PBD kernels → required changes and scratch work.

Per-kernel correspondence (cluster-sorted, keep recoils):

1) `update_bboxes_rigid` @XPDB_new.cl#122-159  
- Copy as-is. Already cluster-sorted, LDS reduction OK.  
- Optional minor change: if we add per-cluster node counts, can carry through but not required.

2) `build_local_topology_rigid` @XPDB_new.cl#162-253  
- Base to copy, but needs extension:  
  - Also map exclusion lists (1st/2nd neighbors) to LOCAL indices for collisions; store in a new buffer `excl_local` (int4 or int8) so collision kernel avoids global lookups.  
  - Possibly cap ghosts/slots; make MAX_GHOSTS configurable via define.  
  - Keep remap of neighs_global → neighs_local for ports and bkSlots correctness.

3) `compute_collision_cluster_rigid` @XPDB_new.cl#256-335  
- Copy with minimal change: swap bonded-exclusion check to use `excl_local` from topology instead of re-reading neighs_local 1st neighbors each time (faster, supports 2nd neighbors).  
- Ensure it reads cluster+ghost from LDS; preserve existing double-counting halving and mass-weighted push.  
- Keep `dpos_coll` output untouched.

4) `compute_ports_cluster_rigid` @XPDB_new.cl#411-518  
- Copy with minor changes:  
  - Inputs should accept `excl_local`? Not needed; ports use neighs_local already.  
  - If we adopt “node-first per cluster, node buffer indexed by inode” keep same nnode_per_group handling.  
  - Optionally allow accumulate flag to start from collision deltas for nodes (current code supports accumulate_dpos).  
  - Keep recoil `dpos_neigh` writing; keep `dpos_node`/`drot_node`.

5) `apply_corrections_rigid_ports` @XPDB_new.cl#521-569  
- Copy with minimal change:  
  - If we add collision exclusions it doesn’t change apply.  
  - Keep cap compensation factor (dx*=2 for caps), bkSlots gather, relaxation factor.  
  - Ensure it matches node-first layout (inode = grp*nnode_per_group + lid).

6) `compute_corrections_rigid_ports` (non-cluster) @XPDB_new.cl#338-408  
- Optional: can omit from RRsp3 if we only keep cluster path. If kept, unchanged.

New/from-scratch items:

A) Topology exclusion mapper (new): either extend `build_local_topology_rigid` or add `build_local_exclusions_rigid` to convert global exclusion lists (1st/2nd neighbors) into local indices aligned with the ghost list (store -1 if absent). This is needed to get branch-free exclusion in collisions.

B) (Optional) Fused kernel (experimental) `solve_cluster_fused` if you still want to try single-kernel collisions+ports; otherwise skip. Keeping recoils and cross-cluster sync means we likely stay split per discussion.

C) (Optional) Small zeroing kernel for cluster deltas if we want device-side clears instead of host fills (not strictly needed).

Layout/sorting implications to preserve:
- Global arrays packed cluster-contiguously; within cluster nodes first then caps; padding to GROUP_SIZE.
- nnode_per_group passed to ports/apply; bkSlots built with clustered node indexing (use make_bk_slots_clustered from XPDB_new.py).
- Ghost list length per cluster; local indices: 0..63 cluster, 64..64+g ghosts.
- Recoils `dpos_neigh` use inode*4+k indexing; bkSlots must be consistent with node-only indexing.

Summary of copy vs modify vs new:
- Copy unchanged: `update_bboxes_rigid`; likely `compute_corrections_rigid_ports` (if kept).
- Copy with minor edits: `build_local_topology_rigid` (add local exclusions); `compute_collision_cluster_rigid` (consume local exclusions); `compute_ports_cluster_rigid` (maybe minor args cleanup); `apply_corrections_rigid_ports` (layout sanity).
- New: exclusion-mapping helper (or merged into topology). Optional fused solver and zeroing helper.


---

# USER

OK, now we should test it properly. I saw in my prevous test were bing problems with proper collisions between atoms, atoms were not coliding as they should because they were not mapped properly into clusters or the mapping of local indexes (ghost atoms, neibhbor exclusion) was not do ne probperly

for that prupose we should add several debug potions into the kernels
1) detailed debug prints gated by #define macros and compile-flags
2) add option to selectively turn on/of the collision and bond evaluation selectively. 
3) then we should run synthetic tests e.g. for H2O checking if atoms which should colide really colide and which should not colide does not colide. We should make clearly defined list of interaction which should be seen in deteild debug print writen e.g. as python dictionary or python list of list. And then we should have post-processing function which analyse the debug prints from the kernels (only one iteration) which analyse it. We should think deeply about how to make such regorous and general testing system, that is quete essential and we can reuse that in other projects

the other think we should do is visial debugging where user (me) can pick and pull atoms, this will allow me to test many features quicky and in flaxyble way and gain some intiution what is the problem. I see we already have such functionality in @test_rigid_XPBD_molecules.py and in @XPTB_utils.py but the matplotlib bakcend is slow and problematic. And I wouldprefer to migrate to pyQt and Vispy for fster full-featured 3D interface. Look how we did it in @MolecularPlacerVisPy.py @SequencePlacerVisPy.py 

The mouse picking in 3D may be tricky (especially when using zoom and camera rotation) , so maybe we should test that separately firt. But we should also have alternative pseudo-2D mode where we still use Vispy D rendered but with camera being in top view so we operate only in x-y plane. Note that by any means we should use othographi projection, not perspective, as it is easier to navigate and also the mouse picking works better in that view.

now read the referenced files and make detailed notes what to implement and how to run the test then, Make detailed notes and plan.

---

there are some ideas about the verification / debugging strategy they to analyse them and thin how to integrate them into our apporach

@RRsp3_XPBD_verification_strategy.md 

---

Another synthetic test which we mist add is preservation of linear and angular momentum over the whole system during relaxation, this is common pitfall and the reason why we added the explicit recoil dpos_neigh which are a bit stupid (slow down, inclreased complexity, architectural limitations), but they were necessary to ensure momentum conservtion. Later we maye tery to get rid of them, but only after we debug these kenrels and ensure momentum conservation with them. 

make anlos detaile dnotes about this

---

 # GPT-5.1-Codex



Below is a detailed plan and notes for adding robust debugging + synthetic tests and a Vispy/Qt interactive picker.

### A) Kernel-level debugging instrumentation (RRsp3.cl)
- Add compile-time gates:
  - `#ifndef DBG_COLL` / `#ifndef DBG_PORT` macros to wrap collision/bond debug prints.
  - `#define DBG_COLL_PRINT 0`, `#define DBG_PORT_PRINT 0` to compile but keep off by default; allow enabling via build options `-DDBG_COLL_PRINT=1`.
  - Optional per-cluster/atom filters: `DBG_COLL_GRP`, `DBG_COLL_LID`, `DBG_PORT_GRP`, `DBG_PORT_LID` (defaults -1 for “all”). This avoids flooding.
- Add selective disable toggles:
  - `#define ENABLE_COLL 1`, `#define ENABLE_PORT 1`; guards in `compute_collision_cluster_rigid` and `compute_ports_cluster_rigid` so we can run one without the other.
- Debug content to print (one iteration):
  - Collisions: `grp`, `lid`, `jloc`, `ri/rj`, `d2`, `rsum`, whether excluded, resulting `dl`/`n*dl*w_i`.
  - Ports: `grp`, `lid`, `k`, `jloc`, `K`, `r_arm`, `tip`, `diff`, [dist](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:249:4-250:53), `impulse_mag`, `P`, recoil slot idx.
  - Ghost info (once per group): `ghost_count`, first few `ghost_indices`.
  - Exclusion mapping (once per atom when printing collisions): `excl1_local`, `excl2_local`, and `neighs_local` so we can validate mapping.
- Guard prints with `if (DBG_COLL_PRINT && (DBG_COLL_GRP<0 || grp==DBG_COLL_GRP) && (DBG_COLL_LID<0 || lid==DBG_COLL_LID))` etc.
- Keep recoils intact.

### B) Synthetic test harness design (Python)
- Reuse current [test_RRsp3_smoke.py](cci:7://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/test_RRsp3_smoke.py:0:0-0:0) structure; extend to a “debug mode”:
  - Build a small scene (H2O or two H2Os) with controlled geometry:
    - Case 1: overlap H–H between molecules → expect collision.
    - Case 2: bonded pairs → must be excluded.
    - Case 3: second neighbors excluded by `excl2` → ensure not colliding.
  - Pack with [pack_molecules_contiguous](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:31:0-196:5), node-first, pad to 64.
  - Build neighs/bkSlots as before; build exclusions via [make_exclusions_1st_2nd](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/RRsp3.py:57:0-98:23) (or custom list to enforce a known exclusion set).
  - Enable kernel debug with build options:
    - `build(options=["-DDBG_COLL_PRINT=1","-DDBG_COLL_GRP=0","-DDBG_COLL_LID=0","-DENABLE_PORT=0"])` for collision-only, or similarly for ports.
  - Run one iteration (coll-only or port-only) to keep logs small.
  - Capture stdout/stderr; parse with a simple Python parser that builds a dict of observed interactions:
    - Expected list/dict (e.g., `[('coll', i, j), ...]`) vs. observed; assert presence/absence.
- Post-processing function:
  - Regex on debug lines → tuples `(grp, lid, jloc, status, dl, ... )`.
  - Compare against expected allowed pairs and forbidden pairs; raise descriptive assertion.
- Make tests deterministic: set positions explicitly (e.g., place two H atoms at 0.8 Å overlap; set radii=1.0 to force penetration).

### C) Vispy/Qt interactive picker plan
- Reuse patterns from `MolecularPlacerVisPy.py` and [SequencePlacerVisPy.py](cci:7://file:///home/prokop/git/FireCore/pyBall/SequencePlacerVisPy.py:0:0-0:0):
  - Scene setup: `SceneCanvas`, orthographic camera for pseudo-2D top view (no perspective). Use Turntable but set fov=0 or use `scene.cameras.PanZoomCamera` in ortho mode.
  - Picking: Vispy provides `visuals.Markers` with `visual.set_gl_state`; for simple picking, do CPU-side ray-plane since ortho top view:
    - With ortho top view (camera looking down -Z), picking reduces to nearest atom in XY screen coords; map mouse pos to data coords via `canvas.transforms`.
  - Controls:
    - Left-click select nearest atom; drag to move in XY; constrain Z=0 or plane.
    - Keys to toggle collision/port enable flags, trigger one RRsp3 step, and display resulting deltas.
  - Rendering:
    - Markers for atoms colored by element; lines for bonds (reuse colors/sizes from `MolecularPlacerVisPy.py` / [SequencePlacerVisPy.py](cci:7://file:///home/prokop/git/FireCore/pyBall/SequencePlacerVisPy.py:0:0-0:0)).
  - Data pipeline:
    - Use the same packed buffers already built for RRsp3.
    - On interaction: update host pos, re-upload pos/radius/quats, run [step_cluster](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/RRsp3.py:265:4-343:16), download pos/quats, update visuals.
  - Optional: “debug overlay” showing ghost count per group and drawing ghost atoms in a different color when debug mode is on.

### D) Steps to implement next
1) Add debug macros/toggles to [RRsp3.cl](cci:7://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/RRsp3.cl:0:0-0:0) (collision/port enable, print gates, filters).
2) Extend [RRsp3.py](cci:7://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/RRsp3.py:0:0-0:0) to accept build options for debug macros and pass them to `cl.Program(...).build(options=...)`.
3) Extend [test_RRsp3_smoke.py](cci:7://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/test_RRsp3_smoke.py:0:0-0:0) into `test_RRsp3_debug.py`:
   - Deterministic geometries for collide/no-collide cases.
   - Expected interaction list and parser to validate kernel prints.
   - Command-line flags to choose collision-only vs port-only, and to enable debug build flags.
4) Prototype Vispy picker:
   - Start with top-view ortho (PanZoomCamera), XY-only moves.
   - Use the existing molecular loading/placement patterns from `MolecularPlacerVisPy.py`/[SequencePlacerVisPy.py](cci:7://file:///home/prokop/git/FireCore/pyBall/SequencePlacerVisPy.py:0:0-0:0).
   - Wire a button to run one RRsp3 step and update visuals.

If you want, I can proceed to implement items 1–3 (kernel debug toggles + debug test parser) before moving to the Vispy picker.

---

Here’s how to fold the verification-strategy ideas into our RRsp3 debugging and testing approach:

### Kernel instrumentation (RRsp3.cl)
- Keep/extend the macros already added (`ENABLE_DEBUG_PRINTS`, `LOG_TOPOLOGY`, `LOG_MAPPING`, `LOG_COLLISION_CHECK`), plus:
  - Filters: `DEBUG_GID_START/END`, optional `DBG_COLL_GRP/LID`, `DBG_PORT_GRP/LID`.
  - Toggles: `ENABLE_COLL`, `ENABLE_PORT` to run collision-only or port-only.
  - Standardized collision log: `COLL: MeG=%d OtherL=%d OtherG=%d Dist=%.4f Action=%s`.
  - Standardized topology log: ghosts found, and global→local mappings for neighs/exclusions.
- Ensure ghost global IDs are available in collision logs (for ghosts, reconstruct `OtherG = ghost_indices[grp*MAX_GHOSTS + (j-GROUP_SIZE)]`).

### “Interaction truth” definition (Python)
- A small module, e.g. `verification_truth.py`, containing tuples `(a, b, expected_action)` with actions: `EXCLUDE_BOND`, `EXCLUDE_ANGLE`, `COLLIDE`, `IGNORE_FAR`.
- Start with two H2O molecules in two clusters (close enough to collide), mirroring the example in the strategy doc.

### Debug build/run path (Python harness)
- Build options: `["-DENABLE_DEBUG_PRINTS", "-DDEBUG_GID_START=0", "-DDEBUG_GID_END=16", "-DENABLE_COLL=1", "-DENABLE_PORT=0"]` etc.
- Run only topology + collision (or topology + ports) for a single iteration to keep logs small.
- Capture stdout/stderr (e.g., subprocess run) to a log string/list.

### Log parser + validator
- Regex for `COLL:` lines → `(me_g, other_g, action)`.
- Map actions vs truth:
  - Expected `EXCLUDE_*` → must see `Action=SKIP_EXCL`.
  - Expected `COLLIDE` → must see `Action=COLLIDE`.
  - `IGNORE_FAR` → no collision logged or `Action=TOO_FAR`.
- Fail fast with descriptive mismatches; emit PASS summary.

### Deterministic synthetic geometries
- Use [pack_molecules_contiguous](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:31:0-196:5) as before, but place molecules deliberately:
  - Case A: overlap inter-molecule pairs to force collisions.
  - Case B: spaced to avoid collisions; verify `TOO_FAR`.
- Radii and positions set explicitly; exclusions built via [make_exclusions_1st_2nd](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/RRsp3.py:57:0-98:23) or explicit lists to match the truth table.

### Vispy/Qt interactive debugger (next step)
- Reuse `MolecularPlacerVisPy.py` / [SequencePlacerVisPy.py](cci:7://file:///home/prokop/git/FireCore/pyBall/SequencePlacerVisPy.py:0:0-0:0) patterns.
- Ortho top-view (PanZoomCamera or Turntable with fov=0); XY-only dragging to simplify picking.
- Picking: nearest atom in view coordinates; update host pos → reupload → run one RRsp3 step → refresh markers/lines.
- Add UI toggles for `ENABLE_COLL/PORT` and a “debug step” button; optionally display ghost counts and exclusions in an overlay.

### Concrete next actions (if you want them implemented)
1) Add the remaining debug filters/toggles and standardized logs to RRsp3.cl.
2) Add build-options plumbing in RRsp3.py to pass debug macros.
3) Create `verification_truth.py` and a new `test_RRsp3_debug.py`:
   - constructs deterministic two-H2O scene,
   - runs collision-only with debug prints on,
   - parses logs and validates against the truth table.
4) Prototype Vispy top-view picker for manual poke-and-pull checks.

---

Notes: momentum conservation synthetic test (validate recoil dpos_neigh)

Goal
- Confirm global linear and angular momentum are preserved during relaxation when collisions + ports run with recoils.
- Use small, controlled systems (2 molecules packed into 2 clusters) so we can compute expected totals and verify after one or a few iterations.

Metrics to check
- Linear momentum P = Σ m v (we infer v ≈ (x_new - x_old)/dt if velocities not explicit; or store prev_pos to compute).
- Angular momentum L = Σ r × m v + I ω (for nodes; caps have only translational term).
- Check before and after a single iteration; difference should be ~ numerical roundoff.

Test setup ideas
1) Two H2O molecules (same as collision debug) with zero initial velocity, small random perturbations to positions and orientations to induce forces but no net initial momentum.
2) Another variant: assign opposite small velocities to two atoms so net P=0, net L=0, and see if stays zero.
3) Ensure masses/inertias:
   - m from elements; invM for nodes/caps; Iiso = 0.4 m R^2 for simple sphere (as used previously).
4) Radii chosen so there is light collision interaction (e.g., shift molecules closer) but not exploding.

What to log/compare
- Before step: positions, quats, previous positions (or explicit velocity array), masses.
- After one iteration (collision+ports, with recoils enabled): recompute P_total, L_total.
- Assert |ΔP| < eps_lin, |ΔL| < eps_ang (set small tolerances, e.g., 1e-5–1e-6 relative to typical magnitudes).
- Also record per-cluster P/L to see if any cluster is leaking momentum.

Kernel toggles
- Keep recoils dpos_neigh ON (this is the mechanism under test).
- Optionally run with `ENABLE_PORT=0` then `ENABLE_COLL=0` to isolate which part breaks conservation if it fails.
- Enable debug prints (selective) for momentum test only if mismatch is detected; otherwise keep off.

Harness steps (Python)
1) Build packed scene (pack_molecules_contiguous) with known geometry and padding; upload masses/radius/quat.
2) Save initial pos/quats (and prev_pos if needed).
3) Run [step_cluster](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/RRsp3.py:265:4-343:16) once.
4) Download pos/quats; compute P/L before and after.
   - v ≈ (x_new - x_old)/dt
   - L_trans = Σ r × m v
   - L_rot (nodes): use ω if available; otherwise infer from quaternion delta: ω ≈ 2 * imag(dq)/dt for small rotations.
5) Compare totals; print diagnostics if above tolerance (include per-cluster P/L and max per-atom contribution).
6) Optional: rerun with collisions-only and ports-only to isolate source.

Edge cases to include
- Padding atoms (mass=0) must not contribute.
- Ghost mapping should not break conservation: ensure ghosts don’t double-count; only real atoms included in P/L.
- Caps recoil: confirm bkSlots indexing matches dpos_neigh slots; if wrong, L will drift.

Outputs expected
- PASS summary: ΔP, ΔL magnitudes, per-cluster P/L.
- On FAIL: detailed print of atoms with largest momentum change, cluster IDs, and whether collision or port-only run caused it.

If you want, I can implement this as a new test script (e.g., `test_RRsp3_momentum.py`) that reuses the packing utilities and computes P/L with the above checks.