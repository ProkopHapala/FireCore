# USER

look we have this two forcefields
/home/prokop/git/FireCore/pyBall/XPDB_AVBD
/home/prokop/git/FireCore/pyBall/XPBD_2D
there we have test scripts like 
@test_XPBD_2D.py 
@test_rigid_XPBD_molecules.py 
@test_TiledJacobi_molecules.py 

and visualization / IO /  testing usitilites in
@XPBD_2D_utils.py 
@XPTB_utils.py 

Because I plan to test many more forcefields and avoid code duplication, I think we should consolidate all the unitilites, visualizers, interface to @AtomicSystem.py @atomicUtils.py , and also some refactored functions (generalized) from inside those scripts  into one place. We wnat to have both 2D and 3D visualizer with functional mouse picking and with callbacks for visualization of buffers downloaded from OpenCL (or other posiion/force buffers) writen in some general reusable way in sngle module /home/prokop/git/FireCore/pyBall/FFutil.py place in

Note - ssumption is that this FFutils will be then comined with arbitrary pyOpenCL forcefield wrpaer such as @XPBD_2D.py and @XPDB_new.py 
which means it should NOT have any direct thependence on thos, it should be general with repsect to replacement of theese, they shoudl comunicate over general numpyarrays and callbacks...

 Now your job is to make detailed notes where is which usefule function noew, and where are some sematically common or similar blocks in the script which make seanse to refactor into some common function. And them weite detailed plan what to do.

---

# GPT-5.1-Codex

Below is a map of the current utilities/visualizers/hooks and overlapping functionality, followed by a consolidation plan for a new shared `pyBall/FFutil.py`. (Sources cited.)

### Where things live now

#### XPBD_2D utilities (2D focus, matplotlib)
- Visualization:
  - [LiveViz2D](cci:2://file:///home/prokop/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:30:0-176:29): live scatter/line plot with optional labels, ports, bond lines, dynamic marker sizing, persistent view. @pyBall/XPBD_2D/XPBD_2D_utils.py#31-178
  - Handles drawing bonds, port markers, atom→port and port→neighbor segments, info text. Uses stored `_last_pos` for picking.
- Picking/interaction:
  - [attach_picker_2d(viz, sim, pick_radius, verbose)](cci:1://file:///home/prokop/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:386:0-432:15): mpl callbacks for mouse press/move/release; stores pick state; expects `viz._last_pos` and `sim` with setters. @pyBall/XPBD_2D/XPBD_2D_utils.py#387-433
- Geometry/topology prep:
  - [_make_permutation_from_groups](cci:1://file:///home/prokop/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:179:0-194:41), [_invert_permutation](cci:1://file:///home/prokop/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:197:0-202:14). @pyBall/XPBD_2D/XPBD_2D_utils.py#180-200
  - [reorder_system_by_groups(...)](cci:1://file:///home/prokop/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:205:0-249:64): reorder atoms by groups, optional node-first, returns perm/inv and remapped bonds/types/charges. @pyBall/XPBD_2D/XPBD_2D_utils.py#206-250
  - [setup_from_mol(sim, mol, ...)](cci:1://file:///home/prokop/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:252:0-350:5) and [setup_from_xyz](cci:1://file:///home/prokop/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:353:0-355:134): build 2D topology (neighs/bk/stiffness/ports), reorder nodes-first, upload topology+state, optional perturbations. @pyBall/XPBD_2D/XPBD_2D_utils.py#253-357
  - Chain helpers: [create_chain_topology](cci:1://file:///home/prokop/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:435:0-438:16), [setup_chain_molecule](cci:1://file:///home/prokop/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:441:0-484:5). @pyBall/XPBD_2D/XPBD_2D_utils.py#436-485
- Metrics:
  - [compute_port_error(pos, rot, neighs, bks, port_local, nnode, port_n)](cci:1://file:///home/prokop/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:363:0-383:23): max/RMS port closure error. @pyBall/XPBD_2D/XPBD_2D_utils.py#364-384
  - Momentum helpers exist elsewhere in file (not yet cited in snippets; present above line 200).
- Rotations:
  - [rot_from_angles(ang)](cci:1://file:///home/prokop/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:358:0-360:75): convert angles to 2D complex rot vectors. @pyBall/XPBD_2D/XPBD_2D_utils.py#359-363
- Verbosity flag `VERBOSE` + setter. @pyBall/XPBD_2D/XPBD_2D_utils.py#23-28

#### XPBD_2D test usage patterns
- `test_XPBD_2D.py` uses [LiveViz2D](cci:2://file:///home/prokop/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:30:0-176:29) and [attach_picker_2d](cci:1://file:///home/prokop/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:386:0-432:15) plus pick-driven mass override and direct `sim.set_atom_pos/vel/omega` for dragging; downloads buffers every 10 steps for diagnostics; uses [compute_port_error](cci:1://file:///home/prokop/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:363:0-383:23). @pyBall/XPBD_2D/test_XPBD_2D.py#80-160

#### XPDB_AVBD utilities (3D rigid focus)
- Permutation/helpers:
  - [invert_permutation](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:5:0-11:14), [apply_permutation_to_bonds](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:14:0-28:27). @pyBall/XPDB_AVBD/XPTB_utils.py#6-30
  - [pack_molecules_contiguous(...)](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:31:0-196:5): packs multiple molecules into workgroup-sized chunks, node-first reorder (or explicit perm), pads to group, returns combined pos/elems/group ids/bonds and per-mol perms. @pyBall/XPDB_AVBD/XPTB_utils.py#32-198
- Simple geometry transforms:
  - [as_unit](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:199:0-204:16), [deform_shift_atom](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:207:0-213:14), [deform_scale_along_direction](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:216:0-232:14). @pyBall/XPDB_AVBD/XPTB_utils.py#200-233
- Small molecule builder:
  - [make_h2o_geometry(add_angle, k_oh, k_hh)](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:235:0-262:21), [bonds_to_max_L0](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:265:0-270:49). @pyBall/XPDB_AVBD/XPTB_utils.py#236-272
- IO:
  - [ensure_outdir](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:273:0-275:15), [print_run_header](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:278:0-281:29), [plot_residual_series](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:284:0-302:15). @pyBall/XPDB_AVBD/XPTB_utils.py#274-304
  - [load_xyz](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:305:0-330:24) (with optional q), [masses_from_elems](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:333:0-346:12) (H/C/N/O only). @pyBall/XPDB_AVBD/XPTB_utils.py#306-347
  - [write_xyz_with_ports](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:385:0-399:23), [write_pdb_trajectory](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:402:0-411:24). @pyBall/XPDB_AVBD/XPTB_utils.py#386-413
- Quaternions and perturbations:
  - [quat_mul](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:349:0-357:14), [normalize_quat](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:360:0-363:29), [perturb_state(pos, quat, pos_scale, rot_scale, rng)](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:366:0-382:28). @pyBall/XPDB_AVBD/XPTB_utils.py#350-383
- Visualization (3D):
  - [LivePortViz](cci:2://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:414:0-467:29): 3D scatter + bonds via ports, keeps labels, optional force quiver; stores `_last_pos`. @pyBall/XPDB_AVBD/XPTB_utils.py#415-470
  - [attach_picker_3d(viz, sim, pick_radius_px, verbose)](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:470:0-537:15): projects 3D to screen, finds closest within pixel radius, tracks mouse3 position at fixed depth; minimal coupling to sim (just returns pick dict). @pyBall/XPDB_AVBD/XPTB_utils.py#471-538
  - [plot_state_with_ports](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:540:0-544:61) singleton wrapper to reuse [LivePortViz](cci:2://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:414:0-467:29). @pyBall/XPDB_AVBD/XPTB_utils.py#541-545

#### XPDB_AVBD test usage patterns
- [test_rigid_XPBD_molecules.py](cci:7://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/test_rigid_XPBD_molecules.py:0:0-0:0):
  - Loads xyz via [load_xyz](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:305:0-330:24) or [load_molecule_any](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/test_rigid_XPBD_molecules.py:57:0-73:35) (AtomicSystem) and reorders nodes-first. @pyBall/XPDB_AVBD/test_rigid_XPBD_molecules.py#22-75
  - Builds rigid topology, ports (local via quat conj), uploads to sim. @pyBall/XPDB_AVBD/test_rigid_XPBD_molecules.py#88-141
  - Visualization toggling: [ensure_viz2d](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/test_rigid_XPBD_molecules.py:205:4-213:35) uses [LiveViz2D](cci:2://file:///home/prokop/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:30:0-176:29)/[attach_picker_2d](cci:1://file:///home/prokop/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:386:0-432:15); [ensure_viz3d](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/test_rigid_XPBD_molecules.py:215:4-222:35) uses [LivePortViz](cci:2://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:414:0-467:29)/[attach_picker_3d](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:470:0-537:15); view switch on key 'v'. @pyBall/XPDB_AVBD/test_rigid_XPBD_molecules.py#185-224
  - Emits frames with port visualization and optional picking-driven interaction (similar to 2D test pattern).
- [test_TiledJacobi_molecules.py](cci:7://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/test_TiledJacobi_molecules.py:0:0-0:0):
  - Bond array builders ([build_bond_arrays_with_angles](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/test_TiledJacobi_molecules.py:74:0-177:69), [_bonds_to_adj](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/test_TiledJacobi_molecules.py:194:0-206:20)) producing fixed-stride arrays; IO dump ([dump_xpdb_inputs_text](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/test_TiledJacobi_molecules.py:34:0-71:37)). @pyBall/XPDB_AVBD/test_TiledJacobi_molecules.py#35-179
  - Uses `pyBall.atomicUtils` and MMFFL; general molecule loading functions further down.

### Overlaps & common patterns ripe for refactor

- Permutations/reordering:
  - [_invert_permutation](cci:1://file:///home/prokop/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:197:0-202:14) (2D) vs [invert_permutation](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:5:0-11:14) (3D) and [apply_permutation_to_bonds](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:14:0-28:27). Merge into a single permutation module.
  - Node-first reordering logic appears in [reorder_system_by_groups](cci:1://file:///home/prokop/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:205:0-249:64), [setup_from_mol](cci:1://file:///home/prokop/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:252:0-350:5), [reorder_nodes_first](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/test_rigid_XPBD_molecules.py:32:0-54:31), [pack_molecules_contiguous](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:31:0-196:5) (node-first with degree>1). These can be unified with a common helper that accepts `node_mask` or `perm` override and returns perm/new bonds and mapping.
- Topology build:
  - 2D uses `build_neighs_bk_from_bonds_2d` and port construction; 3D rigid builds neighs/bLs/bKs and ports via quaternion rotation. Shared scaffolding: convert bonds -> neighbor matrices, build port_local (2D vector / 3D vector) from positions and orientation, generate bk slots/stiffness.
- Visualization:
  - [LiveViz2D](cci:2://file:///home/prokop/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:30:0-176:29) and [LivePortViz](cci:2://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:414:0-467:29) share the same idea: persistent figure, scatter, lines for bonds/ports, labels. Different dimensionality and projection/picking.
  - Picking helpers [attach_picker_2d](cci:1://file:///home/prokop/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:386:0-432:15) and [attach_picker_3d](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:470:0-537:15) share structure: event handlers, distance threshold (world vs screen pixels), storing `pick` dict with `idx`, `active`, `mouse`, optional world coord, and `_last_pos` from viz.
  - Both tests perform pick-driven constraints (set pos/vel/omega or mass) — could be a generic “apply picking to sim via callbacks” hook.
- IO:
  - [load_xyz](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:305:0-330:24) and `AtomicSystem`-based loader coexist; there’s custom [write_xyz_with_ports](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:385:0-399:23), [write_pdb_trajectory](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:402:0-411:24).
  - [print_run_header](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:278:0-281:29), [ensure_outdir](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:273:0-275:15), [plot_residual_series](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:284:0-302:15) small general helpers.
- Quaternion/rotation utilities vs 2D [rot_from_angles](cci:1://file:///home/prokop/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:358:0-360:75).
- Geometry perturbation (positions + rotations) exists in 3D; 2D perturbation is ad-hoc in [setup_from_mol](cci:1://file:///home/prokop/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:252:0-350:5).
- Buffer download/diagnostics:
  - 2D test downloads state, computes momentum and [compute_port_error](cci:1://file:///home/prokop/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:363:0-383:23). Similar diagnostics could be generalized (port closure error, NaN checks, momentum) for both 2D/3D.

### Proposed FFutil.py structure (general, numpy-first, no hard dependency on specific sims)

Organize into cohesive sections; each function pure and data-oriented, minimal deps (numpy, matplotlib when needed):

1. **Imports & Verbosity**
   - `VERBOSE`, [set_verbose](cci:1://file:///home/prokop/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:25:0-27:24).

2. **Permutation / Reordering**
   - [invert_permutation(perm)](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:5:0-11:14).
   - [apply_permutation_to_bonds(bonds, perm_inv)](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:14:0-28:27).
   - [reorder_nodes_first(pos, bonds, *, node_mask=None, nodes_first=True, types=None, charges=None)](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/test_rigid_XPBD_molecules.py:32:0-54:31) returning reordered arrays, perm, inv; works for 2D/3D pos.
   - [reorder_system_by_groups](cci:1://file:///home/prokop/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:205:0-249:64) retained/merged to support explicit groups (from 2D utils).

3. **Topology builders (pure CPU)**
   - `build_neighs_from_bonds(bonds, n_atoms, max_deg=4)` for generic; thin wrappers call 2D/3D-specific builders if needed.
   - `make_bk_slots_generic(neighs, nnode, natoms, max_deg=4)` (wrap existing 2D/3D logic).
   - Port generation:
     - `ports_from_positions(pos, neighs, nnode, dim=2 or 3, rot=None)` returning `port_local`, `port_n`. For 3D, uses optional quats to rotate world offsets into body frame; for 2D, uses complex rot or identity.
   - Stiffness helpers (bond-based to per-node slot) parameterized by k_bond and dim.

4. **Rotation utilities**
   - `rot_from_angles_2d(ang)`.
   - [quat_mul](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:349:0-357:14), [normalize_quat](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:360:0-363:29), [quat_rotate](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/test_rigid_XPBD_molecules.py:176:4-182:43), `perturb_state_3d(pos, quat, ...)` (from XPTB_utils).

5. **Geometry utilities**
   - [deform_shift_atom](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:207:0-213:14), [deform_scale_along_direction](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:216:0-232:14), [make_h2o_geometry](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:235:0-262:21), [bonds_to_max_L0](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:265:0-270:49).
   - Chain/regular topology generators (2D chain) generalized to 3D option.

6. **Visualization**
   - Common base interface? Keep two classes:
     - [LiveViz2D](cci:2://file:///home/prokop/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:30:0-176:29) (largely as-is).
     - `LiveViz3D` (rename from [LivePortViz](cci:2://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:414:0-467:29)), parameterized for optional forces/labels.
   - Shared picking helpers:
     - [attach_picker_2d(viz, *, pick_radius=..., on_pick=None)](cci:1://file:///home/prokop/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:386:0-432:15) returns pick dict; optional callback hook signature `(pick_state)`.
     - [attach_picker_3d(viz, *, pick_radius_px=..., on_pick=None)](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:470:0-537:15) with projection code.
   - Lightweight wrappers:
     - `plot_state_with_ports_2d(...)` and `_3d(...)` reusing singleton viz instances.
   - Optional “pick application” helper:
     - `apply_pick_to_sim(pick, sim, *, set_pos=True, zero_vel=True, zero_omega=True, pin_mass=None, get_state=None)` where caller passes lambdas to set pos/vel/omega/mass because sim API differs.

7. **IO**
   - [load_xyz](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:305:0-330:24) (reuse).
   - [write_xyz_with_ports](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:385:0-399:23) (keep).
   - [write_pdb_trajectory](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:402:0-411:24) (keep).
   - [ensure_outdir](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:273:0-275:15), [print_run_header](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:278:0-281:29).
   - Possibly a thin [masses_from_elems](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:333:0-346:12) with optional element table override.

8. **Diagnostics**
   - [compute_port_error](cci:1://file:///home/prokop/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:363:0-383:23) generalized to 2D/3D (supports rot optional).
   - `check_nans(pos, rot, vel, omega)` returning counts.
   - `compute_momentum_linear/angular` thin wrappers expecting arrays and masses (reuse existing formulas from tests).

9. **Packing**
   - [pack_molecules_contiguous](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:31:0-196:5) kept (for GPU workgroup packing), but optional dependency on group_size; make it independent of specific FF.

10. **AtomicSystem hooks**
    - Generic loader:
      - `load_from_atomic_system(mol: AtomicSystem, *, bAllNodes=False, dim=3, infer_bonds_if_missing=True, reorder_nodes=True)` returning elems, pos (2D/3D slice), bonds, nnode, node_mask, perms.

### Refactor targets (what to move into FFutil)

- From [XPBD_2D_utils.py](cci:7://file:///home/prokop/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:0:0-0:0): [LiveViz2D](cci:2://file:///home/prokop/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:30:0-176:29), [attach_picker_2d](cci:1://file:///home/prokop/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:386:0-432:15), [compute_port_error](cci:1://file:///home/prokop/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:363:0-383:23), permutation helpers, `setup_from_mol/setup_from_xyz` (split into AtomicSystem loader + topology builder), chain topology helpers, [rot_from_angles](cci:1://file:///home/prokop/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:358:0-360:75).
- From [XPTB_utils.py](cci:7://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:0:0-0:0): permutation helpers, [pack_molecules_contiguous](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:31:0-196:5), quat/perturb utilities, IO functions, [LivePortViz](cci:2://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:414:0-467:29) + [attach_picker_3d](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:470:0-537:15), geometry transforms, [make_h2o_geometry](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:235:0-262:21), [bonds_to_max_L0](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPTB_utils.py:265:0-270:49), `print_run_header/ensure_outdir/plot_residual_series`.
- From tests (to avoid duplication in future tests):
  - Pick-driven sim mutation patterns (mass pinning, set pos/vel/omega) — turn into optional callback helper.
  - View toggling patterns (2D/3D) — small helper to register key handler.

### Interfaces to keep generic (no sim coupling)
- Visualizers store `_last_pos`; pickers should only emit pick dict (idx, active, mouse, mouse3, zproj) and leave sim interaction to caller via passed callbacks.
- Topology builders should output numpy arrays (neighs, bkSlots, stiffness, port_local, port_n, perms) without touching OpenCL.
- IO should operate on numpy arrays and plain lists (elems, bonds).
- Packing should accept arbitrary molecule dicts (elems, pos, bonds, nnode, optional perm) and not assume a specific kernel beyond group_size/nnode.

### Suggested next steps
1. Create `pyBall/FFutil.py` and move the shared, pure helpers there, preserving existing behavior; leave thin shim imports in old modules to maintain backward compatibility temporarily.
2. Split sim-specific uploads from data prep: e.g., `build_topology_from_mol(...)` returns dict; [setup_from_mol](cci:1://file:///home/prokop/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:252:0-350:5) in 2D can wrap it to call `sim.upload_*`.
3. Unify pickers with a small common core; expose optional callbacks for “on_press/on_release/on_drag” so tests can inject pin-mass or position override logic.
4. Add minimal tests or adapt existing ones to import from FFutil instead of module-specific utils; start with `test_XPBD_2D.py` and [test_rigid_XPBD_molecules.py](cci:7://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/test_rigid_XPBD_molecules.py:0:0-0:0) by switching imports to the shared functions.
5. Document expected array shapes and dim handling in FFutil docstrings; keep loud validation (ValueError) consistent with current style.

Let me know if you want me to draft the initial `FFutil.py` skeleton and update the two test scripts to consume it.

---

# USER

with this in mind, we should think about consolidation of deubing tools, especially the visual debugging (including loading data, eddicient searching and comparing in output files, reading buffers, making sure debug print have same format, but most importaintly unoform visualization interface avoiding code duplicity, with raycasting into 3D scene for efficient mouse picking, and ploting indexes of atoms and bonds (to easily compare them with debug prints or printed tolpology buffers), and fast BLIT based matplitlib animation interface) we should all this consodidate to single module FFutils.py

there is comprehensive but preminary analysis of this problem or use in XPDB methods @FFutils_visual_debugging_consolidation.md 
now we should dig deeper and analyse how to use it alos for MMFF and UFF from both pyOpenCL,  C++/OpenCL and CPU/C++

what I fogrot to say and what is rather importaint is that we have these utilities to read raw arrays (buffers) from C++ and wrap them into numpy array objects (allowing slicing), but the memory ownership is still on C++ side so we should not allow python Garbage collector to destroy them otherwise it crashes (or there can be some double free).

 see @MMFFmulti_lib.cpp @libUtils.h @MMFF_lib.cpp 
init_buffers
init_buffers_UFF
print_debugs 
setSwitches2
setSwitchesUFF

---

# GPT-5.1-Codex

Here’s a consolidated plan to unify debugging/visualization into a single `pyBall/FFutils.py`, covering CPU/C++ + C++/OpenCL (`MMFF_multi_lib`) + pure pyOpenCL, with safe buffer views and consistent debug formatting.

### A) Safe buffer access from C++ (no double-free)
- Use existing ctypes helpers in [pyBall/MMFF_multi.py](cci:7://file:///home/prokop/git/FireCore/pyBall/MMFF_multi.py:0:0-0:0): `getIBuff/getBuff/getfBuff` wrap raw C++ pointers without transferring ownership. Keep references to the backing ctypes pointer (and to the loaded `lib`) in the Python object to prevent GC from unloading the library.
- Add thin helpers in `FFutils.py`:
  - `as_np_int(ptr, shape)`, `as_np_float(ptr, shape)`, `as_np_double(ptr, shape)` that **only** wrap ctypes pointers; no copies, no frees.
  - A small guard that refuses to `del`/`free`; document “C++ owns memory”.
- For pyOpenCL, add symmetric helpers to download CL buffers to numpy for diff, but keep the above for zero-copy CPU/GPU-host buffers.

### B) Uniform debug print format (CPU & kernels & pyOpenCL)
- Standard lines (fixed width, one line per contribution):
  - `AT sys=%3d i=%4d pos=(% .6f,% .6f,% .6f) f=(% .6f,% .6f,% .6f)`
  - `BOND sys=%3d i=%4d j=%4d k=% .6f l0=% .6f l=% .6f fb=(% .6f,% .6f,% .6f)`
  - `ANGL sys=%3d i=%4d j=%4d k=%4d K=% .6f c0h=(% .6f,% .6f) f1=(% .6f,% .6f,% .6f) f2=(% .6f,% .6f,% .6f)`
  - `PPI sys=%3d i=%4d j=%4d kpp=% .6f dot=% .6f fpi=(% .6f,% .6f,% .6f) fpj=(% .6f,% .6f,% .6f)`
  - `NB sys=%3d i=%4d j=%4d REQ=(% .6f,% .6f,% .6f) r=% .6f f=(% .6f,% .6f,% .6f)`
- Gate by `idebug/verbosity` and by atom list + iteration modulus. Wire C++ side via [setVerbosity](cci:1://file:///home/prokop/git/FireCore/pyBall/MMFF_multi.py:420:0-421:48) / [print_debugs](cci:1://file:///home/prokop/git/FireCore/pyBall/MMFF_multi.py:408:0-410:62) / [setSwitches2](cci:1://file:///home/prokop/git/FireCore/cpp/libs/Molecular/MMFF_lib.cpp:361:0-381:1) / [setSwitchesUFF](cci:1://file:///home/prokop/git/FireCore/pyBall/MMFF_multi.py:472:0-473:120) (see [MMFF_lib.cpp](cci:7://file:///home/prokop/git/FireCore/cpp/libs/Molecular/MMFF_lib.cpp:0:0-0:0), `MMFFmulti_lib.cpp`); mirror same gates in pyOpenCL host wrappers.

### C) Visualization consolidation (FFutils.py)
- Provide two lightweight classes, reusing existing code:
  - [LiveViz2D](cci:2://file:///home/prokop/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:30:0-178:29) (from [XPBD_2D_utils.py](cci:7://file:///home/prokop/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:0:0-0:0)) for 2D projections.
  - `LiveViz3D` (rename/wrap `LivePortViz` from `XPTB_utils.py`) for 3D.
- Common API:
  - [update(pos, bonds=None, ports=None, forces=None, labels=None, title=None, info=None)](cci:1://file:///home/prokop/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:86:4-178:29)
  - `set_picker(on_pick=None, on_drag=None, on_release=None, radius_px=8)` returning a `pick` dict `{active, idx, mouse, world}`; no sim coupling.
  - Optional overlays: atom indices, bond indices (for matching debug logs).
- Provide “adapter” helpers:
  - `make_bonds_from_neighs(neighs)` to turn adjacency into unique edge list.
  - `ports_from_pi(pipos, scale=1.0)` to visualize pi vectors; can also show `hneigh` for UFF.
  - `plot_state_with_ports_2d/3d` single-call convenience.
- Keep Matplotlib fast path: blit-based animation (as in existing LiveViz2D), minimal redraw.

### D) Topology/geometry helpers to share
- Permutations/reordering: unify `invert_permutation`, `apply_permutation_to_bonds`, `reorder_nodes_first` (2D/3D masks).
- Neigh/back-neigh builders from bond lists; bkSlots/stiffness map.
- Port builders: from positions (2D/3D) or from stored pipos/quats.
- Rotation utilities: `rot_from_angles_2d`, `quat_mul/normalize/rotate`, perturbations.
- IO: `load_xyz`, `write_xyz_with_ports`, `write_pdb_trajectory`, `ensure_outdir`, `print_run_header`.
- Diagnostics: [compute_port_error](cci:1://file:///home/prokop/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:365:0-385:23) (2D/3D), NaN checks, momentum checks.

### E) Parity workflow (fits into FFutils)
1) **Init both stacks** with same XYZ and params.
2) **Wrap C++ buffers** via FFutils safe views; gather pyOpenCL arrays.
3) **Static diff**: topology/params (`neighs`, `REQs`, `apars/bLs/bKs/Ksp/Kpp` or UFF `bon/ang/dih/inv`).
4) **Single-step force diff**: run one pass (no integrator), compare `fapos`; if mismatch, enable gated debug prints with the uniform format.
5) **Multi-step drift**: run small N; compare norms/energies.
6) **Visual spot check**: feed positions/bonds/pi/forces into `LiveViz2D/3D`; enable index overlays for matching printed topology.

### F) Action items to implement in code (when in Code mode)
- Add `pyBall/FFutils.py` with:
  - Zero-copy numpy wrappers for C++ buffers.
  - Shared permutation/topology/port builders.
  - Unified visualizers (2D/3D) + pickers.
  - Debug-print format constants and helper to format lines identically across hosts.
- Add thin shims in [XPBD_2D_utils.py](cci:7://file:///home/prokop/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:0:0-0:0) / `XPTB_utils.py` importing from FFutils to maintain backward compatibility.
- Patch pyOpenCL UFF path to ensure `assembleForces_UFF` is called and `fint` sized per force-piece.
- Add index overlays option to visualizers (atom/bond labels).
- Optional: small parity runner script calling both `MMFF_multi` and pyOpenCL, dumping npz for diffs.

If you want, I can draft the `FFutils.py` skeleton (imports, safe buffer wrappers, shared viz classes, adapters) and show how to swap the existing tests to use it.

---

# Claude Sonet 4.5 Thinking

# FFutils.py Consolidation Design

Based on my analysis of the buffer wrapping utilities, visual debugging tools, and the memory safety concerns you raised, here's a comprehensive design for consolidating everything into a unified `FFutils.py` module.

## 1) Memory-Safe Buffer Wrapping

### Current State Analysis

**From `@/home/prokop/.windsurf/worktrees/FireCore/FireCore-1b118e0e/cpp/common/libUtils.h`:**
- C++ maintains global hashmaps: `buffers` (double*), `fbuffers` (float*), `ibuffers` (int*), `bbuffers` (bool*)
- [getBuff(name)](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-1b118e0e/pyBall/MMFF_multi.py:196:0-200:48) / [getIBuff(name)](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-1b118e0e/pyBall/MMFF.py:794:0-798:48) / [getfBuff(name)](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-1b118e0e/pyBall/MMFF_multi.py:214:0-218:49) return **raw pointers**

**From `@/home/prokop/.windsurf/worktrees/FireCore/FireCore-1b118e0e/pyBall/MMFF_multi.py:192`:**
```python
def getBuff(name, sh):
    ptr = lib.getBuff(name)
    return np.ctypeslib.as_array(ptr, shape=sh)
```

**Critical issue:** `np.ctypeslib.as_array()` creates a NumPy array that **views** C++ memory without copying. The returned array does NOT own the memory, so:
1. ✅ **Safe**: If C++ buffer lifetime ≥ Python session (which is true for [MolWorld_sp3](cci:2://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-1b118e0e/cpp/common/molecular/MolWorld_sp3.h:76:0-3258:1) globals)
2. ⚠️ **Unsafe**: If Python GC tries to free the array (shouldn't happen with views, but depends on NumPy version)
3. ❌ **Double-free**: If same buffer is wrapped multiple times and one view is freed

### Proposed Solution: `BufferRegistry` Class

```python
# pyBall/FFutils.py

import numpy as np
import ctypes
from typing import Dict, Tuple, Optional, Literal
import weakref

class BufferRegistry:
    """
    Thread-safe registry for C++-backed NumPy arrays with proper lifetime management.
    
    Memory safety rules:
    - Buffers are VIEWS into C++ memory (no copy)
    - C++ owns the memory (must outlive Python references)
    - Registry keeps weak references to prevent double-wrapping
    - User must ensure C++ buffer lifetime (e.g., don't call dealloc while holding Python refs)
    """
    
    def __init__(self, lib):
        self.lib = lib
        self._cache: Dict[str, weakref.ref] = {}  # name -> weakref to numpy array
        self._shapes: Dict[str, Tuple] = {}       # name -> shape (for re-wrapping if cache expired)
        
    def get_buffer(self, name: str, shape: Tuple[int, ...], dtype: Literal['double','float','int','bool'] = 'double') -> np.ndarray:
        """
        Get NumPy view of C++ buffer. Returns cached array if still alive, else re-wraps.
        
        Args:
            name: Buffer name registered via init_buffers()
            shape: NumPy shape tuple
            dtype: 'double', 'float', 'int', or 'bool'
        
        Returns:
            NumPy array VIEW (no copy, C++ owns memory)
        """
        # Check cache
        if name in self._cache:
            arr = self._cache[name]()  # Dereference weakref
            if arr is not None:
                return arr
        
        # Wrap C++ pointer
        if dtype == 'double':
            ptr = self.lib.getBuff(name.encode('utf8'))
            arr = np.ctypeslib.as_array(ptr, shape=shape)
        elif dtype == 'float':
            ptr = self.lib.getfBuff(name.encode('utf8'))
            arr = np.ctypeslib.as_array(ptr, shape=shape)
        elif dtype == 'int':
            ptr = self.lib.getIBuff(name.encode('utf8'))
            arr = np.ctypeslib.as_array(ptr, shape=shape)
        elif dtype == 'bool':
            ptr = self.lib.getBBuff(name.encode('utf8'))
            arr = np.ctypeslib.as_array(ptr, shape=shape)
        else:
            raise ValueError(f"Unknown dtype: {dtype}")
        
        if arr is None or arr.ctypes.data == 0:
            raise ValueError(f"Buffer '{name}' not found or NULL pointer")
        
        # Cache + store shape for potential re-wrap
        self._cache[name] = weakref.ref(arr)
        self._shapes[name] = (shape, dtype)
        
        # Mark array as non-writeable IF user wants safety (optional)
        # arr.flags.writeable = False
        
        return arr
    
    def clear_cache(self):
        """Clear cache (useful after C++ realloc)."""
        self._cache.clear()
```

**Usage:**
```python
import pyBall.MMFF_multi as uff
buffers = FFutils.BufferRegistry(uff.lib)

# First call: wraps C++ pointer
apos = buffers.get_buffer('apos', (natoms, 3), dtype='double')

# Second call: returns cached array (same object)
apos2 = buffers.get_buffer('apos', (natoms, 3), dtype='double')
assert apos is apos2  # Same object

# After C++ realloc, clear cache
buffers.clear_cache()
```

## 2) Unified Visualization Interface

### Goals
- Single class for 2D/3D molecular visualization
- BLIT-based matplotlib for performance (avoid full canvas redraw)
- Mouse picking with raycasting (3D) or distance-based (2D)
- Atom/bond index labels
- Support for ports, pi-orbitals, force quivers

### Class Hierarchy

```python
class MolViz:
    """Base class for molecular visualization."""
    def __init__(self, elems, dim=2, figsize=(10,8), show_labels=True):
        self.dim = dim
        self.elems = elems
        self.show_labels = show_labels
        self._setup_figure(figsize)
        self._picker = None
        
    def update(self, pos, **kwargs):
        """Update positions and redraw. Subclass implements."""
        raise NotImplementedError
        
    def attach_picker(self, sim=None, **kwargs):
        """Attach mouse picking. Returns pick dict."""
        raise NotImplementedError

class MolViz2D(MolViz):
    """2D molecular visualization with BLIT."""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, dim=2, **kwargs)
        self._setup_artists()
        
    def _setup_artists(self):
        # Scatter for atoms
        self.scatter_atoms = self.ax.scatter([], [], s=500, ...)
        self.scatter_atoms.set_animated(True)  # Enable BLIT
        
        # LineCollection for bonds
        self.lc_bonds = LineCollection([], colors='g', linewidths=1.5)
        self.lc_bonds.set_animated(True)
        self.ax.add_collection(self.lc_bonds)
        
        # Text labels for atom indices
        self._atom_labels = []
        
        # Background cache for BLIT
        self._bg = None
        
    def update(self, pos, bonds=None, neighs=None, forces=None, **kwargs):
        # Clear old labels
        for lbl in self._atom_labels:
            lbl.remove()
        self._atom_labels.clear()
        
        # Update atom positions
        self.scatter_atoms.set_offsets(pos[:, :2])
        
        # Update bonds
        if bonds is not None:
            segs = [[pos[i], pos[j]] for (i,j) in bonds]
            self.lc_bonds.set_segments(segs)
        elif neighs is not None:
            segs = []
            for i, neigh_row in enumerate(neighs):
                for j in neigh_row:
                    if j >= 0:
                        segs.append([pos[i], pos[j]])
            self.lc_bonds.set_segments(segs)
        
        # Draw atom index labels
        if self.show_labels:
            for i, (x, y) in enumerate(pos[:, :2]):
                lbl = self.ax.text(x, y, str(i), fontsize=9, ha='center', va='center', zorder=10)
                lbl.set_animated(True)
                self._atom_labels.append(lbl)
        
        # BLIT update
        if self._bg is None:
            self._bg = self.fig.canvas.copy_from_bbox(self.ax.bbox)
        
        self.fig.canvas.restore_region(self._bg)
        self.ax.draw_artist(self.scatter_atoms)
        self.ax.draw_artist(self.lc_bonds)
        for lbl in self._atom_labels:
            self.ax.draw_artist(lbl)
        self.fig.canvas.blit(self.ax.bbox)
        self.fig.canvas.flush_events()

class MolViz3D(MolViz):
    """3D molecular visualization with raycasting picker."""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, dim=3, **kwargs)
        self.ax = self.fig.add_subplot(111, projection='3d')
        self._setup_artists()
        
    def _project_to_screen(self, pos3d):
        """Project 3D positions to 2D screen coords + depth."""
        from mpl_toolkits.mplot3d import proj3d
        M = self.ax.get_proj()
        x2, y2, z2 = proj3d.proj_transform(pos3d[:, 0], pos3d[:, 1], pos3d[:, 2], M)
        pts2 = np.vstack([x2, y2]).T
        pts_px = self.ax.transData.transform(pts2)
        return pts_px, z2
```

**Consolidation from existing tools:**
- Reuse [LiveViz2D](cci:2://file:///home/prokop/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:30:0-178:29) BLIT logic from `@/home/prokop/.windsurf/worktrees/FireCore/FireCore-1b118e0e/pyBall/XPBD_2D/XPBD_2D_utils.py:31`
- Reuse 3D projection from `@/home/prokop/.windsurf/worktrees/FireCore/FireCore-1b118e0e/pyBall/XPDB_AVBD/XPTB_utils.py:481`
- **Avoid duplicating** these classes; import from FFutils

## 3) Buffer Comparison Utilities

### Unified Comparison Function

```python
def compare_buffers(buf1: dict, buf2: dict, tol=1e-6, max_print=10, label1='CPU', label2='GPU') -> bool:
    """
    Compare two buffer dictionaries with detailed diff reporting.
    
    Args:
        buf1, buf2: Dict[str, np.ndarray]
        tol: Absolute tolerance for floating-point
        max_print: Max number of mismatched elements to print per buffer
        label1, label2: Labels for the two implementations
    
    Returns:
        True if all buffers match within tolerance
    """
    all_match = True
    keys1 = set(buf1.keys())
    keys2 = set(buf2.keys())
    
    if keys1 != keys2:
        print(f"[ERROR] Buffer name mismatch:")
        print(f"  {label1} only: {keys1 - keys2}")
        print(f"  {label2} only: {keys2 - keys1}")
        all_match = False
    
    for name in sorted(keys1 & keys2):
        a = buf1[name]
        b = buf2[name]
        
        if a.shape != b.shape:
            print(f"[ERROR] {name}: shape mismatch {label1}={a.shape} {label2}={b.shape}")
            all_match = False
            continue
        
        if np.issubdtype(a.dtype, np.floating):
            diff = np.abs(a - b)
            max_diff = np.max(diff)
            rms_diff = np.sqrt(np.mean(diff**2))
            
            if max_diff > tol:
                print(f"[DIFF] {name}: max|Δ|={max_diff:.3e} rms|Δ|={rms_diff:.3e}")
                # Find worst offenders
                idx_sorted = np.argsort(diff.ravel())[::-1]
                for k in range(min(max_print, len(idx_sorted))):
                    idx_flat = idx_sorted[k]
                    idx = np.unravel_index(idx_flat, a.shape)
                    print(f"  idx={idx}: {label1}={a[idx]:.6f} {label2}={b[idx]:.6f} Δ={diff[idx]:.3e}")
                all_match = False
        else:  # Integer or bool
            mask = (a != b)
            n_diff = np.sum(mask)
            if n_diff > 0:
                print(f"[DIFF] {name}: {n_diff}/{a.size} elements differ")
                idx_diff = np.nonzero(mask.ravel())[0]
                for k in range(min(max_print, len(idx_diff))):
                    idx_flat = idx_diff[k]
                    idx = np.unravel_index(idx_flat, a.shape)
                    print(f"  idx={idx}: {label1}={a[idx]} {label2}={b[idx]}")
                all_match = False
    
    if all_match:
        print(f"[OK] All buffers match within tol={tol:.0e}")
    
    return all_match
```

## 4) Debug Print Formatting

### Gated Debug Prints (Python Side)

```python
class DebugPrinter:
    """Formatted debug printing for force evaluation."""
    
    def __init__(self, enabled=False, atom_filter=None, iter_filter=None, components=None):
        self.enabled = enabled
        self.atom_filter = set(atom_filter) if atom_filter else None
        self.iter_filter = set(iter_filter) if iter_filter else None
        self.components = set(components) if components else {'bond','angle','dih','inv','nb'}
        
    def should_print(self, ia, iter_idx, component):
        if not self.enabled:
            return False
        if self.atom_filter and ia not in self.atom_filter:
            return False
        if self.iter_filter and iter_idx not in self.iter_filter:
            return False
        if component not in self.components:
            return False
        return True
    
    def print_bond(self, ia, jn, r, l0, k, f, iter_idx=0):
        if self.should_print(ia, iter_idx, 'bond'):
            # Fixed-width format for easy diff
            print(f"UFF_BOND  ia={ia:3d} jn={jn:3d}  r={r:12.6f} l0={l0:12.6f} k={k:12.6f}  f=({f[0]:12.6f},{f[1]:12.6f},{f[2]:12.6f})")
    
    def print_angle(self, ia, jn1, jn2, cos_ang, cos0, k, f, iter_idx=0):
        if self.should_print(ia, iter_idx, 'angle'):
            print(f"UFF_ANGLE ia={ia:3d} jn1={jn1:3d} jn2={jn2:3d}  cos={cos_ang:12.6f} cos0={cos0:12.6f} k={k:12.6f}  f=({f[0]:12.6f},{f[1]:12.6f},{f[2]:12.6f})")
    
    # ... similar for dihedral, inversion, nonbonded
```

**C++ side** (in `@/home/prokop/.windsurf/worktrees/FireCore/FireCore-1b118e0e/cpp/common/molecular/UFF.h`):
```cpp
// Add at top of file
extern int g_debug_level;
extern int g_debug_atom;
extern int g_debug_iter;

// In evalBond():
if (g_debug_level > 0 && ia == g_debug_atom && iter == g_debug_iter) {
    printf("UFF_BOND  ia=%3d jn=%3d  r=%12.6f l0=%12.6f k=%12.6f  f=(%12.6f,%12.6f,%12.6f)\n",
           ia, jn, r, l0, k, f.x, f.y, f.z);
}
```

**Python interface:**
```python
def set_debug_verbosity(level, atom_idx=None, iter_idx=None):
    """Set C++ debug verbosity (for CPU/C++/OpenCL paths)."""
    uff.lib.setDebugVerbosity(level, atom_idx or -1, iter_idx or 0)
```

## 5) Mouse Picking with Raycasting

### Unified Picker Interface

```python
def attach_picker(viz: MolViz, sim=None, pick_radius=0.5, pick_radius_px=20, verbose=0):
    """
    Attach mouse picking to visualization.
    
    For 2D: distance-based picking in data coords
    For 3D: raycasting with projection/unprojection
    
    Returns:
        pick dict with keys: 'idx', 'active', 'mouse' (2D) or 'mouse3' (3D)
    """
    if viz.dim == 2:
        return _attach_picker_2d(viz, sim, pick_radius, verbose)
    elif viz.dim == 3:
        return _attach_picker_3d(viz, sim, pick_radius_px, verbose)
    else:
        raise ValueError(f"Unsupported dimension: {viz.dim}")

def _attach_picker_3d(viz, sim, pick_radius_px, verbose):
    """3D picking with raycasting (reuse from XPTB_utils.py:471)."""
    pick = {"idx": None, "active": False, "mouse3": np.array([0,0,0], dtype=np.float32), "_zproj": None}
    
    def on_press(event):
        if event.inaxes != viz.ax or viz._last_pos is None:
            return
        
        # Project 3D positions to 2D screen
        pts_px, z2 = viz._project_to_screen(viz._last_pos)
        mouse_px = np.array([event.x, event.y], dtype=np.float32)
        
        # Find nearest atom within pick_radius_px
        d2 = np.sum((pts_px - mouse_px[None, :]) ** 2, axis=1)
        i_min = int(np.argmin(d2))
        
        if d2[i_min] <= float(pick_radius_px) ** 2:
            pick["idx"] = i_min
            pick["active"] = True
            pick["_zproj"] = float(z2[i_min])
            pick["mouse3"] = viz._last_pos[i_min].copy()
            if verbose > 0:
                print(f"[PICK3D] press idx={i_min} dpx={np.sqrt(d2[i_min]):.2f}")
    
    def on_motion(event):
        if not pick["active"] or pick["_zproj"] is None:
            return
        
        # Unproject mouse to 3D at fixed depth
        from mpl_toolkits.mplot3d import proj3d
        mouse_px = np.array([event.x, event.y], dtype=np.float32)
        xy_proj = viz.ax.transData.inverted().transform(mouse_px)
        x, y, z = proj3d.inv_transform(xy_proj[0], xy_proj[1], pick["_zproj"], viz.ax.get_proj())
        pick["mouse3"] = np.array([x, y, z], dtype=np.float32)
    
    viz.fig.canvas.mpl_connect('button_press_event', on_press)
    viz.fig.canvas.mpl_connect('motion_notify_event', on_motion)
    return pick
```

## 6) Proposed Module Structure

```
pyBall/
├── FFutils.py                    # NEW: Consolidated utilities
│   ├── BufferRegistry           # Memory-safe C++ buffer wrapping
│   ├── MolViz / MolViz2D / MolViz3D  # Unified visualization
│   ├── compare_buffers()        # Buffer comparison with detailed diff
│   ├── DebugPrinter             # Formatted debug printing
│   ├── attach_picker()          # Unified 2D/3D mouse picking
│   └── format_topology()        # Pretty-print topology buffers
├── MMFF.py                      # Keep for legacy, import from FFutils
├── MMFF_multi.py                # Keep for legacy, import from FFutils
├── XPBD_2D/
│   └── XPBD_2D_utils.py         # Deprecate LiveViz2D, import from FFutils
└── XPDB_AVBD/
    └── XPTB_utils.py             # Deprecate LivePortViz, import from FFutils
```

**Migration strategy:**
1. Create `FFutils.py` with all new classes
2. Update [XPBD_2D_utils.py](cci:7://file:///home/prokop/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:0:0-0:0) to import: `from pyBall.FFutils import MolViz2D as LiveViz2D`
3. Update [XPTB_utils.py](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-1b118e0e/pyBall/XPDB_AVBD/XPTB_utils.py:0:0-0:0) to import: `from pyBall.FFutils import MolViz3D as LivePortViz`
4. Update test scripts to use `FFutils` directly

## 7) Implementation Adapters

For the three implementation paths, we need adapters:

```python
class FFInterface:
    """Abstract interface for forcefield implementations."""
    def get_positions(self) -> np.ndarray:
        raise NotImplementedError
    
    def get_forces(self) -> np.ndarray:
        raise NotImplementedError
    
    def get_topology_buffers(self) -> dict:
        raise NotImplementedError
    
    def eval_step(self):
        raise NotImplementedError

class CPP_UFF(FFInterface):
    """Adapter for C++ UFF (via MMFF_multi or MMFF)."""
    def __init__(self, lib_module):
        self.lib = lib_module
        self.buffers = BufferRegistry(lib_module.lib)
        
    def get_positions(self):
        return self.buffers.get_buffer('apos', (self.lib.natoms, 3), 'double')
    
    # ...

class PyOCL_UFF(FFInterface):
    """Adapter for pyOpenCL UFF."""
    def __init__(self, uff_cl):
        self.uff_cl = uff_cl
        
    def get_positions(self):
        return self.uff_cl.download_buf('apos')  # Already NumPy, pyOpenCL manages memory
    
    # ...
```

## Summary

The consolidation achieves:
✅ **Memory safety**: `BufferRegistry` with weakref caching prevents double-wrapping  
✅ **No code duplication**: Single `MolViz` base class for 2D/3D  
✅ **BLIT performance**: Animated artists + background caching  
✅ **Uniform debug prints**: `DebugPrinter` with consistent formatting  
✅ **Robust mouse picking**: Raycasting (3D) + distance-based (2D)  
✅ **Atom/bond labels**: Integrated into [MolViz.update()](cci:1://file:///home/prokop/git/FireCore/pyBall/XPBD_2D/XPBD_2D_utils.py:86:4-178:29)  
✅ **Buffer comparison**: Detailed diff reporting with tolerances  
✅ **Cross-implementation**: Adapters for pyOpenCL, C++/OpenCL, CPU/C++  

The key insight on memory management: **C++ owns the data, Python creates views**. As long as the C++ [MolWorld_sp3](cci:2://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-1b118e0e/cpp/common/molecular/MolWorld_sp3.h:76:0-3258:1) object outlives the Python session, this is safe. The `BufferRegistry` prevents accidentally creating multiple wrappers of the same pointer.

