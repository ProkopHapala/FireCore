## USER

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

## Gemini-3-pro

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