### XPDB WebGPU – Current Status Report

#### Repo + Files in Play
- **WebGUI / WebGPU implementation**
  - `web/molgui_webgpu/main.js` – MolGUIApp (XPDB init, UI wiring, stepping, logging, disable toggles)
  - `web/molgui_webgpu/GUI.js` – sidebar controls (XPDB enable, angle toggle, dump buttons, pause, verbosity, disable bonds/collisions)
  - `web/molgui_webgpu/XPDBTopology.js` – builds bond adjacency (real vs angle-derived, merged)
  - `web/molgui_webgpu/XPDB_WebGPU.js` – buffer setup, bond uploads, compute dispatch
  - `web/molgui_webgpu/xpdb.wgsl` – compute shaders (bbox, local topology, Jacobi solve)
  - `web/molgui_webgpu/MMParams.js` – parameter tables, verbose lookup logging
  - `web/molgui_webgpu/dump_xpdb_topology.mjs` – Node script to dump topology arrays
  - `web/molgui_webgpu/XPDB_webgpu_status.md` – (empty placeholder; status provided below in chat as requested)

- **Reference / comparison**
  - `pyBall/XPDB_AVBD/XPDB_new.py` + `.cl` – Python/OpenCL solver (correct behavior reference)
  - `pyBall/XPDB_AVBD/test_TiledJacobi_molecules.py` – prints reference buffers (`OUT-test_TiledJacobi_molecules.txt`)
  - `cpp/common_resources/cl/relax_multi_mini.cl` – legacy OpenCL kernels (for bbox / topology logic)

#### What Works Now
1. **No double angle addition**: [XPDBTopology.buildXPDBTopology](cci:1://file:///home/prokophapala/git/FireCore/web/molgui_webgpu/XPDBTopology.js:2:0-201:1) now builds first-order bonds and angle-derived bonds separately, merges them once, and keeps original neighbor lists immutable.
2. **MMParams-driven L0/K** with verbose logging:
   - [getBondL0](cci:1://file:///home/prokophapala/git/FireCore/web/molgui_webgpu/MMParams.js:371:4-389:5) / [getAngleParams](cci:1://file:///home/prokophapala/git/FireCore/web/molgui_webgpu/MMParams.js:411:4-426:5) log every lookup (`VERBOSITY_LEVEL >= 3`) and throw if a pair/triple is missing (no silent fallback).
3. **Topology + buffer dumps**:
   - GUI buttons “Dump topology” and “Dump buffers” print per-atom adjacency (with `i: (count) j:L/K`) and GPU positions.
   - Node script [dump_xpdb_topology.mjs](cci:7://file:///home/prokophapala/git/FireCore/web/molgui_webgpu/dump_xpdb_topology.mjs:0:0-0:0) produces `OUT-dump_xpdb_topology.txt` for 1:1 comparison with the OpenCL reference.
4. **Bond length logging**:
   - After every “Step once” / “Relax N steps”, [logRealBondLengths](cci:1://file:///home/prokophapala/git/FireCore/web/molgui_webgpu/main.js:650:4-679:5) computes actual distances for real bonds from [EditableMolecule](cci:2://file:///home/prokophapala/git/FireCore/web/molgui_webgpu/EditableMolecule.js:240:0-1055:1) (no angle-derived bonds) and prints min/avg/max/sample.
5. **Interactive controls**:
   - Pause toggle, Step once, Relax N steps, Verbosity level, Disable bonds/collisions toggles (UI wiring is in place).
6. **Partial-group collision guard (new)**:
   - `xpdb.wgsl` now zeros `l_rad` for inactive workgroup lanes (`num_atoms < GROUP_SIZE`), eliminating garbage radii. Verified the shader compiles + runs, but this alone did **not** resolve the bonded atoms colliding issue (see "Still Problematic").

#### Still Problematic / Not Working
1. **Bond L0/K still off vs reference (sp2 vs sp3 mix-up)**:
   - Pentacene carbons should be `sp2` (∠≈120°, L0≈1.40 Å). Our dumps show `sp3` parameters (∠≈109°, L0≈1.53 Å), causing tetrahedral geometry. Need to verify which atom-type names are being assigned during MMParams lookup and whether the browser loads the same `.dat` set as the OpenCL reference. Without fixing this, WebGPU cannot converge to the planar reference structure.
2. **Bond disable toggle incomplete**:
   - GUI + MolGUIApp set a flag and call `xpdb.setBondStiffnessScale`, but [XPDB_WebGPU](cci:2://file:///home/prokophapala/git/FireCore/web/molgui_webgpu/XPDB_WebGPU.js:0:0-356:1) does not yet store a host copy of `bondLenStiff` or implement `setBondStiffnessScale`. Disabling bonds currently has no effect.
3. **Collision rejection still broken (root cause unresolved)**:
   - After fixing uninitialized `l_rad`, bonded atoms continue to collide. This implies `bond_indices_local` lacks some bonded neighbors. Likely culprits: (a) `build_local_topology` remap bug for shared workgroups, (b) ghost list truncation when `MAX_GHOSTS` fills, or (c) mismatch vs OpenCL search order. Need to dump and compare the remapped buffer with the Python/OpenCL reference output.
4. **Rendering after relax**:
   - Step/Relax now calls `renderers.molecule.update()` and [requestRender()](cci:1://file:///home/prokophapala/git/FireCore/web/molgui_webgpu/GUI.js:36:4-38:5), so geometry refresh works. (Previously missing; fixed.)

#### Remaining TODO / Next Steps
1. **Fix bond parameter source / atom typing**:
   - Confirm the `.dat` files loaded in-browser match the OpenCL reference (`ElementTypes`, `AtomTypes`, `BondTypes`, `AngleTypes`).
   - Add diagnostics to `XPDBTopology` / `MMParams` to log the atom-type name chosen per atom (expect `C_sp2` for pentacene). Cross-check with `test_TiledJacobi_molecules.py` output to ensure identical naming and parameter IDs.
   - Once WebGPU reports L0/K exactly matching `OUT-test_TiledJacobi_molecules.txt`, re-run dumps for parity.
2. **Implement `XPDB_WebGPU.setBondStiffnessScale(scale)`**:
   - Keep a host-side `Float32Array` copy of `bondLenStiff` after upload.
   - Method multiplies every stiffness (y component) by `scale` and writes updated buffer (or uses compute shader to scale).
   - Used when “Disable bonds” toggle is checked.
3. **Systematic remap debug (`bond_indices_local`)**:
   - After each `build_local_topology` dispatch, enqueue `device.queue.readBuffer` for `bond_indices_local` (and `ghost_packed`).
   - Post-process on CPU to produce lines of the form `atom 12 (global 77): slot0 -> local 5 (global 80)`, `slot1 -> ghost 3 (global 142)`, etc.
   - Compare with OpenCL’s debug output (from `test_TiledJacobi_molecules.py`) to confirm every bonded neighbor is accounted for. If differences exist, inspect:
     1. Off-by-one errors when `num_atoms` not divisible by `GROUP_SIZE`.
     2. Ghost overflow (when `MAX_GHOSTS`=128 fills, entries at the tail vanish — look for warnings and count mismatches).
     3. Divergent search order (WGSL vs OpenCL) causing early termination when `bond_indices_global` contains `-1` mid-array.
4. **Optional debugging kernel**:
   - Temporarily add a WGSL path that short-circuits reindexing: directly read global positions and skip collisions when `bond_indices_global` says so. Use this to confirm collision behavior matches reference when remap is bypassed.

#### Outstanding Challenges
- **Parameter parity**: Need to ensure WebGPU side uses the same MM force-field parameters as Python/OpenCL (atom types vs Z, same `.dat` files).
- **Remapping correctness**: When group ghost buffers overflow (MAX_GHOSTS=128), bonds to atoms outside the local+ghost set disappear, reintroducing bonded collisions.
- **UI/UX**: more logging/dump buttons add console noise; guard with verbosity levels.

#### Session Update – 2026-01-24 (Raw WebGPU bonds/halos + label investigation)
**What we accomplished today**
1. **Raw WebGPU bond + halo rendering validated**
   - `MolGUIApp.requestRender()` now uploads bond indices and selection sets every frame even when running the pure WebGPU backend; visually confirmed that cylinders (bonds) and halos follow camera + XPDB updates.
   - The selection halo UBO is refreshed through `RawWebGPUAtomsRenderer.updateSelection(selIdx, n, 1.3)` (selection radius scale) and the pipeline binding logs prove the buffers remain in sync after XPDB relax steps.
2. **Automatic startup molecule + labels**
   - `installMoleculeIOMethods(EditableMolecule)` is invoked in `main.js`, `_autoLoadDefaultMolecule()` fetches `../common_resources/mol/H2O.mol2`, clears the system, appends parsed atoms/bonds, forces Atom-ID label mode, and calls `raw.setLabelsVisible(true)` so debugging starts immediately without GUI clicks.
3. **Font-atlas instrumentation**
   - `_ensureFontAtlas()` now probes CPU + GPU copies and `_refreshLabelBindGroup()` recreates the label bind group whenever the atlas texture/sampler becomes available; this removed stale texture bindings as a possible root cause.

**Still struggling (labels / texture atlas)**
1. **Label fragments still show UV gradients instead of glyphs**
   - Despite the atlas texture sampling correctly in the debug overlay pipeline, the label fragment shader `fs_labels` returns only the fallback gradient. Bind group snapshots confirm buffers + sampler bindings are valid at draw time.
2. **Probable causes**
   1. The label pipeline might still reference an old texture view because the bind group was created before `_ensureFontAtlas()` finished (mitigated by `_refreshLabelBindGroup`, but needs confirmation via GPU readback).
   2. WGSL UV math could be mis-indexing the atlas (e.g., using normalized indices before multiplying by glyph size), resulting in coordinates that land outside populated texels and sample pure zeros.
   3. Less likely but possible: Dawn validation warning indicates swapchain texture usage flags; if Dawn silently drops the label bind group when the warning fires, the shader could end up sampling an uninitialized default texture.

3. **Mouse picking / drag-box selection inaccurate**
   - The editor consistently selects atoms that are offset from the cursor and marquee rectangles; likely caused by a mismatch between the orthographic camera transforms used for screen-space ray construction and the Raw WebGPU render path (camera matrices live in `MolGUIApp` while selection logic still assumes the older Three.js projection). Needs a full audit of `Editor.pickRayFromScreen()` vs the matrices passed into `RawWebGPUAtomsRenderer.updateCamera()` to ensure identical view/projection transforms and consistent handling of rotations.

**Next debugging steps for labels**
1. **Add an on-demand sampler probe** that reuses the label bind group/sampler and writes sampled RGBA values for a fixed glyph into a staging buffer. Compare with the atlas overlay probe to prove whether the fragment shader sees glyph data.
2. **Instrument the vertex shader outputs** (either via storage buffer or debug color encoding) to log the computed glyph UV rectangles per character; verify they match the atlas layout `(cols=16, rows=6, cw/ch derived from atlas dims)`.
3. **Temporarily bypass UV math** by hardcoding UVs for a known glyph (e.g., ASCII '0') inside `fs_labels`. If glyph pixels appear, the issue is our per-character UV calculation; if not, the pipeline binding still references the wrong texture view.
4. **Silence/resolve the Dawn CopyDst warning** by ensuring the swapchain texture is created with both `RENDER_ATTACHMENT | COPY_DST` (if we control creation) or by avoiding `copyTextureToTexture` on it. Eliminating the warning will clarify whether Dawn is invalidating our render pass mid-frame.
5. **Consider isolating labels into their own render pass** with a dedicated encoder + bind group to ensure the atlas sampler state cannot be clobbered by the atlas-debug overlay.

#### References
- WebGPU implementation:
  - [web/molgui_webgpu/main.js](cci:7://file:///home/prokophapala/git/FireCore/web/molgui_webgpu/main.js:0:0-0:0), [GUI.js](cci:7://file:///home/prokophapala/git/FireCore/web/molgui_webgpu/GUI.js:0:0-0:0), [XPDBTopology.js](cci:7://file:///home/prokophapala/git/FireCore/web/molgui_webgpu/XPDBTopology.js:0:0-0:0), [XPDB_WebGPU.js](cci:7://file:///home/prokophapala/git/FireCore/web/molgui_webgpu/XPDB_WebGPU.js:0:0-0:0), [xpdb.wgsl](cci:7://file:///home/prokophapala/git/FireCore/web/molgui_webgpu/xpdb.wgsl:0:0-0:0), [MMParams.js](cci:7://file:///home/prokophapala/git/FireCore/web/molgui_webgpu/MMParams.js:0:0-0:0)
- Debug scripts: [web/molgui_webgpu/dump_xpdb_topology.mjs](cci:7://file:///home/prokophapala/git/FireCore/web/molgui_webgpu/dump_xpdb_topology.mjs:0:0-0:0)
- Reference OpenCL/Python:
  - `pyBall/XPDB_AVBD/XPDB_new.py`, [pyBall/XPDB_AVBD/XPDB_new.cl](cci:7://file:///home/prokophapala/git/FireCore/pyBall/XPDB_AVBD/XPDB_new.cl:0:0-0:0), `pyBall/XPDB_AVBD/test_TiledJacobi_molecules.py`, `pyBall/XPDB_AVBD/OUT-test_TiledJacobi_molecules.txt`
- Legacy kernels: `cpp/common_resources/cl/relax_multi_mini.cl`

Let me know if you want me to start implementing the bond stiffness scaling or the `bond_indices_local` dump next.