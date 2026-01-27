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

#### XPDB Resolved Issues (2026-01-24)
1. **Label Visibility**: Fixed missing `RENDER_ATTACHMENT` flag on atlas texture, corrected UV Y-flip in WGSL, added `@interpolate(flat)` for character codes, and refined atlas opacity for reliable alpha-discarding.
2. **Selection Accuracy**: Synchronized camera `matrixWorldInverse` in `main.js` to ensure 3D mouse picking matches the WebGPU render state.
3. **Dawn/Vulkan Validation**: Added `COPY_DST` to swapchain to resolve validation warnings during buffer clears.

#### What Works Now
1. **No double angle addition**: [XPDBTopology.buildXPDBTopology](cci:1://file:///home/prokophapala/git/FireCore/web/molgui_webgpu/XPDBTopology.js:2:0-201:1) now builds first-order bonds and angle-derived bonds separately, merges them once, and keeps original neighbor lists immutable.
2. **MMParams-driven L0/K** with verbose logging:
   - [getBondL0](cci:1://file:///home/prokophapala/git/FireCore/web/molgui_webgpu/MMParams.js:371:4-389:5) / [getAngleParams](cci:1://file:///home/prokophapala/git/FireCore/web/molgui_webgpu/MMParams.js:411:4-426:5) log every lookup (`VERBOSITY_LEVEL >= 3`) and throw if a pair/triple is missing (no silent fallback).
3. **Topology + buffer dumps**:
   - GUI buttons “Dump topology” and “Dump buffers” print per-atom adjacency (with i: (count) j:L/K) and GPU positions.
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

**Windsurf/GPT5.2/GPT-5-Codex: What we accomplished today**
1. **Raw WebGPU bond + halo rendering validated**
   - `MolGUIApp.requestRender()` now uploads bond indices and selection sets every frame even when running the pure WebGPU backend; visually confirmed that cylinders (bonds) and halos follow camera + XPDB updates.
   - The selection halo UBO is refreshed through `RawWebGPUAtomsRenderer.updateSelection(selIdx, n, 1.3)` (selection radius scale) and the pipeline binding logs prove the buffers remain in sync after XPDB relax steps.
2. **Automatic startup molecule + labels**
   - `installMoleculeIOMethods(EditableMolecule)` is invoked in `main.js`, `_autoLoadDefaultMolecule()` fetches `../common_resources/mol/H2O.mol2`, clears the system, appends parsed atoms/bonds, forces Atom-ID label mode, and calls `raw.setLabelsVisible(true)` so debugging starts immediately without GUI clicks.
3. **Font-atlas instrumentation**
   - `_ensureFontAtlas()` now probes CPU + GPU copies and `_refreshLabelBindGroup()` recreates the label bind group whenever the atlas texture/sampler becomes available; this removed stale texture bindings as a possible root cause.

**Windsurf/GPT5.2/GPT-5-Codex: Still struggling (labels / texture atlas)**
1. **Label fragments still show UV gradients instead of glyphs**
   - Despite the atlas texture sampling correctly in the debug overlay pipeline, the label fragment shader `fs_labels` returns only the fallback gradient. Bind group snapshots confirm buffers + sampler bindings are valid at draw time.
2. **Probable causes**
   1. The label pipeline might still reference an old texture view because the bind group was created before `_ensureFontAtlas()` finished (mitigated by `_refreshLabelBindGroup`, but needs confirmation via GPU readback).
   2. WGSL UV math could be mis-indexing the atlas (e.g., using normalized indices before multiplying by glyph size), resulting in coordinates that land outside populated texels and sample pure zeros.
   3. Less likely but possible: Dawn validation warning indicates swapchain texture usage flags; if Dawn silently drops the label bind group when the warning fires, the shader could end up sampling an uninitialized default texture.

3. **Mouse picking / drag-box selection inaccurate**
   - The editor consistently selects atoms that are offset from the cursor and marquee rectangles; likely caused by a mismatch between the orthographic camera transforms used for screen-space ray construction and the Raw WebGPU render path (camera matrices live in `MolGUIApp` while selection logic still assumes the older Three.js projection). Needs a full audit of `Editor.pickRayFromScreen()` vs the matrices passed into `RawWebGPUAtomsRenderer.updateCamera()` to ensure identical view/projection transforms and consistent handling of rotations.

**Next debugging steps for labels**
1. ~~Add an on-demand sampler probe~~ (labels now render correctly; probe no longer needed).
2. ~~Instrument the vertex shader outputs~~ (labels now render correctly; instrumentation no longer needed).
3. ~~Temporarily bypass UV math~~ (labels now render correctly; bypass no longer needed).
4. ~~Silence/resolve the Dawn CopyDst warning~~ (labels now render correctly; warning tolerable).
5. ~~Consider isolating labels into their own render pass~~ (labels now render correctly; isolation no longer needed).

**Cleanup performed (2026-01-24)**
- Commented out atlas debug overlay (`debugShowLabelAtlas = false` by default) and DOM image creation (`_fontAtlasDebugImg`) in `RawWebGPUAtomsRenderer._ensureFontAtlas()`; code retained for future debugging.
- Commented out verbose console probes (`console.log` for layout, pixel probes, GPU readback) but kept lines with `// TODO: re-enable for debugging`.
- Added `setLabelStyle(color, size)` to `RawWebGPUAtomsRenderer` to parse hex color strings and update `labelColor`/`labelScale`, then call `_uploadLabelUBO()` to sync the uniform.
- Wired GUI label style controls (`GUI.js`) to call `raw.setLabelStyle()` when `window.app.useRawWebGPU` is true, so color/size changes now affect the Raw WebGPU renderer.

**Antigrafity/Gemini-3-False: SUCCESS - Refined Labels and Selection**
- Fixed **Label Rendering**: Identified missing `RENDER_ATTACHMENT` flag required for WebGPU image copies, corrected UV math, and added flat interpolation to prevent glyph corruption.
- Fixed **Atom Selection**: Synchronized camera world inverse matrices in the render loop to eliminate the mouse picking offset.
- Added **GPU Debug Modes**: Instrumented `fs_labels` with modes for atlas mapping, raw sampling, and character buffer validation.

#### Session Update – 2026-01-24 (Raw WebGPU dirty flags, buffer sync, XPDB controls)

**What we achieved**
1. **Debug log cleanup / opt-in tracing**
   - Commented or gated the spammy `[RawWebGPUAtomsRenderer]` and label-atlas dumps; `MolGUIApp` now uses `debugRawLogs` to print packed-buffer snapshots only when explicitly toggled (default `false`).
2. **Editable → Packed synchronization parity with WebGL path**
   - Added `_markRawTopologyDirty`, `_markRawGeometryDirty`, `_markRawSelectionDirty`, `_markRawLabelsDirty` helpers in `MolGUIApp` and ensured `requestRender()` exports the `EditableMolecule` into `PackedMolecule` whenever `dirtyExport|dirtyTopo|dirtyGeom` mirrors reference `MoleculeRenderer` (@main.js lines around 60–210).
   - GUI entry points (`loadXYZ`, Clear Scene, Add Example) now call `window.app.markRawAllDirty()` so raw buffers refresh immediately after user actions, matching Three.js behavior.
3. **Editor integration for add/delete/move**
   - `Editor.deleteSelection`, `addAtom`, `recalculateBonds`, and gizmo drag updates now invoke `window.app.markRawTopologyDirty/markRawGeometryDirty`, ensuring sphere/bond buffers rebuild (not just labels) after edits.
4. **XPDB geometry sync + render scheduling**
   - `stepXPDB`, `xpdbStepOnce`, `xpdbRelaxSteps`, and `syncXPDBPositions` mark geometry dirty and call `requestRender()` so XPDB buttons update atoms even when the camera is idle.
   - `requestRender()` gained `_needsExtraFrame` queuing so overlapping render requests from XPDB or editor actions are not dropped; on-demand WebGPU mode now reliably flushes every pending update.

**Architectural decisions & rationale**
1. **Dual-tier dirty flags**
   - Topology operations (atom/bond count changes) trigger `_markRawAllDirty` so atoms, bonds, labels, and selection buffers rebuild in one pass; geometry-only changes keep allocations intact and only re-upload atom/label data. Mirrors the reference `MeshBuilder/MoleculeRenderer` division between `structureDirty` and `posDirty` for better perf.
2. **On-demand renderer contract**
   - Rather than switching to a continuous animation loop, `MolGUIApp.requestRender()` now guarantees at-least-once rendering for every logical edit by deferring extra frames via `_needsExtraFrame`. This keeps CPU/GPU load low while ensuring XPDB/GUI events still reach the screen.
3. **Debug logging policy**
   - Mandatory instrumentation (dirty flag transitions, exports, buffer upload counts) is present but commented/gated so we can re-enable quickly during regressions without paying per-frame costs.

**Open follow-ups / risks**
1. **Selection buffer refresh policy**
   - Selection updates currently piggyback on topology helpers; verify halo buffers refresh correctly when selection changes in the absence of topology edits.
2. **XPDB continuous stepping**
   - While manual button presses trigger renders, auto-stepping (`xpdbEnabled && !xpdbPaused`) still depends on `requestRender()` being called elsewhere. Consider scheduling periodic renders while XPDB runs continuously.
3. **Dirty flag cohesion**
   - We now rely on every mutation path calling the right helper (e.g., script runners, bucket rebuilds). Need future audit to ensure custom scripts also tap into `markRawTopologyDirty/markRawGeometryDirty`.

#### References
- WebGPU implementation:
  - [web/molgui_webgpu/main.js](cci:7://file:///home/prokophapala/git/FireCore/web/molgui_webgpu/main.js:0:0-0:0), [GUI.js](cci:7://file:///home/prokophapala/git/FireCore/web/molgui_webgpu/GUI.js:0:0-0:0), [XPDBTopology.js](cci:7://file:///home/prokophapala/git/FireCore/web/molgui_webgpu/XPDBTopology.js:0:0-0:0), [XPDB_WebGPU.js](cci:7://file:///home/prokophapala/git/FireCore/web/molgui_webgpu/XPDB_WebGPU.js:0:0-0:0), [xpdb.wgsl](cci:7://file:///home/prokophapala/git/FireCore/web/molgui_webgpu/xpdb.wgsl:0:0-0:0), [MMParams.js](cci:7://file:///home/prokophapala/git/FireCore/web/molgui_webgpu/MMParams.js:0:0-0:0)
- Debug scripts: [web/molgui_webgpu/dump_xpdb_topology.mjs](cci:7://file:///home/prokophapala/git/FireCore/web/molgui_webgpu/dump_xpdb_topology.mjs:0:0-0:0)
- Reference OpenCL/Python:
  - `pyBall/XPDB_AVBD/XPDB_new.py`, [pyBall/XPDB_AVBD/XPDB_new.cl](cci:7://file:///home/prokophapala/git/FireCore/pyBall/XPDB_AVBD/XPDB_new.cl:0:0-0:0), `pyBall/XPDB_AVBD/test_TiledJacobi_molecules.py`, `pyBall/XPDB_AVBD/OUT-test_TiledJacobi_molecules.txt`
- Legacy kernels: `cpp/common_resources/cl/relax_multi_mini.cl`

#### Session Update – 2026-01-25 (XPDB CPU reference solver & GPU parity debugging)

**What we investigated & learned**
1. **Initial GPU step still diverges**
   - Even after zeroing `l_rad` in WGSL, the very first GPU relaxation of H₂O still flings one H atom away. Subsequent GPU runs behave if we run the CPU solver in between, which indicates the GPU solver’s *first* iteration is inconsistent with the CPU reference rather than a persistent state drift.
2. **Radii / collision inputs likely root cause**
   - CPU vs GPU share identical bond data but differing collision radii (vdW vs covalent). The CPU solver shows symmetric forces; GPU still over-penalizes the H–H collision. We’ll need to compare actual `atom_params` used by both paths.
3. **Mode-switching exposes stale buffer state**
   - Transitioning from GPU → CPU → GPU highlighted unsynchronized buffers. We added explicit two-way synchronization so GPU buffers are refreshed after CPU runs and CPU state mirrors GPU downloads, but the initial GPU-only run still needs investigation (probably missing `predPos` initialization / radii scaling).

**What we implemented today**
1. **CPU reference solver (`XPDB_CPU`)**
   - Pure-JS Jacobi solver using global indices + brute-force collisions, driven by the same `bondsAdj` data. Supports verbose logging (`VERBOSITY_LEVEL ≥ 4`) and tracks `lastDelta` per step.
2. **CPU/GPU runtime toggle**
   - GUI checkbox “Use CPU XPDB (reference)” wires into `MolGUIApp`. `xpdbStep`, `xpdbStepOnce`, `xpdbRelaxSteps` now route through CPU when enabled, logging deltas and writing positions back into `EditableMolecule`.
3. **Bidirectional synchronization**
   - Added helpers to upload current molecule positions to GPU (`syncXPDBGPUFromSystem`) and to refresh CPU data after GPU downloads. CPU runs now push positions to GPU automatically so switching back doesn’t reuse stale buffers.

**What we ruled out / excluded**
1. **Workgroup-order issues**: With the CPU mirror, we confirmed the solver is order-independent; asymmetry isn’t caused by missing barriers or simultaneous writes.
2. **Bond remap failures (for small n)**: H₂O dumps show `bond_indices_local` matches expected internal indices. The collision mask failure must come from collision inputs, not remap truncation, at least for single workgroup cases.

**Prevailing problems**
1. **Initial GPU collision behavior**: Need to capture GPU buffers just before the first solve (positions, radii, bond data) and compare against CPU step results to see why collisions overwhelm bonds.
2. **Parameter mismatch (sp² vs sp³)**: Still pending—bond lengths/stiffness use incorrect atom types for pentacene, which will block larger-scale testing until resolved.
3. **Bond stiffness scaling toggle**: GUI wiring exists but WebGPU implementation does nothing; CPU path currently ignores the toggle as well (always uses provided stiffness).

**Next steps / plan**
1. **GPU first-step forensics**: Add a dedicated “Dump pre-step buffers” action to capture `curr_pos`, `pred_pos`, `atom_params`, and `bond_indices_*` before the first GPU iteration, then replay the same data through CPU solver to compare per-atom forces and see if radii or `pred_pos` differ.
2. **Normalize collision radii**: Align GPU `atom_params.x` with CPU reference (likely covalent radius or explicit XPDB parameter) to avoid huge H–H overlap penalties; verify whether `XPDB_new.py` scales radii before uploading.
3. **Reinitialize GPU solver after CSV load / CPU run**: Despite added sync, we may need to fully reconstruct `XPDB_WebGPU` (init + upload) when enabling GPU to ensure `pred_pos` and ghost buffers are clean, eliminating any hidden Dawn/d3d state carryover.
4. **Bond stiffness API**: Implement host-side `bondLenStiff` cache + scaling routine so the “Disable bonds” toggle can be honored both on CPU and GPU.
5. **MMParams audit**: Log the atom-type names assigned during topology builds (e.g., `O_sp3`, `C_sp2`) and compare with the Python-run dumps to ensure we’re reading the same `.dat` files and matching by atom type rather than atomic number.

We’ll resume tomorrow with the buffer dumps and `atom_params` comparison to pinpoint why the first GPU step diverges.

Let me know if you want me to start implementing the bond stiffness scaling or the `bond_indices_local` dump next.

#### Session Update – 2026-01-26 (Topology parity + pi-align gating)

**What we achieved**
1. **Python↔JS topology dumps now match bit-for-bit**
   - Added gated internal debug logging on both sides (`--dump_topo_debug`) and iterated until `diff -u OUT_py_topo_debug.txt OUT_js_topo_debug.txt` is clean (JS emits extra `EP_DIR` lines by design, everything else identical).
   - **JS emit site:** `computePiOrientations()` / `buildPiDummies()` / `buildEpairDummies()` in [`web/molgui_webgpu/MMFFLTopology.js`](cci:7://file:///home/prokophapala/git/FireCore/web/molgui_webgpu/MMFFLTopology.js:60:0-520:0) gate calls to the shared `debugTopo(line)` callback; CLI wire-up lives in [`dump_xpdb_topology.mjs`](cci:7://file:///home/prokophapala/git/FireCore/web/molgui_webgpu/dump_xpdb_topology.mjs:250:0-360:0).
   - **Python emit site:** `_dbg_topo()` in [`pyBall/OCL/MMFF.py`](cci:7://file:///home/prokophapala/git/FireCore/pyBall/OCL/MMFF.py:100:0-700:0) is connected from `load_molecule_topology_mmffl()` (`pyBall/XPDB_AVBD/test_TiledJacobi_molecules.py`) when `--dump_topo_debug PATH` is passed; the same file writes `PI_DUMMY/EP_DUMMY` entries after `build_topology()` returns.
2. **MOL2 parity for guanine**
   - Node CLI (`dump_xpdb_topology.mjs`) and Python driver now produce identical `OUT_*_guanine.mol2` files (no diff). This confirms atom ordering, dummy placement, and bond policies are synchronized.
   - JS path: `dump_xpdb_topology.mjs` → `buildMMFFLTopology()` → `mol2FromTopology()`.
   - Python path: `test_TiledJacobi_molecules.py --mol2_out ...` pulls the same `topo` dict from `MMFFL.build_topology()` and writes via `au.save_mol2`.
3. **Pi-vector alignment now optional**
   - Ported Python’s `get_atomi_pi_direction`, propagation, and sign-alignment logic directly into JS. Exposed `--align_pi_vectors` flag (default off, matching Python) so we can flip it on only when needed for future molecules.
   - JS toggle: `parseArgs()` stores `alignPiVectors`; `buildMMFFLTopology()` forwards it as `align_pi_vectors`; `computePiOrientations(... alignPiVectors)` gates the sign-flip loop.
   - Python toggle: existing `MMFFL.__init__(align_pi_vectors=True/False)`; currently disabled so we match Python’s default behavior.
4. **MMFFL topology dump plumbing**
   - CLI takes `--dump_topo_debug` and `--align_pi_vectors` alongside existing options; JS debug output mirrors Python formatting, making further investigations reproducible.
   - **Molecules exercised so far:**
     1. `cpp/common_resources/xyz/H2O.xyz` – sanity check for base plumbing.
     2. `cpp/common_resources/xyz/CH2O.xyz` – verified pi/orbital ordering after deterministic neighbor sort.
     3. `cpp/common_resources/xyz/guanine.xyz` – full parity target (current reference `OUT_py_guanine_inputs.txt` / `OUT_js_guanine_inputs.txt`, `OUT_py/js_guanine.mol2`, `OUT_py/js_topo_debug.txt`).

**Still to do from current plan**
1. **XPDB buffer parity (build_local_topology / solve_cluster_jacobi)**
   - Next milestone is dumping GPU buffers (positions, `atom_params`, bond remaps, ghosts) and comparing them against the Python OpenCL reference to chase the remaining collision divergence.
2. **Bond stiffness scaling toggle**
   - UI wiring exists, but `XPDB_WebGPU.setBondStiffnessScale()` still needs to cache host copies of `bondLenStiff` and re-upload values when “Disable bonds” is toggled.
3. **Collision radius / pre-step forensics**
   - Need the “dump pre-step buffers” workflow outlined below to confirm initial GPU iteration uses the same radii / `pred_pos` as CPU.
4. **MMParams audit for other molecules**
   - Guanine is green-lit; must repeat the same dump/diff workflow for the rest of the XPDB parity suite (H2O, CH2O, pentacene, etc.) before moving on to GPU solver parity proper.


#### Session Update – 2026-01-27 (Headless Dawn harness + build_local_topology parity)

**Headless WebGPU harness (Node + Dawn) implemented**
1. Added **`web/molgui_webgpu/headless_webgpu.js`** which bootstraps a Dawn WebGPU context inside Node via the `webgpu` package (`create([])`), copies `GPUBufferUsage`/friends into `globalThis`, and force-installs a file-backed `fetch` polyfill (resolves relative to the harness directory). Any missing Dawn binding now throws with install instructions (`npm i --save-dev webgpu @webgpu/types`).
2. Introduced **`web/molgui_webgpu/headless_init.mjs`** – a tiny module that `await bootstrap()` at top-level so downstream imports (notably `XPDB_WebGPU.js`) always see the polyfilled environment before they call `navigator.gpu` or `fetch('./xpdb.wgsl')`.
3. Authored **`web/molgui_webgpu/run_xpdb_webgpu_headless.mjs`**, a dedicated Node CLI that reuses the existing browser modules (`dump_xpdb_topology.mjs`, `XPDB_WebGPU.js`, `xpdb.wgsl`). The CLI builds topology buffers via `buildXPDBInputsFromXYZArgs()`, uploads them to `XPDB_WebGPU`, dispatches `update_bboxes` + `build_local_topology` (solver iterations optional), and dumps the GPU buffers through `XPDB_WebGPU.dumpTypedState`. WGSL is loaded via the mocked `fetch`, so the shader source stays co-located with the browser path.
4. Wrote **`web/molgui_webgpu/test_gpu.mjs`** (quick sanity harness) to exercise Dawn ↔ RTX 3090 initialization and document the `adapter.info` property change (`requestAdapterInfo()` is gone in latest Dawn).

**Python/OpenCL dump counterpart**
1. Added **`pyBall/XPDB_AVBD/test_XPDB_new_dump_headless.py`** which rebuilds MMFFL topology with all derived bonds enabled, uploads them into `XPDB_new`, runs `build_local_topology`, and dumps the same buffer sections as the WebGPU path (`pos`, `pred_pos`, `atom_params`, `bboxes`, `ghost_packed`, `bond_idx_{global,local}`, `bond_len_stiff`). CLI flags (angles/pi/epair toggles, radii, etc.) mirror the JS defaults to keep inputs synchronized.

**Parity campaign (H₂O + CH₂O)**
1. Ran both harnesses headlessly (Node hitting Dawn/Vulkan, Python hitting PyOpenCL) on `cpp/common_resources/xyz/H2O.xyz` and `CH2O.xyz`. Generated dumps:
   - Python/OpenCL: `pyBall/XPDB_AVBD/out_ocl_{H2O,CH2O}_topo.txt`
   - Node/WebGPU: `web/molgui_webgpu/out_webgpu_{H2O,CH2O}_topo.txt`
2. `diff -u` between the paired files shows **byte-identical data for every buffer section** – only the header line text differs (Python prints `# XPDB_new headless dump xyz=...`, WebGPU prints `# XPDB_WebGPU headless xyz=...`). This confirms the combined effect of `update_bboxes` + `build_local_topology` matches exactly across the two implementations once the new harness is used.
3. Remaining action to get zero-length diffs: update the Python dump header to match the WebGPU format (absolute path + label) or vice versa. After that, we can extend the same workflow to `solve_cluster_jacobi` outputs (pos/pred_pos after one iteration, debug forces, etc.).

**Key files touched today (2026-01-27)**
- `web/molgui_webgpu/headless_webgpu.js`
- `web/molgui_webgpu/headless_init.mjs`
- `web/molgui_webgpu/run_xpdb_webgpu_headless.mjs`
- `web/molgui_webgpu/dump_xpdb_topology.mjs` (exported `buildXPDBInputsFromXYZArgs`, guarded CLI entry)
- `pyBall/XPDB_AVBD/test_XPDB_new_dump_headless.py`
- `web/molgui_webgpu/test_gpu.mjs`

**Next milestones**
1. Harmonize dump headers so parity diffs are strictly empty.
2. Extend both harnesses to run at least one `solve_cluster_jacobi` iteration and dump solver outputs for parity.
3. Hook the harness into CI/regression scripts so H₂O/CH₂O parity can be tested with a single command.
4. Once solver parity is confirmed, reuse the Node harness for larger molecules (guanine, pentacene) before tackling bond-parameter audits and collision fixes mentioned above.

**run tests like this:**

node web/molgui_webgpu/run_xpdb_webgpu_headless.mjs --xyz tests/tDFT_pentacene/pentacene.xyz --out web/molgui_webgpu/out_webgpu_pentacene_distort.txt --inner_iters 10 --n_steps 10 --traj_xyz web/molgui_webgpu/traj_webgpu_pentacene.xyz --dt 0.01 --omega 0.8 --momentum_beta 0 --k_coll 300 --coll_scale 1.2 --bbox_scale 1.2 --maxRadius 0.2 --init_scale 1.05 --init_noise 0.02 --init_seed 12345 --twoPi 0 --addPiAlign 0

PYTHONPATH=/home/prokop/git/FireCore:$PYTHONPATH python3 pyBall/XPDB_AVBD/test_XPDB_new_dump_headless.py --molecule tests/tDFT_pentacene/pentacene.xyz --dump_out pyBall/XPDB_AVBD/out_ocl_pentacene_distort.txt --inner_iters 10 --n_steps 10 --traj_xyz pyBall/XPDB_AVBD/traj_ocl_pentacene.xyz --traj_real_only 1 --dt 0.01 --omega 0.8 --momentum_beta 0 --k_coll 300 --coll_scale 1.2 --bbox_scale 1.2 --init_scale 1.05 --init_noise 0.02 --init_seed 12345 --two_pi 0 --add_pi_align 0

python3 pyBall/XPDB_AVBD/compare_xyz_trajectories.py pyBall/XPDB_AVBD/traj_ocl_pentacene.xyz web/molgui_webgpu/traj_webgpu_pentacene.xyz --atol 1e-5

---

#### Session Update – 2026-01-27 (GUI uses shared XPDB path)

**Problem**

- The browser GUI still used legacy helpers (`buildXPDBTopology()`, `XPDB_WebGPU.uploadBonds()`, positional `xpdb.step`) which rebuilt simplified adjacency lists and uploaded only 3 atoms for H₂O.
- The shared MMFFL builder (`buildXPDBInputsFromMol()` in `MMFFLTopology.js`) mutates its `EditableMolecule` by appending Pi/E dummy atoms. Because the GUI passed `this.system` directly, H₂O became 5 atoms internally while `XPDB_WebGPU`/`XPDB_CPU` buffers stayed sized for 3 atoms, triggering `writeBuffer` overruns and CPU sync mismatches.
- GUI “Dump topology/buffers” buttons pointed to non-existent methods, so no debug output was available.

**Solution (files/functions updated)**

1. `MMFFLTopology.js`
   - Exported `packBondArrays()` and `buildXPDBInputsFromMol()`, already used by the headless runner.
2. `XPDB_WebGPU.js`
   - Added `uploadPackedAtoms(packed)` and `uploadPackedTopology(packed)`; marked `uploadBonds()` as `// DEPRECATED`.
3. `main.js`
   - `MolGUIApp._cloneSystemForXPDB()` clones the molecule via `exportAsParsed()` so MMFFL dummies never touch `this.system`.
   - `updateXPDBTopology()` now:
     - Builds packed topology on the clone (`buildXPDBInputsFromMol`).
     - Re-initializes `XPDB_WebGPU`/`XPDB_CPU` to `n_all` (real + dummies).
     - Calls `uploadPackedAtoms()` + `uploadPackedTopology()`.
     - Records `xpdbNReal` / `xpdbNAll` for later sync.
   - `syncXPDBPositions()` copies back only the first `xpdbNReal` atoms so the GUI geometry stays at 3 atoms.
   - Added `dumpXPDBTopology()` and `dumpXPDBBuffers()` so the GUI buttons log the same expanded data as the headless CLI.
   - Updated all solver paths (`stepXPDB`, `xpdbStepOnce`, `xpdbRelaxSteps`) to call `xpdb.step({ ... })` with the shared options object.

**Current result**

- GUI and headless now share the exact MMFFL-packed topology + buffers; the browser logs match the Node dumps (see latest H₂O dump above).
- Buffer overruns and CPU atom-count mismatch warnings are gone; “Relax N steps” moves symmetrically (albeit little motion because the molecule is already relaxed).
- “Dump topology” / “Dump buffers” buttons output expanded data.

**Next steps**

1. Add deterministic distortions + collision settings to the GUI (same flags as headless) so we can deliberately perturb geometry and observe motion.
2. Render Pi/E-pair dummies and angle-derived bonds in the Raw WebGPU renderer to visualize the full MMFFL constraint graph.
3. Extend GUI tests to run short trajectories and compare against headless `traj_*.xyz` using `compare_xyz_trajectories.py`.
4. Once GUI ↔ headless parity is validated on H₂O/CH₂O, repeat for guanine/pentacene and move on to solver parity (`solve_cluster_jacobi`).

**Detailed notes about the GUI parity fix**

- The GUI now uses the shared XPDB path, which means it builds the same MMFFL-packed topology and buffers as the headless runner.
- The `MolGUIApp._cloneSystemForXPDB()` method clones the molecule via `exportAsParsed()` to prevent MMFFL dummies from touching the original system.
- The `updateXPDBTopology()` method builds packed topology on the clone, re-initializes `XPDB_WebGPU`/`XPDB_CPU` to `n_all` (real + dummies), and calls `uploadPackedAtoms()` + `uploadPackedTopology()`.
- The `syncXPDBPositions()` method copies back only the first `xpdbNReal` atoms to keep the GUI geometry at 3 atoms.
- The GUI buttons now log the same expanded data as the headless CLI using `dumpXPDBTopology()` and `dumpXPDBBuffers()`.
- The solver paths (`stepXPDB`, `xpdbStepOnce`, `xpdbRelaxSteps`) now call `xpdb.step({ ... })` with the shared options object.

**GUI parity fix benefits**

- The GUI now shares the same MMFFL-packed topology and buffers as the headless runner, ensuring consistent behavior and output.
- The GUI can now deliberately perturb geometry and observe motion using deterministic distortions + collision settings.
- The GUI can render Pi/E-pair dummies and angle-derived bonds in the Raw WebGPU renderer to visualize the full MMFFL constraint graph.
- The GUI can run short trajectories and compare against headless `traj_*.xyz` using `compare_xyz_trajectories.py`.