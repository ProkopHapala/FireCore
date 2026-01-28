## USER

is it easy to initialize WebGL with Three.js and conditionally inititalize WebGPU if available? 

I have program in WebGL for molecular visulization and it is quite cosrvative for back compactibility (e.g. on linux suppoer of WebGPU is still experimental and need user to allor some special option vulkan and other, I need run browser like this
`google-chrome --enable-unsafe-webgpu --ozone-platform=x11 --use-angle=vulkan --enable-features=Vulkan,VulkanFromANGLE --ignore-gpu-blocklist`

therefore I do not want to fully rely on WebGPU, but for some things like running interactive molecular dynamics, I really want to use compute shaders as they map perfectly from my existing OpenCL kernels. So I was thinking to make just the simulation module need webgpu and always report if it is avaialbe or not. But I wo not want to double my rendering systems or generally I wo not want to touch my WebGL rendering system. Question is if I have WebGPU and WebGL system initialized at the same time, if I can share buffwers, for example if particle position buffer which is writen by WebGPU simulation, can be rendered using WebGL system and GLSL shaders.

---

## Gemini

Initializing WebGL with a conditional WebGPU fallback is technically easy in Three.js, but your specific requirement—**sharing buffers between a WebGPU simulation and a WebGL renderer**—is currently impossible without a significant performance penalty.

### 1. Can you share buffers between WebGL and WebGPU?
**No, not directly.** WebGL and WebGPU are separate APIs with isolated memory contexts. Even if they are running on the same physical GPU, a buffer created in a WebGPU context (`GPUBuffer`) cannot be bound or accessed by a WebGL context (`WebGLBuffer`) without moving the data through the CPU.

*   **The Cost:** To share data, you would have to:
    1. Run your WebGPU simulation.
    2. Use `gpuBuffer.mapAsync()` to read the data back to a JavaScript `TypedArray`.
    3. Use `gl.bufferSubData()` to upload that same array to your WebGL vertex buffer.
*   **Result:** This "GPU $\to$ CPU $\to$ GPU" round-trip is a major bottleneck. For a high-performance molecular dynamics simulation with thousands of atoms, this will likely negate the speed gains you get from using compute shaders.

### 2. The Three.js Solution: `WebGPURenderer`
Three.js has introduced a new renderer, `THREE.WebGPURenderer`, specifically designed to solve the "WebGPU with WebGL fallback" problem. 

*   **How it works:** You initialize one renderer (`new THREE.WebGPURenderer()`). It automatically detects if WebGPU is available. If not, it falls back to a WebGL 2 backend.
*   **The "One Buffer" Magic:** When using this renderer, you use `StorageBufferAttribute`. Three.js manages the underlying memory. If WebGPU is active, it uses a real GPU storage buffer. If it falls back to WebGL, it emulates the behavior (often using **Transform Feedback** or standard attribute updates).

### 3. The Catch: Your Existing GLSL Shaders
The biggest hurdle is that `WebGPURenderer` **does not support `ShaderMaterial` or `RawShaderMaterial`**.
*   **TSL (Three Shading Language):** To work with the new renderer, you must use the **NodeMaterial** system and **TSL**. TSL is a JavaScript-based language that Three.js transpiles into **WGSL** (for WebGPU) or **GLSL** (for the WebGL fallback).
*   **Porting effort:** Since you mentioned your molecular visualization is "quite conservative" and uses existing GLSL/OpenCL, porting to TSL would be a requirement to keep a single, shared system.

### Summary Comparison

| Approach | Buffer Sharing | Ease of Setup | Compatibility |
| :--- | :--- | :--- | :--- |
| **Two Separate Systems** | Needs CPU Copy (Slow) | Moderate | Keeps existing GLSL system. |
| **Three.js `WebGPURenderer`** | Native GPU Sharing (Fast) | High | Requires porting GLSL to TSL. |

### Recommendation
If you want to use WebGPU for simulation while keeping your WebGL rendering system "untouched," you will be forced to use the slow **CPU copy method**. 

However, if performance is the priority:
1.  **Switch to `WebGPURenderer`:** This handles the conditional initialization and the Linux experimental flags gracefully.
2.  **Report WebGPU Status:** You can check `renderer.backend.isWebGPU` to decide whether to enable your "Interactive MD" module.
3.  **Port Simulation only:** You can write your MD simulation in TSL (which maps very well to OpenCL logic). If the browser falls back to WebGL, Three.js will attempt to run it via Transform Feedback, though it will be slower than the native WebGPU compute path.

---

## USER

OK, here is rendering system from WebGL, can you rewrite it to the new system compactible with WebGPU which should work even if WebGPU is not available?

- [main.js (whole file)](#web-molgui-web-js-main-js)
- [MoleculeRenderer.js (whole file)](#web-molgui-web-js-moleculerenderer-js)
- [MeshRenderer.js (whole file)](#web-common-js-meshrenderer-js)
- [Draw3D.js (whole file)](#web-common-js-draw3d-js)
- [atom.glslf (whole file)](#cpp-common-resources-shaders-atom-glslf)
- [atom.glslv (whole file)](#cpp-common-resources-shaders-atom-glslv)
- [bond_color.glslf (whole file)](#cpp-common-resources-shaders-bond-color-glslf)
- [bond.glslv (whole file)](#cpp-common-resources-shaders-bond-glslv)
- [label.glslf (whole file)](#cpp-common-resources-shaders-label-glslf)
- [label.glslv (whole file)](#cpp-common-resources-shaders-label-glslv)

---

## Gemini

To modernize your system for **Three.js + WebGPU** (with automatic WebGL fallback) while maintaining a **single codebase**, you need to migrate from raw GLSL strings (`ShaderMaterial`) to the **Three.js Shading Language (TSL)** using the **Node Material** system.

This approach solves your problem:
1.  **Single System**: You write materials in TSL (JavaScript). Three.js automatically compiles this to **WGSL** (if WebGPU is active) or **GLSL** (if falling back to WebGL).
2.  **Future Proof**: TSL is required to share data efficiently with WebGPU Compute Shaders (which you want for your MD simulation).

Here is the rewritten rendering system.

### Prerequisites
You must be using a recent version of Three.js (r167+) and have the `three/addons/` mapped correctly in your import map.

### 1. main.js
**Changes**:
- Removed manual `.glsl` file fetching (shaders are now JS/TSL).
- Replaced `WebGLRenderer` with `WebGPURenderer`.

### 2. MeshRenderer.js
**Changes**:
- Replaced `ShaderMaterial` with TSL Nodes.
- Implemented the custom logic (texture fetch, impostors) using `three/tsl`.
- This file contains the logic previously found in `atom.glsl`, `bond.glsl`, etc.

### 3. Draw3D.js
**Changes**:
- Removed the creation of `ShaderMaterial` inside the helpers, or modified them to be compatible.
- `createLabelInstancedMesh` needs to return the geometry/mesh structure but let `MeshRenderer` attach the `NodeMaterial`.

*(Modify `createLabelInstancedMesh` to accept null shaders)*

### 4. MoleculeRenderer.js
**Changes**:
- Minor: Pass `null` for shaders in constructor call (already updated in `main.js` and handled in `MeshRenderer.js`).

### Summary of Benefits for your Goal

1.  **WebGL + WebGPU Compatibility**: The `WebGPURenderer` will attempt to use WebGPU. If your Linux setup flags allow it (`--enable-unsafe-webgpu`), it works. If not, it falls back to WebGL2 backend automatically. The TSL code works on both.
2.  **Buffers**: By using `DataTexture` for positions (as you did in WebGL), we keep the CPU-side logic identical.
3.  **Compute Readiness**: When you are ready for MD, you can create a `StorageInstancedBufferAttribute` for positions. You can then write a TSL Compute Node (`Fn().compute()`) to update those positions on the GPU. The `MeshRenderer` can then simply read that storage buffer instead of the `DataTexture`, keeping the rendering logic almost identical.

This refactor removes the need for conditional initialization logic in your application code—Three.js handles it internally.


---

## USER

now we should accuretely copy this geometry assignment to our javacript implementation of the same problem. which is here
/home/prokop/git/FireCore/web/molgui_webgpu

see the relevant classes
@XPDBTopology.js should be equvalent to @MMFFL.py 
@XPDB_WebGPU.js equvalent to @XPDB_new.py 
@xpdb.wgsl equvalent to @XPDB_new.cl 
@EditableMolecule.js equvalent to @AtomicSystem.py 
and @MMParams.js equivalent to @MMparams.py 

find where are dicrepancies more specifically what is missing in the javascript version.

our goal is to make partity and check it by exporting .mol2 form both, writing the same topology information to console (in same format so it can be excatly compared), and also adding same visualization of the molecule with added dummy atoms, angle/pi/epair-derived bonds like plot_topology_3d from 
@test_TiledJacobi_molecules.py 

perhps for testing it would be good if you implement standalone script which you run in node.js instead of browser which prints thosse outputs and emit .mol2 file, so you can test it headless partity exactly with respect to the python referecne without browser GUI. Can you do it. Then you should run the test end check the parity of python/javascript outputs
make detailed plan how to do that

see the current format of ouput from @test_TiledJacobi_molecules.py which we want to reimplement
``
TYPE ia=   0 ename=O type=O        nbond=1 npi=1 nepair=2
TYPE ia=   1 ename=C type=C        nbond=3 npi=1 nepair=0
TYPE ia=   2 ename=H type=H        nbond=1 npi=0 nepair=0
TYPE ia=   3 ename=H type=H        nbond=1 npi=0 nepair=0
read_AtomAndElementTypes() path /home/prokop/git/FireCore/pyBall/OCL/../../cpp/common_resources/
path1 /home/prokop/git/FireCore/pyBall/OCL/../../cpp/common_resources/ElementTypes.dat
path2 /home/prokop/git/FireCore/pyBall/OCL/../../cpp/common_resources/AtomTypes.dat
read_element_types() filepath= /home/prokop/git/FireCore/pyBall/OCL/../../cpp/common_resources/ElementTypes.dat
DEBUG initAtomProperties: atom 0 O element=O marked node
DEBUG initAtomProperties: atom 1 C element=C marked node
DEBUG MMFF: ia=0 O nbond=1 conf_npi=1 conf_ne=2
DEBUG MMFF: init pi_vec ia=0 -> [0. 0. 0.] (raw=[0. 0. 0.])
DEBUG MMFF: ia=1 C nbond=3 conf_npi=1 conf_ne=0
DEBUG MMFF: init pi_vec ia=1 -> [1. 0. 0.] (raw=[0. 0. 0.])
DEBUG MMFF: propagate pi_vec ia=0 -> [1. 0. 0.] using neighbors [1]
DEBUG MMFF: final pi_vec ia=0 -> [1. 0. 0.] (norm=1.0)
DEBUG MMFF: final pi_vec ia=1 -> [1. 0. 0.] (norm=1.0)
[TOPO] n_real= 4  n_all= 10  n_primary= 3  n_derived= 19
[TOPO] type_names:
     0 O 
     1 C 
     2 H 
     3 H 
     4 Pi DUMMY
     5 Pi DUMMY
     6 Pi DUMMY
     7 Pi DUMMY
     8 E DUMMY
     9 E DUMMY
[TOPO] bonds_primary: [(0, 1), (1, 2), (1, 3)]
[TOPO] bonds_angle: [(0, 2), (0, 3), (2, 3)]
[TOPO] bonds_pi: [(0, 4), (0, 5), (0, 6), (0, 7), (1, 4), (1, 5), (1, 6), (1, 7), (2, 6), (2, 7), (3, 6), (3, 7)]
[TOPO] bonds_epair: [(0, 8), (0, 9), (1, 8), (1, 9)]
[TOPO] wrote /tmp/OCH2_topo.mmffl.mol2
``

now make detailed analysis of the referenced files and detailed notes what should be implemented and where

---

## GPT-5.1-Codex

Parity gaps between the Python reference (MMFFL + XPDB driver) and the current WebGPU stack:

1. **Atom typing / reporting**
   - Python: [MMFFL.assign_type_names(...)](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/MMFFL.py:73:4-95:18) + [report_types_and_valence(...)](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/MMFFL.py:105:4-123:125) print `TYPE ia=... nbond npi nepair` and enforce octet rule using MMparams data.
   - JS: [MMParams.getAtomTypeForAtom()](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/MMParams.js:123:4-136:5) picks types lazily, but [EditableMolecule](cci:2://file:///home/prokop/git/FireCore/web/molgui_webgpu/EditableMolecule.js:240:0-1055:1) never emits the per-atom report nor asserts valence consistency. No CLI/UI hook to select between “table” vs “UFF” typing.
   - Missing: equivalent reporting/logging in WebGUI, plus typed outputs for comparison.

2. **Topology builder scope**
   - Python [MMFFL.build_topology()](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/MMFFL.py:125:4-218:9) delivers primary bonds + derived angle, pi-host, pi-align, epair-host, epair-orth, epair-pair lists, dummy bookkeeping, `is_dummy` mask, aggregated `apos`, MOL2-ready names, etc.
   - JS [buildXPDBTopology()](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/XPDBTopology.js:2:0-201:1) only:
     - reads real bonds from EditableMolecule,
     - derives angle bonds via [MMParams.getAngleParams](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/MMParams.js:411:4-426:5) and [convertAngleToDistance](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/MMParams.js:428:4-444:5),
     - merges them immediately into a single adjacency per atom.
   - Missing features:
     1. **Pi dummies**: no generation of pi positions, no host↔dummy bonds, no dummy↔neighbor springs, no pi-align bonds, no dummy count mask.
     2. **Electron-pair dummies**: not created automatically; [EditableMolecule.addExplicitEPairs()](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/EditableMolecule.js:417:4-484:5) exists but is manual (uses VSEPR heuristics, not MMFFL directions), and even then XPDBTopology ignores “atoms” with `Z=200` when building adjacency (it doesn’t treat them specially).
     3. **Separation of first- vs second-order bonds**: JS merges angle bonds into the same list immediately, so downstream code can’t distinguish classes (needed for selective uploads or plotting).
     4. **Dummy bookkeeping**: no `n_real`, `n_dummy`, `is_dummy[]`, or dummy masks for collision filters.
     5. **Pi / epair stiffness parameters**: no `L_pi`, `L_epair`, or `Kpi_* / Kep_*` equivalents; MMParams has Kss/Ksp but XPDBTopology never uses them.

3. **MMParams data parity**
   - Python `MMparams` loads ElementTypes/AtomTypes/BondTypes/AngleTypes from files on disk, uses them everywhere.
   - JS [MMParams](cci:2://file:///home/prokop/git/FireCore/web/molgui_webgpu/MMParams.js:62:0-470:1) parser mirrors the files but is missing:
     - [convertAngleToDistance](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/MMParams.js:428:4-444:5) behavior identical to [_law_of_cosines](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/MMFFL.py:7:0-11:42)? Need to confirm; if not, parity may break lengths.
     - [getAngleParams](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/MMParams.js:411:4-426:5) in JS only matches by type names; the Python version uses center type’s `Ass` and `Kss`. If WebGUI uses raw angle table while Python uses per-atom-type data, results diverge.
     - Electron pair counts (`nepair`) exist in JS atom types, but the topology builder never uses them.

4. **EditableMolecule vs AtomicSystem**
   - Python [AtomicSystem](cci:2://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:16:0-975:78) stores `mol.enames`, `mol.bonds` (list of pairs), neighbor lists, `npi_list`, `nep_list`, and [add_electron_pairs()](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:874:4-887:52) (mutating). The XPDB script builds bonds, prints topology, saves MOL2.
   - JS [EditableMolecule](cci:2://file:///home/prokop/git/FireCore/web/molgui_webgpu/EditableMolecule.js:240:0-1055:1) tracks atoms/bonds, provides VSEPR helpers, but lacks:
     - `findAngles/neighbor lists` akin to [mol.findAngles](cci:1://file:///home/prokop/git/FireCore/pyBall/atomicUtils.py:335:0-349:21).
     - Mechanism to compute `npi/nepair` per atom automatically (only manual selection-based functions).
     - [.save_mol2()](cci:1://file:///home/prokop/git/FireCore/pyBall/atomicUtils.py:1710:0-1769:95) equivalent. No path to export combined real+dummy positions/bonds; WebGUI can dump JSON but not MOL2 for offline comparison.
     - Matplotlib-like visualization: Web GUI can render in 3D, but there’s no quick debug view replicating [plot_topology_3d](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/test_TiledJacobi_molecules.py:276:0-324:14) (color-coded bonds, dummy scatter).

5. **XPDB driver parity**
   - Python [test_TiledJacobi_molecules.py](cci:7://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/test_TiledJacobi_molecules.py:0:0-0:0) prints `[TOPO]` stats, writes `.mol2`, optional 3D plot, and pipelines topology-only mode for debugging.
   - Web stack lacks any CLI/GUI command to:
     - dump `[TOPO]` blocks in the same format,
     - export `.mol2`,
     - show dummy-specific visualization (colors, separate legends),
     - run topology-only validation vs CPU reference.
   - `dump_xpdb_topology.mjs` prints adjacency, but format differs and doesn’t include dummy info because dummies aren’t generated.

6. **Target parity requirements (based on request)**
   - Need to “accurately copy” geometry assignment:
     - Implement pi and epair dummy creation in JS (likely inside `XPDBTopology` or a new `LinearizedTopology` module mirroring MMFFL). This includes orientation vectors, lengths (`L_pi`, `L_epair`), stiffness constants, optional pi-align / epair-pair bonds.
     - Maintain separate bond lists by class, plus aggregated adjacency for GPU.
     - Track `n_real` vs `n_dummy`, `type_names`, `is_dummy` mask.
   - Matching output:
     - Console output identical to Python: `TYPE ...` lines, `[TOPO]` summary, bond lists.
     - `.mol2` export from Web path (maybe Node script) to compare with Python output.
     - Visualization: add a debug view (either via WebGUI overlay or Node+three.js) matching [plot_topology_3d](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/test_TiledJacobi_molecules.py:276:0-324:14)—distinct colors for primary/angle/pi/epair bonds and dummy types.

7. **Suggested implementation plan**
   1. **Shared topology module**: Port MMFFL logic into JS (`LinearizedMMFF`), using [MMParams](cci:2://file:///home/prokop/git/FireCore/web/molgui_webgpu/MMParams.js:62:0-470:1) for data. Steps:
      - Assign type names per atom (`type_source` parameter).
      - Compute neighbor counts, `npi`, `nepair` (via [initAtomProperties](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/MMFF.py:22:0-59:37) analog).
      - Build angle bonds, pi dummies (with host/neighbor bonds, optional two-sided, align bonds), epair dummies (host/neighbor bonds, optional pair bonds).
      - Produce topology object with same schema as Python.
   2. **EditableMolecule integration**:
      - Provide `calcNeighbors()` to mimic [mol.neighs()](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:150:4-154:23) (adjacency arrays).
      - Provide [findAngles()](cci:1://file:///home/prokop/git/FireCore/pyBall/atomicUtils.py:335:0-349:21) to yield (ja, ib, jc) triplets.
      - Add methods to append dummies (positions + metadata) without interfering with existing atoms (maybe store separate arrays).
   3. **Logging / assertions**:
      - Mirror [report_types_and_valence](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/MMFFL.py:105:4-123:125) output (including octet rule) in JS console (maybe gate by verbosity).
   4. **Dump/export tooling**:
      - Implement `.mol2` writer in JS (shared util) to dump combined atoms/bonds.
      - Add a debug command/button to run topology-only mode: compute topology, log `[TOPO]` blocks, produce `.mol2`, render 3D overlay with colored bonds/dummies (matching Python colors).
   5. **Parity verification**:
      - Run both Python and JS on same molecule (H2O, OCH2, pentacene). Compare console logs and `.mol2` files.

These steps cover the missing functionality and set the stage for GPU parity (later hooking into XPDB_WebGPU.js and shaders).



Here’s a concrete plan for the headless Node.js parity harness:

### 1. Assemble a Node-side “linearized topology” module
1.1 **Port MMFFL logic to JS or wrap existing pieces**
- Create `web/molgui_webgpu/linearized_topology.mjs` (ESM for Node) that instantiates [MMParams](cci:2://file:///home/prokop/git/FireCore/web/molgui_webgpu/MMParams.js:62:0-470:1), loads AtomTypes/ElementTypes/BondTypes/AngleTypes from disk, and mimics [MMFFL.build_topology](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/MMFFL.py:125:4-218:9).
- Reuse/extend [XPDBTopology.js](cci:7://file:///home/prokop/git/FireCore/web/molgui_webgpu/XPDBTopology.js:0:0-0:0) helpers but keep angle, pi, epair builders separate so we can produce the same `topo` dict as Python (`apos`, `type_names`, `is_dummy`, bond lists by class).
- Add per-atom type assignment using the exact `AtomTypes.dat` lookup (table mode) plus optional “UFF” path later.

1.2 **EditableMolecule in Node**
- Add an importable `EditableMoleculeNode` wrapper that can:
  - Load `.xyz` files (reuse logic from [atomicUtils.load_xyz](cci:1://file:///home/prokop/git/FireCore/pyBall/atomicUtils.py:1189:0-1231:36) or implement minimal parser).
  - Build neighbor lists / angle triplets.
  - Append dummy atoms without interfering with real indices (track `n_real`, `n_dummy`).

### 2. Implement the CLI script
2.1 **Script skeleton (`tools/xpdb_topology_cli.mjs`)**
```
node xpdb_topology_cli.mjs --xyz path --type_source table --add_pi ... --mol2 out.mol2
```
- Parse args (yargs or manual).
- Load molecule via new EditableMoleculeNode.
- Call linearized topology builder with flags mirroring Python CLI.
- Capture result `topo`.

2.2 **Output parity**
- Print `TYPE ia=...` lines exactly like Python (pad widths, include `nbond`, `npi`, `nepair`).
- Print `[TOPO] ...` block plus lists (`bonds_primary`, `bonds_angle`, etc.) sorted.
- Write `.mol2` using a shared helper (implement `saveMol2` reusing Python format: `@<TRIPOS>` sections, 1-based indices, include dummy atoms).

2.3 **Optional visualization stub**
- For headless parity, skip 3D rendering; but keep hooks to dump JSON (positions + bond categories) so the browser UI can reuse the same data later.

### 3. Verification workflow
3.1 **Reference run (Python)**
```
python pyBall/XPDB_AVBD/test_TiledJacobi_molecules.py --molecule ... --topology_only 1 --mol2_out ref.mol2 --print_buffers 1 > ref.log
```

3.2 **Node run**
```
node web/molgui_webgpu/tools/xpdb_topology_cli.mjs --xyz ... --mol2_out js.mol2 > js.log
```

3.3 **Diffs**
- `diff -u ref.log js.log` (expect identical `[TYPE]/[TOPO]` sections).
- `diff -u ref.mol2 js.mol2` (positions/dummies/bonds should match up to float formatting).

### 4. Integration steps
4.1 **Shared utilities**
- Put `.xyz` loader and `.mol2` writer under `web/molgui_webgpu/util/` so both Node script and browser code can import them (browser bundler must tree-shake Node-only parts).
- Expose the topology builder as both ESM (Node) and standard JS module for browser.

4.2 **Testing in CI**
- Add npm script `npm run topo-test -- --xyz cpp/common_resources/xyz/OCH2.xyz`.
- Optionally integrate with existing test harness to ensure parity doesn’t regress.

This plan yields a headless parity testbed: run Node CLI → compare logs/MOL2 with Python output; iterate until diffs vanish.

---

## USER

Now we should go furter to step 5 checking exactly the solvers, namely @xpdb.wgsl with respect to @XPDB_new.cl 

In particular we want to check exactly line-by-line component by component 

1) output of build_local_topology, especially the local reindexing bond_indices_local, but also the other, I a bit forgot what cotains which arrrey please refresh my memory

2) solve_cluster_jacobi - we should compere the output positions after one iteration. resp. output forces (resp displaceent, because XPDB does ot compute forces but displacements as far as I know, the name out_force_bond and out_force_coll is missleading I think.

for both we should compare input and output buffers with the same name for the kernels with the same name, with the same inputs. Starting from simple molecules like H2O and CH2O, and then moving to something more complicated like Guainne @guanine.xyz 

---

now make detailed plan what we should implement

the output of the simulation can be printed to console, perhaps you can also produce some output log file? Can you?
we cannot run the WebGPU using node.js without browser right?

The think is I need detailed logfiles so I can give it to you in text or in feedback...


anyway, make dtailed notes what to impleemnt and where,and then think about how to produce text outputs which we can compare linie by line. 

OK, I'm not sure we can debug output buffers without browser, but we can ate least check if all input buffers are the same, you sould make function print buffers in the respective .js file like in @MMFFLTopology.js @XPDB_WebGPU.js and then call these function from @dump_xpdb_topology.mjs and the same functions/datastructures must be also in our web application when we run the XPBD simulations (there should be no code duplications, they shoud import same files and functions)

The output uffers I will then check from console termina, just make sure they are writen in same fomat between python and javascript console, and I do not need to unfold any objects in console (must be writen explicitly line by line)

make these detailed implementtion notes expanding my commands 

---

## GPT-5.1-Codex

Plan for solver parity deep-dive:

1. **Baseline mapping of buffers (done here, no code yet)**
   - `build_local_topology`
     - Inputs shared by CL/WGSL: `curr_pos`, `bboxes` (WGSL packs min/max pairs; CL keeps separate min/max arrays), `bond_indices_global`, group count, cluster size.
     - Outputs to compare:
       - `ghost list`: CL writes `ghost_indices_flat` + `ghost_counts`; WGSL stores `[indices..., count]` per cluster inside `ghost_packed`. Need helper to dump both forms into identical JSON/CSV for diffing.
       - `bond_indices_local`: both flatten to `[atom * N_MAX_BONDED + slot]`. Must log per atom the internal index chosen (0–63 internal, 64+ ghost slot). This includes the “local reindexing” you mentioned.
     - Internal scratch (CL’s `__local` arrays vs WGSL `var<workgroup>`). For logging we only care about final global buffers.

   - `solve_cluster_jacobi`
     - Inputs to ensure identical:
       - `curr_pos` (initial), `pred_pos`, `atom_params` (.x radius, .w mass), `bond_indices_local`, `bond_lengths`, `bond_stiffness`, `ghost data`.
       - Solver constants: `inner_iters`, `dt`, `k_coll`, `omega`, `momentum_beta`.
     - Outputs:
       - Updated `curr_pos` (after one Jacobi iteration).
       - CL also writes `out_force_bond`/`out_force_coll`; WGSL currently doesn’t—but for parity we can add optional storage buffer in WGSL and fill with the equivalent `force` accumulations (or rename to displacement if that’s the intended meaning). That will let us diff “forces/displacements” per atom.

2. **Instrumentation/logging strategy**
   - Add optional “debug dump” paths in both pipelines.
     - For CL side (Python driver [test_TiledJacobi_molecules.py](cci:7://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/test_TiledJacobi_molecules.py:0:0-0:0)), we can already print numpy buffers. Extend the script with `--dump_topology foo.json` and `--dump_solver foo.json` to serialize:
       - `ghost_indices_flat`, `ghost_counts`, `bond_indices_local`.
       - After running a single solver iteration, dump `curr_pos`, `out_force_bond`, `out_force_coll`.
       - Use `json.dump` of plain lists so it’s easy to diff line-by-line.
     - For WGSL/WebGPU side (Node CLI [dump_xpdb_topology.mjs](cci:7://file:///home/prokop/git/FireCore/web/molgui_webgpu/dump_xpdb_topology.mjs:0:0-0:0) or a new `dump_xpdb_solvers.mjs`):
       - Hook after compute pass: map the storage buffers back to CPU (`queue.readBuffer`) and write to disk (`fs.writeFileSync`) with the same schema as Python’s JSON.
       - Produce per-molecule logs: e.g., `out_js_H2O_topology.json`, `out_js_H2O_solver.json`.
       - Console print summaries (counts, first few entries), but the JSON files are the authoritative diffable artifacts you can share.
   - Because true WebGPU compute still needs Chrome (or experimental Node/WebGPU flag), we’ll keep using the existing [dump_xpdb_topology.mjs](cci:7://file:///home/prokop/git/FireCore/web/molgui_webgpu/dump_xpdb_topology.mjs:0:0-0:0) runner (it already uses Dawn/WebGPU via Node 22). So yes, we can read buffers and write text files without a browser.

3. **Verification workflow**
   - Molecules: `H2O`, `CH2O`, `guanine` (path `/cpp/common_resources/xyz/<name>.xyz`). Possibly also `pyridine`/`pyrrol` if time.
   - Procedure per molecule:
     1. Run Python XPDB (`test_TiledJacobi_molecules.py --topology_only 1 --dump_topology out_py_<mol>_top.json` etc.) to dump reference buffers.
     2. Run WebGPU CLI to build topology and immediately dump the same buffers (`--dumpTopology out_js_<mol>_top.json`).
     3. Diff the JSON files (simple `diff -u`). Investigate mismatches in `ghost` or `bond_indices_local`.
     4. Once topology matches, run single Jacobi iteration on both sides with identical solver parameters, dumping `curr_pos` and “forces/displacements” to `*_solver.json`. Compare again.
     5. Scale up to guanine after small cases pass.

4. **Implementation notes (“where to add what”)**
   - [dump_xpdb_topology.mjs](cci:7://file:///home/prokop/git/FireCore/web/molgui_webgpu/dump_xpdb_topology.mjs:0:0-0:0): add optional CLI flags `--dumpTopology file`, `--dumpSolver file`, `--innerIters 1`, etc. After compute pass, call helper `writeBufferJSON(buffer, schema, filename)`.
   - `XPDB_WebGPU.js` (or whichever module dispatches kernels): expose functions to map storage buffers and return typed arrays for logging.
   - [test_TiledJacobi_molecules.py](cci:7://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/test_TiledJacobi_molecules.py:0:0-0:0): add `--dump_topology` and `--dump_solver` flags; after calling CL kernels, use `tolist()` to serialize relevant numpy arrays.
   - Both sides should include metadata in the JSON (molecule name, cluster count, constants) so logs are self-describing.

Once we add those dumping hooks, you’ll be able to run the two pipelines for H₂O/CH₂O/guanine, collect the logs, and we can diff them line-by-line to pinpoint any solver divergence. Let me know when you’re ready to implement the instrumentation; we’ll start with the Python script (easier to modify) before wiring up the WebGPU logging.

Implementation outline (zero code yet):

1. Shared buffer-print utilities
   - Add a tiny module (e.g., `web/molgui_webgpu/debugBuffers.js`) exporting pure functions like `formatGhostList(ghostIndices, ghostCounts)` and `formatBondLocal(bondIndicesLocal, groupSize)`. They should output plain strings (one line per entry, e.g., `atom 05 slot 02 -> 071`).
   - [MMFFLTopology.js](cci:7://file:///home/prokop/git/FireCore/web/molgui_webgpu/MMFFLTopology.js:0:0-0:0) and `XPDB_WebGPU.js` both import these helpers so printing logic stays single-sourced for CLI + browser builds.

2. XPDB-WebGPU instrumentation
   - In `XPDB_WebGPU.js`, expose methods `dumpInputsToConsole(label, buffers)` and `dumpOutputsToConsole(...)`. Each method:
     1. Maps the relevant GPUBuffer to CPU (`device.queue.readBuffer` already supported in Node build). Convert to typed arrays.
     2. Pass arrays through the shared formatting helpers.
     3. `console.log` each line explicitly—no console objects, no JSON folding.
   - Cover at least: `curr_pos`, `pred_pos`, `atom_params`, `bond_indices_global`, `bond_indices_local`, `ghost_packed`. Later we can add solver outputs once WebGPU write-back is working.

3. CLI wiring
   - [dump_xpdb_topology.mjs](cci:7://file:///home/prokop/git/FireCore/web/molgui_webgpu/dump_xpdb_topology.mjs:0:0-0:0):
     - Add CLI flags `--dumpInputs` and `--dumpOutputs`.
     - After loading/building topology, call the XPDB dump helpers if flags are set.
     - The CLI should pass file paths too (e.g., `--dumpInputs out_js_H2O_inputs.txt`). The helper writes to both console and optional files via `fs.writeFileSync`.
   - Ensure the same exported functions are used by the browser app so we’re not duplicating debugging code.

4. Python parity
   - [pyBall/XPDB_AVBD/test_TiledJacobi_molecules.py](cci:7://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/test_TiledJacobi_molecules.py:0:0-0:0) already prints buffers via numpy. Add `--dump_inputs file` and `--dump_outputs file` that emit the same line format (matching the JS helpers). This way diffs can be run directly between the two text files.

5. Workflow reminder
   - Start with H₂O and CH₂O: run Python script with `--dump_inputs out_py_H2O_inputs.txt`, then run Node CLI with `--dumpInputs out_js_H2O_inputs.txt`. Diff the text files. Only after inputs match do we move on to solver outputs and larger molecules like guanine.
   - No browser needed: Node 22 + Dawn backend already powers [dump_xpdb_topology.mjs](cci:7://file:///home/prokop/git/FireCore/web/molgui_webgpu/dump_xpdb_topology.mjs:0:0-0:0), so buffer reads/logging happen entirely in the CLI.

Once you’re ready to implement, we’ll add the shared helper module, expose dump functions in `XPDB_WebGPU.js`, extend the CLI flags, and mirror the formatting in the Python driver.

---

## USER

the discrepancies you found are superficial. Instead we must make sure that we call the same functions for topology builder and biffer preparations, and that we use exactly the same parameters (adding electron pairs, dummy pi, and the angle-derived bonds), only that way we can achieve parity. yes, we should apply same distortion for debugging purposes, but the problem is that just the loaded water geometry without any distorion relax a lot, and asymetrically (one hydrogen atoms goes very far), which is not at all what we see in the headless script. we need to undrstand why is this so, but analysing the code. 

---

## GPT-5.1-Codex



Current parity gap is due to deeper differences than distortion:

1. **Different topology builders**  
   - **Python headless** uses [MMFFL.build_topology()](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/MMFFL.py:125:4-218:9) ([pyBall/OCL/MMFFL.py](cci:7://file:///home/prokop/git/FireCore/pyBall/OCL/MMFFL.py:0:0-0:0)) which emits `bonds_primary`, angle-derived constraints, π dummies, epair dummies, etc., then `_pack_fixed_bonds()` to fill dense `(nAtoms, nMaxBonded)` arrays.  
   - **Node headless** calls [buildXPDBInputsFromXYZArgs()](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/dump_xpdb_topology.mjs:234:0-380:1) (under the hood it invokes the same Python MMFFL via Pyodide or pre-generated JSON), so it gets the exact same topology as the Python runner.  
   - **Browser GUI** instead calls `buildXPDBTopology(this.system, this.mmParams, { includeAngleConstraints, maxBonds:16 })` from `web/molgui_webgpu/XPDBTopology.js`, which is a lightweight JS rebuild and only supports real bonds + optional angles. It doesn’t add π/epair dummies at all, nor does it export the same dense bond matrix. So it’s expected that H₂O relaxes differently: you’re driving the solver with a different constraint set.

2. **Buffer preparation**  
   - Headless pipelines upload precalculated arrays (`bond_idx_global`, `bond_len_stiff`, etc.) matching the Fortran/MMFF data.  
   - GUI’s [uploadBonds(bondsAdj, defaultL, defaultK, angleDeg, enableAngles)](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/XPDB_WebGPU.js:194:4-267:5) simply loops through adjacency lists and writes lengths/stiffness from `bondsAdj`. Default L/K (0.4, 200) aren’t tied to the MMFF parameters. So even the rest lengths differ from what the headless run uses.

Given that, the plan should be:

- **Reuse the same topology generator in all three paths.** Export the Python/MMFFL topology (via JSON) and feed it to both the Node headless script and the browser. Alternatively, port the full MMFFL logic (including π/epair/dummy building) to JS and ensure both headless and GUI call that single implementation. Right now we have two independent builders.

- **Share the same bond packing routine.** Instead of the GUI’s `buildXPDBTopology` + [uploadBonds](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/XPDB_WebGPU.js:194:4-267:5), have the UI load the packed arrays (indices/lengths/stiffness) produced by MMFFL or by a shared JS module, so the GPU sees identical data. That includes honoring CLI-style flags (`twoPi`, `addPiAlign`, `addEpair`, etc.) and even `n_max_bonded`.

- **Only after topology parity is achieved** does it make sense to add distortion toggles or logging.

So the next concrete steps:

1. Extract the packed topology from the headless JSON (or a new API) and wire the GUI to ingest it instead of rebuilding locally.
2. Ensure [XPDB_WebGPU.step](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/XPDB_WebGPU.js:320:4-406:5) is always called with the options object (browser still uses positional args).
3. Once the GUI uses the same buffers, adding console dumps/download buttons becomes meaningful.

Without unifying the topology builder, the browser will keep showing different relaxations no matter how much we instrument it.

---

## USER

aha, again we have some stupid code redudancy making a mess.
like the reduduancy of buildXPDBTopology() and buildXPDBInputsFromXYZArgs(), are there any other similar casses of redudancy
 as I stressed several times, these test make only sense if we are checking same thing, so we should reuse the  same code, maximum part of the execution path MUST be shader between the headless script and the GUI, now when we see that headless is correct as it achive parity to pyOpenCL reference, we must modify the GUI version to reuse same execution path and function as the headless. we should make the old legacy GUI code deprecated (ad there at leas deprecation comments to be sure, we do not want to delete it immadiately just for case something breaks). So the whole XPDBTopology.js is now deprecatec / useless an insted of it we should use what?
we should also consider deprecated UI’s uploadBonds(bondsAdj, defaultL, defaultK, angleDeg, enableAngles) s and anythin like that and replace it by the alterantive function used by headless script. NOTE haowever, and this is very importaint, we should not reimplement anything (creating agaon code duplication) we must make sure that both headless and gui call same function from shader modules. If some functions are not in shared modules now, we should refactroe them and muve them there now. I wanted make them shared between GUI and headless test sript from the start, that was the purpose of testing, to test headlessly the functions used by the GUI, final goal is to make parity of the GUI to pythn, just testing it headlessly is more conveniet. Hope you undrstant the philosophy now.

make detailed notes and plan

---

## GPT-5.1-Codex

Observations on code duplication:

1. **Topology generation**
   - `buildXPDBTopology()` (GUI, JS) builds adjacency lists with minimal features (real bonds + optional angles) and default lengths.
   - [buildXPDBInputsFromXYZArgs()](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/dump_xpdb_topology.mjs:234:0-380:1) (headless Node) calls the same MMFFL Python code that [test_XPDB_new_dump_headless.py](cci:7://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/test_XPDB_new_dump_headless.py:0:0-0:0) uses, producing packed `bond_idx`, `bond_len`, `bond_k`, π/epair dummies, etc.
   - Result: browser uses a different constraint set than both headless paths.

2. **Bond uploads**
   - GUI calls [XPDB_WebGPU.uploadBonds(bondsAdj, defaultL, defaultK, angleDeg, enableAngles)](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/XPDB_WebGPU.js:204:4-277:5), which re-derives lengths/stiffness from defaults.
   - Headless uploads precomputed arrays via `writeBuffer` (exactly matching Python outputs).

3. **Solver parameter invocation**
   - Headless uses the new [xpdb.step({ ... })](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/XPDB_WebGPU.js:310:4-382:5) options object (so coll/bbox scales, dt, etc., stay in sync).
   - GUI still calls the positional signature, relying on legacy default margins.

Plan to eliminate redundancy and keep GUI/headless paths aligned:

1. **Centralize topology build in a shared module**  
   - Move the MMFFL-based packing logic into a reusable JS module (e.g., `XPDBSharedTopology.js`) that is imported by both [run_xpdb_webgpu_headless.mjs](cci:7://file:///home/prokop/git/FireCore/web/molgui_webgpu/run_xpdb_webgpu_headless.mjs:0:0-0:0) and [main.js](cci:7://file:///home/prokop/git/FireCore/web/molgui_webgpu/main.js:0:0-0:0).  
   - If maintaining Python MMFFL is required, expose a JSON export step that the GUI loads (e.g., precompute `topology.json` for each molecule).  
   - Mark `XPDBTopology.js` with `// DEPRECATED: use XPDBSharedTopology` and leave the old code commented with a TODO to remove once the new flow is stable.

2. **Replace [uploadBonds(bondsAdj, ...)](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/XPDB_WebGPU.js:204:4-277:5) with a shared buffer upload helper**  
   - Extract the headless upload routine (the part that writes packed arrays) into `XPDBSharedBuffers.js`.  
   - GUI code calls the same helper, passing the packed arrays generated by the shared topology builder. Deprecate the old [uploadBonds(bondsAdj, ...)](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/XPDB_WebGPU.js:204:4-277:5) by annotating it and routing it through the new helper internally until removal.

3. **Unify solver invocation**  
   - Update [main.js](cci:7://file:///home/prokop/git/FireCore/web/molgui_webgpu/main.js:0:0-0:0) to call [this.xpdb.step({ ... })](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/XPDB_WebGPU.js:310:4-382:5) with the same object used in the headless script.  
   - Store XPDB parameters (dt, inner_iters, k_coll, coll_scale, bbox_scale, init distortion) in a shared config module so both GUI controls and headless CLI read from the same defaults.

4. **Optional initial distortion**  
   - Implement [applyInitialDistortion()](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/run_xpdb_webgpu_headless.mjs:117:0-134:1) in a shared JS helper (`XPDBDebugUtils.js`). GUI can enable it via a checkbox for debugging parity, headless keeps the CLI flags.

5. **Add buffer dump parity**  
   - Expose [XPDB_WebGPU.dumpTypedState](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/XPDB_WebGPU.js:14:4-37:5) in GUI controls (with a download link) so browser output matches the headless `.txt` format. This helps verify that shared topology/upload code is behaving identically.

By moving topology packing, buffer uploads, and XPDB parameter defaults into shared modules (imported by both GUI and headless scripts), we ensure “maximum path sharing” without reimplementing anything. Deprecation comments should be added to `XPDBTopology.js`, [MolGUIApp.updateXPDBTopology](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/main.js:602:4-645:5)’s old code path, and [XPDB_WebGPU.uploadBonds](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/XPDB_WebGPU.js:204:4-277:5) so we remember to remove them once the shared modules are fully adopted.

---

## USER

make clear list which function are not in shared modules and need to be muved there and refactroed.
I hope you do not plan to make any new modules. We already have a lot of modules, and we should consolidate it, putting it to existing modules. Also mention which shared XPBD related function already are in shared modules.

---

## GPT-5.1-Codex



### XPDB functions that still live in one-off modules (must be moved/refactored into shared code)

| Area | Current function / file | Issue | Target shared home |
| --- | --- | --- | --- |
| Topology generation | `buildXPDBTopology(system, mmParams, …)` in `web/molgui_webgpu/XPDBTopology.js` | GUI-only light builder (real bonds + optional angles) reimplements what MMFF already gives us; no π / epair dummies, different rest lengths | Move logic to the same module that headless already uses ([dump_xpdb_topology.mjs](cci:7://file:///home/prokop/git/FireCore/web/molgui_webgpu/dump_xpdb_topology.mjs:0:0-0:0) → extract reusable `buildPackedTopology()` and import it from GUI); mark `XPDBTopology.js` as `// DEPRECATED` immediately |
| Bond upload | [XPDB_WebGPU.uploadBonds(bondsAdj, defaultL, defaultK, …)](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/XPDB_WebGPU.js:204:4-277:5) in [XPDB_WebGPU.js](cci:7://file:///home/prokop/git/FireCore/web/molgui_webgpu/XPDB_WebGPU.js:0:0-0:0) | Re-derives lengths/stiffness from defaults instead of consuming the packed arrays we already compute for headless | Replace GUI call sites with `XPDB_WebGPU.uploadPackedTopology(packed)` (new helper extracted from headless uploader). Keep [uploadBonds](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/XPDB_WebGPU.js:204:4-277:5) only as deprecated wrapper calling the shared helper |
| Solver invocation | GUI calls [this.xpdb.step(dt, iters, …)](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/XPDB_WebGPU.js:320:4-406:5) (positional) in [main.js](cci:7://file:///home/prokop/git/FireCore/web/molgui_webgpu/main.js:0:0-0:0) | Uses legacy defaults (coll/bbox scale 2.0) instead of the shared options object introduced for headless | Refactor GUI to call the same [step({dt,…})](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/XPDB_WebGPU.js:320:4-406:5) object used in headless; remove positional signature after transition |
| Initial distortion utilities | GUI has none; headless uses [applyInitialDistortion()](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/run_xpdb_webgpu_headless.mjs:113:0-130:1) defined inside [run_xpdb_webgpu_headless.mjs](cci:7://file:///home/prokop/git/FireCore/web/molgui_webgpu/run_xpdb_webgpu_headless.mjs:0:0-0:0) | Duplicate math will appear if we need it in the browser | Move the function into an existing shared helper (e.g., `web/molgui_webgpu/XPDBDebugUtils.js` or [dump_xpdb_topology.mjs](cci:7://file:///home/prokop/git/FireCore/web/molgui_webgpu/dump_xpdb_topology.mjs:0:0-0:0)) and import from both headless + GUI |

### XPDB functions already living in shared modules (no action needed)

- [XPDB_WebGPU.dumpTypedState()](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/XPDB_WebGPU.js:14:4-37:5) in [XPDB_WebGPU.js](cci:7://file:///home/prokop/git/FireCore/web/molgui_webgpu/XPDB_WebGPU.js:0:0-0:0): used by both Node runner and (once wired) the GUI to emit identical text dumps.
- [XPDB_WebGPU.readBuffersAsTyped()](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/XPDB_WebGPU.js:288:4-308:5) in [XPDB_WebGPU.js](cci:7://file:///home/prokop/git/FireCore/web/molgui_webgpu/XPDB_WebGPU.js:0:0-0:0): single implementation already shared; both GUI and headless call it.
- [buildXPDBInputsFromXYZArgs()](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/dump_xpdb_topology.mjs:234:0-380:1) in [web/molgui_webgpu/dump_xpdb_topology.mjs](cci:7://file:///home/prokop/git/FireCore/web/molgui_webgpu/dump_xpdb_topology.mjs:0:0-0:0): centralizes MMFFL topology build + buffer packing for headless (we just need GUI to reuse it instead of reimplementing).
- [XPDB_WebGPU.step(opts)](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/XPDB_WebGPU.js:320:4-406:5) in [XPDB_WebGPU.js](cci:7://file:///home/prokop/git/FireCore/web/molgui_webgpu/XPDB_WebGPU.js:0:0-0:0): new options-object path is shared; only GUI still uses the deprecated positional overload.

Once `buildXPDBTopology()` and [uploadBonds()](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/XPDB_WebGPU.js:204:4-277:5) are replaced with the shared builder/uploader, the GUI and headless paths will go through the identical modules, and the older functions can be marked `// DEPRECATED` and eventually removed.

---

## USER

you are an idiot. I just clearly said that buildXPDBTopology is deprecated and we should replace it by buildXPDBInputsFromXYZArgs()  used by the headless script. 
No need to refactor deprecated function, just mark them clearly as depecated in comments, we need to refactor into shareble modules only those which we wnat to use, which are those used by the headless which has partity with pyOpenCL. Do you readme or are you total ignoratn ?

---

Understood: we will stop touching `buildXPDBTopology()` and [XPDB_WebGPU.uploadBonds()](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/XPDB_WebGPU.js:194:4-267:5) except to add explicit `// DEPRECATED – use buildXPDBInputsFromXYZArgs()` comments, and we’ll focus solely on wiring the GUI through the same shared code path the headless runner already uses (i.e., [buildXPDBInputsFromXYZArgs()](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/dump_xpdb_topology.mjs:234:0-380:1) + packed buffers). No new modules—just relocate the existing headless helpers into the current shared files and have both GUI and CLI call them.


### XPDB Parity Migration Script (GUI ⇄ Headless ⇄ PyOpenCL)

1. **Declare legacy JS builder deprecated**
   - In `web/molgui_webgpu/XPDBTopology.js`, add a banner comment:
     ```js
     // DEPRECATED: topology is now generated by buildXPDBInputsFromXYZArgs().
     // This file exists only as a fallback while the GUI migrates to the shared path.
     ```
   - Remove all call sites except the GUI (temporarily); new code must not import this module.

2. **Reuse the headless topology helper**
   - Extract the already-correct helper from [web/molgui_webgpu/dump_xpdb_topology.mjs](cci:7://file:///home/prokop/git/FireCore/web/molgui_webgpu/dump_xpdb_topology.mjs:0:0-0:0):
     - [buildXPDBInputsFromXYZArgs({ ... })](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/dump_xpdb_topology.mjs:234:0-380:1) already produces packed arrays.
     - Export a pure function (e.g., `export async function buildPackedTopology(opts)`) from the same file—no new module needed.
   - In [main.js](cci:7://file:///home/prokop/git/FireCore/web/molgui_webgpu/main.js:0:0-0:0) (MolGUIApp), replace `buildXPDBTopology()` call with:
     ```js
     const packed = await buildXPDBInputsFromXYZArgs({
         xyz: currentXYZPath,
         // pass same flags GUI already exposes (pi, epairs, etc.)
     });
     ```
     Ensure the GUI supplies the same options as the headless CLI (twoPi, addPiAlign, enableAngles, etc.).

3. **Upload the exact packed buffers**
   - Add a helper inside [XPDB_WebGPU.js](cci:7://file:///home/prokop/git/FireCore/web/molgui_webgpu/XPDB_WebGPU.js:0:0-0:0):
     ```js
     uploadPackedTopology({ bondIndices, bondLenStiff }) {
         this.device.queue.writeBuffer(this.buffers.bondIdxGlobal, 0, bondIndices);
         this.device.queue.writeBuffer(this.buffers.bondLenStiff, 0, bondLenStiff);
     }
     ```
   - Mark [uploadBonds()](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/XPDB_WebGPU.js:204:4-277:5) with `// DEPRECATED: use uploadPackedTopology()` and keep it only for old code paths.
   - GUI now does `xpdb.uploadPackedTopology(packed)` exactly like the headless runner.

4. **Share solver invocation**
   - Update all GUI calls to [this.xpdb.step({ ... })](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/XPDB_WebGPU.js:320:4-406:5), passing the same object the headless script uses (`dt`, `iterations`, `k_coll`, `omega`, `momentum_beta`, `maxRadius`, `coll_scale`, `bbox_scale`).
   - Remove positional [step(dt, …)](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/XPDB_WebGPU.js:320:4-406:5) calls once no longer referenced.

5. **Optional debugging hooks**
   - Move [applyInitialDistortion()](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/run_xpdb_webgpu_headless.mjs:117:0-134:1) from [run_xpdb_webgpu_headless.mjs](cci:7://file:///home/prokop/git/FireCore/web/molgui_webgpu/run_xpdb_webgpu_headless.mjs:0:0-0:0) into an existing shared helper file (e.g., [dump_xpdb_topology.mjs](cci:7://file:///home/prokop/git/FireCore/web/molgui_webgpu/dump_xpdb_topology.mjs:0:0-0:0)) so GUI and CLI can import the same function when needed.
   - Add a “Dump XPDB state” button in the GUI that calls [XPDB_WebGPU.dumpTypedState](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/XPDB_WebGPU.js:14:4-37:5) and offers the `.txt` download; this matches the headless output format for diffing.

6. **Track remaining legacy code**
   - Add `TODO(DEPRECATED_XPDB)` tags to:
     - `XPDBTopology.js`
     - [XPDB_WebGPU.uploadBonds](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/XPDB_WebGPU.js:204:4-277:5)
     - Any GUI-specific helper that still calls them
   - Once the GUI is fully wired to the shared helpers, delete the deprecated implementations.

This script keeps everything inside the existing modules, ensures GUI + headless call the same functions, and makes the parity path explicit for future maintenance.





Remaining XPDB refactor work (shared vs legacy):

| Category | Legacy code (to be replaced / marked `// DEPRECATED`) | Shared code already correct |
| --- | --- | --- |
| Topology build | `web/molgui_webgpu/XPDBTopology.js` (`buildXPDBTopology`) — GUI-only builder using default lengths; must be replaced by calling [buildXPDBInputsFromXYZArgs()](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/dump_xpdb_topology.mjs:234:0-380:1) | [web/molgui_webgpu/dump_xpdb_topology.mjs](cci:7://file:///home/prokop/git/FireCore/web/molgui_webgpu/dump_xpdb_topology.mjs:0:0-0:0) ([buildXPDBInputsFromXYZArgs](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/dump_xpdb_topology.mjs:234:0-380:1)) — MMFFL-based packed topology used by headless + PyOpenCL |
| Bond upload | [XPDB_WebGPU.uploadBonds(bondsAdj, defaultL, defaultK, …)](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/XPDB_WebGPU.js:204:4-277:5) — rebuilds lengths per call; should be deprecated once GUI calls `uploadPackedTopology` | The headless flow already writes buffers via packed arrays (same helper logic we’ll expose as `uploadPackedTopology`) |
| Solver invocation | GUI still calls positional [xpdb.step(dt, iterations, …)](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/XPDB_WebGPU.js:320:4-406:5); needs to switch to the shared options-object call | [XPDB_WebGPU.step({ dt, iterations, k_coll, … })](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/XPDB_WebGPU.js:320:4-406:5) — used by Node headless and matches PyOpenCL settings |
| Initial distortion utilities | GUI lacks [applyInitialDistortion](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/run_xpdb_webgpu_headless.mjs:117:0-134:1); headless defines it inline (cannot be reused) | [applyInitialDistortion()](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/run_xpdb_webgpu_headless.mjs:117:0-134:1) in [run_xpdb_webgpu_headless.mjs](cci:7://file:///home/prokop/git/FireCore/web/molgui_webgpu/run_xpdb_webgpu_headless.mjs:0:0-0:0) — move into existing shared helper so both sides can import |

Once the GUI consumes [buildXPDBInputsFromXYZArgs()](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/dump_xpdb_topology.mjs:234:0-380:1) + `uploadPackedTopology()` and calls the shared [step({ … })](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/XPDB_WebGPU.js:320:4-406:5), all XPDB-critical code will live in the shared modules already proven against the PyOpenCL reference.



---

## USER


OK, the rendering of the topology looks good, but what is very bad is how the topology is recomputed
1) we are rebuilding all bonds. That should not be done. rebuilding bolds from geometry should be separate function, than building the axuliary dummy atoms and angle-derived bonds for the XPDB forcefied. Now it is mixed up and rebulding everything :-( this we need to repair and rationalize. When I load a .mol file which already contain bonds I do not want to rebuild bond every time I hit "Relax N steps" it should only rebuild the angle derived ones. 

2) Also when I deflect an atom and it rebuild the XPDB tpology, it does not seem to relax to something which loos chemically meaningfull. And I'm not sure why, but it can be because the angle-derived bonds are created from the actuall initial geometry, rather then from @MMParams.js ?

I mean in @MMParams.js should be optimal bond lengh for CH, OH etc, (every combination of types) and also stiffnesses read from @BondTypes.dat @AngleTypes.dat @AtomTypes.dat 

so for example any sp3 should be ~107deg, sp2 ~120deg etc, and C-C ~1.5A, and so we can compute by cosine law what should be distance of the two 2nd neighbor hydrogen atom in  H-C-H angle just form that params.

But it seems to me that instead what we take is the actuall geometry and distance between atoms in that actuall geometry. Which is pretty stupid because it will never correct initial distorted gometry.

Investigate this how it is done in the code and confirme aor disproof my hypothesis.   

---

## GPT-5.1-Codex

Key findings (matching the memories on “fail loudly” and “minimal dependencies”):

1. **We currently rebuild *all* bonds every time XPDB topology is rebuilt.**  
   - [MolGUIApp._cloneSystemForXPDB()](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/main.js:162:4-172:5) calls `m.recalculateBonds(this.mmParams)` unconditionally, so even MOL2 files that already contain bonds get re-derived from geometry every time you hit Relax/Step.  
   - [buildXPDBInputsFromMol()](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/MMFFLTopology.js:30:0-96:1) always calls [buildMMFFLTopology()](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/MMFFLTopology.js:589:0-719:1), which in turn:
     - Re-derives primary bond adjacency (`bondsAdj1`) using [mmParams.getBondParams()](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/MMParams.js:422:4-443:5) and the current molecule geometry.
     - Recomputes angle/auxiliary structures on top of that.  
   ⇒ You’re right: there’s no separation between “primary bond typing” and “auxiliary XPDB structures”; every rebuild retriggers the expensive, geometry-based bond finding.

2. **Angle-derived bond lengths currently ignore MMFF tables once packing happens.**  
   - Within [buildAngleBonds()](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/MMFFLTopology.js:154:0-184:1) we do use [MMParams](cci:2://file:///home/prokop/git/FireCore/web/molgui_webgpu/MMParams.js:62:0-522:1): `theta = atype.Ass`, `K = tB.Kss`, and `r_ab`, `r_bc` come from `bondsAdj1` whose `l0` was sourced from [mmParams.getBondParams()](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/MMParams.js:422:4-443:5). So the *initial* `linear_bonds` entries are consistent with AtomTypes + BondTypes + AngleTypes.  
   - However [buildXPDBInputsFromMol()](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/MMFFLTopology.js:30:0-96:1) throws away those `l0/k` values. It deduplicates `topo.bonds_*`, then reassigns **L = geometric distance** and **K = opts.defaultK** for *every* bond (even angle/π ones) @web/molgui_webgpu/MMFFLTopology.js#48-78. That means XPDB relaxations chase whatever distorted geometry you started with, not the canonical MM parameters—exactly the misbehavior you’re seeing.

3. **Python MMFFL reference keeps canonical lengths.**  
   - `pyBall/OCL/MMFFL._build_angle_bonds()` mirrors the JS logic but later, when packing for XPDB, Python keeps each `(l0, k)` from `linear_bonds` (see `_pack_fixed_bonds` in `pyBall/XPDB_AVBD/test_TiledJacobi_molecules.py`).  
   - JS diverged at the packing stage.

### Refactor direction

1. **Separate primary bond management vs XPDB augmentation:**
   - Persist the bonds that come from loaders (MOL/MOL2) and only call `recalculateBonds()` when the user explicitly requests it (e.g., “Recalculate Bonds” button).
   - Split [buildMMFFLTopology()](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/MMFFLTopology.js:589:0-719:1) into:
     1. A “primary snapshot” that takes existing `mol.bonds` (no mutation/recalc).
     2. An augmentation phase (angles, pi, epairs) that assumes those bonds are already correct.

2. **Preserve per-bond `l0`/`k` from MM tables:**
   - When building adjacency for XPDB, copy `(l0,k)` directly from:
     - `bondsAdj1` for real bonds (already typed via [mmParams.getBondParams()](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/MMParams.js:422:4-443:5)).
     - `linear_bonds` for angle/π/epair entries (they already store `l0`/`k`).
   - Only fall back to geometric distance + default `K` when MM data is missing.

3. **Optional next steps:**
   1. Add a checksum/version so XPDB rebuilds only when `mol.dirtyTopo` or aux options change.
   2. Surface a “Rebuild bonds from geometry” button for the rare case when the file truly has no bonds.

Let me know if you want me to draft the concrete JS changes (I’ll need Code mode to edit).


---

## USER

OK, so I was right, it is a mess. So correct that! we must take taget geometric params (angles, bond lenghts) from params not from actuall distorted geometry.
And we should not recalculate real bonds in the forcefield xPDB topology builder (we should recalculate bonds only in EditableBolecule when button "Recalculate Bonds" is pressed 
see recalculateBonds in @EditableMolecule.js , we should not have other redudcant bond recalculation.

the philospohy is that the @EditableMolecule.js should be source of real bonding topology and @XPDBTopology.js (or that is deprecated right? we should delete it?) or @MMFFLTopology.js should only add the axuliary pi and epair particles and angle-derived bonds on top of it so we can run XPBD.

Hope this is clearly understanble, Please make a clear plan.

---

## GPT-5.1-Codex

Plan to rationalize XPDB topology rebuilds (based on findings + user rules):

1. **EditableMolecule as sole primary-bond authority**
   - Stop calling `m.recalculateBonds()` inside [_cloneSystemForXPDB()](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/main.js:162:4-172:5) / [buildMMFFLTopology](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/MMFFLTopology.js:589:0-719:1). These routines should trust whatever `EditableMolecule` already has (MOL/MOL2 import, manual edits, or an explicit GUI “Recalculate Bonds” request).
   - Ensure all editor operations that mutate bonds (`addBond`, `deleteSelection`, script ops) already mark topology dirty—that’s where bond recalculation belongs, not in XPDB builders.

2. **Split MMFFL builder into two layers**
   1. **Primary snapshot**: read existing `mol.bonds` and copy their typed parameters (`l0`, `k`) directly via [mmParams.getBondParams()](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/MMParams.js:422:4-443:5); no geometry-based fallback unless a bond is missing MM data (only then use `opts.defaultL/K`, with a loud warn).
   2. **Auxiliary augmentation**: reuse the existing functions ([buildAngleBonds](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/MMFFLTopology.js:154:0-184:1), [buildPiDummies](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/MMFFLTopology.js:443:0-487:1), `buildEpair*`) but make them consume the primary adjacency that already contains canonical `l0/k`. They should append their own `(l0,k)` (using AtomTypes + AngleTypes) into `linear_bonds` without touching real bonds.

3. **Preserve canonical parameters during packing**
   - Change [buildXPDBInputsFromMol()](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/MMFFLTopology.js:30:0-96:1) so it doesn’t overwrite `l0/k` with geometric distances/default stiffness. Instead:
     - Seed adjacency with the typed primary bonds (`bondsAdj1`).
     - For each entry in `linear_bonds`, push its `l0/k` into adjacency as-is.
     - Only if a bond lacks a defined `l0` (e.g., missing data) compute/assign the geometric fallback once.
   - This mirrors Python’s `_pack_fixed_bonds`, letting XPDB relax distorted geometries back toward MMFF targets.

4. **Rationalize rebuild triggers**
   - XPDB rebuild should trigger only when:
     - EditableMolecule topology dirty (user changed bonds or pressed “Recalculate Bonds”).
     - XPDB aux options change (toggle pi/epair/angles, dummy lengths, stiffness sliders).
   - Implement a version counter or hash (`this.system.topoVersion` already exists) and store it alongside `_xpdbTopoLast`. If topoVersion hasn’t changed, skip the expensive rebuild and just reuse cached topology unless aux settings differ.
   - Ensure `markRawAllDirty/markRawTopologyDirty` only clears `_xpdbTopoLast` when actual molecule topology changes; pure geometry edits (moving atoms) shouldn’t invalidate MMFF topology.

5. **Cleanup / consistency**
   - Confirm `EditableMolecule.recalculateBonds()` is the *only* path that derives bonds from geometry. Update GUI hints: “Recalculate Bonds” should be the user-controlled action.
   - Remove the deprecated `XPDBTopology.js` or, at minimum, make it a thin compatibility shim that forwards to `MMFFLTopology` so there’s no confusion about multiple builders.
   - Add loud warnings when bond data is missing (`console.error` or thrown errors per “fail loudly” memory) so we don’t silently fall back to geometry.

6. **Implementation order**
   1. Refactor [_cloneSystemForXPDB](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/main.js:162:4-172:5) + [buildMMFFLTopology](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/MMFFLTopology.js:589:0-719:1) to stop recalc and to reuse existing bonds.
   2. Rework [buildXPDBInputsFromMol](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/MMFFLTopology.js:30:0-96:1) packing logic to preserve `l0/k`.
   3. Add cache invalidation keyed on `EditableMolecule.topoVersion` + aux-options signature.
   4. Delete or deprecate `XPDBTopology.js` to avoid drift.
   5. Add regression logging/tests ensuring primary bonds only change via GUI button.

Let me know when you want me to start editing (need Code mode).

---

## USER

"no geometry-based fallback" period !, I hate when you add stupid fallbacks which hide erros , whe have the rule to fail explicitly, we must throw exaction print error to console, never ad stupid fallback. we do scientific simulations here, not some stupid user webpage, we need to see errors when they happen ! This falback philosophy is a scourge of debugging.

1) when you implement these things add debug print to console (can be gated by verbosity) without them we are blind and can debung nothing !!!

we should also add these rules into memory, I keep saying theat but you still do not follow it enoguh !

update the plan

---

## GPT-5.1-Codex

Updated plan (no hidden fallbacks, with debug visibility):

1. **EditableMolecule owns real bonds exclusively**
   - Strip all implicit `recalculateBonds()` calls from [_cloneSystemForXPDB()](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/main.js:174:4-184:5)/[buildMMFFLTopology()](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/MMFFLTopology.js:589:0-719:1). XPDB code must trust whatever bonds already exist; only the GUI “Recalculate Bonds” button may rebuild them.
   - Ensure any editor action that mutates bonds keeps the existing [markRawTopologyDirty](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/main.js:1242:4-1245:5) flow, so XPDB sees a new topoVersion instead of recomputing bonds itself.

2. **MMFFL builder = augmentation layer only, with loud errors**
   - While ingesting `mol.bonds`, resolve `(l0,k)` via [mmParams.getBondParams()](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/MMParams.js:422:4-443:5); if lookup fails, immediately throw (or `console.error` + throw) rather than guessing. No geometry fallback.
   - Angle/pi/epair helpers must also error when required AtomType/AngleType data is missing. Every such path gets verbosity-gated debug prints (e.g., `if ((window.VERBOSITY_LEVEL|0)>=2) console.log('[XPDB/MMFFL] ...')`) so we can trace inputs when it fails.

3. **XPDB packing preserves canonical parameters**
   - [buildXPDBInputsFromMol()](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/MMFFLTopology.js:30:0-96:1) seeds adjacency directly from the typed primary bonds plus the `linear_bonds` list. Each entry carries its original `l0`/`k`; if either is NaN/<=0, throw with a descriptive message (atom indices, tag).
   - No geometry-derived lengths unless explicitly provided by MM tables; instead, fail loudly so missing data surfaces in dev logs.

4. **Controlled rebuild triggers + diagnostics**
   - Cache `_xpdbTopoLast` with `{topoVersion, auxSignature}`. Rebuild only when either changes; otherwise reuse the cached topo, logging at high verbosity.
   - [markRawTopologyDirty](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/main.js:1242:4-1245:5) should clear the cache only when `EditableMolecule.topoVersion` increments. Pure geometry moves shouldn’t nuke the XPDB topo.
   - Add debug prints (verbosity ≥2) in [_rebuildXPDBTopoCPU](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/main.js:87:4-133:5), [_refreshXPDBTopoOverlay](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/main.js:135:4-172:5), and XPDB packing so we always see reason codes, n_real/n_all, bond counts.

5. **Cleanup + documentation**
   - Deprecate `XPDBTopology.js` (either delete or make it delegate to the new path with a loud warning) to avoid parallel implementations.
   - Document these “no fallback / throw on missing data / debug prints” rules in [migration_to_webgpu.md](cci:7://file:///home/prokop/git/FireCore/web/molgui_webgpu/migration_to_webgpu.md:0:0-0:0) so future contributors follow them.

All future changes will follow the “fail loudly” memory (beab66ed-e401-467a-a73d-25e9a5186c6f) and add verbosity-gated diagnostics per your instruction. When you’re ready for code changes, switch to Code mode.

---

## USER

OK, is seems quite good now,seems to work, but we have probably problem with types, we need to make clear debug ouput in console what types and parameters are assinge, I opened petacene which has all sp2 carbons (aromatics) but the result seems more like sp3 (aliphatic carbons). Also the pi-orbitals are created just on one side. We should have somewhere in gui clickable options which axulary topology (pi-orbitals, electrob pairs, and derived bonds from angles/pi/epair) should be generated, nd if pi-orbitals should be two-sided, we have similar options in the headless scrips like @dump_xpdb_topology.mjs @run_xpdb_webgpu_headless.mjs 

also we should ad option to visialuze bounding boxes

make clear plan.

---

## GPT-5.1-Codex

### Updated plan

1. **High-verbosity type/parameter reporting**  
   - In [buildMMFFLTopology()](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/MMFFLTopology.js:626:0-773:1) (and headless parity via [dump_xpdb_topology.mjs](cci:7://file:///home/prokop/git/FireCore/web/molgui_webgpu/dump_xpdb_topology.mjs:0:0-0:0)), log per-atom assignments when `VERBOSITY_LEVEL>=2`: atom id, symbol, AtomTypes.dat key, valence, `Ass`, `Kss`, bond params used (`l0/k`).  
   - For pentacene we’ll see whether carbons get `C_R`/`C_ene` (sp2) vs `C_3` (sp3). Missing entries will throw already; now we also see what was chosen.

2. **GUI + headless toggles for auxiliary topology & pi dummies**  
   - Add a small XPDB options panel (checkboxes) in the Web GUI mirroring the CLI flags (`addAngle`, `addPi`, `twoPi`, `addPiAlign`, `addEpair`, `addEpairPairs`).  
   - Each toggle updates an internal “XPDB opts” struct and forces `_xpdbTopoLast` rebuild.  
   - Expose the same options through [run_xpdb_webgpu_headless.mjs](cci:7://file:///home/prokop/git/FireCore/web/molgui_webgpu/run_xpdb_webgpu_headless.mjs:0:0-0:0) CLI (already parsed) and propagate to [buildXPDBInputsFromMol](cci:1://file:///home/prokop/git/FireCore/web/molgui_webgpu/MMFFLTopology.js:30:0-123:1).  
   - Add a “show bounding boxes” checkbox to Visualization: when enabled, call the existing bounding-box helper (or add one) in the renderer.

3. **Files touched so far (to stage later)**  
   - [web/molgui_webgpu/main.js](cci:7://file:///home/prokop/git/FireCore/web/molgui_webgpu/main.js:0:0-0:0)  
   - [web/molgui_webgpu/MMFFLTopology.js](cci:7://file:///home/prokop/git/FireCore/web/molgui_webgpu/MMFFLTopology.js:0:0-0:0)

Let me know when you want me to start editing (switch to Code mode).