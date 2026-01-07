# MolGUI Web – Design & TODO

This document is a **clean design + TODO** summary for the MolGUI web editor.
The older `MolGUI_web.md` file contains the original Gemini discussion and is
kept for reference only.

## 1. High‑Level TODO Checklist

### 1.1 Core Viewer/Editor (Implemented)

- [x] **Infrastructure & Logging**
- [x] **Data & Renderer** – `MoleculeSystem` + `MeshRenderer` (instanced impostors, shared position texture)
- [x] **Input / Output (XYZ)** – load/save + inline editor
- [x] **Selection** – click + box, additive/subtractive, highlight
- [x] **Gizmo & Manipulation** – `THREE.TransformControls` on proxy object
- [x] **Labels** – instanced text labels (ID / element / type placeholder)
- [x] **MM Parameters** – `ElementTypes.dat` + `AtomTypes.dat` integration
- [x] **Two-system model (molecule + substrate)** – independent `EditableMolecule` + `MoleculeRenderer` instances, no `main` system
- [x] **ScriptRunner handles** – direct `molecule`/`substrate`/`mol` handles for `clear()`, `load()`, `build_substrate()`, `move/rotate/replicate()`
- [x] **EditableMolecule.lvec + replicate** – per-instance lattice vectors, physical PBC replication (atoms+bonds+lvec)
- [x] **Instanced rendering** – atoms as instanced impostor quads, bonds as line segments, shared shaders; replication clones visual meshes per-lattice without duplicating data

### 1.2 Topology / Chemistry (Planned or Partial)

- [ ] **Connected Components / Fragments**
    - Compute connected components on the bond graph (`findConnectedComponents`).
    - Expose operations: select / hide / delete fragment.
- [ ] **Ring / Cycle Detection**
    - Port `atomicUtils.find_cycles` → JS.
    - Optional visualization of rings.
- [ ] **Bridges / Cut Bonds**
    - Detect bonds whose removal splits a component.
    - Optional highlighting.
- [ ] **Group Substitution / Functional Groups**
- [ ] **Valence / Electron‑Pair Aids**

### 1.3 Interaction & Editing Features (Planned)

- [ ] **Selection Modes: Atom / Molecule / Group**
    - Use connected‑components result to define "molecule" selections.
    - Mode switch: single atom vs whole molecule vs named group.
- [ ] **Soft Selection / Falloff Editing**
- [ ] **Explicit Pivot Control**
- [ ] **Axis‑Based Alignment (by 2–4 picked atoms)**
- [ ] **PCA / Inertia‑Tensor Alignment**
- [ ] **Savable Selections / Groups**

### 1.4 Polymer-on-Surface Editor (Planned)

- [ ] **Crystallography Engine** – unit cell import, surface cleavage, step/terrace generator.
    - [x] Unit cell import from CIF (preset + file load)
    - [x] Optional symmetry expansion (CIF `_symmetry_equiv_pos_as_xyz`)
    - [x] Periodic crystal replication (nx,ny,nz) + optional Miller alignment `(h,k,l)->z`
    - [x] Periodic bonds from `BondTypes.dat` for full crystals (basis + 26 neighbor images)
    - [x] Miller slab cut by `cmin/cmax` in Å (cell-overlap pruning)
    - [ ] Slab bonds (slab cut currently disables bond generation)
    - [ ] Optional in-plane twist control ("up" vector) after aligning normal to Z
    - [ ] Bond order/type propagation from `BondTypes.dat`
- [ ] **Lattice Matcher** – commensurability solver and visual debugger along chosen surface directions.
- [ ] **Monomer Library & Sequencer** – monomer prefabs with head/tail/up anchors, sequence editor for heterogeneous chains.
    - [x] Example monomer preset set (from `tests/tAttach/polymerize.py`) + sequence parser/tokenizer
    - [x] Polymer sequence build (append monomers + join bonds)
    - [x] Repeat-vector sign auto-selection (choose ±`lvec[1]` to minimize head↔tail distance; fixes long wrap-around join bonds)
    - [x] Debug helpers in GUI: save MOL2 + log longest bonds
    - [ ] Full monomer/endgroup/backbone library JSON workflow (curated, not only presets)
- [ ] **Curve & Frame System** – spline paths + robust frame transport along curves/step edges.
- [ ] **Polymer Assembler** – place and orient monomers along path (twist/tilt controls, instancing for performance).

### 1.5 Polymer attachment / substitution (Implemented for MVP)

- [x] **Marker-based group attachment** (ported from legacy, implemented in `EditableMolecule.attachGroupByMarker`)
    - supports distinct marker pairs for backbone vs group (e.g. `Se/Cl` backbone + `Al/Cl` group)
- [x] **Direction-based attachment** (cap+back atom with up-vector + twist) via `EditableMolecule.attachParsedByDirection`
- [x] **Examples UI** (from `tests/tAttach/*`) to load backbones/endgroups and run marker attach in one click

### 2. Current Architecture Overview

### 2.1 HTML & Entry Point

- `web/molgui_web/index.html`
  - Declares `#canvas-container`, help overlay, loads Three.js + controls.
  - Single ES module entrypoint: `js/main.js`.

- `web/molgui_web/js/main.js`
  - Sets up scene, orthographic camera, renderer, `OrbitControls`.
  - Loads shaders from `../common_resources/shaders/` (atoms, bonds, selection, labels).
  - Instantiates two authoritative systems only: `molecule` and `substrate` (each with its own `EditableMolecule`, `PackedMolecule`, `MoleculeRenderer`); no `main`.
  - `Editor` (selection + gizmo), `GUI` (sidebar), `ShortcutManager`, `ScriptRunner`.

### 2.1.1 Controls / UX (Target Scheme)

- **Mouse**
  - **LMB**: selection and gizmo interaction (no camera on left click).
  - **RMB**: rotate camera (OrbitControls RIGHT = ROTATE).
  - **MMB**: pan camera (OrbitControls MIDDLE = PAN), with optional Shift+RMB pan fallback.
  - **Wheel**: zoom.
- **Selection modifiers**
  - No modifier: replace selection.
  - Shift: add to selection.
  - Ctrl/Alt: subtract or toggle.
- **Keyboard shortcuts** (current and planned)
  - `G`: toggle gizmo, `T/R/S`: gizmo mode (translate/rotate/scale).
  - `Esc`: clear selection, `Del/Backspace`: delete selection.
  - `B`: recalc bonds, `A`: add atom.
  - Future: keys to cycle selection mode (atom/molecule/group), toggle soft selection, store/recall groups.

### 2.2 Shared JS Modules (`web/common_js`)

- `Logger.js`
  - Global `logger` with console + DOM output.

- `Draw3D.js`
  - Generic GPU helpers:
    - Create instanced impostor meshes (`createTextureBasedInstancedMesh`).
    - Create line segments driven by position texture.
    - Label system: font atlas, instanced label mesh, `updateLabelBuffers`.

- `MeshRenderer.js`
  - Generic renderer that owns:
    - Position texture (`uPosTex`) as **single source of truth**.
    - Atom mesh (instanced quads, sphere impostor shaders).
    - Selection mesh (instanced highlights).
    - Bond lines (texture‑driven line segments).
    - Label mesh (instanced glyph quads).
  - API:
    - `updatePositions(posArray, count)`
    - `updateParticles(count, colorGetter, scaleGetter)`
    - `updateSelection(indices)`
    - `updateBonds(pairs)`
    - `updateLabels(stringGetter, count)`
    - `setLabelStyle`, `setNodeScale`, `setSelectionScale`, visibility toggles.

- `Selection.js`
  - Generic selection containers (`Selection`, `SelectionBanks`) used for sets of indices.

- `GUIutils.js`, `MeshesUV.js`, `SDfuncs.js`, `Vec3.js`
  - Utility math / mesh helpers (available for future extensions).

### 2.3 MolGUI‑Specific JS (`web/molgui_web/js`)

- `MoleculeSystem.js`
  - SoA storage: `pos: Float32Array`, `types: Uint8Array`.
  - `bonds: [i,j]` array + `neighborList` (adjacency).
  - Selection state: `selection: Set<int>`.
  - Key methods:
    - `addAtom(x,y,z,type)`
    - `addBond(i,j)`
    - `updateNeighborList()`
    - `recalculateBonds(mmParams)` – distance cutoffs from MM covalent radii.
    - `deleteSelectedAtoms()` – rebuilds arrays and bond indices.

- `MoleculeRenderer.js`
  - Thin wrapper around `MeshRenderer` for molecules.
  - Maps atom types → colors/radii using `MMParams`.
  - Controls label mode (ID / element / type).

- `MMParams.js`
  - Parses `ElementTypes.dat`, `AtomTypes.dat`, and `BondTypes.dat` from `common_resources`.
  - Builds lookups:
    - `elementTypes` (by name), `byAtomicNumber[iZ]`, `atomTypes`.
  - Provides bond-length lookup from `BondTypes.dat` (used for crystal bond building).

- `Editor.js`
  - Handles pointer events for:
    - Click picking (ray–sphere intersection per atom).
    - Box selection (project to screen, test in rectangle).
    - Managing `TransformControls` gizmo (centroid pivot, drag updates positions).
  - Coordinates with `MoleculeSystem` + `MoleculeRenderer`.

- `IO.js`
  - Deprecated (XYZ parsing/export lives in `EditableMolecule.js`).

- `GUI.js`
  - Builds sidebar:
    - Selection info, view controls, gizmo toggles, structure controls, geometry load/save, parameter editors, log panel.

## 2.4 Recent changes (Dec 2025)

- **Attachment functions migrated to `EditableMolecule`**
  - marker attach + direction attach implemented (Vec3/Mat3, stable atom IDs)
- **Polymer builder stabilized**
  - sequence parsing fixed (`DDDD_DDDD` expands to letters; `PNA10` works)
  - repeat translation sign auto-chosen by head↔tail distance (fixes wrong long join bonds)
- **Debuggability improved**
  - added MOL2 export
  - added “log longest bonds” helper to detect topology mistakes

- **Faster bond rebuild via spatial buckets (with robust topology handling)**
  - bond modes: `brute`, `bucketNeighbors`, `bucketAllPairsAABB` (GUI: Structure → Bond mode)
  - bucket graph lives in `window.app.lastBucketGraph` (built by crystal/substrate generators)
  - bucket debug overlay:
    - boxes + atom→bucket-center lines (GUI toggles in Structure)
    - auto-update after topology changes (e.g. deletion)
  - robust storage across swap-remove deletions:
    - bucket graph supports flipping atom storage between indices and stable atom IDs: `BucketGraph.toIds(mol)` / `BucketGraph.toInds(mol)`
    - after deletions, empty buckets are removed: `BucketGraph.pruneEmptyBuckets()`
    - bounds are recomputed for visualization: `BucketGraph.recalcBounds(mol)`
  - bond rebuild logs now include `nAtoms` and `nBuckets`
  - key files:
    - `web/common_js/Buckets.js` (`BucketGraph`)
    - `web/molgui_web/js/Editor.js` (deleteSelection + recalcBonds integration)
    - `web/molgui_web/js/main.js` (bucket overlay + refreshBucketDebug)

- **Valence-based passivation (caps + explicit electron pairs)**
  - `EditableMolecule.addCappingAtoms(mmParams, cap='H', opts)`
    - adds missing sigma-bond caps based on `AtomTypes.dat` valence (`at.valence`) and current coordination
    - robust with explicit epair dummy atoms (`Z==200`): epairs count as geometry domains but do not consume sigma valence for capping
  - `EditableMolecule.addExplicitEPairs(mmParams, opts)`
    - inserts explicit lone-pair dummy atoms using `at.nepair` and `at.epair_name`
  - GUI: Structure → Passivation (Cap type textbox + buttons Add Caps / Add EPairs)
  - key file: `web/molgui_web/js/EditableMolecule.js`
  - design justification:
    - keep chemistry metadata in `MMParams` (`ElementTypes.dat`, `AtomTypes.dat`), avoid hardcoded valence tables
    - keep edits in the authoritative graph model (`EditableMolecule`) and re-export to renderer

- **Selection by chemical environment (AND-only query language)**
  - compiled query (parse once) + fast apply over existing per-atom bond lists
  - API:
    - `EditableMolecule.compileSelectQuery(q, mmParams)`
    - `EditableMolecule.applySelectQuery(compiled, {mode:'replace'|'add'|'subtract'})`
  - syntax:
    - atom set: `N|C` (elements) or `C_3|O_OH` (atom types, `_` implies atom-type)
    - neighbor counts: `n{F|Br|Cl}={1,2}` and wildcard `n{*}={1,2}`
    - shorthand: `deg={1,2}`
  - GUI: Selection → Selection Query (Replace/Add/Subtract)
  - design justification:
    - no separate neighbor list is needed: `Atom.bonds[]` already encodes adjacency; we only compile the string into predicates once

- `ShortcutManager.js`
  - Global keyboard shortcuts (gizmo toggle/mode, delete, add atom, recalc bonds).

### 2.4.1 Recent changes (Jan 2026)

- Switched to **two-system model** (`molecule`, `substrate`), removed `main` system; each has its own renderer.
- **ScriptRunner handles**: `molecule`, `substrate`, alias `mol`; methods `clear`, `load`, `build_substrate`, `move/translate`, `rotate/roate`, `replicate/replication`. Commands accept optional `system` without changing global state.
- **EditableMolecule**: per-instance `lvec` and `replicate(nrep,lvec)` clone atoms/bonds and scale lattice; rotation uses `Mat3.fromAxisAngle` and guards missing positions.
- **Default user script**: clears both systems, builds substrate, loads/positions molecule, replicates each; no `main`/merge.
- **Replica rendering refactor**:
    - Retained original GPU-friendly draw-call replication (shared geometry) but added `replicaClones` tracking so every cloned `THREE.InstancedMesh` mirrors the source `count`/matrices.
    - `_syncReplicaInstancedMeshes()` runs after structure changes or replica rebuilds; fixes “ghost” replicas after `clear()` or resizing.
    - Bucket/lattice visuals stay encapsulated in `MoleculeRenderer`; no stray groups in `main.js`.
- **Script + GUI lattice pipeline**:
    - Added `ScriptRunner.addLvec()` / `setViewReplicas()` commands, accepting flexible array/object inputs (`_vecFromInput`).
    - `MolGUI_web/js/main.js::updateReplicas` now logs every invocation and targets `this.renderers[name]`, guaranteeing molecule/substrate parity.
    - GUI sample (“PTCDA on NaCl step”) demonstrates periodic molecule setup via scripts (lvec + view replicas).
- **GUI defaults**: Replica / Lattice panel now defaults to the *molecule* system so scripts and manual edits affect the molecule unless explicitly switched.
- **Rendering/replication**: atoms are instanced impostor quads, bonds are line segments; replication clones visual meshes per lattice shift without duplicating data.

#### Pitfalls & takeaways
- Pure CPU duplication of atom positions (InstancedMesh-of-InstancedMesh) hurt performance and broke user expectations—stick to single draw-call + transform groups when the scene already uses instanced impostors.
- Forgetting to sync cloned instanced meshes (`count`, `instanceMatrix`) leaves “ghost” copies after `molecule.clear()`; always track clones explicitly.
- ScriptRunner commands must be registered in the whitelist **and** support ergonomic inputs; otherwise the script editor silently fails. Use helper parsers to normalize user-provided arrays/objects before touching `Vec3.setV`.
- System defaults matter: default the Replica/Lattice GUI dropdown to *molecule* so script-driven updates and manual edits align with user expectations.

### 2.4.1 Recent changes (Jan 2026)

- **ShortcutManager.js**
  - Global keyboard shortcuts (gizmo toggle/mode, delete, add atom, recalc bonds).

### 3. Planned Topology / Chemistry Features (Details)

This section expands 1.2 with more implementation hints.

### 3.1 Connected Components / Fragments

- Core algorithm: graph traversal (DFS/BFS) over `bonds` to label components.
- Store component ID per atom; cache until bonds change.
- Integrate with selection:
  - In "molecule" selection mode, clicking any atom/bond selects all atoms with the same component ID.

### 3.2 Rings, Bridges, Groups, Valence

- **Rings:** port `atomicUtils.find_cycles`.
- **Bridges:** use graph algorithms on the bond graph (edge‑connectivity).
- **Groups / substitution:** define templates and use alignment utilities (see 4.3) to place them.
- **Valence aids:** quick check of coordination vs expected valence (from `VALENCE_DICT` etc.).

## 4. Planned Interaction Features (Details)

### 4.1 Selection Modes & Savable Groups

- Mode flag in editor: `mode = 'atom' | 'molecule' | 'group'`.
- **Atom mode:** current behavior.
- **Molecule mode:**
  - Require connected‑components (3.1).
  - Picking atom/bond → replace selection with that component.
- **Group mode:**
  - Named groups: maps `name -> Set<atomID>`.
  - GUI list for creating / renaming / deleting groups.
  - Operations: select group, add/remove current selection to group.

### 4.2 Soft Selection / Falloff

- Weight function `w(i)` based on:
  - Geometric distance from pivot atom, or
  - Graph distance (number of bond steps).
- During gizmo drag, apply weighted displacement/rotation.
- GUI parameters:
  - Radius / max graph distance.
  - Falloff curve (e.g. Gaussian / linear).

### 4.3 Pivot & Alignment

- **Pivot:**
  - Maintain explicit pivot object (separate from selection centroid).
  - UI to set pivot to: atom, bond midpoint, arbitrary coordinates.

- **Axis alignment by atoms:**
  - User picks: origin atom, forward pair, up pair.
  - Build orthonormal basis and transform selection.

- **PCA alignment:**
  - Compute covariance of selected atom positions.
  - Eigenvectors = principal axes.
  - Rotate so axes align with global XYZ; use XY for planar molecules.

## 5. Polymer-on-Surface Editor Features (Details)

These features extend the editor towards polymer assembly along crystalline surfaces and step edges.
For full derivations and pseudocode see:

- `doc/Editor_for_Polymers_on_Surfaces_and_Edges.md`

### 5.1 Crystallography Engine

- Unit cell import (`.cif`, `.xyz`): lattice vectors and basis atoms.
  - CIF: preset + file load (symmetry optional)
- Surface cleavage by Miller indices `(h,k,l)` to generate slabs.
  - Current implementation: slab cut by `cmin/cmax` in Å along `nHat` (from reciprocal lattice)
  - Current limitation: slab cut cannot build bonds yet
- Step & terrace (vicinal) surface generator:
  - User controls terrace width and step height in unit cells.
  - Output keeps semantic info (terrace, step edge, lower terrace) for later selection.

### 5.2 Lattice Matcher (Commensurability)

- User picks a substrate direction vector along a surface/step.
- Compute integer pairs `(N, M)` such that `N * |V_sub| ≈ M * L_mono` within strain tolerance.
- Provide a simple visual debugger:
  - Ruler along the chosen direction.
  - Tick marks for repeating positions, ghost preview of monomer backbone.
  - For algorithmic details, see existing 2D lattice matching tools:
    - `cpp/common/molecular/LatticeMatch2D.h`
    - `cpp/libs/Molecular/Lattice2D_lib.cpp`
    - `pyBall/Lattice2D.py`
    - `tests/tLattice2D/run.py`

### 5.3 Monomer Library & Sequencer

- Define monomer prefabs with:
  - Head (entry) anchor, tail (exit) anchor, and up vector.
- Sequence editor UI for heterogeneous chains (e.g. `A-A-A-B-A`), independent of biopolymer conventions.

### 5.4 Curve & Frame System

- User defines control points along step edges or arbitrary paths.
- Generate smooth curves (Catmull–Rom / cubic Bézier).
- Robust frame transport along the curve:
  - Avoid unstable Frenet frames; use rotation-minimizing (parallel transport) or user-up interpolation.
  - Output a transform (position + rotation) per sample along the path.

### 5.5 Polymer Assembler

- Iterate along the curve and sequence:
  - Fetch monomer prefab for each position.
  - Align monomer head→tail vector to curve tangent.
  - Align monomer up vector to curve/surface normal.
  - Apply user twist/tilt controls per monomer or globally.
- Use instanced rendering when repeating many monomers for performance.

## 6. Performance & Coding Guidelines

These are the key rules to keep the implementation scalable while remaining debuggable:

- **No allocations in hot loops**
  - Avoid `new THREE.Vector3()`, `new Float32Array()`, string building, etc. inside per-frame or per-atom loops.
  - Reuse temporaries (`this.tempVec`, shared arrays) and update in place.

- **TypedArrays everywhere for bulk data**
  - Atom positions, types, and other large arrays should live in `Float32Array` / `Uint8Array` buffers.
  - Use SoA layout and resize strategies similar to C++ `std::vector` (double capacity on overflow).

- **Single source of truth on GPU**
  - Keep atomic positions in a single DataTexture (`uPosTex`).
  - Atoms, bonds, selection, labels all fetch from this texture instead of duplicating position buffers.

- **Dirty flags and minimal updates**
  - Use a dirty flag on `MoleculeSystem` to decide when to push data to GPU.
  - For small edits, update the minimal range; avoid full rebuilds unless bonds/topology change.

- **String / DOM work off hot paths**
  - Logger and GUI updates should not run in tight loops over atoms.
  - Use `Logger` with verbosity controls to disable heavy debug output in production.

## 7. Python Reference Toolset (pyBall)

Many topology and editing operations already exist in Python and should be **ported or mirrored**, not re‑invented:

- `pyBall/atomicUtils.py`
- `pyBall/AtomicSystem.py`

Key areas to mirror:

- Bond finding (`findBondsNP`, etc.).
- Neighbor lists and graph utilities.
- Group / fragment detection.
- Cycle/ring detection.
- Orientation and alignment helpers.
- Valence and electron‑pair heuristics.

Keeping MolGUI Web aligned with these tools makes it easier to share workflows
between Python (desktop) and the browser version.
