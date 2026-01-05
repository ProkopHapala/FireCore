## USER 1

I have 3D editor of molecules in javascript / webGL. What is specific about my use-case is that I focus on simulation of molecules on surfaces of crystaline solids, mostly ionc crystals. I would like to have there functions allowing me easily visual generation of these surfaces, and placement of molecule onto it, aligned with respect to the unit-cell r supercell of the crystal. In particular I'm interested in step edges. I should select some crystalocaphic plane or orientation ad make edge on the substrate along it.  I also want to do polymerization along these step edges. So I then want to chose this direction (or other direction) and compute some lattice mathich (small multiples where the lattice match, to make it comensurable). I want to debug this visually. Then I want to take my molecules, align rotate them by two angles with respect to that line, and repat them along that edge forming polymer. This is the task I want to be easy to create in this editor. Also I want to create heterogenouse structures - where the polymer ic copolymer of several different monomers (like sequence of letters, nucleobases in DNA or aminoacids in protein (but I do not work with biopolymers I sesing my own which different backbone and different sidechains)). I want to have general system allowing me to build arbitrary sequence of molecules into such chain, and then orient them. As the guide-line can be used even some cureve - where I have the control points for position (e.g. using cubic spline, and also the up-vectro). So basically I will generate array of points, and upvectros (resp. rotation matrix at each such point), and than translate and rotate the monomers to follow that line. 

Now can you take my chaotic ideas, order them, and make them ito more systematic design document? (listing all the fatures we want to implement, how they should be connected and what mathematicall of computer-sciece tools or algoruthms we can use for it ?)


---


Here is a structured Design Document based on your requirements. This transforms your ideas into a software architecture suitable for a JavaScript/WebGL environment (e.g., using Three.js or Babylon.js).

---

# Design Document: Surface-Supported Polymer Assembly System

## 1. Project Overview
**Goal:** To extend a 3D web-based molecular editor to facilitate the simulation of heterogeneous copolymers on crystalline surfaces, specifically focusing on step-edge alignment, lattice commensurability, and curve-guided placement.

**Core Philosophy:**
1.  **Substrate-First:** The crystal geometry dictates the layout.
2.  **Parametric Assembly:** Polymers are defined by sequences and paths, not just manual atom dragging.
3.  **Visual Debugging:** Immediate visual feedback for geometric constraints (lattice matching, steric clashes).

## 1.1 Implementation decomposition (agreed plan)

We split implementation into 3 mostly-independent subtasks. This is important because only the last step requires chemically correct topology editing.

### Subtask 1: Surface / lattice replication (pure geometry, no topology)

Goal: generate a substrate slab/surface by copying unit-cell atoms along lattice vectors.

- **Inputs**
  - lattice vectors (unit cell) `a,b,c` (3D)
  - basis atoms in the unit cell (positions + element/type)
  - replication counts `(na, nb, nc)` (integers)
- **Output**
  - atom list in world Cartesian coordinates
  - (optional) per-atom tag/group, e.g. `"substrate"`, to allow locking selection and separate rendering
- **Notes**
  - this step does **not** require bonds or neighbor lists
  - start with simple replication; Miller-plane cuts and step-terraces can come later

### Subtask 2: Rigid-body placement/orientation (rotation + translation)

Goal: place molecules/monomers on the substrate and align them along a chosen edge direction.

- **Core operation**
  - compute a rotation matrix (or quaternion) from:
    - **forward vector** (e.g. step-edge direction, or a bond direction)
    - **up vector** (e.g. surface normal)
  - apply rotation + translation to a set of atoms (a monomer prefab, an end-group, or a selected fragment)
- **Reuse**
  - the same math is used for:
    - orienting monomers along the step edge
    - attaching end-groups (chemical residues) to a bond: bond vector = forward vector
    - curve/spline frame alignment later

### Subtask 3: Topology operations (bond editing + atom replacement)

Goal: correctly form polymer chains and attach groups by editing bonds and possibly removing atoms (e.g. replacing H/Cl with a group).

- **Operations needed**
  - add bonds
  - remove bonds
  - delete atoms (e.g. remove leaving group / hydrogen)
  - keep neighbor lists consistent
- **Data structure trade-off (explicit decision)**
  - **arrays with integer indices** are fast (GPU-friendly, good for rendering)
  - but **index stability is fragile** when deleting atoms (reindexing required)
  - a dictionary of unique IDs would simplify logic but is slower and more complex
  - therefore: keep arrays as the core representation, but design helpers for reindexing + stable labels

---

## 2. Functional Modules

### Module A: The Crystallography Engine
Responsible for generating the "stage" upon which chemistry happens.

*   **Unit Cell Import:** Parsers for `.cif` or `.xyz` files to define lattice vectors ($\vec{a}, \vec{b}, \vec{c}$) and basis atoms.
*   **MVP: Replication-only surface generation (Subtask 1):**
    *   Input: lattice vectors + basis atoms + `(na,nb,nc)`.
    *   Output: replicated atoms in Cartesian coordinates.
    *   No bonds/neighbor topology is needed.
    *   This can be used immediately as a substrate for on-surface placement.
*   **Surface Cleavage:**
    *   **Input:** Miller indices $(h, k, l)$.
    *   **Algorithm:** Create a slab geometry by rotating the basis so $(h, k, l)$ is the Z-normal, then replicate in X/Y.
*   **Step & Terrace Generator:**
    *   **Feature:** Create "vicinal surfaces" (high-index planes that look like steps).
    *   **Implementation:** instead of just cutting a flat plane, generate a "staircase" geometry.
    *   **Control:** User defines "Terrace width" (in unit cells) and "Step height" (usually one atomic layer or unit cell height).
    *   **Output:** A mesh representing the stepped surface with semantic data (identifying which atoms belong to the "upper terrace," "step edge," and "lower terrace").

### Module B: The Lattice Matcher (Geometric Solver)
Responsible for calculating how the polymer periodicity fits the crystal periodicity.

*   **Vector Selection:** User selects two points or a vector on the crystal surface (e.g., along the step edge). Let this be vector $\vec{V}_{sub}$.
*   **Commensurability Calculator:**
    *   **Input:** Substrate repeat length $|\vec{V}_{sub}|$ and Monomer backbone length $L_{mono}$.
    *   **Math:** Diophantine approximation. Find integers $N$ (substrate units) and $M$ (monomer units) such that:
        $$ N \times |\vec{V}_{sub}| \approx M \times L_{mono} $$
    *   **Tolerance:** Allow for a specific strain percentage (e.g., $\pm 3\%$).
*   **Visual Debugger:**
    *   Render a "ruler" along the step edge.
    *   Highlight tick marks where the lattice perfectly repeats.
    *   Ghost/Phantom rendering of the monomer backbone to visualize the fit before placing atoms.

### Module C: The Sequencer (Monomer Library)
Responsible for defining the heterogeneous composition of the polymer.

*   **Monomer Definition:**
    *   Each monomer is a pre-fab object (prefab).
    *   **Anchors:** Must define an "Entry Point" (Head), "Exit Point" (Tail), and "Up Vector" (for orientation relative to the surface).
*   **Sequence Editor:**
    *   A UI component (like a string editor) to define the chain: e.g., `A - B - A - C - A`.
    *   Support for arbitrary user-defined libraries (non-bio backbones).

### Module D: The Curve & Frame System
Responsible for the spatial path and orientation of the polymer.

*   **Spline Generation:**
    *   **Input:** Control points set by the user (snapped to crystal lattice points or free-hand).
    *   **Algorithm:** Catmull-Rom Splines (passes through points) or Cubic Bézier curves.
*   **Frame Generation (The "Up-Vector" Problem):**
    *   Standard Frenet-Serret frames are unstable (flip at inflection points).
    *   **Solution:** Use **Parallel Transport Frames** (Rotation Minimizing Frames) or explicitly interpolate user-defined Up-Vectors at each control point.
    *   **Output:** For every point $t$ on the curve, generate a $4\times4$ Transformation Matrix (Position + Rotation).

### Module E: The Assembler
Responsible for instantiating the molecules into the scene.

*   **Conformer Logic:**
    *   Iterate through the *Sequence*.
    *   Fetch the *Monomer Prefab*.
    *   Align Monomer [Head $\to$ Tail] vector with the Curve Tangent vector.
    *   Align Monomer [Up] vector with the Curve/Surface Normal vector.
    *   Apply local rotations (user inputs: "Twist" and "Tilt" angles).
*   **Batching:** Use **Hardware Instancing** (WebGL `InstancedMesh`) if the polymer contains hundreds of repeated units to maintain 60FPS.

---

## 3. Mathematical & Algorithmic Toolset

To implement this efficiently in JavaScript, you will need the following concepts:

### 1. Linear Algebra (glMatrix / Three.js Math)
*   **Basis Change:** Mapping coordinates from Crystallographic Basis (Fractional) to Cartesian Basis (World).
    $$ V_{cart} = [Matrix_{lattice}] \times V_{frac} $$
*   **Quaternions:** Essential for smooth interpolation of rotations along the curve (SLERP) to avoid Gimbal lock when aligning molecules.

### 2. Commensurability Algorithm
Pseudo-code for the Lattice Matcher:
```javascript
function findMatches(crystalPeriod, monomerPeriod, maxStrain) {
  let matches = [];
  for (let i = 1; i < 10; i++) { // Crystal multiples
    for (let j = 1; j < 10; j++) { // Monomer multiples
      let lenCrystal = i * crystalPeriod;
      let lenPolymer = j * monomerPeriod;
      let strain = Math.abs(lenPolymer - lenCrystal) / lenCrystal;
      
      if (strain < maxStrain) {
        matches.push({ N: i, M: j, strain: strain });
      }
    }
  }
  return matches.sort((a,b) => a.strain - b.strain);
}
```

### 3. Curve Frame Transport
To generate the orientation along the curve:
1.  Calculate Tangent $\vec{T}$ at point $P$.
2.  User provides an approximate Up vector $\vec{U}_{approx}$ (usually the surface normal).
3.  Calculate Binormal $\vec{B} = \vec{T} \times \vec{U}_{approx}$ (normalize it).
4.  Recalculate true Normal $\vec{N} = \vec{B} \times \vec{T}$.
5.  Construct Rotation Matrix from $[\vec{T}, \vec{N}, \vec{B}]$.

---

## 4. User Workflow (The "Story")

1.  **Step 1: The Stage.**
    *   User loads `NaCl.cif`.
    *   User opens "Slab Tool," selects Miller Index `(1 0 0)`.
    *   User enables "Step Edge Mode," requests a step every 5 unit cells.
    *   *System renders the terraced surface.*

2.  **Step 2: The Path.**
    *   User clicks "Polymer Path."
    *   User clicks two points along the bottom of a step edge (the "corner").
    *   *System draws a guide line.*
    *   User adds a third point curving away from the step onto the terrace.
    *   *System draws a spline.*

3.  **Step 3: The Match.**
    *   User selects a "Monomer A" from the library.
    *   System highlights that 3 units of Monomer A fit almost perfectly into 2 units of the Crystal Lattice along the straight edge.
    *   User confirms the scaling/stretching.

4.  **Step 4: The Sequence.**
    *   User defines sequence: `A-A-A-B-A-A-A`.
    *   User tweaks the "Tilt" angle because the sidechains are hitting the step wall.
    *   *System updates the visual instantly.*

5.  **Step 5: Export.**
    *   User exports `.xyz` or `.pdb` file containing both the fixed crystal atoms and the generated polymer coordinates.

---

## 4b. GUI Widgets / Panels (WebGL editor)

This section specifies the *concrete UI panels* needed to construct periodic substrates and polymer sequences directly in the web GUI. The goal is to keep the editor logic modular:

- **UI layer**: collects parameters, validates input, calls pure functions.
- **Data layer**: `MoleculeSystem.js` implements geometry/topology/transform operations.
- **Rendering layer**: redraws after `MoleculeSystem` mutations.

---

## 4c. Development status / changelog (what was done, what remains)

### What we achieved (Dec 2025)

- **GUI sidebar regression fixed**
  - **Symptom**: side panel was empty (no widgets).
  - **Cause**: `GUI.createSection()` was accidentally left empty (sections never appended).
  - **Fix**: implemented `createSection()` to append title + content, plus collapsible behavior.

- **GUI init crash fixed**
  - **Symptom**: no atoms rendered, console error `ReferenceError: mkInt is not defined`.
  - **Cause**: missing `mkInt` helper in the attachment-by-direction UI.
  - **Fix**: added local `mkInt` and populated the corresponding row.

- **Selection callback crash fixed**
  - **Symptom**: selecting atoms caused `TypeError: this.gui.updateSelectionCount is not a function`.
  - **Cause**: `main.js` calls `GUI.updateSelectionCount()` but the method did not exist.
  - **Fix**: implemented `GUI.updateSelectionCount()`.

- **Large system rendering bug fixed (critical)**
  - **Symptom**: slabs looked like partial “step edges” / missing layers, e.g. `nz=3` looked like 2 layers + fragmented extra layers.
  - **Cause**: `MoleculeSystem` auto-resized (e.g. capacity 1000 -> 4000) but the GPU renderer (`MeshRenderer`) stayed allocated for the old capacity, so only a prefix of atoms was correctly drawn.
  - **Fix**:
    - added `MeshRenderer.ensureCapacity(newCapacity)` which rebuilds internal textures/meshes
    - `MoleculeRenderer.update()` now calls `ensureCapacity(system.capacity)`
    - `MeshRenderer.updatePositions/updateParticles` now throw if `count > capacity` (fail loudly instead of silently corrupting rendering)

### Current known issues

- **Selection visualization**
  - current selection marker renders as large filled yellow “atoms” (efficient but occludes the real atoms).
  - desired: selection highlight should be translucent / non-occluding (alpha blended halo / outline).

- **CaF2 (fluorite) appears with Ca:F = 1:1**
  - could be a genuine basis/generation error, or could be overlapping positions / visual confusion.
  - next step is to debug by printing element counts and/or exporting generated crystals to `.xyz`.

- **CaCO3 (calcite) preset not implemented**
  - current code explicitly throws.
  - implementation needs a robust unit-cell storage/loading method for non-trivial bases (30 atoms/cell).

### Agreed next steps (implementation plan)

- **Add debugging outputs for crystal generation**
  - per-element counts (e.g. Ca vs F)
  - optional “export generated crystal as XYZ” (download as file) to verify geometry externally.

- **Fix CaF2 and implement CaCO3**
  - keep simple cubic presets hardcoded (rocksalt / fluorite)
  - for complex crystals, support loading unit cells from text (extended XYZ with lattice, or CIF).

- **Refactor for maintainability**
  - move specialized crystal/polymer creation helpers out of `MoleculeSystem.js` into:
    - `web/molgui_web/js/CrystalUtils.js`
    - `web/molgui_web/js/PolymerUtils.js`
  - keep `MoleculeSystem.js` focused on storage + generic operations (append, delete, transform, bonds).

- **Vector math unification**
  - currently `MoleculeSystem.js` contains ad-hoc vector helpers.
  - target: use `web/common_js/Vec3.js` consistently (larger change; do after correctness issues are fixed).

### Panel A: Substrate / Crystal Builder (foldable)

**A.1 Inputs (explicit lattice)**

- **Lattice vectors (3x3)**
  - input boxes for `a⃗, b⃗, c⃗` in Cartesian coordinates (Å)
  - layout: 3 rows (a,b,c) x 3 cols (x,y,z)
- **Replication counts**
  - integer inputs: `nx, ny, nz`
- **Lattice constant `a` (optional convenience)**
  - numeric input (Å)
  - if a preset uses `a` scaling, updating `a` updates `a⃗,b⃗,c⃗` (and/or basis)

**A.2 Presets (crystal type)**

- selection box with presets (initial list)
  - `NaCl` (rocksalt)
  - `KBr` (rocksalt)
  - `MgO` (rocksalt)
  - `CaF2` (fluorite)
  - `CaCO3` (calcite/aragonite; can start with one structure)

Each preset defines:

- `lvec` (or `a` + conventional cell vectors)
- `basis` (fractional or Cartesian)
- default charges (optional)
- default surface orientation and/or step-edge mode (optional)

**A.3 Surface orientation (Miller index slab)**

We need an option to orient the slab by a crystallographic plane (Miller indices) before replication.

- Miller index inputs: integer `(h, k, l)`
- optional: choice of whether the slab normal is `+n` or `-n`
- optional: thickness in unit cells along slab normal

**Intended behavior (high level):**

- For a given lattice (direct basis vectors), compute reciprocal lattice vectors.
- Compute surface normal `n⃗ ∥ h b⃗1 + k b⃗2 + l b⃗3`.
- Construct a rotation matrix that maps this normal to world `+z` (or user-chosen axis).
- Rotate lattice vectors and basis atom positions.
- Then replicate in the in-plane directions.

**Note:** Miller-index slab generation can be implemented incrementally:

1) v0: only allow axis-aligned surfaces via a dropdown (e.g. `100`, `110`, `111`).
2) v1: allow arbitrary `(h,k,l)` and compute orientation.

**A.4 Step-edge / terrace mode (optional, specific presets)**

- For `NaCl`-like surfaces, include a checkbox `step-edge mode` and inputs matching the existing generator:
  - `nx, ny, nz, a, Q0`
  - (later) step period and step direction

### Panel B: Attachment Tool (end-group / capping) (foldable)

This tool attaches a selected fragment (end-group) to a selected atom on the current system. This is distinct from marker-based attachment; it is user-driven picking.

**B.1 User selection**

- User loads/chooses an end-group fragment from the library.
- User **mouse-picks** an atom `A` on the backbone to attach to (typically a hydrogen to replace).
- Optional: user also picks a second atom `B` bonded to `A` (if not, we infer it from bond graph / nearest heavy atom).

**B.2 Geometry definition**

- **Forward direction** `f⃗`
  - default: `f⃗ = normalize( A - B )` (direction outward from backbone)
  - user option: invert direction
- **Bond length**
  - numeric input `bondLen` (Å)
  - modes:
    - `auto`: use covalent radii (approx)
    - `manual`: user overrides `bondLen`
- **Up-vector** `u⃗0` in world coords
  - 3 numeric inputs (default `0,0,1`)
  - will be orthogonalized against `f⃗`
- **Twist angle** `phi` (deg)
  - rotate around `f⃗` after orthogonalization

**B.3 Transform math (exactly defined)**

Given `f⃗` and `u⃗0`:

1) `f⃗ = normalize(f⃗)`
2) `u⃗ = normalize( u⃗0 - f⃗ * dot(u⃗0, f⃗) )`
3) `l⃗ = normalize( cross(u⃗, f⃗) )`
4) build rotation frame `R = [f⃗, u⃗, l⃗]` (column or row convention must match code)
5) apply twist rotation by `phi` in the `{u⃗,l⃗}` plane around `f⃗`

**B.4 Topology update**

- remove the selected capping atom (e.g. H) from the backbone
- attach end-group anchor atom to backbone atom `B` (or to `A` depending on definition)
- add bond of order 1 (initially)

**Note:** this “pick-based attach” is *different* from marker-based attach. Both should exist:

- pick-based: for interactive editing
- marker-based: for robust batch processing and library-driven assembly

### Panel C: Polymer Sequence Builder (foldable)

**C.1 Sequence language (human input)**

We want a compact language:

- Each monomer token starts with an uppercase letter and can contain additional letters (e.g. `D`, `Gly`, `A`, `PNA`)
- An optional integer repeat count after the token repeats it

Examples:

- `D3Gly6A` means `DDD + GlyGlyGlyGlyGlyGly + A`
- `PNA10` means 10 repeats of `PNA`

**Parsing rule (simple, deterministic):**

- token = `[A-Z][a-zA-Z]*`
- count = `[0-9]+` (optional; default=1)

**C.2 Monomer definitions (library)**

We need a monomer library that can be loaded in the GUI and used both for:

- sequence building
- end-group attachment

The recommended approach is a **JSON library file** plus separate geometry files.

#### Proposed JSON schema (v0)

```json
{
  "version": 0,
  "name": "Example monomer library",
  "basePath": "./monomers/",
  "monomers": [
    {
      "id": "D",
      "aliases": ["DANA"],
      "file": "DANA_CG_2inv.mol2",
      "format": "mol2",
      "repeatVec": "lvec1",
      "anchors": {
        "head": { "type": "index", "value": 1 },
        "tail": { "type": "index", "value": 4 }
      },
      "markers": {
        "X": "Se",
        "Y": "Cl"
      }
    }
  ],
  "endGroups": [
    {
      "id": "Guanine",
      "file": "guanine-SeCl.mol2",
      "format": "mol2",
      "markers": { "X": "Se", "Y": "Cl" }
    }
  ]
}
```

**Interpretation:**

- `repeatVec`:
  - `"lvec1"` means translate by the monomer file lattice vector `lvec[1]` (as in `tests/tAttach/polymerize.py`)
  - later we can allow explicit vector: `{ "type":"vec3", "value":[x,y,z] }`
- `anchors.head/tail`:
  - `index` is 1-based atom index in the geometry file
  - alternative selector: by element + occurrence

#### Anchor selector alternative: by symbol occurrence

For convenience in human-made files:

- `{"type":"symbolOrdinal", "symbol":"Cl", "ordinal":1}` selects the first chlorine atom (counting only Cl)
- `{"type":"atomName", "value":"Cl1"}` if input format provides unique atom names (mol2 does)

The editor should support at least `index` first (v0) because it is unambiguous.

**C.3 Build action**

- User selects/loads library JSON.
- User types sequence string.
- User presses `Build`.
- System generates polymer as a new `MoleculeSystem` fragment (or adds into current scene) using `MoleculeSystem.assemblePolymerSequence(...)`.

**C.4 Heterogeneous sequences**

For heterogeneous sequences, the join rule remains:

- connect previous monomer `tail` to new monomer `head`

but placement needs a well-defined repeat translation. Options:

- simplest: use each monomer’s `repeatVec` (e.g. its `lvec[1]`)
- later: allow explicit path-based placement (curve tool) where translation is derived from curve frames

---

## 5. Technical Stack Recommendations

*   **Core Library:** **Three.js** (Standard, huge community) or **Babylon.js** (Great support).
*   **Math:** **Math.js** (for heavy matrix/algebra if Three.js isn't enough) or **numeric.js**.
*   **Crystal Math:** **Crystallography.js** (doesn't exist as a major lib, likely needs custom implementation using standard transformation formulas).
*   **UI Overlay:** **Dat.GUI** or **Leva** (React-based) for tweaking angles/integers.
*   **Performance:** Use **Web Workers** for the "Lattice Matcher" calculation so the UI doesn't freeze while checking integer multiples.

---

# Implementation inspirations / references in FireCore (do not re-invent)

This project already contains working C++ and Python code for lattice matching, replication, and assembly/attachment operations. The WebGL editor should reuse these ideas.

## A) Lattice / surface replication (pure geometry)

- **Reference: periodic replication**
  - `pyBall/AtomicSystem.py`
    - `AtomicSystem.clonePBC(nPBC=(1,1,1))`  # replicates atoms by lattice vectors, no chemistry assumptions
    - note: it shifts by `shift = lvec[0]*ix + lvec[1]*iy + lvec[2]*iz`

## B) Lattice matching / commensurability (2D supercell solver)

- **Core algorithm (C++)**
  - `cpp/common/molecular/LatticeMatch2D.h`
    - `LatticeMatch2D.walk2D(Rmax, dmax)`  # enumerates integer combinations on lat0 and matches lengths vs lat1
    - `LatticeMatch2D.matchAngles(dangMax)`  # pairs candidates by angle compatibility
    - `LatticeMatch2D.exportVecMatch(...)` / `LatticeMatch2D.exportMatch(...)`  # ranking by distortion costs

- **C interface wrapper**
  - `cpp/libs/Molecular/Lattice2D_lib.cpp`
    - `walk2D(double* lat0, double* lat1, ...)`
    - `match(double* lat0, double* lat1, ...)`
    - `getVecMatch(...)`, `getMatches(...)`

- **Python wrapper + usage example**
  - `pyBall/Lattice2D.py`
    - `walk2D(lat0, lat1, Rmax, dRmax)`
    - `match(lat0, lat1, Rmax, dRmax, dAngMax)`
    - `getVecMatch(...)`, `getMatches(...)`
  - `tests/tLattice2D/run.py`
    - example parameters for "Polymer @ NaCl" and retrieval of best matches

## C) Polymer sequencing and fragment assembly (topology + transforms)

- **Sequence building by repeating a monomer cell (translation + bond join)**
  - `tests/tAttach/polymerize.py`
    - uses a `monomers` dictionary mapping letters to `(file, (head,tail))`
    - advances placement by `pos += B.lvec[1]` (monomer repeat vector)
    - joins fragments by adding a bond when merging: `A.addSystems(..., added_bonds=[(...)])`

- **Generic merging (no atom deletion, just append + optional bond connect)**
  - `pyBall/AtomicSystem.py`
    - `AtomicSystem.addSystems(other, pos=None, rot=None, added_bonds=None, _0=1)`
      - applies optional rotation/translation
      - appends arrays and remaps bond indices by `offset`

- **Attach end-group by explicit anchor vectors (forward/up) and delete leaving atom**
  - `pyBall/AtomicSystem.py`
    - `AtomicSystem.attach_group(...)`
      - uses forward vector from a specified bond and an up-vector to build an orientation frame

- **Attach by marker atoms (robust anchor detection, less manual index bookkeeping)**
  - `pyBall/AtomicSystem.py`
    - `AtomicSystem.attach_group_by_marker(G, markerX="Xe", markerY="He", _0=1)`
    - uses `find_marker_pairs(markerX, markerY)` which returns `(iX, iY, iAnchor)` triplets
  - `tests/tAttach/attach_new3.py`  # example: attach guanine endgroup via markers `Se/Cl`

- **Rotation matrix from points / forward+up frame construction**
  - `pyBall/atomicUtils.py`
    - `rotmat_from_points(ps, ifw, iup, fw, up, _0=1)`
    - `makeRotMat(fw, up)`
    - `makeVectros(...)` (build forward/up from chosen indices)

## D) WebGL editor codebase (current state)

- **Web app entry**
  - `web/molgui_web/js/main.js`
- **Data model (atoms + bonds + selection)**
  - `web/molgui_web/js/MoleculeSystem.js`
- **Picking + transform gizmo (selection transform = Subtask 2 primitive)**
  - `web/molgui_web/js/Editor.js`
- **Rendering (instanced atoms/bonds/selection overlays)**
  - `web/molgui_web/js/MoleculeRenderer.js`
  - `web/common_js/MeshRenderer.js`

