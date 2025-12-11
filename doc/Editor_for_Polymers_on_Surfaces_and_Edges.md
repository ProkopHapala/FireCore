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

---

## 2. Functional Modules

### Module A: The Crystallography Engine
Responsible for generating the "stage" upon which chemistry happens.

*   **Unit Cell Import:** Parsers for `.cif` or `.xyz` files to define lattice vectors ($\vec{a}, \vec{b}, \vec{c}$) and basis atoms.
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

## 5. Technical Stack Recommendations

*   **Core Library:** **Three.js** (Standard, huge community) or **Babylon.js** (Great support).
*   **Math:** **Math.js** (for heavy matrix/algebra if Three.js isn't enough) or **numeric.js**.
*   **Crystal Math:** **Crystallography.js** (doesn't exist as a major lib, likely needs custom implementation using standard transformation formulas).
*   **UI Overlay:** **Dat.GUI** or **Leva** (React-based) for tweaking angles/integers.
*   **Performance:** Use **Web Workers** for the "Lattice Matcher" calculation so the UI doesn't freeze while checking integer multiples.
