
## USER 1

I wan to make simple crystal / surface genrator in javascript, something a bit like in ASE (atomistic simulation environment)
It should be based on lattice vectors (or angles), few atoms in elementary unit cell, and number of copies. we have then number of replicas nx,ny,nz, and we can choose some miller indexes for the surface.
Notice that for NaCl which if body-centred of face centered cubic, it is possible to find crystalography direction where Na-Cl ions alternate in checkerboead pattern on the surface, which is the most stable. This is not what ASE on my javascript generate, but I generated it previously with python.


I dont wat webpage or final javascript code from you because it does not reuse modules I already have it does not integrate to the larger system. The crystal unit-cells should be in plain text or like .xyz format.
I want from you clear description of the crystalogrpahy in general and of the systems, so that the LLM can implement it accordingly. You are expert on crystalography the other LLM is expert of java/web-dev
https://en.wikipedia.org/wiki/Miller_index
https://next-gen.materialsproject.org/materials/mp-3953
https://en.wikipedia.org/wiki/Calcium_carbonate
https://en.wikipedia.org/wiki/Calcium_fluoride
https://next-gen.materialsproject.org/materials/mp-2741


---

## Gemini

Here is the technical specification you should feed to the "Coder LLM."

It separates the **Crystallography Data** (the definitions) from the **Transformation Logic** (how to cut/rotate/step).

---

## Prompt for the Coder LLM

**Context:** We are building a crystal generator. You are responsible for the JavaScript implementation.
**Critical Constraint:** Do NOT use "Primitive Unit Cells" (which are often rhombohedral or slanted). You must use **"Conventional Unit Cells"** (which are usually Cubic or Hexagonal boxes). This is required to ensure that standard surfaces (like the NaCl 100 face) appear as stable, neutral checkerboards rather than charged polar slices.

### Part 1: Unit Cell Definitions (The Database)

Store these exactly as defined. These are the **Conventional Cells**.

#### 1. NaCl (Rock Salt Structure)
*Use this to guarantee the checkerboard surface.*
*   **System:** Cubic
*   **Lattice Constant:** $a = 5.64$ Å (approx)
*   **Lattice Vectors:**
    *   $\vec{a} = [5.64, 0, 0]$
    *   $\vec{b} = [0, 5.64, 0]$
    *   $\vec{c} = [0, 0, 5.64]$
*   **Basis (8 atoms - Fractional Coordinates):**
    *   Na: `(0,0,0)`, `(0.5, 0.5, 0)`, `(0.5, 0, 0.5)`, `(0, 0.5, 0.5)`
    *   Cl: `(0.5, 0, 0)`, `(0, 0.5, 0)`, `(0, 0, 0.5)`, `(0.5, 0.5, 0.5)`
*   **Note:** Na is +1, Cl is -1.

#### 2. CaF2 (Fluorite)
*   **System:** Cubic
*   **Lattice Constant:** $a = 5.46$ Å
*   **Lattice Vectors:** orthogonal, length 5.46.
*   **Basis (12 atoms - Fractional Coordinates):**
    *   Ca (FCC positions): `(0,0,0)`, `(0.5, 0.5, 0)`, `(0.5, 0, 0.5)`, `(0, 0.5, 0.5)`
    *   F (Tetrahedral sites):
        *   `(0.25, 0.25, 0.25)`, `(0.75, 0.75, 0.75)`
        *   `(0.75, 0.25, 0.25)`, `(0.25, 0.75, 0.75)`
        *   `(0.25, 0.75, 0.25)`, `(0.75, 0.25, 0.75)`
        *   `(0.25, 0.25, 0.75)`, `(0.75, 0.75, 0.25)`

#### 3. CaCO3 (Calcite)
*Calcite is Trigonal but best handled in the Hexagonal setting for surface generation.*
*   **System:** Hexagonal
*   **Parameters:** $a = 4.99$ Å, $c = 17.06$ Å
*   **Lattice Vectors (Cartesian approximation):**
    *   $\vec{a} = [4.99, 0, 0]$
    *   $\vec{b} = [-2.495, 4.321, 0]$  *(This is $-0.5a, \frac{\sqrt{3}}{2}a, 0$)*
    *   $\vec{c} = [0, 0, 17.06]$
*   **Basis (Simplified for Hexagonal Cell - 30 atoms):**
    *   *Note to Coder:* For CaCO3, allow the user to load a `.xyz` or `.cif` string for the unit cell, as the basis list is too long to hardcode manually.

---

### Part 2: Surface Generation Logic (Miller Indices)

Do not just replicate the unit cell. You must orient the crystal so the user's desired face points UP (Z-axis).

**Algorithm:**
1.  Read inputs: Target Miller indices $(h, k, l)$.
2.  If $(h,k,l) = (0,0,1)$ (and orthogonal cell), no rotation needed.
3.  Otherwise, calculate a **Rotation Matrix** $R$ such that the vector normal to the $(hkl)$ plane becomes the new Z-axis.
    *   *Step A:* Calculate normal vector $\vec{n} = h\vec{b}\times\vec{c} + k\vec{c}\times\vec{a} + l\vec{a}\times\vec{b}$ (Reciprocal lattice vector).
    *   *Step B:* Normalize $\vec{n}$.
    *   *Step C:* Define new Z basis vector $\vec{z'} = \vec{n}$.
    *   *Step D:* Construct arbitrary orthogonal vectors $\vec{x'}$ and $\vec{y'}$ perpendicular to $\vec{z'}$.
4.  Apply this rotation $R$ to all atoms and lattice vectors *before* replication.

---

### Part 3: Vicinal Steps (The "Python Script" Logic)

The user requires a specific "Vicinal Surface" mode that mimics a Python reference script. This creates atomic steps (terraces) by shearing the crystal.

**Implementation Logic:**
1.  First, generate the regular supercell grid (positions $x, y, z$).
2.  Apply the following coordinate transformation to every atom:
    *   Calculate shear factor: `slope = -1.0 / nx` (where `nx` is the number of replicas in X).
    *   `z_new = z + slope * x`
    *   `x_new = x - slope * z` (Small angle approximation rotation/shear).
3.  **Periodic Wrap (The Step Creation):**
    *   Define the step boundary `Lx_step = TotalLengthX * 0.5`.
    *   If `x_new > Lx_step`:
        *   Shift atom up: `z_new += UnitCell_Z_Height`.
        *   (This creates the physical "step up" in the middle of the terrace).

---

### Summary of Inputs/Outputs for the Function

**Function Signature:**
`generateCrystal(latticeVectors, basisAtoms, replicas[nx,ny,nz], miller[h,k,l], applyVicinalShear)`

**Process Flow:**
1.  **Load Conventional Cell** (Vectors + Basis).
2.  **Rotate** Cell based on Miller Indices so $(hkl)$ is up.
3.  **Replicate** atoms $n_x, n_y, n_z$ times.
4.  If `applyVicinalShear == true`: Apply the Z/X shear and the step-jump logic.
5.  **Output:** List of atoms (Symbol, Cartesian X, Y, Z).


----

## USER 2

OK, I think we should generalize this software a bit more so we can create 3D crystals. One topic I also study are silicon and diamond nanocrystals which have basically the same structure. We should be able to generate them in 3D.

But I was thinking that while generating it we can apply some cut-planes which can create some faces, which we can define by miller indexes like 100, 110, 111 and all symmetric permutations. This we can apply as post-processing operator (first generating the atoms and them deleting some, or apply it during generation. perhaps we should try both option, as it will result in different indexing.

We should also thing about making bonds between the atoms and pasivation (e.g. hydrogen capping). Perhaps searching bonds by bond-length in post-processing is fine (althout it is O(n^2) operation unless we use some spacial indexing datastructure). But still it would be interesting to make the bonds explicitly - it should not be so hard if we know the indexing in 3D crystal (perhaps it would be simpler if we aply cutting plane as post procesing operator because than the atom index correspond to its position on 3D grid.

Any other ideas what would be good to implement? I was thinking about some defects and surface reconstruction. But that is more thing for bonding topology, so we will implement this in different module (not the crystalogrpahy module we are discussing now)

please think if through, and write it systematically as a desing doc where you fill in more technical details how to solve it. But pleas do not start to write the code yet.

---

This is a significant architectural expansion. Moving from 2D slabs to 3D nanocrystals with specific faceting (Wulff construction) and topology (bonding/passivation) requires a more robust data structure than a simple list of atoms.

Here is the **Design Document** for the **3D Crystal & Nanocrystal Generator Module**.

---

# Design Document: 3D Nanocrystal & Surface Generator

## 1. System Overview
The system is divided into three distinct pipeline stages:
1.  **Lattice Definition:** Static definitions of unit cells (Si, Diamond, etc.).
2.  **Volumetric Generation & Shaping:** Generating a "supercell" and applying geometric cuts (Miller planes) to define the nanocrystal shape.
3.  **Topology & Passivation:** Post-processing to determine connectivity (bonds) and cap dangling bonds chemically (e.g., Hydrogen).

## 2. Core Crystallography (The "Database")
To support Silicon and Diamond, we must define the **Diamond Cubic** structure strictly.

*   **Structure:** Face-Centered Cubic (FCC) with a 2-atom basis.
*   **Conventional Cell:** Cube with side $a$.
*   **Lattice Vectors:** $\vec{v}_1=[a,0,0], \vec{v}_2=[0,a,0], \vec{v}_3=[0,0,a]$.
*   **Basis Atoms (Fractional):**
    1.  $(0, 0, 0)$
    2.  $(0.25, 0.25, 0.25)$
    3.  *(Plus FCC translations: face centers)*
*   **Parameters:**
    *   **Silicon (Si):** $a \approx 5.431$ Å. Bond length $d \approx 2.35$ Å.
    *   **Diamond (C):** $a \approx 3.567$ Å. Bond length $d \approx 1.54$ Å.

---

## 3. Volumetric Generation & Shaping Strategy

We need to support defining shapes via **Intersecting Half-Spaces** (Planes).

### 3.1. The Shaping Algorithm (Wulff Construction Logic)
Instead of growing the crystal, we will use a **"Carving" strategy** (Subtractive Manufacturing). This preserves the integer lattice indexing ($i,j,k$) which is crucial for $O(N)$ bonding later.

1.  **Define Bounding Box:** Calculate the maximum extent of the crystal based on the requested planes and center point.
2.  **Generate Grid:** Create a dense grid of atoms filling this box (a standard rectangular supercell).
3.  **Apply Cut Planes:**
    *   A plane is defined by Miller index $(hkl)$ and a distance $D$ from the center.
    *   Plane Normal $\vec{n} = (h, k, l)$ (normalized).
    *   Equation: An atom at position $\vec{r}$ is **kept** if $\vec{r} \cdot \vec{n} \le D$.
    *   Iterate through all generated atoms. If an atom violates *any* active cut plane, mark it as `deleted`.

### 3.2. Symmetry Handling (Permutations)
Users usually say "I want {100} facets." The system must expand this into all symmetric equivalents for the cubic system:
*   **{100} Family (6 planes):** $(1,0,0), (-1,0,0), (0,1,0)...$
*   **{110} Family (12 planes):** $(1,1,0), (-1,1,0), (1,0,1)...$
*   **{111} Family (8 planes):** $(1,1,1), (-1,-1,-1)...$

**Design Choice:** The UI should allow adding a "Family" of planes (e.g., `{111}`) with a single slider for distance $D_{111}$.

---

## 4. Topology & Bonding (The Connectivity Module)

You correctly identified that $O(N^2)$ distance searching is inefficient. Since we generate the atoms on a rigid lattice, we can use **Lattice Lookahead** (O(1) per atom, O(N) total).

### 4.1. Explicit Bonding Strategy
In a Diamond Cubic lattice, every atom has exactly 4 neighbors at specific relative vectors.

1.  **Lookup Table:** Define the 4 vectors $\vec{\delta}_1, \vec{\delta}_2, \vec{\delta}_3, \vec{\delta}_4$ connecting a basis atom to its neighbors.
2.  **Grid Indexing:**
    *   Store atoms in a 3D array (or 1D array mapped to $i,j,k$).
    *   For an atom at index $(i,j,k)$, calculating the index of a potential neighbor is integer math (e.g., neighbor is at $i, j, k+1$ or $i+1, j+1, k$, etc.).
3.  **Verification:**
    *   Calculate expected neighbor index.
    *   Check if that neighbor exists AND is not marked `deleted` by the shaping step.
    *   If yes -> Create Bond.
    *   If no -> Increment `dangling_bond_count` and store the `missing_vector`.

### 4.2. Surface Reconstruction & Passivation (Hydrogen Capping)
Once the bulk is carved, we handle the surface.

*   **Passivation Logic:**
    *   Iterate all atoms.
    *   If `bond_count < 4` (for Group IV elements):
        *   We have `4 - bond_count` dangling bonds.
        *   The direction of the dangling bond is exactly the `missing_vector` from the lattice lookup.
    *   **Placement:** Place a Hydrogen atom at $\vec{r}_{atom} + \text{normalize}(\vec{missing\_vector}) \times d_{X-H}$.
    *   *Note:* $d_{Si-H} \approx 1.48$ Å, $d_{C-H} \approx 1.09$ Å.

*   **Steric Hindrance check (Optional but recommended):**
    *   On complex facets (high index), placing H atoms on every dangling bond might cause collisions.
    *   Simple check: If generated H is too close to another H ($< 0.7$ Å), delete one or bridge them (advanced).

---

## 5. Data Structures (Technical Specs)

To make this efficient in JavaScript (and translatable to C++/WASM later if needed), we should avoid arrays of objects `[{x,y,z}, ...]` for the core processing.

**Proposed Structure: Structure-of-Arrays (SoA)**

```javascript
class CrystalSystem {
    // 1. Configuration
    float latticeConstant;
    
    // 2. The Grid (Raw Data)
    // Flattened 1D arrays representing the 3D grid
    Float32Array positions_x;
    Float32Array positions_y;
    Float32Array positions_z;
    Uint8Array   status;      // 0=Deleted, 1=Bulk, 2=Surface, 3=PassivationAgent
    
    // 3. Topology
    Int32Array   neighbor_indices; // Stride 4 (4 neighbors per atom). -1 if missing.
    
    // 4. Output Lists
    // Generated only when exporting to .XYZ or rendering
    List<Atom> finalAtomList;
    List<Bond> finalBondList;
}
```

## 6. Workflow Implementation Plan

1.  **Step 1: Grid Generator:**
    *   Input: `Approximate Radius`.
    *   Action: Create integer grid $N \times N \times N$ sufficient to hold the radius. Populate with basis atoms.

2.  **Step 2: Shape Carver:**
    *   Input: List of `{h, k, l, distance}`.
    *   Action: Loop through atoms. Compute $dist = \vec{r} \cdot \vec{n}_{(hkl)}$. If $dist > distance$, set `status[i] = 0`.

3.  **Step 3: Topology Builder:**
    *   Action: Loop through atoms where `status[i] == 1`.
    *   Check theoretical neighbors using integer offsets.
    *   If neighbor is `status[0]`, record as missing bond vector.
    *   If neighbor is `status[1]`, record bond.

4.  **Step 4: Passivator:**
    *   Action: For every missing bond vector on a `status[1]` atom, create a new atom (Hydrogen) at the correct distance.
    *   Add this new atom to the atom list (with `status[3]`).

5.  **Step 5: Output:**
    *   Format as `.XYZ` or JSON for the visualization engine.

## 7. Future Extensibility (Defects/Reconstruction)

*   **Defects (Vacancies):** Easy. Randomly select `status[1]` atoms and flip them to `status[0]`. Re-run Topology Builder to update bonds.
*   **Surface Reconstruction ($2 \times 1$ dimer):**
    *   This is hard to do geometrically. It requires finding pairs of surface atoms with dangling bonds and moving them closer.
    *   This should be a separate "Relaxation" step, possibly using a simple spring force model (simple molecular mechanics) to let the surface atoms snap into reconstructed positions.

---

### Comparison of User's Options

| Strategy | Pros | Cons | Recommendation |
| :--- | :--- | :--- | :--- |
| **Post-Process (Carving)** | Indexing is trivial ($i,j,k$). Neighbor finding is O(1). Easy to implement on GPU/TypedArrays. | Wastes memory on "deleted" atoms (empty air). | **Recommended** for nanocrystals < 1M atoms. |
| **Generative (Growing)** | Memory efficient. Only stores real atoms. | Neighbor finding requires spatial hash/k-d tree (slower). Harder to handle complex cuts. | Best only if simulating massive billion-atom systems. |