
## User

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