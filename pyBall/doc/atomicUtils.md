# AtomicSystem Class Documentation

The `AtomicSystem` class provides a framework for representing and manipulating molecular systems. It stores atomic positions, types, element names, charges, bonds, lattice vectors, and auxiliary labels. The class also offers a suite of functions for geometry manipulation, bond finding, rotation, cloning with periodic boundaries, and attachment of end–groups.

> **Note:** Most functions assume that per–atom properties (such as `qs`, `Rs`, and `aux_labels`) are pre‐initialized. See the separate pre–initialization documentation for details.

---

## Contents

- [Initialization and I/O](#initialization-and-io)
- [Geometry and Orientation Methods](#geometry-and-orientation-methods)
- [Bond and Neighborhood Functions](#bond-and-neighborhood-functions)
- [System Manipulation and Cloning](#system-manipulation-and-cloning)
- [Attachment of End–Groups](#attachment-of-end-groups)
- [Additional Utility Methods](#additional-utility-methods)
- [Examples](#examples)

---

## Initialization and I/O

### `__init__(fname=None, apos=None, atypes=None, enames=None, lvec=None, qs=None, Rs=None, bonds=None, ngs=None, bReadN=True)`
- **Description:**  
  Initializes an `AtomicSystem` instance. The system can be built from a file (for example, an `.xyz` or `.mol` file) or directly from provided arrays.
- **Parameters:**
  - `fname` (str): Filename to load (determines file format by extension).
  - `apos` (ndarray): Array of atomic positions (N×3).
  - `atypes` (array): Array of atomic numbers.
  - `enames` (list/array): List of element symbols.
  - `lvec` (ndarray): Lattice vectors (3×3 array).
  - `qs` (array): Atomic charges.
  - `Rs` (array): Atomic radii.
  - `bonds` (array): Bonds (each bond is a tuple of two indices).
  - `ngs` (any): Neighbors data (if available).
  - `bReadN` (bool): Whether to use the first line as number-of-atoms.
- **Example:**
  ```python
  from pyBall.atomicUtils import AtomicSystem
  system = AtomicSystem(fname='molecule.xyz')
  system.print()
  ```

### `saveXYZ(fname, mode="w", blvec=True, comment="", ignore_es=None, bQs=True, other_lines=None)`
- **Description:**  
  Saves the current system in XYZ format. If lattice vectors (`lvec`) exist and `blvec` is True, they are prepended to the comment line.
- **Parameters:**  
  See signature.
- **Example:**
  ```python
  system.saveXYZ("output.xyz", comment="My molecule")
  ```

### `toLines()`
- **Description:**  
  Returns the system as a list of strings (one per atom) in a human–readable format.
- **Example:**
  ```python
  lines = system.toLines()
  for line in lines:
      print(line)
  ```

### `toXYZ(fout, comment="#comment", ignore_es=None, other_lines=None, bHeader=False)`
- **Description:**  
  Writes the system to an open file-like object in XYZ format.
- **Example:**
  ```python
  with open("output.xyz", "w") as f:
      system.toXYZ(f, comment="XYZ output")
  ```

---

## Geometry and Orientation Methods

### `print()`
- **Description:**  
  Prints a summary of the system to the console including index, atom type, element, and position. If `aux_labels` are available, they are also printed.
- **Example:**
  ```python
  system.print()
  ```

### `getValenceElectrons()`
- **Description:**  
  Returns an array with the number of valence electrons for each atom. This is determined by looking up the element in the elements dictionary.
- **Example:**
  ```python
  valence = system.getValenceElectrons()
  print(valence)
  ```

### `subtractValenceE(f0=-1.0, f=+1.0)`
- **Description:**  
  Adjusts the system’s atomic charges (`qs`) by subtracting a fraction of the valence electrons. This function multiplies the current charges by `f0` and then adds `f` times the valence electron count.
- **Example:**
  ```python
  system.subtractValenceE(f0=-1.0, f=+1.0)
  ```

### `findBonds(Rcut=3.0, RvdwCut=1.5, RvdWs=None, byRvdW=True)`
- **Description:**  
  Finds bonds in the system using inter–atomic distances and the sum of van der Waals radii (if available).
- **Returns:**  
  A tuple `(bonds, rs)` where `bonds` is a list of atom index pairs and `rs` contains bond lengths.
- **Example:**
  ```python
  bonds, rs = system.findBonds()
  system.printBonds()
  ```

### `findHBonds(Rb=1.5, Rh=2.5, angMax=60.0, typs1={"H"}, typs2={"O","N"}, bPrint=False, bHbase=False)`
- **Description:**  
  Finds hydrogen bonds in the system based on distance and angle criteria.
- **Example:**
  ```python
  hbonds, hb_rs = system.findHBonds()
  print("Found hydrogen bonds:", hbonds)
  ```

### `neighs(bBond=True)`
- **Description:**  
  Computes the neighbor list for each atom based on the bonds.
- **Returns:**  
  A list (or dictionary) of neighbor information.
- **Example:**
  ```python
  neighbors = system.neighs()
  print(neighbors)
  ```

### `find_groups()`
- **Description:**  
  Identifies groups (clusters) of atoms that are bonded together. Typically used to identify functional groups.
- **Example:**
  ```python
  groups = system.find_groups()
  print(groups)
  ```

### `select_by_ename(elist)`
- **Description:**  
  Returns a list of indices for atoms whose element name is in `elist`.
- **Example:**
  ```python
  carbons = system.select_by_ename(["C"])
  print("Carbon indices:", carbons)
  ```

### `getNeighsOfType(selection, typ='N')`
- **Description:**  
  For each atom in `selection`, returns the neighbors of a specific element type (default: nitrogen).
- **Example:**
  ```python
  n_neigh = system.getNeighsOfType([0,1,2], typ='N')
  print(n_neigh)
  ```

### `select_by_neighType(neighs, typ='N', neighTyps={'H':(1,2)})`
- **Description:**  
  Selects atoms based on the type and number of neighboring atoms.
- **Example:**
  ```python
  selected = system.select_by_neighType(system.neighs(), typ='N', neighTyps={'H': (1,2)})
  print(selected)
  ```

### `findAngles(select=None, ngs=None)`
- **Description:**  
  Computes all bond angles (in radians) for atoms (optionally, a subset given by `select`) using neighbor information.
- **Example:**
  ```python
  angles, angle_indices = system.findAngles()
  print("Angles (radians):", angles)
  ```

### `findDihedral(select=None, ngs=None, neighTyp={'H'})`
- **Description:**  
  Computes dihedral angles in the system.
- **Example:**
  ```python
  dihedrals, indices = system.findDihedral()
  print("Dihedral angles:", dihedrals)
  ```

### `findCOG(apos, byBox=False)`
- **Description:**  
  Returns the center-of-geometry (COG) for the provided positions. If `byBox` is True, it computes the COG based on the bounding box.
- **Example:**
  ```python
  center = system.findCOG(system.apos)
  print("Center of geometry:", center)
  ```

### `projectAlongBondDir(i0, i1)`
- **Description:**  
  Projects all atom positions along the bond direction defined by atoms `i0` and `i1`.
- **Example:**
  ```python
  projections = system.projectAlongBondDir(1, 2)
  print(projections)
  ```

---

## System Manipulation and Cloning

### `store_bond_lengths()`
- **Description:**  
  Computes and stores the lengths of all bonds in the system. The results are stored in a dictionary (`bond_lengths`) as keys with a tuple of atom indices.
- **Example:**
  ```python
  bond_lengths = system.store_bond_lengths()
  print("Bond lengths:", bond_lengths)
  ```

### `restore_bond_length(ij, L=None)`
- **Description:**  
  Adjusts the position of atoms in a given bond (specified by tuple `ij`) so that the bond length is restored to a given value or to its original length.
- **Example:**
  ```python
  system.restore_bond_length((0, 1))
  ```

### `clonePBC(nPBC=(1,1,1))`
- **Description:**  
  Creates a periodic boundary condition (PBC) clone of the system. The system is replicated according to `nPBC` (number of copies in x, y, and z directions).
- **Returns:**  
  A new `AtomicSystem` instance with replicated positions and updated lattice vectors.
- **Example:**
  ```python
  system_clone = system.clonePBC(nPBC=(2,2,1))
  system_clone.print()
  ```

### `symmetrized(d=0.1)`
- **Description:**  
  Returns a symmetrized version of the system by replicating atoms near the cell boundaries. Also returns weighting factors.
- **Example:**
  ```python
  sym_sys, ws = system.symmetrized(d=0.1)
  sym_sys.print()
  ```

### `selectSubset(inds)`
- **Description:**  
  Returns a new `AtomicSystem` that is a subset of the current system containing only the atoms whose indices are in `inds`.
- **Example:**
  ```python
  subset = system.selectSubset([0, 2, 4, 6])
  subset.print()
  ```

### `selectBondedCluster(s)`
- **Description:**  
  Given a starting set `s` (as a set of atom indices), returns two lists: one of indices in the cluster and one of the remaining indices.
- **Example:**
  ```python
  ins, outs = system.selectBondedCluster({0})
  print("Cluster:", ins, "Outside:", outs)
  ```

---

## Attachment of End–Groups

### `makeRotMat(ip1, ip2, _0=1)`
- **Description:**  
  Computes a rotation matrix from the forward vector (from atoms in `ip1`) and the up vector (from atoms in `ip2`).
- **Example:**
  ```python
  rotmat = system.makeRotMat((1,2), (3,4))
  print("Rotation matrix:\n", rotmat)
  ```

### `orient_mat(rot, p0=None, bCopy=False)`
- **Description:**  
  Applies a rotation matrix (`rot`) to the atomic positions. Optionally subtracts `p0` before rotation.
- **Example:**
  ```python
  new_positions = system.orient_mat(rotmat, p0=system.apos[0], bCopy=True)
  print(new_positions)
  ```

### `orient_vs(fw, up, p0=None, trans=None, bCopy=False)`
- **Description:**  
  Computes a rotation matrix from the forward (`fw`) and up (`up`) vectors, and applies it to the positions.
- **Example:**
  ```python
  new_positions = system.orient_vs(fw=[1,0,0], up=[0,0,1], p0=system.apos[0])
  print(new_positions)
  ```

### `orient(i0, b1, b2, _0=1, trans=None, bCopy=False)`
- **Description:**  
  A convenience function that computes the pivot position (`i0`), the forward vector (from indices in `b1`), and the up vector (from indices in `b2`) and then orients the system accordingly.
- **Example:**
  ```python
  oriented_positions = system.orient(1, (1,2), (3,4), _0=1)
  print(oriented_positions)
  ```

### `orientPCA(perm=None)`
- **Description:**  
  Reorients the system using Principal Component Analysis (PCA). An optional permutation `perm` can be provided.
- **Example:**
  ```python
  system.orientPCA()
  system.print()
  ```

### `shift(vec, sel=None)`
- **Description:**  
  Translates the system by the vector `vec`. If `sel` is provided, only those atoms are shifted.
- **Example:**
  ```python
  system.shift([1.0, 0.0, 0.0])
  ```

### `rotate_ax(ang, ax=(0,1), p0=None)`
- **Description:**  
  Rotates the system by angle `ang` (in radians) about the specified axis. If `p0` is provided, rotation is performed about that point.
- **Example:**
  ```python
  system.rotate_ax(np.pi/4, ax=(0,1))
  ```

### `delete_atoms(lst)`
- **Description:**  
  Deletes atoms specified by indices in `lst` from all per–atom arrays (positions, types, charges, etc.).
- **Example:**
  ```python
  system.delete_atoms([0, 1])
  system.print()
  ```

### `append_atoms(B, pre="A")`
- **Description:**  
  Appends the atoms from another `AtomicSystem` (B) to the current system. It also concatenates auxiliary labels using a given prefix.
- **Example:**
  ```python
  # Assuming B is another AtomicSystem (e.g., an end–group)
  system.append_atoms(B, pre="X")
  system.print()
  ```

### `remap(lst)`
- **Description:**  
  Remaps a list of auxiliary labels to their new indices after operations such as attachment.
- **Example:**
  ```python
  new_indices = system.remap(["X0", "X1", "X2"])
  print("Remapped indices:", new_indices)
  ```

### `attach_group(G, i0, i1, iup, bond, up=(0.,0.,1.), _0=1, pre="A")`
- **Description:**  
  Attaches an end–group (`G`) to the backbone (self) at a specified bond. The procedure:
  1. Orients the group in its own reference frame using:
     - **Pivot (i0):** The atom that will be placed at the attachment site.
     - **Forward (i1):** The atom used (together with the pivot) to define the forward vector. (This atom is then removed.)
     - **Up (iup):** A pair of indices defining the up vector.
  2. Computes a rotation matrix from the backbone bond (`bond`) and a provided backbone up vector (`up`).
  3. Applies the rotation and translates the group so that its pivot coincides with the backbone’s attachment site.
  4. Appends the group’s atoms to the backbone.
- **Example:**
  ```python
  # For example, attach group G using:
  #   Pivot: atom 1 in G,
  #   Forward: bond from atom 1 to atom 12 (atom 12 is removed),
  #   Up: defined by atoms (10, 11),
  #   Backbone bond: (18,8),
  #   and backbone up vector (0,0,1).
  backbone.attach_group(G, i0=1, i1=12, iup=(10,11), bond=(18,8), up=(0,0,1), _0=1, pre="X")
  backbone.print()
  ```

---

## Additional Utility Methods

### `store_bond_lengths()` and `restore_bond_length(ij, L=None)`
- **Description:**  
  These functions compute and then restore (or adjust) bond lengths. They are useful when bond lengths have been perturbed by a transformation.
- **Example:**
  ```python
  bond_lengths = system.store_bond_lengths()
  system.restore_bond_length((0,1))
  ```

---

## Examples

### 1. Loading a Molecule and Printing Its Data

```python
from pyBall.atomicUtils import AtomicSystem
system = AtomicSystem(fname='molecule.xyz')
system.print()
```

### 2. Cloning a System Under Periodic Boundary Conditions

```python
clone = system.clonePBC(nPBC=(2,2,1))
clone.print()
```

### 3. Attaching an End–Group

```python
# Load backbone and end–group systems.
backbone = AtomicSystem(fname='backbone.xyz')
group = AtomicSystem(fname='endgroups/adenine.xyz')

# (Assume per–atom properties are pre-initialized for both systems.)
# Attach the group using the defined indices:
# - Pivot atom in group: 1
# - Forward defined by bond (1, 12) in group (atom 12 will be deleted)
# - Up defined by atoms (10, 11) in group.
# - Backbone bond used for attachment: (18,8)
backbone.attach_group(group, i0=1, i1=12, iup=(10,11), bond=(18,8), up=(0,0,1), _0=1, pre="X")
backbone.print()
```

### 4. Finding Bonds and Angles

```python
bonds, distances = system.findBonds(Rcut=3.0, RvdwCut=1.5)
print("Bonds:", bonds)
angles, angle_indices = system.findAngles()
print("Angles (radians):", angles)
```

