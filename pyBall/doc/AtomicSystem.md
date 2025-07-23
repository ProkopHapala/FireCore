# `AtomicSystem.py` Documentation

This document provides a comprehensive overview of the `AtomicSystem.py` module, which defines the `AtomicSystem` class. This class serves as a fundamental data structure for representing and manipulating atomic and molecular systems within the `FireCore` project. It leverages many utility functions from `atomicUtils.py` (documented in `atomicUtils.md`).

## 1. `AtomicSystem` Class

### 1.1. Purpose

The `AtomicSystem` class encapsulates all relevant data for an atomic system, including atomic coordinates, types, element names, charges, radii, bond connectivity, neighbor lists, and simulation cell lattice vectors. It provides a rich set of methods for loading, saving, querying, analyzing, and manipulating these systems.

### 1.2. Initialization (`__init__`)

```python
AtomicSystem(fname=None, apos=None, atypes=None, enames=None, lvec=None, qs=None, Rs=None, bonds=None, ngs=None, bReadN=True, bPreinit=True)
```

- **`fname` (str, optional):** Path to a file (`.mol`, `.mol2`, `.xyz`) to load the atomic system from.
- **`apos` (np.ndarray, optional):** (N, 3) array of atomic positions.
- **`atypes` (np.ndarray, optional):** (N,) array of atomic types (integer indices).
- **`enames` (list of str, optional):** List of element names (e.g., `['C', 'H', 'O']`).
- **`lvec` (np.ndarray, optional):** (3, 3) array representing lattice vectors for periodic boundary conditions.
- **`qs` (np.ndarray, optional):** (N,) array of atomic charges.
- **`Rs` (np.ndarray, optional):** (N,) array of atomic radii.
- **`bonds` (np.ndarray, optional):** (M, 2) array of bond indices `[[i, j], ...]`. If not provided, bonds can be found later using `findBonds()`.
- **`ngs` (list of dict, optional):** List of dictionaries representing neighbor lists for each atom.
- **`bReadN` (bool, optional):** If `True`, reads the number of atoms from the file header (for some formats). Default is `True`.
- **`bPreinit` (bool, optional):** If `True` (default), calls `preinitialize_atomic_properties()` after loading/initialization to set up `iZs` (atomic numbers), `masses`, and `Rs` (radii) based on `enames`.

**Note:** When initialized from a file, `AtomicSystem` internally uses functions from `atomicUtils.py` (e.g., `au.loadMol`, `au.load_xyz`).

### 1.3. Attributes

- **`apos` (np.ndarray):** Atomic positions (N, 3).
- **`atypes` (np.ndarray):** Atomic types (N,).
- **`enames` (list of str):** Element names (N,).
- **`qs` (np.ndarray):** Atomic charges (N,).
- **`Rs` (np.ndarray):** Atomic radii (N,).
- **`bonds` (np.ndarray):** Bond connectivity (M, 2).
- **`ngs` (list of dict):** Neighbor lists for each atom.
- **`lvec` (np.ndarray):** Lattice vectors (3, 3).
- **`iZs` (np.ndarray):** Atomic numbers (N,).
- **`masses` (np.ndarray):** Atomic masses (N,).
- **`aux_labels` (list, optional):** Auxiliary labels for atoms.

### 1.4. File I/O Methods

These methods facilitate saving the current atomic system to various file formats.

- **`saveXYZ(fname, mode="w", blvec=True, comment="", ignore_es=None, bQs=True, other_lines=None)`:**
    Saves the system to an XYZ file. Can include lattice vectors in the comment and optionally save charges and radii. Internally uses `atomicUtils.saveXYZ`.
- **`save_mol(fname, title="Avogadro")`:**
    Saves the system to a MOL file. Internally uses `atomicUtils.save_mol`.
- **`save_mol2(fname, comment="")`:**
    Saves the system to a MOL2 file. Internally uses `atomicUtils.save_mol2`.
- **`toLines()`:**
    Returns a list of strings, each representing an atom in a format suitable for geometry files (e.g., `"C   0.000   0.000   0.000"`). Internally uses `atomicUtils.geomLines`.
- **`toXYZ(fout, comment="#comment", ignore_es=None, other_lines=None, bHeader=False)`:**
    Writes the system data to an already open file object `fout` in XYZ format. Internally uses `atomicUtils.writeToXYZ`.

### 1.5. Information & Query Methods

These methods provide ways to inspect and retrieve information about the atomic system.

- **`print()`:**
    Prints a summary of the system, including the number of atoms and details for each atom (index, type, element name, position, and auxiliary labels if present).
- **`getValenceElectrons()`:**
    Returns a NumPy array containing the number of valence electrons for each atom in the system.
- **`subtractValenceE(f0=-1.0, f=+1.0)`:**
    Adjusts the `qs` (charges) array by subtracting or adding valence electron counts, scaled by `f0` and `f`.
- **`printBonds()`:**
    Prints a list of all defined bonds in the system.
- **`printNeighs()`:**
    Prints the neighbor list for each atom.
- **`findBonds(Rcut=3.0, RvdwCut=1.5, RvdWs=None, byRvdW=True)`:**
    Identifies and stores bonds within the system based on interatomic distances and optional Van der Waals radii. Internally uses `atomicUtils.findBondsNP`.
- **`findHBonds(Rb=1.5, Rh=2.5, angMax=60.0, typs1={"H"}, typs2=au.neg_types_set, bPrint=False, bHbase=False)`:**
    Identifies hydrogen bonds within the system based on distance and angular criteria. Internally uses `atomicUtils.findHBondsNP`.
- **`findBondsOfAtom(ia, bAtom=False)`:**
    Returns a list of bonds connected to atom `ia`. If `bAtom` is `True`, returns atom indices, otherwise bond indices.
- **`neighs(bBond=True)`:**
    Generates and returns the neighbor list (`self.ngs`). If `bBond` is `True`, neighbors are determined by existing bonds; otherwise, by proximity (using `atomicUtils.neigh_atoms`).
- **`find_groups()`:**
    Identifies and returns a list of atom indices for each disconnected molecular fragment (bonded cluster) in the system. Internally uses `atomicUtils.selectBondedCluster`.
- **`select_by_ename(elist)`:**
    Returns a boolean mask or list of indices for atoms whose element names are present in `elist`.
- **`getNeighsOfType(selection, typ='N')`:**
    Given a `selection` of atom indices, returns the indices of their neighbors that match the specified `typ`. Internally uses `atomicUtils.findNeighsOfType`.
- **`select_by_neighType(neighs, typ='N', neighTyps={'H':(1,2)})`:**
    Selects atoms based on the types of their neighbors. Internally uses `atomicUtils.findTypeNeigh_`.
- **`findAngles(select=None, ngs=None)`:**
    Finds and returns bond angles within the system. `select` can specify a subset of atoms.
- **`findDihedral(select=None, ngs=None, neighTyp={'H'})`:**
    Finds and returns dihedral angles within the system.
- **`findCOG(apos, byBox=False)`:**
    Calculates the center of geometry (COG) for a given set of positions `apos`. Internally uses `atomicUtils.findCOG`.
- **`projectAlongBondDir(i0, i1)`:**
    Projects atomic positions along the direction defined by the bond between atoms `i0` and `i1`. Internally uses `atomicUtils.projectAlongBondDir`.
- **`store_bond_lengths()`:**
    Calculates and stores the current lengths of all bonds in `self.bond_lengths`.
- **`restore_bond_length(ij, L=None)`:**
    Restores the length of a specific bond `ij` to its stored value or a provided `L`.

### 1.6. System Manipulation & Transformation Methods

These methods allow for modification, transformation, and assembly of atomic systems.

- **`clonePBC(nPBC=(1,1,1))`:**
    Clones the current system across periodic boundary conditions specified by `nPBC` (e.g., `(2,2,1)` for a 2x2x1 supercell).
- **`symmetrized(d=0.1)`:**
    Applies symmetry operations to the system (details depend on internal implementation).
- **`selectSubset(inds)`:**
    Returns a new `AtomicSystem` object containing only the atoms specified by `inds`.
- **`selectBondedCluster(s)`:**
    Selects a bonded cluster of atoms starting from a seed set `s` and returns a new `AtomicSystem` object for that cluster. Internally uses `atomicUtils.selectBondedCluster`.
- **`makeRotMat(ip1, ip2, _0=1)`:**
    Creates a rotation matrix based on two atom indices `ip1` and `ip2`. Internally uses `atomicUtils.makeRotMat`.
- **`orient_mat(rot, p0=None, bCopy=False)`:**
    Orients the system by applying a given rotation matrix `rot`. `p0` specifies the origin of rotation. If `bCopy` is `True`, returns a new `AtomicSystem`.
- **`orient_vs(fw, up, p0=None, trans=None, bCopy=False)`:**
    Orients the system such that the `fw` (forward) and `up` vectors align with the system's axes. Internally uses `atomicUtils.orient_vs`.
- **`orient(i0, b1, b2, _0=1, trans=None, bCopy=False)`:**
    Orients the molecule by aligning three specified atom indices (`i0`, `b1`, `b2`) to a target orientation. Internally uses `atomicUtils.orient`.
- **`orientPCA(perm=None)`:**
    Orients the system by aligning its principal components of inertia with the coordinate axes. Internally uses `atomicUtils.orientPCA`.
- **`shift(vec, sel=None)`:**
    Shifts all or a `sel`ected subset of atoms by a given `vec`tor.
- **`rotate_ax(ang, ax=(0,1), p0=None)`:**
    Rotates all or a selected subset of atoms by `ang` around a specified `ax`is passing through `p0`.
- **`delete_atoms(lst)`:**
    Deletes atoms specified by the list of indices `lst` from the system. It reindexes bonds and neighbor lists to maintain consistency. Internally uses `atomicUtils.reindex_bonds` and `atomicUtils.make_reindex`.
- **`add_atom(pos, ename, atype=-1, q=0.0, R=1.0)`:**
    Adds a new atom to the system with specified `pos`ition, `ename` (element name), `atype` (type), `q` (charge), and `R` (radius).
- **`add_bond(i, j)`:**
    Adds a new bond between atoms `i` and `j`.
- **`merge(other_system, rot=None, trans=None)`:**
    Merges another `AtomicSystem` object (`other_system`) into the current system. The `other_system` can optionally be rotated and translated before merging. It adjusts bond indices accordingly.
- **`attach_group(G, i0, i1, iup, bond, up=(0., 0., 1.), _0=1, pre="A")`:**
    Attaches an end-group (`G`, which is another `AtomicSystem` object) to the current system's backbone at a specified bond. Internally uses `atomicUtils.attach_group`.
- **`attach_group_by_marker(G, markerX="Xe", markerY="He", _0=1, pre="X")`:**
    Attaches an end-group (`G`) to the current system using marker atoms (`markerX`, `markerY`) to define the attachment point and orientation. Internally uses `atomicUtils.attach_group_by_marker`.

### 1.7. Electron Pair (`Epair`) Geometry Methods

These methods are used for calculating and placing electron pairs (e.g., lone pairs, pi-electron pairs) based on atomic bonding configurations.

- **`get_atomi_pi_direction(i)`:**
    Determines the direction of a pi-orbital for atom `i`, typically used for atoms involved in double or triple bonds.
- **`make_epair_geom(i, npi, nb)`:**
    Calculates and places electron pairs around atom `i` based on its number of pi-bonds (`npi`) and number of neighbors (`nb`). It handles different geometries (e.g., like NH3, H2O, =N-, =O).
- **`place_electron_pair(i, direction, distance=0.5, ename='E', atype=200, qs=0.0, Rs=1.0)`:**
    Adds a new "electron pair" to the system. This is represented as a new atom (typically with `ename='E'`) at a calculated position relative to atom `i` and updates `apos`, `atypes`, `enames`, `qs`, and `Rs` arrays. It also adds a bond between atom `i` and the new electron pair.

## 2. Usage Example (Conceptual)

```python
import numpy as np
from pyBall.AtomicSystem import AtomicSystem

# 1. Initialize from a file
sys1 = AtomicSystem(fname="my_molecule.xyz")
sys1.findBonds() # Find bonds if not loaded from file
sys1.print()

# 2. Initialize from arrays
positions = np.array([
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0]
], dtype=np.float32)
e_names = ['C', 'H', 'H']
sys2 = AtomicSystem(apos=positions, enames=e_names)
sys2.add_bond(0, 1)
sys2.add_bond(0, 2)
sys2.printBonds()

# 3. Manipulate the system
sys1.shift(np.array([10.0, 0.0, 0.0]))
sys1.rotate_ax(np.pi/2, ax=(0,1)) # Rotate around Z-axis

# 4. Merge systems
merged_sys = AtomicSystem(apos=np.array([]), enames=[]) # Create an empty system
merged_sys.merge(sys1)
merged_sys.merge(sys2, trans=np.array([5.0, 5.0, 0.0])) # Merge sys2 with translation
merged_sys.saveXYZ("merged_system.xyz")

# 5. Find and place electron pairs (example for a specific atom)
# Assuming sys1 has an oxygen atom at index 0 and bonds are defined
# sys1.make_epair_geom(0, npi=0, nb=2) # Example for H2O-like oxygen
# sys1.saveXYZ("molecule_with_epairs.xyz")
```

This documentation provides a solid foundation for understanding and utilizing the `AtomicSystem` class for various atomic and molecular system manipulations. For detailed parameter descriptions and specific behaviors, always refer to the source code of `AtomicSystem.py` and `atomicUtils.py`.
