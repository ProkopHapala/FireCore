# `atomicUtils.py` Documentation

This document provides a concise overview of the `atomicUtils.py` module.

## Free Functions

### Topology Operations

*   `findAllBonds(atoms, Rcut=3.0, RvdwCut=0.7)`: Determines all bonds based on distance and Van der Waals radii.
*   `filterBonds(bonds, enames, ignore)`: Filters a list of bonds, excluding those that involve atoms with element names in the `ignore` set.
*   `neigh_bonds(natoms, bonds)`: Creates a dictionary-based neighbor list where neighbors are related to atom by bond index.
*   `neigh_atoms(natoms, bonds)`: Creates a simple set-based neighbor list of atom indexes, where neighbors are defined as atoms that are within chemical bond distance
*   `addGroup(base, group, links)`: Adds a group of atoms and bonds to an existing system, connecting them based on a specified list of links.
*   `addBond(base, link, bNew=True)`: Adds a bond to a base system.
*   `disolveAtom(base, ia)`: Dissolves (removes) an atom from a base system, redistributing its bonds.
*   `removeGroup(base, remove)`: Removes a group of atoms from a base system.
*   `selectBondedCluster(s, bonds)`: Selects a bonded cluster of atoms starting from a seed set.

### Finding and Identifying Atoms/Bonds

*   `findBondsNP(apos, atypes=None, Rcut=3.0, RvdwCut=1.5, RvdWs=None, byRvdW=True)`: Efficiently determines bonds based on interatomic distances, optionally considering Van der Waals radii for more accurate bond identification.
*   `findHBondsNP(apos, atypes=None, Rb=1.5, Rh=2.5, angMax=60.0, typs1={"H"}, typs2=neg_types_set, bPrint=False, bHbase=False)`: Identifies hydrogen bonds based on distance and angular criteria involving donor, hydrogen, and acceptor atoms.
*   `findNeighOfType(ia, atypes, neighs, typ='N')`: Given an atom index, return indexes of neighbors with given type.
*   `findNeighsOfType(selection, atypes, neighs, typ='N')`: Find specific type neighbors for many atoms in a given selection
*   `findTypeNeigh(atoms, neighs, typ, neighTyps=[(1, 2, 2)])`: Identifies atoms of one type close to neighbors with set of different types (uses simple array as neigh list)
*   `findTypeNeigh_(types, neighs, typ='N', neighTyps={'H': (1, 2)})`: Identifies atoms of one type close to neighbors with set of different types (uses dictionary to store atom neighbors)
*   `getAllNeighsOfSelected(selected, neighs, atoms, typs={1})`: Gets the neighbors of a selected atom that have the requested type.
*   `findPairs(select1, select2, atoms, Rcut=2.0)`: Returns indexes of pairs of atoms between selected sets that are in the given distance
*   `findPairs_one(select1, atoms, Rcut=2.0)`: Find pairs of atoms in one set that are within the given distance.
*   `countTypeBonds(atoms, ofAtoms, rcut)`: Counts the number of atoms of given type in neighborhood.
*   `findBondsTo(atoms, typ, ofAtoms, rcut)`: Finds Bonds To atoms Of selected Type typ. Returns list of atom indexes and vectors.
*   `findNearest(p, ps, rcut=1e+9)`: Given point p and radius rcut searches nearest neighbor.
*   `pairsNotShareNeigh(pairs, neighs)`: Selects pairs that do not share neighbors, useful to look at distant correlation (or lack of it).

### Geometry Manipulation

*   `rotMatPCA(ps, bUnBorm=False)`: Calculates a rotation matrix based on principal component analysis (PCA) of a set of points.
*   `makeRotMatAng(ang, ax=(0, 1))`: Generates a rotation matrix to perform rotations by given angle along selected axis.
*   `makeRotMat(fw, up)`: Returns rotation matrix that maps selected vectors fw, up to the axis system.
*   `makeRotMatAng2(fw, up, ang)`: Generates rotated forward and up vectors. Returns rotation matrix to the axis system.
*   `rotmat_from_points(ps, ifw=None, iup=None, fw=None, up=None, _0=1)`: Creates rotation matrix from points in system.
*   `mulpos(ps, rot)`: Transforms positions by rotation.
*   `orient_vs(p0, fw, up, apos, trans=None)`: Orients the system so that fw and up point along system axis.
*   `orient(i0, ip1, ip2, apos, _0=1, trans=None, bCopy=True)`: Orients the molecule by giving 3 indexes.
*   `orientPCA(ps, perm=None)`: aligns system to it's principal axis.
*   `groupToPair(p1, p2, group, up, up_by_cog=False)`: Group given atoms so it looks along selected positions
*   `projectAlongBondDir(apos, i0, i1)`: Projects a set of positions along the direction of a specified bond.

### System Assembly

*   `replacePairs(pairs, atoms, group, up_vec=(np.array((0.0, 0.0, 0.0)), 1))`: replace atoms with a given group. Orient the groups properly
*   `replace(atoms, found, to=17, bond_length=2.0, radial=0.0, prob=0.75)`: replace type of atoms which has the bond. This is typically use to replace groups of atoms with predefined bond lengths.
*   `extract_marker_pairs(self, markerX, markerY, remove=True)`: Extracts pairs of atoms in this system based on element types.. This is usefull to attach more molecules to the system. It finds a pair of markers, based on this pair it defines new bond which will be built.
*   `attach_group(self, G, i0, i1, iup, bond, up=(0., 0., 1.), _0=1, pre="A")`: Attaches an end-group (G) to the backbone (self) at a specified bond.
*   `attach_group_by_marker(self, G, markerX="Xe", markerY="He", _0=1, pre="X")`: Attaches an end-group G to this backbone using marker atoms and connectivity.

### File I/O

*   `saveAtoms(atoms, fname, xyz=True)`: Saves a list of atoms to a file.
*   `writeToXYZ(fout, es, xyzs, qs=None, Rs=None, comment="#comment", bHeader=True, ignore_es=None, other_lines=None)`: Writes atomic data to an XYZ file format.
*   `saveXYZ(es, xyzs, fname, qs=None, Rs=None, mode="w", comment="#comment", ignore_es=None, other_lines=None)`: Saves atomic data to an XYZ file.
*   `makeMovie(fname, n, es, func)`: Creates a movie file by writing a series of XYZ frames.
*   `loadAtomsNP(fname=None, fin=None, bReadN=False, nmax=10000, comments=None)`: Loads atomic data from a file.
*   `load_xyz(fname=None, fin=None, bReadN=False, bReadComment=True, nmax=10000)`: Loads atomic data from an XYZ file.
*   `loadMol(fname=None, fin=None, bReadN=False, nmax=10000)`: Loads atomic data and bond information from a MOL file.
*   `loadMol2(fname, bReadN=True, bExitError=True)`: Loads atomic data and bond information from a MOL2 file.
*   `readAtomsXYZ(fin, na)`: Reads atomic data (coordinates and element symbols) from an already opened XYZ file.
*   `read_lammps_lvec(fin)`: Reads lattice vectors from LAMMPS trajectory.
*   `readLammpsTrj(fname=None, fin=None, bReadN=False, nmax=100, selection=None)`: Reads LAMMPS trajectory file.
*   `loadAtoms(name)`: Loads atomic data from file.
*   `load_xyz_movie(fname)`: Loads a series of frames from an XYZ movie file.
*   `psi4frags2string(enames, apos, frags=None)`: Converts atomic data and fragment information into a string format suitable for Psi4.
*   `geomLines(apos, enames)`: Converts atomic data into a list of lines suitable for writing to a geometry file.
*   `scan_xyz(fxyzin, callback=None, kwargs=None)`: Scans an XYZ file, applying a callback function to each frame.

### Auxiliary / Support

*   `string_to_matrix(s, nx=3, ny=3, bExactSize=False)`: Converts string of numbers to matrix.
*   `tryAverage(ip, apos, _0=1)`: Tries to average a set of positions, handling both single atom and multi-atom references.
*   `makeVectros(apos, ip0, b1, b2, _0=1)`: Makes forward and upward vectors.
*   `loadElementTypes(fname='ElementTypes.dat', bDict=False)`: Loads element type parameters from a file.
*   `getVdWparams(iZs, etypes=None, fname='ElementTypes.dat')`: Gets Van der Waals parameters from a list of atomic numbers.
*   `iz2enames(iZs)`: Converts a list of atomic numbers to a list of element names.
*    `histR(ps, dbin=None, Rmax=None, weights=None)`: Calculates a radial distribution histogram of a set of points, optionally weighting each point.
*   `build_frame(forward, up)`: Build an orthonormal frame (a 3×3 rotation matrix) from two non–colinear vectors.
*   `find_attachment_neighbor(system, marker_index, markerX, markerY)`: Find the non-marker attachment neighbor of a marker atom.
*   `compute_attachment_frame_from_indices(ps, iX, iY, system, bFlipFw=False, _0=1)`: Compute the attachment frame for a system from given indices.
*   `loadCoefs(characters=['s'])`: Loads coefficients from files (specific to certain quantum chemistry outputs).
*   `findCOG(ps, byBox=False)`: Calculates the center of geometry (COG) of a set of points.
*   `convert_to_adjacency_list(graph)`: Converts an adjacency matrix graph representation to an adjacency list.
*   `preprocess_graph(graph)`: Simplifies a graph by iteratively removing leaf nodes.
*   `find_cycles(graph, max_length=7)`: Locates all cycles within a graph up to a specified length.
*   `atoms_symmetrized(atypes, apos, lvec, qs=None, REQs=None, d=0.1)`: Symmetrizes atoms in a unit cell by replicating atoms near the cell boundaries.

## Class `AtomicSystem`

### Data Members

*   `apos`: NumPy array of shape `(N, 3)` representing atom positions.
*   `atypes`: NumPy array of shape `(N,)` representing atomic numbers.
*   `enames`: List of element names corresponding to each atom.
*   `qs`: NumPy array of shape `(N,)` representing atomic charges.
*   `Rs`: NumPy array of shape `(N,)` representing atomic radii.
*   `bonds`: NumPy array of shape `(M, 2)` representing bond indices.
*   `ngs`: Neighbor list (as created by `neigh_bonds` or `neigh_atoms`).
*   `lvec`: NumPy array of shape `(3, 3)` representing lattice vectors.
*   `aux_labels`: List of auxiliary labels for each atom.

### Methods

### Construction & File I/O

*   `__init__(self, fname=None, apos=None, atypes=None, enames=None, lvec=None, qs=None, Rs=None, bonds=None, ngs=None, bReadN=True)`: Initializes an `AtomicSystem` object, populating the object with data from a file or provided arrays. The init determines file type by the extention and automatically reads the data.
*   `saveXYZ(self, fname, mode="w", blvec=True, comment="", ignore_es=None, bQs=True, other_lines=None)`: Saves the current atomic structure to an XYZ file, optionally including lattice vectors, comments, charges, and radii. The versatility is controled by arguments.
*   `save_mol(self, fname, title="Avogadro")`: Saves atomic structure to MDL MOL V2000 format.
*   `save_mol2(self, fname, comment="")`: Saves atomic structure to MOL2 format.
*   `toLines(self)`: returns structure as lines of text.
*   `toXYZ(self, fout, comment="#comment", ignore_es=None, other_lines=None, bHeader=False)`: Writes XYZ file to an opened file.

### Topology Operations

*   `findBonds(self, Rcut=3.0, RvdwCut=1.5, RvdWs=None, byRvdW=True)`: Finds bonds based on distance criteria and stores them in the `bonds` attribute. The identified bonds are essential for various topology related analyzes such as ring finding or neighbor finding..
*   `findHBonds(self, Rb=1.5, Rh=2.5, angMax=60.0, typs1={"H"}, typs2=neg_types_set, bPrint=False, bHbase=False)`: Finds hydro bonds based on distance and angle criteria and stores them in the `bonds` attribute. It returns hydroen bonding patters in the system.
*   `findBondsOfAtom(self, ia, bAtom=False)`: Finds Bonds of given Atom
*   `neighs(self, bBond=True)`: Creates a neighbor list based on bonds.
*   `find_groups(self)`: find connected groups.
*   `delete_atoms(self, lst)`: Delete selected atoms. It simply deletes selected atoms.
*   `append_atoms(self, B, pre="A")`: Appends the atoms from another `AtomicSystem` to the current one. This can be used to append group of atoms to another. The groups can also be re-labeled.  It combines two atomic systems into one system

### Property and Type Assignment

*   `getValenceElectrons(self)`: gets vector of valence electrons.
*   `subtractValenceE(self, f0=-1.0, f=+1.0 )`: adjusts atomic charges by subtrating valence electrons to given atoms.
*   `preinitialize_atomic_properties(self)`: Preinitializes per-atom arrays (qs, Rs, aux_labels). It initializes the auxilary properies for the atoms. This function initializes default values and labels, ensuring required fields are properly filled in the class for downstream calculations.
*   `check_atomic_properties(atomicSystem)`: Check atomic Properties in system. Ensures that all atoms have all properties assigned to them

### Selection

*   `select_by_ename(self, elist)`: selects atoms of given ename.
*   `getNeighsOfType(self, selection, typ='N')`: select neighbors by type.
*   `select_by_neighType(self, neighs, typ='N', neighTyps={'H': (1, 2)})`: select by neighbor type.
*   `selectSubset(self, inds)`: create subset of system by indexes.
*   `selectBondedCluster(self, s)`: Selects bonded atoms from system by giving seed set s

### Geometry

*   `findAngles(self, select=None, ngs=None)`: Finds angle for system.
*   `findDihedral(self, select=None, ngs=None, neighTyp={'H'})`: Finds dihedrals for system.
*   `findCOG(self, apos, byBox=False)`: Calculates the center of geometry (COG)
*   `projectAlongBondDir(self, i0, i1)`: Projects a set of positions along the direction of a specified bond.
*   `makeRotMat(self, ip1, ip2, _0=1)`: Makes rotation matrix based on selected points for the system.
*   `orient_mat(self, rot, p0=None, bCopy=False)`: Orient system to specific matrix and origin point p0.
*   `orient_vs(self, fw, up, p0=None, trans=None, bCopy=False)`: Orients the `AtomicSystem` to a specified forward and up direction, relative to point p0.
*   `orient(self, i0, b1, b2, _0=1, trans=None, bCopy=False)`: Orients the `AtomicSystem`, b1, b2 specify to Bonds.
*   `orientPCA(self, perm=None)`: Orients the system by principal component analysis (PCA)
*   `shift(self, vec, sel=None)`: shifts position.
*   `rotate_ax(self, ang, ax=(0, 1), p0=None)`: rotates coordinates to given point

### System Assembly

*   `store_bond_lengths(self)`: Stores the bond lengths of all bonds in the system.
*   `restore_bond_length(self, ij, L=None)`: Restores the bond length of a bond.
*   `attach_group(self, G, i0, i1, iup, bond, up=(0., 0., 1.), _0=1, pre="A")`: Attaches an end-group (G) to the backbone (self) at a specified bond.
*   `extract_marker_pairs(self, markerX, markerY, remove=True)`: Extracts pairs of atoms in this system based on element types..
*   `attach_group_by_marker(self, G, markerX="Xe", markerY="He", _0=1, pre="X")`: Attaches an end-group G to this backbone using marker atoms and connectivity. The function enables joining the given fragment to selected marker, by marker are meant specifically designed types of atoms which are specified with markerX and markerY, which defines anchor point and direction of connection.

### Periodic Boundary Conditions

*   `clonePBC(self, nPBC=(1, 1, 1))`: Clones the `AtomicSystem` to create a periodic supercell.
*   `symmetrized(self, d=0.1)`: Symmetrizes atoms in a unit cell by replicating atoms near the cell boundaries.

### Other

*   `remap(self, lst)`: remap atom's label to its index