# `MMFFBuilder.h` Documentation

This document provides a comprehensive overview of the `MMFFBuilder.h` header file, focusing on the `MM::Builder` class and related structures. This class is designed for building, editing, and analyzing molecular topologies, and for preparing molecular systems for classical molecular mechanics simulations.

## 1. Introduction

The `MMFFBuilder.h` header file defines classes and structures necessary for creating and manipulating molecular systems for molecular mechanics simulations. The central component is the `MM::Builder` class, which provides a rich set of functionalities for:

*   Building molecular topologies from scratch or by loading from files.
*   Modifying molecular structures (adding/removing atoms, bonds, etc.).
*   Assigning atom types and force field parameters.
*   Managing molecular fragments.
*   Preparing the system for specific force fields.

This header file is part of a larger molecular mechanics software package, and it relies on other components such as force field parameter definitions and molecular geometry libraries.

## Helper Classes

## class `Atom`

* Properties:
    *   `id`: Atom Id
    *   `type`:  Atomic type (index into the force field parameter table).
    *   `frag`:  Fragment ID to which this atom belongs.
    *   `iconf`: Index of the `AtomConf` object associated with this atom.
    *   `pos`:   3D coordinates of the atom.
    *   `REQ`:   Van der Waals radius, well depth, charge (stored as a Quaternion, but used as a Vec3)
* Methods:
    *   `print`:  Debugging method printing important properties.
    *   Constructors: Default, and with specified positions, and REQ or all properties.

## class `AtomConf`

Represents a configuration of an atom, including its bonding topology (neighbors and bonding topology).

* Properties:
    *   `iatom`: Index of the atom this configuration belongs to.
    *   `n`:     Total number of neighbors (sigma bonds + pi orbitals + electron pairs + capping atoms)
    *   `nbond`: Number of sigma bonds.
    *   `npi`:   Number of pi orbitals (bonds).
    *   `ne`:    Number of electron pairs.
    *   `nH`:    Number of capping atoms (e.g., hydrogen).
    *   `neighs`: Array of bond indices to neighboring atoms. **NOTE: these are bond *indices*, not atom indices**.  Negative values indicate non-bonding neighbors such as pi orbitals, electron pairs, or capping atoms.
    *   `pi_dir`: direction of pi-orbital
* Methods:
    *   `fillIn_npi_sp3`: calculate npi for sp3 atoms
    *   `rebase`: shifts the bond indexes after insertion of bonds.
    *   `findNeigh`: searches for the given bond in the neighborlist.
    *   `addNeigh`: adds bond (or pi bond, or electron pair, or capping atom) to neighbor list
    *   `replaceNeigh`: replace a neighbor in a list by index
    *   `countBonds`: counts the number of sigma bonds
    *   `findMinBondIndex`: search the smallest bond index starting with given index
    *   `sortBonds`: sort the indexes for sigma bonds.
    *   `addBond`: add sigma bond to neighborlist.
    *   `addH`: add hydrogen to the neighborlist.
    *   `addPi`: add Pi-bond to neighborlist.
    *   `addEpair`: add e-pair to neighborlist
    *   `updateNtot`: update totla count of neighborlist.
    *   `updateNeighs`: shortcut to update neighbor properties
    *   `clearNonBond`: sets all non-bond values to 0
    *   `clearBond`: sets sigma bond to 0
    *   `setNonBond`: initializes number of non-bond neighbors
    *   `init0`: initializes object with zero values and null neighborlist
    *   `checkNeighsRepeat`: looks for repeating bonds in neighborlist.
    *   Constructor: Defauilt and constructor with number of npi and ne

### class `Bond`

Represents a bond between two atoms in a molecular structure. This class is used to store information about the bond, such as its type, the atoms it connects, and its properties.

* Properties:
    *   `type`: Type of the bond (e.g., single, double, triple). Represented as an integer.
    *   `atoms`: Indices of the two atoms connected by this bond.
    *   `l0`: Equilibrium bond length.
    *   `k`: Bond stiffness (force constant).
    *   `kpp`: pi-pi stiffness parameter - Used for describing the allignment of pi-orbitals
    *   `ipbc`: Index of the cell image in periodic boundary conditions.  This allows bonds to "wrap around" the simulation box.
    *   `order`: Bond order which can be non-integer in resonant systems.
* Methods:
    *   `getNeighborAtom`: Given an atom index involved in this bond, returns the index of the other atom. Returns -1 if the given atom is not part of the bond. This is a convenience function for accessing bonded neighbors.
    *   `print`: Prints bond information for debugging.
    *   Constructors:  Default and with specified type, atoms, l0, k, kpp.

### class `Angle`

Represents an angle between three atoms in a molecular structure. This class is used to store information about the angle, such as its type, the atoms it connects, and its properties.

* Properties:
    *   `type`: Type of the angle (for indexing parameters).
    *   `bonds`: Indices of the two bonds forming the angle.
    *   `a0`: Equilibrium angle (in radians).
    *   `k`: Angle stiffness (force constant).
    *   `C0, C1, C2, C3`: Coefficients of the expansion of the angle potential.
    *   `atoms`: Indices of the three atoms forming the angle.

### class `Dihedral`

Represents a dihedral angle.  Dihedral angles are defined by four atoms, and are used to describe the torsional stiffness of a molecule.

* Properties:
    *   `type`: Dihedral type (for parameter indexing).
    *   `bonds`: Indices of the three bonds forming the dihedral angle.
    *   `atoms`: Indices of the four atoms involved in the dihedral angle.
    *   `d`: Parameter 'd' for the dihedral angle
    *   `n`: Multiplicity of the dihedral angle. The number of repeating units over a 360-degree rotation.
    *   `k`: Torsional stiffness (force constant).
    *   `a0`: Equilibrium dihedral angle.

### class `Inversion`

Represents an improper dihedral angle.  Improper dihedrals are defined by four atoms, and are used to describe the inversion stiffness of a molecule perpendiculat to some plane.

* Properties:
    *   `type`: Inversion type (for parameter indexing).
    *   `bonds`: Indices of the three bonds connected to the central atom.
    *   `atoms`: Indices of the four atoms involved in the improper dihedral.
    *   `k`: Inversion stiffness (force constant).
    *   `C0, C1, C2`: Coefficients of the expansion of the inversion potential.

### Fragment

* Properties:
    *   `imolType`: Integer that keeps the ID for the molecule type.
    *   `atomRange`: Range of atom indices belonging to this fragment.
    *   `confRange`: Range of `AtomConf` indices belonging to this fragment.
    *   `bondRange`: Range of bond indices belonging to this fragment.
    *   `angRange`:  Range of angle indices belonging to this fragment.
    *   `dihRange`:  Range of dihedral indices belonging to this fragment.
    *   `pos`:     Position of the fragment (e.g., center of mass).  Used for rigid-body dynamics.
    *   `rot`:     Rotation of the fragment (as a quaternion).  Used for rigid-body dynamics.
    *   `pos0s`:   pointer to the array of initial positions of atoms
* Methods:
    *   `color`: Color used for visualization (e.g. in GUIs)
    *   `finish`: helper method for setting the end values for ranges

## Free Functions

### `splitGraphs()`

```cpp
int splitGraphs( int nb, Vec2i* bonds, int b0, std::unordered_set<int>& innodes );
```

This function performs a graph traversal to identify all atoms connected to a given set of starting atoms (`innodes`), *excluding* a specific bond (`b0`).  It's used for segmenting a molecule based on a bond being broken, creating disconnected fragments. This helps e.g. to select shorter segment when editing

* Arguments:
    *   `nb`: Total number of bonds in the system.
    *   `bonds`: Array of `Vec2i` representing the atom indices connected by each bond.
    *   `b0`:  Index of the bond to exclude from the graph traversal.
    *   `innodes`: Set of initial atom indices to start the traversal from.
    *   Return: Number of nodes (atoms) in resulting graph

### `splitByBond()`

```cpp
int splitByBond( int ib, int nb, Vec2i* bond2atom, Vec3d* apos, int* selection, Vec3d& ax, Vec3d& p0 );
```

This function utilizes `splitGraphs` to divide a molecular system into two fragments by severing a specific bond. It determines the atoms belonging to each fragment and populates a selection array. This helps e.g. to select shorter segment when editing

*   `ib`: Index of the bond to split the system.
*   `nb`: Total number of bonds.
*   `bond2atom`: Array mapping bond indices to atom pairs.
*   `apos`: Array of atom positions.
*   `selection`: Integer array to store the atom indices belonging to the two fragments.
*   `ax`: Reference to a `Vec3d` that will store the normalized vector along the bond being broken.
*   `p0`: Reference to a `Vec3d` that will store position of first atom for splitting.





## Main Class `MM::Builder Class`

The `MM::Builder` class is the central component of this header file. It provides a comprehensive set of methods for building, editing, and analyzing molecular topologies.

### Member Variables

#### Topology Data

*   `atoms`:  `std::vector<Atom>` - Stores the atoms in the system.
*   `bonds`:  `std::vector<Bond>` - Stores the bonds in the system.
*   `angles`: `std::vector<Angle>` - Stores the angles in the system.
*   `dihedrals`: `std::vector<Dihedral>` - Stores the dihedral angles in the system.
*   `inversions`: `std::vector<Inversion>` - Stores the improper dihedral angles (inversions) in the system.
*   `confs`: `std::vector<AtomConf>` - Stores the atom configurations (bonding topology).
*   `selection`: `std::unordered_set<int>` - Keeps selected atoms in set.
*   `atom_permut`:  `std::vector<int>` - Stores atom permutation.

#### PBC (Periodic Boundary Conditions)

*   `bPBC`: `bool` - Flag indicating whether periodic boundary conditions are enabled.
*   `lvec`: `Mat3d` - Lattice vectors defining the simulation box.

#### Force Field Parameters

* `params`: `MMFFparams*` - Pointer to the force field parameter object.
* `atomTypeNames`: `std::vector<std::string>*` - Vector of atom type names (used for I/O).
* `atomTypeDict`: `std::unordered_map<std::string,int>*` - Map from atom type names to integer indices.
* `Lepair`: `double` - Distance of lone pair
* `Kepair`: `double` - Stiffness for lone pair
* `Ksp_default`: `double` - parameter for sp
* `Kpp_default`: `double` - parameter for pp
* `Kpp_e_default`: `double` - Parameter for pp_e
* `ignoreType`: `int` - Atom type to ignore.
* `defaultREQ`: `Quat4d` - Default REQ
* `defaultBond`: `Bond` - Default Bond
* `defaultAngle`: `Angle` - Default Angle
* `bondBrush`: `Bond` - Bond type for brush mode.
* `itypHcap`: `int` - The ID for H atom
* `itypEpair`: `int` - The ID for e-pair

#### Fragment Management

*   `frags`: `std::vector<Fragment>` - Stores the molecular fragments in the system.

#### Capping and Dummy Atom Parameters

*   `capping_types`: `std::unordered_set<int>` - Set of atom types that are considered capping atoms (e.g., hydrogen).
*   `sp3types`: `std::unordered_set<int>` - Set of atom types that are considered sp3 hybridized.

#### Defaults

*   `capAtom`: `Atom` - Default parameters for capping hydrogen atoms.
*   `capAtomEpair`: `Atom` - Default parameters for capping electron pair "atoms".
*   `capAtomPi`: `Atom` - Default parameters for capping Pi "atoms".
*   `capBond`: `Bond` - Default parameters for bonds to capping atoms.

#### Flags

*   `bDummyPi`:  `bool` - Add a bond (or pi bond, or electron pair, or capping atom) into an atom neighbor list
*   `bDummyEpair`: `bool` - Enable adding dummy electron pairs.
*   `bAutoTypes`: `bool` - Flag to enable automatic atom type assignment.
*   `bAddCaps`:   `bool` - Flag to enable automatic addition of capping atoms.
*   `bAddECaps`:   `bool` - Flag to enable automatic addition of epairs
*   `capUp`: `Vec3d` - Default up direction vector.






### Member Functions


*   `makeNeighs( int*& neighs, int perAtom )`: makes Neighbor list, allocates memory and handles it also

#### Topology Editing

* **Atoms**
    *   `removeAtom(int i, bool checkBrute=true)`: Removes an atom from the system.
    *   `insertAtoms( int na, int* atypes, Vec3d* apos, Quat4d* REQs=0, double* qs=0, int* npis=0, const Vec3d& pos=Vec3dZero, const Mat3d& rot=Mat3dIdentity, const Vec3d& pos0=Vec3dZero )`: insert Atoms with atomtypes and bonds.
    *   `insertAtomsBonds( int nb, Vec2i* b2as, int* btypes, int na, int* atypes,  Vec3d* apos, Quat4d* REQs=0, double* qs=0, int* npis=0, const Vec3d& pos=Vec3dZero, const Mat3d& rot=Mat3dIdentity )`: Combine insert atoms and bonds for molecular assembly.
* **Bonds**
    *   `clearBonds()`: Clears all bonds, angles, dihedrals, and configurations from the system, while preserving the atoms.
    *   `insertBond(const Bond& bond )`: insert bond
    *   `insertBond( Vec2i ias, int order )`: inserts bond with certain bond order
    *   `insertBonds( int nb, Vec2i* b2as, int* btypes )`: Insert bonds into molecule.
    *   `removeBondFromConfs(int ib)`: remove bonds from configuration
    *   `autoBonds( double R=-1.20, int i0=0, int imax=-1 )`: Finds and creates bonds between atoms based on a distance criterion. Uses a simple cutoff radius.
    *   `autoBondsPBC( double R=-1.35, int i0=0, int imax=-1, Vec3i npbc=Vec3iOne )`: Finds and creates bonds between atoms considering periodic boundary conditions. This function is essential for simulating extended systems.
    *   `touchingAtoms( int i0, int imax, const Vec3d& p, double R0, double Rfac, std::vector<int>& found )`: finds the set of touching Atoms for given search point.
* **Neihbors** (conectivity graph)
    *   `findHighestValence()`: Finds the atom with the highest valence (number of bonds).
    *   `findNthAtomOfType( int ityp, int n)`: Finds the index of the *n*-th atom of a given type.
* **Angles**
    *   `addAnglesToBond( int ib, int n,const int* neighs, double a0, double k )`: add angles
    *   `addAnglesUpToN( int n, const int* neighs, double a0, double k )`: add angles
    *   `addAnglesToAtom( int ia, double ksigma, double kpi )`: add angles for selected atom
* **Capping Atoms**
    *   `addEpairsToAtoms(int ia, double l=0.5 )`: Deprecated function for adding epairs
    *   `addCappingTypesByIz( int iZ )`: add atom types by atomic number Z to cap type
    *   `addCaps( int ia, int ncap, int ne, int nb, const Vec3d* hs )`: Adds capping atoms (e.g., hydrogen) to an atom.
    *   `addCapTopo(int ia)`: adds capping
* **Other Topology operations**
    *   `clear()`: Clears all data from the builder, including atoms, bonds, angles, dihedrals, configurations, and fragments.
    
#### Geometry Manipulation

* **Transforms** (rotation, translation)
    *   `move_atoms( Vec3d dshift, int i0=0, int imax=-1 )`:  Translates a range of atoms by a given shift vector.
    *   `rotationFromAtoms( int i0, int i1, int i2 )`: compute orientation of the given three atom
    *   `transform_atoms( Mat3d M, Vec3d orig_old=Vec3dZero, Vec3d orig_new=Vec3dZero, int i0=0, int imax=-1 )`: Transforms a range of atoms using a 3x3 matrix `M`, with optional origin translation.
    *   `rotate_atoms( double angle, Vec3d axis=Vec3dZ, Vec3d orig_old=Vec3dZero, Vec3d orig_new=Vec3dZero, int i0=0, int imax=-1 )`: Rotates a range of atoms around a given axis by a specified angle, with optional origin translation.
    *   `orient_atoms( Vec3d fw, Vec3d up, Vec3d orig_old=Vec3dZero, Vec3d orig_new=Vec3dZero, int i0=0, int imax=-1 )`: Orients a range of atoms to a specified forward and up direction, with optional origin translation.
* **Crystal lattice operations**
    *   `pbcShift( Vec3i G )`: calculate PBC
    *   `changeCell( const Mat3d& lvs, Vec3d orig_old=Vec3dZero, Vec3d orig_new=Vec3dZero, int i0=0, int n=-1 )`: Changes the simulation cell (lattice vectors) and transforms the atom coordinates accordingly for PBC.
    *   `updatePBC( Vec3d* pbcShifts, Mat3d* M=0 )`: Update PBC to new configuration
* **Other Geometry operations**
    *   `vecBetweenAtoms(int i, int j)`: Returns a vector pointing from atom `i` to atom `j`.
    *   `bbox( Vec3d& pmin, Vec3d& pmax, int i0=0, int n=-1, bool bInit=true)`:  Calculates the bounding box of a set of atoms.
    *   `findMainAxes(int i0=0,int imax=-1, const bool bRemoveCog=true, const bool bRot=true, Vec3i permut=Vec3i{2,1,0} )`: Calculates and aligns the molecule along its principal axes of inertia.
    *   `findSymmetry( int* found, int i0=0,int imax=-1, double tol=0.1 )`: Finds mirror and point symmetry operations.

#### File I/O

*   `bindParams( MMFFparams* params_ )`: Binds a force field parameter object to the builder.  This allows the builder to access force field parameters for atom types, bonds, angles, etc.
*   `initDefaultAtomTypeDict()`: Creates a default atom type dictionary (mapping atom type names to indices).
*   `load_xyz( const char * fname, bool noH=false, bool bConf=true )`: Loads atom coordinates and types from an XYZ file.
*   `write2xyz( FILE* pfile, const char* comment="#comment" )`: Writes atom coordinates and types to an XYZ file.
*   `save2xyz( const char * fname, const char* comment="#comment" )`: Saves atom coordinates and types to an XYZ file.
*   `saveMol( const char* fname )`: save molecule as MOL file
*   `loadMolTypeXYZ(const char* fname, MMFFparams* params_=0 )`: loads specific molecules types from file.
*   `loadXYZ_Atoms(const char* fname, MMFFparams* params_=0, int iH=-1, bool bCOG=false, const Vec3d& pos=Vec3dZero, const Mat3d& rot=Mat3dIdentity )`: loads atoms from XYZ

#### Re-ordering (re-numbering) of Atoms and Bonds etc.

*   `sortConfAtomsFirst()`: moves all configured atoms (iconf>=0) to the beginnig of array, non-configured (capping) to the end.
*   `sortAtomsOfBonds()`: sorts atoms in bond
*   `sortBonds()`: sorts bonds
*   `setup_atom_permut( bool bPrint=false )`: build atom_permut table
*   `permutAtoms(int* permut, bool doBonds=false )`: Reorder atoms within the class
*   `reindexConfs(int* ias, int* ibs)`: Reindex the atoms to a different array.
*   `reindexBonds(int* ias, bool bUpdateConfs=true )`: Reindex Bonds.
*   `numberAtoms()`: assign id of atoms by number in the atoms array

#### Check Topology Consistency

*   `checkAllAtomsBonded( bool bPrint=true, bool bExit=true, int nbmax=N_NEIGH_MAX, int nbmin=1 )`: Check all atoms bonded
*   `checkBondsOrdered( bool bOrder, bool bPrint )`: check if b.a<b.b
*   `checkAtomHasBond(int ia, int ib, bool bDefault=true)const`: check if given bond is in neighborhood.
*   `checkNeighsRepeat( bool bPrint=true )`: check repeating neighbors
*   `checkBondsInNeighs( bool bPrint=true )`: Check if every atom in b.atoms is in the neighborhood
*   `checkBondsSorted( int iPrint=0 )const`: Check if Bonds are sorted

#### Type and ForceField assignment

* **Type Assignment**
    *   `assignSp3Params( int ityp, int nb, int npi, int ne, int npi_neigh, Quat4d& par )`: Assigns parameters for sp3-hybridized atoms based on their local environment.
    *   `assignSp3Type( int ityp_old, int nb, int npi, int ne, int npi_neigh )`: assign sp3 type
    *   `assignAllSp3Types()`: Automatically assigns sp3 hybridization states to atoms based on their bonding topology.
    *   `assignSpecialTypes( int* neighs )`: Assigns atom types based on connectivity and surrounding atoms (using rules for common functional groups).
    *   `assignSpecialTypesLoop( int nmax, int* neighs )`: Runs the `assignSpecialTypes` function iteratively until convergence or a maximum number of iterations is reached.
    *   `assignTypes( int* neighs=0, int niterMax=10, bool bDeallocNeighs=true )`: Performs full atom type assignment, including sp3 hybridization and special type assignment.
* **Parameter assignment**
    *   `assignBondParams( int ib )`: Assign bond parameters (equilibrium length and stiffness) based on the atom types of the connected atoms and the bond order.
    *   `assignAllBondParams()`: Assigns bond parameters to all bonds in the system.
    *   `assignBondParamsUFF( int ib )`: Assigns UFF bond parameters (equilibrium length and stiffness) using the UFF force field rules.
    *   `assignTorsions( bool bNonPi=false, bool bNO=true )`: Assign Torsions
    *   `autoAngles(double ksigma, double kpi)`: add angles to the whole molecule
    *   `addTorsionsToBond( int ibond )`: add Torsions to selected Bond
    *   `chargeByNeighbors( bool bClean, double factor=0.05, int niters=1, double damp=0.5 )`: Assign charges using neighbor

#### Selection

*   `selectAll( )`: add all atons to the list
*   `selectInverse()`: inverts atoms in list
*   `selectRect( const Vec3d& p0, const Vec3d& p1, const Mat3d& rot )`: add atoms into set which are in box
*   `selectCaping()`: select atoms with capping type
*   `selectBondsBetweenTypes( int imin, int imax, int it1, int it2, bool byZ=false, bool bOnlyFirstNeigh=false )`: Selects bonds where the connected atoms have specific types.

#### Mouse Picking

*   `rayBonds( const Vec3d& ro, const Vec3d& rd, double R )`: Finds the bond that intersects with a ray defined by origin `ro` and direction `rd`, within a radius `R`. This is useful for picking bonds with a mouse cursor in a GUI.
*   `rayPickBond( const Vec3d& ro, const Vec3d& rh, double Rmax )`: select bonds with the mouse by tracing the array
*   `pickBond( const Vec3d& ro, const Vec3d& rh, double Rmax )`: use ray to select bonds

#### Visualization

*   `randomFragmentCollors()`: Assigns random colors to the fragments in the system for visualization purposes.




