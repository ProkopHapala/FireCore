# UFF.h

## Includes

- `<omp.h>`: OpenMP library for parallel programming.
- `"fastmath.h"`: Provides fast mathematical operations.
- `"Vec2.h"`: Defines 2D vector structures and operations.
- `"Vec3.h"`: Defines 3D vector structures and operations.
- `"quaternion.h"`: Implements quaternion mathematics, useful for rotations.
- `"Forces.h"`: Contains various physical interaction forces like Lennard-Jones potential.
- `"SMat3.h"`: Symmetric matrix class for storing and manipulating symmetric matrices.
- `"molecular_utils.h"`: Utility functions related to molecular structures.
- `"NBFF.h"`: Non-Bonded Force Field implementation.
- `"GOpt.h"`: Global optimization utilities.
- `"Buckets.h"`: Data structure for efficient neighbor list management.

## Free Functions

### `bool checkVec3Match (Vec3d f, Vec3d f_, const char* label, int iPrint=1 )`

Checks if two 3D vectors are approximately equal. Prints the cosine similarity and relative length ratio of the vectors if `iPrint > 0`.

- **Parameters:**
  - `Vec3d f`: First vector to compare.
  - `Vec3d f_`: Second vector to compare.
  - `const char* label`: Label for debugging output.
  - `int iPrint`: Verbosity level (default is 1).
  
- **Returns:** 
  - `bool`: True if vectors are approximately equal, false otherwise.

### `bool checkVec3Matches (int n, Vec3d* v, Vec3d* v_, const char* label, int iPrint=1 )`

Checks multiple pairs of 3D vectors for approximate equality. Prints detailed information about each pair if `iPrint > 0`.

- **Parameters:**
  - `int n`: Number of vector pairs to check.
  - `Vec3d* v`: Array of first vectors.
  - `Vec3d* v_`: Array of second vectors.
  - `const char* label`: Label for debugging output.
  - `int iPrint`: Verbosity level (default is 1).
  
- **Returns:** 
  - `bool`: True if all pairs are approximately equal, false otherwise.

## Types (Classes and Structs)

### class `UFF`

The UFF class implements the Universal Force Field (UFF) for molecular mechanics. It provides a framework to evaluate intramolecular interactions like bonds, angles, dihedrals, and inversions using various algorithms.

**Inheritance**

- Inherits from `NBFF` which handles non-bonded forces in molecular systems.

#### Properties

- **public: double Etot**: Total energy of the system.
- **int nbonds**: Number of bond interactions.
- **int i0dih**: Index offset for dihedral interactions.
- **Mat3d invLvec**: Inverse lattice vectors used for periodic boundary conditions.
- **double SubNBTorsionFactor**: Factor to subtract torsion energy from non-bonded energy if greater than zero.
- **Vec3d * fbon**: Forces on atoms due to bonds (not currently in use).
- **Vec3d * fang**: Forces on atoms due to angles (not currently in use).
- **Vec3d * fdih**: Forces on atoms due to dihedrals.
- **Vec3d * finv**: Forces on atoms due to inversions.
- **Buckets a2f**: Data structure for efficient neighbor list management.

#### Methods

##### `void realloc (int natoms_, int nbonds_, int nangles_, int ndihedrals_, int ninversions_ )`

Reallocates memory and initializes arrays based on the number of atoms, bonds, angles, dihedrals, and inversions.

- **Parameters:**
  - `int natoms_`: Number of atoms.
  - `int nbonds_`: Number of bond interactions.
  - `int nangles_`: Number of angle interactions.
  - `int ndihedrals_`: Number of dihedral interactions.
  - `int ninversions_`: Number of inversion interactions.

##### `void dealloc ()`

Deallocates all memory used by the UFF class to prevent memory leaks.

##### `void setLvec (const Mat3d& lvec_)`

Sets the lattice vectors for periodic boundary conditions.

- **Parameters:**
  - `const Mat3d& lvec_`: Lattice vectors defining the unit cell.

##### `void mapAtomInteractions ()`

Maps interactions between atoms to their respective force pieces in the `a2f` buckets structure.

##### `void makeNeighBs ()`

Creates neighbor bond indices for each atom, indicating which bonds it participates in.

##### `void bakeDihedralNeighs ()`

Precomputes neighbor indices for dihedral interactions based on periodic boundary conditions.

##### `void bakeAngleNeighs ()`

Precomputes neighbor indices for angle interactions based on periodic boundary conditions.

##### `void bakeInversionNeighs ()`

Precomputes neighbor indices for inversion interactions based on periodic boundary conditions.

##### `void cleanForce ()`

Sets all forces to zero, preparing the system for a new evaluation or optimization step.

##### `void makeNeighCells (const Vec3i nPBC_ )`

Creates a list of neighbors' cell indices in periodic boundary conditions by iterating over all possible images.

- **Parameters:**
  - `const Vec3i nPBC_`: Number of periodic boundary cells along each axis.

##### `void makeNeighCells (int npbc, Vec3d* pbc_shifts )`

Creates a list of neighbors' cell indices using precomputed PBC shifts.

- **Parameters:**
  - `int npbc`: Number of periodic boundary images.
  - `Vec3d* pbc_shifts`: Precomputed shift vectors for each image.

##### `void assembleForcesDebug (bool bbonds, bool bangles, bool bdihedrals, bool binversions)`

Prints detailed information about the forces assembled from bonds, angles, dihedrals, and inversions. Useful for debugging purposes.

- **Parameters:**
  - `bool bbonds`: Flag to print bond forces.
  - `bool bangles`: Flag to print angle forces.
  - `bool bdihedrals`: Flag to print dihedral forces.
  - `bool binversions`: Flag to print inversion forces.

##### `void assembleForces ()`

Assembles the total force on each atom by summing up contributions from bonds, angles, dihedrals, and inversions. This function is not parallelized due to dependencies between atoms.

##### `void printForcePieces ()`

Prints detailed information about the assembled forces for debugging purposes.

##### `void assembleAtomForce (const int ia)`

Assembles the total force on a single atom by summing up contributions from bonds, angles, dihedrals, and inversions. This function is not parallelized due to dependencies between atoms.

- **Parameters:**
  - `const int ia`: Index of the atom whose forces are being assembled.

##### `void assembleAtomsForces ()`

Iterates over all atoms in the system and calls `assembleAtomForce` for each one, assembling the total force on every atom.

##### `double evalAtomBonds (const int ia, const double R2damp, const double Fmax2)`

Evaluates bond interactions between an atom and its neighbors. Computes both energy and forces acting on the atoms involved in bonds.

- **Parameters:**
  - `const int ia`: Index of the central atom.
  - `const double R2damp`: Damping parameter for Lennard-Jones potential.
  - `const double Fmax2`: Maximum force magnitude squared to clamp non-bonded interactions.

##### `double evalBonds ()`

Evaluates all bond interactions in the system by calling `evalAtomBonds` for each atom. Accumulates total energy and forces from bonds.

##### `double evalAngle_Prokop (const int ia, const double R2damp, const double Fmax2)`

Evaluates angle interactions using Prokop's method. Computes both energy and forces acting on the atoms involved in angles.

- **Parameters:**
  - `const int ia`: Index of the central atom.
  - `const double R2damp`: Damping parameter for Lennard-Jones potential.
  - `const double Fmax2`: Maximum force magnitude squared to clamp non-bonded interactions.

##### `double evalAngle_Paolo (const int ia, const double R2damp, const double Fmax2)`

Evaluates angle interactions using Paolo's method. Computes both energy and forces acting on the atoms involved in angles.

- **Parameters:**
  - `const int ia`: Index of the central atom.
  - `const double R2damp`: Damping parameter for Lennard-Jones potential.
  - `const double Fmax2`: Maximum force magnitude squared to clamp non-bonded interactions.

##### `double evalAngles ()`

Evaluates all angle interactions in the system by calling either `evalAngle_Prokop` or `evalAngle_Paolo` depending on which is enabled. Accumulates total energy and forces from angles.

##### `double evalDihedral_Prokop (const int id, const bool bSubNonBond, const double R2damp, const double Fmax2)`

Evaluates dihedral interactions using Prokop's method. Computes both energy and forces acting on the atoms involved in dihedrals.

- **Parameters:**
  - `const int id`: Index of the central dihedral.
  - `const bool bSubNonBond`: Flag to subtract non-bonded interactions from torsion energy.
  - `const double R2damp`: Damping parameter for Lennard-Jones potential.
  - `const double Fmax2`: Maximum force magnitude squared to clamp non-bonded interactions.

##### `double evalDihedral_Prokop_Old (const int id, const bool bSubNonBond, const double R2damp, const double Fmax2)`

Evaluates dihedral interactions using an older version of Prokop's method. Computes both energy and forces acting on the atoms involved in dihedrals.

##### `double evalDihedral_Paolo (const int id, const bool bSubNonBond, const double R2damp, const double Fmax2)`

Evaluates dihedral interactions using Paolo's method. Computes both energy and forces acting on the atoms involved in dihedrals.

- **Parameters:**
  - `const int id`: Index of the central dihedral.
  - `const bool bSubNonBond`: Flag to subtract non-bonded interactions from torsion energy.
  - `const double R2damp`: Damping parameter for Lennard-Jones potential.
  - `const double Fmax2`: Maximum force magnitude squared to clamp non-bonded interactions.

##### `double evalDihedrals ()`

Evaluates all dihedral interactions in the system by calling either `evalDihedral_Prokop` or `evalDihedral_Paolo` depending on which is enabled. Accumulates total energy and forces from dihedrals.

##### `double evalInversion_Prokop (const int ii)`

Evaluates inversion interactions using Prokop's method. Computes both energy and forces acting on the atoms involved in inversions.

- **Parameters:**
  - `const int ii`: Index of the central inversion.

##### `double evalInversion_Paolo (const int ii)`

Evaluates inversion interactions using Paolo's method. Computes both energy and forces acting on the atoms involved in inversions.

- **Parameters:**
  - `const int ii`: Index of the central inversion.

##### `double evalInversions ()`

Evaluates all inversion interactions in the system by calling either `evalInversion_Prokop` or `evalInversion_Paolo` depending on which is enabled. Accumulates total energy and forces from inversions.

##### `double eval (bool bClean=true)`

Evaluates the full UFF intramolecular force-field, including bonds, angles, dihedrals, and inversions. Optionally cleans forces before evaluation.

- **Parameters:**
  - `bool bClean`: Flag to clean forces before evaluating interactions.
  
##### `double eval_omp_old (bool bClean=true)`

Evaluates the full UFF intramolecular force-field using OpenMP for parallelization. This function is deprecated and not recommended for use in new code.

##### `double eval_omp (bool bClean=true)`

Evaluates the full UFF intramolecular force-field using OpenMP for parallelization, optimizing performance by minimizing dependencies between atoms.

- **Parameters:**
  - `bool bClean`: Flag to clean forces before evaluating interactions.

##### `int run (int niter, double dt, double Fconv, double Flim, double damping=0.1)`

Runs the optimization algorithm to minimize the total energy of the system by moving atoms according to a given force field and constraints.

- **Parameters:**
  - `int niter`: Number of iterations.
  - `double dt`: Time step for each iteration.
  - `double Fconv`: Convergence criterion for forces.
  - `double Flim`: Maximum allowed force magnitude.
  - `double damping`: Damping factor to control the movement of atoms.

##### `int run_t (int niter, double dt, double Fconv, double Flim, double damping=0.1)`

Template version of the optimization algorithm that allows for different behaviors based on template parameters.

- **Parameters:**
  - `int niter`: Number of iterations.
  - `double dt`: Time step for each iteration.
  - `double Fconv`: Convergence criterion for forces.
  - `double Flim`: Maximum allowed force magnitude.
  - `double damping`: Damping factor to control the movement of atoms.

##### `int run_omp (int niter, double dt, double Fconv, double Flim, double damping=0.1)`

Runs the optimization algorithm using OpenMP for parallelization, optimizing performance by minimizing dependencies between atoms.

- **Parameters:**
  - `int niter`: Number of iterations.
  - `double dt`: Time step for each iteration.
  - `double Fconv`: Convergence criterion for forces.
  - `double Flim`: Maximum allowed force magnitude.
  - `double damping`: Damping factor to control the movement of atoms.

##### `void printSizes ()`

Prints the sizes and dimensions of various data structures used in the UFF class, useful for debugging and performance tuning.