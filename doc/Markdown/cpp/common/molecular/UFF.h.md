# UFF.h

The `UFF.h` file implements a universal force field (UFF) for molecular systems based on the article "J. Am. Chem. Soc. 1992, 114, 25, 10024–10035". It provides functionalities to evaluate and optimize interatomic forces in molecules using bonds, angles, proper dihedrals (torsions), and improper dihedrals (plane inversions). The class is designed for efficient evaluation of the force field both serially and in parallel using OpenMP.

## Includes

- `<omp.h>`: OpenMP library for parallel programming.
- `fastmath.h`: Provides fast mathematical operations optimized for performance.
- `Vec2.h`: Defines a 2D vector class.
- `Vec3.h`: Defines a 3D vector class used extensively in the UFF calculations.
- `quaternion.h`: Implements quaternion mathematics, useful for rotations and orientations.
- `Forces.h`: Contains various physical interaction forces like Lennard-Jones (LJ) and Coulomb interactions.
- `SMat3.h`: Represents symmetric matrices which are used to store certain parameters or results.
- `molecular_utils.h`: Provides utility functions related to molecular structures and calculations.
- `NBFF.h`: Non-Bonded Force Field class, from which UFF inherits properties and methods.
- `GOpt.h`: Global optimization utilities for the force field evaluation.
- `Buckets.h`: Data structure used for efficient mapping of atom interactions.

---

## Free Functions

- `checkVec3Match`: Checks if two 3D vectors are approximately equal.
- `checkVec3Matches`: Checks multiple pairs of 3D vectors for approximate equality.

---

## Types (Classes and Structs)

---

### Class `UFF`

The `UFF` class implements the Universal Force Field (UFF) for molecular simulations. It extends the `NBFF` class to include specific functionalities and properties required for UFF calculations.

#### Inheritance

Inherits from `NBFF` to leverage non-bonded force field functionalities.

#### Properties

##### Energy Components
- `double Etot`: Total energy of the system.
- `double Eb`: Energy due to bonds.
- `double Ea`: Energy due to angles.
- `double Ed`: Energy due to dihedrals (proper torsions).
- `double Ei`: Energy due to improper dihedrals (plane inversions).

##### System Configuration
- `int nbonds`: Number of bond interactions in the system.
- `int i0dih`: Index offset for proper dihedral forces.
- `double SubNBTorsionFactor`: Factor to subtract torsion energy from non-bonded energy if greater than zero.

##### Data Structures
- `Vec3d* fint`: Temporary storage of forces on atoms before assembling them into final force arrays.
- `Vec3d* fbon`: Forces due to bonds (before assembly).
- `Vec3d* fang`: Forces due to angles (before assembly).
- `Vec3d* fdih`: Forces due to proper dihedrals (torsions) (before assembly).
- `Vec3d* finv`: Forces due to improper dihedrals (plane inversions) (before assembly).
- `Buckets a2f`: Data structure for fast force assembling.

##### Periodic Boundary Conditions
- `Mat3d invLvec`: Inverse lattice vectors used for periodic boundary conditions.

##### Optimization
- `GOpt* go`: Pointer to global optimization utilities.

#### Methods

##### Memory Management
- `realloc`: Reallocates memory and initializes properties for the UFF object.
- `dealloc`: Deallocates memory used by the UFF object.

##### System Configuration
- `setLvec`: Sets lattice vectors for periodic boundary conditions and computes their inverse.
- `makeNeighCells`: Creates a list of neighbor cell indices for periodic boundary conditions.

##### Neighbor Lists
- `mapAtomInteractions`: Maps atom interactions to force pieces using a buckets structure.
- `makeNeighBs`: Initializes neighbor bond indices for each atom.
- `bakeDihedralNeighs`: Bakes dihedral neighbor information into the system.
- `bakeAngleNeighs`: Bakes angle neighbor information into the system.
- `bakeInversionNeighs`: Bakes improper dihedral (inversion) neighbor information into the system.

##### Force Assembly
- `cleanForce`: Resets all forces and other DOFs to zero.
- `assembleForcesDebug`: Debug function to print out forces assembled from different interactions.
- `assembleForces`: Assembles total forces on atoms.
- `printForcePieces`: Prints out the force pieces for debugging purposes.
- `assembleAtomForce`: Assembles forces on a single atom.
- `assembleAtomsForces`: Iterates over all atoms and assembles their forces.

##### Energy Evaluation
- `evalAtomBonds`: Evaluates bond interactions for a single atom.
- `evalBonds`: Evaluates all bond interactions in the system.
- `evalAngle_Prokop`: Evaluates angle interactions using Prokop's method.
- `evalAngle_Paolo`: Evaluates angle interactions using Paolo's method.
- `evalAngles`: Evaluates all angle interactions in the system.
- `evalDihedral_Prokop`: Evaluates proper dihedral (torsion) interactions using Prokop's method.
- `evalDihedral_Prokop_Old`: Evaluates proper dihedral (torsion) interactions using an older version of Prokop's method.
- `evalDihedral_Paolo`: Evaluates proper dihedral (torsion) interactions using Paolo's method.
- `evalDihedrals`: Evaluates all proper dihedral (torsion) interactions in the system.
- `evalInversion_Prokop`: Evaluates improper dihedral (plane inversion) interactions using Prokop's method.
- `evalInversion_Paolo`: Evaluates improper dihedral (plane inversion) interactions using Paolo's method.
- `evalInversions`: Evaluates all improper dihedral (plane inversion) interactions in the system.
- `eval`: Evaluates the total energy of the UFF force field.
- `eval_omp_old`: Evaluates the total energy using OpenMP with an older implementation.
- `eval_omp`: Evaluates the total energy using OpenMP with a more recent implementation.

##### Simulation Control
- `run`: Runs molecular dynamics or optimization to minimize the system's potential energy.
- `run_t`: Template version of the `run` method with different parameters.
- `run_omp`: Runs molecular dynamics or optimization using OpenMP.

##### Debugging
- `printSizes`: Prints out sizes and dimensions of various data structures.