# MMFFsp3_loc.h

This file defines the `MMFFsp3_loc` class, which implements a localized version of the MMFF force field for classical molecular mechanics. The class focuses on optimizing memory access patterns by storing all parameters per atom, allowing efficient parallel evaluation of bonding interactions.

## Includes

- `<omp.h>`: OpenMP library for parallelization.
- `fastmath.h`: Provides fast mathematical operations.
- `Vec2.h` & `Vec3.h`: Classes for 2D and 3D vectors.
- `quaternion.h`: Defines quaternion structures used in rotations.
- `constants.h`: Contains physical constants.
- `Forces.h`: Includes physical interactions like Lennard-Jones potential.
- `SMat3.h`: Symmetric matrix class.
- `molecular_utils.h`: Molecular structure utilities.
- `NBFF.h`: Non-Bonded Force Field interface.

---

## Types (Classes and Structs)

---

### Class `MMFFsp3_loc`

This class implements a localized version of the MMFF force field for classical molecular mechanics. It focuses on optimizing memory access patterns by storing all parameters per atom, allowing efficient parallel evaluation of bonding interactions.

**Inheritance**

- **NBFF**: Inherits properties and methods from the Non-Bonded Force Field class to handle non-bonded interactions like Lennard-Jones potential.

#### Properties

- `static`:`public:`
  - `nDOFs`: Number of degrees of freedom.
  - `Etot`: Total energy of the system.
  - `DOFs`: Degrees of freedom array for atoms and pi-orbitals.
  - `fDOFs`: Forces corresponding to DOFs.
  - `doBonds`, `doNeighs`, `doPiPiI`, `doPiPiT`, `doPiSigma`, `doAngles`: Flags controlling which interactions are computed.
  - `colDamp`: Collision damping parameters for non-bonded interactions.
  - `cvf`: Vector representing the sum of squared forces in each direction.
  - `bEachAngle`: Boolean flag indicating whether to compute angle energy separately for each angle.
  - `bTorsion`: Boolean flag controlling torsion energy computation.

- **Dynamic Variables**:
  - `pipos`, `fpipos`: Pointers to vectors storing pi-orbitals and their forces, respectively.
  - `fneigh`, `fneighpi`: Temporary storage for neighbor bond and pi-vector forces before assembly into global force arrays.
  - `bkneighs`: Back-neighbors array used in periodic boundary conditions.
  - `angles`: Array of angles between bonds.
  - `invLvec`: Inverse lattice vectors for handling periodic boundary conditions.

#### Methods

##### Memory Management

- `realloc`: Allocates memory for class members.
- `clone`: Clones another `MMFFsp3_loc` object.
- `dealloc`: Deallocates dynamically allocated memory.

##### System Configuration

- `setLvec`: Sets the lattice vectors and computes their inverse.

##### Energy and Force Evaluation

- `eval_atom`: Evaluates energy and forces for a single atom.
- `eval_atoms`: Iterates over all atoms to compute total energy and forces.
- `evalKineticEnergy`: Computes the kinetic energy of the system.
- `eval`: Evaluates the total energy and forces for the system.
- `eval_check`: Performs checks on the system without evaluation cleanup.

##### Pi Orbital Handling

- `initPi`: Initializes pi orbitals.
- `relax_pi`: Refines pi orbital directions.
- `normalizePis`: Normalizes pi orbitals.
- `flipPis`: Rotates pi orbitals to align with a reference vector.

##### Force Assembly

- `assemble_atom`: Assembles neighbor recoil forces into the global force array.
- `asseble_forces`: Iterates over all atoms to assemble forces.

##### Atom Position Manipulation

- `addjustAtomCapLenghs`: Adjusts capping atom positions.
- `addjustCapLenghs`: Adjusts all capping lengths.
- `constrainAtom`: Constrains an atom's position.
- `cleanForce`: Clears all forces.
- `cleanVelocity`: Clears all velocities.
- `shiftBack`: Shifts atomic coordinates back to their primary unit cell.
- `setFromRef`: Sets positions from reference.
- `rotateNodes`: Rotate selected atoms and their caps around a specified axis by a given angle.

##### Movement of Atoms

- `move_atom_GD`: Update atom positions using gradient descent.
- `move_atom_Langevin`: Update atom positions using langevin dynamics.
- `move_Langevin`: Update atom positions using langevin dynamics.
- `move_atom_MD`: Update atom positions using molecular dynamics.
- `move_atom_FIRE`: Update atom positions using FIRE algorithm.
- `move_atom_kvaziFIRE`: Update atom positions using kvazi FIRE algorithm.

##### Periodic Boundary Conditions

- `makeBackNeighs`: Constructs lists of back-neighbors.
- `makeNeighCells`: Constructs lists of neighbor cell indices.

##### Running the Simulation

- `run`: Iteratively evaluates the force field and moves atoms.
- `run_omp`: Iteratively evaluates the force field and moves atoms with OpenMP.
- `optimalTimeStep`: Calculates an optimal time step for optimization.

##### Print Utility Functions

- `printSizes`: Prints sizes of data structures.
- `printAtomParams`: Prints atom parameters.
- `printNeighs`: Prints neighbor information.
- `printBKneighs`: Prints back-neighbor information.
- `print_pipos`: Prints pi orbital positions.
- `print_apos`: Prints atom positions.
- `print_pbc_shifts`: Prints PBC shift vectors.
- `print_constrains`: Prints constraint information.
- `printAngles`: Prints angles between bonds.
- `printAngles`: Prints angles between bonds for a specific atom.

##### Error Checking

- `checkNans`: Checks for NaN values.

##### Charge Distribution

- `chargeToEpairs`: Redistribute charges from atoms to electron pairs.

##### Measurement Functions

- `measureCosPiPi`: Measure the cosine between pi-orbitals.
- `measureAnglePiPi`: Measure the angle between pi-orbitals.
- `measureCosSigmaPi`: Measure the cosine between a sigma bond and a pi-orbital.
- `measureAngleSigmaPi`: Measure the angle between a sigma bond and a pi-orbital.

##### Debugging

- `eval_atom_debug`: A debug version of `eval_atom`.