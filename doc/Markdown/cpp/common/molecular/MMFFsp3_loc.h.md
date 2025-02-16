# MMFFsp3_loc.h

## Purpose

`MMFFsp3_loc.h` is a header file that defines a class `MMFFsp3_loc`, which implements a localized version of the MMFF (Molecular Mechanics Force Field) for evaluating bonding interactions. This force field focuses on atoms with pi-orbitals, optimizing memory access and parallel evaluation by storing parameters and computing energy and forces locally.

## Includes

- `<omp.h>`: For OpenMP parallelization.
- `fastmath.h`: Provides fast mathematical operations.
- `Vec2.h` & `Vec3.h`: Define 2D and 3D vector classes for geometric calculations.
- `quaternion.h`: Handles quaternion operations, useful for representing orientations of pi-orbitals.
- `constants.h`: Contains physical constants used in the force field.
- `Forces.h`: Defines various physical interactions like Lennard-Jones (LJ) and Coulombic forces.
- `SMat3.h`: Symmetric 3x3 matrix class, useful for storing bond parameters.
- `molecular_utils.h`: Utility functions for molecular operations.
- `NBFF.h`: Non-Bonded Force Field base class.

## Types (Classes and Structs)

### Class `MMFFsp3_loc`

`MMFFsp3_loc` is a localized version of the MMFF force field, optimized for parallel evaluation by storing parameters per atom. It evaluates bonds, angles, and pi-orbitals separately to optimize memory access and computation speed.

**Inheritance**

- Inherits from `NBFF`, which provides basic non-bonded interactions like Lennard-Jones (LJ) forces.

#### Properties

- **`nDOFs`**: Number of degrees of freedom.
- **`Etot`**: Total energy of the system.
- **`DOFs`**: Degrees of freedom array for atoms and pi-orbitals.
- **`fDOFs`**: Forces corresponding to `DOFs`.
- **`doBonds`, `doNeighs`, `doPiPiI`, `doPiPiT`, `doPiSigma`, `doAngles`**: Flags controlling which interactions are computed.
- **`colDamp`**: Collision damping parameters for non-bonded interactions.
- **`cvf`**: Vector representing the sum of squared forces in each direction.
- **`bEachAngle`, `bTorsion`**: Flags to control angle and torsion energy evaluation.
- **`pipos`, `fpipos`**: Positions and forces on pi-orbitals.
- **`fneigh`, `fneighpi`**: Temporary storage for neighbor forces and pi-vector forces.
- **`bkneighs`**: Back-neighbors of each atom, used in periodic boundary conditions.
- **`angles`**: Array storing angles between bonds.
- **`invLvec`**: Inverse lattice vectors for periodic boundary conditions.
- **`bAngleCosHalf`**: Flag to use half-angle cosine evaluation.

#### Methods

- **`realloc`**: Reallocates memory for the class members based on new system dimensions.
- **`clone`**: Clones another `MMFFsp3_loc` object, optionally reallocating memory and deep copying parameters.
- **`dealloc`**: Deallocates all allocated memory.
- **`setLvec`**: Sets lattice vectors for periodic boundary conditions.
- **`optimalTimeStep`**: Calculates an optimal time step for the FIRE optimization algorithm.
- **`eval_atom`**: Evaluates energy and forces for a single atom considering its neighbors.
- **`eval_atom_t`**: Template version of `eval_atom`, allowing conditional evaluation based on parameters.
- **`eval_atom_opt`**: Optimized version of `eval_atom`.
- **`eval_atoms`**: Iteratively evaluates all atoms in the system.
- **`evalKineticEnergy`**: Calculates kinetic energy of the system.
- **`addjustAtomCapLenghs`**: Adjusts capping atom lengths to equilibrium distances from node atoms.
- **`addjustCapLenghs`**: Adjusts all capping atom lengths.
- **`initPi`**: Initializes pi-orbitals based on orthogonality of sigma bonds.
- **`relax_pi`**: Iteratively relaxes pi-orbitals until convergence.
- **`normalizePis`**: Normalizes pi-orbitals to unit length.
- **`constrainAtom`**: Constrains an atom's position using a spring force.
- **`cleanForce`**: Clears forces on all atoms and other degrees of freedom.
- **`cleanVelocity`**: Resets velocities of all atoms.
- **`assemble_atom`**: Assembles neighbor recoil forces into the global force array for a single atom.
- **`asseble_forces`**: Assembles forces for all atoms in the system.
- **`eval`**: Evaluates the total energy and forces of the system, optionally cleaning forces before evaluation.
- **`eval_check`**: Checks the consistency of the energy calculation by evaluating and checking for NaNs.
- **`run`**: Iteratively evaluates and moves atoms to minimize energy using gradient descent or MD.
- **`run_omp`**: Parallel version of `run`, utilizing OpenMP for parallelization.
- **`flipPis`**: Flips pi-orbitals to be in the same half-space as a reference vector.
- **`move_atom_GD`**: Updates atom positions using gradient descent.
- **`move_atom_Langevin`**: Updates atom positions using Langevin dynamics.
- **`move_Langevin`**: Wrapper for `move_atom_Langevin`.
- **`move_atom_MD`**: Updates atom positions using molecular dynamics (MD).
- **`move_atom_FIRE`**: Updates atom positions using the FIRE algorithm.
- **`move_atom_kvaziFIRE`**: Updates atom positions using a modified FIRE algorithm.
- **`move_GD`**: Wrapper for `move_atom_GD`.
- **`shiftBack`**: Shifts atom positions back to the first Brillouin zone in periodic boundary conditions.
- **`makeBackNeighs`**: Creates a list of back-neighbors for each atom, used in periodic boundary conditions.
- **`makeNeighCells`**: Sets neighbor cell indices based on lattice vectors or precomputed shifts.
- **`printSizes`**: Prints the sizes and dimensions of the system.
- **`printAtomParams`**: Prints parameters related to a specific atom.
- **`printNeighs`**: Prints neighbors for each atom.
- **`printBKneighs`**: Prints back-neighbors for each atom.
- **`print_pipos`**: Prints positions of pi-orbitals.
- **`print_apos`**: Prints atomic positions.
- **`print_pbc_shifts`**: Prints periodic boundary condition shifts.
- **`print_constrains`**: Prints constraints on atoms.
- **`printAngles`**: Prints angles between bonds for a specific atom.
- **`printTorsions`**: Prints torsion parameters.
- **`printAtomsConstrains`**: Prints constrained atoms.
- **`printDebug`**: Prints detailed information about the system state during debugging.
- **`checkNans`**: Checks for NaNs in various arrays and exits if found.
- **`setFromRef`**: Sets atom positions relative to a reference frame.
- **`rotateNodes`**: Rotates selected node atoms around a specified axis by a given angle.
- **`chargeToEpairs`**: Distributes charges from an atom to its capping electron pairs.
- **`measureCosPiPi`, `measureAnglePiPi`, `measureCosSigmaPi`, `measureAngleSigmaPi`**: Measure cosine and angles between pi-orbitals or sigma bonds with respect to a node atom.
- **`eval_atom_debug`**: Debug version of `eval_atom`, providing more detailed output.

## Implementation Details

The class is designed for efficient parallel evaluation by storing parameters per atom, avoiding synchronization issues in the global force array. It uses OpenMP for parallelization where possible and provides various methods to adjust and evaluate forces on atoms and their pi-orbitals. The localized approach ensures that memory access patterns are optimized for parallel execution, making it suitable for large-scale molecular dynamics simulations.