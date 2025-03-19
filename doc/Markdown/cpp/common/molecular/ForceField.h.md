# ForceField.h

This header file defines the `ForceField` class, which is responsible for managing and calculating forces in a molecular dynamics simulation. It includes various utility functions and classes to handle atomic interactions, force calculations, and other related operations.

## Includes

- `fastmath.h` - Provides mathematical utilities for fast computations.
- `Vec3.h` - Defines a 3-dimensional vector class used throughout the codebase.
- `quaternion.h` - Handles quaternion operations, which are useful in representing rotations and orientations.
- `Atoms.h` - Contains definitions related to atomic structures and properties.
- `Forces.h` - Manages force calculations for molecular dynamics simulations.
- `<functional>` - Allows the use of function objects (functors) as parameters.

---

## Free Functions

- `checkLimits` - Checks if a set of values are within specified limits. If any value is out of bounds, it prints a message indicating which value and its corresponding limits.

- `cos_damp_lin` **Purpose:** Calculates damping coefficients based on the cosine similarity between two vectors. This function is used in collision damping for velocity and force vectors.

---

## Types (Classes and Structs)

---

### class `ForceField`

This is base class for any force field. It provides a common interface for calculating forces and energies in a molecular system. It is designed to be extended by specific force field implementations.

**Inheritance**

- Inherits from `Atoms`.

#### Properties

- `cvf`: `Vec3d` - Represents the dot products between velocity and force vectors.
- `nPBC`: `Vec3i` - Stores the number of periodic boundary condition images along each axis.
- `bPBC`: `bool` - Indicates whether periodic boundary conditions are active.
- `npbc`: `int` - Total number of periodic boundary condition images.
- `bNonBonded`: `bool` - Determines if non-bonded interactions should be calculated.
- `bNonBondNeighs`: `bool` - Specifies the strategy for handling non-bonded neighbors.
- `bSubtractBondNonBond`: `bool` - Controls whether bond energy is subtracted from non-bonded energy.
- `bSubtractAngleNonBond`: `bool` - Controls whether angle energy is subtracted from non-bonded energy.
- `bClampNonBonded`: `bool` - Determines if non-bonded forces should be clamped to a maximum value.
- `FmaxNonBonded`: `double` - Maximum force magnitude for clamping non-bonded interactions.
- `time_per_iter`: `double` - Time taken per iteration of the simulation.
- `colDamp`: `CollisionDamping` - Manages collision damping parameters and behaviors.

#### Methods

- `setNonBondStrategy` - Sets the strategy for handling non-bonded interactions based on an integer mode. Adjusts various flags accordingly.
- `move_atom_MD` - Updates atomic positions using a damped leap-frog method, which is part of molecular dynamics simulations.
- `move_atom_Langevin` - Updates atomic positions using the Langevin thermostat, simulating Brownian motion with friction and random forces.
- `move_atom_FIRE` - Updates atomic positions using the Fast Inertial Relaxation Engine (FIRE) algorithm for efficient relaxation.
- `move_atom_kvaziFIRE` - Updates atomic positions using a modified FIRE algorithm that includes damping based on cosine similarity.
- `move_atom_GD` - Updates atomic positions using gradient descent, minimizing potential energy by adjusting forces.
- `move_GD` - Iterates over all atoms and applies the gradient descent method to update their positions.
- `move_MD` - Iterates over all atoms and applies the molecular dynamics (MD) method with damping to update their positions.
- `cleanForce` - Resets force vectors to zero, ensuring no residual forces are present in the system.
- `cleanVelocity` - Resets velocity vectors to zero, setting all atomic velocities to rest.
- `copyForcesTo` - Copies force vectors from one array to another, useful for transferring data between different parts of the simulation.
- `copyPosTo` - Copies position vectors from one array to another, facilitating the transfer of atomic positions.

---

### class `CollisionDamping`

This class manages collision damping parameters and behaviors. It is used to control the damping of atomic velocities during molecular dynamics simulations.

#### Properties

- `bool`: Public - Flags indicating whether collision damping should be applied for bonds, angles, and non-bonded interactions.
- `nstep`: `int` - Number of steps required to decay velocity to 1/e of its initial value during collision damping.
- `medium`: `double` - Damping coefficient for the medium effect in collision damping.
- `bond`: `double` - Collision damping coefficient for bond interactions.
- `ang`: `double` - Collision damping coefficient for angle interactions.
- `nonB`: `double` - Collision damping coefficient for non-bonded interactions.
- `dRcut1`: `double` - Defines the start of the non-covalent collision damping interaction range.
- `dRcut2`: `double` - Defines the end of the non-covalent collision damping interaction range.
- `cdampB`: `double` - Collision damping coefficient for bond interactions after updating.
- `cdampAng`: `double` - Collision damping coefficient for angle interactions after updating.
- `cdampNB`: `double` - Collision damping coefficient for non-bonded interactions after updating.
- `cos_vf_acc`: `double` - Cosine value used to determine if acceleration can occur in the collision damping process.
- `nstep_acc`: `int` - Number of steps that have been counted towards potential acceleration.
- `nstep_acc_min`: `int` - Minimum number of steps required before considering acceleration.

#### Methods

- `canAccelerate` - Checks if the current step count allows for acceleration in collision damping.
- `tryAccel` - Returns a value indicating whether to attempt acceleration based on the current step count and cosine similarity.
- `update` - Updates the collision damping coefficients based on the specified parameters and time step.
- `set` - Sets various properties of the collision damping object, including flags and damping coefficients.
- `setup_accel` - Configures the minimum number of steps required for acceleration and the cosine value used in this process.
- `update_acceleration` - Updates the acceleration state based on the current velocity-force vector similarity.