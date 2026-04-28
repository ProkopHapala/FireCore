# ForceField.h

## Includes

- "fastmath.h"
- "Vec3.h"
- "quaternion.h"
- "Atoms.h"
- "Forces.h"
- `<functional>`


## Free functions

### `bool checkLimits (int n, int m, const double* vals, const double* vmin, const double* vmax, const char* message, bool bPrint=true )`

#### Purpose
Check if the values in an array are within specified limits.

#### Parameters
- `int n`: Number of elements to check.
- `int m`: Number of dimensions for each element.
- `const double* vals`: Pointer to the array of values to be checked.
- `const double* vmin`: Pointer to the array of minimum allowed values.
- `const double* vmax`: Pointer to the array of maximum allowed values.
- `const char* message`: A string describing what is being checked (used for debugging).
- `bool bPrint`: Flag indicating whether to print messages if limits are exceeded.

#### Return Value
Returns true if any value exceeds its limit, otherwise false.


### `double cos_damp_lin (double c, double& cv, double D, double cmin, double cmax  )`

#### Purpose
Calculate a damping factor based on the cosine of an angle and update velocity accordingly.

#### Parameters
- `double c`: Cosine value of the angle.
- `double& cv`: Reference to the velocity component to be updated.
- `double D`: Damping coefficient.
- `double cmin`: Minimum allowed cosine value.
- `double cmax`: Maximum allowed cosine value.

#### Return Value
Returns a damping factor used in updating the velocity.


## Types (classes and structs)

### class `ForceField`

[Short description what is the purpose of the class, what its role in the bigger context]

**Inheritance**

- Atoms

#### Properties

- `Vec3d cvf`: A vector containing dot products of forces and velocities.
- `Vec3i nPBC`: Number of periodic boundary condition images in each direction.
- `bool bPBC`: Flag indicating whether periodic boundary conditions are used.
- `int npbc`: Total number of periodic images.
- `Vec3d* shifts`: Array of bond vector shifts for periodic boundary conditions.
- `bool bNonBonded`: Flag to indicate if non-bonded interactions should be calculated.
- `bool bNonBondNeighs`: Flag to indicate if non-bonded neighbors are used.
- `bool bSubtractBondNonBond`: Flag to indicate if bond energy is subtracted from non-bonded energy.
- `bool bSubtractAngleNonBond`: Flag to indicate if angle energy is subtracted from non-bonded energy.
- `bool bClampNonBonded`: Flag to indicate if non-bonded energy should be clamped to zero.
- `double FmaxNonBonded`: Maximum force allowed for non-bonded interactions.
- `std::function<double(int, const Vec3d&, Vec3d&)> atomForceFunc`: Function pointer or lambda function for calculating atomic forces.
- `double time_per_iter`: Time per iteration in the simulation.

#### Methods

- **`void setNonBondStrategy (int imode=0 )`**
  - Sets the strategy for non-bonded interactions based on the mode provided. Adjusts flags and settings accordingly to optimize calculations.

- **`Vec3d move_atom_MD (int i, const double dt, const double Flim, const double cdamp=0.9 )`**
  - Moves an atom using the damped leap-frog method of molecular dynamics.
  - Parameters:
    - `int i`: Index of the atom to be moved.
    - `const double dt`: Time step for the simulation.
    - `const double Flim`: Maximum force limit.
    - `const double cdamp`: Damping coefficient (default is 0.9).
  - Returns: A vector containing dot products of forces and velocities.

- **`Vec3d move_atom_Langevin (int i, const float dt, const double Flim, const double gamma_damp=0.1, double T=300 )`**
  - Moves an atom using the Langevin thermostat method.
  - Parameters:
    - `int i`: Index of the atom to be moved.
    - `const float dt`: Time step for the simulation.
    - `const double Flim`: Maximum force limit.
    - `const double gamma_damp`: Damping coefficient (default is 0.1).
    - `double T`: Temperature in Kelvin (default is 300 K).
  - Returns: A vector containing dot products of forces and velocities.

- **`double move_atom_FIRE (int i, const double dt, const double Flim, const double cv, const double cf )`**
  - Moves an atom using the FIRE algorithm.
  - Parameters:
    - `int i`: Index of the atom to be moved.
    - `const double dt`: Time step for the simulation.
    - `const double Flim`: Maximum force limit.
    - `const double cv`: Velocity damping coefficient.
    - `const double cf`: Force damping coefficient.
  - Returns: The squared norm of the force.

- **`double move_atom_kvaziFIRE (int i, const double dt, const double Flim )`**
  - Moves an atom using a modified FIRE algorithm.
  - Parameters:
    - `int i`: Index of the atom to be moved.
    - `const double dt`: Time step for the simulation.
    - `const double Flim`: Maximum force limit.
  - Returns: The squared norm of the force.

- **`double move_atom_GD (int i, double dt, double Flim)`**
  - Moves an atom using gradient descent.
  - Parameters:
    - `int i`: Index of the atom to be moved.
    - `double dt`: Time step for the simulation.
    - `double Flim`: Maximum force limit.
  - Returns: The squared norm of the force.

- **`double move_GD (float dt, double Flim=100.0 )`**
  - Moves all atoms using gradient descent.
  - Parameters:
    - `float dt`: Time step for the simulation.
    - `double Flim`: Maximum force limit (default is 100.0).
  - Returns: The total squared norm of forces.

- **`double move_MD (float dt, const double Flim=100.0, const double cdamp=0.9 )`**
  - Moves all atoms using molecular dynamics.
  - Parameters:
    - `float dt`: Time step for the simulation.
    - `const double Flim`: Maximum force limit (default is 100.0).
    - `const double cdamp`: Damping coefficient (default is 0.9).
  - Returns: The total squared norm of forces.

- **`void cleanForce ()`**
  - Resets the forces on all atoms to zero.
  
- **`void cleanVelocity ()`**
  - Resets the velocities of all atoms to zero.

- **`void copyForcesTo (Vec3d* fapos_)`**
  - Copies the current forces from `fapos` to a provided array.

- **`void copyPosTo (Vec3d* apos_)`**
  - Copies the current positions from `apos` to a provided array.


### class `CollisionDamping`

[Short description what is the purpose of the class, what its role in the bigger context]

#### Properties

- `public: bool bBond`: Flag indicating if collision damping for bonds should be used.
- `bool bAng`: Flag indicating if angular collisions are considered.
- `bool bNonB`: Flag indicating if non-bonded interactions are included in collision damping.
- `int nstep`: Number of steps to decay velocity by a factor of 1/e.
- `double medium`: Damping coefficient for the medium.
- `double bond`: Collision damping coefficient for bonds.
- `double ang`: Collision damping coefficient for angles.
- `double nonB`: Collision damping coefficient for non-bonded interactions.
- `double dRcut1`: Lower limit of the distance range for non-covalent collision damping.
- `double dRcut2`: Upper limit of the distance range for non-covalent collision damping.
- `double cdampB`: Collision damping coefficient for bonds.
- `double cdampAng`: Collision damping coefficient for angles.
- `double cdampNB`: Collision damping coefficient for non-bonded interactions.
- `double cos_vf_acc`: Cosine value used to determine if acceleration can be applied.
- `int nstep_acc`: Number of steps the velocity has been accelerated.
- `int nstep_acc_min`: Minimum number of steps required before acceleration is considered.


#### Methods

- **`bool canAccelerate ()`**
  - Checks if the current step count exceeds the minimum threshold for acceleration.

- **`double tryAccel ()`**
  - Returns a damping factor that allows for velocity acceleration if conditions are met.

- **`double update (double dt )`**
  - Updates collision damping coefficients based on time step and other parameters.
  
- **`void set (int nstep_, double medium_, double bond_, double ang_, double nonB_, double dRcut1_, double dRcut2_ )`**
  - Sets the properties of the `CollisionDamping` object with new values.

- **`void setup_accel (int nstep_acc_min_, double cos_vf_acc_ )`**
  - Configures the parameters for velocity acceleration based on minimum steps and cosine value.

- **`double update_acceleration (Vec3d cvf )`**
  - Updates the number of steps that can be accelerated based on the current state of velocities.