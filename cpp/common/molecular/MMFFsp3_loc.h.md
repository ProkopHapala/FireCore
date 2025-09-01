# MMFFsp3_loc.h

## Includes

- `<omp.h>`: OpenMP header for parallelization.
- `"fastmath.h"`: Provides fast mathematical operations.
- `"Vec2.h"`: Defines 2D vector structures.
- `"Vec3.h"`: Defines 3D vector structures.
- `"quaternion.h"`: Implements quaternion mathematics.
- `"constants.h"`: Contains physical constants used in the force field.
- `"Forces.h"`: Provides various physical interaction forces.
- `"SMat3.h"`: Symmetric matrix operations for internal calculations.
- `"molecular_utils.h"`: Utility functions related to molecular structures and interactions.
- `"NBFF.h"`: Non-Bonded Force Field class, which `MMFFsp3_loc` inherits from.

## Types (classes and structs)

### class `MMFFsp3_loc`

**Purpose**

The `MMFFsp3_loc` class is a localized version of the MMFF force field designed for efficient evaluation of bonding interactions in molecular systems. It optimizes memory access patterns by storing parameters and forces per atom, which allows for parallelized evaluation without synchronization overhead.

**Inheritance**

- Inherits from `NBFF`, providing additional functionalities specific to localized atomic representations.

#### Properties

- **`public: static int nDOFs`**: Number of degrees of freedom in the system.
- **`double Etot`**: Total energy of the system.
- **`double * DOFs`**: Array holding all degrees of freedom (positions and velocities).
- **`double * fDOFs`**: Array for forces corresponding to each degree of freedom.
- **`bool doBonds`**: Flag indicating whether bond interactions should be computed.
- **`bool doNeighs`**: Flag indicating whether neighbor interactions should be computed.
- **`bool doPiPiI`**: Flag indicating whether pi-pi interaction energy should be computed.
- **`bool doPiPiT`**: Flag indicating whether pi-pi torsion energy should be computed.
- **`bool doPiSigma`**: Flag indicating whether pi-sigma interactions should be computed.
- **`bool doAngles`**: Flag indicating whether angle interactions should be computed.
- **`CollisionDamping colDamp`**: Object handling collision damping parameters for non-bonded interactions.
- **`Vec3d cvf`**: Vector representing the sum of squared forces in different directions.
- **`bool bEachAngle`**: Flag indicating if each angle is evaluated separately with its own parameters.
- **`bool bTorsion`**: Flag indicating whether torsional (proper dihedral) interactions should be computed.

#### Methods

##### `void realloc(int nnode_, int ncap_, int ntors_=0)`

Reallocates memory for the system based on new dimensions (`nnode`, `ncap`, and optionally `ntors`). Ensures that all necessary arrays are resized appropriately to accommodate changes in the number of atoms, capping atoms, and torsions.

##### `void clone(MMFFsp3_loc& from, bool bRealloc, bool bREQsDeep=true)`

Clones another instance (`from`) into this object. This includes copying parameters and forces while optionally deep-coping non-bonded interaction parameters (`REQs`).

##### `void dealloc()`

Deallocates all memory allocated for the system, ensuring that resources are freed when no longer needed.

##### `void setLvec(const Mat3d& lvec_)`

Sets the lattice vectors of the system. This is crucial for periodic boundary conditions and ensures correct distance calculations in a crystal-like structure.

##### `double optimalTimeStep(double m=1.0)`

Calculates an optimal time step for the dynamics algorithm based on mass (`m`) and stiffness parameters, ensuring stable integration over time steps.

##### `double eval_atom(const int ia)`

Evaluates the energy and forces for a single atom (`ia`). This method handles bond stretching, angle bending, pi-pi interactions, and other local interactions. Forces are accumulated in temporary arrays to be assembled later.

##### `double eval_atom_t(const int ia)`

Template version of `eval_atom` that allows for conditional evaluation based on various flags passed as template parameters. Useful for optimizing performance by enabling or disabling specific interaction types at compile time.

##### `double eval_atom_opt(const int ia)`

Optimized version of `eval_atom` with reduced branching and simplified logic, enhancing execution speed during runtime.

##### `double eval_atoms(bool bDebug=false, bool bPrint=false)`

Evaluates energy for all atoms in the system. If debugging is enabled, it evaluates each atom individually; otherwise, it uses optimized batch evaluation.

##### `double evalKineticEnergy()`

Calculates and returns the kinetic energy of the system based on atomic velocities. This function provides insights into the thermal state of the molecular system.

##### `void addjustAtomCapLenghs(int ia)`

Adjusts the positions of capping atoms to maintain equilibrium distances from node atoms, ensuring accurate representation of molecular structures.

##### `void addjustCapLenghs()`

Iterates over all node atoms and adjusts their capping atoms' positions using `addjustAtomCapLenghs`.

##### `void initPi(Vec3d* pbc_shifts, double Kmin=0.0001, double r2min=1e-4, bool bCheck=true)`

Initializes the directions of pi orbitals on node atoms based on orthogonality and strength of sigma bonds. This step is crucial for setting up the molecular structure correctly.

##### `int relax_pi(int niter, double dt, double Fconv, double Flim=1000.0)`

Iteratively refines the positions of pi orbitals to minimize energy by adjusting them according to specified convergence criteria and time steps.

##### `void normalizePis()`

Normalizes all pi orbital vectors to unit length, ensuring they are correctly oriented in space.

##### `void constrainAtom(int ia, double Kfix=1.0)`

Constrains the position of an atom (`ia`) to a fixed point using Hooke's law with spring constant `Kfix`.

##### `void cleanForce()`

Clears all forces from the system, preparing it for new evaluations.

##### `void cleanVelocity()`

Sets velocities of all atoms to zero, resetting their motion state.

##### `void assemble_atom(int ia)`

Assembles the recoil forces on an atom (`ia`) by summing up contributions from its neighbors. This step ensures that the total force acting on each atom is correctly computed before moving it according to the dynamics algorithm.

##### `void asseble_forces()`

Iterates over all atoms and calls `assemble_atom` for each, ensuring that forces are assembled correctly across the entire system.

##### `double eval(bool bClean=true, bool bCheck=true)`

Evaluates the total energy of the system by calling `eval_atoms`, then assembles forces using `asseble_forces`. This method is used to minimize the potential energy and move atoms accordingly.

##### `double eval_check()`

A debug version of `eval` that checks for NaN values in various arrays, ensuring numerical stability during evaluation.

##### `int run(int niter, double dt, double Fconv, double Flim, double damping=0.1)`

Iteratively evaluates the force field and moves atoms to minimize energy using a gradient descent algorithm with optional damping. This method is used for molecular dynamics simulations.

##### `int run_omp(int niter, double dt, double Fconv, double Flim, double damping=0.1)`

Parallelized version of `run` that utilizes OpenMP for parallel evaluation and force assembly, enhancing performance on multi-core systems.

##### `void flipPis(Vec3d ax)`

Flips the direction of pi orbitals to align with a given reference vector (`ax`), ensuring they are in the same half-space.

##### `double move_atom_GD(int i, float dt, double Flim)`

Updates atom positions using gradient descent. This method minimizes potential energy by adjusting velocities and positions iteratively.

##### `Vec3d move_atom_Langevin(int i, const float dt, const double Flim, const double gamma_damp=0.1, double T=300)`

Updates atom positions using the Langevin dynamics algorithm, which includes frictional damping and random forces to simulate thermal motion.

##### `Vec3d move_Langevin(const float dt, const double Flim, const double gamma_damp=0.1, double T=300)`

Applies the Langevin dynamics algorithm to all atoms in the system simultaneously, returning a summary of force changes.

##### `Vec3d move_atom_MD(int i, const float dt, const double Flim, const double cdamp=0.9)`

Updates atom positions using molecular dynamics with damping, ensuring stable and accurate motion updates.

##### `double move_atom_FIRE(int i, float dt, double Flim, double cv, double cf)`

Updates atom positions using the Fast Inertial Relaxation Engine (FIRE), a highly efficient algorithm for minimizing potential energy in molecular systems.

##### `double move_atom_kvaziFIRE(int i, float dt, double Flim)`

A modified version of FIRE that includes additional damping factors to improve convergence and stability during minimization processes.

##### `double move_GD(float dt, double Flim=100.0)`

Updates atom positions using gradient descent for all atoms in the system simultaneously, returning a summary of force changes.

##### `Vec3d shiftBack(bool bPBC=false)`

Shifts atomic coordinates back to the primary unit cell if periodic boundary conditions are active, ensuring that all distances and angles are correctly calculated within the simulation box.

##### `void makeBackNeighs(bool bCapNeighs=true)`

Generates a list of back-neighbors for each atom, which can be used in subsequent calculations. This step is particularly useful when dealing with complex molecular structures or large systems.

##### `void makeNeighCells(const Vec3i nPBC_)`

Sets the cell indices for neighbors based on periodic boundary conditions (`nPBC_`). This ensures that all interactions are correctly calculated across different unit cells of a crystal-like structure.

##### `void printSizes()`

Prints out the sizes and dimensions of various arrays used in the system, providing insights into memory usage and performance characteristics.

##### `void printAtomParams(int ia)`

Prints detailed information about an atom (`ia`), including its type, neighbor indices, bond lengths, stiffness constants, and other relevant parameters. This function is useful for debugging and validation purposes.

##### `void printNeighs(int ia)`

Prints the list of neighbors for a given atom (`ia`), helping to visualize connectivity within the molecular structure.

##### `void printBKneighs(int ia)`

Prints back-neighbors for an atom (`ia`), which can be useful in understanding complex interactions or constraints within the system.

##### `void printAtomParams()`

Prints detailed information about all atoms in the system, providing a comprehensive overview of their properties and interactions.

##### `void printNeighs()`

Prints the list of neighbors for all atoms in the system, offering insights into connectivity patterns across the entire molecular structure.

##### `void printBKneighs()`

Prints back-neighbors for all atoms in the system, aiding in understanding complex interaction networks or constraints within the molecular framework.

##### `void print_pipos()`

Prints out the positions of pi orbitals for node atoms, which are crucial for representing electronic structures accurately.

##### `void print_apos()`

Prints out the atomic positions for all atoms in the system, providing a clear visualization of the molecular structure.

##### `void print_pbc_shifts()`

Prints out periodic boundary condition shifts used during calculations, ensuring that interactions across different unit cells are correctly accounted for.

##### `void print_constrains()`

Prints out constraints applied to certain atoms, which can be useful in understanding fixed or restricted motion within the system.

##### `void printAngles(int ia)`

Prints detailed information about angles involving an atom (`ia`), including types of atoms involved and parameters used in angle calculations. This function is helpful for debugging and validation purposes.

##### `void printAngles()`

Prints detailed information about all angles in the system, offering a comprehensive overview of angular interactions within the molecular structure.

##### `void printTorsions()`

Prints out torsional (dihedral) parameters for all atoms in the system, providing insights into rotational constraints and flexibility within the molecular framework.

##### `void printAtomsConstrains(bool bWithOff=false)`

Prints out constraints applied to atoms, including those with zero or non-zero spring constants. This function is useful for debugging and validation purposes.

##### `void printDebug(bool bNg=true, bool bPi=true, bool bA=true)`

Prints detailed information about various aspects of the system (neighbors, pi orbitals, atomic positions), providing a comprehensive overview suitable for debugging and analysis.

##### `bool checkNans(bool bExit=true, bool bNg=true, bool bPi=true, bool bA=true)`

Checks for NaN values in various arrays used during calculations. If any are found, it prints out relevant information and optionally exits the program to prevent further errors.

##### `void setFromRef(Vec3d* aref, Vec3d* piref, Vec3d dp=Vec3dZero, Mat3d rot=Mat3dIdentity)`

Sets atomic positions based on a reference configuration (`aref`), allowing for easy alignment and comparison of molecular structures.

##### `void rotateNodes(int n, int* sel, Vec3d p0, Vec3d ax, double phi)`

Rotates selected node atoms around an axis (`ax`) by a specified angle (`phi`). This function is useful for aligning or transforming specific parts of the molecular structure.

##### `void chargeToEpairs(Quat4d* REQs, int* atypes, double cQ=-0.2, int etyp=-1)`

Redistributes charges from atoms to electron pairs based on specified parameters (`REQs`, `atypes`). This function is useful for modeling electronic properties of molecules accurately.

##### `void chargeToEpairs(double cQ=-0.2, int etyp=-1)`

A simplified version of `chargeToEpairs` that uses default values for redistribution charges and electron pair types.

##### `double measureCosPiPi(int ia, int ib, bool bRenorm=true)`

Measures the cosine of the angle between two pi orbitals (`ia`, `ib`). This function is useful for assessing electronic interactions within molecular structures.

##### `double measureAnglePiPi(int ia, int ib, bool bRenorm=true)`

Calculates the angle between two pi orbitals (`ia`, `ib`) in radians. This function provides insights into the spatial orientation of electron pairs within molecules.

##### `double measureCosSigmaPi(int ipi, int ia, int ib)`

Measures the cosine of the angle between a sigma bond and a pi orbital. This function is useful for understanding hybridization and electronic interactions at atomic scales.

##### `double measureAngleSigmaPi(int ipi, int ia, int ib)`

Calculates the angle between a sigma bond and a pi orbital in radians. This function offers detailed information about spatial arrangements within molecular structures.

##### `double eval_atom_debug(const int ia, bool bPrint=true)`

A debug version of `eval_atom` that prints out intermediate results for atom (`ia`). This method is useful during development or debugging sessions to ensure correct calculations and interactions are occurring as expected.