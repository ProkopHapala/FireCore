# NBFF.h

The `NBFF.h` file defines a class that implements non-bonded force-field interactions between particles using potentials such as Lenard-Jones, Morse, and Coulomb. It supports periodic boundary conditions (PBC) to handle long-range interactions correctly. This class is designed for efficient computation of non-covalent interactions in molecular mechanics simulations.

## Includes

- `fastmath.h`
- `Vec3.h`
- `quaternion.h`
- `Atoms.h`
- `Buckets.h`
- `Forces.h`
- `ForceField.h`
- `simd.h`

---

## Free Functions

- `fitAABB` - Fits an axis-aligned bounding box (AABB) around a set of particles. This function is used to determine the minimum and maximum coordinates of the particles, which can be useful for various spatial partitioning algorithms.
- `makePBCshifts_` - Generates all possible periodic boundary condition (PBC) shifts for a given lattice vector and number of PBC images in each direction. This is used to account for particles that are outside the primary simulation box but within the periodic boundary conditions.
- `evalPointCoulPBC` - Evaluates the Coulomb potential for a given point in space considering periodic boundary conditions (PBC). This function calculates the electrostatic interactions between all particles within the simulation box and its periodic images.
- `sampleCoulombPBC` - Samples the Coulomb potential over a set of points in space, considering periodic boundary conditions. This function is useful for generating a grid of Coulomb potentials or forces at specific positions within the simulation box.

---

## Types (Classes and Structs)

---

### Class `NBFF`

The Non-Bonded Force-Field (`NBFF`) class implements non-covalent interactions between particles using potentials such as Lenard-Jones, Morse, and Coulomb. It supports periodic boundary conditions to handle long-range interactions correctly. This class is a subclass of the `ForceField` class.

**Inheritance**

- **ForceField**: Inherits properties and methods related to force calculations from the base class.

#### Properties

- `nBBs`: Number of bounding boxes used for spatial partitioning.
- `BBs`: Array of bounding box objects (`Vec6d`) that define spatial regions for efficient collision detection.
- `pointBBs`: A `Buckets` object that manages the distribution of particles into bounding boxes.
- `drSR`: Short-range interaction cutoff radius minus a small value to avoid numerical issues.
- `ampSR`: Amplitude factor for short-range repulsion interactions.
- `alphaMorse`: Parameter controlling the shape of the Morse potential.
- `Rdamp`: Damping radius used in Coulomb potential calculations.
- `bPBC`: Boolean flag indicating whether periodic boundary conditions are enabled.
- `npbc`: Total number of PBC shifts.

#### Methods

- **`torq`**: Calculates the total torque on a molecule with respect to a given point. This method is useful for determining rotational forces in molecular dynamics simulations.
  
- **`bindShifts`**: Binds or reallocates an array of PBC shift vectors. This function ensures that the `shifts` array is properly initialized and allocated.

- **`makePBCshifts`**: Generates all possible PBC shift vectors for a given lattice vector and number of images in each direction. This method can be used to initialize the `shifts` array if it is not already allocated.

- **`evalPLQs`**: Evaluates PLQ (Pauli-London-Charge) parameters from REQ (Radius-Energy-Charge) parameters for faster evaluation using a factorized form, especially useful when using a grid-based approach. This method allocates memory and evaluates the PLQ parameters.

- **`makePLQs`**: Allocates memory and evaluates PLQ parameters as in `evalPLQs`.

- **`evalPLQd`**: Similar to `evalPLQs`, but uses double precision (`Quat4f`) for the PLQ parameters. This method alloculates memory and evaluates the PLQ parameters using double precision.

- **`makePLQd`**: Allocates memory and evaluates PLQ parameters using double precision, similar to `makePLQs`.

- **`updatePointBBs`**: Updates the bounding boxes based on the current positions of particles. This method ensures that the bounding boxes accurately represent the spatial distribution of particles.

- **`selectInBox`**: Selects particles within a given bounding box and returns their indices, positions, and non-bonding interaction parameters. This function is useful for efficient collision detection and force calculations.

- **`evalSortRange_BBs`**: Evaluates non-bonded interactions between particles in different bounding boxes. This method is useful for efficiently calculating long-range interactions using spatial partitioning techniques.

- **`evalLJQs`**: Evaluates the total energy of all non-bonded interactions using Lenard-Jones and Coulomb potentials without periodic boundary conditions. This function calculates both the potential energy and forces between particles.

- **`evalLJQs_ng4_atom`**: Evaluates non-bonded interactions for a single atom, excluding bonded atoms. This method is optimized to handle cases where each atom has at most four bonds.

- **`evalLJQs_ng4_omp`**: Similar to `evalLJQs_ng4_atom`, but uses OpenMP parallelization for better performance on multi-core systems.

- **`evalLJQs_ng4`**: Evaluates all non-bonded interactions excluding bonded atoms, using a single-threaded approach. This method is useful when the number of atoms is small or when parallelization overhead outweighs benefits.

- **`addMorseQH_PBC_omp`**: Adds Morse and Coulomb forces to a given atom considering periodic boundary conditions. This function uses OpenMP SIMD for efficient computation.

- **`getMorseQH_PBC_omp`**: Similar to `addMorseQH_PBC_omp`, but returns the total force without modifying it in place.

- **`getLJQs_PBC_omp`**: Evaluates non-bonded interactions using Lenard-Jones and Coulomb potentials, considering periodic boundary conditions. This function uses OpenMP SIMD for efficient computation.

- **`evalLJQs_PBC_atom_omp`**: Similar to `addMorseQH_PBC_omp`, but specifically evaluates the total energy of non-bonded interactions for a single atom using OpenMP SIMD.

- **`evalLJQs_PBC_simd`**: Evaluates all non-bonded interactions using Lenard-Jones and Coulomb potentials, considering periodic boundary conditions. This function uses SIMD instructions for efficient computation.

- **`evalLJQs_ng4_PBC_atom_omp`**: Similar to `evalLJQs_PBC_atom_omp`, but specifically evaluates the total energy of non-bonded interactions for a single atom with up to four bonds, using OpenMP SIMD.

- **`evalLJQs_ng4_PBC_simd`**: Evaluates all non-bonded interactions excluding bonded atoms, considering periodic boundary conditions. This function uses SIMD instructions for efficient computation.

- **`evalLJQs_atom_omp`**: Similar to `evalLJQs_PBC_atom_omp`, but evaluates the total energy of non-bonded interactions for a single atom using OpenMP SIMD without considering periodic boundary conditions.

- **`evalLJQs_simd`**: Evaluates all non-bonded interactions using Lenard-Jones and Coulomb potentials, considering periodic boundary conditions. This function uses SIMD instructions for efficient computation.

- **`evalLJQs_ng4_atom_omp`**: Similar to `evalLJQs_atom_omp`, but specifically evaluates the total energy of non-bonded interactions for a single atom with up to four bonds using OpenMP SIMD.

- **`evalLJQs_ng4_simd`**: Evaluates all non-bonded interactions excluding bonded atoms, considering periodic boundary conditions. This function uses SIMD instructions for efficient computation.

- **`evalLJQs_atom_avx`**: Similar to `evalLJQs_PBC_atom_omp`, but specifically evaluates the total energy of non-bonded interactions for a single atom using AVX SIMD instructions.

- **`evalCollisionDamp_atom_omp`**: Adds collision damping forces to a given atom. This function uses OpenMP SIMD for efficient computation.

- **`evalCollisionDamp_omp`**: Evaluates all collision damping forces considering periodic boundary conditions. This method is useful for simulating friction or other dissipative forces in the system.

- **`evalLJQs_ng4_PBC_atom`**: Similar to `evalLJQs_PBC_atom_omp`, but evaluates non-bonded interactions excluding bonded atoms, considering periodic boundary conditions.

- **`evalLJQs_ng4_PBC_omp`**: Evaluates all non-bonded interactions excluding bonded atoms, considering periodic boundary conditions. This method uses OpenMP parallelization for better performance on multi-core systems.

- **`evalLJQs_ng4_PBC`**: Similar to `evalLJQs_ng4_PBC_omp`, but evaluates the total energy of all non-bonded interactions excluding bonded atoms, considering periodic boundary conditions. This method is useful when the number of atoms is large and parallelization benefits outweigh overhead.

- **`evalLJQ`**: Evaluates the total energy of all non-bonded interactions between two `NBFF` objects using Lenard-Jones potentials without periodic boundary conditions.

- **`evalMorse`**: Evaluates the total energy of all Morse potential interactions between two `NBFF` objects, optionally including recoil forces. This method is useful for simulating systems with Morse potentials instead of Lennard-Jones potentials.

- **`evalMorsePLQ`**: Similar to `evalMorse`, but uses PLQ parameters for faster evaluation using a factorized form, especially when using a grid-based approach.

- **`evalR`**: Evaluates the total repulsion energy between two `NBFF` objects. This method is useful for simulating systems where only repulsive forces are considered.

- **`makePBCshifts`**: Similar to `makePBCshifts_`, but with additional parameters and a different return type, used internally by other methods.

- **`print_nonbonded`**: Prints detailed information about non-bonded interactions for debugging purposes. This function is useful for verifying the correctness of force calculations during development or testing.

- **`bindOrRealloc`**: Binds or reallocates memory for various arrays such as positions (`apos`), forces (`fapos`), non-bonding interaction parameters (`REQs`), and atomic types (`atypes`). This method ensures that the class can handle changes in the number of atoms without losing data.

- **`dealloc`**: Deallocates all resources associated with the `NBFF` object, ensuring proper memory management.