# GridFF.h

`GridFF.h` defines the `GridFF` class, which provides functionality for evaluating force fields on a grid structure. This class is used to compute electrostatic potentials and forces in molecular systems using various interpolation methods such as tricubic, Hermite, and B-spline interpolations. It also supports periodic boundary conditions (PBC) and can handle different types of force field models.

## Includes

- `fastmath.h`
- `Vec2.h`
- `Vec3.h`
- `quaternion.h`
- `Grid.h`
- `Forces.h`
- `MMFFparams.h`
- `Multipoles.h`
- `InterpolateTricubic.h`
- `InterpolateTrilinear.h`
- `Bspline.h`
- `Bspline_fit.h`
- `Bspline_fit_2D.h`
- `Bspline_fit_3D.h`
- `VecN.h`
- `NBFF.h`
- `EwaldGrid.h`
- `IO_utils.h`

---

## Free Functions

- `sum`: Calculates the sum of a vector using an iterator.
- `autoNPBC`: Automatically calculates the number of periodic boundary conditions (nPBC) based on minimum length.
- `evalDipole`: Evaluates the dipole moment around a center point that minimizes the size of the dipole.

---

## Types (Classes and Structs)

---

### Class `GridFF`

#### Purpose

`GridFF` is a class that extends from `NBFF`, providing methods to evaluate force fields on a grid. It supports various interpolation techniques like tricubic, Hermite, and B-spline interpolations.

**Inheritance**

- Inherits from `NBFF`.

#### Properties

##### Grid Definition
- `GridShape grid`: Represents the shape of the grid.
- `Vec3i gridN`: Grid dimensions with additional padding.

##### Force Field Data
- `Quat4f* FFPaul`: Pointer to Pauli force field data (float).
- `Quat4f* FFLond`: Pointer to London force field data (float).
- `Quat4f* FFelec`: Pointer to electric force field data (float).
- `Quat4d* FFPaul_d`: Pointer to Pauli force field data (double).
- `Quat4d* FFLond_d`: Pointer to London force field data (double).
- `Quat4d* FFelec_d`: Pointer to electric force field data (double).
- `Quat4d FEscale`: Scaling factor for the force fields.
- `Quat4f* VPLQH`: Pointer to interpolated force field data.
- `double* HHermite_d`: Pointer to Hermite interpolation coefficients (double).
- `double* Bspline_Pauli`: Pointer to B-spline Pauli force field data (double).
- `double* Bspline_London`: Pointer to B-spline London force field data (double).
- `double* Bspline_Coulomb`: Pointer to B-spline Coulomb force field data (double).
- `Vec3d* Bspline_PLQ`: Array of B-spline PLQ values (double).
- `Quat4f* Bspline_PLQf`: Array of B-spline PLQ values as `Quat4f` (float).

##### Ewald Summation
- `bool bUseEwald`: Flag indicating whether Ewald summation is used.
- `EwaldGrid* ewald`: Pointer to the Ewald grid object for electrostatic calculations.

##### Debugging
- `double* V_debug`: Debugging variable for storing intermediate values.
- `int iDBG`: Debugging flag.

##### Interpolation
- `GridFFmod mode`: Current interpolation method used (`Direct`, `LinearFloat`, etc.).
- `int perVoxel`: Number of data points per voxel.

#### Methods

##### System Setup
- `bindSystem`: Binds a system to the grid by setting up atomic positions and types.
- `allocateAtoms`: Allocates memory for atomic positions and quaternion moments.
- `setAtoms`: Sets up atomic positions and quaternion moments for evaluation.
- `setAtomsSymetrized`: Sets up symmetrized atomic positions by duplicating atoms within specified symmetry bounds.
- `loadCell`: Loads the cell dimensions from a file.
- `init`: Initializes the grid and sets up periodic boundary conditions.
- `initGridFF`: Initializes the grid and sets up atomic positions, quaternion moments, and other parameters for evaluation.

##### Memory Management
- `allocateFFs`: Allocates memory for force field data based on the grid size and whether double precision is used.
- `clear`: Deallocates all allocated memory related to force fields.

##### Grid Evaluation
- `evalAtPoints_REQ`: Evaluates the force field at specified points using the requested interpolation method.
- `evalAtPoints`: Evaluates the force field at specified points using the requested interpolation method and grid dimensions.
- `evalAtPoints_Split`: Evaluates the force field at specified points using separate methods for Coulomb and Morse interactions.
- `evalBsplineRef`: Evaluates the B-spline reference for the force field.
- `evalMorsePBC`: Evaluates the Morse potential with periodic boundary conditions.
- `evalMorsePBC_sym`: Evaluates the Morse potential with periodic boundary conditions using symmetrized atoms.
- `evalMorsePBCatoms_sym`: Evaluates the Morse potential with periodic boundary conditions over a set of atoms.
- `evalGridR`: Evaluates the grid force field based on radial distances from atomic positions.
- `evalCombindGridFF`: Combines Pauli, London, and Coulomb forces into a single output vector.
- `evalGridFFs_symetrized`: Evaluates the force field on the grid with symmetrized atoms.

##### Grid Generation
- `makeGridFF_omp`: Computes the force field on the grid using OpenMP parallelization.
- `makeGridFF`: Computes the force field on the grid using OpenMP parallelization with single precision.
- `makeGridFF_omp_d`: Computes the force field on the grid using OpenMP parallelization with double precision.
- `makeGridFF_d`: Computes the force field on the grid using double precision.
- `makeGridFF_Hherm_d`: Computes the force field on the grid using Hermite interpolation in double precision.
- `makeVPLQHeval`: Evaluates the interpolated force field data and stores it in `VPLQH`.
- `makeVPLQH`: Evaluates and stores interpolated force field data in `VPLQH`.
- `makeGridFF_Bspline_HH_d`: Computes the force field on the grid using B-spline interpolation in double precision with Hermite coefficients.
- `makeGridFF_Bspline_d`: Computes and fits B-spline coefficients for force field potentials on the grid.

##### B-Spline Fitting
- `FitBsplines`: Fits B-spline coefficients to the grid data using an iterative fitting algorithm.

##### Data Copying
- `copyPitch`: Copies data from one pitch to another with specified strides.
- `copyPitchTransp`: Copies data with transposed dimensions for interpolation.
- `copyPBC`: Copies data with periodic boundary conditions applied.
- `pack_Bspline_d`: Packs B-spline coefficients into a single array for easier handling.
- `Bspline_to_f4`: Converts B-spline data to `Quat4f` format and allocates memory if necessary.

##### Electrostatic Potential
- `makeCoulombEwald`: Calculates electrostatic potentials using Ewald summation if available.

##### File I/O
- `tryLoadGridFF_potentials`: Attempts to load precomputed force field data from files or compute it if necessary.
- `tryLoad`: Attempts to load precomputed force field data from files or compute it if necessary with single precision.
- `tryLoad_new`: Attempts to load precomputed force field data from files or compute it if necessary, supporting different interpolation methods and options.
- `saveXsfDebug`: Saves grid data in XSF format for debugging purposes.

##### Electric Field Profile
- `getEFprofile`: Gets the electric field profile along a line segment between two points.
- `getEFprofileToFile`: Saves the electric field profile to a file.
- `checkEFProfileVsNBFF`: Checks if the electric field profile computed by GridFF matches that computed by NBFF within a specified tolerance.
- `checkZProfilesOverAtom`: Checks if the z-profiles over an atom match within a specified tolerance.
- `log_z`: Logs z-coordinates and corresponding force field values to a file.

##### Dipole Evaluation
- `evalCellDipole`: Evaluates the dipole moment for the current grid cell.

##### Other
- `findTop`: Finds the maximum z-coordinate among all atomic positions.
- `checkSum`: Checks the sum of force field components to ensure they are consistent.