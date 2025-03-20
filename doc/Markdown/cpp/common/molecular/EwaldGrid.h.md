# EwaldGrid.h

`EwaldGrid.h` implements the Ewald summation method for calculating electrostatic potentials in molecular dynamics simulations. It projects atomic charges onto a grid and solves Poisson's equation using FFT-based techniques.

## Includes

- `fastmath.h`
- `Vec2.h`
- `Vec3.h`
- `quaternion.h`
- `Grid.h`
- `Bspline.h`
- `IO_utils.h`
- `<omp.h>`
- `<fftw3.h>`

---

## Free functions

### Index Manipulation
- `pbc_ifw`: Calculates the forward periodic boundary condition index.
- `pbc_ibk`: Calculates the backward periodic boundary condition index.

### FFTW Array Conversion
- `array2fftc`: Converts a real array to an FFTW complex array.
- `fftc2array`: Converts an FFTW complex array back to a real array.
- `fftc2array_mul`: Multiplies an FFTW complex array by a scalar and converts it back to a real array.

---

## Types (classes and structs)

---

### class `EwaldGrid`

`EwaldGrid` is a class that handles the Ewald summation method for electrostatic potentials, particularly useful in molecular dynamics simulations. It provides methods to project atomic charges onto a grid and solve Poisson's equation using FFT-based techniques.

**Inheritance**

- GridShape

#### Properties

##### Charge Information
- `Qtot`: Total charge on the grid.
- `Qabs`: Absolute value of total charge on the grid.
- `Qtot_g`: Total charge after applying periodic boundary conditions (PBC).
- `Qabs_g`: Absolute value of total charge after applying PBC.
- `dipole`: Dipole moment vector for slab correction and debugging.

##### FFTW Data
- `V_work`: Temporary array used for FFT operations.
- `vV_work`: Another temporary array used for FFT operations.
- `ifft_plan`: FFTW plan for inverse Fourier transform.

##### PBC Indexing
- `bCubicPBCIndexesDone`: Boolean flag indicating whether cubic PBC indices have been computed.
- `xqs_o3`, `yqs_o3`, `zqs_o3`: Arrays of periodic boundary condition indices for cubic B-spline basis functions.
- `bQuinticPBCIndexesDone`: Boolean flag indicating whether quintic PBC indices have been computed.
- `xqs_o5`, `yqs_o5`, `zqs_o5`: Arrays of periodic boundary condition indices for quintic B-spline basis functions.

##### Debugging
- `iDBG`: Debugging index used to track the current atom being processed.

#### Methods

##### Atom Projection
- `split_atoms_parallel`: Splits atoms into non-overlapping regions to prevent memory write collisions during projection.
- `project_atoms_on_grid_linear`: Projects atomic charges onto a linear interpolation grid.
- `project_atoms_on_grid_cubic`: Projects atomic charges onto a cubic interpolation grid.
- `project_atoms_on_grid_quintic`: Projects atomic charges onto a quintic interpolation grid.
- `projectAtoms`: Projects atomic charges onto a specified interpolation grid based on the B-spline order.

##### System Setup
- `setup`: Initializes the EwaldGrid object with given parameters.

##### Laplace Solver (Real Space)
- `laplace_real`: Solves Poisson's equation using finite differences in real space (no PBC).
- `laplace_real_pbc`: Solves Poisson's equation using finite differences in real space with periodic boundary conditions.
- `laplace_real_loop`: Iteratively solves Laplace's equation using finite differences in real space.
- `laplace_real_loop_inert`: Iteratively solves Laplace's equation using finite differences in real space with inertial relaxation.

##### Slab Correction
- `slabPotential`: Applies a correction potential due to slab geometry in the z-direction.

##### Laplace Solver (Reciprocal Space)
- `laplace_reciprocal_kernel`: Applies a reciprocal space kernel to the FFTW complex arrays.
- `prepare_laplace`: Prepares the necessary data structures and plans for solving Laplace's equation using FFT-based methods.
- `solve_laplace`: Solves Poisson's equation in reciprocal space using FFT-based methods.
- `destroy_laplace`: Frees memory allocated for FFTW plans and arrays.

##### Macro-Level Solver
- `solve_laplace_macro`: Combines multiple steps of the Ewald summation process into a single macro-level function.

##### Potential Calculation
- `potential_of_atoms`: Calculates the potential due to atomic charges on a grid, handling both FFT-based and real-space methods.