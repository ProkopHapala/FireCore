# EwaldGrid.h

## Includes

- "fastmath.h"
- "Vec2.h"
- "Vec3.h"
- "quaternion.h"
- "Grid.h"
- "Bspline.h"
- "IO_utils.h"
- <omp.h>
- <fftw3.h>

## Purpose

The `EwaldGrid` class is designed to handle the Ewald summation method for electrostatic potentials in periodic systems. It provides functions to project atomic charges onto a grid, solve Poisson's equation using FFT and iterative methods, and compute the potential due to these charges.

## Free Functions

### `void array2fftc (int n, const double* in, fftw_complex* out)`

Converts an array of real numbers into a complex array suitable for FFTW operations. This function initializes the imaginary part of each complex number to zero.

- **Parameters:**
  - `n`: Number of elements in the input array.
  - `in`: Pointer to the input array of real numbers.
  - `out`: Pointer to the output array of complex numbers, where the imaginary parts are set to zero.

### `void fftc2array (int n, const fftw_complex* in, double* out)`

Converts a complex FFTW array back into an array of real numbers. This function extracts only the real part of each complex number.

- **Parameters:**
  - `n`: Number of elements in the input array.
  - `in`: Pointer to the input array of complex numbers.
  - `out`: Pointer to the output array of real numbers.

### `void fftc2array_mul (int n, const fftw_complex* in, double* out, double f)`

Multiplies each element of a complex FFTW array by a scalar and converts it back into an array of real numbers. This function extracts only the real part after multiplication.

- **Parameters:**
  - `n`: Number of elements in the input array.
  - `in`: Pointer to the input array of complex numbers.
  - `out`: Pointer to the output array of real numbers.
  - `f`: Scalar value to multiply each element by before conversion.

### `int pbc_ifw (int i, int n)`

Calculates the forward periodic boundary condition index for a given index and grid size. This function wraps an index around if it exceeds the grid size.

- **Parameters:**
  - `i`: Original index.
  - `n`: Grid size.
  
- **Return Value:** 
  - The wrapped index within the valid range of the grid.

### `int pbc_ibk (int i, int n)`

Calculates the backward periodic boundary condition index for a given index and grid size. This function wraps an index around if it is less than zero.

- **Parameters:**
  - `i`: Original index.
  - `n`: Grid size.
  
- **Return Value:** 
  - The wrapped index within the valid range of the grid.

## Types (Classes and Structs)

### class `EwaldGrid`

The `EwaldGrid` class is responsible for handling Ewald summation in a periodic system. It extends the `GridShape` class to provide additional functionality related to charge distribution, potential calculation, and iterative methods.

**Inheritance**

- **Inherits from:** GridShape

#### Properties

- `public: double Qtot`: Total charge on the grid.
- `double Qabs`: Absolute value of total charge.
- `double Qtot_g`: Total charge in reciprocal space.
- `double Qabs_g`: Absolute value of total charge in reciprocal space.
- `Vec3d dipole`: Dipole moment vector for slab correction and debugging.
- `double* V_work`: Pointer to working array for potential calculations.
- `double* vV_work`: Pointer to another working array for potential calculations.
- `bool bCubicPBCIndexesDone`: Flag indicating if cubic PBC indices are computed.
- `Quat4i xqs_o3`: Precomputed cubic PBC indices for the x-axis.
- `Quat4i yqs_o3`: Precomputed cubic PBC indices for the y-axis.
- `Quat4i zqs_o3`: Precomputed cubic PBC indices for the z-axis.
- `bool bQuinticPBCIndexesDone`: Flag indicating if quintic PBC indices are computed.
- `Vec6i xqs_o5`: Precomputed quintic PBC indices for the x-axis.
- `Vec6i yqs_o5`: Precomputed quintic PBC indices for the y-axis.
- `Vec6i zqs_o5`: Precomputed quintic PBC indices for the z-axis.
- `int iDBG`: Debugging index used to track atoms during computation.

#### Methods

- **`void split_atoms_parallel (int na, Vec3d* apos, double Rcut)`**
  
  Splits atomic positions into non-overlapping regions within a specified radius (`Rcut`) to ensure safe projection onto the grid without memory write conflicts.

- **`void project_atoms_on_grid_linear (int na, const Vec3d* apos, const double* qs, double* dens, bool bPBC=true)`**
  
  Projects atomic charges onto a linear interpolation grid. This method handles both periodic and non-periodic boundary conditions based on the `bPBC` flag.

- **`void project_atoms_on_grid_cubic (int na, const Vec3d* apos, const double* qs, double* dens, bool bPBC=true)`**
  
  Projects atomic charges onto a cubic interpolation grid. This method handles both periodic and non-periodic boundary conditions based on the `bPBC` flag.

- **`void project_atoms_on_grid_quintic (int na, const Vec3d* apos, const double* qs, double* dens, bool bPBC=true)`**
  
  Projects atomic charges onto a quintic interpolation grid. This method handles both periodic and non-periodic boundary conditions based on the `bPBC` flag.

- **`int setup (Vec3d pos0_, Mat3d dCell_, Vec3i ns_, bool bPrint=false)`**
  
  Sets up the EwaldGrid with initial parameters such as grid position, cell dimensions, and number of grid points. This method initializes internal variables and prepares for further computations.

- **`void projectAtoms (int na, Vec3d* apos, double* qs, double* dens, int order)`**
  
  Projects atomic charges onto the grid using specified interpolation order (`order`). This function calls one of the projection methods based on the order provided.

- **`double laplace_real (double* Vin, double* Vout, double cSOR)`**
  
  Solves Poisson's equation in real space. This method applies a relaxation factor `cSOR` to iteratively solve for the potential field.

- **`double laplace_real_pbc (double* Vin, double* Vout, double cSOR=0.0)`**
  
  Solves Poisson's equation with periodic boundary conditions in real space. This method handles grid points that are close to the boundaries by wrapping them around using PBC indices.

- **`int laplace_real_loop (double* V, int nmaxiter=1000, double tol=1e-6, bool bPBC=true, double cSOR=0.0)`**
  
  Iteratively solves Poisson's equation in real space until convergence is achieved. This method handles both periodic and non-periodic boundary conditions based on the `bPBC` flag.

- **`int laplace_real_loop_inert (double* V, int nmaxiter=1000, double tol=1e-6, bool bPBC=true, double cSOR=0.0, double cV=0.5)`**
  
  Solves Poisson's equation with an inertial term to accelerate convergence. This method handles both periodic and non-periodic boundary conditions based on the `bPBC` flag.

- **`void slabPotential (int nz_slab, double* Vin, double* Vout)`**
  
  Adds a correction potential for systems with slab geometry. This function accounts for the dipole moment of the system in the z-direction.

- **`void laplace_reciprocal_kernel (fftw_complex* VV)`**
  
  Applies the reciprocal space kernel to the FFTW complex array `VV`. This method modifies the array based on the wave numbers and their inverses.

- **`void prepare_laplace (int flags=-1)`**
  
  Prepares the grid for solving Poisson's equation using FFT. This function allocates memory and sets up FFT plans with specified flags.

- **`void solve_laplace (const double* dens, double* Vout=0)`**
  
  Solves Poisson's equation in reciprocal space using FFTW. This method applies the reciprocal kernel to the density array `dens`.

- **`void destroy_laplace ()`**
  
  Cleans up resources allocated for solving Poisson's equation.

- **`void solve_laplace_macro (double* dens, int nz_slab, double* Vout, bool bPrepare=true, bool bDestroy=true, int flags=-1, bool bOMP=false, int nBlur=0, double cSOR=0, double cV=0.95)`**
  
  Solves Poisson's equation with slab potential correction and optional iterative relaxation. This function handles both FFTW-based and real-space methods.

- **`void potential_of_atoms (int nz, double* VCoul, int natoms, Vec3d* apos, double* qs, int order=3, int nBlur=4, double cV=0.95, bool bClearCharge=true)`**
  
  Computes the electrostatic potential due to atomic charges on the grid. This function handles both FFTW-based and real-space methods.

## Implementation Details

- **FFT Operations:** The class uses FFTW for solving Poisson's equation in reciprocal space. It includes functions to prepare, execute, and destroy FFT plans.
- **Iterative Methods:** For systems where FFT is not feasible or efficient, the class provides iterative methods using relaxation factors (`cSOR`) and inertial terms (`cV`).
- **Slab Correction:** The `slabPotential` function accounts for dipole moments in slab geometries by adding a correction potential.
- **Debugging:** The class includes debugging features such as saving charge densities to files for further analysis.