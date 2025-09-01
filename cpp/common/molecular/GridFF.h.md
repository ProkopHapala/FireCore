# GridFF.h

## Purpose

`GridFF.h` is a header file that defines a class `GridFF`, which provides functionality for evaluating force fields on a grid using various interpolation methods. This class is designed to handle complex interactions between atoms in a system, such as Coulombic and Morse potential forces, and can be used with different force field models like linear, Hermite, or B-spline representations.

## Includes

- "fastmath.h"
- "Vec2.h"
- "Vec3.h"
- "quaternion.h"
- "Grid.h"
- "Forces.h"
- "MMFFparams.h"
- "Multipoles.h"
- "InterpolateTricubic.h"
- "InterpolateTrilinear.h"
- "Bspline.h"
- "Bspline_fit.h"
- "Bspline_fit_2D.h"
- "Bspline_fit_3D.h"
- "VecN.h"
- "NBFF.h"
- "EwaldGrid.h"
- "IO_utils.h"

## Free Functions

### `T sum (int n, T* data, T t )`

Sums the first `n` elements of an array `data` and adds it to a variable `t`.

#### Parameters
- `int n`: Number of elements in the array.
- `T* data`: Pointer to the array containing the data.
- `T t`: Variable to which the sum is added.

#### Return Value
- Returns the updated value of `t`.

### `void autoNPBC (const Mat3d& cell, Vec3i& nPBC, double Lmin=30.0 )`

Automatically calculates the number of periodic boundary conditions based on minimum length.

#### Parameters
- `const Mat3d& cell`: The cell dimensions represented by a 3x3 matrix.
- `Vec3i& nPBC`: Reference to an integer vector that will store the calculated number of periodic boundary conditions in each direction (x, y, z).
- `double Lmin=30.0`: Minimum length for calculating the number of periodic boundary conditions.

#### Return Value
- None

### `double evalDipole (int n, Vec3d* ps, Quat4d* REQs, Vec3d& Dout, Vec3d& p0out )`

Calculates the dipole moment and center of charge for a set of particles.

#### Parameters
- `int n`: Number of particles.
- `Vec3d* ps`: Array of particle positions.
- `Quat4d* REQs`: Array of quaternion moments (charge distribution).
- `Vec3d& Dout`: Output vector to store the calculated dipole moment.
- `Vec3d& p0out`: Output vector to store the center of charge.

#### Return Value
- Returns the total charge `Q`.

## Types (Classes and Structs)

### class `GridFF`

`GridFF` is a class that manages force field calculations on a grid. It inherits from `NBFF`, which provides basic functionality for handling atomic positions, types, and forces.

**Inheritance**

- NBFF

#### Properties

- **`public: GridShape grid`**: Represents the shape of the grid.
- **`bool bUseEwald`**: Indicates whether to use Ewald summation for long-range interactions.
- **`EwaldGrid* ewald`**: Pointer to an `EwaldGrid` object used for Ewald summation if enabled.
- **`Quat4i cubic_yqis[4]`**: Array of quaternions representing cubic indices in the y-direction.
- **`Quat4i cubic_xqis[4]`**: Array of quaternions representing cubic indices in the x-direction.
- **`int iDBG=-1`**: Debugging flag for printing debug information.
- **`GridFFmod mode = GridFFmod::BsplineDouble`**: Specifies the type of force field model to use (e.g., linear, Hermite, B-spline).
- **`int perVoxel = 4`**: Number of data points stored in each voxel.
- **`std::vector<Vec3d> apos_`**: Vector containing atomic positions.
- **`std::vector<Quat4d> REQs_`**: Vector containing quaternion moments for each atom.
- **`std::vector<int> atypes_`**: Vector containing atomic types.
- **`Vec3d dip_p0`**: Center of the dipole moment.
- **`Vec3d dip`**: Dipole moment vector.
- **`double Q`**: Total charge in the system.
- **`double Mpol[10]`**: Array to store multipole moments up to quadrupole.
- **`int iDebugEvalR = 0`**: Debugging flag for evaluating radial forces.
- **`bool bCellSet = false`**: Indicates whether the cell has been set.
- **`bool bSymetrized = false`**: Indicates whether atoms have been symmetrized.

#### Methods

##### `double findTop ()`

Finds the maximum z-coordinate of all atomic positions.

##### `void bindSystem (int natoms_, int* atypes_, Vec3d* apos_, Quat4d* REQs_ )`

Binds an external system to the grid by setting its properties.

##### `void allocateFFs (bool bDouble=false )`

Allocates memory for force field data based on whether double precision is used.

##### `void clear ()`

Deallocates all allocated memory and resets internal variables.

##### `void allocateAtoms (int natoms_)`

Allocates memory for atomic positions and quaternion moments.

##### `int loadCell (const char * fname )`

Loads the cell dimensions from a file.

##### `void evalCellDipole ()`

Evaluates the dipole moment of the current cell.

##### `void init (Vec3i n, Mat3d cell, Vec3d pos0, bool bDouble=false )`

Initializes the grid with given dimensions and position.

##### `void setAtoms (int natoms_, Vec3d * apos_, Quat4d * REQs_ )`

Sets atomic positions and quaternion moments for the system.

##### `void evalAtPoints_REQ (int n, Vec3d* ps, Quat4d* FFout, Quat4d REQ, int natoms_, Vec3d * apos_, Quat4d * REQs_, Vec3i* nPBC_=0 )`

Evaluates the force field at specified points using a given quaternion.

##### `void evalAtPoints_REQ (int n, Vec3d* ps, Quat4d* FFout, Quat4d REQH, Vec3i* nPBC=0 )`

Evaluates the force field at specified points using a given quaternion and periodic boundary conditions.

##### `void evalAtPoints (int n, const Vec3d* ps, Quat4d* FFout, Quat4d PLQH, int natoms_, const Vec3d * apos_, const Quat4d * REQs_, Vec3i* nPBC=0 )`

Evaluates the force field at specified points using a given linear combination of force field components.

##### `void evalAtPoints (int n, const Vec3d* ps, Quat4d* FFout, Quat4d PLQH, Vec3i* nPBC=0 )`

Evaluates the force field at specified points using a given linear combination of force field components and periodic boundary conditions.

##### `void evalAtPoints_Split (int n, Vec3d* ps, Quat4d* FFout, Quat4d PLQH, int natoms_, Vec3d * apos_, Quat4d * REQs_, Vec3i* nPBC=0 )`

Evaluates the force field at specified points using a given linear combination of force field components and splits the evaluation into Coulombic and Morse parts.

##### `void evalAtPoints_Split (int n, Vec3d* ps, Quat4d* FFout, Quat4d PLQH, Vec3i* nPBC=0 )`

Evaluates the force field at specified points using a given linear combination of force field components and splits the evaluation into Coulombic and Morse parts with periodic boundary conditions.

##### `void makeGridFF_omp (int natoms_, Vec3d * apos_, Quat4d * REQs_ )`

Generates the grid force field data in parallel using OpenMP.

##### `void makeGridFF ()`

Generates the grid force field data sequentially.

##### `void makeGridFF_omp_d (int natoms_, Vec3d * apos_, Quat4d * REQs_ )`

Generates double-precision grid force field data in parallel using OpenMP.

##### `void makeGridFF_d ()`

Generates double-precision grid force field data sequentially.

##### `void makeGridFF_Hherm_d (int natoms_, Vec3d * apos_, Quat4d * REQs_ )`

Generates Hermite double-precision grid force field data in parallel using OpenMP.

##### `void makeGridFF_Hherm_d ()`

Generates Hermite double-precision grid force field data sequentially.

##### `void evalBsplineRef (int natoms_, Vec3d * apos_, Quat4d * REQs_, double* VPaul, double* VLond, double* VCoul )`

Evaluates the B-spline reference potential for Coulombic and Morse interactions.

##### `void makeVPLQHeval (int natoms_, Vec3d * apos_, Quat4d * REQs_ )`

Generates evaluation data for force field components using cubic interpolation.

##### `void copyPitch (int n, T* dst, int i0dst, int mdst, const T* src, int i0src, int msrc )`

Copies elements from one array to another with a specified pitch.

##### `void copyPitchTransp (Vec3i ndst, Vec3i transp, T* dst, int i0dst, int mdst, const T* src, int i0src, int msrc )`

Copies elements from one array to another with a specified pitch and transpose order.

##### `void copyPBC (Vec3i nsrc, T* src, int i0src, int msrc, Vec3i ndst, T* dst, int i0dst, int mdst )`

Copies elements from one array to another while applying periodic boundary conditions.

##### `void makeVPLQH ()`

Generates evaluation data for force field components using cubic interpolation with periodic boundary conditions.

##### `void FitBsplines (double Ftol=1e-8, int nmaxiter=1000, double dt=0.1 )`

Fits B-spline coefficients to the grid force field data.

##### `void makeGridFF_Bspline_HH_d (double Ftol=1e-8, int nmaxiter=1000, double dt=0.3 )`

Generates B-spline reference potential using Hermite interpolation in parallel with OpenMP.

##### `void makeCoulombEwald (int natoms_, Vec3d * apos_, Quat4d * REQs_, double* VCoul=0 )`

Calculates the Coulombic interaction energy using Ewald summation.

##### `void tryLoadGridFF_potentials (int natoms_, Vec3d * apos_, Quat4d * REQs_, int ntot, double*& VPaul, double*& VLond, double*& VCoul, bool bSaveNPY=false, bool bSaveXSF=false, bool bForceRecalc=false )`

Attempts to load pre-calculated grid force field potentials from files.

##### `bool makeGridFF_Bspline_d (int natoms_, Vec3d * apos_, Quat4d * REQs_, bool bSaveNPY=false, bool bSaveXSF=false, bool bFit=true, bool bRefine=false )`

Generates B-spline reference potential and fits it to the grid force field data.

##### `void pack_Bspline_d ()`

Packs B-spline coefficients into a single array for easier handling.

##### `int Bspline_to_f4 (bool bAlloc )`

Converts B-spline coefficients to float format, optionally allocating memory if needed.

##### `double evalMorsePBC (Vec3d pi, Quat4d REQi, Vec3d& fi, int natoms, Vec3d * apos, Quat4d * REQs )`

Evaluates the Morse potential for a given point and atom set with periodic boundary conditions.

##### `double evalMorsePBC_sym (Vec3d  pi, Quat4d  REQi, Vec3d& fi     )`

Evaluates the Morse potential for a given point using symmetrized atoms.

##### `double evalMorsePBCatoms_sym (int na, Vec3d* ps, Quat4d* REQs, Vec3d* forces )`

Evaluates the Morse potential for multiple atom sets with periodic boundary conditions and symmetrization.

##### `void evalGridR (int natoms, Vec3d * apos, Quat4d * REQs )`

Evaluates radial force contributions to the grid force field data.

##### `void evalCombindGridFF (Quat4d REQ, Quat4f * FF )`

Combines different components of the force field into a single output vector.

##### `void setAtomsSymetrized (int n, int* atypes, Vec3d* apos, Quat4d* REQs, double d=0.1 )`

Symmetrizes atomic positions and moments by adding replicas within a specified distance.

##### `void getEFprofile (int n, Vec3d p0, Vec3d p1, Quat4d REQ, Quat4d* fes, bool bPrint=false)` 

Generates an electric field profile between two points for debugging purposes.

##### `void getEFprofileToFile (const char* fname, int n, Vec3d p0, Vec3d p1, Quat4d REQ )`

Saves the electric field profile to a file.

##### `double checkEFProfileVsNBFF (int n, Vec3d p0, Vec3d p1, const Quat4d& REQ, double tol=1e-2, bool bExit=false, bool bPrint=false, bool bWarn=true, const char* logfiflename="checkEFProfileVsNBFF.log" )`

Checks the accuracy of force field evaluations against a reference implementation.

##### `bool checkZProfilesOverAtom (int ia, int n, double zmin, double zmax, const Quat4d& REQ, double tol=1e-2, bool bExit=true, bool bPrint=true )`

Checks the accuracy of force field evaluations along the z-axis for a specific atom.

##### `void log_z (const char* fname, int ix=0, int iy=0)` 

Logs grid data to a file for debugging purposes.

##### `void checkSum (bool bDouble )`

Checks the sum of all components in the force field grid.

##### `void initGridFF (const char * name, double z0=NAN, bool bAutoNPBC=true, bool bSymetrize=false, double rAutoPBC=20.0 )`

Initializes the grid with specified parameters and optionally symmetrizes atoms.

##### `bool tryLoad (const char* fname_Coul, const char* fname_Paul, const char* fname_Lond, bool recalcFF=false, bool bDouble=false )`

Attempts to load pre-calculated force field potentials from files or recalculates them if necessary.

##### `bool tryLoad_new (bool bSymetrize=true, bool bFit=true, bool bRefine=true, bool bPrint=false )`

Attempts to load or generate new force field potentials based on specified parameters.