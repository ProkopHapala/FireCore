# Fortran DFT Implementation (Fireball)

This directory contains the Fortran implementation of the Fireball Density Functional Tight Binding (DFTB) code. **Note**: This is an optional component and is NOT the core part of the FireCore project.

## Overview

The Fireball code implements:
- **Density Functional Tight Binding (DFTB)** - Semi-empirical quantum method
- **Self-Consistent Field (SCF)** calculations
- **Real-space grid operations** - Electron density and potential on grids
- **Molecular orbital analysis** - Orbital visualization and analysis
- **QM/MM integration** - Interface with classical force fields

## Directory Structure

```
fortran/
├── MAIN/                # Main programs and SCF drivers
│   ├── fireball.f90    # Main Fireball program
│   ├── libFireCore.f90 # Library interface for C++ integration
│   ├── mixer.f90       # SCF mixing algorithms
│   ├── solveH.f90      # Hamiltonian diagonalization
│   └── anderson2.f90   # Anderson mixing for SCF acceleration
├── MODULES/             # Fortran modules and data structures
│   ├── configuration.f90  # System configuration
│   ├── interactions.f90   # Interaction parameters
│   ├── density.f90        # Electron density handling
│   ├── energy.f90         # Energy components
│   └── forces.f90         # Force calculations
├── ASSEMBLERS/          # Hamiltonian and matrix assembly
│   ├── assemble_mcweda.f90 # Main Hamiltonian assembly
│   ├── assemble_h.f90      # Hamiltonian matrix assembly
│   ├── assemble_S.f90      # Overlap matrix assembly
│   └── buildh.f90          # Hamiltonian construction
├── INTERACTIONS/        # Two-center and three-center interactions
│   ├── doscentros.f90      # Two-center interactions
│   ├── trescentros.f90     # Three-center interactions
│   ├── unocentros.f90      # One-center interactions
│   └── get_ewald.f90       # Ewald summation
├── GRID/                # Real-space grid operations
│   ├── initgrid.f90        # Grid initialization
│   ├── project_dens.f90    # Density projection to grid
│   ├── project_orb.f90     # Orbital projection to grid
│   ├── laplace_fft.f90     # FFT-based Laplace solver
│   └── writeout_xsf.f90    # XSF file output
├── NEIGHBORS/           # Neighbor list construction
│   ├── neighbors.f90       # Main neighbor finding
│   ├── common_neighbors.f90 # Common neighbor utilities
│   └── find_neigh_max.f90  # Maximum neighbor finding
├── READFILES/           # Input file parsing
│   ├── readdata.f90        # Main data reading
│   ├── readbasis.f90       # Basis set reading
│   └── readparam.f90       # Parameter reading
├── INTERPOLATERS/       # Spline interpolation
│   ├── interpolate_1d.f90  # 1D interpolation
│   ├── interpolate_2d.f90  # 2D interpolation
│   └── buildspline_1d.f90  # Spline construction
├── ALLOCATIONS/         # Memory management
│   ├── allocate_h.f90      # Hamiltonian allocation
│   ├── allocate_rho.f90    # Density allocation
│   └── allocate_neigh.f90  # Neighbor allocation
├── DASSEMBLERS/         # Force assembly (derivatives)
│   ├── getforces.f90       # Main force calculation
│   ├── Dassemble_2c.f90    # Two-center force derivatives
│   └── Dassemble_3c.f90    # Three-center force derivatives
├── ROTATIONS/           # Coordinate transformations
│   ├── rotate.f90          # Rotation matrices
│   ├── chooser.f90         # Direction cosine matrices
│   └── epsilon.f90         # Rotation utilities
├── MATH/                # Mathematical utilities
│   ├── factorial.f90       # Factorial functions
│   └── cross.f90           # Vector cross products
├── INITIALIZERS/        # System initialization
│   ├── initbasics.f90      # Basic initialization
│   ├── initcharges.f90     # Charge initialization
│   └── initneighbors.f90   # Neighbor initialization
└── doc/                 # Fortran-specific documentation
    ├── Fireball_general.md
    ├── Fireball_code_structure.md
    └── Fireball_equations.md
```

## Key Components

### Main Program (MAIN/)
- **`fireball.f90`** - Standalone Fireball executable
- **`libFireCore.f90`** - Library interface for C++ integration
- **SCF drivers** - Self-consistent field iteration control
- **Matrix solvers** - Eigenvalue problem solvers

### Core Algorithms (ASSEMBLERS/)
- **Hamiltonian assembly** - Construction of quantum mechanical Hamiltonian
- **Overlap matrix** - Atomic orbital overlap calculations
- **Energy evaluation** - Total energy computation
- **Force calculation** - Analytical force derivatives

### Interaction Calculations (INTERACTIONS/)
- **Two-center integrals** - Bond interactions between atom pairs
- **Three-center integrals** - Three-body interactions
- **Ewald summation** - Long-range electrostatic interactions
- **Spline interpolation** - Efficient integral evaluation

### Grid Operations (GRID/)
- **Real-space grids** - 3D grids for density and potential
- **FFT operations** - Fast Fourier transforms for Poisson solving
- **Density projection** - Mapping orbitals to real space
- **Visualization output** - XSF format for visualization

## Build System

### Compilation
```bash
# From FireCore root directory
./make.sh
```

This script:
1. Creates `Build/` directory
2. Compiles all Fortran sources
3. Links with Intel MKL libraries
4. Creates `fireball.x` executable

### Dependencies
- **Fortran compiler** - `gfortran` or `ifort`
- **Intel MKL** - Math Kernel Library for LAPACK/BLAS
- **FFTW** - Fast Fourier Transform library (optional)

### Installation
```bash
# Install dependencies (Ubuntu)
sudo apt-get install gfortran intel-mkl libmkl-intel-lp64 libmkl-core libmkl-intel-thread

# Compile
./make.sh

# Test
cd tests/Fireball/t02_CH4
./run.sh
```

## Usage

### Standalone Execution
```bash
cd tests/Fireball/t02_CH4
./run.sh  # Runs fireball.x with proper setup
```

### Python Interface
```python
import pyBall.FireCore as fc

# Initialize calculation
fc.init_system("molecule.xyz")
fc.setup_calculation()

# Run SCF
fc.run_scf()

# Get results
energy = fc.get_total_energy()
forces = fc.get_forces()
```

### C++ Integration (QM/MM)
```cpp
#include "FireCore_interface.h"

// Initialize Fireball
firecore_init("molecule.xyz");
firecore_setup();

// Run calculation
firecore_scf();

// Get forces for MM integration
double* forces = firecore_get_forces();
```

## Input Files

### Required Files
- **`input.bas`** - Molecular geometry in basis format
- **`fireball.in`** - Calculation parameters
- **`Fdata/`** - Basis set and interaction data

### Optional Files
- **`answer.bas`** - Initial guess geometry
- **`CHARGES`** - Initial charge distribution

### Example fireball.in
```
&OPTION
itheory = 1          ! DFTB theory level
max_scf_iterations = 100
scf_tolerance_rms = 1.0d-5
&END

&OUTPUT
iwrtxyz = 1          ! Write XYZ trajectory
iwrthop = 1          ! Write Hamiltonian
&END
```

## Output Files

### Standard Output
- **`answer.bas`** - Optimized geometry
- **`answer.xyz`** - XYZ format geometry
- **Energy and force information**

### Grid Output
- **`*.xsf`** - XCrySDen format for visualization
- **Electron density grids**
- **Molecular orbital grids**

### Analysis Output
- **`CHARGES`** - Atomic charges
- **`*.out`** - Detailed calculation logs

## Integration with FireCore

### QM/MM Interface
The Fortran code integrates with C++ classical force fields:
1. **Geometry optimization** - Quantum region optimization
2. **Force calculation** - QM forces for MM integration
3. **Charge distribution** - Electrostatic embedding
4. **Energy evaluation** - QM contribution to total energy

### Grid Interface
Real-space grids are shared with C++ visualization:
1. **Density grids** - Electron density for AFM simulation
2. **Potential grids** - Electrostatic potential for GridFF
3. **Orbital grids** - Molecular orbital visualization

## Performance Considerations

### Parallelization
- **OpenMP** - Shared memory parallelization
- **MKL threading** - Optimized linear algebra
- **Grid operations** - Parallel FFT operations

### Memory Management
- **Dynamic allocation** - Efficient memory usage
- **Large systems** - Careful memory planning
- **Grid storage** - Compressed grid formats

### Optimization
- **Spline interpolation** - Fast integral evaluation
- **Neighbor lists** - Efficient interaction screening
- **SCF acceleration** - Anderson mixing and DIIS

## Troubleshooting

### Compilation Issues
- **MKL linking** - Check Intel MKL installation
- **Fortran compiler** - Ensure compatible compiler version
- **Library paths** - Verify LD_LIBRARY_PATH

### Runtime Issues
- **SCF convergence** - Adjust mixing parameters
- **Memory errors** - Check system memory limits
- **File I/O** - Verify input file formats

### Integration Issues
- **C++ linking** - Check library compatibility
- **Python bindings** - Verify f2py installation
- **Grid alignment** - Check coordinate system consistency
