# C++ Core Implementation

This directory contains the computational core of FireCore - high-performance C++ libraries and applications for molecular mechanics, visualization, and GPU acceleration.

## Directory Structure

```
cpp/
├── common/              # Core computational libraries
│   ├── molecular/      # Molecular mechanics algorithms
│   ├── math/           # Mathematical utilities and algorithms
│   ├── dataStructures/ # Data structures and containers
│   └── OpenCL/         # OpenCL integration utilities
├── common_SDL/          # SDL-based GUI and visualization
│   ├── SDL2/           # SDL2 utilities and wrappers
│   └── SDL2OGL/        # SDL2 + OpenGL integration
├── common_resources/    # Shared resources and data
│   ├── cl/             # OpenCL kernels
│   ├── mol/            # Molecular structure files
│   ├── xyz/            # XYZ coordinate files
│   └── *.dat           # Force field parameters
├── apps/                # Interactive applications
│   ├── MolecularEditor/ # Main molecular editor/simulator
│   ├── EFF/            # Electron Force Field applications
│   └── RFF/            # Reactive Force Field applications
├── apps_OCL/            # OpenCL accelerated applications
├── apps_CUDA/           # CUDA accelerated applications
├── libs/                # Compiled libraries
│   ├── Molecular/      # Molecular mechanics libraries
│   └── quadrature_lib.cpp # Numerical integration
├── libs_SDL/            # SDL-based libraries
├── libs_OCL/            # OpenCL libraries
├── sketches_SDL/        # Development sketches and prototypes
└── Build/               # CMake build directory
```

## Key Components

### Core Libraries (`common/`)

#### Molecular Mechanics (`common/molecular/`)
- **`MMFFsp3.h`** - Molecular mechanics force field with sp3 hybridization
- **`GridFF.h`** - Grid-based force fields for substrate interactions
- **`MMFFBuilder.h`** - Automatic topology generation and charge assignment
- **`MolWorld_sp3.h`** - Molecular world simulation environment
- **`Atoms.h`** - Atomic data structures and manipulation

#### Mathematical Utilities (`common/math/`)
- **`Vec3.h`, `Vec2.h`** - Vector mathematics (Vec3d, Vec2d types)
- **`Mat3.h`** - 3x3 matrix operations
- **Spline interpolation** - B-splines, Hermite splines
- **Numerical integration** - Quadrature methods

#### Data Structures (`common/dataStructures/`)
- **Grid structures** - 3D grids for force fields and density
- **Neighbor lists** - Efficient neighbor finding algorithms
- **Molecular graphs** - Bond topology representation

### Applications (`apps/`)

#### MolecularEditor (`apps/MolecularEditor/`)
The main interactive molecular editor and simulator:
- **Real-time molecular dynamics**
- **Force field visualization**
- **QM/MM hybrid calculations**
- **AFM simulation**
- **Interactive molecule building**

#### Specialized Applications
- **EFF applications** - Electron Force Field methods
- **RFF applications** - Reactive Force Field implementations

### GPU Acceleration

#### OpenCL (`apps_OCL/`, `libs_OCL/`)
- **Parallel force evaluation**
- **Grid-based calculations**
- **Multi-molecule simulations**

#### CUDA (`apps_CUDA/`)
- **CUDA-accelerated molecular mechanics**
- **GPU memory management**
- **Kernel optimization**

## Build System

### CMake Configuration

```bash
cd cpp/Build
cmake .. [OPTIONS]
make -j$(nproc)
```

### Build Options
- **`-DWITH_SDL=ON/OFF`** - Enable/disable SDL GUI support
- **`-DWITH_OPENCL=ON/OFF`** - Enable/disable OpenCL acceleration
- **`-DWITH_CUDA=ON/OFF`** - Enable/disable CUDA acceleration
- **`-DWITH_FFTW=ON/OFF`** - Enable/disable FFTW support
- **`-DWITH_OPENMP=ON/OFF`** - Enable/disable OpenMP parallelization

### Dependencies

#### Required:
- **CMake** (>= 3.10)
- **C++ compiler** (g++, clang++)

#### Optional:
- **SDL2** - For GUI applications (`libsdl2-dev libsdl2-image-dev`)
- **OpenCL** - For GPU acceleration (`nvidia-opencl-dev libclfft-dev`)
- **CUDA** - For NVIDIA GPU acceleration
- **FFTW** - For FFT operations (`libfftw3-dev`)
- **OpenMP** - For CPU parallelization

## Usage Patterns

### Direct Compilation (Not Recommended)
```bash
cd cpp/Build/apps/MolecularEditor
make MolecularEditor
./MolecularEditor [options]
```

### Recommended: Use Test Scripts
```bash
cd tests/tMolGUIapp
./run.sh  # Automatically compiles and runs
```

## Key Libraries and Their Functions

### MMFFsp3 (Molecular Mechanics)
- **Bond stretching** - Harmonic and Morse potentials
- **Angle bending** - Harmonic angular potentials  
- **Torsional rotation** - Cosine series potentials
- **Non-bonded interactions** - Lennard-Jones + Coulomb
- **Charge equilibration** - Electronegativity-based charge assignment

### GridFF (Grid-based Force Fields)
- **Substrate interactions** - Rigid surface potentials
- **Electrostatic grids** - Pre-computed Coulomb potentials
- **Van der Waals grids** - Lennard-Jones interaction grids
- **Pauli repulsion** - Short-range repulsive potentials

### Visualization (SDL)
- **3D molecular rendering** - Ball-and-stick, space-filling models
- **Real-time dynamics** - Interactive molecular dynamics
- **Force visualization** - Vector field display
- **Grid visualization** - Isosurface and slice rendering

## Development Guidelines

### Code Organization
- **Header-only libraries** - Most functionality in `.h` files
- **Template-based** - Extensive use of C++ templates for performance
- **Modular design** - Clear separation of concerns
- **Minimal dependencies** - Avoid external library dependencies where possible

### Performance Considerations
- **Cache-friendly data structures** - Contiguous memory layouts
- **SIMD optimization** - Vectorized operations where applicable
- **GPU offloading** - Computationally intensive operations on GPU
- **Memory management** - Careful allocation/deallocation patterns

### Debugging and Testing
- **Debug prints** - Use `printf` for debugging output
- **Assertions** - Liberal use of `assert()` for invariants
- **Memory sanitizers** - Use AddressSanitizer for memory debugging
- **Profiling** - Built-in timing and performance counters

## Common Issues and Solutions

### Compilation Problems
- **Missing dependencies** - Install required libraries
- **CMake configuration** - Check build options and paths
- **Compiler compatibility** - Ensure C++11 or later support

### Runtime Issues
- **OpenCL errors** - Check GPU drivers and OpenCL installation
- **Memory leaks** - Use AddressSanitizer (`-fsanitize=address`)
- **Performance issues** - Check OpenMP thread count and GPU utilization

### Integration Issues
- **Python bindings** - Ensure shared libraries are built correctly
- **Resource paths** - Check symbolic links to `common_resources`
- **Fortran integration** - Verify library linking for QM/MM
