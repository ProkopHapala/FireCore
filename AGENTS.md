# FireCore Coding Guide
  
# FireCore Project Structure and Development Guide

## Project Overview

FireCore is a comprehensive computational chemistry and physics repository containing multiple subprojects, programs, and test scripts. It implements various numerical methods for:

- **Density Functional Theory (DFT)** - implemented in Fortran, can be called from C++ and Python via library
- **Classical Molecular Dynamics** - implemented in C++, and OpenCL
- **Force Field Development and Fitting** - including reactive and non-reactive force fields
- **Quantum-Classical Hybrid Methods (QM/MM)** - combining DFT with classical force fields (C++ call fortran library)
- **Numerical Method Demonstrations** - standalone Python scripts for various computational techniques

**Important**: FireCore is NOT a single build target but a repository of interconnected subprojects with different compilation and execution workflows.

## Directory Structure and Components

### 1. Core C++ Implementation (`/cpp/`)

The computational core of FireCore, containing high-performance libraries and applications:

#### standalone visual applications (`/cpp/apps/`)
- **Visual programs using SDL**: Interactive molecular editors and simulators
- **GPU-accelerated variants**: 
  - `/cpp/apps_OCL/` - OpenCL accelerated applications
  - `/cpp/apps_CUDA/` - CUDA accelerated applications

#### Libraries (`/cpp/libs/`)
- **`/cpp/libs/Molecular/`** - Core molecular mechanics libraries
- **`/cpp/libs_SDL/`** - SDL-based visualization libraries  
- **`/cpp/libs_OCL/`** - OpenCL accelerated libraries

#### Core Implementation
- **`/cpp/common/`** - Computational core (algorithms, data structures, math)
- **`/cpp/common_SDL/`** - SDL-based GUI and visualization components
- **`/cpp/common_resources/`** - Shared resources (force field parameters, molecular data)

### 2. Fortran DFT Implementation (`/fortran/` and `/fortran2/`)

**Note**: The Fortran module implements Density Functional Tight Binding (DFTB), but it is NOT the core part of this project.

- **`/fortran/`** - Current Fireball DFTB implementation
- **`/fortran2/`** - Reorganized/refactored version (work in progress)

Key subdirectories in `/fortran/`:
- `MAIN/` - Main program and SCF drivers
- `MODULES/` - Fortran modules and data structures
- `ASSEMBLERS/` - Hamiltonian and matrix assembly routines
- `INTERACTIONS/` - Two-center and three-center interaction calculations
- `GRID/` - Real-space grid operations and FFT
- `NEIGHBORS/` - Neighbor list construction

### 3. Python Interface (`/pyBall/`)

Python bindings and utilities providing access to both C++ and Fortran components:

#### Core Python Modules
- **`FireCore.py`** - Main interface to Fireball DFT
- **`MMFF.py`, `MMFFsp3.py`** - Molecular mechanics force fields
- **`AtomicSystem.py`** - Molecular structure manipulation
- **`Forces.py`, `Forces_cpp.py`** - Force calculation interfaces

#### Specialized Modules
- **`/pyBall/OCL/`** - Pure pyOpenCL implementations (independent of C++ libraries)
- **`/pyBall/DFT/`** - DFT-related utilities and high-level interfaces
- **`/pyBall/GUI/`** - Python-based GUI components

#### Utilities
- **`atomicUtils.py`** - Atomic data and utilities
- **`plotUtils.py`** - Visualization and plotting helpers
- **`buildUtils.py`** - Build system utilities (to recompile C++ libraries when imported from python modules)

### 4. OpenCL Implementation (`/pyBall/OCL/`)

Standalone pyOpenCL implementations for GPU acceleration:
- Independent of C++ libraries
- Direct OpenCL kernel implementations
- Alternative to C++ OpenCL bindings

### 5. Tests and Examples (`/tests/`)

**This is the best place to understand what functionality is implemented and how to use it.**

#### Test Categories
- **`/tests/Fireball/`** - Fireball DFT tests (H2, CH4, pentacene, etc.)
- **`/tests/tMMFF*/`** - Molecular mechanics tests (using python bindings and scripts, no GUI)
- **`/tests/tMolGUIapp*/`** - GUI application tests
- **`/tests/tFitREQ*/`** - Non-covalent Force field fitting (both using python+C++ and pyOpenCL bindings)
- **`/tests/tDFT*/`** - DFT calculation tests
- **`/tests/tEFF*/`** - Electron Force Field tests
- **`/tests/pyFireball/`** - Python DFT interface tests
- **`/tests/pySCF/`** - PySCF integration tests
- **`/tests/dftb/`** - DFTB+ integration tests

### 6. Documentation (`/doc/`)

- **`/doc/Markdown/`** - Technical documentation (GUI, force fields, etc.)
- **`/doc/DevNotes/`** - Development notes and TODO items
- **`/doc/py/`** - Python sketches and demonstrations
- **`/doc/Julia/`** - Julia implementations of algorithms
- **`/doc/Maxima/`** - Mathematical derivations and symbolic computations

### 7. Build System

- **`/cpp/Build/`** - CMake build directory for C++ components
- **`/build/`** - Fortran build directory
- **`make.sh`** - Main build script for Fortran components
- **`/cpp/CMakeLists.txt`** - CMake configuration for C++ components

## Critical Development Rules

### Compilation and Execution Workflow

**⚠️ NEVER compile manually using `make` - always use project bash scripts!**

#### For C++ Applications and Libraries:
1. **Use provided `run.sh` scripts** in test directories
2. Scripts automatically:
   - Recompile the program/library
   - Set up paths to inputs and libraries  
   - Run the program with proper arguments
   - Output stdout to logfiles
3. **Edit `run.sh` scripts if needed** for different configurations
4. **DO NOT** try to compile using `make` directly - paths may be incorrect

#### For Fortran (Fireball):
1. Use `make.sh` in project root for initial compilation
2. Use `run.sh` scripts in `/tests/Fireball/` subdirectories
3. Ensure `Fdata_HC_minimal` is available (download from fireball-qmd.github.io)

#### For Python:
1. Set `PYTHONPATH` environment variable: `export PYTHONPATH=/path/to/FireCore:$PYTHONPATH`
2. Use `run.sh` scripts in test directories
3. Python scripts typically call compiled C++ libraries through bindings

### Example Workflow:
```bash
# Test molecular mechanics
cd tests/tMMFF
./run.sh

# Test GUI application  
cd tests/tMolGUIapp
./run.sh

# Test Fireball DFT
cd tests/Fireball/t02_CH4
./run.sh
```

## Global Coding Rules

### Core Principles
- **Minimal dependencies**: stdlib + NumPy/Matplotlib unless explicitly requested
- **Concise, modular code**: small, testable, step-by-step changes; divide-and-conquer
- **Reuse existing functions**: avoid duplication (check first, report if adaptation needed)
- **Pure, data-oriented functions**: explicit inputs/outputs; default named args to avoid long call strings
- **Fail loudly**: no silent handling; assertions for invariants; crashes with stack trace preferred to masking
- **Comment out deprecated/experimental code** instead of deleting; mark with TODO/DEBUG

### NEVER DO THIS
- Never delete, override, or rearrange existing code without explicit permission
- Never perform random aesthetic/style edits unrelated to the task
- Never apply "quick-fixes" that hide root causes (e.g., hard-coded outputs)

### Debugging First
- **Debuggability > UX**: do not hide issues
- Initially add debug prints for key variables and flow; remove later after debugging
- Avoid broad try/except as they mask bugs; prefer loud crashes with stack trace
- Make small, testable changes and run after every change
- Log function entries and key conditions when helpful to track flow
- Mark unfinished/experimental code clearly (e.g., # TODO, # DEBUG)

### Performance Guidelines
- Prefer data-oriented code that is cache-friendly and avoids overheads
- Preallocate and reuse buffers; avoid repeated allocation in hot paths
- Be explicit about dtypes/shapes; prefer contiguous memory where possible

### Style Guidelines (Cross-Language)
- Prefer concise/compact code; avoid bloated structures and unnecessary empty lines
- Short variable names OK (math/physics symbols like E, T, m) when locally clear
- Prefer one-liner expressions, assume unlimited line-width
- Avoid line wrapping that hurts readability of expressions
- Inline comments behind the code line for rationale
- **Doxygen**: use `///`; avoid `/* ... */`
- **C++**: use `printf` for debugging over `std::cout`; prefer plain C arrays (double*) in hot paths
- **Vector math in C/C++**: use `Vec3.h`, `Vec2.h` with `Vec3d`, `Vec2d`, and helpers like `dot()`, `norm()`, `cross()`

### Visualization
- Separate compute vs plotting; no plotting in core algorithms
- Plotting optional via flags (e.g., --noPlot, --saveFig)
- `plt.show()` only in CLI/main, never in library code
- Prefer shared plotting helpers (e.g., plot_utils.py) to avoid duplication

## Getting Started

1. **Explore functionality**: Browse `/tests/` directory to see what's implemented
2. **Run examples**: Use `run.sh` scripts in test directories
3. **Check documentation**: Look in `/doc/` for technical details
4. **Follow the rules**: Always use provided build scripts, never compile manually
5. **Test changes**: Write tests and run them to verify correctness

## Key Integration Points

- **QM/MM**: Fireball DFT + C++ classical force fields
- **GPU Acceleration**: OpenCL/CUDA implementations of force fields
- **Python Bindings**: Access to both C++ and Fortran components
- **Visualization**: SDL-based interactive applications
- **Force Field Fitting**: Tools for parameterizing classical potentials against QM data

## Agent Operating Guidelines

These rules apply to any autonomous assistant working inside this repo:

- **No fallbacks or silent fixes**: if required data/geometry is missing, throw a descriptive error immediately. Never synthesize “best-effort” defaults—the failure must be visible.
- **Debug-first mindset**: prioritize traceable crashes, add targeted logging only while investigating.
- **Preserve context**: never remove comments during active debugging; if you must replace logic, comment out the old block instead of deleting so we can revert instantly if needed.
- **Respect existing structure**: extend modules in-place rather than rewriting or rearranging unless explicitly authorized.
- **Use official scripts**: when commands are required, rely on the provided `run.sh`/`make.sh` helpers; never invoke `make` directly.
- **Document parity work**: when mirroring Python ↔ JS features, cite the reference file/function in comments so future maintainers can diff implementations quickly.
