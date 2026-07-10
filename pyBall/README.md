# Python Interface (pyBall)

This directory contains Python bindings and utilities that provide access to both C++ and Fortran components of FireCore, along with standalone Python implementations of various computational methods.

## Directory Structure

```
pyBall/
├── Core Interfaces
│   ├── FireCore.py          # Main interface to Fireball DFT
│   ├── MMFF.py             # Molecular mechanics force fields
│   ├── MMFFsp3.py          # sp3 hybridization force fields
│   ├── Forces.py           # Force calculation interfaces
│   ├── Forces_cpp.py       # C++ force field bindings
│   ├── FFfit.py            # C++ wrapper + gauge-invariant hybrid Hessian fitting and bounded regularized LSQ
│   ├── FFfit_utils.py      # FF types/topology, hierarchy rows, torsions, analytic valence cross sensitivities
│   └── FFfit_plots.py      # Spectra, equilibrium distributions, and diagnostic Wilson-indicator visualizations
├── Molecular Tools
│   ├── AtomicSystem.py     # Molecular structure manipulation
│   ├── atomicUtils.py      # Atomic data and utilities
│   ├── elements.py         # Periodic table data
│   └── moleculeManagement.py # Molecule loading/saving
├── Specialized Modules
│   ├── OCL/                # Pure pyOpenCL implementations
│   ├── DFT/                # DFT utilities and interfaces
│   ├── GUI/                # Python GUI components
│   └── tests/              # Python-specific tests
├── Utilities
│   ├── plotUtils.py        # Visualization and plotting
│   ├── buildUtils.py       # Build system utilities
│   ├── grid_utils.py       # Grid manipulation utilities
│   └── cpp_utils_.py       # C++ integration utilities
├── External Integrations
│   ├── psi4_utils.py       # Psi4 quantum chemistry integration
│   ├── pyscf_utils.py      # PySCF integration
│   ├── dftb_utils.py       # DFTB+ integration
│   └── lammps_utils.py     # LAMMPS integration
└── Build System
    ├── Makefile_*.py       # Makefile generation scripts
    └── gen_makefile.py     # Build system generator
```

## Core Modules

### FireCore.py - Fireball DFT Interface
Main interface to the Fortran Fireball DFT code:
```python
import pyBall.FireCore as fc

# Initialize DFT calculation
fc.init_system("molecule.xyz")
fc.setup_calculation()
fc.run_scf()

# Get results
energy = fc.get_total_energy()
forces = fc.get_forces()
density = fc.get_electron_density()
```

### MMFF.py - Molecular Mechanics
Interface to C++ molecular mechanics libraries:
```python
import pyBall.MMFF as mmff

# Setup force field
world = mmff.MMFFWorld()
world.loadMolecule("molecule.xyz")
world.setupForceField()

# Run dynamics
world.run_dynamics(nsteps=1000, dt=0.001)
```

### AtomicSystem.py - Molecular Manipulation
High-level molecular structure manipulation:
```python
import pyBall.AtomicSystem as ats

# Load and manipulate molecules
mol = ats.AtomicSystem()
mol.loadXYZ("input.xyz")
mol.center_molecule()
mol.optimize_geometry()
mol.saveXYZ("output.xyz")
```

### FFfit.py - Force-Field Parameter Fitting
The FFfit modules fit transferable bonded parameters to QM Hessians. The least-norm Wilson force matrix is diagnostic only; `fit_hybrid_hessian()` instead compares the dynamical matrix in an orthonormal Wilson row space, alongside mode and local-Hessian residuals. High-level topology, subtype hierarchy, and optional stretch--stretch/stretch--bend terms are in `FFfit_utils.py`.
```python
from pyBall import FFfit

# Fit stiffness parameters to a reference Hessian
fitter = FFfit.FFfit()
fitter.set_geometry(positions)
fitter.set_symbols(['Si', 'H', 'H', 'H', 'H'])
fitter.add_bond(0, 1, 1.48)
fitter.auto_assign_types()
k = fitter.fit_lsq(H_ref)

# Graph algorithms (CSR bond-graph, bounded BFS)
dist = FFfit.FFfit.bond_graph_distances_cpp(bond_pairs, natoms)
mask = FFfit.FFfit.local_hessian_mask_cpp(bonds, natoms, max_graph_distance=2)

# Batch dihedral sensitivity (replaces Python loop)
A = fitter.dihedral_dHdk_batch_typed_cpp(dihedrals, type_idx, n_types)
```
See [`doc/Topics/FFfit/HessianFitting_Theory.md`](../doc/Topics/FFfit/HessianFitting_Theory.md) for the physical model and [`tests/tSiNCs/test_fffit_hybrid.py`](../tests/tSiNCs/test_fffit_hybrid.py) for numerical regression tests.

## Specialized Modules

### OCL/ - Pure pyOpenCL Implementations
Standalone GPU implementations independent of C++ libraries:

#### Key Files:
- **`OpenCLBase.py`** - Base OpenCL setup and utilities
- **`MMFF.py`** - OpenCL molecular mechanics
- **`GridFF.py`** - OpenCL grid-based force fields
- **`UFF.py`** - Universal Force Field implementation
- **`eFF_ocl.py`** - Electron Force Field on GPU

```python
import pyBall.OCL.MMFF as ocl_mmff

# GPU-accelerated molecular mechanics
gpu_world = ocl_mmff.MMFFWorld_OCL()
gpu_world.loadMolecule("molecule.xyz")
gpu_world.run_on_gpu(nsteps=10000)
```

### DFT/ - DFT Utilities
High-level DFT interfaces and utilities:
- **`high_level.py`** - Simplified DFT workflows
- **`utils.py`** - DFT utility functions
- **`jobs.py`** - Batch DFT job management
- **`oclfft.py`** - OpenCL FFT for DFT grids

### GUI/ - Python GUI Components
Python-based graphical interfaces:
- **`GLGUI.py`** - OpenGL-based GUI components
- **`shaders/`** - GLSL shader programs

## Utility Modules

### plotUtils.py - Visualization
Comprehensive plotting and visualization utilities:
```python
import pyBall.plotUtils as plt

# Molecular visualization
plt.plot_molecule("molecule.xyz")
plt.plot_forces(atoms, forces)
plt.plot_energy_surface(grid, energies)

# Scientific plotting
plt.plot_convergence(energies)
plt.plot_spectrum(frequencies, intensities)
```

### grid_utils.py - Grid Operations
3D grid manipulation and analysis:
```python
import pyBall.grid_utils as grid

# Grid operations
density_grid = grid.load_cube("density.cube")
potential_grid = grid.poisson_solve(density_grid)
grid.save_xsf("potential.xsf", potential_grid)
```

## External Integrations

### Quantum Chemistry Packages
- **`psi4_utils.py`** - Psi4 integration for high-level quantum calculations
- **`pyscf_utils.py`** - PySCF integration for advanced DFT methods
- **`dftb_utils.py`** - DFTB+ integration for semi-empirical methods

### Molecular Dynamics
- **`lammps_utils.py`** - LAMMPS integration for large-scale MD

## Setup and Usage

### Environment Setup
```bash
# Set Python path
export PYTHONPATH=/path/to/FireCore:$PYTHONPATH

# Install dependencies
pip install numpy matplotlib scipy
```

### Basic Usage Pattern
```python
# Import pyBall modules
import pyBall.FireCore as fc
import pyBall.MMFF as mmff
import pyBall.plotUtils as plt

# Load molecule
mol = "molecule.xyz"

# Option 1: DFT calculation
fc.init_system(mol)
energy_dft = fc.run_scf()

# Option 2: Classical force field
world = mmff.MMFFWorld()
world.loadMolecule(mol)
energy_mm = world.eval_energy()

# Visualization
plt.plot_molecule(mol)
```

## Build Integration

### C++ Library Compilation
The Python interface automatically triggers C++ library compilation when needed:
```python
# This will compile C++ libraries if needed
import pyBall.MMFF as mmff
world = mmff.MMFFWorld()  # Triggers compilation
```

### Makefile Generation
Automated build system for C++ components:
```python
import pyBall.gen_makefile as gm

# Generate makefiles for current system
gm.generate_makefiles()
gm.compile_libraries()
```

## Testing

### Running Python Tests
```bash
cd pyBall/tests
python3 -m pytest

# Or run specific tests
python3 test_mmff.py
python3 test_dft.py
```

### Integration Tests
```bash
cd tests/pyFireball
./run.sh  # Tests Python-Fortran integration

cd tests/tMMFF
python3 run.py  # Tests Python-C++ integration
```

## Development Guidelines

### Code Organization
- **Modular design** - Clear separation between interfaces and implementations
- **Consistent APIs** - Similar function signatures across modules
- **Error handling** - Proper exception handling and error messages
- **Documentation** - Docstrings for all public functions

### Performance Considerations
- **NumPy arrays** - Use NumPy for numerical data
- **Memory management** - Careful handling of large arrays
- **GPU utilization** - Leverage OpenCL for computationally intensive tasks
- **Caching** - Cache expensive computations when possible

### Integration Best Practices
- **C++ bindings** - Use ctypes or pybind11 for C++ integration
- **Fortran bindings** - Use f2py or direct library calls
- **Resource management** - Proper cleanup of C++/Fortran resources
- **Error propagation** - Translate C++/Fortran errors to Python exceptions

## Common Issues

### Import Errors
- **Missing libraries** - Ensure C++ libraries are compiled
- **Path issues** - Check PYTHONPATH environment variable
- **Dependencies** - Install required Python packages

### Performance Issues
- **Slow calculations** - Check if GPU acceleration is enabled
- **Memory usage** - Monitor memory consumption for large systems
- **Threading** - Be aware of GIL limitations for CPU-bound tasks

### Integration Problems
- **C++ compilation** - Check compiler and dependency installation
- **Fortran linking** - Verify Fortran libraries are available
- **OpenCL issues** - Check GPU drivers and OpenCL installation
