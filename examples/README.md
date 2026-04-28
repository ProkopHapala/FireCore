
# FireCore Tests and Examples

This directory contains examples and test scripts that demonstrate FireCore functionality. **This is the best place to understand what is implemented and how to use it.**

## 🚀 Quick Start

All test directories contain `run.sh` scripts that automatically compile and run the programs:

```bash
cd tests/tMMFF        # Molecular mechanics
./run.sh

cd tests/tMolGUIapp   # GUI applications
./run.sh

cd tests/Fireball/t02_CH4  # DFT calculations
./run.sh
```

## Test Categories

### Molecular Mechanics and Force Fields
- **`tMMFF/`** - Basic molecular mechanics force field tests using python-bining to C++ library
- **`tMMFFmulti/`** - Multi-molecule force field tests using python-bining to C++ library with OpenCL
- **`tMMFFsp3/`** - sp3 hybridization specific tests
- **`tUFF/`** - Universal Force Field tests
- **`tFitFF/`** - Force field parameter fitting
- **`tFitREQ/`** - non-covalent interaction parameter fitting (REQH = {Rvdw,Evdw,Q,Hbond})

### GUI Applications
- **`tMolGUIapp/`** - Basic molecular GUI application
- **`tMolGUIapp_multi/`** - Multi-molecule GUI tests
- **`tMolGUIapp_QMMM/`** - QM/MM hybrid method GUI
- **`tMolGUIapp_QMMM_multi/`** - Multi-molecule QM/MM

### DFT and Quantum Methods
- **`Fireball/`** - Fireball DFT tests (H2, CH4, pentacene, etc.)
- **`tDFT/`** - General DFT calculation tests
- **`tDFT_CO/`** - Carbon monoxide DFT tests
- **`tDFT_pentacene/`** - Pentacene molecule DFT tests
- **`pyFireball/`** - Python interface to Fireball DFT
- **`pySCF/`** - PySCF integration tests
- **`dftb/`** - DFTB+ integration tests

### Specialized Methods
- **`tEFF/`** - Electron Force Field tests
- **`tEFFapp/`** - EFF application tests
- **`tSchroedinger1D/`** - 1D Schrödinger equation solver
- **`tSchroedinger2D/`** - 2D Schrödinger equation solver
- **`tQuadrature/`** - Numerical integration tests

### GPU and Parallel Computing
- **`tCUDA/`** - CUDA acceleration tests
- **`NonBondSampling/`** - Non-bonded interaction sampling

### Molecular Tools and Utilities
- **`tAttach/`** - Molecular attachment and building tools
- **`tKekule/`** - Kekulé structure generation
- **`tLattice2D/`** - 2D lattice generation
- **`tFF2D/`** - 2D force field tests
- **`tLammpsTrj/`** - LAMMPS trajectory analysis
- **`tPsi4resp/`** - Psi4 RESP charge fitting

### Specialized Applications
- **`tQMMM_diacetylene/`** - QM/MM study of diacetylene
- **`blender/`** - Blender integration scripts
- **`pyutils/`** - Python utility scripts

## File Structure in Test Directories

Most test directories contain:
- **`run.sh`** - Main execution script (compiles and runs)
- **`run.py`** - Python execution script (if applicable)
- **Input files** - Molecular structures (.xyz, .mol2), parameters, etc.
- **`data/` or `common_resources/`** - Symbolic links to shared resources
- **Output files** - Results, logs, trajectories (generated after running)

## Usage Patterns

### C++ Applications:
```bash
cd tests/tMolGUIapp
./run.sh                    # Compiles and runs GUI application
./run.sh no                 # Runs without recompiling
```

### Python Scripts:
```bash
cd tests/tMMFF
python3 run.py              # Direct Python execution
./run.sh                    # May also compile C++ libraries first
```

### Fortran DFT:
```bash
cd tests/Fireball/t02_CH4
./run.sh                    # Requires Fortran compilation
```

## Important Notes

1. **Always use `run.sh` scripts** - they handle compilation, paths, and arguments correctly
2. **Check for symbolic links** - many tests create `data/` and `common_resources/` links
3. **Resource requirements** - Some tests require `Fdata_HC_minimal` (download from fireball-qmd.github.io)
4. **GPU tests** - CUDA/OpenCL tests require appropriate hardware and drivers

## Troubleshooting

- **Missing resources**: Check if symbolic links to `common_resources` exist
- **Compilation errors**: Ensure dependencies are installed (see main README.md)
- **Fortran tests failing**: Verify `Fdata_HC_minimal` is available
- **GPU tests failing**: Check OpenCL/CUDA installation and hardware support

## Adding New Tests

When adding new tests:
1. Create a new directory following naming convention (`t<TestName>/`)
2. Include a `run.sh` script that handles compilation and execution
3. Add symbolic links to required resources
4. Document the test purpose and expected outputs
5. Update this README.md with the new test category
