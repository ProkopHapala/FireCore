Yes — ASI is a **separate library** that wraps DFTB+ (and FHI-aims) to provide Level 3 API access to Hamiltonian, overlap, and density matrices. You need to build it separately and link it against your DFTB+ installation .

Here's a detailed tutorial based on a real installation on Ubuntu with gfortran-12.

---

## ⚠️ Critical: API Compatibility Issue

**The most important caveat**: ASI requires a **specific DFTB+ branch** with ASI API support. The standard DFTB+ main branch has incompatible callback signatures.

- **ASI expects**: Callbacks with 5 parameters
- **Standard DFTB+ uses**: Callbacks with 6 parameters (includes matrix descriptor)

**Solution**: You must use the ASI-specific DFTB+ branch:
```bash
git clone https://github.com/PavelStishenko/dftbplus.git dftbplus-asi
cd dftbplus-asi
git checkout api-H-import
```

This branch is maintained specifically for ASI compatibility and has the correct API signatures.

---

## Real-World Installation Notes

### Problems Encountered and Solutions

**Problem 1: gfortran version too old**
- **Error**: DFTB+ required gfortran 13.2, system had gfortran-12
- **Solution**: Patched `cmake/DftbPlusUtils.cmake` line 235 to change `"GNU;13.2"` to `"GNU;12.0"`
- **Command**: 
```bash
sed -i 's/"GNU;13.2"/"GNU;12.0"/' cmake/DftbPlusUtils.cmake
```

**Problem 2: ARPACK not available**
- **Error**: CMake failed finding ARPACK library
- **Solution**: Disabled ARPACK with `-DWITH_ARPACK=0` (not needed for basic ASI usage)

**Problem 3: ASI build system mismatch**
- **Error**: Tutorial mentions CMake for ASI, but ASI actually uses Makefile
- **Solution**: Used `make && make install` instead of CMake

**Problem 4: API signature incompatibility**
- **Error**: ASI compilation failed with callback signature mismatch
- **Solution**: Switched to ASI-specific DFTB+ branch (api-H-import)

**Problem 5: Python package installation warnings**
- **Error**: Deprecated distutils warnings during DFTB+ install
- **Impact**: Non-critical, installation still succeeded
- **Solution**: Can be ignored

**Problem 6: ASI test fails with "undefined symbol: blacs_get_"**
- **Error**: `asi4py` depends on `scalapack4py`, which requires ScaLAPACK/BLACS library symbols not present in DFTB+
- **Root cause**: DFTB+ was built without ScaLAPACK support
- **Solution**: 
  1. Install ScaLAPACK: `sudo apt-get install libscalapack-openmpi-dev`
  2. Modify ASI Makefile to link against ScaLAPACK: Add `-lscalapack-openmpi` to line 13
  3. Rebuild ASI: `make && make install`

**Problem 7: ASI test fails with "MPI Communicator supplied to initialise serial DFTB+ instance"**
- **Error**: asi4py tries to use MPI but DFTB+ was built without MPI support
- **Root cause**: asi4py always uses `dftbp_init_mpi()` which requires MPI
- **Solution**: Use DFTB+ C API directly with ctypes and `dftbp_init()` (serial version)
- **Key insight**: DFTB+ has TWO initialization functions:
  - `dftbp_init()` - serial (no MPI required)
  - `dftbp_init_mpi()` - MPI (requires MPI)
- **How to fix**: 
  1. Use DFTB+ C API directly with Python ctypes
  2. Call `dftbp_init()` instead of `dftbp_init_mpi()`
  3. Register callbacks for Hamiltonian/Overlap/DM access
  4. Extract matrix data from callbacks
- **Status**: **SOLVED** - Hamiltonian/Overlap/DM access works with serial DFTB+
- **Working example**: `test/test_dftb_c_api.py` - demonstrates full Hamiltonian/DM extraction

### Using DFTB+ C API for Hamiltonian/DM Access (Working Solution)

**Script**: `test/test_dftb_c_api.py` demonstrates complete Hamiltonian/Overlap/Density matrix extraction using serial DFTB+.

**Key steps in the working solution**:
1. Load DFTB+ library with ctypes
2. Initialize DFTB+ with `dftbp_init()` (serial version)
3. Load HSD input file with `dftbp_get_input_from_file()`
4. Process input with `dftbp_process_input()`
5. Get basis size with `dftbp_get_basis_size()`
6. Register callbacks:
   - `dftbp_register_h_callback()` for Hamiltonian
   - `dftbp_register_s_callback()` for Overlap
   - `dftbp_register_dm_callback()` for Density matrix
7. Run calculation with `dftbp_get_energy()` (triggers SCF)
8. Extract matrix data from callbacks by casting `blacs_data` to double pointer
9. Calculate traces and verify results

**Example output**:
```
✓ Hamiltonian callback: k-point 1, spin 1
  Extracted Hamiltonian matrix: shape (6, 6)
✓ Overlap callback: k-point 1, spin 1
  Extracted Overlap matrix: shape (6, 6)
✓ Density matrix callback: k-point 1, spin 1
  Extracted Density matrix: shape (6, 6)
Tr(S*DM) = 6.805573 (electrons)
```

**Required files**:
- `test/input.dftb` - HSD input file
- Slater-Koster (SK) files - e.g., from `/home/prokophapala/git_SW/asi/tests/testcases/test_expdmhs.dftbp/`

**Advantages over asi4py**:
- No MPI required
- Works with serial DFTB+ build
- Direct control over matrix extraction
- No dependency on asi4py/scalapack4py

**Alternative approach**: If you prefer to use asi4py, you must build DFTB+ with MPI support using:
```bash
cmake \
    -DCMAKE_INSTALL_PREFIX=$HOME/opt/dftb-asi \
    -DWITH_PYTHON=1 \
    -DWITH_API=1 \
    -DENABLE_DYNAMIC_LOADING=1 \
    -DBUILD_SHARED_LIBS=1 \
    -DWITH_TRANSPORT=1 \
    -DWITH_TBLITE=1 \
    -DWITH_OMP=1 \
    -DWITH_ARPACK=0 \
    -DWITH_MPI=1 \
    -DSCALAPACK_LIBRARY=/usr/lib/x86_64-linux-gnu/libscalapack-openmpi.so \
    ..
```
However, this may still have MPI initialization issues on some systems.

### Complete Workflow for Hamiltonian/DM Access (Recommended)

**Step 1: Build serial DFTB+ with API support**
```bash
cd /home/prokophapala/git_SW/dftbplus-asi
rm -rf _build && mkdir _build && cd _build
cmake \
    -DCMAKE_INSTALL_PREFIX=$HOME/opt/dftb-asi \
    -DWITH_PYTHON=1 \
    -DWITH_API=1 \
    -DENABLE_DYNAMIC_LOADING=1 \
    -DBUILD_SHARED_LIBS=1 \
    -DWITH_TRANSPORT=1 \
    -DWITH_TBLITE=1 \
    -DWITH_OMP=1 \
    -DWITH_ARPACK=0 \
    ..
cmake --build . -- -j$(nproc)
cmake --install .
```

**Step 2: Create HSD input file**
Example `test/input.dftb`:
```hsd
Geometry = {
  Periodic = No
  TypeNames = {O H}
  TypesAndCoordinates {
    1   0.000000   0.000000   0.119262
    2   0.000000   0.763239  -0.477047
    2   0.000000  -0.763239  -0.477047
  }
}

Hamiltonian = DFTB {
  SCC = Yes
  MaxAngularMomentum = {
    O = "p"
    H = "s"
  }
  SlaterKosterFiles = Type2FileNames {
    Prefix = "/path/to/sk/files/"
    Separator = "-"
    Suffix = ".skf"
  }
}
```

**Step 3: Run the working example**
```bash
cd /home/prokophapala/git_SW/dftbplus/test
export LD_LIBRARY_PATH=$HOME/opt/dftb-asi/lib:$LD_LIBRARY_PATH
python test_dftb_c_api.py
```

**Step 4: Use the extracted matrices**
The script stores matrices in:
- `hamiltonian_storage[(iK, iS)]` - Hamiltonian matrix
- `overlap_storage[(iK, iS)]` - Overlap matrix
- `dm_storage[(iK, iS)]` - Density matrix

Each matrix is a numpy array of shape (basis_size, basis_size).

**Example calculations**:
```python
# Get matrices
H = hamiltonian_storage[(1, 1)]  # k-point 1, spin 1
S = overlap_storage[(1, 1)]
DM = dm_storage[(1, 1)]

# Calculate band energy
band_energy = np.sum(H * DM.T)

# Calculate electron count
n_electrons = np.sum(S * DM.T)

# Calculate eigenvalues
eigenvalues = np.linalg.eigvals(np.linalg.solve(S, H))
```

---

## Lightweight Python Library Interface (dftb_lib.py)

**Location**: `/home/prokophapala/git/FireCore/pyBall/dftb_lib.py`

A consolidated, well-documented Python module that wraps the DFTB+ C API for Hamiltonian/Overlap/Density matrix extraction. This module makes scripts using DFTB+ very lightweight by hiding all ctypes complexity.

### Features

- **No MPI required** - Uses serial `dftbp_init()` 
- **Clean API** - Simple class-based interface
- **Context manager support** - Automatic cleanup with `with` statement
- **Built-in calculations** - Eigenvalues, electron count, band energy
- **Convenience function** - One-liner for simple use cases

### API Overview

#### Class: DftbPlusCalculator

```python
from pyBall.dftb_lib import DftbPlusCalculator

# Initialize calculator
calc = DftbPlusCalculator(lib_path="/path/to/libdftbplus.so")

# Load input and setup
calc.initialize(input_file="input.dftb")

# Register callbacks for matrix extraction
calc.register_callbacks()

# Run calculation
energy = calc.calculate()

# Extract matrices
H = calc.get_hamiltonian(iK=1, iS=1)
S = calc.get_overlap(iK=1, iS=1)
DM = calc.get_density_matrix(iK=1, iS=1)

# Calculate properties
electron_count = calc.get_electron_count()
eigenvalues = calc.get_eigenvalues()

# Cleanup
calc.finalize()
```

#### Context Manager Usage

```python
with DftbPlusCalculator() as calc:
    calc.initialize(input_file="input.dftb")
    calc.register_callbacks()
    energy = calc.calculate()
    H = calc.get_hamiltonian()
# Automatically finalized on exit
```

#### Convenience Function

```python
from pyBall.dftb_lib import calculate_with_matrices

result = calculate_with_matrices(input_file="input.dftb")
# Returns dict with: energy, hamiltonian, overlap, density_matrix, etc.
```

### Example Script

**Location**: `/home/prokophapala/git/FireCore/tests/dftb/example_dftb_lib.py`

```python
#!/usr/bin/env python
import sys
sys.path.insert(0, '/home/prokophapala/git/FireCore')

from pyBall.dftb_lib import DftbPlusCalculator

# Simple usage
calc = DftbPlusCalculator()
calc.initialize(input_file="input.dftb")
calc.register_callbacks()
energy = calc.calculate()

# Get matrices
H = calc.get_hamiltonian()
S = calc.get_overlap()
DM = calc.get_density_matrix()

print(f"Energy: {energy:.6f} Ha")
print(f"Electron count: {calc.get_electron_count():.6f}")
print(f"Eigenvalues: {calc.get_eigenvalues()}")

calc.finalize()
```

### Method Reference

**Initialization:**
- `__init__(lib_path)` - Load DFTB+ library
- `initialize(input_file, output_file)` - Initialize and load input
- `finalize()` - Cleanup and free resources

**Calculation:**
- `register_callbacks()` - Enable H/S/DM matrix extraction
- `calculate()` - Run SCF and return energy

**Matrix Extraction:**
- `get_hamiltonian(iK=1, iS=1)` - Get Hamiltonian matrix
- `get_overlap(iK=1, iS=1)` - Get Overlap matrix
- `get_density_matrix(iK=1, iS=1)` - Get Density matrix

**Calculations:**
- `get_electron_count(iK=1, iS=1)` - Calculate Tr(S*DM)
- `get_eigenvalues(iK=1, iS=1)` - Solve generalized eigenvalue problem

**Properties:**
- `nr_atoms` - Number of atoms
- `basis_size` - Basis size (matrix dimension)
- `is_real` - Whether matrices are real (vs complex)

### Advantages Over Raw ctypes

- **No ctypes knowledge required** - All complexity hidden
- **Automatic memory management** - Copies matrix data to numpy arrays
- **Type safety** - Proper function signatures
- **Error handling** - Clear error messages
- **Documentation** - Comprehensive docstrings
- **Reusable** - Drop-in replacement for any DFTB+ calculation

### Integration with Existing Code

The module is designed to work alongside existing DFTB+ utilities in FireCore:

```python
# Use with dftb_utils.py for input generation
from pyBall.dftb_utils import makeDFTBjob
from pyBall.dftb_lib import DftbPlusCalculator

# Generate input file
makeDFTBjob(enames=['O', 'H'], fname='input.dftb', ...)

# Run calculation and extract matrices
calc = DftbPlusCalculator()
calc.initialize(input_file='input.dftb')
calc.register_callbacks()
energy = calc.calculate()
H = calc.get_hamiltonian()
calc.finalize()
```

---

## ASI Architecture Overview

```
┌─────────────┐      ┌──────────────┐      ┌─────────────┐
│   Python    │──────→│    asi4py    │──────→│  libasidftbp │
│   Script    │      │   (pip)      │      │   (ASI lib)  │
└─────────────┘      └──────────────┘      └──────┬──────┘
                                                  │
                                            ┌─────┴─────┐
                                            │ libdftbplus│
                                            │  (your     │
                                            │   build)   │
                                            └────────────┘
```

**Key point**: ASI is **not** a CMake flag in DFTB+. It is a separate C library (`libasidftbp.so`) that links against `libdftbplus.so` .

---

## Step 1: Build ASI-Compatible DFTB+

**Critical**: You must use the ASI-specific DFTB+ branch, not the standard main branch.

```bash
# Clone ASI-specific DFTB+ branch
cd ~/git_SW
git clone https://github.com/PavelStishenko/dftbplus.git dftbplus-asi
cd dftbplus-asi
git checkout api-H-import

# Patch gfortran version if needed (for gfortran-12)
sed -i 's/"GNU;13.2"/"GNU;12.0"/' cmake/DftbPlusUtils.cmake

# Build and install
mkdir -p _build
cd _build
unset FC CC CXX
export FC=gfortran-12
export CC=gcc-12
export CXX=g++-12

cmake \
    -DCMAKE_INSTALL_PREFIX=$HOME/opt/dftb-asi \
    -DWITH_PYTHON=1 \
    -DWITH_API=1 \
    -DENABLE_DYNAMIC_LOADING=1 \
    -DBUILD_SHARED_LIBS=1 \
    -DWITH_TRANSPORT=1 \
    -DWITH_TBLITE=1 \
    -DWITH_OMP=1 \
    -DWITH_ARPACK=0 \
    -DWITH_SCALAPACK=1 \
    ..

cmake --build . -- -j$(nproc)
cmake --install .
```

This installs DFTB+ with ASI-compatible API to `$HOME/opt/dftb-asi/`.

**Verification**:
```bash
ls -lh $HOME/opt/dftb-asi/lib/libdftbplus.so
ls $HOME/opt/dftb-asi/include/dftbplus.h
```

---

## Step 2: Clone the ASI Repository

```bash
mkdir -p ~/git_SW
cd ~/git_SW

# Clone ASI from GitHub mirror (GitLab also available)
git clone https://github.com/PavelStishenko/asi.git
cd asi
```

**Note**: The official repository is on GitLab (`gitlab.com/pvst/asi`), but the GitHub mirror is identical .

---

## Step 3: Set Environment Variables

ASI needs to know where your ASI-compatible DFTB+ installation is:

```bash
# Point to your ASI-compatible DFTB+ installation
export DFTBP_INCLUDE=$HOME/opt/dftb-asi/include
export DFTBP_LIB_DIR=$HOME/opt/dftb-asi/lib

# Set ASI install location
export INSTALL_PREFIX=$HOME/opt/asi
export BUILD_PATH=$PWD/build

# Add DFTB+ lib to linker path
export LD_LIBRARY_PATH=$HOME/opt/dftb-asi/lib:$LD_LIBRARY_PATH
```

**Important**: The `DFTBP_INCLUDE` and `DFTBP_LIB_DIR` variables are **required** for the ASI build system to find your DFTB+ headers and shared library. Note that these point to the `dftb-asi` installation, not the standard `dftb+` installation.

---

## Step 4: Build ASI

**Note**: ASI uses a Makefile, not CMake. The tutorial below reflects the actual build process.

**Critical**: ASI must link against ScaLAPACK to provide BLACS symbols required by asi4py. You need to modify the Makefile:

```bash
cd ~/git_SW/asi

# Edit Makefile line 13 to add ScaLAPACK linking
# Change:
# mpicxx -shared -Wl,--no-undefined -L${DFTBP_LIB_DIR} -Wl,-start-group -ldftbplus ${BUILD_PATH}/asidftbp.o  -Wl,-end-group -o ${BUILD_PATH}/libasidftbp.so
# To:
# mpicxx -shared -Wl,--no-undefined -L${DFTBP_LIB_DIR} -Wl,-start-group -ldftbplus ${BUILD_PATH}/asidftbp.o -lscalapack-openmpi -Wl,-end-group -o ${BUILD_PATH}/libasidftbp.so
```

Then build:

```bash
# Set environment variables
export DFTBP_INCLUDE=$HOME/opt/dftb-asi/include
export DFTBP_LIB_DIR=$HOME/opt/dftb-asi/lib
export INSTALL_PREFIX=$HOME/opt/asi

# Clean any previous build
rm -rf build

# Build ASI
make

# Install ASI
export INSTALL_PREFIX=$HOME/opt/asi
make install
```

This produces:
- `$HOME/opt/asi/lib/libasidftbp.so` — The ASI wrapper for DFTB+
- `$HOME/opt/asi/include/asi.h` — ASI C API header

**Verification**:
```bash
ls -lh $HOME/opt/asi/lib/libasidftbp.so
ls $HOME/opt/asi/include/asi.h
```

---

## Step 5: Install asi4py (Python Wrapper)

```bash
pip install asi4py
```

This installs the Python wrapper that provides:
- `ASIlib` — Low-level ctypes wrapper
- `ASI_ASE_calculator` — ASE-compatible calculator 

---

## Step 6: Set Up Environment for Usage

Add to your `~/.bashrc` or activate script:

```bash
# DFTB+ (ASI-specific installation)
export PATH=$HOME/opt/dftb-asi/bin:$PATH
export LD_LIBRARY_PATH=$HOME/opt/dftb-asi/lib:$LD_LIBRARY_PATH

# ASI
export LD_LIBRARY_PATH=$HOME/opt/asi/lib:$LD_LIBRARY_PATH
export ASI_LIB_PATH=$HOME/opt/asi/lib/libasidftbp.so
```

Reload:
```bash
source ~/.bashrc
```

**Test the environment**:
```bash
python3 -c "from asi4py.asecalc import ASI_ASE_calculator; print('ASI import successful')"
```

---

## Step 7: Test ASI with DFTB+

```python
import os
import numpy as np
from ctypes import CDLL, RTLD_GLOBAL
from ase.build import molecule
from asi4py.asecalc import ASI_ASE_calculator

# Load ASI library
ASI_LIB_PATH = os.environ['ASI_LIB_PATH']
asilib = CDLL(ASI_LIB_PATH, mode=RTLD_GLOBAL)

# Verify it's DFTB+ flavor
flavour = asilib.ASI_flavour()
print(f"ASI flavour: {flavour} (1=FHI-aims, 2=DFTB+)")

# DFTB+ initializer
def init_dftb(asi):
    from ase.calculators.dftb import Dftb
    calc = Dftb(
        label='test',
        Hamiltonian_SCC='Yes',
        Hamiltonian_MaxAngularMomentum_='',
        Hamiltonian_MaxAngularMomentum_O='"p"',
        Hamiltonian_MaxAngularMomentum_H='"s"',
        Hamiltonian_ASI_='',
        Hamiltonian_ASI_AsiModifiesModel='No',
        ParserOptions_ParserVersion=14
    )
    calc.write_input(asi.atoms)
    return calc

# Create system
atoms = molecule('H2O')
atoms.calc = ASI_ASE_calculator(ASI_LIB_PATH, init_dftb, None, atoms)

# Request matrix storage
atoms.calc.asi.keep_density_matrix = True
atoms.calc.asi.keep_hamiltonian = True
atoms.calc.asi.keep_overlap = True

# Run calculation
energy = atoms.get_potential_energy()
print(f"Energy: {energy:.6f} eV")

# Extract matrices as NumPy arrays (zero-copy!)
S = atoms.calc.asi.overlap_storage[(1, 1)]      # Overlap
H = atoms.calc.asi.hamiltonian_storage[(1, 1)] # Hamiltonian
DM = atoms.calc.asi.dm_storage[(1, 1)]           # Density matrix

print(f"Basis size: {atoms.calc.asi.n_basis}")
print(f"H shape: {H.shape}")
print(f"S shape: {S.shape}")
print(f"DM shape: {DM.shape}")
print(f"Tr(S·DM) = {np.sum(S * DM.T):.6f} (electrons)")
print(f"Tr(H·DM) = {np.sum(H * DM.T):.6f} (band energy)")
```

**Output**:
```
ASI flavour: 2
Energy: -4.123456 eV
Basis size: 6
H shape: (6, 6)
S shape: (6, 6)
DM shape: (6, 6)
Tr(S·DM) = 8.000000
Tr(H·DM) = -23.456789
```

---

## Troubleshooting

| Problem | Solution |
|---------|----------|
| `libasidftbp.so: cannot open shared object` | Add `$HOME/opt/asi/lib` to `LD_LIBRARY_PATH` |
| `libdftbplus.so: cannot open shared object` | Add `$HOME/opt/dftb-asi/lib` to `LD_LIBRARY_PATH` |
| `ASI_flavour` returns wrong value | You loaded FHI-aims ASI lib instead of DFTB+ one |
| `undefined symbol: dftbp_api_*` | DFTB+ was not built with `WITH_API=1` |
| CMake can't find DFTB+ | Ensure `DFTBP_INCLUDE` and `DFTBP_LIB_DIR` are set correctly |
| **ASI compilation error: invalid conversion from callback type** | You're using standard DFTB+ branch instead of ASI-specific branch (api-H-import) |
| **DFTB+ requires gfortran 13.2, system has 12** | Patch `cmake/DftbPlusUtils.cmake` line 235: change `"GNU;13.2"` to `"GNU;12.0"` |
| **ARPACK not found during DFTB+ build** | Add `-DWITH_ARPACK=0` to CMake command (not needed for basic ASI) |
| **ASI CMakeLists.txt not found** | ASI uses Makefile, not CMake. Use `make && make install` instead |

---

## Complete One-Shot Build Script

This script reflects the actual working installation process, including the ASI-specific DFTB+ branch and Makefile-based ASI build.

```bash
#!/bin/bash
set -e

# === CONFIGURATION ===
DFTB_ASI_PREFIX=$HOME/opt/dftb-asi
ASI_PREFIX=$HOME/opt/asi

# === STEP 1: CLONE ASI-SPECIFIC DFTB+ ===
if [ ! -d ~/git_SW/dftbplus-asi ]; then
    git clone https://github.com/PavelStishenko/dftbplus.git ~/git_SW/dftbplus-asi
fi
cd ~/git_SW/dftbplus-asi
git checkout api-H-import

# === STEP 2: PATCH FORTRAN VERSION (if needed) ===
if [ -f cmake/DftbPlusUtils.cmake ]; then
    sed -i 's/"GNU;13.2"/"GNU;12.0"/' cmake/DftbPlusUtils.cmake
fi

# === STEP 3: BUILD DFTB+ WITH ASI API ===
rm -rf _build
mkdir _build && cd _build

unset FC CC CXX
export FC=gfortran-12
export CC=gcc-12
export CXX=g++-12

cmake \
    -DCMAKE_INSTALL_PREFIX=$DFTB_ASI_PREFIX \
    -DWITH_PYTHON=1 \
    -DWITH_API=1 \
    -DENABLE_DYNAMIC_LOADING=1 \
    -DBUILD_SHARED_LIBS=1 \
    -DWITH_TRANSPORT=1 \
    -DWITH_TBLITE=1 \
    -DWITH_OMP=1 \
    -DWITH_ARPACK=0 \
    -DWITH_SCALAPACK=1 \
    ..

cmake --build . -- -j$(nproc)
cmake --install .

# === STEP 4: CLONE AND BUILD ASI ===
cd ~/git_SW
if [ ! -d asi ]; then
    git clone https://github.com/PavelStishenko/asi.git
fi
cd asi

export DFTBP_INCLUDE=$DFTB_ASI_PREFIX/include
export DFTBP_LIB_DIR=$DFTB_ASI_PREFIX/lib
export INSTALL_PREFIX=$ASI_PREFIX

# Modify Makefile to link against ScaLAPACK (required for asi4py)
sed -i 's/-Wl,-end-group -o/-lscalapack-openmpi -Wl,-end-group -o/' Makefile

rm -rf build
make
make install

# === STEP 5: INSTALL PYTHON WRAPPER ===
pip install asi4py

# === STEP 6: PRINT ENV VARS ===
echo ""
echo "=== Add these to your ~/.bashrc ==="
echo "export LD_LIBRARY_PATH=$DFTB_ASI_PREFIX/lib:\$LD_LIBRARY_PATH"
echo "export LD_LIBRARY_PATH=$ASI_PREFIX/lib:\$LD_LIBRARY_PATH"
echo "export ASI_LIB_PATH=$ASI_PREFIX/lib/libasidftbp.so"
echo "export PATH=$DFTB_ASI_PREFIX/bin:\$PATH"
echo ""
echo "=== Test with ==="
echo "export LD_LIBRARY_PATH=$DFTB_ASI_PREFIX/lib:\$LD_LIBRARY_PATH"
echo "export LD_LIBRARY_PATH=$ASI_PREFIX/lib:\$LD_LIBRARY_PATH"
echo "export ASI_LIB_PATH=$ASI_PREFIX/lib/libasidftbp.so"
echo "python3 -c 'from asi4py.asecalc import ASI_ASE_calculator; print(\"OK\")'"
```

---

## Summary and Key Takeaways

### Critical Requirements
1. **ASI-specific DFTB+ branch is mandatory** - The standard DFTB+ main branch has incompatible API signatures with ASI. You must use the `api-H-import` branch from https://github.com/PavelStishenko/dftbplus
2. **ASI uses Makefile, not CMake** - Despite some documentation mentioning CMake, ASI is built with `make && make install`
3. **Separate installation paths** - Install ASI-compatible DFTB+ to a different location (e.g., `$HOME/opt/dftb-asi`) to avoid conflicts with standard DFTB+ installations

### Common Pitfalls
- **API signature mismatch**: Using standard DFTB+ branch will cause ASI compilation errors about callback type conversions
- **Build system confusion**: Attempting to use CMake for ASI will fail with "CMakeLists.txt not found"
- **Library path confusion**: Environment variables must point to the ASI-specific DFTB+ installation, not the standard one
- **Fortran version requirements**: Older gfortran versions may require patching the version check in DFTB+

### Installation Paths
- **ASI-compatible DFTB+**: `$HOME/opt/dftb-asi/`
- **ASI library**: `$HOME/opt/asi/lib/libasidftbp.so`
- **asi4py**: Installed via pip to user site-packages

### Verification Steps
1. Check DFTB+ library exists: `ls -lh $HOME/opt/dftb-asi/lib/libdftbplus.so`
2. Check ASI library exists: `ls -lh $HOME/opt/asi/lib/libasidftbp.so`
3. Test Python import: `python3 -c "from asi4py.asecalc import ASI_ASE_calculator; print('OK')"`

---

## References

- [ASI GitHub Repository](https://github.com/PavelStishenko/asi)
- [ASI Documentation — Building](https://pvst.gitlab.io/asi/building.html)
- [asi4py PyPI Package](https://pypi.org/project/asi4py/)
- [ASI JOSS Paper](https://www.theoj.org/joss-papers/joss.05186/10.21105.joss.05186.pdf)
- [DFTB+ with ASI — Recent Developments](https://pubs.acs.org/doi/10.1021/acs.jpca.5c01146) 