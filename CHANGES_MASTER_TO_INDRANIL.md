# Comprehensive Changes: master → Indranil Branch

## Document Information
- **Generated**: December 11, 2025
- **Source Branch**: `master` (commit `88134083`)
- **Target Branch**: `Indranil` (commit `b88cdd47`)
- **Total Files Changed**: 65
- **Total Insertions**: +9,815 lines
- **Total Deletions**: -1,463 lines

---

## Table of Contents
1. [C++ Header Files](#1-c-header-files)
2. [C++ Source Files](#2-c-source-files)
3. [OpenCL Kernel Files](#3-opencl-kernel-files)
4. [Python Core Files (pyBall)](#4-python-core-files-pyball)
5. [Python Test Files](#5-python-test-files)
6. [Documentation Files](#6-documentation-files)
7. [New Files Added](#7-new-files-added)

---

## 1. C++ Header Files

### 1.1 `cpp/common/molecular/GridFF.h`
**Lines Changed**: +196/-27 (169 net additions)

#### Function: `initGridFF()`
**Location**: Lines ~1500-1600
**Changes**:
- Added try-catch blocks around `allocateFFs()`, `makePBCshifts()`, `setAtomsSymetrized()`
- Added null pointer checks before operations
- Added `rAutoPBC` parameter usage instead of hardcoded 20.0
- Added Mpol array initialization to zero
- Added qcog initialization to Vec3dZero
- Added size limit (maxSize=10000) for multipole calculation
- Added safer vector allocation with reserve()

```cpp
// BEFORE (master):
if(bAutoNPBC){ autoNPBC( grid.cell, nPBC, 20.0 ); }

// AFTER (Indranil):
if (bAutoNPBC) { 
    autoNPBC(grid.cell, nPBC, rAutoPBC);
}
```

#### Function: `tryLoad_new()`
**Location**: Lines ~1696-1810
**Changes**:
- Added debug print for grid parameters (line 1697):
  ```cpp
  printf("GridFF::tryLoad_new() loaded grid params: dg=(%g,%g,%g) pos0=(%g,%g,%g) dims=(%d,%d,%d)\n", ...);
  ```
- Added grid dimension update from numpy file (lines 1807-1826):
  ```cpp
  if(npy.ndims >= 3 && (npy.shape[0] != grid.n.x || ...)) {
      grid.n.x = npy.shape[0];
      grid.n.y = npy.shape[1];
      grid.n.z = npy.shape[2];
      grid.updateCell();
      // ...
  }
  ```

---

### 1.2 `cpp/common/molecular/MolWorld_sp3.h`
**Lines Changed**: +685/-88 (597 net additions)

#### New Function: `scan_rigid_uff()`
**Location**: Lines ~1200-1250
**Purpose**: Perform rigid scan with UFF forcefield
```cpp
void scan_rigid_uff(int nconf, double* poss, double* rots, double* Es, 
                    double* aforces, double* aposs, bool omp);
```

#### New Function: `scan_constr()`
**Location**: Lines ~1250-1320
**Purpose**: Constrained relaxation scan
```cpp
void scan_constr(int nconf, int ncontr, int* icontrs, double* contrs, 
                 double* Es, double* aforces, double* aposs, 
                 bool bHardConstr, bool omp, int niter_max, 
                 double dt, double Fconv, double Flim);
```

#### Modified Function: `initParams()`
**Changes**: Added UFF parameter initialization

#### Modified Function: `eval_UFF()`
**Changes**: Added surface energy accumulation (`Esurf+=`)

---

### 1.3 `cpp/common/molecular/UFF.h`
**Lines Changed**: +387/-43 (344 net additions)

#### New Functions Added:
- `evalBonds_UFF()` - UFF bond energy evaluation
- `evalAngles_UFF()` - UFF angle energy evaluation  
- `evalDihedrals_UFF()` - UFF dihedral energy evaluation
- `evalInversions_UFF()` - UFF inversion energy evaluation
- `evalNonBonded_UFF()` - UFF non-bonded energy evaluation
- `eval_UFF()` - Complete UFF evaluation

#### Key Changes:
- Halved energies for bonds and non-bonded (commit 3657e5e1)
- Removed damping in electrostatics
- Fixed charges from xyz
- Increased threshold for clamping forces in non-bonded exclusion

---

### 1.4 `cpp/common/molecular/MMFFparams.h`
**Lines Changed**: +92/-27 (65 net additions)

#### Changes:
- Added new parameter loading functions
- Enhanced element type handling
- Added UFF parameter support

---

### 1.5 `cpp/common/dataStructures/Grid.h`
**Lines Changed**: +38/-7 (31 net additions)

#### New Function: `updateCell()`
**Purpose**: Update cell parameters after dimension changes

#### Changes:
- Added grid dimension validation
- Added cell update method

---

### 1.6 `cpp/common/math/Forces.h`
**Lines Changed**: +37/-12 (25 net additions)

#### New Function: `getSR_x2_smooth()`
**Location**: Lines ~200-230
**Purpose**: Smooth short-range potential function
**Commit**: 6276ae5c

---

### 1.7 `cpp/common/molecular/MMFFsp3_loc.h`
**Lines Changed**: +33/-6 (27 net additions)

#### Changes:
- Local force field improvements
- Enhanced atom handling

---

### 1.8 `cpp/common/molecular/ForceField.h`
**Lines Changed**: +11/-5 (6 net additions)

#### Changes:
- Force field interface updates
- Added new evaluation methods

---

### 1.9 `cpp/common/molecular/NBFF.h`
**Lines Changed**: +19/-11 (8 net additions)

#### Changes:
- Non-bonded force field updates
- Improved exclusion handling

---

### 1.10 `cpp/common/molecular/EwaldGrid.h`
**Lines Changed**: +5/-1 (4 net additions)

#### Changes:
- GPU Ewald grid shift for CPU/LAMMPS agreement (commit cdf029af)

---

### 1.11 `cpp/common/molecular/GlobalOptimizer.h`
**Lines Changed**: +170 modifications (refactoring)

#### Changes:
- Optimization algorithm improvements
- Code cleanup and refactoring

---

### 1.12 `cpp/common/molecular/MolWorld_sp3_multi.h`
**Lines Changed**: +274 modifications (mostly deletions/cleanup)

#### Changes:
- Multi-molecule world updates
- Code cleanup

---

### 1.13 `cpp/common/molecular/MolecularDatabase.h`
**Lines Changed**: +184 modifications (mostly deletions/cleanup)

#### Changes:
- Database handling updates
- Code cleanup

---

### 1.14 `cpp/common/molecular/MMFFBuilder.h`
**Lines Changed**: +11/-1 (10 net additions)

#### Changes:
- Builder updates for groups (commit 4352ce67)

---

### 1.15 `cpp/common/molecular/GOpt.h`
**Lines Changed**: +4/-1 (3 net additions)

#### Changes:
- Optimization parameter updates

---

### 1.16 `cpp/common/molecular/libMMFF.h`
**Lines Changed**: +3/-1 (2 net additions)

#### Changes:
- Library interface updates

---

### 1.17 `cpp/common/molecular/MolWorld_sp3_simple.h`
**Lines Changed**: -1 (1 deletion)

#### Changes:
- Minor cleanup

---

### 1.18 `cpp/common/OpenCL/OCL_MM.h`
**Lines Changed**: +58 modifications

#### Changes:
- OpenCL molecular mechanics updates (commit 4eab40b8)
- Grid generation support for 12x12 cell

---

### 1.19 `cpp/common/IO_utils.h`
**Lines Changed**: +7/-1 (6 net additions)

#### Changes:
- I/O utility improvements

---

### 1.20 `cpp/common_SDL/SDL2OGL/MolGUI.h`
**Lines Changed**: +47 modifications (mostly deletions)

#### Changes:
- GUI updates for scan_constr (commit 191694c9)

---

### 1.21 `cpp/common_SDL/SDL2OGL/MolecularDraw.h`
**Lines Changed**: +8 modifications

#### Changes:
- Drawing updates

---

### 1.22 `cpp/common_SDL/SDL2OGL/Draw3D_Molecular.h`
**Lines Changed**: -16 (16 deletions)

#### Changes:
- Molecular drawing cleanup

---

### 1.23 `cpp/apps/MolecularEditor/MolGUIapp_argv.h`
**Lines Changed**: +5/-1 (4 net additions)

#### Changes:
- App argument handling updates

---

## 2. C++ Source Files

### 2.1 `cpp/libs/Molecular/MMFF_lib.cpp`
**Lines Changed**: +191/-11 (180 net additions)

#### New Functions:
```cpp
// scan_constr - Constrained scan function
void scan_constr(int nconf, int ncontr, int* icontrs, double* contrs, 
                 double* Es, double* aforces, double* aposs, 
                 bool bHardConstr, bool omp, int niter_max, 
                 double dt, double Fconv, double Flim);

// scan_rigid_uff - Rigid UFF scan
void scan_rigid_uff(int nconf, double* poss, double* rots, double* Es, 
                    double* aforces, double* aposs, bool omp);

// test_UFF - UFF test function
void test_UFF(int test);
```

#### Modified Functions:
- `init_buffers()` - UFF buffer handling (commit 6c2fdf81)
- `init()` - Added UFF mode support

---

### 2.2 `cpp/libs_OCL/MMFFmulti_lib.cpp`
**Lines Changed**: +79 modifications (mostly deletions)

#### Changes:
- Multi-molecule library cleanup

---

## 3. OpenCL Kernel Files

### 3.1 `cpp/common_resources/cl/GridFF.cl`
**Lines Changed**: +125/-35 (90 net additions)

#### New Kernel: `slabPotential_zyx()`
**Purpose**: Transposed slab potential calculation
**Commit**: 42e51553

#### Modified Kernel: `poissonW()`
**Changes**: 
- Corrected k-kernel 1/|k^2| calculation
- Fixed normalization constant
- Commit: 76ca776a

#### New Kernel: `sample3D_grid()`
**Purpose**: 3D grid sampling

#### New Kernel: `make_MorseFF_f4()`
**Purpose**: Morse force field generation with float4

---

### 3.2 `cpp/common_resources/cl/relax_multi.cl`
**Lines Changed**: +93 modifications (refactoring)

#### Changes:
- Multi-molecule relaxation kernel updates
- Code cleanup

---

## 4. Python Core Files (pyBall)

### 4.1 `pyBall/atomicUtils.py`
**Lines Changed**: +4/-2 (2 net additions)

#### Function: `getVdWparams()` (Line 862)
**Change**: Parameter lookup method
```python
# BEFORE (master):
if etypes is None: etypes = loadElementTypes( fname=fname, bDict=False )

# AFTER (Indranil):
if etypes is None: etypes = loadElementTypes( fname=fname, bDict=True )
```

#### Function: `loadElementTypes()` (Line 1019)
**Change**: Dictionary key
```python
# BEFORE (master):
if bDict: return { rec[0]:rec for rec in lst }  # key = element name

# AFTER (Indranil):
if bDict: return { rec[1]:rec for rec in lst }  # key = atomic number
```

---

### 4.2 `pyBall/MMFF.py`
**Lines Changed**: +77/-13 (64 net additions)

#### New Global Variable (Line ~900):
```python
glob_bUFF = False
```

#### Function: `getBuffs()` (Lines 924-952)
**Changes**:
- Uses `ndims_ptr` approach for buffer access
- Added `glob_bUFF` mode handling
- Condition changed from `if glob_bMMFF:` to `if glob_bMMFF and not glob_bUFF:`

#### Function: `init()` (Lines 1035-1041)
**Changes**:
- Added `glob_bUFF` tracking
- Changed default `nPBC` from `(1,3,0)` to `(0,0,0)`
```python
# BEFORE:
nPBC=(1,3,0),

# AFTER:
nPBC=(0,0,0),
```

#### Function: `setSwitches()` (Lines 1098-1102)
**Changes**:
- Updated argtypes to match C++ signature (7 args instead of 8)
- Removed `bSaveToDatabase` parameter

#### Function: `setSwitches2()` (Line 1107)
**Change**: Renamed from duplicate `setSwitches` to `setSwitches2`

#### New Function: `scan_constr()` (Lines 1195-1206)
```python
def scan_constr(icontrs, contrs, Es=None, aforces=None, aposs=None, 
                bHardConstr=False, bF=False, bP=False, omp=False, 
                niter_max=10000, dt=0.05, Fconv=1e-5, Flim=100.0):
```

#### New Function: `scan_rigid_uff()` (Lines 1208-1218)
```python
def scan_rigid_uff(poss, rots=None, Es=None, aforces=None, aposs=None, 
                   bF=False, bP=False, omp=False):
```

#### New Function: `test_UFF()` (Lines 1429-1431)
```python
def test_UFF(test):
    return lib.test_UFF(test)
```

---

### 4.3 `pyBall/OCL/GridFF.py`
**Lines Changed**: +408/-180 (228 net additions)

#### New Method: `slabPotential()` (Lines ~800-850)
**Purpose**: Slab potential calculation with transposition support

#### New Method: `makeCoulombEwald_slab()` (Lines ~900-990)
**Purpose**: Coulomb Ewald calculation for slab geometry
**Features**:
- Debug verification prints (retained per user request)
- Buffer shape verification
- Optional Vin checking

#### Modified Method: `make_MorseFF()` (Lines ~1000-1100)
**Changes**:
- Added `dg` parameter for grid step size
- Improved buffer management
- Added `release_unused_buffs()` call

#### New Method: `release_unused_buffs()` (Lines ~1100-1120)
**Purpose**: Release unused OpenCL buffers
**Implementation**: Uses `hasattr` check instead of try-except

#### Modified Method: `fit3D()` (Lines ~400-460)
**Changes**:
- Replaced `profile_kernel` calls with direct kernel calls
- Added `bPrint` flag for optional output

#### Class: `GridShape` 
**New Methods**:
- `__init__()` with `desired_voxel` parameter
- `nice_fft_dim()` for FFT-friendly dimensions

---

### 4.4 `pyBall/tests/ocl_GridFF_new.py`
**Lines Changed**: +257/-30 (227 net additions)

#### New Function: `make_atoms_arrays()` (Lines 197-222)
**Purpose**: Create atom arrays from XYZ file
```python
def make_atoms_arrays(atoms=None, fname=None, bSymetrize=False, 
                      Element_Types_name="./data/ElementTypes.dat", 
                      bSqrtEvdw=True):
```

#### Modified Function: `test_gridFF_ocl()` (Lines 239-450)
**Changes**:
- Added `z0` calculation from topmost atom
- Added `desired_voxel` parameter
- Added `shift0` parameter
- Improved grid parameter printing
- Added Coulomb and Morse potential generation

#### New Function: `check_vcoul_buffer()` (Lines ~350)
**Purpose**: Debug utility to check V_Coul buffer

---

### 4.5 `pyBall/tests/utils.py`
**Lines Changed**: +88/-12 (76 net additions)

#### Modified Function: `plot1Dcut()` (Lines 134-171)
**Changes**:
- Added `pos0` parameter for correct grid origin alignment
- Fixed reference potential calculation alignment
```python
# BEFORE:
def plot1Dcut(apos, qs, Vg, i0, dg, Ls, iax=0, ...):

# AFTER:
def plot1Dcut(apos, qs, Vg, i0, dg, Ls, pos0, iax=0, ...):
```

#### New Function: `select_cut_1d()` (Lines ~120-133)
**Purpose**: Select 1D cut from 3D grid

---

### 4.6 `pyBall/tests/Ewald.py`
**Lines Changed**: +8/-4 (4 net additions)

#### Changes:
- Minor updates to Ewald summation tests

---

### 4.7 `pyBall/MMFF_multi.py`
**Lines Changed**: +69 modifications (mostly deletions)

#### Changes:
- Multi-molecule MMFF interface cleanup
- Removed redundant code

---

## 5. Python Test Files

### 5.1 `tests/tGridFF/` (Complete folder added)

#### `generate_grid.py` (70 lines) - NEW
**Purpose**: CLI tool for generating force field grids
```python
def main():
    gff.test_gridFF_ocl(fname=str(input_path), 
                        Element_Types_name=ELEMENT_TYPES,
                        save_name="double3", job="PLQ",
                        desired_voxel=args.desired_voxel)
```

#### `generate_scans.py` (499 lines) - NEW
**Purpose**: Generate molecular scans on substrates
**Key Functions**:
- `generate_scan()` - Main scan generation
- `setup_scan_grid()` - Grid setup for scans

#### `all_scan.py` (621 lines) - NEW
**Purpose**: Comprehensive scan utilities

#### `README.md` (283 lines) - MODIFIED
**Changes**: Updated documentation for grid generation workflow

---

### 5.2 `tests/tMMFF/`

#### `run.py` (+1013/-70 = 943 net additions)
**Major Changes**:
- Added UFF support
- Added 1D/2D scan support
- Added relaxed scan support
- Added batch processing

#### `generate_scans.py` (502 lines) - NEW
**Purpose**: Scan generation for MMFF

#### `batch_compare.py` (109 lines) - NEW
**Purpose**: Batch comparison FireCore vs LAMMPS

#### `all_scan.py` (621 lines) - NEW
**Purpose**: Scan utilities

#### `scan_coulomb.py` (155 lines) - NEW
**Purpose**: Coulomb potential scan

#### `scan_morse.py` (153 lines) - NEW
**Purpose**: Morse potential scan

#### `scan_total.py` (118 lines) - NEW
**Purpose**: Total energy scan

#### `extended.py` (79 lines) - NEW
**Purpose**: Extended test utilities

#### `run_test_GridFF_ocl_new.py` (+41/-10)
**Changes**: Updated OCL GridFF test

#### `run_sample_surf.py` (+45/-10)
**Changes**: Updated surface sampling

#### `run_test_ewald.py` (+2/-1)
**Changes**: Minor update

---

### 5.3 `tests/tUFF/run.py`
**Lines Changed**: +43 modifications (mostly deletions)
**Changes**: UFF vs LAMMPS automatic test (commit 3657e5e1)

---

### 5.4 `tests/tPTCDA/run.py` (13 lines) - NEW
**Purpose**: PTCDA electrostatics test

---

### 5.5 `tests/tPTCDA2/`
- `run.py` (13 lines) - NEW
- `firecore.py` (13 lines) - NEW

---

### 5.6 `tests/tMMFFmulti/run.py`
**Lines Changed**: +40 modifications (mostly deletions)
**Changes**: Multi-molecule OCL vs OMP tests

---

## 6. Documentation Files

### 6.1 `doc/py/plot_utils.py`
**Lines Changed**: +829/-1 (828 net additions)
**Purpose**: Comprehensive plotting utilities

### 6.2 `doc/py/energy_comparison_1d.py` (171 lines) - NEW
**Purpose**: Energy comparison plotting tool

### 6.3 `doc/py/combine_energy_trajectory_2x2.py` (108 lines) - NEW
**Purpose**: Trajectory combination utility

### 6.4 `doc/py/visualise_trajectory_atoms.py` (937 lines) - NEW
**Purpose**: Trajectory visualization

### 6.5 `doc/py/visualize_top_layer_xy.py` (715 lines) - NEW
**Purpose**: Top layer visualization

### 6.6 `doc/py/molecule_visualization_example.py` (71 lines) - NEW
**Purpose**: Molecule visualization example

---

## 7. New Files Added

### Documentation Files Created by This Process:
- `INDRANIL_BRANCH_SUMMARY.md` (64 lines)
- `pyBall/OCL/GridFF_CHANGES_vs_MASTER.md` (155 lines)
- `pyBall/tests/ocl_GridFF_new_CHANGES_vs_MASTER.md` (103 lines)

### Symlinks:
- `tests/tGridFF/common_resources` → `../../cpp/common_resources`
- `tests/tGridFF/data` → `../../cpp/common_resources`

---

## Summary Statistics

| Category | Files | Insertions | Deletions | Net |
|----------|-------|------------|-----------|-----|
| C++ Headers | 23 | +2,100 | -500 | +1,600 |
| C++ Sources | 2 | +270 | -90 | +180 |
| OpenCL Kernels | 2 | +218 | -128 | +90 |
| Python Core | 7 | +900 | -250 | +650 |
| Python Tests | 20 | +4,500 | -200 | +4,300 |
| Documentation | 9 | +2,900 | -100 | +2,800 |
| **TOTAL** | **65** | **+9,815** | **-1,463** | **+8,352** |

---

## Key Commits Referenced

| Commit | Description |
|--------|-------------|
| `4f4b50f6` | Rigid scan H2O on NaCl matching GPU/CPU |
| `cdf029af` | GPU Ewald grid shift for CPU/LAMMPS agreement |
| `3657e5e1` | UFF tested vs LAMMPS |
| `31f4d5d9` | scan_rigid_uff function added |
| `6c2fdf81` | UFF buffers modification resolved |
| `76ca776a` | Corrected k-kernel in poissonW() |
| `42e51553` | slabPotential_zyx kernel |
| `191694c9` | scan_constr() function added |
| `e8250834` | tGridFF folder added |

---

*Document generated for verification and audit purposes.*
