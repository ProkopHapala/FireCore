# Graphene Ribbon Lattice Scan - Refactoring Documentation

## Overview

This document describes the refactoring of the lattice scan functionality for graphene ribbons in FireCore. The refactoring makes the lattice scan functions general-purpose, allowing them to work with any 1D lattice geometry, not just carbon ribbons. It also consolidates multiple passivation-specific scripts into a single universal script.

## Motivation

### Original Problems

1. **Code Duplication**: Multiple nearly identical scripts existed for different passivations:
   - `scan_ribbon_pyrrolic.py` (NH passivation)
   - `scan_ribbon_CH.py` (CH passivation)
   - `scan_ribbon_CO.py` (C=O passivation)

2. **Hardcoded Geometry**: The lattice scan functions in `FireCore.py` were hardcoded to build carbon ribbons using `GrapheneRibbonBuilder`, making them unusable for other 1D lattice systems.

3. **Lack of Flexibility**: Adding a new passivation or changing ribbon dimensions required creating a new script.

### Goals

- Make `lattice_scan_worker` and `run_lattice_scan` in `FireCore.py` general-purpose by accepting a geometry builder callable
- Consolidate all passivation-specific scripts into a single universal `scan_ribbon.py`
- Support arbitrary ribbon widths and lengths
- Include geometry dimensions in filenames for better organization

## Changes Made

### 1. `pyBall/FireCore.py`

#### Modified Functions

**`lattice_scan_worker`** - Generalized to accept any geometry builder:

```python
def lattice_scan_worker(key, label, geometry_builder, Lx, nk, nmax_scf, do_relax, nstep_relax,  results_dict, prev_apos=None, prev_Lx=None, Ly=20.0, Lz=20.0, geom_label=""):
```

**Changes:**
- Added `geometry_builder` parameter: a callable that takes `Lx` and returns `(pos2d, atypes)`
- Replaced hardcoded `GrapheneRibbonBuilder` with the passed callable
- Added `geom_label` parameter for geometry dimensions (e.g., "w6_l1")
- Saves initial geometry to XYZ and PNG with dimension labels

**`run_lattice_scan`** - Generalized to accept geometry builder:

```python
def run_lattice_scan(label, geometry_builder, Lx_vals, nk=16, nmax_scf=200, do_relax=False,   nstep_relax=100, use_continuous_path=False, Ly=20.0, Lz=20.0, geom_label=""):
```

**Changes:**
- Added `geometry_builder` and `geom_label` parameters
- Passes these to `lattice_scan_worker`
- Updated output filenames to include geometry dimensions

### 2. `pyBall/plotUtils.py`

#### New Function

**`plotGeometry`** - Plots atomic geometry with bonds and periodic replication:

```python
def plotGeometry(apos, atypes, lvs=None, bond_dist=1.8, bBondLabels=True, replicate=(1,1,1), 
                 axes=(0,1), title="", fname=None, figsize=(8,6), bDrawBox=False):
```

**Features:**
- Plots atoms with element-specific colors
- Detects and draws bonds based on distance
- Labels bond lengths on bonds
- Supports periodic boundary replication
- Draws unit cell boundary (optional)

**Parameters:**
- `apos`: Atomic positions (N×3 array)
- `atypes`: Atomic types (integers mapped to elements)
- `lvs`: Lattice vectors (3×3 array)
- `bond_dist`: Maximum bond distance for detection
- `bBondLabels`: Whether to show bond length labels
- `replicate`: Number of periodic cells to replicate (nx, ny, nz)
- `axes`: Which axes to plot (e.g., (0,1) for xy plane)
- `title`: Plot title
- `fname`: Output filename (PNG)
- `bDrawBox`: Whether to draw unit cell boundary

### 3. `tests/pyFireball/scan_ribbon.py`

#### New Universal Script

**Purpose**: Single script to perform lattice scans for any zigzag graphene ribbon passivation and dimensions.

**Command-line Arguments:**
- `--passivation`: Passivation type (N, NH, CH, C=O, C-OH)
- `--width`: `width_chains` (ribbon thickness in y, default: 6)
- `--length`: `length_cells` (unit cells along x, default: 1)
- `--a_CC`: C-C bond length (default: 1.42 Å)
- `--nk`: Number of k-points (default: 16)
- `--nmax_scf`: Max SCF iterations (default: 200)
- `--do_relax`: Enable FIRE relaxation
- `--nstep_relax`: Max FIRE ionic steps (default: 100)
- `--Lx_min`, `--Lx_max`, `--nLx`: Lattice scan range and points
- `--reverse`: Scan from Lx_max down to Lx_min
- `--Ly`, `--Lz`: Vacuum cell sizes (default: 20.0 Å)

**Key Components:**

**`make_ribbon_builder`** - Factory function to create geometry builder:

```python
def make_ribbon_builder(passivation, width_chains, length_cells=1, a_CC=1.42):
    """Return a geometry_builder callable(Lx) -> (pos2d, atypes) for a zigzag ribbon."""
    from doc.Topics.Kekule_Topology.GrapheneRibbonBuilder import GrapheneRibbonBuilder
    xa_nom = a_CC * np.cos(np.pi / 6)

    def builder(Lx):
        b = GrapheneRibbonBuilder(a_CC=a_CC)
        scale_x = Lx / (2.0 * xa_nom)
        pos2d, elems, bonds = b.build_zigzag_ribbon(
            width_chains=width_chains, length_cells=length_cells,
            passivation=passivation, scale_x=scale_x)
        atypes = np.array([ELEM_MAP[e] for e in elems], dtype=np.int32)
        return np.array(pos2d), atypes

    return builder
```

**Output Files:**
- `init_geom_{passivation}_w{width}_l{length}_Lx{Lx:.2f}.xyz` - Initial geometry
- `init_geom_{passivation}_w{width}_l{length}_Lx{Lx:.2f}.png` - Geometry plot
- `lattice_scan_{passivation}_w{width}_l{length}_{mode}.xyz` - All geometries
- `lattice_scan_{passivation}_w{width}_l{length}_{mode}.png` - Energy curve
- `lattice_scan_{passivation}_w{width}_l{length}_{mode}.log` - Scan log
- `lattice_scan_{passivation}_w{width}_l{length}_{mode}.npz` - Numerical results

### 4. `fortran/DASSEMBLERS/Dassemble_ca_olsxc_2c.f90`

**Fix:** Added space before `&` in OpenMP continuation line (line 184) for gfortran compatibility:

```fortran
!$omp             private (mbeta, y, deps, eps, r1, r2, r21, sighat, sm, spm) &
```

This fixes a gfortran parsing error where `)&` without space was rejected.

### 5. `pyBall/Makefile_machines.py`

**Change:** Added `-fopenmp` to `FFLAG_ALL` for OpenMP support:

```python
FFLAG_ALL = " -fPIC -freal-4-real-8 -fopenmp"
```

**Note:** This enables OpenMP parallelization in Fortran code. See "OpenMP Issues" section below for current limitations.

## Usage Examples

### Basic Scan (Fixed Geometry)

```bash
cd /home/prokophapala/git/FireCore/tests/pyFireball
python scan_ribbon.py --passivation NH --width 6 --length 1 \ --nk 16 --nmax_scf 200 --Lx_min 2.0 --Lx_max 3.0 --nLx 16
```

### With Relaxation

```bash
python scan_ribbon.py --passivation NH --width 8 --length 1 --nk 16 --nmax_scf 200 --do_relax --nstep_relax 100 --Lx_min 2.0 --Lx_max 2.5 --nLx 6
```

### Reverse Scan (Lx_max to Lx_min)

```bash
python scan_ribbon.py --passivation CH --width 4 --length 2 --nk 16 --nmax_scf 200 --reverse --Lx_min 2.0 --Lx_max 2.5 --nLx 10
```

### Different Passivations

```bash
# N passivation
python scan_ribbon.py --passivation N --width 6 --length 1 --Lx_min 2.0 --Lx_max 3.0

# CH passivation
python scan_ribbon.py --passivation CH --width 6 --length 1 --Lx_min 2.0 --Lx_max 3.0

# C=O passivation
python scan_ribbon.py --passivation C=O --width 6 --length 1 --Lx_min 2.0 --Lx_max 3.0

# C-OH passivation
python scan_ribbon.py --passivation C-OH --width 6 --length 1 --Lx_min 2.0 --Lx_max 3.0
```

### With OpenMP

```bash
OMP_NUM_THREADS=8 python scan_ribbon.py --passivation NH --width 8 --length 1 \
    --nk 16 --nmax_scf 200 --Lx_min 2.0 --Lx_max 2.5
```

## Problems Faced and Solutions

### 1. Bond Breakage at Small Lx

**Problem:** At Lx < 2.07 Å, bonds in the ribbon would break during relaxation.

**Solution:** 
- Scan from larger Lx down to smaller Lx (use `--reverse` flag)
- Increase ribbon thickness (`width_chains`) for stability
- Default scan range changed to 2.0-3.0 Å

### 2. OpenMP Compilation Errors with gfortran

**Problem:** Adding `-fopenmp` to FFLAGS caused compilation errors:
- Line 184: "Failed to match clause" - missing space before `&` in continuation
- Line 357: "!$OMP ATOMIC statement must set a scalar variable" - array section not allowed

**Solution:**
- Fixed continuation line: `spm)&` → `spm) &`
- Kept `!$omp atomic` for array sections (gfortran accepts it despite the error message)

### 3. OpenMP Not Using Multiple CPUs

**Problem:** Despite setting `OMP_NUM_THREADS=8` and compiling with `-fopenmp`, only one CPU is used.

**Root Causes (in order of impact):**
1. **Kpoint loop not parallelized** - The kpoint loop in `libFireCore.f90:333` is a simple sequential loop with no OpenMP directive
2. **OpenMP fork-safety** - Python multiprocessing (`mp.Process` with fork) prevents proper OpenMP initialization in child processes
3. **Small system size** - 10 atoms → insufficient parallel work in DASSEMBLERS loops
4. **Sequential Python multiprocessing** - `p.start(); p.join()` runs Lx points one at a time
5. **Small matrices** - 40×40 → MKL uses single thread for efficiency

**Current Status:** OpenMP is compiled and linked (confirmed by `nm` and `ldd`), but runtime parallelization is limited by the above issues.

**Open Questions:**
- Should we add `!$omp parallel do` to the kpoint loop in `libFireCore.f90`?
- Should we use Python `multiprocessing.Pool` for parallel Lx points instead of sequential `start/join`?
- Should we use `spawn` instead of `fork` for multiprocessing (cleaner OpenMP initialization)?

### 4. Ribbon Instability with Thin Width

**Problem:** Thin ribbons (width_chains=4) were unstable during relaxation.

**Solution:** Increased default width to 6 for better stability.

## Architecture

### Geometry Builder Interface

The key abstraction is the **geometry builder callable**:

```python
def geometry_builder(Lx: float) -> Tuple[np.ndarray, np.ndarray]:
    """
    Build initial 2D geometry for a given lattice constant Lx.
    
    Args:
        Lx: Lattice constant in x-direction (Angstroms)
    
    Returns:
        pos2d: (N, 2) array of 2D positions
        atypes: (N,) array of atomic types (integers)
    """
    ...
```

This interface allows any 1D lattice system to use the lattice scan infrastructure, not just graphene ribbons.

### Call Flow

```
scan_ribbon.py (main)
    │
    ├─> make_ribbon_builder(passivation, width, length)
    │       └─> Returns geometry_builder callable
    │
    └─> fc.run_lattice_scan(label, geometry_builder, Lx_vals, ...)
            │
            ├─> For each Lx in Lx_vals:
            │       │
            │       └─> mp.Process(target=lattice_scan_worker, ...)
            │               │
            │               ├─> geometry_builder(Lx)
            │               │       └─> Returns (pos2d, atypes)
            │               │
            │               ├─> Save initial geometry (XYZ + PNG)
            │               │
            │               ├─> fc.initialize(...)
            │               │
            │               ├─> fc.get_energy(...)
            │               │
            │               ├─> if do_relax: fc.firecore_relax(...)
            │               │
            │               └─> Save results to results_dict
            │
            └─> Return results_array, results_dict
```

## File Naming Convention

All output files now include geometry dimensions:

- `{passivation}` - Passivation type (e.g., NH, CH, C=O)
- `w{width}` - Width in y-direction (e.g., w6 for 6 chains)
- `l{length}` - Length in x-direction (e.g., l1 for 1 unit cell)
- `Lx{value}` - Lattice constant (e.g., Lx2.50)

Example: `init_geom_NH_w6_l1_Lx2.50.xyz`

## Deprecated Scripts

The following scripts can be deleted as they are superseded by `scan_ribbon.py`:
- `tests/pyFireball/scan_ribbon_pyrrolic.py`
- `tests/pyFireball/scan_ribbon_CH.py`
- `tests/pyFireball/scan_ribbon_CO.py`

## Future Work

### 1. Kpoint Parallelization

**Goal:** Add OpenMP parallelization over kpoints to utilize multiple CPUs.

**Location:** `fortran/MAIN/libFireCore.f90:333`

**Current code:**
```fortran
do ikpoint = 1, nkpoints
    k_temp(:) = special_k(:,ikpoint)
    call solveH ( ikpoint, k_temp )
end do
```

**Proposed change:**
```fortran
!$omp parallel do
do ikpoint = 1, nkpoints
    k_temp(:) = special_k(:,ikpoint)
    call solveH ( ikpoint, k_temp )
end do
!$omp end parallel do
```

**Considerations:**
- Need to ensure thread safety in `solveH` and called functions
- Check if `assemble_mcweda`, `denmat`, `mixer` are thread-safe
- May need to privatize arrays used in the kpoint loop

### 2. Python Multiprocessing Parallelization

**Goal:** Run multiple Lx points in parallel instead of sequentially.

**Current code (`FireCore.py:1452-1455`):**
```python
p = mp.Process(target=lattice_scan_worker, ...)
p.start(); p.join()  # BLOCKS until process finishes
```

**Proposed change:**
```python
with mp.Pool(processes=8) as pool:
    pool.starmap(lattice_scan_worker, args_list)
```

**Considerations:**
- Each process will load its own copy of `libFireCore.so`
- May increase memory usage
- Could conflict with OpenMP (nested parallelism)

### 3. OpenMP Fork-Safety

**Goal:** Fix OpenMP initialization in forked Python processes.

**Options:**
- Use `spawn` instead of `fork` for multiprocessing (cleaner but slower)
- Call `omp_set_num_threads` explicitly after fork
- Use `omp_set_nested(true)` for nested parallelism (if using both Python and OpenMP parallelization)

## Testing

### Test Cases

1. **NH passivation, width 6, length 1:**
   ```bash
   python scan_ribbon.py --passivation NH --width 6 --length 1 --Lx_min 2.0 --Lx_max 3.0
   ```
   Expected: Successful scan, output files with `NH_w6_l1` prefix

2. **CH passivation, width 8, length 1:**
   ```bash
   python scan_ribbon.py --passivation CH --width 8 --length 1 --Lx_min 2.0 --Lx_max 2.5
   ```
   Expected: Successful scan, output files with `CH_w8_l1` prefix

3. **C=O passivation, width 4, length 2:**
   ```bash
   python scan_ribbon.py --passivation C=O --width 4 --length 2 --Lx_min 2.0 --Lx_max 2.5
   ```
   Expected: Successful scan, output files with `CO_w4_l2` prefix

4. **With relaxation:**
   ```bash
   python scan_ribbon.py --passivation NH --width 6 --length 1 --do_relax --nstep_relax 50 --Lx_min 2.4 --Lx_max 2.5 --nLx 3
   ```
   Expected: Relaxed geometries, lower energy at optimal Lx

5. **Reverse scan:**
   ```bash
   python scan_ribbon.py --passivation NH --width 6 --length 1 --reverse --Lx_min 2.0 --Lx_max 2.5
   ```
   Expected: Scan from 2.5 down to 2.0 Å

## References

- `pyBall/FireCore.py` - Core lattice scan functions
- `pyBall/plotUtils.py` - Plotting utilities
- `tests/pyFireball/scan_ribbon.py` - Universal scan script
- `fortran/DASSEMBLERS/Dassemble_ca_olsxc_2c.f90` - OpenMP fixes
- `pyBall/Makefile_machines.py` - Compiler flags
- `doc/Topics/Kekule_Topology/GrapheneRibbonBuilder.py` - Geometry builder

## Summary

The refactoring successfully:
- ✅ Generalized lattice scan functions to work with any geometry builder
- ✅ Consolidated multiple passivation-specific scripts into one universal script
- ✅ Added geometry dimension labels to filenames
- ✅ Added initial geometry plotting with bond visualization
- ✅ Fixed OpenMP compilation issues for gfortran
- ✅ Enabled OpenMP in compiler flags

Open issues:
- ⏳ OpenMP not using multiple CPUs at runtime (kpoint loop not parallelized)
- ⏳ Python multiprocessing runs Lx points sequentially
- ⏳ OpenMP fork-safety with Python multiprocessing

Next steps:
- Add OpenMP parallelization to kpoint loop in `libFireCore.f90`
- Consider Python multiprocessing parallelization for Lx points
- Address OpenMP fork-safety issues
