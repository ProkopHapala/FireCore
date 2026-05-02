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

## Kpoint Parallelization Attempt (2026-05-02)

### Attempted Implementation

**Goal:** Add OpenMP parallelization over kpoints in `fortran/MAIN/libFireCore.f90` to utilize multiple CPUs during SCF iterations.

**Changes Made:**

1. **`fortran/MAIN/libFireCore.f90`** - Added OpenMP directive to kpoint loop:
   ```fortran
   !$omp parallel do private(k_temp)
   do ikpoint = 1, nkpoints
       k_temp(:) = special_k(:,ikpoint)
       call solveH ( ikpoint, k_temp )
   end do
   !$omp end parallel do
   ```

2. **`fortran/MODULES/workmat.f90`** - Made work arrays thread-private:
   ```fortran
   !$omp threadprivate(work, rwork, iwork, slam, xxxx, zzzz)
   ```

3. **`fortran/MAIN/solveH.f90`** - Added critical section for `sm12_save` allocation:
   ```fortran
   !$omp critical
   if (.not. allocated(sm12_save)) then
       allocate (sm12_save(norbitals,norbitals,nkpoints))
   end if
   !$omp end critical
   ```

4. **Fixed OpenMP compilation errors** in multiple files for gfortran compatibility:
   - `Dassemble_ca_olsxc_2c.f90`: Fixed continuation line and variable names
   - `buildh.f90`: Commented out `!$omp atomic` on array sections
   - `assemble_olsxc_off.f90`: Removed empty private clause, commented out atomic
   - `find_neigh_max.f90`: Commented out `!$ volatile` directive
   - `neighbors.f90`: Added `reduction(max:rc_max)` clause
   - `find_neighPP_max.f90`: Commented out `!$ volatile` directive
   - `project_orb.f90`: Removed undeclared variables from private/shared clauses

5. **Commented out conflicting existing OpenMP directives** to isolate kpoint parallelization:
   - `Dassemble_ca_olsxc_2c.f90`: Commented out parallel do
   - `assemble_olsxc_off.f90`: Commented out parallel do
   - `neighbors.f90`: Commented out parallel do
   - `project_orb.f90`: Commented out parallel region

### Results

**Problem:** Kpoint parallelization produced incorrect results:
- Without OpenMP: E_min = -1463.717945 eV at Lx=2.4 Å
- With OpenMP: E_min = -2048.542136 eV at Lx=2.4 Å (difference of ~585 eV)

**Root Cause:** Thread-safety issues with shared module-level arrays. The kpoint loop shares state (density matrix, charges, etc.) that must be accumulated after all kpoints are processed. Many module-level arrays in the `charges` module (Qin, Qout, Qneutral, Qinmixer, Qoutmixer, dq, QLowdin_TOT, QMulliken_TOT) are shared across threads and modified during SCF iterations, causing race conditions.

**Fundamental Issue:** The SCF loop is outside the kpoint loop in `firecore_SCF`:
```fortran
do Kscf = 1, max_scf_iterations
    call assemble_mcweda ()
    do ikpoint = 1, nkpoints  ! <-- kpoint loop
        call solveH ( ikpoint, k_temp )
    end do
    call denmat ()  ! <-- accumulates results from all kpoints
    call mixer ()
end do
```

Kpoints share state that requires careful synchronization. Making all shared state thread-safe would require extensive modifications throughout the codebase.

### Strategy for Future Implementation

**Option 1: Deep refactoring (cleanest but extensive)**
1. Create a `context` structure containing all thread-local data
2. Pass context explicitly to `solveH` and all its callees
3. Refactor functions to use context instead of global module state
4. Parallelize kpoint loop with each thread using its own context

**Option 2: Thread-local module copies (simpler but fragile)**
1. Preallocate N copies of all shared module arrays (N = number of threads)
2. Use OpenMP `threadprivate` on ALL module arrays (not just workmat)
3. This is fragile because it requires identifying every shared variable across the entire codebase

**Option 3: Hybrid approach (practical compromise)**
1. Preallocate work arrays per kpoint (as suggested by user)
2. Keep the density/charge accumulation serial (only parallelize the H(k) construction and diagonalization)
3. After parallel kpoint loop, accumulate results serially via `denmat`

**Key Insight:** `solveH` does two things:
- Builds H(k) and diagonalizes (thread-safe if work arrays are local)
- Potentially modifies global state (not thread-safe)

If we can separate these, we could parallelize just the diagonalization part. However, the current code structure tightly couples them.

### Variables That Need Thread-Local Copies

**Work arrays (already identified):**
- `workmat` module: work, rwork, iwork, slam, xxxx, zzzz
- `sm12_save` (indexed by ikpoint)

**Shared module arrays modified during SCF:**
- `charges` module: Qin, Qout, Qneutral, Qinmixer, Qoutmixer, dq, QLowdin_TOT, QMulliken_TOT
- `density` module: rho (and related arrays)
- `forces` module: if force calculation is enabled
- Other modules accessed by `solveH` callees

### Current Status

The kpoint parallelization attempt has been reverted. All OpenMP-related changes have been removed to restore the code to its working state. The calculation now produces correct results without OpenMP kpoint parallelization.

**Note:** The lattice scan is already parallelized at the Python level (over Lx points) using multiprocessing, but this is not currently being used (the `use_continuous_path` parameter controls it, and it defaults to False). This is the appropriate parallelization strategy for the current use case.

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

## Two-Ribbon H-Transfer NEB Calculation (2026-05-02)

### Overview

New functionality for studying hydrogen transfer between two graphene ribbons with different passivations (N and NH) using Nudged Elastic Band (NEB) method with DFTB+ as the solver via ASE.

### Scripts

#### 1. `tests/pyFireball/build_two_ribbons.py`

**Purpose:** Build a unit cell containing two graphene ribbons (N-passivated and NH-passivated) stacked along the y-direction with controlled hydrogen bond spacing.

**Key Functions:**

**`build_two_ribbon_cell(width_chains, length_cells, Lx, a_CC, L_Hb, shift_x)`**
- Builds N-passivated ribbon (bottom) and NH-passivated ribbon (top)
- Positions ribbons with specified hydrogen bond length `L_Hb`
- Calculates cell size: `Ly = y_span_N + y_span_NH + 2 * L_Hb`
- Supports optional x-shift for one ribbon (in factors of Lx) for registry control
- Returns atomic positions, types, elements, and lattice vectors

**`plot_geometry(apos, atypes, elems, lvs, title)`**
- Plots geometry using `plotUtils.plotGeometry`
- Shows ribbon arrangement and hydrogen bond spacing

**Parameters:**
- `width_chains`: Ribbon width in chains (default: 4)
- `length_cells`: Length in unit cells (default: 1)
- `Lx`: Lattice constant along x (default: 2.4 Å)
- `a_CC`: C-C bond length (default: 1.42 Å)
- `L_Hb`: Hydrogen bond length between ribbons (default: 2.0 Å)
- `shift_x`: x-shift for second ribbon in factors of Lx (default: 0.0)

**Output Files:**
- `two_ribbons_geometry.png` - Geometry visualization
- `two_ribbons.xyz` - XYZ format geometry
- `two_ribbons.gen` - GenFormat for DFTB+

**Usage:**
```bash
cd /home/prokophapala/git/FireCore/tests/pyFireball
python build_two_ribbons.py
```

#### 2. `tests/pyFireball/neb_h_transfer.py`

**Purpose:** Perform NEB calculation for hydrogen transfer between two ribbons using ASE's DFTB+ calculator.

**Key Functions:**

**`read_genformat(fname)`**
- Reads GenFormat geometry file
- Returns ASE Atoms object with lattice vectors

**`create_initial_final_states(atoms)`**
- Identifies H and N atoms involved in transfer
- Creates initial state (H on NH ribbon) and final state (H moved to N ribbon)

**`main()`**
- Reads GenFormat geometry from `two_ribbons.gen`
- Sets up ASE DFTB+ calculator with Slater-Koster parameters
- Creates 7-image NEB path between initial and final states
- Runs CI-NEB optimization
- Saves NEB trajectory and energy profile plot

**Parameters:**
- DFTB+ parameters: Slater-Koster `3ob-3-1`, nk_x=16
- NEB parameters: 7 images, spring force 5.0 eV/Å
- Convergence: fmax=0.05 eV/Å

**Output Files:**
- `neb_trajectory.traj` - ASE trajectory file
- `neb_energy_profile.png` - Energy vs reaction coordinate plot

**Usage:**
```bash
cd /home/prokophapala/git/FireCore/tests/pyFireball
python neb_h_transfer.py
```

#### 3. `tests/pyFireball/scan_LHb.py`

**Purpose:** Scan energy vs hydrogen bond length L_Hb for two-ribbon system, performing both rigid (single-point) and relaxed (geometry optimization) calculations.

**Key Functions:**

**`build_two_ribbon_cell(width_chains, length_cells, Lx, a_CC, L_Hb, shift_x)`**
- Same as in `build_two_ribbons.py` (duplicated for self-contained script)

**`run_dftb_calculation(apos, atypes, lvs, basis_path, do_relax, nk_x, Temperature, MixingParameter, MaxScc, SCCTolerance)`**
- Runs DFTB+ calculation with improved convergence parameters
- Supports both rigid (single-point) and relaxed (geometry optimization) modes
- Uses higher Fermi smearing (600 K) and lower mixing (0.1) for stability
- Returns energy, final geometry, and element names

**`save_xyz_movie(results, fname, label, lvs)`**
- Saves XYZ movie with properties in comment lines
- Comment format: `L_Hb=... E=... eV Lx=... Ly=... Lz=...`
- Includes lattice vectors for full reconstruction

**`main()`**
- Accepts command-line arguments for width, Lx, shift_x
- Scans L_Hb from 1.5 to 2.5 Å with step 0.05 (21 points)
- Runs rigid scan (single-point calculations)
- Runs relaxed scan (geometry optimization)
- Plots energy vs L_Hb for both cases
- Saves XYZ movies with lattice vectors in comments

**Command-line Arguments:**
- `--width`: Ribbon width in chains (default: 4)
- `--Lx`: Lattice constant Lx (default: 2.4)
- `--shift_x`: x-shift for second ribbon in factors of Lx (default: 0.0)

**Convergence Parameters:**
- Temperature: 600 K (higher Fermi smearing)
- MixingParameter: 0.1 (lower mixing for stability)
- MaxScc: 500 (more iterations)
- SCCTolerance: 1e-4 (looser tolerance)

**Output Files:**
- `LHb_scan_w{width}_Lx{Lx}[sx{shift}]_rigid.xyz` - Rigid geometries
- `LHb_scan_w{width}_Lx{Lx}[sx{shift}]_relax.xyz` - Relaxed geometries
- `LHb_scan_w{width}_Lx{Lx}.png` - Energy vs L_Hb plot
- `LHb_scan_w{width}_Lx{Lx}.log` - Scan data log

**Usage Examples:**
```bash
# w4 ribbon, no shift
python scan_LHb.py --width 4

# w6 ribbon with half-cell shift
python scan_LHb.py --width 6 --shift_x 0.5

# w8 ribbon
python scan_LHb.py --width 8
```

**Results:**
- **w4**: Minimum at L_Hb=1.85 Å (rigid), 1.80 Å (relaxed)
- **w6 (shift_x=0.5)**: Minimum at L_Hb=1.85 Å (both rigid and relaxed)
- **w8**: Minimum at L_Hb=1.85 Å (rigid), 1.70 Å (relaxed)

### Modified Modules

#### `pyBall/dftb_utils.py`

**Change:** Added `MixingParameter` parameter to `makeDFTBjob_pbc` function.

**Purpose:** Allow control of Broyden mixing parameter for improved DFTB+ convergence.

**New Parameter:**
- `MixingParameter`: Broyden mixing parameter (default: 0.2, lower values more stable)

**HSD Output:**
```
Mixer = Broyden {
  MixingParameter = 0.2
}
```

### DFTB+ Documentation

#### `tests/dftb/DFTB_docs.md`

**Purpose:** Documentation for DFTB+ integration in FireCore, including Slater-Koster parameters and input file formats.

**Contents:**
- DFTB+ installation and setup
- Slater-Koster parameter sets (3ob-3-1, mio-1-1, etc.)
- GenFormat geometry input
- HSD input file structure
- k-point sampling
- Convergence parameters
- Common issues and solutions

### Key Features

1. **Controlled Ribbon Spacing:** L_Hb parameter controls hydrogen bond distance between ribbons
2. **Registry Control:** shift_x parameter allows x-direction offset for different ribbon registries
3. **Convergence Stability:** Improved DFTB+ parameters (higher temperature, lower mixing) for difficult geometries
4. **XYZ with Lattice Vectors:** XYZ files include Lx, Ly, Lz in comments for full reconstruction
5. **Flexible Widths:** Supports w4, w6, w8 ribbons (and others via command-line)
6. **Both Rigid and Relaxed:** Provides energy curves for both fixed and optimized geometries

### References

- `tests/pyFireball/build_two_ribbons.py` - Two-ribbon geometry builder
- `tests/pyFireball/neb_h_transfer.py` - NEB calculation script
- `tests/pyFireball/scan_LHb.py` - L_Hb scan script
- `tests/pyFireball/scan_constrained.py` - Constrained H-transfer scan (new)
- `pyBall/dftb_utils.py` - DFTB+ utilities (modified)
- `tests/dftb/DFTB_docs.md` - DFTB+ documentation
- `pyBall/plotUtils.py` - Plotting utilities

## Constrained H-Transfer Scan (2026-05-02)

### Overview

After attempting to use NEB for hydrogen transfer between two graphene ribbons (which produced non-smooth paths and excessive temporary files), we implemented a manual constrained scan approach. This provides better control over the reaction coordinate and cleaner output.

### Scripts

#### 1. `tests/pyFireball/neb_h_transfer.py` (NEB Approach - Not Recommended)

**Purpose:** Attempted to use ASE's NEB implementation with DFTB+ as the solver for hydrogen transfer.

**Issues:**
- **Non-smooth path**: The NEB path was "horrible" and not smooth at all due to linear interpolation between initial and final states
- **Excessive temporary files**: Generated numerous temporary directories (`temp_neb_N_stepM`) cluttering the workspace
- **ASE DFTB calculator issues**: Generated malformed HSD input files for periodic systems
- **Complex debugging**: Multiple ASE-specific errors (module imports, calculator sharing, etc.) required extensive debugging

**Status:** The script is kept unchanged but is not recommended for use. It serves as a cautionary example of why a manual constrained scan is preferable for this use case.

#### 2. `tests/pyFireball/scan_constrained.py` (Constrained Scan - Recommended)

**Purpose:** Manual relaxed scan for hydrogen transfer with explicit atom constraints and enhanced visualization.

**Key Features:**
- Fixes specific atoms (N atoms + scanning H) during relaxation
- Scans H position along the N...N axis as the reaction coordinate
- Supports both rigid (single-point) and relaxed (geometry optimization) scans
- Monitors residual forces on fixed atoms
- Generates XYZ movie, energy profile, and per-step geometry plots with force arrows
- Visual indicators for donor N, acceptor N, and scanning H

**Key Functions:**

**`identify_scan_atoms(apos, enames, verbose=True)`**
- Automatically finds the two closest N atoms to a given H atom
- **Donor N**: Closest N to H (covalent N-H bond, ~1.0 Å)
- **Acceptor N**: Second closest N to H (hydrogen bond H...N, ~2.0 Å)
- This robust method works for any H-bond configuration, not just the specific test system
- Returns: `h_scan_idx, n_donor_idx, n_acceptor_idx, fixed_idx`

**`set_h_position(apos, enames, h_scan_idx, n_donor_idx, n_acceptor_idx, L_H)`**
- Moves the scanning H atom to a specified distance `L_H` from the donor N along the N...N axis
- Used to define the reaction coordinate for the scan

**`run_dftb(apos, enames, lvs, temp_dir, do_relax, fixed_atoms, ...)`**
- Runs DFTB+ calculation with atom constraints
- Supports custom optimizer parameters for relaxed scans
- Parses energy, relaxed geometry (from GenFormat), and forces

**`plot_step(apos, atypes, lvs, forces, fixed_idx, h_scan_idx, n_donor_idx, n_acceptor_idx, ...)`**
- Plots geometry with:
  - Bond lengths (red text labels)
  - Fixed atom indicators (red circles on N atoms)
  - Force arrows (red, scaled by magnitude)
  - Highlighted donor N (green), acceptor N (orange), scanning H (magenta)
  - Scan path line (green dashed, showing N-N axis)

**`run_scan(width_chains, Lx, scan_range, n_steps, do_relax, ...)`**
- Main scan loop
- For relaxed scans: uses previous step's relaxed geometry as starting point
- Accumulates results and generates output files

### Critical Implementation Details

#### Atom Selection (The "Which Atoms?" Problem)

**Initial Mistakes:**
1. Selected wrong H atom (H10 on NH ribbon instead of H9 at synaptic junction)
2. Selected wrong N atoms (both on same ribbon, causing H to pass through benzene ring)
3. Used heuristic based on y-coordinates which failed for different geometries

**Correct Approach:**
```python
# Find the two closest N atoms to the H
dists_n_to_h = [(n, np.linalg.norm(apos[h_scan_idx] - apos[n])) for n in n_idx]
dists_n_to_h.sort(key=lambda x: x[1])

n_donor_idx = dists_n_to_h[0][0]      # Closest N (~1.0 Å) - covalent bond
n_acceptor_idx = dists_n_to_h[1][0]   # Second closest N (~2.0 Å) - H-bond
```

This distance-based approach is robust and works for any H-bond configuration.

#### Direction Selection (The "Which Way?" Problem)

**Initial Mistakes:**
1. Selected N atoms on the same ribbon (e.g., N4 and N5 on N-ribbon)
2. Caused H to scan through the benzene ring instead of between ribbons
3. Used y-coordinate heuristics that failed

**Correct Approach:**
The distance-based selection naturally identifies the correct donor-acceptor pair because:
- Covalent N-H bond is always shortest (~1.0 Å)
- Hydrogen bond is second shortest (~2.0 Å)
- These are typically on different ribbons for synaptic junctions

**For the test system:**
- H-bond: N5-H9...N4 (N5 at y=0.92, H9 at y=-0.09, N4 at y=-2.09)
- Scan moves H9 from N5 toward N4 (proton transfer)
- N-N distance: 3.01 Å

#### DFTB+ Atom Fixing (The "How to Fix?" Problem)

**Initial Implementation:**
Used `MovedAtoms` in the `GeometryOptimization` driver to exclude fixed atoms:
```hsd
Driver = GeometryOptimization {
    MovedAtoms = {list of movable atoms}  # All atoms except fixed ones
    ...
}
```

**Critical Issue:**
Initially only fixed N atoms, but the scanning H atom was free to move during relaxation. This caused H to relax back to its equilibrium position regardless of where it was initially placed.

**Solution:**
Fix both N atoms AND the scanning H atom:
```python
# Atoms to fix: all N atoms + the scanning H atom
fixed_idx = n_idx + [h_scan_idx]
```

This ensures:
- N atoms remain fixed (structural constraint)
- Scanning H stays at the scanned position (reaction coordinate constraint)
- Only C atoms and the other H (H10) can relax

#### Relaxation Stability (The "Blowing Up" Problem)

**Initial Issues:**
1. LBFGS optimizer was unstable for constrained relaxation
2. Geometries exploded to positions like `C 43437137920.000000 -35525370790.000000`
3. Energies became unphysical (e.g., E=+96.9 eV instead of -467 eV)

**Solutions:**

**1. Use FIRE optimizer instead of LBFGS:**
```python
relax_params = dftbu.default_params.copy()
relax_params["Optimizer"] = "FIRE{StepSize = 0.1}"  # More stable for constraints
relax_params["MaxSteps"] = 500
relax_params["GradElem"] = 1e-3
```

**2. Tighten SCC convergence:**
```python
Temperature=300, MixingParameter=0.02, MaxScc=500, SCCTolerance=1e-6
```

**3. SteepestDescent also tested but failed** (calculation didn't converge)

**Note:** FIRE is generally more stable than LBFGS for constrained systems because it's designed for molecular dynamics-style optimization.

#### Modified `dftb_utils.py`

**Changes:**

1. **Added `fixed_atoms` parameter to `makeDFTBjob_pbc`:**
```python
def makeDFTBjob_pbc(..., fixed_atoms=None):
    ...
    if opt and fixed_atoms:
        fixed_1based = sorted([i+1 for i in fixed_atoms])
        all_idx = set(range(1, natoms+1))
        moved_idx = sorted(all_idx - set(fixed_1based))
        moved_str = ' '.join(str(i) for i in moved_idx) if moved_idx else "1:-1"
```

2. **Added `Analysis { CalculateForces = Yes }`** to enable force parsing for monitoring constraints

3. **Added `params` parameter to `run_dftb`** to allow custom optimizer settings for relaxed scans

### Usage Examples

#### Rigid Scan (Single-Point Calculations)

```bash
cd /home/prokophapala/git/FireCore/tests/pyFireball
python scan_constrained.py --width 4 --n_steps 13 --L_min 0.8 --L_max 2.5 --rigid --outdir scan_constrained_w4_rigid
```

#### Relaxed Scan (Geometry Optimization)

```bash
python scan_constrained.py --width 4 --n_steps 13 --L_min 0.8 --L_max 2.5 --outdir scan_constrained_w4_relax
```

**Note:** Relaxed scan uses FIRE optimizer with tight SCC convergence for stability.

### Output Files

**Per-Step Geometry Plots:**
- `geom_{tag}_step{N}.png` - Geometry with bonds, force arrows, fixed atom indicators
- Shows donor N (green), acceptor N (orange), scanning H (magenta)
- Green dashed line shows scan path (N-N axis)

**Summary Files:**
- `scan_{tag}.xyz` - XYZ movie with properties in comments
- `scan_{tag}.png` - Energy profile + residual force on fixed atoms

**XYZ Comment Format:**
```
L_H=1.5000 E=-466.143700 eV F_fixed=277.682 Lx=2.40 Ly=11.70 Lz=20.00
```

### Results for w4 Ribbon

**Rigid Scan:**
- Energy minimum at L_H ≈ 1.1-1.2 Å (equilibrium N-H bond)
- Energy rises monotonically as H moves toward acceptor
- Barrier: ~0.8 eV
- Residual forces on fixed atoms: ~411 eV/Å (expected for rigid positions)

**Relaxed Scan:**
- Energy minimum at L_H ≈ 1.0-1.2 Å
- Barrier: ~0.8 eV (similar to rigid, but with C relaxation)
- Residual forces on fixed atoms: 278-476 eV/Å (increases as H moves from equilibrium)
- C atoms relax to accommodate H position

### Key Lessons

1. **Atom selection must be distance-based, not heuristic**: Using y-coordinates or other heuristics fails for different geometries. The "two closest N atoms to H" rule is robust.

2. **Fix the scanning atom during relaxation**: If you're scanning a reaction coordinate by moving an atom, that atom must also be fixed during relaxation, otherwise it will relax back to equilibrium.

3. **Use FIRE optimizer for constrained systems**: LBFGS can be unstable for constrained relaxation. FIRE is more stable for molecular systems.

4. **Tighten SCC convergence for difficult geometries**: Lower temperature (300 K), lower mixing (0.02), and tighter tolerance (1e-6) help stability.

5. **Visual debugging is essential**: The geometry plots with force arrows and atom labels were crucial for identifying the wrong atom selections and directions.

6. **Manual scan > NEB for simple reactions**: For single-bond transfers with a clear reaction coordinate, a manual constrained scan gives better control and cleaner output than NEB.

### References

- `tests/pyFireball/scan_constrained.py` - Constrained scan implementation
- `tests/pyFireball/neb_h_transfer.py` - NEB attempt (not recommended)
- `pyBall/dftb_utils.py` - DFTB+ utilities with atom fixing support
- `pyBall/plotUtils.py` - Plotting utilities (extended for force arrows)
