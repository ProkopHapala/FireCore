# GridFF.py Changes vs Master Branch

**Branch:** Indranil  
**File:** `pyBall/OCL/GridFF.py`  
**Date:** 2025-12-10  
**Master lines:** 936 | **Indranil lines:** 1050

---

## Summary of Changes

This document describes all changes made to `GridFF.py` in the Indranil branch compared to the master branch. These changes originate from the `debug/grid_geration_and_test` branch and have been cleaned (debug prints minimized, try/except removed where possible, commented plotting code removed).

---

## New Functions Added

### 1. `fit3D_with_buffer()` (Lines 310-370)
**Purpose:** Wrapper around `fit3D` that handles buffer size mismatches.

**Commit context:** 
- "corrected serious error in OpenCL fitting (in kernel BsplineConv3D in GridFF.cl)"
- "found optimal parameters for Bspline fitting (dt=0.5,damp=0.15)"

```python
def fit3D_with_buffer(self, buffer, nPerStep=10, nmaxiter=300, dt=0.5, Ftol=1e-16, 
                      damp=0.15, bConvTrj=False, bReturn=True, bPrint=False, 
                      bTime=True, bDebug=True):
```

### 2. `release_unused_buffs()` (Lines 990-1003)
**Purpose:** Release all OpenCL buffers except specified ones to free GPU memory between operations.

**Commit context:**
- "unused buffers are released during grid generation before morse after coulomb, and can run all sized of substrate scan in run.sh"

```python
def release_unused_buffs(self, keep_names=set()):
```

### 3. `make_MorseFF()` (Lines 1005-1050)
**Purpose:** New implementation of Morse force field generation with improved parameter handling.

**Commit context:**
- "GridFF Morse is now generated consistently with CPU in GridFF.h::makeGridFF_Bspline_d() and GPU"

```python
def make_MorseFF(self, atoms, REQs, nPBC=(4, 4, 0), dg=None, ng=None, lvec=None, 
                 g0=None, GFFParams=(0.1, 1.5, 0.0, 0.0), bTime=True, bReturn=True):
```

---

## Modified Functions

### 1. `__init__()` (Lines 27-39)
**Changes:**
- Removed commented-out code blocks
- Improved kernel path resolution using `os.path.dirname(__file__)` for robustness
- Removed try/except block for kernel loading (cleaner, fails fast if kernel missing)

**Before (master):** Uses relative path `../../cpp/common_resources/cl/GridFF.cl`  
**After (Indranil):** Uses absolute path derived from script location

### 2. `try_make_buff()` (Lines 41-53)
**Changes:**
- Added logic to release outdated buffers before reallocating (prevents GPU memory leaks)
- Debug prints guarded with `# if DEBUG:` comments

**Commit context:**
- "unused buffers are released during grid generation"

### 3. `makeCoulombEwald_slab()` (Lines 886-986)
**Changes:**
- Changed `niter` default from 4 to 16
- Added verification buffer logic for debugging
- Returns `V_after_slab` instead of direct `slabPotential()` call
- Added debug prints for buffer verification

**Commit context:**
- "GridFF_cl::makeCoulombEwald_slab() now works fine on GPU"
- "finished real-space smoothing of Ewald_slab"
- "implemented calculation of full GridFF PLQ (Pauli,London,Coulomb) on GPU"

### 4. `slabPotential()` (Lines 730-800)
**Changes:**
- Added debug prints for buffer dimensions and kernel parameters
- Added verification buffer copies for debugging

**Commit context:**
- "debuged PLQ_lin in ocl_GridFF_new.py, main problem was that grid-shape was not set and grid ordering of Vcoul was [iz,iy,ix]"

### 5. `laplace_real_loop_inert()` (Lines 801-850)
**Changes:**
- Added debug prints for grid shape and buffer analysis

**Commit context:**
- "Trying to make work also GridFF_cl::laplace_real_loop_inert() on GPU"

---

## Removed from Master

### 1. `import matplotlib.pyplot as plt`
- Removed as plotting is commented out

---

## Key Commit Messages from Debug Branch

1. **49980dbf** - "unused buffers are released during grid generation before morse after coulomb"
2. **cdf029af** - "Fixed grid origin calculation in ocl_GridFF_new.py"
3. **4f4b50f6** - "rigid scan of H2O molecule on NaCl_1x1_L1 is matching between GPU and CPU"
4. **76ca776a** - "corrected the way k-kernel 1/|k^2| and normalization constant is calculated"
5. **b8acdd80** - "Checked that min,max are the same for Pauli,London for NaCl_1x1 and 8x8"
6. **42e51553** - "debuged PLQ_lin in ocl_GridFF_new.py"
7. **88cfef72** - "implemented calculation of full GridFF PLQ on GPU"
8. **d21302ce** - "finished real-space smoothing of Ewald_slab"
9. **113a7412** - "GridFF_cl::makeCoulombEwald_slab() now works fine on GPU"
10. **3f952318** - "GridFF Morse is now generated consistently with CPU and GPU"

---

## Debug Prints Retained (Minimal)

The following debug prints are retained for verification purposes:

1. **Line 950-957:** Debug before slabPotential (V_Coul_buff state)
2. **Line 977-978:** Verification buffer shape and range
3. **Line 984:** Vin min/max (only if `bCheckVin=True`)
4. **Lines 433-435:** fit3D convergence (only if `bPrint=True`)

---

## Verification Checklist

- [x] Syntax valid (python3 -m py_compile passes)
- [x] All new functions from debug branch included
- [x] Buffer release logic preserved
- [x] Verification code in makeCoulombEwald_slab preserved
- [x] try/except blocks removed or replaced with hasattr checks
- [x] Commented plotting code removed
- [x] Kernel path resolution improved for robustness

---

## Merge Recommendation

**Status:** Ready for merge to master

**Notes:**
- The Indranil branch contains all essential functionality from the debug branch
- Debug prints are minimal but sufficient for verification
- No try/except blocks that could hide errors
- Buffer management improved to prevent GPU memory leaks
