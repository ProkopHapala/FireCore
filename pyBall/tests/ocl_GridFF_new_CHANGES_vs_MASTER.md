# ocl_GridFF_new.py Changes vs Master Branch

**Branch:** Indranil  
**File:** `pyBall/tests/ocl_GridFF_new.py`  
**Date:** 2025-12-10  
**Master lines:** 701 | **Indranil lines:** 884

---

## Summary of Changes

This document describes all changes made to `ocl_GridFF_new.py` in the Indranil branch compared to the master branch.

---

## Key New Features

### 1. Improved z0 Calculation (Lines 251-260)
**Purpose:** Consistent grid origin calculation based on topmost atom layer.

**Commit context:**
- "Fixed grid origin calculation in ocl_GridFF_new.py"
- "Change algorithm for searching z0 to find the first full layer"

```python
# Group atoms by their z-coordinate and find the most populated z-plane
z_coords = xyzq[:, 2]
unique_zs, counts = np.unique(z_coords, return_counts=True)
max_count = counts.max()
most_frequent_zs = unique_zs[counts == max_count]
z0 = most_frequent_zs.max()
```

### 2. Desired Voxel Parameter (Line 239-240)
**Purpose:** Allow specifying desired voxel size for grid generation.

```python
def test_gridFF_ocl(..., desired_voxel=0.1):
    grid = GridShape(desired_voxel=desired_voxel, lvec=atoms.lvec)
```

### 3. Buffer Verification Function (Lines 327-347)
**Purpose:** Quick check of V_Coul_buff contents for debugging.

```python
def check_vcoul_buffer(clgff):
    # Verifies buffer shape and range
```

### 4. Import Changes (Lines 1-16)
- Added `import pyopencl as cl`
- Changed `from ..AtomicSystem import AtomicSystem` to `au.AtomicSystem`

---

## Modified Functions

### 1. `make_atoms_arrays()` (Lines 197-222)
**Changes:**
- Removed verbose per-atom debug prints
- Kept essential `Qtot` print
- Uses `au.AtomicSystem` instead of direct import

### 2. `test_gridFF_ocl()` (Lines 239-882)
**Changes:**
- Added `desired_voxel` parameter
- Improved z0 calculation algorithm
- Added buffer verification in PLQ job
- Cleaned up excessive debug prints
- Added grid parameter saving

---

## Key Commit Messages from Debug Branch

1. **b68be3ac** - "Change algorithm for searching z0 to find the first full layer"
2. **caf4870b** - "Fixed grid origin calculation in ocl_GridFF_new.py"
3. **4f4b50f6** - "rigid scan of H2O molecule on NaCl_1x1_L1 is matching between GPU and CPU"
4. **42e51553** - "debuged PLQ_lin in ocl_GridFF_new.py"
5. **88cfef72** - "implemented calculation of full GridFF PLQ on GPU"

---

## Debug Prints Retained (Essential)

1. **Line 241:** `test_gridFF_ocl() START`
2. **Line 261:** `z0` value
3. **Line 264:** Grid parameters (ns, dg)
4. **Line 350:** Coulomb calculation start
5. **Line 378:** Morse calculation start
6. **Line 380:** autoPBC result
7. **Lines 472-475:** Final min/max values for Paul, Lond, Coul

---

## Verification Checklist

- [x] Syntax valid
- [x] z0 calculation algorithm preserved
- [x] desired_voxel parameter added
- [x] Buffer verification function preserved
- [x] Excessive debug prints removed
- [x] Essential prints retained
