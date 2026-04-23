# STM Implementation Review Notes

## Overview

This document reviews the STM (Scanning Tunneling Microscopy) current computation framework that combines Fireball DFTB Hamiltonian with GPU-accelerated transport calculations.

## What We Are Trying To Do

Compute STM current between a tip and sample molecule using:

1. **Fireball Fortran DFTB** - Provides the electronic structure (Hamiltonian H, Overlap S, wavefunction coefficients C)
2. **GPU-accelerated transport** - Fast evaluation of tunneling current using perturbative approaches
3. **QM/MM framework** - Classical MD provides geometry, quantum transport computed at each snapshot

## Theoretical Framework

### Full NEGF Approach (Exact)

The total system is partitioned into: Left Lead (L), Tip Molecule (T), Surface Molecule (S), Right Lead (R)

Total Hamiltonian:
```
H(t) = H_L + H_R + H_T(R_T) + H_S(R_S) + V_LT + V_SR + V_TS(R_TS)
```

Current via Caroli formula (Landauer-Büttiker):
```
I(V) = (2e/h) ∫ T(E,V) [f(E-μ_L) - f(E-μ_R)] dE
T(E) = Tr[Γ_L G^r Γ_R G^a]
```

Where G^r = [(E+iη)S - H - Σ_L - Σ_R]^{-1}

### Fast Perturbative Approach (GPU)

For weak coupling V_TS, transmission simplifies to:
```
T(E) ≈ Tr[A_T(E) V_TS A_S(E) V_ST]
```

Where A(E) = G^r Γ G^a is the spectral function (projected DOS matrix).

**Diagonal approximation (Tersoff-Hamann limit):**
```
T(E) ≈ Σ_{μ∈T} Σ_{ν∈S} PDOS_μ^T(E) · |V_μν|² · PDOS_ν^S(E)
```

This is O(N_T × N_S) - trivial for GPU.

## Hamiltonian Mapping from Fireball (CRITICAL PITFALLS)

### Fireball Sparse Format

Fireball exports H and S in **blocked sparse format**:
```
H[i, ineigh, :nj, :ni]  # [natoms, neigh_max, numorb_max, numorb_max]
S[i, ineigh, :nj, :ni]
```

Key data structures:
- `neigh_j[i, ineigh]` - neighbor atom index (1-based Fortran indexing)
- `neigh_b[i, ineigh]` - boundary flag (0 = same cell, else periodic image)
- `neigh_self[i]` - index where jatom == iatom AND mbeta == 0
- `iatyp[i]` - atomic number Z for atom i
- `num_orb[ispec]` - orbitals per species
- `nzx` - mapping from atom Z to species index

### Critical Mapping Rules

**1. NEVER assume self-block is last neighbor slot**
```python
# WRONG:
self_block = H[i, -1, :, :]

# CORRECT:
for ineigh in range(neighn[i]):
    if neigh_j[i, ineigh] == i+1 and neigh_b[i, ineigh] == 0:
        self_block = H[i, ineigh, :, :]
```

**2. Fortran uses 1-based indexing, Python uses 0-based**
```python
j_atom_python = neigh_j[i, ineigh] - 1  # Convert to 0-based
```

**3. Block layout is transposed**
```python
# H[i, ineigh, :nj, :ni] has shape (n_orb_j, n_orb_i)
# When assembling dense matrix, need transpose:
M[i0:i0+ni, j0:j0+nj] += blk.T
```

**4. Orbital indexing varies by species**
- Carbon: 4 orbitals (s, px, py, pz)
- Hydrogen: 1 orbital (s)
- Oxygen: 4 orbitals (s, px, py, pz)
- Must map atom index → species index → orbital count

### Dense Conversion Function

The `_blocked_to_dense` function correctly handles this:

```python
def _blocked_to_dense(sparse_data, H_blocks, natoms):
    # 1. Compute orbital layout per atom
    n_orb_atom, offs = _orbital_layout(sparse_data, natoms)
    
    # 2. Find self-neighbor slots
    neigh_self = np.full(natoms, -1, dtype=np.int32)
    for i in range(natoms):
        ii = i + 1  # Fortran 1-based
        for ineigh in range(int(neighn[i])):
            if int(neigh_j[i, ineigh]) == ii and int(neigh_b[i, ineigh]) == 0:
                neigh_self[i] = ineigh
                break
    
    # 3. Assemble dense matrix
    for i in range(natoms):
        ni = int(n_orb_atom[i])
        i0 = int(offs[i])
        for ineigh in range(int(neighn[i])):
            if ineigh == int(neigh_self[i]):
                j = i
            else:
                j = int(neigh_j[i, ineigh]) - 1
            nj = int(n_orb_atom[j])
            j0 = int(offs[j])
            blk = H_blocks[i, ineigh, :nj, :ni]
            M[i0:i0+ni, j0:j0+nj] += blk.T
    return M
```

## Current Implementation in test_stm_orbital_projection.py

### Test Configuration

- **Sample molecule**: PTCDA (38 atoms)
- **Tip**: Single H atom (1 s-orbital)
- **Scan heights**: z = 2.0 Å and 10.0 Å
- **MO range**: HOMO-4 to LUMO+4 (9 orbitals)

### Three Methods Compared

1. **Fortran orb2points** - Reference orbital projection from Fireball
2. **OpenCL true basis** - GPU projection with finite cutoff from Fdata basis
3. **OpenCL exp tail** - GPU projection with exponential vacuum extension

### Response Amplitude Method (NEW)

The test also implements a response amplitude method:

```python
# Combined system: tip + sample
A = [[E+iη-E_tip,  (E+iη)S_ts - H_ts],
     [(E+iη)S_st - H_st,  (E+iη)S_s - H_s]]

# Solve A x = b where b = [1, 0, ..., 0]^T (tip excited)
# Response x_s = -G0 * a_st * x_tip

# STM amplitude (MO-projected):
t_resp(r; E) = |C_MO^T · x_s|²
```

This uses block Gaussian elimination (Dyson equation) to avoid O(N³) inversion.

**Precomputation on CPU (once per MO energy):**
```python
G0 = inv((E+iη)S_s - H_s)  # Sample Green's function
v = C_MO^T @ G0
```

**GPU kernel (per grid point):**
- Build a_st = (E+iη)S_ts - H_ts (sparse, atoms within rcut)
- Compute scalar products: v·a_st, a_st^T·G0·a_st
- Compute t_resp = |x_tip|² |v·a_st|²

## Fdata Linking Issue (RESOLVED)

**Problem**: 
- `tests/pyFireball/Fdata` symlink pointed to `/home/prokop/git/FireCore/tests/Fireball/Fdata_HCNOS` (wrong username)
- Actual path should be `/home/prokophapala/git/FireCore/tests/Fireball/Fdata_HCNOS`

**Solution**:
```bash
cd tests/pyFireball
unlink Fdata
ln -s ../Fireball/Fdata_HCNOS Fdata
```

**Chain of symlinks**:
- `tests/pyFireball/Fdata` → `../Fireball/Fdata_HCNOS`
- `tests/Fireball/Fdata_HCNOS` → `/home/prokophapala/Fireball_Data/Fdata_HCNO/`

The actual data lives in `/home/prokophapala/Fireball_Data/Fdata_HCNO/` which contains:
- Integral tables (2-center, 3-center)
- Basis functions (`.wf` files)
- Slater-Koster parameters

## Key Functions in Fireball Interface

### Hamiltonian Export
- `firecore_get_HS_dims()` - Get array dimensions
- `firecore_get_HS_neighs()` - Get neighbor lists, orbital indices, species info
- `firecore_get_HS_sparse()` - Export sparse H and S matrices

### Wavefunction Export
- `firecore_get_eigen(ikp, norb)` - Get orbital energies
- `firecore_getPointer_wfcoef()` - Get pointer to wavefunction coefficients

### Sampling Functions
- `firecore_orb2points(iband, ikpoint, npoints, points)` - Sample MO at arbitrary points
- `firecore_dens2points(npoints, points, f_den, f_den0)` - Sample electron density

## GPU Infrastructure

### Grid Projection (TESTED & WORKING)
- **Location**: `pyBall/FireballOCL/Grid.py`
- **Kernel**: `project_density_sparse_tiled` in `pyBall/FireballOCL/cl/Grid.cl`
- **Features**:
  - Tiled voxel processing (8³ blocks)
  - Catmull-Rom spline interpolation for radial functions
  - Real spherical harmonics for angular part
  - Sparse neighbor lists for efficiency

### Hamiltonian (PARTIALLY IMPLEMENTED)
- **Location**: `pyBall/FireballOCL/OCL_Hamiltonian.py`
- **Kernel**: `pyBall/FireballOCL/cl/hamiltonian.cl`
- **Features**:
  - 2-center spline interpolation
  - Slater-Koster rotation for s+p basis
  - 3-center recovery and assembly
  - Status: Infrastructure exists but not fully tested

## Next Steps

### Phase 1: Rigid Scan (Validation)
- Reproduce orbital-symmetry-resolved STM images
- Validate against full Fireball calculations
- Test with H2O (small system, 6x6 matrices) for detailed debugging

### Phase 2: Static PDOS MD
- Simulate junction at finite temperature
- Rigid-molecule approximation
- Fast dynamics on GPU

### Phase 3: Soft-Mode Expansion
- Capture molecular deformation effects
- Normal mode projection
- Taylor expansion of PDOS

## Critical Implementation Notes

1. **Always use Fireball H_μν for coupling** - NOT overlap S_μν
   - H_μν contains kinetic energy operator which gives correct decay
   - Overlap decay is too slow and misses nodal structure

2. **Use Slater-Koster tables for GPU** - Same parameters as Fireball
   - Ensures consistency between CPU (Fireball) and GPU (perturbative)
   - Accounts for s, p_x, p_y, p_z orientations correctly

3. **Include Γ broadening from leads** - Without this, DOS becomes delta spikes
   - Add -iγ to diagonal of contact atoms
   - γ typically 0.01-0.1 eV

4. **Check Hamiltonian mapping carefully** - This is the #1 source of bugs
   - Verify orbital counts per atom
   - Verify self-neighbor detection
   - Verify block transposition
   - Print small matrices (H2, CH4) to validate

## Test Results

The test `test_stm_orbital_projection.py` ran successfully, generating:
- 2×3 panel figures for each MO (HOMO-4 to LUMO+4)
- At z = 2.0 Å and z = 10.0 Å
- Comparing: Fortran reference, OpenCL true basis, OpenCL exp tail, response methods

Output directory: `tests/pyFireball/export/stm_orbital_projection/`
