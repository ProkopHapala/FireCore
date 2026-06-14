# Plan: STM Response Function with Tip Atom

## Objective
Implement rigorous STM response function where the tip atom is explicitly included as an additional atom in the quantum system, coupled via off-diagonal hopping to the sample atoms.

## 1. Physics and Mathematics

### System Definition
Combined basis: sample (s, n_s orbitals) + tip (t, n_t orbitals)

```
H = [ H_s    H_st ]      S = [ S_s    S_st ]
    [ H_ts   H_t  ]          [ S_ts   S_t  ]
```

- H_s, S_s: from Fireball (sample only)
- H_t, S_t: tip atom (single H 1s, S_t=1, H_t=E_tip)
- H_st, S_st: tip-sample coupling (computed per tip position)

### Full Green's Function
```
G(E) = [(E + iη)S - H]^{-1} = [ a_tt    a_ts^T ]^{-1}
                                [ a_st    A_ss    ]
```

Where:
- a_tt = (E + iη)S_t - H_t
- A_ss = (E + iη)S_s - H_s
- a_st = (E + iη)S_ts - H_ts

### Block Elimination (Dyson)
Precompute once per energy:
```
G_0 = A_ss^{-1} = [(E + iη)S_s - H_s]^{-1}
v = C^T G_0   (row vector, size n_s)
```

For each tip position r:
```
a_st = (E + iη) S_ts(r) - H_ts(r)    # coupling vector (size n_s)
s1 = v · a_st^H                        # complex scalar
s2 = a_st^T · G_0 · a_st^H           # complex scalar
x_tip = 1 / (a_tt - s2)
resp(E, r) = |x_tip|^2 · |s1|^2
```

### Coupling Approximation
For vacuum tunneling, S_ts ≈ 0, so:
```
a_st ≈ -H_ts(r)
```

H_ts can be computed via:
1. **Fdata path**: Fireball radial integrals + SK angular factors (exact, slow)
2. **Exponential path**: H_ts = exp(-β(r-r0)) × SK angular factors (fast, approximate)

## 2. Relevant Code Modules and Functions

### STM_utils.py (shared utilities)
- `get_orbital_layout(data, natoms)` - orbital count per atom
- `sparse_to_dense(data, H_blocks, natoms)` - sparse → dense matrix
- `dense_to_sparse_blocks(data, M_dense, natoms, numorb_max)` - dense → sparse blocks
- `build_plane_grid(atomPos, z, size, n)` - build XY projection grid

### STM.py (core physics)
- `response_amplitude_map(...)` - CPU reference implementation
  - Computes G_0 = A_ss^{-1}
  - For each tip position: builds H_ts, S_ts via coupling_builder
  - Block elimination: resp = |x_tip|^2 · |v·a_st^H|^2
- `build_inter_system_blocks_fdata(...)` - exact Fdata coupling
- `build_inter_system_blocks_exp_sk(...)` - exponential SK coupling
- `_blocks_to_dense_vector(...)` - convert block list to dense matrix
- `_atom_orb_starts(norb_per)` - orbital offset array

### Grid.py / GridProjector (GPU kernels)
- `response_amplitude_exp(...)` - OpenCL kernel for fast GPU evaluation
  - Precompute G_0, v on CPU
  - GPU builds a_st per grid point using exponential coupling
  - Computes resp = |s1|^2 / |a_tt - s2|^2

### test_stm_orbital_projection.py (current test)
- Uses `response_amplitude_exp` GPU kernel
- Only computes response for exponential coupling
- Missing: true Fdata coupling, full block elimination on GPU

## 3. Implementation Plan

### Phase 1: CPU Reference with Full Augmented Matrix
**Goal**: Validate the physics by explicitly building and inverting the full (sample+tip) matrix on CPU.

1. **Build augmented H and S**
   - Start with H_s, S_s from Fireball (dense, n_s × n_s)
   - Add tip atom: H_t = [[E_tip]], S_t = [[1.0]] (1×1 for H 1s)
   - Build H_full (n_s+1 × n_s+1), S_full (n_s+1 × n_s+1)

2. **Compute coupling H_ts, S_ts per tip position**
   - For each grid point r_tip:
     - Compute distances to all sample atoms
     - Compute SK angular factors (direction cosines l, m, n)
     - Compute radial hopping: H_ts = f(r) × SK(l,m,n)
     - Set S_ts = 0 (or small value for vacuum)

3. **Compute full Green's function**
   ```
   A_full = (E + iη) S_full - H_full
   G_full = A_full^{-1}
   ```

4. **Extract response amplitude**
   - For single-orbital tip: resp = |G_full[tip, sample] · C|^2
   - Or equivalently: resp = |C^T x_s|^2 where x solves A x = b

5. **Compare with block elimination**
   - Verify that full inversion matches Dyson block elimination
   - This validates the mathematical formulation

### Phase 2: CPU Block Elimination (Optimized)
**Goal**: Efficient CPU implementation using precomputed G_0.

1. **Precompute once per energy**
   ```python
   A_ss = (E + iη) * S_s - H_s
   G_0 = np.linalg.inv(A_ss)
   v = C[:, mo].T @ G_0
   ```

2. **Per grid point loop**
   ```python
   for r_tip in grid_points:
       H_ts, S_ts = compute_coupling(r_tip)
       a_st = (E + iη) * S_ts - H_ts   # shape (1, n_s)
       a_vec = a_st[0, :]               # shape (n_s,)
       s1 = np.dot(v, a_vec.conj())
       s2 = np.dot(a_vec, G_0 @ a_vec.conj())
       a_tt = (E + iη) * 1.0 - E_tip
       x_tip = 1.0 / (a_tt - s2)
       resp[ip] = abs(x_tip)**2 * abs(s1)**2
   ```

3. **Implement both coupling models**
   - Fdata path: exact radial integrals
   - Exponential path: exp(-β(r-r0))

4. **Validate against Phase 1**
   - Compare full matrix inversion vs block elimination
   - Both should give identical results (numerical tests)

### Phase 3: GPU OpenCL Kernel with Full Block Elimination
**Goal**: Fast GPU implementation matching the CPU reference.

1. **Extend existing `response_amplitude_exp` kernel**
   - Current kernel builds a_st from exponential coupling
   - Already computes s1, s2, x_tip correctly
   - Need to add Fdata radial lookup path

2. **New kernel: `response_amplitude_fdata`**
   - Same structure as `response_amplitude_exp`
   - But builds a_st using Fdata radial integrals (spline lookup)
   - Requires uploading Fdata basis buffers to GPU

3. **Verify GPU vs CPU parity**
   - Run both paths on same grid
   - Compare max|resp_GPU - resp_CPU|
   - Target: < 1e-3 relative error

### Phase 4: Integration with test_stm_orbital_projection.py
**Goal**: Update the test script to show all four panels properly.

1. **Panel layout (2×3)**
   - Row 1: ψ(r) Fortran, ψ(r) OpenCL true basis, ψ(r) OpenCL exp tail
   - Row 2: Response (true coupling), Response (exp coupling), LDOS

2. **Compute all quantities consistently**
   - Use same tip position grid for all panels
   - Same energy E = ε_MO for all
   - Same broadening η

3. **Add LDOS panel for comparison**
   - Compute LDOS(r; E=ε_MO) from sample Green's function
   - Shows how response differs from simple LDOS

### Phase 5: Generalize to test_mo_vs_ldos.py
**Goal**: Add response function computation to the general MO vs LDOS test.

1. **Add `--tip` flag**
   - `--tip-type 1` (H atom)
   - `--tip-height 2.0` (z above molecule)
   - `--tip-energy 0.0` (E_tip)

2. **Compute response alongside LDOS**
   - For each MO: compute both LDOS and response
   - Add response panel to output figure

3. **Compare LDOS vs Response**
   - In Tersoff-Hamann limit: response ∝ LDOS
   - With tip orbital structure: response shows orbital selectivity (Chen's rules)

## 4. Key Equations Summary

**Full system matrix:**
```
A = (E + iη) S - H
```

**Block elimination:**
```
G_0 = A_ss^{-1}
x_tip = 1 / (a_tt - a_st^T G_0 a_st^H)
x_s = -G_0 a_st^H x_tip
resp = |C^T x_s|^2 = |x_tip|^2 |v · a_st^H|^2
```

**Coupling (exponential approximation):**
```
H_ts[μ, ν] = -exp(-β(|r_μ - R_ν| - r0)) × SK(l, m, n)
S_ts = 0
```

**Coupling (Fdata path):**
```
H_ts[μ, ν] = Fireball radial integral × SK(l, m, n)
S_ts[μ, ν] = Fireball overlap integral × SK(l, m, n)
```

## 5. Files to Modify

**New functions in STM.py:**
- `build_augmented_HS(H_s, S_s, H_t, S_t, H_ts, S_ts)` - stack matrices
- `compute_response_full_matrix(H_full, S_full, C, E, eta)` - explicit inversion
- `compute_response_block_elimination(H_s, S_s, H_ts_func, C, E, eta)` - optimized

**New GPU kernels in Grid.cl:**
- `response_amplitude_fdata` - Fdata-based coupling on GPU

**Modified tests:**
- `test_stm_orbital_projection.py` - add true coupling response
- `test_mo_vs_ldos.py` - add response function panel

## 6. Validation Checklist

- [ ] CPU full matrix = CPU block elimination (exact match)
- [ ] CPU block elimination (exp) = GPU kernel (exp) (parity < 1e-3)
- [ ] CPU block elimination (fdata) = GPU kernel (fdata) (parity < 1e-3)
- [ ] Response symmetry matches molecule symmetry
- [ ] Response for s-tip ≈ LDOS (Tersoff-Hamann limit)
- [ ] Response for p_z-tip shows dz/dz enhancement (Chen's rule)

## 7. Optimization: Avoiding Full Matrix Inversion per Grid Point

### 7.1 Problem Statement

For each tip position **r**, solve **A(r)·x = b** where **A = (E+iη)S(r) - H(r)**. Only the tip-sample coupling **H_st(r)** changes; **A_ss = (E+iη)S_s - H_s** is constant across the grid.

### 7.2 LU/Cholesky Precomputation of A_ss

Precompute **A_ss = L·U** once. Per grid point:
```
  y = U^{-1}(L^{-1} · a_st(r))    # forward/back substitution, O(n_s^2)
  s2 = a_st^T · y                  # O(n_s)
  x_tip = 1 / (a_tt - s2)          # O(1)
  x_s = -x_tip · y                 # O(n_s)
  resp = |C^T · x_s|^2
```
Complexity: **O(n_s^2)** per point (vs **O(n_s^3)** for full inversion). LU is more numerically stable than explicit G_0 = A_ss^{-1} for large systems.

### 7.3 Sparse Iterative Solvers

For large sparse A_ss (tight-binding Hamiltonian):
```
  A_ss_sparse = sparse CSR/CSC matrix
  M = ILU(A_ss)  # incomplete LU preconditioner (precomputed once)

  For each r:
    a_st_sparse = sparse vector (only atoms within rcut)
    y = gmres(A_ss_sparse, a_st_sparse, M=M, x0=x_prev, tol=1e-12)
```
**Advantages**: Memory O(nnz), cost O(nnz) per iteration. **Disadvantages**: Complex arithmetic needs BiCGSTAB/QMR; convergence depends on preconditioner quality.

### 7.4 GPU Warm-Start Iterative Solver

**Physical insight**: x_s(r) changes smoothly across the grid. Previous solution is an excellent initial guess.

**Algorithm (row-major scan)**:
```
  x_prev = 0
  for r in grid_points:
    a_st = build_coupling_gpu(r)          # GPU kernel
    y = sparse_matvec_iterative_gpu(A_ss, a_st, x0=x_prev)
    resp = compute_response_gpu(y, C)    # GPU kernel
    x_prev = y                            # warm start for next point
```
**Benefits**: Converges in ~2-5 iterations instead of 50-100 (for row-major order). Can use Jacobi or SSOR preconditioner computed once on GPU.

### 7.5 Schur Complement + Sparse Direct Solver (Recommended for Medium Systems)

For systems where n_s ~ 1000-10000:
1. Precompute sparse LU of A_ss (e.g., UMFPACK, SuperLU)
2. Per point: sparse forward/back solve with a_st(r)
3. Cost: ~O(nnz) per point after O(nnz) symbolic + O(nnz^1.5) numeric factorization

This is the **optimal CPU path** for systems too large for dense matrices but small enough for direct solvers.

### 7.6 Comparison Table

| Method | Precompute | Per Point | Best For |
|--------|-----------|-----------|----------|
| Dense inversion (G_0) | O(n_s^3) | O(n_s^2) | n_s < 500 |
| LU of A_ss | O(n_s^3) | O(n_s^2) | n_s < 2000 |
| Sparse direct (UMFPACK) | O(nnz^1.5) | O(nnz) | n_s < 10000 |
| Iterative + warm start | O(nnz) for M | O(k·nnz) | n_s > 10000 |
| GPU block elimination (current) | O(n_s^3) CPU | O(n_s^2) GPU | Dense on GPU |
| GPU sparse iterative | O(nnz) CPU | O(k·nnz) GPU | Large sparse on GPU |

### 7.7 Recommended GPU Path for FireCore

For DFTB with sparsity pattern from Fireball (localized orbitals, cutoff ~4-6 Å):

1. **Precompute on CPU**:
   - Sparse LU factorization of A_ss (symbolic + numeric)
   - Upload L, U factors to GPU as sparse CSR
   - Or: upload A_ss in CSR and use GPU sparse iterative library (cuSPARSE, clSPARSE)

2. **Per grid point on GPU**:
   - Build a_st(r) via coupling kernel (exponential or Fdata)
   - Solve A_ss · y = a_st using GPU sparse triangular solve
   - Compute response |C^T·y|^2 via reduction kernel

3. **Optimization**: Batch multiple grid points per kernel launch to amortize overhead.

### 7.8 Next Immediate Steps

1. **Phase 1**: CPU dense reference (full augmented matrix inversion) — validates math
2. **Phase 2**: CPU dense block elimination with explicit G_0 — baseline for parity
3. **Phase 3**: CPU LU-based block elimination — compare timing vs G_0 approach
4. **Phase 4**: GPU block elimination (current kernel) — verify parity with CPU
5. **Phase 5**: GPU sparse iterative with warm start — for large systems (future)
