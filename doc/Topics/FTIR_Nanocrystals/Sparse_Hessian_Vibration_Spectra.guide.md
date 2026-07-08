---
description: Sparse Hessian computation and vibration spectrum analysis for nanocrystals
---

# Sparse Hessian and Vibration Spectrum Analysis

**Docs:** [`README.md`](README.md) · [`doc/topical_audit/Nanocrystal_Vibrations.md`](../../topical_audit/Nanocrystal_Vibrations.md) · [`tests/tSiNCs/README.md`](../../../tests/tSiNCs/README.md)

Reference for computing sparse Hessian matrices from neighbor shells and analyzing vibration spectra of nanocrystals using multiple methods (dense eigenmode, sparse frequency scanning, and comparison).

## Overview

Force field Hessians are naturally sparse because interactions are short-ranged. Instead of computing the full dense 3N×3N matrix (which scales as O(N²) in storage and O(N²) per force evaluation), we extract only the 3×3 blocks for atoms within a few neighbor shells. This reduces both the force evaluation cost (fewer forces read per perturbation) and memory.

We then use the sparse blocks for vibration spectrum analysis via three methods:

1. **Dense eigenmode method** (`vibration_spectrum_from_modes`): Diagonalize once, then evaluate response at all frequencies analytically. Fast for systems up to ~3000 DOF.
2. **Sparse frequency scanning** (`mechanical_greens_probing_sparse`): Build sparse BSR matrix, solve per frequency using `scipy.sparse.linalg.spsolve`. Suitable for medium systems up to ~10000 DOF.
3. **Both modes**: Run both methods and compare results for validation.

## Key Components

### C++: Sparse Hessian Block Computation

**File:** `cpp/libs/Molecular/MMFF_lib.cpp`

```cpp
void getHessianSparseBlocks(int natoms, int max_neigh, int* neigh_idx, 
                            int* neigh_counts, double* out_blocks, double dx)
```

- Computes 3×3 force constant blocks for each atom-neighbor pair using central finite differences
- Reads forces only for listed neighbors (not all atom pairs)
- Temporarily disables nonbonded interactions for speed
- Writes to Python-allocated contiguous array: shape `(natoms, max_neigh, 3, 3)`

### Python: Sparse Hessian Infrastructure

**File:** `pyBall/MMFF.py`

- `buildNeighShells(neighs, n_shells=2, max_neigh=64, include_self=True)`: BFS shell builder from 1st-neighbor list with deduplication, handles 0/1-based indexing
- `getHessianSparseBlocks(neigh_idx, neigh_count, dx=1e-4)`: ctypes wrapper to C++ function
- `getBuffs()`: Exports `neighs` buffer as `(natoms, 4)` for full atom set including capping atoms

**File:** `pyBall/FTIR.py`

- `build_sparse_hessian_from_blocks(neigh_idx, neigh_count, blocks, symmetrize=True)`: Builds scipy.sparse BSR matrix from sparse blocks
- `mechanical_greens_probing_sparse(H_sparse, M, omegas, ...)`: Per-frequency sparse linear solves using `spsolve`
- `vibration_spectrum_from_modes(K, M, omegas, ...)`: Fast eigenmode-based spectrum (dense diagonalization + analytical response)
- `mechanical_greens_probing(K, M, omegas, ...)`: Original dense per-frequency solve (slow, kept for reference)

### Test Scripts

**File:** `tests/tMMFF/test_nanocrystal_sparse_hessian.py`

Main test script for hydrogen-passivated diamond nanocrystals:

- Generates diamond nanoparticle with spherical cut
- Relaxes geometry with bonds+angles only (nonbonded/pi disabled)
- Computes sparse Hessian blocks via neighbor shells
- Reconstructs dense Hessian for diagonalization (optional)
- Computes vibration spectrum via chosen method
- Saves sparse blocks, spectrum data, and optional plots

**CLI options:**
```bash
--R <float>                 # Nanoparticle radius (default 10.0)
--nrep <int>                # Replication count (default 8)
--shells <int>              # Neighbor shell depth (default 2)
--max-neigh <int>           # Max neighbors per atom (default 64)
--dx <float>                # Finite difference step (default 1e-4)
--spectrum                  # Compute vibration spectrum
--spectrum-method {dense,sparse,both}  # Computation method
--fmin <float>              # Min frequency (default 0.01)
--fmax <float>              # Max frequency (default 50.0)
--nfreq <int>               # Number of frequency points (default 2000)
--eta <float>               # Damping parameter (default 0.05)
--diag                      # Diagonalize dense Hessian
--plot                      # Plot Hessian visualization
--plot-spectrum             # Plot spectrum
--outdir <path>             # Output directory (default OUT_nanoparticle)
--no-plot                   # Skip plotting (avoid ASan+matplotlib crashes)
```

**File:** `tests/tMMFF/plot_sparse_hessian.py`

Standalone plotting script for saved sparse Hessian `.npz` files:
- Plots sparse rectangular block heatmap
- Optionally plots dense reconstructed Hessian heatmap with log scale

## Performance Characteristics

| Method | Complexity | Best For | Notes |
|--------|------------|----------|-------|
| Dense eigenmode | O(N³) + O(nfreq×N) | N < 3000 DOF | Fastest for small systems |
| Sparse spsolve | O(nfreq×nnz) | N < 10000 DOF | Scales with sparsity, per-ω solve |
| Dense per-ω solve | O(nfreq×N³) | None | Too slow, kept for reference |

For 558-atom diamond nanocrystal (1674 DOF):
- Dense eigenmode: ~1-2 seconds
- Sparse spsolve (100 freq points): ~10-20 seconds
- Dense per-ω solve (800 freq points): ~minutes (appears frozen)

## Usage Examples

### Generate nanoparticle and compute spectrum (dense eigenmode)
```bash
cd /home/prokop/git/FireCore/tests/tMMFF
bash run.sh test_nanocrystal_sparse_hessian.py \
  --R 8.0 --nrep 6 --shells 2 --max-neigh 48 --dx 1e-4 \
  --spectrum --spectrum-method dense \
  --fmin 0.01 --fmax 10.0 --nfreq 100 --eta 0.05
```

### Compute spectrum with sparse solver
```bash
bash run.sh test_nanocrystal_sparse_hessian.py \
  --R 8.0 --nrep 6 --shells 2 --max-neigh 48 --dx 1e-4 \
  --spectrum --spectrum-method sparse \
  --fmin 0.01 --fmax 10.0 --nfreq 100 --eta 0.05
```

### Compare both methods
```bash
bash run.sh test_nanocrystal_sparse_hessian.py \
  --R 8.0 --nrep 6 --shells 2 --max-neigh 48 --dx 1e-4 \
  --spectrum --spectrum-method both \
  --fmin 0.01 --fmax 10.0 --nfreq 100 --eta 0.05 \
  --plot-spectrum
```

### Plot saved sparse Hessian
```bash
python plot_sparse_hessian.py \
  --npz OUT_nanoparticle/diamond_nc_R8.0_sparse_hessian.npz \
  --plot-dense
```

## Output Files

- `diamond_nc_R{R:.1f}_init.xyz`: Initial nanoparticle structure
- `diamond_nc_R{R:.1f}_relaxed.xyz`: Relaxed structure after MMFF optimization
- `diamond_nc_R{R:.1f}_sparse_hessian.npz`: Sparse Hessian blocks (neigh_idx, neigh_count, blocks)
- `diamond_nc_R{R:.1f}_spectrum_dense.npz`: Spectrum from dense eigenmode method
- `diamond_nc_R{R:.1f}_spectrum_sparse.npz`: Spectrum from sparse solver method
- `diamond_nc_R{R:.1f}_spectrum.png`: Spectrum plot (if --plot-spectrum)

## Implementation Notes

### Sparse Hessian Block Layout

Blocks are stored as `(natoms, max_neigh, 3, 3)` where:
- `blocks[p, j]` = force constant matrix for atom `p` with neighbor `neigh_idx[p, j]`
- Block `H[o,p] = -dF_o/du_p` (negative force derivative)
- Self-neighbor (`include_self=True`): `blocks[p, self_idx]` contains `H[p,p]`, the self-response block
- Symmetric: `H[p,o] = H[o,p].T`

### Central Finite Difference Formula

For each atom `p` and coordinate `k`:
```
H[o,p][l,k] = -(F_o.l(+dx) - F_o.l(-dx)) / (2*dx)
```
where `F_o.l(+dx)` is the `l`-th force component on neighbor `o` when atom `p` is displaced by `+dx` in direction `k`.

### Neighbor Shell Construction

- BFS extends from 1st-neighbor list (MMFF.neighs)
- `n_shells=2` includes 1st and 2nd neighbors
- Deduplicates to avoid repeated entries
- Handles 1-based indexing from MMFF (converts to 0-based internally)

### Reconstructing Dense Hessian from Blocks

The test script includes a `reconstruct_dense_from_blocks()` helper that assembles the full 3N×3N matrix from sparse blocks. This is useful for validation and for methods that require dense matrices (eigenmode-based spectrum). The reconstruction:
1. Places each block `H[o,p]` at the correct `(o*3:o*3+3, p*3:p*3+3)` position
2. Adds transposed blocks for symmetry
3. For diagonal blocks (self-neighbor), places `H[p,p]` on the diagonal

### Rigid Mode Projection

The `project_rigid_modes()` function shifts eigenvalues of 6 rigid body modes (3 translations + 3 rotations) to high frequency to prevent spurious low-energy response.

### Why the Eigenmode Method is So Much Faster

Dense diagonalization is O(N³) **once**. Per-ω dense solves are O(N³) **per frequency point**. For 800 frequency points, that's 800× slower. The eigenmode method diagonalizes once, then evaluates the analytical response at all frequencies in O(nfreq×N). The sparse method avoids dense solves but still requires O(nfreq) sparse linear solves.

### Mass Matrix

Mass matrix is diagonal: `M = diag(m_i, m_i, m_i)` for each atom `i`. Masses are taken from MMFF atom types.

## Dependencies

- NumPy (always required)
- SciPy (required for sparse methods: `scipy.sparse`, `scipy.sparse.linalg.spsolve`)
- Matplotlib (optional, for plotting)

## Known Issues

- **ASan + Matplotlib crashes**: When compiling with Address Sanitizer (ASan), importing matplotlib can cause crashes. Use `--no-plot` to skip plotting, or use the standalone `plot_sparse_hessian.py` script separately after the computation finishes.
- **Sparse spsolve on complex matrices**: `scipy.sparse.linalg.spsolve` may have numerical issues near resonances (where `det(A(ω)) ≈ 0`). The `stabilize` parameter (default 1e-6) adds a small identity shift to improve conditioning.

## Known Bugs and Limitations

### Hydrogen Placement in Generated Nanocrystals

The initial hydrogen passivation in generated SiNC structures produces **incorrect H positions**. The current placement logic does not properly account for surface orientation and bond angles, leading to clashing or poorly positioned capping atoms. This must be fixed before reliable spectra can be computed from generated structures. For now, structures should be visually inspected and manually corrected, or generated using an external tool with proper passivation.

## Future Extensions

- **Move neighbor shell calculation to C++**: The `buildNeighShells()` BFS extension is currently in Python. For large systems, this should be moved to C++ to avoid Python overhead and large array copies.
- **Fast non-bonded repulsion for passivation clashes**: When surface passivation groups clash, add a fast non-bonded interaction using **bounding-box overlap detection** rather than slow long-range electrostatics. This provides Pauli repulsion without full MMFF non-bonded evaluation cost.
- **Harmonic spring alternative to non-bonded repulsion**: For linear response spectra, instead of full non-bonded interactions, detect atom pairs that are too close and add temporary harmonic springs (bonds) to enforce minimum separation. This mimics Pauli repulsion (LJ/Morse wall) at much lower computational cost and keeps the Hessian sparse.
- Partial eigenspectrum via `scipy.sparse.linalg.eigsh` for larger systems (>10000 DOF)
- GPU-accelerated sparse linear solves for very large systems
- Preconditioned iterative solvers for better scaling
- Mode analysis (IR/Raman activity from dipole/derivative tensors)
