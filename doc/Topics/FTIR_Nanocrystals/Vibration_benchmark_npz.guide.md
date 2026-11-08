---
description: NPZ format for MMFF vibration Hessian benchmarks (external sparse solver repos)
---

# Vibration Hessian Benchmark `.npz` Format

Frozen **MMFF** stiffness matrices \(K\) and diagonal mass matrices \(M\) for benchmarking sparse linear solvers, eigenvalue methods, and frequency-domain probes \(A(\omega)=K-(\omega+i\eta)^2 M\).

**Generate locally (gitignored):**

```bash
python3 tests/tSiNCs/export_vibration_benchmarks.py
```

Output directory: `tests/tSiNCs/fixtures/vibration_benchmarks/` (see `README.md` there after export).

**Do not commit** generated `.npz` files; copy to your external solver benchmark repository manually.

---

## Visualization (no Hessian recomputation)

```bash
python3 tests/tSiNCs/plot_vibration_benchmark_sparsity.py
```

Writes alongside each `.npz`:

| PNG | Matrix |
|-----|--------|
| `{name}_Kproj_logabs.png` | Projected solver matrix `K_csr` / `H_dense_projected` — `imshow(log10\|H_ij\|)` |
| `{name}_Hblocks_logabs.png` | Shell-block Hessian assembled from `blocks` only (pre rigid-mode shift) |

**Note:** `Kproj` is **structurally dense** (100% fill) because rigid-mode projection couples all DOFs. The **physical** 1–2 / 1–3 shell truncation is visible in `Hblocks` (~2–30% fill, decreasing with system size).

---

1. **Positive-definite vibrational subspace:** Raw cluster Hessians have 6 zero eigenvalues (rigid body). Exported `K` includes **rigid-mode penalty** (`rigid_shift=1e6`) via `FTIR.apply_rigid_mode_shift`. Check `n_negative_projected == 0` before shipping.
2. **Sparse primary format:** Block list from MMFF `getHessianSparseBlocks` plus assembled **CSR** of projected \(K\).
3. **Dense reference (small systems):** `H_dense_projected` included when `ndof ≤ 3000` for exact `numpy.linalg.eigh` validation.
4. **Self-documenting:** Every file contains `meta_json` (UTF-8 JSON string) with human-readable description.

---

## Size guidance (dense `eigh` on CPU)

Run calibration:

```bash
python3 tests/tSiNCs/export_vibration_benchmarks.py --calibrate-only
```

Typical on a modern CPU (order of magnitude):

| DOF (3N) | Atoms N | Dense `eigh` time | Role |
|----------|---------|-------------------|------|
| 78 | 26 | &lt;0.01 s | Adamantane — exact reference |
| 300–800 | 100–270 | 0.05–0.3 s | Small NC |
| ~810 | 270 | ~0.02 s | `nc_C_R6` |
| ~1242 | 414 | ~0.06 s | `nc_C_R7` |
| ~1674 | 558 | ~0.12 s | **Largest with dense `eigh` &lt;1 s on reference CPU** |
| &gt;3000 | &gt;1000 | &gt;1 s | Sparse solvers only; omit dense arrays in `.npz` |

Target: largest system with **`eigh_seconds < 1 s`** is the ceiling for bundled dense matrices.

---

## File naming

| File | System |
|------|--------|
| `adamantane.npz` | C₁₀H₁₆, no relax |
| `nc_C_R4.npz` | C diamond NC, R=4 Å sphere |
| `nc_C_R5.npz` | R=5 Å |
| `nc_C_R6.npz` | R=6 Å |
| `nc_C_R7.npz` | R=7 Å |
| `nc_C_R8.npz` | R=8 Å (largest dense-reference case) |
| `si_G1_caps_only.npz` | Si NC from JS generator (if fixture exists) |

---

## NPZ array reference

### Metadata

| Key | Type | Description |
|-----|------|-------------|
| `meta_json` | str | Full `BenchmarkMeta` JSON (name, description, units, notes) |
| `name` | str | Short identifier |

### Geometry & mass

| Key | Shape | Description |
|-----|-------|-------------|
| `natoms` | scalar int32 | N |
| `ndof` | scalar int32 | 3N |
| `pos` | (N, 3) float64 | Å, MMFF geometry used for Hessian |
| `elements` | (N,) object str | Element symbols |
| `mass_diag` | (3N,) float64 | Diagonal of M (amu repeated per xyz) |

### Sparse block format (MMFF native)

Force constant blocks \(H_{o,p} = -\partial F_o / \partial u_p\) (3×3):

| Key | Shape | Description |
|-----|-------|-------------|
| `neigh_idx` | (N, max_neigh) int32 | Neighbor atom index per shell slot (-1 pad) |
| `neigh_count` | (N,) int32 | Active neighbors per atom |
| `blocks` | (N, max_neigh, 3, 3) float64 | Raw FD blocks **before** rigid shift |
| `n_shells` | scalar | BFS shell depth used to build `neigh_idx` |
| `max_neigh` | scalar | Second dimension of `neigh_idx` |
| `dx` | scalar | Finite-difference step for Hessian |

**Note:** Rigid-body projection is **non-local**; use CSR arrays below for the projected stiffness \(K\).

### Projected stiffness CSR (use this for solvers)

| Key | Description |
|-----|-------------|
| `K_csr_data` | float64 nonzeros |
| `K_csr_indices` | int32 column indices |
| `K_csr_indptr` | int32 row pointers |
| `K_csr_shape` | (3N, 3N) |

Rebuild in SciPy:

```python
import numpy as np, scipy.sparse as sp, json
d = np.load("adamantane.npz", allow_pickle=True)
meta = json.loads(str(d["meta_json"]))
K = sp.csr_matrix((d["K_csr_data"], d["K_csr_indices"], d["K_csr_indptr"]), shape=tuple(d["K_csr_shape"]))
M = sp.diags(d["mass_diag"])
```

### Reference spectrum (validation)

| Key | Description |
|-----|-------------|
| `omegas_modes_projected` | √(eigvals) of \(M^{-1/2} K_{proj} M^{-1/2}\) |
| `eigvals_mw_projected` | Eigenvalues (internal units) |
| `omega_min_vib` | Smallest ω above `1e-4` threshold (skips rigid penalty modes) |
| `n_negative_raw` | Negative eigvals before rigid shift (expect ≥6 for clusters) |
| `n_negative_projected` | Must be **0** |
| `eigh_seconds` | Time to compute reference `eigh` on this machine |
| `recon_vs_full_max_abs` | ‖K_shell_blocks − K_full‖∞ (shell completeness check) |

### Dense optional (small systems)

| Key | Shape |
|-----|-------|
| `H_dense_full` | (3N, 3N) raw MMFF FD |
| `H_dense_projected` | (3N, 3N) with rigid shift |
| `H_dense_recon_blocks` | (3N, 3N) assembled from `blocks` only |

---

## Dynamic stiffness probe (external repo)

For frequency \(\omega\) and damping \(\eta\):

\[
A(\omega) = K_{proj} - (\omega + i\eta)^2 M
\]

Solve \(A(\omega) u = f\) with your iterative method; compare dipole response or \(\|u\|\) to FireCore `FTIR.mechanical_greens_probing_sparse` on the same fixtures.

Suggested test RHS: unit charge on all atoms, direction \(\hat{x}\):

```python
f = np.zeros(3*N)
for i in range(N):
    f[3*i] = 1.0
```

---

## Units

Internal MMFF units (documented in `meta_json.units`):

- Hessian \(K\): eV/Å²  
- Mass: amu  
- \(\omega\) from `omegas_modes_projected`: \(\sqrt{\text{eV}/(\text{amu}·Å²)}}\) — **not** cm⁻¹  

Convert to cm⁻¹ only if you apply a documented conversion factor for your benchmark paper.

---

## Related

- Export code: `pyBall/vibration_benchmark.py`, `tests/tSiNCs/export_vibration_benchmarks.py`
- Rigid modes: `pyBall/FTIR.py` → `apply_rigid_mode_shift`
- H cap geometry (avoid 53° bug): use `pyBall/nanocrystal_gen.py` VSEPR, **not** legacy `test_nanocrystal_sparse_hessian` bisector formula

---

*Last updated: 2026-06-15*
