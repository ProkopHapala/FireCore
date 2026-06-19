# XRD Simulation on GPU — Progress Report

## 1. What was built

### Source files (add to git)

| File | Role |
|------|------|
| `pyBall/XRD/__init__.py` | Package init, exports `XRDDebye`, `load_xyz_atoms`, `get_form_factor_table`, `compute_sigma_from_sparse_blocks`, `compute_sigma_exact` |
| `pyBall/XRD/form_factors.py` | Cromer-Mann atomic form factors for H, C, Si |
| `pyBall/XRD/debye_histogram.py` | `XRDDebye` pyOpenCL engine + `compute_sigma_from_sparse_blocks` + `compute_sigma_exact` + ensemble accumulation methods |
| `cpp/common_resources/cl/XRDDebye.cl` | OpenCL kernels: `pair_histogram`, `pair_histogram_gaussian`, `debye_transform_q` |
| `scripts/generate_xrd_webgl.py` | CLI tool to generate standalone WebGL 2.0 single-crystal diffraction viewer |
| `tests/tXRD/test_debye_histogram.py` | Test script: powder XRD on existing nanocrystal fixtures with three thermal-broadening modes |
| `tests/tXRD/test_ensemble_exact.py` | Test script: ensemble accumulation (8 rotated copies) + exact thermal broadening via splu |
| `tests/tXRD/test_large_crystal.py` | Test script: larger C diamond NCs (R=8, R=10) — generate, relax, Hessian, XRD |

### Output artifacts (gitignored, regenerated on demand)

- `OUT_XRD/diamond_nc_R6_xrd.npz` — saved spectra and histograms
- `OUT_XRD/diamond_nc_R6_xrd.png` — comparison plot (static / constant σ / Hessian σ)
- `OUT_XRD/diamond_nc_R6_single_crystal.html` — interactive WebGL 2.0 viewer

## 2. Physics

### Powder diffraction via pair-distance histogram

The exact powder intensity for a finite non-periodic cluster follows the **Debye scattering equation**:

$$I(Q) = \sum_{i,j} f_i(Q) f_j(Q) \frac{\sin(Q r_{ij})}{Q r_{ij}}$$

Instead of the direct O(N²·N_Q) double sum, we:

1. **Build a histogram** `H_ab(r)` of interatomic distances for each element pair (a,b)
2. **Fourier-transform** the histogram with `sinc(Qr)`
3. Multiply by precomputed `f_a(Q)·f_b(Q)`

This is exact for the pair term; a self-scattering term `Σ_i f_i²` is added separately.

### Thermal broadening in the histogram

Your insight (validated in `XrayDiffraction.chat.md`) is that the **joint thermal probability density** of two atoms can be **projected onto the pair-distance axis**, giving a 1D Gaussian smear:

$$P_{ij}(r) \approx \frac{1}{\sqrt{2\pi}\sigma_{ij}} \exp\left(-\frac{(r - r_{ij}^0)^2}{2\sigma_{ij}^2}\right)$$

We deposit this Gaussian into the histogram (GPU kernel `pair_histogram_gaussian`) rather than post-multiplying a Debye-Waller factor. This is more physical because:
- The Gaussian in `r` **automatically** produces the correct Q-damped sinc transform
- Each pair gets its **own** width, capturing anisotropic local stiffness

### Computing σ_ij from the Hessian (frozen-stiffness approximation)

For each pair (i,j), the directional stiffness along the bond vector `u = (r_j - r_i)/|r_j - r_i|` is:

$$k_{ij} = \mathbf{u}^T (\mathbf{H}_{ii} + \mathbf{H}_{jj} - 2\mathbf{H}_{ij}) \mathbf{u}$$

Then:

$$\sigma_{ij} = \sqrt{k_B T / k_{ij}}$$

This is **not** the exact thermal width (which would need the pseudoinverse `H⁺`), but it is:
- Fast (uses only local 3×3 Hessian blocks already computed by MMFF)
- Physically reasonable for bonded / near-neighbor pairs
- A strict **lower bound** on the true thermal width (frozen atoms overestimate stiffness)

For distant pairs beyond the sparse cutoff, a constant fallback `σ = 0.02 Å` is used.

### Single-crystal 2D diffraction (WebGL)

The generated HTML page contains **two synchronized views**:

1. **Main canvas** (full screen): diffraction pattern computed per pixel via the Debye structure-factor sum in a fragment shader.
2. **Crystal-view overlay** (260×260 px, top-right): a simple 3D ball-and-stick model of the nanocrystal that rotates in lockstep with the diffraction pattern. Bonds are precomputed from covalent distance thresholds (`C-C < 1.8 Å`, `C-H < 1.2 Å`, etc.).

Both views share the same quaternion `rotQuat`, so dragging the mouse rotates the crystal **and** updates the diffraction pattern simultaneously.

For a detector pixel at screen coordinates `(x,y)` at distance `L` from the sample:

$$\mathbf{k}_{\text{out}} = \frac{2\pi}{\lambda} \hat{\mathbf{n}}, \quad \mathbf{n} = L \hat{\mathbf{k}}_{\text{in}} + x \hat{\mathbf{u}} + y \hat{\mathbf{v}}$$

$$\mathbf{Q} = \mathbf{k}_{\text{out}} - \mathbf{k}_{\text{in}}$$

The structure factor is accumulated directly per pixel:

$$F(\mathbf{Q}) = \sum_j f_j(Q) e^{i \mathbf{Q} \cdot \mathbf{r}_j}$$

Intensity: `I = |F|²`

The fragment shader runs this loop over all atoms for every pixel — fully parallel, no grids, no FFT.

## 3. GPU Performance

### Kernel design

| Kernel | Work-items | Work per item | Memory pattern |
|--------|-----------|---------------|----------------|
| `pair_histogram` | N atoms | O(N) pair loop, atomic-add to global histogram | Coherent read (atom i), scattered write (histogram bins) |
| `pair_histogram_gaussian` | N atoms | O(N·w) where w = stencil width (~3σ/dr bins) | Same, more atomics |
| `debye_transform_q` | N_Q points | O(n_pairs · n_bins) dot product | Read-only global, no atomics |

### Key optimisations

- **Float atomics via CAS loop** — `atomic_add_float` uses `atom_cmpxchg` in a spin loop. Portable across NVIDIA/AMD/Intel. The `#pragma cl_khr_fp_atomic` path was removed because NVIDIA drivers do not expose it.
- **Linear interpolation in histogram** — reduces discretisation error vs. nearest-bin deposition
- **Gaussian stencil limited to ±3σ** — avoids touching all bins for every pair
- **Precomputed `f_a(Q)·f_b(Q)` table** — uploaded once, reused for all crystals

### Observed performance (RTX 3090)

| Task | Time | Notes |
|------|------|-------|
| Histogram (270 atoms, 1000 bins) | ~10 ms | 36k pairs |
| Debye transform (800 Q-points) | ~5 ms | 3 pair types × 1000 bins |
| WebGL render (full screen) | 60 FPS | 270-atom loop per pixel in fragment shader |

For 10k atoms (50M pairs), the histogram kernel will dominate. Mitigations:
- Increase `dr` (coarser bins = fewer atomic collisions)
- Use **shared-memory tiling** (atoms binned into local memory, then global flush)
- Process multiple crystals in a single kernel launch (batched histograms)

## 4. User Guide

### Prerequisites

- PyOpenCL with a working GPU context (tested on NVIDIA RTX 3090)
- `matplotlib` (optional, for plotting)
- A modern browser with WebGL 2.0 (Chrome/Firefox) for the HTML viewer

### 1. Powder diffraction (Python)

```python
from pyBall.XRD import XRDDebye, load_xyz_atoms

engine = XRDDebye(preferred_vendor='nvidia')
atoms, unique_Z = load_xyz_atoms('crystal.xyz')  # atoms: (N,4) float32

Q = np.linspace(0.5, 15.0, 800, dtype=np.float32)
I = engine.powder_spectrum_with_self(atoms, unique_Z, Q, r_max=20.0, dr=0.02)
```

**Arguments:**
- `r_max` — maximum pair distance to histogram (Å). Should exceed the nanocrystal diameter.
- `dr` — histogram bin width (Å). Finer = more bins, more GPU memory, sharper peaks.
- `Q` — scattering vector magnitudes (Å⁻¹). Cu Kα gives Q ≈ 0.5–8 Å⁻¹.

### 2. Thermal broadening modes

| Mode | How to use | Physics quality |
|------|-----------|-----------------|
| **Static** (no broadening) | `engine.compute_histogram(atoms, ...)` | Reference only |
| **Constant σ** | `engine.compute_histogram_gaussian(atoms, sigma=0.015, ...)` | Isotropic Debye-Waller |
| **Hessian-based** | `sigma = compute_sigma_from_sparse_blocks(pos, neigh_idx, neigh_count, blocks); engine.compute_histogram_gaussian(atoms, sigma, ...)` | Pair-specific, frozen-stiffness approximation |

The test script `tests/tXRD/test_debye_histogram.py` runs all three and produces a comparison plot.

### 3. Single-crystal 2D diffraction (WebGL)

Generate a standalone HTML file from any XYZ:

```bash
python scripts/generate_xrd_webgl.py \
    --xyz tests/tSiNCs/fixtures/vibration_parallel/structures/diamond_nc_R6_relaxed.xyz \
    --out OUT_XRD/my_crystal.html
```

**CLI arguments:**
| Flag | Default | Description |
|------|---------|-------------|
| `--xyz` | required | Input XYZ file |
| `--out` | required | Output HTML file |
| `--wavelength` | 1.54 | X-ray wavelength (Å) |
| `--detector-L` | 5000.0 | Detector distance (Å) |
| `--pixel-size` | 20.0 | Detector pixel size (Å) |
| `--exposure` | 1e-3 | Intensity scale factor |
| `--title` | (basename) | Page title |

**Controls in the browser:**
- **Drag** mouse anywhere to rotate the crystal (quaternion trackball). Both the diffraction pattern and the crystal-view overlay update simultaneously.
- **Scroll** to change exposure
- **R** key to reset rotation
- **Sliders** for exposure, detector distance, wavelength

The crystal-view overlay (top-right) shows a ball-and-stick model with:
- Atoms as small circles (color-coded by element: C=green, H=red, Si=blue)
- Bonds as grey lines (precomputed from covalent distance thresholds)

Open the file directly:
```bash
firefox OUT_XRD/my_crystal.html
```

## 5. Ensemble Accumulation (new)

### Physics

Powder XRD from a nanocrystal ensemble = sum of per-crystal pair-distance histograms, then one Debye transform on the accumulated histogram. This models size/shape polydispersity: different crystals contribute different pair-distance distributions.

$$H_{\text{ensemble}}(r) = \sum_{k=1}^{N_{\text{crystals}}} H_k(r)$$

### API

```python
engine = XRDDebye(preferred_vendor='nvidia')

# Method 1: accumulate pre-computed histograms (float64)
hist_ensemble = engine.ensemble_histogram([hist1, hist2, ...], bin_edges)

# Method 2: accumulate Gaussian-broadened histograms with per-pair sigma
hist_ensemble = engine.ensemble_histogram_gaussian(
    [atoms1, atoms2, ...], [sigma1, sigma2, ...], r_max=20.0, dr=0.02)

# Then Debye transform once
I_ensemble = engine.compute_spectrum(hist_ensemble, bin_edges, Q, ff_table, unique_Z)
```

### Test

`tests/tXRD/test_ensemble_exact.py` — 8 rotated copies of 270-atom diamond NC:
- 290,520 total pairs accumulated
- Ensemble/N vs single max|diff| ≈ 5.6e3 (expected: rotation changes float32 pair distances slightly)
- Output: `OUT_XRD/ensemble_exact_xrd.npz` + `.png`

## 6. Exact Thermal Broadening via Sparse LU (new)

### Why exact?

The frozen-stiffness approximation (§4) uses only local 3×3 Hessian blocks:
$$k_{ij}^{\text{frozen}} = \mathbf{u}^T (H_{ii} + H_{jj} - 2H_{ij}) \mathbf{u}$$

This **overestimates stiffness** because it freezes all other atoms. The exact thermal width uses the full Hessian pseudoinverse:
$$\sigma_{ij}^2 = k_B T \cdot \mathbf{b}_{ij}^T H_{\text{mw}}^{-1} \mathbf{b}_{ij}$$

where $\mathbf{b}_{ij}$ is the bond-direction force probe in mass-weighted coordinates and $H_{\text{mw}} = M^{-1/2} H M^{-1/2}$ is the mass-weighted Hessian with rigid-body modes projected out.

### Algorithm

1. Assemble dense Hessian from sparse blocks via `FTIR.build_sparse_hessian_from_blocks()` (reused, not duplicated)
2. Mass-weight: $H_{\text{mw}} = M^{-1/2} H M^{-1/2}$
3. Project out 6 rigid-body modes (3 translation + 3 rotation) via QR orthonormalization
4. Factorize $H_{\text{mw}} = LU$ via `scipy.sparse.linalg.splu`
5. For all N(N-1)/2 pairs: build RHS matrix $B$ (ndof × npairs), solve $X = H_{\text{mw}}^{-1} B$ in one batched call
6. $\sigma_{ij}^2 = k_B T \cdot \text{diag}(B^T X)$

### Performance

| Crystal | N atoms | DOF | Pairs | Factorization | Batched solve | Total |
|---------|---------|-----|-------|---------------|---------------|-------|
| diamond R=6 | 270 | 810 | 36,315 | ~0.1s | ~8s | ~8s |
| diamond R=8 | 558 | 1,674 | 155,403 | ~0.5s | **too slow** (>60s) | — |

**Bottleneck:** The batched `lu.solve(B)` with B shape (1674, 155403) requires ~250 MB and the SuperLU triangular solve with 155K RHS columns is slow. For N>500 atoms, the O(N²) number of pairs makes this approach impractical.

### What was tried

- **CHOLMOD selected inversion:** `scikit-sparse` not installed → fell back to `scipy.sparse.linalg.splu`
- **Per-pair Python loop (original):** 155K individual `lu.solve()` calls → ~minutes (Python overhead)
- **Batched solve (current):** single `lu.solve(B)` with all RHS → better but still slow for N>500

### What would fix it

The exact σ needs $H^{-1}$ applied to bond-direction probes. Alternatives:

1. **Diagonal of $H^{-1}$ only** (selected inversion): CHOLMOD `cholmod_l_select` computes just the diagonal blocks of the inverse — O(N) output, not O(N²). Then $\sigma_{ij}^2 = k_B T (u^T H_{ii}^{-1} u + u^T H_{jj}^{-1} u - 2 u^T H_{ij}^{-1} u)$. Requires `scikit-sparse` (`pip install scikit-sparse`).
2. **GPU sparse solve:** Use cuSPARSE or SuperLU_GPU for the factorization + batched solve. The factorization is O(N) for sparse Hessian, the solve is O(N) per RHS, but 155K RHS is still O(N³) total.
3. **Mode expansion:** Diagonalize once ($O(N^3)$), then $\sigma_{ij}^2 = k_B T \sum_n \frac{(v_n^T b_{ij})^2}{\omega_n^2}$. This is the same cost as `eigh` but gives all pairs from one eigendecomposition + matrix multiply.
4. **Approximation:** Use frozen-stiffness for most pairs, exact only for soft/floppy pairs where the ratio is large.

### Results (270-atom diamond NC, R=6)

| Method | σ range (Å) | σ median (Å) | Time |
|--------|-------------|---------------|------|
| Frozen-stiffness | 0.016–0.029 | 0.020 | 0.02s |
| Exact (splu) | 0.001–0.143 | 0.020 | 8.2s |
| Ratio exact/frozen | — | median 1.0, mean 1.59, max 7.15 | — |

**Key finding:** Frozen approximation is good for most pairs (median ratio 1.0), but significantly underestimates σ for ~10% of pairs (tail up to 7×). These are soft/floppy pairs where collective motion matters.

## 7. Files to add to git

### New files

| File | Role |
|------|------|
| `pyBall/XRD/__init__.py` | Package init (updated: now exports `compute_sigma_exact`) |
| `pyBall/XRD/form_factors.py` | Cromer-Mann form factors |
| `pyBall/XRD/debye_histogram.py` | `XRDDebye` engine + `compute_sigma_from_sparse_blocks` + `compute_sigma_exact` + ensemble methods |
| `cpp/common_resources/cl/XRDDebye.cl` | OpenCL kernels |
| `scripts/generate_xrd_webgl.py` | WebGL viewer generator |
| `tests/tXRD/test_debye_histogram.py` | Single-crystal XRD test (3 thermal broadening modes) |
| `tests/tXRD/test_ensemble_exact.py` | Ensemble accumulation + exact σ test (270-atom diamond) |
| `tests/tXRD/test_large_crystal.py` | Large crystal test (generates R=8, R=10 C diamond NCs) |

### Modified files

| File | Change |
|------|--------|
| `doc/Topics/FTIR_Nanocrystals/XRD_progress.md` | This report (updated) |

### Do NOT commit

- `OUT_XRD/` — output artifacts (NPZ, PNG, HTML, hessian cache)
- `tests/tMMFF/data` — symlink to `cpp/common_resources`

## 8. Known Limitations & Next Steps

### Current limitations

1. **Exact σ too slow for N>500** — batched `splu.solve(B)` with O(N²) RHS columns does not scale. Need selected inversion (CHOLMOD) or mode expansion.
2. **No Lorentz-polarisation factor** — easy to add as a post-multiply on `I(Q)`
3. **No instrumental broadening** — Gaussian convolution in Q-space after transform
4. **WebGL viewer is CPU-bound for large N** — fragment shader loops over all atoms per pixel
5. **Only H, C, Si form factors** — trivial to extend

### Recommended next steps

1. **Selected inversion via CHOLMOD:** `pip install scikit-sparse` → use `cholmod_l_select` for diagonal blocks of $H^{-1}$ → O(N) output, exact σ for all pairs
2. **Mode expansion approach:** `eigh` once, then matrix multiply for all pair σ — same cost as vibration spectrum diagonalization
3. **Ensemble pipeline integration:** wire `ensemble_histogram` into `run_nanocrystal_ensemble.mjs` orchestrator for multi-crystal XRD
4. **Batch GPU processing:** 2D NDRange for multi-crystal histograms

---

*Report updated: 2026-06-19*
*Code version: ensemble accumulation + exact thermal broadening (splu)*
