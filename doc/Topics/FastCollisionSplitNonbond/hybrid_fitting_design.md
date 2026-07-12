# Hybrid Potential Fitting: Design Document

## Goal

Approximate the exact reference potential `V_ref(r)` (Morse+Coulomb from all atoms, large supercell/Ewald)
by a model potential:

```
V_model(r) = V_bspline(r) + V_atomic(r)
```

where:
- `V_bspline(r)` = tricubic B-spline on a regular 3D grid (smooth, captures long-range)
- `V_atomic(r) = Σ_j f_j(|r - r_j|)` = sum of per-atom radial basis functions (short-range, localized)

We want to **simultaneously fit** both the B-spline control points and the atomic basis parameters
to minimize `||V_ref - V_model||²` over the evaluation grid.

---

## Current State (Problems)

1. **Gaussian smoothing proxy**: `fit_bspline_3d` uses `scipy.ndimage.gaussian_filter` as a crude
   proxy for B-spline fitting. This is NOT a fit — it's a low-pass filter. It cannot adapt to capture
   features at the B-spline resolution scale.

2. **No simultaneous fitting**: The atomic basis (folded Morse) is computed analytically and subtracted
   from V_ref, then the residual is smoothed. The atomic basis parameters (D_e, R_eq, alpha) are fixed
   from the force field, not fitted. The B-spline has no knowledge of the atomic basis.

3. **No optimization**: There is no least-squares or optimization step anywhere. The "fit" is just
   `V_bspline = smooth(V_ref - V_atomic)` with fixed atomic parameters.

---

## Proposed Architecture

### Step 1: Define the model

```
V_model(r_i) = Σ_k c_k * B_k(r_i)  +  Σ_j Σ_m p_{jm} * φ_{jm}(|r_i - r_j|)
```

where:
- `B_k(r)` = tricubic B-spline basis function centered at grid node k
- `c_k` = B-spline control point values (to be fitted)
- `φ_{jm}(r)` = m-th radial basis function for atom j (e.g. Morse×window, exponential, polynomial)
- `p_{jm}` = atomic basis coefficients (to be fitted)
- `r_i` = i-th evaluation grid point

### Step 2: Build the design matrix

The model is **linear in all parameters** (c_k and p_{jm}):

```
V_model = A @ x
```

where `x = [c_0, c_1, ..., c_{N_bspline}, p_{0,0}, p_{0,1}, ..., p_{N_atom,N_basis}]`

and `A` is the design matrix with columns:
- First N_bspline columns: B-spline basis values `B_k(r_i)` at each grid point
- Next columns: atomic basis values `φ_{jm}(|r_i - r_j|)` at each grid point

### Step 3: Solve the linear least-squares problem

```
minimize ||A @ x - V_ref||²  +  λ * ||L @ x||²
```

where:
- `||A @ x - V_ref||²` = data misfit
- `λ * ||L @ x||²` = regularization (e.g. smoothness penalty on B-spline, or Tikhonov)
- Solve via `scipy.sparse.linalg.lsqr` or `scipy.linalg.lstsq` (sparse A for large grids)

### Step 4: Evaluate model and compute error

```
V_model = A @ x
error = V_ref - V_model
```

---

## Implementation Details

### B-spline basis construction

Use `scipy.interpolate.BSpline` or `scipy.ndimage` to construct tricubic B-spline basis functions
on a regular 3D grid with node spacing `Δs` (e.g. 1.0 Å).

For each grid point `r_i` and each B-spline node `k`, compute `B_k(r_i)`:
- A tricubic B-spline has compact support: only nodes within 2 cells affect each point
- The basis matrix is **sparse** (each point depends on ~4³=64 nodes)
- Use `scipy.sparse` for efficiency

### Atomic basis functions

For each atom j (top layer only), use radial basis functions with compact support:

**Option A: Windowed Morse (current)**
```
φ_j(r) = Morse_j(r) * w(r)    where w(r) = (1 - r/R_fold)^(2n)
```
Parameters: D_e, R_eq, alpha (from force field, or fitted)

**Option B: Exponential basis (more general)**
```
φ_{jm}(r) = A_m * exp(-α_m * r) * w(r)    for m = 1..M
```
Parameters: A_m, α_m (fitted), with multiple exponentials per atom type

**Option C: Polynomial basis**
```
φ_{jm}(r) = r^m * w(r)    for m = 0..M
```
Parameters: coefficients (fitted)

### Regularization

- **B-spline smoothness**: penalize `∇²c_k` (Laplacian of control points)
- **Atomic basis**: optional Tikhonov to prevent overfitting
- Or use cross-validation to select λ

### PBC handling

- B-spline: use periodic boundary conditions in x,y (wrap mode)
- Atomic basis: sum over periodic images within R_fold cutoff

---

## Algorithm

```
1. Compute V_ref on 3D grid (all atoms, large supercell/Ewald)
2. Select atomic basis atoms (top layer) and basis functions
3. Build sparse design matrix A:
   - B-spline columns: B_k(r_i) for each node k, grid point i
   - Atomic columns: φ_{jm}(|r_i - r_j|) for each atom j, basis m, grid point i
4. Build regularization matrix L (optional)
5. Solve: x = lsqr(A, V_ref, λ, L)  (sparse least-squares)
6. Extract: c_k (B-spline coefficients), p_{jm} (atomic coefficients)
7. Evaluate: V_model = A @ x, error = V_ref - V_model
8. Plot: reference, model, error, B-spline component, atomic component, 1D scans
```

---

## Key Differences from Current Code

| Current | Proposed |
|---------|----------|
| `gaussian_filter` proxy | Real B-spline basis matrix |
| Fixed atomic params | Fitted atomic params (simultaneous) |
| Sequential (subtract, then smooth) | Simultaneous (one LSQ solve) |
| No optimization | Linear least-squares with regularization |
| ~0.87 eV RMSE | Expected much lower (proper fit) |

---

## Suggested Implementation Path

1. **Phase 1**: Replace `gaussian_filter` with real B-spline fitting using `scipy.interpolate.BSpline`
   - Build basis matrix for 1D, extend to 3D (tensor product)
   - Solve least-squares for control points only (no atomic basis yet)
   - Verify error decreases vs Gaussian proxy

2. **Phase 2**: Add atomic basis columns to the design matrix
   - Start with fixed Morse parameters (just add as known component)
   - Solve for B-spline + atomic scale factors

3. **Phase 3**: Fit atomic basis parameters simultaneously
   - Use variable projection (VARPRO): for nonlinear params (α_m, R_eq),
     solve inner linear LSQ for B-spline + linear atomic coeffs,
     outer nonlinear optimization for nonlinear params
   - Or use full nonlinear optimization (scipy.optimize)

4. **Phase 4**: Add regularization, cross-validation, error analysis

---

## File Structure

```
pyBall/OCL/Surface_utils.py:
  - fit_hybrid_3d(V_ref, atoms, basis_config, grid_config) → (V_model, coeffs, error)
    - Builds design matrix (sparse)
    - Solves regularized LSQ
    - Returns fitted model and error
  - eval_hybrid_model(coeffs, atoms, eval_points) → V_model
    - Evaluates fitted model at arbitrary points

tests/tMMFF/test_hybrid_potential.py:
  - Calls fit_hybrid_3d
  - Plots results (same 4-figure layout)
```

---

## Expected Outcome

With proper simultaneous fitting:
- B-spline captures smooth long-range (Coulomb-like) part
- Atomic basis captures short-range (Morse-like) part near atoms
- Error should be dominated by:
  - B-spline resolution (1.0 Å nodes → ~0.01-0.1 eV for smooth parts)
  - Atomic basis expressiveness (how well Morse×window represents the true short-range)
  - Near-atom singularity (1/r Coulomb diverges — atomic basis must handle this)
