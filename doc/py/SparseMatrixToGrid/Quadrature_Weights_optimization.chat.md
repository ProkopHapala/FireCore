https://aistudio.google.com/app/prompts?state=%7B%22ids%22:%5B%221uS4jr2NwrA-QJIimtJND-_vEWi4jWSlr%22%5D,%22action%22:%22open%22,%22userId%22:%22100958146796876347936%22,%22resourceKeys%22:%7B%7D%7D&usp=sharing

## USER

OK, this estimate of weights is good starting point. What I want to do now to find optimial weights for quadrature (intragral) of bread range of function on this grid. My idea is to generate set of functions (eigher exp(-r/w) or gauss exp(-(r/w)^2) with anular part as sp2 comonation (a + b*x + c*y) ) and compare the result with high-resolution regular rectangular grid reference.

Then I should define miminization problem that sum of square errora summed over all these test function is minimal. This should be linear optimization problem solvable by simple algrabra  - please write ow exactly. I would like to add some regularization - e.g. on difference of the new weights from initial guess of the weights. 

first analyse the problem and write math equations how to solve it, then write python script how to do it, which imports the  actuall script 
triangle_plug_gemini_tri2.py you implemented previously so you do not have to copy paste all functions


---

## Gemini

### Mathematical Formulation

We define the problem as a **Regularized Least Squares** optimization.

**1. Definitions:**
*   Let $\mathbf{r}_i$ be the coordinates of the $N$ grid points.
*   Let $w^{guess}_i$ be the initial area-based geometric weights we calculated previously.
*   Let $\{f_j(\mathbf{r})\}_{j=1}^M$ be a set of test functions (Gaussian/Slater orbitals with s, p symmetries).
*   Let $I^{ref}_j = \int f_j(\mathbf{r}) d\mathbf{r}$ be the "true" integral of the $j$-th function, calculated on a very high-resolution reference grid.

**2. The Linear System:**
We want to find new weights $\mathbf{w}$ such that the discrete sum matches the true integral:
$$ \sum_{i=1}^N w_i f_j(\mathbf{r}_i) \approx I^{ref}_j $$

This can be written in matrix form:
$$ \mathbf{F} \mathbf{w} \approx \mathbf{I}^{ref} $$
Where $\mathbf{F}$ is an $M \times N$ matrix with elements $F_{ji} = f_j(\mathbf{r}_i)$.

**3. The Cost Function:**
We want to minimize the integration error while keeping the weights "physical" (close to the geometric area estimate). We introduce a regularization parameter $\lambda$.

$$ J(\mathbf{w}) = \underbrace{\| \mathbf{F}\mathbf{w} - \mathbf{I}^{ref} \|_2^2}_{\text{Error}} + \underbrace{\lambda \| \mathbf{w} - \mathbf{w}^{guess} \|_2^2}_{\text{Regularization}} $$

**4. The Analytic Solution:**
To minimize $J$, we take the derivative with respect to $\mathbf{w}$ and set it to zero. This yields the normal equations (Ridge Regression form):

$$ \nabla_{\mathbf{w}} J = 2 \mathbf{F}^T (\mathbf{F}\mathbf{w} - \mathbf{I}^{ref}) + 2 \lambda (\mathbf{w} - \mathbf{w}^{guess}) = 0 $$

Rearranging to solve for $\mathbf{w}$:

$$ (\mathbf{F}^T \mathbf{F} + \lambda \mathbf{I}) \mathbf{w} = \mathbf{F}^T \mathbf{I}^{ref} + \lambda \mathbf{w}^{guess} $$

This is a linear system $\mathbf{A}\mathbf{x} = \mathbf{b}$ solvable by standard linear algebra:
*   $\mathbf{A} = \mathbf{F}^T \mathbf{F} + \lambda \mathbf{I}$
*   $\mathbf{b} = \mathbf{F}^T \mathbf{I}^{ref} + \lambda \mathbf{w}^{guess}$
*   $\mathbf{w} = \mathbf{A}^{-1} \mathbf{b}$

---

### Python Implementation

Save the previous code as `triangle_plug_gemini_tri2.py` in the same directory. This script imports it, generates the basis set, computes the reference integrals, and solves the linear system.

```python
#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import sys

# Import the geometry generator from the previous script
try:
    import triangle_plug_gemini_tri2 as gen
except ImportError:
    print("Error: Please save the previous code as 'triangle_plug_gemini_tri2.py'")
    sys.exit(1)

# ==========================================
# 1. Basis Function Definitions
# ==========================================

def get_basis_functions():
    """
    Returns a list of dictionaries defining the test functions.
    Types: 
      - Radial: Gauss (exp(-r^2/w^2)) or Slater (exp(-r/w))
      - Angular: s (1), px (x), py (y)
    """
    funcs = []
    
    # We span a range of widths to ensure the grid handles 
    # both sharp core peaks and diffuse tails.
    widths = [0.3, 0.5, 0.8, 1.2, 2.0] 
    
    for w in widths:
        for rad_type in ['gauss', 'slater']:
            # s-orbital (symmetric)
            funcs.append({'type': rad_type, 'w': w, 'ang': (1,0,0), 'name': f'{rad_type}_s_{w}'})
            # px-orbital (antisymmetric x)
            funcs.append({'type': rad_type, 'w': w, 'ang': (0,1,0), 'name': f'{rad_type}_px_{w}'})
            # py-orbital (antisymmetric y)
            funcs.append({'type': rad_type, 'w': w, 'ang': (0,0,1), 'name': f'{rad_type}_py_{w}'})
            
    return funcs

def eval_func(f_def, x, y):
    """Evaluate a single function at x, y."""
    r = np.sqrt(x*x + y*y)
    
    # Radial part
    if f_def['type'] == 'gauss':
        rad = np.exp(-(r/f_def['w'])**2)
    elif f_def['type'] == 'slater':
        rad = np.exp(-r/f_def['w'])
        
    # Angular part (1, x, or y)
    cx, cy, cz = f_def['ang'] # cz unused in 2D
    ang = cx + cy*x + cz*y
    
    return ang * rad

# ==========================================
# 2. Reference Integrator
# ==========================================

def compute_reference_integrals(funcs, half, margin_cells, h_coarse, resolution_factor=10):
    """
    Compute 'True' integrals using a very dense grid.
    resolution_factor: how much finer than the coarse grid.
    """
    # Define bounds roughly covering the coarse grid + margin
    limit = half + (margin_cells + 1) * h_coarse
    
    # Create dense grid
    h_ref = h_coarse / resolution_factor
    x = np.arange(-limit, limit + 1e-9, h_ref)
    y = np.arange(-limit, limit + 1e-9, h_ref)
    xv, yv = np.meshgrid(x, y)
    
    # Area element
    dA = h_ref * h_ref
    
    I_refs = []
    print(f"Computing reference integrals on {len(x)}x{len(y)} grid...")
    
    for f in funcs:
        vals = eval_func(f, xv, yv)
        total = np.sum(vals) * dA
        I_refs.append(total)
        
    return np.array(I_refs)

# ==========================================
# 3. Optimization Solver
# ==========================================

def optimize_weights(coords, w_guess, funcs, I_ref, lambda_reg=1e-4):
    """
    Solve (F.T @ F + lambda*I) @ w = F.T @ I_ref + lambda * w_guess
    """
    N = len(coords)
    M = len(funcs)
    
    print(f"Building Matrix F ({M} functions x {N} points)...")
    
    # 1. Build F matrix
    # Rows = Functions, Cols = Grid Points
    F = np.zeros((M, N))
    
    # To keep order consistent, we convert dict keys to a list
    point_keys = list(coords.keys())
    
    for j, f in enumerate(funcs):
        for i, k in enumerate(point_keys):
            pt = coords[k]
            F[j, i] = eval_func(f, pt[0], pt[1])
            
    # 2. Setup Linear System
    # w_guess vector aligned with point_keys
    w_guess_vec = np.array([w_guess[k] for k in point_keys])
    
    # A = F^T F + lambda * I
    print("Constructing normal equations...")
    FT = F.T
    A = FT @ F
    
    # Regularization: Add lambda to diagonal
    # Scaling lambda by trace ensures it's relative to matrix magnitude
    # or just use fixed small number if normalized.
    # Let's use simple fixed lambda.
    A_reg = A + lambda_reg * np.eye(N)
    
    # b = F^T I_ref + lambda * w_guess
    b = FT @ I_ref + lambda_reg * w_guess_vec
    
    # 3. Solve
    print("Solving linear system...")
    w_opt_vec = np.linalg.solve(A_reg, b)
    
    # 4. Map back to dictionary
    w_opt = {k: w_opt_vec[i] for i, k in enumerate(point_keys)}
    
    return w_opt, point_keys, F, w_guess_vec, w_opt_vec

# ==========================================
# 4. Main Routine
# ==========================================

def main():
    # --- Parameters ---
    n_sub = 4
    half = 1.0
    margin = 2
    dmin = 0.2
    dmax = 1.0
    reg_strength = 1.0e-2 # Strength of regularization (keep close to geometric)
    
    print("--- 1. Generating Grid Geometry ---")
    # Call functions from the previous script
    # Onion coords
    onion_map = gen.get_onion_coords(n_sub, dmin, dmax, half)
    onion_tris = gen.get_onion_topology(n_sub, onion_map)
    
    # Coarse coords
    coarse_quads = gen.get_coarse_grid(n_sub, half, margin)
    
    # Merge and get initial geometric weights
    unique_pts, w_geom, elements = gen.merge_and_analyze(onion_tris, coarse_quads)
    
    print(f"Total Grid Points: {len(unique_pts)}")
    
    # --- 2. Setup Optimization Problem ---
    funcs = get_basis_functions()
    print(f"Generated {len(funcs)} test functions.")
    
    # Calculate Coarse Grid Spacing (for reference resolution)
    h_coarse = (2.0 * half) / n_sub
    
    # Compute Reference Integrals (The "Truth")
    I_ref = compute_reference_integrals(funcs, half, margin, h_coarse, resolution_factor=8)
    
    # --- 3. Solve ---
    w_opt, ordered_keys, F, w_vec_geom, w_vec_opt = optimize_weights(
        unique_pts, w_geom, funcs, I_ref, lambda_reg=reg_strength
    )
    
    # --- 4. Analysis & Plotting ---
    
    # Calculate Errors
    I_geom = F @ w_vec_geom
    I_opt  = F @ w_vec_opt
    
    err_geom = np.abs(I_geom - I_ref)
    err_opt  = np.abs(I_opt - I_ref)
    
    print("\n--- Optimization Results ---")
    print(f"Mean Absolute Error (Geometric): {np.mean(err_geom):.6e}")
    print(f"Mean Absolute Error (Optimized): {np.mean(err_opt):.6e}")
    print(f"Error Reduction Factor:          {np.mean(err_geom)/np.mean(err_opt):.2f}x")
    
    # Plotting
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    # Plot 1: Weight Correction Map
    # (Optimized - Geometric)
    ax = axes[0]
    xs, ys, diffs = [], [], []
    for k in ordered_keys:
        p = unique_pts[k]
        xs.append(p[0])
        ys.append(p[1])
        diffs.append(w_opt[k] - w_geom[k])
        
    sc = ax.scatter(xs, ys, c=diffs, cmap='coolwarm', s=20, edgecolor='none')
    plt.colorbar(sc, ax=ax, label='Weight Correction (Opt - Geom)')
    ax.set_title(f"Weight Corrections\nlambda={reg_strength}")
    ax.set_aspect('equal')
    
    # Plot 2: Optimized Weights Map
    ax = axes[1]
    vals = [w_opt[k] for k in ordered_keys]
    sc2 = ax.scatter(xs, ys, c=vals, cmap='viridis', s=20)
    plt.colorbar(sc2, ax=ax, label='Optimized Weight')
    ax.set_title("Final Optimized Weights")
    ax.set_aspect('equal')
    
    # Plot 3: Error distribution (Test Functions)
    ax = axes[2]
    ax.plot(err_geom, label='Geometric Weights', marker='o', ls='none', alpha=0.5)
    ax.plot(err_opt, label='Optimized Weights', marker='x', ls='none', alpha=0.8)
    ax.set_yscale('log')
    ax.set_xlabel('Test Function Index')
    ax.set_ylabel('Absolute Integration Error')
    ax.set_title("Integration Error per Function")
    ax.legend()
    
    # --- Check Positivity ---
    min_w = np.min(w_vec_opt)
    print(f"\nMinimum Weight: {min_w:.6f} (Should be > 0)")
    if min_w < 0:
        print("WARNING: Some weights are negative. Increase lambda_reg.")
        
    # --- Check Symmetry ---
    # Check if weights on right edge match top edge
    mid_right = (half, 0.0)
    mid_top   = (0.0, half)
    k_r = tuple(np.round(mid_right, 10))
    k_t = tuple(np.round(mid_top, 10))
    
    if k_r in w_opt and k_t in w_opt:
        print(f"Symmetry Check (Opt): R={w_opt[k_r]:.5f}, T={w_opt[k_t]:.5f}")
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
```

### Explanation of the Solution

1.  **Objective Function**:
    The script minimizes $\| \mathbf{F}\mathbf{w} - \mathbf{I}_{ref} \|^2 + \lambda \| \mathbf{w} - \mathbf{w}_{geom} \|^2$.
    *   The first term forces the weights to integrate Gaussian/Slater functions correctly.
    *   The second term forces the weights to stay close to the "Area-based" weights we derived geometrically. This prevents "overfitting" where weights might become huge positive/negative numbers just to satisfy the equations.

2.  **Regularization ($\lambda$)**:
    I set `reg_strength = 1.0e-2`.
    *   If you make this smaller (e.g., $10^{-6}$), the error on the test functions will drop to almost zero, but the weights map (Plot 1) will look "noisy" or "checkerboarded".
    *   If you make this larger (e.g., $1.0$), the weights will hardly change from the geometric area calculation.
    *   $10^{-2}$ is usually a sweet spot where you get smoother weights and significant error reduction.

3.  **Basis Set**:
    I included `gauss` and `slater` types with `s`, `px`, and `py` symmetries.
    *   Including `px` and `py` is crucial. If you only optimized for `s` (spherical) functions, the solver might produce an asymmetric grid (e.g., giving all weight to the left side) because it wouldn't be penalized for getting $x$-dependent integrals wrong.

4.  **Results**:
    The output plot will show:
    *   **Left**: Where the optimization added or subtracted weight. You will likely see corrections at the "seam" (boundary of plug) and potentially near the core.
    *   **Right**: Comparison of errors. The orange crosses (Optimized) should be consistently lower on the log scale than the blue dots (Geometric).

---

## USER

OK this is pretty good ! thank you!  now can please try to consider inherent symmetry ? I mean we have square symmetry, therefore some points map to each other (you may find which by construction, or by applying symmetry operation and searching which overlap within certain limit. Note that all the original "coarse" grid regular grid points have same weights by definition (so they are equvalent by tranalation symmetry. So we have rather small number of unique points camming only from the "plug" and these are the coefs we want to optimize (we do not optimize the regular grid coefs). So the symetry search should be applied only on the plug. Can you do that please? 

Also in the end we sholed print info for each unique point - original weight, new weight. 

---

## Gemini

This is a significant improvement to the formulation. By enforcing symmetry constraints directly, we achieve two things:
1.  **Robustness**: The solution is guaranteed to be symmetric (rotation/reflection), eliminating numerical noise that might make the grid slightly asymmetric.
2.  **Efficiency**: We reduce the number of unknown variables by a factor of ~8 (since one coefficient controls 8 symmetric points).
3.  **Stability**: By fixing the "far-field" coarse points, the optimizer focuses entirely on correcting the plug and the seam.

### Mathematical Formulation with Symmetry & Constraints

We divide the grid points into two sets:
1.  **Fixed Set ($\mathcal{S}_{fix}$)**: Points strictly outside the plug (the regular grid). Their weights are fixed to the geometric area $w_i = h_{coarse}^2$.
2.  **Optimizable Set ($\mathcal{S}_{opt}$)**: Points belonging to the plug (internal + boundary/seam).

We group the optimizable points into **Symmetry Groups**. A group $k$ contains all points $\{r_i\}$ that map to each other via rotation or reflection (e.g., $(x,y)$, $(-x,y)$, $(y,x)$, etc.). All points in group $k$ share a single weight variable $x_k$.

**The Integral Equation:**
$$ \sum_{i \in \mathcal{S}_{fix}} w_i^{fix} f_j(\mathbf{r}_i) + \sum_{k=1}^{K} x_k \left( \sum_{i \in \text{group } k} f_j(\mathbf{r}_i) \right) \approx I_j^{ref} $$

Let's define a **Reduced Matrix** $\tilde{\mathbf{F}}$ and a **Modified Target** $\tilde{\mathbf{b}}$:

1.  **Reduced Basis Matrix** ($\tilde{\mathbf{F}}$):
    Element $\tilde{F}_{jk}$ is the sum of the basis function evaluations for all points in symmetry group $k$:
    $$ \tilde{F}_{jk} = \sum_{i \in \text{group } k} f_j(\mathbf{r}_i) $$

2.  **Modified Target Vector** ($\tilde{\mathbf{b}}$):
    Subtract the known contribution of the fixed points from the reference integral:
    $$ \tilde{b}_j = I_j^{ref} - \sum_{i \in \mathcal{S}_{fix}} w_i^{fix} f_j(\mathbf{r}_i) $$

**The Optimization Problem:**
Find vector $\mathbf{x}$ (size $K$) that minimizes:
$$ J(\mathbf{x}) = \| \tilde{\mathbf{F}}\mathbf{x} - \tilde{\mathbf{b}} \|^2 + \lambda \| \mathbf{x} - \mathbf{x}_{guess} \|^2 $$

**The Solution:**
$$ \mathbf{x} = (\tilde{\mathbf{F}}^T \tilde{\mathbf{F}} + \lambda \mathbf{I})^{-1} (\tilde{\mathbf{F}}^T \tilde{\mathbf{b}} + \lambda \mathbf{x}_{guess}) $$

---

### Python Implementation

Save this as `triangle_plug_optimization.py`. It requires your previous geometry script saved as `triangle_plug_gemini_tri2.py`.

```python
#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import sys

# Import geometry generator
try:
    import triangle_plug_gemini_tri2 as gen
except ImportError:
    print("Error: Save the geometry code as 'triangle_plug_gemini_tri2.py'")
    sys.exit(1)

# ==========================================
# 1. Symmetry & Grouping Logic
# ==========================================

def get_canonical_key(x, y, decimals=9):
    """
    Map any point (x,y) to its fundamental symmetry domain 
    (First octant: x >= 0, y >= 0, y <= x).
    """
    # 1. Reflection symmetry (abs)
    ax, ay = abs(x), abs(y)
    # 2. Diagonal reflection (sort)
    mn, mx = min(ax, ay), max(ax, ay)
    return (round(mn, decimals), round(mx, decimals))

def group_points(unique_pts, half_size, fixed_margin_threshold=1e-9):
    """
    Classify points into:
    1. Fixed Points (strictly outside plug).
    2. Symmetry Groups (inside or on boundary of plug).
    
    Returns:
        fixed_indices: list of point keys
        opt_groups: list of dicts {'key': canonical, 'indices': [pt_keys], 'radius': r}
        pt_to_group_map: dict mapping point_key -> group_index (or -1 if fixed)
    """
    fixed_indices = []
    groups_map = {} # canonical_key -> list of point_keys
    
    for k, pt in unique_pts.items():
        x, y = pt
        
        # Check if strictly outside the plug square [-half, half]
        # Using a small epsilon for float tolerance
        if max(abs(x), abs(y)) > half_size + fixed_margin_threshold:
            fixed_indices.append(k)
        else:
            # It belongs to the plug (or seam) -> Optimize it
            ckey = get_canonical_key(x, y)
            if ckey not in groups_map:
                groups_map[ckey] = []
            groups_map[ckey].append(k)
            
    # Convert groups map to a sorted list (e.g., sorted by radius)
    # This ensures reproducible ordering of coefficients
    opt_groups = []
    for ckey, indices in groups_map.items():
        # Representative point for radius
        r = np.sqrt(ckey[0]**2 + ckey[1]**2)
        opt_groups.append({
            'key': ckey,
            'indices': indices,
            'radius': r,
            'count': len(indices)
        })
        
    # Sort groups by radius (center outward)
    opt_groups.sort(key=lambda g: g['radius'])
    
    return fixed_indices, opt_groups

# ==========================================
# 2. Basis Functions & Reference (Reused)
# ==========================================

def get_basis_functions():
    funcs = []
    widths = [0.3, 0.5, 0.8, 1.2, 2.0]
    for w in widths:
        for rad_type in ['gauss', 'slater']:
            # s, px, py needed to constrain symmetry
            funcs.append({'type': rad_type, 'w': w, 'ang': (1,0,0)}) # s
            funcs.append({'type': rad_type, 'w': w, 'ang': (0,1,0)}) # px
            funcs.append({'type': rad_type, 'w': w, 'ang': (0,0,1)}) # py
    return funcs

def eval_func(f_def, x, y):
    r = np.sqrt(x*x + y*y)
    if f_def['type'] == 'gauss': val = np.exp(-(r/f_def['w'])**2)
    else: val = np.exp(-r/f_def['w'])
    return val * (f_def['ang'][0] + f_def['ang'][1]*x + f_def['ang'][2]*y)

def compute_reference_integrals(funcs, half, margin, h_coarse):
    # Dense grid numerical integration
    limit = half + (margin + 0.5) * h_coarse
    h_ref = h_coarse / 8.0 # 8x subsampling
    x = np.arange(-limit, limit + 1e-9, h_ref)
    xv, yv = np.meshgrid(x, x)
    dA = h_ref**2
    
    I_refs = []
    print(f"Computing Ref Integrals (Grid {len(x)}x{len(x)})...")
    for f in funcs:
        vals = eval_func(f, xv, yv)
        I_refs.append(np.sum(vals) * dA)
    return np.array(I_refs)

# ==========================================
# 3. Symmetry-Constrained Optimization
# ==========================================

def optimize_symmetric(unique_pts, w_geom, funcs, I_ref, half, lambda_reg=1e-2):
    
    # 1. Group Points
    print("Grouping points by symmetry...")
    fixed_keys, opt_groups = group_points(unique_pts, half)
    
    n_groups = len(opt_groups)
    n_fixed = len(fixed_keys)
    n_funcs = len(funcs)
    
    print(f"  Fixed Points (Coarse): {n_fixed}")
    print(f"  Optimizable Groups:    {n_groups} (covering {len(unique_pts) - n_fixed} points)")
    
    # 2. Build Reduced System
    # We need F_reduced (Rows=Funcs, Cols=Groups)
    # And we need to subtract fixed point contributions from I_ref
    
    F_tilde = np.zeros((n_funcs, n_groups))
    contribution_fixed = np.zeros(n_funcs)
    
    # Pre-compute fixed contribution
    print("Computing fixed point contributions...")
    for i_f, k in enumerate(fixed_keys):
        pt = unique_pts[k]
        w = w_geom[k]
        for j, f in enumerate(funcs):
            val = eval_func(f, pt[0], pt[1])
            contribution_fixed[j] += w * val
            
    # Pre-compute reduced matrix columns
    print("Building reduced matrix...")
    for k_grp, group in enumerate(opt_groups):
        indices = group['indices']
        # Sum basis functions for all points in this group
        for idx_pt in indices:
            pt = unique_pts[idx_pt]
            for j, f in enumerate(funcs):
                F_tilde[j, k_grp] += eval_func(f, pt[0], pt[1])

    # 3. Prepare Optimization Targets
    b_tilde = I_ref - contribution_fixed
    
    # Initial guess for the groups (average of geometric weights in group)
    x_guess = np.zeros(n_groups)
    for k_grp, group in enumerate(opt_groups):
        ws = [w_geom[k] for k in group['indices']]
        x_guess[k_grp] = np.mean(ws)
        
    # 4. Solve Linear System (Ridge Regression)
    # (F.T @ F + lam * I) x = F.T @ b + lam * x_guess
    print("Solving linear system...")
    A = F_tilde.T @ F_tilde
    A_reg = A + lambda_reg * np.eye(n_groups)
    rhs = F_tilde.T @ b_tilde + lambda_reg * x_guess
    
    x_opt = np.linalg.solve(A_reg, rhs)
    
    # 5. Distribute weights back to full grid
    w_final = w_geom.copy() # Start with all geometric
    
    # Weights for fixed points stay as they were (area-based)
    # Weights for optimized points are updated
    for k_grp, group in enumerate(opt_groups):
        val = x_opt[k_grp]
        for k in group['indices']:
            w_final[k] = val
            
    return w_final, opt_groups, x_guess, x_opt, F_tilde, b_tilde

# ==========================================
# 4. Main & Reporting
# ==========================================

def main():
    # --- Config ---
    n_sub = 4
    half = 1.0
    margin = 3 # Ensure enough coarse points to stabilize integration
    dmin, dmax = 0.2, 1.0
    reg_strength = 0.005
    
    # --- Generate Geometry ---
    print("--- 1. Generating Geometry ---")
    onion_map = gen.get_onion_coords(n_sub, dmin, dmax, half)
    onion_tris = gen.get_onion_topology(n_sub, onion_map)
    coarse_quads = gen.get_coarse_grid(n_sub, half, margin)
    unique_pts, w_geom, _ = gen.merge_and_analyze(onion_tris, coarse_quads)
    
    # Calculate coarse area h^2 for reference
    h = (2.0*half)/n_sub
    h2 = h*h
    
    # --- Physics & Optimization ---
    print("\n--- 2. Setting up Optimization ---")
    funcs = get_basis_functions()
    I_ref = compute_reference_integrals(funcs, half, margin, h)
    
    w_final, groups, x_geom, x_opt, F_red, b_mod = optimize_symmetric(
        unique_pts, w_geom, funcs, I_ref, half, lambda_reg=reg_strength
    )
    
    # --- Error Analysis ---
    # Recompute errors
    # Error = || F_red * x + contribution_fixed - I_ref ||
    # Which is equivalent to || F_red * x - b_tilde ||
    err_geom = np.abs(F_red @ x_geom - b_mod)
    err_opt  = np.abs(F_red @ x_opt - b_mod)
    
    print("\n--- 3. Results ---")
    print(f"Error Reduction: {np.mean(err_geom)/np.mean(err_opt):.2f}x")
    print(f"Mean Abs Error (Geom): {np.mean(err_geom):.6e}")
    print(f"Mean Abs Error (Opt):  {np.mean(err_opt):.6e}")
    
    # --- Detailed Table ---
    print("\n" + "="*85)
    print(f"| {'Radius':<8} | {'Canonical (x, y)':<20} | {'Mult':<4} | {'Old W':<10} | {'New W':<10} | {'Change %':<8} |")
    print("="*85)
    
    # Groups are already sorted by radius
    for i, g in enumerate(groups):
        k = g['key']
        count = g['count']
        old_w = x_geom[i]
        new_w = x_opt[i]
        chg = (new_w - old_w)/old_w * 100
        
        # Helper to denote if it's seam or core
        note = ""
        if g['radius'] < 1e-9: note = "(CORE)"
        elif abs(max(k) - half) < 1e-9: note = "(SEAM)"
        
        print(f"| {g['radius']:<8.4f} | {str(k):<20} | {count:<4} | {old_w:<10.6f} | {new_w:<10.6f} | {chg:>+7.2f}% | {note}")

    print("="*85)
    
    # --- Visual Check ---
    fig, ax = plt.subplots(figsize=(8,8))
    
    # Extract data for plotting
    xs, ys, cs = [], [], []
    
    # Plot corrections
    for k, p in unique_pts.items():
        xs.append(p[0])
        ys.append(p[1])
        # Correction amount
        cs.append(w_final[k] - w_geom[k])
    
    # Plot Seam Box
    rect = plt.Rectangle((-half, -half), 2*half, 2*half, fill=False, ec='black', ls='--')
    ax.add_patch(rect)
    
    sc = ax.scatter(xs, ys, c=cs, cmap='coolwarm', s=30)
    plt.colorbar(sc, label='Weight Correction (Opt - Geom)')
    
    ax.set_aspect('equal')
    ax.set_title(f"Symmetry-Constrained Weight Corrections\nRegularization={reg_strength}")
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
```

### Interpretation of Output

1.  **The Table**:
    *   **Radius**: Distance from center. Helps identify layers.
    *   **Canonical (x,y)**: The representative coordinate for the symmetry group (always $0 \le y \le x$).
    *   **Mult**: Multiplicity.
        *   Core = 1.
        *   Axis points $(x, 0)$ = 4 (rotational copies).
        *   Diagonals $(x, x)$ = 4.
        *   General points $(x, y)$ = 8 (4 rotations $\times$ 2 reflections).
    *   **Old W / New W**: Compares the area-based estimate to the optimized quadrature weight.
    *   **Change %**: You will likely see:
        *   Small changes ($<1\%$) for internal plug points.
        *   Larger changes (maybe $5-15\%$) for **SEAM** points. This is expected because the simple area summation doesn't perfectly capture the overlap of the infinite tails of the Gaussian/Slater functions across the geometry change.

2.  **Robustness**:
    Notice that the **Regular Grid points are not in the table**. They were subtracted from the reference integral beforehand. The optimizer *only* adjusted the plug and seam to match the remaining integral deficit. This guarantees that your far-field matches the standard grid definition exactly.