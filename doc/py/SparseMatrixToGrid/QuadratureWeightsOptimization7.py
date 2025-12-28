#!/usr/bin/env python3
"""
Quadrature weight optimization on a hybrid plug + coarse grid.

Goal
----
Find quadrature weights w_i that match high-res reference integrals of a diverse
test basis (Gauss/Slater orbitals with s/px/py angular parts) while staying close
to geometric area guesses and respecting square symmetry grouping.

Mathematical problem
--------------------
For functions f_j(r):
    sum_i w_i f_j(r_i) ≈ I_ref,j

Group symmetry-related plug points into K groups with coefficients x_k; coarse
points are fixed to h_coarse^2. After subtracting fixed coarse contribution:
    F̃ x ≈ b̃
with F̃_{jk} = sum_{i in group k} f_j(r_i), b̃_j = I_ref,j - sum_{i in coarse} w_i^fix f_j(r_i)

We solve a regularized least squares (ridge):
    x* = argmin ||F̃ x - b̃||² + λ ||x - x_guess||²
which yields normal equations:
    (F̃ᵀ F̃ + λ I) x = F̃ᵀ b̃ + λ x_guess
We also report coarse-only errors (using the original coarse grid without plug replacement) and
geometric vs optimized errors.

Reference integrals
-------------------
Computed on a dense rectangular grid (spacing h_ref = h_coarse / ref_factor):
    I_ref,j ≈ Σ f_j(x,y) * h_ref²

Test basis
----------
Probes: Gaussian bumps centered at representative plug points, widths scaled by nearest-neighbor
distance × probe_multipliers.
Random orbitals: Gauss or Slater with widths sampled log-uniform in [w_min, w_max]; angular part
uses (1, x/r, y/r) for Slaters (normalized) to avoid excessive spread.

Notes on pitfalls we hit and fixes
----------------------------------
- Slater angular part must be normalized (x/r, y/r); using raw x,y caused huge tails and errors.
- Width sampling must be bounded: very small widths were under-resolved by the grid; very large widths
  extended beyond the integration domain. CLI w_min/w_max now set the log-uniform range explicitly.
- Probe multiplicity: using several scaled widths per symmetry point improves local coverage; controlled
  by --probe-mults.
- Regularization is essential: ridge term (lambda-reg) keeps weights near geometric guesses and prevents
  ill-conditioned spikes; without it errors blew up.
- Coarse-only vs plug: we report coarse-only errors to reveal the baseline deficit before plug replacement.

"""
import argparse
import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.optimize import lsq_linear

# Import geometry generator
try:
    import triangle_plug_gemini_tri2 as gen
except ImportError:
    print("Error: Save the geometry code as 'triangle_plug_gemini_tri2.py'")
    sys.exit(1)

# ==========================================
# 1. Symmetry & Grouping
# ==========================================

def get_canonical_key(x, y, decimals=9):
    ax, ay = abs(x), abs(y)
    mn, mx = min(ax, ay), max(ax, ay)
    return (round(mn, decimals), round(mx, decimals))

def group_points(unique_pts, half_size, fixed_margin_threshold=1e-9):
    fixed_indices = []
    groups_map = {} 
    
    for k, pt in unique_pts.items():
        x, y = pt
        if max(abs(x), abs(y)) > half_size + fixed_margin_threshold:
            fixed_indices.append(k)
        else:
            ckey = get_canonical_key(x, y)
            if ckey not in groups_map: groups_map[ckey] = []
            groups_map[ckey].append(k)
            
    opt_groups = []
    for ckey, indices in groups_map.items():
        r = np.sqrt(ckey[0]**2 + ckey[1]**2)
        opt_groups.append({
            'key': ckey, 'indices': indices, 'radius': r, 'count': len(indices)
        })
    opt_groups.sort(key=lambda g: g['radius'])
    
    point_to_gid = {}
    for gid, g in enumerate(opt_groups):
        for k in g['indices']: point_to_gid[k] = gid
    for k in fixed_indices: point_to_gid[k] = -1
    
    return fixed_indices, opt_groups, point_to_gid

def fmt_key(k):
    return f"({float(k[0]):.6f}, {float(k[1]):.6f})"

def get_nearest_neighbor_dist(pt, all_pts_dict):
    min_d2 = 1e9
    x, y = pt
    for k, p2 in all_pts_dict.items():
        dx = x - p2[0]; dy = y - p2[1]
        d2 = dx*dx + dy*dy
        if d2 > 1e-12 and d2 < min_d2:
            min_d2 = d2
    return np.sqrt(min_d2)

# ==========================================
# 2. Basis Functions (Density Mode)
# ==========================================

def generate_basis_functions(opt_groups, unique_pts, n_random=100, probe_mults=None, w_min=0.2, w_max=2.5):
    funcs = []
    if probe_mults is None or len(probe_mults) == 0:
        probe_mults = [1.0]
    
    # 1. Probes (One per symmetry group)
    for i, g in enumerate(opt_groups):
        rep_key = g['indices'][0]
        rep_pt = unique_pts[rep_key]
        d_nn = get_nearest_neighbor_dist(rep_pt, unique_pts)
        
        for m in probe_mults:
            funcs.append({
                'type'     : 'probe', 
                'w'        : d_nn * m, 
                'center'   : rep_pt,
                'ang'      : (1,0,0), 
                'rad_type' : 'gauss',
                'gid'      : i
            })
        
    # 2. Random Orbitals
    for i in range(n_random):
        # Log-uniform width [w_min, w_max]
        w = np.exp(np.random.uniform(np.log(w_min), np.log(w_max)))
        # Random angular vector
        vec = np.random.normal(0, 1, 3)
        vec /= np.linalg.norm(vec)
        rtype = 'gauss' if np.random.rand() > 0.5 else 'slater'
        
        funcs.append({
            'type': 'orbital',
            'rad_type': rtype,
            'w': w,
            'ang': tuple(vec),
            'center': (0.0, 0.0)
        })
    return funcs

def eval_func_raw(f, x, y):
    """Evaluates the raw wavefunction value (can be negative)."""
    cx, cy = f['center']
    dx = x - cx; dy = y - cy
    r = np.sqrt(dx*dx + dy*dy)
    
    if f['rad_type'] == 'gauss' or f['type'] == 'probe':
        rad = np.exp(-(r/f['w'])**2)
    else:
        rad = np.exp(-r/f['w'])
        
    if f['type'] == 'probe':
        return rad
    else:
        # Use centered coordinates; for slater use normalized angular (x/r, y/r) to avoid excessive spread
        if f['rad_type'] == 'slater':
            inv_r = 1.0 / (r + 1e-12)
            ax = dx * inv_r
            ay = dy * inv_r
        else:
            ax = dx
            ay = dy
        ang = f['ang'][0] + f['ang'][1]*ax + f['ang'][2]*ay
        return rad * ang

def eval_density(f, x, y):
    """Evaluates the DENSITY (square of wavefunction). Always Positive."""
    val = eval_func_raw(f, x, y)
    return val * val

def compute_normalization_factors(funcs, half, margin, h_coarse, ref_factor=10):
    """
    Computes integral of psi^2 for each function.
    Returns: NormFactors (= 1/Integral)
    """
    limit = half + (margin + 0.5) * h_coarse
    h_ref = h_coarse / ref_factor
    
    # Check resolution
    min_w = min(f['w'] for f in funcs)
    if min_w < 2.0 * h_ref:
        print(f"WARNING: h_ref={h_ref:.4f} too coarse for probe w={min_w:.4f}.")

    x = np.arange(-limit, limit + 1e-9, h_ref)
    xv, yv = np.meshgrid(x, x)
    dA = h_ref**2
    
    norm_factors = np.zeros(len(funcs))
    print(f"Computing Norms for {len(funcs)} functions (Grid {len(x)}x{len(x)})...")
    
    for i, f in enumerate(funcs):
        # We integrate the DENSITY (val^2)
        vals = eval_density(f, xv, yv)
        total_integral = np.sum(vals) * dA
        
        if total_integral < 1e-12:
            norm_factors[i] = 0.0 # Should not happen with Gaussians
        else:
            norm_factors[i] = 1.0 / total_integral
            
    return norm_factors

# ==========================================
# 3. Constrained Optimization (Normalized)
# ==========================================

def optimize_constrained(unique_pts, w_geom, funcs, norm_factors, half, 
                         lambda_reg=1e-3, area_weight=1e4, probe_weight=100.0):
    
    fixed_keys, opt_groups, point_to_gid = group_points(unique_pts, half)
    n_groups = len(opt_groups)
    n_funcs = len(funcs)
    
    print("Building Normalized Matrices...")
    
    A_rows = []
    b_rows = []
    
    # --- 1. Physics Rows ---
    # Equation: Sum(w_i * Norm_j * rho_j(r_i)) = 1.0
    
    # Calculate Fixed Point Contribution (Normalized)
    # Fixed_Contr_j = Norm_j * Sum(w_fixed * rho_j(fixed_pt))
    contrib_fixed_norm = np.zeros(n_funcs)
    for k in fixed_keys:
        pt = unique_pts[k]; w = w_geom[k]
        for i, f in enumerate(funcs):
            val = eval_density(f, pt[0], pt[1])
            contrib_fixed_norm[i] += w * val * norm_factors[i]
            
    # Target is 1.0 for everyone (since we normalized)
    # Modified Target b = 1.0 - Fixed_Contribution
    b_physics = np.ones(n_funcs) - contrib_fixed_norm
    
    # Matrix F (Normalized)
    # F_ij = Norm_i * Sum(rho_i(pt)) for pt in group j
    F = np.zeros((n_funcs, n_groups))
    for j, grp in enumerate(opt_groups):
        pts = [unique_pts[k] for k in grp['indices']]
        for i, f in enumerate(funcs):
            raw_sum = sum(eval_density(f, p[0], p[1]) for p in pts)
            F[i, j] = raw_sum * norm_factors[i]
            
    # Apply Weights (Probes vs Randoms)
    weights = np.ones(n_funcs)
    for i, f in enumerate(funcs):
        if f['type'] == 'probe':
            weights[i] = probe_weight
            
    F_weighted = F * weights[:, None]
    b_weighted = b_physics * weights
    
    A_rows.append(F_weighted)
    b_rows.append(b_weighted)
    
    # --- 2. Regularization Rows ---
    # x ~ x0
    x0 = np.array([np.mean([w_geom[k] for k in g['indices']]) for g in opt_groups])
    reg_mat = np.eye(n_groups) * np.sqrt(lambda_reg)
    reg_rhs = x0 * np.sqrt(lambda_reg)
    
    A_rows.append(reg_mat)
    b_rows.append(reg_rhs)
    
    # --- 3. Area Constraint ---
    # Sum(w) = Area
    counts = np.array([g['count'] for g in opt_groups])
    target_area = np.sum(x0 * counts)
    area_row = counts.reshape(1, -1) * area_weight
    area_rhs = np.array([target_area]) * area_weight
    
    A_rows.append(area_row)
    b_rows.append(area_rhs)
    
    # Stack
    A_final = np.vstack(A_rows)
    b_final = np.hstack(b_rows)
    
    # --- Bounds (Strictly Positive) ---
    lb = x0 * 0.1 # Allow dropping to 10%
    ub = x0 * 3.0 # Allow rising to 300%
    
    print("Solving Constrained Least Squares...")
    res = lsq_linear(A_final, b_final, bounds=(lb, ub), method='trf', verbose=0)
    
    x_opt = res.x
    
    # Distribute
    w_final = w_geom.copy()
    for j, grp in enumerate(opt_groups):
        val = x_opt[j]
        for k in grp['indices']: w_final[k] = val
            
    return w_final, opt_groups, x0, x_opt, F, b_physics, contrib_fixed_norm, point_to_gid

# ==========================================
# 4. Main
# ==========================================

def parse_args():
    p = argparse.ArgumentParser(description="Normalized Density Optimization")
    p.add_argument("--n_sub",    type=int,   default=4,     help="Subdivision count along plug edge.")
    p.add_argument("--half",     type=float, default=1.0,   help="Half-size of square plug.")
    p.add_argument("--margin",   type=int,   default=5,     help="Margin in coarse grid steps beyond plug.")
    p.add_argument("--dmin",     type=float, default=0.2,   help="Smallest segment length ratio near core.")
    p.add_argument("--dmax",     type=float, default=1.0,   help="Largest segment length ratio near boundary.")
    p.add_argument("--n-random", type=int,   default=1000,  help="Number of random orbital basis functions.")
    p.add_argument("--w-min",    type=float, default=1.0,   help="Min width for random orbitals (log-uniform).")
    p.add_argument("--w-max",    type=float, default=3.0,   help="Max width for random orbitals (log-uniform).")
    
    # Default params
    p.add_argument("--lambda-reg",   type=float, default=0.0, help="Regularization strength toward geometric weights.")
    p.add_argument("--area-weight",  type=float, default=0.0, help="Weight for total area preservation constraint.")
    p.add_argument("--probe-weight", type=float, default=10.0, help="Relative weight for probe rows in least squares.")
    p.add_argument("--probe-mults",  type=str,   default="2.5,3.5,4.5", help="Comma-separated multipliers for probe widths (NN distance scaled).")
    
    p.add_argument("--ref-factor", type=int, default=10,  help="Refinement factor for reference grid (h_ref = h_coarse/ref_factor).")
    p.add_argument("--label-fontsize", type=float, default=8.0, help="Font size for point annotations in plots.")
    p.add_argument("--top-n", type=int, default=5, help="How many highest-error functions to visualize.")
    return p.parse_args()

def main():
    args = parse_args()
    
    print("--- 1. Generating Geometry ---")
    onion_map = gen.get_onion_coords(args.n_sub, args.dmin, args.dmax, args.half)
    onion_tris = gen.get_onion_topology(args.n_sub, onion_map)
    coarse_quads = gen.get_coarse_grid(args.n_sub, args.half, args.margin)
    unique_pts, w_total, _, w_coarse, _ = gen.merge_and_analyze(onion_tris, coarse_quads)
    w_geom = w_total
    
    h_coarse = (2.0*args.half)/args.n_sub
    
    # Normalize coarse perimeter weights: outside the plug, use uniform coarse area h^2
    h2_coarse = h_coarse * h_coarse
    for k, p in unique_pts.items():
        if max(abs(p[0]), abs(p[1])) > args.half + 1e-9:  # strictly outside plug box
            w_geom[k] = h2_coarse
    
    print("\n--- 2. Setting up Optimization ---")
    _, temp_groups, _ = group_points(unique_pts, args.half)
    probe_mults = [float(s) for s in args.probe_mults.split(",") if s.strip()!=""]
    funcs = generate_basis_functions(temp_groups, unique_pts, n_random=args.n_random, probe_mults=probe_mults, w_min=args.w_min, w_max=args.w_max)
    
    # Compute Normalization Factors (Integrals of psi^2)
    norm_factors = compute_normalization_factors(funcs, args.half, args.margin, h_coarse, ref_factor=args.ref_factor)
    
    # Coarse-only normalized integrals (original coarse grid, no plug replacement)
    I_coarse_norm = np.zeros(len(funcs))
    for k, w in w_coarse.items():
        if w == 0.0: continue
        pt = unique_pts[k]
        for i, f in enumerate(funcs):
            I_coarse_norm[i] += w * eval_density(f, pt[0], pt[1]) * norm_factors[i]
    
    w_final, groups, x_geom, x_opt, F_norm, b_physics, contrib_fixed, point_to_gid = optimize_constrained(
        unique_pts, w_geom, funcs, norm_factors, args.half, 
        lambda_reg=args.lambda_reg, 
        area_weight=args.area_weight,
        probe_weight=args.probe_weight
    )
    
    # 3. Results Analysis (Unscaled Errors)
    # We report the error in the Normalized integral (Target = 1.0)
    I_geom_norm = contrib_fixed + F_norm @ x_geom
    I_opt_norm  = contrib_fixed + F_norm @ x_opt
    err_coarse  = np.abs(I_coarse_norm - 1.0)
    
    # Since Target is 1.0, the error is effectively relative error
    err_geom = np.abs(I_geom_norm - 1.0)
    err_opt  = np.abs(I_opt_norm - 1.0)
    
    print("\n--- 3. Results (Normalized Density Integrals) ---")
    print(f"Error Reduction (geom->opt): {np.mean(err_geom)/np.mean(err_opt):.2f}x")
    print(f"Mean Error (Geom): {np.mean(err_geom):.6f} (Relative to 1.0)")
    print(f"Mean Error (Opt):  {np.mean(err_opt):.6f} (Relative to 1.0)")
    print(f"Mean Error (Coarse-only): {np.mean(err_coarse):.6f} (Relative to 1.0)")
    
    print("\nPer-function Normalized errors (Top 15):")
    for j, f in enumerate(funcs):
        if j > 15 and f['type'] != 'probe': continue 
        
        if f['type'] == 'probe':
            desc = f"PROBE GID={f['gid']} w={f['w']:.4f}"
        else:
            a = f['ang']
            astr = f"({a[0]:.2f},{a[1]:.2f},{a[2]:.2f})"
            desc = f"{f['rad_type']} w={f['w']:.3f} ang={astr}"
            
        print(f"{j:03d}: {desc:<40} | coarse={err_coarse[j]:.6f} geom={err_geom[j]:.6f} -> opt={err_opt[j]:.6f}")

    print("\n" + "="*85)
    print(f"| {'Radius':<8} | {'Canonical (x, y)':<20} | {'Mult':<4} | {'Old W':<10} | {'New W':<10} | {'Change %':<8} | {'Note'}")
    print("="*85)
    
    for i, g in enumerate(groups):
        k = g['key']
        count = g['count']
        old_w = x_geom[i]
        new_w = x_opt[i]
        chg = (new_w - old_w)/old_w * 100
        
        note = ""
        if g['radius'] < 1e-9: note = "(CORE)"
        elif abs(max(k) - args.half) < 1e-9: note = "(SEAM)"
        elif abs(new_w - old_w*0.1) < 1e-6: note = "LO BOUND"
        elif abs(new_w - old_w*3.0) < 1e-6: note = "HI BOUND"
        
        print(f"| {g['radius']:<8.4f} | {fmt_key(k):<20} | {count:<4} | {old_w:<10.6f} | {new_w:<10.6f} | {chg:>+7.2f}% | {note}")
    print("="*85)
    
    # 4. Plot
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    xs, ys = [], []
    for k in unique_pts:
        p = unique_pts[k]
        xs.append(p[0]); ys.append(p[1])
    
    # background reference grid (high-res) for context
    limit = args.half + (args.margin + 0.5) * h_coarse
    h_ref = h_coarse / args.ref_factor
    ref_x = np.arange(-limit, limit + 1e-9, h_ref)
    ref_y = ref_x
    ref_xv, ref_yv = np.meshgrid(ref_x, ref_y)
    ref_flat_x = ref_xv.ravel()
    ref_flat_y = ref_yv.ravel()
        
    ax = axes[0]
    ax.scatter(ref_flat_x, ref_flat_y, c='lightgray', s=5, alpha=0.2, edgecolor='none', zorder=0)
    diffs = [w_final[k] - w_geom[k] for k in unique_pts]
    sc0 = ax.scatter(xs, ys, c=diffs, cmap='coolwarm', s=30, edgecolor='none', zorder=1)
    plt.colorbar(sc0, ax=ax, label="Correction")
    for k in unique_pts:
        gid = point_to_gid.get(k, -1)
        if gid >= 0:
            ax.annotate(f"{gid}", (unique_pts[k][0], unique_pts[k][1]), 
                        fontsize=args.label_fontsize, ha="center", va="center", alpha=0.5)
    ax.set_aspect('equal'); ax.set_title("Corrections & GIDs")
    
    ax = axes[1]
    ax.scatter(ref_flat_x, ref_flat_y, s=0.5, alpha=0.2, color='k', edgecolor='none', zorder=0)
    vals = [w_final[k] for k in unique_pts]
    sc1 = ax.scatter(xs, ys, c=vals, cmap='viridis', s=30, edgecolor='none', zorder=1)
    plt.colorbar(sc1, ax=ax, label="Weight")
    for k in unique_pts:
        ax.annotate(f"{w_final[k]:.3g}", (unique_pts[k][0], unique_pts[k][1]), 
                    fontsize=args.label_fontsize, ha="center", va="center")
    ax.set_aspect('equal'); ax.set_title("Final Weights (Positive)")
    
    ax = axes[2]
    # Filter out fixed points for clearer histogram? No, standard error plot
    ax.plot(err_coarse, label='Coarse', marker='.', ls='none', alpha=0.5)
    ax.plot(err_geom, label='Geom', marker='o', ls='none', alpha=0.5)
    ax.plot(err_opt, label='Opt', marker='x', ls='none', alpha=0.8)
    ax.set_yscale('log')
    ax.set_title("Normalized Errors (Target=1.0)")
    ax.set_ylabel("Abs Error (Relative to Area)")
    ax.legend()
    
    for ax in axes[:2]:
        rect = plt.Rectangle((-args.half, -args.half), 2*args.half, 2*args.half, fill=False, ec='black', ls='--')
        ax.add_patch(rect)
    
    # --- Plot Top-N worst errors (heatmap of density with sample dots) ---
    top_n = max(1, args.top_n)
    worst_idx = np.argsort(err_opt)[::-1][:top_n]
    n_cols = min(3, top_n)
    n_rows = int(np.ceil(top_n / n_cols))
    fig_top, axes_top = plt.subplots(n_rows, n_cols, figsize=(6*n_cols, 5*n_rows))
    axes_top = np.array(axes_top).reshape(n_rows, n_cols)
    
    extent = (-limit, limit, -limit, limit)
    for ax_idx, j in enumerate(worst_idx):
        r = ax_idx // n_cols
        c = ax_idx % n_cols
        axh = axes_top[r, c]
        vals = eval_density(funcs[j], ref_xv, ref_yv)
        im = axh.imshow(vals, extent=extent, origin='lower', cmap='magma')
        fig_top.colorbar(im, ax=axh, fraction=0.046, pad=0.04, label='Density')
        axh.scatter(xs, ys, s=8, c='white', alpha=0.7, edgecolor='none')
        f = funcs[j]
        if f['type'] == 'probe':
            desc = f"PROBE gid={f['gid']} w={f['w']:.3f}"
        else:
            a = f['ang']; desc = f"{f['rad_type']} w={f['w']:.3f} ang=({a[0]:.2f},{a[1]:.2f},{a[2]:.2f})"
        axh.set_title(f"#{j} err {err_opt[j]:.3e} ({desc})")
        axh.set_xlim(-limit, limit); axh.set_ylim(-limit, limit); axh.set_aspect('equal')
    # Hide unused axes if any
    for k in range(len(worst_idx), n_rows*n_cols):
        r = k // n_cols; c = k % n_cols
        axes_top[r, c].axis('off')
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()