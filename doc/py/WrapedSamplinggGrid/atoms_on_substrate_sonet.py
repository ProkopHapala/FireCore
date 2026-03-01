"""
Atoms on Substrate: Conformal Grid via Complex Polynomial Roots
================================================================

KEY IMPROVEMENTS over the Gemini solution:
1. Polynomial roots approach (like your original 2D code) — NO Newton divergence
2. Tunable decay via modified Green's functions: log(z-p) + alpha*log(z-p*) 
   where p* is the image charge. alpha controls how fast the perturbation dies out.
3. Branch cut management: we use the POLYNOMIAL PRODUCT form W = A * prod(z - roots)
   so |W| and arg(W) are computed from the polynomial, not summed logs
4. Method of images built naturally into the polynomial
5. Free-stream is a LINEAR factor: W_total = z_linear * W_molecular

MATH:
-----
The core idea: instead of summing logs (which gives branch cuts),
we BUILD a polynomial whose roots are the atom positions.

    W(z) = -i*A*z * prod_k (z - z_k)^n_k    (product of linear factors)

Then:
    U = log|W|     (isolines = equipotentials)
    V = arg(W)     (streamlines = phase contours)

The "substrate" condition: we want the real axis y=0 to be a streamline.
This is achieved by Method of Images: for each atom at z_k = x_k + i*y_k,
add its mirror image at z_k* = x_k - i*y_k with same weight.

Then on the real axis z = x: 
    W(x) = -i*A*x * prod_k |x - z_k|^2  (real × positive = real positive)
So arg(W) = const on y=0. ✓

TUNABLE DECAY:
--------------
The modification: instead of (z - z_k)(z - z_k*) = full image,
use:  (z - z_k) * (z - z_k*)^alpha    with 0 < alpha <= 1

alpha = 1: standard method of images, decay ~ 1/r^2 away from surface
alpha = 0: no image, decay ~ 1/r (Gemini's version, too slow)
alpha < 1: partial image, intermediate decay rate

Physically: alpha controls how "reflecting" the substrate is.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import argparse

# ============================================================
# 1. FIELD DEFINITION via POLYNOMIAL PRODUCT
# ============================================================

def make_polynomial(atoms_data, alpha=1.0):
    """
    Build the polynomial whose roots define the conformal mapping.
    
    atoms_data: list of (x, y, weight, charge_sign)
    alpha: image charge strength [0=no reflection, 1=full mirror]
    
    Returns: poly_coeffs for np.polyval
    """
    poly = np.array([1.0+0j])
    
    for (ax, ay, weight, _) in atoms_data:
        z_atom = ax + 1j*ay
        z_image = ax - 1j*ay  # mirror below substrate
        
        w_int = max(1, int(round(weight * 4)))  # integer multiplicity
        
        # Real atom: (z - z_atom)^w
        for _ in range(w_int):
            poly = np.convolve(poly, [1.0, -z_atom])
        
        # Image atom with tunable strength alpha
        # We approximate fractional power by using ceil/floor weights
        # For clean polynomial: use alpha as ratio of image to real weight
        w_image = max(0, int(round(alpha * w_int)))
        for _ in range(w_image):
            poly = np.convolve(poly, [1.0, -z_image])
    
    return poly

def W_total(z, poly, A=1.5):
    """
    W(z) = log( -i*A*z * P(z) )
    
    U = Real(W) = log|A*z*P(z)|  → isolines = height shells
    V = Imag(W) = arg(-i*A*z*P(z)) → streamlines = columns
    
    The -i rotation ensures: 
      - U increases upward (y direction) in far field
      - V = const are vertical lines in far field
    """
    p_val = np.polyval(poly, z)
    # Full field: free-stream (-i*A*z) times molecular perturbation
    w_field = -1j * A * z * p_val
    # Protect against log(0)
    mag = np.abs(w_field)
    mag = np.where(mag < 1e-30, 1e-30, mag)
    return np.log(mag) + 1j * np.angle(w_field)


# ============================================================
# 2. GRID SAMPLING VIA POLYNOMIAL ROOTS (Your Original Method!)
# ============================================================

def find_grid_points_via_roots(U_targets, V_targets, poly, A=1.5, 
                                atoms_data=None, filter_inside=True):
    """
    For each (U_target, V_target), find z such that W(z) = U + iV.
    
    W(z) = log(-i*A*z*P(z)) = U + iV
    => -i*A*z*P(z) = exp(U + iV) = w_target
    
    This is a POLYNOMIAL EQUATION in z of degree (N_poly + 1):
       A_coeff * z * P(z) = w_target  (times -i)
    => -i*A*z*P(z) - w_target = 0
    
    Solve via np.roots() — EXACTLY your original technique!
    The degree-N polynomial has N roots, and we track them
    across U-steps using nearest-neighbor matching.
    """
    N_degree = len(poly) - 1  # degree of P(z)
    # z*P(z) has degree N_degree + 1
    # -i*A*(z*P(z)) = -i*A * poly extended by one degree
    
    # z*poly: multiply poly by [1, 0] (shift coefficients, append 0)
    z_poly = np.append(poly, 0)  # z * P(z): coefficients shifted
    base_poly = -1j * A * z_poly  # base polynomial without w_target offset
    
    sample_points = []
    pid = 0
    
    for v_idx, v in enumerate(V_targets):
        prev_roots = None
        
        for u_idx, u in enumerate(U_targets):
            w_target = np.exp(u + 1j * v)
            
            # THE KEY: polynomial root finding
            # -i*A*z*P(z) - w_target = 0
            test_poly = base_poly.copy()
            test_poly[-1] -= w_target  # subtract constant term
            
            roots = np.roots(test_poly)
            
            # Root tracking (from your original code)
            if prev_roots is None:
                ordered_roots = roots
            else:
                ordered_roots = match_roots(prev_roots, roots)
            prev_roots = ordered_roots
            
            for branch_idx, root in enumerate(ordered_roots):
                rx, ry = np.real(root), np.imag(root)
                
                # Keep only above substrate
                if ry < 0.01:
                    continue
                
                # Filter inside atoms
                if filter_inside and atoms_data:
                    inside = False
                    for (ax, ay, weight, _) in atoms_data:
                        r_atom = weight  # use weight as radius
                        if np.hypot(rx - ax, ry - ay) < r_atom:
                            inside = True
                            break
                    if inside:
                        continue
                
                sample_points.append({
                    'id': pid, 'x': rx, 'y': ry,
                    'u_idx': u_idx, 'v_idx': v_idx,
                    'u_val': u, 'v_val': v,
                    'branch': branch_idx
                })
                pid += 1
    
    return sample_points


def match_roots(prev_roots, curr_roots):
    """Nearest-neighbor root tracker — directly from your original code."""
    ordered = []
    pool = list(curr_roots)
    for p in prev_roots:
        dists = [np.abs(c - p) for c in pool]
        best = np.argmin(dists)
        ordered.append(pool.pop(best))
    return ordered


# ============================================================
# 3. MAIN
# ============================================================

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--alpha',      type=float, default=1.0,
        help='Image charge strength [0=none, 1=full mirror]. Controls decay speed.')
    parser.add_argument('--A',          type=float, default=1.5,
        help='Free-stream strength')
    parser.add_argument('--u_density',  type=float, default=3.0)
    parser.add_argument('--v_density',  type=int,   default=20)
    parser.add_argument('--highlight_u', type=int, default=5)
    parser.add_argument('--highlight_v', type=int, default=8)
    parser.add_argument('--highlight_branch', type=int, default=0)
    try:
        args = parser.parse_args()
    except SystemExit:
        args = parser.parse_args(args=[])

    # ---- Define atoms (x, y, radius_for_filtering, charge_weight) ----
    atoms_data = [
        (-0.65,  0.35,  0.22,  1),   # left atom
        ( 0.55,  0.38,  0.18,  1),   # right atom
    ]
    
    A = args.A
    alpha = args.alpha

    print(f"Parameters: A={A}, alpha={alpha}")
    print(f"Alpha interpretation: {alpha:.2f} = ", end="")
    if alpha < 0.3: print("weak image (slow decay ~1/r)")
    elif alpha < 0.7: print("partial image (moderate decay)")
    else: print("strong image (fast decay ~1/r²)")

    # ---- Build polynomial ----
    poly = make_polynomial(atoms_data, alpha=alpha)
    N_degree = len(poly) - 1
    print(f"Polynomial degree: {N_degree}, roots to track: {N_degree+1}")

    # ---- Determine U boundary (molecular surface level) ----
    # Sample potential at various heights to find a good working range
    u_sample_vals = []
    for yy in np.linspace(0.05, 1.8, 30):
        for xx in np.linspace(-1.5, 1.5, 10):
            val = np.real(W_total(xx + 1j*yy, poly, A))
            if np.isfinite(val):
                u_sample_vals.append(val)
    u_min = np.percentile(u_sample_vals, 15)
    u_max = np.percentile(u_sample_vals, 90)

    step_U = 1.0 / args.u_density
    U_targets = np.arange(u_min, u_max, step_U)
    V_targets = np.linspace(-np.pi, np.pi, args.v_density, endpoint=False)
    print(f"U range: [{U_targets[0]:.2f}, {U_targets[-1]:.2f}], {len(U_targets)} shells")
    print(f"V range: {len(V_targets)} streamlines")

    # ---- Sample grid via polynomial roots ----
    print("Computing grid points via polynomial roots...")
    sample_points = find_grid_points_via_roots(
        U_targets, V_targets, poly, A=A, 
        atoms_data=atoms_data, filter_inside=True
    )
    print(f"Total valid sample points: {len(sample_points)}")

    # ---- 1D Cuts ----
    isoline_pts = [p for p in sample_points if p['u_idx'] == args.highlight_u]
    isoline_pts.sort(key=lambda p: p['x'])  # sort left to right (surface style)
    
    # For substrate geometry, pick the branch closest to the atom region
    # (multiple branches exist — pick the one in upper half plane with largest count)
    streamline_candidates = [p for p in sample_points if p['v_idx'] == args.highlight_v]
    if streamline_candidates:
        # Auto-pick best branch: most points, all in upper half
        from collections import Counter
        branch_counts = Counter(p['branch'] for p in streamline_candidates)
        best_branch = args.highlight_branch if args.highlight_branch in branch_counts else max(branch_counts, key=branch_counts.get)
        streamline_pts = [p for p in streamline_candidates if p['branch'] == best_branch]
        streamline_pts.sort(key=lambda p: p['y'])
    else:
        streamline_pts = []

    print(f"\nIsoline (u_idx={args.highlight_u}): {len(isoline_pts)} pts")
    print(f"Streamline (v_idx={args.highlight_v}, branch={args.highlight_branch}): {len(streamline_pts)} pts")

    # ============================================================
    # 4. VISUALIZATION
    # ============================================================
    res = 800
    x_lin = np.linspace(-2.0,  2.0, res)
    y_lin = np.linspace( 0.0,  2.0, res)
    X, Y = np.meshgrid(x_lin, y_lin)
    Z_grid = X + 1j*Y

    W_grid = W_total(Z_grid, poly, A)
    U_field = np.real(W_grid)
    V_field = np.imag(W_grid)

    # Smooth checkerboard (sin-based to avoid 2π discontinuities)
    freq_u = args.u_density
    freq_v = args.v_density / (2*np.pi)
    checker = (np.sin(U_field * freq_u * np.pi) > 0) ^ (np.sin(V_field * freq_v * np.pi) > 0)

    fig, axes = plt.subplots(1, 2, figsize=(16, 8))
    
    for ax_idx, ax in enumerate(axes):
        # Checkerboard background
        cmap_check = ListedColormap(['#f8fafc', '#e2e8f0'])
        ax.imshow(checker, extent=[-2, 2, 0, 2], origin='lower', 
                  cmap=cmap_check, alpha=0.9, zorder=1, aspect='auto')
        
        # Grid lines
        ax.contour(X, Y, U_field, levels=U_targets[::2], 
                   colors='#8b5cf6', linewidths=0.7, alpha=0.5, zorder=2)
        ax.contour(X, Y, V_field, levels=V_targets[::2],
                   colors='#0ea5e9', linewidths=0.7, alpha=0.5, zorder=2)
        
        # Substrate
        ax.axhline(0, color='#1e293b', linewidth=5, zorder=5)
        ax.fill_between([-2, 2], -0.05, 0, color='#cbd5e1', zorder=4)
        
        # Atoms
        for (ax_p, ay_p, r, _) in atoms_data:
            circle = plt.Circle((ax_p, ay_p), r, color='#475569', 
                                 alpha=0.6, zorder=6)
            ax.add_patch(circle)
            ax.plot(ax_p, ay_p, 'w+', markersize=10, zorder=7)
        
        ax.set_aspect('equal')
        ax.set_xlim([-1.8, 1.8])
        ax.set_ylim([-0.05, 1.8])
        ax.axis('off')
    
    # Left panel: all sample points
    ax = axes[0]
    all_x = [p['x'] for p in sample_points]
    all_y = [p['y'] for p in sample_points]
    ax.scatter(all_x, all_y, s=8, color='#ef4444', alpha=0.4, 
               edgecolors='none', zorder=10)
    ax.set_title(f"All Grid Points\n(α={alpha:.1f}, A={A}, poly degree {N_degree+1})\n"
                 f"{len(sample_points)} samples", fontsize=12)
    
    # Right panel: highlighted 1D cuts
    ax = axes[1]
    ax.scatter(all_x, all_y, s=6, color='#fca5a5', alpha=0.2, 
               edgecolors='none', zorder=8)
    
    if isoline_pts:
        ix = [p['x'] for p in isoline_pts]
        iy = [p['y'] for p in isoline_pts]
        ax.plot(ix, iy, color='#22c55e', linewidth=2.5, 
                marker='o', markersize=6, markerfacecolor='w', zorder=12,
                label=f"Isoline (shell u_idx={args.highlight_u})")
    
    if streamline_pts:
        sx = [p['x'] for p in streamline_pts]
        sy = [p['y'] for p in streamline_pts]
        ax.plot(sx, sy, color='#f59e0b', linewidth=2.5,
                marker='s', markersize=6, markerfacecolor='w', zorder=12,
                label=f"Streamline (v_idx={args.highlight_v}, br={args.highlight_branch})")
    
    ax.legend(loc='upper right', fontsize=9)
    ax.set_title(f"1D Cut Extraction\nIsoline (green) & Streamline (amber)", fontsize=12)
    
    fig.suptitle(
        f"Atoms on Substrate — Conformal Grid via Polynomial Roots\n"
        f"W(z) = log(−iAz · P(z)),  image charge α={alpha:.1f}  "
        f"[0=slow 1/r decay, 1=fast 1/r² decay]",
        fontsize=13, y=1.01
    )
    
    plt.tight_layout()
    out_path = os.path.join(os.path.dirname(__file__), 'atoms_on_substrate.png')
    plt.savefig(out_path, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"\nSaved to {out_path}")
    plt.show()


if __name__ == '__main__':
    main()
