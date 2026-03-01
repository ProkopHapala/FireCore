import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

def match_roots(prev_roots, curr_roots):
    """
    Nearest-neighbor tracker to match polynomial roots across U-steps.
    This prevents the N roots from shuffling, keeping streamlines perfectly continuous!
    """
    ordered_roots = []
    curr_pool = list(curr_roots)
    for p in prev_roots:
        # Find the root in current pool closest to the previous root
        dists = [np.abs(c - p) for c in curr_pool]
        best_idx = np.argmin(dists)
        ordered_roots.append(curr_pool.pop(best_idx))
    return ordered_roots

def main():
    parser = argparse.ArgumentParser(description="Analytical Conformal Grid: 1D Cuts")
    parser.add_argument('--u_density',        type=float, default=2.5)
    parser.add_argument('--v_density',        type=int,   default=16)
    parser.add_argument('--filter_mode',      type=str,   default='isoline', choices=['isoline','geometric','both'])
    parser.add_argument('--highlight_u_idx',  type=int,   default=4,   help="Index of the Isoline (Shell) to extract")
    parser.add_argument('--highlight_v_idx',  type=int,   default=0,   help="Index of the Streamline (Ray) to extract")
    parser.add_argument('--highlight_branch', type=int,   default=0,   help="Which of the N roots (branches) to track for the streamline")
    parser.add_argument('--u_boundary',       type=float, default=1.5, help='U_boundary multiplier. Default=2.5')

    try:
        args = parser.parse_args()
    except SystemExit:
        args = parser.parse_args(args=[])

    # 1. Define Molecule
    atoms = [
    #     x     y     r    multiplicity
        ( 0.0,  0.0, 0.40, 2       ),  # small
        (-0.5, -0.0, 0.20, 1       ),  # medium
        ( 0.0, -0.5, 0.20, 1       )   # large
    ]
    N_degree = sum([a[3] for a in atoms]) # Total polynomial degree (N=6)

    # Center of Mass (used for sorting the isoline loop later)
    com_x = np.mean([a[0] for a in atoms])
    com_y = np.mean([a[1] for a in atoms])

    base_poly = [1.0]
    for (cx, cy, r, weight) in atoms:
        weight_int = int(round(weight))
        assert abs(weight - weight_int) < 1e-8, f"Atom weight must be integer-like for multiplicity, got {weight}"
        center = cx + 1j * cy
        for _ in range(weight_int):
            base_poly = np.convolve(base_poly, [1.0, -center])

    # 2. Calculate Molecular Surface
    u_surface_vals = []
    for (cx, cy, r, _) in atoms:
        for angle in [0, np.pi/2, np.pi, 3*np.pi/2]:
            test_z = (cx + r * np.cos(angle)) + 1j * (cy + r * np.sin(angle))
            val = np.log(np.abs(np.polyval(base_poly, test_z)))
            u_surface_vals.append(val)
    U_boundary = max(u_surface_vals)*args.u_boundary

    # 3. Setup Grid Targets
    step_U = 1.0 / args.u_density
    freq_U = args.u_density
    freq_V = args.v_density / (2 * np.pi)

    # Checkerboard background (visualizes mapping continuity)
    res = 1000
    extent = [-2.0, 2.0, -2.0, 2.0]
    x_lin = np.linspace(extent[0], extent[1], res)
    y_lin = np.linspace(extent[2], extent[3], res)
    X, Y = np.meshgrid(x_lin, y_lin)
    Z_grid = X + 1j * Y
    W = np.polyval(base_poly, Z_grid)
    U = np.log(np.abs(W))
    V = np.angle(W)
    checker = (np.floor(U * freq_U) + np.floor(V * freq_V)) % 2

    if args.filter_mode in ['isoline', 'both']:
        start_idx = int(np.ceil(U_boundary / step_U))
        U_targets = np.arange(start_idx * step_U, 6.0, step_U)
    else:
        U_targets = np.arange(-8.0, 6.0, step_U)

    V_targets = np.linspace(-np.pi, np.pi, args.v_density, endpoint=False)

    # ==========================================
    # 4. STRUCTURED DATA COLLECTION (The Database)
    # ==========================================
    sample_points = []
    point_id_counter = 0

    # We loop V first, then U. This lets us trace a streamline outwards from the atoms!
    for v_idx, v in enumerate(V_targets):
        prev_roots = None
        
        for u_idx, u in enumerate(U_targets):
            w_target = np.exp(u + 1j * v)
            poly = np.copy(base_poly)
            poly[-1] -= w_target
            roots = np.roots(poly)
            
            # Root Tracking to maintain streamline continuity
            if prev_roots is None:
                ordered_roots = roots # Initial arbitrary order at the surface
            else:
                ordered_roots = match_roots(prev_roots, roots)
            
            prev_roots = ordered_roots
            
            # Store structured data points
            for branch_idx, root in enumerate(ordered_roots):
                rx, ry = np.real(root), np.imag(root)
                
                point_data = {
                    'id': point_id_counter,
                    'x': rx, 'y': ry,
                    'u_idx': u_idx, 'v_idx': v_idx,
                    'u_val': u, 'v_val': v,
                    'branch': branch_idx
                }
                sample_points.append(point_data)
                point_id_counter += 1

    # ==========================================
    # 5. EXTRACT 1D CUTS (Queries)
    # ==========================================
    
    # Query A: Extract a single continuous ISOLINE (Shell)
    isoline_pts = [p for p in sample_points if p['u_idx'] == args.highlight_u_idx]
    # To draw a line through them, we sort them circularly around the Center of Mass
    isoline_pts.sort(key=lambda p: np.arctan2(p['y'] - com_y, p['x'] - com_x))
    
    # Query B: Extract a single continuous STREAMLINE (Ray)
    # Note: We must specify v_idx AND the branch, otherwise we get N radiating rays!
    streamline_pts = [p for p in sample_points if p['v_idx'] == args.highlight_v_idx and p['branch'] == args.highlight_branch]
    # Already sorted outward automatically because we looped U in increasing order!

    # Console Output
    print(f"--- 1D CUT EXTRACTION ---")
    print(f"Total points generated: {len(sample_points)}")
    
    print(f"\nISOLINE (u_idx = {args.highlight_u_idx}):")
    iso_ids = [p['id'] for p in isoline_pts]
    print(f"Point IDs: {iso_ids}")
    
    print(f"\nSTREAMLINE (v_idx = {args.highlight_v_idx}, branch = {args.highlight_branch}):")
    stream_ids = [p['id'] for p in streamline_pts]
    print(f"Point IDs: {stream_ids}")

    # ==========================================
    # 6. PLOTTING
    # ==========================================
    fig, ax = plt.subplots(figsize=(10, 10))

    # Checkerboard background
    cmap = ListedColormap(['#ffffff', '#e2e8f0'])
    ax.imshow(checker, extent=extent, origin='lower', cmap=cmap, alpha=1.0, zorder=1)
    
    # Background points (Dim Red)
    all_x = [p['x'] for p in sample_points]
    all_y = [p['y'] for p in sample_points]
    ax.plot(all_x, all_y, linestyle='none', marker='.', color='#e63946', alpha=0.3, markersize=5, zorder=2)

    # Highlight 1: The Isoline (Neon Green)
    if isoline_pts:
        iso_x = [p['x'] for p in isoline_pts] + [isoline_pts[0]['x']] # append first point to close the loop
        iso_y = [p['y'] for p in isoline_pts] + [isoline_pts[0]['y']]
        ax.plot(iso_x, iso_y, color='#00ff00', linewidth=2.5, marker='o', markersize=6, markerfacecolor='w', zorder=4, label=f"1D Isoline Cut (u_idx={args.highlight_u_idx})")

    # Highlight 2: The Streamline (Bright Cyan)
    if streamline_pts:
        str_x = [p['x'] for p in streamline_pts]
        str_y = [p['y'] for p in streamline_pts]
        ax.plot(str_x, str_y, color='#00e5ff', linewidth=2.5, marker='s', markersize=6, markerfacecolor='w', zorder=5, label=f"1D Streamline Cut (v_idx={args.highlight_v_idx}, br={args.highlight_branch})")

    # Draw Atom Boundaries
    for (cx, cy, r, _) in atoms:
        circle = plt.Circle((cx, cy), r, color='#1e293b', fill=False, linestyle='--', linewidth=1.5, zorder=3)
        ax.add_patch(circle)
        ax.plot(cx, cy, 'k+', markersize=8)

    ax.set_aspect('equal')
    ax.set_xlim([-2.0, 2.0])
    ax.set_ylim([-2.0, 2.0])
    ax.set_title("1D Manifold Extraction from Conformal Grid\nContinuous Isoline & Streamline Tracking", fontsize=14)
    ax.legend(loc='upper right')
    
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    main()