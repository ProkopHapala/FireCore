#!/usr/bin/env python3
"""
Corrected 2D plug visualizer.

Fixes:
1. Symmetry: Uses a structured quad-grid topology that collapses to the core.
   This ensures the number of points is constant per layer, preventing "sawtooth" asymmetry.
2. Weighting: Splits every resulting quad into 4 triangles using a centroid.
   This ensures weights are distributed symmetrically and smoothly.
"""

import numpy as np
import matplotlib.pyplot as plt
import math
import argparse
import sys


def mapped_t_values(n_layers, dmin=1.0, dmax=1.0, weight_func=None):
    """
    Return t-levels (size n_layers+1) for radial spacing.
    """
    if n_layers == 0:
        return np.array([0.0, 1.0])
    u = np.linspace(0.0, 1.0, n_layers)
    if weight_func is None:
        # u=0 -> boundary (dmax), u=1 -> core (dmin)
        weights = dmax - (dmax - dmin) * u
    else:
        weights = weight_func(u)
    weights = np.maximum(weights, 1e-12)
    weights = weights / weights.sum()
    levels = np.concatenate([[0.0], np.cumsum(weights)])
    levels[-1] = 1.0
    return levels


def triangle_subdivision(v0, v1, core, n_sub, dmin=0.2, dmax=1.0, weight_func=None, decimals=12):
    """
    Generate a symmetric mesh for the triangle (v0, v1, core).
    
    Strategy: 
    - Create layers with constant point count (n_sub + 1).
    - Layer 0 is the edge v0-v1.
    - Layer n_sub is collapsed to 'core'.
    - This forms a grid of Quads.
    - Each Quad is split into 4 triangles using its centroid to preserve symmetry.
    """
    t_levels = mapped_t_values(n_sub, dmin=dmin, dmax=dmax, weight_func=weight_func)
    points = {}
    triangles = []

    # 1. Generate Points
    # grid[l][m] where l is radial (0..n_sub), m is lateral (0..n_sub)
    for l in range(n_sub + 1):
        t = t_levels[l]
        # Endpoints of the chord at this radial level
        row_start = (1 - t) * v0 + t * core
        row_end   = (1 - t) * v1 + t * core
        
        # If at core (l == n_sub), collapse all points to core
        if l == n_sub:
            for m in range(n_sub + 1):
                points[(l, m)] = core
        else:
            # Interpolate uniformly along the chord
            for m in range(n_sub + 1):
                s = m / n_sub
                pt = (1 - s) * row_start + s * row_end
                points[(l, m)] = pt

    # 2. Generate Segments and Triangles (Splitting Quads)
    segments = []
    
    # We iterate through the "cells" of the grid
    for l in range(n_sub):
        for m in range(n_sub):
            # The four corners of the quad/cell
            p00 = points[(l, m)]       # Top-Left
            p10 = points[(l, m+1)]     # Top-Right
            p11 = points[(l+1, m+1)]   # Bot-Right
            p01 = points[(l+1, m)]     # Bot-Left
            
            # Check for degeneracy (at the core, p01 and p11 are the same point)
            # However, for robustness, we just treat everything as 4 triangles 
            # or 1 triangle if the bottom edge is collapsed.
            
            # Compute Centroid
            center = (p00 + p10 + p11 + p01) / 4.0
            
            # Store centroid point (keyed by tuple to avoid conflict with grid points)
            # We use a float key or distinct tuple structure
            c_key = (l, m, 'c') 
            
            # Note: We don't need to export the centroid to the 'points' dict for 
            # external usage unless we want it annotated. 
            # We'll just use the coordinate for triangles.
            
            # Create 4 triangles: (p00, p10, c), (p10, p11, c), (p11, p01, c), (p01, p00, c)
            # This ensures perfect symmetry.
            
            tris_local = [
                (p00, p10, center),
                (p10, p11, center),
                (p11, p01, center),
                (p01, p00, center)
            ]
            
            # Filter out degenerate triangles (zero area) which happen at the core
            # At the core, p11 == p01 == core.
            # Triangle (p11, p01, c) -> (core, core, c) -> Area 0.
            valid_tris = []
            for a, b, c_pt in tris_local:
                # Simple check: if vertices are numerically identical, skip
                d1 = np.linalg.norm(a - b)
                d2 = np.linalg.norm(b - c_pt)
                d3 = np.linalg.norm(c_pt - a)
                if d1 > 1e-12 and d2 > 1e-12 and d3 > 1e-12:
                    triangles.append((a, b, c_pt))
                    # Add internal segments for visualization
                    segments.append((a, c_pt))
                    segments.append((b, c_pt))

            # Add boundary segments of the quad (only once)
            # Horizontal (l)
            segments.append((p00, p10))
            # Vertical (m)
            segments.append((p00, p01))
            # The other two are handled by neighbors, except at boundaries
            if m == n_sub - 1: segments.append((p10, p11))
            if l == n_sub - 1: segments.append((p01, p11))

    # Deduplicate triangles
    tri_keys = set()
    unique_triangles = []
    for a, b, c in triangles:
        ka = tuple(np.round(a, decimals=decimals))
        kb = tuple(np.round(b, decimals=decimals))
        kc = tuple(np.round(c, decimals=decimals))
        key = tuple(sorted([ka, kb, kc]))
        if key in tri_keys:
            continue
        tri_keys.add(key)
        unique_triangles.append((a, b, c))

    # Clean points dict (remove duplicates if core was added multiple times)
    # The original logic expects points to be returned.
    return points, segments, unique_triangles


def build_square_plug(center=np.array([0.0, 0.0]), half=1.0):
    cx, cy = center
    verts = [
        np.array([cx - half, cy - half]),
        np.array([cx + half, cy - half]),
        np.array([cx + half, cy + half]),
        np.array([cx - half, cy + half]),
    ]
    return verts, center


def plot_segments(ax, segments, color="C0", lw=1.0, alpha=0.8):
    for a, b in segments:
        xs = [a[0], b[0]]
        ys = [a[1], b[1]]
        ax.plot(xs, ys, color=color, lw=lw, alpha=alpha)


def dedupe_segments(raw_segments, unique_pts, decimals=12):
    seg_keys = set()
    unique_segs = []
    for a, b in raw_segments:
        ka = tuple(np.round(a, decimals=decimals))
        kb = tuple(np.round(b, decimals=decimals))
        # Ensure we only keep segments between known unique points
        if ka in unique_pts and kb in unique_pts:
            key = tuple(sorted([ka, kb]))
            if key in seg_keys:
                continue
            seg_keys.add(key)
            unique_segs.append((unique_pts[ka], unique_pts[kb]))
    return unique_segs


def triangle_weights(triangles, decimals=12):
    """Accumulate area/3 per vertex for each triangle."""
    weights = {}
    for a, b, c in triangles:
        # Cross product for 2D area
        area = 0.5 * abs((b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0]))
        contrib = area / 3.0
        for p in (a, b, c):
            k = tuple(np.round(p, decimals=decimals))
            weights[k] = weights.get(k, 0.0) + contrib
    return weights


def generate_plug(verts, core, n_sub, dmin, dmax, weight_func=None, decimals=12):
    unique_pts = {}
    all_segments = []
    all_triangles = []
    for i in range(4):
        v0 = verts[i]
        v1 = verts[(i + 1) % 4]
        pts, segs, tris = triangle_subdivision(v0, v1, core, n_sub=n_sub, dmin=dmin, dmax=dmax, weight_func=weight_func, decimals=decimals)
        
        # Flatten points
        for p in pts.values():
            key = tuple(np.round(p, decimals=decimals))
            if key not in unique_pts:
                unique_pts[key] = p
        
        all_segments.extend(segs)
        all_triangles.extend(tris)

    unique_segs = dedupe_segments(all_segments, unique_pts, decimals=decimals)
    return unique_pts, unique_segs, all_triangles


def generate_regular_grid(half, n_sub, margin, remove_center=False, decimals=12):
    h = 2 * half / n_sub
    n_side = n_sub + 2 * margin + 1
    coords = np.linspace(-half - margin * h, half + margin * h, n_side)
    pts = {}
    
    # Generate points
    for x in coords:
        for y in coords:
            # remove_center logic: strictly inside the square
            if remove_center and (abs(x) < half - 1e-9 and abs(y) < half - 1e-9):
                continue
            key = tuple(np.round([x, y], decimals=decimals))
            pts[key] = np.array([x, y])

    raw_segments = []
    weights = {}
    
    # Generate Cells and Weights
    for ix, x in enumerate(coords):
        for iy, y in enumerate(coords):
            # Skip if cell is strictly inside hole (approx check based on center)
            if ix + 1 >= len(coords) or iy + 1 >= len(coords):
                continue
                
            x2 = coords[ix + 1]
            y2 = coords[iy + 1]
            cx, cy = 0.5*(x+x2), 0.5*(y+y2)
            
            # If we are removing center, skip cells that are effectively inside
            # A cell is inside if its center is inside
            if remove_center and (abs(cx) < half and abs(cy) < half):
                continue

            keys_cell = [
                tuple(np.round([x, y], decimals=decimals)),
                tuple(np.round([x2, y], decimals=decimals)),
                tuple(np.round([x, y2], decimals=decimals)),
                tuple(np.round([x2, y2], decimals=decimals)),
            ]
            
            # Only process if all corners exist in our point set
            if all(k in pts for k in keys_cell):
                # Add segments
                p1, p2, p3, p4 = [pts[k] for k in keys_cell] # BL, BR, TL, TR order depends on grid
                # Just add edges
                raw_segments.append((pts[keys_cell[0]], pts[keys_cell[1]])) # Bottom
                raw_segments.append((pts[keys_cell[0]], pts[keys_cell[2]])) # Left
                raw_segments.append((pts[keys_cell[1]], pts[keys_cell[3]])) # Right
                raw_segments.append((pts[keys_cell[2]], pts[keys_cell[3]])) # Top
                
                # Add Weight (Area/4)
                area_quarter = (h * h) * 0.25
                for kc in keys_cell:
                    weights[kc] = weights.get(kc, 0.0) + area_quarter

    segs = dedupe_segments(raw_segments, pts, decimals=decimals)
    return pts, segs, h, weights, None


def merge_grids(coarse_pts, coarse_segs, plug_pts, plug_segs, decimals=12):
    merged_pts = dict(coarse_pts)
    merged_pts.update(plug_pts)
    raw_segments = []
    raw_segments.extend(coarse_segs)
    raw_segments.extend(plug_segs)
    merged_segs = dedupe_segments(raw_segments, merged_pts, decimals=decimals)
    return merged_pts, merged_segs


def make_basis(radial, angular, width=1.0):
    def radial_fn(x, y):
        r2 = x * x + y * y
        r = np.sqrt(r2)
        if radial == "gaus":
            return np.exp(-r2 / (width * width))
        if radial == "slater":
            return np.exp(-r / width)
        raise ValueError(f"Unknown radial '{radial}'")

    def angular_fn(x, y, cs=[1.0, 0.0, 0.0]):
        return 1.0*cs[0] + x*cs[1] + y*cs[2]
    return lambda x, y: angular_fn(x, y, cs=angular) * radial_fn(x, y)


def evaluate_integral_with_weights(points_dict, weights_dict, basis_fn, decimals=12):
    total = 0.0
    for p in points_dict.values():
        k = tuple(np.round(p, decimals=decimals))
        w = weights_dict.get(k, 0.0)
        v = basis_fn(p[0], p[1])
        total += v * w
    return float(total)


def generate_reference_grid(half, margin, h_coarse, ref_mult, basis_fn):
    h_ref = h_coarse / ref_mult
    domain_half = half + margin * h_coarse
    coords = np.arange(-domain_half, domain_half + 1e-9, h_ref)
    xv, yv = np.meshgrid(coords, coords, indexing="xy")
    vals = basis_fn(xv, yv)
    integral = float(vals.sum() * (h_ref * h_ref))
    return integral


def parse_args():
    p = argparse.ArgumentParser(description="Corrected Plug Debugger")
    p.add_argument("--n_sub", type=int, default=4)
    p.add_argument("--dmin", type=float, default=0.2)
    p.add_argument("--dmax", type=float, default=1.0)
    p.add_argument("--half", type=float, default=1.0)
    p.add_argument("--margin", type=int, default=1)
    return p.parse_args()


if __name__== "__main__":
    args = parse_args()
    
    # 1. Setup
    n_sub = args.n_sub
    h_coarse = 2 * args.half / n_sub
    basis_fn = make_basis("slater", [1.0, 0.0, 0.0], width=1.0)
    
    verts, core = build_square_plug(half=args.half)
    
    # 2. Generate Grids
    # Plug
    plug_pts, plug_segs, plug_tris = generate_plug(verts, core, n_sub, args.dmin, args.dmax)
    plug_weights = triangle_weights(plug_tris)
    
    # Coarse (Hole)
    coarse_pts_no_plug, coarse_segs_no_plug, _, coarse_weights_no_plug, _ = generate_regular_grid(
        args.half, n_sub, args.margin, remove_center=True
    )
    
    # Full Coarse (for reference of boundary values)
    _, _, _, coarse_weights_full, _ = generate_regular_grid(
        args.half, n_sub, args.margin, remove_center=False
    )

    # 3. Merge Weights
    merged_pts, merged_segs = merge_grids(coarse_pts_no_plug, coarse_segs_no_plug, plug_pts, plug_segs)
    merged_weights = dict(coarse_weights_no_plug)
    
    # Add plug weights
    for k, v in plug_weights.items():
        merged_weights[k] = merged_weights.get(k, 0.0) + v
        
    # Recovery check: For points on the boundary, they should have roughly area = h_coarse^2.
    # Coarse side gives ~0.5 * h^2. Plug side gives ~0.5 * h^2.
    
    # 4. Plot
    fig, ax = plt.subplots(figsize=(7, 7))
    plot_segments(ax, merged_segs, color="black", lw=0.5, alpha=0.3)
    
    # Annotate weights
    for k, p in merged_pts.items():
        w = merged_weights.get(k, 0.0)
        # Color code: Boundary vs Inner
        col = "red" if w > 0.26 or w < 0.24 else "blue" # arbitrary threshold for visual inspection
        ax.annotate(f"{w:.3f}", (p[0], p[1]), fontsize=7, color=col, ha="center", va="center")
        
    ax.set_aspect("equal")
    ax.set_title(f"Symmetric Plug (n={n_sub})\nCorner weights should be consistent")
    plt.tight_layout()
    plt.show()

    # 5. Sanity Check Numbers
    print("=== Weight Analysis ===")
    p_corner = (args.half, args.half)
    k_corner = tuple(np.round(p_corner, 12))
    print(f"Corner {k_corner} Weight: {merged_weights.get(k_corner, 0.0):.4f}")
    
    p_edge = (args.half, 0.0)
    k_edge = tuple(np.round(p_edge, 12))
    print(f"Edge Mid {k_edge} Weight: {merged_weights.get(k_edge, 0.0):.4f}")
    
    # Verify Symmetry on Edge
    w_vals = []
    # Check weights along the right edge x=half
    for y in np.linspace(-args.half, args.half, n_sub + 1):
        k = tuple(np.round([args.half, y], 12))
        if k in merged_weights:
            w_vals.append(merged_weights[k])
    
    print(f"Weights along edge: {[float(f'{x:.4f}') for x in w_vals]}")
    is_symmetric = np.allclose(w_vals, w_vals[::-1])
    print(f"Edge Symmetry Preserved: {is_symmetric}")