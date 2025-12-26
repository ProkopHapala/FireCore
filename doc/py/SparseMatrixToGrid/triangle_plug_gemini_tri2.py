#!/usr/bin/env python3
"""
Integrated Plug + Coarse Grid Generator (Gemini-based geometry, FireCore CLI/debug).

Features:
1. "Onion" Triangular Plug (Triangles inside) with symmetric topology.
2. Coarse Grid (Quads outside) with configurable margin layers.
3. Seamless merging at the boundary; area/3 (tris) and area/4 (quads) weighting.
4. CLI parity with triangle_plug_debug.py (angular coeffs, weight-func, ref grid, labels).
5. Diagnostic prints: weight sums, plug/coarse/hybrid integrals vs. reference.
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse
import math
import sys

# ==========================================
# 0. Utilities and basis
# ==========================================

def make_weight_func(name, dmin, dmax):
    if name is None or name == "linear":
        return None
    if name == "cosine":
        return lambda u: dmax - (dmax - dmin) * (0.5 - 0.5 * np.cos(math.pi * u))
    raise ValueError(f"Unknown weight_func '{name}'")


def mapped_t_values(n_layers, dmin=1.0, dmax=1.0, weight_func=None):
    if n_layers == 0:
        return np.array([0.0, 1.0])
    u = np.linspace(0.0, 1.0, n_layers)
    if weight_func is None:
        weights = dmax - (dmax - dmin) * u
    else:
        weights = weight_func(u)
    weights = np.maximum(weights, 1e-12)
    weights = weights / weights.sum()
    levels = np.concatenate([[0.0], np.cumsum(weights)])
    levels[-1] = 1.0
    return levels


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
        return 1.0 * cs[0] + x * cs[1] + y * cs[2]
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
    coords = np.arange(-domain_half, domain_half + 1e-12, h_ref)
    xv, yv = np.meshgrid(coords, coords, indexing="xy")
    vals = basis_fn(xv, yv)
    integral = float(vals.sum() * (h_ref * h_ref))
    return integral


# ==========================================
# 1. Plug Generation ("Onion" / Triangular)
# ==========================================

def get_onion_coords(n_sub, dmin, dmax, half, weight_func=None):
    """
    Generates coordinates for the plug.
    Layer 0 = Boundary (Linear spacing to match Coarse Grid).
    Layer N = Core.
    """
    logical_map = {}
    t_levels = mapped_t_values(n_sub, dmin=dmin, dmax=dmax, weight_func=weight_func)
    center = np.array([0.0, 0.0])
    
    corners = [
        np.array([-half, -half]),
        np.array([ half, -half]),
        np.array([ half,  half]),
        np.array([-half,  half])
    ]
    
    for l in range(n_sub + 1):
        t = t_levels[l]
        
        if l == n_sub:
            logical_map[("CORE",)] = center
            continue
            
        for q in range(4):
            v_start = corners[q]
            v_end   = corners[(q + 1) % 4]
            
            p_start = (1 - t) * v_start + t * center
            p_end   = (1 - t) * v_end   + t * center
            
            n_seg = n_sub - l
            for i in range(n_seg):
                s = i / n_seg
                pt = (1 - s) * p_start + s * p_end
                logical_map[(q, l, i)] = pt

    return logical_map

def get_onion_topology(n_sub, coord_map):
    """
    Connects the dots to form Triangles.
    Returns list of coords: [(p1, p2, p3), ...]
    """
    tris = []
    
    # Helper to lookup coord
    def get_pt(q, l, i):
        # Normalize
        q_norm = q % 4
        
        # Max index in this layer
        max_idx = n_sub - l
        
        # Wrap corners
        if i == max_idx:
            q_norm = (q_norm + 1) % 4
            i = 0
            
        if l == n_sub:
            return coord_map[("CORE",)]
            
        return coord_map[(q_norm, l, i)]

    for q in range(4):
        for l in range(n_sub):
            n_seg = n_sub - l
            
            for i in range(n_seg):
                # Current layer points
                p_out_1 = get_pt(q, l, i)
                p_out_2 = get_pt(q, l, i + 1)
                
                # Next layer points
                if l == n_sub - 1:
                    # Connect to core
                    p_core = coord_map[("CORE",)]
                    tris.append((p_out_1, p_out_2, p_core))
                else:
                    # Standard strip
                    p_in_1 = get_pt(q, l + 1, i)
                    
                    # Triangle 1 (Base on outer)
                    tris.append((p_out_1, p_out_2, p_in_1))
                    
                    # Triangle 2 (Inverted, Base on inner)
                    # Only if not the last segment of the strip
                    if i < n_seg - 1:
                        p_in_2 = get_pt(q, l + 1, i + 1)
                        tris.append((p_out_2, p_in_2, p_in_1))
                        
    return tris

# ==========================================
# 2. Coarse Grid Generation
# ==========================================

def get_coarse_grid(n_sub, half, margin_cells):
    """
    Generates Quads for the area surrounding the plug.
    Spacing h is determined by n_sub and half to match plug boundary.
    """
    h = (2.0 * half) / n_sub
    
    # Extent
    min_val = -half - margin_cells * h
    max_val =  half + margin_cells * h
    
    # Number of points along one side of the total domain
    # Plug takes n_sub segments. Margin takes 2*margin_cells segments.
    N_total = n_sub + 2 * margin_cells + 1
    
    coords = np.linspace(min_val, max_val, N_total)
    
    quads = []
    
    for i in range(len(coords) - 1):
        for j in range(len(coords) - 1):
            x0, x1 = coords[i], coords[i+1]
            y0, y1 = coords[j], coords[j+1]
            
            cx, cy = (x0+x1)/2, (y0+y1)/2
            
            # Check if this quad is strictly inside the plug area
            # If center is inside [-half, half], skip it
            if abs(cx) < half - 1e-9 and abs(cy) < half - 1e-9:
                continue
            
            # Create Quad (Counter-clockwise)
            p00 = np.array([x0, y0])
            p10 = np.array([x1, y0])
            p11 = np.array([x1, y1])
            p01 = np.array([x0, y1])
            
            quads.append((p00, p10, p11, p01))
            
    return quads

# ==========================================
# 3. Merging & Weight Calculation
# ==========================================

def merge_and_analyze(plug_tris, coarse_quads):
    """
    Merge geometries by coordinate rounding.
    Returns:
      unique_pts: {coord_key: np.array}
      weights_total: {coord_key: float}  (plug + coarse)
      weights_plug:  {coord_key: float}
      weights_coarse:{coord_key: float}
      elements: list of (type_str, list_of_keys, area, centroid)
    """
    unique_pts = {}
    weights_total = {}
    weights_plug = {}
    weights_coarse = {}
    elements = []

    def get_key(pt):
        return tuple(np.round(pt, 10))

    def reg(pt):
        k = get_key(pt)
        if k not in unique_pts:
            unique_pts[k] = pt
            weights_total[k] = 0.0
            weights_plug[k] = 0.0
            weights_coarse[k] = 0.0
        return k

    # Plug triangles
    for tri in plug_tris:
        keys = [reg(p) for p in tri]
        pts = [unique_pts[k] for k in keys]
        p0, p1, p2 = pts
        area = 0.5 * abs((p1[0]-p0[0])*(p2[1]-p0[1]) - (p1[1]-p0[1])*(p2[0]-p0[0]))
        c = np.mean(pts, axis=0)
        w_share = area / 3.0
        for k in keys:
            weights_total[k] += w_share
            weights_plug[k] += w_share
        elements.append(("tri", keys, area, c))

    # Coarse quads
    for quad in coarse_quads:
        keys = [reg(p) for p in quad]
        pts_arr = np.array([unique_pts[k] for k in keys])
        c = np.mean(pts_arr, axis=0)
        w = np.linalg.norm(pts_arr[1] - pts_arr[0])
        h = np.linalg.norm(pts_arr[3] - pts_arr[0])
        area = w * h
        w_share = area / 4.0
        for k in keys:
            weights_total[k] += w_share
            weights_coarse[k] += w_share
        elements.append(("quad", keys, area, c))

    return unique_pts, weights_total, weights_plug, weights_coarse, elements

# ==========================================
# 4. Visualization
# ==========================================

def list_points_and_duplicates(unique_pts, decimals=12):
    print("\n=== Points (rounded) ===")
    for k in sorted(unique_pts.keys()):
        print(tuple(float(f"{x:.12g}") for x in k))
    print("\n=== Duplicate keys (count>1) ===")
    print("None")
    return False


def annotate_points(ax, unique_pts, weights=None, fontsize=8, color="red", decimals=12):
    sorted_items = sorted(unique_pts.items())
    for idx, (k, p) in enumerate(sorted_items):
        text = f"{weights.get(k, 0.0):.3g}" if weights is not None else f"{idx}"
        ax.annotate(text, (p[0], p[1]), color=color, fontsize=fontsize, ha="center", va="center")


def print_boundary_weights(unique_pts, weights, half, tol=1e-12):
    items = []
    for k, p in unique_pts.items():
        x, y = p
        if abs(abs(x) - half) < tol or abs(abs(y) - half) < tol:
            items.append((float(x), float(y), float(weights.get(k, 0.0))))
    items.sort(key=lambda t: (t[0], t[1]))
    print("=== Plug boundary weights (|x|=half or |y|=half) ===")
    for x, y, w in items:
        print(f"({x:.6g}, {y:.6g}) -> {w:.6g}")
    print("")


def print_element_areas(elements, unique_pts, max_items=50):
    print("=== Element areas (centroid) ===")
    for idx, (etype, keys, area, cent) in enumerate(elements[:max_items]):
        cx, cy = cent
        print(f"{idx:04d} {etype} area={area:.6g} centroid=({cx:.6g},{cy:.6g})")
    if len(elements) > max_items:
        print(f"... {len(elements)-max_items} more elements not shown")
    print("")


# ==========================================
# 4. Main
# ==========================================

def parse_args():
    p = argparse.ArgumentParser(description="Hybrid grid with symmetric plug (Gemini geometry) + FireCore CLI.")
    p.add_argument("--weight-func", dest="weight_func_name", choices=["linear", "cosine"], default="linear", help="Shape of segment length ramp from boundary to core.")
    p.add_argument("--n_sub",       type=int,   default=4,    help="Subdivision count along plug edge.")
    p.add_argument("--dmin",        type=float, default=0.2,  help="Relative segment length near core.")
    p.add_argument("--dmax",        type=float, default=1.0,  help="Relative segment length near boundary.")
    p.add_argument("--half",        type=float, default=1.0,  help="Half-size of square plug.")
    p.add_argument("--list-points", type=int,   default=0,    help="Print rounded unique points and duplicate counts, then exit.")
    p.add_argument("--font-size",   type=float, default=12.0, help="Font size for point labels.")
    p.add_argument("--margin",      type=int,   default=3,    help="Margin (in coarse-grid steps) around plug.")
    p.add_argument("--radial",      choices=["gaus", "slater"], default="slater", help="Radial part of basis.")
    p.add_argument("--angular",     type=str, default="1.0,0.0,0.0", help="Angular coefficients for s,px,py (e.g., '1.0,0.0,0.0').")
    p.add_argument("--width",       type=float, default=1.0, help="Width parameter for basis (sigma for gaus, 1/beta for slater).")
    p.add_argument("--ref-mult",    type=int,   default=4, help="Reference grid refinement factor over coarse spacing.")
    p.add_argument("--label",       choices=["id", "weight"], default="weight", help="Annotate points by index or by computed weight.")
    return p.parse_args()


def main():
    args = parse_args()
    n_sub = args.n_sub
    weight_func = make_weight_func(args.weight_func_name, args.dmin, args.dmax)

    # Basis
    angular_coeffs = [float(c) for c in args.angular.split(",")]
    if len(angular_coeffs) != 3:
        raise ValueError(f"Expected 3 angular coefficients, got {len(angular_coeffs)}")
    basis_fn = make_basis(args.radial, angular_coeffs, width=args.width)

    h_coarse = 2 * args.half / n_sub

    # Plug
    p_map = get_onion_coords(n_sub, args.dmin, args.dmax, args.half, weight_func=weight_func)
    plug_tris = get_onion_topology(n_sub, p_map)

    # Coarse grid
    coarse_quads = get_coarse_grid(n_sub, args.half, args.margin)

    # Merge + weights
    unique_pts, weights_total, weights_plug, weights_coarse, elements = merge_and_analyze(plug_tris, coarse_quads)

    # Integrals
    coarse_int = evaluate_integral_with_weights(unique_pts, weights_coarse, basis_fn)
    plug_int = evaluate_integral_with_weights(unique_pts, weights_plug, basis_fn)
    hybrid_int = evaluate_integral_with_weights(unique_pts, weights_total, basis_fn)
    ref_int = generate_reference_grid(args.half, args.margin, h_coarse, args.ref_mult, basis_fn)

    # Areas
    domain_half = args.half + args.margin * h_coarse
    ref_area = (2 * domain_half) * (2 * domain_half)
    print(f"=== Grid Statistics ===")
    print(f"Coarse grid quads:        {len(coarse_quads)}")
    print(f"Plug triangles:           {len(plug_tris)}")
    print(f"Unique points:            {len(unique_pts)}")
    print("")

    print(f"=== Area / weight sums ===")
    print(f"Reference area:           {ref_area:.12g}")
    print(f"Coarse weights sum:       {sum(weights_coarse.values()):.12g}")
    print(f"Plug weights sum:         {sum(weights_plug.values()):.12g}")
    print(f"Hybrid weights sum:       {sum(weights_total.values()):.12g}")
    print("")
    print_boundary_weights(unique_pts, weights_total, args.half)
    print_element_areas(elements, unique_pts)

    print(f"=== Integral Comparison ===")
    print(f"Coarse only:              {coarse_int:.12g}")
    print(f"Plug only:                {plug_int:.12g}")
    print(f"Hybrid (coarse+plug):     {hybrid_int:.12g}")
    print(f"High-resolution reference:{ref_int:.12g}")
    print(f"Coarse error: abs={abs(coarse_int - ref_int):.6g}, rel={(abs(coarse_int - ref_int)/(abs(ref_int)+1e-15)):.6g}")
    print(f"Hybrid error: abs={abs(hybrid_int - ref_int):.6g}, rel={(abs(hybrid_int - ref_int)/(abs(ref_int)+1e-15)):.6g}")

    # Optional point list
    if args.list_points:
        list_points_and_duplicates(unique_pts)
        sys.exit(0)

    # Plot
    fig, ax = plt.subplots(figsize=(8, 8))
    # Draw elements (light)
    for etype, keys, area, cent in elements:
        color = "green" if etype == "tri" else "blue"
        alpha = 0.08
        pts = [unique_pts[k] for k in keys]
        poly = plt.Polygon(pts, closed=True, fill=True, facecolor=color, edgecolor='k', alpha=alpha, lw=0.5)
        ax.add_patch(poly)
        for k in keys:
            p_v = unique_pts[k]
            ax.plot([p_v[0], cent[0]], [p_v[1], cent[1]], color='red', lw=0.3, alpha=0.3)
        # annotate element area at centroid
        ax.annotate(f"{area:.3g}", (cent[0], cent[1]), color="darkslategray", fontsize=6, ha="center", va="center")

    ann_weights = weights_total if args.label == "weight" else None
    annotate_points(ax, unique_pts, weights=ann_weights, fontsize=args.font_size, color="black")

    rect = plt.Rectangle((-args.half, -args.half), 2*args.half, 2*args.half, fill=False, edgecolor='k', lw=2, linestyle='--')
    ax.add_patch(rect)

    ax.set_aspect('equal')
    ax.set_title(f"Hybrid Grid (n_sub={n_sub}, margin={args.margin})")
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()