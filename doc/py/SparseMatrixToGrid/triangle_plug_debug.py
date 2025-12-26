#!/usr/bin/env python3
"""
Quick visual debugger for the 2D plug (square split into 4 triangles).

Now: start from a **linear tessellation** of each triangle into similar small
triangles. Keep the outer edge uniformly subdivided; for each base point on
the outer edge, the points along the segment to the core are shifted using an
arbitrary mapping f:[0,1]->[0,1]. This lets us cluster nodes toward/away from
the core without moving the outer boundary.

Usage (default): python3 triangle_plug_debug.py
"""

import numpy as np
import matplotlib.pyplot as plt
import math
import argparse
import sys


# DEBUG helper: keep explicit prints to trace parameters and generated grids
def debug_print(label, arr):
    print(f"{label}: shape={arr.shape}, min={arr.min():.4f}, max={arr.max():.4f}, first={arr.flatten()[:5]}")


def mapped_t_values(n_layers, dmin=1.0, dmax=1.0, weight_func=None):
    """
    Return t-levels (size n_layers+1) so that segment lengths increase from
    core (small ~ dmin) to boundary (large ~ dmax) after normalization.

    weights[i] is the desired length of segment i (boundary outward to core).
    They are normalized to sum to 1 so total length equals the edge length.
    """
    if n_layers == 0:
        return np.array([0.0, 1.0])
    u = np.linspace(0.0, 1.0, n_layers)
    if weight_func is None:
        # u=0 -> boundary segment (largest), u=1 -> near core (smallest)
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
    Linear tessellation with uniform base AB and mapped core edges AC, BC.
    For each layer l (0=base, n_sub=core):
      t = mapped_t_values(...) controls position along AC/BC; AC and BC get denser toward C.
      The line between AC(t) and BC(t) is uniformly subdivided within the layer.
    """
    t_levels = mapped_t_values(n_sub, dmin=dmin, dmax=dmax, weight_func=weight_func)
    points = {}
    triangles = []

    # Generate points layer by layer
    for l in range(n_sub + 1):
        t = t_levels[l]
        a_l = (1 - t) * v0 + t * core  # point on AC
        b_l = (1 - t) * v1 + t * core  # point on BC
        mmax = n_sub - l
        for m in range(mmax + 1):
            s = 0.0 if mmax == 0 else m / mmax  # uniform within the layer
            pt = (1 - s) * a_l + s * b_l
            points[(l, m)] = pt

    # Build edge list for wireframe triangles and triangle connectivity
    segments = []
    for l in range(n_sub + 1):
        mmax = n_sub - l
        for m in range(mmax + 1):
            p = points[(l, m)]
            # within-layer edge
            if m + 1 <= mmax:
                q = points[(l, m + 1)]
                segments.append((p, q))
            # edge to next layer (vertical)
            if l + 1 <= n_sub and m <= n_sub - (l + 1):
                q = points[(l + 1, m)]
                segments.append((p, q))
            # diagonal edge to next layer
            if l + 1 <= n_sub and m - 1 >= 0 and m - 1 <= n_sub - (l + 1):
                q = points[(l + 1, m - 1)]
                segments.append((p, q))
            # small triangles for weights (two per rhombus where applicable)
            if l + 1 <= n_sub and m + 1 <= mmax:
                # Triangle (l,m)-(l+1,m)-(l,m+1)
                a = points[(l, m)]
                b = points[(l + 1, m)]
                c = points[(l, m + 1)]
                triangles.append((a, b, c))
                # Triangle (l+1,m)-(l+1,m-1)-(l,m+1) when valid
                if m - 1 >= 0 and m - 1 <= n_sub - (l + 1):
                    d = points[(l + 1, m - 1)]
                    triangles.append((b, d, c))
    # Deduplicate triangles by rounding (to avoid overlaps on shared edges)
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

    return points, segments, unique_triangles


def build_square_plug(center=np.array([0.0, 0.0]), half=1.0):
    """Return square vertices ordered ccw and center."""
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


def dedupe_points(all_points, decimals=12):
    unique = {}
    for p in all_points:
        key = tuple(np.round(p, decimals=decimals))
        if key not in unique:
            unique[key] = p
    return unique


def dedupe_segments(raw_segments, unique_pts, decimals=12):
    seg_keys = set()
    unique_segs = []
    for a, b in raw_segments:
        ka = tuple(np.round(a, decimals=decimals))
        kb = tuple(np.round(b, decimals=decimals))
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
        area = 0.5 * abs((b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0]))
        contrib = area / 3.0
        for p in (a, b, c):
            k = tuple(np.round(p, decimals=decimals))
            weights[k] = weights.get(k, 0.0) + contrib
    return weights


def generate_plug(verts, core, n_sub, dmin, dmax, weight_func=None, decimals=12):
    """Generate all points/segments/triangles for 4 triangles, dedup shared edge nodes/segments."""
    unique_pts = {}
    all_segments = []
    all_triangles = []
    for i in range(4):
        v0 = verts[i]
        v1 = verts[(i + 1) % 4]
        pts, segs, tris = triangle_subdivision(v0, v1, core, n_sub=n_sub, dmin=dmin, dmax=dmax, weight_func=weight_func, decimals=decimals)
        for p in pts.values():
            key = tuple(np.round(p, decimals=decimals))
            if key not in unique_pts:
                unique_pts[key] = p
        all_segments.extend(segs)
        all_triangles.extend(tris)

    # dedupe segments by rounded endpoints
    unique_segs = dedupe_segments(all_segments, unique_pts, decimals=decimals)
    return unique_pts, unique_segs, all_triangles


def generate_regular_grid(half, n_sub, margin, remove_center=False, decimals=12):
    """
    Generate square grid with spacing matching plug boundary.
    Total points per side = n_sub + 2*margin + 1 (so boundary aligns with plug).
    If remove_center=True, points inside plug square [-half, half]^2 are skipped.
    """
    h = 2 * half / n_sub  # spacing
    n_side = n_sub + 2 * margin + 1
    coords = np.linspace(-half - margin * h, half + margin * h, n_side)
    pts = {}
    for x in coords:
        for y in coords:
            # keep boundary; remove only strict interior of plug
            if remove_center and (abs(x) < half - 1e-12 and abs(y) < half - 1e-12):
                continue
            key = tuple(np.round([x, y], decimals=decimals))
            pts[key] = np.array([x, y])

    raw_segments = []
    weights = {}
    cells = []  # list of (keys_cell, center)
    for ix, x in enumerate(coords):
        for iy, y in enumerate(coords):
            if remove_center and (abs(x) < half - 1e-12 and abs(y) < half - 1e-12):
                continue
            key = tuple(np.round([x, y], decimals=decimals))
            if key not in pts:
                continue
            # right neighbor
            if ix + 1 < len(coords):
                x2 = coords[ix + 1]
                if not (remove_center and (abs(x2) < half - 1e-12 and abs(y) < half - 1e-12)):
                    key2 = tuple(np.round([x2, y], decimals=decimals))
                    if key2 in pts:
                        raw_segments.append((pts[key], pts[key2]))
            # up neighbor
            if iy + 1 < len(coords):
                y2 = coords[iy + 1]
                if not (remove_center and (abs(x) < half - 1e-12 and abs(y2) < half - 1e-12)):
                    key2 = tuple(np.round([x, y2], decimals=decimals))
                    if key2 in pts:
                        raw_segments.append((pts[key], pts[key2]))
            # square cell contribution (only if all four corners exist)
            if ix + 1 < len(coords) and iy + 1 < len(coords):
                x2 = coords[ix + 1]
                y2 = coords[iy + 1]
                keys_cell = [
                    tuple(np.round([x, y], decimals=decimals)),
                    tuple(np.round([x2, y], decimals=decimals)),
                    tuple(np.round([x, y2], decimals=decimals)),
                    tuple(np.round([x2, y2], decimals=decimals)),
                ]
                if all(k in pts for k in keys_cell):
                    area_quarter = (h * h) * 0.25
                    for kc in keys_cell:
                        weights[kc] = weights.get(kc, 0.0) + area_quarter
                    cx = 0.5 * (x + x2)
                    cy = 0.5 * (y + y2)
                    cells.append((keys_cell, cx, cy))

    segs = dedupe_segments(raw_segments, pts, decimals=decimals)
    return pts, segs, h, weights, cells


def merge_grids(coarse_pts, coarse_segs, plug_pts, plug_segs, decimals=12):
    merged_pts = dict(coarse_pts)
    merged_pts.update(plug_pts)  # plug supersedes (shares) boundary coords
    raw_segments = []
    raw_segments.extend(coarse_segs)
    raw_segments.extend(plug_segs)
    merged_segs = dedupe_segments(raw_segments, merged_pts, decimals=decimals)
    return merged_pts, merged_segs


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
        if weights is not None:
            text = f"{weights.get(k, 0.0):.3g}"
        else:
            text = f"{idx}"
        ax.annotate(text, (p[0], p[1]), color=color, fontsize=fontsize, ha="center", va="center")
    xs = [p[0] for _, p in sorted_items]
    ys = [p[1] for _, p in sorted_items]
    #ax.scatter(xs, ys, s=10, color="red", zorder=3, alpha=0.8)


def make_weight_func(name, dmin, dmax):
    if name is None or name == "linear":
        return None  # use built-in linear ramp from dmin->dmax
    if name == "cosine":
        return lambda u: dmax - (dmax - dmin) * (0.5 - 0.5 * np.cos(math.pi * u))
    raise ValueError(f"Unknown weight_func '{name}'")


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
        return 1.0*cs[0] + x*cs[1] + y*cs[2]  # s + px + py combination
    return lambda x, y: angular_fn(x, y, cs=angular) * radial_fn(x, y)


def evaluate_integral(points_dict, basis_fn, weight):
    """Legacy function for backward compatibility."""
    vals = []
    for p in points_dict.values():
        vals.append(basis_fn(p[0], p[1]))
    vals = np.array(vals)
    return float(vals.sum() * weight), vals


def evaluate_integral_with_weights(points_dict, weights_dict, basis_fn, decimals=12):
    """Evaluate integral using per-point weights."""
    total = 0.0
    vals = []
    for p in points_dict.values():
        k = tuple(np.round(p, decimals=decimals))
        w = weights_dict.get(k, 0.0)
        v = basis_fn(p[0], p[1])
        vals.append(v)
        total += v * w
    return float(total), np.array(vals)


def generate_reference_grid(half, margin, h_coarse, ref_mult, basis_fn):
    h_ref = h_coarse / ref_mult
    domain_half = half + margin * h_coarse
    coords = np.arange(-domain_half, domain_half + 1e-12, h_ref)
    xv, yv = np.meshgrid(coords, coords, indexing="xy")
    vals = basis_fn(xv, yv)
    integral = float(vals.sum() * (h_ref * h_ref))
    return integral


def generate_reference_grid_for_area(half, h_coarse, ref_mult, basis_fn):
    """Generate reference integral for a specific area (e.g., plug area)."""
    h_ref = h_coarse / ref_mult
    coords = np.linspace(-half, half, int(2*half/h_ref) + 1)
    xv, yv = np.meshgrid(coords, coords, indexing="xy")
    vals = basis_fn(xv, yv)
    integral = float(vals.sum() * (h_ref * h_ref))
    return integral


def parse_args():
    p = argparse.ArgumentParser(description="Visualize square plug triangle subdivision with non-uniform AC/BC spacing.")
    p.add_argument("--weight-func", dest="weight_func_name", choices=["linear", "cosine"], default="linear", help="Shape of segment length ramp from boundary to core.")
    p.add_argument("--n_sub",       type=int,   default=4,    help="Subdivision count along base edge AB (uniform).")
    p.add_argument("--dmin",        type=float, default=0.2,  help="Relative segment length near core along AC/BC.")
    p.add_argument("--dmax",        type=float, default=1.0,  help="Relative segment length near boundary along AC/BC.")
    p.add_argument("--half",        type=float, default=1.0,  help="Half-size of square plug.")
    p.add_argument("--list-points", type=int,   default=0,    help="Print rounded unique points and duplicate counts, then exit.")
    p.add_argument("--font-size",   type=float, default=12.0, help="Font size for point labels.")
    p.add_argument("--margin",      type=int,   default=3,    help="Margin (in coarse-grid steps) around plug.")
    p.add_argument("--radial",      choices=["gaus", "slater"], default="slater", help="Radial part of basis.")
    p.add_argument("--angular",     type=str, default="1.0,0.0,0.0", help="Angular coefficients for s,px,py (e.g., '1.0,0.0,0.0' for pure s, '1.0,1.0,1.0' for sp2).")
    p.add_argument("--width",       type=float, default=1.0, help="Width parameter for basis (sigma for gaus, 1/beta for slater).")
    p.add_argument("--ref-mult",    type=int,   default=4, help="Reference grid refinement factor over coarse spacing.")
    p.add_argument("--label",       choices=["id", "weight"], default="weight", help="Annotate points by index or by computed weight.")
    return p.parse_args()


if __name__== "__main__":
    args = parse_args()
    n_sub = args.n_sub
    dmin = args.dmin
    dmax = args.dmax
    weight_func = make_weight_func(args.weight_func_name, dmin, dmax)
    
    # Parse angular coefficients
    angular_coeffs = [float(c) for c in args.angular.split(',')]
    if len(angular_coeffs) != 3:
        raise ValueError(f"Expected 3 angular coefficients, got {len(angular_coeffs)}")
    
    basis_fn = make_basis(args.radial, angular_coeffs, width=args.width)
    h_coarse = 2 * args.half / n_sub

    verts, core = build_square_plug(half=args.half)
    fig, ax = plt.subplots(figsize=(6, 6))

    t_levels = mapped_t_values(n_sub, dmin=dmin, dmax=dmax, weight_func=weight_func)
    debug_print("t_levels (0=boundary -> 1=core)", t_levels)

    plug_pts, plug_segs, plug_tris = generate_plug(verts, core, n_sub=n_sub, dmin=dmin, dmax=dmax, weight_func=weight_func)
    plug_weights = triangle_weights(plug_tris)
    
    # Generate two versions of coarse grid
    coarse_pts_full, coarse_segs_full, _, coarse_weights_full, _ = generate_regular_grid(half=args.half, n_sub=n_sub, margin=args.margin, remove_center=False)
    coarse_pts_no_plug, coarse_segs_no_plug, _, coarse_weights_no_plug, _ = generate_regular_grid(half=args.half, n_sub=n_sub, margin=args.margin, remove_center=True)
    
    merged_pts, merged_segs = merge_grids(coarse_pts_no_plug, coarse_segs_no_plug, plug_pts, plug_segs)
    # merged per-point weights: start from coarse grid without plug, then add plug contributions,
    # and finally add back coarse contributions for plug boundary points that still exist
    merged_weights = dict(coarse_weights_no_plug)
    for k, v in plug_weights.items():
        merged_weights[k] = merged_weights.get(k, 0.0) + v
    # add coarse weights for points that were removed from coarse_no_plug but are present in plug (boundary)
    for k, v in coarse_weights_full.items():
        if k in merged_weights and k not in coarse_weights_no_plug:
            merged_weights[k] += v
    duplicates = 0

    # Integrals - compute all three cases (using per-point weights)
    # 1. Coarse grid alone (full grid)
    coarse_int, _ = evaluate_integral_with_weights(coarse_pts_full, coarse_weights_full, basis_fn)
    
    # 2. Hybrid grid (coarse + plug) using merged weights
    hybrid_int, _ = evaluate_integral_with_weights(merged_pts, merged_weights, basis_fn)
    
    # 3. High-resolution reference grid
    ref_int = generate_reference_grid(args.half, args.margin, h_coarse, args.ref_mult, basis_fn)
    
    # Plug diagnostics
    plug_ref_int = generate_reference_grid_for_area(args.half, h_coarse, args.ref_mult, basis_fn)
    plug_coarse_pts = {}
    plug_half = args.half
    for k, p in coarse_pts_full.items():
        if abs(p[0]) <= plug_half + 1e-12 and abs(p[1]) <= plug_half + 1e-12:
            plug_coarse_pts[k] = p
    plug_coarse_int, _ = evaluate_integral_with_weights(plug_coarse_pts, coarse_weights_full, basis_fn)
    plug_part, _ = evaluate_integral_with_weights(plug_pts, plug_weights, basis_fn)

    # Weight (area) sums for sanity check
    coarse_w_sum = sum(coarse_weights_full.values())
    hybrid_w_sum = sum(merged_weights.values())
    plug_w_sum = sum(plug_weights.values())
    plug_coarse_w_sum = sum(coarse_weights_full.get(k, 0.0) for k in plug_coarse_pts.keys())
    domain_half = args.half + args.margin * h_coarse
    ref_area = (2 * domain_half) * (2 * domain_half)
    
    print(f"=== Grid Statistics ===")
    print(f"Coarse grid points:        {len(coarse_pts_full)}")
    print(f"Coarse grid (no plug):     {len(coarse_pts_no_plug)}")
    print(f"Plug points:               {len(plug_pts)}")
    print(f"Plug area coarse points:   {len(plug_coarse_pts)}")
    print(f"Hybrid grid points:        {len(merged_pts)}")
    print(f"")
    
    print(f"=== Area / weight sums ===")
    print(f"Reference area:            {ref_area:.12g}")
    print(f"Coarse weights sum:        {coarse_w_sum:.12g}")
    print(f"Hybrid weights sum:        {hybrid_w_sum:.12g}")
    print(f"Plug weights sum:          {plug_w_sum:.12g}")
    print(f"Plug coarse weights sum:   {plug_coarse_w_sum:.12g}")
    print(f"")
    
    print(f"=== Area Contributions ===")
    print(f"Plug area (coarse grid):   {plug_coarse_int:.12g}")
    print(f"Plug area (plug weights):  {plug_part:.12g}")
    print(f"Plug area (reference):     {plug_ref_int:.12g}")
    print(f"Coarse error:              {abs(plug_coarse_int - plug_ref_int):.6g}")
    print(f"Plug error:                {abs(plug_part - plug_ref_int):.6g}")
    print(f"")
    
    print(f"=== Integral Comparison ===")
    print(f"Coarse grid only:           {coarse_int:.12g}")
    print(f"Hybrid grid (coarse+plug):  {hybrid_int:.12g}")
    print(f"High-resolution reference:  {ref_int:.12g}")
    print(f"")
    print(f"Coarse error:   abs={abs(coarse_int - ref_int):.6g}, rel={(abs(coarse_int - ref_int)/(abs(ref_int)+1e-15)):.6g}")
    print(f"Hybrid error:   abs={abs(hybrid_int - ref_int):.6g}, rel={(abs(hybrid_int - ref_int)/(abs(ref_int)+1e-15)):.6g}")
    print(f"Improvement:    {(abs(coarse_int - ref_int) - abs(hybrid_int - ref_int))/abs(coarse_int - ref_int)*100:.2f}% reduction in error")

    if args.list_points:
        dup = list_points_and_duplicates(merged_pts)
        if dup:
            sys.exit(1)

    plot_segments(ax, merged_segs, color="C0", lw=1.0, alpha=0.2)
    ann_weights = merged_weights if args.label == "weight" else None
    annotate_points(ax, merged_pts, weights=ann_weights, fontsize=args.font_size, color="red")

    # Outline square and core
    outline = verts + [verts[0]]
    ax.plot([p[0] for p in outline], [p[1] for p in outline], "k--", lw=1.0, label="plug boundary")
    ax.scatter([core[0]], [core[1]], color="k", s=40, label="core")

    ax.set_aspect("equal", "box")
    ax.set_title(f"Square plug subdivided into linear grids (n_sub={n_sub})")
    ax.legend()
    ax.grid(True, alpha=0.2)
    plt.tight_layout()
    plt.show()
