#!/usr/bin/env python3
"""
Triangular Plug Generator with Symmetric "Onion" Layering.

Features:
- Generates a grid that decreases resolution from Boundary -> Core.
- Uses Index-based topology (Quadrant, Layer, Index) to ensure perfect stitching 
  without messy float deduplication on the boundary.
- Visualizes the "Weight Attribution" by drawing lines from vertices to triangle centroids.
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse

# ==========================================
# 1. Topology & Geometry Helpers
# ==========================================

def get_point_key(quadrant, layer, index, n_sub):
    """
    Normalize point indices to handle shared boundaries and the single core point.
    
    Layers go from 0 (boundary) to n_sub (core).
    Points in layer L range from index 0 to (n_sub - L).
    
    Rules:
    1. Core: Layer == n_sub -> Key is unique 'CORE'.
    2. Boundary sharing: The last point of Quadrant Q is the first point of Quadrant Q+1.
       (Q, L, max_idx) -> (Q+1, L, 0)
    """
    # 1. Check for Core
    if layer == n_sub:
        return ("CORE",)
    
    # Max index for this layer
    max_idx = n_sub - layer
    
    # 2. Normalize Quadrant/Index
    q = quadrant % 4
    idx = index
    
    # If we are at the end of a sector, wrap to the beginning of the next
    if idx == max_idx:
        q = (q + 1) % 4
        idx = 0
        
    return (q, layer, idx)

def get_coordinates(n_sub, dmin, dmax, half, center_pos):
    """
    Calculate coordinates for all logical points.
    Returns: dict { point_key : np.array([x, y]) }
    """
    coords = {}
    cx, cy = center_pos
    
    # Define the 4 corners of the box relative to center
    corners = [
        np.array([cx - half, cy - half]), # BL
        np.array([cx + half, cy - half]), # BR
        np.array([cx + half, cy + half]), # TR
        np.array([cx - half, cy + half])  # TL
    ]
    
    # Radial spacing function (t=0 boundary, t=1 core)
    # Using linear mapping for simplicity, can be swapped for cosine/etc.
    # To match your previous logic: weights scale from dmax (boundary) to dmin (core).
    
    # Generate radial t values
    u_vals = np.linspace(0.0, 1.0, n_sub)
    weights = dmax - (dmax - dmin) * u_vals
    weights = weights / weights.sum()
    t_levels = np.concatenate([[0.0], np.cumsum(weights)])
    t_levels[-1] = 1.0 # Force exact 1.0
    
    # Store Core
    coords[("CORE",)] = np.array(center_pos)
    
    for q in range(4):
        v_start = corners[q]
        v_end   = corners[(q + 1) % 4]
        
        for l in range(n_sub):
            t = t_levels[l]
            
            # The segment for this layer moves towards center
            # Start/End points of the "Arc" (straight line in this case)
            p_layer_start = (1 - t) * v_start + t * center_pos
            p_layer_end   = (1 - t) * v_end   + t * center_pos
            
            # Number of segments in this layer = n_sub - l
            # Number of points = n_sub - l + 1
            n_segments = n_sub - l
            
            # Generate points along the line
            # We only generate indices 0 to n_segments-1 for this quadrant.
            # The last point (index n_segments) is handled by the next quadrant's index 0.
            for i in range(n_segments):
                s = i / n_segments
                pt = (1 - s) * p_layer_start + s * p_layer_end
                
                key = get_point_key(q, l, i, n_sub)
                coords[key] = pt
                
    return coords

def generate_mesh(n_sub):
    """
    Generate the topology (list of triangles).
    Each triangle is a tuple of keys (k1, k2, k3).
    """
    triangles = []
    
    for q in range(4):
        for l in range(n_sub):
            # Layer l has (N-l) segments -> (N-l+1) points
            # Layer l+1 has (N-l-1) segments -> (N-l) points
            
            # Number of "base" segments in current layer
            n_seg = n_sub - l
            
            for i in range(n_seg):
                # Points in current layer (Outer)
                p_out_1 = get_point_key(q, l, i, n_sub)
                p_out_2 = get_point_key(q, l, i + 1, n_sub)
                
                # Points in next layer (Inner)
                # Note: Next layer has 1 fewer segment.
                # If we are at the very last layer (l = n_sub - 1), 
                # next layer is just the CORE.
                
                if l == n_sub - 1:
                    # Special case: Tip of the onion, connects to core
                    p_in = get_point_key(0, n_sub, 0, n_sub) # Core
                    triangles.append((p_out_1, p_out_2, p_in))
                else:
                    # Standard strip
                    # To allow decreasing resolution, the topology is:
                    # Outer: i, i+1
                    # Inner: i
                    
                    # However, strictly connecting i to i creates a slant. 
                    # For a generic "Sector" with N segments outer and N-1 inner:
                    # We usually generate:
                    # 1. Triangle (Out[i], Out[i+1], In[i])
                    # 2. Triangle (Out[i+1], In[i+1], In[i]) -> Inverted
                    
                    # Check bounds for inner points
                    # Inner layer has n_seg - 1 segments, so indices 0..(n_seg-1)
                    
                    # Triangle 1 (Base on outer): Always valid
                    # Maps Outer[i]--Outer[i+1] to Inner[i]
                    if i < n_seg:
                         p_in_1 = get_point_key(q, l + 1, i, n_sub)
                         # But wait, does Inner[i] exist?
                         # Layer L+1 has size (N-(L+1)) segments = N-L-1 segments.
                         # Indices 0 to N-L-1.
                         # i goes from 0 to N-L-1.
                         # So yes, valid unless i == n_seg-1 (last segment)
                         
                         if i < (n_seg): 
                             # Wait, logic check. 
                             # Outer: 0..K. Inner: 0..K-1.
                             # If we connect Out[i], Out[i+1], In[i], we cover the whole strip 
                             # BUT we miss the inverted triangles.
                             
                             # Correct strip topology for reduction:
                             # Triangle UP: (Out[i], Out[i+1], In[i])  <-- fills "gaps" near base? No.
                             # Let's visualize:
                             # O0--O1--O2
                             #   \/  \/
                             #   I0--I1
                             
                             # Triangles: (O0, O1, I0), (O1, O2, I1)
                             # Inverted:  (O1, I1, I0)
                             
                             p_in_left = get_point_key(q, l + 1, i, n_sub)
                             
                             # Limit: Inner layer has fewer points.
                             # If i < n_seg (which is max indices of inner?), let's be careful.
                             # n_seg is current layer segments. Inner has n_seg-1 segments.
                             
                             # Case 1: The "V" shape pointing in (O_i, O_i+1, I_i)
                             # This is valid if I_i exists.
                             # Inner indices go 0..(n_seg-1).
                             # So valid for i < n_seg. 
                             # EXCEPT the last outer segment (i = n_seg-1) connects to...
                             # In the standard reduction, the last triangle is just (O_last, O_last+1, I_last)? 
                             # No, indices don't match up 1:1 at the end.
                             
                             # Let's map explicitly:
                             # We need to cover the rectangle defined by Outer[0..N] and Inner[0..N-1].
                             # We generate 2*N - 1 triangles.
                             
                             # 1. "Up" Triangle: (Out[i], Out[i+1], In[i])
                             # Valid for i = 0 to n_seg - 1.
                             # But wait, for the last one (i=n_seg-1), In[i] is the LAST point of inner.
                             # That works.
                             if i < n_seg: 
                                 # We need to handle the "Index out of bounds" for In?
                                 # get_point_key handles wrapping, but here we want strictly the points in this sector.
                                 # Actually, the wrapping handles it. 
                                 # In the visualization O0..O_n connects to I0..I_{n-1}.
                                 pass

                             # Constructing the strip:
                             # For k in 0 to n_seg - 1:
                             #    Tri (Out[k], Out[k+1], In[k])
                             #    If k < n_seg - 1:
                             #        Tri (Out[k+1], In[k+1], In[k])
                             
                             p_in_1 = get_point_key(q, l+1, i, n_sub)
                             
                             # "Pointing In" Triangle
                             triangles.append((p_out_1, p_out_2, p_in_1))
                             
                             # "Pointing Out" Triangle (filling the gap)
                             if i < n_seg - 1:
                                 p_in_2 = get_point_key(q, l+1, i+1, n_sub)
                                 triangles.append((p_out_2, p_in_2, p_in_1))

    return triangles

# ==========================================
# 2. Physics / Weight Logic
# ==========================================

def compute_weights_and_centroids(coords, triangles):
    """
    Compute weights for each point and centroids for each triangle.
    Weight[p] += Area[tri] / 3 for all tris containing p.
    """
    weights = {k: 0.0 for k in coords}
    centroids = [] # List of (cx, cy, [(px,py), ...]) for debug lines
    tri_data = []  # List of (p1, p2, p3, area)
    
    for tri_keys in triangles:
        pts = [coords[k] for k in tri_keys]
        p0, p1, p2 = pts
        
        # Area = 0.5 * |(x1-x0)(y2-y0) - (y1-y0)(x2-x0)|
        area = 0.5 * abs((p1[0]-p0[0])*(p2[1]-p0[1]) - (p1[1]-p0[1])*(p2[0]-p0[0]))
        
        # Centroid
        c = (p0 + p1 + p2) / 3.0
        
        # Accumulate Weights
        w_contrib = area / 3.0
        for k in tri_keys:
            weights[k] += w_contrib
            
        # Store debug info
        # Store lines from Vertex -> Centroid
        lines = [(p, c) for p in pts]
        centroids.append((c, lines, area))
        tri_data.append(pts)
        
    return weights, centroids, tri_data

# ==========================================
# 3. Visualization
# ==========================================

def visualize(n_sub, dmin, dmax, half):
    # Setup
    coords = get_coordinates(n_sub, dmin, dmax, half, np.array([0.0, 0.0]))
    tri_keys = generate_mesh(n_sub)
    weights, debug_data, tris_pts = compute_weights_and_centroids(coords, tri_keys)
    
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # 1. Plot Triangles (Wireframe)
    for pts in tris_pts:
        poly = plt.Polygon(pts, fill=None, edgecolor='k', lw=0.5, alpha=0.3)
        ax.add_patch(poly)

    # 2. Plot "Fair Share" Lines (Vertex to Centroid)
    # This visualizes the Voronoi-dual-like area attribution
    for cent, lines, area in debug_data:
        # Plot lines from vertex to centroid
        for p_v, p_c in lines:
            ax.plot([p_v[0], p_c[0]], [p_v[1], p_c[1]], color='orange', lw=0.5, alpha=0.6)
        
        # Optional: Plot centroid dot
        # ax.scatter([cent[0]], [cent[1]], s=1, color='orange', alpha=0.5)

    # 3. Annotate Points with Weights
    # Filter points to avoid clutter? No, user wants to debug values.
    xs, ys = [], []
    for k, pos in coords.items():
        w = weights[k]
        xs.append(pos[0])
        ys.append(pos[1])
        
        # Check symmetry visually by color
        # Red = Corner-ish, Blue = Edge-ish
        color = 'black'
        
        # Label offset for readability
        ax.annotate(f"{w:.4f}", (pos[0], pos[1]), color='blue', fontsize=12, ha='center', va='center', weight='bold',  bbox=dict(boxstyle='round,pad=0.1', fc='white', alpha=0.7, ec='none'))

    #ax.scatter(xs, ys, s=15, color='black', zorder=10)
    
    # 4. Check Symmetry numerically
    print("--- Symmetry Check ---")
    # Get weights of the 4 main corners
    corners = [
        get_point_key(0, 0, 0, n_sub), # BL (conceptually, depending on start index)
        get_point_key(1, 0, 0, n_sub),
        get_point_key(2, 0, 0, n_sub),
        get_point_key(3, 0, 0, n_sub)
    ]
    corner_weights = [weights[k] for k in corners]
    print(f"Corner Weights (Should be identical): {[float(f'{x:.5f}') for x in corner_weights]}")
    
    # Get weights of the midpoints of the edges
    # Midpoint of layer 0 is index n_sub // 2 ?
    mid_idx = n_sub // 2
    if n_sub % 2 == 0:
        # If even segments, we have an odd point in the middle
        mids = [get_point_key(q, 0, mid_idx, n_sub) for q in range(4)]
        mid_weights = [weights[k] for k in mids]
        print(f"Edge Midpoint Weights (Should be identical): {[float(f'{x:.5f}') for x in mid_weights]}")

    print(f"Core Weight: {weights[('CORE',)]:.5f}")

    ax.set_aspect('equal')
    ax.set_title(f"Triangular Plug (n_sub={n_sub})\nOrange Lines connect Vertices to Triangle Centroids (Visualizing Area Share)")
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=4, help="Subdivisions")
    parser.add_argument("--dmin", type=float, default=0.2, help="Core density factor")
    parser.add_argument("--dmax", type=float, default=1.0, help="Boundary density factor")
    args = parser.parse_args()
    
    visualize(args.n, args.dmin, args.dmax, 1.0)