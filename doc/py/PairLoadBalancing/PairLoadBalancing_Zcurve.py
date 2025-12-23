
"""
This script is designed as a didactic "textbook" lecture on solving the two most difficult problems in GPU-accelerated molecular simulations: **Memory Bandwidth Bottlenecks** and **Workload Imbalance.**

### The Problem Summary
1.  **Memory Bottleneck:** Reading atom data from Global VRAM is 100x slower than math. To speed up the code, a GPU "work-group" must load a set of atoms into **Local Memory (LDS)** once and reuse them for many calculations. If interactions are scattered randomly, we reload the same data constantly, killing performance.
2.  **Load Balancing:** Molecules are uneven. Some regions are dense (3D clusters), others are sparse (1D chains). If we assign "10 atoms per work-group," the "dense" work-groups will take 100x longer to finish than the "sparse" ones.

### The Reasoning for the Solution
*   **Spatial Sorting (Z-Curve):** We map 3D space to a 1D line. By sorting atoms along this line, atoms that are close in space become close in the array index. This clusters the "1s" in our Density Matrix near the diagonal.
*   **Linearized Pair-Lists:** Instead of partitioning by **atoms**, we partition by **interactions (pairs)**. By taking a fixed-size "chunk" of pairs from the sorted list, we guarantee every GPU work-group does exactly the same amount of math (Perfect Load Balance), while the Z-order ensures they reuse the same atoms (Efficient Memory Use).

### Critical Observations for your integration:

1.  **Midpoint Strategy:** For the Density Matrix, sorting by `i` (the first atom) is usually enough. But for **Grid Projection**, I implemented the **Midpoint Z-order** in the script. This ensures that the work-groups processing "bond charges" are grouped by where that bond exists in 3D space, which is critical for caching the grid voxels.
2.  **The LDS Limit:** Notice the "Avg unique atoms per WG" in the printout. This number tells you how many `float4` atomic positions and `float4` orbital coefficients you need to fit into the GPU's **Local Memory (LDS)**. At ~30 atoms, you are well within the 32KB-64KB limit of modern GPUs.
3.  **Load Balancing:** Because we partition the **array of pairs** (the `pairs_i` and `pairs_j` arrays), every work-group is active for the same number of clock cycles. There is no "tail" of work where one thread is waiting for another.


### Python Implementation: Z-Curve Analysis for Fireball
"""

import argparse
import math
import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# LECTURE PART 1: The Z-Curve (Morton Order)
# Why? We need to turn 3D coordinates into 1D indices while keeping 
# neighboring atoms close together in the array.
# =============================================================================

def interleave_3d(x, y, z):
    """
    Spreads bits of x, y, z so they can be interleaved.
    Interleaving [x2 x1 x0], [y2 y1 y0], [z2 z1 z0] -> [x2 y2 z2 x1 y1 z1 x0 y0 z0]
    """
    def part_1d(n):
        n = int(n) # ensure Python int to avoid numpy overflow
        n &= 0x000003ff # Use 10 bits (total 30 bits)
        n = (n | (n << 16)) & 0xff0000ff
        n = (n | (n << 8))  & 0x0300f00f
        n = (n | (n << 4))  & 0x030c30c3
        n = (n | (n << 2))  & 0x09249249
        return n
    return (part_1d(x) << 2) | (part_1d(y) << 1) | part_1d(z)

def hilbert_index_3d(x, y, z, bits=10):
    """Simple 3D Hilbert index (Skilling-inspired). bits <= 16 recommended."""
    x, y, z = int(x), int(y), int(z)
    n_side = 1 << bits
    mask = 1 << (bits - 1)
    xi, yi, zi = x, y, z
    index = 0
    for _ in range(bits):
        rx = 1 if (xi & mask) else 0
        ry = 1 if (yi & mask) else 0
        rz = 1 if (zi & mask) else 0
        digit = (rx << 2) | (ry << 1) | rz
        index = (index << 3) | digit
        # rotate/reflection following 3D Hilbert rules
        if ry == 0:
            if rz == 1:
                xi = (n_side - 1) ^ xi
                zi = (n_side - 1) ^ zi
            xi, zi = zi, xi
        mask >>= 1
    return index

def get_indices(pts, res=1024, curve="morton", bits=10):
    """ Maps 3D points to 1D curve integers (morton or hilbert). """
    p_min, p_max = pts.min(axis=0), pts.max(axis=0)
    q = ((pts - p_min) / (p_max - p_min + 1e-9) * (res - 1)).astype(np.int32)
    if curve == "morton":
        return np.array([interleave_3d(p[0], p[1], p[2]) for p in q])
    if curve == "hilbert":
        return np.array([hilbert_index_3d(p[0], p[1], p[2], bits=bits) for p in q])
    raise ValueError(f"Unknown curve {curve}")

# alias for backward compatibility in comments
get_z_indices = get_indices

# =============================================================================
# LECTURE PART 2: Data Generation
# Two strategies:
#  - randomwalk: polymer-like folded chain
#  - box: uniform random points in a cube
# =============================================================================

def make_atoms(strategy, n_atoms, box, walk_sigma=0.5):
    if strategy == "randomwalk":
        steps = np.random.normal(0, walk_sigma, (n_atoms, 3))
        return np.cumsum(steps, axis=0)
    if strategy == "box":
        low, high = -box*0.5, box*0.5
        return np.random.uniform(low, high, (n_atoms, 3))
    raise ValueError(f"Unknown strategy {strategy}")

def build_system(args):
    np.random.seed(args.seed)
    pts = make_atoms(args.strategy, args.n_atoms, args.box, args.walk_sigma)
    z_codes = get_indices(pts, res=args.res, curve=args.curve, bits=args.bits)
    z_sort_idx = np.argsort(z_codes)
    pts_sorted = pts[z_sort_idx]
    return pts_sorted

def get_indices_2d(i, j):
    # Simple Morton (Z-curve) interleave for 2D indices (uses 16 bits safe range)
    def part_1d(n):
        n = n & 0x0000ffff
        n = (n | (n << 8)) & 0x00ff00ff
        n = (n | (n << 4)) & 0x0f0f0f0f
        n = (n | (n << 2)) & 0x33333333
        n = (n | (n << 1)) & 0x55555555
        return n
    return (part_1d(i) << 1) | part_1d(j)

def sort_pairs(pairs_i, pairs_j, pts, args):
    if args.pair_sort == "midpoint":
        midpoints = (pts[pairs_i] + pts[pairs_j]) * 0.5
        sort_val = get_indices(midpoints, res=args.res, curve=args.curve, bits=args.bits)
    elif args.pair_sort == "lex":
        sort_val = pairs_i * (pts.shape[0] + 1) + pairs_j
    elif args.pair_sort == "hilbert2d":
        sort_val = get_indices_2d(pairs_i, pairs_j)
    else:  # none
        return pairs_i, pairs_j
    idx = np.argsort(sort_val)
    return pairs_i[idx], pairs_j[idx]

def find_pairs(pts, rcut, res, curve, bits, args):
    diff = pts[:, np.newaxis, :] - pts[np.newaxis, :, :]
    dist = np.sqrt(np.sum(diff**2, axis=-1))
    adj = dist < rcut
    pairs_i, pairs_j = np.where(np.triu(adj, k=1))
    pairs_i, pairs_j = sort_pairs(pairs_i, pairs_j, pts, args)
    return adj, pairs_i, pairs_j

def analyze_chunks(pairs_i, pairs_j, chunk_size, n_atoms):
    n_pairs = len(pairs_i)
    n_chunks = n_pairs // chunk_size
    stats_unique_atoms = []
    workload_map = np.zeros((n_atoms, n_atoms))
    for c in range(n_chunks):
        idx = slice(c*chunk_size, (c+1)*chunk_size)
        ci, cj = pairs_i[idx], pairs_j[idx]
        unique_atoms = len(np.unique(np.concatenate([ci, cj])))
        stats_unique_atoms.append(unique_atoms)
        workload_map[ci, cj] = (c % 10) + 1
    return n_pairs, n_chunks, np.array(stats_unique_atoms), workload_map

def neighbor_stats(adj):
    adj_no_self = adj.copy()
    np.fill_diagonal(adj_no_self, False)
    counts = np.sum(adj_no_self, axis=1)
    return counts.min(), counts.mean(), counts.max()

def plot_all(pts, adj, workload_map, n_atoms, chunk_size, plots):
    # plots: list of tokens from CLI, e.g. ["spatial","workload","adj","xyz"]
    plot_funcs = []
    if "spatial" in plots:
        plot_funcs.append("spatial")
    if "adj" in plots:
        plot_funcs.append("adj")
    if "workload" in plots:
        plot_funcs.append("workload")
    if "xyz" in plots:
        plot_funcs.append("xyz")
    n = len(plot_funcs)
    if n == 0:
        return
    cols = 2 if n > 1 else 1
    rows = (n + cols - 1) // cols
    fig = plt.figure(figsize=(8*cols, 5*rows))
    for idx, key in enumerate(plot_funcs):
        ax = fig.add_subplot(rows, cols, idx+1, projection='3d' if key == "spatial" else None)
        if key == "spatial":
            ax.plot(pts[:,0], pts[:,1], pts[:,2], c='gray', alpha=0.3)
            ax.scatter(pts[:,0], pts[:,1], pts[:,2], c=np.arange(n_atoms), cmap='jet')
            ax.set_title("Atoms Sorted by Z-Curve\n(Color = Array Index)")
        elif key == "adj":
            ax.imshow(adj, cmap='Greys', interpolation='nearest')
            ax.set_title("Adjacency Matrix\n(Z-sorting concentrates pairs near diagonal)")
        elif key == "workload":
            cmap_chunks = plt.colormaps.get_cmap('tab10').copy()
            cmap_chunks.set_bad('white')
            workload_map_masked = np.ma.masked_equal(workload_map, 0)
            ax.imshow(workload_map_masked, cmap=cmap_chunks, interpolation='nearest', vmin=1, vmax=10)
            ax.set_title(f"GPU Work-groups (Chunks of {chunk_size})\n(Every color = Equal Math load)")
        elif key == "xyz":
            t = np.arange(n_atoms)
            ax.plot(t, pts[:,0], label='x')
            ax.plot(t, pts[:,1], label='y')
            ax.plot(t, pts[:,2], label='z')
            ax.set_xlabel("Atom index (Z-ordered)")
            ax.set_ylabel("Coordinate")
            ax.set_title("Coordinates along Z-curve ordering")
            ax.legend()
    plt.tight_layout()

def print_stats(n_atoms, n_pairs, chunk_size, stats_unique_atoms, neigh_min, neigh_avg, neigh_max):
    avg_atoms = np.mean(stats_unique_atoms)
    min_atoms = np.min(stats_unique_atoms)
    max_atoms = np.max(stats_unique_atoms)
    efficiency = 1.0 - (avg_atoms / (2 * chunk_size))
    fill_factor = n_pairs / float(n_atoms * n_atoms)
    print(f"--- RESULTS ---")
    print(f"Total Atoms:             {n_atoms}")
    print(f"Total Pairs to calculate: {n_pairs}")
    print(f"Pairs per Work-group:    {chunk_size} (Perfectly Balanced)")
    print(f"Unique atoms per WG (min/avg/max): {min_atoms:.0f}/{avg_atoms:.1f}/{max_atoms:.0f} (out of {2*chunk_size} possible)")
    print(f"Neighbors per atom (min/avg/max): {neigh_min:.0f}/{neigh_avg:.1f}/{neigh_max:.0f}")
    print(f"Fill factor (pairs / N^2): {fill_factor:.4f}")
    print(f"Memory Reuse Efficiency: {efficiency:.1%}")
    print(f"\nREASONING:")
    print(f"1. Every GPU Work-group calculates exactly {chunk_size} interactions.")
    print(f"2. Because of Z-sorting, those {chunk_size} interactions only involve ~{avg_atoms:.1f} atoms.")
    print(f"3. We load {avg_atoms:.1f} atoms into Local Memory and perform {chunk_size} x (Basis Evals).")
    print(f"4. This minimizes Global Memory access, which is the main bottleneck.")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n_atoms",    type=int,   default=400, help="Number of atoms")
    ap.add_argument("--rcut",       type=float, default=2.0, help="Cutoff distance")
    ap.add_argument("--chunk",      type=int,   default=32, help="Pairs per work-group")
    ap.add_argument("--seed",       type=int,   default=42, help="RNG seed")
    ap.add_argument("--strategy",   choices=["randomwalk", "box"], default="box", help="Atom placement strategy")
    ap.add_argument("--box",        type=float, default=10.0, help="Box size for uniform strategy")
    ap.add_argument("--walk_sigma", type=float, default=0.5, help="Step stddev for randomwalk")
    ap.add_argument("--curve",      choices=["morton", "hilbert"], default="hilbert", help="Space-filling curve for sorting")
    ap.add_argument("--res",        type=int,   default=1024, help="Quantization resolution per axis (points snapped to [0,res-1])")
    ap.add_argument("--bits",       type=int,   default=10, help="Bits per axis (must match res<=2^bits)")
    ap.add_argument("--plots",      type=str,   default="spatial,workload,xyz", help="Comma-separated plots: spatial, workload, adj, xyz")
    ap.add_argument("--pair_sort",  choices=["midpoint", "lex", "hilbert2d", "none"], default="midpoint", help="How to order pairs before chunking")
    args = ap.parse_args()
    pts = build_system(args)
    adj, pairs_i, pairs_j = find_pairs(pts, args.rcut, args.res, args.curve, args.bits, args)
    neigh_min, neigh_avg, neigh_max = neighbor_stats(adj)
    n_pairs, n_chunks, stats_unique_atoms, workload_map = analyze_chunks(pairs_i, pairs_j, args.chunk, len(pts))
    plots = [p.strip() for p in args.plots.split(",") if p.strip()]
    plot_all(pts, adj, workload_map, len(pts), args.chunk, plots)
    print_stats(len(pts), n_pairs, args.chunk, stats_unique_atoms, neigh_min, neigh_avg, neigh_max)
    plt.show()