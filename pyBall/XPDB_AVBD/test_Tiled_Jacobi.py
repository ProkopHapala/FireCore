import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation   import FuncAnimation
from matplotlib.collections import LineCollection
import argparse
import sys
import os

# Ensure we can import modules from the current directory
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from XPDB import XPDB
# Reuse utilities from the original test script
from test_Jacobi_Convergence import (
    make_dense_block, 
    make_random_bonds, 
    _scatter_sizes_from_radius,
    plot_violation_history,
    init_plot,
    update_plot
)

def step_iteration_tiled(sim, args, Rmax, host_pos, iter_count, coll_overlaps, update_plot_fn):
    """One animation frame: rebuild topology and run tiled solver."""
    # Rebuild tiled topology (ghosts and local bond indices)
    sim.build_tiled_topology(Rmax=Rmax, coll_scale=args.coll_scale, bbox_scale=args.bbox_scale)
    
    # Run the tiled solver
    sim.solve_tiled_jacobi(
        dt=args.dt,
        iterations=args.iters_per_frame,
        k_coll=args.k_coll,
        omega=args.omega,
        momentum_beta=args.momentum_beta
    )
    
    # Diagnostics (using global positions for simplicity)
    diag = sim.diagnostics()
    iter_count[0] += args.iters_per_frame
    coll_violation = max(-diag["coll_min"], 0.0)
    bond_violation = max(abs(diag["bond_min"]), abs(diag["bond_max"]))
    coll_overlaps.append((iter_count[0], coll_violation, bond_violation))
    
    # Download for visualization
    sim.get_positions(out=host_pos)

    # Optional: Download debug info
    if args.debug_viz:
        import pyopencl as cl
        f_bond = np.empty((sim.num_atoms, 4), dtype=np.float32)
        f_coll = np.empty((sim.num_atoms, 4), dtype=np.float32)
        cl.enqueue_copy(sim.queue, f_bond, sim.cl_debug_force_bond)
        cl.enqueue_copy(sim.queue, f_coll, sim.cl_debug_force_coll).wait()
        f_bond *= args.force_scale
        f_coll *= args.force_scale
        
        # We can store these in host_pos or return them to update_plot_fn
        # For now let's just use them if update_plot_fn can handle them
        return update_plot_fn(host_pos, f_bond, f_coll)
    
    return update_plot_fn(host_pos)

def update_plot_tiled(scatter, scatter_sizes, bonds, line_segments, force_lines, p, f_bond=None, f_coll=None):
    has_bonds = any(bonds)
    scatter.set_offsets(p[:, :2])
    scatter.set_sizes(scatter_sizes)
    
    if has_bonds:
        segs = []
        for i, blist in enumerate(bonds):
            for j, _, _ in blist:
                if j > i:
                    segs.append([p[i, :2], p[j, :2]])
        line_segments.set_segments(segs)
    else:
        line_segments.set_segments([])

    # Plot forces as lines
    if f_bond is not None and f_coll is not None:
        f_segs = []
        f_colors = []
        scale = 0.1 # Scale for force visualization
        for i in range(len(p)):
            # Bond force (green)
            if np.linalg.norm(f_bond[i, :2]) > 1e-6:
                f_segs.append([p[i, :2], p[i, :2] + f_bond[i, :2] * scale])
                f_colors.append((0, 1, 0, 0.8))
            # Collision force (red)
            if np.linalg.norm(f_coll[i, :2]) > 1e-6:
                f_segs.append([p[i, :2], p[i, :2] + f_coll[i, :2] * scale])
                f_colors.append((1, 0, 0, 0.8))
        force_lines.set_segments(f_segs)
        force_lines.set_colors(f_colors)

    return (scatter, line_segments, force_lines)

def verify_reindexing(sim, bonds, Rmax, coll_scale, bbox_scale):
    """
    Verify that build_local_topology correctly reindexed bonds.
    Compares global bond indices with local bond indices + ghost list.
    """
    # Use current parameters; caller should have uploaded radius already.
    sim.build_tiled_topology(Rmax=Rmax, coll_scale=coll_scale, bbox_scale=bbox_scale)
    
    # Download results
    import pyopencl as cl
    local_bonds = np.empty((sim.num_atoms, 4), dtype=np.int32)
    cl.enqueue_copy(sim.queue, local_bonds, sim.cl_local_bond_indices).wait()
    
    ghost_indices = np.empty(sim.num_groups * sim.max_ghosts, dtype=np.int32)
    cl.enqueue_copy(sim.queue, ghost_indices, sim.cl_ghost_indices).wait()
    
    ghost_counts = np.empty(sim.num_groups, dtype=np.int32)
    cl.enqueue_copy(sim.queue, ghost_counts, sim.cl_ghost_counts).wait()
    
    print("\n--- Verifying Tiled Re-indexing ---")
    errors = 0
    for i in range(sim.num_atoms):
        grp = i // sim.group_size
        g_start = grp * sim.group_size
        g_count = ghost_counts[grp]
        g_list = ghost_indices[grp*sim.max_ghosts : grp*sim.max_ghosts + g_count]
        
        orig_neighs = [b[0] for b in bonds[i]]
        loc_b = local_bonds[i]
        
        for k, ln in enumerate(loc_b):
            if ln == -1: continue
            
            # Resolve local index ln back to global
            if ln < sim.group_size:
                global_target = g_start + ln
            else:
                ghost_idx = ln - sim.group_size
                if ghost_idx >= len(g_list):
                    print(f"Error: Atom {i} has invalid local bond index {ln} (ghost {ghost_idx} out of {len(g_list)})")
                    errors += 1
                    continue
                global_target = g_list[ghost_idx]
                
            if global_target not in orig_neighs:
                print(f"Error: Atom {i} bond {k} points to {global_target}, which is not in original bonds {orig_neighs}")
                errors += 1
    
    if errors == 0:
        print("Success: Tiled re-indexing verified for all atoms.")
    else:
        print(f"Failure: Found {errors} re-indexing errors.")

def print_cluster_mapping(sim, Rmax, coll_scale, bbox_scale):
    """Print per-cluster mapping: internal atoms, ghosts, and bonds with local indices."""
    import pyopencl as cl
    sim.build_tiled_topology(Rmax=Rmax, coll_scale=coll_scale, bbox_scale=bbox_scale)

    # Download buffers
    ghost_indices = np.empty(sim.num_groups * sim.max_ghosts, dtype=np.int32)
    ghost_counts = np.empty(sim.num_groups, dtype=np.int32)
    local_bonds = np.empty((sim.num_atoms, 4), dtype=np.int32)
    cl.enqueue_copy(sim.queue, ghost_indices, sim.cl_ghost_indices).wait()
    cl.enqueue_copy(sim.queue, ghost_counts, sim.cl_ghost_counts).wait()
    cl.enqueue_copy(sim.queue, local_bonds, sim.cl_local_bond_indices).wait()

    for grp in range(sim.num_groups):
        g_start = grp * sim.group_size
        g_end = min(g_start + sim.group_size, sim.num_atoms)
        g_count = ghost_counts[grp]
        g_list = ghost_indices[grp*sim.max_ghosts : grp*sim.max_ghosts + g_count]
        print(f"\nCluster {grp}: internal [{g_start}:{g_end}), ghosts {g_count}")
        print("  Internal (local->global):", [f"{li}->{g_start+li}" for li in range(g_end - g_start)])
        print("  Ghosts   (local->global):", [f"{sim.group_size+gi}->{g_list[gi]}" for gi in range(g_count)])
        print("  Bonds (global: local targets):")
        for gid in range(g_start, g_end):
            lbs = local_bonds[gid]
            mapped = []
            for ln in lbs:
                if ln == -1:
                    mapped.append("-1")
                elif ln < sim.group_size:
                    mapped.append(f"{ln}(int)")
                else:
                    mapped.append(f"{ln}(ghost)")
            print(f"    {gid}: {mapped}")

def report_overstretched_bonds(sim, bonds, Rmax, coll_scale, bbox_scale, tol):
    """
    After relaxation, report bonds whose length exceeds (1+tol) * rest_len.
    Also prints local mapping for each offending bond.
    """
    import pyopencl as cl
    # Rebuild topology to get current ghosts/local indices
    sim.build_tiled_topology(Rmax=Rmax, coll_scale=coll_scale, bbox_scale=bbox_scale)

    # Download positions and mapping
    pos = sim.get_positions(out=None)  # (num_atoms,3)
    ghost_indices = np.empty(sim.num_groups * sim.max_ghosts, dtype=np.int32)
    ghost_counts = np.empty(sim.num_groups, dtype=np.int32)
    local_bonds = np.empty((sim.num_atoms, 4), dtype=np.int32)
    cl.enqueue_copy(sim.queue, ghost_indices, sim.cl_ghost_indices).wait()
    cl.enqueue_copy(sim.queue, ghost_counts, sim.cl_ghost_counts).wait()
    cl.enqueue_copy(sim.queue, local_bonds, sim.cl_local_bond_indices).wait()

    print("\n--- Overstretched Bonds Report ---")
    count = 0
    for i, b_list in enumerate(bonds):
        grp = i // sim.group_size
        g_start = grp * sim.group_size
        g_count = ghost_counts[grp]
        g_list = ghost_indices[grp*sim.max_ghosts : grp*sim.max_ghosts + g_count]

        for slot, (j, rest, k_b) in enumerate(b_list[:4]):
            p_i = pos[i]; p_j = pos[j]
            dist = float(np.linalg.norm(p_i - p_j))
            if dist > (1.0 + tol) * rest:
                # Find local index used in mapping
                ln = local_bonds[i, slot] if slot < 4 else -1
                loc_str = ""
                if ln == -1:
                    loc_str = "absent"
                elif ln < sim.group_size:
                    loc_str = f"{ln}(int)"
                else:
                    gidx = ln - sim.group_size
                    if 0 <= gidx < len(g_list):
                        loc_str = f"{ln}(ghost->{g_list[gidx]})"
                    else:
                        loc_str = f"{ln}(ghost?out-of-range)"
                print(f"Bond {i}->{j} slot{slot}: dist={dist:.4f} rest={rest:.4f} over={(dist/rest-1)*100:.1f}% "
                      f"[grp {grp}, local {loc_str}]")
                count += 1
    if count == 0:
        print("No overstretched bonds beyond tolerance.")

def run_tiled_animation():
    parser = argparse.ArgumentParser(description="Tiled Jacobi Solver Test")
    parser.add_argument("--nx",              type=int,   default=10)
    parser.add_argument("--ny",              type=int,   default=10)
    parser.add_argument("--spacing",         type=float, default=0.35)
    parser.add_argument("--jitter",          type=float, default=0.1)
    parser.add_argument("--atom_rad",        type=float, default=0.2)
    parser.add_argument("--atom_rad_max",    type=float, default=None, help="If set, draw radii uniformly in [atom_rad, atom_rad_max]")
    parser.add_argument("--omega",           type=float, default=0.5)
    parser.add_argument("--momentum_beta",   type=float, default=0.7)
    parser.add_argument("--dt",              type=float, default=1.0)
    parser.add_argument("--seed",            type=int,   default=587)
    parser.add_argument("--bond_prob",       type=float, default=0.25)
    parser.add_argument("--k_coll",          type=float, default=100.0)
    parser.add_argument("--bond_k",          type=float, default=200.0)
    parser.add_argument("--bond_len",        type=float, default=0.4)
    parser.add_argument("--iters_per_frame", type=int,   default=1, help="Jacobi iterations per animation frame")
    parser.add_argument("--frames",          type=int,   default=100, help="Number of animation frames")
    parser.add_argument("--coll_scale",      type=float, default=2.0, help="Collision search margin multiplier of Rmax")
    parser.add_argument("--bbox_scale",      type=float, default=2.0, help="BBox padding multiplier of Rmax")
    parser.add_argument("--group_size",      type=int,   default=64, help="Workgroup size (tiles)")
    parser.add_argument("--max_ghosts",      type=int,   default=128, help="Max ghosts per cluster")
    parser.add_argument("--debug_viz",       type=int,   default=1, help="Visualize debug forces")
    parser.add_argument("--force_scale",     type=float, default=0.1, help="Scale factor applied to debug forces when plotting")
    parser.add_argument("--print_mapping",   type=int,   default=0, help="Print cluster->local/global mapping after topology build")
    parser.add_argument("--label_atoms",     type=int,   default=1, help="Draw global atom indices on the plot")
    parser.add_argument("--stretch_tol",     type=float, default=0.1, help="Relative tolerance for bond overstretch reporting (e.g., 0.1 = 10%)")
    args = parser.parse_args()

    np.random.seed(args.seed)
    pos0, vel0, radius, mass = make_dense_block(args.nx, args.ny, args.spacing, args.jitter, args.atom_rad)
    if args.atom_rad_max is not None:
        r_hi = max(args.atom_rad, args.atom_rad_max)
        r_lo = min(args.atom_rad, args.atom_rad_max)
        radius = np.random.uniform(r_lo, r_hi, size=radius.shape).astype(np.float32)
    Rmax = float(np.max(radius))
    num_atoms = len(pos0)
    
    # XPDB setup (using WG_SIZE=64 for NVIDIA optimization as discussed)
    sim = XPDB(num_atoms, group_size=args.group_size)
    sim.upload_data(pos0, vel0, radius, mass)
    bonds = make_random_bonds(pos0, prob=args.bond_prob, rest_len=args.bond_len, k_bond=args.bond_k)
    
    # Initialize tiled structures
    sim.init_tiled(max_ghosts=args.max_ghosts)
    sim.upload_global_bonds(bonds)
    
    # Upload bond params for solver (standard CSR used for diagnostics and potentially solver)
    if any(bonds):
        sim.upload_bonds(bonds)
    else:
        sim.set_empty_bonds()

    # 1. Verification Step
    verify_reindexing(sim, bonds, Rmax=Rmax, coll_scale=args.coll_scale, bbox_scale=args.bbox_scale)
    if args.print_mapping:
        print_cluster_mapping(sim, Rmax=Rmax, coll_scale=args.coll_scale, bbox_scale=args.bbox_scale)

    # 2. Animation Step
    fig, ax = plt.subplots(figsize=(8, 8))
    margin = 1.0
    ax.set_xlim(np.min(pos0[:, 0]) - margin, np.max(pos0[:, 0]) + margin)
    ax.set_ylim(np.min(pos0[:, 1]) - margin, np.max(pos0[:, 1]) + margin)
    ax.set_aspect('equal')
    fig.canvas.draw()
    
    scatter_sizes = _scatter_sizes_from_radius(ax, fig, radius)
    cmap = plt.cm.get_cmap("tab20")
    cluster_colors = cmap((np.arange(sim.num_groups) % cmap.N) / cmap.N)
    atom_colors = cluster_colors[np.arange(num_atoms) // sim.group_size]
    scatter = ax.scatter(pos0[:, 0], pos0[:, 1], s=scatter_sizes, c=atom_colors, alpha=0.7, edgecolors='k', linewidths=0.3)
    line_segments = LineCollection([], colors='gray', linewidths=0.5)
    ax.add_collection(line_segments)
    
    force_lines = LineCollection([], linewidths=1.0)
    ax.add_collection(force_lines)

    # Optional labels
    text_labels = []
    if args.label_atoms:
        for i in range(num_atoms):
            t = ax.text(pos0[i, 0], pos0[i, 1], str(i), fontsize=7, ha='center', va='center')
            text_labels.append(t)
    
    # Setup iteration state
    iter_count = [0]
    coll_overlaps = []
    host_pos = np.empty((num_atoms, 4), dtype=np.float32)
    
    # Seed iteration buffers (A, prev) for momentum
    sim.reset_iteration_buffers(src_buffer=sim.cl_pos, reset_pred=True)

    def update_plot_fn(p, fb=None, fc=None):
        artists = update_plot_tiled(scatter, scatter_sizes, bonds, line_segments, force_lines, p, fb, fc)
        if args.label_atoms:
            for i, t in enumerate(text_labels):
                t.set_position((p[i, 0], p[i, 1]))
            artists = artists + tuple(text_labels)
        return artists

    step_fn = lambda frame: step_iteration_tiled(sim, args, Rmax, host_pos, iter_count, coll_overlaps, update_plot_fn)
    init_fn = lambda: init_plot(scatter, scatter_sizes, pos0, bonds, line_segments)

    ani = FuncAnimation(fig, step_fn, frames=args.frames, init_func=init_fn, interval=30, blit=False, repeat=False)
    plt.show()

    if coll_overlaps:
        plot_violation_history(coll_overlaps)
        plt.show()

    # Post-run bond stretch report
    report_overstretched_bonds(sim, bonds, Rmax=Rmax, coll_scale=args.coll_scale, bbox_scale=args.bbox_scale, tol=args.stretch_tol)

if __name__ == "__main__":
    run_tiled_animation()
