import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.collections import LineCollection
import argparse
import sys
import os
import itertools
import pyopencl as cl

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from XPDB_new import XPDB_new


def make_dense_block(nx, ny, spacing, jitter, atom_rad):
    """Create a dense block of atoms"""
    pos = []
    for i in range(nx):
        for j in range(ny):
            x = (i - nx/2) * spacing + np.random.uniform(-jitter, jitter)
            y = (j - ny/2) * spacing + np.random.uniform(-jitter, jitter)
            pos.append([x, y, 0.0])
    pos = np.array(pos, dtype=np.float32)
    vel = np.zeros_like(pos)
    radius = np.full(len(pos), atom_rad, dtype=np.float32)
    mass = np.ones(len(pos), dtype=np.float32)
    return pos, vel, radius, mass


def make_random_bonds(pos, prob=0.25, rest_len=0.4, k_bond=200.0):
    """Create random bonds between nearby atoms"""
    num_atoms = len(pos)
    bonds_adj = [[] for _ in range(num_atoms)]
    
    for i in range(num_atoms):
        for j in range(i+1, num_atoms):
            dist = np.linalg.norm(pos[i] - pos[j])
            if dist < rest_len * 1.5 and np.random.random() < prob:
                bonds_adj[i].append([j, rest_len, k_bond])
                bonds_adj[j].append([i, rest_len, k_bond])
    
    return bonds_adj


def update_plot(scatter, scatter_sizes, bonds, line_segments, force_lines, p, f_bond=None, f_coll=None):
    """Update plot for animation"""
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
        scale = 0.1
        for i in range(len(p)):
            if np.linalg.norm(f_bond[i, :2]) > 1e-6:
                f_segs.append([p[i, :2], p[i, :2] + f_bond[i, :2] * scale])
                f_colors.append((0, 1, 0, 0.8))
            if np.linalg.norm(f_coll[i, :2]) > 1e-6:
                f_segs.append([p[i, :2], p[i, :2] + f_coll[i, :2] * scale])
                f_colors.append((1, 0, 0, 0.8))
        force_lines.set_segments(f_segs)
        force_lines.set_colors(f_colors)

    return (scatter, line_segments, force_lines)


def run_test():
    parser = argparse.ArgumentParser(description="Tiled Jacobi Solver Test (Simplified)")
    parser.add_argument("--nx", type=int, default=10)
    parser.add_argument("--ny", type=int, default=10)
    parser.add_argument("--spacing", type=float, default=0.35)
    parser.add_argument("--jitter", type=float, default=0.1)
    parser.add_argument("--atom_rad", type=float, default=0.2)
    parser.add_argument("--omega", type=float, default=0.5)
    parser.add_argument("--momentum_beta", type=float, default=0.7)
    parser.add_argument("--dt", type=float, default=1.0)
    parser.add_argument("--seed", type=int, default=587)
    parser.add_argument("--bond_prob", type=float, default=0.25)
    parser.add_argument("--k_coll", type=float, default=100.0)
    parser.add_argument("--bond_k", type=float, default=200.0)
    parser.add_argument("--bond_len", type=float, default=0.4)
    parser.add_argument("--iters_per_frame", type=int, default=1)
    parser.add_argument("--frames", type=int, default=-1, help="Number of frames; <0 to run until window close")
    parser.add_argument("--coll_scale", type=float, default=2.0)
    parser.add_argument("--bbox_scale", type=float, default=2.0)
    parser.add_argument("--group_size", type=int, default=64)
    parser.add_argument("--max_ghosts", type=int, default=128)
    parser.add_argument("--debug_viz", type=int, default=1)
    parser.add_argument("--force_scale", type=float, default=0.1)
    parser.add_argument("--pick_radius", type=float, default=0.4, help="Max distance for picking (world units)")
    args = parser.parse_args()

    np.random.seed(args.seed)
    pos0, vel0, radius, mass = make_dense_block(args.nx, args.ny, args.spacing, args.jitter, args.atom_rad)
    Rmax = float(np.max(radius))
    num_atoms = len(pos0)
    
    # Initialize XPDB_new
    sim = XPDB_new(num_atoms, group_size=args.group_size)
    sim.upload_data(pos0, vel0, radius, mass)
    
    bonds = make_random_bonds(pos0, prob=args.bond_prob, rest_len=args.bond_len, k_bond=args.bond_k)
    sim.upload_bonds_fixed(bonds)
    
    # Setup plot
    fig, ax = plt.subplots(figsize=(8, 8))
    margin = 1.0
    ax.set_xlim(np.min(pos0[:, 0]) - margin, np.max(pos0[:, 0]) + margin)
    ax.set_ylim(np.min(pos0[:, 1]) - margin, np.max(pos0[:, 1]) + margin)
    ax.set_aspect('equal')
    fig.canvas.draw()
    
    scatter_sizes = (radius * 1000).astype(int)
    cmap = plt.cm.get_cmap("tab20")
    cluster_colors = cmap((np.arange(sim.num_groups) % cmap.N) / cmap.N)
    atom_colors = cluster_colors[np.arange(num_atoms) // sim.group_size]
    scatter = ax.scatter(pos0[:, 0], pos0[:, 1], s=scatter_sizes, c=atom_colors, alpha=0.7, edgecolors='k', linewidths=0.3)
    line_segments = LineCollection([], colors='gray', linewidths=0.5)
    ax.add_collection(line_segments)
    
    force_lines = LineCollection([], linewidths=1.0)
    ax.add_collection(force_lines)

    host_pos = np.empty((num_atoms, 4), dtype=np.float32)
    host_pred = np.empty((num_atoms, 4), dtype=np.float32)
    # Seed host_pos with initial positions for picking
    sim.get_positions(out=host_pos)

    pick_state = {
        "idx": None,
        "active": False,
        "mouse": np.array([0.0, 0.0], dtype=np.float32),
    }

    def on_press(event):
        if event.inaxes != ax or event.xdata is None or event.ydata is None:
            return
        mouse_xy = np.array([event.xdata, event.ydata], dtype=np.float32)
        # Compute distance to all atoms (2D)
        d2 = np.sum((host_pos[:num_atoms, :2] - mouse_xy) ** 2, axis=1)
        i_min = int(np.argmin(d2))
        if d2[i_min] <= args.pick_radius * args.pick_radius:
            pick_state["idx"] = i_min
            pick_state["active"] = True
            pick_state["mouse"] = mouse_xy
            print(f"[DEBUG] on_press hit: idx={i_min} d={np.sqrt(d2[i_min]):.4f} radius={args.pick_radius} mouse={mouse_xy}")
        else:
            print(f"[DEBUG] on_press miss: closest idx={i_min} d={np.sqrt(d2[i_min]):.4f} radius={args.pick_radius} mouse={mouse_xy}")
            print(f"[DEBUG] pos closest={host_pos[i_min, :2]}")

    def on_release(event):
        if pick_state["active"]:
            print(f"[DEBUG] on_release: idx={pick_state['idx']}")
        pick_state["active"] = False
        pick_state["idx"] = None

    def on_motion(event):
        if not pick_state["active"]:
            return
        if event.xdata is None or event.ydata is None:
            return
        pick_state["mouse"] = np.array([event.xdata, event.ydata], dtype=np.float32)
        print(f"[DEBUG] on_motion: idx={pick_state['idx']} mouse={pick_state['mouse']}")

    cid_press = fig.canvas.mpl_connect('button_press_event', on_press)
    cid_release = fig.canvas.mpl_connect('button_release_event', on_release)
    cid_motion = fig.canvas.mpl_connect('motion_notify_event', on_motion)
    
    def update(frame):
        # Build pred targets from current positions, override picked atom to mouse
        host_pred[:, :3] = host_pos[:, :3]
        host_pred[:, 3] = 0.0
        if pick_state["active"] and pick_state["idx"] is not None:
            idx = pick_state["idx"]
            host_pred[idx, 0] = pick_state["mouse"][0]
            host_pred[idx, 1] = pick_state["mouse"][1]
            host_pred[idx, 2] = 0.0
            # Force the current position to the mouse location for a strong drag
            host_pos[idx, 0] = pick_state["mouse"][0]
            host_pos[idx, 1] = pick_state["mouse"][1]
            host_pos[idx, 2] = 0.0
            cl.enqueue_copy(sim.queue, sim.cl_pos, host_pos)
        # Push pred targets to device
        cl.enqueue_copy(sim.queue, sim.cl_pred_pos, host_pred)

        # Build topology
        sim.build_local_topology(Rmax=Rmax, coll_scale=args.coll_scale, bbox_scale=args.bbox_scale)

        # Solve
        sim.solve_cluster_jacobi(
            dt=args.dt,
            iterations=args.iters_per_frame,
            k_coll=args.k_coll,
            omega=args.omega,
            momentum_beta=args.momentum_beta
        )
        
        # Download positions
        sim.get_positions(out=host_pos)
        
        # Download forces if debug viz enabled
        if args.debug_viz:
            f_bond, f_coll = sim.get_debug_forces()
            f_bond *= args.force_scale
            f_coll *= args.force_scale
            return update_plot(scatter, scatter_sizes, bonds, line_segments, force_lines, host_pos, f_bond, f_coll)
        
        return update_plot(scatter, scatter_sizes, bonds, line_segments, force_lines, host_pos)
    
    frames_iter = itertools.count() if args.frames < 0 else range(args.frames)
    ani = FuncAnimation(fig, update, frames=frames_iter, interval=30, blit=False, repeat=False)
    plt.show()


if __name__ == "__main__":
    run_test()
