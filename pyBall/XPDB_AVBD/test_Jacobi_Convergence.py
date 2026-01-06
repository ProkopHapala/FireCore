import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.collections import LineCollection
import argparse
import sys
import os
# Silence PyOpenCL compiler output warnings unless explicitly enabled
os.environ.setdefault("PYOPENCL_COMPILER_OUTPUT", "0")
import pyopencl as cl
import warnings
warnings.filterwarnings("ignore", category=cl.CompilerWarning)

# Ensure we can import XPDB from the current directory
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from XPDB import XPDB

def _scatter_sizes_from_radius(ax, fig, radius):
    """Compute scatter areas (points^2) so marker radius matches data-space radius.
    Accepts scalar or array radius; returns array matching input."""
    p0 = ax.transData.transform((0, 0))
    p1 = ax.transData.transform((1, 1))
    sx = abs(p1[0] - p0[0])
    sy = abs(p1[1] - p0[1])
    sp = min(sx, sy)
    dpi = fig.dpi
    pr_points = np.asarray(radius) * sp * 72.0 / dpi
    return np.pi * (pr_points**2)

def plot_violation_series(results, num_iters, title_suffix=""):
    """Plot collision and bond violation series for multiple momentum betas."""
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharex=True)
    axc, axb = axes
    iters = np.arange(1, num_iters + 1)
    for label, cs, bs in results:
        axc.semilogy(iters, np.maximum(cs, 1e-12), label=str(label))
        axb.semilogy(iters, np.maximum(bs, 1e-12), label=str(label))
    axc.set_title(f"Collision violation{title_suffix}")
    axb.set_title(f"Bond violation{title_suffix}")
    for ax in axes:
        ax.set_xlabel("Iteration")
        ax.grid(True, alpha=0.3)
    axc.set_ylabel("Violation magnitude")
    axc.legend(title="beta")
    plt.tight_layout()
    return fig, axes

def plot_violation_history(coll_overlaps):
    """Plot violation history from animation run (iteration, coll, bond)."""
    data = np.array(coll_overlaps, dtype=np.float32)
    iters = data[:, 0]
    coll_viol = data[:, 1]
    bond_viol = data[:, 2]

    plt.figure()
    plt.semilogy(iters, np.maximum(coll_viol, 1e-12), label="collision |min|")
    plt.semilogy(iters, np.maximum(bond_viol, 1e-12), label="bond |max|")
    plt.xlabel("Iteration")
    plt.ylabel("Violation magnitude")
    plt.title("Constraint violation vs iterations (log scale)")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

def prepare_sim(sim, pos, vel, radius, mass, bonds, neighbor_margin=1.0):
    sim.upload_data(pos, vel, radius, mass)
    if any(bonds):
        sim.upload_bonds(bonds)
    else:
        sim.set_empty_bonds()
    sim.ensure_neighbors(margin=neighbor_margin)
    sim.reset_iteration_buffers(src_buffer=sim.cl_pos, reset_pred=True)

def make_dense_block(nx, ny, spacing, jitter, radius_val):
    """Create a dense block of particles with small random jitter."""
    xs = np.linspace(0, (nx-1)*spacing, nx)
    ys = np.linspace(0, (ny-1)*spacing, ny)
    X, Y = np.meshgrid(xs, ys)
    X = X.flatten()
    Y = Y.flatten()
    
    num_atoms = len(X)
    pos = np.zeros((num_atoms, 3), dtype=np.float32)
    pos[:, 0] = X + np.random.uniform(-jitter, jitter, num_atoms)
    pos[:, 1] = Y + np.random.uniform(-jitter, jitter, num_atoms)
    
    radius = np.full(num_atoms, radius_val, dtype=np.float32)
    mass = np.ones(num_atoms, dtype=np.float32)
    vel = np.zeros((num_atoms, 3), dtype=np.float32)
    
    return pos, vel, radius, mass

def make_random_bonds(pos, prob, rest_len, k_bond, max_neighbors=4):
    """Create symmetric bonds between nearest neighbors with probability `prob`."""
    n = len(pos)
    bonds = [[] for _ in range(n)]
    if prob <= 0.0 or k_bond <= 0.0:
        return bonds

    # Pairwise distances (small systems; OK to use full matrix)
    dists = np.linalg.norm(pos[:, None, :] - pos[None, :, :], axis=2)
    for i in range(n):
        nn_indices = np.argsort(dists[i])
        added = 0
        for j in nn_indices:
            if i == j or added >= max_neighbors:
                continue
            if np.random.rand() <= prob:
                bonds[i].append((int(j), float(rest_len), float(k_bond)))
                bonds[j].append((int(i), float(rest_len), float(k_bond)))
            added += 1
    return bonds


def init_plot(scatter, scatter_sizes, pos, bonds, line_segments):
    has_bonds = any(bonds)
    scatter.set_offsets(pos[:, :2])
    scatter.set_sizes(scatter_sizes)
    if has_bonds:
        segs = []
        for i, blist in enumerate(bonds):
            for j, _, _ in blist:
                if j > i:
                    segs.append([pos[i, :2], pos[j, :2]])
        line_segments.set_segments(segs)
    else:
        line_segments.set_segments([])
    return (scatter, line_segments) if has_bonds else (scatter,)


def update_plot(scatter, scatter_sizes, bonds, line_segments, p):
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
    return (scatter, line_segments) if has_bonds else (scatter,)


def step_iteration(sim, args, cheby_state, host_pos, iter_count, coll_overlaps, update_plot_fn):
    diag, _ = sim.jacobi_iteration(
        dt=args.dt,
        k_coll=args.k_coll,
        omega=args.omega,
        pd_scale=args.pd_scale,
        momentum_beta=args.momentum_beta,
        cheby_state=cheby_state,
        pred_buffer=sim.cl_iter_pos_A,  # use A as predictor for simple relaxation
        collect_diag=True,
        out_pos_host=host_pos,
    )
    iter_count[0] += 1
    coll_violation = max(-diag["coll_min"], 0.0)
    bond_violation = max(abs(diag["bond_min"]), abs(diag["bond_max"]))
    coll_overlaps.append((iter_count[0], coll_violation, bond_violation))
    return update_plot_fn(host_pos)

def run_momentum_sweep():
    mvals = [float(x) for x in args.momentum_list.split(',') if x.strip() != ""]
    results = []
    for mbeta in mvals:
        prepare_sim(sim, pos0, vel0, radius, mass, bonds, neighbor_margin=1.0)
        coll_series, bond_series = sim.run_jacobi_iterations(
            iterations=args.iters,
            dt=args.dt,
            k_coll=args.k_coll,
            omega=args.omega,
            pd_scale=args.pd_scale,
            momentum_beta=mbeta,
            cheby_enable=bool(args.cheby_enable),
            cheby_rho=args.cheby_rho,
            cheby_delay=args.cheby_delay,
            collect_diagnostics=True,
        )
        results.append((f"beta={mbeta}", coll_series, bond_series))
    plot_violation_series(results, args.iters)
    plt.show()
    sys.exit(0)

def run_animation():
    pos = pos0.copy()
    vel = vel0.copy()
    prepare_sim(sim, pos, vel, radius, mass, bonds, neighbor_margin=1.0)

    fig, ax = plt.subplots(figsize=(8, 8))
    margin = 1.0
    ax.set_xlim(np.min(pos[:, 0]) - margin, np.max(pos[:, 0]) + margin)
    ax.set_ylim(np.min(pos[:, 1]) - margin, np.max(pos[:, 1]) + margin)
    ax.set_aspect('equal')
    
    fig.canvas.draw()
    scatter_sizes = _scatter_sizes_from_radius(ax, fig, radius)
    scatter = ax.scatter(pos[:, 0], pos[:, 1], s=scatter_sizes, alpha=0.6)

    line_segments = LineCollection([], colors='gray', linewidths=0.5)
    ax.add_collection(line_segments)
    has_bonds = any(bonds)

    iter_count = [0]
    coll_overlaps = []
    cheby_state = {"rho": args.cheby_rho, "delay": args.cheby_delay, "omega_k": 1.0, "iter": 0} if args.cheby_enable else None
    host_pos = np.empty((num_atoms, 4), dtype=np.float32)

    update_plot_fn = lambda p: update_plot(scatter, scatter_sizes, bonds, line_segments, p)
    step_fn = lambda frame: step_iteration(sim, args, cheby_state, host_pos, iter_count, coll_overlaps, update_plot_fn)
    init_fn = lambda: init_plot(scatter, scatter_sizes, pos, bonds, line_segments)

    ani = FuncAnimation(fig, step_fn, frames=1000, init_func=init_fn, interval=30, blit=False, repeat=False)
    plt.show()

    if coll_overlaps:
        plot_violation_history(coll_overlaps)
        plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Jacobi Solver Convergence Test")
    parser.add_argument("--nx",       type=int,    default=10, help="number of atoms in x")
    parser.add_argument("--ny",       type=int,    default=10, help="number of atoms in y")
    parser.add_argument("--spacing",  type=float, default=0.35, help="spacing between atoms (radius=0.2, so <0.4 is overlapping)")
    parser.add_argument("--jitter",   type=float, default=0.02, help="random jitter")
    parser.add_argument("--atom_rad", type=float, default=0.2, help="atom radius")
    parser.add_argument("--k_coll",   type=float, default=100.0, help="collision stiffness")
    parser.add_argument("--omega",    type=float, default=0.5, help="relaxation factor")
    parser.add_argument("--pd_scale", type=float, default=0.0, help="projective inertia scale (0=pure Jacobi)")
    parser.add_argument("--momentum_beta", type=float, default=0.7, help="heavy-ball momentum coefficient")
    parser.add_argument("--cheby_enable", type=int, default=0, help="1=enable Chebyshev accel (two-step)")
    parser.add_argument("--cheby_rho", type=float, default=0.99, help="spectral radius estimate for Chebyshev")
    parser.add_argument("--cheby_delay", type=int, default=5, help="delay iterations before Chebyshev kicks in")
    parser.add_argument("--dt",       type=float, default=0.05, help="dummy dt for inertia term")
    parser.add_argument("--seed",     type=int,   default=587,    help="rng seed for reproducibility")
    parser.add_argument("--bond_prob",type=float, default=0.2,  help="probability to create a bond to a nearest neighbor")
    parser.add_argument("--bond_k",   type=float, default=200.0,help="bond stiffness")
    parser.add_argument("--bond_len", type=float, default=0.35, help="bond rest length (near spacing)")
    parser.add_argument("--iters",    type=int,   default=40,  help="solver iterations")
    parser.add_argument("--momentum_list", type=str, default="0.3,0.5,0.6,0.7,0.8,0.9", help="comma-separated momentum_beta values to sweep (disables animation)")
    #parser.add_argument("--momentum_list", type=str, default=None, help="comma-separated momentum_beta values to sweep (disables animation)")
    args = parser.parse_args()

    np.random.seed(args.seed)

    pos0, vel0, radius, mass = make_dense_block(args.nx, args.ny, args.spacing, args.jitter, args.atom_rad)
    num_atoms = len(pos0)

    sim = XPDB(num_atoms)
    sim.upload_data(pos0, vel0, radius, mass)
    bonds = make_random_bonds(pos0, prob=args.bond_prob, rest_len=args.bond_len, k_bond=args.bond_k)
    prepare_sim(sim, pos0, vel0, radius, mass, bonds, neighbor_margin=1.0)

    # Main entry: choose sweep vs animation (either/or)
    if args.momentum_list is not None:
        run_momentum_sweep()
    else:
        run_animation()
