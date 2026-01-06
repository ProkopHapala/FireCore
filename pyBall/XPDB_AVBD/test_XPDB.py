# test_XPDB.py
import argparse
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.animation import FuncAnimation
from XPDB import XPDB

# PYOPENCL_COMPILER_OUTPUT=1
os.environ['PYOPENCL_COMPILER_OUTPUT'] = os.environ.get('PYOPENCL_COMPILER_OUTPUT', '0')

# Scene bounds / external forces (constants)
BOX_MIN = (-10.0, 0.0, -5.0)
BOX_MAX = (5.0, 5.0, 5.0)
GRAVITY = (0.0, -9.81, 0.0)


def make_scene(args, no_bonds):
    """Create initial positions, velocities, radii, masses, bonds based on args."""
    num_atoms = args.chains * args.atoms_per_chain
    rng = np.random.default_rng(args.seed)

    pos = np.zeros((num_atoms, 3), dtype=np.float32)
    vel = np.zeros((num_atoms, 3), dtype=np.float32)
    radius = np.full(num_atoms, args.atom_rad, dtype=np.float32)
    mass = np.full(num_atoms, args.atom_mass, dtype=np.float32)
    bonds = [[] for _ in range(num_atoms)]

    span_x = BOX_MAX[0] - BOX_MIN[0]
    span_y = BOX_MAX[1] - BOX_MIN[1]
    base_y = BOX_MIN[1] + 0.3 * span_y
    # place chains evenly within box width
    stride_x = span_x * 0.6 / max(1, args.chains - 1) if args.chains > 1 else 0.0
    start_x = BOX_MIN[0] + 0.2 * span_x

    for c in range(args.chains):
        offset_x = start_x + c * stride_x
        start_idx = c * args.atoms_per_chain

        for i in range(args.atoms_per_chain):
            idx = start_idx + i
            jitter = rng.normal(scale=args.spread, size=3).astype(np.float32)
            pos[idx] = [offset_x, base_y + i * args.bond_len * 0.9, 0.0] + jitter

            if i > 0 and not no_bonds:
                prev = idx - 1
                bonds[idx].append((prev, args.bond_len, args.stiffness))
                bonds[prev].append((idx, args.bond_len, args.stiffness))

        if not no_bonds:
            mass[start_idx + args.atoms_per_chain - 1] = 1e6

    return pos, vel, radius, mass, bonds

def _scatter_sizes_from_radius(ax, fig, radius):
    """Compute scatter areas (points^2) so marker radius matches data-space radius."""
    p0 = ax.transData.transform((0, 0))
    p1 = ax.transData.transform((1, 1))
    sx = abs(p1[0] - p0[0])
    sy = abs(p1[1] - p0[1])
    sp = min(sx, sy)  # pixels per data unit, conservative for both axes
    dpi = fig.dpi
    pr_points = radius * sp * 72.0 / dpi  # convert pixel radius to points
    return np.pi * (pr_points**2)  # scatter size expects area in points^2

def make_animation(sim, bonds, radius, args, no_bonds):
    fig, ax = plt.subplots()
    # Slightly larger viewport than box for visibility
    margin = 2.0
    ax.set_xlim(BOX_MIN[0] - margin, BOX_MAX[0] + margin)
    ax.set_ylim(BOX_MIN[1] - margin, BOX_MAX[1] + margin)
    scatter       = ax.scatter([], [], s=1)
    line_segments = LineCollection([], colors='gray', linewidths=0.5)
    ax.add_collection(line_segments)
    line_segments.set_visible(not no_bonds)
    # Draw box outline
    box_x = [BOX_MIN[0], BOX_MAX[0], BOX_MAX[0], BOX_MIN[0], BOX_MIN[0]]
    box_y = [BOX_MIN[1], BOX_MIN[1], BOX_MAX[1], BOX_MAX[1], BOX_MIN[1]]
    (box_line,) = ax.plot(box_x, box_y, 'k--', linewidth=0.8)
    # Draw parabolic bottom curve
    px = np.linspace(BOX_MIN[0]-margin, BOX_MAX[0]+margin, 200)
    py = (args.bottom_y if args.bottom_y is not None else BOX_MIN[1]) + args.bottom_a * (px - 0.5*(BOX_MIN[0]+BOX_MAX[0]))**2
    (bottom_line,) = ax.plot(px, py, 'r:', linewidth=1.0)

    fig.canvas.draw()  # ensure transforms are ready
    scatter_sizes = _scatter_sizes_from_radius(ax, fig, radius)

    def update_plot(p):
        scatter.set_offsets(p[:, :2])
        scatter.set_sizes(scatter_sizes)

        if no_bonds:
            line_segments.set_segments([])
        else:
            segs = []
            for i in range(NUM_ATOMS):
                for j, _, _ in bonds[i]:
                    if j > i:  # draw each bond once
                        segs.append([p[i, :2], p[j, :2]])
            line_segments.set_segments(segs)
        return scatter, line_segments, box_line, bottom_line,

    return fig, update_plot

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="XPDB chain drop demo")
    parser.add_argument("--chains",          type=int,   default=1,     help="number of chains")
    parser.add_argument("--atoms_per_chain", type=int,   default=64,     help="atoms per chain")
    parser.add_argument("--group_size",      type=int,   default=64,    help="OpenCL group size")
    parser.add_argument("--bond_len",        type=float, default=0.5,   help="bond rest length")
    parser.add_argument("--atom_rad",        type=float, default=0.2,   help="atom radius")
    parser.add_argument("--atom_mass",       type=float, default=1.0,   help="atom mass")
    parser.add_argument("--stiffness",       type=float, default=0.0,   help="bond stiffness")
    parser.add_argument("--dt",              type=float, default=0.05,  help="timestep")
    parser.add_argument("--bonds",           type=int,   default=0,     help="1=use bonds, 0=free particles")
    parser.add_argument("--spread",          type=float, default=0.2,   help="random spread for initial pos")
    parser.add_argument("--seed",            type=int,   default=0,     help="rng seed")
    parser.add_argument("--k_coll",          type=float, default=80.0,  help="collision stiffness")
    parser.add_argument("--box_k",           type=float, default=40.0,  help="box wall stiffness")
    parser.add_argument("--bottom_c",        type=float, default=-5.0,  help="parabolic bottom coefficient (0 disables)")
    parser.add_argument("--bottom_y",        type=float, default=None,  help="parabolic bottom reference height (defaults to box min y)")
    parser.add_argument("--bottom_a",        type=float, default=0.05,  help="parabolic bottom curvature coefficient")
    parser.add_argument("--substeps",        type=int,   default=2,     help="physics substeps per animation frame (predict+solve)")
    parser.add_argument("--solver_iters",    type=int,   default=5,     help="constraint solver iterations per substep")
    parser.add_argument("--verlet_every",    type=int,   default=1,     help="rebuild neighbors every N frames")
    parser.add_argument("--diag_every",      type=int,   default=1,     help="diagnostics print cadence (frames)")
    parser.add_argument("--clamp_z",         type=int,   default=1,     help="force z=0 each step for 2D debug")
    parser.add_argument("--pd_scale",        type=float, default=0.0,   help="projective inertia diagonal scale (mass/dt^2 term)")
    args = parser.parse_args()

    # --------------------------
    # SCENE SETUP
    # --------------------------
    GROUP_SIZE = args.group_size
    NUM_CHAINS = args.chains
    ATOMS_PER_CHAIN = args.atoms_per_chain
    NUM_ATOMS = NUM_CHAINS * ATOMS_PER_CHAIN

    # Physics Params
    BOND_LEN = args.bond_len
    ATOM_RAD = args.atom_rad
    ATOM_MASS = args.atom_mass
    STIFFNESS = args.stiffness
    DT = args.dt
    SUBSTEPS = args.substeps          # physics substeps per animation frame (predict + solver)
    SOLVER_ITERS = args.solver_iters  # Jacobi constraint iterations per substep
    VERLET_EVERY = args.verlet_every  # rebuild neighbor lists every N frames
    DIAG_EVERY = args.diag_every      # print diagnostics every N frames

    # Initialize scene data
    no_bonds = (args.bonds == 0)
    pos, vel, radius, mass, bonds = make_scene(args, no_bonds)

    # --------------------------
    # SETTINGS SUMMARY
    # --------------------------
    print(f"chains={NUM_CHAINS}, atoms_per_chain={ATOMS_PER_CHAIN}, bonds={'on' if not no_bonds else 'off'}, "
          f"dt={DT}, substeps={SUBSTEPS}, solver_iters={SOLVER_ITERS}, verlet_every={VERLET_EVERY}, diag_every={DIAG_EVERY}")
    print(f"box_min={BOX_MIN}, box_max={BOX_MAX}, bottom_y={args.bottom_y if args.bottom_y is not None else BOX_MIN[1]}, "
          f"bottom_c={args.bottom_c}, bottom_a={args.bottom_a}")
    print("initial positions (x,y,z,R):")
    pos_print = np.concatenate([pos, radius[:, None]], axis=1)
    np.set_printoptions(precision=3, suppress=True)
    print(pos_print)

    # --------------------------
    # INIT SIMULATION
    # --------------------------
    sim = XPDB(NUM_ATOMS, GROUP_SIZE)
    sim.upload_data(pos, vel, radius, mass)
    if not no_bonds:
        sim.upload_bonds(bonds)
        sim.update_verlet(margin_sq=(2.0*ATOM_RAD + 1.0)**2)
    else:
        sim.set_empty_bonds()

    # --------------------------
    # ANIMATION LOOP
    # --------------------------
    fig, update_plot = make_animation(sim, bonds, radius, args, no_bonds)
    step_id = {"value": 0}  # mutable holder

    def step_and_draw(frame):
        sid = step_id["value"]
        if sid % VERLET_EVERY == 0:
            sim.update_verlet(margin_sq=(2.0*ATOM_RAD + 0.5)**2)
        for _ in range(SUBSTEPS):
            sim.step(
                dt=DT,
                iterations=SOLVER_ITERS,
                k_coll=args.k_coll,
                gravity=GRAVITY,
                box_min=BOX_MIN,
                box_max=BOX_MAX,
                box_k=args.box_k,
                bottom_y=args.bottom_y if args.bottom_y is not None else BOX_MIN[1],
                bottom_c=args.bottom_c,
                bottom_x0=0.5*(BOX_MIN[0]+BOX_MAX[0]),
                bottom_a=args.bottom_a,
                clamp_z=args.clamp_z,
            )
        p = sim.get_positions()
        if sid % DIAG_EVERY == 0:
            diag = sim.diagnostics()
            if diag:
                print(f"[step {sid:6d}] bond: {diag['bond_min']:+8.4f} {diag['bond_max']:+8.4f} | coll: {diag['coll_min']:+8.4f} {diag['coll_max']:+8.4f}")
        step_id["value"] = sid + 1
        return update_plot(p)

    ani = FuncAnimation(fig, step_and_draw, frames=20000, interval=20, blit=True, repeat=False)
    plt.show()