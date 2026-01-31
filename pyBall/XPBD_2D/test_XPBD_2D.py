#!/usr/bin/env python3
"""
test_XPBD_2D.py - Lightweight test script for 2D XPBD simulator

Usage:
    python test_XPBD_2D.py --method xpbd --n_atoms 10 --iters 200
    python test_XPBD_2D.py --xyz /path/to/H2O.xyz --iters 500
    python test_XPBD_2D.py --xyz /path/to/guanine.xyz --bAllNodes
"""

import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import itertools

# Add repo root for pyBall.AtomicSystem relative imports
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from XPBD_2D import XPBD_2D
from XPBD_2D_utils import (
    LiveViz2D, setup_from_mol,
    compute_momentum_2d, compute_port_error, attach_picker_2d
)

np.set_printoptions(linewidth=np.inf, threshold=np.inf)


def run_simulation(sim, method, nnode, iters, *, dt=0.01, inner_iters=10, damp_vel=0.98, topology=None, noshow=False, viz_every=1, interval=30, pick_radius=0.5, viz_substeps=False):
    """Run simulation.

    Performance notes:
    - We only download state and redraw every `viz_every` iterations.
    - For interactive mode we use matplotlib blitting.
    """
    if topology is None:
        raise ValueError('run_simulation: topology is None')

    neighs = topology.get('neighs')
    bks = topology.get('bks')
    port_local = topology.get('port_local')
    port_n = topology.get('port_n')

    if noshow:
        for it in range(iters):
            if method == 'force':
                sim.step_explicit_force(nnode=nnode, dt=dt, damp=float(damp_vel), nsteps=1)
            elif method == 'xpbd_md':
                sim.step_xpbd(nnode=nnode, dt=dt, iterations=int(inner_iters), reset_lambda=(it == 0))
            elif method == 'xpbd_relax':
                sim.step_xpbd(nnode=nnode, dt=dt, iterations=int(inner_iters), reset_lambda=(it == 0))
            else:
                raise ValueError(f"Unknown method {method}")

            if it % 20 == 0 or it == iters - 1:
                pos, rot, vel, omega = sim.download_state()
                P, L = compute_momentum_2d(pos, vel, omega)
                max_port, rms_port = compute_port_error(pos, rot, neighs, bks, port_local, nnode)
                print(f"iter {it:4d}: |P|={np.linalg.norm(P):.4e}, |L|={abs(L):.4e} port_max={max_port:.3e} port_rms={rms_port:.3e}")
        return

    viz = LiveViz2D()
    pick = attach_picker_2d(viz, sim, pick_radius=float(pick_radius))

    def apply_pick():
        if pick.get("active") and pick.get("idx") is not None:
            ia = int(pick["idx"])
            mouse_xy = np.asarray(pick.get("mouse", [0.0, 0.0]), dtype=np.float32)
            sim.set_atom_pos(ia, mouse_xy)
            sim.set_atom_vel(ia, [0.0, 0.0])
            sim.set_atom_omega(ia, 0.0)

    state = {"it": 0, "sub": 0, "reset_lambda": True}

    def step_once():
        it = int(state["it"])
        apply_pick()
        if method == 'force':
            sim.step_explicit_force(nnode=nnode, dt=dt, damp=float(damp_vel), nsteps=1)
            state["it"] = it + 1
        elif method == 'xpbd_md':
            sim.step_xpbd(nnode=nnode, dt=dt, iterations=int(inner_iters), reset_lambda=bool(state["reset_lambda"]))
            state["reset_lambda"] = False
            state["it"] = it + 1
        elif method == 'xpbd_relax':
            if viz_substeps:
                sim.step_xpbd(nnode=nnode, dt=dt, iterations=1, reset_lambda=bool(state["reset_lambda"]))
                state["reset_lambda"] = False
                state["sub"] += 1
                if state["sub"] >= int(inner_iters):
                    state["sub"] = 0
                    state["it"] = it + 1
            else:
                sim.step_xpbd(nnode=nnode, dt=dt, iterations=int(inner_iters), reset_lambda=bool(state["reset_lambda"]))
                state["reset_lambda"] = False
                state["it"] = it + 1
        else:
            raise ValueError(f"Unknown method {method}")

    def anim_init():
        pos, rot, vel, omega = sim.download_state()
        max_port, rms_port = compute_port_error(pos, rot, neighs, bks, port_local, nnode)
        info = f'Iter: 0/{iters}\n|P|: 0\n|L|: 0\nport max: {max_port:.3e}  rms: {rms_port:.3e}'
        return viz.update(pos, neighs=neighs, nnode=nnode, title=f'XPBD_2D - {method}', info=info, port_local=port_local, port_n=port_n, rot=rot)

    ani_ref = {"ani": None}

    def anim_update(frame):
        # advance simulation by viz_every steps
        for _ in range(int(viz_every)):
            if int(state["it"]) >= int(iters):
                break
            step_once()
            apply_pick()

        pos, rot, vel, omega = sim.download_state()
        P, L = compute_momentum_2d(pos, vel, omega)
        max_port, rms_port = compute_port_error(pos, rot, neighs, bks, port_local, nnode)
        tag = f"{int(state['it'])}/{iters}" if not viz_substeps else f"{int(state['it'])}/{iters} sub={int(state['sub'])}"
        info = f'Iter: {tag}\n|P|: {np.linalg.norm(P):.4e}\n|L|: {abs(L):.4e}\nport max: {max_port:.3e}  rms: {rms_port:.3e}'
        if int(state["it"]) % 20 == 0 or int(state["it"]) == int(iters) - 1:
            print(f"iter {int(state['it']):4d}: |P|={np.linalg.norm(P):.4e}, |L|={abs(L):.4e} port_max={max_port:.3e} port_rms={rms_port:.3e}")

        artists = viz.update(pos, neighs=neighs, nnode=nnode, title=f'XPBD_2D - {method}', info=info, port_local=port_local, port_n=port_n, rot=rot)

        if int(state["it"]) >= int(iters) and ani_ref["ani"] is not None:
            ani_ref["ani"].event_source.stop()
        return artists

    frames = itertools.count()  # run until we stop the event source
    ani = FuncAnimation(viz.fig, anim_update, init_func=anim_init, frames=frames, interval=int(interval), blit=False, repeat=False, cache_frame_data=False)
    ani_ref["ani"] = ani
    # keep a reference to avoid GC
    viz._ani = ani
    plt.show(block=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Test XPBD_2D simulator')
    #parser.add_argument('--method',     default='force', choices=['force', 'xpbd_md', 'xpbd_relax'])
    parser.add_argument('--method',     default='xpbd_relax', choices=['force', 'xpbd_md', 'xpbd_relax'])

    parser.add_argument("--molecule",   type=str, default="../../cpp/common_resources/xyz/pentacene.xyz", help="Path to molecule file (.xyz, .mol, .mol2)")
    parser.add_argument('--iters',      type=int,   default=20000)
    parser.add_argument('--dt',         type=float, default=0.01)
    parser.add_argument('--inner_iters',type=int,   default=10)
    parser.add_argument('--perturb',    type=float, default=0.1)
    parser.add_argument('--perturb_rot',type=float, default=0.1)
    parser.add_argument('--seed',       type=int,   default=0)
    parser.add_argument('--bAllNodes',  action='store_true')
    parser.add_argument('--noshow',     action='store_true')
    parser.add_argument('--viz_every',  type=int,   default=1)
    parser.add_argument('--interval',   type=int,   default=30)
    parser.add_argument('--pick_radius',type=float, default=0.5)
    parser.add_argument('--viz_substeps',type=int,  default=0)
    parser.add_argument('--damp_vel',   type=float, default=0.9)
    
    args = parser.parse_args()
    
    print(f"XPBD_2D Test: method={args.method}, molecule={args.molecule}")

    if args.noshow:
        plt.ioff()

    from pyBall.AtomicSystem import AtomicSystem
    mol = AtomicSystem(fname=args.molecule)
    n_atoms = len(mol.apos)
    sim = XPBD_2D(num_atoms=n_atoms)

    topology = setup_from_mol(sim, mol, k_bond=200.0,
                              perturbation=args.perturb, perturb_rot=args.perturb_rot,
                              bAllNodes=args.bAllNodes, seed=args.seed)
    nnode = topology['nnode']
    
    print(f"Running {args.iters} iterations with {args.method} solver... n_atoms={n_atoms} nnode={nnode}")
    run_simulation(sim, args.method, nnode, args.iters,
                   dt=args.dt, inner_iters=args.inner_iters, damp_vel=args.damp_vel,
                   topology=topology, noshow=bool(args.noshow), viz_every=args.viz_every, interval=args.interval,
                   pick_radius=args.pick_radius, viz_substeps=bool(args.viz_substeps))
    
    pos, rot, vel, omega = sim.download_state()
    P, L = compute_momentum_2d(pos, vel, omega)
    print(f"\nFinal: |P|={np.linalg.norm(P):.6e}, |L|={abs(L):.6e}")
