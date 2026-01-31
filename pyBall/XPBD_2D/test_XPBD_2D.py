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


def run_simulation(sim, method, nnode, iters, *, dt=0.01, inner_iters=10, damp_vel=0.98, relax_dt=None, relax=0.7, bmix=0.0, topology=None, noshow=False, viz_every=1, interval=30, pick_radius=0.5, viz_substeps=False, view_scale=None, verbose=0):
    """Run simulation.

    Performance notes:
    - We only download state and redraw every `viz_every` iterations.
    - For interactive mode we use matplotlib blitting.
    - Real velocity (for MD) is computed as v = (x_new - x_old)/dt
    - Heavy-ball momentum (for relaxation) is: x_new = x' + dp*bmix where dp = x_new - x_prev
    """
    if topology is None:
        raise ValueError('run_simulation: topology is None')

    neighs = topology.get('neighs')
    bks = topology.get('bks')
    port_local = topology.get('port_local')
    port_n = topology.get('port_n')

    if noshow:
        if relax_dt is None:
            relax_dt = float(dt)
        for it in range(iters):
            if method == 'force':
                sim.step_explicit_force(nnode=nnode, dt=dt, damp=float(damp_vel), nsteps=1)
            elif method == 'xpbd_md':
                sim.step_xpbd(nnode=nnode, dt=dt, iterations=int(inner_iters), reset_lambda=(it == 0), bmix=float(bmix), reset_hb=(it == 0))
            elif method == 'xpbd_relax':
                # Use debug version for first few iterations to see lambda behavior
                if (int(verbose) > 0) and (it < 3):
                    sim.step_xpbd_debug(nnode=nnode, dt=float(relax_dt), iterations=int(inner_iters), reset_lambda=True, bmix=float(bmix), reset_hb=True, max_debug_steps=3)
                    # Capture detailed diagnostics after first inner iteration
                    diag = sim.download_constraint_diagnostics(nnode)
                    print(f"\n=== CONSTRAINT DIAGNOSTICS (iter {it}, after first inner) ===")
                    # Print all active constraints
                    for ni in range(nnode):
                        for pi in range(4):
                            C_val = diag['C'][ni, pi]
                            if abs(C_val) > 1e-12:  # Active constraint
                                lam = diag['lambda'][ni, pi]
                                K = diag['K'][ni, pi]
                                dtheta = diag['dtheta'][ni, pi]
                                dpos_i = diag['dpos_i'][ni, pi]
                                dpos_j = diag['dpos_j'][ni, pi]
                                r_world = diag['r_world'][ni, pi]
                                n = diag['n'][ni, pi]
                                print(f"  node[{ni}].port[{pi}]: C={C_val:12.4e} lam={lam:12.4e} K={K:8.2f} dtheta={dtheta:10.4e} r_world=({r_world[0]:10.4e},{r_world[1]:10.4e}) n=({n[0]:10.4e},{n[1]:10.4e}) dpos_i=({dpos_i[0]:10.4e},{dpos_i[1]:10.4e}) dpos_j=({dpos_j[0]:10.4e},{dpos_j[1]:10.4e})")
                    print("=== END DIAGNOSTICS ===\n")
                else:
                    sim.step_xpbd(nnode=nnode, dt=float(relax_dt), iterations=int(inner_iters), reset_lambda=True, bmix=float(bmix), reset_hb=True)
            elif method == 'pbd_md':
                sim.step_pbd_md(nnode=nnode, dt=float(dt), iterations=int(inner_iters), relax=float(relax), damp=float(damp_vel))
            elif method == 'pbd_relax':
                sim.step_pbd(nnode=nnode, iterations=int(inner_iters), relax=float(relax), bmix=float(bmix), reset_hb=(it == 0))
            else:
                raise ValueError(f"Unknown method {method}")

            if it % 10 == 0 or it == iters - 1:
                pos, rot, vel, omega = sim.download_state()
                npos_nan = int(np.isnan(pos).sum())
                nrot_nan = int(np.isnan(rot).sum())
                nvel_nan = int(np.isnan(vel).sum())
                nomega_nan = int(np.isnan(omega).sum())
                if (int(verbose) > 0) and ((npos_nan + nrot_nan + nvel_nan + nomega_nan) > 0):
                    print(f"[DEBUG] NaNs detected: pos={npos_nan} rot={nrot_nan} vel={nvel_nan} omega={nomega_nan}")
                P, L = compute_momentum_2d(pos, vel, omega)
                max_port, rms_port = compute_port_error(pos, rot, neighs, bks, port_local, nnode)
                # Download lambda for diagnostic
                if method == 'xpbd_relax' or method == 'xpbd_md':
                    lam = sim.download_lambda(nnode)
                    lam_min = np.min(lam)
                    lam_max = np.max(lam)
                    lam_mean = np.mean(lam)
                    print(f"iter {it:4d}: |P|={np.linalg.norm(P):.4e}, |L|={abs(L):.4e} port_max={max_port:.3e} port_rms={rms_port:.3e}  lambda[min={lam_min:.3e} max={lam_max:.3e} mean={lam_mean:.3e}]")
                else:
                    print(f"iter {it:4d}: |P|={np.linalg.norm(P):.4e}, |L|={abs(L):.4e} port_max={max_port:.3e} port_rms={rms_port:.3e}")
        return

    relax_dt_ref = {"dt": float(dt) if (relax_dt is None) else float(relax_dt)}

    viz = LiveViz2D(view_scale=view_scale)
    pick = attach_picker_2d(viz, sim, pick_radius=float(pick_radius), verbose=int(verbose))

    state = {"it": 0, "sub": 0, "reset_lambda": True}

    pick_state_prev = {"active": False, "idx": None}
    PIN_MASS = 1e8  # Very high mass to effectively pin atom
    NORMAL_MASS = 1.0  # Normal mass

    def apply_pick():
        # Detect pick start/end to set/restore mass
        curr_active = pick.get("active") and pick.get("idx") is not None
        prev_active = pick_state_prev["active"]
        
        # If pick just started, set high mass on the picked atom
        if curr_active and not prev_active:
            ia = int(pick["idx"])
            sim.set_atom_mass(ia, PIN_MASS)
            if int(verbose) > 0:
                print(f"[DEBUG] pick press idx={ia} mass set to {PIN_MASS}")
        
        # If pick just ended, restore normal mass on the previously picked atom
        if not curr_active and prev_active and pick_state_prev["idx"] is not None:
            ia = int(pick_state_prev["idx"])
            sim.set_atom_mass(ia, NORMAL_MASS)
            if int(verbose) > 0:
                print(f"[DEBUG] pick release idx={ia} mass restored to {NORMAL_MASS}")
        
        # Update previous state
        pick_state_prev["active"] = curr_active
        pick_state_prev["idx"] = pick.get("idx")
        
        # Apply position/velocity constraints if actively picking
        if curr_active:
            ia = int(pick["idx"])
            mouse_xy = np.asarray(pick.get("mouse", [0.0, 0.0]), dtype=np.float32)
            sim.set_atom_pos(ia, mouse_xy)
            sim.set_atom_vel(ia, [0.0, 0.0])
            sim.set_atom_omega(ia, 0.0)

    def inner_callback(itr):
        """Callback for inner iterations - visualize every viz_every Jacobi steps"""
        apply_pick()
        if itr % int(viz_every) == 0:
            pos, rot, vel, omega = sim.download_state()
            P, L = compute_momentum_2d(pos, vel, omega)
            max_port, rms_port = compute_port_error(pos, rot, neighs, bks, port_local, nnode)
            # Print positions for debugging divergence
            if int(verbose) > 0:
                print(f"\n=== Inner iteration {itr} ===")
                print(f"  port_max={max_port:.3e} port_rms={rms_port:.3e}")
                print(f"  Positions (first 5 nodes): {pos[:5]}")
            tag = f"{int(state['it'])}/{iters} inner={itr}"
            info = f'Iter: {tag}\n|P|: {np.linalg.norm(P):.4e}\n|L|: {abs(L):.4e}\nport max: {max_port:.3e}  rms: {rms_port:.3e}'
            viz.update(pos, neighs=neighs, nnode=nnode, title=f'XPBD_2D - {method}', info=info, port_local=port_local, port_n=port_n, rot=rot)
            viz.fig.canvas.draw()
            viz.fig.canvas.flush_events()
            # Use plt.pause for delay (interval is in ms, convert to seconds)
            # NOTE: even with interval==0, we need a tiny pause to process GUI events (mouse picking)
            if interval > 0:
                plt.pause(interval / 1000.0)
            else:
                plt.pause(0.001)

        return

    # For xpbd_relax, use callback to visualize inner iterations
    if method == 'xpbd_relax':
        for outer_it in range(int(iters)):
            apply_pick()
            sim.step_xpbd(nnode=nnode, dt=float(relax_dt_ref["dt"]), iterations=int(inner_iters), 
                         reset_lambda=True, bmix=float(bmix), reset_hb=True, callback=inner_callback)
            state["it"] = outer_it + 1
        plt.show(block=True)
        return

    # For pbd_md, use animation loop with picking
    if method == 'pbd_md':
        for outer_it in range(int(iters)):
            apply_pick()
            sim.step_pbd_md(nnode=nnode, dt=float(dt), iterations=int(inner_iters), relax=float(relax), damp=float(damp_vel))
            state["it"] = outer_it + 1
            # Download and visualize every viz_every steps
            if outer_it % int(viz_every) == 0:
                pos, rot, vel, omega = sim.download_state()
                P, L = compute_momentum_2d(pos, vel, omega)
                max_port, rms_port = compute_port_error(pos, rot, neighs, bks, port_local, nnode)
                tag = f"{int(state['it'])}/{iters}"
                info = f'Iter: {tag}\n|P|: {np.linalg.norm(P):.4e}\n|L|: {abs(L):.4e}\nport max: {max_port:.3e}  rms: {rms_port:.3e}\nvel_mean: {np.linalg.norm(vel):.3e}'
                viz.update(pos, neighs=neighs, nnode=nnode, title=f'XPBD_2D - {method}', info=info, port_local=port_local, port_n=port_n, rot=rot)
                viz.fig.canvas.draw()
                viz.fig.canvas.flush_events()
                if interval > 0:
                    plt.pause(interval / 1000.0)
        plt.show(block=True)
        return

    # For pbd_relax, use callback to visualize inner iterations
    if method == 'pbd_relax':
        for outer_it in range(int(iters)):
            apply_pick()
            sim.step_pbd(nnode=nnode, iterations=int(inner_iters), 
                        relax=float(relax), bmix=float(bmix), reset_hb=(outer_it == 0), callback=inner_callback)
            state["it"] = outer_it + 1
        plt.show(block=True)
        return

    def step_once():
        it = int(state["it"])
        apply_pick()
        if method == 'force':
            sim.step_explicit_force(nnode=nnode, dt=dt, damp=float(damp_vel), nsteps=1)
            state["it"] = it + 1
        elif method == 'xpbd_md':
            sim.step_xpbd(nnode=nnode, dt=dt, iterations=int(inner_iters), reset_lambda=bool(state["reset_lambda"]), bmix=float(bmix), reset_hb=bool(state["reset_lambda"]))
            state["reset_lambda"] = False
            state["it"] = it + 1
        elif method == 'xpbd_relax':
            if viz_substeps:
                sim.step_xpbd(nnode=nnode, dt=float(relax_dt_ref["dt"]), iterations=1, reset_lambda=True, bmix=float(bmix), reset_hb=True)
                state["sub"] += 1
                if state["sub"] >= int(inner_iters):
                    state["sub"] = 0
                    state["it"] = it + 1
            else:
                sim.step_xpbd(nnode=nnode, dt=float(relax_dt_ref["dt"]), iterations=int(inner_iters), reset_lambda=True, bmix=float(bmix), reset_hb=True)
                state["it"] = it + 1
        elif method == 'pbd_relax':
            sim.step_pbd(nnode=nnode, iterations=int(inner_iters), relax=float(relax), bmix=float(bmix), reset_hb=bool(state["reset_lambda"]))
            state["reset_lambda"] = False
            state["it"] = it + 1
        elif method == 'pbd_md':
            sim.step_pbd_md(nnode=nnode, dt=float(dt), iterations=int(inner_iters), relax=float(relax), damp=float(damp_vel))
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

'''
python test_XPBD_2D.py --method xpbd_relax --molecule ../../cpp/common_resources/xyz/pentacene.xyz --interval 200 --viz_scale 40 --perturb 1.0 --iters 10

python test_XPBD_2D.py --method xpbd_relax --molecule ../../cpp/common_resources/xyz/pentacene.xyz --interval 0.01 --viz_scale 100 --perturb 0.1 perturb_rot 0.0 --iters 1 --inner_iters 10

python test_XPBD_2D.py --method xpbd_relax --molecule ../../cpp/common_resources/xyz/pyrrol.xyz --interval 200 --viz_scale 10 --perturb 0.1 --iters 1 --inner_iters 10

python test_XPBD_2D.py --method xpbd_relax --molecule ../../cpp/common_resources/xyz/pyrrol.xyz --interval 200 --viz_scale 10 --perturb 0.1 --iters 1 --inner_iters 100 --relax_dt 0.01 


python test_XPBD_2D.py --method pbd_relax --molecule ../../cpp/common_resources/xyz/pyrrol.xyz --interval 0 --viz_scale 10 --perturb 0.1 --iters 1 --inner_iters 5000 --relax_dt 0.1 | tee OUT-pbd

python test_XPBD_2D.py --method pbd_relax --molecule ../../cpp/common_resources/xyz/pyrrol.xyz --interval 0 --viz_scale 10 --perturb 0.1 --iters 10000 --relax_dt 0.1 --method pbd_md

python test_XPBD_2D.py --method pbd_md    --molecule ../../cpp/common_resources/xyz/pyrrol.xyz --interval 0 --viz_scale 10 --perturb 0.1 --iters 10000 --relax_dt 0.1

'''

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Test XPBD_2D simulator')
    #parser.add_argument('--method',     default='force', choices=['force', 'xpbd_md', 'xpbd_relax'])
    #parser.add_argument('--method',     default='xpbd_relax', choices=['force', 'xpbd_md', 'xpbd_relax', 'pbd_relax'])
    parser.add_argument('--method',     default='pbd_relax', choices=['force', 'xpbd_md', 'xpbd_relax', 'pbd_relax', 'pbd_md'])

    parser.add_argument("--molecule",   type=str, default="../../cpp/common_resources/xyz/pentacene.xyz", help="Path to molecule file (.xyz, .mol, .mol2)")
    parser.add_argument('--iters',      type=int,   default=20000)
    parser.add_argument('--dt',         type=float, default=0.01)
    parser.add_argument('--relax_dt',   type=float, default=0.2)
    parser.add_argument('--inner_iters',type=int,   default=10)
    parser.add_argument('--relax',      type=float, default=0.7)
    parser.add_argument('--bmix',       type=float, default=0.2)
    parser.add_argument('--perturb',    type=float, default=0.1)
    parser.add_argument('--perturb_rot',type=float, default=0.1)
    parser.add_argument('--seed',       type=int,   default=0)
    parser.add_argument('--bAllNodes',  action='store_true')
    parser.add_argument('--noshow',     action='store_true')
    parser.add_argument('--viz_every',  type=int,   default=1)
    parser.add_argument('--interval',   type=int,   default=30, help='Animation frame delay in milliseconds (default: 30)')
    parser.add_argument('--viz_scale',  type=float, default=None, help='Fixed viewport scale (e.g., 20 for 20x20 view). If not set, auto-scales to molecule size.')
    parser.add_argument('--pick_radius',type=float, default=0.5)
    parser.add_argument('--viz_substeps',type=int,  default=0)
    parser.add_argument('--damp_vel',   type=float, default=0.9)
    parser.add_argument('--verbose',    type=int,   default=0)
    
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
    
    # Print initialization arrays for debugging
    print("\n=== INITIALIZATION ARRAYS ===")
    print(f"nnode={nnode}, n_atoms={n_atoms}")
    print(f"\nbkSlots ({len(topology['bkSlots'])} atoms):")
    for i, bk in enumerate(topology['bkSlots']): print(f"  atom[{i:2d}]: bkSlots=({bk[0]:3d}, {bk[1]:3d}, {bk[2]:3d}, {bk[3]:3d})")
    print(f"\nneighs ({len(topology['neighs'])} nodes):")
    for i, ng in enumerate(topology['neighs']):  print(f"  node[{i:2d}]: neighs=({ng[0]:3d}, {ng[1]:3d}, {ng[2]:3d}, {ng[3]:3d})")
    print(f"\nport_local (first 5 nodes, all ports):")
    for i in range(min(5, nnode)):
        for k in range(4):
            pl = topology['port_local'][i, k]
            if np.any(pl != 0):
                print(f"  node[{i}].port[{k}]: ({pl[0]:10.6f}, {pl[1]:10.6f})")
    print(f"\nstiffness ({len(topology['stiffness'])} nodes):")
    for i, st in enumerate(topology['stiffness']): print(f"  node[{i:2d}]: K=({st[0]:8.2f}, {st[1]:8.2f}, {st[2]:8.2f}, {st[3]:8.2f})")
    print("=== END INITIALIZATION ===\n")
    
    print(f"Running {args.iters} iterations with {args.method} solver... n_atoms={n_atoms} nnode={nnode}")
    run_simulation(sim, args.method, nnode, args.iters,
                   dt=args.dt, inner_iters=args.inner_iters, damp_vel=args.damp_vel,
                   relax_dt=args.relax_dt, relax=args.relax, bmix=args.bmix,
                   topology=topology, noshow=bool(args.noshow), viz_every=args.viz_every, interval=args.interval,
                   pick_radius=args.pick_radius, viz_substeps=bool(args.viz_substeps), view_scale=args.viz_scale, verbose=args.verbose)
    
    pos, rot, vel, omega = sim.download_state()
    P, L = compute_momentum_2d(pos, vel, omega)
    print(f"\nFinal: |P|={np.linalg.norm(P):.6e}, |L|={abs(L):.6e}")
