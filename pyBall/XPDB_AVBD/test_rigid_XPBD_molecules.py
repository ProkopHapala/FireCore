import numpy as np
import argparse
import os
import sys

from numpy.random import default_rng

import matplotlib.pyplot as plt

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

_XPBD2D_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "XPBD_2D"))
if _XPBD2D_DIR not in sys.path:
    sys.path.insert(0, _XPBD2D_DIR)

from XPDB_new import XPDB_new, build_neighs_bk_from_bonds, make_bLs_bKs_from_neighs, make_bk_slots, linear_momentum, angular_momentum, run_RRsp3_PD, run_RRsp3_force
from XPTB_utils import load_xyz, masses_from_elems, perturb_state, write_xyz_with_ports, write_pdb_trajectory, plot_state_with_ports


def bonds_for_molecule(name, elems):
    # Debug test only: hardcoded bonds for the specific xyz files
    if name.lower().startswith('h2o'):
        # O(0)-H(1), O(0)-H(2)
        return [(0, 1), (0, 2)], 1  # nnode=1 (O is node, H are caps)
    if name.lower().startswith('ch2nh'):
        # C(0)=N(1), C(0)-H(2), C(0)-H(3), N(1)-H(4)
        return [(0, 1), (0, 2), (0, 3), (1, 4)], 2  # nnode=2 (C,N)
    raise ValueError(f"No bond template for '{name}'")


def reorder_nodes_first(elems, xyz, bonds, *, nnode, is_node=None):
    elems = list(elems)
    xyz = np.asarray(xyz, dtype=np.float32)
    n = int(len(elems))
    if xyz.shape != (n, 3):
        raise ValueError(f"reorder_nodes_first: xyz.shape={xyz.shape} expected ({n},3)")
    if is_node is None:
        is_node = np.zeros((n,), dtype=bool)
        is_node[:int(nnode)] = True
    is_node = np.asarray(is_node, dtype=bool)
    if is_node.shape != (n,):
        raise ValueError(f"reorder_nodes_first: is_node.shape={is_node.shape} expected ({n},)")
    if int(np.sum(is_node)) != int(nnode):
        raise ValueError(f"reorder_nodes_first: sum(is_node)={int(np.sum(is_node))} != nnode={int(nnode)}")

    order = np.concatenate([np.nonzero(is_node)[0], np.nonzero(~is_node)[0]]).astype(np.int32)
    perm_inv = np.empty((n,), dtype=np.int32)
    perm_inv[order] = np.arange(n, dtype=np.int32)

    elems2 = [elems[i] for i in order]
    xyz2 = xyz[order, :].copy()
    bonds2 = [(int(perm_inv[int(i)]), int(perm_inv[int(j)])) for (i, j) in bonds]
    return elems2, xyz2, bonds2


def load_molecule_any(xyz_path, *, bAllNodes=False):
    from pyBall.AtomicSystem import AtomicSystem
    mol = AtomicSystem(fname=xyz_path)
    if mol.bonds is None or len(mol.bonds) == 0:
        mol.findBonds()
    bonds = [(int(b[0]), int(b[1])) for b in mol.bonds]
    elems = list(mol.enames)
    xyz = np.asarray(mol.apos[:, :3], dtype=np.float32)
    if bool(bAllNodes):
        is_node = np.ones((len(elems),), dtype=bool)
    else:
        is_node = np.array([e != 'H' for e in elems], dtype=bool)
    nnode = int(np.sum(is_node))
    if nnode <= 0:
        raise ValueError(f"load_molecule_any: nnode={nnode} (bAllNodes={bAllNodes})")
    elems, xyz, bonds = reorder_nodes_first(elems, xyz, bonds, nnode=nnode, is_node=is_node)
    return elems, xyz, bonds, nnode


def run_case(xyz_path, *, method='force_explicit', dt=0.1, dt_force=0.01, iters=50, iters_force=2000, k_bond=200.0, k_rot=50.0, perturb_pos=0.1, perturb_rot=0.1, seed=0, damp_force=0.98,
             dump_xyz=None, dump_every=10, viz_force=False, viz_every=100, dump_pdb=None, pdb_every=10, plot_conv=True,
             pbd_relax=0.5, xpbd_reset_lambda=True, xpbd_variant='scalar', bAllNodes=False):
    name = os.path.basename(xyz_path)
    if name.lower() in ('h2o.xyz', 'ch2nh.xyz'):
        elems, xyz0, _q = load_xyz(xyz_path)
        bonds, nnode = bonds_for_molecule(name, elems)
    else:
        elems, xyz0, bonds, nnode = load_molecule_any(xyz_path, bAllNodes=bool(bAllNodes))
    m = masses_from_elems(elems)

    neighs, bks = build_neighs_bk_from_bonds(len(elems), bonds, max_deg=4)
    bLs, bKs = make_bLs_bKs_from_neighs(xyz0, neighs, k_bond=k_bond)

    quat0 = np.zeros((len(elems), 4), dtype=np.float32)
    quat0[:, 3] = 1.0

    vel0 = np.zeros((len(elems), 3), dtype=np.float32)
    omega0 = np.zeros((len(elems), 3), dtype=np.float32)

    rng = default_rng(seed)
    pos_init, quat_init = perturb_state(xyz0, quat0, perturb_pos, perturb_rot, rng)

    sim = XPDB_new(len(elems))
    sim.upload_rigid_state(pos_init, m, quat=quat_init, vel=vel0, omega=omega0)
    sim.upload_rigid_topology(neighs, bks, bLs, bKs)

    # atom port types: 0 none, 1 sp1, 2 sp2, 3 sp3
    atom_types = np.zeros((len(elems),), dtype=np.uint8)
    for i, e in enumerate(elems):
        if e in ('C', 'N'):
            atom_types[i] = 2  # sp2
        else:
            atom_types[i] = 0
    sim.upload_rigid_atom_types(atom_types)

    bkSlots = make_bk_slots(neighs, nnode=int(nnode))
    sim.upload_rigid_bk_slots(bkSlots)

    # Explicit-force local ports: node-only buffers (shape nnode)
    # For this simple test we use the initial geometry as body frame (quat=identity).
    port_local = np.zeros((int(nnode), 4, 4), dtype=np.float32)
    port_n = np.zeros((int(nnode),), dtype=np.uint8)
    for ia in range(int(nnode)):
        nn = 0
        for k in range(4):
            j = int(neighs[ia, k])
            if j < 0:
                continue
            v = xyz0[j] - xyz0[ia]
            port_local[ia, k, :3] = v
            nn += 1
        port_n[ia] = nn
    sim.upload_rigid_ports_local(port_local, port_n, nnode=int(nnode))
    sim.upload_rigid_node_stiffness_flat(bKs, nnode=int(nnode))

    Iiso = 0.4 * m * 1.0 * 1.0

    pos4, q4, v4, om4 = sim.download_rigid_state()
    p0 = pos4[:, :3].copy()

    # 1) projective (position-based) single iteration
    sim.rigid_projective_step(nnode=nnode, dt=dt, iterations=1)
    pos4, q4, v4, om4 = sim.download_rigid_state()
    p1 = pos4[:, :3].copy()
    v1 = (p1 - p0) / float(dt)

    P1 = linear_momentum(v1, m)
    L1 = angular_momentum(p1, v1, m, om4[:, :3], Iiso)

    print(f"\n=== {name} ===")
    print(f"nAtoms={len(elems)} nnode={nnode} dt={dt} iters={iters}")
    print(f"P(after 1 iter) = {P1}")
    print(f"L(after 1 iter) = {L1}")

    f_hist = []
    it_hist = []
    frames_pdb = []

    def emit_frame_3d(it, title=None):
        if not viz_force:
            return
        if (it % int(viz_every)) != 0 and it != int(iters) - 1:
            return
        pos4_it, q4_it, *_ = sim.download_rigid_state()
        pneigh_full = np.zeros((len(elems), 4, 4), dtype=np.float32)
        port_n_full = np.zeros((len(elems),), dtype=np.uint8)
        pneigh_full[:int(nnode), :, :] = port_local
        port_n_full[:int(nnode)] = port_n
        plot_state_with_ports(elems, pos4_it[:, :3], pneigh_full, port_n_full, title=title or f"{name} it={it}")

    def quat_rotate(q, v):
        q = np.asarray(q, dtype=np.float32)
        v = np.asarray(v, dtype=np.float32)
        qv = q[..., :3]
        qw = q[..., 3:4]
        t = 2.0 * np.cross(qv, v)
        return v + qw * t + np.cross(qv, t)

    viz2d = None
    viz3d = None
    pick2d = None
    pick3d = None
    active_view = getattr(run_case, '_view', 'none')

    def hook_toggle_view(fig):
        nonlocal active_view
        def on_key(event):
            nonlocal active_view
            if event.key == 'v':
                active_view = '3d' if active_view == '2d' else '2d'
        fig.canvas.mpl_connect('key_press_event', on_key)

    pneigh_full0 = np.zeros((len(elems), 4, 4), dtype=np.float32)
    port_n_full = np.zeros((len(elems),), dtype=np.uint8)
    pneigh_full0[:int(nnode), :, :] = port_local
    port_n_full[:int(nnode)] = port_n
    pneigh2d_full0 = pneigh_full0[:, :, :2].copy()

    def ensure_viz2d():
        nonlocal viz2d, pick2d
        if viz2d is not None:
            return
        import XPBD_2D_utils as utils2d
        from XPBD_2D_utils import LiveViz2D, attach_picker_2d
        viz2d = LiveViz2D(show_labels=True)
        pick2d = attach_picker_2d(viz2d, sim, pick_radius=0.5, verbose=1)
        hook_toggle_view(viz2d.fig)

    def ensure_viz3d():
        nonlocal viz3d, pick3d
        if viz3d is not None:
            return
        from XPTB_utils import LivePortViz, attach_picker_3d
        viz3d = LivePortViz(elems)
        pick3d = attach_picker_3d(viz3d, sim, pick_radius_px=20, verbose=1)
        hook_toggle_view(viz3d.fig)

    if viz_force and active_view == '2d':
        ensure_viz2d()
    if viz_force and active_view == '3d':
        ensure_viz3d()

    if viz_force and active_view != 'none':
        def emit_frame(it, title=None):
            if (it % int(viz_every)) != 0 and it != int(iters) - 1:
                return
            pos4_it, q4_it, *_ = sim.download_rigid_state()

            pneigh_full = pneigh_full0.copy()
            if int(nnode) > 0:
                vloc = pneigh_full[:int(nnode), :, :3]
                qn = q4_it[:int(nnode), :]
                pneigh_full[:int(nnode), :, :3] = quat_rotate(qn[:, None, :], vloc)

            if active_view == '2d':
                ensure_viz2d()
                pos2 = pos4_it[:, :2]
                viz2d.update(pos2, neighs=neighs, nnode=int(nnode), title=title or f"{name} it={it}", info=f"it={it}", port_local=pneigh_full[:, :, :2], port_n=port_n_full, rot=None)
                viz2d.fig.canvas.draw(); viz2d.fig.canvas.flush_events(); plt.pause(0.001)
            elif active_view == '3d':
                ensure_viz3d()
                viz3d.update(pos4_it[:, :3], pneigh_full, port_n_full, title=title or f"{name} it={it}")
                viz3d.fig.canvas.draw(); viz3d.fig.canvas.flush_events(); plt.pause(0.001)
    else:
        def emit_frame(it, title=None):
            return

    PIN_MASS = 1e8
    pick_prev = {"active": False, "idx": None}

    def apply_pick():
        nonlocal active_view
        pick = pick2d if active_view == '2d' else pick3d if active_view == '3d' else None
        if pick is None:
            return
        curr_active = bool(pick.get('active', False)) and (pick.get('idx', None) is not None)
        prev_active = bool(pick_prev['active'])

        if curr_active and (not prev_active):
            ia = int(pick['idx'])
            M_prev, I_prev = sim.get_atom_mass(ia)
            sim.set_atom_mass(ia, (PIN_MASS, I_prev))
        if (not curr_active) and prev_active and (pick_prev['idx'] is not None):
            ia = int(pick_prev['idx'])
            M_prev, I_prev = sim.get_atom_mass(ia)
            sim.set_atom_mass(ia, (1.0, I_prev))
        pick_prev['active'] = curr_active
        pick_prev['idx'] = int(pick['idx']) if curr_active else None

        if curr_active:
            ia = int(pick['idx'])
            if active_view == '2d':
                mouse = np.asarray(pick.get('mouse', [0.0, 0.0]), dtype=np.float32)
                sim.set_atom_pos(ia, [float(mouse[0]), float(mouse[1]), 0.0])
                sim.set_atom_vel(ia, [0.0, 0.0, 0.0])
                sim.set_atom_omega(ia, [0.0, 0.0, 0.0])
            elif active_view == '3d':
                mouse3 = np.asarray(pick.get('mouse3', [0.0, 0.0, 0.0]), dtype=np.float32)
                sim.set_atom_pos(ia, [float(mouse3[0]), float(mouse3[1]), float(mouse3[2])])
                sim.set_atom_vel(ia, [0.0, 0.0, 0.0])
                sim.set_atom_omega(ia, [0.0, 0.0, 0.0])

    def hook_toggle_view(fig):
        nonlocal active_view
        def on_key(event):
            nonlocal active_view
            if event.key == 'v':
                active_view = '3d' if active_view == '2d' else '2d'
        fig.canvas.mpl_connect('key_press_event', on_key)

    # Reset to initial state
    sim.upload_rigid_state(pos_init, m, quat=quat_init, vel=vel0, omega=omega0)
    sim.upload_rigid_atom_types(atom_types)

    if method == 'projective':
        sim.upload_rigid_state(pos_init, m, quat=quat_init, vel=vel0, omega=omega0)
        for it in range(int(iters)):
            apply_pick()
            sim.rigid_projective_step(nnode=int(nnode), dt=dt, iterations=1)
            emit_frame(it, title=f"{name} projective it={it}")
        pos4, q4, v4, om4 = sim.download_rigid_state()
        pos_prev = pos4[:, :3].copy()
    elif method == 'pbd_ports':
        sim.upload_rigid_state(pos_init, m, quat=quat_init, vel=vel0, omega=omega0)
        for it in range(int(iters)):
            apply_pick()
            sim.rigid_ports_pbd_step(nnode=int(nnode), iterations=1, relaxation=float(pbd_relax))
            emit_frame(it, title=f"{name} pbd it={it}")
        pos4, q4, v4, om4 = sim.download_rigid_state()
        pos_prev = pos4[:, :3].copy()
    elif method.startswith('xpbd_ports'):
        sim.upload_rigid_state(pos_init, m, quat=quat_init, vel=vel0, omega=omega0)
        for it in range(int(iters)):
            apply_pick()
            sim.rigid_ports_xpbd_step(
                nnode=int(nnode), dt=float(dt), iterations=1,
                reset_lambda=bool(xpbd_reset_lambda) if it == 0 else False,
                variant=xpbd_variant
            )
            emit_frame(it, title=f"{name} {method} it={it}")
        pos4, q4, v4, om4 = sim.download_rigid_state()
        pos_prev = pos4[:, :3].copy()
    else:
        pos_prev, max_err = run_RRsp3_PD(
            sim,
            nnode=int(nnode), dt=dt, iters=iters,
            pos_init=pos_init, m=m, quat_init=quat_init, vel0=vel0, omega0=omega0,
            neighs=neighs, bks=bks, bLs=bLs, bKs=bKs,
            atom_types=atom_types,
            verbose=True
        )
        print(f"max bond |d-L0| after proj iters = {max_err:.6e}")

    if method == 'force_explicit':
        nnode_force = int(nnode)
        sim.upload_rigid_state(pos_init, m, quat=quat_init, vel=vel0, omega=omega0)
        sim.upload_rigid_atom_types(atom_types)
        sim.upload_rigid_bk_slots(bkSlots)
        sim.upload_rigid_ports_local(port_local, port_n, nnode=int(nnode))
        pos_prev = pos_init.copy()
        print_every = max(1, int(iters_force // 20))
        dump_xyz_i = None
        if dump_xyz is not None:
            base, ext = os.path.splitext(dump_xyz)
            if ext == '':
                ext = '.xyz'
            dump_xyz_i = f"{base}_{name}{ext}"
            open(dump_xyz_i, 'w').close()
        def _on_force(it, fnorm, f4):
            f_hist.append(fnorm)
            it_hist.append(it)
            print(f"force it={it:6d} |F|={fnorm:.6e} log10|F|={np.log10(fnorm+1e-30): .3f}")

        def _on_frame(it):
            if dump_xyz_i is not None and ((it % int(dump_every)) == 0 or (it == int(iters_force) - 1)):
                pos4_it, q4_it, v4_it, om4_it = sim.download_rigid_state()
                pneigh_node = sim.download_buffer(sim.cl_rpneigh, (nnode_force * 4, 4)).reshape(nnode_force, 4, 4)
                pneigh_full = np.zeros((len(elems), 4, 4), dtype=np.float32)
                port_n_full = np.zeros((len(elems),), dtype=np.uint8)
                pneigh_full[:nnode_force, :, :] = pneigh_node
                port_n_full[:nnode_force] = port_n
                write_xyz_with_ports(dump_xyz_i, elems, pos4_it[:, :3], pneigh_full, port_n_full)

            if dump_pdb is not None and ((it % int(pdb_every)) == 0 or (it == int(iters_force) - 1)):
                pos4_it, q4_it, v4_it, om4_it = sim.download_rigid_state()
                frames_pdb.append(pos4_it[:, :3].copy())

            if viz_force and ((it % int(viz_every)) == 0 or (it == int(iters_force) - 1)):
                pos4_it, q4_it, v4_it, om4_it = sim.download_rigid_state()
                pneigh_node = sim.download_buffer(sim.cl_rpneigh, (nnode_force * 4, 4)).reshape(nnode_force, 4, 4)
                f4 = sim.download_buffer(sim.cl_rforce, (sim.num_atoms, 4))
                pneigh_full = np.zeros((len(elems), 4, 4), dtype=np.float32)
                port_n_full = np.zeros((len(elems),), dtype=np.uint8)
                pneigh_full[:nnode_force, :, :] = pneigh_node
                port_n_full[:nnode_force] = port_n
                plot_state_with_ports(elems, pos4_it[:, :3], pneigh_full, port_n_full, force=f4, title=f"{name} it={it}")

        run_RRsp3_force(
            sim,
            nnode=nnode_force, dt_force=dt_force, iters_force=iters_force, damp_force=damp_force,
            pos_init=pos_init, m=m, quat_init=quat_init, vel0=vel0, omega0=omega0,
            atom_types=atom_types,
            bkSlots=bkSlots,
            port_local=port_local, port_n=port_n,
            on_frame=_on_frame,
            on_force=_on_force
        )

    pos4, q4, v4, om4 = sim.download_rigid_state()
    if not np.all(np.isfinite(pos4)):
        raise RuntimeError(f"{method}: non-finite pos detected")
    if not np.all(np.isfinite(q4)):
        raise RuntimeError(f"{method}: non-finite quat detected")
    p = pos4[:, :3]
    if method == 'force_explicit':
        v = v4[:, :3]
        if not np.all(np.isfinite(v4)):
            raise RuntimeError(f"{method}: non-finite vel detected")
        if not np.all(np.isfinite(om4)):
            raise RuntimeError(f"{method}: non-finite omega detected")
        P = linear_momentum(v, m)
        L = angular_momentum(p, v, m, om4[:, :3], Iiso)
        print(f"force final |P|={np.linalg.norm(P):.6e} |L|={np.linalg.norm(L):.6e}")
    else:
        print(f"{method}: vel/omega not updated (PBD/XPBD step), skipping momentum print")
    pos_prev = p.copy()
    max_err = 0.0
    for i in range(len(elems)):
        for k in range(4):
            j = int(neighs[i, k])
            if j < 0:
                continue
            d = float(np.linalg.norm(pos_prev[j] - pos_prev[i]))
            err = abs(d - float(bLs[i, k]))
            if err > max_err:
                max_err = err
    print(f"max bond |d-L0| after {method} = {max_err:.6e}")

    if (method == 'force_explicit') and (len(f_hist) > 1) and plot_conv:
        plt.figure(figsize=(6, 4))
        plt.plot(it_hist, np.log10(np.array(f_hist) + 1e-30))
        plt.xlabel('iter')
        plt.ylabel('log10 |F|')
        plt.title(f'force convergence: {name}')
        plt.tight_layout()
        if not getattr(run_case, '_noshow', False):
            plt.show()

    if dump_pdb is not None:
        mol = os.path.splitext(os.path.basename(xyz_path))[0]
        if '{mol}' in dump_pdb:
            pdb_path = dump_pdb.format(mol=mol)
        else:
            base, ext = os.path.splitext(dump_pdb)
            if ext == '':
                ext = '.pdb'
            pdb_path = f"{base}_{mol}{ext}"
        if len(frames_pdb) == 0:
            raise RuntimeError(f"dump_pdb requested but no frames were collected (pdb_every={pdb_every})")
        write_pdb_trajectory(pdb_path, frames_pdb, elems, bonds)


'''
python3 pyBall/XPDB_AVBD/test_rigid_XPBD_molecules.py  --dt_force 0.001 --iters_force 1000 --perturb_pos 0.1 --perturb_rot 0.1 --damp_force 0.98 --viz_force --viz_every 50

python3 pyBall/XPDB_AVBD/test_rigid_XPBD_molecules.py --dt_force 0.001 --iters_force 2000 --perturb_pos 0.1 --perturb_rot 0.1 --damp_force 0.98 --dump_xyz traj.xyz --dump_every 10 --noshow


python test_rigid_XPBD_molecules.py --molecule ch2nh --dt_force 0.001 --iters_force 1000 --perturb_pos 0.1 --perturb_rot 0.1 --damp_force 0.98 --viz_force --viz_every 50 --dump_pdb traj.pdb


python test_rigid_XPBD_molecules.py --molecule ch2nh --dt_force 0.001 --iters_force 1000 --perturb_pos 0.1 --perturb_rot 0.1 --damp_force 0.98 --dump_pdb traj.pdb --pdb_every 10 --noshow

python test_rigid_XPBD_molecules.py --molecule ch2nh --dt_force 0.001 --iters_force 1000 --perturb_pos 0.1 --perturb_rot 0.1 --damp_force 0.98 --viz_force --viz_every 50 --dump_pdb traj_{mol}.pdb --pdb_every 10

python test_rigid_XPBD_molecules.py --molecule ch2nh --dt_force 0.001 --iters_force 1000 --perturb_pos 0.1 --perturb_rot 0.1 --damp_force 0.98 --viz_force --viz_every 10 --dump_pdb traj.pdb

'''

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--molecule', type=str, default='all')
    ap.add_argument('--xyz', type=str, default=None)
    ap.add_argument('--preset', type=str, default=None, choices=['h2o', 'ch2nh', 'hcooh', 'pyrrole', 'guanine', 'pentacene'])
    ap.add_argument('--bAllNodes', action='store_true')
    ap.add_argument('--method', type=str, default='force_explicit', choices=['force_explicit', 'projective', 'pbd_ports', 'xpbd_ports_scalar', 'xpbd_ports_vector'])
    ap.add_argument('--view', type=str, default='none', choices=['none', '2d', '3d'])
    ap.add_argument('--dt', type=float, default=0.1)
    ap.add_argument('--iters', type=int, default=50)
    ap.add_argument('--dt_force', type=float, default=0.01)
    ap.add_argument('--iters_force', type=int, default=2000)
    ap.add_argument('--k_bond', type=float, default=200.0)
    ap.add_argument('--k_rot', type=float, default=50.0)
    ap.add_argument('--perturb_pos', type=float, default=0.1)
    ap.add_argument('--perturb_rot', type=float, default=0.1)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--damp_force', type=float, default=0.98)
    ap.add_argument('--dump_xyz', type=str, default=None)
    ap.add_argument('--dump_every', type=int, default=10)
    ap.add_argument('--dump_pdb', type=str, default=None)
    ap.add_argument('--pdb_every', type=int, default=10)
    ap.add_argument('--viz_force', action='store_true')
    ap.add_argument('--viz_every', type=int, default=100)
    ap.add_argument('--plot_conv', action='store_true')
    ap.add_argument('--noshow', action='store_true')
    ap.add_argument('--pbd_relax', type=float, default=0.5)
    ap.add_argument('--xpbd_reset_lambda', type=int, default=1)
    args = ap.parse_args()

    base = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'cpp', 'common_resources', 'xyz'))
    xyz_h2o = os.path.join(base, 'H2O.xyz')
    xyz_ch2nh = os.path.join(base, 'CH2NH.xyz')
    presets = {
        'h2o': xyz_h2o,
        'ch2nh': xyz_ch2nh,
        'hcooh': os.path.join(base, 'HCOOH.xyz'),
        'pyrrole': os.path.join(base, 'pyrrol.xyz'),
        'guanine': os.path.join(base, 'guanine.xyz'),
        'pentacene': os.path.join(base, 'pentacene.xyz'),
    }
    run_case._noshow = bool(args.noshow)
    run_case._view = str(args.view)

    def resolve_xyz(args):
        if args.preset is not None:
            key = str(args.preset).lower()
            if key not in presets:
                raise ValueError(f"unknown preset '{args.preset}'")
            return presets[key]
        if args.xyz is not None:
            p = os.path.abspath(str(args.xyz))
            if not os.path.exists(p):
                raise FileNotFoundError(f"--xyz not found: {p}")
            return p
        if args.molecule.lower() in presets:
            return presets[str(args.molecule).lower()]
        if args.molecule.lower() in ('all',):
            return None
        p = os.path.abspath(str(args.molecule))
        if not os.path.exists(p):
            raise FileNotFoundError(f"--molecule not found: {p}")
        return p

    xyz_one = resolve_xyz(args)

    if (xyz_one is None) and (args.molecule in ('all', 'h2o')):
        run_case(xyz_h2o, dt=args.dt, iters=args.iters, dt_force=args.dt_force, iters_force=args.iters_force, k_bond=args.k_bond, k_rot=args.k_rot,
                 perturb_pos=args.perturb_pos, perturb_rot=args.perturb_rot, seed=args.seed + 0, damp_force=args.damp_force,
                 method=args.method if args.method != 'xpbd_ports_vector' else 'xpbd_ports',
                 xpbd_variant='vector' if args.method == 'xpbd_ports_vector' else 'scalar',
                 xpbd_reset_lambda=bool(args.xpbd_reset_lambda),
                 pbd_relax=args.pbd_relax,
                 bAllNodes=bool(args.bAllNodes), dump_xyz=args.dump_xyz, dump_every=args.dump_every, viz_force=args.viz_force, viz_every=args.viz_every,
                 dump_pdb=args.dump_pdb, pdb_every=args.pdb_every, plot_conv=bool(args.plot_conv))
    if (xyz_one is None) and (args.molecule in ('all', 'ch2nh')):
        run_case(xyz_ch2nh, dt=args.dt, iters=args.iters, dt_force=args.dt_force, iters_force=args.iters_force, k_bond=args.k_bond, k_rot=args.k_rot,
                 perturb_pos=args.perturb_pos, perturb_rot=args.perturb_rot, seed=args.seed + 1, damp_force=args.damp_force,
                 method=args.method if args.method != 'xpbd_ports_vector' else 'xpbd_ports',
                 xpbd_variant='vector' if args.method == 'xpbd_ports_vector' else 'scalar',
                 xpbd_reset_lambda=bool(args.xpbd_reset_lambda),
                 pbd_relax=args.pbd_relax,
                 bAllNodes=bool(args.bAllNodes), dump_xyz=args.dump_xyz, dump_every=args.dump_every, viz_force=args.viz_force, viz_every=args.viz_every,
                 dump_pdb=args.dump_pdb, pdb_every=args.pdb_every, plot_conv=bool(args.plot_conv))

    if xyz_one is not None:
        run_case(xyz_one, dt=args.dt, iters=args.iters, dt_force=args.dt_force, iters_force=args.iters_force, k_bond=args.k_bond, k_rot=args.k_rot,
                 perturb_pos=args.perturb_pos, perturb_rot=args.perturb_rot, seed=args.seed + 0, damp_force=args.damp_force,
                 method=args.method if args.method != 'xpbd_ports_vector' else 'xpbd_ports',
                 xpbd_variant='vector' if args.method == 'xpbd_ports_vector' else 'scalar',
                 xpbd_reset_lambda=bool(args.xpbd_reset_lambda),
                 pbd_relax=args.pbd_relax,
                 bAllNodes=bool(args.bAllNodes), dump_xyz=args.dump_xyz, dump_every=args.dump_every, viz_force=args.viz_force, viz_every=args.viz_every,
                 dump_pdb=args.dump_pdb, pdb_every=args.pdb_every, plot_conv=bool(args.plot_conv))
