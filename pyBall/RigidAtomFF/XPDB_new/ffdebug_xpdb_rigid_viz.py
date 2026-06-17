#!/usr/bin/env python3
import os, sys, argparse
import numpy as np

_THIS_DIR = os.path.abspath(os.path.dirname(__file__))
_SHARED_DIR = os.path.abspath(os.path.join(_THIS_DIR, '..', 'shared'))
_REPO_ROOT = os.path.abspath(os.path.join(_THIS_DIR, os.pardir))
for _p in (_THIS_DIR, _SHARED_DIR, _REPO_ROOT):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from pyBall.FFutils import MolViz2D, MolViz3D, attach_picker_2d, attach_picker_3d
from XPDB_new import XPDB_new, build_neighs_bk_from_bonds, make_bLs_bKs_from_neighs, make_bk_slots
import XPTB_utils as xu


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--xyz', default=os.path.join(_ROOT, 'cpp/common_resources/xyz_mini/H2O.xyz'))
    ap.add_argument('--iters', type=int, default=500)
    ap.add_argument('--dt', type=float, default=0.01)
    ap.add_argument('--damp', type=float, default=0.98)
    ap.add_argument('--view', default='3d', choices=['2d', '3d'])
    ap.add_argument('--pick_radius_px', type=int, default=20)
    ap.add_argument('--pick_radius', type=float, default=0.5)
    ap.add_argument('--verbose', type=int, default=0)
    args = ap.parse_args()

    elems, xyz, q = xu.load_xyz(args.xyz)
    # bonds+ nnode heuristics: non-H nodes
    bonds = []
    # Use AtomicSystem for bonds if possible
    try:
        from pyBall.AtomicSystem import AtomicSystem
        mol = AtomicSystem(fname=args.xyz)
        if mol.bonds is None or len(mol.bonds) == 0:
            mol.findBonds()
        bonds = [(int(b[0]), int(b[1])) for b in mol.bonds]
    except Exception:
        # fallback: empty
        bonds = []
    if len(bonds) == 0:
        raise RuntimeError('ffdebug_xpdb_rigid_viz: no bonds found; provide xyz with bonds or use AtomicSystem bond inference')

    is_node = np.array([e != 'H' for e in elems], dtype=bool)
    nnode = int(np.sum(is_node))
    if nnode <= 0: raise RuntimeError('nnode<=0')

    # reorder nodes first
    order = np.concatenate([np.nonzero(is_node)[0], np.nonzero(~is_node)[0]]).astype(np.int32)
    perm_inv = np.empty((len(elems),), dtype=np.int32); perm_inv[order] = np.arange(len(elems), dtype=np.int32)
    elems2 = [elems[i] for i in order]
    xyz2 = xyz[order, :].copy()
    bonds2 = [(int(perm_inv[i]), int(perm_inv[j])) for (i, j) in bonds]

    neighs, bks = build_neighs_bk_from_bonds(len(elems2), bonds2, max_deg=4)
    bLs, bKs = make_bLs_bKs_from_neighs(xyz2, neighs, k_bond=200.0)
    bkSlots = make_bk_slots(neighs, nnode=nnode)

    m = xu.masses_from_elems(elems2)
    quat = np.zeros((len(elems2), 4), dtype=np.float32); quat[:, 3] = 1.0
    vel = np.zeros((len(elems2), 3), dtype=np.float32)
    omg = np.zeros((len(elems2), 3), dtype=np.float32)

    sim = XPDB_new(len(elems2))
    sim.upload_rigid_state(xyz2, m, quat=quat, vel=vel, omega=omg)
    sim.upload_rigid_topology(neighs, bks, bLs, bKs)
    sim.upload_rigid_bk_slots(bkSlots)

    # build simple ports in local frame from initial geometry (forcefield-specific)
    port_local = np.zeros((nnode, 4, 4), dtype=np.float32)
    port_n = np.zeros((nnode,), dtype=np.uint8)
    for ia in range(nnode):
        nn = 0
        for k in range(4):
            j = int(neighs[ia, k])
            if j < 0: continue
            v = xyz2[j] - xyz2[ia]
            port_local[ia, k, :3] = v
            nn += 1
        port_n[ia] = nn
    sim.upload_rigid_ports_local(port_local, port_n, nnode=nnode)
    sim.upload_rigid_node_stiffness_flat(bKs, nnode=nnode)

    b3d = (str(args.view).lower() == '3d')
    if b3d:
        viz = MolViz3D(elems=elems2, show_labels=True)
        pick = attach_picker_3d(viz, pick_radius_px=int(args.pick_radius_px), verbose=int(args.verbose))
    else:
        viz = MolViz2D(elems=elems2, show_labels=True, use_blit=True)
        pick = attach_picker_2d(viz, pick_radius=float(args.pick_radius), verbose=int(args.verbose))
        pick['_z_keep'] = None

    for it in range(int(args.iters)):
        if pick.get('active') and pick.get('idx') is not None:
            ia = int(pick['idx'])
            if b3d:
                sim.set_atom_pos(ia, pick['mouse3'])
            else:
                if pick.get('_z_keep', None) is None:
                    p4, _, _, _ = sim.download_rigid_state()
                    pick['_z_keep'] = float(p4[ia, 2])
                sim.set_atom_pos(ia, [float(pick['mouse'][0]), float(pick['mouse'][1]), float(pick['_z_keep'])])
            sim.set_atom_vel(ia, [0.0, 0.0, 0.0])
            sim.set_atom_omega(ia, [0.0, 0.0, 0.0])
        if (not pick.get('active', False)) and (not b3d):
            pick['_z_keep'] = None

        sim.rigid_force_explicit_step(nnode=nnode, dt=float(args.dt), nsteps=1, damp=float(args.damp))
        pos4, q4, v4, om4 = sim.download_rigid_state()
        pos = pos4[:, :3]
        f4 = sim.download_buffer(sim.cl_rforce, (sim.num_atoms, 4))
        if b3d:
            viz.update(pos, bonds=bonds2, force=f4[:, :3], title=f'XPDB rigid it={it}')
        else:
            viz.update(pos[:, :2], bonds=bonds2, title=f'XPDB rigid it={it}', info='(2D view: drag changes x,y only)')

    viz.plt.show(block=True)


if __name__ == '__main__':
    main()
