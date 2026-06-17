#!/usr/bin/env python3
import os, sys, argparse
import numpy as np

_THIS_DIR = os.path.abspath(os.path.dirname(__file__))
_REPO_ROOT = os.path.abspath(os.path.join(_THIS_DIR, '..'))
for _p in (_THIS_DIR, _REPO_ROOT):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from pyBall.FFutils import MolViz2D, MolViz3D, attach_picker_2d, attach_picker_3d
from XPBD_2D import XPBD_2D
import XPBD_2D_utils as u2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--xyz', default=os.path.join(_ROOT, 'cpp/common_resources/xyz_mini/H2O.xyz'))
    ap.add_argument('--bAllNodes', action='store_true')
    ap.add_argument('--iters', type=int, default=200)
    ap.add_argument('--inner', type=int, default=10)
    ap.add_argument('--dt', type=float, default=0.1)
    ap.add_argument('--bmix', type=float, default=0.0)
    ap.add_argument('--view_scale', type=float, default=None)
    ap.add_argument('--view', default='2d', choices=['2d', '3d'])
    ap.add_argument('--pick_radius', type=float, default=0.5)
    ap.add_argument('--viz_every', type=int, default=1)
    ap.add_argument('--interval', type=float, default=0.001)
    ap.add_argument('--verbose', type=int, default=0)
    args = ap.parse_args()

    sim = XPBD_2D(num_atoms=1)  # will be resized by setup? no, XPBD_2D allocates by constructor; use natoms from mol
    # Build topology first to know natoms
    from pyBall.AtomicSystem import AtomicSystem
    mol = AtomicSystem(fname=args.xyz)
    natoms = len(mol.apos)
    sim = XPBD_2D(num_atoms=natoms)
    topo = u2.setup_from_mol(sim, mol, k_bond=200.0, perturbation=0.0, perturb_rot=0.0, bAllNodes=bool(args.bAllNodes), seed=0)

    nnode = int(topo['nnode'])
    neighs = topo['neighs']
    bks = topo['bks']
    port_local = topo['port_local']

    b3d = (str(args.view).lower() == '3d')
    if b3d:
        viz = MolViz3D(elems=topo.get('elems', None), show_labels=True)
        pick = attach_picker_3d(viz, pick_radius_px=int(max(2, float(args.pick_radius) * 100.0)), verbose=int(args.verbose))
    else:
        viz = MolViz2D(elems=topo.get('elems', None), show_labels=True, view_scale=args.view_scale, use_blit=True)
        pick = attach_picker_2d(viz, pick_radius=float(args.pick_radius), verbose=int(args.verbose))
    pick_prev = {'active': False, 'idx': None}
    PIN_MASS = 1e8
    NORMAL_MASS = 1.0

    def apply_pick(pos, rot):
        curr_active = pick.get('active') and pick.get('idx') is not None
        prev_active = pick_prev['active']
        if curr_active and not prev_active:
            ia = int(pick['idx'])
            M_prev, I_prev = sim.get_atom_mass(ia)
            sim.set_atom_mass(ia, (PIN_MASS, I_prev))
            if int(args.verbose) > 0: print(f"[pick] press idx={ia} mass M {M_prev}->{PIN_MASS}")
        if (not curr_active) and prev_active and (pick_prev['idx'] is not None):
            ia = int(pick_prev['idx'])
            M_prev, I_prev = sim.get_atom_mass(ia)
            sim.set_atom_mass(ia, (NORMAL_MASS, I_prev))
            if int(args.verbose) > 0: print(f"[pick] release idx={ia} mass M {M_prev}->{NORMAL_MASS}")
        pick_prev['active'] = bool(curr_active)
        pick_prev['idx'] = pick.get('idx')
        if curr_active:
            ia = int(pick['idx'])
            if b3d:
                sim.set_atom_pos(ia, np.asarray(pick['mouse3'][:2], dtype=np.float32))
            else:
                sim.set_atom_pos(ia, np.asarray(pick['mouse'], dtype=np.float32))
            sim.set_atom_vel(ia, [0.0, 0.0])
            sim.set_atom_omega(ia, 0.0)

    for it in range(int(args.iters)):
        sim.step_xpbd(nnode=nnode, dt=float(args.dt), iterations=int(args.inner), reset_lambda=True, bmix=float(args.bmix), reset_hb=True)
        if (it % int(args.viz_every)) == 0:
            pos, rot, vel, omega = sim.download_state()
            apply_pick(pos, rot)
            max_port, rms_port = u2.compute_port_error(pos, rot, neighs, bks, port_local, nnode)
            extra_pts = []
            extra_segs = []
            # Ports overlay (forcefield-specific mapping)
            for i in range(nnode):
                zi = rot[i]
                for k in range(4):
                    j = int(neighs[i, k])
                    if j < 0: continue
                    p = port_local[i, k]
                    pr = np.array([zi[0]*p[0] - zi[1]*p[1], zi[1]*p[0] + zi[0]*p[1]], dtype=np.float32)
                    tip = pos[i] + pr
                    extra_pts.append(tip)
                    extra_segs.append([pos[i], tip])
                    extra_segs.append([tip, pos[j]])
            extra_pts = np.array(extra_pts, dtype=np.float32) if len(extra_pts) else None
            info = f"it={it}  port_max={max_port:.3e} port_rms={rms_port:.3e}"
            if b3d:
                pos3 = np.zeros((pos.shape[0], 3), dtype=np.float32)
                pos3[:, :2] = pos
                extra_segments3 = []
                if extra_segs is not None:
                    for (a, b) in extra_segs:
                        aa = np.array([a[0], a[1], 0.0], dtype=np.float32)
                        bb = np.array([b[0], b[1], 0.0], dtype=np.float32)
                        extra_segments3.append(([aa[0], bb[0]], [aa[1], bb[1]], [aa[2], bb[2]]))
                viz.update(pos3, bonds=None, title='XPBD_2D ffdebug (3D view)', extra_segments=extra_segments3)
            else:
                viz.update(pos, bonds=None, title='XPBD_2D ffdebug', info=info, extra_points=extra_pts, extra_segments=extra_segs)
                viz.plt.pause(float(args.interval))

    viz.plt.show(block=True)


if __name__ == '__main__':
    main()
