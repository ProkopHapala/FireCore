#!/usr/bin/env python3
import os, sys, argparse
import numpy as np
import pyopencl as cl

_THIS_DIR = os.path.abspath(os.path.dirname(__file__))
_REPO_ROOT = os.path.abspath(os.path.join(_THIS_DIR, os.pardir))
for _p in (_THIS_DIR, _REPO_ROOT):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from pyBall.FFutils import MolViz2D, MolViz3D, attach_picker_2d, attach_picker_3d
from XPDB_new import XPDB_new
from test_TiledJacobi_molecules import load_molecule_topology_mmffl


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mol', default=os.path.join(_ROOT, 'cpp/common_resources/mol/formic_acid.mol2'))
    ap.add_argument('--iters', type=int, default=200)
    ap.add_argument('--inner', type=int, default=10)
    ap.add_argument('--dt', type=float, default=0.01)
    ap.add_argument('--k_coll', type=float, default=0.0)
    ap.add_argument('--Rmax', type=float, default=2.0)
    ap.add_argument('--view', default='3d', choices=['2d', '3d'])
    ap.add_argument('--pick_radius_px', type=int, default=20)
    ap.add_argument('--pick_radius', type=float, default=0.5)
    ap.add_argument('--verbose', type=int, default=0)
    args = ap.parse_args()

    topo, apos, bond_idx, bond_L, bond_K, bonds_adj, bond_type_mask = load_molecule_topology_mmffl(
        args.mol,
        n_max_bonded=16,
        default_K=200.0,
        add_angle=False,
        add_pi=False,
        add_epair=False,
        print_topology=False,
    )
    apos = np.asarray(apos, dtype=np.float32)
    bond_idx = np.asarray(bond_idx, dtype=np.int32)
    bond_L = np.asarray(bond_L, dtype=np.float32)
    bond_K = np.asarray(bond_K, dtype=np.float32)
    if apos.ndim != 2 or apos.shape[1] != 3:
        raise ValueError(f"apos shape {apos.shape} expected (N,3)")
    if bond_idx.ndim != 2 or bond_idx.shape[0] != apos.shape[0]:
        raise ValueError(f"bond_idx shape {bond_idx.shape} expected (N,M)")
    rad = np.ones((apos.shape[0],), dtype=np.float32) * 0.2
    mass = np.ones((apos.shape[0],), dtype=np.float32)

    n = apos.shape[0]
    sim = XPDB_new(n)
    sim.upload_data(apos, np.zeros_like(apos), rad, mass)
    sim.upload_bonds_fixed_from_arrays(bond_idx, bond_L, bond_K)

    elems = topo.get('type_names', None)
    bonds = []
    # derive bonds list from bond_idx (unique)
    for i in range(n):
        for k in range(bond_idx.shape[1]):
            j = int(bond_idx[i, k])
            if j >= 0 and j > i:
                bonds.append((i, j))

    b3d = (str(args.view).lower() == '3d')
    if b3d:
        viz = MolViz3D(elems=elems, show_labels=True)
        pick = attach_picker_3d(viz, pick_radius_px=int(args.pick_radius_px), verbose=int(args.verbose))
    else:
        viz = MolViz2D(elems=elems, show_labels=True, use_blit=True)
        pick = attach_picker_2d(viz, pick_radius=float(args.pick_radius), verbose=int(args.verbose))
        pick['_z_keep'] = None

    def _set_atom_pos_tiled(ia, xyz):
        ia = int(ia)
        xyz = np.asarray(xyz, dtype=np.float32).reshape(3,)
        p4 = np.zeros((1, 4), dtype=np.float32)
        p4[0, :3] = xyz
        # Keep solver stable: overwrite current + prev + vel + momentum for this atom
        cl.enqueue_copy(sim.queue, sim.cl_pos, p4, device_offset=ia * 16).wait()
        cl.enqueue_copy(sim.queue, sim.cl_prev_pos, p4, device_offset=ia * 16).wait()
        cl.enqueue_copy(sim.queue, sim.cl_vel, np.zeros((1, 4), dtype=np.float32), device_offset=ia * 16).wait()
        cl.enqueue_fill_buffer(sim.queue, sim.cl_momentum, np.float32(0.0), ia * 16, 16).wait()

    for it in range(int(args.iters)):
        if pick.get('active') and pick.get('idx') is not None:
            ia = int(pick['idx'])
            if b3d:
                _set_atom_pos_tiled(ia, np.asarray(pick['mouse3'], dtype=np.float32).reshape(3,))
            else:
                if pick.get('_z_keep', None) is None:
                    pos4 = sim.download_buffer(sim.cl_pos, (sim.num_atoms, 4))
                    pick['_z_keep'] = float(pos4[ia, 2])
                _set_atom_pos_tiled(ia, [float(pick['mouse'][0]), float(pick['mouse'][1]), float(pick['_z_keep'])])
        if (not pick.get('active', False)) and (not b3d):
            pick['_z_keep'] = None

        sim.build_local_topology(float(args.Rmax))
        sim.solve_cluster_jacobi_local(dt=float(args.dt), inner_iters=int(args.inner), k_coll=float(args.k_coll), omega=0.8)
        pos4 = sim.download_buffer(sim.cl_pos, (sim.num_atoms, 4))
        pos = pos4[:, :3]
        if b3d:
            viz.update(pos, bonds=bonds, title=f'XPDB tiled it={it}')
        else:
            viz.update(pos[:, :2], bonds=bonds, title=f'XPDB tiled it={it}', info='(2D view: drag changes x,y only)')

    viz.plt.show(block=True)


if __name__ == '__main__':
    main()
