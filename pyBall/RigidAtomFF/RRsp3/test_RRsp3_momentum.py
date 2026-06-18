import sys, os
import argparse
import numpy as np

_THIS_DIR = os.path.abspath(os.path.dirname(__file__))
_SHARED_DIR = os.path.abspath(os.path.join(_THIS_DIR, '..', 'shared'))
for _p in (_THIS_DIR, _SHARED_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from RRsp3 import RRsp3, build_neighs_bk_from_bonds, make_bk_slots_clustered, make_exclusions_1st_2nd
from XPTB_utils import pack_molecules_contiguous, make_h2o_geometry, masses_from_elems, perturb_state
from RRsp3_utils import make_ports_from_neighs, write_xyz_frame

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--steps', type=int, default=1)
    ap.add_argument('--traj', type=str, default='')
    ap.add_argument('--traj_every', type=int, default=1)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--pos_scale', type=float, default=0.05)
    ap.add_argument('--rot_scale', type=float, default=0.0)
    args = ap.parse_args()

    group_size = 64
    dt = 0.1

    pos_h2o, _ = make_h2o_geometry(add_angle=False)
    elems = ['O','H','H']
    bonds = [(0,1),(0,2)]
    nnode = 1

    # Put two molecules close enough for collisions but not extremely overlapping
    shift = np.array([2.1, 0.0, 0.0], dtype=np.float32)

    mols = [
        {'elems': elems, 'pos': pos_h2o, 'bonds': bonds, 'nnode': nnode},
        {'elems': elems, 'pos': pos_h2o + shift, 'bonds': bonds, 'nnode': nnode},
    ]

    packed = pack_molecules_contiguous(mols, group_size=group_size, nodes_first=True, pad_to_group=True)
    natoms = int(packed['natoms_total'])
    pos = packed['pos']
    elems_all = packed['elems']
    is_pad = packed['is_padding']
    real = ~is_pad

    # masses / inv masses
    elems_real = [e for e in elems_all if e != 'X']
    m_real = masses_from_elems(elems_real)
    m = np.zeros((natoms,), dtype=np.float32)
    m[real] = m_real
    invm = np.zeros_like(m)
    invm[real] = 1.0 / m[real]

    # neighs/exclusions
    neighs, _ = build_neighs_bk_from_bonds(natoms, packed['bonds'], max_deg=4)
    excl1, excl2 = make_exclusions_1st_2nd(neighs)

    nnode_per_group = int(packed['nnode_group'][0])
    bkSlots = make_bk_slots_clustered(neighs, group_size=group_size, nnode_per_group=nnode_per_group, natoms=natoms)

    # ports
    port_local_atoms, K_atoms = make_ports_from_neighs(pos, neighs, K=200.0)

    # perturb state to induce corrections while keeping system symmetric-ish
    quat = np.zeros((natoms, 4), dtype=np.float32)
    quat[:, 3] = 1.0
    rng = np.random.default_rng(int(args.seed))
    pos_p, quat_p = perturb_state(pos, quat, pos_scale=float(args.pos_scale), rot_scale=float(args.rot_scale), rng=rng)

    # radius
    rad = np.zeros((natoms,), dtype=np.float32)
    rad[real] = 1.0

    sim = RRsp3(natoms, group_size=group_size, prefer_gpu=True)
    sim.upload_state(pos_p, invm, quat=quat_p)
    fixmask = np.zeros((natoms,), dtype=np.int32)
    fixmask[is_pad] |= (1 | 2 | 4)
    sim.upload_fixmask(fixmask)
    sim.upload_radius(rad)
    sim.upload_neighs_and_exclusions(neighs, excl1, excl2)
    sim.upload_cluster_ports(port_local_atoms, K_atoms, nnode_per_group=nnode_per_group)
    sim.upload_bkSlots(bkSlots)

    relax = 0.5
    epsP = 1e-5
    epsL = 1e-5

    elems_real_xyz = [e for e in elems_all if e != 'X']
    fxyz = None
    if str(args.traj):
        fxyz = open(str(args.traj), 'w')

    for istep in range(int(args.steps)):
        pos4, quat4 = sim.download_pos_quat()
        pos0 = pos4[:, :3].copy()

        dpos_coll, dpos_node, drot_node, dpos_neigh = sim.compute_cluster_deltas(nnode_per_group=nnode_per_group, dt=dt, k_coll=50.0, bbox_margin=0.5)

        # Reconstruct per-atom dx_port by gathering node/cap contributions via bkSlots
        dx_port = np.zeros((natoms, 3), dtype=np.float32)

        for ig in range(natoms // group_size):
            abase = ig * group_size
            inode_base = ig * nnode_per_group
            for il in range(nnode_per_group):
                ia = abase + il
                inode = inode_base + il
                dx_port[ia, :] += dpos_node[inode, :3]

        for ia in range(natoms):
            for k in range(4):
                slot = int(bkSlots[ia, k])
                if slot >= 0:
                    dx_port[ia, :] += dpos_neigh[slot, :3]

        dx = (dpos_coll[:, :3] + dx_port) * float(relax)
        mdx = dx * m[:, None]
        dP = np.sum(mdx, axis=0)

        com = np.sum(pos0[real] * m[real, None], axis=0) / (np.sum(m[real]) + 1e-12)
        r = pos0[real] - com[None, :]
        dL_trans = np.sum(np.cross(r, mdx[real]), axis=0)

        I = 0.4 * m
        dtheta = np.zeros((natoms, 3), dtype=np.float32)
        for ig in range(natoms // group_size):
            abase = ig * group_size
            inode_base = ig * nnode_per_group
            for il in range(nnode_per_group):
                ia = abase + il
                inode = inode_base + il
                dtheta[ia, :] = drot_node[inode, :3] * float(relax)
        dL_rot = np.sum((I[:, None] * dtheta), axis=0)
        dL = dL_trans + dL_rot

        nP = float(np.linalg.norm(dP))
        nL = float(np.linalg.norm(dL))
        print(f"step {istep:4d} |sum(m*dx)|={nP:.3e} |sum(rxm*dx)+sum(I*dtheta)|={nL:.3e}")
        if nP > epsP or nL > epsL:
            raise RuntimeError(f"PBD correction conservation failed at step={istep}: |sum(m*dx)|={nP:.3e} |sum(rxm*dx)+sum(I*dtheta)|={nL:.3e}")

        if (fxyz is not None) and ((istep % int(args.traj_every)) == 0):
            comment = f"step={istep} sum_m_dx=({dP[0]:.3e},{dP[1]:.3e},{dP[2]:.3e}) |sum_m_dx|={nP:.3e} sum_Lcorr=({dL[0]:.3e},{dL[1]:.3e},{dL[2]:.3e}) |sum_Lcorr|={nL:.3e}"
            pos_xyz = pos0[real]
            if len(elems_real_xyz) != pos_xyz.shape[0]:
                raise RuntimeError(f"XYZ export mismatch: len(elems_real_xyz)={len(elems_real_xyz)} pos_xyz.shape={pos_xyz.shape}")
            write_xyz_frame(fxyz, elems_real_xyz, pos_xyz, comment=comment)

        sim.apply_cluster_corrections(nnode_per_group=nnode_per_group, relaxation=relax)

    if fxyz is not None:
        fxyz.close()

    print("[PASS] PBD correction conservation within tolerance")


if __name__ == '__main__':
    main()
