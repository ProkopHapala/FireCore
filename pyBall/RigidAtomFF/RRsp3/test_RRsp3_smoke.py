import sys, os
import numpy as np

_THIS_DIR = os.path.abspath(os.path.dirname(__file__))
_SHARED_DIR = os.path.abspath(os.path.join(_THIS_DIR, '..', 'shared'))
for _p in (_THIS_DIR, _SHARED_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from RRsp3 import RRsp3, build_neighs_bk_from_bonds, make_bk_slots_clustered, make_exclusions_1st_2nd
from XPTB_utils import pack_molecules_contiguous, make_h2o_geometry, masses_from_elems
from RRsp3_utils import make_ports_from_neighs

def check_local_ranges(neighs_local, excl1_local, excl2_local, ghost_counts, *, group_size=64):
    neighs_local = np.asarray(neighs_local, dtype=np.int32)
    excl1_local = np.asarray(excl1_local, dtype=np.int32)
    excl2_local = np.asarray(excl2_local, dtype=np.int32)
    ghost_counts = np.asarray(ghost_counts, dtype=np.int32)
    natoms = int(neighs_local.shape[0])
    if neighs_local.shape != (natoms, 4):
        raise ValueError(f"check_local_ranges: neighs_local.shape={neighs_local.shape} expected ({natoms},4)")

    ng = natoms // int(group_size)
    for ig in range(ng):
        g = int(ghost_counts[ig])
        hi = int(group_size) + g
        abase = ig * int(group_size)
        arrs = (neighs_local[abase:abase+group_size], excl1_local[abase:abase+group_size], excl2_local[abase:abase+group_size])
        for A in arrs:
            bad = (A >= hi) & (A != -1)
            if np.any(bad):
                w = np.where(bad)
                raise RuntimeError(f"local index out of range in group {ig}: hi={hi} example atom={w[0][0]} k={w[1][0]} val={A[w[0][0], w[1][0]]}")


def main():
    group_size = 64

    pos_h2o, _bonds_adj = make_h2o_geometry(add_angle=False)
    elems = ['O', 'H', 'H']
    bonds = [(0, 1), (0, 2)]
    nnode = 1

    mols = [
        {'elems': elems, 'pos': pos_h2o, 'bonds': bonds, 'nnode': nnode},
        {'elems': elems, 'pos': pos_h2o + np.array([4.0, 0.0, 0.0], dtype=np.float32), 'bonds': bonds, 'nnode': nnode},
    ]

    packed = pack_molecules_contiguous(mols, group_size=group_size, nodes_first=True, pad_to_group=True)
    natoms = int(packed['natoms_total'])
    ng = natoms // group_size
    if ng != 2:
        raise RuntimeError(f"expected 2 groups, got ng={ng} natoms={natoms}")

    pos = packed['pos']
    elems_all = packed['elems']
    is_pad = packed['is_padding']

    # masses / inv masses
    elems_real = [e for e in elems_all if e != 'X']
    m_real = masses_from_elems(elems_real)
    m = np.ones((natoms,), dtype=np.float32)
    real = ~is_pad
    m[real] = m_real
    invm = np.zeros_like(m)
    invm[real] = 1.0 / m[real]

    # neighs on packed indexing
    neighs, _bks = build_neighs_bk_from_bonds(natoms, packed['bonds'], max_deg=4)

    # cluster invariant for this smoke test: constant nnode per group
    nnode_per_group = int(packed['nnode_group'][0])
    if not np.all(packed['nnode_group'] == nnode_per_group):
        raise RuntimeError(f"smoke test assumes constant nnode_per_group across groups, got {packed['nnode_group']}")

    # exclusions
    excl1, excl2 = make_exclusions_1st_2nd(neighs)

    # ports + stiffness from neighs
    port_local_atoms, K_atoms = make_ports_from_neighs(pos, neighs, K=200.0)

    # radius
    rad = np.zeros((natoms,), dtype=np.float32)
    rad[real] = 1.0

    sim = RRsp3(natoms, group_size=group_size, prefer_gpu=True)

    quat = np.zeros((natoms, 4), dtype=np.float32)
    quat[:, 3] = 1.0
    sim.upload_state(pos, invm, quat=quat)
    sim.upload_radius(rad)
    sim.upload_neighs_and_exclusions(neighs, excl1, excl2)

    sim.upload_cluster_ports(port_local_atoms, K_atoms, nnode_per_group=nnode_per_group)

    bkSlots = make_bk_slots_clustered(neighs, group_size=group_size, nnode_per_group=nnode_per_group, natoms=natoms)
    sim.upload_bkSlots(bkSlots)

    sim.step_cluster(nnode_per_group=nnode_per_group, dt=0.1, k_coll=50.0, relaxation=0.5, bbox_margin=0.5)

    ghost_counts = sim.download(sim.cl_ghost_counts, (ng,), np.int32)
    ghost_indices = sim.download(sim.cl_ghost_indices, (ng * sim.max_ghosts,), np.int32)
    neighs_local = sim.download(sim.cl_neighs_local, (natoms, 4), np.int32)
    excl1_local = sim.download(sim.cl_excl1_local, (natoms, 4), np.int32)
    excl2_local = sim.download(sim.cl_excl2_local, (natoms, 4), np.int32)
    dpos_coll = sim.download(sim.cl_dpos_coll, (natoms, 4), np.float32)
    dpos_node = sim.download(sim.cl_dpos_node, (sim._nnode_tot, 4), np.float32)
    dpos_neigh = sim.download(sim.cl_dpos_neigh, (sim._nnode_tot * 4, 4), np.float32)

    check_local_ranges(neighs_local, excl1_local, excl2_local, ghost_counts, group_size=group_size)

    print("ghost_counts", ghost_counts.tolist())
    for ig in range(ng):
        g = int(ghost_counts[ig])
        print(f"group {ig} ghosts {g} indices (first {min(g,8)}):", ghost_indices[ig * sim.max_ghosts:ig * sim.max_ghosts + min(g, 8)].tolist())

    print("neighs_local group0 atoms0..2:")
    print(neighs_local[0:3])
    print("excl1_local group0 atoms0..2:")
    print(excl1_local[0:3])
    print("excl2_local group0 atoms0..2:")
    print(excl2_local[0:3])

    print("dpos_coll group0 atoms0..2:")
    print(dpos_coll[0:3])

    print("dpos_node (nodes) first 2:")
    print(dpos_node[:2])

    print("dpos_neigh (node recoil slots) first 8:")
    print(dpos_neigh[:8])

    pos4, quat4 = sim.download_pos_quat()
    print("pos after step group0 atoms0..2:")
    print(pos4[0:3])
    print("quat after step group0 atoms0..0:")
    print(quat4[0])


if __name__ == "__main__":
    main()
