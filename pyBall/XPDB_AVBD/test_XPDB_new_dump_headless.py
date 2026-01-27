#!/usr/bin/env python3
import argparse
import numpy as np
import pyopencl as cl

from pyBall import AtomicSystem
from pyBall.OCL.MMFFL import MMFFL
from pyBall.XPDB_AVBD.XPDB_new import XPDB_new
from pyBall.XPDB_AVBD.test_TiledJacobi_molecules import _bonds_to_adj, _pack_fixed_bonds


def _pad_int(i, w=4):
    s = str(int(i))
    if len(s) >= w:
        return s
    return (' ' * (w - len(s))) + s


def dump_xpdb_gpu_state_text(path_out, *, xyz_path, fixed, num_atoms, n_max_bonded, group_size, max_ghosts,
                            pos4, pred4, params4, bboxes_min, bboxes_max, ghost_indices, ghost_counts,
                            bond_idx_global, bond_len, bond_k, bond_idx_local):
    if path_out is None:
        raise ValueError('dump_xpdb_gpu_state_text: path_out is None')

    n = int(num_atoms)
    ng = int((n + group_size - 1) // group_size)
    ffmt = '{:.' + str(int(fixed)) + 'f}'

    out = []
    out.append(f"# XPDB_new headless dump xyz={xyz_path}")
    out.append(f"# numAtoms={n} numGroups={ng} nMaxBonded={int(n_max_bonded)} maxGhosts={int(max_ghosts)}")

    out.append(f"# pos n={n} stride=4 cols=3")
    for i in range(n):
        out.append(f"pos[{_pad_int(i,4)}] {ffmt.format(float(pos4[i,0]))} {ffmt.format(float(pos4[i,1]))} {ffmt.format(float(pos4[i,2]))}")

    out.append(f"# pred_pos n={n} stride=4 cols=3")
    for i in range(n):
        out.append(f"pred_pos[{_pad_int(i,4)}] {ffmt.format(float(pred4[i,0]))} {ffmt.format(float(pred4[i,1]))} {ffmt.format(float(pred4[i,2]))}")

    out.append(f"# atom_params n={n} layout=vec4(radius,0,0,mass)")
    for i in range(n):
        out.append(f"atom_params[{_pad_int(i,4)}] r={ffmt.format(float(params4[i,0]))} m={ffmt.format(float(params4[i,3]))}")

    out.append(f"# bboxes n={ng*2} stride=4 cols=3")
    for g in range(ng):
        out.append(f"bboxes[{_pad_int(g*2+0,4)}] {ffmt.format(float(bboxes_min[g,0]))} {ffmt.format(float(bboxes_min[g,1]))} {ffmt.format(float(bboxes_min[g,2]))}")
        out.append(f"bboxes[{_pad_int(g*2+1,4)}] {ffmt.format(float(bboxes_max[g,0]))} {ffmt.format(float(bboxes_max[g,1]))} {ffmt.format(float(bboxes_max[g,2]))}")

    out.append(f"# ghost_packed numGroups={ng} maxGhosts={int(max_ghosts)} layout=[idx... idx][count@maxGhosts]")
    for g in range(ng):
        c = int(ghost_counts[g])
        out.append(f"ghost_packed[grp={_pad_int(g,3)}] count={c}")
        base = g * max_ghosts
        for k in range(c):
            out.append(f"ghost_packed[grp={_pad_int(g,3)}][{_pad_int(k,3)}] {int(ghost_indices[base + k])}")

    out.append(f"# bond_idx_global nAtoms={n} nMaxBonded={int(n_max_bonded)}")
    for i in range(n):
        row = ' '.join(str(int(bond_idx_global[i,k])) for k in range(n_max_bonded))
        out.append(f"bond_idx_global[{_pad_int(i,4)}] {row}")

    out.append(f"# bond_idx_local nAtoms={n} nMaxBonded={int(n_max_bonded)}")
    for i in range(n):
        row = ' '.join(str(int(bond_idx_local[i,k])) for k in range(n_max_bonded))
        out.append(f"bond_idx_local[{_pad_int(i,4)}] {row}")

    out.append(f"# bond_len_stiff nAtoms={n} nMaxBonded={int(n_max_bonded)} layout=(L,K) per slot")
    for i in range(n):
        parts = []
        for k in range(n_max_bonded):
            parts.append(f"{ffmt.format(float(bond_len[i,k]))},{ffmt.format(float(bond_k[i,k]))}")
        out.append(f"bond_len_stiff[{_pad_int(i,4)}] {' '.join(parts)}")

    out.append('')
    with open(path_out, 'w') as f:
        f.write('\n'.join(out))


def main():
    p = argparse.ArgumentParser(description='Headless XPDB_new(OpenCL) dump of build_local_topology outputs')
    p.add_argument('--molecule', type=str, required=True)
    p.add_argument('--type_source', type=str, default='table')

    p.add_argument('--enable_angles', type=int, default=1)
    p.add_argument('--add_pi', type=int, default=1)
    p.add_argument('--two_pi', type=int, default=1)
    p.add_argument('--add_pi_align', type=int, default=1)
    p.add_argument('--add_epair', type=int, default=1)
    p.add_argument('--add_epair_pairs', type=int, default=1)
    p.add_argument('--align_pi_vectors', type=int, default=0)

    p.add_argument('--L_pi', type=float, default=1.0)
    p.add_argument('--L_epair', type=float, default=0.5)
    p.add_argument('--k_angle', type=float, default=100.0)
    p.add_argument('--k_pi', type=float, default=50.0)
    p.add_argument('--k_pi_orth', type=float, default=30.0)
    p.add_argument('--k_pi_align', type=float, default=15.0)
    p.add_argument('--k_ep', type=float, default=40.0)
    p.add_argument('--k_ep_orth', type=float, default=25.0)
    p.add_argument('--k_ep_pair', type=float, default=10.0)

    p.add_argument('--bond_k', type=float, default=200.0)
    p.add_argument('--bond_len', type=float, default=1.3)
    p.add_argument('--atom_rad', type=float, default=0.2)
    p.add_argument('--atom_mass', type=float, default=1.0)
    p.add_argument('--max_bonded', type=int, default=16)

    p.add_argument('--dump_out', type=str, required=True)
    p.add_argument('--dump_fixed', type=int, default=6)

    p.add_argument('--group_size', type=int, default=64)
    p.add_argument('--max_ghosts', type=int, default=128)
    p.add_argument('--Rmax', type=float, default=None)
    p.add_argument('--coll_scale', type=float, default=2.0)
    p.add_argument('--bbox_scale', type=float, default=2.0)

    args = p.parse_args()

    mol = AtomicSystem.AtomicSystem(fname=args.molecule)
    if mol.bonds is None or len(mol.bonds) == 0:
        mol.findBonds()
    mol.neighs()

    ff = MMFFL(
        L_pi=args.L_pi,
        two_pi_dummies=bool(args.two_pi),
        Kang=args.k_angle,
        Kpi_host=args.k_pi,
        Kpi_orth=args.k_pi_orth,
        Kpi_align=args.k_pi_align,
        L_epair=args.L_epair,
        Kep_host=args.k_ep,
        Kep_orth=args.k_ep_orth,
        Kep_pair=args.k_ep_pair,
        align_pi_vectors=bool(args.align_pi_vectors),
    )

    topo = ff.build_topology(
        mol,
        type_source=args.type_source,
        add_angle=bool(args.enable_angles),
        add_pi=bool(args.add_pi),
        two_pi_dummies=bool(args.two_pi),
        add_pi_align=bool(args.add_pi_align),
        add_epair=bool(args.add_epair),
        add_epair_pairs=bool(args.add_epair_pairs),
    )

    apos_all = np.asarray(ff.apos[:ff.natoms, :3], dtype=np.float32)
    bonds_primary = topo.get('bonds_primary', [])
    bonds_derived = topo.get('bonds_angle', []) + topo.get('bonds_pi', []) + topo.get('bonds_pi_align', []) + topo.get('bonds_epair', []) + topo.get('bonds_epair_pair', [])
    bonds_all = sorted(set(bonds_primary + bonds_derived))

    n_all = int(apos_all.shape[0])
    bonds_adj = _bonds_to_adj(n_all, bonds_all, default_L=None, default_K=args.bond_k, apos=apos_all)
    bond_idx, bond_len, bond_k, _ = _pack_fixed_bonds(n_all, bonds_adj, n_max_bonded=args.max_bonded)

    radius = np.full(n_all, args.atom_rad, dtype=np.float32)
    mass = np.full(n_all, args.atom_mass, dtype=np.float32)

    sim = XPDB_new(n_all, group_size=args.group_size)
    sim.max_ghosts = int(args.max_ghosts)

    vel = np.zeros((n_all, 3), dtype=np.float32)
    sim.upload_data(pos=apos_all, vel=vel, radius=radius, mass=mass)
    sim.upload_bonds_fixed_from_arrays(bond_idx, bond_len, bond_k)

    Rmax = float(args.Rmax) if args.Rmax is not None else float(np.max(radius))
    sim.build_local_topology(Rmax=Rmax, coll_scale=float(args.coll_scale), bbox_scale=float(args.bbox_scale))

    # download buffers
    bmin = np.empty((sim.num_groups, 4), dtype=np.float32)
    bmax = np.empty((sim.num_groups, 4), dtype=np.float32)
    cl.enqueue_copy(sim.queue, bmin, sim.cl_bboxes_min).wait()
    cl.enqueue_copy(sim.queue, bmax, sim.cl_bboxes_max).wait()

    ghost_idx = np.empty(sim.num_groups * sim.max_ghosts, dtype=np.int32)
    ghost_cnt = np.empty(sim.num_groups, dtype=np.int32)
    cl.enqueue_copy(sim.queue, ghost_idx, sim.cl_ghost_indices).wait()
    cl.enqueue_copy(sim.queue, ghost_cnt, sim.cl_ghost_counts).wait()

    bond_loc = np.empty((n_all, sim.n_max_bonded), dtype=np.int32)
    cl.enqueue_copy(sim.queue, bond_loc, sim.cl_bond_indices_local).wait()

    pos4 = np.zeros((n_all, 4), dtype=np.float32)
    pred4 = np.zeros((n_all, 4), dtype=np.float32)
    params4 = np.zeros((n_all, 4), dtype=np.float32)
    cl.enqueue_copy(sim.queue, pos4, sim.cl_pos).wait()
    cl.enqueue_copy(sim.queue, pred4, sim.cl_pred_pos).wait()
    cl.enqueue_copy(sim.queue, params4, sim.cl_params).wait()

    dump_xpdb_gpu_state_text(
        args.dump_out,
        xyz_path=args.molecule,
        fixed=args.dump_fixed,
        num_atoms=n_all,
        n_max_bonded=sim.n_max_bonded,
        group_size=sim.group_size,
        max_ghosts=sim.max_ghosts,
        pos4=pos4,
        pred4=pred4,
        params4=params4,
        bboxes_min=bmin,
        bboxes_max=bmax,
        ghost_indices=ghost_idx,
        ghost_counts=ghost_cnt,
        bond_idx_global=bond_idx,
        bond_len=bond_len,
        bond_k=bond_k,
        bond_idx_local=bond_loc,
    )
    print(f"[DUMP] wrote {args.dump_out}")


if __name__ == '__main__':
    main()
