#!/usr/bin/env python3
import argparse
import math
import os
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
                            bond_idx_global, bond_len, bond_k, bond_idx_local,
                            force_bond=None, force_coll=None):
    if path_out is None:
        raise ValueError('dump_xpdb_gpu_state_text: path_out is None')

    n = int(num_atoms)
    ng = int((n + group_size - 1) // group_size)
    ffmt = '{:.' + str(int(fixed)) + 'f}'

    out = []
    out.append(f"# XPDB_WebGPU headless xyz={os.path.abspath(xyz_path)} numAtoms={n} numGroups={ng} nMaxBonded={int(n_max_bonded)} maxGhosts={int(max_ghosts)}")

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

    if force_bond is not None:
        out.append(f"# force_bond n={n} stride=4 cols=3")
        for i in range(n):
            fx, fy, fz = force_bond[i]
            out.append(f"force_bond[{_pad_int(i,4)}] {ffmt.format(float(fx))} {ffmt.format(float(fy))} {ffmt.format(float(fz))}")

    if force_coll is not None:
        out.append(f"# force_coll n={n} stride=4 cols=3")
        for i in range(n):
            fx, fy, fz = force_coll[i]
            out.append(f"force_coll[{_pad_int(i,4)}] {ffmt.format(float(fx))} {ffmt.format(float(fy))} {ffmt.format(float(fz))}")

    out.append('')
    with open(path_out, 'w') as f:
        f.write('\n'.join(out))


def write_xyz_frame(f, symbols, pos4, *, step=0, dt=0.01, title='XPDB'):
    n = int(len(symbols))
    f.write(f"{n}\n")
    f.write(f"{title} step={int(step)} time={float(step)*float(dt):.6f}\n")
    for i, s in enumerate(symbols):
        x = float(pos4[i, 0]); y = float(pos4[i, 1]); z = float(pos4[i, 2])
        f.write(f"{s} {x:.8f} {y:.8f} {z:.8f}\n")


def _deterministic_noise(idx, axis, seed):
    # Simple portable hash -> [-1,1] using sin fract noise; matches JS implementation.
    x = math.sin(float(seed) * 12.9898 + float(idx) * 78.233 + float(axis) * 37.719) * 43758.5453
    frac = x - math.floor(x)
    return (frac * 2.0) - 1.0


def apply_initial_distortion(apos, scale=1.0, noise_amp=0.0, seed=1337):
    if scale != 1.0:
        apos[:, :3] *= float(scale)
    if noise_amp != 0.0:
        amp = float(noise_amp)
        for idx in range(apos.shape[0]):
            for axis in range(3):
                apos[idx, axis] += amp * _deterministic_noise(idx, axis, seed)


def main():
    p = argparse.ArgumentParser(description='Headless XPDB_new(OpenCL) dump of build_local_topology + solve outputs')
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
    p.add_argument('--dt', type=float, default=0.01)
    p.add_argument('--inner_iters', type=int, default=0)
    p.add_argument('--k_coll', type=float, default=200.0)
    p.add_argument('--omega', type=float, default=0.0)
    p.add_argument('--momentum_beta', type=float, default=0.0)

    p.add_argument('--n_steps', type=int, default=0)
    p.add_argument('--traj_xyz', type=str, default=None)
    p.add_argument('--traj_real_only', type=int, default=1)
    p.add_argument('--init_scale', type=float, default=1.0, help='Uniform scale applied to initial coordinates before upload')
    p.add_argument('--init_noise', type=float, default=0.0, help='Additive jitter amplitude (Angstrom) applied per-component before upload')
    p.add_argument('--init_seed', type=int, default=1337, help='Seed for deterministic initial jitter (shared with WebGPU harness)')

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

    apos_ref = np.array(ff.apos[:ff.natoms, :3], dtype=np.float32, copy=True)
    apos_all = apos_ref.copy()
    apply_initial_distortion(apos_all, scale=float(args.init_scale), noise_amp=float(args.init_noise), seed=int(args.init_seed))
    bonds_primary = topo.get('bonds_primary', [])
    bonds_derived = topo.get('bonds_angle', []) + topo.get('bonds_pi', []) + topo.get('bonds_pi_align', []) + topo.get('bonds_epair', []) + topo.get('bonds_epair_pair', [])
    bonds_all = sorted(set(bonds_primary + bonds_derived))

    n_all = int(apos_all.shape[0])
    bonds_adj = _bonds_to_adj(n_all, bonds_all, default_L=None, default_K=args.bond_k, apos=apos_ref)
    bond_idx, bond_len, bond_k, _ = _pack_fixed_bonds(n_all, bonds_adj, n_max_bonded=args.max_bonded)

    radius = np.full(n_all, args.atom_rad, dtype=np.float32)
    mass = np.full(n_all, args.atom_mass, dtype=np.float32)

    sim = XPDB_new(n_all, group_size=args.group_size)
    sim.max_ghosts = int(args.max_ghosts)

    vel = np.zeros((n_all, 3), dtype=np.float32)
    sim.upload_data(pos=apos_all, vel=vel, radius=radius, mass=mass)
    sim.upload_bonds_fixed_from_arrays(bond_idx, bond_len, bond_k)

    Rmax = float(args.Rmax) if args.Rmax is not None else float(np.max(radius))

    # Optional trajectory run: rebuild local topology each step to match WebGPU path.
    if int(args.n_steps) > 0:
        if not args.traj_xyz:
            raise ValueError('--n_steps > 0 requires --traj_xyz PATH')

        type_names = topo.get('type_names', None)
        if type_names is None or len(type_names) != n_all:
            raise ValueError('Missing topo.type_names for trajectory symbol export')

        # Use element-ish symbols for visualization.
        # For many molecules types are simple (C,O,H), and dummies are Pi/E.
        symbols_all = []
        for t in type_names:
            s = str(t)
            if s in ('Pi', 'E'):
                symbols_all.append(s)
            else:
                s2 = s.split('_')[0]
                if len(s2) > 2:
                    s2 = s2[0]
                symbols_all.append(s2)

        n_real = int(topo.get('n_real', 0))
        if int(args.traj_real_only) != 0:
            symbols = symbols_all[:n_real]
        else:
            symbols = symbols_all

        pos4 = np.zeros((n_all, 4), dtype=np.float32)
        with open(args.traj_xyz, 'w') as ftraj:
            for istep in range(int(args.n_steps)):
                sim.build_local_topology(Rmax=Rmax, coll_scale=float(args.coll_scale), bbox_scale=float(args.bbox_scale))
                # Match WebGPU step(): pred_pos := curr_pos before solve
                cl.enqueue_copy(sim.queue, sim.cl_pred_pos, sim.cl_pos)
                if args.inner_iters > 0:
                    sim.solve_cluster_jacobi(
                        dt=float(args.dt),
                        iterations=int(args.inner_iters),
                        k_coll=float(args.k_coll),
                        omega=float(args.omega),
                        momentum_beta=float(args.momentum_beta),
                    )
                cl.enqueue_copy(sim.queue, pos4, sim.cl_pos).wait()
                if int(args.traj_real_only) != 0:
                    write_xyz_frame(ftraj, symbols, pos4[:len(symbols)], step=istep, dt=float(args.dt), title='XPDB_OCL')
                else:
                    write_xyz_frame(ftraj, symbols, pos4, step=istep, dt=float(args.dt), title='XPDB_OCL')

        # After trajectory run, continue with final-state dump below.

    # One-shot topology build (or final refresh before dump)
    sim.build_local_topology(Rmax=Rmax, coll_scale=float(args.coll_scale), bbox_scale=float(args.bbox_scale))
    if args.inner_iters > 0 and int(args.n_steps) == 0:
        sim.solve_cluster_jacobi(
            dt=float(args.dt),
            iterations=int(args.inner_iters),
            k_coll=float(args.k_coll),
            omega=float(args.omega),
            momentum_beta=float(args.momentum_beta),
        )

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

    f_bond = np.zeros((n_all, 3), dtype=np.float32)
    f_coll = np.zeros((n_all, 3), dtype=np.float32)
    if args.inner_iters > 0:
        fb, fc = sim.get_debug_forces()
        f_bond[:, :] = fb
        f_coll[:, :] = fc

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
        force_bond=f_bond if args.inner_iters > 0 else None,
        force_coll=f_coll if args.inner_iters > 0 else None,
    )
    print(f"[DUMP] wrote {args.dump_out}")


if __name__ == '__main__':
    main()
