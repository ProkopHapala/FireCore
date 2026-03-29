#!/usr/bin/env python3
import os, sys, time, argparse
import numpy as np
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(BASE, os.pardir, os.pardir))
if ROOT not in sys.path:
    sys.path.append(ROOT)
if BASE not in sys.path:
    sys.path.append(BASE)

from pyBall.AtomicSystem import AtomicSystem
from pyBall.OCL.MMFF import MMFF as MMFF_pyocl
from pyBall.OCL.MolecularDynamics import MolecularDynamics
import pyopencl as cl

np.set_printoptions(linewidth=200, suppress=True)


def write_xyz_frame(fout, enames, pos3, comment=""):
    nat = pos3.shape[0]
    fout.write(f"{nat}\n")
    fout.write(f"{comment}\n")
    for i in range(nat):
        fout.write(f"{enames[i]:3s} {pos3[i,0]:12.6f} {pos3[i,1]:12.6f} {pos3[i,2]:12.6f}\n")


def download(md, natoms, nvecs):
    apos_buf = np.empty((nvecs, 4), dtype=np.float32)
    af_buf   = np.empty((nvecs, 4), dtype=np.float32)
    cl.enqueue_copy(md.queue, apos_buf, md.buffer_dict['apos'])
    cl.enqueue_copy(md.queue, af_buf,   md.buffer_dict['aforce'])
    md.queue.finish()
    return apos_buf[:natoms, :3].copy(), af_buf[:natoms, :3].copy(), af_buf[:natoms, 3].copy()


def make_h2o():
    mol = AtomicSystem()
    mol.enames = np.array(["O", "H", "H"], dtype=object)
    mol.apos = np.array([
        [0.0, 0.0, 0.0],
        [0.757, 0.586, 0.0],
        [-0.757, 0.586, 0.0]
    ], dtype=np.float32)
    mol.qs = np.array([-0.82, +0.41, +0.41], dtype=np.float32)
    mol.natoms = 3
    mol.atypes = [8, 1, 1]
    mol.bonds = [(0, 1), (0, 2)]
    mol.neighs()
    mol.atom_types_mmff = ["O_3", "H_", "H_"]
    return mol


def bonds_from_mol(mol):
    b = getattr(mol, 'bonds', None)
    if b is None:
        return []
    out = []
    for ij in b:
        i, j = int(ij[0]), int(ij[1])
        if i == j:
            continue
        if i < 0 or j < 0:
            continue
        out.append((i, j))
    return out


def eval_bond_stats(pos, bonds, ref=None):
    if not bonds:
        return None
    ls = np.array([np.linalg.norm(pos[i] - pos[j]) for (i, j) in bonds], dtype=float)
    out = {
        'l_min': float(ls.min()),
        'l_max': float(ls.max()),
        'l_mean': float(ls.mean()),
    }
    if ref is not None:
        l0 = np.array([np.linalg.norm(ref[i] - ref[j]) for (i, j) in bonds], dtype=float)
        rat = ls / (l0 + 1e-12)
        out.update({
            'stretch_max': float(rat.max()),
            'stretch_mean': float(rat.mean()),
        })
    return out


def relax_at_tip(md, mm, enames, tip_pos, anchor_idx, dt, damp, Fconv, nstep_max, out_stride=50, bonds=None, bonds_ref=None, verbose=True):
    nat = mm.natoms

    # Set constraint for anchor atom
    constr  = np.zeros((nat, 4), dtype=np.float32)
    constrK = np.zeros((nat, 4), dtype=np.float32)
    constr[anchor_idx, :3] = np.asarray(tip_pos, dtype=np.float32)
    constr[anchor_idx, 3]  = 1.0
    constrK[anchor_idx, 0] = 2000.0
    constrK[anchor_idx, 1] = 2000.0
    constrK[anchor_idx, 2] = 2000.0
    md.toGPU('constr',  constr)
    md.toGPU('constrK', constrK)

    vfac = np.float32(1.0 - damp)
    md.toGPU('MDparams', np.array([dt, 1e+6, vfac, 0.0], dtype=np.float32), byte_offset=0)

    fmax_hist = []

    for istep in range(1, nstep_max + 1):
        md.run_cleanForceMMFFf4()

        md.run_getNonBond_GridFF_Bspline_ex2()
        md.run_getMMFFf4()
        md.run_updateAtomsMMFFf4()

        want = (istep % out_stride == 0) or (istep == 1) or (istep == nstep_max)
        if want:
            pos, F, _ = download(md, nat, mm.nvecs)
            fmag = np.sqrt(np.sum(F * F, axis=1))
            fmax = float(fmag.max())
            fmax_hist.append((istep, fmax))

            if not np.isfinite(pos).all():
                raise RuntimeError(f"NaN/Inf positions at relax step {istep}")

            bstat = eval_bond_stats(pos, bonds, ref=bonds_ref) if bonds is not None else None
            if verbose:
                if bstat is None:
                    print(f"    relax step {istep:6d} fmax={fmax:12.6e}  pos_min={pos.min(axis=0)} pos_max={pos.max(axis=0)}")
                else:
                    print(f"    relax step {istep:6d} fmax={fmax:12.6e}  l[min,max,mean]=({bstat['l_min']:.3f},{bstat['l_max']:.3f},{bstat['l_mean']:.3f})  stretch_max={bstat.get('stretch_max',np.nan):.3f}")

            if fmax < Fconv:
                return pos, F, istep, fmax_hist

    pos, F, _ = download(md, nat, mm.nvecs)
    return pos, F, nstep_max, fmax_hist


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mol',        type=str, default='', help='Input molecule .mol2 or .xyz. If empty, uses built-in H2O.')
    ap.add_argument('--anchor',     type=int, default=1, help='Anchor atom index constrained to tip position.')
    ap.add_argument('--grid', default=os.path.join(ROOT, 'tests/tUFF/data/NaCl_1x1_L3/Bspline_PLQd.npy'), help='GridFF Bspline_PLQd.npy path (set empty to run vacuum like debugging runs)')
    ap.add_argument('--grid_p0',    type=float, nargs=3, default=[-2.0, -2.0, 0.0])
    ap.add_argument('--grid_step',  type=float, nargs=3, default=[0.1, 0.1, 0.1])
    ap.add_argument('--grid_alpha', type=float, default=2.2)
    ap.add_argument('--z_start',    type=float, default=4.0)

    ap.add_argument('--dx',         type=float, default=0.1, help='Tip step in x (A)')
    ap.add_argument('--nscan',      type=int,   default=40, help='Number of relax points; total path = dx*nscan')

    ap.add_argument('--dt',         type=float, default=0.02)
    ap.add_argument('--damp',       type=float, default=0.01)
    ap.add_argument('--Fconv',      type=float, default=1e-4)
    ap.add_argument('--nstep_max',  type=int,   default=2000)
    ap.add_argument('--zero_vel',   type=int,   default=1, help='If 1, initialize velocities to zero (recommended for deterministic relaxation). If 0, use random velocities.')
    ap.add_argument('--cold_start', type=int,   default=0, help='If 1, reset geometry (and optionally velocities) at each scan point to the same rigid-shifted initial configuration. This makes the scan directly comparable to PathOpt band-off replicas.')
    ap.add_argument('--out',        type=str,   default='out_relaxed_scan')
    ap.add_argument('--stride',     type=int,   default=200, help='Print relax stats every stride steps')

    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)

    if args.mol:
        mol = AtomicSystem(args.mol)
    else:
        mol = make_h2o()

    if getattr(mol, 'enames', None) is None:
        raise RuntimeError('Molecule has no enames')
    if getattr(mol, 'apos', None) is None:
        raise RuntimeError('Molecule has no apos')

    mol.enames = np.array([e.split('_')[0] if isinstance(e, str) else str(e) for e in mol.enames], dtype=object)

    mol.apos = np.asarray(mol.apos, dtype=np.float32)
    mol.apos[:, 2] += np.float32(args.z_start)

    anchor = int(args.anchor)
    if anchor < 0 or anchor >= mol.apos.shape[0]:
        raise ValueError(f"anchor out of range: {anchor}")

    bonds = bonds_from_mol(mol)
    bonds_ref = mol.apos.copy()

    # Build MMFF (same style as test_ditetraceno_surface.py, but with no PBC for molecule-molecule)
    mm = MMFF_pyocl(bTorsion=False, verbosity=0)
    mm.capping_atoms = set()
    mm.reorder_nodes_first = False
    mm.nPBC = (0, 0, 0)
    mm.lvec = (np.eye(3, dtype=np.float32) * 100.0)
    mm.toMMFFsp3_loc(mol=mol, atom_types=mm.atom_types, bRealloc=True, bEPairs=False, bUFF=False)
    mm.make_back_neighs(b_cap_neighs=False)
    mm.excl = mm._make_excl_1_2_3(mm.neighs, neighCell=mm.neighCell, npbc=mm.npbc, EXCL_MAX=16)

    print(f"MMFF built: natoms={mm.natoms} nnode={mm.nnode} nvecs={mm.nvecs} npbc={mm.npbc} nPBC={mm.nPBC}")

    # Setup GPU
    md = MolecularDynamics(nloc=32, debug_build_options='-DDBG_UFF=0', enable_nonbond=True)
    md.realloc(mm, nSystems=1)
    md.setup_kernels()
    md.pack_system(0, mm)

    grid_ok = False
    V4 = None
    if args.grid:
        grid_path = args.grid
        if not os.path.isabs(grid_path):
            grid_path = os.path.abspath(os.path.join(os.getcwd(), grid_path))
        if not os.path.exists(grid_path):
            print(f"[WARN] Grid file not found: {grid_path}. Running in vacuum (no GridFF).")
        else:
            V = np.load(grid_path)
            if V.ndim != 4:
                raise ValueError(f"Grid Bspline must be 4D [nx,ny,nz,nch]; got {V.shape}")
            if V.shape[3] == 3:
                V4 = np.zeros((V.shape[0], V.shape[1], V.shape[2], 4), dtype=np.float32)
                V4[:, :, :, :3] = V.astype(np.float32)
            elif V.shape[3] == 4:
                V4 = V.astype(np.float32)
            else:
                raise ValueError(f"Grid Bspline channels must be 3 or 4; got {V.shape}")
            grid_ok = True

    if grid_ok:
        md.initGridFF(grid_shape=tuple(int(x) for x in V4.shape[:3]), bspline_data=np.ascontiguousarray(V4), grid_p0=tuple(args.grid_p0), grid_step=tuple(args.grid_step), use_texture=False, r_damp=0.0, alpha_morse=float(args.grid_alpha), bKernels=True)
        if getattr(md, 'kernel_args_getNonBond_GridFF_Bspline_ex2', None) is None:
            raise RuntimeError('kernel_args_getNonBond_GridFF_Bspline_ex2 is None; ensure relax_multi.cl provides getNonBond_GridFF_Bspline_ex2 and exclusions buffer exists')

    if int(args.zero_vel) != 0:
        avel = np.zeros((mm.nvecs, 4), dtype=np.float32)
        md.toGPU('avel', avel)
    else:
        np.random.seed(42)
        avel = np.random.normal(0, 0.1, (mm.nvecs, 4)).astype(np.float32)
        md.toGPU('avel', avel)

    # Scan trajectory
    tip0 = mol.apos[anchor].copy()
    tip0[2] += 0.0

    # For cold-start scans we need an immutable reference geometry for resets
    apos0_ref = None
    if int(args.cold_start) != 0:
        apos0_ref = mol.apos.copy()

    out_xyz = os.path.join(args.out, 'scan_relaxed.xyz')
    out_npz = os.path.join(args.out, 'scan_relaxed.npz')
    out_png = os.path.join(args.out, 'scan_relaxed_xz.png')

    A_scan = np.zeros((args.nscan + 1, mm.natoms, 3), dtype=np.float32)
    tip_scan = np.zeros((args.nscan + 1, 3), dtype=np.float32)
    iters = np.zeros((args.nscan + 1,), dtype=np.int32)
    fmaxs = np.zeros((args.nscan + 1,), dtype=np.float32)

    with open(out_xyz, 'w') as fout:
        for iscan in range(args.nscan + 1):
            tip = tip0.copy()
            tip[0] += np.float32(args.dx * iscan)
            tip_scan[iscan] = tip

            if int(args.cold_start) != 0:
                # Reset positions to the same reference geometry rigidly shifted by tip displacement (like init_trajectory_shift)
                shift = tip - tip0
                pos = apos0_ref.copy()
                pos[:, :3] += shift[None, :3]
                # Upload to GPU (positions only)
                apos_buf = np.zeros((mm.nvecs, 4), dtype=np.float32)
                apos_buf[:mm.natoms, :3] = pos
                md.toGPU('apos', apos_buf)
                if int(args.zero_vel) != 0:
                    md.toGPU('avel', np.zeros((mm.nvecs, 4), dtype=np.float32))

            print(f"\n=== scan {iscan:4d}/{args.nscan} tip=({tip[0]:.3f},{tip[1]:.3f},{tip[2]:.3f}) ===")

            pos, F, niter_used, fmax_hist = relax_at_tip(
                md, mm, mol.enames, tip_pos=tip, anchor_idx=anchor,
                dt=float(args.dt), damp=float(args.damp), Fconv=float(args.Fconv), nstep_max=int(args.nstep_max),
                out_stride=int(args.stride), bonds=bonds, bonds_ref=bonds_ref, verbose=True,
            )

            fmag = np.sqrt(np.sum(F * F, axis=1))
            fmax = float(fmag.max())

            iters[iscan] = niter_used
            fmaxs[iscan] = fmax
            A_scan[iscan] = pos

            bstat = eval_bond_stats(pos, bonds, ref=bonds_ref)
            comment = f"# scan={iscan} tip=({tip[0]:.6f},{tip[1]:.6f},{tip[2]:.6f}) niter={niter_used} fmax={fmax:.6e} lmax={bstat['l_max']:.6f} stretch_max={bstat.get('stretch_max',np.nan):.3f}"
            write_xyz_frame(fout, mol.enames, pos, comment=comment)

            if (not np.isfinite(pos).all()) or (bstat is not None and bstat.get('stretch_max', 1.0) > 1.8):
                raise RuntimeError(f"Unphysical geometry detected at scan {iscan}: stretch_max={bstat.get('stretch_max',np.nan)}")

    np.savez(out_npz, A_scan=A_scan, tip_scan=tip_scan, iters=iters, fmaxs=fmaxs)

    # Plot x-z traces
    plt.figure(figsize=(6, 6))
    for ia in range(mm.natoms):
        plt.plot(A_scan[:, ia, 0], A_scan[:, ia, 2], '-', lw=1)
    plt.plot(tip_scan[:, 0], tip_scan[:, 2], 'k--', lw=2, label='tip')
    plt.xlabel('x (A)')
    plt.ylabel('z (A)')
    plt.axis('equal')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(out_png)
    plt.close()

    print("\nDONE")
    print(f"- xyz: {out_xyz}")
    print(f"- npz: {out_npz}")
    print(f"- plot: {out_png}")


if __name__ == '__main__':
    main()
