#!/usr/bin/env python3

import os, sys, argparse
import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(BASE)

from pyBall import MMFF_multi as mm_cpp
from pyBall.OCL.MMFF import MMFF as MMFF_pyocl
from pyBall.OCL.MolecularDynamics import MolecularDynamics
from pyBall.AtomicSystem import AtomicSystem

np.set_printoptions(linewidth=np.inf, suppress=False)

DATA_DIR = os.path.join(BASE, 'cpp/common_resources')
XYZ_DIR  = os.path.join(DATA_DIR, 'xyz')


def _print_arr(name, A, n=0):
    if A is None:
        print(f"{name}: None")
        return
    A = np.asarray(A)
    if (n is not None) and (n > 0) and (A.ndim >= 1) and (A.shape[0] > n):
        print(f"{name} shape={A.shape} (head {n})\n", A[:n])
    else:
        print(f"{name} shape={A.shape}\n", A)


def _absmax(A):
    A = np.asarray(A)
    return float(np.max(np.abs(A))) if A.size else 0.0


def _minmax(A):
    A = np.asarray(A)
    if A.size == 0:
        return 0.0, 0.0
    return float(np.min(A)), float(np.max(A))


def compare(name, A, B, tol=0.0):
    A = np.asarray(A)
    B = np.asarray(B)
    if A.shape != B.shape:
        print(f"{name:<10}: FAIL  shape {A.shape} != {B.shape}")
        return False
    d = B - A
    m = _absmax(d)
    ok = m <= tol
    amin, amax = _minmax(A)
    bmin, bmax = _minmax(B)
    status = 'PASS' if ok else 'FAIL'
    print(f"{name:<10}: {status}  max|Δ|={m:.3e}  tol={tol:.1e}  C++[{amin:.3e},{amax:.3e}]  PyOCL[{bmin:.3e},{bmax:.3e}]")
    if not ok:
        idx = np.unravel_index(int(np.argmax(np.abs(d))), d.shape)
        print(f"    worst idx={idx} A={A[idx]} B={B[idx]} Δ={d[idx]}")
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--xyz', type=str, default=os.path.join(XYZ_DIR, 'HCOOH.xyz'))
    ap.add_argument('--tolBuf', type=float, default=0.0)
    ap.add_argument('--tolF', type=float, default=1e-4)
    ap.add_argument('--nonbond', type=int, default=0)
    ap.add_argument('--angles', type=int, default=1)
    ap.add_argument('--pisigma', type=int, default=1)
    ap.add_argument('--pipii', type=int, default=1)
    ap.add_argument('--force-node-all', type=int, default=1, help='1: force PyOCL builder to treat all atoms as nodes (nnode=natoms) to match C++/OCL path')
    ap.add_argument('--dump', type=int, default=0, help='If >0, print full input/output buffers (use --dump-n to limit rows)')
    ap.add_argument('--dump-n', type=int, default=0, help='If >0, print only first N rows of long buffers')
    ap.add_argument('--save-npz', type=str, default='', help='If non-empty, save compared buffers/forces to NPZ at this path')
    ap.add_argument('--fast-exit', type=int, default=0, help='Exit via os._exit(0) to bypass known C++ double-free at interpreter shutdown')
    ap.add_argument('--label', type=str, default='', help='Optional label/preset name for summary printing')
    args = ap.parse_args()

    # ---- C++/OpenCL path (reference for this parity step)
    mm_cpp.init(
        nSys_=1,
        xyz_name=args.xyz,
        surf_name=None,
        smile_name=None,
        sElementTypes=os.path.join(DATA_DIR, 'ElementTypes.dat'),
        sAtomTypes=os.path.join(DATA_DIR, 'AtomTypes.dat'),
        sBondTypes=os.path.join(DATA_DIR, 'BondTypes.dat'),
        sAngleTypes=os.path.join(DATA_DIR, 'AngleTypes.dat'),
        sDihedralTypes=os.path.join(DATA_DIR, 'DihedralTypes.dat'),
        bMMFF=True,
        bEpairs=False,
        nPBC=(0,0,0),
        gridStep=0.0,
        bUFF=False,
        b141=True,
        bSimple=False,
        bConj=True,
        bCumulene=True,
    )

    mm_cpp.setSwitches2(
        CheckInvariants=-1,
        PBC=-1,
        NonBonded=(1 if args.nonbond else -1),
        NonBondNeighs=(1 if args.nonbond else -1),
        SurfAtoms=-1,
        GridFF=-1,
        MMFF=1,
        Angles=(1 if args.angles else -1),
        PiSigma=(1 if args.pisigma else -1),
        PiPiI=(1 if args.pipii else -1),
    )

    mm_cpp.getBuffs()

    # Evaluate on C++ OpenCL
    mm_cpp.eval_force_MMFF(2)
    # IMPORTANT: compare against the actual OpenCL output buffers (gpu_aforces), not stale CPU arrays
    F_cpp = np.array(mm_cpp.gpu_aforces[0, :mm_cpp.natoms, :3], copy=True, dtype=np.float64)
    E_cpp = float(np.sum(mm_cpp.gpu_aforces[0, :mm_cpp.natoms, 3]))

    if args.dump:
        n = args.dump_n
        _print_arr('CPP.neighs(int4)',   mm_cpp.neighs[:mm_cpp.natoms, :], n=n)
        _print_arr('CPP.gpu_bkNeighs',   mm_cpp.gpu_bkNeighs[0, :mm_cpp.nvecs, :], n=n)
        _print_arr('CPP.gpu_REQs',       mm_cpp.gpu_REQs[0, :mm_cpp.natoms, :4], n=n)
        _print_arr('CPP.gpu_MMpars',     mm_cpp.gpu_MMpars[0, :mm_cpp.nnode, :4], n=n)
        _print_arr('CPP.gpu_BLs',        mm_cpp.gpu_BLs[0, :mm_cpp.nnode, :4], n=n)
        _print_arr('CPP.gpu_BKs',        mm_cpp.gpu_BKs[0, :mm_cpp.nnode, :4], n=n)
        _print_arr('CPP.gpu_Ksp',        mm_cpp.gpu_Ksp[0, :mm_cpp.nnode, :4], n=n)
        _print_arr('CPP.gpu_Kpp',        mm_cpp.gpu_Kpp[0, :mm_cpp.nnode, :4], n=n)
        _print_arr('CPP.gpu_atoms',      mm_cpp.gpu_atoms[0, :mm_cpp.nvecs, :4], n=n)
        _print_arr('CPP.gpu_aforces',    mm_cpp.gpu_aforces[0, :mm_cpp.natoms, :4], n=n)

    # ---- PyOpenCL path
    mol = AtomicSystem(fname=args.xyz)
    if getattr(mol, 'ngs', None) is None: mol.neighs()

    mm = MMFF_pyocl(bTorsion=False, verbosity=0)
    if args.force_node_all:
        mm.capping_atoms = set()   # treat all atoms as nodes to match C++/OCL builder layout
        mm.reorder_nodes_first = False

    mm.toMMFFsp3_loc(mol=mol, atom_types=mm.atom_types, bRealloc=True, bEPairs=False, bUFF=False)
    mm.make_back_neighs(b_cap_neighs=not args.force_node_all)

    md = MolecularDynamics(nloc=32)
    md.realloc(mm, nSystems=1)
    md.setup_kernels()
    md.set_pack_system_debug(False)
    md.pack_system(0, mm)

    # force-only: clean -> getMMFFf4 -> updateAtomsMMFFf4 with dt=0 for recoil assembly
    md.kernel_params['bSubtractVdW'] = np.int32(0 if not args.nonbond else 1)

    # set dt=0 in MDparams buffer (dt, damp, damp)
    md.toGPU('MDparams', np.array([0.0, mm.damp, mm.damp, 0.0], dtype=np.float32), byte_offset=0)

    md.run_cleanForceMMFFf4()
    md.run_getMMFFf4()
    md.run_updateAtomsMMFFf4()

    # download forces
    af = np.empty((mm.nvecs, 4), dtype=np.float32)
    import pyopencl as cl
    cl.enqueue_copy(md.queue, af, md.buffer_dict['aforce'])
    md.queue.finish()

    F_py = np.asarray(af[:mm.natoms, :3], dtype=np.float64)
    E_py = float(np.sum(af[:mm.natoms, 3]))

    if args.dump:
        n = args.dump_n
        _print_arr('PY.neighs(int4)',       mm.neighs[:mm.natoms, :], n=n)
        _print_arr('PY.back_neighs',        mm.back_neighs[:mm.nvecs, :], n=n)
        _print_arr('PY.REQs(float4)',       mm.REQs[:mm.natoms, :4], n=n)
        _print_arr('PY.apars(float4)',      mm.apars[:mm.nnode, :4], n=n)
        _print_arr('PY.bLs(float4)',        mm.bLs[:mm.nnode, :4], n=n)
        _print_arr('PY.bKs(float4)',        mm.bKs[:mm.nnode, :4], n=n)
        _print_arr('PY.Ksp(float4)',        mm.Ksp[:mm.nnode, :4], n=n)
        _print_arr('PY.Kpp(float4)',        mm.Kpp[:mm.nnode, :4], n=n)
        _print_arr('PY.apos(float4)',       mm.apos[:mm.nvecs, :4], n=n)
        _print_arr('PY.aforce(float4)',     af[:mm.natoms, :4], n=n)

    # ---- Input parity (only if shapes match)
    ok_buf = True
    print("========= MMFF PyOpenCL vs C++ OpenCL =========")
    print(f"molecule: {args.xyz}")
    label = args.label if args.label else f"nonbond={args.nonbond} angles={args.angles} pisigma={args.pisigma} pipii={args.pipii}"
    print(f"preset: {label}")
    print(f"C++ dims natoms={mm_cpp.natoms} nnode={mm_cpp.nnode} nvecs={mm_cpp.nvecs} | PyOCL natoms={mm.natoms} nnode={mm.nnode} nvecs={mm.nvecs}")

    if (mm_cpp.natoms == mm.natoms) and (mm_cpp.nnode == mm.nnode) and (mm_cpp.nvecs == mm.nvecs):
        print("\n--- Buffer parity (builder outputs) ---")
        ok_buf &= compare('neighs', mm_cpp.neighs[:mm_cpp.natoms, :], mm.neighs[:mm.natoms, :], tol=args.tolBuf)
        ok_buf &= compare('bkNeighs', mm_cpp.gpu_bkNeighs[0, :mm_cpp.nvecs, :], mm.back_neighs[:mm.nvecs, :], tol=0.0)

        # C++ path exposes GPU-packed float4 buffers; compare all 4 components (w can matter)
        ok_buf &= compare('REQs',   mm_cpp.gpu_REQs[0, :mm_cpp.natoms, :4], mm.REQs[:mm.natoms, :4], tol=1e-6)
        ok_buf &= compare('MMpars', mm_cpp.gpu_MMpars[0, :mm_cpp.nnode, :4], mm.apars[:mm.nnode, :4], tol=1e-6)
        ok_buf &= compare('BLs',    mm_cpp.gpu_BLs[0, :mm_cpp.nnode, :4], mm.bLs[:mm.nnode, :4], tol=1e-6)
        ok_buf &= compare('BKs',    mm_cpp.gpu_BKs[0, :mm_cpp.nnode, :4], mm.bKs[:mm.nnode, :4], tol=1e-5)
        ok_buf &= compare('Ksp',    mm_cpp.gpu_Ksp[0, :mm_cpp.nnode, :4], mm.Ksp[:mm.nnode, :4], tol=1e-6)
        ok_buf &= compare('Kpp',    mm_cpp.gpu_Kpp[0, :mm_cpp.nnode, :4], mm.Kpp[:mm.nnode, :4], tol=1e-6)

        # positions: compare full vectors buffer (atoms + pi) including w
        ok_buf &= compare('gpu_atoms', mm_cpp.gpu_atoms[0, :mm_cpp.nvecs, :4], mm.apos[:mm.nvecs, :4], tol=1e-6)
    else:
        print('SKIP input buffer compare: natoms/nnode differ')
        ok_buf = False

    # ---- Force parity
    dF = F_py - F_cpp
    max_abs = _absmax(dF)
    rms = float(np.sqrt(np.mean(dF*dF)))
    idx = np.unravel_index(int(np.argmax(np.abs(dF))), dF.shape)
    cmin, cmax = _minmax(F_cpp)
    pmin, pmax = _minmax(F_py)

    print("\n--- Force parity ---")
    print(f"Max |ΔF| = {max_abs:.6e}")
    print(f"RMS  ΔF  = {rms:.6e}")
    okF = max_abs <= args.tolF
    print(f"{'PASS' if okF else 'FAIL'}: forces within tol {args.tolF:.1e} (worst idx={idx} Δ={dF[idx]})")
    print(f"Force ranges: C++ [{cmin:.3e}, {cmax:.3e}]  PyOCL [{pmin:.3e}, {pmax:.3e}]")

    print('\nPASS' if (ok_buf and okF) else '\nFAIL')

    if args.save_npz:
        np.savez(
            args.save_npz,
            xyz=args.xyz,
            neighs_cpp=mm_cpp.neighs[:mm_cpp.natoms, :],
            bkNeighs_cpp=mm_cpp.gpu_bkNeighs[0, :mm_cpp.nvecs, :],
            REQs_cpp=mm_cpp.gpu_REQs[0, :mm_cpp.natoms, :4],
            MMpars_cpp=mm_cpp.gpu_MMpars[0, :mm_cpp.nnode, :4],
            BLs_cpp=mm_cpp.gpu_BLs[0, :mm_cpp.nnode, :4],
            BKs_cpp=mm_cpp.gpu_BKs[0, :mm_cpp.nnode, :4],
            Ksp_cpp=mm_cpp.gpu_Ksp[0, :mm_cpp.nnode, :4],
            Kpp_cpp=mm_cpp.gpu_Kpp[0, :mm_cpp.nnode, :4],
            atoms_cpp=mm_cpp.gpu_atoms[0, :mm_cpp.nvecs, :4],
            aforce_cpp=mm_cpp.gpu_aforces[0, :mm_cpp.natoms, :4],
            neighs_py=mm.neighs[:mm.natoms, :],
            bkNeighs_py=mm.back_neighs[:mm.nvecs, :],
            REQs_py=mm.REQs[:mm.natoms, :4],
            MMpars_py=mm.apars[:mm.nnode, :4],
            BLs_py=mm.bLs[:mm.nnode, :4],
            BKs_py=mm.bKs[:mm.nnode, :4],
            Ksp_py=mm.Ksp[:mm.nnode, :4],
            Kpp_py=mm.Kpp[:mm.nnode, :4],
            atoms_py=mm.apos[:mm.nvecs, :4],
            aforce_py=af[:mm.natoms, :4],
        )

    if args.fast_exit:
        sys.stdout.flush(); sys.stderr.flush(); os._exit(0)
    return 0 if (ok_buf and okF) else 1


if __name__ == '__main__':
    raise SystemExit(main())
