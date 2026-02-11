#!/usr/bin/env python3

import os, sys, argparse
import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(BASE)

from pyBall import MMFF_multi as mm

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


def compare(name, A, B, tol=0.0):
    A = np.asarray(A)
    B = np.asarray(B)
    if A.shape != B.shape:
        print(f"BUF MISMATCH {name}: shape {A.shape} != {B.shape}")
        return False
    d = B - A
    m = _absmax(d)
    ok = m <= tol
    print(f"BUF {name}: max|d|={m:.6e} tol={tol:.1e} {'OK' if ok else 'FAIL'}")
    if not ok:
        idx = np.unravel_index(int(np.argmax(np.abs(d))), d.shape)
        print(f"  worst idx={idx} A={A[idx]} B={B[idx]} d={d[idx]}")
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--xyz', type=str, default=os.path.join(XYZ_DIR, 'HCOOH.xyz'))
    ap.add_argument('--nsys', type=int, default=1)
    ap.add_argument('--tolF', type=float, default=1e-4)
    ap.add_argument('--tolBufI', type=float, default=0.0)
    ap.add_argument('--tolBufF', type=float, default=0.0)
    ap.add_argument('--nonbond', type=int, default=0)
    ap.add_argument('--angles', type=int, default=1)
    ap.add_argument('--pisigma', type=int, default=1)
    ap.add_argument('--pipii', type=int, default=1)
    ap.add_argument('--print-debugs', type=int, default=1)
    ap.add_argument('--dump', type=int, default=0, help='If >0, print full input/output buffers (use --dump-n to limit rows)')
    ap.add_argument('--dump-n', type=int, default=0, help='If >0, print only first N rows of long buffers')
    ap.add_argument('--save-npz', type=str, default='', help='If non-empty, save compared buffers/forces to NPZ at this path')
    ap.add_argument('--fast-exit', type=int, default=0, help='Exit via os._exit(0) to bypass known C++ double-free at interpreter shutdown')
    args = ap.parse_args()

    # IMPORTANT: run from tests/tUFF so data/ symlinks resolve
    # init uses relative default paths (data/*.dat)
    mm.init(
        nSys_=args.nsys,
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

    # Force switches (0=keep, >0 on, <0 off)
    mm.setSwitches2(
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

    # expose buffers
    mm.getBuffs()

    if args.print_debugs:
        mm.print_debugs(True, True, False, False)

    # --- INPUT PARITY: compare CPU packed vs GPU packed buffers (system 0)
    ok_buf = True
    # neighs: CPU uses nnode x 4; GPU uses nvecs x 4 (atoms+pi). Compare only node part.
    ok_buf &= compare('neighs(node)', mm.neighs[:mm.nnode,:], mm.gpu_neighs[0,:mm.nnode,:], tol=args.tolBufI)

    if args.dump:
        n = args.dump_n
        _print_arr('CPU.neighs(node)', mm.neighs[:mm.nnode, :], n=n)
        _print_arr('GPU.neighs(node)', mm.gpu_neighs[0, :mm.nnode, :], n=n)
        for nm in ('gpu_REQs','gpu_MMpars','gpu_BLs','gpu_BKs','gpu_Ksp','gpu_Kpp','gpu_atoms'):
            if hasattr(mm, nm):
                _print_arr(f'GPU.{nm}', getattr(mm, nm)[0], n=n)

    # compare packed param arrays on GPU (float32) vs CPU-derived buffers are not directly exposed here;
    # best-effort check: make sure GPU arrays are finite and reasonable
    for nm in ('gpu_REQs','gpu_MMpars','gpu_BLs','gpu_BKs','gpu_Ksp','gpu_Kpp'):
        A = getattr(mm, nm)
        if not np.isfinite(A).all():
            print(f"BUF {nm}: contains NaN/Inf")
            ok_buf = False
        else:
            print(f"BUF {nm}: finite OK shape={A.shape} min={A.min():.3e} max={A.max():.3e}")

    # --- FORCE PARITY
    # Both CPU and GPU paths copy forces to nbmol.fapos; read from mm.fapos[:natoms,:3]
    na = mm.natoms

    E_cpu = mm.eval_force_MMFF(0)
    F_cpu = np.array(mm.fapos[:na,:3], copy=True, dtype=np.float64)

    E_gpu = mm.eval_force_MMFF(2)
    F_gpu = np.array(mm.fapos[:na,:3], copy=True, dtype=np.float64)

    if args.dump:
        n = args.dump_n
        _print_arr('F_cpu', F_cpu, n=n)
        _print_arr('F_gpu', F_gpu, n=n)
        _print_arr('dF', (F_gpu - F_cpu), n=n)

    dF = F_gpu - F_cpu
    max_abs = _absmax(dF)
    rms = float(np.sqrt(np.mean(dF*dF)))
    idx = np.unravel_index(int(np.argmax(np.abs(dF))), dF.shape)

    print(f"\nE_cpu={E_cpu:+.8e} E_gpu={E_gpu:+.8e} dE={E_gpu-E_cpu:+.3e}")
    print(f"Max|dF|={max_abs:.6e} RMS={rms:.6e} tolF={args.tolF:.1e} worst={idx} dF={dF[idx]}")
    # per-atom forces
    for ia in range(na):
        d = F_gpu[ia] - F_cpu[ia]
        print(f"  atom[{ia}] CPU=({F_cpu[ia,0]:+.6f},{F_cpu[ia,1]:+.6f},{F_cpu[ia,2]:+.6f}) GPU=({F_gpu[ia,0]:+.6f},{F_gpu[ia,1]:+.6f},{F_gpu[ia,2]:+.6f}) dF=({d[0]:+.3e},{d[1]:+.3e},{d[2]:+.3e})")

    okF = max_abs <= args.tolF
    print('\nPASS' if (ok_buf and okF) else '\nFAIL')

    if args.save_npz:
        np.savez(
            args.save_npz,
            xyz=args.xyz,
            neighs_cpu=mm.neighs[:mm.nnode, :],
            neighs_gpu=mm.gpu_neighs[0, :mm.nnode, :],
            gpu_REQs=mm.gpu_REQs[0],
            gpu_MMpars=mm.gpu_MMpars[0],
            gpu_BLs=mm.gpu_BLs[0],
            gpu_BKs=mm.gpu_BKs[0],
            gpu_Ksp=mm.gpu_Ksp[0],
            gpu_Kpp=mm.gpu_Kpp[0],
            gpu_atoms=mm.gpu_atoms[0],
            E_cpu=E_cpu,
            E_gpu=E_gpu,
            F_cpu=F_cpu,
            F_gpu=F_gpu,
        )

    if args.fast_exit:
        sys.stdout.flush(); sys.stderr.flush(); os._exit(0)
    return 0 if (ok_buf and okF) else 1


if __name__ == '__main__':
    raise SystemExit(main())
