#!/usr/bin/env python3

import os, sys, argparse
import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(BASE)

from pyBall.MD_test_utils import configure_md, evaluate_mmff_gpu, evaluate_mmff_cpp, mmff_cpp_init_from_atoms, force_parity_report

np.set_printoptions(linewidth=np.inf)

DATA_DIR = os.path.join(BASE, 'cpp/common_resources')
DEFAULT_MOL = os.path.join(DATA_DIR, 'mol', 'formic_acid.mol2')


def _print_arr(name, A, n=0):
    if A is None:
        print(f"{name}: None")
        return
    A = np.asarray(A)
    if (n is not None) and (n > 0) and (A.ndim >= 1) and (A.shape[0] > n):
        print(f"{name} shape={A.shape} (head {n})\n", A[:n])
    else:
        print(f"{name} shape={A.shape}\n", A)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--molecule', type=str, default=DEFAULT_MOL)
    ap.add_argument('--tol', type=float, default=5e-5)
    ap.add_argument('--do-nb', type=int, default=0)
    ap.add_argument('--angles', type=int, default=1)
    ap.add_argument('--pisigma', type=int, default=1)
    ap.add_argument('--pipii', type=int, default=1)
    ap.add_argument('--print-forces', type=int, default=0)
    ap.add_argument('--print-energy', type=int, default=1)
    ap.add_argument('--dump', type=int, default=0, help='If >0, print full input/output buffers (use --dump-n to limit rows)')
    ap.add_argument('--dump-n', type=int, default=0, help='If >0, print only first N rows of long buffers')
    ap.add_argument('--save-npz', type=str, default='', help='If non-empty, save compared buffers/forces to NPZ at this path')
    ap.add_argument('--fast-exit', type=int, default=0, help='Exit via os._exit(0) to bypass known C++ double-free at interpreter shutdown')
    args = ap.parse_args()

    # --- GPU side (PyOpenCL)
    mol, mm, md = configure_md(args.molecule, dt=0.01, damp=1.0, flim=10.0, subtract_vdw=False, drive_temp=0.0, drive_gamma=0.0, drive_seed=0.0, print_params=False)

    if args.dump:
        n = args.dump_n
        _print_arr('PY.apos(float4)',  mm.apos[:mm.nvecs, :4], n=n)
        _print_arr('PY.neighs(int4)',  mm.neighs[:mm.natoms, :], n=n)
        if hasattr(mm, 'back_neighs') and (mm.back_neighs is not None):
            _print_arr('PY.back_neighs', mm.back_neighs[:mm.nvecs, :], n=n)
        if hasattr(mm, 'REQs') and (mm.REQs is not None):
            _print_arr('PY.REQs(float4)',  mm.REQs[:mm.natoms, :4], n=n)
        if hasattr(mm, 'apars') and (mm.apars is not None):
            _print_arr('PY.apars(float4)', mm.apars[:mm.nnode, :4], n=n)
        if hasattr(mm, 'bLs') and (mm.bLs is not None):
            _print_arr('PY.bLs(float4)',   mm.bLs[:mm.nnode, :4], n=n)
        if hasattr(mm, 'bKs') and (mm.bKs is not None):
            _print_arr('PY.bKs(float4)',   mm.bKs[:mm.nnode, :4], n=n)
        if hasattr(mm, 'Ksp') and (mm.Ksp is not None):
            _print_arr('PY.Ksp(float4)',   mm.Ksp[:mm.nnode, :4], n=n)
        if hasattr(mm, 'Kpp') and (mm.Kpp is not None):
            _print_arr('PY.Kpp(float4)',   mm.Kpp[:mm.nnode, :4], n=n)

    # --- CPU side (C++ reference via pyBall.MMFF)
    mmff_cpp_init_from_atoms(mm.apos[:mm.natoms, :3], mol.enames[:mm.natoms], DATA_DIR, nPBC=(0,0,0), b141=True, bSimple=False, bConj=True, bCumulene=True)

    # evaluate at same geometry
    E_cpu, F_cpu = evaluate_mmff_cpp(mm.apos[:mm.natoms, :3], do_nonbond=bool(args.do_nb), do_angles=bool(args.angles), do_pisigma=bool(args.pisigma), do_pipii=bool(args.pipii), do_pbc=False)
    E_gpu, F_gpu = evaluate_mmff_gpu(mm, md, mm.apos, do_clean=True, do_nb=bool(args.do_nb), do_mmff=True, mode='none')

    rep = force_parity_report(F_cpu, F_gpu, tol=args.tol, label_ref='CPU(C++)', label_tst='GPU(PyOCL)')

    if args.print_energy:
        print(f"E_cpu={E_cpu:+.8e}  E_gpu={E_gpu:+.8e}  dE={E_gpu-E_cpu:+.3e}")

    print(f"Max |dF| = {rep['max_abs']:.6e}  RMS dF = {rep['rms']:.6e}  tol={rep['tol']:.2e}")
    ia = rep['idx'][0]
    print(f"worst atom {ia}: CPU={rep['ref']}  GPU={rep['tst']}  dF={rep['diff']}")

    if args.print_forces:
        print('CPU forces:\n', F_cpu)
        print('GPU forces:\n', F_gpu)

    if args.dump:
        n = args.dump_n
        _print_arr('F_cpu(C++)', F_cpu, n=n)
        _print_arr('F_gpu(PyOCL)', F_gpu, n=n)
        _print_arr('dF', (F_gpu - F_cpu), n=n)

    if args.save_npz:
        np.savez(
            args.save_npz,
            molecule=args.molecule,
            apos=mm.apos[:mm.nvecs, :4],
            neighs=mm.neighs[:mm.natoms, :],
            REQs=getattr(mm, 'REQs', None),
            apars=getattr(mm, 'apars', None),
            bLs=getattr(mm, 'bLs', None),
            bKs=getattr(mm, 'bKs', None),
            Ksp=getattr(mm, 'Ksp', None),
            Kpp=getattr(mm, 'Kpp', None),
            E_cpu=E_cpu,
            E_gpu=E_gpu,
            F_cpu=F_cpu,
            F_gpu=F_gpu,
        )

    if rep['ok']:
        print('PASS')
        if args.fast_exit:
            sys.stdout.flush(); sys.stderr.flush(); os._exit(0)
        return 0
    else:
        print('FAIL')
        if args.fast_exit:
            sys.stdout.flush(); sys.stderr.flush(); os._exit(0)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
