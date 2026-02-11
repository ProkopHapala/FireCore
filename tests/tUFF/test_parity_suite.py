#!/usr/bin/env python3

import os, sys, argparse, subprocess, shlex
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent
BASE     = THIS_DIR.parent.parent
DATA_DIR = THIS_DIR / 'data'
XYZ_DIR  = DATA_DIR / 'xyz'
MOL_DIR  = DATA_DIR / 'mol'


def _resolve_mol(name: str) -> str:
    p = Path(name)
    if p.is_absolute() or str(p).startswith('./') or '/' in str(p):
        return str(p)
    ext = p.suffix.lower()
    if ext in ('.mol', '.mol2', '.sdf'):
        return str(MOL_DIR / name)
    return str(XYZ_DIR / name)


def _run(cmd, *, cwd, out_path=None):
    cmd_str = cmd if isinstance(cmd, str) else ' '.join(shlex.quote(c) for c in cmd)
    print(f"\n=== RUN: {cmd_str} ===")
    if out_path is not None:
        out_path = Path(out_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with out_path.open('w') as f:
            p = subprocess.run(cmd_str, cwd=str(cwd), shell=True, stdout=f, stderr=subprocess.STDOUT)
        print(f"log: {out_path}  rc={p.returncode}")
        return p.returncode
    p = subprocess.run(cmd_str, cwd=str(cwd), shell=True)
    return p.returncode


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mols', type=str, default='methanol.mol2,HCONH2.xyz,uracil.xyz,xylitol.mol2,guanine.xyz,Si10_H.xyz')
    ap.add_argument('--outdir', type=str, default='OUT_parity_suite')
    ap.add_argument('--nonbond', type=int, default=0)
    ap.add_argument('--angles', type=int, default=1)
    ap.add_argument('--pisigma', type=int, default=1)
    ap.add_argument('--pipii', type=int, default=1)
    ap.add_argument('--tolF', type=float, default=1e-4)
    ap.add_argument('--tolBuf', type=float, default=0.0)
    ap.add_argument('--rerun-fails', type=int, default=1)
    ap.add_argument('--dump-n', type=int, default=0)
    ap.add_argument('--uff-components', type=str, default='bonds,angles,dihedrals,inversions', help='UFF components to test (comma-separated). Suite will stage up to this set to localize failures.')
    args = ap.parse_args()

    outdir = THIS_DIR / args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    mols = [m.strip() for m in args.mols.split(',') if m.strip()]

    results = []

    # UFF staged sets to localize failures
    want = [c.strip() for c in args.uff_components.split(',') if c.strip()]
    stages = []
    if 'bonds' in want: stages.append('bonds')
    if 'angles' in want: stages.append('bonds,angles')
    if 'dihedrals' in want: stages.append('bonds,angles,dihedrals')
    if 'inversions' in want: stages.append('bonds,angles,dihedrals,inversions')
    if not stages:
        stages = ['bonds,angles,dihedrals,inversions']

    for m in mols:
        mp = _resolve_mol(m)
        tag = Path(m).stem

        # --- MMFF: C++ CPU vs C++ OpenCL
        mmff_cpu_gpu = [
            'python3 -u test_MMFF_multi_parity.py',
            f"--xyz {shlex.quote(mp)}",
            f"--tolF {args.tolF}",
            f"--tolBufI {args.tolBuf}",
            f"--nonbond {args.nonbond}",
            f"--angles {args.angles}",
            f"--pisigma {args.pisigma}",
            f"--pipii {args.pipii}",
            '--print-debugs 0',
            '--fast-exit 1',
        ]
        cmd1 = ' '.join(mmff_cpu_gpu)
        log1 = outdir / f"{tag}__mmff_cpu_vs_ocl.log"
        rc1 = _run(cmd1, cwd=THIS_DIR, out_path=log1)

        # --- MMFF: C++ OpenCL vs PyOpenCL
        mmff_ocl_pyocl = [
            'python3 -u test_MMFF_cpp_vs_pyocl_parity.py',
            f"--xyz {shlex.quote(mp)}",
            f"--tolF {args.tolF}",
            f"--tolBuf {args.tolBuf}",
            f"--nonbond {args.nonbond}",
            f"--angles {args.angles}",
            f"--pisigma {args.pisigma}",
            f"--pipii {args.pipii}",
            '--fast-exit 1',
        ]
        cmd2 = ' '.join(mmff_ocl_pyocl)
        log2 = outdir / f"{tag}__mmff_ocl_vs_pyocl.log"
        rc2 = _run(cmd2, cwd=THIS_DIR, out_path=log2)

        # --- UFF: C++ CPU vs OpenCL
        rc3 = 0
        for st in stages:
            uff_cpu_gpu = [
                'python3 -u test_UFF_ocl.py',
                f"--molecule {shlex.quote(mp)}",
                f"--tolerance {args.tolF}",
                '--gpu 1',
                f"--components {shlex.quote(st)}",
                '--fast-exit 1',
            ]
            cmd3 = ' '.join(uff_cpu_gpu)
            log3 = outdir / f"{tag}__uff_cpu_vs_ocl__{st.replace(',','_')}.log"
            rc3 = _run(cmd3, cwd=THIS_DIR, out_path=log3)
            if args.rerun_fails and (rc3 != 0):
                _run(cmd3 + f" --dump 1 --dump-n {args.dump_n}", cwd=THIS_DIR, out_path=outdir / f"{tag}__uff_cpu_vs_ocl__{st.replace(',','_')}.DUMP.log")
                break

        results.append((m, rc1, rc2, rc3))

        if args.rerun_fails:
            if rc1 != 0:
                _run(cmd1 + f" --dump 1 --dump-n {args.dump_n} --save-npz {shlex.quote(str(outdir / (tag+'__mmff_cpu_vs_ocl.npz')))}", cwd=THIS_DIR, out_path=outdir / f"{tag}__mmff_cpu_vs_ocl.DUMP.log")
            if rc2 != 0:
                _run(cmd2 + f" --dump 1 --dump-n {args.dump_n} --save-npz {shlex.quote(str(outdir / (tag+'__mmff_ocl_vs_pyocl.npz')))}", cwd=THIS_DIR, out_path=outdir / f"{tag}__mmff_ocl_vs_pyocl.DUMP.log")

    print("\n==================== SUMMARY ====================")
    print("molecule | mmff cpu-vs-ocl | mmff ocl-vs-pyocl | uff cpu-vs-ocl")
    ok_all = True
    for m, rc1, rc2, rc3 in results:
        print(f"{m} | {rc1} | {rc2} | {rc3}")
        ok_all &= (rc1 == 0 and rc2 == 0 and rc3 == 0)

    print(f"\nlogs in: {outdir}")
    return 0 if ok_all else 1


if __name__ == '__main__':
    raise SystemExit(main())
