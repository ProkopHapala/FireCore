#!/usr/bin/env python3

import os
import sys
import argparse
import numpy as np

base_path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(base_path)

data_dir = os.path.join(base_path, "cpp/common_resources")

from pyBall import MMFF_multi as uff_cppocl
from pyBall.AtomicSystem import AtomicSystem
from pyBall.OCL import UFF as uff_pyocl


def _parse_components(s: str):
    comp_set = {c.strip() for c in (s or "").split(',') if c.strip()}
    need_bonds = ('bonds' in comp_set) or ('angles' in comp_set) or ('dihedrals' in comp_set) or ('inversions' in comp_set)
    return {
        'bonds': need_bonds,
        'angles': 'angles' in comp_set,
        'dihedrals': 'dihedrals' in comp_set,
        'inversions': 'inversions' in comp_set,
    }


def _pad_int_cols(a, ncol, fill=-1):
    a = np.asarray(a)
    assert a.ndim == 2
    if a.shape[1] == ncol:
        return a
    assert a.shape[1] < ncol
    out = np.full((a.shape[0], ncol), fill, dtype=a.dtype)
    out[:, :a.shape[1]] = a
    return out


def _pyocl_download_buf(uff_cl, name, shape, dtype):
    import pyopencl as cl
    buf = uff_cl.buffer_dict.get(name, None)
    if buf is None:
        return np.zeros(shape, dtype=dtype)
    host = np.zeros(int(np.prod(shape)), dtype=dtype)
    cl.enqueue_copy(uff_cl.queue, host, buf)
    return host.reshape(shape)


def _set_cpp_switches(components):
    # Disable all noncovalent paths (semantics: 0=keep, <0=force off, >0=force on)
    uff_cppocl.setSwitches2(PBC=-1, NonBonded=-1, NonBondNeighs=-1, SurfAtoms=-1, GridFF=-1)

    DoBond      = 1 if components['bonds'] else -1
    DoAngle     = 1 if components['angles'] else -1
    DoDihedral  = 1 if components['dihedrals'] else -1
    DoInversion = 1 if components['inversions'] else -1
    DoAssemble  = 1 if (DoAngle > 0 or DoDihedral > 0 or DoInversion > 0) else -1

    uff_cppocl.setSwitchesUFF(
        DoBond=DoBond,
        DoAngle=DoAngle,
        DoDihedral=DoDihedral,
        DoInversion=DoInversion,
        DoAssemble=DoAssemble,
        SubtractBondNonBond=-1,
        ClampNonBonded=-1,
    )
    return DoBond, DoAngle, DoDihedral, DoInversion, DoAssemble


def _run_cpp_opencl_forces(mol_path: str, components: dict):
    os.chdir(os.path.join(base_path, "cpp"))
    uff_cppocl.init(
        nSys_=1,
        xyz_name=mol_path,
        surf_name=None,
        sElementTypes=os.path.join(data_dir, "ElementTypes.dat"),
        sAtomTypes=os.path.join(data_dir, "AtomTypes.dat"),
        sBondTypes=os.path.join(data_dir, "BondTypes.dat"),
        sAngleTypes=os.path.join(data_dir, "AngleTypes.dat"),
        sDihedralTypes=os.path.join(data_dir, "DihedralTypes.dat"),
        bMMFF=False,
        bUFF=True,
        b141=True,
        bSimple=True,
        bConj=True,
        bCumulene=True,
        nPBC=(0,0,0),
    )

    uff_cppocl.getBuffs_UFF()
    sw = _set_cpp_switches(components)

    # Use scan() to evaluate forces at fixed geometry (no integration)
    confs = uff_cppocl.apos.copy()[None, :, :]
    F = uff_cppocl.scan(confs, iParalel=2)[0]

    # Inputs/buffers from C++ builder (host-side)
    bufs = {
        'bonAtoms':  uff_cppocl.bonAtoms.copy(),
        'angAtoms':  _pad_int_cols(uff_cppocl.angAtoms.copy(), 4, fill=-1),
        'dihAtoms':  uff_cppocl.dihAtoms.copy(),
        'invAtoms':  uff_cppocl.invAtoms.copy(),
        'neighs':    uff_cppocl.neighs.copy(),
        'neighBs':   uff_cppocl.neighBs.copy(),
        'bonParams': uff_cppocl.bonParams.copy(),
        'angParams': uff_cppocl.angParams.copy(),
        'dihParams': uff_cppocl.dihParams.copy(),
        'invParams': uff_cppocl.invParams.copy(),
        'angNgs':    uff_cppocl.angNgs.copy(),
        'dihNgs':    _pad_int_cols(uff_cppocl.dihNgs.copy(), 4, fill=-1),
        'invNgs':    uff_cppocl.invNgs.copy(),
        'hneigh':    (uff_cppocl.gpu_hneigh.copy() if (getattr(uff_cppocl,'gpu_hneigh',None) is not None) else uff_cppocl.hneigh.copy()),
        'a2f_offsets': uff_cppocl.a2f_offsets.copy(),
        'a2f_counts':  uff_cppocl.a2f_counts.copy(),
        'a2f_indices': uff_cppocl.a2f_indices.copy(),
    }
    return F, bufs, sw


def _run_pyopencl_forces(mol_path: str, components: dict, debug_opts=None):
    mol = AtomicSystem(fname=mol_path)

    uff_cl = uff_pyocl.UFF_CL(nloc=32, debug_build_options=debug_opts)
    uff_data = uff_cl.toUFF(mol)

    uff_cl.upload_topology_params(uff_data)
    uff_cl.upload_positions(mol.apos)

    uff_cl.bDoBonds      = components['bonds']
    uff_cl.bDoAngles     = components['angles']
    uff_cl.bDoDihedrals  = components['dihedrals']
    uff_cl.bDoInversions = components['inversions']

    uff_cl.run_eval_step(bClearForce=True)
    F = uff_cl.get_forces(iSys=0)

    # Normalize some layout differences for comparison
    bufs = {
        'bonAtoms':  np.ascontiguousarray(uff_data.get('bonAtoms', np.zeros((0, 2), np.int32)), dtype=np.int32),
        'angAtoms':  np.ascontiguousarray(uff_data.get('angAtoms', np.zeros((0, 4), np.int32)), dtype=np.int32),
        'dihAtoms':  np.ascontiguousarray(uff_data.get('dihAtoms', np.zeros((0, 4), np.int32)), dtype=np.int32),
        'invAtoms':  np.ascontiguousarray(uff_data.get('invAtoms', np.zeros((0, 4), np.int32)), dtype=np.int32),
        'neighs':    np.ascontiguousarray(uff_data.get('neighs', np.zeros((0, 4), np.int32)), dtype=np.int32),
        'neighBs':   np.ascontiguousarray(uff_data.get('neighBs', np.zeros((0, 4), np.int32)), dtype=np.int32),
        'bonParams': np.ascontiguousarray(uff_data.get('bonParams', np.zeros((0, 2), np.float32)), dtype=np.float32),
        'angParams': np.ascontiguousarray(uff_data.get('angParams', np.zeros((0, 5), np.float32)), dtype=np.float32),
        'dihParams': np.ascontiguousarray(uff_data.get('dihParams', np.zeros((0, 3), np.float32)), dtype=np.float32),
        'invParams': np.ascontiguousarray(uff_data.get('invParams', np.zeros((0, 4), np.float32)), dtype=np.float32),
        'angNgs':    np.ascontiguousarray(uff_data.get('angNgs', np.zeros((0, 2), np.int32)), dtype=np.int32),
        'dihNgs':    np.ascontiguousarray(uff_data.get('dihNgs', np.zeros((0, 4), np.int32)), dtype=np.int32),
        'invNgs':    np.ascontiguousarray(uff_data.get('invNgs', np.zeros((0, 3), np.int32)), dtype=np.int32),
        # torsion-driving buffers
        'hneigh':    np.ascontiguousarray(_pyocl_download_buf(uff_cl, 'hneigh', shape=(uff_cl.natoms*4, 4), dtype=np.float32), dtype=np.float32),
        'a2f_offsets': np.ascontiguousarray(uff_data.get('a2f_offsets', np.zeros((0,), np.int32)), dtype=np.int32),
        'a2f_counts':  np.ascontiguousarray(uff_data.get('a2f_counts',  np.zeros((0,), np.int32)), dtype=np.int32),
        'a2f_indices': np.ascontiguousarray(uff_data.get('a2f_indices', np.zeros((0,), np.int32)), dtype=np.int32),
    }

    return F, bufs


def _max_abs_diff(a, b):
    d = np.asarray(a) - np.asarray(b)
    if d.size == 0:
        return 0.0
    return float(np.max(np.abs(d)))


def _min_max(arr):
    arr = np.asarray(arr)
    if arr.size == 0:
        return (0.0, 0.0, True)
    return float(np.min(arr)), float(np.max(arr)), False


def _canon_bonAtoms(a):
    a = np.asarray(a)
    if a.size == 0:
        return a
    assert a.ndim == 2 and a.shape[1] == 2
    aa = np.array(a, copy=True)
    i = np.minimum(aa[:, 0], aa[:, 1])
    j = np.maximum(aa[:, 0], aa[:, 1])
    aa[:, 0] = i
    aa[:, 1] = j
    key = aa[:, 0].astype(np.int64) * (np.int64(1) << 32) + aa[:, 1].astype(np.int64)
    return aa[np.argsort(key)]


def _compare_bufs(cpp_bufs, py_bufs, tolF=0.0):
    ok = True
    diffs = {}
    for k, a in cpp_bufs.items():
        b = py_bufs.get(k, None)
        if b is None:
            diffs[k] = (np.nan, 'missing_py')
            ok = False
            continue
        # Some builders represent empty buffers differently (e.g. (0,4) vs (0,)).
        # If both are empty, consider it a match regardless of shape.
        if (np.asarray(a).size == 0) and (np.asarray(b).size == 0):
            diffs[k] = (0.0, '')
            continue
        if a.shape != b.shape:
            diffs[k] = (np.nan, f"shape {a.shape} vs {b.shape}")
            ok = False
            continue
        tol = 0.0
        if a.dtype.kind == 'f' or b.dtype.kind == 'f':
            tol = tolF
        md = _max_abs_diff(a, b)
        diffs[k] = (md, '')
        if md > tol:
            ok = False
    return ok, diffs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--molecule', type=str, required=True)
    ap.add_argument('--components', type=str, default='bonds,angles,dihedrals,inversions')
    ap.add_argument('--tolerance', type=float, default=5e-5)
    ap.add_argument('--tol-buf-f', type=float, default=1e-6)
    ap.add_argument('--dump', type=int, default=0)
    ap.add_argument('--dump-n', type=int, default=0)
    ap.add_argument('--fast-exit', type=int, default=1)
    ap.add_argument('--ocl-dbg', type=int, default=0)
    ap.add_argument('--ocl-dbg-sys', type=int, default=0)
    ap.add_argument('--ocl-dbg-dih', type=int, default=0)
    ap.add_argument('--check-torsion-inputs', type=int, default=1)
    args = ap.parse_args()

    mol_path = os.path.abspath(args.molecule)
    os.chdir(base_path)

    comps = _parse_components(args.components)

    debug_opts = None
    if args.ocl_dbg > 0:
        debug_opts = [
            f"-DDBG_UFF={args.ocl_dbg}",
            f"-DIDBG_SYS={args.ocl_dbg_sys}",
            f"-DIDBG_DIH={args.ocl_dbg_dih}",
        ]

    run_err = None
    F_cpp = bufs_cpp = sw = None
    F_py  = bufs_py  = None
    try:
        F_cpp, bufs_cpp, sw = _run_cpp_opencl_forces(mol_path, comps)
        F_py,  bufs_py      = _run_pyopencl_forces(mol_path, comps, debug_opts=debug_opts)
    except Exception as e:
        run_err = e

    if run_err is not None:
        print("\n========= UFF PyOpenCL vs C++ OpenCL =========")
        print(f"molecule: {mol_path}")
        print(f"components: {args.components}")
        print("FAIL: exception during run")
        import traceback
        traceback.print_exception(type(run_err), run_err, run_err.__traceback__)
        if args.fast_exit:
            sys.stdout.flush(); sys.stderr.flush(); os._exit(1)
        raise

    # Canonicalize bond ordering to avoid false FAILs due to list order
    if 'bonAtoms' in bufs_cpp: bufs_cpp['bonAtoms'] = _canon_bonAtoms(bufs_cpp['bonAtoms'])
    if 'bonAtoms' in bufs_py:  bufs_py ['bonAtoms'] = _canon_bonAtoms(bufs_py ['bonAtoms'])

    if args.check_torsion_inputs <= 0:
        for k in ('hneigh','a2f_offsets','a2f_counts','a2f_indices'):
            bufs_cpp.pop(k, None)
            bufs_py.pop(k, None)

    # For a2f_indices compare only meaningful prefix (sum of counts). Tail is unused capacity in C++ packing.
    if ('a2f_counts' in bufs_cpp) and ('a2f_indices' in bufs_cpp) and ('a2f_counts' in bufs_py) and ('a2f_indices' in bufs_py):
        nref_cpp = int(np.sum(bufs_cpp['a2f_counts']))
        nref_py  = int(np.sum(bufs_py ['a2f_counts']))
        nref = min(nref_cpp, nref_py)
        bufs_cpp['a2f_indices'] = np.ascontiguousarray(bufs_cpp['a2f_indices'][:nref])
        bufs_py ['a2f_indices'] = np.ascontiguousarray(bufs_py ['a2f_indices'][:nref])

    # Compare forces (fail loudly on NaN/Inf)
    F_cpp64 = np.asarray(F_cpp, dtype=np.float64)
    F_py64  = np.asarray(F_py,  dtype=np.float64)
    minF_cpp = float(np.min(F_cpp64)) if F_cpp64.size else 0.0
    maxF_cpp = float(np.max(F_cpp64)) if F_cpp64.size else 0.0
    minF_py  = float(np.min(F_py64 )) if F_py64.size  else 0.0
    maxF_py  = float(np.max(F_py64 )) if F_py64.size  else 0.0
    bad_cpp = ~np.isfinite(F_cpp64)
    bad_py  = ~np.isfinite(F_py64)
    has_bad = bool(bad_cpp.any() or bad_py.any())
    if has_bad:
        n_bad_cpp = int(bad_cpp.sum())
        n_bad_py  = int(bad_py.sum())
        dF = np.full_like(F_py64, np.nan)
        max_dF = float('inf')
        rms_dF = float('inf')
    else:
        dF = F_py64 - F_cpp64
        max_dF = float(np.max(np.abs(dF))) if dF.size else 0.0
        rms_dF = float(np.sqrt(np.mean(dF*dF))) if dF.size else 0.0

    # Compare buffers (host builder outputs)
    ok_bufs, diffs = _compare_bufs(bufs_cpp, bufs_py, tolF=args.tol_buf_f)

    zero_force = ((maxF_cpp == 0.0 and minF_cpp == 0.0) or (maxF_py == 0.0 and minF_py == 0.0))
    ok_force = (max_dF <= args.tolerance) and (not has_bad) and (not zero_force)
    ok_all = bool(ok_force) and bool(ok_bufs)

    print("\n========= UFF PyOpenCL vs C++ OpenCL =========")
    print(f"molecule: {mol_path}")
    print(f"components: {args.components}")
    print(f"C++ switches (DoBond,DoAngle,DoDihedral,DoInversion,DoAssemble) = {sw}")
    print("\n--- Buffer parity (builder outputs) ---")
    for k, (md, note) in diffs.items():
        if np.isnan(md):
            print(f"{k:10s} : FAIL {note}")
        else:
            tol = args.tol_buf_f if (bufs_cpp[k].dtype.kind == 'f' or bufs_py[k].dtype.kind == 'f') else 0.0
            status = "PASS" if md <= tol else "FAIL"
            min_cpp, max_cpp, empty_cpp = _min_max(bufs_cpp[k])
            min_py,  max_py,  empty_py  = _min_max(bufs_py[k])
            extra = " (empty)" if empty_cpp and empty_py else ""
            print(f"{k:10s} : {status}  max|Δ|={md:.3e}  tol={tol:.1e}  C++[{min_cpp:.3e},{max_cpp:.3e}]  PyOCL[{min_py:.3e},{max_py:.3e}]{extra}")

    print("\n--- Force parity ---")
    if has_bad:
        print(f"FAIL: non-finite forces detected  nBad(C++)={n_bad_cpp}  nBad(PyOCL)={n_bad_py}")
    elif zero_force:
        print(f"FAIL: all forces are zero  C++[min,max]=({minF_cpp:.3e},{maxF_cpp:.3e})  PyOCL[min,max]=({minF_py:.3e},{maxF_py:.3e})")
    else:
        print(f"Max |ΔF| = {max_dF:.6e}")
        print(f"RMS  ΔF  = {rms_dF:.6e}")
        if not ok_force:
            print(f"FAIL: max|ΔF| {max_dF:.6e} > tol {args.tolerance:.2e}")
        else:
            print(f"PASS: forces within tol {args.tolerance:.2e}")

    print(f"Force ranges: C++ [{minF_cpp:.3e}, {maxF_cpp:.3e}]  PyOCL [{minF_py:.3e}, {maxF_py:.3e}]")

    if args.dump or (not ok_all):
        n = args.dump_n
        sl = slice(None) if n <= 0 else slice(0, n)
        print("\n--- FORCES (C++ OCL) ---")
        print(F_cpp[sl])
        print("\n--- FORCES (PyOCL) ---")
        print(F_py[sl])
        print("\n--- AbsDiff(F) ---")
        print(np.abs(dF)[sl])

        print("\n--- BUFFERS (C++ builder) ---")
        for k, v in bufs_cpp.items():
            print(k)
            print(v[sl] if v.ndim >= 2 else v)

        print("\n--- BUFFERS (PyOCL builder) ---")
        for k, v in bufs_py.items():
            print(k)
            print(v[sl] if v.ndim >= 2 else v)

    if args.fast_exit:
        sys.stdout.flush(); sys.stderr.flush(); os._exit(0 if ok_all else 1)

    raise SystemExit(0 if ok_all else 1)


if __name__ == '__main__':
    main()
