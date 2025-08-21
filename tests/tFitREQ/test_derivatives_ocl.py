#!/usr/bin/env python3

import sys, os
import argparse
import numpy as np
import matplotlib.pyplot as plt

# keep imports minimal and uniform (match tests/tFitREQ/test_derivatives.py style)
# IMPORTANT: do NOT clear sys.path or we lose site-packages (e.g., pyopencl)
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../"))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)
from pyBall.OCL.NonBondFitting import FittingDriver, extract_macro_block

np.set_printoptions(linewidth=200)

# --- utils --------------------------------------------------------------

def parse_dof_minmax(fname):
    mins, maxs, names = [], [], []
    comps = "REQH"
    with open(fname,'r') as f:
        for line in f:
            if not line.strip() or line[0]=='#':
                continue
            w = line.split()
            t, ic = w[0], int(w[1])
            mn, mx = float(w[2]), float(w[3])
            names.append(f"{t}.{comps[ic]}")
            mins.append(mn); maxs.append(mx)
    return np.array(mins), np.array(maxs), names

def num_deriv(xs, Es):
    # central difference on interior points
    xs = np.asarray(xs)
    Es = np.asarray(Es)
    dE = (Es[2:] - Es[:-2]) / (xs[2:] - xs[:-2])
    xc = 0.5*(xs[2:] + xs[:-2])
    return dE, xc

# --- main ---------------------------------------------------------------

if __name__ == '__main__':
    # run like this:
    #   python3 test_derivatives_ocl.py --xyz HHalogens/porcessed/HF-A1_HF-D1.xyz --dof dofSelection_MorseSR.dat   --points 50 --tol 1e-6 
    #   python3 test_derivatives_ocl.py --xyz all.xyz --dof dofSelection_MorseSR.dat --model MODEL_MorseQ_PAIR --energy_model ENERGY_MorseQ_PAIR --points 100 --eps 1e-6 --tol 1e-6 --save out_  --show --regularize 0 

    
    p = argparse.ArgumentParser(description='OCL: Scan DOFs and compare analytical vs numerical derivatives')
    p.add_argument('--xyz', default='all.xyz')
    p.add_argument('--dof', default='dofSelection_MorseSR.dat')
    p.add_argument('--model', default='MODEL_MorseQ_PAIR')
    p.add_argument('--energy_model', default='ENERGY_MorseQ_PAIR')
    p.add_argument('--points', type=int, default=100)
    p.add_argument('--eps', type=float, default=1e-6)
    p.add_argument('--tol', type=float, default=1e-6)
    p.add_argument('--show', action='store_true')
    p.add_argument('--save', default='')
    p.add_argument('--exclude', nargs='*', default=[])
    p.add_argument('--regularize', type=int, default=1, help='1 to include regularization; 0 to disable (zero stiffness in regParams)')
    p.add_argument('--verbose', type=int, default=0)
    args = p.parse_args()

    # Resolve project root to locate common resources
    this_dir = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.abspath(os.path.join(this_dir, '..', '..'))
    atom_types_file = os.path.join(repo_root, 'cpp', 'common_resources', 'AtomTypes.dat')
    forces_path = os.path.join(repo_root, 'cpp', 'common_resources', 'cl', 'Forces.cl')

    # Setup driver
    drv = FittingDriver(verbose=args.verbose)
    drv.load_atom_types(atom_types_file)
    drv.load_data(os.path.join(this_dir, args.xyz) if not os.path.isabs(args.xyz) else args.xyz)
    drv.load_dofs(os.path.join(this_dir, args.dof) if not os.path.isabs(args.dof) else args.dof)
    drv.init_and_upload()

    # Compile templated derivative and energy kernels
    macro_der = extract_macro_block(forces_path, args.model)
    macro_en  = extract_macro_block(forces_path, args.energy_model)
    drv.compile_with_model(macros={'MODEL_PAIR_ACCUMULATION': macro_der,
                                   'MODEL_PAIR_ENERGY':       macro_en}, bPrint=False)
    drv.set_energy_kernel_args()

    # Optionally disable regularization by zeroing stiffness terms and re-uploading
    if int(args.regularize) == 0 and hasattr(drv, 'host_regParams'):
        # Columns: [min,max, xlo,xhi, Klo,Khi, K0,x0] -> zero Klo,Khi,K0
        rp = drv.host_regParams.copy()
        rp[:, 4:7] = 0.0
        drv.host_regParams = rp
        if hasattr(drv, 'regParams_buff'):
            drv.toGPU_(drv.regParams_buff, drv.host_regParams)

    # DOF names and ranges
    mins, maxs, names_from_file = parse_dof_minmax(os.path.join(this_dir, args.dof) if not os.path.isabs(args.dof) else args.dof)
    dof_names = [f"{d['typename']}." + "REQH"[int(d['comp'])] for d in drv.dof_definitions]
    assert len(dof_names) == len(names_from_file)

    # select dofs
    exclude_set = set(args.exclude)
    iDOFs = [i for i,n in enumerate(dof_names) if n not in exclude_set]

    # baseline vector
    x0 = np.array([d['xstart'] for d in drv.dof_definitions], dtype=np.float32)

    diffs = []
    for i in iDOFs:
        mn, mx = mins[i]+args.eps, maxs[i]-args.eps
        if not np.isfinite(mn) or not np.isfinite(mx) or mn==mx:
            continue
        xs = np.linspace(mn, mx, args.points)
        Es = np.zeros_like(xs)  # J(x)
        Fs = np.zeros_like(xs)
        for j, x in enumerate(xs):
            xv = x0.copy(); xv[i] = x
            Es[j] = drv.evaluate_objective(xv)  # scalar objective J
            g = drv.get_forces(xv)
            Fs[j] = g[i]  # analytic dJ/dx
        Fnum, xs_num = num_deriv(xs, Es)
        # Compare analytic dJ/dx to numerical dJ/dx
        Fs_inner = Fs[1:-1]
        diff = float(np.max(np.abs(Fs_inner - Fnum))) if Fnum.size>0 else np.nan
        diffs.append((diff, i, dof_names[i]))

        # Text diagnostics per DOF
        Es_min, Es_max = float(np.nanmin(Es)), float(np.nanmax(Es))
        Fa_min, Fa_max = float(np.nanmin(Fs_inner)), float(np.nanmax(Fs_inner)) if Fs_inner.size>0 else (np.nan, np.nan)
        Fn_min, Fn_max = float(np.nanmin(Fnum)), float(np.nanmax(Fnum)) if Fnum.size>0 else (np.nan, np.nan)
        print(f"[{i:3d} {dof_names[i]:12s}] regularize={int(args.regularize)} | J(x): min {Es_min:.3e} max {Es_max:.3e} NaN? {np.isnan(Es).any()} | dJ/dx ana: min {Fa_min:.3e} max {Fa_max:.3e} NaN? {np.isnan(Fs_inner).any()} | dJ/dx num: min {Fn_min:.3e} max {Fn_max:.3e} NaN? {np.isnan(Fnum).any()} | Linf |Δ| {diff:.3e}")

        # Plots per DOF: (1) J(x), (2) dJ/dx analytic vs numerical
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
        ax1.plot(xs, Es, '-', lw=1.2)
        ax1.set_xlabel('DOF value')
        ax1.set_ylabel('J(x)')
        ax1.grid(True)
        ax2.plot(xs_num, Fs_inner, '-', lw=1.0, label='analytic')
        ax2.plot(xs_num, Fnum,     ':', lw=1.5, label='numeric')
        ax2.set_xlabel('DOF value')
        ax2.set_ylabel('dJ/dx')
        ax2.grid(True)
        ax2.legend(fontsize=8)
        fig.suptitle(f"{dof_names[i]} | regularize={int(args.regularize)}")

        if args.save:
            base = args.save if args.save.endswith('.') else (args.save+'_' if args.save else 'test_derivatives_ocl_')
            safe_name = dof_names[i].replace('/', '_').replace(' ', '_')
            fig.savefig(f"{base}{i:03d}_{safe_name}.png", dpi=150, bbox_inches='tight')

    # checks
    failed = [(d,i,n) for d,i,n in diffs if (np.isfinite(d) and d>args.tol) or (not np.isfinite(d))]
    for d,i,n in diffs:
        print(f"{i:3d} {n:12s} max|F_ana-F_num| = {d:.3e}")
    if failed:
        print("Failures:")
        for d,i,n in failed:
            print(f"  {i:3d} {n:12s} diff {d:.3e} > tol {args.tol}")

    if args.show or not args.save:
        plt.show()

    # non-zero exit to act as a test
    if failed:
        sys.exit(1)

