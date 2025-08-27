#!/usr/bin/python3 -u

import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time

#sys.path.append("/home/niko/work/FIRECORE/FireCore/")
#sys.path.append("/home/prokop/git/FireCore-fitREQH")

# not line limit for numpy, not text wrap
np.set_printoptions(linewidth=1000, threshold=np.inf)


sys.path.append("../../")
from pyBall import FitREQ as fit
from pyBall.OCL.NonBondFitting import FittingDriver, extract_macro_block
from pyBall import atomicUtils as au

fit.plt = plt

# =====================
#      Functions
# =====================

def setup_cpu_fit(xyz_file, dof_file, morse=1, verbosity=1):
    """Initialize CPU-side FitREQ, load data, set weights, and return essentials."""
    fit.setVerbosity(verbosity, PrintDOFs=1, PrintfDOFs=1, PrintBeforReg=-1, PrintAfterReg=-1)
    fit.loadTypes()
    fit.loadDOFSelection(dof_file)
    dof_names, dof_specs = fit.loadDOFnames(dof_file, return_specs=True)
    nbatch = fit.loadXYZ(xyz_file, bAddEpairs=True, bOutXYZ=False)
    Erefs, x0s = fit.read_xyz_data(xyz_file)
    fit.setGlobalParams(kMorse=1.8, Lepairs=0.7)
    if morse:
        imodel = 2
        weights0, lens = fit.split_and_weight_curves(Erefs, x0s, n_before_min=100, weight_func=lambda E: fit.exp_weight_func(E, a=1.0, alpha=4.0))
    else:
        imodel = 3
        weights0, lens = fit.split_and_weight_curves(Erefs, x0s, n_before_min=2, weight_func=lambda E: fit.exp_weight_func(E, a=1.0, alpha=4.0))
    fit.setup(imodel=imodel, EvalJ=1, WriteJ=1, Regularize=-1)
    fit.setWeights(weights0)
    fit.getBuffs()
    fit.setFilter(EmodelCutStart=0.0, EmodelCut=0.5, PrintOverRepulsive=-1, DiscardOverRepulsive=-1, SaveOverRepulsive=-1, ListOverRepulsive=-1)
    # All DOFs by default
    exclude = set([])
    iDOFs = list(range(len(dof_names)))
    iDOFs_ = [i for i in iDOFs if dof_names[i] not in exclude]
    return {
        'imodel': imodel,
        'dof_names': dof_names,
        'dof_specs': dof_specs,
        'iDOFs': iDOFs_,
    }

def setup_gpu_driver(xyz_file, dof_file, model_macro, hb_gate=1, regularize=False, verbose=0, charge_from_type=False):
    """Create and initialize the OpenCL fitting driver with a selected model macro.
    charge_from_type selects runtime charge source (False=atoms.w, True=tREQHs[:,2]).
    """
    this_dir   = os.path.dirname(os.path.abspath(__file__))
    repo_root  = os.path.abspath(os.path.join(this_dir, '..', '..'))
    atom_types_file = os.path.join(repo_root, 'cpp', 'common_resources', 'AtomTypes.dat')
    forces_path     = os.path.join(repo_root, 'cpp', 'common_resources', 'cl', 'Forces.cl')

    fit_ocl = FittingDriver(verbose=verbose, use_type_charges=bool(charge_from_type))
    fit_ocl.load_atom_types(atom_types_file)
    fit_ocl.load_data(xyz_file if os.path.isabs(xyz_file) else os.path.join(this_dir, xyz_file))
    fit_ocl.load_dofs(dof_file if os.path.isabs(dof_file) else os.path.join(this_dir, dof_file))
    fit_ocl.init_and_upload()

    macro_der = extract_macro_block(forces_path, model_macro)
    macros = {
        'MODEL_PAIR_ACCUMULATION': macro_der,
        'HBOND_GATE_DEFINE': f"#define HBOND_GATE {int(hb_gate)}",
    }
    fit_ocl.compile_with_model(macros=macros, bPrint=False)
    fit_ocl.set_regularization_enabled(enabled=bool(regularize))
    return fit_ocl

def build_xs_for_dof(i, dof_specs, npts):
    dspec = dof_specs[i]
    xmin = dspec['min'] if dspec['min'] is not None else 0.0
    xmax = dspec['max'] if dspec['max'] is not None else 1.0
    return np.linspace(xmin+1e-6, xmax-1e-6, npts)

def scan_dof(backend, iDOF, xs, dof_names=None, fit_ocl=None, debug=False, print_every=0):
    """Unified DOF scan for both CPU and GPU backends.
    Returns {'xs','Es','Fs'}.
    """
    Es = np.zeros_like(xs, dtype=np.float64)
    Fs = np.zeros_like(xs, dtype=np.float64)
    if backend == 'cpu':
        y_backup = fit.DOFs[iDOF]
        Es[:], Fs[:] = fit.scanParam(iDOF, xs, bEvalSamples=True)
        fit.DOFs[iDOF] = y_backup
        return {'xs': xs, 'Es': Es, 'Fs': Fs}
    # GPU path
    assert fit_ocl is not None, 'fit_ocl is required for GPU backend'
    x0 = np.array([d['xstart'] for d in fit_ocl.dof_definitions], dtype=np.float32)
    for j, x in enumerate(xs):
        xv = x0.copy(); xv[iDOF] = float(x)
        if debug and (j == 0 or j == len(xs)//2 or j == len(xs)-1 or (print_every and (j % int(print_every) == 0))):
            print(f"[GPU] scan iDOF={iDOF} name={dof_names[iDOF] if dof_names else iDOF} x={x:.6g}")
        J, g = fit_ocl.getErrorDerivs(xv)
        Es[j] = J
        Fs[j] = g[iDOF]
    return {'xs': xs, 'Es': Es, 'Fs': Fs}

def run_scans_cpu(iDOFs, dof_specs, npts):
    scans = {}
    for i in iDOFs:
        xs = build_xs_for_dof(i, dof_specs, npts)
        scans[i] = scan_dof('cpu', iDOF=i, xs=xs)
    return scans

def run_scans_gpu(iDOFs, dof_specs, npts, fit_ocl, dof_names, debug=False):
    scans = {}
    for i in iDOFs:
        xs = build_xs_for_dof(i, dof_specs, npts)
        pe = max(1, len(xs)//2)
        scans[i] = scan_dof('gpu', iDOF=i, xs=xs, dof_names=dof_names, fit_ocl=fit_ocl, debug=debug, print_every=pe)
    return scans

def plot_overlays(scans_cpu, scans_gpu, dof_names, iDOFs):
    for i in iDOFs:
        fig = plt.figure(figsize=(8,10.0))
        axE = plt.subplot(2,1,1)
        axF = plt.subplot(2,1,2)
        print('-- DOF', i, dof_names[i])
        fit.plotDOFscan_one(i, DOFname=f"{dof_names[i]} [CPU]", bEs=True, bFs=True, verb=1, axE=axE, axF=axF, color='C0', data=scans_cpu[i])
        fit.plotDOFscan_one(i, DOFname=f"{dof_names[i]} [GPU]", bEs=True, bFs=True, verb=1, axE=axE, axF=axF, color='C1', data=scans_gpu[i])
        axE.legend(); axE.set_xlabel('DOF value'); axE.set_ylabel('E [kcal/mol]');   axE.grid(alpha=0.2)
        axF.legend(); axF.set_xlabel('DOF value'); axF.set_ylabel('F [kcal/mol/A]'); axF.grid(alpha=0.2)
        plt.suptitle(f"DOF scan: {dof_names[i]}")

if __name__ == '__main__':

    import argparse
    p = argparse.ArgumentParser(description='Compare CPU vs GPU DOF scans for FitREQ with optional runtime charge source switch')
    p.add_argument('--xyz', default='H2O_single.xyz')
    p.add_argument('--dof', default='dofSelection_H2O_Morse.dat')
    p.add_argument('--morse', type=int, default=1, help='1=Morse model, 0=LJ model')
    p.add_argument('--use_type_charges', type=int, default=1, help='GPU runtime charge source: 0=per-atom atoms.w, 1=type-based tREQHs[:,2]')
    p.add_argument('--npts', type=int, default=50)
    p.add_argument('--hb_gate', type=int, default=1)
    p.add_argument('--regularize', type=int, default=0)
    p.add_argument('--verbose', type=int, default=0)
    p.add_argument('--debug', action='store_true')
    args = p.parse_args()

    this_dir = os.path.dirname(os.path.abspath(__file__))
    xyz_abs = args.xyz if os.path.isabs(args.xyz) else os.path.join(this_dir, args.xyz)
    dof_abs = args.dof if os.path.isabs(args.dof) else os.path.join(this_dir, args.dof)
    # CPU
    cpu = setup_cpu_fit(xyz_abs, dof_abs, morse=int(args.morse), verbosity=args.verbose)
    # GPU
    model_macro = 'MODEL_MorseQ_PAIR' if int(args.morse) else 'MODEL_LJQH2_PAIR'
    fit_ocl = setup_gpu_driver(xyz_abs, dof_abs, model_macro=model_macro, hb_gate=int(args.hb_gate), regularize=bool(args.regularize), verbose=int(args.verbose), charge_from_type=bool(args.use_type_charges))
    # Scans
    fit.setVerbosity(0)
    scans_cpu = run_scans_cpu(cpu['iDOFs'], cpu['dof_specs'], int(args.npts))
    scans_gpu = run_scans_gpu(cpu['iDOFs'], cpu['dof_specs'], int(args.npts), fit_ocl, cpu['dof_names'], debug=bool(args.debug))
    # Plot
    plot_overlays(scans_cpu, scans_gpu, cpu['dof_names'], cpu['iDOFs'])
    plt.show()

