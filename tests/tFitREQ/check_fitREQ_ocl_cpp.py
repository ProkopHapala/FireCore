#!/usr/bin/python3 -u

import sys
import numpy as np
import os
import time

#sys.path.append("/home/niko/work/FIRECORE/FireCore/")
#sys.path.append("/home/prokop/git/FireCore-fitREQH")

# not line limit for numpy, not text wrap
np.set_printoptions(linewidth=1000, threshold=np.inf)


sys.path.append("../../")
from pyBall import FitREQ as fit_cpp
#from pyBall.OCL.NonBondFitting import FittingDriver
from pyBall.OCL.FittingDriver import FittingDriver
from pyBall.OCL.OpenCLBase import OpenCLBase
from pyBall import atomicUtils as au

fit_cpp.plt = None  # Disable plotting

# =====================
#      Functions
# =====================

def setup_cpu_fit(args):
    """Initialize CPU-side FitREQ from argparse args, load data, and return config + handle."""
    # Resolve paths
    this_dir = os.path.dirname(os.path.abspath(__file__))
    xyz_file = args.xyz if os.path.isabs(args.xyz) else os.path.join(this_dir, args.xyz)
    dof_file = args.dof if os.path.isabs(args.dof) else os.path.join(this_dir, args.dof)

    # Core config
    fit_cpp.setVerbosity(int(args.verbose), PrintDOFs=1, PrintfDOFs=1, PrintBeforReg=-1, PrintAfterReg=-1)
    fit_cpp.loadTypes()
    fit_cpp.loadDOFSelection(dof_file)
    dof_names, dof_specs = fit_cpp.loadDOFnames(dof_file, return_specs=True)
    nbatch = fit_cpp.loadXYZ(xyz_file, bAddEpairs=False, bOutXYZ=False)
    Erefs, x0s = fit_cpp.read_xyz_data(xyz_file)
    fit_cpp.setGlobalParams(kMorse=1.8, Lepairs=0.7)

    # Model and weights
    morse = int(args.morse)
    if morse:
        imodel = 2
        weights0, lens = fit_cpp.split_and_weight_curves(Erefs, x0s, n_before_min=100, weight_func=lambda E: fit_cpp.exp_weight_func(E, a=1.0, alpha=4.0))
    else:
        imodel = 3
        weights0, lens = fit_cpp.split_and_weight_curves(Erefs, x0s, n_before_min=2,   weight_func=lambda E: fit_cpp.exp_weight_func(E, a=1.0, alpha=4.0))

    fit_cpp.setup(imodel=imodel, EvalJ=1, WriteJ=1, Regularize=-1, useTypeQ=args.use_type_charges*2-1)
    fit_cpp.setWeights(weights0)
    fit_cpp.getBuffs()
    fit_cpp.setFilter(EmodelCutStart=0.0, EmodelCut=0.5, PrintOverRepulsive=-1, DiscardOverRepulsive=-1, SaveOverRepulsive=-1, ListOverRepulsive=-1)

    # Match GPU charge-source option
    try:
        fit_cpp.useTypeQ = int(args.use_type_charges)
    except Exception:
        pass

    # All DOFs by default
    exclude = set([])
    iDOFs = list(range(len(dof_names)))
    iDOFs_ = [i for i in iDOFs if dof_names[i] not in exclude]
    return {
        'fit_cpp': fit_cpp,
        'imodel': imodel,
        'dof_names': dof_names,
        'dof_specs': dof_specs,
        'iDOFs': iDOFs_,
        'xyz_file': xyz_file,
        'dof_file': dof_file,
    }

def setup_gpu_driver(args):
    """Create and initialize the OpenCL fitting driver using argparse args."""
    # Resolve paths
    this_dir   = os.path.dirname(os.path.abspath(__file__))
    repo_root  = os.path.abspath(os.path.join(this_dir, '..', '..'))
    xyz_file   = args.xyz if os.path.isabs(args.xyz) else os.path.join(this_dir, args.xyz)
    dof_file   = args.dof if os.path.isabs(args.dof) else os.path.join(this_dir, args.dof)
    atom_types_file = os.path.join(repo_root, 'cpp', 'common_resources', 'AtomTypes.dat')
    forces_path     = os.path.join(repo_root, 'cpp', 'common_resources', 'cl', 'Forces.cl')

    # Driver
    fit_ocl = FittingDriver(verbose=int(args.verbose), use_type_charges=bool(args.use_type_charges))
    fit_ocl.load_atom_types(atom_types_file)
    fit_ocl.load_data(xyz_file)
    fit_ocl.load_dofs(dof_file)
    fit_ocl.init_and_upload()
    fit_ocl.serial_mode = (args.parallel == 0)


    # Model selection
    model_macro   = 'MODEL_MorseQ_PAIR' if int(args.morse) else 'MODEL_LJQH2_PAIR'
    model_macro_E = 'ENERGY_MorseQ_PAIR' if int(args.morse) else 'ENERGY_LJQH2_PAIR'
    
    # Parse CL file to get macros
    cl_parser = OpenCLBase()
    cl_content = cl_parser.parse_cl_lib(forces_path)
    macro_der   = cl_content['macros'].get(model_macro)
    macro_der_E = cl_content['macros'].get(model_macro_E)
    
    if macro_der is None or macro_der_E is None:
        raise RuntimeError(f"Could not find required macros in {forces_path}")
    
    macros = {
        'MODEL_PAIR_ACCUMULATION': macro_der,
        'MODEL_PAIR_ENERGY':       macro_der_E,
        'HBOND_GATE_DEFINE': f"#define HBOND_GATE {1}",
    }
    for k,v in macros.items():  print(f"{k}: {v}")
    #exit()
    output_path = os.path.join(this_dir, 'FitREQ_preprocessed.cl')
    fit_ocl.compile_with_model(macros=macros, bPrint=True, output_path=output_path)
    fit_ocl.load_program(output_path)
    return fit_ocl

def compare_objectives(cpu_fit, gpu_driver, dof_values, verbose=True):
    """Compare objective function values between CPU and GPU implementations."""
    
    # Set CPU DOF values
    for i, val in enumerate(dof_values):
        fit_cpp.DOFs[i] = float(val)
    
    # CPU objective: evaluate exactly once
    # Use scanParam to get the objective (it sets DOFs internally)
    cpu_result = fit_cpp.scanParam(0, np.array([dof_values[0]]), bEvalSamples=True)
    cpu_J = cpu_result[0][0]
    if verbose:
        print(f"CPU: evalFitError returned J: {cpu_J:.8e}")

    # GPU objective: evaluate exactly once via derivative path (applies DOFs to tREQHs)
    J, g = gpu_driver.getErrorDerivs(dof_values.astype(np.float32))
    gpu_J = J
    if verbose:
        print(f"GPU: J={J:.8e}, this should be sum over samples of 0.5*(Emol-Eref)^2")
    
    # Compute difference
    diff = abs(cpu_J - gpu_J)
    rel_diff = diff / (abs(cpu_J) + 1e-12)
    
    if verbose:
        print(f"DOF values: {dof_values}")
        print(f"CPU J: {cpu_J:.8e}")
        print(f"GPU J: {gpu_J:.8e}")
        print(f"Abs diff: {diff:.2e}")
        print(f"Rel diff: {rel_diff:.2e}")
    
    return {
        'cpu_J': cpu_J,
        'gpu_J': gpu_J,
        'abs_diff': diff,
        'rel_diff': rel_diff,
        'dof_values': dof_values.copy()
    }

if __name__ == '__main__':

    import argparse
    p = argparse.ArgumentParser(description='Compare CPU vs GPU objective function consistency for FitREQ')
    p.add_argument('--xyz', default='H2O_single.xyz')
    p.add_argument('--dof', default='dofSelection_H2O_Morse.dat')
    p.add_argument('--morse', type=int, default=1, help='1=Morse model, 0=LJ model')
    p.add_argument('--use_type_charges', type=int, default=1, help='GPU runtime charge source: 0=per-atom atoms.w, 1=type-based tREQHs[:,2]')
    p.add_argument('--uniform_weights', action='store_true', help='Use uniform weights (all 1.0) instead of exponential weights')
    p.add_argument('--verbose', type=int, default=1)
    p.add_argument('--parallel', type=int, default=0, help='Use serial mode (single lane per workgroup)')
    p.add_argument('--test_starting', action='store_true', help='Test using starting DOF values')
    p.add_argument('--test_custom', nargs='+', type=float, help='Test using custom DOF values')
    args = p.parse_args()

    this_dir = os.path.dirname(os.path.abspath(__file__))

    print("="*70)
    print("GPU vs CPU Consistency Test for FitREQ")
    print("="*70)
    print(f"XYZ file: {args.xyz}")
    print(f"DOF file: {args.dof}")
    print(f"Model: {'Morse' if args.morse else 'LJ'}")
    print(f"GPU charge source: {'type-based' if args.use_type_charges else 'per-atom'}")
    print("="*70)

    # CPU setup
    job_cpp = setup_cpu_fit(args)
    job_ocl = setup_gpu_driver(args)

    # Determine test DOF values
    if args.test_custom:
        test_dofs = np.array(args.test_custom, dtype=np.float64)
        print(f"Testing custom DOF values: {test_dofs}")
    elif args.test_starting:
        # Get starting values from DOF specs
        test_dofs = np.array([job_cpp['dof_specs'][i]['xstart'] for i in job_cpp['iDOFs']], dtype=np.float64)
        print(f"Testing starting DOF values: {test_dofs}")
    else:
        # Test with small random perturbations from starting values
        start_vals = np.array([job_cpp['dof_specs'][i]['xstart'] for i in job_cpp['iDOFs']], dtype=np.float64)
        np.random.seed(42)  # For reproducible results
        perturbations = np.random.normal(0, 0.1, len(start_vals))
        test_dofs = start_vals + perturbations
        print(f"Testing perturbed DOF values: {test_dofs}")
        print(f"  (base: {start_vals}, perturbations: {perturbations})")

    print("-"*70)

    # Run comparison
    print("Computing objective functions...")
    result = compare_objectives(job_cpp, job_ocl, test_dofs, verbose=True)

    print("-"*70)
    print("SUMMARY:")
    if result['rel_diff'] < 1e-6:
        print("✅ EXCELLENT: CPU and GPU results match within 1 ppm!")
    elif result['rel_diff'] < 1e-4:
        print("✅ GOOD: CPU and GPU results match within 0.01%!")
    elif result['rel_diff'] < 1e-2:
        print("⚠️  FAIR: CPU and GPU results differ by ~1%")
    else:
        print("❌ POOR: CPU and GPU results differ significantly")

    print(f"Relative error: {result['rel_diff']:.2e}")
    print("="*70)
