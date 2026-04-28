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

def setup_cpu_fit(args, xyz_file, dof_file ):
    """Initialize CPU-side FitREQ from argparse args, load data, and return config + handle."""
    # Resolve paths
    #this_dir = os.path.dirname(os.path.abspath(__file__))
    #xyz_file = args.xyz if os.path.isabs(args.xyz) else os.path.join(this_dir, args.xyz)
    #dof_file = args.dof if os.path.isabs(args.dof) else os.path.join(this_dir, args.dof)

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
    fit_cpp.useTypeQ = int(args.use_type_charges)
    return fit_cpp



def setup_gpu_driver(args, xyz_file, dof_file, ):
    """Create and initialize the OpenCL fitting driver using argparse args."""
    # Resolve paths
    #this_dir   = os.path.dirname(os.path.abspath(__file__))
    #repo_root  = os.path.abspath(os.path.join(this_dir, '..', '..'))
    #xyz_file   = args.xyz if os.path.isabs(args.xyz) else os.path.join(this_dir, args.xyz)
    #dof_file   = args.dof if os.path.isabs(args.dof) else os.path.join(this_dir, args.dof)

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
    cl_content = fit_ocl.parse_cl_lib(forces_path)
    macro_der   = cl_content['macros'].get(model_macro)
    macro_der_E = cl_content['macros'].get(model_macro_E)
    
    if macro_der is None or macro_der_E is None:
        raise RuntimeError(f"Could not find required macros in {forces_path}")
    
    macros = {
        'MODEL_PAIR_ACCUMULATION': macro_der,
        'MODEL_PAIR_ENERGY':       macro_der_E,
        'HBOND_GATE_DEFINE': f"#define HBOND_GATE {1}",
    }
    #for k,v in macros.items():  print(f"{k}: {v}")
    #exit()
    output_path = os.path.join(this_dir, 'FitREQ_preprocessed.cl')
    fit_ocl.compile_with_model(macros=macros, bPrint=True, output_path=output_path)
    fit_ocl.load_program(output_path)
    return fit_ocl

def compare_objectives(fit_cpp, fit_ocl, dof_values, verbose=True):
    """Compare objective function values between CPU and GPU implementations."""
    
    # Set CPU DOF values
    for i, val in enumerate(dof_values):
        fit_cpp.DOFs[i] = float(val)
    
    # CPU objective: evaluate exactly once
    # Use scanParam to get the objective (it sets DOFs internally)
    J_cpu, F_cpu = fit_cpp.scanParam(0, np.array([dof_values[0]]), bEvalSamples=True)
    print(f"py.CPU: evalFitError returned J: {J_cpu}  F: {F_cpu}")

    # GPU objective: evaluate exactly once via derivative path (applies DOFs to tREQHs)
    J_gpu, F_gpu = fit_ocl.getErrorDerivs(dof_values.astype(np.float32))
    print(f"py.GPU: evalFitError returned J: {J_gpu}  F: {F_gpu}")
    
    # Compute difference
    diff = abs(J_cpu - J_gpu)
    rel_diff = diff / (abs(J_cpu) + 1e-12)
    
    if verbose:
        print(f"DOF values: {dof_values}")
        print(f"CPU J: {J_cpu}  F {F_cpu}")
        print(f"GPU J: {J_gpu}  F {F_gpu}")
        #print(f"Abs diff:   {diff:.2e}")
        #print(f"Rel diff:   {rel_diff:.2e}")
        
    return {
        'cpu_J': J_cpu,
        'gpu_J': J_gpu,
        'abs_diff': diff,
        'rel_diff': rel_diff,
        'dof_values': dof_values.copy()
    }

def check_derivatives(cpu_fit_config, gpu_driver, dof_values, delta=1e-6, verbose=True, compare_cpu_gpu=True):
    """Check derivatives using numerical differentiation."""
    
    # Get objective and derivatives at the original point
    result0 = compare_objectives(cpu_fit_config, gpu_driver, dof_values, verbose=False)
    J0 = result0['gpu_J']
    g0 = gpu_driver.getErrorDerivs(dof_values.astype(np.float32))[1]

    if verbose:
      print(f"Analytical Derivatives GPU: {g0}")

    # CPU objective function setup
    fit_cpp = cpu_fit_config['fit_cpp']
    for i, val in enumerate(dof_values):
        fit_cpp.DOFs[i] = float(val)
    cpu_result0 = fit_cpp.scanParam(0, np.array([dof_values[0]]), bEvalSamples=True)
    cpu_J0 = cpu_result0[0][0]
    

    # Perturb each DOF and compute numerical derivative
    numerical_derivatives_gpu = np.zeros_like(dof_values)
    numerical_derivatives_cpu = np.zeros_like(dof_values)
    
    for i in range(len(dof_values)):
        dof_values_plus = dof_values.copy()
        dof_values_plus[i] += delta

        # CPU side
        fit_cpp_plus = cpu_fit_config['fit_cpp']
        for j, val in enumerate(dof_values_plus):
            fit_cpp_plus.DOFs[j] = float(val)
        cpu_result_plus = fit_cpp_plus.scanParam(0, np.array([dof_values_plus[0]]), bEvalSamples=True)
        cpu_J_plus = cpu_result_plus[0][0]
        numerical_derivatives_cpu[i] = (cpu_J_plus - cpu_J0) / delta

        # GPU side
        result_plus = compare_objectives(cpu_fit_config, gpu_driver, dof_values_plus, verbose=False)
        J_plus = result_plus['gpu_J']
        numerical_derivatives_gpu[i] = (J_plus - J0) / delta

    # Compare numerical and analytical derivatives
    diff = numerical_derivatives_gpu - g0 #g0 analytical deriv GPU
    rel_diff = diff / (np.abs(g0) + 1e-12)

    if verbose:
        print("--------------------------------------------------")
        print("Derivative Check:")

        for i in range(len(dof_values)):
            print(f"  DOF {i}:")
            print(f"    Numerical GPU: {numerical_derivatives_gpu[i]:.6e}")
            print(f"    Numerical CPU: {numerical_derivatives_cpu[i]:.6e}")
            print(f"    Analytical: {g0[i]:.6e}")

            print(f"    Abs Diff (GPU): {diff[i]:.2e}")
            print(f"    Rel Diff: {rel_diff[i]:.2e}")

            diff_cpu = numerical_derivatives_cpu[i] - g0[i]
            print(f"    Abs Diff (CPU): {diff_cpu:.2e}")

            if compare_cpu_gpu:
                diff_gpu_cpu = numerical_derivatives_gpu[i] - numerical_derivatives_cpu[i]
                print(f"    Diff (GPU-CPU): {diff_gpu_cpu:.2e}")
        print("--------------------------------------------------")

    return {'numerical': numerical_derivatives, 'analytical': g0, 
            'abs_diff': diff, 'rel_diff': rel_diff}

if __name__ == '__main__':

    import argparse
    p = argparse.ArgumentParser(description='Compare CPU vs GPU objective function consistency for FitREQ')
    p.add_argument('--xyz', default='H2O_single.xyz')
    p.add_argument('--dof',   default='dofSelection_H2O_Morse.dat')
    p.add_argument('--morse', type=int, default=1, help='1=Morse model, 0=LJ model')
    p.add_argument('--use_type_charges', type=int, default=1, help='GPU runtime charge source: 0=per-atom atoms.w, 1=type-based tREQHs[:,2]')
    p.add_argument('--uniform_weights', action='store_true', help='Use uniform weights (all 1.0) instead of exponential weights')
    p.add_argument('--verbose',   type=int, default=1 )
    p.add_argument('--derivs',    type=int, default=0,  help='Check derivatives using numerical differentiation')
    p.add_argument('--parallel',  type=int, default=1, help='Use serial mode (single lane per workgroup)')
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

    this_dir   = os.path.dirname(os.path.abspath(__file__))
    repo_root  = os.path.abspath(os.path.join(this_dir, '..', '..'))
    xyz_file   = args.xyz if os.path.isabs(args.xyz) else os.path.join(this_dir, args.xyz)
    dof_file   = args.dof if os.path.isabs(args.dof) else os.path.join(this_dir, args.dof)

    fit_cpp = setup_cpu_fit   (args, xyz_file, dof_file )
    fit_ocl = setup_gpu_driver(args, xyz_file, dof_file )


    dof_names, dof_specs = fit_cpp.loadDOFnames(dof_file, return_specs=True)
    exclude = set([])
    iDOFs = list(range(len(dof_names)))
    iDOFs_ = [i for i in iDOFs if dof_names[i] not in exclude]

    # Determine test DOF values
    if args.test_custom:
        test_dofs = np.array(args.test_custom, dtype=np.float64)
        print(f"Testing custom DOF values: {test_dofs}")
    elif args.test_starting:
        # Get starting values from DOF specs
        test_dofs = np.array([dof_specs[i]['xstart'] for i in iDOFs_], dtype=np.float64)
        print(f"Testing starting DOF values: {test_dofs}")
    else:
        # Test with small random perturbations from starting values
        start_vals = np.array([dof_specs[i]['xstart'] for i in iDOFs_], dtype=np.float64)
        np.random.seed(42)  # For reproducible results
        perturbations = np.random.normal(0, 0.1, len(start_vals))
        test_dofs = start_vals + perturbations
        print(f"Testing perturbed DOF values: {test_dofs}")
        print(f"  (base: {start_vals}, perturbations: {perturbations})")

    print("-"*70)

    # Run comparison
    print("Computing objective functions...")
    result = compare_objectives(fit_cpp, fit_ocl, test_dofs, verbose=True)

    # Optional derivative check
    if args.derivs > 0:
        check_derivatives(fit_cpp, fit_ocl, test_dofs, verbose=bool(args.verbose), compare_cpu_gpu=True)

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
