#!/usr/bin/env python3
"""
Test script for comparing CPU and GPU implementations of the Universal Force Field (UFF)
using a unified library interface.
"""

import sys
import os
import argparse
import numpy as np

# Add FireCore to the Python path
base_path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(base_path)

from pyBall import MMFF_multi as uff

# Define the path to common resource files
data_dir = os.path.join(base_path, "cpp/common_resources")

# ==================
#  Buffer Specs (mirrored from test_UFF_ocl.py)
# ==================
# Unified buffer specification to collect and print CPU buffers.
# Format:
# 'buffer_name': {
#     'stride': canonical columns (GPU layout),
#     'cpu_stride': optional CPU columns if different,
#     'type': 'int' (topology) or 'float' (parameters)
# }
BUF_SPECS = {
    # --- Topology Buffers (Integers) ---
    'bonAtoms':  {'stride': 2, 'type': 'int'},
    'angAtoms':  {'stride': 4, 'cpu_stride': 3, 'type': 'int'},
    'dihAtoms':  {'stride': 4, 'type': 'int'},
    'invAtoms':  {'stride': 4, 'type': 'int'},
    'neighs':    {'stride': 4, 'type': 'int'},
    'neighBs':   {'stride': 4, 'type': 'int'},
    # --- Parameter Buffers (Floats) ---
    'bonParams': {'stride': 2, 'type': 'float'},
    'angParams': {'stride': 5, 'type': 'float'},
    'dihParams': {'stride': 3, 'type': 'float'},
    'invParams': {'stride': 4, 'type': 'float'},
}

TOPOLOGY_SPECS = {k: v for k, v in BUF_SPECS.items() if v['type'] == 'int'}
PARAMS_SPECS   = {k: v for k, v in BUF_SPECS.items() if v['type'] == 'float'}

# ==================
#  Helper Functions
# ==================
def get_cpu_bufs(uff_obj, specs):
    """Collect CPU buffers from the given UFF object and pad to canonical stride when needed."""
    bufs = {}
    for name, spec in specs.items():
        if hasattr(uff_obj, name):
            buf = getattr(uff_obj, name)
            canonical_stride = spec['stride']
            cpu_stride = spec.get('cpu_stride', canonical_stride)
            if cpu_stride != canonical_stride:
                # e.g. angAtoms have 3 columns on CPU, 4 on GPU
                import numpy as _np
                padded = _np.full((buf.shape[0], canonical_stride), -1, dtype=_np.int32)
                padded[:, :cpu_stride] = buf
                bufs[name] = padded
            else:
                bufs[name] = buf
    return bufs

def print_bufs(bufs, title="Buffers"):
    print(f"--- {title} ---")
    if not bufs:
        print(" (No buffers to print)")
        return
    for name, buf in bufs.items():
        print(f"{name}\n", buf)

def run_uff(use_gpu, components, bPrintBufs=False, nPrintSetup=False):
    """
    Run a single UFF evaluation step on either CPU or GPU.

    Args:
        use_gpu (bool): If True, run on GPU (OpenCL), otherwise run on CPU.
        components (dict): A dictionary of flags to enable/disable UFF components.

    Returns:
        tuple: A tuple containing the calculated energy (placeholder) and forces.
    """
    print(f"\n--- Running UFF on {'GPU' if use_gpu else 'CPU'} ---")
    
    uff.setSwitches2(NonBonded=-1, SurfAtoms=-1, GridFF=-1)
    # Set which components to evaluate for this run
    uff.setSwitchesUFF(
        DoBond      =components.get( 'bonds',      -1 ),
        DoAngle     =components.get( 'angles',     -1 ),
        DoDihedral  =components.get( 'dihedrals',  -1 ),
        DoInversion =components.get( 'inversions', -1 ),
        DoAssemble          =  1,
        SubtractBondNonBond = -1,
        ClampNonBonded      = -1
    )

    print("py.DEBUG 4")
    
    # Select CPU or GPU path via the iParalel flag expected by C++ run() switch
    # CPU: 0 (serial) or 1 (OpenMP); GPU (OpenCL UFF): 2
    iParalel = 2 if use_gpu else 0

    print("py.DEBUG 5")
    
    #uff.run(nstepMax=1, iParalel=iParalel)

    uff.setTrjName("trj_multi.xyz", savePerNsteps=1)
    uff.run( nstepMax=1, dt=0.00, Fconv=1e-6, ialg=2, damping=0.1, iParalel=iParalel )
    # uff.run( nstepMax=1000, dt=0.02, Fconv=1e-6, ialg=2, damping=0.1, iParalel=iParalel )
    #uff.run( nstepMax=10000, dt=0.01, Fconv=1e-6, ialg=2, damping=0.1, iParalel=iParalel )
    #uff.run( nstepMax=1, dt=0.02, Fconv=1e-6, ialg=2, damping=0.1, iParalel=iParalel )
    print("py.DEBUG 6")
    energy = 0 
    if use_gpu:
        # Download results from GPU to host buffers; UFF exposes fapos (double) via init_buffers_UFF
        uff.download(bForces=True)
        forces = uff.fapos.copy()
    else:
        forces = uff.fapos.copy()
    
    return energy, forces

def compare_results(cpu_energy, cpu_forces, gpu_energy, gpu_forces, tol=1e-5, component_name=""):
    """Compare energy and forces from CPU and GPU, and report differences."""
    print(f"\n--- {component_name} Comparison ---")
    print(f"CPU Energy: {cpu_energy:.6f} | GPU Energy: {gpu_energy:.6f} (NOTE: GPU energy not fully implemented for comparison)")
    
    if cpu_forces.shape != gpu_forces.shape:
        print(f"ERROR: Force array shapes mismatch! CPU: {cpu_forces.shape}, GPU: {gpu_forces.shape}")
        return False

    force_diff = np.abs(cpu_forces - gpu_forces)
    max_force_diff = np.max(force_diff) if force_diff.size > 0 else 0.0
    
    print(f"Max force component difference: {max_force_diff:.6f}")
    
    if max_force_diff > tol:
        print(f"Validation FAILED for {component_name}: Forces differ more than tolerance.")
        print("CPU Forces:\n", cpu_forces)
        print("GPU Forces:\n", gpu_forces)
        print("Difference:\n", force_diff)
        return False
    else:
        print(f"Validation PASSED for {component_name}: Forces are consistent.")
        return True

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='UFF CPU vs. GPU Validation Test')
    default_mol = os.path.join(data_dir, 'mol', 'formic_acid.mol2')
    parser.add_argument('-m', '--molecule',  type=str, default=default_mol, help='Molecule file (.mol2, .xyz)')
    parser.add_argument('-t', '--tolerance', type=float, default=1e-5, help='Numerical tolerance for comparison')
    parser.add_argument('-p', '--print-buffers', type=int, default=0, help='Print buffer contents before run')
    parser.add_argument('-v', '--verbose',       type=int, default=0, help='Verbose output')
    args = parser.parse_args()

    # --- Initialize the library once ---
    print("--- Initializing MMFF_multi library ---")
    uff.init(
        nSys_=1,
        xyz_name=args.molecule,
        sElementTypes  = os.path.join(data_dir, "ElementTypes.dat"),
        sAtomTypes     = os.path.join(data_dir, "AtomTypes.dat"),
        sBondTypes     = os.path.join(data_dir, "BondTypes.dat"),
        sAngleTypes    = os.path.join(data_dir, "AngleTypes.dat"),
        sDihedralTypes = os.path.join(data_dir, "DihedralTypes.dat"),
        bMMFF=True, # to use UFF MMFF should be True!
        bUFF=True
    )
    # NOTE: setSwitches() does not accept the 'NonBonded' keyword; non-bonded handling for UFF
    # is controlled via setSwitchesUFF (e.g., ClampNonBonded, SubtractBondNonBond). We already
    # set ClampNonBonded/SubtractBondNonBond inside run_uff(), so we disable this invalid call.  # TODO
    # uff.setSwitches(NonBonded=-1)

    print("py.DEBUG 3")

    if args.print_buffers:
        print("--- Buffers Before Run ---")
        print("TOPOLOGY_SPECS ", TOPOLOGY_SPECS)
        topo_bufs = get_cpu_bufs(uff, TOPOLOGY_SPECS)
        print("PARAMS_SPECS ", PARAMS_SPECS)
        print_bufs(topo_bufs, "Topology Buffers")
        param_bufs = get_cpu_bufs(uff, PARAMS_SPECS)
        print_bufs(param_bufs, "Parameter Buffers")

    if args.verbose:
        uff.print_debugs() 
        uff.print_setup()

    uff.getBuffs_UFF()

    uff.print_debugs(bParams=False, bNeighs=False, bShifts=False, bAtoms=True)
    uff.apos[:,:] += 0.2 * np.random.randn(uff.natoms, 3)
    uff.print_debugs(bParams=False, bNeighs=False, bShifts=False, bAtoms=True)

    # Run CPU then GPU with the same initialization and switches
    #components = ['bonds', 'angles', 'dihedrals', 'inversions']
    #components = ['bonds']
    #components = ['bonds', 'angles']
    components = ['bonds', 'angles', 'dihedrals']
    #components = ['bonds', 'inversions']
    component_flags = {key: 1 for key in components }

    cpu_energy, cpu_forces = run_uff(use_gpu=False, components=component_flags)
    gpu_energy, gpu_forces = run_uff(use_gpu=True,  components=component_flags)

    # Single comparison across all enabled components
    passed = compare_results(cpu_energy, cpu_forces, gpu_energy, gpu_forces, tol=args.tolerance, component_name="ALL")
    print("\n================= SUMMARY =================")
    if passed:
        print("CPU vs GPU UFF comparison PASSED")
    else:
        print("CPU vs GPU UFF comparison FAILED")
    print("=========================================")
