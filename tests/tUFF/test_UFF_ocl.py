#!/usr/bin/env python3
"""Simple UFF validation test comparing CPU and GPU implementations"""

import sys
import os
import argparse

# Add FireCore to path
base_path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(base_path)

data_dir = os.path.join(base_path, "cpp/common_resources")

import numpy as np
from pyBall.AtomicSystem import AtomicSystem
from pyBall import     MMFF as uff_cpp
from pyBall.OCL import UFF  as uff_ocl

# ==================
#  Buffer Specs
# ==================
# A single, unified definition for all buffers used in the comparison.
# This makes the code DRY and easier to maintain.
#
# Format:
# 'buffer_name': {
#     'stride': The canonical number of columns for this buffer (usually matches GPU).
#     'cpu_stride': (Optional) The number of columns in the CPU buffer if different.
#     'gpu_attr': (Optional) The attribute name(s) in the UFF_CL class.
#                   Can be a string for a single buffer or a tuple for composed buffers.
#     'type': 'int' for topology buffers, 'float' for parameter buffers.
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
    'angParams': {'stride': 5, 'gpu_attr': ('angParams1', 'angParams2_w'), 'type': 'float'},
    'dihParams': {'stride': 3, 'type': 'float'},
    'invParams': {'stride': 4, 'type': 'float'},
}

TOPOLOGY_SPECS = {k: v for k, v in BUF_SPECS.items() if v['type'] == 'int'}
PARAMS_SPECS   = {k: v for k, v in BUF_SPECS.items() if v['type'] == 'float'}

# ==================
#  Helper Functions
# ==================

def get_cpu_bufs(uff_cpp, specs):
    """
    Automatically collect all CPU buffers and pad them to match GPU layout where necessary.
    """
    bufs = {}
    for name, spec in specs.items():
        attr = name
        if hasattr(uff_cpp, attr):
            buf = getattr(uff_cpp, attr)
            # The canonical stride for comparison is the GPU stride
            canonical_stride = spec['stride']
            cpu_stride = spec.get('cpu_stride', canonical_stride)

            # If CPU stride is different, we need to pad it for comparison
            if cpu_stride != canonical_stride:
                # This is the case for angAtoms
                padded_buf = np.full((buf.shape[0], canonical_stride), -1, dtype=np.int32)
                padded_buf[:, :cpu_stride] = buf
                bufs[name] = padded_buf
            else:
                bufs[name] = buf
    return bufs

def get_gpu_bufs(uff_cl, specs):
    """Automatically collect all GPU buffers and reshape them."""
    bufs = {}
    for name, spec in specs.items():
        gpu_attr = spec.get('gpu_attr', name)
        stride = spec['stride']
        shape = (-1, stride)

        if isinstance(gpu_attr, tuple):  # Special case for composed buffers like angParams
            # This part is inherently specific to the composition
            buf1 = uff_cl.download_buf(gpu_attr[0]).reshape(-1, 4)
            buf2 = uff_cl.download_buf(gpu_attr[1]).reshape(-1, 1)
            bufs[name] = np.hstack((buf1, buf2))
        else:
            try:
                downloaded = uff_cl.download_buf(gpu_attr)
                if downloaded is not None:
                    bufs[name] = downloaded.reshape(shape)
            except (KeyError, AttributeError, ValueError):
                # Buffer might not exist or have wrong size, which is a valid state to check
                pass
    return bufs

def compare_bufs(cpu_bufs, gpu_bufs, tol=1e-6, buf_type='Buffers'):
    """
    Compare CPU and GPU buffers for consistency.
    If verbose, prints detailed side-by-side comparisons for any mismatched buffers.
    """
    error_messages = []
    has_errors = False

    # Use cpu_bufs keys as the reference set of buffers to check
    for name in cpu_bufs:
        if name not in gpu_bufs:
            # This is not necessarily an error if GPU doesn't need the buffer,
            # but for UFF validation it is.
            error_messages.append(f"--- ERROR in buffer: {name} ---\nGPU buffer not found or failed to download.\n")
            has_errors = True
            continue

        cpu = cpu_bufs[name]
        gpu = gpu_bufs[name]

        is_error = False
        error_msg = ""

        if cpu.shape != gpu.shape:
            is_error = True
            error_msg = f"Shape mismatch in '{name}': CPU {cpu.shape} vs GPU {gpu.shape}"
        elif np.issubdtype(cpu.dtype, np.floating):
            diff = np.abs(cpu - gpu)
            max_diff = np.max(diff) if diff.size > 0 else 0.0
            if max_diff > tol:
                is_error = True
                error_msg = f"Numerical diff in '{name}': max_diff={max_diff:.2e}"
        else:  # Integer types
            if not np.array_equal(cpu, gpu):
                is_error = True
                error_msg = f"Value mismatch in '{name}'"

        if is_error:
            has_errors = True
            details = f"--- ERROR in buffer: {name} ---\n{error_msg}\n"
            details += f"CPU buffer '{name}':\n{cpu}\n"
            details += f"GPU buffer '{name}':\n{gpu}\n"
            error_messages.append(details)

    if has_errors:
        print(f"\n========= {buf_type} Comparison Errors =========")
        for msg in error_messages:
            print(msg)
        print("==========================================")
        return False
    else:
        print(f"\nAll {buf_type} buffers match between CPU and GPU")
        return True

def print_bufs(bufs, title="Buffers"):
    """Generic function to print a dictionary of buffers with a title."""
    print(f"--- {title} ---")
    if not bufs:
        print(" (No buffers to print)")
        return
    for name, buf in bufs.items():
        print(f"{name}\n", buf)


def run_uff_cpp( args ):
    uff_cpp.init(
        xyz_name=args.molecule,
        sElementTypes=os.path.join(data_dir, "ElementTypes.dat"),
        sAtomTypes=os.path.join(data_dir, "AtomTypes.dat"),
        sBondTypes=os.path.join(data_dir, "BondTypes.dat"),
        sAngleTypes=os.path.join(data_dir, "AngleTypes.dat"),
        sDihedralTypes=os.path.join(data_dir, "DihedralTypes.dat"),
        bMMFF=True,
        bUFF=True,
        b141=True,
    )
    uff_cpp.setSwitches(NonBonded=-1, SurfAtoms=-1, GridFF=-1)
    #uff_cpp.setSwitchesUFF( DoAssemble=+1, DoBond=1, DoAngle=-1, DoDihedral=-1, DoInversion=-1,  SubtractBondNonBond=-1, ClampNonBonded=-1 )
    uff_cpp.setSwitchesUFF( DoAssemble=+1, DoBond=1, DoAngle=1, DoDihedral=1, DoInversion=1,  SubtractBondNonBond=-1, ClampNonBonded=-1 )
    print("-----------\n uff_cpp.getBuffs_UFF()  ")
    uff_cpp.getBuffs_UFF()
    print("-----------\n uff_cpp.setTrjName()  ")

    print("--- CPU Buffers Before Run ---")
    cpu_topo_bufs = get_cpu_bufs(uff_cpp, TOPOLOGY_SPECS)
    print_bufs(cpu_topo_bufs, "CPU Topology Buffers")
    cpu_param_bufs = get_cpu_bufs(uff_cpp, PARAMS_SPECS)
    print_bufs(cpu_param_bufs, "CPU Parameter Buffers")


    uff_cpp.print_debugs()
    uff_cpp.print_setup()

    uff_cpp.setTrjName("trj.xyz", savePerNsteps=1)
    #print("-----------\n uff_cpp.run()  ")
    #uff_cpp.run( nstepMax=1, dt=0.02, Fconv=1e-6, ialg=2, damping=0.1 )
    #uff_cpp.run( nstepMax=10000, dt=0.02, Fconv=1e-6, ialg=2, damping=0.1 )
    uff_cpp.run( nstepMax=10000, dt=0.01, Fconv=1e-6, ialg=2, damping=0.1 )
    energy = uff_cpp.Es[0]
    forces = uff_cpp.fapos.copy()
    return energy, forces, uff_cpp

def run_uff_ocl(args):
    mol = AtomicSystem(fname=args.molecule)

    # Initialize the OpenCL runner, which now handles the builder internally
    uff_cl = uff_ocl.UFF_CL(nloc=32)

    # Build the UFF topology, parameters, and allocate buffers
    uff_data = uff_cl.toUFF(mol)

    # Upload data to the GPU
    uff_cl.upload_topology_params(uff_data)
    uff_cl.upload_positions(mol.apos)

    # Parse components string into dictionary
    components = {
        'bonds': 'bonds' in args.components,
        'angles': 'angles' in args.components,
        'dihedrals': 'dihedrals' in args.components,
        'inversions': 'inversions' in args.components
    }

    uff_cl.bDoBonds      = components['bonds']
    uff_cl.bDoAngles     = components['angles']
    uff_cl.bDoDihedrals  = components['dihedrals']
    uff_cl.bDoInversions = components['inversions']

    print("--- GPU Buffers Before Run ---")
    gpu_topo_bufs = get_gpu_bufs(uff_cl, TOPOLOGY_SPECS)
    print_bufs(gpu_topo_bufs, "GPU Topology Buffers")
    gpu_param_bufs = get_gpu_bufs(uff_cl, PARAMS_SPECS)
    print_bufs(gpu_param_bufs, "GPU Parameter Buffers")

    uff_cl.run_eval_step(bClearForce=True)

    energy = uff_cl.get_total_energy()
    forces = uff_cl.get_forces()
    return energy, forces, uff_cl

def compare_results(cpu_energy, cpu_forces, gpu_energy, gpu_forces):
    print(f"CPU Energy: {cpu_energy}")
    print(f"GPU Energy: {gpu_energy}")
    print(f"Energy diff: {abs(cpu_energy - gpu_energy)}")

    force_diff = np.linalg.norm(cpu_forces - gpu_forces, axis=1)
    print(f"Max force diff: {np.max(force_diff)}")
    print(f"Mean force diff: {np.mean(force_diff)}")


'''
NOTE:
/home/prokop/git/FireCore/cpp/common_resources/mol/formic_acid.mol2
@<TRIPOS>MOLECULE
*****
 5 4 0 0 0
SMALL
GASTEIGER

@<TRIPOS>ATOM
      1 C           1.0554   -0.6718    0.0204 C.2     1  UNL1        0.2919
      2 O          -0.0556   -0.1659    0.0050 O.2     1  UNL1       -0.2548
      3 H           1.1556   -1.7510    0.0532 H       1  UNL1        0.1510
      4 O           2.1590    0.1028   -0.0031 O.3     1  UNL1       -0.4828
      5 H           2.0965    1.0681   -0.0324 H       1  UNL1        0.2948
@<TRIPOS>BOND
     1     1     2    2
     2     1     3    1
     3     1     4    1
     4     4     5    1
------
this means that bonds should be like:

bonAtoms
 [[0 1]
  [0 2]
  [0 3]
  [3 4]]



'''


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='UFF Validation Test Script')
    # Build absolute default path to avoid CWD dependence
    default_mol = os.path.join(data_dir, 'mol', 'formic_acid.mol2')
    parser.add_argument('-m','--molecule',      type=str,     default=default_mol,             help='Molecule file (.mol, .xyz)')
    parser.add_argument('-g','--gpu',           type=int,     default=1,                     help='Use GPU/OpenCL implementation' )
    parser.add_argument('-c','--components',    default='bonds,angles,dihedrals,inversions', help='Comma-separated list of UFF components to test')
    parser.add_argument('-v','--verbose',       action='count', default=1,                   help='Increase verbosity level'    )
    parser.add_argument('-t','--tolerance',     type=float,   default=1e-6,                  help='Energy comparison tolerance' )
    args = parser.parse_args()

    # --- we skip this for the moment
    cpu_energy, cpu_forces, uff_cpp = run_uff_cpp(args)  # DEBUG: uncomment this when GPU runs without error


    exit()

    if args.gpu:
        gpu_energy, gpu_forces, uff_cl = run_uff_ocl(args)


        #exit() # DEBUG: comment this when GPU runs without error
        compare_results(cpu_energy, cpu_forces, gpu_energy, gpu_forces)

        # New two-phase comparison
        cpu_topo_bufs = get_cpu_bufs(uff_cpp, TOPOLOGY_SPECS)
        gpu_topo_bufs = get_gpu_bufs(uff_cl,  TOPOLOGY_SPECS)

        # First, compare topology
        if compare_bufs(cpu_topo_bufs, gpu_topo_bufs, buf_type='Topology'):
            # If topology matches, compare parameters
            cpu_param_bufs = get_cpu_bufs(uff_cpp, PARAMS_SPECS)
            gpu_param_bufs = get_gpu_bufs(uff_cl,  PARAMS_SPECS)
            compare_bufs(cpu_param_bufs, gpu_param_bufs, buf_type='Parameters', tol=args.tolerance)
    else:
        print(f"CPU Energy: {cpu_energy}")
        print(f"CPU Forces:\n{cpu_forces}")

    print("=========\nALL DONE")
