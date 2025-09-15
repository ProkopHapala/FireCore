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

    uff.fapos[:, :] = 0.0 # make sure it's initialized

    uff.setTrjName("trj_multi.xyz", savePerNsteps=1)
    uff.run( nstepMax=1, dt=0.00, Fconv=1e-6, ialg=2, damping=0.1, iParalel=iParalel )
    # uff.run( nstepMax=1000, dt=0.02, Fconv=1e-6, ialg=2, damping=0.1, iParalel=iParalel )
    #uff.run( nstepMax=10000, dt=0.01, Fconv=1e-6, ialg=2, damping=0.1, iParalel=iParalel )
    #uff.run( nstepMax=1, dt=0.02, Fconv=1e-6, ialg=2, damping=0.1, iParalel=iParalel )
    print("py.DEBUG 6")
    energy = 0 
    if use_gpu:
        # Download results from GPU to host buffers; UFF exposes fapos (double) via init_buffers_UFF # make sure it's initialized
        forces = uff.fapos.copy()
    else:
        forces = uff.fapos.copy()
    
    return energy, forces


def scan_uff(nconf, nsys, components, tol=1e-3, seed=123):
    """Run multi-configuration CPU vs GPU force comparison using C++ scan().

    - Generates nconf configurations by perturbing the current uff.apos deterministically
    - CPU path evaluates sequentially
    - GPU path evaluates batched over nsys replicas
    """
    print(f"\n--- Running UFF scan for nconf={nconf} with nsys={nsys} ---")

    # Switches as in single run
    uff.setSwitches2(NonBonded=-1, SurfAtoms=-1, GridFF=-1)
    uff.setSwitchesUFF(
        DoBond      =components.get('bonds',      -1),
        DoAngle     =components.get('angles',     -1),
        DoDihedral  =components.get('dihedrals',  -1),
        DoInversion =components.get('inversions', -1),
        DoAssemble          =  1,
        SubtractBondNonBond = -1,
        ClampNonBonded      = -1,
    )

    # Build configurations from current geometry
    base = uff.apos.copy()
    rng = np.random.default_rng(seed)
    confs = np.repeat(base[None, :, :], nconf, axis=0)
    confs += 0.1 * rng.standard_normal(confs.shape)

    # CPU forces
    F_cpu = uff.scan(confs, iParalel=0)
    # GPU forces (batched)
    F_gpu = uff.scan(confs, iParalel=2)

    # Compare
    ok = True
    diffs = F_gpu - F_cpu
    max_abs = np.max(np.abs(diffs)) if diffs.size else 0.0
    idx = np.unravel_index(np.argmax(np.abs(diffs)), diffs.shape) if diffs.size else (0, 0, 0)
    print(f"scan(): max |ΔF| = {max_abs:.3e} at index {idx}")
    if max_abs > tol:
        ok = False
        print("scan(): Differences exceed tolerance. Showing summary per-config:")
        for ic in range(nconf):
            Fi = F_cpu[ic]
            Gi = F_gpu[ic]
            di = diffs[ic]
            ma = float(np.max(np.abs(di)))
            cpu_min = float(np.min(Fi)) if Fi.size else 0.0
            cpu_max = float(np.max(Fi)) if Fi.size else 0.0
            gpu_min = float(np.min(Gi)) if Gi.size else 0.0
            gpu_max = float(np.max(Gi)) if Gi.size else 0.0
            cpu_l2  = float(np.linalg.norm(Fi))
            gpu_l2  = float(np.linalg.norm(Gi))
            print("\n############################################")
            print(f"########### system {ic}  E_GPU=NA  E_CPU=NA")
            print(f"CPU_force stats: min={cpu_min:.6e} max={cpu_max:.6e} ||F||_2={cpu_l2:.6e}")
            print(Fi)
            print(f"GPU_force stats: min={gpu_min:.6e} max={gpu_max:.6e} ||F||_2={gpu_l2:.6e}")
            print(Gi)
            print(f"AbsDiff max|ΔF|={ma:.3e}")
            print("--------------------------------------------")
    else:
        print(f"scan(): PASSED within tol={tol:.2e}")
    return ok, F_cpu, F_gpu

def compare_results(cpu_energy, cpu_forces, gpu_energy, gpu_forces, tol=1e-5, component_name=""):
    """Rich comparison of CPU vs GPU forces and energy.

    - Prints min/max and norms for both CPU/GPU forces to catch all‑zero cases
    - Prints max absolute component difference and its location
    - On failure, prints full CPU/GPU force matrices and the absolute diff
    """
    print(f"\n--- {component_name} Comparison ---")
    print(f"CPU Energy: {cpu_energy:.6f} | GPU Energy: {gpu_energy:.6f} (NOTE: GPU energy not fully implemented for comparison)")

    if cpu_forces.shape != gpu_forces.shape:
        print(f"ERROR: Force array shapes mismatch! CPU: {cpu_forces.shape}, GPU: {gpu_forces.shape}")
        return False

    def stats(name, F):
        F = np.asarray(F)
        finite = np.isfinite(F)
        if not np.all(finite):
            print(f"WARNING: {name} has non‑finite values: NaN={np.isnan(F).sum()}, Inf={np.isinf(F).sum()}")
        minv = np.min(F) if F.size else 0.0
        maxv = np.max(F) if F.size else 0.0
        maxabs = np.max(np.abs(F)) if F.size else 0.0
        l2 = float(np.linalg.norm(F)) if F.size else 0.0
        print(f"{name} stats: min={minv:.6e} max={maxv:.6e} max|.|={maxabs:.6e} ||F||_2={l2:.6e}")
        return maxabs, l2

    cpu_maxabs, cpu_l2 = stats("CPU Forces", cpu_forces)
    gpu_maxabs, gpu_l2 = stats("GPU Forces", gpu_forces)

    # Detect all‑zero (or near‑zero) force arrays
    tiny = max(1e-12, tol*0.1)
    cpu_all_zero = cpu_maxabs < tiny
    gpu_all_zero = gpu_maxabs < tiny
    if cpu_all_zero and gpu_all_zero:
        print("NOTE: Both CPU and GPU forces are (near) zero; differences may be meaningless. Check inputs/switches.")

    # Differences
    diff = np.asarray(gpu_forces) - np.asarray(cpu_forces)
    adiff = np.abs(diff)
    max_force_diff = np.max(adiff) if adiff.size else 0.0
    max_idx = np.unravel_index(np.argmax(adiff), adiff.shape) if adiff.size else (0, 0)
    print(f"Max force component difference: {max_force_diff:.6e} at index {max_idx}")

    # Also print per‑atom max norms if over tolerance
    if max_force_diff > tol:
        atom_norm_diff = np.linalg.norm(diff.reshape(-1, 3), axis=1) if diff.ndim == 2 and diff.shape[1] >= 3 else np.linalg.norm(diff, axis=-1)
        worst_atom = int(np.argmax(atom_norm_diff)) if atom_norm_diff.size else -1
        if worst_atom >= 0:
            print(f"Worst atom index: {worst_atom}  |ΔF|={atom_norm_diff[worst_atom]:.6e}")

    # Decide pass/fail
    passed = max_force_diff <= tol

    if not passed:
        print(f"Validation FAILED for {component_name}: Forces differ more than tolerance = {tol:.2e}.")
        print("CPU Forces (N,3):\n", cpu_forces)
        print("GPU Forces (N,3):\n", gpu_forces)
        print("Abs Diff (N,3):\n", adiff)
        return False

    print(f"Validation PASSED for {component_name}: Forces are consistent within tol={tol:.2e}.")

    # If one side is all‑zero, highlight that despite PASS
    if cpu_all_zero or gpu_all_zero:
        which = "CPU" if cpu_all_zero else "GPU"
        print(f"WARNING: {which} forces are (near) zero: max|F| < {tiny:g}. Investigate why fapos is zero.")
        print("Hints: ensure DoAssemble=1, selected components are enabled, and the system/topology isn’t empty.")
    return True

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='UFF CPU vs. GPU Validation Test')
    default_mol = os.path.join(data_dir, 'mol', 'formic_acid.mol2')
    #default_mol = os.path.join(data_dir, 'mol', 'xylitol.mol2')
    parser.add_argument('-m', '--molecule',      type=str,   default=default_mol, help='Molecule file (.mol2, .xyz)')
    parser.add_argument('-t', '--tolerance',     type=float, default=1e-3, help='Numerical tolerance for comparison')
    parser.add_argument('-p', '--print-buffers', type=int,   default=0, help='Print buffer contents before run')
    parser.add_argument('-v', '--verbose',       type=int,   default=0, help='Verbose output')
    parser.add_argument('--nsys',     type=int, default=2, help='Number of GPU replicas (systems)')
    parser.add_argument('--nconf',    type=int, default=2, help='Number of configurations for scan(); 0 disables scan')
    parser.add_argument('--use-scan', type=int, default=1, help='Use scan() path (1) or single-step run() (0)')
    args = parser.parse_args()

    # --- Initialize the library once ---
    print("--- Initializing MMFF_multi library ---")
    uff.init(
        nSys_=args.nsys,
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
    #components = ['bonds',  'dihedrals']
    #components = ['bonds', 'inversions']
    #components = ['bonds', 'angles', 'dihedrals']
    components = ['bonds', 'angles', 'dihedrals', 'inversions']
    component_flags = {key: 1 for key in components }

    if args.use_scan or args.nconf>0:
        nconf = args.nconf if args.nconf>0 else args.nsys
        ok, F_cpu, F_gpu = scan_uff(nconf=nconf, nsys=args.nsys, components=component_flags, tol=args.tolerance)
        passed = ok
    else:
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
