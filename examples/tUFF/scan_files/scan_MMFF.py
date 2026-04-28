# this file should emulate interactive scans    uff.init(
import sys
import os
import argparse

import numpy as np
import matplotlib.pyplot as plt

# Add FireCore to the Python path
base_path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.append(base_path)

from pyBall import MMFF_multi as uff
from functions import *

# Define the path to common resource files
data_dir = os.path.join(base_path, "cpp/common_resources")



def components_to_switches(components, *, bNonBonded=False, bGridFF=False):
    def resolve(name: str) -> int:
        val = components.get(name, -1)
        if val is None:
            return -1
        return int(val)

    DoBond      = resolve('bonds')
    DoAngle     = resolve('angles')
    DoDihedral  = resolve('dihedrals')
    DoInversion = resolve('inversions')
    DoAssemble  = 1 if any(v > 0 for v in (DoBond, DoAngle, DoDihedral, DoInversion)) else -1
    SubtractBondNonBond = 1
    ClampNonBonded      = 1
    uff.setSwitches2(
        NonBonded=1 if bNonBonded else -1,
        SurfAtoms=1 if bGridFF else -1,
        GridFF=1 if bGridFF else -1,
    )
    uff.setSwitches(
        doAngles=DoAngle,
        doPiPiI=DoDihedral,
        doPiSigma=DoInversion,
        doBonded=DoBond
    )
    return DoBond, DoAngle, DoDihedral, DoInversion, DoAssemble, SubtractBondNonBond, ClampNonBonded


def print_bufs(bufs, title="Buffers"):
    print(f"--- {title} ---")
    if not bufs:
        print(" (No buffers to print)")
        return
    for name, buf in bufs.items():
        print(f"{name}\n", buf)

def scan_uff(confs, components, bNonBonded=False, bGridFF=False):
    components_to_switches(components, bNonBonded=bNonBonded, bGridFF=bGridFF)
    print("components_to_switches() done")
    F_gpu = uff.scan(confs, iParalel=3)
    return F_gpu

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='UFF CPU vs. GPU Validation Test')
    #default_mol = os.path.join(data_dir, 'xyz', 'H2O.xyz')
    #default_mol = os.path.join(data_dir, 'mol', 'formic_acid.mol2')
    default_mol = os.path.join(data_dir, 'xyz', 'xylitol_WO_gridFF.xyz')
    #default_mol = os.path.join(data_dir, 'xyz', 'xylitol_for_gridFF.xyz')

    parser.add_argument('-m', '--molecule',      type=str,      default=default_mol, help='Molecule file (.mol2, .xyz)')
    parser.add_argument('-t', '--tolerance',     type=float,    default=1e-5,        help='Numerical tolerance for comparison')
    parser.add_argument('-p', '--print-buffers', type=int,      default=0,           help='Print buffer contents before run')
    parser.add_argument('-v', '--verbose',       type=int,      default=0,           help='Verbose output')
    parser.add_argument('--nsys',                type=int,      default=100,           help='Number of GPU replicas (systems)')
    parser.add_argument('--nconf',               type=int,      default=10000,           help='Number of configurations for scan(); 0 disables scan')
    parser.add_argument('--niter',               type=int,      default=3000,         help='Relaxation steps per configuration for scan_relaxed()')
    parser.add_argument('--fconv',               type=float,    default=1e-6,        help='Force convergence threshold for scan_relaxed()')
    parser.add_argument('--dt',                  type=float,    default=0.01,        help='Integration timestep for relaxation')
    parser.add_argument('--damping',             type=float,    default=0.1,         help='Damping for relaxation')
    parser.add_argument('--flim',                type=float,    default=1000.0,      help='Force limit for relaxation')
    parser.add_argument('--preset', choices=['all', 'bonded', 'bonded-only', 'none', 'grid-only', 'nonbonded-only'], default='all', help='Component preset before applying --enable/--disable/--comps.')
    parser.add_argument('--comps', nargs='*', choices=UFF_COMPONENTS, default=None, help='Explicit list of bonded UFF components to enable (overrides preset).')
    parser.add_argument('--enable', nargs='*', choices=UFF_COMPONENTS, default=None, help='Additional UFF components to enable.')
    parser.add_argument('--disable', nargs='*', choices=UFF_COMPONENTS, default=None, help='UFF components to disable.')
    parser.add_argument('--non-bonded', action='store_true', help='Enable non-bonded interactions.')
    parser.add_argument('--grid-ff', action='store_true', help='Enable GridFF interactions.')
    parser.add_argument('--surf', type=str, default=None, help='Surface file (.xyz)')
    parser.add_argument('--scan_atom', type=int, default=5, help='Atom to scan')
    parser.add_argument('--scan_dim', type=str, default='z', choices=['z', 'xy'], help='Scan dimension')
    parser.add_argument('--z_range', type=str, default='-8.0,-7.0,0.1', help='Scan range for z: "start,end,step" ')
    parser.add_argument('--xy_range', type=str, default='0.0,10.0,1.0,0.0,10.0,1.0', help='Scan range for xy: "x_start,x_end,x_step,y_start,y_end,y_step" ')
    parser.add_argument('--scan_pos', type=str, default='0.0,0.0', help='Scan position (e.g., for z-scan: "x,y"; for xy-scan: "z")')
    
    args = parser.parse_args()

    bNonBonded = args.non_bonded
    bGridFF    = args.grid_ff

    base_components = None if args.comps is None else list(args.comps)
    component_flags = build_component_flags(
        base_components,
        enable=args.enable,
        disable=args.disable,
        preset=args.preset,
    )

    # Cleanup old files before running
    cleanup_xyz_files(["scan_relaxed_cpu_*.xyz", "scan_relaxed_gpu_*.xyz", "trj_multi.xyz", "relaxed_gpu_*.xyz"])
    surf_name = os.path.join(data_dir, args.surf) if args.surf else None

    # --- Initialize the library once ---
    print("--- Initializing MMFF_multi library ---")
    uff.init(
        nSys_=args.nsys,
        xyz_name=args.molecule,
        surf_name=surf_name,
        sElementTypes  = os.path.join(data_dir, "ElementTypes.dat"),
        sAtomTypes     = os.path.join(data_dir, "AtomTypes.dat"),
        sBondTypes     = os.path.join(data_dir, "BondTypes.dat"),
        sAngleTypes    = os.path.join(data_dir, "AngleTypes.dat"),
        sDihedralTypes = os.path.join(data_dir, "DihedralTypes.dat"),
        bMMFF=True, # to use UFF MMFF should be True!
        bUFF=False,
        nExplore=1000,
        nRelax=500,
        T=300,
        gamma=0.1
    )


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

    uff.getBuffs()

    uff.print_debugs(bParams=False, bNeighs=False, bShifts=False, bAtoms=True)
    #uff.apos[:,:] += 0.2 * np.random.randn(uff.natoms, 3)
    uff.print_debugs(bParams=False, bNeighs=False, bShifts=False, bAtoms=True)

    if args.nconf > 0:
        print("scan_mmff()")
        confs, z_scan, xy_scan, nconf = generate_confs(
            uff.apos,
            args.nconf,
            123,
            args.scan_atom,
            args.scan_dim,
            args.z_range,
            args.xy_range,
            args.scan_pos
        )
        F_gpu = scan_uff(
            confs,
            components=component_flags,
            bNonBonded=bNonBonded,
            bGridFF=bGridFF,
        )

        if z_scan:
            scan_params = [float(x) for x in args.z_range.split(',')]
            start, end, step = scan_params
            Z = np.arange(start, end, step)
            Z_forces = F_gpu[:, args.scan_atom, 2]
            plt.figure(figsize=(8, 5))
            plt.plot(Z, Z_forces, 'o-')
            plt.xlabel('Z [Ang]')
            plt.ylabel('Fz [eV/Ang]')
            plt.grid(True)
            plt.tight_layout()
            plt.show()

        if xy_scan:
            scan_params = [float(x) for x in args.xy_range.split(',')]
            x_start, x_end, x_step, y_start, y_end, y_step = scan_params
            x_values = np.arange(x_start, x_end, x_step)
            y_values = np.arange(y_start, y_end, y_step)
            X, Y = np.meshgrid(x_values, y_values)
            Z_forces = F_gpu[:, args.scan_atom, 2]
            plt.figure(figsize=(8, 5))
            plt.contourf(X, Y, Z_forces.reshape(len(y_values), len(x_values)))
            plt.xlabel('X [Å]')
            plt.ylabel('Y [Å]')
            plt.colorbar(label='Fz [eV/Å]')
            plt.grid(True)
            plt.tight_layout()
            plt.show()

