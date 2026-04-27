import sys
import os
import numpy as np
import argparse

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import MMFF_multi as mmff

def load_constraints(filename="constraints.txt"):
    constraints = []

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            # Skip comments and empty lines
            if line.startswith('#') or not line:
                continue

            parts = line.split()
            if len(parts) != 7:
                print(f"Warning: Invalid line in constraints file: {line}")
                continue

            atom_idx = int(parts[0])
            initial_pos = [float(parts[1]), float(parts[2]), float(parts[3])]
            final_pos = [float(parts[4]), float(parts[5]), float(parts[6])]
            constraints.append((initial_pos, final_pos))


    return constraints if constraints else None

# Thermodynamic Integration for Entropic Spring
def main():
    parser = argparse.ArgumentParser(description="Thermodynamic Integration for Entropic Spring")
    parser.add_argument("--nSys", type=int, default=10, help="Number of parallel systems")
    parser.add_argument("--xyz_name", type=str, default="../../cpp/common_resources/entropic_spring_20.xyz", help="Path to the molecule file")
    parser.add_argument("--nLambda", type=int, default=20, help="Number of Lambda windows")
    parser.add_argument("--nMDsteps", type=int, default=100000, help="Number of MD steps per window")
    parser.add_argument("--nEQsteps", type=int, default=5000, help="Number of equilibration steps per window")
    parser.add_argument("--Fconv", type=float, default=1e-6, help="Force convergence criterion")
    parser.add_argument("--constraints", type=str, default="constraints.txt", help="Path to constraints file")
    parser.add_argument("--surf_name", type=str, default="../../cpp/common_resources/xyz/NaCl_3x3_L3", help="Path to the surface file; use 'none' to disable the surface")
    parser.add_argument("--mode", type=str, default="TI", choices=["TI", "JE", "BOTH"], help="Mode of calculation: TI, JE, or BOTH")
    parser.add_argument("--K", type=float, default=5.0, help="Jarzynski Equality force constant")
    parser.add_argument("-T", "--temperature", type=float, default=300.0, help="Temperature in Kelvin")
    parser.add_argument("--dt", type=float, default=0.05, help="Time step")
    parser.add_argument("--t_damp", type=float, default=100.0, help="Damping time constant")
    parser.add_argument("--hard_atoms", action="store_true", help="Use hard atom constraints")
    parser.add_argument("--soft_atoms", action="store_true", help="Use soft atom constraints")
    parser.add_argument("--hard_dist",  action="store_true", help="Use hard distance constraints")
    parser.add_argument("--soft_dist",  action="store_true", help="Use soft distance constraints")

    args = parser.parse_args()
    surf_name = None if args.surf_name.lower() in ("none", "null", "off") else args.surf_name

    # Initialize MMFF_multi
    mmff.init(
        nSys_=args.nSys,
        xyz_name=args.xyz_name,
        surf_name=surf_name,
        sElementTypes="../../cpp/common_resources/ElementTypes.dat",
        sAtomTypes="../../cpp/common_resources/AtomTypes.dat",
        sBondTypes="../../cpp/common_resources/BondTypes.dat",
        sAngleTypes="../../cpp/common_resources/AngleTypes.dat",
        bMMFF=True,
        bEpairs=False,
        T=args.temperature,
        gamma=1/(args.t_damp*args.dt),
        GridFF=6
    )
    mmff.set_opt(
        dt_max=args.dt,
        dt_min=max(args.dt * 0.1, 1.0e-6),
        damp_max=0.2,
        finc=1.1,
        fdec=0.5,
        falpha=0.8,
        minLastNeg=5,
        cvf_min=-0.1,
        cvf_max=+0.1,
    )

    constraints = load_constraints(args.constraints)

    if constraints is None or len(constraints) == 0:
        print("ERROR: No constraints loaded!")
        sys.exit(1)

    nCVs = len(constraints)
    print(f"Loaded {nCVs} constraint(s) from {args.constraints}")

    # Flatten constraints into separate lists for initial and final positions
    initial_positions = []
    final_positions = []
    for i, (init_pos, final_pos) in enumerate(constraints):
        initial_positions.extend(init_pos)
        final_positions.extend(final_pos)
        print(f"  CV {i+1}: ({init_pos[0]:.1f}, {init_pos[1]:.1f}, {init_pos[2]:.1f}) → ({final_pos[0]:.1f}, {final_pos[1]:.1f}, {final_pos[2]:.1f})")

    print(f"\nParameters: nLambda={args.nLambda}, nMDsteps={args.nMDsteps}, nEQsteps={args.nEQsteps}, Mode={args.mode}\n")

    # Map mode to integer
    mode_map = {"TI": 0, "JE": 1, "BOTH": 2}
    mode_int = mode_map[args.mode.upper()]

    # Run free energy calculation
    result = mmff.computeFreeEnergy(
        nCVs=nCVs,
        initial_positions=initial_positions,
        final_positions=final_positions,
        nLambda=args.nLambda,
        nMDsteps=args.nMDsteps,
        nEQsteps=args.nEQsteps,
        Fconv=args.Fconv,
        mode=mode_int,
        K=args.K,
        hardAtoms = 1 if args.hard_atoms else -1,
        softAtoms = 1 if args.soft_atoms else -1,
        hardDist  = 1 if args.hard_dist  else -1,
        softDist  = 1 if args.soft_dist  else -1
    )
    print(f"Raw free energy change: {result:.6f} eV")

    out_base = os.path.splitext(os.path.basename(args.xyz_name))[0]
    out_data = f"{out_base}_free_energy.dat"
    print(f"\n{'=' * 60}")
    print(f"  Free energy change: {result:.6f} eV")
    print(f"  Results saved to: {out_data}")
    print(f"{'=' * 60}\n")

if __name__ == "__main__":
    main()
