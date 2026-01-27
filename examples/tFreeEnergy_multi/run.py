import sys
import numpy as np
import matplotlib.pyplot as plt
import argparse

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import MMFF_multi as mmff

# Example usage of computeFreeEnergy function
def main():
    parser = argparse.ArgumentParser(description="Free Energy Calculation Example")
    parser.add_argument("--nSys", type=int, default=50, help="Number of systems")
    parser.add_argument("--xyz_name", type=str, default="data/DA.mol2", help="Path to the molecule file")
    parser.add_argument("--nCV", type=float, default=1.0, help="Number of Collective Variables")
    parser.add_argument("--nLambda", type=int, default=10, help="Number of Lambda steps")
    parser.add_argument("--nbStep", type=int, default=100, help="Number of steps")
    parser.add_argument("--nMDsteps", type=int, default=100000, help="Number of MD steps")
    parser.add_argument("--nEQsteps", type=int, default=10000, help="Number of equilibration steps")
    parser.add_argument("--dt", type=float, default=0.5, help="Time step")
    parser.add_argument("--tdamp", type=float, default=100.0, help="Temperature damping")
    parser.add_argument("--T", type=float, default=300.0, help="Temperature in Kelvin")

    args = parser.parse_args()

    print("=== Free Energy Calculation Example ===")

    print(f"Initializing MMFF_multi with {args.nSys} systems")
    print(f"Loading molecule from: {args.xyz_name}")

    # Initialize MMFF_multi with the DA.mol2 file
    mmff.init(
        nSys_=args.nSys,
        xyz_name=args.xyz_name,
        sElementTypes="../../cpp/common_resources/ElementTypes.dat",
        sAtomTypes="../../cpp/common_resources/AtomTypes.dat",
        sBondTypes="../../cpp/common_resources/BondTypes.dat",
        sAngleTypes="../../cpp/common_resources/AngleTypes.dat",
        bMMFF=True,
        bEpairs=False,
        gamma = 1 / (args.dt * args.tdamp),     # Temperature damping
        T = args.T         # Temperature in Kelvin

    )
    print("Initialization complete")

    # Define positions for Si
    initial_pos_1 = [-10.0, 0.0, 0.0]
    final_pos_1   = [-60.0, 0.0, 0.0]
    initial_pos_2 = [+1.0, 0.0, 0.0]
    final_pos_2   = [+6.0, 0.0, 0.0]

    # Call computeFreeEnergy
    print(f"Computing free energy with parameters:")
    print(f"  nCV = {args.nCV}")
    print(f"  nLambda = {args.nLambda}")
    print(f"  nbStep = {args.nbStep}")
    print(f"  nMDsteps = {args.nMDsteps}")
    print(f"  nEQsteps = {args.nEQsteps}")
    print(f"  dt = {args.dt}")
    print(f"  T = {args.T}")

    result = mmff.computeFreeEnergy(
        nCV=args.nCV,
        initial_pos_1=initial_pos_1,
        final_pos_1=final_pos_1,
        initial_pos_2=initial_pos_2,
        final_pos_2=final_pos_2,
        nLambda=args.nLambda,
        nbStep=args.nbStep,
        nMDsteps=args.nMDsteps,
        nEQsteps=args.nEQsteps,
        tdamp=args.tdamp,
        T=args.T,
        dt=args.dt
    )

    print(f"\nFree energy result: {result}")

if __name__ == "__main__":
    main()
