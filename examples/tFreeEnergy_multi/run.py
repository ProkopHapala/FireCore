import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import MMFF_multi as mmff

# Example usage of computeFreeEnergy function
def main():
    print("=== Free Energy Calculation Example ===")

    # Initialize the system
    nSys = 10
    xyz_name = "data/DA.mol2"  # Load the DA molecule pair

    print(f"Initializing MMFF_multi with {nSys} systems")
    print(f"Loading molecule from: {xyz_name}")

    # Parameters for free energy calculation
    lamda1 = 0.0      # Initial lambda value
    lamda2 = 1.0      # Final lambda value
    dc = [0, 1, 2]    # Distance constraint indices (example)
    nbStep = 100      # Number of steps
    nMDsteps = 100000 # Number of MD steps
    nEQsteps = 10000  # Number of equilibration steps
    dt = 0.5          # Time step
    tdamp = 100.0     # Temperature damping

    # Initialize MMFF_multi with the DA.mol2 file
    mmff.init(
        nSys_=nSys,
        xyz_name=xyz_name,
        sElementTypes="../../common_resources/ElementTypes.dat",
        sAtomTypes="../../common_resources/AtomTypes.dat",
        sBondTypes="../../common_resources/BondTypes.dat",
        sAngleTypes="../../common_resources/AngleTypes.dat",
        bMMFF=True,
        bEpairs=False,
        gamma = 1 / (dt * tdamp),     # Temperature damping
        T = 300.0         # Temperature in Kelvin

    )
    print("Initialization complete")



    # Call computeFreeEnergy
    print(f"Computing free energy with parameters:")
    print(f"  lamda1 = {lamda1}")
    print(f"  lamda2 = {lamda2}")
    print(f"  dc = {dc}")
    print(f"  nbStep = {nbStep}")
    print(f"  nMDsteps = {nMDsteps}")
    print(f"  nEQsteps = {nEQsteps}")
    print(f"  dt = {dt}")

    result = mmff.computeFreeEnergy(
        lamda1=lamda1,
        lamda2=lamda2,
        dc=dc,
        nbStep=nbStep,
        nMDsteps=nMDsteps,
        nEQsteps=nEQsteps,
        dt=dt
    )

    print(f"\nFree energy result: {result}")

if __name__ == "__main__":
    main()
