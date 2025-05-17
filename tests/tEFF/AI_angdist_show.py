from ctypes import wstring_at
import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time

sys.path.append("../../")
from pyBall import eFF as eff
elementPath = "export/scan_data/angdistscan_CH4.xyz"
fileToSavePath = "results/AI/result_1.txt"

def plot_energy_landscape( Xs, Ys, Es, Espan=None):
    """Plot energy landscape from XYZ file using imshow (simple and robust)."""
    #params, nrec = extract_blocks(xyz_file)
    xs, idx = np.unique ( Xs, return_inverse=True )
    ys, idy  = np.unique( Ys, return_inverse=True )
    energy_grid = np.full((len(xs), len(ys)), np.nan)
    # Fill the grid
    for i in range(len(Xs)):
        energy_grid[ idx[i], idy[i]] = Es[i]
    # Plot
    plt.figure(figsize=(10,8))
    if Espan is not None: 
        vmin = np.nanmin(energy_grid)
        vmax = vmin + Espan
    else:
        vmin = None
        vmax = None
    plt.imshow(energy_grid, extent=[ys.min(), ys.max(), xs.max(), xs.min()], aspect='auto', cmap='inferno', vmin=vmin, vmax=vmax)
    plt.colorbar(label='Total Energy ')
    plt.xlabel('Distance (Å)')
    plt.ylabel('Angle (rad)')
    #plt.title('Potential Energy Surface')     


def read_simulation(fileToReadPath):
    with open(fileToReadPath, "r") as f:
        lines = f.readlines()

    # Strip whitespace
    lines = [line.strip() for line in lines if line.strip() != ""]

    # Read angleArr and distArr
    angleArr = np.array([float(x) for x in lines[0].split()])
    distArr = np.array([float(x) for x in lines[1].split()])

    # The rest are blocks of 3 lines each: flexVar, variance, allEtot
    flexVar = []
    variance = []
    allEtot = []

    # Start after angle and dist (which were at line 0 and 1)
    i = 2
    while i < len(lines):
        flex_line = np.array([float(x) for x in lines[i].split()])
        var_line = np.array([float(x) for x in lines[i + 1].split()])
        etot_line = np.array([float(x) for x in lines[i + 2].split()])

        flexVar.append(flex_line)
        variance.append(var_line)
        allEtot.append(etot_line)

        i += 3  # advance by 3 lines per block

    return angleArr, distArr, flexVar, variance, allEtot


if __name__ == "__main__":
    print("#=========== RUN /home/gabriel/git/FireCore/tests/tEFF/AI_angdist_show.py")
    fileToReadPath = "results/AI/result_1.txt"
    angleArr, distArr, flexVar, variance, allEtot = read_simulation(fileToReadPath)

    plot_energy_landscape( angleArr, distArr, allEtot[0])
    plt.title(f"Initial for theta = { flexVar[0]} and variance = {variance[0]}")

    plot_energy_landscape( angleArr, distArr, allEtot[0])
    plt.title(f"Final for theta = { flexVar[-1]} and variance = {variance[-1]}")
    plt.savefig("finalDiff.png")

    print("#=========== DONE /home/gabriel/git/FireCore/tests/tEFF/AI_angdist_show.py")
    plt.show()

