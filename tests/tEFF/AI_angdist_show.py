from ctypes import wstring_at
import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time
import argparse

sys.path.append("../../")
from pyBall import eFF as eff
elementPath = "export/scan_data/angdistscan_CH4.xyz"
fileToReadPath = "results/AI/result_3.txt"
maxBars = 20

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

def extract_blocks(xyz_file):
    """Extract parameters from XYZ file comments (lines starting with #)
    Returns:
        dict: Dictionary of extracted parameters with NaNs for missing values
        (e.g. {'ang': [...], 'dist': [...], 'Etot': [...]})
    """
    all_keys = set()
    records = []
    # First pass: collect all keys and raw records
    with open(xyz_file) as f:
        for line in f:
            if line.startswith('#'):
                parts = line[1:].strip().split()
                record = {}
                for i in range(0, len(parts)-1, 2):
                    key = parts[i]
                    val = float(parts[i+1])
                    record[key] = val
                    all_keys.add(key)
                records.append(record)
    # Initialize params with NaN-filled arrays
    nrec= len(records)
    params = {key: np.full(nrec, np.nan) for key in all_keys}
    # Second pass: fill values
    for i, record in enumerate(records):
        for key, val in record.items():
            params[key][i] = val
    return params, nrec

def allVal():
    print("#=========== RUN /home/gabriel/git/FireCore/tests/tEFF/AI_angdist_show.py, all values")
    print(f"Loading from file {fileToReadPath}")
    angleArr, distArr, flexVar, variance, allEtot = read_simulation(fileToReadPath)
    print(flexVar)
    variance = [x[0] for x in variance]
    labels = [f"({v[0]:.3f}, {v[1]:.3f}, {v[2]:.3f})" for v in flexVar]
    print(f"Total simulations: {len(flexVar)}")
    print(f"Minimum varieance: {min(variance)}")
    print(f"For KSrho: {flexVar[variance.index(min(variance))]}")
    print(f"In position {variance.index(min(variance))}")
    # Plot
    plt.figure(figsize=(15, 6))
    plt.bar(labels[-maxBars:], variance[-maxBars:])
    plt.xticks(rotation=90)
    plt.xlabel("Vector")
    plt.ylabel("Value")
    plt.title("Bar Plot of Values for Each Vector")
    plt.tight_layout()    
    print("#=========== DONE /home/gabriel/git/FireCore/tests/tEFF/AI_angdist_show.py")
    plt.show()


def minVal():
    print("#=========== RUN /home/gabriel/git/FireCore/tests/tEFF/AI_angdist_show.py, all values")
    print(f"Loading from file {fileToReadPath}")
    angleArr, distArr, flexVar, variance, allEtot = read_simulation(fileToReadPath)
    print(flexVar)
    variance = [x[0] for x in variance]
    index = variance.index(min(variance))
    KRSrho = flexVar[index]

    eff.setVerbosity(0,0)
    print("verbos")
    atomParams = np.array([
    #  Q   sQ   sP   cP
    [ 0.,  1.0, 1.00, 0.0 ], # 0
    [ 1.,  0.0, 0.00, 0.0 ], # 1 H
    [ 0.,  1.0, 1.00, 1.0 ], # 2 He
    [ 1.,  0.0, 0.10, 1.0 ], # 3 Li
    [ 2.,  0.0, 0.10, 1.0 ], # 4 Be
    [ 3.,  0.0, 0.10, 1.0 ], # 5 B
    [ 4.,  0.0, 0.10, 1.0 ], # 6 C
    [ 5.,  0.0, 0.10, 1.0 ], # 7 N
    [ 6.,  0.0, 0.2, 1.0 ], # 8 O
    [ 7.,  0.0, 0.10, 1.0 ], # 9 F
    ], dtype=np.float64)
    eff.setAtomParams( atomParams )
    print("set atom par")
    params, nrec = extract_blocks("export/scan_data/angdistscan_CH4.xyz")
    plot_energy_landscape( params['ang'], params['dist'], params['Etot'], Espan=5.0 )
    plt.title("Before relaxetion")
    plt.savefig("map2D_referece.png")
    print("plt.savefig")
    outEs = np.zeros((nrec,5))
    # apos = np.zeros((nrec,,3))
    # epos = np.zeros((nrec,4))
    print("============================================================================")
    with open("processXYZ.xyz", "w") as f: f.write("")
    #eff.processXYZ( "export/scan_data/distscan_H2O.xyz", bOutXYZ=True, outEs );
    print("opend process XYZ")
    eff.initOpt( dt=0.005, damping=0.005, f_limit=1000.0)

    bCoreElectrons = False
    eff.setSwitches( coreCoul=1 )
    #eff.setSwitches( coreCoul=0 )
    eff.preAllocateXYZ("export/scan_data/angdistscan_CH4.xyz", Rfac=-1.35, bCoreElectrons=bCoreElectrons )
    print("preallocate")
    eff.getBuffs()
    eff.info()
    print("get buffs")
    #eff.aPars[0,2]=1
    eff.esize[:]=0.7
    print("kRSrho: ", KRSrho)
    eff.setKRSrho(KRSrho)
    # eff.processXYZ( "export/scan_data/angdistscan_CH4.xyz", outEs=outEs, bCoreElectrons=bCoreElectrons, bChangeCore=False, bChangeEsize=True, nstepMax=0, dt=0.005, Fconv=1e-3, ialg=2 ) #, KRSrho=KRSrho 
    eff.processXYZ_e( "export/scan_data/angdistscan_CH4_e2.xyz", outEs=outEs, nstepMax=10000, dt=0.005, Fconv=1e-3) #, KRSrho=KRSrho 
    print("processXYZ")

    plot_energy_landscape( params['ang'], params['dist'], outEs[:,0] )
    plt.title("After relaxetion")
    plt.savefig("map2d_eFF.png")

    print("#=========== DONE /home/gabriel/git/FireCore/tests/tEFF/AI_angdist_show.py, all values")
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--variant', choices=['all', 'min'], required=True)
    args = parser.parse_args()

    if args.variant == 'all':
        allVal()
    else:
        minVal()    
