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

def count_mask_lines( fgo_file, mask='#iconf' ):
    line = None
    nline = 0
    with open(fgo_file) as f:
        for l in f:
            if mask in l:
                line = l
                nline += 1
    return nline, line
                

def extract_nae( xyz_file ):
    with open(xyz_file) as f:
        ws = f.readline().strip().split()
        na= int(ws[0])
        ne= int(ws[1])
        return na, ne

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

def save_simulation(angleArr, distArr, allEtot, flexVar, variance, fileToSavePath):
    os.remove(fileToSavePath)
    with open(fileToSavePath, "w") as f:
        # Convert and flatten arrays (to ensure 1D), then write manually
        angleArrnp = np.array(angleArr).flatten()
        distArrnp = np.array(distArr).flatten()

        # Write angle and dist arrays in one line each
        f.write(" ".join(f"{val:.6f}" for val in angleArrnp) + "\n")
        f.write(" ".join(f"{val:.6f}" for val in distArrnp) + "\n")
        f.write("\n\n\n")  # 3 blank lines
        print(flexVar)

        for i in range(len(flexVar)):
            # Ensure each input is 1D
            print(flexVar)

            flex_row = np.array(flexVar[i]).flatten()
            print(flex_row)
            var_row = np.array(variance[i]).flatten()
            print(var_row)
            etot_row = np.array(allEtot[i]).flatten()
            print(etot_row)


            # Write each on a single line
            f.write(" ".join(f"{val:.6f}" for val in flex_row) + "\n")
            f.write(" ".join(f"{val:.6f}" for val in var_row) + "\n")
            f.write(" ".join(f"{val:.6f}" for val in etot_row) + "\n")
            f.write("\n")  # blank line between each block

if __name__ == "__main__":
    print("#=========== RUN /home/gabriel/git/FireCore/tests/tEFF/AI_rough.py")

    eff.setVerbosity(1,0)
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
    params, nrec = extract_blocks(elementPath)
    eff.initOpt( dt=0.005, damping=0.005, f_limit=1000.0)

    bCoreElectrons = False
    eff.setSwitches( coreCoul=1 )
    eff.preAllocateXYZ(elementPath, Rfac=-1.35, bCoreElectrons=bCoreElectrons )
    eff.getBuffs()
    eff.info()
    eff.esize[:]=0.7

    angleArr = params['ang']
    distArr = params['dist']
    # using lists for better performance
    allEtot = []
    flexVar = []
    variance = []

    for i in range(2):
        outEs = np.zeros((nrec,5))
        eff.processXYZ(elementPath, bOutXYZ=True, outEs=outEs, bCoreElectrons=bCoreElectrons, bChangeCore=False, bChangeEsize=True, nstepMax=10000, dt=0.005, Fconv=1e-3, ialg=2 );
        outEsdiff = np.zeros((nrec,5))
        outEsdiff[:,0] = params["Etot"] - outEs[:,0]
        outEsdiff[:,0] -= sum(outEsdiff[:,0])/len(outEsdiff[:,0])
        outEsdiff[:,0][params['dist'] == 0.7] = 0 #getting rid of unwanted valuess 
        var = np.var(outEsdiff[:,0])

        allEtot.append(outEsdiff[:,0])
        theta = np.array([i, i+1])
        print(theta)
        flexVar.append(theta)
        variance.append(var)
        print("======================================================================================================")
        print(i)
        print("======================================================================================================")


    save_simulation(angleArr, distArr, allEtot, flexVar, variance, fileToSavePath)

    print(angleArr)
    print(distArr)
    print(allEtot[0])
    print(len(allEtot[0]))
    print(len(angleArr))
    print(len(distArr))
    

    plot_energy_landscape( angleArr, distArr, allEtot[0])
    plt.title("Initial for theta = {flexVar[0]} and variance = {variance[0]}")

    # plot_energy_landscape( params['ang'], params['dist'], outEsdiff[:,0] )

    # plt.title("Difference between DFT and eFF")
    # plt.savefig("map2d_diff.png")

    print("#=========== DONE /home/gabriel/git/FireCore/tests/tEFF/AI_angdist_sim.py")
    plt.show()