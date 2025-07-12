from ctypes import wstring_at
import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time
import argparse

sys.path.append("../../")
from pyBall import eFF as eff
elementPath = "export/scan_data/single_CH4.xyz"
elementPath_e = "export/scan_data/single_CH4_e.xyz"
fileToSaveProcess = "processXYZ.xyz"
maxBars = 20

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

if __name__ == "__main__":
    print("#=========== RUN /home/gabriel/git/FireCore/tests/tEFF/Single_relax.py, all values")
    eff.setTrjName(fileToSaveProcess, savePerNsteps=1)
    if os.path.exists(fileToSaveProcess):
        os.remove(fileToSaveProcess) # deleting useless information

    KRSrho = np.array( [ 1.66487924,  1.61387271, -2.22015201])
    KRSrho = np.array([1.125, 0.9, -0.2])
    eff.setKRSrho(KRSrho)

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
    # params, nrec = extract_blocks("export/scan_data/angdistscan_CH4.xyz")
    params, nrec = extract_blocks(elementPath)
    print("nrec: ", nrec)
    outEs = np.zeros((nrec,5))
    # apos = np.zeros((nrec,,3))
    # epos = np.zeros((nrec,4))
    print("============================================================================")
    # with open("processXYZ.xyz", "w") as f: f.write("")
    #eff.processXYZ( "export/scan_data/distscan_H2O.xyz", bOutXYZ=True, outEs );
    print("opend process XYZ")
    eff.initOpt( dt=0.005, damping=0.005, f_limit=1000.0)

    bCoreElectrons = False
    eff.setSwitches( coreCoul=1 )
    #eff.setSwitches( coreCoul=0 )
    eff.preAllocateXYZ(elementPath, Rfac=-1.35, bCoreElectrons=bCoreElectrons )
    print("preallocate")
    eff.getBuffs()
    eff.info()
    print("get buffs")
    #eff.aPars[0,2]=1
    eff.esize[:]=0.7
    # eff.processXYZ( "export/scan_data/angdistscan_CH4.xyz", outEs=outEs, bCoreElectrons=bCoreElectrons, bChangeCore=False, bChangeEsize=True, nstepMax=0, dt=0.005, Fconv=1e-3, ialg=2 ) #, KRSrho=KRSrho 
    eff.processXYZ_e( elementPath_e, outEs=outEs, nstepMax=10000, dt=0.005, Fconv=1e-3) #, KRSrho=KRSrho 
    print("#=========== DONE /home/gabriel/git/FireCore/tests/tEFF/Single_relax.py, all values")
