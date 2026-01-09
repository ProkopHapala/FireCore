from ctypes import wstring_at
import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time
import argparse

sys.path.append("../../")
from pyBall import eFF as eff
elementPath = "export/scan_data/single_H2O.xyz"
# elementPath_e = "export/scan_data/single_H2O_ee.xyz"
# elementPath_e = "export/scan_data/H2O_10rnd_be.xyz"
elementPath_e = "export/scan_data/single_H2O_ee2.xyz"

# elementPath_e = "export/scan_data/angdistscan_CH4_forces_fixed.xyz"
# elementPath_e = "export/scan_data/CH4_10rnd_be.xyz" #be stands for better electrons. They are positioned at hydrogen atoms


# elementPath_e = "H2O_spins_fc.xyz"
# elementPath = "export/scan_data/single_H2O.xyz"
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
    eff.setTrjName(fileToSaveProcess, savePerNsteps=10)
    if os.path.exists(fileToSaveProcess):
        os.remove(fileToSaveProcess) # deleting useless information

    eff.setVerbosity(0,0)

    # KRSrho = np.array( [ 1.66487924,  1.61387271, -2.22015201])
    # KRSrho = np.array([1.125, 0.9, -0.2])
    # ECandPS =  [1.69444666e+02, 7.36370000e-02, 1.06690000e-02 ,1.57720000e-02] #wrong energy, right geometry
    # ECandPS = [26.183748553593603, 6.99412639825118, 3.457650280614723, 0.05897200521166494]
    # ECandPS = [66.4180837113872, 0.03500394691838782, 0.08530704650416196, 0.02056517014046732, 0.02986361488528383, 0.07165776124962407] # maybe right geometry
    ECandPS =  [1.835013 ,5.296729, 3.639857, 5.897393] #right energy, wrong geometry
    # ECandPS =  [2.764537 ,10.640034, 4.246970, 0.159267]
    # ECandPS =  [9.215005620624611, 0.29407744427742133, 0.34942230631929727, 0.32340019994012675]
    # ECandPS = [0,0,0,0]
    # KRSrho = np.array([ 1. ,  1. , -0.3])
    ECandPS =  [0.01503911370072064, 0.21527243058716336, 0.18234018105197047, 7.171216751018456]
    ECandPS =     [4651.225601 ,0.003989 ,0.009395 ,0.001288]

    # eff.setKRSrho(KRSrho)
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
    # params, nrec = extract_blocks("export/scan_data/angdistscan_CH4.xyz")
    # params, nrec = extract_blocks(elementPath)
    # print("nrec: ", nrec)
    # outEs = np.zeros((nrec,5))
    # apos = np.zeros((nrec,,3))
    # epos = np.zeros((nrec,4))
    # with open("processXYZ.xyz", "w") as f: f.write("")
    #eff.processXYZ( "export/scan_data/distscan_H2O.xyz", bOutXYZ=True, outEs );
    eff.initOpt( dt=0.005, damping=0.005, f_limit=1000.0)

    bCoreElectrons = False
    eff.setSwitches( coreCoul=1 )
    #eff.setSwitches( coreCoul=0 )
    eff.preAllocateXYZ(elementPath, Rfac=-1.35, bCoreElectrons=bCoreElectrons )
    eff.getBuffs()
    eff.info()
    #eff.aPars[0,2]=1
    eff.esize[:]=0.7
    nPos = 10
    na = 3
    apos  = np.zeros( (nPos, na, 3) )
    fapos = np.zeros( (nPos, na, 3) )
    outEs = np.zeros((nPos,5))
    # eff.processXYZ( "export/scan_data/angdistscan_CH4.xyz", outEs=outEs, bCoreElectrons=bCoreElectrons, bChangeCore=False, bChangeEsize=True, nstepMax=0, dt=0.005, Fconv=1e-3, ialg=2 ) #, KRSrho=KRSrho 
    
    eff.setFixedAtoms(False)
    eff.setParsECandPS(ECandPS) 
    eff.processXYZ_e( elementPath_e, outEs=outEs, apos=apos, fapos=fapos, nstepMax=10000, dt=0.005, Fconv=1e-3) #, KRSrho=KRSrho 

    for snglApos in apos:
        print("New molecule")
        angls = []
        bonds = []
        hatoms = len(snglApos) -1
        for i in range(hatoms):
            bond = np.linalg.norm(snglApos[i+1] - snglApos[0])
            bonds.append(bond)
            print("bond: ", bond)
            for j in range(i+1, hatoms):
                OH1 = snglApos[i+1] - snglApos[0]
                OH2 = snglApos[j+1] - snglApos[0]
                norm1 = np.linalg.norm(OH1)
                norm2 = np.linalg.norm(OH2)
                if norm1 > 1e-6 and norm2 > 1e-6:
                    cos_angle = np.dot(OH1, OH2) / (norm1 * norm2)
                    # Clip to handle floating point inaccuracies
                    # cos_angle = np.clip(cos_angle, -1.0, 1.0)
                    angl = np.arccos(cos_angle)
                else:
                    angl = 3.14 # if the norm is zero, set angle to pi (180 degrees)
                print("angle: ", angl)
                angls.append(angl)
    CH4angl = 1.9106 # in radians
    anglErr = np.sum((a-CH4angl)**2 for a in angls)

    forces = np.linalg.norm(fapos)
    print("Forces: ", forces)
    print("Angle error: ", anglErr)
    print("Bonds: ", bonds)
    print("#=========== DONE /home/gabriel/git/FireCore/tests/tEFF/Single_relax.py, all values")
