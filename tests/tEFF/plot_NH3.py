from ctypes import wstring_at
import sys
import numpy as np
import os
import time
import argparse
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import eFF as eff



evalfile = "EMolecules/NH3/.mol2"
fileToSaveProcess = "processMol2.xyz"

if __name__ == "__main__":
    end = 0.0
    _start = 0.5
    iterations = 10
    eff.setVerbosity(0,0)
    if os.path.exists(fileToSaveProcess):
        os.remove(fileToSaveProcess) # deleting useless information
    eff.setFixedAtoms(True)
    aforce = np.zeros((5, 3))
    outf = np.empty([0])
    apos  = np.zeros((5, 3))
    outpos = np.empty([0])
    for i in range(0, iterations):
        R_eff = _start + (end - _start)*i/iterations
        atomParams = np.array([
        # Z_nuc, R_eff, Zcore_eff,   PA,        PB,        PC,        PD,        PE
        [ 0.0,   1.0,      0.0,      0.0,       0.0,       0.0,       0.0,       0.0      ], # 0: Dummy
        [ 1.0,   0.0,      0.0,      0.0,       0.0,       0.0,       0.0,       0.0      ], # 1: H (Bare nucleus)
        [ 2.0,   0.1,      2.0,      0.0,       0.0,       0.0,       0.0,       0.0      ], # 2: He (Simple core: 1 pair (2e), radius 0.3. sQ=0.3, sP=0.3, cP=1.0)
        [ 3.0,   0.1,      2.0,      0.0,       0.0,       0.0,       0.0,       0.0      ], # 3: Li (Simple core: Z=3, sQ=0.5, sP=0.5, cP=1.0 for 1s2)
        [ 4.0,   0.1,      2.0,      0.0,       0.0,       0.0,       0.0,       0.0      ], # 4: Be (Simple core: Z=4, sQ=0.4, sP=0.4, cP=1.0 for 1s2)
        [ 5.0,   0.1,      2.0,      0.0,       0.0,       0.0,       0.0,       0.0      ], # 5: B  (Simple core: Z=5, sQ=0.35,sP=0.35,cP=1.0 for 1s2)
        [ 6.0,   0.179   , 2.0,      22.721015, 0.728733,  1.103199,  17.695345, 6.693621 ], # 6: C (ECP: Z_nuc=6, R_core=0.621, Z_core=2. p-type)
        [ 7.0,   R_eff,    2.0,      0.0,       0.0,       0.0,       0.0,       0.0      ], # 7: N (ECP: Z_nuc=7, R_core=0.0,   Z_core=2. p-type)
        [ 8.0,   0.3,      2.0,      25.080199, 0.331574,  1.276183,  12.910142, 3.189333 ], # 8: O (ECP: Z_nuc=8, R_core=0.167813, Z_core=2. p-type)
        [ 9.0,   0.3,      2.0,      0.0,       0.0,       0.0,       0.0,       0.0      ]  # 9: F (Simple core: Z=9, sQ=0.3, sP=0.3, cP=1.0 for 1s2)
        ], dtype=np.float64)
        eff.setAtomParams( atomParams , mode=2)
        if (i == iterations): 
            eff.setTrjName(fileToSaveProcess, savePerNsteps=1)
        eff.processMol2(evalfile, nstepMax=170000, dt=0.00002, Fconv=1e-7, ialg=2, aforce=aforce, apos=apos)
        fsum = np.sum(np.abs(aforce))
        outf = np.append(outf, fsum)
        dif  =  apos[0] - apos[4]  
        outpos = np.append(outpos, np.sqrt(dif.dot(dif)))

        
    #print("H1", H2[1])
    fig, axs = plt.subplots(2)
    x = np.linspace(_start, end, iterations)
    axs[0].plot(  x, outpos, label="Distance" )
    axs[1].plot(  x, outf, label="Force" )
    #plt.plot( x, H2[1], label="e1" )   
    #plt.plot( x, H2[2], label = "e2" )
    axs[0].set(xlabel='R_eff', ylabel='distance')
    axs[1].set(xlabel='R_eff', ylabel='force')
    axs[0].axhline(y=1.0925, color='r', linestyle='--')
    plt.grid() 
    plt.legend()
    plt.show()