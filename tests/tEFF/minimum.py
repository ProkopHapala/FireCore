from ctypes import wstring_at
import sys
import numpy as np
import os
import time
import argparse

sys.path.append("../../")
from pyBall import eFF as eff

#elementPath = "export/scan_data/single_H2O-2.xyz"            #for H2O
#elementPath_e = "export/scan_data/single_H2O_ee2_relaxed.xyz"
elementPathMol2 = "EMolecules/H2O/.mol2"

#elementPath_e = "EMolecules/CO/e.xyz"          #for CO Does not work because of triple bond
#elementPath_e = "EMolecules/CO/ee2.xyz"

#elementPath_e = "EMolecules/CO2/e_relaxed.xyz"          #for CO2 does not work because of double bond
#elementPath_e = "EMolecules/CO2/ee2.xyz"

#elementPath_e = "export/scan_data/single_H2O_ee.xyz" #something older
#elementPath_e = "export/scan_data/H2O_10rnd_be.xyz"

#elementPath = "export/scan_data/single_CH4.xyz"             #for CH4
#elementPath_e = "export/scan_data/single_CH4_e.xyz" #1 CH4
#elementPath = "EMolecules/CH3/.mol2" 

#elementPath_e = "export/scan_data/CH4_10rnd_be.xyz" #10 CH4       

#elementPath = "EMolecules/Experiment.xyz"             #CH3CH3
#elementPath_e = "EMolecules/CH3CH3/E_ee2_relaxed.xyz"
#elementPath_e = "EMolecules/NH3/_e.xyz"

# elementPath_e = "export/scan_data/angdistscan_CH4_forces_fixed.xyz"
# elementPath_e = "export/scan_data/CH4_10rnd_be.xyz" #be stands for better electrons. They are positioned at hydrogen atoms

#elementPath = "testofMol2.mol2"
#elementPath_e = "testofMol2.xyz"
elementPath_e = "export/scan_data/double_H2O_ee2.xyz"                   #HydrogenBonds
#elementPath = "export/scan_data/double_H2O.xyz"
fileToSaveProcess = "processXYZ.xyz"


import shutil

if __name__ == "__main__":
    
    

    eff.setTrjName(fileToSaveProcess, savePerNsteps=100)
    if os.path.exists(fileToSaveProcess):
        os.remove(fileToSaveProcess) # deleting useless information

    eff.setVerbosity(3,2)
    #eff.processXYZ( elementPath, nstepMax=100000,ialg=2, dt=0.0001, Fconv=1e-7, bCoreElectrons=True) 
    #eff.processMol2( elementPathMol2,nstepMax=170000, dt=0.005,ialg=2, Fconv=1e-7, xyz_out=fileToSaveProcess)
    apos_ = np.zeros((8,3))
    eff.processXYZ_e( elementPath_e, nstepMax=200000, dt=0.01, optAlg = 2, Fconv=1e-7, apos=apos_, xyz_out=fileToSaveProcess) #, KRSrho=KRSrho 
    print(apos_)
    eff.info()
    dif = (apos_[0]-apos_[1])
    distance = np.sqrt(dif.dot(dif))
    print(distance)
    #exec(open("xyz_view_new-----.py").read())