from ctypes import wstring_at
import sys
import numpy as np
import os
import time
import argparse

sys.path.append("../../")
from pyBall import eFF as eff

elementPath = "export/scan_data/single_H2O-2.xyz"            #for H2O
#elementPath_e = "export/scan_data/single_H2O_ee2_relaxed.xyz"

#elementPath_e = "EMolecules/CO/e.xyz"          #for CO Does not work because of triple bond
#elementPath_e = "EMolecules/CO/ee2_relaxed.xyz"

#elementPath_e = "EMolecules/CO2/e_relaxed.xyz"          #for CO2 does not work because of double bond
#elementPath_e = "EMolecules/CO2/ee2_relaxed.xyz"

#elementPath_e = "export/scan_data/single_H2O_ee.xyz" #something older
#elementPath_e = "export/scan_data/H2O_10rnd_be.xyz"

#elementPath = "export/scan_data/single_CH4.xyz"             #for CH4
#elementPath_e = "export/scan_data/single_CH4_e.xyz" #1 CH4

#elementPath_e = "export/scan_data/CH4_10rnd_be.xyz" #10 CH4       

#elementPath = "EMolecules/Experiment.xyz"             #TomsExperiment
#elementPath_e = "EMolecules/CH3CH3/E_ee2.xyz"
#elementPath_e = "EMolecules/Amoniak/_e.xyz"

# elementPath_e = "export/scan_data/angdistscan_CH4_forces_fixed.xyz"
# elementPath_e = "export/scan_data/CH4_10rnd_be.xyz" #be stands for better electrons. They are positioned at hydrogen atoms


# elementPath_e = "H2O_spins_fc.xyz"
# elementPath = "export/scan_data/single_H2O.xyz"
fileToSaveProcess = "processXYZ.xyz"


import shutil

if __name__ == "__main__":
    
    

    eff.setTrjName(fileToSaveProcess, savePerNsteps=100)
    if os.path.exists(fileToSaveProcess):
        os.remove(fileToSaveProcess) # deleting useless information

    eff.setVerbosity(3,2)
    #eff.processXYZ( elementPath, nstepMax=0,ialg=2, dt=0.0001, Fconv=1e-7, bCoreElectrons=True) 
    #eff.processXYZ( elementPath, nstepMax=0,ialg=2, dt=0.0001, Fconv=1e-7); shutil.copy("processXYZ.xyz","processXYZ-.xyz")

    eff.processXYZ_e( "processXYZ-.xyz", nstepMax=5000000,optAlg=2, dt=0.0001, Fconv=1e-7) #, KRSrho=KRSrho 
    exec(open("xyz_view_new-----.py").read())