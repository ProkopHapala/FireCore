

import sys
import os

sys.path.append('../../')
from pyBall import atomicUtils as au
#from pyBall import FFFit as fff
#from pyBall import dftb_utils as dftbu
import subprocess

import numpy as np
import matplotlib.pyplot as plt
#from functools import partial

# ============= SETUP 

names=[
 "2_ammonia_ammonia",
 "2_formamide_formamide",
 "2_formicAcid_formicAcid", 
 "2_HCN_HCN",
 "2_water_water",
 "2_carbonicAcid_carbonicAcid",
 "2_cyanidenitrate_urea",
 "2_formicAcid_formamide",
 "2_nitro_diamine",
 "2_urea_urea",
 "2_water_ammonia",
 "2_water_carbonicAcid",
 "2_water_urea",
]

dir_h5d3   = "/home/prokophapala/git/FireCore/tests/dftb/HBsmall_scan"
dir_3ob    = "/home/prokophapala/git/FireCore/tests/dftb/HBsmall_scan_3ob"
dir_b3lyp  = "/home/prokophapala/git/FireCore/tests/tPsi4resp/HBsmall_scan"

# ============ Functions



# ============ MAIN

wd = os.getcwd()
for name in names:
    plt.figure()
    print( name )

    dat_h5d3  = np.genfromtxt( dir_h5d3  + ("/%s.dat" %name ) )
    dat_3ob  = np.genfromtxt( dir_3ob   + ("/%s.dat" %name ) )
    dat_b3lyp = np.genfromtxt( dir_b3lyp + ("/%s.dat" %name ) )

    plt.plot( dat_h5d3 [:,0], dat_h5d3 [:,1]-dat_h5d3 [-1,1], "g-", label="dftb3_h5d3" )
    plt.plot( dat_3ob  [:,0], dat_3ob  [:,1]-dat_3ob  [-1,1], "b-", label="dftb3_3ob"  )
    plt.plot( dat_b3lyp[:,0], dat_b3lyp[:,1]-dat_b3lyp[-1,1], "k-", label="b3lyp+d3/pVDZ, CP"   )
    plt.plot( dat_b3lyp[:,0], dat_b3lyp[:,2]-dat_b3lyp[-1,2], "k--",label="b3lyp+d3/pVDZ, noCP" )

    plt.ylim(-20.0,5.0)
    plt.xlabel( "x [A]" )
    plt.ylabel( "E [kcal/mol]" )
    plt.title( name )
    plt.legend()
    plt.grid()
    plt.savefig( name+".png", bbox_inches='tight')
plt.show()


