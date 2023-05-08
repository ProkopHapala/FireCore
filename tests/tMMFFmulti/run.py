import sys
import numpy as np
import os
import time
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import MMFF_multi as mmff

#======== Body


mmff.setVerbosity( verbosity=1, idebug=0 )
mmff.init( nSys_=40, xyz_name="data/polymer-2_new", surf_name="data/NaCl_1x1_L2"  )         # without MMFF
mmff.getBuffs()

#print(mmff.gpu_atoms)
#print(mmff.gpu_lvecs)
#print(mmff.gpu_lvecs.dtype)
#exit()

t0 = time.time_ns()
for i in range(5):
    mmff.change_lvec( [[0.2,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]], bAdd=True )
    mmff.run(10000,iParalel=-1)
    #mmff.run(10000,iParalel=0)
    #mmff.run(10000,iParalel=1)
    #mmff.run(1000,iParalel=2)
t = time.time_ns()-t0;  print( "Py: time(optimizeLattice_1d) %g[s]" %(t*1e-9) )

#20,20, Mat3d{   0.0,0.5,0.0,    0.0,0.0,0.0,    0.0,0.0,0.0  }
#mmff.optimizeLattice_1d( [[0.0,0.5,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]], n1=20, n2=20, initMode=0, tol=1e-6)

#mmff.eval()
#mmff.relax(1000)
# print( "FORCES:\n mmff.fapos:\n ", mmff.fapos )
# mmff.plot(bForce=True, Fscale=10.0 )
# plt.show()
# exit(0)