import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import MMFF_multi as mmff

#======== Body

#./$name -m $nsys -x common_resources/BPBA   -g common_resources/NaCl_1x1_L2
#./$name -m 40 -x common_resources/polymer-2_new                  -g common_resources/NaCl_1x1_L2

mmff.setVerbosity( verbosity=1, idebug=0 )
mmff.init( xyz_name="data/polymer-2_new", surf_name="data/NaCl_1x1_L2"  )         # without MMFF
mmff.getBuffs()

print(mmff.gpu_atoms)
print(mmff.gpu_lvecs)
#print(mmff.gpu_lvecs.dtype)
exit()

#20,20, Mat3d{   0.0,0.5,0.0,    0.0,0.0,0.0,    0.0,0.0,0.0  }

mmff.optimizeLattice_1d( [[0.0,0.5,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]], n1=20, n2=20, initMode=0, tol=1e-6)

#mmff.eval()
#mmff.relax(1000)
# print( "FORCES:\n mmff.fapos:\n ", mmff.fapos )
# mmff.plot(bForce=True, Fscale=10.0 )
# plt.show()
# exit(0)