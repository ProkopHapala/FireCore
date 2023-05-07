import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import MMFF_multi as mmff

#======== Body

mmff.setVerbosity( verbosity=1, idebug=0 )
#mmff.init( xyz_name="data/pyridine", surf_name="data/NaCl_1x1_L2r" )                             # all
mmff.init( xyz_name="data/pyridine", surf_name="data/NaCl_1x1_L2", bMMFF=False  )              # without MMFF
#mmff.init( xyz_name="data/pyridine", surf_name="data/NaCl_1x1_L2", bMMFF=False, gridStep=-1 )  # without gridFF
mmff.getBuffs()
mmff.eval()
#mmff.relax(1000)
print( "FORCES:\n mmff.fapos:\n ", mmff.fapos )
mmff.plot(bForce=True, Fscale=10.0 )
plt.show()
exit(0)