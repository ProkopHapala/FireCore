import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import MMFF as mmff


#======== Body

mmff.setVerbosity( verbosity=1, idebug=0 )
mmff.init( xyz_name="data/propandiol", bMMFF=True  )              # without MMFF
mmff.getBuffs()
mmff.eval()
mmff.setTrjName("relax.xyz",1)
mmff.relax(niter=100, Ftol=1e-6, bWriteTrj=False )

#plt.show()