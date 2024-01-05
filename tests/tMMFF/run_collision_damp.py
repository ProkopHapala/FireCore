import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import MMFF as mmff


#======== Body

mmff.setVerbosity( verbosity=1, idebug=0 )
#mmff.init( xyz_name="data/pyridine", surf_name="data/NaCl_1x1_L2" )                             # all
#mmff.init( xyz_name="data/pyridine", surf_name="data/NaCl_1x1_L2", bMMFF=False  )              # without MMFF
#mmff.init( xyz_name="data/pyridine", surf_name="data/NaCl_1x1_L2", bMMFF=False, gridStep=-1 )  # without gridFF

mmff.init( xyz_name="data/nHexadecan" )

#mmff.getBuffs()

nstepMax = 1000
outE  = np.zeros( nstepMax )
outF  = np.zeros( nstepMax )
outV = np.zeros( nstepMax )
outVF = np.zeros( nstepMax )

#mmff.eval()
#mmff.run(nstepMax=1000, dt=-1, Fconv=1e-6, ialg=2, outE=None, outF=None, omp=False)
nitr = mmff.run(nstepMax=1000, dt=-1, Fconv=1e-6, ialg=2, outE=outE, outF=outF, outV=outV, outVF=outVF, omp=False)

plt.plot( outE[:nitr], "-k", label="E"  )
plt.plot( outF[:nitr], "-r", label="F"  )
plt.plot( outV[:nitr], "-b", label="VF" )
plt.plot( outVF[:nitr],"-g", label="VF" )

#plt.xlim(0, 100)
plt.ylim(1e-6, 1e+6)
plt.yscale('log')
plt.legend()
plt.show()


#print( "FORCES:\n mmff.fapos:\n ", mmff.fapos )
#mmff.plot(bForce=True, Fscale=10.0 )
#plt.show()
exit(0)

