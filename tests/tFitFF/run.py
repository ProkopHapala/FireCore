import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import MMFF as mmff
from pyBall import FFFit as fff
from pyBall import FireCore as fc
fff.fc = fc

#======== Body

apos,iZs,enames,qs = au.loadAtomsNP( "C2H4.xyz" )
print( "Before relaxation: apos=\n", apos )

fc.initialize( atomType=iZs, atomPos=apos, verbosity=3 )
fc.relax( apos, forces=None, fixPos=[0,2,3], nstepf=1000, nmax_scf=100, Es=None )
fc.relax( apos, forces=None, fixPos=[0,2,3], nstepf=1000, nmax_scf=100, Es=None )
fc.relax( apos, forces=None, fixPos=[-1   ], nstepf=1000, nmax_scf=100, Es=None )
print( "After relaxation: apos=\n", apos )

sel=[1,4,5]

#print( "##### LINEAR SCAN " )
#apos[sel,0]+=-0.5
#Es = fff.makeLinearScan_firecore( 20, sel, [0.1,0.0,0.0], apos, nmax_scf=200 )

print( "##### ROTATION SCAN " )
Es = fff.makeRotationScan_firecore( 20, sel, fff.makeRotMat( np.pi/20 ), apos[0], apos, nmax_scf=200 )

print( Es )
plt.plot(Es, '.-')
plt.grid()
plt.show()