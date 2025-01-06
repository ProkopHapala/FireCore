import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import atomicUtils as au
#from pyBall import MMFF as mmff
from pyBall import FFFit as fff
from pyBall import FireCore as fc
fff.fc = fc

#======== Body

path='./molecules'

# list all molecules in path
molecules = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]
print(molecules)

for name in molecules[:2]:
    print("#===========", name)
    #apos,iZs,enames,qs 
    #mol = au.AtomicSystem( os.path.join(path,name) )
    apos,iZs,enames,qs = au.loadAtomsNP( os.path.join(path,name) )
    #apos,iZs,enames,qs = au.loadAtomsNP( "C2H4.xyz" )
    #iZs = mol.atypes

    #iZs  = np.array( iZs,  dtype=np.int32,   order='C' ) # Ensure correct type and contiguous
    #apos = np.array( apos, dtype=np.float64, order='C' ) # Ensure correct type and contiguous

    print( "Before relaxation: apos=\n", apos )
    fc.initialize( atomType=iZs, atomPos=apos, verbosity=3 )
    fc.relax( apos, nstepf=1000, nmax_scf=100 )
    #mol.saveXYZ( "relax_%s.xyz" %name )
    print( "After relaxation: apos=\n", apos )
    fc.reload()