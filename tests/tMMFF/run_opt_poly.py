import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall  import atomicUtils as au
from pyBall.atomicUtils import AtomicSystem
from pyBall  import plotUtils   as plu


from pyBall import MMFF as mmff

#======== Functions

def getFromCommnet( fname, varname="Hbonds" ):
    f = open(fname); f.readline(); comment = f.readline()
    i0 = comment.find( varname )
    #print( comment[i0:] )
    return eval( comment[i0+len(varname)+1:] )
    #print( Hbonds['X'], Hbonds['Y'] )

def setHBondConstrains( fname ):
    Hbonds = getFromCommnet( fname+".xyz" )
    hbXs = Hbonds['X']
    hbYs = Hbonds['Y']
    if( len(hbXs) != len(hbYs) ): 
        print("ERROR: len(hbXs)(%i) != len(hbYs)(%i) => exit() " %(en(hbXs),len(hbYs)))
        exit()
    print( Hbonds )

    for i in range(len(hbXs)):
        mmff.addDistConstrain( hbXs[i], hbYs[i], lmin=1.5, lmax=1.7, kmin=0.0, kmax=1., flim=10.0, l=None, k=None, shift=(1.,0.,0.), bOldIndex=True )

#fname = "out/BB.HNH-h.NHO-hh"

#======== Main Body



#mmff.setVerbosity( verbosity=1, idebug=0 )
mmff.setVerbosity( verbosity=0, idebug=0 )
#mmff.init( xyz_name="data/pyridine", surf_name="data/NaCl_1x1_L2", bMMFF=False  )              # without MMFF

names = [ os.path.splitext(fname)[0] for fname in os.listdir("out") ]
print(names)
nstepMax=2000
outE = np.zeros(nstepMax)
outF = np.zeros(nstepMax)

for name in names:
    print("########### " + name )
    mmff.init( xyz_name="out/"+name  )              # without MMFF
    setHBondConstrains( "out/"+name )
    mmff.change_lvec( [[-5.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]], bAdd=True )

    #mmff.getBuffs()
   
    mmff.setTrjName( "relax_trjs/"+name+".xyz" )
    outE[:]=0;outF[:]=0
    mmff.run( nstepMax=nstepMax, outE=outE, outF=outF )

    mmff.clear()



    '''
    print( "outE ", outE )
    print( "outF ", outF )
    plt.plot( outE )
    plt.plot( np.log10(outF) )
    plt.show()
    '''


