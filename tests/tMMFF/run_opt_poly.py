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

def setHBondConstrains( fname, k=5., lmax=1.7, flim=10.0 ):
    Hbonds = getFromCommnet( fname+".xyz" )
    hbXs = Hbonds['X']
    hbYs = Hbonds['Y']
    if( len(hbXs) != len(hbYs) ): 
        print("ERROR: len(hbXs)(%i) != len(hbYs)(%i) => exit() " %(en(hbXs),len(hbYs)))
        exit()
    print( Hbonds )

    for i in range(len(hbXs)):
        mmff.addDistConstrain( hbXs[i], hbYs[i], lmin=1.5, lmax=lmax, kmin=0.0, kmax=k, flim=flim, l=None, k=None, shift=(1.,0.,0.), bOldIndex=True )

#fname = "out/BB.HNH-h.NHO-hh"

#======== Main Body



#mmff.setVerbosity( verbosity=2, idebug=0 )
#mmff.setVerbosity( verbosity=1, idebug=0 )
mmff.setVerbosity( verbosity=0, idebug=0 )
#mmff.init( xyz_name="data/pyridine", surf_name="data/NaCl_1x1_L2", bMMFF=False  )              # without MMFF

names = [ name for name in os.listdir("out") ]
names = [ name for name in names if ( "2x2" not in name ) ]
names = [ os.path.splitext(name)[0] for name in names if ( os.path.splitext(name)[1] != '.sh' ) ]
print(names)
#nstepMax=0
nstepMax=10000
outE = np.zeros(nstepMax)
outF = np.zeros(nstepMax)

#names = sys.argv[1].split(',')

'''
names = [
    'BB.HNH-hp.OHO-h_1',
    #'BB.HNH-hh.NHO-hp',
    'BB.HNH-hh.NHO-hp',
]
'''

#names = [  'BB.HNH-hp.OHO-h_1', 'BB.HNH-hh.NHO-hp' ]

#names = [ 'BB.HNH-hh.NHO-hh' ]
#names = [ 'BB.HNH-hp.OHO-h_1', 'BB.HNH-hh.NHO-hh' ]


for name in names:
    print("########### (%s)" %name )
    mmff.init( xyz_name="out/"+name  )              # without MMFF
    setHBondConstrains( "out/"+name )
    mmff.change_lvec( [[-3.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]], bAdd=True )

    #mmff.getBuffs()
   
    mmff.setTrjName( "relax_trjs_omp/"+name+".xyz", nPBC=(2,2,1), savePerNsteps=1000000 )
    #mmff.setTrjName( "relax_trjs_omp/"+name+".xyz", nPBC=(2,2,1), savePerNsteps=1 )
    outE[:]=0;outF[:]=0

    #mmff.print_debugs()
    #nsteps = mmff.run( nstepMax=nstepMax, outE=outE, outF=outF )  ;print("relaxation took %i  steps of(%i)" %(nsteps,nstepMax) )
    nsteps = mmff.run( nstepMax=nstepMax, outE=outE, outF=outF, omp=True )  ;print("relaxation took %i  steps of(%i)" %(nsteps,nstepMax) )

    mmff.clear()


