import sys
import numpy as np
import os
import time
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import MMFF_multi as mmff

#======== Body
def scanPlot( nscan = 1000, span=(0.0,8.0), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0.0), label="E" ):
    ts = np.linspace( span[0],span[1], nscan, endpoint=False)
    poss  = np.zeros( (nscan,3) )
    poss[:,0] = p0[0] + ts*dir[0]
    poss[:,1] = p0[1] + ts*dir[1]
    poss[:,2] = p0[2] + ts*dir[2]

    Es,Fs,Ps = mmff.scan( poss, bF=True, bP=True, Fconv=1e-5, niter_max=100000 )
    #print( "Es.shape ", Es.shape )
    plt.plot( ts, Es, '-', lw=0.5, label=label  )


#======== Body

mmff.setVerbosity( verbosity=1, idebug=1 )

#mmff.init( xyz_name="data/xyz/pyridine", surf_name="data/NaCl_1x1_L2" )    
#mmff.init( xyz_name="data/xyz/nHexadecan_dicarboxylic", bMMFF=True  )     
#mmff.init( xyz_name="data/xyz/O", surf_name="data/xyz/NaCl_1x1_L3" )  
mmff.init( xyz_name="data/xyz/H2O", surf_name="data/xyz/NaCl_1x1_L3", nSys_=10 )    
#mmff.init( xyz_name="data/xyz/PTCDA", surf_name="data/xyz/NaCl_1x1_L3" )    
mmff.getBuffs()
print("natoms=", mmff.natoms )
#print( "ffflags ", mmff.ffflags )

#mmff.setSwitches( NonBonded=-1, MMFF=-1, SurfAtoms=0, GridFF=1 )

#mmff.PLQs[:,2 ] = 0.0 # delete Coulomb (charges)
#mmff.PLQs[:,:2] = 0.0 # delete Morse (EvdW)
scanPlot( nscan=10, span=(0.0,8.0), dir=(1.0,0.0,0.0), p0=(1.0,0.0,0.0),  label="E_x" )
# scanPlot( nscan=1000, span=(0.0,8.0), dir=(0.0,1.0,0.0), p0=(0.0,0.0,0.0),  label="E_y" )
# #scanPlot( nscan=1000, span=(-5.0,5.0), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0.0), label="E_z" )

plt.legend()
plt.grid()
plt.show()
exit(0)
'''

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
    #mmff.run(10000,iParalel=-1)
    #mmff.run(10000,iParalel=0)
    #mmff.run(10000,iParalel=1)
    mmff.run(10000,iParalel=2)
t = time.time_ns()-t0;  print( "Py: time(optimizeLattice_1d) %g[s]" %(t*1e-9) )

#20,20, Mat3d{   0.0,0.5,0.0,    0.0,0.0,0.0,    0.0,0.0,0.0  }
#mmff.optimizeLattice_1d( [[0.0,0.5,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]], n1=20, n2=20, initMode=0, tol=1e-6)

#mmff.eval()
#mmff.relax(1000)
# print( "FORCES:\n mmff.fapos:\n ", mmff.fapos )
# mmff.plot(bForce=True, Fscale=10.0 )
# plt.show()
# exit(0)'
'''