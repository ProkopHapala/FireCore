import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time

sys.path.append("../../")
from pyBall import eFF as eff


eff.load_fgo("data/H2O_shoot.fgo", True )
#eff.loadFromFile_fgo( "data/H2O_shoot.fgo", true );    # ifixe.push_back(ff.ne-1);      bMD=true; bRelaxed=true;
#eff.init()
eff.initOpt(0.1,0.1)
eff.getBuffs( )  #;exit()

print( "evel ", eff.evel  )

print( "aPars ", eff.aPars  )

#exit()
#eff.run(100, 0.001, ialg=-1 )

eff.setTrjName("my_new_beautiful_trj.xyz")

nsamp = 100
ie0 = -1
trj_e0 = np.zeros( (nsamp,3) )

for i in range(100):
    eff.run(10, 0.001, ialg=-1 )
    trj_e0[i,:] = eff.epos[ie0,:]
    #print( i,"\n", eff.epos, eff.evel )

#print(trj_e0)

plt.plot( trj_e0[:,1], trj_e0[:,2], '.-' )
plt.axis('equal'); plt.grid()
plt.show()

