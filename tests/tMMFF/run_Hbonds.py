import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import MMFF as mmff


# =========== Setting

nmaxiter = 10000
#xyz_name = "butandiol-2"
xyz_name = "polymer-2_new-OH"


#ialg  = 2
#ialg  = 3
#alg_names=["GD","MDdamp","FIRE","FIREsmooth"]

#======== Body

#alg_name = alg_names[ialg]
mmff.setVerbosity( verbosity=2, idebug=0 )
mmff.init( xyz_name="data/"+xyz_name, bMMFF=True )     
#mmff.getBuffs()
#mmff.eval()

clrs = [ 'r','g','b','m','c' ]

rs = np.linspace(1.5,5.5, 400 )
nb = mmff.findHbonds( Rcut=4.0, angMax=30.0 );
for ib in range( nb ):
    c = clrs[ib]
    Es,Fs = mmff.sampleHbond( ib, rs );               plt.plot( rs, Es - Es[-1], c=c, ls='-',  label="Hb #"+str(ib)         )
    Es,Fs = mmff.sampleHbond( ib, rs, maskH=0.0 );    plt.plot( rs, Es - Es[-1], c=c, ls='--', label="Hb #"+str(ib)+" noHb" )
    Es,Fs = mmff.sampleHbond( ib, rs, maskQ=0.0 );    plt.plot( rs, Es - Es[-1], c=c, ls=':',  label="Hb #"+str(ib)+" noQ"  )

plt.legend()
plt.grid()
plt.ylim(-0.1,0.1)

plt.show()


print("ALL DONE")
#plt.show()