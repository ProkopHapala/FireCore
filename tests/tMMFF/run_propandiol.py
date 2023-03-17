import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import MMFF as mmff




xyz_name = "butandiol-2"

#ialg  = 2
ialg  = 3

alg_names=["GD","MDdamp","FIRE","FIREsmooth"]

#======== Body

alg_name = alg_names[ialg]
mmff.setVerbosity( verbosity=1, idebug=0 )
mmff.init( xyz_name="data/"+xyz_name, bMMFF=True  )     
mmff.getBuffs()
mmff.eval()

nmaxiter = 10000
cos,f,v,dt,damp = mmff.setOptLog( nmaxiter )
mmff.setTrjName(xyz_name+"."+alg_name+".xyz",1)

#mmff.relax(niter=nmaxiter, Ftol=1e-6, bWriteTrj=True )
#nsteps = mmff.run(nmaxiter, ialg=3 )  # run with FIRE_smooth 
nsteps = mmff.run(nmaxiter, ialg=ialg )  # run with FIRE

plt.figure( figsize=(10,10))
plt.subplot(2,1,1)
plt.plot( cos [:nsteps-1], label="cos(v,f)" );
plt.legend(); plt.grid();  plt.ylim( -1., 1. );
plt.subplot(2,1,2)
plt.plot( f   [:nsteps-1], label="|f|"  );
plt.plot( v   [:nsteps-1], label="|v|"  );
plt.plot( dt  [:nsteps-1], label="dt"   );
plt.plot( damp[:nsteps-1], label="damp" );
plt.legend(); plt.grid(); plt.yscale('log');  plt.ylim( 1e-8, 1e+4 );
plt.savefig(xyz_name+"."+alg_name+".png", bbox_inches='tight')
plt.show()

#print(f,v,damp,cos,dt)

print("ALL DONE")
#plt.show()