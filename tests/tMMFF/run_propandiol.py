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
xyz_name = "propandiol"


#ialg  = 2
ialg  = 3

alg_names=["GD","MDdamp","FIRE","FIREsmooth"]

#======== Body

alg_name = alg_names[ialg]
mmff.setVerbosity( verbosity=1, idebug=0 )
mmff.init( xyz_name="data/"+xyz_name, bMMFF=True  )     
mmff.getBuffs()
mmff.eval()

apos_bak = mmff.DOFs.copy()

damp_maxs = np.linspace( 0.05,0.3,50 )
ns = np.zeros( len(damp_maxs) )

T0 = time.time_ns()
for i,damp_max in enumerate(damp_maxs):
    #mmff.setTrjName("trj_%03i.xyz" %i ,100)
    mmff.DOFs [:,:] = apos_bak[:,:]
    mmff.vDOFs[:,:] = 0
    #mmff.set_opt( dt_max=0.1, dt_min=0.02, damp_max=0.2, finc=1.1, fdec=0.5, falpha=0.8, minLastNeg=5, cvf_min=-0.1, cvf_max=+0.1 )
    mmff.set_opt( damp_max=damp_max )
    ns[i] = mmff.run(nmaxiter, ialg=ialg )
T=time.time_ns()-T0; print(  "Time ", T*1.e-9, "[s]" )

plt.plot( damp_maxs, ns, '.-' )
plt.show()


'''

# ---------- Plot

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
'''

print("ALL DONE")
#plt.show()