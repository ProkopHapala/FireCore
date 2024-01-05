import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import MMFF as mmff


#======== Body

mmff.setVerbosity( verbosity=1, idebug=0 )

#mmff.init( xyz_name="data/nHexadecan" ); mmff.setTrjName( "opt_nHexadecan.xyz" );
mmff.init( xyz_name="data/nHexadecan_fold" ); mmff.setTrjName( "opt_nHexadecan_fold.xyz", savePerNsteps=10 );

fconv = 1e-3
nds  =  1
cdB  = -1.0
#cdB  = 0.1
cdNB = -1.0
cmd  =  0.01
dt   =  0.05

#mmff.getBuffs()
mmff.setupCollisionDamping( nds, cmd, cdB, cdNB )

nstepMax = 20000
outE  = np.zeros( nstepMax )
outF  = np.zeros( nstepMax )
outV  = np.zeros( nstepMax )
outVF = np.zeros( nstepMax )


#mmff.eval()

#nitr = mmff.run(nstepMax=1000, dt=-1, Fconv=1e-6, outE=outE, outF=outF, outV=outV, outVF=outVF)
nitr = mmff.run(nstepMax=nstepMax, dt=dt, Fconv=fconv, outE=outE, outF=outF, outV=outV, outVF=outVF)

fig, ax1 = plt.subplots()

Emin = outE.min()

#plt.subplot(2,1,1)
#ax1.plot( outE [:nitr]-Emin,"-k", label="E"  )
ax1.plot( outF [:nitr],"-r", label="|f|" , lw=0.5 )
ax1.plot( outV [:nitr],"-b", label="|v|" , lw=0.5 )
ax1.set_xlabel('MD step')
ax1.set_ylabel('E, F, V')
ax1.set_ylim(1e-5, 1e+2)
ax1.set_yscale('log')
plt.grid()
plt.legend()

ax2 = ax1.twinx()
ax2.plot( range(nitr), outVF[:nitr],"-g", label="cos(v,f)", lw=0.5 )
ax2.set_ylim(-1.1,+4.1)

plt.legend( loc='upper left')

#print( "outE ", outE [:nitr] )

#plt.xlim(0, 100)

plt.savefig("opt_nds%i_cdB%6.3f_cdNB%6.3f_cmd%6.3f_dt%6.3f.png" %(nds,cdB,cdNB,cmd,dt), bbox_inches='tight')
plt.show()


#print( "FORCES:\n mmff.fapos:\n ", mmff.fapos )
#mmff.plot(bForce=True, Fscale=10.0 )
#plt.show()
exit(0)

