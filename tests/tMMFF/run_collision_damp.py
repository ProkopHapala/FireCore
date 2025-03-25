import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import MMFF as mmff


#======== Body

mmff.setVerbosity( verbosity=1, idebug=0 )

#fname = "nHexadecan_fold"
#fname = "pentacene_dimer"
#fname = "pentacene_cross2"
#fname = "pentacene_cross3"
#fname = "pentacene_cross_preopt"

#fname = "pentacene_cross_preopt"

#fname = "diamin_and_diether_C10"
fname = "hydropentacene_cross"

#mmff.init( xyz_name="data/nHexadecan" ); mmff.setTrjName( "opt_nHexadecan.xyz" );
mmff.init( xyz_name="data/xyz/"+fname ); mmff.setTrjName( "opt_"+fname+".xyz", savePerNsteps=10 );


nds  =  1
cdB=-1.0;cdA=-1.0;cdNB=-1.0;cmd=0.0; dt=-1.0 
#fconv = 1e-3
fconv = 1e-4
#fconv = 1e-6

#cdB  =  0.1
#cdA  =  0.1
#cdB  =  0.5
#cdA  =  0.3

#cdA  =  0.2
#cdB  =  0.5

cdA  =  0.2
cdB  =  0.7

#cmd = 0.01

#cmd = -0.02   
#cmd = -0.05   
#cmd = -0.01
#cmd = -0.01


#dt   =  0.05
#dt   =  0.08
dt    =  0.10
#dt   =  0.12
#dt   =  0.15
#dt   =  0.20

#cmd  = 0.0   # acceleration
nmin = 50
cos  = 0.5


#mmff.getBuffs()
mmff.setupCollisionDamping( nstep=nds, medium=cmd, cB=cdB, cA=cdA, cNB=cdNB )
mmff.setup_accel( nstep_acc_min=nmin, cos_vf_acc=cos )

nstepMax = 10000
#nstepMax = 20000
outE  = np.zeros( nstepMax )
outF  = np.zeros( nstepMax )
outV  = np.zeros( nstepMax )
outVF = np.zeros( nstepMax )


#mmff.eval()

#nitr = mmff.run(nstepMax=1000, dt=-1, Fconv=1e-6, outE=outE, outF=outF, outV=outV, outVF=outVF)
nitr = mmff.run(nstepMax=nstepMax, dt=dt, Fconv=fconv, outE=outE, outF=outF, outV=outV, outVF=outVF)

fig, ax1 = plt.subplots( figsize=(10,5) )

Emin = outE[:nitr].min(); print( "Emin ", Emin )

#plt.subplot(2,1,1)
#ax1.plot( outE [:nitr]-Emin,"-k", label="E"  )
ax1.plot( outF [:nitr],"-r", label="|f|" , lw=0.5 )
ax1.plot( outV [:nitr],"-b", label="|v|" , lw=0.5 )
ax1.set_xlabel('MD step')
ax1.set_ylabel('E, F, V')
ax1.set_ylim(1e-6, 1e+2)
ax1.set_xlim(0,nstepMax)
ax1.set_yscale('log')
plt.grid()
plt.legend()

ax2 = ax1.twinx()
ax2.plot( range(nitr), outVF[:nitr],"-g", label="cos(v,f)", lw=0.5 )
ax2.axhline( y=1.0, color='gray', linestyle='--', lw=0.5 )
ax2.axhline( y=0.5, color='gray', linestyle='--', lw=0.5 )
ax2.axhline( y=0.0, color='gray', linestyle='--', lw=0.5 )
ax2.axhline( y=-1.0, color='gray', linestyle='--', lw=0.5 )
ax2.set_ylim(-1.1,+4.1)



plt.legend( loc='upper left')

#print( "outE ", outE [:nitr] )
#plt.xlim(0, 100)

#spar = "_nds%i_cdB%6.3f_cdA%6.3f_cdNB%6.3f_cmd%6.3f_dt%6.3f"  %(nds,cdB,cdA,cdNB,cmd,dt)
spar = "_cdB%5.2f_cdA%5.2f_cmd%6.3f_nmin%i_cos%5.2f_dt%6.3f"  %(cdB,cdA,cmd,nmin,cos,dt)
plt.title( fname+spar )

#plt.savefig("opt_new_"+fname+spar+".png", bbox_inches='tight')
plt.savefig("opt_FIRE_"+fname+spar+".png", bbox_inches='tight')
plt.show()


#print( "FORCES:\n mmff.fapos:\n ", mmff.fapos )
#mmff.plot(bForce=True, Fscale=10.0 )
#plt.show()
exit(0)

