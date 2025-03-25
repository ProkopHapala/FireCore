import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import MMFF as mmff


#======= Functions
def make_scan(  dct, key, vals, bPrint=True ):
    print( dct )
    results = []
    for i,val in enumerate( vals ):
        dct[key] = val
        if(bPrint): print("\n\n#============= run %i %s=%g =============\n\n" %(i,key,val) )
        if i==2: mmff.setVerbosity( verbosity=0, idebug=1 )
        mmff.init( xyz_name="data/xyz/"+fname );
        mmff.setup_accel( nstep_acc_min=dct['nmin'], cos_vf_acc=dct['cos'] )
        mmff.setupCollisionDamping( nstep=1, medium=dct['cmd'], cB=dct['cdB'], cA=dct['cdA'], cNB=dct['cdNB'] )
        nitr = mmff.run(nstepMax=nstepMax, dt=dct['dt'], Fconv=dct['fconv'] )
        if(bPrint): print( "dt = ", dt, "nitr = ", nitr )
        if nitr<0: nitr = np.nan
        results.append( nitr )
        mmff.clear()
    return results

def makeParams():
    params={
    'cdA'  : cdA,
    'cdB'  : cdB,
    'cdNB' : cdNB,
    'cmd'  : cmd,
    'dt'   : dt,
    'fconv': fconv,
    'nmin' : nmin,
    'cos'  : cos
    }   
    return params

#======== Body

mmff.setVerbosity( verbosity=0, idebug=0 )

#fname = "nHexadecan_fold"
#fname = "pentacene_dimer"
#fname = "pentacene_cross2"
#fname = "pentacene_cross3"
#fname = "pentacene_cross_preopt"

#fname = "pentacene_cross_preopt"

#fname = "diamin_and_diether_C10"
fname = "hydropentacene_cross"

#fname = "H2O"

#mmff.init( xyz_name="data/nHexadecan" ); mmff.setTrjName( "opt_nHexadecan.xyz" );
#mmff.init( xyz_name="data/"+fname ); mmff.setTrjName( "opt_"+fname+".xyz", savePerNsteps=10 );

cdB=-1.0;cdA=-1.0;cdNB=-1.0;cmd=0.0; dt=-1.0; nmin=50; cos=0.5
#fconv = 1e-3
fconv = 1e-4
#fconv = 1e-6
#cdB  =  0.1
#cdA  =  0.1
#cdB  =  0.5
#cdA  =  0.3

#cdA  =  0.2
#cdB  =  0.5
#cdA  =  0.3
#cdB  =  0.5
#cmd = 0.01
#cmd = -0.03   # acceleration
#cmd = -0.02   # acceleration
#cmd = -0.05   # acceleration
#cmd = -0.01*0 # acceleration
#cmd = -0.01*0 # acceleration
#dt   =  0.05
#dt   =  0.08
#dt   =  0.10
#dt   =  0.12
#dt   =  0.15
#dt   =  0.20

fig, ax1 = plt.subplots()

#dts     = [0.05,0.06,0.07,0.08,0.09,0.10,0.11,0.12,0.13,0.14,0.15]
#dts     = [0.06,0.07,0.08,0.09,0.10,0.11,0.12,0.13,0.14,0.15]
#dts     = [0.07,0.08,0.09,0.10,0.11,0.12,0.13,0.14,0.15]
#dts     = [0.08,0.09,0.10,0.11,0.12,0.13,0.14,0.15]
#dts = np.arange(0.05,0.15,0.005)

dts = np.arange(0.10,0.15,0.005)

'''
xs=dts
nstepMax = 10000
cdB  =  0.5
cdA  =  0.1
cmd  = -0.01 # acceleration
nmax = 0
#for cdB in [0.5,0.1,0.02]:
#for cdA in [0.5,0.1,0.02]:
#for cmd in [0.000,-0.005,-0.010,-0.015,-0.020]:
for cos in [0.3,0.4,0.5,0.6,0.7]:
    res = make_scan(  makeParams(), "dt", xs, bPrint=True )
    #label="cdB=%6.3f cdA=%6.3f" %(cdB,cdA)
    #label="cmd=%6.3f" %(cmd)
    label="cos=%6.3f" %(cos)
    ax1.plot( xs, res, 'o-', label=label )
    nmax = max( nmax, np.nanmax(res) )
'''

'''
nstepMax = 10000
dt   = 0.10
cdB  =  0.5
cdA  =  0.1
cmd  = -0.01
xs = [0.3,0.4,0.5,0.6,0.7]
res = make_scan(  makeParams(), "cos", xs, bPrint=True )
ax1.plot( xs, res, 'o-' )
nmax = np.nanmax(res)
'''

#cs    = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
nmins = [10,20,30,40,50]
cmds = [0.000,-0.005,-0.010,-0.015,-0.020] ; xs = cmds
#cdBs = [0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.05,0.02]
#cdBs = [0.3,0.4,0.5,0.6]; xs = cdBs
nmin = 20
cos  = 0.5 
dt   = 0.1 
nstepMax = 10000
cdB  =  0.5
cdA  =  0.3
cmd  = -0.01 # acceleration
nmax = 0
#for cdB in [0.5,0.1,0.02]:
#for cdA in [0.5,0.1,0.02]:
#for cmd in [0.000,-0.005,-0.010,-0.015,-0.020]:
#for nmin in nmins:
for cos in [0.3,0.4,0.5,0.6,0.7]:
    #for cdA in [0.1,0.2,0.3,0.4,0.5]:
    res = make_scan( makeParams(), "cmd", xs, bPrint=True )
    #label="nmin=%6.3f" %(nmin)
    #label="cdA=%6.3f" %(cdA)
    #label="nmin=%6.3f" %(nmin)
    label="cos=%6.3f" %(cos)
    ax1.plot( xs, res, 'o-', label=label )
    nmax = max( nmax, np.nanmax(res) )

# #cdA  =  0.3
# cdB  =  0.1
# res2 = make_scan(  makeParams(), "dt", xs, bPrint=True )
# ax1.plot( xs, res2, 'o-', label="cdB=%6.3f cdA=%6.3f" %(cdB,cdA)  );

plt.legend()
plt.ylim(0,nmax*1.1)
#plt.ylim(0,nstepMax)
#ax1.axhline( y=nstepMax, color='r', linestyle='--' )
plt.savefig("opt_scan.png", bbox_inches='tight')
plt.show()

#plt.savefig("opt_"+fname+spar+".png", bbox_inches='tight')
#plt.show()


