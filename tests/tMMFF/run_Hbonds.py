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


ref_path = "/home/prokop/Desktop/CARBSIS/PEOPLE/Mithun/HBond_Fit_Radial_Scans/energy_angle_0/"

ref_O_NH = np.genfromtxt(ref_path+'ref-CH2O_NH3-000.dat')
ref_N_NH = np.genfromtxt(ref_path+'ref-NH3_NH3-000.dat')

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

rs = np.linspace(1.5,20.0, 100 )
nb = mmff.findHbonds( Rcut=4.0, angMax=30.0 );
for ib in range( nb ):
    c = clrs[ib]
    Es,Fs,s = mmff.sampleHbond( ib, rs                       );    plt.plot( rs, Es - Es[-1], c=c, ls='-',  label="Hb #"+str(ib)         )
    #Es,Fs,s = mmff.sampleHbond( ib, rs,           maskH=0.0  );    plt.plot( rs, Es - Es[-1], c=c, ls='--', label="Hb #"+str(ib)+" noHb" )
    #Es,Fs,s = mmff.sampleHbond( ib, rs, maskQ=0.0            );    plt.plot( rs, Es - Es[-1], c=c, ls='-',  label="Hb #"+str(ib)+" noQ |"+s  )
    #Es,Fs,s = mmff.sampleHbond( ib, rs, maskQ=0.0, maskH=0.0 );    plt.plot( rs, Es - Es[-1], c=c, ls=':',  label="Hb #"+str(ib)+" noQ |"+s  )
    #print ,s

plt.plot(ref_O_NH[:,0], ref_O_NH[:,1], label='ref_O_NH' )
plt.plot(ref_N_NH[:,0], ref_N_NH[:,1], label='ref_N_NH' )

plt.legend()
plt.grid()
plt.ylim(-0.3,0.3)
#plt.ylim(-0.1,0.1)
#plt.ylim(-0.03,0.03)

plt.show()


print("ALL DONE")
#plt.show()