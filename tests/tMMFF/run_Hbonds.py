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
#xyz_name = "polymer-2_new-OH"
xyz_name = "formic_dimer"



ref_path = "/home/prokop/Desktop/CARBSIS/PEOPLE/Mithun/HBond_Fit_Radial_Scans/energy_angle_0/"

#  Charges :
#  NH3    +0.28363  -0.84922      
#  H2O    +0.34068  -0.68135
#  HF     +0.41277  -0.41277    



#ref_O_NH   = np.genfromtxt(ref_path+'ref-CH2O_NH3-000.dat')
#ref_N_NH   = np.genfromtxt(ref_path+'ref-NH3_NH3-000.dat')
ref_H2O_H2O =  np.genfromtxt(ref_path+'ref-H2O_H2O-000.dat' )
ref_HF_HF   =  np.genfromtxt(ref_path+'ref-HF_HF-000.dat'   )
ref_NH3_NH3 =  np.genfromtxt(ref_path+'ref-NH3_NH3-000.dat' )


#ialg  = 2
#ialg  = 3
#alg_names=["GD","MDdamp","FIRE","FIREsmooth"]

#======== Body

#alg_name = alg_names[ialg]
mmff.setVerbosity( verbosity=2, idebug=0 )
mmff.init( xyz_name="data/"+xyz_name, bMMFF=True )     
#mmff.getBuffs()
#mmff.eval()
rs = np.linspace(1.5,20.0, 100 )

Es,Fs = mmff.sampleNonBondTypes( "H_F   F"  , rs, qH=+0.41277, qX=-0.41277 );    plt.plot( rs, Es - Es[-1], '-g',  label="H_F   F"    )
Es,Fs = mmff.sampleNonBondTypes( "H_OH  O_3", rs, qH=+0.34068, qX=-0.34068 );    plt.plot( rs, Es - Es[-1], '-r',  label="H_OH  O_3"  )
Es,Fs = mmff.sampleNonBondTypes( "H_NH2 N_3", rs, qH=+0.28363, qX=-0.28363 );    plt.plot( rs, Es - Es[-1], '-b',  label="H_NH2 N_3"  )

'''
clrs = [ 'r','g','b','m','c' ]

nb = mmff.findHbonds( Rcut=4.0, angMax=30.0 );
for ib in range( nb )[:1]:
    c = clrs[ib]
    Es,Fs,s = mmff.sampleHbond( ib, rs                   , kind=1            );    plt.plot( rs, Es - Es[-1], '-k',  label="Hb #"+str(ib) + s         )
    Es,Fs,s = mmff.sampleHbond( ib, rs                   , kind=2, dcomp=3.0 );    plt.plot( rs, Es - Es[-1], '-k',  label="Hb #"+str(ib) + s         )
    Es,Fs,s = mmff.sampleHbond( ib, rs                   , kind=2, dcomp=1.5 );    plt.plot( rs, Es - Es[-1], '--k', label="Hb #"+str(ib) + s         )
    #Es,Fs,s = mmff.sampleHbond( ib, rs,           maskH=0.0  );    plt.plot( rs, Es - Es[-1], c=c, ls='--k', label="Hb #"+str(ib)+" noHb" )
    Es,Fs,s = mmff.sampleHbond( ib, rs, maskQ=0.0            );    plt.plot( rs, Es - Es[-1], ':k',  label="Hb #"+str(ib)+" noQ |"+s  )
    #Es,Fs,s = mmff.sampleHbond( ib, rs, maskQ=0.0, maskH=0.0 );    plt.plot( rs, Es - Es[-1], c=c, ls=':k',  label="Hb #"+str(ib)+" noQ |"+s  )
    #print ,s
'''

#plt.plot(ref_O_NH[:,0]-0.3,    ref_O_NH[:,1],   label='ref_O_NH' )
#plt.plot(ref_N_NH[:,0]-0.3,    ref_N_NH[:,1],   label='ref_N_NH' )
plt.plot(ref_HF_HF[:,0]-0.3,   ref_HF_HF[:,1],  ".:g", label='ref HF_HF'   )
plt.plot(ref_H2O_H2O[:,0]-0.3, ref_H2O_H2O[:,1],".:r", label='ref H2O_H2O' )
plt.plot(ref_NH3_NH3[:,0]-0.3, ref_NH3_NH3[:,1],".:b", label='ref NH3_NH3' )

plt.legend()

plt.ylim(1.5,20.0);   # plt.xticks(np.arange(1,5, 20.5, 0.2)) 
plt.ylim(-0.3,0.3)
plt.minorticks_on()
plt.grid(which='both')
#plt.ylim(-0.1,0.1)
#plt.ylim(-0.03,0.03)

plt.show()


print("ALL DONE")
#plt.show()