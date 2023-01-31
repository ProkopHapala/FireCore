import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time

sys.path.append("../../")
from pyBall import eFF as eff

from pyBall import eFF_terms as pyeff

'''
rho=-0.2; r=1.2; si=0.3; sj=0.8
rho=-0.2; r=0.7; si=1.3; sj=0.15

Euuref,DTref, Sref = pyeff.pyeff_E_up_up   ( r, si, sj, rho,)
Eudref             = pyeff.pyeff_E_up_down ( r, si, sj, rho,)
Euu,Eud, DT, S     = pyeff.pyeff_EPaul     ( r, si, sj, rho,)

print( "Eupdown ", Eudref,  Eud  )
print( "Eupup   ", Euuref,  Euu  )
print( "S       ", Sref,  S  )
print( "DT      ", DTref, DT )

eff.eval_ee( r, si, sj )
'''

#eff.check_DerivsPauli( r0=0.0,r1=2.5,   s0=0.5,s1=0.5, n=100 , spin=1 )
#eff.check_DerivsPauli( r0=0.0,r1=2.5,   s0=0.5,s1=0.5, n=100 , spin=-1 )
#eff.check_DerivsPauli( r0=0.7,r1=0.7,   s0=0.2,s1=2.5,  sj=0.5, n=100 , spin=1)
#eff.check_DerivsPauli( r0=0.7,r1=0.7,   s0=0.2,s1=2.5,  sj=0.5, n=100 , spin=-1)
#exit()

#eff.setSwitches( pauli=1,   kinetic=-1, coulomb=-1, AA=-1, AE=-1, AECoulomb=-1, AEPauli=-1  )
#eff.setSwitches( pauli=-1,   kinetic=1, coulomb=-1, AA=-1, AE=-1, AECoulomb=-1, AEPauli=-1  )
#eff.setSwitches( pauli=-1,   kinetic=-1, coulomb=1, AA=1, AE=1, AECoulomb=1, AEPauli=1  )

#eff.check_Derivs_ie( "e2_1g_2o_triplet", ie=0, r0=0.5,r1=1.5, s0=0.5,s1=0.5, n=100 )
#eff.check_Derivs_ie( "H2_eFF_asym", ie=0, r0=0.5,r1=1.5, s0=0.5,s1=0.5, n=100 )
#eff.check_Derivs_ie( "H2_eFF_asym", ie=0, r0=1.0,r1=1.0, s0=0.5,s1=1.5, n=100 )

def run_H2O_vs_ebullet( ie0 = -1, nsamp=100, bBsize=False ):
    eff.load_fgo("data/H2O_shoot.fgo", bVel_=True )    # load H2O moleule with electron-bullet in .fgo format (i.e. floating-gaussian-orbital) including initial velocities
    eff.getBuffs( )                                    # share internal data arrays form eFF_lib.so as numpy arrays
    eff.initOpt(dt=0.001,damping=0.0, bMass=True)      # initialize optimizer/propagator
    if not bBsize: eff.invSmass[:] = 0                 # fix size of electrons by setting invMass=0 (mass=innfinity)  
    #print( "aPars\n", eff.aPars  )                    # print atom parameters (parameters of nuclei potential, charge, core electron pauli repulsion and radius)
    eff.setTrjName("trj_bullet_vs_H2O.xyz")            # setup output .xyz file to save trajectory of all atoms and electrons at each timespep (comment-out to ommit .xyz and improve performance ) 
    print( eff.invMasses )
    trj_e0 = np.zeros( (nsamp,3) )                     # array to store trajectroy electron bullet 
    for i in range(100):
        eff.run( 10, ialg=-1 )                         # run eFF for 10 iterations, using leap-frog progrgator (ialg=-1) 
        trj_e0[i,:] = eff.epos[ie0,:]                  # store current position of electron bullet
        #print( "E[%i] %g [eV]" %(i,eff.Es[0]) )        # Print total energy

    # plot trajectroy of electron bullet
    plt.plot( trj_e0[:,1], trj_e0[:,2], '.-' )    
    plt.axis('equal'); plt.grid()
    plt.show()

eff.setVerbosity(1)

const_bohr = 0.5291772105638411


#eff.test_Hatom()        ;exit()
#eff.test_Hatom(True)    ;exit()
#eff.check_H2(False)     ;exit()
#eff.check_H2()          ;exit()

# limit Euu(r=0,si=sj) = 1/s^2 = 1/1.7007535132532356^2 = 0.34571422244 [Hartree] = 9.407362451147032 [eV]
#   !!! But this should be infinite ( two same-spin electrons cannot occupy same orbital !!!


#Euu,_,_,_    = pyeff.E_up_up__( 1.125e-7, 1.7007535132532354, 1.7007535132532354, rho=-0.2 )  
#Euu,DT,S    = pyeff.pyeff_E_up_up(  1.125e-7, 1.7007535132532354, 1.7007535132532354, rho=-0.2 ); print( "!!!!! pyeff_E_up_up(): Euu",Euu, "DT",DT, "S",S )
#Euu,DT,S     = pyeff.pyeff_E_up_up(  1.125e-7, 0.9/const_bohr, 0.9/const_bohr, rho=-0.2 );     print( "!!!!! pyeff_E_up_up(): Euu",Euu, "DT",DT, "S",S )
#Euu,Eud,DT,S = pyeff.pyeff_EPaul( 1e-5, 0.9/const_bohr, 0.9/const_bohr, rho=-0.2 );            print( "!!!!! pyeff_EPaul():   Euu",Euu, "Eud",Eud, "DT",DT, "S",S )

#eff.check_H2( False, "H2_eFF_upup")     ;exit()
#eff.check_H2( False, "H2_eFF_asym")     ;exit()
#eff.check_H2( False, "H2_eFF_sym")     ;exit()
#eff.check_H2( True, "H2_eFF")     ;exit()

#eff.relax_mol("H_eFF")
#eff.relax_mol("H2_eFF")
#eff.relax_mol("H2O")
#eff.relax_mol("CH4")

#eff.setVerbosity(3)
eff.eval_mol("CH4_lmps", fUnits=const_bohr)
#eff.eval_mol("CH4_lmps_fixcore", fUnits=const_bohr)

exit()

nmax=2000
outE=np.zeros(nmax)
eff.relax_mol("CH4_lmps", dt=0.005, nMaxIter=nmax, outE=outE )

plt.plot(outE); plt.show()

#run_H2O_vs_ebullet()