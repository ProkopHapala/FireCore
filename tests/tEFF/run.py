import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time

sys.path.append("../../")
from pyBall import eFF as eff

from pyBall import eFF_terms as pyeff


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



exit()

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

#eff.test_Hatom()       ;exit()
#eff.test_Hatom(True)   ;exit()
eff.check_H2(False)     ;exit()
#eff.check_H2()          ;exit()

#eff.relax_mol("H_eFF")
#eff.relax_mol("H2_eFF")
#eff.relax_mol("H2O")
#run_H2O_vs_ebullet()