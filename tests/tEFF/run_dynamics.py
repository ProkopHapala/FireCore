import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time

sys.path.append("../../")
from pyBall import eFF as eff

def run_H2O_vs_ebullet( ie0 = -1, nsamp=100, nsubStep=10, bBsize=False ):
    eff.load_fgo("data/H2O_shoot.fgo", bVel_=True )    # load H2O moleule with electron-bullet in .fgo format (i.e. floating-gaussian-orbital) including initial velocities
    eff.getBuffs( )                                    # share internal data arrays form eFF_lib.so as numpy arrays
    eff.initOpt(dt=0.001,damping=0.0, bMass=True)      # initialize optimizer/propagator dt=0.001 [fs]
    if not bBsize: eff.invSmass[:] = 0                 # fix size of electrons by setting invMass=0 (mass=innfinity)  
    #print( "aPars\n", eff.aPars  )                    # print atom parameters (parameters of nuclei potential, charge, core electron pauli repulsion and radius)
    eff.setTrjName("trj_bullet_vs_H2O.xyz")            # setup output .xyz file to save trajectory of all atoms and electrons at each timespep (comment-out to ommit .xyz and improve performance ) 
    print( 1/eff.invMasses )                           # print masses of electrons (for debugging)
    trj_e0 = np.zeros( (nsamp,3) )                     # array to store electron bullet trajectory
    for i in range(nsamp):
        eff.run( nsubStep, ialg=-1 )                   # run eFF for 10 iterations, using leap-frog progrgator (ialg=-1) 
        trj_e0[i,:] = eff.epos[ie0,:]                  # store current position of electron bullet
        print( "E[%i] %g [eV]" %(i,eff.Es[0]) )        # Print total energy

    # plot trajectroy of electron bullet
    plt.plot( trj_e0[:,1], trj_e0[:,2], '.-' )    
    plt.axis('equal'); plt.grid()
    plt.show()

#eff.setVerbosity(1)
eff.setVerbosity(0)
#eff.eval_mol( "H2O_shoot" )
run_H2O_vs_ebullet()

