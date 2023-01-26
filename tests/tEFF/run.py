import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time

sys.path.append("../../")
from pyBall import eFF as eff



eff.test_Hatom()

exit()





def relax_mol(name):
    eff.load_fgo("data/"+name+".fgo" )                 # load molecule in  .fgo format (i.e. floating-gaussian-orbital)
    eff.setVerbosity(verbosity=1, idebug=0)             # set verbosity of simulation (defualt verbosity=0)
    eff.initOpt(dt=0.03,damping=0.1 )                   # initialize optimizer/propagator
    #eff.setPauliModel(1)                               # Pauli Repulsion model from eFF paper (http://aip.scitation.org/doi/10.1063/1.3272671) ... (Default)
    #eff.setPauliModel(0);  eff.setKPauli(30.0);        # Pauli Repulsion model proportional to overlap    EPauli = KPauli*<psi|psi>^2 for same spin electrons (=0 for oposite sipns)
    #eff.setPauliModel(2)                               # Pauli Repulsion model using only Valence-Bond theory (Eq.2 in  http://doi.wiley.com/10.1002/jcc.21637)
    eff.setTrjName(name+"_relax.xyz", savePerNsteps=1 ) # setup output .xyz file to save trajectory of all atoms and electrons at each timespep (comment-out to ommit .xyz and improve performance ) 
    eff.run( 10000, Fconv=1e-3, ialg=2 )                # run simuation for maximum 1000 time steps intil it converge to |F|<1e-3, ialg=2 is FIRE http://users.jyu.fi/~pekkosk/resources/pdf/FIRE.pdf   https://www.sciencedirect.com/science/article/pii/S0927025620300756
    eff.save_fgo(name+"_relaxed.fgo"  )                  # save final relaxed geometry to .fgo format (i.e. floating-gaussian-orbital).

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

#relax_mol("H_eFF")
#relax_mol("H2_eFF")
#relax_mol("H2O")
run_H2O_vs_ebullet()