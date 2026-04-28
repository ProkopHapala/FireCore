import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time

sys.path.append("../../")
from pyBall import eFF as eff
#from pyBall import eFF_terms as pyeff

# ----- Fucntions

def plot_energy_and_force(energies, forces, steps=None,  lable1='Energ [eV]', lable2='Force eV/Å', figsize=(8,5) ):
    """
    Plot energy and force on shared x-axis with separate y-axes.
    
    Parameters:
    -----------
    steps : array-like
        Step numbers or simulation time points
    energies : array-like
        Total energy values
    forces : array-like
        Force magnitude values (will be log10 transformed)
    lable1 : str, optional
        Unit label for energy axis (default: 'eV')
    lable2 : str, optional
        Unit label for force axis (default: 'eV/Å')
    """
    # Create figure and primary axis
    fig, ax1 = plt.subplots(figsize=figsize)
    if steps is None: steps=range(len(energies))

    # Plot energy on primary axis (left)
    ax1.plot(steps, energies, ".-k", label='Total energy', ms=1.5, lw=0.5)
    ax1.set_xlabel('Step')
    ax1.set_ylabel( lable1, color='k')
    ax1.tick_params(axis='y', labelcolor='k')
    ax1.grid(True, alpha=0.3)
    
    # Create secondary axis (right)
    ax2 = ax1.twinx()
    
    # Plot log10 force on secondary axis
    ax2.plot(steps, np.log10(forces), ".-r", label='log10|Force|', ms=1.5, lw=0.5)
    ax2.set_ylabel( lable2, color='r')
    ax2.tick_params(axis='y', labelcolor='r')
    
    # Combine legends and adjust layout
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right')
    
    plt.title('Energy and Force Convergence')
    plt.tight_layout()
    return fig, (ax1, ax2)

# ------ Body

#eff.setVerbosity(1) # If 1: it will  write more stuff to console; If 0: it will wirte less stuff to console
#eff.setVerbosity(2)
eff.setVerbosity(0)

eff.load_fgo("data/H2_eFF.fgo" )  
#eff.load_fgo("data/H2_far.fgo", bVel_=True )   
#eff.load_fgo("data/H2O.fgo", bVel_=True )   
eff.getBuffs()

# --- setup constraints
fixed_inds = np.array([
#    index, binary mask    
(0, 0b111), # 1st atom, x,y,z
(1, 0b111), # 2nd atom, x,y,z
], dtype=np.int32)

# # ---- optimization for single configuration ( investigate best optimizer settings to find optimal dt, damping, f_limit for fast convergence)
# fixed_poss = np.array([ [0,0,0,0], [0,0,1.0,0] ], dtype=np.float64)
# eff.setTrjName("trj.xyz", savePerNsteps=1, bDel=True)
# eff.initOpt(0.1,0.001,1000.0)
# eff.set_constrains(2, fixed_poss, fixed_inds)
# outE, outF = eff.run( nstepMax=1000, dt=0.02, Fconv=1e-6, ialg=2, bOutF=True, bOutE=True )
# fig, (ax1, ax2) = plot_energy_and_force( energies=outE, forces=outF )
# plt.show()
#exit()

nconf = 27
# --- setup positions of fixed atoms (1st atoms at (0,0,0), 2nd at (x,0,0) )
xs = np.linspace( 3.0, 0.3, nconf,endpoint=False);    
#xs = np.linspace( 0.3, 3.0, nconf,endpoint=False);    

fixed_poss = np.zeros((nconf, eff.na, 4 ))
fixed_poss[:,1,0] = xs   # set x coordinate of 2nd atom
print( "xs ", xs )

#exit()

with open("scan.xyz", "w") as f: f.write( "" )  # clear file, since we are appending it inside eff.relaxed_scan()

eff.initOpt( dt=0.02, damping=0.001, f_limit=1000.0)
#apos, epos, Es = 1eff.relaxed_scan( fixed_poss, fixed_inds, nstepMax=100, dt=1e-2, Fconv=1e-6, ialg=0 )
apos, epos, Es = eff.relaxed_scan( fixed_poss, fixed_inds, nstepMax=10000, dt=0.02, Fconv=1e-6, ialg=2, scan_trj_name="scan.xyz" )


plt.figure(figsize=(5,15))
# plt.plot(radiusSpace, Etot_value,"o--", label='Etot', color='blue')
plt.subplot(311)
plt.plot( xs, Es[:,0],".-k", label='Etot')
plt.xlabel('Radius [A]')
plt.ylabel('Energy [eV]')
plt.legend()
plt.grid()

plt.subplot(312)
plt.plot( xs, Es[:,1],".:", label='Ek')
plt.plot( xs, Es[:,2],".:", label='Eee')
plt.plot( xs, Es[:,3],".:", label='EeePaul')
plt.plot( xs, Es[:,4],".:", label='EeeExch')
plt.plot( xs, Es[:,5],".:", label='Eae')
plt.plot( xs, Es[:,6],".:", label='EaePaul')
plt.plot( xs, Es[:,7],".:", label='Eaa')
plt.xlabel('Radius [A]')
plt.ylabel('Energy [eV]')
plt.legend()
plt.grid()

plt.subplot(313)
plt.plot( xs, epos[:,0,3],".-", label='e1 size')
plt.plot( xs, epos[:,1,3],".-", label='e2 size')
plt.plot( xs, epos[:,0,0],".-", label='e1 x')
plt.plot( xs, epos[:,1,0],".-", label='e2 x')
plt.plot( xs, apos[:,0,0],".--", label='a1 x')
plt.plot( xs, apos[:,1,0],".--", label='a2 x')
plt.xlabel('Radius   [A]')
plt.ylabel('Position [A]')
plt.legend()
plt.grid()
plt.savefig("scan.png")
plt.show()