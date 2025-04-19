import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time

sys.path.append("../../")
from pyBall import eFF as eff
#from pyBall import eFF_terms as pyeff

# Set line length to unlimited for numpy array printing
np.set_printoptions(linewidth=np.inf)


'''
Eis[0] = ff.Etot;
Eis[1] = ff.Ek;
Eis[2] = ff.Eee;
Eis[3] = ff.EeePaul;
Eis[4] = ff.EeeExch;
Eis[5] = ff.Eae;
Eis[6] = ff.EaePaul;
Eis[7] = ff.Eaa;
'''

def build_constraints(xs, fixed_inds=None, fixed_poss=None, vidx=-1 ):
    eff.getBuffs()                                          # refresh buffers after loading FGO
    if fixed_inds is None: fixed_inds = np.array([[vidx, 0b1111]], dtype=np.int32)  # constrain that electron's DOFs
    if fixed_poss is None: fixed_poss = np.zeros((len(xs),1,4))                    # positions for each scan point
    fixed_poss[:,0,0] = xs                                  # x-coordinate varying
    fixed_poss[:,0,1] = eff.epos[vidx,1]                    # y,z,size static
    fixed_poss[:,0,2] = eff.epos[vidx,2]
    fixed_poss[:,0,3] = eff.esize[vidx]
    return fixed_inds, fixed_poss


jobs=[
    ['full','data/Oe_full.fgo',"-"],
    ['ecp', 'data/Oe_ecp.fgo',"--"]
]

# Turn off stdout buffering for immediate output
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', buffering=1)
eff.setVerbosity(1,0)

columns=[
    (0,"Etot","k"),
    #(1,"Ek",  "r"),
    #(2,"Eee", "g"),
    #(3,"EeePaul","cyan"),
    #(4,"EeeExch","k"),
    #(5,"Eae","r"),
    #(6,"EaePaul","orange"),
    #(7,"Eaa","k")
]

iEtot=0
iEk=1 
iEee=2 
iEeePaul=3 
iEeeExch=4 
iEae=5 
iEaePaul=6 
iEaa=7 




# ------ Body
nconf = 50
xs = np.linspace(5.0, 0.0, nconf)
energies = {}
plt.figure()
for label, fgo, ls in jobs:
    print( " ======= ", label )
    eff.load_fgo(fgo, bVel_=True)
    fixed_inds, fixed_poss = build_constraints(xs)
    #print("fixed_poss ", fixed_poss)
    #print("fixed_inds ", fixed_inds)
    eff.initOpt(dt=0.02, damping=0.001, f_limit=1000.0)
    apos, epos_scan, Es = eff.relaxed_scan(fixed_poss, fixed_inds, nstepMax=0, dt=0.02, Fconv=1e-6, ialg=2 )
    #energies[label] = Es[:,0]
    #plt.plot(xs, Es[:,0], '.-', label=label)
    #plt.plot(xs, Es[:,1], '.-', label=label)
    for icol, name, clr in columns:
        plt.plot(xs, Es[:,icol], ls=ls, c=clr, label=label+"_"+name)
    plt.plot(xs, Es[:,iEee]+Es[:,iEae], ls=ls, c='#FF00FF', label=label+"_Ecoul")
    #print("energies (python)\n", Es)

plt.xlabel('Distance [A]')
plt.ylabel('Total Energy [eV]')
plt.legend( loc='upper right' )
plt.grid()
plt.savefig("scan_Oe_ECP.png")
plt.show()
