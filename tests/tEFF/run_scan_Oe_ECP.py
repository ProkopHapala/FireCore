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
    ['full'  ,'data/Oe_full.fgo',"-"          , 0.5,0 ],
    ['ecp_q8', 'data/Oe_ecp_coreCoul.fgo' ,":" , 2.0, 1 ],    # (1) cQ=-8.0e =>  coreCoul=True
    ['ecp_q6', 'data/Oe_ecp.fgo' ,"--"         , 1.5,-1 ],    # (2) cQ=-6.0e =>  coreCoul=False
]

# Turn off stdout buffering for immediate output
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', buffering=1)
eff.setVerbosity(0,0)

#eff.setSwitches(coreCoul=-1)   # (1) for Oe_ecp_coreCoul.fgo se to 1   (2) or Oe_ecp.fgo set to 0

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
#xs = np.linspace(5.0, 0.0, nconf)
xs = np.array([0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,
               0.6,0.7,0.8,0.9,1.0,
               1.2,1.4,1.6,1.8,2.0,
               2.5,3.0,4.5,5.0,6.0
               ])
energies = {}
plt.figure()
for label, fgo, ls, lw, coreCoul in jobs:
    print( " ======= ", label )
    eff.load_fgo(fgo, bVel_=True)
    print("loaded fgo")
    print(coreCoul)
    eff.setSwitches(coreCoul=coreCoul)
    print("set switches")
    eff.info()
    print("info")
    fixed_inds, fixed_poss = build_constraints(xs)
    #print("fixed_poss ", fixed_poss)
    #print("fixed_inds ", fixed_inds)
    print("Built constraints")
    eff.initOpt(dt=0.02, damping=0.001, f_limit=1000.0)
    apos, epos_scan, Es = eff.relaxed_scan(fixed_poss, fixed_inds, nstepMax=0, dt=0.02, Fconv=1e-6, ialg=2 )
    #energies[label] = Es[:,0]
    #plt.plot(xs, Es[:,0], '.-', label=label)
    #plt.plot(xs, Es[:,1], '.-', label=label)
    print("scanned")
    for icol, name, clr in columns:
        plt.plot(xs, Es[:,icol]-Es[-1,icol], ls=ls, c=clr, lw=lw, label=label+"_"+name)
    #plt.plot(xs, Es[:,iEee]+Es[:,iEae], ls=ls, c='#FF00FF', label=label+"_Ecoul")
    plt.plot(xs, Es[:,iEeePaul]+Es[:,iEaePaul], ls=ls, c='#FF00FF', lw=lw, label=label+"_EPauli")
    plt.plot(xs, Es[:,iEee    ]+Es[:,iEae    ] - Es[-1,iEee    ]-Es[-1,iEae    ], ls=ls, c='#00FFFF', lw=lw, label=label+"_Ecoul")
    #print("energies (python)\n", Es)

plt.xlabel('Distance [A]')
plt.ylabel('Total Energy [eV]')
plt.legend( loc='upper right' )
plt.grid()
plt.savefig("scan_Oe_ECP.png")
plt.show()
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
    ['full'  ,'data/Oe_full.fgo',"-"          , 0.5,-1 ],
    ['ecp_q8', 'data/Oe_ecp_coreCoul.fgo' ,":" , 2.0, 1 ],    # (1) cQ=-8.0e =>  coreCoul=True
    ['ecp_q6', 'data/Oe_ecp.fgo' ,"--"         , 1.5,-1 ],    # (2) cQ=-6.0e =>  coreCoul=False
]

# Turn off stdout buffering for immediate output
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', buffering=1)
eff.setVerbosity(0,0)

#eff.setSwitches(coreCoul=-1)   # (1) for Oe_ecp_coreCoul.fgo se to 1   (2) or Oe_ecp.fgo set to 0

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
#xs = np.linspace(5.0, 0.0, nconf)
xs = np.array([0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,
               0.6,0.7,0.8,0.9,1.0,
               1.2,1.4,1.6,1.8,2.0,
               2.5,3.0,4.5,5.0,6.0
               ])
energies = {}
plt.figure()
for label, fgo, ls, lw, coreCoul in jobs:
    print( " ======= ", label )
    eff.load_fgo(fgo, bVel_=True)
    eff.setSwitches(coreCoul=coreCoul)
    eff.info()
    fixed_inds, fixed_poss = build_constraints(xs)
    #print("fixed_poss ", fixed_poss)
    #print("fixed_inds ", fixed_inds)
    eff.initOpt(dt=0.02, damping=0.001, f_limit=1000.0)
    apos, epos_scan, Es = eff.relaxed_scan(fixed_poss, fixed_inds, nstepMax=0, dt=0.02, Fconv=1e-6, ialg=2 )
    #energies[label] = Es[:,0]
    #plt.plot(xs, Es[:,0], '.-', label=label)
    #plt.plot(xs, Es[:,1], '.-', label=label)
    for icol, name, clr in columns:
        plt.plot(xs, Es[:,icol]-Es[-1,icol], ls=ls, c=clr, lw=lw, label=label+"_"+name)
    #plt.plot(xs, Es[:,iEee]+Es[:,iEae], ls=ls, c='#FF00FF', label=label+"_Ecoul")
    plt.plot(xs, Es[:,iEeePaul]+Es[:,iEaePaul], ls=ls, c='#FF00FF', lw=lw, label=label+"_EPauli")
    plt.plot(xs, Es[:,iEee    ]+Es[:,iEae    ] - Es[-1,iEee    ]-Es[-1,iEae    ], ls=ls, c='#00FFFF', lw=lw, label=label+"_Ecoul")
    #print("energies (python)\n", Es)

plt.xlabel('Distance [A]')
plt.ylabel('Total Energy [eV]')
plt.legend( loc='upper right' )
plt.grid()
plt.savefig("scan_Oe_ECP.png")
plt.show()
