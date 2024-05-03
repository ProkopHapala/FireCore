import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time

sys.path.append("../../")
from pyBall import FitREQ as fit
from pyBall import atomicUtils as au

# ============== Setup

Hartree2eV   = 27.2114079527

# --------- Load Reference energy
fname = "scan_H2O_b3lyp_cc-pvdz.xyz"

r0 = 1.91
Eref,xs=fit.EnergyFromXYZ(fname)    # load reference energy and position from xyz file
xs+=r0
Eref*=Hartree2eV

# ------ Set strength of hydrogen bond correction
#Hfactor = 0.99
#Hfactor = 0.9
Hfactor = 0.95
#Hfactor = 0.99
#Hfactor = 1.0

typeMask = np.array([ [0,0,0,0], [0,1,1,1], ], dtype=np.int32 )
typREQs  = np.array([ 
    #     R0vdW,    E0vdW,            Q      HBcorrection           
    # ---- imodel=1
    #[    1.487 ,    np.sqrt(0.0006808),     +0.35,         +0.22                 ],     # H
    #[    1.661 ,    np.sqrt(0.0091063),     -0.7 ,         -0.22                 ],     # O
    # ---- imodel=2
    [     1.487 ,    np.sqrt(0.0006808),     +0.35,    +np.sqrt(0.0006808)*Hfactor ],    # H
    [     1.661 ,    np.sqrt(0.0091063),     -0.7 ,    -np.sqrt(0.0091063)*Hfactor ],    # O
])   
fit.init_types( typeMask, typREQs, bCopy=True ) 
#fit.loadXYZ( fname, [0,1,2], [3,4,5] )
fit.loadXYZ( fname, [3,4,5], [0,1,2], types0=[0,1,0], testtypes=[0,1,0]  )     # load reference geometry


# ------ obtain energy profile from classical model (fit)
#Es     = fit.getEs( imodel=2, bRigid=False)
Es     = fit.getEs( imodel=2 )

# ------ obtain energy progile with HBond correction set to zero
typREQs[0,3] = 0.0              # set HBond correction to zero
fit.setType(0, typREQs[0,:] )   # set atom type 0
fit.setType(1, typREQs[1,:] )   # set atom type 1
#Es_noH = fit.getEs( imodel=1, bRigid=False) # get energy profile
Es_noH = fit.getEs( imodel=1 ) # get energy profile

print( Es )

# -------- Plot graphs
Emin = Eref.min()
plt.plot(xs,Eref  , '-',  label="E_ref", lw=3, c='grey' )
plt.plot(xs,Es_noH, '-r', label="E_noH" )
plt.plot(xs,Es    , '-g', label="E_fit" )
plt.ylim( Emin*1.2, -Emin )

plt.xlabel("r(O-H) [A]"); plt.ylabel("E[eV]");
plt.legend(); plt.grid()
plt.savefig("Hbond_correction.png", bbox_inches='tight')
plt.show()