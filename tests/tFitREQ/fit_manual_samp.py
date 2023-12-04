import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time

sys.path.append("../../")
from pyBall import FitREQ as fit
from pyBall import atomicUtils as au

# ============== Setup

fname = "samp"
imodel = 3

# ============== Functions

fit.setVerbosity(1)
#fit.loadXYZ_new( fname+".xyz", fname_AtomTypes="atypes.dat" )     # load reference geometry

#fit.loadXYZ_new( fname+".xyz", fname_AtomTypes="atypes.dat", bAddEpairs=True, bOutXYZ=True )     # load reference geometry
fit.loadXYZ_new( fname+".xyz", bAddEpairs=True, bOutXYZ=True )     # load reference geometry

#fit.getBuffs()
#print( "typToREQ\n" , fit.typToREQ )
#print( "typeREQs\n" , fit.typeREQs )
#print( "types1\n"   , fit.types1   )
#print( "types2"   , fit.types2   )
#print( "types3"   , fit.types3   )

# ------ obtain energy profile from classical model (fit)
Es     = fit.getEs( imodel=imodel, isampmode=2 )     #;print( "Es_noH", Es     )
print( "Es", Es )



'''
# ------ obtain energy progile with HBond correction set to zero
#typREQs[0,3] = 0.0              # set HBond correction to zero
#fit.setType(0, typREQs[0,:] )   # set atom type 0
#Es_noH = fit.getEs( imodel=imodel, bRigid=False)     #;print( "Es_noH", Es_noH )

print( "rs.shape", rs.shape, "Es.shape", Es.shape, "Eref.shape", Eref.shape ) 
plt.plot( rs, Eref  , '-k', label="E_ref", lw=3 )
plt.plot( rs, Es    , '-g', label="E",     lw=2 )
#plt.plot( rs, Es_noH, '-r', label="E_noH", lw=2 )
plt.axhline(0.0, c='k', ls='--')
plt.grid()
plt.legend()
plt.ylim( -0.03, 0.03 )
#plt.ylim( -2.0, 2.0 )
'''

plt.show()
