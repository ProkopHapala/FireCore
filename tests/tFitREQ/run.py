import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time

sys.path.append("../../")
from pyBall import FitREQ as fit
from pyBall import atomicUtils as au

# ============== Setup

h2o = au.AtomiSystem( "common_resources/H2O.xyz")
h2o.print()
h2o.findBonds()
h2o.printBonds()
exit()

# fit Q
#typeMask = np.array([ [0,0,0], [0,0,1], ], dtype=np.int32 )
#typREQs  = np.array([ [ 1.0,0.0, -0.2],[ 1.0,0.0, +0.2] ] )   # charge = 0.2
# fit E0
#typeMask = np.array([ [0,0,0], [0,1,0], ], dtype=np.int32 )
#typREQs  = np.array([ [ 1.0,0.1, 0.0],[ 1.0,0.1, 0.0] ] )   # charge = 0.2
# fit R
#typeMask = np.array([ [0,0,0], [1,0,0], ], dtype=np.int32 )
#typREQs  = np.array([ [ 1.0,0.1, 0.0],[ 1.2,0.1, 0.0] ] )   # charge = 0.2

# fit E0 Q R
typeMask = np.array([ [0,0,0], [0,1,1], ], dtype=np.int32 )
typREQs  = np.array([ [ 1.0,0.1,-0.2],[ 1.2,0.1, 0.2] ] )   # charge = 0.2

ne = 50
poses        = np.zeros   ( (ne,3,3) )
poses[:,0,0] = np.linspace( 0,6.0, ne )
poses[:,1,0] = 1.0
poses[:,2,1] = 1.0
Es           = np.zeros( ne )

m1_types=np.array([0], dtype=np.int32)
m1_ps   =np.array([-1.5,0.0,0.0])

m2_types=np.array([1], dtype=np.int32)
m2_ps   =np.array([0.0,0.0,0.0])

# =========== Run

#print( "poses \n", poses )

fit.init_types( typeMask, typREQs,  bCopy=True )               
fit.setRigidSamples( Es, poses,     bCopy=True, bAlloc=True )  
fit.setSystem( -1, m1_types, m1_ps, bCopy=True )               
fit.setSystem( -2, m2_types, m2_ps, bCopy=True )             
fit.setSystem( -3, m2_types, m2_ps, bCopy=True )              

fit.getBuffs()  
Es = fit.getEs(bRigid=True)     # ;print(Es)
#fit.DOFs[:] = 0
fit.DOFs[0] = 0.0
fit.DOFs[1] = 0.0
#fit.DOFs[2] = 0.0
fit.Es  [:] = Es[:]

DOFs0  = fit.DOFs.copy() 

#print(fit.poses)

#fit.run(100,1e-4,0.1,True, 0)      
fit.run(300,1e-4,0.2,True, 1)  

print( "BEFOR fit.REQ", typREQs[1]   )
print( "BEFOR fit.DOFs", DOFs0     )
print( "AFTER fit.DOFs", fit.DOFs  )
#print( "fit.typeREQs", fit.typeREQs )


Es_after = fit.getEs(bRigid=True)      #;print(Es)

xs= fit.poses[:,0,0]
plt.plot( xs, fit.Es  ,'.-',label="E_ref" )
plt.plot( xs, Es_after,'.-',label="E_fit" )
plt.show()

print("ALL DONE")
