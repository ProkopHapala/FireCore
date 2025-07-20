import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import MMFF as mmff



mmff.setVerbosity( verbosity=1, idebug=1 )


#mmff.loadXYZ("data/xyz/H2O.xyz")  # Replace with your molecule file

#mmff.init( xyz_name="data/xyz/H2O", bEpairs=True )  

mmff.init( xyz_name="data/xyz/CH2O.xyz", bEpairs=True )  



mmff.getBuffs()
mmff.setSwitches( NonBonded=-1, MMFF=1, SurfAtoms=0, GridFF=0 )

# Get initial positions
natoms    = mmff.natoms


# Select atoms for Hessian calculation (e.g., first 5 atoms)
inds = np.arange(5, dtype=np.int32)

# Calculate Hessians
h3x3 = mmff.getHessian3x3( [3,4] )
#H   = mmff.getHessian3Nx3N_selected(inds, positions)

print("3x3 Hessian blocks shape:", h3x3.shape)
for i in range(h3x3.shape[0]): print(f"Block atom {i}\n", h3x3[i])
#print("Full Hessian shape:", H.shape)