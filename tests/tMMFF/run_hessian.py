import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import MMFF as mmff
from pyBall import grid_utils as gu



mmff.setVerbosity( verbosity=1, idebug=1 )


#mmff.loadXYZ("data/xyz/H2O.xyz")  # Replace with your molecule file

#mmff.init( xyz_name="data/xyz/H2O", bEpairs=True )  ; selected_atoms = [3,4]

mmff.init( xyz_name="data/xyz/CH2O.xyz", bEpairs=True )  ; selected_atoms = [4,5]



mmff.getBuffs()
mmff.setSwitches( NonBonded=-1, MMFF=1, SurfAtoms=0, GridFF=0 )

# Get initial positions
natoms    = mmff.natoms

extent=[-1.0,1.0, -1.0,1.0]


X,Y = gu.makeGrid2D( extent=extent, nxy=[100,100], endpoint=False, dpix=0.5 )
ps = np.zeros( (X.shape[0]*X.shape[1],3) )
ps[:,0] = X.flatten()
ps[:,1] = Y.flatten()
ps[:,2] = 0.0

for ia in selected_atoms: 
    E,F = mmff.scan_atoms_rigid( [ia], ps, bRelative=True )
    plt.imshow(E.reshape(X.shape), origin='lower', cmap='viridis')
    plt.colorbar()
    plt.title(f'Atom {ia}')
    plt.show()

# Select atoms for Hessian calculation (e.g., first 5 atoms)
inds = np.arange(5, dtype=np.int32)

# Calculate Hessians
h3x3 = mmff.getHessian3x3( selected_atoms )
#H   = mmff.getHessian3Nx3N_selected(inds, positions)

print("3x3 Hessian blocks shape:", h3x3.shape)
for i in range(h3x3.shape[0]): print(f"Block atom {i}\n", h3x3[i])
#print("Full Hessian shape:", H.shape)