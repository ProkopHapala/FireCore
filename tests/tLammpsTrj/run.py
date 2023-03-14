import sys
import matplotlib.pyplot as plt
sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import plotUtils   as plu

#selection=[1,2,3,4]
selection=range(100)

# ---- read trajectory
trj = au.readLammpsTrj( fname="traj.lammpstrj", nmax=100, selection=set(selection) )

# ---- multiply trajectroy in periodic boundary conditions (PBC)
trj_221 = [ sys.clonePBC( (2,2,1) ) for sys in trj  ]

# ---- save PBC-mutiplicated trajcetory
#open("trj.xyz"    ,"w").close(); 
#for i,sys in enumerate(trj): sys.saveXYZ( "trj.xyz", mode="a" )
open("trj_221.xyz","w").close(); 
for i,sys in enumerate(trj_221): sys.saveXYZ( "trj_221.xyz", mode="a" )

'''
# ---- save trajectroy to PBC
#plu.plotTrj( trj, bBonds=True, sz=50., numbers=None, axes=(0,1), extent=None, prefix="mol_" )
plu.plotTrj( trj_221, bBonds=True, sz=50., numbers=selection, axes=(0,1), prefix="mol_", extent=(-10.,10.,0.0,20.) )
'''