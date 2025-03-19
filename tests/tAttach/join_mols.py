
import sys
import numpy as np
sys.path.append('../../')
from pyBall.AtomicSystem import AtomicSystem
import pyBall.plotUtils as plu
import matplotlib.pyplot as plt

#A = AtomicSystem(fname='backbone_subs_1.mol2')
#B = AtomicSystem(fname='backbone_subs_1.mol2')
#A.addSystems(B, pos=(0,5.0,0), added_bonds=[(1,4)], _0=1 )

A = AtomicSystem(fname='sequence_PPPPPPP.mol2')
#B = AtomicSystem(fname='porph_Cytosine/porph2s_C.mol2')
#A.addSystems(B, pos=(25.0,10.0,0), rot=[(-1,0,0),(0,-1,0),(0,0,1)] )

B = AtomicSystem(fname='porphs_Guanine/porph2s_G.mol2')
A.addSystems(B, pos=(-23.0,9.0,0) )
A.addSystems(B, pos=(-23.0,20.0,0) )

plu.plotSystem( A, bLabels=False )
plt.savefig('join_mols.png')
plt.grid()
plt.show()

print( "A.addSystems(B) save files: A_B.mol2, A_B.xyz" )
A.saveXYZ  ('A_B.xyz')
A.save_mol2('A_B.mol2')

