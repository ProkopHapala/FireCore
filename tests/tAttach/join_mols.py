
import sys
import numpy as np
sys.path.append('../../')
from pyBall.AtomicSystem import AtomicSystem

A = AtomicSystem(fname='backbone_subs_1.mol2')
B = AtomicSystem(fname='backbone_subs_1.mol2')

A.addSystems(B, pos=(0,5.0,0), added_bonds=[(1,4)], _0=1 )

print( "A.addSystems(B) save files: A_B.mol2, A_B.xyz" )
A.saveXYZ  ('A_B.xyz')
A.save_mol2('A_B.mol2')

