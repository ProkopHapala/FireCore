#!/python3
import os
import sys
import numpy as np
#import matplotlib.pyplot as plt

deg2rad = np.pi/180.0

sys.path.append("../../")
import pyBall.atomicUtils as au

#fname = sys.argv[1]

#npbc = 16
npbc = 8
path="../../cpp/common_resources/xyz/"

mol_name="PTCDA.xyz"
mol = au.AtomicSystem( fname=path+mol_name )
# cog = mol.apos.sum(axis=0)/len(mol.apos); #print(cog)
# mol.apos-=cog
# mol.orientPCA()


surf_name="NaCl_1x1_L3.xyz"
surf = au.AtomicSystem( fname=path+surf_name )
s8x8 = surf.clonePBC(nPBC=(npbc,npbc,1) )
s8x8.saveXYZ('NaCl_8x8_L3.xyz')

mol1 = mol.clonePBC()
mol1.shift( (4.0*2.5, 4.0*2.5, 0.0) )
s8x8.append_atoms( mol1 )

s8x8.saveXYZ('NaCl_8x8_L3_PTCDA.xyz')

mol2 = mol.clonePBC()
mol2.rotate_ax( 90.0*deg2rad );
mol2.shift( (4.0*5.5, 4.0*2.5, 0.0) )
s8x8.append_atoms( mol2 )

s8x8.saveXYZ('NaCl_8x8_L3_2PTCDA.xyz')

mol3 = mol.clonePBC()
mol3.rotate_ax( 90.0*deg2rad );
mol3.shift( (4.0*2.5, 4.0*5.5, 0.0) )
s8x8.append_atoms( mol3 )

mol4 = mol.clonePBC()
mol4.shift( (4.0*5.5, 4.0*5.5, 0.0) )
s8x8.append_atoms( mol4 )

s8x8.saveXYZ('NaCl_8x8_L3_4PTCDA.xyz')

