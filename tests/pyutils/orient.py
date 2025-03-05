#!/usr/bin/python
import sys
sys.path.append("../../")
#import pyBall.atomicUtils as au
from pyBall.AtomicSystem import AtomicSystem
fname = sys.argv[1]
mol = AtomicSystem( fname=fname )
#mol.orient( 2, (5,2), (5,6), trans=(2,1,0)  )
cog = mol.apos.sum(axis=0)/len(mol.apos); #print(cog)
mol.apos-=cog
mol.orientPCA()
mol.saveXYZ(fname+'-oriented.xyz')
