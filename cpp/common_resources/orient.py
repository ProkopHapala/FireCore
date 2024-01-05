#!/usr/bin/python
import sys
sys.path.append("../../")
import pyBall.atomicUtils as au
fname = sys.argv[1]
mol = au.AtomicSystem( fname=fname )
#mol.orient( 2, (5,2), (5,6), trans=(2,1,0)  )
mol.orientPCA()
mol.saveXYZ(fname+'-oriented.xyz')
