#!/usr/bin/python

import sys
sys.path.append("../../")

from pyBall.DFT import oclfft as ocl
from pyBall.DFT import jobs
from pyBall import atomicUtils as au

#jobs.Test_projectDens( )

#xyzs,Zs,enames,qs = au.loadAtomsNP( "answer.xyz")
xyzs,Zs,enames,qs = au.loadAtomsNP( "input.bas")
#jobs.Test_projectDens( atomType=Zs, atomPos=xyzs )
jobs.projectDens0_new()