#!/usr/bin/python

import sys
sys.path.append("../")

from pyOCL import oclfft as ocl
from pyOCL import jobs
from pyOCL import atomicUtils as au

#jobs.Test_projectDens( )

xyzs,Zs,enames,qs = au.loadAtomsNP( "answer.xyz")
jobs.Test_projectDens( atomType=Zs, atomPos=xyzs )