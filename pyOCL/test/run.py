#!/usr/bin/python

import sys
sys.path.append("../")

from pyOCL import oclfft as ocl
from pyOCL import jobs
from pyOCL import atomicUtils as au


xyzs,Zs,enames,qs = au.loadAtomsNP( "answer.xyz")

print( xyzs )

#jobs.Test_projectDens( )