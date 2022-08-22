#!/usr/bin/python


# INTEROPERABLE GLOBAL VARIABLES
# https://gcc.gnu.org/onlinedocs/gfortran/Interoperable-Global-Variables.html#Interoperable-Global-Variables


import sys
import numpy as np

sys.path.append("../")
from pyOCL import oclfft as ocl
from pyOCL import utils as oclu
from pyOCL import high_level as oclh
from pyOCL import jobs
from pyOCL import atomicUtils as au


#ocl.setErrorCheck( 0 )
ocl.setErrorCheck( 1 )

#xyzs,Zs,enames,qs = au.loadAtomsNP( "answer.xyz")
#xyzs,Zs,enames,qs = au.loadAtomsNP( "pentacene.xyz")
#xyzs,Zs,enames,qs = au.loadAtomsNP( "CH4.xyz")
#jobs.Test_projectDens( atomType=Zs, atomPos=xyzs )

Zs = np.array([ 8,6 ], np.int32)

apos = np.array([
[ 0.0,0.0, 0.1   ],
[ 0.0,0.0, 1.228 ],
])

ngrid=(128,64,32)
dcell = [0.2,0.2,0.2,0.2]
iA=0; iC=1


#ocl.tryInitFFT( ngrid)           ;print( "DEBUG poisson 1 " )

#jobs.projectDens ( iOutBuff=iA, atomType=Zs, atomPos=apos, iMO0=1, iMO1=8, ngrid=ngrid, dcell=dcell, bSaveXsf=True, bSaveBin=True )
jobs.projectDens ( iOutBuff=iA, atomType=Zs, atomPos=apos, iMO0=1, iMO1=8, ngrid=ngrid, dcell=dcell, bSaveXsf=False, bSaveBin=True )
jobs.projectDens0( iOutBuff=iA, atomType=Zs, atomPos=apos, iMO0=0, iMO1=8, ngrid=ngrid, dcell=dcell, bSaveXsf=True, bSaveBin=True )


#jobs.projectDens( iOutBuff=iA, atomType=Zs, atomPos=apos, iMO0=1, iMO1=102//2, ngrid=ngrid, dcell=dcell, bSaveXsf=False, bSaveBin=True )
#ocl.initFFT( ngrid )
#ocl.tryInitFFT( ngrid)           ;print( "DEBUG poisson 1 " )
#ocl.runfft (iA )
#ocl.poisson   (  iA,iC, dcell )  ;print( "DEBUG poisson 3 " )
