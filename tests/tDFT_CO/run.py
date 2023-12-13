#!/usr/bin/python


# INTEROPERABLE GLOBAL VARIABLES
# https://gcc.gnu.org/onlinedocs/gfortran/Interoperable-Global-Variables.html#Interoperable-Global-Variables


import sys
import numpy as np

sys.path.append("../../")
from pyBall.DFT import oclfft as ocl
from pyBall.DFT import utils as oclu
from pyBall.DFT import high_level as oclh
from pyBall.DFT import jobs
from pyBall     import atomicUtils as au


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

def export_density_shifted():
    shift = [ngrid[0]//2,ngrid[1]//2,ngrid[2]//2]
    iA=0; iB=1
    print( "# ---- total SCF density" )
    jobs.projectDens( iOutBuff=iA, atomType=Zs, atomPos=apos, iMO0=0, ngrid=ngrid, dcell=dcell, bSaveXsf=False, bSaveBin=False, bSCF=True, bDen0diff=False )
    ocl.roll( iA, iB, shift )
    ocl.saveToXsf( "dens_scf.xsf", iB )
    ocl.saveToBin( "dens_scf.bin", iB )
    
    exit(0)

    print( "# ---- density difference" )
    jobs.projectDens0_new( iOutBuff=iA,  atomType=Zs, atomPos=apos, ngrid=ngrid, dcell=dcell, bSaveXsf=False, bSaveBin=False, acumCoef=[1.0,-1.0] )
    ocl.roll( iA, iB, shift )
    ocl.saveToXsf( "dens_diff.xsf", iB )
    ocl.saveToBin( "dens_diff.bin", iB )


#ocl.tryInitFFT( ngrid)           ;print( "DEBUG poisson 1 " )

#jobs.projectDens ( iOutBuff=iA, atomType=Zs, atomPos=apos, iMO0=1, iMO1=8, ngrid=ngrid, dcell=dcell, bSaveXsf=True, bSaveBin=True )
#jobs.projectDens ( iOutBuff=iA, atomType=Zs, atomPos=apos, iMO0=1, iMO1=8, ngrid=ngrid, dcell=dcell, bSaveXsf=True, bSaveBin=True, saveName="dens_scf"  )
#jobs.projectDens0( iOutBuff=iA, atomType=Zs, atomPos=apos,                 ngrid=ngrid, dcell=dcell, bSaveXsf=True, bSaveBin=True, saveName="dens_diff"  )


export_density_shifted()



#jobs.projectDens( iOutBuff=iA, atomType=Zs, atomPos=apos, iMO0=1, iMO1=102//2, ngrid=ngrid, dcell=dcell, bSaveXsf=False, bSaveBin=True )
#ocl.initFFT( ngrid )
#ocl.tryInitFFT( ngrid)           ;print( "DEBUG poisson 1 " )
#ocl.runfft (iA )
#ocl.poisson   (  iA,iC, dcell )  ;print( "DEBUG poisson 3 " )
