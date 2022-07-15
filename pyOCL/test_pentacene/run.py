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
xyzs,Zs,enames,qs = au.loadAtomsNP( "pentacene.xyz")
#xyzs,Zs,enames,qs = au.loadAtomsNP( "CH4.xyz")
#jobs.Test_projectDens( atomType=Zs, atomPos=xyzs )
ngrid=(128,64,32)
dcell = [0.2,0.2,0.2,0.2]
iA=0; iC=1


#ocl.tryInitFFT( ngrid)           ;print( "DEBUG poisson 1 " )
jobs.projectDens( iOutBuff=iA, atomType=Zs, atomPos=xyzs, iMO0=1, ngrid=ngrid, dcell=dcell )
#ocl.poisson   (  iA,iC, dcell )  ;print( "DEBUG poisson 3 " )

ibuff_MCOconv = 2
ibuff_densCO  = ocl.newFFTbuffer( "dens_CO" )
ibuff_MCOconv = ocl.newFFTbuffer( "MCOconv" )
ocl.loadFromBin( "../test_CO/dens.bin", ibuff_densCO )
ocl.convolve( iA,ibuff_densCO, ibuff_MCOconv  )
ocl.saveToXsf( "MCOconv.xsf", ibuff_MCOconv )





'''
ocl.tryInitFFT( ngrid ) 
V     = oclu.poisson( rho, dcell,  iA=0, iC=1 )
#ocl.tryInitFFT( A.shape )     ;print( "DEBUG poisson 1 " )
#print( "print Vmin Vmax ", np.min(V), np.max(V) )
'''

exit()


#jobs.Test_fft2d(n=23)   # WORKS for non-2^N arrays
#exit()

sys.path.append("../../")
import pyBall as pb
from   pyBall import FireCore as fc

atomPos,atomType,enames,qs = au.loadAtomsNP( "input.bas")

iMO = 1

fc.run_nonSCF( atomType, atomPos )

'''
ngrid, dCell, lvs = fc.setupGrid()
#fc.orb2xsf(iMO); #exit()
#fc.dens2xsf(); #exit()
#wf = fc.getGridMO( iMO,ngrid=ngrid)
rho = fc.getGridDens( ngrid=ngrid[::-1] )
print( " !!!!!! ngrid ", ngrid )
ocl.setGridShape_dCell( ngrid, dCell )
'''

ngrid=(128,64,32)     # we can only do multiples of 2^N
dcell = [0.1,0.1,0.1,0.1]
dCell = np.array([[0.1,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]],dtype=np.int32)

ocl.init()         
ocl.initFFT( ngrid )                  
ocl.loadWfBasis( [1,6], Rcuts=[4.5,4.5] )

ocl.setGridShape_dCell( ngrid, dCell )
ocl.tryInitFFT( ngrid )     ;print( "DEBUG poisson 1 " )

i0orb  = oclu.countOrbs( atomType )                  ;print("i0orb ", i0orb)  
wfcoef = fc.get_wfcoef(norb=i0orb[-1])               ;print("DEBUG 1");
ocl.convCoefsC( atomType, [2,1,1,1,1], atomPos, wfcoef,  iorb0=1, iorb1=2 , bInit=True )    ;print("DEBUG 2");
ocl.projectAtomsDens( 0 )            ;print("DEBUG 3");
ocl.saveToXsf( "rho_gpu.xsf", 1 )    ;print("DEBUG 4");

exit()

V     = oclu.poisson( rho, dcell,  iA=0, iC=1 )
ocl.tryInitFFT( A.shape )     ;print( "DEBUG poisson 1 " )
print( "print Vmin Vmax ", np.min(V), np.max(V) )
#ocl.saveToXsf( "V.xsf", 1 );


#ewfaux = fc.getGridMO( iMO,ngrid=ngrid)   ;print( "ewfaux.min(),ewfaux.max() ", ewfaux.min(),ewfaux.max() )
#sh = ewfaux.shape                         ;print( "ewfaux.shape ", sh )

'''
#i0orb  = jobs.countOrbs( atomType )           ;print("i0orb ", i0orb)  
#wfcoef = fc.get_wfcoef(norb=i0orb[-1])


#print( xyzs )

#jobs.Test_projectDens( )
'''