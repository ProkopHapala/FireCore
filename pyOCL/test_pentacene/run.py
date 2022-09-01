#!/usr/bin/python

# INTEROPERABLE GLOBAL VARIABLES
# https://gcc.gnu.org/onlinedocs/gfortran/Interoperable-Global-Variables.html#Interoperable-Global-Variables

import sys
import numpy as np
import os

sys.path.append("../")
from pyOCL import oclfft as ocl
from pyOCL import utils as oclu
from pyOCL import high_level as oclh
from pyOCL import jobs
from pyOCL import atomicUtils as au

def job_make_Eelec_Epauli():
    ocl.setErrorCheck( 1 )
    print( "# --- Preparation" )
    apos,Zs,enames,qs = au.loadAtomsNP( "pentacene.xyz")
    ngrid=(128,64,32)
    dcell = [0.2,0.2,0.2,0.2]
    iA=0; iC=1
    print( "# --- Allocations")

    #jobs.projectDensFireball( atomType=Zs, atomPos=apos, bSCF=True, saveXsf=1, f_den0=-1.0 )
    #jobs.density_from_firecore_to_xsf( atomType=Zs, atomPos=apos, bSCF=False, saveXsf="dens_check.xsf", Cden=1.0, Cden0=-1.0  )
    jobs.projectDens( iOutBuff=iA, atomType=Zs, atomPos=apos, iMO0=1, ngrid=ngrid, dcell=dcell, bSaveXsf=True, bSaveBin=False, saveName="dens_scf", bSCF=True )
    #jobs.check_density_projection(  atomType=Zs, atomPos=apos, ngrid=ngrid, dcell=dcell, bSCF=False, iOutBuff=0, Cden=1.0, Cden0=-1.0 )

    exit(0)


    print( "# --- SCF density")
    jobs.projectDens( iOutBuff=iA, atomType=Zs, atomPos=apos, iMO0=1, ngrid=ngrid, dcell=dcell, bSaveXsf=False, bSaveBin=False, saveName="dens_scf", bSCF=True )
    
    ibuff_densCO  = ocl.newFFTbuffer( "dens_CO" )
    ibuff_ConvOut = ocl.newFFTbuffer( "MCOconv" )
    ibuff_DensBak = ocl.newFFTbuffer( "DensBak" )

    ocl.copy( iA, ibuff_DensBak )
    ocl.saveToXsf( "Dens_bak.xsf",  ibuff_DensBak )
    ocl.saveToXsf( "Dens_orig.xsf", iA)

    print( "# ==== E_Pauli ( density convolution )")
    ocl.loadFromBin( "../test_CO/dens_scf.bin", ibuff_densCO )
    ocl.convolve( ibuff_DensBak,ibuff_densCO, ibuff_ConvOut )
    ocl.saveToXsf( "Epaul.xsf",    ibuff_ConvOut )

    #print( "# === E_elec ( density and potential convolution )")
    #jobs.projectDens0( iOutBuff=ibuff_DensBak, atomType=Zs, atomPos=apos, ngrid=ngrid, dcell=dcell, bSaveXsf=False,  bSaveBin=False, saveName="dens_diff" )   ;print( "DEBUG 2.1 " )
    #ocl.saveToXsf( "dens_diff.xsf", ibuff_DensBak )         ;print( "DEBUG 2.2 " )
    
    iBuffDens0 = iA
    #jobs.projectDens0( iOutBuff=iBuffDens0, atomType=Zs, atomPos=apos, ngrid=ngrid, dcell=dcell, bSaveXsf=False,  bSaveBin=False, saveName="dens_diff" )   ;print( "DEBUG 2.1 " )
    jobs.projectDens0_new( iOutBuff=iBuffDens0, atomPos=apos, atomType=Zs, ngrid=ngrid, dcell=dcell, bSaveXsf=False, bSaveBin=False, saveName="dens_diff" )
    ocl.saveToXsf( "dens_diff.xsf", iBuffDens0 )
    print( "# --- Poisson (rho->V)")
    ocl.poisson( iA=iBuffDens0, iOut=iC, dcell=dcell )
    ocl.saveToXsf( "Vout.xsf", iC )

    exit()
    print( "# --- E_elec = convolution( rho, V )  " )
    ocl.loadFromBin( "../test_CO/dens_diff.bin", ibuff_densCO )
    ocl.convolve( iC,ibuff_densCO, ibuff_ConvOut )
    ocl.saveToXsf( "E_elec.xsf",   ibuff_ConvOut )

def job_convolve_density_with_CO_orig( ):
    ocl.setErrorCheck( 1 )
    xyzs,Zs,enames,qs = au.loadAtomsNP( "pentacene.xyz")
    ngrid=(128,64,32)
    dcell = [0.2,0.2,0.2,0.2]
    iA=0; iC=1
    jobs.projectDens( iOutBuff=iA, atomType=Zs, atomPos=xyzs, iMO0=1, ngrid=ngrid, dcell=dcell )
    ibuff_MCOconv = 2
    ibuff_densCO  = ocl.newFFTbuffer( "dens_CO" )
    ibuff_MCOconv = ocl.newFFTbuffer( "MCOconv" )
    ocl.loadFromBin( "../test_CO/dens.bin", ibuff_densCO )
    ocl.convolve( iA,ibuff_densCO, ibuff_MCOconv  )
    ocl.saveToXsf( "MCOconv.xsf", ibuff_MCOconv )

def project_or_load_density( ngrid, iBuff=0, dcell=None ):
    if dcell is None: dcell = [0.2,0.2,0.2,0.2]
    if not os.path.exists( "dens.bin" ):
        print("!!!!! job_convolve_density_with_CO :  PROJECT+SAVE ./dens.bin ")
        xyzs,Zs,enames,qs = au.loadAtomsNP( "pentacene.xyz")
        jobs.projectDens( iOutBuff=iBuff, atomType=Zs, atomPos=xyzs, iMO0=1, ngrid=ngrid, dcell=dcell, bSaveXsf=False, bSaveBin=True )
    else:
        print("!!!!! job_convolve_density_with_CO :  LOAD ./dens.bin ")
        Ns = (ngrid[0],ngrid[1],ngrid[2])
        ocl.initFFTgrid( Ns, dcell=dcell )
        ocl.loadFromBin( "./dens.bin", iBuff )

def job_convolve_density_with_CO( iA=0 ):
    print( "JOB: convolve_density_with_CO() " )
    ocl.setErrorCheck( 1 )
    ngrid=(128,64,32)
    project_or_load_density( ngrid, iBuff=iA )
    ibuff_densCO  = ocl.newFFTbuffer( "dens_CO" )
    ibuff_MCOconv = ocl.newFFTbuffer( "MCOconv" )
    ocl.loadFromBin( "../test_CO/dens.bin", ibuff_densCO )
    ocl.convolve( iA,ibuff_densCO, ibuff_MCOconv  )
    ocl.saveToXsf( "MCOconv.xsf", ibuff_MCOconv )

def job_poisson_equation( iA=0, iC=1 ):
    ocl.setErrorCheck( 1 )
    ngrid = (128,64,32)
    dcell = [0.2,0.2,0.2,0.2]
    project_or_load_density( ngrid, iBuff=iA, dcell=dcell )
    ocl.poisson( iA=iA, iOut=iC, dcell=dcell )
    ocl.saveToXsf( "Vout.xsf", iC )

#job_convolve_density_with_CO_orig()
#job_convolve_density_with_CO()
#job_poisson_equation()
job_make_Eelec_Epauli()

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
