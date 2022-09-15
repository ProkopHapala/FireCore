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
from pyOCL import PP as pp

def test_PP_scan_LJQ(iZPP=8, nx=200, ny=100, nz=20, dtip=0.1, bUseBuffFE=False):
    ocl.setErrorCheck( 1 )
    print( "# --- Preparation" )
    apos,Zs,enames,qs = au.loadAtomsNP( "pentacene.xyz")
    ngrid=(128,64,64)
    dcell = [0.2,0.2,0.2,0.2]
    dCell = np.array([[dcell[0],0.0,0.0],[0.0,dcell[1],0.0],[0.0,0.0,dcell[2]]] )
    cell = [ngrid[0]*dcell[0],ngrid[1]*dcell[1],ngrid[2]*dcell[2]]
    p0   = np.array([cell[0]*-0.5,cell[1]*-0.5,cell[2]*-0.5 + 3.0])
    p0shift = p0*-1.; 
    p0shift[0]-=10.0
    p0shift[1]-=5.0
    p0shift[2]+=8.0
    iA=0; iC=1
    Ns= (ngrid[0],ngrid[1],ngrid[2])

    print( "#======= Allocate GPU memory " );
    itex_FE  = ocl.initPP( Ns )
    iBuffOut = ocl.newFFTbuffer( "OutFE", nfloat=4, ntot=nx*ny*nz )

    print( "#======= Make Force-Field " );
    ocl.setGridShapePP ( dCell, p0=p0 )
    #ocl.setGridShapePP ( dCell, p0=[.0,.0,.0] )
    typeParams = pp.loadSpecies('atomtypes.ini')
    cLJs       = pp.getAtomsLJ( iZPP, Zs, typeParams )
    if bUseBuffFE:
        ibuff_FE = ocl.newFFTbuffer( "FE" , 4 )
        ocl.evalLJC_QZs( ibuff_FE, apos, cLJs, Qs=qs )
        ocl.saveToXsf( "FE_w.xsf", ibuff_FE, stride=4, offset=3 )
        ocl.copyBuffToImage( ibuff_FE, itex_FE, ngrid[0],ngrid[1],ngrid[2] )    #;print("DEBUG 8 ")
    else:
        ocl.evalLJC_QZs_toImg( apos, cLJs, Qs=qs )


    print( "#======= Relaxed Scan " );
    ocl.makeStartPointGrid( nx, ny, p0shift, [0.1,0.0,0.0], [0.0,0.1,0.0] ) #;print("DEBUG 9 ")
    #ocl.getFEinStrokes ( iBuffOut, nz, [0.0,0.0,0.2] )                     #;print("DEBUG 10")
    ocl.relaxStrokesTilted( iBuffOut )                                      #;print("DEBUG 10")
    #ocl.getFEinStrokes ( iBuffOut, nz, [0.0,0.0,0.2] )                    #;print("DEBUG 10")
    OutFE = ocl.download  ( iBuffOut, Ns=(nx,ny,nz,4), dtype=np.float32 )   #;print("DEBUG 11")
    #print( "OutFE.shape ", OutFE.shape )

    print( "#======= Plot " );
    import matplotlib.pyplot as plt
    #print( "OutFE[:,:,:,0] \n", OutFE[:,:,:,0] )
    #print( "OutFE[:,:,:,1] \n", OutFE[:,:,:,1] )
    #print( "OutFE[:,:,:,2] \n", OutFE[:,:,:,2] )
    nnz = 5
    dnz = nz//nnz
    plt.figure(figsize=(5*4,4*nnz))
    cmap = "gray"
    for i in range(nnz):
        iz = i*dnz
        plt.subplot(5,4,i*4+1); plt.imshow( OutFE[:,:,iz,0], cmap=cmap ); plt.colorbar() ;plt.title("iz=%i" %iz)
        plt.subplot(5,4,i*4+2); plt.imshow( OutFE[:,:,iz,1], cmap=cmap ); plt.colorbar()
        plt.subplot(5,4,i*4+3); plt.imshow( OutFE[:,:,iz,2], cmap=cmap ); plt.colorbar()
        plt.subplot(5,4,i*4+4); plt.imshow( OutFE[:,:,iz,3], cmap=cmap ); plt.colorbar()
    plt.show()

def test_PP_makeFF_LJQ(iZPP=8):
    ocl.setErrorCheck( 1 )
    print( "# --- Preparation" )
    apos,Zs,enames,qs = au.loadAtomsNP( "pentacene.xyz")
    ngrid=(128,64,32)
    dcell = [0.2,0.2,0.2,0.2]
    dCell = np.array([[dcell[0],0.0,0.0],[0.0,dcell[1],0.0],[0.0,0.0,dcell[2]]] )
    p0 = [ngrid[0]*dcell[0]*-0.5,ngrid[1]*dcell[1]*-0.5,ngrid[2]*dcell[2]*-0.5 + 3.0]
    iA=0; iC=1
    Ns= (ngrid[0],ngrid[1],ngrid[2])

    itex_FE  = ocl.initPP( Ns )    #;print("DEBUG 1 itex_FE ", itex_FE )  
    ibuff_FE = ocl.newFFTbuffer( "FE" , 4 )   #;print("DEBUG 3 ")
    ocl.setGridShapePP ( dCell, p0=p0 )
    #ocl.setGridShapePP ( dCell, p0=[.0,.0,.0] )
    typeParams = pp.loadSpecies('atomtypes.ini')
    cLJs       = pp.getAtomsLJ( iZPP, Zs, typeParams )
    ocl.evalLJC_QZs( ibuff_FE, apos, cLJs, Qs=qs )
    ocl.saveToXsf( "FE_w.xsf", ibuff_FE, stride=4, offset=3 )

def test_PP_sampleFF():
    ocl.setErrorCheck( 1 )
    print( "# --- Preparation" )
    apos,Zs,enames,qs = au.loadAtomsNP( "pentacene.xyz")
    ngrid=(128,64,32)
    dcell = [0.2,0.2,0.2,0.2]
    dCell = np.array([[dcell[0],0.0,0.0],[0.0,dcell[1],0.0],[0.0,0.0,dcell[2]]] )
    iA=0; iC=1
    Ns= (ngrid[0],ngrid[1],ngrid[2])

    itex_FE = ocl.initPP( Ns )    #;print("DEBUG 1 itex_FE ", itex_FE )  
    
    print( "# --- SCF density")
    ibuff_rho  = ocl.newFFTbuffer( "rho", 4 )   #;print("DEBUG 2 ")
    ibuff_FE   = ocl.newFFTbuffer( "FE" , 4 )   #;print("DEBUG 3 ")
    jobs.projectDens( iOutBuff=ibuff_rho, atomType=Zs, atomPos=apos, iMO0=0, ngrid=ngrid, dcell=dcell, bSaveXsf=False, bSaveBin=False, bSCF=True, bDen0diff=False )
    ocl.gradient( ibuff_rho,  ibuff_FE, dcell)
    ocl.saveToXsf( "F_x.xsf", ibuff_FE, stride=4, offset=0 )
    ocl.saveToXsf( "F_y.xsf", ibuff_FE, stride=4, offset=1 )
    ocl.saveToXsf( "F_z.xsf", ibuff_FE, stride=4, offset=2 )
    ocl.saveToXsf( "F_w.xsf", ibuff_FE, stride=4, offset=3 )
    #ocl.ocl.saveToXsf( "test.xsf", ibuff_FE, stride=4, offset=0 )

    #nx=5;ny=5;nz=5
    nx=100;ny=100;nz=5
    iBuffOut = ocl.newFFTbuffer( "OutFE", nfloat=4, ntot=nx*ny*nz )      #;print("DEBUG 6 ")
    ocl.copyBuffToImage( ibuff_FE, itex_FE, ngrid[0],ngrid[1],ngrid[2] ) #;print("DEBUG 7 ")
    #ocl.copyBuffToImage( ibuff_FE, itex_FE, ngrid[2],ngrid[1],ngrid[0] ) ;print("DEBUG 7 ")
    ocl.setGridShapePP ( dCell, p0=None )                                #;print("DEBUG 8 ")    
    #ocl.makeStartPointGrid( nx, ny, [2.0,2.0,6.0], [0.1,0.0,0.0], [0.0,0.1,0.0] ) ;print("DEBUG 9 ")
    ocl.makeStartPointGrid( nx, ny, [0.0,0.0,4.0], [1.0/nx,0.0,0.0], [0.0,1.0/ny,0.0] ) #;print("DEBUG 9 ")
    #ocl.getFEinStrokes ( iBuffOut, nz )                                  #;print("DEBUG 10")
    ocl.getFEinStrokes ( iBuffOut, nz, [0.0,0.0,0.2] )                    #;print("DEBUG 10")
    OutFE = ocl.download( iBuffOut, Ns=(nx,ny,nz,4), dtype=np.float32 )   #;print("DEBUG 11")

    import matplotlib.pyplot as plt
    #print( "OutFE[:,:,:,0] \n", OutFE[:,:,:,0] )
    #print( "OutFE[:,:,:,1] \n", OutFE[:,:,:,1] )
    #print( "OutFE[:,:,:,2] \n", OutFE[:,:,:,2] )
    plt.figure(figsize=(5*4,4*nz))
    for i in range(nz):
        plt.subplot(5,4,i*4+1); plt.imshow( OutFE[:,:,i,0] ); plt.colorbar()
        plt.subplot(5,4,i*4+2); plt.imshow( OutFE[:,:,i,1] ); plt.colorbar()
        plt.subplot(5,4,i*4+3); plt.imshow( OutFE[:,:,i,2] ); plt.colorbar()
        plt.subplot(5,4,i*4+4); plt.imshow( OutFE[:,:,i,3] ); plt.colorbar()
    plt.show()

def test_job_Density_Gradient():
    ocl.setErrorCheck( 1 )
    print( "# --- Preparation" )
    apos,Zs,enames,qs = au.loadAtomsNP( "pentacene.xyz")
    ngrid=(128,64,32)
    dcell = [0.2,0.2,0.2,0.2]
    iA=0; iC=1

    print( "# --- SCF density")
    jobs.projectDens( iOutBuff=iA, atomType=Zs, atomPos=apos, iMO0=0, ngrid=ngrid, dcell=dcell, bSaveXsf=False, bSaveBin=False, bSCF=True, bDen0diff=False )
    ibuff_FE  = ocl.newFFTbuffer( "FE", 4 )
    ocl.gradient( iA, ibuff_FE, dcell)
    ocl.saveToXsf( "F_x.xsf", ibuff_FE, stride=4, offset=0 )
    ocl.saveToXsf( "F_y.xsf", ibuff_FE, stride=4, offset=1 )
    ocl.saveToXsf( "F_z.xsf", ibuff_FE, stride=4, offset=2 )
    ocl.saveToXsf( "F_w.xsf", ibuff_FE, stride=4, offset=3 )
    #ocl.ocl.saveToXsf( "test.xsf", ibuff_FE, stride=4, offset=0 )

def job_make_Eelec_Epauli():
    ocl.setErrorCheck( 1 )
    print( "# --- Preparation" )
    apos,Zs,enames,qs = au.loadAtomsNP( "pentacene.xyz")
    ngrid=(128,64,32)
    dcell = [0.2,0.2,0.2,0.2]
    iA=0; iC=1
    print( "# --- Allocations")
    iMO1 = 51
    #jobs.projectDensFireball( atomType=Zs, atomPos=apos, bSCF=True, saveXsf=1, f_den0=-1.0 )
    #jobs.density_from_firecore_to_xsf( atomType=Zs, atomPos=apos, bSCF=False, saveXsf="dens_check.xsf", Cden=1.0, Cden0=-1.0  )
    #jobs.projectDens( iOutBuff=iA, atomType=Zs, atomPos=apos, iMO0=1, iMO1=iMO1, ngrid=ngrid, dcell=dcell, bSaveXsf=True, bSaveBin=False, saveName="dens_scf", bSCF=True, bDen0diff=False )
    #jobs.projectDens( iOutBuff=iA, atomType=Zs, atomPos=apos, iMO0=1, iMO1=iMO1, ngrid=ngrid, dcell=dcell, bSaveXsf=True, bSaveBin=False, saveName="dens_scf", bSCF=True, bDen0diff=True )
    #jobs.check_density_projection( atomType=Zs, atomPos=apos, ngrid=ngrid, dcell=dcell, bSCF=False, iOutBuff=0, Cden=1.0, Cden0=0.0, iMO1=iMO1 )
    #dCell = np.array([[dcell[0],0.0,0.0],[0.0,dcell[1],0.0],[0.0,0.0,dcell[2]]],dtype=np.float32)
    #jobs.orbitals_from_firecore( orbitals=[50,51,52,53], atomType=Zs, atomPos=apos, ngrid=ngrid, dCell=dCell, g0=None )
    #exit(0)

    print( "# --- SCF density")
    jobs.projectDens( iOutBuff=iA, atomType=Zs, atomPos=apos, iMO0=0, ngrid=ngrid, dcell=dcell, bSaveXsf=False, bSaveBin=False, bSCF=True, bDen0diff=False )
    ibuff_densCO  = ocl.newFFTbuffer( "dens_CO" )
    ibuff_ConvOut = ocl.newFFTbuffer( "MCOconv" )
    ibuff_DensBak = ocl.newFFTbuffer( "DensBak" )
    ocl.copy( iA, ibuff_DensBak )
    #ocl.roll( iA, ibuff_ConvOut, [ngrid[0]//2,ngrid[1]//2,ngrid[2]//2] )
    #ocl.saveToXsf( "Dens_bak.xsf",  ibuff_DensBak )
    #ocl.saveToXsf( "Dens_orig.xsf", iA)
    #ocl.saveToXsf( "Dens_roll.xsf", ibuff_ConvOut )

    print( "# ==== E_Pauli ( density convolution )")
    ocl.loadFromBin( "../test_CO/dens_scf.bin", ibuff_densCO )
    ocl.convolve( ibuff_DensBak,ibuff_densCO, ibuff_ConvOut )
    ocl.saveToXsf( "Epaul.xsf",    ibuff_ConvOut )
    
    print( "# ==== E_Electrostatic ( density convolution )")
    iBuffDens0 = iA
    jobs.projectDens0_new( iOutBuff=iBuffDens0, atomPos=apos, atomType=Zs, ngrid=ngrid, dcell=dcell, bSaveXsf=False, bSaveBin=False, acumCoef=[1.0,-1.0] )
    Ns = (ngrid[0],ngrid[1],ngrid[2])
    data = ocl.download( iBuffDens0, data=None, Ns=Ns )
    dvol = dcell[0]*dcell[1]*dcell[2]
    Q_diff = np.sum(data  )*dvol; print( "Q_diff ", Q_diff )
    #ocl.saveToXsfAtomsData( "dens_diff.xsf",  data, Zs, apos )
    ocl.saveToXsf( "dens_diff.xsf", iBuffDens0 )
    print( "# --- Poisson (rho->V)")
    ocl.poisson( iA=iBuffDens0, iOut=iC, dcell=dcell )
    ocl.saveToXsf( "Vout.xsf", iC )

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
#job_make_Eelec_Epauli()
#test_job_Density_Gradient()
#test_PP_sampleFF()
#test_PP_makeFF_LJQ()
test_PP_scan_LJQ()

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
