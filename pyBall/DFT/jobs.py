#!/usr/bin/python

import numpy as np
#import os
import sys
from . import oclfft as ocl
from . import utils  as oclu

def projectAtoms__( atoms, acoefs):
    import matplotlib.pyplot as plt
    import time
    print( "# ========= TEST   projectAtoms__()  " )
    ocl.init()
    Ns=(100,100,100)
    ocl.initFFT( Ns  )
    ocl.loadWfBasis( [1,6], Rcuts=[4.5,4.5] )
    ocl.initAtoms( len(atoms) )
    #initBasisTable( basis.shape[0], basis.shape[1], basis )

    ocl.setGridShape( )
    t0 = time.clock()
    ocl.projectAtoms( atoms, acoefs, 0 )
    ocl.saveToXsf( "test.xsf", 0 )
    arrA = np.zeros(Ns,dtype=np.csingle  )
    ocl.download ( 0, arrA )    
    t = time.clock()-t0; print( "projection time ", t )
    plt.figure(); 
    #print( arrA[10,10].real )
    plt.imshow( arrA[10].real ) 
    #plt.imshow( np.log( np.abs(arrA[10])) ) 
    plt.grid()
    plt.show(); 
    ocl.cleanup()

def Test_projectAtoms(n=64, natoms=1000):
    import matplotlib.pyplot as plt
    import time
    print( "# ========= TEST   Test_projectAtoms()  " )
    
    acs=[
    [[0.0,2.0,0.0,1.001],    [0.0,1.0,0.0,0.0]],  
    [[6.0,2.0,0.0,0.999],    [0.0,0.0,0.0,1.0]],
    #[[2.0, 2.0,0.0,1.0],    [1.0,0.0,0.0,1.0]], 
    #[[-1.0, 1.0,0.0,1.0],    [1.0,0.0,0.0,1.0]], 
    #[[-1.0, 1.0,0.0,1.0],    [1.0,0.0,0.0,1.0]], 
    #[[-1.0,-1.0,0.0,1.0],    [1.0,0.0,0.0,1.0]], 
    #[[2.0,-1.0,0.0,1.0],    [1.0,0.0,0.0,1.0]], 
    #[[2.0,-2.0,0.0,1.0],    [1.0,0.0,0.0,1.0]], 
    ]
    atoms = np.array([ a[0] for a in acs ], dtype=np.float32)
    coefs = np.array([ a[1] for a in acs ], dtype=np.float32)
    
    #atoms = np.random.rand(natoms,4).astype(np.float32);    atoms[:,:3]*=10.0; atoms[:,3]=3.0;
    #coefs = np.random.rand(natoms,4).astype(np.float32);    coefs[:,:3] = 0.0; coefs[:,3]=1.0; 

    '''
    ny=32
    nx=128
    basis        = np.zeros( (ny,nx,4) )  # RGBA   resp {x,y,z,s}
    basis[0:8 ,:3] = 1.0
    basis[8:32,:3] = 1.0
    basis[0:8,10:20] = 1.0
    basis[0:8,40:50] = 1.0
    basis[8:32,30:40] = 1.0
    basis[8:32,40:50] = -1.0
    basis[8:32,50:60] = 1.0
    basis=basis.astype(np.float32)
    plt.imshow(basis[:,:,0], interpolation='nearest'); #plt.show()
    '''

    #print( "atoms ", atoms, atoms.dtype )
    #print( "coefs ", coefs )
    ocl.init()
    #Ns=(128,128,128)
    Ns=(100,100,100)
    ocl.initFFT( Ns  )
    ocl.loadWfBasis( [1,6], Rcuts=[4.5,4.5] )
    ocl.initAtoms( len(atoms) )
    #initBasisTable( basis.shape[0], basis.shape[1], basis )

    ocl.setGridShape()
    t0 = time.clock()
    ocl.projectAtoms( atoms, coefs, 0 )
    ocl.saveToXsf( "test.xsf", 0 )
    arrA = np.zeros(Ns,dtype=np.csingle  )
    ocl.download ( 0, arrA )    
    t = time.clock()-t0; print( "projection time ", t )
    plt.figure(); 
    #print( arrA[10,10].real )
    plt.imshow( arrA[10].real ) 
    #plt.imshow( np.log( np.abs(arrA[10])) ) 
    plt.grid()
    plt.show(); 
    ocl.cleanup()

def Test_projectAtomPosTex():
    import matplotlib.pyplot as plt
    import time
    print( "# ========= Test_projectAtomPosTex  " )
    
    xref = np.linspace(0.0,4.50*ocl.const_Bohr_Radius,1000)
    yref = ocl.loadWf_C( "Fdata/basis/001_450.wf1", n=1000 )
    #exit()

    # 4.5 * const_Bohr_Radius = 2.38129744906  #[A]

    acs=[
    [1,  [0.0,0.0,0.0,0.001],    [0.0,0.0,0.0,1.0]],  
    #[6, [6.0,2.0,0.0,1.001],    [0.0,0.0,0.0,1.0]],
    ]

    atomType = np.array( [ a[0] for a in acs ], dtype=np.int32 )
    atomPos  = np.array( [ a[1] for a in acs ] )
    #apos     = np.array( [ a[1] for a in acs ], dtype=np.float32 )
    #coefs    = np.array( [ a[2] for a in acs ], dtype=np.float32) 
    
    npos = 100   # must be multiple of local group size
    xs = np.linspace(0.0,5.0,npos)
    poss  = np.zeros((npos,3)                 ); poss [:,0]= xs
    poss_ = np.zeros((npos,4),dtype=np.float32); poss_[:,0]= xs

    # ============= CPU Fortran
    fc.preinit()
    fc.init( atomType, atomPos )
    # =========== Electron Density
    fc.assembleH( atomPos)
    fc.solveH()
    sigma=fc.updateCharges() ; print( sigma )
    ngrid, dCell, lvs = fc.setupGrid()

    ys = fc.getpsi( poss, in1=1, issh=1, l=0, m=1 )
    #print( "Test_projectAtomPosTex ys ", ys )

    wfcoef = fc.get_wfcoef(norb=1)
    #print( "Test_projectAtomPosTex wfcoef: \n", wfcoef )
    apos_,wfcoef_ = oclu.convCoefs( atomType, wfcoef, atomPos )
    # ========== GPU
    ocl.init()                                   #;print("DEBUG py 2")
    ocl.loadWfBasis( [1,6], Rcuts=[4.5,4.5] )         #;print("DEBUG py 3") 
    ocl.initAtoms( len(apos_) )                  #;print("DEBUG py 4")
    #initBasisTable( basis.shape[0], basis.shape[1], basis )
    out = ocl.projectAtomPosTex(apos_,wfcoef_,poss_) #;print("DEBUG py 5")
    plt.plot( xref, yref*ocl.pref_s,    label="load" ) 
    plt.plot( poss_[:,0], out[:,0], label="GPU" ) 
    plt.plot( poss [:,0], ys,       label="CPU", ls="--" ) 
    #plt.imshow( np.log( np.abs(arrA[10])) ) 
    plt.legend()
    plt.grid()
    plt.show(); 
    ocl.cleanup()

def Test_fft2d(n=64):
    import matplotlib.pyplot as plt
    print( "# ========= TEST   runFFT()  " )
    xs    = np.linspace(-2.0,2.0,n)
    Xs,Ys = np.meshgrid(xs,xs)
    #arrA   = np.sin(Xs*3 + (Xs**2)*5.1 )*np.cos(Ys*20 + (Xs**2)*3.0)   + np.random.rand(n,n)*0.1
    arrA   = 1/( np.sin(Xs*2)**2 + np.cos(Ys*3 )**2 + 0.1 ) + np.random.rand(n,n)*0.1
    arrA   = arrA.astype( np.csingle )
    arrB   = ( np.sin(Xs*6)**2 + np.cos(Ys*3 )**2 + 0.1 )/np.exp( -0.1*(Xs**2 + Ys**2) )  + np.random.rand(n,n)*0.1
    arrB   = arrB.astype( np.csingle )
    arrC = arrA.copy(); arrC[:,:]=0
    ocl.init()
    ocl.initFFT(    arrA.shape )
    ocl.upload ( 0, arrA )    ;plt.figure(); plt.imshow( arrA.real )
    ocl.upload ( 1, arrB )    ;plt.figure(); plt.imshow( arrB.real )
    #convolve( 0,1,   2 )
    #arrC = np.zeros(arrA.shape)
    #download( 2, arrC)    ;plt.figure(); plt.imshow( arrC.real ) 
    ocl.runfft (0 ); ocl.download ( 0, arrA )    ;plt.figure(); plt.imshow( np.log( np.abs(arrA)) ) 
    ocl.runfft (1 ); ocl.download ( 1, arrB )    ;plt.figure(); plt.imshow( np.log( np.abs(arrB)) ) 
    plt.show(); 
    ocl.cleanup()

def Test_Convolution2d( n=1024):
    import matplotlib.pyplot as plt
    print( "# ========= TEST   convolve()  " )
    xs    = np.linspace(-2.0,2.0,n)
    Xs,Ys = np.meshgrid(xs,xs)
    #arrA   = np.sin(Xs*3 + (Xs**2)*5.1 )*np.cos(Ys*20 + (Xs**2)*3.0)   + np.random.rand(n,n)*0.1
    arrA   = 1/( np.sin(Xs*2)**2 + np.cos(Ys*3 )**2 + 0.1 )
    arrA   = arrA.astype( np.csingle )
    arrB   = (Xs*Ys)/np.exp( -0.55*(Xs**2 + Ys**2) )  + np.random.rand(n,n)*0.1
    arrB   = arrB.astype( np.csingle )
    arrC = arrA.copy(); arrC[:,:]=0
    ocl.init()
    ocl.initFFT(    arrA.shape )
    ocl.upload ( 0, arrA )    ;plt.figure(); plt.imshow( arrA.real )
    ocl.upload ( 1, arrB )    ;plt.figure(); plt.imshow( arrB.real )
    ocl.upload ( 2, arrC ) 
    ocl.convolve( 0,1,2 )
    ocl.download( 2, arrC)    ;plt.figure(); plt.imshow( arrC.real )
    plt.show(); 
    ocl.cleanup()

def Test_projectWf( iMO=1 ):
    sys.path.append("../../")
    import pyBall as pb
    from pyBall import FireCore as fc
    # ----- Make CH4 molecule
    atomType = np.array([6,1,1,1,1]).astype(np.int32)
    atomPos  = np.array([
        [ 0.01,     0.0,     0.0],
        [-1.0,     +1.0,    -1.0],
        [+1.0,     -1.0,    -1.0],
        [-1.0,     -1.0,    +1.0],
        [+1.0,     +1.0,    +1.0],
    ])
    print("# ======== FireCore Run " )
    #print ("atomType ", atomType)
    #print ("atomPos  ", atomPos)
    fc.preinit()
    norb = fc.init( atomType, atomPos )
    # --------- Electron Density
    fc.assembleH( atomPos)
    fc.solveH()
    sigma=fc.updateCharges() ; print( sigma )
    ngrid, dCell, lvs = fc.setupGrid()
    ewfaux = fc.getGridMO( iMO,ngrid=ngrid)   #;print( "ewfaux.min(),ewfaux.max() ", ewfaux.min(),ewfaux.max() )
    sh = ewfaux.shape                         #;print( "ewfaux.shape ", sh )
    fc.orb2xsf(iMO); #exit()

    i0orb  = oclu.countOrbs( atomType )       #;print("i0orb ", i0orb)  
    wfcoef = fc.get_wfcoef(norb=i0orb[-1])
    print("# ========= PyOCL Wave-Function Projection " )
    #wfcoef = wfcoef[1,:]
    #print( "wfcoef: \n", wfcoef )
    #print( "coefs for MO#%i " %iMO, wfcoef[iMO-1,:] )
    typeDict={ 1:0, 6:1 }
    #apos_,wfcoef_ = convCoefs( atomType, wfcoef[iMO,:], atomPos )
    apos_,wfcoef_ = oclu.convCoefs( atomType, wfcoef[iMO-1,:], atomPos, typeDict )
    #wfcoef_ = np.array( [ 0.0,0.0,0.0,1.0,   0.0,0.0,0.0,1.0,   0.0,0.0,0.0,1.0,   0.0,0.0,0.0,1.0,   0.0,0.0,0.0,1.0 ], dtype=np.float32)
    ocl.init()                           
    Ns = (ngrid[0],ngrid[1],ngrid[2])
    ocl.initFFT( Ns  )                  
    ocl.loadWfBasis( [1,6], Rcuts=[4.5,4.5] )    
    ocl.initAtoms( len(apos_) )          
    ocl.setGridShape_dCell( Ns, dCell )
    ocl.projectAtoms( apos_, wfcoef_, 0 ) 
    #saveToXsf( "test.xsf", 0 )
    ocl.saveToXsfAtoms( "orb_%03i.xsf" %iMO, 0,    atomType, atomPos  )
    arrA = np.zeros(Ns+(2,),dtype=np.float32)  
    ocl.download(0,arrA)                  

    print("# ========= Plot GPU vs CPU comparison " )
    #print( arrA  [16,16,:,0] )
    #print( ewfaux[16,16,:] )
    import matplotlib.pyplot as plt
    ix0=ngrid[0]//2
    iy0=ngrid[1]//2
    plt.plot( arrA  [ix0,iy0,:,0], label="GPU" )
    plt.plot( ewfaux[ix0,iy0,:  ], label="CPU" )
    #print( "arrA  .shape ", arrA  .shape )
    #print( "ewfaux.shape ", ewfaux.shape )
    #plt.figure(); plt.imshow( arrA  [ngrid[0]//2 ,:,:,0] ); plt.title("GPU")
    #plt.figure(); plt.imshow( ewfaux[ngrid[0]//2 ,:,:  ] ); plt.title("CPU")
    plt.legend()
    plt.show()

def Test_projectDens( iMO0=1, iMO1=2 ):
    sys.path.append("../../")
    import pyBall as pb
    from pyBall import FireCore as fc
    # ----- Make CH4 molecule
    atomType = np.array([6,1,1,1,1]).astype(np.int32)
    atomPos  = np.array([
        [ 0.01,     0.0,     0.0],
        [-1.0,     +1.0,    -1.0],
        [+1.0,     -1.0,    -1.0],
        [-1.0,     -1.0,    +1.0],
        [+1.0,     +1.0,    +1.0],
    ])
    
    print("# ======== FireCore Run " )
    print ("atomType ", atomType)
    print ("atomPos  ", atomPos)
    fc.preinit()
    norb = fc.init( atomType, atomPos )
    # --------- Electron Density
    fc.assembleH( atomPos)
    fc.solveH()
    sigma=fc.updateCharges() ; print( sigma )
    # ======== Project Grid using FireCore "
    ngrid, dCell, lvs = fc.setupGrid()
    #ewfaux = fc.getGridMO( iMO,ngrid=ngrid)   ;print( "ewfaux.min(),ewfaux.max() ", ewfaux.min(),ewfaux.max() )
    #sh = ewfaux.shape                       ;print( "ewfaux.shape ", sh )
    #fc.orb2xsf(iMO); #exit()

    i0orb  = ocl.countOrbs( atomType )           ;print("i0orb ", i0orb)  
    wfcoef = fc.get_wfcoef(norb=i0orb[-1])

    
    print("# ========= PyOCL Density-Function Projection " )
    ocl.init()            
    Ns = (ngrid[0],ngrid[1],ngrid[2])
    ocl.initFFT( Ns  )                  
    ocl.loadWfBasis( [1,6], Rcuts=[4.5,4.5] )    
    #initAtoms( len(apos_) )          
    ocl.setGridShape_dCell( Ns, dCell )
    print( "wfcoef \n", wfcoef )
    ocl.convCoefsC( atomType, [2,1,1,1,1], atomPos, wfcoef,  iorb0=iMO0, iorb1=iMO1 , bInit=True ) 
    ocl.projectAtomsDens( 0 ) 

    print( "DEBUG before saveToXsfAtoms " )
    ocl.saveToXsfAtoms( "dens_%03i_%03i.xsf" %(iMO0,iMO1), 0,    atomType, atomPos  )

    exit(0)

    print( "DEBUG after saveToXsfAtoms " )
    arrA = np.zeros(Ns+(2,),dtype=np.float32)  
    download(0,arrA) 
                     

    print("# ========= Plot GPU vs CPU comparison " )
    #print( arrA  [16,16,:,0] )
    #print( ewfaux[16,16,:] )
    import matplotlib.pyplot as plt
    ix0=ngrid[0]//2
    iy0=ngrid[1]//2
    plt.plot( arrA  [ix0,iy0,:,0], label="GPU" )
    #plt.plot( ewfaux[ix0,iy0,:  ], label="CPU" )
    print( "arrA  .shape ", arrA  .shape )
    print( "ewfaux.shape ", ewfaux.shape )
    #plt.figure(); plt.imshow( arrA  [ngrid[0]//2 ,:,:,0] ); plt.title("GPU")
    #plt.figure(); plt.imshow( ewfaux[ngrid[0]//2 ,:,:  ] ); plt.title("CPU")
    plt.legend()
    plt.show()

def iZs2dict(iZs, dr="./Fdata/basis"):
    import glob
    s = { i for i in iZs }
    elems = sorted( s )   ;print( "elems ", elems )
    dct   = { k:(v+1) for v,k in enumerate(elems) }    
    ords  = [ dct[i] for i in iZs ]                   
    Rcuts = [ ]
    for i in elems:
        path = '%s/%03i_*.na0' %(dr,i) 
        lstr = glob.glob( path )[0].split('_')[1]
        rcut = float(lstr[:3])/100.0 
        print( i, path,  lstr, rcut )
        Rcuts.append( rcut )
    return elems, dct, ords, Rcuts

def prepare_fireball(  atomType=None, atomPos=None, bSCF=False ):
    sys.path.append("../../")
    import pyBall as pb
    from pyBall import FireCore as fc
    fc.setVerbosity(verbosity=0)
    fc.preinit()
    norb = fc.init( atomType, atomPos )
    # --------- Electron Density
    if bSCF:
        fc.setVerbosity(verbosity=1)
        fc.SCF( atomPos, iforce=0, nmax_scf=200 )
    else:
        fc.assembleH( atomPos)
        fc.solveH()
        sigma=fc.updateCharges() ; print( sigma )
    return fc

def orbitals_from_firecore( atomType=None, atomPos=None, bSCF=False, orbitals=[0,], ngrid=None, dCell=None, g0=None ):
    fc = prepare_fireball( atomType=atomType, atomPos=atomPos, bSCF=bSCF )
    ngrid, dCell, lvs = fc.setupGrid( ngrid=ngrid, dCell=dCell, g0=g0 )
    for iMO in orbitals:
        fc.orb2xsf( iMO )
        #ewfaux = fc.getGridDens( ngrid=ngrid, Cden=Cden, Cden0=Cden0 )

def density_from_firecore( atomType=None, atomPos=None, bSCF=False, Cden=1.0, Cden0=0.0, bGetGrid=False, bGetCoefs=False,  ngrid=None, dCell=None, g0=None ):
    fc = prepare_fireball( atomType=atomType, atomPos=atomPos, bSCF=bSCF )
    grid=None
    if bGetGrid:
        ngrid, dCell, lvs = fc.setupGrid( ngrid=ngrid, dCell=dCell, g0=g0 )
        ewfaux = fc.getGridDens( ngrid=ngrid, Cden=Cden, Cden0=Cden0 )
        grid = ( ewfaux, ngrid, dCell, lvs )
    coefs=None
    if bGetCoefs:
        i0orb  = oclu.countOrbs( atomType ) 
        wfcoef = fc.get_wfcoef(norb=i0orb[-1])
        coefs = (wfcoef,i0orb)
    
    return coefs, grid

def density_from_firecore_to_xsf( atomType=None, atomPos=None, bSCF=False, saveXsf="dens_check.xsf", Cden=1.0, Cden0=0.0 ):
    print( "DEBUG density_from_firecore_to_xsf()" )
    _, (ewfaux,ngrid,dCell,lvs) = density_from_firecore( atomType=atomType, atomPos=atomPos, bSCF=bSCF, Cden=Cden, Cden0=Cden0, bGetGrid=True, bGetCoefs=False )
    if saveXsf is not None:
        #ngrid, dCell, lvs = fc.setupGrid()
        #ewfaux = fc.getGridDens( ngrid=ngrid, Cden=1.0, Cden0=0.0 )
        #ewfaux = fc.getGridDens( ngrid=ngrid, Cden=Cden, Cden0=Cden0 )
        Ns = (ngrid[0],ngrid[1],ngrid[2])
        ocl.setGridShape_dCell( Ns, dCell, pos0=[0.0,0.0,0.0] )
        pos0 = ocl.getCellHalf( Ns, dCell );   print( "!!!!!!!! pos0 ", pos0  )
        atomPos[:,0]-=pos0[2]
        atomPos[:,1]-=pos0[1]
        atomPos[:,2]-=pos0[0]
        ocl.saveToXsfAtomsData( saveXsf, ewfaux, atomType, atomPos )

def project_dens_GPU( wfcoef, atomType=None, atomPos=None, ngrid=(64,64,64), dcell=[0.2,0.2,0.2,0.2], iOutBuff=0, iMO0=0, iMO1=None, i0orb=None, bDownalod=False, bDen0diff=False ):
    print( "#======= project_dens_GPU | ngrid: ", ngrid," dcell: ", dcell )
    if iMO1 is None:
        iMO1 = i0orb[-1]//2
        print( "project_dens_GPU iMO1 ", iMO1 )
    Ns = (ngrid[0],ngrid[1],ngrid[2])
    dCell = np.array([[dcell[0],0.0,0.0],[0.0,dcell[1],0.0],[0.0,0.0,dcell[2]]],dtype=np.float32)
    elems, dct, ords, Rcuts = iZs2dict(atomType)
    ocl.tryInitFFT( Ns )
    ocl.tryLoadWfBasis( elems, Rcuts=Rcuts )
    ocl.setGridShape_dCell( Ns, dCell )
    ocl.convCoefsC( atomType, ords, atomPos, wfcoef, bInit=True )
    ocl.projectAtomsDens( iOutBuff, iorb0=iMO0, iorb1=iMO1, acumCoef=[0.0,2.0] )     #  2.0 electrions per atom
    if bDen0diff:
        ords=np.array( ords, dtype=np.int32)
        #ocl.setTypes( [1,4], [[1.0,0.0],[1.0,3.0]] )  ;print( "WARRNING: project_dens_GPU(): ocl.setTypes([4,4],[[1.0,3.0],[1.0,5.0]]) is works just for Hydrocarbons !!!!" )
        ocl.setTypesZ( sorted(list(set(atomType))) )
        ocl.projectAtomsDens0( iOutBuff, apos=atomPos, atypes=ords, acumCoef=[1.0,-1.0] ) 
    if bDownalod:
        data = ocl.download( iOutBuff, data=None, Ns=Ns )
        data = data.real.astype(np.float)
        return data

def check_density_projection( atomType=None, atomPos=None, ngrid=(64,64,64), dcell = [0.2,0.2,0.2,0.2], bSCF=False, iOutBuff=0, Cden=1.0, Cden0=-1.0, iMO1=None ):
    print( "#======= check_density_projection()" )
    dCell = np.array([[dcell[0],0.0,0.0],[0.0,dcell[1],0.0],[0.0,0.0,dcell[2]]],dtype=np.float32)
    #g0 = [-12.0,-6.0,-3.0]
    g0 = None
    dvol = dcell[0]*dcell[1]*dcell[2]
    (wfcoef,i0orb), (ewfaux,ngrid_,dCell,lvs) = density_from_firecore( atomType=atomType, atomPos=atomPos, bSCF=bSCF, Cden=Cden, Cden0=Cden0, bGetGrid=True, bGetCoefs=True, ngrid=ngrid, dCell=dCell, g0=g0 )
    data = project_dens_GPU( wfcoef, atomType=atomType, atomPos=atomPos, ngrid=ngrid, dcell=dcell, iOutBuff=iOutBuff, iMO0=0, iMO1=iMO1, i0orb=i0orb, bDownalod=True )
    #ocl.saveToXsfAtoms    ( "dens_GPU_.xsf", iOutBuff, atomType, atomPos )
    Q_gpu = np.sum(data  )*dvol
    Q_cpu = np.sum(ewfaux)*dvol
    print( "Qtot gpu ",Q_gpu, "cpu", Q_cpu )
    #error = ewfaux - data
    #Ns = (ngrid[0],ngrid[1],ngrid[2])
    #ocl.setGridShape_dCell( Ns, dCell, pos0=[0.0,0.0,0.0] )
    #pos0 = ocl.getCellHalf( Ns, dCell );   print( "!!!!!!!! pos0 ", pos0  )
    #atomPos[:,0]-=pos0[2]
    #atomPos[:,1]-=pos0[1]
    #atomPos[:,2]-=pos0[0]
    ocl.saveToXsfAtomsData( "dens_GPU.xsf",  data    , atomType, atomPos )
    ocl.saveToXsfAtomsData( "dens_CPU.xsf", ewfaux, atomType, atomPos )
    error = ewfaux - data
    ocl.saveToXsfAtomsData( "dens_err.xsf", error , atomType, atomPos )
    print( "DEBUG check_density_projection() DONE" )

def projectDens( iMO0=1, iMO1=None, atomType=None, atomPos=None, ngrid=(64,64,64), dcell = [0.2,0.2,0.2,0.2], p0=None, iOutBuff=0, Rcuts=[4.5,4.5], bSCF=False, bSaveXsf=False, bSaveBin=False, saveName="dens", bDen0diff=False ):
    print("# ========= projectDens() bSCF ", bSCF )
    (wfcoef,i0orb),_ = density_from_firecore( atomType=atomType, atomPos=atomPos, bSCF=bSCF, bGetGrid=False, bGetCoefs=True )
    project_dens_GPU( wfcoef, atomType=atomType, atomPos=atomPos, ngrid=ngrid, dcell=dcell, iOutBuff=iOutBuff, iMO0=iMO0, iMO1=iMO1, i0orb=i0orb, bDen0diff=bDen0diff )
    if bSaveXsf:
        ocl.saveToXsfAtoms( saveName+".xsf", iOutBuff,    atomType, atomPos  )
    if bSaveBin:
        ocl.saveToBin(      saveName+".bin", iOutBuff )

def projectDens0_new( atomType=None, atomPos=None, ngrid=(64,64,64), dcell=[0.2,0.2,0.2,0.2], iOutBuff=0, Rcuts=[4.5,4.5], bSaveXsf=False, bSaveBin=False, saveName="dens0", acumCoef=[1.0,-1.0]  ):
    print("# ========= projectDens0_new() " )
    sys.path.append("../../")
    import pyBall as pb
    elems, dct, ords, Rcuts = iZs2dict(atomType)
    #print("# ======== FireCore Run " )
    #print ("atomType ", atomType)
    #print ("atomPos  ", atomPos)
    #print ("ords  ", ords)
    dCell = np.array([[dcell[0],0.0,0.0],[0.0,dcell[1],0.0],[0.0,0.0,dcell[2]]],dtype=np.float32)
    Ns = (ngrid[0],ngrid[1],ngrid[2])    
    ocl.tryInitFFT( Ns )
    ocl.tryLoadWfBasis( elems, Rcuts=Rcuts )
    ocl.setGridShape_dCell( Ns, dCell )
    #ocl.setTypes( [1,4], [[1.0,0.0],[1.0,3.0]] )
    ocl.setTypesZ( sorted(list(set(atomType))) )
    ords=np.array( ords, dtype=np.int32)
    ocl.projectAtomsDens0( iOutBuff, apos=atomPos, atypes=ords, acumCoef=acumCoef) 
    if bSaveXsf:
        ocl.saveToXsfAtoms( saveName+".xsf", iOutBuff, atomType, atomPos  )
    if bSaveBin:
        ocl.saveToBin( saveName+".bin", iOutBuff )

if __name__ == "__main__":
    #Test_projectWf( iMO=2 )
    #Test_projectAtomPosTex2()
    Test_projectDens( )



    #exit()
    #Test_projectAtoms(n=64)
    