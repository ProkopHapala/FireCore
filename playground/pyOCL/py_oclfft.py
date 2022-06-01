#!/usr/bin/python

import numpy as np
from   ctypes import c_int, c_double, c_long, c_bool, c_float, c_char_p, c_bool, c_void_p, c_char_p
import ctypes as ct
import os
import sys
#import pyopencl as cl

c_float_p  = ct.POINTER(c_float)
c_double_p = ct.POINTER(c_double)
c_int_p    = ct.POINTER(c_int)
c_long_p   = ct.POINTER(c_long)

def fortReal(a):
    ct.byref(ct.c_double(a))

def _np_as(arr,atype):
    if arr is None:
        return None
    else: 
        return arr.ctypes.data_as(atype)

const_Bohr_Radius = 0.529177210903
pref_s = 0.28209479177 # sqrt(1.0f/(4.0f*M_PI));
pref_p = 0.4886025119 # sqrt(3.0f/(4.0f*M_PI));
pref_d = 1.09254843059 # sqrt(15.0f/(4.0f*M_PI));

array1ui = np.ctypeslib.ndpointer(dtype=np.uint32, ndim=1, flags='CONTIGUOUS')
array1i  = np.ctypeslib.ndpointer(dtype=np.int32,  ndim=1, flags='CONTIGUOUS')
array1l  = np.ctypeslib.ndpointer(dtype=np.int64,  ndim=1, flags='CONTIGUOUS')
array2i  = np.ctypeslib.ndpointer(dtype=np.int32,  ndim=2, flags='CONTIGUOUS')
array1d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS')
array2d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=2, flags='CONTIGUOUS')
array3d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=3, flags='CONTIGUOUS')

'''
    void init(){
    int  upload   ( int i, const float* cpu_data )
    int  download ( int i,       float* cpu_data )
    void initFFT  ( int ndim, size_t* Ns_ ){
    void run_fft( int ibuff, bool fwd, float* data )
    void convolve ( int ibuffA, int ibuffB, int ibuff_result )
    void projectAtoms( float4* atoms, float4* coefs, int ibuff_result )
    void projectAtomPosTex( float* atoms, float* coefs, int nPos, float* poss, float* out )
    void setGridShape( float* pos0, float* dA, float* dB, float* dC ){
    int  initBasisTable( int nx, int ny, float* data );
    void approx( int npoints, int npolys, double* xs, double* ys, double* ws ){
    void loadWf(const char* fname, double* out){
    void loadWfBasis( double Rcut, int nsamp, int ntmp, int nZ, int* iZs ){
    void saveToXsf(const char* fname, int ibuff){
    void initFireBall( int natoms, int* atypes, double* apos ){
'''

#mode=ct.RTLD_GLOBAL
mode=ct.RTLD_LOCAL

#lib1 = ct.CDLL( name, ct.RTLD_GLOBAL )
lib2 = ct.CDLL(  "/usr/lib/x86_64-linux-gnu/libOpenCL.so", mode )
#lib2 = ct.CDLL(  "/usr/local/lib64/libclFFT.so", mode )
lib  = ct.CDLL(  "./libOCLfft.so",               mode )

lib.init.argtypes  = [ ] 
lib.init.restype   =  None
def init():
    lib.init( )

lib.cleanup.argtypes  = [ ] 
lib.cleanup.restype   =  None
def cleanup():
    lib.cleanup( )

#void loadData( float* data_ );
lib.upload.argtypes  = [ c_int, c_float_p ] 
lib.upload.restype   =  c_int
def loadData( i, data):
    return lib.loadData( i, _np_as( data, c_float_p ) )

#void loadData( float* data_ );
lib.download.argtypes  = [ c_int, c_float_p ] 
lib.download.restype   =  c_int
def download( i, data):
    return lib.download( i, _np_as( data, c_float_p ) )

#void loadData( float* data_ );
lib.upload.argtypes  = [ c_int, c_float_p ] 
lib.upload.restype   =  None
def upload( i, data):
    lib.upload( i, _np_as( data, c_float_p ) )

#void initFFT ( int ndim, int* Ns );
lib.initFFT.argtypes  = [ c_int, array1l ] 
lib.initFFT.restype   =  None
def initFFT( Ns ):
    Ns=np.array(Ns,dtype=np.int64)
    ndim=len(Ns)
    lib.initFFT( ndim, Ns )
    print( "DEBUG py: initFFT DONE" )

#void initFFT ( int ndim, int* Ns );
lib.initAtoms.argtypes  = [ c_int ] 
lib.initAtoms.restype   =  None
def initAtoms( N ):
    lib.initAtoms( N )

#int  initBasisTable( int nx, int ny, float* data );
lib.initBasisTable.argtypes  = [ c_int, c_int, c_float_p ] 
lib.initBasisTable.restype   =  c_int
def initBasisTable( nx,ny, data ):
    return lib.initBasisTable( nx,ny, _np_as( data, c_float_p ) )

lib.convolve.argtypes  = [ c_int, c_int, c_int ] 
lib.convolve.restype   =  None
def convolve(iA,iB,iOut):
    lib.convolve( iA,iB,iOut )

#projectAtoms( float* atoms, float4* coefs, int ibuff_result )
lib.projectAtoms.argtypes  = [ c_float_p, c_float_p, c_int ] 
lib.projectAtoms.restype   =  None
def projectAtoms(atoms,coefs,iOut):
    lib.projectAtoms( _np_as( atoms, c_float_p ),_np_as( coefs, c_float_p ),iOut )


#void projectAtomPosTex( float* atoms, float* coefs, int nPos, float* poss, float* out )
lib.projectAtomPosTex.argtypes  = [ c_float_p, c_float_p, c_int, c_float_p, c_float_p ] 
lib.projectAtomPosTex.restype   =  None
def projectAtomPosTex(atoms,coefs,poss,out=None):
    nPos=len(poss)
    if(out is None):
        out = np.zeros((nPos,2), dtype=np.float32)
    lib.projectAtomPosTex( _np_as( atoms, c_float_p ),_np_as( coefs, c_float_p ),nPos, _np_as( poss, c_float_p ), _np_as( out, c_float_p ) )
    return out

#setGridShape( float* pos0, float* dA, float* dB, float* dC ){
lib.setGridShape.argtypes  = [ c_float_p, c_float_p, c_float_p, c_float_p ] 
lib.setGridShape.restype   =  None
def setGridShape( pos0=[0.,0.,0.,0.],dA=[0.1,0.,0.,0.],dB=[0.,0.1,0.,0.],dC=[0.,0.,0.1,0.0]):
    pos0=np.array(pos0,dtype=np.float32)
    dA=np.array(dA,dtype=np.float32)
    dB=np.array(dB,dtype=np.float32)
    dC=np.array(dC,dtype=np.float32)
    print( "dA ", dA);
    print( "dB ", dB);
    print( "dC ", dC);
    lib.setGridShape( _np_as( pos0, c_float_p ),_np_as( dA, c_float_p ), _np_as( dB, c_float_p ), _np_as( dC, c_float_p ) )

def setGridShape_dCell( Ns, dCell ):
    pos0=[  dCell[0,0]*(-Ns[0]//2),  dCell[1,1]*(-Ns[1]//2),  dCell[2,2]*(-Ns[1]//2),  0.0 ]
    print( " setGridShape_dCell() pos0 ", pos0 )
    setGridShape( 
        pos0=pos0,
        dA=dCell[0]+[0.0],
        dB=dCell[1]+[0.0],
        dC=dCell[2]+[0.0]
    )


'''
# void run_fft( int ibuff, bool fwd, float* data )
lib.runfft.argtypes  = [ c_int, c_bool, c_float_p ] 
lib.runfft.restype   =  None
def runfft(ibuff, outbuff, fwd=True ):
    lib.runfft( ibuff, fwd,  _np_as( outbuff, c_float_p ) )
'''


#void initFireBall( int natoms, int* atypes, double* apos ){
lib.initFireBall.argtypes  = [ c_int, array1i, array2d ] 
lib.initFireBall.restype   =  None
def initFireBall( atypes, apos ):
    natoms = len(atypes)
    lib.initFireBall( natoms, atypes, apos )

# void run_fft( int ibuff, bool fwd, float* data )
lib.runfft.argtypes  = [ c_int, c_bool ] 
lib.runfft.restype   =  None
def runfft(ibuff, fwd=True ):
    lib.runfft( ibuff, fwd )

#void approx( int npoints, int npolys, double* xs, double* ys, double* ws ){
lib.approx.argtypes  = [ c_int, c_int, c_double_p, c_double_p, c_double_p ] 
lib.approx.restype   =  None
def approx( xs, ys, ws=None, npoly=14 ):
    n=len(xs)
    if ws is None:
        ws = np.ones(n)
    return lib.approx( n, npoly, _np_as( xs, c_double_p ),_np_as( ys, c_double_p ), _np_as( ws, c_double_p ) )

# loadWf(const char* fname, double* out){
lib.loadWf.argtypes  = [ c_char_p, c_float_p ] 
lib.loadWf.restype   =  None
def loadWf_C( fname, n=1000 ):
    data=np.zeros(n, dtype=np.float32)
    #fname = c_char_p( fname) 
    fname = fname.encode('utf-8')
    lib.loadWf( fname, _np_as( data, c_float_p ) )
    #print( data )
    return data

# void saveToXsf(const char* fname, int ibuff){
lib.saveToXsf.argtypes  = [ c_char_p, c_int ] 
lib.saveToXsf.restype   =  None
def saveToXsf( fname, ibuff ):
    fname = fname.encode('utf-8')
    lib.saveToXsf( fname, ibuff )

#void saveToXsfAtoms(const char* fname, int ibuff, int natoms, int* atypes, double* apos )
lib.saveToXsfAtoms.argtypes  = [ c_char_p, c_int,  c_int,c_int_p,c_double_p ] 
lib.saveToXsfAtoms.restype   =  None
def saveToXsfAtoms( fname, ibuff, atypes, apos ):
    fname = fname.encode('utf-8')
    na = len(atypes)
    lib.saveToXsfAtoms( fname, ibuff,   na, _np_as( atypes, c_int_p ), _np_as( apos, c_double_p ) )

#    void loadWfBasis( const char* path, float RcutSamp, int nsamp, int ntmp, int nZ, int* iZs, float* Rcuts ){
lib.loadWfBasis.argtypes  = [  c_char_p, c_float, c_int, c_int, c_int, c_int_p, c_float_p ] 
lib.loadWfBasis.restype   =  None
def loadWfBasis( iZs, nsamp=100, ntmp=1000, RcutSamp=5.0, path="Fdata/basis/", Rcuts=None, RcutDef=4.5 ):
    # NOTE : RcutSamp is in [Angstroem] while Rcuts and RcutDef are in [bohr_radius];    RcutSamp should not be changed without chaning "wf_tiles_per_angstroem" in myprog.cl

    nZ=len(iZs)
    iZs=np.array(iZs,dtype=np.int32)
    path = path.encode('utf-8')
    if Rcuts is None:
        Rcuts=np.ones(nZ,dtype=np.float32)*RcutDef
    else:
        Rcuts=np.array(Rcuts,dtype=np.float32)
    return lib.loadWfBasis( path, RcutSamp, nsamp, ntmp, nZ, _np_as( iZs, c_int_p ), _np_as( Rcuts, c_float_p ) )


def projectAtoms__( atoms, acoefs):
    import matplotlib.pyplot as plt
    import time
    print( "# ========= TEST   runFFT()  " )

    #initFireBall( atypes, atoms )

    init()
    Ns=(100,100,100)
    initFFT( Ns  )
    loadWfBasis( [1,6], Rcuts=[4.5,4.5] )
    initAtoms( len(atoms) )
    #initBasisTable( basis.shape[0], basis.shape[1], basis )

    setGridShape( )
    t0 = time.clock()
    projectAtoms( atoms, acoefs, 0 )
    saveToXsf( "test.xsf", 0 )
    arrA = np.zeros(Ns,dtype=np.csingle  )
    download ( 0, arrA )    
    t = time.clock()-t0; print( "projection time ", t )
    plt.figure(); 
    #print( arrA[10,10].real )
    plt.imshow( arrA[10].real ) 
    #plt.imshow( np.log( np.abs(arrA[10])) ) 
    plt.grid()
    plt.show(); 
    cleanup()










def Test_projectAtoms(n=64, natoms=1000):
    import matplotlib.pyplot as plt
    import time
    print( "# ========= TEST   runFFT()  " )
    
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
    init()
    #Ns=(128,128,128)
    Ns=(100,100,100)
    initFFT( Ns  )
    loadWfBasis( [1,6], Rcuts=[4.5,4.5] )
    initAtoms( len(atoms) )
    #initBasisTable( basis.shape[0], basis.shape[1], basis )

    setGridShape()
    t0 = time.clock()
    projectAtoms( atoms, coefs, 0 )
    saveToXsf( "test.xsf", 0 )
    arrA = np.zeros(Ns,dtype=np.csingle  )
    download ( 0, arrA )    
    t = time.clock()-t0; print( "projection time ", t )
    plt.figure(); 
    #print( arrA[10,10].real )
    plt.imshow( arrA[10].real ) 
    #plt.imshow( np.log( np.abs(arrA[10])) ) 
    plt.grid()
    plt.show(); 
    cleanup()



def Test_projectAtomPosTex():
    import matplotlib.pyplot as plt
    import time
    print( "# ========= Test_projectAtomPosTex  " )
    
    xref = np.linspace(0.0,4.50*const_Bohr_Radius,1000)
    yref = loadWf_C( "Fdata/basis/001_450.wf1", n=1000 )
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
    apos_,wfcoef_ = convCoefs( atomType, wfcoef, atomPos )
    # ========== GPU
    init()                                   #;print("DEBUG py 2")
    loadWfBasis( [1,6], Rcuts=[4.5,4.5] )         #;print("DEBUG py 3") 
    initAtoms( len(apos_) )                  #;print("DEBUG py 4")
    #initBasisTable( basis.shape[0], basis.shape[1], basis )
    out = projectAtomPosTex(apos_,wfcoef_,poss_) #;print("DEBUG py 5")
    plt.plot( xref, yref*pref_s,    label="load" ) 
    plt.plot( poss_[:,0], out[:,0], label="GPU" ) 
    plt.plot( poss [:,0], ys,       label="CPU", ls="--" ) 
    #plt.imshow( np.log( np.abs(arrA[10])) ) 
    plt.legend()
    plt.grid()
    plt.show(); 
    cleanup()



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
    init()
    initFFT(    arrA.shape )
    upload ( 0, arrA )    ;plt.figure(); plt.imshow( arrA.real )
    upload ( 1, arrB )    ;plt.figure(); plt.imshow( arrB.real )
    #convolve( 0,1,   2 )
    #arrC = np.zeros(arrA.shape)
    #download( 2, arrC)    ;plt.figure(); plt.imshow( arrC.real ) 
    runfft (0 ); download ( 0, arrA )    ;plt.figure(); plt.imshow( np.log( np.abs(arrA)) ) 
    runfft (1 ); download ( 1, arrB )    ;plt.figure(); plt.imshow( np.log( np.abs(arrB)) ) 
    plt.show(); 
    cleanup()

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
    init()
    initFFT(    arrA.shape )
    upload ( 0, arrA )    ;plt.figure(); plt.imshow( arrA.real )
    upload ( 1, arrB )    ;plt.figure(); plt.imshow( arrB.real )
    upload ( 2, arrC ) 
    convolve( 0,1,2 )
    download( 2, arrC)    ;plt.figure(); plt.imshow( arrC.real )
    plt.show(); 
    cleanup()

def countOrbs(atypes):
    norbs = [ 1 if (x == 1) else 4 for x in atypes ]
    return np.cumsum(norbs)

def convCoefs( atypes, oCs, oatoms, typeDict ):
    na = len( atypes)
    norbs = [ 1 if (x == 1) else 4 for x in atypes ]
    norb = sum(norbs)
    atoms = np.zeros( (na,4), dtype=np.float32)
    coefs = np.zeros( (na,4), dtype=np.float32)
    #print( "atoms.shape ", atoms.shape, "oatoms.shape ", oatoms.shape, " oCs.shape ", oCs.shape  )
    atoms[:,:3] = oatoms[:,:3]
    atoms[:,3]=0.1
    io=0
    for ia,no in enumerate(norbs):
        coefs[ia,3]=oCs[io]; io+=1
        if(no>1):
            coefs[ia,0]=oCs[io+0]
            coefs[ia,1]=oCs[io+1]
            coefs[ia,2]=oCs[io+2]
            io+=3
        atoms[ia,3] += typeDict[atypes[ia]] + 0.1
    #print( "atoms ", atoms )
    #print( "coefs ", coefs )
    print( " atomTypes CPU ", atypes[ia]  )
    print( " atomTypes GPU ", atoms[:,3]  )
    return atoms, coefs

def Test_projectWf( ):
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
    ngrid, dCell, lvs = fc.setupGrid()
    ewfaux = fc.getGridMO( 1,ngrid=ngrid)   ;print( "ewfaux.min(),ewfaux.max() ", ewfaux.min(),ewfaux.max() )
    sh = ewfaux.shape                       ;print( "ewfaux.shape ", sh )
    fc.orb2xsf(1); #exit()

    i0orb  = countOrbs( atomType )           ;print("i0orb ", i0orb)  
    wfcoef = fc.get_wfcoef(norb=i0orb[-1])

    print("# ========= PyOCL Wave-Function Projection " )
    iMO=0
    #wfcoef = wfcoef[1,:]
    print( "wfcoef: \n", wfcoef )
    apos_,wfcoef_ = convCoefs( atomType, wfcoef[iMO,:], atomPos )
    #wfcoef_ = np.array( [ 0.0,0.0,0.0,1.0,   0.0,0.0,0.0,1.0,   0.0,0.0,0.0,1.0,   0.0,0.0,0.0,1.0,   0.0,0.0,0.0,1.0 ], dtype=np.float32)
    init()                           
    Ns = (ngrid[0],ngrid[1],ngrid[2])
    initFFT( Ns  )                  
    loadWfBasis( [1,6], Rcuts=[4.5,4.5] )    
    initAtoms( len(apos_) )          
    setGridShape_dCell( Ns, dCell )
    projectAtoms( apos_, wfcoef_, 0 ) 
    #saveToXsf( "test.xsf", 0 )
    saveToXsfAtoms( "test.xsf", 0,    atomType, atomPos  )
    arrA = np.zeros(Ns+(2,),dtype=np.float32)  
    download(0,arrA)                  

    print("# ========= Plot GPU vs CPU comparison " )
    #print( arrA  [16,16,:,0] )
    #print( ewfaux[16,16,:] )
    import matplotlib.pyplot as plt
    ix0=ngrid[0]//2
    iy0=ngrid[1]//2
    plt.plot( arrA  [ix0,iy0,:,0], label="GPU" )
    plt.plot( ewfaux[ix0,iy0,:  ], label="CPU" )
    print( "arrA  .shape ", arrA  .shape )
    print( "ewfaux.shape ", ewfaux.shape )
    #plt.figure(); plt.imshow( arrA  [ngrid[0]//2 ,:,:,0] ); plt.title("GPU")
    #plt.figure(); plt.imshow( ewfaux[ngrid[0]//2 ,:,:  ] ); plt.title("CPU")
    plt.legend()
    plt.show()






def Test_projectAtomPosTex2():
    sys.path.append("../../")
    import pyBall as pb
    from pyBall import FireCore as fc
    #atomType = np.array( [ a[0] for a in acs ], dtype=np.int32 )
    #atomPos  = np.array( [ a[1] for a in acs ] )
    #apos     = np.array( [ a[1] for a in acs ], dtype=np.float32 )
    #coefs    = np.array( [ a[2] for a in acs ], dtype=np.float32) 
    # ----- Sampling Points
    npos = 100   # must be multiple of local group size
    xs = np.linspace(0.0,5.0,npos)
    poss  = np.zeros((npos,3)                 ); poss [:,0]= xs
    poss_ = np.zeros((npos,4),dtype=np.float32); poss_[:,0]= xs

    # ----- Make CH4 molecule
    
    atomType = np.array([6,1,1,1,1]).astype(np.int32)     # CH4
    #atomType = np.array([1,1,1,1,1]).astype(np.int32)    # HH4
    atomPos  = np.array([
        [ 0.01,     0.0,     0.0],
        [-1.0,     +1.0,    -1.0],
        [+1.0,     -1.0,    -1.0],
        [-1.0,     -1.0,    +1.0],
        [+1.0,     +1.0,    +1.0],
    ])
    
    '''
    # ----- Make H-atom
    atomType = np.array([1]).astype(np.int32)
    atomPos  = np.array([[ 0.01,     0.0,     0.0],])
    '''
    
    print("# ======== FireCore Run " )
    # ============= CPU Fortran
    fc.preinit()
    fc.init( atomType, atomPos )
    # =========== Electron Density
    fc.assembleH( atomPos)
    fc.solveH()
    sigma=fc.updateCharges() ; print( sigma )
    ngrid, dCell, lvs = fc.setupGrid()

    wfcoef = np.array( [ 1.0,0.0,0.0,0.0,  0.0,0.0,0.0,0.0 ] )  # CH4
    #wfcoef = np.array( [ 1.0,  0.0,0.0,0.0,0.0 ] )             # HH4
    fc.set_wfcoef(wfcoef,iMO=1,ikp=1)

    ys = fc.orb2points( poss, ys=None, iMO=1,  ikpoint=1 )     # ;print( " ys ", ys) 

    print("# ========= PyOCL Wave-Function Projection " )
    #i0orb  = countOrbs( atomType ) 
    #wfcoef = fc.get_wfcoef( norb=i0orb[-1] )
    typeDict={ 1:0, 6:1 }
    apos_,wfcoef_ = convCoefs( atomType, wfcoef, atomPos, typeDict )
    init()                                   #;print("DEBUG py 2")
    loadWfBasis( [1,6], Rcuts=[4.5,4.5] )    #;print("DEBUG py 3") 
    initAtoms( len(apos_) )                  #;print("DEBUG py 4")
    out = projectAtomPosTex(apos_,wfcoef_,poss_) #;print("DEBUG py 5")

    print("# ========= Plot GPU vs CPU comparison " )
    import matplotlib.pyplot as plt
    #plt.plot( xref, yref*pref_s,    label="load" ) 
    plt.plot( poss_[:,0], out[:,0], label="GPU" ) 
    plt.plot( poss [:,0], ys,       label="CPU", ls="--" ) 
    plt.legend()
    plt.grid()
    plt.show(); 
    cleanup()






if __name__ == "__main__":
    #Test_projectWf( )
    Test_projectAtomPosTex2()




    #exit()
    #Test_projectAtoms(n=64)
    