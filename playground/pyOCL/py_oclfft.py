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

#void loadWfBasis( double Rcut, int nsamp, int ntmp, int nZ, int* iZs ){
lib.loadWfBasis.argtypes  = [  c_char_p, c_float, c_int, c_int, c_int, c_int_p ] 
lib.loadWfBasis.restype   =  None
def loadWfBasis( iZs, nsamp=100, ntmp=1000, Rcut=4.5, path="Fdata/basis/" ):
    nZ=len(iZs)
    iZs=np.array(iZs,dtype=np.int32)
    path = path.encode('utf-8')
    return lib.loadWfBasis( path, Rcut, nsamp, ntmp, nZ, _np_as( iZs, c_int_p ) )


def projectAtoms__( atoms, acoefs):
    import matplotlib.pyplot as plt
    import time
    print( "# ========= TEST   runFFT()  " )

    #initFireBall( atypes, atoms )

    init()
    Ns=(100,100,100)
    initFFT( Ns  )
    loadWfBasis( [1,6], Rcut=4.5 )
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
    loadWfBasis( [1,6], Rcut=4.5 )
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

'''
n=256
xs       = np.linspace(-2.0,2.0,n)
Xs,Ys,Zs = np.meshgrid(xs,xs,xs)
arrA      = np.sin(Xs*3)*np.cos(Ys*20)*np.cos(Zs*20)
arrB      = 1/( 1 + Xs**2 + Ys**2 + Zs**2)
'''

#data = loadWf_C( "basis/001_450.wf1"  )
#print( data )
#exit()

'''

def loadWf( fname="basis/001_450.wf1" ):
    txt = open( fname, 'r' ).readlines()
    txt = [ l.replace('D','e') for l in txt ]
    data = np.genfromtxt( txt, skip_header=5, skip_footer=1 )
    data = data.flatten()
    print( len(data) )
    return data

import matplotlib.pyplot as plt
xs = np.linspace( 0.0,4.5,1000 )
wf_Hs = loadWf( "basis/001_450.wf1" )  ;plt.plot( xs, wf_Hs, label="Hs" )
wf_Cs = loadWf( "basis/006_450.wf1" )  ;plt.plot( xs, wf_Cs, label="Cs" )
wf_Cp = loadWf( "basis/006_450.wf2" )  ;plt.plot( xs, wf_Cp, label="Cp" )

approx( xs, wf_Hs, npoly=15 )
#approx( xs, wf_Cs, npoly=15 )
#approx( xs, wf_Cp, npoly=15 )

plt.grid()
plt.show()

'''

#Test_fft2d()
#Test_Convolution2d()


sys.path.append("../../")


def convCoefs( atypes, oCs, oatoms, iMO=0 ):
    na = len( atypes)
    norbs = [ 1 if (x == 1) else 4 for x in atypes ]
    norb = sum(norbs)
    atoms = np.zeros( (na,4), dtype=np.float32)
    coefs = np.zeros( (na,4), dtype=np.float32)
    print( "atoms.shape ", atoms.shape, "oatoms.shape ", oatoms.shape  )
    atoms[:,:3] = oatoms[:,:3]
    atoms[:,3]=0.1
    j=0
    for i,no in enumerate(norbs):
        coefs[i,3]=oCs[iMO,j]; j+=1
        if(no>1):
            coefs[i,0]=oCs[iMO,j+0]
            coefs[i,1]=oCs[iMO,j+1]
            coefs[i,2]=oCs[iMO,j+2]
            j+=3
    print( "atoms ", atoms )
    print( "coefs ", coefs )
    return atoms, coefs


if __name__ == "__main__":
    
    import pyBall as pb
    from pyBall import FireCore as fc

    #atomType = np.random.randint(6, size=natoms).astype(np.int32)
    #atomPos  = np.random.random((3,natoms))
    '''
    atomType = np.array([6,1,1,1,1]).astype(np.int32)
    atomPos  = np.array([
        [ 0.1,      0.0,     0.0],
        [-1.0,     +1.0,    -1.0],
        [+1.0,     -1.0,    -1.0],
        [-1.0,     -1.0,    +1.0],
        [+1.0,     +1.0,    +1.0],
    ])
    '''
    
    atomType = np.array([1]).astype(np.int32)
    atomPos  = np.array([
        [ 0.1,      0.0,     0.0],
    ])
    

    #initFireBall( atomType, atomPos )

    print ("atomType ", atomType)
    print ("atomPos  ", atomPos)
    fc.preinit()
    fc.init( atomType, atomPos )
    
    # =========== Electron Density
    fc.assembleH( atomPos)
    fc.solveH()
    sigma=fc.updateCharges() ; print( sigma )
    ngrid, dCell, lvs = fc.setupGrid()

    '''
    pref_s = np.sqrt( 1.0/(4.0*np.pi))
    pref_p = np.sqrt( 3.0/(4.0*np.pi))
    pref_d = np.sqrt(15.0/(4.0*np.pi))
    ys_2 = loadWf_C( "Fdata/basis/001_450.wf1" )*pref_s
    ys = fc.getpsi(in1=1, issh=1, n=250, dx=0.01, x0=0.0, ys=None, theta=0.0, phi=0.0)
    import matplotlib.pyplot as plt
    plt.plot( np.linspace(0,3.5,len(ys)),   ys,   label="FireCore" )
    plt.plot( np.linspace(0,3.5,len(ys_2)), ys_2, label="pyOCL" )
    plt.legend()
    plt.show()
    exit(0)
    '''

    ewfaux = fc.getGridMO( 1,ngrid=ngrid)
    #ewfaux = fc.getGridDens( ngrid=ngrid )
    print( ewfaux.min(),ewfaux.max() )
    
    sh = ewfaux.shape
    print( "ewfaux.shape ", ewfaux.shape )
    #plt.figure(); plt.imshow( ewfaux[ sh[0]//2+5,:,: ] )
    #plt.figure(); plt.imshow( ewfaux[ sh[0]//2  ,:,: ] )
    #plt.figure(); plt.imshow( ewfaux[ sh[0]//2-5,:,: ] )
    #plt.show()


    wfcoef = fc.firecore_get_wfcoef(norb=1)
    print( "wfcoef: \n", wfcoef )
    apos_,wfcoef_ = convCoefs( atomType, wfcoef, atomPos )

    init()                            ;print("DEBUG 0 ")
    Ns=(32,32,32)
    initFFT( Ns  )                    ;print("DEBUG 1 ")
    loadWfBasis( [1,6], Rcut=4.5 )    ;print("DEBUG 2 ")
    initAtoms( len(apos_) )           ;print("DEBUG 3 ")
    setGridShape( pos0=[-2.5812965+0.16133103125*1.5,-2.5812965+0.16133103125*1.5*0,-2.5812965+0.16133103125*1.5*0,0.],dA=[0.16133103125,0.,0.,0.],dB=[0.,0.16133103125,0.,0.],dC=[0.,0.,0.16133103125,0.0] )
    projectAtoms( apos_, wfcoef_, 0 ) ;print("DEBUG 4 ") 
    #saveToXsf( "test.xsf", 0 )
    arrA = np.zeros(Ns+(2,),dtype=np.float32)  
    download(0,arrA)                  ;print("DEBUG 5 ")

    print( arrA  [16,16,:,0] )
    print( ewfaux[16,16,:] )

    import matplotlib.pyplot as plt
    plt.plot( arrA  [16,16,:,0], label="GPU" )
    plt.plot( ewfaux[16,16,:],   label="CPU" )
    #plt.figure(); plt.imshow( arrA  [16,:,:,0] )
    #plt.figure(); plt.imshow( ewfaux[16,:,:  ] )
    plt.legend()
    plt.show()



    #exit()
    #Test_projectAtoms(n=64)
    