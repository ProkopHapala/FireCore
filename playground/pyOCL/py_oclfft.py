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
def setGridShape( pos0=[0.,0.,0.,0.],dA=[0.1,0.,0.,0.],dB=[0.,0.1,0.,0.],dC=[0.,0.,0.,0.1]):
    pos0=np.array(pos0,dtype=np.float32)
    dA=np.array(dA,dtype=np.float32)
    dB=np.array(dB,dtype=np.float32)
    dC=np.array(dC,dtype=np.float32)
    lib.setGridShape( _np_as( pos0, c_float_p ),_np_as( dA, c_float_p ), _np_as( dB, c_float_p ), _np_as( dC, c_float_p ) )


'''
# void run_fft( int ibuff, bool fwd, float* data )
lib.runfft.argtypes  = [ c_int, c_bool, c_float_p ] 
lib.runfft.restype   =  None
def runfft(ibuff, outbuff, fwd=True ):
    lib.runfft( ibuff, fwd,  _np_as( outbuff, c_float_p ) )
'''

# void run_fft( int ibuff, bool fwd, float* data )
lib.runfft.argtypes  = [ c_int, c_bool ] 
lib.runfft.restype   =  None
def runfft(ibuff, fwd=True ):
    lib.runfft( ibuff, fwd )


def Test_projectAtoms(n=64, natoms=1000):
    import matplotlib.pyplot as plt
    import time
    print( "# ========= TEST   runFFT()  " )
    
    acs=[
    [[1.0,2.0,0.0,1.0],    [0.0,0.0,0.0,1.0]],  
    #[[2.0,3.0,0.0,1.0],    [0.0,0.0,0.0,1.0]],
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

    xs = np.linspace(0.0,10.0,100)*-1.0
    ys = np.ones(16)
    basis,_=np.meshgrid(xs,ys); #basis=np.exp(-0.3*(basis-2.0)**2)*-1+np.exp(-0.3*(basis)**2)*1;  
    basis=basis.astype(np.float32)
    plt.imshow(basis); plt.show()

    #print( "atoms ", atoms, atoms.dtype )
    #print( "coefs ", coefs )
    init()
    #Ns=(128,128,128)
    Ns=(100,100,100)
    initFFT( Ns  )
    initAtoms( len(atoms) )

    initBasisTable( basis.shape[0], basis.shape[1], basis )
    setGridShape()
    t0 = time.clock()
    projectAtoms( atoms, coefs, 0 )
    arrA = np.zeros(Ns,dtype=np.csingle  )
    download ( 0, arrA )    
    t = time.clock()-t0; print( "projection time ", t )
    plt.figure(); 
    #print( arrA[10,10].real )
    plt.imshow( arrA[10].real ) 
    #plt.imshow( np.log( np.abs(arrA[10])) ) 
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

#Test_fft2d()
#Test_Convolution2d()

Test_projectAtoms(n=64)