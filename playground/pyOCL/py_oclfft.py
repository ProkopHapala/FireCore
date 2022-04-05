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
    void convolve ( int ibuffA, int ibuffB, int ibuff_result ){
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
def initFFT( arr ):
    ndim=len(arr.shape)
    Ns = np.array(arr.shape,dtype=np.int64)
    lib.initFFT( ndim, Ns )


lib.convolve.argtypes  = [ c_int, c_int, c_int ] 
lib.convolve.restype   =  None
def convolve(iA,iB,iOut):
    lib.convolve( iA,iB,iOut )

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
    initFFT(    arrA )
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
    initFFT(    arrA )
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
Test_Convolution2d()
