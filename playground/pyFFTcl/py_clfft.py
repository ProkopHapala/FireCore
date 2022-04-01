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
void initOCL();
void initFFT ( int ndim, int* Ns );
void loadData( float* data_ );
void runFFT();
void cleanup();
void runAll( );
'''

#mode=ct.RTLD_GLOBAL
mode=ct.RTLD_LOCAL

#lib1 = ct.CDLL( name, ct.RTLD_GLOBAL )
lib2 = ct.CDLL(  "/usr/lib/x86_64-linux-gnu/libOpenCL.so", mode )
#lib2 = ct.CDLL(  "/usr/local/lib64/libclFFT.so", mode )
lib  = ct.CDLL(  "./libfft3d.so",               mode )

lib.initOCL.argtypes  = [ ] 
lib.initOCL.restype   =  None
def initOCL():
    lib.initOCL( )

lib.runFFT.argtypes  = [ ] 
lib.runFFT.restype   =  None
def runFFT():
    lib.runFFT( )

lib.cleanup.argtypes  = [ ] 
lib.cleanup.restype   =  None
def cleanup():
    lib.cleanup( )

lib.runAll.argtypes  = [ ] 
lib.runAll.restype   =  None
def runAll():
    lib.runAll( )

#void loadData( float* data_ );
lib.loadData.argtypes  = [ c_float_p ] 
lib.loadData.restype   =  None
def loadData( data):
    lib.loadData( _np_as( data, c_float_p ) )

#void initFFT ( int ndim, int* Ns );
lib.initFFT.argtypes  = [ c_int, array1l ] 
lib.initFFT.restype   =  None
def initFFT( arr ):
    ndim=len(arr.shape)
    Ns = np.array(arr.shape,dtype=np.int64)
    lib.initFFT( ndim, Ns )


# ======================= RUN ===============

#runAll()
#initOCL()

'''
print(  "============ RUN FFT 2D ==========" )
xs       = np.linspace(-2.0,2.0,4)
Xs,Ys,Zs = np.meshgrid(xs,xs,xs)
arr      = np.sin(Xs*3)*np.cos(Ys*20)*np.cos(Zs*20)/( 1 + Xs**2 + Ys**2 + Zs**2)
arr = arr.astype( np.csingle )
print( arr )
initOCL()         ;print( "DEBUG 1")          
initFFT (arr)     ;print( "DEBUG 2")                  
loadData(arr)     ;print( "DEBUG 3")                 
runFFT()          ;print( "DEBUG 4")   
print( arr )
'''

print(  "============ RUN FFT 3D ==========" )
n=1024
xs    = np.linspace(-2.0,2.0,n)
Xs,Ys = np.meshgrid(xs,xs)
arr   = np.sin(Xs*3 + (Xs**2)*5.1 )*np.cos(Ys*20 + (Xs**2)*3.0)/(1+Xs**2+Ys**2)   + np.random.rand(n,n)*0.1
arr   = arr.astype( np.csingle )
#arr.imag = np.cos(Xs*3)*np.sin(Ys*20)/(1+Xs**2+Ys**2)

import matplotlib.pyplot as plt
plt.figure(); plt.imshow( arr.real )
#plt.figure(); plt.imshow( arr.imag )


initOCL()         ;print( "DEBUG 1")          
initFFT(arr)      ;print( "DEBUG 2")                  
loadData(arr)     ;print( "DEBUG 3")                 
runFFT()          ;print( "DEBUG 4")   

print( arr )
#plt.figure(); plt.imshow( arr.real ); plt.colorbar(); 
plt.figure(); plt.imshow( arr.real, vmin=-100.0, vmax=100.0 ); plt.colorbar(); 
plt.show()

cleanup()         ;print( "DEBUG 5")   