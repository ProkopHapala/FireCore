#!/usr/bin/python

import numpy as np
from   ctypes import c_int, c_double, c_bool, c_float, c_char_p, c_bool, c_void_p, c_char_p
import ctypes as ct
import os
import sys
#import pyopencl as cl

c_float_p  = ct.POINTER(c_float)
c_double_p = ct.POINTER(c_double)
c_int_p    = ct.POINTER(c_int)

def fortReal(a):
    ct.byref(ct.c_double(a))

def _np_as(arr,atype):
    if arr is None:
        return None
    else: 
        return arr.ctypes.data_as(atype)


array1ui = np.ctypeslib.ndpointer(dtype=np.uint32, ndim=1, flags='CONTIGUOUS')
array1i  = np.ctypeslib.ndpointer(dtype=np.int32,  ndim=1, flags='CONTIGUOUS')
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
lib.initFFT.argtypes  = [ c_int, array1i] 
lib.initFFT.restype   =  None
def initFFT( arr ):
    ndim=len(arr.shape)
    Ns = np.array(arr.shape,dtype=np.int32)
    lib.initFFT( ndim, Ns )

runAll()
#initOCL()

print( lib2 )