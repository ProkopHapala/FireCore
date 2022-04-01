#!/usr/bin/python

import numpy as np
from   ctypes import c_int, c_double, c_bool, c_float, c_char_p, c_bool, c_void_p, c_char_p
import ctypes as ct
import os
import sys
#import pyopencl as cl

'''
void makeData();
void printData( float* X );
void printDataX();
void initOCL();
void initFFT();
void planFFT();
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


lib.runAll.argtypes  = [ ] 
lib.runAll.restype   =  None
def runAll():
    lib.runAll( )

#lib.initOCL.argtypes  = [ ] 
#lib.initOCL.restype   =  None
#def initOCL():
#    lib.initOCL( )

runAll()
#initOCL()

print( lib2 )