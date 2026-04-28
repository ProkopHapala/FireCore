
import numpy as np
from   ctypes import c_int, c_double, c_bool, c_float, c_char_p, c_bool, c_void_p, c_char_p
import ctypes
import os
import sys

#sys.path.append('../')
#from pyMeta import cpp_utils 
from . import cpp_utils_ as cpp_utils
#import cpp_utils_ as cpp_utils


cpp_utils.BUILD_PATH = os.path.normpath( cpp_utils.PACKAGE_PATH + '../../cpp/Build/libs/Molecular' ) 
lib = cpp_utils.loadLib('Forces_lib', recompile=True)

c_double_p = ctypes.POINTER(c_double)
c_float_p  = ctypes.POINTER(c_float)
c_int_p    = ctypes.POINTER(c_int)
c_bool_p   = ctypes.POINTER(c_bool)

def _np_as(arr,atype):
    if arr is None:
        return None
    else: 
        return arr.ctypes.data_as(atype)

cpp_utils.s_numpy_data_as_call = "_np_as(%s,%s)"


# int init() {
lib.init.argtypes  = []
lib.init.restype   =  c_int
def init():
    return lib.init()

# void evalForce( char* funcName, int n, double* xs, Vec2d* FEs, double* params ){
lib.evalForce.argtypes  = [ c_char_p, c_int, c_double_p, c_double_p, c_double_p ]
lib.evalForce.restype   =  None
def evalForce( funcName, params, FEs=None, xs=None, x0=0.0, dx=0.1, n=100 ):
    if xs is not None: 
        n = len(xs)
    else:
        xs=np.linspace(x0,x0+n*dx,n, endpoint=False)
    if FEs is None: FEs=np.zeros((n,2))
    params = np.array( params )
    funcName = funcName.encode('utf-8')
    lib.evalForce( funcName, n, _np_as(xs,c_double_p), _np_as(FEs,c_double_p), _np_as(params,c_double_p) )
    return FEs, xs, funcName

def diff( xs, Es ):
    return -(Es[2:]-Es[:-2])/(xs[2:]-xs[:-2]), xs[1:-1]

def plotAll( funcName, params, x0=2.0, dx=0.01, n=1000 ):
    FEs, xs, funcName = evalForce( funcName, params, x0=x0, dx=dx, n=n )
    Fnum, xs_ = diff( xs, FEs[:,1] )
    plt.figure(figsize=(5,10))
    plt.subplot(211)
    plt.plot( xs,  FEs[:,0], '-',lw=0.5, label="E" ); Emin=FEs[:,0].min(); plt.ylim(Emin,-Emin)
    plt.subplot(212)
    plt.plot( xs,  -FEs[:,1], '-',lw=0.5, label="F" ); Fmin=FEs[:,1].min();
    plt.plot( xs_, -Fnum    , ':',lw=1.5, label="Fnum" ); plt.ylim(Fmin,-Fmin)
    plt.legend()


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    init()
    
    # funcName = "getLJQ"
    # funcName = "getMorseQH"
    # params = [ 3.0, 1.0, 0.0 ]
    # params = [ 3.0, 1.0, 0.0, 0.5, 1.0, 1.0 ]

    tests=[
        ("getLJQ",     [ 3.0, 1.0, 0.0, 1e-9 ] ),
        ("getMorseQH", [ 3.0, 1.0, 0.0, 0.5, 1.0, 1e-9 ] ),
        ("getMorseP4", [ 3.0, 1.0, 1.5 ] ),
    ]

    for funcName, params in tests:
        plotAll( funcName, params)

    plt.show()
    