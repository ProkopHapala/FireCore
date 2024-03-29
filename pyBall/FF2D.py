
import numpy as np
from   ctypes import c_int, c_double, c_bool, c_float, c_char_p, c_bool, c_void_p, c_char_p
import ctypes
import os
import sys

sys.path.append('../')
from . import cpp_utils_ as cpp_utils

c_double_p = ctypes.POINTER(c_double)
c_int_p    = ctypes.POINTER(c_int)

def _np_as(arr,atype):
    if arr is None:
        return None
    else: 
        return arr.ctypes.data_as(atype)

cpp_utils.s_numpy_data_as_call = "_np_as(%s,%s)"

verbosity = 0
dt_glob   = 0.1
bVel      = False

# ===== To generate Interfaces automatically from headers call:
header_strings = [
"int init( char* str ){",
"void toArrays( int* types, double* apos, int* neighs ){",
"int run( int n, double dt, double damp, double F2conv, bool bCleanForce){",
]
#cpp_utils.writeFuncInterfaces( header_strings );        exit()     #   uncomment this to re-generate C-python interfaces

cpp_utils.BUILD_PATH = os.path.normpath( cpp_utils.PACKAGE_PATH + '../../cpp/Build/libs/Molecular' ) 
lib = cpp_utils.loadLib('FF2D_lib', recompile=False)
array1ui = np.ctypeslib.ndpointer(dtype=np.uint32, ndim=1, flags='CONTIGUOUS')
array1i  = np.ctypeslib.ndpointer(dtype=np.int32,  ndim=1, flags='CONTIGUOUS')
array2i  = np.ctypeslib.ndpointer(dtype=np.int32,  ndim=2, flags='CONTIGUOUS')
array1d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS')
array2d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=2, flags='CONTIGUOUS')
array3d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=3, flags='CONTIGUOUS')
# ========= C functions

def cstr( s ):
    if s is None: return None
    return s.encode('utf8')

#  int init( char* str ){
lib.init.argtypes  = [c_char_p, c_int] 
lib.init.restype   =  c_int
def init(s,  seed=454 ):
    return lib.init( cstr(s), seed )

#  void toArrays( int* types, double* apos, int* neighs ){
lib.toArrays.argtypes  = [c_int_p, c_double_p, c_int_p] 
lib.toArrays.restype   =  None
def toArrays(types=None, apos=None, neighs=None):
    return lib.toArrays(_np_as(types,c_int_p), _np_as(apos,c_double_p), _np_as(neighs,c_int_p))

#  double step( double dt, double damp ){
lib.step.argtypes  = [c_double, c_double] 
lib.step.restype   =  c_double
def step( dt=0.1, damp=0.1 ):
    return lib.step( dt, damp )

#  int run( int n, double dt, double damp, double F2conv, bool bCleanForce){
lib.run.argtypes  = [c_int, c_double, c_double, c_double, c_bool] 
lib.run.restype   =  c_int
def run(n=1000, dt=0.1, damp=0.1, Fconv=1e-2, bCleanForce=False):
    return lib.run(n, dt, damp, Fconv, bCleanForce)

# =========  Tests
