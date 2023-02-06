
import numpy as np
from   ctypes import c_int, c_double, c_bool, c_float, c_char_p, c_bool, c_void_p, c_char_p, c_void_p
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
"void init( void* W, int nframes ){",
]
#cpp_utils.writeFuncInterfaces( header_strings );        exit()     #   uncomment this to re-generate C-python interfaces

cpp_utils.BUILD_PATH = os.path.normpath( cpp_utils.PACKAGE_PATH + '../../cpp/Build/libs_SDL' ) 
lib = cpp_utils.loadLib('MolGUIlib', recompile=False)
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
lib.init.argtypes  = [c_void_p] 
lib.init.restype   =  None
def init( W, n=1000000 ):
    return lib.init( W )

#  int init( int nframes ){
lib.run.argtypes  = [c_int] 
lib.run.restype   =  None
def run( n=1000000 ):
    return lib.run(  n )

# =========  Tests
