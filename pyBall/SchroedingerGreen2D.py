
import numpy as np
from   ctypes import c_int, c_double, c_bool, c_float, c_char_p, c_bool, c_void_p, c_char_p
import ctypes
import os
import sys

#sys.path.append('../')
#from pyMeta import cpp_utils 
from . import cpp_utils_ as cpp_utils
#import cpp_utils_ as cpp_utils

c_double_p = ctypes.POINTER(c_double)
c_int_p    = ctypes.POINTER(c_int)

def _np_as(arr,atype):
    if arr is None:
        return None
    else: 
        return arr.ctypes.data_as(atype)

cpp_utils.s_numpy_data_as_call = "_np_as(%s,%s)"

# ===== To generate Interfaces automatically from headers call:
header_strings = [
"void init(int nx, int ny){",
"double step( double E0, double dt ){",
]
#cpp_utils.writeFuncInterfaces( header_strings );        exit()     #   uncomment this to re-generate C-python interfaces

cpp_utils.BUILD_PATH = os.path.normpath( cpp_utils.PACKAGE_PATH + '../../cpp/Build/libs/Molecular' ) 
lib = cpp_utils.loadLib('SchroedingerGreen2D_lib', recompile=False)
array1ui = np.ctypeslib.ndpointer(dtype=np.uint32, ndim=1, flags='CONTIGUOUS')
array1i  = np.ctypeslib.ndpointer(dtype=np.int32,  ndim=1, flags='CONTIGUOUS')
array2i  = np.ctypeslib.ndpointer(dtype=np.int32,  ndim=2, flags='CONTIGUOUS')
array1d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS')
array2d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=2, flags='CONTIGUOUS')
array3d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=3, flags='CONTIGUOUS')

# ====================================
# ========= Globals
# ====================================

isInitialized = False

# ====================================
# ========= C functions
# ====================================

#double* getBuff(const char* name){ 
lib.getBuff.argtypes = [c_char_p]
lib.getBuff.restype  = c_double_p 
def getBuff(name,sh):
    if not isinstance(sh, tuple): sh=(sh,)
    name=name.encode('utf8')
    ptr = lib.getBuff(name)
    return np.ctypeslib.as_array( ptr, shape=sh)

def getBuffs():
    lib.init_buffers()
    global EQF
    EQF    = getBuff( "EQF",    (6,) ) # [ E,Q,F2sum; ]
    global V,source,psi,fpsi
    V      = getBuff( "V",      (ny,nx) )
    source = getBuff( "source", (ny,nx) )
    psi    = getBuff( "psi",    (ny,nx) )
    fpsi   = getBuff( "fpsi",   (ny,nx) )

#  void init(int nx_, int ny_){
lib.init.argtypes  = [c_int, c_int] 
lib.init.restype   =  None
def init(nx_, ny_):
    global nx,ny
    nx=nx_;ny=ny_
    return lib.init(nx, ny)

#  double step( double E0, double dt ){
lib.step.argtypes  = [c_double, c_double] 
lib.step.restype   =  c_double
def step(E0=0.0, dt=0.1):
    return lib.step(E0, dt)

#  double step( double E0, double dt ){
lib.step_Green.argtypes  = [] 
lib.step_Green.restype   =  c_double
def step_Green():
    return lib.step_Green()

# ====================================
# ========= Test Functions
# ====================================

# ====================================
# ========= Python Functions
# ====================================

# ====================================
# ========= MAIN
# ====================================

# if __name__ == "__main__":
#     import matplotlib.pyplot as plt
#     plt.show()