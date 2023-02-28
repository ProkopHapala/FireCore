
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
"void init_types(int nbatch, int ntyp, int* typeMask, double* typREQs ){",
"void setSystem( int isys, int na, int* types, double* ps, bool bCopy=false ){",
"void setRigidSamples( int n, double* Es_, Mat3d* poses_, bool bCopy ){",
"double run( int nstep, double ErrMax, double dt, bool bRigid ){",
]
#cpp_utils.writeFuncInterfaces( header_strings );        exit()     #   uncomment this to re-generate C-python interfaces

cpp_utils.BUILD_PATH = os.path.normpath( cpp_utils.PACKAGE_PATH + '../../cpp/Build/libs/Molecular' ) 
lib = cpp_utils.loadLib('FitREQ_lib', recompile=False)
array1ui = np.ctypeslib.ndpointer(dtype=np.uint32, ndim=1, flags='CONTIGUOUS')
array1i  = np.ctypeslib.ndpointer(dtype=np.int32,  ndim=1, flags='CONTIGUOUS')
array2i  = np.ctypeslib.ndpointer(dtype=np.int32,  ndim=2, flags='CONTIGUOUS')
array1d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS')
array2d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=2, flags='CONTIGUOUS')
array3d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=3, flags='CONTIGUOUS')

# ====================================
# ========= Globals
# ====================================

#isInitialized = False
#nfound = -1

# ====================================
# ========= C functions
# ====================================


#  void init_types(int ntyp, int* typeMask, double* typREQs ){
lib.init_types.argtypes  = [c_int, c_int_p, c_double_p] 
lib.init_types.restype   =  None
def init_types(typeMask, typREQs=None):
    ntyp = len(typeMask)
    return lib.init_types( ntyp, _np_as(typeMask,c_int_p), _np_as(typREQs,c_double_p))

#  void setSystem( int isys, int na, int* types, double* ps, bool bCopy=false ){
lib.setSystem.argtypes  = [c_int, c_int, c_int_p, c_double_p, c_bool] 
lib.setSystem.restype   =  None
def setSystem(isys, na, types, ps, bCopy=False):
    return lib.setSystem(isys, na, _np_as(types,c_int_p), _np_as(ps,c_double_p), bCopy)

#  void setRigidSamples( int n, double* Es_, Mat3d* poses_, bool bCopy ){
lib.setRigidSamples.argtypes  = [c_int, c_double_p, c_double_p, c_bool, c_bool] 
lib.setRigidSamples.restype   =  None
def setRigidSamples(Es, poses, bCopy=False, bAlloc=False):
    if Es    is not None: n = len(Es)
    if poses is not None: n = len(poses)
    return lib.setRigidSamples(n, _np_as(Es,c_double_p), _np_as(poses,c_double_p), bCopy, bAlloc)

#  double run( int nstep, double ErrMax, double dt, bool bRigid ){
lib.run.argtypes  = [c_int, c_double, c_double, c_bool] 
lib.run.restype   =  c_double
def run(nstep, ErrMax, dt, bRigid):
    return lib.run(nstep, ErrMax, dt, bRigid)
