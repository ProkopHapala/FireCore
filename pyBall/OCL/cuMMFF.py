import ctypes as ct
from ctypes import c_int, c_float, c_char, c_void_p
import ctypes as ct
import numpy as np
import os
from .. import cpp_utils_ as cpp_utils

# Define the ctypes for the arrays we'll be working with
c_float_p = ct.POINTER(c_float)
c_int_p   = ct.POINTER(c_int)
c_char_p  = ct.POINTER(c_char)


os.environ['LD_PRELOAD'] = '/usr/lib/x86_64-linux-gnu/libstdc++.so.6'

stdcpp = ct.CDLL('/usr/lib/x86_64-linux-gnu/libstdc++.so.6', ct.RTLD_GLOBAL)
cudart = ct.CDLL('/usr/local/cuda/lib64/libcudart.so',       ct.RTLD_GLOBAL)

CUDA_BUILD_PATH = os.path.normpath("../../cpp/Build/apps_CUDA/")  
lib = cpp_utils.loadLib('cuMMFF_lib', recompile=False, is_cuda=True, path=CUDA_BUILD_PATH, mode=ct.RTLD_GLOBAL)
array1ui = np.ctypeslib.ndpointer(dtype=np.uint32, ndim=1, flags='CONTIGUOUS')
array1i  = np.ctypeslib.ndpointer(dtype=np.int32,  ndim=1, flags='CONTIGUOUS')
array2i  = np.ctypeslib.ndpointer(dtype=np.int32,  ndim=2, flags='CONTIGUOUS')
array1d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS')
array2d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=2, flags='CONTIGUOUS')
array3d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=3, flags='CONTIGUOUS')


# upload function
lib.upload.argtypes = [c_char_p, c_void_p, c_int]
lib.upload.restype  = c_int
def upload(name, data, elem_size=-1):
    name_    = name.encode('utf-8')
    return lib.upload(name_, data.ctypes.data_as(c_void_p), elem_size)

# download function
lib.download.argtypes = [c_char_p, c_void_p, c_int]
lib.download.restype  = c_int
def download(name, data, elem_size=-1):
    name_    = name.encode('utf-8')
    return lib.download(name_, data.ctypes.data_as(c_void_p), elem_size)


#void init(int nSystems_, int nAtoms_, int nnode_, int npbc_, int nMaxSysNeighs_) {
lib.init.argtypes = [c_int, c_int, c_int, c_int, c_int]
lib.init.restype  = None
def init( nAtoms, nnode, npbc=1, nMaxSysNeighs=4, nSystems=1):   
    lib.init(nSystems, nAtoms, nnode, npbc, nMaxSysNeighs)

# ======= Kernels

# run_cleanForceMMFFf4 function
lib.run_cleanForceMMFFf4.argtypes = []
lib.run_cleanForceMMFFf4.restype  = c_int
def run_cleanForceMMFFf4():
    return lib.run_cleanForceMMFFf4()

# run_getNonBond function
lib.run_getNonBond.argtypes = []
lib.run_getNonBond.restype  = c_int
def run_getNonBond():
    return lib.run_getNonBond()

# run_getMMFFf4 function
lib.run_getMMFFf4.argtypes = []
lib.run_getMMFFf4.restype  = c_int
def run_getMMFFf4():
    return lib.run_getMMFFf4()

def run_updateAtomsMMFFf4():
    return lib.run_updateAtomsMMFFf4()
# run_updateAtomsMMFFf4 function
lib.run_updateAtomsMMFFf4.argtypes = []
lib.run_updateAtomsMMFFf4.restype  = c_int

# run_printOnGPU function
lib.run_printOnGPU.argtypes = [c_int, c_int_p]
lib.run_printOnGPU.restype  = c_int
def run_printOnGPU(sys_index=0, mask=(1,1,1,1)):
    mask_cuda = (c_int * 4)(*mask)
    return lib.run_printOnGPU(sys_index, mask_cuda)

# synchronize function
lib.synchronize.argtypes = []
lib.synchronize.restype  = c_int
def synchronize():
    return lib.synchronize()
