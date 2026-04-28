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


#int getBufferIndex(const char* name) {
lib.getBufferIndex.argtypes = [c_char_p]
lib.getBufferIndex.restype  = c_int
def getBufferIndex(name):
    name_    = name.encode('utf-8')
    return lib.getBufferIndex(name_)

#int uploadId  (int id, const void* h_data, int nbyte, int offset=0) { 
lib.uploadId.argtypes = [c_int, c_void_p, c_int, c_int]
lib.uploadId.restype  = c_int
def uploadId(id, data, nbyte, offset=0):
    return lib.uploadId(id, data.ctypes.data_as(c_void_p), nbyte, offset)

#int downloadId(int id, void*       h_data, int nbyte, int offset=0) {
lib.downloadId.argtypes = [c_int, c_void_p, c_int, c_int]
lib.downloadId.restype  = c_int
def downloadId(id, data, nbyte, offset=0):
    return lib.downloadId(id, data.ctypes.data_as(c_void_p), nbyte, offset)

# upload function
lib.upload.argtypes = [c_char_p, c_void_p, c_int, c_int]
lib.upload.restype  = c_int
def upload(name, data, nbyte=-1, offset=0):
    name_    = name.encode('utf-8')
    return lib.upload(name_, data.ctypes.data_as(c_void_p), nbyte, offset)

# download function
lib.download.argtypes = [c_char_p, c_void_p, c_int, c_int]
lib.download.restype  = c_int
def download(name, data, nbyte=-1, offset=0):
    name_    = name.encode('utf-8')
    return lib.download(name_, data.ctypes.data_as(c_void_p), nbyte, offset)


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

# int run_MD(int nstep){
lib.run_MD.argtypes = [c_int, c_int]
lib.run_MD.restype  = c_int
def run_MD(nstep, mask=0b111):
    return lib.run_MD(nstep, mask)

# ====== Python

def pack_system(mmff, iSys=0):
    #print("pack_system() iSys=%d" % iSys)
    nvecs   = mmff.nvecs
    natoms  = mmff.natoms
    nnode   = mmff.nnode
    float4_size = 4 * np.float32().itemsize
    int4_size   = 4 * np.int32().itemsize
    
    # Calculate offsets for different buffer types
    offset_atoms = iSys * nvecs * float4_size
    upload("apos",   mmff.apos,  nbyte=mmff.apos.nbytes,  offset=offset_atoms)
    upload("aforce", mmff.fapos, nbyte=mmff.fapos.nbytes, offset=offset_atoms)
    
    # REQs offset calculation
    offset_REQs = iSys * natoms * float4_size
    upload("REQs",   mmff.REQs, nbyte=mmff.REQs.nbytes, offset=offset_REQs)
    
    # neighs offset calculation
    offset_neighs = iSys * natoms * int4_size
    upload("neighs",    mmff.neighs,    nbyte=mmff.neighs.nbytes, offset=offset_neighs)
    upload("neighCell", mmff.neighCell, nbyte=mmff.neighCell.nbytes, offset=offset_neighs)
    
    # Parameters offset calculation
    offset_apars = iSys * nnode * float4_size
    upload("MMpars", mmff.apars, nbyte=mmff.apars.nbytes, offset=offset_apars)
    upload("BLs",   mmff.bLs,    nbyte=mmff.bLs.nbytes, offset=offset_apars)
    upload("BKs",   mmff.bKs,    nbyte=mmff.bKs.nbytes, offset=offset_apars)
    upload("Ksp",   mmff.Ksp,    nbyte=mmff.Ksp.nbytes, offset=offset_apars)
    upload("Kpp",   mmff.Kpp,    nbyte=mmff.Kpp.nbytes, offset=offset_apars)
    
    # Upload MD parameters
    upload("MDpars", np.array([mmff.dt, mmff.damp, mmff.Flimit], dtype=np.float32), nbyte=float4_size, offset=iSys*float4_size)