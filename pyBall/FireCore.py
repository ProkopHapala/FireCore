
import numpy as np
from   ctypes import c_int, c_double, c_bool, c_float, c_char_p, c_bool, c_void_p, c_char_p
import ctypes as ct
import os
import sys

#sys.path.append('../')
#from pyMeta import cpp_utils 
import cpp_utils

c_double_p = ct.POINTER(c_double)
c_int_p    = ct.POINTER(c_int)

def fortReal(a):
    ct.byref(ct.c_double(a))

def _np_as(arr,atype):
    if arr is None:
        return None
    else: 
        return arr.ct.data_as(atype)

cpp_utils.s_numpy_data_as_call = "_np_as(%s,%s)"

# ===== To generate Interfaces automatically from headers call:
header_strings = [
"void init_buffers(){",
]
#cpp_utils.writeFuncInterfaces( header_strings );        exit()     #   uncomment this to re-generate C-python interfaces

# LFLAGS   = " -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lfftw3xf_gnu -lm -Bdynamic -lmpi "

libs = {}
def loadLib(name, dic=libs):
    lib = ct.CDLL( name, ct.RTLD_GLOBAL )
    dic[name] = lib
    return lib

mkldir = "/home/prokop/SW/intel/compilers_and_libraries_2019.5.281/linux/mkl/lib/intel64_lin/"

cpp_utils.BUILD_PATH = os.path.normpath( cpp_utils.PACKAGE_PATH + '/../build/' ) 
lib = cpp_utils.loadLib('FireCore', recompile=False, mode=ct.RTLD_GLOBAL )

'''
loadLib( mkldir+"libmkl_def.so"        )
loadLib( mkldir+"libmkl_core.so"       )
loadLib( mkldir+"libmkl_avx2.so", )
'''

array1ui = np.ctypeslib.ndpointer(dtype=np.uint32, ndim=1, flags='CONTIGUOUS')
array1i  = np.ctypeslib.ndpointer(dtype=np.int32,  ndim=1, flags='CONTIGUOUS')
array2i  = np.ctypeslib.ndpointer(dtype=np.int32,  ndim=2, flags='CONTIGUOUS')
array1d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS')
array2d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=2, flags='CONTIGUOUS')
array3d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=3, flags='CONTIGUOUS')
# ========= C functions

#  subroutine hello( )
lib.firecore_hello.argtypes  = [ ] 
lib.firecore_hello.restype   =  None
def firecore_hello():
    return lib.firecore_hello()

#  subroutine firecore_init( natoms_, atomTypes, atomsPos )
lib.sum2.argtypes  = [c_double_p ] 
lib.sum2.restype   =  c_double
def sum2( a ):
    a = ct.c_double(a)
    return lib.sum2( ct.byref(a) )

#  subroutine firecore_init( natoms_, atomTypes, atomsPos )
lib.sum2val.argtypes  = [c_double,c_double ] 
lib.sum2val.restype   =  c_double
def sum2val( a, b ):
    return lib.sum2val( a, b )

'''
#  subroutine firecore_init( natoms_, atomTypes, atomsPos )
lib.firecore_init.argtypes  = [c_int, c_int_p, c_double_p ] 
lib.firecore_init.restype   =  None
def firecore_init(natoms, atomTypes, atomPos ):
    atomTypes_ = atomTypes.ctypes. data_as(c_int_p)
    atomPos_   = atomPos  .ctypes.data_as(c_double_p)
    return lib.firecore_init(natoms, atomTypes_, atomPos_ )
'''

#  subroutine firecore_init( natoms_, atomTypes, atomsPos )
lib.firecore_init.argtypes  = [c_int, array1i, array2d ] 
lib.firecore_init.restype   =  None
def firecore_init(natoms, atomTypes, atomPos ):
    return lib.firecore_init(natoms, atomTypes, atomPos )

#  subroutine firecore_evalForce( nmax_scf, forces_ )
lib.firecore_evalForce.argtypes  = [c_int, array2d ] 
lib.firecore_evalForce.restype   =  None
def firecore_evalForce( nmax_scf=100, forces=None, natom=5 ):
    if forces is None:
        forces = np.zeros( (3,natom) )
    return lib.firecore_evalForce( nmax_scf, forces )

# ========= Python Functions

if __name__ == "__main__":

    # b = sum2   ( 5.0 ); print "sum2   ", b
    # b = sum2val( 5.0,6.0 ); print "sum2val", b

    #a = ct.c_double(5)
    #b = lib.sum2(ct.byref(a))
    #print b

    
    firecore_hello()
    
    natoms = 5
    #atomType = np.random.randint(6, size=natoms).astype(np.int32)
    #atomPos  = np.random.random((3,natoms))
    atomType = np.array([6,1,1,1,1]).astype(np.int32)
    atomPos  = np.array([
        [ 0.0,0.0,0.1  ],
        [ -1.0,+1.0,+1.0  ],
        [ +1.0,-1.0,+1.0  ],
        [ +1.0,+1.0,-1.0  ],
        [ -1.0,-1.0,-1.0  ],
    ])
    print ("atomType ", atomType)
    print ("atomPos  ", atomPos)
    firecore_init( natoms, atomType, atomPos )
    forces = firecore_evalForce(nmax_scf=100)
    print( "Python: Forces", forces )
    

