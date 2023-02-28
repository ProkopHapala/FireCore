
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
"void getEs( double* Es, bool bRigid ){",
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
lib.init_types.argtypes  = [c_int, c_int_p, c_double_p, c_bool ] 
lib.init_types.restype   =  None
def init_types(typeMask, typREQs=None, bCopy=False ):
    ntyp = len(typeMask)
    return lib.init_types( ntyp, _np_as(typeMask,c_int_p), _np_as(typREQs,c_double_p), bCopy)

#  void setSystem( int isys, int na, int* types, double* ps, bool bCopy=false ){
lib.setSystem.argtypes  = [c_int, c_int, c_int_p, c_double_p, c_bool] 
lib.setSystem.restype   =  None
def setSystem(isys, types, ps, bCopy=False):
    na=len(types)
    return lib.setSystem(isys, na, _np_as(types,c_int_p), _np_as(ps,c_double_p), bCopy)

#  void setRigidSamples( int n, double* Es_, Mat3d* poses_, bool bCopy ){
lib.setRigidSamples.argtypes  = [c_int, c_double_p, c_double_p, c_bool, c_bool] 
lib.setRigidSamples.restype   =  None
def setRigidSamples(Es, poses, bCopy=False, bAlloc=False):
    if Es    is not None: n = len(Es)
    if poses is not None: n = len(poses)
    return lib.setRigidSamples(n, _np_as(Es,c_double_p), _np_as(poses,c_double_p), bCopy, bAlloc)

#  double run( int nstep, double ErrMax, double dt, bool bRigid ){
lib.run.argtypes  = [c_int, c_double, c_double, c_bool, c_int, c_bool ] 
lib.run.restype   =  c_double
def run(nstep, ErrMax, dt, bRigid, ialg=1, bRegularize=False, bClamp=False ):
    return lib.run(nstep, ErrMax, dt, bRigid, ialg, bRegularize )

#void getEs( double* Es, bool bRigid ){
lib.getEs.argtypes  = [c_double_p,  c_bool] 
lib.getEs.restype   =  None
def getEs( Es=None, bRigid=True):
    if Es is None: Es = np.zeros( nbatch )
    lib.getEs( _np_as(Es,c_double_p), bRigid)
    return Es

# =============== Buffers

#printBuffNames(){
lib.printBuffNames.argtypes = []
lib.printBuffNames.restype  = None
def printBuffNames():
    lib.printBuffNames()

#int* getIBuff(const char* name){ 
lib.getIBuff.argtypes = [c_char_p]
lib.getIBuff.restype  = c_int_p
def getIBuff(name,sh):
    if not isinstance(sh, tuple): sh=(sh,)
    name=name.encode('utf8')
    ptr = lib.getIBuff(name)
    return np.ctypeslib.as_array( ptr, shape=sh)

#double* getBuff(const char* name){ 
lib.getBuff.argtypes = [c_char_p]
lib.getBuff.restype  = c_double_p 
def getBuff(name,sh):
    if not isinstance(sh, tuple): sh=(sh,)
    name=name.encode('utf8')
    ptr = lib.getBuff(name)
    return np.ctypeslib.as_array( ptr, shape=sh)

def getBuffs():
    init_buffers()
    global ndims,nDOFs,ntype,nbatch,n0,n1
    ndims = getIBuff( "ndims", (6,) )  # [nDOFs,natoms,nnode,ncap,npi,nbonds]
    nDOFs=ndims[0]; ntype=ndims[1]; nbatch=ndims[2];n0=ndims[3];n1=ndims[4]; 
    print( "getBuffs(): nDOFs %i ntype %i nbatch %i n0 %i n1 %i" %(nDOFs,ntype,nbatch,n0,n1) )

    global DOFs,fDOFs,typeREQs,typeREQsMin,typeREQsMax,typeREQs0,typeKreg,typToREQ,weights,Es,poses,ps1,ps2,ps3, types1,types2,types3
    DOFs     = getBuff ( "DOFs",     nDOFs  )
    fDOFs    = getBuff ( "fDOFs",    nDOFs  )
    typToREQ = getIBuff( "typToREQ", (ntype,3)  )
    typeREQs0= getBuff ( "typeREQs0",(ntype,3)  )
    typeREQsMin= getBuff ( "typeREQsMin",(ntype,3)  )
    typeREQsMax= getBuff ( "typeREQsMax",(ntype,3)  )
    typeKreg = getBuff ( "typeKreg", (ntype,3)  )
    #weights = getBuff ( "weights",  nbatch )
    Es       = getBuff ( "Es",       nbatch ) 
    poses    = getBuff ( "poses",    (nbatch,3,3) )
    ps1      = getBuff ( "ps1",  (n0,3)    )
    ps2      = getBuff ( "ps2",  (n1,3)   )
    ps3      = getBuff ( "ps3",  (n1,3)   )
    types1   = getIBuff( "types1", n0 )
    types2   = getIBuff( "types2", n1 )
    types3   = getIBuff( "types3", n1 )

#  void init_buffers()
lib.init_buffers.argtypes  = []
lib.init_buffers.restype   =  None
def init_buffers():
    return lib.init_buffers()
