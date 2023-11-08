
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
"double loadXYZ( char* fname, int n0, int* i0s, int ntest, int* itests, int* types0=0, int testtypes=0 ){",
"void setType(int i, double* REQ )",
"void getType(int i, double* REQ )",
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

def cstr( s ):
    if s is None: return None
    return s.encode('utf8')

#  void setVerbosity( int verbosity_, int idebug_ ){
lib.setVerbosity.argtypes  = [c_int, c_int] 
lib.setVerbosity.restype   =  None
def setVerbosity( verbosity=1, idebug=0 ):
    return lib.setVerbosity( verbosity, idebug )

#  void init_types(int ntyp, int* typeMask, double* typREQs ){
lib.init_types.argtypes  = [c_int, c_int_p, c_double_p, c_bool ] 
lib.init_types.restype   =  None
def init_types(typeMask, typREQs=None, bCopy=True ):
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
lib.run.argtypes  = [c_int, c_int, c_double, c_double, c_bool, c_int, c_bool ] 
lib.run.restype   =  c_double
def run( nstep, ErrMax=1e-6, dt=0.1, bRigid=False, imodel=1, ialg=1, bRegularize=False, bClamp=False ):
    return lib.run(imodel, nstep, ErrMax, dt, bRigid, ialg, bRegularize )

#void getEs( double* Es, bool bRigid ){
lib.getEs.argtypes  = [ c_int, c_double_p,  c_bool] 
lib.getEs.restype   =  c_double
def getEs( imodel=1, Es=None, bRigid=True):
    if Es is None: Es = np.zeros( nbatch )
    Eerr = lib.getEs( imodel, _np_as(Es,c_double_p), bRigid)
    return Es

#  double loadXYZ( char* fname, int n0, int* i0s, int ntest, int* itests, int* types0, int testtypes ){
lib.loadXYZ.argtypes  = [c_char_p, c_int, c_int_p, c_int, c_int_p, c_int_p, c_int_p ] 
lib.loadXYZ.restype   =  c_int
def loadXYZ( fname,  i0s, itests, types0=None, testtypes=None, fname_AtomTypes="data/AtomTypes.dat" ):
    global nbatch
    n0     = len( i0s    )
    ntest  = len( itests )
    i0s    = np.array(i0s   ,np.int32)
    itests = np.array(itests,np.int32)
    if(types0    is not None): types0    = np.array(types0   ,np.int32)
    if(testtypes is not None): testtypes = np.array(testtypes,np.int32)
    nbatch = lib.loadXYZ( cstr(fname), n0, _np_as(i0s,c_int_p), ntest, _np_as(itests,c_int_p), _np_as(types0,c_int_p), _np_as(testtypes,c_int_p), cstr(fname_AtomTypes) )
    return nbatch


#  void setType(int i, double* REQ )
lib.setType.argtypes  = [c_int, c_double_p] 
lib.setType.restype   =  None
def setType(i, REQ):
    REQ=np.array(REQ)
    return lib.setType(i, _np_as(REQ,c_double_p))

#  void getType(int i, double* REQ )
lib.getType.argtypes  = [c_int, c_double_p] 
lib.getType.restype   =  None
def getType(i, REQ=None):
    if(REQ is None): REQ=np.zeros(4)
    lib.getType(i, _np_as(REQ,c_double_p))
    return REQ

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
    if(ptr):
        return np.ctypeslib.as_array( ptr, shape=sh)
    else:
        return None
    
#double* getBuff(const char* name){ 
lib.getBuff.argtypes = [c_char_p]
lib.getBuff.restype  = c_double_p 
def getBuff(name,sh):
    if not isinstance(sh, tuple): sh=(sh,)
    name=name.encode('utf8')
    ptr = lib.getBuff(name)
    if(ptr):
        return np.ctypeslib.as_array( ptr, shape=sh)
    else:
        return None

def getBuffs():
    init_buffers()
    global ndims,nDOFs,ntype,nbatch,n0,n1, params
    ndims = getIBuff( "ndims", (6,) )  # [nDOFs,natoms,nnode,ncap,npi,nbonds]
    nDOFs=ndims[0]; ntype=ndims[1]; nbatch=ndims[2];n0=ndims[3];n1=ndims[4]; 
    print( "getBuffs(): nDOFs %i ntype %i nbatch %i n0 %i n1 %i" %(nDOFs,ntype,nbatch,n0,n1) )
    params     = getBuff( "params", 4 )

    global DOFs,fDOFs,typeREQs,typeREQsMin,typeREQsMax,typeREQs0,typeKreg,typToREQ,weights,Es,poses,ps1,ps2,ps3, types1,types2,types3
    DOFs       = getBuff ( "DOFs",     nDOFs  )
    fDOFs      = getBuff ( "fDOFs",    nDOFs  )
    typToREQ   = getIBuff( "typToREQ",    (ntype,4)  )
    typeREQs   = getBuff ( "typeREQs",    (ntype,4)  )
    typeREQs0  = getBuff ( "typeREQs0",   (ntype,4)  )
    typeREQsMin= getBuff ( "typeREQsMin", (ntype,4)  )
    typeREQsMax= getBuff ( "typeREQsMax", (ntype,4)  )
    typeKreg   = getBuff ( "typeKreg",    (ntype,4)  )
    
    #weights = getBuff ( "weights",  nbatch )
    Es       = getBuff ( "Es",       nbatch ) 
    poses    = getBuff ( "poses",  (nbatch,3,3) )
    ps1      = getBuff ( "ps1",    (n0,3)    )
    #ps2      = getBuff ( "ps2",  (n1,3)   )
    #ps3      = getBuff ( "ps3",  (n1,3)   )
    types1   = getIBuff( "types1", n0 )
    #types2   = getIBuff( "types2", n1 )
    #types3   = getIBuff( "types3", n1 )

#  void init_buffers()
lib.init_buffers.argtypes  = []
lib.init_buffers.restype   =  None
def init_buffers():
    return lib.init_buffers()

################## Python ###############

def EnergyFromXYZ(fname):
    fin = open(fname,'r')
    il=0
    nl=0
    Es = []
    xs = []
    for line in fin:
        if(il==0):
            nl=2+int(line.split()[0])
        else:
            if( il%nl==1 ):
                ws = line.split()
                Es.append(float(ws[3]))
                xs.append(float(ws[5]))
        il+=1
    Es = np.array(Es)
    xs = np.array(xs)
    return Es,xs






