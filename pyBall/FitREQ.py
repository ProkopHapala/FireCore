
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
"void loadTypes( const char* fname_ElemTypes, const char* fname_AtomTypes ){",
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
def setVerbosity(verbosity=1, idebug=0):
    return lib.setVerbosity(verbosity, idebug)

# void setSwitches( int EvalJ, int WriteJ, int CheckRepulsion, int Regularize, int Epairs){
lib.setSwitches.argtypes  = [c_int, c_int, c_int, c_int, c_int]
lib.setSwitches.restype   =  None    
def setSwitches(EvalJ=0, WriteJ=0, CheckRepulsion=0, Regularize=0, Epairs=0):
    return lib.setSwitches(EvalJ, WriteJ, CheckRepulsion, Regularize, Epairs)

# int export_Erefs( double* Erefs ){ 
lib.export_Erefs.argtypes  = [c_double_p]
lib.export_Erefs.restype   =  c_int
def export_Erefs(Erefs=None, n=None):
    if Erefs is None:
        if n is None: n = lib.export_Erefs(_np_as(None,c_double_p))
        Erefs = np.zeros(n)
    lib.export_Erefs(_np_as(Erefs,c_double_p))
    return Erefs

#void setWeights( int n, double* weights ){
lib.setWeights.argtypes  = [c_int, c_double_p]
lib.setWeights.restype   =  None
def setWeights(weights):
    n = len(weights)
    lib.setWeights(n, _np_as(weights,c_double_p))
    
#  double run( int nstep, double Fmax, double dt, int imodel_, int ialg, int iparallel, bool bClamp, double max_step ){
lib.run.argtypes  = [c_int, c_double, c_double, c_int, c_int, c_int, c_bool, c_double]
lib.run.restype   =  c_double
def run(nstep=1000, Fmax=1e-8, dt=0.01, imodel=0, ialg=2, iparallel=2, bClamp=False, max_step=0.05):
    return lib.run(nstep, Fmax, dt, imodel, ialg, iparallel, bClamp, max_step)

# double getEs( int imodel, double* Es, double* Fs, bool bOmp, bool bDOFtoTypes ){
lib.getEs.argtypes  = [c_int, c_double_p,  c_double_p, c_bool, c_bool]
lib.getEs.restype   =  c_double
def getEs(imodel=0, Es=None, Fs=None, bOmp=False, bDOFtoTypes=False, bEs=True, bFs=False ):
    if bEs and (Es is None): Es = np.zeros( nbatch )
    if bFs and (Fs is None): Fs = np.zeros( nDOFs  )
    Eerr = lib.getEs(imodel, _np_as(Es,c_double_p), _np_as(Es,c_double_p), bOmp, bDOFtoTypes)
    return Eerr, Es, Fs

# void scanParam( int iDOF, int imodel,  int n, double* xs,  double* Es, double* Fs, bool bRegularize ){
lib.scanParam.argtypes  = [c_int, c_int, c_int, c_double_p, c_double_p, c_double_p, c_bool]
lib.scanParam.restype   = None
def scanParam( iDOF, xs, Es=None, Fs=None, imodel=2, bRegularize=False ):
    n = len(xs)
    if Es is None: Es = np.zeros( n )
    if Fs is None: Fs = np.zeros( n )
    lib.scanParam(iDOF, imodel, n, _np_as(xs,c_double_p), _np_as(Es,c_double_p), _np_as(Fs,c_double_p), bRegularize )
    return Es,Fs

#void scanParam2D( int iDOFx, int iDOFy, int imodel, int nx, int ny, double* xs, double* ys,  double* Es, double* Fx, double* Fy, bool bRegularize ){
lib.scanParam2D.argtypes  = [c_int, c_int, c_int, c_int, c_int, c_double_p, c_double_p, c_double_p, c_double_p, c_double_p, c_bool]
lib.scanParam2D.restype   = None
def scanParam2D(iDOFx, iDOFy, xs, ys, Es=None, Fx=None, Fy=None, imodel=2, bRegularize=False):
    nx, ny = len(xs), len(ys)
    if Es is None: Es = np.zeros((ny,nx))
    if Fx is None: Fx = np.zeros((ny,nx))
    if Fy is None: Fy = np.zeros((ny,nx))
    lib.scanParam2D(iDOFx, iDOFy, imodel, nx, ny, _np_as(xs,c_double_p), _np_as(ys,c_double_p),    _np_as(Es,c_double_p), _np_as(Fx,c_double_p), _np_as(Fy,c_double_p), bRegularize)
    return Es, Fx, Fy

# void loadTypes_new( const char* fname_ElemTypes, const char* fname_AtomTypes ){
lib.loadTypes.argtypes  = [c_char_p, c_char_p]
lib.loadTypes.restype   =  None
def loadTypes(fEtypes="data/ElementTypes.dat", fAtypes="data/AtomTypes.dat"):
    return lib.loadTypes(cstr(fEtypes), cstr(fAtypes))

# int loadTypeSelection_walls( const char* fname ){
lib.loadTypeSelection.argtypes  = [c_char_p]
lib.loadTypeSelection.restype   =  c_int
def loadTypeSelection(fname="typeSelection.dat"):
    return lib.loadTypeSelection(cstr(fname))

#void loadWeights( const char* fname ){
lib.loadWeights.argtypes  = [c_char_p]
lib.loadWeights.restype   =  c_int
def loadWeights(fname="weights.dat"):
    return lib.loadWeights(cstr(fname))

# int loadXYZ_new( const char* fname, const char* fname_AtomTypes  ){
lib.loadXYZ.argtypes  = [c_char_p, c_bool, c_bool]
lib.loadXYZ.restype   =  c_int
def loadXYZ(fname, bAddEpairs=False, bOutXYZ=False):
    global nbatch
    nbatch = lib.loadXYZ(cstr(fname), bAddEpairs, bOutXYZ)
    return nbatch

#  void setTypeToDOFs(int i, double* REQ )
lib.setTypeToDOFs.argtypes  = [c_int, c_double_p] 
lib.setTypeToDOFs.restype   =  None
def setTypeToDOFs(i, REQ):
    REQ=np.array(REQ)
    return lib.setTypeToDOFs(i, _np_as(REQ,c_double_p))

#  void getTypeFromDOFs(int i, double* REQ )
lib.getTypeFromDOFs.argtypes  = [c_int, c_double_p] 
lib.getTypeFromDOFs.restype   =  None
def getTypeFromDOFs(i, REQ=None):
    if(REQ is None): REQ=np.zeros(4)
    lib.getTypeFromDOFs(i, _np_as(REQ,c_double_p))
    return REQ

# =============== Buffers

#  void init_buffers()
lib.init_buffers.argtypes  = []
lib.init_buffers.restype   =  None
def init_buffers():
    #print( "init_buffers()" )
    return lib.init_buffers()

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
        return np.ctypeslib.as_array(ptr, shape=sh)
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
        return np.ctypeslib.as_array(ptr, shape=sh)
    else:
        return None

#def getBuffs( nnode, npi, ncap, nbond, NEIGH_MAX=4 ):
def getBuffs( ):
    #print( "getBuffs()" )
    init_buffers()
    #print( "getBuffs().1" )
    global ndims,typToREQ
    ndims = getIBuff( "ndims", (6,) )  # [nDOFs,ntype,nbatch,n0,n1,imodel]
    global nDOFs,ntype,nbatch,n0,n1,imodel
    nDOFs=ndims[0]; ntype=ndims[1]; nbatch=ndims[2]; n0=ndims[3]; n1=ndims[4]; imodel=ndims[5]
    typToREQ      = getIBuff( "typToREQ",    (ntype,) )
    #print( "getBuffs().2" )
    global DOFs,fDOFs,vDOFs
    DOFs          = getBuff ( "DOFs",   (nDOFs,)  )
    fDOFs         = getBuff ( "fDOFs",  (nDOFs,)  ) 
    vDOFs         = getBuff ( "vDOFs",  (nDOFs,)  ) 
    global typeREQs,  typeREQsMin,   typeREQsMax
    typeREQs      = getBuff ( "typeREQs",    (ntype,) )
    typeREQsMin   = getBuff ( "typeREQsMin", (ntype,) )
    typeREQsMax   = getBuff ( "typeREQsMax", (ntype,) )
    global typeREQs0, typeREQs0_low, typeREQs0_high
    typeREQs0     = getBuff ( "typeREQs0",   (ntype,) )
    typeREQs0_low = getBuff ( "typeREQs0_low",   (ntype,) )
    typeREQs0_high= getBuff ( "typeREQs0_high",  (ntype,) )
    global typeKreg,  typeKreg_low,  typeKreg_high
    typeKreg      = getBuff ( "typeKreg",    (ntype,) )
    typeKreg_low  = getBuff ( "typeKreg_low", (ntype,) )
    typeKreg_high = getBuff ( "typeKreg_high", (ntype,) )
    #print( "getBuffs().3" )

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






