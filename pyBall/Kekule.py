
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
'void setVerbosity( int verbosity_, int idebug_ ){',
'void init_buffers(){',
'void init( int natom, int nbond, int* bond2atom, double* atomValence0 ){',
'void setDefaultBondOrders( double min, double max ){',
'void pinBondOrders(int n, int* ibonds, int* target ){',
'double relax( double dt, double F2conv=1e-6, int maxIter=1000 ){',
]
#cpp_utils.writeFuncInterfaces( header_strings );        exit()     #   uncomment this to re-generate C-python interfaces

cpp_utils.BUILD_PATH = os.path.normpath( cpp_utils.PACKAGE_PATH + '../../cpp/Build/libs/Molecular' ) 
lib = cpp_utils.loadLib('Kekule_lib', recompile=False)
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

def getBuffs( ):
    global ndims,Es,Ks,natoms,nbonds
    ndims = getIBuff( "ndims", (2,) ) # [natom,nbond]
    Es    = getBuff ( "Es",    (3,) ) # [ Etot,Eb,Ea; ]
    Ks    = getBuff ( "Ks",    (4,) ) # [ Etot,Eb,Ea; ]
    natoms=ndims[0]; nbonds=ndims[1]  ; print(natoms,nbonds, ndims )
    global bond2atom,  bondOrder,bondOrderMin,bondOrderMax,   atomValence,atomValenceMin,atomValenceMax

    bondOrder      = getBuff( "bondOrder",    nbonds )
    bond2atom      = getIBuff( "bond2atom",(nbonds,2) )
    

    bondOrderMin   = getBuff( "bondOrderMin", nbonds )
    bondOrderMax   = getBuff( "bondOrderMax", nbonds )
    
    atomValence    = getBuff( "atomValence",    natoms )
    atomValenceMin = getBuff( "atomValenceMin", natoms )
    atomValenceMax = getBuff( "atomValenceMax", natoms )


#  void setVerbosity( int verbosity_, int idebug_ ){
lib.setVerbosity.argtypes  = [c_int, c_int] 
lib.setVerbosity.restype   =  None
def setVerbosity(verbosity=0, idebug=0):
    return lib.setVerbosity(verbosity, idebug)

#  void init_buffers(){
lib.init_buffers.argtypes  = [] 
lib.init_buffers.restype   =  None
def init_buffers():
    return lib.init_buffers()

#  void init( int natom, int nbond, int* bond2atom, double* atomValence0 ){
lib.init.argtypes  = [c_int, c_int, c_int_p, c_double_p, c_int] 
lib.init.restype   =  None
def init(natom, nbond, bond2atom, atomValence0, seed=4445454 ):
    return lib.init(natom, nbond, _np_as(bond2atom,c_int_p), _np_as(atomValence0,c_double_p), seed)

#  void setDefaultBondOrders( double min, double max ){
lib.setDefaultBondOrders.argtypes  = [c_double, c_double] 
lib.setDefaultBondOrders.restype   =  None
def setDefaultBondOrders(min=1, max=2):
    return lib.setDefaultBondOrders(min, max)

#  void pinBondOrders(int n, int* ibonds, int* target ){double relax( double dt, double F2conv=1e-6, int maxIter=1000 ){
lib.pinBondOrders.argtypes  = [c_int, c_int_p, c_int_p] 
lib.pinBondOrders.restype   =  None
def pinBondOrders(n, ibonds, target):
    return lib.pinBondOrders(n, _np_as(ibonds,c_int_p), _np_as(target,c_int_p))

#  double eval(){
lib.eval.argtypes  = [] 
lib.eval.restype   =  c_double
def eval():
    return lib.eval()

#  double relax( double dt, double F2conv=1e-6, int maxIter=1000 ){
lib.relax.argtypes  = [c_double, c_double, c_int, c_bool, c_int] 
lib.relax.restype   =  c_double
def relax(dt=0.01, F2conv=1e-4, maxIter=1000, bRandStart=True, ialg=1):
    return lib.relax(dt, F2conv, maxIter, bRandStart, ialg)

#  double relax( double dt, double F2conv=1e-6, int maxIter=1000 ){
lib.relaxComb.argtypes  = [ c_bool] 
lib.relaxComb.restype   =  c_double
def relaxComb(bRandStart=True):
    return lib.relaxComb( bRandStart)

# ====================================
# ========= Test Functions
# ====================================

# ====================================
# ========= Python Functions
# ====================================

def testBond( n=100, d=0.1, Eout=None, ib=0, val0=None ):
    if(Eout is None):
        Eout=np.zeros(n)
    else:
        n=len(Eout)
    xs=np.zeros(n)
    if val0 is not None: bondOrder[ib]=val0
    for i in range(n):
        Eout[i] = lib.eval()
        bondOrder[ib]+=d
        xs[i]=bondOrder[ib]
    return Eout, xs

def runSystem( AOs, b2a , bPrintResults=True):
    AOs = np.array(AOs)
    b2a = np.array(b2a,dtype=np.int32)
    init( len(AOs), len(b2a), b2a, AOs, seed=np.random.randint(15454) )
    init_buffers()
    getBuffs()
    setDefaultBondOrders()
    #relax(maxIter=100, dt=0.5,  ialg=1)
    #Ks[0]=1;   # Katom
    #Ks[1]=1;   # Kbond
    #Ks[2]=0.0;   # KatomInt
    #Ks[3]=0.1;   # KbondInt
    #kek.Ks[0]=0;   # Katom
    #kek.Ks[1]=0;   # Kbond
    #kek.Ks[2]=0;   # KatomInt
    #kek.Ks[3]=0;   # KbondInt
    relax(maxIter=200, dt=0.5,  ialg=1)
    if(bPrintResults):
        print("============= Relaxed: ")
        print( "atomValence    ", atomValence    )
        print( "atomValenceMin ", atomValenceMin )
        print( "atomValenceMax ", atomValenceMax )
        print( "bondOrder      ", bondOrder      )
        print( "bondOrderMin   ", bondOrderMin   )
        print( "bondOrderMax   ", bondOrderMax   )

# ====================================
# ========= MAIN
# ====================================

# if __name__ == "__main__":
#     import matplotlib.pyplot as plt
#     plt.show()