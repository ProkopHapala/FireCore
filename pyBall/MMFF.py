
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
"void init_buffers()",
"void init_params(char* fatomtypes, char* fbondtypes)",
"void init_nonbond()",
"void buildFF( bool bNonBonded_, bool bOptimizer_ )",
"int loadmol( char* fname_mol )",
"void initWithMolFile(char* fname_mol, bool bNonBonded_, bool bOptimizer_ )",
"void setSwitches( int doAngles, int doPiPiT, int  doPiSigma, int doPiPiI, int doBonded_, int PBC, int CheckInvariants )",
"bool checkInvariants( double maxVcog, double maxFcog, double maxTg )",
"double eval()",
"bool relax( int niter, double Ftol )",
]
#cpp_utils.writeFuncInterfaces( header_strings );        exit()     #   uncomment this to re-generate C-python interfaces

#libSDL = ctypes.CDLL( "/usr/lib/x86_64-linux-gnu/libSDL2.so", ctypes.RTLD_GLOBAL )
#libGL  = ctypes.CDLL( "/usr/lib/x86_64-linux-gnu/libGL.so",   ctypes.RTLD_GLOBAL )


#cpp_name='CombatModels'
#cpp_utils.make(cpp_name)
#LIB_PATH      = os.path.dirname( os.path.realpath(__file__) )
#LIB_PATH_CPP  = os.path.normpath(LIB_PATH+'../../../'+'/cpp/Build/libs/'+cpp_name )
#lib = ctypes.CDLL( LIB_PATH_CPP+("/lib%s.so" %cpp_name) )

cpp_utils.BUILD_PATH = os.path.normpath( cpp_utils.PACKAGE_PATH + '../../cpp/Build_OCL/libs/Molecular' ) 
lib = cpp_utils.loadLib('MMFF_lib', recompile=False)
array1ui = np.ctypeslib.ndpointer(dtype=np.uint32, ndim=1, flags='CONTIGUOUS')
array1i  = np.ctypeslib.ndpointer(dtype=np.int32,  ndim=1, flags='CONTIGUOUS')
array2i  = np.ctypeslib.ndpointer(dtype=np.int32,  ndim=2, flags='CONTIGUOUS')
array1d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS')
array2d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=2, flags='CONTIGUOUS')
array3d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=3, flags='CONTIGUOUS')

# ====================================
# ========= C functions
# ====================================

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

#def getBuffs( nnode, npi, ncap, nbond, NEIGH_MAX=4 ):
def getBuffs( NEIGH_MAX=4 ):
    #natom=nnode+ncap
    #nvecs=natom+npi
    #nDOFs=nvecs*3
    global ndims,Es
    ndims = getIBuff( "ndims", (6,) )  # [nDOFs,natoms,nnode,ncap,npi,nbonds]
    Es    = getBuff ( "Es",    (6,) )  # [ Etot,Eb,Ea, Eps,EppT,EppI; ]
    global nDOFs,natoms,nnode,ncap,nbonds,npi,nvecs
    nDOFs=ndims[0]; natoms=ndims[1]; nnode=ndims[2];ncap=ndims[3];npi=ndims[4];nbonds=ndims[5];
    nvecs=natoms+npi
    print( "nbonds %i nvecs %i npi %i natoms %i nnode %i ncap %i" %(nbonds,nvecs,npi,natoms,nnode,ncap) )
    global DOFs,fDOFs,apos,fapos,pipos,fpipos,bond_l0,bond_k, Kneighs,bond2atom,aneighs
    #Ebuf     = getEnergyTerms( )
    DOFs      = getBuff ( "DOFs",     (nvecs,3)  )
    fDOFs     = getBuff ( "fDOFs",    (nvecs,3)  ) 
    apos      = getBuff ( "apos",     (natoms,3) )
    fapos     = getBuff ( "fapos",    (natoms,3) )
    pipos     = getBuff ( "pipos",    (npi,3)    )
    fpipos    = getBuff ( "fpipos",   (npi,3)    )
    bond_l0   = getBuff ( "bond_l0",  (nbonds)   )
    bond_k    = getBuff ( "bond_k",   (nbonds)   )
    Kneighs   = getBuff ( "Kneighs",  (nnode,NEIGH_MAX) )
    bond2atom = getIBuff( "bond2atom",(nbonds,2) )
    aneighs   = getIBuff( "aneighs",  (nnode,NEIGH_MAX) )

#  void init_buffers()
lib.init_buffers.argtypes  = []
lib.init_buffers.restype   =  None
def init_buffers():
    return lib.init_buffers()

#  void init_params(const char* fatomtypes, const char* fbondtypes)
lib.init_params.argtypes  = [c_char_p, c_char_p] 
lib.init_params.restype   =  None
def init_params(fatomtypes, fbondtypes):
    fatomtypes = fatomtypes.encode('utf8')
    fbondtypes = fbondtypes.encode('utf8')
    return lib.init_params(fatomtypes,fbondtypes)

#  void init_nonbond()
lib.init_nonbond.argtypes  = [] 
lib.init_nonbond.restype   =  None
def init_nonbond():
    return lib.init_nonbond()

#  void buildFF( bool bNonBonded_, bool bOptimizer_ )
lib.buildFF.argtypes  = [c_bool, c_bool] 
lib.buildFF.restype   =  None
def buildFF(bNonBonded=True, bOptimizer=True):
    return lib.buildFF(bNonBonded, bOptimizer)

#  int loadmol( const char* fname_mol )
lib.loadmol.argtypes  = [c_char_p] 
lib.loadmol.restype   =  c_int
def loadmol(fname_mol):
    fname_mol=fname_mol.encode('utf8')
    return lib.loadmol(fname_mol)

#  void initWithMolFile(char* fname_mol, bool bNonBonded_, bool bOptimizer_ )
lib.initWithMolFile.argtypes  = [c_char_p, c_bool, c_bool] 
lib.initWithMolFile.restype   =  None
def initWithMolFile(fname_mol, bNonBonded=True, bOptimizer=True):
    fname_mol = fname_mol.encode('utf8')
    return lib.initWithMolFile(fname_mol, bNonBonded, bOptimizer)

#  void setSwitches( int doAngles, int doPiPiT, int  doPiSigma, int doPiPiI, int doBonded_, int PBC, int CheckInvariants )
lib.setSwitches.argtypes  = [c_int, c_int, c_int , c_int, c_int, c_int, c_int] 
lib.setSwitches.restype   =  None
def setSwitches(doAngles=0, doPiPiT=0, doPiSigma=0, doPiPiI=0, doBonded=0, PBC=0, CheckInvariants=0):
    return lib.setSwitches(doAngles, doPiPiT, doPiSigma, doPiPiI, doBonded, PBC, CheckInvariants)

#  bool checkInvariants( double maxVcog, double maxFcog, double maxTg )
lib.checkInvariants.argtypes  = [c_double, c_double, c_double] 
lib.checkInvariants.restype   =  c_bool
def checkInvariants(maxVcog=1e-9, maxFcog=1e-9, maxTg=1e-1):
    return lib.checkInvariants(maxVcog, maxFcog, maxTg)

#  double eval()
lib.eval.argtypes  = [] 
lib.eval.restype   =  c_double
def eval():
    return lib.eval()

#  bool relax( int niter, double Ftol )
lib.relax.argtypes  = [c_int, c_double] 
lib.relax.restype   =  c_bool
def relax(niter, Ftol):
    return lib.relax(niter, Ftol)

# ====================================
# ========= Python Functions
# ====================================


# ====================================
# ========= MAIN
# ====================================

# if __name__ == "__main__":
#     import matplotlib.pyplot as plt
#     plt.show()