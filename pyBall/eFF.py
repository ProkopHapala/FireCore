import numpy as np
from   ctypes import c_int, c_double, c_bool, c_float, c_char_p, c_bool, c_void_p, c_char_p
import ctypes
import os
import sys

sys.path.append('../')
from . import cpp_utils_ as cpp_utils

c_double_p = ctypes.POINTER(c_double)
c_int_p    = ctypes.POINTER(c_int)

def _np_as(arr,atype):
    if arr is None:
        return None
    else: 
        return arr.ctypes.data_as(atype)

cpp_utils.s_numpy_data_as_call = "_np_as(%s,%s)"

verbosity = 0
dt_glob   = 0.1
bVel      = False

default_path = "data/xyz/"
default_path = "data/"

# ===== To generate Interfaces automatically from headers call:
header_strings = [
#"void init_buffers(){",
#"bool load_xyz( const char* fname ){",
#"void init( int na, int ne ){",
#"void eval(){"
#"void info(){",
#"double* getEnergyPointer(){",
#"int*    getDimPointer   (){",
#"double* getBuff(const char* name){",
#"void setBuff(const char* name, double* buff){",
#"int* getIBuff(const char* name){",
#"void setIBuff(const char* name, int* buff){", 
#"void setPauliModel(int i){",
#"void setKPauli( double KPauli ){",
#"void initOpt( double dt, double damping, double f_limit ){",
#"int  run( int nstepMax, double dt, double Fconv=1e-6, int ialg=0 ){",
#"void writeTo_fgo( char const* filename, bool bVel, bool bAppend ){",
]
#cpp_utils.writeFuncInterfaces( header_strings );        exit()     #   uncomment this to re-generate C-python interfaces

#libSDL = ctypes.CDLL( "/usr/lib/x86_64-linux-gnu/libSDL2.so", ctypes.RTLD_GLOBAL )
#libGL  = ctypes.CDLL( "/usr/lib/x86_64-linux-gnu/libGL.so",   ctypes.RTLD_GLOBAL )


#cpp_name='CombatModels'
#cpp_utils.make(cpp_name)
#LIB_PATH      = os.path.dirname( os.path.realpath(__file__) )
#LIB_PATH_CPP  = os.path.normpath(LIB_PATH+'../../../'+'/cpp/Build/libs/'+cpp_name )
#lib = ctypes.CDLL( LIB_PATH_CPP+("/lib%s.so" %cpp_name) )

cpp_utils.BUILD_PATH = os.path.normpath( cpp_utils.PACKAGE_PATH + '../../cpp/Build/libs/Molecular' ) 
lib = cpp_utils.loadLib('eFF_lib', recompile=False)
array1ui = np.ctypeslib.ndpointer(dtype=np.uint32, ndim=1, flags='CONTIGUOUS')
array1i  = np.ctypeslib.ndpointer(dtype=np.int32,  ndim=1, flags='CONTIGUOUS')
array2i  = np.ctypeslib.ndpointer(dtype=np.int32,  ndim=2, flags='CONTIGUOUS')
array1d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS')
array2d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=2, flags='CONTIGUOUS')
array3d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=3, flags='CONTIGUOUS')
# ========= C functions

def cstr( s ):
    if s is None: return None
    return s.encode('utf8')

#  void setVerbosity( int verbosity_, int idebug_ ){
lib.setVerbosity.argtypes  = [c_int, c_int] 
lib.setVerbosity.restype   =  None
def setVerbosity(verbosity_=0, idebug=0):
    global verbosity
    verbosity = verbosity_
    return lib.setVerbosity(verbosity, idebug)

#  void init_buffers(){
lib.init_buffers.argtypes  = [] 
lib.init_buffers.restype   =  None
def init_buffers():
    return lib.init_buffers() 

#  void load_xyz( const char* fname ){
lib.load_xyz.argtypes  = [c_char_p] 
lib.load_xyz.restype   =  c_bool
def load_xyz(fname):
    return lib.load_xyz(cstr(fname))
    #return lib.load_xyz(_np_as(fname,c_char_p)) 

def getBuffs( ):
    init_buffers()
    global ne,na,nDOFs, ndims, Es
    ndims = getIBuff( "ndims", (3,) ); 
    Es    = getBuff ( "Es",    (8,) ) # [0Etot,1Ek,2Eee,3EeePaul,4EeeExch,5Eae,6EaePaul,7Eaa]
    ne=ndims[0]; na=ndims[1]; nDOFs=ndims[2]     #;print("ne,na,nDOFs ", ne,na,nDOFs)
    global pDOFs, fDOFs, apos, aforce, epos, eforce, esize, fsize, aPars, espin
    pDOFs  = getBuff ( "pDOFs",  nDOFs )
    fDOFs  = getBuff ( "fDOFs",  nDOFs )
    apos   = getBuff ( "apos",   (na,3) )
    aforce = getBuff ( "aforce", (na,3) )
    epos   = getBuff ( "epos",   (ne,3) )
    eforce = getBuff ( "eforce", (ne,3) )
    esize  = getBuff ( "esize",   ne )
    fsize  = getBuff ( "fsize",   ne )
    aPars  = getBuff ( "aPars",   (na,4) )
    espin  = getIBuff( "espin",   ne )
    if(bVel):
        global vDOFs, avel, evel, vsize, invMasses, invAmass,invEmass,invSmass
        vDOFs     = getBuff ( "vDOFs",    nDOFs  )
        avel      = getBuff ( "avel",     (na,3) )
        evel      = getBuff ( "evel",     (ne,3) )
        vsize     = getBuff ( "vsize",     ne    )
        invMasses = getBuff ( "invMasses",nDOFs  )
        invAmass = getBuff ( "invAmass",  (na,3) )
        invEmass = getBuff ( "invEmass",  (ne,3) )
        invSmass = getBuff ( "invSmass",   ne    )

#  void load_xyz( const char* fname ){
lib.load_fgo.argtypes  = [c_char_p, c_bool, c_double] 
lib.load_fgo.restype   =  c_bool
def load_fgo(fname, bVel_=True, fUnits=1.):
    global bVel
    bVel=bVel_
    return lib.load_fgo( cstr(fname), bVel, fUnits)

#  void save_fgo( char const* filename, bool bVel, bool bAppend ){
lib.save_fgo.argtypes  = [c_char_p, c_bool, c_bool] 
lib.save_fgo.restype   =  None
def save_fgo(filename, bVel=False, bAppend=False):
    return lib.save_fgo( cstr(filename), bVel, bAppend)

#  void save_xyz( char const* filename, bool bVel, bool bAppend ){
lib.save_xyz.argtypes  = [c_char_p, c_bool] 
lib.save_xyz.restype   =  None
def save_xyz(filename, bAppend=False):
    return lib.save_xyz( cstr(filename), bAppend)

#  void setTrjName( char* trj_fname_ ){ 
lib.setTrjName.argtypes  = [c_char_p, c_int ] 
lib.setTrjName.restype   =  c_bool
def setTrjName(trj_fname_="trj.xyz", savePerNsteps=1, bDel=True ):
    if bDel: open(trj_fname_,"w").close()
    global trj_fname
    trj_fname=cstr(trj_fname_)
    return lib.setTrjName( trj_fname, savePerNsteps )

#  void init( int na, int ne ){
lib.init.argtypes  = [c_int, c_int] 
lib.init.restype   =  None
def init(na, ne, bVel_=False):
    global bVel
    bVel=bVel
    return lib.init(na, ne) 

#  void eval(){
lib.eval.argtypes  = [] 
lib.eval.restype   =  c_double
def eval():
    return lib.eval() 

#void evalFuncDerivs( int ie, int n, double* r, double* s, double* Es, double* Fr, double* Fs ){
lib.evalFuncDerivs.argtypes = [ c_int, c_int, array1d, array1d, array1d, array1d, array1d ]
lib.evalFuncDerivs.restype  = None
def evalFuncDerivs( r, s, Es=None, Fs=None, Fr=None, ie=0 ):
    r = r + s*0
    s = s + r*0
    n = len(r)
    if Es is None: Es=np.zeros(n)
    if Fr is None: Fr=np.zeros(n)
    if Fs is None: Fs=np.zeros(n) 
    lib.evalFuncDerivs( ie, n, r, s, Es, Fr, Fs )
    return Es,Fr,Fs

#  void info(){

lib.printSwitches.argtypes  = [] 
lib.printSwitches.restype   =  None
def printSwitches():
    return lib.printSwitches() 

lib.info.argtypes  = [] 
lib.info.restype   =  None
def info():
    return lib.info() 

#  double* getEnergyPointer(){
lib.getEnergyPointer.argtypes  = [] 
lib.getEnergyPointer.restype   = c_double_p
def getEnergyTerms( sh=(7,) ):
    # Ek=0, Eee EeePaul EeeExch Eae EaePaul Eaa
    ptr = lib.getEnergyPointer()
    return  np.ctypeslib.as_array( ptr, shape=sh )

#int*    getDimPointer   (){
lib.getDimPointer.argtypes  = [] 
lib.getDimPointer.restype   = c_int_p
def getDimPointer( sh=(3,) ):
    # ne=0 na=0 nDOFs=0
    ptr = lib.getDimPointer()
    return  np.ctypeslib.as_array( ptr, shape=sh )

#  double* getBuff(const char* name){
lib.getBuff.argtypes  = [c_char_p] 
lib.getBuff.restype   =  c_double_p
def getBuff( name, sh ):
    ptr = lib.getBuff(cstr(name))
    if not isinstance(sh, tuple): sh=(sh,)
    #sh_ = (natom,)
    #if sh is not None:
    #    sh_ = sh_ + sh
    #print "DEBUG type( ptr ) ", type( ptr ), sh
    return np.ctypeslib.as_array( ptr, shape=sh)

#  void setBuff(const char* name, double* buff){
lib.setBuff.argtypes  = [c_char_p, c_double_p] 
lib.setBuff.restype   =  None
def setBuff(name, buff):
    return lib.setBuff( cstr(name), _np_as(buff,c_double_p)) 
    #return lib.setBuff(_np_as(name,c_char_p), _np_as(buff,c_double_p)) 

#  int* getIBuff(const char* name){
lib.getIBuff.argtypes  = [c_char_p] 
lib.getIBuff.restype   =  c_int_p
def getIBuff(name,sh):
    ptr = lib.getIBuff(cstr(name))
    if not isinstance(sh, tuple): sh=(sh,)
    return np.ctypeslib.as_array( ptr, shape=sh)
    #return lib.getIBuff(_np_as(name,c_char_p)) 

#  void setIBuff(const char* name, int* buff){
lib.setIBuff.argtypes  = [c_char_p, c_int_p] 
lib.setIBuff.restype   =  None
def setIBuff(name, buff):
    return lib.setIBuff(name, _np_as(buff,c_int_p)) 
    #return lib.setIBuff(_np_as(name,c_char_p), _np_as(buff,c_int_p)) 

#  void setPauliModel(int i){
lib.setPauliModel.argtypes  = [c_int] 
lib.setPauliModel.restype   =  None
def setPauliModel(i):
    return lib.setPauliModel(i) 

#  void setKPauli( double KPauli ){
lib.setKPauli.argtypes  = [c_double] 
lib.setKPauli.restype   =  None
def setKPauli(KPauli):
    return lib.setKPauli(KPauli) 

# void setAtomParams( int n, const double* params_, bool bCopy=true, int mode=1 ){ 
lib.setAtomParams.argtypes = [ c_int, c_double_p, c_bool, c_int ]
lib.setAtomParams.restype  = None
def setAtomParams( params, bCopy=True, mode=1 ):
    n = len(params)
    if bCopy: params = np.array(params, dtype=np.double)
    return lib.setAtomParams( n, _np_as(params,c_double_p), bCopy, mode )

#void setSwitches( int bEvalKinetic, int bEvalCoulomb, int  bEvalPauli, int bEvalAA, int bEvalAE, int bEvalAECoulomb, int bEvalAEPauli, int bCoreCoul, int bEvalCoreCorect ){
lib.setSwitches.argtypes = [ c_int, c_int, c_int, c_int, c_int, c_int, c_int, c_int, c_int ]
lib.setSwitches.restype  = None
def setSwitches( kinetic=0, coulomb=0, pauli=0, AA=0, AE=0, AECoulomb=0, AEPauli=0, coreCoul=0, coreCorect=0 ):
    lib.setSwitches( kinetic, coulomb, pauli, AA, AE, AECoulomb, AEPauli, coreCoul, coreCorect )

#void setup( int isetup ){
lib.setup.argtypes = [ c_int ]
lib.setup.restype  = None
def setup( isetup ):
    lib.setup( isetup )

#  void initOpt( double dt, double damping, double f_limit ){
lib.initOpt.argtypes  = [c_double, c_double, c_double, c_bool ] 
lib.initOpt.restype   =  None
def initOpt(dt=0.1, damping=0.1, f_limit=1000.0, bMass=False ):
    global dt_glob
    dt_glob = dt
    return lib.initOpt(dt, damping, f_limit, bMass)

#  int  run( int nstepMax, double dt, double Fconv=1e-6, int ialg=0 ){
lib. run.argtypes  = [c_int, c_double, c_double, c_int, c_double_p, c_double_p] 
lib. run.restype   =  c_int
def  run(nstepMax=1000, dt=None, Fconv=1e-6, ialg=0, outE=None, outF=None, bOutE=True, bOutF=True):
    if dt is None: dt=dt_glob
    if (outE is None) and bOutE: outE=np.full(nstepMax, np.nan)
    if (outF is None) and bOutF: outF=np.full(nstepMax, np.nan)
    lib.run(nstepMax, dt, Fconv, ialg, _np_as(outE,c_double_p), _np_as(outF,c_double_p) )
    return  outE, outF


# void set_constrains( int nfix, Quat4d* fixed_poss, Vec2i* fixed_inds, bool bRealloc=true  ){
lib.set_constrains.argtypes  = [c_int, c_double_p, c_int_p, c_bool ]
lib.set_constrains.restype   =  None
def set_constrains( nfix, fixed_poss, fixed_inds, bRealloc=True ):
    return lib.set_constrains( nfix, _np_as(fixed_poss, c_double_p), _np_as(fixed_inds, c_int_p), bRealloc )

#void relaxed_scan( int nconf, int nfix, double* fixed_poss, int* fixed_inds_, double* outEs, double* apos_, double* epos_, int nstepMax, double dt, double Fconv, int ialg, char* scan_trj_name ){
lib.relaxed_scan.argtypes  = [c_int, c_int, c_double_p, c_int_p, c_double_p, c_double_p, c_double_p, c_int, c_double, c_double, c_int, c_char_p ]
lib.relaxed_scan.restype   =  None
def relaxed_scan( fixed_poss, fixed_inds, outEs=None, apos=None, epos=None, nstepMax=1000, dt=1e-2, Fconv=1e-6, ialg=0, scan_trj_name="scan.xyz" ):
    nconf, nfix, _ = fixed_poss.shape
    if apos is None: apos = np.zeros( (nconf, na, 3) )
    if epos is None: epos = np.zeros( (nconf, ne, 4) )
    if outEs is None: outEs = np.zeros( (nconf,8) )
    lib.relaxed_scan( nconf, nfix, _np_as(fixed_poss, c_double_p), _np_as(fixed_inds, c_int_p), _np_as(outEs, c_double_p), _np_as(apos, c_double_p), _np_as(epos, c_double_p), nstepMax, dt, Fconv, ialg, cstr(scan_trj_name) )
    return apos, epos, outEs

# void evalNumDerivs( double* Fnum, double d ){
lib. evalNumDerivs.argtypes  = [array1d, c_double] 
lib. evalNumDerivs.restype   =  None
def evalNumDerivs( Fnum=None, d=0.01):
    if Fnum is None: Fnum=np.zeros(nDOFs)
    lib.evalNumDerivs( Fnum, d )
    return Fnum


#void sample_ee( int n, double* RSs_, double* FEout_, int spin, double* KRSrho_, bool bEvalCoulomb, bool bEvalPauli, int iPauliModel ){
lib.sample_ee.argtypes  = [c_int, array2d, array2d, c_int, array1d, c_bool, c_bool, c_int ]
lib.sample_ee.restype   =  None
def sample_ee( RSs, spin, FEout=None, KRSrho=[1.125,0.9,-0.2], bEvalCoulomb=True, bEvalPauli=True, iPauliModel=1 ):
    n = len(RSs)

    FEout  = np.zeros((n,4)) #TOHLE VYKRESLIT
    KRSrho = np.array(KRSrho)
    lib.sample_ee(n, RSs, FEout, spin, KRSrho, bEvalCoulomb, bEvalPauli, iPauliModel )
    print(n)
    print("sample_ee eFF.py")
#     if FEout is None: FEout = np.zeros((n, 4), dtype=np.float64, order='C')  # Quat4d: [fx, fy, fz, E]
#     #print("RSs.shape", RSs.shape)
#     #print("FEout.shape", FEout.shape)
#     KRSrho = np.array(KRSrho, dtype=np.float64)
    return FEout

#void sample_EA( int n, double* RSs_, double* FEout_, double* KRSrho_,  double* aPar_,  bool bEvalAECoulomb, bool bCoreCoul, bool bEvalAEPauli ){
lib.sample_EA.argtypes  = [c_int, array2d, array2d, array1d, array1d, c_bool, c_bool, c_bool ]
lib.sample_EA.restype   =  None
def sample_EA( RSs, FEout=None, KRSrho=[1.125,0.9,-0.2], aPar=[4.,0.1,0.1,2.0], bEvalAECoulomb=True, bCoreCoul=True, bEvalAEPauli=True ):
    n = len(RSs)
    if FEout is None: FEout = np.zeros((n, 3), dtype=np.float64, order='C')
    #print("RSs.shape", RSs.shape)
    #print("FEout.shape", FEout.shape)
    KRSrho = np.array(KRSrho, dtype=np.float64)
    aPar = np.array(aPar, dtype=np.float64)
    lib.sample_EA(n, RSs, FEout, KRSrho, aPar, bEvalAECoulomb, bCoreCoul, bEvalAEPauli)
    return FEout

#int processXYZ( const char* fname, double Rfac=-0.5, double* outEs=0, double* apos_=0, double* epos_=0, int nstepMax=1000, double dt=0.001, double Fconv=1e-3, int ialg=2, bool bAddEpairs=false, bool bCoreElectrons=true, bool bChangeCore=true, bool bChangeEsize=true, const char* xyz_out="processXYZ.xyz", const char* fgo_out="processXYZ.fgo" ){
lib.processXYZ.argtypes  = [c_char_p, c_double, c_double_p, c_double_p, c_double_p, c_int, c_double, c_double, c_int, c_bool, c_bool,c_bool, c_bool, c_char_p, c_char_p ]
lib.processXYZ.restype   =  c_int # The C++ function returns an int (number of configurations processed)
def processXYZ( fname, Rfac=-1.35, outEs=None, apos=None, epos=None, nstepMax=1000, dt=0.5e-2, Fconv=1e-3, ialg=2, bAddEpairs=False, bCoreElectrons=False, bChangeCore=True, bChangeEsize=True, xyz_out="processXYZ.xyz", fgo_out="processXYZ.fgo",  bOutputs=(0,0,0) ):
    if bOutputs[0] and outEs is None: outEs = np.zeros(8, dtype=np.float64)
    if bOutputs[1] and apos  is None: apos  = np.zeros( (na, 3) )
    if bOutputs[2] and epos  is None: epos  = np.zeros( (ne, 4) )    
    print("xyz_out ", xyz_out )
    print("fgo_out ", fgo_out )
    lib.processXYZ( cstr(fname), Rfac, _np_as(outEs, c_double_p), _np_as(apos, c_double_p), _np_as(epos, c_double_p), nstepMax, dt, Fconv, ialg, bAddEpairs, bCoreElectrons, bChangeCore, bChangeEsize, cstr(xyz_out), cstr(fgo_out))
    return outEs, apos, epos

#int preAllocateXYZ(const char* fname, double Rfac=-0.5, bool bCoreElectrons=true )
lib.preAllocateXYZ.argtypes = [c_char_p, c_double, c_bool]
lib.preAllocateXYZ.restype  = c_int
def preAllocateXYZ(fname, Rfac=-0.5, bCoreElectrons=True):
    """Pre-initialize eFF from a single-config XYZ without dynamics"""
    return lib.preAllocateXYZ(cstr(fname), Rfac, bCoreElectrons)

lib.processXYZ_e.argtypes  = [c_char_p, c_double_p, c_double_p, c_double_p, c_int, c_double, c_double, c_int, c_char_p, c_char_p ]
lib.processXYZ_e.restype   =  c_int
def processXYZ_e( fname, outEs=None, apos=None, epos=None, nstepMax=0, dt=0.001, Fconv=1e-3, optAlg=2, xyz_out="processXYZ.xyz", fgo_out="processXYZ.fgo", bOutputs=(0,0,0) ):
    """
    Process XYZ file with electrons
    Returns: outEs, apos, epos
    """    
    # Get number of atoms and electrons from file (first and second line)
    if (bOutputs[1] and apos  is None) or (bOutputs[2] and epos  is None):
        with open(fname) as f:
            nae = int(f.readline().strip().split()[0])
            ws  = f.readline().strip().split()
            na  = int(ws[1])
            ne  = int(ws[2])
            if nae != na + ne: raise Exception(f"nae({nae}) != na({na}) + ne({ne}) while reading `{fname}`" )
            nconf = 1
            for line in f:
                ws = line.strip().split()
                if len(ws) == 1: nconf += 1
    if bOutputs[0] and outEs is None: outEs = np.zeros( (nconf,5) )
    if bOutputs[1] and apos  is None: apos  = np.zeros( (nconf, na, 3) )
    if bOutputs[2] and epos  is None: epos  = np.zeros( (nconf, ne, 4) )
    lib.processXYZ_e( cstr(fname), _np_as(outEs, c_double_p), _np_as(apos, c_double_p), _np_as(epos, c_double_p), nstepMax, dt, Fconv, optAlg, cstr(xyz_out), cstr(fgo_out) )
    return outEs, apos, epos

# void setKRSrho(double* _KRSrho){
lib.setKRSrho.argtypes  = [c_double_p] 
lib.setKRSrho.restype   =  None
def setKRSrho( KRSrho ):
    _KRSrho = np.array( KRSrho, dtype=np.float64 )
    return lib.setKRSrho( _np_as(_KRSrho, c_double_p) )




#int preAllocateFGO(const char* fname, bool bVel, double fUnits)
lib.preAllocateFGO.argtypes = [c_char_p, c_bool, c_double]
lib.preAllocateFGO.restype  = c_int
def preAllocateFGO(fname, bVel_=True, fUnits=1.):
    global bVel
    bVel = bVel_
    return lib.preAllocateFGO(cstr(fname), bVel, fUnits)

#int processFGO(const char* fname, bool bVel, double fUnits, double* outEs, double* apos, Quat4d* epos, int nstepMax, double dt, double Fconv, int ialg, bool bOutXYZ, bool bOutFGO)
lib.processFGO.argtypes = [c_char_p, c_double, c_double_p, c_double_p, c_double_p, c_int, c_double, c_double, c_int, c_bool, c_bool]
lib.processFGO.restype  = c_int
def processFGO(fname, fUnits=1., outEs=None, apos=None, epos=None, nstepMax=1000, dt=0.001, Fconv=1e-3, ialg=2, bOutXYZ=False, bOutFGO=False):
    return lib.processFGO(cstr(fname), fUnits, _np_as(outEs, c_double_p), _np_as(apos, c_double_p), _np_as(epos, c_double_p), nstepMax, dt, Fconv, ialg, bOutXYZ, bOutFGO)

# =========  Tests

def printEs():
    #print( " Etot %g Ek %g Eee %g EeePaul %g Eae %g EaePaul %g Eaa %g [eV]" %(Es[0],Es[1],Es[2],Es[3],Es[5],Es[6],Es[7]) )  # [0Etot,1Ek,2Eee,3EeePaul,4EeeExch,5Eae,6EaePaul,7Eaa]
    print("eFF.py::printEs() Etot ",Es[0],"Ek",Es[1],"Eee",Es[2],"EeePaul",Es[3],"Eae",Es[5],"EaePaul",Es[6],"Eaa", Es[7], "[eV]" )  # [0Etot,1Ek,2Eee,3EeePaul,4EeeExch,5Eae,6EaePaul,7Eaa]

def printAtoms():
    for i in range(na):
        print( "atom[%i] Q %g apos(%g,%g,%g)" %(i, aPars[i,0], apos[i,0],apos[i,1],apos[i,2]) )

def printElectrons():
    for i in range(ne):
        print( "electron[%i] apos(%g,%g,%g) size %g spin %i" %(i, epos[i,0],epos[i,1],epos[i,2], esize[i], espin[i]) )

def getNearestAtom( p, apos, ino=-1 ):
    r2s = np.sum( (apos[:,:]-p[None,:])**2, axis=1)
    if ino >=0: r2s[ino]=1e+300
    imin = np.argmin(r2s)
    rmin = np.sqrt( r2s[imin] )
    return imin,rmin

def getNearestAtoms( apos, bPrint=False ):
    n=len(apos)
    imins=np.zeros(n,np.int32)
    rmins=np.zeros(n)
    for i in range(n):
        imins[i],rmins[i] = getNearestAtom( apos[i,:], apos, ino=i )
        if bPrint:
            print( "bond[%i,%i] L=%g" %(i,imins[i],rmins[i]) )
    return imins, rmins
    
def eval_mol(name, fUnits=1., bPrint=True ):
    load_fgo(default_path+name+".fgo" )                               # load molecule in  .fgo format (i.e. floating-gaussian-orbital)
    eval()
    if bPrint:
        getBuffs()
        print("\n BAF")
        printEs()

def eval_ee( r, si, sj, spin=-1 ):
    init( 0, 2 )
    getBuffs()
    epos[:,:]=0; epos[0,0]=r
    esize[0]=si; esize[1]=sj
    espin[0]=1;  espin[1]=spin
    #setSwitches( kinetic=0, coulomb=0, pauli=0 )
    eval()
    print( "Eee", Es[2], "Epaul", Es[3] )  # [0Etot,1Ek,2Eee,3EeePaul,4EeeExch,5Eae,6EaePaul,7Eaa]

def check_DerivsPauli( r0=0.5,r1=1.5, s0=0.5,s1=0.5,   sj=1.0, n=10, spin=1 ):
    init( 0, 2 )
    setPauliModel(1)
    setSwitches( kinetic=-1, coulomb=-1, pauli=1 )
    getBuffs()
    espin[0]=1;  espin[1]=spin
    epos[:,:]=0; esize[:]=sj;
    rs = np.linspace(r0,r1,n,endpoint=False); dr=rs[1]-rs[0]      #; print( "rs \n", rs, r0, r1 );
    ss = np.linspace(s0,s1,n,endpoint=False); ds=ss[1]-ss[0]      #; print( "ss \n", ss, s0, s1 );
    Es,Fr,Fs = evalFuncDerivs( rs, ss ); #print( Es, Fr, Fs )
    dE   = (Es[2:]-Es[:-2])*0.5
    #print( rs, ss, Es, Fr, Fs )
    import matplotlib.pyplot as plt 
    #F= Fr*dr + Fs*ds; plt.figure(); plt.plot( dE   ,':', label="dE" );   plt.plot( F,  label="F" );   plt.legend(); plt.grid()
    plt.figure(); plt.plot( rs[1:-1], dE/dr,':',lw=2, label="dEdr" ); plt.plot( rs, Fr*-1, label="Fr" ); plt.legend(); plt.grid()
    plt.figure(); plt.plot( ss[1:-1], dE/ds,':',lw=2, label="dEds" ); plt.plot( ss, Fs*-1, label="Fs" ); plt.legend(); plt.grid()
    plt.show()
    #return Es,Fr,Fs,


def checkNumDerivs(name, d=0.001, bDetail=False, bPrint=True):
    load_fgo(default_path+name+".fgo" )
    getBuffs()
    eval()
    Fana = fDOFs.copy()
    Fnum = evalNumDerivs( d=d)*-1
    Ferr = Fnum-Fana
    maxError    = np.max( np.abs(Ferr                     ) )
    maxErrorRel = np.max( np.abs(Ferr)/(np.abs(Fana)+1e-8 ) )
    if bDetail: print (Fnum); print( Fana); print( Ferr)
    if bPrint:  print ( "maxError [eV/A]", maxError, " [%]", maxErrorRel*100. )
    #for i in range(nDOFs):
    #    print( Fnum )
    return maxError, maxErrorRel

'''
def eval_Derivs( name, iderivs=None, d=0.01 ):
    load_fgo(default_path+name+".fgo" )
    if iderivs is None: iderivs = range(nDOFs)
    n=len(iderivs)
    Fana = np.zeros(n)
    Fnum = np.zeros(n)
    for i,id in enumerate(iderivs):
        o=pDOFs[id]
        pDOFs[id] = o-d; E1 = eval()
        pDOFs[id] = o+d; E2 = eval()
        Fnum[i]   = (E2-E1)/(2*d)
        pDOFs[id] = o
        eval()
        Fana[i]   = fDOFs[id] 
    return Fana,Fnum
'''

def check_Derivs_ie( name, ie=0, r0=0.5,r1=1.5, s0=0.5,s1=0.5, n=10 ):
    load_fgo(default_path+name+".fgo" ) 
    getBuffs()
    rs = np.linspace(r0,r1,n,endpoint=False); dr=rs[1]-rs[0]      ; print( "rs \n", rs, r0, r1 );
    ss = np.linspace(s0,s1,n,endpoint=False); ds=ss[1]-ss[0]      ; print( "ss \n", ss, s0, s1 );
    Es,Fr,Fs = evalFuncDerivs( rs, ss,  ie=ie ); #print( Es, Fr, Fs )
    dE   = (Es[2:]-Es[:-2])*0.5
    #print( rs, ss, Es, Fr, Fs )
    import matplotlib.pyplot as plt 
    #F= Fr*dr + Fs*ds; plt.figure(); plt.plot( dE   ,':', label="dE" );   plt.plot( F,  label="F" );   plt.legend(); plt.grid()
    plt.figure(); plt.plot( rs[1:-1], dE/dr,':',lw=2, label="dEdr" ); plt.plot( rs, Fr*-1, label="Fr" ); plt.legend(); plt.grid()
    plt.figure(); plt.plot( ss[1:-1], dE/ds,':',lw=2, label="dEds" ); plt.plot( ss, Fs*-1, label="Fs" ); plt.legend(); plt.grid()
    plt.show()


def relax_mol(name, dt=0.03,damping=0.1, bTrj=True, bResult=True, perN=1, bPrintLbonds=True, nMaxIter=10000, outE=None, outF=None, fUnits=1., bFixNuclei=False ):
    if outE==True: outE=np.zeros(nMaxIter)
    load_fgo(default_path+name+".fgo", fUnits=fUnits , bVel_=bFixNuclei)                               # load molecule in  .fgo format (i.e. floating-gaussian-orbital)
    initOpt(dt=dt,damping=damping )                              # initialize optimizer/propagator
    if(bTrj): setTrjName(name+"_relax.xyz", savePerNsteps=perN ) # setup output .xyz file to save trajectory of all atoms and electrons at each timespep (comment-out to ommit .xyz and improve performance ) 
    getBuffs()
    print("Got buffs")
    #print("invMasses", invMasses )
    if(bFixNuclei): invAmass[:]=0 
    #print("invMasses", invMasses )
    nstep=run( nMaxIter, Fconv=1e-3, ialg=2, outE=outE, outF=outF ) # run simuation for maximum 1000 time steps until it converge to |F|<1e-3, ialg=2 is FIRE http://users.jyu.fi/~pekkosk/resources/pdf/FIRE.pdf   https://www.sciencedirect.com/science/article/pii/S0927025620300756
    if(verbosity>0):printEs()
    if bPrintLbonds:
        getNearestAtoms( apos, bPrint=False ) # bPrint defines, whether to print bond lengths
    if(bResult): 
        result_name=name+"_relaxed.fgo"
        if(verbosity>0): print("Optimized molecule saved to ", result_name)
        save_fgo( result_name )                 # save final relaxed geometry to .fgo format (i.e. floating-gaussian-orbital).
    print(outE)
    print("nstep", nstep)
    if outE is not None: 
        return outE[:nstep] #HUH?
    else:
        return nstep

def scan_core_size( name, core_sizes, dt=0.03,damping=0.1, nMaxIter=10000, fUnits=1., ia=0 ):
    load_fgo(default_path+name+".fgo", fUnits=fUnits )                  # load molecule in  .fgo format (i.e. floating-gaussian-orbital)
    initOpt(dt=dt,damping=damping )                                # initialize optimizer/propagator
    getBuffs()
    bondLengths = np.zeros(len(core_sizes));  bondLengths[:] = np.NaN
    result_name= name+"_corescan.fgo" ; open(result_name,'w').close()
    xyz_name   = name+"_corescan.xyz" ; open(xyz_name   ,'w').close()
    setVerbosity(0)
    for i,csize in enumerate(core_sizes): 
        aPars[ia,2]=csize
        nstep=run( nMaxIter, Fconv=1e-3, ialg=2 )    # run simuation for maximum 1000 time steps intil it converge to |F|<1e-3, ialg=2 is FIRE http://users.jyu.fi/~pekkosk/resources/pdf/FIRE.pdf   https://www.sciencedirect.com/science/article/pii/S0927025620300756
        if(nstep>=nMaxIter): 
            print( "cannot coverge csize ", csize ); break;
        imin,bondLengths[i] = getNearestAtom( apos[0,:], apos, ino=0 )
        print(csize, bondLengths[i], Es[0], nstep, imin ); # eff.printEs()
        save_fgo( result_name, bAppend=True )                 # save final relaxed geometry to .fgo format (i.e. floating-gaussian-orbital).
        save_xyz( xyz_name, bAppend=True )
    print( "bondLengths", bondLengths )
    return bondLengths

def check_H2(bRelax=True, name="H2_eFF", bPyeff=True):
    '''
    # limit Euu(r=0,si=sj) = 1/s^2 = 1/1.7007535132532356^2 = 0.34571422244 [Hartree] = 9.407362451147032 [eV]
    #   !!! But this should be infinite ( two same-spin electrons cannot occupy same orbital !!!
    '''
    print( "#======= check_H2(%s)" %name )
    load_fgo(default_path+name+".fgo" )  
    getBuffs()
    if bRelax:
        #relax_mol(name)
        initOpt(dt=0.03,damping=0.01 ) 
        run( 10000, Fconv=1e-3, ialg=2 )  
        bond_length = np.sqrt( ( (apos[0,:] - apos[1,:])**2 ).sum() ) ; print("apos\n",apos)
        Etot        = Es[0]
        print( "check_H2 E %g [eV] lbond %g [A]" %(Etot, bond_length) )
    else:
        if bPyeff:
            from pyBall import eFF_terms as effpy
            print(  "effpy.run_H2_molecule.__doc__:\n", effpy.run_H2_molecule.__doc__ )
            r = np.sqrt(((epos[0]-epos[1])**2).sum())
            effpy.pyeff_E_up_up( 1.125*np.sqrt(r**2+1e-8)/0.5291772105638411, 0.9*esize[0]/0.5291772105638411, 0.9*esize[1]/0.5291772105638411, rho=-0.2)
            Euu,Eud,DT,S = effpy.pyeff_EPaul( r, esize[0], esize[1], rho=-0.2 );           # print( "!!!!! pyeff_EPaul():   Euu",Euu, "Eud",Eud, "DT",DT, "S",S )
        eval()
        #eval_mol(name)
        #eval_mol("H2_eFF_relaxed")
    printAtoms()
    printElectrons()
    printEs()
    
def init_eff( natom_=0, nelec_=1, s=0.5,  aQ=1.0,aQs=0.0,aP=0.0,aPs=0.1 ):
    global natom,nelec
    natom=natom_; nelec=nelec_; 
    init( natom, nelec )
    aPar  = getBuff( "aPars",(natom,4) )
    apos  = getBuff( "apos",(natom,3) )
    epos  = getBuff( "epos",(nelec,3) )
    esize = getBuff( "esize",(nelec)  )
    aPar [:,0]=aQ;aPar[:,1]=aQs;aPar[:,2]=aPs;aPar[:,3]=aP;
    apos [:,:] = 0
    epos [:,:] = 0
    esize[:]   = s
    '''
    epos [:,:,:]= 0              + (np.random.rand(norb,perOrb,3)-0.5)*rnd_pos
    esize[:,:  ]=sz              + (np.random.rand(norb,perOrb  )-0.5)*rnd_size
    ecoef[:,:  ]=1               + (np.random.rand(norb,perOrb  )-0.5)*rnd_coef
    rhoP [:,:,:]=0               + (np.random.rand(norb,nqOrb,3 )-0.5)*rnd_pos
    rhoS [:,:  ]=sz*np.sqrt(0.5) + (np.random.rand(norb,nqOrb   )-0.5)*rnd_size
    rhoQ [:,:  ]=1               + (np.random.rand(norb,nqOrb   )-0.5)*rnd_coef
    '''

def scan_size( ss, ie0 ):
    Escan = np.zeros(   (len(ss),len(Es)) )
    for i,s in enumerate(ss):
        esize[ie0] = s
        lib.eval()
        Escan[i,:] =  Es[:]
    return Escan

def test_Hatom( bDerivs=False ):
    from . import eFF_terms as effpy
    import matplotlib.pyplot as plt 
    init_eff( natom_=1, nelec_=1, s=0.5 )
    ss = np.arange( 0.3,1.0, 0.01 )

    print(  "effpy.run_Hatom.__doc__:\n", effpy.run_Hatom.__doc__ )
    Ek_ref,Eae_ref = effpy.Hatom( ss );  E_ref=Ek_ref+Eae_ref 
    if bDerivs:
        E,dE   = evalFuncDerivs(1,ss)
        xs=ss
        plt.plot(xs      ,dE                               ,'-',label="F_ana")
        plt.plot(xs[1:-1],(E[2:]-E[:-2])/(-2*(xs[1]-xs[0])),':',label="F_num")
    else:
        getBuffs()
        Es = scan_size( ss, 0 );   E=Es[:,0]; Ek=Es[:,1]; Eae=Es[:,5]     # [0Etot,1Ek,2Eee,3EeePaul,4EeeExch,5Eae,6EaePaul,7Eaa]
        plt.plot( ss, Eae,    '-r', label="Eae" )
        plt.plot( ss, Ek,     '-b', label="Ek"  )     
        plt.plot( ss, Eae_ref,':r', lw=2, label="Eae_ref" )
        plt.plot( ss, Ek_ref, ':b', lw=2, label="Ek_ref"  )   
    
    i0ref,E0ref,x0ref = effpy.getExmin1D(E_ref,ss) ; print( "Hatom opt Reference: E %g [eV] s %g [A]" %(E0ref,x0ref) )
    i0   ,E0   ,x0    = effpy.getExmin1D(E    ,ss) ; print( "Hatom opt Numerical: E %g [eV] s %g [A]" %(E0   ,x0   ) )

    plt.plot( ss, E_ref,':k',lw=3, label="E_ref" )
    plt.plot( ss, E,    'grey'    ,lw=3,label="E"     )

    plt.xlabel('size [A]')
    plt.legend()
    plt.grid()
    plt.show()

# ========= Python Functions

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    load_xyz("../../cpp/sketches_SDL/Molecular/data/e2_eFF.xyz")
    #load_xyz("../../cpp/sketches_SDL/Molecular/data/H2O_eFF.xyz")
    info()
    eval()

    plt.show()