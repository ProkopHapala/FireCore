
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
c_float_p  = ctypes.POINTER(c_float)
c_int_p    = ctypes.POINTER(c_int)

def _np_as(arr,atype):
    if arr is None:
        return None
    else: 
        return arr.ctypes.data_as(atype)

cpp_utils.s_numpy_data_as_call = "_np_as(%s,%s)"

# ===== To generate Interfaces automatically from headers call:
header_strings = [
#"void init_buffers()",
#"void init_params(char* fatomtypes, char* fbondtypes)",
#"void init_nonbond()",
#"void buildFF( bool bNonBonded_, bool bOptimizer_ )",
#"int loadmol( char* fname_mol )",
#"void initWithMolFile(char* fname_mol, bool bNonBonded_, bool bOptimizer_ )",
#"void setSwitches( int doAngles, int doPiPiT, int  doPiSigma, int doPiPiI, int doBonded_, int PBC, int CheckInvariants )",
#"bool checkInvariants( double maxVcog, double maxFcog, double maxTg )",
#"double eval()",
#"bool relax( int niter, double Ftol )",
#"void shift_atoms_ax( int n, int* selection, double* d                  )",
#"void shift_atoms   ( int n, int* selection, int ia0, int ia1, double l )",
#"void rotate_atoms_ax( int n, int* selection, double* p0, double* ax, double phi      )",
#"void rotate_atoms   ( int n, int* selection, int ia0, int iax0, int iax1, double phi )",
#"int  splitAtBond( int ib, int* selection )",
#"void scanTranslation_ax( int n, int* selection, double* vec, int nstep, double* Es, bool bWriteTrj )",
#"void scanTranslation( int n, int* selection, int ia0, int ia1, double l, int nstep, double* Es, bool bWriteTrj )",
#"void sampleNonBond(int n, double* rs, double* Es, double* fs, int kind, double*REQi_,double*REQj_, double K ){",
#"#void sampleSurf(char* name, int n, double* rs, double* Es, double* fs, int kind, double*REQ_, double K, double Rdamp ){",
#"void init( char* xyz_name, char* surf_name, char* smile_name, bool bMMFF=false, int* nPBC, double gridStep, char* sAtomTypes, char* sBondTypes, char* sAngleTypes ){",
#"void scanRotation_ax( int n, int* selection, double* p0, double* ax, double phi, int nstep, double* Es, bool bWriteTrj )",
#"void scanRotation( int n, int* selection,int ia0, int iax0, int iax1, double phi, int nstep, double* Es, bool bWriteTrj )",
#"void set_opt( double dt_max,  double dt_min, double damp_max, double finc,    double fdec,   double falpha, int minLastNeg, double cvf_min, double cvf_max){",
#"void sample_evalAngleCos( double K, double c0, int n, double* angles, double* Es, double* Fs ){",
#"void sample_DistConstr( double lmin, double lmax, double kmin, double kmax, double flim , int n, double* xs, double* Es, double* Fs ){",
#"void addDistConstrain(  int i0,int i1, double lmin,double lmax,double kmin,double kmax,double flim, double k ){",
#"void addAngConstrain(  int i0,int i1,int i2, double ang0, double k ){",
#"void change_lvec( double* lvec, bool bAdd, bool  ){",
#"void optimizeLattice_1d( double* dlvec, int n1, int n2, int initMode, double tol ){",
#"void optimizeLattice_2d_multi( double* dlvec, int nstesp, int initMode, double tol ){",
#"virtual void upload_pop( const char* fname ){",
# "void pack_system  ( int isys, bool bParams, bool bForces, bool bVel, bool blvec ){",
# "void unpack_system( int isys, MMFFsp3_loc& ff, bool bForces=false, bool bVel=false ){",
# "void upload_sys   ( int isys, bool bParams, bool bForces, bool bVel ){",
# "void download_sys ( int isys, bool bForces, bool bVel ){",
# "void upload       ( bool bParams, bool bForces, bool bVel ){",
# "void download     ( bool bForces, bool bVel ){",
]
#cpp_utils.writeFuncInterfaces( header_strings );        exit()     #   uncomment this to re-generate C-python interfaces

#libSDL = ctypes.CDLL( "/usr/lib/x86_64-linux-gnu/libSDL2.so", ctypes.RTLD_GLOBAL )
#libGL  = ctypes.CDLL( "/usr/lib/x86_64-linux-gnu/libGL.so",   ctypes.RTLD_GLOBAL )


#cpp_name='CombatModels'
#cpp_utils.make(cpp_name)
#LIB_PATH      = os.path.dirname( os.path.realpath(__file__) )
#LIB_PATH_CPP  = os.path.normpath(LIB_PATH+'../../../'+'/cpp/Build/libs/'+cpp_name )
#lib = ctypes.CDLL( LIB_PATH_CPP+("/lib%s.so" %cpp_name) )

cpp_utils.BUILD_PATH = os.path.normpath( cpp_utils.PACKAGE_PATH + '../../cpp/Build/libs_OCL' ) 
lib = cpp_utils.loadLib('MMFFmulti_lib', recompile=False)
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
glob_bMMFF    = True

# ====================================
# ========= C functions
# ====================================

#  void sample_DistConstr( double lmin, double lmax, double kmin, double kmax, double flim , int n, double* xs, double* Es, double* Fs ){
lib.sample_DistConstr.argtypes  = [c_double, c_double, c_double, c_double, c_double, c_int, c_double_p, c_double_p, c_double_p] 
lib.sample_DistConstr.restype   =  None
def sample_DistConstr( xs, lmin=1, lmax=1, kmin=1, kmax=1, flim=1e+300, Es=None, Fs=None):
    n = len(xs)
    if Es is None: Es=np.zeros(n)
    if Fs is None: Fs=np.zeros(n)
    lib.sample_DistConstr(lmin, lmax, kmin, kmax, flim, n, _np_as(xs,c_double_p), _np_as(Es,c_double_p), _np_as(Fs,c_double_p))
    return Es,Fs

#  void sample_evalPiAling( double K, double c0, int n, double* angles, double* Es, double* Fs ){
lib.sample_evalPiAling.argtypes  = [c_double, c_double, c_double, c_double, c_int, c_double_p, c_double_p, c_double_p] 
lib.sample_evalPiAling.restype   =  None
def sample_evalPiAling( angles, K=1.0, ang0=0.0, r1=1.,r2=1., Es=None, Fs=None):
    n = len(angles)
    if Es is None: Es=np.zeros(n)
    if Fs is None: Fs=np.zeros(n)
    lib.sample_evalPiAling(K, ang0, r1, r2, n, _np_as(angles,c_double_p), _np_as(Es,c_double_p), _np_as(Fs,c_double_p))
    return Es,Fs

#  void sample_evalAngleCos( double K, double c0, int n, double* angles, double* Es, double* Fs ){
lib.sample_evalAngleCos.argtypes  = [c_double, c_double, c_double, c_double, c_int, c_double_p, c_double_p, c_double_p] 
lib.sample_evalAngleCos.restype   =  None
def sample_evalAngleCos( angles, K=1.0, ang0=0.0, r1=1.,r2=1., Es=None, Fs=None):
    n = len(angles)
    if Es is None: Es=np.zeros(n)
    if Fs is None: Fs=np.zeros(n)
    lib.sample_evalAngleCos(K, ang0, r1, r2, n, _np_as(angles,c_double_p), _np_as(Es,c_double_p), _np_as(Fs,c_double_p))
    return Es,Fs

#  void sample_evalAngleCosHalf( double K, double c0, int n, double* angles, double* Es, double* Fs ){
lib.sample_evalAngleCosHalf.argtypes  = [c_double, c_double, c_double, c_double, c_int, c_double_p, c_double_p, c_double_p] 
lib.sample_evalAngleCosHalf.restype   =  None
def sample_evalAngleCosHalf( angles, K=1.0, ang0=0.0, r1=1.,r2=1., Es=None, Fs=None):
    n = len(angles)
    if Es is None: Es=np.zeros(n)
    if Fs is None: Fs=np.zeros(n)
    lib.sample_evalAngleCosHalf(K, ang0, r1, r2, n, _np_as(angles,c_double_p), _np_as(Es,c_double_p), _np_as(Fs,c_double_p))
    return Es,Fs

#  void sampleNonBond(int n, double* rs, double* Es, double* fs, int kind, double*REQi_,double*REQj_, double K ){
lib.sampleNonBond.argtypes  = [c_int, array1d, array1d, array1d, c_int, array1d, array1d, c_double, c_double ] 
lib.sampleNonBond.restype   =  None
def sampleNonBond( rs, Es=None, fs=None, kind=1, REQi=(1.487,0.0006808,0.0), REQj=(1.487,0.0006808,0.0), K=-1.0, Rdamp=1.0 ):
    n =len(rs)
    if Es is None: Es=np.zeros(n)
    if Fs is None: Fs=np.zeros(n)
    rs=np.array(rs)
    REQi=np.array(REQi)
    REQj=np.array(REQj) 
    lib.sampleNonBond(n, rs, Es, Fs, kind, REQi, REQj, K, Rdamp)
    return Es,fs

# void sampleSurf(char* name, int n, double* rs, double* Es, double* fs, int kind, double*REQ_, double K, double Rdamp ){
lib.sampleSurf.argtypes  = [c_char_p, c_int, array1d, array1d, array1d, c_int, c_int, c_double, c_double, c_double, array1d, c_bool] 
lib.sampleSurf.restype   =  None
def sampleSurf( name, rs, Es=None, fs=None, kind=1, atyp=0, Q=0.0, K=-1.0, Rdamp=1.0, pos0=(0.,0.,0.), bSave=False ):
    if name is not None: name=name.encode('utf8')
    n =len(rs)
    if Es is None: Es=np.zeros(n)
    if fs is None: fs=np.zeros(n)
    rs=np.array(rs)
    pos0=np.array(pos0)
    lib.sampleSurf( name, n, rs, Es, fs, kind, atyp, Q, K, Rdamp, pos0, bSave )
    return Es,fs

# void sampleSurf_vecs(char* name, int n, double* rs, double* Es, double* fs, int kind, double*REQ_, double K, double Rdamp ){
lib.sampleSurf_vecs.argtypes  = [c_char_p, c_int, array2d, array1d, array2d, c_int, c_int, c_double, c_double, c_double, array1d, c_bool] 
lib.sampleSurf_vecs.restype   =  None
def sampleSurf_vecs( name, poss, Es=None, fs=None, kind=1, atyp=0, Q=0.0, K=-1.0, Rdamp=1.0, pos0=(0.,0.,0.), bSave=False ):
    if name is not None: name=name.encode('utf8')
    n =len(poss)
    if Es is None: Es=np.zeros(n)
    if fs is None: fs=np.zeros((n,3))
    pos0=np.array(pos0)
    lib.sampleSurf_vecs( name, n, poss, Es, fs, kind, atyp, Q, K, Rdamp, pos0, bSave )
    return Es,fs

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

#double* getBuff(const char* name){ 
lib.getfBuff.argtypes = [c_char_p]
lib.getfBuff.restype  = c_float_p 
def getfBuff(name,sh):
    if not isinstance(sh, tuple): sh=(sh,)
    name=name.encode('utf8')
    ptr = lib.getfBuff(name)
    return np.ctypeslib.as_array( ptr, shape=sh )

#def getBuffs( nnode, npi, ncap, nbond, NEIGH_MAX=4 ):
def getBuffs( NEIGH_MAX=4 ):
    # int  nDOFs=0,nnode=0,ncap=0,nvecs=0;
    # double Etot,Eb,Ea, Eps,EppT,EppI;
    global ndims,Es
    ndims = getIBuff( "ndims", (6,) )  # [nDOFs,natoms,nnode,ncap,npi,nbonds]
    Es    = getBuff ( "Es",    (6,) )  # [ Etot,Eb,Ea, Eps,EppT,EppI; ]
    global nDOFs,natoms,nnode,ncap,npi,nvecs
    nDOFs=ndims[0]; nnode=ndims[1]; ncap=ndims[2];nvecs=ndims[3]; natoms=nnode+ncap; npi=nnode
    print( "getBuffs(): nDOFs %i nvecs %i  natoms %i nnode %i ncap %i npi %i" %(nDOFs,nvecs,natoms,nnode,ncap,npi) )
    global DOFs,fDOFs,vDOFs,apos,fapos,pipos,fpipos,bond_l0,bond_k, bond2atom,neighs,selection
    #Ebuf     = getEnergyTerms( )
    apos      = getBuff ( "apos",     (natoms,3) )
    fapos     = getBuff ( "fapos",    (natoms,3) )
    if glob_bMMFF:
        DOFs      = getBuff ( "DOFs",     (nvecs,3)  )
        fDOFs     = getBuff ( "fDOFs",    (nvecs,3)  ) 
        vDOFs     = getBuff ( "vDOFs",    (nvecs,3)  ) 
        pipos     = getBuff ( "pipos",    (npi,3)    )
        fpipos    = getBuff ( "fpipos",   (npi,3)    )
        #bond_l0   = getBuff ( "bond_l0",  (nbonds)   )
        #bond_k    = getBuff ( "bond_k",   (nbonds)   )
        #Kneighs   = getBuff ( "Kneighs",  (nnode,NEIGH_MAX) )
        #bond2atom = getIBuff( "bond2atom",(nbonds,2) )
        neighs   = getIBuff( "neighs",  (nnode,NEIGH_MAX) )
        selection = getIBuff( "selection",  (natoms) )
    
    # --- GPU buffers (multi-replicas)
    global gpu_neighs,gpu_neighCell,gpu_bkNeighs,  gpu_atoms, gpu_aforces, gpu_avel, gpu_constr,  gpu_REQs, gpu_MMpars, gpu_BLs, gpu_BKs, gpu_Ksp, gpu_Kpp, gpu_lvecs, gpu_ilvecs
    gpu_neighs    = getIBuff ( "gpu_neighs",    (nSys,nvecs,4)  )
    gpu_neighCell = getIBuff ( "gpu_neighCell", (nSys,nvecs,4)  ) 
    gpu_bkNeighs  = getIBuff ( "gpu_bkNeighs",  (nSys,nvecs,4)  ) 

    gpu_atoms     = getfBuff ( "gpu_atoms",     (nSys,nvecs,4)  )
    gpu_aforces   = getfBuff ( "gpu_aforces",   (nSys,nvecs,4)  ) 
    gpu_avel      = getfBuff ( "gpu_avel",      (nSys,nvecs,4)  ) 
    gpu_constr    = getfBuff ( "gpu_constr",    (nSys,natoms,4) )

    gpu_REQs      = getfBuff ( "gpu_REQs",      (nSys,natoms,3) )
    gpu_MMpars    = getfBuff ( "gpu_MMpars",    (nSys,nnode,3)  ) 
    gpu_BLs       = getfBuff ( "gpu_BLs",       (nSys,nnode,3)  ) 
    gpu_BKs       = getfBuff ( "gpu_BKs",       (nSys,nnode,3)  )
    gpu_Ksp       = getfBuff ( "gpu_Ksp",       (nSys,nnode,3)  )
    gpu_Kpp       = getfBuff ( "gpu_Kpp",       (nSys,nnode,3)  )

    gpu_lvecs     = getfBuff ( "gpu_lvecs",     (nSys,3,4)    )
    gpu_ilvecs    = getfBuff ( "gpu_ilvecs",    (nSys,3,4)    )
    #gpu_pbcshifts = getfBuff ( "gpu_pbcshifts", (nSys,npi,4)    )



#  void init_buffers()
lib.init_buffers.argtypes  = []
lib.init_buffers.restype   =  None
def init_buffers():
    return lib.init_buffers()

def cstr( s ):
    if s is None: return None
    return s.encode('utf8')

#  void setVerbosity( int verbosity_, int idebug_ ){
lib.setVerbosity.argtypes  = [c_int, c_int] 
lib.setVerbosity.restype   =  None
def setVerbosity( verbosity=1, idebug=0 ):
    return lib.setVerbosity( verbosity, idebug )

#  void* init( int nSys, char* xyz_name, char* surf_name, char* smile_name, bool bMMFF, bool bEpairs, int* nPBC, double gridStep, char* sAtomTypes, char* sBondTypes, char* sAngleTypes, double T, double gamma, int nExplore, int nRelax, double pos_kick, double vel_kick ){
lib.init.argtypes = [
    c_int,        # nSys
    c_char_p,     # xyz_name
    c_char_p,     # surf_name
    c_char_p,     # smile_name
    c_bool,       # bMMFF
    c_bool,       # bEpairs
    array1i,      # nPBC
    array1i,      # grid_nPBC
    c_double,     # gridStep
    c_char_p,     # sElementTypes
    c_char_p,     # sAtomTypes
    c_char_p,     # sBondTypes
    c_char_p,     # sAngleTypes
    c_double,     # T
    c_double,     # gamma
    c_int,        # nExplore
    c_int,        # nRelax
    c_double,     # pos_kick
    c_double,      # vel_kick
    c_int         # GridFF
]
lib.init.restype   =  c_void_p
def init(
        nSys_=10,
        xyz_name  ="input.xyz", 
        surf_name =None, 
        smile_name=None, 
        sElementTypes = "data/ElementTypes.dat",
        sAtomTypes = "data/AtomTypes.dat",
        sBondTypes = "data/BondTypes.dat",
        sAngleTypes= "data/AngleTypes.dat",
        bMMFF=True, bEpairs=False, nPBC=(1,1,0), gridnPBC=(1,1,0), gridStep=0.1,
        T = -1, gamma = -1,
        nExplore=0, nRelax=0, pos_kick=0.0, vel_kick=0.0,
        GridFF=5
    ):
    global glob_bMMFF, nSys
    nSys=nSys_
    glob_bMMFF = bMMFF
    # Convert integer tuples to numpy arrays for C compatibility
    nPBC = np.array(nPBC, dtype=np.int32)
    gridnPBC = np.array(gridnPBC, dtype=np.int32)
    return lib.init( nSys, cstr(xyz_name), cstr(surf_name), cstr(smile_name), bMMFF, bEpairs, nPBC, gridnPBC, gridStep, cstr(sElementTypes), cstr(sAtomTypes), cstr(sBondTypes), cstr(sAngleTypes), T, gamma, nExplore, nRelax, pos_kick, vel_kick, GridFF )

def tryInit():
    if not isInitialized:
        init()
        #getBuffs( NEIGH_MAX=4 )

#  void insertSMILES(char* s)
lib.insertSMILES.argtypes  = [c_char_p,c_bool] 
lib.insertSMILES.restype   =  None
def insertSMILES(s ):
    s = s.encode('utf8')
    return lib.insertSMILES(s)

#  void setSwitches_multi( int doAngles, int doPiPiT, int  doPiSigma, int doPiPiI, int doBonded_, int PBC, int CheckInvariants )
lib.setSwitches_multi.argtypes  = [c_int, c_int, c_int , c_int, c_int, c_int, c_int, c_int, c_int] 
lib.setSwitches_multi.restype   =  None
def setSwitches(doAngles=0, doPiPiT=0, doPiSigma=0, doPiPiI=0, doBonded=0, PBC=0, CheckInvariants=0, dovdW=0, bSaveToDatabase=0):
    return lib.setSwitches_multi(doAngles, doPiPiT, doPiSigma, doPiPiI, doBonded, PBC, CheckInvariants, dovdW, bSaveToDatabase)

#  bool checkInvariants( double maxVcog, double maxFcog, double maxTg )
lib.checkInvariants.argtypes  = [c_double, c_double, c_double] 
lib.checkInvariants.restype   =  c_bool
def checkInvariants(maxVcog=1e-9, maxFcog=1e-9, maxTg=1e-1):
    return lib.checkInvariants(maxVcog, maxFcog, maxTg)

#  void set_opt( double dt_max,  double dt_min, double damp_max, double finc,    double fdec,   double falpha, int minLastNeg, double cvf_min, double cvf_max){
lib.set_opt.argtypes  = [c_double, c_double, c_double, c_double, c_double, c_double, c_int, c_double, c_double] 
lib.set_opt.restype   =  None
def set_opt( dt_max=0.1, dt_min=0.02, damp_max=0.2, finc=1.1, fdec=0.5, falpha=0.8, minLastNeg=5, cvf_min=-0.1, cvf_max=+0.1 ):
    return lib.set_opt(dt_max, dt_min, damp_max, finc, fdec, falpha, minLastNeg, cvf_min, cvf_max)

#  void setOptLog( int n, double* cos, double* f, double* v, double* dt, double* damp ){
lib.setOptLog.argtypes  = [c_int, c_double_p, c_double_p, c_double_p, c_double_p, c_double_p] 
lib.setOptLog.restype   =  None
def setOptLog( n ):
    cos  = np.zeros(n) 
    f    = np.zeros(n) 
    v    = np.zeros(n) 
    dt   = np.zeros(n) 
    damp = np.zeros(n)
    lib.setOptLog(n, _np_as(cos,c_double_p), _np_as(f,c_double_p), _np_as(v,c_double_p), _np_as(dt,c_double_p), _np_as(damp,c_double_p))
    return cos,f,v,dt,damp

#  void setTrjName( char* trj_fname_ ){ 
lib.setTrjName.argtypes  = [c_char_p, c_int ] 
lib.setTrjName.restype   =  c_bool
def setTrjName(trj_fname_="trj.xyz", savePerNsteps=1, bDel=True ):
    if bDel: open(trj_fname_,"w").close()
    global trj_fname
    trj_fname=cstr(trj_fname_)
    return lib.setTrjName( trj_fname, savePerNsteps )

#  int saveXYZ( const char* fname, const char* comment)
lib.saveXYZ.argtypes  = [c_char_p, c_char_p, c_int ] 
lib.saveXYZ.restype   =  c_int
def saveXYZ(fname, comment, imod=1 ):
    return lib.saveXYZ( cstr(fname), cstr(comment), imod )

# char* getTypeName( int ia, bool fromFF){
lib.getTypeName.argtypes  = [c_int, c_bool ] 
lib.getTypeName.restype   =  c_char_p
def getTypeName( ia, fromFF=True ):
    s = lib.getTypeName( ia, fromFF )
    ss = s.decode()
    return ss 

#double eval()
lib.eval.argtypes  = [] 
lib.eval.restype   =  c_double
def eval():
    return lib.eval()

#  int  run( int nstepMax, double dt, double Fconv=1e-6, int ialg=0 ){
lib. run.argtypes  = [c_int, c_double, c_double, c_int, c_double_p, c_double_p, c_int ] 
lib. run.restype   =  c_int
def  run(nstepMax=1000, dt=-1, Fconv=1e-3, ialg=2, outE=None, outF=None, iParalel=1 ):
    return lib.run(nstepMax, dt, Fconv, ialg, _np_as(outE,c_double_p), _np_as(outF,c_double_p), iParalel )

# void MDloop( int nIter, double Ftol = -1, int iParalel = 3 )
lib.MDloop.argtypes  = [c_int, c_double, c_int, c_int ] 
lib.MDloop.restype   =  None
def MDloop( perframe=100, Ftol=-1, iParalel=3,  perVF=100 ):
    return lib.MDloop( perframe, Ftol, iParalel, perVF )

# ========= GPU Replicas management

#  void pack_system  ( int isys, bool bParams, bool bForces, bool bVel, bool blvec ){
lib.pack_system.argtypes = [c_int, c_bool, c_bool, c_bool, c_bool] 
lib.pack_system.restype  =  None
def pack_system(isys, bParams=False, bForces=False, bVel=False, blvec=True ):
    return lib.pack_system(isys, bParams, bForces, bVel, blvec)

#  void unpack_system( int isys, MMFFsp3_loc& ff, bool bForces=false, bool bVel=false ){
lib.unpack_system.argtypes = [c_int, c_bool, c_bool] 
lib.unpack_system.restype  =  None
def unpack_system(isys, bForces=True, bVel=True ):
    return lib.unpack_system(isys, bForces, bVel)

#  void upload_sys( int isys, bool bParams, bool bForces, bool bVel ){
lib.upload_sys.argtypes  = [c_int, c_bool, c_bool, c_bool, c_bool ] 
lib.upload_sys.restype   =  None
def upload_sys(isys, bParams=True, bForces=False, bVel=True, blvec=True ):
    return lib.upload_sys(isys, bParams, bForces, bVel, blvec)

#  void download_sys ( int isys, bool bForces, bool bVel ){
lib.download_sys .argtypes  = [c_int, c_bool, c_bool] 
lib.download_sys .restype   =  None
def download_sys (isys, bForces=True, bVel=True):
    return lib.download_sys (isys, bForces, bVel)

#  void upload( bool bParams, bool bForces, bool bVel ){
lib.upload.argtypes  = [c_bool, c_bool, c_bool, c_bool] 
lib.upload.restype   =  None
def upload(bParams=True, bForces=False, bVel=True, blvec=True):
    return lib.upload(bParams, bForces, bVel, blvec)

#  void download( bool bForces, bool bVel ){
lib.download.argtypes = [c_bool, c_bool] 
lib.download.restype  =  None
def download(bForces=True, bVel=True):
    return lib.download(bForces, bVel)

# ========= Lattice Optimization

#  void change_lvec( double* lvec, bool bAdd, bool  ){
lib.change_lvec.argtypes  = [c_double_p, c_bool, c_bool] 
lib.change_lvec.restype   =  None
def change_lvec( lvec, bAdd=False, bUpdatePi=False ):
    lvec=np.array( lvec,dtype=np.float64)
    return lib.change_lvec(_np_as(lvec,c_double_p), bAdd, bool)

#  void optimizeLattice_1d( double* dlvec, int n1, int n2, int initMode, double tol ){
lib.optimizeLattice_1d.argtypes  = [c_double_p, c_int, c_int, c_int, c_double] 
lib.optimizeLattice_1d.restype   =  None
def optimizeLattice_1d(dlvec, n1=5, n2=5, initMode=0, tol=1e-6):
    dlvec=np.array( dlvec,dtype=np.float64)
    return lib.optimizeLattice_1d(_np_as(dlvec,c_double_p), n1, n2, initMode, tol)

#  void optimizeLattice_2d_multi( double* dlvec, int nstesp, int initMode, double tol ){
lib.optimizeLattice_2d_multi.argtypes  = [c_double_p, c_int, c_int, c_double] 
lib.optimizeLattice_2d_multi.restype   =  None
def optimizeLattice_2d_multi(dlvec, nstesp=10, initMode=0, tol=1e-6):
    dlvec=np.array( dlvec,dtype=np.float64)
    return lib.optimizeLattice_2d_multi(_np_as(dlvec,c_double_p), nstesp, initMode, tol)

#  virtual void upload_pop( const char* fname ){
lib.upload_pop.argtypes  = [c_char_p] 
lib.upload_pop.restype   =  None
def upload_pop(fname):
    return lib.upload_pop( cstr(fname) )

# ========= Constrains

#  void addDistConstrain(  int i0,int i1, double lmin,double lmax,double kmin,double kmax,double flim, double k, double* shift ){
lib.addDistConstrain.argtypes  = [c_int, c_int, c_double, c_double, c_double, c_double, c_double, c_double_p ] 
lib.addDistConstrain.restype   =  None
def addDistConstrain( i0, i1, lmin=1, lmax=1, kmin=1, kmax=1, flim=1e+300, l=None, k=None, shift=(0.,0.,0.) ):
    if l is not None: 
        lmin=l; lmax=l
    if k is not None: 
        kmin=k; kmax=k
    shift=np.array(shift)
    return lib.addDistConstrain(i0, i1, lmin, lmax, kmin, kmax, flim, shift )

#  void addAngConstrain(  int i0,int i1,int i2, double ang0, double k ){
lib.addAngConstrain.argtypes  = [c_int, c_int, c_int, c_double, c_double] 
lib.addAngConstrain.restype   =  None
def addAngConstrain(i0, i1, i2, ang0=0.0, k=1.0):
    return lib.addAngConstrain(i0, i1, i2, ang0, k)

# ============= Manipulation

#void shift_atoms_ax( int n, int* sel, double* d                  )
lib.shift_atoms_ax.argtypes  = [c_int, c_int_p, c_double_p] 
lib.shift_atoms_ax.restype   =  None
def shift_atoms_ax(d, sel=None):
    n=0
    if sel is not None:
        n=len(sel)
        sel=np.array(sel,np.int32)
    if d is not None: d=np.array(d)
    return lib.shift_atoms_ax(n, _np_as(sel,c_int_p), _np_as(d,c_double_p))

#void shift_atoms( int n, int* sel, int ia0, int ia1, double l )
lib.shift_atoms.argtypes  = [c_int, c_int_p, c_int, c_int, c_double] 
lib.shift_atoms.restype   =  None
def shift_atoms( ia0, ia1, l, sel=None):
    n=0
    if sel is not None:
        n=len(sel)
        sel=np.array(sel,np.int32)
    return lib.shift_atoms(n, _np_as(sel,c_int_p), ia0, ia1, l)

#  rotate_atoms_ax( int n, int* sel, double* p0, double* ax, double phi      )
lib.rotate_atoms_ax.argtypes  = [c_int, c_int_p, c_double_p, c_double_p, c_double] 
lib.rotate_atoms_ax.restype   =  None
def rotate_atoms_ax( phi, sel=None, p0=None, ax=None ):
    n=0
    if sel is not None:
        n=len(sel)
        sel=np.array(sel,np.int32)
    if p0 is not None: p0=np.array(p0)
    if ax is not None: ax=np.array(ax)
    return lib.rotate_atoms_ax(n, _np_as(sel,c_int_p), _np_as(p0,c_double_p), _np_as(ax,c_double_p), phi)

#  rotate_atoms( int n, int* sel, int ia0, int iax0, int iax1, double phi )
lib.rotate_atoms_ax.argtypes  = [c_int, c_int_p, c_int, c_int, c_int, c_double] 
lib.rotate_atoms_ax.restype   =  None
def rotate_atoms_ax(ia0, iax0, iax1, phi, sel=None):
    n=0
    if sel is not None:
        n=len(sel)
        sel=np.array(sel,np.int32)
    return lib.rotate_atoms_ax(n, _np_as(sel,c_int_p), ia0, iax0, iax1, phi)

#  void scanTranslation_ax( int n, int* sel, double* vec, int nstep, double* Es, bool bWriteTrj )
lib.scanTranslation_ax.argtypes  = [c_int, c_int_p, c_double_p, c_int, c_double_p, c_bool] 
lib.scanTranslation_ax.restype   =  None
def scanTranslation_ax(nstep, Es, bWriteTrj, sel=None, vec=None, ):
    n=0
    if sel is not None:
        n=len(sel)
        sel=np.array(sel,np.int32)
    if vec is not None: vec=np.array(vec)
    return lib.scanTranslation_ax(n, _np_as(sel,c_int_p), _np_as(vec,c_double_p), nstep, _np_as(Es,c_double_p), bWriteTrj)

#  void scanTranslation( int n, int* sel, int ia0, int ia1, double l, int nstep, double* Es, bool bWriteTrj )
lib.scanTranslation.argtypes  = [c_int, c_int_p, c_int, c_int, c_double, c_int, c_double_p, c_bool] 
lib.scanTranslation.restype   =  None
def scanTranslation( ia0, ia1, l, nstep, Es=None, sel=None, bWriteTrj=False):
    if Es is None: Es=np.zeros(nstep)
    return lib.scanTranslation(n, _np_as(sel,c_int_p), ia0, ia1, l, nstep, _np_as(Es,c_double_p), bWriteTrj)

#  void scanRotation_ax( int n, int* sel, double* p0, double* ax, double phi, int nstep, double* Es, bool bWriteTrj ){
lib.scanRotation_ax.argtypes  = [c_int, c_int_p, c_double_p, c_double_p, c_double, c_int, c_double_p, c_bool] 
lib.scanRotation_ax.restype   =  None
def scanRotation_ax( phi, nstep, sel=0, p0=None, ax=None,  Es=None, bWriteTrj=False):
    n=0
    if sel is not None:
        n=len(sel)
        sel=np.array(sel,np.int32)
    if p0 is not None: p0=np.array(p0)
    if ax is not None: ax=np.array(ax)
    lib.scanRotation_ax(n, _np_as(sel,c_int_p), _np_as(p0,c_double_p), _np_as(ax,c_double_p), phi, nstep, _np_as(Es,c_double_p), bWriteTrj)
    return Es

#  void scanRotation( int n, int* sel,int ia0, int iax0, int iax1, double phi, int nstep, double* Es, bool bWriteTrj ){
lib.scanRotation.argtypes  = [c_int, c_int_p, c_int, c_int, c_int, c_double, c_int, c_double_p, c_bool] 
lib.scanRotation.restype   =  None
def scanRotation( ia0, iax0, iax1, phi, nstep, sel=None, Es=None, bWriteTrj=False, _0=0):
    n=0
    if _0!=0: ia0 -=_0; iax0-=_0; iax1-=_0;
    if sel is not None:
        n=len(sel)
        sel=np.array(sel,np.int32)
        sel-=_0
    if Es is None: Es=np.zeros(nstep)
    lib.scanRotation(n, _np_as(sel,c_int_p), ia0, iax0, iax1, phi, nstep, _np_as(Es,c_double_p), bWriteTrj)
    return Es

def scanBondRotation( ib, phi, nstep, Es=None, bWriteTrj=False, bPrintSel=False):
    nsel = splitAtBond(ib) 
    if bPrintSel: print( "split to:\n", selection[:nsel],"\n", selection[nsel:] )
    ias = bond2atom[ib,:]
    return scanRotation( ias[0], ias[0], ias[1], phi, nstep, sel=None, Es=Es, bWriteTrj=bWriteTrj)


#lib.scan.argtypes  = [c_int, array2d, array2d, array1d, array2d, array2d, c_bool, c_bool, c_int, c_double, c_double, c_double]
# int  scan( int nconf, double* poss, double* rots, double* Es, double* aforces, double* aposs, bool omp, bool bRelax, int niter_max, double dt, double Fconv, double Flim ){
lib.scan.argtypes = [ c_int, c_double_p, c_double_p, c_double_p, c_double_p, c_double_p, c_double_p, c_bool, c_bool, c_int, c_double, c_double, c_double ]
lib.scan.restype   =  None
def scan(poss, rots=None, dirs=None,  Es=None, aforces=None, aposs=None,  bF=False,bP=False, omp=False, bRelax=False, niter_max=10000, dt=0.05, Fconv=1e-5, Flim=100.0 ):
    nconf=len(poss)
    if Es is None: Es=np.zeros(nconf)
    if (aforces is None) and bF: aforces=np.zeros( (nconf,natoms,3) )
    if (aposs is None) and bP:   aposs=np.zeros(   (nconf,natoms,3) )
    #lib.scan(nconf, poss, rots, Es, aforces, aposs, omp, bRelax, niter_max, dt, Fconv, Flim )
    lib.scan( nconf, _np_as(poss,c_double_p), _np_as(rots,c_double_p), _np_as(dirs,c_double_p), _np_as(Es,c_double_p), _np_as(aforces,c_double_p), _np_as(aposs,c_double_p), omp, bRelax, niter_max, dt, Fconv, Flim )
    return Es, aforces, aposs 


# ====================================
# ========= Python Functions
# ====================================

def plot(b2as=None,ax1=0,ax2=1,ps=None, fs=None, bForce=False, Fscale=1.0 ):
    import matplotlib.pyplot as plt
    from matplotlib import collections  as mc
    if ps is None: ps=apos
    if fs is None: fs=fapos
    plt.plot( ps[:,ax1], ps[:,ax2],'o', c='#8080FF' )
    if(bForce): plt.quiver( ps[:,ax1], ps[:,ax2], fs[:,ax1], fs[:,ax2], scale = 1/Fscale )
    if glob_bMMFF:
        if b2as  is None: b2as=bond2atom
        lines = [  ((ps[b[0],ax1],ps[b[0],ax2]),(ps[b[1],ax1],ps[b[1],ax2])) for b in b2as ]
        lc = mc.LineCollection(lines, colors='#808080', linewidths=2)
        ax=plt.gca()
        ax.add_collection(lc)

        for i,p in enumerate(ps):    ax.annotate( "%i"%i , (p[ax1],p[ax2]), color='b' )
        for i,l in enumerate(lines): 
            p= ((l[0][0]+l[1][0])*0.5,(l[0][1]+l[1][1])*0.5)
            ax.annotate( "%i"%i , p, color='k')
    plt.axis('equal')
    plt.grid()

def plot_selection(sel=None,ax1=0,ax2=1,ps=None, s=100):
    import matplotlib.pyplot as plt
    if sel is None: sel=selection
    if ps is None: ps=apos
    asel=ps[sel]
    plt.scatter( asel[:,ax1], asel[:,ax2], s=s, facecolors='none', edgecolors='r' )



# ====================================
# ========= Test Functions
# ====================================


# ====================================
# ========= Python Functions
# ====================================





# ====================================
# ========= MAIN
# ====================================

# if __name__ == "__main__":
#     import matplotlib.pyplot as plt
#     plt.show()