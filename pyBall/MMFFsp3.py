
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
#"void init_buffers()",
#"void init_params(char* fatomtypes, char* fbondtypes)",
#"void init_nonbond()",
#"void initParams       ( const char* sAtomTypes, const char* sBondTypes, const char* sAngleTypes )",
#"int  buildMolecule_xyz( const char* xyz_name )",
#"void makeMMFF         (                      )",
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
#"void init( char* xyz_name, char* surf_name, char* smile_name, bool bMMFF=false, int* nPBC, double gridStep, char* sAtomTypes, char* sBondTypes, char* sAngleTypes ){",
#"void scanRotation_ax( int n, int* selection, double* p0, double* ax, double phi, int nstep, double* Es, bool bWriteTrj )",
#"void scanRotation( int n, int* selection,int ia0, int iax0, int iax1, double phi, int nstep, double* Es, bool bWriteTrj )",
#"void scanAngleToAxis_ax( int n, int* selection, double r, double R, double* p0, double* ax, int nstep, double* angs, double* Es, const char* trjName )",
#"void setOptLog( int n, double* cos, double* f, double* v, double* dt, double* damp ){",
#"void printAtomConfs ( bool bOmmitCap, bool bNeighs )",
#"void printAtomTypes ( )",
#"void printBonds     ( )",
#"void printBondParams( )",
#"int saveXYZ( const char* fname, const char* comment)",
#"int selectBondsBetweenTypes( int imin, int imax, int it1, int it2, bool byZ, bool bOnlyFirstNeigh, int* atoms_ ){"
#"int getFrament( int ifrag, int* bounds_, double* pose ){",
#"void scanHBond( const char* fname, int n, double d,  int ifrag1, int ifrag2, int i1a,int i1b, int i2a,int i2b ){",
#"void orient( int fw1,int fw2,  int up1,int up2,  int i0,  int imin, int imax ){",
#"void findMainAxes( double* rot, ifrag=-1, int imin=0,int imax=-1, bool bRot=true){",
#"void findSymmetry( int* found, int i0=0,int imax=-1, double tol=0.1 ){",
#"double measureBond(int ia, int ib          )",
#"double measureAngle(int ic, int ia, int ib )",
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
lib = cpp_utils.loadLib('MMFFsp3_lib', recompile=False)
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

# #  void sample_DistConstr( double lmin, double lmax, double kmin, double kmax, double flim , int n, double* xs, double* Es, double* Fs ){
# lib.sample_DistConstr.argtypes  = [c_double, c_double, c_double, c_double, c_double, c_int, c_double_p, c_double_p, c_double_p] 
# lib.sample_DistConstr.restype   =  None
# def sample_DistConstr( xs, lmin=1, lmax=1, kmin=1, kmax=1, flim=1e+300, Es=None, Fs=None):
#     n = len(xs)
#     if Es is None: Es=np.zeros(n)
#     if Fs is None: Fs=np.zeros(n)
#     lib.sample_DistConstr(lmin, lmax, kmin, kmax, flim, n, _np_as(xs,c_double_p), _np_as(Es,c_double_p), _np_as(Fs,c_double_p))
#     return Es,Fs

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
    if fs is None: fs=np.zeros(n)
    rs=np.array(rs)
    REQi=np.array(REQi)
    REQj=np.array(REQj) 
    lib.sampleNonBond(n, rs, Es, fs, kind, REQi, REQj, K, Rdamp)
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
    ptr = lib.getIBuff( cstr(name) )
    #print(name," : ",ptr)
    return np.ctypeslib.as_array( ptr, shape=sh)

#double* getBuff(const char* name){ 
lib.getBuff.argtypes = [c_char_p]
lib.getBuff.restype  = c_double_p 
def getBuff(name,sh):
    if not isinstance(sh, tuple): sh=(sh,)
    ptr = lib.getBuff( cstr(name) )
    #print(name," : ",ptr)
    return np.ctypeslib.as_array( ptr, shape=sh)

def getBuffs( ):
    init_buffers()
    #natom=nnode+ncap
    #nvecs=natom+npi
    #nDOFs=nvecs*3
    global ndims,Es
    ndims = getIBuff( "ndims", (5,) )  # [nDOFs,natoms,nnode,ncap,npi,nbonds]
    Es    = getBuff ( "Es",    (6,) )  # [ Etot,Eb,Ea, Eps,EppT,EppI; ]
    global nDOFs,natoms,nnode,ncap,nvecs
    #nDOFs=0,natoms=0,nnode=0,ncap=0,nvecs=0;
    nDOFs=ndims[0]; nnode=ndims[1]; ncap=ndims[2];nvecs=ndims[3]; natoms=nnode+ncap;
    print( "getBuffs(): natoms %i nnode %i ncap %i nvecs %i"  %(natoms,nnode,ncap,nvecs) )
    global DOFs,fDOFs,apos,fapos,pipos,fpipos, neighs #,selection
    #Ebuf     = getEnergyTerms( )
    apos      = getBuff ( "apos",     (natoms,3) )
    fapos     = getBuff ( "fapos",    (natoms,3) )
    if glob_bMMFF:
        DOFs      = getBuff ( "DOFs",     (nvecs,3)  )
        fDOFs     = getBuff ( "fDOFs",    (nvecs,3)  ) 
        pipos     = getBuff ( "pipos",    (nnode,3)  )
        fpipos    = getBuff ( "fpipos",   (nnode,3)  )
        neighs   = getIBuff( "neighs",  (natoms,4) )
        #selection = getIBuff( "selection",  (natoms) )

#  void init_buffers()
lib.init_buffers.argtypes  = []
lib.init_buffers.restype   =  None
def init_buffers():
    return lib.init_buffers()

#  void init_nonbond()
#lib.init.argtypes  = [] 
#lib.init.restype   =  None
#def init():
#    isInitialized = True
#    lib.init()

def cstr( s ):
    if s is None: return None
    return s.encode('utf8')

#  void setVerbosity( int verbosity_, int idebug_ ){
lib.setVerbosity.argtypes  = [c_int, c_int] 
lib.setVerbosity.restype   =  None
def setVerbosity( verbosity=1, idebug=0 ):
    return lib.setVerbosity( verbosity, idebug )

# void* init( char* xyz_name, char* smile_name, int* nPBC, char* sAtomTypes, char* sBondTypes, char* sAngleTypes, const char* sDihedralTypes ){
lib.init.argtypes  = [c_char_p, c_char_p, array1i, c_char_p, c_char_p, c_char_p, c_char_p, c_char_p ] 
lib.init.restype   =  c_void_p
def init(
        xyz_name  ="input.xyz", 
        smile_name=None, 
        sElementTypes = "common_resources/ElementTypes.dat",
        sAtomTypes = "common_resources/AtomTypes.dat", 
        sBondTypes = "common_resources/BondTypes.dat", 
        sAngleTypes = "common_resources/AngleTypes.dat",
        sDihedralTypes = "common_resources/DihedralTypes.dat",
        nPBC=(1,1,0)
    ):
    nPBC=np.array(nPBC,dtype=np.int32)
    return lib.init( cstr(xyz_name),  cstr(smile_name), nPBC, cstr(sElementTypes), cstr(sAtomTypes), cstr(sBondTypes), cstr(sAngleTypes), cstr(sDihedralTypes) )

def tryInit():
    if not isInitialized:
        init()
        #getBuffs( NEIGH_MAX=4 )

#  void initParams       ( const char* sAtomTypes, const char* sBondTypes, const char* sAngleTypes )
lib.initParams.argtypes  = [c_char_p, c_char_p, c_char_p, c_char_p , c_char_p ] 
lib.initParams.restype   =  None
def initParams( sElementTypes="common_resources/ElementTypes.dat", sAtomTypes="common_resources/AtomTypes.dat", sBondTypes = "common_resources/BondTypes.dat", sAngleTypes= "common_resources/AngleTypes.dat", sDihedralTypes = "common_resources/DihedralTypes.dat",  ):
    return lib.initParams( cstr(sElementTypes), cstr(sAtomTypes), cstr(sBondTypes), cstr(sAngleTypes),  cstr(sDihedralTypes) )

#  int  buildMolecule_xyz( const char* xyz_name )
lib. buildMolecule_xyz.argtypes  = [c_char_p, c_bool, c_double, c_bool, c_bool ] 
lib. buildMolecule_xyz.restype   =  c_int
def  buildMolecule_xyz(xyz_name, bEpairs=False, fAutoCharges=-1, bAutoTypes=True, bRelaxPi=False ):
    return lib.buildMolecule_xyz( cstr(xyz_name), bEpairs, fAutoCharges, bAutoTypes, bRelaxPi )

#  void makeFFs         (                      )
lib.makeFFs.argtypes  = [] 
lib.makeFFs.restype   =  None
def makeFFs():
    return lib.makeFFs()

#  void clear         (                      )
lib.clear.argtypes  = [] 
lib.clear.restype   =  None
def clear():
    return lib.clear()

#  void insertSMILES(char* s)
lib.insertSMILES.argtypes  = [c_char_p,c_bool] 
lib.insertSMILES.restype   =  None
def insertSMILES(s ):
    s = s.encode('utf8')
    return lib.insertSMILES(s)

# #  void initWithMolFile(char* fname_mol, bool bNonBonded_, bool bOptimizer_ )
# lib.initWithSMILES.argtypes  = [c_char_p, c_bool, c_bool, c_bool, c_bool] 
# lib.initWithSMILES.restype   =  None
# def initWithSMILES(fname_mol, bPrint=True, bCap=True, bNonBonded=True, bOptimizer=True):
#     fname_mol = fname_mol.encode('utf8')
#     return lib.initWithSMILES(fname_mol, bPrint, bCap, bNonBonded, bOptimizer)

#  void setSwitches( int CheckInvariants, int PBC, int NonBonded, int MMFF, int Angles, int PiSigma, int PiPiI  ){
lib.setSwitches.argtypes  = [c_int, c_int, c_int, c_int, c_int, c_int, c_int] 
lib.setSwitches.restype   =  None
def setSwitches( CheckInvariants=0, PBC=0, NonBonded=0, MMFF=0, Angles=0, PiSigma=0, PiPiI=0 ):
    return lib.setSwitches( CheckInvariants, PBC, NonBonded, MMFF, Angles, PiSigma, PiPiI )

lib.printSwitches.argtypes  = [] 
lib.printSwitches.restype   =  None
def printSwitches():
    return lib.printSwitches()

#  bool checkInvariants( double maxVcog, double maxFcog, double maxTg )
lib.checkInvariants.argtypes  = [c_double, c_double, c_double] 
lib.checkInvariants.restype   =  c_bool
def checkInvariants(maxVcog=1e-9, maxFcog=1e-9, maxTg=1e-1):
    return lib.checkInvariants(maxVcog, maxFcog, maxTg)

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
lib. run.argtypes  = [c_int, c_double, c_double, c_int, c_double_p, c_double_p] 
lib. run.restype   =  c_int
def  run(nstepMax=1000, dt=-1, Fconv=1e-6, ialg=2, outE=None, outF=None):
    return lib.run(nstepMax, dt, Fconv, ialg, _np_as(outE,c_double_p), _np_as(outF,c_double_p) )

#bool relax( int niter, double Ftol )
#lib.relax.argtypes  = [c_int, c_double, c_bool ] 
#lib.relax.restype   =  c_bool
#def relax(niter=100, Ftol=1e-6, bWriteTrj=False ):
#    return lib.relax(niter, Ftol, bWriteTrj )

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

#  int splitAtBond( int ib, int* sel )
# lib.splitAtBond.argtypes  = [c_int, c_int_p] 
# lib.splitAtBond.restype   =  c_int
# def splitAtBond(ib, sel=None):
#     return lib.splitAtBond(ib, _np_as(sel,c_int_p))


#  double measureBond(int ia, int ib          )
lib.measureBond.argtypes  = [c_int, c_int] 
lib.measureBond.restype   =  c_double
def measureBond(ia, ib):
    return lib.measureBond(ia, ib)

#  double measureAngle(int ic, int ia, int ib )
lib.measureAngle.argtypes  = [c_int, c_int, c_int] 
lib.measureAngle.restype   =  c_double
def measureAngle(ic, ia, ib):
    return lib.measureAngle(ic, ia, ib)

#  double measureAnglePiPi(int ia, int ib          )
lib.measureAnglePiPi.argtypes  = [c_int, c_int] 
lib.measureAnglePiPi.restype   =  c_double
def measureAnglePiPi(ia, ib):
    return lib.measureAnglePiPi(ia, ib)

#  double measureAngleSigmaPi(int ic, int ia, int ib )
lib.measureAngleSigmaPi.argtypes  = [c_int, c_int, c_int] 
lib.measureAngleSigmaPi.restype   =  c_double
def measureAngleSigmaPi(ipi, ia, ib):
    return lib.measureAngleSigmaPi(ipi, ia, ib)

#  void scanTranslation_ax( int n, int* sel, double* vec, int nstep, double* Es, bool bWriteTrj )
lib.scanTranslation_ax.argtypes  = [c_int, c_int_p, c_double_p, c_int, c_double_p, c_char_p, c_bool] 
lib.scanTranslation_ax.restype   =  None
def scanTranslation_ax(nstep, Es=None, trjName=None, sel=None, vec=None, bAddjustCaps=False ):
    n=0
    if sel is not None:
        n=len(sel)
        sel=np.array(sel,np.int32)
    if vec is not None: vec=np.array(vec)
    return lib.scanTranslation_ax(n, _np_as(sel,c_int_p), _np_as(vec,c_double_p), nstep, _np_as(Es,c_double_p), cstr(trjName) )

#  void scanTranslation( int n, int* sel, int ia0, int ia1, double l, int nstep, double* Es, bool bWriteTrj )
lib.scanTranslation.argtypes  = [c_int, c_int_p, c_int, c_int, c_double, c_int, c_double_p, c_char_p, c_bool ] 
lib.scanTranslation.restype   =  None
def scanTranslation( ia0, ia1, l, nstep, Es=None, sel=None, trjName=None, ):
    if Es is None: Es=np.zeros(nstep)
    return lib.scanTranslation(n, _np_as(sel,c_int_p), ia0, ia1, l, nstep, _np_as(Es,c_double_p), cstr(trjName), bAddjustCaps=False )

#  void scanRotation_ax( int n, int* sel, double* p0, double* ax, double phi, int nstep, double* Es, bool bWriteTrj ){
lib.scanRotation_ax.argtypes  = [c_int, c_int_p, c_double_p, c_double_p, c_double, c_int, c_double_p, c_char_p ] 
lib.scanRotation_ax.restype   =  None
def scanRotation_ax( phi, nstep, sel=0, p0=None, ax=None,  Es=None, trjName=None, bEs=True):
    n=0
    if sel is not None:
        n=len(sel)
        sel=np.array(sel,np.int32)
    if p0 is not None: p0=np.array(p0)
    if ax is not None: ax=np.array(ax)
    if bEs and ( Es is None ): Es = np.zeros(nstep)
    lib.scanRotation_ax(n, _np_as(sel,c_int_p), _np_as(p0,c_double_p), _np_as(ax,c_double_p), phi, nstep, _np_as(Es,c_double_p), cstr(trjName) )
    return Es

#  void scanRotation( int n, int* sel,int ia0, int iax0, int iax1, double phi, int nstep, double* Es, bool bWriteTrj ){
lib.scanRotation.argtypes  = [c_int, c_int_p, c_int, c_int, c_int, c_double, c_int, c_double_p, c_char_p ] 
lib.scanRotation.restype   =  None
def scanRotation( ia0, iax0, iax1, phi, nstep, sel=None, Es=None, trjName=None,  _0=0):
    n=0
    if _0!=0: ia0 -=_0; iax0-=_0; iax1-=_0;
    if sel is not None:
        n=len(sel)
        sel=np.array(sel,np.int32)
        sel-=_0
    if Es is None: Es=np.zeros(nstep)
    lib.scanRotation(n, _np_as(sel,c_int_p), ia0, iax0, iax1, phi, nstep, _np_as(Es,c_double_p), str(trjName)  )
    return Es

#  void scanAngleToAxis_ax( int n, int* selection, double r, double R, double* p0, double* ax, int nstep, double* angs, double* Es, const char* trjName )
lib.scanAngleToAxis_ax.argtypes  = [c_int, c_int_p, c_double, c_double, c_double_p, c_double_p, c_int, c_double_p, c_double_p, c_char_p] 
lib.scanAngleToAxis_ax.restype   =  None
def scanAngleToAxis_ax(angs, sel, r=1.0, p0=[0.,0.,0.], ax=[0.,0.,1.], R=0, Es=None, trjName=None, bEs=True):
    n=0
    if sel is not None:
        n=len(sel)
        sel=np.array(sel,np.int32)
    if p0 is not None: p0=np.array(p0)
    if ax is not None: ax=np.array(ax)
    nstep = len(angs)
    if bEs and ( Es is None ): Es = np.zeros(nstep)
    lib.scanAngleToAxis_ax(n, _np_as(sel,c_int_p), r, R, _np_as(p0,c_double_p), _np_as(ax,c_double_p), nstep, _np_as(angs,c_double_p), _np_as(Es,c_double_p), cstr(trjName)  )
    return Es

def scanBondRotation( ib, phi, nstep, Es=None, bWriteTrj=False, bPrintSel=False):
    nsel = splitAtBond(ib) 
    if bPrintSel: print( "split to:\n", selection[:nsel],"\n", selection[nsel:] )
    ias = bond2atom[ib,:]
    return scanRotation( ias[0], ias[0], ias[1], phi, nstep, sel=None, Es=Es, bWriteTrj=bWriteTrj)

#  int selectBondsBetweenTypes( int imin, int imax, int it1, int it2, bool byZ, bool bOnlyFirstNeigh, int* atoms_ ){
lib.selectBondsBetweenTypes.argtypes  = [c_int, c_int, c_int, c_int, c_bool, c_bool, c_int_p ] 
lib.selectBondsBetweenTypes.restype   =  c_int
def selectBondsBetweenTypes( imin, imax, it1, it2, byZ=False, bOnlyFirstNeigh=False, atoms=None ):
    if(atoms is None): atoms=np.zeros((100,2), dtype=np.int32 )
    n = lib.selectBondsBetweenTypes(imin, imax, it1, it2, byZ, bOnlyFirstNeigh, _np_as(atoms ,c_int_p) )
    return atoms[:n,:]

#  void findMainAxes( double* rot, ifrag=-1, int imin=0,int imax=-1, bool bRot=true){
lib.findMainAxes.argtypes  = [c_double_p, c_int, c_int, c_int, c_int_p, c_bool] 
lib.findMainAxes.restype   =  None
def findMainAxes(ifrag=-1, imin=0, imax=-1, rot=None, bRot=True, bSaveRot=True, permut=None ):
    if(permut is not None): permut=np.array(permut,dtype=np.int23)
    if( bSaveRot and (rot is None) ): rot=np.zeros((3,3))
    lib.findMainAxes( _np_as(rot,c_double_p), ifrag, imin, imax, permut, bRot )
    return rot

#  void findSymmetry( int* found, int i0=0,int imax=-1, double tol=0.1 ){
lib.findSymmetry.argtypes  = [c_int_p, c_int, c_int, c_int, c_double] 
lib.findSymmetry.restype   =  None
def findSymmetry(found=None, ifrag=-1, i0=0, imax=-1, tol=0.1):
    if(found is None): found=np.zeros(4, dtype=np.int32)
    lib.findSymmetry(_np_as(found,c_int_p), ifrag, i0, imax, tol)
    return found


#  void orient( int fw1,int fw2,  int up1,int up2,  int i0,  int imin, int imax ){
lib.orient.argtypes  = [c_char_p, c_int, c_int, c_int, c_int, c_int, c_int, c_int] 
lib.orient.restype   =  None
def orient(fw, up, i0=None, imin=0, imax=-1, fname="oriented.xyz" ):
    if i0 is None: i0=fw[0]
    return lib.orient( cstr(fname), fw[0],fw[1], up[0],up[1], i0, imin, imax)

#  void scanHBond( const char* fname, int n, double d,  int ifrag1, int ifrag2, int i1a,int i1b, int i2a,int i2b ){
lib.scanHBond.argtypes  = [c_char_p, c_int, c_double, c_double, c_int, c_int, c_int, c_int, c_int, c_int, c_bool,c_bool, c_double_p ] 
lib.scanHBond.restype   =  None
def scanHBond( b1, b2, ifrag1=0, ifrag2=1, fname="scanHBond.xyz", n=10, dl=0.2, l0=5.0, isDonor=(True,False), ups=None ):
    if(ups is not None): ups=np.array(ups)
    return lib.scanHBond( cstr(fname), n, dl, l0, ifrag1, ifrag2, b1[0],b1[1], b2[0],b2[1],  isDonor[0],isDonor[1],  _np_as(ups,c_double_p)  )

#  int getFrament( int ifrag, int* bounds_, double* pose ){
lib.getFrament.argtypes  = [c_int, c_int_p, c_double_p] 
lib.getFrament.restype   =  c_int
def getFrament(ifrag, bounds=None, pose=None, bBounds=True, bPose=False):
    if( bBounds and (bounds is None) ): bounds=np.zeros((3,2), dtype=np.int32 )
    if( bPose   and (bPose  is None) ): pose  =np.zeros((4,3))
    imol = lib.getFrament(ifrag, _np_as(bounds,c_int_p), _np_as(pose,c_double_p))
    return bounds, pose, imol

#  void printTypes ( )
lib.printTypes.argtypes  = [] 
lib.printTypes.restype   =  None
def printTypes():
    lib.printTypes()

#  void printAtomConfs ( bool bOmmitCap, bool bNeighs )
lib.printAtomConfs.argtypes  = [c_bool, c_bool] 
lib.printAtomConfs.restype   =  None
def printAtomConfs(bOmmitCap=False, bNeighs=True):
    lib.printAtomConfs(bOmmitCap, bNeighs)

#  void printAtomTypes ( )
lib.printAtomTypes.argtypes  = [] 
lib.printAtomTypes.restype   =  None
def printAtomTypes():
    lib.printAtomTypes()

#  void printBonds     ( )
lib.printBonds.argtypes  = [] 
lib.printBonds.restype   =  None
def printBonds():
    lib.printBonds()

#  void printBondParams( )
lib.printBondParams.argtypes  = [] 
lib.printBondParams.restype   =  None
def printBondParams():
    lib.printBondParams()

#  void printAtomParams( )
lib.printAtomParams.argtypes  = [] 
lib.printAtomParams.restype   =  None
def printAtomParams():
    lib.printAtomParams()

# ====================================
# ========= Python Functions
# ====================================

def scanAllHBonds( path1, path2, t1=8,t2=8, ups=None, bOnlyFirstNeigh1=True, bOnlyFirstNeigh2=True ):
    clear()
    f1 = buildMolecule_xyz( xyz_name=path1, bEpairs=True )
    f2 = buildMolecule_xyz( xyz_name=path2, bEpairs=True )
    b1,_,_ = getFrament(f1);                                                        print( "bonds_1", b1[2] )
    b2,_,_ = getFrament(f2);                                                        print( "bonds_2", b2[2] )
    donors_1    = selectBondsBetweenTypes( b1[2,0], b1[2,1], t1, 1,   True, bOnlyFirstNeigh1 ); print( "donors_1    ", donors_1    )
    donors_2    = selectBondsBetweenTypes( b2[2,0], b2[2,1], t2, 1,   True, bOnlyFirstNeigh1 ); print( "donors_2    ", donors_2    )
    acceptors_1 = selectBondsBetweenTypes( b1[2,0], b1[2,1], t1, 200, True, bOnlyFirstNeigh2 ); print( "acceptors_1 ", acceptors_1 )
    acceptors_2 = selectBondsBetweenTypes( b2[2,0], b2[2,1], t2, 200, True, bOnlyFirstNeigh2 ); print( "acceptors_2 ", acceptors_2 )

    fname1 = os.path.split(path1)[1];  fname1=os.path.splitext(fname1)[0]
    fname2 = os.path.split(path2)[1];  fname2=os.path.splitext(fname2)[0]

    for i,b1 in enumerate(donors_1):
        for j,b2 in enumerate(acceptors_2):
            fname =  "scanHBond_"+fname1+"-H"+str(i)+"_vs_"+fname2+"-e"+str(j)+".xyz"
            print( "scanHBond ",  b1, b2, fname )
            scanHBond( b1, b2, l0=2.0, fname=fname, isDonor=(True,False), ups=ups )
    for i,b1 in enumerate(acceptors_1):
        for j,b2 in enumerate(donors_2):
            fname = "scanHBond_"+fname1+"-e"+str(i)+"_vs_"+fname2+"-H"+str(j)+".xyz"
            scanHBond( b1, b2, l0=2.0, fname=fname, isDonor=(False,True), ups=ups )
            print( "scanHBond ", b1, b2, fname )



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
# ========= MAIN
# ====================================

# if __name__ == "__main__":
#     import matplotlib.pyplot as plt
#     plt.show()