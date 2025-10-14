import numpy as np
from   ctypes import c_int, c_double, c_bool, c_float, c_char_p, c_bool, c_void_p, c_char_p, POINTER
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
c_bool_p   = ctypes.POINTER(c_bool)

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
#"int  sampleHbond( int ib, int n, double* rs, double* Es, double* fs, int kind, double K, double Rdamp ){",
#"#void sampleSurf(char* name, int n, double* rs, double* Es, double* fs, int kind, double*REQ_, double K, double Rdamp ){",
#"void init( char* xyz_name, char* surf_name, char* smile_name, bool bMMFF=false, int* nPBC, double gridStep, char* sAtomTypes, char* sBondTypes, char* sAngleTypes ){",
#"void scanRotation_ax( int n, int* selection, double* p0, double* ax, double phi, int nstep, double* Es, bool bWriteTrj )",
#"void scanRotation( int n, int* selection,int ia0, int iax0, int iax1, double phi, int nstep, double* Es, bool bWriteTrj )",
#"void set_opt( double dt_max,  double dt_min, double damp_max, double finc,    double fdec,   double falpha, int minLastNeg, double cvf_min, double cvf_max){",
#"void sample_evalAngleCos( double K, double c0, int n, double* angles, double* Es, double* Fs ){",

#void sample_SplineHermite( double x0, double dx, int np, double* Eps, int n, double* xs, double* Es, double* Fs )
#void sample_SplineHermite2D( double* g0, double* dg, int* ng, double* Eg, int n, double* ps, double* fes )
#void sample_SplineHermite3D( double* g0, double* dg, int* ng, double* Eg, int n, double* ps, double* fes )

# void sample_SplineConstr( double lmin, double lmax, double dx, int n, double* Eps, double* Es, double* Fs ){
#"void sample_DistConstr( double lmin, double lmax, double kmin, double kmax, double flim , int n, double* xs, double* Es, double* Fs ){",
#"void addDistConstrain(  int i0,int i1, double lmin,double lmax,double kmin,double kmax,double flim, double k ){",
#"void addAngConstrain(  int i0,int i1,int i2, double ang0, double k ){",
#"int substituteMolecule( const char* fname, int ib, double* up, int ipivot, bool bSwapBond )",
#'void print_debugs( bool bParams, bool bNeighs, bool bShifts ){',
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
lib = cpp_utils.loadLib('MMFF_lib', recompile=False)






array1ui = np.ctypeslib.ndpointer(dtype=np.uint32, ndim=1, flags='CONTIGUOUS')
array1i  = np.ctypeslib.ndpointer(dtype=np.int32,  ndim=1, flags='CONTIGUOUS')
array2i  = np.ctypeslib.ndpointer(dtype=np.int32,  ndim=2, flags='CONTIGUOUS')
array1d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS')
array2d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=2, flags='CONTIGUOUS')
array3d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=3, flags='CONTIGUOUS')


class AtomType(ctypes.Structure):
    _fields_ = [
        ("name",      ctypes.c_char * 8),
        ("iZ",        ctypes.c_uint8),
        ("valence",   ctypes.c_uint8),
        ("nepair",    ctypes.c_uint8),
        ("npi",       ctypes.c_uint8),
        ("sym",       ctypes.c_uint8),
        ("color",     ctypes.c_uint32),
        ("Ruff",      ctypes.c_double),
        ("RvdW",      ctypes.c_double),
        ("EvdW",      ctypes.c_double),
        ("Qbase",     ctypes.c_double),
        ("Hb",        ctypes.c_double),
        ("parrent",   ctypes.c_int),
        ("element",   ctypes.c_int),
        ("ePairType", ctypes.c_int),
        ("bMMFF",     ctypes.c_bool),
        ("Ass",       ctypes.c_double),
        ("Asp",       ctypes.c_double),
        ("Kss",       ctypes.c_double),
        ("Ksp",       ctypes.c_double),
        ("Kep",       ctypes.c_double),
        ("Kpp",       ctypes.c_double),
        ("subTyp_x",  ctypes.c_int),
        ("subTyp_y",  ctypes.c_int),
        ("subTyp_z",  ctypes.c_int),
        #Vec3i     subTypes=Vec3iZero;  // sp1 sp2 sp3    // Q1 Q2 Q3 (polarized)
    ]

p_AtomType = ctypes.POINTER(AtomType)

lib.getAtomTypes.restype     = p_AtomType
lib.getAtomTypeCount.restype = c_int

def getAtomTypes():
    ptr = lib.getAtomTypes()
    n   = lib.getAtomTypeCount()
    array_type = AtomType * n
    arr = ctypes.cast(ptr, ctypes.POINTER(array_type)).contents
    return arr, n

# ====================================
# ========= Globals
# ====================================

isInitialized       = False
glob_bMMFF          = True
glob_bUFF           = False
bBuffersInitialized = False

# ====================================
# ========= C functions
# ====================================


# void samplePBCindexes( int n, int* inds, int ng, int* iout, int order ){
lib.samplePBCindexes.argtypes  = [ c_int, array1i, c_int, array2i, c_int ]
lib.samplePBCindexes.restype   =  None
def samplePBCindexes( inds, ng, iout=None, order=3 ):
    n = len(inds)
    if iout is None: iout = np.zeros( (n,(order+1)), dtype=np.int32 )
    lib.samplePBCindexes( n, inds, ng, iout, order )
    return iout

# void projectBspline1D( int nx, double* xs, double* ws, double g0, double dg, int ng, double* ys, int order ){
lib.projectBspline1D.argtypes  = [ c_int, c_double_p, c_double_p, c_double, c_double, c_int, c_double_p, c_int ]
lib.projectBspline1D.restype   =  None
def projectBspline1D( xs, g0, dg, ng, ys=None, ws=None, order=3 ):
    n = len(xs)
    xs = np.array( xs )
    if ys is None: ys = np.zeros( ng )
    if ws is None: ws = np.ones( n )
    lib.projectBspline1D( n, _np_as(xs,c_double_p), _np_as(ws,c_double_p), g0, dg, ng, _np_as(ys,c_double_p), order )
    return ys

# void projectBspline2D( int nx, double* ps_, double* ws, double* g0_, double* dg_, int* ng_, double* ys, int order ){
lib.projectBspline2D.argtypes  = [ c_int, c_double_p, c_double_p, c_double_p, c_double_p, c_int_p, c_double_p, c_int ]
lib.projectBspline2D.restype   =  None
def projectBspline2D( xs, g0, dg, ng, ys=None, ws=None, order=3 ):
    n = len(xs)
    xs = np.array( xs )
    dg = np.array( dg )
    g0 = np.array( g0 )
    ng = np.array( ng, dtype=np.int32 )
    if ys is None:
        ys = np.zeros( ng )
    else:
        ys = np.array( ys )
    if ws is None:
        ws = np.ones( n )
    else:
        ws = np.array( ws )
    lib.projectBspline2D( n, _np_as(xs,c_double_p), _np_as(ws,c_double_p), _np_as(g0,c_double_p), _np_as(dg,c_double_p), _np_as(ng,c_int_p), _np_as(ys,c_double_p), order )
    return ys

#  int fit_Bspline( const int n, double* Gs, double* Es, double* Ws, double Ftol, int nmaxiter, double dt, double Kreg, bool bPBC, bool bHalf  ){
lib.fit_Bspline.argtypes  = [ c_int, c_double_p, c_double_p, c_double_p, c_double, c_int, c_double, c_double,  c_bool,c_bool  ]
lib.fit_Bspline.restype   =  None
def fit_Bspline( Es, Gs=None, Ws=None, Ftol=1e-6, nmaxiter=100, dt=0.1, Kreg=-1.0, bPBC=True, bHalf=False ):
    n = len(Es)
    if(bHalf):
        n = n//2
        if Ws is None: Ws = np.ones( n*2 )
        if Gs is None: Gs = Es[::2].copy()
    else:
        if Ws is None: Ws = np.ones( n )
        if Gs is None: Gs = Es.copy()
    lib.fit_Bspline( n, _np_as(Gs,c_double_p), _np_as(Es,c_double_p), _np_as(Ws,c_double_p), Ftol, nmaxiter, dt, Kreg, bPBC, bHalf )
    return Gs, Ws

# int fitEF_Bspline( const int n, const double* Gs, double* fes, double* Ws, double Ftol, int nmaxiter, double dt ){
lib.fitEF_Bspline.argtypes  = [ c_double,  c_int, c_double_p, c_double_p, c_double_p, c_double, c_int, c_double ]
lib.fitEF_Bspline.restype   =  None
def fitEF_Bspline( dg, Fes, Gs=None, Ws=None, Ftol=1e-6, nmaxiter=100, dt=0.1 ):
    n = len(Fes)
    if Ws is None: Ws = np.ones( (n,2) )
    if Gs is None: Gs = Fes[:,0].copy()
    lib.fitEF_Bspline( dg, n, _np_as(Gs,c_double_p), _np_as(Fes,c_double_p), _np_as(Ws,c_double_p), Ftol, nmaxiter, dt )
    return Gs, Ws

# int fit2D_Bspline( const int* ns, double* Gs, double* Es, double* Ws, double Ftol, int nmaxiter, double dt ){
lib.fit2D_Bspline.argtypes  = [ c_int_p, c_double_p, c_double_p, c_double_p, c_double, c_int, c_double, c_bool ]
lib.fit2D_Bspline.restype   =  None
def fit2D_Bspline( Es, Gs=None, Ws=None, Ftol=1e-6, nmaxiter=100, dt=0.1, bPBC=False ):
    ns = Es.shape
    if Ws is None: Ws = np.ones( ns )
    if Gs is None: Gs = Es.copy()
    ns = np.array( ns[::-1], dtype=np.int32 )
    lib.fit2D_Bspline( _np_as(ns,c_int_p) , _np_as(Gs,c_double_p), _np_as(Es,c_double_p), _np_as(Ws,c_double_p), Ftol, nmaxiter, dt, bPBC )
    return Gs, Ws

# int fit3D_Bspline( const int* ns, double* Gs, double* Es, double* Ws, double Ftol, int nmaxiter, double dt, bool bOMP ){
lib.fit3D_Bspline.argtypes  = [ c_int_p, c_double_p, c_double_p, c_double_p, c_double, c_int, c_double, c_bool, c_bool ]
lib.fit3D_Bspline.restype   =  None
def fit3D_Bspline( Es, Gs=None, Ws=None, Ftol=1e-6, nmaxiter=100, dt=0.1, bPBC=False, bOMP=True ):
    ns = np.array( Es.shape[::-1], dtype=np.int32 )
    n  = Es.size
    if Ws is None: Ws = np.ones( ns )
    if Gs is None: Gs = Es[:].copy()
    lib.fit3D_Bspline( _np_as(ns,c_int_p) , _np_as(Gs,c_double_p), _np_as(Es,c_double_p), _np_as(Ws,c_double_p), Ftol, nmaxiter, dt, bPBC, bOMP )
    return Gs

#void setupGrid( int* ns, double* cell, bool bAlloc ){
lib.setupGrid.argtypes  = [ c_int_p, c_double_p, c_bool ]
lib.setupGrid.restype   =  c_double_p
def setupGrid( ns, cell=None, bAlloc=True ):
    ns = np.array( (ns[2],ns[1],ns[0]) , dtype=np.int32 )
    if cell is not None: cell = np.array( cell )
    data_ptr =  lib.setupGrid( _np_as(ns,c_int_p), _np_as(cell,c_double_p), bAlloc )
    ff = None
    if( bAlloc ):
        ff   = np.ctypeslib.as_array(data_ptr, (ns[2],ns[1],ns[0]) )
    return ff

#void loadBin_d( const char* fname, double* data ){
lib.loadBin_f.argtypes  = [ c_char_p, c_float_p ]
lib.loadBin_f.restype   = c_int
def loadBin_f( name , data=None, ns=None):
    print( "ns ", ns )
    if data is None: data=np.zeros( ns[::-1], dtype=np.float32 )
    print( "data.shape ", data.shape )
    name=name.encode('utf8')
    ret = lib.loadBin_f( name, _np_as(data,c_float_p) )
    return data

#void loadBin_d( const char* fname, double* data ){
lib.loadBin_d.argtypes  = [ c_char_p, c_double_p ]
lib.loadBin_d.restype   = c_int
def loadBin_d( name , data=None, ns=None, chan=None ):
    print( "ns ", ns )
    ns[::-1]
    if chan is not None: ns = ns + (chan,)
    if data is None:     data=np.zeros( ns )
    print( "data.shape ", data.shape )
    name=name.encode('utf8')
    ret = lib.loadBin_d( name, _np_as(data,c_double_p) )
    return data

#void saveBin_d( const char* fname, double* data ){
lib.saveBin_d.argtypes  = [ c_char_p, c_double_p ]
lib.saveBin_d.restype   =  c_int
def saveBin_d( name, data):
    name=name.encode('utf8')
    return lib.saveBin_d( name, _np_as(data,c_double_p) )

# double* loadXSF( const char* fname, int* ns, double* cell ){
lib.loadXSF.argtypes  = [ c_char_p, c_int_p, c_double_p ]
lib.loadXSF.restype   =  c_double_p
def loadXSF( name ):
    if name is not None: name=name.encode('utf8')
    cell = np.zeros((3,3))
    ns   = np.zeros(3, dtype=np.int32 )
    data_ptr = lib.loadXSF( name, _np_as(ns,c_int_p), _np_as(cell,c_double_p)  )
    ff       = np.ctypeslib.as_array(data_ptr, (ns[2],ns[1],ns[0]) )
    return ff, cell

# double* loadXSF( const char* fname, int* ns, double* cell ){
lib.saveXSF.argtypes  = [ c_char_p, c_double_p, c_int_p, c_double_p ]
lib.saveXSF.restype   =  c_double_p
def saveXSF( name, FF, cell=None ):
    ns = FF.shape
    ns   = np.array( (ns[2],ns[1],ns[0]) , dtype=np.int32 )
    #ns   = np.array( ns, dtype=np.int32 )
    if cell is not None: cell = np.array( cell )
    if name is not None: name = name.encode('utf8')
    lib.saveXSF( name, _np_as(FF,c_double_p), _np_as(ns,c_int_p), _np_as(cell,c_double_p) )

# # New wrapper for saving grid plus geometry.
# # We assume that the C++ library now exposes a function called saveXSF_all
# # with signature similar to:
# #   int saveXSF_all(const char* fname, const double* FF, int pitch, int offset,
# #                     int natoms, int* atypes, Vec3d* apos, bool bPrimCoord);
# # (Vec3d is assumed to be represented as three doubles.)

# # import ctypes
# # from ctypes import c_int, c_double, c_bool, c_char_p, POINTER
# lib.saveXSF_all.argtypes = [ c_char_p, 
#                              ctypes.POINTER(c_double), 
#                              c_int, c_int, c_int,
#                              POINTER(c_int), 
#                              ctypes.POINTER(c_double), 
#                              c_bool ]
# lib.saveXSF_all.restype  = c_int

# def saveXSF_geometry(name, FF, pitch=1, offset=0, natoms=0, atypes=None, fapos=None, bPrimCoord=True):
#     """
#     Save grid data (FF) along with geometry to an XSF file.

#     Parameters:
#       name       : filename (str)
#       FF         : a numpy array with grid data
#       pitch      : grid pitch (int, default 1)
#       offset     : offset (int, default 0)
#       natoms     : number of atoms to write
#       atypes     : atom types array (list or numpy array of int)
#       fapos      : atomic positions (numpy array of shape (natoms,3))
#       bPrimCoord : bool flag; if True, positions are in the primitive coordinate system
#     """
#     # Prepare grid dimensions (reorder: [z,y,x]) for saving
#     ns = FF.shape
#     ns_arr = np.array((ns[2], ns[1], ns[0]), dtype=np.int32)
#     # For our function we don’t pass ns_arr explicitly because the C++ code
#     # likely derives it from the grid dimensions.
    
#     if atypes is not None:
#         atypes = np.array(atypes, dtype=np.int32)
#     if fapos is not None:
#         fapos = np.ascontiguousarray(fapos, dtype=np.float64)
        
#     return lib.saveXSF_all(name.encode('utf8'),
#                            _np_as(FF, ctypes.POINTER(c_double)),
#                            pitch, offset, natoms,
#                            _np_as(atypes, POINTER(c_int)) if atypes is not None else None,
#                            _np_as(fapos, ctypes.POINTER(c_double)) if fapos is not None else None,
#                            bPrimCoord)


#void makeGridFF( const char* name, int* ffshape, int mode, double z0, double* cel0, bool bSymmetrize, bool bAutoNPBC, bool bFit, bool bRefine ){
lib.makeGridFF.argtypes  = [ c_char_p, c_int_p, c_int, c_double, c_double_p, c_bool, c_bool, c_bool, c_bool ]
lib.makeGridFF.restype   =  None
def makeGridFF( name, mode=1, z0=np.nan, cel0=[-0.5,-0.5,0.0], bSymmetrize=False, bAutoNPBC=True, bFit=True, bRefine=True ):
    name=name.encode('utf8')
    cel0 = np.array( cel0 )
    ffshape = np.zeros( 4, dtype=np.int32 )   #;print( "ffshape ", ffshape )
    #print("mmff.makeGridFF: ", " bSymmetrize=", bSymmetrize, " bAutoNPBC=", bAutoNPBC, " bFit=", bFit  )
    lib.makeGridFF( name,  _np_as(ffshape,c_int_p), mode, float(z0), _np_as(cel0,c_double_p), bSymmetrize, bAutoNPBC, bFit, bRefine )
    #ffshape = ffshape[::-1]
    #print( "ffshape ", ffshape )
    #ff_ = np.ctypeslib.as_array(ff, ffshape )
    #print( "makeGridFF() DONE" )
    #return ff_

#double* getArrayPointer( const char* name, int* shape  ){
lib.getArrayPointer.argtypes  = [ c_char_p, c_int_p ]
lib.getArrayPointer.restype   =  c_double_p
def getArrayPointer( name ):
    name=name.encode('utf8')
    shape = np.zeros( 4, dtype=np.int32 )
    ptr = lib.getArrayPointer( name,  _np_as(shape,c_int_p) )
    if ptr:
        valid_shape = [dim for dim in shape if dim > 0]
        return np.ctypeslib.as_array( ptr, valid_shape )
    return None


# int setupEwaldGrid( double* pos0, double* dCell, int* ns, bool bPrint ){
lib.setupEwaldGrid.argtypes  = [ c_double_p, c_double_p, c_int_p, c_bool ]
lib.setupEwaldGrid.restype   =  c_int
def setupEwaldGrid( ns, pos0=[0.0,0.0,0.0], dCell=None, dg=[0.1,0.1,0.1], bPrint=False ):
    if dCell is None: dCell = [[dg[0],0.,0.],[0.,dg[1],0.],[0.,0.,dg[2]]]
    pos0  = np.array( pos0 )
    dCell = np.array( dCell )
    ns = np.array( ns, dtype=np.int32 )
    return lib.setupEwaldGrid( _np_as(pos0,c_double_p), _np_as(dCell,c_double_p), _np_as(ns,c_int_p), bPrint )


#void projectAtomsEwaldGrid( int na, double* apos, double* qs, double* dens, int order ){
lib.projectAtomsEwaldGrid.argtypes  = [ c_int, c_double_p, c_double_p, c_double_p, c_int ]
lib.projectAtomsEwaldGrid.restype   =  None
def projectAtomsEwaldGrid( apos, qs, dens=None, ns=None, order=2 ):
    na = len(apos)
    apos = np.array( apos )
    qs   = np.array( qs )
    if dens is None: dens = np.zeros( ns[::-1], dtype=np.float64 )
    lib.projectAtomsEwaldGrid( na, _np_as(apos,c_double_p), _np_as(qs,c_double_p), _np_as(dens,c_double_p), order )
    return dens


#void EwaldGridSolveLaplace( double* dens, int nz, double* Vout, bool bPrepare, bool bDestroy, int flags, bool bOMP, int nBlur, double cSOR, double cV ){
lib.EwaldGridSolveLaplace.argtypes  = [ c_double_p, c_int, c_double_p, c_bool, c_bool,  c_int, c_bool, c_int, c_double, c_double ]
lib.EwaldGridSolveLaplace.restype   =  None
def EwaldGridSolveLaplace( dens, nz_slab=-1, Vout=None, bPrepare=True, bDestroy=True, flags=-1, bOMP=False, nBlur=0, cSOR=0.0, cV=0.5 ):
    if Vout is None: Vout = np.zeros( dens.shape, dtype=np.float64 )
    lib.EwaldGridSolveLaplace( _np_as(dens,c_double_p), nz_slab, _np_as(Vout,c_double_p), bPrepare, bDestroy, flags, bOMP, nBlur, cSOR, cV )
    return Vout

# void EwaldGridSolveLaplaceDebug( double* dens, double* Vout, double* densw, double* kerw, double* VwKer ){
lib.EwaldGridSolveLaplaceDebug.argtypes  = [ c_double_p, c_double_p, c_double_p, c_double_p, c_double_p ]
lib.EwaldGridSolveLaplaceDebug.restype   =  None
def EwaldGridSolveLaplaceDebug( dens, Vout=None, densw=None, kerw=None, VwKer=None ):
    if Vout  is None: Vout  = np.zeros( dens.shape, dtype=np.float64 )
    if densw is None: densw = np.zeros( dens.shape, dtype=np.float64 )
    if kerw  is None: kerw  = np.zeros( dens.shape, dtype=np.float64 )
    if VwKer is None: VwKer = np.zeros( dens.shape, dtype=np.float64 )
    lib.EwaldGridSolveLaplaceDebug( _np_as(dens,c_double_p), _np_as(Vout,c_double_p), _np_as(densw,c_double_p), _np_as(kerw,c_double_p), _np_as(VwKer,c_double_p) )
    return Vout, densw, kerw, VwKer

#void projectMultiPole( double* p0, int n, double* ps, double* Qs, int order, double* cs ){
lib.projectMultiPole.argtypes  = [ c_double_p, c_int, c_double_p, c_double_p, c_int, c_double_p ]
lib.projectMultiPole.restype   =  None
def projectMultiPole( ps, Qs, order=2, p0=None, cs=None ):
    n = len(ps)
    ps = np.array( ps )
    Qs = np.array( Qs )
    if p0 is None: p0 = np.zeros( 3, dtype=np.float64 )
    if cs is None: cs = np.zeros( n, dtype=np.float64 )
    lib.projectMultiPole( _np_as(p0,c_double_p), n, _np_as(ps,c_double_p), _np_as(Qs,c_double_p), order, _np_as(cs,c_double_p) )
    return cs, p0

#void sampleMultipole( int n, double* ps_, double* fe_, double* p0_, int order, double* cs ){
lib.sampleMultipole.argtypes  = [ c_int, c_double_p, c_double_p, c_double_p, c_int, c_double_p ]
lib.sampleMultipole.restype   =  None
def sampleMultipole( ps, p0, cs, order=2, fe=None ):
    n = len(ps)
    ps = np.array( ps )
    #p0 = np.array( p0 )
    if fe is None: fe = np.zeros( n, dtype=np.float64 )
    lib.sampleMultipole( n, _np_as(fe,c_double_p), _np_as(ps,c_double_p), _np_as(p0,c_double_p), order, _np_as(cs,c_double_p) )

# void evalGridFFAtPoints( int n, double* ps, double* FFout, double* PLQH, int* nPBC ){
lib.evalGridFFAtPoints.argtypes  = [ c_int, c_double_p, c_double_p, c_double_p, c_bool, c_int_p  ]
lib.evalGridFFAtPoints.restype   =  None
def evalGridFFAtPoints( ps, FFout=None, PLQH=[0.0,0.0,1.0,0.0], bSplit=True, nPBC=None ):
    n = len(ps)
    if FFout is None: FFout=np.zeros( (n,4) )
    PLQH = np.array( PLQH )
    if nPBC is not None: nPBC = np.array( nPBC, dtype=np.int32 )
    lib.evalGridFFAtPoints( n, _np_as(ps,c_double_p), _np_as(FFout,c_double_p), _np_as(PLQH,c_double_p), bSplit, _np_as(nPBC,c_int_p) )
    return FFout

# void sample_func( int n, double* xs, double* ys, int kind ){
lib.sample_func.argtypes  = [c_int, c_double_p, c_double_p, c_int, c_double_p]
lib.sample_func.restype   =  None
def sample_func( xs, ys=None, kind=0, params=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0] ):
    n = len(xs)
    params = np.array( params )
    if ys is None: ys=np.zeros(n)
    lib.sample_func( n, _np_as(xs,c_double_p), _np_as(ys,c_double_p), kind, _np_as(params,c_double_p) )
    return ys

# void sample_func( int n, double* xs, double* ys, int kind ){
lib.sample_func.argtypes  = [c_int, c_double_p, c_double_p, c_int, c_double_p]
lib.sample_func.restype   =  None
def sample_funcEF( xs, EFs=None, kind=0, params=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0] ):
    n = len(xs)
    params = np.array( params )
    if EFs is None: EFs=np.zeros((n,2))
    lib.sample_funcEF( n, _np_as(xs,c_double_p), _np_as(EFs,c_double_p), kind, _np_as(params,c_double_p) )
    return EFs

# void sample_Bspline( double g0, double dg, int ng, double* Gs, int n, double* xs, double* fes , int order, bool bPBC ){
lib.sample_Bspline.argtypes  = [c_double, c_double, c_int, c_double_p, c_int, c_double_p, c_double_p, c_int, c_bool ]
lib.sample_Bspline.restype   =  None
def sample_Bspline( xs, Gs, x0=0.0, dx=1.0, fes=None, order=3, bPBC=True ):
    n = len(xs)
    if fes is None: fes=np.zeros((n,2))
    lib.sample_Bspline(x0, dx, len(Gs), _np_as(Gs,c_double_p), n, _np_as(xs,c_double_p), _np_as(fes,c_double_p), order, bPBC )
    return fes

# void sample_NURBS( double g0, double dg, int ng, double* Gs, double* Ws, int n, double* xs, double* fes ){
lib.sample_NURBS.argtypes  = [c_double, c_double, c_int, c_double_p, c_double_p, c_int, c_double_p, c_double_p ]
lib.sample_NURBS.restype   =  None
def sample_NURBS( xs, Gs, Ws, x0=0.0, dx=1.0, fes=None ):
    n = len(xs)
    if fes is None: fes=np.zeros((n,2))
    lib.sample_NURBS(x0, dx, len(Gs), _np_as(Gs,c_double_p), _np_as(Ws,c_double_p), n, _np_as(xs,c_double_p), _np_as(fes,c_double_p) )
    return fes


# void sample_Bspline2D( double* g0, double* dg, int* ng, double* G, int n, double* ps, double* fes ){
lib.sample_Bspline2D.argtypes  = [c_double_p, c_double_p, c_int_p, c_double_p,                c_int, c_double_p, c_double_p]
lib.sample_Bspline2D.restype   =  None
def sample_Bspline2D( ps, Gs, g0, dg, fes=None):
    n = len(ps)
    g0 = np.array(g0)
    dg = np.array(dg)
    ng = np.array( Gs.shape[::-1], np.int32 )
    if fes is None: fes=np.zeros((n,3))
    lib.sample_Bspline2D( _np_as(g0,c_double_p), _np_as(dg,c_double_p), _np_as(ng,c_int_p), _np_as(Gs,c_double_p), n, _np_as(ps,c_double_p), _np_as(fes,c_double_p) )
    return fes

#void sample_Bspline2D_comb3( double* g0, double* dg, int* ng, double* G, int n, double* ps, double* fes, double* Cs ){
lib.sample_Bspline2D_comb3.argtypes  = [c_double_p, c_double_p, c_int_p, c_double_p,                c_int, c_double_p, c_double_p, c_double_p]
lib.sample_Bspline2D_comb3.restype   =  None
def sample_Bspline2D_comb3( ps, Gs, g0, dg, fes=None, Cs=[1.0,1.0]  ):
    n  = len(ps)
    g0 = np.array(g0)
    dg = np.array(dg)
    ng = np.array( Gs.shape[:2][::-1], np.int32 )
    Cs = np.array( Cs )
    if fes is None: fes=np.zeros((n,3))
    lib.sample_Bspline2D_comb3( _np_as(g0,c_double_p), _np_as(dg,c_double_p), _np_as(ng,c_int_p), _np_as(Gs,c_double_p), n, _np_as(ps,c_double_p), _np_as(fes,c_double_p), _np_as(Cs,c_double_p) )
    return  fes


#void sample_Bspline3D_comb3( double* g0, double* dg, int* ng, double* G, int n, double* ps, double* fes, double* Cs ){
lib.sample_Bspline3D_comb3.argtypes  = [c_double_p, c_double_p, c_int_p, c_double_p,                c_int, c_double_p, c_double_p, c_double_p ]
lib.sample_Bspline3D_comb3.restype   =  None
def sample_Bspline3D_comb3( ps, Gs, g0, dg, fes=None, Cs=[1.0,1.0,1.0] ):
    n  = len(ps)
    g0 = np.array(g0)
    dg = np.array(dg)
    #ng = np.array( Gs.shape[:3][::-1], np.int32 )
    ng = np.array( Gs.shape[:3], np.int32 )
    Cs = np.array( Cs )
    if fes is None: fes=np.zeros((n,4))
    lib.sample_Bspline3D_comb3( _np_as(g0,c_double_p), _np_as(dg,c_double_p), _np_as(ng,c_int_p), _np_as(Gs,c_double_p), n, _np_as(ps,c_double_p), _np_as(fes,c_double_p), _np_as(Cs,c_double_p) )
    return fes


#void sample_Bspline3D            ( double* g0, double* dg, int* ng, double* G,               int n, double* ps, double* fes ){
#void sample_SplineHermite3D_deriv( double* g0, double* dg, int* ng, double* Eg, double* dEg, int n, double* ps, double* fes ){
lib.sample_Bspline3D.argtypes  = [c_double_p, c_double_p, c_int_p, c_double_p,                c_int, c_double_p, c_double_p]
lib.sample_Bspline3D.restype   =  None
def sample_Bspline3D( ps, Eg, g0, dg, fes=None):
    n = len(ps)
    g0 = np.array(g0)
    dg = np.array(dg)
    ng = np.array( Eg.shape[::-1], np.int32 )
    if fes is None: fes=np.zeros((n,4))
    lib.sample_Bspline3D( _np_as(g0,c_double_p), _np_as(dg,c_double_p), _np_as(ng,c_int_p), _np_as(Eg,c_double_p), n, _np_as(ps,c_double_p), _np_as(fes,c_double_p) )
    return fes

#void sample_SplineHermite( double g0, double dg, int ng, double* Eg, int n, double* xs, double* fes ){
lib.sample_SplineHermite.argtypes  = [c_double, c_double, c_int, c_double_p, c_int, c_double_p, c_double_p ]
lib.sample_SplineHermite.restype   =  None
def sample_SplineHermite( xs, Eps, g0=0.0, dg=1.0, fes=None ):
    n = len(xs)
    if fes is None: fes=np.zeros((n,2))
    lib.sample_SplineHermite(g0, dg, len(Eps), _np_as(Eps,c_double_p), n, _np_as(xs,c_double_p), _np_as(fes,c_double_p) )
    return fes

#void sample_SplineHermite_comb( double g0, double dg, int ng, double* Eg, int n, double* xs, double* fes, int ncomb, double* Cs ){
lib.sample_SplineHermite_comb.argtypes  = [c_double, c_double, c_int, c_double_p, c_int, c_double_p, c_double_p, c_int, c_double_p ]
lib.sample_SplineHermite_comb.restype   =  None
def sample_SplineHermite_comb( xs, Eps, Cs, ncomb=2, g0=0.0, dg=1.0, fes=None ):
    n = len(xs)
    if fes is None: fes=np.zeros((n,2))
    Cs = np.array(Cs,dtype=np.float64)
    lib.sample_SplineHermite_comb(g0, dg, len(Eps), _np_as(Eps,c_double_p), n, _np_as(xs,c_double_p), _np_as(fes,c_double_p), ncomb, _np_as(Cs,c_double_p) )
    return fes

#void sample1D_deriv( const double g0, const double dg, const int ng, const Vec2d* FE, const int n, const double* ps, Vec2d* fes ){
lib.sample_SplineHermite1D_deriv.argtypes  = [c_double, c_double, c_int, c_double_p, c_int, c_double_p, c_double_p ]
lib.sample_SplineHermite1D_deriv.restype   =  None
def sample_SplineHermite1D_deriv( ps, EFg, g0, dg, fes=None):
    n = len(ps)
    ng= len(EFg)
    if fes is None: fes=np.zeros((n,2))
    lib.sample_SplineHermite1D_deriv( g0, dg, ng, _np_as(EFg,c_double_p), n, _np_as(ps,c_double_p), _np_as(fes,c_double_p) )
    return fes

#void sample_SplineHermite2D_deriv( double* g0, double* dg, int* ng, double* Eg, double* dEg, int n, double* ps, double* fes ){
lib.sample_SplineHermite2D_deriv.argtypes  = [c_double_p, c_double_p, c_int_p, c_double_p,  c_double_p, c_int, c_double_p, c_double_p]
lib.sample_SplineHermite2D_deriv.restype   =  None
def sample_SplineHermite2D_deriv( ps, Eg, dEg, g0, dg, fes=None):
    n = len(ps)
    g0 = np.array(g0)
    dg = np.array(dg)
    ng = np.array( Eg.shape, np.int32 )
    if fes is None: fes=np.zeros((n,4))
    lib.sample_SplineHermite2D_deriv( _np_as(g0,c_double_p), _np_as(dg,c_double_p), _np_as(ng,c_int_p), _np_as(Eg,c_double_p), _np_as(dEg,c_double_p),  n, _np_as(ps,c_double_p), _np_as(fes,c_double_p) )
    return fes


# void sample_SplineHermite2D_comb( double* g0, double* dg, int* ng, double* Gs, int n, double* ps, double* fes, int ncomb, double* Cs ){
lib.sample_SplineHermite2D_comb.argtypes  = [c_double_p, c_double_p, c_int_p, c_double_p,  c_int, c_double_p, c_double_p, c_int, c_double_p ]
lib.sample_SplineHermite2D_comb.restype   =  None
def sample_SplineHermite2D_comb( ps, Gs, g0, dg, ncomb=2, fes=None, Cs=[1.0,1.0] ):
    n = len(ps)
    g0 = np.array(g0)
    dg = np.array(dg)
    #ng = np.array( Gs.shape[:2][::-1], np.int32 )
    ng = np.array( Gs.shape[:2], np.int32 )    # y-axis should be fastest
    C = np.array(Cs,dtype=np.float64)
    if fes is None: fes=np.zeros((n,3))
    lib.sample_SplineHermite2D_comb( _np_as(g0,c_double_p), _np_as(dg,c_double_p), _np_as(ng,c_int_p), _np_as(Gs,c_double_p), n,_np_as(ps,c_double_p), _np_as(fes,c_double_p), ncomb, _np_as(C,c_double_p) )
    return fes

#void sample_SplineHermite3D_deriv( double* g0, double* dg, int* ng, double* Eg, double* dEg, int n, double* ps, double* fes ){
lib.sample_SplineHermite3D_deriv.argtypes  = [c_double_p, c_double_p, c_int_p, c_double_p, c_double_p,  c_int, c_double_p, c_double_p]
lib.sample_SplineHermite3D_deriv.restype   =  None
def sample_SplineHermite3D_deriv( ps, Eg, dEg, g0, dg, fes=None):
    n = len(ps)
    g0 = np.array(g0)
    dg = np.array(dg)
    ng = np.array( Eg.shape, np.int32 )
    if fes is None: fes=np.zeros((n,4))
    lib.sample_SplineHermite3D_deriv( _np_as(g0,c_double_p), _np_as(dg,c_double_p), _np_as(ng,c_int_p), _np_as(Eg,c_double_p), _np_as(dEg,c_double_p), n, _np_as(ps,c_double_p), _np_as(fes,c_double_p) )
    return fes

#void sample_SplineHermite3D_comb3( double* g0, double* dg, int* ng, double* EFg, int n, double* ps, double* fes, double* Cs ){
lib.sample_SplineHermite3D_comb3.argtypes  = [c_double_p, c_double_p, c_int_p, c_double_p, c_int, c_double_p, c_double_p, c_double_p ]
lib.sample_SplineHermite3D_comb3.restype   =  None
def sample_SplineHermite3D_comb3( ps, EFg, g0, dg, fes=None, Cs=[1.0,1.0,1.0] ):
    n = len(ps)
    g0 = np.array(g0)
    dg = np.array(dg)
    ng = np.array( EFg.shape, np.int32 )
    C = np.array(Cs,dtype=np.float64)
    if fes is None: fes=np.zeros((n,4))
    lib.sample_SplineHermite3D_comb3( _np_as(g0,c_double_p), _np_as(dg,c_double_p), _np_as(ng,c_int_p), _np_as(EFg,c_double_p), n, _np_as(ps,c_double_p), _np_as(fes,c_double_p), _np_as(C,c_double_p) )
    return fes

#void sample_SplineHermite2D( double* g0, double* dg, int* ng, double* Eg, int n, double* ps, double* fes )
lib.sample_SplineHermite2D.argtypes  = [c_double_p, c_double_p, c_int_p, c_double_p, c_int, c_double_p, c_double_p]
lib.sample_SplineHermite2D.restype   =  None
def sample_SplineHermite2D( ps, Eg, g0, dg, fes=None):
    n = len(ps)
    g0 = np.array(g0)
    dg = np.array(dg)
    ng = np.array( Eg.shape[::-1], np.int32 )
    if fes is None: fes=np.zeros((n,3))
    lib.sample_SplineHermite2D( _np_as(g0,c_double_p), _np_as(dg,c_double_p), _np_as(ng,c_int_p), _np_as(Eg,c_double_p), n, _np_as(ps,c_double_p), _np_as(fes,c_double_p) )
    return fes

#void sample_SplineHermite3D( double* g0, double* dg, int* ng, double* Eg, int n, double* ps, double* fes )
lib.sample_SplineHermite3D.argtypes  = [c_double_p, c_double_p, c_int_p, c_double_p, c_int, c_double_p, c_double_p]
lib.sample_SplineHermite3D.restype   =  None
def sample_SplineHermite3D( ps, Eg, g0, dg, fes=None):
    n = len(ps)
    g0 = np.array(g0)
    dg = np.array(dg)
    if fes is None: fes=np.zeros((n,4))
    ng = np.array( Eg.shape, np.int32 )
    lib.sample_SplineHermite3D( _np_as(g0,c_double_p), _np_as(dg,c_double_p), _np_as(ng,c_int_p), _np_as(Eg,c_double_p), n, _np_as(ps,c_double_p), _np_as(fes,c_double_p) )
    return fes

#void sample_SplineHermite3D_f( float* g0, float* dg, int* ng, float* Eg, int n, float* ps, float* fes ){
lib.sample_SplineHermite3D_f.argtypes  = [c_float_p, c_float_p, c_int_p, c_float_p, c_int, c_float_p, c_float_p]
lib.sample_SplineHermite3D_f.restype   =  None
def sample_SplineHermite3D_f( ps, Eg, g0, dg, fes=None):
    n = len(ps)
    g0 = np.array(g0, dtype=np.float32)
    dg = np.array(dg, dtype=np.float32)
    ng = np.array( Eg.shape, np.int32 )
    if fes is None: fes=np.zeros((n,4), dtype=np.float32)
    lib.sample_SplineHermite3D_f( _np_as(g0,c_float_p), _np_as(dg,c_float_p), _np_as(ng,c_int_p), _np_as(Eg,c_float_p), n, _np_as(ps,c_float_p), _np_as(fes,c_float_p) )
    return fes

#void sampleCoulombPBC( int nps, double* ps, double* fe, int natom,  double* apos, double* Qs, double* lvec, int* nPBC, double Rdamp ){
lib.sampleCoulombPBC.argtypes  = [c_int, c_double_p, c_double_p, c_int, c_double_p, c_double_p, c_double_p, c_int_p, c_double]
lib.sampleCoulombPBC.restype   =  None
def sampleCoulombPBC( ps, apos, Qs, lvec=[[10.0,0.0,0.0],[10.0,0.0,0.0],[10.0,0.0,0.0]], nPBC=[1,1,1], fe=None, Rdamp=1e-32 ):
    nps   = len(ps)
    natom = len(apos)
    ps    = np.array(ps,   dtype=np.float64)
    apos  = np.array(apos, dtype=np.float64)
    Qs    = np.array(Qs,   dtype=np.float64)
    if fe is None: fe=np.zeros( (nps,4) )
    lvec  = np.array(lvec, dtype=np.float64)
    nPBC = np.array(nPBC,  dtype=np.int32)
    lib.sampleCoulombPBC( nps, _np_as(ps,c_double_p), _np_as(fe,c_double_p), natom, _np_as(apos,c_double_p), _np_as(Qs,c_double_p), _np_as(lvec,c_double_p), _np_as(nPBC,c_int_p), Rdamp )
    return fe


# void sample_SplineConstr( double x0, double dx, int np, double* Eps, int n, double* xs, double* Es, double* Fs ){
lib.sample_SplineConstr.argtypes  = [c_double, c_double, c_int, c_double_p, c_int, c_double_p, c_double_p, c_double_p]
lib.sample_SplineConstr.restype   =  None
def sample_SplineConstr( xs, Eps, x0=0.0, dx=1.0, Es=None, Fs=None):
    n = len(xs)
    if Es is None: Es=np.zeros(n)
    if Fs is None: Fs=np.zeros(n)
    lib.sample_SplineConstr(x0, dx, len(Eps), _np_as(Eps,c_double_p), n, _np_as(xs,c_double_p), _np_as(Es,c_double_p), _np_as(Fs,c_double_p))
    return Es,Fs

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
def sampleNonBond( rs, Es=None, Fs=None, kind=1, REQi=(1.487,0.0006808,0.0), REQj=(1.487,0.0006808,0.0), K=-1.0, Rdamp=1.0 ):
    n =len(rs)
    if Es is None: Es=np.zeros(n)
    if Fs is None: Fs=np.zeros(n)
    rs=np.array(rs)
    REQi=np.array(REQi)
    REQj=np.array(REQj)
    lib.sampleNonBond(n, rs, Es, Fs, kind, REQi, REQj, K, Rdamp)
    return Es,Fs

#  int findHbonds( double Rcut, double Hcut, double angMax ){
lib.findHbonds.argtypes  = [ c_double, c_double, c_double]
lib.findHbonds.restype   =  c_int
def findHbonds( Rcut=4.0, Hcut=0.0001, angMax=30.0 ):
    return lib.findHbonds( Rcut, Hcut, angMax )

#  int sampleHbond( int ib, int n, double* rs, double* Es, double* fs, int kind, Vec2d mask, double K, double Rdamp ){
lib.sampleHbond.argtypes  = [c_int,c_int, array1d, array1d, array1d, c_int, c_double, c_double, c_double, c_double, c_double, c_char_p ]
lib.sampleHbond.restype   =  c_int
def sampleHbond( ib, rs, Es=None, Fs=None, kind=1, maskQ=1.0, maskH=1.0, K=-1.0, Rdamp=1.0, dcomp=1.0 ):
    n =len(rs)
    if Es is None: Es=np.zeros(n)
    if Fs is None: Fs=np.zeros(n)
    rs  =np.array(rs)
    s = ctypes.create_string_buffer(1024)
    lib.sampleHbond(ib, n, rs, Es, Fs, kind, maskQ, maskH, K, Rdamp, dcomp, s )
    s = s.value.decode('utf-8')
    return Es,Fs,s

#void sampleNonBondTypes( int n, double* rs, double* Es, double* fs, int kind, double qH, double qX, double K, double Rdamp, double dcomp, char* type_str ){
lib.sampleNonBondTypes.argtypes  = [ c_int, array1d, array1d, array1d, c_int, c_double, c_double, c_double, c_double, c_double, c_char_p ]
lib.sampleNonBondTypes.restype   =  c_int
def sampleNonBondTypes( type_str, rs, Es=None, Fs=None, kind=1, qH=1.0, qX=1.0, K=-1.0, Rdamp=1.0, dcomp=1.0 ):
    n =len(rs)
    if Es is None: Es=np.zeros(n)
    if Fs is None: Fs=np.zeros(n)
    rs  =np.array(rs)
    s = type_str.encode('utf8')
    lib.sampleNonBondTypes( n, rs, Es, Fs, kind, qH, qX, K, Rdamp, dcomp, s )
    return Es,Fs

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

# void sampleSurf_vecs(int n, double* poss_, double* FEs_, int kind, int ityp, double RvdW, double EvdW, double Q, double K, double RQ, int npbc, bool bSave){
lib.sampleSurf_vecs.argtypes  = [ c_int, array2d, array2d, c_int, c_int,  c_double, c_double, c_double, c_double, c_double, c_int, c_bool]
lib.sampleSurf_vecs.restype   =  None
def sampleSurf_vecs(ps, FEs=None, kind=1, ityp=-1, RvdW=1.487, EvdW=0.0006808, Q=0.0, K=+1.5, Rdamp=0.1, npbc=1, bSave=False ):
    n =len(ps)
    if FEs is None: FEs=np.zeros((n,4))
    lib.sampleSurf_vecs( n, ps, FEs, kind, ityp, RvdW, EvdW, Q, K, Rdamp, npbc, bSave )
    return FEs


# void sampleSurf_new( int n, double* ps_, double* FEout_, int kind, double* REQ_, double K, double RQ ){
lib.sampleSurf_new.argtypes  = [ c_int, array2d, array2d, c_int, array1d, c_double, c_double ]
lib.sampleSurf_new.restype   =  None
def sampleSurf_new( ps, PLQH, FEout=None, mode=1, K=-1.0, Rdamp=1.0 ):
    n =len(ps)
    if FEout is None: FEout=np.zeros( (n,4) )
    PLQH=np.array(PLQH)
    lib.sampleSurf_new( n, ps, FEout, mode, PLQH, K, Rdamp )
    return FEout

# # void sampleSurf_vecs(char* name, int n, double* rs, double* Es, double* fs, int kind, double*REQ_, double K, double Rdamp ){
# lib.sampleSurf_vecs.argtypes  = [c_char_p, c_int, array2d, array1d, array2d, c_int, c_int, c_double, c_double, c_double, array1d, c_bool]
# lib.sampleSurf_vecs.restype   =  None
# def sampleSurf_vecs( name, poss, Es=None, fs=None, kind=1, atyp=0, Q=0.0, K=-1.0, Rdamp=1.0, pos0=(0.,0.,0.), bSave=False ):
#     if name is not None: name=name.encode('utf8')
#     n =len(poss)
#     if Es is None: Es=np.zeros(n)
#     if fs is None: fs=np.zeros((n,3))
#     pos0=np.array(pos0)
#     lib.sampleSurf_vecs( name, n, poss, Es, fs, kind, atyp, Q, K, Rdamp, pos0, bSave )
#     return Es,fs

# void setupCollisionDamping( int ndampstep, double damping_medium, double collisionDamping, double collisionDamping_NB, double col_damp_dRcut ){
lib.setupCollisionDamping.argtypes  = [c_int, c_double, c_double, c_double, c_double, c_double, c_double]
lib.setupCollisionDamping.restype   =  None
def setupCollisionDamping( nstep=10, medium=0.02, cB=-1.0, cA=-1.0, cNB=-1.0, dRcut1=-0.2, dRcut2=0.3 ):
    #print( "setupCollisionDamping(): ",nstep,medium,cB,cA,cNB,dRcut1,dRcut2 )
    lib.setupCollisionDamping( nstep, medium, cB, cA, cNB, dRcut1, dRcut2 )

#void setup_accel(int nstep_acc_min_, double cos_vf_acc_ ){
lib.setup_accel.argtypes  = [c_int, c_double]
lib.setup_accel.restype   =  None
def setup_accel( nstep_acc_min=10, cos_vf_acc=0.5 ):
    #print( "setup_accel(): ",nstep_acc_min,cos_vf_acc )
    lib.setup_accel( nstep_acc_min, cos_vf_acc )

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
    print("getBuff(): %s ptr=%s" %(name,ptr))
    if ptr is None: print("getBuff(): %s not found" %name); return None
    return np.ctypeslib.as_array( ptr, shape=sh)

#double* getBuff(const char* name){
lib.getBBuff.argtypes = [c_char_p]
lib.getBBuff.restype  = c_bool_p
def getBBuff(name,sh):
    if not isinstance(sh, tuple): sh=(sh,)
    name=name.encode('utf8')
    ptr = lib.getBBuff(name)
    return np.ctypeslib.as_array( ptr, shape=sh)


# Add the new function definition for get_molecule_natoms
lib.get_molecule_natoms.argtypes = []
lib.get_molecule_natoms.restype = c_int

def get_molecule_natoms():
    """
    Retrieves the number of atoms in the molecule.

    Returns:
        int: The number of atoms in the molecule.
    """
    return lib.get_molecule_natoms()

# Add the new function definitions
lib.get_gridFF_info.argtypes = [c_int_p, c_double_p]
lib.get_gridFF_info.restype = None

def get_gridFF_info():
    """
    Retrieves grid FF information.

    Returns:
        tuple: A tuple containing shift0, pos0, natoms, and atoms_.size.
    """
    global gff_shift0, gff_pos0, gff_cell, gff_dCell, gff_natoms, gff_natoms_
    int_data = np.zeros(2, dtype=np.int32)
    float_data = np.zeros(24, dtype=np.float64)
    lib.get_gridFF_info(int_data.ctypes.data_as(c_int_p),
                        float_data.ctypes.data_as(c_double_p))
    gff_shift0 = float_data[0:3]
    gff_pos0 = float_data[3:6]
    gff_cell = float_data[6:15].reshape((3, 3))      # Reshape to 3x3 matrix
    gff_dCell = float_data[15:24].reshape((3, 3))     # Reshape to 3x3 matrix
    gff_natoms = int_data[0]
    gff_natoms_ = int_data[1]
    print("GridFF info -> shift0:", gff_shift0, " pos0:", gff_pos0,
          " substrate natoms:", gff_natoms, " atoms_.size:", gff_natoms_)
    return gff_shift0, gff_pos0,gff_cell,gff_dCell, gff_natoms, gff_natoms_


# Add the new function definitions
lib.get_atom_positions.argtypes = [c_double_p, c_double_p]
lib.get_atom_positions.restype = None


def get_atom_positions():
    """
    Retrieves the atom positions of the substrate and molecule.

    Returns:
        tuple: A tuple containing two NumPy arrays:
            - substrate_apos (numpy.ndarray): Substrate atom positions (natoms_substrate, 3).
            - molecule_apos (numpy.ndarray): Molecule atom positions (natoms_molecule, 3).
    """
    global gff_natoms
    natoms_molecule = get_molecule_natoms()
    substrate_apos = np.zeros((gff_natoms, 3), dtype=np.float64)  # Use gff_natoms
    molecule_apos = np.zeros((natoms_molecule, 3), dtype=np.float64)  # Use natoms_molecule
    lib.get_atom_positions(substrate_apos.ctypes.data_as(c_double_p),
                            molecule_apos.ctypes.data_as(c_double_p))
    return substrate_apos, molecule_apos


#def getBuffs( nnode, npi, ncap, nbond, NEIGH_MAX=4 ):
def getBuffs( NEIGH_MAX=4 ):
    print("getBuffs()")
    init_buffers( False )
    #natom=nnode+ncap
    #nvecs=natom+npi
    #nDOFs=nvecs*3
    global ffflags
    ffflags = getBBuff( "ffflags" , (14,) )
    global ndims,Es
    ndims = getIBuff( "ndims", (9,) )  # [nDOFs,natoms,nnode,ncap,npi,nbonds]
    global nDOFs,natoms,nnode,ncap,npi,nvecs,nbonds,ne,ie0
    # MFF_lib.cpp::init_buffers() ndims{nDOFs=9,natoms=3,nnode=1,ncap=2,npi=0,nbonds=2,nvecs=3,ne=0,ie0=3}
    nDOFs=ndims[0]; natoms=ndims[1];  nnode=ndims[2]; ncap=ndims[3]; npi=ndims[4]; nbonds=ndims[5]; nvecs=ndims[6]; ne=ndims[7]; ie0=ndims[8]
    print( "getBuffs(): nDOFs=%i nvecs=%i  natoms=%i nnode=%i ncap=%i npi=%i nbonds=%i nvecs=%i ne=%i ie0=%i " %(nDOFs,nvecs,natoms,nnode,ncap,npi,nbonds,nvecs,ne,ie0) )
    Es    = getBuff ( "Es",    (6,) )  # [ Etot,Eb,Ea, Eps,EppT,EppI; ]
    global DOFs,fDOFs,vDOFs,apos,fapos,REQs,PLQs,pipos,fpipos,bond_l0,bond_k, bond2atom,neighs,selection
    #Ebuf     = getEnergyTerms( )
    apos      = getBuff ( "apos",     (natoms,3) )
    fapos     = getBuff ( "fapos",    (natoms,3) )
    REQs      = getBuff ( "REQs",     (natoms,4) )
    PLQs      = getBuff ( "PLQs",     (natoms,4) )
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
    neighs    = getIBuff( "neighs",  (nnode,NEIGH_MAX) )
    selection = getIBuff( "selection",  (natoms) )
    print( "getBuffs DONE" )

def getBuffs_UFF( NEIGH_MAX=4 ):
    print("getBuffs_UFF()")
    init_buffers( True )
    # ndims :  int _natoms, nbonds, nangles, ndihedrals, ninversions, nf; // 5
    # ndims :  int i0dih,i0inv,i0ang,i0bon;                               // 4
    # Es: double Etot, Eb, Ea, Ed, Ei;                                    // 5
    global ffflags
    ffflags = getBBuff( "ffflags" , (14,) )
    global ndims,Es
    ndims = getIBuff( "ndims", (10,) )  # [_natoms, nbonds, nangles, ndihedrals, ninversions, nf, i0dih,i0inv,i0ang,i0bon]
    global natoms, nbonds, nangles, ndihedrals, ninversions, nf, i0dih,i0inv,i0ang,i0bon
    natoms=ndims[0]; nbonds=ndims[1]; nangles=ndims[2]; ndihedrals=ndims[3]; ninversions=ndims[4]; nf=ndims[5]; i0dih=ndims[6]; i0inv=ndims[7]; i0ang=ndims[8]; i0bon=ndims[9]
    print( "getBuffs(): natoms=%i nbonds=%i nangles=%i ndihedrals=%i ninversions=%i nf=%i i0dih=%i i0inv=%i i0ang=%i i0bon=%i " %(natoms,nbonds,nangles,ndihedrals,ninversions,nf,i0dih,i0inv,i0ang,i0bon) )
    Es    = getBuff ( "Es",    (5,) )  # [ Etot,Eb,Ea,Ed,Ei ]
    global apos,fapos,REQs,hneigh,fint,bonAtoms,angAtoms,dihAtoms,invAtoms,neighs,neighBs,bonParams,angParams,dihParams,invParams,angNgs,dihNgs,invNgs
    #Ebuf     = getEnergyTerms( )
    apos      = getBuff ( "apos",     (natoms,3) )
    fapos     = getBuff ( "fapos",    (natoms,3) )
    REQs      = getBuff ( "REQs",     (natoms,4) )
    # ------ UFF
    hneigh    = getBuff ( "hneigh",    (natoms*NEIGH_MAX,4) )
    fint      = getBuff ( "fint",      (nf,3) )

    bonParams = getBuff ( "bonParams", (nbonds,2)      )
    angParams = getBuff ( "angParams", (nangles,5)     )
    dihParams = getBuff ( "dihParams", (ndihedrals,3)  )
    invParams = getBuff ( "invParams", (ninversions,4) )

    bonAtoms  = getIBuff( "bonAtoms",  (nbonds,2)      )
    angAtoms  = getIBuff( "angAtoms",  (nangles,3)     )
    dihAtoms  = getIBuff( "dihAtoms",  (ndihedrals,4)  )
    invAtoms  = getIBuff( "invAtoms",  (ninversions,4) )

    angNgs    = getIBuff( "angNgs",    (nangles,2)     )
    dihNgs    = getIBuff( "dihNgs",    (ndihedrals,3)  )
    invNgs    = getIBuff( "invNgs",    (ninversions,3) )

    neighs    = getIBuff( "neighs",    (natoms,4)      )
    neighBs   = getIBuff( "neighBs",   (natoms,4)      )

    # // Quat4i *  neighBs   __attribute__((aligned(64))) = 0; // [natoms]      bond indices for each neighbor
    # // Vec2i  *  bonAtoms  __attribute__((aligned(64))) = 0; // [nbonds]      bonds atoms
    # // Vec2d  *  bonParams __attribute__((aligned(64))) = 0; // [nbonds]      bonds parameters
    # // Vec3i  *  angAtoms  __attribute__((aligned(64))) = 0; // [nangles]     angles atoms
    # // double5*  angParams __attribute__((aligned(64))) = 0; // [nangles]     angles parameters
    # // Quat4i *  dihAtoms  __attribute__((aligned(64))) = 0; // [ndihedrals]  dihedrals atoms
    # // Vec3d  *  dihParams __attribute__((aligned(64))) = 0; // [ndihedrals]  dihedrals parameters
    # // Quat4i *  invAtoms  __attribute__((aligned(64))) = 0; // [ninversions] inversions atoms
    # // Quat4d *  invParams __attribute__((aligned(64))) = 0; // [ninversions] inversions parameters

    # // Vec2i * angNgs __attribute__((aligned(64))) = 0; // [nangles]     angles neighbor index
    # // Vec3i * dihNgs __attribute__((aligned(64))) = 0; // [ndihedrals]  dihedrals neighbor index
    # // Vec3i * invNgs __attribute__((aligned(64))) = 0; // [ninversions] inversions neighbor index

    

    print( "getBuffs DONE" )

#  void init_buffers()
lib.init_buffers.argtypes  = []
lib.init_buffers.restype   =  None
def init_buffers( bUFF=False ):
    global bBuffersInitialized
    if bBuffersInitialized: 
        print("init_buffers():Buffers already initialized => return")
        return
    bBuffersInitialized = True
    #print( "init_buffers()" )
    if bUFF:
        lib.init_buffers_UFF()
    else:
        lib.init_buffers()


'''
#  void init_params(const char* fatomtypes, const char* fbondtypes)
lib.init_params.argtypes  = [c_char_p, c_char_p, c_char_p]
lib.init_params.restype   =  None
def init_params(fatomtypes, fbondtypes, fbondangles ):
    fatomtypes = fatomtypes.encode('utf8')
    fbondtypes = fbondtypes.encode('utf8')
    fbondangles = fbondangles.encode('utf8')
    return lib.init_params(fatomtypes,fbondtypes,fbondangles)
'''

'''
#  void init_nonbond()
lib.init_nonbond.argtypes  = []
lib.init_nonbond.restype   =  None
def init_nonbond():
    return lib.init_nonbond()
'''

# void print_setup(){
lib.print_setup.argtypes  = []
lib.print_setup.restype   =  None
def print_setup():
    return lib.print_setup()


#  void print_debugs( bool bParams, bool bNeighs, bool bShifts ){
lib.print_debugs.argtypes  = [c_bool, c_bool, c_bool]
lib.print_debugs.restype   =  None
def print_debugs(bParams=True, bNeighs=True, bShifts=False):
    return lib.print_debugs(bParams, bNeighs, bShifts)

#  void clear         (                      )
lib.clear.argtypes  = []
lib.clear.restype   =  None
def clear():
    return lib.clear()

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

# interface to init on the C++ side (MMFF_lib.cpp)

#void* init( char* xyz_name, char* surf_name, char* smile_name, bool bMMFF, bool bEpairs, bool bUFF, bool b141, bool bSimple, bool bConj, bool bCumulene, int* nPBC, double gridStep, char* sElementTypes, char* sAtomTypes, char* sBondTypes, char* sAngleTypes, char* sDihedralTypes ){
#lib.init(      cstr(xyz_name), cstr(surf_name), cstr(smile_name),  bMMFF, bEpairs,   bUFF,   b141, bSimple,  bConj, bCumulene,    nPBC, gridStep, cstr(sElementTypes), cstr(sAtomTypes), cstr(sBondTypes), cstr(sAngleTypes), cstr(sDihedralTypes) )
lib.init.argtypes  = [c_char_p,        c_char_p,         c_char_p, c_bool,  c_bool, c_bool, c_bool,  c_bool, c_bool,    c_bool, array1i, c_double,            c_char_p,         c_char_p,         c_char_p,          c_char_p,             c_char_p]
lib.init.restype   =  c_void_p
def init(
        xyz_name  =None,
        surf_name =None,
        smile_name=None,
        sElementTypes  = "data/ElementTypes.dat",
        sAtomTypes     = "data/AtomTypes.dat",
        sBondTypes     = "data/BondTypes.dat",
        sAngleTypes    = "data/AngleTypes.dat",
        sDihedralTypes = "data/DihedralTypes.dat",
        bMMFF=True,
        bEpairs=False,
        nPBC=(1,3,0),
        gridStep=0.1,
        bUFF=False,
        b141=True,
        bSimple=False,
        bConj=True,
        bCumulene=True
    ):
    global glob_bMMFF, glob_bUFF
    glob_bMMFF = bMMFF
    glob_bUFF = bUFF
    nPBC=np.array(nPBC,dtype=np.int32)
    return lib.init( cstr(xyz_name), cstr(surf_name), cstr(smile_name), bMMFF, bEpairs, bUFF, b141, bSimple, bConj, bCumulene, nPBC, gridStep, cstr(sElementTypes), cstr(sAtomTypes), cstr(sBondTypes), cstr(sAngleTypes), cstr(sDihedralTypes) )


#void initParams          ( const char* sElementTypes, const char* sAtomTypes, const char* sBondTypes, const char* sAngleTypes, const char* sDihedralTypes ){
lib.initParams.argtypes  = [c_char_p, c_char_p, c_char_p, c_char_p, c_char_p]
lib.initParams.restype   =  None
def initParams(
        sElementTypes  = "data/ElementTypes.dat",
        sAtomTypes     = "data/AtomTypes.dat",
        sBondTypes     = "data/BondTypes.dat",
        sAngleTypes    = "data/AngleTypes.dat",
        sDihedralTypes = "data/DihedralTypes.dat",
    ):
    return lib.initParams(  cstr(sElementTypes), cstr(sAtomTypes), cstr(sBondTypes), cstr(sAngleTypes), cstr(sDihedralTypes) )

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

'''
#  void buildFF( bool bNonBonded_, bool bOptimizer_ )
lib.buildFF.argtypes  = [c_bool, c_bool]
lib.buildFF.restype   =  None
def buildFF(bNonBonded=True, bOptimizer=True):
    return lib.buildFF(bNonBonded, bOptimizer)
'''

# #  int loadmol( const char* fname_mol )
# lib.loadmol.argtypes  = [c_char_p]
# lib.loadmol.restype   =  c_int
# def loadmol(fname_mol):
#     fname_mol=fname_mol.encode('utf8')
#     return lib.loadmol(fname_mol)

# #  void initWithMolFile(char* fname_mol, bool bNonBonded_, bool bOptimizer_ )
# lib.initWithMolFile.argtypes  = [c_char_p, c_bool, c_bool]
# lib.initWithMolFile.restype   =  None
# def initWithMolFile(fname_mol, bNonBonded=True, bOptimizer=True):
#     fname_mol = fname_mol.encode('utf8')
#     return lib.initWithMolFile(fname_mol, bNonBonded, bOptimizer)

# #  void initWithMolFile(char* fname_mol, bool bNonBonded_, bool bOptimizer_ )
# lib.initWithSMILES.argtypes  = [c_char_p, c_bool, c_bool, c_bool, c_bool]
# lib.initWithSMILES.restype   =  None
# def initWithSMILES(fname_mol, bPrint=True, bCap=True, bNonBonded=True, bOptimizer=True):
#     fname_mol = fname_mol.encode('utf8')
#     return lib.initWithSMILES(fname_mol, bPrint, bCap, bNonBonded, bOptimizer)

#  void setSwitches( int doAngles, int doPiPiT, int  doPiSigma, int doPiPiI, int doBonded_, int PBC, int CheckInvariants )
lib.setSwitches.argtypes  = [c_int, c_int, c_int , c_int, c_int, c_int, c_int]
lib.setSwitches.restype   =  None
def setSwitches(doAngles=0, doPiPiT=0, doPiSigma=0, doPiPiI=0, doBonded=0, PBC=0, CheckInvariants=0, bSaveToDatabase=0):
    return lib.setSwitches(doAngles, doPiPiT, doPiSigma, doPiPiI, doBonded, PBC, CheckInvariants, bSaveToDatabase)

# void setSwitches2( int CheckInvariants, int PBC, int NonBonded, int NonBondNeighs,  int SurfAtoms, int GridFF, int MMFF, int Angles, int PiSigma, int PiPiI ){
lib.setSwitches2.argtypes  = [c_int, c_int, c_int, c_int, c_int, c_int, c_int, c_int, c_int, c_int]
lib.setSwitches2.restype   =  None
def setSwitches( CheckInvariants=0, PBC=0, NonBonded=0, NonBondNeighs=0, SurfAtoms=0, GridFF=0, MMFF=0, Angles=0, PiSigma=0, PiPiI=0):
    return lib.setSwitches2(CheckInvariants, PBC, NonBonded, NonBondNeighs, SurfAtoms, GridFF, MMFF, Angles, PiSigma, PiPiI)

# void setSwitchesUFF( int DoBond, int DoAngle, int DoDihedral, int DoInversion, int DoAssemble, int SubtractBondNonBond, int ClampNonBonded )
lib.setSwitchesUFF.argtypes  = [c_int, c_int, c_int, c_int, c_int, c_int, c_int]
lib.setSwitchesUFF.restype   =  None
def setSwitchesUFF( DoBond=0, DoAngle=0, DoDihedral=0, DoInversion=0, DoAssemble=0, SubtractBondNonBond=0, ClampNonBonded=0):
    return lib.setSwitchesUFF(DoBond, DoAngle, DoDihedral, DoInversion, DoAssemble, SubtractBondNonBond, ClampNonBonded)

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
def setTrjName(trj_fname_="trj.xyz", savePerNsteps=1, bDel=True, nPBC=(1,1,1) ):
    if bDel: open(trj_fname_,"w").close()
    global trj_fname
    trj_fname=cstr(trj_fname_)
    nPBC=np.array(nPBC, np.int32)
    return lib.setTrjName( trj_fname, savePerNsteps, _np_as(nPBC,c_int_p)  )

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

# #bool relax( int niter, double Ftol )
# lib.relax.argtypes  = [c_int, c_double, c_bool ]
# lib.relax.restype   =  c_bool
# def relax(niter=100, Ftol=1e-6, bWriteTrj=False ):
#     return lib.relax(niter, Ftol, bWriteTrj )

# #  int  run( int nstepMax, double dt, double Fconv=1e-6, int ialg=0 ){
# lib. run.argtypes  = [c_int, c_double, c_double, c_int ]
# lib. run.restype   =  c_int
# def  run(nstepMax=1000, dt=-1, Fconv=1e-6, ialg=2 ):
#     return lib.run(nstepMax, dt, Fconv, ialg )

#  int  run( int nstepMax, double dt, double Fconv=1e-6, int ialg=0 ){
lib. run.argtypes  = [c_int, c_double, c_double, c_int, c_double, c_double_p, c_double_p, c_double_p, c_double_p, c_bool ]
lib. run.restype   =  c_int
def run(nstepMax=1000, dt=-1, Fconv=1e-6, ialg=2, damping=-1.0, outE=None, outF=None, outV=None, outVF=None, omp=False):
    return lib.run(nstepMax, dt, Fconv, ialg, damping, _np_as(outE,c_double_p), _np_as(outF,c_double_p), _np_as(outV,c_double_p), _np_as(outVF,c_double_p), omp )

#lib.scan.argtypes  = [c_int, array2d, array2d, array1d, array2d, array2d, c_bool, c_bool, c_int, c_double, c_double, c_double]
# int  scan( int nconf, double* poss, double* rots, double* Es, double* aforces, double* aposs, bool omp, bool bRelax, int niter_max, double dt, double Fconv, double Flim ){
lib.scan.argtypes = [ c_int, c_double_p, c_double_p, c_double_p, c_double_p, c_double_p, c_bool, c_bool, c_int, c_double, c_double, c_double ]
lib.scan.restype   =  None
def scan(poss, rots=None, Es=None, aforces=None, aposs=None,  bF=False,bP=False, omp=False, bRelax=False, niter_max=10000, dt=0.05, Fconv=1e-5, Flim=100.0 ):
    nconf=len(poss)
    if Es is None: Es=np.zeros(nconf)
    if (aforces is None) and bF: aforces=np.zeros( (nconf,natoms,3) )
    if (aposs is None) and bP:   aposs=np.zeros(   (nconf,natoms,3) )
    #lib.scan(nconf, poss, rots, Es, aforces, aposs, omp, bRelax, niter_max, dt, Fconv, Flim )
    lib.scan( nconf, _np_as(poss,c_double_p), _np_as(rots,c_double_p), _np_as(Es,c_double_p), _np_as(aforces,c_double_p), _np_as(aposs,c_double_p), omp, bRelax, niter_max, dt, Fconv, Flim )
    return Es, aforces, aposs

# int getHessian3x3( int n_atoms, int* inds, double* out_hessians, double dx, bool bDiag )
lib.getHessian3x3.argtypes = [c_int, c_int_p, c_double_p, c_double, c_bool]
lib.getHessian3x3.restype  = None
def getHessian3x3(inds, dx=1e-4, bDiag=True):
    n = len(inds)
    inds_c = np.ascontiguousarray(inds, dtype=np.int32)
    out    = np.zeros((n,4,3), dtype=np.float64)
    lib.getHessian3x3( n, _np_as(inds_c, c_int_p), _np_as(out,    c_double_p), dx, bDiag )
    return out

# int getHessian3Nx3N_selected( int n_atoms, int* inds, double* out_hessian_full, double dx )
lib.getHessian3Nx3N.argtypes = [c_int, c_int_p,  c_double_p, c_double]
lib.getHessian3Nx3N.restype  = None
def getHessian3Nx3N(inds, dx=1e-4):
    n = len(inds)
    inds_c = np.array(inds, dtype=np.int32)
    dim    = 3 * n
    out    = np.zeros((dim,dim), dtype=np.float64)
    lib.getHessian3Nx3N( n, _np_as(inds_c, c_int_p), _np_as(out,    c_double_p), dx )
    return out

# void scan_atoms_rigid(int nscan, int nsel, int* inds, double* scan_pos, double* out_forces, double* out_Es, bool bRelative){
lib.scan_atoms_rigid.argtypes = [ c_int, c_int, c_int_p,  c_double_p, c_double_p, c_double_p, c_bool]
lib.scan_atoms_rigid.restype = None
def scan_atoms_rigid(inds, poss, out_forces=None, out_Es=None, bRelative=True, bForce=True, bEnergy=True ):
    nscan = len(poss)
    nsel  = len(inds)
    if (out_forces is None) and bForce:  out_forces = np.zeros((nscan,nsel,3), dtype=np.float64)
    if (out_Es     is None) and bEnergy: out_Es     = np.zeros( nscan,         dtype=np.float64)
    inds_c = np.ascontiguousarray(inds, dtype=np.int32)
    pos_c  = np.ascontiguousarray(poss, dtype=np.float64)
    lib.scan_atoms_rigid( nscan, nsel, _np_as(inds_c, c_int_p), _np_as(pos_c, c_double_p), _np_as(out_forces, c_double_p), _np_as(out_Es, c_double_p), bRelative )
    return  out_Es, out_forces

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
    print("DEBUG_3")
    return lib.optimizeLattice_1d(_np_as(dlvec,c_double_p), n1, n2, initMode, tol)

# ========= Constrains

#  void addDistConstrain(  int i0,int i1, double lmin,double lmax,double kmin,double kmax,double flim, double k, double* shift ){
lib.addDistConstrain.argtypes  = [c_int, c_int, c_double, c_double, c_double, c_double, c_double, c_double_p ]
lib.addDistConstrain.restype   =  None
def addDistConstrain( i0, i1, lmin=1, lmax=1, kmin=1, kmax=1, flim=1e+300, l=None, k=None, shift=(0.,0.,0.), bOldIndex=True ):
    if l is not None:
        lmin=l; lmax=l
    if k is not None:
        kmin=k; kmax=k
    shift=np.array(shift)
    return lib.addDistConstrain(i0, i1, lmin, lmax, kmin, kmax, flim, _np_as(shift,c_double_p) )

#  void addAngConstrain(  int i0,int i1,int i2, double ang0, double k ){
lib.addAngConstrain.argtypes  = [c_int, c_int, c_int, c_double, c_double]
lib.addAngConstrain.restype   =  None
def addAngConstrain(i0, i1, i2, ang0=0.0, k=1.0):
    return lib.addAngConstrain(i0, i1, i2, ang0, k)

#  int substituteMolecule( const char* fname, int ib, double* up, int ipivot, bool bSwapBond )
lib.substituteMolecule.argtypes  = [c_char_p, c_int, c_double_p, c_int, c_bool]
lib.substituteMolecule.restype   =  c_int
def substituteMolecule(fname, ib, ipivot, up=(0,0,1), bSwapBond=False):
    up = np.array(up, np.int23)
    return lib.substituteMolecule( cstr(fname), ib, _np_as(up,c_double_p), ipivot, bSwapBond)

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

# #  int splitAtBond( int ib, int* sel )
# lib.splitAtBond.argtypes  = [c_int, c_int_p]
# lib.splitAtBond.restype   =  c_int
# def splitAtBond(ib, sel=None):
#     return lib.splitAtBond(ib, _np_as(sel,c_int_p))

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


# =====================================
# ========= Database
# =====================================

def addSnapshot(ifNew=False, fname=None):
    lib.addSnapshot(ifNew, cstr(fname))

def printDatabase():
    lib.printDatabase()

def computeDistance(ia,ib,dist=None):
    if dist is None: dist=np.zeros(1)
    return lib.computeDistance(ia,ib,_np_as(dist,c_double_p))


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
