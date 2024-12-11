
import numpy as np
from   ctypes import c_int, c_double, c_bool, c_float, c_char_p, c_bool, c_void_p, c_char_p
import ctypes
import os
import sys

#sys.path.append('../')
#from pyMeta import cpp_utils 
from . import cpp_utils_ as cpp_utils
#import cpp_utils_ as cpp_utils

plt = None

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

bWeightsSet = False
#isInitialized = False
#nfound = -1

# ====================================
# ========= C functions
# ====================================

def cstr( s ):
    if s is None: return None
    return s.encode('utf8')

#  void setVerbosity( int verbosity_, int idebug_, int PrintDOFs, int PrintfDOFs, int PrintBeforReg, int PrintAfterReg ){
lib.setVerbosity.argtypes  = [c_int, c_int, c_int, c_int, c_int, c_int]
lib.setVerbosity.restype   =  None
def setVerbosity(verbosity=1, idebug=0, PrintDOFs=0, PrintfDOFs=0, PrintBeforReg=0, PrintAfterReg=0):
    return lib.setVerbosity(verbosity, idebug, PrintDOFs, PrintfDOFs, PrintBeforReg, PrintAfterReg)

# void setup( int imodel, int EvalJ, int WriteJ, int CheckRepulsion, int Regularize, int AddRegError, int Epairs, int BroadcastFDOFs, int UdateDOFbounds){
lib.setup.argtypes  = [c_int,c_int, c_int, c_int, c_int, c_int, c_int, c_int, c_int]
lib.setup.restype   =  None    
def setup(imodel=1, EvalJ=0, WriteJ=0, CheckRepulsion=0, Regularize=0, AddRegError=0, Epairs=0, BroadcastFDOFs=0, UdateDOFbounds=0):
    return lib.setup( imodel,EvalJ, WriteJ, CheckRepulsion, Regularize, AddRegError, Epairs, BroadcastFDOFs, UdateDOFbounds)

#void setFilter( double EmodelCut, double EmodelCutStart, int iWeightModel, int ListOverRepulsive, int SaveOverRepulsive, int PrintOverRepulsive, int DiscardOverRepulsive, int WeightByEmodel ){
lib.setFilter.argtypes  = [c_double, c_double, c_int, c_int, c_int, c_int, c_int, c_int]
lib.setFilter.restype   =  None
def setFilter( EmodelCut=1.0, EmodelCutStart=None, EmodelCutFactor=0.75, iWeightModel=1, ListOverRepulsive=0, SaveOverRepulsive=0, PrintOverRepulsive=0, DiscardOverRepulsive=0, WeightByEmodel=0 ):
    if EmodelCutStart is None: EmodelCutStart = EmodelCut*EmodelCutFactor
    return lib.setFilter( EmodelCut, EmodelCutStart, iWeightModel, ListOverRepulsive, SaveOverRepulsive, PrintOverRepulsive, DiscardOverRepulsive, WeightByEmodel )

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
    global bWeightsSet
    bWeightsSet = True
    n = len(weights)
    lib.setWeights(n, _np_as(weights,c_double_p))
    
# double run( int ialg, int iparallel, int nstep, double Fmax, double dt, double max_step, double damping, bool bClamp ){
lib.run.argtypes  = [c_int, c_int, c_int, c_double, c_double, c_double, c_double, c_bool]
lib.run.restype   =  c_double
def run(ialg=2, iparallel=1, nstep=100, Fmax=1e-8, dt=0.01, max_step=0.05, damping=0.0, bClamp=False):
    return lib.run( ialg, iparallel, nstep, Fmax, dt, max_step, damping, bClamp )

# double getEs( double* Es, double* Fs, bool bOmp, bool bDOFtoTypes ){
lib.getEs.argtypes  = [ c_double_p,  c_double_p, c_bool, c_bool]
lib.getEs.restype   =  c_double
def getEs(Es=None, Fs=None, bOmp=False, bDOFtoTypes=True, bEs=True, bFs=False ):
    if bEs and (Es is None): Es = np.zeros( nbatch )
    if bFs and (Fs is None): Fs = np.zeros( nDOFs  )
    Eerr = lib.getEs(_np_as(Es,c_double_p), _np_as(Fs,c_double_p), bOmp, bDOFtoTypes)
    #print("getEs(): Es", Es)
    return Eerr, Es, Fs

# void scanParam( int iDOF, int n, double* xs,  double* Es, double* Fs, bool bRegularize, bool bEvalSamples ){
lib.scanParam.argtypes  = [c_int, c_int, c_double_p, c_double_p, c_double_p, c_bool]
lib.scanParam.restype   = None
def scanParam( iDOF, xs, Es=None, Fs=None, bEvalSamples=True ):
    n = len(xs)
    if Es is None: Es = np.zeros( n )
    if Fs is None: Fs = np.zeros( n )
    lib.scanParam(iDOF, n, _np_as(xs,c_double_p), _np_as(Es,c_double_p), _np_as(Fs,c_double_p), bEvalSamples )
    return Es,Fs

# void scanParam2D( int iDOFx, int iDOFy, int nx, int ny, double* xs, double* ys,  double* Es, double* Fx, double* Fy, bool bRegularize, bool bEvalSamples ){
lib.scanParam2D.argtypes  = [c_int, c_int, c_int, c_int, c_double_p, c_double_p, c_double_p, c_double_p, c_double_p, c_bool]
lib.scanParam2D.restype   = None
def scanParam2D(iDOFx, iDOFy, xs, ys, Es=None, Fx=None, Fy=None, bEvalSamples=True):
    nx, ny = len(xs), len(ys)
    if Es is None: Es = np.zeros((ny,nx))
    if Fx is None: Fx = np.zeros((ny,nx))
    if Fy is None: Fy = np.zeros((ny,nx))
    lib.scanParam2D(iDOFx, iDOFy, nx, ny, _np_as(xs,c_double_p), _np_as(ys,c_double_p),    _np_as(Es,c_double_p), _np_as(Fx,c_double_p), _np_as(Fy,c_double_p), bEvalSamples)
    return Es, Fx, Fy

# void loadTypes_new( const char* fname_ElemTypes, const char* fname_AtomTypes ){
lib.loadTypes.argtypes  = [c_char_p, c_char_p]
lib.loadTypes.restype   =  None
def loadTypes(fEtypes="data/ElementTypes.dat", fAtypes="data/AtomTypes.dat"):
    return lib.loadTypes(cstr(fEtypes), cstr(fAtypes))

#int loadDOFSelection( const char* fname ){
lib.loadDOFSelection.argtypes  = [c_char_p]
lib.loadDOFSelection.restype   =  c_int
def loadDOFSelection(fname="DOFSelection.dat"):
    return lib.loadDOFSelection(cstr(fname))

# int loadTypeSelection_walls( const char* fname ){
# lib.loadTypeSelection.argtypes  = [c_char_p]
# lib.loadTypeSelection.restype   =  c_int
# def loadTypeSelection(fname="typeSelection.dat"):
#     return lib.loadTypeSelection(cstr(fname))

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
    global DOFs,fDOFs,vDOFs, fDOFbounds
    DOFs          = getBuff ( "DOFs",   (nDOFs,)  )
    fDOFs         = getBuff ( "fDOFs",  (nDOFs,)  ) 
    vDOFs         = getBuff ( "vDOFs",  (nDOFs,)  ) 
    fDOFbounds    = getBuff ( "fDOFbounds", (nDOFs,2)  ) 
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
    if bWeightsSet:
        global weights
        weights = getBuff ( "weights",    (nbatch,) )
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



def check_array_difference(arr1, arr2, name, max_error=1e-8, err_message="arrays differs" ):
    dmax = (arr1-arr2).max()
    print(f"{name} dmax={dmax}")
    if not np.allclose(arr1, arr2, atol=max_error):
        print(f"{name} arrays differ:")
        for i, (v1, v2) in enumerate(zip(arr1, arr2)):
            if not np.isclose(v1, v2): print(f"{i}\t{v1}\t{v2}")
        assert False, f"{name} "+err_message

def test_getEs_openmp():
    E1, Es1, Fs1 = fit.getEs(imodel=imodel, bOmp=False, bEs=True, bFs=True)
    E2, Es2, Fs2 = fit.getEs(imodel=imodel, bOmp=True,  bEs=True, bFs=True)
    check_array_difference(Es1, Es2, "Es w/o OpenMP")
    check_array_difference(Fs1, Fs2, "Fs w/o OpenMP")
    print( "test_getEs_openmp() E1,E2", E1, E2 )
    assert np.isclose(E1, E2), f"E values differ: E1={E1}, E2={E2}"
    print( "test_getEs_openmp() passed " )

def read_xyz_data(fname="input_all.xyz"):
    """Read XYZ file and extract Etot and x0 values from comment lines"""
    #print("read_xyz_data()\n")
    #print("Reading XYZ file:", fname)
    Erefs = []
    x0s = []
    with open(fname, 'r') as f:
        while True:
            line = f.readline()
            #print(line)
            if not line: break
            if line.startswith('# n0'):
                #print(line)
                # Parse line like "# n0 5 Etot .70501356708840164618 x0 1.40"
                parts = line.split()
                Etot  = float(parts[4])
                x0    = float(parts[6])
                Erefs.append(Etot)
                x0s.append(x0)
            # Skip the rest of the xyz structure
            #natoms = int(line) if line[0].isdigit() else 0
            #for _ in range(natoms):
            #    f.readline()
    return np.array(Erefs), np.array(x0s)

def split_and_weight_curves(Erefs, x0s, n_before_min=4, weight_func=None ):
    """
    Split energy curves based on x0 discontinuities and assign weights.
    
    Args:
        Erefs: numpy array of total energies
        x0s: numpy array of x0 values (monotonic within each curve)
        n_before_min: number of points before minimum to keep with positive weight
    
    Returns:
        weights: numpy array of weights (0.0 or 1.0)
    """
    weights = np.zeros_like(Erefs)
    
    # Find where x0 values reset (non-monotonic changes)
    dx0 = np.diff(x0s)
    curve_starts = np.where(dx0 < 0)[0] + 1
    
    # Add start and end indices to process all segments
    all_splits = np.concatenate(([0], curve_starts, [len(x0s)]))
    
    # Process each curve segment
    for start, end in zip(all_splits[:-1], all_splits[1:]):
        segment = Erefs[start:end]
        if len(segment) == 0:
            continue
        imin = np.argmin(segment) + start
        icut = imin - n_before_min
        weight_start = max(icut, start)
        if weight_func is None:
            weights[weight_start:end] = 1.0
        else:
            weights[weight_start:end] = weight_func(Erefs[weight_start:end])
    
    return weights

def exp_weight_func(Erefs, a=1.0, alpha=3.0 ):
    Emin = np.min(Erefs)
    return np.exp( alpha*(Erefs-Emin)/Emin )*a

def genWeights(Erefs, Ecut ):
    mask = Erefs<Ecut
    weights = np.zeros( len(Erefs) )
    weights[mask] = 1.0
    return weights

def plotEWs(Erefs=None,Emodel=None,  weights=None,  weights0=None, bLimEref=True, Emin=None, EminFac=1.2):
    plt.figure( figsize=(15,5) )
    if( Erefs    is not None): plt.plot( Erefs   ,'.k' , lw=0.5, ms=1.0, label="E_ref")
    if( Emodel   is not None): plt.plot( Emodel  ,'.-r', lw=0.5, ms=1.0, label="E_model")
    if( weights  is not None): plt.plot( weights ,'-g' , lw=0.5,         label="weights")
    if( weights0 is not None): plt.plot( weights0,'--c', lw=0.5,         label="weights0")
    plt.legend()
    plt.xlabel("#sample(conf)")
    plt.ylabel("E [kcal/mol]")
    plt.grid()
    if Emin is not None:
        plt.ylim( Emin*EminFac, -Emin*EminFac )
    elif bLimEref:
        Emin =  Erefs.min()
        plt.ylim( Emin*EminFac, -Emin*EminFac )

def plotDOFscans( iDOFs, xs, DOFnames, bEs=True, bFs=False,  label="plotDOFscans", bEvalSamples=False ):
    plt.figure()
    for iDOF in iDOFs:
        y = DOFs[iDOF]    # store backup value of this DOF
        Es,Fs = scanParam( iDOF, xs, bEvalSamples=bEvalSamples )   # do 1D scan
        #print( "iDOF", iDOF, DOFnames[iDOF], "Es", Es )
        if bEs: plt.plot(xs,Es, '-', label=DOFnames[iDOF] )       # plot 1D scan
        if bFs: plt.plot(xs,Fs, '-', label=DOFnames[iDOF] )       # plot 1D scan
        DOFs[iDOF] = y    # restore
    plt.legend()
    plt.xlabel("DOF value")
    plt.ylabel("E [kcal/mol]")    
    plt.title( label )
    plt.grid()
    plt.show()






