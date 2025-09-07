#from nt import write
import numpy as np
from   ctypes import c_int, c_double, c_bool, c_float, c_char_p, c_bool, c_void_p, c_char_p
import ctypes
import os
import sys
#import glob

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
ev2kcal = 23.060548
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

# void setGlobalParams( double kMorse, double Lepairs, double EijMax, double softClamp_start, double softClamp_max ){
lib.setGlobalParams.argtypes  = [c_double, c_double, c_double, c_double, c_double]
lib.setGlobalParams.restype   =  None
def setGlobalParams(kMorse=1.6, Lepairs=0.5, EijMax=5.0, softClamp_start=4.0, softClamp_max=6.0):
    return lib.setGlobalParams(kMorse, Lepairs, EijMax, softClamp_start, softClamp_max)

# void setup( int imodel, int EvalJ, int WriteJ, int CheckRepulsion, int Regularize, int AddRegError, int Epairs, int BroadcastFDOFs, int UdateDOFbounds, int EvalOnlyCorrections, int SaveJustElementXYZ){
lib.setup.argtypes  = [c_int,c_int, c_int, c_int, c_int, c_int, c_int, c_int, c_int, c_int]
lib.setup.restype   =  None    
def setup(imodel=1, EvalJ=0, WriteJ=0, CheckRepulsion=0, Regularize=0, AddRegError=0, Epairs=0, BroadcastFDOFs=0, UdateDOFbounds=0, EvalOnlyCorrections=0, SaveJustElementXYZ=0):
    return lib.setup( imodel,EvalJ, WriteJ, CheckRepulsion, Regularize, AddRegError, Epairs, BroadcastFDOFs, UdateDOFbounds, EvalOnlyCorrections, SaveJustElementXYZ)

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
    
#void setTrjBuffs( double* trj_E, double* trj_F, double* trj_DOFs, double* trj_fDOFs){
lib.setTrjBuffs.argtypes  = [c_double_p, c_double_p, c_double_p, c_double_p]
lib.setTrjBuffs.restype   =  None
def setTrjBuffs( niter, trj_E=None, trj_F=None, trj_DOFs=None, trj_fDOFs=None, nDOFs_=None, bE=True, bF=True, bDOFs=True, bfDOFs=False ):
    if nDOFs_     is None: nDOFs_ = nDOFs
    print("setTrjBuffs(): niter=%i nDOFs=%i bE=%i bF=%i bDOFs=%i bfDOFs=%i" % (niter, nDOFs_, bE, bF, bDOFs, bfDOFs))
    if (trj_E     is None) and bE     : trj_E     = np.zeros( niter )
    if (trj_F     is None) and bF     : trj_F     = np.zeros( niter )
    if (trj_DOFs  is None) and bDOFs  : trj_DOFs  = np.zeros( (niter, nDOFs_) )
    if (trj_fDOFs is None) and bfDOFs : trj_fDOFs = np.zeros( (niter, nDOFs_) )
    lib.setTrjBuffs(_np_as(trj_E,c_double_p), _np_as(trj_F,c_double_p), _np_as(trj_DOFs,c_double_p), _np_as(trj_fDOFs,c_double_p))
    return trj_E, trj_F, trj_DOFs, trj_fDOFs

# double run( int ialg, int iparallel, int nstep, double Fmax, double dt, double max_step, double damping, bool bClamp ){
lib.run.argtypes  = [c_int, c_int, c_int, c_double, c_double, c_double, c_double, c_bool]
lib.run.restype   =  c_double
def run(ialg=2, iparallel=1, nstep=100, Fmax=1e-8, dt=0.01, max_step=0.05, damping=0.0, bClamp=False):
    return lib.run( ialg, iparallel, nstep, Fmax, dt, max_step, damping, bClamp )

# double getEs( double* Es, double* Fs, bool bOmp, bool bDOFtoTypes, char* xyz_name){
lib.getEs.argtypes  = [ c_double_p,  c_double_p, c_bool, c_bool, c_char_p ]
lib.getEs.restype   =  c_double
def getEs(Es=None, Fs=None, bOmp=False, bDOFtoTypes=True, bEs=True, bFs=False, xyz_name=None  ):
    if xyz_name is not None: 
        with open(xyz_name,'w') as f: f.write("")
    if bEs and (Es is None): Es = np.zeros( nbatch )
    if bFs and (Fs is None): Fs = np.zeros( nDOFs  )
    Eerr = lib.getEs(_np_as(Es,c_double_p), _np_as(Fs,c_double_p), bOmp, bDOFtoTypes, cstr(xyz_name) )
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

# ------------------------------------
# Sampling interface for damped Coulomb functions (C++ -> Python)
# void sample_funcEF( int n, double* xs, double* EFs_, int kind, double* params )
lib.sample_funcEF.argtypes = [c_int, c_double_p, c_double_p, c_int, c_double_p]
lib.sample_funcEF.restype  = None
def sample_funcEF(xs, kind=0, params=None):
    """
    Sample damped Coulomb functions from C++ for given radii xs.

    Args:
      xs: 1D array of radii r
      kind:
        0   -> bare 1/r
        10  -> Boys exact erf(r)/r
        11  -> Boys C1 Hermite Cubic
        12  -> Boys C2 Hermite Quintic
        13  -> Boys C1 Quartic Even
        14  -> Boys C2 Sextic Even
        20  -> Soft clamp (positive)
        21  -> Smooth clamp (positive)
        22  -> Soft clamp (negative)
        23  -> Smooth clamp (negative)
      params:
        for Boys kinds: [rmin]
        for clamp kinds: [y1, y2]

    Returns:
      E, F arrays where E=y(r) and F=-dE/dr (radial force)
    """
    xs = np.asarray(xs, dtype=np.double)
    n  = xs.size
    EFs = np.zeros((n,2), dtype=np.double)
    if params is not None:
        params = np.asarray(params, dtype=np.double)
        pptr = _np_as(params, c_double_p)
    else:
        pptr = _np_as(None, c_double_p)
    lib.sample_funcEF(n, _np_as(xs,c_double_p), _np_as(EFs,c_double_p), kind, pptr)
    return EFs[:,0].copy(), EFs[:,1].copy()

# void loadTypes_new( const char* fname_ElemTypes, const char* fname_AtomTypes ){
lib.loadTypes.argtypes  = [c_char_p, c_char_p]
lib.loadTypes.restype   =  None
def loadTypes(fEtypes="data/ElementTypes.dat", fAtypes="data/AtomTypes.dat"):
    return lib.loadTypes(cstr(fEtypes), cstr(fAtypes))

#int loadDOFSelection( const char* fname ){
lib.loadDOFSelection.argtypes  = [c_char_p]
lib.loadDOFSelection.restype   =  c_int
def loadDOFSelection(fname="dofSelection.dat"):
    global nDOFs
    nDOFs = lib.loadDOFSelection(cstr(fname))
    return nDOFs

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

# int loadXYZ( const char* fname, bool bAddEpairs, bool bOutXYZ, bool bEvalOnlyCorrections, bool bAppend ){ 
lib.loadXYZ.argtypes  = [c_char_p, c_bool, c_bool, c_bool, c_bool]
lib.loadXYZ.restype   =  c_int
def loadXYZ(fname, bAddEpairs=False, bOutXYZ=False, bEvalOnlyCorrections=False, bAppend=False ):
    global nbatch
    nbatch = lib.loadXYZ(cstr(fname), bAddEpairs, bOutXYZ, bEvalOnlyCorrections, bAppend)
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

def find_all_dirs(base_path):
    return [ item for item in os.listdir(base_path) if os.path.isdir(os.path.join(base_path, item)) ]

def extract_fragment_names( directories=None, base_path='./'):
    if directories is None:
        directories  = find_all_dirs(base_path)
    first_molecules  = set()
    second_molecules = set()
    for directory in directories:
        parts = directory.split("_")
        if len(parts) == 2:
            first_molecules .add(parts[0])
            second_molecules.add(parts[1])
        else:
            print(f"Warning: Directory name {directory} does not contain two molecules")
    sorted_first  = sorted(list(first_molecules))
    sorted_second = sorted(list(second_molecules))
    return sorted_first, sorted_second

def combine_fragments( frags1, frags2, path=None, ext=""):
    dirs = [ f"{f1}_{f2}" for f1 in frags1 for f2 in frags2 ]
    if path is not None:
        dirs_ = [ d for d in dirs if os.path.exists( os.path.join(path,d+ext) ) ]
        return dirs_
    else:
        return dirs

def concatenate_xyz_files(directories=None, base_path='./', fname="all.xyz", output_file="all.xyz", mode='w'):
    """
    Concatenates all.xyz files from specified directories, adding the directory name
    to the comment lines (starting with '#'), and writes directly to an output file.

    Args:
        directories (list): A list of directory names to process.
        base_path (str): The base path to the directories containing the all.xyz files.
        output_file (str): The path to the output file.
    """
    if directories is None:
        directories = find_all_dirs(base_path)
    marks = []
    with open(output_file, mode) as outfile:
        i=0
        for directory in directories:
            dir_path  = os.path.join(base_path, directory)
            file_path = os.path.join(dir_path, fname)
            if not os.path.exists(file_path):
                print(f"Warning in concatenate_xyz_files(): file {fname} not found in {directory}. Skipping it...")
                continue
            with open(file_path, 'r') as infile:
                i0 = i
                for line in infile:
                    if line.startswith("#"):
                        new_line = f"{line.strip()} {directory}\n"
                        outfile.write(new_line)
                        i+=1
                    else:
                        outfile.write(line)
                i1 = i
            marks.append( (i0,i1) )
    return marks

def concatenate_xyz_files_flat(names=None, base_path='./', output_file="all.xyz", mode='w'):
    """
    Concatenates xyz files from a flat directory structure, where xyz files are directly in the path.
    Adds the filename (without extension) to the comment lines.

    Args:
        names (list): A list of names to process. If None, will find all xyz files in base_path.
        base_path (str): The path containing the xyz files.
        output_file (str): The path to the output file.
        mode (str): Write mode for the output file ('w' for write, 'a' for append).
    """
    if names is None:
        # Find all xyz files in the directory
        names = []
        for f in os.listdir(base_path):
            if f.endswith('.xyz') and f != output_file:
                names.append(os.path.splitext(f)[0])
    
    marks = []
    with open(output_file, mode) as outfile:
        i = 0
        for name in names:
            file_path = os.path.join(base_path, f"{name}.xyz")
            if not os.path.exists(file_path):
                print(f"Warning in concatenate_xyz_files_flat(): file {name}.xyz not found in {base_path}. Skipping it...")
                continue
            with open(file_path, 'r') as infile:
                i0 = i
                for line in infile:
                    if line.startswith("#"):
                        new_line = f"{line.strip()} {name}\n"
                        outfile.write(new_line)
                        i += 1
                    else:
                        outfile.write(line)
                i1 = i
            marks.append((i0, i1))
    return marks

def read_file_comments(fname, comment_sign='#'):
    with open(fname, 'r') as f:
        return [ line.strip() for line in f if line.startswith(comment_sign) ]

def extract_comments_and_types(fname, comment_sign='#'):
    type_names = set()  # Use a set to avoid duplicates
    comment_lines = []
    i=0
    with open(fname, 'r') as file:
        for line in file:
            if line.startswith(comment_sign):
                comment_lines.append(line)
            else:
                parts = line.split()
                if parts: 
                    if "_" in parts[0]:
                        #parts = parts[0].split("_")
                        type_names.add(parts[0])
                #print(parts)
                #i+=1
            #if i>10: break
    return type_names, comment_lines

def add_epair_types(types):
    epair_types = []
    for t1 in types:
        # copy type name without "_" (N_3 -> N3)
        epair_types.append( "E_"+t1.replace("_","") )
    # add to set
    types.update(epair_types)
    return types


def comment_non_matching_lines( type_names, fname_in, bWriteAsComment=False, fname_out="dofSelection.dat"):
    with open(fname_in, 'r') as file:
        lines = file.readlines()
    with open(fname_out, 'w') as file:
        for line in lines: 
            if line.startswith('#'):
                file.write(line)
                continue
            parts = line.split()
            if parts:
                typename = parts[0]  # The first column is the type name
                if typename not in type_names:
                    if bWriteAsComment:
                        file.write(f"# {line.strip()}\n")  # Comment out the line
                else:
                    file.write(line)  # Keep the line as is

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

def loadDOFnames( fname, comps="REQH" ):
    names = []
    with open(fname, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split()
            icomp = int(parts[1])
            names.append( parts[0]+"."+ comps[icomp] )
    return names

def mark_molecule_blocks(lines, ):
    marks = []
    lens  = []
    #current_mark = 0
    prev_angle  = None
    prev_mol    = None
    angle_data  = []
    nxs   = []   #store number of nx for each angle (iy)
    i0 = 0
    i1 = -1
    nx = 0
    for i,line in enumerate(lines):
        i1 = i
        parts = line.split()
        #E     = float(parts[4]) 
        #r     = float(parts[6]) 
        angle = float(parts[8])  # Extract angle after 'z'
        mol   = parts[-1].strip()  # Last word is the molecule pair
        #print( i, E, angle, r, mol, line )
        if prev_mol is None: prev_mol = mol
        if mol != prev_mol:
            angle_data.append( nxs )
            nxs = []
            marks.append( (i0,i1) )
            i0 = i1
            #print( len(marks), (i0,i1), mol, prev_mol, angle )
            prev_mol = mol
        if prev_angle is None: prev_angle = angle
        if angle != prev_angle:
            nxs.append( nx )
            nx=0
            prev_angle = angle
        nx += 1
    # Append the last molecule pair
    if prev_angle is not None:
        marks.append( (i0,i1) )
        angle_data.append( nxs )
    return marks, angle_data

def slice_and_reshape(Es, marks, angle_data):
    Eplots = []
    nmol   = len(marks)
    if nmol != len(angle_data):
        print(f"Error: len(marks={nmol}) != len(angle_data={len(angle_data)})")
        exit(0)
    for i in range( nmol ):
        i0,i1    = marks[i]
        nxs      = angle_data[i]
        nx       = np.max(nxs)
        ny       = len(nxs) 
        Eplot    = np.empty( (ny,nx) )
        Eplot[:,:] = np.nan
        for iy in range(ny):
            nxi             = nxs[iy]
            #print(  len(Es), "i0:i0+nxi ", i0,i0+nxi )
            Eplot[iy,0:nxi] = Es[i0:i0+nxi]
            i0+=nxi
        Eplots.append(Eplot)
    return Eplots

def split_and_weight_curves(Erefs, x0s, n_before_min=4, weight_func=None, EminMin=-0.02 ):
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

    lens = []
    # Process each curve segment
    for start, end in zip(all_splits[:-1], all_splits[1:]):
        segment = Erefs[start:end]
        n = len(segment)
        lens.append(n)
        #print(  )
        if len(segment) == 0: continue
        imin = np.argmin(segment) + start

        Emin = np.min(segment)
        if Emin < EminMin:
            icut = imin - n_before_min
            weight_start = max(icut, start)
            if weight_func is None:
                weights[weight_start:end] = 1.0
            else:
                weights[weight_start:end] = weight_func(Erefs[weight_start:end])
        else:
            print("Emin=", Emin)
            if weight_func is None:
                weights[start:end] = 1.0
            else:
                weights[start:end] = weight_func( Erefs[start:end] )
                #weights[start:end] = 1.0
    
    return weights, lens

def exp_weight_func(Erefs, a=1.0, alpha=3.0, Emin0=0.1 ):
    Emin = np.min(Erefs)
    return np.exp( -alpha*(Erefs-Emin)/( np.abs(Emin) + Emin0 ) )*a

def genWeights(Erefs, Ecut ):
    mask = Erefs<Ecut
    weights = np.zeros( len(Erefs) )
    weights[mask] = 1.0
    return weights

def plotEWs(Erefs=None,Emodel=None,  weights=None,  weights0=None, bLimEref=True, Emin=None, EminFac=1.2 , bKcal=True):
    E_units = 1.0
    if bKcal: 
        E_units = ev2kcal
        units="[kcal/mol]"
    else:
        units="[eV]" 
    plt.figure( figsize=(15,5) )
    if( Erefs    is not None): plt.plot( Erefs*E_units   ,'.k' , lw=0.5, ms=1.0, label="E_ref")
    #if( Emodel   is not None): plt.plot( Emodel*E_units  ,'.-r', lw=0.5, ms=1.0, label="E_model")
    if( Emodel   is not None): plt.plot( Emodel*E_units  ,'-r', lw=0.5, ms=1.0, label="E_model")
    if( weights  is not None): plt.plot( weights ,'-g' , lw=0.5,         label="weights")
    if( weights0 is not None): plt.plot( weights0,'--c', lw=0.5,         label="weights0")
    plt.legend()
    plt.xlabel("#sample(conf)")
    plt.ylabel("E "+units)
    plt.grid()
    if Emin is not None:
        plt.ylim( Emin*EminFac, -Emin*EminFac )
    elif bLimEref:
        Emin =  Erefs.min()*EminFac
        plt.ylim( Emin*EminFac, -Emin*EminFac )

def numDeriv( x, y):
    dx = x[2:]-x[:-2]
    dy = y[2:]-y[:-2]
    xs = x[1:-1]
    return dy/dx, xs

def plotDOFscans( iDOFs, xs, DOFnames, bEs=True, bFs=False,  title="plotDOFscans", bEvalSamples=True, bPrint=False ):
    plt.figure(figsize=(8,10.0))
    color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
    ncol = len(color_cycle)
    for iDOF in iDOFs:
        y = DOFs[iDOF]    # store backup value of this DOF
        #print( f"#======= DOF[{iDOF}]: {xs}" )
        Es,Fs = scanParam( iDOF, xs, bEvalSamples=bEvalSamples )   # do 1D scan
        #print( f"#======= fDOF[{iDOF}]: {Efs}" )
        #print( "iDOF", iDOF, DOFnames[iDOF], "Es", Es )
        # take color from standard matplotlib color cycle
        c = color_cycle[iDOF % ncol]
        if bEs: 
            plt.subplot(2,1,1); plt.plot(xs,Es, '-', color=c, label="E "+DOFnames[iDOF] )       # plot 1D scan
        if bFs: 
            Fs_num, xs_num = numDeriv(xs,Es)
            plt.subplot(2,1,2); plt.plot(xs,Fs,    '-', lw=1.0, color=c, label="F "+DOFnames[iDOF] )       # This is error in the E_O3 charge derivative
            plt.subplot(2,1,2); plt.plot(xs_num,-Fs_num, ':', lw=1.5, color=c, label="F "+DOFnames[iDOF] ) 
            if bPrint:
                print ( "# plotDOFscans DOF ", iDOF, DOFnames[iDOF], " dx= ", xs[1]-xs[0] ); 
                print ( "#  i          x              E            F_ana=-dE/dDOF     F_num          F_ana-F_num          // F_num=-(E[i+1]-E[i-1])/(x[i+1]-x[i-1])" ); 
                for i in range(1, len(xs)-1 ): print( " %3i %15.5f %15.5f %15.5f %15.5f %15.5f"%(i, xs[i], Es[i], Fs[i], -Fs_num[i-1], Fs[i]+Fs_num[i-1]) )
        DOFs[iDOF] = y    # restore

    plt.subplot(2,1,1);
    plt.legend()
    plt.xlabel("DOF value")
    plt.ylabel("E [kcal/mol]")   
    plt.grid()

    plt.subplot(2,1,2);
    plt.xlabel("DOF value")
    plt.ylabel("F [kcal/mol/A]")   
    plt.grid()

    plt.suptitle( title )


def checkDOFderiv( iDOF, x0=0.5, d=0.001, bEvalSamples=True ):
        xs = np.array([x0-d,x0,x0+d])
        Es,Fs = scanParam( iDOF, xs, bEvalSamples=bEvalSamples )   # do 1D scan
        F_num = numDeriv(xs,Es)[0][0]
        print( f"checkDOFderiv(iDOF={iDOF}) x={x0} d={d} F_ana=-dE/dx={Fs[1]} F_num=-(E(x+d)-E(x-d))/(2d)={-F_num} F_ana-F_num={Fs[1]+F_num}" )

def plot_Epanels(Eplots, ref_dirs, bColorbar=True, Emin=-5.0, bKcal=False ):
    E_units = 1.0
    if bKcal: 
        E_units = ev2kcal
    nmols = len(Eplots)
    if nmols != len(ref_dirs):
        print(f"Error: len(Eplots={nmols}) != len(ref_dirs={len(ref_dirs)})")
        return
    fig, axs = plt.subplots( 1,nmols, figsize=(20, 3))  # 3 rows for Erefs, Emodel, Ediff
    for i in range(nmols):
        # Plot Reference Energies
        im = axs[i].imshow(Eplots[i].T*E_units, aspect='auto', origin='lower', vmin=Emin, vmax=-Emin, cmap='bwr' )
        axs[i].set_title(f"Ref: {ref_dirs[i]}")
        axs[i].set_ylabel('Reference Energies')
        if bColorbar:  plt.colorbar(im, ax=axs[i])
    plt.tight_layout()
    return fig

def plot_Epanels_diff(Emodels, Erefs, ref_dirs, bColorbar=True, Emin=-5.0, bKcal=False ):
    E_units = 1.0
    if bKcal: 
        E_units = ev2kcal
    nmols = len(Erefs)
    if nmols != len(ref_dirs):
        print(f"Error: len(Erefs={nmols}) != len(ref_dirs={len(ref_dirs)})")
        return
    fig, axs = plt.subplots( 3,nmols, figsize=(nmols*4, 6))  # 3 rows for Erefs, Emodel, Ediff
    for i in range(nmols):
        # Check if axs is 1D or 2D
        if nmols == 1:
            # Use single index for axs
            ax_ref   = axs[0]
            ax_model = axs[1]
            ax_diff  = axs[2]
        else:
            # Use double indexing for axs
            ax_ref   = axs[0, i]
            ax_model = axs[1, i]
            ax_diff  = axs[2, i]

        Emodel = Emodels[i]
        Eref   = Erefs[i]
        
        # Plot Reference Energies
        im = ax_ref.imshow(Eref.T*E_units, aspect='auto', origin='lower', vmin=Emin, vmax=-Emin, cmap='bwr' )
        ax_ref.set_ylabel('Reference Energies')
        if bColorbar:  plt.colorbar(im, ax=ax_ref)
        ax_ref.set_title(f"Ref: {ref_dirs[i]}")
    
        # Plot Model Energies
        im = ax_model.imshow(Emodel.T*E_units, aspect='auto', origin='lower', vmin=Emin, vmax=-Emin, cmap='bwr' )
        ax_model.set_ylabel('Model Energies')
        if bColorbar: plt.colorbar(im, ax=ax_model)

        # Plot Difference
        Ediff = Emodel - Eref  # Calculate difference
        im = ax_diff.imshow(Ediff.T*E_units, aspect='auto', origin='lower', vmin=Emin, vmax=-Emin, cmap='bwr' )
        ax_diff.set_ylabel('Error ')
        if bColorbar: plt.colorbar(im, ax=ax_diff)
    plt.tight_layout()
    return fig

def plot2Dlong_diff( Erefs, Es, lens  ):
    lens = np.array(lens)
    nmax = np.max(lens) 
    nseg = len(lens)
    Eplot  = np.zeros( (nseg, nmax)  ); Eplot[:,:]  = np.nan
    Eplot_ = np.zeros( (nseg, nmax)  ); Eplot_[:,:] = np.nan
    ii = 0
    for i in range( len(lens) ):
        ni = lens[i]
        #Eplot[i,0:lens[i]] = Es[i]
        Eplot [i,0:ni] = Erefs[ii:ii+ni]
        Eplot_[i,0:ni] = Es   [ii:ii+ni]
        ii+=ni
    dEplot = Eplot - Eplot_
    plt.figure(figsize=(20,12))
    Emin = np.min(Erefs)
    dEmax = max( -np.nanmin(dEplot),np.nanmax(dEplot) )
    dEmax = 0.1
    print( "dEmax: ", dEmax, "Emin ", Emin )
    #Emax = np.max(Eplot)
    plt.subplot(6,1,1); plt.imshow( Eplot [:nseg//2,:].T, origin='lower', vmin=Emin,   vmax=-Emin, cmap='bwr', extent=[ 0, len(lens), 0, 180.0 ] )
    plt.subplot(6,1,2); plt.imshow( Eplot_[:nseg//2,:].T, origin='lower', vmin=Emin,   vmax=-Emin, cmap='bwr', extent=[ 0, len(lens), 0, 180.0 ] )
    plt.subplot(6,1,3); plt.imshow( dEplot[:nseg//2,:].T, origin='lower', vmin=-dEmax, vmax=dEmax, cmap='bwr', extent=[ 0, len(lens), 0, 180.0 ] )
    plt.subplot(6,1,4); plt.imshow( Eplot [nseg//2:,:].T, origin='lower', vmin=Emin,   vmax=-Emin, cmap='bwr', extent=[ 0, len(lens), 0, 180.0 ] )
    plt.subplot(6,1,5); plt.imshow( Eplot_[nseg//2:,:].T, origin='lower', vmin=Emin,   vmax=-Emin, cmap='bwr', extent=[ 0, len(lens), 0, 180.0 ] )
    plt.subplot(6,1,6); plt.imshow( dEplot[nseg//2:,:].T, origin='lower', vmin=-dEmax, vmax=dEmax, cmap='bwr', extent=[ 0, len(lens), 0, 180.0 ] )
    #plt.colorbar()
    plt.xlabel("DOF")
    plt.ylabel("segment")
    plt.tight_layout()

def plot_Epanels_diff_separate(Emodels, Erefs, ref_dirs, save_prefix=None, bColorbar=True, Emin=-5.0, bKcal=False, bClose=False):
    """
    Plot each donor-acceptor pair in a separate figure with three panels (reference, model, difference).
    
    Args:
        Emodels: List of model energy arrays
        Erefs: List of reference energy arrays
        ref_dirs: List of reference directory names (donor-acceptor pairs)
        save_prefix: If provided, save each figure as '{save_prefix}_{ref_dir}.png'
        bColorbar: Whether to show colorbar
        Emin: Minimum energy for colormap scale
        bKcal: If True, convert energies to kcal/mol
        bClose: If True, close each figure after saving
    
    Returns:
        List of figure handles
    """
    E_units = 1.0
    if bKcal: 
        E_units = ev2kcal
    nmols = len(Erefs)
    if nmols != len(ref_dirs):
        print(f"Error: len(Erefs={nmols}) != len(ref_dirs={len(ref_dirs)})")
        return []
    
    figs = []
    for i in range(nmols):
        fig, axs = plt.subplots(3, 1, figsize=(4,6))  # 3 rows, 1 column
        
        Emodel = Emodels[i]
        Eref   = Erefs[i]
        
        # Plot Reference Energies
        im = axs[0].imshow(Eref.T*E_units, aspect='auto', origin='lower', vmin=Emin, vmax=-Emin, cmap='bwr')
        axs[0].set_ylabel('Reference Energies')
        if bColorbar: plt.colorbar(im, ax=axs[0])
        axs[0].set_title(f"{ref_dirs[i]}")
        
        # Plot Model Energies
        im = axs[1].imshow(Emodel.T*E_units, aspect='auto', origin='lower', vmin=Emin, vmax=-Emin, cmap='bwr')
        axs[1].set_ylabel('Model Energies')
        if bColorbar: plt.colorbar(im, ax=axs[1])
        
        # Plot Difference
        Ediff = Emodel - Eref
        im = axs[2].imshow(Ediff.T*E_units, aspect='auto', origin='lower', vmin=Emin, vmax=-Emin, cmap='bwr')
        axs[2].set_ylabel('Error')
        if bColorbar: plt.colorbar(im, ax=axs[2])
        
        plt.tight_layout()
        
        if save_prefix is not None:
            # Clean ref_dir name for filename (remove special characters)
            clean_name = ref_dirs[i].replace('/', '_').replace('\\', '_')
            filename = f"{save_prefix}_{clean_name}.png"
            fig.savefig(filename)
            print(f"Saved figure to: {filename}")
        
        if bClose:
            plt.close(fig)
        else:
            figs.append(fig)
    
    return figs


# ================== From opt_2D_new.py



# ============== Plotting Helper ==============

def plot_energy_2d_from_xyz(
    xyz_path,
    distances=None,
    angles=None,
    title=None,
    cmap='bwr',
    vmin=None,
    vmax=None,
    save_path=None,
):
    """Read an .xyz trajectory with comment lines carrying Etot/x0 and y|z angle tokens, build a 2D grid and plot it.

    - distances: list of distance values (floats) to define the x-axis grid order
    - angles: list of angle values (ints) to define the y-axis grid order
    Missing samples are filled with NaN.
    The function detects whether the file uses 'y' or 'z' angles and labels the axis accordingly.
    """
    import numpy as _np

    # Defaults from user's provided grids
    if distances is None:
        distances = [
            1.40, 1.45, 1.50, 1.55, 1.60, 1.65, 1.70, 1.75, 1.80, 1.85,
            1.90, 1.95, 2.00, 2.05, 2.10, 2.15, 2.20, 2.25, 2.30, 2.35,
            2.40, 2.45, 2.50, 2.60, 2.70, 2.80, 2.90, 3.00, 3.50, 4.00,
            4.50, 5.00, 6.00, 8.00, 10.00, 15.00, 20.00,
        ]
    if angles is None:
        angles = list(range(-90, 100, 10))  # -90..90 step 10, includes 0

    # Fast lookup maps
    dist_to_ix = {float(d): i for i, d in enumerate(distances)}
    ang_to_iy = {int(a): i for i, a in enumerate(angles)}

    # Prepare grid filled with NaN (shape: distances x angles)
    G = np.empty((len(distances), len(angles)), dtype=float)
    G[:] = np.nan

    direction_found = None  # 'y' or 'z'

    def _parse_comment(line: str):
        # Example: "# n0 3 Etot .4284 x0 01.80 z -90"
        toks = line.strip().split()
        vals = {"Etot": None, "x0": None, "angle": None, "axis": None}
        for i, t in enumerate(toks):
            if t == 'Etot' and i + 1 < len(toks):
                try:
                    vals["Etot"] = float(toks[i + 1])
                except Exception:
                    pass
            elif t == 'x0' and i + 1 < len(toks):
                try:
                    vals["x0"] = float(toks[i + 1])
                except Exception:
                    pass
            elif t in ('y', 'z', 'Y', 'Z') and i + 1 < len(toks):
                vals["axis"] = t.lower()
                try:
                    # angle tokens may be like '00' or '-90'
                    vals["angle"] = int(float(toks[i + 1]))
                except Exception:
                    pass
        return vals

    # Iterate frames: count line, comment line, then skip N atom lines
    with open(xyz_path, 'r', encoding='utf-8', errors='replace') as f:
        while True:
            count = f.readline()
            if not count:
                break
            count = count.strip()
            if not count:
                continue
            try:
                n = int(count)
            except Exception:
                # Not a proper count line; attempt to continue
                continue
            comment = f.readline()
            if not comment:
                break

            vals = _parse_comment(comment)
            if vals["axis"] is not None and direction_found is None:
                direction_found = vals["axis"]

            # Map to indices
            if vals["x0"] is not None and vals["angle"] is not None and vals["Etot"] is not None:
                # Map to (distance index, angle index)
                idist = dist_to_ix.get(round(vals["x0"], 2))
                if idist is None:
                    idist = dist_to_ix.get(vals["x0"])  # fallback exact float
                iang = ang_to_iy.get(int(vals["angle"]))
                if idist is not None and iang is not None:
                    G[idist, iang] = vals["Etot"]

            # Skip atom lines
            for _ in range(n):
                _ = f.readline()

    # Shift by asymptotic reference (last distance row minimum)
    ref_row = G[-1, :]
    ref = np.nanmin(ref_row[np.isfinite(ref_row)]) if np.any(np.isfinite(ref_row)) else 0.0
    GS = G - ref
    # Color scaling: vmin from global min (<=0), vmax = -vmin
    if vmin is None:
        try:
            vmin = float(np.nanmin(GS))
            if np.isfinite(vmin) and vmin > 0: vmin = 0.0
        except Exception:
            vmin = None
    if (vmin is not None) and (vmax is None):
        vmax = -vmin

    # Plot with angle on x-axis, distance on y-axis
    axlabel = f"angle ({direction_found})" if direction_found in ('y', 'z') else "angle"
    mesh = plt.imshow(GS, origin='lower', aspect='auto', cmap=cmap, vmin=vmin, vmax=vmax, interpolation='nearest')
    # label ticks by physical values while keeping uniform spacing
    xt = np.linspace(0, GS.shape[1]-1, min(6, GS.shape[1])).astype(int)
    yt = np.linspace(0, GS.shape[0]-1, min(6, GS.shape[0])).astype(int)
    plt.xticks(xt, [f"{angles[i]:.0f}" for i in xt])
    plt.yticks(yt, [f"{distances[i]:.2f}" for i in yt])
    plt.xlabel(axlabel + ' [deg]')
    plt.ylabel('distance x0 [Å]')
    if title is None:
        title = os.path.basename(xyz_path)
    plt.title(title)
    plt.colorbar(mesh, label='Etot shifted [a.u.]')
    if save_path is not None:
        print("Saving plot to:", save_path)
        plt.savefig(save_path, dpi=200, bbox_inches='tight')
    return GS


def parse_xyz_mapping(xyz_path, distances=None, angles=None):
    """Parse .xyz to build:
    - ref_grid (angles x distances) filled with Etot where present, NaN otherwise
    - sequence list of (iy, ix) per frame in file order for mapping model energies
    - detected axis label ('y' or 'z' or None)
    """
    if distances is None:
        distances = [
            1.40, 1.45, 1.50, 1.55, 1.60, 1.65, 1.70, 1.75, 1.80, 1.85,
            1.90, 1.95, 2.00, 2.05, 2.10, 2.15, 2.20, 2.25, 2.30, 2.35,
            2.40, 2.45, 2.50, 2.60, 2.70, 2.80, 2.90, 3.00, 3.50, 4.00,
            4.50, 5.00, 6.00, 8.00, 10.00, 15.00, 20.00,
        ]
    if angles is None:
        angles = list(range(-90, 100, 10))
    dist_to_ix = {float(d): i for i, d in enumerate(distances)}
    ang_to_iy = {int(a): i for i, a in enumerate(angles)}

    # grid shape: (distances, angles)
    Gref = np.empty((len(distances), len(angles)), dtype=float)
    Gref[:] = np.nan
    seq = []
    axis = None

    def _parse_comment(line: str):
        toks = line.strip().split()
        vals = {"Etot": None, "x0": None, "angle": None, "axis": None}
        for i, t in enumerate(toks):
            if t == 'Etot' and i + 1 < len(toks):
                try: vals["Etot"] = float(toks[i + 1])
                except Exception: pass
            elif t == 'x0' and i + 1 < len(toks):
                try: vals["x0"] = float(toks[i + 1])
                except Exception: pass
            elif t in ('y', 'z', 'Y', 'Z') and i + 1 < len(toks):
                vals["axis"] = t.lower()
                try: vals["angle"] = int(float(toks[i + 1]))
                except Exception: pass
        return vals

    with open(xyz_path, 'r', encoding='utf-8', errors='replace') as f:
        while True:
            count = f.readline()
            if not count: break
            count = count.strip()
            if not count: continue
            try:
                n = int(count)
            except Exception:
                continue
            comment = f.readline()
            if not comment: break
            vals = _parse_comment(comment)
            if vals["axis"] is not None and axis is None:
                axis = vals["axis"]
            idist = iang = None
            if vals["x0"] is not None:
                idist = dist_to_ix.get(round(vals["x0"], 2)) or dist_to_ix.get(vals["x0"]) 
            if vals["angle"] is not None:
                iang = ang_to_iy.get(int(vals["angle"]))
            if (iang is not None) and (idist is not None):
                seq.append((idist, iang))
                if vals["Etot"] is not None:
                    Gref[idist, iang] = vals["Etot"]
            for _ in range(n):
                _ = f.readline()
    return Gref, seq, axis, distances, angles


def compute_model_grid(xyz_path, seq, shape, do_fit=False, bAddEpairs=True, run_params=None, bOutXYZ=False):
    """Load the .xyz into FitREQ, optionally run fitting, compute model energies for each frame and map to grid.
    - seq: list of (iy, ix) in file order from parse_xyz_mapping
    - shape: (ny, nx)
    - run_params: dict with keys like nstep, ErrMax, dt, max_step, damping, bClamp
    """
    loadXYZ(xyz_path, bAddEpairs=bAddEpairs, bOutXYZ=bOutXYZ, bEvalOnlyCorrections=False)
    if do_fit:
        if run_params is None: run_params = dict(Fmax=1e-8, dt=0.5, damping=0.1, max_step=-1, bClamp=True, iparallel=0, ialg=1,)
        #run(**run_params)
        run( iparallel=0, ialg=1, nstep=run_params["nstep"], Fmax=1e-8, dt=0.5, damping=0.1,   max_step=-1,  bClamp=True )
    # Compute model energies for all frames
    _, Es, _ = getEs(bOmp=False, bEs=True, bFs=False)
    Gm = np.empty(shape, dtype=float); Gm[:] = np.nan
    nmap = min(len(Es), len(seq))
    for i in range(nmap):
        idist, iang = seq[i]
        Gm[idist, iang] = Es[i]
    return Gm


def shift_grid(G):
    import numpy as _np
    if G.size == 0:
        return G, 0.0, 0.0
    last = G[-1, :] if (G.ndim == 2 and G.shape[0] > 0) else np.array([0.0])
    ref = float(np.nanmin(last[np.isfinite(last)])) if np.any(np.isfinite(last)) else 0.0
    GS = G - ref
    mloc = float(np.nanmin(GS)) if np.any(np.isfinite(GS)) else 0.0
    if np.isfinite(mloc) and mloc > 0: mloc = 0.0
    return GS, ref, mloc


def extract_min_curves(angles, distances, G, rmax=None):
    import numpy as _np
    nA = len(angles)
    rmin = np.full(nA, np.nan)
    emin = np.full(nA, np.nan)
    for j in range(nA):
        col = G[:, j]
        if np.any(np.isfinite(col)):
            idx = int(np.nanargmin(col))
            emin[j] = col[idx]
            rmin[j] = distances[idx]
            if (rmax is not None) and np.isfinite(rmin[j]) and (rmin[j] > rmax):
                rmin[j] = np.nan; emin[j] = np.nan
    return rmin, emin


def plot_compare(Gref, Gmodel, angles, distances, title, save_prefix=None, vmin=None, vmax=None, line=False, kcal=False):
    import matplotlib.pyplot as plt

    # Shift each by its own baseline and determine symmetric limits from reference
    GRS, refR, mlocR = shift_grid(Gref)
    GMS = None
    if Gmodel is not None: GMS, refM, mlocM = shift_grid(Gmodel)
        
    if line:
        # Reuse plot_min_lines_pair by building panel-shaped inputs (angles x distances)
        Epanel_ref = GRS.T
        Epanel_mod = GMS.T if GMS is not None else None
        # Xpanel: replicate distances across all angles
        Xpanel = np.tile(np.asarray(distances, dtype=float), (len(angles), 1))
        fig2 = plot_min_lines_pair(Epanel_ref, Epanel_mod, Xpanel, angles, title=None, save_path=None, to_kcal=kcal)
        if save_prefix:
            p2 = save_prefix.replace('.png','') + '_lines.png'
            print("Saving plot to:", p2)
            fig2.savefig(p2, dpi=150, bbox_inches='tight')

    if kcal:
        if GRS is not None: GRS *= ev2kcal
        if GMS is not None: GMS *= ev2kcal

    vminR = float(np.nanmin(GRS)) if np.any(np.isfinite(GRS)) else None
    if vmin is None:
        vmin = vminR if (vminR is not None) else None
        if (vmin is not None) and (vmin > 0): vmin = 0.0
    if (vmin is not None) and (vmax is None):
        vmax = -vmin

    # 3 rows: Ref, Model, Diff (if model present); else only Ref
    rows = 3 if (GMS is not None) else 1
    fig, axs = plt.subplots(rows, 1, figsize=(4, 2.2*rows), constrained_layout=True)
    if rows == 1:
        axs = [axs]
    # Ref
    im0 = axs[0].imshow(GRS, origin='lower', aspect='auto', cmap='bwr', vmin=vmin, vmax=vmax)
    xt = np.linspace(0, GRS.shape[1]-1, min(6, GRS.shape[1])).astype(int)
    yt = np.linspace(0, GRS.shape[0]-1, min(6, GRS.shape[0])).astype(int)
    axs[0].set_xticks(xt); axs[0].set_yticks(yt)
    x_ticks=[f"{angles[i]:.0f}" for i in xt]
    y_ticks=[f"{distances[i]:.2f}" for i in yt]
    axs[0].set_xticklabels(x_ticks); 
    axs[0].set_yticklabels(y_ticks)
    axs[0].set_ylabel('Distance (Å)')
    axs[0].set_title(title)
    plt.colorbar(im0, ax=axs[0])
    if GMS is not None:
        im1 = axs[1].imshow(GMS, origin='lower', aspect='auto', cmap='bwr', vmin=vmin, vmax=vmax)
        axs[1].set_xticks(xt); axs[1].set_yticks(yt)
        axs[1].set_xticklabels(x_ticks); 
        axs[1].set_yticklabels(y_ticks)
        axs[1].set_ylabel('Distance (Å)')
        plt.colorbar(im1, ax=axs[1])
        # Diff
        D = GMS - GRS
        dmax = max(-np.nanmin(D), np.nanmax(D)) if np.any(np.isfinite(D)) else 1.0
        im2 = axs[2].imshow(D, origin='lower', aspect='auto', cmap='bwr', vmin=-dmax, vmax=dmax)
        axs[2].set_xticks(xt); axs[2].set_yticks(yt)
        axs[2].set_xticklabels(x_ticks); 
        axs[2].set_yticklabels(y_ticks)
        axs[2].set_xlabel('Angle (deg)')
        axs[2].set_ylabel('Distance (Å)')
        plt.colorbar(im2, ax=axs[2])
    else:
        axs[0].set_xlabel('Angle (deg)')

    if save_prefix:
        fname = f"{save_prefix}.png" if not save_prefix.endswith('.png') else save_prefix
        print("Saving plot to:", fname)
        fig.savefig(fname, dpi=150, bbox_inches='tight')

# ===== Reusable helpers for Rmin/Emin from panel-shaped data (angles x distances)

def _distances_from_Xpanel(Xpanel):
    """Given Xpanel shaped (ny_angles, nx_distances) with distances per sample,
    build a single distances vector of length nx by taking the first finite value in each column.
    Fallback to NaN if none is finite.
    """
    ny, nx = Xpanel.shape
    dists = np.full(nx, np.nan)
    for j in range(nx):
        col = Xpanel[:, j]
        mask = np.isfinite(col)
        if np.any(mask):
            dists[j] = col[np.where(mask)[0][0]]
    return dists


def compute_min_lines_from_panel(Epanel, Xpanel, angles, rmax=None, do_shift=True):
    """Compute Rmin(angle) and Emin(angle) from panel-shaped data.

    Inputs:
      - Epanel: 2D array (ny_angles, nx_distances) of energies
      - Xpanel: 2D array (ny_angles, nx_distances) of distances x0 corresponding to columns
      - angles: list/array of angle values of length ny_angles
      - rmax: optional cutoff; if provided, entries with rmin > rmax are set to NaN
      - do_shift: if True, shift grid by asymptotic baseline before extracting

    Returns:
      rmin, emin arrays of length ny_angles
    """
    # Convert to grid shaped (distances, angles) expected by shift_grid/extract_min_curves
    G = np.asarray(Epanel, dtype=float).T  # (nx, ny)
    distances = _distances_from_Xpanel(np.asarray(Xpanel, dtype=float))
    if do_shift:
        G, _, _ = shift_grid(G)
    rmin, emin = extract_min_curves(angles=np.asarray(angles), distances=distances, G=G, rmax=rmax)
    return rmin, emin


def plot_min_lines_pair(Epanel_ref, Epanel_mod, Xpanel, angles, title=None, save_path=None, to_kcal=False, ms=2, lw=0.5):
    """Plot Rmin(angle) and Emin(angle) lines for a ref/model pair using panel-shaped inputs.

    - Epanel_ref/Epanel_mod: (ny_angles, nx_distances)
    - Xpanel: distances panel (ny_angles, nx_distances)
    - angles: sequence of ny angle values
    - to_kcal: if True, convert energy lines to kcal/mol using ev2kcal
    """
    rR, eR = compute_min_lines_from_panel(Epanel_ref, Xpanel, angles)
    rM = eM = None
    if Epanel_mod is not None:
        rM, eM = compute_min_lines_from_panel(Epanel_mod, Xpanel, angles)
    import matplotlib.pyplot as plt
    fig, (axR, axE) = plt.subplots(1, 2, figsize=(6.5, 2.8), constrained_layout=True)
    # ------ Distance
    axR.plot(angles, rR, '.-k', lw=lw,ms=ms, label='Ref')
    if rM is not None: 
        axR.plot(angles, rM, '.-r', lw=lw,ms=ms, label='Model')
    axR.set_title('r_min(angle)');  
    axR.set_xlabel('Angle (deg)');   
    axR.set_ylabel('Distance x0 [Å]');
    axR.set_ylim(None,3.0)
    axR.grid(alpha=0.3, linestyle='--'); 
    axR.legend(fontsize=8);
    # ------ Energy
    Er = eR * (ev2kcal if to_kcal else 1.0)
    axE.plot(angles, Er, '.-k', lw=lw,ms=ms, label='Ref')
    if eM is not None:
        Em = eM * (ev2kcal if to_kcal else 1.0)
        axE.plot(angles, Em, '.-r', lw=lw,ms=ms, label='Model')
    axE.set_title('E_min(angle)')
    axE.set_xlabel('Angle (deg)')
    axE.set_ylabel('Energy [{}]'.format('kcal/mol' if to_kcal else 'eV'))
    axE.grid(alpha=0.3, linestyle='--')
    axE.legend(fontsize=8)
    if title:     fig.suptitle(title, fontsize=10)
    if save_path: 
        print("Saving plot to:", save_path)
        fig.savefig(save_path, bbox_inches='tight')
    return fig


def plot_trj_dofs(trj_DOFs, DOFnames=None, lss=None, clrs=None, title=None, save_path=None):
    """Plot trajectories of DOFs over optimization iterations.

    Args:
        trj_DOFs: ndarray (niter, nDOFs)
        DOFnames: optional list of names per DOF; if None, uses generic labels
        lss: list of line styles to cycle
        clrs: list of colors to cycle
        title: optional title for the plot

    Returns:
        matplotlib.figure.Figure or None
    """
    import matplotlib.pyplot as plt
    if trj_DOFs is None: return None
    niter, nDOFs_ = trj_DOFs.shape
    if DOFnames is None: DOFnames = [f"DOF {i}" for i in range(nDOFs_)]
    if lss      is None: lss = ['-','--',":",'-','--',":",'-','-']
    if clrs     is None: clrs = ['b','b','b','r','r','r','c','m']
    fig = plt.figure()
    for i in range(min(nDOFs_, len(DOFnames))):
        ls = lss[i % len(lss)]
        c  = clrs[i % len(clrs)]
        plt.plot(trj_DOFs[:, i], ls, c=c, lw=1.0, ms=2.0, label=DOFnames[i])
    if title:
        plt.title(title)
    plt.legend(fontsize=8)
    plt.grid(alpha=0.2)
    if save_path: 
        print("Saving plot to:", save_path)
        fig.savefig(save_path, bbox_inches='tight')
    return fig
