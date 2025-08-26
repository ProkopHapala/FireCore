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

# void setGlobalParams( double kMorse, double Lepairs ){
lib.setGlobalParams.argtypes  = [c_double, c_double]
lib.setGlobalParams.restype   =  None
def setGlobalParams(kMorse=1.6, Lepairs=0.5):
    return lib.setGlobalParams(kMorse, Lepairs)

# void setup( int imodel, int EvalJ, int WriteJ, int CheckRepulsion, int Regularize, int AddRegError, int Epairs, int BroadcastFDOFs, int UdateDOFbounds, int EvalOnlyCorrections){
lib.setup.argtypes  = [c_int,c_int, c_int, c_int, c_int, c_int, c_int, c_int, c_int, c_int]
lib.setup.restype   =  None    
def setup(imodel=1, EvalJ=0, WriteJ=0, CheckRepulsion=0, Regularize=0, AddRegError=0, Epairs=0, BroadcastFDOFs=0, UdateDOFbounds=0, EvalOnlyCorrections=0):
    return lib.setup( imodel,EvalJ, WriteJ, CheckRepulsion, Regularize, AddRegError, Epairs, BroadcastFDOFs, UdateDOFbounds, EvalOnlyCorrections)

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

# void loadTypes_new( const char* fname_ElemTypes, const char* fname_AtomTypes ){
lib.loadTypes.argtypes  = [c_char_p, c_char_p]
lib.loadTypes.restype   =  None
def loadTypes(fEtypes="data/ElementTypes.dat", fAtypes="data/AtomTypes.dat"):
    return lib.loadTypes(cstr(fEtypes), cstr(fAtypes))

#int loadDOFSelection( const char* fname ){
lib.loadDOFSelection.argtypes  = [c_char_p]
lib.loadDOFSelection.restype   =  c_int
def loadDOFSelection(fname="dofSelection.dat"):
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

# int loadXYZ( const char* fname, bool bAddEpairs, bool bOutXYZ, bool bEvalOnlyCorrections ){ 
lib.loadXYZ.argtypes  = [c_char_p, c_bool, c_bool, c_bool]
lib.loadXYZ.restype   =  c_int
def loadXYZ(fname, bAddEpairs=False, bOutXYZ=False, bEvalOnlyCorrections=False ):
    global nbatch
    nbatch = lib.loadXYZ(cstr(fname), bAddEpairs, bOutXYZ, bEvalOnlyCorrections)
    return nbatch

# void exportSampleToXYZ(int i, const char* fname)
lib.exportSampleToXYZ.argtypes = [c_int, c_char_p]
lib.exportSampleToXYZ.restype  = None
def exportSampleToXYZ(i, fname):
    return lib.exportSampleToXYZ(i, cstr(fname))

# void exportAllSamplesToXYZ(const char* fname)
lib.exportAllSamplesToXYZ.argtypes = [c_char_p]
lib.exportAllSamplesToXYZ.restype  = None
def exportAllSamplesToXYZ(fname):
    return lib.exportAllSamplesToXYZ(cstr(fname))

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

def loadDOFnames( fname="dofSelection.dat", comps="REQH", return_specs=False ):
    """Load DOF names from a selection file.

    If return_specs=True, also parse and return specification dicts per DOF with keys:
      name, comp, comp_char, min, max, xlo, xhi, Klo, Khi, K0, xstart

    Notes:
    - Lines starting with '#' or empty after stripping inline comments are skipped
    - Inline comments after '#' are ignored
    - Missing numeric fields are returned as None
    """
    names = []
    specs = []
    with open(fname, 'r') as f:
        for line in f:
            # Strip inline comments and whitespace
            line_clean = line.split('#', 1)[0].strip()
            if not line_clean:
                continue
            parts = line_clean.split()
            if len(parts) < 2:
                continue
            typename = parts[0]
            icomp = int(parts[1])
            comp_char = comps[icomp] if (0 <= icomp < len(comps)) else str(icomp)
            names.append(f"{typename}.{comp_char}")
            if return_specs:
                # Safely parse remaining numeric fields if present
                def fget(i):
                    return float(parts[i]) if len(parts) > i else None
                spec = {
                    'name': typename,
                    'comp': icomp,
                    'comp_char': comp_char,
                    'min': fget(2),
                    'max': fget(3),
                    'xlo': fget(4),
                    'xhi': fget(5),
                    'Klo': fget(6),
                    'Khi': fget(7),
                    'K0':  fget(8),
                    'xstart': fget(9),
                }
                specs.append(spec)
    if return_specs:
        return names, specs
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

def plotDOFscan_one(iDOF, xs=None, DOFname=None, bEs=True, bFs=False, bEvalSamples=True, axE=None, axF=None, color=None, verb=1, npts=100, xrange=(0.0,1.0), data=None ):
    """Plot a single DOF scan.

    - If axE/axF are provided, plot into them; else create a new 2x1 figure.
    - If xs is None and dof_spec has 'min'/'max', use that range with npts samples.
    - If data is provided, it must be either a dict with keys 'xs','Es','Fs' or a tuple (xs,Es,Fs);
      in that case, values are plotted directly without calling scanParam().
    - Returns (xs, Es, Fs, Fs_num, xs_num)
    """
    if plt is None:
        raise RuntimeError("Matplotlib handle 'plt' not set. Assign fit.plt = plt in the caller.")
    make_fig = (axE is None) or (axF is None)
    if make_fig:
        fig = plt.figure(figsize=(8,10.0))
        axE = plt.subplot(2,1,1)
        axF = plt.subplot(2,1,2)
    # Prepare data (either provided or computed)
    if data is not None:
        # Accept dict with keys or a tuple
        if isinstance(data, dict):
            xs = data.get('xs', xs)
            Es = data.get('Es', None)
            Fs = data.get('Fs', None)
        else:
            try:
                xs, Es, Fs = data
            except Exception as e:
                raise ValueError("data must be dict with keys 'xs','Es','Fs' or tuple (xs,Es,Fs)") from e
        if (xs is None) or (Es is None) or (Fs is None):
            raise ValueError("plotDOFscan_one(): incomplete data provided; require xs, Es, Fs")
    else:
        # Determine xs if not provided, then compute by CPU scan
        if xs is None: xs = np.linspace(xrange[0]+1e-6, xrange[1]-1e-6, npts)
        y_backup = DOFs[iDOF]
        Es, Fs = scanParam(iDOF, xs, bEvalSamples=bEvalSamples)
        DOFs[iDOF] = y_backup
    # Choose color
    label = ("E "+DOFname) if DOFname else f"E DOF {iDOF}"
    if bEs:
        axE.plot(xs, Es, '-', lw=0.5, color=color, label=label)
    xs_num = None
    Fs_num = None
    if bFs:
        Fs_num, xs_num = numDeriv(xs, Es)
        Fs_num = -Fs_num  # F = -dE/dx
        labF = ("F "+DOFname) if DOFname else f"F DOF {iDOF}"
        axF.plot(xs,     Fs,     '-', lw=0.5, color=color, label=labF+"_ana")
        axF.plot(xs_num, Fs_num, ':', lw=1.5, color=color, label=labF+"_num")
        if verb > 0:
            Eminmax = (np.min(Es), np.max(Es))
            Fminmax = (np.min(Fs), np.max(Fs))
            Fnumminmax = (np.min(Fs_num), np.max(Fs_num)) if Fs_num is not None else (None, None)
            Ferr = Fs[1:-1] + (-Fs_num) if Fs_num is not None else None
            Ferrmax = (np.min(Ferr), np.max(Ferr)) if Ferr is not None else (None, None)
            relFerrmax = np.max(np.abs(Ferr)/(np.max(Fs)+1e-8))
            nm = DOFname if DOFname else f"{iDOF}"
            print(f"plotDOFscan_one(iDOF={iDOF} : {nm})  relFerrmax: {relFerrmax}  (F_ana-F_num)[min,max]: {Ferrmax}  |  E[min,max]: {Eminmax} F_ana[min,max]: {Fminmax} F_num[min,max]: {Fnumminmax}  ")
    if make_fig:
        axE.legend(); axE.set_xlabel("DOF value"); axE.set_ylabel("E [kcal/mol]");   axE.grid( alpha=0.2)
        axF.legend(); axF.set_xlabel("DOF value"); axF.set_ylabel("F [kcal/mol/A]"); axF.grid( alpha=0.2)
    return xs, Es, Fs, Fs_num, xs_num

def plotDOFscans( iDOFs, xs, DOFnames, bEs=True, bFs=False,  title="plotDOFscans", bEvalSamples=True, verb=1 ):
    plt.figure(figsize=(8,10.0))
    axE = plt.subplot(2,1,1)
    axF = plt.subplot(2,1,2)
    color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
    ncol = len(color_cycle)
    for iDOF in iDOFs:
        c = color_cycle[iDOF % ncol]
        plotDOFscan_one(iDOF, xs=xs, DOFname=DOFnames[iDOF], bEs=bEs, bFs=bFs, bEvalSamples=bEvalSamples, axE=axE, axF=axF, color=c, verb=verb)
    axE.legend(); axE.set_xlabel("DOF value"); axE.set_ylabel("E [kcal/mol]"); axE.grid( alpha=0.2)
    axF.set_xlabel("DOF value"); axF.set_ylabel("F [kcal/mol/A]"); axF.grid( alpha=0.2)
    plt.suptitle(title)
    # plt.show()  # Let caller decide

def checkDOFderiv( iDOF, x0=0.5, d=0.001, bEvalSamples=True, dof_names=None ):
    xs = np.array([x0-d,x0,x0+d])
    Es, Fs = scanParam( iDOF, xs, bEvalSamples=bEvalSamples )   # do 1D scan
    F_num, xs_num = numDeriv(xs, Es)
    F_num = -F_num  # F = -dE/dx
    if dof_names is not None:
        Eminmax=(np.min(Es),np.max(Es))
        Fminmax=(np.min(Fs),np.max(Fs))
        Fnumminmax=(np.min(F_num),np.max(F_num))
        Ferr = Fs[1:-1]-F_num  # analytic minus numeric
        Ferrmax=(np.min(Ferr),np.max(Ferr))
        print( f"checkDOFderiv(iDOF={iDOF} : {dof_names[iDOF]}) x={x0} d={d} E[min,max]={Eminmax} F_ana[min,max]={Fminmax} F_num[min,max]={Fnumminmax} (F_ana-F_num)[min,max]={Ferrmax}" )
    return Es, Fs, F_num

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
