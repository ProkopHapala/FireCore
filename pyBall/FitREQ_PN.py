import numpy as np
from   ctypes import c_int, c_double, c_bool, c_float, c_char_p, c_bool, c_void_p, c_char_p
import ctypes
import os
from . import cpp_utils_ as cpp_utils

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
lib = cpp_utils.loadLib('FitREQ_PN_lib', recompile=False)
array1ui = np.ctypeslib.ndpointer(dtype=np.uint32, ndim=1, flags='CONTIGUOUS')
array1i  = np.ctypeslib.ndpointer(dtype=np.int32,  ndim=1, flags='CONTIGUOUS')
array2i  = np.ctypeslib.ndpointer(dtype=np.int32,  ndim=2, flags='CONTIGUOUS')
array1d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS')
array2d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=2, flags='CONTIGUOUS')
array3d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=3, flags='CONTIGUOUS')

# ====================================
# ========= Globals
# ====================================
ev2kcal = 23.060547831
bWeightsSet = False
#isInitialized = False
#nfound = -1

# ====================================
# ========= C functions
# ====================================

def cstr( s ):
    if s is None: return None
    return s.encode('utf8')

#void setVerbosity( int verbosity_, int idebug_, int PrintDOFs, int PrintfDOFs, int PrintBeforReg, int PrintAfterReg, int PrintOverRepulsive ){
lib.setVerbosity.argtypes  = [c_int, c_int, c_int, c_int, c_int, c_int, c_int]
lib.setVerbosity.restype   =  None
def setVerbosity(verbosity=0, idebug=0, PrintDOFs=0, PrintfDOFs=0, PrintBeforReg=0, PrintAfterReg=0, PrintOverRepulsive=0):
    return lib.setVerbosity(verbosity, idebug, PrintDOFs, PrintfDOFs, PrintBeforReg, PrintAfterReg, PrintOverRepulsive)

#void setModel( int ivdW, int iCoul, int iHbond, int Epairs, int iEpairs, double kMorse, double Lepairs, bool bPN ){
lib.setModel.argtypes  = [c_int, c_int, c_int, c_int, c_int, c_double, c_double, c_bool]
lib.setModel.restype   =  None
def setModel(ivdW=4, iCoul=1, iHbond=0, Epairs=0, iEpairs=0, kMorse=1.6, Lepairs=0.5, bPN=True):
    return lib.setModel(ivdW, iCoul, iHbond, Epairs, iEpairs, kMorse, Lepairs, bPN)

#void loadTypes( const char* fname_ElemTypes, const char* fname_AtomTypes ){
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

#int loadXYZ( const char* fname, bool bAddEpairs, bool bOutXYZ, bool bEvalOnlyCorrections, bool bAppend ){ 
lib.loadXYZ.argtypes  = [c_char_p, c_bool, c_bool, c_bool, c_char_p, c_bool, c_bool]
lib.loadXYZ.restype   =  c_int
def loadXYZ(fname, bAddEpairs=False, bOutXYZ=False, bSaveJustElementXYZ=False, OutXYZ_fname='out_epairs.xyz', bEvalOnlyCorrections=False, bAppend=False ):
    global nbatch
    nbatch = lib.loadXYZ(cstr(fname), bAddEpairs, bOutXYZ, bSaveJustElementXYZ, cstr(OutXYZ_fname), bEvalOnlyCorrections, bAppend)
    return nbatch

lib.setPenalty.argtypes  = [c_int, c_int, c_int, c_int, c_int, c_double, c_double]
lib.setPenalty.restype   =  None
def setPenalty(Clamp=1, Regularize=0, AddRegError=0, RegCountWeight=0, SoftClamp=0, softClamp_start=0.17, softClamp_max=0.26):
    lib.setPenalty(Clamp, Regularize, AddRegError, RegCountWeight, SoftClamp, softClamp_start, softClamp_max)

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
def setTrjBuffs( niter, trj_E=None, trj_F=None, trj_DOFs=None, trj_fDOFs=None, nDOFs_=None, bE=True, bF=True, bDOFs=True, bfDOFs=True ):
    if nDOFs_     is None: nDOFs_ = nDOFs
    print("setTrjBuffs(): niter=%i nDOFs=%i bE=%i bF=%i bDOFs=%i bfDOFs=%i" % (niter, nDOFs_, bE, bF, bDOFs, bfDOFs))
    #if (trj_E     is None) and bE     : trj_E     = np.zeros( niter )
    #if (trj_F     is None) and bF     : trj_F     = np.zeros( niter )
    #if (trj_DOFs  is None) and bDOFs  : trj_DOFs  = np.zeros( (niter, nDOFs_) )
    #if (trj_fDOFs is None) and bfDOFs : trj_fDOFs = np.zeros( (niter, nDOFs_) )
    if (trj_E     is None) and bE     : trj_E     = np.full( niter, np.nan )
    if (trj_F     is None) and bF     : trj_F     = np.full( niter, np.nan )
    if (trj_DOFs  is None) and bDOFs  : trj_DOFs  = np.full( (niter, nDOFs_), np.nan )
    if (trj_fDOFs is None) and bfDOFs : trj_fDOFs = np.full( (niter, nDOFs_), np.nan )
    lib.setTrjBuffs(_np_as(trj_E,c_double_p), _np_as(trj_F,c_double_p), _np_as(trj_DOFs,c_double_p), _np_as(trj_fDOFs,c_double_p))
    return trj_E, trj_F, trj_DOFs, trj_fDOFs

#double run_PN( int ialg, int iparallel, int nstep, double Fmax, double dt, double max_step, double damping ){
lib.run_PN.argtypes  = [c_int, c_int, c_int, c_double, c_double, c_double, c_double]
lib.run_PN.restype   =  c_double
def run_PN(ialg=2, iparallel=0, nstep=1000, Fmax=1e-8, dt=0.01, max_step=0.05, damping=0.01):
    return lib.run_PN( ialg, iparallel, nstep, Fmax, dt, max_step, damping )

#double getEs_components( double* Es, double* Es_Coul, double* Es_vdW, double* Es_Epairs, double* Es_Hbond){
lib.getEs_components.argtypes  = [ c_double_p, c_double_p, c_double_p, c_double_p, c_double_p ]
lib.getEs_components.restype   =  c_double
def getEs_components( Es=None, Es_Coul=None, Es_vdW=None, Es_Hcorr=None, Es_Epairs=None ):
    Es       = np.zeros( nbatch )
    Es_Coul  = np.zeros( nbatch )
    Es_vdW   = np.zeros( nbatch )
    Es_Hcorr = np.zeros( nbatch )
    Es_Epairs= np.zeros( nbatch )
    lib.getEs_components( _np_as(Es,c_double_p), _np_as(Es_Coul,c_double_p), _np_as(Es_vdW,c_double_p), _np_as(Es_Hcorr,c_double_p), _np_as(Es_Epairs,c_double_p) )
    return Es, Es_Coul, Es_vdW, Es_Hcorr, Es_Epairs

#double getError( int iparallel ){
lib.getError.argtypes  = [ c_int ]
lib.getError.restype   =  c_double
def getError( iparallel=0 ):
    return lib.getError( iparallel )

# void scanParam( int iDOF, int n, double* xs,  double* Es, double* Fs, bool bRegularize, bool bEvalSamples ){
lib.scanParam.argtypes  = [c_int, c_int, c_double_p, c_double_p, c_double_p, c_bool]
lib.scanParam.restype   = None
def scanParam( iDOF, xs, Es=None, Fs=None, bEvalSamples=True ):
    n = len(xs)
    if Es is None: Es = np.zeros( n )
    if Fs is None: Fs = np.zeros( n )
    lib.scanParam(iDOF, n, _np_as(xs,c_double_p), _np_as(Es,c_double_p), _np_as(Fs,c_double_p), bEvalSamples )
    return Es,Fs

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


# =============== Buffers
    
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

#void init_buffers(){
lib.init_buffers.argtypes  = []
lib.init_buffers.restype   =  None
def init_buffers():
    #print( "init_buffers()" )
    return lib.init_buffers()

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
    ptr = lib.getBuff(name) # from libUtils.h (I guess)
    if(ptr):
        return np.ctypeslib.as_array(ptr, shape=sh)
    else:
        return None

################## Python ###############

def read_xyz_data(fname="input_all.xyz"):
    """Read XYZ file and extract Etot and x0 values from comment lines"""
    #print(f"DEBUG: read_xyz_data() reading file: {fname}", flush=True)
    Erefs = []
    x0s = []
    line_count = 0
    matched_count = 0
    sample_lines = []
    
    with open(fname, 'r') as f:
        while True:
            line = f.readline()
            line_count += 1
            if not line: break
            
            # Keep first 5 comment lines as samples
            if line.startswith('#') and len(sample_lines) < 5:
                sample_lines.append(line.strip())
            
            if line.startswith('# n0'):
                matched_count += 1
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
    
    #print(f"DEBUG: Total lines read: {line_count}", flush=True)
    #print(f"DEBUG: Lines matched with '# n0': {matched_count}", flush=True)
    #print(f"DEBUG: Sample comment lines:", flush=True)
    #for sl in sample_lines:
    #    print(f"  {sl}", flush=True)
    #print(f"DEBUG: Erefs array length: {len(Erefs)}, x0s array length: {len(x0s)}", flush=True)
    #if len(Erefs) > 0:
    #    print(f"DEBUG: First Eref: {Erefs[0]}, First x0: {x0s[0]}", flush=True)
    
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

#def save_trj_dofs(trj_DOFs, DOFnames=None, folder=None):
def save_trj_dofs(trj_DOFs, DOFnames=None, folder=None, n_steps=None):
    """
    Save each DOF trajectory as a separate .dat file in a folder.

    Each file will contain two columns:
        iteration_index   DOF_value

    Example output file (for DOF_0.dat):
        # iter   DOF_0
        0  0.12345
        1  0.12789
        2  0.13001
    """
    if trj_DOFs is None or folder is None:
        return

    niter, nDOFs_ = trj_DOFs.shape
    
    if n_steps is None:
        nans = np.isnan(trj_DOFs[:,0])
        if np.any(nans):
            n_steps = np.argmax(nans)
        else:
            n_steps = niter
    n_save = min(niter, int(n_steps))

    if DOFnames is None:
        DOFnames = [f"DOF_{i}" for i in range(nDOFs_)]

    # Ensure the folder exists
    try:
        os.makedirs(folder, exist_ok=True)
    except Exception:
        pass

    # Iteration indices for first column
    #iters = np.arange(niter)
    iters = np.arange(n_save)

    # Write each DOF file
    for i in range(nDOFs_):
        filename = os.path.join(folder, f"{DOFnames[i]}.dat")
        #data = np.column_stack((iters, trj_DOFs[:, i]))
        data = np.column_stack((iters, trj_DOFs[:n_save, i]))
        header = f"# iter  {DOFnames[i]}"
        np.savetxt(filename, data, header=header, comments="", fmt="%.12g")
        #print(f"Saved {DOFnames[i]}.dat trajectory file with iteration indices to: {folder}")
    print(f"Saved {nDOFs_} trajectory files with iteration indices to: {folder}")

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
                    #print("Etot", vals["Etot"], " at x0=", vals["x0"], " angle=", vals["angle"])
            for _ in range(n):
                _ = f.readline()
    return Gref, seq, axis, distances, angles

def save_data(Gref, Gmodel, angles, distances, save_data_prefix=None, save_fmt="both", kcal=False, line=False):

    # Shift each by its own baseline and determine symmetric limits from reference
    GRS, refR, mlocR = shift_grid(Gref)
    GMS = None
    if Gmodel is not None: 
        GMS, refM, mlocM = shift_grid(Gmodel)
        print( "save_data() Gmodel min,max", Gmodel.min(), Gmodel.max() )
        print( "save_data() GMS    min,max", GMS.min(), GMS.max() )

    if kcal:
        if GRS is not None: GRS *= ev2kcal
        if GMS is not None: GMS *= ev2kcal

    if save_data_prefix is None:
        print("No save_data_prefix specified, nothing written.")
        return

    # Ensure output folder exists
    folder = os.path.dirname(save_data_prefix)
    if folder:
        os.makedirs(folder, exist_ok=True)

    if GRS is not None:
        if save_fmt in ("both", "npz"):     save_grid_npz(angles, distances, GRS, save_data_prefix + "__ref.npz")
        if save_fmt in ("both", "gnuplot"): save_grid_gnuplot(angles, distances, GRS, save_data_prefix + "__ref.dat")
    if GMS is not None:
        if save_fmt in ("both", "npz"):     save_grid_npz(angles, distances, GMS, save_data_prefix + "__model.npz")
        if save_fmt in ("both", "gnuplot"): save_grid_gnuplot(angles, distances, GMS, save_data_prefix + "__model.dat")
        D = GMS - GRS
        if save_fmt in ("both", "npz"):     save_grid_npz(angles, distances, D, save_data_prefix + "__diff.npz")
        if save_fmt in ("both", "gnuplot"): save_grid_gnuplot(angles, distances, D, save_data_prefix + "__diff.dat")

    if line:
        # Reuse plot_min_lines_pair by building panel-shaped inputs (angles x distances)
        Epanel_ref = GRS.T
        if GMS is not None:
            Epanel_mod = GMS.T 
        # Xpanel: replicate distances across all angles
        Xpanel = np.tile(np.asarray(distances, dtype=float), (len(angles), 1))
        rR, eR = compute_min_lines_from_panel(Epanel_ref, Xpanel, angles)
        if GMS is not None:
            rM, eM = compute_min_lines_from_panel(Epanel_mod, Xpanel, angles)
        if save_fmt in ("both","npz"):         save_min_lines_npz(angles, rR, eR, save_data_prefix + "__ref_lines.npz")
        if save_fmt in ("both","gnuplot"):     save_min_lines_gnuplot(angles, rR, eR, save_data_prefix + "__ref_lines.dat")
        if GMS is not None:
            if save_fmt in ("both","npz"):     save_min_lines_npz(angles, rM, eM, save_data_prefix + "__model_lines.npz")
            if save_fmt in ("both","gnuplot"): save_min_lines_gnuplot(angles, rM, eM, save_data_prefix + "__model_lines.dat")

def shift_grid(G):
    import numpy as _np
    if G.size == 0:
        return G, 0.0, 0.0
    last = G[-1, :] if (G.ndim == 2 and G.shape[0] > 0) else np.array([0.0])
    ref = float(np.nanmin(last[np.isfinite(last)])) if np.any(np.isfinite(last)) else 0.0
    #GS = G - ref
    GS = G
    mloc = float(np.nanmin(GS)) if np.any(np.isfinite(GS)) else 0.0
    if np.isfinite(mloc) and mloc > 0: mloc = 0.0
    return GS, ref, mloc

def save_grid_npz(angles, distances, grid, filepath):
    """Save 2D map to NPZ with keys a, d, g."""
    folder = os.path.dirname(filepath)
    if folder:
        try: os.makedirs(folder, exist_ok=True)
        except Exception: pass
    print("save_grid_npz(): saving to", filepath)
    np.savez_compressed(filepath, a=np.asarray(angles), d=np.asarray(distances), g=np.asarray(grid))

def save_grid_gnuplot(angles, distances, grid, filepath):
    """Save 2D map as gnuplot-friendly triplets: angle distance energy per row. NaNs -> 'nan'."""
    folder = os.path.dirname(filepath)
    if folder:
        try: os.makedirs(folder, exist_ok=True)
        except Exception: pass
    print("save_grid_gnuplot(): saving to", filepath)
    with open(filepath, 'w') as f:
        f.write("# angle  distance  energy\n")
        for j, a in enumerate(angles):
            for i, d in enumerate(distances):
                e = grid[i, j]
                if np.isnan(e):
                    f.write(f"{a:<12.6g} {d:<12.6g} nan\n")
                else:
                    f.write(f"{a:<12.6g} {d:<12.6g} {e:<12.12g}\n")
            f.write(f"\n")

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

def extract_min_curves(angles, distances, G, rmax=None):
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

def save_min_lines_npz(angles, rmin, emin, filepath):
    folder = os.path.dirname(filepath)
    if folder:
        try: os.makedirs(folder, exist_ok=True)
        except Exception: pass
    np.savez_compressed(filepath, a=np.asarray(angles), r=np.asarray(rmin), e=np.asarray(emin))

def save_min_lines_gnuplot(angles, rmin, emin, filepath):
    folder = os.path.dirname(filepath)
    if folder:
        try: os.makedirs(folder, exist_ok=True)
        except Exception: pass
    with open(filepath, 'w') as f:
        f.write("# angle  rmin  emin\n")
        for a, r, e in zip(angles, rmin, emin):
            rs = "nan" if (r is None or not np.isfinite(r)) else f"{r:.12g}"
            es = "nan" if (e is None or not np.isfinite(e)) else f"{e:.12g}"
            f.write(f"{a:.6g} {rs} {es}\n")

def save_data_single(G, angles, distances, save_data_prefix=None, save_fmt="both", kcal=False, suffix=None):
    if kcal:
        if G is not None: G *= ev2kcal

    if save_data_prefix is None:
        print("No save_data_prefix specified, nothing written.")
        return

    if suffix is None:
        print("No suffix specified, nothing written.")
        return

    # Ensure output folder exists
    folder = os.path.dirname(save_data_prefix)
    if folder:
        os.makedirs(folder, exist_ok=True)

    if G is not None:
        if save_fmt in ("both", "npz"):     save_grid_npz(    angles, distances, G, save_data_prefix + "_" + suffix + ".npz")
        if save_fmt in ("both", "gnuplot"): save_grid_gnuplot(angles, distances, G, save_data_prefix + "_" + suffix + ".dat")

def exp_weight_func(Erefs, a=1.0, alpha=3.0, Emin0=0.1 ):
    Emin = np.min(Erefs)
    return np.exp( -alpha*(Erefs-Emin)/( np.abs(Emin) + Emin0 ) )*a

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

    # Process each curve segment
    lens = []
    #segment_counter = 0
    for start, end in zip(all_splits[:-1], all_splits[1:]):
        segment = Erefs[start:end]
        n = len(segment)
        lens.append(n)
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
            #print("Emin=", Emin)
            if weight_func is None:
                weights[start:end] = 1.0
            else:
                weights[start:end] = weight_func( Erefs[start:end] )
        
        #for i in range(start, end):
        #    print(f"    {segment_counter:6d} {i:6d} {x0s[i]:12.6f} {Erefs[i]:15.8f} {weights[i]:10.6f}", flush=True)

        #segment_counter += 1
    
    return weights, lens

def plotDOFscans( iDOFs, xs, DOFnames, bEs=True, bFs=False,  title="plotDOFscans", bEvalSamples=True, bPrint=False, outfile=None ):
    plt.figure(figsize=(8,10.0))
    color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
    ncol = len(color_cycle)
    for iDOF in iDOFs:
        y = DOFs[iDOF]    # store backup value of this DOF
        #print( f"#======= DOF[{iDOF}]: {xs}" )
        Es,Fs = scanParam( iDOF, xs, bEvalSamples=bEvalSamples )   # do 1D scan
        #print( f"#======= fDOF[{iDOF}]: {Efs}" )

        #print( "iDOF", iDOF, DOFnames[iDOF], " Fitness       \n", Es )
        #print( "iDOF", iDOF, DOFnames[iDOF], " dFitness/dDOF \n", Fs )
        Fs_num, xs_num = numDeriv(xs,Es)
        if outfile is not None:
            fname = f"{outfile}_DOF{iDOF}_ana.dat"
            np.savetxt( fname, np.array([xs, Es,Fs]).T, fmt="%23.15g", header=f"DOF {iDOF} {DOFnames[iDOF]} x[a.u.]  E[eV]  F[eV/a.u.] " )
            fname = f"{outfile}_DOF{iDOF}_num.dat"
            np.savetxt( fname, np.array([xs_num, -Fs_num]).T, fmt="%23.15g", header=f"DOF {iDOF} {DOFnames[iDOF]} x[a.u.]  F[eV/a.u.] " )
        # take color from standard matplotlib color cycle
        c = color_cycle[iDOF % ncol]
        if bEs: 
            plt.subplot(2,1,1); plt.plot(xs,Es, '-', color=c, label="E "+DOFnames[iDOF] )       # plot 1D scan
        if bFs: 
            Fs_num, xs_num = numDeriv(xs,Es)
            plt.subplot(2,1,2); plt.plot(xs,Fs,    '-', lw=1.0, color=c, label="F_ana "+DOFnames[iDOF] )       # This is error in the E_O3 charge derivative
            plt.subplot(2,1,2); plt.plot(xs_num,-Fs_num, ':', lw=1.5, color=c, label="F_num "+DOFnames[iDOF] ) 
            if bPrint:
                print ( "# plotDOFscans DOF ", iDOF, DOFnames[iDOF], " dx= ", xs[1]-xs[0] ); 
                print ( "#  i          x              E            F_ana=-dE/dDOF     F_num          F_ana-F_num          // F_num=-(E[i+1]-E[i-1])/(x[i+1]-x[i-1])" ); 
                for i in range(1, len(xs)-1 ): print( " %3i %15.5f %15.5f %15.5f %15.5f %15.5f"%(i, xs[i], Es[i], Fs[i], -Fs_num[i-1], Fs[i]+Fs_num[i-1]) )
        DOFs[iDOF] = y    # restore

    plt.subplot(2,1,1);
    plt.legend()
    plt.xlabel("DOF value")
    plt.ylabel("E [eV]")   
    plt.grid()

    plt.subplot(2,1,2);
    plt.xlabel("DOF value")
    plt.ylabel("F [eV/a.u.]")   
    plt.legend()
    plt.grid()

    plt.suptitle( title )

def numDeriv( x, y):
    dx = x[2:]-x[:-2]
    dy = y[2:]-y[:-2]
    xs = x[1:-1]
    return dy/dx, xs
