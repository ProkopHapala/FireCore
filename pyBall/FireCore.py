from ast import arg
import numpy as np
from ctypes import *

# needed for export_bcna_table wrapper
from ctypes import byref
c_int, c_double, c_bool, c_float, c_char_p, c_bool, c_void_p, c_char_p
import ctypes as ct
import os
import sys

#sys.path.append('../')
#from pyMeta import cpp_utils 
#import cpp_utils
from . import cpp_utils_ as cpp_utils
#cpp_utils = cpp_utils_

lib = None

c_double_p = ct.POINTER(c_double)
c_int_p    = ct.POINTER(c_int)

def fortReal(a):
    ct.byref(ct.c_double(a))

def _np_as(arr,atype):
    if arr is None:
        return None
    else: 
        return arr.ctypes.data_as(atype)

cpp_utils.s_numpy_data_as_call = "_np_as(%s,%s)"

# ===== To generate Interfaces automatically from headers call:
header_strings = [
"void firecore_getCharges( double* charges )",
"void firecore_preinit( ) ",
"void firecore_set_lvs( double* lvs )" ,
"void firecore_init( int natoms, int* atomTypes, double* atomsPos )",
"void firecore_assembleH( int iforce, int Kscf, double* positions )",
"void firecore_solveH( double* k_temp, int ikpoint )",
"void firecore_updateCharges( double sigmatol, double* sigma )",
"void firecore_SCF( int nmax_scf, double* positions, int iforce  )",
"void firecore_evalForce( int nmax_scf, double* positions_, double* forces_, double* energies )",
"void firecore_relax( nmax_scf, positions_, forces_, fixPos, energies )",
"void firecore_setupGrid( double Ecut, it ifixg0, doube* g0,  int ngrid, double* dCell  )", 
"void firecore_getGridMO( int iMO, double* ewfaux )",
"void firecore_getGridDens( int imo0, int imo1, double* ewfaux )",
"void firecore_getPointer_wfcoef( double* bbnkre )",
"void firecore_get_wfcoef( int ikp, double* wfcoefs )",
"void firecore_set_wfcoef( int iMO, int ikp, double* wfcoefs )",
"void firecore_getpsi( int in1, int issh, int n, double x0, double dx, double ys )",
"void firecore_MOtoXSF( int iMO )",
"void firecore_orb2points( int iband, int ikpoint, int npoints, double* points, double* ewfaux )",
"void firecore_get_HS_dims( int* natoms_out, int* norbitals_out, int* nspecies_out, int* neigh_max_out, int* numorb_max_out, int* nsh_max_out, int* ME2c_max_out, int* max_mu_dim1_out, int* max_mu_dim2_out, int* max_mu_dim3_out, int* mbeta_max_out, int* nspecies_fdata_out, int* nelec_out)",
    "void firecore_get_HS_neighs( int* num_orb_out, int* degelec_out, int* iatyp_out, int* lssh_out, int* mu_out, int* nu_out, int* mvalue_out, int* nssh_out, int* nzx_out, int* neighn_out, int* neigh_j_out, int* neigh_b_out, double* xl_out )",
    "void firecore_get_HS_sparse( double* h_mat_out, double* s_mat_out )",
    "void firecore_get_rho_sparse( double* rho_out )",
"void firecore_get_HS_k( double* kpoint_vec, void* Hk_out, void* Sk_out )",
"void firecore_get_nspecies( int* nspecies_out )",
"void firecore_get_eigen( int ikp, double* eigen_out )",
"void firecore_set_options( int ioff_S, int ioff_T, int ioff_Vna, int ioff_Vnl, int ioff_Vxc, int ioff_Vca, int ioff_Vxc_ca )",
"void firecore_set_export_mode( int export_mode )",
"void firecore_scanHamPiece2c( int interaction, int isub, int in1, int in2, int in3, double* dR, int applyRotation, double* sx_out )",
"void firecore_scanHamPiece3c( int interaction, int isorp, int in1, int in2, int indna, double* dRj, double* dRk, int applyRotation, double* bcnax_out )",
]
#cpp_utils.writeFuncInterfaces( header_strings );        exit()     #   uncomment this to re-generate C-python interfaces

# LFLAGS   = " -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lfftw3xf_gnu -lm -Bdynamic -lmpi "


mkldir = "/home/prokop/SW/intel/compilers_and_libraries_2019.5.281/linux/mkl/lib/intel64_lin/"

argDict={}

def reload():
    global lib
    bReaload = False
    if lib is not None: 
        cpp_utils.unload_lib(lib)
        bReaload = True
    cpp_utils.BUILD_PATH = os.path.normpath( cpp_utils.PACKAGE_PATH + '/../build/' ) 
    #lib = cpp_utils.loadLib('FireCore', recompile=False, mode=ct.RTLD_GLOBAL )
    lib = cpp_utils.loadLib('FireCore', recompile=False, mode=ct.RTLD_LOCAL )
    if bReaload:
        cpp_utils.set_args_dict(lib, argDict)
    return lib

lib = reload()

'''
loadLib( mkldir+"libmkl_def.so"        )
loadLib( mkldir+"libmkl_core.so"       )
loadLib( mkldir+"libmkl_avx2.so", )
'''

array1ui = np.ctypeslib.ndpointer(dtype=np.uint32, ndim=1, flags='CONTIGUOUS')
array1i  = np.ctypeslib.ndpointer(dtype=np.int32,  ndim=1, flags='CONTIGUOUS')
array2i  = np.ctypeslib.ndpointer(dtype=np.int32,  ndim=2, flags='CONTIGUOUS')
array1d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS')
array2d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=2, flags='CONTIGUOUS')
array3d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=3, flags='CONTIGUOUS')
array3i  = np.ctypeslib.ndpointer(dtype=np.int32,  ndim=3, flags='CONTIGUOUS')
array4d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=4, flags='CONTIGUOUS')
array2cd = np.ctypeslib.ndpointer(dtype=np.complex128, ndim=2, flags='CONTIGUOUS') # 
# ========= C functions



#  void firecore_setVerbosity( int verbosity_, int idebugWrite_ )
# lib.firecore_setVerbosity.argtypes  = [c_int, c_int ] 
# lib.firecore_setVerbosity.restype   =  None
argDict["firecore_setVerbosity"]=( None, [c_int, c_int ] )
def setVerbosity( verbosity=0, idebugWrite=0 ):
    return lib.firecore_setVerbosity(verbosity, idebugWrite ) 

#  subroutine firecore_init( natoms_, atomTypes, atomsPos )
# lib.firecore_init.argtypes  = [c_int, array1i, array2d ] 
# lib.firecore_init.restype   =  None
argDict["firecore_init"]=( None, [c_int, array1i, array2d ]  )
def init(atomTypes, atomPos ):
    #atomTypes = np.array( atomTypes , np.int32   )
    #atomPos   = np.array( atomPos   , np.float64 )
    #print("firecore_init() atomTypes: ", type(atomTypes), atomTypes.dtype )
    #print("firecore_init() atomPos  : ", type(atomPos),   atomPos.dtype )
    #print("firecore_init() \n atomTypes:" , atomTypes, "\n atomPos:", atomPos )
    natoms = len(atomTypes)
    return lib.firecore_init(natoms, atomTypes, atomPos )

#  subroutine firecore_evalForce( nmax_scf, forces_ )
#lib.firecore_evalForce.argtypes  = [c_int, array2d, array2d, array1d, c_int ] 
#lib.firecore_evalForce.restype   =  None
argDict["firecore_evalForce"]=( None, [c_int, array2d, array2d, array1d, c_int ] )
def evalForce( pos, forces=None, nmax_scf=100, Es=None, ixyz=-1 ):
    if Es     is None: Es     = np.zeros(8)
    if forces is None: forces = np.zeros( pos.shape )
    lib.firecore_evalForce( nmax_scf, pos, forces, Es, ixyz )
    return forces, Es

# "void firecore_relax( nmax_scf, positions_, forces_, fixPos, energies )",
#lib.firecore_relax.argtypes  = [c_int, c_int, array2d, array2d, array2i, array1d ] 
#lib.firecore_relax.restype   =  None
argDict["firecore_relax"]=( None, [c_int, c_int, array2d, array2d, array2i, array1d ] )
def relax( pos, forces=None, fixPos=None, nstepf=1000, nmax_scf=100, Es=None, ):
    if Es     is None: Es     = np.zeros(8)
    if forces is None: forces = np.zeros( pos.shape )
    if fixPos is None: 
        fixPos = np.zeros( pos.shape, np.int32 )
    else:
        if( isinstance(fixPos[0],int) ):
            fixPos_=np.zeros(pos.shape,np.int32)
            if(fixPos[0]>=0):
                fixPos_[fixPos,:]=1; 
            fixPos=fixPos_
        else:
            fixPos=np.array(fixPos, np.int32)
    lib.firecore_relax( nstepf, nmax_scf, pos, forces, fixPos, Es )
    return forces

#  void getCharges( double* charges )
#lib.firecore_getCharges.argtypes  = [array2d] 
#lib.firecore_getCharges.restype   =  None
argDict["firecore_getCharges"]=( None, [array2d] )
def getCharges(charges):
    return lib.firecore_getCharges(charges) 

#  "void firecore_getPointer_wfcoef( double* bbnkre )
#lib.firecore_getPointer_wfcoef.argtypes  = [array2d] 
#lib.firecore_getPointer_wfcoef.restype   =  None
argDict["firecore_getPointer_wfcoef"]=( None, [array2d] )
def getPointer_wfcoef(charges):
    return lib.firecore_getPointer_wfcoef(charges) 

#"void firecore_get_wfcoef( int ikp, double* wfcoefs )",
#  void getCharges( double* charges )
#lib.firecore_get_wfcoef.argtypes  = [c_int, array2d] 
#lib.firecore_get_wfcoef.restype   =  None
argDict["firecore_get_wfcoef"]=( None, [c_int, array2d] )
def get_wfcoef(wfcoef=None,norb=None, ikp=1):
    if(wfcoef is None):
        wfcoef=np.zeros( (norb,norb) )
    lib.firecore_get_wfcoef(ikp,wfcoef) 
    return wfcoef

#"void firecore_set_wfcoef( int iMO, int ikp, double* wfcoefs )",
#lib.firecore_set_wfcoef.argtypes  = [c_int,c_int, array1d] 
#lib.firecore_set_wfcoef.restype   =  None
argDict["firecore_set_wfcoef"]=( None, [c_int,c_int, array1d] )
def set_wfcoef(wfcoef,iMO=1,ikp=1):
    return lib.firecore_set_wfcoef(iMO,ikp,wfcoef)

#  void preinit( ) 
#lib.firecore_preinit.argtypes  = [] 
#lib.firecore_preinit.restype   =  None
argDict["firecore_preinit"]=( None, [] )
def preinit():
    return lib.firecore_preinit() 

#  void set_lvs( double* lvs )
#lib.firecore_set_lvs.argtypes  = [array2d] 
#lib.firecore_set_lvs.restype   =  None
argDict["firecore_set_lvs"]=( None, [array2d] )
def set_lvs(lvs):
    return lib.firecore_set_lvs(lvs) 

#  void assembleH( int iforce, int Kscf, double* positions )
#lib.firecore_assembleH.argtypes  = [c_int, c_int, array2d ] 
#lib.firecore_assembleH.restype   =  None
argDict["firecore_assembleH"]=( None, [c_int, c_int, array2d ] )
def assembleH( positions, iforce=0, Kscf=1 ):
    return lib.firecore_assembleH(iforce, Kscf, positions) 

#  void solveH( double* k_temp, int ikpoint )
#lib.firecore_solveH.argtypes  = [ array1d, c_int] 
#lib.firecore_solveH.restype   =  None
argDict["firecore_solveH"]=( None, [ array1d, c_int] )
def solveH(k_temp=None, ikpoint=1):
    if k_temp is None:
        k_temp = np.array([0.0,0.0,0.0])
    return lib.firecore_solveH(k_temp, ikpoint ) 

#  void updateCharges( double sigmatol, double* sigma )
#lib.firecore_updateCharges.argtypes  = [c_double, c_double_p] 
#lib.firecore_updateCharges.restype   =  None
argDict["firecore_updateCharges"]=( None, [c_double, c_double_p] )
def updateCharges( sigmatol=1e-6, sigma=None):
    if sigma is None:
        sigma = np.zeros(1)
    lib.firecore_updateCharges(sigmatol, _np_as(sigma,c_double_p)) 
    return sigma

#  void SCF( int nmax_scf, double* positions, int iforce  )
#lib.firecore_SCF.argtypes  = [c_int, array2d, c_int] 
#lib.firecore_SCF.restype   =  None
argDict["firecore_SCF"]=( None, [c_int, array2d, c_int] )
def SCF( positions, iforce=0, nmax_scf=200 ):
    return lib.firecore_SCF( nmax_scf, positions, iforce ) 

#  void setupGrid( double Ecut, int ifixg0, doube* g0,  int ngrid, double* dCell  )
#lib.firecore_setupGrid.argtypes  = [c_double, c_int, array1d, array1i, array2d ] 
#lib.firecore_setupGrid.restype   =  None
argDict["firecore_setupGrid"]=( None, [c_double, c_int, array1d, array1i, array2d ] )
def setupGrid(Ecut=100, g0=None, ngrid=None, dCell=None):
    if g0 is None:
        ifixg0=0
        g0=np.array([0.0,0.0,0.0])
    else:
        ifixg0=1
        g0 = np.array( g0 )
    if ngrid is None:
        ngrid = np.zeros(3,dtype=np.int32); ngrid[:]=-1
    else:
        ngrid = np.array(ngrid,dtype=np.int32)
    if dCell is None:
        dCell = np.zeros( (3,3) )
    else:
        dCell = np.array( dCell, dtype=np.float )
    lib.firecore_setupGrid(Ecut, ifixg0, g0, ngrid, dCell )
    lvs = np.array( [ [0.0,0.0,0.0], dCell[0]*ngrid[0], dCell[1]*ngrid[1], dCell[2]*ngrid[2], ] )
    #ngrid = ngrid[::-1].copy()
    return ngrid, dCell, lvs

#  void getGridMO( int iMO, double* ewfaux )
#lib.firecore_getGridMO.argtypes  = [ c_int, array3d ] 
#lib.firecore_getGridMO.restype   =  None
argDict["firecore_getGridMO"]=( None, [ c_int, array3d ] )
def getGridMO(iMO, ewfaux=None, ngrid=None):
    if ewfaux is None:
        ewfaux = np.zeros(ngrid)
    lib.firecore_getGridMO(iMO, ewfaux )
    return ewfaux

#  void getGridDens( int imo0, int imo1, double* ewfaux )
#lib.firecore_getGridDens.argtypes  = [ array3d, c_double, c_double ] 
#lib.firecore_getGridDens.restype   =  None
argDict["firecore_getGridDens"]=( None, [ array3d, c_double, c_double ] )
def getGridDens(ewfaux=None, ngrid=None, Cden = 1.0, Cden0 = 0.0 ):
    ngrid=ngrid[::-1]
    print( " getGridDens() ngrid ", ngrid ); #exit(0)
    if ewfaux is None:
        ewfaux = np.zeros(ngrid)
    lib.firecore_getGridDens( ewfaux, Cden, Cden0  ) 
    return ewfaux 

#  void firecore_MOtoXSF( int iMO )
#lib.firecore_orb2xsf.argtypes  = [ c_int ] 
#lib.firecore_orb2xsf.restype   =  None
argDict["firecore_orb2xsf"]=( None, [ c_int ] )
def orb2xsf(iMO):
    lib.firecore_orb2xsf(iMO )

#  void firecore_dens2xsf( int iMO )
#lib.firecore_dens2xsf.argtypes  = [ c_double ] 
#lib.firecore_dens2xsf.restype   =  None
argDict["firecore_dens2xsf"]=( None, [ c_double ] )
def dens2xsf( f_den0=0.0 ):
    print( "DEBUG FireCore.py:: dens2xsf()" )
    lib.firecore_dens2xsf( f_den0 )
    print( "DEBUG FireCore.py:: dens2xsf() END" )

'''
#"void firecore_getpsi( int in1, int issh, int n, double x0, double dx, double ys )",
lib.firecore_getpsi.argtypes  = [c_int, c_int, c_double, c_double,  c_int, c_int, c_int, c_double, c_double, array1d ] 
lib.firecore_getpsi.restype   =  None
def getpsi(in1=1, issh=1, n=50, dx=0.1, x0=0.0, ys=None, l=0, m=1, theta=0.0, phi=0.0):
    if ys is None:
        ys = np.zeros(n)
    lib.firecore_getpsi(l,m,theta,phi,in1, issh, n, x0, dx, ys )
    return ys
'''

#"void firecore_getpsi( int in1, int issh, int n, double* poss, double* ys )",
#lib.firecore_getpsi.argtypes  = [c_int, c_int, c_int, c_int,    c_int, array2d, array1d ] 
#lib.firecore_getpsi.restype   =  None
argDict["firecore_getpsi"]=( None, [c_int, c_int, c_int, c_int,    c_int, array2d, array1d ] )
def getpsi( poss, ys=None, in1=1, issh=1, l=0, m=1 ):
    n = len(poss)
    if ys is None:
        ys = np.zeros(n)
    lib.firecore_getpsi( l, m, in1, issh, n, poss, ys )
    return ys

#void firecore_orb2points( int iband, int ikpoint, int npoints, double* points, double* ewfaux )
#lib.firecore_orb2points.argtypes  = [c_int, c_int, c_int, array2d, array1d ] 
#lib.firecore_orb2points.restype   =  None
argDict["firecore_orb2points"]=( None, [c_int, c_int, c_int, array2d, array1d ] )
def orb2points( poss, ys=None, iMO=1,  ikpoint=1 ):
    n = len(poss)
    if ys is None:
        ys = np.zeros(n)
    lib.firecore_orb2points( iMO, ikpoint, n, poss, ys )
    return ys

# subroutine firecore_dens2points(npoints, points, f_den, f_den0, ewfaux_out) bind(c, name='firecore_dens2points')
argDict["firecore_dens2points"]=( None, [c_int, array2d, c_double, c_double, array1d ] )
def dens2points( points, f_den=1.0, f_den0=0.0, ewfaux_out=None ):
    n = len(points)
    if ewfaux_out is None:
        ewfaux_out = np.zeros(n)
    lib.firecore_dens2points( n, points, f_den, f_den0, ewfaux_out )
    return ewfaux_out

# eigenvalues export: firecore_get_eigen( int ikp, double* eigen_out )
argDict["firecore_get_eigen"]=( None, [c_int, array1d] )
def get_eigen( ikp=1, norb=None ):
    """Return array of orbital energies eigen_k(:, ikp) for given k-point index (1-based)."""
    if norb is None:
        dims = get_HS_dims()
        norb = dims.norbitals
    eigen = np.zeros(norb, dtype=np.float64)
    lib.firecore_get_eigen( ikp, eigen )
    return eigen

# --- Export H and S matrices ---

# Define classes to hold the structured data
class FireballDims:
    def __init__(self, natoms, norbitals, nspecies, neigh_max, numorb_max, nsh_max, ME2c_max,
                 max_mu_dim1, max_mu_dim2, max_mu_dim3, mbeta_max, nspecies_fdata, nelec):
        self.natoms = natoms
        self.norbitals = norbitals
        self.nspecies = nspecies  # Distinct species in current calculation
        self.neigh_max = neigh_max
        self.numorb_max = numorb_max
        self.nsh_max = nsh_max
        self.ME2c_max = ME2c_max
        self.max_mu_dim1 = max_mu_dim1
        self.max_mu_dim2 = max_mu_dim2
        self.max_mu_dim3 = max_mu_dim3
        self.mbeta_max = mbeta_max
        self.nspecies_fdata = nspecies_fdata # Total species in info.dat
        self.nelec = nelec  # Total electron count (from ztot)

class FireballData:
    # def __init__(self, dims: FireballDims):
    #     if not isinstance(dims, FireballDims):
    #         raise TypeError("dims argument must be an instance of FireballDims")

    #     # Fortran: rho(numorb_max, numorb_max, neigh_max, natoms)
    #     # Python:  rho[imu, inu, ineigh, iatom]
    #     self.rho     = np.zeros((dims.numorb_max, dims.numorb_max, dims.neigh_max, dims.natoms), dtype=np.float64, order='F')
    #     self.h_mat   = np.zeros((dims.numorb_max, dims.numorb_max, dims.neigh_max, dims.natoms), dtype=np.float64, order='F')
    #     self.s_mat   = np.zeros((dims.numorb_max, dims.numorb_max, dims.neigh_max, dims.natoms), dtype=np.float64, order='F')
    #     self.num_orb = np.zeros(dims.nspecies, dtype=np.int32, order='F')
    #     self.degelec = np.zeros(dims.natoms, dtype=np.int32, order='F')
    #     self.iatyp   = np.zeros(dims.natoms, dtype=np.int32, order='F')
    #     # Fortran: lssh(nsh_max, nspecies) -> Python: lssh[ish, ispec]
    #     self.lssh    = np.zeros((dims.nsh_max, dims.nspecies), dtype=np.int32, order='F')
    #     # Fortran: mu(d1,d2,d3) -> Python: mu[d1,d2,d3]
    #     self.mu      = np.zeros((dims.max_mu_dim1, dims.max_mu_dim2, dims.max_mu_dim3), dtype=np.int32, order='F')
    #     self.nu      = np.zeros((dims.max_mu_dim1, dims.max_mu_dim2, dims.max_mu_dim3), dtype=np.int32, order='F')
    #     self.mvalue  = np.zeros((dims.max_mu_dim1, dims.max_mu_dim2, dims.max_mu_dim3), dtype=np.int32, order='F')
    #     self.nssh    = np.zeros(dims.nspecies, dtype=np.int32, order='F')
    #     self.nzx     = np.zeros(dims.nspecies_fdata, dtype=np.int32, order='F')
    #     self.neighn  = np.zeros(dims.natoms, dtype=np.int32, order='F')
    #     # Fortran: neigh_j(neigh_max, natoms) -> Python: neigh_j[ineigh, iatom]
    #     self.neigh_j = np.zeros((dims.neigh_max, dims.natoms), dtype=np.int32, order='F')
    #     self.neigh_b = np.zeros((dims.neigh_max, dims.natoms), dtype=np.int32, order='F')
    #     # Fortran: xl(3, mbeta_max) -> Python: xl[i, mbeta]
    #     self.xl      = np.zeros((3, dims.mbeta_max), dtype=np.float64, order='F')

    # NOTE: We intentionally keep Python arrays in default C-order (row-major)
    # and reverse axis order so the linear layout matches Fortran column-major.
    # Example: Fortran A(i,j,k,l) -> Python shape (l,k,j,i); last index is fastest.
    # Do NOT use order='F' because ctypes ndpointer requires C_CONTIGUOUS.
    def __init__(self, dims: FireballDims):
        if not isinstance(dims, FireballDims):
            raise TypeError("dims argument must be an instance of FireballDims")

        # Fortran: rho(numorb_max, numorb_max, neigh_max, natoms)
        # Python (axis-reversed for matching linear layout): rho[iatom, ineigh, imu, inu]
        self.neigh_max = dims.neigh_max
        self.numorb_max = dims.numorb_max
        self.rho     = np.zeros((dims.natoms, dims.neigh_max, dims.numorb_max, dims.numorb_max), dtype=np.float64)
        self.h_mat   = np.zeros((dims.natoms, dims.neigh_max, dims.numorb_max, dims.numorb_max), dtype=np.float64)
        self.s_mat   = np.zeros((dims.natoms, dims.neigh_max, dims.numorb_max, dims.numorb_max), dtype=np.float64)
        # NOTE: iatyp returned by Fortran is the atomic number (Z), not a compact
        # species index. To avoid out-of-bounds when iatyp is used directly, size
        # num_orb by the full nspecies_fdata table (max species defined in Fdata).
        self.num_orb = np.zeros(dims.nspecies_fdata, dtype=np.int32)
        self.degelec = np.zeros(dims.natoms, dtype=np.int32)
        self.iatyp   = np.zeros(dims.natoms, dtype=np.int32)
        # Fortran: lssh(nsh_max, nspecies) -> Python: lssh[ish, ispec]
        self.lssh    = np.zeros((dims.nspecies, dims.nsh_max), dtype=np.int32)
        # Fortran: mu(d1,d2,d3) -> Python: mu[d1,d2,d3]
        self.mu      = np.zeros((dims.max_mu_dim3, dims.max_mu_dim2, dims.max_mu_dim1), dtype=np.int32)
        self.nu      = np.zeros((dims.max_mu_dim3, dims.max_mu_dim2, dims.max_mu_dim1), dtype=np.int32)
        self.mvalue  = np.zeros((dims.max_mu_dim3, dims.max_mu_dim2, dims.max_mu_dim1), dtype=np.int32)
        self.nssh    = np.zeros(dims.nspecies, dtype=np.int32)
        self.nzx     = np.zeros(dims.nspecies_fdata, dtype=np.int32)
        self.neighn  = np.zeros(dims.natoms, dtype=np.int32)
        # Fortran: neigh_j(neigh_max, natoms) -> Python: neigh_j[ineigh, iatom]
        self.neigh_j = np.zeros((dims.natoms, dims.neigh_max), dtype=np.int32)
        self.neigh_b = np.zeros((dims.natoms, dims.neigh_max), dtype=np.int32)
        # Fortran: xl(3, mbeta_max) -> Python reversed: xl[mbeta, 3]
        self.xl      = np.zeros((dims.mbeta_max, 3), dtype=np.float64)

# void firecore_get_HS_dims( int* natoms_out, int* norbitals_out, int* nspecies_out, int* neigh_max_out, int* numorb_max_out, int* nsh_max_out, int* ME2c_max_out, int* max_mu_dim1_out, int* max_mu_dim2_out, int* max_mu_dim3_out, int* mbeta_max_out, int* nspecies_fdata_out, int* nelec_out)
argDict["firecore_get_HS_dims"]=( None, [c_int_p, c_int_p, c_int_p, c_int_p, c_int_p, c_int_p, c_int_p, c_int_p, c_int_p, c_int_p, c_int_p, c_int_p, c_int_p] ) # 13 args
def get_HS_dims():
    natoms_c         = ct.c_int()
    norbitals_c      = ct.c_int()
    nspecies_c       = ct.c_int()
    neigh_max_c      = ct.c_int()
    numorb_max_c     = ct.c_int()
    nsh_max_c        = ct.c_int()
    ME2c_max_c       = ct.c_int()
    max_mu_dim1_c    = ct.c_int()
    max_mu_dim2_c    = ct.c_int()
    max_mu_dim3_c    = ct.c_int()
    mbeta_max_c      = ct.c_int()
    nspecies_fdata_c = ct.c_int()
    nelec_c          = ct.c_int()
    lib.firecore_get_HS_dims(
        ct.byref(natoms_c), ct.byref(norbitals_c), ct.byref(nspecies_c),
        ct.byref(neigh_max_c), ct.byref(numorb_max_c),
        ct.byref(nsh_max_c), ct.byref(ME2c_max_c),
        ct.byref(max_mu_dim1_c), ct.byref(max_mu_dim2_c), ct.byref(max_mu_dim3_c),
        ct.byref(mbeta_max_c), ct.byref(nspecies_fdata_c), ct.byref(nelec_c)
    )
    return FireballDims(
        natoms_c.value, norbitals_c.value,
        nspecies_c.value,
        neigh_max_c.value, numorb_max_c.value,
        nsh_max_c.value, ME2c_max_c.value,
        max_mu_dim1_c.value, max_mu_dim2_c.value, max_mu_dim3_c.value,
        mbeta_max_c.value, nspecies_fdata_c.value, nelec_c.value
    )

# void firecore_get_HS_neighs( int* num_orb_out, int* degelec_out, int* iatyp_out, int* lssh_out, int* mu_out, int* nu_out, int* mvalue_out, int* nssh_out, int* nzx_out, int* neighn_out, int* neigh_j_out, int* neigh_b_out, double* xl_out )
argDict["firecore_get_HS_neighs"]=( None, [array1i, array1i, array1i, array2i, array3i, array3i, array3i, array1i, array1i, array1i, array2i, array2i, array2d] )
def get_HS_neighs(dims):
    if not isinstance(dims, FireballDims):
        raise TypeError("dims argument must be an instance of FireballDims")
    data = FireballData(dims)
    # Ensure C-contiguous arrays for ctypes ndpointer (flags='C_CONTIGUOUS')
    def c_array(a, name):
        arr = np.ascontiguousarray(a, dtype=a.dtype)
        assert arr.flags['C_CONTIGUOUS'], f"{name} not C-contiguous"
        return arr
    data.num_orb = c_array(data.num_orb, "num_orb")
    data.degelec = c_array(data.degelec, "degelec")
    data.iatyp   = c_array(data.iatyp,   "iatyp")
    data.lssh    = c_array(data.lssh,    "lssh")
    data.mu      = c_array(data.mu,      "mu")
    data.nu      = c_array(data.nu,      "nu")
    data.mvalue  = c_array(data.mvalue,  "mvalue")
    data.nssh    = c_array(data.nssh,    "nssh")
    data.nzx     = c_array(data.nzx,     "nzx")
    data.neighn  = c_array(data.neighn,  "neighn")
    data.neigh_j = c_array(data.neigh_j, "neigh_j")
    data.neigh_b = c_array(data.neigh_b, "neigh_b")
    data.xl      = c_array(data.xl,      "xl")
    lib.firecore_get_HS_neighs(
        data.num_orb, data.degelec, data.iatyp,
        data.lssh, data.mu, data.nu, data.mvalue, data.nssh, data.nzx,
        data.neighn, data.neigh_j, data.neigh_b, data.xl
    )
    return data

# void firecore_get_HS_sparse( double* h_mat_out, double* s_mat_out )
argDict["firecore_get_HS_sparse"]=( None, [array4d, array4d] )
def get_HS_sparse(dims, data=None):
    if not isinstance(dims, FireballDims):
        raise TypeError("dims argument must be an instance of FireballDims")
    if data is None: data = FireballData(dims)
    lib.firecore_get_HS_sparse( data.h_mat, data.s_mat )
    return data

# void firecore_get_rho_sparse( double* rho_out )
argDict["firecore_get_rho_sparse"]=( None, [array4d] )
def get_rho_sparse(dims, data=None):
    if not isinstance(dims, FireballDims):
        raise TypeError("dims argument must be an instance of FireballDims")
    if data is None: data = FireballData(dims)
    lib.firecore_get_rho_sparse( data.rho )
    return data

# void firecore_get_HS_k( double* kpoint_vec, void* Hk_out, void* Sk_out ) # Using void* for complex arrays
argDict["firecore_get_HS_k"]=( None, [array1d, array2cd, array2cd] ) # array2cd for complex double
def get_HS_k(kpoint_vec, norbitals):
    Hk_out = np.zeros((norbitals, norbitals), dtype=np.complex128)
    Sk_out = np.zeros((norbitals, norbitals), dtype=np.complex128)
    kpoint_vec_np = np.array(kpoint_vec, dtype=np.float64)
    lib.firecore_get_HS_k(kpoint_vec_np, Hk_out, Sk_out)
    return Hk_out, Sk_out

cpp_utils.set_args_dict(lib, argDict)

#===================================
# ========= Python Functions
#===================================

_cached_norb = None
_cached_dims = None

def _get_norb(norb=None):
    global _cached_norb
    if norb is not None:
        _cached_norb = norb
        return norb
    if _cached_norb is None:
        _cached_norb = get_HS_dims().numorb_max
    return _cached_norb

def run_nonSCF( atomType, atomPos ):
    preinit()
    norb = init( atomType, atomPos )
    assembleH( atomPos)
    solveH()
    sigma=updateCharges(); #print( sigma )
    return norb, sigma

def _load_HS_dims():
    natoms_c         = ct.c_int()
    norbitals_c      = ct.c_int()
    nspecies_c       = ct.c_int()
    neigh_max_c      = ct.c_int()
    numorb_max_c     = ct.c_int()
    nsh_max_c        = ct.c_int()
    ME2c_max_c       = ct.c_int()
    max_mu_dim1_c    = ct.c_int()
    max_mu_dim2_c    = ct.c_int()
    max_mu_dim3_c    = ct.c_int()
    mbeta_max_c      = ct.c_int()
    nspecies_fdata_c = ct.c_int()
    nelec_c          = ct.c_int()
    lib.firecore_get_HS_dims(
        ct.byref(natoms_c), ct.byref(norbitals_c), ct.byref(nspecies_c),
        ct.byref(neigh_max_c), ct.byref(numorb_max_c),
        ct.byref(nsh_max_c), ct.byref(ME2c_max_c),
        ct.byref(max_mu_dim1_c), ct.byref(max_mu_dim2_c), ct.byref(max_mu_dim3_c),
        ct.byref(mbeta_max_c), ct.byref(nspecies_fdata_c), ct.byref(nelec_c)
    )
    return FireballDims(
        natoms_c.value, norbitals_c.value,
        nspecies_c.value,
        neigh_max_c.value, numorb_max_c.value,
        nsh_max_c.value, ME2c_max_c.value,
        max_mu_dim1_c.value, max_mu_dim2_c.value, max_mu_dim3_c.value,
        mbeta_max_c.value, nspecies_fdata_c.value, nelec_c.value
    )

def get_HS_dims(force_refresh=False):
    global _cached_dims
    if force_refresh or (_cached_dims is None):
        _cached_dims = _load_HS_dims()
    return _cached_dims


def set_options(ioff_S=1, ioff_T=1, ioff_Vna=1, ioff_Vnl=1, ioff_Vxc=1, ioff_Vca=1, ioff_Vxc_ca=1):
    lib.firecore_set_options(ioff_S, ioff_T, ioff_Vna, ioff_Vnl, ioff_Vxc, ioff_Vca, ioff_Vxc_ca)

argDict["firecore_set_export_mode"]=( None, [c_int] )
def set_export_mode(export_mode=0):
    return lib.firecore_set_export_mode(int(export_mode))

# void firecore_scanHamPiece2c( int interaction, int isub, int in1, int in2, int in3, double* dR, int applyRotation, double* sx_out )
argDict["firecore_scanHamPiece2c"]=( None, [c_int, c_int, c_int, c_int, c_int, c_double_p, c_int, c_double_p] )
argDict["firecore_scanHamPiece2c_batch"]=( None, [c_int, c_int, c_int, c_int, c_int, c_int, c_double_p, c_int, c_double_p] )
def scanHamPiece2c( interaction, isub, in1, in2, in3, dR, applyRotation=True, norb=None, sx_out=None ):
    if sx_out is None:
        norb = _get_norb(norb)
        sx_out = np.zeros( (norb, norb) )
    lib.firecore_scanHamPiece2c( interaction, isub, in1, in2, in3, _np_as(dR,c_double_p), int(applyRotation), _np_as(sx_out,c_double_p) )
    return sx_out

def scanHamPiece2c_batch( interaction, isub, in1, in2, in3, dRs, applyRotation=True, norb=None, sx_out=None ):
    dRs_np = np.ascontiguousarray(dRs, dtype=np.float64)
    npoints = dRs_np.shape[0]
    if sx_out is None:
        norb = _get_norb(norb)
        sx_out = np.zeros( (npoints, norb, norb) )
    lib.firecore_scanHamPiece2c_batch( interaction, isub, in1, in2, in3, npoints, _np_as(dRs_np.ravel(), c_double_p), int(applyRotation), _np_as(sx_out, c_double_p) )
    return sx_out

# void firecore_scanHamPiece3c( int interaction, int isorp, int in1, int in2, int indna, double* dRj, double* dRk, int applyRotation, double* bcnax_out )
argDict["firecore_scanHamPiece3c"]=( None, [c_int, c_int, c_int, c_int, c_int, c_double_p, c_double_p, c_int, c_double_p] )
argDict["firecore_scanHamPiece3c_batch"]=( None, [c_int, c_int, c_int, c_int, c_int, c_int, c_double_p, c_double_p, c_int, c_double_p] )
# raw 3c (no Legendre/recover/rotate)
argDict["firecore_scanHamPiece3c_raw"]=( None, [c_int, c_int, c_int, c_int, c_int, c_double_p, c_double_p, c_double_p] )
argDict["firecore_scanHamPiece3c_raw_batch"]=( None, [c_int, c_int, c_int, c_int, c_int, c_int, c_double_p, c_double_p, c_double_p] )
# export raw bcna table (itheta 1..5), flattened x-major (ix fastest) then y
argDict["firecore_export_bcna_table"]=( None, [c_int, c_int, c_int, c_int, c_int, c_int, c_int, c_int, c_int_p, c_int_p, c_double_p, c_double_p, c_double_p, c_int_p] )
def scanHamPiece3c( interaction, isorp, in1, in2, indna, dRj, dRk, applyRotation=True, norb=None, bcnax_out=None ):
    if bcnax_out is None:
        norb = _get_norb(norb)
        bcnax_out = np.zeros( (norb, norb) )
    lib.firecore_scanHamPiece3c( interaction, isorp, in1, in2, indna, _np_as(dRj,c_double_p), _np_as(dRk,c_double_p), int(applyRotation), _np_as(bcnax_out,c_double_p) )
    return bcnax_out

def scanHamPiece3c_raw( interaction, isorp, in1, in2, indna, dRj, dRk, out=None ):
    """
    Raw bcna interpolation: no Legendre, no mvalue*sint, no recover_3c, no rotate.
    Returns array shape (5, max_mu_dim3).
    """
    dims = get_HS_dims()
    max_me = dims.max_mu_dim3
    dRj_c = _np_as(dRj, c_double_p)
    dRk_c = _np_as(dRk, c_double_p)
    if out is None:
        out = np.zeros((5, max_me), dtype=np.float64, order='F')
    lib.firecore_scanHamPiece3c_raw(
        c_int(interaction), c_int(isorp), c_int(in1), c_int(in2), c_int(indna),
        dRj_c, dRk_c, _np_as(out, c_double_p)
    )
    return out

def scanHamPiece3c_raw_batch( interaction, isorp, in1, in2, indna, dRjs, dRks ):
    """
    Batched raw bcna interpolation over an array of geometries.
    Returns array shape (npoints, 5, max_mu_dim3).
    """
    dims = get_HS_dims()
    max_me = dims.max_mu_dim3
    dRjs_np = np.ascontiguousarray(dRjs, dtype=np.float64)
    dRks_np = np.ascontiguousarray(dRks, dtype=np.float64)
    npoints = dRjs_np.shape[0]
    if dRks_np.shape[0] != npoints:
        raise ValueError("dRjs and dRks must have the same length")
    dRjs_c = _np_as(dRjs_np.ravel(), c_double_p)
    dRks_c = _np_as(dRks_np.ravel(), c_double_p)
    out_raw = np.zeros((5, max_me, npoints), dtype=np.float64, order='F')
    lib.firecore_scanHamPiece3c_raw_batch(
        c_int(interaction), c_int(isorp), c_int(in1), c_int(in2), c_int(indna),
        c_int(npoints), dRjs_c, dRks_c, _np_as(out_raw, c_double_p)
    )
    return np.moveaxis(out_raw, 2, 0)

def export_bcna_table( interaction, isorp, in1, in2, indna, itheta, iME, maxsize ):
    """
    Export raw bcna_0{itheta} table from Fortran; returns (arr[ny,nx], hx, hy, status).
    Data layout in arr is [y,x] because Fortran filled buffer with ix fastest, then iy.
    """
    buf = np.zeros(maxsize, dtype=np.float64, order='C')
    nx_out = c_int(0); ny_out = c_int(0)
    hx_out = c_double(0.0); hy_out = c_double(0.0)
    status = c_int(0)
    lib.firecore_export_bcna_table(
        c_int(interaction), c_int(isorp), c_int(in1), c_int(in2), c_int(indna),
        c_int(itheta), c_int(iME), c_int(maxsize),
        byref(nx_out), byref(ny_out), byref(hx_out), byref(hy_out),
        _np_as(buf, c_double_p), byref(status)
    )
    if status.value != 0:
        raise RuntimeError(f"firecore_export_bcna_table failed: status={status.value}")
    nxy = nx_out.value * ny_out.value
    arr = buf[:nxy].reshape((ny_out.value, nx_out.value))
    return arr, float(hx_out.value), float(hy_out.value), status.value

def scanHamPiece3c_batch( interaction, isorp, in1, in2, indna, dRjs, dRks, applyRotation=True, norb=None, bcnax_out=None ):
    dRjs_np = np.ascontiguousarray(dRjs, dtype=np.float64)
    dRks_np = np.ascontiguousarray(dRks, dtype=np.float64)
    npoints = dRjs_np.shape[0]
    if bcnax_out is None:
        norb = _get_norb(norb)
        bcnax_out = np.zeros( (npoints, norb, norb) )
    lib.firecore_scanHamPiece3c_batch( interaction, isorp, in1, in2, indna, npoints, _np_as(dRjs_np.ravel(), c_double_p), _np_as(dRks_np.ravel(), c_double_p), int(applyRotation), _np_as(bcnax_out,c_double_p) )
    return bcnax_out

def initialize( atomType=None, atomPos=None, verbosity=1 ):
    global norb
    #if verbosity is not None: 
    setVerbosity(verbosity=verbosity)  # NOTE: verbosity must be set here, because previous values would be overwritten by preinit()
    preinit()
    norb = init( atomType, atomPos )
    _get_norb(norb)
    return norb


#===================================
# ========= Main
#===================================

if __name__ == "__main__":

    # b = sum2   ( 5.0 ); print "sum2   ", b
    # b = sum2val( 5.0,6.0 ); print "sum2val", b

    #a = ct.c_double(5)
    #b = lib.sum2(ct.byref(a))
    #print b

    #firecore_hello()
    
    natoms = 5
    #atomType = np.random.randint(6, size=natoms).astype(np.int32)
    #atomPos  = np.random.random((3,natoms))
    atomType = np.array([6,1,1,1,1]).astype(np.int32)
    atomPos  = np.array([
        [ 0.1,      0.0,     0.0],
        [-1.0,     +1.0,    -1.0],
        [+1.0,     -1.0,    -1.0],
        [-1.0,     -1.0,    +1.0],
        [+1.0,     +1.0,    +1.0],
    ])
    print ("atomType ", atomType)
    print ("atomPos  ", atomPos)
    preinit()
    init( natoms, atomType, atomPos )

    '''
    # =========== Molecular Orbital
    assembleH( atomPos)
    solveH()
    ngrid, dCell = setupGrid()
    ewfaux = getGridMO( 2, ngrid=ngrid )
    #ewfaux = getGridDens( ngrid=ngrid )
    print( ewfaux.min(),ewfaux.max() )
    import matplotlib.pyplot as plt
    sh = ewfaux.shape
    plt.figure(); plt.imshow( ewfaux[ sh[0]/2+5,:,: ] )
    plt.figure(); plt.imshow( ewfaux[ sh[0]/2  ,:,: ] )
    plt.figure(); plt.imshow( ewfaux[ sh[0]/2-5,:,: ] )
    plt.show()
    '''

    # =========== Electron Density
    assembleH( atomPos)
    solveH()
    sigma=updateCharges() ; print( sigma )
    ngrid, dCell, lvs = setupGrid()
    ewfaux = getGridDens( ngrid=ngrid )
    print( ewfaux.min(),ewfaux.max() )
    
    import matplotlib.pyplot as plt
    sh = ewfaux.shape
    plt.figure(); plt.imshow( ewfaux[ sh[0]/2+5,:,: ] )
    plt.figure(); plt.imshow( ewfaux[ sh[0]/2  ,:,: ] )
    plt.figure(); plt.imshow( ewfaux[ sh[0]/2-5,:,: ] )
    plt.show()

    # ============= Forces
    #forces = evalForce(atomPos, nmax_scf=100)
    #print( "Python: Forces", forces.transpose() )
    
    # Test new HS export functions
    print("\nTesting HS export functions:")
    dims = get_HS_dims()
    print("Dimensions retrieved (type):", type(dims))
    if dims.natoms > 0 and dims.norbitals > 0 : # Ensure system is initialized
        sparse_data = get_HS_sparse(dims)
        print("iatyp:", sparse_data.iatyp)
        print("h_mat shape:", sparse_data.h_mat.shape)
        print("s_mat shape:", sparse_data.s_mat.shape)

        kvec = [0.0, 0.0, 0.0]
        Hk, Sk = get_HS_k(kvec, dims.norbitals)
        print("Hk shape:", Hk.shape)
        print("Sk shape:", Sk.shape)
