
from ast import arg
import numpy as np
from   ctypes import c_int, c_double, c_bool, c_float, c_char_p, c_bool, c_void_p, c_char_p
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
    #ngrid=ngrid[::-1]
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


cpp_utils.set_args_dict(lib, argDict)

#===================================
# ========= Python Functions
#===================================

def run_nonSCF( atomType, atomPos ):
    preinit()
    norb = init( atomType, atomPos )
    assembleH( atomPos)
    solveH()
    sigma=updateCharges(); #print( sigma )
    return norb, sigma

def initialize( atomType=None, atomPos=None, verbosity=0 ):
    global norb
    setVerbosity(verbosity=verbosity)
    preinit()
    norb = init( atomType, atomPos )
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
    

