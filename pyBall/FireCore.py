
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
"void getCharges( double* charges )",
"void preinit( ) ",
"void set_lvs( double* lvs )" ,
"void init( int natoms, int* atomTypes, double* atomsPos )",
"void assembleH( int iforce, int Kscf, double* positions )",
"void solveH( double* k_temp, int ikpoint )",
"void updateCharges( double sigmatol, double* sigma )",
"void SCF( int nmax_scf, double* positions, int iforce  )",
"void evalForce( int nmax_scf, double* positions_, double* forces )",
"void setupGrid( double Ecut, it ifixg0, doube* g0,  int ngrid, double* dCell  )", 
"void getGridMO( int iMO, double* ewfaux )",
"void getGridDens( int imo0, int imo1, double* ewfaux )",
"void firecore_getPointer_wfcoef( double* bbnkre )",
"void firecore_get_wfcoef( int ikp, double* wfcoefs )",
"void firecore_getpsi( int in1, int issh, int n, double x0, double dx, double ys )",
"void firecore_MOtoXSF( int iMO )",
]
#cpp_utils.writeFuncInterfaces( header_strings );        exit()     #   uncomment this to re-generate C-python interfaces

# LFLAGS   = " -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lfftw3xf_gnu -lm -Bdynamic -lmpi "

libs = {}
def loadLib(name, dic=libs):
    lib = ct.CDLL( name, ct.RTLD_GLOBAL )
    dic[name] = lib
    return lib

mkldir = "/home/prokop/SW/intel/compilers_and_libraries_2019.5.281/linux/mkl/lib/intel64_lin/"

cpp_utils.BUILD_PATH = os.path.normpath( cpp_utils.PACKAGE_PATH + '/../build/' ) 
lib = cpp_utils.loadLib('FireCore', recompile=False, mode=ct.RTLD_GLOBAL )

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


#  subroutine firecore_init( natoms_, atomTypes, atomsPos )
lib.firecore_init.argtypes  = [c_int, array1i, array2d ] 
lib.firecore_init.restype   =  None
def init(atomTypes, atomPos ):
    natoms = len(atomTypes)
    return  lib.firecore_init(natoms, atomTypes, atomPos )

#  subroutine firecore_evalForce( nmax_scf, forces_ )
lib.firecore_evalForce.argtypes  = [c_int, array2d, array2d ] 
lib.firecore_evalForce.restype   =  None
def evalForce( pos, forces=None, natom=5, nmax_scf=100 ):
    if forces is None:
        forces = np.zeros( (3,natom) )
    lib.firecore_evalForce( nmax_scf, pos, forces )
    return forces

#  void getCharges( double* charges )
lib.firecore_getCharges.argtypes  = [array2d] 
lib.firecore_getCharges.restype   =  None
def getCharges(charges):
    return lib.firecore_getCharges(charges) 


#  "void firecore_getPointer_wfcoef( double* bbnkre )
lib.firecore_getPointer_wfcoef.argtypes  = [array2d] 
lib.firecore_getPointer_wfcoef.restype   =  None
def getPointer_wfcoef(charges):
    return lib.firecore_getPointer_wfcoef(charges) 




"void firecore_get_wfcoef( int ikp, double* wfcoefs )",
#  void getCharges( double* charges )
lib.firecore_get_wfcoef.argtypes  = [c_int, array2d] 
lib.firecore_get_wfcoef.restype   =  None
def get_wfcoef(wfcoef=None,norb=None, ikp=1):
    if(wfcoef is None):
        wfcoef=np.zeros( (norb,norb) )
    lib.firecore_get_wfcoef(ikp,wfcoef) 
    return wfcoef

#  void preinit( ) 
lib.firecore_preinit.argtypes  = [] 
lib.firecore_preinit.restype   =  None
def preinit():
    return lib.firecore_preinit() 



#  void set_lvs( double* lvs )
lib.firecore_set_lvs.argtypes  = [array2d] 
lib.firecore_set_lvs.restype   =  None
def set_lvs(lvs):
    return lib.firecore_set_lvs(lvs) 

#  void assembleH( int iforce, int Kscf, double* positions )
lib.firecore_assembleH.argtypes  = [c_int, c_int, array2d ] 
lib.firecore_assembleH.restype   =  None
def assembleH( positions, iforce=0, Kscf=1 ):
    return lib.firecore_assembleH(iforce, Kscf, positions) 



#  void solveH( double* k_temp, int ikpoint )
lib.firecore_solveH.argtypes  = [ array1d, c_int] 
lib.firecore_solveH.restype   =  None
def solveH(k_temp=None, ikpoint=1):
    if k_temp is None:
        k_temp = np.array([0.0,0.0,0.0])
    return lib.firecore_solveH(k_temp, ikpoint ) 



#  void updateCharges( double sigmatol, double* sigma )
lib.firecore_updateCharges.argtypes  = [c_double, c_double_p] 
lib.firecore_updateCharges.restype   =  None
def updateCharges( sigmatol=1e-6, sigma=None):
    if sigma is None:
        sigma = np.zeros(1)
    lib.firecore_updateCharges(sigmatol, _np_as(sigma,c_double_p)) 
    return sigma


#  void SCF( int nmax_scf, double* positions, int iforce  )
lib.firecore_SCF.argtypes  = [c_int, array2d, c_int] 
lib.firecore_SCF.restype   =  None
def SCF( positions, iforce=0, nmax_scf=200 ):
    return lib.firecore_SCF( nmax_scf, positions, iforce ) 


#  void setupGrid( double Ecut, int ifixg0, doube* g0,  int ngrid, double* dCell  )
lib.firecore_setupGrid.argtypes  = [c_double, c_int, array1d, array1i, array2d ] 
lib.firecore_setupGrid.restype   =  None
def setupGrid(Ecut=100, ifixg0=0, g0=None, ngrid=None, dCell=None):
    if g0 is None:
        g0=np.array([0.0,0.0,0.0])
    if ngrid is None:
        ngrid = np.zeros(3,dtype=np.int32)
    if dCell is None:
        dCell = np.zeros( (3,3) )
    lib.firecore_setupGrid(Ecut, ifixg0, g0, ngrid, dCell )
    lvs = np.array( [ [0.0,0.0,0.0], dCell[0]*ngrid[0], dCell[1]*ngrid[1], dCell[2]*ngrid[2], ] )
    ngrid = ngrid[::-1].copy()
    return ngrid, dCell, lvs



#  void getGridMO( int iMO, double* ewfaux )
lib.firecore_getGridMO.argtypes  = [ c_int, array3d ] 
lib.firecore_getGridMO.restype   =  None
def getGridMO(iMO, ewfaux=None, ngrid=None):
    if ewfaux is None:
        ewfaux = np.zeros(ngrid)
    lib.firecore_getGridMO(iMO, ewfaux )
    return ewfaux

#  void firecore_MOtoXSF( int iMO )
lib.firecore_orb2xsf.argtypes  = [ c_int ] 
lib.firecore_orb2xsf.restype   =  None
def orb2xsf(iMO, ewfaux=None, ngrid=None):
    if ewfaux is None:
        ewfaux = np.zeros(ngrid)
    lib.firecore_orb2xsf(iMO )
    return ewfaux

#  void getGridDens( int imo0, int imo1, double* ewfaux )
lib.firecore_getGridDens.argtypes  = [c_int, c_int, array3d ] 
lib.firecore_getGridDens.restype   =  None
def getGridDens(imo0=0, imo1=1000, ewfaux=None, ngrid=None):
    #ngrid=ngrid[::-1]
    print( ngrid ); #exit(0)
    if ewfaux is None:
        ewfaux = np.zeros(ngrid)
    lib.firecore_getGridDens(imo0, imo1, ewfaux ) 
    return ewfaux 


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
lib.firecore_getpsi.argtypes  = [c_int, c_int, c_int, c_int,    c_int, array2d, array1d ] 
lib.firecore_getpsi.restype   =  None
def getpsi( poss, ys=None, in1=1, issh=1, l=0, m=1 ):
    n = len(poss)
    if ys is None:
        ys = np.zeros(n)
    lib.firecore_getpsi( l, m, in1, issh, n, poss, ys )
    return ys

# ========= Python Functions

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
    

