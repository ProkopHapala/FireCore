
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
'int match( double* lat0, double* lat1, double Rmax, double dRmax, double dAngMax ){',
'void getMatches( int* inds, double* errs, bool bSort, double* Ks ){',
]
#cpp_utils.writeFuncInterfaces( header_strings );        exit()     #   uncomment this to re-generate C-python interfaces

cpp_utils.BUILD_PATH = os.path.normpath( cpp_utils.PACKAGE_PATH + '../../cpp/Build/libs/Molecular' ) 
lib = cpp_utils.loadLib('Lattice2D_lib', recompile=False)
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

nfound = -1

# ====================================
# ========= C functions
# ====================================

#  int walk2D( double* lat0, double* lat1, double Rmax, double dRmax ){
lib.walk2D.argtypes  = [c_double_p, c_double_p, c_double, c_double] 
lib.walk2D.restype   =  None
def walk2D(lat0, lat1, Rmax=10.0, dRmax=0.1 ):
    lat0 = np.array(lat0)
    lat1 = np.array(lat1)
    lib.walk2D(_np_as(lat0,c_double_p), _np_as(lat1,c_double_p), Rmax, dRmax )

#  int match( double* lat0, double* lat1, double Rmax, double dRmax, double dAngMax ){
lib.match.argtypes  = [c_double_p, c_double_p, c_double, c_double, c_double] 
lib.match.restype   =  c_int
def match(lat0, lat1, Rmax=10.0, dRmax=0.1, dAngMax=0.1  ):
    global nfound
    lat0 = np.array(lat0)
    lat1 = np.array(lat1)
    nfound = lib.match(_np_as(lat0,c_double_p), _np_as(lat1,c_double_p), Rmax, dRmax, dAngMax)
    return nfound

#  int getVecMatch( int nmax, int* inds, int* ns, double* errs, bool bSort, bool bV ){
lib.getVecMatch.argtypes  = [c_int, c_int_p, c_int_p, c_double_p, c_bool, c_bool, c_bool, c_bool ] 
lib.getVecMatch.restype   =  c_int
def getVecMatch(inds=None,ns=None, errs=None, bSort=True, bV=False, nmax=10, bCleans=True, bPrint=True ):
    if inds is None: inds = np.zeros( (nmax,2), dtype=np.int32)
    if ns   is None: ns   = np.zeros( (nmax),   dtype=np.int32)
    if errs is None: errs = np.zeros( (nmax)                  )
    #print( "getMatches nfound ", nfound, inds.shape, ns.shape, errs.shape )
    n=lib.getVecMatch(nmax, _np_as(inds,c_int_p), _np_as(ns,c_int_p), _np_as(errs,c_double_p), bSort, bV, bCleans, bPrint )
    return inds[:n,:],ns[:n],errs[:n]

#  void getMatches( int* inds, double* errs, bool bSort, double* Ks ){
lib.getMatches.argtypes  = [c_int_p, c_int_p, c_double_p, c_bool, c_double_p] 
lib.getMatches.restype   =  None
def getMatches(inds=None,ns=None, errs=None, bSort=True, Ks=(1.,1.,1.,0.)):
    Ks = np.array(Ks)
    if inds is None: inds = np.zeros( (nfound,4), dtype=np.int32)
    if ns   is None: ns   = np.zeros( (nfound,2), dtype=np.int32)
    if errs is None: errs = np.zeros( (nfound,4)                )
    #print( "getMatches nfound ", nfound, inds.shape, ns.shape, errs.shape )
    lib.getMatches(_np_as(inds,c_int_p), _np_as(ns,c_int_p), _np_as(errs,c_double_p), bSort, _np_as(Ks,c_double_p))
    return inds,ns, errs, Ks

# ====================================
# ========= Test Functions
# ====================================

# ====================================
# ========= Python Functions
# ====================================

def plotLattice(lat, plt, ls='k', label=None ):
    plt.plot( [lat[0,0],0.0,lat[1,0], lat[1,0]+lat[0,0], lat[0,0] ],  [lat[0,1],0.0,lat[1,1],lat[1,1]+lat[0,1], lat[0,1]], ls, label=label )

def plotSuperLattices( lat, inds, plt, n=2, ls=':' ):
    # ----- Plot best found latticle matches
    for i in range(n):
        u=lat[0,:]*inds[i,0] + lat[1,:]*inds[i,1] 
        v=lat[0,:]*inds[i,2] + lat[1,:]*inds[i,3]
        plt.plot( [u[0],0.,v[0]], [u[1],0.,v[1]], ls, label="a(%i,%i) b(%i,%i)" %(inds[i,0],inds[i,1],inds[i,2],inds[i,3])  )

# ====================================
# ========= MAIN
# ====================================

# if __name__ == "__main__":
#     import matplotlib.pyplot as plt
#     plt.show()