#!/usr/bin/python

import numpy as np
from   ctypes import c_int, c_double, c_long, c_bool, c_float, c_char_p, c_bool, c_void_p, c_char_p
import ctypes as ct
import os
import sys
#from . import utils as oclu
#import pyopencl as cl

# ========= Flags - These flags inform about state of the C/C++ library
b_Init    = False
b_initFFT = False
fft_shape = None

def tryInitFFT( sh ):
    if not b_Init: init()
    if not b_initFFT:
        fft_shape = sh
        initFFT( fft_shape )

c_float_p  = ct.POINTER(c_float)
c_double_p = ct.POINTER(c_double)
c_int_p    = ct.POINTER(c_int)
c_long_p   = ct.POINTER(c_long)

def fortReal(a):
    ct.byref(ct.c_double(a))

def _np_as(arr,atype):
    if arr is None:
        return None
    else: 
        return arr.ctypes.data_as(atype)

const_Bohr_Radius = 0.529177210903
pref_s = 0.28209479177 # sqrt(1.0f/(4.0f*M_PI));
pref_p = 0.4886025119 # sqrt(3.0f/(4.0f*M_PI));
pref_d = 1.09254843059 # sqrt(15.0f/(4.0f*M_PI));

array1ui = np.ctypeslib.ndpointer(dtype=np.uint32, ndim=1, flags='CONTIGUOUS')
array1i  = np.ctypeslib.ndpointer(dtype=np.int32,  ndim=1, flags='CONTIGUOUS')
array1l  = np.ctypeslib.ndpointer(dtype=np.int64,  ndim=1, flags='CONTIGUOUS')
array2i  = np.ctypeslib.ndpointer(dtype=np.int32,  ndim=2, flags='CONTIGUOUS')
array1d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS')
array2d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=2, flags='CONTIGUOUS')
array3d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=3, flags='CONTIGUOUS')

'''
    void init(){
    int  upload   ( int i, const float* cpu_data )
    int  download ( int i,       float* cpu_data )
    void initFFT  ( int ndim, size_t* Ns_ ){
    void run_fft( int ibuff, bool fwd, float* data )
    void convolve ( int ibuffA, int ibuffB, int ibuff_result )
    void projectAtoms( float4* atoms, float4* coefs, int ibuff_result )
    void projectAtomPosTex( float* atoms, float* coefs, int nPos, float* poss, float* out )
    void setGridShape( float* pos0, float* dA, float* dB, float* dC ){
    int  initBasisTable( int nx, int ny, float* data );
    void approx( int npoints, int npolys, double* xs, double* ys, double* ws ){
    void loadWf(const char* fname, double* out){
    void loadWfBasis( double Rcut, int nsamp, int ntmp, int nZ, int* iZs ){
    void saveToXsf(const char* fname, int ibuff){
    void initFireBall( int natoms, int* atypes, double* apos ){
    void convCoefs( int natoms, int* iZs, int* ityps, double* ocoefs, double* oatoms, int iorb0, int iorb1 ){
'''


#mode=ct.RTLD_GLOBAL
mode=ct.RTLD_LOCAL

#lib1 = ct.CDLL( name, ct.RTLD_GLOBAL )
lib2 = ct.CDLL(  "/usr/lib/x86_64-linux-gnu/libOpenCL.so", mode )
#lib2 = ct.CDLL(  "/usr/local/lib64/libclFFT.so", mode )
lib  = ct.CDLL(  "../cpp_build/libOCLfft.so",               mode )

lib.init.argtypes  = [ c_char_p ] 
lib.init.restype   =  None
def init( cl_src_dir='../cl' ):
    b_Init = True
    cl_src_dir = cl_src_dir.encode('utf8')
    lib.init( cl_src_dir )

lib.cleanup.argtypes  = [ ] 
lib.cleanup.restype   =  None
def cleanup():
    lib.cleanup( )

#void loadData( float* data_ );
#lib.upload.argtypes  = [ c_int, c_float_p ] 
#lib.upload.restype   =  c_int
#def loadData( i, data):
#    return lib.loadData( i, _np_as( data, c_float_p ) )

#void loadData( float* data_ );
lib.download.argtypes  = [ c_int, c_float_p ] 
lib.download.restype   =  c_int
def download( i, data=None, Ns=-1,dtype=np.csingle):
    if data is None:
        data = np.zeros(Ns,dtype=dtype)
    return lib.download( i, _np_as( data, c_float_p ) )

#void upload( float* data_ );
lib.upload.argtypes  = [ c_int, c_float_p ] 
lib.upload.restype   =  None
def upload( i, data):
    lib.upload( i, _np_as( data, c_float_p ) )

#void upload( float* data_ );
lib.upload_d.argtypes  = [ c_int, c_double_p, c_bool ] 
lib.upload_d.restype   =  None
def upload_d( i, data, bComplex=False ):
    lib.upload_d( i, _np_as( data, c_double_p ), bComplex )

#void initFFT ( int ndim, int* Ns );
lib.initFFT.argtypes  = [ c_int, array1l ] 
lib.initFFT.restype   =  None
def initFFT( Ns ):
    b_initFFT = True
    Ns=np.array(Ns,dtype=np.int64)
    ndim=len(Ns)
    lib.initFFT( ndim, Ns )
    print( "DEBUG py: initFFT DONE" )

#void initFFT ( int ndim, int* Ns );
lib.setErrorCheck.argtypes  = [ c_int ] 
lib.setErrorCheck.restype   =  None
def setErrorCheck( ierr ):
    lib.setErrorCheck( ierr )

#void initFFT ( int ndim, int* Ns );
lib.initAtoms.argtypes  = [ c_int ] 
lib.initAtoms.restype   =  None
def initAtoms( nAtoms, nOrbs=1 ):
    lib.initAtoms( nAtoms, nOrbs )

#int  initBasisTable( int nx, int ny, float* data );
lib.initBasisTable.argtypes  = [ c_int, c_int, c_float_p ] 
lib.initBasisTable.restype   =  c_int
def initBasisTable( nx,ny, data ):
    return lib.initBasisTable( nx,ny, _np_as( data, c_float_p ) )

lib.convolve.argtypes  = [ c_int, c_int, c_int ] 
lib.convolve.restype   =  None
def convolve(iA,iB,iOut):
    lib.convolve( iA,iB,iOut )

lib.poisson.argtypes  = [ c_int, c_int, c_float_p ] 
lib.poisson.restype   =  None
def poisson(iA,iOut, dcell):
    #if len(dcell) < 4: dcell +=  
    dcell=np.array(dcell, np.float32)
    lib.poisson( iA,iOut, _np_as( dcell, c_float_p )  )

#projectAtoms( float* atoms, float4* coefs, int ibuff_result )
lib.projectAtoms.argtypes  = [ c_float_p, c_float_p, c_int ] 
lib.projectAtoms.restype   =  None
def projectAtoms(atoms,coefs,iOut):
    lib.projectAtoms( _np_as( atoms, c_float_p ),_np_as( coefs, c_float_p ),iOut )

#projectAtomsDens( float* atoms, float* coefs, int ibuff_result, int iorb1, int iorb2 )
lib.projectAtomsDens.argtypes  = [ c_float_p, c_float_p, c_int, c_int, c_int ] 
lib.projectAtomsDens.restype   =  None
def projectAtomsDens( iOut, atoms=None,coefs=None,  iorb0=1, iorb1=2 ):
    lib.projectAtomsDens( _np_as( atoms, c_float_p ),_np_as( coefs, c_float_p ),iOut, iorb0, iorb1 )

#void projectAtomPosTex( float* atoms, float* coefs, int nPos, float* poss, float* out )
lib.projectAtomPosTex.argtypes  = [ c_float_p, c_float_p, c_int, c_float_p, c_float_p ] 
lib.projectAtomPosTex.restype   =  None
def projectAtomPosTex(atoms,coefs,poss,out=None):
    nPos=len(poss)
    if(out is None):
        out = np.zeros((nPos,2), dtype=np.float32)
    lib.projectAtomPosTex( _np_as( atoms, c_float_p ),_np_as( coefs, c_float_p ),nPos, _np_as( poss, c_float_p ), _np_as( out, c_float_p ) )
    return out

#void convCoefs( int natoms, int* iZs, int* ityps, double* ocoefs, double* oatoms, int iorb0, int iorb1 ){
lib.convCoefs.argtypes  = [ c_int, c_int_p, c_int_p, c_double_p, c_double_p, c_int, c_int, c_bool ] 
lib.convCoefs.restype   =  None
def convCoefsC( iZs, ityps, apos, wfcoefs,  iorb0=1, iorb1=2, bInit=False   ):
    natoms=len(iZs)
    iZs   = np.array(iZs,   dtype=np.int32)
    ityps = np.array(ityps, dtype=np.int32)
    lib.convCoefs( natoms, _np_as(iZs,c_int_p),_np_as(ityps,c_int_p), _np_as(wfcoefs,c_double_p), _np_as(apos,c_double_p), iorb0, iorb1, bInit )

#setGridShape( float* pos0, float* dA, float* dB, float* dC ){
lib.setGridShape.argtypes  = [ c_float_p, c_float_p, c_float_p, c_float_p ] 
lib.setGridShape.restype   =  None
def setGridShape( pos0=[0.,0.,0.,0.],dA=[0.1,0.,0.,0.],dB=[0.,0.1,0.,0.],dC=[0.,0.,0.1,0.0]):
    pos0=np.array(pos0,dtype=np.float32)
    dA=np.array(dA,dtype=np.float32)
    dB=np.array(dB,dtype=np.float32)
    dC=np.array(dC,dtype=np.float32)
    print( "dA ", dA);
    print( "dB ", dB);
    print( "dC ", dC);
    lib.setGridShape( _np_as( pos0, c_float_p ),_np_as( dA, c_float_p ), _np_as( dB, c_float_p ), _np_as( dC, c_float_p ) )

def setGridShape_dCell( Ns, dCell ):
    #pos0=[  dCell[0,0]*(-Ns[0]//2),  dCell[1,1]*(-Ns[1]//2),  dCell[2,2]*(-Ns[1]//2),  0.0 ]
    pos0=[  dCell[0,0]*(1-Ns[0]//2),  dCell[1,1]*(-Ns[1]//2),  dCell[2,2]*(-Ns[1]//2),  0.0 ]   # NOT sure why this fits best
    #print( " setGridShape_dCell() pos0 ", pos0 )
    setGridShape( 
        pos0=pos0,
        dA=dCell[0]+[0.0],
        dB=dCell[1]+[0.0],
        dC=dCell[2]+[0.0]
    )

#void initFireBall( int natoms, int* atypes, double* apos ){
lib.initFireBall.argtypes  = [ c_int, array1i, array2d ] 
lib.initFireBall.restype   =  None
def initFireBall( atypes, apos ):
    natoms = len(atypes)
    lib.initFireBall( natoms, atypes, apos )

# void run_fft( int ibuff, bool fwd, float* data )
lib.runfft.argtypes  = [ c_int, c_bool ] 
lib.runfft.restype   =  None
def runfft(ibuff, fwd=True ):
    lib.runfft( ibuff, fwd )

#void approx( int npoints, int npolys, double* xs, double* ys, double* ws ){
lib.approx.argtypes  = [ c_int, c_int, c_double_p, c_double_p, c_double_p ] 
lib.approx.restype   =  None
def approx( xs, ys, ws=None, npoly=14 ):
    n=len(xs)
    if ws is None:
        ws = np.ones(n)
    return lib.approx( n, npoly, _np_as( xs, c_double_p ),_np_as( ys, c_double_p ), _np_as( ws, c_double_p ) )

# loadWf(const char* fname, double* out){
lib.loadWf.argtypes  = [ c_char_p, c_float_p ] 
lib.loadWf.restype   =  None
def loadWf_C( fname, n=1000 ):
    data=np.zeros(n, dtype=np.float32)
    #fname = c_char_p( fname) 
    fname = fname.encode('utf-8')
    lib.loadWf( fname, _np_as( data, c_float_p ) )
    #print( data )
    return data

# void saveToXsf(const char* fname, int ibuff){
lib.saveToXsf.argtypes  = [ c_char_p, c_int ] 
lib.saveToXsf.restype   =  None
def saveToXsf( fname, ibuff ):
    fname = fname.encode('utf-8')
    lib.saveToXsf( fname, ibuff )

#void saveToXsfAtoms(const char* fname, int ibuff, int natoms, int* atypes, double* apos )
lib.saveToXsfAtoms.argtypes  = [ c_char_p, c_int,  c_int,c_int_p,c_double_p ] 
lib.saveToXsfAtoms.restype   =  None
def saveToXsfAtoms( fname, ibuff, atypes, apos ):
    fname = fname.encode('utf-8')
    na = len(atypes)
    lib.saveToXsfAtoms( fname, ibuff,   na, _np_as( atypes, c_int_p ), _np_as( apos, c_double_p ) )

#    void loadWfBasis( const char* path, float RcutSamp, int nsamp, int ntmp, int nZ, int* iZs, float* Rcuts ){
lib.loadWfBasis.argtypes  = [  c_char_p, c_float, c_int, c_int, c_int, c_int_p, c_float_p ] 
lib.loadWfBasis.restype   =  None
def loadWfBasis( iZs, nsamp=100, ntmp=1000, RcutSamp=5.0, path="Fdata/basis/", Rcuts=None, RcutDef=4.5 ):
    # NOTE : RcutSamp is in [Angstroem] while Rcuts and RcutDef are in [bohr_radius];    RcutSamp should not be changed without chaning "wf_tiles_per_angstroem" in myprog.cl

    nZ=len(iZs)
    iZs=np.array(iZs,dtype=np.int32)
    path = path.encode('utf-8')
    if Rcuts is None:
        Rcuts=np.ones(nZ,dtype=np.float32)*RcutDef
    else:
        Rcuts=np.array(Rcuts,dtype=np.float32)
    return lib.loadWfBasis( path, RcutSamp, nsamp, ntmp, nZ, _np_as( iZs, c_int_p ), _np_as( Rcuts, c_float_p ) )