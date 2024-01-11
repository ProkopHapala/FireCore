#!/usr/bin/python

import numpy as np
from   ctypes import c_int, c_double, c_long, c_bool, c_float, c_char_p, c_bool, c_void_p, c_char_p
import ctypes as ct
import os
import sys
#from . import utils as oclu
#import pyopencl as cl

# ========= Flags - These flags inform about state of the C/C++ library
b_Init        = False
b_initFFT     = False
b_basisLoaded = False
fft_shape = None

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

#mode=ct.RTLD_LAZY
#mode=ct.RTLD_NOW
mode=ct.RTLD_LOCAL

#lib1 = ct.CDLL( name, ct.RTLD_GLOBAL )
lib2 = ct.CDLL(  "/usr/lib/x86_64-linux-gnu/libOpenCL.so", mode )
#lib2 = ct.CDLL(  "/usr/local/lib64/libclFFT.so", mode )

#bUseLocal = True
bUseLocal = False
if bUseLocal:
    lib  = ct.CDLL(  "../cpp_build/libOCLfft.so",                      mode )
else:
    # === NOTE: if we don't compile in delendencies to libOCL_GridFF.so we need to load them before with RTLD_GLOBAL
    #mode=ct.RTLD_GLOBAL
    #lib1 = ct.CDLL(  "/usr/lib/x86_64-linux-gnu/libOpenCL.so",  mode )
    #lib2 = ct.CDLL(  "/usr/lib/x86_64-linux-gnu/libclFFT.so",   mode )
    #lib  = ct.CDLL(  "../../cpp/Build_OCL/libs_OCL/libOCL_GridFF.so",  mode )
    lib  = ct.CDLL(  "../../cpp/Build/libs_OCL/libOCL_GridFF.so",  mode )


lib.init.argtypes  = [ c_char_p ] 
lib.init.restype   =  None
def init( cl_src_dir='./data/cl' ):
    global b_Init; b_Init = True
    cl_src_dir = cl_src_dir.encode('utf8')
    lib.init( cl_src_dir )

# void release( bool bReleaseOCL, bool bReleaseOCLfft);
lib.release.argtypes  = [ c_bool, c_bool ]
lib.release.restype   =  None
def release( bReleaseOCL=False, bReleaseOCLfft=True ):
    lib.release( bReleaseOCL, bReleaseOCLfft )

# void printDeviceInfo( bool bDetails )
lib.printDeviceInfo.argtypes  = [ c_bool ]
lib.printDeviceInfo.restype   =  None
def printDeviceInfo( bDetails=False ):
    lib.printDeviceInfo( bDetails )

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
        print( "download() allocated %i bytes ", data.nbytes )
    lib.download( i, _np_as( data, c_float_p ) )
    return data

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

# int copy( int iBufFrom, int iBufTo )
lib.copy.argtypes  = [ c_int, c_int,  c_int,  c_int,  c_int ] 
lib.copy.restype   =  c_int
def copy( iBufFrom, iBufTo, nbytes=-1, src_offset=0, dst_offset=0   ):
    return lib.copy( iBufFrom, iBufTo, nbytes, src_offset, dst_offset )

#int copyBuffToImage( int iBuff, int itex, int nbytes=-1, int src_offset=0                   ){
lib.copyBuffToImage.argtypes  = [ c_int, c_int,  c_int,  c_int, c_int ] 
lib.copyBuffToImage.restype   =  c_int
def copyBuffToImage( iBuff, itex, nx,ny,nz ):
    return lib.copyBuffToImage( iBuff, itex, nx,ny,nz )

# void roll_buf( int ibuffA, int ibuffB, int* shift )
lib.roll_buf.argtypes  = [ c_int, c_int, c_int_p ] 
lib.roll_buf.restype   =  None
def roll( iBufFrom, iBufTo, shift ):
    shift = np.array( shift, np.int32 )
    return lib.roll_buf( iBufFrom, iBufTo,  _np_as( shift, c_int_p )  )

#void initFFT ( int ndim, int* Ns );
lib.initFFT.argtypes  = [ c_int, array1l ] 
lib.initFFT.restype   =  None
def initFFT( Ns ):
    global b_initFFT; b_initFFT = True
    Ns=np.array(Ns,dtype=np.int64)
    ndim=len(Ns)
    lib.initFFT( ndim, Ns )

#void initFFT ( int ndim, int* Ns );
lib.setErrorCheck.argtypes  = [ c_int ] 
lib.setErrorCheck.restype   =  None
def setErrorCheck( ierr ):
    lib.setErrorCheck( ierr )

#void setVerbosity(int verbosity_ )
lib.setVerbosity.argtypes  = [ c_int ] 
lib.setVerbosity.restype   =  None
def setVerbosity( verbosity ):
    lib.setVerbosity( verbosity )

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
def poisson( iA, iOut, dcell):
    #if len(dcell) < 4: dcell +=  
    dcell=np.array(dcell, np.float32)
    lib.poisson( iA,iOut, _np_as( dcell, c_float_p )  )

lib.gradient.argtypes  = [ c_int, c_int, c_float_p ] 
lib.gradient.restype   =  None
def gradient(iA,iOut, dcell):
    #if len(dcell) < 4: dcell +=  
    dcell=np.array(dcell, np.float32)
    lib.gradient( iA,iOut, _np_as( dcell, c_float_p )  )

# void evalVpointChargesPBC( int na, double* apos, double* aQs, int np, double* ps, double* Vps, int* nPBC, double* cell ){
lib.evalVpointChargesPBC.argtypes  = [ c_int, c_double_p, c_double_p, c_int, c_double_p, c_double_p, c_int_p, c_double_p ]
lib.evalVpointChargesPBC.restype   =  None
def evalVpointChargesPBC( apos, aQs, ps, cell, Vps=None, nPBC=[0,0,0] ):
    na = len(apos)
    n = len(ps)
    apos = np.array(apos, dtype=np.float64)
    aQs  = np.array(aQs,  dtype=np.float64)
    ps   = np.array(ps,   dtype=np.float64)
    cell = np.array(cell, dtype=np.float64)
    if Vps is None:
        Vps = np.zeros(n, dtype=np.float64)
    nPBC = np.array(nPBC, dtype=np.int32)
    lib.evalVpointChargesPBC( na, _np_as(apos,c_double_p), _np_as(aQs,c_double_p), n, _np_as(ps,c_double_p), _np_as(Vps,c_double_p), _np_as(nPBC,c_int_p), _np_as(cell,c_double_p) )
    return Vps

#projectAtoms( float* atoms, float4* coefs, int ibuff_result )
lib.projectAtoms.argtypes  = [ c_float_p, c_float_p, c_int ] 
lib.projectAtoms.restype   =  None
def projectAtoms(atoms,coefs,iOut):
    lib.projectAtoms( _np_as( atoms, c_float_p ),_np_as( coefs, c_float_p ),iOut )

#projectAtomsDens( float* atoms, float* coefs, int ibuff_result, int iorb1, int iorb2 )
lib.projectAtomsDens.argtypes  = [ c_float_p, c_float_p, c_int, c_int, c_int, c_float_p ] 
lib.projectAtomsDens.restype   =  None
def projectAtomsDens( iOut, atoms=None,coefs=None,  iorb0=0, iorb1=1, acumCoef=[0.0,1.0] ):
    acumCoef = np.array(acumCoef,dtype=np.float32)
    lib.projectAtomsDens( _np_as( atoms, c_float_p ),_np_as( coefs, c_float_p ),iOut, iorb0, iorb1, _np_as( acumCoef, c_float_p ) )

#void projectAtomsDens0( int ibuff_result, float* acumCoef, int natoms=0, int* ityps=0, Vec3d* oatoms=0 )
lib.projectAtomsDens0.argtypes  = [ c_int, c_float_p, c_int, c_int_p, c_double_p, c_float_p ] 
lib.projectAtomsDens0.restype   =  None
def projectAtomsDens0( iOut, atypes=None, apos=None, acumCoef=[0.0,1.0], coefs=None ):
    natom=0
    if atypes is not None: natom=len(atypes)
    acumCoef = np.array(acumCoef,dtype=np.float32)
    #print( "type(atypes)    ", type(atypes)   )
    #print( "type(apos)      ", type(apos)     )
    #print( "type(acumCoef)  ", type(acumCoef) )
    lib.projectAtomsDens0( iOut, _np_as( acumCoef, c_float_p ), natom, _np_as( atypes, c_int_p ), _np_as( apos, c_double_p ), _np_as( coefs, c_float_p ) )

#void projectAtomPosTex( float* atoms, float* coefs, int nPos, float* poss, float* out )
lib.projectAtomPosTex.argtypes  = [ c_float_p, c_float_p, c_int, c_float_p, c_float_p ] 
lib.projectAtomPosTex.restype   =  None
def projectAtomPosTex(atoms,coefs,poss,out=None):
    nPos=len(poss)
    if(out is None):
        out = np.zeros((nPos,2), dtype=np.float32)
    lib.projectAtomPosTex( _np_as( atoms, c_float_p ),_np_as( coefs, c_float_p ),nPos, _np_as( poss, c_float_p ), _np_as( out, c_float_p ) )
    return out

# projectDenmat( int natoms, int* iZs, int* ityps, double* ocoefs, double* apos, int iorb0, int iorb1, double Rcut, bool bInit ){ 
lib.projectDenmat.argtypes  = [ c_int, c_int_p, c_int_p, c_double_p, c_double_p, c_int, c_int, c_double, c_bool ]
lib.projectDenmat.restype   =  None
def projectDenmat( iZs, ityps, apos, iorb0=0, iorb1=1, Rcut=5.0, bInit=False ):
    natoms=len(iZs)
    iZs   = np.array(iZs,   dtype=np.int32)
    ityps = np.array(ityps, dtype=np.int32)-1
    lib.projectDenmat( natoms, _np_as(iZs,c_int_p),_np_as(ityps,c_int_p), _np_as(apos,c_double_p), iorb0, iorb1, Rcut, bInit )

#void setTypes( int* atype_nOrb_, float* atype_Qconfs_ ){
lib.setTypes.argtypes  = [ c_int, c_int_p, c_float_p, c_bool ] 
lib.setTypes.restype   =  None
def setTypes( atype_nOrb, atype_Qconfs, bInternal=True ):
    atype_nOrb  =np.array(atype_nOrb, dtype=np.int32)
    atype_Qconfs=np.array(atype_Qconfs, dtype=np.float32)
    ntyp = len(atype_nOrb)
    lib.setTypes( ntyp, _np_as(atype_nOrb,c_int_p),_np_as(atype_Qconfs,c_float_p), bInternal )

#void convCoefs( int natoms, int* iZs, int* ityps, double* ocoefs, double* oatoms, int iorb0, int iorb1 ){
lib.convCoefs.argtypes  = [ c_int, c_int_p, c_int_p, c_double_p, c_double_p, c_bool, c_bool ] 
lib.convCoefs.restype   =  c_int
def convCoefsC( iZs, ityps, apos, wfcoefs, bInit=False, bDiagonal=False ):
    natoms=len(iZs)
    iZs   = np.array(iZs,   dtype=np.int32)
    ityps = np.array(ityps, dtype=np.int32)-1
    return lib.convCoefs( natoms, _np_as(iZs,c_int_p),_np_as(ityps,c_int_p), _np_as(wfcoefs,c_double_p), _np_as(apos,c_double_p), bInit, bDiagonal )

#setGridShape( float* pos0, float* dA, float* dB, float* dC ){
lib.setGridShape.argtypes  = [ c_float_p, c_float_p, c_float_p, c_float_p ] 
lib.setGridShape.restype   =  None
def setGridShape( pos0=[0.,0.,0.,0.],dA=[0.1,0.,0.,0.],dB=[0.,0.1,0.,0.],dC=[0.,0.,0.1,0.0]):
    pos0=np.array(pos0,dtype=np.float32)
    dA=np.array(dA,dtype=np.float32)
    dB=np.array(dB,dtype=np.float32)
    dC=np.array(dC,dtype=np.float32)
    #print( "dA ", dA); print( "dB ", dB); print( "dC ", dC);
    lib.setGridShape( _np_as( pos0, c_float_p ),_np_as( dA, c_float_p ), _np_as( dB, c_float_p ), _np_as( dC, c_float_p ) )

def getCellHalf( Ns, dCell ):
    return [  dCell[0,0]*(1-Ns[0]//2),  dCell[1,1]*(-Ns[1]//2),  dCell[2,2]*(-Ns[2]//2),  0.0 ]

def setGridShape_dCell( Ns, dCell, pos0=None ):
    if pos0 is None:
        pos0 = getCellHalf( Ns, dCell )  # NOT sure why this fits best
    setGridShape( pos0=pos0, dA=dCell[0]+[0.0], dB=dCell[1]+[0.0], dC=dCell[2]+[0.0] )

# void run_fft( int ibuff, bool fwd, float* data )
lib.runfft.argtypes  = [ c_int, c_bool ] 
lib.runfft.restype   =  None
def runfft(ibuff, fwd=True ):
    lib.runfft( ibuff, fwd )

#int newFFTbuffer( char* name ){ oclfft.newFFTbuffer( name ); }
lib.newFFTbuffer.argtypes  = [ c_char_p, c_int, c_int ] 
lib.newFFTbuffer.restype   =  c_int
def newFFTbuffer( fname, nfloat=2, ntot=-1 ):
    fname = fname.encode('utf-8')
    return lib.newFFTbuffer( fname, nfloat, ntot )

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

#void saveToBin(const char* fname, int ibuff){
lib.saveToBin.argtypes  = [ c_char_p, c_int ] 
lib.saveToBin.restype   =  None
def saveToBin( fname, ibuff ):
    fname = fname.encode('utf-8')
    lib.saveToBin( fname, ibuff )

#void loadFromBin(const char* fname, int ibuff){
lib.loadFromBin.argtypes  = [ c_char_p, c_int ] 
lib.loadFromBin.restype   =  None
def loadFromBin( fname, ibuff ):
    fname = fname.encode('utf-8')
    #print( "DEBUG loadFromBin | fname ", fname  )
    #print( "DEBUG loadFromBin | ibuff ", ibuff  )
    lib.loadFromBin( fname, ibuff )

# void saveToXsf(const char* fname, int ibuff){
lib.saveToXsf.argtypes  = [ c_char_p, c_int, c_int, c_int ] 
lib.saveToXsf.restype   =  None
def saveToXsf( fname, ibuff, stride=2, offset=0 ):
    fname = fname.encode('utf-8')
    lib.saveToXsf( fname, ibuff, stride, offset )

#void saveToXsfAtoms(const char* fname, int ibuff, int natoms, int* atypes, double* apos )
lib.saveToXsfAtoms.argtypes  = [ c_char_p, c_int, c_int, c_int, c_int,c_int_p,c_double_p ] 
lib.saveToXsfAtoms.restype   =  None
def saveToXsfAtoms( fname, ibuff, atypes, apos, stride=2, offset=0 ):
    fname = fname.encode('utf-8')
    na = len(atypes)
    lib.saveToXsfAtoms( fname, ibuff, stride, offset,  na, _np_as( atypes, c_int_p ), _np_as( apos, c_double_p ) )

#void saveToXsfAtomsData(const char* fname, double* data, int natoms, int* atypes, double* apos )
lib.saveToXsfAtomsData.argtypes  = [ c_char_p, c_int_p, c_double_p, c_int,c_int_p,c_double_p ] 
lib.saveToXsfAtomsData.restype   =  None
def saveToXsfAtomsData( fname, data, atypes, apos ):
    fname = fname.encode('utf-8')
    na = len(atypes)
    ngrid = np.array( data.shape, dtype=np.int32 )
    lib.saveToXsfAtomsData( fname, _np_as(ngrid, c_int_p ),  _np_as( data, c_double_p ),  na, _np_as( atypes, c_int_p ), _np_as( apos, c_double_p ) )

#    void loadWfBasis( const char* path, float RcutSamp, int nsamp, int ntmp, int nZ, int* iZs, float* Rcuts ){
lib.loadWfBasis.argtypes  = [  c_char_p, c_float, c_int, c_int, c_int, c_int_p, c_float_p, c_bool ] 
lib.loadWfBasis.restype   = c_float_p
def loadWfBasis( iZs, nsamp=100, ntmp=1000, RcutSamp=5.0, path="Fdata/basis/", Rcuts=None, RcutDef=4.5, bDelete=True ):
    # NOTE : RcutSamp is in [Angstroem] while Rcuts and RcutDef are in [bohr_radius];    RcutSamp should not be changed without chaning "wf_tiles_per_angstroem" in myprog.cl
    global b_basisLoaded; b_basisLoaded=True; 
    nZ=len(iZs)
    iZs=np.array(iZs,dtype=np.int32)
    path = path.encode('utf-8')
    if Rcuts is None:
        Rcuts=np.ones(nZ,dtype=np.float32)*RcutDef
    else:
        Rcuts=np.array(Rcuts,dtype=np.float32)
    wf_data = lib.loadWfBasis( path, RcutSamp, nsamp, ntmp, nZ, _np_as( iZs, c_int_p ), _np_as( Rcuts, c_float_p ), bDelete )
    # convert wf_data from c_float_p to np.array, if it is not null pointer
    if wf_data: # check if wf_data is not null pointer
        print( "loadWfBasis ", wf_data )
        wf_data = np.ctypeslib.as_array(wf_data, shape=(nZ,nsamp,2))
    #else:
    #    print( "loadWfBasis  bDelete, wf_data ", bDelete, wf_data )
    return wf_data

# =================== OCL_PP

#void initPP( const char* cl_src_dir, size_t* Ns_ ){
lib.initPP.argtypes  = [ c_char_p , c_long_p ] 
lib.initPP.restype   =  c_int
def initPP( Ns, path='./data/cl' ):
    global b_Init, b_initFFT
    b_Init    = True
    b_initFFT = True
    path = path.encode('utf-8')
    Ns=np.array(Ns,dtype=np.int64)
    return lib.initPP(path, _np_as(Ns,c_long_p) )

#void makeStartPointGrid( int nx, int ny, double* p0, double* da, double* db )
lib.makeStartPointGrid.argtypes  = [ c_int,c_int, array1d,array1d,array1d ] 
lib.makeStartPointGrid.restype   =  None
def makeStartPointGrid( nx,ny, p0, da, db ):
    p0=np.array(p0)
    da=np.array(da)
    db=np.array(db)
    lib.makeStartPointGrid( nx,ny, p0,da,db )

#void setGridShapePP    ( double* p0, double* dCell                          )
lib.setGridShapePP.argtypes  = [ c_double_p, array2d ] 
lib.setGridShapePP.restype   =  None
def setGridShapePP( dCell, p0=None ):
    dCell=np.array(dCell)
    if p0 is not None:
        p0=np.array(p0)
    lib.setGridShapePP( _np_as(p0,c_double_p), dCell )

#void relaxStrokesTilted( int ibuff_out, int np=0, float* points=0 )
lib.relaxStrokesTilted.argtypes  = [ c_int, c_int, c_float, c_int, c_float_p ] 
lib.relaxStrokesTilted.restype   =  None
def relaxStrokesTilted( iBuffOut, nz=20, dtip=-0.1 ):
    lib.relaxStrokesTilted( iBuffOut, nz, dtip, 0,None )

#void getFEinStrokes    ( int ibuff_out, int np=0, float* points=0 )
lib.getFEinStrokes.argtypes  = [ c_int,  c_int, array1d,  c_int, c_float_p ] 
lib.getFEinStrokes.restype   =  None
def getFEinStrokes( iBuffOut, nz=10, dTip=[0.0,0.0,-0.1] ):
    dTip = np.array(dTip)
    lib.getFEinStrokes( iBuffOut, nz, dTip,    0, None )

#void evalLJC_QZs( int ibuff_out, int na=0, float* atoms=0, float* coefs=0 ){
lib.evalLJC_QZs.argtypes  = [ c_int, c_int, c_float_p, c_float_p ] 
lib.evalLJC_QZs.restype   =  None
def evalLJC_QZs( iBuffOut, apos, cLJs, Qs=None ):
    na = len(apos)
    atoms = np.zeros( (na,4), np.float32 ); atoms[:,:3] = apos[:,:]; atoms[:,3]=Qs;
    coefs = np.zeros( (na,2), np.float32 ); coefs[:,:]  = cLJs[:,:]
    lib.evalLJC_QZs( iBuffOut, na, _np_as(atoms,c_float_p), _np_as(coefs,c_float_p) )

#void evalLJC_QZs_toImg(  int na=0, float* atoms=0, float* coefs=0 ){
lib.evalLJC_QZs_toImg.argtypes  = [ c_int, c_float_p, c_float_p ] 
lib.evalLJC_QZs_toImg.restype   =  None
def evalLJC_QZs_toImg( apos, cLJs, Qs=None ):
    na = len(apos)
    atoms = np.zeros( (na,4), np.float32 ); atoms[:,:3] = apos[:,:]; atoms[:,3]=Qs;
    coefs = np.zeros( (na,2), np.float32 ); coefs[:,:]  = cLJs[:,:]
    lib.evalLJC_QZs_toImg( na, _np_as(atoms,c_float_p), _np_as(coefs,c_float_p) )

'''
# ===================== COMMENTED / DEPRECATED FUNCTIONS

#void approx( int npoints, int npolys, double* xs, double* ys, double* ws ){
lib.approx.argtypes  = [ c_int, c_int, c_double_p, c_double_p, c_double_p ] 
lib.approx.restype   =  None
def approx( xs, ys, ws=None, npoly=14 ):
    n=len(xs)
    if ws is None:
        ws = np.ones(n)
    return lib.approx( n, npoly, _np_as( xs, c_double_p ),_np_as( ys, c_double_p ), _np_as( ws, c_double_p ) )

#void initFireBall( int natoms, int* atypes, double* apos ){
lib.initFireBall.argtypes  = [ c_int, array1i, array2d ] 
lib.initFireBall.restype   =  None
def initFireBall( atypes, apos ):
    natoms = len(atypes)
    lib.initFireBall( natoms, atypes, apos )
'''

# ===================== PYTHON

def saveBuff( iBuff, fname ):
    if   fname[-4:] == ".xsf":
        saveToXsf( fname, iBuff )
    elif fname[-4:] == ".bin":  
        saveToBin( fname, iBuff )

def tryInitFFT( sh ):
    #print( "!!!!------------------------- tryInitFFT: b_Init, b_initFFT ", b_Init, b_initFFT )
    if not b_Init: 
        #print( "!!!!------------------------- tryInitFFT: init() " )
        init()
    if not b_initFFT:
        #print( "!!!!------------------------- initFFT: init() " )
        fft_shape = sh
        initFFT( fft_shape )

def tryLoadWfBasis( iZs, nsamp=100, ntmp=1000, RcutSamp=5.0, path="Fdata/basis/", Rcuts=None, RcutDef=4.5, bDelete=True ):
    if not b_basisLoaded:
        return loadWfBasis( iZs, nsamp=nsamp, ntmp=ntmp, RcutSamp=RcutSamp, path=path, Rcuts=Rcuts, RcutDef=RcutDef, bDelete=bDelete )

def initFFTgrid( Ns, dcell = [0.2,0.2,0.2,0.2], dCell=None ):
    init() 
    initFFT( Ns  )
    if dCell is None:
        dCell = np.array([[dcell[0],0.,0.],[0.,dcell[1],0.],[0.,0.,dcell[2]]])
    setGridShape_dCell( Ns, dCell )

default_AtomTypeConfs=[ 
#  norb,   Qshell,   Name
( 1, (1.0,0.0), "H"  ),
( 4, (2.0,0.0), "He" ),
( 4, (1.0,0.0), "Li" ),
( 4, (1.0,1.0), "Be" ),
( 4, (1.0,2.0), "B"  ),
( 4, (1.0,3.0), "C"  ), 
( 4, (2.0,3.0), "N"  ), 
( 4, (2.0,4.0), "O"  ),
( 4, (2.0,5.0), "F"  ), 
]

def setTypesZ( iZs, bInternal=True, AtomTypeConfs=default_AtomTypeConfs ):
    atype_nOrb   = []
    atype_Qconfs = []
    if AtomTypeConfs is None:
        AtomTypeConfs=default_AtomTypeConfs
    for iZ in iZs:
        typ=AtomTypeConfs[iZ-1]
        print( "iZ, typ ", iZ, typ )
        atype_nOrb  .append(typ[0])
        atype_Qconfs.append(typ[1])
    return setTypes( atype_nOrb, atype_Qconfs )