
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
#"void init_buffers()",
#"void init_params(char* fatomtypes, char* fbondtypes)",
#"void init_nonbond()",
#"void buildFF( bool bNonBonded_, bool bOptimizer_ )",
#"int loadmol( char* fname_mol )",
#"void initWithMolFile(char* fname_mol, bool bNonBonded_, bool bOptimizer_ )",
#"void setSwitches( int doAngles, int doPiPiT, int  doPiSigma, int doPiPiI, int doBonded_, int PBC, int CheckInvariants )",
#"bool checkInvariants( double maxVcog, double maxFcog, double maxTg )",
#"double eval()",
#"bool relax( int niter, double Ftol )",
"void shift_atoms_ax( int n, int* selection, double* d                  )",
"void shift_atoms   ( int n, int* selection, int ia0, int ia1, double l )",
"void rotate_atoms_ax( int n, int* selection, double* p0, double* ax, double phi      )",
"void rotate_atoms   ( int n, int* selection, int ia0, int iax0, int iax1, double phi )",
"int  splitAtBond( int ib, int* selection )",
"void scanTranslation_ax( int n, int* selection, double* vec, int nstep, double* Es, bool bWriteTrj )",
"void scanTranslation( int n, int* selection, int ia0, int ia1, double l, int nstep, double* Es, bool bWriteTrj )",
#"void scanRotation_ax( int n, int* selection, double* p0, double* ax, double phi, int nstep, double* Es, bool bWriteTrj )",
#"void scanRotation( int n, int* selection,int ia0, int iax0, int iax1, double phi, int nstep, double* Es, bool bWriteTrj )",
]
#cpp_utils.writeFuncInterfaces( header_strings );        exit()     #   uncomment this to re-generate C-python interfaces

#libSDL = ctypes.CDLL( "/usr/lib/x86_64-linux-gnu/libSDL2.so", ctypes.RTLD_GLOBAL )
#libGL  = ctypes.CDLL( "/usr/lib/x86_64-linux-gnu/libGL.so",   ctypes.RTLD_GLOBAL )


#cpp_name='CombatModels'
#cpp_utils.make(cpp_name)
#LIB_PATH      = os.path.dirname( os.path.realpath(__file__) )
#LIB_PATH_CPP  = os.path.normpath(LIB_PATH+'../../../'+'/cpp/Build/libs/'+cpp_name )
#lib = ctypes.CDLL( LIB_PATH_CPP+("/lib%s.so" %cpp_name) )

cpp_utils.BUILD_PATH = os.path.normpath( cpp_utils.PACKAGE_PATH + '../../cpp/Build/libs/Molecular' ) 
lib = cpp_utils.loadLib('MMFF_lib', recompile=False)
array1ui = np.ctypeslib.ndpointer(dtype=np.uint32, ndim=1, flags='CONTIGUOUS')
array1i  = np.ctypeslib.ndpointer(dtype=np.int32,  ndim=1, flags='CONTIGUOUS')
array2i  = np.ctypeslib.ndpointer(dtype=np.int32,  ndim=2, flags='CONTIGUOUS')
array1d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS')
array2d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=2, flags='CONTIGUOUS')
array3d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=3, flags='CONTIGUOUS')

# ====================================
# ========= C functions
# ====================================

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
    return np.ctypeslib.as_array( ptr, shape=sh)


#double* getBuff(const char* name){ 
lib.getBuff.argtypes = [c_char_p]
lib.getBuff.restype  = c_double_p 
def getBuff(name,sh):
    if not isinstance(sh, tuple): sh=(sh,)
    name=name.encode('utf8')
    ptr = lib.getBuff(name)
    return np.ctypeslib.as_array( ptr, shape=sh)

#def getBuffs( nnode, npi, ncap, nbond, NEIGH_MAX=4 ):
def getBuffs( NEIGH_MAX=4 ):
    #natom=nnode+ncap
    #nvecs=natom+npi
    #nDOFs=nvecs*3
    global ndims,Es
    ndims = getIBuff( "ndims", (6,) )  # [nDOFs,natoms,nnode,ncap,npi,nbonds]
    Es    = getBuff ( "Es",    (6,) )  # [ Etot,Eb,Ea, Eps,EppT,EppI; ]
    global nDOFs,natoms,nnode,ncap,nbonds,npi,nvecs
    nDOFs=ndims[0]; natoms=ndims[1]; nnode=ndims[2];ncap=ndims[3];npi=ndims[4];nbonds=ndims[5];
    nvecs=natoms+npi
    print( "getBuffs(): nbonds %i nvecs %i npi %i natoms %i nnode %i ncap %i" %(nbonds,nvecs,npi,natoms,nnode,ncap) )
    global DOFs,fDOFs,apos,fapos,pipos,fpipos,bond_l0,bond_k, Kneighs,bond2atom,aneighs,selection
    #Ebuf     = getEnergyTerms( )
    DOFs      = getBuff ( "DOFs",     (nvecs,3)  )
    fDOFs     = getBuff ( "fDOFs",    (nvecs,3)  ) 
    apos      = getBuff ( "apos",     (natoms,3) )
    fapos     = getBuff ( "fapos",    (natoms,3) )
    pipos     = getBuff ( "pipos",    (npi,3)    )
    fpipos    = getBuff ( "fpipos",   (npi,3)    )
    bond_l0   = getBuff ( "bond_l0",  (nbonds)   )
    bond_k    = getBuff ( "bond_k",   (nbonds)   )
    Kneighs   = getBuff ( "Kneighs",  (nnode,NEIGH_MAX) )
    bond2atom = getIBuff( "bond2atom",(nbonds,2) )
    aneighs   = getIBuff( "aneighs",  (nnode,NEIGH_MAX) )
    selection = getIBuff( "selection",  (natoms) )

#  void init_buffers()
lib.init_buffers.argtypes  = []
lib.init_buffers.restype   =  None
def init_buffers():
    return lib.init_buffers()

#  void init_params(const char* fatomtypes, const char* fbondtypes)
lib.init_params.argtypes  = [c_char_p, c_char_p, c_char_p] 
lib.init_params.restype   =  None
def init_params(fatomtypes, fbondtypes, fbondangles ):
    fatomtypes = fatomtypes.encode('utf8')
    fbondtypes = fbondtypes.encode('utf8')
    fbondangles = fbondangles.encode('utf8')
    return lib.init_params(fatomtypes,fbondtypes,fbondangles)

#  void init_nonbond()
lib.init_nonbond.argtypes  = [] 
lib.init_nonbond.restype   =  None
def init_nonbond():
    return lib.init_nonbond()

#  void insertSMILES(char* s)
lib.insertSMILES.argtypes  = [c_char_p,c_bool] 
lib.insertSMILES.restype   =  None
def insertSMILES(s, bPrint=False, bCap=False ):
    s = s.encode('utf8')
    return lib.insertSMILES(s, bPrint, bCap )

#  void buildFF( bool bNonBonded_, bool bOptimizer_ )
lib.buildFF.argtypes  = [c_bool, c_bool] 
lib.buildFF.restype   =  None
def buildFF(bNonBonded=True, bOptimizer=True):
    return lib.buildFF(bNonBonded, bOptimizer)

#  int loadmol( const char* fname_mol )
lib.loadmol.argtypes  = [c_char_p] 
lib.loadmol.restype   =  c_int
def loadmol(fname_mol):
    fname_mol=fname_mol.encode('utf8')
    return lib.loadmol(fname_mol)

#  void initWithMolFile(char* fname_mol, bool bNonBonded_, bool bOptimizer_ )
lib.initWithMolFile.argtypes  = [c_char_p, c_bool, c_bool] 
lib.initWithMolFile.restype   =  None
def initWithMolFile(fname_mol, bNonBonded=True, bOptimizer=True):
    fname_mol = fname_mol.encode('utf8')
    return lib.initWithMolFile(fname_mol, bNonBonded, bOptimizer)

#  void setSwitches( int doAngles, int doPiPiT, int  doPiSigma, int doPiPiI, int doBonded_, int PBC, int CheckInvariants )
lib.setSwitches.argtypes  = [c_int, c_int, c_int , c_int, c_int, c_int, c_int] 
lib.setSwitches.restype   =  None
def setSwitches(doAngles=0, doPiPiT=0, doPiSigma=0, doPiPiI=0, doBonded=0, PBC=0, CheckInvariants=0):
    return lib.setSwitches(doAngles, doPiPiT, doPiSigma, doPiPiI, doBonded, PBC, CheckInvariants)

#  bool checkInvariants( double maxVcog, double maxFcog, double maxTg )
lib.checkInvariants.argtypes  = [c_double, c_double, c_double] 
lib.checkInvariants.restype   =  c_bool
def checkInvariants(maxVcog=1e-9, maxFcog=1e-9, maxTg=1e-1):
    return lib.checkInvariants(maxVcog, maxFcog, maxTg)

#double eval()
lib.eval.argtypes  = [] 
lib.eval.restype   =  c_double
def eval():
    return lib.eval()

#bool relax( int niter, double Ftol )
lib.relax.argtypes  = [c_int, c_double, c_bool ] 
lib.relax.restype   =  c_bool
def relax(niter=100, Ftol=1e-6, bWriteTrj=False ):
    return lib.relax(niter, Ftol, bWriteTrj )

# ============= Manipulation

#void shift_atoms_ax( int n, int* sel, double* d                  )
lib.shift_atoms_ax.argtypes  = [c_int, c_int_p, c_double_p] 
lib.shift_atoms_ax.restype   =  None
def shift_atoms_ax(d, sel=None):
    n=0
    if sel is not None:
        n=len(sel)
        sel=np.array(sel,np.int32)
    if d is not None: d=np.array(d)
    return lib.shift_atoms_ax(n, _np_as(sel,c_int_p), _np_as(d,c_double_p))

#void shift_atoms( int n, int* sel, int ia0, int ia1, double l )
lib.shift_atoms.argtypes  = [c_int, c_int_p, c_int, c_int, c_double] 
lib.shift_atoms.restype   =  None
def shift_atoms( ia0, ia1, l, sel=None):
    n=0
    if sel is not None:
        n=len(sel)
        sel=np.array(sel,np.int32)
    return lib.shift_atoms(n, _np_as(sel,c_int_p), ia0, ia1, l)

#  rotate_atoms_ax( int n, int* sel, double* p0, double* ax, double phi      )
lib.rotate_atoms_ax.argtypes  = [c_int, c_int_p, c_double_p, c_double_p, c_double] 
lib.rotate_atoms_ax.restype   =  None
def rotate_atoms_ax( phi, sel=None, p0=None, ax=None ):
    n=0
    if sel is not None:
        n=len(sel)
        sel=np.array(sel,np.int32)
    if p0 is not None: p0=np.array(p0)
    if ax is not None: ax=np.array(ax)
    return lib.rotate_atoms_ax(n, _np_as(sel,c_int_p), _np_as(p0,c_double_p), _np_as(ax,c_double_p), phi)

#  rotate_atoms( int n, int* sel, int ia0, int iax0, int iax1, double phi )
lib.rotate_atoms_ax.argtypes  = [c_int, c_int_p, c_int, c_int, c_int, c_double] 
lib.rotate_atoms_ax.restype   =  None
def rotate_atoms_ax(ia0, iax0, iax1, phi, sel=None):
    n=0
    if sel is not None:
        n=len(sel)
        sel=np.array(sel,np.int32)
    return lib.rotate_atoms_ax(n, _np_as(sel,c_int_p), ia0, iax0, iax1, phi)

#  int splitAtBond( int ib, int* sel )
lib.splitAtBond.argtypes  = [c_int, c_int_p] 
lib.splitAtBond.restype   =  c_int
def splitAtBond(ib, sel=None):
    return lib.splitAtBond(ib, _np_as(sel,c_int_p))

#  void scanTranslation_ax( int n, int* sel, double* vec, int nstep, double* Es, bool bWriteTrj )
lib.scanTranslation_ax.argtypes  = [c_int, c_int_p, c_double_p, c_int, c_double_p, c_bool] 
lib.scanTranslation_ax.restype   =  None
def scanTranslation_ax(nstep, Es, bWriteTrj, sel=None, vec=None, ):
    n=0
    if sel is not None:
        n=len(sel)
        sel=np.array(sel,np.int32)
    if vec is not None: vec=np.array(vec)
    return lib.scanTranslation_ax(n, _np_as(sel,c_int_p), _np_as(vec,c_double_p), nstep, _np_as(Es,c_double_p), bWriteTrj)

#  void scanTranslation( int n, int* sel, int ia0, int ia1, double l, int nstep, double* Es, bool bWriteTrj )
lib.scanTranslation.argtypes  = [c_int, c_int_p, c_int, c_int, c_double, c_int, c_double_p, c_bool] 
lib.scanTranslation.restype   =  None
def scanTranslation( ia0, ia1, l, nstep, Es=None, sel=None, bWriteTrj=False):
    if Es is None: Es=np.zeros(nstep)
    return lib.scanTranslation(n, _np_as(sel,c_int_p), ia0, ia1, l, nstep, _np_as(Es,c_double_p), bWriteTrj)

#  void scanRotation_ax( int n, int* sel, double* p0, double* ax, double phi, int nstep, double* Es, bool bWriteTrj ){
lib.scanRotation_ax.argtypes  = [c_int, c_int_p, c_double_p, c_double_p, c_double, c_int, c_double_p, c_bool] 
lib.scanRotation_ax.restype   =  None
def scanRotation_ax( phi, nstep, sel=0, p0=None, ax=None,  Es=None, bWriteTrj=False):
    n=0
    if sel is not None:
        n=len(sel)
        sel=np.array(sel,np.int32)
    if p0 is not None: p0=np.array(p0)
    if ax is not None: ax=np.array(ax)
    lib.scanRotation_ax(n, _np_as(sel,c_int_p), _np_as(p0,c_double_p), _np_as(ax,c_double_p), phi, nstep, _np_as(Es,c_double_p), bWriteTrj)
    return Es

#  void scanRotation( int n, int* sel,int ia0, int iax0, int iax1, double phi, int nstep, double* Es, bool bWriteTrj ){
lib.scanRotation.argtypes  = [c_int, c_int_p, c_int, c_int, c_int, c_double, c_int, c_double_p, c_bool] 
lib.scanRotation.restype   =  None
def scanRotation( ia0, iax0, iax1, phi, nstep, sel=None, Es=None, bWriteTrj=False, _0=0):
    n=0
    if _0!=0: ia0 -=_0; iax0-=_0; iax1-=_0;
    if sel is not None:
        n=len(sel)
        sel=np.array(sel,np.int32)
        sel-=_0
    if Es is None: Es=np.zeros(nstep)
    lib.scanRotation(n, _np_as(sel,c_int_p), ia0, iax0, iax1, phi, nstep, _np_as(Es,c_double_p), bWriteTrj)
    return Es

def scanBondRotation( ib, phi, nstep, Es=None, bWriteTrj=False, bPrintSel=False):
    nsel = splitAtBond(ib) 
    if bPrintSel: print( "split to:\n", selection[:nsel],"\n", selection[nsel:] )
    ias = bond2atom[ib,:]
    return scanRotation( ias[0], ias[0], ias[1], phi, nstep, sel=None, Es=Es, bWriteTrj=bWriteTrj)





# ====================================
# ========= Python Functions
# ====================================

def plot(b2as=None,ax1=0,ax2=1,ps=None):
    import matplotlib.pyplot as plt
    from matplotlib import collections  as mc
    if ps is None: ps=apos
    if b2as  is None: b2as=bond2atom
    lines = [  ((ps[b[0],ax1],ps[b[0],ax2]),(ps[b[1],ax1],ps[b[1],ax2])) for b in b2as ]
    lc = mc.LineCollection(lines, colors='#808080', linewidths=2)
    ax=plt.gca()
    ax.add_collection(lc)
    plt.plot( ps[:,ax1], ps[:,ax2],'o', c='#8080FF' )
    for i,p in enumerate(ps):    ax.annotate( "%i"%i , (p[ax1],p[ax2]), color='b' )
    for i,l in enumerate(lines): 
        p= ((l[0][0]+l[1][0])*0.5,(l[0][1]+l[1][1])*0.5)
        ax.annotate( "%i"%i , p, color='k')
    plt.axis('equal')
    plt.grid()

def plot_selection(sel=None,ax1=0,ax2=1,ps=None, s=100):
    import matplotlib.pyplot as plt
    if sel is None: sel=selection
    if ps is None: ps=apos
    asel=ps[sel]
    plt.scatter( asel[:,ax1], asel[:,ax2], s=s, facecolors='none', edgecolors='r' )


# ====================================
# ========= MAIN
# ====================================

# if __name__ == "__main__":
#     import matplotlib.pyplot as plt
#     plt.show()