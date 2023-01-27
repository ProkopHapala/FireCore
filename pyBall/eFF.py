
import numpy as np
from   ctypes import c_int, c_double, c_bool, c_float, c_char_p, c_bool, c_void_p, c_char_p
import ctypes
import os
import sys

sys.path.append('../')
from . import cpp_utils_ as cpp_utils

c_double_p = ctypes.POINTER(c_double)
c_int_p    = ctypes.POINTER(c_int)

def _np_as(arr,atype):
    if arr is None:
        return None
    else: 
        return arr.ctypes.data_as(atype)

cpp_utils.s_numpy_data_as_call = "_np_as(%s,%s)"

verbosity = 0
dt_glob   = 0.1
bVel      = False

# ===== To generate Interfaces automatically from headers call:
header_strings = [
#"void init_buffers(){",
#"bool load_xyz( const char* fname ){",
#"void init( int na, int ne ){",
#"void eval(){",
#"void info(){",
#"double* getEnergyPointer(){",
#"int*    getDimPointer   (){",
#"double* getBuff(const char* name){",
#"void setBuff(const char* name, double* buff){",
#"int* getIBuff(const char* name){",
#"void setIBuff(const char* name, int* buff){", 
#"void setPauliModel(int i){",
#"void setKPauli( double KPauli ){",
#"void initOpt( double dt, double damping, double f_limit ){",
#"int  run( int nstepMax, double dt, double Fconv=1e-6, int ialg=0 ){",
#"void writeTo_fgo( char const* filename, bool bVel, bool bAppend ){",
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
lib = cpp_utils.loadLib('eFF_lib', recompile=False)
array1ui = np.ctypeslib.ndpointer(dtype=np.uint32, ndim=1, flags='CONTIGUOUS')
array1i  = np.ctypeslib.ndpointer(dtype=np.int32,  ndim=1, flags='CONTIGUOUS')
array2i  = np.ctypeslib.ndpointer(dtype=np.int32,  ndim=2, flags='CONTIGUOUS')
array1d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS')
array2d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=2, flags='CONTIGUOUS')
array3d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=3, flags='CONTIGUOUS')
# ========= C functions

def cstr( s ):
    if s is None: return None
    return s.encode('utf8')

#  void setVerbosity( int verbosity_, int idebug_ ){
lib.setVerbosity.argtypes  = [c_int, c_int] 
lib.setVerbosity.restype   =  None
def setVerbosity(verbosity_=0, idebug=0):
    global verbosity
    verbosity = 1
    return lib.setVerbosity(verbosity, idebug)

#  void init_buffers(){
lib.init_buffers.argtypes  = [] 
lib.init_buffers.restype   =  None
def init_buffers():
    return lib.init_buffers() 

#  void load_xyz( const char* fname ){
lib.load_xyz.argtypes  = [c_char_p] 
lib.load_xyz.restype   =  c_bool
def load_xyz(fname):
    return lib.load_xyz(cstr(fname))
    #return lib.load_xyz(_np_as(fname,c_char_p)) 

def getBuffs( ):
    init_buffers()
    global ne,na,nDOFs, ndims, Es
    ndims = getIBuff( "ndims", (3,) ); 
    Es    = getBuff ( "Es",    (8,) ) # [0Etot,1Ek,2Eee,3EeePaul,4EeeExch,5Eae,6EaePaul,7Eaa]
    ne=ndims[0]; na=ndims[1]; nDOFs=ndims[2]     #;print("ne,na,nDOFs ", ne,na,nDOFs)
    global pDOFs, fDOFs, apos, aforce, epos, eforce, esize, fsize, aPars, espin
    pDOFs  = getBuff ( "pDOFs",  nDOFs )
    fDOFs  = getBuff ( "fDOFs",  nDOFs )
    apos   = getBuff ( "apos",   (na,3) )
    aforce = getBuff ( "aforce", (na,3) )
    epos   = getBuff ( "epos",   (ne,3) )
    eforce = getBuff ( "eforce", (ne,3) )
    esize  = getBuff ( "esize",   ne )
    fsize  = getBuff ( "fsize",   ne )
    aPars  = getBuff ( "aPars",   (na,4) )
    espin  = getIBuff( "espin",   ne )
    if(bVel):
        global vDOFs, avel, evel, vsize, invMasses, invAmass,invEmass,invSmass
        vDOFs     = getBuff ( "vDOFs",    nDOFs  )
        avel      = getBuff ( "avel",     (na,3) )
        evel      = getBuff ( "evel",     (ne,3) )
        vsize     = getBuff ( "vsize",     ne    )
        invMasses = getBuff ( "invMasses",nDOFs  )
        invAmass = getBuff ( "invAmass",  (na,3) )
        invEmass = getBuff ( "invEmass",  (ne,3) )
        invSmass = getBuff ( "invSmass",   ne    )

#  void load_xyz( const char* fname ){
lib.load_fgo.argtypes  = [c_char_p, c_bool] 
lib.load_fgo.restype   =  c_bool
def load_fgo(fname, bVel_=False):
    global bVel
    bVel=bVel_
    return lib.load_fgo( cstr(fname), bVel)

#  void save_fgo( char const* filename, bool bVel, bool bAppend ){
lib.save_fgo.argtypes  = [c_char_p, c_bool, c_bool] 
lib.save_fgo.restype   =  None
def save_fgo(filename, bVel=False, bAppend=False):
    return lib.save_fgo( cstr(filename), bVel, bAppend)

#  void setTrjName( char* trj_fname_ ){ 
lib.setTrjName.argtypes  = [c_char_p, c_int ] 
lib.setTrjName.restype   =  c_bool
def setTrjName(trj_fname_="trj.xyz", savePerNsteps=1, bDel=True ):
    if bDel: open(trj_fname_,"w").close()
    global trj_fname
    trj_fname=cstr(trj_fname_)
    return lib.setTrjName( trj_fname, savePerNsteps )

#  void init( int na, int ne ){
lib.init.argtypes  = [c_int, c_int] 
lib.init.restype   =  None
def init(na, ne, bVel_=False):
    global bVel
    bVel=bVel
    return lib.init(na, ne) 

#  void eval(){
lib.eval.argtypes  = [] 
lib.eval.restype   =  c_double
def eval():
    return lib.eval() 

#void evalFuncDerivs( int n, double* r, double* s, double* Es, double* Fs ){
lib.evalFuncDerivs.argtypes = [ c_int, array1d, array1d, array1d, array1d ]
lib.evalFuncDerivs.restype  = None
def evalFuncDerivs( r, s, Es=None, Fs=None ):
    r = r + s*0
    s = s + r*0
    n = len(r)
    if Es is None: Es=np.zeros(n)
    if Fs is None: Fs=np.zeros(n) 
    lib.evalFuncDerivs( n, r, s, Es, Fs )
    return Es,Fs

#  void info(){
lib.info.argtypes  = [] 
lib.info.restype   =  None
def info():
    return lib.info() 

#  double* getEnergyPointer(){
lib.getEnergyPointer.argtypes  = [] 
lib.getEnergyPointer.restype   = c_double_p
def getEnergyTerms( sh=(7,) ):
    # Ek=0, Eee EeePaul EeeExch Eae EaePaul Eaa
    ptr = lib.getEnergyPointer()
    return  np.ctypeslib.as_array( ptr, shape=sh )

#int*    getDimPointer   (){
lib.getDimPointer.argtypes  = [] 
lib.getDimPointer.restype   = c_int_p
def getDimPointer( sh=(3,) ):
    # ne=0 na=0 nDOFs=0
    ptr = lib.getDimPointer()
    return  np.ctypeslib.as_array( ptr, shape=sh )

#  double* getBuff(const char* name){
lib.getBuff.argtypes  = [c_char_p] 
lib.getBuff.restype   =  c_double_p
def getBuff( name, sh ):
    ptr = lib.getBuff(cstr(name))
    if not isinstance(sh, tuple): sh=(sh,)
    #sh_ = (natom,)
    #if sh is not None:
    #    sh_ = sh_ + sh
    #print "DEBUG type( ptr ) ", type( ptr ), sh
    return np.ctypeslib.as_array( ptr, shape=sh)

#  void setBuff(const char* name, double* buff){
lib.setBuff.argtypes  = [c_char_p, c_double_p] 
lib.setBuff.restype   =  None
def setBuff(name, buff):
    return lib.setBuff( cstr(name), _np_as(buff,c_double_p)) 
    #return lib.setBuff(_np_as(name,c_char_p), _np_as(buff,c_double_p)) 

#  int* getIBuff(const char* name){
lib.getIBuff.argtypes  = [c_char_p] 
lib.getIBuff.restype   =  c_int_p
def getIBuff(name,sh):
    ptr = lib.getIBuff(cstr(name))
    if not isinstance(sh, tuple): sh=(sh,)
    return np.ctypeslib.as_array( ptr, shape=sh)
    #return lib.getIBuff(_np_as(name,c_char_p)) 

#  void setIBuff(const char* name, int* buff){
lib.setIBuff.argtypes  = [c_char_p, c_int_p] 
lib.setIBuff.restype   =  None
def setIBuff(name, buff):
    return lib.setIBuff(name, _np_as(buff,c_int_p)) 
    #return lib.setIBuff(_np_as(name,c_char_p), _np_as(buff,c_int_p)) 

#  void setPauliModel(int i){
lib.setPauliModel.argtypes  = [c_int] 
lib.setPauliModel.restype   =  None
def setPauliModel(i):
    return lib.setPauliModel(i) 

#  void setKPauli( double KPauli ){
lib.setKPauli.argtypes  = [c_double] 
lib.setKPauli.restype   =  None
def setKPauli(KPauli):
    return lib.setKPauli(KPauli) 

#void setSwitches_(int bNormalize, int bNormForce, int bEvalKinetic, int bEvalCoulomb, int  bEvalExchange, int  bEvalPauli, int bEvalAA, int bEvalAE, int bEvalAECoulomb, int bEvalAEPauli ){
lib.setSwitches.argtypes = [ c_int, c_int, c_int, c_int, c_int, c_int, c_int ]
lib.setSwitches.restype  = None
def setSwitches( kinetic=0, coulomb=0, pauli=0, AA=0, AE=0, AECoulomb=0, AEPauli=0 ):
    lib.setSwitches( kinetic, coulomb, pauli, AA, AE, AECoulomb, AEPauli )


#  void initOpt( double dt, double damping, double f_limit ){
lib.initOpt.argtypes  = [c_double, c_double, c_double, c_bool ] 
lib.initOpt.restype   =  None
def initOpt(dt=0.1, damping=0.1, f_limit=1000.0, bMass=False ):
    global dt_glob
    dt_glob = dt
    return lib.initOpt(dt, damping, f_limit, bMass)

#  int  run( int nstepMax, double dt, double Fconv=1e-6, int ialg=0 ){
lib. run.argtypes  = [c_int, c_double, c_double, c_int] 
lib. run.restype   =  c_int
def  run(nstepMax=1000, dt=None, Fconv=1e-6, ialg=0):
    if dt is None: dt=dt_glob
    return lib.run(nstepMax, dt, Fconv, ialg)



# =========  Tests


def eval_mol(name ):
    load_fgo("data/"+name+".fgo" )                               # load molecule in  .fgo format (i.e. floating-gaussian-orbital)
    eval()

def relax_mol(name, dt=0.03,damping=0.1, bTrj=True, bResult=True, perN=1 ):
    load_fgo("data/"+name+".fgo" )                               # load molecule in  .fgo format (i.e. floating-gaussian-orbital)
    initOpt(dt=dt,damping=damping )                              # initialize optimizer/propagator
    if(bTrj): setTrjName(name+"_relax.xyz", savePerNsteps=perN ) # setup output .xyz file to save trajectory of all atoms and electrons at each timespep (comment-out to ommit .xyz and improve performance ) 
    run( 10000, Fconv=1e-3, ialg=2 )                             # run simuation for maximum 1000 time steps intil it converge to |F|<1e-3, ialg=2 is FIRE http://users.jyu.fi/~pekkosk/resources/pdf/FIRE.pdf   https://www.sciencedirect.com/science/article/pii/S0927025620300756
    if(bResult): 
        result_name=name+"_relaxed.fgo"
        if(verbosity>0): print("Optimized molecule saved to ", result_name)
        save_fgo( result_name )                 # save final relaxed geometry to .fgo format (i.e. floating-gaussian-orbital).

def printEs():
    print( " Etot %g Ek %g Eee %g EeePaul %g Eae%g EaePaul %g Eaa %g [A]" %(Es[0],Es[1],Es[2],Es[3],Es[5],Es[6],Es[7]) )  # [0Etot,1Ek,2Eee,3EeePaul,4EeeExch,5Eae,6EaePaul,7Eaa]

def printAtoms():
    for i in range(na):
        print( "atom[%i] Q %g apos(%g,%g,%g)" %(i, aPars[i,0], apos[i,0],apos[i,1],apos[i,2]) )

def printElectrons():
    for i in range(ne):
        print( "electron[%i] apos(%g,%g,%g) size %g spin %i" %(i, epos[i,0],epos[i,1],epos[i,2], esize[i], espin[i]) )

def check_H2(bRelax=True):
    from pyBall import eFF_terms as effpy
    if bRelax:
        relax_mol("H2_eFF")
    else:
        eval_mol("H2_eFF")
    getBuffs()
    bond_length = np.sqrt( ( (apos[0,:] - apos[1,:])**2 ).sum() ) ; print("apos\n",apos)
    Etot        = Es[0]
    printAtoms()
    printElectrons()
    print(  "effpy.run_H2_molecule.__doc__:\n", effpy.run_H2_molecule.__doc__ )
    print( "check_H2 E %g [eV] lbond %g [A]" %(Etot, bond_length) )
    printEs()
    
def init_eff( natom_=0, nelec_=1, s=0.5,  aQ=1.0,aQs=0.0,aP=0.0,aPs=0.1 ):
    global natom,nelec
    natom=natom_; nelec=nelec_; 
    init( natom, nelec )
    aPar  = getBuff( "aPars",(natom,4) )
    apos  = getBuff( "apos",(natom,3) )
    epos  = getBuff( "epos",(nelec,3) )
    esize = getBuff( "esize",(nelec)  )
    aPar [:,0]=aQ;aPar[:,1]=aQs;aPar[:,2]=aPs;aPar[:,3]=aP;
    apos [:,:] = 0
    epos [:,:] = 0
    esize[:]   = s
    '''
    epos [:,:,:]= 0              + (np.random.rand(norb,perOrb,3)-0.5)*rnd_pos
    esize[:,:  ]=sz              + (np.random.rand(norb,perOrb  )-0.5)*rnd_size
    ecoef[:,:  ]=1               + (np.random.rand(norb,perOrb  )-0.5)*rnd_coef
    rhoP [:,:,:]=0               + (np.random.rand(norb,nqOrb,3 )-0.5)*rnd_pos
    rhoS [:,:  ]=sz*np.sqrt(0.5) + (np.random.rand(norb,nqOrb   )-0.5)*rnd_size
    rhoQ [:,:  ]=1               + (np.random.rand(norb,nqOrb   )-0.5)*rnd_coef
    '''

def scan_size( ss, ie0 ):
    Escan = np.zeros(   (len(ss),len(Es)) )
    for i,s in enumerate(ss):
        esize[ie0] = s
        lib.eval()
        Escan[i,:] =  Es[:]
    return Escan

def test_Hatom( bDerivs=False ):
    from . import eFF_terms as effpy
    import matplotlib.pyplot as plt 
    init_eff( natom_=1, nelec_=1, s=0.5 )
    ss = np.arange( 0.3,1.0, 0.01 )

    print(  "effpy.run_Hatom.__doc__:\n", effpy.run_Hatom.__doc__ )
    Ek_ref,Eae_ref = effpy.Hatom( ss );  E_ref=Ek_ref+Eae_ref 
    if bDerivs:
        E,dE   = evalFuncDerivs(1,ss)
        xs=ss
        plt.plot(xs      ,dE                               ,'-',label="F_ana")
        plt.plot(xs[1:-1],(E[2:]-E[:-2])/(-2*(xs[1]-xs[0])),':',label="F_num")
    else:
        getBuffs()
        Es = scan_size( ss, 0 );   E=Es[:,0]; Ek=Es[:,1]; Eae=Es[:,5]     # [0Etot,1Ek,2Eee,3EeePaul,4EeeExch,5Eae,6EaePaul,7Eaa]
        plt.plot( ss, Eae,    '-r', label="Eae" )
        plt.plot( ss, Ek,     '-b', label="Ek"  )     
        plt.plot( ss, Eae_ref,':r', lw=2, label="Eae_ref" )
        plt.plot( ss, Ek_ref, ':b', lw=2, label="Ek_ref"  )   
    
    i0ref,E0ref,x0ref = effpy.getExmin1D(E_ref,ss) ; print( "Hatom opt Reference: E %g [eV] s %g [A]" %(E0ref,x0ref) )
    i0   ,E0   ,x0    = effpy.getExmin1D(E    ,ss) ; print( "Hatom opt Numerical: E %g [eV] s %g [A]" %(E0   ,x0   ) )

    plt.plot( ss, E_ref,':k',lw=3, label="E_ref" )
    plt.plot( ss, E,    'grey'    ,lw=3,label="E"     )

    plt.xlabel('size [A]')
    plt.legend()
    plt.grid()
    plt.show()

# ========= Python Functions

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    load_xyz("../../cpp/sketches_SDL/Molecular/data/e2_eFF.xyz")
    #load_xyz("../../cpp/sketches_SDL/Molecular/data/H2O_eFF.xyz")
    info()
    eval()

    plt.show()