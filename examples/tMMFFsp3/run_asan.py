import ctypes; from ctypes import c_int, c_double, c_bool, c_float, c_char_p, c_bool, c_void_p, c_char_p    ;print("DEBUG 1 ")
#asan = ctypes.CDLL(  "/usr/lib/gcc/x86_64-linux-gnu/11/libasan.so", mode=ctypes.RTLD_LOCAL )                ;print("DEBUG 2 ")
lib  = ctypes.CDLL(  "../../cpp/Build/libs/Molecular/libMMFFsp3_lib.so", mode=ctypes.RTLD_LOCAL )           ;print("DEBUG 3 ")


'''
array1i  = np.ctypeslib.ndpointer(dtype=np.int32,  ndim=1, flags='CONTIGUOUS')
def cstr( s ):
    if s is None: return None
    return s.encode('utf8')

# void* init( char* xyz_name, char* smile_name, int* nPBC, char* sAtomTypes, char* sBondTypes, char* sAngleTypes ){
lib.init.argtypes  = [c_char_p, c_char_p, array1i, c_char_p, c_char_p, c_char_p] 
lib.init.restype   =  c_void_p
def init( xyz_name  ="input.xyz", smile_name=None, sAtomTypes = "common_resources/AtomTypes.dat", sBondTypes = "common_resources/BondTypes.dat", sAngleTypes= "common_resources/AngleTypes.dat", nPBC=(1,1,0) ):
    nPBC=np.array(nPBC,dtype=np.int32)
    return lib.init( cstr(xyz_name),  cstr(smile_name), nPBC, cstr(sAtomTypes), cstr(sBondTypes), cstr(sAngleTypes) )


init( xyz_name="data/pyridine", surf_name="data/NaCl_sym-center", bMMFF=False  )  
'''
