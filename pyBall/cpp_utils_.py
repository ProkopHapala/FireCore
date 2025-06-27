import os
import ctypes

loaded_libs = {} 

clean_build    = True 
recompile_glob = True
lib_ext        ='.so'
s_numpy_data_as_call = "_np_as(%s,%s)"

def work_dir( v__file__ ): 
    return os.path.dirname( os.path.realpath( v__file__ ) )

PACKAGE_PATH = work_dir( __file__ )

# Base root path (FireCore repository root)
BASE_PATH = os.path.normpath(os.path.join(PACKAGE_PATH, '../'))
# Default library path
BUILD_PATH = os.path.normpath(os.path.join(BASE_PATH, 'cpp/Build/libs/'))

#print (" PACKAGE_PATH : ", PACKAGE_PATH)
#print (" BUILD_PATH   : ", BUILD_PATH)

def set_args_dict( lib, argDict):
    for k,v in  argDict.items():
        f = lib.__getattr__(k)
        f.restype  = v[0]
        f.argtypes = v[1]
        

def compile_lib( name,
        #FFLAGS = "-std=c++11 -Og -g -Wall",
        #FFLAGS = "-std=c++11 -O3 -ftree-vectorize -unroll-loops -ffast-math",
        FFLAGS = "-std=c++11 -Ofast",
        LFLAGS = "-I/usr/local/include/SDL2 -lSDL2",
        path   = BUILD_PATH,
        clean  = True,
        is_cuda=False
    ):
    lib_name = name+lib_ext
    print (" COMPILATION OF : "+name)
    if path is not None:
        dir_bak = os.getcwd()
        os.chdir( path )
    print (os.getcwd())
    if clean:
        try:
            os.remove( lib_name  )
            os.remove( name+".o" ) 
        except:
            pass 

    if is_cuda:
        # Allow unsupported compiler versions and add header include paths
        cpp_root = os.path.normpath(os.path.dirname(path) + "/../..")
        include_paths = f"-I. -I{path} -I{cpp_root} -I{cpp_root}/common -I{cpp_root}/common/dataStructures"
        cuda_flags = f"-allow-unsupported-compiler {include_paths} -x cu -std=c++11"
        linker_flags = "--linker-options -lstdc++,-rpath,$ORIGIN"
        
        print(f"CPP root path: {cpp_root}")
        print(f"Include paths: {include_paths}")
        print(f"CUDA flags: {cuda_flags}")
        
        # Compile step with detailed output
        compile_cmd = f"nvcc {cuda_flags} -Xcompiler \"-fPIC -fvisibility=default\" -c {name}.cpp -o {name}.o"
        print(f"Compile command: {compile_cmd}")
        compile_result = os.system(compile_cmd)
        if compile_result != 0:
            print(f"Error: Compilation failed with code {compile_result}")
        
        # Link step with more aggressive C++ standard library linking
        link_cmd = f"nvcc -shared {linker_flags} -o {lib_name} {name}.o -lcudart -lstdc++"
        print(f"Link command: {link_cmd}")
        link_result = os.system(link_cmd)
        if link_result != 0:
            print(f"Error: Linking failed with code {link_result}")
    else:
        # Regular C++ compilation
        os.system(f"g++ {FFLAGS} -c -fPIC {name}.cpp -o {name}.o {LFLAGS}")
        os.system(f"g++ {FFLAGS} -shared -Wl,-soname,{lib_name} -o {lib_name} {name}.o -lstdc++ {LFLAGS}")

    if path is not None:
        os.chdir( dir_bak )

def make( what="" ):
    current_directory = os.getcwd()
    os.chdir ( BUILD_PATH          )
    print ("CPP_PATH " + BUILD_PATH)
    if clean_build:
        os.system("make clean")
    os.system( "make "+what      )
    os.chdir ( current_directory )

def loadLib( cpp_name, recompile=True, mode=ctypes.RTLD_LOCAL, is_cuda=False, path=None ):
    if recompile and recompile_glob:  
        if is_cuda:
            compile_lib(cpp_name, is_cuda=True, path=path)
        else:
            make(cpp_name)
    # If path is None, use the default BUILD_PATH
    if path is None:
        path = BUILD_PATH
        
    # First try the provided path
    lib_path = os.path.join(path, "lib" + cpp_name + lib_ext)
    
    # If file doesn't exist, try to find it in common locations
    if not os.path.exists(lib_path):
        # Try Molecular subdirectory
        molecular_path = os.path.join(path, 'Molecular')
        molecular_lib_path = os.path.join(molecular_path, "lib" + cpp_name + lib_ext)
        
        if os.path.exists(molecular_lib_path):
            lib_path = molecular_lib_path
        else:
            # Try using absolute path from BASE_PATH
            base_lib_path = os.path.join(BASE_PATH, 'cpp/Build/libs/Molecular', "lib" + cpp_name + lib_ext)
            if os.path.exists(base_lib_path):
                lib_path = base_lib_path
    print(f"Loading library from: {lib_path}")
    if lib_path in loaded_libs: 
        unload_lib_by_path(lib_path)  # Unload if already loaded
    lib = ctypes.CDLL(lib_path, mode) #nacita c knihovnu
    loaded_libs[lib_path] = lib  # Store the loaded library
    return lib

def unload_lib_by_path(lib_path):
    if lib_path in loaded_libs:
        lib = loaded_libs.pop(lib_path)  # Remove from dictionary
        try:
            # Attempt to find and use dlclose (more common on Unix-like systems)
            unload_lib(lib)
            print(f"Library {lib_path} unloaded successfully.")
        except (AttributeError, OSError) as e:
            print(f"Warning: Could not unload library {lib_path} properly.")
            print(f"Error: {e}")

def unload_lib(lib):
    dlclose_func          = ctypes.CDLL(None).dlclose
    dlclose_func.argtypes = [ctypes.c_void_p]
    dlclose_func(lib._handle)


# ============ automatic C-python interface generation

def _np_as(arr,atype):
    if arr is None:
        return None
    elif isinstance( arr, str ):
        #arr = arr.encode('utf-8')
        #print "arr = ", arr
        #return arr
        return arr.encode('utf-8')
    else: 
        return arr.ctypes.data_as(atype)

def parseFuncHeader( s ):
    #args = []
    arg_types = []
    arg_names = []
    i0 = s.find(' ')
    i1 = s.find('(')
    i2 = s.find(')')
    ret_type = s[:i0]
    name     = s[i0+1:i1]
    sargs = s[i1+1:i2]
    largs = sargs.split(',')
    for sarg in largs:
        sarg = sarg.strip()
        i    = sarg.rfind(' ')
        #arg_name=sarg[i+1:]
        #arg_type=sarg[  :i]
        #args.append( (arg_type, arg_name) )
        arg_types.append(sarg[  :i])
        arg_names.append(sarg[i+1:])
    return (name,ret_type,arg_types,arg_names)

def translateTypeName( tname ):
    np = tname.count('*')
    if np > 1 :
        print ("Cannot do pointer-to-pointer (**) ", s)
        print ("=> exit() ")
        exit()
    else:
        if(np==1):
            i  = tname.find ('*')
            return "c_"+tname[:i]+"_p", True
        else:
            return "c_"+tname     , False

s_numpy_data_as_call = "%s.ctypes.data_as(%s)"

def writePointerCall( name, ttype ):
    if ttype[1]:
        return s_numpy_data_as_call %(name,ttype[0])
    else:
        return name

def writeFuncInterface( parsed ):
    name,ret_type,arg_types,arg_names = parsed
    arg_types_ = [ ]
    #arg_names = [ ]
    if ret_type=="void" : 
        ret_type="None"
    else:
        ret_type=translateTypeName(ret_type)[0]
    sdef_args  = ", ".join( [                    translateTypeName(t)[0] for   t in arg_types                ] )
    scall_args = ", ".join( [ writePointerCall(n,translateTypeName(t))   for n,t in zip(arg_names,arg_types) ] )
    lines = [
        "lib."+name+".argtypes  = ["+ sdef_args + "] ",
        "lib."+name+".restype   =  " + ret_type          ,
        "def "+name+"("+ ", ".join(arg_names) +"):"    ,
        "    return lib."+name+"("+scall_args+")"         ,
    ]
    return "\n".join( lines )

def writeFuncInterfaces( func_headers, debug=False ):
    for s in func_headers:
        parsed = parseFuncHeader( s ); 
        if debug : print ("parsed :\n", parsed)
        sgen   = writeFuncInterface( parsed )
        print ("\n# ", s)
        #print (sgen,"\n\n")
        print (sgen)
