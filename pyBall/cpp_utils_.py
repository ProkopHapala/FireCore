
import os
import ctypes

clean_build = True 
recompile_glob = True
lib_ext   ='.so'
s_numpy_data_as_call = "_np_as(%s,%s)"

def work_dir( v__file__ ): 
    return os.path.dirname( os.path.realpath( v__file__ ) )

PACKAGE_PATH = work_dir( __file__ )
BUILD_PATH   = os.path.normpath( PACKAGE_PATH + '../../../cpp/Build/libs/CombatModels' )

#print (" PACKAGE_PATH : ", PACKAGE_PATH)
#print (" BUILD_PATH   : ", BUILD_PATH)

def compile_lib( name,
        #FFLAGS = "-std=c++11 -Og -g -Wall",
        FFLAGS = "-std=c++11 -O3 -ftree-vectorize -unroll-loops -ffast-math",
        LFLAGS = "-I/usr/local/include/SDL2 -lSDL2",
        path   = BUILD_PATH,
        clean  = True,
    ):
    lib_name = name+lib_ext
    print (" COMPILATION OF : "+name)
    if path is not None:
        dir_bak = os.getcwd()
        os.chdir( path );
    print (os.getcwd())
    if clean:
        try:
            os.remove( lib_name  )
            os.remove( name+".o" ) 
        except:
            pass 
    os.system("g++ "+FFLAGS+" -c -fPIC "+name+".cpp -o "+name+".o "+LFLAGS )
    os.system("g++ "+FFLAGS+" -shared -Wl,-soname,"+lib_name+" -o "+lib_name+" "+name+".o "+LFLAGS)
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

def loadLib( cpp_name, recompile=True, mode=ctypes.RTLD_LOCAL ):
    if recompile and recompile_glob:
        make(cpp_name)
    return ctypes.CDLL(  BUILD_PATH + "/lib" + cpp_name + lib_ext, mode ) 


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






