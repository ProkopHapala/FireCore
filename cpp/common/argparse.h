
#ifndef  argparse_h
#define  argparse_h

#include <string>
#include <unordered_map>
#include <functional>

//typedef int MyInt;
class ArgFunc{ public:
    int n;
    std::function<void(const char**)> func;
}; 
using LambdaDict = std::unordered_map<std::string, ArgFunc>;

void process_args( int argc, char *argv[], const LambdaDict& funcs ){
    for(int i=0; i<argc; i++ ){
        //printf( ".process_args[%i] %s \n", i, argv[i] );
        auto found = funcs.find(argv[i]);
        if( found != funcs.end() ){
            const ArgFunc& argf = found->second;
            //printf( "---process_args[%i] (%s) argf.n %i \n", i, found->first.c_str(), argf.n );
            const char** args = (const char** )(argv+i+1); 
            argf.func( args );
            i+=argf.n;
        }
    }
}

#endif
