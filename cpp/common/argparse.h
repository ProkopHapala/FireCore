
#ifndef  argparse_h
#define  argparse_h

#include <string>
#include <unordered_map>
#include <functional>

//typedef int MyInt;
using LambdaDict = std::unordered_map<std::string, std::function<void(const char*)>>;

void process_args( int argc, char *argv[], const LambdaDict& funcs ){
    for(int i=0; i<argc; i++ ){
        //printf( "process_args[%i] %s \n", i, argv[i] );
        auto found = funcs.find(argv[i]);
        if( found != funcs.end() ){
            //printf( "process_args[%i] (%s)->(%s) \n", i, found->first.c_str(), argv[i] );
            i++;
            found->second( argv[i] );
        }
    }
}

#endif
