
#ifndef  script_parser_h
#define  script_parser_h

#include <any>
#include <functional>
#include <string>
#include <unordered_map>
#include "Tree.h"

// see also cpp/common/repl.h and cpp/common/argparse.h
#include <repl.h>
#include <argparse.h>

// from argparse.h
// class ArgFunc{ public:
//     int n;
//     std::function<void(const char**)> func;
// }; 
// using LambdaDict = std::unordered_map<std::string, ArgFunc>;

class code_segment{ public:
    int level;
    char* code;
    code_segment* parent;
    code_segment* children;
    code_segment* next;
    code_segment* prev;
};




class ScripParser{ public:
    LambdaDict functions;
    int nmaxchar = 1024*16;
    int nlevels  = 5;
    char* separators = " ,;:|\n";
    //char* sepMap     = new char[256];

    char* program_buff = new char[nmaxchar];
    //char* line    = new char[nmaxchar];

    void parse_file_inRAM(FILE* file) {  
        fseek(file, 0, SEEK_END);    // got to end of file
        long fileSize = ftell(file); // get current position (which is the size)
        fseek(file, 0, SEEK_SET);    // go back to beginning of file
        char* program_buff = (char*)malloc(fileSize + 1); // WARRNING: this may be too much memory if file is large !!!!
        if (program_buff) {
            fread(program_buff, fileSize, 1, file);       // WARRNING:  Read the entire file into the buffer, it can crash if file is too large
            program_buff[fileSize] = '\0'; // Null-terminate the buffer
            int ntot=0;
            while(ntot<fileSize){
                const int n += parse_line(&program_buff[ntot]);
                ntot += n;
                if(n==0){break;} // Prevent infinite loop if parse_line returns 0
            }
            //free(program_buff);
        } else {
            // Handle error: malloc failed
        }
    }
    int parse_line(char* line){
        if( line[0] == '#' ){ return 0; } // comment line
        for(int i=0; i<nmaxchar; i++ ){
            const char c = line[i];
            switch(c){
                // explicit fallthrough
                case '\n': [[fallthrough]];
                case '\0': {
                    return i;
                }
                case '\t':
                //case '\t': 
                case ' ':   // level 0
                case ',':   // level 1
                case ';':   // level 2
                case ':':   // level 3
                //case '|': // level 4
                if(i>0){ }
                break;
            }
            //if( (c=='\n')||(c=='\0')) ){ return i; } // end of line
        }
    }

};

#endif
