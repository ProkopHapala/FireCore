
#ifndef libUtils_h
#define libUtils_h

#include <vector>
#include <unordered_map>
#include <string>

std::unordered_map<std::string, double*>  buffers;
std::unordered_map<std::string, float*>   fbuffers;
std::unordered_map<std::string, int*>     ibuffers;

extern "C"{

void printBuffNames(){
    printf( "#=== Buffers:  \n"); for (auto& it:  buffers) { printf( "%s: %li\n", it.first.c_str(), it.second ); };
    printf( "#=== fBuffers: \n"); for (auto& it: fbuffers) { printf( "%s: %li\n", it.first.c_str(), it.second ); };
    printf( "#=== IBuffers: \n"); for (auto& it: ibuffers) { printf( "%s: %li\n", it.first.c_str(), it.second ); };
}

double* getBuff(const char* name){ 
    auto got = buffers.find( name );
    if( got==buffers.end() ){ printf("ERROR: getBuff(%s) NOT FOUND\n", name); printBuffNames(); return 0;        }
    else                    { return got->second; }
}

void setBuff(const char* name, double* buff){ 
    buffers[name] = buff;
}

float* getfBuff(const char* name){ 
    auto got = fbuffers.find( name );
    if( got==fbuffers.end() ){ printf("ERROR: getfBuff(%s) NOT FOUND\n", name); printBuffNames(); return 0;        }
    else                    { return got->second; }
}

void setfBuff(const char* name, float* buff){ 
    fbuffers[name] = buff;
}

int* getIBuff(const char* name){ 
    auto got = ibuffers.find( name );
    if( got == ibuffers.end() ){ printf("ERROR: getIBuff(%s) NOT FOUND\n", name); printBuffNames(); return 0;        }
    else                       { return got->second; }
}

void setIBuff(const char* name, int* buff){ 
    ibuffers[name] = buff;
    //auto got = buffers.find( name );
    //if( got==buffers.end() ){ return null;        }
    //else                    { return got->second; }
}

};

#endif
