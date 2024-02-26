
#ifndef libUtils_h
#define libUtils_h

#include <vector>
#include <unordered_map>
#include <string>

std::unordered_map<std::string, double*>  buffers;
std::unordered_map<std::string, float*>   fbuffers;
std::unordered_map<std::string, int*>     ibuffers;
std::unordered_map<std::string, bool*>    bbuffers;

extern "C"{

void printBuffNames(){
    printf( "#=== Buffers:  \n"); for (auto& it:  buffers) { printf( "%s: %li\n", it.first.c_str(), it.second ); };
    printf( "#=== fBuffers: \n"); for (auto& it: fbuffers) { printf( "%s: %li\n", it.first.c_str(), it.second ); };
    printf( "#=== IBuffers: \n"); for (auto& it: ibuffers) { printf( "%s: %li\n", it.first.c_str(), it.second ); };
    printf( "#=== bBuffers: \n"); for (auto& it: bbuffers) { printf( "%s: %li\n", it.first.c_str(), it.second ); };
}

double* getBuff(const char* name){ 
    auto got = buffers.find( name );
    if( got==buffers.end() )[[unlikely]]{ printf("ERROR: getBuff(%s) NOT FOUND\n", name); printBuffNames(); return 0;        }
    else                    { return got->second; }
}

void setBuff(const char* name, double* buff){ 
    buffers[name] = buff;
}

float* getfBuff(const char* name){ 
    auto got = fbuffers.find( name );
    if( got==fbuffers.end() )[[unlikely]]{ printf("ERROR: getfBuff(%s) NOT FOUND\n", name); printBuffNames(); return 0;        }
    else                    { return got->second; }
}

void setfBuff(const char* name, float* buff){ 
    fbuffers[name] = buff;
}

int* getIBuff(const char* name){ 
    auto got = ibuffers.find( name );
    if( got == ibuffers.end() )[[unlikely]]{ printf("ERROR: getIBuff(%s) NOT FOUND\n", name); printBuffNames(); return 0;        }
    else                       { return got->second; }
}

void setIBuff(const char* name, int* buff){ 
    ibuffers[name] = buff;
    //auto got = buffers.find( name );
    //if( got==buffers.end() ){ return null;        }
    //else                    { return got->second; }
}

bool* getBBuff(const char* name){ 
    auto got = bbuffers.find( name );
    if( got == bbuffers.end() )[[unlikely]]{ printf("ERROR: getBBuff(%s) NOT FOUND\n", name); printBuffNames(); return 0;        }
    else                       { return got->second; }
}

void setBBuff(const char* name, bool* buff){  bbuffers[name] = buff; }

};

#endif
