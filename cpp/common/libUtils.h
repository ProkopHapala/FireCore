
#ifndef libUtils_h
#define libUtils_h

#include <vector>
#include <unordered_map>
#include <string>

std::unordered_map<std::string, double*>  buffers;
std::unordered_map<std::string, int*>     ibuffers;

extern "C"{

double* getBuff(const char* name){ 
    auto got = buffers.find( name );
    if( got==buffers.end() ){ printf( "buffer[%s] NOT FOUNT !!! \n", name ); return 0;        }
    else                    { return got->second; }
}

void setBuff(const char* name, double* buff){ 
    buffers[name] = buff;
}

int* getIBuff(const char* name){ 
    auto got = ibuffers.find( name );
    if( got == ibuffers.end() ){ return 0;        }
    else                    { return got->second; }
}

void setIBuff(const char* name, int* buff){ 
    ibuffers[name] = buff;
    //auto got = buffers.find( name );
    //if( got==buffers.end() ){ return null;        }
    //else                    { return got->second; }
}

};

#endif
