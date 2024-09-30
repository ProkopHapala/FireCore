
#ifndef  test_globals_h
#define  test_globals_h

#include <unordered_map>
#include <string>
#include "quaternion.h"

struct NDArray{
    double* data=0;
    Quat4i  dims=Quat4i{-1,-1,-1,-1};
};
static std::unordered_map<std::string,NDArray> golbal_array_dict;

#endif




