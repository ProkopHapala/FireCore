
#ifndef  ScriptLang_h
#define  ScriptLang_h

//#include <vector>
//#include <string>
//#include <unordered_map>

#include "contains.h"

class Variable{ public:
    std::string name;
    int   type;
    void* data=0;
};

class Function{ public:
    std::string name;
    std::vector<int> args;
    int returnType;
    std::vector<Variable> args;
    std::string body;
};

class ScriptLang{ public:
    Dict<Function> functions;
    Dict<Variable> variables;    

    

};

#endif
