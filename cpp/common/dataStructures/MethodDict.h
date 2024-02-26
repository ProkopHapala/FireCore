
#ifndef  MethodDict_h
#define  MethodDict_h

#include "containers.h"

template<typename MyClass>
class MethodDict : public Dict< void (MyClass::*)() >{ public:
    using Method = void (MyClass::*)();

    bool actionDispatch( MyClass* instance, int id ){
        Method method = Dict<void (MyClass::*)()>::vec[id];
        if( method ){ 
            (instance->*method)();
            return true; 
        }else{
            printf( "MethodDict::actionDispatch(%i): method not found. \n", id );
            return false;
        }
    }

};

#endif
