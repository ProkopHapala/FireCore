#ifndef FireCore_h
#define FireCore_h

#include <dlfcn.h>

namespace FireCore{

// subroutine firecore_evalForce( nmax_scf, forces_ )  bind(c, name='firecore_evalForce')
//  subroutine firecore_init( natoms_, atomTypes, atomsPos ) bind(c, name='firecore_init')
//typedef void (*Pprocedure)();
//typedef void (*Pfirecore_evalForce)(int,double*);

typedef void (*P_evalForce  )(int,double*,double*);
typedef void (*P_init       )(int,int*,double*);
typedef void (*P_getCharges )(double*);

class Lib{ public:
    void        *lib_handle = 0;
    P_evalForce  evalForce =0;
    P_init       init      =0;
    P_getCharges getCharges=0;

    inline int loadLib( const char* fullname ){
        char *error;
        if(lib_handle){ dlclose(lib_handle); lib_handle=0; };
        //char fullname[1024] = "/home/prokop/git/FireCore/build/libFireCore.so";
        printf("loading %s :\n", fullname );
        void* plib = dlopen( fullname, RTLD_LAZY | RTLD_GLOBAL );
        if (plib){
            lib_handle=plib;
        }else{ printf( "%s\n", dlerror()); return -1; }
        init      = (P_init)     dlsym(lib_handle, "firecore_init");
        if ((error = dlerror())){ printf( "%s\n", error); init    =0; exit(0); }
        evalForce = (P_evalForce)dlsym(lib_handle, "firecore_evalForce");
        if ((error = dlerror())){ printf("%s\n", error); evalForce=0; exit(0); }
        getCharges = (P_getCharges)dlsym(lib_handle, "firecore_getCharges");
        if ((error = dlerror())){ printf("%s\n", error); getCharges=0; exit(0); }
        return 0;
    }

};

}

#endif