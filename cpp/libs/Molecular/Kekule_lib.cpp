

constexpr int ntmpstr=2048;
char tmpstr[ntmpstr];

#include <globals.h>
//int verbosity = 1;
//int idebug    = 0;
//double tick2second=1e-9;

#include "Kekule.h"

Kekule W;

//============================

//const std::vector<std::string>* atomTypeNames = 0;
#include "libUtils.h"

extern "C"{

void setVerbosity( int verbosity_, int idebug_ ){
    verbosity = verbosity_;
    idebug    = idebug_;
}

void init_buffers(){
    // double* atomValence0 =0; // [natom] how many bonds atom can have?
    // double* atomValenceK =0; // [natom] stiffness of that bond order
    // Vec2i*  bond2atom    =0; // [nbond] atoms coresponding to the bond;
    // double* bondOrder0   =0; // [nbond] wanted obnd order of that bond 
    // double* bondOrderK   =0; // [nbond] stiffness of that bond order
    // // dynamical variables
    // double* atomValence  =0; // [natom]
    // double* atomValenceF =0; // [natom]
    // double* bondOrder    =0; // [natom]
    // double* bondOrderF   =0; // [nbond]
    buffers .insert( { "Es",                     &W.Etot           } );
    buffers .insert( { "Ks",                     &W.Katom          } );
    ibuffers.insert( { "ndims",                  &W.natom          } );
    ibuffers.insert( { "bond2atom",      (int*   )W.bond2atom      } );
    buffers .insert( { "atomValence",    (double*)W.atomValence    } );
    buffers .insert( { "atomValenceMin", (double*)W.atomValenceMin } );
    buffers .insert( { "atomValenceMax", (double*)W.atomValenceMax } );
    buffers .insert( { "bondOrder",      (double*)W.bondOrder        } );
    buffers .insert( { "bondOrderMin",   (double*)W.bondOrderMin   } );
    buffers .insert( { "bondOrderMax",   (double*)W.bondOrderMax   } );
}

void init( int natom, int nbond, int* bond2atom, double* atomValence0, int seed ){
    srand(seed);
    W.realloc( natom, nbond, (Vec2i*)bond2atom, 0,0, atomValence0, atomValence0 );
}

void setDefaultBondOrders( double min, double max ){
    W.setDefaultBondOrders( min, max );
}

void pinBondOrders(int n, int* ibonds, int* target ){
    W.pinBondOrders( n, ibonds, target  );
}

double eval(){ return W.eval(); }

double relax( double dt, double F2conv, int maxIter, bool bRandStart, int ialg ){
    return W.relax( dt, F2conv, maxIter, bRandStart, ialg );
};

double relaxComb( bool bRandStart ){
    if(verbosity>0) printf("// ---- Relax phase# 1\n");
    W.Katom=1;    W.Kbond=1;
    W.KatomInt=0; W.KbondInt=0;
    W.relax( 0.1, 1e-2, 10, bRandStart , 0);

    if(verbosity>0) printf("// ---- Relax phase# 2\n");
    W.KatomInt=0.1; W.KbondInt=0.1;
    W.relax( 0.5, 1e-2, 100, false , 1);

    if(verbosity>0) printf("// ---- Relax phase# 3\n");
    W.KatomInt=0.3; W.KbondInt=0.3;
    return W.relax( 0.5, 1e-4, 100, false , 1);

};

} // extern "C"
