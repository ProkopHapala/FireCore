
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <vector>
#include <unordered_map>
#include <string>

#include <globals.h>
//int verbosity=0;
int idebug   =0;

#include "FF2D.h"
//#include "DynamicOpt.h"
#include "Vec3Utils.h"


FF2D ff;

// ============ Global Variables

extern "C"{

/*
void init_buffers(){
    // atoms (ions)
    buffers.insert( { "pDOFs",           ff.pDOFs  } );
    buffers.insert( { "fDOFs",           ff.fDOFs  } );
    buffers.insert( { "apos",   (double*)ff.apos   } );
    buffers.insert( { "aforce", (double*)ff.aforce } );
    buffers.insert( { "epos",   (double*)ff.epos   } );
    buffers.insert( { "eforce", (double*)ff.eforce } );
    buffers.insert( { "esize",  (double*)ff.esize  } );
    buffers.insert( { "fsize",  (double*)ff.fsize  } );
    buffers.insert ( { "aPars", (double*)ff.aPars  } );
    buffers.insert ( { "Es",            &ff.Etot   } );
    ibuffers.insert( { "espin",          ff.espin  } );
    ibuffers.insert( { "ndims",         &ff.ne     } );
}
*/

int init( char* str, int seed ){
    ff.insertString( str );
    ff.try_realloc();
    ff.cleanVelocity();
    if(seed>0){
        srand(seed);
        ff.setPosRandom(5.0);
    }
    return ff.atoms.size();
}

void toArrays( int* types, double* apos, int* neighs ){
    ff.toArrays( types, (Vec2d*)apos, (int4*)neighs );
}

double step( double dt, double damp ){ return ff.step(dt,damp); }    

int run( int n, double dt, double damp, double Fconv, bool bCleanForce){
    printf( "Fconv %g\n", Fconv);
    double F2conv=Fconv*Fconv;
    ff.try_realloc();
    if(bCleanForce)ff.cleanVelocity();
    for(int itr=0; itr<n; itr++){
        double F2sum = ff.step( dt, damp);
        printf( "run[%i] E %g |F| %g \n", itr, ff.Etot, sqrt(F2sum) );
        if(F2sum<F2conv)return itr;
    }
    return 0;
}

} // extern "C"
