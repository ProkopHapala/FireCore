
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
//int idebug   =0;

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

void getBonds( int* bonds ){
    ff.getBonds( (Vec2i*)bonds );
}

double step( double dt, double damp ){ return ff.step(dt,damp); }    

int run( int niter, double dt, double damp, double Fconv, bool bCleanForce){
    return ff.run( niter, dt, damp, Fconv, bCleanForce );
}

int getAtomNumber(){ return ff.atoms.size(); }
int getBondNumber(){ return ff.bonds.size(); }

int addAtom( double x, double y, int type=0, int* neighs=0 ){ return ff.addAtom( Vec2d{x,y}, type, neighs );}

int addBond( int ia, int ja ){ return ff.addBond(ia,ja); }

bool removeAtom(int i){ 
    //return ff.removeAtom(i); 
    //ff.swapAtoms( i, ff.atoms.size()-1 );
    //ff.atoms.pop_back();
    ff.removeAtomSwap( i );
    ff.bondsFromNeighs();
    return true;
}
bool removeBond(int i){ return ff.removeBond(i); }
int findBondAt( double x, double y, double R ){ return ff.findBondAt( {x,y},R); }
int findAtomAt( double x, double y, double R ){ return ff.findAtomAt( {x,y},R); }


void print_atoms(){ ff.print_atoms(); }

} // extern "C"
