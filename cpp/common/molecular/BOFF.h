/*
Bond-Order Forcefield
https://github.com/ProkopHapala/FireCore/wiki/Bond-Order-potentials
*/

#ifndef BOFF_h
#define BOFF_h

#include "fastmath.h"
#include "Vec3.h"
//#include "quaternion.h"
#include "Forces.h"


class BOFF{ public:

    int     natom   =0;
    Vec3d*  apos    =0;
    double* aO      =0;  // atoms bond order
    int*    atype   =0;
    Quat4i* aneighs =0;  // indexes of 4 most bonded atoms 
    Quat4i* Oneighs =0;  // bond order contribution of 4 most bonded atoms


bool insertSorted( int ja, double Oij, int* aneighs, double* Oneighs ){
    return false;
};

double BOfunc( int ti, int tj, double r, double& dOij ){
    // this function evaluate Bond order from interactomic distance and atom type 
    //  * perhaps by spline interpolation ?
    //  * or we can fit it to exponential function
    return 0;
}

double eval_atom(const int ia){
    Vec3d   pi    = apos   [ia];
    int     ti    = atype  [ia];
    //int*    ngs   = (int*   )(aneighs + ia); // do we need this ?
    //double* ngOs  = (double*)(Oneighs + ia);
    //int ngs   = Quat4i{-1,-1,-,1,-1};

    int    ngs[nMaxNeigh];
    double ngO[nMaxNeigh];
    
    // ------- sum bond order over all atoms, find most bonded neighbors
    int nneigh = 0;
    double Oi  = 0;
    for(int ja=0; ja<na; ja++){
        if( ja == ia   )continue; // self-interactions
        Vec3d  d   = apos[ja]-pa;
        double r2  = d.norm2();
        if( r2 > R2cut )continue; // cutoff
        double r   = sqrt(r2);
        double dOij;
        double Oij  = BOfunc( ti, atype[ja], r, dOij );
        Oi += Oij; 
        if(Oij > min_Oij ){ // consider this as candidate for the bond
            //insertSorted( ja, Oij, aneighs, Oneighs );
            // maybe we can insert un-sorted and sort later?
            ngs[nng]=ja;
            ngO[nng]=ngO;
            nng++;
            if(nng>=nMaxNeigh){ printf("atom[%i] has > maximum number of neighbors %i => exit() \n", nMaxNeigh ); exit(0); }
        }
    }

    // ----- sort bonded atoms?

    // ------- Evaluate bond order energy and its derivative
    double dOi = Oi - Oi0s[ti];
    double FOi = K*dOi;
    double EOi = 0.5*FOi*dOi;

    // ------ Apply forces due to bond order 
    for(int ja=0; ja<na; ja++){
        if(ja==ia)continue; // self-interactions
        Vec3d  d   = apos[ja]-pa;
        double r2  = d.norm2();
        if( r2 > R2cut )continue; // cutoff
        double r   = sqrt(r2);
        double dOij;
        BOfunc( ti, atype[ja], r, dOij );
        fa.add_mul( d, FOi*dOij/r );
    }

    // ----- Angular forces betwee bonded atoms
    for(int i=0; i<4; i++){
        int ig = ngs[i];
        if(ig<0) break;
        for(int j=0; j<4; j++){
            int jg = ngs[j];
            if(ig<0) break;
            // bond betwee bonded atoms
        }
    }

    return E;
}

double eval_atoms(){
    double E=0;
    for(int ia=0; ia<nnode; ia++){ 
        E+=eval_atom(ia); 
    }
    return E;
}


}

#endif

