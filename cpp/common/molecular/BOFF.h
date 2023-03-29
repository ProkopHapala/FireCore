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

#include "arrayAlgs.h"


class BOFF{ public:

    int     natom   =0;
    Vec3d*  apos    =0;
    Vec3d*  aforce  =0;
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

double BOAngfunc( int ti, int tj, double r, double& dOij ){
    // this function evaluate Bond order from interactomic distance and atom type 
    //  * perhaps by spline interpolation ?
    //  * or we can fit it to exponential function
    return 0;
}


struct BObond{
    int       i;
    double  val;
    double dval;
    double    r;
};


double eval_atom(const int ia){
    Vec3d   pi    = apos   [ia];
    int     ti    = atype  [ia];
    //int*    ngs   = (int*   )(aneighs + ia); // do we need this ?
    //double* ngOs  = (double*)(Oneighs + ia);
    //int ngs   = Quat4i{-1,-1,-,1,-1};

    //int    ngs[nMaxNeigh];
    //double ngO[nMaxNeigh];

    BObond bonds[nMaxNeigh];
    
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
        if( Oij > 0 ){ // consider this as candidate for the bond
            ngs[nng]=BObond{ ja, Oij, dOij, r };
            nng++;
            if(nng>=nMaxNeigh){ printf("atom[%i] has > maximum number of neighbors %i => exit() \n", nMaxNeigh ); exit(0); }
        }
    }

    // ----- sort bonded atoms?    - is it worth it ?
    //int insertSort( nneigh, ngs );   // we sort it because we want to apply angular forces only on first few atoms

    // ------- Evaluate bond order energy and its derivative
    double dOi = Oi - Oi0s[ti];
    double FOi = K*dOi;
    double EOi = 0.5*FOi*dOi;

    const double OangMin = 0.5; 
    Quat4   hs[nneigh];
    double  AK[nneigh];
    int nang = 0;
    // ------ Apply forces due to bond order 
    for(int j=0; j<nneigh; j++){
        const BObond& B = ngs[j];
        Vec3d h = apos[B.i] - pi;
        double ivr = 1/B.r;
        h.mul( invr );
        
        fa.add_mul( h, FOi*B.dval );
        if(Bi.val<OangMin){
            hs[j].f = h;
            hs[j].e = ivr;
            double dAK;
            AK[j]   = BOAngfunc( int ti, atype[B.i], B.r, dAK );
            nang++;
        }
    }
    
    // ----- Angular forces betwee bonded atoms
    //       Cutoff for angular forces is shorter than cutoff for BO forces
    for(int i=0; i<nang; i++){
        const BObond& Bi = ngs[i];
        const Vec3d&  hi = hs [i];
        double Ki = AK[i]*Kang;
        for(int j=0; j<nang; j++){
            const BObond& Bj = ngs[j];
            const Vec3d&  hj = hs [j];
            double K = Ki*AK[j];
            E += evalAngleCos    ( hi.f, hj.f, hi.e, hj.e,  ssK*damp, ssC0, f1, f2 );  
            fa          .sub( f1  );
            fa          .sub( f2  );
            aforce[Bi.i].add( f1  ); // ToDo: This may be problem for paralelization
            aforce[Bj.i].add( f2  );
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

