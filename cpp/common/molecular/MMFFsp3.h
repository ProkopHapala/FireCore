
#ifndef MMFFsp3_h
#define MMFFsp3_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"

// ======================
// ====   MMFFsp3
// ======================

inline void pairs2triple( const  Vec2i b1, const Vec2i b2, Vec3i& tri, bool& flip1, bool& flip2 ){
    if     ( b1.y == b2.x ){ tri.set( b1.x, b1.y, b2.y ); flip1=false; flip2=false;  }
    else if( b1.y == b2.y ){ tri.set( b1.x, b1.y, b2.x ); flip1=false; flip2=true;   }
    else if( b1.x == b2.x ){ tri.set( b1.y, b1.x, b2.y ); flip1=true;  flip2=false;  }
    else if( b1.x == b2.y ){ tri.set( b1.y, b1.x, b2.x ); flip1=true;  flip2=true;   }
}

class MMFFsp3{ public:
    constexpr const int nneigh_max = 4;
    int  nDOFs=0,natoms=0, nbonds=0, nang=0, ncap0=0, npi=0;
    bool bPBC=false;

    // 
    double * DOFs  = 0;   // degrees of freedom
    double * fDOFs = 0;   // forces

    bool doPi=0;
    int ipi0;
    Vec3d * apos=0;
    Vec3d * pipos=0;
    //Vec3d * cappos=0;

    Vec3d * fapos=0;
    Vec3d * fpipos=0;
    //Vec3d * fcappos=0;

    // --- Parameters
    Vec2i  * bond2atom = 0;
    double * bond_l0   = 0;  // [A]
    double * bond_k    = 0;  // [eV/A] ?
    Vec3d  * pbcShifts = 0;  // [A]

    int*    aneighs = 0;  // [natom*nneigh_max] neigbors of atom
    double* Kneighs = 0;  // [natom*nneigh_max] neigbors of atom

    // --- Axuliary variables
    double * lbond  = 0;   // bond lengths
    Vec3d  * hbond  = 0;   // normalized bond unitary vectors

void realloc( int natoms_, int nbonds_, int nang_, int ntors_ ){
    natoms=natoms_; nbonds=nbonds_; nang=nang_; ntors=ntors_;
    nDOFs = (natoms + npis)*3;
    _realloc( DOFs     , nDOFs );
    _realloc( fDOFs    , nDOFs );

    _realloc( lbond     , nbonds );
    _realloc( hbond     , nbonds );

    _realloc( bond2atom , nbonds );
    _realloc( bond_l0   , nbonds );
    _realloc( bond_k    , nbonds );

}

void cleanAtomForce(){ for(int i=0; i<natoms; i++){ aforce[i].set(0.0); } }

// ============== Evaluation


int i_DEBUG = 0;


inline double evalPi( const Vec2i& at, int ipi, int jpi, double K ){  // interaction between two pi-bonds
    // E = -K*<pi|pj>^2
    Vec3d d =  
}

inline double evalAngle_dist( int ia, int ing, int jng, double K ){
    //  E = -K*|pi-pj|
    // ToDo: This may be made more efficient if we store hbonds
    Vec3d d  = apos[ing] - apos[jng];   d .normalize();
    Vec3d hi = apos[ing] - apos[ia];    hi.normalize();
    Vec3d hj = apos[jng] - apos[ia];    hj.normalize();
    Vec3d fi  = d* K;
    Vec3d fj  = d*-K;
    Vec3d fiT = hi* hi.dot(fi);   fi.sub(fiT);
    Vec3d fjT = hj* hj.dot(fj);   fj.sub(fjT);
    aforce[ia] .add(fiT);
    aforce[ia] .add(fiT);
    aforce[ing].add(fi );
    aforce[jng].add(fj );
}

void evalAngle_cos(  i,  ia, ja, ib, jb, double K){
    Vec2i ib = ang2bond[ig];
    Vec3i ia = ang2atom[ig];
    Vec3d h1,h2;

    if(ib.x&SIGN_MASK){ ib.x&=0xFFFF; h1 = hbond[ib.x]; h1.mul(-1.0d); }else{ h1 = hbond[ib.x]; };
    if(ib.y&SIGN_MASK){ ib.y&=0xFFFF; h2 = hbond[ib.y]; h2.mul(-1.0d); }else{ h2 = hbond[ib.y]; };

    double c = h1.dot(h2);

    Vec3d hf1,hf2; // unitary vectors of force
    hf1 = h2 - h1*c;
    hf2 = h1 - h2*c;

    double fang = -K/(1.02-c);
    hf1.mul( fang/lbond[ib] );
    hf2.mul( fang/lbond[jb] );
    aforce[ia.x].add( hf1     );
    aforce[ia.y].add( hf2     );
    aforce[ia.z].sub( hf1+hf2 );
}

double eval_bond(int ib){
    //printf( "bond %i\n", ib );
    Vec2i iat = bond2atom[ib];
    Vec3d f; f.set_sub( apos[iat.y], apos[iat.x] );
    if(pbcShifts)f.add( pbcShifts[ib] );
    //printf( "bond[%i|%i,%i] (%g,%g,%g) (%g,%g,%g) \n", ib,iat.a,iat.b, apos[iat.x].x, apos[iat.x].y, apos[iat.x].z, apos[iat.y].x, apos[iat.y].y, apos[iat.y].z );
    double l = f.normalize();
    //printf( "bond[%i|%i,%i] (%g,%g,%g) %g \n", ib,iat.a,iat.b, f.x, f.y, f.z, l );
    lbond [ib] = l;
    hbond [ib] = f;
    if(bondMasked)if(bondMasked[ib])return 0;
    const double k = bond_k[ib];
    double dl = (l-bond_l0[ib]);
    f.mul( dl*k*2 );
    aforce[iat.x].add( f );
    aforce[iat.y].sub( f );
    double E = k*dl*dl;
    double Epi = 0;
    if(doPi){ // interaction between pi-bonds of given atom
        int i0=iat.a+1*nneigh_max;
        int j0=iat.b+1*nneigh_max;
        for(int i=-1;i>=-2;i--){
            int ipi=aneighs[i0+i];
            if(ipi<ipi0) continue;
            for(int j=-1;j>=-2;j--){
                int jpi=aneighs[j0+j];
                if(jpi<ipi0) continue;
                Epi += evalPi( at, ipi,jpi, Kneighs[i0+i]*Kneighs[j0+j] );
            }
        }
    }
    Epis += Epi;
    return E+Epi;
}

double eval_neighs(int ia){
    double E=0;
    for(int i=0; i<nneigh_max; i++){
        int ing = aneighs[i];
        if( ing<0 ){ break; // empty
        }else if( ing<i0pi ){ // signa-bonds
            for(int j=0; j<i; j++){
                int jng = aneighs[i];
                evalAngle( ia, ing, jng  );
            }
        }
    }
    return E;
}

double eval_bonds(){
    double E=0;
    for(int ib=0; ib<nbonds; ib++){  E+=eval_bond(ib); }
    return E;
}

double eval_neighs(){
    double E=0;
    for(int ia=0; ia<natoms; ia++){ E+=eval_neighs(ia); }
    return E;
}

double eval( bool bClean=true ){
    if(bClean)cleanAtomForce();
    Eb = eval_bonds();
    Ea = eval_angles();
    return Eb+Ea+Et;
};

//double getFmax(){ double Fmax=0; for(int i=0; i<natoms;i++){ _max( Fmax, aforce[i].norm2() ); }; return Fmax; };
//void moveGD(double dt){ for(int i=0; i<natoms;i++){ apos[i].add_mul( aforce[i],dt); } };

}


}; // MMFF

#endif
