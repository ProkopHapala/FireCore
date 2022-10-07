
#ifndef MMFFsp3_h
#define MMFFsp3_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "quaternion.h"
#include "molecular_utils.h"

// ======================
// ====   MMFFsp3_loc
// ======================

/*
 MMFFsp3_loc is atempt to rewrite MMFFsp3 with more memory-local datastructures
 this should make it more suitable for paralelization and implementation on GPU using OpenCL with local memory cache. 
*/


#define CAP_PI -1

class MMFFsp3_loc{ public:
    static constexpr const int nneigh_max = 4;
    int  nDOFs=0,nnode=0,ncap=0,nps=0;
    
    bool bPBC=false;
    double Kpipi = 0.25;

    double * DOFs  = 0;   // degrees of freedom
    double * fDOFs = 0;   // forces

    // Node atoms
    int   * atype=0;
    Vec3d  *  apos=0;  // node atom positions
    Vec3d  * fapos=0;  //
    double * aKs=0;    // angular strengh of atoms
    double * a0s=0;    // equilibirum angle

    // Neighbors (cap-atoms, pi-vectors)[natom*nneighmax];
    Vec3d *  cpos=0;     // capping atom poistions
    Vec3d * fcpos=0;     // 
    double* lcpos=0;      // length of bond to neighbor;

    double* l0cpos=0;    // length of bond to neighbor;
    double* kcpos =0;   
    int*    neighs = 0;  // neigbor indexes, [+] index of node atoms [-] type of cap
    
    // ----- DEBUG / Tweaks / Modifiees
    double Etot,Eb,Ess,Esp,EppT,EppI;
    bool bPPI =true;
    bool bPPT =true;
    bool bSP  =true;
    bool bSS  =true;
    int  nPPI =0;
    int  nPPT =0;
    int  nSP  =0;
    int  nSS  =0;
    int  nB   =0;
    bool bDEBUG_plot=0;
    int iDEBUG_n=0,iDEBUG_pick=0;

void realloc( int nnode_ ){
    nnode =nnode_;
    ncap=nnode*nneigh_max;
    nps =nnode_+ncap; 
    nDOFs = nps*3;
    //printf( "MMFFsp3::realloc() natom(%i,nnode=%i,ncap=%i), npi=%i, nbond=%i \n", natoms, nnode, ncap, npi, nbonds );
    _realloc( DOFs  , nDOFs );
    _realloc( fDOFs , nDOFs );

    _realloc( atype, nnode );
    _realloc( aKs,   nnode );
    _realloc( a0s,   nnode );

    _realloc( neighs, ncap );  // indexes of neighbors
    _realloc( lcpos,  ncap );  // lengh of bonds to neighbors

    _realloc( l0cpos, ncap );  // bond lengh (equilibrium)
    _realloc( kcpos,  ncap );  // bond stiffness
    
    apos   = (Vec3d*)DOFs ;
    fapos  = (Vec3d*)fDOFs;
    cpos   = (Vec3d*)(DOFs +nnode);
    fcpos  = (Vec3d*)(fDOFs+nnode);
}

void cleanEnergy   (){ Etot=0;Eb=0;Ess=0;Esp=0;EppT=0;EppI=0; };
void cleanAtomForce(){ for(int i=0; i<nnode; i++){ fapos[i].set(0.0); } }
void cleanPiForce  (){ for(int i=0; i<ncap;  i++){ fcpos[i].set(0.0); } }
//void normalizePi   (){ for(int i=0; i<npi;   i++){ pipos [i].normalize(); } }

void cleanAll(){
    cleanEnergy();
    cleanAtomForce();
    //cleanPiForce();
}

// ============== Evaluation

int i_DEBUG = 0;

// interaction responsible for angle between sigma-bonds on given atom
inline double evalSS(  int ia, int ing, int jng, double K ){
    nSS++;
    Vec3d h1 = cpos[ing]; double ir1 = lcpos[ing];
    Vec3d h2 = cpos[jng]; double ir2 = lcpos[ing];
    double c = h1.dot(h2);
    Vec3d hf1,hf2;
    hf1 = h2 - h1*c;
    hf2 = h1 - h2*c;
    double c_ = c+1.0;
    double E    =  K*c_*c_;
    double fang = -K*c_*2;
    hf1.mul( fang*ir1 );
    hf2.mul( fang*ir2 );
    fapos[ing].add( hf1     );
    fapos[jng].add( hf2     );
    fapos[ia ].sub( hf1+hf2 );
    return E;
}

// interactin makes pi-orbital orthogonal on sigma bonds
inline double evalSP( int ia, int ing, int jng, double K ){
    nSP++;
    Vec3d h1 = cpos[ing]; double ir1 = lcpos[ing];
    Vec3d h2 = cpos[jng];
    double c = h1.dot(h2);
    Vec3d hf1,hf2;
    hf1 = h2 - h1*c;
    hf2 = h1 - h2*c;
    double E = K*c*c;
    double fang = -K*c*2;
    hf1.mul( fang*ir1 );
    hf2.mul( fang     );
    fcpos [ing].add( hf1 );
    fcpos [jng].add( hf2 );
    fapos [ia ].sub( hf1 );
    return E;
}

// interaction makes two p-orbitals on the same atom orthogonal
inline double evalPPT( int ia, int ing, int jng, double K ){
    nPPT++;
    Vec3d h1 = cpos[ing];
    Vec3d h2 = cpos[jng];
    double c = h1.dot(h2);
    Vec3d hf1,hf2;
    hf1 = h2 - h1*c;
    hf2 = h1 - h2*c;
    double E = K*c*c;
    double fang = -K*c*2;
    hf1.mul( fang );
    hf2.mul( fang );
    fcpos [ing].add( hf1 );
    fcpos [jng].add( hf2 );
    return E;
}

// interaction makes two p-orbitals on neighboring atoms co-linear ( e.g. planarize aromatic rings, cis-trans izomerization barrier in ethylene)
inline double evalPPI( int ing, int jng, double K ){  
    nPPI++;
    Vec3d h1 = cpos[ing];
    Vec3d h2 = cpos[jng];
    double c = h1.dot(h2);
    Vec3d hf1,hf2;
    hf1 = h2 - h1*c;
    hf2 = h1 - h2*c;
    bool sign = c<0; if(sign) c=-c;
    double E    = -K*c;
    double fang =  K;
    if(sign)fang=-fang;
    hf1.mul( fang );
    hf2.mul( fang );
    fcpos [ing].add( hf1 );
    fcpos [jng].add( hf2 );
    return E;
}

// interaction responsible for bond-length (bond strech) 
inline double eval_bond( int ia, int ja, int ing, Vec3d d, double l0, double k ){
    nB++;
    double l=d.normalize();
    lcpos[ing]=l;
     cpos[ing]=d;
    double dl = (l-l0);
    d.mul( dl*k*2 );
    fapos[ia ].add( d );
    fapos[ing].sub( d );
    return k*dl*dl;
}

double eval_neighs_bonds(int ia){
    double E=0;
    int ioff = ia*nneigh_max;
    for(int i=0; i<nneigh_max; i++){
        int ing=ioff+i;
        Vec3d  pi = apos  [ia];
        int    ja = neighs[ing];
        double l0 = l0cpos[ing];
        double k  =  kcpos[ing];
        bool bipi=ja<0;
        Vec3d  d;
        if  (ja==CAP_PI){ cpos[ing].normalize(); }  //-- pi-orbital
        else{ if(ja<0){ d=cpos[ing];             }  //-- cap  atom
              else    { d=apos[ja]-pi;           }  //-- node atom
              E+= eval_bond( ia,ja,ing, d, l0, k );
        }
    }
    return E;
}

// consider neighbors like 
// [s,s,p,p] and fact that (i<j) 
//     => angle interactions are allways like   (s-p) never (p-s) 

double eval_neighs_angles(int ia){
    double E=0;
    int ioff = ia*nneigh_max;
    const bool bPPT_ = bPPT; 
    const bool bPPI_ = bPPI;
    const bool bSS_  = bSS;
    const bool bSP_  = bSP;
    double K  = aKs[ia]; // K0
    double c0 = a0s[ia]; // angle 0
    for(int i=0; i<nneigh_max; i++){
        int   ing = ioff+i;
        int    ib = neighs[ing];
        bool bi=(ib!=CAP_PI);
        for(int j=i+1; j<nneigh_max; j++){
            int jng = ioff+j;
            int jb  = neighs[jng];
            bool bj=(jb!=CAP_PI);
            if(bi){ if( bj){ if(bSS_ ){ double e=evalSS ( ia, ing, jng, K ); Ess +=e; E+=e; } } // sigma-sigma angle
                    else   { if(bSP_ ){ double e=evalSP ( ia, ing, jng, K ); Esp +=e; E+=e; } } // sigma-pi orthogonality
            }else{  if(!bj){ if(bPPT_){ double e=evalPPT( ia, ing, jng, K ); EppT+=e; E+=e; } } // pi-pi orthogonality
                    else   { if(bPPI_){                                                         // pi-pi colinearity on neighbors
                        int joff = (jb+1)*nneigh_max;   
                        for(int ii=1; ii<3; ii++){
                            int kng=joff-ii;
                            int kb =neighs[kng];
                            if(kb!=CAP_PI) continue;
                            double e=evalPPI( ing, kng, K ); EppI  +=e; E+=e; 
                        }
                    }
                }
            }    
        }
    }
    return E;
}

double eval(){
    double E=0;
    for(int ia=0; ia<nnode; ia++){ E+=eval_neighs_bonds (ia); } // bonds
    for(int ia=0; ia<nnode; ia++){ E+=eval_neighs_angles(ia); } // neighbors/angles
    for(int i =0; i <ncap;  i++ ){ int ia=neighs[i]; if(ia>=0){ fapos[i].add(fcpos[i]); } }  // apply force back to node atoms
    return E;
}

};


#endif
