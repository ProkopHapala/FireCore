
#ifndef MMFFf4_h
#define MMFFf4_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "quaternion.h"
#include "molecular_utils.h"

// ======================
// ====   MMFFsp3
// ======================

inline float evalBond( const Vec3f& h, float dl, float k, Vec3f& f ){
    float fr = dl*k;
    f.set_mul ( h, fr );
    //f1.add( h );
    //f2.sub( h );
    return fr*dl*0.5;
}

inline float evalAngleCos( const Vec3f& h1, const Vec3f& h2, float ir1, float ir2, float K, float c0, Vec3f& f1, Vec3f& f2 ){
    float c = h1.dot(h2);
    //f1 = h2 - h1*c;
    //f2 = h1 - h2*c;
    f1.set_add_mul( h2,h1,-c );
    f2.set_add_mul( h1,h2,-c );
    float c_   = c-c0;
    float E    =  K*c_*c_;
    float fang = -K*c_*2;
    f1.mul( fang*ir1 );
    f2.mul( fang*ir2 );
    return E;
}

inline float evalPiAling( const Vec3f& h1, const Vec3f& h2, float ir1, float ir2, float K, Vec3f& f1, Vec3f& f2 ){  // interaction between two pi-bonds
    float c = h1.dot(h2);
    //f1 = h2 - h1*c;
    //f2 = h1 - h2*c;
    f1.set_add_mul( h2,h1,-c );
    f2.set_add_mul( h1,h2,-c );
    bool sign = c<0; if(sign) c=-c;
    float E    = -K*c;
    float fang =  K;
    if(sign)fang=-fang;
    f1.mul( fang );
    f2.mul( fang );
    return E;
}

class MMFFf4{ public:
    static constexpr const int nneigh_max = 4;
    int  nDOFs=0,natoms=0,nnode=0,ncap=0,nvecs=0;
    bool bPBC=false;
    float Etot,Eb,Ea, Eps,EppT,EppI;

    float *  DOFs = 0;   // degrees of freedom
    float * fDOFs = 0;   // forces
    float * vDOFs = 0;   // velocities
    
    //                           c0     Kss    Ksp    c0_e
    Quat4d default_NeighParams{ -1.0,   1.0,   1.0,   -1.0 };

    // Dynamical Varaibles;
    Quat4f *   apos=0;   // [natom]
    Quat4f *  fapos=0;   // [natom]
    Quat4f *  pipos=0;   // [nnode]
    Quat4f * fpipos=0;   // [nnode]
    // Aux Dynamil
    Quat4f * fneigh  =0;  // [nnode*4]     temporary store of forces on atoms form neighbors (before assembling step)
    Quat4f * fneighpi=0;  // [nnode*4]     temporary store of forces on pi    form neighbors (before assembling step)

    // Params
    Quat4i*  aneighs =0; // [nnode*4]   index of neighboring atoms
    Quat4i*  bkneighs=0; // [natoms*4]  inverse neighbors

    Quat4f*  apars=0;  // [nnode] per atom forcefield parametrs
    Quat4f*  REQs =0;  // [nnode] parameters of non-covalent interactions
    Quat4f*  bLs  =0;  // [nnode] bond lengths
    Quat4f*  bKs  =0;  // [nnode] bond stiffness
    Quat4f*  Ksp  =0;  // [nnode] stiffness of pi-alignment
    Quat4f*  Kpp  =0;  // [nnode] stiffness of pi-planarization


    bool    bSubtractAngleNonBond=false;
    bool    bPBCbyLvec  =false;
    Mat3d   invLvec, lvec;

// =========================== Functions

void realloc( int nnode_, int ncap_ ){
    nnode=nnode_; ncap=ncap_;
    natoms= nnode  + ncap; 
    nvecs = natoms+nnode;  // each atom as also pi-orientiation (like up-vector)
    nDOFs = nvecs*4;
    //printf( "MMFFsp3::realloc() natom(%i,nnode=%i,ncap=%i), npi=%i, nbond=%i \n", natoms, nnode, ncap, npi, nbonds );
    int ipi0=natoms;

    // ----- Dynamical
    _realloc( DOFs     , nDOFs );
    _realloc( fDOFs    , nDOFs );
    _realloc( vDOFs    , nDOFs );
    apos   = (Quat4f*) DOFs ;
    fapos  = (Quat4f*)fDOFs;
    pipos  = apos  + ipi0;
    fpipos = fapos + ipi0;
    // ---- Aux
    _realloc( fneigh  , nnode*4 );
    _realloc( fneighpi, nnode*4 );
    // ----- Params [natom]
    _realloc( aneighs , nnode  );
    _realloc( bkneighs, natoms );
    _realloc( apars  , nnode );
    _realloc( bLs    , nnode );
    _realloc( bKs    , nnode );
    _realloc( Ksp    , nnode );
    _realloc( Kpp    , nnode );


}

// ============== Evaluation

float eval_atom(int ia){

    float E=0;
    const Vec3f pa  = apos [ia].f; 
    const Vec3f hpi = pipos[ia].f; 

    //Vec3f& fa  = fapos [ia].f; 
    //Vec3f& fpi = fpipos[ia].f; 
    Vec3f fa  = Vec3fZero; 
    Vec3f fpi = Vec3fZero; 
    
    //--- array aliases
    const int*   ings = aneighs[ia].array;
    const float* bK   = bKs    [ia].array;
    const float* bL   = bLs    [ia].array;
    const float* Kspi = Ksp    [ia].array;
    const float* Kppi = Kpp    [ia].array;
    Quat4f* fbs  = fneigh   +ia*4;
    Quat4f* fps  = fneighpi +ia*4;

    // --- settings
    float  ssC0 = apars[ia].x;
    float  ssK  = apars[ia].y;
    float  piC0 = apars[ia].z;
    bool   bPi  = ings[3]<0;

    //--- Aux Variables 
    Quat4f hs[4];
    Vec3f  f1,f2;

    // --------- Bonds Step
    for(int i=0; i<4; i++){
 int ing = ings[i];
        //printf( "bond[%i|%i=%i]\n", ia,i,ing );
        //fbs[i]=Quat4fOnes; fps[i]=Quat4fOnes;
        fbs[i]=Quat4fZero; fps[i]=Quat4fZero;
        if(ing<0) break;
        Quat4f h; 
        h.f.set_sub( apos[ing].f, pa );
        //if(idebug)printf( "bond[%i|%i=%i] l=%g pj[%i](%g,%g,%g) pi[%i](%g,%g,%g)\n", ia,i,ing, h.f.norm(), ing,apos[ing].x,apos[ing].y,apos[ing].z, ia,pa.x,pa.y,pa.z  );
        float  l = h.f.normalize();
        h.e      = 1/l;
        hs [i]   = h;

        if(ia<ing){   // we should avoid double counting because otherwise node atoms would be computed 2x, but capping only once

            //E+= evalBond( h.f, l-bL[i], bK[i], f1 ); fbs[i].add(f1);  fa.sub(f1);    // bond length force
            E+= evalBond( h.f, l-bL[i], bK[i], f1 );  fbs[i].f.sub(f1);  fa.add(f1);    

            double kpp = Kppi[i];
            if( (ing<nnode) && (kpp>1e-6) ){   // Only node atoms have pi-pi alignemnt interaction
                E += evalPiAling( hpi, pipos[ing].f, 1., 1.,   kpp,       f1, f2 );   fpi.add(f1);  fps[i].f.add(f2);    //   pi-alignment     (konjugation)
                if(idebug)printf( "pi-pi[%i|%i] kpp=%g c=%g l(%g,%g) f1(%g,%g,%g) f2(%g,%g,%g) \n", ia,ing, kpp, hpi.dot(pipos[ing].f),1.,1., f1.x,f1.y,f1.z,  f2.x,f2.y,f2.z  );
            }
            // ToDo: triple bonds ?

        } 

        // pi-sigma 
        //if(bPi){    
        double ksp = Kspi[i];
        if(ksp>1e-6){  
            E += evalAngleCos( hpi, h.f      , 1., h.e, ksp, piC0, f1, f2 );   fpi.add(f1);  fbs[i].f.add(f2);    //   pi-planarization (orthogonality)
            if(idebug)printf( "pi-sigma[%i|%i] ksp=%g c=%g l(%g,%g) f1(%g,%g,%g) f2(%g,%g,%g) \n", ia,ing, ksp, hpi.dot(h.f),1.,h.e, f1.x,f1.y,f1.z,  f2.x,f2.y,f2.z  );
        }
        //}

    }

    // --------- Angle Step
    for(int i=0; i<4; i++){
        int ing = ings[i];
        if(ing<0) break;
        const Quat4f& hi = hs[i];
        for(int j=i+1; j<4; j++){
            int jng  = ings[j];
            if(jng<0) break;
            const Quat4f& hj = hs[j];
            E += evalAngleCos( hi.f, hj.f, hi.e, hj.e, ssK, ssC0, f1, f2 );     // angles between sigma bonds
            //if(idebug)printf( "ang[%i|%i,%i] kss=%g c=%g l(%g,%g) f1(%g,%g,%g) f2(%g,%g,%g)\n", ia,ing,jng, ssK, hi.f.dot(hj.f),hi.e,hj.e, f1.x,f1.y,f1.z,  f2.x,f2.y,f2.z  );
            fbs[i].f.add( f1     );
            fbs[j].f.add( f2     );
            fa.      sub( f1+f2  );
            // ToDo: subtract non-covalent interactions
        }
    }

    //fapos [ia].add(fa ); 
    //fpipos[ia].add(fpi);
    fapos [ia].f=fa; 
    fpipos[ia].f=fpi;
    return E;
}

double eval_atoms(){
    double E=0;
    for(int ia=0; ia<nnode; ia++){ E+=eval_atom(ia); }
    return E;
}

void normalizePis(){ 
    for(int i=0; i<nnode; i++){ pipos[i].f.normalize(); } 
}

void cleanForce(){ 
    for(int i=0; i<nDOFs; i++){ fDOFs[i]=0;  } 
}

void asseble_forces(){
    for(int ia=0; ia<natoms; ia++){
        Quat4f fa=Quat4fZero,fp=Quat4fZero;
        const int* ings = bkneighs[ia].array;
        bool bpi = ia<nnode;
        for(int i=0; i<4; i++){
            int j = ings[i];
            if(j<0) break;
            //if(j>=(nnode*4)){ printf("ERROR bkngs[%i|%i] %i>=4*nnode(%i)\n", ia, i, j, nnode*4 ); exit(0); }
            fa.add(fneigh  [j]);
            if(bpi)fp.add(fneighpi[j]);
        }
        fa.e=0;
        fp.e=0;
        fapos [ia].add( fa ); 
        if(bpi){
            fpipos[ia].add( fp );
            fpipos[ia].f.makeOrthoU( pipos[ia].f );  // subtract force component which change pi-vector size
        }
    }
}

float eval( bool bClean=true, bool bCheck=true ){
    //if(bClean){ cleanAll(); }
    //printf( "print_apos() BEFORE\n" );print_apos();
    cleanForce();
    normalizePis();
    //printf( "print_apos() AFTER \n" ); print_apos();
    eval_atoms();
    asseble_forces();
    //Etot = Eb + Ea + Eps + EppT + EppI;
    return Etot;
}


void move_GD(float dt){
    //for(int i=0; i<0; i++){}
    //Quat4f vs = (Quat4f*) vDOFs; 
    for(int i=0; i<nvecs; i++){
        //vs[i].f.mul( f, dt );
        apos[i].f.add_mul( fapos[i].f, dt ); 
    }
}

void makeBackNeighs( ){
    for(int i=0; i<natoms; i++){ bkneighs[i]=Quat4i{-1,-1,-1,-1}; };
    for(int ia=0; ia<nnode; ia++){
        for(int j=0; j<4; j++){        // 4 neighbors
            int ja = aneighs[ia].array[j];
            if( ja<0 )continue;
            //NOTE: We deliberately ignore back-neighbors from caping atoms 
            bool ret = addFirstEmpty( bkneighs[ja].array, 4, ia*4+j, -1 );
            if(!ret){ printf("ERROR in MMFFf4_loc::makeBackNeighs(): Atom #%i has >4 back-Neighbors (while adding atom #%i) \n", ja, ia ); exit(0); }
        };
    }
}

void printAtomParams(int ia){ printf("atom[%i] ngs{%3i,%3i,%3i,%3i} par(%5.3f,%5.3f,%5.3f)  bL(%5.3f,%5.3f,%5.3f,%5.3f) bK(%6.3f,%6.3f,%6.3f,%6.3f)  Ksp(%5.3f,%5.3f,%5.3f,%5.3f) Kpp(%5.3f,%5.3f,%5.3f,%5.3f) \n", ia, aneighs[ia].x,aneighs[ia].y,aneighs[ia].z,aneighs[ia].w,    apars[ia].x,apars[ia].y,apars[ia].z,    bLs[ia].x,bLs[ia].y,bLs[ia].z,bLs[ia].w,   bKs[ia].x,bKs[ia].y,bKs[ia].z,bKs[ia].w,     Ksp[ia].x,Ksp[ia].y,Ksp[ia].z,Ksp[ia].w,   Kpp[ia].x,Kpp[ia].y,Kpp[ia].z,Kpp[ia].w  ); };
void printAtomParams(){for(int ia=0; ia<nnode; ia++){ printAtomParams(ia); }; };

void printBKneighs(int ia){ printf("atom[%i] bkngs{%3i,%3i,%3i,%3i} \n", ia, bkneighs[ia].x,bkneighs[ia].y,bkneighs[ia].z,bkneighs[ia].w ); };
void printBKneighs(){for(int ia=0; ia<natoms; ia++){ printBKneighs(ia); }; };

void print_apos(){
    for(int ia=0;ia<natoms;ia++){ printf( "print_apos[%i](%g,%g,%g)\n", ia, apos[ia].x,apos[ia].y,apos[ia].z ); }
}

};


#endif
