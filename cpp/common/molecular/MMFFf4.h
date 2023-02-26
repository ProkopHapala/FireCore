
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
    f1.set_add_mul( h2,h1,c );
    f2.set_add_mul( h1,h2,c );
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
    f1.set_add_mul( h2,h1,c );
    f2.set_add_mul( h1,h2,c );
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
    
    //                           c0     Kss    Ksp    c0_e
    Quat4d default_NeighParams{ -1.0,   1.0,   1.0,   -1.0 };

    // Dynamical Varaibles;
    Quat4f *   apos=0;   // [natom]
    Quat4f *  fapos=0;   // [natom]
    Quat4f *  pipos=0;   // [nnode]
    Quat4f * fpipos=0;   // [nnode]
    // Aux Dynamil
    Quat4f * fneih  =0;  // [nnode*4]     temporary store of forces on atoms form neighbors (before assembling step)
    Quat4f * fneihpi=0;  // [nnode*4]     temporary store of forces on pi    form neighbors (before assembling step)

    // Params
    Quat4i*  aneighs=0; // [nnode*4]   index of neighboring atoms
    Quat4i*  ineighs=0; // [natoms*4]  inverse neighbors

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
    apos   = (Quat4f*) DOFs ;
    fapos  = (Quat4f*)fDOFs;
    pipos  = apos  + ipi0;
    fpipos = fapos + ipi0;
    // ---- Aux
    _realloc( fneih  , nnode*4 );
    _realloc( fneihpi, nnode*4 );
    // ----- Params [natom]
    _realloc( aneighs, nnode );
    _realloc( ineighs, natoms );
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
    const float* Kppi = Ksp    [ia].array;
    const float* Kspi = Kpp    [ia].array;
    Quat4f* fbs  = fneih   +ia*4;
    Quat4f* fps  = fneihpi +ia*4;

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
        if(ing<0) break;
        Quat4f h; 
        h.f.set_sub( apos[ing].f, pa );
        float l = h.f.normalize();
        h.e    = 1/l;
        hs [i] = h;
        // bond length force
        E+= evalBond( h.f, l-bL[i], bK[i], f1 ); fbs[i].f.set(f1);  fa.sub(f1);      // bond length force
        // pi-force
        if(bPi){    
            double ksp = Kspi[i];
            if(ksp>1e-6)  
            E += evalAngleCos( hpi, h.f         , 1., h.e, ksp, piC0, f1, f2 );   fpi.add(f1);  fbs[i].f.add(f2);    //   pi-planarization (orthogonality)
            double kpp = Kppi[i];
            if(kpp>1e-6)
            E += evalPiAling( hpi, pipos[ing].f, 1., 1,   kpp,       f1, f2 );   fpi.add(f1);  fps[i].f.add(f2);    //   pi-alignment     (konjugation)
            // ToDo: triple bonds ?
        }
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

float eval_atoms(){
    float E=0;
    for(int ia=0; ia<nnode; ia++){ E+=eval_atom(ia); }
    return E;
}

void asseble_forces(){
    for(int ia=0; ia<nnode; ia++){
        int io=ia*4;
        Quat4f fa=Quat4fZero,fp=Quat4fZero;
        const int* ings = ineighs[ia].array;
        for(int i=0; i<4; i++){
            int j = ings[i];
            if(j<0) break;
            fa.add(fneih  [j]);
            fp.add(fneihpi[j]);
        }
        fapos [ia].add( fa ); 
        fpipos[ia].add( fp );
    }
}

float eval( bool bClean=true, bool bCheck=true ){
    //if(bClean){ cleanAll(); }
    //eval_atoms();
    asseble_forces();
    //Etot = Eb + Ea + Eps + EppT + EppI;
    return Etot;
}


};


#endif
