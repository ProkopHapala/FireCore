
#ifndef MMFFsp3_loc_h
#define MMFFsp3_loc_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Forces.h"
#include "quaternion.h"
#include "molecular_utils.h"

// ======================
// ====   MMFFsp3
// ======================

class MMFFsp3_loc{ public:
    static constexpr const int nneigh_max = 4;
    int  nDOFs=0,natoms=0,nnode=0,ncap=0,nvecs=0;
    bool bPBC=false;
    double Etot,Eb,Ea, Eps,EppT,EppI;

    double *  DOFs = 0;   // degrees of freedom
    double * fDOFs = 0;   // forces
    
    //                           c0     Kss    Ksp    c0_e
    Quat4d default_NeighParams{ -1.0,   1.0,   1.0,   -1.0 };

    // Dynamical Varaibles;
    Vec3d *   apos=0;   // [natom]
    Vec3d *  fapos=0;   // [natom]
    Vec3d *  pipos=0;   // [nnode]
    Vec3d * fpipos=0;   // [nnode]
    // Aux Dynamil
    Vec3d * fneih  =0;  // [nnode*4]     temporary store of forces on atoms form neighbors (before assembling step)
    Vec3d * fneihpi=0;  // [nnode*4]     temporary store of forces on pi    form neighbors (before assembling step)

    // Params
    Quat4i*  aneighs =0; // [nnode*4]   index of neighboring atoms
    Quat4i*  bkneighs=0; // [natoms*4]  inverse neighbors

    Quat4d*  apars=0;  // [nnode] per atom forcefield parametrs
    Quat4d*  bLs  =0;  // [nnode] bond lengths
    Quat4d*  bKs  =0;  // [nnode] bond stiffness
    Quat4d*  Ksp  =0;  // [nnode] stiffness of pi-alignment
    Quat4d*  Kpp  =0;  // [nnode] stiffness of pi-planarization
    Vec3d*  REQs =0;   // [nnode] parameters of non-covalent interactions

    bool    bSubtractAngleNonBond=false;
    Mat3d   invLvec, lvec;

// =========================== Functions

void realloc( int nnode_, int ncap_ ){
    nnode=nnode_; ncap=ncap_;
    natoms= nnode + ncap; 
    nvecs = natoms+nnode;  // each atom as also pi-orientiation (like up-vector)
    nDOFs = nvecs*3;
    //printf( "MMFFsp3::realloc() natom(%i,nnode=%i,ncap=%i), npi=%i, nbond=%i \n", natoms, nnode, ncap, npi, nbonds );
    int ipi0=natoms;

    // ----- Dynamical
    _realloc(  DOFs    , nDOFs );
    _realloc( fDOFs    , nDOFs );
    apos   = (Vec3d*) DOFs ;
    fapos  = (Vec3d*)fDOFs;
    pipos  = apos  + ipi0;
    fpipos = fapos + ipi0;
    // ---- Aux
    _realloc( fneih  , nnode*4 );
    _realloc( fneihpi, nnode*4 );
    // ----- Params [natom]
    _realloc( aneighs , nnode );
    _realloc( bkneighs, natoms );
    _realloc( apars   , nnode );
    _realloc( bLs     , nnode );
    _realloc( bKs     , nnode );
    _realloc( Ksp     , nnode );
    _realloc( Kpp     , nnode );


}

// ============== Evaluation

double eval_atom(const int ia){
    //printf( "MMFFsp3_loc::eval_atom(%i)\n", ia );
    double E=0;
    const Vec3d pa  = apos [ia]; 
    const Vec3d hpi = pipos[ia]; 

    //printf( "apos[%i](%g,%g,%g)\n",ia,apos[ia].x,apos[ia].y,apos[ia].z );
    //return E;

    //Vec3d& fa  = fapos [ia]; 
    //Vec3d& fpi = fpipos[ia];
    Vec3d fa   = Vec3dZero;
    Vec3d fpi  = Vec3dZero; 
    
    //--- array aliases
    const int*    ings = aneighs[ia].array;
    const double* bK   = bKs    [ia].array;
    const double* bL   = bLs    [ia].array;
    const double* Kspi = Ksp    [ia].array;
    const double* Kppi = Kpp    [ia].array;
    Vec3d* fbs  = fneih   +ia*4;
    Vec3d* fps  = fneihpi +ia*4;

    // --- settings
    double  ssC0 = apars[ia].x;
    double  ssK  = apars[ia].y;
    double  piC0 = apars[ia].z;
    //bool    bPi  = ings[3]<0;   we distinguish this by Ksp, otherwise it would be difficult for electron pairs e.g. (-O-C=)

    //--- Aux Variables 
    Quat4d  hs[4];
    Vec3d   f1,f2;

    //if(idebug)printf( "atom[%i] Ksp{%5.3f,%5.3f,%5.3f,%5.3f} \n", ia, Kspi[0],Kspi[1],Kspi[2],Kspi[3] );
    // --------- Bonds Step
    for(int i=0; i<4; i++){
        int ing = ings[i];
        //printf( "bond[%i|%i=%i]\n", ia,i,ing );
        //fbs[i]=Vec3dOne; fps[i]=Vec3dOne;
        fbs[i]=Vec3dZero; fps[i]=Vec3dZero;
        if(ing<0) break;
        Quat4d h; 
        h.f.set_sub( apos[ing], pa );
        //if(idebug)printf( "bond[%i|%i=%i] l=%g pj[%i](%g,%g,%g) pi[%i](%g,%g,%g)\n", ia,i,ing, h.f.norm(), ing,apos[ing].x,apos[ing].y,apos[ing].z, ia,pa.x,pa.y,pa.z  );
        double l = h.f.normalize();
        h.e    = 1/l;
        hs [i] = h;
        // bond length force
        //continue; 
        if(ia<ing){   // we should avoid double counting because otherwise node atoms would be computed 2x, but capping only once
            //E+= evalBond( h.f, l-bL[i], bK[i], f1 ); fbs[i].add(f1);  fa.sub(f1);    // bond length force
            E+= evalBond( h.f, l-bL[i], bK[i], f1 ); 
            fbs[i].sub(f1);  
            fa.add(f1);    
            //if(idebug)printf( "bond[%i|%i=%i] l %g dl0,k(%g,%g) h(%g,%g,%g) f(%g,%g,%g)\n", ia,i,ing, l,bL[i], bK[i], h.f.x,h.f.y,h.f.z, f1.x,f1.y,f1.z  );

            double kpp = Kppi[i];
            if( (ing<nnode) && (kpp>1e-6) ){   // Only node atoms have pi-pi alignemnt interaction
                E += evalPiAling( hpi, pipos[ing], 1., 1,   kpp,       f1, f2 );   fpi.add(f1);  fps[i].add(f2);    //   pi-alignment     (konjugation)
            }
            // ToDo: triple bonds ?
        } 
        
        // pi-sigma 
        //if(bPi){    
        double ksp = Kspi[i];
        if(ksp>1e-6){  
            E += evalAngleCos( hpi, h.f      , 1., h.e, ksp, piC0, f1, f2 );   fpi.add(f1);  fbs[i].add(f2);    //   pi-planarization (orthogonality)
            //if(idebug)printf( "pi-sigma[%i|%i=%i] ksp=%g e=%g \n", ia,i,ing, ksp, e  );
        }
        //}
        
    }
    
    // --------- Angle Step
    for(int i=0; i<4; i++){
        int ing = ings[i];
        if(ing<0) break;
        const Quat4d& hi = hs[i];
        for(int j=i+1; j<4; j++){
            int jng  = ings[j];
            if(jng<0) break;
            const Quat4d& hj = hs[j];
            E += evalAngleCos( hi.f, hj.f, hi.e, hj.e, ssK, ssC0, f1, f2 );     // angles between sigma bonds
            fbs[i].add( f1     );
            fbs[j].add( f2     );
            fa    .sub( f1+f2  );
            // ToDo: subtract non-covalent interactions
        }
    }
    
    //fapos [ia].add(fa ); 
    //fpipos[ia].add(fpi);
    fapos [ia]=fa; 
    fpipos[ia]=fpi;
    return E;
}

double eval_atoms(){
    double E=0;
    for(int ia=0; ia<nnode; ia++){ E+=eval_atom(ia); }
    return E;
}

void normalizePis(){ 
    for(int i=0; i<nnode; i++){ pipos[i].normalize(); } 
}

void cleanForce(){ 
    //for(int i=0; i<natoms; i++){ fapos [i].set(0.0);  } 
    //for(int i=0; i<nnode;  i++){ fpipos[i].set(0.0);  } 
    for(int i=0; i<nDOFs; i++){ fDOFs[i]=0;  } 
}

void asseble_forces(){
    for(int ia=0; ia<natoms; ia++){
        Vec3d fa=Vec3dZero,fp=Vec3dZero;
        const int* ings = bkneighs[ia].array;
        bool bpi = ia<nnode;
        for(int i=0; i<4; i++){
            int j = ings[i];
            if(j<0) break;
            //if(j>=(nnode*4)){ printf("ERROR bkngs[%i|%i] %i>=4*nnode(%i)\n", ia, i, j, nnode*4 ); exit(0); }
            fa.add(fneih  [j]);
            if(bpi)fp.add(fneihpi[j]);
        }
        fapos [ia].add( fa ); 
        if(bpi){
            fpipos[ia].add( fp );
            fpipos[ia].makeOrthoU( pipos[ia] );  // subtract force component which change pi-vector size
        }
    }
}

double eval( bool bClean=true, bool bCheck=true ){
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

void makeBackNeighs( ){
    for(int i=0; i<natoms; i++){ bkneighs[i]=Quat4i{-1,-1,-1,-1}; };
    for(int ia=0; ia<nnode; ia++){
        for(int j=0; j<4; j++){        // 4 neighbors
            int ja = aneighs[ia].array[j];
            if( ja<0 )continue;
            //NOTE: We deliberately ignore back-neighbors from caping atoms 
            bool ret = addFirstEmpty( bkneighs[ja].array, 4, ia*4+j, -1 );
            if(!ret){ printf("ERROR in MMFFsp3_loc::makeBackNeighs(): Atom #%i has >4 back-Neighbors (while adding atom #%i) \n", ja, ia ); exit(0); }
        };
    }
    //for(int i=0; i<natoms; i++){printf( "bkneigh[%i] (%i,%i,%i,%i) \n", i, bkneighs[i].x, bkneighs[i].y, bkneighs[i].z, bkneighs[i].w );}
    //checkBkNeighCPU();
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
