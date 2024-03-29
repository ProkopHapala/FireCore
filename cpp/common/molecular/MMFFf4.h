
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

inline void combineREQ(const Quat4f& a, const Quat4f& b, Quat4f& out){
    out.x=a.x+b.x; // radius
    out.y=a.y*b.y; // epsilon
    out.z=a.z*b.z; // q*q
    out.w=out.w*out.w;
}

inline float evalAngleCosHalf( const Vec3f& h1, const Vec3f& h2, float ir1, float ir2, const Vec2f& cs0, float k, Vec3f& f1, Vec3f& f2 ){
    //printf( " ir1 %g ir2 %g \n", ir1, ir2 );
    // This is much better angular function than evalAngleCos() with just a little higher computational cost ( 2x sqrt )
    Vec3f h; h.set_add( h1, h2 );
    float c2 = h.norm2()*0.25;     // cos(a/2) = |ha+hb|
    float s2 = 1-c2 + 1e-7;        // s2 must be positive number
    float c  = sqrt(c2);
    float s  = sqrt(s2);
    Vec2f cs  = cs0;
    cs.udiv_cmplx({c,s});
    //Vec2f cs{c,s};
    //cs.mul_cmplx(cs0);
    float E         =  k*( 1 - cs.x );  // just for debug ?
    float fr        = -k*(     cs.y );
    c2 *=-2;
    fr /= 4*c*s;   //    |h - 2*c2*a| =  1/(2*s*c) = 1/sin(a)
    float fr1    = fr*ir1;
    float fr2    = fr*ir2;
    f1.set_lincomb( fr1, h,  fr1*c2, h1 );  //fa = (h - 2*c2*a)*fr / ( la* |h - 2*c2*a| );
    f2.set_lincomb( fr2, h,  fr2*c2, h2 );  //fb = (h - 2*c2*b)*fr / ( lb* |h - 2*c2*b| );
    return E;
}

inline float addAtomicForceLJQ( const Vec3f& dp, Vec3f& f, const Quat4f& REQ ){
    //Vec3f dp; dp.set_sub( p2, p1 );
    const float COULOMB_CONST_ = 14.3996448915;  //  [V*A/e] = [ (eV/A) * A^2 /e^2]
    float ir2  = 1/( dp.norm2() + 1e-4 );
    float ir   = sqrt(ir2);
    float ir2_ = ir2*REQ.x*REQ.x;
    float ir6  = ir2_*ir2_*ir2_;
    //float fr   = ( ( 1 - ir6 )*ir6*12*REQ.b + ir*REQ.c*-COULOMB_CONST )*ir2;
    float Eel  = ir*REQ.z*COULOMB_CONST_;
    float vdW  = ir6*REQ.y;
    float fr   = ( ( 1 - ir6 )*12*vdW - Eel )*ir2;
    //printf( " (%g,%g,%g) r %g fr %g \n", dp.x,dp.y,dp.z, 1/ir, fr );
    f.add_mul( dp, fr );
    return  ( ir6 - 2 )*vdW + Eel;
}

//==================================================================
//    class      MMFFf4
///================================================================

class MMFFf4{ public:
    static constexpr const int nneigh_max = 4;
    int  nDOFs=0,natoms=0,nnode=0,ncap=0,nvecs=0;
    bool bPBC=false;
    float Etot,Eb,Ea, Eps,EppT,EppI;

    float *  DOFs = 0;   // degrees of freedom
    float * fDOFs = 0;   // forces
    float * vDOFs = 0;   // velocities
    
    //                            c0     Kss    Ksp     c0_e
    Quat4f default_NeighParams{ -1.0f,   1.0f,   1.0f, -1.0f };

    // Dynamical Varaibles;
    Quat4f *   apos=0;   // [natom]
    Quat4f *  fapos=0;   // [natom]
    Quat4f *  pipos=0;   // [nnode]
    Quat4f * fpipos=0;   // [nnode]
    // Aux Dynamil
    Quat4f * fneigh  =0;  // [nnode*4]     temporary store of forces on atoms form neighbors (before assembling step)
    Quat4f * fneighpi=0;  // [nnode*4]     temporary store of forces on pi    form neighbors (before assembling step)

    // Params
    Quat4i*  neighs     =0; // [natoms]   index of neighboring atoms
    Quat4i*  neighCell  =0; // [natoms]   index of neighboring atoms
    Quat4i*  bkneighs    =0; // [natoms]  inverse neighbors

    Quat4f*  apars=0;  // [nnode] per atom forcefield parametrs
    Quat4f*  REQs =0;  // [nnode] parameters of non-covalent interactions
    Quat4f*  bLs  =0;  // [nnode] bond lengths
    Quat4f*  bKs  =0;  // [nnode] bond stiffness
    Quat4f*  Ksp  =0;  // [nnode] stiffness of pi-alignment
    Quat4f*  Kpp  =0;  // [nnode] stiffness of pi-planarization

    bool    bAngleCosHalf         = true;
    bool    bSubtractAngleNonBond = false;
    bool    bPBCbyLvec  =false;
    Mat3f   invLvec, lvec;
    Vec3i   nPBC;
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
    //_realloc( neighs , nnode  );   // We need neighs for all atoms because of Non-Bonded
    _realloc( neighs   , natoms );   // We need neighs for all atoms because of Non-Bonded
    _realloc( neighCell, natoms );   // We need neighs for all atoms because of Non-Bonded
    _realloc( bkneighs  , natoms );
    _realloc( apars  , nnode );
    _realloc( bLs    , nnode );
    _realloc( bKs    , nnode );
    _realloc( Ksp    , nnode );
    _realloc( Kpp    , nnode );

}

void dealloc(){
    _dealloc(DOFs );
    _dealloc(fDOFs);
    _dealloc(vDOFs);
    apos   = 0;
    fapos  = 0;
    pipos  = 0;
    fpipos = 0;
    _dealloc(neighs);
    _dealloc(neighCell);
    _dealloc(bkneighs);
    _dealloc(apars);
    _dealloc(bLs);
    _dealloc(bKs);
    _dealloc(Ksp);
    _dealloc(Kpp);
    //_dealloc(angles);
    //_dealloc(tors2atom  );
    //_dealloc(torsParams );
    nnode=0; ncap=0; natoms=0; nvecs=0; nDOFs=0;
}

void setLvec(const Mat3f& lvec_){ lvec=lvec_; lvec.invert_T_to( invLvec ); }

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
    const int*   ings = neighs[ia].array;
    const float* bK   = bKs   [ia].array;
    const float* bL   = bLs   [ia].array;
    const float* Kspi = Ksp   [ia].array;
    const float* Kppi = Kpp   [ia].array;
    Quat4f* fbs  = fneigh   +ia*4;
    Quat4f* fps  = fneighpi +ia*4;

    // --- settings
    // float  ssC0 = apars[ia].x;
    // float  ssK  = apars[ia].y;
    // float  piC0 = apars[ia].z;

    const Quat4f& apar   = apars[ia];
    const float   ssK    = apar.z;
    const float   piC0   = apar.w;
    const Vec2f   cs0_ss = Vec2f{apar.x,apar.y};
    const float   ssC0   = cs0_ss.x*cs0_ss.x - cs0_ss.y*cs0_ss.y;   // cos(2x) = cos(x)^2 - sin(x)^2, because we store cos(ang0/2) to use in  evalAngleCosHalf

    bool   bPi  = ings[3]<0;

    //--- Aux Variables 
    Quat4f hs[4];
    Vec3f  f1,f2;

    //const int iG_DBG = 0;
    const int ia_DBG = 1;
    //if(ia==ia_DBG)printf( "CPU[%i] neighs(%i,%i,%i,%i) \n", ia, ings[0],ings[1],ings[2],ings[3] );

    for(int i=0; i<4; i++){ fbs[i]=Quat4fZero; fps[i]=Quat4fZero; } // we initialize it here because of the break

    // --------- Bonds Step
    for(int i=0; i<4; i++){
        int ing = ings[i];
        //printf( "bond[%i|%i=%i]\n", ia,i,ing );
        fbs[i]=Quat4fZero; fps[i]=Quat4fZero;
        //fbs[i]=Quat4fOnes; fps[i]=Quat4fOnes;
        //fbs[i]=Quat4f{1,0,1,0}; fps[i]=Quat4f{1,2,1,2};
        if(ing<0) break;
        Quat4f h; 
        h.f.set_sub( apos[ing].f, pa );
        //if(idebug)printf( "bond[%i|%i=%i] l=%g pj[%i](%g,%g,%g) pi[%i](%g,%g,%g)\n", ia,i,ing, h.f.norm(), ing,apos[ing].x,apos[ing].y,apos[ing].z, ia,pa.x,pa.y,pa.z  );
        float  l = h.f.normalize();
        h.e      = 1/l;
        hs [i]   = h;

        if(ia<ing){   // we should avoid double counting because otherwise node atoms would be computed 2x, but capping only once
            E+= evalBond( h.f, l-bL[i], bK[i], f1 );  fbs[i].f.sub(f1);  fa.add(f1);    
            //if(ia==ia_DBG)printf( "CPU bond[%i|%i] kpb=%g l0=%g l=%g h(%g,%g,%g) f(%g,%g,%g) \n", ia,ing, bK[i],bL[i], l, h.x,h.y,h.z,  f1.x,f1.y,f1.z  );
            
            float kpp = Kppi[i];
            if( (ing<nnode) && (kpp>1e-6) ){   // Only node atoms have pi-pi alignemnt interaction
                E += evalPiAling( hpi, pipos[ing].f, 1., 1.,   kpp,       f1, f2 );   fpi.add(f1);  fps[i].f.add(f2);    //   pi-alignment     (konjugation)
                //if(ia==ia_DBG)printf( "CPU:pipi[%i|%i] kpp=%g c=%g f1(%g,%g,%g) f2(%g,%g,%g)\n", ia,ing, kpp, hpi.dot(pipos[ing].f), f1.x,f1.y,f1.z,  f2.x,f2.y,f2.z  );
                //if(idebug)printf( "pi-pi[%i|%i] kpp=%g c=%g l(%g,%g) f1(%g,%g,%g) f2(%g,%g,%g) \n", ia,ing, kpp, hpi.dot(pipos[ing].f),1.,1., f1.x,f1.y,f1.z,  f2.x,f2.y,f2.z  );
            }
            // ToDo: triple bonds ?
            
        } 
        
        // pi-sigma 
        //if(bPi){    
        float ksp = Kspi[i];
        if(ksp>1e-6){  
            E += evalAngleCos( hpi, h.f      , 1., h.e, ksp, piC0, f1, f2 );   fpi.add(f1);  fa.sub(f2); fbs[i].f.add(f2);    //   pi-planarization (orthogonality)
            //if(idebug)printf( "pi-sigma[%i|%i] ksp=%g c=%g l(%g,%g) f1(%g,%g,%g) f2(%g,%g,%g) \n", ia,ing, ksp, hpi.dot(h.f),1.,h.e, f1.x,f1.y,f1.z,  f2.x,f2.y,f2.z  );
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
            if( bAngleCosHalf ){
                E += evalAngleCosHalf( hi.f, hj.f,  hi.e, hj.e,  cs0_ss,  ssK, f1, f2 );
            }else{             
                E += evalAngleCos( hi.f, hj.f, hi.e, hj.e, ssK, ssC0, f1, f2 );     // angles between sigma bonds
            }
            //if(idebug)printf( "ang[%i|%i,%i] kss=%g c=%g l(%g,%g) f1(%g,%g,%g) f2(%g,%g,%g)\n", ia,ing,jng, ssK, hi.f.dot(hj.f),hi.e,hj.e, f1.x,f1.y,f1.z,  f2.x,f2.y,f2.z  );
            //if(ia==ia_DBG)printf( "CPU:ang[%i|%i,%i] kss=%g c0=%g c=%g l(%g,%g) f1(%g,%g,%g) f2(%g,%g,%g)\n", ia,ing,jng, ssK, ssC0, hi.f.dot(hj.f),hi.e,hj.e, f1.x,f1.y,f1.z,  f2.x,f2.y,f2.z  );
            fa.sub( f1+f2  );
            if(bSubtractAngleNonBond){
                Vec3f fij=Vec3fZero;
                Quat4f REQij; combineREQ( REQs[ing],REQs[jng], REQij );
                Vec3f dij; dij.set_lincomb( -1./hi.e, hi.f, 1./hj.e, hj.f );  // method without reading from global buffer
                //Vec3f dij = apos[jng] - apos[ing];                          // method with    reading from global buffer
                E -= addAtomicForceLJQ( dij, fij, REQij );
                f1.sub(fij);
                f2.add(fij);
            }
            fbs[i].f.add( f1     );
            fbs[j].f.add( f2     );
            
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

void normalizePis(){ 
    for(int i=0; i<nnode; i++){ pipos[i].f.normalize(); pipos[i].e=0; } 
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
        fapos [ia].e=0; 
        if(bpi){
            fpipos[ia].add( fp );
            fpipos[ia].e=0;
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
    //if(idebug){printf("CPU BEFORE assemble() \n"); printDEBUG();} 
    asseble_forces();
    //Etot = Eb + Ea + Eps + EppT + EppI;
    return Etot;
}

double eval_check(){
    if(verbosity>0){
        printf(" ============ check MMFFf4 START\n " );
        printSizes();
    }
    eval();
    checkNans();
    if(verbosity>0)printf(" ============ check MMFFf4 DONE\n " );
    return Etot;
}

void flipPis( Vec3f ax ){
    for(int i=0; i<nnode; i++){
        float c = pipos[i].f.dot(ax);
        if( c<0 ){ pipos[i].mul(-1); } 
    }
}

void move_GD(float dt, float Flim=100.0 ){
    //for(int i=0; i<0; i++){}
    //Quat4f vs = (Quat4f*) vDOFs; 
    float F2lim=Flim*Flim;
    for(int i=0; i<nvecs; i++){
        //vs[i].f.mul( f, dt );
        Vec3f f = fapos[i].f;
        float fr2 = f.norm2();
        if(fr2>F2lim){ f.mul(Flim/sqrt(fr2)); };
        apos[i].f.add_mul( f, dt ); 
    }
}

void makeBackNeighs( bool bCapNeighs=true ){
    for(int i=0; i<natoms; i++){ bkneighs[i]=Quat4i{-1,-1,-1,-1}; };
    for(int ia=0; ia<nnode; ia++){
        for(int j=0; j<4; j++){        // 4 neighbors
            int ja = neighs[ia].array[j];
            if( ja<0 )continue;
            //NOTE: We deliberately ignore back-neighbors from caping atoms 
            bool ret = addFirstEmpty( bkneighs[ja].array, 4, ia*4+j, -1 );
            if(!ret){ printf("ERROR in MMFFf4_loc::makeBackNeighs(): Atom #%i has >4 back-Neighbors (while adding atom #%i) \n", ja, ia ); exit(0); }
        };
    }
    if(bCapNeighs){   // set neighbors for capping atoms
        for(int ia=nnode; ia<natoms; ia++){ neighs[ia]=Quat4i{-1,-1,-1,-1}; neighs[ia].x = bkneighs[ia].x/4;  }
    }
}

void makeNeighCells( const Vec3i nPBC_ ){ 
    //printf( "makeNeighCells() nPBC_(%i,%i,%i) lvec (%g,%g,%g) (%g,%g,%g) (%g,%g,%g)\n", nPBC.x,nPBC.y,nPBC.z, lvec.a.x,lvec.a.y,lvec.a.z,  lvec.b.x,lvec.b.y,lvec.b.z,   lvec.c.x,lvec.c.y,lvec.c.z );
    nPBC=nPBC_;
    for(int ia=0; ia<natoms; ia++){
        for(int j=0; j<4; j++){
            //printf("ngcell[%i,j=%i] \n", ia, j);
            int ja = neighs[ia].array[j];
            //printf("ngcell[%i,ja=%i] \n", ia, ja);
            if( ja<0 )continue;
            Vec3f d = apos[ja].f - apos[ia].f;
            int ipbc=0;
            int imin=-1;
            float r2min = 1e+30;
            for(int ic=-nPBC.z; ic<=nPBC.z; ic++){ for(int ib=-nPBC.y; ib<=nPBC.y; ib++){ for(int ia=-nPBC.x; ia<=nPBC.x; ia++){   
                Vec3f shift= (lvec.a*ia) + (lvec.b*ib) + (lvec.c*ic); 
                shift.add(d);
                float r2 = shift.norm2();
                if(r2<r2min){   // find nearest distance
                    r2min=r2;
                    imin=ipbc;
                }
                ipbc++; 
            }}}
            //printf("ngcell[%i,%i] imin=%i \n", ia, ja, imin);
            neighCell[ia].array[j] = imin;
            //printf("ngcell[%i,%i] imin=%i ---- \n", ia, ja, imin);
        }
    }
}

void makeNeighCells( int npbc, Vec3d* pbc_shifts ){ 
    for(int ia=0; ia<natoms; ia++){
        for(int j=0; j<4; j++){
            //printf("ngcell[%i,j=%i] \n", ia, j);
            int ja = neighs[ia].array[j];
            //printf("ngcell[%i,ja=%i] \n", ia, ja);
            if( ja<0 )continue;
            Vec3f d = apos[ja].f - apos[ia].f;
            int imin=-1;
            float r2min = 1e+30;
            for( int ipbc=0; ipbc<npbc; ipbc++ ){   
                Vec3f shift= (Vec3f)pbc_shifts[ipbc]; 
                shift.add(d);
                float r2 = shift.norm2();
                if(r2<r2min){   // find nearest distance
                    r2min=r2;
                    imin=ipbc;
                }
            }
            //printf("ngcell[%i,%i] imin=%i \n", ia, ja, imin);
            neighCell[ia].array[j] = imin;
            //printf("ngcell[%i,%i] imin=%i ---- \n", ia, ja, imin);
        }
    }
}

void printSizes(){ printf( "MMFFf4::printSizes(): nDOFs(%i) natoms(%i) nnode(%i) ncap(%i) nvecs(%i) \n", nDOFs,natoms,nnode,ncap,nvecs ); };
void printAtomParams(int ia){ printf("atom[%i] ngs{%3i,%3i,%3i,%3i} par(%5.3f,%5.3f,%5.3f)  bL(%5.3f,%5.3f,%5.3f,%5.3f) bK(%6.3f,%6.3f,%6.3f,%6.3f)  Ksp(%5.3f,%5.3f,%5.3f,%5.3f) Kpp(%5.3f,%5.3f,%5.3f,%5.3f) \n", ia, neighs[ia].x,neighs[ia].y,neighs[ia].z,neighs[ia].w,    apars[ia].x,apars[ia].y,apars[ia].z,    bLs[ia].x,bLs[ia].y,bLs[ia].z,bLs[ia].w,   bKs[ia].x,bKs[ia].y,bKs[ia].z,bKs[ia].w,     Ksp[ia].x,Ksp[ia].y,Ksp[ia].z,Ksp[ia].w,   Kpp[ia].x,Kpp[ia].y,Kpp[ia].z,Kpp[ia].w  ); };
void printAtomParams(){for(int ia=0; ia<nnode; ia++){ printAtomParams(ia); }; };

void printBKneighs(int ia){ printf("atom[%i] bkngs{%3i,%3i,%3i,%3i} \n", ia, bkneighs[ia].x,bkneighs[ia].y,bkneighs[ia].z,bkneighs[ia].w ); };
void printBKneighs(){for(int ia=0; ia<natoms; ia++){ printBKneighs(ia); }; };

void print_apos(){
    for(int ia=0;ia<natoms;ia++){ printf( "print_apos[%i](%g,%g,%g)\n", ia, apos[ia].x,apos[ia].y,apos[ia].z ); }
}

bool checkNans( bool bExit=true, bool bNg=true, bool bPi=true, bool bA=true ){
    bool ret = false;
    if(bA)  ret |= ckeckNaN_f(natoms,  4, (float*)  apos,   "apos"    );
    if(bA)  ret |= ckeckNaN_f(natoms,  4, (float*) fapos,  "fapos"    );
    if(bPi) ret |= ckeckNaN_f(nnode,   4, (float*) pipos,  "pipos"    );
    if(bPi) ret |= ckeckNaN_f(nnode,   4, (float*)fpipos, "fpipos"    );
    if(bNg) ret |= ckeckNaN_f(nnode*4, 4, (float*)fneigh,  "fneigh"   );
    if(bNg) ret |= ckeckNaN_f(nnode*4, 4, (float*)fneighpi,"fneighpi" );
    if(bExit&&ret){ printf("ERROR: NaNs detected in %s in %s => exit(0)\n", __FUNCTION__, __FILE__ ); exit(0); };
    return ret;
}

void printDEBUG(  bool bNg=true, bool bPi=true, bool bA=true ){
    printf( "MMFFf4::printDEBUG() \n" );
    if(bA)for(int i=0; i<natoms; i++){
        printf( "CPU[%i] ", i );
        //printf( "bkngs{%2i,%2i,%2i,%2i} ",         bkNeighs[i].x, bkNeighs[i].y, bkNeighs[i].z, bkNeighs[i].w );
        printf( "fapos{%6.3f,%6.3f,%6.3f,%6.3f} ", fapos[i].x, fapos[i].y, fapos[i].z, fapos[i].w );
        //printf(  "avel{%6.3f,%6.3f,%6.3f,%6.3f} ", avel[i].x, avel[i].y, avel[i].z, avel[i].w );
        //printf(  "apos{%6.3f,%6.3f,%6.3f,%6.3f} ", apos[i].x, apos[i].y, apos[i].z, apos[i].w );
        printf( "\n" );
    }
    if(bPi)for(int i=0; i<nnode; i++){
        int i1=i+natoms;
        printf( "CPU[%i] ", i1 );
        printf(  "fpipos{%6.3f,%6.3f,%6.3f,%6.3f} ", fapos[i1].x, fapos[i1].y, fapos[i1].z, fapos[i1].w );
        //printf(  "vpipos{%6.3f,%6.3f,%6.3f,%6.3f} ", avel[i1].x, avel[i1].y, avel[i1].z, avel[i1].w );
        //printf(   "pipos{%6.3f,%6.3f,%6.3f,%6.3f} ", apos[i1].x, apos[i1].y, apos[i1].z, apos[i1].w );
        printf( "\n" );
    }
    if(bNg)for(int i=0; i<nnode; i++){ for(int j=0; j<4; j++){
        int i1=i*4+j;
        //int i2=(i+natoms)*4+j;
        printf( "CPU[%i,%i] ", i, j );
        printf( "fneigh  {%6.3f,%6.3f,%6.3f,%6.3f} ", fneigh  [i1].x, fneigh  [i1].y, fneigh  [i1].z, fneigh  [i1].w );
        printf( "fneighpi{%6.3f,%6.3f,%6.3f,%6.3f} ", fneighpi[i1].x, fneighpi[i1].y, fneighpi[i1].z, fneighpi[i1].w );
        printf( "\n" );
    }}
}

};


#endif
