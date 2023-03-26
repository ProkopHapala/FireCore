
#ifndef MMFFsp3_loc_h
#define MMFFsp3_loc_h

#define ANG_HALF_COS  1


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

    bool doBonds  =true;
    bool doNeighs =true;
    bool doPiPiI  =true;
    bool doPiPiT  =true;
    bool doPiSigma=true;
    bool doAngles =true;
    bool doEpi    =true; 
    
    //                           c0     Kss    Ksp    c0_e
    Quat4d default_NeighParams{ -1.0,   1.0,   1.0,   -1.0 };

    // Dynamical Varaibles;
    Vec3d *   apos=0;   // [natom]
    Vec3d *  fapos=0;   // [natom]
    Vec3d *  pipos=0;   // [nnode]
    Vec3d * fpipos=0;   // [nnode]
    // Aux Dynamil
    Vec3d * fneigh  =0;  // [nnode*4]     temporary store of forces on atoms form neighbors (before assembling step)
    Vec3d * fneighpi=0;  // [nnode*4]     temporary store of forces on pi    form neighbors (before assembling step)

    // Params
    int   *  atypes  =0;
    Quat4i*  aneighs =0;   // [natoms]  index of neighboring atoms
    Quat4i*  bkneighs=0;   // [natoms]  inverse neighbors
    Quat4i*  aneighCell=0; // [natoms]  cell index for neighbors

    Quat4d*  apars=0;  // [nnode] per atom forcefield parametrs
    Quat4d*  bLs  =0;  // [nnode] bond lengths
    Quat4d*  bKs  =0;  // [nnode] bond stiffness
    Quat4d*  Ksp  =0;  // [nnode] stiffness of pi-alignment
    Quat4d*  Kpp  =0;  // [nnode] stiffness of pi-planarization
    Vec3d*  REQs =0;   // [nnode] parameters of non-covalent interactions

    bool    bAngleCosHalf         = true;
    bool    bSubtractAngleNonBond = false;
    double  Rdamp  = 1.0;
    Mat3d   invLvec, lvec;
    Vec3i   nPBC;

// =========================== Functions

void realloc( int nnode_, int ncap_ ){
    nnode=nnode_; ncap=ncap_;
    natoms= nnode + ncap; 
    nvecs = natoms+nnode;  // each atom as also pi-orientiation (like up-vector)
    nDOFs = nvecs*3;
    //printf( "MMFFsp3::realloc() natom(%i,nnode=%i,ncap=%i), npi=%i, nbond=%i \n", natoms, nnode, ncap, npi, nbonds );
    int ipi0=natoms;

    /*
    // ----- Dynamical
    _realloc(  DOFs    , nDOFs );
    _realloc( fDOFs    , nDOFs );
    apos   = (Vec3d*) DOFs ;
    fapos  = (Vec3d*)fDOFs;
    pipos  = apos  + ipi0;
    fpipos = fapos + ipi0;
    // ---- Aux
    _realloc( fneigh  , nnode*4 );
    _realloc( fneighpi, nnode*4 );
    // ----- Params [natom]
    _realloc( atypes    , natoms );
    _realloc( aneighs   , natoms );
    _realloc( aneighCell, natoms );
    _realloc( bkneighs  , natoms );
    _realloc( apars     , nnode );
    _realloc( bLs       , nnode );
    _realloc( bKs       , nnode );
    _realloc( Ksp       , nnode );
    _realloc( Kpp       , nnode );
    */
    
    
    /*
    _realloc0(  DOFs    , nDOFs , 0. );
    _realloc0( fDOFs    , nDOFs , 0. );
    apos   = (Vec3d*) DOFs ;
    fapos  = (Vec3d*)fDOFs;
    pipos  = apos  + ipi0;
    fpipos = fapos + ipi0;
    // ---- Aux
    _realloc0( fneigh  , nnode*4, Vec3dZero );
    _realloc0( fneighpi, nnode*4, Vec3dZero );
    // ----- Params [natom]
    _realloc0( atypes    , natoms, -1 );
    _realloc0( aneighs   , natoms, Quat4iMinusOnes );
    _realloc0( aneighCell, natoms, Quat4iMinusOnes );
    _realloc0( bkneighs  , natoms, Quat4iMinusOnes);
    _realloc0( apars     , nnode, Quat4dZero );
    _realloc0( bLs       , nnode, Quat4dZero );
    _realloc0( bKs       , nnode, Quat4dZero );
    _realloc0( Ksp       , nnode, Quat4dZero );
    _realloc0( Kpp       , nnode, Quat4dZero );
    */
    
    _realloc0(  DOFs    , nDOFs , (double)NAN );
    _realloc0( fDOFs    , nDOFs , (double)NAN );
    apos   = (Vec3d*) DOFs ;
    fapos  = (Vec3d*)fDOFs;
    pipos  = apos  + ipi0;
    fpipos = fapos + ipi0;
    // ---- Aux
    _realloc0( fneigh  , nnode*4, Vec3dNAN );
    _realloc0( fneighpi, nnode*4, Vec3dNAN );
    // ----- Params [natom]
    _realloc0( atypes    , natoms, -1 );
    _realloc0( aneighs   , natoms, Quat4iMinusOnes );
    _realloc0( aneighCell, natoms, Quat4iMinusOnes );
    _realloc0( bkneighs  , natoms, Quat4iMinusOnes);
    _realloc0( apars     , nnode, Quat4dNAN );
    _realloc0( bLs       , nnode, Quat4dNAN );
    _realloc0( bKs       , nnode, Quat4dNAN );
    _realloc0( Ksp       , nnode, Quat4dNAN );
    _realloc0( Kpp       , nnode, Quat4dNAN );

}

void dealloc(){
    _dealloc(DOFs );
    _dealloc(fDOFs);
    apos   = 0;
    fapos  = 0;
    pipos  = 0;
    fpipos = 0;
    _dealloc(atypes);
    _dealloc(aneighs);
    _dealloc(aneighCell);
    _dealloc(bkneighs);
    _dealloc(apars);
    _dealloc(bLs);
    _dealloc(bKs);
    _dealloc(Ksp);
    _dealloc(Kpp);
}

void setLvec(const Mat3d& lvec_){ lvec=lvec_; lvec.invert_T_to( invLvec ); }

double optimalTimeStep(double m=1.0){
    double Kmax = 1.0;
    for(int i=0; i<nnode; i++){ 
        Kmax=fmax(Kmax, bKs[i].x ); 
        Kmax=fmax(Kmax, bKs[i].y ); 
        Kmax=fmax(Kmax, bKs[i].z ); 
        Kmax=fmax(Kmax, bKs[i].w ); 
    }
    return M_PI*2.0*sqrt(m/Kmax)/10.0;  // dt=T/10;   T = 2*pi/omega = 2*pi*sqrt(m/k)
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
    Vec3d* fbs  = fneigh   +ia*4;
    Vec3d* fps  = fneighpi +ia*4;

    // // --- settings
    // double  ssC0 = apars[ia].x;
    // double  ssK  = apars[ia].y;
    // double  piC0 = apars[ia].z;
    // //bool    bPi  = ings[3]<0;   we distinguish this by Ksp, otherwise it would be difficult for electron pairs e.g. (-O-C=)

    const Quat4d& apar  = apars[ia];
    const double  ssK  = apar.z;
    const double  piC0 = apar.w;
    const Vec2d cs0_ss = Vec2d{apar.x,apar.y};
    const double  ssC0 = cs0_ss.x*cs0_ss.x - cs0_ss.y*cs0_ss.y;   // cos(2x) = cos(x)^2 - sin(x)^2, because we store cos(ang0/2) to use in  evalAngleCosHalf

    //--- Aux Variables 
    Quat4d  hs[4];
    Vec3d   f1,f2;
    
    const int ia_DBG = 0;
    //if(ia==ia_DBG)printf( "ffl[%i] neighs(%i,%i,%i,%i) \n", ia, ings[0],ings[1],ings[2],ings[3] );

    //for(int i=0; i<4; i++){ fneigh[ia*4+i]=Vec3dZero; fneighpi[ia*4+i]=Vec3dZero; }
    for(int i=0; i<4; i++){ fbs[i]=Vec3dZero; fps[i]=Vec3dZero; } // we initialize it here because of the break

    // --------- Bonds Step
    for(int i=0; i<4; i++){
        int ing = ings[i];
        //printf( "bond[%i|%i=%i]\n", ia,i,ing );
        //fbs[i]=Vec3dOne; fps[i]=Vec3dOne;
        //fbs[i]=Vec3dZero; fps[i]=Vec3dZero; // NOTE: wee need to initialize it before, because of the break
        if(ing<0) break;
        Quat4d h; 
        h.f.set_sub( apos[ing], pa );
        //if(idebug)printf( "bond[%i|%i=%i] l=%g pj[%i](%g,%g,%g) pi[%i](%g,%g,%g)\n", ia,i,ing, h.f.norm(), ing,apos[ing].x,apos[ing].y,apos[ing].z, ia,pa.x,pa.y,pa.z  );
        double l = h.f.normalize();
        h.e    = 1/l;
        hs [i] = h;
        // bond length force
        //continue; 

        //if(ia==ia_DBG) printf( "ffl:h[%i|%i=%i] l %g h(%g,%g,%g) pj(%g,%g,%g) pa(%g,%g,%g) \n", ia,i,ing, l, h.x,h.y,h.z, apos[ing].x,apos[ing].y,apos[ing].z, pa.x,pa.y,pa.z );

        if(ia<ing){   // we should avoid double counting because otherwise node atoms would be computed 2x, but capping only once
            E+= evalBond( h.f, l-bL[i], bK[i], f1 ); fbs[i].sub(f1);  fa.add(f1);    
            //if(ia==ia_DBG)printf( "ffl:bond[%i|%i=%i] kb=%g l0=%g l=%g h(%g,%g,%g) f(%g,%g,%g) \n", ia,i,ing, bK[i],bL[i], l, h.x,h.y,h.z,  f1.x,f1.y,f1.z  );
            
            double kpp = Kppi[i];
            if( (ing<nnode) && (kpp>1e-6) ){   // Only node atoms have pi-pi alignemnt interaction
                E += evalPiAling( hpi, pipos[ing], 1., 1.,   kpp,       f1, f2 );   fpi.add(f1);  fps[i].add(f2);    //   pi-alignment     (konjugation)
                //if(ia==ia_DBG)printf( "ffl:pp[%i|%i] kpp=%g c=%g f1(%g,%g,%g) f2(%g,%g,%g)\n", ia,ing, kpp, hpi.dot(pipos[ing]), f1.x,f1.y,f1.z,  f2.x,f2.y,f2.z  );
            }
            // ToDo: triple bonds ?
            
        } 
        
        // pi-sigma 
        //if(bPi){    
        double ksp = Kspi[i];
        if(ksp>1e-6){  
            E += evalAngleCos( hpi, h.f      , 1., h.e, ksp, piC0, f1, f2 );   fpi.add(f1); fa.sub(f2);  fbs[i].add(f2);    //   pi-planarization (orthogonality)
            //if(ia==ia_DBG)printf( "ffl:sp[%i|%i] ksp=%g piC0=%g c=%g hp(%g,%g,%g) h(%g,%g,%g)\n", ia,ing, ksp,piC0, hpi.dot(h.f), hpi.x,hpi.y,hpi.z,  h.x,h.y,h.z  );
            //if(ia==ia_DBG)printf( "ffl:sp[%i|%i] ksp=%g piC0=%g c=%g f1(%g,%g,%g) f2(%g,%g,%g)\n", ia,ing, ksp,piC0, hpi.dot(h.f), f1.x,f1.y,f1.z,  f2.x,f2.y,f2.z  );
        }
        //}
        
    }

    
    // --------- Angle Step
    const double R2damp=Rdamp*Rdamp;
    for(int i=0; i<4; i++){
        int ing = ings[i];
        if(ing<0) break;
        const Quat4d& hi = hs[i];
        for(int j=i+1; j<4; j++){
            int jng  = ings[j];
            if(jng<0) break;
            const Quat4d& hj = hs[j];    
            // #if ANG_HALF_COS
            //     E += evalAngleCosHalf( hi.f, hj.f,  hi.e, hj.e,  cs0_ss,  ssK, f1, f2 );
            // #else             
            //     E += evalAngleCos( hi.f, hj.f, hi.e, hj.e, ssK, ssC0, f1, f2 );     // angles between sigma bonds
            // #endif

            if( bAngleCosHalf ){
                E += evalAngleCosHalf( hi.f, hj.f,  hi.e, hj.e,  cs0_ss,  ssK, f1, f2 );
            }else{             
                E += evalAngleCos( hi.f, hj.f, hi.e, hj.e, ssK, ssC0, f1, f2 );     // angles between sigma bonds
            }
            //if(ia==ia_DBG)printf( "ffl:ang[%i|%i,%i] kss=%g cs0(%g,%g) c=%g l(%g,%g) f1(%g,%g,%g) f2(%g,%g,%g)\n", ia,ing,jng, ssK, cs0_ss.x,cs0_ss.y, hi.f.dot(hj.f),hi.w,hj.w, f1.x,f1.y,f1.z,  f2.x,f2.y,f2.z  );

            fa    .sub( f1+f2  );
            if(bSubtractAngleNonBond){
                Vec3d fij=Vec3dZero;
                //printf( "non-bond[%i|%i=%i,%i,=%i] REQs=%li \n", ia, i,ing,j,jng, REQs  );
                //printf( "non-bond[%i|%i=%i,%i,=%i] REQi(%g,%g,%g) REQj(%g,%g,%g) \n", ia, i,ing,j,jng, REQs[ing].x,REQs[ing].y,REQs[ing].z,   REQs[jng].x,REQs[jng].y,REQs[jng].z  );
                Vec3d REQij; combineREQ( REQs[ing],REQs[jng], REQij );
                Vec3d dp; dp.set_lincomb( -1./hi.e, hi.f, 1./hj.e, hj.f );  // method without reading from global buffer
                //Vec3d dp = apos[jng] - apos[ing];                          // method with    reading from global buffer
                //E -= addAtomicForceLJQ( dp, fij, REQij );
                E -= getLJQ( dp, REQij, R2damp, fij );
                //if(ia==ia_DBG)printf( "ffl:LJQ[%i|%i,%i] r=%g REQ(%g,%g,%g) fij(%g,%g,%g)\n", ia,ing,jng, dp.norm(), REQij.x,REQij.y,REQij.z, fij.x,fij.y,fij.z );
                f1.sub(fij);
                f2.add(fij);
            }
            fbs[i].add( f1     );
            fbs[j].add( f2     );
            //if(ia==ia_DBG)printf( "ffl:ANG[%i|%i,%i] fa(%g,%g,%g) fbs[%i](%g,%g,%g) fbs[%i](%g,%g,%g)\n", ia,ing,jng, fa.x,fa.y,fa.z, i,fbs[i].x,fbs[i].y,fbs[i].z,   j,fbs[j].x,fbs[j].y,fbs[j].z  );
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
    // NOTE: We do not need clean fneigh,fneighpi because they are set in eval_atoms 
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
            fa.add(fneigh  [j]);
            if(bpi){
                //printf( "assemble[%i,%i|%i] pi(%g,%g,%g) fp(%g,%g,%g) fpng(%g,%g,%g) \n", ia,i,j, pipos[ia].x,pipos[ia].y,pipos[ia].z, fpipos[ia].x,fpipos[ia].y,fpipos[ia].z, fneighpi[j].x,fneighpi[j].y,fneighpi[j].z );
                fp.add(fneighpi[j]);
            }
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
    //if(idebug){printf("CPU BEFORE assemble() \n"); printDEBUG();} 
    asseble_forces();
    //Etot = Eb + Ea + Eps + EppT + EppI;
    return Etot;
}

double eval_check(){
    printf(" ============ check MMFFsp3_loc START\n " );
    printSizes();
    eval();
    checkNans();
    printf(" ============ check MMFFsp3_loc DONE\n " );
    return Etot;
}

void flipPis( Vec3d ax ){
    for(int i=0; i<nnode; i++){
        double c = pipos[i].dot(ax);
        if( c<0 ){ pipos[i].mul(-1); } 
    }
}


void move_GD(float dt, double Flim=100.0 ){
    double F2lim=Flim*Flim;
    for(int i=0; i<nvecs; i++){
        Vec3d  f   = fapos[i];
        double fr2 = f.norm2();
        if(fr2>F2lim){ f.mul(Flim/sqrt(fr2)); };
        apos[i].add_mul( f, dt ); 
    }
}

void makeBackNeighs( bool bCapNeighs=true ){
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
    if(bCapNeighs){   // set neighbors for capping atoms
        for(int ia=nnode; ia<natoms; ia++){ aneighs[ia]=Quat4i{-1,-1,-1,-1};  aneighs[ia].x = bkneighs[ia].x/4;  }
    }
}

void makeNeighCells( const Vec3i nPBC_ ){ 
    nPBC=nPBC_;
    for(int ia=0; ia<natoms; ia++){
        aneighCell[ia]=Quat4i{0,0,0,0};
        for(int j=0; j<4; j++){
            //printf("ngcell[%i,j=%i] \n", ia, j);
            int ja = aneighs[ia].array[j];
            //printf("ngcell[%i,ja=%i] \n", ia, ja);
            if( ja<0 )continue;
            Vec3d d = apos[ja] - apos[ia];
            int ipbc=0;
            int imin=-1;
            double r2min = 1e+300;
            for(int ia=-nPBC.x; ia<=nPBC.x; ia++){ for(int ib=-nPBC.y; ib<=nPBC.y; ib++){ for(int ic=-nPBC.z; ic<=nPBC.z; ic++){ 
                Vec3d shift= (lvec.a*ia) + (lvec.b*ib) + (lvec.c*ic); 
                shift.add(d);
                double r2 = shift.norm();
                if(r2<r2min){   // find nearest distance
                    r2min=r2;
                    imin=ipbc;
                }
                ipbc++; 
            }}}
            //printf("ngcell[%i,%i] imin=%i \n", ia, ja, imin);
            aneighCell[ia].array[j] = imin;
            //printf("ngcell[%i,%i] imin=%i ---- \n", ia, ja, imin);
        }
    }
    //printf("MMFFsp3_loc::makeNeighCells() DONE \n");
}

void printSizes(){ printf( "MMFFf4::printSizes(): nDOFs(%i) natoms(%i) nnode(%i) ncap(%i) nvecs(%i) \n", nDOFs,natoms,nnode,ncap,nvecs ); };
void printAtomParams(int ia){ printf("atom[%i] ngs{%3i,%3i,%3i,%3i} par(%5.3f,%5.3f,%5.3f)  bL(%5.3f,%5.3f,%5.3f,%5.3f) bK(%6.3f,%6.3f,%6.3f,%6.3f)  Ksp(%5.3f,%5.3f,%5.3f,%5.3f) Kpp(%5.3f,%5.3f,%5.3f,%5.3f) \n", ia, aneighs[ia].x,aneighs[ia].y,aneighs[ia].z,aneighs[ia].w,    apars[ia].x,apars[ia].y,apars[ia].z,    bLs[ia].x,bLs[ia].y,bLs[ia].z,bLs[ia].w,   bKs[ia].x,bKs[ia].y,bKs[ia].z,bKs[ia].w,     Ksp[ia].x,Ksp[ia].y,Ksp[ia].z,Ksp[ia].w,   Kpp[ia].x,Kpp[ia].y,Kpp[ia].z,Kpp[ia].w  ); };
void printAtomParams(){for(int ia=0; ia<nnode; ia++){ printAtomParams(ia); }; };

void printBKneighs(int ia){ printf("atom[%i] bkngs{%3i,%3i,%3i,%3i} \n", ia, bkneighs[ia].x,bkneighs[ia].y,bkneighs[ia].z,bkneighs[ia].w ); };
void printBKneighs(){for(int ia=0; ia<natoms; ia++){ printBKneighs(ia); }; };


void print_apos(){
    for(int ia=0;ia<natoms;ia++){ printf( "print_apos[%i](%g,%g,%g)\n", ia, apos[ia].x,apos[ia].y,apos[ia].z ); }
}

void printDEBUG(  bool bNg=true, bool bPi=true, bool bA=true ){
    printf( "MMFFsp3_loc::printDEBUG() \n" );
    if(bA)for(int i=0; i<natoms; i++){
        printf( "CPU[%i] ", i );
        //printf( "bkngs{%2i,%2i,%2i,%2i,%2i} ",         bkNeighs[i].x, bkNeighs[i].y, bkNeighs[i].z, bkNeighs[i].w );
        printf( "fapos{%6.3f,%6.3f,%6.3f} ", fapos[i].x, fapos[i].y, fapos[i].z );
        //printf(  "avel{%6.3f,%6.3f,%6.3f} ", avel[i].x, avel[i].y, avel[i].z);
        printf(  "apos{%6.3f,%6.3f,%6.3f} ", apos[i].x, apos[i].y, apos[i].z );
        printf( "\n" );
    }
    if(bPi)for(int i=0; i<nnode; i++){
        int i1=i+natoms;
        printf( "CPU[%i] ", i1 );
        printf(  "fpipos{%6.3f,%6.3f,%6.3f} ", fapos[i1].x, fapos[i1].y, fapos[i1].z );
        //printf(  "vpipos{%6.3f,%6.3f,%6.3f} ", avel[i1].x, avel[i1].y, avel[i1].z );
        printf(   "pipos{%6.3f,%6.3f,%6.3f} ", apos[i1].x, apos[i1].y, apos[i1].z );
        printf( "\n" );
    }
    if(bNg)for(int i=0; i<nnode; i++){ for(int j=0; j<4; j++){
        int i1=i*4+j;
        //int i2=(i+natoms)*4+j;
        printf( "CPU[%i,%i] ", i, j );
        printf( "fneigh  {%6.3f,%6.3f,%6.3f} ", fneigh  [i1].x, fneigh  [i1].y, fneigh  [i1].z );
        printf( "fneighpi{%6.3f,%6.3f,%6.3f} ", fneighpi[i1].x, fneighpi[i1].y, fneighpi[i1].z );
        printf( "\n" );
    }}
}

bool checkNans( bool bExit=true, bool bNg=true, bool bPi=true, bool bA=true ){
    bool ret = false;
    if(bA)  ret |= ckeckNaN_d(natoms,  3, (double*) apos,   "apos"  );
    if(bA)  ret |= ckeckNaN_d(natoms,  3, (double*)fapos,  "fapos"  );
    if(bPi) ret |= ckeckNaN_d(nnode,   3, (double*) pipos,  "pipos" );
    if(bPi) ret |= ckeckNaN_d(nnode,   3, (double*)fpipos, "fpipos" );
    if(bNg) ret |= ckeckNaN_d(nnode*4, 3, (double*)fneigh,  "fneigh"   );
    if(bNg) ret |= ckeckNaN_d(nnode*4, 3, (double*)fneighpi,"fneighpi" );
    if(bExit&&ret){ printf("ERROR: NaNs detected in %s in %s => exit(0)\n", __FUNCTION__, __FILE__ ); exit(0); };
    return ret;
}

void rotateNodes(int n, int* sel, Vec3d p0, Vec3d ax, double phi ){
    ax.normalize();
    double ca=cos(phi);
    double sa=sin(phi);
    for(int i=0;i<n; i++){
        int ia = sel[i];
        if(ia>=nnode)continue;
        apos [ia].rotate_csa( ca, sa, ax, p0 );
        pipos[ia].rotate_csa( ca, sa, ax     );
        int* ngs=aneighs[ia].array; 
        for(int j=0;j<4;j++){
            int ja = ngs[j];
            if(ja>=0){ if(ja>nnode) apos[ ja  ].rotate_csa( ca, sa, ax, p0 ); }
        }
    }
}

void chargeToEpairs( Vec3d* REQs, int* atypes, double cQ=-0.2, int etyp=-1 ){
    for( int ia=0; ia<nnode; ia++ ){
        int* ngs=aneighs[ia].array; 
        for( int j=0; j<4; j++ ){
            int ja = ngs[j]; 
            if(ja<0) continue;
            if( atypes[ja]==etyp ){ REQs[ja].z+=cQ; REQs[ia].z-=cQ; };
        }
    }
}

};


#endif
