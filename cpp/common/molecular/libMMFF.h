
#include "Forces.h"
#include "InterpolateTricubic.h"
#include "Bspline.h"

extern "C"{

// ================ COMMONS

void setVerbosity( int verbosity_, int idebug_ ){
    verbosity = verbosity_;
    idebug    = idebug_;
}

// ================ INITIALIZATION

//int loadmol(char* fname_mol ){ return W.loadmol(fname_mol ); }
void insertSMILES( char* s ){  W.insertSMILES(s); };

void initParams          ( const char* sElementTypes, const char* sAtomTypes, const char* sBondTypes, const char* sAngleTypes, const char* sDihedralTypes ){ W.initParams(sElementTypes,sAtomTypes,sBondTypes,sAngleTypes,sDihedralTypes); }

//int  buildMolecule_xyz   ( const char* xyz_name, bool bEpairs, double fAutoCharges ){  W.builder.bDummyEpair=bEpairs; W.bEpairs=bEpairs; W.fAutoCharges=fAutoCharges;  return W.buildMolecule_xyz( xyz_name );  }

int  buildMolecule_xyz( const char* xyz_name, bool bEpairs, double fAutoCharges, bool bAutoTypes, bool bRelaxPi ){  
    W.builder.bDummyEpair=bEpairs; W.builder.bAutoTypes=bAutoTypes; W.bRelaxPi=bRelaxPi; W.fAutoCharges=fAutoCharges;  
    //printf( "buildMolecule_xyz W.builder.bDummyEpair=%i bEpairs=%i \n", W.builder.bDummyEpair, bEpairs );
    return W.buildMolecule_xyz( xyz_name ); 
}

int  buildMolecule_SMILES( const char* xyz_name, bool bEpairs, double fAutoCharges ){  W.builder.bDummyEpair=bEpairs; W.bEpairs=bEpairs; W.fAutoCharges=fAutoCharges;  return W.buildMolecule_xyz( xyz_name );  }
void makeFFs          ( ){ W.makeFFs(); }
void clear            ( ){ W.clear();   }

void addDistConstrain( int i0,int i1, double lmin,double lmax,double kmin,double kmax,double flim, double* shift, bool bOldIndex  ){
    //W.constrs.bonds.push_back( DistConstr{ {i0,i1}, {lmax,lmin}, {kmax,kmin}, flim, *(Vec3d*)shift } );
    W.addDistConstrain( i0,i1, lmin,lmax,kmin,kmax,flim, *(Vec3d*)shift, bOldIndex );
    W.bConstrains=true;
}
void addAngConstrain( int i0,int i1,int i2, double ang0, double k ){
    W.constrs.angles.push_back( AngleConstr{ {i0,i1,i2}, {cos(ang0/2.),sin(ang0/2.)}, k, } );
}

// ================ UTILS

bool checkInvariants( double maxVcog, double maxFcog, double maxTg ){ return W.checkInvariants( maxVcog, maxFcog, maxTg ); }
void setTrjName     ( const char* trj_fname_, int savePerNsteps_, int* nPBC ){ W.trj_fname=trj_fname_; W.savePerNsteps=savePerNsteps_; if(verbosity>0)printf( "setTrjName(%s)\n", W.trj_fname ); W.nPBC_save=*(Vec3i*)nPBC; }
int toXYZ           (const char* comment="#comment"){ return W.toXYZ(comment); }
int saveXYZ( const char* fname, const char* comment, int imod){ 
    int ret=-1;
    switch (imod){
        case 0: { ret=W.builder.save2xyz(fname,comment ); } break;
        case 1: { ret=W        .saveXYZ (fname,comment ); } break;
    }
    return ret;  
}


AtomType* getAtomTypes()    { return &W.params.atypes[0];    }
int       getAtomTypeCount(){ return W.params.atypes.size(); }

// ================ RUN / EVAL

void setupCollisionDamping( int nstep, double medium, double bond, double ang, double nonB, double dRcut1, double dRcut2 ){
    W.ffl.colDamp.set( nstep, medium, bond, ang, nonB, dRcut1, dRcut2 );
    // bool    bCollisionDamping        = false; // if true we use collision damping
    // bool    bCollisionDampingNonBond = false;  // if true we use collision damping for non-bonded interactions
    // double  damping_medium           = 1.0;   // cdamp       = 1 -(damping_medium     /ndampstep     )
    // double  collisionDamping         = 0.1;   // col_damp    =     collisionDamping   /(dt*ndampstep )
    // double  collisionDamping_NB      = 0.1;   // col_damp_NB =     collisionDamping_NB/(dt*ndampstep )
    // int     ndampstep                = 10;    // how many steps it takes to decay velocity to to 1/e of the initial value
    // double  col_damp_dRcut           = 0.5;   // non-covalent collision damping interaction goes between 1.0 to 0.0 on interval  [ Rvdw , Rvdw+col_damp_dRcut ]
    // double col_damp      = 0.0;  //  collisionDamping   /(dt*ndampstep );
    // double col_damp_NB   = 0.0;  //  collisionDamping_NB/(dt*ndampstep );
    // W.ffl.bCollisionDamping        = collisionDamping   >0;
    // W.ffl.bCollisionDampingNonBond = collisionDamping_NB>0;
    // W.ffl.damping_medium           = damping_medium;
    // W.ffl.collisionDamping         = fmax( collisionDamping   ,0 );
    // W.ffl.collisionDamping_NB      = fmax( collisionDamping_NB,0 );
    // W.ffl.ndampstep                = ndampstep;
    // W.ffl.col_damp_dRcut1          = col_damp_dRcut1;
    // W.ffl.col_damp_dRcut2          = col_damp_dRcut2;
}


void setup_accel(int nstep_acc_min_, double cos_vf_acc_ ){
    W.ffl.colDamp.setup_accel( nstep_acc_min_, cos_vf_acc_ );
}



double eval (){ return W.eval(); };
//int    run  ( int nstepMax, double dt=-1, double Fconv=1e-6, int ialg=2, double* outE=0, double* outF=0 ){ return W.run(nstepMax,dt,Fconv,ialg,outE,outF);  }
//bool   relax( int niter,    double Ftol, bool bWriteTrj ){ return W.relax( niter, Ftol, bWriteTrj );}

// ================ PRINT STATE

void printTypes     ( ){ W.params.printAtomTypes(true); }
void printAtomConfs ( bool bOmmitCap, bool bNeighs ){ W.builder.printAtomConfs(bOmmitCap,bNeighs); }
void printAtomTypes ( ){ W.builder.printAtomTypes ( ); }
void printBonds     ( ){ W.builder.printBonds     ( ); }
void printBondParams( ){ W.builder.printBondParams( ); }
void printAtomParams( ){ W.ffl.printAtomParams( );     }
void printSwitches  ( ){ W.printSwitches( );           }

const char* getTypeName( int ia, bool fromFF){
    int it;
    if(fromFF){ if(ia>=W.nbmol.natoms)           return "?"; it = W.nbmol.atypes[ia];       }
    else      { if(ia>=W.builder.atoms.size())   return "?"; it = W.builder.atoms[ia].type; }
    if( (it<0) || (it>=W.params.atypes.size()) ) return "?";
    return W.params.atypes[it].name;
}

// ========= Manipulation with the molecule

void shift_atoms_ax ( int n, int* selection, double* d                               ){ W.shift_atoms ( n, selection,*(Vec3d*)d);                     };
void rotate_atoms_ax( int n, int* selection, double* p0, double* ax, double phi      ){ W.rotate_atoms( n, selection,*(Vec3d*)p0, *(Vec3d*)ax, phi ); };
void shift_atoms    ( int n, int* selection, int ia0, int ia1, double l              ){ W.shift_atoms ( n, selection, ia0, ia1, l );                  };
void rotate_atoms   ( int n, int* selection, int ia0, int iax0, int iax1, double phi ){ W.rotate_atoms( n, selection, ia0, iax0, iax1, phi );         };

double measureBond        (int ia, int ib         ){ return W.ffl.measureBondLegth(ia,ib); }
double measureAngle       (int ic, int ia, int ib ){ return W.ffl.measureAngle (ic,ia,ib); }
double measureAnglePiPi   (int ia, int ib         ){ return W.ffl.measureAnglePiPi(ia, ib, true ); }
double measureAngleSigmaPi(int ipi, int ia, int ib){ return W.ffl.measureAngleSigmaPi(ipi, ia, ib ); }

// ========= Force-Field Scanning

void scanTranslation_ax( int n, int* selection, double* vec, int nstep, double* Es, const char* trjName, bool bAddjustCaps ){
    if(selection==0){ selection=W.manipulation_sel; n=W.manipulation_nsel; }
    W.scanTranslation_ax( n, selection, *(Vec3d*)vec, nstep, Es, trjName, bAddjustCaps );
}
void scanTranslation( int n, int* selection, int ia0, int ia1, double l, int nstep, double* Es, const char* trjName, bool bAddjustCaps ){ 
    W.scanTranslation( n, selection, ia0, ia1, l, nstep, Es, trjName, bAddjustCaps ); 
}
void scanRotation_ax( int n, int* selection, double* p0, double* ax, double phi, int nstep, double* Es, const char* trjName ){
    if(p0==0) p0=(double*)&W.manipulation_p0;
    if(ax==0) ax=(double*)&W.manipulation_ax;
    if(selection==0){selection=W.manipulation_sel; n=W.manipulation_nsel; }
    W.scanRotation_ax( n, selection, *(Vec3d*)p0, *(Vec3d*)ax, phi, nstep, Es, trjName );
}
void scanRotation( int n, int* selection,int ia0, int iax0, int iax1, double phi, int nstep, double* Es, const char* trjName ){ 
    W.scanRotation( n, selection,ia0, iax0, iax1, phi, nstep, Es, trjName );
}

void scanAngleToAxis_ax( int n, int* selection, double r, double R, double* p0, double* ax, int nstep, double* angs, double* Es, const char* trjName ){
    W.scanAngleToAxis_ax( n, selection, r, R, *(Vec3d*)p0, *(Vec3d*)ax, nstep, angs, Es, trjName );
}

// ========= Force-Field Component Sampling  

int fit_Bspline( const int n, double* Gs, double* Es, double* Ws, double Ftol, int nmaxiter, double dt ){
    return Bspline::fit1D( n, Gs, Es, Ws, Ftol, nmaxiter, dt );
}
int fitEF_Bspline( double dg, const int n, double* Gs, double* fes, double* Ws, double Ftol, int nmaxiter, double dt ){
    return Bspline::fit1D_EF( dg, n, Gs,  (Vec2d*)fes, (Vec2d*)Ws, Ftol, nmaxiter, dt );
}

int fit2D_Bspline( const int* ns, double* Gs, double* Es, double* Ws, double Ftol, int nmaxiter, double dt ){
    return Bspline::fit2D( *(Vec2i*)ns, Gs,  Es, Ws, Ftol, nmaxiter, dt );
}

int fit3D_Bspline( const int* ns, double* Gs, double* Es, double* Ws, double Ftol, int nmaxiter, double dt ){
    return Bspline::fit3D( *(Vec3i*)ns, Gs,  Es, Ws, Ftol, nmaxiter, dt );
}

void sample_Bspline( double g0, double dg, int ng, double* Gs, int n, double* xs, double* fes ){
    Bspline::sample1D( g0,dg,ng,Gs, n, xs, (Vec2d*)fes );
}

void sample_Bspline3D( double* g0, double* dg, int* ng, double* G, int n, double* ps, double* fes ){
    long t0 = getCPUticks();
    Bspline::sample3D( *(Vec3d*)g0, *(Vec3d*) dg, *(Vec3i*)ng, G, n, (Vec3d*)ps, (Quat4d*)fes );    // sample3D(n=10000) time=2490.29[kTick] 249.029[tick/point]
    double t = (getCPUticks()-t0); printf( "sample_Bspline3D(n=%i) time=%g[kTick] %g[tick/point]\n", n, t*(1.e-3), t/n );
}



void sample_SplineHermite( double g0, double dg, int ng, double* Eg, int n, double* xs, double* fes ){
    Spline_Hermite::sample1D( g0,dg,ng,Eg, n, xs, (Vec2d*)fes );
}

void sample_SplineHermite2D( double* g0, double* dg, int* ng, double* Eg, int n, double* ps, double* fes ){
    long t0 = getCPUticks();
    Spline_Hermite::sample2D( *(Vec2d*)g0, *(Vec2d*)dg, *(Vec2i*)ng, Eg, n, (Vec2d*)ps, (Vec3d*)fes );     //    54 [tick/point] with -Ofast
    //Spline_Hermite::sample2D_avx( *(Vec2d*)g0, *(Vec2d*)dg, *(Vec2i*)ng, Eg, n, (Vec2d*)ps, (Vec3d*)fes );   //   43 [tick/point] with -Ofast  
    double t = (getCPUticks()-t0); printf( "sample2D(n=%i) time=%g[kTick] %g[tick/point]\n", n, t*(1.e-3), t/n );
}

void sample_SplineHermite3D( double* g0, double* dg, int* ng, double* Eg, int n, double* ps, double* fes ){
    long t0 = getCPUticks();
    Spline_Hermite::sample3D    ( *(Vec3d*)g0, *(Vec3d*) dg, *(Vec3i*)ng, Eg, n, (Vec3d*)ps, (Quat4d*)fes );    // sample3D(n=10000) time=2490.29[kTick] 249.029[tick/point]
    //Spline_Hermite::sample3D_avx( *(Vec3d*)g0, *(Vec3d*) dg, *(Vec3i*)ng, Eg, n, (Vec3d*)ps, (Quat4d*)fes );  // sample3D(n=10000) time=2578.11[kTick] 257.811[tick/point]
    double t = (getCPUticks()-t0); printf( "sample3D(n=%i) time=%g[kTick] %g[tick/point]\n", n, t*(1.e-3), t/n );
}

void sample_SplineHermite3D_f( float* g0, float* dg, int* ng, float* Eg, int n, float* ps, float* fes ){
    long t0 = getCPUticks();
    Spline_Hermite::sample3D( *(Vec3f*)g0, *(Vec3f*) dg, *(Vec3i*)ng, Eg, n, (Vec3f*)ps, (Quat4f*)fes );    // sample3D(n=10000) time=2490.29[kTick] 249.029[tick/point]
    double t = (getCPUticks()-t0); printf( "sample3D_f(n=%i) time=%g[kTick] %g[tick/point]\n", n, t*(1.e-3), t/n );
}

void sample_SplineHermite3D_deriv( double* g0, double* dg, int* ng, double* Eg, double* dEg, int n, double* ps, double* fes ){
    Spline_Hermite::sample3D_deriv( *(Vec3d*)g0, *(Vec3d*) dg, *(Vec3i*)ng, (Quat4d*)Eg, (Quat4d*)dEg, n, (Vec3d*)ps, (Quat4d*)fes ); 
}

void sample_SplineHermite2D_deriv( double* g0, double* dg, int* ng, double* Eg, double* dEg, int n, double* ps, double* fes ){
    Spline_Hermite::sample2D_deriv( *(Vec2d*)g0, *(Vec2d*)dg, *(Vec2i*)ng, (Quat4d*)Eg, (Quat4d*)dEg, n, (Vec2d*)ps, (Quat4d*)fes );  
}

void sample_SplineHermite1D_deriv( double g0, double dg, int ng, double* EFg, int n, double* ps, double* fes ){
    Spline_Hermite::sample1D_deriv( g0, dg, ng, (Vec2d*)EFg, n, ps, (Vec2d*)fes );  
}



inline double periodic_x2(double x, int nmax=100){
    //int ix    = (int)(x+nmax) - nmax; // fast floor
    int ix=(int)x; if(x<0)ix--;
    double dx = (x-ix-0.5)*2.0;
    double y = 1-dx*dx;
    if( ix&1 ){ y*=-1; } 
    return y;
}


void sample_func( int n, double* xs, double* ys, int kind ){
    for(int i=0; i<n; i++ ){
        double x = xs[i];
        switch (kind){
            case 0:  { ys[i]= sin(x);         } break;
            case 1:  { ys[i]= periodic_x2(x); } break;
            default: { ys[i]= NAN;            } break;
        }
    }
}


void sample_SplineConstr( double x0, double dx, int np, double* Eps, int n, double* xs, double* Es, double* Fs ){
    SplineConstr C( {0,1}, x0, dx, np, Eps );
    Vec3d ps[2]{{.0,.0,.0},{.0,.0,.0}};
    Vec3d fs[2];
    for(int i=0; i<n; i++ ){
        ps[1]={xs[i],0.0,0.0};
        fs[0]=Vec3dZero;
        fs[1]=Vec3dZero;
        Es[i] = C.apply( ps, fs );
        Fs[i] = fs[0].x;
    }
}

void sample_DistConstr( double lmin, double lmax, double kmin, double kmax, double flim , int n, double* xs, double* Es, double* Fs ){
    DistConstr C( {0,1}, {lmax,lmin}, {kmax,kmin}, flim  );
    Vec3d ps[2]{{.0,.0,.0},{.0,.0,.0}};
    Vec3d fs[2];
    for(int i=0; i<n; i++ ){
        ps[1]={xs[i],0.0,0.0};
        fs[0]=Vec3dZero;
        fs[1]=Vec3dZero;
        Es[i] = C.apply( ps, fs );
        Fs[i] = fs[0].x;
    }
}

void sample_evalPiAling( double k, double ang0, double r1, double r2, int n, double* angles, double* Es, double* Fs ){
    Vec3d h1={1,0,0};
    Vec3d f1,f2;
    double c0 = cos(ang0);
    for(int i=0; i<n; i++ ){
        double a = angles[i];
        Vec3d h2={cos(a),sin(a),0.0};
        Es[i] = evalPiAling( h1, h2, 1./r1, 1./r2, k, f1, f2 );  
        Fs[i] = f1.y;
    }
}

void sample_evalAngleCos( double k, double ang0, double r1, double r2, int n, double* angles, double* Es, double* Fs ){
    Vec3d h1={1,0,0};
    Vec3d f1,f2;
    double c0 = cos(ang0);
    for(int i=0; i<n; i++ ){
        double a = angles[i];
        Vec3d h2={cos(a),sin(a),0.0};
        Es[i] = evalAngleCos( h1, h2, 1./r1, 1./r2, k, c0, f1, f2 );
        Fs[i] = f1.y;
    }
}

void sample_evalAngleCosHalf( double k, double ang0, double r1, double r2, int n, double* angles, double* Es, double* Fs ){
    Vec3d h1={1,0,0};
    Vec3d f1,f2;
    Vec2d cs0; cs0.fromAngle( ang0/2. );
    for(int i=0; i<n; i++ ){
        double a = angles[i];
        Vec3d h2={cos(a),sin(a),0.0};
        Es[i] = evalAngleCosHalf( h1, h2, 1./r1, 1./r2, cs0, k, f1, f2 );
        Fs[i] = f1.y;
    }
}

void sampleNonBond(int n, double* rs, double* Es, double* fs, int kind, double*REQi_,double*REQj_, double K, double Rdamp ){
    Quat4d REQi = *(Quat4d*)REQi_;
    Quat4d REQj = *(Quat4d*)REQj_;
    Quat4d REQij; combineREQ( REQi, REQj, REQij );
    REQij.y = sqrt(REQij.y);
    Vec3d pi=Vec3dZero;
    Vec3d pj=Vec3dZero;
    double R2damp=Rdamp*Rdamp;
    for(int i=0; i<n; i++){
        double E;
        Vec3d  f=Vec3dZero;
        pj.x=rs[i];
        switch(kind){
            case 1: E=addAtomicForceMorseQ( pj-pi, f, REQij.x, REQij.y, REQij.z, K, R2damp ); break;  // Morse
            case 2: E=addAtomicForceLJQ   ( pj-pi, f, REQij );                                break;  // Lenard Jones
            case 3: double fr; E=erfx_e6( pj.x, K, fr ); f.x=fr; break;  // gauss damped electrostatics
            case 4: E=repulsion_R4( pj-pi, f, REQij.x-Rdamp, REQij.x, K );
        }
        //printf( "i %i r %g E %g f %g \n", i, pj.x, E, f.x );
        fs[i]=f.x;
        Es[i]=E;
    }
}

void sampleNonBondTypes( int n, double* rs, double* Es, double* fs, int kind, double qH, double qX, double K, double Rdamp, double dcomp, char* type_str ){
    char nameH[16];
    char nameX[16];
    sscanf( type_str, "%s %s", nameH, nameX );
    //int iH = W.params.getAtomType( nameH );
    //int iX = W.params.getAtomType( nameX );
    const AtomType* tH = W.params.getAtomTypeObj( nameH );
    const AtomType* tX = W.params.getAtomTypeObj( nameX );
    int itE = tX->ePairType;
    if( (itE<0)||(itE>W.params.atypes.size())){ printf("ERROR: type(%s).ePairType=%i => exit()\n", nameX, itE ); exit(0); };
    const AtomType* tE = &W.params.atypes[ tX->ePairType ];
    Quat4d REQi = tH->assignREQH();   REQi.z=qH;
    Quat4d REQj = tX->assignREQH();   REQj.z=qX;
    Quat4d REQe = tE->assignREQH();
    Quat4d REQdi{ 1.0,0.0,-REQi.z, 0.0 };
    Quat4d REQdj{ 1.0,0.0,-REQj.z, 0.0 };
    REQj.z -= tE->Qbase;
    double de = tE->Ruff;
    //printf( "REQ_H %-8s R=%5.3f E=%7.5f Q=%6.3f H=%6.3f \n",  tH->name, REQi.x,REQi.y,REQi.z,REQi.w );
    //printf( "REQ_X %-8s R=%5.3f E=%7.5f Q=%6.3f H=%6.3f \n",  tX->name, REQj.x,REQj.y,REQj.z,REQj.w );
    //printf( "REQ_E %-8s R=%5.3f E=%7.5f Q=%6.3f H=%6.3f \n",  tE->name, REQe.x,REQe.y,REQe.z,REQe.w );
    //sprintf( s, "%s[%i]-%s[%i] (%4.2f,%5.4f,%4.2f,%4.2f) (%4.2f,%5.4f,%4.2f,%4.2f)", W.params.atypes[ W.ffl.atypes[b.x]].name, b.x, W.params.atypes[ W.ffl.atypes[b.y]].name, b.y, REQi.x,REQi.y,REQi.z,REQi.w,  REQj.x,REQj.y,REQj.z,REQj.w  );
    double R2damp=Rdamp*Rdamp;
    Vec3d d{dcomp,0.0,0.0};
    for(int i=0; i<n; i++){
        double E = 0;
        Vec3d  f = Vec3dZero;
        double r = rs[i];
        Quat4d REQ;
        switch(kind){
            case 1:{
                Vec3d f_=Vec3dZero;
                combineREQ( REQi, REQj,  REQ ); E += getLJQH( Vec3d{0.0,0.0,r         }, f_, REQ, R2damp );  f.add(f_);
                combineREQ( REQi, REQdj, REQ ); E += getLJQH( Vec3d{0.0,0.0,r+dcomp   }, f_, REQ, R2damp );  f.add(f_);
                combineREQ( REQj, REQdi, REQ ); E += getLJQH( Vec3d{0.0,0.0,r+dcomp   }, f_, REQ, R2damp );  f.add(f_);
                combineREQ( REQdi,REQdj, REQ ); E += getLJQH( Vec3d{0.0,0.0,r+dcomp*2 }, f_, REQ, R2damp );  f.add(f_);

                combineREQ( REQe, REQdi, REQ ); E += getLJQH( Vec3d{0.0,0.0,r+dcomp-de}, f_, REQ, R2damp );  f.add(f_);
                combineREQ( REQe, REQi,  REQ ); E += getLJQH( Vec3d{0.0,0.0,r-de      }, f_, REQ, R2damp );  f.add(f_);

                } break;
        }
        //printf( "i %i r %g E %g f %g \n", i, pj.x, E, f.x );
        fs[i]=f.x;
        Es[i]=E;
    }
    //return nb;
}

int selectBondsBetweenTypes( int imin, int imax, int it1, int it2, bool byZ, bool bOnlyFirstNeigh, int* atoms_ ){
    W.builder.selectBondsBetweenTypes( imin, imax, it1, it2, byZ, bOnlyFirstNeigh );
    Vec2i* atoms = (Vec2i*)atoms_;
    int i=0;
    for( int ib : W.builder.selection ){
        Vec2i b = W.builder.bonds[ib].atoms;
        //int t1 = W.builder.atoms[b.a].type;
        int t2 = W.builder.atoms[b.b].type;
        //if( byZ ){ t1=W.params.atypes[t1].iZ; t2=W.params.atypes[t2].iZ; };
        if( byZ ){ t2=W.params.atypes[t2].iZ; };
        //printf( "b(%i,%i) types(%i,%i) its(%i,%i)\n",  b.a,b.b,  t1,t2, it1,it2 );
        //if(t2==it1){ b=b.swaped(); printf( "bond swap b(%i,%i)\n", b.a, b.b ); };
        if(t2==it1){ b=b.swaped(); };
        //printf( "selectBondsBetweenTypes[%i] %i,%i \n", i, b.a,b.b );
        atoms[i]=b;
        i++;
    }
    return i;
}

int getFrament( int ifrag, int* bounds_, double* pose ){
    const MM::Fragment& frag = W.builder.frags[ifrag];
    if(bounds_){
        Vec2i* bounds=(Vec2i*)bounds_;
        bounds[0] = frag.atomRange;
        bounds[1] = frag.confRange;
        bounds[2] = frag.bondRange;
    }
    if( pose ){
        const double* ds = &(frag.pos.x);
        for(int i=0; i<7; i++){ pose[i] = ds[i]; }
    }
    return frag.imolType;
    //angRange;
    //dihRange;
}

void findMainAxes( double* rot, int ifrag, int imin, int imax, int* permut_, bool bRot){
    Vec3i permut{2,1,0};
    if(permut_){ permut=*((Vec3i*)permut_); };
    if(ifrag>=0){ const MM::Fragment& frag = W.builder.frags[ifrag]; imin=frag.atomRange.a; imax=frag.atomRange.b; }
    if(rot){ *((Mat3d*)rot)= W.builder.findMainAxes(imin,imax,true,bRot,permut); }
    else   {                 W.builder.findMainAxes(imin,imax,true,bRot,permut); }    
}

void findSymmetry( int* found, int ifrag, int imin=0,int imax=-1, double tol=0.1 ){
    if(ifrag>=0){ const MM::Fragment& frag = W.builder.frags[ifrag]; imin=frag.atomRange.a; imax=frag.atomRange.b; }
    W.builder.findSymmetry( found, imin,imax, tol );
}

void orient( const char* fname, int fw1,int fw2,  int up1,int up2,  int i0,  int imin, int imax ){
    FILE* fout = fopen(fname,"w");    
    Vec3d fw = W.builder.atoms[fw2].pos - W.builder.atoms[fw1].pos;  fw.normalize();
    Vec3d up = W.builder.atoms[up2].pos - W.builder.atoms[up1].pos;  up.normalize();
    W.builder.orient_atoms( fw, up, W.builder.atoms[i0].pos,  Vec3dZero,   imin,imax );
    W.builder.write2xyz(fout, "#scanHBond[0]" );
    fclose(fout);
}

//

void setSwitches( int CheckInvariants, int PBC, int NonBonded, int MMFF, int Angles, int PiSigma, int PiPiI ){
    #define _setbool(b,i) { if(i>0){b=true;}else if(i<0){b=false;} }
    _setbool( W.bCheckInvariants, CheckInvariants  );
    _setbool( W.bPBC         , PBC       );
    _setbool( W.bNonBonded   , NonBonded );
    _setbool( W.bMMFF        , MMFF      );
    _setbool( W.ffl.doAngles , Angles    );
    _setbool( W.ffl.doPiSigma, PiSigma   );
    _setbool( W.ffl.doPiPiI  , PiPiI     );
    W.ffl.bSubtractAngleNonBond = W.bNonBonded;
    #undef _setbool
}

void setOptLog( int n, double* cos, double* f, double* v, double* dt, double* damp ){
    W.opt_log.n    = n;
    W.opt_log.cos  = cos;
    W.opt_log.f    = f;
    W.opt_log.v    = v;
    W.opt_log.dt   = dt;
    W.opt_log.damp = damp;
}

} // extern "C"
