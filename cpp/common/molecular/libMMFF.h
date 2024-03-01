
#include "Forces.h"
#include "InterpolateTricubic.h"

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

void sample_SplineHermite( double g0, double dg, int ng, double* Eg, int n, double* xs, double* Es, double* Fs ){
    Spline_Hermite::sample1D( g0,dg,ng,Eg, n, xs, Es, Fs );
}

void sample_SplineHermite2D( double* g0, double* dg, int* ng, double* Eg, int n, double* ps, double* fes ){
    long t0 = getCPUticks();
    //Spline_Hermite::sample2D( *(Vec2d*)g0, *(Vec2d*) dg, *(Vec2i*)ng, Eg, n, (Vec2d*)ps, (Vec3d*)fes );         //    82 [tick/point] with -Ofast
    Spline_Hermite::sample2D_avx( *(Vec2d*)g0, *(Vec2d*) dg, *(Vec2i*)ng, Eg, n, (Vec2d*)ps, (Vec3d*)fes );   //   111 [tick/point] with -Ofast  
    double t = (getCPUticks()-t0); printf( "sample2D(n=%i) time=%g[kTick] %g[tick/point]\n", n, t*(1.e-3), t/n );
}

void sample_SplineHermite3D( double* g0, double* dg, int* ng, double* Eg, int n, double* ps, double* fes ){
    Spline_Hermite::sample3D( *(Vec3d*)g0, *(Vec3d*) dg, *(Vec3i*)ng, Eg, n, (Vec3d*)ps, (Quat4d*)fes );
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
