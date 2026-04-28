#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <vector>
#include <unordered_map>
#include <string>

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
//#include "VecN.h"

#include "testUtils.h"
#include "libUtils.h"

//#include "integration.h"
//#include "AOIntegrals.h"
//#include "AOrotations.h"

Vec3d DEBUG_dQdp;
int DEBUG_iter     = 0;
int DEBUG_log_iter = 0;
int i_DEBUG=0;

#include <globals.h>
//int verbosity=0;
//int idebug   =0;

const char* prefix = "#Epiece";

#include "InteractionsGauss.h"
#include "eFF.h"
#include "DynamicOpt.h"
#include "Vec3Utils.h"

#define WITH_BUILDER 1
#if WITH_BUILDER
#include "MMFFparams.h"
#include "MMFFBuilder.h"
#include "Atoms.h"
MMFFparams params;
std::vector<Atoms*> samples;
MM::Builder builder;
#endif

// ============ Global Variables

EFF         ff;
DynamicOpt opt;

bool opt_initialized=false;

char* trj_fname = 0;
int   savePerNsteps = 1;

bool isPreAllocated=false;
//bool bChangeCore=true;

int measure_i0=0;
int measure_n=0;
int nDist = 0;
int nAng  = 0;
Vec2i* dist_inds   = 0;
Vec3i* ang_inds    = 0;
double* dists_vals = 0;
double* ang_vals   = 0;

// initialize measurement parameters: start index, count, and pointers
void setup_measurements(int i0_, int n_, int nd, Vec2i* d_inds, int na_, Vec3i* a_inds, double* d_vals, double* a_vals){
    measure_i0 = i0_;
    measure_n  = n_;
    nDist      = nd;
    dist_inds  = d_inds;
    nAng       = na_;
    ang_inds   = a_inds;
    dists_vals = d_vals;
    ang_vals   = a_vals;
}

extern "C"{
// ========= Grid initialization

void setVerbosity( int verbosity_, int idebug_ ){
    // Disable buffering completely on stdout and stderr
    setvbuf(stdout, NULL, _IONBF, 0);
    setvbuf(stderr, NULL, _IONBF, 0);
    verbosity = verbosity_;
    idebug    = idebug_;
}

double eval(){ return ff.eval(); };

void evalFuncDerivs( int ie, int n, double* r, double* s, double* Es, double* Fr, double* Fs ){
    double fr,fs;
    for(int i=0; i<n; i++){
        ff.epos [ie].x = r[i];
        ff.esize[ie]   = s[i];
        //ff.esize[1]   = s[i]; 
        Es[i]         = ff.eval();
        Fr[i]         = ff.eforce[ie].x; 
        Fs[i]         = ff.fsize [ie];
        //Fs[i]         = ff.fsize [1];
        printf("electron[%i] x %g s %g Es %g Fr %g Fs %g \n", ie, r[i], s[i], Es[i], Fr[i], Fs[i] );
        //printf( "[%i] Es %g Fr %g Fs %g \n", i, Es[i], ff.eforce[0].x, ff.fsize[0] );
    }
}

void evalNumDerivs( double* Fnum, double d ){
    double denom=1/(2.*d);
    for(int i=0; i<ff.nDOFs; i++){
        double o = ff.pDOFs[i];
        ff.pDOFs[i]=o-d;  double E1=ff.eval();
        ff.pDOFs[i]=o+d;  double E2=ff.eval();  
        Fnum[i] = (E2-E1)*denom;
        ff.pDOFs[i]=o;
    }
}

void init_buffers(){
    // atoms (ions)
    buffers.insert( { "pDOFs",           ff.pDOFs  } );
    buffers.insert( { "fDOFs",           ff.fDOFs  } );
    buffers.insert( { "apos",   (double*)ff.apos   } );
    buffers.insert( { "aforce", (double*)ff.aforce } );
    buffers.insert( { "epos",   (double*)ff.epos   } );
    buffers.insert( { "eforce", (double*)ff.eforce } );
    buffers.insert( { "esize",  (double*)ff.esize  } );
    buffers.insert( { "fsize",  (double*)ff.fsize  } );
    buffers.insert ( { "aPars", (double*)ff.aPars  } );
    buffers.insert ( { "Es",            &ff.Etot   } );
    ibuffers.insert( { "espin",          ff.espin  } );
    ibuffers.insert( { "ndims",         &ff.ne     } );
    if(ff.vDOFs){
        buffers.insert( { "vDOFs",              ff.vDOFs  } );
        buffers.insert( { "avel",      (double*)ff.avel   } );
        buffers.insert( { "evel",      (double*)ff.evel   } );
        buffers.insert( { "vsize",              ff.vsize  } );
        buffers.insert( { "invMasses",          ff.invMasses           } );
        buffers.insert( { "invAmass",           ff.invMasses           } );
        buffers.insert( { "invEmass",          (ff.invMasses+ff.na*3)         } );
        buffers.insert( { "invSmass",          (ff.invMasses+(ff.na+ff.ne)*3) } );
    }
}


void printShortestBondLengths(){
    for(int ia=0; ia<ff.na; ia++){
        int ifound; double r2 = findNearest( ff.na, ff.apos, ff.apos[ia], ifound, ia );
        printf( "bond[%i,%i] L=%g [A]\n", ia, ifound, sqrt(r2) ); 
    }
}

void setTrjName( char* trj_fname_, int savePerNsteps_ ){ trj_fname=trj_fname_; if(verbosity>0)printf( "setTrjName(%s)\n", trj_fname ); savePerNsteps=savePerNsteps_;  }

bool load_xyz( const char* fname ){ 
    //printf( "load_xyz \n" );
    bool b = ff.loadFromFile_xyz( fname );
    init_buffers();
    return b; 
}

bool load_fgo( const char* fname, bool bVel, double fUnits ){ 
    //printf( "load_xyz \n" );
    bool b = ff.loadFromFile_fgo( fname, bVel, fUnits );
    init_buffers();
    return b; 
}

void init( int na, int ne ){
    ff.realloc ( na, ne, true );
    //ff.autoAbWs( default_aAbWs, default_eAbWs );
    init_buffers();
}

void info(){ ff.info(); }
void printSwitches(){ ff.printSwitches(); }

double* getEnergyPointer(){ return &ff.Ek; }
int*    getDimPointer   (){ return &ff.ne; }

void initOpt( double dt, double damping, double f_limit, bool bMass ){
    if(ff.vDOFs){ if(bMass){ff.makeMasses(ff.invMasses);}else{ff.makeMasses(ff.invMasses,1.0);} }
    printf("initOpt()  ff.nDOFs %i, ff.pDOFs %p, ff.vDOFs %p, ff.fDOFs %p, ff.invMasses %p\n", ff.nDOFs, ff.pDOFs, ff.vDOFs, ff.fDOFs, ff.invMasses ); 
    opt.bindOrAlloc( ff.nDOFs, ff.pDOFs, ff.vDOFs, ff.fDOFs, ff.invMasses );
    //opt.cleanVel( ); // this is already inside initOpt
    //opt.initOpt( dt, damping );
    opt.setTimeSteps(dt );
    opt.setDamping  (damping);
    opt.f_limit = f_limit;
    opt_initialized=true;
    opt.fixmask = ff.fixmask;
    opt.bfixmask=true;
};

//int run( int nstepMax, double dt, double Fconv=1e-6, int ialg=0, double* outE, double* outF ){ 
int run( int nstepMax, double dt, double Fconv, int ialg, double* outE, double* outF ){ 
    double F2conv=Fconv*Fconv;
    double F2 = 1.0;
    double Etot;
    int itr=0;
    if( (ialg!=0)&(!opt_initialized) ){ printf("ERROR ialg(%i)>0 but optimizer not initialized => call initOpt() first !"); exit(0); };
    opt.setTimeSteps(dt);
    //printf( "verbosity %i nstepMax %i Fconv %g dt %g \n", verbosity, nstepMax, Fconv, dt );
    //printf( "trj_fname %s \n", trj_fname );
    if(ialg>0){
        //printf("run() opt.vel %p\n", opt.vel); 
        opt.cleanVel( ); 
    }
    bool bConv=false;
    //if(ff.nfix>0){ ff.apply_hard_fix(); }
    for(itr=0; itr<nstepMax; itr++ ){
        ff.clearForce();
        Etot = ff.eval();
        if( ff.bNegativeSizes & (verbosity>0) ){ printf( "negative electron sizes in step #%i => perhaps decrease relaxation time step dt=%g[fs]? \n", itr, opt.dt ); }
        //if(ff.nfix>0){ if(ialg>0){ ff.clear_fixed_dynamics(); }else{ ff.clear_fixed_force(); } }
        switch(ialg){
            case  0: F2 = ff .move_GD      (dt);      break;
            case -1: F2 = opt.move_LeapFrog(dt);      break;
            case  1: F2 = opt.move_MD (dt,opt.damping); break;
            case  2: F2 = opt.move_FIRE();       break;
        }
        //if(ff.nfix>0){ ff.apply_hard_fix(); }
        if(outE){ outE[itr]=Etot;     }
        if(outF){ outF[itr]=sqrt(F2); }
        if(verbosity>2){ printf("itr: %6i Etot[eV] %16.8f |F|[eV/A] %16.8f \n", itr, Etot, sqrt(F2) ); };
        if(F2<F2conv){
            //if(verbosity>0){ printf("Converged in %i iteration Etot %g[eV] |F| %g[eV/A] <(Fconv=%g) \n", itr, Etot, sqrt(F2), Fconv ); };
            bConv=true;
            break;
        }
        if( (trj_fname) && (itr%savePerNsteps==0) )  ff.save_xyz( trj_fname, "a" );
    }
    if(verbosity>1){ printf("run() %s in %6i iterations Etot[eV] %16.8f |F|[eV/A] %16.8f (Fconv=%g) \n", bConv ? "    CONVERGED" : "NOT-CONVERGED", itr, Etot, sqrt(F2), Fconv ); };
    //printShortestBondLengths();
    return itr;
}

void set_constrains( int nfix, Quat4d* fixed_poss, Vec2i* fixed_inds, bool bRealloc=true  ){
    // Set following variables of eFF object
    // int nfix=0;
    // Quat4d* fixed_poss=0; // [nfix], {x,y,z,w} position of fixed particles
    // Vec2i*  fixed_inds=0; // [nfix]  {ia/-ie, bitmask{x|y|z|w}}, ia/-ie is index of atom/ or negative index of electron, bitmax indicate which component is fixed
    
    if(bRealloc) ff.realloc_fixed(nfix);
    for(int i=0; i<nfix; i++){
        ff.fixed_poss[i] = fixed_poss[i];
        ff.fixed_inds[i] = fixed_inds[i];
    }

    // for(int i=0; i<nfix; i++){
    //     Quat4i ip = fixed_poss[i];
    //     if(fixed_inds[i].x>=0){
    //         ff.apos_fix[i] = ip.f;
    //     }else{
    //         ff.epos_fix [i].x = ip.f;
    //         ff.esize_fix[i] = fixed_inds[i].e;
    //     }
    // }
}

void relaxed_scan( int nconf, int nfix, double* fixed_poss, int* fixed_inds_, double* outEs, double* apos_, double* epos_, int nstepMax, double dt, double Fconv, int ialg, char* scan_trj_name ){
    //printf( "relaxed_scan() nconf %i nfix %i @fixed_poss=%p @fixed_inds=%p @outEs=%p @apos_=%p @epos_=%p @nstepMax %i @dt %g @Fconv %g @ialg %i \n", nconf, nfix, fixed_poss, fixed_inds_, outEs, apos_, epos_, nstepMax, dt, Fconv, ialg );
    //printf( "relaxed_scan() ialg: %i @vsize=%p @evel=%p @avel=%p \n", ialg, ff.vsize, ff.evel, ff.avel );
    Vec2i* fixed_inds = (Vec2i*)fixed_inds_;
    //for(int i=0; i<nfix; i++){ printf( "fixed[%3i] %3i %3i \n", i, fixed_inds[i].x, fixed_inds[i].y ); }
    //exit(0);   
    ff.realloc_fixed(nfix);
    for(int iconf=0; iconf<nconf; iconf++){
        //if((verbosity>0)&&(nstepMax>0))printf( "relaxed_scan() iconf %3i \n", iconf );
        if(verbosity>0)printf( "--- relaxed_scan() iconf %3i \n", iconf );
        if( nfix>0 ){ set_constrains( nfix, ((Quat4d*)fixed_poss)+iconf*nfix, fixed_inds, false ); }
        if( nstepMax > 0 ){
            run( nstepMax, dt, Fconv, ialg, 0, 0 );
        }else{
            ff.apply_hard_fix();
            ff.eval();
        }
        if(outEs){        
            // double Etot=0,Ek=0, Eee=0,EeePaul=0,EeeExch=0,  Eae=0,EaePaul=0,  Eaa=0; 
            //printf( "Etot %g Ek %g Eee %g EeePaul %g EeeExch %g Eae %g EaePaul %g Eaa %g \n", ff.Etot, ff.Ek, ff.Eee, ff.EeePaul, ff.EeeExch, ff.Eae, ff.EaePaul, ff.Eaa );
            double* Eis = outEs + iconf*8;
            Eis[0] = ff.Etot;
            Eis[1] = ff.Ek;
            Eis[2] = ff.Eee;
            Eis[3] = ff.EeePaul;
            Eis[4] = ff.EeeExch;
            Eis[5] = ff.Eae;
            Eis[6] = ff.EaePaul;
            Eis[7] = ff.Eaa;
        }
        if( apos_ ){
            Vec3d*  apos = ((Vec3d* )apos_)+iconf*ff.na;
            for(int j=0; j<ff.na; j++){ apos[j]   = ff.apos[j]; }
        }
        if( epos_ ){
            Quat4d* epos = ((Quat4d*)epos_)+iconf*ff.ne;
            for(int j=0; j<ff.ne; j++){ epos[j].f = ff.epos[j]; epos[j].w = ff.esize[j]; }
        }
        // record distances and angles in specified interval
        if(dist_inds){ int ic=iconf-measure_i0; if(ic>=0 && ic<measure_n){
            ff.analyse_distances( dist_inds, nDist,  dists_vals + ic*nDist);
        }}
        if(ang_inds){ int ic=iconf-measure_i0; if(ic>=0 && ic<measure_n){
            ff.analyse_angles( ang_inds, nAng, ang_vals + ic*nAng);
        }}
        if( scan_trj_name )  ff.save_xyz( scan_trj_name, "a" );
        fflush(stdout); 
    }
}


void sample_ee( int n, double* RSs_, double* FEout_, int spin, double* KRSrho_, bool bEvalCoulomb, bool bEvalPauli, int iPauliModel ){
    Quat4d* FEout =(Quat4d*)FEout_;
    Vec3d* RSs    =(Vec3d*)RSs_;
    Vec3d  KRSrho =*(Vec3d*)KRSrho_;
    //using namespace std;
    //auto [x, y, z] = RSs[0];
    for(int i=0; i<n; i++){
        double ri = RSs[i].x;
        double si = RSs[i].y;
        double sj = RSs[i].z;
        Vec3d f=Vec3dZero;
        Quat4d EFi=Quat4dZero;
        Vec3d dR{0.0,0.0,ri};
        if(bEvalCoulomb){
            EFi.w += addCoulombGauss( dR, si, sj, f, EFi.y, EFi.z, 1.0 );
        }
        if(bEvalPauli){
            // Pauli repulsion form this eFF paper http://aip.scitation.org/doi/10.1063/1.3272671  
            // iPauliModel==0 Pauli repulasion from Valence-Bond theory  // Pauli repulsion only for electrons with same spin
            // iPauliModel==0 Pauli repulasion as overlap between same spin orbitals 
            if     ( iPauliModel == 1 ){              EFi.w += addPauliGauss_New    ( dR, si, sj, f, EFi.y, EFi.z, spin, KRSrho ); }
            else if( iPauliModel == 2 ){  if(spin>0){ EFi.w += addPauliGaussVB      ( dR, si, sj, f, EFi.y, EFi.z ); } }
            else if( iPauliModel == 0 ){  if(spin>0){ EFi.w += addDensOverlapGauss_S( dR, si*M_SQRT2, sj*M_SQRT2, ff.KPauliOverlap, f,  EFi.y, EFi.z ); } }
       }
       EFi.x=f.x;
       FEout[i]=EFi;
    }
}

void sample_EA( int n, double* RSs_, double* FEout_, double* KRSrho_,  double* aPar_,  bool bEvalAECoulomb, bool bCoreCoul, bool bEvalAEPauli ){
    Vec3d* FEout  =(Vec3d*)FEout_;
    Vec2d* RSs    =(Vec2d*)RSs_;
    Vec3d  KRSrho =*(Vec3d*)KRSrho_;
    Quat4d  aPar  =*(Quat4d*)aPar_;
    for(int i=0; i<n; i++){
        double ri = RSs[i].x;
        double si = RSs[i].y;
        Vec3d f=Vec3dZero;
        Vec3d EFi=Vec3dZero;
        Vec3d dR{0.0,0.0,ri};
        double fs_junk=0;
        if(bEvalAECoulomb){ EFi.z += addCoulombGauss  ( dR,  aPar.y, si, f, fs_junk, EFi.y, aPar.x );  }
        if( bEvalAEPauli ){ EFi.z += addPauliGauss_New( dR, si, aPar.z,  f, EFi.y, fs_junk, 0, KRSrho, aPar.w       );    
            if(bCoreCoul ){ EFi.z += addCoulombGauss  ( dR, si, aPar.z,  f, EFi.y, fs_junk,            aPar.w*2.0  ); }
        }
        EFi.x=f.x;
        FEout[i]=EFi;
    }

}

void save_fgo( char const* filename, bool bVel, bool bAppend ){ 
    if(bAppend){ ff.writeTo_fgo( filename, bVel, "a" ); }
    else       { ff.writeTo_fgo( filename, bVel, "w" ); } 
};

void save_xyz( char const* filename, bool bAppend ){ 
    if(bAppend){ ff.save_xyz( filename, "a" ); }
    else       { ff.save_xyz( filename, "w" ); } 
};


//#define NEWBUFF(name,N)   double* name = new double[N]; buffers.insert( {#name, name} );

void setPauliModel(int i  ){ ff.iPauliModel = i; }
void setKPauli( double KPauli ){
    if( ff.iPauliModel == 0 ){ ff.KPauliOverlap=KPauli; }
    else                  { ff.KPauliKin=KPauli;     }
}

void setSwitches( int bEvalKinetic, int bEvalCoulomb, int  bEvalPauli, int bEvalAA, int bEvalAE, int bEvalAECoulomb, int bEvalAEPauli, int bCoreCoul, int bEvalCoreCorect ){
    //printf( "\n\n\n\n#### setSwitches_ bEvalAEPauli %i \n", bEvalAEPauli );
#define _setbool(b,i) { if(i>0){b=true;}else if(i<0){b=false;} }
    //_setbool( ff.bNormalize     , bNormalize     );
    //_setbool( ff.bNormForce     , bNormForce     );
    _setbool( ff.bEvalKinetic   , bEvalKinetic   );
    _setbool( ff.bEvalCoulomb   , bEvalCoulomb   );
    //_setbool( ff.bEvalExchange  , bEvalExchange  );
    _setbool( ff.bEvalPauli     , bEvalPauli     );
    _setbool( ff.bEvalAA        , bEvalAA        );
    _setbool( ff.bEvalAE        , bEvalAE        );
    _setbool( ff.bEvalAECoulomb , bEvalAECoulomb );
    _setbool( ff.bEvalAEPauli   , bEvalAEPauli   );
    _setbool( ff.bCoreCoul      , bCoreCoul      );
    _setbool( ff.bEvalCoreCorect, bEvalCoreCorect);
#undef _setbool
}

void setup( int isetup ){
    switch(isetup){
        case 0:
            ff.bEvalAECoulomb = true;
            ff.bEvalAEPauli   = true;
            ff.iECPmodel      = 0;
            ff.iPauliModel    = 0;
            break;
        case 1:
            ff.bEvalAECoulomb = true;
            ff.bEvalAEPauli   = true;
            ff.iECPmodel      = 1;
            ff.iPauliModel    = 0;
            break;
        default:
            ff.bEvalAECoulomb = true;
            ff.bEvalAEPauli   = true;
            ff.iECPmodel      = 2;
            ff.iPauliModel    = 0;
            break;
    }
}

void setAtomParams( int n, const double* params, bool bCopy=true, int mode=1 ){ 
    if(bCopy){ // if we worry that the source pointer may be deleted
        if(mode==1){
            Quat4d* atom_params = new Quat4d[n];
            for(int i=0; i<n; i++){ atom_params[i] = ((Quat4d*)params)[i]; }
            ff.atom_params = atom_params;
        }else if(mode==2){
            double8* atom_params = new double8[n];
            for(int i=0; i<n; i++){ atom_params[i] = ((double8*)params)[i]; }
            ff.atom_params2 = atom_params;
        }
    }else{
        if     (mode==1){ ff.atom_params  = (Quat4d* )params; }
        else if(mode==2){ ff.atom_params2 = (double8*)params; }
    }
}


int preAllocateFGO(const char* fname, double fUnits=1.){
    if(!ff.loadFromFile_fgo(fname,true,fUnits)) return -1;
    //init_buffers();
    printf("preAllocateFGO()  ff.nDOFs %i, ff.pDOFs %p, ff.vDOFs %p, ff.fDOFs %p, ff.invMasses %p\n", ff.nDOFs, ff.pDOFs, ff.vDOFs, ff.fDOFs, ff.invMasses ); 
    opt.bindOrAlloc(ff.nDOFs, ff.pDOFs, ff.vDOFs, ff.fDOFs, ff.invMasses);
    isPreAllocated=true;
    return 0;
}

int processFGO(const char* fname, double fUnits=1., double* outEs=0, double* apos_=0, Quat4d* epos_=0, int nstepMax=1000, double dt=0.001, double Fconv=1e-3, int ialg=2, bool bOutXYZ=false, bool bOutFGO=false){
    FILE* f = fopen(fname,"r"); if(!f) return 0;
    //init_buffers();
    int iconf=0;
    while(true){
        if( ff.readSingleFGO(f,true,fUnits, !isPreAllocated) < 2 ) break;
        if(!isPreAllocated){ 
            printf("processFGO()  ff.nDOFs %i, ff.pDOFs %p, ff.vDOFs %p, ff.fDOFs %p, ff.invMasses %p\n", ff.nDOFs, ff.pDOFs, ff.vDOFs, ff.fDOFs, ff.invMasses ); 
            opt.bindOrAlloc(ff.nDOFs, ff.pDOFs, ff.vDOFs, ff.fDOFs, ff.invMasses); isPreAllocated=true; 
        }
        run(nstepMax,dt,Fconv,ialg,0,0);
        ff.eval();
        if(verbosity>0){ printf( "processFGO() iconf %i Etot %g Ek %g Eee %g Eae %g Eaa %g \n", iconf, ff.Etot, ff.Ek, ff.Eee, ff.Eae, ff.Eaa ); }
        if(verbosity>1){ ff.info(); }
        if(bOutXYZ) ff.save_xyz("processFGO.xyz","a");
        if(bOutFGO) ff.writeTo_fgo("processFGO.fgo", false, "a", iconf);
        if(apos_){ Vec3d*  ap=((Vec3d* )apos_)+iconf*ff.na; for(int i=0;i<ff.na;i++){ ap[i]=ff.apos[i];                        }     }
        if(epos_){ Quat4d* ep=((Quat4d*)epos_)+iconf*ff.ne; for(int i=0;i<ff.ne;i++){ ep[i].f=ff.epos[i]; ep[i].w=ff.esize[i]; }     }
        if(outEs){ double* outEi=outEs+iconf*5; outEi[0]=ff.Etot; outEi[1]=ff.Ek; outEi[2]=ff.Eee; outEi[3]=ff.Eae; outEi[4]=ff.Eaa; }
        iconf++;
    }
    fclose(f);
    return iconf;
}

// #define WITH_BUILDER 1
#if WITH_BUILDER


void builder2EFF(EFF& ff, const MM::Builder& builder, bool bRealloc=true ){
    int ne=0;
    int na=0;
    int epairtyp = params.getAtomType("E");
    // --- loop over atoms and electron pairs
    for(int i=0; i<builder.atoms.size(); i++){
        int it = builder.atoms[i].type;
        if( it!=epairtyp ){
            int iZ = params.atypes[it].iZ;
            if(!bRealloc){
                //ff.aPars [na] = EFF::default_AtomParams[iZ];
                ff.aPars [na] = ff.atom_params[iZ];
                ff.apos  [na] = builder.atoms[i].pos;
            }
            na++;
        }
    }
    Vec3d* epos  = (bRealloc) ? 0 : ff.epos;
    int*   espin = (bRealloc) ? 0 : ff.espin;
    ne = builder.exportElectronPairs( epos, espin, true, true );
    if(bRealloc){
        ff.realloc( na, ne, true );
    }
}

int builder2EFFstatic( EFF* ff, MM::Builder& builder, bool bCoreElectrons=true, bool bChangeCore=true, bool bChangeEsize=true, double esize0=0.5, double le=-0.5 ){
    //printf( "builder2EFFstatic() builder.atoms.size() %i builder.bonds.size() %i ff=%p \n", builder.atoms.size(), builder.bonds.size(), ff );
    //if(ff){ printf( "builder2EFFstatic() ff.na %i ff.ne %i \n", ff->na, ff->ne ); }
    int ie=0;
    int ia=0;
    Vec3d epos[4];
    for(int ib=0; ib<builder.bonds.size(); ib++){
        if( ff ){ 
            const MM::Bond& b = builder.bonds[ib];
            Vec3d pa = builder.atoms[b.atoms.i].pos;
            Vec3d pb = builder.atoms[b.atoms.j].pos;
            double c = 0.6;
            Vec3d p1 = pb*(1-c) + pa*c;
            Vec3d p2 = pa*(1-c) + pb*c;

            ff->epos[ie]=p1; ff->espin[ie]= 1; if(bChangeEsize){ff->esize[ie]=esize0;}; ie++;
            ff->epos[ie]=p2; ff->espin[ie]=-1; if(bChangeEsize){ff->esize[ie]=esize0;}; ie++;
            //ff->set_electron( ie, p, esize0,  1 );  ie++;
            //ff->set_electron( ie, p, esize0, -1 );  ie++;
        }else{ ie+=2; }
    }
    for(int ia=0; ia<builder.atoms.size(); ia++){
        //printf( "builder2EFFstatic() ia=%i na=%i \n", ia, builder.atoms.size() );
        Vec3d pa = builder.atoms[ia].pos;
        int   it = builder.atoms[ia].type;
        int   iZ = params.atypes[it].iZ;
        bool bCore = bCoreElectrons && (iZ>1);
        if( ff ){
            if( ia>=ff->na ){ printf( "builder2EFFstatic() ia(%i) >= ff->na(%i) => exit \n", ia, ff->na ); exit(0); }
            ff->apos [ia] = pa;
            if(bCore){
                ff->aPars[ia] = (Quat4d){ (double)iZ, 0.0, 0.0, 0.0 };
                //double ecsize = EFF::default_AtomParams[iZ].z;
                double ecsize = ff->atom_params[iZ].z;
                ff->epos[ie]=pa; ff->espin[ie]= 1; if(bChangeEsize){ff->esize[ie]=ecsize;}; ie++;
                ff->epos[ie]=pa; ff->espin[ie]=-1; if(bChangeEsize){ff->esize[ie]=ecsize;}; ie++;
                //ff->set_electron( ie, pa, ecsize,  1 ); ie++;
                //ff->set_electron( ie, pa, ecsize, -1 ); ie++;
            }else if (bChangeCore){
                Quat4d& apar = ff->aPars[ia];
                //apar = EFF::default_AtomParams[iZ];
                apar = ff->atom_params[iZ];
                if   (ff->bCoreCoul){ apar.x = iZ    ; }
                else                { apar.x = iZ-2.0; }
            }
        }else if(bCore){ ie+=2; }
        int nei = builder.addEpairsByPi(ia, le, epos, false );
        //printf( "builder2EFFstatic() ia=%i nei=%i \n", ia, nei );
        for(int i=0; i<nei; i++){
            if( ff ){
                ff->epos[ie]=epos[i]; ff->espin[ie]= 1; if(bChangeEsize){ff->esize[ie]=esize0;}; ie++;
                ff->epos[ie]=epos[i]; ff->espin[ie]=-1; if(bChangeEsize){ff->esize[ie]=esize0;}; ie++;
                //ff->set_electron( ie, epos[i], esize0,  1 ); ie++;
                //ff->set_electron( ie, epos[i], esize0, -1 ); ie++;
            }else{ ie+=2; }
        }
    }
    return ie;
}

int preAllocateXYZ(const char* fname, double Rfac=-0.5, bool bCoreElectrons=true){
    printf( "preAllocateXYZ() fname=%s Rfac=%g bCoreElectrons=%i\n", fname, Rfac, bCoreElectrons );
    if(params.atypes.size()==0){ 
        const char* sElementTypes ="common_resources/ElementTypes.dat"; 
        const char* sAtomTypes    ="common_resources/AtomTypes.dat"; 
        const char* sBondTypes    ="common_resources/BondTypes.dat"; 
        const char* sAngleTypes   ="common_resources/AngleTypes.dat"; 
        const char* sDihedralTypes="common_resources/DihedralTypes.dat"; 
        params.init(sElementTypes,sAtomTypes,sBondTypes,sAngleTypes,sDihedralTypes); 
    }
    FILE* fin=fopen(fname,"r"); if(fin==0){ printf("cannot open '%s' \n",fname); exit(0);} 
    const int nline=1024; 
    char line[1024], at_name[8]; 
    int il=0; 
    Atoms atoms; 
    builder.params=&params;
    while(fgets(line,nline,fin)){
        if(il==0){ 
            int na=-1; 
            sscanf(line,"%i",&na); 
            atoms.allocNew(na); 
            //printf( "preAllocateXYZ() na=%i atoms.apos=%p atoms.atypes=%p \n", na, atoms.apos, atoms.atypes ); 
        }
        else if( (il>1) && (il<atoms.natoms+2) ){ 
            double x,y,z,q; 
            int nret=sscanf(line,"%s %lf %lf %lf %lf",at_name,&x,&y,&z,&q); 
            if(nret<5) q=0; 
            int i=il-2; 
            //printf( "preAllocateXYZ() i=%i at_name=%s x=%g y=%g z=%g q=%g | %s \n", i, at_name, x, y, z, q, line ); 
            atoms.apos[i].set(x,y,z); 
            atoms.atypes[i]=params.getAtomType(at_name); 
        }
        if(il==atoms.natoms+2){ 
            builder.insertAtoms(atoms); 
            builder.tryAddConfsToAtoms(0,-1); 
            builder.autoBonds(Rfac); 
            int ne=builder2EFFstatic(0,builder,bCoreElectrons); 
            ff.realloc(atoms.natoms,ne,true); 
            opt.bindOrAlloc(ff.nDOFs,ff.pDOFs,ff.vDOFs,ff.fDOFs,ff.invMasses); 
            builder.load_atom_pos(atoms.apos,0); 
            builder2EFFstatic(&ff,builder,bCoreElectrons); 
            isPreAllocated=true; 
            break; 
        }
        il++;
    }
    fclose(fin);
    return 1;
}

int processXYZ( const char* fname, double Rfac=-0.5, double* outEs=0, double* apos_=0, double* epos_=0, int nstepMax=1000, double dt=0.001, double Fconv=1e-3, int ialg=2, bool bAddEpairs=false, bool bCoreElectrons=true, bool bChangeCore=true, bool bChangeEsize=true, const char* xyz_out="processXYZ.xyz", const char* fgo_out="processXYZ.fgo" ){
    setvbuf(stdout, NULL, _IONBF, 0);
    printf( "processXYZ(%s) bAddEpairs=%i Rfac %g xyz_out=%s fgo_out=%s \n", fname, bAddEpairs, Rfac, xyz_out, fgo_out );

    if(params.atypes.size()==0){
        const char* sElementTypes  = "common_resources/ElementTypes.dat";
        const char* sAtomTypes     = "common_resources/AtomTypes.dat"; 
        const char* sBondTypes     = "common_resources/BondTypes.dat"; 
        const char* sAngleTypes    = "common_resources/AngleTypes.dat";
        const char* sDihedralTypes = "common_resources/DihedralTypes.dat";
        params.init( sElementTypes, sAtomTypes, sBondTypes, sAngleTypes, sDihedralTypes );
        //exit(0);
    }
    
    FILE* fin = fopen( fname, "r" );
    if(fin==0){ printf("cannot open '%s' \n", fname ); exit(0);}
    const int nline=1024;
    char line[nline];
    char comment[nline];
    char at_name[8];
    // --- Open output file
    int il   = 0;
    Atoms*  atoms=0;
    builder.params = &params;
    bool bOnlyFirst = true;
    int iconf = 0; 

    while( fgets(line, nline, fin) ){
        if ( il==0 ){               // --- Read number of atoms
            int na=-1;
            sscanf( line, "%i", &na );
            if( (na>0)&&(na<10000) ){
                if( bOnlyFirst && (iconf==0) ){
                atoms = new Atoms(na);
            }
            }else{ printf( "ERROR in FitREQ::loadXYZ() Suspicious number of atoms (%i) while reading `%s`  => Exit() \n", na, fname ); exit(0); }
        }else if( il==1 ){               // --- Read comment line ( read reference energy )
            //printf("comment_line=%s\n",line);
            sscanf( line, "%*s %i %*s %lf ", &(atoms->n0), &(atoms->Energy) );
            sprintf(comment,"%s",line);
        }else if( il<atoms->natoms+2 ){  // --- Road atom line (type, position, charge)
            double x,y,z,q;
            int nret = sscanf( line, "%s %lf %lf %lf %lf", at_name, &x, &y, &z, &q );
            //printf( ".xyz[%i] %s %lf %lf %lf %lf\n", il, at_name, x, y, z, q );
            if(nret<5){q=0;}
            int i=il-2;
            atoms->apos[i].set(x,y,z);
            atoms->atypes[i]=params.getAtomType(at_name);
            //atoms->charge[i]=q;
        }
        if( il==atoms->natoms+2 ){
            if( (!isPreAllocated) && iconf==0 ){
                builder.insertAtoms( *atoms );
                builder.tryAddConfsToAtoms( 0, -1 );
                builder.printAtomConfs();
                builder.autoBonds( Rfac ); 
                int ne = builder2EFFstatic( 0, builder, bCoreElectrons, bChangeCore, bChangeEsize );
                if(verbosity>0)printf("processXYZ() iconf=%i natoms=%i ne=%i builder.atoms.size()=%i builder.bonds.size()=%i\n", iconf, atoms->natoms, ne, builder.atoms.size(), builder.bonds.size() );
                ff.realloc( atoms->natoms, ne, true );
                opt.bindOrAlloc( ff.nDOFs, ff.pDOFs, ff.vDOFs, ff.fDOFs, ff.invMasses );
            }
            builder.load_atom_pos( atoms->apos, 0 );
            builder2EFFstatic( &ff, builder, bCoreElectrons, bChangeCore, bChangeEsize );    

            if(verbosity>1)ff.info();
            if( nstepMax>0 ){
                { // constrain
                    int nfix=ff.na;
                    if(iconf==0)ff.realloc_fixed(nfix);
                    for(int i=0; i<nfix; i++){
                        ff.fixed_poss[i].f = ff.apos[i];
                        //ff.fixed_poss[i].w = 1;
                        ff.fixed_inds[i] = Vec2i{i,7};
                    }
                    //set_constrains( ff.na, fix, fixed_inds, true  );
                }
                //printf("processXYZ() iconf=%i natoms=%i na=%i ne=%i \n", iconf, atoms->natoms, ff.na, ff.ne );
                //builder2EFF( ff, builder );
                //run( int nstepMax, double dt, double Fconv, int ialg, double* outE, double* outF );
                run( nstepMax, dt, Fconv, ialg, 0, 0 );
            }
            ff.eval();
            //snprintf(comment, nline, "na,ne %i %i Etot(%.3f)=T(%.2f)+ee(%.3f)+ea(%.3f)+aa(%.3f)", ff.na, ff.ne, ff.Etot, ff.Ek, ff.Eee, ff.Eae, ff.Eaa);
            if(verbosity>0){
                //printf("processXYZ() iconf=%i natoms=%i | %s \n", iconf, atoms->natoms, comment);
                printf("processXYZ() iconf=%i na: %i ne: %i Etot(%.3f)=T(%.2f)+ee(%.3f)+ea(%.3f)+aa(%.3f)\n", iconf, ff.na, ff.ne, ff.Etot, ff.Ek, ff.Eee, ff.Eae, ff.Eaa );
            }
            if(xyz_out)ff.save_xyz(xyz_out, "a", comment);
            if(fgo_out)ff.writeTo_fgo(fgo_out, false, "a", iconf);
            ff.copyEnergies         (outEs, iconf);
            ff.copyAtomPositions    ((Vec3d*)apos_, iconf);
            ff.copyElectronPositions((Quat4d*)epos_, iconf);
            il=0;
            iconf++;
        }
        il++;
    }
    fclose(fin);
    return iconf;
}



// This functions takes .xyz file with full electron description
int processXYZ_e( const char* fname, double* outEs=0, double* apos_=0, double* epos_=0, int nstepMax=1000, double dt=0.001, double Fconv=1e-3, int optAlg=2, const char* xyz_out="processXYZ.xyz", const char* fgo_out="processXYZ.fgo" ){
    printf("processXYZ_e(%s)\n", fname);
    //double8* apars = ff.atom_params2;
    FILE* fin = fopen(fname, "r");
    if(!fin){ printf("Cannot open '%s'\n", fname); exit(0); }
    char line[1024];
    int il = 0, iconf = 0;
    int ie=0,ia=0;
    int nae=0,na=0,ne=0;
    while(fgets(line, sizeof(line), fin)){
        if      (il == 0){
            sscanf(line, "%d", &nae);
        }else if(il == 1){
            char coreMode;
            sscanf(line, "na,ne,core %i %i %c ", &na, &ne, &coreMode );
            ff.setCoreMode(coreMode);
            //na=nae-ne;
            if(nae!=na+ne){ printf("ERROR in processXYZ_e() nae(%i) != na(%i) + ne(%i) while reading `%s`  => Exit() \n", nae, na, ne, fname ); exit(0); }
            if(iconf == 0){
                ff.realloc(na, ne, true);
                //opt.bindOrAlloc(ff.nDOFs, ff.pDOFs, ff.vDOFs, ff.fDOFs, ff.invMasses);
                initOpt( dt, 0.1, 100.0, false );
            }
        }else{
            //printf( "particle_line[%i]: %s ", il, line );
            char type = ff.from_xyz_line(line, ie, ia);
        }
        il++;
        if (il>=nae+2){
            printf( "----- conf: %i \n", iconf );
            ff.info();
            if(nstepMax > 0) run(nstepMax, dt, Fconv, optAlg, 0, 0);   
            ff.eval();
            ff.copyEnergies         (         outEs, iconf );
            ff.copyAtomPositions    ((Vec3d* )apos_, iconf );
            ff.copyElectronPositions((Quat4d*)epos_, iconf );
            if(verbosity>0) printf(" processXYZ_e() iconf: %3i na: %3i ne: %3i Etot: %16.8f\n", iconf, na, ne, ff.Etot);
            if(xyz_out) ff.save_xyz(xyz_out, "a", line);
            il = 0;
            ia = 0;
            ie = 0;
            iconf++;
        }        
    }
    fclose(fin);
    return iconf;
}




// void eval_xyz_movie(){
//     Atoms::atomsFromXYZ();
// }


#endif // WITH_BUILDER




} // extern "C"
