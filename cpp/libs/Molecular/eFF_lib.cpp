
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
int idebug   =0;

const char* prefix = "#Epiece";

#include "InteractionsGauss.h"
#include "eFF.h"
#include "DynamicOpt.h"
#include "Vec3Utils.h"

// ============ Global Variables

EFF         ff;
DynamicOpt opt;

bool opt_initialized=false;

char* trj_fname = 0;
int   savePerNsteps = 1;


extern "C"{
// ========= Grid initialization

void setVerbosity( int verbosity_, int idebug_ ){
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
    ff.realloc ( na, ne );
    //ff.autoAbWs( default_aAbWs, default_eAbWs );
    init_buffers();
}

void info(){ ff.info(); }

double* getEnergyPointer(){ return &ff.Ek; }
int*    getDimPointer   (){ return &ff.ne; }

void initOpt( double dt, double damping, double f_limit, bool bMass ){
    if(ff.vDOFs){ if(bMass){ff.makeMasses(ff.invMasses);}else{ff.makeMasses(ff.invMasses,1.0);} }
    opt.bindOrAlloc( ff.nDOFs, ff.pDOFs, ff.vDOFs, ff.fDOFs, ff.invMasses );
    //opt.cleanVel( ); // this is already inside initOpt
    //opt.initOpt( dt, damping );
    opt.setTimeSteps(dt );
    opt.setDamping  (damping);
    opt.f_limit = f_limit;
    opt_initialized=true;
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
    if(ialg>0){ opt.cleanVel( ); }
    for(itr=0; itr<nstepMax; itr++ ){
        ff.clearForce();
        Etot = ff.eval();
        if( ff.bNegativeSizes & (verbosity>0) ){ printf( "negative electron sizes in step #%i => perhaps decrease relaxation time step dt=%g[fs]? \n", itr, opt.dt ); }
        switch(ialg){
            case  0: ff .move_GD      (dt);      break;
            case -1: opt.move_LeapFrog(dt);      break;
            case  1: F2 = opt.move_MD (dt,opt.damping); break;
            case  2: F2 = opt.move_FIRE();       break;
        }
        if(outE){ outE[itr]=Etot; }
        if(outF){ outF[itr]=F2;   }
        if(verbosity>0){ printf("[%i] Etot %g[eV] |F| %g [eV/A] \n", itr, Etot, sqrt(F2) ); };
        if(F2<F2conv){
            if(verbosity>0){ printf("Converged in %i iteration Etot %g[eV] |F| %g[eV/A] <(Fconv=%g) \n", itr, Etot, sqrt(F2), Fconv ); };
            break;
        }
        if( (trj_fname) && (itr%savePerNsteps==0) )  ff.save_xyz( trj_fname, "a" );
    }
    //printShortestBondLengths();
    return itr;
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

void setSwitches( int bEvalKinetic, int bEvalCoulomb, int  bEvalPauli, int bEvalAA, int bEvalAE, int bEvalAECoulomb, int bEvalAEPauli ){
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
#undef _setbool
}

} // extern "C"
