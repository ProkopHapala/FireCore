
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
MMFFparams* params=0;
std::vector<Atoms*> samples;
MM::Builder builder;
#endif

// ============ Global Variables

EFF         ff;
DynamicOpt opt;

bool opt_initialized=false;

char* trj_fname = 0;
int   savePerNsteps = 1;


extern "C"{
// ========= Grid initialization

void setVerbosity( int verbosity_, int idebug_ ){
    //setbuf(stdout, NULL);  // Disable buffering completely
    setvbuf(stdout, NULL, _IOLBF, 0);  // Line buffering (flushed on newline)
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
    bool bConv=false;
    if(ff.nfix>0){ ff.apply_hard_fix(); }
    for(itr=0; itr<nstepMax; itr++ ){
        ff.clearForce();
        Etot = ff.eval();
        if( ff.bNegativeSizes & (verbosity>0) ){ printf( "negative electron sizes in step #%i => perhaps decrease relaxation time step dt=%g[fs]? \n", itr, opt.dt ); }

        if(ff.nfix>0){ if(ialg>0){ ff.clear_fixed_dynamics(); }else{ ff.clear_fixed_force(); } }
        switch(ialg){
            case  0: F2 = ff .move_GD      (dt);      break;
            case -1: F2 = opt.move_LeapFrog(dt);      break;
            case  1: F2 = opt.move_MD (dt,opt.damping); break;
            case  2: F2 = opt.move_FIRE();       break;
        }
        if(ff.nfix>0){ ff.apply_hard_fix(); }
        if(outE){ outE[itr]=Etot;     }
        if(outF){ outF[itr]=sqrt(F2); }
        if(verbosity>1){ printf("itr: %6i Etot[eV] %16.8f |F|[eV/A] %16.8f \n", itr, Etot, sqrt(F2) ); };
        if(F2<F2conv){
            //if(verbosity>0){ printf("Converged in %i iteration Etot %g[eV] |F| %g[eV/A] <(Fconv=%g) \n", itr, Etot, sqrt(F2), Fconv ); };
            bConv=true;
            break;
        }
        if( (trj_fname) && (itr%savePerNsteps==0) )  ff.save_xyz( trj_fname, "a" );
    }
    if(verbosity>0){ printf("%s in %6i iterations Etot[eV] %16.8f |F|[eV/A] %16.8f (Fconv=%g) \n", bConv ? "    CONVERGED" : "NOT-CONVERGED", itr, Etot, sqrt(F2), Fconv ); };
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
}

void relaxed_scan( int nconf, int nfix, double* fixed_poss, int* fixed_inds_, double* outEs, double* apos_, double* epos_, int nstepMax, double dt, double Fconv, int ialg, char* scan_trj_name ){
    //printf( "relaxed_scan() nconf %i nfix %i @fixed_poss=%p @fixed_inds=%p @outEs=%p @apos_=%p @epos_=%p @nstepMax %i @dt %g @Fconv %g @ialg %i \n", nconf, nfix, fixed_poss, fixed_inds_, outEs, apos_, epos_, nstepMax, dt, Fconv, ialg );
    //printf( "relaxed_scan() ialg: %i @vsize=%p @evel=%p @avel=%p \n", ialg, ff.vsize, ff.evel, ff.avel );
    Vec2i* fixed_inds = (Vec2i*)fixed_inds_;
    //for(int i=0; i<nfix; i++){ printf( "fixed[%3i] %3i %3i \n", i, fixed_inds[i].x, fixed_inds[i].y ); }
    //exit(0);   
    ff.realloc_fixed(nfix);
    for(int iconf=0; iconf<nconf; iconf++){
        if(verbosity>0)printf( "relaxed_scan() iconf %3i ", iconf );
        set_constrains( nfix, ((Quat4d*)fixed_poss)+iconf*nfix, fixed_inds, false );
        run( nstepMax, dt, Fconv, ialg, 0, 0 );
        outEs[iconf] = ff.Etot;
        // double Etot=0,Ek=0, Eee=0,EeePaul=0,EeeExch=0,  Eae=0,EaePaul=0,  Eaa=0; 
        {
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
        Vec3d*  apos = ((Vec3d* )apos_)+iconf*ff.na;
        Quat4d* epos = ((Quat4d*)epos_)+iconf*ff.ne;
        for(int j=0; j<ff.na; j++){ apos[j]   = ff.apos[j]; }
        for(int j=0; j<ff.ne; j++){ epos[j].f = ff.epos[j]; epos[j].w = ff.esize[j]; }
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


// #define WITH_BUILDER 1
#if WITH_BUILDER


void builder2EFF(EFF& ff, const MM::Builder& builder, bool bRealloc=true ){
    int ne=0;
    int na=0;
    int epairtyp = params->getAtomType("E");
    // --- loop over atoms and electron pairs
    for(int i=0; i<builder.atoms.size(); i++){
        int it = builder.atoms[i].type;
        if( it!=epairtyp ){
            int iZ = params->atypes[it].iZ;
            if(!bRealloc){
                ff.aPars [na] = EFF::default_AtomParams[iZ];
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

// void from_Atoms(EFF& ff, Atoms* atoms, bool bRealloc=true ){
//     int ne=0;
//     int na=0;
//     int epairtyp = params->getAtomType("E");
//     // --- loop over atoms and electron pairs
//     for(int i=0; i<atoms->natoms; i++){
//         int it = atoms->atypes[i];
//         if( it==epairtyp ){
//             ff.epos[ne] = atoms->apos[i];
//             ff.esize[ne] = 0.5;
//             ff.espin[ne] = 1;
//             ne++;
//             ff.epos[ne] = atoms->apos[i];
//             ff.esize[ne] = 0.5;
//             ff.espin[ne] = 1;
//             ne++;
//         }else{
//             int iZ = params->atypes[it].iZ;
//             ff.aPars [na] = EFF::default_AtomParams[iZ];
//             ff.apos  [na] = atoms->apos[i];
//             na++;
//         }
//     }
// }

int loadXYZ( const char* fname, bool bAddEpairs=false, bool bOutXYZ=false ){
    printf( "FitREQ::loadXYZ(%s) bAddEpairs=%i bOutXYZ=%i\n", fname, bAddEpairs, bOutXYZ );
    FILE* fin = fopen( fname, "r" );
    if(fin==0){ printf("cannot open '%s' \n", fname ); exit(0);}
    const int nline=1024;
    char line[1024];
    char at_name[8];
    // --- Open output file
    FILE* fout=0;
    if(bAddEpairs && bOutXYZ){
        sprintf(line,"%s_Epairs.xyz", fname );
        fout = fopen(line,"w");
    }
    int il   = 0;
    Atoms* atoms=0;
    while( fgets(line, nline, fin) ){
        if      ( il==0 ){               // --- Read number of atoms
            int na=-1;
            sscanf( line, "%i", &na );
            if( (na>0)&&(na<10000) ){
                atoms = new Atoms(na);
            }else{ printf( "ERROR in FitREQ::loadXYZ() Suspicious number of atoms (%i) while reading `%s`  => Exit() \n", na, fname ); exit(0); }
        }else if( il==1 ){               // --- Read comment line ( read reference energy )
            sscanf( line, "%*s %*s %i %*s %lf ", &(atoms->n0), &(atoms->Energy) );
        }else if( il<atoms->natoms+2 ){  // --- Road atom line (type, position, charge)
            double x,y,z,q;
            int nret = sscanf( line, "%s %lf %lf %lf %lf", at_name, &x, &y, &z, &q );
            //printf( ".xyz[%i] %s %lf %lf %lf %lf\n", il, at_name, x, y, z, q );
            if(nret<5){q=0;}
            int i=il-2;
            atoms->apos[i].set(x,y,z);
            atoms->atypes[i]=params->getAtomType(at_name);
            //atoms->charge[i]=q;
        }
        il++;
        if( bAddEpairs){    // ---- store sample atoms to batch
            builder.params = params;
            //atoms = builder.buildBondsAndEpairs(atoms);
            builder2EFF( ff, builder );

            //builder.exportAtoms( Vec3d* apos, int* atypes, int i0=0, int n=-1 );
            //builder.exportElectronPairs( ff.epos, ff.espin, true, true );

            // int nbad = atoms->checkTypeInRange( params->atypes.size()-1, 1, true );
            // if( nbad>0 ){ 
            //     printf("ERROR in FitREQ::loadXYZ() samples[%5i]->checkTypeInRange() return nbad=%i => exit()\n", nbatch, nbad ); 
            //     params->printTypesOfAtoms(atoms->natoms, atoms->atypes);
            //     exit(0); 
            // }
            //samples.push_back( atoms );
            //if(bEvalOnlyCorrections){ printFittedAdata( samples.size()-1 ); }
            il=0;
        }
    }
    if(fout)fclose(fout);
    fclose(fin);
    return samples.size();
}




// void eval_xyz_movie(){
//     Atoms::atomsFromXYZ();
// }


#endif // WITH_BUILDER




} // extern "C"
