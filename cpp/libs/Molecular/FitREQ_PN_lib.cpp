#include "globals.h"
#include "testUtils.h"
#include "FitREQ_PN.h"

MMFFparams params;
FitREQ_PN W;

#include "OptRandomWalk.h"
OptRandomWalk ropt;

//============================

#include "libUtils.h"

extern "C"{


// =========================
// Sampling interface for testing C++ dampers from Python
// kind codes:
//  0  : bare coulomb 1/r
// 10-14 : Boys family (exact, C1 cubic, C2 quintic, C1 even quartic, C2 even sextic) with params: [rmin]
// 20-23 : Soft-clamp family (soft+, smooth+, soft-, smooth-) with params: [y1, y2]
// Output EFs_ is an array of Vec2d (E,F) where F = - dE/dr (radial force)
void setVerbosity( int verbosity_, int idebug_, int PrintDOFs, int PrintfDOFs, int PrintBeforReg, int PrintAfterReg, int PrintOverRepulsive ){
    verbosity = verbosity_;
    idebug    = idebug_;
    // no buffering
    setbuf(stdout, NULL);
    setbuf(stderr, NULL);
    #define _setbool(name) { if(name>0){W.b##name=true;}else if(name<0){W.b##name=false;} }
    _setbool( PrintDOFs     );
    _setbool( PrintfDOFs    );
    _setbool( PrintBeforReg );
    _setbool( PrintAfterReg );
    _setbool( PrintOverRepulsive );
    #undef _setbool
    printf( "setVerbosity() verbosity %i idebug %i PrintDOFs %i PrintfDOFs %i PrintBeforReg %i PrintAfterReg %i PrintOverRepulsive %i\n", verbosity, idebug, W.bPrintDOFs, W.bPrintfDOFs, W.bPrintBeforReg, W.bPrintAfterReg, W.bPrintOverRepulsive );
}

void setModel( int ivdW, int iCoul, int iHbond, int Epairs, int iEpairs, double kMorse, double Lepairs, bool bPN ){
    W.ivdW    = ivdW;
    W.iCoul   = iCoul;
    W.iHbond  = iHbond;
    if(Epairs>0){W.bEpairs=true;}else if(Epairs<0){W.bEpairs=false;}
    W.iEpairs = iEpairs;
    W.kMorse  = kMorse;
    W.Lepairs = Lepairs;
    W.bPN     = bPN;
    printf( "setModel() ivdW %i iCoul %i iHbond %i iEpairs %i kMorse %f Lepairs %f bPN %i \n", W.ivdW, W.iCoul, W.iHbond, W.iEpairs, W.kMorse, W.Lepairs, W.bPN );
}

void loadTypes( const char* fname_ElemTypes, const char* fname_AtomTypes ){
    params.loadElementTypes( fname_ElemTypes );
    params.loadAtomTypes( fname_AtomTypes ); 
    W.params=&params;
}

int loadDOFSelection( const char* fname ){
    return W.loadDOFSelection( fname );
}

int loadXYZ( const char* fname, bool bAddEpairs, bool bOutXYZ, bool bSaveJustElementXYZ, char* OutXYZ_fname, bool bEvalOnlyCorrections, bool bAppend ){ 
    W.bSaveJustElementXYZ = bSaveJustElementXYZ;
    W.bEvalOnlyCorrections=bEvalOnlyCorrections; 
    //printf( "loadXYZ(fname=%s, bAddEpairs=%i, bOutXYZ=%i, bSaveJustElementXYZ=%i, bEvalOnlyCorrections=%i, bAppend=%i ) \n", fname, bAddEpairs, bOutXYZ, bSaveJustElementXYZ, bEvalOnlyCorrections, bAppend );
    return W.loadXYZ( fname, bAddEpairs, bOutXYZ, OutXYZ_fname, bAppend );
}

void init_buffers(){
    //printf( "init_buffers() \n" );

    ibuffers.insert( { "ndims", &W.nDOFs  } );
    ibuffers.insert( { "typToREQ",      (int*)W.typToREQ  } );

    buffers.insert( { "DOFs",  (double*)W.DOFs  } );
    buffers.insert( { "fDOFs", (double*)W.fDOFs } );
    buffers.insert( { "vDOFs", (double*)W.fDOFs } );

    buffers.insert( { "fDOFbounds",     (double*)W.fDOFbounds } );
    
    buffers.insert( { "typeREQs",       (double*)W.typeREQs  } );
    buffers.insert( { "typeREQsMin",    (double*)W.typeREQsMin  } );
    buffers.insert( { "typeREQsMax",    (double*)W.typeREQsMax  } );

    buffers.insert( { "typeREQs0",      (double*)W.typeREQs0 } );
    buffers.insert( { "typeREQs0_low",  (double*)W.typeREQs0_low } );
    buffers.insert( { "typeREQs0_high", (double*)W.typeREQs0_high } );

    buffers.insert( { "typeKreg",       (double*)W.typeKreg  } );
    buffers.insert( { "typeKreg_low",   (double*)W.typeKreg_low  } );
    buffers.insert( { "typeKreg_high",  (double*)W.typeKreg_high  } );

    if(W.weights)buffers.insert( { "weights", (double*)W.weights  } );

    //ibuffers.insert( { "vDOFs", (double*)W.fDOFs } );
    //printBuffNames();
}

void setPenalty( int Clamp, int Regularize, int AddRegError, int RegCountWeight, int SoftClamp, double softClamp_start, double softClamp_max ){
    //if(Clamp>0){W.bClamp=true;}else if(Clamp<0){W.bClamp=false;}
    #define _setbool(name) { if(name>0){W.b##name=true;}else if(name<0){W.b##name=false;} }
    _setbool( Clamp );
    _setbool( Regularize     );
    _setbool( AddRegError    );
    _setbool( RegCountWeight );
    _setbool( SoftClamp );
    #undef _setbool
    W.softClamp_start = softClamp_start;
    W.softClamp_max   = softClamp_max;
    printf( "setPenalty() Clamp %i Regularize %i AddRegError %i RegCountWeight %i SoftClamp %i softClamp_start=%f softClamp_max=%f\n", W.bClamp, W.bRegularize, W.bAddRegError, W.bRegCountWeight, W.bSoftClamp, W.softClamp_start, W.softClamp_max );
}

void setWeights( int n, double* weights ){
   _realloc( W.weights, n );
   for(int i=0; i<n; i++){ W.weights[i]=weights[i]; }     
}

void setTrjBuffs( double* trj_E, double* trj_F, double* trj_DOFs, double* trj_fDOFs){
    W.trj_E = trj_E;
    W.trj_F = trj_F;
    W.trj_DOFs = trj_DOFs;
    W.trj_fDOFs = trj_fDOFs;
}

double run_PN( int ialg, int iparallel, int nstep, double Fmax, double dt, double max_step, double damping ){
    printf( "run(ialg=%i,iparallel=%i,imodel=%i,nstep=%6i,nsamp=%6i)\n", ialg, iparallel, W.imodel, nstep, W.samples.size() );
    long t0 = getCPUticks();
    double Err=0;
    switch (iparallel){
        case 0:{ Err=W.run_PN ( ialg, nstep, Fmax, dt, max_step, damping, false ); } break;
        case 1:{ Err=W.run_PN ( ialg, nstep, Fmax, dt, max_step, damping, true  ); } break;
    }
    double T = (getCPUticks()-t0);
    printf( "Time: run(nstep=%6i,nsamp=%6i,iparallel=%i) T= %8.3f [MTicks] %8.3f [ticks/conf]\n", nstep, W.samples.size(), iparallel, T*1e-6, T/(W.samples.size()*nstep) );
    return Err;
}

double getEs( double* Es, double* Fs, bool bOmp, bool bDOFtoTypes, char* xyz_name){
    printf( "getEs() imodel %i nDOFs %i bOmp %i bDOFtoTypes %i bEvalOnlyCorrections=%i xyz_name=%s @DOFtoTyp=%p  @Es=%p @Fs=%p \n", W.imodel, W.nDOFs, bOmp, bDOFtoTypes, W.bEvalOnlyCorrections, xyz_name, W.DOFtoTyp, Es, Fs );
    //W.imodel=imodel;
    if(xyz_name){ W.xyz_out=xyz_name; W.bSaveSampleToXYZ=true; }else{ W.bSaveSampleToXYZ=false; }
    if(bDOFtoTypes)W.DOFsToTypes(); 
    W.clean_fDOFs();
    double E = 0;
    if( bOmp ){ E = W.evalSamples_omp  ( Es ); }
    else      { E = W.evalSamples_serial( Es ); }
    //for(int i=0; i<W.samples.size(); i++){ printf( "getEs() sample[%i] E: %20.10f \n", i, Es[i] ); }
    if( Fs ){ for(int i=0; i<W.nDOFs; i++){ Fs[i] = W.fDOFs[i]; } }
    W.bSaveSampleToXYZ=false; 
    return E;
}

void getEs_components( double* Es, double* Es_Coul, double* Es_vdW, double* Es_Epairs, double* Es_Hbond ){
    W.DOFsToTypes(); 
    int nsamp = W.samples.size();
    for(int isamp=0; isamp<nsamp; isamp++){
        Atoms* atoms  = W.samples[isamp];
        const AddedData* adata = (const AddedData*)(atoms->userData);
        alignas(32) double Qs  [atoms->natoms];
        alignas(32) Vec3d  apos[atoms->natoms];   // atomic positions
        W.fillTempArrays( atoms, apos, Qs );
        int     nj = atoms->n0;
        int     j0 = 0; 
        int     ni = atoms->natoms - atoms->n0;
        int     i0 = atoms->n0;
        W.evalEnergyComponents( i0, ni, j0, nj, atoms->atypes, apos, W.typeREQs, Qs, adata->host, Es_Coul[isamp], Es_vdW[isamp], Es_Epairs[isamp], Es_Hbond[isamp] );
        Es[isamp] = Es_Coul[isamp] + Es_vdW[isamp] + Es_Epairs[isamp] + Es_Hbond[isamp];
        //printf( "getEs_components() sample[%i] E: %20.10f E_Coul: %20.10f E_vdW: %20.10f E_Epairs: %20.10f E_Hbond: %20.10f\n", isamp, Es[isamp], Es_Coul[isamp], Es_vdW[isamp], Es_Epairs[isamp], Es_Hbond[isamp] );
    }
}

double getError( int iparallel ){
    W.bSaveSampleToXYZ=false; 
    if(W.weights){W.updateWeightsSum();}
    double Err=0.0;
    W.clear_fDOFbounds();
    switch (iparallel){
        case 0:{ Err = W.evalFitError( 0, false ); } break;
        case 1:{ Err = W.evalFitError( 0, true  ); } break;
    }
    return Err;
}

void scanParam( int iDOF, int n, double* xs,  double* Es, double* Fs, bool bEvalSamples ){
    //W.bRegularize=bRegularize;
    bool bOmp = W.iparallel>0;
    printf( "scanParam() iDOF: %i imodel: %i n: %i bOmp: %i\n", iDOF, W.imodel, n, bOmp );
    W.clear_fDOFbounds();
    if(bOmp){ W.bBroadcastFDOFs=true; W.realloc_sample_fdofs();  }
    for(int i=0; i<n; i++){
        W.DOFs[iDOF] = xs[i];
        //printf( "\n##### scanParam()[%3i] W.DOFs[%3i]: %20.10f      %8s.%c \n", i, iDOF,     W.DOFs[iDOF],    W.params->atypes[W.DOFtoTyp[iDOF].x].name,  "REQH"[W.DOFtoTyp[iDOF].y]  );
        double E = W.evalFitError( i, bOmp, bEvalSamples );
        if(Fs)Es[i] = E;
        if(Fs)Fs[i] = W.fDOFs[iDOF];
        //if( verbosity>1){ printf( "scanParam()[%3i] W.DOFs[%3i]: %20.10f E: %20.10f  F: %20.10f \n", i, iDOF, W.DOFs[iDOF], E, W.fDOFs[iDOF] ); }
    }
}

} // extern "C"
