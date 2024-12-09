

#include "globals.h"

#include "testUtils.h"
#include "FitREQ.h"

MMFFparams params;
FitREQ W;

#include "OptRandomWalk.h"
OptRandomWalk ropt;

//============================

//const std::vector<std::string>* atomTypeNames = 0;
#include "libUtils.h"

extern "C"{

void setVerbosity( int verbosity_, int idebug_ ){
    verbosity = verbosity_;
    idebug    = idebug_;
}

void setSwitches( int EvalJ, int WriteJ, int CheckRepulsion, int Regularize, int Epairs){
    #define _setbool(name) { if(name>0){W.b##name=true;}else if(name<0){W.b##name=false;} }
    _setbool( EvalJ          );
    _setbool( WriteJ         );
    _setbool( CheckRepulsion );
    _setbool( Regularize     );
    _setbool( Epairs         );
    #undef _setbool
}

int loadTypeSelection( const char* fname ){
    return W.loadTypeSelection( fname );
}

int loadWeights( const char* fname ){
    return W.loadWeights( fname );
}

void loadTypes( const char* fname_ElemTypes, const char* fname_AtomTypes ){
    params.loadElementTypes( fname_ElemTypes );
    params.loadAtomTypes( fname_AtomTypes ); 
    W.params=&params;
}

int loadXYZ( const char* fname, bool bAddEpairs, bool bOutXYZ ){  return W.loadXYZ( fname, bAddEpairs, bOutXYZ );}

void setWeights( int n, double* weights ){
   _realloc( W.weights, n );
   for(int i=0; i<n; i++){ W.weights[i]=weights[i]; }     
}

int export_Erefs( double* Erefs ){ return W.export_Erefs( Erefs ); }

double run( int nstep, double Fmax, double dt, int imodel_, int ialg, int iparallel, bool bClamp, double max_step ){
    printf( "run(nstep=%6i,nsamp=%6i,iparallel=%i) bEvalJ=%i bWriteJ=%i bJ=%i \n", nstep, W.samples.size(), iparallel, W.bEvalJ, W.bWriteJ, (W.bEvalJ&&(!W.bWriteJ)) );
    long t0 = getCPUticks();
    double Err=0;
    switch (iparallel){
        case 0:{ Err=W.run    ( nstep, Fmax, dt, imodel_, ialg, false, bClamp, max_step ); } break;
        case 1:{ Err=W.run    ( nstep, Fmax, dt, imodel_, ialg, true,  bClamp, max_step ); } break;
        case 2:{ Err=W.run_omp( nstep, Fmax, dt, imodel_, ialg,        bClamp, max_step ); } break;
    }
    double T = (getCPUticks()-t0);
    printf( "Time: run(nstep=%6i,nsamp=%6i,iparallel=%i) T= %8.3f [MTicks] %8.3f [ticks/conf]\n", nstep, W.samples.size(), iparallel, T*1e-6, T/(W.samples.size()*nstep) );
    return Err;
}

static double evalFitError( int n, double * Xs ){
    for(int i=0; i<W.nDOFs; i++) { W.DOFs[i] = Xs[i]; }
    return W.evalFitError( -1, W.iparallel>0 );
};

double optimize_random(int nstep, double stepSize=0.1) {
    ropt.getEnergy = evalFitError;
    ropt.realloc(W.nDOFs, W.DOFs);
    ropt.run(nstep, stepSize);
    ropt.dealloc();
    return ropt.Ebest;
}

void setTypeToDOFs  (int i, double* REQ ){ W.setTypeToDOFs  ( i, *(Quat4d*)REQ ); }
void getTypeFromDOFs(int i, double* REQ ){ W.getTypeFromDOFs( i, *(Quat4d*)REQ ); }

double getEs( int imodel, double* Es, double* Fs, bool bOmp, bool bDOFtoTypes ){
    //printf( "getEs() imodel %i bOmp %i bDOFtoTypes %i \n", imodel, bOmp, bDOFtoTypes );
    W.imodel=imodel;
    if(bDOFtoTypes)W.DOFsToTypes(); 
    W.clean_fDOFs();
    double E = 0;
    if( bOmp ){ E = W.evalSamples_omp  ( Es ); }
    else      { E = W.evalSamples_noOmp( Es ); }
    if( Fs ){ for(int i=0; i<W.nDOFs; i++){ Fs[i] = W.fDOFs[i]; } }
    return E;
}

void scanParam( int iDOF, int imodel,  int n, double* xs,  double* Es, double* Fs, bool bRegularize ){
    //printf( "scanParam() iDOF %i imodel %i n %i \n", iDOF, imodel, n );
    W.imodel=imodel;
    for(int i=0; i<n; i++){
        W.DOFs[iDOF] = xs[i];
        double E = W.evalFitError( i, W.iparallel>0 );
        if(Fs)Es[i] = E;
        if(Fs)Fs[i] = W.fDOFs[iDOF];
        //if( verbosity>1){ printf( "scanParam()[%3i] W.DOFs[%3i]: %20.10f E: %20.10f  F: %20.10f \n", i, iDOF, W.DOFs[iDOF], E, W.fDOFs[iDOF] ); }
    }
}

void scanParam2D( int iDOFx, int iDOFy, int imodel, int nx, int ny, double* xs, double* ys,  double* Es, double* Fx, double* Fy, bool bRegularize ){
    //printf( "scanParam() iDOF %i imodel %i nx %i ny %i \n", iDOFx, imodel, nx, ny );
    W.imodel=imodel;
    for(int iy=0; iy<ny; iy++){
        W.DOFs[iDOFy] = ys[iy];
        for(int ix=0; ix<nx; ix++){
            int i = iy*nx + ix; 
            //printf( "scanParam2D() ix %i iy %i i %i \n", ix, iy, i );
            W.DOFs[iDOFx] = xs[ix];
            double E = W.evalFitError( i, W.iparallel>0 );
            if(Es)Es[i] = E;
            if(Fx)Fx[i] = W.fDOFs[iDOFx];
            if(Fy)Fy[i] = W.fDOFs[iDOFy];
        }
    }
}

void init_buffers(){
    //printf( "init_buffers() \n" );

    ibuffers.insert( { "ndims", &W.nDOFs  } );
    ibuffers.insert( { "typToREQ",      (int*)W.typToREQ  } );

    buffers.insert( { "DOFs",  (double*)W.DOFs  } );
    buffers.insert( { "fDOFs", (double*)W.fDOFs } );
    buffers.insert( { "vDOFs", (double*)W.fDOFs } );

    buffers.insert( { "typeREQs",       (double*)W.typeREQs  } );
    buffers.insert( { "typeREQsMin",    (double*)W.typeREQsMin  } );
    buffers.insert( { "typeREQsMax",    (double*)W.typeREQsMax  } );

    buffers.insert( { "typeREQs0",      (double*)W.typeREQs0 } );
    buffers.insert( { "typeREQs0_low",  (double*)W.typeREQs0_low } );
    buffers.insert( { "typeREQs0_high", (double*)W.typeREQs0_high } );

    buffers.insert( { "typeKreg",       (double*)W.typeKreg  } );
    buffers.insert( { "typeKreg_low",   (double*)W.typeKreg_low  } );
    buffers.insert( { "typeKreg_high",  (double*)W.typeKreg_high  } );

    //ibuffers.insert( { "vDOFs", (double*)W.fDOFs } );
    //printBuffNames();
}

} // extern "C"
