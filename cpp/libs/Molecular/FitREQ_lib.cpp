

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


// =========================
// Sampling interface for testing C++ dampers from Python
// kind codes:
//  0  : bare coulomb 1/r
// 10-14 : Boys family (exact, C1 cubic, C2 quintic, C1 even quartic, C2 even sextic) with params: [rmin]
// 20-23 : Soft-clamp family (soft+, smooth+, soft-, smooth-) with params: [y1, y2]
// Output EFs_ is an array of Vec2d (E,F) where F = - dE/dr (radial force)
void sample_funcEF( int n, double* xs, double* EFs_, int kind, double* params ){
    Vec2d* EFs = (Vec2d*)EFs_;
    for(int i=0; i<n; i++ ){
        double x = xs[i];
        double y=0.0, dydr=0.0;
        switch (kind){
            case 0:{ // bare
                if(x<1e-16){ y=1e16; dydr=-1e32; }
                else{ y = 1.0/x; dydr = -1.0/(x*x); }
            }break;
            // Boys family
            case 10: case 11: case 12: case 13: case 14:{
                double rmin = params? params[0] : W.boys_rmin;
                int mode = 0;
                if(kind==11) mode=1; else if(kind==12) mode=2; else if(kind==13) mode=3; else if(kind==14) mode=4; // 10=exact
                dampCoulomb_Boys_val_deriv( x, rmin, mode, y, dydr );
            }break;
            // Soft clamp family (positive/negative)
            case 20: case 21: case 22: case 23:{
                double y1 = params? params[0] : W.clamp_y1;
                double y2 = params? params[1] : W.clamp_y2;
                int mode;
                if(kind==20) mode=1; else if(kind==21) mode=2; else if(kind==22) mode=3; else mode=4;
                dampCoulomb_SoftClamp_val_deriv( x, y1, y2, mode, y, dydr );
            }break;
            default:{ y=NAN; dydr=NAN; }break;
        }
        EFs[i].x = y;
        EFs[i].y = -dydr; // force is -dE/dr = -dy/dr (since charge scaling can be applied in Python if needed)
    }
}

void setVerbosity( int verbosity_, int idebug_, int PrintDOFs, int PrintfDOFs, int PrintBeforReg, int PrintAfterReg ){
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
    #undef _setbool
    printf( "setVerbosity() verbosity %i idebug %i \n", verbosity, idebug );
    printf( "setVerbosity() PrintDOFs %i PrintfDOFs %i PrintBeforReg %i PrintAfterReg %i \n", W.bPrintDOFs, W.bPrintfDOFs, W.bPrintBeforReg, W.bPrintAfterReg );
}

// bool  bEvalJ          = false;    // Should we evaluate variational derivatives on Fregment J 
// bool  bWriteJ         = false;    // Should we write variational derivatives to Fregment J ( inner loop over j )
// bool  bCheckRepulsion = false;    // Should we check maximum repulsion (EijMax) inside inner loop over j for each sample atoms ?
// bool  bRegularize     = true;     // Should we apply additional regularization forces to otimizer ( beside the true variational forces from inter-atomic forcefield ? )
// bool  bAddRegError    = true;     // Should we add regularization error to total error ?
// bool  bEpairs         = true;     // Should we add electron pairs to the molecule ?
// //bool  bOptEpR = false;          // Should we optimize electron pair distance (from host atom) ?
// bool  bBroadcastFDOFs = false;    // Should we broadcast fDOFs (each sample to its own chunk of memory) to prevent atomic-write conflicts ?
// bool  bUdateDOFbounds = true;     // Should we update fDOFbounds after each sample ?
// bool  bEvalOnlyCorrections = false;  // Split evaluation and optimization to Emodel0 and Ecorrection (where only Ecorrection is updated every iteration)
void setup( int imodel, int EvalJ, int WriteJ, int CheckRepulsion, int Regularize, int RegCountWeight, int AddRegError, int Epairs, int BroadcastFDOFs, int UdateDOFbounds, int EvalOnlyCorrections, int SaveJustElementXYZ, int SoftClamp){
    printf("C++ setup()  SoftClamp %i \n", SoftClamp );
    W.imodel = imodel;
    #define _setbool(name) { if(name>0){W.b##name=true;}else if(name<0){W.b##name=false;} }
    _setbool( EvalJ          );
    _setbool( WriteJ         );
    _setbool( CheckRepulsion );
    _setbool( Regularize     );
    _setbool( RegCountWeight );
    _setbool( AddRegError    );
    _setbool( Epairs         );
    _setbool( BroadcastFDOFs );
    _setbool( UdateDOFbounds );
    _setbool( EvalOnlyCorrections );
    _setbool( SaveJustElementXYZ );
    _setbool( SoftClamp );
    printf( "setup() imodel %i EvalJ %i WriteJ %i CheckRepulsion %i Regularize %i RegCountWeight %i AddRegError %i Epairs %i BroadcastFDOFs %i UdateDOFbounds %i EvalOnlyCorrections %i SaveJustElementXYZ %i SoftClamp %i \n", 
                     W.imodel, W.bEvalJ, W.bWriteJ, W.bCheckRepulsion, W.bRegularize, W.bRegCountWeight, W.bAddRegError, W.bEpairs, W.bBroadcastFDOFs, W.bUdateDOFbounds, W.bEvalOnlyCorrections, W.bSaveJustElementXYZ, W.bSoftClamp );
    #undef _setbool
}

void setGlobalParams( double kMorse, double Lepairs, double EijMax, double softClamp_start, double softClamp_max ){
    W.kMorse   = kMorse;
    W.Lepairs  = Lepairs;
    W.EijMax   = EijMax;
    W.softClamp_start = softClamp_start;
    W.softClamp_max   = softClamp_max;
    printf( "setGlobalParams() kMorse %g Lepairs %g EijMax %g softClamp_start %g softClamp_max %g \n", kMorse, Lepairs, EijMax, softClamp_start, softClamp_max );
}

//bool bListOverRepulsive    = true;   // Should we list overrepulsive samples? 
//bool bSaveOverRepulsive    = false;  // Should we save overrepulsive samples to .xyz file?
//bool bPrintOverRepulsive   = true;   // Should we print overrepulsive samples? 
//bool bDiscardOverRepulsive = true;   // Should we discard overrepulsive samples? ( i.e. ignore them as training examples )

// int    iWeightModel    = 1;    // weight of model energy 1=linear, 2=cubic_smooth_step  
// double EmodelCut       = 10.0; // sample model energy when we consider it too repulsive and ignore it during fitting
// double EmodelCutStart  = 5.0;  // sample model energy when we start to decrease weight in the fitting  

void setFilter( double EmodelCut, double EmodelCutStart, int iWeightModel, int ListOverRepulsive, int SaveOverRepulsive, int PrintOverRepulsive, int DiscardOverRepulsive, int WeightByEmodel ){
    W.EmodelCut      = EmodelCut;
    W.EmodelCutStart = EmodelCutStart;
    W.iWeightModel   = iWeightModel;
    #define _setbool(name) { if(name>0){W.b##name=true;}else if(name<0){W.b##name=false;} }
    _setbool( ListOverRepulsive );
    _setbool( SaveOverRepulsive );
    _setbool( PrintOverRepulsive);
    _setbool( DiscardOverRepulsive );
    _setbool( WeightByEmodel );
    #undef _setbool
    printf( "setFilter(): EmodelCut=%g EmodelCutStart=%g iWeightModel=%i ListOverRepulsive=%i SaveOverRepulsive=%i PrintOverRepulsive=%i DiscardOverRepulsive=%i \n",  W.EmodelCut, W.EmodelCutStart, W.iWeightModel, W.bListOverRepulsive, W.bSaveOverRepulsive, W.bPrintOverRepulsive, W.bDiscardOverRepulsive );
    { printf( "setFilter() weight samples:\n #i     E            weight \n" ); int n=10; for(int i=0; i<=n; i++){  double wi,E=W.EmodelCutStart+i*(W.EmodelCut-W.EmodelCutStart)/n; W.smoothWeight( E, wi ); printf( "%3i %16.8e %16.8e \n", i, E, wi ); } }
}


int loadDOFSelection( const char* fname ){
    return W.loadDOFSelection( fname );
}

// int loadTypeSelection( const char* fname ){
//     return W.loadTypeSelection( fname );
// }

int loadWeights( const char* fname ){
    return W.loadWeights( fname );
}

void loadTypes( const char* fname_ElemTypes, const char* fname_AtomTypes ){
    params.loadElementTypes( fname_ElemTypes );
    params.loadAtomTypes( fname_AtomTypes ); 
    W.params=&params;
}

int loadXYZ( const char* fname, bool bAddEpairs, bool bOutXYZ, bool bEvalOnlyCorrections, bool bAppend ){ W.bEvalOnlyCorrections=bEvalOnlyCorrections; return W.loadXYZ( fname, bAddEpairs, bOutXYZ, bAppend );}

void setWeights( int n, double* weights ){
   _realloc( W.weights, n );
   for(int i=0; i<n; i++){ W.weights[i]=weights[i]; }     
}

int export_Erefs( double* Erefs ){ return W.export_Erefs( Erefs ); }


void setTrjBuffs( double* trj_E, double* trj_F, double* trj_DOFs, double* trj_fDOFs){
    W.trj_E = trj_E;
    W.trj_F = trj_F;
    W.trj_DOFs = trj_DOFs;
    W.trj_fDOFs = trj_fDOFs;
}

double run( int ialg, int iparallel, int nstep, double Fmax, double dt, double max_step, double damping, bool bClamp ){
    printf( "run(ialg=%i,iparallel=%i,imodel=%i,nstep=%6i,nsamp=%6i) bEvalJ=%i bWriteJ=%i bJ=%i \n", ialg, iparallel, W.imodel, nstep, W.samples.size(),  W.bEvalJ, W.bWriteJ, (W.bEvalJ&&(!W.bWriteJ)) );
    long t0 = getCPUticks();
    double Err=0;
    switch (iparallel){
        case 0:{ Err=W.run    ( ialg, nstep, Fmax, dt, max_step, damping, bClamp, false ); } break;
        case 1:{ Err=W.run    ( ialg, nstep, Fmax, dt, max_step, damping, bClamp, true  ); } break;
        case 2:{ Err=W.run_omp( ialg, nstep, Fmax, dt, max_step, damping, bClamp        ); } break;
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

double getEs( double* Es, double* Fs, bool bOmp, bool bDOFtoTypes, char* xyz_name){
    printf( "getEs() imodel %i nDOFs %i bOmp %i bDOFtoTypes %i bEvalOnlyCorrections=%i xyz_name=%s @DOFtoTyp=%p  @Es=%p @Fs=%p \n", W.imodel, W.nDOFs, bOmp, bDOFtoTypes, W.bEvalOnlyCorrections, xyz_name, W.DOFtoTyp, Es, Fs );
    //W.imodel=imodel;
    if(xyz_name){ W.xyz_out=xyz_name; W.bSaveSampleToXYZ=true; }else{ W.bSaveSampleToXYZ=false; }
    if(bDOFtoTypes)W.DOFsToTypes(); 
    W.clean_fDOFs();
    double E = 0;
    if(W.bEvalOnlyCorrections){
       int isamp=0;
       W.printSampleFitSplit(isamp);
       W.printSampleFittedAtoms(isamp);
       double Ecorr=0, Efull=0;  
       W.evalSample_uncorr   ( isamp );
       W.evalSampleError_corr( isamp, Ecorr );
       W.evalSampleError     ( isamp, Efull );
       exit(0);
       W.evalSamples_uncorr();
       E = W.evalSamples_corr( Es ); 
       
    }else
    if( bOmp ){ E = W.evalSamples_omp  ( Es ); }
    else      { E = W.evalSamples_noOmp( Es ); }
    //for(int i=0; i<W.samples.size(); i++){ printf( "getEs() sample[%i] E: %20.10f \n", i, Es[i] ); }
    if( Fs ){ for(int i=0; i<W.nDOFs; i++){ Fs[i] = W.fDOFs[i]; } }
    W.bSaveSampleToXYZ=false; 
    return E;
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

void scanParam2D( int iDOFx, int iDOFy, int nx, int ny, double* xs, double* ys,  double* Es, double* Fx, double* Fy, bool bEvalSamples ){
    //printf( "scanParam() iDOF %i imodel %i nx %i ny %i \n", iDOFx, imodel, nx, ny );
    W.clear_fDOFbounds();
    //W.bRegularize=bRegularize;
    bool bOmp = W.iparallel>0;
    if(bOmp){ W.bBroadcastFDOFs=true; W.realloc_sample_fdofs();  }
    for(int iy=0; iy<ny; iy++){
        W.DOFs[iDOFy] = ys[iy];
        for(int ix=0; ix<nx; ix++){
            int i = iy*nx + ix; 
            //printf( "scanParam2D() ix %i iy %i i %i \n", ix, iy, i );
            W.DOFs[iDOFx] = xs[ix];
            double E = W.evalFitError( i, bOmp );
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

} // extern "C"
