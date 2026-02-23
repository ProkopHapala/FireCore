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

void setModel( int ivdW, int iCoul, int iHbond, int Epairs, int iEpairs, double kMorse, double Lepairs, bool bPN, double svdW, double sCoul, double sHcorr, double sEpairs ){
    W.ivdW    = ivdW;
    W.iCoul   = iCoul;
    W.iHbond  = iHbond;
    if(Epairs>0){W.bEpairs=true;}else if(Epairs<0){W.bEpairs=false;}
    W.iEpairs = iEpairs;
    W.kMorse  = kMorse;
    W.Lepairs = Lepairs;
    W.bPN     = bPN;
    W.svdW    = svdW;
    W.sCoul   = sCoul;
    W.sHcorr  = sHcorr;
    W.sEpairs  = sEpairs;
    printf( "setModel() ivdW %i iCoul %i iHbond %i iEpairs %i kMorse %f Lepairs %f bPN %i svdW %f sCoul %f sHcorr %f sEpairs %f\n", W.ivdW, W.iCoul, W.iHbond, W.iEpairs, W.kMorse, W.Lepairs, W.bPN, W.svdW, W.sCoul, W.sHcorr, W.sEpairs );
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

int getSampleGeom(int isamp, double* pos, int* types, double* Qs, int* host, int* n0){
    if(isamp < 0 || isamp >= W.samples.size()) return 0;
    Atoms* atoms = W.samples[isamp];
    const AddedData* adata = (const AddedData*)(atoms->userData);
    
    W.fillTempArrays( atoms, (Vec3d*)pos, Qs );
    
    for(int i=0; i<atoms->natoms; i++){
        types[i] = atoms->atypes[i];
        if(host) host[i] = adata->host[i];
    }
    if(n0) *n0 = atoms->n0;
    return atoms->natoms;
}

void evalSamplePairs(int isamp, double* pair_out){
    if(isamp < 0 || isamp >= W.samples.size()) return;
    W.DOFsToTypes(); 
    Atoms* atoms = W.samples[isamp];
    const AddedData* adata = (const AddedData*)(atoms->userData);
    int natoms = atoms->natoms;
    alignas(32) double Qs[natoms];
    alignas(32) Vec3d apos[natoms];
    W.fillTempArrays( atoms, apos, Qs );
    
    int i0 = atoms->n0;
    int ni = natoms - atoms->n0;
    int j0 = 0;
    int nj = atoms->n0;

    for(int i=0; i<natoms*natoms*4; i++) pair_out[i] = 0.0;

    for(int ii=0; ii<ni; ii++){
        const int i = i0+ii;
        const int ih = adata->host[i];
        const bool bEpi = ih>=0;
        const Vec3d& pi = apos[i];
        const double Qi = Qs[i];
        const int ti = atoms->atypes[i];
        const Quat4d& REQi = W.typeREQs[ti];

        for(int jj=0; jj<nj; jj++){
            const int j = j0+jj;
            const int jh = adata->host[j];
            const bool bEpj = jh>=0;
            const double Qj = Qs[j];
            const Vec3d dij = apos[j] - pi;
            const int tj = atoms->atypes[j];
            const Quat4d& REQj = W.typeREQs[tj];
            const double R0 = REQi.x + REQj.x;
            const double eps = REQi.y * REQj.y;
            const double Q = Qi * Qj;
            double H = REQi.w * REQj.w;
            const double sH = (H<0.0) ? 1.0 : 0.0;
            const double r = dij.norm();

            double Eij_Coul = 0.0, Eij_vdW = 0.0, Eij_Hcorr = 0.0, Eij_Epairs = 0.0;
            double dE_dH=0.0, dE_dR=0.0;
            double fA=0.0, fR=0.0, fH1=0.0, fH2=0.0;

            if( bEpi ){
                if(bEpj) continue;
                if(W.iEpairs==1)      Eij_Epairs = getSR_PN( r, H, REQi.x, dE_dH, dE_dR );
                else if(W.iEpairs==2) Eij_Epairs = getSR2_PN( r, H, REQi.x, dE_dH, dE_dR );
                else if(W.iEpairs==3) Eij_Epairs = getSR3_PN( r, H, REQi.x, dE_dH, dE_dR );
            }else if( bEpj ){
                if(W.iEpairs==1)      Eij_Epairs = getSR_PN( r, H, REQj.x, dE_dH, dE_dR );
                else if(W.iEpairs==2) Eij_Epairs = getSR2_PN( r, H, REQj.x, dE_dH, dE_dR );
                else if(W.iEpairs==3) Eij_Epairs = getSR3_PN( r, H, REQj.x, dE_dH, dE_dR );
            }else{
                if(W.iCoul==1) Eij_Coul = Q * COULOMB_CONST / r;
                else if(W.iCoul==2) Eij_Coul = Q * dampCoulomb_SoftClamp(r, W.clamp_y1, W.clamp_y2) * COULOMB_CONST;
                else if(W.iCoul>9) Eij_Coul = Q * dampCoulomb_Boys(r, W.boys_rmin, W.iCoul-10) * COULOMB_CONST;

                if(W.ivdW==1){
                    double u = R0/r; double u3 = u*u*u; fA = u3*u3; fR = fA*fA; fH1=1.0; fH2=2.0;
                }else if(W.ivdW==2){
                    double u = R0/r; double u2 = u*u; fA = u2*u2*u2; fR = fA*u2; fH1=3.0; fH2=4.0;
                }else if(W.ivdW==3){
                    double u = R0/r; double u3 = u*u*u; fA = u3*u3; fR = fA*u3; fH1=2.0; fH2=3.0;
                }else if(W.ivdW==4){
                    double alpha = W.kMorse;
                    if(alpha < 0.0) alpha = 6.0/R0;
                    fA = exp(-alpha*(r-R0)); fR = fA*fA; fH1=1.0; fH2=2.0;
                }else if(W.ivdW==5){
                    double alpha = W.kMorse;
                    if(alpha < 0.0) alpha = 6.0/R0;
                    double u = R0/r; double u3 = u*u*u; double e = exp(-alpha*(r-R0));
                    fA = u3*u3; fR = e*e; fH1=1.0; fH2=2.0;
                }
                Eij_vdW = eps * (fH1*fR - fH2*fA);

                if(sH>0.0){
                    if(W.iHbond==1 || W.iHbond==3) Eij_Hcorr += eps * H * fH1 * fR;
                    if(W.iHbond==2 || W.iHbond==3) Eij_Hcorr += eps * H * fH2 * fA;
                }
            }

            int idx1 = (i * natoms + j) * 4;
            int idx2 = (j * natoms + i) * 4;
            pair_out[idx1+0] = Eij_vdW;
            pair_out[idx1+1] = Eij_Coul;
            pair_out[idx1+2] = Eij_Hcorr;
            pair_out[idx1+3] = Eij_Epairs;
            pair_out[idx2+0] = Eij_vdW;
            pair_out[idx2+1] = Eij_Coul;
            pair_out[idx2+2] = Eij_Hcorr;
            pair_out[idx2+3] = Eij_Epairs;
        }
    }
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
