

#include "globals.h"

#include "testUtils.h"
#include "FitREQ.h"

FitREQ W;

MMFFparams params;

//============================

//const std::vector<std::string>* atomTypeNames = 0;
#include "libUtils.h"

extern "C"{

void setVerbosity( int verbosity_, int idebug_ ){
    verbosity = verbosity_;
    idebug    = idebug_;
}



// void init_types(int ntyp, int* typeMask, double* typREQs, bool bCopy ){
//     W.init_types( ntyp, (Quat4i*)typeMask, (Quat4d*)typREQs, bCopy );
// }

//int loadTypeSelection( const char* fname, int imodel ){    return W.loadTypeSelection( fname, imodel );}

int loadTypeSelection_walls( const char* fname ){
    return W.loadTypeSelection_walls( fname );
}

int loadWeights( const char* fname ){
    return W.loadWeights( fname );
}

// void setSystem( int isys, int na, int* types, double* ps, bool bCopy=false ){
//     W.setSystem( isys, na, types, (Vec3d*)ps, bCopy );
// }

// void setRigidSamples( int n, double* Es_, double* poses_, bool bCopy, bool bAlloc ){
//     W.setRigidSamples( n, Es_, (Mat3d*)poses_, bCopy );
// }

// int loadXYZ( const char* fname, int n0, int* i0s, int ntest, int* itests, int* types0, int* testtypes, const char* fname_AtomTypes  ){
//     bool bReadTypes = !(types0 && testtypes);
//     if( bReadTypes && !W.params ){ params.loadAtomTypes( fname_AtomTypes ); W.params=&params; }
//     return W.loadXYZ( fname, n0, i0s, ntest, itests, types0, testtypes );
// }

// void loadTypes( const char* fname_ElemTypes, const char* fname_AtomTypes ){
//     params.loadElementTypes( fname_ElemTypes );
//     params.loadAtomTypes( fname_AtomTypes ); 
//     W.params=&params;
//     W.init_types_par();
// }

void loadTypes_new( const char* fname_ElemTypes, const char* fname_AtomTypes ){
    params.loadElementTypes( fname_ElemTypes );
    params.loadAtomTypes( fname_AtomTypes ); 
    W.params=&params;
}

int loadXYZ_new( const char* fname, bool bAddEpairs, bool bOutXYZ ){
    return W.loadXYZ_new( fname, bAddEpairs, bOutXYZ );
    //return W.loadXYZ_new_bak( fname, bAddEpairs, bOutXYZ );
}

void setWeights( int n, double* weights ){
   _realloc( W.weights, n );
   for(int i=0; i<n; i++){ W.weights[i]=weights[i]; }
   //for(int i=0; i<n; i++){ printf("init.weights[%i]=%g\n",i,W.weights[i]); }
   //W.renormWeights(n);   // Prokop: it was not working because W.nbatch was =0 (modified loadXYZ_new to set nbatch=samples.size() ) 
   //for(int i=0; i<n; i++){ printf("lib.weights[%i]=%g\n",i,W.weights[i]); }       
}

int export_Erefs( double* Erefs ){ return W.export_Erefs( Erefs ); }

double run( int nstep, double Fmax, double dt, int imodel, int isampmode, int ialg, bool bRegularize, bool bClamp, double max_step, bool bEpairs ){
    return W.run( nstep, Fmax, dt, imodel, isampmode, ialg, bRegularize, bClamp, max_step, bEpairs );
}

void setTypeToDOFs  (int i, double* REQ ){ W.setTypeToDOFs  ( i, *(Quat4d*)REQ ); }
void getTypeFromDOFs(int i, double* REQ ){ W.getTypeFromDOFs( i, *(Quat4d*)REQ ); }

double getEs( int imodel, double* Es, int isampmode, bool bEpairs ){
    //printf("USED DOFs= ");for(int j=0;j<W.nDOFs;j++){ printf("%g ",W. DOFs[j]); };printf("\n");
    W.imodel=imodel;
    W.bEpairs=bEpairs;
    switch(isampmode){
        //case 0: return W.evalDerivsRigid( Es ); break;
        //case 1: return W.evalDerivs     ( Es ); break;
        case 2: return W.evalDerivsSamp ( Es ); break;
    }
    //printf("FINAL DOFs= ");for(int j=0;j<W.nDOFs;j++){ printf("%g ",W. DOFs[j]); };printf("\n");
    return 0;
}

void scanParam( int iDOF, int imodel,  int n, double* xs,  double* Es, double* Fs, bool bRegularize ){
    //printf( "scanParam() iDOF %i imodel %i n %i \n", iDOF, imodel, n );
    W.imodel=imodel;
    for(int i=0; i<n; i++){
        W.DOFs[iDOF] = xs[i];
        W.DOFsToTypes();
        double E = W.evalDerivsSamp();
        //if(bRegularize){ W.regularization_force_walls(); }
        if(bRegularize){ W.regularizeDOFs(); }
        if(Fs)Es[i] = E;
        if(Fs)Fs[i] = W.fDOFs[iDOF];
        //if( verbosity>1){ printf( "scanParam()[%3i] W.DOFs[%3i]: %20.10f E: %20.10f  F: %20.10f \n", i, iDOF, W.DOFs[iDOF], E, W.fDOFs[iDOF] ); }
    }
}

void scanParam2D( int iDOFx, int iDOFy, int imodel, int nx, int ny, double* xs, double* ys,  double* Es, double* Fx, double* Fy, bool bRegularize ){
    printf( "scanParam() iDOF %i imodel %i nx %i ny %i \n", iDOFx, imodel, nx, ny );
    W.imodel=imodel;
    for(int iy=0; iy<ny; iy++){
        W.DOFs[iDOFy] = ys[iy];
        for(int ix=0; ix<nx; ix++){
            int i = iy*nx + ix; 
            //printf( "scanParam2D() ix %i iy %i i %i \n", ix, iy, i );
            W.DOFs[iDOFx] = xs[ix];
            W.DOFsToTypes();
            double E = W.evalDerivsSamp();
            if(bRegularize){ W.regularizeDOFs(); }
            if(Es)Es[i] = E;
            if(Fx)Fx[i] = W.fDOFs[iDOFx];
            if(Fy)Fy[i] = W.fDOFs[iDOFy];
        }
    }
}


// void init_buffers(){
//     ibuffers.insert( { "ndims", (int*)&W.nDOFs } );
//     buffers .insert( { "DOFs",   (double*)W.DOFs  } );
//     buffers .insert( { "fDOFs",  (double*)W.fDOFs } );
    
//     ibuffers.insert( { "typToREQ",   (int*)   W.typToREQ  } );
//     buffers .insert( { "typeREQs",   (double*)W.typeREQs  } );
//     buffers .insert( { "typeREQs0",  (double*)W.typeREQs0 } );
//     buffers .insert( { "typeREQsMin",(double*)W.typeREQsMin } );
//     buffers .insert( { "typeREQsMax",(double*)W.typeREQsMax } );
//     buffers .insert( { "typeKreg",   (double*)W.typeKreg  } );
//     //buffers .insert( { "weights",          W.weights  } );
//     buffers .insert( { "Es",               W.Es       } );
//     //if(W.poses)
//     buffers .insert( { "poses",   (double*)W.poses  } );
//     buffers .insert( { "params",  &W.Kneutral       } );
//     //printf( "init_buffers() @Es %li @poses %li \n", (long)W.Es, (long)W.poses  );
//     //printf( "init_buffers() W.system0 @ %li %li %li \n", (long)W.system0, (long)W.systemTest0, (long)W.systemTest );
//     buffers .insert( { "ps1",    (double*)W.system0    ->apos  } );
//     //buffers .insert( { "ps2",    (double*)W.systemTest0->apos  } );
//     //buffers .insert( { "ps3",    (double*)W.systemTest ->apos  } );
//     ibuffers.insert( { "types1", (int*)W.system0    ->atypes } );
//     //ibuffers.insert( { "types2", (int*)W.systemTest0->atypes } );
//     //ibuffers.insert( { "types3", (int*)W.systemTest ->atypes } );
// }



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
