

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

double run( int nstep, double Fmax, double dt, int imodel, int isampmode, int ialg, bool bRegularize, bool bClamp, double max_step, bool bEpairs ){
    return W.run( nstep, Fmax, dt, imodel, isampmode, ialg, bRegularize, bClamp, max_step, bEpairs );
}


// double run( int nstep, double Fmax, double dt, int imodel, int isampmode, int ialg, bool bRegularize, bool bClamp, double max_step, bool bEpairs ){
//     W.imodel=imodel;
//     W.bEpairs=bEpairs;
//     double Err=0;
//     printf( "run( nstep %i Fmax %g dt %g isamp %i )\n", nstep, Fmax, dt, isampmode  );
//     double F2max=Fmax*Fmax;
//     double F2;
//     for(int i=0; i<nstep; i++){
//         //printf("[%i]  DOFs=", i);for(int j=0;j<W.nDOFs;j++){ printf("%g ",W. DOFs[j]); };printf("\n");
//         W.DOFsToTypes(); 
//         W.clean_derivs();
//         //printf("[%i]  DOFs=", i);for(int j=0;j<W.nDOFs;j++){ printf("%g ",W. DOFs[j]); };printf("\n");
//         switch(isampmode){
//             case 0: Err = W.evalDerivsRigid(); break;
//             case 1: Err = W.evalDerivs     (); break;
//             case 2: Err = W.evalDerivsSamp (); break;
//         }   
//         printf("step= %i DOFs= ", i);for(int j=0;j<W.nDOFs;j++){ printf("%g ",W.DOFs[j]); };printf("\n");
//         //if(bRegularize){ W.regularization_force(); }
//         printf("step= %i fDOFs= ", i);for(int j=0;j<W.nDOFs;j++){ printf("%g ",W.fDOFs[j]); };printf("\n");
//         if(bRegularize){ W.regularization_force_walls(); }
//         printf("step= %i after_reg fDOFs= ", i);for(int j=0;j<W.nDOFs;j++){ printf("%g ",W.fDOFs[j]); };printf("\n");
// //exit(0);        
//         switch(ialg){
//             case 0: F2 = W.move_GD( dt, max_step ); break;
//             case 1: F2 = W.move_MD( dt, max_step ); break;
//             case 2: F2 = W.move_GD_BB_short( i, dt, max_step ); break;
//             case 3: F2 = W.move_GD_BB_long( i, dt, max_step ); break;
//             case 4: F2 = W.move_MD_nodamp( dt, max_step ); break;
//         }
//         // regularization must be done before evaluation of derivatives
//         if(bClamp     ){ W.limit_params();         }
//         //printf("step= %i dt= %g\n", i, dt );
//         printf("step= %i RMSE= %g |F|= %g\n", i, sqrt(Err), sqrt(abs(F2)) );
//         printf("[%i]\n", i );
//         if( F2<0.0 ){ printf("DYNAMICS STOPPED after %i iterations \n", i); printf("VERY FINAL DOFs= ");for(int j=0;j<W.nDOFs;j++){ printf("%.15g ",W. DOFs[j]); };printf("\n"); return Err; }
//         if( F2<F2max ){ printf("CONVERGED in %i iterations \n", i); printf("VERY FINAL DOFs= ");for(int j=0;j<W.nDOFs;j++){ printf("%.15g ",W. DOFs[j]); };printf("\n"); return Err; }
//     }
//     printf("step= %i DOFs= ", nstep);for(int j=0;j<W.nDOFs;j++){ printf("%g ",W. DOFs[j]); };printf("\n");
//     printf("VERY FINAL DOFs= ");for(int j=0;j<W.nDOFs;j++){ printf("%.15g ",W. DOFs[j]); };printf("\n");
//     return Err;
// }

void setType(int i, double* REQ ){ W.setType( i, *(Quat4d*)REQ ); }
void getType(int i, double* REQ ){ W.getType( i, *(Quat4d*)REQ ); }

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

void init_buffers(){

    ibuffers.insert( { "ndims", (int*)&W.nDOFs } );
    
    buffers .insert( { "DOFs",   (double*)W.DOFs  } );
    buffers .insert( { "fDOFs",  (double*)W.fDOFs } );
    
    ibuffers.insert( { "typToREQ",   (int*)   W.typToREQ  } );
    buffers .insert( { "typeREQs",   (double*)W.typeREQs  } );
    buffers .insert( { "typeREQs0",  (double*)W.typeREQs0 } );
    buffers .insert( { "typeREQsMin",(double*)W.typeREQsMin } );
    buffers .insert( { "typeREQsMax",(double*)W.typeREQsMax } );
    buffers .insert( { "typeKreg",   (double*)W.typeKreg  } );

    //buffers .insert( { "weights",          W.weights  } );
    buffers .insert( { "Es",               W.Es       } );
    //if(W.poses)
    buffers .insert( { "poses",   (double*)W.poses  } );
    buffers .insert( { "params",  &W.Kneutral       } );

    //printf( "init_buffers() @Es %li @poses %li \n", (long)W.Es, (long)W.poses  );

    //printf( "init_buffers() W.system0 @ %li %li %li \n", (long)W.system0, (long)W.systemTest0, (long)W.systemTest );
    buffers .insert( { "ps1",    (double*)W.system0    ->apos  } );
    //buffers .insert( { "ps2",    (double*)W.systemTest0->apos  } );
    //buffers .insert( { "ps3",    (double*)W.systemTest ->apos  } );
    
    ibuffers.insert( { "types1", (int*)W.system0    ->atypes } );
    //ibuffers.insert( { "types2", (int*)W.systemTest0->atypes } );
    //ibuffers.insert( { "types3", (int*)W.systemTest ->atypes } );
    
}

} // extern "C"
