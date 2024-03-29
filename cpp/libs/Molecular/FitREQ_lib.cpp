﻿

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

void init_types(int ntyp, int* typeMask, double* typREQs, bool bCopy ){
    W.init_types( ntyp, (Quat4i*)typeMask, (Quat4d*)typREQs, bCopy );
}

void setSystem( int isys, int na, int* types, double* ps, bool bCopy=false ){
    W.setSystem( isys, na, types, (Vec3d*)ps, bCopy );
}

void setRigidSamples( int n, double* Es_, double* poses_, bool bCopy, bool bAlloc ){
    W.setRigidSamples( n, Es_, (Mat3d*)poses_, bCopy );
}

int loadXYZ( const char* fname, int n0, int* i0s, int ntest, int* itests, int* types0, int* testtypes, const char* fname_AtomTypes  ){
    bool bReadTypes = !(types0 && testtypes);
    if( bReadTypes && !W.params ){ params.loadAtomTypes( fname_AtomTypes ); W.params=&params; }
    return W.loadXYZ( fname, n0, i0s, ntest, itests, types0, testtypes );
}

void loadTypes( const char* fname_ElemTypes, const char* fname_AtomTypes ){
    params.loadElementTypes( fname_ElemTypes );
    params.loadAtomTypes( fname_AtomTypes ); 
    W.params=&params;
    W.init_types_par();
}

int loadXYZ_new( const char* fname, bool bAddEpairs, bool bOutXYZ ){
    return W.loadXYZ_new( fname, bAddEpairs, bOutXYZ );
}

double run( int imodel,  int nstep, double Fmax, double dt, bool bRigid , int ialg, bool bRegularize, bool bClamp){
    W.imodel=imodel;
    double Err=0;
    printf( "run( nstep %i Fmax %g dt %g bRigid %i )\n", nstep, Fmax, dt, bRigid  );
    double F2max=Fmax*Fmax;
    for(int i=0; i<nstep; i++){
        W.DOFsToTypes(); 
        W.clean_derivs();
        if(bRigid){ Err= W.evalDerivsRigid(); }
        else      { Err= W.evalDerivs     (); }
        if(bRegularize){ W.regularization_force(); }
        if(bClamp     ){ W.limit_params();         }
        //printf(" DOFs=");for(int j=0;j<W.nDOFs;j++){ printf("%g ",W. DOFs[j]); };printf("\n");
        //printf("fDOFs=");for(int j=0;j<W.nDOFs;j++){ printf("%g ",W.fDOFs[j]); };printf("\n");
        double F2=1;
        switch(ialg){
            case 0: F2 = W.move_GD( dt ); break;
            case 1: F2 = W.move_MD( dt ); break;
        }
        printf("[%i] RMS=%g |F|=%g\n", i, sqrt(Err), sqrt(F2) );
        if( F2<F2max ){ printf("CONVERGED in %i iterations \n", i); break; }
    }
    return Err;
}

void setType(int i, double* REQ ){ W.setType( i, *(Quat4d*)REQ ); }
void getType(int i, double* REQ ){ W.getType( i, *(Quat4d*)REQ ); }

double getEs( int imodel, double* Es, int isampmode ){
    W.imodel=imodel;
    switch(isampmode){
        case 0: return W.evalDerivsRigid( Es ); break;
        case 1: return W.evalDerivs     ( Es ); break;
        case 2: return W.evalDerivsSamp ( Es ); break;
    }
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
