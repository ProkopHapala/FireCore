﻿
#include  "globals.h"

#include "testUtils.h"
#include "MolWorld_sp3.h"

// ============ Global Variables

MolWorld_sp3 W;

//============================

#include "libMMFF.h"
#include "libUtils.h"

extern "C"{

void init_buffers(){
    //printf( "init_buffers() \n" );
    buffers .insert( { "apos",   (double*)W.nbmol.apos } );
    buffers .insert( { "fapos",  (double*)W.nbmol.fapos } );
    if(W.bMMFF){
        buffers .insert( { "DOFs",      W.ffl.DOFs  } );
        buffers .insert( { "fDOFs",     W.ffl.fDOFs } );
        buffers .insert( { "vDOFs",     W.opt.vel  } );
        if(!W.bUFF){
            buffers .insert( { "pipos",  (double*)W.ffl.pipos   } );
            buffers .insert( { "fpipos", (double*)W.ffl.fpipos } );
            ibuffers.insert( { "neighs",      (int*)W.ffl.neighs  } );
        }
    }else{
        W.ff.natoms=W.nbmol.natoms;
    }
    ibuffers.insert( { "ndims",    &W.ff.nDOFs } );
    buffers .insert( { "Es",       &W.ff.Etot  } );
    ibuffers.insert( { "selection", W.manipulation_sel  } );
    bbuffers.insert( { "ffflags", &W.doBonded  } );
    //printBuffNames();
}

// int loadmol(char* fname_mol ){ return W.loadmol(fname_mol ); }
//lib.init( cstr(xyz_name), cstr(surf_name), cstr(smile_name),      bMMFF,      bEpairs,      bUFF,      b141,      bSimple,      bConj,      bCumulene,      nPBC,        gridStep, cstr(sElementTypes), cstr(sAtomTypes), cstr(sBondTypes), cstr(sAngleTypes), cstr(sDihedralTypes) )
void* init( char* xyz_name, char* surf_name, char* smile_name, bool bMMFF, bool bEpairs, bool bUFF, bool b141, bool bSimple, bool bConj, bool bCumulene, int* nPBC, double gridStep, char* sElementTypes, char* sAtomTypes, char* sBondTypes, char* sAngleTypes, char* sDihedralTypes ){
	W.smile_name = smile_name;
	W.xyz_name   = xyz_name;
	W.surf_name  = surf_name;
	W.bMMFF      = bMMFF;
    W.bEpairs    = bEpairs;
    W.gridStep   = gridStep;
    W.nPBC       = *(Vec3i*)nPBC;
    W.bUFF       = bUFF; 
    W.b141       = b141;
    W.bSimple    = bSimple;
    W.bConj      = bConj;
    W.bCumulene  = bCumulene;
    // read and store parameters from tables
    // TBD pass bUFF to MMFFparams::init so that if true, no need to read bonds, angles nor dihedrals...
    //W.params.verbosity = verbosity;
    W.params.init( sElementTypes, sAtomTypes, sBondTypes, sAngleTypes, sDihedralTypes );
    // bring names of atom types into builder (H is capping atom, E is electron pair)
	W.builder.bindParams(&W.params);
    bool bGrid = gridStep>0;
    // initialize the main
    //W.init( bGrid, bUFF );
    W.bGridFF=bGrid;
    W.bUFF   =bUFF;
    W.init();
    init_buffers();
    return &W;
}

int    run( int nstepMax, double dt, double Fconv, int ialg, double damping, double* outE, double* outF, double* outV, double* outVF, bool omp ){
    //printf( "bOpenMP = %i \n", omp );
    //W.rum_omp_ocl( nstepMax, dt, Fconv, 1000.0, 1000 ); 
    // run_omp( int niter_max, double dt, double Fconv=1e-6, double Flim=1000, double timeLimit=0.02, double* outE=0, double* outF=0 ){
    if(omp){ return W.run_omp   (nstepMax,dt,Fconv,   10.0, -1.0, outE, outF, outV, outVF ); }
    else   { return W.run_no_omp(nstepMax,dt,Fconv, 1000.0,  damping, outE, outF, outV, outVF ); }
    //else   { return W.run       (nstepMax,dt,Fconv,ialg,       outE, outF, outV, outVF ); }
}

int substituteMolecule( const char* fname, int ib, double* up, int ipivot, bool bSwapBond ){
    return W.substituteMolecule( fname, ib, *(Vec3d*)up, ipivot, bSwapBond );
}

void set_opt( 
        double dt_max,  double dt_min, double damp_max, 
        double finc,    double fdec,   double falpha, int minLastNeg,
        double cvf_min, double cvf_max
    ){

    W.opt.dt_max  = dt_max;
    W.opt.dt_min  = dt_min;
    W.opt.dt      = dt_max;

    W.opt.damp_max   = damp_max;
    W.opt.damping    = damp_max;

    W.opt.cvf_min    = cvf_min;
    W.opt.cvf_max    = cvf_max;

    W.opt.minLastNeg =  minLastNeg;
    W.opt.finc       =  finc;
    W.opt.fdec       =  fdec;
    W.opt.falpha     =  falpha;
    
    //W.opt.f_limit  =  f_limit  ;
    //W.opt.v_limit  =  v_limit  ;
    //W.opt.dr_limit =  dr_limit ;

}

void sampleSurf(char* name, int n, double* rs, double* Es, double* fs, int kind, int atyp, double Q, double K, double RQ, double* pos0_, bool bSave){
    if(name){
        W.ff.realloc( 1,0,0,0, true );
        W.ff.apos [0] = *(Vec3d*)pos0_;
        W.ff.atype[0] = atyp;
        bool bGrid=(kind>10);
        if( kind==10 ) W.gridFF.iDebugEvalR=1;
        W.gridFF.alphaMorse = K;
        W.gridFF.Rdamp = RQ;
        W.loadSurf( name, bGrid, bSave );
        W.nbmol.REQs[0].z = Q;
        if(bSave){
            Quat4f* FFtot = new Quat4f[W.gridFF.grid.getNtot()];
            W.gridFF.evalCombindGridFF ( W.nbmol.REQs[0], FFtot );
            W.gridFF.grid.saveXSF<float>( "FFtot_E.xsf", (float*)FFtot, 4, 3, W.surf.natoms, W.surf.atypes, W.surf.apos );
            delete [] FFtot;
        }
    }
    Quat4d REQ=W.nbmol.REQs[0];
    Quat4f PLQ = REQ2PLQ( REQ, K );
    printf( "DEBUG sampleSurf REQ(%g,%g,%g) \n", REQ.x, REQ.y, REQ.z );
    printf( "DEBUG sampleSurf PLQ(%g,%g,%g) \n", PLQ.x, PLQ.y, PLQ.z );
    //exit(0);
    double R2Q=RQ*RQ;
    for(int i=0; i<n; i++){
        Quat4f fe=Quat4fZero;
        W.nbmol.apos[0].z=rs[i];
        W.ff.cleanAtomForce();
        switch(kind){
            case  0: fe.e=   W.nbmol.evalR         (W.surf ); break; 
            case  1: fe.e=   W.nbmol.evalMorse     (W.surf, false,                           K,RQ  ); fe.f=(Vec3f)W.nbmol.fapos[0]; break; 
            //case  5: fe.e=   W.nbmol.evalMorsePLQ  (W.surf, PLQ, W.gridFF.grid.cell, {1,1,0},K,R2Q ); fe.f=(Vec3f)W.nbmol.fapos[0]; break; 
            case 10:         W.gridFF.addForce_surf(W.nbmol.apos[0], {1.,0.,0.}, fe );  break;
            case 11:         W.gridFF.addForce_surf(W.nbmol.apos[0], PLQ, fe );  break;
            case 12:         W.gridFF.addForce     (W.nbmol.apos[0], PLQ, fe );  break;
        }
        fs[i]=fe.z;
        Es[i]=fe.e;
    }
}

void print_debugs( bool bParams, bool bNeighs, bool bShifts ){
    W.ffl.printSizes();
    if( bParams ) W.ffl.printAtomParams();
    if( bNeighs ) W.ffl.printNeighs();
    if( bShifts ) W.ffl.print_pbc_shifts();
}

void sampleSurf_vecs(char* name, int n, double* poss_, double* Es, double* fs_, int kind, int atyp, double Q, double K, double RQ, double* pos0_, bool bSave){
    Vec3d* poss =(Vec3d*)poss_;
    Vec3d* fs   =(Vec3d*)fs_;
    if(name){
        W.ff.realloc( 1,0,0,0, true );
        W.ff.apos [0] = *(Vec3d*)pos0_;
        W.ff.atype[0] = atyp;
        bool bGrid=(kind>=10);
        if( kind==10 ) W.gridFF.iDebugEvalR=1;
        W.gridFF.alphaMorse = K;
        W.gridFF.Rdamp = RQ;
        W.loadSurf( name, bGrid, bSave );
        W.nbmol.REQs[0].z = Q;
        if(bSave){
            Quat4f* FFtot = new Quat4f[W.gridFF.grid.getNtot()];
            W.gridFF.evalCombindGridFF ( W.nbmol.REQs[0], FFtot );
            W.gridFF.grid.saveXSF<float>( "FFtot_E.xsf", (float*)FFtot, 4, 3, W.surf.natoms, W.surf.atypes, W.surf.apos );
            printf( "DEBUG saveXSF() DONE \n" );
            delete [] FFtot;
        }
    }
    printf( "DEBUG start sampling kind=%i \n", kind );
    Quat4f PLQ = REQ2PLQ( W.nbmol.REQs[0], K );
    //printf( "PLQ(%g,%g,%g) \n", PLQ.x, PLQ.y, PLQ.z );
    double R2Q=RQ*RQ;
    for(int i=0; i<n; i++){
        Quat4f fe=Quat4fZero;
        W.nbmol.apos[0]=poss[i];
        //printf( "[%i] (%g,%g,%g)\n", i, W.nbmol.apos[0].x,W.nbmol.apos[0].y,W.nbmol.apos[0].z );
        W.ff.cleanAtomForce();
        switch(kind){
            case  0: fe.e=   W.nbmol.evalR         (W.surf ); break; 
            case  1: fe.e=   W.nbmol.evalMorse     (W.surf, false,                           K,RQ  ); fe.f=(Vec3f)W.nbmol.fapos[0]; break; 
            //case  5: fe.e=   W.nbmol.evalMorsePLQ  (W.surf, PLQ, W.gridFF.grid.cell, {1,1,0},K,R2Q ); fe.f=(Vec3f)W.nbmol.fapos[0]; break; 
            case 10:         W.gridFF.addForce_surf(W.nbmol.apos[0], {1.,0.,0.}, fe );  break;
            case 11:         W.gridFF.addForce_surf(W.nbmol.apos[0], PLQ, fe );  break;
            case 12:         W.gridFF.addForce     (W.nbmol.apos[0], PLQ, fe );  break;
        }
        fs[i]=(Vec3d)(fe.f);
        Es[i]=fe.e;
    }
}

void change_lvec( double* lvec, bool bAdd, bool  ){
    if(bAdd){ W.add_to_lvec( *(Mat3d*)lvec ); }
    else    { W.change_lvec( *(Mat3d*)lvec ); }
}

void optimizeLattice_1d( double* dlvec, int n1, int n2, int initMode, double tol ){
    printf("MMFF_lib::optimizeLattice_1d(n1=%i,n2=%i,initMode=%i,tol=%g) \n", n1, n2, initMode, tol );
    W.gopt.tolerance=tol;
    W.gopt.initMode =initMode; 
    W.optimizeLattice_1d( n1, n2, *(Mat3d*)dlvec );
}

void addSnapshot(){
    W.addSnapshotIfNew();
}

void printDatabase(){
    W.printDatabase();
}

} // extern "C"
