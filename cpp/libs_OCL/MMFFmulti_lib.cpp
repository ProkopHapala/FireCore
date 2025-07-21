
#include <globals.h>

#include "testUtils.h"

//#include "MolWorld_sp3.h"
#include "MolWorld_sp3_multi.h"

#include <thread>
#include <chrono>

// ============ Global Variables

MolWorld_sp3_multi W;

//============================

#include "libMMFF.h"
#include "libUtils.h"

extern "C"{

void init_buffers(){
    buffers .insert( { "apos",   (double*)W.nbmol.apos } );
    buffers .insert( { "fapos",  (double*)W.nbmol.fapos } );
    if(W.bMMFF){
        buffers .insert( { "DOFs",      W.ffl.DOFs  } );
        buffers .insert( { "fDOFs",     W.ffl.fDOFs } );
        buffers .insert( { "vDOFs",     W.opt.vel  } );
        //buffers .insert( { "apos",   (double*)W.ff.apos   } );
        //buffers .insert( { "fapos",  (double*)W.ff.fapos } );
        buffers .insert( { "apos",   (double*)W.nbmol.apos } );
        buffers .insert( { "fapos",  (double*)W.nbmol.fapos } );
        buffers .insert( { "pipos",  (double*)W.ffl.pipos   } );
        buffers .insert( { "fpipos", (double*)W.ffl.fpipos } );
        //buffers .insert( { "bond_l0",   (double*)W.ffl.bond_l0   } );
        //buffers .insert( { "bond_k",    (double*)W.ffl.bond_k    } );
        //buffers .insert( { "pbcShifts", (double*)W.ff.pbcShifts } );
        //buffers .insert( { "Kneighs",   (double*)W.ff.Kneighs   } );
        //ibuffers.insert( { "bond2atom",    (int*)W.ff.bond2atom  } );
        ibuffers.insert( { "neighs",       (int*)W.ffl.neighs  } );

        ibuffers.insert( { "gpu_neighs",   (int*)W.neighs    } );
        ibuffers.insert( { "gpu_neighCell",(int*)W.neighCell } );
        ibuffers.insert( { "gpu_bkNeighs", (int*)W.bkNeighs  } );

        fbuffers.insert( { "gpu_atoms",    (float*)W.atoms   } );
        fbuffers.insert( { "gpu_aforces",  (float*)W.aforces } );
        fbuffers.insert( { "gpu_avel",     (float*)W.avel    } );
        fbuffers.insert( { "gpu_constr",   (float*)W.constr  } );
        
        fbuffers.insert( { "gpu_REQs",     (float*)W.REQs   } );
        fbuffers.insert( { "gpu_MMpars",   (float*)W.MMpars } );
        fbuffers.insert( { "gpu_BLs",      (float*)W.BLs    } );
        fbuffers.insert( { "gpu_BKs",      (float*)W.BKs    } );
        fbuffers.insert( { "gpu_Ksp",      (float*)W.Ksp    } );
        fbuffers.insert( { "gpu_Kpp",      (float*)W.Kpp    } );

        fbuffers.insert( { "gpu_lvecs",    (float*)W.lvecs     } );
        fbuffers.insert( { "gpu_ilvecs",   (float*)W.ilvecs    } );
        fbuffers.insert( { "gpu_pbcshifts",(float*)W.pbcshifts } );

        //float* gpu_lvecs= getfBuff("gpu_lvecs"); 
        //for(int i=0; i<12; i++){ printf( "gpu_lvecs[%i]=%g\n", i, gpu_lvecs[i] ); };
    }else{
        W.ff.natoms=W.nbmol.natoms;
    }
    ibuffers.insert( { "ndims",    &W.ffl.nDOFs } );
    buffers .insert( { "Es",       &W.ffl.Etot  } );
    ibuffers.insert( { "selection", W.manipulation_sel  } );
}

// int loadmol(char* fname_mol ){ return W.loadmol(fname_mol ); }

void* init( int nSys, char* xyz_name, char* surf_name, char* smile_name, bool bMMFF, bool bEpairs, int* nPBC, int* grid_nPBC, double gridStep, char* sElementTypes, char* sAtomTypes, char* sBondTypes, char* sAngleTypes, double T, double gamma, int nExplore, int nRelax, double pos_kick, double vel_kick, int GridFF ){
    printf( "MMFFmulti_lib::init() nSys=%i xyz_name(%s) surf_name(%s) bMMFF=%i bEpairs=%i T=%g gamma=%g nExplore=%i nRelax=%i pos_kick=%g vel_kick=%g \n", nSys, xyz_name, surf_name, bMMFF, bEpairs, T, gamma, nExplore, nRelax, pos_kick, vel_kick );
    W.smile_name = smile_name;
    W.xyz_name   = xyz_name;
    W.surf_name  = surf_name;
    W.bMMFF      = bMMFF;
    W.bEpairs    = bEpairs;
    W.gridStep   = gridStep;
    W.nPBC       = *(Vec3i*)nPBC;
    W.gridFF.nPBC= *(Vec3i*)grid_nPBC;
//    W.params.init( sElementTypes, sAtomTypes, sBondTypes, sAngleTypes );
//	W.builder.bindParams(&W.params);
    W.nSystems=nSys;
    bool bGrid = GridFF>0;
//    
    W.bGopt = true;
    W.go.T_target = T;
    W.go.gamma_damp = gamma;
    W.go.nExplore = nExplore;
    W.go.nRelax = nRelax;
    W.go.pos_kick = pos_kick;
    W.go.vel_kick = vel_kick;

    long T0=getCPUticks();
    float nseconds = 0.1;
    std::this_thread::sleep_for(std::chrono::milliseconds((int)(1000*0.1)) );
    tick2second = nseconds/(getCPUticks()-T0);

    if(surf_name){
        printf("surf_name(%s)\n", surf_name);
    }
    else{
        printf("No surface file specified\n");
    }
    W.init();
    W.bGridFF = bGrid;
    printf("GridFF=%i\n", W.bGridFF);
    init_buffers();
    return &W;
}

int run( int nstepMax, double dt, double Fconv, int ialg, double* outE, double* outF, int iParalel ){
    W.bOcl= iParalel > 0;
    Mat3d lvec = W.ffls[0].lvec;
    //printf( "run Fconv=%g lvec{{%6.3f,%6.3f,%6.3f}{%6.3f,%6.3f,%6.3f}{%6.3f,%6.3f,%6.3f}}\n", Fconv, lvec.a.x,lvec.a.y,lvec.a.z, lvec.b.x,lvec.b.y,lvec.b.z, lvec.c.x,lvec.c.y,lvec.c.z );
    int nitrdione=0;
    switch(iParalel){
        case -1: nitrdione = W.run_multi_serial( nstepMax, Fconv, 1000.0, 1000 ); break; 
        case  0:
        case  1: nitrdione = W.run_omp_ocl     ( nstepMax, Fconv, 1000.0, 1000 ); break; 
        case  2: nitrdione = W.run_ocl_opt( nstepMax, Fconv    ); break; 
        case  3: nitrdione = W.run_ocl_loc( nstepMax, Fconv, 1 ); break; 
        case  4: nitrdione = W.run_ocl_loc( nstepMax, Fconv, 2 ); break; 
    }
    return nitrdione;
    //return W.rum_omp_ocl( nstepMax, dt, Fconv, 1000.0, 1000 ); 
    //return W.run(nstepMax,dt,Fconv,ialg,outE,outF);  
}

void MDloop( int perframe, double Ftol = -1, int iParalel = 3, int perVF = 100 ){
    W.iParalel = iParalel;
    W.nPerVFs = perVF;
    W.iterPerFrame = perframe;
    W.MDloop( perframe, Ftol );
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


void pack_system( int isys, bool bParams, bool bForces, bool bVel, bool blvec ){
    W.pack_system( isys, W.ffls[isys], bParams, bForces, bVel, blvec );
}
void unpack_system( int isys, bool bForces=false, bool bVel=false ){
    W.unpack_system( isys, W.ffls[isys], bForces, bVel );
}

void upload_sys  ( int isys, bool bParams, bool bForces, bool bVel, bool blvec ){
    W.upload_sys  ( isys, bParams, bForces, bVel, blvec );
}
void download_sys( int isys, bool bForces, bool bVel ){
    W.download_sys( isys, bForces, bVel );
}

void upload(  bool bParams, bool bForces, bool bVel, bool blvec ){
    W.upload( bParams, bForces, bVel, blvec );
}
void download( bool bForces, bool bVel ){
    W.download( bForces, bVel );
}


void change_lvec( double* lvec, bool bAdd, bool bUpdatePi ){
    if(bAdd){ W.add_to_lvec( *(Mat3d*)lvec ); }
    else    { W.change_lvec( *(Mat3d*)lvec ); }
    //if(bUpdatePi){ W.ffl.initPi( W.pbc_shifts ); }
}

void upload_pop( const char* fname ){
    printf("MolWorld_sp3::upload_pop(%s)\n", fname );
    W.gopt.loadPopXYZ( fname );
    int nmult=W.gopt.population.size(); 
    int npara=W.paralel_size(); if( nmult!=npara ){ printf("ERROR in GlobalOptimizer::lattice_scan_1d_multi(): (imax-imin)=(%i) != solver.paralel_size(%i) => Exit() \n", nmult, npara ); exit(0); }
    W.gopt.upload_multi(nmult,0, true, true );
}

void optimizeLattice_1d( double* dlvec, int n1, int n2, int initMode, double tol ){
    printf("MMFFmulti_lib::optimizeLattice_1d(n1=%i,n2=%i,initMode=%i,tol=%g) \n", n1, n2, initMode, tol );
    W.optimizeLattice_1d( n1, n2, *(Mat3d*)dlvec );
}

void optimizeLattice_2d_multi( double* dlvec, int nstesp, int initMode, double tol ){
    //int initMode=1;
    //int initMode = 0;
    //int nstesp   = 40;
    //int nstesp = 2;
    //gopt.tolerance = 0.02;
    //W.gopt.tolerance = 0.01;
    //Mat3d dlvec = Mat3d{   0.2,0.0,0.0,    0.0,0.0,0.0,    0.0,0.0,0.0  };
    W.gopt.lattice_scan_2d_multi( nstesp, *(Mat3d*)dlvec, initMode, "lattice_scan_2d_multi.xyz" );
}

void  scan( int nconf, double* poss, double* rots, double* dirs, double* Es, double* aforces, double* aposs, bool omp, bool bRelax, int niter_max, double dt, double Fconv, double Flim ){
    // if(bRelax){
    //     if(omp){ printf("ERROR: scan_relaxed() not implemented witht OMP\n"); exit(0); } 
    //     else   { W.scan_relaxed( nconf, (Vec3d*)poss, (Mat3d*)rots, Es, (Vec3d*)aforces, (Vec3d*)aposs, omp, niter_max, dt, Fconv, Flim );  }
    // }else{
    //     if(omp){ printf("ERROR: scan_rigid() not implemented witht OMP\n"); exit(0); } 
    //     else   { W.scan_rigid( nconf, (Vec3d*)poss, (Mat3d*)rots, Es, (Vec3d*)aforces, (Vec3d*)aposs, omp ); }
    // }
    W.scan_relaxed( nconf, (Vec3d*)poss, (Mat3d*)rots, (Vec3d*)dirs, Es, (Vec3d*)aforces, (Vec3d*)aposs, omp, niter_max, dt, Fconv, Flim );
}

void setSwitches_multi( int CheckInvariants, int PBC, int NonBonded, int MMFF, int Angles, int PiSigma, int PiPiI, int dovdW){
    #define _setbool(b,i) { if(i>0){b=true;}else if(i<0){b=false;} }
    _setbool( W.bCheckInvariants, CheckInvariants  );
    _setbool( W.bPBC         , PBC       );
    _setbool( W.bNonBonded   , NonBonded );
    _setbool( W.bMMFF        , MMFF      );
    _setbool( W.ffl.doAngles , Angles    );
    _setbool( W.ffl.doPiSigma, PiSigma   );
    _setbool( W.ffl.doPiPiI  , PiPiI     );
    _setbool( W.dovdW, dovdW );
    W.ffl.bSubtractAngleNonBond = W.bNonBonded;
    #undef _setbool
}


} // extern "C"
