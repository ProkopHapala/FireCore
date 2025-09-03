
#include <globals.h>

#include "testUtils.h"

//#include "MolWorld_sp3.h"
#include "MolWorld_sp3_multi.h"

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

void* init( int nSys, char* xyz_name, char* surf_name, char* smile_name, bool bMMFF, bool bEpairs, bool bUFF, bool b141, bool bSimple, bool bConj, bool bCumulene, int* nPBC, double gridStep, char* sElementTypes, char* sAtomTypes, char* sBondTypes, char* sAngleTypes, char* sDihedralTypes ){
    printf( "MMFFmulti_lib::init() nSys=%i xyz_name(%s) surf_name(%s) bMMFF=%i bEpairs=%i bUFF=%i \n", nSys, xyz_name, surf_name, bMMFF, bEpairs, bUFF );
	W.smile_name = smile_name;
	W.xyz_name   = xyz_name;
	W.surf_name  = surf_name;
	W.bMMFF      = bMMFF;
    W.bUFF       = bUFF;
    W.bEpairs    = bEpairs;
    W.gridStep   = gridStep;
    W.nPBC       = *(Vec3i*)nPBC;
    W.params.init( sElementTypes, sAtomTypes, sBondTypes, sAngleTypes, sDihedralTypes );
	W.builder.bindParams(&W.params);
    W.nSystems=nSys;
    bool bGrid = gridStep>0;
    W.bGridFF = bGrid;
    W.init();
    init_buffers();
    return &W;
}

void setSwitchesUFF( int DoBond, int DoAngle, int DoDihedral, int DoInversion, int DoAssemble, int SubtractBondNonBond, int ClampNonBonded ){
    #define _setbool(b,i) { if(i>0){b=true;}else if(i<0){b=false;} }
    if(W.uff_ocl){
        _setbool( W.uff_ocl->bUFF_bonds,      DoBond );
        _setbool( W.uff_ocl->bUFF_angles,     DoAngle );
        _setbool( W.uff_ocl->bUFF_dihedrals,  DoDihedral );
        _setbool( W.uff_ocl->bUFF_inversions, DoInversion );
        _setbool( W.uff_ocl->bSubtractNB,     SubtractBondNonBond );
        _setbool( W.uff_ocl->bClampNonBonded, ClampNonBonded );
    }
    _setbool( W.ffu.bDoBond,              DoBond );
    _setbool( W.ffu.bDoAngle,             DoAngle );
    _setbool( W.ffu.bDoDihedral,          DoDihedral );
    _setbool( W.ffu.bDoInversion,         DoInversion );
    _setbool( W.ffu.bDoAssemble,          DoAssemble );
    _setbool( W.ffu.bSubtractBondNonBond, SubtractBondNonBond );
    _setbool( W.ffu.bClampNonBonded,      ClampNonBonded );
    #undef _setbool
}

int run( int nstepMax, double dt, double Fconv, int ialg, double* outE, double* outF, int iParalel ){
    W.bOcl= iParalel > 0;
    Mat3d lvec = W.ffls[0].lvec;
    //printf( "run Fconv=%g lvec{{%6.3f,%6.3f,%6.3f}{%6.3f,%6.3f,%6.3f}{%6.3f,%6.3f,%6.3f}}\n", Fconv, lvec.a.x,lvec.a.y,lvec.a.z, lvec.b.x,lvec.b.y,lvec.b.z, lvec.c.x,lvec.c.y,lvec.c.z );
    int nitrdione=0;
    if(W.bUFF){
        // UFF branch: CPU (serial/OpenMP) vs GPU (OpenCL) based on iParalel
        switch(iParalel){
            // CPU paths (use UFF::run on host)
            case -1:
            case  0: {
                // Run UFF on CPU for system 0 (host-side integrator)
                // Note: UFF::run updates ffu.apos/fapos internally
                nitrdione = W.ffu.run( nstepMax, dt, Fconv, 1000.0, 0.1, outE, outF, nullptr, nullptr );
                // Reflect updated positions/forces into host buffers for system 0
                W.pack_uff_system( 0, W.ffu, /*bParams=*/false, /*bForces=*/true, /*bVel=*/false, /*blvec=*/false );
                break;
            }

            // GPU paths (evaluate UFF via OpenCL kernels)
            default: {
                // Sync current UFF state into host-side GPU arrays, then upload
                W.pack_uff_system( 0, W.ffu, /*bParams=*/true, /*bForces=*/false, /*bVel=*/false, /*blvec=*/true );
                // Ensure data are uploaded (safe to call multiple times)
                W.upload( /*bParams=*/true, /*bForces=*/false, /*bVel=*/false, /*blvec=*/true );
                double Etot = W.eval_UFF_ocl( nstepMax );
                // Download forces and positions back to host for inspection/use
                W.download( /*bForces=*/true, /*bVel=*/false );
                if(outE){ outE[0] = Etot; }
                nitrdione = nstepMax; // we executed nstepMax evaluation steps on GPU
                break;
            }
        }
    }else{
        switch(iParalel){
            case -1: nitrdione = W.run_multi_serial( nstepMax, Fconv, 1000.0, 1000 ); break; 
            case  0:
            case  1: nitrdione = W.run_omp_ocl     ( nstepMax, Fconv, 1000.0, 1000 ); break; 
            case  2: nitrdione = W.run_ocl_opt( nstepMax, Fconv    ); break; 
            case  3: nitrdione = W.run_ocl_loc( nstepMax, Fconv, 1 ); break; 
            case  4: nitrdione = W.run_ocl_loc( nstepMax, Fconv, 2 ); break; 
        }
    }
    return nitrdione;
    //return W.rum_omp_ocl( nstepMax, dt, Fconv, 1000.0, 1000 ); 
    //return W.run(nstepMax,dt,Fconv,ialg,outE,outF);  
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
    if(W.bUFF) W.pack_uff_system( isys, W.ffu, bParams, bForces, bVel, blvec ); else W.pack_system( isys, W.ffls[isys], bParams, bForces, bVel, blvec );
}
void unpack_system( int isys, bool bForces=false, bool bVel=false ){
    if(W.bUFF) W.unpack_uff_system( isys, W.ffu, bForces, bVel ); else W.unpack_system( isys, W.ffls[isys], bForces, bVel );
}

void upload_sys  ( int isys, bool bParams, bool bForces, bool bVel, bool blvec ){
    if(W.bUFF) W.upload_uff_sys( isys, bParams, bForces, bVel, blvec ); else W.upload_sys( isys, bParams, bForces, bVel, blvec );
}
void download_sys( int isys, bool bForces, bool bVel ){
    if(W.bUFF) W.download_uff_sys( isys, bForces, bVel ); else W.download_mmff_sys( isys, bForces, bVel );
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

} // extern "C"
