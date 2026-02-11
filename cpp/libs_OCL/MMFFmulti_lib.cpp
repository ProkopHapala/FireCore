
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
    printf( "init_buffers() \n" );
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

void init_buffers_UFF(){
    printf( "init_buffers_UFF() \n" );
    // Common buffers
    buffers.insert( { "apos",   (double*)W.nbmol.apos  } );
    buffers.insert( { "fapos",  (double*)W.nbmol.fapos } );
    buffers.insert( { "REQs",   (double*)W.nbmol.REQs  } );
    //buffers.insert( { "PLQs",   (double*)W.nbmol.PLQs  } );
    // UFF-specific buffers
    if(W.bUFF){ // UFF-specific buffers
        buffers.insert(  { "hneigh",    (double*)W.ffu.hneigh } );
        // GPU-computed hneigh (float4) downloaded from OpenCL into MolWorld_sp3_multi::hneigh_ocl
        if(W.hneigh_ocl) fbuffers.insert( { "gpu_hneigh", (float*)W.hneigh_ocl } );
        buffers.insert(  { "fint",      (double*)W.ffu.fint   } );
        buffers.insert(  { "bonParams", (double*)W.ffu.bonParams } );
        buffers.insert(  { "angParams", (double*)W.ffu.angParams } );
        buffers.insert(  { "dihParams", (double*)W.ffu.dihParams } );
        buffers.insert(  { "invParams", (double*)W.ffu.invParams } );

        ibuffers.insert( { "neighs",    (int*)W.ffu.neighs    } );
        ibuffers.insert( { "neighBs",   (int*)W.ffu.neighBs   } );
        ibuffers.insert( { "bonAtoms",  (int*)W.ffu.bonAtoms  } );
        ibuffers.insert( { "angAtoms",  (int*)W.ffu.angAtoms  } );
        ibuffers.insert( { "dihAtoms",  (int*)W.ffu.dihAtoms  } );
        ibuffers.insert( { "invAtoms",  (int*)W.ffu.invAtoms  } );

        // neighbor indices for angles, dihedrals, inversions
        ibuffers.insert( { "angNgs",    (int*)W.ffu.angNgs    } );
        ibuffers.insert( { "dihNgs",    (int*)W.ffu.dihNgs    } );
        ibuffers.insert( { "invNgs",    (int*)W.ffu.invNgs    } );

        // Atom-to-force-piece assembly map packed for OpenCL (MolWorld_sp3_multi::pack_uff_system)
        // NOTE: host_a2f_indices has length nA2F == ffu.nf; only first sum(a2f_counts) entries are meaningful,
        // remaining entries should be -1.
        if(W.host_a2f_offsets.size()>0) ibuffers.insert( { "a2f_offsets", (int*)W.host_a2f_offsets.data() } );
        if(W.host_a2f_counts .size()>0) ibuffers.insert( { "a2f_counts",  (int*)W.host_a2f_counts .data() } );
        if(W.host_a2f_indices.size()>0) ibuffers.insert( { "a2f_indices", (int*)W.host_a2f_indices.data() } );


        // ---- TODO: GPU buffers for UFF

    }
    // UFF-specific dimensions
    if(W.bUFF){
        ibuffers.insert( { "ndims", &W.ffu._natoms } ); //
        buffers.insert ( { "Es",    &W.ffu.Etot    } );
    }
    ibuffers.insert( { "selection", W.manipulation_sel  } );
    bbuffers.insert( { "ffflags",  &W.doBonded  } );
    // int _natoms, nbonds, nangles, ndihedrals, ninversions, nf; // 5
    // int i0dih,i0inv,i0ang,i0bon;                               // 4
    // double Etot, Eb, Ea, Ed, Ei;                               // 5
    printf( "MMFF_lib.cpp::init_buffers_UFF() ndims{natoms=%i, nbonds=%i, nangles=%i, ndihedrals=%i, ninversions=%i, nf=%i, i0dih=%i,i0inv=%i,i0ang=%i,i0bon=%i, }\n", W.ffu._natoms, W.ffu.nbonds, W.ffu.nangles, W.ffu.ndihedrals, W.ffu.ninversions, W.ffu.nf, W.ffu.i0dih, W.ffu.i0inv, W.ffu.i0ang, W.ffu.i0bon );
    printf( "MMFF_lib.cpp::init_buffers_UFF() Es{ Etot=%f, Eb=%f, Ea=%f, Ed=%f, Ei=%f, }\n", W.ffu.Etot, W.ffu.Eb, W.ffu.Ea, W.ffu.Ed, W.ffu.Ei );
}

void print_debugs( bool bParams, bool bNeighs, bool bShifts, bool bAtoms ){
    printf("print_debugs() W.bUFF=%i, bParams=%i, bNeighs=%i, bShifts=%i, bAtoms=%i \n", W.bUFF, bParams, bNeighs, bShifts, bAtoms);
    if(W.bUFF){
        W.ffu.printSizes();
        if(bParams) W.ffu.printAllParams(true, true, true, true, true);
        if(bAtoms ) W.ffu.print();
    } else {
        W.ffl.printSizes();
        if( bParams ) W.ffl.printAtomParams();
        if( bNeighs ) W.ffl.printNeighs();
        if( bShifts ) W.ffl.print_pbc_shifts();
        if( bAtoms  ) W.ffl.print();
    }
}

void print_setup(){
    if(W.bUFF){
        W.ffu.printSimulationSetup();
    }else{
        printf("MMFFmulti_lib::print_setup() bUFF=false\n");
       //W.ffl.printSimulationSetup();
    }
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
    // unbuffered printf()
    setbuf(stdout, NULL);
    setbuf(stderr, NULL);
    W.params.init( sElementTypes, sAtomTypes, sBondTypes, sAngleTypes, sDihedralTypes );
	W.builder.bindParams(&W.params);
    W.nSystems=nSys;
    bool bGrid = gridStep>0;
    W.bGridFF = bGrid;
    W.init();
    //init_buffers();
    return &W;
}

void setSwitches2( int CheckInvariants, int PBC, int NonBonded, int NonBondNeighs,  int SurfAtoms, int GridFF, int MMFF, int Angles, int PiSigma, int PiPiI ){
    #define _setbool(b,i) { if(i>0){b=true;}else if(i<0){b=false;} }
    _setbool( W.bCheckInvariants, CheckInvariants  );
    _setbool( W.bPBC           , PBC       );

    _setbool( W.bNonBonded     , NonBonded );
    _setbool( W.bNonBondNeighs , NonBondNeighs );

    _setbool( W.bSurfAtoms   , SurfAtoms );
    _setbool( W.bGridFF      , GridFF    );

    _setbool( W.bMMFF        , MMFF      );
    _setbool( W.ffl.doAngles , Angles    );
    _setbool( W.ffl.doPiSigma, PiSigma   );
    _setbool( W.ffl.doPiPiI  , PiPiI     );

    // Keep subtraction/clamp consistent with nonbond evaluation.
    // NOTE: switch semantics are 0=keep, >0=force true, <0=force false.
    if(NonBonded!=0){
        bool b = (NonBonded>0);
        W.ffl.bSubtractBondNonBond  = b;
        W.ffl.bSubtractAngleNonBond = b;
        W.ffl.bClampNonBonded       = b;
    }

    printf( "setSwitches2() W.bCheckInvariants==%i bPBC=%i | bNonBonded=%i bNonBondNeighs=%i | bSurfAtoms=%i bGridFF=%i | bMMFF=%i doAngles=%i doPiSigma=%i doPiPiI=%i | subNB(bond=%i,angle=%i,clamp=%i) \n",
        W.bCheckInvariants, W.bPBC,  W.bNonBonded, W.bNonBondNeighs, W.bSurfAtoms, W.bGridFF, W.bMMFF, W.ffl.doAngles, W.ffl.doPiSigma, W.ffl.doPiPiI,
        W.ffl.bSubtractBondNonBond, W.ffl.bSubtractAngleNonBond, W.ffl.bClampNonBonded
    );

    #undef _setbool
}

double eval_force_MMFF( int iParalel ){
    if(W.bUFF){ printf("eval_force_MMFF(): ERROR W.bUFF=true; this is MMFF-only\n"); return 0.0; }
    const int isys = 0;
    if(iParalel==2){
        int err=0;
        W.bOcl=true;
        if( W.task_MMFF==0 ){ W.setup_MMFFf4_ocl(); }
        // Keep GPU host-side flat buffers in sync with current CPU state
        // blvec=false: lvec already uploaded during init; avoids 0-size pbcshifts upload when npbc=0
        W.pack_system  ( isys, W.ffls[isys], true, false, false, false );
        W.upload_sys   ( isys, true, false, false, false );
        // Set bSubtractVdW consistent with NonBonded switch (same as run_ocl_opt)
        W.ocl.bSubtractVdW = W.bNonBonded ? (!W.bExclusion2) : 0;
        if(W.task_MMFF){
            int err_set = clSetKernelArg( W.ocl.kernels[W.task_MMFF->ikernel], 16, sizeof(int), &W.ocl.bSubtractVdW );
            OCL_checkError(err_set, "eval_force_MMFF.set_bSubtractVdW");
        }
        // Force evaluation: getMMFFf4 computes per-atom forces and stores recoil in fneigh;
        // updateAtomsMMFFf4 assembles recoil via bkNeighs. With dt=0 it skips dynamics.
        Quat4f MDpars_save = W.MDpars[isys];
        W.MDpars[isys] = Quat4f{0.0f, 0.0f, 0.0f, 0.0f};  // dt=0 → assembly only, no move
        err |= W.ocl.upload( W.ocl.ibuff_MDpars, W.MDpars );
        err |= W.task_cleanF->enque_raw();
        err |= W.task_MMFF  ->enque_raw();
        err |= W.task_move  ->enque_raw();  // assembles recoil forces, skips dynamics (dt=0)
        err |= W.ocl.finishRaw();
        OCL_checkError(err, "eval_force_MMFF(GPU).finishRaw");
        W.MDpars[isys] = MDpars_save;  // restore original MDpars
        // Download assembled forces (aforces now includes recoil)
        W.download_mmff_sys( isys, true, false );
        W.unpack_system( isys, W.ffls[isys], true, false );
        // Etot is not guaranteed to be updated on GPU path; compute explicitly from aforces.w (atoms only)
        {
            const int i0v = isys * W.ocl.nvecs;
            double Etot = 0.0;
            for(int i=0; i<W.ffl.natoms; i++){
                Etot += (double)W.aforces[i0v+i].w;
            }
            W.ffls[isys].Etot = Etot;
        }
        // Copy forces back to shared nbmol arrays so Python can read them
        for(int i=0; i<W.ffl.natoms; i++){ W.nbmol.fapos[i] = W.ffls[isys].fapos[i]; }
        return W.ffls[isys].Etot;
    }else{
        W.bOcl=false;
        W.setNonBond( W.bNonBonded );
        // Propagate switches from template ffl to per-system copy ffls[isys]
        W.ffls[isys].doBonds              = W.ffl.doBonds;
        W.ffls[isys].doAngles             = W.ffl.doAngles;
        W.ffls[isys].doPiSigma            = W.ffl.doPiSigma;
        W.ffls[isys].doPiPiI              = W.ffl.doPiPiI;
        W.ffls[isys].bSubtractBondNonBond = W.ffl.bSubtractBondNonBond;
        W.ffls[isys].bSubtractAngleNonBond= W.ffl.bSubtractAngleNonBond;
        W.ffls[isys].bClampNonBonded      = W.ffl.bClampNonBonded;
        printf("eval_force_MMFF(CPU) doAngles=%i doPiSigma=%i doPiPiI=%i subBondNB=%i subAngNB=%i clampNB=%i\n",
            W.ffls[isys].doAngles, W.ffls[isys].doPiSigma, W.ffls[isys].doPiPiI,
            W.ffls[isys].bSubtractBondNonBond, W.ffls[isys].bSubtractAngleNonBond, W.ffls[isys].bClampNonBonded);
        // Force-only evaluation on CPU
        W.ffls[isys].cleanForce();
        double E = W.ffls[isys].eval(true);
        // Copy forces back to shared nbmol arrays so Python can read them
        for(int i=0; i<W.ffl.natoms; i++){ W.nbmol.fapos[i] = W.ffls[isys].fapos[i]; }
        return E;
    }
}


void setSwitchesUFF( int DoBond, int DoAngle, int DoDihedral, int DoInversion, int DoAssemble, int SubtractBondNonBond, int ClampNonBonded ){
    #define _setbool(b,i) { if(i>0){b=true;}else if(i<0){b=false;} }
    if(W.uff_ocl){
        _setbool( W.uff_ocl->bUFF_bonds,      DoBond );
        _setbool( W.uff_ocl->bUFF_angles,     DoAngle );
        _setbool( W.uff_ocl->bUFF_dihedrals,  DoDihedral );
        _setbool( W.uff_ocl->bUFF_inversions, DoInversion );
        _setbool( W.uff_ocl->bUFF_assemble,   DoAssemble );
        _setbool( W.uff_ocl->bSubtractNB,     SubtractBondNonBond );
        _setbool( W.uff_ocl->bClampNonBonded, ClampNonBonded );
        // If any higher-order interactions are enabled, force assembly on GPU
        // Angle/dihedral/inversion kernels are interaction-centric and do not safely accumulate to fapos directly.
        // if( (W.uff_ocl->bUFF_angles || W.uff_ocl->bUFF_dihedrals || W.uff_ocl->bUFF_inversions) ){
        //     W.uff_ocl->bUFF_assemble = true;
        // }
    }
    _setbool( W.ffu.bDoBond,              DoBond );
    _setbool( W.ffu.bDoAngle,             DoAngle );
    _setbool( W.ffu.bDoDihedral,          DoDihedral );
    _setbool( W.ffu.bDoInversion,         DoInversion );
    _setbool( W.ffu.bDoAssemble,          DoAssemble );
    _setbool( W.ffu.bSubtractBondNonBond, SubtractBondNonBond );
    _setbool( W.ffu.bClampNonBonded,      ClampNonBonded );
    // Sync UFF non-bond master switch from world; ForceField::setNonBondStrategy() depends on this
    W.ffu.bNonBonded = W.bNonBonded;
    // Enforce a consistent non-bond strategy on CPU to avoid UFF::run() flipping flags implicitly
    // Map: SubtractBondNonBond==1 -> imode<0 (no neighbor culling, subtract+clamp)
    //      SubtractBondNonBond==0 -> imode>0 (neighbor-based, no subtract, no clamp)
    W.ffu.setNonBondStrategy( !W.ffu.bSubtractBondNonBond );
    // Propagate the final CPU flags to GPU so kernels match CPU behavior
    W.uff_ocl->bSubtractNB        = W.ffu.bSubtractBondNonBond && W.ffu.bNonBonded;
    W.uff_ocl->bSubtractNB_angle  = W.ffu.bSubtractAngleNonBond && W.ffu.bNonBonded;
    W.uff_ocl->bClampNonBonded    = W.ffu.bClampNonBonded    && W.ffu.bNonBonded;
    printf( "setSwitchesUFF() DoBond=%i DoAngle=%i DoDihedral=%i DoInversion=%i DoAssemble=%i SubtractBondNonBond=%i ClampNonBonded=%i \n", W.ffu.bDoBond, W.ffu.bDoAngle, W.ffu.bDoDihedral, W.ffu.bDoInversion, W.ffu.bDoAssemble, W.ffu.bSubtractBondNonBond, W.ffu.bClampNonBonded );
    if(W.uff_ocl){
        printf("GPU flags: bUFF_bonds=%i bUFF_angles=%i bUFF_dihedrals=%i bUFF_inversions=%i bUFF_assemble=%i bSubtractNB=%i bSubtractNB_angle=%i bClampNB=%i\n",
            W.uff_ocl->bUFF_bonds, W.uff_ocl->bUFF_angles, W.uff_ocl->bUFF_dihedrals, W.uff_ocl->bUFF_inversions, W.uff_ocl->bUFF_assemble,
            W.uff_ocl->bSubtractNB, W.uff_ocl->bSubtractNB_angle, W.uff_ocl->bClampNonBonded );
        // Re-bind kernel args with the updated flags and UFF parameters to keep GPU in sync with CPU
        // Use current UFF parameters from W.ffu
        float Rdamp = (float)W.ffu.Rdamp;
        float FmaxNB = (float)W.ffu.FmaxNonBonded;
        float SubNBTorsion = (float)W.ffu.SubNBTorsionFactor;
        W.uff_ocl->setup_kernels( Rdamp, FmaxNB, SubNBTorsion );
    }
    #undef _setbool
}

//int run( int nstepMax, double dt, double Fconv, int ialg, double* outE, double* outF, int iParalel ){
int run( int nstepMax, double dt, double Fconv, int ialg, double damping, double* outE, double* outF, double* outV, double* outVF, int iParalel ){
    W.bOcl= iParalel > 0;
    Mat3d lvec = W.ffls[0].lvec;
    printf( "MMFFmulti_lib.cpp run() iParalel=%i W.bUFF=%i ialg=%i dt=%g Fconv=%g nstepMax=%i  outE=%p outF=%p \n", iParalel, W.bUFF, ialg, dt, Fconv, nstepMax, outE, outF );
    //printf( "run Fconv=%g lvec{{%6.3f,%6.3f,%6.3f}{%6.3f,%6.3f,%6.3f}{%6.3f,%6.3f,%6.3f}}\n", Fconv, lvec.a.x,lvec.a.y,lvec.a.z, lvec.b.x,lvec.b.y,lvec.b.z, lvec.c.x,lvec.c.y,lvec.c.z );
    int nitrdione=0;
    if(W.bUFF){
        // UFF branch: CPU (serial/OpenMP) vs GPU (OpenCL) based on iParalel
        switch(iParalel){
            // CPU paths (use UFF::run on host)
            case -1:
            case  0: { nitrdione = W.ffu.run    ( nstepMax, dt, Fconv, 1000.0, damping, outE, outF, outV, outVF ); break; }
            case 1:  { nitrdione = W.ffu.run_omp( nstepMax, dt, Fconv, 1000.0, damping, outE, outF, outV, outVF ); break; }
            case 2:  { // GPU paths (evaluate UFF via OpenCL kernels)
                //                           bParams bForces  bVel    bLvec
                W.pack_uff_system( 0, W.ffu, true,   false,  false,   true   );
                W.upload_uff     (           true,   false,  false,   true   );
                double Etot = W.eval_UFF_ocl( nstepMax );
                W.download_uff( true, false );            // fills W.atoms and W.aforces from GPU
                W.unpack_uff_system( 0, W.ffu, true, false ); // propagate W.aforces -> W.ffu.fapos (and positions)
                if(outE){ outE[0] = Etot; }
                nitrdione = nstepMax;
                break;
            }
            default: printf("run() iParalel=%i not implemented\n", iParalel); break;
        }
    }else{
        switch(iParalel){
            case -1: nitrdione = W.run_multi_serial( nstepMax, Fconv, 1000.0, 1000 ); break;
            case  0:
            case  1: nitrdione = W.run_omp_ocl( nstepMax, Fconv, 1000.0, 1000 ); break;
            case  2: nitrdione = W.run_ocl_opt( nstepMax, Fconv    ); break;
            case  3: nitrdione = W.run_ocl_loc( nstepMax, Fconv, 1 ); break;
            case  4: nitrdione = W.run_ocl_loc( nstepMax, Fconv, 2 ); break;
            default: printf("run() iParalel=%i not implemented\n", iParalel); break;
        }
    }
    return nitrdione;
    //return W.rum_omp_ocl( nstepMax, dt, Fconv, 1000.0, 1000 );
    //return W.run(nstepMax,dt,Fconv,ialg,outE,outF);
}

// Evaluate forces for multiple configurations
// confs: [nConf, natoms, 3] (row-major doubles)
// outF : [nConf, natoms, 3]
// iParalel: CPU serial/OMP uses same mapping as run(); GPU==2 uses OpenCL UFF batched over W.nSystems replicas
int scan( int nConf, double* confs_, double* outF_, int iParalel ){
    if(!W.bUFF){ printf("scan(): ERROR bUFF=false; only UFF supported for now\n"); return 0; }
    const int natoms = W.ffu._natoms;
    auto confs = (Vec3d*)confs_;
    auto outF  = (Vec3d*)outF_;
    int nDone = 0;
    const int dbg_sys = 0;

    if( (iParalel==2) && W.uff_ocl ){
        if(!W.uff_ocl->bKernelPrepared){ W.uff_ocl->setup_kernels( (float)W.ffu.Rdamp, (float)W.ffu.FmaxNonBonded, (float)W.ffu.SubNBTorsionFactor ); }
        for(int isys=0; isys<W.nSystems; ++isys){ W.pack_uff_system( isys, W.ffu, true, false, false, true ); }
        W.upload_uff( true, false, false, true );

        for(int ib=0; ib<nConf; ){
            int nBatch = W.nSystems; if(ib+nBatch>nConf) nBatch = nConf-ib;
            for(int i=0;i<nBatch;i++){
                int isys = i;
                idebug = (isys==dbg_sys) ? 4 : 0;
                Vec3d* apos = W.ffu.apos;
                for(int ia=0; ia<natoms; ia++){ apos[ia] = confs[(ib+i)*natoms + ia]; }
                W.pack_uff_system( isys, W.ffu, false, false, false, false );
            }
            W.upload_uff( false, false, false, false );
            W.eval_UFF_ocl( 1 );
            W.download_uff( true, false );
            for(int i=0;i<nBatch;i++){
                int isys = i;
                W.unpack_uff_system( isys, W.ffu, true, false );
                Vec3d* fapos = W.ffu.fapos;
                for(int ia=0; ia<natoms; ia++){ outF[(ib+i)*natoms + ia] = fapos[ia]; }
            }
            ib += nBatch; nDone += nBatch;
        }
    }else{
        for(int ic=0; ic<nConf; ic++){
            idebug = (ic==dbg_sys) ? 4 : 0;
            for(int ia=0; ia<natoms; ia++){ W.ffu.apos[ia] = confs[ic*natoms + ia]; }
            for(int ia=0; ia<natoms; ia++){ W.ffu.fapos[ia] = Vec3dZero; }
            if(iParalel==1){ W.ffu.run_omp( 1, 0.0, 1e-6, 1000.0, 0.1, nullptr, nullptr, nullptr, nullptr ); }
            else            { W.ffu.run    ( 1, 0.0, 1e-6, 1000.0, 0.1, nullptr, nullptr, nullptr, nullptr ); }
            for(int ia=0; ia<natoms; ia++){ outF[ic*natoms + ia] = W.ffu.fapos[ia]; }
            nDone++;
        }
    }
    return nDone;
}

int scan_relaxed( int nConf, double* confs_, double* outF_, int niter, double dt, double damping, double Fconv, double Flim, int iParalel ){
    if(!W.bUFF){ printf("scan_relaxed(): ERROR bUFF=false; only UFF supported\n"); return 0; }
    const int natoms = W.ffu._natoms;
    auto confs = (Vec3d*)confs_;
    auto outF  = (Vec3d*)outF_;
    int nDone = 0;
    const int dbg_sys = 0;
    if( (iParalel==2) && W.uff_ocl ){
        W.setup_UFF_ocl();
        if(!W.uff_ocl->bKernelPrepared){ W.uff_ocl->setup_kernels( (float)W.ffu.Rdamp, (float)W.ffu.FmaxNonBonded, (float)W.ffu.SubNBTorsionFactor ); }
        for(int isys=0; isys<W.nSystems; ++isys){ W.pack_uff_system( isys, W.ffu, true, false, false, true ); }
        W.upload_uff( true, false, false, true );
        setTrjName("scan_relaxed_gpu", 1, (int*)&W.nPBC);
        for(int ib=0; ib<nConf; ){
            int nBatch = W.nSystems; if(ib+nBatch>nConf) nBatch = nConf-ib;
            for(int i=0;i<nBatch;i++){
                int isys = i;
                idebug = (isys==dbg_sys) ? 4 : 0;
                for(int ia=0; ia<natoms; ia++){ W.ffu.apos[ia] = confs[(ib+i)*natoms + ia]; }
                W.pack_uff_system( isys, W.ffu, false, false, false, false );
            }
            W.upload_uff( false, false, false, false );
            // Zero velocities on host and upload to GPU to ensure a clean start for each batch
            memset( W.avel, 0, W.uff_ocl->nAtomsTot * sizeof(Quat4f) );
            W.uff_ocl->upload( W.uff_ocl->ibuff_avel, (float*)W.avel );

            W.run_uff_ocl( niter, dt, damping, Fconv, Flim );

            W.download_uff( true, true ); // Download final forces and positions
            for(int i=0;i<nBatch;i++){
                int isys = i;
                W.unpack_uff_system( isys, W.ffu, true, false ); // Unpack final forces
                Vec3d* fapos = W.ffu.fapos;
                for(int ia=0; ia<natoms; ia++){ outF[(ib+i)*natoms + ia] = fapos[ia]; }
                // Trajectory saving is now handled inside run_uff_ocl
            }
            ib += nBatch; nDone += nBatch;
        }
    }else{
        for(int ic=0; ic<nConf; ic++){
            idebug = (ic==dbg_sys) ? 4 : 0;
            for(int ia=0; ia<natoms; ia++){ W.ffu.apos[ia] = confs[ic*natoms + ia]; }
            W.ffu.cleanVelocity(); // Zero velocities to ensure a clean start for each configuration
            char fname[256]; sprintf(fname, "scan_relaxed_cpu_%03i.xyz", ic);
            setTrjName(fname, 1, (int*)&W.nPBC);
            if(iParalel==1){ W.ffu.run_omp( niter, dt, Fconv, Flim, damping, nullptr, nullptr, nullptr, nullptr ); } // fapos is cleaned inside run()
            else           { W.ffu.run    ( niter, dt, Fconv, Flim, damping, nullptr, nullptr, nullptr, nullptr ); } // fapos is cleaned inside run()
            setTrjName(nullptr, -1, (int*)&W.nPBC ); // Disable trajectory saving
            for(int ia=0; ia<natoms; ia++){ outF[ic*natoms + ia] = W.ffu.fapos[ia]; }
            nDone++;
        }
    }
    return nDone;
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
