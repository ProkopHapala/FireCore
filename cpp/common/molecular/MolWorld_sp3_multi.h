
#ifndef MolWorld_sp3_ocl_h
#define MolWorld_sp3_ocl_h

#include "MolWorld_sp3.h"
//#include "OCL_DFT.h"
//#include "OCL_PP.h"
#include "OCL_MM.h"
#include "datatypes_utils.h"

#include "Confs.h"

// ======================================
// class:        MolWorld_sp3_ocl
// ======================================

class MolWorld_sp3_multi : public MolWorld_sp3 { public:
    OCL_MM     ocl;

    int nSystems    = 1;
    //int iSystemCur  = 0;    // currently selected system replica
    int iSystemCur  = 5;    // currently selected system replica
    //int iSystemCur  = 8;    // currently selected system replica
    bool bGPU_MMFF = true;

    Quat4f* atoms      =0;
    Quat4f* aforces    =0;
    Quat4f* avel       =0;

    Quat4f* constr     =0;

    Quat4i* neighs     =0;
    Quat4i* neighCell  =0;
    Quat4i* bkNeighs   =0;
    //Quat4f* neighForce =0;

    Quat4f* REQs       =0;
    Quat4f* MMpars     =0;
    Quat4f* BLs        =0;
    Quat4f* BKs        =0;
    Quat4f* Ksp        =0;
    Quat4f* Kpp        =0;

    cl_Mat3*  lvecs    =0;
    cl_Mat3* ilvecs    =0;

    OCLtask* task_cleanF=0;
    OCLtask* task_NBFF=0;
    OCLtask* task_MMFF=0;
    OCLtask* task_move=0;
    OCLtask* task_print=0;

// ======== Functions

void realloc( int nSystems_ ){
    printf("MolWorld_sp3_multi::realloc() \n");
    nSystems=nSystems_;
    printf( "MolWorld_sp3_multi::realloc() Systems %i nAtoms %i nnode %i \n", nSystems, ffl.natoms,  ffl.nnode );
    ocl.initAtomsForces( nSystems, ffl.natoms,  ffl.nnode );
    //printf( "MolWorld_sp3_multi::realloc() Systems %i nAtoms %i nnode %i nvecs %i \n", nSystems, ocl.nAtoms, ocl.nnode, ocl.nvecs );
    // --- dynamical
    _realloc( atoms,     ocl.nvecs*nSystems  );
    _realloc( aforces,   ocl.nvecs*nSystems  );
    //_realloc( avel,    ocl.nvecs*nSystems  );
    _realloc( constr,    ocl.nAtoms*nSystems  );
    // --- params
    _realloc( neighs,    ocl.nAtoms*nSystems );
    _realloc( neighCell, ocl.nAtoms*nSystems );
    //_realloc( bkNeighs,  ocl.nAtoms*nSystems );
    _realloc( bkNeighs,  ocl.nvecs*nSystems );
    _realloc( atoms,     ocl.nvecs*nSystems  );
    _realloc( atoms,     ocl.nvecs*nSystems  );
    _realloc( atoms,     ocl.nvecs*nSystems  );
    _realloc( REQs,      ocl.nAtoms*nSystems );
    _realloc( MMpars,    ocl.nnode*nSystems  );
    _realloc( BLs,       ocl.nnode*nSystems  );
    _realloc( BKs,       ocl.nnode*nSystems  );
    _realloc( Ksp,       ocl.nnode*nSystems  );
    _realloc( Kpp,       ocl.nnode*nSystems  );

    _realloc( lvecs,     nSystems  );
    _realloc( ilvecs,    nSystems  );

    // ToDo : it may be good to bind buffer directly in p_cpu buffer inside   OCLsystem::newBuffer()
}

void pack_system( int isys, MMFFsp3_loc& ff, bool bParams=0, bool bForces=0 , float l_rnd=-1 ){
    printf("MolWorld_sp3_multi::pack_system(%i) \n", isys);
    //ocl.nvecs;
    int i0n = isys * ocl.nnode;
    int i0a = isys * ocl.nAtoms;
    int i0v = isys * ocl.nvecs;
    pack( ff.nvecs, ff.apos, atoms+i0v );
    for(int i=0; i<ocl.nAtoms; i++){ Quat4f a=atoms[i+i0v]; a.w=-1.0; constr[i+i0a] = a; }  // contrains
    /*
    if(l_rnd>0){ 
        printf("WARRNING: random noise added to apos \n");
        for(int i=0; i<ff.nvecs; i++){ 
            atoms[i+i0v].f.addRandomCube(l_rnd); 
            if(i>=ff.natoms){ atoms[i+i0v].f.normalize(); } // normalize pi-bonds 
        } 
        atoms[0+i0v].x=isys;
    }
    */
    if(bForces){
        pack( ff.nvecs, ff.fapos, aforces+i0v );
        //pack( f.nvecs, ff.avel, avel+i0v );
    }
    if(bParams){
        Mat3_to_cl( ff.   lvec,  lvecs[isys] );
        Mat3_to_cl( ff.invLvec, ilvecs[isys] );
        copy    ( ff.natoms, ff.aneighCell, neighCell+i0a );
        copy    ( ff.natoms, ff.aneighs,    neighs   +i0a );
        //copy_add( ff.natoms, ff.aneighs,    neighs   +i0a,           0      );
        copy_add( ff.natoms, ff.bkneighs,   bkNeighs +i0v,           i0n*8  );
        copy_add( ff.nnode , ff.bkneighs,   bkNeighs +i0v+ff.natoms, i0n*8 + 4*ff.nnode );
        pack    ( ff.nnode , ff.apars,      MMpars   +i0n );
        pack    ( ff.nnode , ff.bLs,        BLs      +i0n );
        pack    ( ff.nnode , ff.bKs,        BKs      +i0n );
        pack    ( ff.nnode , ff.Ksp,        Ksp      +i0n );
        pack    ( ff.nnode , ff.Kpp,        Kpp      +i0n );
        pack    ( nbmol.n  , nbmol.REQs, REQs        +i0a );
    }

    //if(isys==5)
    {
        //for(int i=0; i<ocl.nAtoms; i++){ Quat4i ng = neighs[i0a+i]; printf( "ngs[%4i] (%4i,%4i,%4i,%4i) \n", i, ng.x,ng.y,ng.z,ng.w ); }
        //for(int i=0; i<ocl.nvecs; i++){ Quat4f p = atoms[i0v+i]; printf( "apos[%4i] (%5.3f,%5.3f,%5.3f) pi %i \n", i, p.x,p.y,p.z, i>=ocl.nAtoms ); }
        //if(isys==0)for(int i=0; i<ocl.nvecs; i++){Quat4i bkng = bkNeighs[i0v+i]; printf( "bkng[%4i] (%4i,%4i,%4i,%4i) pi %i\n", i, bkng.x,bkng.y,bkng.z,bkng.w, i>=ocl.nAtoms ); }
    }
}

void unpack_system(  int isys, MMFFsp3_loc& ff, bool bForces=0 ){
    printf("MolWorld_sp3_multi::unpack_system(%i) \n", isys);
    int i0n = isys * ocl.nnode;
    int i0a = isys * ocl.nAtoms;
    int i0v = isys * ocl.nvecs;
    unpack( ff.nvecs, ff.apos, atoms+i0v );
    if(bForces){
        unpack( ff.nvecs, ff.fapos, aforces+i0v );
        //unpack( ff.nvecs, ff.fapos, avel+i0v );
    }
}

void upload(  bool bParams=0, bool bForces=0 ){
    printf("MolWorld_sp3_multi::upload() \n");
    ocl.upload( ocl.ibuff_atoms,  atoms  );
    ocl.upload( ocl.ibuff_constr, constr );
    if(bForces){
        ocl.upload( ocl.ibuff_aforces, aforces );
        //ocl.upload( ocl.ibuff_avel,    avel    );
    }
    if(bParams){
        ocl.upload( ocl.ibuff_lvecs,   lvecs );
        ocl.upload( ocl.ibuff_ilvecs,  lvecs );
        ocl.upload( ocl.ibuff_neighs,     neighs    );
        ocl.upload( ocl.ibuff_neighCell,  neighCell );
        ocl.upload( ocl.ibuff_bkNeighs,   bkNeighs  );
        //ocl.upload( ocl.ibuff_neighForce, neighForce  );
        ocl.upload( ocl.ibuff_REQs,   REQs   );
        ocl.upload( ocl.ibuff_MMpars, MMpars );
        ocl.upload( ocl.ibuff_BLs, BLs );
        ocl.upload( ocl.ibuff_BKs, BKs );
        ocl.upload( ocl.ibuff_Ksp, Ksp );
        ocl.upload( ocl.ibuff_Kpp, Kpp ); 
    }
}

void download( bool bForces=0  ){
    printf("MolWorld_sp3_multi::download() \n");
    ocl.download( ocl.ibuff_atoms, atoms );
    if(bForces){
        ocl.download( ocl.ibuff_aforces, aforces );
        //ocl.download( ocl.ibuff_avel,    avel    );
    }
}

virtual void init( bool bGrid ) override {
    printf("# ========== MolWorld_sp3_multi::init() START\n");
    ocl.print_devices(true);
    ocl.init();
    ocl.makeKrenels_MM("common_resources/cl" );
    MolWorld_sp3::init(bGrid);

    // ----- init systems
    realloc( nSystems );
    float random_init = 0.5;
    for(int i=0; i<nSystems; i++){
        pack_system( i, ffl, true, false, random_init );
    }
    upload( true, false );
    bGridFF=false;
    bOcl   =false;
    setup_MMFFf4_ocl();
    printf("# ========== MolWorld_sp3_multi::init() DONE\n");
}

// ======================== SETUP OCL KERNEL functions

void setup_MMFFf4_ocl(){
    printf("MolWorld_sp3_multi::setup_MMFFf4_ocl() \n");
    ocl.nDOFs.x=ff.natoms;
    ocl.nDOFs.y=ff.nnode;
    if(!task_move  )task_move   = ocl.setup_updateAtomsMMFFf4( ff4.natoms, ff4.nnode       );
    if(!task_print )task_print  = ocl.setup_printOnGPU       ( ff4.natoms, ff4.nnode       );    /// Print on GPU 
    if(!task_MMFF  )task_MMFF   = ocl.setup_getMMFFf4        ( ff4.natoms, ff4.nnode, bPBC );
    if(!task_NBFF  )task_NBFF   = ocl.setup_getNonBond       ( ff4.natoms, ff4.nnode, nPBC, gridFF.Rdamp  );
    if(!task_cleanF)task_cleanF = ocl.setup_cleanForceMMFFf4 ( ff4.natoms, ff4.nnode       );
}

void setup_NBFF_ocl(){
    printf("MolWorld_sp3_multi::setup_NBFF_ocl() \n");
    ocl.nDOFs.x=ff.natoms;
    ocl.nDOFs.y=ff.nnode;
    if(!task_cleanF )task_cleanF = ocl.setup_cleanForceMMFFf4 ( ff4.natoms, ff4.nnode        );
    if(!task_move   )task_move   = ocl.setup_updateAtomsMMFFf4( ff4.natoms, ff4.nnode        ); 
    if(!task_print  )task_print  = ocl.setup_printOnGPU       ( ff4.natoms, ff4.nnode        );    /// Print on GPU 
    if(!task_NBFF   )task_NBFF   = ocl.setup_getNonBond       ( ff4.natoms, ff4.nnode, nPBC, gridFF.Rdamp  );
}


void evaluateSubstrate(){

}



// ======================== EVAL OCL KERNEL functions

double eval_MMFFf4_ocl( int niter ){ 
    printf("MolWorld_sp3_multi::eval_MMFFf4_ocl() \n");
    int err=0;
    if( task_MMFF==0 )setup_MMFFf4_ocl();
    
    {  // DEBUG --- evaluate on CPU
        //bool bEval_ffl = false;
        bool bEval_ffl = true;
        unpack_system( iSystemCur, ffl );
        ffl.Rdamp = gridFF.Rdamp;
        ffl.cleanForce();
        ffl.eval();    //for(int i=0; i<ffl.nvecs; i++){  printf("ffl[%4i] f(%10.5f,%10.5f,%10.5f) p(%10.5f,%10.5f,%10.5f) pi %i \n", i,  ffl.fapos[i].x,ffl.fapos[i].y,ffl.fapos[i].z,   ffl.apos[i].x,ffl.apos[i].y,ffl.apos[i].z,  i>=ffl.natoms ); }
        nbmol.evalLJQs_ng4_PBC( ffl.aneighs, ffl.aneighCell, ffl.lvec, ffl.nPBC, gridFF.Rdamp );
        fcog  = sum ( ffl.natoms, ffl.fapos   );
        tqcog = torq( ffl.natoms, ffl.apos, ffl.fapos );
        if(  fcog.norm2()>1e-8 ){ printf("WARRNING: ffl.eval_MMFFf4_ocl() CPU |fcog| =%g; fcog=(%g,%g,%g) bEval_ffl %i \n", bEval_ffl, fcog.norm(),  fcog.x, fcog.y, fcog.z, bEval_ffl ); exit(0); }
    }
    
    /*
    { // DEBUG --- evaluate on CPU
        unpack_system( iSystemCur, ffl );
        ffl  .cleanForce();
        nbmol.evalLJQs_ng4_PBC( ffl.aneighs, ffl.aneighCell, ffl.lvec, ffl.nPBC, gridFF.Rdamp );
        fcog  = sum ( ffl.natoms, ffl.fapos   );
        tqcog = torq( ffl.natoms, ffl.apos, ffl.fapos );
        if(  fcog.norm2()>1e-8 ){ printf("WARRNING: ffl.eval_MMFFf4_ocl() CPU |fcog| =%g; fcog=(%g,%g,%g)\n", fcog.norm(),  fcog.x, fcog.y, fcog.z ); exit(0); }
    }
    */

    // evaluate on GPU
    for(int i=0; i<niter; i++){
        err |= task_cleanF->enque_raw();  // DEBUG: this should be solved inside  task_move->enque_raw();
        err |= task_MMFF  ->enque_raw();
        err |= task_NBFF  ->enque_raw();
        err |= task_print ->enque_raw(); // DEBUG: just printing the forces before assempling
        err |= task_move  ->enque_raw(); 
        OCL_checkError(err, "eval_MMFFf4_ocl_1");
    }
    //err |= ocl.finishRaw();
    //OCL_checkError(err, "eval_MMFFf4_ocl_2");

    //printf( "ocl.download(n=%i) \n", n );
    //ocl.download( ocl.ibuff_aforces, ff4.fapos, ff4.nvecs );
    //ocl.download( ocl.ibuff_atoms,   ff4.apos , ff4.nvecs );
    ocl.download( ocl.ibuff_aforces, ff4.fapos, ff4.nvecs, ff4.nvecs*iSystemCur );
    ocl.download( ocl.ibuff_atoms,   ff4.apos , ff4.nvecs, ff4.nvecs*iSystemCur );
    //for(int i=0; i<ff4.nvecs; i++){  printf("CPU[%i] p(%g,%g,%g) f(%g,%g,%g) pi %i \n", i, ff4.apos[i].x,ff4.apos[i].y,ff4.apos[i].z,  ff4.fapos[i].x,ff4.fapos[i].y,ff4.fapos[i].z, i>=ff4.natoms ); }
    err |= ocl.finishRaw();
    OCL_checkError(err, "eval_MMFFf4_ocl_2");

    for(int i=0; i<ff4.nvecs; i++){  printf("OCL[%4i] f(%10.5f,%10.5f,%10.5f) p(%10.5f,%10.5f,%10.5f) pi %i \n", i, ff4.fapos[i].x,ff4.fapos[i].y,ff4.fapos[i].z,  ff4.apos[i].x,ff4.apos[i].y,ff4.apos[i].z,  i>=ff4.natoms ); }

    //ffl.eval();   // We already computed this above
    bool ret=false;
    printf("### Compare ffl.fapos,  GPU.fapos  \n"); ret |= compareVecs( ff4.natoms, ffl.fapos,  ff4.fapos,  1e-4, true );
    printf("### Compare ffl.fpipos, GPU.fpipos \n"); ret |= compareVecs( ff4.nnode,  ffl.fpipos, ff4.fpipos, 1e-4, true ); 
    if(ret){ printf("ERROR: GPU.eval() and ffl.eval() produce different results => exit() \n"); exit(0); }else{ printf("CHECKED: GPU task_MMFF.eval() == CPU ffl.eval() \n"); }

    //printf("GPU AFTER assemble() \n"); ff4.printDEBUG( false,false );
    unpack( ff4.natoms, ffl.  apos, ff4.  apos );
    unpack( ff4.natoms, ffl. fapos, ff4. fapos );
    unpack( ff4.nnode,  ffl. pipos, ff4. pipos );
    unpack( ff4.nnode,  ffl.fpipos, ff4.fpipos );

    // ---- Check Invariatns
    fcog  = sum ( ffl.natoms, ffl.fapos   );
    tqcog = torq( ffl.natoms, ffl.apos, ffl.fapos );
    if(  fcog.norm2()>1e-8 ){ printf("WARRNING: eval_MMFFf4_ocl |fcog| =%g; fcog=(%g,%g,%g)\n", fcog.norm(),  fcog.x, fcog.y, fcog.z ); exit(0); }
    //if( tqcog.norm2()>1e-8 ){ printf("WARRNING: eval_MMFFf4_ocl |torq| =%g; torq=(%g,%g,%g)\n", tqcog.norm(),tqcog.x,tqcog.y,tqcog.z ); exit(0); }   // NOTE: torq is non-zero because pi-orbs have inertia
    
    exit(0);

    return 0;
}


double eval_NBFF_ocl( int niter ){ 
    printf("MolWorld_sp3_multi::eval_NBFF_ocl() \n");
    int err=0;
    if( task_NBFF==0 ){ setup_NBFF_ocl(); }

    { // DEBUG --- evaluate on CPU
        unpack_system( iSystemCur, ffl );
        ffl  .cleanForce();
        nbmol.evalLJQs_ng4_PBC( ffl.aneighs, ffl.aneighCell, ffl.lvec, ffl.nPBC, gridFF.Rdamp );
        fcog  = sum ( ffl.natoms, ffl.fapos   );
        tqcog = torq( ffl.natoms, ffl.apos, ffl.fapos );
        if(  fcog.norm2()>1e-8 ){ printf("WARRNING: eval_NBFF_ocl() CPU |fcog| =%g; fcog=(%g,%g,%g)\n", fcog.norm(),  fcog.x, fcog.y, fcog.z ); exit(0); }
    }
    // evaluate on GPU
    for(int i=0; i<niter; i++){
        err |= task_cleanF->enque_raw(); // DEBUG: this should be solved inside  task_move->enque_raw();
        err |= task_NBFF  ->enque_raw();
        err |= task_print ->enque_raw(); // DEBUG: just printing the forces before assempling
        err |= task_move  ->enque_raw();
        
    }
    
    ocl.download( ocl.ibuff_aforces, ff4.fapos, ff4.nvecs, ff4.nvecs*iSystemCur );
    ocl.download( ocl.ibuff_atoms,   ff4.apos , ff4.nvecs, ff4.nvecs*iSystemCur );
    err |=  ocl.finishRaw();                              //DEBUG
    OCL_checkError(err, "eval_MMFFf4_ocl_2");

    for(int i=0; i<ff4.nvecs; i++){  printf("OCL[%4i] f(%10.5f,%10.5f,%10.5f) p(%10.5f,%10.5f,%10.5f) pi %i \n", i, ff4.fapos[i].x,ff4.fapos[i].y,ff4.fapos[i].z,  ff4.apos[i].x,ff4.apos[i].y,ff4.apos[i].z,  i>=ff4.natoms ); }

    // ---- Compare to ffl
    bool ret=false;
    printf("### Compare ffl.fapos,  ff4.fapos   \n"); ret |= compareVecs( ff4.natoms, ffl.fapos,  ff4.fapos,  1e-4, true );
    if(ret){ printf("ERROR: GPU task_NBFF.eval() != ffl.nbmol.evalLJQs_ng4_PBC() => exit() \n"); exit(0); }else{ printf("CHECKED: GPU task_NBFF.eval() == ffl.nbmol.evalLJQs_ng4_PBC() \n"); }

    //printf("GPU AFTER assemble() \n"); ff4.printDEBUG( false,false );
    unpack( ff4.natoms, ffl.  apos, ff4.  apos );
    unpack( ff4.natoms, ffl. fapos, ff4. fapos );
    //opt.move_FIRE();
    
    
    //ff4.printDEBUG( false, false );
    // ============== CHECKS
    // ---- Check Invariatns
    fcog  = sum ( ffl.natoms, ffl.fapos   );
    tqcog = torq( ffl.natoms, ffl.apos, ffl.fapos );
    if(  fcog.norm2()>1e-8 ){ printf("WARRNING: eval_MMFFf4_ocl |fcog| =%g; fcog=(%g,%g,%g)\n", fcog.norm(),  fcog.x, fcog.y, fcog.z ); exit(0); }
    //if( tqcog.norm2()>1e-8 ){ printf("WARRNING: eval_MMFFf4_ocl |torq| =%g; torq=(%g,%g,%g)\n", tqcog.norm(),tqcog.x,tqcog.y,tqcog.z ); exit(0); }   // NOTE: torq is non-zero because pi-orbs have inertia
    
    exit(0);
    return 0;
}

// ======================== OTHER

void eval(){
    //printf("#======= MDloop[%i] \n", nloop );
    double E=0;
    setNonBond( bNonBonded );
    if(bGPU_MMFF){
        eval_MMFFf4_ocl( 1 );
        //eval_NBFF_ocl  ( 1 ); 
    }else{
        //printf( " ### CPU \n" );
        if(bMMFF){ E += ff.eval();  } 
        else     { VecN::set( nbmol.n*3, 0.0, (double*)nbmol.fs );  }
        if(bSurfAtoms){ 
            if  (bGridFF){ 
                //if(bOcl){ E+= eval_gridFF_ocl(nbmol.n, nbmol.ps,             nbmol.fs ); } 
                //else 
                { 
                    E+= gridFF.eval    (nbmol.n, nbmol.ps, nbmol.PLQs, nbmol.fs ); 
                    E+= nbmol.evalNeighs();
                }
            }else { 
                E+= nbmol.evalNeighs();   // Non-bonded interactions between atoms within molecule
              //E+= nbmol.evalMorse   (surf, false,                   gridFF.alpha, gridFF.Rdamp );
                E+= nbmol.evalMorsePBC( surf, gridFF.grid.cell, nPBC, gridFF.alpha, gridFF.Rdamp );
              //E+= nbmol.evalMorsePLQ( surf, gridFF.grid.cell, nPBC, gridFF.alpha, gridFF.Rdamp ); 
            }
        }
    }
    //for(int i=0; i<ff.natoms; i++){ printf("atom[%i] f(%g,%g,%g) \n", i, ff.fapos[i].x ,ff.fapos[i].y ,ff.fapos [i].z ); };
    //for(int i=0; i<ff.npi   ; i++){ printf("pvec[%i] f(%g,%g,%g) \n", i, ff.fpipos[i].x,ff.fpipos[i].y,ff.fpipos[i].z ); };
}

virtual void MDloop( int nIter, double Ftol = 1e-6 ) override {
    //printf( "MolWorld_sp3_ocl::MDloop(%i) bGridFF %i bOcl %i bMMFF %i \n", nIter, bGridFF, bOcl, bMMFF );
    //bMMFF=false;
    if(bMMFF)ff.cleanAll();
    for(int itr=0; itr<nIter; itr++){
        eval();
        //for(int i=0; i<nbmol.n; i++){ printf("atom[%i] f(%g,%g,%g)\n", i, nbmol.fs[i].x,nbmol.fs[i].y,nbmol.fs[i].z ); }
        ckeckNaN_d( nbmol.n, 3, (double*)nbmol.fs, "nbmol.fs" );
        if(ipicked>=0){
             float K = -2.0;
             Vec3d f = getForceSpringRay( ff.apos[ipicked], pick_hray, pick_ray0, K );
             ff.fapos[ipicked].add( f );
        };
        if( !bGPU_MMFF){ // update atomic positions
            //ff.fapos[  10 ].set(0.0); // This is Hack to stop molecule from moving
            //opt.move_GD(0.001);
            //opt.move_LeapFrog(0.01);
            //opt.move_MDquench();
            opt.move_FIRE();
        }
        double f2=1;
        if(f2<sq(Ftol)){
            bConverged=true;
        }
        nloop++;
    }
    bChargeUpdated=false;
}

virtual void swith_method()override{ 
    bGPU_MMFF=!bGPU_MMFF;    bOcl=bGPU_MMFF;
}

virtual char* info_str   ( char* str=0 ){ if(str==0)str=tmpstr; sprintf(str,"bGridFF %i bOcl %i \n", bGridFF, bOcl ); return str; }

// ==================================
//       Grid evaluation
// ==================================

void surf2ocl(Vec3i nPBC, bool bSaveDebug=false){
    int err=0;
    printf( "surf2ocl() na(%i) = ncell(%i) * natom(%i)\n", gridFF.natoms, nPBC );
    long T0=getCPUticks();
    Quat4f* atoms_surf = new Quat4f[gridFF.natoms];
    Quat4f* REQs_surf  = new Quat4f[gridFF.natoms];
    double R2damp = sq(gridFF.Rdamp);
    int ii=0;
    pack    ( gridFF.natoms, gridFF.apos,  atoms_surf );
    pack    ( gridFF.natoms, gridFF.aREQs, REQs_surf  );
    printf( ">>time(surf_to_GPU) %g \n", (getCPUticks()-T0)*tick2second );
    long T1=getCPUticks();
    ocl.makeGridFF( gridFF.grid, nPBC, gridFF.natoms, (float4*)atoms, (float4*)REQs, true );
    //ocl.addDipoleField( gridFF.grid, (float4*)dipole_ps, (float4*), true );
    printf( ">>time(ocl.makeGridFF() %g \n", (getCPUticks()-T1)*tick2second );
    delete [] atoms_surf;
    delete [] REQs_surf;
    if(bSaveDebug){
        gridFF.allocateFFs();
        ocl.download( ocl.itex_FE_Paul, gridFF.FFPauli  );  
        ocl.download( ocl.itex_FE_Lond, gridFF.FFLondon );
        ocl.download( ocl.itex_FE_Coul, gridFF.FFelec   );
        err =  ocl.finishRaw();    OCL_checkError(err, "surf2ocl.download.finish");
        printf( ">>time(ocl.surf2ocl.download() %g \n", (getCPUticks()-T1)*tick2second );
        saveGridXsfDebug();
    }
    ocl.buffers[ocl.ibuff_atoms_surf].release();
    ocl.buffers[ocl.ibuff_REQs_surf ].release();
}

virtual void initGridFF( const char * name, bool bGrid=true, bool bSaveDebugXSFs=false, double z0=NAN, Vec3d cel0={-0.5,-0.5,0.0}, bool bAutoNPBC=true )override{
    printf( "MolWorld_sp3_multi::initGridFF() \n");
    if(verbosity>0)printf("MolWorld_sp3_multi::initGridFF(%s,bGrid=%i,z0=%g,cel0={%g,%g,%g})\n",  name, bGrid, z0, cel0.x,cel0.y,cel0.z  );
    gridFF.grid.center_cell( cel0 );
    bGridFF=true;
    //gridFF.bindSystem      (surf.n, surf.atypes, surf.ps, surf.REQs );
    gridFF.setAtomsSymetrized(surf.n, surf.atypes, surf.ps, surf.REQs );
    if( isnan(z0) ){ z0=gridFF.findTop();   if(verbosity>0) printf("GridFF::findTop() %g \n", z0);  };
    gridFF.grid.pos0.z=z0;
    if(verbosity>1)gridFF.grid.printCell();
    if(bAutoNPBC){ autoNPBC( gridFF.grid.cell, nPBC, 30.0 ); }    
    long T0 = getCPUticks();
    surf2ocl( nPBC, bSaveDebugXSFs );
    printf( ">>time(init_ocl;GridFF_ocl): %g [s] \n", (getCPUticks()-T0)*tick2second  );
    bGridFF   =true; 
    //bSurfAtoms=false;
}

}; // class MolWorld_sp3_ocl

#endif
