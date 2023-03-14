
#ifndef MolWorld_sp3_ocl_h
#define MolWorld_sp3_ocl_h

#include "MolWorld_sp3.h"
#include "OCL_DFT.h"
//#include "OCL_PP.h"
#include "OCL_MM.h"
#include "datatypes_utils.h"

// ======================================
// class:        MolWorld_sp3_ocl
// ======================================

class MolWorld_sp3_multi : public MolWorld_sp3 { public:
    OCL_MM     ocl;

    int nSystems    = 1;
    int iSystemCur  = 0;    // currently selected system replica
    bool bGPU_MMFF = true;

    Quat4f* atoms      =0;
    Quat4f* aforces    =0;
    Quat4f* avel       =0;

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
    // --- params
    _realloc( neighs,    ocl.nAtoms*nSystems );
    _realloc( neighCell, ocl.nAtoms*nSystems );
    _realloc( bkNeighs,  ocl.nAtoms*nSystems );
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

void pack_system( int isys, MMFFsp3_loc& ff, bool bParams=0, bool bForces=0 ){
    printf("MolWorld_sp3_multi::pack_system(%i) \n", isys);
    //ocl.nvecs;
    int i0n = isys * ocl.nnode;
    int i0a = isys * ocl.nAtoms;
    int i0v = isys * ocl.nvecs;
    pack( ff.nvecs, ff.apos, atoms+i0v );
    if(bForces){
        pack( ff.nvecs, ff.fapos, aforces+i0v );
        //pack( f.nvecs, ff.avel, avel+i0v );
    }
    if(bParams){
        Mat3_to_cl( ff.   lvec,  lvecs[isys] );
        Mat3_to_cl( ff.invLvec, ilvecs[isys] );
        copy( ff.natoms, ff.aneighs,      neighs   +i0a );
        copy( ff.natoms, ff.aneighCell,   neighCell+i0a );
        copy( ff.natoms, ff.bkneighs,     bkNeighs +i0a );
        pack( nbmol.n,  nbmol.REQs,       REQs     +i0a );
        pack( ff.nnode, ff.apars, MMpars+i0n );
        pack( ff.nnode, ff.bLs,   BLs   +i0n );
        pack( ff.nnode, ff.bKs,   BKs   +i0n );
        pack( ff.nnode, ff.Ksp,   Ksp   +i0n );
        pack( ff.nnode, ff.Kpp,   Kpp   +i0n );
    }
}

void unpack_system(  int isys, MMFFsp3_loc& ff, bool bForces=0 ){
    printf("MolWorld_sp3_multi::unpack_system(%i) \n", isys);
    int i0n = isys * ocl.nnode;
    int i0a = isys * ocl.nAtoms;
    int i0v = isys * ocl.nvecs;
    unpack( ff.nvecs, ff.fapos, atoms+i0v );
    if(bForces){
        unpack( ff.nvecs, ff.fapos, aforces+i0v );
        //unpack( ff.nvecs, ff.fapos, avel+i0v );
    }
}

void upload(  bool bParams=0, bool bForces=0 ){
    printf("MolWorld_sp3_multi::upload() \n");
    ocl.upload( ocl.ibuff_atoms, atoms );
    if(bForces){
        ocl.upload( ocl.ibuff_aforces, aforces );
        //ocl.upload( ocl.ibuff_avel,    avel    );
    }
    if(bParams){
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
    ocl.init();
    ocl.makeKrenels_MM("common_resources/cl" );
    MolWorld_sp3::init(bGrid);

    // ----- init systems
    realloc( nSystems );
    for(int i=0; i<nSystems; i++){
        pack_system( i, ffl, true, false );
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
    task_move   = ocl.setup_updateAtomsMMFFf4( ff4.natoms, ff4.nnode       );   
    task_MMFF   = ocl.setup_getMMFFf4        ( ff4.natoms, ff4.nnode, bPBC );
    task_NBFF   = ocl.setup_getNonBond       ( ff4.natoms, ff4.nnode, nPBC, gridFF.Rdamp  );
    task_cleanF = ocl.setup_cleanForceMMFFf4 ( ff4.natoms, ff4.nnode       );
}

// ======================== EVAL OCL KERNEL functions

double eval_MMFFf4_ocl( int niter ){ 
    printf("MolWorld_sp3_multi::eval_MMFFf4_ocl() \n");
    if( task_MMFF==0 )setup_MMFFf4_ocl();
    for(int i=0; i<niter; i++){
        task_cleanF->enque_raw();  // DEBUG: this should be solved inside  task_move->enque_raw();
        task_MMFF  ->enque_raw();
        //task_NBFF  ->enque_raw();
        task_move->enque_raw(); //DEBUG
    }
    //printf( "ocl.download(n=%i) \n", n );
    ocl.download( ocl.ibuff_aforces, ff4.fapos, ff4.nvecs );
    ocl.download( ocl.ibuff_atoms,   ff4.apos , ff4.nvecs );
    //for(int i=0; i<ff4.natoms; i++){  printf("CPU[%i] p(%g,%g,%g) f(%g,%g,%g) \n", i, ff4.apos[i].x,ff4.apos[i].y,ff4.apos[i].z,  ff4.fapos[i].x,ff4.fapos[i].y,ff4.fapos[i].z ); }
    ocl.finishRaw();                              //DEBUG
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

    return 0;
}

// ======================== OTHER

void eval(){
    //printf("#======= MDloop[%i] \n", nloop );
    double E=0;
    //bSurfAtoms= false;
    if(bGPU_MMFF){
        eval_MMFFf4_ocl( 1 );
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

}; // class MolWorld_sp3_ocl

#endif
