﻿
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

// ==================================
//         Initialization
// ==================================

void realloc( int nSystems_ ){
    printf("MolWorld_sp3_multi::realloc() \n");
    nSystems=nSystems_;
    printf( "MolWorld_sp3_multi::realloc() Systems %i nAtoms %i nnode %i \n", nSystems, ffl.natoms,  ffl.nnode );
    ocl.initAtomsForces( nSystems, ffl.natoms,  ffl.nnode );
    //printf( "MolWorld_sp3_multi::realloc() Systems %i nAtoms %i nnode %i nvecs %i \n", nSystems, ocl.nAtoms, ocl.nnode, ocl.nvecs );
    // --- dynamical
    _realloc( atoms,     ocl.nvecs*nSystems  );
    _realloc( aforces,   ocl.nvecs*nSystems  );
    _realloc( avel,      ocl.nvecs*nSystems  );    for(int i=0; i<ocl.nvecs*nSystems; i++){ avel[i]=Quat4fZero; }  
    _realloc( constr,    ocl.nAtoms*nSystems );
    // --- params
    _realloc( neighs,    ocl.nAtoms*nSystems );
    _realloc( neighCell, ocl.nAtoms*nSystems );
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

virtual void init( bool bGrid ) override {
    int err = 0;
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
    //bGridFF=false;
    //bOcl   =false;
    bOcl   =true;
    setup_MMFFf4_ocl();
    int4 mask{1,1,0,0};
    //ocl.printOnGPU( 0,mask );
    //ocl.printOnGPU( 1,mask );
    //ocl.printOnGPU( 2,mask );
    printf("# ========== MolWorld_sp3_multi::init() DONE\n");
}

// ==================================
//         GPU <-> CPU I/O
// ==================================

void pack_system( int isys, MMFFsp3_loc& ff, bool bParams=0, bool bForces=0, bool bVel=false, float l_rnd=-1 ){
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
    if(bForces){ pack( ff.nvecs, ff.fapos, aforces+i0v ); }
    //if(bVel   ){ pack( ff.nvecs, opt.vel,  avel   +i0v ); }
    if(bParams){
        Mat3_to_cl( ff.   lvec,  lvecs[isys] );
        Mat3_to_cl( ff.invLvec, ilvecs[isys] );
        copy    ( ff.natoms, ff.neighCell, neighCell+i0a );
        copy    ( ff.natoms, ff.neighs,    neighs   +i0a );
        //copy_add( ff.natoms, ff.neighs,    neighs   +i0a,           0      );
        copy_add( ff.natoms, ff.bkneighs,   bkNeighs +i0v,           i0n*8  );
        copy_add( ff.nnode , ff.bkneighs,   bkNeighs +i0v+ff.natoms, i0n*8 + 4*ff.nnode );
        pack    ( ff.nnode , ff.apars,      MMpars   +i0n );
        pack    ( ff.nnode , ff.bLs,        BLs      +i0n );
        pack    ( ff.nnode , ff.bKs,        BKs      +i0n );
        pack    ( ff.nnode , ff.Ksp,        Ksp      +i0n );
        pack    ( ff.nnode , ff.Kpp,        Kpp      +i0n );
        pack    ( nbmol.natoms, nbmol.REQs, REQs        +i0a, fabs(gridFF.alpha) );
    }

    //if(isys==5)
    {
        //for(int i=0; i<ocl.nAtoms; i++){ Quat4i ng = neighs[i0a+i]; printf( "ngs[%4i] (%4i,%4i,%4i,%4i) \n", i, ng.x,ng.y,ng.z,ng.w ); }
        //for(int i=0; i<ocl.nvecs; i++){ Quat4f p = atoms[i0v+i]; printf( "apos[%4i] (%5.3f,%5.3f,%5.3f) pi %i \n", i, p.x,p.y,p.z, i>=ocl.nAtoms ); }
        //if(isys==0)for(int i=0; i<ocl.nvecs; i++){Quat4i bkng = bkNeighs[i0v+i]; printf( "bkng[%4i] (%4i,%4i,%4i,%4i) pi %i\n", i, bkng.x,bkng.y,bkng.z,bkng.w, i>=ocl.nAtoms ); }
    }
}

void unpack_system(  int isys, MMFFsp3_loc& ff, bool bForces=0, bool bVel=false ){
    printf("MolWorld_sp3_multi::unpack_system(%i) \n", isys);
    int i0n = isys * ocl.nnode;
    int i0a = isys * ocl.nAtoms;
    int i0v = isys * ocl.nvecs;
    unpack( ff.nvecs, ff.apos, atoms+i0v );
    if(bForces){ unpack( ff.nvecs, ff.fapos, aforces+i0v ); }
    if(bVel   ){ unpack( ff.nvecs, ff.fapos, avel   +i0v ); }
}

void upload(  bool bParams=0, bool bForces=0, bool bVel=true ){
    printf("MolWorld_sp3_multi::upload() \n");
    ocl.upload( ocl.ibuff_atoms,  atoms  );
    ocl.upload( ocl.ibuff_constr, constr );
    if(bForces){ ocl.upload( ocl.ibuff_aforces, aforces ); }
    if(bVel   ){ ocl.upload( ocl.ibuff_avel,    avel    ); }
    if(bParams){
        ocl.upload( ocl.ibuff_lvecs,   lvecs );
        ocl.upload( ocl.ibuff_ilvecs, ilvecs );
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

void download( bool bForces=0, bool bVel=false ){
    printf("MolWorld_sp3_multi::download() \n");
    ocl.download( ocl.ibuff_atoms, atoms );
    if(bForces){ ocl.download( ocl.ibuff_aforces, aforces ); }
    if(bVel   ){ ocl.download( ocl.ibuff_avel,    avel    ); }
}

// ==================================
//          SETUP OCL KERNELS 
// ==================================

void setup_MMFFf4_ocl(){
    printf("MolWorld_sp3_multi::setup_MMFFf4_ocl() \n");
    ocl.nDOFs.x=ff.natoms;
    ocl.nDOFs.y=ff.nnode;
    if(!task_move  )   task_move  = ocl.setup_updateAtomsMMFFf4( ff4.natoms, ff4.nnode       );
    if(!task_print )   task_print = ocl.setup_printOnGPU       ( ff4.natoms, ff4.nnode       );
    if(!task_MMFF  )   task_MMFF  = ocl.setup_getMMFFf4        ( ff4.natoms, ff4.nnode, bPBC );
    if(!task_NBFF  ) { 
        if( bGridFF ){ task_NBFF  = ocl.setup_getNonBond_GridFF( ff4.natoms, ff4.nnode, nPBC, gridFF.Rdamp ); } 
        else         { task_NBFF  = ocl.setup_getNonBond       ( ff4.natoms, ff4.nnode, nPBC, gridFF.Rdamp ); }
    }
    if(!task_cleanF)task_cleanF = ocl.setup_cleanForceMMFFf4 ( ff4.natoms, ff4.nnode       );
}

void setup_NBFF_ocl(){
    printf("MolWorld_sp3_multi::setup_NBFF_ocl() \n");
    ocl.nDOFs.x=ff.natoms;
    ocl.nDOFs.y=ff.nnode;
    if(!task_cleanF )task_cleanF = ocl.setup_cleanForceMMFFf4 ( ff4.natoms, ff4.nnode        );
    if(!task_move   )task_move   = ocl.setup_updateAtomsMMFFf4( ff4.natoms, ff4.nnode        ); 
    if(!task_print  )task_print  = ocl.setup_printOnGPU       ( ff4.natoms, ff4.nnode        );
    if(!task_NBFF   )task_NBFF   = ocl.setup_getNonBond       ( ff4.natoms, ff4.nnode, nPBC, gridFF.Rdamp  );
}

void picked2GPU( int ipick,  double K ){
    //printf( "picked2GPU() ipick %i iSystemCur %i \n", ipick, iSystemCur );
    int i0a = ocl.nAtoms*iSystemCur;
    int i0v = ocl.nvecs *iSystemCur;
    if(ipick>=0){
        Quat4f& acon = constr[i0a + ipick];
        Vec3f hray   = (Vec3f)pick_hray;
        Vec3f ray0   = (Vec3f)pick_ray0;
        const Quat4f& atom = atoms [i0v + ipick];
        float c = hray.dot( atom.f ) - hray.dot( ray0 );
        acon.f = ray0 + hray*c;
        acon.w = K;
    }else{
        for(int i=0; i<ocl.nAtoms; i++){   constr[i0a + i].w=-1.0;  };
    }
    //for(int i=0; i<ocl.nAtoms; i++){ printf( "CPU:constr[%i](%7.3f,%7.3f,%7.3f |K= %7.3f) \n", i, constr[i0a+i].x,constr[i0a+i].y,constr[i0a+i].z,  constr[i0a+i].w   ); }
    ocl.upload( ocl.ibuff_constr, constr );   // ToDo: instead of updating the whole buffer we may update just relevant part?
}

// ==================================
//           eval @GPU 
// ==================================

double eval_MMFFf4_ocl( int niter, bool bForce=false ){ 
    //printf("MolWorld_sp3_multi::eval_MMFFf4_ocl() niter=%i \n", niter );
    picked2GPU( ipicked,  1.0 );
    int err=0;
    if( task_MMFF==0 )setup_MMFFf4_ocl();
    // evaluate on GPU
    for(int i=0; i<niter; i++){
        err |= task_cleanF->enque_raw();  // DEBUG: this should be solved inside  task_move->enque_raw();
        err |= task_MMFF  ->enque_raw();
        err |= task_NBFF  ->enque_raw();
        //err |= task_print ->enque_raw(); // DEBUG: just printing the forces before assempling
        err |= task_move  ->enque_raw(); 
        //OCL_checkError(err, "eval_MMFFf4_ocl_1");
    }
    if(bForce)ocl.download( ocl.ibuff_aforces, ff4.fapos, ff4.nvecs, ff4.nvecs*iSystemCur );
    ocl          .download( ocl.ibuff_atoms,   ff4.apos , ff4.nvecs, ff4.nvecs*iSystemCur );
    //for(int i=0; i<ff4.nvecs; i++){  printf("OCL[%4i] f(%10.5f,%10.5f,%10.5f) p(%10.5f,%10.5f,%10.5f) pi %i \n", i, ff4.fapos[i].x,ff4.fapos[i].y,ff4.fapos[i].z,  ff4.apos[i].x,ff4.apos[i].y,ff4.apos[i].z,  i>=ff4.natoms ); }

    err |= ocl.finishRaw();
    OCL_checkError(err, "eval_MMFFf4_ocl");
    unpack( ff4.natoms, ffl.  apos, ff4.  apos );
    unpack( ff4.nnode,  ffl. pipos, ff4. pipos );
    if(bForce){
        unpack( ff4.natoms, ffl. fapos, ff4. fapos );
        unpack( ff4.nnode,  ffl.fpipos, ff4.fpipos );
        // ---- Check Invariatns   - this works only if we unpack forces
        fcog  = sum ( ffl.natoms, ffl.fapos   );
        tqcog = torq( ffl.natoms, ffl.apos, ffl.fapos );
        if(  fcog.norm2()>1e-8 ){ printf("WARRNING: eval_MMFFf4_ocl |fcog| =%g; fcog=(%g,%g,%g)\n", fcog.norm(),  fcog.x, fcog.y, fcog.z ); exit(0); }
        //if( tqcog.norm2()>1e-8 ){ printf("WARRNING: eval_MMFFf4_ocl |torq| =%g; torq=(%g,%g,%g)\n", tqcog.norm(),tqcog.x,tqcog.y,tqcog.z ); exit(0); }   // NOTE: torq is non-zero because pi-orbs have inertia
    }
    return 0;
}

double eval_NBFF_ocl( int niter, bool bForce=false ){ 
    //printf("MolWorld_sp3_multi::eval_NBFF_ocl() \n");
    int err=0;
    if( task_NBFF==0 ){ setup_NBFF_ocl(); }
    // evaluate on GPU
    for(int i=0; i<niter; i++){
        //err |= task_cleanF->enque_raw(); // DEBUG: this should be solved inside  task_move->enque_raw();
        err |= task_NBFF  ->enque_raw();
        //err |= task_print ->enque_raw(); // DEBUG: just printing the forces before assempling
        err |= task_move  ->enque_raw();
        //OCL_checkError(err, "eval_NBFF_ocl.1");
    }
    if(bForce)ocl.download( ocl.ibuff_aforces, ff4.fapos, ff4.nvecs, ff4.nvecs*iSystemCur );
    ocl.download( ocl.ibuff_atoms,   ff4.apos , ff4.nvecs, ff4.nvecs*iSystemCur );
    err |=  ocl.finishRaw();                              //DEBUG
    OCL_checkError(err, "eval_NBFF_ocl.2");
    if(bForce){
        unpack( ff4.natoms, ffl.  apos, ff4.  apos );
        unpack( ff4.natoms, ffl. fapos, ff4. fapos );
        fcog  = sum ( ffl.natoms, ffl.fapos   );
        tqcog = torq( ffl.natoms, ffl.apos, ffl.fapos );
        if(  fcog.norm2()>1e-8 ){ printf("WARRNING: eval_NBFF_ocl |fcog| =%g; fcog=(%g,%g,%g)\n", fcog.norm(),  fcog.x, fcog.y, fcog.z ); exit(0); }
        //if( tqcog.norm2()>1e-8 ){ printf("WARRNING: eval_MMFFf4_ocl |torq| =%g; torq=(%g,%g,%g)\n", tqcog.norm(),tqcog.x,tqcog.y,tqcog.z ); exit(0); }   // NOTE: torq is non-zero because pi-orbs have inertia
    }
    return 0;
}

// ==================================
//                 eval
// ==================================

double eval( ){
    double E=0;
    setNonBond( bNonBonded );  // Make sure ffl subtracts non-covalent interction for angles
    if(bMMFF){ 
        //E += ff .eval();
        E += ffl.eval();
        //E += eval_f4();
        //printf( "atom[0] nbmol(%g,%g,%g) ff(%g,%g,%g) ffl(%g,%g,%g) \n", nbmol.apos[0].x,nbmol.apos[0].y,nbmol.apos[0].z,  ff.apos[0].x,ff.apos[0].y,ff.apos[0].z,  ffl.apos[0].x,ffl.apos[0].y,ffl.apos[0].z );
        
    }else{ VecN::set( nbmol.natoms*3, 0.0, (double*)nbmol.fapos );  }
    if(bNonBonded){
        if(bMMFF){    
            //if  (bPBC){ E += nbmol.evalLJQs_ng4_PBC( ffl.neighs, ffl.neighCell, ff.lvec, {1,1,0} );             }   // atoms outside cell
            if  (bPBC){ E += nbmol.evalLJQs_ng4_PBC( ffl.neighs, ffl.neighCell, npbc, pbc_shifts, gridFF.Rdamp ); }  
            else      { E += nbmol.evalLJQs_ng4    ( ffl.neighs );                                                }   // atoms in cell ignoring bondede neighbors       
        }else{
            if  (bPBC){ E += nbmol.evalLJQs_PBC    ( ff.lvec, {1,1,0} ); }   // atoms outside cell
            else      { E += nbmol.evalLJQs        ( );                  }   // atoms in cell ignoring bondede neighbors    
        }
    }
    if(bConstrains)constrs.apply( nbmol.apos, nbmol.fapos );
    /*
    if(bSurfAtoms){ 
        if   (bGridFF){ E+= gridFF.eval(nbmol.natoms, nbmol.apos, nbmol.PLQs, nbmol.fapos ); }
        //else        { E+= nbmol .evalMorse   ( surf, false,                  gridFF.alpha, gridFF.Rdamp );  }
        else          { E+= nbmol .evalMorsePBC( surf, gridFF.grid.cell, nPBC, gridFF.alpha, gridFF.Rdamp );  }
    }
    */
    //printf( "eval() bSurfAtoms %i bGridFF %i \n", bSurfAtoms, bGridFF );
    //for(int i=0; i<nbmol.natoms; i++){ printf("atom[%i] f(%g,%g,%g)\n", i, nbmol.fapos[i].x,nbmol.fapos[i].y,nbmol.fapos[i].z ); }
    return E;
}

// ==================================
//                 MDloop
// ==================================

virtual void MDloop( int nIter, double Ftol = 1e-6 ) override {
    //printf( "MolWorld_sp3_ocl::MDloop(%i) bGridFF %i bOcl %i bMMFF %i \n", nIter, bGridFF, bOcl, bMMFF );
    //bMMFF=false;
    if( bOcl ){
        //printf( "GPU frame[%i] -- \n", nIter );
        if( (iSystemCur<0) || (iSystemCur>=nSystems) ){  printf("ERROR: iSystemCur(%i) not in range [ 0 .. nSystems(%i) ] => exit() \n", iSystemCur, nSystems ); exit(0); }
        eval_MMFFf4_ocl( nIter );
        //eval_NBFF_ocl  ( 1 ); 
    }else{
        printf( "CPU frame[%i] \n", nIter );
        for(int itr=0; itr<nIter; itr++){
            double E = eval();
            //double E = MolWorld_sp3::eval();
            ckeckNaN_d( nbmol.natoms, 3, (double*)nbmol.fapos, "nbmol.fapos" );
            //if( bPlaneSurfForce )for(int i=0; i<ff.natoms; i++){ ff.fapos[i].add( getForceMorsePlane( ff.apos[i], {0.0,0.0,1.0}, -5.0, 0.0, 0.01 ) ); }
            //printf( "apos(%g,%g,%g) f(%g,%g,%g)\n", ff.apos[0].x,ff.apos[0].y,ff.apos[0].z,   ff.fapos[0].x,ff.fapos[0].y,ff.fapos[0].z );
            //if(bCheckInvariants){ checkInvariants(maxVcog,maxFcog,maxTg); }
            if(ipicked>=0){ pullAtom( ipicked );  }; // printf( "pullAtom(%i) E=%g\n", ipicked, E ); };
            //ff.fapos[  10 ].set(0.0); // This is Hack to stop molecule from moving
            //opt.move_GD(0.001);
            //opt.move_LeapFrog(0.01);
            //opt.move_MDquench();
            double f2=opt.move_FIRE();   
            //printf( "[%i] E= %g [eV] |F|= %g [eV/A]\n", nloop, E, sqrt(f2) );
            //double f2=1;
            if(f2<sq(Ftol)){
                bConverged=true;
            }
            nloop++;
        }
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
    printf( "surf2ocl() gridFF.natoms=%i nPBC(%i,%i,%i)\n", gridFF.natoms, nPBC.x,nPBC.y,nPBC.z );
    long T0=getCPUticks();
    Quat4f* atoms_surf = new Quat4f[gridFF.natoms];
    Quat4f* REQs_surf  = new Quat4f[gridFF.natoms];
    pack( gridFF.natoms, gridFF.apos,  atoms_surf,   sq(gridFF.Rdamp) );
    pack( gridFF.natoms, gridFF.REQs, REQs_surf , fabs(gridFF.alpha) );
    long T1=getCPUticks();
    ocl.makeGridFF( gridFF.grid, nPBC, gridFF.natoms, (float4*)atoms_surf, (float4*)REQs_surf, true );
    //ocl.addDipoleField( gridFF.grid, (float4*)dipole_ps, (float4*), true );
    printf( ">>time(ocl.makeGridFF() %g \n", (getCPUticks()-T1)*tick2second );
    delete [] atoms_surf;
    delete [] REQs_surf;
    bSaveDebug=true;
    if(bSaveDebug){
        gridFF.allocateFFs();
        ocl.download( ocl.itex_FE_Paul, gridFF.FFPauli  );  
        ocl.download( ocl.itex_FE_Lond, gridFF.FFLondon );
        ocl.download( ocl.itex_FE_Coul, gridFF.FFelec   );
        err =  ocl.finishRaw();    OCL_checkError(err, "surf2ocl.download.finish");
        printf( ">>time(surf2ocl.download() %g \n", (getCPUticks()-T1)*tick2second );
        saveGridXsfDebug();
    }
    ocl.buffers[ocl.ibuff_atoms_surf].release();
    ocl.buffers[ocl.ibuff_REQs_surf ].release();
    exit(0);
}

virtual void initGridFF( const char * name, bool bGrid=true, bool bSaveDebugXSFs=false, double z0=NAN, Vec3d cel0={-0.5,-0.5,0.0}, bool bAutoNPBC=true )override{
    printf( "MolWorld_sp3_multi::initGridFF() \n");
    if(verbosity>0)printf("MolWorld_sp3_multi::initGridFF(%s,bGrid=%i,z0=%g,cel0={%g,%g,%g})\n",  name, bGrid, z0, cel0.x,cel0.y,cel0.z  );
    gridFF.grid.center_cell( cel0 );
    bGridFF=true;
    gridFF.bindSystem      (surf.natoms, surf.atypes, surf.apos, surf.REQs );
    //gridFF.setAtomsSymetrized(surf.natoms, surf.atypes, surf.apos, surf.REQs );
    gridFF.evalCellDipole();
    if( ( fabs(gridFF.Q)>1e-6 ) || (gridFF.dip.norm2()>1e-8) ){ printf("ERROR: GridFF has dipole and dipole correction not yet implemented => exit() \n"); exit(0); }
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

// ##############################################################
// ##############################################################
//                 Debugging versions 
// ##############################################################
// ##############################################################


// ==============================================================
//                   eval_MMFFf4_ocl_debug    
// ==============================================================

double eval_MMFFf4_ocl_debug( int niter ){ 
    printf("MolWorld_sp3_multi::eval_MMFFf4_ocl_debug() \n");
    int err=0;
    if( task_MMFF==0 )setup_MMFFf4_ocl();
    
    {  // DEBUG --- evaluate on CPU
        //bool bEval_ffl = false;
        bool bEval_ffl = true;
        unpack_system( iSystemCur, ffl );
        ffl.Rdamp = gridFF.Rdamp;
        ffl.cleanForce();
        ffl.eval();    //for(int i=0; i<ffl.nvecs; i++){  printf("ffl[%4i] f(%10.5f,%10.5f,%10.5f) p(%10.5f,%10.5f,%10.5f) pi %i \n", i,  ffl.fapos[i].x,ffl.fapos[i].y,ffl.fapos[i].z,   ffl.apos[i].x,ffl.apos[i].y,ffl.apos[i].z,  i>=ffl.natoms ); }
        //nbmol.evalLJQs_ng4_PBC( ffl.neighs, ffl.neighCell, ffl.lvec, ffl.nPBC, gridFF.Rdamp );
        nbmol.evalLJQs_ng4_PBC( ffl.neighs, ffl.neighCell, npbc, pbc_shifts, gridFF.Rdamp ); 
        fcog  = sum ( ffl.natoms, ffl.fapos   );
        tqcog = torq( ffl.natoms, ffl.apos, ffl.fapos );
        if(  fcog.norm2()>1e-8 ){ printf("WARRNING: ffl.eval_MMFFf4_ocl_debug() CPU |fcog| =%g; fcog=(%g,%g,%g) bEval_ffl %i \n", bEval_ffl, fcog.norm(),  fcog.x, fcog.y, fcog.z, bEval_ffl ); exit(0); }
    }
    
    /*
    { // DEBUG --- evaluate on CPU
        unpack_system( iSystemCur, ffl );
        ffl  .cleanForce();
        //nbmol.evalLJQs_ng4_PBC( ffl.neighs, ffl.neighCell, ffl.lvec, ffl.nPBC, gridFF.Rdamp );
        nbmol.evalLJQs_ng4_PBC( ffl.neighs, ffl.neighCell, npbc, pbc_shifts, gridFF.Rdamp ); } 
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
        OCL_checkError(err, "eval_MMFFf4_ocl_debug.1");
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
    OCL_checkError(err, "eval_MMFFf4_ocl_debug.2");

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

// ==============================================================
//                   eval_NBFF_ocl_debug    
// ==============================================================

double eval_NBFF_ocl_debug( int niter ){ 
    printf("MolWorld_sp3_multi::eval_NBFF_ocl_debug() \n");
    int err=0;
    if( task_NBFF==0 ){ setup_NBFF_ocl(); }

    { // DEBUG --- evaluate on CPU
        unpack_system( iSystemCur, ffl );
        ffl  .cleanForce();
        //nbmol.evalLJQs_ng4_PBC( ffl.neighs, ffl.neighCell, ffl.lvec, ffl.nPBC, gridFF.Rdamp );
        nbmol.evalLJQs_ng4_PBC( ffl.neighs, ffl.neighCell, npbc, pbc_shifts, gridFF.Rdamp ); 
        fcog  = sum ( ffl.natoms, ffl.fapos   );
        tqcog = torq( ffl.natoms, ffl.apos, ffl.fapos );
        if(  fcog.norm2()>1e-8 ){ printf("WARRNING: eval_NBFF_ocl_debug() CPU |fcog| =%g; fcog=(%g,%g,%g)\n", fcog.norm(),  fcog.x, fcog.y, fcog.z ); exit(0); }
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
    OCL_checkError(err, "eval_NBFF_ocl_debug");

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
    if(  fcog.norm2()>1e-8 ){ printf("WARRNING: eval_NBFF_ocl_debug |fcog| =%g; fcog=(%g,%g,%g)\n", fcog.norm(),  fcog.x, fcog.y, fcog.z ); exit(0); }
    //if( tqcog.norm2()>1e-8 ){ printf("WARRNING: eval_MMFFf4_ocl |torq| =%g; torq=(%g,%g,%g)\n", tqcog.norm(),tqcog.x,tqcog.y,tqcog.z ); exit(0); }   // NOTE: torq is non-zero because pi-orbs have inertia
    
    exit(0);
    return 0;
}


}; // class MolWorld_sp3_ocl

#endif
