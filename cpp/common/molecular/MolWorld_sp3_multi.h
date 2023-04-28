﻿
#ifndef MolWorld_sp3_ocl_h
#define MolWorld_sp3_ocl_h

#include "MolWorld_sp3.h"
//#include "OCL_DFT.h"
//#include "OCL_PP.h"
#include "OCL_MM.h"
#include "datatypes_utils.h"

#include "MultiSolverInterface.h"
#include "Confs.h"

struct FIRE_setup{
    float cvf_min   = -0.1f;  // minimum cosine for velocity damping in  move_FIRE_smooth()
    float cvf_max   = +0.1f;  // maximum cosine for velocity damping in  move_FIRE_smooth()
    float cv_kill   =  0.0f;
    float ff_safety = 1e-32;

    int   minLastNeg  = 3;
    float finc         = 1.1;
    float fdec         = 0.5;
    float falpha       = 0.8;

    float dt_max       = 0;
    float dt_min       = 0;
    float damp_max     = 0;

    inline void setup( float dt_, float damp_=0.1 ){
        dt_max=dt_; dt_min=dt_*0.1; damp_max=damp_;
    }

    FIRE_setup()=default;
    FIRE_setup(float dt_, float damp_=0.1){ setup(dt_,damp_); };
};

struct FIRE{
    constexpr static float ff_safety = 1.e-32f;
    FIRE_setup* par=0;

    int   lastNeg;
    float ff,vv,vf;
    float cos_vf;
    float renorm_vf;
    float cf,cv;
    float dt;
    float damping;

    float damp_func_FIRE( float c, float& cv ){ // Original FIRE damping
        double cf;
        if( c < 0.0 ){
            cv = par->cv_kill;
            cf = 0.;
        }else{  
            cv = 1-damping;
            cf = damping;
            //cv = 1.;
            //cf = 0.;
        }
        return cf;
    }

    void update_params(){
        float f_len      = sqrt(ff);
        float v_len      = sqrt(vv);
        cos_vf     = vf    / ( f_len*v_len + ff_safety );
        renorm_vf  = v_len / ( f_len       + ff_safety );
        cf  = renorm_vf * damp_func_FIRE( cos_vf, cv );
        if( cos_vf < 0.0 ){
            //dt       = fmax( dt * par->fdec, par->dt_min );
            lastNeg  = 0;
            damping  = par->damp_max;
        }else{
            if( lastNeg > par->minLastNeg ){
                //dt        = fmin( dt * par->finc, par->dt_max );
                damping   = damping  * par->falpha;
            }
            lastNeg++;
        }
    }

    void print(){ printf( "dt %g damping %g cos_vf %g cv %g cf %g \n", dt, damping, cos_vf, cv, cf ); }
    void bind_params( FIRE_setup* par_ ){ par=par_; dt=par->dt_max; damping=par->damp_max; };

    FIRE()=default;
    FIRE( FIRE_setup* par_ ){ bind_params(par_); }

};



// ======================================
// class:        MolWorld_sp3_ocl
// ======================================

class MolWorld_sp3_multi : public MolWorld_sp3, public MultiSolverInterface { public:
    OCL_MM     ocl;
    FIRE_setup fire_setup{0.1,0.1};

    int nSystems    = 1;
    //int iSystemCur  = 5;    // currently selected system replica
    //int iSystemCur  = 8;    // currently selected system replica
    bool bGPU_MMFF = true;

    Quat4f* atoms      =0;
    Quat4f* aforces    =0;
    Quat4f* avel       =0;

    FIRE*   fire         =0;
    Quat4f* MDpars       =0;

    Quat4f* constr     =0;

    Quat4i* neighs     =0;
    Quat4i* neighCell  =0;
    Quat4i* bkNeighs   =0;
    Quat4i* bkNeighs_new=0;
    //Quat4f* neighForce =0;

    Quat4f* REQs       =0;
    Quat4f* MMpars     =0;
    Quat4f* BLs        =0;
    Quat4f* BKs        =0;
    Quat4f* Ksp        =0;
    Quat4f* Kpp        =0;

    cl_Mat3*  lvecs    =0;
    cl_Mat3* ilvecs    =0;
    Quat4f*   pbcshifts =0;

    OCLtask* task_cleanF=0;
    OCLtask* task_NBFF=0;
    OCLtask* task_NBFF_Grid=0;
    OCLtask* task_MMFF=0;
    OCLtask* task_move=0;
    OCLtask* task_print=0;
    OCLtask* task_MMFFloc1=0;
    OCLtask* task_MMFFloc2=0;
    OCLtask* task_MMFFloc_test=0;

// ==================================
//         Initialization
// ==================================

void realloc( int nSystems_ ){
    printf("MolWorld_sp3_multi::realloc() \n");
    nSystems=nSystems_;
    printf( "MolWorld_sp3_multi::realloc() Systems %i nAtoms %i nnode %i \n", nSystems, ffl.natoms,  ffl.nnode );
    ocl.initAtomsForces( nSystems, ffl.natoms,  ffl.nnode, npbc );
    //printf( "MolWorld_sp3_multi::realloc() Systems %i nAtoms %i nnode %i nvecs %i \n", nSystems, ocl.nAtoms, ocl.nnode, ocl.nvecs );
    // --- dynamical
    _realloc ( atoms,     ocl.nvecs*nSystems  );
    _realloc0( aforces,   ocl.nvecs*nSystems  ,Quat4fZero);
    _realloc0( avel,      ocl.nvecs*nSystems  ,Quat4fZero);
    _realloc ( constr,    ocl.nAtoms*nSystems );
    // --- params
    _realloc( neighs,    ocl.nAtoms*nSystems );
    _realloc( neighCell, ocl.nAtoms*nSystems );
    _realloc( bkNeighs,    ocl.nvecs*nSystems );
    _realloc( bkNeighs_new,ocl.nvecs*nSystems );
    _realloc( REQs,      ocl.nAtoms*nSystems );
    _realloc( MMpars,    ocl.nnode*nSystems  );
    _realloc( BLs,       ocl.nnode*nSystems  );
    _realloc( BKs,       ocl.nnode*nSystems  );
    _realloc( Ksp,       ocl.nnode*nSystems  );
    _realloc( Kpp,       ocl.nnode*nSystems  );

    _realloc( lvecs,     nSystems  );
    _realloc( ilvecs,    nSystems  );
    _realloc( MDpars,    nSystems  );

    _realloc( pbcshifts, ocl.npbc*nSystems );

    _realloc( fire,      nSystems  );

    // ToDo : it may be good to bind buffer directly in p_cpu buffer inside   OCLsystem::newBuffer()
}

virtual void init( bool bGrid ) override {
    int err = 0;
    printf("# ========== MolWorld_sp3_multi::init() START\n");
    gopt.msolver = this;
    int i_nvidia = ocl.print_devices(true);
    ocl.init(i_nvidia);
    ocl.makeKrenels_MM("common_resources/cl" );
    MolWorld_sp3::init(bGrid);
    // ----- init systems
    realloc( nSystems );
    //if(bGridFF) evalCheckGridFF_ocl();  // this must be after we make buffers but before we fill them
    float random_init = 0.5;
    for(int i=0; i<nSystems; i++){
        pack_system( i, ffl, true, false, random_init );
        fire[i].bind_params( &fire_setup );
    }
    for(int i=0;i<nSystems; i++){ printMSystem(i); }
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
    //printf("MolWorld_sp3_multi::pack_system(%i) \n", isys);
    //ocl.nvecs;
    int i0n = isys * ocl.nnode;
    int i0a = isys * ocl.nAtoms;
    int i0v = isys * ocl.nvecs;
    int i0pbc = isys*ocl.npbc;
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

    Mat3d lvec0=ff.lvec;  // DEBUG

    if(bForces){ pack( ff.nvecs, ff.fapos, aforces+i0v ); }
    //if(bVel   ){ pack( ff.nvecs, opt.vel,  avel   +i0v ); }
    if(bParams){
        
        //Mat3d lvec=lvec0; lvec.a.addRandomCube(0.5); ffl.setLvec(lvec);  // DEBUG
        //Mat3d lvec=lvec0; lvec.b.x += randf(-0.5,0.5); ffl.setLvec(lvec);  // DEBUG
        //Mat3d lvec=lvec0; lvec.b.x += isys*0.00163; ffl.setLvec(lvec);  // DEBUG
        //Mat3d lvec=lvec0; lvec.b.x += isys*0.15; ffl.setLvec(lvec);  // DEBUG

        Mat3d lvec=lvec0; lvec.a.y += isys*2.25; ffl.setLvec(lvec);  // DEBUG

        Mat3_to_cl( ff.   lvec,  lvecs[isys] );
        Mat3_to_cl( ff.invLvec, ilvecs[isys] );
        evalPBCshifts( nPBC, ff.lvec, pbcshifts + i0pbc );

        copy    ( ff.natoms, ff.neighCell, neighCell+i0a );
        copy    ( ff.natoms, ff.neighs,    neighs   +i0a );
        //copy_add( ff.natoms, ff.neighs,    neighs   +i0a,           0      );
        copy    ( ff.natoms, ff.bkneighs,  bkNeighs_new +i0v           );
        copy    ( ff.nnode,  ff.bkneighs,  bkNeighs_new +i0v+ff.natoms );
        copy_add( ff.natoms, ff.bkneighs,   bkNeighs +i0v,           i0n*8  );
        copy_add( ff.nnode , ff.bkneighs,   bkNeighs +i0v+ff.natoms, i0n*8 + 4*ff.nnode );
        pack    ( ff.nnode , ff.apars,      MMpars   +i0n );
        pack    ( ff.nnode , ff.bLs,        BLs      +i0n );
        pack    ( ff.nnode , ff.bKs,        BKs      +i0n );
        pack    ( ff.nnode , ff.Ksp,        Ksp      +i0n );
        pack    ( ff.nnode , ff.Kpp,        Kpp      +i0n );
        pack    ( nbmol.natoms, nbmol.REQs, REQs     +i0a );
    }

    //if(isys==5)
    {
        //for(int i=0; i<ocl.nAtoms; i++){ Quat4i ng = neighs[i0a+i]; printf( "ngs[%4i] (%4i,%4i,%4i,%4i) \n", i, ng.x,ng.y,ng.z,ng.w ); }
        //for(int i=0; i<ocl.nvecs; i++){ Quat4f p = atoms[i0v+i]; printf( "apos[%4i] (%5.3f,%5.3f,%5.3f) pi %i \n", i, p.x,p.y,p.z, i>=ocl.nAtoms ); }
        //if(isys==0)for(int i=0; i<ocl.nvecs; i++){Quat4i bkng = bkNeighs[i0v+i]; printf( "bkng[%4i] (%4i,%4i,%4i,%4i) pi %i\n", i, bkng.x,bkng.y,bkng.z,bkng.w, i>=ocl.nAtoms ); }
    }
}

void evalVF( int n, Quat4f* fs, Quat4f* vs, FIRE& fire, Quat4f& MDpar ){
    double vv=0,ff=0,vf=0;
    for(int i=0; i<n; i++){
        Vec3f f = fs[i].f;
        Vec3f v = vs[i].f;
        ff += f.dot(f);
        vv += v.dot(v);
        vf += f.dot(v);
    }
    fire.vv=vv;
    fire.ff=ff;
    fire.vf=vf;
    fire.update_params();
    MDpar.x = fire.dt;
    //MDpar.x = fire.dt * 0.5; // DEBUG
    //MDpar.x = 0.01; // DEBUG

    MDpar.y = 1 - fire.damping;
    MDpar.z = fire.cv;
    MDpar.w = fire.cf;
}

double evalVFs(){
    int err=0;
    ocl.download( ocl.ibuff_aforces , aforces );
    ocl.download( ocl.ibuff_avel    , avel    );
    err |= ocl.finishRaw();  OCL_checkError(err, "evalVFs()");
    //printf("MolWorld_sp3_multi::evalVFs(%i) \n", isys);
    double F2max = 0;
    for(int isys=0; isys<nSystems; isys++){
        int i0v = isys * ocl.nvecs;
        evalVF( ocl.nvecs, aforces+i0v, avel   +i0v, fire[isys], MDpars[isys] );
        F2max = fmax( F2max, fire[isys].ff );
    }
    //printf( "MDpars{%g,%g,%g,%g}\n", MDpars[0].x,MDpars[0].y,MDpars[0].z,MDpars[0].w );
    err |= ocl.upload( ocl.ibuff_MDpars, MDpars );
    //err |= ocl.finishRaw();  
    OCL_checkError(err, "evalVFs()");
    return F2max;
}

void unpack_system(  int isys, MMFFsp3_loc& ff, bool bForces=false, bool bVel=false ){
    //rintf("MolWorld_sp3_multi::unpack_system(%i) \n", isys);
    int i0n = isys * ocl.nnode;
    int i0a = isys * ocl.nAtoms;
    int i0v = isys * ocl.nvecs;
    unpack( ff.nvecs, ff.apos, atoms+i0v );
    if(bForces){ unpack( ff.nvecs, ff.fapos, aforces+i0v ); }
    if(bVel   ){ unpack( ff.nvecs, ff.vapos, avel   +i0v ); }
}

void upload(  bool bParams=false, bool bForces=0, bool bVel=true ){
    //printf("MolWorld_sp3_multi::upload() \n");
    int err=0;
    err|= ocl.upload( ocl.ibuff_atoms,  atoms  );
    err|= ocl.upload( ocl.ibuff_constr, constr );
    if(bForces){ err|= ocl.upload( ocl.ibuff_aforces, aforces ); }
    if(bVel   ){ err|= ocl.upload( ocl.ibuff_avel,    avel    ); }
    if(bParams){
        err|= ocl.upload( ocl.ibuff_lvecs,   lvecs );
        err|= ocl.upload( ocl.ibuff_ilvecs, ilvecs );
        err|= ocl.upload( ocl.ibuff_pbcshifts, pbcshifts );
        err|= ocl.upload( ocl.ibuff_neighs,     neighs    );
        err|= ocl.upload( ocl.ibuff_neighCell,  neighCell );
        err|= ocl.upload( ocl.ibuff_bkNeighs,   bkNeighs  );
        err|= ocl.upload( ocl.ibuff_bkNeighs_new,   bkNeighs_new  );
        //err|= ocl.upload( ocl.ibuff_neighForce, neighForce  );
        err|= ocl.upload( ocl.ibuff_REQs,   REQs   );
        err|= ocl.upload( ocl.ibuff_MMpars, MMpars );
        err|= ocl.upload( ocl.ibuff_BLs, BLs );
        err|= ocl.upload( ocl.ibuff_BKs, BKs );
        err|= ocl.upload( ocl.ibuff_Ksp, Ksp );
        err|= ocl.upload( ocl.ibuff_Kpp, Kpp );
    }
    err |= ocl.finishRaw(); 
    OCL_checkError(err, "MolWorld_sp2_multi::upload().finish");
}

void download( bool bForces=false, bool bVel=false ){
    int err=0;
    //printf("MolWorld_sp3_multi::download() \n");
    ocl.download( ocl.ibuff_atoms, atoms );
    if(bForces){ ocl.download( ocl.ibuff_aforces, aforces ); }
    if(bVel   ){ ocl.download( ocl.ibuff_avel,    avel    ); }
    err |= ocl.finishRaw(); 
    OCL_checkError(err, "MolWorld_sp2_multi::upload().finish");
}

void printMSystem( int isys, bool blvec=true, bool bilvec=false, bool bNg=true, bool bNgC=true, bool bPos=true, bool bF=false, bool bV=false ){
    int i0n = isys * ocl.nnode;
    int i0a = isys * ocl.nAtoms;
    int i0v = isys * ocl.nvecs;
    printf( "### System[%i] \n", isys );
    if(blvec ){printf("lvec\n");  printMat(lvecs [isys]); }
    if(bilvec){printf("ilvec\n"); printMat(ilvecs[isys]); }
    for(int i=0; i<ocl.nAtoms; i++){
        printf("[%i]", i );
        if(bNg )printf("ng (%3i,%3i,%3i,%3i)",  neighs   [i0a+i].x, neighs   [i0a+i].y, neighs   [i0a+i].z, neighs   [i0a+i].w );
        if(bNgC)printf("ngC(%2i,%2i,%2i,%2i)",  neighCell[i0a+i].x, neighCell[i0a+i].y, neighCell[i0a+i].z, neighCell[i0a+i].w );
        if(bPos)printf("pa(%6.3f,%6.3f,%6.3f)", atoms    [i0v+i].x, atoms    [i0v+i].y, atoms    [i0v+i].z );
        if(bF  )printf("va(%6.3f,%6.3f,%6.3f)", aforces  [i0v+i].x, aforces  [i0v+i].y, aforces  [i0v+i].z );
        if(bV  )printf("fa(%6.3f,%6.3f,%6.3f)", avel     [i0v+i].x, avel     [i0v+i].y, avel     [i0v+i].z );
        printf("\n" );
    }
}



// ===============================================
//       Implement    MultiSolverInterface
// ===============================================

virtual int paralel_size( )override{ return nSystems; }

virtual double solve_multi ( int nmax, double tol )override{
    return eval_MMFFf4_ocl( nmax, tol );
}

virtual void setGeom( int isys, Vec3d* ps, Mat3d *lvec, bool bPrepared )override{
    int i0n = isys * ocl.nnode;
    int i0a = isys * ocl.nAtoms;
    int i0v = isys * ocl.nvecs;


    // we cannot pack directly ps because  pipos is not defined in setGeom
    for(int i=0; i<ffl.natoms; i++){ ffl.apos[i]=ps[i]; } 
    ffl.initPi();         
    pack( ocl.nvecs, ffl.apos, atoms+i0v );

    if(lvec){ ffl.setLvec(*lvec); };
    Mat3_to_cl( ffl.lvec,    lvecs [isys] );
    Mat3_to_cl( ffl.invLvec, ilvecs[isys] );
    printf( "setGeom[%i] lvec{%6.3f,%6.3f,%6.3f}{%6.3f,%6.3f,%6.3f}{%6.3f,%6.3f,%6.3f}\n",  isys,  ffl.lvec.a.x,ffl.lvec.a.y,ffl.lvec.a.z,   ffl.lvec.b.x,ffl.lvec.b.y,ffl.lvec.b.z,   ffl.lvec.c.x,ffl.lvec.c.y,ffl.lvec.c.z  ); 
    // if( ! bPrepared ){
    //     set ( ocl.nvecs,     avel+i0v  );
    //     // ---- upload to GPU
    //     ocl.upload( ocl.ibuff_atoms,   atoms,   ocl.nvecs, i0v  );
    //     ocl.upload( ocl.ibuff_avel,    avel,    ocl.nvecs, i0v  );
    //     ocl.upload( ocl.ibuff_lvecs,   lvecs,   1,         isys );
    //     ocl.upload( ocl.ibuff_ilvecs,  ilvecs,  1,         isys );
    //     //for(int i=0; i<ocl.nAtoms; i++){ Quat4f a=atoms[i+i0v]; a.w=-1.0; constr[i+i0a] = a; }  // contrains
    // }
}

virtual double getGeom     ( int isys, Vec3d* ps, Mat3d *lvec, bool bPrepared )override{
    int i0n = isys * ocl.nnode;
    int i0a = isys * ocl.nAtoms;
    int i0v = isys * ocl.nvecs;
    //if( ! bPrepared ) ocl.download( ocl.ibuff_atoms,   atoms,   ocl.nvecs, i0v  );
    pack( ocl.nvecs, ps, atoms+i0v );
    //if(lvec){ lvec.a lvecs };
    return 0;
}

virtual void downloadPop()override{
    ocl.download( ocl.ibuff_atoms,  atoms  );
}

virtual void uploadPop  ()override{
    int ntot = nSystems*ocl.nvecs;
    set ( ntot, avel);
    ocl.upload( ocl.ibuff_avel , avel   );
    ocl.upload( ocl.ibuff_atoms, atoms  );
    ocl.upload( ocl.ibuff_lvecs,   lvecs  );
    ocl.upload( ocl.ibuff_ilvecs,  ilvecs );
    //ocl.upload( ocl.ibuff_constr, constr );
}

virtual void upload_pop( const char* fname ){
    printf("MolWorld_sp3::upload_pop(%s)\n", fname );
    gopt.loadPopXYZ( fname );
    int nmult=gopt.population.size(); 
    int npara=paralel_size(); if( nmult!=npara ){ printf("ERROR in GlobalOptimizer::lattice_scan_1d_multi(): (imax-imin)=(%i) != solver.paralel_size(%i) => Exit() \n", nmult, npara ); exit(0); }
    gopt.upload_multi(nmult,0);

    for(int i=0;i<nSystems; i++){ printMSystem(i); }

    /*
    //int initMode=1;
    int initMode=0;
    //int nstesp = 40;
    int nstesp = 2;
    Mat3d dlvec =  Mat3d{   0.2,0.0,0.0,    0.0,0.0,0.0,    0.0,0.0,0.0  };
    gopt.lattice_scan_2d_multi( nstesp, dlvec, initMode, "lattice_scan_2d_multi.xyz" );
    */
}

virtual void setSystemReplica (int i){ 
    int err=0;
    iSystemCur = i;   
    int i0v = iSystemCur * ocl.nvecs;
    int i0a = iSystemCur * ocl.nAtoms;
    ocl.download( ocl.ibuff_atoms,   atoms,  ocl.nvecs, i0v );
    //ocl.download( ocl.ibuff_aforces, aforces, ocl.nvecs, i0v );
    //ocl.download( ocl.ibuff_avel,    avel,    ocl.nvecs, i0v );
    err|=ocl.finishRaw();  OCL_checkError(err, "setSystemReplica()");
    unpack_system(iSystemCur, ffl, true, true); 

    copy( ffl.natoms, neighs+i0a, ffl.neighs );
    //for(int j=0; j<ffl.natoms; j++){ printf( "setSystemReplica[%i] ffl.neighs[%i](%3i,%3i,%3i,%3i)\n", i,  j,  ffl.neighs[j].x,ffl.neighs[j].y,ffl.neighs[j].z,ffl.neighs[j].w ); };

    Mat3_from_cl( builder.lvec, lvecs[iSystemCur] );
    ffl.setLvec( builder.lvec );
    makePBCshifts( nPBC, ffl.lvec );
}
virtual int countSystemReplica(     ){ return nSystems; }

void saveDebugXYZreplicas( int itr, double F ){
    int err=0;
    char fname  [100];
    char commnet[100];
    ocl.download( ocl.ibuff_atoms   , atoms   );
    err|=ocl.finishRaw();  OCL_checkError(err, "saveDebugXYZreplicas()");
    for(int isys=0; isys<nSystems; isys++){
        sprintf( fname  , "DEBUG_sys%03i.xyz", isys );
        sprintf( commnet, "# itr %i |F| %g ", itr, F );
        unpack_system( isys, ffl, false, false );
        if( ckeckNaN_d( ffl.nvecs, 3, (double*)ffl.apos, fname, false ) ){ printf( "saveDebugXYZreplicas(itr=%i) NaNs in replica[%i] => Exit() \n", itr, isys  ); exit(0); };
        saveXYZ( fname, commnet, false, "a" );
    }
}

// ==================================
//          SETUP OCL KERNELS 
// ==================================

void setup_MMFFf4_ocl(){
    printf("MolWorld_sp3_multi::setup_MMFFf4_ocl() \n");
    ocl.nDOFs.x=ff.natoms;
    ocl.nDOFs.y=ff.nnode;
    if(!task_cleanF)   task_cleanF = ocl.setup_cleanForceMMFFf4( ff4.natoms, ff4.nnode       );
    if(!task_move  )   task_move  = ocl.setup_updateAtomsMMFFf4( ff4.natoms, ff4.nnode       );
    if(!task_print )   task_print = ocl.setup_printOnGPU       ( ff4.natoms, ff4.nnode       );
    if(!task_MMFF  )   task_MMFF  = ocl.setup_getMMFFf4        ( ff4.natoms, ff4.nnode, bPBC );

    Vec3i nPBC_=nPBC; if(!bPBC){ nPBC_=Vec3iZero; }; printf( "MolWorld_sp3_multi::setup_MMFFf4_ocl() bPBC=%i nPBC(%i,%i,%i) n", bPBC, nPBC_.x,nPBC_.y,nPBC_.z );
    if((!task_NBFF_Grid)&&bGridFF ){ task_NBFF_Grid = ocl.setup_getNonBond_GridFF( ff4.natoms, ff4.nnode, nPBC_ ); } 
    if(!task_NBFF                 ){ task_NBFF      = ocl.setup_getNonBond       ( ff4.natoms, ff4.nnode, nPBC_ ); }
    //exit(0);

    // if(!task_NBFF  ) { 
    //     if( bGridFF ){ task_NBFF  = ocl.setup_getNonBond_GridFF( ff4.natoms, ff4.nnode, nPBC ); } 
    //     else         { task_NBFF  = ocl.setup_getNonBond       ( ff4.natoms, ff4.nnode, nPBC ); }
    // }
    if(!task_cleanF)task_cleanF = ocl.setup_cleanForceMMFFf4 ( ff4.natoms, ff4.nnode       );
}

void setup_NBFF_ocl(){
    printf("MolWorld_sp3_multi::setup_NBFF_ocl() \n");
    ocl.nDOFs.x=ff.natoms;
    ocl.nDOFs.y=ff.nnode;
    if(!task_cleanF )task_cleanF = ocl.setup_cleanForceMMFFf4 ( ff4.natoms, ff4.nnode        );
    if(!task_move   )task_move   = ocl.setup_updateAtomsMMFFf4( ff4.natoms, ff4.nnode        ); 
    if(!task_print  )task_print  = ocl.setup_printOnGPU       ( ff4.natoms, ff4.nnode        );
    if(!task_NBFF   )task_NBFF   = ocl.setup_getNonBond       ( ff4.natoms, ff4.nnode, nPBC  );
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

int eval_MMFFf4_ocl( int niter, double Fconv=1e-6, bool bForce=false ){ 
    //printf("MolWorld_sp3_multi::eval_MMFFf4_ocl() niter=%i \n", niter );
    //for(int i=0;i<npbc;i++){ printf( "CPU ipbc %i shift(%7.3g,%7.3g,%7.3g)\n", i, pbc_shifts[i].x,pbc_shifts[i].y,pbc_shifts[i].z ); }
    double F2conv = Fconv*Fconv;
    picked2GPU( ipicked,  1.0 );
    int err=0;
    if( task_MMFF    ==0 )setup_MMFFf4_ocl();
    //if( task_MMFFloc ==0 )task_MMFFloc=ocl.setup_evalMMFFf4_local( niter );
    // evaluate on GPU
    //long T0 = getCPUticks();
    //if(task_MMFFloc){
    //    task_MMFFloc->enque_raw();

    //niter=10;

    int nPerVFs = _min(10,niter);
    int nVFs    = niter/nPerVFs;


    long T0 = getCPUticks();
    int niterdone=0;
    if( itest != 0 ){
        //niter=1;
        //niter=10;
        niter=100;
        if(itest==1){ 
            if( task_MMFFloc1 ==0 )task_MMFFloc1=ocl.setup_evalMMFFf4_local1( niter );
            task_MMFFloc1     ->enque_raw(); 
        }else if(itest==2){
            if( task_MMFFloc2 ==0 )task_MMFFloc2=ocl.setup_evalMMFFf4_local2( niter );
            task_MMFFloc2     ->enque_raw(); 
        }else if(itest==3){
            if( task_MMFFloc_test==0  )task_MMFFloc_test=ocl.setup_evalMMFFf4_local_test( niter );
            task_MMFFloc_test->enque_raw();
        }
        //exit(0);    
    }else for(int i=0; i<nVFs; i++){
        //long T0 = getCPUticks();
        for(int j=0; j<nPerVFs; j++){
            err |= task_cleanF->enque_raw();      // this should be solved inside  task_move->enque_raw();   if we do not need to output force 
            if(bGridFF){ err |= task_NBFF_Grid ->enque_raw(); }
            else       { err |= task_NBFF      ->enque_raw(); }
            err |= task_MMFF->enque_raw();
            //ocl.printOnGPU( iSystemCur, int4{1,1,0,0}  );
            //OCL_checkError(err, "eval_MMFFf4_ocl.task_NBFF_Grid");
            //err |= task_print   ->enque_raw();    // just printing the forces before assempling
            err |= task_move      ->enque_raw(); 
            //OCL_checkError(err, "eval_MMFFf4_ocl_1");
            niterdone++;

            { //DEBUG
                download( true, true );
                err|=ocl.finishRaw();  OCL_checkError(err, "eval_MMFFf4_ocl().DEBUG.download");
                for(int isys=0; isys<nSystems; isys++){
                    unpack_system( isys, ffl, false, false );

                    double frange=100.0;
                    bool berr=false;
                    berr|= ckeckRange( ffl.nvecs, 3, (double*)ffl.apos,  -frange, frange, "apos",  true );
                    berr|= ckeckRange( ffl.nvecs, 3, (double*)ffl.fapos, -frange, frange, "fapos", true );
                    berr|= ckeckRange( ffl.nvecs, 3, (double*)ffl.vapos, -frange, frange, "vapos", true );
                    if(berr){ printf( "ERROR eval_MMFFf4_ocl().DEBUG outOfRange(%g) in replica[%i].apos  => Exit() \n", frange, isys ); exit(0); };

                    berr=false;
                    berr|= ckeckNaN_d( ffl.nvecs, 3, (double*)ffl.apos,  "apos",  true );
                    berr|= ckeckNaN_d( ffl.nvecs, 3, (double*)ffl.fapos, "fapos", true );
                    berr|= ckeckNaN_d( ffl.nvecs, 3, (double*)ffl.vapos, "vapos", true );
                    if(berr){ printf( "ERROR eval_MMFFf4_ocl().DEBUG NaNs in replica[%i].apos  => Exit() \n", isys ); exit(0); };
                }
                //saveDebugXYZreplicas( niterdone, 0.0 );
            }

        }
        err |= ocl.finishRaw();  OCL_checkError(err, "evalVFs()");
        double F2 = evalVFs();
        //saveDebugXYZreplicas( niterdone, sqrt(F2) );
        // if( F2<F2conv  ){ 
        //     printf( "eval_MMFFf4_ocl() all %i replicas relaxed in %i steps, |F|(%g)<%g \n", nSystems, niterdone, sqrt(F2), Fconv ); 
        //     return niterdone; 
        // }
        
        //long T2 =  getCPUticks();
        //printf("eval_MMFFf4_ocl(),evalVFs() time.tot=%g[ms] time.download=%g[ms] niter=%i \n", ( T2-T0 )*tick2second*1000, ( T2-T1)*tick2second*1000, niter );
        //printf("eval_MMFFf4_ocl(),evalVFs() time=%g[ms] niter=%i \n", ( getCPUticks()-T0 )*tick2second*1000 , niter );
        //fire[iSystemCur].print();
    }
    err |= ocl.finishRaw(); printf("eval_MMFFf4_ocl() time=%7.3f[ms] niter=%i \n", ( getCPUticks()-T0 )*tick2second*1000 , niterdone );

    // ===============================================================================================================================================================================
    //     Performance Measurements ( 10 replicas of polymer-2_new )
    //     for niter = 100         
    //         time=  5.0018 [ms]   task_cleanF; task_MMFF;            task_move;
    //         time= 77.3187 [ms]   task_cleanF; task_MMFF; task_NBFF; task_move;
    //   => Need to optimize task_NBFF ... i.e.  kernel getNonBond(), 
    //          1) first try by using work_group_size = get_local_size(0) = 32  =     __local float4 LATOMS[32],LCLJS [32];
    //              * with task->local.x=32   time=14.4335[ms]   but system explodes 
    //              * with nsys=50 rather than 10 time goes to time=28.9872[ms] rather than time=14.4335[ms]
    //              * with task->local.x=64 and LATOMS[64],LCLJS[64]; we fit the whole system on one work_group => system does not explode !!!!!!! :DDDD
    //          2) then try to optimize inner most loop over pbc_shifts 
    //          3) remove IF condition for vdw ? ( better use force limit )
    // ===============================================================================================================================================================================


    return niterdone;
}

double eval_NBFF_ocl( int niter, bool bForce=false ){ 
    //printf("MolWorld_sp3_multi::eval_NBFF_ocl() \n");
    int err=0;
    if( task_NBFF==0 ){ setup_NBFF_ocl(); }
    // evaluate on GPU
    for(int i=0; i<niter; i++){
        //err |= task_cleanF->enque_raw(); // this should be solved inside  task_move->enque_raw();
        err |= task_NBFF  ->enque_raw();
        //err |= task_print ->enque_raw(); // just printing the forces before assempling
        err |= task_move  ->enque_raw();
        //OCL_checkError(err, "eval_NBFF_ocl.1");
    }
    if(bForce)ocl.download( ocl.ibuff_aforces, ff4.fapos, ff4.nvecs, ff4.nvecs*iSystemCur );
    ocl.download( ocl.ibuff_atoms,   ff4.apos , ff4.nvecs, ff4.nvecs*iSystemCur );
    err |=  ocl.finishRaw();
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
        bool err=0;
        //E += ff .eval();
        E += ffl.eval();
        //E += eval_f4();
        //printf( "atom[0] nbmol(%g,%g,%g) ff(%g,%g,%g) ffl(%g,%g,%g) \n", nbmol.apos[0].x,nbmol.apos[0].y,nbmol.apos[0].z,  ff.apos[0].x,ff.apos[0].y,ff.apos[0].z,  ffl.apos[0].x,ffl.apos[0].y,ffl.apos[0].z );
        
        unpack_system(iSystemCur, ffl, true, true); 

        /*
        bool bForce = true;
        if(bForce)ocl.download( ocl.ibuff_aforces, ff4.fapos, ff4.nvecs, ff4.nvecs*iSystemCur );
        ocl          .download( ocl.ibuff_atoms,   ff4.apos , ff4.nvecs, ff4.nvecs*iSystemCur );
        //for(int i=0; i<ff4.nvecs; i++){  printf("OCL[%4i] f(%10.5f,%10.5f,%10.5f) p(%10.5f,%10.5f,%10.5f) pi %i \n", i, ff4.fapos[i].x,ff4.fapos[i].y,ff4.fapos[i].z,  ff4.apos[i].x,ff4.apos[i].y,ff4.apos[i].z,  i>=ff4.natoms ); }
        err |= ocl.finishRaw(); //printf("eval_MMFFf4_ocl() time=%g[ms] niter=%i \n", ( getCPUticks()-T0 )*tick2second*1000 , niter );
        OCL_checkError(err, "eval_MMFFf4_ocl");
        unpack( ff4.natoms, ffl.  apos, ff4.  apos );
        unpack( ff4.nnode,  ffl. pipos, ff4. pipos );
        if(bForce){
            unpack( ff4.natoms, ffl. fapos, ff4. fapos );
            unpack( ff4.nnode,  ffl.fpipos, ff4.fpipos );
            // ---- Check Invariatns   - this works only if we unpack forces
            //fcog  = sum ( ffl.natoms, ffl.fapos   );
            //tqcog = torq( ffl.natoms, ffl.apos, ffl.fapos );
            //if(  fcog.norm2()>1e-8 ){ printf("WARRNING: eval_MMFFf4_ocl |fcog| =%g; fcog=(%g,%g,%g)\n", fcog.norm(),  fcog.x, fcog.y, fcog.z ); exit(0); }
            //if( tqcog.norm2()>1e-8 ){ printf("WARRNING: eval_MMFFf4_ocl |torq| =%g; torq=(%g,%g,%g)\n", tqcog.norm(),tqcog.x,tqcog.y,tqcog.z ); exit(0); }   // NOTE: torq is non-zero because pi-orbs have inertia
        }
        */

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
        //else        { E+= nbmol .evalMorse   ( surf, false,                  gridFF.alphaMorse, gridFF.Rdamp );  }
        else          { E+= nbmol .evalMorsePBC( surf, gridFF.grid.cell, nPBC, gridFF.alphaMorse, gridFF.Rdamp );  }
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
        //nIter = 100;
        nIter = 1;
        eval_MMFFf4_ocl( nIter );
        unpack_system(iSystemCur, ffl, true, true); 
        //eval_NBFF_ocl  ( 1 ); 
        //eval_NBFF_ocl_debug(1); //exit(0);
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

virtual void scanSurfFF( int n, Quat4f* ps, Quat4f* REQs, Quat4f* fs )override{
    ocl.sampleGridFF( n, fs, ps, REQs, true  );
}

bool checkSampleGridFF( int n, Vec3d p0, Vec3d p1, Quat4d REQ=Quat4d{ 1.487, 0.02609214441, +0.1, 0.}, double tol=1e-2, bool bExit=false, bool bPrint=false, bool bWarn=true, const char* logfiflename="checkSampleGridFF.log" ){
    if(bPrint){ printf("MolWorld_sp3_multi::checkSampleGridFF(np=%i,p0{%6.3f,%6.3f,%6.3f},p1{%6.3f,%6.3f,%6.3f}REQ{%6.3f,%10.7f,%6.3f,%10.7f}) \n", n, ocl.nAtoms, p0.x,p0.y,p0.z,  p1.x,p1.y,p1.z, REQ.x,REQ.y,REQ.z,REQ.w ); };
    //if((ocl.nAtoms*ocl.nSystems)<n){ printf("ERROR in MolWorld_sp3_multi::checkSampleGridFF() natom(%i)<n(%i) => Exit()\n", ocl.nAtoms*ocl.nSystems, n ); exit(0); } // This is no-longer needed
    Quat4f* samp_fs   = new Quat4f[n]; 
    Quat4f* samp_ps   = new Quat4f[n];
    Quat4f* samp_REQs = new Quat4f[n];
    Vec3d dp = (p1-p0)*(1./n);
    for(int i=0; i<n; i++){
        v2f4(p0+dp*i, *((float4*)samp_ps+i) );  
        samp_REQs [i] = (Quat4f)REQ;
        //samp_fs =  
        //samp_ps   = 
        //samp_REQs    = 
        //printf( "checkSampleGridFF[%i] p(%g,%g,%g) REQ(%g,%g,%g,%g)\n", i,  atoms[i].x,atoms[i].y,atoms[i].z,    REQs[i].x,REQs[i].y,REQs[i].z,REQs[i].w  );
    }
    ocl.grid_shift0.f=(Vec3f)gridFF.shift0;
    ocl.sampleGridFF( n, samp_fs, samp_ps, samp_REQs, true   );
    FILE * logf=0;
    if(logfiflename){ 
        logf = fopen(logfiflename,"w");
        fprintf( logf, "GridFF::checkEFProfileVsNBFF(np=%i,natoms=%i,npbc=%i,p2{%6.3f,%6.3f,%6.3f},p1{,%6.3f,%6.3f,%6.3f}REQ{%g,%g,%g,%g}) \n", n, gridFF.natoms,gridFF.npbc, p0.x,p0.y,p0.z,  p1.x,p1.y,p1.z, REQ.x,REQ.y,REQ.z,REQ.w );
        fprintf(     logf, "i   x y z     E  Eref     fx fx_ref      fy fy_ref     fz  fz_ref\n");
    }
    if(bPrint){     printf("i   x y z     E  Eref     fx fx_ref      fy fy_ref     fz  fz_ref\n"); }
    double tol2=tol*tol;
    //Quat4f PLQ = REQ2PLQ( REQ, alphaMorse );   //printf( "PLQ %6.3f %10.7f %6.3f \n", PLQ.x,PLQ.y,PLQ.z   );
    //bool err = false;
    bool bErr=false;
    double err2Max=0;
    int imax = -1;
    for(int i=0; i<n; i++){
        Vec3d fref;
        Vec3d  p    = p0 + dp*i;
        Quat4f fe   = samp_fs[i];
        float  Eref = gridFF.getMorseQH_PBC_omp( p, REQ, fref );
        float  dE   = fe.e-Eref;
        Vec3f  df   = fe.f - (Vec3f)fref;
        double e2e = dE*dE     /( Eref*Eref + fe.e*fe.e        + 1 );         
        double e2f = df.norm2()/( fe.f.norm2() + fref.norm2() + 1 ); 
        if( e2e>err2Max ){ err2Max=e2e; imax=i; }
        if( e2f>err2Max ){ err2Max=e2f; imax=i; }
        if( (e2e>tol2) || (e2f>tol2) ){
            bErr=true;
            if(bWarn)printf("WARRNING[%i/%i] dE=%g |dF|=%g p(%6.3f,%6.3f,%6.3f) GridFF(%g,%g,%g|%g)  NBFF(%g.%g,%g|%g)\n", i, n, dE, df.norm(), p.x,p.y,p.z,   fe.x,fe.y,fe.z,fe.w,   fref.x,fref.y,fref.z,Eref  );
            if(bExit){ printf("ERROR in GridFF::checkEFProfileVsNBFF() - GridFF force does not match NBFF reference at test point %i MaxRelativeError=%g => Exit()\n", i, sqrt(err2Max) ); exit(0); }
        } 
        if(bPrint){ printf(       "%i    %6.3f %6.3f %6.3f    %g %g   %g %g    %g %g    %g %g\n", i, p.x, p.y, p.z,    fe.e, Eref,    fe.x,fref.x,    fe.y,fref.y,    fe.z,fref.z ); }
        if(logf  ){fprintf( logf, "%3i    %6.3f %6.3f %6.3f    %14.6f %14.6f   %14.6f %14.6f    %14.6f %14.6f    %14.6f %14.6f\n", i, p.x, p.y, p.z,    fe.e, Eref,    fe.x,fref.x,    fe.y,fref.y,    fe.z,fref.z ); }
    }
    if(logf){ fclose(logf); }
    if(bWarn && bErr ){
        //printf("WARRNING GridFF MaxRelativeError=%g at point[%i]\n",  sqrt(err2Max), imax );
        Vec3d fref;
        Vec3d  p    = p0 + dp*imax;
        Quat4f fe   = samp_fs[imax];
        float  Eref = gridFF.getMorseQH_PBC_omp( p, REQ, fref );
        float  dE   = fe.e-Eref;
        Vec3f  df   = fe.f - (Vec3f)fref;
        double e2e = dE*dE     /(Eref*Eref + fe.e*fe.e + 1);          err2Max=fmax(err2Max,e2e);
        double e2f = df.norm2()/( fe.f.norm2() + fref.norm2() + 1 );  err2Max=fmax(err2Max,e2f);
        printf("WARRNING GridFF MaxError=%g at[%i/%i] dE=%g |dF|=%g p(%6.3f,%6.3f,%6.3f) GridFF(%g.%g,%g|%g)  NBFF(%g.%g,%g|%g)\n",  sqrt(err2Max), imax, n, dE, df.norm(), p.x,p.y,p.z, fe.x,fe.y,fe.z,fe.w,   fref.x,fref.y,fref.z,Eref  );
    }
    delete [] samp_fs; 
    delete [] samp_ps;
    delete [] samp_REQs;
    return bErr;
}

bool evalCheckGridFF_ocl( int imin=0, int imax=1, bool bExit=true, bool bPrint=true, double tol=1e-2, Quat4d REQ=Quat4d{ 1.487, 0.02609214441, +0.1, 0.}, double dz=0.05 ){
    REQ=Quat4d{ 1.487, 0.02609214441,-0.1, 0.};
    printf( "MolWorld_sp3::evalCheckGridFF_ocl() natoms=%i npbc=%i apos=%li REQs=%li shifts=%li \n", gridFF.natoms, gridFF.npbc, gridFF.apos, gridFF.REQs, gridFF.shifts );
    _checkNull(gridFF.shifts)
    _checkNull(gridFF.REQs)
    _checkNull(gridFF.apos)
    bool err = false;
    double zmin=gridFF.grid.pos0.z+gridFF.shift0.z+1.0;
    double zmax=gridFF.grid.pos0.z+gridFF.grid.cell.c.z*0.5;
    zmax=fmin(zmax,10);
    int nz = ((int)((zmax-zmin)/dz)) + 1;
    for( int ia=imin; ia<imax; ia++ ){
        Vec3d p0=gridFF.apos[ia]; p0.z=zmin;
        Vec3d p1=p0;       p1.z=zmax;
        //double err =  checkEFProfileVsNBFF( n, p0, p1, REQ, tol, false,false );
        err |= checkSampleGridFF( nz, p0,p1, REQ, tol, false, false, true );
        if(err>tol){    
            //if(bPrint)checkEFProfileVsNBFF( n, p0, p1, REQ, tol, false, true, false );
            if(bPrint)checkSampleGridFF   ( nz, p0,p1, REQ, tol, false, true, false );
            if(bExit){ printf("ERROR in GridFF::checkZProfilesOverAtom(%i) - GridFF force does not match NBFF reference, MaxRelativeError=%g => Exit()\n", ia, err ); exit(0); }
            break;
        }
        //checkEFProfileVsNBFF( n, p0, p1, REQ, tol,  false, true, false );
        //double err =  checkEFProfileVsNBFF( n, p0, p1, REQ, tol, false,false );
        //err |= checkZProfilesOverAtom( ia, nz, zmin, zmax, REQ, tol, bExit, bPrint );
    }
    return err;
}

void surf2ocl(Vec3i nPBC, bool bSaveDebug=false){
    int err=0;
    printf( "surf2ocl() gridFF.natoms=%i nPBC(%i,%i,%i)\n", gridFF.natoms, nPBC.x,nPBC.y,nPBC.z );
    long T0=getCPUticks();
    Quat4f* atoms_surf = new Quat4f[gridFF.natoms];
    Quat4f* REQs_surf  = new Quat4f[gridFF.natoms];
    pack( gridFF.natoms, gridFF.apos,  atoms_surf, sq(gridFF.Rdamp) );
    pack( gridFF.natoms, gridFF.REQs, REQs_surf                     );   // ToDo: H-bonds should be here
    long T1=getCPUticks();
    ocl.GFFparams.x = gridFF.Rdamp;
    ocl.GFFparams.y = gridFF.alphaMorse;
    //v2f4( gridFF.grid.pos0, ocl.grid_p0 );
    ocl.grid_p0.f     = (Vec3f)gridFF.grid.pos0;
    ocl.grid_shift0.f = (Vec3f)gridFF.shift0;
    ocl.makeGridFF( gridFF.grid, nPBC, gridFF.natoms, (float4*)atoms_surf, (float4*)REQs_surf, true );
    err |=  ocl.finishRaw();    OCL_checkError(err, "surf2ocl.makeGridFF.finish");
    //ocl.addDipoleField( gridFF.grid, (float4*)dipole_ps, (float4*), true );
    printf( ">>time(ocl.makeGridFF() %g \n", (getCPUticks()-T1)*tick2second );
    bool bDownload =true; // WARRNING : if I set this OFF it crashes unexpectedly
    //bSaveDebug=true;
    bSaveDebug=false; 
    if(bDownload){
        gridFF.allocateFFs();
        ocl.download( ocl.itex_FE_Paul, gridFF.FFPaul );  
        ocl.download( ocl.itex_FE_Lond, gridFF.FFLond );
        ocl.download( ocl.itex_FE_Coul, gridFF.FFelec );    
        err |=  ocl.finishRaw();    OCL_checkError(err, "surf2ocl.download.finish");
        printf( ">>time(surf2ocl.download() %g \n", (getCPUticks()-T1)*tick2second );
        if(bSaveDebug){ saveGridXsfDebug(); }
    }
    delete [] atoms_surf;
    delete [] REQs_surf;
    ocl.buffers[ocl.ibuff_atoms_surf].release();
    ocl.buffers[ocl.ibuff_REQs_surf ].release();
    err |=  ocl.finishRaw();    OCL_checkError(err, "surf2ocl.makeGridFF.atoms_surf.release");
    //exit(0);
}

virtual void initGridFF( const char * name, bool bGrid=true, bool bSaveDebugXSFs=false, double z0=NAN, Vec3d cel0={-0.5,-0.5,0.0}, bool bAutoNPBC=true )override{
    printf( "MolWorld_sp3_multi::initGridFF() \n");
    if(verbosity>0)printf("MolWorld_sp3_multi::initGridFF(%s,bGrid=%i,z0=%g,cel0={%g,%g,%g})\n",  name, bGrid, z0, cel0.x,cel0.y,cel0.z  );
    if(gridFF.grid.n.anyEqual(0)){ printf("ERROR in MolWorld_sp3_multi::initGridFF() zero grid.n(%i,%i,%i) => Exit() \n", gridFF.grid.n.x,gridFF.grid.n.y,gridFF.grid.n.z ); exit(0); };
    gridFF.grid.center_cell( cel0 );
    bGridFF=true;
    gridFF.bindSystem      (surf.natoms, surf.atypes, surf.apos, surf.REQs );
    gridFF.setAtomsSymetrized( gridFF.natoms, gridFF.atypes, gridFF.apos, gridFF.REQs, 0.1 );
    //gridFF.setAtomsSymetrized(surf.natoms, surf.atypes, surf.apos, surf.REQs );
    gridFF.evalCellDipole();
    if( ( fabs(gridFF.Q)>1e-6 ) || (gridFF.dip.norm2()>1e-8) ){ printf("ERROR: GridFF has dipole and dipole correction not yet implemented => exit() \n"); exit(0); }
    if( isnan(z0) ){ z0=gridFF.findTop();   if(verbosity>0) printf("GridFF::findTop() %g \n", z0);  };
    gridFF.grid.pos0.z=z0;
    gridFF.lvec = gridFF.grid.cell;
    //if(verbosity>1)
    //gridFF.grid.printCell();
    gridFF.nPBC=Vec3i{1,1,0};
    if(bAutoNPBC){ autoNPBC( gridFF.grid.cell, gridFF.nPBC, 20.0 ); }
    gridFF.makePBCshifts ( gridFF.nPBC, gridFF.lvec );
    long T0 = getCPUticks();
    //bSurfAtoms=false;
    //gridFF.shift0 = Vec3d{0.,0., 0.0};
    gridFF.shift0 = Vec3d{0.,0.,-2.0};
    //gridFF.shift0 = Vec3d{0.,0.,-6.0};
    surf2ocl( gridFF.nPBC, bSaveDebugXSFs );
    printf( ">>time(init_ocl;GridFF_ocl): %g [s] \n", (getCPUticks()-T0)*tick2second  );
    bGridFF   =true; 
    // evalCheckGridFF_ocl();   // here are not initialized buffers atoms.aforce,REQs, so it will crash.
}


virtual int getMultiSystemPointers( int*& M_neighs,  int*& M_neighCell, Quat4f*& M_apos, int& nvec ) override {
    nvec        = ocl.nvecs;
    M_neighs    = (int*)neighs;
    M_neighCell = (int*)neighCell;
    M_apos      = atoms;
    return nSystems;
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
    
    {  // --- evaluate on CPU
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
    { // --- evaluate on CPU
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
        err |= task_cleanF->enque_raw(); // this should be solved inside  task_move->enque_raw();
        err |= task_MMFF  ->enque_raw();
        err |= task_NBFF  ->enque_raw();
        err |= task_print ->enque_raw(); // just printing the forces before assempling
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
    
    //exit(0);
    return 0;
}

// ==============================================================
//                   eval_NBFF_ocl_debug    
// ==============================================================

double eval_NBFF_ocl_debug( int niter ){ 
    printf("MolWorld_sp3_multi::eval_NBFF_ocl_debug() \n");
    int err=0;
    if( task_NBFF==0 ){ setup_NBFF_ocl(); }

    { // --- evaluate on CPU
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
        err |= task_cleanF->enque_raw(); // this should be solved inside  task_move->enque_raw();
        err |= task_NBFF  ->enque_raw();
        err |= task_print ->enque_raw(); // just printing the forces before assempling
        err |= task_move  ->enque_raw();
    }
    
    ocl.download( ocl.ibuff_aforces, ff4.fapos, ff4.nvecs, ff4.nvecs*iSystemCur );
    ocl.download( ocl.ibuff_atoms,   ff4.apos , ff4.nvecs, ff4.nvecs*iSystemCur );
    err |=  ocl.finishRaw();
    OCL_checkError(err, "eval_NBFF_ocl_debug");

    for(int i=0; i<ff4.nvecs; i++){  printf("OCL[%4i] f(%10.5f,%10.5f,%10.5f) p(%10.5f,%10.5f,%10.5f) pi %i \n", i, ff4.fapos[i].x,ff4.fapos[i].y,ff4.fapos[i].z,  ff4.apos[i].x,ff4.apos[i].y,ff4.apos[i].z,  i>=ff4.natoms ); }

    if( ckeckNaN_f( ff4.nvecs, 4, (float*)ff4.fapos, "gpu.fapos" ) ){ printf("ERROR in MolWorld_sp3_multi::eval_NBFF_ocl_debug() fapos contain NaNs \n"); exit(0); };

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
    
    //exit(0);
    return 0;
}


}; // class MolWorld_sp3_ocl

#endif
