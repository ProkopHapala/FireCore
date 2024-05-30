
#ifndef MolWorld_sp3_ocl_h
#define MolWorld_sp3_ocl_h

#include "MMFFf4.h"
//#include "Molecule.h"
//#include "MMFFBuilder.h"
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
    int id;

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
        float f_len = sqrt(ff);
        float v_len = sqrt(vv);
        cos_vf      = vf    / ( f_len*v_len + ff_safety );
        //renorm_vf   = v_len / ( f_len       + ff_safety );
        renorm_vf   = v_len / ( f_len       + v_len );
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
        //printf( "FIRE[%i].update_params() |f|=%g,|v|=%g,cos_vf=%g,renorm_vf=%g,cf=%g,cv=%g,damping=%g,dt=%g,lastNeg=%i \n", id, f_len,v_len,cos_vf, renorm_vf, cf, cv, damping, dt, lastNeg );
    }

    void print(){ printf( "dt %g damping %g cos_vf %g cv %g cf %g \n", dt, damping, cos_vf, cv, cf ); }
    void bind_params( FIRE_setup* par_ ){ par=par_; dt=par->dt_max; damping=par->damp_max; lastNeg=0; };

    FIRE()=default;
    FIRE( FIRE_setup* par_ ){ bind_params(par_); }

};



// ======================================
// class:        MolWorld_sp3_ocl
// ======================================

class MolWorld_sp3_multi : public MolWorld_sp3, public MultiSolverInterface { public:
    OCL_MM     ocl;
    FIRE_setup fire_setup{0.1,0.1};

    //int nSystems    = 1;
    int iSysFMax=-1;
    //int iSystemCur  = 5;    // currently selected system replica
    //int iSystemCur  = 8;    // currently selected system replica
    bool bGPU_MMFF = true;

    Quat4f* atoms      =0;
    Quat4f* aforces    =0;
    Quat4f* avel       =0;
    Quat4f* cvfs       =0;

    FIRE*   fire       =0;  // FIRE-relaxation state
    Quat4f* MDpars     =0;  // Molecular dynamics params
    Quat4f* TDrive     =0;  // temperature and drived dynamics

    Quat4f* constr     =0;
    Quat4f* constrK    =0;
    cl_Mat3*  bboxes   =0;

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


    // ---- Groups of atoms
    Vec2i*  granges  = 0; // [ nSystems*nGroup ]    
    int*    a2g      = 0; // [ nSystems*nAtoms ]  
    int*    g2a      = 0; // [ nSystems*nAtoms ]       
    Quat4f* gforces  = 0; // [ nSystems*nGroup ]    
    Quat4f* gtorqs   = 0; // [ nSystems*nGroup ]   
    Quat4f* gcenters = 0; // [ nSystems*nGroup ] 

    Quat4f* gfws      = 0; // [ nSystems*nGroup ] 
    Quat4f* gups      = 0; // [ nSystems*nGroup ]    
    Quat4f* gweights  = 0; // [ nSystems*nGroup ]   
    Vec2f*  gfweights = 0; // [ nSystems*nGroup ] 



    cl_Mat3*  lvecs    =0;
    cl_Mat3* ilvecs    =0;
    Quat4f*   pbcshifts =0;

    OCLtask* task_cleanF=0;
    OCLtask* task_NBFF=0;
    OCLtask* task_NBFF_Grid=0;
    OCLtask* task_SurfAtoms=0;
    OCLtask* task_MMFF=0;
    OCLtask* task_move=0;
    OCLtask* task_print=0;

    OCLtask* task_GroupUpdate=0;
    OCLtask* task_GroupForce=0;

    OCLtask* task_MMFFloc =0;
    OCLtask* task_MMFFloc1=0;
    OCLtask* task_MMFFloc2=0;
    OCLtask* task_MMFFloc_test=0;

    bool bDONE_surf2ocl = false;

    DynamicOpt*   opts=0;
    MMFFsp3_loc*  ffls=0;

    GOpt*         gopts=0;

    MMFFf4     ff4;

    Quat4f* afm_ps=0;

    const char* uploadPopName=0;

virtual MolWorldVersion getMolWorldVersion() const override { return MolWorldVersion::GPU; };

// ==================================
//         Initialization
// ==================================

void realloc( int nSystems_ ){
    printf("MolWorld_sp3_multi::realloc() \n");
    nSystems=nSystems_;
    printf( "MolWorld_sp3_multi::realloc() Systems %i nAtoms %i nnode %i \n", nSystems, ffl.natoms,  ffl.nnode );
    ocl.initAtomsForces( nSystems, ffl.natoms,  ffl.nnode, npbc+1 );
    //printf( "MolWorld_sp3_multi::realloc() Systems %i nAtoms %i nnode %i nvecs %i \n", nSystems, ocl.nAtoms, ocl.nnode, ocl.nvecs );
    // --- dynamical
    _realloc ( atoms,     ocl.nvecs*nSystems  );
    _realloc0( aforces,   ocl.nvecs*nSystems  , Quat4fZero );
    _realloc0( avel,      ocl.nvecs*nSystems  , Quat4fZero );
    _realloc0( cvfs,      ocl.nvecs*nSystems  , Quat4fZero );
    _realloc0( constr,    ocl.nAtoms*nSystems , Quat4fOnes*-1. );
    _realloc0( constrK,   ocl.nAtoms*nSystems , Quat4fOnes*-1. );
    _realloc0( bboxes,   nSystems, cl_Mat3{cl_float4{-1e+8,-1e+8,-1e+8,-1e+8,},cl_float4{+1e+8,+1e+8,+1e+8,+1e+8,}, cl_float4{-1.,-1.,-1.,-1.} }   );
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
    _realloc0( TDrive,   nSystems, Quat4f{0.0,-1.0,0.0,0.0} );

    _realloc( pbcshifts, ocl.npbc*nSystems );

    _realloc( fire,      nSystems  );

    // ToDo : it may be good to bind buffer directly in p_cpu buffer inside   OCLsystem::newBuffer()
}

void initMultiCPU(int nSys){
    _realloc( opts, nSys );
    _realloc( ffls, nSys );
    _realloc( gopts,nSys );

    double dtopt=ff.optimalTimeStep();
    for(int isys=0; isys<nSys; isys++){
        ffls[isys].clone( ffl, true, true );
        ffls[isys].makePBCshifts( nPBC, true );
        ffls[isys].id=isys;

        ffls[isys].PLQs = ffl.PLQs;  // WARNING : maybe we should make own copy of PLQs for each system because we have own copy of REQs for each system

        // ----- optimizer
        opts[isys].bindOrAlloc( ffls[isys].nDOFs, ffls[isys].DOFs, 0, ffls[isys].fDOFs, 0 ); 
        ffls[isys].vapos=(Vec3d*)opts[isys].vel;    
        opts[isys].initOpt( dtopt );
        opts[isys].cleanVel();
        gopts[isys].copy( go );
        //gopts[isys].startExploring();
        gopts[isys].bExploring = false;
        //gopts[isys].print();
        //gopts[isys].constrs.printSizes();  // Debug
        //gopts[isys].constrs.printDrives( );
        gopts[isys].nExplore = 1000;
    }
}

virtual void stopExploring () override { printf("MolWorld_sp3_multi::stopExploring()\n");  go.bExploring=false; for(int isys=0; isys<nSystems; isys++){ gopts[isys].bExploring = false; }; };
virtual void startExploring() override { printf("MolWorld_sp3_multi::startExploring()\n"); go.startExploring(); for(int isys=0; isys<nSystems; isys++){ gopts[isys].startExploring();   }; };

virtual int getMultiConf( float* Fconvs , bool* bExplors )override{ 
    //printf("MolWorld_sp3_multi::getMultiConf()\n");
    for(int isys=0; isys<nSystems; isys++){
        Fconvs  [isys] = sqrt(fire[isys].ff);   
        bExplors[isys] = gopts[isys].bExploring;
    }
    return nSystems; 
};


virtual void init() override {
    int err = 0;
    printf("# ========== MolWorld_sp3_multi::init() START\n");
    gopt.msolver = this;
    int i_nvidia = ocl.print_devices(true);
    ocl.init(i_nvidia);
    ocl.makeKrenels_MM("common_resources/cl" );
    MolWorld_sp3::init();

    // initialization of ff4 is here because parrent MolWorld_sp3 does not contain ff4 anymore 
    builder.toMMFFf4( ff4, true, bEpairs );  //ff4.printAtomParams(); ff4.printBKneighs(); 
    ff4.flipPis( Vec3fOne );
    ff4.setLvec((Mat3f)builder.lvec);
    ff4.makeNeighCells  ( nPBC );

    // ----- init systems
    realloc( nSystems );
    //if(bGridFF) evalCheckGridFF_ocl();  // this must be after we make buffers but before we fill them
    float random_init = 0.5;

    fire_setup.setup(opt.dt_max,opt.damp_max);
    for(int i=0; i<nSystems; i++){
        pack_system( i, ffl, true, false, random_init );
        fire[i].bind_params( &fire_setup );
        fire[i].id=i;
    }
    initMultiCPU(nSystems);
    //for(int i=0;i<nSystems; i++){ printMSystem(i); }

    upload( true, false );
    //bGridFF=false;
    bOcl   =false;
    //bOcl   =true;
    //setup_MMFFf4_ocl();
    int4 mask{1,1,0,0};
    //ocl.printOnGPU( 0,mask );
    //ocl.printOnGPU( 1,mask );
    //ocl.printOnGPU( 2,mask );

    iParalelMax= 4;
    iParalelMin=-1;
    //iParalel=1;
    iParalel=3;
    //iParalel=iParalelMax;

    
    printf( "uploadPopName @ %li", uploadPopName );
    if( uploadPopName ){   printf( "!!!!!!!!!!!!\n UPLOADING POPULATION FROM FILE (%s)", uploadPopName );  upload_pop( uploadPopName ); }


    printf("# ========== MolWorld_sp3_multi::init() DONE\n");
}


virtual int getGroupPose( Quat4f*& gpos, Quat4f*& gfw, Quat4f*& gup ) override { 
    gpos = gcenters + iSystemCur*ocl.nvecs; 
    gfw  = gfws     + iSystemCur*ocl.nvecs;
    gup  = gups     + iSystemCur*ocl.nvecs;
    return ocl.nGroup; 
};

void groups2ocl( int isys, bool bForce=true, bool bPose=false, bool bWeights=true ){
    int nG = ocl.nGroup;
    if( nG != groups.groups.size() ){  printf("ERROR ocl.nGroup(%i)!=groups.size(%i) => exit() \n", ocl.nGroup, groups.groups.size() ); exit(0); }
    int i0g = isys*nG;
    int i0a = isys*ocl.nvecs;
    for(int ig=0; ig<nG; ig++){
        Group& g = groups.groups[ig];
        int i = ig + i0g;
        if(bPose){
            gcenters[i].f = (Vec3f)g.cog;
            gfws    [i].f = (Vec3f)g.fw;
            gups    [i].f = (Vec3f)g.up;
        }
        if(bForce){
            gforces [i].f = (Vec3f)g.force;
            gtorqs  [i].f = (Vec3f)g.torq;
        }
        if(bWeights){
            printf( "CPU [isys=%i,ig=%i]\n", isys, ig );
            for(int j=0; j<g.i0n.y; j++){
                int ia = groups.g2a[ j + g.i0n.x ];
                //int ia = g2a[ i0g + j ];
                int iia = i0a + ia;
                gweights [ iia ] = groups. weights[ia];
                gfweights[ iia ] = groups.fweights[ia];
                printf( "CPU gweights[%i](%g,%g,%g,%g) gfweights(%g,%g)\n", iia, gweights[iia].x,gweights[iia].y,gweights[iia].z,gweights[iia].w,  gfweights[iia].x,gfweights[iia].y );
            }
        }
    }
}


int init_groups(){
    printf("MolWorld_sp3_multi::init_groups()\n" );
    // int ngroup = -1;
    // printf("atom2group.size()==%i\n", atom2group.size() );
    // for(int i=0; i<atom2group.size(); i++){ 
    //     //printf("atom2group[%i]==%i\n", i, atom2group[i]);
    //     ngroup=_max(ngroup,atom2group[i]); 
    // }
    // ngroup+=1;
    int ngroup = groups.groups.size();

    if(ngroup<=0) return -1;

    // int2*   granges  = 0; // [ nSystems*nGroup ]    
    // int*    a2g      = 0; // [ nSystems*nAtoms ]  
    // int*    g2a      = 0; // [ nSystems*nAtoms ]       
    // Quat4f* gforces  = 0; // [ nSystems*nGroup ]    
    // Quat4f* gtorqs   = 0; // [ nSystems*nGroup ]   
    // Quat4f* gcenters = 0; // [ nSystems*nGroup ] 

    int nGroupTot = ngroup*ocl.nSystems;
    int nAtomTot  = ocl.nvecs*ocl.nSystems;
    _realloc0( gforces,  nGroupTot, Quat4fZero );
    _realloc0( gtorqs,   nGroupTot, Quat4fZero );
    //_realloc0( gtorqs,   nGroupTot, Quat4f{1.0f,0.0f,0.0f,0.0f} );

    _realloc0( gcenters, nGroupTot, Quat4fZero );
    _realloc0( granges,  nGroupTot, Vec2iZero  );
    _realloc0( gfws,     nGroupTot, Quat4fZero  );
    _realloc0( gups,     nGroupTot, Quat4fZero  );

    _realloc0( gweights, nAtomTot, Quat4fZero  );
    _realloc0( gfweights,nAtomTot, Vec2fOnes   );
    _realloc0( a2g, nAtomTot, -1 );
    _realloc0( g2a, nAtomTot, -1 );

    //exit(0);
    ocl.initAtomGroups( ngroup );
    ocl.setGroupMapping( &atom2group[0] );
    for(int isys=0; isys<nSystems; isys++){
        groups2ocl( isys, true, false, true );
        Mat3_to_cl( bbox, bboxes[isys] );        // ToDo: this may need to go to different place
    }
    int err=0;
    err|= ocl.upload( ocl.ibuff_gweights,  gweights  ); OCL_checkError(err, "init_groups.upload(gweights)");
    err|= ocl.upload( ocl.ibuff_gfweights, gfweights ); OCL_checkError(err, "init_groups.upload(gfweights)");

    err|= ocl.upload( ocl.ibuff_gforces, gforces );     OCL_checkError(err, "init_groups.upload(gforces)");
    err|= ocl.upload( ocl.ibuff_gtorqs,  gtorqs  );     OCL_checkError(err, "init_groups.upload(gtorqs)");

    err|= ocl.upload( ocl.ibuff_bboxes,  bboxes  );     OCL_checkError(err, "init_groups.upload(bboxes)");  // ToDo: this may need to go to different place
    
    return err;
}

virtual void pre_loop() override {
    //printf("MolWorld_sp3_multi::pre_loop()\n" );
    init_groups();

    
    for(int isys=0; isys<nSystems; isys++){
        for(int ia=0; ia<ffl.natoms; ia++ ){
            int iaa = isys * ocl.nAtoms + ia;
            constr [iaa] = (Quat4f)ffl.constr[ia];
            constrK[iaa] = (Quat4f)ffl.constrK[ia];
        }
    }

    // for(int ic : constrain_list ){
    //     for(int isys=0; isys<nSystems; isys++){
    //         int i0a   = isys * ocl.nAtoms;
    //         int i0v   = isys * ocl.nvecs;
    //         Quat4f a=atoms[ic+i0v]; a.w=1.0; constr[ic+i0a] = a;
    //     }
    // }
    // ocl.upload( ocl.ibuff_constr, constr );     //OCL_checkError(err, "init_groups.upload(constr)");
}

// ==================================
//         GPU <-> CPU I/O
// ==================================

void pack_system( int isys, MMFFsp3_loc& ff, bool bParams=0, bool bForces=false, bool bVel=false, bool blvec=true, float l_rnd=-1 ){
    //printf("MolWorld_sp3_multi::pack_system(%i) \n", isys);
    //ocl.nvecs;
    int i0n   = isys * ocl.nnode;
    int i0a   = isys * ocl.nAtoms;
    int i0v   = isys * ocl.nvecs;
    int i0pbc = isys*ocl.npbc;

    if(blvec){
        if(npbc==0){ pbcshifts[isys].f=Vec3fZero; };
        //evalPBCshifts( nPBC, ff.lvec, pbcshifts + i0pbc );
        //evalPBCshifts( nPBC, ff.lvec, pbc_shifts );  // This must be before pack(ff.apos)
        pack( npbc,  pbc_shifts, pbcshifts+i0pbc );
        //ffl.initPi(pbc_shifts);
        Mat3_to_cl( ff.   lvec,  lvecs[isys] );
        Mat3_to_cl( ff.invLvec, ilvecs[isys] );
    }
    
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

    //Mat3d lvec0=ff.lvec;  // DEBUG
    if(bForces){ pack( ff.nvecs, ff.fapos, aforces+i0v ); }
    //if(bVel   ){ pack( ff.nvecs, opt.vel,  avel   +i0v ); }
    if(bParams){


        //Mat3d lvec=lvec0; lvec.a.addRandomCube(1.5); ffl.setLvec(lvec);  // DEBUG 
        //Mat3d lvec=lvec0; lvec.a=crashed_lvecs_a[isys]; ffl.setLvec(lvec);  // DEBUG 
        //Mat3d lvec=lvec0; lvec.a=crashed_lvecs_a[16]; ffl.setLvec(lvec);  // DEBUG 
        //Mat3d lvec=lvec0; lvec.a=crashed_lvecs_a[3]; ffl.setLvec(lvec);  // DEBUG 
        //Mat3d lvec=lvec0; lvec.a.addRandomCube(0.5); ffl.setLvec(lvec);  // DEBUG
        //Mat3d lvec=lvec0; lvec.a.addRandomCube(0.4); ffl.setLvec(lvec);  // DEBUG
        //Mat3d lvec=lvec0; lvec.a.addRandomCube(0.25); ffl.setLvec(lvec);  // DEBUG
        //Mat3d lvec=lvec0; lvec.b.x += randf(-0.5,0.5); ffl.setLvec(lvec);  // DEBUG
        //Mat3d lvec=lvec0; lvec.b.x += isys*0.00163; ffl.setLvec(lvec);  // DEBUG
        //Mat3d lvec=lvec0; lvec.b.x += isys*0.15; ffl.setLvec(lvec);  // DEBUG
        //Mat3d lvec=lvec0; lvec.a.y += isys*2.25; ffl.setLvec(lvec);  // DEBUG
        //printf( "   Vec3d{%g,%g,%g}, \n",  lvec.a.x, lvec.a.y, lvec.a.z  );

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

void unpack_system(  int isys, MMFFsp3_loc& ff, bool bForces=false, bool bVel=false ){
    //rintf("MolWorld_sp3_multi::unpack_system(%i) \n", isys);
    int i0n = isys * ocl.nnode;
    int i0a = isys * ocl.nAtoms;
    int i0v = isys * ocl.nvecs;
    unpack( ff.nvecs, ff.apos, atoms+i0v );
    if(bForces){ unpack( ff.nvecs, ff.fapos, aforces+i0v ); }
    if(bVel   ){ unpack( ff.nvecs, ff.vapos, avel   +i0v ); }
}

void upload_sys( int isys, bool bParams=false, bool bForces=0, bool bVel=true, bool blvec=true ){
    //printf("MolWorld_sp3_multi::upload() \n");
    int i0v   = isys * ocl.nvecs;
    int i0a   = isys * ocl.nAtoms;
    int err=0;
    err|= ocl.upload( ocl.ibuff_atoms,  atoms,    ocl.nvecs,  i0v );
    err|= ocl.upload( ocl.ibuff_constr,  constr,  ocl.nAtoms, i0a );
    err|= ocl.upload( ocl.ibuff_constrK, constrK, ocl.nAtoms, i0a );
    err|= ocl.upload( ocl.ibuff_bboxes,  bboxes, 1, isys  );
    if(bForces){ err|= ocl.upload( ocl.ibuff_aforces, aforces, ocl.nvecs, i0v ); }
    if(bVel   ){ 
        err|= ocl.upload( ocl.ibuff_avel,    avel,    ocl.nvecs, i0v ); 
        err|= ocl.upload( ocl.ibuff_cvf,     cvfs,    ocl.nvecs, i0v ); 
    }
    if(blvec){
        int i0pbc = isys * ocl.npbc;
        err|= ocl.upload( ocl.ibuff_lvecs,     lvecs     , 1, isys  );
        err|= ocl.upload( ocl.ibuff_ilvecs,    ilvecs    , 1, isys  );
        err|= ocl.upload( ocl.ibuff_pbcshifts, pbcshifts , ocl.npbc, i0pbc );
    }
    if(bParams){
        int i0n   = isys * ocl.nnode;
        int i0bk  = isys * ocl.nbkng;
        err|= ocl.upload( ocl.ibuff_neighs,    neighs    , ocl.nAtoms, i0a  );
        err|= ocl.upload( ocl.ibuff_neighCell, neighCell , ocl.nAtoms, i0a  );
        err|= ocl.upload( ocl.ibuff_bkNeighs,  bkNeighs  , ocl.nbkng, i0bk  );
        err|= ocl.upload( ocl.ibuff_bkNeighs_new, bkNeighs_new, ocl.nbkng, i0bk  );
        //err|= ocl.upload( ocl.ibuff_neighForce, neighForce  );
        err|= ocl.upload( ocl.ibuff_REQs,   REQs,  ocl.nAtoms, i0a  );
        err|= ocl.upload( ocl.ibuff_MMpars, MMpars, ocl.nnode, i0n );
        err|= ocl.upload( ocl.ibuff_BLs, BLs, ocl.nnode, i0n );
        err|= ocl.upload( ocl.ibuff_BKs, BKs, ocl.nnode, i0n );
        err|= ocl.upload( ocl.ibuff_Ksp, Ksp, ocl.nnode, i0n );
        err|= ocl.upload( ocl.ibuff_Kpp, Kpp, ocl.nnode, i0n );
    }
    err |= ocl.finishRaw(); 
    OCL_checkError(err, "MolWorld_sp2_multi::upload().finish");
}

void download_sys( int isys, bool bForces=false, bool bVel=false ){
    //int i0n   = isys * ocl.nnode;
    int i0v   = isys * ocl.nvecs;
    //int i0a   = isys * ocl.nAtoms;
    //int i0bk  = isys * ocl.nbkng;
    //int i0pbc = isys * ocl.npbc;
    int err=0;
    //printf("MolWorld_sp3_multi::download() \n");
    ocl.download             ( ocl.ibuff_atoms,   atoms,    ocl.nvecs, i0v );
    if(bForces){ ocl.download( ocl.ibuff_aforces, aforces , ocl.nvecs, i0v ); }
    if(bVel   ){ ocl.download( ocl.ibuff_avel,    avel ,    ocl.nvecs, i0v ); }
    err |= ocl.finishRaw(); 
    OCL_checkError(err, "MolWorld_sp2_multi::upload().finish");
}

void upload(  bool bParams=false, bool bForces=0, bool bVel=true, bool blvec=true ){
    //printf("MolWorld_sp3_multi::upload() \n");
    int err=0;
    err|= ocl.upload( ocl.ibuff_atoms,  atoms  );
    err|= ocl.upload( ocl.ibuff_constr,  constr );
    err|= ocl.upload( ocl.ibuff_constrK, constrK );
    err|= ocl.upload( ocl.ibuff_bboxes,  bboxes );
    if(bForces){ err|= ocl.upload( ocl.ibuff_aforces, aforces ); }
    if(bVel   ){ err|= ocl.upload( ocl.ibuff_avel,    avel    ); }
    if(blvec){
        err|= ocl.upload( ocl.ibuff_lvecs,   lvecs );
        err|= ocl.upload( ocl.ibuff_ilvecs, ilvecs );
        err|= ocl.upload( ocl.ibuff_pbcshifts, pbcshifts );
    }
    if(bParams){
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

void evalVF( int n, Quat4f* fs, Quat4f* vs, FIRE& fire, Quat4f& MDpar ){
    //printf("evalVF() fire[%i]\n", fire.id );
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

void evalVF_new( int n, Quat4f* cvfs, FIRE& fire, Quat4f& MDpar, bool bExploring ){
    //printf("evalVF() fire[%i]\n", fire.id );
    Vec3f cvf=Vec3fZero;
    for(int i=0; i<n; i++){ 
        cvf.add( cvfs[i].f ); 
        cvfs[i] = Quat4fZero; 
    }
    fire.vv=cvf.x;
    fire.ff=cvf.y;
    fire.vf=cvf.z;
    fire.update_params();
    // if(bExploring){
    //     MDpar.x = fire.par->dt_max;
    //     MDpar.y = 1.0;
    //     MDpar.z = 1.0;
    //     MDpar.w = 0.0;
    // }else{
        MDpar.x = fire.dt;
        MDpar.y = 1 - fire.damping;
        MDpar.z = fire.cv;
        MDpar.w = fire.cf;
    //}
}

bool updateMultiExploring( double Fconv=1e-6, float fsc = 0.02, float tsc = 0.3 ){
    int err=0;
    double F2conv = Fconv*Fconv;
    bool bGroupUpdate=false;
    bool bExploring = false;
    for(int isys=0; isys<nSystems; isys++){
        int i0v = isys * ocl.nvecs;
        double f2 = fire[isys].ff;
        // -------- Global Optimization
        if( ( f2 < F2conv ) && (!gopts[isys].bExploring) ){            // Start Exploring
            gopts[isys].startExploring();
            bGroupUpdate=true;
            //printf("MolWorld_sp3_multi::evalVFs() isys=%3i Start Exploring \n", isys );
            //for(int ig=0; ig<ocl.nGroup; ig++){ setGroupDrive(isys, ig, {fsc,fsc,0.0} ); }   // Shift driver
            //for(int ig=0; ig<ocl.nGroup; ig++){ setGroupDrive(isys, ig, Vec3fZero, {tsc,0.0,0.0} ); }   // group rotate driver
            for(int ig=0; ig<ocl.nGroup; ig++){ setGroupDrive(isys, ig, {fsc,fsc,0.0}, {tsc,0.0,0.0} ); }   // Shift driver
        }
        if( gopts[isys].update() ){ // Stop Exploring
            bGroupUpdate=true;
            //printf("MolWorld_sp3_multi::evalVFs() isys=%3i Stop Exploring \n", isys );
            for(int ig=0; ig<ocl.nGroup; ig++){ setGroupDrive(isys, ig, Vec3fZero, Vec3fZero ); }
        };
        bExploring |= gopts[isys].bExploring;
    }
    //err |= ocl.upload( ocl.ibuff_TDrive, TDrive );
    //printf("MolWorld_sp3_multi::evalVFs() bGroupUpdate=%i \n", bGroupUpdate );
    if(bGroupUpdate){
        //printf("MolWorld_sp3_multi::evalVFs() bGroupUpdate=%i \n", bGroupUpdate );
        err |= ocl.upload( ocl.ibuff_gforces, gforces );
        err |= ocl.upload( ocl.ibuff_gtorqs , gtorqs  );
    }
    return bExploring;
}

double evalVFs( double Fconv=1e-6 ){
    double F2conv = Fconv*Fconv;
    //printf("MolWorld_sp3_multi::evalVFs()\n");
    int err=0;
    //ocl.download( ocl.ibuff_aforces , aforces );
    //ocl.download( ocl.ibuff_avel    , avel    );
    ocl.download( ocl.ibuff_cvf    , cvfs  );
    err |= ocl.finishRaw();  OCL_checkError(err, "evalVFs().1");
    //printf("MolWorld_sp3_multi::evalVFs(%i) \n", isys);
    double F2max = 0;
    iSysFMax=-1;
    bool bGroupUpdate=false;
    for(int isys=0; isys<nSystems; isys++){
        int i0v = isys * ocl.nvecs;
        //evalVF( ocl.nvecs, aforces+i0v, avel   +i0v, fire[isys], MDpars[isys] );
        evalVF_new( ocl.nvecs, cvfs+i0v, fire[isys], MDpars[isys], gopts[isys].bExploring );
        double f2 = fire[isys].ff;
        if(f2>F2max){ F2max=f2; iSysFMax=isys; }
        // // -------- Global Optimization
        // if( ( f2 < F2conv ) && (!gopts[isys].bExploring) ){            // Start Exploring
        //     //printf( "evalVFs() iSys=%i CONVERGED |F|=%g \n", isys, sqrt(f2) );
        //     gopts[isys].startExploring();
        //     bGroupUpdate=true;
        //     //printf("MolWorld_sp3_multi::evalVFs() isys=%3i Start Exploring \n", isys );
        //     //for(int ig=0; ig<ocl.nGroup; ig++){ setGroupDrive(isys, ig, {fsc,fsc,0.0} ); }   // Shift driver
        //     for(int ig=0; ig<ocl.nGroup; ig++){ setGroupDrive(isys, ig, Vec3fZero, {tsc,0.0,0.0} ); }   // group rotate driver
        // }
        // if( gopts[isys].update() ){ // Stop Exploring
        //     bGroupUpdate=true;
        //     //printf("MolWorld_sp3_multi::evalVFs() isys=%3i Stop Exploring \n", isys );
        //     for(int ig=0; ig<ocl.nGroup; ig++){ setGroupDrive(isys, ig, Vec3fZero, Vec3fZero ); }
        // };
        // if( gopts[isys].bExploring ){
        //     TDrive[isys].x = 1000;  // Temperature [K]
        //     TDrive[isys].y = 0.1;  // gamma_damp
        //     TDrive[isys].z = 0;    // ?
        //     TDrive[isys].w = randf(-1.0,1.0); 
        // }else{
        //     TDrive[isys].y = -1.0; // gamma_damp
        // }
        //printf( "evalVFs()[iSys=%i]  bExploring=%i (%i/%i)  |F|=%g \n", isys, gopts[isys].bExploring,   gopts[isys].istep, gopts[isys].nExplore,  sqrt(f2) );
        //printf( "evalF2[sys=%i] |f|=%g MDpars(dt=%g,damp=%g,cv=%g,cf=%g)\n", isys, sqrt(f2), MDpars[isys].x, MDpars[isys].y, MDpars[isys].z, MDpars[isys].w );
        //F2max = fmax( F2max, fire[isys].ff );
    }
    //printf( "MDpars{%g,%g,%g,%g}\n", MDpars[0].x,MDpars[0].y,MDpars[0].z,MDpars[0].w );
    err |= ocl.upload( ocl.ibuff_MDpars, MDpars );
    // err |= ocl.upload( ocl.ibuff_TDrive, TDrive );
    err |= ocl.upload( ocl.ibuff_cvf   , cvfs   );
    // //printf("MolWorld_sp3_multi::evalVFs() bGroupUpdate=%i \n", bGroupUpdate );
    // if(bGroupUpdate){
    //     //printf("MolWorld_sp3_multi::evalVFs() bGroupUpdate=%i \n", bGroupUpdate );
    //     err |= ocl.upload( ocl.ibuff_gforces, gforces );
    //     err |= ocl.upload( ocl.ibuff_gtorqs , gtorqs  );
    // }
    //err |= ocl.finishRaw();  
    OCL_checkError(err, "evalVFs().2");
    return F2max;
}

void setGroupDrive(int isys, int ig, Vec3f fsc=Vec3fZero, Vec3f tsc=Vec3fZero){
    int igs = ig + isys*ocl.nGroup;
    gforces[igs] = Quat4f{randf(-fsc.x,fsc.x),randf(-fsc.y,fsc.y),randf(-fsc.z,fsc.z),0.0};
    gtorqs [igs]   = Quat4f{randf(-tsc.x,tsc.x),randf(-tsc.y,tsc.y),randf(-tsc.z,tsc.z),0.0};
    //gforces[igs] = Quat4fZero;
    //gtorqs [igs] = Quat4fZero;
}

void setGroupDrives(){
    for(int isys=0; isys<nSystems; isys++){
        for(int ig=0; ig<ocl.nGroup; ig++){
           setGroupDrive(isys, ig,  {1.0,1.0,0.0} );
        }
    }
    int err=0;
    err |= ocl.upload( ocl.ibuff_gforces, gforces );
    err |= ocl.upload( ocl.ibuff_gtorqs , gtorqs  );
}

double evalF2(){
    int err=0;
    ocl.download( ocl.ibuff_aforces , aforces );
    err |= ocl.finishRaw();  OCL_checkError(err, "evalF2().1");
    double F2max = 0;
    iSysFMax=-1;
    for(int isys=0; isys<nSystems; isys++){
        double f2 = 0;
        Quat4f* fs = aforces + isys*ocl.nvecs;
        for(int i=0; i<ocl.nvecs; i++){ f2 += fs[i].f.norm2(); }
        //printf( "evalF2[sys=%i] |f|=%g \n", isys, sqrt(f2) );
        if(f2>F2max){ F2max=f2; iSysFMax=isys; }
    }
    OCL_checkError(err, "evalF2().2");
    return F2max;
}

// ===============================================
//       Implement    MultiSolverInterface
// ===============================================

void move_MultiConstrain( Vec3d d, Vec3d di, float Kfixmin=0.001f ){
    printf( "MolWorld_sp3_multi::move_MultiConstrain()\n" );
    for(int isys=0; isys<nSystems; isys++){
        int i0a   = isys * ocl.nAtoms;
        int i0v   = isys * ocl.nvecs;
        for(int ia=0; ia<ffls[isys].natoms; ia++){
            if( ffls[isys].constr[ia].w > Kfixmin ){
                //ffls[isys].constr[ia].w=1.0f;
                ffls[isys].constr[ia].f.add( d + di*isys ); 
                constr[ia+i0a]=(Quat4f)ffls[isys].constr[ia];
            };
        };
    }
    //for(int isys=0; isys<nSystems; isys++){ printf("replica %i : ", isys); ffls[isys].printAtomsConstrains(); }
    //upload( ocl.ibuff_constr, constr );
}

virtual void setConstrains(bool bClear=true, double Kfix=1.0 ){
    printf("MolWorld_sp3_multi::setConstrains()\n");
    MolWorld_sp3::setConstrains( bClear, Kfix );
    for(int isys=0; isys<nSystems; isys++){
        int i0a   = isys * ocl.nAtoms;
        int i0v   = isys * ocl.nvecs;
        // set it to GPU buffers
        for( int i=0; i<ocl.nAtoms; i++ ){ constr[i+i0a].w=-1;                                    }
        for( int i : constrain_list     ){ constr[i+i0a].w=Kfix; constr[i+i0a].f= atoms[i+i0v].f; }
        // set it to CPU ffls
        for( int i=0; i<ffls[isys].natoms; i++ ){ ffls[isys].constr[i].w=-1;                            }
        for( int i : constrain_list            ){ ffls[isys].constr[i].w=Kfix; ffls[isys].constr[i].f=ffls[isys].apos[i]; }
    }
    move_MultiConstrain( Vec3d{0.0, 0.4, 5.0}, Vec3d{0.0, 0.4, 0.0} );
    /*
    // temporary hack to show up paralell manipulation
    double Kfixmin = 0.001; 
    for(int isys=0; isys<nSystems; isys++){
        int i0a   = isys * ocl.nAtoms;
        int i0v   = isys * ocl.nvecs;
        // set it to GPU buffers
        //for( int i : constrain_list     ){ constr[i+i0a].w=Kfix; constr[i+i0a].f= atoms[i+i0v].f; }
        for(int ia=0; ia<ffls[isys].natoms; ia++){
            if( ffls[isys].constr[ia].w > Kfixmin ){
                //ffls[isys].constr[ia].w=1.0f;
                ffls[isys].constr[ia].f.add( {0.0, 0.4*isys, 5.0} ); 
                constr[ia+i0a]=(Quat4f)ffls[isys].constr[ia];
            };
        };
    }
    for(int isys=0; isys<nSystems; isys++){ printf("replica %i : ", isys); ffls[isys].printAtomsConstrains(); }
    */

    upload( ocl.ibuff_constr, constr );
}





virtual int paralel_size( )override{ return nSystems; }

virtual double solve_multi ( int nmax, double tol )override{
    //return eval_MMFFf4_ocl( nmax, tol );
    //return run_ocl_opt( nmax, tol );
    int nitrdione;
    switch(iParalel){
        case -1: bOcl=false;  nitrdione = run_multi_serial(nmax,tol);  break; 
        case  0:
        case  1: bOcl=false;  nitrdione = run_omp_ocl( nmax, tol    ); break; 
        case  2: bOcl=true;   nitrdione = run_omp_ocl( nmax, tol    ); break; 
        case  3: bOcl=true;   nitrdione = run_ocl_opt( nmax, tol    ); break; 
        //case  3: bOcl=true; nitrdione = run_ocl_loc( nmax, tol, 1 ); break; 
        //case  4: bOcl=true; nitrdione = run_ocl_loc( nmax, tol, 2 ); break; 
        //default:
        //    eval_NBFF_ocl_debug();
    }
    return nitrdione;
}

virtual void setGeom( int isys, Vec3d* ps, Mat3d *lvec, bool bPrepared )override{


    /* //   OLD
    
    int i0n = isys * ocl.nnode;
    int i0a = isys * ocl.nAtoms;
    int i0v = isys * ocl.nvecs;
    int i0pbc = isys * ocl.npbc;
    if(lvec){ 
        ffl.setLvec(*lvec); 
        evalPBCshifts( nPBC, ff.lvec, pbc_shifts ); // this must be before   ffl.initPi(pbc_shifts);

        printf( "MolWorld_sp3_multi::setGeom() lvec(a={}) ",  );
    }
    if(ps){
        for(int i=0; i<ffl.natoms; i++){ ffl.apos[i]=ps[i]; } 
        ffl.initPi(pbc_shifts);
    }
    pack_system(isys,ffl);
    */

    {
        if(lvec){
            ffls[isys].setLvec(*lvec);
            int i0pbc = isys*ocl.npbc;
            ffls[isys].setLvec ( *lvec );
            ffls[isys].makePBCshifts( nPBC, true );

            pack( ffls[isys].npbc,  ffls[isys].shifts, pbcshifts+i0pbc);
            Mat3_to_cl( ffls[isys].   lvec,  lvecs[isys] );
            Mat3_to_cl( ffls[isys].invLvec, ilvecs[isys] );  
        }
        
        if(ps){
            for(int i=0; i<ffls[isys].natoms; i++){ ffls[isys].apos[i]=ps[i]; } 
            //ffls[isys].initPi(pbc_shifts);

            //pack_system(isys,ffls[isys]);
        }

    }




    /*
    if(lvec){ ffl.setLvec(*lvec); };
    Mat3_to_cl( ffl.lvec,    lvecs [isys] );
    Mat3_to_cl( ffl.invLvec, ilvecs[isys] );
    //evalPBCshifts( nPBC, ffl.lvec, pbcshifts + i0pbc );
    evalPBCshifts( nPBC, ffl.lvec, pbc_shifts );
    pack( npbc, pbc_shifts, pbcshifts+i0pbc );
    printf( "setGeom[%i] lvec{%6.3f,%6.3f,%6.3f}{%6.3f,%6.3f,%6.3f}{%6.3f,%6.3f,%6.3f}\n",  isys,  ffl.lvec.a.x,ffl.lvec.a.y,ffl.lvec.a.z,   ffl.lvec.b.x,ffl.lvec.b.y,ffl.lvec.b.z,   ffl.lvec.c.x,ffl.lvec.c.y,ffl.lvec.c.z  ); 
    // we cannot pack directly ps because  pipos is not defined in setGeom
    for(int i=0; i<ffl.natoms; i++){ ffl.apos[i]=ps[i]; } 
    ffl.initPi();         
    pack( ocl.nvecs, ffl.apos, atoms+i0v );
    */

}

virtual double getGeom     ( int isys, Vec3d* ps, Mat3d *lvec, bool bPrepared )override{
    /*
    int i0n = isys * ocl.nnode;
    int i0a = isys * ocl.nAtoms;
    int i0v = isys * ocl.nvecs;
    //if( ! bPrepared ) ocl.download( ocl.ibuff_atoms,   atoms,   ocl.nvecs, i0v  );
    pack( ocl.nvecs, ps, atoms+i0v );
    //if(lvec){ lvec.a lvecs };
    return 0;
    */
   for(int ia=0;ia<ffls[isys].natoms;ia++){ ps[ia]=ffls[isys].apos[ia]; };
   //printf( "MolWorld_sp3_multi::getGeom(isys=%i) Etot %g \n", isys, ffls[isys].Etot );
   return ffls[isys].Etot;
}

virtual void downloadPop()override{
    ocl.download( ocl.ibuff_atoms,  atoms  );
}

virtual void uploadPop  ()override{
    int ntot = nSystems*ocl.nvecs;
    set ( ntot, avel);
    ocl.upload( ocl.ibuff_avel ,   avel   );
    ocl.upload( ocl.ibuff_atoms,   atoms  );
    ocl.upload( ocl.ibuff_lvecs,   lvecs  );
    ocl.upload( ocl.ibuff_ilvecs,  ilvecs );
    //ocl.upload( ocl.ibuff_constr, constr );
}

virtual void change_lvec( const Mat3d& lvec )override{
    //printf("MolWold_sp3_multi::change_lvec() iSystemCur=%i \n", iSystemCur);
    int isys  = iSystemCur;
    int i0pbc = isys*ocl.npbc; 
    ffls[isys].setLvec ( lvec );
    ffls[isys].makePBCshifts( nPBC, true );
    pack( ffls[isys].npbc,  ffls[isys].shifts, pbcshifts+i0pbc );
    Mat3_to_cl( ffls[isys].   lvec,  lvecs[isys] );
    Mat3_to_cl( ffls[isys].invLvec, ilvecs[isys] );
    
    builder.lvec = ffls[iSystemCur].lvec;
    ffl.setLvec  ( ffls[isys].lvec );
    evalPBCshifts( nPBC, ffl.lvec, pbc_shifts );

    int err=0;
    err|= ocl.upload( ocl.ibuff_lvecs,     lvecs,  1,  isys  );
    err|= ocl.upload( ocl.ibuff_ilvecs,    ilvecs, 1,  isys  );
    err|= ocl.upload( ocl.ibuff_pbcshifts, pbcshifts , i0pbc );
    err |= ocl.finishRaw(); 
    OCL_checkError(err, "MolWorld_sp3_multi::upload().finish");
}

virtual void add_to_lvec( const Mat3d& dlvec )override{
    //printf("MolWold_sp3_multi::add_to_lvec()\n");
    for(int isys=0; isys<nSystems; isys++){
        int i0pbc = isys*ocl.npbc;
        ffls[isys].setLvec ( ffls[isys].lvec+dlvec );
        ffls[isys].makePBCshifts( nPBC, true );
        pack( ffls[isys].npbc,  ffls[isys].shifts, pbcshifts+i0pbc);
        Mat3_to_cl( ffls[isys].   lvec,  lvecs[isys] );
        Mat3_to_cl( ffls[isys].invLvec, ilvecs[isys] );
    }

    builder.lvec = ffls[iSystemCur].lvec;
    ffl.setLvec ( ffls[iSystemCur].lvec );
    evalPBCshifts( nPBC, ffl.lvec, pbc_shifts );
    
    int err=0;
    err|= ocl.upload( ocl.ibuff_lvecs,     lvecs     );
    err|= ocl.upload( ocl.ibuff_ilvecs,    ilvecs    );
    err|= ocl.upload( ocl.ibuff_pbcshifts, pbcshifts );
    err |= ocl.finishRaw(); 
    OCL_checkError(err, "MolWorld_sp2_multi::upload().finish");
}

virtual void upload_pop( const char* fname ){
    printf("MolWorld_sp3::upload_pop(%s)\n", fname );
    gopt.loadPopXYZ( fname );
    int nmult=gopt.population.size(); 
    int npara=paralel_size(); if( nmult!=npara ){ printf("ERROR in GlobalOptimizer::lattice_scan_1d_multi(): (imax-imin)=(%i) != solver.paralel_size(%i) => Exit() \n", nmult, npara ); exit(0); }
    gopt.upload_multi(nmult,0, true, true );

    // //int initMode=1;
    // int initMode = 0;
    // int nstesp   = 40;
    // //int nstesp = 2;
    // //gopt.tolerance = 0.02;
    // gopt.tolerance = 0.01;
    // Mat3d dlvec =  Mat3d{   0.2,0.0,0.0,    0.0,0.0,0.0,    0.0,0.0,0.0  };
    // gopt.lattice_scan_2d_multi( nstesp, dlvec, initMode, "lattice_scan_2d_multi.xyz" );
}

virtual void optimizeLattice_1d( int n1, int n2, Mat3d dlvec ){
    printf("\n\n\n######### MolWorld_sp3_multi::optimizeLattice_1d(%i,%i) \n", n1, n2 );
    int initMode=1;
    // //int initMode = 0;
    // //int nstesp   = 40;
    // //int nstesp = 2;
    // //gopt.tolerance = 0.02;
    // gopt.tolerance = 0.01;
    ///Mat3d dlvec =  Mat3d{   0.2,0.0,0.0,    0.0,0.0,0.0,    0.0,0.0,0.0  };
    gopt.lattice_scan_2d_multi( n1, dlvec, initMode, "lattice_scan_2d_multi.xyz" );
    
}

int saveSysXYZ( int isys, const char* fname, const char* comment="#comment", bool bNodeOnly=false, const char* mode="w", Vec3i nPBC=Vec3i{1,1,1} ){ 
    char str_tmp[1024];
    MMFFsp3_loc& ffl = ffls[isys];
    if(bPBC){ sprintf( str_tmp, "lvs %8.3f %8.3f %8.3f    %8.3f %8.3f %8.3f    %8.3f %8.3f %8.3f %s", ffl.lvec.a.x, ffl.lvec.a.y, ffl.lvec.a.z, ffl.lvec.b.x, ffl.lvec.b.y, ffl.lvec.b.z, ffl.lvec.c.x, ffl.lvec.c.y, ffl.lvec.c.z, comment ); }
    else    { sprintf( str_tmp, "%s", comment ); }
    return params.saveXYZ( fname, (bNodeOnly ? ffl.nnode : ffl.natoms) , ffl.atypes, ffl.apos, str_tmp, ffl.REQs, mode, true, nPBC, ffl.lvec ); 
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
    evalPBCshifts( nPBC, ffl.lvec, pbc_shifts );
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

virtual void evalAFMscan( GridShape& scan, Quat4f*& OutFE, Quat4f*& OutPos, Quat4f** ps=0, bool bSaveDebug=false )override{ 
    int err=0;
    //Quat4f* OutFE =0;
    //Quat4f* OutPos=0;
    //if( OutFE_  ){ OutFE =OutFE_;  }else{ OutFE =new Quat4f[scan.n.totprod()];   }
    //if( OutPos_ ){ OutPos=OutPos_; }else{ OutPos=new Quat4f[scan.n.totprod()]; }
    if( OutFE ==0 ){ OutFE =new Quat4f[scan.n.totprod()]; }
    if( OutPos==0 ){ OutPos=new Quat4f[scan.n.totprod()]; }
    int nxy = scan.n.x*scan.n.y;
    _realloc( afm_ps, nxy );

    // ToDo: ztop should be probably found and used in  evalAFM_FF()
    int i0v = iSystemCur*ocl.nvecs;
    float ztop= -1e+300;
    for( int i=0; i<ocl.nAtoms; i++ ){ ztop=fmax( atoms[i0v+i].z, ztop ); };
    printf( "MolWrold_sp3_multi::evalAFMscan() ztop=%g \n", ztop );

    scan.pos0=Vec3d{0.0,0.0,ztop+4.0+5.0+1.0};
    for(int iy=0; iy<scan.n.y; iy++) for(int ix=0; ix<scan.n.x; ix++){ afm_ps[iy*scan.n.x+ix].f = (Vec3f)(  scan.pos0 + scan.dCell.a*ix + scan.dCell.b*iy ); afm_ps[iy*scan.n.x+ix].w=0; }
    printf( "MolWrold_sp3_multi::evalAFMscan() scan.ns(%i,%i,%i) nxy=%i \n", scan.n.x,scan.n.y,scan.n.z,nxy );
    long T0=getCPUticks();
    ocl.PPAFM_scan( nxy, scan.n.z, afm_ps, OutFE, OutPos, scan.dCell.c.z );  
    printf( ">>time(surf2ocl.download() %g \n", (getCPUticks()-T0)*tick2second );

    // ToDo: we should probably return vmin,vmax to MolGUI for visualization? 
    int ntot = scan.n.x*scan.n.y*scan.n.z;
    Quat4f vmin=Quat4fmax,vmax=Quat4fmin;
    for(int i=0; i<ntot; i++){ OutFE[i].update_bounds( vmin, vmax ); }
    printf( "MolWrold_sp3_multi::evalAFMscan() vmin{%g,%g,%g,%g} vmax{%g,%g,%g,%g} \n", vmin.x,vmin.y,vmin.z,vmin.w,  vmax.x,vmax.y,vmax.z,vmax.w  );
    //bSaveDebug=true;
    //bSaveDebug=false; 
    if(bSaveDebug){ 
        scan.saveXSF( "AFM_E.xsf",       (float*)OutFE, 4,3 );
        scan.saveXSF( "AFM_Fx.xsf",      (float*)OutFE, 4,0 );
        scan.saveXSF( "AFM_Fy.xsf",      (float*)OutFE, 4,1 );
        scan.saveXSF( "AFM_Fz.xsf",      (float*)OutFE, 4,2 );
        scan.saveXSF( "AFM_PPpos_w.xsf", (float*)OutPos, 4,3 );
        scan.saveXSF( "AFM_PPpos_x.xsf", (float*)OutPos, 4,0 );
        scan.saveXSF( "AFM_PPpos_y.xsf", (float*)OutPos, 4,1 );
        scan.saveXSF( "AFM_PPpos_z.xsf", (float*)OutPos, 4,2 );
    }
    if(ps){ *ps=afm_ps; }else{ delete [] ps;   }
    //if(OutFE==0 ){ delete [] OutFE; }
    //if(OutPos==0){ delete [] OutPos; }
}

virtual void evalAFM_FF( GridShape& grid, Quat4f* data_=0, bool bSaveDebug=false )override{ 
    int err=0;
    Quat4f* data=0;
    bool bDownload = bSaveDebug && (data_==0);
    if( bDownload ){ data=new Quat4f[grid.n.totprod()]; }
    long T0=getCPUticks();
    Vec3i afm_nPBC{1,1,0};
    autoNPBC( grid.cell, afm_nPBC, 30.0 );
    ocl.PPAFM_makeFF( iSystemCur, grid, afm_nPBC );
    err |=  ocl.finishRaw(); OCL_checkError(err, "evalAFM_FF.PPAFM_makeFF"); printf( ">>time(ocl.makeGridFF() %g \n", (getCPUticks()-T0)*tick2second );
    if( bDownload )ocl.download( ocl.itex_FEAFM, data );
    err |=  ocl.finishRaw();  OCL_checkError(err, "surf2ocl.download.finish");   
    printf( ">>time(surf2ocl.finishRaw() %g[s] download=%i \n", (getCPUticks()-T0)*tick2second, bDownload );
    if(bSaveDebug){
        grid.saveXSF( "AFM_FF_E.xsf", (float*)data, 4,3 );
        grid.saveXSF( "AFM_FF_z.xsf", (float*)data, 4,2 );
    }
    if( bDownload && (data_==0) ){ delete [] data; }
    // downalod ?
}

// ==================================
//          SETUP OCL KERNELS 
// ==================================

void setup_MMFFf4_ocl(){
    printf("MolWorld_sp3_multi::setup_MMFFf4_ocl() \n");
    ocl.nDOFs.x=ff.natoms;
    ocl.nDOFs.y=ff.nnode;
    if(!task_cleanF)   task_cleanF = ocl.setup_cleanForceMMFFf4 ( ffl.natoms, ffl.nnode       );
    if(!task_move  )   task_move   = ocl.setup_updateAtomsMMFFf4( ffl.natoms, ffl.nnode       );
    if(!task_print )   task_print  = ocl.setup_printOnGPU       ( ffl.natoms, ffl.nnode       );
    if(!task_MMFF  )   task_MMFF   = ocl.setup_getMMFFf4        ( ffl.natoms, ffl.nnode, bPBC );

    Vec3i nPBC_=nPBC; if(!bPBC){ nPBC_=Vec3iZero; }; printf( "MolWorld_sp3_multi::setup_MMFFf4_ocl() bPBC=%i nPBC(%i,%i,%i) n", bPBC, nPBC_.x,nPBC_.y,nPBC_.z );
    if((!task_NBFF_Grid)&&bGridFF ){ task_NBFF_Grid = ocl.setup_getNonBond_GridFF( ffl.natoms, ffl.nnode, nPBC_ ); } 
    if(!task_NBFF                 ){ 
        task_NBFF      = ocl.setup_getNonBond       ( ffl.natoms, ffl.nnode, nPBC_ ); 
    }

    // OCLtask* getSurfMorse(  Vec3i nPBC_, int na=0, float4* atoms=0, float4* REQs=0, int na_s=0, float4* atoms_s=0, float4* REQs_s=0,  bool bRun=true, OCLtask* task=0 ){
    if(!task_SurfAtoms && bSurfAtoms ){ 
        if( bDONE_surf2ocl==false ){ printf("ERROR in MolWorld_sp3_multi::setup_MMFFf4_ocl() must call MolWorld_sp3_multi::surf2ocl() first\n"); }
        task_SurfAtoms = ocl.getSurfMorse( gridFF.nPBC, gridFF.natoms, 0,0,0,0,0,false ); 
    } 
    // ocl.makeGridFF( gridFF.grid, nPBC, gridFF.natoms, (float4*)atoms_surf, (float4*)REQs_surf, true );
    /// HERE_HERE

    task_GroupUpdate=ocl.setup_updateGroups( );
    task_GroupForce =ocl.setup_groupForce(   );

    //exit(0);

    // if(!task_NBFF  ) { 
    //     if( bGridFF ){ task_NBFF  = ocl.setup_getNonBond_GridFF( ff4.natoms, ff4.nnode, nPBC ); } 
    //     else         { task_NBFF  = ocl.setup_getNonBond       ( ff4.natoms, ff4.nnode, nPBC ); }
    // }
    if(!task_cleanF)task_cleanF = ocl.setup_cleanForceMMFFf4 ( ffl.natoms, ffl.nnode       );

    //exit(0);
}

void setup_NBFF_ocl(){
    printf("MolWorld_sp3_multi::setup_NBFF_ocl() \n");
    ocl.nDOFs.x=ff.natoms;
    ocl.nDOFs.y=ff.nnode;
    if(!task_cleanF )task_cleanF = ocl.setup_cleanForceMMFFf4 ( ffl.natoms, ffl.nnode        );
    if(!task_move   )task_move   = ocl.setup_updateAtomsMMFFf4( ffl.natoms, ffl.nnode        ); 
    if(!task_print  )task_print  = ocl.setup_printOnGPU       ( ffl.natoms, ffl.nnode        );
    if(!task_NBFF   )task_NBFF   = ocl.setup_getNonBond       ( ffl.natoms, ffl.nnode, nPBC  );
}

void picked2GPU( int ipick,  double K ){
    printf( "picked2GPU() ipick %i iSystemCur %i \n", ipick, iSystemCur );
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
        for(int i=0; i<ocl.nAtoms; i++){ constr[i0a + i].w=-1.0;  };
        for(int i: constrain_list     ){ constr[i0a + i].w=Kfix;  };
        // for(int ia=0; ia<ocl.nAtoms; ia++){ 
        //     int iaa = i0a + ia; 
        //     constr [iaa].w=ffls[iSystemCur]->constr [ia];
        //     constrK[iaa].w=ffls[iSystemCur]->constrK[ia];
        // };
    }
    //for(int i=0; i<ocl.nAtoms; i++){ printf( "CPU:constr[%i](%7.3f,%7.3f,%7.3f |K= %7.3f) \n", i, constr[i0a+i].x,constr[i0a+i].y,constr[i0a+i].z,  constr[i0a+i].w   ); }
    ocl.upload( ocl.ibuff_constr,  constr  );   // ToDo: instead of updating the whole buffer we may update just relevant part?
    ocl.upload( ocl.ibuff_constrK, constrK );
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
        double t = ( getCPUticks()-T0 )*tick2second;
        err |= ocl.finishRaw(); OCL_checkError(err, "eval_MMFFf4_ocl.test.finish");
        printf("eval_MMFFf4_ocl.test(itest=%i) time=%7.3f[ms] %7.3f[us/iter] niter=%i na=%i \n", itest, t*1000, t*1e+6, niter, ocl.nAtoms );
        return niter;
        //exit(0);    
    }else{ for(int i=0; i<nVFs; i++){
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

            /*
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
            */

        }
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
    download(false,false);
    err |= ocl.finishRaw(); OCL_checkError(err, "eval_MMFFf4_ocl.finish");
    printf("eval_MMFFf4_ocl() time=%7.3f[ms] niter=%i \n", ( getCPUticks()-T0 )*tick2second*1000 , niterdone );
    }

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

int run_ocl_loc( int niter, double Fconv=1e-6 , int iVersion=1 ){ 
    //printf("MolWorld_sp3_multi::run_ocl_loc() niter=%i \n", niter );
    double F2conv = Fconv*Fconv;
    picked2GPU( ipicked,  1.0 );
    int err=0;
    if( task_MMFF==0)setup_MMFFf4_ocl();
    int nPerVFs = _min(10,niter);
    int nVFs    = niter/nPerVFs;
    int niterdone=0;
    double F2=0;
    if( task_MMFFloc1==0 ){ task_MMFFloc1=ocl.setup_evalMMFFf4_local1( niter ); }
    if( task_MMFFloc2==0 ){ task_MMFFloc2=ocl.setup_evalMMFFf4_local2( niter ); }
    if     ( iVersion==1 ){ task_MMFFloc=task_MMFFloc1; }
    else if( iVersion==2 ){ task_MMFFloc=task_MMFFloc2; }
    long T0 = getCPUticks();
    for(int i=0; i<nVFs; i++){
        for(int j=0; j<nPerVFs; j++){
            task_MMFFloc->enque_raw(); 
            niterdone++;
        }
        ocl.download( ocl.ibuff_atoms,  atoms  );
        F2 = evalVFs();
        if( F2<F2conv  ){ 
            double t=(getCPUticks()-T0)*tick2second;
            printf( "run_ocl_loc(nsys=%i|iPara=%i) CONVERGED in <%i steps, |F|(%g)<%g time %g[ms] %g[us/step] bGridFF=%i \n", nSystems,iParalel, niterdone, sqrt(F2), Fconv, t*1000, t*1e+6/niterdone, bGridFF ); 
            return niterdone; 
        }
    }
    double t=(getCPUticks()-T0)*tick2second;
    printf( "run_ocl_loc(nsys=%i|iPara=%i) NOT CONVERGED in %i steps, |F|(%g)>%g time %g[ms] %g[us/step] bGridFF=%i iSysFMax=%i \n", nSystems,iParalel, niter, sqrt(F2), Fconv, t*1000, t*1e+6/niterdone, bGridFF, iSysFMax ); 
    //err |= ocl.finishRaw(); 
    //printf("eval_MMFFf4_ocl() time=%7.3f[ms] niter=%i \n", ( getCPUticks()-T0 )*tick2second*1000 , niterdone );
    return niterdone;
}

int debug_eval(){
    printf("MolWorld_sp3_multi::debug_eval() GridFF.npbc(%i)  GridFF.nPBC(%i,%i,%i) \n", gridFF.nPBC.x,gridFF.nPBC.y,gridFF.nPBC.z, gridFF.npbc);
    int err  = 0;
    int isys = 0;
    int ia   = 0;
    ffls[isys].cleanForce();
    //gridFF.evalMorsePBC_sym( ffls[isys].apos[ia], ffls[isys].REQs[ia], ffls[isys].fapos[ia] ); 
    gridFF.evalMorsePBCatoms_sym( ffls[isys].natoms, ffls[isys].apos, ffls[isys].REQs, ffls[isys].fapos );
    for(int i=0; i<ffls[isys].natoms; i++){
        Vec3d f = ffls[isys].fapos[i];
        printf( "cpu_aforces[%i] f(%g,%g,%g)\n", i, f.x, f.y, f.z );
    }
    err = task_cleanF->enque_raw();                    OCL_checkError(err, "MolWorld_sp3_multi::debug_eval().task_cleanF()"     ); 
    err = task_SurfAtoms->enque_raw();                 OCL_checkError(err, "MolWorld_sp3_multi::debug_eval().task_SurfAtoms()"  );  
    err = ocl.download( ocl.ibuff_atoms,    atoms   ); OCL_checkError(err, "MolWorld_sp3_multi::debug_eval().download(atoms)"   ); 
    err = ocl.download( ocl.ibuff_aforces,  aforces ); OCL_checkError(err, "MolWorld_sp3_multi::debug_eval().download(aforces)" ); 
    for(int i=0; i<ffls[isys].natoms; i++){
        //Vec3d f = ffls[isys].fapos;
        printf( "gpu_aforces[%i]  f(%g,%g,%g)\n", i, aforces[i].x, aforces[i].y, aforces[i].z );
    }
    exit(0);
}

int run_ocl_opt( int niter, double Fconv=1e-6 ){ 
    //printf("MolWorld_sp3_multi::run_ocl_opt() niter=%i \n", niter );
    //for(int i=0;i<npbc;i++){ printf( "CPU ipbc %i shift(%7.3g,%7.3g,%7.3g)\n", i, pbc_shifts[i].x,pbc_shifts[i].y,pbc_shifts[i].z ); }
    //debug_eval(); return 0;

    double F2conv = Fconv*Fconv;
    picked2GPU( ipicked,  1.0 );
    int err=0;
    if( task_MMFF==0)setup_MMFFf4_ocl();

    //int nPerVFs = _min(1,niter);
    int nPerVFs = _min(10,niter);
    //int nPerVFs = _min(50,niter);
    int nVFs    = niter/nPerVFs;
    long T0     = getCPUticks();
    int niterdone=0;
    double F2=0;

    bool dovdW=true;
    //bool dovdW=false;
    ocl.bSubtractVdW=dovdW;
    
    bool bExplore = false;
    //for(int isys=0; isys<nSystems; isys++){ if(gopts[isys].bExploring) bExplore = true; }

    for(int i=0; i<nVFs; i++){

        // if(bGroups){
        //     ocl.upload( ocl.ibuff_gtorqs,  gtorqs  );
        //     ocl.upload( ocl.ibuff_gforces, gforces );
        // }

        // bExplore = false;
        // for(int isys=0; isys<nSystems; isys++){ 
        //     if(gopts[isys].bExploring) bExplore = true; 
            
        // }

        bExplore = updateMultiExploring( Fconv );

        bool bGroupDrive = bGroups && bExplore;
        //bGroupDrive = false;

        //printf( "CPU::bbox(%g,%g,%g)(%g,%g,%g)(%g,%g,%g)\n", bbox.a.x,bbox.a.y,bbox.a.z,   bbox.b.x,bbox.b.y,bbox.b.z,   bbox.c.x,bbox.c.y,bbox.c.z );
        for(int ia=0; ia<ffl.natoms; ia++){
            if( ffl.constr[ia].w > 0 ) printf( "CPU:atom[%i] constr(%g,%g,%g|%g) constrK(%g,%g,%g|%g)\n", ia, ffl.constr[ia].x,ffl.constr[ia].y,ffl.constr[ia].z,ffl.constr[ia].w,   ffl.constrK[ia].x,ffl.constrK[ia].y,ffl.constrK[ia].z,ffl.constrK[ia].w  );
        }


        //bGroupDrive = true;
        // if(bGroupDrive){
        //     for(int ig=0; ig<ocl.nGroupTot; ig++){
        //         gtorqs[ig] = Quat4f{ 0.5f*sin(nloop*0.02f), 0.0f,0.0f,0.0f  };
        //         //printf( "CPU:gtorqs[%i](%g,%g,%g,%g)\n", ig, gtorqs[ig].x,gtorqs[ig].y,gtorqs[ig].z,gtorqs[ig].w );  
        //     }
        //     err |= ocl.upload( ocl.ibuff_gtorqs, gtorqs );  //OCL_checkError(err, "task_MMFF->enque_raw()");
        //     //for(int ig=0; ig<ocl.nGroupTot; ig++){ printf( "CPU:gtorqs[%i](%g,%g,%g,%g)\n", ig, gtorqs[ig].x,gtorqs[ig].y,gtorqs[ig].z,gtorqs[ig].w );     }
        // }

        for(int j=0; j<nPerVFs; j++){
            //if(bGroupDrive)printf( "bGroupDrive==true\n" );
            {    
                if( bGroupDrive )err |= task_GroupUpdate->enque_raw();
                if(dovdW)[[likely]]{
                    if(bSurfAtoms)[[likely]]{
                        if  (bGridFF)[[likely]]{ 
                            err |= task_NBFF_Grid ->enque_raw();   //OCL_checkError(err, "task_NBFF_Grid->enque_raw(); ");
                        }else { 
                            //printf( "task_NBFF(), task_SurfAtoms() \n" );
                            err |= task_NBFF     ->enque_raw();  //OCL_checkError(err, "MolWorld_sp3_multi::run_ocl_opt().task_NBFF()" ); 
                            err |= task_SurfAtoms->enque_raw();  //OCL_checkError(err, "MolWorld_sp3_multi::run_ocl_opt().task_SurfAtoms()" );
                        }
                    }else{ 
                        err |= task_NBFF      ->enque_raw();     //OCL_checkError(err, "task_NBFF->enque_raw();");
                    }
                }
                err |= task_MMFF->enque_raw();    //OCL_checkError(err, "task_MMFF->enque_raw()");

                if( bGroupDrive ) err |= task_GroupForce->enque_raw();
                err |= task_move->enque_raw();    //OCL_checkError(err, "task_move->enque_raw()");
            }
            niterdone++;
            nloop++;
        }
        //ocl.download( ocl.ibuff_atoms,    atoms   );
        //ocl.download( ocl.ibuff_aforces,  aforces );
        F2 = evalVFs( Fconv );

        // if(bGroups){ // Check Groups
        //     printf("MolWorld_sp3_multi::eval_MMFFf4_ocl()[niterdone=%i/%i]\n", niterdone, niter );
        //     ocl.download( ocl.ibuff_gcenters, gcenters );
        //     for(int i=0; i<ocl.nGroupTot; i++){
        //         int isys = i/ocl.nGroup; int ig = i - isys*ocl.nGroup;
        //         printf( "gcenter[%i|isys=%i,ig=%i] pos(%10.6f,%10.6f,%10.6f)\n", i, isys,ig, gcenters[i].x,gcenters[i].y,gcenters[i].z  ); 
        //     }
        // }

        /*
        { // ======= DEBUG - Check vs CPU
            int iS=iSystemCur;
            //printf( "MDpars[%i](dt=%g,damp=%g,cv=%g,cf=%g)cos_vf=%g\n", iSystemCur, MDpars[iS].x,MDpars[iS].y,MDpars[iS].z,MDpars[iS].w, fire[iS].cos_vf );
            int i0v = iSystemCur * ocl.nvecs; 
            ffl.cleanForce();
            ffl.eval(false);
            ffl.evalLJQs_ng4_simd    ();
            //ffl.evalLJQs_ng4_PBC_simd();
            //ffl.evalLJQs_ng4( ffl.neighs, gridFF.Rdamp );
            // ---- Compare to ffl
            bool ret=false;
            //printf("### run_ocl_opt: Compare GPU vs CPU ffl \n"); 
            ret |= compareVecs( ocl.nAtoms, ffl.fapos, aforces+i0v, 1e-4, 2 );
            if(ret){ printf("ERROR in run_ocl_opt()[nloop=%i] GPU.aforce[iSystemCur=%i] != ffl.evalLJQs_ng4_PBC() => Exit()\n", nloop, iSystemCur ); exit(0); }
            //else { printf("run_ocl_opt[nloop=%i] CHECKED: GPU.eval() == ffl.evalLJQs_ng4_PBC()\n", nloop); }
            unpack( ffl.nvecs, ffl.apos, atoms+i0v );
        }
        */
        //F2 = evalF2();
        if( F2<F2conv  ){ 
            double t=(getCPUticks()-T0)*tick2second;
            //printf( "run_omp_ocl(nSys=%i|iPara=%i) CONVERGED in %i/%i nsteps |F|=%g time=%g[ms]\n", nSystems, iParalel, itr,niter_max, sqrt(F2max), T1*1000 );
            if(verbosity>0)
            printf( "run_ocl_opt(nSys=%i|iPara=%i,bSurfAtoms=%i,bGridFF=%i,bExplore=%i) CONVERGED in %i/%i steps, |F|(%g)<%g time %g [ms]( %g [us/step]) bGridFF=%i \n", nSystems, iParalel, bSurfAtoms, bGridFF, bExplore, niterdone,niter, sqrt(F2), Fconv, t*1000, t*1e+6/niterdone, bGridFF ); 
            return niterdone; 
        }
    }

    if(bGroups){
        err|= ocl.download( ocl.ibuff_gcenters, gcenters );
        err|= ocl.download( ocl.ibuff_gfws, gfws );
        err|= ocl.download( ocl.ibuff_gups, gups );
    }
    err|= ocl.download( ocl.ibuff_atoms,    atoms   );
    err|= ocl.download( ocl.ibuff_aforces,  aforces );
    err|= ocl.finishRaw();  
    OCL_checkError(err, "run_ocl_opt().finishRaw()");

    double t=(getCPUticks()-T0)*tick2second;
    printf( "run_ocl_opt(nSys=%i|iPara=%i,bSurfAtoms=%i,bGridFF=%i,bExplore=%i) NOT CONVERGED in %i steps, |F|(%g)>%g time %g [ms]( %g [us/step]) bGridFF=%i iSysFMax=%i dovdW=%i \n", nSystems, iParalel, bSurfAtoms, bGridFF, bExplore,  niter, sqrt(F2), Fconv, t*1000, t*1e+6/niterdone, bGridFF, iSysFMax, dovdW ); 
    //err |= ocl.finishRaw(); 
    //printf("eval_MMFFf4_ocl() time=%7.3f[ms] niter=%i \n", ( getCPUticks()-T0 )*tick2second*1000 , niterdone );
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
        unpack( ff4.nvecs, ffl. apos, ff4.  apos );
        unpack( ff4.nvecs, ffl.fapos, ff4. fapos );
        fcog  = sum ( ffl.natoms, ffl.fapos   );
        tqcog = torq( ffl.natoms, ffl.apos, ffl.fapos );
        if(  fcog.norm2()>1e-8 ){ printf("WARRNING: eval_NBFF_ocl |fcog| =%g; fcog=(%g,%g,%g)\n", fcog.norm(),  fcog.x, fcog.y, fcog.z ); exit(0); }
        //if( tqcog.norm2()>1e-8 ){ printf("WARRNING: eval_MMFFf4_ocl |torq| =%g; torq=(%g,%g,%g)\n", tqcog.norm(),tqcog.x,tqcog.y,tqcog.z ); exit(0); }   // NOTE: torq is non-zero because pi-orbs have inertia
    }
    return 0;
}
int run_omp( int niter_max, double Fconv=1e-3, double Flim=1000, double timeLimit=0.02 ){
    //printf( "run_omp() niter_max=%i %dt=g Fconv%g \n", niter_max, dt, Fconv  ); 
    double F2conv=Fconv*Fconv;
    double F2max=0;
    int itr=0,niter=niter_max;
    long T0,T00;
    double T1,T2;
    int err=0;
    T00 = getCPUticks();
    //#pragma omp parallel shared(E,F2,ff,vv,vf,ffl) private(itr)
    #pragma omp parallel shared(niter,itr,F2max,T0,T1,T2,err)
    while(itr<niter){
        if(itr<niter){
        #pragma omp for
        for( int isys=0; isys<nSystems; isys++ ){
            //printf( "run_omp[itr=%i] isys=%i @cpu(%i/%i) \n", itr, isys, omp_get_thread_num(), omp_get_num_threads() );
            ffls[isys].cleanForce();
            ffls[isys].eval(false);
            { // Non-Bonded
                if(bPBC){ ffls[isys].Etot += ffls[isys].evalLJQs_ng4_PBC_simd(); }
                else    { ffls[isys].Etot += ffls[isys].evalLJQs_ng4_simd    (); } 
                if   (bGridFF){ gridFF.addForces( ffls[isys].natoms, ffls[isys].apos, ffls[isys].PLQs, ffls[isys].fapos ); }        // GridFF
            }
            if( (ipicked>=0) && (isys==iSystemCur) ){ 
                pullAtom( ipicked, ffls[isys].apos, ffls[isys].fapos );  
            };
            if(bConstrains){
                ffls[isys].Etot += constrs.apply( ffls[isys].apos, ffls[isys].fapos, &ffls[isys].lvec );
            }
        }        
        #pragma omp barrier
        #pragma omp single
        {   F2max=0;}
        #pragma omp for reduction(max:F2max)
        for( int isys=0; isys<nSystems; isys++ ){
            int i0v = isys*ocl.nvecs;
            double F2 = opts[isys].move_FIRE();
            F2max = fmax(F2max,F2);
        }
        #pragma omp single
        { 
            nloop++;
            itr++; 
            if(F2max<F2conv){ 
                niter=0; 
                T1 = (getCPUticks()-T00)*tick2second;
                if(verbosity>0)printf( "run_omp(nSys=%i|iPara=%i,bOcl=%i,bGridFF=%i) CONVERGED in %i/%i nsteps |F|=%g time=%g[ms] %g[us/step] \n", nSystems, iParalel,bOcl,bGridFF, itr,niter_max, sqrt(F2max), T1*1000, T1*1e+6/itr );
                itr--;
            }
        }
        }
    } 
    for(int i=0; i<ffl.nvecs; i++){
        ffl.apos [i]=ffls[iSystemCur].apos [i];
        ffl.fapos[i]=ffls[iSystemCur].fapos[i];
    }
    T1 = (getCPUticks()-T00)*tick2second;
    if(itr>=niter_max)if(verbosity>0)printf( "run_omp(nSys=%i|iPara=%i,bOcl=%i,bGridFF=%i) NOT CONVERGED in %i/%i nsteps |F|=%g time=%g[ms] %g[us/step]\n", nSystems, iParalel,bOcl,bGridFF, itr,niter_max, sqrt(F2max), T1*1000, T1*1e+6/niter_max );
    return itr;
}



int run_omp_ocl( int niter_max, double Fconv=1e-3, double Flim=1000, double timeLimit=0.02 ){
    //printf( "run_omp_ocl() niter_max=%i %dt=g Fconv%g \n", niter_max, dt, Fconv  ); 
    double F2conv=Fconv*Fconv;
    double F2max=0;
    int itr=0,niter=niter_max;
    long T0,T00;
    double T1,T2;
    int err=0;
    T00 = getCPUticks();
    //#pragma omp parallel shared(E,F2,ff,vv,vf,ffl) private(itr)
    #pragma omp parallel shared(niter,itr,F2max,T0,T1,T2,err,gopt_ifound)
    while(itr<niter){
        if(itr<niter){
        #pragma omp single
        {
        if(bOcl){
            T0 = getCPUticks();
            ocl.upload( ocl.ibuff_atoms, atoms  );
            if( bSurfAtoms ){
                if(bGridFF){ 
                    err |= task_NBFF_Grid ->enque_raw(); 
                }else{
                    err |= task_NBFF     ->enque_raw(); 
                    err |= task_SurfAtoms->enque_raw(); 
                }
            }else{ err |= task_NBFF      ->enque_raw(); }
            ocl.download( ocl.ibuff_aforces , aforces );
            //ocl.download( ocl.ibuff_avel    , avel    );

            
        }
        }
        #pragma omp for
        for( int isys=0; isys<nSystems; isys++ ){
            if(bGopt)gopts[isys].update();    
            //printf( "run_omp_ocl[itr=%i] isys=%i @cpu(%i/%i) \n", itr, isys, omp_get_thread_num(), omp_get_num_threads() );
            ffls[isys].cleanForce();
            ffls[isys].eval(false);
            if(!bOcl){
                if(bPBC){ ffls[isys].Etot += ffls[isys].evalLJQs_ng4_PBC_simd(); }
                else    { ffls[isys].Etot += ffls[isys].evalLJQs_ng4_simd    (); } 
                //ffls[isys].evalLJQs_ng4_simd    ();
                if(bSurfAtoms)[[likely]]{ 
                    if   (bGridFF)[[likely]]{ gridFF.addForces            ( ffls[isys].natoms, ffls[isys].apos, ffls[isys].PLQs, ffls[isys].fapos ); }        // GridFF
                    else                    { gridFF.evalMorsePBCatoms_sym( ffls[isys].natoms, ffls[isys].apos, ffls[isys].REQs, ffls[isys].fapos ); }
                }
            }
            if( (ipicked>=0) && (isys==iSystemCur) ){ 
                pullAtom( ipicked, ffls[isys].apos, ffls[isys].fapos );  
            };
            if(bConstrains){ // ToDo: this apply the same constrains to all replicas, we may want to apply different constrains to different replicas
                //printf( "run_omp() constrs[%i].apply()\n", constrs.bonds.size() );
                ffls[isys].Etot += constrs.apply( ffls[isys].apos, ffls[isys].fapos, &ffls[isys].lvec );
            }
            //printf("ffls[%i].Etot(itr=%i) = %g \n", isys, itr, ffls[isys].Etot );
            //if( bGopt && gopts[isys].bExploring ){  ffls[isys].Etot += gopts[isys].constrs.apply( ffls[isys].apos, ffls[isys].fapos, &ffls[isys].lvec ); }
        }        
        #pragma omp barrier
        #pragma omp single
        {   F2max=0;
            if(bOcl){
            //T1 = (getCPUticks()-T0)*tick2second;
            err |= ocl.finishRaw(); OCL_checkError(err, "eval_MMFFf4_ocl.finish");
            //T2 = (getCPUticks()-T0)*tick2second;
            //printf( "CPU.finish %g[us] GPU.finish %g[us] \n", T1*1e+6, T2*1e+6 );
        }}
        #pragma omp for reduction(max:F2max)
        for( int isys=0; isys<nSystems; isys++ ){
            int i0v = isys*ocl.nvecs;
            //if(bOcl)unpack_add( ffls[isys].natoms, ffls[isys].fapos, aforces+i0v );  // APPLY non-covalent forces form GPU
            double F2=1.0;
            if( gopts[isys].bExploring ){
                //printf( "run_omp_ocl() gopts[isys].bExploring=true \n", isys );
                ffls[isys].Etot += gopts[isys].constrs.apply( ffls[isys].apos, ffls[isys].fapos, &ffls[isys].lvec );
                ffls[isys].move_Langevin( opts[isys].dt_max, 10000.0, gopts[isys].gamma_damp, gopts[isys].T_target );
                if(bToCOG){ Vec3d cog=average( ffls[isys].natoms, ffls[isys].apos );  move( ffls[isys].natoms, ffls[isys].apos, cog*-1.0 ); }
            }
            else{
                F2 = opts[isys].move_FIRE();
            }
            //double F2 = opts[isys].move_FIRE();
            // if(isys==iSystemCur){
            //     double cv;
            //     double cf  = opts[isys].renorm_vf * opts[isys].damp_func_FIRE( opts[isys].cos_vf, cv );
            //     printf( "opt[%i](dt=%g/%g,damp=%g,cf=%g,cv=%g)cos_vf=%g\n", iSystemCur, opts[isys].dt,opts[isys].dt_max,opts[isys].damping,opts[isys].cv,opts[isys].cf, opts[isys].cos_vf );
            // }
            F2max = fmax(F2max,F2);
            //if(bOcl)
            pack( ffls[isys].nvecs,  ffls[isys].apos, atoms  +i0v );
        }
        #pragma omp single
        { 
            nloop++;
            itr++; 
            // if(timeLimit>0){
            //     double t = (getCPUticks() - T0)*tick2second;
            //     if(t>0.02){ 
            //         niter=0; 
            //         if(verbosity>0)printf( "run_omp() ended due to time limit after %i nsteps ( %6.3f [s]) \n", itr, t ); 
            //     }
            // }

            if( bGopt ){
                for(int isys=0; isys<nSystems; isys++){
                    if( ( opts[isys].ff < F2conv )  && ( !gopts[isys].bExploring ) ){
                        gopt_ifound++;
                        sprintf(tmpstr,"# %i sys: %i E: %g |F|: %g istep: %i", gopt_ifound, isys, ffls[isys].Etot, sqrt(ffl.cvf.z), gopts[isys].istep );
                        printf( "run_omp_ocl(GOopt).save %s \n", tmpstr );
                        saveSysXYZ( isys, "gopt.xyz", tmpstr, false, "a", nPBC_save );
                        gopts[isys].startExploring();
                        gopts[isys].apply_kick( ffl.natoms, ffl.apos, ffl.vapos );
                        if( gopt_ifound > gopt_nfoundMax ){ printf( "gopt_ifound(%i)>gopt_nfoundMax(%i) => exit \n", gopt_ifound, gopt_nfoundMax ); exit(0); };
                    }
                }
            }else if(F2max<F2conv){ 
                niter=0; 
                T1 = (getCPUticks()-T00)*tick2second;
                //printf( "run_omp_ocl(nSys=%i|iPara=%i) CONVERGED in %i/%i steps, |F|(%g)<%g time %g[ms] %g[us/step] bGridFF=%i \n", nSystems, iParalel, niterdone,niter, sqrt(F2), Fconv, t*1000, t*1e+6/niterdone, bGridFF ); 
                if(verbosity>0)printf( "run_omp_ocl(nSys=%i|iPara=%i,bOcl=%i,bGridFF=%i) CONVERGED in %i/%i nsteps |F|=%g time=%g[ms] %g[us/step] \n", nSystems, iParalel,bOcl,bGridFF, itr,niter_max, sqrt(F2max), T1*1000, T1*1e+6/itr );
                itr--;
            }
            //printf( "step[%i] E %g |F| %g ncpu[%i] \n", itr, E, sqrt(F2), omp_get_num_threads() ); 
            //{printf( "step[%i] dt %g(%g) cv %g cf %g cos_vf %g \n", itr, opt.dt, opt.dt_min, opt.cv, opt.cf, opt.cos_vf );}
            //if(verbosity>2){printf( "step[%i] E %g |F| %g ncpu[%i] \n", itr, E, sqrt(F2), omp_get_num_threads() );}
        }
        }
    } // END while(itr<niter) OPENMP_BLOCK
    //printf( "run_omp_ocl().copy iSystemCur=%i \n", iSystemCur );
    for(int i=0; i<ffl.nvecs; i++){
        ffl.apos [i]=ffls[iSystemCur].apos [i];
        ffl.fapos[i]=ffls[iSystemCur].fapos[i];
    }
    T1 = (getCPUticks()-T00)*tick2second;
    if(itr>=niter_max)if(verbosity>0)printf( "run_omp_ocl(nSys=%i|iPara=%i,bOcl=%i,bGridFF=%i) NOT CONVERGED in %i/%i nsteps |F|=%g time=%g[ms] %g[us/step]\n", nSystems, iParalel,bOcl,bGridFF, itr,niter_max, sqrt(F2max), T1*1000, T1*1e+6/niter_max );
    return itr;
}

int run_multi_serial( int niter_max, double Fconv=1e-3, double Flim=1000, double timeLimit=0.02 ){
    //printf( "run_omp_ocl() niter_max=%i %dt=g Fconv%g \n", niter_max, dt, Fconv  ); 
    double F2conv=Fconv*Fconv;
    double F2max=0;
    int itr=0,niter=niter_max;
    long T0,T00;
    double T1,T2;
    int err=0;
    T00 = getCPUticks();
    while(itr<niter){
        if(itr<niter){
        F2max=0;
        for( int isys=0; isys<nSystems; isys++ ){
            ffls[isys].cleanForce();
            ffls[isys].eval(false);
            if(bPBC){ ffls[isys].evalLJQs_ng4_PBC_simd(); }
            else    { ffls[isys].evalLJQs_ng4_simd    (); } 
            if( (ipicked>=0) && (isys==iSystemCur) ){ 
                pullAtom( ipicked, ffls[isys].apos, ffls[isys].fapos );  
            };
            double F2 = opts[isys].move_FIRE();
            F2max = fmax(F2max,F2);
        }        
        nloop++;
        itr++; 
        if(F2max<F2conv){ 
            niter=0; 
            T1 = (getCPUticks()-T00)*tick2second;
            if(verbosity>0)printf( "run_multi_serial(nSys=%i|iPara=%i) CONVERGED in %i/%i nsteps |F|=%g time=%g[ms] %g[us/step]\n", nSystems,iParalel, itr,niter_max, sqrt(F2max), T1*1000, T1*1e+6/itr );
            itr--;
        }
        }
    } // END while(itr<niter) OPENMP_BLOCK
    //printf( "run_omp_ocl().copy iSystemCur=%i \n", iSystemCur );
    for(int i=0; i<ffl.nvecs; i++){
        ffl.apos [i]=ffls[iSystemCur].apos [i];
        ffl.fapos[i]=ffls[iSystemCur].fapos[i];
    }
    T1 = (getCPUticks()-T00)*tick2second;
    if(itr>=niter_max)if(verbosity>0)printf( "run_multi_serial(nSys=%i|iPara=%i) NOT CONVERGED in %i/%i nsteps |F|=%g time=%g[ms] %g[us/step]\n", nSystems,iParalel, itr,niter_max, sqrt(F2max), T1*1000, T1*1e+6/niter_max );
    return itr;
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
    if(bConstrains)constrs.apply( nbmol.apos, nbmol.fapos, &ffl.lvec );
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

virtual char* getStatusString( char* s, int nmax ) override {
    //double F2 = evalVFs();
    double  F2max=0;
    double  F2min=1e+300;
    for(int isys=0; isys<nSystems; isys++){
        int i0v = isys * ocl.nvecs;
        //evalVF( ocl.nvecs, aforces+i0v, avel   +i0v, fire[isys], MDpars[isys] );
        F2max = fmax( F2max, fire[isys].ff );
        F2min = fmin( F2min, fire[isys].ff );
    }
    s += sprintf(s, "iSystemCur %i/%i \n",  iSystemCur, nSystems );
    s += sprintf(s, "eval_MMFFf4_ocl |F|max=%g |F|min=%g \n", sqrt(F2max), sqrt(F2min) );
    return s;
}

virtual void MDloop( int nIter, double Ftol = -1 ) override {
    if(Ftol<0)  Ftol = Ftol_default;
    //printf( "MolWorld_sp3_ocl::MDloop(%i) bGridFF %i bOcl %i bMMFF %i \n", nIter, bGridFF, bOcl, bMMFF );
    //bMMFF=false;


    if(bAnimManipulation){
        float dir_sign = ((nloop/10000)%2)*2 - 1 ; //printf("AnimManipulation dir_sign=%f \n", dir_sign);
        move_MultiConstrain( Vec3d{ 0.2*dir_sign, 0.0, 0.0}, Vec3dZero );
    }

    // //run_omp_ocl( 1 );
    // run_omp_ocl( 50 );
    // //run_omp_ocl( 100 );
    // return;

    nIter=iterPerFrame;
    //bOcl=iParalel>0;
    int nitrdione=0;
    switch(iParalel){
        case -1: bOcl=false;  nitrdione = run_multi_serial(nIter,Ftol);  break; 
        case  0:
        case  1: bOcl=false;  nitrdione = run_omp_ocl( nIter, Ftol    ); break; 
        case  2: bOcl=true;   nitrdione = run_omp_ocl( nIter, Ftol    ); break; 
        case  3: bOcl=true;   nitrdione = run_ocl_opt( nIter, Ftol    ); break; 
        //case  3: bOcl=true; nitrdione = run_ocl_loc( nIter, Ftol, 1 ); break; 
        //case  4: bOcl=true; nitrdione = run_ocl_loc( nIter, Ftol, 2 ); break; 
        default:
            eval_NBFF_ocl_debug();
    }

    unpack_system( iSystemCur, ffl, true, true );

    /*
    if( bOcl ){
        //printf( "GPU frame[%i] -- \n", nIter );
        if( (iSystemCur<0) || (iSystemCur>=nSystems) ){  printf("ERROR: iSystemCur(%i) not in range [ 0 .. nSystems(%i) ] => exit() \n", iSystemCur, nSystems ); exit(0); }
        nIter = 50;
        //nIter = 1;
        eval_MMFFf4_ocl( nIter );
        unpack_system(iSystemCur, ffl, true, true); 
        //SDL_Delay(1000);
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
            if(ipicked>=0){ pullAtom( ipicked, nbmol.apos, nbmol.fapos );  }; // printf( "pullAtom(%i) E=%g\n", ipicked, E ); };
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
    */
    
    icurIter+=nitrdione;
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
    pack( gridFF.natoms, gridFF.apos, atoms_surf, sq(gridFF.Rdamp) );
    pack( gridFF.natoms, gridFF.REQs, REQs_surf                     );   // ToDo: H-bonds should be here
    long T1=getCPUticks();
    ocl.GFFparams.x = gridFF.Rdamp;
    ocl.GFFparams.y = gridFF.alphaMorse;
    //v2f4( gridFF.grid.pos0, ocl.grid_p0 );
    ocl.grid_p0.f     = (Vec3f)gridFF.grid.pos0;
    ocl.grid_shift0.f = (Vec3f)gridFF.shift0;
    //printf( "grid_p0(%g,%g,%g) grid_shift0(%g,%g,%g)\n", ocl.grid_p0.x,ocl.grid_p0.y,ocl.grid_p0.z,    ocl.grid_shift0.x,ocl.grid_shift0.y,ocl.grid_shift0.z  );
    //exit(0);
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
    //  ToDo: We do no release these buffers because we use them for kernell:: getSurfMorse()
    //ocl.buffers[ocl.ibuff_atoms_surf].release();
    //ocl.buffers[ocl.ibuff_REQs_surf ].release();
    err |=  ocl.finishRaw();    OCL_checkError(err, "surf2ocl.makeGridFF.atoms_surf.release");
    bDONE_surf2ocl = true;
    //exit(0);
}

virtual void initGridFF( const char * name, bool bGrid=true, bool bSaveDebugXSFs=false, double z0=NAN, Vec3d cel0={-0.5,-0.5,0.0}, bool bAutoNPBC=true )override{
    printf( "MolWorld_sp3_multi::initGridFF() \n");
    if(verbosity>0)printf("MolWorld_sp3_multi::initGridFF(%s,bGrid=%i,z0=%g,cel0={%g,%g,%g})\n",  name, bGrid, z0, cel0.x,cel0.y,cel0.z  );
    if(gridFF.grid.n.anyEqual(0)){ printf("ERROR in MolWorld_sp3_multi::initGridFF() zero grid.n(%i,%i,%i) => Exit() \n", gridFF.grid.n.x,gridFF.grid.n.y,gridFF.grid.n.z ); exit(0); };
    gridFF.grid.center_cell( cel0 );
    bGridFF=true;
    gridFF.bindSystem        ( surf  .natoms, surf  .atypes, surf.apos  , surf.REQs        );
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
    gridFF.makePBCshifts   ( gridFF.nPBC,      gridFF.lvec       );
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
    //ocl.download( ocl.ibuff_aforces, ff4.fapos, ff4.nvecs, ff4.nvecs*iSystemCur );
    //ocl.download( ocl.ibuff_atoms,   ff4.apos , ff4.nvecs, ff4.nvecs*iSystemCur );
    //for(int i=0; i<ff4.nvecs; i++){  printf("CPU[%i] p(%g,%g,%g) f(%g,%g,%g) pi %i \n", i, ff4.apos[i].x,ff4.apos[i].y,ff4.apos[i].z,  ff4.fapos[i].x,ff4.fapos[i].y,ff4.fapos[i].z, i>=ff4.natoms ); }
    err |= ocl.finishRaw();
    OCL_checkError(err, "eval_MMFFf4_ocl_debug.2");

    //for(int i=0; i<ff4.nvecs; i++){  printf("OCL[%4i] f(%10.5f,%10.5f,%10.5f) p(%10.5f,%10.5f,%10.5f) pi %i \n", i, ff4.fapos[i].x,ff4.fapos[i].y,ff4.fapos[i].z,  ff4.apos[i].x,ff4.apos[i].y,ff4.apos[i].z,  i>=ff4.natoms ); }

    //ffl.eval();   // We already computed this above
    bool ret=false;
    //printf("### Compare ffl.fapos,  GPU.fapos  \n"); ret |= compareVecs( ffl.natoms, ffl.fapos,  ff4.fapos,  1e-4, true );
    //printf("### Compare ffl.fpipos, GPU.fpipos \n"); ret |= compareVecs( ffl.nnode,  ffl.fpipos, ff4.fpipos, 1e-4, true ); 
    //if(ret){ printf("ERROR: GPU.eval() and ffl.eval() produce different results => exit() \n"); exit(0); }else{ printf("CHECKED: GPU task_MMFF.eval() == CPU ffl.eval() \n"); }

    //printf("GPU AFTER assemble() \n"); ff4.printDEBUG( false,false );
    //unpack( ffl.natoms, ffl.  apos, ff4.  apos );
    //unpack( ffl.natoms, ffl. fapos, ff4. fapos );
    //unpack( ffl.nnode,  ffl. pipos, ff4. pipos );
    //unpack( ffl.nnode,  ffl.fpipos, ff4.fpipos );

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

double eval_NBFF_ocl_debug( bool bCompareToCPU=true, bool bMove=false, bool bPrintOnGPU=false, bool bPrintOnCPU=false ){ 
    printf("MolWorld_sp3_multi::eval_NBFF_ocl_debug() \n");
    int err=0;
    if( task_NBFF==0 ){ setup_NBFF_ocl(); }

    if(bCompareToCPU){ // --- evaluate on CPU
        unpack_system( iSystemCur, ffl );
        ffl  .cleanForce();
        //nbmol.evalLJQs_ng4_PBC( ffl.neighs, ffl.neighCell, ffl.lvec, ffl.nPBC, gridFF.Rdamp );
        //nbmol.evalLJQs_ng4_PBC( ffl.neighs, ffl.neighCell, npbc, pbc_shifts, gridFF.Rdamp ); 
        nbmol.evalLJQs_ng4( ffl.neighs, gridFF.Rdamp );
        fcog  = sum ( ffl.natoms, ffl.fapos   );
        tqcog = torq( ffl.natoms, ffl.apos, ffl.fapos );
        if(  fcog.norm2()>1e-8 ){ printf("WARRNING: eval_NBFF_ocl_debug() CPU |fcog| =%g; fcog=(%g,%g,%g)\n", fcog.norm(),  fcog.x, fcog.y, fcog.z ); exit(0); }
    }
    
    { // --- evaluate on GPU
        err |= task_cleanF->enque_raw(); // this should be solved inside  task_move->enque_raw();
        err |= task_NBFF  ->enque_raw();
        if(bPrintOnGPU) err |= task_print ->enque_raw(); // just printing the forces before assempling
        if(bMove      ) err |= task_move  ->enque_raw();
    }
    
    ocl.download( ocl.ibuff_aforces, ff4.fapos, ff4.nvecs, ff4.nvecs*iSystemCur );
    ocl.download( ocl.ibuff_atoms,   ff4.apos , ff4.nvecs, ff4.nvecs*iSystemCur );
    err |=  ocl.finishRaw();
    OCL_checkError(err, "eval_NBFF_ocl_debug");

    if(bPrintOnCPU)for(int i=0; i<ff4.nvecs; i++){  printf("OCL[%4i] f(%10.5f,%10.5f,%10.5f) p(%10.5f,%10.5f,%10.5f) pi %i \n", i, ff4.fapos[i].x,ff4.fapos[i].y,ff4.fapos[i].z,  ff4.apos[i].x,ff4.apos[i].y,ff4.apos[i].z,  i>=ff4.natoms ); }

    if( ckeckNaN_f( ff4.nvecs, 4, (float*)ff4.fapos, "gpu.fapos" ) ){ printf("ERROR in MolWorld_sp3_multi::eval_NBFF_ocl_debug() fapos contain NaNs \n"); exit(0); };

    if(bCompareToCPU){
        // ---- Compare to ffl
        bool ret=false;
        printf("### Compare ffl.fapos,  ff4.fapos   \n"); ret |= compareVecs( ff4.natoms, ffl.fapos,  ff4.fapos,  1e-6, 1 );
        if(ret){ printf("ERROR: GPU task_NBFF.eval() != ffl.nbmol.evalLJQs_ng4_PBC() => exit() \n"); exit(0); }else{ printf("CHECKED: GPU task_NBFF.eval() == ffl.nbmol.evalLJQs_ng4_PBC() \n"); }
    }

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
