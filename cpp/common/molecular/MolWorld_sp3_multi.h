
#ifndef MolWorld_sp3_ocl_h
#define MolWorld_sp3_ocl_h

#include <string>
#include "MMFFf4.h"
//#include "Molecule.h"
//#include "MMFFBuilder.h"
#include "MolWorld_sp3.h"
//#include "OCL_DFT.h"
//#include "OCL_PP.h"
#include "OCL_MM.h"
#include "OCL_UFF.h"
//#include "OCL_EFF.h"
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
    OCL_UFF*   uff_ocl = nullptr;
   // bool       bUFF    = false;  -- this is inherited from MolWorld_sp3

    FIRE_setup fire_setup{0.1,0.1};

    int iParalel_default=3;
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
    int*    excl       =0;
    Quat4i* bkNeighs   =0;
    Quat4i* bkNeighs_new=0;
    //Quat4f* neighForce =0;

    Quat4f* REQs       =0;
    Quat4f* MMpars     =0;
    Quat4f* BLs        =0;
    Quat4f* BKs        =0;
    Quat4f* Ksp        =0;
    Quat4f* Kpp        =0;

    // ---- UFF Host-side parameter arrays
    // Topology/params per system packed consecutively across systems
    std::vector<int2>    host_bon_atoms;
    std::vector<float2>  host_bon_params;
    std::vector<int4>    host_ang_atoms;
    std::vector<int2>    host_ang_ngs;
    std::vector<float4>  host_ang_params1;
    std::vector<float>   host_ang_params2_w;
    std::vector<int4>    host_dih_atoms;
    std::vector<int4>    host_dih_ngs;
    std::vector<float4>  host_dih_params;
    std::vector<int4>    host_inv_atoms;
    std::vector<int4>    host_inv_ngs;
    std::vector<float4>  host_inv_params;
    std::vector<int4>    host_neighs_UFF;     // [nSystems*nAtoms]
    std::vector<int4>    host_neighCell_UFF;  // [nSystems*nAtoms]
    std::vector<int4>    host_neighBs_UFF;    // [nSystems*nAtoms]
    std::vector<int>     host_a2f_offsets;
    std::vector<int>     host_a2f_counts;
    std::vector<int>     host_a2f_indices;


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
    OCLtask* task_NBFF_ex2=0;
    OCLtask* task_NBFF_Grid=0;
    OCLtask* task_NBFF_Grid_Bspline=0;
    OCLtask* task_NBFF_Grid_Bspline_tex=0;
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

    int nPerVFs = 10;

    GOpt*         gopts=0;

    MMFFf4     ff4;

    Quat4f* afm_ps=0;

    const char* uploadPopName=0;

    bool bMILAN = false;
    bool bSaveToDatabase=false;

    // For deferred GridFF initialization
    bool bGridFF_pending = false;
    const char* gridFF_name = nullptr;
    double gridFF_z0 = NAN;
    Vec3d gridFF_cel0 = {-0.5,-0.5,0.0};
    bool gridFF_bSymetrize = true;
    bool gridFF_bAutoNPBC = true;
    bool gridFF_bCheckEval = true;
    bool gridFF_bUseEwald = true;
    bool gridFF_bFit = true;
    bool gridFF_bRefine = true;

    MolecularDatabase* database = 0;
    long nStepConvSum = 0;
    long nStepNonConvSum = 0;
    long nStepExplorSum = 0;
    int  nbConverged=0;
    int  nbNonConverged=0;
    int  nbEvaluation=0;
    int  nExploring=0;

virtual int getMolWorldVersion() const override { return (int)MolWorldVersion::GPU; };

// ==================================
//         Initialization
// ==================================

void realloc( int nSystems_ ){
    DEBUG;
    printf("MolWorld_sp3_multi::realloc() \n");
    nSystems=nSystems_;

    if(bUFF){
        printf( "MolWorld_sp3_multi::realloc() UFF Systems %i nAtoms %i \n", nSystems, ffu.natoms );
        int nAtoms      = ffu.natoms;
        int nBonds      = ffu.nbonds;
        int nAngles     = ffu.nangles;
        int nDihedrals  = ffu.ndihedrals;
        int nInversions = ffu.ninversions;
        // The size of the atom-to-force mapping array is the total number of force components
        int nA2F        = ffu.nf;
        uff_ocl->realloc( nSystems, ffu.natoms, ffu.nbonds, ffu.nangles, ffu.ndihedrals, ffu.ninversions, 0, nA2F );
        // Allocate host-side arrays
        host_bon_atoms    .resize(nSystems * ffu.nbonds);
        host_bon_params   .resize(nSystems * ffu.nbonds);
        host_ang_atoms    .resize(nSystems * ffu.nangles);
        host_ang_ngs      .resize(nSystems * nAngles);
        host_ang_params1  .resize(nSystems * nAngles);
        host_ang_params2_w.resize(nSystems * nAngles);
        host_dih_atoms    .resize(nSystems * nDihedrals);
        host_dih_ngs      .resize(nSystems * nDihedrals);
        host_dih_params   .resize(nSystems * nDihedrals);
        host_inv_atoms    .resize(nSystems * nInversions);
        host_inv_ngs      .resize(nSystems * nInversions);
        host_inv_params   .resize(nSystems * nInversions);
        host_neighs_UFF   .resize(nSystems * nAtoms);
        host_neighCell_UFF.resize(nSystems * nAtoms);
        host_neighBs_UFF  .resize(nSystems * nAtoms);
        host_a2f_offsets  .resize(nSystems * (nAtoms+1));
        host_a2f_counts   .resize(nSystems * nAtoms);
        host_a2f_indices  .resize(nSystems * nA2F);
        _realloc ( atoms,     uff_ocl->nAtomsTot  );
        _realloc0( aforces,   uff_ocl->nAtomsTot  , Quat4fZero );
        _realloc0( avel,      uff_ocl->nAtomsTot  , Quat4fZero );
        _realloc( REQs,       uff_ocl->nAtomsTot );
        _realloc( pbcshifts, uff_ocl->nSystems * (npbc+1) ); // npbc is from MolWorld_sp3
    }else{
        printf( "MolWorld_sp3_multi::realloc() MMFF Systems %i nAtoms %i nnode %i \n", nSystems, ffl.natoms,  ffl.nnode );
        ocl.initAtomsForces( nSystems, ffl.natoms,  ffl.nnode, npbc+1 );
        _realloc ( atoms,     ocl.nvecs*nSystems  );
        _realloc0( aforces,   ocl.nvecs*nSystems  , Quat4fZero );
        _realloc0( avel,      ocl.nvecs*nSystems  , Quat4fZero );
        _realloc0( cvfs,      ocl.nvecs*nSystems  , Quat4fZero );
        _realloc0( constr,    ocl.nAtoms*nSystems , Quat4fOnes*-1. );
        _realloc0( constrK,   ocl.nAtoms*nSystems , Quat4fOnes*-1. );
        _realloc( neighs,    ocl.nAtoms*nSystems );
        _realloc( neighCell, ocl.nAtoms*nSystems );
        _realloc( excl,      ocl.nAtoms*nSystems*EXCL_MAX );
        _realloc( bkNeighs,    ocl.nvecs*nSystems );
        _realloc( bkNeighs_new,ocl.nvecs*nSystems );
        _realloc( REQs,      ocl.nAtoms*nSystems );
        _realloc( MMpars,    ocl.nnode*nSystems  );
        _realloc( BLs,       ocl.nnode*nSystems  );
        _realloc( BKs,       ocl.nnode*nSystems  );
        _realloc( Ksp,       ocl.nnode*nSystems  );
        _realloc( Kpp,       ocl.nnode*nSystems  );
        _realloc( pbcshifts, ocl.npbc*nSystems );
    }

    // Common buffers for both FF
    cl_Mat3 bbox0{cl_float4{-1e+8,-1e+8,-1e+8,-1e+8,},cl_float4{+1e+8,+1e+8,+1e+8,+1e+8,}, cl_float4{-1.,-1.,-1.,-1.} };
    _realloc0( bboxes,   nSystems, bbox0 );
    _realloc( lvecs,     nSystems  );
    _realloc( ilvecs,    nSystems  );
    _realloc( MDpars,    nSystems  );
    Quat4f TDrive0{0.0,-1.0,0.0,0.0};
    _realloc0( TDrive,   nSystems, TDrive0 );
    _realloc( fire,      nSystems  );
}

void initMultiCPU(int nSys){
    printf("MolWorld_sp3_multi::initMultiCPU(%i)\n", nSys);
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
        gopts[isys].nExplore = 1000/nPerVFs;
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

    // Temporarily set surf_name to null to prevent base class from calling initGridFF
    const char* original_surf_name = surf_name;
    surf_name = nullptr;

    // --- Set flags for base class BEFORE calling its init()
    printf("DEBUG MolWorld_sp3_multi::init() BEFORE base init, this->bUFF=%d, MolWorld_sp3::bUFF=%d\n", this->bUFF, MolWorld_sp3::bUFF);

    MolWorld_sp3::init(); // loading of .xyz and .mol/.mo2 is done here, look what is done hare and respcet it

    printf("DEBUG MolWorld_sp3_multi::init() AFTER base init, &ffu=%p\n", &ffu);

    // Restore surf_name
    surf_name = original_surf_name;

    if(bUFF){
        printf("MolWorld_sp3_multi::init() for UFF_OCL\n");
        ffu.Rdamp = 1.0;
        ffu.FmaxNonBonded = 1000.0;
        ffu.SubNBTorsionFactor = 0.0;
        uff_ocl = new OCL_UFF();
        uff_ocl->init(i_nvidia); // Initialize the uff_ocl object
        uff_ocl->makeKernels("common_resources/cl" );
        // The base class MolWorld_sp3::init() already takes care of building the UFF force field (ffu)
        // when bUFF is set. We just need to initialize the OpenCL part here.
        // The redundant calls below were likely causing inconsistencies.
    } else {
        ocl.init(i_nvidia);
        ocl.makeKrenels_MM("common_resources/cl" );
        // initialization of ff4 is here because parrent MolWorld_sp3 does not contain ff4 anymore
        builder.toMMFFf4( ff4, true, bEpairs );  //ff4.printAtomParams(); ff4.printBKneighs();
        ff4.flipPis( Vec3fOne );
        ff4.setLvec((Mat3f)builder.lvec);
        ff4.makeNeighCells  ( nPBC );
    }

    // Handle surface loading and GridFF initialization if we have a surface
    if(original_surf_name != nullptr){
        printf("MolWorld_sp3_multi::init() - Loading surface %s after OpenCL initialization\n", original_surf_name);
        loadSurf(original_surf_name, true, false, NAN);
    }

    // ----- init systems
    realloc( nSystems );
    // Prepare UFF kernels after sizes are known
    if(bUFF && uff_ocl){uff_ocl->setup_kernels( ffu.Rdamp, ffu.FmaxNonBonded, ffu.SubNBTorsionFactor); }
    //if(bGridFF) evalCheckGridFF_ocl();  // this must be after we make buffers but before we fill them
    float random_init = 0.5;

    fire_setup.setup(opt.dt_max,opt.damp_max);
    for(int i=0; i<nSystems; i++){
        if(bUFF) pack_uff_system( i, ffu, true, false, true, true ); else pack_system( i, ffl, true, false, true, true, random_init );
        fire[i].bind_params( &fire_setup );
        fire[i].id=i;
    }
    initMultiCPU(nSystems);
    //for(int i=0;i<nSystems; i++){ printMSystem(i); }
    upload( true, false, true, true );
    //bGridFF=false;
    bOcl   =false;
    //bOcl   =true;

    // Complete GridFF initialization if it was deferred
    if(bGridFF_pending){
        completeGridFFInit();
    }
    //setup_MMFFf4_ocl();
    int4 mask{1,1,0,0};
    //ocl.printOnGPU( 0,mask );
    //ocl.printOnGPU( 1,mask );
    //ocl.printOnGPU( 2,mask );

    iParalelMax= 4;
    iParalelMin=-1;
    //iParalel=1;
    //iParalel=3;
    //iParalel=iParalelMax;

    if(database)
    if(!database){
        database = new MolecularDatabase();
        database->setDescriptors();
    }

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
            //printf( "CPU [isys=%i,ig=%i]\n", isys, ig );
            for(int j=0; j<g.i0n.y; j++){
                int ia = groups.g2a[ j + g.i0n.x ];
                //int ia = g2a[ i0g + j ];
                int iia = i0a + ia;
                gweights [ iia ] = groups. weights[ia];
                gfweights[ iia ] = groups.fweights[ia];
                //printf( "CPU gweights[%i](%g,%g,%g,%g) gfweights(%g,%g)\n", iia, gweights[iia].x,gweights[iia].y,gweights[iia].z,gweights[iia].w,  gfweights[iia].x,gfweights[iia].y );
            }
        }
    }
}


int init_groups(){
    int ngroup = groups.groups.size();
    printf("MolWorld_sp3_multi::init_groups() ngroup=%i\n", ngroup );
    // int ngroup = -1;
    // printf("atom2group.size()==%i\n", atom2group.size() );
    // for(int i=0; i<atom2group.size(); i++){
    //     //printf("atom2group[%i]==%i\n", i, atom2group[i]);
    //     ngroup=_max(ngroup,atom2group[i]);
    // }
    // ngroup+=1;
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
    err|= ocl.upload( ocl.ibuff_gforces,   gforces );   OCL_checkError(err, "init_groups.upload(gforces)");
    err|= ocl.upload( ocl.ibuff_gtorqs,    gtorqs  );   OCL_checkError(err, "init_groups.upload(gtorqs)");
    err|= ocl.upload( ocl.ibuff_bboxes,    bboxes  );   OCL_checkError(err, "init_groups.upload(bboxes)");  // ToDo: this may need to go to different place
    return err;
}

virtual void pre_loop() override {
    printf("MolWorld_sp3_multi::pre_loop()\n" );
    init_groups();
    for(int isys=0; isys<nSystems; isys++){
        for(int ia=0; ia<ffl.natoms; ia++ ){
            int iaa = isys * ocl.nAtoms + ia;
            //constr [iaa] = (Quat4f)ffl.constr[ia];
            //constrK[iaa] = (Quat4f)ffl.constrK[ia];
            ffls[isys].constr [ia] = ffl.constr [ia];
            ffls[isys].constrK[ia] = ffl.constrK[ia];
            //if(isys==0){    printf( "pre_loop() ffl[ia=%i] constr(%g,%g,%g|%g) constrK(%g,%g,%g) \n", isys, ia, ffl.constr[ia].x,ffl.constr[ia].y,ffl.constr[ia].z,ffl.constr[ia].w, ffl.constrK[ia].x,ffl.constrK[ia].y,ffl.constrK[ia].z ); }
        }
    }
    //printConstrains();
    // for(int ic : constrain_list ){
    //     for(int isys=0; isys<nSystems; isys++){
    //         int i0a   = isys * ocl.nAtoms;
    //         int i0v   = isys * ocl.nvecs;
    //         Quat4f a=atoms[ic+i0v]; a.w=1.0; constr[ic+i0a] = a;
    //     }
    // }
    // ocl.upload( ocl.ibuff_constr, constr );     //OCL_checkError(err, "init_groups.upload(constr)");
    print("MolWorld_sp3_multi::pre_loop() DONE\n");
}

// ==================================
//         GPU <-> CPU I/O  MMFF
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
    for(int i=0; i<ocl.nAtoms; i++){
        Quat4f a=atoms[i+i0v];
        //a.w=-1.0;
        a.w = ffl.constr[i].w;
        constr [i+i0a] = a;
        constrK[i+i0a] = (Quat4f)ffl.constrK[i];
        //printf(  "atom[%3i|sys=%i](%8.5f,%8.5f,%8.5f) constr(%8.5f,%8.5f,%8.5f|%8.5f)\n", i, isys, ff.apos[i].x,ff.apos[i].y,ff.apos[i].z,     a.x,a.y,a.z,a.w );
    }  // contrains
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

    //Mat3d lvec0=ff.lvec;
    if(bForces){ pack( ff.nvecs, ff.fapos, aforces+i0v ); }
    //if(bVel   ){ pack( ff.nvecs, opt.vel,  avel   +i0v ); }
    if(bParams){


        //Mat3d lvec=lvec0; lvec.a.addRandomCube(1.5); ffl.setLvec(lvec);
        //Mat3d lvec=lvec0; lvec.a=crashed_lvecs_a[isys]; ffl.setLvec(lvec);
        //Mat3d lvec=lvec0; lvec.a=crashed_lvecs_a[16]; ffl.setLvec(lvec);
        //Mat3d lvec=lvec0; lvec.a=crashed_lvecs_a[3]; ffl.setLvec(lvec);
        //Mat3d lvec=lvec0; lvec.a.addRandomCube(0.5); ffl.setLvec(lvec);
        //Mat3d lvec=lvec0; lvec.a.addRandomCube(0.4); ffl.setLvec(lvec);
        //Mat3d lvec=lvec0; lvec.a.addRandomCube(0.25); ffl.setLvec(lvec);
        //Mat3d lvec=lvec0; lvec.b.x += randf(-0.5,0.5); ffl.setLvec(lvec);
        //Mat3d lvec=lvec0; lvec.b.x += isys*0.00163; ffl.setLvec(lvec);
        //Mat3d lvec=lvec0; lvec.b.x += isys*0.15; ffl.setLvec(lvec);
        //Mat3d lvec=lvec0; lvec.a.y += isys*2.25; ffl.setLvec(lvec);
        //printf( "   Vec3d{%g,%g,%g}, \n",  lvec.a.x, lvec.a.y, lvec.a.z  );

        copy    ( ff.natoms, ff.neighCell, neighCell+i0a );
        copy    ( ff.natoms, ff.neighs,    neighs   +i0a );
        copy    ( ff.natoms*EXCL_MAX, ff.excl, excl + i0a*EXCL_MAX );
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

void pack_systems( bool bParams=0, bool bForces=false, bool bVel=false, bool blvec=true, float l_rnd=-1 ){
    for(int isys=0; isys<nSystems; isys++){
        pack_system( isys, ffls[isys], bParams, bForces, bVel, blvec, l_rnd );
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
        err|= ocl.upload( ocl.ibuff_excl,      excl      , ocl.nAtoms*EXCL_MAX, i0a*EXCL_MAX  );
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

void upload_mmff(  bool bParams, bool bForces, bool bVel, bool blvec ){
    printf("MolWorld_sp3_multi::upload() nSys %i nAtoms %i nvecs %i bParams %i bForces %i bVel %i blvec %i \n", nSystems, ocl.nAtoms, ocl.nvecs, bParams, bForces, bVel, blvec);
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
        err|= ocl.upload( ocl.ibuff_excl,       excl      );
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

void download_mmff_sys( int isys, bool bForces, bool bVel ){
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

void download_mmff( bool bForces, bool bVel ){
    int err=0;
    //printf("MolWorld_sp3_multi::download() \n");
    ocl.download( ocl.ibuff_atoms, atoms );
    if(bForces){ ocl.download( ocl.ibuff_aforces, aforces ); }
    if(bVel   ){ ocl.download( ocl.ibuff_avel,    avel    ); }
    err |= ocl.finishRaw();
    OCL_checkError(err, "MolWorld_sp2_multi::download().finish");
}

// ==================================
//         GPU <-> CPU I/O  UFF
// ==================================


void pack_uff_system( int isys, const UFF& ff, bool bParams=false, bool bForces=false, bool bVel=false, bool blvec=true ){
    printf("MolWorld_sp3_multi::pack_uff_system(isys=%i) nAtoms=%i nsystems=%i uff.nAtoms=%i\n", isys, ff.natoms, nSystems, ff.natoms);
    int nAtoms      = ff.natoms;
    int i0a = isys * nAtoms;

    // Pack atoms
    for(int i=0; i<nAtoms; ++i){
        atoms[i0a + i] = (Quat4f){(float)ff.apos[i].x, (float)ff.apos[i].y, (float)ff.apos[i].z, 0.f};
    }
    if(bForces){ pack( nAtoms, ff.fapos, aforces+i0a ); }
    if(bVel)   { pack( nAtoms, ff.vapos, avel+i0a );    }

    if(bParams){
        int nBonds      = ff.nbonds;
        int nAngles     = ff.nangles;
        int nDihedrals  = ff.ndihedrals;
        int nInversions = ff.ninversions;
        int nA2F        = ff.nf;

        int i0b = isys * nBonds;
        int i0A = isys * nAngles;
        int i0D = isys * nDihedrals;
        int i0I = isys * nInversions;
        int i0F = isys * nA2F;

        // Pack REQs
        for(int i=0; i<nAtoms; ++i){
            REQs [i0a + i] = (Quat4f){(float)ff.REQs[i].x, (float)ff.REQs[i].y, (float)ff.REQs[i].z, (float)ff.REQs[i].w};
        }

        if(isys == 0 || isys == 1){ // DEBUG: Print REQs for first two systems
            printf("MolWorld_sp3_multi::pack_uff_system(isys=%d).REQs:\n", isys);
            for(int i=0; i<nAtoms; ++i){ Quat4f req = REQs[i0a + i];  printf("REQ [s=%d,a=%d]: %g %g %g %g\n", isys, i, req.x, req.y, req.z, req.w); }
        }

        // Pack per-atom topology used by UFF kernels
        for(int i=0; i<nAtoms; ++i){ host_neighs_UFF   [i0a + i] = (int4)ff.neighs[i]; }
        for(int i=0; i<nAtoms; ++i){ host_neighCell_UFF[i0a + i] = (int4)ff.neighCell[i]; }
        for(int i=0; i<nAtoms; ++i){
            int4 nb = (int4)ff.neighBs[i];
            int* pin = (int*)&nb;
            for(int k=0;k<4;k++){ if(pin[k]>=0) pin[k] += i0b; }
            host_neighBs_UFF[i0a + i] = nb;
        }

        // Pack bonds, angles, dihedrals, inversions
        for(int i=0; i<nBonds;      ++i){ host_bon_atoms [i0b + i] = (int2){ff.bonAtoms[i].x, ff.bonAtoms[i].y}; host_bon_params[i0b + i] = (float2){(float)ff.bonParams[i].x, (float)ff.bonParams[i].y}; }
        int i0h = isys * (nAtoms*4);
        for(int i=0; i<nAngles;     ++i){ host_ang_atoms [i0A + i] = (int4){ff.angAtoms[i].x, ff.angAtoms[i].y, ff.angAtoms[i].z, 0}; host_ang_ngs[i0A + i] = (int2){ff.angNgs[i].x + i0h, ff.angNgs[i].y + i0h}; host_ang_params1[i0A + i] = (float4){(float)ff.angParams[i].c0, (float)ff.angParams[i].c1, (float)ff.angParams[i].c2, (float)ff.angParams[i].c3}; host_ang_params2_w[i0A + i] = (float)ff.angParams[i].k; }
        for(int i=0; i<nDihedrals;  ++i){ host_dih_atoms [i0D + i] = (int4)ff.dihAtoms[i]; host_dih_ngs[i0D + i] = (int4){ff.dihNgs[i].x + i0h, ff.dihNgs[i].y + i0h, ff.dihNgs[i].z + i0h, 0}; host_dih_params[i0D + i] = (float4){(float)ff.dihParams[i].x, (float)ff.dihParams[i].y, (float)ff.dihParams[i].z, 0.f}; }
        for(int i=0; i<nInversions; ++i){ host_inv_atoms [i0I + i] = (int4)ff.invAtoms[i]; host_inv_ngs[i0I + i] = (int4){ff.invNgs[i].x + i0h, ff.invNgs[i].y + i0h, ff.invNgs[i].z + i0h, 0}; host_inv_params[i0I + i] = (float4){(float)ff.invParams[i].x, (float)ff.invParams[i].y, (float)ff.invParams[i].z, (float)ff.invParams[i].w}; }

        // Pack force assembly map
        for(int i=0; i<nAtoms; ++i){
            host_a2f_offsets[i0a + i] = ff.a2f.cellI0s[i];
            host_a2f_counts [i0a + i] = ff.a2f.cellNs[i];
        }
        for(int i=0; i<nA2F; ++i){
            host_a2f_indices[i0F + i] = ff.a2f.cell2obj[i];
        }
    }

    if(blvec){
        // Pack PBC
        Mat3_to_cl( ff.lvec,    lvecs[isys] );
        Mat3_to_cl( ff.invLvec, ilvecs[isys] );
        int i0pbc = isys * (npbc+1);
        pack( npbc, ff.shifts, pbcshifts + i0pbc );
    }
}

void unpack_uff_system( int isys, UFF& ff, bool bForces=false, bool bVel=false ){
    int i0a = isys * ff.natoms;
    unpack( ff.natoms, ff.apos, atoms+i0a );
    if(bForces){ unpack( ff.natoms, ff.fapos, aforces+i0a ); }
    if(bVel)   { unpack( ff.natoms, ff.vapos, avel+i0a );    }
}


void upload_uff_sys( int isys, bool bParams, bool bForces, bool bVel, bool blvec, bool bPrint=false ){
    //bPrint=true;
    if(bPrint)printf("MolWorld_sp3_multi::upload_uff_sys(%i) nAtoms=%i nBonds=%i nAngles=%i nDihedrals=%i nInversions=%i \n", isys, uff_ocl->nAtoms, uff_ocl->nBonds, uff_ocl->nAngles, uff_ocl->nDihedrals, uff_ocl->nInversions);
    int i0a = isys * uff_ocl->nAtoms;
    int err=0;
    err             |= uff_ocl->upload( uff_ocl->ibuff_apos,  (float*)(atoms  +i0a), uff_ocl->nAtoms );
    if(bForces){ err|= uff_ocl->upload( uff_ocl->ibuff_fapos, (float*)(aforces+i0a), uff_ocl->nAtoms ); }
    if(blvec){
        err|= uff_ocl->upload( uff_ocl->ibuff_lvecs, lvecs+isys, 1 );
        // Upload at least one PBC shift record even if npbc==0 (buffer was allocated with safe size 1)
        int nPBC_safe = (npbc>0)? npbc : 1;
        int i0pbc = isys * nPBC_safe;
        err|= uff_ocl->upload( uff_ocl->ibuff_pbcshifts, pbcshifts + i0pbc, nPBC_safe );
    }
    if(bParams){
        int i0b = isys * uff_ocl->nBonds;
        int i0A = isys * uff_ocl->nAngles;
        int i0D = isys * uff_ocl->nDihedrals;
        int i0I = isys * uff_ocl->nInversions;
        int i0F = isys * uff_ocl->nA2F;
        // Topology needed by all kernels: neighbors and neighbor bond indices
        err|= uff_ocl->upload( uff_ocl->ibuff_neighs,           (int*)  (ffu.neighs               + i0a), uff_ocl->nAtoms,      i0a, bPrint );
        err|= uff_ocl->upload( uff_ocl->ibuff_neighCell,        (int*)  (ffu.neighCell            + i0a), uff_ocl->nAtoms,      i0a, bPrint );
        err|= uff_ocl->upload( uff_ocl->ibuff_neighBs,          (int*)  (ffu.neighBs              + i0a), uff_ocl->nAtoms,      i0a, bPrint );
        if(uff_ocl->nBonds>0){
            err|= uff_ocl->upload( uff_ocl->ibuff_bonAtoms,     (int*)  (host_bon_atoms.data()    + i0b), uff_ocl->nBonds,      i0b, bPrint );
            err|= uff_ocl->upload( uff_ocl->ibuff_bonParams,    (float*)(host_bon_params.data()   + i0b), uff_ocl->nBonds,      i0b, bPrint );
        }
        if(uff_ocl->nAngles>0){
            err|= uff_ocl->upload( uff_ocl->ibuff_angAtoms,     (int*)  (host_ang_atoms.data()    + i0A), uff_ocl->nAngles,     i0A, bPrint );
            err|= uff_ocl->upload( uff_ocl->ibuff_angNgs,       (int*)  (host_ang_ngs.data()      + i0A), uff_ocl->nAngles,     i0A, bPrint );
            err|= uff_ocl->upload( uff_ocl->ibuff_angParams1,   (float*)(host_ang_params1.data()  + i0A), uff_ocl->nAngles,     i0A, bPrint );
            err|= uff_ocl->upload( uff_ocl->ibuff_angParams2_w, (float*)(host_ang_params2_w.data()+ i0A), uff_ocl->nAngles,     i0A, bPrint );
        }
        if(uff_ocl->nDihedrals>0){
            err|= uff_ocl->upload( uff_ocl->ibuff_dihAtoms,      (int*)  (host_dih_atoms.data()   + i0D), uff_ocl->nDihedrals,  i0D, bPrint );
            err|= uff_ocl->upload( uff_ocl->ibuff_dihNgs,        (int*)  (host_dih_ngs.data()     + i0D), uff_ocl->nDihedrals,  i0D, bPrint );
            err|= uff_ocl->upload( uff_ocl->ibuff_dihParams,     (float*)(host_dih_params.data()  + i0D), uff_ocl->nDihedrals,  i0D, bPrint );
        }
        if(uff_ocl->nInversions>0){
            err|= uff_ocl->upload( uff_ocl->ibuff_invAtoms,      (int*)  (host_inv_atoms.data()   + i0I), uff_ocl->nInversions, i0I, bPrint );
            err|= uff_ocl->upload( uff_ocl->ibuff_invNgs,        (int*)  (host_inv_ngs.data()     + i0I), uff_ocl->nInversions, i0I, bPrint );
            err|= uff_ocl->upload( uff_ocl->ibuff_invParams,     (float*)(host_inv_params.data()  + i0I), uff_ocl->nInversions, i0I, bPrint );
        }
        err|= uff_ocl->upload( uff_ocl->ibuff_a2f_offsets,       host_a2f_offsets.data()          + i0a,  uff_ocl->nAtoms,      i0a, bPrint );
        err|= uff_ocl->upload( uff_ocl->ibuff_a2f_counts,        host_a2f_counts.data()           + i0a,  uff_ocl->nAtoms  ,    i0a, bPrint );
        err|= uff_ocl->upload( uff_ocl->ibuff_a2f_indices,       host_a2f_indices.data()          + i0F,  uff_ocl->nA2F    ,    i0F, bPrint );
    }
    err |= uff_ocl->finishRaw();
    OCL_checkError(err, "MolWorld_sp3_multi::upload_uff_sys().finish");
}

void upload_uff( bool bParams, bool bForces, bool bVel, bool blvec, bool bPrint=false ){
    //bPrint=true;
    if(bPrint)printf("MolWorld_sp3_multi::upload_uff() nSystems=%i nAtomsTot=%i nBondsTot=%i nAnglesTot=%i nDihedralsTot=%i nInversionsTot=%i nA2FTot=%i \n",  uff_ocl->nSystems, uff_ocl->nAtoms, uff_ocl->nBonds, uff_ocl->nAngles, uff_ocl->nDihedrals, uff_ocl->nInversions, uff_ocl->nA2F);
    int err=0;
    int nAtomsTot      = uff_ocl->nSystems * uff_ocl->nAtoms;
    int nBondsTot      = uff_ocl->nSystems * uff_ocl->nBonds;
    int nAnglesTot     = uff_ocl->nSystems * uff_ocl->nAngles;
    int nDihedralsTot  = uff_ocl->nSystems * uff_ocl->nDihedrals;
    int nInversionsTot = uff_ocl->nSystems * uff_ocl->nInversions;
    int nA2FTot        = uff_ocl->nSystems * uff_ocl->nA2F;
    if(bPrint)printf("MolWorld_sp3_multi::upload_uff() nAtomsTot=%i nBondsTot=%i nAnglesTot=%i nDihedralsTot=%i nInversionsTot=%i nA2FTot=%i \n", nAtomsTot, nBondsTot, nAnglesTot, nDihedralsTot, nInversionsTot, nA2FTot);
    err|=             uff_ocl->upload( uff_ocl->ibuff_apos,  (float*)atoms, nAtomsTot );
    err|=             uff_ocl->upload( uff_ocl->ibuff_REQs,  (float*)REQs,  nAtomsTot );
    if(bForces) err|= uff_ocl->upload( uff_ocl->ibuff_fapos, (float*)aforces );
    if(blvec){
        err|= uff_ocl->upload( uff_ocl->ibuff_lvecs,     lvecs );
        // Upload at least one PBC shift record even if npbc==0 (buffer allocated with safe size 1)
        int nPBC_safe = (npbc>0)? npbc : 1;
        err|= uff_ocl->upload( uff_ocl->ibuff_pbcshifts, (float*)pbcshifts, nSystems * nPBC_safe );
    }
    OCL_checkError(err, "MolWorld_sp3_multi::upload_uff().1");
    if(bParams){
        // Topology needed by all kernels (bulk upload across all systems) from packed host buffers
        err|= uff_ocl->upload( uff_ocl->ibuff_neighs,            (int*)host_neighs_UFF.data(),     nAtomsTot,    0, bPrint );
        err|= uff_ocl->upload( uff_ocl->ibuff_neighCell,         (int*)host_neighCell_UFF.data(),  nAtomsTot,    0, bPrint );
        err|= uff_ocl->upload( uff_ocl->ibuff_neighBs,           (int*)host_neighBs_UFF.data(),    nAtomsTot,    0, bPrint );
        OCL_checkError(err, "MolWorld_sp3_multi::upload_uff().2");
        if(nBondsTot>0){
            err|= uff_ocl->upload( uff_ocl->ibuff_bonAtoms,      (int*)host_bon_atoms.data(),   nBondsTot,    0, bPrint );
            err|= uff_ocl->upload( uff_ocl->ibuff_bonParams,     (float*)host_bon_params.data(),nBondsTot,    0, bPrint );
            OCL_checkError(err, "MolWorld_sp3_multi::upload_uff().3");
        }
        if(nAnglesTot>0){
            err|= uff_ocl->upload( uff_ocl->ibuff_angAtoms,     (int*)host_ang_atoms.data(),   nAnglesTot,    0, bPrint );
            err|= uff_ocl->upload( uff_ocl->ibuff_angNgs,       (int*)host_ang_ngs.data(),     nAnglesTot,    0, bPrint );
            err|= uff_ocl->upload( uff_ocl->ibuff_angParams1,   host_ang_params1.data(),       nAnglesTot,    0, bPrint );
            err|= uff_ocl->upload( uff_ocl->ibuff_angParams2_w, host_ang_params2_w.data(),     nAnglesTot,    0, bPrint );
            OCL_checkError(err, "MolWorld_sp3_multi::upload_uff().4");
        }
        if(nDihedralsTot>0){
            err|= uff_ocl->upload( uff_ocl->ibuff_dihAtoms,     (int*)host_dih_atoms.data(),   nDihedralsTot, 0, bPrint );
            err|= uff_ocl->upload( uff_ocl->ibuff_dihNgs,       (int*)host_dih_ngs.data(),     nDihedralsTot, 0, bPrint );
            err|= uff_ocl->upload( uff_ocl->ibuff_dihParams,    (float*)host_dih_params.data(),nDihedralsTot, 0, bPrint );
            OCL_checkError(err, "MolWorld_sp3_multi::upload_uff().5");
        }
        if(nInversionsTot>0){
            err|= uff_ocl->upload( uff_ocl->ibuff_invAtoms,     (int*)host_inv_atoms.data(),   nInversionsTot, 0, bPrint );
            err|= uff_ocl->upload( uff_ocl->ibuff_invNgs,       (int*)host_inv_ngs.data(),     nInversionsTot, 0, bPrint );
            err|= uff_ocl->upload( uff_ocl->ibuff_invParams,    (float*)host_inv_params.data(),nInversionsTot, 0, bPrint );
            OCL_checkError(err, "MolWorld_sp3_multi::upload_uff().6");
        }
    }
    err|= uff_ocl->upload( uff_ocl->ibuff_a2f_offsets, host_a2f_offsets.data(), nAtomsTot, 0, bPrint );
    err|= uff_ocl->upload( uff_ocl->ibuff_a2f_counts,  host_a2f_counts.data(),  nAtomsTot, 0, bPrint );
    err|= uff_ocl->upload( uff_ocl->ibuff_a2f_indices, host_a2f_indices.data(), nA2FTot,   0, bPrint );
    err |= uff_ocl->finishRaw(); OCL_checkError(err, "MolWorld_sp3_multi::upload_uff().finish");
}

void download_uff_sys( int isys, bool bForces, bool bVel ){
    int i0a = isys * uff_ocl->nAtoms;
    uff_ocl->download( uff_ocl->ibuff_apos, (float*)(atoms+i0a), uff_ocl->nAtoms );
    if(bForces){ uff_ocl->download( uff_ocl->ibuff_fapos, (float*)(aforces+i0a), uff_ocl->nAtoms ); }
    //if(bVel   ){ uff_ocl->download( uff_ocl->ibuff_avel,  (float*)(avel+i0a),    uff_ocl->nAtoms ); } // not implemented yet
    uff_ocl->finishRaw();
}


void download_uff( bool bForces, bool bVel ){
    uff_ocl->download_results( (float*)aforces, 0 ); // energies not handled yet
    uff_ocl->download( uff_ocl->ibuff_apos, (float*)atoms );
    uff_ocl->finishRaw();
}

// ==================================
//         GPU <-> CPU I/O  UFF
// ==================================

void upload(  bool bParams=true, bool bForces=false, bool bVel=false, bool blvec=true ){
    if(bUFF) upload_uff( bParams, bForces, bVel, blvec ); else upload_mmff( bParams, bForces, bVel, blvec );
}

void download( bool bForces=false, bool bVel=false ){
    if(bUFF) download_uff( bForces, bVel ); else download_mmff( bForces, bVel );
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
    //MDpar.x = fire.dt * 0.5;
    //MDpar.x = 0.01;
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

    int nMaxSteps = 10000/nPerVFs;

    for(int isys=0; isys<nSystems; isys++){
        int i0v = isys * ocl.nvecs;
        double f2 = fire[isys].ff;
        // -------- Global Optimization
        if( ( f2 < F2conv ) && (!gopts[isys].bExploring) ){            // Start Exploring
            nbConverged++;
            nStepConvSum+=gopts[isys].istep*nPerVFs;
            gopts[isys].startExploring();
            if(bGroups){
                bGroupUpdate=true;
                //printf("MolWorld_sp3_multi::evalVFs() isys=%3i Start Exploring \n", isys );
                //for(int ig=0; ig<ocl.nGroup; ig++){ setGroupDrive(isys, ig, {fsc,fsc,0.0} ); }   // Shift driver
                //for(int ig=0; ig<ocl.nGroup; ig++){ setGroupDrive(isys, ig, Vec3fZero, {tsc,0.0,0.0} ); }   // group rotate driver
                for(int ig=0; ig<ocl.nGroup; ig++){
                    Vec2f fDrive = groups.groups[ig].fDrive;
                    setGroupDrive(isys, ig, {fDrive.x,fDrive.x,0.0}, {fDrive.y,0.0,0.0} );
                    //setGroupDrive(isys, ig, {fsc,fsc,0.0}, {tsc,0.0,0.0} );
                }   // Shift driver
            }
        }
        else if( ( gopts[isys].istep > nMaxSteps ) && (!gopts[isys].bExploring)){
            nbNonConverged++;
            nStepNonConvSum+=gopts[isys].istep*nPerVFs;
            gopts[isys].startExploring();
            if(bGroups){ bGroupUpdate=true; }
            database->convergedStructure.push_back(false);

            int i0v = isys * ocl.nvecs;
            unpack( ffls[isys].nvecs,  ffls[isys].apos, atoms+i0v);
            if(database->addIfNewDescriptor(&ffls[isys])==-1){
                sprintf(tmpstr,"# %i E %g |F| %g istep=%i", database->getNMembers(), ffls[isys].Etot, sqrt(ffl.cvf.z), gopts[isys].istep );
                saveXYZ( "gopt.xyz", tmpstr, false, "a", nPBC_save );
                database->convergedStructure.push_back(false);
            }

            for(int ig=0; ig<ocl.nGroup; ig++){ setGroupDrive(isys, ig, {fsc,fsc,0.0}, {tsc,0.0,0.0} ); }   // Shift driver
        }
        if( gopts[isys].update() ){ // Stop Exploring
            bGroupUpdate=true;
            nExploring++;
            nStepExplorSum+=gopts[isys].nExplore*nPerVFs;
            //printf("MolWorld_sp3_multi::evalVFs() isys=%3i Stop Exploring \n", isys );
            for(int ig=0; ig<ocl.nGroup; ig++){ setGroupDrive(isys, ig, Vec3fZero, Vec3fZero ); }
        };
        bExploring |= gopts[isys].bExploring;
        //printf("gopts[%i].bExploring=%i\n", isys, gopts[isys].bExploring );
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
        nbEvaluation++;
        int i0v = isys * ocl.nvecs;
        //evalVF( ocl.nvecs, aforces+i0v, avel   +i0v, fire[isys], MDpars[isys] );
        evalVF_new( ocl.nvecs, cvfs+i0v, fire[isys], MDpars[isys], gopts[isys].bExploring );
        double f2 = fire[isys].ff;
        if(f2>F2max){ F2max=f2; iSysFMax=isys; }
        // -------- Global Optimization
        if( ( f2 < F2conv ) && (!gopts[isys].bExploring) ){            // Start Exploring
            int i0v = isys * ocl.nvecs;
            unpack( ffls[isys].nvecs,  ffls[isys].apos, atoms+i0v);

            if(bMILAN){
                int sameMember = database->addIfNewDescriptor(&ffls[isys]);
                if(sameMember==-1){
                    std::vector<double> theta;
                    for(auto g : groups.groups){
                        theta.push_back( g.fw.dot(Vec3dZ)/g.fw.norm()/2/M_PI*360 );
                    }
                    sprintf(tmpstr,"# %i E %g |F| %g istep=%i isys=%i,", database->getNMembers(), 0.5*fire[isys].vv, sqrt(f2), gopts[isys].istep, isys );
                    nbmol.copyOf( ffls[isys] );
                    Vec3d cog = nbmol.getBBcog();
                    for(int i=0; i<nbmol.natoms; i++){
                        nbmol.apos[i].sub( cog );
                    }
                    saveXYZ( "gopt.xyz", tmpstr, false, "a", nPBC_save );
                    database->convergedStructure.push_back(true);
                }
            }

            // //printf( "evalVFs() iSys=%i CONVERGED |F|=%g \n", isys, sqrt(f2) );
            // gopts[isys].startExploring();
            // bGroupUpdate=true;
            // //printf("MolWorld_sp3_multi::evalVFs() isys=%3i Start Exploring \n", isys );
            // //for(int ig=0; ig<ocl.nGroup; ig++){ setGroupDrive(isys, ig, {fsc,fsc,0.0} ); }   // Shift driver
            // for(int ig=0; ig<ocl.nGroup; ig++){ setGroupDrive(isys, ig, Vec3fZero, {tsc,0.0,0.0} ); }   // group rotate driver
        }
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

void checkBordersOfBbox(){
    int err=0;
    Vec3d period = {bbox.b.x-bbox.a.x, bbox.b.y-bbox.a.y, bbox.b.z-bbox.a.z};
    for(int isys=0; isys<nSystems; isys++){
        int i0v = isys * ocl.nvecs;
        for(int i=0; i<ocl.nvecs; i++){
            Vec3f p = atoms[i0v+i].f;
            if( p.x < bbox.a.x - period.a){ for(int j=0; j<ocl.nvecs; j++)atoms[i0v+j].f.x += 2*period.a; break; }
            if( p.x > bbox.b.x + period.a){ for(int j=0; j<ocl.nvecs; j++)atoms[i0v+j].f.x -= 2*period.a; break; }
            if( p.y < bbox.a.y - period.b){ for(int j=0; j<ocl.nvecs; j++)atoms[i0v+j].f.y += 2*period.b; break; }
            if( p.y > bbox.b.y + period.b){ for(int j=0; j<ocl.nvecs; j++)atoms[i0v+j].f.y -= 2*period.b; break; }
            if( p.z < bbox.a.z - period.c){ for(int j=0; j<ocl.nvecs; j++)atoms[i0v+j].f.z += 2*period.c; break; }
            if( p.z > bbox.b.z + period.c){ for(int j=0; j<ocl.nvecs; j++)atoms[i0v+j].f.z -= 2*period.c; break; }
            //atoms[i0v+i].f = p;
        }
    }
    err |= ocl.upload( ocl.ibuff_atoms, atoms );
    err |= ocl.finishRaw();  OCL_checkError(err, "checkBorders().2");

}

void setGroupDrive(int isys, int ig, Vec3f fsc=Vec3fZero, Vec3f tsc=Vec3fZero){
    int igs = ig + isys*ocl.nGroup;
    gforces[igs] = Quat4f{randf(-fsc.x,fsc.x),randf(-fsc.y,fsc.y),randf(-fsc.z,fsc.z),0.0};
    gtorqs [igs] = Quat4f{randf(-tsc.x,tsc.x),randf(-tsc.y,tsc.y),randf(-tsc.z,tsc.z),0.0};
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

void printConstrains(){
    printf( "MolWorld_sp3_multi::printConstrains()\n" );
    for(int isys=0; isys<nSystems; isys++){
        int i0a   = isys * ocl.nAtoms;
        int i0v   = isys * ocl.nvecs;
        for(int ia=0; ia<ffls[isys].natoms; ia++){
            if(ffls[isys].constr[ia].w>0)printf( "ffls[isys=%i,ia=%i] constr(%g,%g,%g|%g) constrK(%g,%g,%g) \n", isys, ia, ffls[isys].constr[ia].x,ffls[isys].constr[ia].y,ffls[isys].constr[ia].z,ffls[isys].constr[ia].w, ffls[isys].constrK[ia].x,ffls[isys].constrK[ia].y,ffls[isys].constrK[ia].z );
        }
    }
}

void move_MultiConstrain( Vec3d d, Vec3d di, float Kfixmin=0.001f ){
    printf( "MolWorld_sp3_multi::move_MultiConstrain()\n" );
    for(int isys=0; isys<nSystems; isys++){
        int i0a   = isys * ocl.nAtoms;
        int i0v   = isys * ocl.nvecs;
        for(int ia=0; ia<ffls[isys].natoms; ia++){
            if( ffls[isys].constr[ia].w > Kfixmin ){
                int iaa = ia+i0a;
                //ffls[isys].constr[ia].w=1.0f;
                //ffls[isys].constr[ia].f.add( d + di*isys );   // Shift Constrains by di per isys
                ffls[isys].constr[ia].f.add( d );   // Shift Constrains by di per isys
                constr [iaa]=(Quat4f)ffls[isys].constr[ia];
                constrK[iaa]=(Quat4f)ffls[isys].constrK[ia];
                //printf( "ffls[isys=%i,ia=%i] constr(%g,%g,%g|%g) constrK(%g,%g,%g) \n", isys, ia, ffls[isys].constr[ia].x,ffls[isys].constr[ia].y,ffls[isys].constr[ia].z,ffls[isys].constr[ia].w, ffls[isys].constrK[ia].x,ffls[isys].constrK[ia].y,ffls[isys].constrK[ia].z );
            };
        };
    }
    //for(int isys=0; isys<nSystems; isys++){ printf("replica %i : ", isys); ffls[isys].printAtomsConstrains(); }
    //upload( ocl.ibuff_constr, constr );
}

virtual void setConstrains(bool bClear=true, double Kfix=1.0 ){
    printf("MolWorld_sp3_multi::setConstrains()\n");
    //printConstrains();
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
    //printConstrains();
    move_MultiConstrain( Vec3d{0.0, 0.0, 5.0}, Vec3d{0.0, 0.4, 0.0} );
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

    ocl.upload( ocl.ibuff_constr, constr );
    //exit(0);
}





virtual int paralel_size( )override{ return nSystems; }

double eval_UFF_ocl( int niter=1 ){
    printf("MolWorld_sp3_multi::eval_UFF_ocl()\n");
    if(!bUFF) return 0.0;
    if(!uff_ocl->bKernelPrepared) uff_ocl->setup_kernels( ffu.Rdamp, ffu.FmaxNonBonded, ffu.SubNBTorsionFactor ); // setup kernels if not done yet
    for(int i=0; i<niter; i++){
        uff_ocl->eval();
    }
    uff_ocl->finishRaw();
    std::vector<float> energies( nSystems * 5 );
    uff_ocl->download( uff_ocl->ibuff_energies, energies.data() );
    uff_ocl->finishRaw();
    return energies[4]; // E_tot for system 0
}

void setup_UFF_ocl(){
    printf("MolWorld_sp3_multi::setup_UFF_ocl()\n");
    if(!uff_ocl){ printf("ERROR: setup_UFF_ocl() called but uff_ocl is null!\n"); return; }

    // Ensure kernels are prepared if they haven't been already.
    if(!uff_ocl->bKernelPrepared){
        float Rdamp        = (float)ffu.Rdamp;
        float FmaxNB       = (float)ffu.FmaxNonBonded;
        float SubNBTorsion = (float)ffu.SubNBTorsionFactor;
        uff_ocl->setup_kernels(Rdamp, FmaxNB, SubNBTorsion);
    }

    Vec3i nPBC_ = nPBC; if(!bPBC) nPBC_ = Vec3iZero;

    if(bGridFF){
        if(!uff_ocl->task_NBFF_Grid_Bspline){
            printf("MolWorld_sp3_multi::setup_UFF_ocl() for UFF GridFF\n");
            // UFF should set grid parameters directly from GridFF object, not copy from MMFF context
            // Set grid parameters directly from GridFF object
            uff_ocl->setGridShape(gridFF.grid);
            uff_ocl->bNonBond=bNonBonded;
            printf("DEBUG setup_UFF_ocl: uff_ocl->ibuff_BsplinePLQ=%d (should be 41)\n", uff_ocl->ibuff_BsplinePLQ);
            uff_ocl->task_NBFF_Grid_Bspline = uff_ocl->setup_getNonBond_GridFF_Bspline(ffu.natoms, nPBC_);
        }
    } else if(bNonBonded){
        if(!uff_ocl->task_NBFF){
            printf("MolWorld_sp3_multi::setup_UFF_ocl() for UFF NonBonded\n");
            uff_ocl->task_NBFF = uff_ocl->setup_getNonBond(ffu.natoms, nPBC_);
        }
    }
    if(!uff_ocl->task_updateAtoms) uff_ocl->setup_updateAtomsMMFFf4( ffu.natoms, 0 );
}

virtual double solve_multi ( int nmax, double tol )override{
    //return eval_MMFFf4_ocl( nmax, tol );
    //return run_ocl_opt( nmax, tol );
    if(bUFF){ return eval_UFF_ocl(nmax); }
    int nitrdione = 0;
    switch(iParalel){
        case -1: bOcl=false;  nitrdione = run_multi_serial(nmax,tol);  break;
        case  0: //break;
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

int saveXYZ_system(int isys, const char* fname, const char* comment="#comment", const char* mode="w"){
    printf("MolWorld_sp3_multi::saveXYZ_system(isys=%i, fname=%s, mode=%s, comment=%s)\n", isys, fname, mode, comment);
    if(isys < 0 || isys >= nSystems) return -1;
    if(bUFF){
        unpack_uff_system(isys, ffu, true, false);
        return params.saveXYZ( fname, ffu.natoms, ffu.atypes, ffu.apos, comment, ffu.REQs, mode, true );
    }else{
        unpack_system(isys, ffls[isys], true, false);
        return params.saveXYZ( fname, ffls[isys].natoms, ffls[isys].atypes, ffls[isys].apos, comment, ffls[isys].REQs, mode, true );
    }
    return 0;
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
        sprintf( fname  , "debug_sys%03i.xyz", isys );
        sprintf( commnet, "# itr %i |F| %g ", itr, F );
        unpack_system( isys, ffl, false, false );
        if( ckeckNaN_d( ffl.nvecs, 3, (double*)ffl.apos, fname, false ) ){ printf( "saveDebugXYZreplicas(itr=%i) NaNs in replica[%i] => Exit() \n", itr, isys  ); exit(0); };
        saveXYZ( fname, commnet, false, "a" );
    }
}


void findAFMGrid(){


}

virtual void evalAFMscan( GridShape& scan, Quat4f*& OutFE, Quat4f*& OutPos, Quat4f** ps=0, bool bSaveDebug=false )override{

    printf( "MolWrold_sp3_multi::evalAFMscan() scan.ns(%i,%i,%i) \n", scan.n.x,scan.n.y,scan.n.z );
    int err=0;
    //Quat4f* OutFE =0;
    //Quat4f* OutPos=0;
    //if( OutFE_  ){ OutFE =OutFE_;  }else{ OutFE =new Quat4f[scan.n.totprod()];   }
    //if( OutPos_ ){ OutPos=OutPos_; }else{ OutPos=new Quat4f[scan.n.totprod()]; }

    // ---- update scan bounds
    // ToDo: ztop should be probably found and used in  evalAFM_FF()
    int i0v = iSystemCur*ocl.nvecs;
    float ztop= -1e+300;
    Vec2f pmin=Vec2fmax,pmax=Vec2fmin;
    for( int i=0; i<ocl.nAtoms; i++ ){
        Quat4f p = atoms[i0v+i];
        p.xy().update_bounds(pmin,pmax);
        ztop=fmax( atoms[i0v+i].z, ztop );
    };
    float margin=3.0;
    pmin.x-=margin,pmin.y-=margin, pmax.x+=margin, pmax.y+=margin;
    printf( "MolWrold_sp3_multi::evalAFMscan() ztop=%g pmin(%g,%g) pmax(%g,%g)\n", ztop, pmin.x,pmin.y, pmax.x, pmax.y );
    scan.pos0=Vec3d{pmin.x,pmin.y,ztop+4.0+5.0+1.0};
    scan.cell.a=Vec3d{pmax.x-pmin.x,0.0,0.0};
    scan.cell.b=Vec3d{pmax.x-pmin.x,0.0,0.0};
    scan.autoN( scan.dCell.xx );
    scan.printCell();

    int nxy = scan.n.x*scan.n.y;
    if( OutFE ==0 ){ OutFE =new Quat4f[scan.n.totprod()]; }
    if( OutPos==0 ){ OutPos=new Quat4f[scan.n.totprod()]; }
    _realloc( afm_ps, nxy );

    for(int iy=0; iy<scan.n.y; iy++) for(int ix=0; ix<scan.n.x; ix++){ afm_ps[iy*scan.n.x+ix].f = (Vec3f)(  scan.pos0 + scan.dCell.a*ix + scan.dCell.b*iy ); afm_ps[iy*scan.n.x+ix].w=0; }
    printf( "MolWrold_sp3_multi::evalAFMscan() scan.ns(%i,%i,%i) nxy=%i \n", scan.n.x,scan.n.y,scan.n.z,nxy );
    long T0=getCPUticks();
    ocl.PPAFM_scan( nxy, scan.n.z, afm_ps, OutFE, OutPos, scan.dCell.c.z );
    printf( ">>time(surf2ocl.download() %g \n", (getCPUticks()-T0)*tick2second );

    // ToDo: we should probably return vmin,vmax to MolGUI for visualization?
    int ntot = scan.n.x*scan.n.y*scan.n.z;
    Quat4f vmin=Quat4fmax,vmax=Quat4fmin;
    for(int i=0; i<ntot; i++){ OutFE[i].update_bounds( vmin, vmax ); }
    printf( "MolWrold_sp3_multi::evalAFMscan() vmin{%g,%g,%g|%g} vmax{%g,%g,%g|%g} \n", vmin.x,vmin.y,vmin.z,vmin.w,  vmax.x,vmax.y,vmax.z,vmax.w  );
    bSaveDebug=true;
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
    printf("MolWorld_sp3_multi::setup_MMFFf4_ocl() bUFF=%d\n", bUFF);
    ocl.nDOFs.x=ff.natoms;
    ocl.nDOFs.y=ff.nnode;
    if(!task_cleanF)   task_cleanF = ocl.setup_cleanForceMMFFf4 ( ffl.natoms, ffl.nnode       );
    if(!task_move  )   task_move   = ocl.setup_updateAtomsMMFFf4( ffl.natoms, ffl.nnode       );
    if(!task_print )   task_print  = ocl.setup_printOnGPU       ( ffl.natoms, ffl.nnode       );
    if(!task_MMFF  )   task_MMFF   = ocl.setup_getMMFFf4        ( ffl.natoms, ffl.nnode, bPBC );

    Vec3i nPBC_=nPBC; if(!bPBC){ nPBC_=Vec3iZero; }; printf( "MolWorld_sp3_multi::setup_MMFFf4_ocl() bPBC=%i nPBC(%i,%i,%i) n", bPBC, nPBC_.x,nPBC_.y,nPBC_.z );

    if(bGridFF) switch( gridFF.mode ){
        case GridFFmod::LinearFloat: [[fallthrough]]
        case GridFFmod::LinearDouble:   if(!task_NBFF_Grid         ){
            printf( "MolWorld_sp3_multi::setup_MMFFf4_ocl() GridFFmod::%i => ocl.setup_getNonBond_GridFF() \n", gridFF.mode );
            task_NBFF_Grid         = ocl.setup_getNonBond_GridFF        ( ffl.natoms, ffl.nnode, nPBC_ ); } break;
        case GridFFmod::BsplineFloat:  [[fallthrough]]
        case GridFFmod::BsplineDouble:
            if(ocl.bUseTexture){
                if(!task_NBFF_Grid_Bspline_tex ){
                    printf( "MolWorld_sp3_multi::setup_MMFFf4_ocl() GridFFmod::%i && bUseTexture=%i => ocl.setup_getNonBond_GridFF_Bspline_tex() \n", gridFF.mode, ocl.bUseTexture );
                    task_NBFF_Grid_Bspline_tex = ocl.setup_getNonBond_GridFF_Bspline_tex( ffl.natoms, ffl.nnode, nPBC_ );
                }
            }else{
                if(!task_NBFF_Grid_Bspline ){
                    printf( "MolWorld_sp3_multi::setup_MMFFf4_ocl() GridFFmod::%i && bUseTexture=%i => ocl.setup_getNonBond_GridFF_Bspline() \n", gridFF.mode, ocl.bUseTexture );
                    task_NBFF_Grid_Bspline = ocl.setup_getNonBond_GridFF_Bspline( ffl.natoms, ffl.nnode, nPBC_ );
                }
            }
        break;
        default: { printf( "MolWorld_sp3_multi::setup_MMFFf4_ocl() GridFFmod::%i not implemented \n", gridFF.mode ); } break;
    }

    if(!task_NBFF    ){  task_NBFF      = ocl.setup_getNonBond    ( ffl.natoms, ffl.nnode, nPBC_ );  }
    if(!task_NBFF_ex2){  task_NBFF_ex2  = ocl.setup_getNonBond_ex2( ffl.natoms, ffl.nnode, nPBC_ );  }

    // OCLtask* getSurfMorse(  Vec3i nPBC_, int na=0, float4* atoms=0, float4* REQs=0, int na_s=0, float4* atoms_s=0, float4* REQs_s=0,  bool bRun=true, OCLtask* task=0 ){
    if(!task_SurfAtoms && bSurfAtoms ){
        if( bDONE_surf2ocl==false ){ printf("ERROR in MolWorld_sp3_multi::setup_MMFFf4_ocl() must call MolWorld_sp3_multi::surf2ocl() first\n"); exit(0); }
        task_SurfAtoms = ocl.getSurfMorse( gridFF.nPBC, gridFF.natoms, 0,0,0,0,0,false );
    }
    // ocl.makeGridFF( gridFF.grid, nPBC, gridFF.natoms, (float4*)atoms_surf, (float4*)REQs_surf, true );
    /// HERE_HERE

    if( ocl.nGroupTot > 0 ){
        task_GroupUpdate=ocl.setup_updateGroups( );
        task_GroupForce =ocl.setup_groupForce(   );
    }
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

void picked2GPU( int ipick,  float K ){
    //printf( "MolWorld_sp3_multi::picked2GPU() ipick %i iSystemCur %i \n", ipick, iSystemCur );
    int i0a = ocl.nAtoms*iSystemCur;
    int i0v = ocl.nvecs *iSystemCur;
    if(ipick>=0){
        Quat4f& acon  = constr [i0a + ipick];
        Quat4f& aconK = constrK[i0a + ipick];
        Vec3f hray   = (Vec3f)pick_hray;
        Vec3f ray0   = (Vec3f)pick_ray0;
        const Quat4f& atom = atoms [i0v + ipick];
        float c = hray.dot( atom.f ) - hray.dot( ray0 );
        acon.f = ray0 + hray*c;
        acon.w = K;
        aconK = Quat4f{K,K,K,K};
    }else{
        for(int i=0; i<ocl.nAtoms; i++){ constr[i0a + i].w=-1.0;  };
        for(int i=0; i<ocl.nAtoms; i++){ constrK[i0a + i]=Quat4fZero; };
        //for(int i: constrain_list     ){ constr[i0a + i].w=Kfix;  };
        for(int ia=0; ia<ocl.nAtoms; ia++){
            int iaa = i0a + ia;
            constr [iaa]=(Quat4f)ffls[iSystemCur].constr [ia];
            constrK[iaa]=(Quat4f)ffls[iSystemCur].constrK[ia];
        };
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
            {
                download( true, true );
                err|=ocl.finishRaw();  OCL_checkError(err, "eval_MMFFf4_ocl().debug.download");
                for(int isys=0; isys<nSystems; isys++){
                    unpack_system( isys, ffl, false, false );

                    double frange=100.0;
                    bool berr=false;
                    berr|= ckeckRange( ffl.nvecs, 3, (double*)ffl.apos,  -frange, frange, "apos",  true );
                    berr|= ckeckRange( ffl.nvecs, 3, (double*)ffl.fapos, -frange, frange, "fapos", true );
                    berr|= ckeckRange( ffl.nvecs, 3, (double*)ffl.vapos, -frange, frange, "vapos", true );
                    if(berr){ printf( "ERROR eval_MMFFf4_ocl().debug outOfRange(%g) in replica[%i].apos  => Exit() \n", frange, isys ); exit(0); };

                    berr=false;
                    berr|= ckeckNaN_d( ffl.nvecs, 3, (double*)ffl.apos,  "apos",  true );
                    berr|= ckeckNaN_d( ffl.nvecs, 3, (double*)ffl.fapos, "fapos", true );
                    berr|= ckeckNaN_d( ffl.nvecs, 3, (double*)ffl.vapos, "vapos", true );
                    if(berr){ printf( "ERROR eval_MMFFf4_ocl().debug NaNs in replica[%i].apos  => Exit() \n", isys ); exit(0); };
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
    this->download(false,false);
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

// Minimal UFF GPU relaxation loop analogous to run_ocl_opt
int run_uff_ocl( int niter, double dt, double damping, double Fconv, double Flim ){
    if(!bUFF){ printf("MolWorld_sp3_multi::run_uff_ocl(): bUFF=false\n"); return 0; }
    if(!uff_ocl){ printf("MolWorld_sp3_multi::run_uff_ocl(): uff_ocl==nullptr\n"); return 0; } // This should not happen, uff_ocl is initialized in init()
    printf("MolWorld_sp3_multi::run_uff_ocl(  bGridFF=%i, bNonBonded=%i | niter=%i, dt=%f, damping=%f, Fconv=%f, Flim=%f) \n",bGridFF, bNonBonded,  niter, dt, damping, Fconv, Flim);
    // --- Prepare MD parameters for all systems, matching the kernel's expectation
    // UFF.cl -> updateAtomsMMFFf4 -> const float4 MDpars  = MDparams[iS]; // (dt,damp,Flimit)
    // The kernel uses MDpars.z as a multiplicative velocity damping factor. CPU uses (1.0-damping).
    for(int i=0; i<nSystems; i++){
        MDpars[i].x = dt;               // dt
        MDpars[i].y = Flim;             // Flimit
        MDpars[i].z = 1.0f - damping;   // velocity damping factor (1.0 - damping)
        MDpars[i].w = 0.0f;             // unused
        TDrive[i]   = (Quat4f){0.0f,0.0f,0.0f,0.0f}; // No thermal drive for relaxation
    }
    uff_ocl->upload( uff_ocl->ibuff_MDpars, MDpars );
    uff_ocl->upload( uff_ocl->ibuff_TDrive, TDrive );
    // --- Setup atom update kernel once before the loop
    // The user has corrected OCL_UFF.h, this signature is now correct.
    uff_ocl->setup_updateAtomsMMFFf4( uff_ocl->nAtoms, 0 ); // For UFF, nNode is 0
    // --- Run relaxation loop
    bBonding = ffu.bDoBond;
    for(int i=0; i<niter; i++){
        if     (bBonding)  { uff_ocl->eval(false);                     }
        else               { uff_ocl->task_clear_fapos->enque();       }
        if     (bGridFF   ){ uff_ocl->task_NBFF_Grid_Bspline->enque(); }
        else if(bNonBonded){ uff_ocl->task_NBFF             ->enque(); }
        uff_ocl->task_updateAtoms->enque(); 
        // --- Save current frame if trajectory is requested
        if((savePerNsteps>0) && (i%savePerNsteps==0) && (trj_fname)){
            download_uff( true, false ); // download positions for this step
            for(int isys=0; isys<nSystems; ++isys){
                char fname[256];
                sprintf(fname, "%s_%03i.xyz", trj_fname, isys);
                unpack_uff_system( isys, ffu, true, false );
                char comment[256];
                sprintf(comment, "# step %d", i+1);
                saveXYZ_system(isys, fname, comment, "a"); // Append frame
            }
        }
    }
    return niter;
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
    //printf("MolWorld_sp3_multi::run_ocl_opt() niter=%i bGroups=%i ocl.nGroupTot=%i \n", niter, bGroups, ocl.nGroupTot );
    //for(int i=0;i<npbc;i++){ printf( "CPU ipbc %i shift(%7.3g,%7.3g,%7.3g)\n", i, pbc_shifts[i].x,pbc_shifts[i].y,pbc_shifts[i].z ); }
    //debug_eval(); return 0;

    double F2conv = Fconv*Fconv;
    picked2GPU( ipicked,  1.0 );

    int err=0;
    if( task_MMFF==0)setup_MMFFf4_ocl();

    //int nPerVFs = 1;
    nPerVFs = _min(10,niter);
    //int nPerVFs = _min(50,niter);
    int nVFs    = niter/nPerVFs;
    long T0     = getCPUticks();
    int niterdone=0;
    double F2=0;

    bool dovdW=true;
    //bool dovdW=false;
    ocl.bSubtractVdW=dovdW;

    if(ocl.nGroupTot<=0){ bGroups = false; };
    bool bExplore = false;
    //for(int isys=0; isys<nSystems; isys++){ if(gopts[isys].bExploring) bExplore = true; }
    bool bBspline = (gridFF.mode==GridFFmod::BsplineFloat) || (gridFF.mode==GridFFmod::BsplineDouble);
    for(int i=0; i<nVFs; i++){

        if(bGroups){
            ocl.upload( ocl.ibuff_gtorqs,  gtorqs  );
            ocl.upload( ocl.ibuff_gforces, gforces );
        }

        bExplore = false;
        for(int isys=0; isys<nSystems; isys++){
            if(gopts[isys].bExploring) bExplore = true;
        }

        if(bAnimManipulation){ animate(); }
        if(bMILAN){
            bExplore = updateMultiExploring( Fconv );
        }
        bool bGroupDrive = bGroups && bExplore;
        //bGroupDrive = false;

        //printf( "CPU::bbox(%g,%g,%g)(%g,%g,%g)(%g,%g,%g)\n", bbox.a.x,bbox.a.y,bbox.a.z,   bbox.b.x,bbox.b.y,bbox.b.z,   bbox.c.x,bbox.c.y,bbox.c.z );
        //for(int ia=0; ia<ffl.natoms; ia++){      if( ffl.constr[ia].w > 0 ) printf( "CPU:atom[%i] constr(%g,%g,%g|%g) constrK(%g,%g,%g|%g)\n", ia, ffl.constr[ia].x,ffl.constr[ia].y,ffl.constr[ia].z,ffl.constr[ia].w,   ffl.constrK[ia].x,ffl.constrK[ia].y,ffl.constrK[ia].z,ffl.constrK[ia].w  ); }
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
            //printf("run_ocl_opt[i=%i,j=%i]\n",  i, j );
            //if(bGroupDrive)printf( "bGroupDrive==true\n" );
            {
                if( bGroupDrive )err |= task_GroupUpdate->enque_raw();
                if(dovdW)[[likely]]{
                    if(bSurfAtoms)[[likely]]{
                        if  (bGridFF)[[likely]]{
                            if(bBspline)[[likely]]{
                                //printf( " task_NBFF_Grid_Bspline ->enque_raw(); \n" );
                                if( ocl.bUseTexture ){
                                    err |= task_NBFF_Grid_Bspline_tex ->enque_raw(); //OCL_checkError(err, "task_NBFF_Grid->enque_raw(); ");
                                }else{
                                    err |= task_NBFF_Grid_Bspline ->enque_raw(); //OCL_checkError(err, "task_NBFF_Grid->enque_raw(); ");
                                }
                            }else{
                                err |= task_NBFF_Grid ->enque_raw();   //OCL_checkError(err, "task_NBFF_Grid->enque_raw(); ");
                            }
                        }else {
                            //printf( "task_NBFF(), task_SurfAtoms() \n" );
                            err |= task_NBFF     ->enque_raw();  //OCL_checkError(err, "MolWorld_sp3_multi::run_ocl_opt().task_NBFF()" );
                            err |= task_SurfAtoms->enque_raw();  //OCL_checkError(err, "MolWorld_sp3_multi::run_ocl_opt().task_SurfAtoms()" );
                        }
                    }else{
                        if(bExclusion2){
                            err |= task_NBFF_ex2->enque_raw();
                        }else{
                            err |= task_NBFF      ->enque_raw();     //OCL_checkError(err, "task_NBFF->enque_raw();");
                        }
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
        { // ======= Check vs CPU
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
            printf( "run_ocl_opt(nSys=%i|iPara=%i,bSurfAtoms=%i,bGridFF=%i,bUseTexture=%i,bExplore=%i,) CONVERGED in %i/%i steps, |F|(%g)<%g time %g [ms]( %g [us/step]) bGridFF=%i \n", nSystems, iParalel, bSurfAtoms, bGridFF, ocl.bUseTexture, bExplore, niterdone,niter, sqrt(F2), Fconv, t*1000, t*1e+6/niterdone, bGridFF );
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

    if(bMILAN){ checkBordersOfBbox(); }
    double t=(getCPUticks()-T0)*tick2second;
    if(verbosity>0)printf( "run_ocl_opt(nSys=%i|iPara=%i,bSurfAtoms=%i,bGridFF=%i,bExplore=%i) NOT CONVERGED in %i steps, |F|(%g)>%g time %g [ms]( %g [us/step]) bGridFF=%i iSysFMax=%i dovdW=%i \n", nSystems, iParalel, bSurfAtoms, bGridFF, bExplore,  niter, sqrt(F2), Fconv, t*1000, t*1e+6/niterdone, bGridFF, iSysFMax, dovdW );
    //if(database->getNMembers()>0)    printf("%i  converged: %s\n", database->getNMembers(), database->convergedStructure.back() ? "true" : "false");
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


virtual void animate(){
    int err=0;
    double phase = sin( icurIter * anim_speed );
    Vec3f d = (Vec3f)( anim_vec * phase );
    for(int isys=0; isys<nSystems; isys++){
        int i0a   = isys * ocl.nAtoms;
        for(int ia=0; ia<ocl.nAtoms; ia++){
            if( constr[ia+i0a].w >0 ) constr[ia+i0a].f.add( d );
            //constrK[i+i0a] = (Quat4f)ffl.constrK[i];
        }
    }
    err|= ocl.upload( ocl.ibuff_constr,  constr );
    //err|= ocl.upload( ocl.ibuff_constrK, constrK );
}

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

bool first=true;
uint64_t zeroT=0;

virtual void MDloop( int nIter, double Ftol = -1 ) override {
    if(iParalel<-100){ iParalel=iParalel_default; };

    if(Ftol<0)  Ftol = Ftol_default;
    //printf( "MolWorld_sp3_multi::MDloop(%i) bGridFF %i bOcl %i bMMFF %i iParalel=%i \n", nIter, bGridFF, bOcl, bMMFF, iParalel );
    //bMMFF=false;
    if(first){ zeroT = getCPUticks(); first=false; }

    // if(bAnimManipulation){
    //     float dir_sign = ((nloop/10000)%2)*2 - 1 ; //printf("AnimManipulation dir_sign=%f \n", dir_sign);
    //     move_MultiConstrain( Vec3d{ 0.2*dir_sign, 0.0, 0.0}, Vec3dZero );
    // }

    // //run_omp_ocl( 1 );
    // run_omp_ocl( 50 );
    // //run_omp_ocl( 100 );
    // return;

    nIter=iterPerFrame;
    //bOcl=iParalel>0;
    int nitrdione=0;
    switch(iParalel){
        case -1: bOcl=false;  nitrdione = run_multi_serial(nIter,Ftol);  break;
        case  0: bOcl=false;  nitrdione = MolWorld_sp3::run_no_omp( nIter, Ftol ); break;
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

    if( bMILAN && bSaveToDatabase  ){ // Milan
        std::ofstream file("minima.dat", std::ios::app);
        if(icurIter%1000 < 200){
            if(icurIter<1000) {
                printf("%15s %20s %20s %20s %20s %20s %20s %20s %20s %20s\n", "Nb Iteration", "Totally", "Unique converged", "Time", "nbConverged", "nbEvaluation", "nStepConvAvg", "nStepNonConvAvg", "nStepExploringAvg", "total evaluation");
                file << "Iterations" << " " << "Totally" << " " << "Unique" << " " << "Time"  << "\n";
            }
            // /((double)(nbConverged+nbNonConverged))"nStepSaveAvg",
            printf("%15d %20d %20d %20g %20d %20d %20g %20g %20g %20d \n", icurIter, database->totalEntries, database->getNMembers(), (getCPUticks()-zeroT)*tick2second, nbConverged, nbEvaluation*nPerVFs, (nStepNonConvSum+nStepConvSum)/(double)(nbConverged+nbNonConverged), nStepNonConvSum/((double)nbNonConverged), nStepExplorSum/((double)nExploring), (nStepConvSum+nStepNonConvSum+nStepExplorSum) );
            file << icurIter << " " << database->totalEntries << " " << database->getNMembers() << " " <<  (getCPUticks()-zeroT)*tick2second << "\n";
        }
        file.close();  // zavření souboru
    }
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

void surf2ocl( Vec3i nPBC ){
    int err=0;
    printf( "MolWorld_sp3_multi::surf2ocl() gridFF.natoms=%i nPBC(%i,%i,%i)\n", gridFF.natoms, nPBC.x,nPBC.y,nPBC.z );
    // Debug: Check if OpenCL context is initialized
    printf( "DEBUG surf2ocl: OCL context = %p, commands = %p\n", (void*)ocl.context, (void*)ocl.commands );
    if(ocl.context == 0){
        printf( "ERROR in MolWorld_sp3_multi::surf2ocl(): OpenCL context not initialized! This will cause surf2ocl to fail. =>exit()\n" );
        exit(0);
    }
    //long T0=getCPUticks();
    Quat4f* atoms_surf = new Quat4f[gridFF.natoms];
    Quat4f* REQs_surf  = new Quat4f[gridFF.natoms];
    pack( gridFF.natoms, gridFF.apos, atoms_surf, sq(gridFF.Rdamp) );
    pack( gridFF.natoms, gridFF.REQs, REQs_surf                     );   // ToDo: H-bonds should be here
    ocl.GFFparams.x = gridFF.Rdamp;
    ocl.GFFparams.y = gridFF.alphaMorse;
    //v2f4( gridFF.grid.pos0, ocl.grid_p0 );
    ocl.grid_p0.f     = (Vec3f)gridFF.grid.pos0;
    ocl.grid_shift0.f = (Vec3f)gridFF.shift0;
    //printf( "grid_p0(%g,%g,%g) grid_shift0(%g,%g,%g)\n", ocl.grid_p0.x,ocl.grid_p0.y,ocl.grid_p0.z,    ocl.grid_shift0.x,ocl.grid_shift0.y,ocl.grid_shift0.z  );
    //exit(0);
    ocl.surf2ocl(  gridFF.grid, nPBC, gridFF.natoms, (float4*)atoms_surf, (float4*)REQs_surf );
    delete [] atoms_surf;
    delete [] REQs_surf;
    bDONE_surf2ocl = true;
}

void surf2ocl_uff( Vec3i nPBC ){
    int err=0;
    printf( "MolWorld_sp3_multi::surf2ocl_uff() gridFF.natoms=%i nPBC(%i,%i,%i)\n", gridFF.natoms, nPBC.x,nPBC.y,nPBC.z );
    // Debug: Check if OpenCL context is initialized
    printf( "DEBUG surf2ocl_uff: OCL context = %p, commands = %p\n", (void*)uff_ocl->context, (void*)uff_ocl->commands );
    if(uff_ocl->context == 0){
        printf( "ERROR in MolWorld_sp3_multi::surf2ocl_uff(): OpenCL context not initialized! This will cause surf2ocl_uff to fail. =>exit()\n" );
        exit(0);
    }
    //long T0=getCPUticks();
    Quat4f* atoms_surf = new Quat4f[gridFF.natoms];
    Quat4f* REQs_surf  = new Quat4f[gridFF.natoms];
    pack( gridFF.natoms, gridFF.apos, atoms_surf, sq(gridFF.Rdamp) );
    pack( gridFF.natoms, gridFF.REQs, REQs_surf                     );   // ToDo: H-bonds should be here
    uff_ocl->GFFparams.x = gridFF.Rdamp;
    uff_ocl->GFFparams.y = gridFF.alphaMorse;
    //v2f4( gridFF.grid.pos0, ocl.grid_p0 );
    uff_ocl->grid_p0.f     = (Vec3f)gridFF.grid.pos0;
    uff_ocl->grid_shift0.f = (Vec3f)gridFF.shift0;
    //printf( "grid_p0(%g,%g,%g) grid_shift0(%g,%g,%g)\n", ocl.grid_p0.x,ocl.grid_p0.y,ocl.grid_p0.z,    ocl.grid_shift0.x,ocl.grid_shift0.y,ocl.grid_shift0.z  );
    //exit(0);
    uff_ocl->surf2ocl(  gridFF.grid, nPBC, gridFF.natoms, (float4*)atoms_surf, (float4*)REQs_surf );
    delete [] atoms_surf;
    delete [] REQs_surf;
    bDONE_surf2ocl = true;
}

void completeGridFFInit(){
    if(!bGridFF_pending) return;

    printf("MolWorld_sp3_multi::completeGridFFInit() - Completing deferred GridFF initialization\n");

    // First initialize the GridFF properly
    printf("Initializing GridFF with name: %s\n", gridFF_name);
    initGridFF(gridFF_name, gridFF_z0, gridFF_cel0, gridFF_bSymetrize);

    // Now call the appropriate surf2ocl with the properly initialized GridFF
    if(bUFF){
        surf2ocl_uff(gridFF.nPBC);
    }else{
        surf2ocl(gridFF.nPBC);
    }

    bGridFF_pending = false;
}

// Override loadSurf to handle deferred GridFF initialization
bool loadSurf(const char* name, bool bGrid=true, bool bSaveDebugXSFs=false, double z0=NAN, Vec3d cel0={-0.5,-0.5,0.0}) {
    printf("MolWorld_sp3_multi::loadSurf(%s) - deferring GridFF initialization\n", name);

    // Store parameters for deferred initialization
    bGridFF_pending = bGrid;
    gridFF_name = name;
    gridFF_z0 = z0;
    gridFF_cel0 = cel0;
    gridFF_bSymetrize = true; // Default value

    // Call base class to load the surface but don't initialize GridFF yet
    bool result = MolWorld_sp3::loadSurf(name, false, bSaveDebugXSFs, z0, cel0);

    return result;
}

 void make_gridFF_ocl( bool bSaveDebug=false ){
    if(bDONE_surf2ocl==false){ printf("ERROR in MolWorld_sp3_multi::make_gridFF_ocl() call surf2ocl() first \n"); exit(0); };
    int err=0;
    long T1=getCPUticks();
    ocl.makeGridFF( true );
    err |=  ocl.finishRaw();    OCL_checkError(err, "make_gridFF_ocl.makeGridFF.finish");
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
        err |=  ocl.finishRaw();    OCL_checkError(err, "make_gridFF_ocl.download.finish");
        printf( ">>time(make_gridFF_ocl.download() %g \n", (getCPUticks()-T1)*tick2second );
        if(bSaveDebug){ gridFF.saveXsfDebug(); }
    }
    //  ToDo: We do no release these buffers because we use them for kernell:: getSurfMorse()
    //ocl.buffers[ocl.ibuff_atoms_surf].release();
    //ocl.buffers[ocl.ibuff_REQs_surf ].release();
    err |=  ocl.finishRaw();    OCL_checkError(err, "make_gridFF_ocl.makeGridFF.atoms_surf.release");
    //bDONE_surf2ocl = true;
    //exit(0);
}

//virtual double* initGridFF( const char * name, bool bGrid=true, bool bSaveDebugXSFs=false, double z0=NAN, Vec3d cel0={-0.5,-0.5,0.0}, bool bAutoNPBC=true, bool bCheckEval=true )override{
virtual void initGridFF( const char * name, double z0=NAN, Vec3d cel0={-0.5,-0.5,0.0}, bool bSymetrize=true, bool bAutoNPBC=true, bool bCheckEval=true, bool bUseEwald=true, bool bFit=true, bool bRefine=true ) override {
    int err=0;
    printf( "MolWorld_sp3_multi::initGridFF(%s) bSymetrize=%i bAutoNPBC=%i bCheckEval=%i bUseEwald=%i bFit=%i bRefine=%i ,cel0={%g,%g,%g}) \n", name,  bSymetrize, bAutoNPBC, bCheckEval, bUseEwald, bFit, bRefine, cel0.x,cel0.y,cel0.z );

    // Call base class implementation first
    MolWorld_sp3::initGridFF(name, z0, cel0, bSymetrize, bAutoNPBC, bCheckEval, bUseEwald, bFit, bRefine);

    // Check if OpenCL is initialized before calling surf2ocl
    if(bUFF){
        if(uff_ocl == 0 || uff_ocl->context == 0){
            printf("WARNING: MolWorld_sp3_multi::initGridFF() UFF OpenCL context not initialized yet. Deferring surf2ocl call.\n");
            bGridFF_pending = true;
        } else {
            // OpenCL is initialized, call surf2ocl immediately
            surf2ocl_uff(gridFF.nPBC);
        }
    }else{
        if(ocl.context == 0){
            printf("WARNING: MolWorld_sp3_multi::initGridFF() MMFF OpenCL context not initialized yet. Deferring surf2ocl call.\n");
            bGridFF_pending = true;
        } else {
            // OpenCL is initialized, call surf2ocl immediately
            surf2ocl(gridFF.nPBC);
        }
    }
    // if(verbosity>0)printf("MolWorld_sp3_multi::initGridFF(%s,bGrid=%i,z0=%g,cel0={%g,%g,%g})\n",  name, z0, cel0.x,cel0.y,cel0.z  );
    // if(gridFF.grid.n.anyEqual(0)){ printf("ERROR in MolWorld_sp3_multi::initGridFF() zero grid.n(%i,%i,%i) => Exit() \n", gridFF.grid.n.x,gridFF.grid.n.y,gridFF.grid.n.z ); exit(0); };
    // gridFF.grid.center_cell( cel0 );
    // bGridFF=true;
    // gridFF.bindSystem        ( surf  .natoms, surf  .atypes, surf.apos  , surf.REQs        );
    // gridFF.setAtomsSymetrized( gridFF.natoms, gridFF.atypes, gridFF.apos, gridFF.REQs, 0.1 );
    // //gridFF.setAtomsSymetrized(surf.natoms, surf.atypes, surf.apos, surf.REQs );
    // gridFF.evalCellDipole();
    // if( ( fabs(gridFF.Q)>1e-6 ) || (gridFF.dip.norm2()>1e-8) ){ printf("ERROR: GridFF has dipole and dipole correction not yet implemented => exit() \n"); exit(0); }
    // if( isnan(z0) ){ z0=gridFF.findTop();   if(verbosity>0) printf("GridFF::findTop() %g \n", z0);  };
    // gridFF.grid.pos0.z=z0;
    // gridFF.lvec = gridFF.grid.cell;
    // //if(verbosity>1)
    // //gridFF.grid.printCell();
    // gridFF.nPBC=Vec3i{1,1,0};
    // if(bAutoNPBC){ autoNPBC( gridFF.grid.cell, gridFF.nPBC, 20.0 ); }
    // gridFF.makePBCshifts   ( gridFF.nPBC,      gridFF.lvec       );
    // long T0 = getCPUticks();
    // //bSurfAtoms=false;
    // //gridFF.shift0 = Vec3d{0.,0., 0.0};
    // gridFF.shift0 = Vec3d{0.,0.,-2.0};
    // //gridFF.shift0 = Vec3d{0.,0.,-6.0};
    // double bSaveDebugXSFs=false;

    bool bOnGPU=false; // Currently we have not implemented B-spline GridFF on GPU

     //if(verbosity>0)
    printf("MolWorld_sp3_multi::initGridFF(%s,bSymetrize=%i,bAutoNPBC=%i bCheckEval=%i bUseEwald=%i bFit=%i bRefine=%i bGridDouble=%i gridStep=%g,z0=%g,cel0={%g,%g,%g} )\n",  name, bSymetrize, bAutoNPBC,bCheckEval,bUseEwald,bFit,bRefine, bGridDouble, gridStep, z0, cel0.x,cel0.y,cel0.z  );
    sprintf(tmpstr, "%s.lvs", name );
    if( file_exist(tmpstr) ){  gridFF.grid.loadCell( tmpstr, gridStep );  gridFF.bCellSet=true; }
    if( !gridFF.bCellSet ){
        bGridFF=false;
        printf( "WARRNING!!! GridFF not initialized because %s not found\n", tmpstr );
        return;
    }
    //double* ffgrid = 0;
    gridFF.grid.center_cell( cel0 );
    bGridFF=true;
    gridFF.bindSystem(surf.natoms, surf.atypes, surf.apos, surf.REQs );
    gridFF.initGridFF( name, z0, bAutoNPBC, bSymetrize );
    char wd0[1024]; getcwd(wd0,1024); printf( "MolWorld_sp3::initGridFF() 1 wd0=`%s`\n", wd0 );
    const char* last_slash = strrchr(name, '/');
    const char* result = (last_slash) ? last_slash + 1 : name;
    char wd[128]; sprintf( wd, "data/%s", result); printf( "MolWorld_sp3::initGridFF() 1 wd=`%s`\n", wd );
    tryMakeDir  ( wd );
    tryChangeDir( wd );
    gridFF.bUseEwald = bUseEwald;
    gridFF.ewald     = &gewald;
    //ffgrid = gridFF.HHermite_d;
    //getcwd(tmpstr, 1024 ); printf( "initGridFF() 3 WD=`%s`\n", tmpstr );
    gridFF.shift0 = Vec3d{0.,0.,-2.0};
    //gridFF.shift0 = Vec3d{0.,0.,0.0};
    //if(bCheckEval)gridFF.evalCheck();    // WARRNING:  CHECK FOR gridFF TURNED OFF !!!!!!!!!!!!!!!!!!!!!!!!!
    //return ffgrid;

    if(bOnGPU){
        bool bSaveDebugXSFs=false;
        make_gridFF_ocl( bSaveDebugXSFs );
        //surf2ocl( gridFF.nPBC, bSaveDebugXSFs );
        //printf( ">>time(init_ocl;GridFF_ocl): %g [s] \n", (getCPUticks()-T0)*tick2second  );
        bGridFF   =true;
        // evalCheckGridFF_ocl();   // here are not initialized buffers atoms.aforce,REQs, so it will crash.
        //return 0;
    }else{
        gridFF.tryLoad_new( bSymetrize, bFit, bRefine );
        bGridFF   =true;
        //  copy to GPU - use the correct OpenCL context based on bUFF flag
        DEBUG
        GridShape& gsh = gridFF.grid;

        // Select the correct OpenCL object based on bUFF flag
        OCLsystem* target_ocl = bUFF ? (OCLsystem*)uff_ocl : (OCLsystem*)&ocl;

        if(bUFF){
            // For UFF, set grid parameters in uff_ocl
            uff_ocl->grid_step    = Quat4f{ (float)gsh.dCell.xx, (float)gsh.dCell.yy, (float)gsh.dCell.zz, 0.0 };
            uff_ocl->grid_invStep.f.set_inv( uff_ocl->grid_step.f );
        }else{
            // For MMFF, set grid parameters in ocl
            ocl.grid_step    = Quat4f{ (float)gsh.dCell.xx, (float)gsh.dCell.yy, (float)gsh.dCell.zz, 0.0 };
            ocl.grid_invStep.f.set_inv( ocl.grid_step.f );
        }

        DEBUG
        int nxyz = gridFF.Bspline_to_f4(true);
        DEBUG
        if(bUFF ? uff_ocl->bUseTexture : ocl.bUseTexture){
            printf("Uploading B-spline data to 3D Texture (itex_BsplinePLQH).\n");
            if(bUFF){
                uff_ocl->itex_BsplinePLQH = uff_ocl->newBufferImage3D( "BsplinePLQH_tex",  gsh.n.x, gsh.n.y, gsh.n.z, sizeof(cl_float4), gridFF.Bspline_PLQf, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, {CL_RGBA, CL_FLOAT} );
                DEBUG
                if(uff_ocl->itex_BsplinePLQH < 0) { printf("Error creating BsplinePLQH_tex texture.\n"); exit(0); }
            }else{
                ocl.itex_BsplinePLQH = ocl.newBufferImage3D( "BsplinePLQH_tex",  gsh.n.x, gsh.n.y, gsh.n.z, sizeof(cl_float4), gridFF.Bspline_PLQf, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, {CL_RGBA, CL_FLOAT} );
                DEBUG
                if(ocl.itex_BsplinePLQH < 0) { printf("Error creating BsplinePLQH_tex texture.\n"); exit(0); }
            }
        }else{
            if(bUFF){
                printf("DEBUG MolWorld_sp3_multi::initGridFF() - Allocating BsplinePLQ buffer for UFF, uff_ocl=%p\n", uff_ocl);
                uff_ocl->ibuff_BsplinePLQ = uff_ocl->newBuffer( "BsplinePLQ", nxyz, sizeof(float4), 0, CL_MEM_READ_ONLY  );
                printf("DEBUG MolWorld_sp3_multi::initGridFF() - After newBuffer, ibuff_BsplinePLQ=%d\n", uff_ocl->ibuff_BsplinePLQ);
                err = uff_ocl->upload( uff_ocl->ibuff_BsplinePLQ, gridFF.Bspline_PLQf ); OCL_checkError(err, "MolWorld_sp3_multi::initGridFF().upload(BsplinePLQ)");
                printf("DEBUG MolWorld_sp3_multi::initGridFF() - After upload, ibuff_BsplinePLQ=%d\n", uff_ocl->ibuff_BsplinePLQ);
            }else{
                ocl.ibuff_BsplinePLQ = ocl.newBuffer( "BsplinePLQ", nxyz, sizeof(float4), 0, CL_MEM_READ_ONLY  );
                err = ocl.upload( ocl.ibuff_BsplinePLQ, gridFF.Bspline_PLQf ); OCL_checkError(err, "MolWorld_sp3_multi::initGridFF().upload(BsplinePLQ)");
            }
            { // Debug
                for(int i=0; i<nxyz; i++){ gridFF.Bspline_PLQf[i]=Quat4fZero; }
                if(bUFF){
                    err = uff_ocl->download( uff_ocl->ibuff_BsplinePLQ, gridFF.Bspline_PLQf);  OCL_checkError(err, "MolWorld_sp3_multi::initGridFF().download(BsplinePLQ)");
                }else{
                    err = ocl.download( ocl.ibuff_BsplinePLQ, gridFF.Bspline_PLQf);  OCL_checkError(err, "MolWorld_sp3_multi::initGridFF().download(BsplinePLQ)");
                }
                Vec3d  min1=Vec3dZero,max1=Vec3dZero;
                Quat4f min2=Quat4fZero,max2=Quat4fZero;
                for(int i=0; i<nxyz; i++){
                    gridFF.Bspline_PLQf[i].update_bounds( min2, max2 );
                    gridFF.Bspline_PLQ [i].update_bounds( min1, max1 );
                }
                printf( "Bspline_PLQ  min(%g,%g,%g   ) max(%g,%g,%g   ) \n", min1.x,min1.y,min1.z,        max1.x,max1.y,max1.z        );
                printf( "Bspline_PLQf min(%g,%g,%g,%g) max(%g,%g,%g,%g) \n", min2.x,min2.y,min2.z,min2.w, max2.x,max2.y,max2.z,max2.w );
            }
        }
        DEBUG
        //exit(0);
    }
    tryChangeDir( wd0 );
    printf( "MolWorld_sp3_multi::initGridFF() DONE\n" );
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

    //printf("GPU AFTER assemble() \n"); ff4.printDebug( false,false );
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

    //printf("GPU AFTER assemble() \n"); ff4.printDebug( false,false );
    unpack( ff4.natoms, ffl.  apos, ff4.  apos );
    unpack( ff4.natoms, ffl. fapos, ff4. fapos );
    //opt.move_FIRE();
    //ff4.printDebug( false, false );
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
