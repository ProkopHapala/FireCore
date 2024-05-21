#ifndef  OCL_MM_h
#define  OCL_MM_h

/*
This is a OpenCL interface for Classical Forcefield such as  MMFFf4 (molecular mechanics with 4 neighbor), GridFF, and PP-AFM (probe particle AFM)
It can simulate multiple replica of systems in parallel in order to exploit GPU parallelism and use all cores of the GPU even when simulating small systems (<100-1000 atoms)
Non-covalent interactions are calculated in periodic boundary conditions (PBC) excluding interactions between bonded atoms which is recognized from list of neighbors and its cell indexes
Non-covalent interactions between dynamical atoms consist of Lennard-Jones (LJ) and Coulomb interactions and Hydrogen bond correction (pseudo-charge H)
Non-covalent with rigid substrate (Pauli repulsion, London dispersion, and Coulomb electrostatics) is calculated using GridFF by 3D texture interpolation
*/

//#include "OCL_DFT.h"
#include "OCL.h"

// coversions between OpenCL and C++ types
void v2f4  ( const Vec3d& v, float4& f4 ){ f4.x=(float)v.x; f4.y=(float)v.y; f4.z=(float)v.z; };                      // pack Vec3d to float4
void f4toV3( const cl_float4& f4, Vec3d& v ){ v.x=f4.s[0]; v.y=f4.s[1]; v.z=f4.s[2]; };                               // unpack float4 to Vec3d
void v2i4( const Vec3i& v, int4& f4   ){ f4.x=(float)v.x; f4.y=(float)v.y; f4.z=(float)v.z; };                        // pack Vec3i to int4
//void v2f4( const Vec3d& v, cl_float4& f4 ){ f4.s[0]=(cl_float)v.x; f4.s[1]=(cl_float)v.y; f4.s[2]=(cl_float)v.z; }; // pack Vec3d to float4
cl_float4 cl_f4( const Vec3d& v ){ return (cl_float4){(cl_float)v.x,(cl_float)v.y,(cl_float)v.z,0.f}; };              // pack Vec3d to float4

// see https://stackoverflow.com/questions/33119233/opencl-using-struct-as-kernel-argument
typedef struct __attribute__ ((packed)) cl_Mat3{
    cl_float4 a;
    cl_float4 b;
    cl_float4 c;
}st_foo;


inline void printMat( const cl_Mat3& mat  ){
	printf( " %f %f %f \n", mat.a.s[0], mat.a.s[1], mat.a.s[2] );
	printf( " %f %f %f \n", mat.b.s[0], mat.b.s[1], mat.b.s[2] );
	printf( " %f %f %f \n", mat.c.s[0], mat.c.s[1], mat.c.s[2] );
}

// pack Mat3d to OpenCL cl_Mat3
void Mat3_to_cl( const Mat3d& m, cl_Mat3& clm ){ 
    //v2f4( m.a, clm.a ); v2f4( m.b, clm.b ); v2f4( m.c, clm.c ); 
    clm.a = cl_f4( m.a );
    clm.b = cl_f4( m.b );
    clm.c = cl_f4( m.c );
}

// unpack OpenCL cl_Mat3 to Mat3d
void Mat3_from_cl( Mat3d& m, const cl_Mat3& clm ){ 
    f4toV3( clm.a, m.a );
    f4toV3( clm.b, m.b );
    f4toV3( clm.c, m.c );
}

//=======================================================================
//=======================================================================

class OCL_MM: public OCLsystem { public:
    cl_program program_relax=0;

    // dimensions
    int nAtoms=0; 
    int nnode=0, nvecs=0, nneigh=0, npi=0, nSystems=0, nbkng=0, ncap=0; // nnode: number of node atoms; nvecs: number of vectors (atoms and pi-orbitals); nneigh: number of neighbors); npi: number of pi orbitals; nSystems: number of system replicas; nbkng: number of back-neighbors; ncap: number of capping atoms

    int4   print_mask{1,1,1,1};  // what to print on GPU
    int4   nDOFs    {0,0,0,0};   // number of DOFs (nAtoms,nnode,nSystems,0)
    int4   nPBC     {0,0,0,0};   // number of PBC images in each direction (nx,ny,nz,0)
    int    npbc=0;               // total number of PBC images
    int    bSubtractVdW=1;       // subtract non-bonded interactions for atoms bonnded to common neighbor ?
    //float4 md_params{0.05,0.9,100.0,0.0};     // (dt,cdamp,forceLimit)
    float4 md_params{0.05,0.98,100.0,0.0};      // (dt,cdamp,forceLimit)
    //float4 md_params{0.05,0.995,100.0,0.0};   // (dt,cdamp,forceLimit)
    
    int    niter;                    // number of iterations in relax
    float4 GFFparams{1.0,1.5,0.,0.}; // (Rdamp, alphaMorse, 0, 0) for GridFF

    cl_Mat3 cl_lvec;      // lattice vectors         - DEPRECATED: but we do not need them anymore, we use arrays of lattice vectors for each system
    cl_Mat3 cl_invLvec;   // inverse lattice vectors - DEPRECATED: but we do not need them anymore, we use arrays of lattice vectors for each system

    // OpenCL buffers and textures ids
    int ibuff_atoms=-1,ibuff_aforces=-1,ibuff_neighs=-1,ibuff_neighCell=-1;
    int ibuff_avel=-1,ibuff_cvf=-1, ibuff_neighForce=-1,  ibuff_bkNeighs=-1, ibuff_bkNeighs_new=-1;
    int ibuff_REQs=-1, ibuff_MMpars=-1, ibuff_BLs=-1,ibuff_BKs=-1,ibuff_Ksp=-1, ibuff_Kpp=-1;   // MMFFf4 params
    int ibuff_lvecs=-1, ibuff_ilvecs=-1,ibuff_MDpars=-1,ibuff_TDrive=-1, ibuff_pbcshifts=-1; 
    int ibuff_constr=-1;
    int ibuff_samp_ps=-1;
    int ibuff_samp_fs=-1;
    int ibuff_samp_REQ=-1;
    // OpenCL buffers and textures ids
    int itex_FE_Paul=-1;
    int itex_FE_Lond=-1;
    int itex_FE_Coul=-1;
    int ibuff_atoms_surf=-1;
    int ibuff_REQs_surf =-1;
    int ibuff_dipole_ps=-1;
    int ibuff_dipoles  =-1;

    // ------- Grid
    //size_t Ns[4]; // = {N0, N1, N2};
    //size_t Ntot;
    int natom_surf=0;  // number of atoms in surface
    int4    grid_n;    // number of grid points in each direction (nx,ny,nz,0)
    //float4  grid_p0       { 0.f, 0.f, 0.f, 0.f };
    //float4  grid_shift0   { 0.f, 0.f, 0.f, 0.f };
    //float4  grid_shift0_p0{ 0.f, 0.f, 0.f, 0.f };

    //Quat4f  surf_p0       { 0.f, 0.f, 0.f, 0.f }; // surface shift0
    Quat4f  grid_p0       { 0.f, 0.f, 0.f, 0.f }; // grid origin
    Quat4f  grid_shift0   { 0.f, 0.f, 0.f, 0.f }; // grid shift
    Quat4f  grid_shift0_p0{ 0.f, 0.f, 0.f, 0.f }; // grid shift + origin

    cl_Mat3 cl_dGrid;       // grid cell step (voxel rhombus)
    cl_Mat3 cl_diGrid;      // inverse grid cell step
    cl_Mat3 cl_grid_lvec;   // grid lattice vectors
    cl_Mat3 cl_grid_ilvec;  // inverse grid lattice vectors
    
    // =================== PP_AFM - Probe Particle AFM variables

        int itex_FEAFM      =-1;     // Forcefield texture with Force and Energy for AFM  
        int ibuff_afm_ps    =-1;     // buffer with AFM probe particle positions at which we sample the forcefield
        int ibuff_afm_ws    =-1;     // buffer with AFM probe particle weights (for convolution)
        int ibuff_afm_FEout =-1;     // buffer with AFM probe particle Force and Energy output
        int ibuff_afm_PPpos =-1;     // buffer with AFM probe particle positions after relaxation

        int4     afm_nPBC{0,0,0,0};        // number of PBC images in each direction (nx,ny,nz,0)
        int4     afm_grid_n{0,0,0,0};      // number of grid points in each direction (nx,ny,nz,0)
        Quat4f   afm_grid_p0{0.,0.,0.,0.}; // grid origin
        cl_Mat3  afm_grid_lvec;            // grid lattice vectors
        cl_Mat3  afm_grid_ilvec;           // inverse grid lattice vectors
        //float4 tipRot[3]={{1.,0.,0.,0.},{0.,1.,0.,0.},{0.,0.,1.,0.1}};  // tip rotation
        cl_Mat3  tipRot;                                          // tip rotation
        float4  afm_relax_params{0.5f,0.1f,0.000001f,0.5f};       // (dt,cdamp,forceLimit,forceLimit2)
        //float4  tipParams{  1.661f, 0.0091063f, -0.1f, 0.0f };  // tip stiffness
        float4  tipParams{  1.661f, 0.0091063f, 0.0f, 0.0f };     // 
        float4  tip_stiffness { -0.03f,-0.03f,-0.03f,-1.00 };  // tip stiffness
        float4  tip_dpos0{0.0f,0.0f,-4.0f, 4.0f};    
        //float4  tip_Qs {-0.0f,-0.4f,-0.0f,0.0f};
        //float4  tip_Qs {-0.4f,+0.4f,-0.0f,0.0f};
        float4  tip_Qs {-2.0f,+4.0f,-2.0f,0.0f};        // tip charge
        float4  tip_QZs{-0.1f, 0.0f,+0.1f,0.0f};        // tip charge z-position
        float4  afm_surfFF{0.f,0.f,0.f,0.f};            // surface forcefield parameters
        int     afm_nz,afm_nzout, afm_nMaxItr=128;      // number of z-slices, number of z-slices in output, max number of iterations

    // ====================== Functions


    void makeKrenels_MM( const char*  cl_src_dir ){
        printf( "makeKrenels_MM() \n" );
        char srcpath[1024];
        sprintf( srcpath, "%s/relax_multi.cl", cl_src_dir );     
        buildProgram( srcpath, program_relax );
        newTask( "getNonBond"             ,program_relax, 2);
        newTask( "getMMFFf4"              ,program_relax, 2);
        newTask( "cleanForceMMFFf4"       ,program_relax, 2);
        newTask( "updateAtomsMMFFf4"      ,program_relax, 2);
        newTask( "printOnGPU"             ,program_relax, 2);
        newTask( "getNonBond_GridFF"      ,program_relax, 2);
        newTask( "getSurfMorse"           ,program_relax, 2);
        newTask( "make_GridFF"            ,program_relax, 1);
        newTask( "sampleGridFF"           ,program_relax, 1);
        newTask( "addDipoleField"         ,program_relax, 1);
        newTask( "evalMMFFf4_local1"      ,program_relax, 2);
        newTask( "evalMMFFf4_local2"      ,program_relax, 2);
        newTask( "evalMMFFf4_local_test"  ,program_relax, 2);
        newTask( "PPAFM_makeFF"           ,program_relax, 1);
        newTask( "PPAFM_scan"             ,program_relax, 1);
        newTask( "PPAFM_scan_df"          ,program_relax, 1);

        newTask( "getSortRangeBuckets"    ,program_relax, 1); 

        //newTask( "write_toImg"     ,program_relax, 3,{0,0,0,0},{1,1,1,0} ); 
        printf( "... makeKrenels_MM() DONE \n" );
    }

    int initAtomsForces( int nSystems_, int nAtoms_, int nnode_, int npbc_ ){
        nSystems=nSystems_;
        nnode  = nnode_;
        nAtoms = nAtoms_;
        npi    = nnode_;
        nvecs  = nAtoms+npi;   // number of vectors (atoms and pi-orbitals)
        nbkng  = nnode*4*2;    // number of back-neighbors (4 neighbors per node, position and pi-orbital)
        ncap   = nAtoms-nnode; // number of capping atoms
        if(npbc_==0){ npbc=1; }; 
        npbc   = npbc_;
        printf( "initAtomsForces() nSystems %i nvecs %i natoms %i nnode %i nbkng %i \n", nSystems, nvecs, nAtoms, nnode, nbkng );
        printf( "initAtomsForces() nS*nvecs %i nS*natoms %i nS*nnode %i nS*nbkng %i \n", nSystems*nvecs,  nSystems*nAtoms, nSystems*nnode, nSystems*nbkng );
        if( (nSystems<=0)||(nAtoms<=0) ){ printf("ERROR in OCL_MM::initAtomsForces() invalid size nSystems=%i nAtoms=%i => Exit() \n", nSystems, nAtoms); exit(0); }
        ibuff_atoms      = newBuffer( "atoms",      nSystems*nvecs , sizeof(float4), 0, CL_MEM_READ_WRITE ); // atoms positions and velocities (x,y,z,m)
        ibuff_aforces    = newBuffer( "aforces",    nSystems*nvecs , sizeof(float4), 0, CL_MEM_READ_WRITE ); // atoms forces
        ibuff_REQs       = newBuffer( "REQs",       nSystems*nAtoms, sizeof(float4), 0, CL_MEM_READ_ONLY  ); // atoms parameters {R0,E0,Q}
        ibuff_neighs     = newBuffer( "neighs",     nSystems*nAtoms, sizeof(int4  ), 0, CL_MEM_READ_ONLY  );
        ibuff_neighCell  = newBuffer( "neighCell" , nSystems*nAtoms, sizeof(int4  ), 0, CL_MEM_READ_ONLY  );

        ibuff_constr     = newBuffer( "constr",    nSystems*nAtoms , sizeof(float4), 0, CL_MEM_READ_WRITE );
        //ibuff_constr0    = newBuffer( "constr0",   nSystems*nAtoms , sizeof(float4), 0, CL_MEM_READ_WRITE );
        //ibuff_constrK    = newBuffer( "constrK",   nSystems*nAtoms , sizeof(float4), 0, CL_MEM_READ_WRITE );

        ibuff_bkNeighs     = newBuffer( "bkNeighs", nSystems*nvecs,  sizeof(int4  ), 0, CL_MEM_READ_ONLY  );     // back neighbors
        ibuff_bkNeighs_new = newBuffer( "bkNeighs_new", nSystems*nvecs,  sizeof(int4  ), 0, CL_MEM_READ_ONLY  );   
        ibuff_avel       = newBuffer( "avel",       nSystems*nvecs,  sizeof(float4), 0, CL_MEM_READ_WRITE );     // atoms velocities (x,y,z,m)
        ibuff_cvf        = newBuffer( "cvf",        nSystems*nvecs , sizeof(float4), 0, CL_MEM_READ_WRITE );
        ibuff_neighForce = newBuffer( "neighForce", nSystems*nbkng,  sizeof(float4), 0, CL_MEM_READ_WRITE );

        ibuff_MMpars     = newBuffer( "MMpars",     nSystems*nnode,  sizeof(int4),   0, CL_MEM_READ_ONLY  );
        ibuff_BLs        = newBuffer( "BLs",        nSystems*nnode,  sizeof(float4), 0, CL_MEM_READ_ONLY  );
        ibuff_BKs        = newBuffer( "BKs",        nSystems*nnode,  sizeof(float4), 0, CL_MEM_READ_ONLY  );
        ibuff_Ksp        = newBuffer( "Ksp",        nSystems*nnode,  sizeof(float4), 0, CL_MEM_READ_ONLY  );
        ibuff_Kpp        = newBuffer( "Kpp",        nSystems*nnode,  sizeof(float4), 0, CL_MEM_READ_ONLY  );

        ibuff_MDpars     = newBuffer( "MDpars",     nSystems,        sizeof(float4),  0, CL_MEM_READ_ONLY  );
        ibuff_TDrive     = newBuffer( "TDrive",     nSystems,        sizeof(float4),  0, CL_MEM_READ_ONLY  );
        ibuff_lvecs      = newBuffer( "lvecs",      nSystems,        sizeof(cl_Mat3), 0, CL_MEM_READ_ONLY  );
        ibuff_ilvecs     = newBuffer( "ilvecs",     nSystems,        sizeof(cl_Mat3), 0, CL_MEM_READ_ONLY  );

        ibuff_pbcshifts  = newBuffer( "pbcshifts",  nSystems*npbc,   sizeof(float4), 0, CL_MEM_READ_ONLY  );

        // int nsamp_max = 1000; DEBUG
        // ibuff_samp_fs   = newBuffer( "samp_fs",   nsamp_max, sizeof(float4), 0, CL_MEM_READ_WRITE );   DEBUG
        // ibuff_samp_ps   = newBuffer( "samp_ps",   nsamp_max, sizeof(float4), 0, CL_MEM_READ_ONLY  );    DEBUG
        // ibuff_samp_REQ  = newBuffer( "samp_REQ",  nsamp_max, sizeof(float4), 0, CL_MEM_READ_ONLY  );    DEBUG

        return ibuff_atoms;
    }


    OCLtask* setup_getNonBond( int na, int nNode, Vec3i nPBC_, OCLtask* task=0){
        printf("setup_getNonBond(na=%i,nnode=%i) \n", na, nNode);
        if(task==0) task = getTask("getNonBond");
        //int nloc = 1;
        //int nloc = 2;
        //int nloc = 4;
        //int nloc = 8;
        //int nloc = 16;
        int nloc = 32;
        //int nloc = 64;
        task->local.x  = nloc;
        task->global.x = na + nloc-(na%nloc); // round up to multiple of nloc
        task->global.y = nSystems;

        useKernel( task->ikernel );
        nDOFs.x=na; 
        nDOFs.y=nNode;
        //nDOFs.x=bPBC; 
        v2i4( nPBC_, nPBC );
        // ------- Maybe We do-not need to do this every frame ?
        int err=0;
        err |= _useArg   ( nDOFs );               // 1 
        // Dynamical
        err |= useArgBuff( ibuff_atoms      ); // 2  
        err |= useArgBuff( ibuff_aforces    ); // 3
        // parameters
        err |= useArgBuff( ibuff_REQs      );  // 4
        err |= useArgBuff( ibuff_neighs    );  // 5
        err |= useArgBuff( ibuff_neighCell );  // 6
        err |= useArgBuff( ibuff_lvecs     );  // 7
        //err |= _useArg( cl_lvec          );  // 7
        err |= _useArg( nPBC               );  // 8
        err |= _useArg( GFFparams          );  // 9
        OCL_checkError(err, "setup_getNonBond");
        return task;
        // const int4 ns,                  // 1
        // // Dynamical
        // __global float4*  atoms,        // 2
        // __global float4*  forces,       // 3
        // // Parameters
        // __global float4*  REQKs,        // 4
        // __global int4*    neighs,       // 5
        // __global int4*    neighCell,    // 6
        // const int4 nPBC,                // 7
        // const cl_Mat3 lvec,             // 8
        // float R2damp                    // 9
    }

    OCLtask* setup_getNonBond_GridFF( int na, int nNode, Vec3i nPBC_, OCLtask* task=0){
        printf("setup_getNonBond_GridFF(na=%i,nnode=%i) itex_FE_Paul=%i itex_FE_Lond=%i itex_FE_Coul=%i\n", na, nNode,  itex_FE_Paul,itex_FE_Lond,itex_FE_Coul );
        if((itex_FE_Paul<0)||(itex_FE_Paul<=0)||(itex_FE_Paul<=0)){ printf( "ERROR in setup_getNonBond_GridFF() GridFF textures not initialized(itex_FE_Paul=%i itex_FE_Lond=%i itex_FE_Coul=%i) => Exit() \n", itex_FE_Paul,itex_FE_Lond,itex_FE_Coul ); exit(0); }
        if(task==0) task = getTask("getNonBond_GridFF");        
        //int nloc = 1;
        //int nloc = 4;
        //int nloc = 8;
        int nloc = 32;
        //int nloc = 64;
        task->local.x  = nloc;
        task->global.x = na + nloc-(na%nloc);
        task->global.y = nSystems;
        grid_shift0_p0 = grid_p0 + grid_shift0;

        useKernel( task->ikernel );
        nDOFs.x=na; 
        nDOFs.y=nNode; 

        v2i4( nPBC_, nPBC );
        // ------- Maybe We do-not need to do this every frame ?
        int err=0;
        err |= _useArg   ( nDOFs );           // 1
        // Dynamical
        err |= useArgBuff( ibuff_atoms     ); // 2
        err |= useArgBuff( ibuff_aforces   ); // 3
        // parameters
        err |= useArgBuff( ibuff_REQs      );  // 4
        err |= useArgBuff( ibuff_neighs    );  // 5
        err |= useArgBuff( ibuff_neighCell );  // 6
        err |= useArgBuff( ibuff_lvecs     );  // 7
        err |= _useArg( nPBC               );  // 8
        err |= _useArg( GFFparams          );  // 9
        err |= useArgBuff( itex_FE_Paul    );  // 10
        err |= useArgBuff( itex_FE_Lond    );  // 11
        err |= useArgBuff( itex_FE_Coul    );  // 12   
        //err |= _useArg( cl_diGrid          );  // 13
        err |= _useArg( cl_grid_ilvec      );  // 13
        err |= _useArg( grid_shift0_p0     );  // 14
        OCL_checkError(err, "setup_getNonBond_GridFF");
        return task;
        // const int4 ns,                  // 1
        // // Dynamical
        // __global float4*  atoms,        // 2
        // __global float4*  forces,       // 3
        // // Parameters
        // __global float4*  REQKs,        // 4
        // __global int4*    neighs,       // 5
        // __global int4*    neighCell,    // 6
        // __global cl_Mat3* lvecs,        // 7
        // const int4 nPBC,                // 8
        // const float4 GFFparams,         // 9
        // // GridFF
        // __read_only image3d_t  FE_Paul, // 10
        // __read_only image3d_t  FE_Lond, // 11
        // __read_only image3d_t  FE_Coul, // 12
        // const cl_Mat3  grid_invd,       // 13
        // const float4   grid_p0          // 14
    }

    OCLtask* setup_getMMFFf4( int na, int nNode, bool bPBC=false, OCLtask* task=0){
        printf("setup_getMMFFf4(na=%i,nnode=%i) \n", na, nNode);
        if(task==0) task = getTask("getMMFFf4");
        task->global.x = nNode;
        task->global.y = nSystems;
        useKernel( task->ikernel );
        nDOFs.x=na; 
        nDOFs.y=nNode; 
        nDOFs.w=bPBC; 
        // ------- Maybe We do-not need to do this every frame ?
        int err=0;
        err |= _useArg   ( nDOFs );            // 1
        // Dynamical
        err |= useArgBuff( ibuff_atoms  );     // 2
        err |= useArgBuff( ibuff_aforces);     // 3
        err |= useArgBuff( ibuff_neighForce ); // 4
        // parameters
        err |= useArgBuff( ibuff_neighs    );  // 5
        err |= useArgBuff( ibuff_neighCell );  // 5
        err |= useArgBuff( ibuff_REQs   );     // 6
        err |= useArgBuff( ibuff_MMpars );     // 7
        err |= useArgBuff( ibuff_BLs    );     // 8
        err |= useArgBuff( ibuff_BKs    );     // 9
        err |= useArgBuff( ibuff_Ksp    );     // 10
        err |= useArgBuff( ibuff_Kpp    );     // 11
        err |= useArgBuff( ibuff_lvecs  );     // 12
        err |= useArgBuff( ibuff_ilvecs );     // 13
        err |= useArgBuff( ibuff_pbcshifts );  // 13
        err |= _useArg   ( npbc         );  
        err |= _useArg   ( bSubtractVdW ); 
        

        //err |= _useArg( cl_lvec    );        // 12
        //err |= _useArg( cl_invLvec );        // 13
        OCL_checkError(err, "setup_getMMFFf4");
        return task;
        // const int4 nDOFs,              // 1   (nAtoms,nnode)
        // // Dynamical
        // __global float4*  apos,        // 2    [natoms]
        // __global float4*  fapos,       // 3    [natoms]     
        // __global float4*  fneigh,      // 4    [nnode*4]
        // // parameters
        // __global int4*    neighs,       // 5  [nnode]  neighboring atoms
        // __global float4*  REQKs,        // 6  [natoms] non-boding parametes {R0,E0,Q} 
        // __global float4*  apars,        // 7 [nnode]  per atom forcefield parametrs {c0ss,Kss,c0sp}
        // __global float4*  bLs,          // 8 [nnode]  bond lengths  for each neighbor
        // __global float4*  bKs,          // 9 [nnode]  bond stiffness for each neighbor
        // __global float4*  Ksp,          // 10 [nnode]  stiffness of pi-alignment for each neighbor
        // __global float4*  Kpp,          // 11 [nnode]  stiffness of pi-planarization for each neighbor
        // const cl_Mat3 lvec,             // 12
        // const cl_Mat3 invLvec           // 13
    }

    OCLtask* setup_updateAtomsMMFFf4( int na, int nNode,  OCLtask* task=0 ){
        printf( "setup_updateAtomsMMFFf4() \n" );
        if(task==0) task = getTask("updateAtomsMMFFf4");
        task->global.x = na+nNode;
        task->global.y = nSystems;
        //task->local .x = 1;
        //task->roundSizes();
        //if(n >=0  ) 
        nDOFs.x=na; 
        nDOFs.y=nNode; 
        int err=0;
        useKernel( task->ikernel  );
        err |= _useArg( nDOFs     );           // 1
        err |= useArgBuff( ibuff_atoms      ); // 2
        err |= useArgBuff( ibuff_avel       ); // 3
        err |= useArgBuff( ibuff_aforces    ); // 4
        err |= useArgBuff( ibuff_cvf        ); // 5
        err |= useArgBuff( ibuff_neighForce ); // 6
        err |= useArgBuff( ibuff_bkNeighs   ); // 7
        err |= useArgBuff( ibuff_constr     ); // 8
        err |= useArgBuff( ibuff_MDpars     ); // 9
        err |= useArgBuff( ibuff_TDrive     ); // 10
        OCL_checkError(err, "setup_updateAtomsMMFFf4");
        return task;
        // const int4        n,            // 1
        // __global float4*  apos,         // 2
        // __global float4*  avel,         // 3
        // __global float4*  aforce,       // 4
        // __global float4*  cvf,          // 5
        // __global float4*  fneigh,       // 6
        // __global int4*    bkNeighs,     // 7
        // __global float4*  constr,       // 8
        // __global float4*  MDparams      // 9
    }

    void printOnGPU( int isys, int4 mask=int4{1,1,1,1}, OCLtask* task=0 ){
        printf( "printOnGPU() \n" );
        if(task==0) task = getTask("printOnGPU");
        printf(  "task %li \n", task  );
        //task->global.x = na+nNode;
        //task->global.y = nSystems;
        task->global.x = 1;
        task->global.y = 1;                 
        nDOFs.x=nAtoms; 
        nDOFs.y=nnode;
        nDOFs.z=isys;                      
        print_mask=mask;
        int err=0;
        useKernel( task->ikernel  );       
        err |= _useArg( nDOFs     );       
        err |= _useArg( mask      );       
        OCL_checkError(err, "printOnGPU().setup");    
        err |= task->enque_raw();                      
        OCL_checkError(err, "printOnGPU().enque"  );  
        err |= finishRaw();                           
        OCL_checkError(err, "printOnGPU().finish" );  
    }


    OCLtask* setup_printOnGPU( int na, int nNode,  OCLtask* task=0 ){
        printf( "setup_printOnGPU() \n" );
        if(task==0) task = getTask("printOnGPU");
        //task->global.x = na+nNode;
        //task->global.y = nSystems;
        task->global.x = 1;
        task->global.y = 1;
        nDOFs.x=na; 
        nDOFs.y=nNode; 
        nDOFs.z=0; 
        print_mask=int4{1,1,1,1};
        int err=0;
        useKernel( task->ikernel  );
        err |= _useArg( nDOFs          );      // 1
        err |= _useArg( print_mask     );      // 2
        err |= useArgBuff( ibuff_atoms      ); // 3
        err |= useArgBuff( ibuff_avel       ); // 4
        err |= useArgBuff( ibuff_aforces    ); // 5
        err |= useArgBuff( ibuff_neighForce ); // 6
        err |= useArgBuff( ibuff_bkNeighs   ); // 7
        err |= useArgBuff( ibuff_constr     ); // 8
        OCL_checkError(err, "setup_printOnGPU");
        return task;
        // const int4        n,            // 1
        // __global float4*  apos,         // 2
        // __global float4*  avel,         // 3
        // __global float4*  aforce,       // 4
        // __global float4*  fneigh,       // 5
        // __global int4*    bkNeighs      // 6
    }

    OCLtask* setup_cleanForceMMFFf4( int na, int nNode,  OCLtask* task=0 ){
        printf( "setup_cleanForceMMFFf4() \n" );
        if(task==0) task = getTask("cleanForceMMFFf4");
        task->global.x = na+nNode;
        task->global.y = nSystems;
        nDOFs.x=na; 
        nDOFs.y=nNode;
        int err=0; 
        useKernel( task->ikernel );
        err |= _useArg( nDOFs     );           // 1
        err |= useArgBuff( ibuff_aforces    ); // 2
        err |= useArgBuff( ibuff_neighForce ); // 3
        OCL_checkError(err, "setup_cleanForceMMFFf4");
        return task;
        // const int4        n,           // 2
        // __global float4*  aforce,      // 5
        // __global float4*  fneigh       // 6
    }


    OCLtask* sampleGridFF( int n, Quat4f* fs=0, Quat4f* ps=0, Quat4f* REQs=0, bool bRun=true, OCLtask* task=0){
        //printf("OCL_MM::sampleGridFF() n=%i bRun=%i fs=%li ps=%li REQs=%li\n", n, bRun, fs,ps,REQs );
        int err=0;
        if((itex_FE_Paul<0)||(itex_FE_Lond<0 )||(itex_FE_Coul<0)){ printf("ERROR in OCL_MM::sampleGridFF() textures not initialized itex_FE_Paul=%i itex_FE_Lond=%i itex_FE_Coul=%i => Exit()\n",  itex_FE_Paul, itex_FE_Lond,  itex_FE_Coul ); }
        if(( ibuff_atoms<0)||(ibuff_aforces<0)||(ibuff_REQs<0  )){ printf("ERROR in OCL_MM::sampleGridFF() buffers  not initialized ibuff_atoms=%i ibuff_aforces=%i ibuff_REQs=%i => Exit()\n",    ibuff_atoms,  ibuff_aforces, ibuff_REQs   ); }
        // DEBUG
        // //if(ibuff_samp_fs <=0)ibuff_samp_fs   = newBuffer( "samp_fs",   n, sizeof(float4), 0, CL_MEM_WRITE_ONLY  );   DEBUG
        
        int nalloc = _max(n,1000);
        if(ibuff_samp_fs <0)ibuff_samp_fs   = newBuffer( "samp_fs",   nalloc, sizeof(float4), 0, CL_MEM_WRITE_ONLY );
        if(ibuff_samp_ps <0)ibuff_samp_ps   = newBuffer( "samp_ps",   nalloc, sizeof(float4), 0, CL_MEM_READ_ONLY  );
        if(ibuff_samp_REQ<0)ibuff_samp_REQ  = newBuffer( "samp_REQ",  nalloc, sizeof(float4), 0, CL_MEM_READ_ONLY  );
        if( buffers[ibuff_samp_ps].n < n ){ printf("ERROR in OCL_MM::sampleGridFF() buffer samp_ps.n(%i)<n(%i) => Exit() \n", buffers[ibuff_samp_ps].n, n ); exit(0); }
        //if(ibuff_atoms_surf<=0) ibuff_atoms_surf = newBuffer( "atoms_surf", na, sizeof(float4), 0, CL_MEM_READ_ONLY );
        //if(ibuff_REQs_surf <=0) ibuff_REQs_surf  = newBuffer( "REQs_surf",  na, sizeof(float4), 0, CL_MEM_READ_ONLY );
        if(ps   )upload( ibuff_samp_ps , ps,   n );
        if(REQs )upload( ibuff_samp_REQ, REQs, n );
        //if(ps   )upload( ibuff_atoms, ps,   n );   
        //if(REQs )upload( ibuff_REQs,  REQs, n );

        OCL_checkError(err, "sampleGridFF().upload" );
        if(task==0) task = getTask("sampleGridFF");
        task->local.x  = 1;
        task->global.x = n;
        useKernel( task->ikernel );
        grid_shift0_p0 = grid_p0 + grid_shift0;  
        nDOFs.x=n; 
        // ------- Maybe We do-not need to do this every frame ?
        err |= _useArg   ( nDOFs );            // 1
        err |= useArgBuff( ibuff_samp_ps   );  // 2
        err |= useArgBuff( ibuff_samp_fs   );  // 3
        err |= useArgBuff( ibuff_samp_REQ  );  // 4
        err |= _useArg   ( GFFparams       );  // 5
        err |= useArgBuff( itex_FE_Paul    );  // 6
        err |= useArgBuff( itex_FE_Lond    );  // 7
        err |= useArgBuff( itex_FE_Coul    );  // 8   
        //err |= _useArg( cl_diGrid        );  // 9    With un-Normalized Coordienates for texture sampler
        err |= _useArg( cl_grid_ilvec      );  // 9      With Normalized Coordienates for texture sampler
        err |= _useArg( grid_shift0_p0     );  // 10
        OCL_checkError(err, "sampleGridFF.setup");
        if(bRun){
            err |= task->enque_raw();                 OCL_checkError(err, "sampleGridFF().enque"    );
            err |= download( ibuff_samp_fs, fs, n );  OCL_checkError(err, "sampleGridFF().downalod" );
            err |= finishRaw();                       OCL_checkError(err, "sampleGridFF().finish"   );
        }
        //exit(0);
        return task;
        // const int4 ns,                  // 1
        // __global float4*  atoms,        // 2
        // __global float4*  forces,       // 3
        // __global float4*  REQKs,        // 4
        // const float4  GFFParams,        // 5
        // __read_only image3d_t  FE_Paul, // 6
        // __read_only image3d_t  FE_Lond, // 7
        // __read_only image3d_t  FE_Coul, // 8
        // const cl_Mat3  diGrid,          // 9
        // const float4   grid_p0          // 10
        
    }

    //void setGridShape( const Mat3d& dCell ){
    void setGridShape( const GridShape& grid ){
        v2i4      ( grid.n      , grid_n       );
        //v2f4      ( grid.pos0   , grid_p0      );
        grid_p0.f = (Vec3f)grid.pos0;
        Mat3_to_cl( grid.dCell  , cl_dGrid     );
        Mat3_to_cl( grid.diCell , cl_diGrid    );
        Mat3_to_cl( grid.cell   , cl_grid_lvec );
        Mat3_to_cl( grid.iCell  , cl_grid_ilvec );
    }

    OCLtask* getSurfMorse(  Vec3i nPBC_, int na=0, float4* atoms=0, float4* REQs=0, int na_s=0, float4* atoms_s=0, float4* REQs_s=0,  bool bRun=true, OCLtask* task=0 ){
        v2i4( nPBC_, nPBC );
        //if(ibuff_atoms_surf<0) ibuff_atoms_surf = newBuffer( "atoms_surf", na, sizeof(float4), 0, CL_MEM_READ_ONLY );
        //if(ibuff_REQs_surf <0) ibuff_REQs_surf  = newBuffer( "REQs_surf",  na, sizeof(float4), 0, CL_MEM_READ_ONLY );
        printf( "!!!!!!!!!! OCL_MM::getSurfMorse() ibuffs: atoms_surf(%i) REQs_surf(%i) atoms(%i) REQs(%i) aforces(%i) \n", ibuff_atoms_surf, ibuff_REQs_surf, ibuff_atoms, ibuff_REQs, ibuff_aforces );
        int err=0;
        err |= finishRaw();       OCL_checkError(err, "getSurfMorse().imgAlloc" );
        //OCLtask* task = tasks[ task_dict["getSurfMorse"] ];
        nDOFs.x = nAtoms;
        nDOFs.y = nnode;
        nDOFs.z = natom_surf; 
        if(task==0) task = getTask("getSurfMorse");
        //int nloc = 1;
        //int nloc = 4;
        //int nloc = 8;
        int nloc  = 32;
        //int nloc = 64;
        task->local.x = nloc;
        task->global.x = nAtoms + nloc-(nAtoms%nloc);
        //task->local.x = 1;
        task->local.y = 1;
        //task->global.x = nAtoms;
        task->global.y = nSystems;
        //if(atoms){ err |= upload( ibuff_atoms_surf, atoms, na ); OCL_checkError(err, "getSurfMorse().upload(atoms)" ); natom_surf = na; }
        //if(REQs ){ err |= upload( ibuff_REQs_surf , REQs , na ); OCL_checkError(err, "getSurfMorse().upload(REQs )" ); }
        useKernel( task->ikernel );

        err |= _useArg   ( nDOFs );            // 1
        err |= useArgBuff( ibuff_atoms      ); // 2
        err |= useArgBuff( ibuff_REQs       ); // 3
        err |= useArgBuff( ibuff_aforces    ); // 4
        //   Still good - no segfault
        err |= useArgBuff( ibuff_atoms_surf ); // 5
        err |= useArgBuff( ibuff_REQs_surf  ); // 6
        err |= _useArg( nPBC            );     // 7       
        err |= _useArg( cl_grid_lvec    );     // 8
        //err |= _useArg( grid_p0         );     // 9
        err |= _useArg( grid_shift0     );     // 9
        err |= _useArg( GFFparams       );     // 10

        // err |= _useArg   ( nDOFs );            OCL_checkError(err, "arg[1]: " );  // 1
        // err |= useArgBuff( ibuff_atoms      ); OCL_checkError(err, "arg[2]: " );// 2
        // err |= useArgBuff( ibuff_REQs       ); OCL_checkError(err, "arg[3]: " );// 3
        // err |= useArgBuff( ibuff_aforces    ); OCL_checkError(err, "arg[4]: " );// 4
        // err |= useArgBuff( ibuff_atoms_surf ); OCL_checkError(err, "arg[5]: " );// 5
        // err |= useArgBuff( ibuff_REQs_surf  ); OCL_checkError(err, "arg[6]: " );// 6
        // err |= _useArg( nPBC            );     OCL_checkError(err, "arg[7]: " );// 7       
        // err |= _useArg( cl_grid_lvec    );     OCL_checkError(err, "arg[8]: " );// 8
        // err |= _useArg( grid_p0         );     OCL_checkError(err, "arg[9]: " );// 9
        // err |= _useArg( GFFparams       );     OCL_checkError(err, "arg[10]: " );// 10
        OCL_checkError(err, "getSurfMorse().setup");
        if(bRun){
            err |= task->enque_raw(); OCL_checkError(err, "getSurfMorse().enque"  );
            err |= finishRaw();       OCL_checkError(err, "getSurfMorse().finish" );
        }
        return task;
        // const int4 ns,                // 1
        // __global float4*  atoms,      // 2
        // __global float4*  REQs,       // 3
        // __global float4*  forces,     // 4
        // __global float4*  atoms_s,    // 5
        // __global float4*  REQ_s,      // 6
        // const int4        nPBC,       // 7
        // const cl_Mat3     lvec,       // 8
        // const float4      pos0,       // 9
        // const float4      GFFParams   // 10
        exit(0);
    }

    OCLtask* makeGridFF( const GridShape& grid, Vec3i nPBC_, int na=0, float4* atoms=0, float4* REQs=0, bool bRun=true, OCLtask* task=0 ){
        setGridShape( grid );
        v2i4( nPBC_, nPBC );
        if(ibuff_atoms_surf<0) ibuff_atoms_surf = newBuffer( "atoms_surf", na, sizeof(float4), 0, CL_MEM_READ_ONLY );
        if(ibuff_REQs_surf <0) ibuff_REQs_surf  = newBuffer( "REQs_surf",  na, sizeof(float4), 0, CL_MEM_READ_ONLY );
        printf( "OCL_MM::makeGridFF() grid_n(%i,%i,%i)\n", grid_n.x,grid_n.y,grid_n.z );
        if(itex_FE_Paul<=0) itex_FE_Paul         = newBufferImage3D( "FEPaul", grid_n.x, grid_n.y, grid_n.z, sizeof(float)*4, 0, CL_MEM_READ_WRITE, {CL_RGBA, CL_FLOAT} );
        if(itex_FE_Lond<=0) itex_FE_Lond         = newBufferImage3D( "FFLond", grid_n.x, grid_n.y, grid_n.z, sizeof(float)*4, 0, CL_MEM_READ_WRITE, {CL_RGBA, CL_FLOAT} );
        if(itex_FE_Coul<=0) itex_FE_Coul         = newBufferImage3D( "FFCoul", grid_n.x, grid_n.y, grid_n.z, sizeof(float)*4, 0, CL_MEM_READ_WRITE, {CL_RGBA, CL_FLOAT} );
        int err=0;
        err |= finishRaw();       OCL_checkError(err, "makeGridFF().imgAlloc" );
        //OCLtask* task = tasks[ task_dict["make_GridFF"] ];
        if(task==0) task = getTask("make_GridFF");

        //int nloc = 1;
        //int nloc = 4;
        //int nloc = 8;
        int nloc  = 32;
        //int nloc = 64;
        int ngrid = grid_n.x*grid_n.y*grid_n.z;
        task->local.x  = nloc;
        task->global.x = ngrid + nloc-(ngrid%nloc);

        //printf( "makeGridFF() na=%i nG=%i(%i,%i,%i) nPBC(%i,%i,%i) \n", na, task->global.x, grid_n.x,grid_n.y,grid_n.z,  nPBC.x,nPBC.y,nPBC.z );
        //printf("ibuff_atoms_surf %li, ibuff_REQs_surf %li \n", ibuff_atoms_surf, ibuff_REQs_surf );
        if(atoms){ err |= upload( ibuff_atoms_surf, atoms, na ); OCL_checkError(err, "makeGridFF().upload(atoms)" ); natom_surf = na; }
        if(REQs ){ err |= upload( ibuff_REQs_surf , REQs , na ); OCL_checkError(err, "makeGridFF().upload(REQs )" ); }
        useKernel( task->ikernel );
        err |= useArg    ( natom_surf       ); // 1
        err |= useArgBuff( ibuff_atoms_surf ); // 2
        err |= useArgBuff( ibuff_REQs_surf  ); // 3
        err |= useArgBuff( itex_FE_Paul );     // 4
        err |= useArgBuff( itex_FE_Lond );     // 5
        err |= useArgBuff( itex_FE_Coul );     // 6
        err |= _useArg( nPBC            );     // 7     
        err |= _useArg( grid_n          );     // 8      
        err |= _useArg( cl_grid_lvec    );     // 9
        err |= _useArg( grid_p0         );     // 10
        err |= _useArg( GFFparams       );     // 11
        OCL_checkError(err, "makeGridFF().setup");
        if(bRun){
            err |= task->enque_raw(); OCL_checkError(err, "makeGridFF().enque"  );
            err |= finishRaw();       OCL_checkError(err, "makeGridFF().finish" );
        }
        return task;
        // const int nAtoms,                // 1
        // __global float4*  atoms,         // 2
        // __global float4*  REQKs,         // 3
        // __write_only image3d_t  FE_Paul, // 4
        // __write_only image3d_t  FE_Lond, // 5
        // __write_only image3d_t  FE_Coul, // 6
        // const int4     nPBC,             // 7
        // const int4     nGrid,            // 8
        // const cl_Mat3  lvec,             // 9
        // const float4   grid_p0           // 10

    }

    OCLtask* addDipoleField( const GridShape& grid, Vec3i nPBC_, int n=0, float4* dipole_ps=0, float4* dipoles=0, bool bRun=true, OCLtask* task=0 ){
        setGridShape( grid );
        v2i4( nPBC_, nPBC );
        if(ibuff_dipole_ps<0) ibuff_dipole_ps = newBuffer( "dipole_ps", n, sizeof(float4), 0, CL_MEM_READ_ONLY );
        if(ibuff_dipoles  <0) ibuff_dipoles   = newBuffer( "dipoles",   n, sizeof(float4), 0, CL_MEM_READ_ONLY );        
        //OCLtask* task = tasks[ task_dict["make_GridFF"] ];
        if(task==0) task = getTask("addDipoleField");
        task->global.x = grid_n.x*grid_n.y*grid_n.z;
        if(dipole_ps)upload( ibuff_dipole_ps, dipole_ps, n );
        if(dipoles  )upload( ibuff_dipoles  , dipoles  , n );
        useKernel( task->ikernel );
        int4 ngrid{ grid_n.x, grid_n.y, grid_n.z,0 };
        int err=0;
        err |= useArg    ( n               ); // 1
        err |= useArgBuff( ibuff_dipole_ps ); // 2
        err |= useArgBuff( ibuff_dipoles   ); // 3
        err |= useArgBuff( itex_FE_Coul    ); // 6   
        err |= _useArg( ngrid        );        // 8      
        err |= _useArg( grid_p0      );        // 9
        err |= _useArg( cl_dGrid     );        // 10
        OCL_checkError(err, "addDipoleField().setup");
        if(bRun){
            err |= task->enque_raw(); OCL_checkError(err, "addDipoleField().enque"  );
            err |= finishRaw();       OCL_checkError(err, "addDipoleField().finish" );
        }
        return task;
        // const int n,                     // 1
        // __global float4*  ps,            // 2
        // __global float4*  dipols,        // 3
        // __write_only image3d_t  FE_Coul, // 4
        // const int4     nGrid,            // 5
        // const cl_Mat3  dGrid,            // 6
        // const float4   grid_p0           // 7
    }


    OCLtask* PPAFM_makeFF( int isys, const GridShape& grid, Vec3i nPBC_, bool bRun=true, OCLtask* task=0 ){
        setGridShape( grid );
        v2i4      ( nPBC_      , afm_nPBC      );
        v2i4      ( grid.n     , afm_grid_n    );
        Mat3_to_cl( grid.cell  , afm_grid_lvec );
        Mat3_to_cl( grid.iCell , afm_grid_ilvec );
        afm_grid_p0.f = (Vec3f)grid.pos0;
        printf( "OCL_MM::PPAFM_makeFF() grid_n(%i,%i,%i)\n",       afm_grid_n.x, afm_grid_n.y, afm_grid_n.z );
        if(itex_FEAFM<=0) itex_FEAFM  = newBufferImage3D( "FEAFM", afm_grid_n.x, afm_grid_n.y, afm_grid_n.z, sizeof(float)*4, 0, CL_MEM_READ_WRITE, {CL_RGBA, CL_FLOAT} );
        printf( "OCL_MM::PPAFM_makeFF() ibuff_atoms=%li ibuff_REQs=%li itex_FEAFM=%li \n", ibuff_atoms, ibuff_REQs, itex_FEAFM );
        int err=0;
        err |= finishRaw();       OCL_checkError(err, "PPAFM_makeFF().alloc(FEAFM)" );
        //OCLtask* task = tasks[ task_dict["make_GridFF"] ];
        if(task==0) task = getTask("PPAFM_makeFF");

        //int nloc = 1;
        //int nloc = 4;
        //int nloc = 8;
        int nloc  = 32;
        //int nloc = 64;
        int ngrid      = afm_grid_n.x*afm_grid_n.y*afm_grid_n.z;
        task->local.x  = nloc;
        task->global.x = ngrid + nloc-(ngrid%nloc);

        //printf( "makeGridFF() na=%i nG=%i(%i,%i,%i) nPBC(%i,%i,%i) \n", na, task->global.x, grid_n.x,grid_n.y,grid_n.z,  nPBC.x,nPBC.y,nPBC.z );
        //printf("ibuff_atoms_surf %li, ibuff_REQs_surf %li \n", ibuff_atoms_surf, ibuff_REQs_surf );
        //if(atoms){ err |= upload( ibuff_atoms_surf, atoms, na ); OCL_checkError(err, "PPAFM_makeFF().upload(atoms)" ); natom_surf = na; }
        //if(REQs ){ err |= upload( ibuff_REQs_surf , REQs , na ); OCL_checkError(err, "PPAFM_makeFF().upload(REQs )" ); }
   
        // int4      afm_grid_n;
        // float4    afm_grid_p0;         
        // cl_Mat3   afm_grid_lvec;
        // //float4  tipRot[3]={{1.,0.,0.,0.},{0.,1.,0.,0.},{0.,0.,1.,0.1}};  // tip rotation
        // cl_Mat3   tipRot;          
        // float4  afm_relax_params{0.5f,0.1f,0.02f,0.5f};
        // float4  tipParams{  1.661f, 0.0091063f, -0.1f, 0.0f };
        // float4  tip_stiffness { -0.03f,-0.03f,-0.03f,-1.00 };  // tip stiffness
        // float4  tip_dpos0{0.0f,0.0f,-4.0f, 4.0f};    
        // float4  tip_Qs {0.f,-0.05f,0.f,0.0f};
        // float4  tip_QZs{0.1f, 0.0f,-0.1f,0.0f};
        // float4  afm_surfFF{0.f,0.f,0.f,0.f};
        // int     afm_nz,afm_nzout; 



        useKernel( task->ikernel );
        err |= _useArg   ( nDOFs        );   // 1
        err |= useArgBuff( ibuff_atoms  );   // 2
        err |= useArgBuff( ibuff_REQs   );   // 3
        err |= useArgBuff( itex_FEAFM   );   // 4
        err |= _useArg( afm_nPBC        );   // 5     
        err |= _useArg( afm_grid_n      );   // 6      
        err |= _useArg( afm_grid_lvec   );   // 7
        err |= _useArg( afm_grid_p0     );   // 8
        err |= _useArg( tipParams       );   // 9
        err |= _useArg( tip_Qs          );   // 10
        err |= _useArg( tip_QZs         );   // 11
        OCL_checkError(err, "PPAFM_makeFF().setup");
        if(bRun){
            err |= task->enque_raw(); OCL_checkError(err, "PPAFM_makeFF().enque"  );
            err |= finishRaw();       OCL_checkError(err, "PPAFM_makeFF().finish" );
        }
        return task;
        // __kernel void PPAFM_makeFF(
        // const int nAtoms,                // 1
        // __global float4*  atoms,         // 2
        // __global float4*  REQs,          // 3
        // __write_only image3d_t  imgOut,  // 4
        // const int4     nPBC,             // 5
        // const int4     nGrid,            // 6
        // const cl_Mat3  dlvec,            // 7
        // const float4   grid_p0,          // 8
        // const float4 tipParams,          // 9
        // const float4 Qs,                 // 10
        // const float4 QZs                 // 11
    }

    OCLtask* PPAFM_scan( int np, int nz, Quat4f* ps, Quat4f* FEout, Quat4f* PPpos, float dTip=0.1, Mat3d tipRot_=Mat3dIdentity, int nMaxItr=128, bool bRun=true, OCLtask* task=0 ){
        //if(ibuff_atoms_surf<=0) ibuff_atoms_surf = newBuffer( "atoms_surf", na, sizeof(float4), 0, CL_MEM_READ_ONLY );
        //if(ibuff_REQs_surf <=0) ibuff_REQs_surf  = newBuffer( "REQs_surf",  na, sizeof(float4), 0, CL_MEM_READ_ONLY );
        printf( "OCL_MM::PPAFM_scan() np=%i nz=%i nw=%i \n", np, nz );
        if(ibuff_afm_ps   <=0) ibuff_afm_ps    = newBuffer( "afm_ps",    np, sizeof(float4),    0, CL_MEM_READ_ONLY  );
        if(ibuff_afm_FEout<=0) ibuff_afm_FEout = newBuffer( "afm_FEout", np*nz, sizeof(float4), 0, CL_MEM_WRITE_ONLY );
        if(ibuff_afm_PPpos<=0) ibuff_afm_PPpos = newBuffer( "afm_PPpos", np*nz, sizeof(float4), 0, CL_MEM_WRITE_ONLY );
        int err=0;
        err |= finishRaw();       OCL_checkError(err, "PPAFM_makeFF().imgAlloc" );;
        if(task==0) task = getTask("PPAFM_scan");
        task->global.x = np;
        if(ps){ err |= upload( ibuff_afm_ps, ps, np ); OCL_checkError(err, "PPAFM_scan().upload(ps)" ); }
        afm_nz = nz;

        Mat3_to_cl( tipRot_, tipRot );
        tipRot.c.s[3]=-dTip;
        afm_nMaxItr=nMaxItr;

        printf( "OCL_MM::PPAFM_scan() tipRot{{%g,%g,%g,%g},{%g,%g,%g,%g},{%g,%g,%g,%g}}\n", tipRot.a.s[0],tipRot.a.s[1],tipRot.a.s[2],tipRot.a.s[3],   tipRot.b.s[0],tipRot.b.s[1],tipRot.b.s[2],tipRot.b.s[3],   tipRot.c.s[0],tipRot.c.s[1],tipRot.c.s[2],tipRot.c.s[3] );

        // int4      afm_grid_n;
        // float4    afm_grid_p0;         
        // cl_Mat3   afm_grid_lvec;
        // //float4  tipRot[3]={{1.,0.,0.,0.},{0.,1.,0.,0.},{0.,0.,1.,0.1}};  // tip rotation
        // cl_Mat3   tipRot;          
        // float4  afm_relax_params{0.5f,0.1f,0.02f,0.5f};
        // float4  tipParams{  1.661f, 0.0091063f, -0.1f, 0.0f };
        // float4  tip_stiffness { -0.03f,-0.03f,-0.03f,-1.00 };  // tip stiffness
        // float4  tip_dpos0{0.0f,0.0f,-4.0f, 4.0f};    
        // float4  tip_Qs {0.f,-0.05f,0.f,0.0f};
        // float4  tip_QZs{0.1f, 0.0f,-0.1f,0.0f};
        // float4  afm_surfFF{0.f,0.f,0.f,0.f};
        // int     afm_nz,afm_nzout; 

        useKernel( task->ikernel );
        err |= useArgBuff( itex_FEAFM      ); // 1
        err |= useArgBuff( ibuff_afm_ps    ); // 2
        err |= useArgBuff( ibuff_afm_FEout ); // 3
        err |= useArgBuff( ibuff_afm_PPpos ); // 4      
        err |= _useArg( afm_grid_ilvec     ); // 5
        err |= _useArg( tipRot             ); // 6
        err |= _useArg( tip_stiffness      ); // 7
        err |= _useArg( tip_dpos0          ); // 8
        err |= _useArg( afm_relax_params   ); // 9
        err |= _useArg( afm_nz             ); // 10
        err |= _useArg( afm_nMaxItr        ); // 11
        OCL_checkError(err, "PPAFM_scan().setup");
        if(bRun){
            err |= task->enque_raw();                            OCL_checkError(err, "PPAFM_scan().enque"    );
            if(FEout)err |= download( ibuff_afm_FEout, FEout );  OCL_checkError(err, "PPAFM_scan().download(FEout)" );
            if(PPpos)err |= download( ibuff_afm_PPpos, PPpos );  OCL_checkError(err, "PPAFM_scan().download(PPpos)" );
            err |= finishRaw();                                  OCL_checkError(err, "PPAFM_scan().finish"   );
            //if(FEout){ err |= upload( ibuff_afm_ws, ws, na ); OCL_checkError(err, "PPAFM_scan().upload(ws)" ); }
        }
        return task;
        // __read_only image3d_t  imgIn,   // 1 
        // __global  float4*      points,  // 2
        // __global  float4*      FEs,     // 3
        // __global  float4*      PPpos,   // 4
        // const cl_Mat3  diGrid,          // 5
        // const cl_Mat3  tipRot,          // 6
        // float4 stiffness,               // 7
        // float4 dpos0,                   // 8
        // float4 relax_params,            // 9
        // const int nz,                   // 10
        // const int nMaxItr               // 11
    }


    OCLtask* PPAFM_scan_df( int np, int nz, int nw, Quat4f* ps, float* ws, Quat4f* FEout, bool bRun=true, OCLtask* task=0 ){
        //if(ibuff_atoms_surf<=0) ibuff_atoms_surf = newBuffer( "atoms_surf", na, sizeof(float4), 0, CL_MEM_READ_ONLY );
        //if(ibuff_REQs_surf <=0) ibuff_REQs_surf  = newBuffer( "REQs_surf",  na, sizeof(float4), 0, CL_MEM_READ_ONLY );
        printf( "OCL_MM::PPAFM_scan_df() np=%i nz=%i nw=%i \n", np, nz, nw );
        if(ibuff_afm_ps   <=0) ibuff_afm_ps    = newBuffer( "afm_ps",    np, sizeof(float4),    0, CL_MEM_READ_ONLY  );
        if(ibuff_afm_ws   <=0) ibuff_afm_ws    = newBuffer( "afm_ws",    nw, sizeof(float4),    0, CL_MEM_READ_ONLY  );
        if(ibuff_afm_FEout<=0) ibuff_afm_FEout = newBuffer( "afm_FEout", np*nz, sizeof(float4), 0, CL_MEM_WRITE_ONLY );
        int err=0;
        err |= finishRaw();       OCL_checkError(err, "PPAFM_scan_df().imgAlloc" );
        //OCLtask* task = tasks[ task_dict["make_GridFF"] ];
        if(task==0) task = getTask("PPAFM_scan_df");
        task->global.x = np;
        //printf( "makeGridFF() na=%i nG=%i(%i,%i,%i) nPBC(%i,%i,%i) \n", na, task->global.x, grid_n.x,grid_n.y,grid_n.z,  nPBC.x,nPBC.y,nPBC.z );
        //printf("ibuff_atoms_surf %li, ibuff_REQs_surf %li \n", ibuff_atoms_surf, ibuff_REQs_surf );
        if(ps){ err |= upload( ibuff_afm_ps, ps, np ); OCL_checkError(err, "PPAFM_scan().upload(ps)" ); }
        if(ws){ err |= upload( ibuff_afm_ws, ws, nw ); OCL_checkError(err, "PPAFM_scan().upload(ws)" ); }
        afm_nz = nz;
        useKernel( task->ikernel );
        err |= useArg ( itex_FEAFM      ); // 1
        err |= useArg ( ibuff_afm_ps    ); // 2
        err |= useArg ( ibuff_afm_ws    ); // 3
        err |= useArg ( ibuff_afm_FEout ); // 4
        err |= _useArg( nPBC            ); // 7       
        err |= _useArg( afm_grid_lvec   ); // 9
        err |= _useArg( tipRot          ); // 10
        err |= _useArg( tip_stiffness   ); // 11
        err |= _useArg( tip_dpos0       ); // 11
        err |= _useArg( afm_surfFF      ); // 11
        err |= _useArg( afm_nz              ); // 11
        err |= _useArg( afm_nzout           ); // 11
        OCL_checkError(err, "PPAFM_scan().setup");
        if(bRun){
            err |= task->enque_raw();                   OCL_checkError(err, "PPAFM_scan_df().enque"    );
            err |= download( ibuff_afm_FEout, FEout );  OCL_checkError(err, "PPAFM_scan_df().download" );
            err |= finishRaw();                         OCL_checkError(err, "PPAFM_scan_df().finish"   );
            //if(FEout){ err |= upload( ibuff_afm_ws, ws, na ); OCL_checkError(err, "PPAFM_scan().upload(ws)" ); }
        }
        return task;
        // __kernel void PPAFM_scan(
        //     __read_only image3d_t  imgIn,   // 1 
        //     __global    float4*    points,  // 2
        //     __constant  float*     weighs,  // 3
        //     __global    float4*    FEs,     // 4
        //     const cl_Mat3  diGrid,          // 5
        //     const cl_Mat3  tipRot,          // 6
        //     float4 stiffness,               // 7
        //     float4 dpos0,                   // 8
        //     float4 relax_params,            // 9
        //     float4 surfFF,                  // 11
        //     const int nz,                   // 12
        //     const int nzout                 // 13
    }




    OCLtask* setup_evalMMFFf4_local1( int niter_, OCLtask* task=0 ){
        printf("OCL_MM::setup_evalMMFFf4_local1()\n");
        if(task==0) task = getTask("evalMMFFf4_local1");
        
        //md_params.y = 0.9;
        //int nloc = 1;
        //int nloc = 4;
        //int nloc = 8;
        //int nloc  = 32;
        //int nloc = 64;
        int nloc  = nAtoms;
        task->local.x  = nloc;
        task->global.x = nnode + nloc-(nnode%nloc);
        task->global.y = nSystems;
        useKernel( task->ikernel );
        niter   = niter_;
        nDOFs.x = nAtoms; 
        nDOFs.y = nnode; 
        // ------- Maybe We do-not need to do this every frame ?
        int err=0;
        err |= _useArg   ( nDOFs           );    OCL_checkError(err, "setup_evalMMFFf4_local.1");   DEBUG
        err |= useArgBuff( ibuff_atoms     );    OCL_checkError(err, "setup_evalMMFFf4_local.2");   DEBUG
        err |= useArgBuff( ibuff_avel      );    OCL_checkError(err, "setup_evalMMFFf4_local.3");   DEBUG
        err |= useArgBuff( ibuff_constr    );    OCL_checkError(err, "setup_evalMMFFf4_local.4");   DEBUG
        err |= useArgBuff( ibuff_neighs    );    OCL_checkError(err, "setup_evalMMFFf4_local.5");   DEBUG
        err |= useArgBuff( ibuff_neighCell );    OCL_checkError(err, "setup_evalMMFFf4_local.6");   DEBUG
        //err |= useArgBuff( ibuff_bkNeighs  );    OCL_checkError(err, "setup_evalMMFFf4_local.7");   DEBUG
        err |= useArgBuff( ibuff_bkNeighs_new  );    OCL_checkError(err, "setup_evalMMFFf4_local.7");   DEBUG
        err |= useArgBuff( ibuff_REQs      );    OCL_checkError(err, "setup_evalMMFFf4_local.8");   DEBUG
        err |= useArgBuff( ibuff_MMpars    );    OCL_checkError(err, "setup_evalMMFFf4_local.9");   DEBUG
        err |= useArgBuff( ibuff_BLs       );    OCL_checkError(err, "setup_evalMMFFf4_local.10");  DEBUG
        err |= useArgBuff( ibuff_BKs       );    OCL_checkError(err, "setup_evalMMFFf4_local.11");  DEBUG
        err |= useArgBuff( ibuff_Ksp       );    OCL_checkError(err, "setup_evalMMFFf4_local.12");  DEBUG
        err |= useArgBuff( ibuff_Kpp       );    OCL_checkError(err, "setup_evalMMFFf4_local.13");  DEBUG
        err |= useArgBuff( ibuff_lvecs     );    OCL_checkError(err, "setup_evalMMFFf4_local.14");  DEBUG
        err |= useArgBuff( ibuff_ilvecs    );    OCL_checkError(err, "setup_evalMMFFf4_local.15");  DEBUG
        err |= _useArg   ( nPBC            );    OCL_checkError(err, "setup_evalMMFFf4_local.16");  DEBUG
        err |= _useArg   ( GFFparams       );    OCL_checkError(err, "setup_evalMMFFf4_local.17");  DEBUG
        err |= _useArg   ( md_params       );    OCL_checkError(err, "setup_evalMMFFf4_local.18");  DEBUG
        err |= _useArg   ( niter           );    OCL_checkError(err, "setup_evalMMFFf4_local.19");  DEBUG
        OCL_checkError(err, "setup_evalMMFFf4_local1"); DEBUG;
        return task;
        // const int4 nDOFs,               // 1  (nAtoms,nnode)
        // __global float4*  apos,         // 2  [natoms]
        // __global float4*  avel,         // 3
        // __global float4*  constr,       // 4
        // __global int4*    neighs,       // 5  [nnode]  neighboring atoms
        // __global int4*    neighCell,    // 6
        // __global int4*    bkneighs,     // 7
        // __global float4*  REQs,         // 8  [natoms] non-boding parametes {R0,E0,Q} 
        // __global float4*  apars,        // 9  [nnode]  per atom forcefield parametrs {c0ss,Kss,c0sp}
        // __global float4*  bLs,          // 10 [nnode]  bond lengths  for each neighbor
        // __global float4*  bKs,          // 11 [nnode]  bond stiffness for each neighbor
        // __global float4*  Ksp,          // 12 [nnode]  stiffness of pi-alignment for each neighbor
        // __global float4*  Kpp,          // 13 [nnode]  stiffness of pi-planarization for each neighbor
        // __global cl_Mat3* lvecs,        // 14
        // __global cl_Mat3* ilvecs,       // 15
        // const int4        nPBC,         // 16
        // const float4      GFFparams,    // 17
        // const float4      MDpars,       // 18
        // const int         niter         // 19
    }

    OCLtask* setup_evalMMFFf4_local2( int niter_, OCLtask* task=0 ){
        printf("OCL_MM::setup_evalMMFFf4_local2()\n");
        if(task==0) task = getTask("evalMMFFf4_local2");
        
        //md_params.y = 0.9;
        //int nloc = 1;
        //int nloc = 4;
        //int nloc = 8;
        //int nloc  = 32;
        //int nloc = 64;
        int nloc  = _max( nnode, ncap );
        task->local.x  = nloc;
        task->global.x = nnode + nloc-(nnode%nloc);
        task->global.y = nSystems;
        useKernel( task->ikernel );
        niter   = niter_;
        nDOFs.x = nAtoms; 
        nDOFs.y = nnode; 
        // ------- Maybe We do-not need to do this every frame ?
        int err=0;
        err |= _useArg   ( nDOFs           );    OCL_checkError(err, "setup_evalMMFFf4_local.1");   DEBUG
        err |= useArgBuff( ibuff_atoms     );    OCL_checkError(err, "setup_evalMMFFf4_local.2");   DEBUG
        err |= useArgBuff( ibuff_avel      );    OCL_checkError(err, "setup_evalMMFFf4_local.3");   DEBUG
        err |= useArgBuff( ibuff_constr    );    OCL_checkError(err, "setup_evalMMFFf4_local.4");   DEBUG
        err |= useArgBuff( ibuff_neighs    );    OCL_checkError(err, "setup_evalMMFFf4_local.5");   DEBUG
        err |= useArgBuff( ibuff_neighCell );    OCL_checkError(err, "setup_evalMMFFf4_local.6");   DEBUG
        //err |= useArgBuff( ibuff_bkNeighs  );    OCL_checkError(err, "setup_evalMMFFf4_local.7");   DEBUG
        err |= useArgBuff( ibuff_bkNeighs_new  );    OCL_checkError(err, "setup_evalMMFFf4_local.7");   DEBUG
        err |= useArgBuff( ibuff_REQs      );    OCL_checkError(err, "setup_evalMMFFf4_local.8");   DEBUG
        err |= useArgBuff( ibuff_MMpars    );    OCL_checkError(err, "setup_evalMMFFf4_local.9");   DEBUG
        err |= useArgBuff( ibuff_BLs       );    OCL_checkError(err, "setup_evalMMFFf4_local.10");  DEBUG
        err |= useArgBuff( ibuff_BKs       );    OCL_checkError(err, "setup_evalMMFFf4_local.11");  DEBUG
        err |= useArgBuff( ibuff_Ksp       );    OCL_checkError(err, "setup_evalMMFFf4_local.12");  DEBUG
        err |= useArgBuff( ibuff_Kpp       );    OCL_checkError(err, "setup_evalMMFFf4_local.13");  DEBUG
        err |= useArgBuff( ibuff_lvecs     );    OCL_checkError(err, "setup_evalMMFFf4_local.14");  DEBUG
        err |= useArgBuff( ibuff_ilvecs    );    OCL_checkError(err, "setup_evalMMFFf4_local.15");  DEBUG
        err |= _useArg   ( nPBC            );    OCL_checkError(err, "setup_evalMMFFf4_local.16");  DEBUG
        err |= _useArg   ( GFFparams       );    OCL_checkError(err, "setup_evalMMFFf4_local.17");  DEBUG
        err |= _useArg   ( md_params       );    OCL_checkError(err, "setup_evalMMFFf4_local.18");  DEBUG
        err |= _useArg   ( niter           );    OCL_checkError(err, "setup_evalMMFFf4_local.19");  DEBUG
        OCL_checkError(err, "setup_evalMMFFf4_local2"); DEBUG;
        return task;
        // const int4 nDOFs,               // 1  (nAtoms,nnode)
        // __global float4*  apos,         // 2  [natoms]
        // __global float4*  avel,         // 3
        // __global float4*  constr,       // 4
        // __global int4*    neighs,       // 5  [nnode]  neighboring atoms
        // __global int4*    neighCell,    // 6
        // __global int4*    bkneighs,     // 7
        // __global float4*  REQs,         // 8  [natoms] non-boding parametes {R0,E0,Q} 
        // __global float4*  apars,        // 9  [nnode]  per atom forcefield parametrs {c0ss,Kss,c0sp}
        // __global float4*  bLs,          // 10 [nnode]  bond lengths  for each neighbor
        // __global float4*  bKs,          // 11 [nnode]  bond stiffness for each neighbor
        // __global float4*  Ksp,          // 12 [nnode]  stiffness of pi-alignment for each neighbor
        // __global float4*  Kpp,          // 13 [nnode]  stiffness of pi-planarization for each neighbor
        // __global cl_Mat3* lvecs,        // 14
        // __global cl_Mat3* ilvecs,       // 15
        // const int4        nPBC,         // 16
        // const float4      GFFparams,    // 17
        // const float4      MDpars,       // 18
        // const int         niter         // 19
    }

   OCLtask* setup_evalMMFFf4_local_test( int niter_, OCLtask* task=0 ){
        
        if(task==0) task = getTask("evalMMFFf4_local_test");
        //md_params.y = 0.9;
        //int nloc = 1;
        //int nloc = 4;
        //int nloc = 8;
        //int nloc  = 32;
        //int nloc = 64;
        int nloc  = _max( nnode, ncap );
        task->local.x  = nloc;
        task->global.x = nnode + nloc-(nnode%nloc);
        task->global.y = nSystems;
        useKernel( task->ikernel );
        niter   = niter_;
        nDOFs.x = nAtoms; 
        nDOFs.y = nnode;  
        printf("OCL_MM::evalMMFFf4_local_test(niter_=%i) nAtoms=%i nnode=%i nglob(%i,%i) nloc(%i,%i) \n", niter_,    nAtoms,nnode,     task->global.x,task->global.y,   task->local.x,task->local.y );
        DEBUG
        // ------- Maybe We do-not need to do this every frame ?
        int err=0;
        err |= _useArg   ( nDOFs           );    OCL_checkError(err, "setup_evalMMFFf4_local_test.1");   DEBUG
        err |= useArgBuff( ibuff_atoms     );    OCL_checkError(err, "setup_evalMMFFf4_local_test.2");   DEBUG
        err |= useArgBuff( ibuff_avel      );    OCL_checkError(err, "setup_evalMMFFf4_local_test.3");   DEBUG
        err |= useArgBuff( ibuff_constr    );    OCL_checkError(err, "setup_evalMMFFf4_local_test.4");   DEBUG
        err |= useArgBuff( ibuff_neighs    );    OCL_checkError(err, "setup_evalMMFFf4_local_test.5");   DEBUG
        err |= useArgBuff( ibuff_neighCell );    OCL_checkError(err, "setup_evalMMFFf4_local_test.6");   DEBUG
        //err |= useArgBuff( ibuff_bkNeighs  );    OCL_checkError(err, "setup_evalMMFFf4_local_test.7");   DEBUG
        err |= useArgBuff( ibuff_bkNeighs_new  );    OCL_checkError(err, "setup_evalMMFFf4_local_test.7");   DEBUG
        err |= useArgBuff( ibuff_REQs      );    OCL_checkError(err, "setup_evalMMFFf4_local_test.8");   DEBUG
        err |= useArgBuff( ibuff_MMpars    );    OCL_checkError(err, "setup_evalMMFFf4_local_test.9");   DEBUG
        err |= useArgBuff( ibuff_BLs       );    OCL_checkError(err, "setup_evalMMFFf4_local_test.10");  DEBUG
        err |= useArgBuff( ibuff_BKs       );    OCL_checkError(err, "setup_evalMMFFf4_local_test.11");  DEBUG
        err |= useArgBuff( ibuff_Ksp       );    OCL_checkError(err, "setup_evalMMFFf4_local_test.12");  DEBUG
        err |= useArgBuff( ibuff_Kpp       );    OCL_checkError(err, "setup_evalMMFFf4_local_test.13");  DEBUG
        err |= useArgBuff( ibuff_lvecs     );    OCL_checkError(err, "setup_evalMMFFf4_local_test.14");  DEBUG
        err |= useArgBuff( ibuff_ilvecs    );    OCL_checkError(err, "setup_evalMMFFf4_local_test.15");  DEBUG
        err |= _useArg   ( nPBC            );    OCL_checkError(err, "setup_evalMMFFf4_local_test.16");  DEBUG
        err |= _useArg   ( GFFparams       );    OCL_checkError(err, "setup_evalMMFFf4_local_test.17");  DEBUG
        err |= _useArg   ( md_params       );    OCL_checkError(err, "setup_evalMMFFf4_local_test.18");  DEBUG
        err |= _useArg   ( niter           );    OCL_checkError(err, "setup_evalMMFFf4_local_test.19");  DEBUG
        OCL_checkError(err, "setup_evalMMFFf4_local_test"); DEBUG;
        return task;
        // const int4 nDOFs,               // 1  (nAtoms,nnode)
        // __global float4*  apos,         // 2  [natoms]
        // __global float4*  avel,         // 3
        // __global float4*  constr,       // 4
        // __global int4*    neighs,       // 5  [nnode]  neighboring atoms
        // __global int4*    neighCell,    // 6
        // __global int4*    bkneighs,     // 7
        // __global float4*  REQs,         // 8  [natoms] non-boding parametes {R0,E0,Q} 
        // __global float4*  apars,        // 9  [nnode]  per atom forcefield parametrs {c0ss,Kss,c0sp}
        // __global float4*  bLs,          // 10 [nnode]  bond lengths  for each neighbor
        // __global float4*  bKs,          // 11 [nnode]  bond stiffness for each neighbor
        // __global float4*  Ksp,          // 12 [nnode]  stiffness of pi-alignment for each neighbor
        // __global float4*  Kpp,          // 13 [nnode]  stiffness of pi-planarization for each neighbor
        // __global cl_Mat3* lvecs,        // 14
        // __global cl_Mat3* ilvecs,       // 15
        // const int4        nPBC,         // 16
        // const float4      GFFparams,    // 17
        // const float4      MDpars,       // 18
        // const int         niter         // 19
    }

};

#endif
