#ifndef  OCL_MM_h
#define  OCL_MM_h

//#include "OCL_DFT.h"
#include "OCL.h"

void v2f4( const Vec3d& v, float4& f4 ){ f4.x=(float)v.x; f4.y=(float)v.y; f4.z=(float)v.z; };
void v2i4( const Vec3i& v, int4& f4   ){ f4.x=(float)v.x; f4.y=(float)v.y; f4.z=(float)v.z; };
//void v2f4( const Vec3d& v, cl_float4& f4 ){ f4.s[0]=(cl_float)v.x; f4.s[1]=(cl_float)v.y; f4.s[2]=(cl_float)v.z; };
cl_float4 cl_f4( const Vec3d& v ){ return (cl_float4){(cl_float)v.x,(cl_float)v.y,(cl_float)v.z,0.f}; };

// see https://stackoverflow.com/questions/33119233/opencl-using-struct-as-kernel-argument
typedef struct __attribute__ ((packed)) cl_Mat3{
    cl_float4 a;
    cl_float4 b;
    cl_float4 c;
}st_foo;


void Mat3_to_cl( const Mat3d& m, cl_Mat3& clm ){ 
    //v2f4( m.a, clm.a ); v2f4( m.b, clm.b ); v2f4( m.c, clm.c ); 
    clm.a = cl_f4( m.a );
    clm.b = cl_f4( m.b );
    clm.c = cl_f4( m.c );
}

//=======================================================================
//=======================================================================
//class OCL_MM: public OCL_DFT { public:
class OCL_MM: public OCLsystem { public:
    cl_program program_relax=0;

    int nAtoms=0;
    int nnode=0, nvecs=0, nneigh=0, npi=0, nSystems=0, nbkng=0;

    int4   print_mask{1,1,1,1};
    int4   nDOFs    {0,0,0,0};
    int4   nPBC     {0,0,0,0};
    float4 md_params{0.05,0.9,100.0,0.0};    // (dt,cdamp,forceLimit)
    //float Rdamp  = 1e-4;
    float Rdamp  = 1.;
    //float R2damp = Rdamp*Rdamp;

    cl_Mat3 cl_lvec;
    cl_Mat3 cl_invLvec;

    int ibuff_atoms=-1,ibuff_aforces=-1,ibuff_neighs=-1,ibuff_neighCell=-1;
    int ibuff_avel=-1, ibuff_neighForce=-1,  ibuff_bkNeighs=-1;
    int ibuff_REQs=-1, ibuff_MMpars=-1, ibuff_BLs=-1,ibuff_BKs=-1,ibuff_Ksp=-1, ibuff_Kpp=-1;   // MMFFf4 params
    int ibuff_lvecs=-1, ibuff_ilvecs=-1; 
    int ibuff_constr=-1;

    // ------- Grid
    //size_t Ns[4]; // = {N0, N1, N2};
    //size_t Ntot;
    int natom_surf=0;
    int4    grid_n;
    float4  grid_p0{0.f,0.f,0.f,0.f};
    cl_Mat3 cl_dGrid;
    cl_Mat3 cl_diGrid;
    cl_Mat3 cl_grid_lvec;


    int itex_FE_Paul=-1;
    int itex_FE_Lond=-1;
    int itex_FE_Coul=-1;
    int ibuff_atoms_surf=-1;
    int ibuff_REQs_surf =-1;
    int ibuff_dipole_ps=-1;
    int ibuff_dipoles  =-1;

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
        newTask( "make_GridFF"            ,program_relax, 1);
        newTask( "addDipoleField"         ,program_relax, 1);
        //newTask( "write_toImg"     ,program_relax, 3,{0,0,0,0},{1,1,1,0} ); 
        printf( "... makeKrenels_MM() DONE \n" );
    }

    int initAtomsForces( int nSystems_, int nAtoms_, int nnode_ ){
        nSystems=nSystems_;
        nnode  = nnode_;
        nAtoms = nAtoms_;
        npi    = nnode_;
        nvecs  = nAtoms+npi;
        nbkng  = nnode*4*2;
        printf( "initAtomsForces() nSystems %i nvecs %i natoms %i nnode %i nbkng %i \n", nSystems, nvecs, nAtoms, nnode, nbkng );
        printf( "initAtomsForces() nS*nvecs %i nS*natoms %i nS*nnode %i nS*nbkng %i \n", nSystems*nvecs,  nSystems*nAtoms, nSystems*nnode, nSystems*nbkng );
        ibuff_atoms      = newBuffer( "atoms",      nSystems*nvecs , sizeof(float4), 0, CL_MEM_READ_WRITE );
        ibuff_aforces    = newBuffer( "aforces",    nSystems*nvecs , sizeof(float4), 0, CL_MEM_READ_WRITE );
        ibuff_REQs       = newBuffer( "REQs",       nSystems*nAtoms, sizeof(float4), 0, CL_MEM_READ_ONLY  );
        ibuff_neighs     = newBuffer( "neighs",     nSystems*nAtoms, sizeof(int4  ), 0, CL_MEM_READ_ONLY  );
        ibuff_neighCell  = newBuffer( "neighCell" , nSystems*nAtoms, sizeof(int4  ), 0, CL_MEM_READ_ONLY  );

        ibuff_constr     = newBuffer( "constr",    nSystems*nAtoms , sizeof(float4), 0, CL_MEM_READ_WRITE );
        //ibuff_constr0    = newBuffer( "constr0",   nSystems*nAtoms , sizeof(float4), 0, CL_MEM_READ_WRITE );
        //ibuff_constrK    = newBuffer( "constrK",   nSystems*nAtoms , sizeof(float4), 0, CL_MEM_READ_WRITE );

        ibuff_bkNeighs   = newBuffer( "bkNeighs",   nSystems*nvecs , sizeof(int4  ), 0, CL_MEM_READ_ONLY  );
        //ibuff_bkNeighs = newBuffer( "bkNeighs",   nSystems*nAtoms, sizeof(int4  ), 0, CL_MEM_READ_ONLY  );
        ibuff_avel       = newBuffer( "avel",       nSystems*nvecs,  sizeof(float4), 0, CL_MEM_READ_WRITE );
        ibuff_neighForce = newBuffer( "neighForce", nSystems*nbkng,  sizeof(float4), 0, CL_MEM_READ_WRITE );

        ibuff_MMpars     = newBuffer( "MMpars",     nSystems*nnode,  sizeof(int4),   0, CL_MEM_READ_ONLY  );
        ibuff_BLs        = newBuffer( "BLs",        nSystems*nnode,  sizeof(float4), 0, CL_MEM_READ_ONLY  );
        ibuff_BKs        = newBuffer( "BKs",        nSystems*nnode,  sizeof(float4), 0, CL_MEM_READ_ONLY  );
        ibuff_Ksp        = newBuffer( "Ksp",        nSystems*nnode,  sizeof(float4), 0, CL_MEM_READ_ONLY  );
        ibuff_Kpp        = newBuffer( "Kpp",        nSystems*nnode,  sizeof(float4), 0, CL_MEM_READ_ONLY  );

        ibuff_lvecs      = newBuffer( "lvecs",      nSystems,        sizeof(cl_Mat3), 0, CL_MEM_READ_ONLY  );
        ibuff_ilvecs     = newBuffer( "ilvecs",     nSystems,        sizeof(cl_Mat3), 0, CL_MEM_READ_ONLY  );

        return ibuff_atoms;
    }

    OCLtask* setup_getNonBond( int na, int nNode, Vec3i nPBC_, float Rdamp_, OCLtask* task=0){
        printf("setup_getNonBond(na=%i,nnode=%i) \n", na, nNode);
        if(task==0) task = getTask("getNonBond");
        task->global.x = na;
        task->global.y = nSystems;
        useKernel( task->ikernel );
        nDOFs.x=na; 
        nDOFs.y=nNode; 
        //nDOFs.x=bPBC; 
        Rdamp = Rdamp_;
        v2i4( nPBC_, nPBC );
        // ------- Maybe We do-not need to do this every frame ?
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
        err |= _useArg( Rdamp              );  // 9
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

    OCLtask* setup_getNonBond_GridFF( int na, int nNode, Vec3i nPBC_, float Rdamp_, OCLtask* task=0){
        printf("setup_getNonBond_GridFF(na=%i,nnode=%i) \n", na, nNode);
        if(task==0) task = getTask("getNonBond_GridFF");
        task->global.x = na;
        task->global.y = nSystems;
        useKernel( task->ikernel );
        nDOFs.x=na; 
        nDOFs.y=nNode; 
        //nDOFs.x=bPBC; 
        Rdamp = Rdamp_;
        v2i4( nPBC_, nPBC );
        // ------- Maybe We do-not need to do this every frame ?
        err |= _useArg   ( nDOFs );               // 1
        // Dynamical
        err |= useArgBuff( ibuff_atoms     ); // 2
        err |= useArgBuff( ibuff_aforces   ); // 3
        // parameters
        err |= useArgBuff( ibuff_REQs      );  // 4
        err |= useArgBuff( ibuff_neighs    );  // 5
        err |= useArgBuff( ibuff_neighCell );  // 6
        err |= useArgBuff( ibuff_lvecs     );  // 7
        err |= _useArg( nPBC               );  // 8
        err |= _useArg( Rdamp              );  // 9
        err |= useArgBuff( itex_FE_Paul    );  // 10
        err |= useArgBuff( itex_FE_Lond    );  // 11
        err |= useArgBuff( itex_FE_Coul    );  // 12    
        err |= _useArg( cl_diGrid          );  // 13
        err |= _useArg( grid_p0            );  // 14
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
        // const float Rdamp,              // 9
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
        err |= _useArg   ( nDOFs );            // 1
        // Dynamical
        err |= useArgBuff( ibuff_atoms  );     // 2
        err |= useArgBuff( ibuff_aforces);     // 3
        err |= useArgBuff( ibuff_neighForce ); // 4
        // parameters
        err |= useArgBuff( ibuff_neighs );     // 5
        err |= useArgBuff( ibuff_REQs   );     // 6
        err |= useArgBuff( ibuff_MMpars );     // 7
        err |= useArgBuff( ibuff_BLs    );     // 8
        err |= useArgBuff( ibuff_BKs    );     // 9
        err |= useArgBuff( ibuff_Ksp    );     // 10
        err |= useArgBuff( ibuff_Kpp    );     // 11
        err |= useArgBuff( ibuff_lvecs  );     // 12
        err |= useArgBuff( ibuff_ilvecs );     // 13
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
        useKernel( task->ikernel  );
        err |= _useArg( md_params );           // 1
        err |= _useArg( nDOFs     );           // 2
        err |= useArgBuff( ibuff_atoms      ); // 3
        err |= useArgBuff( ibuff_avel       ); // 4
        err |= useArgBuff( ibuff_aforces    ); // 5
        err |= useArgBuff( ibuff_neighForce ); // 6
        err |= useArgBuff( ibuff_bkNeighs   ); // 7
        err |= useArgBuff( ibuff_constr     ); // 8
        OCL_checkError(err, "setup_updateAtomsMMFFf4");
        return task;
        // const float4      MDpars,       // 1
        // const int4        n,            // 2
        // __global float4*  apos,         // 3
        // __global float4*  avel,         // 4
        // __global float4*  aforce,       // 5
        // __global float4*  fneigh,       // 6
        // __global int4*    bkNeighs      // 7
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

    //void setGridShape( const Mat3d& dCell ){
    void setGridShape( const GridShape& grid ){
        v2i4      ( grid.n      , grid_n       );
        v2f4      ( grid.pos0   , grid_p0      );
        Mat3_to_cl( grid.dCell  , cl_dGrid     );
        Mat3_to_cl( grid.diCell , cl_diGrid    );
        Mat3_to_cl( grid.cell   , cl_grid_lvec );
    }

    OCLtask* makeGridFF( const GridShape& grid, Vec3i nPBC_, int na=0, float4* atoms=0, float4* REQs=0, bool bRun=true, OCLtask* task=0 ){
        setGridShape( grid );
        v2i4( nPBC_, nPBC );
        if(ibuff_atoms_surf<=0) ibuff_atoms_surf = newBuffer( "atoms_surf", na, sizeof(float4), 0, CL_MEM_READ_ONLY );
        if(ibuff_REQs_surf <=0) ibuff_REQs_surf  = newBuffer( "REQs_surf",  na, sizeof(float4), 0, CL_MEM_READ_ONLY );
        if(itex_FE_Paul<=0) itex_FE_Paul         = newBufferImage3D( "FEPaul", grid_n.x, grid_n.y, grid_n.z, sizeof(float)*4, 0, CL_MEM_READ_WRITE, {CL_RGBA, CL_FLOAT} );
        if(itex_FE_Lond<=0) itex_FE_Lond         = newBufferImage3D( "FFLond", grid_n.x, grid_n.y, grid_n.z, sizeof(float)*4, 0, CL_MEM_READ_WRITE, {CL_RGBA, CL_FLOAT} );
        if(itex_FE_Coul<=0) itex_FE_Coul         = newBufferImage3D( "FFCoul", grid_n.x, grid_n.y, grid_n.z, sizeof(float)*4, 0, CL_MEM_READ_WRITE, {CL_RGBA, CL_FLOAT} );
        err |= finishRaw();       OCL_checkError(err, "makeGridFF().imgAlloc" );
        //OCLtask* task = tasks[ task_dict["make_GridFF"] ];
        if(task==0) task = getTask("make_GridFF");
        task->global.x = grid_n.x*grid_n.y*grid_n.z;
        //printf( "makeGridFF() na=%i nG=%i(%i,%i,%i) nPBC(%i,%i,%i) \n", na, task->global.x, grid_n.x,grid_n.y,grid_n.z,  nPBC.x,nPBC.y,nPBC.z );
        //printf("ibuff_atoms_surf %li, ibuff_REQs_surf %li \n", ibuff_atoms_surf, ibuff_REQs_surf );
        if(atoms){ err = upload( ibuff_atoms_surf, atoms, na ); OCL_checkError(err, "makeGridFF().upload(atoms)" ); natom_surf = na; }
        if(REQs ){ err = upload( ibuff_REQs_surf , REQs , na ); OCL_checkError(err, "makeGridFF().upload(REQs )" ); }
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
        if(ibuff_dipole_ps<=0) ibuff_dipole_ps = newBuffer( "dipole_ps", n, sizeof(float4), 0, CL_MEM_READ_ONLY );
        if(ibuff_dipoles  <=0) ibuff_dipoles   = newBuffer( "dipoles",   n, sizeof(float4), 0, CL_MEM_READ_ONLY );        
        //OCLtask* task = tasks[ task_dict["make_GridFF"] ];
        if(task==0) task = getTask("addDipoleField");
        task->global.x = grid_n.x*grid_n.y*grid_n.z;
        if(dipole_ps)upload( ibuff_dipole_ps, dipole_ps, n );
        if(dipoles  )upload( ibuff_dipoles  , dipoles  , n );
        useKernel( task->ikernel );
        int4 ngrid{ grid_n.x, grid_n.y, grid_n.z,0 };
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

};

#endif
