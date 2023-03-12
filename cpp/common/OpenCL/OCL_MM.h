#ifndef  OCL_MM_h
#define  OCL_MM_h

#include "OCL_DFT.h"

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
class OCL_MM: public OCL_DFT { public:
    cl_program program_relax=0;

    int4   nDOFs    {0,0,0,0};
    int4   nPBC     {0,0,0,0};
    float4 md_params{0.05,0.9,100.0,0.0};    // (dt,cdamp,forceLimit)
    //float Rdamp  = 1e-4;
    float Rdamp  = 1.;
    //float R2damp = Rdamp*Rdamp;

    cl_Mat3 cl_lvec;
    cl_Mat3 cl_invLvec;

    int nz;

    int n_start_point = 0;
    int ibuff_start_point=-1;
    int ibuff_avel=-1, ibuff_pi0s=-1, ibuff_neighForce=-1,  ibuff_bkNeighs=-1;
    int ibuff_bondLK=-1, ibuff_ang0K=-1;
    int itex_FF=-1;
    int ibuff_MMpars=-1, ibuff_BLs=-1,ibuff_BKs=-1,ibuff_Ksp=-1, ibuff_Kpp=-1;   // MMFFf4 params

    int itex_FE_Paul=-1;
    int itex_FE_Lond=-1;
    int itex_FE_Coul=-1;

    void makeKrenels_MM( const char*  cl_src_dir ){
        printf( "makeKrenels_MM() \n" );
        char srcpath[1024];
        sprintf( srcpath, "%s/relax_multi.cl", cl_src_dir );     
        buildProgram( srcpath, program_relax );
        newTask( "getNonBond"             ,program_relax);
        newTask( "getMMFFf4"              ,program_relax);
        newTask( "cleanForceMMFFf4"       ,program_relax);
        newTask( "updateAtomsMMFFf4"      ,program_relax);
        //newTask( "write_toImg"     ,program_relax, 3,{0,0,0,0},{1,1,1,0} ); 
        printf( "... makeKrenels_MM() DONE \n" );
    }

    int initMM( const char*  cl_src_dir ){
        makeKrenels_MM( cl_src_dir );
        printf( "initMM() Ns(%li,%li,%li) \n", Ns[0], Ns[1], Ns[2] );
        return itex_FF;
    }

    int initAtomsForces( int nSystems, int nAtoms_, int nnode, bool bMMFFsp3=false, bool bMMFFf4=false ){
        nAtoms=nAtoms_;
        int npi=nnode;
        int nvecs=nAtoms+npi;
        int nneigh=nnode*4*2;
        ibuff_atoms      = newBuffer( "atoms",      nSystems*nvecs , sizeof(float4), 0, CL_MEM_READ_WRITE );
        ibuff_aforces    = newBuffer( "aforces",    nSystems*nvecs , sizeof(float4), 0, CL_MEM_READ_WRITE );
        ibuff_coefs      = newBuffer( "coefs",      nSystems*nAtoms, sizeof(float4), 0, CL_MEM_READ_ONLY  );
        ibuff_neighs     = newBuffer( "neighs",     nSystems*nAtoms, sizeof(int4  ), 0, CL_MEM_READ_ONLY  ); // need neihgs for all atoms because of Non-Bonded
        ibuff_neighCell  = newBuffer( "neighCell" , nSystems*nAtoms, sizeof(int4  ), 0, CL_MEM_READ_ONLY  );

        ibuff_bkNeighs   = newBuffer( "bkNeighs",   nSystems*nvecs,  sizeof(int4  ), 0, CL_MEM_READ_ONLY  );
        ibuff_avel       = newBuffer( "avel",       nSystems*nvecs,  sizeof(float4), 0, CL_MEM_READ_WRITE );
        ibuff_neighForce = newBuffer( "neighForce", nSystems*nneigh, sizeof(float4), 0, CL_MEM_READ_WRITE );

        ibuff_MMpars     = newBuffer( "MMpars",     nSystems*nnode,  sizeof(int4),   0, CL_MEM_READ_ONLY );
        ibuff_BLs        = newBuffer( "BLs",        nSystems*nnode,  sizeof(float4), 0, CL_MEM_READ_ONLY  );
        ibuff_BKs        = newBuffer( "BKs",        nSystems*nnode,  sizeof(float4), 0, CL_MEM_READ_ONLY  );
        ibuff_Ksp        = newBuffer( "Ksp",        nSystems*nnode,  sizeof(float4), 0, CL_MEM_READ_ONLY );
        ibuff_Kpp        = newBuffer( "Kpp",        nSystems*nnode,  sizeof(float4), 0, CL_MEM_READ_ONLY );

        return ibuff_atoms;
    }


    OCLtask* setup_getNonBond( int na, int nNode, Vec3i nPBC_, float Rdamp_, OCLtask* task=0){
        printf("!!!!! setup_getNonBond(na=%i,nnode=%i) \n", na, nNode);
        if(task==0) task = getTask("getNonBond");
        task->global.x = na;
        useKernel( task->ikernel );
        nDOFs.x=na; 
        nDOFs.y=nNode; 
        //nDOFs.x=bPBC; 
        Rdamp = Rdamp_;
        nPBC.x = nPBC_.x;
        nPBC.y = nPBC_.y;
        nPBC.z = nPBC_.z;
        // ------- Maybe We do-not need to do this every frame ?
        err |= _useArg   ( nDOFs );               // 1
        // Dynamical
        err |= useArgBuff( ibuff_atoms      ); // 2
        err |= useArgBuff( ibuff_aforces    ); // 3
        // parameters
        err |= useArgBuff( ibuff_coefs     );  // 4
        err |= useArgBuff( ibuff_neighs    );  // 5
        err |= useArgBuff( ibuff_neighCell );  // 6
        err |= _useArg( nPBC               );  // 7
        err |= _useArg( cl_lvec            );  // 8
        err |= _useArg( Rdamp              );  // 9
        OCL_checkError(err, "setup_getNonBond");
        return task;
        /*
        __kernel void getNonBond(
            const int4 ns,                  // 1
            // Dynamical
            __global float4*  atoms,        // 2
            __global float4*  forces,       // 3
            // Parameters
            __global float4*  REQKs,        // 4
            __global int4*    neighs,       // 5
            __global int4*    neighCell,    // 6
            const int4 nPBC,                // 7
            const cl_Mat3 lvec,             // 8
            float R2damp                    // 9
        ){
        */
    }

    OCLtask* setup_getMMFFf4( int na, int nNode, bool bPBC=false, OCLtask* task=0){
        //printf("setup_getMMFFsp3(na=%i,nnode=%i) \n", na, nNode);
        if(task==0) task = getTask("getMMFFf4");
        task->global.x = nNode;
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
        err |= useArgBuff( ibuff_coefs  );     // 6
        err |= useArgBuff( ibuff_MMpars );     // 7
        err |= useArgBuff( ibuff_BLs    );     // 8
        err |= useArgBuff( ibuff_BKs    );     // 9
        err |= useArgBuff( ibuff_Ksp    );     // 10
        err |= useArgBuff( ibuff_Kpp    );     // 11
        err |= _useArg( cl_lvec    );          // 12
        err |= _useArg( cl_invLvec );          // 13
        OCL_checkError(err, "setup_getMMFFf4");
        return task;
        /*
        __kernel void getMMFFf4(
            const int4 nDOFs,              // 1   (nAtoms,nnode)
            // Dynamical
            __global float4*  apos,        // 2    [natoms]
            __global float4*  fapos,       // 3    [natoms]     
            __global float4*  fneigh,      // 4    [nnode*4]
            // parameters
            __global int4*    neighs,       // 5  [nnode]  neighboring atoms
            __global float4*  REQKs,        // 6  [natoms] non-boding parametes {R0,E0,Q} 
            __global float4*  apars,        // 7 [nnode]  per atom forcefield parametrs {c0ss,Kss,c0sp}
            __global float4*  bLs,          // 8 [nnode]  bond lengths  for each neighbor
            __global float4*  bKs,          // 9 [nnode]  bond stiffness for each neighbor
            __global float4*  Ksp,          // 10 [nnode]  stiffness of pi-alignment for each neighbor
            __global float4*  Kpp,          // 11 [nnode]  stiffness of pi-planarization for each neighbor
            const cl_Mat3 lvec,             // 12
            const cl_Mat3 invLvec           // 13
        ){
    */
    }

    OCLtask* setup_updateAtomsMMFFf4( int na, int nNode,  OCLtask* task=0 ){
        if(task==0) task = getTask("updateAtomsMMFFf4");
        task->global.x = na+nNode;
        //task->local .x = 1;
        //task->roundSizes();
        //if(n >=0  ) 
        nDOFs.x=na; 
        nDOFs.y=nNode; 
        useKernel( task->ikernel );
        err |= _useArg( md_params );           // 1
        err |= _useArg( nDOFs     );           // 2
        err |= useArgBuff( ibuff_atoms      ); // 3
        err |= useArgBuff( ibuff_avel       ); // 4
        err |= useArgBuff( ibuff_aforces    ); // 5
        err |= useArgBuff( ibuff_neighForce ); // 6
        err |= useArgBuff( ibuff_bkNeighs   ); // 7
        OCL_checkError(err, "setup_updateAtomsMMFFf4");
        return task;
        /*
            const float4      MDpars,       // 1
            const int4        n,            // 2
            __global float4*  apos,         // 3
            __global float4*  avel,         // 4
            __global float4*  aforce,       // 5
            __global float4*  fneigh,       // 6
            __global int4*    bkNeighs      // 7
        */
    }

    OCLtask* setup_cleanForceMMFFf4( int na, int nNode,  OCLtask* task=0 ){
        if(task==0) task = getTask("cleanForceMMFFf4");
        task->global.x = na+nNode;
        nDOFs.x=na; 
        nDOFs.y=nNode; 
        useKernel( task->ikernel );
        err |= _useArg( nDOFs     );           // 1
        err |= useArgBuff( ibuff_aforces    ); // 2
        err |= useArgBuff( ibuff_neighForce ); // 3
        OCL_checkError(err, "setup_cleanForceMMFFf4");
        return task;
        /*
            const int4        n,           // 2
            __global float4*  aforce,      // 5
            __global float4*  fneigh       // 6
        */
    }

};

#endif