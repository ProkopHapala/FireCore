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
//class OCL_MM: public OCL_DFT { public:
class OCL_MM: public OCLsystem { public:
    cl_program program_relax=0;

    int nAtoms=0;
    int nnode=0, nvecs=0, nneigh=0, npi=0, nSystems=0;

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

    void makeKrenels_MM( const char*  cl_src_dir ){
        printf( "makeKrenels_MM() \n" );
        char srcpath[1024];
        sprintf( srcpath, "%s/relax_multi.cl", cl_src_dir );     
        buildProgram( srcpath, program_relax );
        newTask( "getNonBond"             ,program_relax, 2);
        newTask( "getMMFFf4"              ,program_relax, 2);
        newTask( "cleanForceMMFFf4"       ,program_relax, 2);
        newTask( "updateAtomsMMFFf4"      ,program_relax, 2);
        //newTask( "write_toImg"     ,program_relax, 3,{0,0,0,0},{1,1,1,0} ); 
        printf( "... makeKrenels_MM() DONE \n" );
    }

    int initAtomsForces( int nSystems_, int nAtoms_, int nnode_ ){
        nSystems=nSystems_;
        nnode  = nnode_;
        nAtoms = nAtoms_;
        npi    = nnode_;
        nvecs  = nAtoms+npi;
        nneigh = nnode*4*2;
        printf( "initAtomsForces() nSystems*nvecs %i nSystems %i nvecs %i natoms %i nnode %i npi %i \n", nSystems*nvecs, nSystems, nvecs, nAtoms, nnode, npi );
        ibuff_atoms      = newBuffer( "atoms",      nSystems*nvecs , sizeof(float4), 0, CL_MEM_READ_WRITE );
        ibuff_aforces    = newBuffer( "aforces",    nSystems*nvecs , sizeof(float4), 0, CL_MEM_READ_WRITE );
        ibuff_REQs       = newBuffer( "REQs",       nSystems*nAtoms, sizeof(float4), 0, CL_MEM_READ_ONLY  );
        ibuff_neighs     = newBuffer( "neighs",     nSystems*nAtoms, sizeof(int4  ), 0, CL_MEM_READ_ONLY  );
        ibuff_neighCell  = newBuffer( "neighCell" , nSystems*nAtoms, sizeof(int4  ), 0, CL_MEM_READ_ONLY  );

        ibuff_bkNeighs   = newBuffer( "bkNeighs",   nSystems*nvecs,  sizeof(int4  ), 0, CL_MEM_READ_ONLY  );
        ibuff_avel       = newBuffer( "avel",       nSystems*nvecs,  sizeof(float4), 0, CL_MEM_READ_WRITE );
        ibuff_neighForce = newBuffer( "neighForce", nSystems*nneigh, sizeof(float4), 0, CL_MEM_READ_WRITE );

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
        nPBC.x = nPBC_.x;
        nPBC.y = nPBC_.y;
        nPBC.z = nPBC_.z;
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
        //err |= _useArg( cl_lvec    );          // 12
        //err |= _useArg( cl_invLvec );          // 13
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
        /*
            const int4        n,           // 2
            __global float4*  aforce,      // 5
            __global float4*  fneigh       // 6
        */
    }

};

#endif
