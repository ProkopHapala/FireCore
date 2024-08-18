#ifndef  OCL_EFF_h
#define  OCL_EFF_h

#include "OCL.h"

class OCL_EFF : public OCLsystem { public:
    cl_program program_eff=0;

    int nSystems=0, nAtom=0, nElec=0;
    int nAtomTot=0, nElecTot=0;

    int
    ibuff_apos    =-1,
    ibuff_aforce  =-1,
    ibuff_avel    =-1,
    ibuff_aParams =-1,

    ibuff_epos    =-1,
    ibuff_eforce  =-1,
    ibuff_evel    =-1,
    ibuff_espin   =-1;

    float4 KRSrho{ 1.125, 0.9, -0.2, 1.0 };  // 
    

    // ====================== Functions

    void makeKrenels_EFF( const char*  cl_src_dir ){
        printf( "makeKrenels_EFF() \n" );
        char srcpath[1024];
        sprintf( srcpath, "%s/relax_multi.cl", cl_src_dir );     
        buildProgram( srcpath, program_eff );
        newTask( "eval_electrons"         ,program_eff, 2);
        newTask( "eval_ions"              ,program_eff, 2);
        printf( "... makeKrenels_EFF() DONE \n" );
    }

    void initBuffersEFF( int nSystems_, int nAtom_, int nElec_ ){

        nSystems=nSystems_;
        nAtom   = nAtom_;
        nElec   = nElec_;
        nAtomTot = nSystems*nAtom;
        nElecTot = nSystems*nElec;

        ibuff_apos      = newBuffer( "atoms",      nAtomTot , sizeof(float4), 0, CL_MEM_READ_WRITE );
        ibuff_aforce    = newBuffer( "aforces",    nAtomTot , sizeof(float4), 0, CL_MEM_READ_WRITE );
        ibuff_avel      = newBuffer( "avel",       nAtomTot , sizeof(float4), 0, CL_MEM_READ_WRITE );
        ibuff_aParams   = newBuffer( "aParams",    nAtomTot , sizeof(float4), 0, CL_MEM_READ_WRITE );

        ibuff_epos      = newBuffer( "electrons",  nElecTot , sizeof(float4), 0, CL_MEM_READ_WRITE );
        ibuff_eforce    = newBuffer( "eforces",    nElecTot , sizeof(float4), 0, CL_MEM_READ_WRITE );
        ibuff_evel      = newBuffer( "evel",       nElecTot , sizeof(float4), 0, CL_MEM_READ_WRITE );
        ibuff_espin     = newBuffer( "espin",      nElecTot , sizeof(float4), 0, CL_MEM_READ_WRITE );

    }


    OCLtask* setup_eval_electrons(OCLtask* task=0){
        printf("setup_eval_electrons() \n" );
        if(task==0) task = getTask("eval_electrons");
        int nloc = 32;
        //int nloc = 64;
        task->local.x  = nloc;
        task->global.x = nElec;   // round up to multiple of nloc
        task->global.y = nSystems;

        useKernel( task->ikernel );
        int err=0;
        err |= _useArg   ( nAtom    );       // 1 
        err |= _useArg   ( nElec    );       // 2 
        err |= useArgBuff( ibuff_apos    );  // 3  
        err |= useArgBuff( ibuff_aforce  );  // 4
        err |= useArgBuff( ibuff_aParams );  // 5
        err |= useArgBuff( ibuff_epos    );  // 6
        err |= useArgBuff( ibuff_eforce  );  // 7
        err |= useArgBuff( ibuff_espin   );  // 8
        err |= _useArg( KRSrho           );  // 9
        OCL_checkError(err, "setup_eval_electrons");
        return task;

    }
/*
__kernel void eval_electrons(
    int na,                    // 1 
    int ne,                    // 2
    __global float4*  apos,    // 3      
    __global float4*  aforce,  // 4  
    __global float4*  aParams, // 5 
    __global float4*  epos,    // 6     
    __global float4*  eforce,  // 7
    __global int*     espin,   // 8
    float4 KRSrho              // 9
){
*/


    OCLtask* setup_eval_eval_ions(OCLtask* task=0){
        printf("setup_eval_ions() \n" );
        if(task==0) task = getTask("eval_ions");
        int nloc = 32;
        //int nloc = 64;
        task->local.x  = nloc;
        task->global.x = nAtom;   // round up to multiple of nloc
        task->global.y = nSystems;

        useKernel( task->ikernel );
        int err=0;
        err |= _useArg   ( nAtom    );       // 1 
        err |= _useArg   ( nElec    );       // 2 
        err |= useArgBuff( ibuff_apos    );  // 3  
        err |= useArgBuff( ibuff_aforce  );  // 4
        err |= useArgBuff( ibuff_aParams );  // 5
        err |= useArgBuff( ibuff_epos    );  // 6
        err |= useArgBuff( ibuff_eforce  );  // 7
        err |= _useArg( KRSrho           );  // 8
        OCL_checkError(err, "setup_eval_ions");
        return task;

    }
/*
__kernel void eval_ions(
    int na,                    // 1
    int ne,                    // 2
    __global float4*  apos,    // 3       
    __global float4*  aforce,  // 4    
    __global float4*  aParams, // 5  
    __global float4*  epos,    // 6      
    __global float4*  eforce,  // 7 
     float4 KRSrho             // 8
){
*/

};

#endif
