#ifndef  OCL_PP_h
#define  OCL_PP_h

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
class OCL_PP: public OCL_DFT { public:
    //DEFAULT_dTip         = np.array( [ 0.0 , 0.0 , -0.1 , 0.0 ], dtype=np.float32 );
    //DEFAULT_stiffness    = np.array( [-0.03,-0.03, -0.03,-1.0 ], dtype=np.float32 );
    //DEFAULT_dpos0        = np.array( [ 0.0 , 0.0 , -4.0 , 4.0 ], dtype=np.float32 );
    //DEFAULT_relax_params = np.array( [ 0.5 , 0.1 ,  0.02, 0.5 ], dtype=np.float32 );
    cl_program program_relax=0;

    int4   nDOFs    {0,0,0,0};
    int4   nPBC     {0,0,0,0};
    float4 md_params{0.05,0.9,100.0,0.0};    // (dt,cdamp,forceLimit)
    float4 dinv[3];    // grid step
    float4 tipRot[3]={{1.,0.,0.,0.},{0.,1.,0.,0.},{0.,0.,1.,0.1}};  // tip rotation
    //float4 tipRot[3]={{-1.,0.,0.,0.},{0.,-1.,0.,0.},{0.,0.,-1.,0.1}};  // tip rotation
    float4 stiffness { -0.03f,-0.03f,-0.03f,-1.00 };  // tip stiffness
    float4 dTip {0.0f,0.0f,-0.1f,0.0f};       // step of tip approch
    float4 dpos0{0.0f,0.0f,-4.0f, 4.0f};      // shift of initial positions
    //float4 dpos0{0.0f,0.0f, 0.0f,0.0f};      // shift of initial positions
    float4 relax_params{0.5f,0.1f,0.02f,0.5f};
    float4 surfFF;
    float4 tipQs {0.f,-0.05f,0.f,0.0f};
    //float4 tipQs {-10.f,+20.f,-10.f,0.0f};
    float4 tipQZs{0.1f, 0.0f,-0.1f,0.0f};
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

    void setGridShape( const Mat3d& dCell ){
        grid.n = {(int)Ns[0],(int)Ns[1],(int)Ns[2]};
        grid.dCell = dCell;
        grid.updateCell_2();
        v2f4(grid.dCell.a,dA);
        v2f4(grid.dCell.b,dB);
        v2f4(grid.dCell.c,dC);
        v2f4(grid.diCell.a,dinv[0]);
        v2f4(grid.diCell.b,dinv[1]);
        v2f4(grid.diCell.c,dinv[2]);
        //v2f4(grid.diCell.a*(1./grid.n.x),dinv[0]);
        //v2f4(grid.diCell.b*(1./grid.n.y),dinv[1]);
        //v2f4(grid.diCell.c*(1./grid.n.z),dinv[2]);
        //dinv[0].mul( 1/grid.n.x );
        //dinv[1].mul( 1/grid.n.y );
        //dinv[2].mul( 1/grid.n.z );
        //printf( "grid.diCell \n" );
        //println( grid.dCell.a );
        //println( grid.dCell.b );
        //println( grid.dCell.c );
        //printf( "grid.diCell \n" );
        //println( grid.diCell.a );
        //println( grid.diCell.b );
        //println( grid.diCell.c );
        //printf("DEBUG setGridShape dinvA(%g,%g,%g) dinvA(%g,%g,%g) dinvC(%g,%g,%g)\n",   dinv[0].x,dinv[0].y,dinv[0].z,   dinv[1].x,dinv[1].y,dinv[1].z,   dinv[2].x,dinv[2].y,dinv[2].z );
    }

    void makeStartPointGrid( Vec2i ns, Vec3d p0, Vec3d da, Vec3d db ){
        //printf( "makeStartPointGrid() nxy(%i,%i) p0(%g,%g,%g) da(%g,%g,%g) db(%g,%g,%g) \n", ns.x, ns.y,  p0.x,p0.y,p0.z,  da.x,da.y,da.z,  db.x,db.y,db.z );
        n_start_point = ns.x*ns.y;
        float4* ps = new float4[n_start_point];
        for(int ix=0; ix<ns.x; ix++){
            for(int iy=0; iy<ns.y; iy++){
                Vec3f p = (Vec3f)(p0 + (da*ix) + (db*iy));
                //ps[ ix+ ns.x*iy] = (float4){ p.x,p.y,p.z,0.0f };
                ps[iy+ns.y*ix] = (float4){ p.x,p.y,p.z,0.0f };
            }
        }
        //printf( "makeStartPointGrid() 1 \n" );
        if(ibuff_start_point<0){ ibuff_start_point=newBuffer( "StartPoints", n_start_point, sizeof(float4) ); }
        //printf( "makeStartPointGrid() 2 \n" );
        upload( ibuff_start_point, ps, n_start_point );
        //printf( "makeStartPointGrid() 2 \n" );
        delete [] ps;
    }

    void makeKrenels_PP( const char*  cl_src_dir ){
        printf( "makeKrenels_PP() \n" );
        char srcpath[1024];
        sprintf( srcpath, "%s/relax.cl", cl_src_dir );     
        buildProgram( srcpath, program_relax );
        newTask( "getFEinStrokes"         ,program_relax);
        newTask( "relaxStrokesTilted"     ,program_relax);
        newTask( "evalLJC_QZs"            ,program_relax);
        newTask( "evalLJC_QZs_toImg"      ,program_relax);
        newTask( "make_GridFF"            ,program_relax);
        newTask( "getNonBond"             ,program_relax);
        newTask( "getNonBondForce_GridFF" ,program_relax);
        newTask( "getMMFFsp3"             ,program_relax);
        newTask( "getMMFFf4"              ,program_relax);
        newTask( "cleanForceMMFFf4"       ,program_relax);
        newTask( "updateAtomsMMFFf4"      ,program_relax);
        newTask( "gatherForceAndMove"     ,program_relax);
        newTask( "updatePiPos0"           ,program_relax);
        newTask( "evalPiPi"               ,program_relax);
        //newTask( "write_toImg"     ,program_relax, 3,{0,0,0,0},{1,1,1,0} ); 
        printf( "... makeKrenels_PP() DONE \n" );
    }

    int initPP( const char*  cl_src_dir ){
        makeKrenels_PP( cl_src_dir );
        printf( "initPP() Ns(%li,%li,%li) \n", Ns[0], Ns[1], Ns[2] );
        itex_FF = newBufferImage3D( "FF", Ns[0], Ns[1], Ns[2], sizeof(float)*4, 0, CL_MEM_READ_WRITE, {CL_RGBA, CL_FLOAT} );
        return itex_FF;
    }

/*
    float4* debug_gen_FE( ){ 
        float4* data = new float4[Ntot];
        float d=0.02;
        printf( "Ns (%li,%li,%li) \n", Ns[0],Ns[1],Ns[2] );
        for(int ix=0; ix<(int)Ns[0]; ix++){
            //printf( "ix %i \n", ix );
            for(int iy=0; iy<(int)Ns[1]; iy++){
                for(int iz=0; iz<(int)Ns[2]; iz++){
                    //float fx=sin(ix*d);
                    //float fy=sin(iy*d);
                    //float fz=sin(iz*d); 
                    float fx=ix;
                    float fy=iy;
                    float fz=iz; 
                    //int i = iz+Ns[2]*(iy + Ns[1]*ix );
                    int i = ix+Ns[0]*(iy + Ns[1]*iz );
                    if(i>(int)Ntot){ printf("i(%i)>Ntot(%li) \n", i, Ntot ); exit(0); }
                    data[ i ] = (float4){fx,fy,fz, sqrt( fabs(fx*fy*fz)) };
                }
            }
        }
        return data;
    }
*/
    void getFEinStrokes( int ibuff_out, int nz_, Vec3d dTip_, int np=0, float4* points_=0 ){
        nz=nz_;
        dTip=(float4){(float)dTip_.x, (float)dTip_.y, (float)dTip_.z,0};
        OCLtask* task = tasks[ task_dict["getFEinStrokes"] ];
        task->global.x = n_start_point;
        //printf( "DEBUG roll_buf iKernell_roll %i ibuffA %i ibuffB %i \n", iKernell_roll, ibuffA, ibuffB );
        if( points_ != 0){
            upload( ibuff_start_point, points_, np);
        }
        useKernel( task->ikernel );
        int err=0;
        err |= useArgBuff( itex_FF      );      // 1
        err |= useArgBuff( ibuff_start_point ); // 2
        err |= useArgBuff( ibuff_out    );      // 3
        err |= _useArg( dinv[0] );              // 4
        err |= _useArg( dinv[1] );              // 5
        err |= _useArg( dinv[2] );              // 6
        err |= _useArg( dTip  );                // 7
        err |= _useArg( dpos0 );                // 8
        err |= useArg( nz );                // 9
        //printf("DEBUG getFEinStrokes 1 dinvA(%g,%g,%g) dinvA(%g,%g,%g) dinvC(%g,%g,%g)\n",   dinv[0].x,dinv[0].y,dinv[0].z,   dinv[1].x,dinv[1].y,dinv[1].z,   dinv[2].x,dinv[2].y,dinv[2].z );
        //printf("DEBUG getFEinStrokes 1 (%g,%g,%g)\n", dTip.x,dTip.y,dTip.z );
        OCL_checkError(err, "getFEinStrokes_1");
        err = task->enque();
        //err = task->enque( 3, *(size_t4*)&Ns, (size_t4){1,1,1,1} );
        OCL_checkError(err, "getFEinStrokes_2");  
    }

    void relaxStrokesTilted( int ibuff_out, int nz_, float dtip, int np=0, float4* points_=0 ){
        nz=nz_;
        tipRot[2].w = dtip;
        OCLtask* task = tasks[ task_dict["relaxStrokesTilted"] ];
        //printf( "relaxStrokesTilted() n_start_point %i \n", n_start_point );
        task->global.x = n_start_point;
        //printf( "DEBUG roll_buf iKernell_roll %i ibuffA %i ibuffB %i \n", iKernell_roll, ibuffA, ibuffB );
        if( points_ != 0){
            upload( ibuff_start_point, points_, np);
        }
        useKernel( task->ikernel );
        int err=0;
        err |= useArgBuff( itex_FF      );      // 1
        err |= useArgBuff( ibuff_start_point ); // 2
        err |= useArgBuff( ibuff_out    );      // 3
        err |= _useArg( dinv[0] );              // 4
        err |= _useArg( dinv[1] );              // 5
        err |= _useArg( dinv[2] );              // 6
        err |= _useArg( tipRot[0] );            // 7
        err |= _useArg( tipRot[1] );            // 8
        err |= _useArg( tipRot[2] );            // 9
        err |= _useArg( stiffness );            // 10
        err |= _useArg( dpos0 );                // 11
        err |= _useArg( relax_params );         // 12
        err |= _useArg( surfFF );               // 13
        err |= useArg( nz );                    // 14
        OCL_checkError(err, "relaxStrokesTilted_1");
        err = task->enque();
        //err = task->enque( 3, *(size_t4*)&Ns, (size_t4){1,1,1,1} );
        OCL_checkError(err, "relaxStrokesTilted_2");  
    }

    void evalLJC_QZs( int ibuff_out, int na=0, float4* atoms=0, float4* coefs=0 ){
        OCLtask* task = tasks[ task_dict["evalLJC_QZs"] ];
        task->global.x = Ntot;
        //task->global.y = Ns[1];
        //task->global.z = Ns[2];
        //printf( "DEBUG roll_buf iKernell_roll %i ibuffA %i ibuffB %i \n", iKernell_roll, ibuffA, ibuffB );
        //if(ibuff_atoms<0)initAtoms( na, na );
        if(atoms)upload( ibuff_atoms, atoms, na);
        if(coefs)upload( ibuff_coefs, coefs, na);
        useKernel( task->ikernel );
        printf("ibuff_out %i \n", ibuff_out);
        int4 ngrid{ (int)Ns[0],(int)Ns[1],(int)Ns[2],(int)Ns[3] };
        int err=0;
        err |= useArg    ( nAtoms        ); // 1
        err |= useArgBuff( ibuff_atoms   ); // 2
        err |= useArgBuff( ibuff_coefs   ); // 3
        err |= useArgBuff( ibuff_out     ); // 4
        //err |= useArg    ( (int)Ns     ); // 5
        err |= _useArg( ngrid  );         // 5        
        err |= _useArg( pos0 );           // 6
        err |= _useArg( dA );             // 7
        err |= _useArg( dB );             // 8
        err |= _useArg( dC );             // 9
        err |= _useArg( tipQs );             // 10
        err |= _useArg( tipQZs );            // 11
        OCL_checkError(err, "evalLJC_QZs_1");
        err = task->enque();
        OCL_checkError(err, "evalLJC_QZs_2");  
        /*
        __kernel void evalLJC_QZs(
            const int nAtoms,        // 1
            __global float4* atoms,  // 2
            __global float2*  cLJs,  // 3
            __global float4*    FE,  // 4
            int4 nGrid,              // 5
            float4 grid_p0,          // 6 
            float4 grid_dA,          // 7
            float4 grid_dB,          // 8
            float4 grid_dC,          // 9
            float4 Qs,               // 10
            float4 QZs               // 11
        */
    }

    void evalLJC_QZs_toImg( int na=0, float4* atoms=0, float4* coefs=0 ){        
        OCLtask* task = tasks[ task_dict["evalLJC_QZs_toImg"] ];
        task->global.x = Ntot;
        //task->global.y = Ns[1];
        //task->global.z = Ns[2];
        //printf( "DEBUG roll_buf iKernell_roll %i ibuffA %i ibuffB %i \n", iKernell_roll, ibuffA, ibuffB );
        if(ibuff_atoms<0)initAtoms( na, na );
        if(atoms)upload( ibuff_atoms, atoms, na);
        if(coefs)upload( ibuff_coefs, coefs, na);
        useKernel( task->ikernel );
        int4 ngrid{ (int)Ns[0],(int)Ns[1],(int)Ns[2],(int)Ns[3] };
        int err=0;
        err |= useArg    ( nAtoms        ); // 1
        err |= useArgBuff( ibuff_atoms   ); // 2
        err |= useArgBuff( ibuff_coefs   ); // 3
        err |= useArgBuff( itex_FF       ); // 4
        //err |= useArg    ( (int)Ns     ); // 5
        err |= _useArg( ngrid  );          // 5        
        err |= _useArg( pos0 );            // 6
        err |= _useArg( dA );              // 7
        err |= _useArg( dB );              // 8
        err |= _useArg( dC );              // 9
        err |= _useArg( tipQs );           // 10
        err |= _useArg( tipQZs );          // 11
        OCL_checkError(err, "evalLJC_QZs_toImg_1");
        err = task->enque_raw();
        OCL_checkError(err, "evalLJC_QZs_toImg_2");  
        /*
        __kernel void evalLJC_QZs(
            const int nAtoms,        // 1
            __global float4* atoms,  // 2
            __global float2*  cLJs,  // 3
            __write_only image3d_t  imgOut, // 4
            int4 nGrid,              // 5
            float4 grid_p0,          // 6 
            float4 grid_dA,          // 7
            float4 grid_dB,          // 8
            float4 grid_dC,          // 9
            float4 Qs,               // 10
            float4 QZs               // 11
        */
    }

    void makeGridFF( int na=0, float4* atoms=0, float4* coefs=0 ){
        if(itex_FE_Paul<=0) itex_FE_Paul = newBufferImage3D( "FEPaul", Ns[0], Ns[1], Ns[2], sizeof(float)*4, 0, CL_MEM_READ_WRITE, {CL_RGBA, CL_FLOAT} );
        if(itex_FE_Lond<=0) itex_FE_Lond = newBufferImage3D( "FFLond", Ns[0], Ns[1], Ns[2], sizeof(float)*4, 0, CL_MEM_READ_WRITE, {CL_RGBA, CL_FLOAT} );
        if(itex_FE_Coul<=0) itex_FE_Coul = newBufferImage3D( "FFCoul", Ns[0], Ns[1], Ns[2], sizeof(float)*4, 0, CL_MEM_READ_WRITE, {CL_RGBA, CL_FLOAT} );
        //OCLtask* task = tasks[ task_dict["make_GridFF"] ];
        OCLtask* task = getTask("make_GridFF");
        task->global.x = Ntot;
        //task->local.x  = 32;
        //task->roundSizes();
        //task->global.x = roundUpSize(Ntot,task->local.x);
        //task->global.y = Ns[1];
        //task->global.z = Ns[2];
        //printf( "DEBUG roll_buf iKernell_roll %i ibuffA %i ibuffB %i \n", iKernell_roll, ibuffA, ibuffB );
        if(ibuff_atoms<0)initAtoms( na, na );
        if(atoms)upload( ibuff_atoms, atoms, na );
        if(coefs)upload( ibuff_coefs, coefs, na );
        useKernel( task->ikernel );
        int4 ngrid{ (int)Ns[0],(int)Ns[1],(int)Ns[2],(int)Ns[3] };
        int err=0;
        err |= useArg    ( nAtoms       ); // 1
        err |= useArgBuff( ibuff_atoms  ); // 2
        err |= useArgBuff( ibuff_coefs  ); // 3
        err |= useArgBuff( itex_FE_Paul ); // 4
        err |= useArgBuff( itex_FE_Lond ); // 5
        err |= useArgBuff( itex_FE_Coul ); // 6
        err |= _useArg( ngrid  );          // 7        
        err |= _useArg( pos0 );            // 8
        err |= _useArg( dA );              // 9
        err |= _useArg( dB );              // 10
        err |= _useArg( dC );              // 11
        OCL_checkError(err, "makeGridFF.setup");
        err |= task->enque_raw(); OCL_checkError(err, "makeGridFF.enque");
        err |= finishRaw();       OCL_checkError(err, "makeGridFF.finish");
        /*
            const int nAtoms,                // 1
            __global float4*  atoms,         // 2
            __global float4*  REKQs,         // 3
            //__global float4*  FE_Paul,     // 4
            //__global float4*  FE_Lond,     // 5
            //__global float4*  FE_Coul,     // 6
            __write_only image3d_t  FE_Paul, // 4
            __write_only image3d_t  FE_Lond, // 5
            __write_only image3d_t  FE_Coul, // 6
            int4 nGrid,         // 7
            float4 grid_p0,     // 8
            float4 grid_dA,     // 9
            float4 grid_dB,     // 10
            float4 grid_dC      // 11
        */
    }

    int initAtomsForces( int nAtoms_, int npi, int nnode, bool bMMFFsp3=false, bool bMMFFf4=false ){
        nAtoms=nAtoms_;
        int nvecs=nAtoms+npi;
        int nneigh=nnode*4;
        if(bMMFFf4)nneigh*=2;
        ibuff_atoms     =newBuffer( "atoms",     nvecs,  sizeof(float4), 0, CL_MEM_READ_WRITE );
        ibuff_aforces   =newBuffer( "aforces",   nvecs,  sizeof(float4), 0, CL_MEM_READ_WRITE );
        ibuff_coefs     =newBuffer( "coefs",     nAtoms, sizeof(float4), 0, CL_MEM_READ_ONLY  );
        ibuff_neighs    =newBuffer( "neighs",    nAtoms, sizeof(int4  ), 0, CL_MEM_READ_ONLY  ); // need neihgs for all atoms because of Non-Bonded
        ibuff_neighCell =newBuffer( "neighCell" ,nAtoms, sizeof(int4  ), 0, CL_MEM_READ_ONLY  );
        if(bMMFFsp3 || bMMFFf4){
            ibuff_bkNeighs    = newBuffer( "bkNeighs",   nvecs,  sizeof(int4  ), 0, CL_MEM_READ_ONLY  );
            ibuff_avel        = newBuffer( "avel",       nvecs,  sizeof(float4), 0, CL_MEM_READ_WRITE );
            ibuff_neighForce  = newBuffer( "neighForce", nneigh, sizeof(float4), 0, CL_MEM_READ_WRITE );
            if(bMMFFsp3){
                printf( "initAtomsForces bMMFFsp3==true\n" );
                ibuff_bondLK      = newBuffer( "bondLK ",    nAtoms   , sizeof(float8), 0, CL_MEM_READ_ONLY  );
                ibuff_ang0K       = newBuffer( "ang0K",      nnode    , sizeof(float4), 0, CL_MEM_READ_ONLY  );
                ibuff_pi0s        = newBuffer( "pi0s",       npi      , sizeof(float4), 0, CL_MEM_READ_WRITE );
            }
            if(bMMFFf4){ // int ibuff_MMpars=-1, ibuff_BLs=-1,ibuff_BKs=-1,ibuff_Ksp=-1, ibuff_Kpp=-1;   // MMFFf4 params
                printf( "initAtomsForces bMMFFf4==true\n" );
                ibuff_MMpars = newBuffer( "MMpars", nnode, sizeof(int4),   0, CL_MEM_READ_ONLY );
                ibuff_BLs    = newBuffer( "BLs",    nnode, sizeof(float4), 0, CL_MEM_READ_ONLY  );
                ibuff_BKs    = newBuffer( "BKs",    nnode, sizeof(float4), 0, CL_MEM_READ_ONLY  );
                ibuff_Ksp    = newBuffer( "Ksp",    nnode, sizeof(float4), 0, CL_MEM_READ_ONLY );
                ibuff_Kpp    = newBuffer( "Kpp",    nnode, sizeof(float4), 0, CL_MEM_READ_ONLY );
            }
        }
        //printBuffers();
        //exit(0);
        return ibuff_atoms;
    }


    OCLtask* setup_getNonBond( int na, int nNode, Vec3i nPBC_, OCLtask* task=0){
        printf("!!!!! setup_getNonBond(na=%i,nnode=%i) \n", na, nNode);
        if(task==0) task = getTask("getNonBond");
        task->global.x = na;
        useKernel( task->ikernel );
        nDOFs.x=na; 
        nDOFs.y=nNode; 
        //nDOFs.x=bPBC; 
        nPBC.x = nPBC_.x;
        nPBC.y = nPBC_.y;
        nPBC.z = nPBC_.z;
        // ------- Maybe We do-not need to do this every frame ?
        int err=0;
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


    OCLtask* setup_getNonBondForce_GridFF( OCLtask* task=0, int na=-1 ){
        if(task==0) task = getTask("getNonBondForce_GridFF");
        if(na>=0  ) task->global.x = na;
        task->global.x = na;
        useKernel( task->ikernel );
        int err=0;
        err |= useArg    ( nAtoms       ); // 1
        err |= useArgBuff( ibuff_atoms  ); // 2
        err |= useArgBuff( ibuff_coefs  ); // 3
        err |= useArgBuff( ibuff_aforces); // 4
        err |= useArgBuff( ibuff_neighs ); // 5
        err |= useArgBuff( itex_FE_Paul ); // 6
        err |= useArgBuff( itex_FE_Lond ); // 7
        err |= useArgBuff( itex_FE_Coul ); // 8     
        err |= _useArg( pos0 );            // 9
        err |= _useArg( dinv[0] );              // 10
        err |= _useArg( dinv[1] );              // 11
        err |= _useArg( dinv[2] );              // 12
        OCL_checkError(err, "setup_getNonBondForce_GridFF");
        return task;
        /*
            const int nAtoms,               // 1
            __global float4*  atoms,        // 2
            __global float4*  REKQs,        // 3
            __global float4*  forces,       // 4
            __read_only image3d_t  FE_Paul, // 5
            __read_only image3d_t  FE_Lond, // 6
            __read_only image3d_t  FE_Coul, // 7
            float4 pos0,     // 8
            float4 dinvA,    // 9
            float4 dinvB,    // 10
            float4 dinvC     // 11
        */
    }

    void getNonBondForce_GridFF( int na=0, float4* atoms=0, float4* coefs=0, float4* aforces=0, int4* neighs=0 ){
        //printf("getNonBondForce_GridFF(na=%i) \n", na);
        if(ibuff_atoms<0)initAtoms( na, 1 );
        if(atoms  )upload( ibuff_atoms,   atoms,  na); // Note - these are other atoms than used for makeGridFF()
        if(coefs  )upload( ibuff_coefs,   coefs,  na);
        //if(coefs  )upload( ibuff_neighs,  neighs, na);
        //if(aforces)upload( ibuff_aforces, aforces, na);
        OCLtask* task = setup_getNonBondForce_GridFF( 0, na );
        int err=0;
        err = task->enque_raw();
        OCL_checkError(err, "getNonBondForce_GridFF_2");  
        if(aforces)err=download( ibuff_aforces, aforces, na);
        OCL_checkError(err, "getNonBondForce_GridFF_3");  
    }

    




    OCLtask* setup_getMMFFsp3( int na, int nNode, bool bPBC=false, OCLtask* task=0){
        //printf("setup_getMMFFsp3(na=%i,nnode=%i) \n", na, nNode);
        if(task==0) task = getTask("getMMFFsp3");
        task->global.x = na;
        //task->local .x = 1;
        //task->roundSizes();
        //if(na>=0  ) 
        useKernel( task->ikernel );
        nDOFs.x=na; 
        nDOFs.y=nNode; 
        nDOFs.w=bPBC; 
        int err=0;
        // ------- Maybe We do-not need to do this every frame ?
        //err |= useArg    ( nAtoms       );   // 1
        err |= _useArg   ( nDOFs );            // 1
        err |= useArgBuff( ibuff_atoms  );     // 2
        err |= useArgBuff( ibuff_coefs  );     // 3
        err |= useArgBuff( ibuff_aforces);     // 4
        err |= useArgBuff( ibuff_neighs );     // 5
        err |= useArgBuff( ibuff_bondLK );     // 6
        err |= useArgBuff( ibuff_ang0K  );     // 7
        err |= useArgBuff( ibuff_neighForce ); // 8
        err |= useArgBuff( itex_FE_Paul );     // 9
        err |= useArgBuff( itex_FE_Lond );     // 10
        err |= useArgBuff( itex_FE_Coul );     // 11    
        err |= _useArg( pos0 );                // 12
        err |= _useArg( dinv[0] );             // 13
        err |= _useArg( dinv[1] );             // 14
        err |= _useArg( dinv[2] );             // 15
        err |= _useArg( cl_lvec    );          // 16
        err |= _useArg( cl_invLvec );          // 16
        OCL_checkError(err, "setup_getMMFFsp3");
        return task;
        /*
        __kernel void getMMFFsp3(
            const int nAtoms,               // 1
            __global float4*  atoms,        // 2
            __global float4*  REQKs,        // 3
            __global float4*  forces,       // 4
            __global int4*    neighs,       // 5
            __global float8*  bondLK,       // 6 
            __global float2*  ang0K,        // 7
            __global float4*  neighForces,  // 8
            __read_only image3d_t  FE_Paul, // 9
            __read_only image3d_t  FE_Lond, // 10
            __read_only image3d_t  FE_Coul, // 11
            const float4 pos0,     // 12
            const float4 dinvA,    // 13
            const float4 dinvB,    // 14
            const float4 dinvC     // 15
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
        int err=0;
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
        int err=0;
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
        int err=0;
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

    OCLtask* setup_gatherForceAndMove( int n, int natom,  OCLtask* task=0 ){
        //printf("setup_gatherForceAndMove(na=%i) \n", n);
        //if(ibuff_atoms<0)initAtoms( na, 1 );
        //if(atoms  )upload( ibuff_atoms,   atoms,  na); // Note - these are other atoms than used for makeGridFF()
        //if(coefs  )upload( ibuff_coefs,   coefs,  na);
        //if(coefs  )upload( ibuff_neighs,  neighs, na);
        //if(aforces)upload( ibuff_aforces, aforces, na);
        if(task==0) task = getTask("gatherForceAndMove");
        task->global.x = n;
        //task->local .x = 1;
        //task->roundSizes();
        //if(n >=0  ) 
        nDOFs.x=n; 
        nDOFs.y=natom; 
        useKernel( task->ikernel );
        int err=0;
        err |= _useArg( md_params );           // 1
        err |= _useArg( nDOFs     );           // 2
        err |= useArgBuff( ibuff_atoms  );     // 3
        err |= useArgBuff( ibuff_avel   );     // 4
        err |= useArgBuff( ibuff_aforces);     // 5
        err |= useArgBuff( ibuff_neighForce ); // 6
        err |= useArgBuff( ibuff_bkNeighs   );   // 7
        err |= useArgBuff( ibuff_pi0s       );   // 8
        OCL_checkError(err, "setup_gatherForceAndMove");
        //err = task->enque_raw();
        //OCL_checkError(err, "gatherForceAndMove");  
        //if(aforces)err=download( ibuff_aforces, aforces, na);
        //OCL_checkError(err, "gatherForceAndMove");  
        return task;
        /*
            __kernel void gatherForceAndMove(
                const float4      params,       // 1
                const int         nAtoms,       // 2
                __global float4*  apos,         // 3
                __global float4*  avel,         // 4
                __global float4*  aforce,       // 5
                __global float4*  neighForces,  // 6
                __global int4*    bkNeighs      // 7
            ){
        */
    }


    OCLtask* setup_updatePiPos0( int natom, int npi, OCLtask* task=0 ){
        if(task==0) task = getTask("updatePiPos0");
        task->global.x = npi;
        //task->local .x = 1;
        //task->roundSizes();
        //if(n >=0  ) 
        nDOFs.x=natom; 
        nDOFs.y=npi; 
        //printf("setup_updatePiPos0 nDOFs(natom=%i,npi=%i) \n", nDOFs.x,nDOFs.y );
        useKernel( task->ikernel );
        int err=0;
        err |= _useArg( nDOFs     );            // 1
        err |= useArgBuff( ibuff_atoms   );     // 2
        err |= useArgBuff( ibuff_pi0s    );     // 3
        err |= useArgBuff( ibuff_neighs  );     // 4
        OCL_checkError(err, "setup_updatePiPos0");
        return task;
        /*
            __kernel void updatePiPos0(
                const int4        n,           // 1
                __global float4*  apos,        // 2
                __global float4*  pi0s,        // 3
                __global int4*    pi_neighs    // 4
            ){
        */
    }


    OCLtask* setup_evalPiPi( int natom, int npi, OCLtask* task=0 ){
        if(task==0) task = getTask("evalPiPi");
        task->global.x = npi;
        //task->local .x = 1;
        //task->roundSizes();
        //if(n >=0  ) 
        nDOFs.x=natom; 
        nDOFs.y=npi; 
        //printf("setup_updatePiPos0 nDOFs(natom=%i,npi=%i) \n", nDOFs.x,nDOFs.y );
        useKernel( task->ikernel );
        int err=0;
        err |= _useArg( nDOFs    );             // 1
        err |= useArgBuff( ibuff_atoms   );     // 2
        err |= useArgBuff( ibuff_aforces );     // 3
        err |= useArgBuff( ibuff_neighs  );     // 4
        OCL_checkError(err, "setup_evalPiPi");
        return task;
        /*
            __kernel void evalPiPi(
                const int4        n,        // 1
                __global float4*  apos,     // 2
                __global float4*  aforce,   // 2
                __global int4*    neighs    // 4
            ){
        */
    }


};

#endif
