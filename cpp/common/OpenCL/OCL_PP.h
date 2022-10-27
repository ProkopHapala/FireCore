#ifndef  OCL_PP_h
#define  OCL_PP_h

#include "OCL_DFT.h"

//=======================================================================
//=======================================================================
class OCL_PP: public OCL_DFT { public:
    //DEFAULT_dTip         = np.array( [ 0.0 , 0.0 , -0.1 , 0.0 ], dtype=np.float32 );
    //DEFAULT_stiffness    = np.array( [-0.03,-0.03, -0.03,-1.0 ], dtype=np.float32 );
    //DEFAULT_dpos0        = np.array( [ 0.0 , 0.0 , -4.0 , 4.0 ], dtype=np.float32 );
    //DEFAULT_relax_params = np.array( [ 0.5 , 0.1 ,  0.02, 0.5 ], dtype=np.float32 );
    cl_program program_relax=0;

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
    int nz;

    int n_start_point = 0;
    int ibuff_start_point=-1;
    int itex_FF=-1;

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
        //v2f4(grid.diCell.a,dinv[0]);
        //v2f4(grid.diCell.b,dinv[1]);
        //v2f4(grid.diCell.c,dinv[2]);
        v2f4(grid.diCell.a*(1./grid.n.x),dinv[0]);
        v2f4(grid.diCell.b*(1./grid.n.y),dinv[1]);
        v2f4(grid.diCell.c*(1./grid.n.z),dinv[2]);
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
        newTask( "getFEinStrokes"    ,1,{0,0,0,0},{1,0,0,0},program_relax);
        newTask( "relaxStrokesTilted",1,{0,0,0,0},{1,0,0,0},program_relax);
        newTask( "evalLJC_QZs"       ,1,{0,0,0,0},{1,0,0,0},program_relax);
        newTask( "evalLJC_QZs_toImg" ,1,{0,0,0,0},{1,0,0,0},program_relax);
        newTask( "make_GridFF"        ,1,{0,0,0,0},{1,0,0,0},program_relax);
        newTask( "getNonBondForce_GridFF",1,{0,0,0,0},{1,0,0,0},program_relax);
        //newTask( "write_toImg"       ,3,{0,0,0,0},{1,1,1,0},program_relax);
        //tasks[ newTask( "relaxStrokesTilted",1,{0,0,0,0},{1,0,0,0},program_relax) ]->args={}; 
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
        //task->global.y = Ns[1];
        //task->global.z = Ns[2];
        //printf( "DEBUG roll_buf iKernell_roll %i ibuffA %i ibuffB %i \n", iKernell_roll, ibuffA, ibuffB );
        if(ibuff_atoms<0)initAtoms( na, na );
        if(atoms)upload( ibuff_atoms, atoms, na );
        if(coefs)upload( ibuff_coefs, coefs, na );
        useKernel( task->ikernel );
        int4 ngrid{ (int)Ns[0],(int)Ns[1],(int)Ns[2],(int)Ns[3] };
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
        OCL_checkError(err, "makeGridFF_1");
        err = task->enque_raw();
        OCL_checkError(err, "makeGridFF_2");  
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

    int initAtomsForces( int nAtoms_ ){
        nAtoms=nAtoms_;
        ibuff_atoms   =newBuffer( "atoms",    nAtoms, sizeof(float4), 0, CL_MEM_READ_ONLY );
        ibuff_coefs   =newBuffer( "coefs",    nAtoms, sizeof(float4), 0, CL_MEM_READ_ONLY );
        ibuff_aforces =newBuffer( "aforces",  nAtoms, sizeof(float4), 0, CL_MEM_READ_ONLY );
        return ibuff_atoms;
    }

    void getNonBondForce_GridFF( int na=0, float4* atoms=0, float4* coefs=0, float4* aforces=0 ){
        printf("getNonBondForce_GridFF(na=%i) \n", na);
        if(ibuff_atoms<0)initAtoms( na, 1 );
        if(atoms  )upload( ibuff_atoms,   atoms, na); // Note - these are other atoms than used for makeGridFF()
        if(coefs  )upload( ibuff_coefs,   coefs, na);
        DEBUG
        //if(aforces)upload( ibuff_aforces, aforces, na);
        OCLtask* task = getTask("getNonBondForce_GridFF");
        task->global.x = na;
        useKernel( task->ikernel );
        err |= useArg    ( nAtoms       ); // 1
        err |= useArgBuff( ibuff_atoms  ); // 2
        err |= useArgBuff( ibuff_coefs  ); // 3
        err |= useArgBuff( ibuff_aforces); // 4
        err |= useArgBuff( itex_FE_Paul ); // 5
        err |= useArgBuff( itex_FE_Lond ); // 6
        err |= useArgBuff( itex_FE_Coul ); // 7     
        err |= _useArg( pos0 );            // 8
        err |= _useArg( dA );              // 9
        err |= _useArg( dB );              // 10
        err |= _useArg( dC );              // 11
        OCL_checkError(err, "getNonBondForce_GridFF_1");
        err = task->enque_raw();
        OCL_checkError(err, "getNonBondForce_GridFF_2");  
        if(aforces)err=download( ibuff_aforces, aforces, na);
        OCL_checkError(err, "getNonBondForce_GridFF_3");  
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

};

#endif
