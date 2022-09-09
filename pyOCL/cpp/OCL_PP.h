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
    float4 tipRot[3];  // tip rotation
    float4 stiffness { -0.03f,-0.03f,-0.03f,-1.00 };  // tip stiffness
    float4 dTip {0.0f,0.0f,-0.1f,0.0f};       // step of tip approch
    float4 dpos0{0.0f,0.0f,-4.0f,0.0f};      // shift of initial positions
    //float4 dpos0{0.0f,0.0f, 0.0f,0.0f};      // shift of initial positions
    float4 relax_params{0.5f,0.1f,0.02f,0.5f};
    float4 surfFF;
    //QZs             = [ 0.1,  0, -0.1, 0 ],
    //Qs              = [ -10, 20,  -10, 0 ],
    float4 tipQs {10.f,20.f,-10.f,0.0f};
    float4 tipQZs{0.1f,0.0f, 0.1f,0.0f};
    int nz;

    int n_start_point = 0;
    int ibuff_start_point=-1;
    int itex_FF=-1;

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
        for(int iy=0; iy<ns.y; iy++){
            for(int ix=0; ix<ns.x; ix++){
                Vec3f p = (Vec3f)(p0 + (da*ix) + (db*iy));
                ps[ ix+ ns.x*iy] = (float4){ p.x,p.y,p.z,0.0f };
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
        char srcpath[1024];
        sprintf( srcpath, "%s/relax.cl", cl_src_dir );       
        printf( "DEBUG makeKrenels_PP() %s \n", srcpath ); 
        buildProgram( srcpath, program_relax );
        newTask( "getFEinStrokes"    ,1,{0,0,0,0},{1,0,0,0},program_relax);
        newTask( "relaxStrokesTilted",1,{0,0,0,0},{1,0,0,0},program_relax);
        //tasks[ newTask( "relaxStrokesTilted",1,{0,0,0,0},{1,0,0,0},program_relax) ]->args={}; 
    }

    int initPP( const char*  cl_src_dir ){
        makeKrenels_PP( cl_src_dir );
        itex_FF = newBufferImage3D( "FF", Ns[0], Ns[1], Ns[1], sizeof(float)*4, 0, CL_MEM_READ_ONLY, {CL_RGBA, CL_FLOAT} );
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

    void relaxStrokesTilted( int ibuff_out, int np=0, float4* points_=0 ){
        OCLtask* task = tasks[ task_dict["relaxStrokesTilted"] ];
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
        OCL_checkError(err, "getFEinStrokes");
        err = task->enque();
        //err = task->enque( 3, *(size_t4*)&Ns, (size_t4){1,1,1,1} );
        OCL_checkError(err, "getFEinStrokes");  
    }

    void evalLJC_QZs( int ibuff_out, int na=0, float4* atoms=0, float4* coefs=0 ){
        OCLtask* task = tasks[ task_dict["relaxStrokesTilted"] ];
        task->global.x = n_start_point;
        //printf( "DEBUG roll_buf iKernell_roll %i ibuffA %i ibuffB %i \n", iKernell_roll, ibuffA, ibuffB );
        if(atoms)upload( ibuff_atoms, atoms, na);
        if(coefs)upload( ibuff_coefs, coefs, na);
        useKernel( task->ikernel );
        err |= useArg    ( nAtoms        );  // 1
        err |= useArgBuff( ibuff_atoms   );  // 2
        err |= useArgBuff( ibuff_coefs   );  // 3
        err |= useArgBuff( ibuff_out     );  // 4
        err |= useArg    ( (int)Ntot     );  // 5
        err |= _useArg( pos0 );           // 6
        err |= _useArg( dA );             // 7
        err |= _useArg( dB );             // 8
        err |= _useArg( dC );             // 9
        err |= _useArg( tipQs );             // 10
        err |= _useArg( tipQZs );            // 11
        OCL_checkError(err, "getFEinStrokes");
        err = task->enque();
        OCL_checkError(err, "getFEinStrokes");  
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

};

#endif