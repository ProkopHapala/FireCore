

// No need to explicitely include the OpenCL headers 
//#include <clFFT.h>
#include "libOCLfft.h"

#include  "OCL.h"

OCLsystem ocl;

typedef struct{
    cl_uint  dim;
    size_t global[3];
    size_t local[3];
    //cl_int blocksize;
} KernelDims;

class OCLfft : public OCLsystem { public:

    clfftPlanHandle planHandle;
    clfftDim fft_dim = CLFFT_3D;

    //const size_t N0 = 4, N1 = 4, N2 = 4;
    int ndim=0;
    size_t Ns[4]; // = {N0, N1, N2};
    size_t Ntot;
    int4   Nvec;
    //size_t buffer_size;

    //static cl_mem data_cl;
    //static float *data;

    //int N=0;
    //int size = N * N;
    //static cl_mem d_a, d_b, d_c;

    int iKernell_mull=-1;
    int iKernell_project=-1;
    int iKernell_project_tex=-1;
    OCLtask cltask_mul;
    OCLtask cltask_project;
    OCLtask cltask_project_tex;
    int itex_basis=-1;

    int    nAtoms;
    float4 pos0, dA, dB, dC;
    int ibuff_atoms,ibuff_coefs;

    void updateNtot(){
        Ntot=1; for(int i=0; i<ndim; i++){ Ntot=Ns[i]; };
    }

    void planFFT(){
        // Create a default plan for a complex FFT. 
        err = clfftCreateDefaultPlan(&planHandle, context, fft_dim, Ns );
        err = clfftSetPlanPrecision (planHandle, CLFFT_SINGLE);
        err = clfftSetLayout        (planHandle, CLFFT_COMPLEX_INTERLEAVED, CLFFT_COMPLEX_INTERLEAVED);
        err = clfftSetResultLocation(planHandle, CLFFT_INPLACE);
        err = clfftBakePlan         (planHandle, 1, &commands, NULL, NULL);
    }

    void initFFT( int ndim, size_t* Ns_ ){
        //printf("DEBUG initFFT() ndim %i Ns[%li,%li,%li]\n", ndim, Ns[0], Ns[1], Ns[2] );
        if     (ndim==1){ fft_dim=CLFFT_1D; }
        else if(ndim==2){ fft_dim=CLFFT_2D; }
        else if(ndim==3){ fft_dim=CLFFT_3D; };
        //buffer_size  = sizeof(float2);
        Ntot=1; for(int i=0; i<ndim;i++){  Ns[i]=Ns_[i];  Ntot*= Ns[i]; }
        printf( "initFFT ndim %i Ntot %li [%li,%li,%li]\n", ndim, Ntot, Ns[0],Ns[1],Ns[2] );
        clfftSetupData fftSetup;
        err  = clfftInitSetupData(&fftSetup);
        err  = clfftSetup        (&fftSetup);
        //data_cl = clCreateBuffer( context, CL_MEM_READ_WRITE, buffer_size, NULL, &err );
        planFFT(  );
    }

    int newFFTbuffer( char* name ){
        return newBuffer( name, Ntot, sizeof(float2), 0, CL_MEM_READ_WRITE );
    }

    int initAtoms( int nAtoms_ ){
        nAtoms=nAtoms_;
        ibuff_atoms=newBuffer( "atoms", nAtoms, sizeof(float4), 0, CL_MEM_READ_ONLY );
        ibuff_coefs=newBuffer( "coefs", nAtoms, sizeof(float4), 0, CL_MEM_READ_ONLY );
        return ibuff_atoms;
    };

    int initBasisTable( int nx, int ny, float* data ){
        itex_basis = newBufferImage2D( "BasisTable", nx, ny,   sizeof(float),  data, CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR , {CL_RGBA, CL_FLOAT} );
        return itex_basis;
    }

    //void loadData( float* data_ ){
        //printf("DEBUG loadData() \n");
        //data = data_;
        //makeData();
        //printData( data ); 
        //err = clEnqueueWriteBuffer  ( commands, data_cl, CL_TRUE, 0, buffer_size, data, 0, NULL, NULL );
    //}

    void runFFT( int ibuff, bool fwd, float* data=0 ){
        //err = clEnqueueWriteBuffer ( queue, data_cl, CL_TRUE, 0, buffer_size, data, 0, NULL, NULL );
        cl_mem data_cl =  buffers[ibuff].p_gpu;
        if(fwd){
            err = clfftEnqueueTransform( planHandle, CLFFT_FORWARD, 1, &commands, 0, NULL, NULL, &data_cl, NULL, NULL);  // Execute the plan. -> Forward Transform
        }else{
            err = clfftEnqueueTransform( planHandle, CLFFT_BACKWARD, 1, &commands, 0, NULL, NULL, &data_cl, NULL, NULL);  // Execute the plan.  -> Backward Transform      
        }                                                  
        if(data){
            err = clEnqueueReadBuffer  ( commands, data_cl, CL_TRUE, 0, buffers[ibuff].byteSize(), data, 0, NULL, NULL ); // Fetch results of calculations. 
            err = clFinish(commands);    // Wait for calculations to be finished. 
        }
        //printData( data ); 
    }

    void mul_buffs( int ibuffA, int ibuffB, int ibuff_result ){
        KernelDims kdim;
        kdim.dim        = 1;
        kdim.global[0]  = Ntot*2;
        kdim.local [0]  = 16;
        cl_kernel kernel = kernels[iKernell_mull]; 
        //int N2 = Ntot*2;
        err =  clSetKernelArg(kernel, 0, sizeof(int),    &kdim.global        );
        err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &(buffers[0].p_gpu) );
        err |= clSetKernelArg(kernel, 2, sizeof(cl_mem), &(buffers[1].p_gpu) );
        err |= clSetKernelArg(kernel, 3, sizeof(cl_mem), &(buffers[2].p_gpu) );
        //checkError(err, "Setting kernel args");
        //double start_time = wtime();
        printf( "mul_buffs kdim: dim %i global %li local %li \n", kdim.dim, kdim.global[0], kdim.local[0] ); 
        err = clEnqueueNDRangeKernel(  commands, kernel,   kdim.dim, NULL, kdim.global, kdim.local, 0, NULL, NULL);    
        OCL_checkError(err, "Enqueueing kernel");
        //err = clFinish(commands);
        //OCL_checkError(err, "Waiting for kernel to finish");
        //double run_time = wtime() - start_time;
    }

    void initTask_mul( int ibuffA, int ibuffB, int ibuff_result ){
        //cltask_mul.setup( this, iKernell_mull, 1, Ntot*2, 16 );
        cltask_mul.setup( this, iKernell_mull, 1, 8, 1 ); printf( "WARRNING : initTask_mul() IS WRONG !!!! %i %i \n", ibuffA, ibuffB );
        cltask_mul.args = { 
            INTarg (cltask_mul.global[0]),
            BUFFarg(ibuffA),
            BUFFarg(ibuffB),
            BUFFarg(ibuff_result),
        };
    }

    void initTask_project( int ibuffAtoms, int ibuffCoefs, int ibuff_result ){
        Nvec  =(int4){(int)Ns[0],(int)Ns[1],(int)Ns[2],(int)Ns[3]};
        cltask_project.setup( this, iKernell_project, 1, Ntot*2, 16 );
        cltask_project.args = { 
            INTarg (nAtoms),        //1
            BUFFarg(ibuffAtoms),    //2
            BUFFarg(ibuffCoefs),    //3
            BUFFarg(ibuff_result),  //4
            REFarg(Nvec),           //5
            REFarg(pos0),           //6
            REFarg(dA),           //7
            REFarg(dB),           //8
            REFarg(dC)            //9
        };
        //cltask_project.print_arg_list();
    }

    void initTask_project_tex( int ibuffAtoms, int ibuffCoefs, int ibuff_result ){
        Nvec  =(int4){(int)Ns[0],(int)Ns[1],(int)Ns[2],(int)Ns[3]};
        cltask_project_tex.setup( this, iKernell_project_tex, 1, Ntot*2, 16 );
        cltask_project_tex.args = { 
            INTarg (nAtoms),        //1
            BUFFarg(ibuffAtoms),    //2
            BUFFarg(ibuffCoefs),    //3
            BUFFarg(ibuff_result),  //4
            BUFFarg(itex_basis),    //5
            REFarg(Nvec),           //6
            REFarg(pos0),           //7
            REFarg(dA),           //8
            REFarg(dB),           //9
            REFarg(dC)            //10
        };
        cltask_project_tex.print_arg_list();
    }

    void projectAtoms( float4* atoms, float4* coefs, int ibuff_result ){
        //for(int i=0; i<nAtoms;i++){printf( "atom[%i] xyz|e(%g,%g,%g|%g) coefs(%g,%g,%g|%g)\n", i, atoms[i].x,atoms[i].y,atoms[i].z,atoms[i].w,  coefs[i].x,coefs[i].y,coefs[i].z,coefs[i].w  );}
        upload(ibuff_atoms,atoms);
        upload(ibuff_coefs,coefs);
        //ibuff_atoms=0;
        //ibuff_coefs=1;
        //printf("ibuff_atoms %i ibuff_coefs %i \n", ibuff_atoms, ibuff_coefs);
        //err = clEnqueueWriteBuffer  ( queue, data_cl,                   CL_TRUE, 0, buffer_size,           data,  0, NULL, NULL );
        //err = clEnqueueWriteBuffer ( commands, buffers[ibuff_atoms].p_gpu, CL_TRUE, 0, sizeof(float4)*nAtoms, atoms, 0, NULL, NULL ); OCL_checkError(err, "Creating ibuff_atoms");
        //err = clEnqueueWriteBuffer ( commands, buffers[ibuff_coefs].p_gpu, CL_TRUE, 0, sizeof(float4)*nAtoms, coefs, 0, NULL, NULL ); OCL_checkError(err, "Creating ibuff_coefs");
        clFinish(commands); 
        //initTask_project( ibuff_atoms, ibuff_coefs, ibuff_result );
        //cltask_project.enque( );
        initTask_project_tex( ibuff_atoms, ibuff_coefs, ibuff_result );
        cltask_project_tex.enque( );
        //initTask_mul( ibuff_atoms, ibuff_coefs,   ibuff_result   );
        //cltask_mul.enque( );
        clFinish(commands); 
    }
    
    void convolution( int ibuffA, int ibuffB, int ibuff_result ){
        err = clfftEnqueueTransform( planHandle, CLFFT_FORWARD, 1, &commands, 0, NULL, NULL, &buffers[ibuffA].p_gpu, NULL, NULL);
        err = clfftEnqueueTransform( planHandle, CLFFT_FORWARD, 1, &commands, 0, NULL, NULL, &buffers[ibuffB].p_gpu, NULL, NULL);  
        //mul_buffs( ibuffA, ibuffB, ibuff_result );
        initTask_mul( ibuffA, ibuffB, ibuff_result );  //cltask_mul.print_arg_list();
        cltask_mul.enque( );
        err = clfftEnqueueTransform( planHandle, CLFFT_BACKWARD, 1, &commands, 0, NULL, NULL, &buffers[ibuff_result].p_gpu, NULL, NULL);  
    }

    void cleanup(){
        //clReleaseMemObject( data_cl );
        //free( data );
        err = clfftDestroyPlan( &planHandle );
        clfftTeardown( );
        //clReleaseCommandQueue( commands );
        //clReleaseContext( context );
        release_OCL();
    }


    void makeMyKernels(){
        buildProgram( "myprog.cl" );
        iKernell_mull    = newKernel( "mul" );
        iKernell_project = newKernel( "projectAtomsToGrid" );
        iKernell_project_tex = newKernel( "projectAtomsToGrid_texture" );
    };


/*
    void runAll( ){
        //printf("DEBUG Ns[%li,%li,%li]\n", clLengths[0], clLengths[1], clLengths[2] );
        makeData();  printData( data );    printf( "DEBUG 1 \n" );
        initOCL();                         printf( "DEBUG 2 \n" );
        //printf("DEBUG Ns[%li,%li,%li]\n", clLengths[0], clLengths[1], clLengths[2] );
        initFFT( 3, clLengths );           printf( "DEBUG 3 \n" );
        loadData( data );
        runFFT();    printData(data );     printf( "DEBUG 5 \n" );
        cleanup();                         printf( "DEBUG 6 \n" );
    }
*/

};

OCLfft oclfft;

extern "C" {

    void  init(){
        oclfft.init();
        oclfft.makeMyKernels();
    }

    int   upload(int i, const float* cpu_data ){ return oclfft  .upload(i,cpu_data);                        };
    int download(int i,       float* cpu_data ){ return oclfft.download(i,cpu_data);  oclfft.finishRaw();   };

    void initFFT( int ndim, size_t* Ns_ ){
        oclfft.initFFT( ndim, Ns_ );
        oclfft.newFFTbuffer( "inputA" );
        oclfft.newFFTbuffer( "inputB" );
        oclfft.newFFTbuffer( "outputC" );
        //oclfft.initTask_mul( 0, 1, 2 );    // If we know arguments in front, we may define it right now
    }

    int initAtoms( int nAtoms          ){  return oclfft.initAtoms( nAtoms ); };
    void runfft( int ibuff, bool fwd   ){ oclfft.runFFT( ibuff,fwd,0);     };
    //void runfft( int ibuff, bool fwd, float* data ){ oclfft.runFFT( ibuff,fwd, data); };
    void convolve( int ibuffA, int ibuffB, int ibuff_result ){
        oclfft.convolution( ibuffA, ibuffB, ibuff_result  );
        //oclfft.mul_buffs( ibuffA, ibuffB, ibuff_result );
    }
    void projectAtoms( float* atoms, float* coefs, int ibuff_result ){ oclfft.projectAtoms( (float4*)atoms, (float4*)coefs, ibuff_result ); }
    void cleanup(){ oclfft.cleanup(); }

    void setGridShape( float* pos0, float* dA, float* dB, float* dC ){
        oclfft.pos0=*(float4*)pos0;
        oclfft.dA  =*(float4*)dA;
        oclfft.dB  =*(float4*)dB;
        oclfft.dC  =*(float4*)dC;
    }

    int initBasisTable( int nx, int ny, float* data ){  return oclfft.initBasisTable(nx,ny,data ); };

};