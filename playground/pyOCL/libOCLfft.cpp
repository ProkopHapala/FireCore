

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
    size_t Ns[3]; // = {N0, N1, N2};
    size_t Ntot;
    //size_t buffer_size;

    //static cl_mem data_cl;
    //static float *data;

    //int N=0;
    //int size = N * N;
    //static cl_mem d_a, d_b, d_c;

    int iKernell_mull=-1;
    int iKernell_project=-1;
    OCLtask cltask_mul;
    OCLtask cltask_project;

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
        cltask_mul.setup( this, iKernell_mull, 1, Ntot*2, 16 );
        cltask_mul.args = { 
            INTarg (cltask_mul.global[0]),
            BUFFarg(ibuffA),
            BUFFarg(ibuffB),
            BUFFarg(ibuff_result),
        };
    }

    void initTask_project( int ibuffAtoms, int ibuffCoefs, int ibuff_result ){
        pos0=(float4){0.0f,0.0f,0.0f,0.0f};
        dA  =(float4){0.1f,0.0f,0.0f,0.0f};
        dB  =(float4){0.0f,0.1f,0.0f,0.0f};
        dC  =(float4){0.0f,0.0f,0.1f,0.0f};
        cltask_project.setup( this, iKernell_project, 1, Ntot*2, 16 );
        cltask_project.args = { 
            INTarg (nAtoms),
            BUFFarg(ibuffAtoms),
            BUFFarg(ibuffCoefs),
            BUFFarg(ibuff_result),
            REFarg(pos0),
            REFarg(dA),
            REFarg(dB),
            REFarg(dC)
        };
    }

    void projectAtoms( float4* atoms, float4* coefs, int ibuff_result ){
        upload(ibuff_atoms,atoms);
        upload(ibuff_coefs,coefs);
        initTask_project( ibuff_atoms, ibuff_coefs, ibuff_result );
        cltask_project.enque( );
        //upload(ibuff_coefs,coefs);
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

    int initAtoms( int nAtoms_ ){
        oclfft.nAtoms=nAtoms_;
        int ibuff0=oclfft.newBuffer( "atoms", nAtoms_, sizeof(float4), 0, CL_MEM_READ_ONLY );
                   oclfft.newBuffer( "coefs", nAtoms_, sizeof(float4), 0, CL_MEM_READ_ONLY );
        return ibuff0;
    };

    void runfft( int ibuff, bool fwd              ){ oclfft.runFFT( ibuff,fwd,0);     };
    //void runfft( int ibuff, bool fwd, float* data ){ oclfft.runFFT( ibuff,fwd, data); };
    void convolve( int ibuffA, int ibuffB, int ibuff_result ){
        oclfft.convolution( ibuffA, ibuffB, ibuff_result  );
        //oclfft.mul_buffs( ibuffA, ibuffB, ibuff_result );
    }
    void projectAtoms( float4* atoms, float4* coefs, int ibuff_result ){ oclfft.projectAtoms( atoms, coefs, ibuff_result ); }
    void cleanup(){ oclfft.cleanup(); }

};