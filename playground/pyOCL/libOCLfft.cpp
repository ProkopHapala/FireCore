

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
        kdim.dim        = 2;
        //kdim.blocksize  = 16;
        kdim.global[0]  = Ntot*2;
        kdim.local [0]  = 16;
        //cl_uint  dim       = 1;
        //cl_int   blocksize = 16;
        //size_t global[]    = {Ntot*2};
        //size_t local []    = {blocksize};

        cl_kernel kernel = kernels[iKernell_mull]; 
        //int N2 = Ntot*2;
        err =  clSetKernelArg(kernel, 0, sizeof(int),    &kdim.global        );
        err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &(buffers[0].p_gpu) );
        err |= clSetKernelArg(kernel, 2, sizeof(cl_mem), &(buffers[1].p_gpu) );
        err |= clSetKernelArg(kernel, 3, sizeof(cl_mem), &(buffers[2].p_gpu) );
        //checkError(err, "Setting kernel args");
        //double start_time = wtime();
        err = clEnqueueNDRangeKernel(  commands, kernel,   kdim.dim, NULL, kdim.global, kdim.local, 0, NULL, NULL);    
        //checkError(err, "Enqueueing kernel");
        //err = clFinish(queue);
        //checkError(err, "Waiting for kernel to finish");
        //double run_time = wtime() - start_time;
    }

    void convolution( int ibuffA, int ibuffB, int ibuff_result ){
        err = clfftEnqueueTransform( planHandle, CLFFT_FORWARD, 1, &commands, 0, NULL, NULL, &buffers[ibuffA].p_gpu, NULL, NULL);
        err = clfftEnqueueTransform( planHandle, CLFFT_FORWARD, 1, &commands, 0, NULL, NULL, &buffers[ibuffB].p_gpu, NULL, NULL);  
        mul_buffs( ibuffA, ibuffB, ibuff_result );
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
        iKernell_mull = newKernel( "mul" );
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
    }

    void runfft( int ibuff, bool fwd, float* data ){ oclfft.runFFT( ibuff,fwd, data); };
    void convolve( int ibuffA, int ibuffB, int ibuff_result ){
        //oclfft.convolution( ibuffA, ibuffB, ibuff_result  );
        oclfft.mul_buffs( ibuffA, ibuffB, ibuff_result );
    }

    void cleanup(){ oclfft.cleanup(); }

};