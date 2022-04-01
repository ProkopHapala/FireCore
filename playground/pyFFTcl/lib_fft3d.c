
#include <stdio.h>
#include <stdlib.h>

// No need to explicitely include the OpenCL headers 
#include <clFFT.h>
#include "lib_fft3.h"

static    cl_int err;
static    cl_platform_id platform = 0;
static    cl_device_id device = 0;
static    cl_context_properties props[3] = { CL_CONTEXT_PLATFORM, 0, 0 };
static    cl_context ctx = 0;
static    cl_command_queue queue = 0;
static    cl_mem bufX;
static    float *X;
static    cl_event event = NULL;
static    int ret = 0;
static    const size_t N0 = 4, N1 = 4, N2 = 4;
static    char platform_name[128];
static    char device_name[128];

// FFT library realted declarations
static    clfftPlanHandle planHandle;
static    clfftDim dim = CLFFT_3D;
static    size_t clLengths[3] = {N0, N1, N2};
static    size_t buffer_size;


void makeData(){
    buffer_size  = N0 * N1 * N2 * 2 * sizeof(*X);
    X = (float *)malloc(buffer_size);
    printf("\nPerforming fft on an two dimensional array of size N0 x N1 x N2 : %lu x %lu x %lu\n", (unsigned long)N0, (unsigned long)N1, (unsigned long)N2);
    for (size_t i=0; i<N0; ++i) {
        for (size_t j=0; j<N1; ++j) {
            for (size_t k=0; k<N2; ++k) {
                float x = 0.0f;
                float y = 0.0f;
                if (i==0 && j==0 && k==0) { x = y = 0.5f; }
                size_t idx = 2*(k+j*N2+i*N1*N2);
                X[idx  ] = x  + (rand()&0xFF)/256.0;
                X[idx+1] = y;
                //printf("(%f, %f) ", X[idx], X[idx+1]);
            }
            //printf("\n");
        }
        //printf("\n");
    }
}

void printData( float* X ){
    // print output array
    printf("\n\nfft result: \n");
    for (size_t i=0; i<N0; ++i) {
        for (size_t j=0; j<N1; ++j) {
            for (size_t k=0; k<N2; ++k) {
                size_t idx = 2*(k+j*N2+i*N1*N2);
                printf("(%f, %f) ", X[idx], X[idx+1]);
            }
            printf("\n");
        }
        printf("\n");
    }
    printf("\n");
}

void printDataX(){ printData(X); };

void initOCL(){
    printf("DEBUG 1.0 \n");
    err = clGetPlatformIDs( 1, &platform, NULL );      printf("DEBUG 1.1 \n");
    size_t ret_param_size = 0;

    err = clGetPlatformInfo(platform, CL_PLATFORM_NAME,sizeof(platform_name), platform_name, &ret_param_size);   printf("DEBUG 1.2 \n");
    printf("Platform found: %s\n", platform_name);
    err = clGetDeviceIDs( platform, CL_DEVICE_TYPE_DEFAULT, 1, &device, NULL );                        printf("DEBUG 1.3 \n");
    err = clGetDeviceInfo(device, CL_DEVICE_NAME,sizeof(device_name), device_name,&ret_param_size);    printf("DEBUG 1.4 \n");
    printf("Device found on the above platform: %s\n", device_name);
    props[1] = (cl_context_properties)platform;
    ctx   = clCreateContext( props, 1, &device, NULL, NULL, &err );    printf("DEBUG 1.5 \n");
    queue = clCreateCommandQueue( ctx,  device, 0, &err );             printf("DEBUG 1.6 \n");
}

void initFFT(){
    /* Setup clFFT. */
    clfftSetupData fftSetup;
    err  = clfftInitSetupData(&fftSetup);
    err  = clfftSetup(&fftSetup);
    bufX = clCreateBuffer( ctx, CL_MEM_READ_WRITE, buffer_size, NULL, &err );
}

void planFFT(){
    // Create a default plan for a complex FFT. 
    err = clfftCreateDefaultPlan(&planHandle, ctx, dim, clLengths);
    err = clfftSetPlanPrecision (planHandle, CLFFT_SINGLE);
    err = clfftSetLayout        (planHandle, CLFFT_COMPLEX_INTERLEAVED, CLFFT_COMPLEX_INTERLEAVED);
    err = clfftSetResultLocation(planHandle, CLFFT_INPLACE);
    err = clfftBakePlan         (planHandle, 1, &queue, NULL, NULL);
}

void runFFT(){
    err = clEnqueueWriteBuffer ( queue, bufX, CL_TRUE, 0, buffer_size, X, 0, NULL, NULL );
    err = clfftEnqueueTransform( planHandle, CLFFT_FORWARD, 1, &queue, 0, NULL, NULL, &bufX, NULL, NULL);  // Execute the plan.                                                              
    err = clEnqueueReadBuffer  ( queue, bufX, CL_TRUE, 0, buffer_size, X, 0, NULL, NULL ); // Fetch results of calculations. 
    err = clFinish(queue);    // Wait for calculations to be finished. 
}

void cleanup(){
    clReleaseMemObject( bufX );
    free(X);
    err = clfftDestroyPlan( &planHandle );
    clfftTeardown( );
    clReleaseCommandQueue( queue );
    clReleaseContext( ctx );
}

void runAll( ){
    makeData();  printData( X );    printf( "DEBUG 1 \n" );
    initOCL();                            printf( "DEBUG 2 \n" );
    initFFT();                            printf( "DEBUG 3 \n" );
    planFFT();                            printf( "DEBUG 4 \n" );
    runFFT();    printData( X );    printf( "DEBUG 5 \n" );
    cleanup();                            printf( "DEBUG 6 \n" );
}
