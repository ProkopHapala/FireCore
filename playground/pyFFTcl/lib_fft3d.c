
#include <stdio.h>
#include <stdlib.h>

// No need to explicitely include the OpenCL headers 
//#include <clFFT.h>
#include "lib_fft3.h"

static int ret = 0;
static cl_int err;
static char platform_name[128];
static char device_name  [128];
static cl_platform_id platform = 0;
static cl_device_id device     = 0;
static cl_context_properties props[3] = { CL_CONTEXT_PLATFORM, 0, 0 };
static cl_context          ctx   = 0;
static cl_command_queue    queue = 0;

static cl_event event = NULL;
static clfftPlanHandle planHandle;
static clfftDim dim = CLFFT_3D;

static const size_t N0 = 4, N1 = 4, N2 = 4;
static size_t clLengths[3] = {N0, N1, N2};
static size_t buffer_size;
//static int clLengths[3] = {N0, N1, N2};
//static int buffer_size;

static cl_mem data_cl;
static float *data;


int N=0;
int size = N * N;
static cl_mem d_a, d_b, d_c;


// ========= GRENERAL OpenCL ========== 
static cl_program       program;       // compute program
static cl_kernel        kernel;        // compute kernel
static KernelDims kdim;

void makeData(){
    buffer_size  = N0 * N1 * N2 * 2 * sizeof(*data);
    data = (float *)malloc(buffer_size);
    printf("\nPerforming fft on an two dimensional array of size N0 x N1 x N2 : %lu x %lu x %lu\n", (unsigned long)N0, (unsigned long)N1, (unsigned long)N2);
    for (size_t i=0; i<N0; ++i) {
        for (size_t j=0; j<N1; ++j) {
            for (size_t k=0; k<N2; ++k) {
                float x = 0.0f;
                float y = 0.0f;
                if (i==0 && j==0 && k==0) { x = y = 0.5f; }
                size_t idx = 2*(k+j*N2+i*N1*N2);
                data[idx  ] = x  + (rand()&0xFF)/256.0;
                data[idx+1] = y;
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

void printDataX(){ printData(data); };

void initOCL(){
    //printf("DEBUG 1.0 \n");
    err = clGetPlatformIDs( 1, &platform, NULL );      //printf("DEBUG 1.1 \n");
    size_t ret_param_size = 0;

    err = clGetPlatformInfo(platform, CL_PLATFORM_NAME,sizeof(platform_name), platform_name, &ret_param_size);   //printf("DEBUG 1.2 \n");
    printf("Platform found: %s\n", platform_name);
    err = clGetDeviceIDs( platform, CL_DEVICE_TYPE_DEFAULT, 1, &device, NULL );                        //printf("DEBUG 1.3 \n");
    err = clGetDeviceInfo(device, CL_DEVICE_NAME,sizeof(device_name), device_name,&ret_param_size);    //printf("DEBUG 1.4 \n");
    printf("Device found on the above platform: %s\n", device_name);
    props[1] = (cl_context_properties)platform;
    ctx   = clCreateContext( props, 1, &device, NULL, NULL, &err );    //printf("DEBUG 1.5 \n");
    queue = clCreateCommandQueue( ctx,  device, 0, &err );             //printf("DEBUG 1.6 \n");
}


void planFFT(){
    // Create a default plan for a complex FFT. 
    err = clfftCreateDefaultPlan(&planHandle, ctx, dim, clLengths);
    err = clfftSetPlanPrecision (planHandle, CLFFT_SINGLE);
    err = clfftSetLayout        (planHandle, CLFFT_COMPLEX_INTERLEAVED, CLFFT_COMPLEX_INTERLEAVED);
    err = clfftSetResultLocation(planHandle, CLFFT_INPLACE);
    err = clfftBakePlan         (planHandle, 1, &queue, NULL, NULL);
}


void initFFT( int ndim, size_t* Ns ){
    //printf("DEBUG initFFT() ndim %i Ns[%li,%li,%li]\n", ndim, Ns[0], Ns[1], Ns[2] );
    if     (ndim==1){ dim=CLFFT_1D; }
    else if(ndim==2){ dim=CLFFT_2D; }
    else if(ndim==3){ dim=CLFFT_3D; };
    clLengths[0] = Ns[0];
    clLengths[1] = Ns[1];
    clLengths[2] = Ns[2];
    buffer_size  = 2 * sizeof(*data);
    for(int i=0; i<ndim;i++){ buffer_size *= Ns[i]; }

    clfftSetupData fftSetup;
    err  = clfftInitSetupData(&fftSetup);
    err  = clfftSetup(&fftSetup);
    data_cl = clCreateBuffer( ctx, CL_MEM_READ_WRITE, buffer_size, NULL, &err );
    planFFT(  );
}

void loadData( float* data_ ){
    //printf("DEBUG loadData() \n");
    data = data_;
    //makeData();
    //printData( data ); 
    err = clEnqueueWriteBuffer  ( queue, data_cl, CL_TRUE, 0, buffer_size, data, 0, NULL, NULL );
}

void runFFT(){
    //err = clEnqueueWriteBuffer ( queue, data_cl, CL_TRUE, 0, buffer_size, data, 0, NULL, NULL );
    err = clfftEnqueueTransform( planHandle, CLFFT_FORWARD, 1, &queue, 0, NULL, NULL, &data_cl, NULL, NULL);  // Execute the plan. -> Forward Transform
    //err = clfftEnqueueTransform( planHandle, CLFFT_BACKWARD, 1, &queue, 0, NULL, NULL, &data_cl, NULL, NULL);  // Execute the plan.  -> Backward Transform                                                        
    err = clEnqueueReadBuffer  ( queue, data_cl, CL_TRUE, 0, buffer_size, data, 0, NULL, NULL ); // Fetch results of calculations. 
    err = clFinish(queue);    // Wait for calculations to be finished. 
    //printData( data ); 
}

void cleanup(){
    clReleaseMemObject( data_cl );
    //free( data );
    err = clfftDestroyPlan( &planHandle );
    clfftTeardown( );
    clReleaseCommandQueue( queue );
    clReleaseContext( ctx );
}

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




// ======================================================
// ==================== GENERAL OPENCL
// ======================================================
// CODE FROM /home/prokop/Dropbox/MyDevSW/OpenCL/OpenCL_in_C/MatrixMutl

char * getKernelSource(char *filename){
    FILE *file = fopen(filename, "r");
    if (!file){
        printf( "Error: Could not open kernel source file\n");
        exit(EXIT_FAILURE);
    }
    fseek(file, 0, SEEK_END);
    int len = ftell(file) + 1;
    rewind(file);

    char *source = (char *)calloc(sizeof(char), len);
    if (!source){
        printf( "Error: Could not allocate memory for source string\n");
        exit(EXIT_FAILURE);
    }
    fread(source, sizeof(char), len, file);
    fclose(file);
    return source;
}

int makeKernel( char * fname, char * name ){
    char * kernelsource = getKernelSource( fname);
    // Create the comput program from the source buffer
    program = clCreateProgramWithSource(ctx, 1, (const char **) & kernelsource, NULL, &err);
    char tmpstr[1024];
    sprintf(tmpstr,"Creating program with %s", fname);
    checkError(err, tmpstr);
    free(kernelsource);
    // Build the program
    err = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
    if (err != CL_SUCCESS){
        size_t len;
        char buffer[2048];
        printf("Error: Failed to build program executable!\n%s\n", err_code(err));
        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
        printf("%s\n", buffer);
        return EXIT_FAILURE;
    }
    // Create the compute kernel from the program
    kernel = clCreateKernel(program, name, &err);
    checkError(err, tmpstr);
}


int runKernel( int kind, KernelDims kdim ){

    //zero_mat(N, h_C);
    //setArray      ( N*N,  0.0,   h_C );
    //cl_uint  work_dim  = 1;
    //size_t   global[2] = {N,N};
    //size_t   local [2] = {N,N};

    err =  clSetKernelArg(kernel, 0, sizeof(int),    &N   );
    err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &d_a );
    err |= clSetKernelArg(kernel, 2, sizeof(cl_mem), &d_b );
    err |= clSetKernelArg(kernel, 3, sizeof(cl_mem), &d_c );
    if       (kind==1){
        err |= clSetKernelArg(kernel, 4, sizeof(float) * N, NULL);
    }else if (kind==2){
        err |= clSetKernelArg(kernel, 4, sizeof(float) * kdim.blocksize * kdim.blocksize, NULL);
        err |= clSetKernelArg(kernel, 5, sizeof(float) * kdim.blocksize * kdim.blocksize, NULL);
    }
    checkError(err, "Setting kernel args");

    double start_time = wtime();
    err = clEnqueueNDRangeKernel(  queue,   kernel,   kdim.dim, NULL, kdim.global, kdim.local, 0, NULL, NULL);    
    checkError(err, "Enqueueing kernel");

    err = clFinish(queue);
    checkError(err, "Waiting for kernel to finish");
    double run_time = wtime() - start_time;

    //err = clEnqueueReadBuffer(  queue, d_c, CL_TRUE, 0, sizeof(float) * size, h_C,  0, NULL, NULL);
    //checkError(err, "Reading back d_c");
    //results(N, h_C, run_time);
}


void runCLKernel(){
    makeKernel( "cl/C_block_form.cl", "mmul" );
    kdim.dim        = 2;
    kdim.blocksize  = 16;
    kdim.global[0]  = N;              kdim.global[1]  = N;
    kdim.local [0]  = kdim.blocksize; kdim.local [1]  = kdim.blocksize;
    runKernel( 2, kdim );
}

void releasOCL(){
    clReleaseMemObject(d_a);
    clReleaseMemObject(d_b);
    clReleaseMemObject(d_c);
    clReleaseProgram(program);
    clReleaseKernel(kernel);
    clReleaseCommandQueue(queue);
    clReleaseContext(ctx);
}