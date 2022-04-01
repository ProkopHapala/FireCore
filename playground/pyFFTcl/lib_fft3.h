#include <clFFT.h>

typedef struct{
    cl_uint  dim;
    size_t * global;
    size_t * local;
    cl_int   blocksize;
} KernelDims;


void makeData();
void printData( float* X );
void printDataX();
void initOCL();
void initFFT();
void planFFT();
void runFFT();
void cleanup();

//extern "C" {
//    void runAll( );
//}