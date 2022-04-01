/* ************************************************************************
 * Copyright 2013 Advanced Micro Devices, Inc.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * ************************************************************************/

#include <stdio.h>
#include <stdlib.h>

/* No need to explicitely include the OpenCL headers */
//#include <clFFT.h>
#include "lib_fft3.h"

//FFTcl fftcl;

int main( void ){
    makeData();  printData( X );    printf( "DEBUG 1 \n" );
    initOCL();                                  printf( "DEBUG 2 \n" );
    initFFT();                                  printf( "DEBUG 3 \n" );
    planFFT();                                  printf( "DEBUG 4 \n" );
    runFFT();    printData( X );    printf( "DEBUG 5 \n" );
    cleanup();                                  printf( "DEBUG 6 \n" );
    return 0;
}
