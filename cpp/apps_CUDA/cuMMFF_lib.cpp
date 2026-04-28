// File: mmcuda_lib_c_interface.cu

//#include "relax_multi.cuh" // Public C API header (using CUDA types)
#include "cuMM.h"            // The C++ CUDA implementation class

// vector_types.h is included by cuMM.h via relax_multi.cuh

#include <cstdio>             // For error messages (keep minimal)
#include <string>             // For using std::string names

// --- Global Instance ---
static CU_MM ff;

// --- C API Implementation ---
extern "C" {

    // ... mmcuda_init, mmcuda_cleanup, mmcuda_isInitialized ...

    int getBufferIndex(const char* name) { return ff.buffers.getId(name); }
    int uploadId  (int id, const void* h_data, int nbyte, int offset=0) { return ff.uploadById  (id, h_data, nbyte, offset); }
    int downloadId(int id, void*       h_data, int nbyte, int offset=0) { return ff.downloadById(id, h_data, nbyte, offset); }

    int upload(const char* name, const void* h_data, int nbyte=-1, int offset=0) {
        return ff.uploadByName(std::string(name), h_data, nbyte, offset);
    }
    int download(const char* name, void* h_data, int nbyte=-1, int offset=0 ) {
        return ff.downloadByName(std::string(name), h_data, nbyte, offset);
    }

    void init(int nSystems_, int nAtoms_, int nnode_, int npbc_, int nMaxSysNeighs_) {
        ff.init(nSystems_, nAtoms_, nnode_, npbc_, nMaxSysNeighs_);
    }

    // --- Kernel Execution Functions ---
    // Call the public methods of the CU_MM instance
    int run_cleanForceMMFFf4 ()                                  { return ff.run_cleanForceMMFFf4(); }
    int run_getNonBond       ()                                  { return ff.run_getNonBond(); }
    int run_getMMFFf4        ()                                  { return ff.run_getMMFFf4(); }
    int run_updateAtomsMMFFf4()                                  { return ff.run_updateAtomsMMFFf4(); }
    int run_printOnGPU       (int iSys, int* mask)               { return ff.run_printOnGPU(iSys, *(int4*)mask); }
    int synchronize          ()                                  { 
        printf("ff.synchronize()\n");
        return ff.synchronize(); 
    }

    int run_MD(int nstep, int mask){
        printf("cuMMFF_lib.cpp run_MD(): nstep %i\n", nstep);
        for(int i=0; i<nstep; i++){
            //printf("cuMMFF_lib.cpp run_MD(): step %i\n", i);
            //run_cleanForceMMFFf4();
            if(mask&1)run_getNonBond();
            if(mask&2)run_getMMFFf4();
            if(mask&4)run_updateAtomsMMFFf4();
        }
        CUDA_CHECK_MM(cudaDeviceSynchronize());
        return 0;
    }

} // extern "C"