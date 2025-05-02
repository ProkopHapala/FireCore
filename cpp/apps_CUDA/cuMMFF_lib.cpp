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
    int uploadId  (int id, const void* h_data, int nbyte) { return ff.uploadById  (id, h_data, nbyte); }
    int downloadId(int id, void*       h_data, int nbyte) { return ff.downloadById(id, h_data, nbyte); }

    int upload(const char* name, const void* h_data, int nbyte) {
        return ff.uploadByName(std::string(name), h_data);
    }
    int download(const char* name, void* h_data, int nbyte ) {
        return ff.downloadByName(std::string(name), h_data);
    }

    void init(int nSystems_, int nAtoms_, int nnode_, int npbc_, int nMaxSysNeighs_) {
        ff.init(nSystems_, nAtoms_, nnode_, npbc_, nMaxSysNeighs_);
    }


    // --- Kernel Execution Functions ---
    // Call the public methods of the CU_MM instance
    int run_cleanForceMMFFf4() {
        return ff.run_cleanForceMMFFf4();
    }
    int run_getNonBond(int* h_nPBC, float* h_GFFparams) {
        return ff.run_getNonBond( *(int4*)h_nPBC, *(float4*)h_GFFparams);
    }
    int run_getMMFFf4(int bSubtractVdW) {
        return ff.run_getMMFFf4(bSubtractVdW);
    }
    int run_updateAtomsMMFFf4() {
        return ff.run_updateAtomsMMFFf4();
    }
    int run_printOnGPU(int iSys, int* mask) {
        return ff.run_printOnGPU(iSys, *(int4*)mask);
    }
    int synchronize() {
        return ff.synchronize();
    }

} // extern "C"