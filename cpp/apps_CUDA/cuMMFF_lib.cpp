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

    // --- Generic Upload/Download by Name ---
    int mmcuda_upload(const char* name, const void* h_data, size_t elem_size, size_t count) {
        // Use CU_MM's map-based upload for convenience
        return ff.uploadByName(std::string(name), h_data);
        // We need count/elem_size if uploadByName needs them for size check
        // Assuming uploadByName calculates size internally based on metadata stored during init
        // If not, the C interface needs count/elem_size:
        // CudaBuffer bufInfo = ff.getBufferInfo(std::string(name));
        // if (bufInfo.elemSize != elem_size || bufInfo.count != count ) { /* size mismatch error */ return 1; }
        // return ff.uploadByName(std::string(name), h_data);
    }
    int mmcuda_download(const char* name, void* h_data, size_t elem_size, size_t count) {
         // Similar logic as upload
        return ff.downloadByName(std::string(name), h_data);
    }


    // --- Kernel Execution Functions ---
    // Call the public methods of the CU_MM instance
    int mmcuda_run_cleanForceMMFFf4() {
        return ff.run_cleanForceMMFFf4();
    }
    int mmcuda_run_getNonBond(int4 h_nPBC, float4 h_GFFparams) {
        return ff.run_getNonBond(h_nPBC, h_GFFparams);
    }
    int mmcuda_run_getMMFFf4(int bSubtractVdW) {
        return ff.run_getMMFFf4(bSubtractVdW);
    }
    int mmcuda_run_updateAtomsMMFFf4() {
        return ff.run_updateAtomsMMFFf4();
    }
    int mmcuda_run_printOnGPU(int iSys, int4 mask) {
        return ff.run_printOnGPU(iSys, mask);
    }
    int mmcuda_synchronize() {
        return ff.synchronize();
    }

} // extern "C"