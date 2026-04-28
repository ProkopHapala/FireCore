#ifndef CU_MM_H
#define CU_MM_H

#include <cuda_runtime.h>
#include <vector_types.h>
#include <string>
// #include <vector>
// #include <unordered_map>
// #include <cstdio>
// #include <cstdlib>
// #include <stdexcept>

#include "containers.h"


// --- cu_Mat3 Definition ---
// #ifndef cu_Mat3_DEFINED
// #define cu_Mat3_DEFINED
// typedef struct { float4 a; float4 b; float4 c; } cu_Mat3;
// #endif

// --- Kernel Forward Declarations ---
#include "relax_multi.cuh"

// --- Error Checking ---
#ifndef CUDA_CHECK_MM
static void HandleErrorCUMM(cudaError_t err, const char *file, int line) {
    if (err != cudaSuccess) {
        fprintf(stderr, "CUDA Error: %s in %s at line %d\n", cudaGetErrorString(err), file, line);
        exit(EXIT_FAILURE); // Exit on error
    }
}
#define CUDA_CHECK_MM(err) (HandleErrorCUMM(err, __FILE__, __LINE__))
#endif

// Contains pointer and metadata for a single GPU buffer
struct Buffer{
    void*  ptr      = nullptr;   // Raw device pointer
    size_t nBytes   = 0;         // Total allocated size in bytes
    size_t elemSize = 0;         // Size of a single element
    size_t count    = 0;         // Number of elements allocated

    Buffer() = default;
    Buffer(void* p, size_t nb, size_t es, size_t c) : ptr(p), nBytes(nb), elemSize(es), count(c) {}

    // Free the CUDA memory associated with this buffer
    void free() {
        if (ptr) {
            cudaFree(ptr);
            ptr      = nullptr;
            nBytes   = 0;
            elemSize = 0;
            count    = 0;
        }
    }

    // Templated getter for typed pointer
    template <typename T> T* get() const {
        // Optional: Add runtime check: if(sizeof(T) != elemSize && elemSize != 0) { fprintf... }
        return static_cast<T*>(ptr);
    }
};


class CU_MM { public:
    // --- Internal State ---
    bool initialized = false;

    bool bCheckErros = true;

    // Dimensions
    int nSystems = 0;
    int nAtoms = 0;
    int nnode = 0;
    int nvecs = 0;
    int npbc = 0;
    int nMaxSysNeighs = 0;

    // --- Buffer Management ---
    // Map buffer name to its Buffer struct
    // std::unordered_map<std::string, Buffer> buffers;

    Dict<Buffer> buffers;

    // --- Cached Pointers for Kernel Launches (Efficiency) ---
    // These are copies of pointers stored in the `buffers` map,
    // retrieved once after allocation for fast access in run_* methods.
    float4* d_apos = nullptr; // Using d_ prefix for clarity
    float4* d_aforce = nullptr;
    float4* d_avel = nullptr;
    float4* d_cvf = nullptr;
    float4* d_fneigh = nullptr;
    int4*   d_bkNeighs = nullptr;
    float4* d_REQs = nullptr;
    int4*   d_neighs = nullptr;
    int4*   d_neighCell = nullptr;
    float4* d_MMpars = nullptr;
    float4* d_BLs = nullptr;
    float4* d_BKs = nullptr;
    float4* d_Ksp = nullptr;
    float4* d_Kpp = nullptr;
    cu_Mat3*   d_lvecs = nullptr;
    cu_Mat3*   d_ilvecs = nullptr;
    float4* d_pbcshifts = nullptr;
    float4* d_MDpars = nullptr;
    float4* d_TDrive = nullptr;
    float4* d_constr = nullptr;
    float4* d_constrK = nullptr;
    cu_Mat3*   d_bboxes = nullptr;
    int*    d_sysneighs = nullptr;
    float4* d_sysbonds = nullptr;

    int4 nPBC_val{0, 0, 0, 0};
    float4 GFFParams_val{0.0f, 0.0f, 0.0f, 0.0f};
    int bSubtractVdW=0;


    bool isInitialized() const { return initialized; }


    void* allocateBuffer(const std::string& name, size_t count, size_t elem_size) {
        size_t nBytes = count * elem_size;
        if(nBytes<=0){ return nullptr; }
        void* d_ptr   = nullptr;
        int i = buffers.getId( name );
        if( i >= 0 ){
            printf("ERROR in CU_MM::allocateBuffer() buffer '%s' already exists.\n", name.c_str());
            exit(0);
        }
        CUDA_CHECK_MM(cudaMalloc(&d_ptr, nBytes));
        buffers.add( name, Buffer(d_ptr, nBytes, elem_size, count) );
        return d_ptr;
    }

    int uploadById(int id, const void* h_data, int nBytes, int offset=0) {
        //printf("CU_MM::uploadById(): id=%i, nBytes=%i, offset=%i\n", id, nBytes, offset);
        char* ptr = static_cast<char*>(buffers.vec[id].ptr) + offset;
        CUDA_CHECK_MM(cudaMemcpy(ptr, h_data, nBytes, cudaMemcpyHostToDevice));
        return 0;
    }

    int downloadById(int id, void* h_data, int nBytes, int offset=0) {
        char* ptr = static_cast<char*>(buffers.vec[id].ptr) + offset;
        CUDA_CHECK_MM(cudaMemcpy(h_data, ptr, nBytes, cudaMemcpyDeviceToHost));
        return 0;
    }

    void* getDevicePtrByName(const std::string& name) {
        int i = buffers.getId(name);
        if (i < 0){ printf("CU_MM::getDevicePtrByName() Error: Buffer '%s' not found.\n", name.c_str()); return nullptr; }
        return buffers.vec[i].ptr;
    }

    int uploadByName(const std::string& name, const void* h_data, int nBytes=-1, int offset=0) {
        //printf("CU_MM::uploadByName(): name=%s, nBytes=%i, offset=%i\n", name.c_str(), nBytes, offset);
        int id = buffers.getId(name);
        if (id < 0) {
            printf("ERROR in CU_MM::uploadByName() buffer '%s' not found.\n", name.c_str());
            return -1;
        }
        // If nBytes is -1, use the size from the buffer allocation
        if (nBytes < 0) nBytes = buffers.vec[id].nBytes;
        return uploadById(id, h_data, nBytes, offset);
    }

    int downloadByName(const std::string& name, void* h_data, int nBytes=-1, int offset=0) {
        int id = buffers.getId(name);
        if (id < 0){ printf("CU_MM::downloadByName() Error: Buffer '%s' not found.\n", name.c_str()); return -1; }
        // If nBytes is -1, use the size from the buffer allocation
        if (nBytes < 0) nBytes = buffers.vec[id].nBytes;
        return downloadById(id, h_data, nBytes, offset);
    }


    CU_MM() { printf("CU_MM object created (Buffer Map version).\n"); }
    ~CU_MM() { cleanup(); printf("CU_MM object destroyed (Buffer Map version).\n"); }

    // --- Initialization and Cleanup ---

    void init(int nSystems_, int nAtoms_, int nnode_, int npbc_, int nMaxSysNeighs_) {
        if (initialized) { cleanup(); }
        printf("CU_MM::init(nSys=%d, nAtoms=%d, nNode=%d, npbc=%d, nMaxSysN=%d)\n",  nSystems_, nAtoms_, nnode_, npbc_, nMaxSysNeighs_);

        nSystems = nSystems_; nAtoms = nAtoms_; nnode = nnode_;
        npbc = (npbc_ > 0) ? npbc_ : 0;
        nMaxSysNeighs = (nMaxSysNeighs_ > 0) ? nMaxSysNeighs_ : 0;
        nvecs = nAtoms + nnode;
        size_t fneigh_count = (nnode > 0) ? (size_t)nSystems * nnode * 8 : 0;

        if (nSystems <= 0 || nAtoms <= 0 || nnode < 0 ) { exit(EXIT_FAILURE); }

        // Allocate buffers using the map
        d_apos      = (float4*)allocateBuffer("apos",      (size_t)nSystems * nvecs, sizeof(float4));
        d_aforce    = (float4*)allocateBuffer("aforce",    (size_t)nSystems * nvecs, sizeof(float4));
        d_avel      = (float4*)allocateBuffer("avel",      (size_t)nSystems * nvecs, sizeof(float4));
        d_cvf       = (float4*)allocateBuffer("cvf",       (size_t)nSystems * nvecs, sizeof(float4));
        d_bkNeighs  = (int4*  )allocateBuffer("bkNeighs",  (size_t)nSystems * nvecs, sizeof(int4));
        d_fneigh    = (float4*)allocateBuffer("fneigh",    fneigh_count,             sizeof(float4));
        d_REQs      = (float4*)allocateBuffer("REQs",      (size_t)nSystems * nAtoms, sizeof(float4));
        d_neighs    = (int4*  )allocateBuffer("neighs",    (size_t)nSystems * nAtoms, sizeof(int4));
        d_neighCell = (int4*  )allocateBuffer("neighCell", (size_t)nSystems * nAtoms, sizeof(int4));
        d_constr    = (float4*)allocateBuffer("constr",    (size_t)nSystems * nAtoms, sizeof(float4));
        d_constrK   = (float4*)allocateBuffer("constrK",   (size_t)nSystems * nAtoms, sizeof(float4));
        d_MMpars    = (float4*)allocateBuffer("MMpars",    (size_t)nSystems * nnode, sizeof(float4));
        d_BLs       = (float4*)allocateBuffer("BLs",       (size_t)nSystems * nnode, sizeof(float4));
        d_BKs       = (float4*)allocateBuffer("BKs",       (size_t)nSystems * nnode, sizeof(float4));
        d_Ksp       = (float4*)allocateBuffer("Ksp",       (size_t)nSystems * nnode, sizeof(float4));
        d_Kpp       = (float4*)allocateBuffer("Kpp",       (size_t)nSystems * nnode, sizeof(float4));
        d_MDpars    = (float4*)allocateBuffer("MDpars",    (size_t)nSystems, sizeof(float4));
        d_TDrive    = (float4*)allocateBuffer("TDrive",    (size_t)nSystems, sizeof(float4));
        d_lvecs     = (cu_Mat3*  )allocateBuffer("lvecs",     (size_t)nSystems, sizeof(cu_Mat3));
        d_ilvecs    = (cu_Mat3*  )allocateBuffer("ilvecs",    (size_t)nSystems, sizeof(cu_Mat3));
        d_bboxes    = (cu_Mat3*  )allocateBuffer("bboxes",    (size_t)nSystems, sizeof(cu_Mat3));
        d_pbcshifts = (float4*)allocateBuffer("pbcshifts", (npbc > 0) ? (size_t)nSystems * npbc : 0, sizeof(float4));
        d_sysneighs = (int*   )allocateBuffer("sysneighs", (nMaxSysNeighs > 0) ? (size_t)nSystems * nMaxSysNeighs : 0, sizeof(int));
        d_sysbonds  = (float4*)allocateBuffer("sysbonds",  (nMaxSysNeighs > 0) ? (size_t)nSystems * nMaxSysNeighs : 0, sizeof(float4));

        initialized = true;
        printf("CU_MM::init() finished successfully (Buffer Map version).\n");
    }

    void cleanup() {
        if (!initialized) return;
        printf("CU_MM::cleanup() freeing %zu GPU buffers...\n", buffers.vec.size());
        for (auto& buf : buffers.vec) {
            buf.free();
        }
        buffers.vec.clear();
        buffers.map.clear();
        // Reset dimensions and cached pointers
        nSystems = nAtoms = nnode = nvecs = npbc = nMaxSysNeighs = 0;
        initialized = false;
        // Nullify cached pointers
        d_apos = d_aforce = d_avel = d_cvf = d_fneigh = nullptr;
        d_bkNeighs = nullptr; d_REQs = nullptr; d_neighs = nullptr; d_neighCell = nullptr;
        d_MMpars = d_BLs = d_BKs = d_Ksp = d_Kpp = nullptr;
        d_lvecs = d_ilvecs = nullptr; d_pbcshifts = nullptr; d_MDpars = d_TDrive = nullptr;
        d_constr = nullptr;
        d_constrK = nullptr;
        d_bboxes = nullptr;
        d_sysneighs = nullptr;
        d_sysbonds = nullptr;
        printf("CU_MM::cleanup() finished.\n");
    }

    // --- Kernel Execution Functions ---
    // Use cached d_* members for efficiency

    int run_cleanForceMMFFf4() {
        dim3 blockDim(32, 1);
        dim3 gridDim((nvecs + blockDim.x - 1) / blockDim.x, nSystems);
        int4 n_dims = {nAtoms, nnode, 0, 0};
        //printf("run_cleanForceMMFFf4(): nAtoms=%i, nnode=%i, nvecs=%i, nSystems=%i, gridDim=(%i,%i), blockDim=(%i,%i)\n", nAtoms, nnode, nvecs, nSystems, gridDim.x, gridDim.y, blockDim.x, blockDim.y);
        cleanForceMMFFf4<<<gridDim, blockDim>>>(n_dims, d_aforce, d_fneigh);
        if(bCheckErros){ CUDA_CHECK_MM(cudaGetLastError()); }
        // CUDA FINISH ?
        //CUDA_CHECK_MM(cudaDeviceSynchronize());
        return 0;
    }

    int run_getNonBond() {
        dim3 blockDim(32, 1);
        dim3 gridDim((nAtoms + blockDim.x - 1) / blockDim.x, nSystems);
        int4 ns_dims = {nAtoms, nnode, 0, (int)blockDim.x};
        size_t shared_mem_size = blockDim.x * (sizeof(float4) + sizeof(float4));
        getNonBond<<<gridDim, blockDim, shared_mem_size>>>(
            ns_dims, 
            d_apos, 
            d_aforce, 
            d_REQs, 
            d_neighs, 
            d_neighCell,
            d_lvecs, 
            nPBC_val, 
            GFFParams_val
        );
        if(bCheckErros){ CUDA_CHECK_MM(cudaGetLastError()); }
        // CUDA FINISH ?
        //CUDA_CHECK_MM(cudaDeviceSynchronize());
        return 0;
    }

    int run_getMMFFf4() {
        dim3 blockDim(32, 1);
        dim3 gridDim((nnode + blockDim.x - 1) / blockDim.x, nSystems);
        int4 ndofs_val = {nAtoms, nnode, nSystems, 0};
        getMMFFf4<<<gridDim, blockDim>>>(
            ndofs_val, 
            d_apos, 
            d_aforce, 
            d_fneigh,
            d_neighs, 
            d_neighCell, 
            d_REQs, 
            d_MMpars, 
            d_BLs, 
            d_BKs, 
            d_Ksp,
            d_Kpp,
            d_lvecs,
            d_ilvecs, 
            d_pbcshifts, 
            npbc, 
            bSubtractVdW
        );
        if(bCheckErros){ CUDA_CHECK_MM(cudaGetLastError()); }
        // CUDA FINISH ?
        //CUDA_CHECK_MM(cudaDeviceSynchronize());
        return 0;
    }

    int run_updateAtomsMMFFf4() {
        dim3 blockDim(32, 1);
        dim3 gridDim((nvecs + blockDim.x - 1) / blockDim.x, nSystems);
        int4 n_dims = {nAtoms, nnode, nSystems, nMaxSysNeighs};
        updateAtomsMMFFf4<<<gridDim, blockDim>>>(
            n_dims,
            d_apos, 
            d_avel, 
            d_aforce, 
            d_cvf, 
            d_fneigh,
            d_bkNeighs, 
            d_constr, 
            d_constrK, 
            d_MDpars, 
            d_TDrive,
            d_bboxes,
            d_sysneighs, 
            d_sysbonds
        );
        if(bCheckErros){ CUDA_CHECK_MM(cudaGetLastError()); }
        // CUDA FINISH ?
        //CUDA_CHECK_MM(cudaDeviceSynchronize());
        return 0;
    }

    int run_printOnGPU(int iSys, int4 mask) {
        dim3 blockDim(1, 1);
        dim3 gridDim(1, 1);
        int4 n_dims = {nAtoms, nnode, iSys, nSystems};
        printOnGPU<<<gridDim, blockDim>>>(
            n_dims, 
            mask, 
            d_apos, 
            d_avel, 
            d_aforce, 
            d_fneigh, 
            d_bkNeighs, 
            d_constr
        );
        CUDA_CHECK_MM(cudaGetLastError()); // Check launch
        CUDA_CHECK_MM(cudaDeviceSynchronize()); // Sync after print
        return 0;
    }

    int synchronize() {
        CUDA_CHECK_MM(cudaGetLastError());
        CUDA_CHECK_MM(cudaDeviceSynchronize());
        return 0;
    }

}; // END class CU_MM

#endif // CU_MM_H