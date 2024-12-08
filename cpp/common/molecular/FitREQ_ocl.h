#ifndef FitREQ_ocl_h
#define FitREQ_ocl_h

#include <vector>
#include "Vec3.h"
#include "quaternion.h"
//#include "NBFF.h"
#include "Atoms.h"
#include "MMFFparams.h"
//#include "Forces.h"
//#include "MMFFBuilder.h"
//#include "functions.h"
#include "FitREQ.h"
#include "OCLfft_errors.h"
#include "OCL.h"

class FitREQ_ocl : public FitREQ { public:
    // OpenCL context and management
    OCLsystem* cl = nullptr;
    
    // OpenCL buffers for kernel execution
    cl_mem cl_ranges    = nullptr;  // int4:    (i0,ni, j0,nj) start and number of atoms in fragments 1,2
    cl_mem cl_tREQHs    = nullptr;  // float4:  (RvdW,EvdW,Q,Hbond) parameters for each type
    cl_mem cl_atypes    = nullptr;  // int:     atom type indices
    cl_mem cl_ieps      = nullptr;  // int2:    electron pair indices
    cl_mem cl_atoms     = nullptr;  // float4:  (x,y,z,Q) atom positions and charges
    cl_mem cl_dEdREQs   = nullptr;  // float4:  derivatives output
    cl_mem cl_fDOFs     = nullptr;  // float:   DOF derivatives
    cl_mem cl_DOFtoAtom = nullptr; // int:     mapping from DOFs to atoms
    cl_mem cl_DOFnis    = nullptr;  // int2:    (i0,ni) ranges for DOF atom mappings
    
    //const int4 ns,                  
    // __global int4*    ranges,       // [nsample]  (i0,ni, j0,nj) star and number of atoms in fragments 1,2 
    // __global float4*  tREQHs,       // [ntypes]   non-bonded parameters (RvdW,EvdW,QvdW,Hbond)
    // __global int*     atypes,       // [natomTot] atom types
    // //__global int*     hosts,      // [natomTot] index of host atoms (it electron pair), -1 if it is not electron pair
    // __global int2*    ieps,         // [natomTot] {iep1,iep2} index of electron pair type ( used to subtract  charge)
    // __global float4*  atoms,        // [natomTot] positions and REPS charge of each atom (x,y,z,Q) 
    // __global float4*  dEdREQs       // [natomTot] output derivatives of type REQH parameters

    // __global float*   fDOFs,      // [nDOFs]    derivatives of REQH parameters
    // __global int2*    DOFnis,     // [nDOFs]    (i0,ni) star and number of atoms in fragments 1,2    
    // __global int*     DOFtoAtom,  // [nInds]    list of atom indexes relevant for each DOF (non-uniform ranges indexed by DOFnis)
    // __global int*     DOFcofefs,  // [nInds]    factors for update of each DOF from the atom dEdREQH parameters   fDOFi = dot( DOFcofefs[i], dEdREQH[i] )
    // __global float4*  dEdREQs     // [nAtomTot] derivatives of REQH parameters

    // Temporary CPU-side arrays for buffer preparation
    std::vector<Quat4i>  ranges;    // [nsample]  (i0,ni, j0,nj) star and number of atoms in fragments 1,2 
    std::vector<int>     atypes;    // [natomTot]  atom types
    std::vector<Vec2i>   ieps;      // [natomTot] {iep1,iep2} index of electron pair type ( used to subtract  charge)
    std::vector<Quat4f>  atoms;     // [natomTot] positions and REPS charge of each atom (x,y,z,Q)
    std::vector<Quat4f>  dEdREQs;   // [natomTot] output derivatives of type REQH parameters
    std::vector<float>   fDOFs_;     // [nDOFs]    derivatives of REQH parameters
    std::vector<Vec2i>   DOFnis;    // [nDOFs] Ranges for DOF mappings
    std::vector<int>     DOFtoAtom; // [nInds]  list of atom indexes relevant for each DOF
    std::vector<Quat4f>  DOFcofefs; // [nInds] Factors for update of each DOF from the atom dEdREQH parameters   fDOFi = dot( DOFcofefs[i], dEdREQH[i] )
    
    // Buffer sizes and counts
    int natomsTot = 0;     // Total number of atoms across all samples
    int nsamples  = 0;      // Number of samples
    int nDOFsTot  = 0;      // Total number of DOF indices
    //int nDOFs     = 0;      // Number of DOFs [from 
        
// ============= Functions =============

void export_range(Atoms* sample, int offset, Quat4i& range) {
    int n1 = sample->n0;           // Number of atoms in fragment 1
    int n2 = sample->natoms - n1;  // Number of atoms in fragment 2
    range = Quat4i{
        offset,     // i0: Start of fragment 1
        n1,         // ni: Number of atoms in fragment 1
        offset+n1,  // j0: Start of fragment 2
        n2          // nj: Number of atoms in fragment 2
    };
}

void export_epairs(Atoms* sample, int offset) {
    // First set all to no electron pairs
    int na = sample->natoms;
    for(int i=0; i<na; i++){  ieps[offset+i]=Vec2i{-1,-1}; }
    // Then process electron pairs if present
    AddedData* adata = (AddedData*)sample->userData;
    if(adata){
        for(int i=0; i<adata->nep; i++) {
            int ih = adata->bs[i].x;
            int ie = adata->bs[i].y;
            Vec2i& epi = ieps[offset+ih];
            if(epi.x == -1){ epi.x = ie; }
            else           { epi.y = ie; }
        }
    }
}

void prepare_samples() {
    // Count total atoms and samples
    nsamples = samples.size();
    natomsTot = 0;
    for(auto* sample : samples) {  natomsTot += sample->natoms; }
    // Resize all arrays
    ranges.resize(nsamples);
    atypes.resize(natomsTot);
    ieps  .resize(natomsTot);
    atoms .resize(natomsTot);
    // Process all samples in a single loop
    int offset = 0;
    for(int i=0; i<nsamples; i++) {
        Atoms* sample = samples[i];
        for(int j=0; j<sample->natoms; j++) {
            int ia = offset + j; 
            atypes[ia]   = sample->atypes[j];
            atoms [ia].f = (Vec3f)sample->apos[j];
            atoms [ia].w = sample->charge[j];
        }
        export_range (sample, offset, ranges[i]);
        export_epairs(sample, offset);
        offset += sample->natoms;
    }
}

int processDOFmappings(bool bWrite=false) {
    int nfound = 0;
    std::vector<int> nextIndex;
    if(bWrite){
        nextIndex.resize(nDOFs, 0);
    }
    // Process each DOF
    for(int idof=0; idof<nDOFs; idof++) {
        Vec2i typeComp = REQtoTyp[idof];
        int   type     = typeComp.x;
        int   icomp    = typeComp.y;
        
        // Count atoms of this type
        int count = 0;
        for(int i=0; i<natomsTot; i++) {
            if(atypes[i] == type){
                if(bWrite){
                    int baseIdx = DOFnis[idof].x;
                    int idx     = nextIndex[idof]++;
                    DOFtoAtom [baseIdx + idx] = i;
                    // Set coefficient vector (like {0,0,1,0} for Q component)
                    Quat4f coef = Quat4fZero;
                    coef.array[icomp] = 1.0f;
                    DOFcofefs[baseIdx + idx] = coef;
                }
                count++;
            }
        }
        if(bWrite){
            DOFnis[idof] = Vec2i{nfound, count};
        }
        nfound += count;
    }
    return nfound;
}

void prepareDOFMappings() {
    // First pass: count total DOF indices
    nDOFsTot = processDOFmappings(false);
    
    // Allocate arrays
    DOFtoAtom.resize(nDOFsTot);
    DOFcofefs.resize(nDOFsTot);
    DOFnis.resize(nDOFs);
    
    // Second pass: build mappings
    processDOFmappings(true);
}

/*
void allocBuffers() {
    // Calculate total sizes
    nsamples = samples.size();
    natomsTot = 0;
    for(auto* sample : samples) {
        natomsTot += sample->natoms;
    }
    // Allocate OpenCL buffers
    cl_ranges    = clCreateBuffer(cl->context, CL_MEM_READ_ONLY,  sizeof(cl_int4) * nsamples, NULL, NULL);
    cl_tREQHs    = clCreateBuffer(cl->context, CL_MEM_READ_ONLY,  sizeof(cl_float4) * ntype, NULL, NULL);
    cl_atypes    = clCreateBuffer(cl->context, CL_MEM_READ_ONLY,  sizeof(cl_int) * natomsTot, NULL, NULL);
    cl_ieps      = clCreateBuffer(cl->context, CL_MEM_READ_ONLY,  sizeof(cl_int2) * natomsTot, NULL, NULL);
    cl_atoms     = clCreateBuffer(cl->context, CL_MEM_READ_ONLY,  sizeof(cl_float4) * natomsTot, NULL, NULL);
    cl_dEdREQs   = clCreateBuffer(cl->context, CL_MEM_WRITE_ONLY, sizeof(cl_float4) * natomsTot, NULL, NULL);
    cl_fDOFs     = clCreateBuffer(cl->context, CL_MEM_WRITE_ONLY, sizeof(cl_float) * nDOFs, NULL, NULL);
    cl_DOFtoAtom = clCreateBuffer(cl->context, CL_MEM_READ_ONLY,  sizeof(cl_int) * nDOFsTot, NULL, NULL);
    cl_DOFnis    = clCreateBuffer(cl->context, CL_MEM_READ_ONLY,  sizeof(cl_int2) * nDOFs, NULL, NULL);
}

void uploadBuffers() {
    // Upload all prepared data to GPU
    clEnqueueWriteBuffer(cl->queue, cl_ranges,    CL_TRUE, 0, sizeof(cl_int4) * nsamples, ranges.data(), 0, NULL, NULL);
    clEnqueueWriteBuffer(cl->queue, cl_tREQHs,    CL_TRUE, 0, sizeof(cl_float4) * ntype, typeREQs, 0, NULL, NULL);
    clEnqueueWriteBuffer(cl->queue, cl_atypes,    CL_TRUE, 0, sizeof(cl_int) * natomsTot, atypes.data(), 0, NULL, NULL);
    clEnqueueWriteBuffer(cl->queue, cl_ieps,      CL_TRUE, 0, sizeof(cl_int2) * natomsTot, ieps.data(), 0, NULL, NULL);
    clEnqueueWriteBuffer(cl->queue, cl_atoms,     CL_TRUE, 0, sizeof(cl_float4) * natomsTot, atoms.data(), 0, NULL, NULL);
    clEnqueueWriteBuffer(cl->queue, cl_DOFtoAtom, CL_TRUE, 0, sizeof(cl_int) * nDOFsTot, DOFtoAtom.data(), 0, NULL, NULL);
    clEnqueueWriteBuffer(cl->queue, cl_DOFnis,    CL_TRUE, 0, sizeof(cl_int2) * nDOFs, DOFnis.data(), 0, NULL, NULL);
}

void downloadBuffers() {
    // Download results from GPU
    clEnqueueReadBuffer(cl->queue, cl_dEdREQs, CL_TRUE, 0, sizeof(cl_float4) * natomsTot, dEdREQs.data(), 0, NULL, NULL);
    clEnqueueReadBuffer(cl->queue, cl_fDOFs,   CL_TRUE, 0, sizeof(cl_float) * nDOFs, fDOFs, 0, NULL, NULL);
}


void evalSampleDerivatives() {
    // Prepare all data
    prepare_samples();
    prepareDOFMappings();
    
    // Allocate and upload buffers
    allocBuffers();
    uploadBuffers();
    
    // Set kernel arguments and execute
    cl_kernel kernel = cl->getKernel("evalSampleDerivatives");
    clSetKernelArg(kernel, 0, sizeof(cl_mem), &cl_ranges);
    clSetKernelArg(kernel, 1, sizeof(cl_mem), &cl_tREQHs);
    clSetKernelArg(kernel, 2, sizeof(cl_mem), &cl_atypes);
    clSetKernelArg(kernel, 3, sizeof(cl_mem), &cl_ieps);
    clSetKernelArg(kernel, 4, sizeof(cl_mem), &cl_atoms);
    clSetKernelArg(kernel, 5, sizeof(cl_mem), &cl_dEdREQs);
    
    // Launch kernel with appropriate work size
    size_t globalWorkSize = natomsTot;
    size_t localWorkSize = 32;  // As specified in kernel
    clEnqueueNDRangeKernel(cl->queue, kernel, 1, NULL, &globalWorkSize, &localWorkSize, 0, NULL, NULL);
    
    // Download results
    downloadBuffers();
}
*/

}; // class FitREQ_ocl

#endif
