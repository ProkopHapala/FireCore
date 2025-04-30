// mmcuda_kernels.cuh
#ifndef MMCUDA_KERNELS_CUH
#define MMCUDA_KERNELS_CUH

#include <cuda_runtime.h>
#include <vector_types.h> // For float4, int4 etc.

// Local Mat3 definition for kernels
typedef struct { float4 a; float4 b; float4 c; } Mat3;

// Declare all CUDA kernels with proper signatures
__global__ void cleanForceMMFFf4(
    int4 n,                // (natoms,nnode,?,?)
    float4* aforce,        // forces on atoms and pi-orbitals
    float4* fneigh         // recoil forces on neighbors
);
__global__ void getNonBond(
    int4 ns,               // (natoms,nnode,?,buffer elems)
    float4* atoms,         // positions (no pi-orbitals)
    float4* forces,        // accumulated forces
    float4* REQKs,         // non-bond params per atom
    int4* neighs,          // neighbor indices per atom
    int4* neighCell,       // neighbor cell indices per atom
    Mat3* lvecs,           // lattice vectors per system
    const int4 nPBC,       // PBC image counts
    const float4 GFFParams // grid-force-field params
);
__global__ void getMMFFf4(
    int4 nDOFs,            // (nAtoms,nnode,nSystems,?)
    float4* apos,          // positions incl. capping and pi-orbitals
    float4* fapos,         // forces on atoms (node atoms)
    float4* fneigh,        // recoil forces on neighbors and pi
    int4* neighs,          // neighbor indices per node atom
    int4* neighCell,       // neighbor cell index per node atom
    float4* REQKs,         // non-bond params per atom
    float4* apars,         // per-node FF params
    float4* bLs,           // bond lengths per neighbor
    float4* bKs,           // bond stiffness per neighbor
    float4* Ksp,           // sigma-pi orth stiffness
    float4* Kpp,           // pi-pi planar stf
    Mat3* lvecs,           // lattice vectors per system
    Mat3* ilvecs,          // inverse lattice vectors per system
    float4* pbc_shifts,    // precomputed PBC shifts
    const int npbc,        // count of PBC shifts
    const int bSubtractVdW // flag to subtract bonded nonbonded
);
__global__ void updateAtomsMMFFf4(
    int4 n,                // (natoms,nnode,nsys,nMaxSysNeighs)
    float4* apos,          // positions incl. pi-orbitals
    float4* avel,          // velocities
    float4* aforce,        // forces on atoms
    float4* cvf,           // accum |f|^2,|v|^2,<f|v>
    float4* fneigh,        // recoil forces
    int4* bkNeighs,        // back neighbor indices
    float4* constr,        // constraints per atom
    float4* constrK,       // constraint stiffness per atom
    float4* MDparams,      // MD params dt,damp, etc.
    float4* TDrives,       // thermal drive params
    Mat3* bboxes,          // bounding boxes per system
    int* sysneighs,        // sys neighbor lists
    float4* sysbonds       // system-system bond params
);
__global__ void printOnGPU(
    int4 n,                // (natoms,nnode,sysToPrint,?)
    int4 mask,             // print flags
    float4* apos,          // positions
    float4* avel,          // velocities
    float4* aforce,        // forces
    float4* fneigh,        // recoil forces
    int4* bkNeighs,        // back neighbor indices
    float4* constr         // constraints
);

#endif // MMCUDA_KERNELS_CUH