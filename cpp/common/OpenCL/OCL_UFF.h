#ifndef OCL_UFF_h
#define OCL_UFF_h

/*
This is an OpenCL interface for the Universal Force Field (UFF).
It is designed to simulate multiple replicas of a system in parallel on a GPU
to maximize hardware utilization, especially for small molecular systems.

The design follows the pattern established by OCL_MM.h, with separate kernels
for each UFF interaction component (bonds, angles, dihedrals, inversions).
Forces are calculated in pieces and stored in a temporary buffer (`fint`) to
avoid race conditions, and then assembled in a final step.
*/

#include "OCL.h"

class OCL_UFF : public OCLsystem {
public:
    // --- Dimensions
    int nSystems    = 0;
    int nAtoms      = 0;
    int nBonds      = 0;
    int nAngles     = 0;
    int nDihedrals  = 0;
    int nInversions = 0;
    int nPBC        = 0;
    int nAtomsTot   = 0;
    int nA2F        = 0;

    // --- Component Flags
    bool bUFF_bonds      = true;
    bool bUFF_angles     = true;
    bool bUFF_dihedrals  = true;
    bool bUFF_inversions = true;
    bool bSubtractNB     = true; // Subtract 1-3 and 1-4 non-bonded interactions
    bool bClampNonBonded = true;

    // --- OpenCL Buffers (identified by integer index)
    // Particle buffers
    int ibuff_apos    = -1; // Atom positions {x,y,z,w}
    int ibuff_fapos   = -1; // Atom forces    {fx,fy,fz,E}
    int ibuff_REQs    = -1; // Non-bonded parameters {R,E,Q,H}

    // Auxiliary buffers (computed on GPU or temporary)
    int ibuff_hneigh  = -1; // Normalized bond vectors {hx,hy,hz,1/L}
    int ibuff_fint    = -1; // Intermediate force pieces

    // Topology buffers
    int ibuff_neighs    = -1; // Neighbor atom indices
    int ibuff_neighCell = -1; // Neighbor cell indices for PBC
    int ibuff_neighBs   = -1; // Neighbor bond indices

    // Bond buffers
    int ibuff_bonAtoms  = -1; // Bond atom indices {i,j}
    int ibuff_bonParams = -1; // Bond parameters {k, l0}

    // Angle buffers
    int ibuff_angAtoms   = -1; // Angle atom indices {i,j,k}
    int ibuff_angNgs     = -1; // Precomputed angle neighbor indices in hneigh
    int ibuff_angParams1 = -1; // Angle params {k, c0, c1, c2}
    int ibuff_angParams2_w = -1; // Angle param {c3}

    // Dihedral buffers
    int ibuff_dihAtoms  = -1; // Dihedral atom indices {i,j,k,l}
    int ibuff_dihNgs    = -1; // Precomputed dihedral neighbor indices in hneigh
    int ibuff_dihParams = -1; // Dihedral parameters {V, d, n}

    // Inversion buffers
    int ibuff_invAtoms  = -1; // Inversion atom indices {i,j,k,l}
    int ibuff_invNgs    = -1; // Precomputed inversion neighbor indices in hneigh
    int ibuff_invParams = -1; // Inversion parameters {K, c0, c1, c2}

    // Atom-to-Force mapping buffers for assembly
    int ibuff_a2f_offsets = -1;
    int ibuff_a2f_counts  = -1;
    int ibuff_a2f_indices = -1;

    // System-wide buffers
    int ibuff_pbcshifts = -1; // PBC shift vectors
    int ibuff_lvecs     = -1; // Lattice vectors for each system
    int ibuff_energies  = -1; // Buffer to store computed energies

    // --- OpenCL Kernels (identified by OCLtask object)
    OCLtask* task_evalBonds     = nullptr;
    OCLtask* task_evalAngles    = nullptr;
    OCLtask* task_evalDihedrals = nullptr;
    OCLtask* task_evalInversions= nullptr;
    OCLtask* task_assemble      = nullptr;

    // ====================== Functions

    void makeKernels(const char* cl_src_dir) {
        char srcpath[1024];
        sprintf(srcpath, "%s/UFF.cl", cl_src_dir);
        buildProgram(srcpath, program); // Assuming 'program' is the member from OCL base class

        // Create tasks for each kernel
        // TODO: The local work-group sizes (e.g., 32) are hardcoded for now. They should be tuned for optimal performance based on the device and kernel characteristics.
        //                                name                  program  nL nG
        newTask("evalBondsAndHNeigh_UFF", program, 1, (size_t4){0,0,0,0}, (size_t4){32,1,1,1} );
        newTask("evalAngles_UFF",         program, 1, (size_t4){0,0,0,0}, (size_t4){32,1,1,1} );
        newTask("evalDihedrals_UFF",      program, 1, (size_t4){0,0,0,0}, (size_t4){32,1,1,1} );
        newTask("evalInversions_UFF",     program, 1, (size_t4){0,0,0,0}, (size_t4){32,1,1,1} );
        newTask("assembleForces_UFF",     program, 1, (size_t4){0,0,0,0}, (size_t4){32,1,1,1} );

        // Get task pointers
        task_evalBonds      = getTask("evalBondsAndHNeigh_UFF");
        task_evalAngles     = getTask("evalAngles_UFF");
        task_evalDihedrals  = getTask("evalDihedrals_UFF");
        task_evalInversions = getTask("evalInversions_UFF");
        task_assemble       = getTask("assembleForces_UFF");
    }

    void realloc(int nSystems_, int nAtoms_, int nBonds_, int nAngles_, int nDihedrals_, int nInversions_, int nPBC_, int nA2F_) {
        nSystems    = nSystems_;
        nAtoms      = nAtoms_;
        nBonds      = nBonds_;
        nAngles     = nAngles_;
        nDihedrals  = nDihedrals_;
        nInversions = nInversions_;
        nPBC        = nPBC_;
        nA2F        = nA2F_;
        nAtomsTot   = nSystems * nAtoms;

        // Calculate total size of the intermediate force buffer `fint`
        int nf_per_system = (nBonds_ * 2) + (nAngles_ * 3) + (nDihedrals_ * 4) + (nInversions_ * 4);

        // Allocate buffers on the GPU
        ibuff_apos        = newBuffer("apos",        nAtomsTot,       sizeof(cl_float4), 0, CL_MEM_READ_WRITE);
        ibuff_fapos       = newBuffer("fapos",       nAtomsTot,       sizeof(cl_float4), 0, CL_MEM_READ_WRITE);
        ibuff_REQs        = newBuffer("REQs",        nAtomsTot,       sizeof(cl_float4), 0, CL_MEM_READ_ONLY);
        ibuff_hneigh      = newBuffer("hneigh",      nAtomsTot * 4,   sizeof(cl_float4), 0, CL_MEM_READ_WRITE);
        ibuff_fint        = newBuffer("fint",        nSystems * nf_per_system, sizeof(cl_float4), 0, CL_MEM_READ_WRITE);

        ibuff_neighs      = newBuffer("neighs",      nAtomsTot,       sizeof(cl_int4),   0, CL_MEM_READ_ONLY);
        ibuff_neighCell   = newBuffer("neighCell",   nAtomsTot,       sizeof(cl_int4),   0, CL_MEM_READ_ONLY);
        ibuff_neighBs     = newBuffer("neighBs",     nAtomsTot,       sizeof(cl_int4),   0, CL_MEM_READ_ONLY);

        ibuff_bonAtoms    = newBuffer("bonAtoms",    nSystems * nBonds_,      sizeof(cl_int2),   0, CL_MEM_READ_ONLY);
        ibuff_bonParams   = newBuffer("bonParams",   nSystems * nBonds_,      sizeof(cl_float2), 0, CL_MEM_READ_ONLY);

        ibuff_angAtoms    = newBuffer("angAtoms",    nSystems * nAngles_,     sizeof(cl_int4),   0, CL_MEM_READ_ONLY); // Padded to int4
        ibuff_angNgs      = newBuffer("angNgs",      nSystems * nAngles_,     sizeof(cl_int2),   0, CL_MEM_READ_ONLY);
        ibuff_angParams1  = newBuffer("angParams1",  nSystems * nAngles_,     sizeof(cl_float4), 0, CL_MEM_READ_ONLY);
        ibuff_angParams2_w= newBuffer("angParams2_w",nSystems * nAngles_,     sizeof(cl_float),  0, CL_MEM_READ_ONLY);

        ibuff_dihAtoms    = newBuffer("dihAtoms",    nSystems * nDihedrals_,  sizeof(cl_int4),   0, CL_MEM_READ_ONLY);
        ibuff_dihNgs      = newBuffer("dihNgs",      nSystems * nDihedrals_,  sizeof(cl_int4),   0, CL_MEM_READ_ONLY); // Padded to int4
        ibuff_dihParams   = newBuffer("dihParams",   nSystems * nDihedrals_,  sizeof(cl_float4), 0, CL_MEM_READ_ONLY); // Padded to float4

        ibuff_invAtoms    = newBuffer("invAtoms",    nSystems * nInversions_, sizeof(cl_int4),   0, CL_MEM_READ_ONLY);
        ibuff_invNgs      = newBuffer("invNgs",      nSystems * nInversions_, sizeof(cl_int4),   0, CL_MEM_READ_ONLY); // Padded to int4
        ibuff_invParams   = newBuffer("invParams",   nSystems * nInversions_, sizeof(cl_float4), 0, CL_MEM_READ_ONLY);

        ibuff_a2f_offsets = newBuffer("a2f_offsets", nAtomsTot,       sizeof(cl_int),    0, CL_MEM_READ_ONLY);
        ibuff_a2f_counts  = newBuffer("a2f_counts",  nAtomsTot,       sizeof(cl_int),    0, CL_MEM_READ_ONLY);
        ibuff_a2f_indices = newBuffer("a2f_indices", nSystems * nA2F_,sizeof(cl_int),    0, CL_MEM_READ_ONLY);

        ibuff_pbcshifts   = newBuffer("pbcshifts",   nSystems * nPBC, sizeof(cl_float4), 0, CL_MEM_READ_ONLY);
        ibuff_lvecs       = newBuffer("lvecs",       nSystems,        sizeof(cl_Mat3),   0, CL_MEM_READ_ONLY);
        ibuff_energies    = newBuffer("energies",    nSystems * 5,    sizeof(cl_float),  0, CL_MEM_WRITE_ONLY); // E_b, E_a, E_d, E_i, E_tot
    }

    void setup_kernels() {
        // This function would set up the arguments for each kernel.
        // It's called once after realloc and before the main loop.
        // Example for one kernel:
        // task_evalBonds->args = {
        //     BUFFarg(ibuff_apos), BUFFarg(ibuff_fapos), ...
        // };
        // TODO: Implement argument binding for all kernels based on UFF.cl
    }

    void eval(bool bClearForce = true) {
        // This function enqueues all the kernels for a full UFF evaluation.
        if (bClearForce) {
            // Enqueue a kernel to zero the force buffers (fapos, fint)
        }

        if (bUFF_bonds)      task_evalBonds->enque();
        if (bUFF_angles)     task_evalAngles->enque();
        if (bUFF_dihedrals)  task_evalDihedrals->enque();
        if (bUFF_inversions) task_evalInversions->enque();

        task_assemble->enque();
    }

    void download_results(float* fapos_host, float* energies_host = nullptr) {
        download(ibuff_fapos, fapos_host);
        if (energies_host) {
            download(ibuff_energies, energies_host);
        }
    }
};

#endif