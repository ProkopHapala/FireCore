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
    bool bUFF_assemble   = true;
    bool bSubtractNB     = true; // Subtract 1-2 (bond) non-bonded interactions
    bool bSubtractNB_angle = true; // Subtract 1-3 (angle) non-bonded interactions
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

    // Optional per-interaction energy contributions
    int ibuff_Ea = -1; // per-angle energies [nSystems*nAngles]
    int ibuff_Ed = -1; // per-dihedral energies [nSystems*nDihedrals]
    int ibuff_Ei = -1; // per-inversion energies [nSystems*nInversions]

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
    OCLtask* task_clear_fapos   = nullptr;
    OCLtask* task_clear_fint    = nullptr;
    bool bKernelPrepared = false;

    // ====================== Functions

    void makeKernels(const char* cl_src_dir) {
        char srcpath[1024];
        sprintf(srcpath, "%s/UFF.cl", cl_src_dir);
        buildProgram(srcpath, program); // Assuming 'program' is the member from OCL base class
        // Create tasks for each kernel
        // TODO: The local work-group sizes (e.g., 32) are hardcoded for now. They should be tuned for optimal performance based on the device and kernel characteristics.
        //                                name                  program  nL nG
        newTask("clear_fapos_UFF",        program, 1, (size_t4){0,0,0,0}, (size_t4){32,1,1,1} );
        newTask("clear_fint_UFF",         program, 1, (size_t4){0,0,0,0}, (size_t4){32,1,1,1} );
        newTask("evalBondsAndHNeigh_UFF", program, 1, (size_t4){0,0,0,0}, (size_t4){32,1,1,1} );
        newTask("evalAngles_UFF",         program, 1, (size_t4){0,0,0,0}, (size_t4){32,1,1,1} );
        newTask("evalDihedrals_UFF",      program, 1, (size_t4){0,0,0,0}, (size_t4){32,1,1,1} );
        newTask("evalInversions_UFF",     program, 1, (size_t4){0,0,0,0}, (size_t4){32,1,1,1} );
        newTask("assembleForces_UFF",     program, 1, (size_t4){0,0,0,0}, (size_t4){32,1,1,1} );
        bKernelPrepared = false;
    }

    void realloc(int nSystems_, int nAtoms_, int nBonds_, int nAngles_, int nDihedrals_, int nInversions_, int nPBC_, int nA2F_) {
        printf("DEBUG: OCL_UFF::realloc() START\n");
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

        // Debug summary of counts
        printf("OCL_UFF::realloc counts: nSystems=%d nAtoms=%d nBonds=%d nAngles=%d nDihedrals=%d nInversions=%d nPBC=%d nA2F=%d nAtomsTot=%d nf_per_system=%d\n",
               nSystems, nAtoms, nBonds, nAngles, nDihedrals, nInversions, nPBC, nA2F, nAtomsTot, nf_per_system);

        // Safe sizes to avoid zero-sized buffer creation (CL_INVALID_BUFFER_SIZE)
        auto safeN = [](int n){ return (n>0)? n: 1; };
        int nBondsTot      = nSystems * nBonds_;
        int nAnglesTot     = nSystems * nAngles_;
        int nDihedralsTot  = nSystems * nDihedrals_;
        int nInversionsTot = nSystems * nInversions_;
        int nA2FTot        = nSystems * nA2F_;

        // Allocate buffers on the GPU (with debug prints and safe sizes where appropriate)
        ibuff_apos        = newBuffer("apos",        nAtomsTot,                   sizeof(cl_float4), 0, CL_MEM_READ_WRITE);
        ibuff_fapos       = newBuffer("fapos",       nAtomsTot,                   sizeof(cl_float4), 0, CL_MEM_READ_WRITE);
        ibuff_REQs        = newBuffer("REQs",        nAtomsTot,                   sizeof(cl_float4), 0, CL_MEM_READ_ONLY);
        ibuff_hneigh      = newBuffer("hneigh",      nAtomsTot * 4,               sizeof(cl_float4), 0, CL_MEM_READ_WRITE);
        ibuff_fint        = newBuffer("fint",        nSystems * nf_per_system,    sizeof(cl_float4), 0, CL_MEM_READ_WRITE);
        ibuff_neighs      = newBuffer("neighs",      nAtomsTot,                   sizeof(cl_int4),   0, CL_MEM_READ_ONLY);
        ibuff_neighCell   = newBuffer("neighCell",   nAtomsTot,                   sizeof(cl_int4),   0, CL_MEM_READ_ONLY);
        ibuff_neighBs     = newBuffer("neighBs",     nAtomsTot,                   sizeof(cl_int4),   0, CL_MEM_READ_ONLY);
        ibuff_bonAtoms    = newBuffer("bonAtoms",    safeN(nBondsTot),            sizeof(cl_int2),   0, CL_MEM_READ_ONLY);
        ibuff_bonParams   = newBuffer("bonParams",   safeN(nBondsTot),            sizeof(cl_float2), 0, CL_MEM_READ_ONLY);
        ibuff_angAtoms    = newBuffer("angAtoms",    safeN(nAnglesTot),           sizeof(cl_int4),   0, CL_MEM_READ_ONLY); // Padded to int4
        ibuff_angNgs      = newBuffer("angNgs",      safeN(nAnglesTot),           sizeof(cl_int2),   0, CL_MEM_READ_ONLY);
        ibuff_angParams1  = newBuffer("angParams1",  safeN(nAnglesTot),           sizeof(cl_float4), 0, CL_MEM_READ_ONLY);
        ibuff_angParams2_w= newBuffer("angParams2_w",safeN(nAnglesTot),           sizeof(cl_float),  0, CL_MEM_READ_ONLY);
        ibuff_dihAtoms    = newBuffer("dihAtoms",    safeN(nDihedralsTot),        sizeof(cl_int4),   0, CL_MEM_READ_ONLY);
        ibuff_dihNgs      = newBuffer("dihNgs",      safeN(nDihedralsTot),        sizeof(cl_int4),   0, CL_MEM_READ_ONLY); // Padded to int4
        ibuff_dihParams   = newBuffer("dihParams",   safeN(nDihedralsTot),        sizeof(cl_float4), 0, CL_MEM_READ_ONLY); // Padded to float4
        ibuff_invAtoms    = newBuffer("invAtoms",    safeN(nInversionsTot),       sizeof(cl_int4),   0, CL_MEM_READ_ONLY);
        ibuff_invNgs      = newBuffer("invNgs",      safeN(nInversionsTot),       sizeof(cl_int4),   0, CL_MEM_READ_ONLY); // Padded to int4
        ibuff_invParams   = newBuffer("invParams",   safeN(nInversionsTot),       sizeof(cl_float4), 0, CL_MEM_READ_ONLY);

        ibuff_a2f_offsets = newBuffer("a2f_offsets", nAtomsTot,                   sizeof(cl_int),    0, CL_MEM_READ_ONLY);
        ibuff_a2f_counts  = newBuffer("a2f_counts",  nAtomsTot,                   sizeof(cl_int),    0, CL_MEM_READ_ONLY);
        ibuff_a2f_indices = newBuffer("a2f_indices", safeN(nA2FTot),              sizeof(cl_int),    0, CL_MEM_READ_ONLY);

        // Use minimum size of 1 for buffers that could be 0-sized
        int nPBC_safe = (nPBC > 0) ? nPBC : 1;
        ibuff_pbcshifts   = newBuffer("pbcshifts",   nSystems * nPBC_safe,        sizeof(cl_float4), 0, CL_MEM_READ_ONLY);
        ibuff_lvecs       = newBuffer("lvecs",       nSystems,                    sizeof(cl_Mat3),   0, CL_MEM_READ_ONLY);
        ibuff_energies    = newBuffer("energies",    nSystems * 5,                sizeof(cl_float),  0, CL_MEM_WRITE_ONLY); // E_b, E_a, E_d, E_i, E_tot

        // Optional energy contributions per interaction
        if(nAngles_>0)    ibuff_Ea = newBuffer("Ea_contrib", nSystems * nAngles_,    sizeof(cl_float), 0, CL_MEM_WRITE_ONLY);
        if(nDihedrals_>0) ibuff_Ed = newBuffer("Ed_contrib", nSystems * nDihedrals_, sizeof(cl_float), 0, CL_MEM_WRITE_ONLY);
        if(nInversions_>0)ibuff_Ei = newBuffer("Ei_contrib", nSystems * nInversions_,sizeof(cl_float), 0, CL_MEM_WRITE_ONLY);
    }

    void setup_kernels( float Rdamp, float FmaxNonBonded, float SubNBTorsionFactor ) {
        printf("OCL_UFF::setup_kernels()\n");
        // Common temporary scalars; TODO wire from UFF host
        //float Rdamp              = 0.0f;
        //float FmaxNonBonded      = 1e6f;
        //float SubNBTorsionFactor = 0.0f;

        // Compute offsets into fint buffer to MATCH CPU UFF.h layout exactly:
        // CPU: nf = 4*ndihedrals + 4*ninversions + 3*nangles + nbonds;
        //       i0dih=0; i0inv=i0dih+4*ndihedrals; i0ang=i0inv+4*ninversions; i0bon=i0ang+3*nangles;
        int i0dih = 0;
        int i0inv = i0dih + nDihedrals * 4;
        int i0ang = i0inv + nInversions * 4;
        int i0bon = i0ang + nAngles * 3;
        int nf_per_system = (nDihedrals * 4) + (nInversions * 4) + (nAngles * 3) + (nBonds);

        // Disable components with zero interactions
        bUFF_bonds      &= (nBonds     > 0);
        bUFF_angles     &= (nAngles    > 0);
        bUFF_dihedrals  &= (nDihedrals > 0);
        bUFF_inversions &= (nInversions> 0);

        bKernelPrepared = false;
        // Get task pointers
        task_evalBonds      = getTask("evalBondsAndHNeigh_UFF");
        task_evalAngles     = getTask("evalAngles_UFF");
        task_evalDihedrals  = getTask("evalDihedrals_UFF");
        task_evalInversions = getTask("evalInversions_UFF");
        task_assemble       = getTask("assembleForces_UFF");
        task_clear_fapos    = getTask("clear_fapos_UFF");
        task_clear_fint     = getTask("clear_fint_UFF");

        // --- evalBondsAndHNeigh_UFF ---
        if(task_evalBonds){
            printf("OCL_UFF::setup_kernels().task_evalBonds \n");
            int nloc = 32;
            task_evalBonds->local.x  = nloc;
            task_evalBonds->global.x = nAtoms + nloc - (nAtoms % nloc);
            task_evalBonds->global.y = nSystems;
            useKernel( task_evalBonds->ikernel );
            int err=0;
            err |= useArg   ( nAtoms );                    OCL_checkError(err, "evalBonds.arg1 nAtoms");
            err |= useArg   ( nPBC );                      OCL_checkError(err, "evalBonds.arg2 nPBC");
            err |= useArg   ( i0bon );                     OCL_checkError(err, "evalBonds.arg3 i0bon");
            int bSub = bSubtractNB?1:0;
            err |= useArg   ( bSub );                      OCL_checkError(err, "evalBonds.arg4 bSub");
            err |= useArg   ( Rdamp );                     OCL_checkError(err, "evalBonds.arg5 Rdamp");
            err |= useArg   ( FmaxNonBonded );             OCL_checkError(err, "evalBonds.arg6 FmaxNB");
            err |= useArgBuff( ibuff_apos );               OCL_checkError(err, "evalBonds.arg7 apos");
            err |= useArgBuff( ibuff_fapos );              OCL_checkError(err, "evalBonds.arg8 fapos");
            err |= useArgBuff( ibuff_neighs );             OCL_checkError(err, "evalBonds.arg9 neighs");
            err |= useArgBuff( ibuff_neighCell );          OCL_checkError(err, "evalBonds.arg10 neighCell");
            err |= useArgBuff( ibuff_pbcshifts );          OCL_checkError(err, "evalBonds.arg11 pbc_shifts");
            err |= useArgBuff( ibuff_neighBs );            OCL_checkError(err, "evalBonds.arg12 neighBs");
            err |= useArgBuff( ibuff_bonParams );          OCL_checkError(err, "evalBonds.arg13 bonParams");
            err |= useArgBuff( ibuff_REQs );               OCL_checkError(err, "evalBonds.arg14 REQs");
            err |= useArgBuff( ibuff_bonAtoms );           OCL_checkError(err, "evalBonds.arg15 bonAtoms");
            err |= useArgBuff( ibuff_hneigh );             OCL_checkError(err, "evalBonds.arg16 hneigh");
            err |= useArgBuff( ibuff_fint );               OCL_checkError(err, "evalBonds.arg17 fint");
            //int bAssemble = bUFF_assemble ? 1 : 0;
            //err |= useArg   ( bAssemble );                 OCL_checkError(err, "evalBonds.arg18 bAssemble");
        }

        // --- evalAngles_UFF ---
        if(task_evalAngles){
            printf("OCL_UFF::setup_kernels().task_evalAngles \n");
            int nloc = 32;
            task_evalAngles->local.x  = nloc;
            task_evalAngles->global.x = (nAngles>0)? (nAngles + nloc - (nAngles % nloc)) : 0;
            task_evalAngles->global.y = nSystems;
            useKernel( task_evalAngles->ikernel );
            int err=0;
            // Offsets match CPU mapping (angles start after dihedrals and inversions)
            int bSub = bSubtractNB_angle?1:0;
            err |= useArg   ( nAngles );                   OCL_checkError(err, "evalAngles.arg1 nAngles");
            err |= useArg   ( i0ang );                     OCL_checkError(err, "evalAngles.arg2 i0ang");
            err |= useArg   ( bSub );                      OCL_checkError(err, "evalAngles.arg3 bSub");
            err |= useArg   ( Rdamp );                     OCL_checkError(err, "evalAngles.arg4 Rdamp");
            err |= useArg   ( FmaxNonBonded );             OCL_checkError(err, "evalAngles.arg5 FmaxNB");
            err |= useArgBuff( ibuff_angAtoms );           OCL_checkError(err, "evalAngles.arg6 angAtoms");
            err |= useArgBuff( ibuff_angNgs );             OCL_checkError(err, "evalAngles.arg7 angNgs");
            err |= useArgBuff( ibuff_angParams1 );         OCL_checkError(err, "evalAngles.arg8 angParams1");
            err |= useArgBuff( ibuff_angParams2_w );       OCL_checkError(err, "evalAngles.arg9 angParams2_w");
            err |= useArgBuff( ibuff_hneigh );             OCL_checkError(err, "evalAngles.arg10 hneigh");
            err |= useArgBuff( ibuff_REQs );               OCL_checkError(err, "evalAngles.arg11 REQs");
            err |= useArgBuff( ibuff_apos );               OCL_checkError(err, "evalAngles.arg12 apos");
            err |= useArgBuff( ibuff_pbcshifts );          OCL_checkError(err, "evalAngles.arg13 pbc_shifts");
            err |= useArgBuff( ibuff_neighs );             OCL_checkError(err, "evalAngles.arg14 neighs");
            err |= useArgBuff( ibuff_neighCell );          OCL_checkError(err, "evalAngles.arg15 neighCell");
            err |= useArg   ( nPBC );                      OCL_checkError(err, "evalAngles.arg16 nPBC");
            err |= useArgBuff( ibuff_fint );               OCL_checkError(err, "evalAngles.arg17 fint");
            if(ibuff_Ea>=0){ err |= useArgBuff( ibuff_Ea ); OCL_checkError(err, "evalAngles.arg18 Ea_contrib"); }
        }

        // --- evalDihedrals_UFF ---
        if(task_evalDihedrals){
            printf("OCL_UFF::setup_kernels().task_evalDihedrals \n");
            int nloc = 32;
            task_evalDihedrals->local.x  = nloc;
            task_evalDihedrals->global.x = (nDihedrals>0)? (nDihedrals + nloc - (nDihedrals % nloc)) : 0;
            task_evalDihedrals->global.y = nSystems;
            useKernel( task_evalDihedrals->ikernel );
            int err=0;
            // offset after angles
            err |= useArg   ( nDihedrals );                OCL_checkError(err, "evalDihedrals.arg1 nDihedrals");
            err |= useArg   ( i0dih );                     OCL_checkError(err, "evalDihedrals.arg2 i0dih");
            err |= useArg   ( SubNBTorsionFactor );        OCL_checkError(err, "evalDihedrals.arg3 SubNB");
            err |= useArg   ( Rdamp );                     OCL_checkError(err, "evalDihedrals.arg4 Rdamp");
            err |= useArg   ( FmaxNonBonded );             OCL_checkError(err, "evalDihedrals.arg5 FmaxNB");
            err |= useArgBuff( ibuff_dihAtoms );           OCL_checkError(err, "evalDihedrals.arg6 dihAtoms");
            err |= useArgBuff( ibuff_dihNgs );             OCL_checkError(err, "evalDihedrals.arg7 dihNgs");
            err |= useArgBuff( ibuff_dihParams );          OCL_checkError(err, "evalDihedrals.arg8 dihParams");
            err |= useArgBuff( ibuff_hneigh );             OCL_checkError(err, "evalDihedrals.arg9 hneigh");
            err |= useArgBuff( ibuff_REQs );               OCL_checkError(err, "evalDihedrals.arg10 REQs");
            err |= useArgBuff( ibuff_apos );               OCL_checkError(err, "evalDihedrals.arg11 apos");
            err |= useArgBuff( ibuff_pbcshifts );          OCL_checkError(err, "evalDihedrals.arg12 pbc_shifts");
            err |= useArgBuff( ibuff_neighs );             OCL_checkError(err, "evalDihedrals.arg13 neighs");
            err |= useArgBuff( ibuff_neighCell );          OCL_checkError(err, "evalDihedrals.arg14 neighCell");
            err |= useArg   ( nPBC );                      OCL_checkError(err, "evalDihedrals.arg15 nPBC");
            err |= useArgBuff( ibuff_fint );               OCL_checkError(err, "evalDihedrals.arg16 fint");
            if(ibuff_Ed>=0){ err |= useArgBuff( ibuff_Ed ); OCL_checkError(err, "evalDihedrals.arg17 Ed_contrib"); }
        }

        // --- evalInversions_UFF ---
        if(task_evalInversions){
            printf("OCL_UFF::setup_kernels().task_evalInversions \n");
            int nloc = 32;
            task_evalInversions->local.x  = nloc;
            task_evalInversions->global.x = (nInversions>0)? (nInversions + nloc - (nInversions % nloc)) : 0;
            task_evalInversions->global.y = nSystems;
            useKernel( task_evalInversions->ikernel );
            int err=0;
            // offset after dihedrals
            err |= useArg   ( nInversions );               OCL_checkError(err, "evalInversions.arg1 nInversions");
            err |= useArg   ( i0inv );                     OCL_checkError(err, "evalInversions.arg2 i0inv");
            err |= useArgBuff( ibuff_invAtoms );           OCL_checkError(err, "evalInversions.arg3 invAtoms");
            err |= useArgBuff( ibuff_invNgs );             OCL_checkError(err, "evalInversions.arg4 invNgs");
            err |= useArgBuff( ibuff_invParams );          OCL_checkError(err, "evalInversions.arg5 invParams");
            err |= useArgBuff( ibuff_hneigh );             OCL_checkError(err, "evalInversions.arg6 hneigh");
            err |= useArgBuff( ibuff_fint );               OCL_checkError(err, "evalInversions.arg7 fint");
            if(ibuff_Ei>=0){ err |= useArgBuff( ibuff_Ei ); OCL_checkError(err, "evalInversions.arg8 Ei_contrib"); }
        }

        // --- assembleForces_UFF ---
        if(task_assemble){
            printf("OCL_UFF::setup_kernels().task_assemble \n");
            int nloc = 32;
            task_assemble->local.x  = nloc;
            task_assemble->global.x = nAtoms + nloc - (nAtoms % nloc);
            task_assemble->global.y = nSystems;
            useKernel( task_assemble->ikernel );
            int err=0;
            err |= useArg   ( nAtoms );                    OCL_checkError(err, "assemble.arg1 nAtoms");
            err |= useArgBuff( ibuff_fint );               OCL_checkError(err, "assemble.arg2 fint");
            err |= useArgBuff( ibuff_a2f_offsets );        OCL_checkError(err, "assemble.arg3 a2f_offsets");
            err |= useArgBuff( ibuff_a2f_counts );         OCL_checkError(err, "assemble.arg4 a2f_counts");
            err |= useArgBuff( ibuff_a2f_indices );        OCL_checkError(err, "assemble.arg5 a2f_indices");
            err |= useArgBuff( ibuff_fapos );              OCL_checkError(err, "assemble.arg6 fapos");
            int bClearForce = 1;
            err |= useArg   ( bClearForce );               OCL_checkError(err, "assemble.arg7 bClearForce");
        }

        // --- clear_fapos_UFF ---
        if(task_clear_fapos){
            int nloc=32;
            task_clear_fapos->local.x  = nloc;
            task_clear_fapos->global.x = nAtomsTot + nloc - (nAtomsTot % nloc);
            useKernel( task_clear_fapos->ikernel );
            int err=0;
            err |= useArg   ( nAtomsTot ); OCL_checkError(err, "clear_fapos.arg1 nAtomsTot");
            err |= useArgBuff( ibuff_fapos ); OCL_checkError(err, "clear_fapos.arg2 fapos");
        }

        // --- clear_fint_UFF ---
        if(task_clear_fint){
            int nloc=32;
            task_clear_fint->local.x  = nloc;
            task_clear_fint->global.x = nSystems*nf_per_system + nloc - ((nSystems*nf_per_system) % nloc);
            useKernel( task_clear_fint->ikernel );
            int err=0;
            int nTot = nSystems*nf_per_system;
            err |= useArg   ( nTot ); OCL_checkError(err, "clear_fint.arg1 nTot");
            err |= useArgBuff( ibuff_fint ); OCL_checkError(err, "clear_fint.arg2 fint");
        }

        bKernelPrepared = true;
    }

    void eval(bool bClearForce = true) {
        printf("OCL_UFF::eval() bClearForce=%i bUFF_bonds=%i bUFF_angles=%i bUFF_dihedrals=%i bUFF_inversions=%i bUFF_assemble=%i\n", bClearForce, bUFF_bonds, bUFF_angles, bUFF_dihedrals, bUFF_inversions, bUFF_assemble);
        // bUFF_bonds      &= nBonds>0;
        // bUFF_angles     &= nAngles>0;
        // bUFF_dihedrals  &= nDihedrals>0;
        // bUFF_inversions &= nInversions>0;
        // bUFF_assemble   &= nAtoms>0;
        // printf("OCL_UFF::eval() AFTER bClearForce=%i bUFF_bonds=%i bUFF_angles=%i bUFF_dihedrals=%i bUFF_inversions=%i bUFF_assemble=%i\n", bClearForce, bUFF_bonds, bUFF_angles, bUFF_dihedrals, bUFF_inversions, bUFF_assemble);
        // This function enqueues all the kernels for a full UFF evaluation.
        if (bClearForce) {
            if(task_clear_fapos){ task_clear_fapos->enque(); }
            if(task_clear_fint ){ task_clear_fint ->enque(); }
        }
        if (bUFF_bonds      ){ printf("OCL_UFF::eval().task_evalBonds      \n"); task_evalBonds     ->enque(); }
        if (bUFF_angles     ){ printf("OCL_UFF::eval().task_evalAngles     \n"); task_evalAngles    ->enque(); }
        if (bUFF_dihedrals  ){ printf("OCL_UFF::eval().task_evalDihedrals  \n"); task_evalDihedrals ->enque(); }
        if (bUFF_inversions ){ printf("OCL_UFF::eval().task_evalInversions \n"); task_evalInversions->enque(); }
        if (bUFF_assemble   ){ printf("OCL_UFF::eval().task_assemble       \n"); task_assemble      ->enque(); }
        printf("OCL_UFF::eval() DONE");
    }

    void download_results(float* fapos_host, float* energies_host = nullptr) {
        download(ibuff_fapos, fapos_host);
        if (energies_host) {
            download(ibuff_energies, energies_host);
        }
    }
};

#endif