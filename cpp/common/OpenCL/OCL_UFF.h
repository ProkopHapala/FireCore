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
    int  npbc       = 0;
    int4 nPBC       = {0,0,0,0};
    int nAtomsTot   = 0;
    int nA2F        = 0;

    // --- Component Flags
    bool bUFF_bonds      = true;
    bool bUFF_angles     = true;
    bool bUFF_dihedrals  = true;
    bool bUFF_inversions = true;
    bool bUFF_assemble   = true;
    bool bNonBond        = true;
    bool bGridFF         = true;
    bool bSubtractNB     = true; // Subtract 1-2 (bond) non-bonded interactions
    bool bSubtractNB_angle = true; // Subtract 1-3 (angle) non-bonded interactions
    bool bClampNonBonded = true;
    bool bUseTexture     = false;

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
    int ibuff_excl      = -1; // Second-neighbor exclusion table
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
    int ibuff_ilvecs    = -1; // Inverse lattice vectors for each system (optional, kept for parity)
    int ibuff_energies  = -1; // Buffer to store computed energies

    // --- GridFF Buffers and Parameters (copied from OCL_MM)
    int itex_BsplinePLQH = -1;
    int ibuff_BsplinePLQ = -1;
    int4 grid_n;
    Quat4f grid_p0;
    Quat4f grid_shift0;
    Quat4f grid_shift0_p0;
    Quat4f grid_step     { 0.f, 0.f, 0.f, 0.f }; // grid cell step ( assuming it is orghogonal )
    Quat4f grid_invStep  { 0.f, 0.f, 0.f, 0.f }; // inverse grid cell step (assuming it is orghogonal)
    float4 GFFparams{1.0, 1.5, 0., 0.}; // (Rdamp, alphaMorse, 0, 0)
    int ibuff_atoms_surf = -1;
    int ibuff_REQs_surf  = -1;
    int natom_surf=0;
    cl_Mat3 cl_dGrid;       // grid cell step (voxel rhombus)
    cl_Mat3 cl_diGrid;      // inverse grid cell step
    cl_Mat3 cl_grid_lvec;   // grid lattice vectors
    cl_Mat3 cl_grid_ilvec;  // inverse grid lattice vectors

    // --- Common MD/update state (needed for atom updates on GPU)
    int4   nDOFs{0,0,0,0};   // (natoms, nnode, unused, nMaxSysNeighs)

    // --- Buffers for atom updates (shared with MM path naming)
    int ibuff_avel      = -1; // velocities (x,y,z,m)
    int ibuff_cvf       = -1; // accumulators for FIRE damping { |f|^2, |v|^2, <f|v>, 0 }
    int ibuff_constr    = -1; // constraints target positions (x,y,z,K)
    int ibuff_constrK   = -1; // constraint stiffness (kx,ky,kz,unused)
    int ibuff_MDpars    = -1; // per-system MD parameters (dt, damp, Flimit, reserved)
    int ibuff_TDrive    = -1; // per-system thermal driving (T, gamma, seed, reserved)
    int ibuff_bboxes    = -1; // per-system bounding boxes (cl_Mat3)
    int ibuff_sysneighs = -1; // optional inter-system neighbor list
    int ibuff_sysbonds  = -1; // optional inter-system bond parameters

    // --- OpenCL Kernels (identified by OCLtask object)
    OCLtask* task_evalBonds     = nullptr;
    OCLtask* task_evalAngles    = nullptr;
    OCLtask* task_evalDihedrals = nullptr;
    OCLtask* task_evalInversions= nullptr;
    OCLtask* task_assemble      = nullptr;
    OCLtask* task_updateAtoms   = nullptr;
    OCLtask* task_clear_fapos   = nullptr;
    OCLtask* task_clear_fint    = nullptr;
    OCLtask* task_NBFF = nullptr;
    OCLtask* task_NBFF_ex2 = nullptr;
    OCLtask* task_NBFF_Grid_Bspline = nullptr;
    OCLtask* task_NBFF_Grid_Bspline_ex2 = nullptr;
    OCLtask* task_SurfAtoms   = nullptr;
    bool bKernelPrepared = false;
    bool bSetUp = false;

    // ====================== Functions

    void setGridShape( const GridShape& grid ){
        printf("DEBUG OCL_UFF::setGridShape() called\n");
        v2i4      ( grid.n      , grid_n       );
        grid_p0.f = (Vec3f)grid.pos0;
        // Set grid step and inverse step from diagonal components of dCell matrix
        grid_step.f   = Vec3f{ (float)grid.dCell.xx, (float)grid.dCell.yy, (float)grid.dCell.zz };
        grid_invStep.f.set_inv( grid_step.f );
        Mat3_to_cl( grid.dCell  , cl_dGrid     );
        Mat3_to_cl( grid.diCell , cl_diGrid    );
        Mat3_to_cl( grid.cell   , cl_grid_lvec );
        Mat3_to_cl( grid.iCell  , cl_grid_ilvec );
        
        // Debug print to verify parameters are set correctly
        printf("DEBUG OCL_UFF::setGridShape() grid_p0=(%f,%f,%f) grid_step=(%f,%f,%f) grid_invStep=(%f,%f,%f)\n",
               grid_p0.f.x, grid_p0.f.y, grid_p0.f.z,
               grid_step.f.x, grid_step.f.y, grid_step.f.z,
               grid_invStep.f.x, grid_invStep.f.y, grid_invStep.f.z);
        printf("DEBUG OCL_UFF::setGridShape() grid_n=(%d,%d,%d,%d)\n", grid_n.x, grid_n.y, grid_n.z, grid_n.w);
    }

    void surf2ocl( const GridShape& grid, Vec3i nPBC_, int na=0, float4* atoms=0, float4* REQs=0 ){
        int err=0;
        setGridShape( grid );
        v2i4( nPBC_, nPBC );
        if(ibuff_atoms_surf<0) ibuff_atoms_surf = newBuffer( "atoms_surf", na, sizeof(float4), 0, CL_MEM_READ_ONLY );
        if(ibuff_REQs_surf <0) ibuff_REQs_surf  = newBuffer( "REQs_surf",  na, sizeof(float4), 0, CL_MEM_READ_ONLY );
        if(atoms){ err |= upload( ibuff_atoms_surf, atoms, na ); OCL_checkError(err, "makeGridFF().upload(atoms)" ); natom_surf = na; }
        if(REQs ){ err |= upload( ibuff_REQs_surf , REQs , na ); OCL_checkError(err, "makeGridFF().upload(REQs )" ); }
        OCL_checkError(err, "surf2ocl()" );
    }

    OCLtask* getSurfMorse(  Vec3i nPBC_, int na=0, float4* atoms=0, float4* REQs=0, int na_s=0, float4* atoms_s=0, float4* REQs_s=0,  bool bRun=true, OCLtask* task=0 ){
        v2i4( nPBC_, nPBC );
        int err=0;
        err |= finishRaw();       OCL_checkError(err, "getSurfMorse().imgAlloc" );
        nDOFs.x = nAtoms;
        nDOFs.y = 0; // nnode is 0 for UFF
        nDOFs.z = natom_surf;
        if(task==0) task = getTask("getSurfMorse");
        int nloc  = 32;
        task->local.x = nloc;
        task->global.x = nAtoms + nloc-(nAtoms%nloc);
        task->local.y = 1;
        task->global.y = nSystems;
        useKernel( task->ikernel );

        err |= _useArg   ( nDOFs );
        err |= useArgBuff( ibuff_apos      );
        err |= useArgBuff( ibuff_REQs       );
        err |= useArgBuff( ibuff_fapos    );
        err |= useArgBuff( ibuff_atoms_surf );
        err |= useArgBuff( ibuff_REQs_surf  );
        err |= _useArg( nPBC        );
        err |= _useArg( cl_grid_lvec    );
        err |= _useArg( grid_shift0     );
        err |= _useArg( GFFparams       );
        OCL_checkError(err, "getSurfMorse().setup");
        if(bRun){
            err |= task->enque_raw(); OCL_checkError(err, "getSurfMorse().enque"  );
            err |= finishRaw();       OCL_checkError(err, "getSurfMorse().finish" );
        }
        return task;
    }

    void makeKernels(const char* cl_src_dir) {
        char srcpath[1024];
        sprintf(srcpath, "%s/UFF.cl", cl_src_dir);
        // Build with debug flags enabled for pairwise interaction comparison
        // Print interactions for atom 0 with all other atoms
        char build_options[256];
        sprintf(build_options, "-I. -cl-std=CL2.0 ");
        buildProgram(srcpath, build_options); // Assuming 'program' is the member from OCL base class
        // Create tasks for each kernel
        // TODO: The local work-group sizes (e.g., 32) are hardcoded for now. They should be tuned for optimal performance based on the device and kernel characteristics.
        //                                name                  program  nL nG
        newTask("clear_fapos_UFF",        program, 2, (size_t4){0,0,0,0}, (size_t4){32,1,1,1} );
        newTask("clear_fint_UFF",         program, 2, (size_t4){0,0,0,0}, (size_t4){32,1,1,1} );
        newTask("evalBondsAndHNeigh_UFF", program, 2, (size_t4){0,0,0,0}, (size_t4){32,1,1,1} );
        newTask("evalAngles_UFF",         program, 2, (size_t4){0,0,0,0}, (size_t4){32,1,1,1} );
        newTask("evalDihedrals_UFF",      program, 2, (size_t4){0,0,0,0}, (size_t4){32,1,1,1} );
        newTask("evalInversions_UFF",     program, 2, (size_t4){0,0,0,0}, (size_t4){32,1,1,1} );
        newTask("assembleForces_UFF",     program, 2, (size_t4){0,0,0,0}, (size_t4){32,1,1,1} );
        newTask("updateAtomsMMFFf4",      program, 2, (size_t4){0,0,0,0}, (size_t4){32,1,1,1} );
        newTask("getNonBond",             program, 2, (size_t4){0,0,0,0}, (size_t4){32,1,1,1});
        newTask("getNonBond_ex2",         program, 2, (size_t4){0,0,0,0}, (size_t4){32,1,1,1});
        newTask("getNonBond_GridFF_Bspline", program, 2, (size_t4){0,0,0,0}, (size_t4){32,1,1,1});
        newTask("getNonBond_GridFF_Bspline_ex2", program, 2, (size_t4){0,0,0,0}, (size_t4){32,1,1,1});
        newTask("getSurfMorse",           program, 2, (size_t4){0,0,0,0}, (size_t4){32,1,1,1});
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
        nPBC.x      = nPBC_;
        npbc = (nPBC.x*2+1)*(nPBC.y*2+1)*(nPBC.z*2+1);
        nA2F        = nA2F_;
        nAtomsTot   = nSystems * nAtoms;

        // Calculate total size of the intermediate force buffer `fint`
        int nf_per_system = (nBonds_ * 2) + (nAngles_ * 3) + (nDihedrals_ * 4) + (nInversions_ * 4);

        // Debug summary of counts
        printf("OCL_UFF::realloc counts: nSystems=%d nAtoms=%d nBonds=%d nAngles=%d nDihedrals=%d nInversions=%d nPBC=%d nA2F=%d nAtomsTot=%d nf_per_system=%d\n",
               nSystems, nAtoms, nBonds, nAngles, nDihedrals, nInversions, npbc, nA2F, nAtomsTot, nf_per_system);

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
        ibuff_excl        = newBuffer("excl",        nAtomsTot*EXCL_MAX,          sizeof(cl_int ),   0, CL_MEM_READ_ONLY);
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
        int nPBC_safe = (npbc > 0) ? npbc : 1;
        ibuff_pbcshifts   = newBuffer("pbcshifts",   nSystems * nPBC_safe,        sizeof(cl_float4), 0, CL_MEM_READ_ONLY);
        ibuff_lvecs       = newBuffer("lvecs",       nSystems,                    sizeof(cl_Mat3),   0, CL_MEM_READ_ONLY);
        ibuff_ilvecs      = newBuffer("ilvecs",      nSystems,                    sizeof(cl_Mat3),   0, CL_MEM_READ_ONLY);
        ibuff_energies    = newBuffer("energies",    nSystems * 5,                sizeof(cl_float),  0, CL_MEM_WRITE_ONLY); // E_b, E_a, E_d, E_i, E_tot

        // Optional energy contributions per interaction
        if(nAngles_>0)    ibuff_Ea = newBuffer("Ea_contrib", nSystems * nAngles_,    sizeof(cl_float), 0, CL_MEM_WRITE_ONLY);
        if(nDihedrals_>0) ibuff_Ed = newBuffer("Ed_contrib", nSystems * nDihedrals_, sizeof(cl_float), 0, CL_MEM_WRITE_ONLY);
        if(nInversions_>0)ibuff_Ei = newBuffer("Ei_contrib", nSystems * nInversions_,sizeof(cl_float), 0, CL_MEM_WRITE_ONLY);

        // --- Allocate MD/update related buffers (sizes available now)
        // Note: For UFF we do not use pi orbitals; vectors count equals atoms count per system
        ibuff_avel     = newBuffer("avel",      nAtomsTot, sizeof(cl_float4), 0, CL_MEM_READ_WRITE);
        ibuff_cvf      = newBuffer("cvf",       nAtomsTot, sizeof(cl_float4), 0, CL_MEM_READ_WRITE);
        ibuff_constr   = newBuffer("constr",    nAtomsTot, sizeof(cl_float4), 0, CL_MEM_READ_WRITE);
        ibuff_constrK  = newBuffer("constrK",   nAtomsTot, sizeof(cl_float4), 0, CL_MEM_READ_WRITE);
        ibuff_MDpars   = newBuffer("MDpars",    nSystems,  sizeof(cl_float4), 0, CL_MEM_READ_ONLY );
        ibuff_TDrive   = newBuffer("TDrive",    nSystems,  sizeof(cl_float4), 0, CL_MEM_READ_ONLY );
        ibuff_bboxes   = newBuffer("bboxes",    nSystems,  sizeof(cl_Mat3),   0, CL_MEM_READ_ONLY );
        // Inter-system coupling (optional); allocate minimal placeholders
        ibuff_sysneighs= newBuffer("sysneighs", nSystems,  sizeof(cl_int),    0, CL_MEM_READ_ONLY );
        ibuff_sysbonds = newBuffer("sysbonds",  nSystems,  sizeof(cl_float4), 0, CL_MEM_READ_ONLY );
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
        task_updateAtoms     = getTask("updateAtomsMMFFf4");
        task_NBFF           = getTask("getNonBond");
        task_NBFF_ex2       = getTask("getNonBond_ex2");
        task_NBFF_Grid_Bspline = getTask("getNonBond_GridFF_Bspline");
        task_NBFF_Grid_Bspline_ex2 = getTask("getNonBond_GridFF_Bspline_ex2");

        // --- evalBondsAndHNeigh_UFF ---
        if(task_evalBonds && nBonds > 0){
            printf("OCL_UFF::setup_kernels().task_evalBonds \n");
            int nloc = 32;
            task_evalBonds->local.x  = nloc;
            task_evalBonds->global.x = nAtoms + nloc - (nAtoms % nloc);
            task_evalBonds->global.y = nSystems;
            useKernel( task_evalBonds->ikernel );
            int err=0;
            err |= useArg   ( nAtoms );                    OCL_checkError(err, "evalBonds.arg1 nAtoms");
            err |= useArg   ( npbc );                      OCL_checkError(err, "evalBonds.arg2 nPBC");
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
        if(task_evalAngles && nAngles > 0){
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
            err |= useArg   ( npbc );                      OCL_checkError(err, "evalAngles.arg16 nPBC");
            err |= useArgBuff( ibuff_fint );               OCL_checkError(err, "evalAngles.arg17 fint");
            if(ibuff_Ea>=0){ err |= useArgBuff( ibuff_Ea ); OCL_checkError(err, "evalAngles.arg18 Ea_contrib"); }
            err |= useArg   ( nf_per_system );             OCL_checkError(err, "evalAngles.arg19 nf_per_system");
        }

        // --- evalDihedrals_UFF ---
        if(task_evalDihedrals && nDihedrals > 0){
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
            err |= useArg   ( npbc );                      OCL_checkError(err, "evalDihedrals.arg15 nPBC");
            err |= useArgBuff( ibuff_fint );               OCL_checkError(err, "evalDihedrals.arg16 fint");
            if(ibuff_Ed>=0){ err |= useArgBuff( ibuff_Ed ); OCL_checkError(err, "evalDihedrals.arg17 Ed_contrib"); }
            err |= useArg   ( nf_per_system );             OCL_checkError(err, "evalDihedrals.arg18 nf_per_system");
        }

        // --- evalInversions_UFF ---
        if(task_evalInversions && nInversions > 0){
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
            err |= useArg   ( nf_per_system );             OCL_checkError(err, "evalInversions.arg9 nf_per_system");
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
            err |= useArg   ( nf_per_system );             OCL_checkError(err, "assemble.arg8 nf_per_system");
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

        // if(task_NBFF && bNonBond && (nAtoms>0)){
        //     task_NBFF = setup_getNonBond(nAtoms, Vec3i{nPBC.x,nPBC.y,nPBC.z}, task_NBFF);
        // }
        // if(task_NBFF_ex2 && bNonBond && (nAtoms>0)){
        //     task_NBFF_ex2 = setup_getNonBond_ex2(nAtoms, Vec3i{nPBC.x,nPBC.y,nPBC.z}, task_NBFF_ex2);
        // }
        // if(task_NBFF_Grid_Bspline && bGridFF && (nAtoms>0)){
        //     task_NBFF_Grid_Bspline = setup_getNonBond_GridFF_Bspline(nAtoms, Vec3i{nPBC.x,nPBC.y,nPBC.z}, task_NBFF_Grid_Bspline);
        // }
        // if(task_NBFF_Grid_Bspline_ex2 && bGridFF && (nAtoms>0)){
        //     task_NBFF_Grid_Bspline_ex2 = setup_getNonBond_GridFF_Bspline_ex2(nAtoms, Vec3i{nPBC.x,nPBC.y,nPBC.z}, task_NBFF_Grid_Bspline_ex2);
        // }
        bKernelPrepared = true;
    }

    OCLtask* setup_getNonBond(int na, Vec3i nPBC_, OCLtask* task = nullptr) {
        if (!task) task = getTask("getNonBond");
        int nloc = 32;
        task->local.x = nloc;
        task->local.y = 1;  // CRITICAL FIX: Set local size for systems dimension to 1
        task->global.x = na + nloc - (na % nloc);
        task->global.y = nSystems;
        useKernel(task->ikernel);
        nDOFs = (int4){na,0,0,0};
        // nDOFs.y = 0; // nNode is 0 for UFF
        // nDOFs.z = 0; // nNode is 0 for UFF
        // nDOFs.w = 0; // nNode is 0 for UFF
        int4 npbc_int4;
        v2i4(nPBC_, npbc_int4);
        int err = 0;
        err |= _useArg(nDOFs);                     OCL_checkError(err, "setup_getNonBond.arg 1");      // 1
        err |= useArgBuff(ibuff_apos);             OCL_checkError(err, "setup_getNonBond.arg 2");      // 2
        err |= useArgBuff(ibuff_fapos);            OCL_checkError(err, "setup_getNonBond.arg 3");      // 3
        err |= useArgBuff(ibuff_REQs);             OCL_checkError(err, "setup_getNonBond.arg 4");      // 4
        err |= useArgBuff(ibuff_neighs);           OCL_checkError(err, "setup_getNonBond.arg 5");      // 5
        err |= useArgBuff(ibuff_neighCell);        OCL_checkError(err, "setup_getNonBond.arg 6");      // 6
        err |= useArgBuff(ibuff_lvecs);            OCL_checkError(err, "setup_getNonBond.arg 7");      // 7
        err |= _useArg(npbc_int4);                 OCL_checkError(err, "setup_getNonBond.arg 8");      // 9
        err |= _useArg(GFFparams);                 OCL_checkError(err, "setup_getNonBond.arg 9");      // 10
        OCL_checkError(err, "setup_getNonBond_UFF");
        return task;
    }

    OCLtask* setup_getNonBond_ex2(int na, Vec3i nPBC_, OCLtask* task = nullptr) {
        if (!task) task = getTask("getNonBond_ex2");
        int nloc = 32;
        task->local.x = nloc;
        task->local.y = 1;  // CRITICAL FIX: Set local size for systems dimension to 1
        task->global.x = na + nloc - (na % nloc);
        task->global.y = nSystems;
        useKernel(task->ikernel);
        nDOFs = (int4){na,0,0,0};
        // nDOFs.y = 0; // nNode is 0 for UFF
        // nDOFs.z = 0; // nNode is 0 for UFF
        // nDOFs.w = 0; // nNode is 0 for UFF
        int4 npbc_int4;
        v2i4(nPBC_, npbc_int4);
        int err = 0;
        err |= _useArg(nDOFs);                     OCL_checkError(err, "setup_getNonBond_ex2.arg 1");      // 1
        err |= useArgBuff(ibuff_apos);             OCL_checkError(err, "setup_getNonBond_ex2.arg 2");      // 2
        err |= useArgBuff(ibuff_fapos);            OCL_checkError(err, "setup_getNonBond_ex2.arg 3");      // 3
        err |= useArgBuff(ibuff_REQs);             OCL_checkError(err, "setup_getNonBond_ex2.arg 4");      // 4
        err |= useArgBuff(ibuff_excl);             OCL_checkError(err, "setup_getNonBond_ex2.arg 5");      // 5
        err |= useArgBuff(ibuff_lvecs);            OCL_checkError(err, "setup_getNonBond_ex2.arg 6");      // 6
        err |= _useArg(npbc_int4);                 OCL_checkError(err, "setup_getNonBond_ex2.arg 7");      // 7
        err |= _useArg(GFFparams);                 OCL_checkError(err, "setup_getNonBond_ex2.arg 8");      // 8
        OCL_checkError(err, "setup_getNonBond_ex2_UFF");
        return task;
    }

    OCLtask* setup_getNonBond_GridFF_Bspline(int na, Vec3i nPBC_, OCLtask* task = nullptr) {
        if (!task) task = getTask("getNonBond_GridFF_Bspline");
        int nloc = 32;
        task->local.x = nloc;
        task->global.x = na + nloc - (na % nloc);
        task->local.y = 1;  
        task->global.y = nSystems;
        grid_shift0_p0 = grid_p0;
        useKernel(task->ikernel);
        nDOFs = (int4){na,0,0,bNonBond?0:-1};
        //nDOFs.y = 0; // nNode is 0 for UFF
        int4 npbc_int4;
        v2i4(nPBC_, npbc_int4);
        int err = 0;
        
        // Debug check for BsplinePLQ buffer
        printf("DEBUG OCL_UFF::setup_getNonBond_GridFF_Bspline() ibuff_BsplinePLQ=%d\n", ibuff_BsplinePLQ);
        if(ibuff_BsplinePLQ < 0) {
            printf("ERROR: ibuff_BsplinePLQ not initialized (%d)\n", ibuff_BsplinePLQ);
            exit(1);
        }
        
        err |= _useArg(nDOFs);                    OCL_checkError(err, "setup_getNonBond_GridFF_Bspline.arg 1");       // 1
        err |= useArgBuff(ibuff_apos);            OCL_checkError(err, "setup_getNonBond_GridFF_Bspline.arg 2");       // 2
        err |= useArgBuff(ibuff_fapos);           OCL_checkError(err, "setup_getNonBond_GridFF_Bspline.arg 3");       // 3
        err |= useArgBuff(ibuff_REQs);            OCL_checkError(err, "setup_getNonBond_GridFF_Bspline.arg 4");       // 4
        err |= useArgBuff(ibuff_neighs);          OCL_checkError(err, "setup_getNonBond_GridFF_Bspline.arg 5");       // 5
        err |= useArgBuff(ibuff_neighCell);       OCL_checkError(err, "setup_getNonBond_GridFF_Bspline.arg 6");       // 6
        err |= useArgBuff(ibuff_lvecs);           OCL_checkError(err, "setup_getNonBond_GridFF_Bspline.arg 7");       // 7
        err |= _useArg(npbc_int4);                OCL_checkError(err, "setup_getNonBond_GridFF_Bspline.arg 8");       // 8
        err |= _useArg(GFFparams);                OCL_checkError(err, "setup_getNonBond_GridFF_Bspline.arg 9");       // 9
        err |= useArgBuff(ibuff_BsplinePLQ);      OCL_checkError(err, "setup_getNonBond_GridFF_Bspline.arg 10");      // 10
        err |= _useArg(grid_n);                   OCL_checkError(err, "setup_getNonBond_GridFF_Bspline.arg 11");      // 11
        err |= _useArg(grid_invStep);             OCL_checkError(err, "setup_getNonBond_GridFF_Bspline.arg 12");      // 12
        err |= _useArg(grid_p0);                  OCL_checkError(err, "setup_getNonBond_GridFF_Bspline.arg 13");      // 13
        //err |= useArgBuff(ibuff_pbcshifts);       OCL_checkError(err, "setup_getNonBond_GridFF_Bspline.arg 14");      // 14
        
        // Debug print to verify parameters are being passed correctly
        printf("DEBUG OCL_UFF::setup_getNonBond_GridFF_Bspline() passing parameters:\n");
        printf("  grid_n=(%d,%d,%d,%d)\n", grid_n.x, grid_n.y, grid_n.z, grid_n.w);
        printf("  grid_invStep=(%f,%f,%f)\n", grid_invStep.f.x, grid_invStep.f.y, grid_invStep.f.z);
        printf("  grid_p0=(%f,%f,%f)\n", grid_p0.f.x, grid_p0.f.y, grid_p0.f.z);
        
        OCL_checkError(err, "setup_getNonBond_GridFF_Bspline_UFF");
        return task;
    }

    OCLtask* setup_getNonBond_GridFF_Bspline_ex2(int na, Vec3i nPBC_, OCLtask* task = nullptr) {
        if (!task) task = getTask("getNonBond_GridFF_Bspline_ex2");
        int nloc = 32;
        task->local.x = nloc;
        task->global.x = na + nloc - (na % nloc);
        task->local.y = 1;  
        task->global.y = nSystems;
        grid_shift0_p0 = grid_p0;
        useKernel(task->ikernel);
        nDOFs = (int4){na,0,0,bNonBond?0:-1};
        int4 npbc_int4;
        v2i4(nPBC_, npbc_int4);
        int err = 0;

        if(ibuff_BsplinePLQ < 0) {
            printf("ERROR: ibuff_BsplinePLQ not initialized (%d)\n", ibuff_BsplinePLQ);
            exit(1);
        }

        err |= _useArg(nDOFs);                    OCL_checkError(err, "setup_getNonBond_GridFF_Bspline_ex2.arg 1");       // 1
        err |= useArgBuff(ibuff_apos);            OCL_checkError(err, "setup_getNonBond_GridFF_Bspline_ex2.arg 2");       // 2
        err |= useArgBuff(ibuff_fapos);           OCL_checkError(err, "setup_getNonBond_GridFF_Bspline_ex2.arg 3");       // 3
        err |= useArgBuff(ibuff_REQs);            OCL_checkError(err, "setup_getNonBond_GridFF_Bspline_ex2.arg 4");       // 4
        err |= useArgBuff(ibuff_excl);            OCL_checkError(err, "setup_getNonBond_GridFF_Bspline_ex2.arg 5");       // 5
        err |= useArgBuff(ibuff_lvecs);           OCL_checkError(err, "setup_getNonBond_GridFF_Bspline_ex2.arg 6");       // 6
        err |= _useArg(npbc_int4);                OCL_checkError(err, "setup_getNonBond_GridFF_Bspline_ex2.arg 7");       // 7
        err |= _useArg(GFFparams);                OCL_checkError(err, "setup_getNonBond_GridFF_Bspline_ex2.arg 8");       // 8
        err |= useArgBuff(ibuff_BsplinePLQ);      OCL_checkError(err, "setup_getNonBond_GridFF_Bspline_ex2.arg 9");       // 9
        err |= _useArg(grid_n);                   OCL_checkError(err, "setup_getNonBond_GridFF_Bspline_ex2.arg 10");      // 10
        err |= _useArg(grid_invStep);             OCL_checkError(err, "setup_getNonBond_GridFF_Bspline_ex2.arg 11");      // 11
        err |= _useArg(grid_p0);                  OCL_checkError(err, "setup_getNonBond_GridFF_Bspline_ex2.arg 12");      // 12
        OCL_checkError(err, "setup_getNonBond_GridFF_Bspline_ex2_UFF");
        return task;
    }

    // Setup and bind arguments for updateAtomsMMFFf4 (UFF)
    // Must match UFF.cl::updateAtomsMMFFf4 signature exactly:
    // 0: const int4 n
    // 1: __global float4* apos
    // 2: __global float4* avel
    // 3: __global float4* aforce   (we pass ibuff_fapos)
    // 4: __global float4* cvf
    // 5: __global float4* constr
    // 6: __global float4* constrK
    // 7: __global float4* MDparams
    // 8: __global float4* TDrives
    OCLtask* setup_updateAtomsMMFFf4(int natoms, int nNode=0){
        if(!task_updateAtoms){ task_updateAtoms = getTask("updateAtomsMMFFf4"); }
        if(!task_updateAtoms) return nullptr;
        int nloc = 32;
        int nvec = natoms + nNode; // for UFF, nNode==0 typically
        task_updateAtoms->local.x  = nloc;
        task_updateAtoms->global.x = nvec + nloc - (nvec % nloc);
        task_updateAtoms->global.y = nSystems;
        useKernel( task_updateAtoms->ikernel );
        // Pack dimensions (natoms, nnode, 0, nMaxSysNeighs)
        nDOFs.x = natoms; nDOFs.y = nNode; nDOFs.z = 0; nDOFs.w = 0;
        int err=0;
        err |= _useArg   ( nDOFs        ); OCL_checkError(err, "setup_updateAtomsMMFFf4.arg 0"); // 0
        err |= useArgBuff( ibuff_apos   ); OCL_checkError(err, "setup_updateAtomsMMFFf4.arg 1"); // 1 positions
        err |= useArgBuff( ibuff_avel   ); OCL_checkError(err, "setup_updateAtomsMMFFf4.arg 2"); // 2 velocities
        err |= useArgBuff( ibuff_fapos  ); OCL_checkError(err, "setup_updateAtomsMMFFf4.arg 3"); // 3 forces (aforce)
        err |= useArgBuff( ibuff_cvf    ); OCL_checkError(err, "setup_updateAtomsMMFFf4.arg 4"); // 4 cvf accumulators
        err |= useArgBuff( ibuff_constr ); OCL_checkError(err, "setup_updateAtomsMMFFf4.arg 5"); // 5 constraints target
        err |= useArgBuff( ibuff_constrK); OCL_checkError(err, "setup_updateAtomsMMFFf4.arg 6"); // 6 constraint stiffness
        err |= useArgBuff( ibuff_MDpars ); OCL_checkError(err, "setup_updateAtomsMMFFf4.arg 7"); // 7 MD params per system
        err |= useArgBuff( ibuff_TDrive ); OCL_checkError(err, "setup_updateAtomsMMFFf4.arg 8"); // 8 thermal driving per system
        OCL_checkError(err, "OCL_UFF::setup_updateAtomsMMFFf4");
        return task_updateAtoms;
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
        if (bUFF_bonds      ){ task_evalBonds     ->enque(); }
        if (bUFF_angles     ){ task_evalAngles    ->enque(); }
        if (bUFF_dihedrals  ){ task_evalDihedrals ->enque(); }
        if (bUFF_inversions ){ task_evalInversions->enque(); }
        if (bUFF_assemble   ){ task_assemble      ->enque(); }
        //printf("OCL_UFF::eval() DONE\n");
    }

    void download_results(float* fapos_host, float* energies_host = nullptr) {
        download(ibuff_fapos, fapos_host);
        if (energies_host) {
            download(ibuff_energies, energies_host);
        }
    }
};

#endif
