/**
 * @file hubbard_solver.cl
 * @brief OpenCL kernel for brute-force ground state search of a Hubbard/Ising-like Hamiltonian.
 *
 * This kernel finds the minimum energy many-body state for a system of interacting sites on a substrate,
 * influenced by an STM tip at various positions.
 *
 * Parallelization Strategy:
 * - Each work-item (thread) corresponds to a single STM tip position and voltage.
 * - The global work size should be equal to `nTips`.
 *
 * Optimization:
 * - The state of the many-body system (occupancy of `nSingle` sites) is encoded as bits in a single `uint`.
 * - The site-site interaction matrix `W_ij` is pre-calculated collaboratively by the work-group and stored
 *   in fast __local memory to accelerate the main energy calculation loop.
 * - Tip-dependent potentials `U_i` and tunneling factors `T_i` are pre-calculated per thread.
 *
 * Hamiltonian Model Assumed:
 *   E = sum_i{ n_i * (E0_i + V_i) } + sum_{i<j}{ n_i * n_j * W_ij }
 *   where:
 *     n_i:   Occupancy of site i (0 or 1)
 *     E0_i:  On-site energy of site i
 *     V_i:   Potential at site i from the tip -> Vbias * (tipRadius / |r_i - r_tip|)
 *     W_ij:  Interaction energy between sites i and j.
 *
 * Current Model:
 *   I_occ  = sum_i{ n_i * exp(-b * |r_i - r_tip|) }
 *   I_unoc = sum_i{ (1 - n_i) * exp(-b * |r_i - r_tip|) }
 *   Calculated for the found ground state.
 */

// Define max number of sites supported. This must be a compile-time constant.
// 32*32*4 bytes = 4KB, which should fit in the local memory of most GPUs.
#define MAX_WORKGROUP_SIZE 32
#define MAX_SITES 32
#define MAX_WIJ (MAX_SITES * (MAX_SITES - 1) / 2)
#define COULOMB_CONST 14.399644730092272


float Vtip( float3 pSite, float3 pTip, float Vbias, float Rtip ){
    // Tip potential V_i = Vbias * (R_tip / |r_i - r_tip|)
    const float3 d   = pSite - pTip;
    const float  r   = length(d);
    return  Vbias * (Rtip/r); // Tip potential V_i = Vbias * (R_tip / |r_i - r_tip|)
}


float get_Wij( float3 pi, float3 pj ){
    // --- TODO: PHYSICAL MODEL FOR W_ij ---
    // The interaction W_ij depends on your specific physical system (e.g., screening, dipole type).
    // Replace the simple 1/r model below with your actual formula.
    // For dipole-dipole interaction, it might be proportional to 1/r^3.
    // For a simple screened Coulomb interaction, it could be: prefactor * exp(-dist / lambda) / dist.
    const float3  d = pi - pj;
    const float   r = length(d);
    return  COULOMB_CONST/r;
}

float get_Wij_mirror( float3 pi, float3 pj, float zMirror ){
    const float3 pj_mirror = (float3)(pj.x, pj.y, 2.0f * zMirror - pj.z);
    const float  r_rr      = distance(pi, pj);
    const float  r_ri      = distance(pi, pj_mirror);
    return 2.0f * COULOMB_CONST * ( (1.0f / r_rr) - (1.0f / r_ri) );
}


int get_W_idx(int i, int j, int n) {
    // This formula is derived from the sum of an arithmetic series for the preceding rows.
    // Sum of lengths of rows 0..i-1 is: i*n - i*(i+1)/2
    // The offset in the current row i is (j - i - 1).
    return (i*n-i * (i+1)/2) + (j-i-1);
}

/**
 * @file hubbard_solver_on_the_fly.cl
 * @brief Kernel using on-the-fly W_ij calculation to conserve local memory.
 *
 * This kernel uses the "one work-group per tip" model but avoids storing the
 * large W_ij matrix in local memory. Instead, W_ij values are computed on-the-fly
 * as needed inside the main loop. To make this efficient, site positions are
 * pre-loaded into local memory.
 *
 * This strategy is optimal when local memory is too limited to hold the full W_ij matrix.
 *
 * Parallelization Strategy:
 * - tip_idx      = get_group_id(0)
 * - lid          = get_local_id(0)
 * - Global Size  = nTips * workgroup_size
 * - Local Size   = workgroup_size
 */


__kernel void eval_coupling_matrix(
    const int nSingle,
    __global const float4* posE,
    __global float*  Rij,
    __global float*  Wij,
    const float zMirror
){
    const int gid    = get_global_id(0);
    //cnst int lid    = get_local_id(0);
    //const int local_size = get_local_size(0);

    if (gid >= nSingle*nSingle ) return;
    //if (lid >= nSingle         ) return;

    const int i = gid / nSingle;
    const int j = gid % nSingle;

    if (i == j) {
        Rij[gid] = 0.0f;
        Wij[gid] = 0.0f;
        return;
    }

    const float4 site  = posE[i];
    const float4 site2 = posE[j];
    const float3 d     = site.xyz - site2.xyz;
    const float  r     = length(d);
    Rij[gid] = r;
    Wij[gid] = get_Wij_mirror(site.xyz, site2.xyz, zMirror);
}

__kernel void eval_Oriented_Hopping(
    const int nSingle,
    const int nTips,
    __global const float4* posE,   // [nSingle] {x,y,z,?}        site positions
    __global float2*       rots,    // [nSingle] {cos, sin}       site rotation as unit complex number
    __global float4*       orbs, // [nSingle] {s,px,py,decay}  tip wavefunctions
    __global const float4* pTips,  // [nTips]   {x,y,z,?}        tip positions
    __global float*        Tout   // [nTips*nSingle]            output tunneling amplitudes
){
    const int gid    = get_global_id(0);
    if (gid >= nTips*nSingle ) return;
    const int iTip   = gid / nSingle;
    const int iSite  = gid % nSingle;
    const float3 d = pTips[iTip].xyz - posE[iSite].xyz;
    const float  r = sqrt( dot(d,d) + 1e-9f );
    const float2 rot   = rots[iSite];      // {cos(a), sin(a)}
    const float da =  d.x * rot.x + d.y * rot.y;
    const float db = -d.x * rot.y + d.y * rot.x;
    const float inv_r = 1.0f / r;
    const float4 orb   = orbs[iSite];   // {C_s, C_px, C_py, decay_const}
    const float wfa =
          orb.x              // s-wave part (isotropic)
        + orb.y * da * inv_r // px'-wave part
        + orb.z * db * inv_r;// py'-wave part
    const float wfr = exp(-orb.w * r);
    Tout[gid] = wfr * wfa;
}

__kernel void solve_minBrute_fly(
    const int nSingle,
    __global const float4* posE,
    const int nTips,
    __global const float4* pTips,
    const float tipDecay,
    const float tipRadius,
    const float zMirror,
    const float Wcouple,
    __global float*  Emin,
    __global int*    iMin,
    __global float2* Itot
) {
    //=========================================================================
    // 1. Initialization and ID setup
    //=========================================================================
    const int iTip       = get_group_id(0);
    const int lid        = get_local_id(0);
    const int gid        = get_global_id(0);
    const int local_size = get_local_size(0);

    if (iTip >= nTips) return;
    if (nSingle > MAX_SITES) { if (lid == 0) Emin[iTip] = 1e10f; return; }

    bool bNoCrossCoupling = false;
    if( fabs(Wcouple)<1e-9 ) bNoCrossCoupling = true;

    // Shared memory for this work-group. Note the absence of the large 'W' matrix.
    __local float4 Lsites[MAX_SITES];
    //__local float4 Tsite [MAX_SITES];

    #define iDBG 0

    //if (iDBG==gid){  for (int i = 0; i < nSingle; ++i) { printf("GPU posE %4i E: %12.8f  pos( %12.8f , %12.8f , %12.8f ) \n", i, posE[i].w, posE[i].x, posE[i].y, posE[i].z);  } }

    //=========================================================================
    // 2. Collaborative Pre-calculation in Local Memory
    //=========================================================================
    // All threads in the workgroup collaborate to fill the shared arrays.
    const float4 tip = pTips[iTip];
    for (int i = lid; i < nSingle; i += local_size) {
        Lsites[i].xyz = posE[i].xyz;
        Lsites[i].w   = posE[i].w + Vtip(Lsites[i].xyz, tip.xyz, tip.w, tipRadius);
    }
    
    // Synchronize to ensure all shared arrays are fully computed.
    barrier(CLK_LOCAL_MEM_FENCE);

    //if (iDBG==gid){  for (int i = 0; i < nSingle; ++i) { printf("GPU Lsite %4i E: %12.8f  pos( %12.8f , %12.8f , %12.8f ) \n", i, Lsites[i].w, Lsites[i].x, Lsites[i].y, Lsites[i].z); } }

    //=========================================================================
    // 3. Parallel Brute-Force Minimization with On-the-Fly W_ij
    //=========================================================================
    const uint nMany = 1 << nSingle;
    
    float lEmin = 1e10f;
    uint  liMin = 0;

    for (uint iMany = lid; iMany < nMany; iMany += local_size) {
        float E = 0.0f;
        for (int i = 0; i < nSingle; ++i) {  // On-site and interaction energy calculation
            if ((iMany >> i) & 1) {          // If site 'i' is occupied...
                const float4 site = Lsites[i];             // Add its on-site energy
                E += site.w;
                if(bNoCrossCoupling) continue;
                for (int j = i + 1; j < nSingle; ++j) { // ...calculate interaction with other occupied sites j > i
                    if ((iMany >> j) & 1) {  E += get_Wij_mirror(site.xyz, Lsites[j].xyz, zMirror ) * Wcouple; }
                }
            }
        }
        //printf("GPU Emany %4i lid %2i  E: %12.8f  lEmin %12.8f \n", iMany, lid, E, lEmin);
        if (E < lEmin) {
            //printf("GPU Emany-min %4i  lid %2i E: %12.8f ->  lEmin %12.8f \n", iMany, lid, E, lEmin);
            lEmin = E;
            liMin = iMany;
        }
    }

    //=========================================================================
    // 4. Work-group Reduction to find the final minimum
    //=========================================================================
    __local float LEmins[MAX_WORKGROUP_SIZE];
    __local uint  LImins[MAX_WORKGROUP_SIZE];

    LEmins[lid] = lEmin;
    LImins[lid] = liMin;
    barrier(CLK_LOCAL_MEM_FENCE);

    // Thread 0 performs the final reduction and writes the result.
    if (lid == 0) {
        float gEmin = 1e10f;
        uint  gimin = -1;
        int nred = (nMany<MAX_WORKGROUP_SIZE) ? nMany : MAX_WORKGROUP_SIZE;
        for (int i=0; i<nred; ++i) {
            const float LEmin = LEmins[i];
            //if (iDBG==gid){ printf("GPU gMin %4i i %2i LEmin: %12.8f imin: %i \n", iTip, i, LEmin, LImins[i]); }
            if (LEmin<gEmin) {
                gEmin = LEmin;
                gimin = LImins[i];
            }
        }
        //=====================================================================
        // 5. Final Calculations & Output (done only by thread 0)
        //=====================================================================
        //if (iDBG==gid){ printf("GPU final gMin %4i Emin: %12.8f imin: %i \n", iTip, gEmin, gimin); }
        Emin[iTip] = gEmin;
        iMin[iTip] = (int)gimin;

        float I_occ  = 0.0f;
        float I_unoc = 0.0f;
        for (int i=0; i<nSingle; ++i) {  // Check occupancy in the found ground state (min_idx)
            float r = distance(tip.xyz, Lsites[i].xyz);
            float T = exp(-tipDecay * r);
            if ((gimin >> i) & 1){ I_occ += T; }else{ I_unoc += T; }
        }
        Itot[iTip] = (float2)(I_occ, I_unoc);
    }
}


/**
 * @file hubbard_solver_boltzmann_single_pass.cl
 * @brief Highly efficient single-pass kernel for Boltzmann-weighted averaging.
 *
 * This kernel avoids the two-pass approach by using a running correction algorithm.
 * It maintains a running E_min and rescales the accumulated partition function (Z)
 * and weighted currents whenever a new, lower minimum energy is found.
 * This is numerically stable and computationally efficient.
 *
 * Algorithm per Work-Group:
 * 1. Pre-load shared data into local memory.
 * 2. Each thread performs a single-pass search over its subset of states,
 *    maintaining its own local_Emin, local_Z, and local weighted sums.
 * 3. A final reduction step combines the results from all threads. It finds the
 *    true group_Emin and rescales each thread's partial sums to this common
 *    reference energy before summing them to get the final result.
 */

__kernel void solve_minBrute_boltzmann(
    const int nSingle,
    __global const float4* posE,
    const int nTips,
    __global const float4* pTips,
    const float tipDecay,
    const float tipRadius,
    const float zMirror,
    const float Wcouple,
    const float bBoltzmann, // 1/kT [1/eV]
    __global float*  Emin,
    __global int*    iMin,
    __global float2* Itot
) {
     //=========================================================================
    // 1. Initialization and ID setup
    //=========================================================================
    const int tip_idx    = get_group_id(0);
    const int lid        = get_local_id(0);
    const int local_size = get_local_size(0);

    if (tip_idx >= nTips) return;
    if (nSingle > MAX_SITES) { if (lid == 0) Emin[tip_idx] = 1e10f; return; }

    bool bNoCrossCoupling = false;
    if( fabs(Wcouple)<1e-9 ) bNoCrossCoupling = true;

    // Shared memory for this work-group. Note the absence of the large 'W' matrix.
    __local float4 Lsites[MAX_SITES];
    __local float  Tsite [MAX_SITES];

    //=========================================================================
    // 2. Collaborative Pre-calculation in Local Memory
    //=========================================================================
    // All threads in the workgroup collaborate to fill the shared arrays.
    const float4 tip = pTips[tip_idx];
    for (int i = lid; i < nSingle; i += local_size) {
        Lsites[i].xyz = posE[i].xyz;
        Lsites[i].w   = posE[i].w + Vtip(Lsites[i].xyz, tip.xyz, tip.w, tipRadius);
        Tsite[i]      = exp(-tipDecay * distance(tip.xyz, Lsites[i].xyz));
    }
    
    // Synchronize to ensure all shared arrays are fully computed.
    barrier(CLK_LOCAL_MEM_FENCE);

    //=========================================================================
    // 3. SINGLE-PASS Boltzmann-weighted summation
    //=========================================================================
    const uint nMany = 1u << nSingle;
    
    // Each thread finds the minimum and sums in its own subset of states
    float local_Emin   = 1e10f;
    uint  local_imin   = 0;
    float local_Z      = 0.0f;
    float local_Iocc   = 0.0f;
    float local_Iunoc  = 0.0f;

    for (uint iMany = lid; iMany < nMany; iMany += local_size) {
        float E = 0.0f;
        float i_occ_state = 0.0f;
        float i_unoc_state = 0.0f;

        for (int i = 0; i < nSingle; ++i) {
            float T = Tsite[i];
            if ((iMany >> i) & 1) {
                float4 site = Lsites[i];
                E += site.w;
                i_occ_state +=  T;
                if(bNoCrossCoupling) continue;
                for (int j = i + 1; j < nSingle; ++j) {
                    if ((iMany >> j) & 1) {  E += get_Wij_mirror(site.xyz, Lsites[j].xyz, zMirror ) * Wcouple; }
                }
            } else { i_unoc_state += T; }
        }

        // boltzmann-weighted average
        float w = exp( - E * bBoltzmann );
        local_Z         *= w;
        local_Iocc      *= w;
        local_Iunoc     *= w;
        if (E < local_Emin) {
            // New minimum found! Rescale previous sums to the new reference energy.
            local_Emin = E;
            local_imin = iMany;
        }
    }

    //=========================================================================
    // 4. Final Reduction and Output
    //=========================================================================
    __local float reduction_E    [MAX_WORKGROUP_SIZE];
    __local uint  reduction_idx  [MAX_WORKGROUP_SIZE];
    __local float reduction_Z    [MAX_WORKGROUP_SIZE];
    __local float reduction_Iocc [MAX_WORKGROUP_SIZE];
    __local float reduction_Iunoc[MAX_WORKGROUP_SIZE];

    reduction_E[lid]     = local_Emin;
    reduction_idx[lid]   = local_imin;
    reduction_Z[lid]     = local_Z;
    reduction_Iocc[lid]  = local_Iocc;
    reduction_Iunoc[lid] = local_Iunoc;
    barrier(CLK_LOCAL_MEM_FENCE);

    if (lid == 0) {
        // --- First, find the true group-wide minimum energy ---
        float group_Emin   = reduction_E[0];
        uint  group_imin = reduction_idx[0];
        for (int i = 1; i < local_size; ++i) {
            if (reduction_E[i] < group_Emin) {
                group_Emin = reduction_E[i];
                group_imin = reduction_idx[i];
            }
        }
        // --- Now, sum up all contributions, rescaling them to the true group_Emin ---
        float total_Z = 0.0f, total_Iocc = 0.0f, total_Iunoc = 0.0f;
        for (int i = 0; i < local_size; ++i) {
            total_Z     += reduction_Z[i];
            total_Iocc  += reduction_Iocc[i];
            total_Iunoc += reduction_Iunoc[i];
        }

        // Store debug ground state info
        Emin[tip_idx] =      group_Emin;
        iMin[tip_idx] = (int)group_imin;

        // Calculate final averaged currents and store them
        if (total_Z > 1e-30f) { // Avoid division by zero
            Itot[tip_idx] = (float2)(total_Iocc / total_Z, total_Iunoc / total_Z);
        } else {
            Itot[tip_idx] = (float2)(0.0f, 0.0f);
        }
    }
}




/**
 * @file parallel_local_update_solver.cl
 * @brief Advanced solver using work-group collaboration and bitmasked state
 *        for large-scale Monte Carlo simulations.
 *
 * This kernel assigns an entire work-group to a single tip position to accelerate
 * the search for a low-energy state. It is highly memory-efficient due to bitmasking.
 *
 * Parallel Strategy (per work-group):
 * 1. The current state is stored in a __local bitmask array, shared by all threads.
 * 2. In each iteration, all threads in parallel propose a random site-flip.
 * 3. Each thread calculates the energy change (dE) for its proposed flip.
 * 4. A parallel reduction finds the "best" move among all proposals (the one with the
 *    lowest dE).
 * 5. One thread updates the shared state with this single best move.
 * 6. Barriers are used to synchronize the work-group between steps.
 * This massively parallelizes the "search" phase of the Monte Carlo step.
 */


#define MAX_WORKGROUP_SIZE_BIG 16
#define MAX_SITES_BIG 256
#define OCC_BYTES     32
//#define MAX_NEIGHBORS 8

//#define MAX_WORKGROUP_SIZE 256 // For sizing reduction arrays

// --- Bitmask Helper Macros ---
// Calculates the value of the i-th bit in the uchar array 'mask'
#define GET_OCC(i, mask) (((mask)[(i) >> 3] >> ((i) & 7)) & 1)

// Flips the i-th bit in the uchar array 'mask' using XOR
#define FLIP_OCC(i, mask) ((mask)[(i) >> 3] ^= (1 << ((i) & 7)))

// --- RNG Helper Functions ---
uint wang_hash(uint seed) {
    seed = (seed ^ 61u) ^ (seed >> 16u); seed *= 9u;
    seed = seed ^ (seed >> 4u); seed *= 0x27d4eb2du;
    seed = seed ^ (seed >> 15u); return seed;
}

uint wang_hash_uint(uint* seed) { *seed = wang_hash(*seed); return *seed; }

float wang_hash_float(uint* seed) {
    *seed = wang_hash(*seed);
    return ((float)(*seed)) / ((float)UINT_MAX);
}

// Helper function to calculate total energy of a configuration
static float calculate_total_energy(
    __local uchar* occ_mask,
    const int nSite,
    __global const float* Esite,
    __global const float* W_val,
    __global const int* W_idx,
    __global const int* nNeigh,
    const int max_neighs,
    const int tip_offset
) {
    float E = 0.0f;
    for (int i = 0; i < nSite; ++i) {
        if (GET_OCC(i, occ_mask)) {
            E += Esite[tip_offset + i];
            if (max_neighs > 0) {
                const int iw0 = i * max_neighs;
                for (int k = 0; k < nNeigh[i]; ++k) {
                    const int iw = iw0 + k;
                    const int j = W_idx[iw];
                    if (i < j && GET_OCC(j, occ_mask)) {
                        E += W_val[iw];
                    }
                }
            }
        }
    }
    return E;
}



__kernel void solve_local_updates(
    const int nSite, 
    const int nTips, 
    const float4 solver_params,
    __global const float*  Esite,
    __global const float*  Tsite,
    __global const float*  W_val, 
    __global const int*    W_idx, 
    __global const int*    nNeigh,
    __global uchar*        occ_out, 
    __global float*        E_out, 
    __global float2*       Itot_out,
    const int              max_neighs,
    int                    initMode
) {
    //=========================================================================
    // 1. Initialization and ID setup
    //=========================================================================
    const int itip       = get_group_id(0);
    const int lid        = get_local_id(0);
    const int local_size = get_local_size(0);

    const int gid  = get_global_id(0);
    
    if (itip >= nTips) return;
    if (nSite > MAX_SITES_BIG || local_size > MAX_WORKGROUP_SIZE_BIG) return;

    // Unpack solver parameters
    const float kT         = solver_params.x;
    const int   nIter      = (int)(solver_params.y+0.5);
    const int   solverMode = (int)(solver_params.z+0.5);
    uint rng_state         = (uint)solver_params.w + itip * local_size + lid; // Unique seed per thread
    rng_state=wang_hash(rng_state);


    // --- Shared memory for the work-group ---
    // State is stored as a bitmask to save memory
    const int occ_bytes        = min( (nSite+7)/8, OCC_BYTES );
    __local uchar occ_mask      [OCC_BYTES];
    __local float reduction_dE  [MAX_WORKGROUP_SIZE_BIG];
    __local int   reduction_site[MAX_WORKGROUP_SIZE_BIG];
    
    const int tip_offset = itip * nSite;   // Offset for accessing tip-specific global data

    //#define iDBG 612
    // for shape(200,200)  30100 = 200*200*(0.75) + 200/2
    #define iDBG 30100  

    if((itip==iDBG)&&(lid==0)){ 
        int nG = get_global_size(0);
        printf("GPU solve_local_updates() iDBG %i  nSite: %i nTips: %i nIter %i max_neighs %i solverMode %i kT %8.4e initMode %i local_size: %i global_size %i  occ_bytes %i\n", iDBG, nSite, nTips, nIter, max_neighs, solverMode, kT, initMode, local_size, nG, occ_bytes ); 

        // printf("GPU isite,Esite, Tsite: \n");
        // for( int is=0; is<nSite; ++is){ printf("GPU isite %3i Esite %16.8f Tsite %16.8f\n", is, Esite[tip_offset + is], Tsite[tip_offset + is] ); }

        // printf("GPU Wij (sparse inter-site couplings): \n");
        // for( int is=0; is<nSite; ++is){
        //     int nng = nNeigh[is];
        //     printf("GPU site %3i nng %i Wij: ", is, nng );
        //     for( int j=0; j<nng; ++j){
        //         int iw = is * max_neighs + j;
        //         printf(" %3i:%10.6f ", W_idx[iw], W_val[iw] );
        //     }
        //     printf("\n");
        // }
    }


    // Thread 0 initializes the state to all zeros
    if (lid == 0) {
        //initMode = 1; // DEBUG - override initMode
        if     (initMode == 0 ){ for (int i=0;i<occ_bytes;++i) { occ_mask[i] = 0;    } } 
        else if(initMode == 1 ){ for (int i=0;i<occ_bytes;++i) { occ_mask[i] = 0xFF; } } 
        else if(initMode == 2 ){ for (int i=0;i<occ_bytes;++i) { occ_mask[i] = occ_out[ itip*OCC_BYTES + i]; } } 
        else if(initMode == 3 ){ for (int i=0;i<occ_bytes;++i) { occ_mask[i] = wang_hash_uint(&rng_state)&0xFF; } }
        //if(itip==iDBG){ printf("GPU itip %3i occ: ", itip ); for(int i=0;i<occ_bytes;++i){ printf("%02x", occ_mask[i]); } printf("\n"); }
    }
    barrier(CLK_LOCAL_MEM_FENCE); // Ensure all threads see the initialized state

    //=========================================================================
    // 2. Main Collaborative Monte Carlo Loop
    //=========================================================================
    for (int iter = 0; iter < nIter; ++iter) {
        // --- Step A: Each thread proposes a move in parallel ---
        int i_site;
        if (solverMode == 0){ i_site = (iter * local_size + lid ) % nSite; } // Deterministic sweep across threads
        else                { i_site = wang_hash_uint(&rng_state) % nSite; } // Stochastic proposal

        const uchar n_i  = GET_OCC(i_site, occ_mask);     // site occupancy
        float      Ei    = Esite[tip_offset + i_site];    // on-site energy
        
        if(max_neighs>0){ // inter-site coupling - only if there are neighbors
            const int  iw0   = i_site * max_neighs;
            for (int k = 0; k < nNeigh[i_site]; ++k) {        // coulomb interactions with neighbors
                const int iw    = iw0 + k;
                const int jsite = W_idx[iw];
                float occ_k  = (float)GET_OCC(jsite, occ_mask); // neighbor occupation
                Ei          += occ_k * W_val[iw];
            }
        }

        //reduction_dE  [lid] = (2.0f * (float)n_i - 1.0f) * Ei; // Corrected energy change
        if( n_i>0){ reduction_dE[lid] = -Ei; } // if occupied, we make it unoccupied, so we substrat site energy from total sum
        else      { reduction_dE[lid] = +Ei; } // if unoccupied, we make it occupied, so we add site energy to total sum


        // if((itip==iDBG)&&(lid==0)){   printf("GPU itip %3i iter %3i i_site %3i n_i %i Ei %16.8f dE %16.8f \n", itip, iter, i_site, n_i, Ei, reduction_dE[lid] ); }

        reduction_site[lid] = i_site;

        barrier(CLK_LOCAL_MEM_FENCE); // Wait for all proposals to be calculated

        // --- Step B: Thread 0 finds the best move and decides ---
        if (lid == 0) {
            float dEmin = 1e10f; // Start with a high energy
            int   imin  = -1;
            for (int il = 0; il < local_size; ++il) {  // Find the move with the minimum dE among all threads
                //if(itip==iDBG){  printf("GPU ibest? iter %3i itip %3i lid %3i i_site %3i dE %16.8f \n", iter, itip, il, reduction_site[il], reduction_dE[il] ); }
                if (reduction_dE[il] < dEmin) {
                    dEmin = reduction_dE[il];
                    imin  = reduction_site[il];
                }
            }
            if(imin < 0){ continue; }
            bool bDo = false;
            if      (dEmin < 0.0f) { bDo = true; } 
            else if ( solverMode == 2  ) { // Simulated Annealing
                //if (wang_hash_float(&rng_state) < exp(-dEmin / kT)) { bDo = true; }
            }
            //if(itip==iDBG){  printf("GPU flip? iter %3i itip %3i i_site %3i dE %16.8f bDo %i\n", iter, itip, imin, dEmin, bDo ); }
            if (bDo){ FLIP_OCC(imin, occ_mask); }
        }
        barrier(CLK_LOCAL_MEM_FENCE); // Wait for state update before next iteration
    }

    //=========================================================================
    // 3. Post-Loop: Thread 0 calculates final properties and writes output
    //=========================================================================
    if (lid == 0) {
        float E     = 0.0f;
        float Iocc  = 0.0f;
        float Iunoc = 0.0f;
        for (int i = 0; i < nSite; ++i) {
            if (GET_OCC(i, occ_mask)) {
                Iocc += Tsite[tip_offset + i];          // on-site current
                E    += Esite[tip_offset + i];          // on-site energy
                if(max_neighs>0){ // inter-site coupling - only if there are neighbors
                    const int iw0 = i * max_neighs;
                    for (int k = 0; k < nNeigh[i]; ++k) {  // coulomb interactions with neighbors
                        const int iw = iw0 + k;
                        const int j  = W_idx[iw];
                        if (i < j && GET_OCC(j, occ_mask)) {  E += W_val[iw]; }
                    }
                }
            } else {
                Iunoc += Tsite[tip_offset + i];          // on-site current
            }
        }
        // if(itip==iDBG){  
        //     printf("GPU final itip %3i E %16.8f Iocc %16.8f Iunoc %16.8f \n", itip, E, Iocc, Iunoc ); 
        //     for(int i=0;i<nSite;++i){ printf("GPU occ[%3i] %3i E %16.8f \n", i, GET_OCC(i, occ_mask), Esite[tip_offset + i] ); }
        // }

        for (int j = 0; j < occ_bytes; ++j) { occ_out[itip*OCC_BYTES + j] = occ_mask[j]; } // store optimized state configuration
        E_out   [itip]     = E;                  // store energy of optimized state
        Itot_out[itip]  = (float2)(Iocc, Iunoc); // store current of optimized state
    }
}


/**
 * @brief Calculates the total energy of a configuration in parallel within a workgroup.
 * @param occ_mask Local memory buffer with the occupancy bitmask.
 * @param reduction_buffer Scratch local memory for parallel reduction.
 * @return The total energy of the configuration.
 */
static float calculate_total_energy_parallel(
    __local uchar* occ_mask,
    const int nSite,
    __global const float* Esite,
    __global const float* W_val,
    __global const int* W_idx,
    __global const int* nNeigh,
    const int max_neighs,
    const int tip_offset,
    __local float* reduction_buffer
) {
    const int lid = get_local_id(0);
    const int local_size = get_local_size(0);

    // Each thread calculates a partial sum of energy for its subset of sites
    float partial_E = 0.0f;
    for (int i = lid; i < nSite; i += local_size) {
        if (GET_OCC(i, occ_mask)) {
            partial_E += Esite[tip_offset + i];
            if (max_neighs > 0) {
                const int iw0 = i * max_neighs;
                for (int k = 0; k < nNeigh[i]; ++k) {
                    const int iw = iw0 + k;
                    const int j = W_idx[iw];
                    // i < j ensures each interaction is counted only once
                    if (i < j && GET_OCC(j, occ_mask)) {
                        partial_E += W_val[iw];
                    }
                }
            }
        }
    }
    reduction_buffer[lid] = partial_E;

    barrier(CLK_LOCAL_MEM_FENCE);

    // Parallel reduction to sum up the partial energies
    for (int offset = local_size / 2; offset > 0; offset /= 2) {
        if (lid < offset) {
            reduction_buffer[lid] += reduction_buffer[lid + offset];
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    // The final result is in thread 0's element
    return reduction_buffer[0];
}

__kernel void solve_MC(
    const int nSite,
    const int nTips,
    const float4 solver_params,
    const int nLocalIter,
    const float4 prob_params, // NEW: Probabilities {p_keep, p_reload_best, p_random_reset, unused}
    const int nRelaxSteps, // Unused
    __global const float*  Esite,
    __global const float*  Tsite,
    __global const float*  W_val,
    __global const int*    W_idx,
    __global const int*    nNeigh,
    __global uchar*        occ_best,
    __global float*        E_best,
    __global float2*       Itot_out,
    const int              max_neighs,
    int                    initMode
) {
    //=========================================================================
    // 1. Initialization and ID setup
    //=========================================================================
    const int itip       = get_group_id(0);
    const int lid        = get_local_id(0);
    const int local_size = get_local_size(0);

    if (itip >= nTips) return;

    const int nIterTotal   = (int)(solver_params.y + 0.5);
    //const int nIterTotal     = 10000;
    const int solverMode     = (int)(solver_params.z + 0.5);
    const int occ_bytes      = min((nSite + 7) / 8, OCC_BYTES);
    const int tip_offset     = itip * nSite;
    const int tip_occ_offset = itip * OCC_BYTES;

    __local uchar occ_mask         [OCC_BYTES];
    __local float reduction_buffer [MAX_WORKGROUP_SIZE_BIG];
    __local int   reduction_site   [MAX_WORKGROUP_SIZE_BIG];

    float E_best_cached;
    const int nGlobalSteps = max(1, nIterTotal / nLocalIter);

    uint rng_state = (uint)solver_params.w + get_global_id(0);
    rng_state = wang_hash(rng_state);

    if ((itip == iDBG) && (lid == 0)) {
        printf("GPU kernel solve_MC() nSite: %i, nLocalIter: %i, nIterTotal: %i iDBG %i solverMode %i initMode %i OCC_BYTES %i MAX_WORKGROUP_SIZE_BIG %i\n", nSite, nLocalIter, nIterTotal, iDBG, solverMode, initMode, OCC_BYTES, MAX_WORKGROUP_SIZE_BIG);
        printf("GPU kernel solve_MC() Probabilities: Keep=%.2f, Reload=%.2f, Random=%.2f\n", prob_params.x, prob_params.y, prob_params.z);
    }

    //=========================================================================
    // 2. Initialize State
    //=========================================================================
    if (lid == 0) {
        E_best_cached = E_best[itip];
        if (initMode == 2) {
            if (itip == iDBG) printf("iDBG: InitMode 2 - Loading from global occ_best.\n");
            for (int i = 0; i < occ_bytes; ++i) occ_mask[i] = occ_best[tip_occ_offset + i];
        } else {
            if (itip == iDBG) printf("iDBG: InitMode %i - Generating new configuration.\n", initMode);
            if (initMode == 0) { for (int i = 0; i < occ_bytes; ++i) occ_mask[i] = 0; }
            else if (initMode == 1) { for (int i = 0; i < occ_bytes; ++i) occ_mask[i] = 0xFF; }
            else { for (int i = 0; i < occ_bytes; ++i) occ_mask[i] = wang_hash_uint(&rng_state) & 0xFF; }
        }
    }
    barrier(CLK_LOCAL_MEM_FENCE);
    return;

    //=========================================================================
    // 3. Main Hybrid Optimization Loop
    //=========================================================================
    
    for (int iGlobal = 0; iGlobal < nGlobalSteps; ++iGlobal) {
        if ((itip == iDBG) && (lid == 0)) {  printf("\n--- iDBG: Global Step %i/%i, E_best_cached = %e ---\n", iGlobal, nGlobalSteps, E_best_cached); }

        // --- A: Perform a batch of greedy local updates ---
        for (int iLocal = 0; iLocal < nLocalIter; ++iLocal) {
            int i_site = (solverMode == 0) ? ((iGlobal * nLocalIter + iLocal) * local_size + lid) % nSite : (wang_hash_uint(&rng_state) % nSite);
            uchar n_i = GET_OCC(i_site, occ_mask);
            float Ei = Esite[tip_offset + i_site];
            if(max_neighs > 0){
                const int iw0 = i_site * max_neighs;
                for (int k = 0; k < nNeigh[i_site]; ++k) {
                    Ei += (float)GET_OCC(W_idx[iw0 + k], occ_mask) * W_val[iw0 + k];
                }
            }
            reduction_buffer[lid] = (n_i > 0) ? -Ei : Ei; // Use reduction_buffer for dE
            reduction_site[lid] = i_site;
            barrier(CLK_LOCAL_MEM_FENCE);

            if (lid == 0) {
                float dEmin = 1e10f;
                int imin = -1;
                for (int il = 0; il < local_size; ++il) {
                    if (reduction_buffer[il] < dEmin) {
                        dEmin = reduction_buffer[il];
                        imin = reduction_site[il];
                    }
                }
                if (imin >= 0 && dEmin < 0.0f) {
                    FLIP_OCC(imin, occ_mask);
                }
            }
            barrier(CLK_LOCAL_MEM_FENCE);
        }

        // // --- B: Evaluate, Compare, and Decide Next Step ---
        float E_current = calculate_total_energy_parallel(occ_mask, nSite, Esite, W_val, W_idx, nNeigh, max_neighs, tip_offset, reduction_buffer);
        // // TODO: we may elimiate this calculate_total_energy_parallel if we accumulate all the updates of energy in some local variable (not need to be private)

        if (lid == 0) {
            if ((itip == iDBG)) { printf("iDBG: Local optimization finished. E_current = %e, E_best_cached = %e\n", E_current, E_best_cached); }
            if (E_current < E_best_cached) {
                if (itip == iDBG) printf("iDBG: *** NEW BEST FOUND! *** (%e < %e)\n", E_current, E_best_cached);
                E_best_cached = E_current;
                E_best[itip] = E_current;
                for (int i = 0; i < occ_bytes; ++i) {
                    occ_best[tip_occ_offset + i] = occ_mask[i];
                }
            } else {
                if (itip == iDBG) printf("iDBG: No improvement. Deciding next action...\n");
                float r = wang_hash_float(&rng_state);
                if (r < prob_params.x) { // Option 1: Keep current state
                    if (itip == iDBG) printf("iDBG: Action -> KEEP current config (r=%.2f < p=%.2f)\n", r, prob_params.x);
                    // Do nothing, continue with current occ_mask
                } else if (r < prob_params.x + prob_params.y) { // Option 2: Reload global best
                    if (itip == iDBG) printf("iDBG: Action -> RELOAD global best (r=%.2f < p=%.2f)\n", r, prob_params.x + prob_params.y);
                    for (int i = 0; i < occ_bytes; ++i) {
                        occ_mask[i] = occ_best[tip_occ_offset + i];
                    }
                } else { // Option 3: Random reset
                    if (itip == iDBG) printf("iDBG: Action -> RANDOM reset\n");
                    for (int i = 0; i < occ_bytes; ++i) {
                        occ_mask[i] = wang_hash_uint(&rng_state) & 0xFF;
                    }
                }
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    //=========================================================================
    // 4. Finalization: Calculate currents for the best state
    //=========================================================================
    if (lid == 0) {
        if ((itip == iDBG)) printf("====== KERNEL FINISH (iDBG=%i) ======\n\n", iDBG);
        float Iocc = 0.0f, Iunoc = 0.0f;
        for (int j = 0; j < occ_bytes; ++j) {
            occ_mask[j] = occ_best[tip_occ_offset + j];
        }
        for (int i = 0; i < nSite; ++i) {
            if (GET_OCC(i, occ_mask)) Iocc += Tsite[tip_offset + i];
            else Iunoc += Tsite[tip_offset + i];
        }
        Itot_out[itip] = (float2)(Iocc, Iunoc);
    }
}


__kernel void solve_MC_bak(
    const int nSite,
    const int nTips,
    const float4 solver_params,
    // --- New global exploration parameters ---
    const int nLocalIter,
    const float randConfigProb,
    const int nRelaxSteps,
    // --- Global Data Buffers ---
    __global const float*  Esite,
    __global const float*  Tsite,
    __global const float*  W_val,
    __global const int*    W_idx,
    __global const int*    nNeigh,
    // --- Input/Output Buffers for Best State ---
    __global uchar*        occ_best, // Stores the best configuration found globally
    __global float*        E_best,   // Stores the best energy found globally
    __global float2*       Itot_out, // Final current output
    const int              max_neighs,
    int                    initMode
) {
    //=========================================================================
    // 1. Initialization and ID setup
    //=========================================================================
    const int itip       = get_group_id(0);
    const int lid        = get_local_id(0);
    const int local_size = get_local_size(0);

    if (itip >= nTips) return;
    if (nSite > MAX_SITES_BIG || local_size > MAX_WORKGROUP_SIZE_BIG) return;

    // Unpack solver parameters
    const float kT         = solver_params.x;
    const int   nIterTotal = (int)(solver_params.y + 0.5);
    const int   solverMode = (int)(solver_params.z + 0.5);
    uint rng_state         = (uint)solver_params.w + get_global_id(0);
    rng_state = wang_hash(rng_state);

    const int occ_bytes = min((nSite + 7) / 8, OCC_BYTES);
    const int tip_offset = itip * nSite;
    const int tip_occ_offset = itip * OCC_BYTES;

    // --- Shared memory for the work-group ---
    __local uchar occ_mask[OCC_BYTES];      // Current working configuration
    __local uchar occ_best_l[OCC_BYTES];    // Best configuration found in this run
    __local float reduction_dE[MAX_WORKGROUP_SIZE_BIG];
    __local int   reduction_site[MAX_WORKGROUP_SIZE_BIG];

    float E_best_l; // Best energy found in this run

    //=========================================================================
    // 2. Initialize State
    //=========================================================================
    if (lid == 0) {
        //initMode = 1; // DEBUG - override initMode
        if     (initMode == 0 ){ for (int i=0;i<occ_bytes;++i) { occ_mask[i] = 0;    } } 
        else if(initMode == 1 ){ for (int i=0;i<occ_bytes;++i) { occ_mask[i] = 0xFF; } } 
        else if(initMode == 2 ){ for (int i=0;i<occ_bytes;++i) { occ_mask[i] = occ_best[ itip*OCC_BYTES + i]; } } 
        else if(initMode == 3 ){ for (int i=0;i<occ_bytes;++i) { occ_mask[i] = wang_hash_uint(&rng_state)&0xFF; } }
        //if(itip==iDBG){ printf("GPU itip %3i occ: ", itip ); for(int i=0;i<occ_bytes;++i){ printf("%02x", occ_mask[i]); } printf("\n"); }
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    // Perform initial relaxation steps to settle the configuration
    for (int iter = 0; iter < nRelaxSteps; ++iter) {
        int i_site = wang_hash_uint(&rng_state) % nSite;
        uchar n_i = GET_OCC(i_site, occ_mask);
        float Ei = Esite[tip_offset + i_site];
        if(max_neighs > 0){
            const int iw0 = i_site * max_neighs;
            for (int k = 0; k < nNeigh[i_site]; ++k) {
                Ei += (float)GET_OCC(W_idx[iw0 + k], occ_mask) * W_val[iw0 + k];
            }
        }
        float dE = (n_i > 0) ? -Ei : Ei;
        reduction_dE[lid] = dE;
        reduction_site[lid] = i_site;
        barrier(CLK_LOCAL_MEM_FENCE);
        if (lid == 0) {
            if (reduction_dE[0] < 0.0f) { // Simple greedy flip for relaxation
                FLIP_OCC(reduction_site[0], occ_mask);
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    // After relaxation, calculate the initial energy and set it as the best
    if (lid == 0) {
        E_best_l = calculate_total_energy_parallel(occ_mask, nSite, Esite, W_val, W_idx, nNeigh, max_neighs, tip_offset, reduction_dE);
        for (int i = 0; i < occ_bytes; ++i) { occ_best_l[i] = occ_mask[i]; }
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    //=========================================================================
    // 3. Main Hybrid Optimization Loop
    //=========================================================================
    const int nGlobalSteps = max(1, nIterTotal / nLocalIter);
    for (int iGlobal = 0; iGlobal < nGlobalSteps; ++iGlobal) {

        // --- A: Perform a batch of local updates ---
        for (int iLocal = 0; iLocal < nLocalIter; ++iLocal) {
            int i_site = (solverMode == 0) ? ((iLocal * local_size + lid) % nSite) : (wang_hash_uint(&rng_state) % nSite);
            uchar n_i = GET_OCC(i_site, occ_mask);
            float Ei = Esite[tip_offset + i_site];
            if(max_neighs > 0){
                const int iw0 = i_site * max_neighs;
                for (int k = 0; k < nNeigh[i_site]; ++k) {
                    Ei += (float)GET_OCC(W_idx[iw0 + k], occ_mask) * W_val[iw0 + k];
                }
            }
            reduction_dE[lid] = (n_i > 0) ? -Ei : Ei;
            reduction_site[lid] = i_site;
            barrier(CLK_LOCAL_MEM_FENCE);

            if (lid == 0) {
                float dEmin = 1e10f;
                int imin = -1;
                for (int il = 0; il < local_size; ++il) {
                    if (reduction_dE[il] < dEmin) {
                        dEmin = reduction_dE[il];
                        imin = reduction_site[il];
                    }
                }
                if (imin >= 0) {
                    bool bDo = (dEmin < 0.0f);
                    if (!bDo && solverMode == 2 && kT > 1e-9) {
                        if (wang_hash_float(&rng_state) < exp(-dEmin / kT)) { bDo = true; }
                    }
                    if (bDo) { FLIP_OCC(imin, occ_mask); }
                }
            }
            barrier(CLK_LOCAL_MEM_FENCE);
        } // End of local iterations

        // --- B: Evaluate, Compare, and Prepare for Next Global Step ---
        float E_current = calculate_total_energy_parallel(occ_mask, nSite, Esite, W_val, W_idx, nNeigh, max_neighs, tip_offset, reduction_dE);

        if (lid == 0) {
            if (E_current < E_best_l) {
                E_best_l = E_current;
                for (int i = 0; i < occ_bytes; ++i) { occ_best_l[i] = occ_mask[i]; }
            }
            
            // Decide what configuration to start the next local search from
            if (wang_hash_float(&rng_state) < randConfigProb) {
                // Reset to a completely random state
                for (int i = 0; i < occ_bytes; ++i) { occ_mask[i] = wang_hash_uint(&rng_state) & 0xFF; }
            } else {
                // Start from the best known state
                for (int i = 0; i < occ_bytes; ++i) { occ_mask[i] = occ_best_l[i]; }
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE); // Ensure all threads see the new occ_mask
    }

    //=========================================================================
    // 4. Finalization and Output
    //=========================================================================
    if (lid == 0) {
        // Only write to global memory if we found a better solution than what's already there
        if (E_best_l < E_best[itip]) {
            E_best[itip] = E_best_l;
            for (int j = 0; j < occ_bytes; ++j) {
                occ_best[tip_occ_offset + j] = occ_best_l[j];
            }

            float Iocc = 0.0f, Iunoc = 0.0f;
            for (int i = 0; i < nSite; ++i) {
                if (GET_OCC(i, occ_best_l)) {
                    Iocc += Tsite[tip_offset + i];
                } else {
                    Iunoc += Tsite[tip_offset + i];
                }
            }
            Itot_out[itip] = (float2)(Iocc, Iunoc);
        }
    }
}



/**
 * @file precalc_esite_thop.cl
 * @brief OpenCL kernel to pre-calculate on-site energies and hopping amplitudes.
 *
 * This kernel implements a sophisticated physical model including:
 * - Multipole interactions (up to quadrupole order) between the tip and sample sites.
 * - An electrostatic mirror charge effect from a conductive substrate.
 * - An optional linear potential ramp to model a flat capacitor field.
 * - Optional rotation of the sample site's local coordinate frame.
 *
 * It is designed to generate the input `Esite_Thop` array required by the
 * local update Monte Carlo solvers.
 *
 * Parallelism:
 * - Each work-item computes the energy and hopping for one unique (tip, site) pair.
 * - Global work size should be (nTips * nSingle).
 */


// --- Type Definitions ---
typedef struct __attribute__ ((packed)) { float4 a; float4 b; float4 c; } cl_Mat3;

// --- Helper Functions ---
float3 mat3_dot_vec3(cl_Mat3 m, float3 v) {  return (float3)( dot(m.a.xyz, v), dot(m.b.xyz, v), dot(m.c.xyz, v) ); }

/**
 * @brief Computes multipole interaction energy.
 * MODIFIED: Now takes a private pointer to coefficients.
 */
float Emultipole(__private const float* cs, float3 d, int nMulti) {
    float ir2 = 1.0f / dot(d, d);    
    float E = cs[0];  // Assumes cs[0] always exists
    if (nMulti >= 1) {   // Dipole term (requires at least 4 coefficients: c0, px, py, pz)
        E += ir2 * ( cs[1] * d.x + 
                     cs[2] * d.y + 
                     cs[3] * d.z );  }  
    if (nMulti >= 4) {   // diagonal Quadrupole term
        E += ir2 * ir2 * ( cs[4] * d.x*d.x +
                           cs[5] * d.y*d.y +
                           cs[6] * d.z*d.z ); }  
    if (nMulti >= 7) {  // off-diagonal Quadrupole term
        E += ir2 * ir2 * (  cs[7] * d.y*d.z +
                            cs[8] * d.z*d.x +
                            cs[9] * d.x*d.y ); }
    return sqrt(ir2) * E;
}

// --- Main Kernel ---
__kernel void precalc_Esite_Thop(
    // --- System & Tip Configuration
    const int nSingle,
    const int nTips,
    __global const float4*  posE,
    __global const cl_Mat3* rotSite,
    __global const float4*  pTips,

    // --- Physical Model Parameters
    const float  Rtip,
    const float2 zV,
    // MODIFIED: Multipole coefficients are now a site-dependent global array
    __global const float* multipoleCoefs,
    const int    nMulti,             // Max number of multipole coeffs per site
    const int    bMirror,
    const int    bRamp,
    const float  Thop_decay,
    const float  Thop_amp,

    // --- Output Array
    __global float* Esite,
    __global float* Tsite
    
) {
    // 1. ID setup
    const int gid = get_global_id(0);
    if (gid >= (nTips * nSingle)) return;
    const int iTip  = gid / nSingle;
    const int iSite = gid % nSingle;

    // 2. Load Data
    const float4 tip_data  = pTips[iTip];
    const float4 site_data = posE[iSite];

    if( (iTip==0)&&(iSite==0) ) {
        printf("GPU: Tip( %12.8f , %12.8f , %12.8f | %12.8f ) Site( %12.8f , %12.8f , %12.8f | %12.8f ) \n", tip_data.x, tip_data.y, tip_data.z, tip_data.w, site_data.x, site_data.y, site_data.z, site_data.w );
        // print site rotation
        if(rotSite) {
            cl_Mat3 rot = rotSite[iSite];
            printf("GPU: Site %i rotation: \n", iSite);
            printf("GPU: %12.8f %12.8f %12.8f %12.8f \n", rot.a.x, rot.a.y, rot.a.z, rot.a.w);
            printf("GPU: %12.8f %12.8f %12.8f %12.8f \n", rot.b.x, rot.b.y, rot.b.z, rot.b.w);
            printf("GPU: %12.8f %12.8f %12.8f %12.8f \n", rot.c.x, rot.c.y, rot.c.z, rot.c.w);
        } 

    }

    float3       pTip      = tip_data.xyz;
    const float  VBias     = tip_data.w;
    const float3 pSite     = site_data.xyz;
    const float  E0        = site_data.w;

    // 3. On-Site Energy (Esite) Calculation
    const float zV0   = zV.x;
    const float zV1   = pTip.z + zV.y;
    float3 pTipMirror = pTip;
    pTipMirror.z      = 2.0f * zV0 - pTip.z;
    pTip             -= pSite;
    pTipMirror       -= pSite;

    // if (rotSite) {
    //     cl_Mat3 rot = rotSite[iSite];
    //     pTip       = mat3_dot_vec3(rot, pTip);
    //     pTipMirror = mat3_dot_vec3(rot, pTipMirror);
    // }

    // MODIFIED: Copy this site's coefficients from global to private memory. This makes access within Emultipole much faster.
    __private float cs[16]; // Max supported multipole coeffs
    const int nMulti_ = min(nMulti, 16);
    for(int i = 0; i < nMulti_; ++i){  cs[i] = multipoleCoefs[iSite * nMulti + i];  }


    
    float E                 = Emultipole(cs, pTip, nMulti);
    //if (bMirror != 0) {  E -= Emultipole(cs, pTipMirror, nMulti); }
    E *= (VBias * Rtip);

    if( (iTip==0)&&(iSite==0) ) {
        printf("GPU: VBias %12.8f    Rtip %12.8f  zV1 %12.8f  zV0 %12.8f \n", VBias, Rtip, zV1, zV0 );
        printf("GPU: E %12.8f pTip %12.8f %12.8f %12.8f    cs:   ", E, pTip.x, pTip.y, pTip.z);
        for(int i = 0; i < nMulti_; ++i){  printf("%12.8f ", cs[i]);  }
        printf("\n");
    }

    if (bRamp != 0 && (zV1 - zV0) != 0.0f) {
        float ramp = clamp((pSite.z - zV0) / (zV1 - zV0), 0.0f, 1.0f);
        // Ramp term only interacts with the monopole component (cs[0])
        if(nMulti > 0) E += (cs[0] * VBias * ramp);
    }
    const float final_Esite = E + E0;

    // 4. Hopping Amplitude (Thop) Calculation
    const float r          = distance(tip_data.xyz, site_data.xyz);
    const float final_Thop = Thop_amp * exp(-Thop_decay * r);

    // 5. Store Results
    Esite[gid] = final_Esite;
    Tsite[gid] = final_Thop;
}