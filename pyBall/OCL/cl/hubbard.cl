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

#define MAX_SITES_BIG 256
#define MAX_NEIGHBORS 16
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

float wang_hash_float(uint* seed) {
    *seed = wang_hash(*seed);
    return (float)(*seed) / (float)(UINT_MAX);
}

__kernel void solve_local_updates_parallel(
    const int nSite, 
    const int nTips, 
    const float4 solver_params,
    __global const float2* Esite_Thop,
    __global const float*  W_val, 
    __global const int*    W_idx, 
    __global const int*    nNeigh,
    __global ulong*        state_out, 
    __global float*        E_out, 
    __global float2*       Itot_out
) {
    //=========================================================================
    // 1. Initialization and ID setup
    //=========================================================================
    const int tip_idx    = get_group_id(0);
    const int lid        = get_local_id(0);
    const int local_size = get_local_size(0);

    if (tip_idx >= nTips) return;
    if (nSite > 64 || local_size > MAX_WORKGROUP_SIZE) return;

    // Unpack solver parameters
    const float kT    = solver_params.x;
    const int   nIter = (int)solver_params.y;
    const int   mode  = (int)solver_params.z;
    uint rng_state = (uint)solver_params.w + tip_idx * local_size + lid; // Unique seed per thread

    // --- Shared memory for the work-group ---
    // State is stored as a bitmask to save memory
    __local uchar occ_mask[MAX_SITES / 8 + 1];
    // Reduction arrays to find the best move
    __local float reduction_dE  [MAX_WORKGROUP_SIZE];
    __local int   reduction_site[MAX_WORKGROUP_SIZE];

    // Thread 0 initializes the state to all zeros
    if (lid == 0) {
        for (int i = 0; i < (nSite / 8 + 1); ++i) {
            occ_mask[i] = 0;
        }
    }
    barrier(CLK_LOCAL_MEM_FENCE); // Ensure all threads see the initialized state

    // Offset for accessing tip-specific global data
    const int tip_offset = tip_idx * nSite;

    //=========================================================================
    // 2. Main Collaborative Monte Carlo Loop
    //=========================================================================
    for (int iter = 0; iter < nIter; ++iter) {
        // --- Step A: Each thread proposes a move in parallel ---
        int i_site;
        if (mode == 0) { // Deterministic sweep across threads
            i_site = (iter * local_size + lid) % nSite;
        } else { // Stochastic proposal
            i_site = wang_hash(&rng_state) % nSite;
        }

        const uchar n_i = GET_OCC(i_site, occ_mask);

        float h_local          = Esite_Thop[tip_offset + i_site].x;
        const int w_row_offset = i_site * MAX_NEIGHBORS;
        for (int k = 0; k < nNeigh[i_site]; ++k) {
            h_local += (float)GET_OCC(W_idx[w_row_offset + k], occ_mask) * W_val[w_row_offset + k];
        }

        reduction_dE[lid] = (1.0f - 2.0f * (float)n_i) * h_local;
        reduction_site[lid] = i_site;

        barrier(CLK_LOCAL_MEM_FENCE); // Wait for all proposals to be calculated

        // --- Step B: Thread 0 finds the best move and decides ---
        if (lid == 0) {
            float best_dE = 1e10f; // Start with a high energy
            int   best_site_to_flip = -1;

            // Find the move with the minimum dE among all threads
            for (int i = 0; i < local_size; ++i) {
                if (reduction_dE[i] < best_dE) {
                    best_dE           = reduction_dE[i];
                    best_site_to_flip = reduction_site[i];
                }
            }

            bool accept_move = false;
            if (best_dE < 0.0f) {
                accept_move = true;
            } else if (mode == 2 && kT > 1e-9f) { // Simulated Annealing
                if (wang_hash_float(&rng_state) < exp(-best_dE / kT)) {
                    accept_move = true;
                }
            }
            
            // If accepted, flip the bit in the shared state mask
            if (accept_move && best_site_to_flip != -1) {
                FLIP_OCC(best_site_to_flip, occ_mask);
            }
        }
        
        barrier(CLK_LOCAL_MEM_FENCE); // Wait for state update before next iteration
    }

    //=========================================================================
    // 3. Post-Loop: Thread 0 calculates final properties and writes output
    //=========================================================================
    if (lid == 0) {
        float final_E = 0.0f;
        float final_Iocc = 0.0f;
        float final_Iunoc = 0.0f;
        ulong final_state_mask = 0;

        for (int i = 0; i < nSite; ++i) {
            if (GET_OCC(i, occ_mask)) {
                final_state_mask |= (1UL << i);
                final_E += Esite_Thop[tip_offset + i].x;
                final_Iocc += Esite_Thop[tip_offset + i].y;

                const int w_row_offset = i * MAX_NEIGHBORS;
                for (int k = 0; k < nNeigh[i]; ++k) {
                    const int j = W_idx[w_row_offset + k];
                    if (i < j && GET_OCC(j, occ_mask)) {
                        final_E += W_val[w_row_offset + k];
                    }
                }
            } else {
                final_Iunoc += Esite_Thop[tip_offset + i].y;
            }
        }

        state_out[tip_idx] = final_state_mask;
        E_out[tip_idx]     = final_E;
        Itot_out[tip_idx]  = (float2)(final_Iocc, final_Iunoc);
    }
}