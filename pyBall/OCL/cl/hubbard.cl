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
    float3 d   = pSite - pTip;
    float  r   = length(d);
    return  Vbias * (Rtip/r); // Tip potential V_i = Vbias * (R_tip / |r_i - r_tip|)
}


float get_Wij( float3 pi, float3 pj ){
    // --- TODO: PHYSICAL MODEL FOR W_ij ---
    // The interaction W_ij depends on your specific physical system (e.g., screening, dipole type).
    // Replace the simple 1/r model below with your actual formula.
    // For dipole-dipole interaction, it might be proportional to 1/r^3.
    // For a simple screened Coulomb interaction, it could be: prefactor * exp(-dist / lambda) / dist.
    float3  d = pi - pj;
    float   r = length(d);
    return  COULOMB_CONST/r;
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

__kernel void solve_minBrute_fly(
    const int nSingle,
    __global const float4* posE,
    const int nTips,
    __global const float4* pTips,
    const float tipDecay,
    const float tipRadius,
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

    // Shared memory for this work-group. Note the absence of the large 'W' matrix.
    __local float4 Lsites[MAX_SITES];
    //__local float4 Tsite [MAX_SITES];

    //=========================================================================
    // 2. Collaborative Pre-calculation in Local Memory
    //=========================================================================
    // All threads in the workgroup collaborate to fill the shared arrays.
    const float4 tip = pTips[tip_idx];
    for (int i = lid; i < nSingle; i += local_size) {
        Lsites[i].xyz = posE[i].xyz;
        Lsites[i].w   = posE[i].w + Vtip(Lsites[i].xyz, tip.xyz, tip.w, tipRadius);
    }
    
    // Synchronize to ensure all shared arrays are fully computed.
    barrier(CLK_LOCAL_MEM_FENCE);

    //=========================================================================
    // 3. Parallel Brute-Force Minimization with On-the-Fly W_ij
    //=========================================================================
    const uint nMany = 1 << nSingle;
    
    float local_min_E   = 1e10f;
    uint  local_min_idx = 0;

    for (uint iMany = lid; iMany < nMany; iMany += local_size) {
        float E = 0.0f;
        for (int i = 0; i < nSingle; ++i) {  // On-site and interaction energy calculation
            if ((iMany >> i) & 1) {          // If site 'i' is occupied...
                float4 site = Lsites[i];
                E         += site.w;               // Add its on-site energy
                float3 pi  = site.xyz;
                for (int j = i + 1; j < nSingle; ++j) { // ...calculate interaction with other occupied sites j > i
                    if ((iMany >> j) & 1) {  E += get_Wij(pi, Lsites[j].xyz); }
                }
            }
        }
        if (E < local_min_E) {
            local_min_E = E;
            local_min_idx = iMany;
        }
    }

    //=========================================================================
    // 4. Work-group Reduction to find the final minimum
    //=========================================================================
    __local float reduction_E  [MAX_WORKGROUP_SIZE];
    __local uint  reduction_idx[MAX_WORKGROUP_SIZE];

    reduction_E[lid]   = local_min_E;
    reduction_idx[lid] = local_min_idx;
    barrier(CLK_LOCAL_MEM_FENCE);

    // Thread 0 performs the final reduction and writes the result.
    if (lid == 0) {
        float final_min_E   = reduction_E[0];
        uint  final_min_idx = reduction_idx[0];
        for (int i = 1; i < local_size; ++i) {
            if (reduction_E[i] < final_min_E) {
                final_min_E   = reduction_E[i];
                final_min_idx = reduction_idx[i];
            }
        }
        //=====================================================================
        // 5. Final Calculations & Output (done only by thread 0)
        //=====================================================================
        Emin[tip_idx] = final_min_E;
        iMin[tip_idx] = (int)final_min_idx;

        float I_occ  = 0.0f;
        float I_unoc = 0.0f;
        for (int i=0; i<nSingle; ++i) {  // Check occupancy in the found ground state (min_idx)
            float r = distance(tip.xyz, Lsites[i].xyz);
            float T = exp(-tipDecay * r);
            if ((final_min_idx >> i) & 1){ I_occ += T; }else{ I_unoc += T; }
        }
        Itot[tip_idx] = (float2)(I_occ, I_unoc);
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
    const float kT, // Thermal energy (k_B * T) in eV
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
    float local_Emin  = 1e10f;
    uint  local_imin  = 0;
    float local_Z     = 0.0f;
    float local_Iocc  = 0.0f;
    float local_Iunoc = 0.0f;
    const float inv_kT = (kT > 1e-9f) ? 1.0f / kT : 1e9f;

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
                for (int j = i + 1; j < nSingle; ++j) {
                    if ((iMany >> j) & 1) { E += get_Wij(site.xyz, Lsites[j].xyz); }
                }
            } else { i_unoc_state += T; }
        }

        // boltzmann-weighted average
        float w = exp((local_Emin - E) * inv_kT);
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

