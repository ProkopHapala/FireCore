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

__kernel void solve_minBrute(
    // --- Input: System parameters
    const int nSingle,         // Number of sites (must be <= MAX_SITES)
    __global const float4* posE,      // [nSingle] {x, y, z, E0} - Site positions and on-site energies
    // --- Input: Tip parameters (defines the parallel work)
    const int nTips,           // Number of tip positions (total work items)
    __global const float4* pTips,     // [nTips] {x, y, z, Vbias} - Tip positions and bias voltages
    const float tipDecay,      // Decay constant 'b' for tunneling
    const float tipRadius,     // Effective radius of the tip for interaction potential
    // --- Output: Results for each tip position
    __global float*  Emin,      // [nTips] Energy of the many-body ground state
    __global int*    iMin,      // [nTips] Index of the ground state (bitmask of occupancies)
    __global float2* Itot       // [nTips] {I_occ, I_unoc} - Tunneling currents
) {
    //=========================================================================
    // 1. Initialization and ID setup
    //=========================================================================
    const int gid = get_global_id(0);
    if (gid >= nTips){ return; }
    if (nSingle > MAX_SITES) {Emin[gid] = 1e10f; return;} // Ensure nSingle does not exceed our static allocation limit

    // Local memory for the interaction matrix W_ij. Shared by all threads in a work-group.
    __local float W[MAX_WIJ];

    //=========================================================================
    // 2. Collaborative Pre-calculation of W_ij Matrix (using Local Memory)
    //=========================================================================
    // This calculation is independent of the tip position, so all threads in a work-group
    // can collaborate to compute it once. This avoids redundant calculations.
    const int lid        = get_local_id(0);
    const int local_size = get_local_size(0);
    for (int i=lid; i<nSingle; i+=local_size) {
        float3 pi = posE[i].xyz;
        int iw0 = i*(nSingle-i-1)/2;
        for (int j=i+1; j<nSingle; ++j) {
            float3 pj = posE[j].xyz;
            int iw = iw0 + j-i-1;
            W[iw] = get_Wij(pi, pj);
        }
    }
    barrier(CLK_LOCAL_MEM_FENCE); // Synchronize all threads in the work-group to ensure W is fully computed before anyone uses it.

    //=========================================================================
    // 3. Per-Thread Pre-calculation of Tip-Dependent Terms
    //=========================================================================
    
    // Esites must be private ( because they depend on tip position cannot be shared )
    float Esite[MAX_SITES]; // U_i = E0_i + V_i (On-site + tip potential)
    { // memory sanity
    const float4 tip = pTips[gid];
    for (int i = 0; i < nSingle; ++i) {
        Esite[i]      = posE[i].w + Vtip(posE[i].xyz, tip.xyz, tip.w, tipRadius);
    }
    }
    
    //=========================================================================
    // 4. Brute-Force Minimization Loop
    //=========================================================================
    // WARNING: This loop runs 2^nSingle times. It will be very slow for nSingle > 22.
    const uint nMany = 1 << nSingle;
    float min_E   = 1e10f;
    uint  min_idx = 0;
    for (uint iMany=0; iMany<nMany; ++iMany) {
        float E = 0.0f;        
        // Calculate energy for the current state configuration (state_idx)
        int   iW = 0;
        for (int i = 0; i<nSingle; ++i){  
            bool ni = (iMany >> i) & 1; // Check if site 'i' is occupied
            if (ni){ E += Esite[i]; }   // Add on-site and tip potential energy                   
            for (int j=i+1; j<nSingle; ++j) {  //  sum_{i<j} n_i*n_j*W_ij
                if ( ni && ((iMany>>j)&1) ) { E += W[iW]; }
                iW++; 
            }
        }
        if (E < min_E) { min_E = E; min_idx = iMany; } // Check if this state has a lower energy
    }
    Emin[gid] = min_E;
    iMin[gid] = (int)min_idx;

    //=========================================================================
    // 5. Post-Minimization: Calculate Tunneling Current for the Ground State
    //=========================================================================
    float I_occ  = 0.0f;
    float I_unoc = 0.0f;
    for (int i=0; i<nSingle; ++i) {  // Check occupancy in the found ground state (min_idx)
        float dist_tip_site = distance(tip.xyz, posE[i].xyz);
        float T = exp(-tipDecay * dist_tip_site);
        if ((min_idx >> i) & 1){ I_occ += T; }else{ I_unoc += T; }
    }
    
    //=========================================================================
    // 6. Store Final Results in Global Memory
    //=========================================================================

    Itot[gid] = (float2)(I_occ, I_unoc);
}



/**
 * @file hubbard_solver_workgroup_per_tip.cl
 * @brief Alternative kernel using one work-group per tip position for higher parallelism.
 *
 * This kernel assigns an entire work-group to a single tip position. This allows:
 * 1. Storing tip-dependent Esite and T arrays in fast __local memory.
 * 2. Parallelizing the brute-force search over all threads in the work-group.
 *
 * Parallelization Strategy:
 * - tip_idx      = get_group_id(0)
 * - lid          = get_local_id(0)
 * - Global Size  = nTips * workgroup_size
 * - Local Size   = workgroup_size (e.g., 64, 128, 256)
 *
 * All threads in a work-group (e.g., 128 threads) collaborate to find the
 * ground state for a single tip position before moving to the next.
 */

__kernel void solve_minBrute_loc(
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
    if (nSingle > MAX_SITES) {if (lid == 0) Emin[tip_idx] = 1e10f; return; }

    // Shared memory for all threads in this work-group (processing one tip)
    
    float local_min_E = 1e10f;
    uint  local_min_idx = 0;

    //{ // memory sanity  float[MAX_SITES+MAX_WIJ]
    __local float Esite[MAX_SITES];
    __local float W    [MAX_WIJ];

    //=========================================================================
    // 2. Collaborative Pre-calculation in Local Memory
    //=========================================================================
    // All threads in the workgroup collaborate to fill the shared arrays.
    
    // a) Calculate W_ij matrix (tip-independent)
    for (int i = lid; i < nSingle; i += local_size) {
        for (int j = i + 1; j < nSingle; ++j) {
            // This is the correct, robust formula for the flattened index k
            int k = (i * nSingle - i * (i + 1) / 2) + (j - i - 1);
            W[k] = get_Wij(posE[i].xyz, posE[j].xyz);
        }
    }

    // b) Calculate Esite and T arrays (tip-dependent)
    const float4 tip = pTips[tip_idx];
    for (int i = lid; i < nSingle; i += local_size) {
        Esite[i] = posE[i].w + Vtip(posE[i].xyz, tip.xyz, tip.w, tipRadius);
    }

    // Synchronize all threads to ensure W, Esite, and T are fully computed before use.
    barrier(CLK_LOCAL_MEM_FENCE);

    //=========================================================================
    // 3. Parallel Brute-Force Minimization
    //=========================================================================
    const uint nMany = 1 << nSingle;
    
    // Each thread finds the minimum in its own subset of states


    // Each thread (lid) processes states lid, lid+local_size, lid+2*local_size, ...
    for (uint iMany = lid; iMany < nMany; iMany += local_size) {
        float E = 0.0f;
        bool ni = (iMany >> i) & 1; // Check if site 'i' is occupied
        for (int i = 0; i < nSingle; ++i){ if (ni){ E += Esite[i]; } }
        // Interaction energy part
        int iW = 0;
        for (int i = 0; i < nSingle; ++i) {
            for (int j = i + 1; j < nSingle; ++j) {
                // Check if BOTH sites i and j are occupied
                if (ni && ((iMany >> j) & 1)) { E += W[iW]; }
                iW++; // Always advance index
            }
        }
        if (E < local_min_E) {
            local_min_E = E;
            local_min_idx = iMany;
        }
    }

    //} // end // memory sanity  float[MAX_SITES+MAX_WIJ]


    //{ // local memory sanity 

    //=========================================================================
    // 4. Work-group Reduction to find the final minimum
    //=========================================================================
    // We need to find the minimum among all threads in the work-group.
    __local float reduction_E  [MAX_WORKGROUP_SIZE];
    __local uint  reduction_idx[MAX_WORKGROUP_SIZE];

    reduction_E  [lid] = local_min_E;
    reduction_idx[lid] = local_min_idx;
    barrier(CLK_LOCAL_MEM_FENCE);

    // Thread 0 of the work-group performs the final reduction and writes the result.
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
            float r = distance(tip.xyz, posE[i].xyz);
            float T = exp(-tipDecay * r);
            if ((final_min_idx >> i) & 1){ I_occ += T; }else{ I_unoc += T; }
        }
        Itot[tip_idx] = (float2)(I_occ, I_unoc);
    }
    
    //}
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



