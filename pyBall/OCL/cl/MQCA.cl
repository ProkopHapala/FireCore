/**
 * @file MQCA.cl
 * @brief OpenCL kernels for Molecular Quantum Cellular Automata (MQCA) ground-state solver.
 *
 * Hamiltonian:
 *   H = sum_i n_i [ eps_i + sum_{j<i} W_ij n_j ]
 *
 * Key design:
 *  - Max 16 active sites => 2^16 = 65536 states, brute-force is tractable
 *  - Gray code traversal: only 1 bit flips per step => incremental energy update
 *  - One workgroup per system instance, 16 threads per workgroup (one per site)
 *  - Sparse W_ij: max MAX_NEIGH non-zero neighbors per site
 *  - Fixed (boundary) sites pre-folded into eps_i bias on the Python side
 *
 * Gray code:  G(i) = i ^ (i>>1)
 * Flip index: bit that changed = ctz( G(i) ^ G(i-1) )
 */

#define MAX_SITES 16
#define MAX_NEIGH 8

// ============================================================
//  Helper: count trailing zeros (which bit flipped)
// ============================================================
inline int ctz16(uint x) {
    int n = 0;
    if ((x & 0x00FF) == 0) { n += 8; x >>= 8; }
    if ((x & 0x000F) == 0) { n += 4; x >>= 4; }
    if ((x & 0x0003) == 0) { n += 2; x >>= 2; }
    if ((x & 0x0001) == 0) { n += 1; }
    return n;
}

// ============================================================
//  Compute total energy for a given occupancy mask
//  (used for initial E_current = 0 check at occ=0)
// ============================================================
inline float total_energy_mask(
    uint occ_mask, int nSite,
    __global const float* Esite,      // [nSite]
    __global const float* W_val,      // [nSite * MAX_NEIGH]
    __global const int*   W_idx,      // [nSite * MAX_NEIGH]
    __global const int*   nNeigh      // [nSite]
) {
    float E = 0.0f;
    for (int i = 0; i < nSite; i++) {
        if (!((occ_mask >> i) & 1u)) continue;
        E += Esite[i];
        int iw0 = i * MAX_NEIGH;
        for (int k = 0; k < nNeigh[i]; k++) {
            int j = W_idx[iw0 + k];
            if (j < i && ((occ_mask >> j) & 1u)) {  // only j<i to count each pair once
                E += W_val[iw0 + k];
            }
        }
    }
    return E;
}

// ============================================================
//  Kernel: mqca_groundstate
//
//  Finds the ground-state energy and occupancy for one or more
//  independent MQCA instances using Gray code brute-force.
//
//  Parallelisation:
//    get_group_id(0)  = instance index
//    get_local_id(0)  = site index (0..nSite-1, threads nSite..15 are idle)
//
//  Global size  = nInstances * 16
//  Local  size  = 16
//
//  Inputs per instance (instance-major layout, stride = nSite):
//    Esite[inst*nSite + i]  : on-site energy of site i (including fixed-neighbor bias)
//
//  Shared coupling (same for all instances in a batch that use same geometry):
//    W_val[i*MAX_NEIGH + k] : k-th non-zero coupling of site i
//    W_idx[i*MAX_NEIGH + k] : column index j for that coupling
//    nNeigh[i]              : number of neighbors of site i
//
//  Outputs per instance:
//    E_out[inst]   : ground-state energy
//    occ_out[inst] : ground-state 16-bit occupancy mask
// ============================================================
__kernel void mqca_groundstate(
    const int nSite,                       // 1
    const int nInstances,                  // 2
    __global const float* Esite,           // 3  [nInstances * nSite]
    __global const float* W_val,           // 4  [nSite * MAX_NEIGH]
    __global const int*   W_idx,           // 5  [nSite * MAX_NEIGH]
    __global const int*   nNeigh,          // 6  [nSite]
    __global float*       E_out,           // 7  [nInstances]
    __global int*         occ_out          // 8  [nInstances]  (ushort fits in int)
) {
    const int inst = get_group_id(0);
    const int lid  = get_local_id(0);   // 0..15

    if (inst >= nInstances) return;

    // --- local memory shared within workgroup ---
    __local float  dE_contrib[MAX_SITES];   // each thread writes its contribution
    __local float  E_current_l;
    __local uint   occ_mask_l;
    __local float  E_min_l;
    __local uint   occ_min_l;

    const int inst_off = inst * nSite;

    // --- Thread 0 initialises shared state ---
    if (lid == 0) {
        occ_mask_l = 0u;
        E_min_l    = 0.0f;   // E(all-zero) = 0 by definition
        occ_min_l  = 0u;
        E_current_l = 0.0f;
    }
    if (lid < MAX_SITES) dE_contrib[lid] = 0.0f;
    barrier(CLK_LOCAL_MEM_FENCE);

    // ---- Gray code loop over all 2^nSite states ----
    const int nStates = 1 << nSite;
    uint prev_g = 0u;

    for (int step = 1; step < nStates; step++) {
        uint g    = (uint)step ^ ((uint)step >> 1);   // current Gray code
        uint diff = g ^ prev_g;                        // exactly one bit set
        int  flip = ctz16(diff);                       // which site flips (0-based)
        prev_g = g;

        // ---- Each thread computes its site's coupling contribution ----
        float contrib = 0.0f;
        if (lid < nSite) {
            if (lid == flip) {
                // dE for flipping site `flip`:
                //   dE = sign * (eps_flip + sum_{k: neighbor k occupied} W_{flip,k})
                float site_E = Esite[inst_off + flip];
                uint  cur_occ = occ_mask_l;  // read shared (set by thread 0)
                int   iw0 = flip * MAX_NEIGH;
                for (int k = 0; k < nNeigh[flip]; k++) {
                    int   j  = W_idx[iw0 + k];
                    float wv = W_val[iw0 + k];
                    if ((cur_occ >> j) & 1u) site_E += wv;
                }
                // sign: if currently occupied (+1 -> will become 0, so energy goes down)
                int n_flip = (cur_occ >> flip) & 1u;
                contrib = (n_flip ? -1.0f : 1.0f) * site_E;
            }
        }
        dE_contrib[lid] = contrib;
        barrier(CLK_LOCAL_MEM_FENCE);

        // ---- Thread 0 aggregates and tracks minimum ----
        if (lid == 0) {
            float dE = 0.0f;
            for (int s = 0; s < nSite; s++) dE += dE_contrib[s];
            E_current_l += dE;
            occ_mask_l ^= (1u << flip);   // flip the bit
            if (E_current_l < E_min_l) {
                E_min_l   = E_current_l;
                occ_min_l = occ_mask_l;
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    // ---- Write results ----
    if (lid == 0) {
        E_out  [inst] = E_min_l;
        occ_out[inst] = (int)occ_min_l;
    }
}

// ============================================================
//  Kernel: mqca_groundstate_batch_W
//
//  Variant where BOTH Esite AND W_val/W_idx differ per instance.
//  Useful for the geometry-variation genetic algorithm where each
//  individual has its own coupling matrix.
//
//  Layout:
//    Esite [inst*nSite + i]
//    W_val [inst*nSite*MAX_NEIGH + i*MAX_NEIGH + k]
//    W_idx [inst*nSite*MAX_NEIGH + i*MAX_NEIGH + k]
//    nNeigh[inst*nSite + i]
// ============================================================
__kernel void mqca_groundstate_batch_W(
    const int nSite,
    const int nInstances,
    __global const float* Esite,     // [nInstances * nSite]
    __global const float* W_val,     // [nInstances * nSite * MAX_NEIGH]
    __global const int*   W_idx,     // [nInstances * nSite * MAX_NEIGH]
    __global const int*   nNeigh,    // [nInstances * nSite]
    __global float*       E_out,     // [nInstances]
    __global int*         occ_out    // [nInstances]
) {
    const int inst = get_group_id(0);
    const int lid  = get_local_id(0);

    if (inst >= nInstances) return;

    __local float  dE_contrib[MAX_SITES];
    __local float  E_current_l;
    __local uint   occ_mask_l;
    __local float  E_min_l;
    __local uint   occ_min_l;

    const int inst_site_off  = inst * nSite;
    const int inst_neigh_off = inst * nSite * MAX_NEIGH;

    if (lid == 0) {
        occ_mask_l  = 0u;
        E_min_l     = 0.0f;
        occ_min_l   = 0u;
        E_current_l = 0.0f;
    }
    if (lid < MAX_SITES) dE_contrib[lid] = 0.0f;
    barrier(CLK_LOCAL_MEM_FENCE);

    const int nStates = 1 << nSite;
    uint prev_g = 0u;

    for (int step = 1; step < nStates; step++) {
        uint g    = (uint)step ^ ((uint)step >> 1);
        uint diff = g ^ prev_g;
        int  flip = ctz16(diff);
        prev_g = g;

        float contrib = 0.0f;
        if (lid < nSite && lid == flip) {
            float site_E = Esite[inst_site_off + flip];
            uint  cur_occ = occ_mask_l;
            int   iw0 = inst_neigh_off + flip * MAX_NEIGH;
            int   nn  = nNeigh[inst_site_off + flip];
            for (int k = 0; k < nn; k++) {
                int   j  = W_idx[iw0 + k];
                float wv = W_val[iw0 + k];
                if ((cur_occ >> j) & 1u) site_E += wv;
            }
            int n_flip = (cur_occ >> flip) & 1u;
            contrib = (n_flip ? -1.0f : 1.0f) * site_E;
        }
        dE_contrib[lid] = contrib;
        barrier(CLK_LOCAL_MEM_FENCE);

        if (lid == 0) {
            float dE = 0.0f;
            for (int s = 0; s < nSite; s++) dE += dE_contrib[s];
            E_current_l += dE;
            occ_mask_l ^= (1u << flip);
            if (E_current_l < E_min_l) {
                E_min_l   = E_current_l;
                occ_min_l = occ_mask_l;
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    if (lid == 0) {
        E_out  [inst] = E_min_l;
        occ_out[inst] = (int)occ_min_l;
    }
}
