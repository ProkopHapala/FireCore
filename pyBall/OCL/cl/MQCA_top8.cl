/**
 * @file MQCA_top8.cl
 * @brief OpenCL kernel for MQCA ground-state solver tracking top-8 states.
 *
 * This kernel tracks the 8 lowest energy states to detect degeneracy.
 */

#define MAX_SITES 16
#define MAX_NEIGH 8
#define TOP_K 8

// ============================================================
//  Helper: count trailing zeros
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
//  Kernel: mqca_groundstate_top8
//
//  Tracks TOP 8 lowest energy states for detecting degeneracy.
//  Uses simple bubble-sort at end instead of insertion during loop.
//
//  Outputs:
//    E_top8    [nInstances * 8]  : energies (sorted low to high)
//    occ_top8  [nInstances * 8]  : occupancy masks
// ============================================================
__kernel void mqca_groundstate_top8(
    const int nSite,
    const int nInstances,
    __global const float* Esite,
    __global const float* W_val,
    __global const int*   W_idx,
    __global const int*   nNeigh,
    __global float*       E_top8,
    __global int*         occ_top8
) {
    const int inst = get_group_id(0);
    const int lid  = get_local_id(0);

    if (inst >= nInstances) return;

    __local float  dE_contrib[MAX_SITES];
    __local float  E_current_l;
    __local uint   occ_mask_l;

    // Local arrays for top-K tracking
    __local float  E_best[TOP_K];
    __local uint   occ_best[TOP_K];

    const int inst_site_off  = inst * nSite;
    const int inst_neigh_off = inst * nSite * MAX_NEIGH;
    const int out_off        = inst * TOP_K;

    if (lid == 0) {
        occ_mask_l  = 0u;
        E_current_l = 0.0f;
        // Initialize with large values
        for (int k = 0; k < TOP_K; k++) {
            E_best[k] = 1e20f;
            occ_best[k] = 0u;
        }
        // State 0 (all zeros) is first candidate
        E_best[0] = 0.0f;
        occ_best[0] = 0u;
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

            // Simple replacement: if better than worst, replace worst
            if (E_current_l < E_best[TOP_K-1]) {
                E_best[TOP_K-1] = E_current_l;
                occ_best[TOP_K-1] = occ_mask_l;
                // Simple bubble-up to maintain partial order
                for (int k = TOP_K-1; k > 0; k--) {
                    if (E_best[k] < E_best[k-1]) {
                        float tmpE = E_best[k];
                        uint tmpO = occ_best[k];
                        E_best[k] = E_best[k-1];
                        occ_best[k] = occ_best[k-1];
                        E_best[k-1] = tmpE;
                        occ_best[k-1] = tmpO;
                    } else break;
                }
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    if (lid == 0) {
        // Final bubble sort to ensure correct order
        for (int i = 0; i < TOP_K; i++) {
            for (int j = i + 1; j < TOP_K; j++) {
                if (E_best[i] > E_best[j]) {
                    float tmpE = E_best[i];
                    uint tmpO = occ_best[i];
                    E_best[i] = E_best[j];
                    occ_best[i] = occ_best[j];
                    E_best[j] = tmpE;
                    occ_best[j] = tmpO;
                }
            }
        }

        for (int k = 0; k < TOP_K; k++) {
            E_top8[out_off + k]   = E_best[k];
            occ_top8[out_off + k] = (int)occ_best[k];
        }
    }
}
