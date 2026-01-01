
// =========================================================================
// OMM.cl - Order-N Molecular Orbitals Kernels
// Split-Strategy Implementation (3-Phase Design)
// =========================================================================

#define BLOCK_SIZE 256

// -------------------------------------------------------------------------
// KERNEL 1: apply_operator (Generic Convolution: \tilde{\psi} = \hat{O} \psi)
// Used for both \hat{S} and \hat{H}.
// Supports computing result on an expanded support (orbSupportExt).
// -------------------------------------------------------------------------
__kernel void apply_operator(
    const int nOrbs,
    const __global int*   orbSupport,      // Original support [nOrbs, nMaxSupport]
    const __global float* orbCoefs,        // Coefficients [nOrbs, nMaxSupport]
    const __global int*   nSupport,        // Length of original support [nOrbs]
    const __global int*   orbSupportExt,   // Expanded support [nOrbs, nMaxSupportExt]
    const __global int*   nSupportExt,     // Length of expanded support [nOrbs]
    const __global int*   atomNeigh,       // [nAtoms, nMaxNeigh]
    const __global float* atomMatrix,      // H or S values [nAtoms, nMaxNeigh]
    __global float*       outTwiddle,      // [nOrbs, nMaxSupportExt] Resulting \tilde{C}
    const int nMaxSupport,
    const int nMaxSupportExt,
    const int nMaxNeigh
) {
    int orb_idx = get_group_id(0);
    int lid     = get_local_id(0);

    if (orb_idx >= nOrbs) return;
    int len_orig = nSupport[orb_idx];
    int len_ext  = nSupportExt[orb_idx];
    
    // Each thread calculates one element of the transformed vector on the EXPANDED support
    for (int a_ext = lid; a_ext < len_ext; a_ext += get_local_size(0)) {
        int atom_a = orbSupportExt[orb_idx * nMaxSupportExt + a_ext];
        float res_a = 0.0f;
        int neigh_start = atom_a * nMaxNeigh;

        // Optimization: Find which atoms of the original support are neighbors of atom_a
        for (int n = 0; n < nMaxNeigh; n++) {
            int atom_b = atomNeigh[neigh_start + n];
            if (atom_b == -1) break;

            float Mat_ab = atomMatrix[neigh_start + n];

            // Find coefficient of atom_b in the ORIGINAL orbital's support
            // This is the bottleneck: O(len_orig) search. 
            // Since len_orig is small (e.g. 4-64), this is okay for now.
            for (int k = 0; k < len_orig; k++) {
                if (orbSupport[orb_idx * nMaxSupport + k] == atom_b) {
                    res_a += Mat_ab * orbCoefs[orb_idx * nMaxSupport + k];
                    break;
                }
            }
        }
        outTwiddle[orb_idx * nMaxSupportExt + a_ext] = res_a;
    }
}

// -------------------------------------------------------------------------
// KERNEL 2: project_orbitals (Sparse Dot Product: S_ij = <\psi_i | \tilde{\psi}_j>)
// Computes overlap between original support i and expanded support j.
// -------------------------------------------------------------------------
__kernel void project_orbitals(
    const int nPairs,
    const __global int2*  orbPairs,        // Pairs (i, j)
    const __global int*   orbSupport,      // Original [nOrbs, nMaxSupport]
    const __global float* orbCoefs,        // Original C [nOrbs, nMaxSupport]
    const __global int*   nSupport,        // [nOrbs]
    const __global int*   orbSupportExt,   // Expanded [nOrbs, nMaxSupportExt]
    const __global float* orbTwiddle,      // Expanded \tilde{C} [nOrbs, nMaxSupportExt]
    const __global int*   nSupportExt,     // [nOrbs]
    __global float*       outScalars,      // S_ij [nPairs]
    const int nMaxSupport,
    const int nMaxSupportExt
) {
    int pairIdx = get_group_id(0);
    if (pairIdx >= nPairs) return;

    int i_idx = orbPairs[pairIdx].x;
    int j_idx = orbPairs[pairIdx].y;
    int lid = get_local_id(0);

    // Load Orbital J's EXPANDED SUPPORT and TWIDDLE value into SLM
    __local int   supp_ext_j[BLOCK_SIZE];
    __local float twid_ext_j[BLOCK_SIZE];
    __local int   len_ext_j;

    if (lid == 0) len_ext_j = nSupportExt[j_idx];
    
    // Collaborative Load (if len_ext_j > BLOCK_SIZE, needs loop, but assuming <= 256 for now)
    for (int k = lid; k < nMaxSupportExt; k += get_local_size(0)) {
        supp_ext_j[k] = orbSupportExt[j_idx * nMaxSupportExt + k];
        twid_ext_j[k] = orbTwiddle[j_idx * nMaxSupportExt + k];
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    // Perform Sparse Dot Product: Sum_a c_{i,a} * \tilde{c}_{j,a}
    float my_dot = 0.0f;
    int len_orig_i = nSupport[i_idx];

    for (int a_orig = lid; a_orig < len_orig_i; a_orig += get_local_size(0)) {
        int atom_a = orbSupport[i_idx * nMaxSupport + a_orig];
        float c_ia = orbCoefs[i_idx * nMaxSupport + a_orig];

        // Find atom_a in J's EXPANDED support (SLM)
        for (int k = 0; k < len_ext_j; k++) {
            if (supp_ext_j[k] == atom_a) {
                // S_ij += c_ia * \tilde{c}_ja
                float term = c_ia * twid_ext_j[k];
                my_dot += term; 
                break;
            }
        }
    }

    // Reduction
    __local float reduce_buf[BLOCK_SIZE];
    reduce_buf[lid] = my_dot;
    barrier(CLK_LOCAL_MEM_FENCE);

    for (int s = BLOCK_SIZE / 2; s > 0; s >>= 1) {
        if (lid < s) reduce_buf[lid] += reduce_buf[lid + s];
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    if (lid == 0) outScalars[pairIdx] = reduce_buf[0];
}

// -------------------------------------------------------------------------
// KERNEL 3: assemble_gradient (Sparse Mixing: g_i = \tilde{\psi}^H_i - \sum \lambda_{ij} \tilde{\psi}^S_j)
// Note: result is on ORIGINAL support.
// -------------------------------------------------------------------------
__kernel void assemble_gradient(
    const int nOrbs,
    const __global int*   pairRowPtr,      // CSR Index for pairs [nOrbs+1]
    const __global int2*  orbPairs,        // [nTotalPairs]
    const __global float* scalarLambda,    // \lambda_{ij} [nTotalPairs]
    const __global int*   orbSupport,      // Original [nOrbs, nMaxSupport]
    const __global int*   nSupport,        // [nOrbs]
    const __global int*   orbSupportExt,   // Expanded [nOrbs, nMaxSupportExt]
    const __global float* twiddleH,        // Expanded \tilde{C}^H [nOrbs, nMaxSupportExt]
    const __global float* twiddleS,        // Expanded \tilde{C}^S [nOrbs, nMaxSupportExt]
    const __global int*   nSupportExt,     // [nOrbs]
    __global float*       outGrads,        // Resulting Gradient [nOrbs, nMaxSupport]
    const int nMaxSupport,
    const int nMaxSupportExt
) {
    int i_idx = get_group_id(0); // One workgroup handles Orbital I
    int lid   = get_local_id(0);

    if (i_idx >= nOrbs) return;

    int len_orig_i = nSupport[i_idx];
    int len_ext_i  = nSupportExt[i_idx];
    
    // Gradient result is on ORIGINAL support.
    // Initialize with H|i> projection onto original support.
    // We need to find which index in expanded support corresponds to each index in original support.
    float my_grad = 0.0f;
    int atom_a = -1;
    
    if (lid < len_orig_i) {
        atom_a = orbSupport[i_idx * nMaxSupport + lid];
        // Find atom_a in expanded support of i to get \tilde{H}_{ia}
        for (int k = 0; k < len_ext_i; k++) {
            if (orbSupportExt[i_idx * nMaxSupportExt + k] == atom_a) {
                my_grad = twiddleH[i_idx * nMaxSupportExt + k];
                break;
            }
        }
    }

    // Loop over interacting orbitals 'j' to add Orthogonality constraints: -\lambda_{ij} * S|j>
    int start_pair = pairRowPtr[i_idx];
    int end_pair   = pairRowPtr[i_idx + 1];

    __local int   supp_ext_j[BLOCK_SIZE];
    __local float twid_ext_s_j[BLOCK_SIZE];
    __local int   len_ext_j;

    for (int p = start_pair; p < end_pair; p++) {
        int j_idx = orbPairs[p].y; 
        float lambda = scalarLambda[p];

        if (lid == 0) len_ext_j = nSupportExt[j_idx];
        barrier(CLK_LOCAL_MEM_FENCE);
        
        for (int k = lid; k < nMaxSupportExt; k += get_local_size(0)) {
            supp_ext_j[k]   = orbSupportExt[j_idx * nMaxSupportExt + k];
            twid_ext_s_j[k] = twiddleS[j_idx * nMaxSupportExt + k];
        }
        barrier(CLK_LOCAL_MEM_FENCE);

        if (lid < len_orig_i) {
            for (int k = 0; k < len_ext_j; k++) {
                if (supp_ext_j[k] == atom_a) {
                    my_grad -= lambda * twid_ext_s_j[k];
                    break;
                }
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    if (lid < len_orig_i) {
        outGrads[i_idx * nMaxSupport + lid] = my_grad;
    }
}
