/// XRDDebye.cl — Pair-distance histogram + Debye transform for powder diffraction
///
/// Kernel 1: pair_histogram — deposits interatomic distances into binned histograms
/// Kernel 2: pair_histogram_gaussian — same with Gaussian thermal broadening
/// Kernel 3: debye_transform_q — converts histograms to I(Q) via sinc(Qr)
///
/// Float atomics implemented with cmpxchg loop (portable across NVIDIA/AMD/Intel).

/// Portable float atomic add using 32-bit compare-and-swap.
inline void atomic_add_float(__global float* addr, float val) {
    union { uint u; float f; } old_val, new_val;
    do {
        old_val.f = *addr;
        new_val.f = old_val.f + val;
    } while (atom_cmpxchg((__global uint*)addr, old_val.u, new_val.u) != old_val.u);
}

/// Helper: linearly interpolate pair count into two adjacent bins
inline void deposit_pair(__global float* hist, int n_bins, int pair_type, float bin_f, float weight) {
    int ib = (int)floor(bin_f);
    float w = bin_f - (float)ib;
    int off = pair_type * n_bins;
    if (ib >= 0 && ib < n_bins) {
        atomic_add_float(&hist[off + ib], weight * (1.0f - w));
    }
    if (ib + 1 >= 0 && ib + 1 < n_bins) {
        atomic_add_float(&hist[off + ib + 1], weight * w);
    }
}

/// sinc(x) = sin(x)/x with small-x Taylor expansion
inline float sinc(float x) {
    float ax = fabs(x);
    if (ax < 1.0e-4f) {
        float x2 = x * x;
        return 1.0f - x2 * (1.0f / 6.0f - x2 / 120.0f);
    }
    return sin(x) / x;
}

/// Pair-distance histogram with linear interpolation.
/// Each work-item handles atom i and loops over j>i.
///
/// Args:
///   atoms     — float4[n_atoms]: xyz + element_index (as float, cast to int)
///   hist      — float[n_pair_types * n_bins]: global atomic accumulation target
///   n_atoms   — total atoms
///   n_species — number of element types
///   n_bins    — histogram bins
///   r_min     — lower histogram edge (Å)
///   dr        — bin width (Å)
///   n_pair_types — n_species*(n_species+1)/2
__kernel void pair_histogram(
    __global const float4* atoms,
    __global float* hist,
    const int n_atoms,
    const int n_species,
    const int n_bins,
    const float r_min,
    const float dr,
    const int n_pair_types
) {
    int i = get_global_id(0);
    if (i >= n_atoms) return;

    float3 ri = atoms[i].xyz;
    int ti = (int)(atoms[i].w);

    for (int j = i + 1; j < n_atoms; j++) {
        float3 rj = atoms[j].xyz;
        int tj = (int)(atoms[j].w);

        float3 d = ri - rj;
        float r = sqrt(d.x * d.x + d.y * d.y + d.z * d.z);

        // sort species indices so ti <= tj (triangular indexing)
        int a = ti, b = tj;
        if (a > b) { int tmp = a; a = b; b = tmp; }

        // pair_type = a * n_species + b - a*(a+1)/2
        int pair_type = a * n_species + b - (a * (a + 1)) / 2;
        if (pair_type < 0 || pair_type >= n_pair_types) continue;

        float bin_f = (r - r_min) / dr;
        deposit_pair(hist, n_bins, pair_type, bin_f, 1.0f);
    }
}

/// Gaussian-splat histogram (thermal broadening).
/// Each pair contributes a Gaussian centred at r0 with width sigma_ij.
/// Only bins within ±3σ are touched.
__kernel void pair_histogram_gaussian(
    __global const float4* atoms,
    __global const float* sigma,   // per-pair sigma, flat upper-triangle or full N*N
    __global float* hist,
    const int n_atoms,
    const int n_species,
    const int n_bins,
    const float r_min,
    const float dr,
    const int n_pair_types,
    const int sigma_mode  // 0=use sigma buffer, 1=constant sigma per species-pair
) {
    int i = get_global_id(0);
    if (i >= n_atoms) return;

    float3 ri = atoms[i].xyz;
    int ti = (int)(atoms[i].w);

    for (int j = i + 1; j < n_atoms; j++) {
        float3 rj = atoms[j].xyz;
        int tj = (int)(atoms[j].w);

        float3 d = ri - rj;
        float r0 = sqrt(d.x * d.x + d.y * d.y + d.z * d.z);

        int a = ti, b = tj;
        if (a > b) { int tmp = a; a = b; b = tmp; }
        int pair_type = a * n_species + b - (a * (a + 1)) / 2;
        if (pair_type < 0 || pair_type >= n_pair_types) continue;

        float sig = 0.0f;
        if (sigma_mode == 0) {
            sig = sigma[i * n_atoms + j];  // full matrix lookup
        }

        if (sig < 1.0e-6f) {
            // fallback to delta (linear interp)
            float bin_f = (r0 - r_min) / dr;
            deposit_pair(hist, n_bins, pair_type, bin_f, 1.0f);
        } else {
            float inv_s = 1.0f / sig;
            float r_low  = r0 - 3.0f * sig;
            float r_high = r0 + 3.0f * sig;
            int k0 = max(0, (int)floor((r_low  - r_min) / dr));
            int k1 = min(n_bins - 1, (int)ceil ((r_high - r_min) / dr));

            float norm = 0.0f;
            // two-pass: first accumulate weights locally to normalise
            float ws[32];  // stencil width <= 6σ/dr; 32 is safe for σ<=5*dr
            int nst = 0;
            for (int k = k0; k <= k1; k++) {
                float rb = r_min + (k + 0.5f) * dr;
                float x = (rb - r0) * inv_s;
                float w = exp(-0.5f * x * x);
                ws[nst++] = w;
                norm += w;
            }
            if (norm > 1.0e-12f) {
                float inv_norm = 1.0f / norm;
                int off = pair_type * n_bins;
                nst = 0;
                for (int k = k0; k <= k1; k++) {
                    atomic_add_float(&hist[off + k], ws[nst++] * inv_norm);
                }
            }
        }
    }
}

/// Debye transform: histogram -> I(Q)
///
/// Each work-item handles one Q-point.
///
/// Args:
///   hist      — float[n_pair_types * n_bins]
///   ff_prod   — float[n_pair_types * n_Q]: precomputed f_a(Q)*f_b(Q) for each pair type
///   Q_vals    — float[n_Q]: scattering vector magnitudes
///   r_centers — float[n_bins]: bin centre distances
///   I_Q       — float[n_Q]: output intensity
///   n_pair_types
///   n_bins
///   n_Q
__kernel void debye_transform_q(
    __global const float* hist,
    __global const float* ff_prod,
    __global const float* Q_vals,
    __global const float* r_centers,
    __global float* I_Q,
    const int n_pair_types,
    const int n_bins,
    const int n_Q
) {
    int iq = get_global_id(0);
    if (iq >= n_Q) return;

    float Q = Q_vals[iq];
    float I = 0.0f;
    for (int pt = 0; pt < n_pair_types; pt++) {
        float ff = ff_prod[pt * n_Q + iq];
        if (fabs(ff) < 1.0e-12f) continue;
        float sum = 0.0f;
        int off = pt * n_bins;
        for (int k = 0; k < n_bins; k++) {
            float qr = Q * r_centers[k];
            float s = sinc(qr);
            sum += hist[off + k] * s;
        }
        I += ff * sum;
    }
    I_Q[iq] = I;
}

