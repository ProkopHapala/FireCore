
float interpolate_2d(
    float xin, float yin,
    int nx, int ny, float hx, float hy,
    const __global float* data, // [nx, ny, n_nz]
    int i_nz, int n_nz
) {
    float inv_hx = 1.0f / hx;
    float inv_hy = 1.0f / hy;

    int imidx = (int)(xin * inv_hx);
    if (imidx < 1) imidx = 1; else if (imidx > nx - 3) imidx = nx - 3;
    int imidy = (int)(yin * inv_hy);
    if (imidy < 1) imidy = 1; else if (imidy > ny - 3) imidy = ny - 3;

    float px = xin * inv_hx - imidx;
    float py = yin * inv_hy - imidy;

    float g[4];
    for (int k = 0; k < 4; k++) {
        int ik = k + imidx - 1;
        float f1m1 = data[(ik * ny + imidy    ) * n_nz + i_nz];
        float f0p3 = 3.0f * data[(ik * ny + imidy + 1) * n_nz + i_nz];
        float f1p3 = 3.0f * data[(ik * ny + imidy + 2) * n_nz + i_nz];
        float f2p1 = data[(ik * ny + imidy + 3) * n_nz + i_nz];

        float bb3 = -f1m1 + f0p3 - f1p3 + f2p1;
        float bb2 = 3.0f * f1m1 - 2.0f * f0p3 + f1p3;
        float bb1 = -2.0f * f1m1 - f0p3 + 2.0f * f1p3 - f2p1;
        float bb0 = 2.0f * f0p3;
        g[k] = ((bb3 * py + bb2) * py + bb1) * py + bb0;
    }

    float f1m1 = g[0];
    float f0p3 = 3.0f * g[1];
    float f1p3 = 3.0f * g[2];
    float f2p1 = g[4]; // Wait, g is size 4, so g[3]? Yes, indices 0,1,2,3.
    // Fixed:
    f2p1 = g[3];

    float bb3 = -f1m1 + f0p3 - f1p3 + f2p1;
    float bb2 = 3.0f * f1m1 - 2.0f * f0p3 + f1p3;
    float bb1 = -2.0f * f1m1 - f0p3 + 2.0f * f1p3 - f2p1;
    float bb0 = 2.0f * f0p3;
    return (((bb3 * px + bb2) * px + bb1) * px + bb0) / 36.0f;
}

__kernel void assemble_3c(
    const int n_triplets_work,
    const int n_nz_max,
    const int nx_max, const int ny_max,
    const __global float3* ratoms,    // [natoms]
    const __global int4* triplets,    // [n_triplets_work] (i, j, k, type_idx)
    const __global float* data_3c,    // [n_species_triplets, 5, ny, nx, n_nz]
    const __global float2* h_grids_3c, // [n_species_triplets] (hx, hy)
    __global float* results           // [n_triplets_work, n_nz_max]
) {
    int idx = get_global_id(0);
    if (idx >= n_triplets_work) return;

    int4 ijk_type = triplets[idx];
    float3 r_i = ratoms[ijk_type.x];
    float3 r_j = ratoms[ijk_type.y];
    float3 r_k = ratoms[ijk_type.z];
    int type_idx = ijk_type.w;

    float3 Rij = r_j - r_i;
    float3 Rik = r_k - r_i;
    float y = length(Rij);
    float x = length(Rik);
    
    float cost = dot(Rij, Rik) / (x * y + 1e-12f);
    if (cost > 1.0f) cost = 1.0f; else if (cost < -1.0f) cost = -1.0f;
    float cost2 = cost * cost;
    float sint = sqrt(max(0.0f, 1.0f - cost2));

    float2 h = h_grids_3c[type_idx];
    float hx = h.x;
    float hy = h.y;

    float p[5];
    p[0] = 1.0f;
    p[1] = cost;
    p[2] = (3.0f * cost2 - 1.0f) * 0.5f;
    p[3] = (5.0f * cost2 * cost - 3.0f * cost) * 0.5f;
    p[4] = (35.0f * cost2 * cost2 - 30.0f * cost2 + 3.0f) * 0.125f;

    int triplet_stride = 5 * ny_max * nx_max * n_nz_max;
    int theta_stride = ny_max * nx_max * n_nz_max;

    for (int i_nz = 0; i_nz < n_nz_max; i_nz++) {
        float hlist = 0.0f;
        for (int it = 0; it < 5; it++) {
            const __global float* subdata = &data_3c[type_idx * triplet_stride + it * theta_stride];
            float val = interpolate_2d(x, y, nx_max, ny_max, hx, hy, subdata, i_nz, n_nz_max);
            hlist += p[it] * val;
        }
        // Simplified: Fireball checks mvalue here. Assume m=0 for den3.
        // For bcna, m can be 1. We'll need mvalue passed in.
        results[idx * n_nz_max + i_nz] = hlist;
    }
}

typedef struct {
    float4 a, b, c, d; // Spline coefficients for one point
} SplineCoeffs;

__kernel void assemble_2c(
    const int n_pairs,
    const int n_nonzero_max,
    const int numz,
    const __global float3* ratoms,  // [natoms] - careful with float3 vs float4 padding!
    const __global int2* neighbors, // [n_pairs] (i, j)
    const __global float* splines,  // [n_species_pairs, numz, n_nonzero_max, 4]
    const __global int* pair_types, // [n_pairs] species pair index
    const __global float* h_grids,  // [n_species_pairs]
    __global float* blocks          // [n_pairs, 4, 4] block-sparse matrix
) {
    int i_pair = get_global_id(0);
    if (i_pair >= n_pairs) return;

    int2 ij = neighbors[i_pair];
    
    // Use float* to avoid float3 padding issues if not aligned
    // But since cl_float3 is often 16 bytes, we should be careful.
    // For now, let's just use float3 and make sure we pass float4-padded data from Python.
    float3 r_i = ratoms[ij.x];
    float3 r_j = ratoms[ij.y];
    
    float3 dR = r_j - r_i;
    float r = length(dR);
    
    float3 ez = (float3)(0.0f, 0.0f, 0.0f);
    if (r > 1e-10f) {
        ez = dR / r;
    }
    
    int spec_pair = pair_types[i_pair];
    float h_grid = h_grids[spec_pair];
    
    int i_z = (int)(r / h_grid);
    if (i_z < 0) i_z = 0;
    if (i_z >= numz - 1) i_z = numz - 2;
    float dr = r - i_z * h_grid;
    
    int pair_stride = numz * n_nonzero_max * 4;
    int z_stride = n_nonzero_max * 4;
    
    // Interpolate sigma/pi components
    float comps[6]; // Max 6 for s-p
    for (int i_nz = 0; i_nz < n_nonzero_max; i_nz++) {
        int base = spec_pair * pair_stride + i_z * z_stride + i_nz * 4;
        comps[i_nz] = splines[base + 0] + dr*(splines[base + 1] + dr*(splines[base + 2] + dr*splines[base + 3]));
    }
    
    // Map components to local matrix (s, py, pz, px order)
    // 1: ss_sig, 2: sp_sig, 3: ps_sig, 4: pp_pi, 5: pp_sig, 6: pp_pi
    float ss_sig = comps[0];
    float sp_sig = (n_nonzero_max > 1) ? comps[1] : 0.0f;
    float ps_sig = (n_nonzero_max > 2) ? comps[2] : 0.0f;
    float pp_pi  = (n_nonzero_max > 3) ? comps[3] : 0.0f;
    float pp_sig = (n_nonzero_max > 4) ? comps[4] : 0.0f;
    
    // Rotation assembly (Slater-Koster)
    // Orbital indices: 0:s, 1:py, 2:pz, 3:px
    __global float* b = &blocks[i_pair * 16];
    
    // s-s
    b[0*4 + 0] = ss_sig;
    
    // s-p
    b[0*4 + 1] = sp_sig * ez.y; // py
    b[0*4 + 2] = sp_sig * ez.z; // pz
    b[0*4 + 3] = sp_sig * ez.x; // px
    
    // p-s
    b[1*4 + 0] = ps_sig * ez.y;
    b[2*4 + 0] = ps_sig * ez.z;
    b[3*4 + 0] = ps_sig * ez.x;
    
    // p-p: H_ij = (pp_sig - pp_pi) * ez_i * ez_j + pp_pi * delta_ij
    float diff = pp_sig - pp_pi;
    
    // py-py, py-pz, py-px
    b[1*4 + 1] = diff * ez.y * ez.y + pp_pi;
    b[1*4 + 2] = diff * ez.y * ez.z;
    b[1*4 + 3] = diff * ez.y * ez.x;
    
    // pz-py, pz-pz, pz-px
    b[2*4 + 1] = diff * ez.z * ez.y;
    b[2*4 + 2] = diff * ez.z * ez.z + pp_pi;
    b[2*4 + 3] = diff * ez.z * ez.x;
    
    // px-py, px-pz, px-px
    b[3*4 + 1] = diff * ez.x * ez.y;
    b[3*4 + 2] = diff * ez.x * ez.z;
    b[3*4 + 3] = diff * ez.x * ez.x + pp_pi;
}
