
float interpolate_2d(
    float xin, float yin,
    int nx, int ny, float hx, float hy,
    const __global float* data, // [nx, ny, n_nz]
    int i_nz, int n_nz
) {
    float inv_hx = 1.0f / hx;
    float inv_hy = 1.0f / hy;

    // Fortran interpolate_2d:
    //   imidx = int(xin/hx); clamp to [1, nx-1]
    //   imidy = int(yin/hy); clamp to [1, ny-1]
    int imidx = (int)(xin * inv_hx);
    if (imidx < 1) imidx = 1; else if (imidx > nx - 1) imidx = nx - 1;
    int imidy = (int)(yin * inv_hy);
    if (imidy < 1) imidy = 1; else if (imidy > ny - 1) imidy = ny - 1;

    float px = xin * inv_hx - imidx;
    float py = yin * inv_hy - imidy;

    // NOTE: data_3c is packed as [ny, nx, n_nz] i.e. [y,x] in the OpenCL buffer.
    // Therefore linear indexing here is (iy*nx + ix)*n_nz + i_nz.

    float g[4];
    for (int k = 0; k < 4; k++) {
        // Fortran: k=1..4, ikF = k + imidx - 1 (1-based)
        // Here:   k=0..3, ik0 = (k+1 + imidx - 1) - 1 = k + imidx - 1
        int ik = clamp(k + imidx - 1, 0, nx - 1);
        int iy0 = imidy - 1;
        int iy1 = clamp(iy0 + 1, 0, ny - 1);
        int iy2 = clamp(iy0 + 2, 0, ny - 1);
        int iy3 = clamp(iy0 + 3, 0, ny - 1);
        iy0 = clamp(iy0, 0, ny - 1);

        float f1m1 = data[(iy0 * nx + ik) * n_nz + i_nz];
        float f0p3 = 3.0f * data[(iy1 * nx + ik) * n_nz + i_nz];
        float f1p3 = 3.0f * data[(iy2 * nx + ik) * n_nz + i_nz];
        float f2p1 = data[(iy3 * nx + ik) * n_nz + i_nz];

        float bb3 = -f1m1 + f0p3 - f1p3 + f2p1;
        float bb2 = 3.0f * f1m1 - 2.0f * f0p3 + f1p3;
        float bb1 = -2.0f * f1m1 - f0p3 + 2.0f * f1p3 - f2p1;
        float bb0 = 2.0f * f0p3;
        g[k] = ((bb3 * py + bb2) * py + bb1) * py + bb0;
    }

    float f1m1 = g[0];
    float f0p3 = 3.0f * g[1];
    float f1p3 = 3.0f * g[2];
    float f2p1 = g[3];

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
    const int n_isorp_max, const int isorp_idx,
    const __global float3* ratoms,    // [natoms]
    const __global int4* triplets,    // [n_triplets_work] (i, j, k, type_idx)
    const __global float* data_3c,    // [n_triplets, n_isorp, 5, ny, nx, n_nz]
    const __global float2* h_grids_3c,// [n_triplets, n_isorp] (hx, hy)
    const __global int2* dims_3c,     // [n_triplets, n_isorp] (numx, numy)
    __global float* results           // [n_triplets_work, n_nz_max]
    // TODO: mvalue for bcna (when needed)
    // TODO: rotation of hlist into matrix elements
    // TODO: ensure n_nz_max matches index_max3c(in1,in2)
    // TODO: match Fortran interpolation clamping and bounds
) {
    int idx = get_global_id(0);
    if (idx >= n_triplets_work) return;

    int4 ijk_type = triplets[idx];
    float3 r_i = ratoms[ijk_type.x];
    float3 r_j = ratoms[ijk_type.y];
    float3 r_k = ratoms[ijk_type.z];
    int type_idx = ijk_type.w;

    // Fortran geometry: y = |r2 - r1|; x = |rna - 0.5*(r1+r2)| ; cost = dot(sighat, rhat)
    float3 r21 = r_j - r_i;
    float3 rnabc = r_k - 0.5f*(r_i + r_j);
    float y = length(r21);
    float x = length(rnabc);
    
    float3 sighat = (y > 1e-12f) ? r21 / y : (float3)(0.0f,0.0f,1.0f);
    float3 rhat  = (x > 1e-12f) ? rnabc / x : (float3)(0.0f,0.0f,1.0f);
    float cost = dot(sighat, rhat);
    if (cost > 1.0f) cost = 1.0f; else if (cost < -1.0f) cost = -1.0f;
    float cost2 = cost * cost;
    float sint = sqrt(max(0.0f, 1.0f - cost2));

    int grid_idx = type_idx * n_isorp_max + isorp_idx;
    float2 h = h_grids_3c[grid_idx];
    float hx = h.x;
    float hy = h.y;
    int2 dims = dims_3c[grid_idx];
    int nx_use = dims.x;
    int ny_use = dims.y;

    float p[5];
    p[0] = 1.0f;
    p[1] = cost;
    p[2] = (3.0f * cost2 - 1.0f) * 0.5f;
    p[3] = (5.0f * cost2 * cost - 3.0f * cost) * 0.5f;
    p[4] = (35.0f * cost2 * cost2 - 30.0f * cost2 + 3.0f) * 0.125f;

    int triplet_stride = n_isorp_max * 5 * ny_max * nx_max * n_nz_max;
    int isorp_stride = 5 * ny_max * nx_max * n_nz_max;
    int theta_stride = ny_max * nx_max * n_nz_max;

    for (int i_nz = 0; i_nz < n_nz_max; i_nz++) {
        float hlist = 0.0f;
        for (int it = 0; it < 5; it++) {
            const __global float* subdata = &data_3c[type_idx * triplet_stride + isorp_idx * isorp_stride + it * theta_stride];
            float val = interpolate_2d(x, y, nx_use, ny_use, hx, hy, subdata, i_nz, n_nz_max);
            hlist += p[it] * val;
        }
        // Simplified: Fireball checks mvalue here. Assume m=0 for den3.
        // For bcna, m can be 1. We'll need mvalue passed in.
        results[idx * n_nz_max + i_nz] = hlist;
    }
}

// Batch probe of 3c points; one work-item per point
__kernel void scan_3c_points(
    const int npoints,
    const int n_nz_max,
    const int nx_max, const int ny_max,
    const int n_isorp_max, const int isorp_idx,
    const int type_idx,
    const __global float* dRjs,   // [npoints,3]
    const __global float* dRks,   // [npoints,3]
    const __global float* data_3c,
    const __global float2* h_grids_3c,
    const __global int2* dims_3c,
    __global float* results        // [npoints, n_nz_max]
) {
    int idx = get_global_id(0);
    if (idx >= npoints) return;

    float3 r_i = (float3)(0.0f,0.0f,0.0f);
    float3 r_j = vload3(idx, dRjs);
    float3 r_k = vload3(idx, dRks);

    // Fortran geometry (firecore_scanHamPiece3c):
    // r1=(0,0,0), r2=dRj, rna=dRk
    // y = |r2 - r1|, x = |rna - 0.5*(r1+r2)|, cost = dot(sighat, rhat)
    float3 r21 = r_j - r_i;
    float y = length(r21);
    float inv_y = (y > 1e-12f) ? 1.0f / y : 0.0f;

    float3 mid = 0.5f * (r_i + r_j);
    float3 rnabc = r_k - mid;
    float x = length(rnabc);
    float inv_x = (x > 1e-12f) ? 1.0f / x : 0.0f;

    float3 sighat = (inv_y > 0.0f) ? (r21 * inv_y) : (float3)(0.0f,0.0f,1.0f);
    float3 rhat  = (inv_x > 0.0f) ? (rnabc * inv_x) : (float3)(0.0f,0.0f,1.0f);
    float cost = dot(sighat, rhat);
    if (cost > 1.0f) cost = 1.0f; else if (cost < -1.0f) cost = -1.0f;
    float cost2 = cost * cost;

    int grid_idx = type_idx * n_isorp_max + isorp_idx;
    float2 h = h_grids_3c[grid_idx];
    float hx = h.x;
    float hy = h.y;
    int2 dims = dims_3c[grid_idx];
    int nx_use = dims.x;
    int ny_use = dims.y;

    float p[5];
    p[0] = 1.0f;
    p[1] = cost;
    p[2] = (3.0f * cost2 - 1.0f) * 0.5f;
    p[3] = (5.0f * cost2 * cost - 3.0f * cost) * 0.5f;
    p[4] = (35.0f * cost2 * cost2 - 30.0f * cost2 + 3.0f) * 0.125f;

    int triplet_stride = n_isorp_max * 5 * ny_max * nx_max * n_nz_max;
    int isorp_stride = 5 * ny_max * nx_max * n_nz_max;
    int theta_stride = ny_max * nx_max * n_nz_max;

    for (int i_nz = 0; i_nz < n_nz_max; i_nz++) {
        float hlist = 0.0f;
        for (int it = 0; it < 5; it++) {
            const __global float* subdata = &data_3c[type_idx * triplet_stride + isorp_idx * isorp_stride + it * theta_stride];
            float val = interpolate_2d(x, y, nx_use, ny_use, hx, hy, subdata, i_nz, n_nz_max);
            hlist += p[it] * val;
        }
        results[idx * n_nz_max + i_nz] = hlist;
    }
}

// Raw 3c probe: returns the five bcna_0k components (no Legendre sum/rotation)
__kernel void scan_3c_raw_points(
    const int npoints,
    const int n_nz_max,
    const int nx_max, const int ny_max,
    const int n_isorp_max, const int isorp_idx,
    const int type_idx,
    const __global float* dRjs,   // [npoints,3]
    const __global float* dRks,   // [npoints,3]
    const __global float* data_3c,
    const __global float2* h_grids_3c,
    const __global int2* dims_3c,
    __global float* results        // [npoints, 5, n_nz_max]
) {
    int idx = get_global_id(0);
    if (idx >= npoints) return;


    if (idx != 0) return; // For the moment we do it serial

    const int debug_limit = 64;
    for (int idx = 0; idx < npoints; idx++) { // for debugging we iterate over points allowing serial print (better for debugging)

    float3 r_i = (float3)(0.0f,0.0f,0.0f);
    float3 r_j = vload3(idx, dRjs);
    float3 r_k = vload3(idx, dRks);

    float3 r21 = r_j - r_i;
    float y = length(r21);
    float3 rnabc = r_k - 0.5f*(r_i + r_j);
    float x = length(rnabc);

    int grid_idx = type_idx * n_isorp_max + isorp_idx;
    float2 h = h_grids_3c[grid_idx];
    float hx = h.x;
    float hy = h.y;
    int2 dims = dims_3c[grid_idx];
    int nx_use = dims.x;
    int ny_use = dims.y;

    int triplet_stride = n_isorp_max * 5 * ny_max * nx_max * n_nz_max;
    int isorp_stride = 5 * ny_max * nx_max * n_nz_max;
    int theta_stride = ny_max * nx_max * n_nz_max;
    const __global float* triplet_base = data_3c + type_idx * triplet_stride + isorp_idx * isorp_stride;

    if (idx < debug_limit) {
        // Match Fortran format (fields + 1-based ip)
        // Fortran prints: [scanHamPiece3c_raw_batch] ip=   NN      x    y    hx  hy   nx  ny
        printf("[OCL::scan_3c_raw_points] ip=%5d  %12.6f %12.6f  %9.6f %9.6f  %5d %5d\n",
               idx+1, x, y, hx, hy, nx_use, ny_use);
    }
    

    for (int it = 0; it < 5; it++) {
        const __global float* subdata = triplet_base + it * theta_stride;
        for (int i_nz = 0; i_nz < n_nz_max; i_nz++) {
            float val = interpolate_2d(x, y, nx_use, ny_use, hx, hy, subdata, i_nz, n_nz_max);
            results[(idx * 5 + it) * n_nz_max + i_nz] = val;
        }
    }

    }
}

typedef struct {
    float4 a, b, c, d; // Spline coefficients for one point
} SplineCoeffs;

// Batch probe of 2c blocks; one work-item per point (for verification / scanning)
__kernel void scan_2c_points(
    const int npoints,
    const int n_nonzero_max,
    const int numz,
    const int pair_type,
    const int applyRotation,
    const __global float* dRs,     // [npoints,3]
    const __global float* splines, // [n_pairs, numz, n_nonzero_max, 4]
    const __global float* h_grids, // [n_pairs]
    __global float* blocks         // [npoints, 4, 4]
) {
    int idx = get_global_id(0);
    if (idx >= npoints) return;

    float3 dR = vload3(idx, dRs);
    float r = length(dR);
    float3 ez;
    if (applyRotation != 0) {
        if (r > 1e-10f) ez = dR / r; else ez = (float3)(0.0f, 0.0f, 1.0f);
    } else {
        ez = (float3)(0.0f, 0.0f, 1.0f);
        dR = (float3)(0.0f, 0.0f, r);
    }

    float h_grid = h_grids[pair_type];
    int i_z = (int)(r / h_grid);
    if (i_z < 0) i_z = 0;
    if (i_z >= numz - 1) i_z = numz - 2;
    float dr = r - i_z * h_grid;

    int pair_stride = numz * n_nonzero_max * 4;
    int z_stride = n_nonzero_max * 4;
    float comps[6];
    for (int i_nz = 0; i_nz < n_nonzero_max; i_nz++) {
        int base = pair_type * pair_stride + i_z * z_stride + i_nz * 4;
        comps[i_nz] = splines[base + 0] + dr*(splines[base + 1] + dr*(splines[base + 2] + dr*splines[base + 3]));
    }

    float ss_sig = comps[0];
    float sp_sig = (n_nonzero_max > 1) ? comps[1] : 0.0f;
    float ps_sig = (n_nonzero_max > 2) ? comps[2] : 0.0f;
    float pp_pi  = (n_nonzero_max > 3) ? comps[3] : 0.0f;
    float pp_sig = (n_nonzero_max > 4) ? comps[4] : 0.0f;

    float ezx = ez.x;
    float ezy = -ez.y; // match Fortran sign convention (py)
    float ezz = -ez.z; // match Fortran sign convention (pz)
    __global float* b = &blocks[idx * 16];

    b[0*4 + 0] = ss_sig;
    b[0*4 + 1] = sp_sig * ezy;
    b[0*4 + 2] = sp_sig * ezz;
    b[0*4 + 3] = sp_sig * ezx;

    b[1*4 + 0] = ps_sig * ezy;
    b[2*4 + 0] = ps_sig * ezz;
    b[3*4 + 0] = ps_sig * ezx;

    float diff = pp_sig - pp_pi;
    b[1*4 + 1] = diff * ezy * ezy + pp_pi;
    b[1*4 + 2] = diff * ezy * ezz;
    b[1*4 + 3] = diff * ezy * ezx;

    b[2*4 + 1] = diff * ezz * ezy;
    b[2*4 + 2] = diff * ezz * ezz + pp_pi;
    b[2*4 + 3] = diff * ezz * ezx;

    b[3*4 + 1] = diff * ezx * ezy;
    b[3*4 + 2] = diff * ezx * ezz;
    b[3*4 + 3] = diff * ezx * ezx + pp_pi;
}

// Batch probe of 2c blocks; one work-item per pair (for assembly)
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
    // NOTE: Fortran reference uses opposite sign convention for pz; flip ez.z to match.
    float ezx = ez.x;
    float ezy = -ez.y; // DEBUG: align py sign with Fortran epsilon convention
    float ezz = -ez.z; // DEBUG: match Fortran pz sign
    __global float* b = &blocks[i_pair * 16];
    
    // s-s
    b[0*4 + 0] = ss_sig;
    
    // s-p
    b[0*4 + 1] = sp_sig * ezy; // py
    b[0*4 + 2] = sp_sig * ezz; // pz
    b[0*4 + 3] = sp_sig * ezx; // px
    
    // p-s
    b[1*4 + 0] = ps_sig * ezy;
    b[2*4 + 0] = ps_sig * ezz;
    b[3*4 + 0] = ps_sig * ezx;
    
    // p-p: H_ij = (pp_sig - pp_pi) * ez_i * ez_j + pp_pi * delta_ij
    float diff = pp_sig - pp_pi;
    
    // py-py, py-pz, py-px
    b[1*4 + 1] = diff * ezy * ezy + pp_pi;
    b[1*4 + 2] = diff * ezy * ezz;
    b[1*4 + 3] = diff * ezy * ezx;
    
    // pz-py, pz-pz, pz-px
    b[2*4 + 1] = diff * ezz * ezy;
    b[2*4 + 2] = diff * ezz * ezz + pp_pi;
    b[2*4 + 3] = diff * ezz * ezx;
    
    // px-py, px-pz, px-px
    b[3*4 + 1] = diff * ezx * ezy;
    b[3*4 + 2] = diff * ezx * ezz;
    b[3*4 + 3] = diff * ezx * ezx + pp_pi;
}
