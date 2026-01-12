
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

// ---------------------------
// PP rotation utilities (epsilon + twister p-matrix)
// ---------------------------

inline void epsilon_fb(float3 r1, float3 r2, __private float spe[3][3]){
    // Fortran epsilon(R1,R2,spe)
    // Here: r1 = r2 (absolute position of atom2), r2 = sighat (bond direction)
    float r1mag = length(r1);
    float r2mag = length(r2);
    // unit matrix if r2==0
    if(r2mag < 1e-4f){
        spe[0][0]=1.0f; spe[0][1]=0.0f; spe[0][2]=0.0f;
        spe[1][0]=0.0f; spe[1][1]=1.0f; spe[1][2]=0.0f;
        spe[2][0]=0.0f; spe[2][1]=0.0f; spe[2][2]=1.0f;
        return;
    }

    float3 zphat = r2 / r2mag;
    float3 yphat;
    float ypmag;

    // yphat = zphat x r1hat (cross with r1hat)
    if(r1mag > 1e-4f){
        float3 r1hat = r1 / r1mag;
        yphat = cross(zphat, r1hat);
        ypmag = length(yphat);
        if(ypmag > 1e-6f){
            yphat /= ypmag;
            float3 xphat = cross(yphat, zphat);
            // spe(ix,1)=xphat(ix), spe(ix,2)=yphat(ix), spe(ix,3)=zphat(ix) (Fortran 1-based)
            spe[0][0]=xphat.x; spe[0][1]=yphat.x; spe[0][2]=zphat.x;
            spe[1][0]=xphat.y; spe[1][1]=yphat.y; spe[1][2]=zphat.y;
            spe[2][0]=xphat.z; spe[2][1]=yphat.z; spe[2][2]=zphat.z;
            return;
        }
    }

    // fallback if colinear
    if(fabs(zphat.x) > 1e-4f){
        yphat = (float3)(-(zphat.y+zphat.z)/zphat.x, 1.0f, 1.0f);
    }else if(fabs(zphat.y) > 1e-4f){
        yphat = (float3)(1.0f, -(zphat.x+zphat.z)/zphat.y, 1.0f);
    }else{
        yphat = (float3)(1.0f, 1.0f, -(zphat.x+zphat.y)/zphat.z);
    }
    ypmag = length(yphat);
    yphat /= ypmag;
    float3 xphat = cross(yphat, zphat);
    spe[0][0]=xphat.x; spe[0][1]=yphat.x; spe[0][2]=zphat.x;
    spe[1][0]=xphat.y; spe[1][1]=yphat.y; spe[1][2]=zphat.y;
    spe[2][0]=xphat.z; spe[2][1]=yphat.z; spe[2][2]=zphat.z;
}

inline void twister_pmat(const __private float eps[3][3], __private float pmat[3][3]){
    // Fortran twister: pmat - y,z,x ordering
    pmat[0][0] = eps[1][1];
    pmat[0][1] = eps[1][2];
    pmat[0][2] = eps[1][0];

    pmat[1][0] = eps[2][1];
    pmat[1][1] = eps[2][2];
    pmat[1][2] = eps[2][0];

    pmat[2][0] = eps[0][1];
    pmat[2][1] = eps[0][2];
    pmat[2][2] = eps[0][0];
}

inline void rotatePP_sp(
    const __private float pmat[3][3],
    const __private float sm[4][4],
    __private float sx[4][4]
){
    // s+p only, Ortega order (s,py,pz,px)
    // L = diag(1, pmat) ; R = diag(1, pmat)
    // sx = L * sm * R^T
    float L[4][4];
    float R[4][4];
    for(int i=0;i<4;i++){ for(int j=0;j<4;j++){ L[i][j]=0.0f; R[i][j]=0.0f; } }
    L[0][0]=1.0f; R[0][0]=1.0f;
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            L[1+i][1+j] = pmat[i][j];
            R[1+i][1+j] = pmat[i][j];
        }
    }
    // temp = L*sm
    float tmp[4][4];
    for(int i=0;i<4;i++){
        for(int j=0;j<4;j++){
            float s=0.0f;
            for(int k=0;k<4;k++) s += L[i][k]*sm[k][j];
            tmp[i][j]=s;
        }
    }
    // sx = tmp * R^T
    for(int i=0;i<4;i++){
        for(int j=0;j<4;j++){
            float s=0.0f;
            for(int k=0;k<4;k++) s += tmp[i][k]*R[j][k];
            sx[i][j]=s;
        }
    }
}

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
    float ezy = -ez.y;
    float ezz = -ez.z;
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

// PP projector overlap (sVNL) assembly: Fortran doscentrosPP + recover_PP + rotatePP (s+p only)
__kernel void assemble_pp(
    const int n_pairs,
    const int n_nonzero_max,
    const int numz,
    const __global float3* ratoms,
    const __global int2* neighbors,
    const __global float* splines,
    const __global int* pair_types,
    const __global float* h_grids,
    const __global short* muPP_map,   // [n_species_pairs, n_nonzero_max] 1-based indices
    const __global short* nuPP_map,   // [n_species_pairs, n_nonzero_max] 1-based indices
    __global float* blocks            // [n_pairs, 4, 4]
) {
    int i_pair = get_global_id(0);
    if (i_pair >= n_pairs) return;

    // Safety: this kernel currently assumes n_nonzero_max <= 16
    if (n_nonzero_max > 16) return;

    int2 ij = neighbors[i_pair];
    float3 r_i = ratoms[ij.x];
    float3 r_j = ratoms[ij.y];
    float3 dR = r_j - r_i;
    float r = length(dR);
    float3 ez = (r > 1e-10f) ? dR / r : (float3)(0.0f,0.0f,1.0f);
    float3 sighat = ez;

    int spec_pair = pair_types[i_pair];
    float h_grid = h_grids[spec_pair];
    int i_z = (int)(r / h_grid);
    if (i_z < 0) i_z = 0;
    if (i_z >= numz - 1) i_z = numz - 2;
    float dr = r - i_z * h_grid;

    int pair_stride = numz * n_nonzero_max * 4;
    int z_stride = n_nonzero_max * 4;

    // Interpolate pplist(index)
    float vals[16];
    for (int i_nz = 0; i_nz < n_nonzero_max; i_nz++) {
        int base = spec_pair * pair_stride + i_z * z_stride + i_nz * 4;
        vals[i_nz] = splines[base + 0] + dr*(splines[base + 1] + dr*(splines[base + 2] + dr*splines[base + 3]));
    }

    // recover_PP into sm (molecular)
    float sm[4][4];
    for(int a=0;a<4;a++){ for(int b=0;b<4;b++){ sm[a][b]=0.0f; } }
    int mbase = spec_pair * n_nonzero_max;
    for (int idx = 0; idx < n_nonzero_max; idx++) {
        short imu1 = muPP_map[mbase + idx];
        short inu1 = nuPP_map[mbase + idx];
        if(imu1 <= 0 || inu1 <= 0) continue;
        int imu = (int)imu1 - 1;
        int inu = (int)inu1 - 1;
        if(imu < 4 && inu < 4) sm[imu][inu] = vals[idx];
    }

    // rotatePP (s+p)
    // Fortran assemble_sVNL uses epsilon(r2, sighat, eps) with r2 = absolute neighbor position
    float eps[3][3];
    epsilon_fb(r_j, sighat, eps);
    float pmat[3][3];
    twister_pmat(eps, pmat);
    float sx[4][4];
    rotatePP_sp(pmat, sm, sx);

    __global float* b = &blocks[i_pair * 16];
    for(int a=0;a<4;a++){
        for(int c=0;c<4;c++){
            // Match assemble_2c convention: store as (inu,imu)
            // so that Python's dense reconstruction (which transposes blocks) works consistently.
            b[c*4 + a] = sx[a][c];
        }
    }
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
    const __global int* is_vna_pair, // [n_species_pairs] 1 if Vna pair
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
    
    int is_vna = is_vna_pair[spec_pair];
    
    // Map components to local matrix (s, py, pz, px order)
    // 1: ss_sig, 2: sp_sig, 3: ps_sig, 4: pp_pi, 5: pp_sig, 6: pp_pi
    float ss_sig = comps[0];
    float sp_sig = comps[1];
    float ps_sig = (n_nonzero_max > 2) ? comps[2] : 0.0f;
    float pp_pi  = (n_nonzero_max > 3) ? comps[3] : 0.0f;
    float pp_sig = (n_nonzero_max > 4) ? comps[4] : 0.0f;

    // Fortran vna_atom path (interaction=4) uses epsilon(r2, sighat, eps) with r2=absolute neighbor position.
    // That makes the x/y axes depend on absolute r_j, not only on dR. This matters for sign/parity.
    // We implement that strictly for vna_atom pairs (flagged in is_vna_pair buffer).
    if(is_vna != 0){
        // local (bond-frame) matrix sm in Ortega order (s,py,pz,px) where pz is along bond
        float sm[4][4];
        for(int a=0;a<4;a++){ for(int c=0;c<4;c++){ sm[a][c]=0.0f; } }
        sm[0][0] = ss_sig;
        sm[0][2] = sp_sig;
        sm[2][0] = ps_sig;
        sm[1][1] = pp_pi;
        sm[2][2] = pp_sig;
        sm[3][3] = pp_pi;

        float eps[3][3];
        epsilon_fb(r_j, ez, eps);
        float pmat[3][3];
        twister_pmat(eps, pmat);
        float sx[4][4];
        rotatePP_sp(pmat, sm, sx);

        __global float* b = &blocks[i_pair * 16];
        for(int a=0;a<4;a++){
            for(int c=0;c<4;c++){
                b[a*4 + c] = sx[a][c];
            }
        }
        return;
    }
    
    // Rotation assembly (Slater-Koster)
    // Orbital indices: 0:s, 1:py, 2:pz, 3:px
    float ezx = -ez.x;
    float ezy = -ez.y;
    float ezz = -ez.z;
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




// hamiltonian_vnl_shared.cl

// Compile-time constants (passed via -D to PyOpenCL)
#define N_ORB 4     // e.g., 4 for sp, 9 for spd
#define N_PROJ 4    // e.g., 4 for sp, 9 for spd
#define N_ORB_SQ 16 // N_ORB * N_ORB

/**
 * Kernel: contract_vnl_shared
 * 
 * Logic:
 * 1. Map Workgroup -> Neighbor Pair (i, j)
 * 2. Cooperative Load: Threads load S_ii, S_jj, S_ij, S_ji into Local Memory.
 * 3. Compute: H_ij = Term1 (Proj on i) + Term2 (Proj on j)
 * 4. Write: H_out[pair_index] (No atomics)
 */
__kernel void contract_vnl_shared(
    const int n_pairs,
    __global const int2* neighbor_pairs,   // (i, j)
    __global const int* self_offset_map,   // atom_idx -> index of (i, i) in sVNL
    __global const int* reverse_pair_map,  // pair_idx -> index of (j, i) in sVNL
    __global const float* sVNL_global,     // [total_pairs, N_ORB, N_PROJ]
    __global const float* cl_coeffs,       // [n_atoms, N_PROJ]
    __global float* H_vnl_out              // [total_pairs, N_ORB, N_ORB]
) {
    // 1. Identify Workgroup Task
    int pair_id = get_group_id(0); // Each WG handles one pair
    if (pair_id >= n_pairs) return;

    // Thread ID within the block (0..15 for sp, 0..80 for spd)
    int tid = get_local_id(0); 
    
    // Identify dimensions (row/col of the output matrix)
    int row = tid / N_ORB;
    int col = tid % N_ORB;

    // --- LOCAL MEMORY ALLOCATION ---
    // We need 4 matrices of size [N_ORB][N_PROJ]
    // For sp (4x4): 16 floats * 4 matrices = 64 floats (256 bytes) -> Fits easily.
    __local float S_ii[N_ORB * N_PROJ];
    __local float S_jj[N_ORB * N_PROJ];
    __local float S_ij[N_ORB * N_PROJ];
    __local float S_ji[N_ORB * N_PROJ];

    // --- DATA LOADING PHASE ---
    
    int2 atoms = neighbor_pairs[pair_id];
    int atom_i = atoms.x;
    int atom_j = atoms.y;

    int idx_S_ii = self_offset_map[atom_i];
    int idx_S_jj = self_offset_map[atom_j];
    int idx_S_ij = pair_id;               // The current pair
    int idx_S_ji = reverse_pair_map[pair_id]; // The reverse pair

    // Threads cooperate to load data. 
    // Since N_ORB*N_PROJ might be != LocalSize, we loop.
    int total_elements = N_ORB * N_PROJ;
    
    for (int i = tid; i < total_elements; i += get_local_size(0)) {
        // Flat indices for global buffer
        int b_row = i / N_PROJ;
        int b_col = i % N_PROJ;
        
        // sVNL layout: [pair][basis_row][proj_col]
        // Stride = total_elements
        
        S_ii[i] = sVNL_global[idx_S_ii * total_elements + i];
        S_jj[i] = sVNL_global[idx_S_jj * total_elements + i];
        S_ij[i] = sVNL_global[idx_S_ij * total_elements + i];
        
        // Handle case where reverse pair might not exist (boundary conditions)
        if (idx_S_ji >= 0) {
            S_ji[i] = sVNL_global[idx_S_ji * total_elements + i];
        } else {
            S_ji[i] = 0.0f; 
        }
    }

    // Barrier: Wait for all LMEM to be populated
    barrier(CLK_LOCAL_MEM_FENCE);

    // --- COMPUTE PHASE ---
    
    // Check if this thread maps to a valid matrix element
    if (row < N_ORB && col < N_ORB) {
        
        float val = 0.0f;
        
        // Pointers to coefficients
        __global const float* C_i = &cl_coeffs[atom_i * N_PROJ];
        __global const float* C_j = &cl_coeffs[atom_j * N_PROJ];

        // Contraction Loop over Projectors (k)
        for (int k = 0; k < N_PROJ; k++) {
            
            // Term 1: Projector on Atom J
            // <phi_i | Psi_j> * C_j * <Psi_j | phi_j>
            // Matrix: S_ij * diag(C_j) * S_jj^T
            // Element: S_ij[row, k] * C_j[k] * S_jj[col, k] 
            // Note: S_jj[col, k] corresponds to Transpose(S_jj)[k, col]
            
            float term_j = S_ij[row * N_PROJ + k] * C_j[k] * S_jj[col * N_PROJ + k];

            // Term 2: Projector on Atom I
            // <phi_i | Psi_i> * C_i * <Psi_i | phi_j>
            // Matrix: S_ii * diag(C_i) * S_ji^T
            // Note: We need <Psi_i | phi_j>. 
            // In the reverse pair (j->i), sVNL stored <phi_j | Psi_i>.
            // So <Psi_i | phi_j> is the TRANSPOSE of sVNL(j,i).
            // So we need (S_ji)^T. Element at [row, k] of S_ji^T is S_ji[k, row]? 
            // No.
            // Let's look at indices:
            // We want element (row, col) of the result.
            // S_ii[row, k] * C_i[k] * (Something relating i and j)[k, col]
            // The last term is <Psi_i_k | phi_j_col>.
            // The stored block S_ji is <phi_j | Psi_i>.
            // Element S_ji[basis_idx, proj_idx].
            // We want <Psi_i_k | phi_j_col> = Conjugate(<phi_j_col | Psi_i_k>).
            // This is S_ji[col, k].
            
            float term_i = S_ii[row * N_PROJ + k] * C_i[k] * S_ji[col * N_PROJ + k];

            val += (term_j + term_i);
        }

        // --- WRITE PHASE ---
        // Writes are coalesced if N_ORB is multiple of 16/32, 
        // or at least grouped by workgroup.
        int out_idx = pair_id * (N_ORB * N_ORB) + row * N_ORB + col;
        H_vnl_out[out_idx] = val;
    }
}

// NOTE:
// The sketch kernel above assumes that pair_id indexes directly into sVNL_global (idx_S_ij=pair_id).
// In our Python path, sVNL is typically evaluated for a separate list of (phi_atom, pp_atom) pairs
// (often all (i,k)), whose ordering does NOT match the 2-center neighbor list.
// Therefore we provide a map-driven kernel below which explicitly takes indices of the needed sVNL blocks.
//
// Conventions:
// - sVNL blocks are stored in the same orientation as assemble_pp output: block[nu,mu] (row=projector index, col=basis index).
// - Vnl output blocks should follow the same convention as other 2c blocks in this file: block[nu,mu].
// - The contraction for one term is:
//     Vnl(nu,mu) += sum_k sVNL_A(k,mu) * cl(k) * sVNL_B(k,nu)
//   where sVNL_A corresponds to <phi_i|Psi_K> and sVNL_B corresponds to <phi_j|Psi_K>.

__kernel void contract_vnl_pairs_sp(
    const int n_pairs,
    __global const int2* ij_pairs,          // [n_pairs] (i,j)
    __global const int* idx_ii,             // [n_pairs] index into sVNL_global for (i,i)
    __global const int* idx_ij,             // [n_pairs] index into sVNL_global for (i,j)
    __global const int* idx_ji,             // [n_pairs] index into sVNL_global for (j,i)
    __global const int* idx_jj,             // [n_pairs] index into sVNL_global for (j,j)
    __global const float* sVNL_global,      // [n_sVNL_pairs, 16] flattened blocks (nu,mu)
    __global const float* cl_coeffs,        // [n_atoms, 4] expanded per-PP-orb coefficients (s,py,pz,px)
    __global float* vnl_out                 // [n_pairs, 16] output blocks (nu,mu)
){
    // One work-group per pair, local size 16.
    int gid  = get_global_id(0);
    int tid  = get_local_id(0);
    int pair = gid >> 4;                   // /16
    if(pair >= n_pairs) return;

    // Map tid -> (nu,mu) for output element
    int nu = tid >> 2;                     // /4
    int mu = tid & 3;                      // %4

    int2 ij = ij_pairs[pair];
    int i = ij.x;
    int j = ij.y;

    int i_ii = idx_ii[pair];
    int i_ij = idx_ij[pair];
    int i_ji = idx_ji[pair];
    int i_jj = idx_jj[pair];

    __local float Sii[16];
    __local float Sij[16];
    __local float Sji[16];
    __local float Sjj[16];

    // Load all 4 blocks into local memory (each thread loads one element from each)
    // If any index is negative (missing reverse/self), treat as zero.
    Sii[tid] = (i_ii >= 0) ? sVNL_global[(i_ii<<4) + tid] : 0.0f;
    Sij[tid] = (i_ij >= 0) ? sVNL_global[(i_ij<<4) + tid] : 0.0f;
    Sji[tid] = (i_ji >= 0) ? sVNL_global[(i_ji<<4) + tid] : 0.0f;
    Sjj[tid] = (i_jj >= 0) ? sVNL_global[(i_jj<<4) + tid] : 0.0f;

    barrier(CLK_LOCAL_MEM_FENCE);

    __global const float* Ci = cl_coeffs + (i<<2);
    __global const float* Cj = cl_coeffs + (j<<2);

    // Compute Vnl(nu,mu) for this pair
    float val = 0.0f;

    // PP at i: <phi_i|Psi_i> * cl(i) * <phi_j|Psi_i>
    // Uses Sii(k,mu) and Sji(k,nu)
    for(int k=0;k<4;k++){
        val += Sii[(k<<2) + mu] * Ci[k] * Sji[(k<<2) + nu];
    }

    // PP at j: <phi_i|Psi_j> * cl(j) * <phi_j|Psi_j>
    // Uses Sij(k,mu) and Sjj(k,nu)
    for(int k=0;k<4;k++){
        val += Sij[(k<<2) + mu] * Cj[k] * Sjj[(k<<2) + nu];
    }

    vnl_out[(pair<<4) + tid] = val;
}

// General contraction over PP centers k.
// For each pair (i,j), and each (nu,mu) element:
//   Vnl(nu,mu) = sum_k sum_cc sVNL(i,k)[cc,mu] * cl(k)[cc] * sVNL(j,k)[cc,nu]
// Here sVNL(i,k) is stored as block[nu,mu] => nu=cc (projector), mu=basis.
__kernel void contract_vnl_sumk_sp(
    const int n_pairs,
    const int n_k,
    __global const int2* ij_pairs,     // [n_pairs] (i,j)
    __global const int* idx_ik,        // [n_pairs*n_k] index of sVNL(i,k)
    __global const int* idx_jk,        // [n_pairs*n_k] index of sVNL(j,k)
    __global const float* sVNL_global, // [n_sVNL_pairs,16] blocks (cc,mu)
    __global const float* cl_coeffs,   // [n_atoms,4] expanded cl per projector orbital
    __global float* vnl_out            // [n_pairs,16] blocks (nu,mu)
){
    int gid = get_global_id(0);
    int tid = get_local_id(0);
    int pair = gid >> 4;
    if(pair >= n_pairs) return;

    int nu = tid >> 2;
    int mu = tid & 3;

    float val = 0.0f;
    int base = pair * n_k;
    for(int ik=0; ik<n_k; ik++){
        int a = idx_ik[base + ik];
        int b = idx_jk[base + ik];
        if(a < 0 || b < 0) continue;

        // k is encoded implicitly by ik: host should order k list consistently and provide cl in same order.
        // For current usage, host uses k=0..n_atoms-1, so k==ik.
        __global const float* ck = cl_coeffs + (ik<<2);
        int ao = a<<4;
        int bo = b<<4;
        for(int cc=0; cc<4; cc++){
            val += sVNL_global[ao + (cc<<2) + mu] * ck[cc] * sVNL_global[bo + (cc<<2) + nu];
        }
    }
    vnl_out[(pair<<4) + tid] = val;
}



// ================================================================================================
// HELPER FUNCTIONS: 3-CENTER ROTATION & SMOOTHING
// ================================================================================================

// Fireball Smoother (from smoother.f)
// Returns val in [0, 1] and derivative (optional, passed as float*)
inline float smoother_fb(float r, float rcut, float smt_elect) {
    // rcut corresponds to rend in Fortran
    // smt_elect corresponds to the width of the transition
    float rstart = rcut - smt_elect;
    
    if (r < rstart) return 1.0f;
    if (r > rcut) return 0.0f;

    // Polynomial smooth transition
    float x = (r - rstart) / (rcut - rstart);
    // Standard Fireball typically uses (1-x)^3 * (1+3x) or similar variations.
    // Here assuming the standard "sticky" smoother: (1-x*x)^2 or similar.
    // Let's use a standard C2 continuous smoother: 1 - 3x^2 + 2x^3 (Hermite) 
    // or (1-x)^3(1+3x) which is often used in tight-binding.
    // Mapping 0->1, 1->0:
    float y = 1.0f - x;
    return y * y * (3.0f - 2.0f * y); 
}

inline void eval_legendre_0_4(float cost, float* P){
    cost = clamp(cost, -1.0f, 1.0f);
    float c2 = cost*cost;
    P[0]=1.0f;
    P[1]=cost;
    P[2]=0.5f*(3.0f*c2-1.0f);
    P[3]=0.5f*(5.0f*c2*cost - 3.0f*cost);
    P[4]=0.125f*(35.0f*c2*c2 - 30.0f*c2 + 3.0f);
}

// Convert interpolated 3C coefficients (hlist) into a rotated 4x4 submatrix block
// Input: hlist[5] = {ss, sp_sig, ps_sig, pp_sig, pp_pi}
// Output: block[16] (row-major)
inline void recover_and_rotate_3c_sp(
    const float* hlist, 
    const __private float eps[3][3], 
    float* block // 16 elements
){
    // 1. Construct Molecular Frame Matrix M (z-axis along bond)
    // Indices: s=0, py=1, pz=2, px=3 (Ortega ordering)
    float M[4][4];
    for(int i=0; i<16; i++) ((float*)M)[i] = 0.0f;

    M[0][0] = hlist[0]; // s-s
    M[0][2] = hlist[1]; // s-pz
    M[2][0] = hlist[2]; // pz-s
    M[2][2] = hlist[3]; // pz-pz (sigma)
    M[1][1] = hlist[4]; // py-py (pi)
    M[3][3] = hlist[4]; // px-px (pi)

    // 2. Construct Rotation Matrix R (4x4) from eps (3x3)
    // eps rows are x', y', z' axes. Fireball P-ordering is y, z, x.
    float R[4][4];
    // S-orbital part (identity)
    R[0][0]=1.0f; R[0][1]=0.0f; R[0][2]=0.0f; R[0][3]=0.0f;
    R[1][0]=0.0f; R[2][0]=0.0f; R[3][0]=0.0f;

    // P-orbital part (Twister logic)
    // Local y(1) -> Global axes via eps row 1
    // Local z(2) -> Global axes via eps row 2
    // Local x(3) -> Global axes via eps row 0
    
    for(int k=0; k<3; k++) {
        R[1][k+1] = eps[1][k]; 
        R[2][k+1] = eps[2][k];
        R[3][k+1] = eps[0][k];
    }

    // 3. Apply Rotation: B = R * M * R^T
    // Compute tmp = R * M
    float tmp[4][4];
    
    for(int i=0; i<4; i++) {
        for(int j=0; j<4; j++) {
            float val = 0.0f;
            for(int k=0; k<4; k++) val += R[i][k] * M[k][j];
            tmp[i][j] = val;
        }
    }

    // Compute block = tmp * R^T
    for(int i=0; i<4; i++) {
        for(int j=0; j<4; j++) {
            float val = 0.0f;
            for(int k=0; k<4; k++) val += tmp[i][k] * R[j][k]; // R[j][k] is transpose(R)[k][j]
            block[i*4 + j] = val;
        }
    }
}

// ================================================================================================
// KERNEL 1: GATHER AVERAGE DENSITY (compute_avg_rho)
// ================================================================================================
__kernel void gather_density_3c(
    const int n_pairs,
    const int n_nz_max,
    // Geometry & Lists
    const __global float3* ratoms,
    const __global int4* pairs,       // (i, j, type_pair_idx, global_off)
    const __global int* cn_offsets,   // Common Neighbor Offsets (CSR)
    const __global int* cn_indices,   // Common Neighbor Indices (CSR)
    // Data
    const __global float* charges,    // Qin [natoms]
    const __global float* S_mat,      // Overlap matrix [n_pairs, 16] (flattened 4x4)
    const __global float* rho_2c,     // Precomputed 2-center density [n_pairs, 16]
    // Interpolation Tables (Interaction=3)
    const __global float* data_3c,
    const __global float2* h_grids_3c,
    const __global int2* dims_3c,
    const int n_isorp_max, const int nx_max, const int ny_max,
    // Output
    __global float* rho_avg_out       // [n_pairs, 16]
) {
    int pid = get_global_id(0);
    if (pid >= n_pairs) return;

    int4 pair_info = pairs[pid];
    int atom_i = pair_info.x;
    int atom_j = pair_info.y;
    int type_idx = pair_info.z; // Maps to species pair for tables

    float3 r_i = ratoms[atom_i];
    float3 r_j = ratoms[atom_j];
    
    // --- 1. Initialize with 2-center contribution ---
    float block_acc[16];
    const __global float* rho_2c_ptr = &rho_2c[pid * 16];
    for(int i=0; i<16; i++) block_acc[i] = rho_2c_ptr[i];

    // --- 2. Geometry Prep for Pair ---
    float3 r21 = r_j - r_i;
    float y = length(r21);
    float3 sighat = (y > 1e-8f) ? r21 / y : (float3)(0.0f,0.0f,1.0f);
    float eps[3][3];
    epsilon_fb(r_j, sighat, eps); // Setup rotation matrix for this pair

    // Legendre Polynomial buffer
    float P[5];
    float hlist[5]; // Storage for interpolated values (s-s, s-p, p-s, pp-sig, pp-pi)
    float rot_block[16];

    // --- 3. Gather Loop (3-Center) ---
    int start_k = cn_offsets[pid];
    int end_k   = cn_offsets[pid+1];

    for (int k_ptr = start_k; k_ptr < end_k; k_ptr++) {
        int atom_k = cn_indices[k_ptr];
        float q_k = charges[atom_k]; // Neutral atom charge Q_neutral assumed handled in pre-calc or deviation

        // 3C Geometry
        float3 r_k = ratoms[atom_k];
        float3 mid = 0.5f * (r_i + r_j);
        float3 rnabc = r_k - mid;
        float x = length(rnabc);
        float3 rhat = (x > 1e-8f) ? rnabc / x : (float3)(0.0f,0.0f,1.0f);
        float cost = dot(sighat, rhat);
        eval_legendre_0_4(cost, P);

        // Grid info
        // Assuming Interaction 3 (Density) uses isorp=0 (or specific shell mapping passed in uniforms)
        // Here simplifying assuming 1 component per interaction type or handling striding logic
        int isorp_idx = 0; // Simplified. Real code needs species(k) mapping.
        
        int grid_idx = type_idx * n_isorp_max + isorp_idx;
        float hx = h_grids_3c[grid_idx].x;
        float hy = h_grids_3c[grid_idx].y;
        int nx_use = dims_3c[grid_idx].x;
        int ny_use = dims_3c[grid_idx].y;

        // Packed layout: data_3c[triplet, isorp, itheta(0..4), y, x, nz]
        // Correct contraction: hlist[iME] = sum_{itheta} P[itheta] * interp(plane_itheta, nz=iME)
        int triplet_stride = n_isorp_max * 5 * ny_max * nx_max * n_nz_max;
        int isorp_stride   = 5 * ny_max * nx_max * n_nz_max;
        int theta_stride   = ny_max * nx_max * n_nz_max;
        const __global float* base_trip = &data_3c[type_idx * triplet_stride + isorp_idx * isorp_stride];
        for(int iME=0; iME<5; iME++){
            float sum = 0.0f;
            for(int it=0; it<5; it++){
                const __global float* plane = base_trip + it * theta_stride;
                sum += P[it] * interpolate_2d(x, y, nx_use, ny_use, hx, hy, plane, iME, n_nz_max);
            }
            hlist[iME] = sum;
        }

        // Rotate to Crystal Frame
        recover_and_rotate_3c_sp(hlist, eps, rot_block);

        // Accumulate Weighted by Charge
        for(int i=0; i<16; i++) {
            block_acc[i] += rot_block[i] * q_k;
        }
    }

    // --- 4. Divide by Overlap (S) ---
    const __global float* S_ptr = &S_mat[pid * 16];
    
    __global float* out_ptr = &rho_avg_out[pid * 16];
    
    for(int i=0; i<16; i++) {
        float s_val = S_ptr[i];
        // Avoid division by zero
        if (fabs(s_val) < 1e-8f) s_val = (s_val >= 0) ? 1e-8f : -1e-8f;
        out_ptr[i] = block_acc[i] / s_val;
    }
}

// ================================================================================================
// Average Density Assembly ( average_rho )
// ================================================================================================
__kernel void compute_avg_rho(
    const int n_pairs,
    const int n_nz_max, 
    // Geometry
    const __global float3* ratoms,
    const __global int4* pairs,       // [i, j, type_idx, _]
    const __global int* cn_offsets,   // CSR offsets
    const __global int* cn_indices,   // CSR column indices
    // Input Data
    const __global float* S_mat,      // Overlap [n_pairs, 16]
    const __global float* rho_2c,     // Precomputed 2-center density [n_pairs, 16]
    const __global float* Q_neutral,  // Neutral atomic charges [natoms]
    // Tables (Interaction type = 3)
    const __global float* data_3c,
    const __global float2* h_grids_3c,
    const __global int2* dims_3c,
    const int n_isorp_max, const int nx_max, const int ny_max,
    // Output
    __global float* rho_avg_out       // [n_pairs, 16]
) {
    int pid = get_global_id(0);
    if (pid >= n_pairs) return;

    // 1. Setup Pair Geometry
    int4 p = pairs[pid];
    int i = p.x; int j = p.y; int type_idx = p.z;
    float3 r_i = ratoms[i];
    float3 r_j = ratoms[j];
    float3 r21 = r_j - r_i;
    float d_ij = length(r21);
    float3 sighat = (d_ij > 1e-10f) ? r21/d_ij : (float3)(0,0,1);
    
    // Rotation Matrix
    float eps[3][3];
    epsilon_fb(r_j, sighat, eps); 

    // Initialize with 2-center part
    float rho_acc[16];
    for(int k=0; k<16; k++) rho_acc[k] = rho_2c[pid*16 + k];

    // 2. Loop Common Neighbors
    int start = cn_offsets[pid];
    int end   = cn_offsets[pid+1];

    for (int idx = start; idx < end; idx++) {
        int k = cn_indices[idx];
        float q_k = Q_neutral[k]; // Weight by neutral charge
        float3 r_k = ratoms[k];

        // 3C Coords
        float3 mid = 0.5f * (r_i + r_j);
        float3 rnabc = r_k - mid;
        float x = length(rnabc);
        float3 rhat = (x > 1e-10f) ? rnabc/x : (float3)(0,0,1);
        float cost = dot(sighat, rhat);

        // Legendre Polys
        float P[5];
        eval_legendre_0_4(cost, P);

        // Interpolate hlist[0..4] using correct Legendre sum over itheta planes
        float hlist[5];
        int isorp_idx = 0;
        int grid_idx = type_idx * n_isorp_max + isorp_idx;
        float hx = h_grids_3c[grid_idx].x;
        float hy = h_grids_3c[grid_idx].y;
        int nx = dims_3c[grid_idx].x;
        int ny = dims_3c[grid_idx].y;

        int trip_stride  = n_isorp_max * 5 * ny_max * nx_max * n_nz_max;
        int isorp_stride = 5 * ny_max * nx_max * n_nz_max;
        int theta_stride = ny_max * nx_max * n_nz_max;
        const __global float* base_trip = &data_3c[type_idx * trip_stride + isorp_idx * isorp_stride];

        for(int iME=0; iME<5; iME++){
            float sum = 0.0f;
            for(int it=0; it<5; it++){
                const __global float* plane = base_trip + it * theta_stride;
                sum += P[it] * interpolate_2d(x, d_ij, nx, ny, hx, hy, plane, iME, n_nz_max);
            }
            hlist[iME] = sum;
        }

        // Rotate and Accumulate
        float block[16];
        recover_and_rotate_3c_sp(hlist, eps, block);
        
        for(int m=0; m<16; m++) rho_acc[m] += block[m] * q_k;
    }

    // 3. Output accumulated density block (rho_off-like). Normalization (if any) is handled separately.
    __global float* out_ptr = &rho_avg_out[pid*16];
    for(int m=0; m<16; m++) out_ptr[m] = rho_acc[m];
}

// ================================================================================================
// Charged Atom Potential (assemble_ca_3c)
// ================================================================================================


__kernel void assemble_ca_3c(
    const int n_pairs,
    const int n_nz_max,
    const float r_cut_global, // Global smoothing cutoff
    const float smt_elect,    // Smoother width
    // Geometry
    const __global float3* ratoms,
    const __global int4* pairs,
    const __global int* cn_offsets,
    const __global int* cn_indices,
    // Data
    const __global float* delta_charges, // Qin - Qneutral
    const __global float* S_mat,         // Need overlap for long-range term
    // Tables (Interaction=CA)
    const __global float* data_3c,
    const __global float2* h_grids_3c,
    const __global int2* dims_3c,
    const int n_isorp_max, const int nx_max, const int ny_max,
    // Output
    __global float* H_3c_out             // Accumulates into here
) {
    int pid = get_global_id(0);
    if (pid >= n_pairs) return;

    int4 pair_info = pairs[pid];
    int atom_i = pair_info.x;
    int atom_j = pair_info.y;
    int type_idx = pair_info.z;

    float3 r_i = ratoms[atom_i];
    float3 r_j = ratoms[atom_j];
    float3 r21 = r_j - r_i;
    float y = length(r21);
    float3 sighat = (y > 1e-8f) ? r21 / y : (float3)(0.0f,0.0f,1.0f);
    float eps[3][3];
    epsilon_fb(r_j, sighat, eps);

    float block_acc[16];
    for(int i=0; i<16; i++) block_acc[i] = 0.0f;

    const __global float* S_ptr = &S_mat[pid * 16];

    int start_k = cn_offsets[pid];
    int end_k   = cn_offsets[pid+1];

    for (int k_ptr = start_k; k_ptr < end_k; k_ptr++) {
        int atom_k = cn_indices[k_ptr];
        float dq_k = delta_charges[atom_k];
        
        if (fabs(dq_k) < 1e-6f) continue; // Skip if neutral

        float3 r_k = ratoms[atom_k];
        float d_ik = distance(r_i, r_k);
        float d_jk = distance(r_j, r_k);

        // --- Smoother Calculation ---
        // Fireball uses combined smoother stn = stn(ik) * stn(jk)
        float stn1 = smoother_fb(d_ik, r_cut_global, smt_elect);
        float stn2 = smoother_fb(d_jk, r_cut_global, smt_elect);
        float stn = stn1 * stn2;
        float one_minus_stn = 1.0f - stn;

        // --- Short Range Term (Exact 3-Center Interpolation) ---
        float v_short[16];
        if (stn > 1e-6f) {
            // Geometry
            float3 mid = 0.5f * (r_i + r_j);
            float3 rnabc = r_k - mid;
            float x = length(rnabc);
            float3 rhat = (x > 1e-10f) ? rnabc / x : (float3)(0.0f,0.0f,1.0f);
            float cost = clamp(dot(sighat, rhat), -1.0f, 1.0f);
            float cost2 = cost*cost;
            
            float P[5];
            eval_legendre_0_4(cost, P);

            // Interpolate
            float hlist[5];
            // ... setup grid strides ... (simplified for brevity)
            int grid_idx = type_idx * n_isorp_max; // + shell mapping
            float hx = h_grids_3c[grid_idx].x;
            float hy = h_grids_3c[grid_idx].y;
            int nx_u = dims_3c[grid_idx].x;
            int ny_u = dims_3c[grid_idx].y;
            
            int trip_str = n_isorp_max * 5 * ny_max * nx_max * n_nz_max;
            int th_str = ny_max * nx_max * n_nz_max;
            
            for(int c=0; c<5; c++){
                const __global float* sd = &data_3c[type_idx*trip_str + c*th_str];
                float val = interpolate_2d(x, y, nx_u, ny_u, hx, hy, sd, 0, 1);
                hlist[c] = val * P[c];
            }
            recover_and_rotate_3c_sp(hlist, eps, v_short);
        } else {
            for(int i=0; i<16; i++) v_short[i] = 0.0f;
        }

        // --- Long Range Term (Multipole/Ewald) ---
        // V_long = (S_ij / 2) * (1/dik + 1/djk) (Simplified monopole)
        // Fireball includes dipole terms here too, effectively (S/R - D/R^2).
        // Using simplified monopole here:
        float inv_dik = (d_ik > 1e-4f) ? 1.0f/d_ik : 0.0f;
        float inv_djk = (d_jk > 1e-4f) ? 1.0f/d_jk : 0.0f;
        float geo_factor = 0.5f * (inv_dik + inv_djk);

        // Accumulate
        for(int i=0; i<16; i++) {
            float s_val = S_ptr[i];
            float v_long = s_val * geo_factor; // Monopole term
            
            // Final combination: DeltaQ * [ stn*Vshort + (1-stn)*Vlong ]
            float term = dq_k * (stn * v_short[i] + one_minus_stn * v_long);
            block_acc[i] += term;
        }
    }

    // Write to H_3c output
    // Note: This needs to be ADDED to Neutral Atom 3c if they share a buffer, 
    // or kept separate. Here we assume a dedicated buffer or pre-filled buffer.
    __global float* out = &H_3c_out[pid * 16];
    for(int i=0; i<16; i++) {
        // Atomic add usually not needed if we parallelize by pair and this is the only kernel writing to this pair.
        // If NA terms are separate, we just write.
        out[i] += block_acc[i]; 
    }

    return;

#if 0
            float hy = h_grids_3c[grid_idx].y;
            int nx = dims_3c[grid_idx].x;
            int ny = dims_3c[grid_idx].y;

            int trip_stride  = n_isorp_max * 5 * ny_max * nx_max * n_nz_max;
            int isorp_stride = 5 * ny_max * nx_max * n_nz_max;
            int theta_stride = ny_max * nx_max * n_nz_max;
            const __global float* base_trip = &data_3c[type_idx * trip_stride + isorp_idx * isorp_stride];

            for(int iME=0; iME<5; iME++){
                float sum = 0.0f;
                for(int it=0; it<5; it++){
                    const __global float* plane = base_trip + it * theta_stride;
                    sum += P[it] * interpolate_2d(x, d_ij, nx, ny, hx, hy, plane, iME, n_nz_max);
                }
                hlist[iME] = sum;
            }

            recover_and_rotate_3c_sp(hlist, eps, v_short);
        } else {
            for(int m=0; m<16; m++) v_short[m] = 0.0f;
        }

        // --- 3. Long Range (Ewald / Multipole) ---
        // V_lr = 0.5 * S_ij * (1/Rik + 1/Rjk)
        float inv_R = 0.5f * ((d_ik > 1e-3f ? 1.0f/d_ik : 0.0f) + (d_jk > 1e-3f ? 1.0f/d_jk : 0.0f));
        
        // --- 4. Mix and Accumulate ---
        for(int m=0; m<16; m++) {
            float v_long = s_ptr[m] * inv_R;
            float term = stn * v_short[m] + (1.0f - stn) * v_long;
            block_acc[m] += term * dq;
        }
    }

    // Write output
    __global float* out = &V_ca_out[pid*16];
    for(int m=0; m<16; m++) out[m] = block_acc[m];
#endif
}


// ================================================================================================
// Pseudopotential 3C (assemble_3c_PP)
// ================================================================================================

__kernel void assemble_3c_PP(
    const int n_pairs,
    // Lists
    const __global int4* pairs,       // [i, j, type_idx, _]
    // Since K must overlap with I and J, we usually iterate all atoms or use a neighbor list.
    // Efficient approach: iterate K in neighbor list of I, check if K is neighbor of J.
    // OR: Use the same Common Neighbor list as above.
    const __global int* cn_offsets,
    const __global int* cn_indices,
    // Data
    const __global float* sVNL,       // Sparse [n_sVNL_pairs, 16] (basis x projector)
    const __global int* sVNL_map,     // Map (atom_basis, atom_proj) -> index in sVNL
    // sVNL_map is usually too big (NxN). 
    // Alternative: pass sVNL indices in the CN list if pre-calculated in Python.
    // Assuming here we have a way to fetch S(i,k) and S(j,k).
    // Let's assume sVNL is stored linearly per neighbor pair in a separate list 
    // and we passed indices in `cn_indices` (packed struct).
    // FOR SIMPLICITY: Using the previously provided 'contract_vnl_sumk_sp' logic 
    // where we explicitly pass indices of IK and JK pairs.
    const __global int* idx_ik_list,  // [n_pairs * n_k_max] indices into sVNL
    const __global int* idx_jk_list,  // [n_pairs * n_k_max] indices into sVNL
    const int n_k_max,                // Max common neighbors (padding)
    const __global float* cl_coeffs,  // [natoms, 4]
    __global float* V_nl_out
) {
    int pid = get_global_id(0);
    if (pid >= n_pairs) return;

    float acc[16];
    for(int m=0; m<16; m++) acc[m] = 0.0f;

    // Loop over possible K atoms (columns in the index lists)
    for(int k_iter = 0; k_iter < n_k_max; k_iter++) {
        int idx_ik = idx_ik_list[pid * n_k_max + k_iter];
        int idx_jk = idx_jk_list[pid * n_k_max + k_iter];

        if (idx_ik < 0 || idx_jk < 0) continue; // No interaction

        // We also need the atom index of K to fetch Coefficients
        // This should be stored alongside indices or derived.
        // Assuming we look it up from the common neighbor list
        int atom_k = cn_indices[cn_offsets[pid] + k_iter]; 
        
        const __global float* C = &cl_coeffs[atom_k * 4];
        const __global float* Sik = &sVNL[idx_ik * 16];
        const __global float* Sjk = &sVNL[idx_jk * 16];

        // Contraction: Sum_n (Sik[mu, n] * C[n] * Sjk[nu, n])
        // Note: sVNL is typically (basis, proj) or (proj, basis). 
        // Fireball: <phi | alpha>.
        // Let's assume sVNL blocks are row=basis, col=projector.
        
        for(int mu=0; mu<4; mu++) {
            for(int nu=0; nu<4; nu++) {
                float val = 0.0f;
                for(int n=0; n<4; n++) {
                    // C[n] is diagonal epsilon
                    // Sjk[nu, n] is <phi_j_nu | alpha_k_n>
                    // Sik[mu, n] is <phi_i_mu | alpha_k_n>
                    val += Sik[mu*4 + n] * C[n] * Sjk[nu*4 + n];
                }
                acc[mu*4 + nu] += val;
            }
        }
    }

    __global float* out = &V_nl_out[pid * 16];
    for(int m=0; m<16; m++) out[m] = acc[m];
}