
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