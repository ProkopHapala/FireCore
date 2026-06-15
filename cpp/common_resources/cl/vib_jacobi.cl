// Vibration frequency-domain Jacobi solver: A(omega) u = f
// Block Jacobi with heavy-ball momentum; one workgroup per cluster/molecule.
// Pattern ported from LFF.cl lff_jacobi (local gather, private neighbor cache).

#define VIB_WG_SIZE     32
#define VIB_MAX_NEIGH   32

inline float2 cmul(float2 a, float2 b) { return (float2)(a.x*b.x - a.y*b.y, a.x*b.y + a.y*b.x); }
inline float2 cadd(float2 a, float2 b) { return a + b; }

inline void cmat3_vec_reim(
    __global const float* Bre, __global const float* Bim,
    float3 vre, float3 vim, float* out_re, float* out_im
){
    for (int r = 0; r < 3; ++r) {
        float2 acc = (float2)(0.f, 0.f);
        for (int c = 0; c < 3; ++c) {
            float2 vc = (float2)(vre[c], vim[c]);
            float2 br = (float2)(Bre[r*3+c], Bim[r*3+c]);
            acc = cadd(acc, cmul(br, vc));
        }
        out_re[r] = acc.x; out_im[r] = acc.y;
    }
}

__kernel void vib_jacobi_solve(
    __global const int*    row_count,     // [natoms]
    __global const int*    row_neigh,     // [natoms, VIB_MAX_NEIGH]
    __global const float*  row_blk_re,    // [natoms, VIB_MAX_NEIGH, 9]
    __global const float*  row_blk_im,    // [natoms, VIB_MAX_NEIGH, 9]
    __global const float*  diag_inv_re,   // [natoms, 9]
    __global const float*  diag_inv_im,   // [natoms, 9]
    __global const float*  rhs_re,        // [natoms, 3]
    __global const float*  rhs_im,        // [natoms, 3]
    __global       float*  u_re,          // [natoms, 3] in/out
    __global       float*  u_im,          // [natoms, 3] in/out
    const int              natoms,
    const int              n_iter,
    const float            b_mix
){
    const int ia0   = 0; // single cluster per launch
    const int lid   = get_local_id(0);
    const int lsz   = get_local_size(0);
    if (lid >= natoms) return;

    const int idx = ia0 + lid;
    const int neigh_base = idx * VIB_MAX_NEIGH;
    const int blk_stride = VIB_MAX_NEIGH * 9;

    __local float lu_re[VIB_WG_SIZE][3];
    __local float lu_im[VIB_WG_SIZE][3];

    int nneigh = row_count[idx];
    int ng_idx[VIB_MAX_NEIGH];
    for (int k = 0; k < VIB_MAX_NEIGH; ++k) {
        ng_idx[k] = row_neigh[neigh_base + k];
        if (k >= nneigh) ng_idx[k] = -1;
    }

    float3 u_old_re = (float3)(u_re[idx*3+0], u_re[idx*3+1], u_re[idx*3+2]);
    float3 u_old_im = (float3)(u_im[idx*3+0], u_im[idx*3+1], u_im[idx*3+2]);
    float3 pi_re = u_old_re;
    float3 pi_im = u_old_im;

    for (int iter = 0; iter < n_iter; ++iter) {
        lu_re[lid][0] = pi_re.x; lu_re[lid][1] = pi_re.y; lu_re[lid][2] = pi_re.z;
        lu_im[lid][0] = pi_im.x; lu_im[lid][1] = pi_im.y; lu_im[lid][2] = pi_im.z;
        barrier(CLK_LOCAL_MEM_FENCE);

        float acc_re[3] = {rhs_re[idx*3+0], rhs_re[idx*3+1], rhs_re[idx*3+2]};
        float acc_im[3] = {rhs_im[idx*3+0], rhs_im[idx*3+1], rhs_im[idx*3+2]};

        for (int k = 0; k < nneigh; ++k) {
            int j = ng_idx[k];
            if (j < 0) break;
            const int boff = idx * blk_stride + k * 9;
            float sub_re[3], sub_im[3];
            cmat3_vec_reim(row_blk_re + boff, row_blk_im + boff,
                (float3)(lu_re[j][0], lu_re[j][1], lu_re[j][2]),
                (float3)(lu_im[j][0], lu_im[j][1], lu_im[j][2]),
                sub_re, sub_im);
            for (int r = 0; r < 3; ++r) { acc_re[r] -= sub_re[r]; acc_im[r] -= sub_im[r]; }
        }

        const int dioff = idx * 9;
        float sol_re[3], sol_im[3];
        cmat3_vec_reim(diag_inv_re + dioff, diag_inv_im + dioff,
            (float3)(acc_re[0], acc_re[1], acc_re[2]),
            (float3)(acc_im[0], acc_im[1], acc_im[2]),
            sol_re, sol_im);

        float3 pi_new_re = (float3)(sol_re[0], sol_re[1], sol_re[2]);
        float3 pi_new_im = (float3)(sol_im[0], sol_im[1], sol_im[2]);

        // heavy-ball momentum: u = u_new + b_mix * (u - u_old)
        float3 du_re = pi_new_re - u_old_re;
        float3 du_im = pi_new_im - u_old_im;
        pi_re = pi_new_re + b_mix * du_re;
        pi_im = pi_new_im + b_mix * du_im;
        u_old_re = pi_re; u_old_im = pi_im;

        barrier(CLK_LOCAL_MEM_FENCE);
    }

    u_re[idx*3+0] = pi_re.x; u_re[idx*3+1] = pi_re.y; u_re[idx*3+2] = pi_re.z;
    u_im[idx*3+0] = pi_im.x; u_im[idx*3+1] = pi_im.y; u_im[idx*3+2] = pi_im.z;
}
