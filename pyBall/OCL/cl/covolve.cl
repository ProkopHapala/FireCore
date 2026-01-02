// Defaults
#ifndef LS_X
#define LS_X 8
#endif
#ifndef LS_Y
#define LS_Y 8
#endif
#ifndef LS_Z
#define LS_Z 8
#endif

#define TW (LS_X + 2)
#define TH (LS_Y + 2)
#define TD (LS_Z + 2)

inline int pmod(int i, int n) {
    int r = i % n;
    return (r < 0) ? r + n : r;
}

__kernel void Convolution3D_General(
    const int4 ns,            // {nx, ny, nz, 0}
    __global const float* Gs, // Input
    __global const float* G0, // Optional Additive Input
    __global       float* out,// Output
    const float4   weights,   // x=Center, y=Face, z=Edge, w=Corner
    const int4     bPBC,      // Boundary: 1=Periodic, 0=Zero
    const float2   coefs      // x=Multiplicative, y=Additive G0 scale
) {
    const int lx = get_local_id(0);
    const int ly = get_local_id(1);
    const int lz = get_local_id(2);
    
    const int g_start_x = get_group_id(0) * LS_X;
    const int g_start_y = get_group_id(1) * LS_Y;
    const int g_start_z = get_group_id(2) * LS_Z;

    __local float tile[TD][TH][TW];

    const int tid = lz * (LS_Y * LS_X) + ly * LS_X + lx; 
    const int group_nthreads = LS_X * LS_Y * LS_Z;
    const int tile_nelements = TW * TH * TD;

    // --- 1. LOAD TILE ---
    for (int i = tid; i < tile_nelements; i += group_nthreads) {
        int r = i;
        int tx = r % TW; r /= TW;
        int ty = r % TH; r /= TH;
        int tz = r;

        int gx = g_start_x + tx - 1;
        int gy = g_start_y + ty - 1;
        int gz = g_start_z + tz - 1;

        bool inside = true;
        if (gx < 0 || gx >= ns.x) { if (bPBC.x) gx = pmod(gx, ns.x); else inside = false; }
        if (gy < 0 || gy >= ns.y) { if (bPBC.y) gy = pmod(gy, ns.y); else inside = false; }
        if (gz < 0 || gz >= ns.z) { if (bPBC.z) gz = pmod(gz, ns.z); else inside = false; }

        float val = 0.0f;
        if (gx >= 0 && gx < ns.x && gy >= 0 && gy < ns.y && gz >= 0 && gz < ns.z) {
            val = Gs[gz * (ns.x * ns.y) + gy * ns.x + gx];
        }
        
        tile[tz][ty][tx] = val;
    }

    barrier(CLK_LOCAL_MEM_FENCE);

    // --- 2. CHECK OUTPUT BOUNDS ---
    int out_x = g_start_x + lx;
    int out_y = g_start_y + ly;
    int out_z = g_start_z + lz;

    if (out_x >= ns.x || out_y >= ns.y || out_z >= ns.z) return;

    // --- 3. CONVOLUTION ---
    int tx = lx + 1;
    int ty = ly + 1;
    int tz = lz + 1;

    float sum = 0.0f;

    // Center
    sum += tile[tz][ty][tx] * weights.x;

    // Faces (6)
    sum += (tile[tz][ty][tx-1] + tile[tz][ty][tx+1] +
            tile[tz][ty-1][tx] + tile[tz][ty+1][tx] +
            tile[tz-1][ty][tx] + tile[tz+1][ty][tx]) * weights.y;

    // Edges (12)
    sum += (tile[tz][ty-1][tx-1] + tile[tz][ty-1][tx+1] +
            tile[tz][ty+1][tx-1] + tile[tz][ty+1][tx+1] +
            tile[tz-1][ty][tx-1] + tile[tz-1][ty][tx+1] +
            tile[tz+1][ty][tx-1] + tile[tz+1][ty][tx+1] +
            tile[tz-1][ty-1][tx] + tile[tz-1][ty+1][tx] +
            tile[tz+1][ty-1][tx] + tile[tz+1][ty+1][tx]) * weights.z;

    // Corners (8)
    sum += (tile[tz-1][ty-1][tx-1] + tile[tz-1][ty-1][tx+1] +
            tile[tz-1][ty+1][tx-1] + tile[tz-1][ty+1][tx+1] +
            tile[tz+1][ty-1][tx-1] + tile[tz+1][ty-1][tx+1] +
            tile[tz+1][ty+1][tx-1] + tile[tz+1][ty+1][tx+1]) * weights.w;

    // --- 4. WRITE ---
    int g_idx = out_z * (ns.x * ns.y) + out_y * ns.x + out_x;
    
    sum *= coefs.x;
    if (G0) { sum += G0[g_idx] * coefs.y; }

    out[g_idx] = sum;
}