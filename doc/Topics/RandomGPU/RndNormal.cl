// =============================================================================
// RndNormal.cl  –  GPU normal-distribution random number generators
// -----------------------------------------------------------------------------
// Required defines (injected from host via -D or prepended string):
//   BUFFER_SIZE   – number of floats in the pre-generated normal buffer
//   LDS_SIZE      – local work-group size for shared-memory kernel (power of 2)
//   INV_SQRT2     – 0.7071067811865475f
// =============================================================================

// -----------------------------------------------------------------------------
// PRNG helpers
// -----------------------------------------------------------------------------

// Wang hash: fast 32-bit avalanche, good for seeding per-thread state
inline uint wang_hash(uint x) {
    x = (x ^ 61u) ^ (x >> 16u);
    x *= 9u;
    x ^= x >> 4u;
    x *= 0x27d4eb2du;
    x ^= x >> 15u;
    return x;
}

// Xorshift32: cheap 32-bit state advance (period 2^32-1)
inline uint xorshift32(uint s) {
    s ^= s << 13u; s ^= s >> 17u; s ^= s << 5u;
    return s;
}

// SplitMix64: counter-based 64-bit hash, excellent avalanche, no state
inline ulong splitmix64(ulong x) {
    x += 0x9e3779b97f4a7c15UL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9UL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebUL;
    return x ^ (x >> 31);
}

// Convert upper 24 bits of a 32-bit word to float in (0, 1]
// NOTE: guarded with max() so log() never sees 0
inline float uint_to_f01(uint x) {
    return max((float)(x >> 8u) * (1.0f / 16777216.0f), 1e-7f);
}

// Convert upper 24 bits of a 64-bit word to float in (0, 1]
inline float ulong_to_f01(ulong x) {
    return max((float)(x >> 40) * (1.0f / 16777216.0f), 1e-7f);
}

// -----------------------------------------------------------------------------
// 1) Box-Muller  –  32-bit PRNG (wang_hash + xorshift32)
//    Fastest on hardware with weak 64-bit ALU (e.g. older NVIDIA).
//    Uses native_sqrt/log/cos for maximum throughput.
// -----------------------------------------------------------------------------
__kernel void box_muller_32(__global float* out, uint seed) {
    uint gid   = get_global_id(0);
    uint state = wang_hash(gid + seed);
    uint s2    = xorshift32(state);
    float u1   = uint_to_f01(state);
    float u2   = (float)(s2 >> 8u) * (1.0f / 16777216.0f);    // [0,1) ok for theta
    float r    = native_sqrt(-2.0f * native_log(u1));
    out[gid]   = r * native_cos(2.0f * M_PI_F * u2);
}

// -----------------------------------------------------------------------------
// 2) Box-Muller  –  64-bit PRNG (splitmix64)
//    Better statistical quality (two fully independent hashes per sample).
//    Preferred when 64-bit throughput is acceptable.
// -----------------------------------------------------------------------------
__kernel void box_muller_64(__global float* out, uint seed) {
    uint  gid = get_global_id(0);
    ulong h1  = splitmix64((ulong)(gid * 2u + seed));
    ulong h2  = splitmix64((ulong)(gid * 2u + seed + 1u));
    float u1  = ulong_to_f01(h1);
    float u2  = (float)(h2 >> 40) * (1.0f / 16777216.0f);
    float r   = native_sqrt(-2.0f * native_log(u1));
    out[gid]  = r * native_cos(2.0f * M_PI_F * u2);
}

// -----------------------------------------------------------------------------
// 3) Inverse-CDF (probit) via Beasley-Springer-Moro rational approximation
//    Branch-free central region + log-log tail. No rejection, deterministic.
//    Good quality, moderate cost (rational poly + 2x native_log in tail).
// -----------------------------------------------------------------------------
__kernel void inverse_cdf(__global float* out, uint seed) {
    uint  gid = get_global_id(0);
    float u   = ulong_to_f01(splitmix64((ulong)(gid + seed)));
    u = min(u, 1.0f - 1e-7f);
    float p   = u - 0.5f;
    float r, x;
    if (fabs(p) < 0.42f) {
        r = p * p;
        x = p * ((((-25.44106049637f * r + 41.39119773534f) * r - 18.61500062529f) * r + 2.50662823884f)
              / ((((3.13082909833f * r - 21.06224101826f) * r + 23.08336743743f) * r - 8.47351093090f) * r + 1.0f));
    } else {
        r = (p > 0.0f) ? (1.0f - u) : u;
        r = native_log(-native_log(r));
        x = 0.3374754822726869f + r*(0.9761690190917186f + r*(0.1607979714918209f
          + r*(0.0276438810333863f + r*(0.0038405729373609f + r*(0.0003951896511349f
          + r*(0.0000321767881768f + r*(0.0000002888167364f + r*0.0000003960315187f)))))));
        if (p < 0.0f) x = -x;
    }
    out[gid] = x;
}

// -----------------------------------------------------------------------------
// 4) Marsaglia polar  –  rejection sampling
//    NOTE: GPU threads diverge on rejection → slower than Box-Muller on GPU.
//    Included for comparison / completeness. Iteration capped at 32 to
//    prevent GPU hangs (probability of hitting cap: ~(1-pi/4)^32 ≈ 1e-8).
// -----------------------------------------------------------------------------
__kernel void marsaglia_polar(__global float* out, uint seed) {
    uint  gid  = get_global_id(0);
    uint  i    = gid;
    float u1, u2, s;
    int   iter = 0;
    do {
        ulong h1 = splitmix64((ulong)(i * 2u + seed));
        ulong h2 = splitmix64((ulong)(i * 2u + seed + 1u));
        u1 = 2.0f * (float)(h1 >> 40) * (1.0f / 16777216.0f) - 1.0f;
        u2 = 2.0f * (float)(h2 >> 40) * (1.0f / 16777216.0f) - 1.0f;
        s  = u1*u1 + u2*u2;
        i++; iter++;
    } while ((s >= 1.0f || s == 0.0f) && iter < 32);
    if (s <= 0.0f) s = 1e-7f;
    out[gid] = u1 * native_sqrt(-2.0f * native_log(s) / s);
}

// -----------------------------------------------------------------------------
// 5) Global-memory buffer  –  random lookup into pre-generated normal pool
//    Bottleneck: random (non-coalesced) global memory reads → high latency.
//    Index uses power-of-2 mask when BUFFER_SIZE is power-of-2 (fast),
//    otherwise modulo (slow). Host should pass a power-of-2 BUFFER_SIZE.
//    Two independent reads + average preserves N(0,1) (sum of two N(0,1)/sqrt(2)).
// -----------------------------------------------------------------------------
__kernel void global_buffer(__global float* out, __global const float* buf, uint seed) {
    uint gid  = get_global_id(0);
    uint idx1 = (uint)(splitmix64((ulong)(gid + seed))              % BUFFER_SIZE);
    uint idx2 = (uint)(splitmix64((ulong)(gid + seed + 0x9e3779b9u)) % BUFFER_SIZE);
    out[gid]  = (buf[idx1] + buf[idx2]) * INV_SQRT2;
}

// -----------------------------------------------------------------------------
// 6) Shared-memory buffer  –  cooperative load into LDS, then random lookup
//    Each workgroup picks a random contiguous block from the global buffer
//    (coalesced reads → high bandwidth), loads it to fast LDS, then each
//    thread draws two random indices within that block.
//    LDS_SIZE should be >= local_work_size and a power of 2.
//    local_work_size must equal LDS_SIZE when launching this kernel.
// -----------------------------------------------------------------------------
__kernel void shared_buffer(__global float* out, __global const float* buf, uint seed) {
    uint lid        = get_local_id(0);
    uint gid        = get_global_id(0);
    uint group_id   = get_group_id(0);
    uint local_size = get_local_size(0);   // == LDS_SIZE at launch

    __local float lds[LDS_SIZE];

    // Whole workgroup picks a random starting offset (coalesced load)
    uint group_state = wang_hash(group_id + seed);
    uint offset      = group_state % (BUFFER_SIZE - LDS_SIZE);

    // Cooperatively load LDS_SIZE floats; each thread loads (LDS_SIZE/local_size) elements
    for (uint i = lid; i < LDS_SIZE; i += local_size)
        lds[i] = buf[offset + i];
    barrier(CLK_LOCAL_MEM_FENCE);

    // Per-thread independent random indices within LDS (power-of-2 mask if possible)
    uint ts    = wang_hash(gid + seed + 0x12345u);
    ts         = xorshift32(ts);
    uint idx1  = ts % LDS_SIZE;
    ts         = xorshift32(ts);
    uint idx2  = ts % LDS_SIZE;

    out[gid] = (lds[idx1] + lds[idx2]) * INV_SQRT2;
}
