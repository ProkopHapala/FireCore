
// From Google Gemini: 
//  https://aistudio.google.com/app/prompts?state=%7B%22ids%22:%5B%221iMLknKwbiMjTYf1iMgVCodYwCWJuVfaP%22%5D,%22action%22:%22open%22,%22userId%22:%22100958146796876347936%22,%22resourceKeys%22:%7B%7D%7D&usp=sharing



// We assume support for float atomic add is available, either via
// cl_khr_fp_atomic extension or OpenCL 2.0+.
// If not available, a cmpxchg version for float must be used.
#pragma OPENCL EXTENSION cl_khr_fp_atomic : enable // Explicitly request float atomics

// Define Coulomb constant (adjust value and units as needed)
#define K_COULOMB 138.935458f // Use 'f' suffix for float literals

//----------------------------------------------------------------------
// Helper: Atomic add for float3 (component-wise)
// Assumes atomic_addf is available
//----------------------------------------------------------------------
void atomic_add_float3(__global float *addr, float3 val) {
    atomic_addf(&addr[0], val.x);
    atomic_addf(&addr[1], val.y);
    atomic_addf(&addr[2], val.z);
}

// If atomic_addf is NOT supported on your target device/OpenCL version,
// use this alternative cmpxchg version instead:
/*
inline void atomic_add_float(__global float *addr, float val) {
    union {
        uint u;
        float f;
    } old_val, new_val;
    do {
        old_val.f = *addr;
        new_val.f = old_val.f + val;
    } while (atom_cmpxchg((__global uint *)addr, old_val.u, new_val.u) != old_val.u);
}
void atomic_add_float3(__global float *addr, float3 val) {
    atomic_add_float(&addr[0], val.x);
    atomic_add_float(&addr[1], val.y);
    atomic_add_float(&addr[2], val.z);
}
*/


//======================================================================
// Interaction Force Calculation Functions (Inline)
// These calculate the force contribution ON bucket B due to bucket A.
// F_b = Grad_{R_ab}(U_ab)
// Inputs:
//   M_a, M_b: Packed float16 moments for buckets A and B
//   Rhat: Normalized vector R_b - R_a
//   RinvN: Appropriate power of 1/|R_b - R_a| (e.g., Rinv2 = 1/R^2)
//======================================================================

//----------------------------------------------------------------------
// Extract Dipole Vector from packed data
//----------------------------------------------------------------------
inline float3 get_dipole(float16 M) {
    return (float3)(M.s4, M.s5, M.s6);
}

//----------------------------------------------------------------------
// Calculate Quadrupole Matrix (Theta) * Vector (v) product
//----------------------------------------------------------------------
inline float3 quadrupole_mul_vec(float16 M_quad, float3 v) {
    // M_quad holds moments of the bucket defining the quadrupole tensor Theta
    // Access Theta components: s7=xx, s8=xy, s9=xz, sa=yx, sb=yy, sc=yz, sd=zx, se=zy, sf=zz
    return (float3)(
        M_quad.s7 * v.x + M_quad.s8 * v.y + M_quad.s9 * v.z, // Row 1: xx, xy, xz
        M_quad.sa * v.x + M_quad.sb * v.y + M_quad.sc * v.z, // Row 2: yx, yy, yz
        M_quad.sd * v.x + M_quad.se * v.y + M_quad.sf * v.z  // Row 3: zx, zy, zz
    );
    // Assumes host stored symmetric tensor correctly (e.g., M_quad.s8 == M_quad.sa)
}

//----------------------------------------------------------------------
// Calculate Rhat . Theta . Rhat = Rhat^T * Theta * Rhat
//----------------------------------------------------------------------
inline float rhat_dot_theta_dot_rhat(float16 M_quad, float3 Rhat) {
    float3 TRhat = quadrupole_mul_vec(M_quad, Rhat);
    return dot(Rhat, TRhat);
}

//----------------------------------------------------------------------
// p=0: Monopole-Monopole (MM) Force
// U = Q_a * Q_b / R
// Grad(U) = - Q_a * Q_b * Rhat / R^2
//----------------------------------------------------------------------
inline float3 calc_force_MM(
    float16 M_a, float16 M_b,
    float3 Rhat, float Rinv2)
{
    float Q_a = M_a.s3;
    float Q_b = M_b.s3;
    return -Q_a * Q_b * Rinv2 * Rhat;
}

//----------------------------------------------------------------------
// p=1: Monopole(a)-Dipole(b) (MD) Force
// U = Q_a * (D_b . Rhat) / R^2
// Grad(U) = Q_a * [ D_b - 3 * (D_b.Rhat) * Rhat ] / R^3
//----------------------------------------------------------------------
inline float3 calc_force_MD(
    float16 M_a, float16 M_b,
    float3 Rhat, float Rinv3)
{
    float Q_a = M_a.s3;
    float3 D_b = get_dipole(M_b);
    float D_b_dot_Rhat = dot(D_b, Rhat);
    return Q_a * Rinv3 * (D_b - 3.0f * D_b_dot_Rhat * Rhat);
}

//----------------------------------------------------------------------
// p=1: Dipole(a)-Monopole(b) (DM) Force
// U = - Q_b * (D_a . Rhat) / R^2
// Grad(U) = -Q_b * [ D_a - 3 * (D_a.Rhat) * Rhat ] / R^3
//----------------------------------------------------------------------
inline float3 calc_force_DM(
    float16 M_a, float16 M_b,
    float3 Rhat, float Rinv3)
{
    float3 D_a = get_dipole(M_a);
    float Q_b = M_b.s3;
    float D_a_dot_Rhat = dot(D_a, Rhat);
    return -Q_b * Rinv3 * (D_a - 3.0f * D_a_dot_Rhat * Rhat);
}

//----------------------------------------------------------------------
// p=1: Dipole-Dipole (DD) Force
// U = [ D_a.D_b - 3 * (D_a.Rhat) * (D_b.Rhat) ] / R^3
// Grad(U) = -3/R^4 * [ (D_a.R)D_b + (D_b.R)D_a + (D_a.D_b)Rhat - 5(D_a.R)(D_b.R)Rhat ]
//----------------------------------------------------------------------
inline float3 calc_force_DD(
    float16 M_a, float16 M_b,
    float3 Rhat, float Rinv4)
{
    float3 D_a = get_dipole(M_a);
    float3 D_b = get_dipole(M_b);
    float D_a_dot_Rhat = dot(D_a, Rhat);
    float D_b_dot_Rhat = dot(D_b, Rhat);
    float D_a_dot_D_b = dot(D_a, D_b);
    float3 term1 = D_a_dot_Rhat * D_b;
    float3 term2 = D_b_dot_Rhat * D_a;
    float3 term3 = D_a_dot_D_b * Rhat;
    float3 term4 = -5.0f * D_a_dot_Rhat * D_b_dot_Rhat * Rhat;
    return -3.0f * Rinv4 * (term1 + term2 + term3 + term4);
}

//----------------------------------------------------------------------
// p=2: Monopole(a)-Quadrupole(b) (MQ) Force
// U = Q_a * (Rhat . Theta_b . Rhat) / (2 * R^3)
// Grad(U) = (Q_a / 2R^4) * [ 2 * Theta_b . Rhat - 3 * (Rhat.Theta_b.Rhat) * Rhat ]
//----------------------------------------------------------------------
inline float3 calc_force_MQ(
    float16 M_a, float16 M_b,
    float3 Rhat, float Rinv4)
{
    float Q_a = M_a.s3;
    float3 Thet_b_Rhat = quadrupole_mul_vec(M_b, Rhat);
    float Rhat_Thet_b_Rhat = dot(Rhat, Thet_b_Rhat); // = Rhat . Theta_b . Rhat
    return (Q_a * 0.5f * Rinv4) * (2.0f * Thet_b_Rhat - 3.0f * Rhat_Thet_b_Rhat * Rhat);
}

//----------------------------------------------------------------------
// p=2: Quadrupole(a)-Monopole(b) (QM) Force
// U = Q_b * (Rhat . Theta_a . Rhat) / (2 * R^3)
// Grad(U) = (Q_b / 2R^4) * [ 2 * Theta_a . Rhat - 3 * (Rhat.Theta_a.Rhat) * Rhat ]
//----------------------------------------------------------------------
inline float3 calc_force_QM(
    float16 M_a, float16 M_b,
    float3 Rhat, float Rinv4)
{
    float Q_b = M_b.s3;
    float3 Thet_a_Rhat = quadrupole_mul_vec(M_a, Rhat);
    float Rhat_Thet_a_Rhat = dot(Rhat, Thet_a_Rhat); // = Rhat . Theta_a . Rhat
     return (Q_b * 0.5f * Rinv4) * (2.0f * Thet_a_Rhat - 3.0f * Rhat_Thet_a_Rhat * Rhat);
}

//----------------------------------------------------------------------
// p=2: Dipole(a)-Quadrupole(b) (DQ) Force
// ** NOTE: Formulas for DQ/QD/QQ forces are complex and require careful validation. **
// ** The version here is based on common forms but should be checked against literature. **
//----------------------------------------------------------------------
inline float3 calc_force_DQ(
    float16 M_a, float16 M_b,
    float3 Rhat, float Rinv5)
{
    float3 D_a = get_dipole(M_a);
    float D_a_dot_Rhat = dot(D_a, Rhat);

    float3 Thet_b_Rhat = quadrupole_mul_vec(M_b, Rhat);
    float Rhat_Thet_b_Rhat = dot(Rhat, Thet_b_Rhat); // R.T_b.R
    float D_a_dot_Thet_b_Rhat = dot(D_a, Thet_b_Rhat); // D_a.T_b.R

    // Based on gradient of U_DQ ~ Rinv4 (potential)
    // Force F_b = Grad_{R_ab}(U_DQ) ~ Rinv5
    // Using the simplified form from previous answer (requires validation):
    float3 force = Rinv5 * (-2.0f * D_a_dot_Thet_b_Rhat * Rhat
                             + 5.0f * D_a_dot_Rhat * Rhat_Thet_b_Rhat * Rhat // Factor 5 might be different based on potential definition
                             - D_a * Rhat_Thet_b_Rhat
                             - 2.0f * D_a_dot_Rhat * Thet_b_Rhat
                            );
    return force;
}

//----------------------------------------------------------------------
// p=2: Quadrupole(a)-Dipole(b) (QD) Force
// ** See NOTE for DQ **
//----------------------------------------------------------------------
inline float3 calc_force_QD(
    float16 M_a, float16 M_b,
    float3 Rhat, float Rinv5)
{
    float3 D_b = get_dipole(M_b);
    float D_b_dot_Rhat = dot(D_b, Rhat);

    float3 Thet_a_Rhat = quadrupole_mul_vec(M_a, Rhat);
    float Rhat_Thet_a_Rhat = dot(Rhat, Thet_a_Rhat); // R.T_a.R
    float D_b_dot_Thet_a_Rhat = dot(D_b, Thet_a_Rhat); // D_b.T_a.R

    // Based on gradient of U_QD ~ Rinv4 (potential)
    // Force F_b = Grad_{R_ab}(U_QD) ~ Rinv5
    // Using symmetric form based on DQ (requires validation):
    float3 force = Rinv5 * (+2.0f * D_b_dot_Thet_a_Rhat * Rhat
                             - 5.0f * D_b_dot_Rhat * Rhat_Thet_a_Rhat * Rhat // Factor 5 might be different
                             + D_b * Rhat_Thet_a_Rhat
                             + 2.0f * D_b_dot_Rhat * Thet_a_Rhat
                            );
    return force;
}

//----------------------------------------------------------------------
// p=2: Quadrupole-Quadrupole (QQ) Force
// ** Omitted due to extreme complexity of the gradient formula. **
//----------------------------------------------------------------------
// inline float3 calc_force_QQ(...) { ... return (float3)(0.0f); }


//======================================================================
// Main M2M Force Kernel (Packed Float Data Version)
//======================================================================
__kernel void compute_m2m_force_packed_float(
    // --- Input Buffers ---
    __global const int2* restrict bucket_pairs,          // Array of (idx_a, idx_b) pairs to compute
    __global const float16* restrict bucket_moments,    // Packed moments (pos, Q, D, Theta) for each bucket
    // Particle Data (needed for force distribution)
    __global const float* restrict particle_charges,     // Charge of each particle
    __global const int* restrict bucket_particle_indices,// Global index of particles per bucket (flattened/offset)
    __global const int* restrict bucket_particle_counts, // Number of particles per bucket
    __global const int* restrict bucket_particle_offsets,// Starting offset for each bucket in `bucket_particle_indices`
    // Accuracy control
    const int expansion_order, // 0:MM, 1:up to DD, 2:up to QD (QQ omitted)
    // --- Output Buffer ---
    __global float* restrict particle_forces            // Global buffer (Nx3 floats) for forces (accumulates!)
) {
    int gid = get_global_id(0);

    // Get bucket indices for this work item
    int idx_a = bucket_pairs[gid].x;
    int idx_b = bucket_pairs[gid].y;

    // --- Load Packed Moments ---
    float16 M_a = bucket_moments[idx_a];
    float16 M_b = bucket_moments[idx_b];

    // --- Geometric Factors ---
    float3 R_a = (float3)(M_a.s0, M_a.s1, M_a.s2);
    float3 R_b = (float3)(M_b.s0, M_b.s1, M_b.s2);
    float3 R_ab_vec = R_b - R_a;
    float R2 = dot(R_ab_vec, R_ab_vec);

    // Use a suitable small float threshold for singularity
    if (R2 < 1.0e-6f) return;

    float R = sqrt(R2);
    float Rinv = 1.0f / R;
    float Rinv2 = Rinv * Rinv;
    float Rinv3 = Rinv2 * Rinv;
    float Rinv4 = Rinv3 * Rinv;
    float Rinv5 = Rinv4 * Rinv;
    // float Rinv6 = Rinv5 * Rinv; // If QQ were included

    float3 Rhat = R_ab_vec * Rinv;

    // --- Calculate Total Force on B due to A ---
    float3 force_b = (float3)(0.0f, 0.0f, 0.0f);

    // Order p=0
    force_b += calc_force_MM(M_a, M_b, Rhat, Rinv2);

    // Order p=1
    if (expansion_order >= 1) {
        force_b += calc_force_MD(M_a, M_b, Rhat, Rinv3);
        force_b += calc_force_DM(M_a, M_b, Rhat, Rinv3);
        force_b += calc_force_DD(M_a, M_b, Rhat, Rinv4);
    }

    // Order p=2
    if (expansion_order >= 2) {
        force_b += calc_force_MQ(M_a, M_b, Rhat, Rinv4);
        force_b += calc_force_QM(M_a, M_b, Rhat, Rinv4);
        // ** Use DQ/QD with caution - validation needed **
        force_b += calc_force_DQ(M_a, M_b, Rhat, Rinv5);
        force_b += calc_force_QD(M_a, M_b, Rhat, Rinv5);
        // QQ term omitted
        // force_b += calc_force_QQ(...);
    }

    // Apply Coulomb Constant
    force_b *= K_COULOMB;

    // --- Force Distribution (Atomic Accumulation) ---
    float3 force_a = -force_b; // Action-Reaction

    // Distribute force_a to particles in bucket a
    int count_a = bucket_particle_counts[idx_a];
    int offset_a = bucket_particle_offsets[idx_a];
    float Q_a = M_a.s3; // Get charge for distribution logic
    // Use a suitable small float threshold for charge check
    if (count_a > 0) {
        if (fabs(Q_a) > 1.0e-5f) { // Charge-proportional
            float inv_Q_a = 1.0f / Q_a;
            for (int i = 0; i < count_a; ++i) {
                int p_idx = bucket_particle_indices[offset_a + i];
                float q_i = particle_charges[p_idx];
                float3 f_contrib = force_a * (q_i * inv_Q_a);
                atomic_add_float3(&particle_forces[p_idx * 3], f_contrib);
            }
        } else { // Even distribution for neutral buckets
            float3 f_per_particle = force_a / (float)count_a;
            for (int i = 0; i < count_a; ++i) {
                int p_idx = bucket_particle_indices[offset_a + i];
                atomic_add_float3(&particle_forces[p_idx * 3], f_per_particle);
            }
        }
    }

    // Distribute force_b to particles in bucket b
    int count_b = bucket_particle_counts[idx_b];
    int offset_b = bucket_particle_offsets[idx_b];
    float Q_b = M_b.s3; // Get charge for distribution logic
    if (count_b > 0) {
         if (fabs(Q_b) > 1.0e-5f) { // Charge-proportional
            float inv_Q_b = 1.0f / Q_b;
             for (int i = 0; i < count_b; ++i) {
                int p_idx = bucket_particle_indices[offset_b + i];
                float q_i = particle_charges[p_idx];
                float3 f_contrib = force_b * (q_i * inv_Q_b);
                atomic_add_float3(&particle_forces[p_idx * 3], f_contrib);
            }
        } else { // Even distribution for neutral buckets
            float3 f_per_particle = force_b / (float)count_b;
             for (int i = 0; i < count_b; ++i) {
                int p_idx = bucket_particle_indices[offset_b + i];
                atomic_add_float3(&particle_forces[p_idx * 3], f_per_particle);
            }
        }
    }
}




#pragma OPENCL EXTENSION cl_khr_fp_atomic : enable // Or use cmpxchg fallback

// Define constants if needed
#define SQRT_PI 1.77245385f

//----------------------------------------------------------------------
// Helper: Atomic add for float
//----------------------------------------------------------------------
#ifdef cl_khr_fp_atomic
inline void atomic_addf_local(__global float *addr, float val) {
    atomic_addf(addr, val);
}
#else
inline void atomic_addf_local(__global float *addr, float val) {
    union { uint u; float f; } old_val, new_val;
    do { old_val.f = *addr; new_val.f = old_val.f + val; }
    while (atom_cmpxchg((__global uint *)addr, old_val.u, new_val.u) != old_val.u);
}
#endif

//----------------------------------------------------------------------
// Helper: Calculate Regular Solid Harmonics R_lm(v) for p=2 (UNNORMALIZED REAL BASIS)
// Basis: 1, y, z, x, xy, yz, 2z^2-x^2-y^2, xz, x^2-y^2 -> s0..s8
//----------------------------------------------------------------------
inline float16 calculate_regular_harmonics_p2(float3 v) {
    float16 Rlm = (float16)(0.0f);
    float x = v.x; float y = v.y; float z = v.z;
    float x2 = x*x; float y2 = y*y; float z2 = z*z;

    Rlm.s0 = 1.0f;                     // H_0
    Rlm.s1 = y;                        // H_1
    Rlm.s2 = z;                        // H_2
    Rlm.s3 = x;                        // H_3
    Rlm.s4 = x*y;                      // H_4
    Rlm.s5 = y*z;                      // H_5
    Rlm.s6 = 2.0f * z2 - x2 - y2;      // H_6
    Rlm.s7 = x*z;                      // H_7
    Rlm.s8 = x2 - y2;                  // H_8
    return Rlm;
}

//----------------------------------------------------------------------
// Helper: Calculate Irregular Solid Harmonics O_lm(v) for p=2 (using same basis)
// O_lm ~ R_lm(v) / |v|^(2l+1)
//----------------------------------------------------------------------
inline float16 calculate_irregular_harmonics_p2(float3 v) {
    float16 Olm = (float16)(0.0f);
    float R2 = dot(v, v);
    if (R2 < 1e-12f) { return Olm; } // Avoid division by zero
    float Rinv = rsqrt(R2); // Use reciprocal sqrt for potential speedup
    float Rinv2 = Rinv * Rinv;
    float Rinv3 = Rinv2 * Rinv;
    float Rinv5 = Rinv3 * Rinv2;

    // Get regular harmonics R_lm(v)
    float16 Rlm = calculate_regular_harmonics_p2(v);

    // O_lm = R_lm / R^(2l+1) (UNNORMALIZED)
    Olm.s0 = Rlm.s0 * Rinv;   // l=0 -> Rinv^1
    Olm.s1 = Rlm.s1 * Rinv3;  // l=1 -> Rinv^3
    Olm.s2 = Rlm.s2 * Rinv3;
    Olm.s3 = Rlm.s3 * Rinv3;
    Olm.s4 = Rlm.s4 * Rinv5;  // l=2 -> Rinv^5
    Olm.s5 = Rlm.s5 * Rinv5;
    Olm.s6 = Rlm.s6 * Rinv5;
    Olm.s7 = Rlm.s7 * Rinv5;
    Olm.s8 = Rlm.s8 * Rinv5;
    return Olm;
}


//======================================================================
// Main M2L Translation Kernel (p=2)
//======================================================================
__kernel void m2l_translate_p2(
    __global const int2* restrict m2l_pairs,
    __global const float3* restrict bucket_centers,
    __global const float16* restrict M_lm_coeffs, // Source M_lm (basis H_i)
    __global float* restrict L_lm_coeffs_global   // Target L_lm (basis H_i) accumulator
) {
    int gid = get_global_id(0);
    int idx_a = m2l_pairs[gid].x; // Source index
    int idx_b = m2l_pairs[gid].y; // Target index

    // --- Load Source Multipole Moments M_lm(a) ---
    float16 Mlm_a = M_lm_coeffs[idx_a];

    // --- Calculate Displacement Vector & Translation Operators ---
    float3 R_a = bucket_centers[idx_a];
    float3 R_b = bucket_centers[idx_b];
    float3 R_ab = R_b - R_a;
    // O_lm(-R_ab): Irregular harmonics evaluated at vector from target to source
    float16 Olm_neg_Rab = calculate_irregular_harmonics_p2(-R_ab);

    // --- Prepare source M and operator O arrays for easier indexing ---
    float M[9];
    M[0] = Mlm_a.s0; M[1] = Mlm_a.s1; M[2] = Mlm_a.s2;
    M[3] = Mlm_a.s3; M[4] = Mlm_a.s4; M[5] = Mlm_a.s5;
    M[6] = Mlm_a.s6; M[7] = Mlm_a.s7; M[8] = Mlm_a.s8;

    float O[9];
    O[0] = Olm_neg_Rab.s0; O[1] = Olm_neg_Rab.s1; O[2] = Olm_neg_Rab.s2;
    O[3] = Olm_neg_Rab.s3; O[4] = Olm_neg_Rab.s4; O[5] = Olm_neg_Rab.s5;
    O[6] = Olm_neg_Rab.s6; O[7] = Olm_neg_Rab.s7; O[8] = Olm_neg_Rab.s8;

    // --- Perform M2L Translation into L ---
    float L[9] = {0.0f}; // Holds the 9 L_lm contributions

    // L[0] (l=0) = M0*O0
    L[0] = M[0] * O[0];

    // L[1..3] (l=1) = M0*O1 + M1*O0
    L[1] = M[0] * O[1] + M[1] * O[0]; // Ly = M0*Oy + My*O0
    L[2] = M[0] * O[2] + M[2] * O[0]; // Lz = M0*Oz + Mz*O0
    L[3] = M[0] * O[3] + M[3] * O[0]; // Lx = M0*Ox + Mx*O0

    // L[4..8] (l=2) = M0*O2 + M1*O1 + M2*O0
    // L[4] (xy) = M0*Oxy + Mxy*O0 + (Mx*Oy + My*Ox)
    L[4] =  M[0] * O[4] + M[4] * O[0]
         + (M[3] * O[1] + M[1] * O[3]);

    // L[5] (yz) = M0*Oyz + Myz*O0 + (My*Oz + Mz*Oy)
    L[5] =  M[0] * O[5] + M[5] * O[0]
         + (M[1] * O[2] + M[2] * O[1]);

    // L[6] (2z2-x2-y2) = M0*O(2z2...) + M(2z2...)*O0 + 2*(Mz*Oz) - (Mx*Ox) - (My*Oy)
    L[6] =  M[0] * O[6] 
         +  M[6] * O[0]
         +  M[2] * O[2] * 2.0f 
         -  M[3] * O[3] 
         -  M[1] * O[1];

    // L[7] (xz) = M0*Oxz + Mxz*O0 + (Mx*Oz + Mz*Ox)
    L[7] = M[0] * O[7] + M[7] * O[0]
         +(M[3] * O[2] + M[2] * O[3]);

    // L[8] (x2-y2) = M0*O(x2-y2) + M(x2-y2)*O0 + (Mx*Ox - My*Oy)
    L[8] = M[0] * O[8] + M[8] * O[0]
         +(M[3] * O[3] - M[1] * O[1]);


    // --- Atomically Accumulate Results into Global L_lm Buffer ---
    __global float* L_lm_target_ptr = L_lm_coeffs_global + idx_b * 16; // Base pointer for target L_lm(b)

    // Add the computed contributions L[0..8]
    for (int i = 0; i < 9; ++i) {
         if (fabs(L[i]) > 1e-30f) { // Use a small float threshold
            atomic_addf_local(&L_lm_target_ptr[i], L[i]);
         }
    }
}