// InteractionEnergy.cl - Molecule-substrate interaction energy kernel
// Computes LJ/Morse + Coulomb + H-bond energy for rigid-body scan poses
// Uses MMFF REQ parameter layout: float4 REQ = (RvdW, sqrt(EvdW), Q, Hbond)
// Mixing rules match relax_multi.cl: R_ij = R_i + R_j, E_ij = E_i * E_j, Q_ij = Q_i * Q_j
// Each workgroup evaluates one scan pose; each thread handles one molecule atom

// Rigid-body transform: 3x3 rotation (rows packed in float4.xyz) + translation in .w columns
typedef struct { float4 row0; float4 row1; float4 row2; } Transform;

// ======== Interaction energy functions (matching MMFF conventions) ========

// LJ energy with combined parameters: E_LJ = eps * ( (r0/r)^12 - 2*(r0/r)^6 )
// r0 = R_i + R_j (sum of vdW radii), eps = E_i * E_j (product of sqrt(EvdW))
inline float energy_LJ(float r, float r0, float eps) {
    float u  = r0 / r;
    float u6 = u * u * u; u6 = u6 * u6;
    return eps * u6 * (u6 - 2.0f);
}

// Morse energy: E_Morse = De * (exp(-2*a*(r-r0)) - 2*exp(-a*(r-r0)))
inline float energy_Morse(float r, float r0, float De, float alpha) {
    float ex = exp(-alpha * (r - r0));
    return De * (ex * ex - 2.0f * ex);
}

// Coulomb energy: E_Coul = COULOMB_CONST * Q_i * Q_j / r
inline float energy_Coulomb(float r, float Qij, float Coulomb_const) {
    return Coulomb_const * Qij / r;
}

// H-bond correction energy (only for donor-acceptor pairs where H_i * H_j < 0)
inline float energy_HBond(float r, float r0, float Hij) {
    float u  = r0 / r;
    float u2 = u * u;
    float u4 = u2 * u2;
    return Hij * u4 * (u4 - 2.0f);  // 8-4 well
}

// Apply rigid transform: p_rot = R * p0 + t
inline float3 apply_transform(Transform T, float4 p0) {
    float3 p_rot;
    p_rot.x = T.row0.x * p0.x + T.row0.y * p0.y + T.row0.z * p0.z + T.row0.w;
    p_rot.y = T.row1.x * p0.x + T.row1.y * p0.y + T.row1.z * p0.z + T.row1.w;
    p_rot.z = T.row2.x * p0.x + T.row2.y * p0.y + T.row2.z * p0.z + T.row2.w;
    return p_rot;
}

inline float macro_phi_rect_dipole(float3 p, float4 P, float4 AB) {
    // Potential of a uniformly polarized rectangle centered at origin, spanning x∈[-Ax,+Ax], y∈[-By,+By] in plane z=0.
    // P = (Px,Py,Pz,0) is polarization (dipole moment per area) in [e/Ang].
    // Returns potential in units [e/Ang]. Multiply by Coulomb_const*q to get energy [eV].
    float Ax = AB.x;
    float By = AB.y;
    float x = p.x;
    float y = p.y;
    float z = p.z;
    float sumOmega = 0.0f;
    float sumLogY  = 0.0f;
    float sumLogX  = 0.0f;
    float xs[2] = {-Ax, Ax};
    float ys[2] = {-By, By};
    for (int ix=0; ix<2; ix++) {
        for (int iy=0; iy<2; iy++) {
            float X = x - xs[ix];
            float Y = y - ys[iy];
            float R = sqrt(X*X + Y*Y + z*z);
            float s = ((ix==0)?-1.0f:1.0f) * ((iy==0)?-1.0f:1.0f); // (-,-)=+; (+,-)=-; (-,+)=-; (+,+)=+
            // solid angle (for Pz)
            float denom = z * R;
            float ang = 0.0f;
            if (fabs(denom) > 1e-8f) {
                ang = atan2( X*Y, denom );
            }
            sumOmega += s * ang;
            // log terms (for Px,Py)
            sumLogY += s * log( Y + R );
            sumLogX += s * log( X + R );
        }
    }
    return (P.z * sumOmega) - (P.x * sumLogY) - (P.y * sumLogX);
}

// ======== Main kernel ========

__kernel void evaluate_interaction(
    __global const float4*     mol_atoms,      // molecule base atoms (x,y,z,0)
    __global const float4*     mol_REQs,       // molecule REQ params (R, E, Q, H)
    int                        nmol,            // number of molecule atoms
    __global const float4*     sub_atoms,       // substrate atoms (x,y,z,0)
    __global const float4*     sub_REQs,        // substrate REQ params (R, E, Q, H)
    int                        nsub,            // number of substrate atoms in unit cell
    int3                       nPBC,            // (nx, ny, nz) max PBC shifts
    float4                     lvec_a,          // lattice vector a
    float4                     lvec_b,          // lattice vector b
    float4                     lvec_c,          // lattice vector c
    __global const Transform*  transforms,      // one per scan pose
    int                        enable_LJ,       // flag: enable LJ/Morse
    int                        enable_Coulomb,   // flag: enable Coulomb
    int                        enable_HBond,     // flag: enable H-bond
    int                        enable_Morse,     // flag: Morse instead of LJ
    int                        enable_macro,     // subtract macro-potential (dipole sheet)
    float4                     macro_P,          // (Px,Py,Pz,0) polarization in [e/Ang]
    float4                     macro_AB,         // (Ax,By,0,0) rectangle half-sizes [Ang]
    float                      Coulomb_const,    // 14.3996 eV*Ang
    float                      Morse_alpha,      // Morse alpha parameter
    __local  float*            local_scores,     // [nloc] local reduction
    __global float*            results_total,    // [nconf] total energy per pose
    __global float*            results_LJ,       // [nconf] LJ component
    __global float*            results_Coul,     // [nconf] Coulomb component
    __global float*            results_HB        // [nconf] H-bond component
) {
    int iloc  = get_local_id(0);
    int nloc  = get_local_size(0);
    int iconf = get_group_id(0);
    Transform T = transforms[iconf];

    float E_lj   = 0.0f;
    float E_coul = 0.0f;
    float E_hb   = 0.0f;

    for (int im = iloc; im < nmol; im += nloc) {
        float4 p0   = mol_atoms[im];
        float4 REQm = mol_REQs[im];
        float3 p_rot = apply_transform(T, p0);

        if (enable_macro) {
            float phi = macro_phi_rect_dipole(p_rot, macro_P, macro_AB);
            E_coul -= Coulomb_const * REQm.z * phi;
        }

        float3 start_shift = lvec_a.xyz * (float)(-nPBC.x) + 
                             lvec_b.xyz * (float)(-nPBC.y) + 
                             lvec_c.xyz * (float)(-nPBC.z);

        float3 shift_c = start_shift;
        for (int ic = -nPBC.z; ic <= nPBC.z; ic++) {
            float3 shift_b = shift_c;
            for (int ib = -nPBC.y; ib <= nPBC.y; ib++) {
                float3 shift_a = shift_b;
                for (int ia = -nPBC.x; ia <= nPBC.x; ia++) {
                    
                    for (int js = 0; js < nsub; js++) {
                        float4 ps   = sub_atoms[js];
                        float4 REQs = sub_REQs[js];
            
                        // p_mol - (p_sub + shift)
                        float3 dp = p_rot - (ps.xyz + shift_a);
                        float r2 = dot(dp, dp);
                        float r  = sqrt(r2);
                        if (r < 0.1f) r = 0.1f;
            
                        float Rij = REQm.x + REQs.x;
                        float Eij = REQm.y * REQs.y;
                        float Qij = REQm.z * REQs.z;
                        float Hij = REQm.w * REQs.w;
            
                        if (enable_Morse) { E_lj += energy_Morse(r, Rij, Eij, Morse_alpha); }
                        else if (enable_LJ) { E_lj += energy_LJ(r, Rij, Eij); }
                        if (enable_Coulomb) { E_coul += energy_Coulomb(r, Qij, Coulomb_const); }
                        if (enable_HBond && (Hij < 0.0f)) { E_hb += energy_HBond(r, Rij, Hij); }
                    }
                    shift_a += lvec_a.xyz;
                }
                shift_b += lvec_b.xyz;
            }
            shift_c += lvec_c.xyz;
        }
    }

    // ---- Parallel reduction for each component ----
    local_scores[iloc] = E_lj;
    barrier(CLK_LOCAL_MEM_FENCE);
    for (int s = nloc / 2; s > 0; s >>= 1) { if (iloc < s) local_scores[iloc] += local_scores[iloc + s]; barrier(CLK_LOCAL_MEM_FENCE); }
    if (iloc == 0) results_LJ[iconf] = local_scores[0];

    barrier(CLK_LOCAL_MEM_FENCE);
    local_scores[iloc] = E_coul;
    barrier(CLK_LOCAL_MEM_FENCE);
    for (int s = nloc / 2; s > 0; s >>= 1) { if (iloc < s) local_scores[iloc] += local_scores[iloc + s]; barrier(CLK_LOCAL_MEM_FENCE); }
    if (iloc == 0) results_Coul[iconf] = local_scores[0];

    barrier(CLK_LOCAL_MEM_FENCE);
    local_scores[iloc] = E_hb;
    barrier(CLK_LOCAL_MEM_FENCE);
    for (int s = nloc / 2; s > 0; s >>= 1) { if (iloc < s) local_scores[iloc] += local_scores[iloc + s]; barrier(CLK_LOCAL_MEM_FENCE); }
    if (iloc == 0) {
        results_HB[iconf] = local_scores[0];
        results_total[iconf] = results_LJ[iconf] + results_Coul[iconf] + results_HB[iconf];
    }
}

// ======== Constrained relaxation kernel ========

__kernel void relax_constrained(
    __global const float4*     mol_atoms,
    __global const float4*     mol_REQs,
    int                        nmol,
    __global const float4*     sub_atoms,
    __global const float4*     sub_REQs,
    int                        nsub,
    int3                       nPBC,
    float4                     lvec_a,
    float4                     lvec_b,
    float4                     lvec_c,
    __global const Transform*  transforms,
    int                        enable_LJ,
    int                        enable_Coulomb,
    int                        enable_HBond,
    int                        enable_Morse,
    int                        enable_macro,
    float4                     macro_P,
    float4                     macro_AB,
    float                      Coulomb_const,
    float                      Morse_alpha,
    float                      spring_k,        // tether stiffness (like constrK)
    float                      dt,              // step size for steepest descent
    int                        nsteps,          // number of relaxation steps
    __local  float*            local_scores,
    __global float*            results_total,
    __global float*            results_LJ,
    __global float*            results_Coul,
    __global float*            results_HB,
    __global float4*           relaxed_pos      // [nconf * nmol] output relaxed positions
) {
    int iloc  = get_local_id(0);
    int nloc  = get_local_size(0);
    int iconf = get_group_id(0);
    Transform T = transforms[iconf];

    for (int im = iloc; im < nmol; im += nloc) {
        float4 p0   = mol_atoms[im];
        float4 REQm = mol_REQs[im];
        float3 anchor = apply_transform(T, p0);
        float3 pos = anchor;

        for (int step = 0; step < nsteps; step++) {
            float3 force = (float3)(0.0f, 0.0f, 0.0f);
            // Tether force: F = -k*(pos - anchor)
            force -= spring_k * (pos - anchor);

            if (enable_macro) {
                float phi = macro_phi_rect_dipole(pos, macro_P, macro_AB);
                // E = -C*q*phi  => F = -dE/dr = +C*q*grad(phi)
                // Numerical grad would be too expensive; macro correction here is energy-level only.
            }

            float3 start_shift = lvec_a.xyz * (float)(-nPBC.x) + 
                                 lvec_b.xyz * (float)(-nPBC.y) + 
                                 lvec_c.xyz * (float)(-nPBC.z);
            
            float3 shift_c = start_shift;
            for (int ic = -nPBC.z; ic <= nPBC.z; ic++) {
                float3 shift_b = shift_c;
                for (int ib = -nPBC.y; ib <= nPBC.y; ib++) {
                    float3 shift_a = shift_b;
                    for (int ia = -nPBC.x; ia <= nPBC.x; ia++) {
                        
                        for (int js = 0; js < nsub; js++) {
                            float4 ps   = sub_atoms[js];
                            float4 REQs = sub_REQs[js];
                            float3 dp = pos - (ps.xyz + shift_a);
                            float r2 = dot(dp, dp);
                            float r  = sqrt(r2);
                            if (r < 0.1f) r = 0.1f;
                            float invr = 1.0f / r;
            
                            float Rij = REQm.x + REQs.x;
                            float Eij = REQm.y * REQs.y;
                            float Qij = REQm.z * REQs.z;
                            float Hij = REQm.w * REQs.w;
                            float fr = 0.0f; // dE/dr (positive = repulsive)
            
                            if (enable_Morse) {
                                float ex = exp(-Morse_alpha * (r - Rij));
                                fr += 2.0f * Eij * Morse_alpha * ex * (ex - 1.0f);
                            } else if (enable_LJ) {
                                float u = Rij * invr;
                                float u6 = u*u*u; u6 = u6*u6;
                                fr += 12.0f * Eij * invr * u6 * (u6 - 1.0f);
                            }
                            if (enable_Coulomb) {
                                fr += -Coulomb_const * Qij * invr * invr;
                            }
                            if (enable_HBond && (Hij < 0.0f)) {
                                float u = Rij * invr;
                                float u4 = u*u; u4 = u4*u4;
                                fr += 8.0f * Hij * invr * u4 * (u4 - 1.0f);
                            }
                            force -= fr * dp * invr; // F = -dE/dr * r_hat
                        }
                        shift_a += lvec_a.xyz;
                    }
                    shift_b += lvec_b.xyz;
                }
                shift_c += lvec_c.xyz;
            }
            pos += dt * force;
        }
        relaxed_pos[iconf * nmol + im] = (float4)(pos.x, pos.y, pos.z, 0.0f);
    }

    // Compute final energy at relaxed positions
    float E_lj = 0.0f, E_coul = 0.0f, E_hb = 0.0f;
    for (int im = iloc; im < nmol; im += nloc) {
        float4 rp   = relaxed_pos[iconf * nmol + im];
        float4 REQm = mol_REQs[im];

        if (enable_macro) {
            float phi = macro_phi_rect_dipole(rp.xyz, macro_P, macro_AB);
            E_coul -= Coulomb_const * REQm.z * phi;
        }
        
        float3 start_shift = lvec_a.xyz * (float)(-nPBC.x) + 
                             lvec_b.xyz * (float)(-nPBC.y) + 
                             lvec_c.xyz * (float)(-nPBC.z);
        float3 shift_c = start_shift;
        for (int ic = -nPBC.z; ic <= nPBC.z; ic++) {
            float3 shift_b = shift_c;
            for (int ib = -nPBC.y; ib <= nPBC.y; ib++) {
                float3 shift_a = shift_b;
                for (int ia = -nPBC.x; ia <= nPBC.x; ia++) {
                    for (int js = 0; js < nsub; js++) {
                        float4 ps   = sub_atoms[js];
                        float4 REQs = sub_REQs[js];
                        float3 dp = rp.xyz - (ps.xyz + shift_a);
                        float r = sqrt(dot(dp, dp));
                        if (r < 0.1f) r = 0.1f;
                        float Rij = REQm.x + REQs.x;
                        float Eij = REQm.y * REQs.y;
                        float Qij = REQm.z * REQs.z;
                        float Hij = REQm.w * REQs.w;
                        if (enable_Morse) { E_lj += energy_Morse(r, Rij, Eij, Morse_alpha); }
                        else if (enable_LJ) { E_lj += energy_LJ(r, Rij, Eij); }
                        if (enable_Coulomb) { E_coul += energy_Coulomb(r, Qij, Coulomb_const); }
                        if (enable_HBond && (Hij < 0.0f)) { E_hb += energy_HBond(r, Rij, Hij); }
                    }
                    shift_a += lvec_a.xyz;
                }
                shift_b += lvec_b.xyz;
            }
            shift_c += lvec_c.xyz;
        }
        
        float3 anchor = apply_transform(T, mol_atoms[im]);
        float3 dd = rp.xyz - anchor;
        E_lj += 0.5f * spring_k * dot(dd, dd);
    }

    // Reduction
    local_scores[iloc] = E_lj;   barrier(CLK_LOCAL_MEM_FENCE);
    for (int s = nloc/2; s > 0; s >>= 1) { if (iloc < s) local_scores[iloc] += local_scores[iloc+s]; barrier(CLK_LOCAL_MEM_FENCE); }
    if (iloc == 0) results_LJ[iconf] = local_scores[0];
    barrier(CLK_LOCAL_MEM_FENCE);
    local_scores[iloc] = E_coul;  barrier(CLK_LOCAL_MEM_FENCE);
    for (int s = nloc/2; s > 0; s >>= 1) { if (iloc < s) local_scores[iloc] += local_scores[iloc+s]; barrier(CLK_LOCAL_MEM_FENCE); }
    if (iloc == 0) results_Coul[iconf] = local_scores[0];
    barrier(CLK_LOCAL_MEM_FENCE);
    local_scores[iloc] = E_hb;    barrier(CLK_LOCAL_MEM_FENCE);
    for (int s = nloc/2; s > 0; s >>= 1) { if (iloc < s) local_scores[iloc] += local_scores[iloc+s]; barrier(CLK_LOCAL_MEM_FENCE); }
    if (iloc == 0) {
        results_HB[iconf] = local_scores[0];
        results_total[iconf] = results_LJ[iconf] + results_Coul[iconf] + results_HB[iconf];
    }
}
