// XPBD_2D.cl - 2D Position Based Dynamics with Complex Number Rotation
// In 2D, rotation is represented by complex numbers z = (cosθ, sinθ)
// Rotating a vector v by z: z * v = (z.x*v.x - z.y*v.y, z.y*v.x + z.x*v.y)

// ------------------------------------------------------------------
// CONFIGURATION
// ------------------------------------------------------------------
#define MAX_DEGREE     4      // Max ports per node

// ------------------------------------------------------------------
// 2D COMPLEX NUMBER HELPERS (rotation representation)
// ------------------------------------------------------------------

// Rotate vector v by complex rotation z = (cosθ, sinθ)
inline float2 cmplx_rotate(float2 z, float2 v) {
    return (float2)(
        z.x * v.x - z.y * v.y,  // real part
        z.y * v.x + z.x * v.y   // imag part
    );
}

// Complex multiplication: z1 * z2
inline float2 cmplx_mul(float2 z1, float2 z2) {
    return (float2)(
        z1.x * z2.x - z1.y * z2.y,
        z1.y * z2.x + z1.x * z2.y
    );
}

// Get complex rotation from angle (radians)
inline float2 cmplx_from_angle(float angle) {
    return (float2)(cos(angle), sin(angle));
}

// Get angle from complex rotation
inline float cmplx_angle(float2 z) {
    return atan2(z.y, z.x);
}

inline float2 cmplx_normalize(float2 z) {
    float r2 = z.x*z.x + z.y*z.y;
    if (r2 < 1e-20f) return (float2)(1.0f, 0.0f);
    float invr = rsqrt(r2);
    return z * invr;
}

// ------------------------------------------------------------------
// KERNEL: Initialize heavy-ball buffers
// ------------------------------------------------------------------
__kernel void init_hb_pos_2d(
    const int natoms,
    __global const float2* pos,
    __global float2* hb_pos
) {
    int i = get_global_id(0);
    if (i >= natoms) return;
    hb_pos[i] = pos[i];
}

__kernel void init_hb_rot_2d(
    const int nnode,
    __global const float2* rot,
    __global float2* hb_rot
) {
    int i = get_global_id(0);
    if (i >= nnode) return;
    hb_rot[i] = cmplx_normalize(rot[i]);
}

// ------------------------------------------------------------------
// KERNEL: Clear forces and node buffers
// ------------------------------------------------------------------
__kernel void clear_2d_forces(
    const int natoms,
    __global float2* force
) {
    int i = get_global_id(0);
    if (i >= natoms) return;
    force[i] = (float2)(0.0f, 0.0f);
}

__kernel void clear_2d_node_buffers(
    const int nnode,
    __global float2* fneigh,   // size nnode*4
    __global float2* pneigh    // size nnode*4 (lever arms)
) {
    int i = get_global_id(0);
    int n = nnode * MAX_DEGREE;
    if (i >= n) return;
    fneigh[i] = (float2)(0.0f, 0.0f);
    pneigh[i] = (float2)(0.0f, 0.0f);
}

// ------------------------------------------------------------------
// KERNEL: Gather port forces (2D version)
// Each node gathers forces from its ports connected to neighbors
// ------------------------------------------------------------------
__kernel void gather_port_forces_2d(
    const int nnode,
    __global const float2* pos,           // w component stores invMass in .x (reuse)
    __global const float2* rot,           // complex rotation per node
    __global const int4*   neighs,        // neighbor indices (up to 4)
    __global const float4* stiffness,     // K per port (only .x used)
    __global const float2* port_local,    // local port offsets [nnode*4]
    __global const uchar*  port_n,        // number of active ports per node
    __global float2*       force,         // output force per atom
    __global float2*       fneigh,        // recoil force per port slot
    __global float2*       pneigh         // lever arm per port slot
) {
    int i = get_global_id(0);
    if (i >= nnode) return;

    float2 p_i = pos[i];
    float2 z_i = rot[i];

    int4 ng = neighs[i];
    int* neighbors = (int*)&ng;
    int npi = (int)port_n[i];

    float2 fi = (float2)(0.0f, 0.0f);
    int i4 = i * MAX_DEGREE;

    for (int k = 0; k < MAX_DEGREE; k++) {
        int j = neighbors[k];
        if (j < 0) continue;

        float Kij = (k == 0) ? stiffness[i].x : (k == 1) ? stiffness[i].y : (k == 2) ? stiffness[i].z : stiffness[i].w;
        if (Kij <= 0.0f) continue;

        float2 p_j = pos[j];

        // Rotate local port offset to global frame
        float2 r_arm = (float2)(0.0f, 0.0f);
        if (k < npi) {
            r_arm = cmplx_rotate(z_i, port_local[i4 + k]);
        }

        float2 tip = p_i + r_arm;
        float2 f = (p_j - tip) * Kij;

        pneigh[i4 + k] = r_arm;
        fi += f;
        fneigh[i4 + k] = -f;
    }

    force[i] = fi;
}

// ------------------------------------------------------------------
// KERNEL: Explicit integration with force gather (2D)
// ------------------------------------------------------------------
__kernel void integrate_2d_explicit(
    const int natoms,
    const int nnode,
    __global float2*       pos,        // .x = invMass (stored separately for capping atoms)
    __global float2*       vel,
    __global float2*       rot,        // complex rotation
    __global float*        omega,      // angular velocity (scalar in 2D)
    __global const int4*   bkSlots,    // back-slot indices for recoil gather
    __global float2*       force,
    __global float2*       fneigh,     // recoil forces from nodes
    __global float2*       pneigh,     // lever arms
    const float dt,
    const float damp
) {
    int i = get_global_id(0);
    if (i >= natoms) return;

    // For simplicity: uniform mass. invMass passed via separate buffer or constant
    float invMass = 1.0f;  // Could be passed as additional buffer
    float mass = 1.0f;
    float invI = 1.0f / (0.4f * mass);  // Simplified inertia for 2D

    float2 f = force[i];

    // Gather recoil forces from connected nodes
    int4 bk = bkSlots[i];
    if (bk.x >= 0) f += fneigh[bk.x];
    if (bk.y >= 0) f += fneigh[bk.y];
    if (bk.z >= 0) f += fneigh[bk.z];
    if (bk.w >= 0) f += fneigh[bk.w];

    // Linear integration
    float2 v = vel[i];
    v = v * damp + f * invMass * dt;
    pos[i] = pos[i] + v * dt;
    vel[i] = v;

    // Rotational integration (only for nodes)
    if (i < nnode) {
        float2 z = rot[i];

        // Compute torque from recoil forces
        float tau = 0.0f;
        int i4 = i * MAX_DEGREE;
        for (int k = 0; k < MAX_DEGREE; k++) {
            float2 r_arm = pneigh[i4 + k];
            float2 Fk = -fneigh[i4 + k];
            // 2D cross product: tau = r.x * F.y - r.y * F.x
            tau += r_arm.x * Fk.y - r_arm.y * Fk.x;
        }

        // Integrate angular velocity
        float w = omega[i];
        w = w * damp + tau * invI * dt;
        omega[i] = w;

        // Update rotation: integrate angle then convert to complex
        float angle = w * dt;
        float2 dz = cmplx_from_angle(angle);
        rot[i] = cmplx_mul(dz, z);
    }
}

// ------------------------------------------------------------------
// KERNEL: Compute real velocities from position changes (for MD mode)
// v = (x_new - x_old) / dt
// ------------------------------------------------------------------
__kernel void compute_velocities_from_positions(
    const int natoms,
    __global const float2* pos_new,
    __global const float2* pos_old,
    __global float2* vel,
    const float dt,
    const float damp
) {
    int i = get_global_id(0);
    if (i >= natoms) return;
    
    float2 dx = pos_new[i] - pos_old[i];
    float inv_dt = 1.0f / (dt + 1e-20f);
    vel[i] = dx * inv_dt * damp;
}

// ------------------------------------------------------------------
// KERNEL: Compute angular velocities from rotation changes (for MD mode)
// omega = (angle_new - angle_old) / dt
// ------------------------------------------------------------------
__kernel void compute_angular_velocities_from_rotations(
    const int nnode,
    __global const float2* rot_new,
    __global const float2* rot_old,
    __global float* omega,
    const float dt,
    const float damp
) {
    int i = get_global_id(0);
    if (i >= nnode) return;
    
    float2 z_new = rot_new[i];
    float2 z_old = rot_old[i];
    
    // Compute relative rotation: z_rel = conj(z_old) * z_new
    // In 2D: conj(z) = (z.x, -z.y)
    float2 z_old_conj = (float2)(z_old.x, -z_old.y);
    float2 z_rel = cmplx_mul(z_old_conj, z_new);
    
    // Extract angle from z_rel = (cos(dtheta), sin(dtheta))
    float dtheta = atan2(z_rel.y, z_rel.x);
    
    float inv_dt = 1.0f / (dt + 1e-20f);
    omega[i] = dtheta * inv_dt * damp;
}

// ------------------------------------------------------------------
// KERNEL: PBD constraint solve (2D version)
// Computes position corrections for port constraints with mass (M/dt^2) diagonal term
// xi_cor = (xi_pred * a + sum_j Kij * xj) / (a + sum_j Kij) where a = M/dt^2
// ------------------------------------------------------------------
__kernel void compute_corrections_2d(
    const int nnode,
    __global const float2* pos,
    __global const float2* rot,
    __global const int4*   neighs,
    __global const float2* port_local,
    __global const float*  stiffness_flat,  // size nnode*4
    __global const float*  mass,            // atom masses [natoms]
    __global float2*       dpos_node,       // accumulated position correction
    __global float*        dtheta_node,     // accumulated rotation correction (scalar)
    __global float2*       dpos_neigh,      // corrections for neighbors
    const float dt
) {
    int i = get_global_id(0);
    if (i >= nnode) return;

    float2 p_i = pos[i];
    float2 z_i = rot[i];
    
    // Projective Dynamics diagonal term: a = M/dt^2
    float dt2 = dt * dt;
    float a_i = mass[i] / (dt2 + 1e-20f);

    int4 ng = neighs[i];
    int* neighbors = (int*)&ng;
    float2 dpos = (float2)(0.0f, 0.0f);
    float dtheta = 0.0f;

    int i4 = i * MAX_DEGREE;

    for (int k = 0; k < MAX_DEGREE; k++) {
        int j = neighbors[k];
        if (j < 0) continue;

        float K = stiffness_flat[i4 + k];
        if (K <= 0.0f) continue;

        float2 r_local = port_local[i4 + k];
        float2 r_world = cmplx_rotate(z_i, r_local);
        float2 tip = p_i + r_world;

        float2 p_j = pos[j];

        float2 diff = tip - p_j;  // constraint violation (port->atom)
        float C2 = dot(diff, diff);
        if (C2 < 1e-12f) continue;

        float C = sqrt(C2);
        float2 n = diff / C;

        // Projective Dynamics with inertia term a_i = M/dt^2
        // Denominator: a_i + sum_j Kij
        float I_i = 0.4f;
        float r_cross_n = r_world.x * n.y - r_world.y * n.x;
        float w_i_trans = 1.0f;
        float w_i_rot   = (I_i > 1e-12f) ? (r_cross_n * r_cross_n) / I_i : 0.0f;
        float w_j       = 1.0f;
        float w_total   = a_i + K * (w_i_trans + w_i_rot + w_j);
        float delta_lambda = (-C * K) / (w_total + 1e-12f);

        float2 dp_i = delta_lambda * n * w_i_trans;
        float2 dp_j = -delta_lambda * n * w_j;

        float dtheta_k = (I_i > 1e-12f) ? (delta_lambda * r_cross_n / I_i) : 0.0f;

        dpos += dp_i;
        dtheta += dtheta_k;
        dpos_neigh[i4 + k] = dp_j;
    }

    dpos_node[i] = dpos;
    dtheta_node[i] = dtheta;
}

// ------------------------------------------------------------------
// KERNEL: Initialize momentum buffers to zero
// ------------------------------------------------------------------
__kernel void init_mom_pos_2d(
    const int natoms,
    __global float2* mom_pos
) {
    int i = get_global_id(0);
    if (i >= natoms) return;
    mom_pos[i] = (float2)(0.0f, 0.0f);
}

__kernel void init_mom_rot_2d(
    const int nnode,
    __global float2* mom_rot
) {
    int i = get_global_id(0);
    if (i >= nnode) return;
    mom_rot[i] = (float2)(0.0f, 0.0f);
}

// ------------------------------------------------------------------
// KERNEL: Apply corrections with relaxation and proper momentum mixing (2D)
// ------------------------------------------------------------------
__kernel void apply_corrections_2d(
    const int natoms,
    const int nnode,
    __global float2*       pos,
    __global float2*       rot,
    __global const int4*   bkSlots,
    __global const float2* dpos_node,
    __global const float*  dtheta_node,
    __global const float2* dpos_neigh,
    const float relaxation,
    __global float2*       hb_pos,
    __global float2*       hb_rot,
    __global float2*       mom_pos,
    __global float2*       mom_rot,
    const float bmix_pos,
    const float bmix_rot
) {
    int i = get_global_id(0);
    if (i >= natoms) return;

    // Apply position correction
    float2 corr = (float2)(0.0f, 0.0f);

    if (i < nnode) {
        corr += dpos_node[i];
    }

    // Gather corrections from neighbors
    int4 bk = bkSlots[i];
    if (bk.x >= 0) corr += dpos_neigh[bk.x];
    if (bk.y >= 0) corr += dpos_neigh[bk.y];
    if (bk.z >= 0) corr += dpos_neigh[bk.z];
    if (bk.w >= 0) corr += dpos_neigh[bk.w];

    // Relaxation step
    float2 p_corr = pos[i] + corr * relaxation;

    // Heavy-ball momentum mixing for position (C++ style)
    if (bmix_pos > 1e-6f) {
        float2 d_prev = mom_pos[i];        // momentum direction from previous iter
        float2 p_prev = hb_pos[i];         // previous iterate
        float2 p_next = p_corr + d_prev * bmix_pos;
        mom_pos[i] = p_next - p_prev;      // store new momentum direction
        hb_pos[i] = p_next;
        pos[i] = p_next;
    } else {
        pos[i] = p_corr;
        if (bmix_pos > -0.5f) {            // Only update hb if not explicitly disabled
            hb_pos[i] = p_corr;
        }
    }

    // Apply rotation correction (nodes only)
    if (i < nnode) {
        float dtheta = dtheta_node[i] * relaxation;
        float2 dz = cmplx_from_angle(dtheta);
        float2 r_corr = cmplx_mul(dz, rot[i]);
        
        // Heavy-ball momentum mixing for rotation
        if (bmix_rot > 1e-6f) {
            float2 d_prev = mom_rot[i];
            float2 r_prev = hb_rot[i];
            // For rotation, we apply momentum in the tangent space
            // r_next = normalize(r_corr + bmix * momentum direction)
            float2 r_next = r_corr + d_prev * bmix_rot;
            mom_rot[i] = r_next - r_prev;
            hb_rot[i] = r_next;
            rot[i] = cmplx_normalize(r_next);
        } else {
            rot[i] = cmplx_normalize(r_corr);
            if (bmix_rot > -0.5f) {
                hb_rot[i] = cmplx_normalize(r_corr);
            }
        }
    }
}

// ------------------------------------------------------------------
// KERNEL: XPBD constraint solve with lambda accumulation (2D)
// ------------------------------------------------------------------
__kernel void compute_xpbd_corrections_2d(
    const int nnode,
    __global const float2* pos,
    __global const float2* rot,
    __global const int4*   neighs,
    __global const float4* stiffness,     // per-node stiffness
    __global const float2* port_local,
    __global float*        lambda,        // accumulated multipliers [nnode*4]
    __global float2*       dpos_neigh,
    __global float2*       dpos_node,
    __global float*        dtheta_node,
    const float dt
) {
    int i = get_global_id(0);
    if (i >= nnode) return;

    float2 p_i = pos[i];
    float2 z_i = rot[i];

    int4 ng = neighs[i];
    int* neighbors = (int*)&ng;
    float2 dpos = (float2)(0.0f, 0.0f);
    float dtheta = 0.0f;

    int i4 = i * MAX_DEGREE;

    for (int k = 0; k < MAX_DEGREE; k++) {
        int j = neighbors[k];
        if (j < 0) continue;

        float K = (k == 0) ? stiffness[i].x : (k == 1) ? stiffness[i].y : (k == 2) ? stiffness[i].z : stiffness[i].w;
        if (K <= 0.0f) continue;

        int idx = i4 + k;
        float lambda_prev = lambda[idx];

        float2 r_local = port_local[idx];
        float2 r_world = cmplx_rotate(z_i, r_local);
        float2 tip = p_i + r_world;

        float2 p_j = pos[j];

        float2 diff = tip - p_j;
        float C2 = dot(diff, diff);
        if (C2 < 1e-16f) continue;

        float C = sqrt(C2);
        float2 n = diff / C;

        // Generalized mass (translation + rotation coupling) + XPBD compliance, like XPDB_new.cl
        float I_i = 0.4f;
        float dt2 = dt * dt;

        float r_cross_n = r_world.x * n.y - r_world.y * n.x;

        float w_i_trans = 1.0f;
        float w_i_rot   = (I_i > 1e-12f) ? (r_cross_n * r_cross_n) / I_i : 0.0f;
        float w_j       = 1.0f;

        float alpha = 1.0f / (K + 1e-12f);
        float alpha_tilde = alpha / (dt2 + 1e-20f);

        float w_total = w_i_trans + w_i_rot + w_j + alpha_tilde;
        float delta_lambda = (-C - alpha_tilde * lambda_prev) / (w_total + 1e-12f);
        lambda[idx] = lambda_prev + delta_lambda;

        float2 dp_i = delta_lambda * n * w_i_trans;
        float2 dp_j = -delta_lambda * n * w_j;
        float dtheta_k = (I_i > 1e-12f) ? (delta_lambda * r_cross_n / I_i) : 0.0f;

        dpos += dp_i;
        dtheta += dtheta_k;
        dpos_neigh[idx] = dp_j;
    }

    dpos_node[i] = dpos;
    dtheta_node[i] = dtheta;
}

// ------------------------------------------------------------------
// KERNEL: Reset lambda for XPBD
// ------------------------------------------------------------------
__kernel void reset_lambda_2d(
    const int n,
    __global float* lambda
) {
    int i = get_global_id(0);
    if (i >= n) return;
    lambda[i] = 0.0f;
}

// ------------------------------------------------------------------
// KERNEL: XPBD Corrections with DEBUG PRINTS (diagnostic version)
// ------------------------------------------------------------------
__kernel void compute_xpbd_corrections_2d_debug(
    const int nnode,
    __global const float2* pos,
    __global const float2* rot,
    __global const int4*   neighs,
    __global const float4* stiffness,     // per-node stiffness
    __global const float2* port_local,
    __global float*        lambda,        // accumulated multipliers [nnode*4]
    __global float2*       dpos_neigh,
    __global float2*       dpos_node,
    __global float*        dtheta_node,
    const float dt,
    const int debug_step,                 // current step for gating output
    const int max_debug_steps,            // max steps to output (e.g., 5)
    // Diagnostic output buffers (only written for debug_step < max_debug_steps)
    __global float*        dbg_data,      // [nnode*4*5] = C, lambda, dtheta, K, alpha per constraint
    __global float2*       dbg_dpos_i,    // [nnode*4] position correction on node i
    __global float2*       dbg_dpos_j,    // [nnode*4] position correction on neighbor j
    __global float2*       dbg_r_world,   // [nnode*4] rotated port vector
    __global float2*       dbg_n          // [nnode*4] constraint normal
) {
    int i = get_global_id(0);
    if (i >= nnode) return;

    float2 p_i = pos[i];
    float2 z_i = rot[i];

    int4 ng = neighs[i];
    int* neighbors = (int*)&ng;
    float2 dpos = (float2)(0.0f, 0.0f);
    float dtheta = 0.0f;

    int i4 = i * MAX_DEGREE;
    int base_offset = i4;  // For diagnostic arrays

    for (int k = 0; k < MAX_DEGREE; k++) {
        int j = neighbors[k];
        int idx = i4 + k;

        // Initialize diagnostics to zero (for inactive constraints)
        if (debug_step < max_debug_steps) {
            dbg_data[base_offset + k] = 0.0f;           // C
            dbg_data[base_offset + k + nnode*4] = 0.0f; // lambda
            dbg_data[base_offset + k + nnode*8] = 0.0f; // dtheta
            dbg_data[base_offset + k + nnode*12] = 0.0f; // K
            dbg_data[base_offset + k + nnode*16] = 0.0f; // alpha
            dbg_dpos_i[idx] = (float2)(0.0f, 0.0f);
            dbg_dpos_j[idx] = (float2)(0.0f, 0.0f);
            dbg_r_world[idx] = (float2)(0.0f, 0.0f);
            dbg_n[idx] = (float2)(0.0f, 0.0f);
        }

        if (j < 0) continue;

        float K = (k == 0) ? stiffness[i].x : (k == 1) ? stiffness[i].y : (k == 2) ? stiffness[i].z : stiffness[i].w;
        if (K <= 0.0f) continue;

        float lambda_prev = lambda[idx];

        float2 r_local = port_local[idx];
        float2 r_world = cmplx_rotate(z_i, r_local);
        float2 tip = p_i + r_world;

        float2 p_j = pos[j];

        float2 diff = tip - p_j;
        float C2 = dot(diff, diff);
        if (C2 < 1e-16f) continue;

        float C = sqrt(C2);
        float2 n = diff / C;

        // Generalized mass (translation + rotation coupling) + XPBD compliance
        float I_i = 0.4f;
        float dt2 = dt * dt;

        float r_cross_n = r_world.x * n.y - r_world.y * n.x;

        float w_i_trans = 1.0f;
        float w_i_rot   = (I_i > 1e-12f) ? (r_cross_n * r_cross_n) / I_i : 0.0f;
        float w_j       = 1.0f;

        float alpha = 1.0f / (K + 1e-12f);
        float alpha_tilde = alpha / (dt2 + 1e-20f);

        float w_total = w_i_trans + w_i_rot + w_j + alpha_tilde;
        float delta_lambda = (-C - alpha_tilde * lambda_prev) / (w_total + 1e-12f);
        lambda[idx] = lambda_prev + delta_lambda;

        float2 dp_i = delta_lambda * n * w_i_trans;
        float2 dp_j = -delta_lambda * n * w_j;
        float dtheta_k = (I_i > 1e-12f) ? (delta_lambda * r_cross_n / I_i) : 0.0f;

        // Store diagnostics for active constraints
        if (debug_step < max_debug_steps) {
            dbg_data[base_offset + k] = C;
            dbg_data[base_offset + k + nnode*4] = lambda[idx];
            dbg_data[base_offset + k + nnode*8] = dtheta_k;
            dbg_data[base_offset + k + nnode*12] = K;
            dbg_data[base_offset + k + nnode*16] = alpha;
            dbg_dpos_i[idx] = dp_i;
            dbg_dpos_j[idx] = dp_j;
            dbg_r_world[idx] = r_world;
            dbg_n[idx] = n;
        }

        dpos += dp_i;
        dtheta += dtheta_k;
        dpos_neigh[idx] = dp_j;
    }

    dpos_node[i] = dpos;
    dtheta_node[i] = dtheta;
}

// ------------------------------------------------------------------
// KERNEL: Gather and apply XPBD corrections (2D) with momentum buffers
// ------------------------------------------------------------------
__kernel void gather_and_apply_xpbd_2d(
    const int natoms,
    const int nnode,
    __global float2*       pos,
    __global float2*       rot,
    __global const int4*   bkSlots,
    __global const float2* dpos_neigh,
    __global const float2* dpos_node,
    __global const float*  dtheta_node,
    __global float2*       hb_pos,
    __global float2*       hb_rot,
    __global float2*       mom_pos,
    __global float2*       mom_rot,
    const float bmix_pos,
    const float bmix_rot
) {
    int i = get_global_id(0);
    if (i >= natoms) return;

    float2 corr = (float2)(0.0f, 0.0f);

    if (i < nnode) {
        corr += dpos_node[i];
    }

    int4 bk = bkSlots[i];
    if (bk.x >= 0) corr += dpos_neigh[bk.x];
    if (bk.y >= 0) corr += dpos_neigh[bk.y];
    if (bk.z >= 0) corr += dpos_neigh[bk.z];
    if (bk.w >= 0) corr += dpos_neigh[bk.w];

    // For XPBD, no relaxation factor - corrections are already scaled
    float2 p_corr = pos[i] + corr;

    // Heavy-ball momentum mixing for position
    if (bmix_pos > 1e-6f) {
        float2 d_prev = mom_pos[i];
        float2 p_prev = hb_pos[i];
        float2 p_next = p_corr + d_prev * bmix_pos;
        mom_pos[i] = p_next - p_prev;
        hb_pos[i] = p_next;
        pos[i] = p_next;
    } else {
        pos[i] = p_corr;
        if (bmix_pos > -0.5f) {
            hb_pos[i] = p_corr;
        }
    }

    if (i < nnode) {
        float dtheta = dtheta_node[i];
        float2 dz = cmplx_from_angle(dtheta);
        float2 r_corr = cmplx_mul(dz, rot[i]);
        
        // Heavy-ball momentum mixing for rotation
        if (bmix_rot > 1e-6f) {
            float2 d_prev = mom_rot[i];
            float2 r_prev = hb_rot[i];
            float2 r_next = r_corr + d_prev * bmix_rot;
            mom_rot[i] = r_next - r_prev;
            hb_rot[i] = r_next;
            rot[i] = cmplx_normalize(r_next);
        } else {
            rot[i] = cmplx_normalize(r_corr);
            if (bmix_rot > -0.5f) {
                hb_rot[i] = cmplx_normalize(r_corr);
            }
        }
    }
}
