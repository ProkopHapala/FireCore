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
    int npi = (int)port_n[i];

    float2 fi = (float2)(0.0f, 0.0f);
    int i4 = i * MAX_DEGREE;

    for (int k = 0; k < MAX_DEGREE; k++) {
        int j = (k == 0) ? ng.x : (k == 1) ? ng.y : (k == 2) ? ng.z : ng.w;
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
// KERNEL: PBD constraint solve (2D version)
// Computes position corrections for port constraints
// ------------------------------------------------------------------
__kernel void compute_corrections_2d(
    const int nnode,
    __global const float2* pos,
    __global const float2* rot,
    __global const int4*   neighs,
    __global const float2* port_local,
    __global const float*  stiffness_flat,  // size nnode*4
    __global float2*       dpos_node,       // accumulated position correction
    __global float*        dtheta_node,     // accumulated rotation correction (scalar)
    __global float2*       dpos_neigh       // corrections for neighbors
) {
    int i = get_global_id(0);
    if (i >= nnode) return;

    float2 p_i = pos[i];
    float2 z_i = rot[i];

    int4 ng = neighs[i];
    float2 dpos = (float2)(0.0f, 0.0f);
    float dtheta = 0.0f;

    int i4 = i * MAX_DEGREE;

    for (int k = 0; k < MAX_DEGREE; k++) {
        int j = (k == 0) ? ng.x : (k == 1) ? ng.y : (k == 2) ? ng.z : ng.w;
        if (j < 0) continue;

        float K = stiffness_flat[i4 + k];
        if (K <= 0.0f) continue;

        float2 r_local = port_local[i4 + k];
        float2 r_world = cmplx_rotate(z_i, r_local);
        float2 tip = p_i + r_world;
        float2 p_j = pos[j];

        float2 diff = tip - p_j;  // constraint violation
        float C2 = dot(diff, diff);
        if (C2 < 1e-12f) continue;

        // Simplified mass weights (uniform mass)
        float W_i = 1.0f;
        float W_j = 1.0f;

        // Rotation correction: dtheta contribution from r x diff
        // In 2D: dtheta ~ (r.x * diff.y - r.y * diff.x) / |r|^2
        float r2 = dot(r_world, r_world);
        float r_cross_diff = r_world.x * diff.y - r_world.y * diff.x;
        float dtheta_k = (r2 > 1e-12f) ? (r_cross_diff / r2) : 0.0f;

        dpos += diff * W_i;
        dtheta += dtheta_k * W_i;
        dpos_neigh[i4 + k] = -diff * W_j;
    }

    dpos_node[i] = dpos;
    dtheta_node[i] = dtheta;
}

// ------------------------------------------------------------------
// KERNEL: Apply corrections with relaxation (2D)
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
    const float relaxation
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

    pos[i] += corr * relaxation;

    // Apply rotation correction (nodes only)
    if (i < nnode) {
        float dtheta = dtheta_node[i] * relaxation;
        float2 dz = cmplx_from_angle(dtheta);
        rot[i] = cmplx_mul(dz, rot[i]);
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
    float2 dpos = (float2)(0.0f, 0.0f);
    float dtheta = 0.0f;

    int i4 = i * MAX_DEGREE;

    for (int k = 0; k < MAX_DEGREE; k++) {
        int j = (k == 0) ? ng.x : (k == 1) ? ng.y : (k == 2) ? ng.z : ng.w;
        if (j < 0) continue;

        float K = (k == 0) ? stiffness[i].x : (k == 1) ? stiffness[i].y : (k == 2) ? stiffness[i].z : stiffness[i].w;
        if (K <= 0.0f) continue;

        float alpha = 1.0f / (K * dt * dt + 1e-12f);
        int idx = i4 + k;
        float lambda_prev = lambda[idx];

        float2 r_local = port_local[idx];
        float2 r_world = cmplx_rotate(z_i, r_local);
        float2 tip = p_i + r_world;
        float2 p_j = pos[j];

        float2 diff = tip - p_j;
        float C = sqrt(dot(diff, diff));
        if (C < 1e-8f) continue;

        float2 n = diff / C;

        // Simplified compliance-weighted delta
        float W_total = 2.0f + alpha;  // simplified: node + neighbor + compliance
        float delta_lambda = (-C - alpha * lambda_prev) / W_total;

        lambda[idx] += delta_lambda;

        float2 dp_i = n * delta_lambda;
        float2 dp_j = -n * delta_lambda;

        dpos += dp_i;

        // Rotation correction
        float r2 = dot(r_world, r_world);
        float r_cross_n = r_world.x * n.y - r_world.y * n.x;
        float dtheta_k = (r2 > 1e-12f) ? (r_cross_n * delta_lambda / r2) : 0.0f;
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
// KERNEL: Gather and apply XPBD corrections (2D)
// ------------------------------------------------------------------
__kernel void gather_and_apply_xpbd_2d(
    const int natoms,
    const int nnode,
    __global float2*       pos,
    __global float2*       rot,
    __global const int4*   bkSlots,
    __global const float2* dpos_neigh,
    __global const float2* dpos_node,
    __global const float*  dtheta_node
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

    pos[i] += corr;

    if (i < nnode) {
        float dtheta = dtheta_node[i];
        float2 dz = cmplx_from_angle(dtheta);
        rot[i] = cmplx_mul(dz, rot[i]);
    }
}
