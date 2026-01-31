## USER

I have classical forcefield and I'm having special way how to represent angles as axillary bonds. But I'm thinking how to best represent planarization of SP3 Carbon (aromatic) without torsions and dihedrals. I wanted to apply some term like shift atom I so that it is in plane of the other atoms. But it should work as position based dynamics and I worry that it will not preserve linear and angular momentum.
I mean in position based dynamics every Jacobi update only move the central atom ignoring the recoil forces on neighbors. In bonds this is fine as update of pi, a d update of pj although independent are symetric compensating the recoil but in this planarization moves (where pi is moved into plane of it's neighbors) I'm not sure about this symmetry and recoil compensation.
I realized that doing it properly (tetrahednon volume constrians, point-to plane contrains) is rather costly and require project recoils back. Instead I was trying to add axuliary particle (dummy atom) representing pi-orbital. That produce more lean forcefield. That axuliary particle can have distance constrains with the other atoms (neighbors of p0,   call them p1,p2,p3) and in that case my jacobi update can update onnly pi and the momentum will be conserved. But I'm thinking if I want to do it different way where I use angle constrain between pi-dummy and neighbors of p0, then I do not see why to write that jacobi update properly without writing forces to the other atoms (?).
would something change if the pi-dummy instead of being partilce (position) is just normalized orientation vectro like in this article
Stable Cosserat Rods
https://graphics.cs.utah.edu/research/projects/stable-cosserat-rods/ 
https://github.com/jerry060599/StableCosseratRods 

OK, finally I modified the force-based solver like this

notice few details and desing decisions.
in `gather_port_forces` we
1) we project ports from local frame to world frame using qaternion rotation of atom i q[i]
2) then we compute the distance force fij = Kij*(pj - port[k]), where j=neigh[k], and store them to slots located under atom i, fneigh[i*4+k], k is index of neighbor
3) we also write the force on the central atom i = sum_j fij

important detail is that this kernel run only over node atoms, and many arrays like fneigh qrot etc are allocated only for nneigh not for all natom (i.e. capping atoms bonded to single ohost, like hydrogens are considered as point particles without rotation)

then in `integrate_rigid_explicit` we 
1) gather the recoil forces forces from the atoms using  bkNeigh array which store not index of neihgbor atom but index of slot on neibhor atom from which we should accumulate the recoild force dedicated to this atom (if it is understanable) we build this map on host.
2) we evaluate torque due to recoil forces on ports (from neighbors) and undate the rotation, we also accumulate this recoil force and update the translation, both using leapfrog. 

please analyse the code think it through and tell me if you undrestand, try to explain you in maximaly undertanbale didacti way, and try to coment on efficiency, physical correctnes, conservation of angular linear momentum, efficient use of the global/buffers


```python
def make_bk_slots(neighs, *, nnode, natoms=None):
    if natoms is None:
        natoms = int(neighs.shape[0])
    bkSlots = np.full((natoms, 4), -1, dtype=np.int32)
    bkCount = np.zeros((natoms,), dtype=np.int32)
    for ia in range(int(nnode)):
        for k in range(4):
            ja = int(neighs[ia, k])
            if ja < 0:
                continue
            s = int(bkCount[ja])
            if s >= 4:
                raise RuntimeError(f"bkSlots overflow: atom {ja} has >4 back slots (from node {ia})")
            bkSlots[ja, s] = ia * 4 + k
            bkCount[ja] += 1
    return bkSlots
```


```OpenCL
__kernel void gather_port_forces(
    const int nnode,
    __global const float4* pos,
    __global const float4* quat,
    __global int4*   neighs,
    __global float4* bKs,
    __global const float4* port_local,
    __global const uchar*  port_n,
    __global float4* force,
    __global float4* fneigh,
    __global float4* pneigh
) {
    int i = get_global_id(0);
    if (i >= nnode) return;

    float3 p_i = pos[i].xyz;
    float4 q_i = quat[i];

    int4 ng = neighs[i];
    int* neighbors = (int*)&ng;

    int npi = (int)port_n[i];

    float3 fi = (float3)(0.0f);

    int i4 = i * 4;
    for (int k = 0; k < 4; k++) {
        int j = neighbors[k];
        if (j < 0) break;

        float Kij = fetch4(bKs, i, k);
        float3 p_j = pos[j].xyz;

        float3 r_arm = (float3)(0.0f);
        if (k < npi) r_arm = quat_rotate(q_i, port_local[i4 + k].xyz);

        float3 tip = p_i + r_arm;
        float3 f   = (p_j - tip) * Kij;
        pneigh[i4 + k] = (float4)(r_arm, 0.0f);

        fi += f;
        fneigh[i4 + k] = (float4)(-f, 0.0f);
    }

    force[i]  = (float4)(fi, 0.0f);
}


__kernel void integrate_rigid_explicit(
    const int natoms,
    const int nnode,
    __global float4* pos,
    __global float4* vel,
    __global float4* quat,
    __global float4* omega,
    __global int4*   bkSlots,
    __global float4* force,
    __global float4* fneigh,
    __global float4* pneigh,
    const float dt,
    const float damp
) {
    int i = get_global_id(0);
    if (i >= natoms) return;

    float invMass = pos[i].w;
    float mass    = (invMass > 1e-12f) ? (1.0f / invMass) : 1e12f;
    float invI    = 1.0f / (0.4f * mass + 1e-12f);

    float3 f = force[i].xyz;

    int4 bk = bkSlots[i];
    int* islots = (int*)&bk;
    for (int k = 0; k < 4; k++) {
        int islot = islots[k];
        if (islot >= 0) {
            f += fneigh[islot].xyz;
        }
    }

    float3 v = vel[i].xyz;
    v = v * damp + (f * invMass) * dt;
    float3 p = pos[i].xyz + v * dt;
    vel[i] = (float4)(v, 0.0f);
    pos[i] = (float4)(p, invMass);

    if (i < nnode) {
        float4 q = quat[i];
        float3 tau = (float3)(0.0f);
        int i4 = i * 4;
        for (int k = 0; k < 4; k++) {
            float3 r_arm = pneigh[i4 + k].xyz;
            float3 Fk    = -fneigh[i4 + k].xyz;
            tau += cross(r_arm, Fk);
        }

        float3 w = omega[i].xyz;
        w = w * damp + (tau * invI) * dt;
        omega[i] = (float4)(w, 0.0f);
        float4 dq = quat_from_vec(w * dt);
        quat[i] = normalize(quat_mul(dq, q));
    }
}
```

yes to clarify:

1) ad pneigh - actually it does not matter if we read from memory ports in local frame, or pneigh (ports in global frame) so I think we can gain nothing bu eliminating this as we would need ports and rotate them so there would be more global memory reads and more random, but it is true that this cots us few more writes of pneigh in  gather_port_forces, but that is coalesced (packed under on atom i)
2) npi refers to the number of neighbors, it is just bad name

---

now I would like to make the position-based projective dynamics (instead of force based) version of this.

Although you said then im projective dynamics we do not need to pass recoil forces to neighbors, I'm still not hundert percent convinced about that. Also what is a bit problematic for me - that we are doing rotation dynamics here, and I'm not sure if this is fully compactible with the idea that we solve some linear stiffness matrix problem. I mean rotation is not linear when the rotation angle is high, it can be only approximated as linear for small angles by 1st order taylor. So can we really say thais is Jacobi linar matrix solver?

So try to figure out how to make position based version of this solver with split between node and caping atoms (only node have rotation quaternions an fneighs) and with explicit recoils.

I still thin we should do the update in this framework

pi' = ( pi *alpha + sum_j fij ) /( alpha + sum_j Kij )
ti = sum_j torq(fij,pij); // torg
qi' = rotate( qi, ti ) /( alpha + sum_j Kijrij )

or mabybe we should do slerp (quaternion interpolation) because as I understand postion based dynamivs, it is basically interpolation (weighted average) between original position of particle pi weighted by the interia alpha=M/dt^2 and the position demanded by the contrains (bonds to neighbors)

so if I tranfer/generalize this to rotation we should to also interpolation but interpolation on surface of sphare (rotation manifold), which is slerp used for tracball rotation? Or not? What do you think?


what the fuck, there is again some confusion and overcomplication, what if the forward and backward neibhot. There is only one force between the ports (where atom i wants to have the neibhors) and the real position on nenihbor atoms. I do not see why you split it in two. I see this mistake did the GPT-5.2 while coding, and it produce disfuctional code (did not converged, did not conserved angular momentum), and too me quite some time to beat it and fore it to do it properly as I want whth the fneigh recoils. I think you have the same error in your metal picture. 

now you write explicit formula for the linear PD update for the postion of points (tranlation, linear position), but you do not write explict formula for the rotational update, which if bad, What I miss in your formula for rotation update (inconsistent with respect to translation update) is missing alpha*qi (current rotation) in the numerator of the interpolation. I man if you have p_new = (a*p_pred + sum K pj )/(a + sum K) why you do not have also for rotation something like 
q_new = (a*q_pred + sum K qj )/(a + sum K)
and instead you have just 
q_new = (sum K qj )/(a + sum K)

anyway, we should not use just plus and minis becasue the rotation is not really linear but we should use rotation, I think but we can try both and see which converge better and which preserve angular momentum better

---

## Kimi 2.5

 You are absolutely right about the **"one constraint"** principle. Splitting the interaction into forward/backward updates breaks Newton's third law and destroys momentum conservation. In PBD (Position-Based Dynamics), this manifests as asymmetric corrections: if you independently project the constraint from $i$'s side and then from $j$'s side without proper weighting, you create spurious forces. Your `fneigh`/`bkSlots` pattern is the correct way to handle this—you compute the constraint violation once, distribute the correction to both bodies according to their effective masses, and accumulate.

Here is the didactic breakdown of how to convert your force-based rigid-body solver to **Projective Dynamics** (PBD with a direct solver/Jacobi structure), while preserving momentum, handling rotations correctly, and avoiding the pitfalls you identified.

### 1. The Physics: Single Constraint, Distributed Correction

Consider one bond between a **port** on rigid body $i$ and a neighbor atom $j$.

*   **Port position:** $\mathbf{p}_i = \mathbf{x}_i + \mathbf{R}_i \mathbf{a}_i$ (where $\mathbf{a}_i$ is the local offset, $\mathbf{R}_i$ the rotation matrix from quaternion $\mathbf{q}_i$)
*   **Neighbor position:** $\mathbf{x}_j$
*   **Constraint:** $C = |\mathbf{p}_i - \mathbf{x}_j| - d_0 = 0$

When we project this constraint, we need to move $\mathbf{x}_i$, rotate $\mathbf{R}_i$, and move $\mathbf{x}_j$. The key is the **generalized inverse mass** $W$, which determines how much each degree of freedom moves.

**Generalized Inverse Masses:**
*   **Translation of $i$:** $W_{\text{trans},i} = 1/m_i$
*   **Translation of $j$:** $W_{\text{trans},j} = 1/m_j$  
*   **Rotation of $i$:** $W_{\text{rot},i} = (\mathbf{r}_i \times \mathbf{n})^T \mathbf{I}_i^{-1} (\mathbf{r}_i \times \mathbf{n})$  
    where $\mathbf{r}_i = \mathbf{p}_i - \mathbf{x}_i$ (the arm) and $\mathbf{n}$ is the constraint direction.

**Total compliance:** $W_{\text{total}} = W_{\text{trans},i} + W_{\text{trans},j} + W_{\text{rot},i}$

**The Corrections (Jacobi Update):**
For a violation $\Delta d = |\mathbf{p}_i - \mathbf{x}_j| - d_0$, the position-based correction is:
$$
\Delta \mathbf{x}_i = + \frac{W_{\text{trans},i}}{W_{\text{total}}} \Delta d \, \mathbf{n} \quad \text{(move COM toward target)}
$$
$$
\Delta \mathbf{x}_j = - \frac{W_{\text{trans},j}}{W_{\text{total}}} \Delta d \, \mathbf{n} \quad \text{(move neighbor away)}
$$
$$
\Delta \boldsymbol{\theta}_i = + \frac{W_{\text{rot},i}}{W_{\text{total}}} \Delta d \, \mathbf{I}_i^{-1} (\mathbf{r}_i \times \mathbf{n}) \quad \text{(rotate body)}
$$

**Why this conserves momentum:** The linear impulses $m_i \Delta \mathbf{x}_i$ and $m_j \Delta \mathbf{x}_j$ are equal and opposite. The angular impulse $\mathbf{I}_i \Delta \boldsymbol{\theta}_i$ exactly matches the moment of the linear impulse applied at the port ($\mathbf{r}_i \times (m_i \Delta \mathbf{x}_{\text{from rot}})$). This is preserved because we solve the constraint **once** and distribute, rather than solving twice.

### 2. Rotation Update: Why NO SLERP and NO $\alpha \mathbf{q}$

You asked why the formula isn't $\mathbf{q}_{\text{new}} = (\alpha \mathbf{q}_{\text{old}} + \sum K \mathbf{q}_{\text{target}}) / \dots$

**The fundamental reason:** PBD does **not** interpolate between orientations. It applies an **infinitesimal rotation** (a torque impulse) derived from the constraint violation. 

*   **Position space:** $p_{\text{new}} = p_{\text{old}} + \Delta \mathbf{x}$ (additive)
*   **Rotation space:** $\mathbf{q}_{\text{new}} = \delta \mathbf{q} \otimes \mathbf{q}_{\text{old}}$ (multiplicative, or "additive" in the tangent space)

SLERP is for blending between two known orientations (e.g., animation keyframes). Here, we have a **torque** (angular displacement vector $\Delta \boldsymbol{\theta}$) that we want to apply to the current orientation.

**The Correct Update:**
```opencl
// Delta_theta is the axis-angle vector (small rotation)
float3 dtheta = ...; // From constraint projection above

// Convert to quaternion delta: [axis * sin(|dtheta|/2), cos(|dtheta|/2)]
// For small angles: approximately [dtheta/2, 1]
float angle = length(dtheta);
if (angle > 1e-6f) {
    float4 dq = (float4)(normalize(dtheta) * sin(angle * 0.5f), cos(angle * 0.5f));
    q[i] = normalize(quat_mul(dq, q[i]));  // Apply on the left (local rotation)
}
```

**Why no $\alpha \mathbf{q}$ term:** In your position formula $p' = (\alpha p + \sum K p_{\text{target}})/(\alpha + \sum K)$, the $\alpha p$ term represents the **inertia** (resistance to change). In rotation, this resistance is embedded in $W_{\text{rot},i}$ (the moment of inertia). The "interpolation" happens in the **calculation of $\Delta \boldsymbol{\theta}$** via the weight $W_{\text{rot},i} / W_{\text{total}}$. You don't interpolate the quaternions; you rotate the current quaternion by a weighted amount.

### 3. The OpenCL Implementation Structure

This mirrors your force-based design but exchanges forces for **position corrections**.

#### Kernel 1: Compute Corrections (Replaces `gather_port_forces`)
This runs only on **nodes** ($i < \text{nnode}$). For each bond, it computes the correction needed and stores the **neighbor's portion** in `delta_neigh` (equivalent to your `fneigh` but storing position deltas).

```opencl
__kernel void compute_corrections(
    const int nnode,
    __global const float4* pos,        // w = invMass
    __global const float4* quat,
    __global const int4*   neighs,
    __global const float4* port_local, // Local offsets (xyz, 0)
    __global const float*  rest_len,
    __global const float*  stiffness,  // Or compliance alpha
    __global float4* dpos_node,        // Accumulated COM correction for node i
    __global float4* drot_node,        // Accumulated rotation vector for node i
    __global float4* dpos_neigh        // Correction to write to neighbor's slot (recoils)
) {
    int i = get_global_id(0);
    if (i >= nnode) return;

    float3 xi = pos[i].xyz;
    float4 qi = quat[i];
    float invMi = pos[i].w;
    float mi = (invMi > 1e-12f) ? 1.0f/invMi : 1e12f;
    float3 invIi_vec = (float3)(2.5f * invMi); // Simplified: 1/(0.4*m) for sphere

    float3 dx_accum = (float3)(0.0f);
    float3 dtheta_accum = (float3)(0.0f);
    int i4 = i * 4;

    for (int k = 0; k < 4; k++) {
        int j = neighs[i][k];
        if (j < 0) break;

        float3 aj = port_local[i4 + k].xyz;
        float3 ri = quat_rotate(qi, aj); // Arm in world space
        float3 pi = xi + ri;             // Port position
        
        float3 xj = pos[j].xyz;
        float invMj = pos[j].w;
        
        float3 diff = pi - xj;
        float dist = length(diff);
        if (dist < 1e-6f) continue;
        
        float3 n = diff / dist;
        float C = dist - rest_len[i4 + k];
        float K = stiffness[i4 + k]; // This acts as the 'w' weight in XPBD
        
        // Generalized inverse masses
        float W_lin_i = invMi;
        float W_lin_j = invMj;
        
        float3 rxn = cross(ri, n);
        // W_rot = (r x n)^T * I^-1 * (r x n). Using diagonal approx:
        float W_rot_i = dot(rxn*rxn, invIi_vec); 
        
        float W_total = W_lin_i + W_lin_j + W_rot_i + 1e-12f;
        
        // Lagrange multiplier (signed)
        float delta_lambda = -C / W_total;
        
        // --- Correction for node i (accumulated) ---
        // Translational part
        dx_accum += (W_lin_i / W_total) * delta_lambda * n * dist; 
        // Note: above is simplified; strictly it's delta_lambda * n * (1/mi factor already in W)
        // Actually cleaner:
        float3 dx_i = (invMi / W_total) * -C * n;
        dx_accum += dx_i;
        
        // Rotational part: dtheta = I^-1 * (r x n) * scale
        float3 dtheta_i = (invIi_vec * rxn) * (W_rot_i / W_total) * delta_lambda * dist;
        // Or directly:
        dtheta_i = invIi_vec * rxn * (-C / W_total);
        dtheta_accum += dtheta_i;
        
        // --- Correction for neighbor j (recoil) ---
        // This is stored in the slot so j can pick it up
        float3 dx_j = -(invMj / W_total) * -C * n; // Moves opposite to n
        dpos_neigh[i4 + k] = (float4)(dx_j, 0.0f);
    }
    
    dpos_node[i] = (float4)(dx_accum, 0.0f);
    drot_node[i] = (float4)(dtheta_accum, 0.0f);
}
```

#### Kernel 2: Apply Corrections (Replaces `integrate_rigid_explicit`)
This runs on **all atoms**. It gathers the recoils via `bkSlots` (exactly like your force gathering), then updates positions and quaternions. There is no explicit velocity here (PBD is direct position update), but you can update velocities post-projection via `v = (x_new - x_old)/dt` if needed for output.

```opencl
__kernel void apply_corrections(
    const int natoms,
    const int nnode,
    __global float4* pos,
    __global float4* quat,
    __global const int4* bkSlots,      // Back-mapping: which slots affect me?
    __global const float4* dpos_node,  // From kernel 1 (nodes only)
    __global const float4* drot_node,  // From kernel 1 (nodes only)
    __global const float4* dpos_neigh, // From kernel 1 (all bonds)
    const float relaxation             // 0.2 - 0.8 for Jacobi stability
) {
    int i = get_global_id(0);
    if (i >= natoms) return;

    float3 dx = (float3)(0.0f);
    
    // If node, start with its own COM correction from its ports
    if (i < nnode) {
        dx += dpos_node[i].xyz;
    }
    
    // Gather recoils from all neighbors that wrote to my slots
    int4 bk = bkSlots[i];
    for (int k = 0; k < 4; k++) {
        int slot = bk[k];
        if (slot >= 0) {
            dx += dpos_neigh[slot].xyz;
        }
    }
    
    // Apply position update (Jacobi step with under-relaxation)
    float3 xi = pos[i].xyz;
    xi += dx * relaxation;
    pos[i].xyz = xi;
    
    // Apply rotation update (only nodes)
    if (i < nnode) {
        float3 dtheta = drot_node[i].xyz * relaxation;
        
        float angle = length(dtheta);
        if (angle > 1e-8f) {
            float3 axis = dtheta / angle;
            float4 dq = (float4)(axis * sin(angle * 0.5f), cos(angle * 0.5f));
            float4 qi = quat[i];
            // Apply rotation: q_new = dq * q_old (local rotation applied to current)
            quat[i] = normalize(quat_mul(dq, qi));
        }
    }
}
```

### 4. Why This Conserves Momentum

*   **Linear:** For each constraint $i$-$j$, the code computes $\Delta \mathbf{x}_i$ and $\Delta \mathbf{x}_j$ such that $m_i \Delta \mathbf{x}_i = -m_j \Delta \mathbf{x}_j$ (equal and opposite impulses). When you sum over the system in `apply_corrections`, the internal forces cancel, and only external constraints (anchors) change total momentum.
*   **Angular:** The term $\Delta \boldsymbol{\theta}_i = \mathbf{I}_i^{-1} (\mathbf{r}_i \times \mathbf{n}) \cdot \text{scale}$ ensures that the angular momentum change $\mathbf{I}_i \Delta \boldsymbol{\theta}_i$ equals the moment of the linear correction applied at the port ($\mathbf{r}_i \times (m_i \Delta \mathbf{x}_{i,\text{rot}})$). Since the neighbor feels an equal and opposite linear impulse at the same point in space (the port location), the total angular momentum is conserved.

### 5. Regarding Auxiliary Particles vs. Orientation Vectors

If you use an **orientation vector** (Cosserat/director) instead of a dummy particle:
*   You store a unit vector $\mathbf{d}$ (the pi-orbital direction) instead of a position $\mathbf{x}_{\text{dummy}}$.
*   The constraint becomes: $\mathbf{d} \cdot (\mathbf{p}_1 - \mathbf{p}_2) = 0$ (orthogonality to bond) or similar.
*   The projection updates $\mathbf{d}$ via parallel transport: $\mathbf{d}' = \text{normalize}(\mathbf{d} + \Delta \boldsymbol{\theta} \times \mathbf{d})$.
*   **Advantage:** No extra position DOF, no risk of the dummy drifting away.
*   **Disadvantage:** You must implement the constraint gradient on the manifold $S^2$ (slightly more complex cross-product logic).

However, your **port-based rigid body approach** is superior to both dummy particles and Cosserat rods for this specific problem because:
1.  It avoids the non-physical inertia of a dummy particle.
2.  It avoids the drift/numerical normalization issues of Cosserat directors.
3.  It directly maps to your existing force-field architecture (bonds from ports).

### 6. Final Notes on Convergence

*   **Jacobi vs Gauss-Seidel:** The scheme above is Jacobi (updates based on state at iteration start). It requires `relaxation < 1.0` (typically 0.2-0.5) to avoid divergence.
*   **Rotation Linearization:** You were worried about rotation non-linearity. PBD handles this by iterating. Each iteration linearizes the rotation (small angle approx). For large rotations, you need more iterations, but the scheme remains stable because the quaternion re-normalization projects back to the manifold exactly.
*   **Efficiency:** The memory traffic pattern (`dpos_neigh` writes in kernel 1, reads via `bkSlots` in kernel 2) is identical to your force-based solver and is coalesced (each node writes 4 float4s contiguously).

This design preserves your `fneigh`/`bkSlots` philosophy, conserves momentum exactly, and correctly handles the non-Euclidean nature of rotations without expensive SLERP operations.

---

## Caude-4.5-Sonet

I'll analyze this code systematically and help you develop a position-based dynamics version.

## Current Force-Based Implementation Analysis

Your current design is quite clever. Let me break it down:

### Architecture
1. **Split atom types**: Node atoms (with rotation) vs capping atoms (point particles)
2. **Port-based bonding**: Each node has up to 4 ports in local frame, rotated to world frame
3. **Explicit recoil tracking**: `fneigh` stores forces that neighbors exert back on the central atom

### Physical Correctness

**Linear momentum**: ✅ Conserved
- For each bond i→j: force `f` on neighbor j, force `-f` stored in `fneigh[i*4+k]`
- In integration, atom j gathers this `-f` via `bkSlots` mapping
- Net force on system = 0

**Angular momentum**: ✅ Conserved
- Torques computed from port positions (`r_arm`) and recoil forces
- Each action-reaction pair creates equal-opposite torques about CoM

### Efficiency Observations

**Good**:
- Coalesced writes to `fneigh`, `pneigh` (packed per atom i)
- `bkSlots` precomputed on host (one-time cost)
- Only node atoms do rotation math

**Potential issues**:
- `pneigh` write/read could be eliminated (as you noted) but probably not worth complexity
- Random access in `bkSlots` gathering (unavoidable)

---

## Position-Based Dynamics Version

You're absolutely right to be concerned about rotations in PBD. Here's the key insight:

### The Rotation Manifold Problem

Standard PBD solves: **minimize** `½α||p' - p||² + ½∑Kᵢⱼ||p' - constraint||²`

For rotations, we need: **minimize** `½α·d(q', q)² + ½∑Kᵢⱼ·d(q', target)²`

where `d(q₁, q₂)` is the **geodesic distance** on SO(3), not Euclidean distance!

### Your Interpolation Intuition is Correct

For positions: `p_new = (α·p_old + ∑Kⱼ·p_target) / (α + ∑Kⱼ)`

For rotations, the analog is **NOT** `q_new = (α·q_old + ∑Kⱼ·q_target) / (α + ∑Kⱼ)` (this doesn't preserve unit norm!)

Instead: **weighted SLERP** or exponential map averaging.

### Proposed PBD Implementation

```python
# Kernel 1: Compute constraint violations (no forces, just deltas)
__kernel void compute_constraints(
    const int nnode,
    __global const float4* pos,
    __global const float4* quat,
    __global int4* neighs,
    __global float4* bKs,
    __global const float4* port_local,
    __global const uchar* port_n,
    __global float4* delta_p,      // translation corrections
    __global float4* delta_w,      // angular velocity corrections (tangent space)
    __global float4* pneigh        // cached port positions
) {
    int i = get_global_id(0);
    if (i >= nnode) return;
    
    float3 p_i = pos[i].xyz;
    float4 q_i = quat[i];
    
    int4 ng = neighs[i];
    int* neighbors = (int*)&ng;
    int npi = (int)port_n[i];
    
    float3 dp_i = (float3)(0.0f);
    float3 dw_i = (float3)(0.0f);  // torque accumulator
    
    int i4 = i * 4;
    for (int k = 0; k < 4; k++) {
        int j = neighbors[k];
        if (j < 0) break;
        
        float Kij = fetch4(bKs, i, k);
        float3 p_j = pos[j].xyz;
        
        // Rotate port to world frame
        float3 r_arm = (float3)(0.0f);
        if (k < npi) r_arm = quat_rotate(q_i, port_local[i4 + k].xyz);
        
        float3 tip = p_i + r_arm;
        float3 violation = p_j - tip;  // how far is neighbor from port?
        
        // Accumulate translation correction (weighted by stiffness)
        dp_i += violation * Kij;
        
        // Accumulate torque correction
        dw_i += cross(r_arm, violation) * Kij;
        
        // Cache for neighbor updates
        pneigh[i4 + k] = (float4)(r_arm, 0.0f);
    }
    
    delta_p[i] = (float4)(dp_i, 0.0f);
    delta_w[i] = (float4)(dw_i, 0.0f);
}

# Kernel 2: Apply corrections with inertia
__kernel void apply_corrections(
    const int natoms,
    const int nnode,
    __global float4* pos,
    __global float4* pos_prev,
    __global float4* quat,
    __global float4* quat_prev,
    __global int4* bkSlots,
    __global float4* delta_p,
    __global float4* delta_w,
    __global float4* pneigh,
    const float dt,
    const float alpha_scale  // = 1/(compliance * dt^2)
) {
    int i = get_global_id(0);
    if (i >= natoms) return;
    
    float invMass = pos[i].w;
    float mass = (invMass > 1e-12f) ? (1.0f / invMass) : 1e12f;
    float alpha_p = mass / (dt * dt);  // positional inertia
    
    // --- TRANSLATION UPDATE ---
    float3 p_old = pos[i].xyz;
    float3 dp = delta_p[i].xyz;
    
    // Gather recoil corrections from neighbors
    int4 bk = bkSlots[i];
    int* islots = (int*)&bk;
    float total_K = 0.0f;
    for (int k = 0; k < 4; k++) {
        int islot = islots[k];
        if (islot >= 0) {
            int i_node = islot / 4;
            int k_port = islot % 4;
            // The neighbor pushed us, we push back
            float3 r_arm = pneigh[islot].xyz;
            float3 violation = p_old - (pos[i_node].xyz + r_arm);
            float Kij = /* need to store this */;  // TODO: cache K values
            dp += violation * Kij;
            total_K += Kij;
        }
    }
    
    // PBD update: interpolate between inertia and constraint
    float3 p_new = (alpha_p * p_old + dp) / (alpha_p + total_K + 1e-12f);
    pos[i].xyz = p_new;
    
    // --- ROTATION UPDATE (only for nodes) ---
    if (i < nnode) {
        float I = 0.4f * mass;  // moment of inertia
        float alpha_r = I / (dt * dt);  // rotational inertia
        
        float3 dw = delta_w[i].xyz;
        float4 q_old = quat[i];
        
        // Method 1: Small angle approximation (fast but less accurate)
        // float3 delta_angle = dw / (alpha_r + total_K + 1e-12f);
        // float4 dq = quat_from_vec(delta_angle);
        // quat[i] = normalize(quat_mul(dq, q_old));
        
        // Method 2: Proper geodesic interpolation (recommended)
        // Convert torque to target rotation
        float3 axis = normalize(dw + (float3)(1e-12f));
        float angle = length(dw) / (alpha_r + total_K + 1e-12f);
        float4 q_target = quat_from_axis_angle(axis, angle);
        float4 q_delta = quat_mul(q_target, quat_conj(q_old));
        
        // SLERP with inertia weight
        float t = total_K / (alpha_r + total_K + 1e-12f);
        quat[i] = quat_slerp(q_old, quat_mul(q_delta, q_old), t);
    }
}
```

### Key Differences from Your Force Version

1. **No `fneigh` storage** - we compute violations, not forces
2. **Direct position updates** - no velocity integration
3. **Proper rotation interpolation** - SLERP or exponential map, not linear addition
4. **Missing α·q_old term** - You're right! It should be:
   ```
   q_new = SLERP(q_old, q_target, K/(α+K))
   ```
   Where `α` is rotational inertia. When `K→∞`, `q_new→q_target` (stiff constraint). When `α→∞`, `q_new→q_old` (heavy particle).

### Momentum Conservation in PBD

**Your concern is valid**: Standard PBD Jacobi iterations do NOT conserve momentum per-iteration. However:

- **Gauss-Seidel PBD**: Conserves momentum if you update recoils immediately
- **Extended PBD (XPBD)**: Adds compliance, recovers momentum conservation in limit
- **Your explicit recoil approach**: Should work but needs careful derivation

For your case, I'd recommend:
1. Use **symmetric updates**: When updating atom i, immediately update its neighbors
2. Or use **multi-pass**: Gather phase (compute all deltas) → Scatter phase (apply symmetrically)

Would you like me to work out the full symmetric update version with proper momentum conservation guarantees?

---

## USER

define quat_slerp

explain in detail this

* Gauss-Seidel PBD: Conserves momentum if you update recoils immediately
* Extended PBD (XPBD): Adds compliance, recovers momentum conservation in limit

ad "Use symmetric updates: When updating atom i, immediately update its neighbors"
* no this is exactly what I do not want to do because that lead to race-condigtion or require atomic writes at random places which is slow on GPU
ad "Or use multi-pass: Gather phase (compute all deltas) → Scatter phase (apply symmetrically)"
* yes, that is what I'm doing or not? My strategy is all about how to make this properly so that it conserve the momenrum and use only gather operations.

---

## Claude-4.5-Sonet

## Quaternion SLERP Definition

```c
float4 quat_slerp(float4 q0, float4 q1, float t) {
    // Spherical Linear intERPolation
    // Interpolates along great circle on 4D unit sphere
    
    float dot_prod = dot(q0, q1);
    
    // If quaternions are nearly opposite, flip one
    if (dot_prod < 0.0f) {
        q1 = -q1;
        dot_prod = -dot_prod;
    }
    
    // If very close, use linear interpolation (avoid division by ~0)
    if (dot_prod > 0.9995f) {
        return normalize(q0 + t * (q1 - q0));
    }
    
    // Spherical interpolation
    float theta = acos(dot_prod);
    float sin_theta = sin(theta);
    float w0 = sin((1.0f - t) * theta) / sin_theta;
    float w1 = sin(t * theta) / sin_theta;
    
    return w0 * q0 + w1 * q1;
}
```

**Why SLERP for rotations?**
- Linear interpolation `(1-t)·q0 + t·q1` doesn't stay on unit sphere
- SLERP maintains constant angular velocity along geodesic
- For rotations: interpolating **angles**, not coordinates

## Momentum Conservation in PBD - Detailed Explanation

### Standard Jacobi PBD: Does NOT Conserve Momentum

```
for each constraint (i,j):
    δp_i = compute_correction_for_i()
    δp_j = compute_correction_for_j()
    
for each particle i:
    p_i += δp_i  // sum all corrections for i
```

**Problem**: When particle `i` is updated based on constraint with `j`, the correction `δp_i` accounts for the current position of `j`. But when `j` is updated later in the same iteration, it moves based on the **old** position of `i`. The action-reaction symmetry is broken **within** the iteration.

**Result**: `δp_i + δp_j ≠ 0` in general → momentum drift per iteration (though it may converge to conserving solution after many iterations).

### Gauss-Seidel PBD: CAN Conserve Momentum

```
for each constraint (i,j) in sequential order:
    // Compute corrections based on CURRENT positions
    violation = compute_constraint_violation(p_i, p_j)
    
    // Update BOTH particles IMMEDIATELY
    δp_i = -w_i * violation  // w_i = 1/(m_i * (1/m_i + 1/m_j))
    δp_j = +w_j * violation  // w_j = 1/(m_j * (1/m_i + 1/m_j))
    
    p_i += δp_i  // update NOW
    p_j += δp_j  // update NOW
    
    // Next constraint sees updated positions
```

**Key difference**: Immediate update means:
- When constraint `(i,j)` is processed: `δp_i + δp_j = 0` (by construction)
- When constraint `(i,k)` is processed next: it sees the **new** `p_i`, maintaining consistency

**Momentum conservation**: Each constraint update preserves momentum **exactly** because updates are symmetric and immediate.

**Why you can't use this on GPU**: Race conditions! If constraints `(i,j)` and `(i,k)` both update `p_i` simultaneously → undefined behavior.

### Extended PBD (XPBD): The Compliance Trick

XPBD adds **Lagrange multipliers** that track constraint forces across iterations:

```
Initialize: λ = 0 for each constraint

for iteration in range(num_iterations):
    for each constraint (i,j):
        violation = C(p_i, p_j)  // constraint function
        
        // Compute Lagrange multiplier update
        Δλ = -violation / (w_i + w_j + α/dt²)  // α = compliance
        λ += Δλ
        
        // Apply correction
        δp_i = -Δλ * w_i * ∇C_i
        δp_j = -Δλ * w_j * ∇C_j
        
        p_i += δp_i
        p_j += δp_j
```

**The magic**: 
- `α = 0`: infinitely stiff constraint (standard PBD)
- `α > 0`: soft constraint (like spring with stiffness `k = 1/α`)
- `λ` accumulates the "total force" applied by constraint across iterations

**Momentum conservation**: 
- Over **multiple iterations**, `λ` converges so that `∑δp_i = 0` globally
- In the limit of convergence: `λ = λ_true` (the actual constraint force)
- The accumulated corrections become equivalent to solving the implicit system exactly

**Why it helps**: Even with Jacobi updates, XPBD's Lagrange multipliers "remember" past corrections, allowing momentum to equilibrate across iterations.

## Your Gather-Only Strategy: The Challenge

You want:
```
Phase 1 (Gather): Each atom i computes what IT wants
    δp_i = f(p_i, p_j, p_k, ...)  // based on neighbors

Phase 2 (Apply): Each atom updates independently
    p_i += δp_i
```

**The fundamental issue**: 

Consider bond between atoms `i` and `j`:
- Atom `i` computes: `δp_i = Δ·(p_j - p_i) / (m_i + m_j)` (pulls toward j)
- Atom `j` computes: `δp_j = Δ·(p_i - p_j) / (m_i + m_j)` (pulls toward i)

For momentum conservation: `m_i·δp_i + m_j·δp_j = 0`

This requires: `δp_i = -(m_j/m_i)·δp_j`

But if both are computed **independently** in parallel, they're both based on the **same** positions `(p_i, p_j)`, so the symmetry works!

**Wait, does standard Jacobi conserve momentum then?**

Let me recalculate:

```
Constraint: |p_j - p_i| = L (bond length)
Current violation: C = |p_j - p_i| - L

Standard PBD correction (symmetric masses):
δp_i = -0.5 * C * n̂  // n̂ = (p_j - p_i)/|p_j - p_i|
δp_j = +0.5 * C * n̂

Check: δp_i + δp_j = 0 ✓
```

**So Jacobi DOES conserve momentum for pairwise constraints!**

The issue arises with:
1. **Multi-way constraints** (angle, dihedral) - harder to make symmetric
2. **Rotations** - torques need careful handling
3. **Convergence** - Jacobi is slower than Gauss-Seidel

## Your Port-Based System: The Rotation Complication

Here's the tricky part:

```
Atom i (node) has ports at positions: p_i + R(q_i)·r_local_k
Atom j (neighbor) is at: p_j

Constraint: atom j should be at port k of atom i
Violation: v_k = p_j - (p_i + R(q_i)·r_local_k)
```

**Updates needed**:
- Translation of i: `δp_i`
- Rotation of i: `δq_i` (or `δω_i`)
- Translation of j: `δp_j`

**For momentum conservation**:
```
Linear: m_i·δp_i + m_j·δp_j = 0
Angular: I_i·δω_i + r_k × (m_j·δp_j) = 0  // torque balance about CoM
```

**Your current approach (force-based)**:
```
f_k = K·v_k  // force at port k

On atom i:
  F_i = ∑_k f_k                    // net force
  τ_i = ∑_k (r_k × f_k)            // net torque
  
On atom j:
  F_j = -f_k  // stored in fneigh, gathered via bkSlots
  
Update:
  v_i += F_i/m_i · dt
  ω_i += τ_i/I_i · dt
  v_j += F_j/m_j · dt
```

This conserves momentum because `F_i = -∑F_j` by construction (every `f_k` has recoil `-f_k`).

## Position-Based Version with Gather-Only

**The key insight**: You need to compute corrections such that:

```
m_i·δp_i + ∑_{j∈neighbors} m_j·δp_j = 0
I_i·δω_i + ∑_{j} (r_k × m_j·δp_j) = 0
```

**Strategy**: Treat each port constraint independently (like your force version), compute symmetric corrections:

```c
__kernel void compute_and_apply_corrections(
    // ... parameters ...
) {
    int i = get_global_id(0);
    if (i >= nnode) return;
    
    float3 p_i = pos[i].xyz;
    float4 q_i = quat[i];
    float m_i = 1.0f / pos[i].w;
    float I_i = 0.4f * m_i;
    float alpha_p = m_i / (dt * dt);
    float alpha_r = I_i / (dt * dt);
    
    float3 dp_i_total = (float3)(0.0f);
    float3 dw_i_total = (float3)(0.0f);
    
    // Process each port constraint
    for (int k = 0; k < 4; k++) {
        int j = neighs[i][k];
        if (j < 0) break;
        
        float K_ij = bKs[i][k];
        float m_j = 1.0f / pos[j].w;
        
        // Port position in world frame
        float3 r_k = quat_rotate(q_i, port_local[i*4 + k].xyz);
        float3 port_world = p_i + r_k;
        
        // Constraint violation
        float3 v_k = pos[j].xyz - port_world;
        float C = length(v_k);
        if (C < 1e-8f) continue;
        float3 n = v_k / C;
        
        // Compute weights (effective mass accounting for rotation)
        // This is the TRICKY part!
        float3 dr_dq = cross(r_k, n);  // how rotation affects constraint
        float w_i_trans = 1.0f / m_i;
        float w_i_rot = dot(dr_dq, dr_dq) / I_i;
        float w_j = 1.0f / m_j;
        
        float w_total = w_i_trans + w_i_rot + w_j + 1.0f/(K_ij + 1e-12f);
        
        // Corrections (PBD style)
        float λ = -C / w_total;  // Lagrange multiplier
        
        float3 dp_i = λ * n * w_i_trans;
        float3 dw_i = λ * cross(r_k, n) / I_i;  // angular impulse
        float3 dp_j = -λ * n * w_j;
        
        // Accumulate corrections for atom i
        dp_i_total += dp_i;
        dw_i_total += dw_i;
        
        // Store correction for atom j (to be gathered later)
        delta_neigh[i*4 + k] = (float4)(dp_j, 0.0f);
    }
    
    // Apply corrections to atom i
    pos_delta[i] = (float4)(dp_i_total, 0.0f);
    omega_delta[i] = (float4)(dw_i_total, 0.0f);
}

__kernel void gather_and_apply(
    // ... parameters ...
) {
    int i = get_global_id(0);
    if (i >= natoms) return;
    
    float3 dp = pos_delta[i].xyz;
    
    // Gather recoil corrections from neighbors
    int4 bk = bkSlots[i];
    for (int k = 0; k < 4; k++) {
        int islot = bk[k];
        if (islot >= 0) {
            dp += delta_neigh[islot].xyz;
        }
    }
    
    // Apply
    pos[i].xyz += dp;
    
    if (i < nnode) {
        float3 dw = omega_delta[i].xyz;
        float angle = length(dw);
        if (angle > 1e-8f) {
            float3 axis = dw / angle;
            float4 dq = quat_from_axis_angle(axis, angle);
            quat[i] = normalize(quat_mul(dq, quat[i]));
        }
    }
}
```

**Does this conserve momentum?**

Check for port k connecting atoms i and j:
```
m_i·dp_i + m_j·dp_j 
= m_i·(λ·n/m_i) + m_j·(-λ·n/m_j)
= λ·n - λ·n = 0 ✓

I_i·dω_i + r_k × (m_j·dp_j)
= I_i·(λ·cross(r_k,n)/I_i) + r_k × (m_j·(-λ·n/m_j))
= λ·cross(r_k,n) - λ·cross(r_k,n) = 0 ✓
```

**Yes! Momentum conserved per constraint, per iteration!**

The key was computing `λ` using the **generalized mass** `w_total` that accounts for both translation AND rotation coupling.

Does this clarify the approach?