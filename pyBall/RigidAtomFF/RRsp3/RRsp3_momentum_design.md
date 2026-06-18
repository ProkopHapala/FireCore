# RRsp3: Momentum-Corrected Rigid Body Dynamics Design

## 1. Executive Summary
This document details the design and implementation of a momentum-conserving, accelerated solver for the `RRsp3` rigid body molecular dynamics module. The system uses **Extended Position Based Dynamics (XPBD)** with **Heavy-Ball Momentum** (solver acceleration) to ensure fast convergence while maintaining physical correctness (linear and angular momentum conservation).

## 2. Problem Analysis

### 2.1 Conservation Laws in Iterative Solvers
*   **Linear Momentum**: In PBD/XPBD, linear momentum is conserved if every positional correction $\Delta \mathbf{x}_{ij}$ applied to particle $i$ due to a constraint with $j$ is matched by an equal and opposite correction $\Delta \mathbf{x}_{ji} = -\Delta \mathbf{x}_{ij}$ to $j$.
    *   *Implementation*: The `RRsp3` solver achieves this by computing the constraint impulse once, applying it to $i$, and writing the recoil to $j$'s accumulator (`dpos_neigh`). The accumulation phase sums these symmetrically.
*   **Angular Momentum**: 
    *   Rigid bodies must rotate in response to off-center impulses. The update $\Delta \boldsymbol{\theta} \propto \mathbf{I}^{-1} (\mathbf{r} \times \mathbf{p})$ ensures that the change in angular momentum equals the torque applied.
    *   *Issue*: Exact conservation requires that the sum of torques vanishes. In PBD, this holds if the constraint gradients are symmetric. The current implementation uses a linearized angular update which is accurate for small steps but may exhibit drift for large rotations if not carefully regularized.

### 2.2 Projective Dynamics (PD) vs XPBD
The user noted confusion regarding $M/\Delta t^2$ vs inverse mass weighting.
*   **Projective Dynamics**: Formulates the problem as an energy minimization:
    $$ E(\mathbf{x}) = \frac{1}{2 \Delta t^2} (\mathbf{x} - \mathbf{y})^T \mathbf{M} (\mathbf{x} - \mathbf{y}) + \sum W_c(\mathbf{x}) $$
    The solver step involves solving a linear system $(\mathbf{M}/\Delta t^2 + \mathbf{L}) \Delta \mathbf{x} = \text{forces}$. The diagonal term is $\mathbf{M}/\Delta t^2$.
*   **XPBD**: Formulates constraints $C(\mathbf{x})=0$ with compliance $\alpha = 1/(K \Delta t^2)$.
    The positional update is:
    $$ \Delta \mathbf{x} = \frac{-C(\mathbf{x})}{\alpha + \sum w_i \|\nabla C\|^2} w_i \nabla C $$
    where $w_i = 1/m_i$.
    *   *Equivalence*: Dividing the PD equation by $M/\Delta t^2$ roughly yields the XPBD form. The term $w_i$ in XPBD corresponds to the inverse of the inertial term in PD. 
    *   *Decision*: `RRsp3` uses the **XPBD formulation**. It is explicit, easier to implement on GPU (no global matrix solve), and handles non-linear constraints (like collisions) naturally. The compliance $\alpha = 1/(K \Delta t^2)$ correctly regularizes the constraint stiffness against the inertia $M/\Delta t^2$.

### 2.3 Solver Momentum vs Physical Momentum
*   **Physical Momentum**: The velocity $\mathbf{v}$ carried between time steps. Conserved by the symplectic Euler integration (or Verlet).
*   **Solver Momentum (Heavy-Ball)**: An acceleration technique for the Jacobi iterations *within* a time step.
    $$ \mathbf{x}_{k+1} = \text{Jacobi}(\mathbf{x}_k) + \beta (\mathbf{x}_k - \mathbf{x}_{k-1}) $$
    This is purely algorithmic and does not affect the physical momentum of the system, but speeds up convergence to the constraint manifold.

## 3. Implementation Details

### 3.1 Kernel: `apply_corrections_rigid_ports`
The core logic was updated to support Heavy-Ball momentum for both translation and rotation.

**Translation Update:**
```c
float3 d_mom = dpos_mom[i].xyz;       // Previous step delta
float3 dx_total = dx_coll + dx_port;  // Constraint correction
float3 move = dx_total * relaxation + d_mom * beta;

pos[i].xyz += move;
dpos_mom[i] = (float4)(move, 0.0f);   // Store for next step
```

**Rotation Update:**
Rotation requires care because quaternions live on a hypersphere. We use a linear approximation for the momentum term, which is valid for the small adjustments typical in PBD iterations, followed by renormalization.
```c
// 1. Apply Jacobi Correction
float3 dtheta = drot_node[i].xyz * relaxation;
float4 q_jacobi = quat_rotate(q_old, dtheta);

// 2. Apply Solver Momentum
float4 dq_mom = dquat_mom[i];         // Previous quaternion delta
float4 q_new = q_jacobi + dq_mom * beta;
q_new = normalize(q_new);

// 3. Update State
quat[i] = q_new;
dquat_mom[i] = q_new - q_old;
```

### 3.2 Python Wrapper (`RRsp3.py`)
*   Added persistent buffers `cl_dpos_mom` and `cl_dquat_mom` to store the step differences between iterations.
*   Added `reset_momentum()` to clear these buffers at the start of a new physics time step (crucial, as solver momentum should not persist across time steps).

## 4. Testing & Verification

### 4.1 Convergence Test (`test_RRsp3_convergence.py`)
A headless batch script was created to rigorously test convergence.
*   **Setup**: Loads `backbone_pasivated-H.xyz`.
*   **Distortion**: Applies random noise to atom positions.
*   **Relaxation**: Runs 100 iterations of the Jacobi solver with varying $\beta$ (momentum) parameters: 0.0, 0.5, 0.8, 0.9.
*   **Metric**: Calculates mean constraint violation (bond length error) at each step.
*   **Result**: 
    *   $\beta=0.0$ (Standard Jacobi): Slow, monotonic convergence.
    *   $\beta=0.9$ (High Momentum): significantly faster convergence (orders of magnitude lower error for the same iteration count), potentially with slight oscillations.
*   **Output**: Generates `convergence_noise.png` (Log-Linear plot) and saves the final trajectory.

### 4.2 Reference comparison
The implementation aligns with `cpp/common/math/ProjectiveDynamics_d.cpp` in terms of the momentum update structure ($x_{new} = x_{corrected} + \beta \Delta x_{prev}$), adapted for the parallel GPU architecture of `RRsp3`.

## 5. Conclusion
The `RRsp3` module now possesses a fully functional, momentum-accelerated rigid body solver. The XPBD formulation ensures stability and physical plausibility (via compliance), while the Heavy-Ball momentum drastically reduces the computational cost (iterations) required to reach a given error tolerance.

---

## 6. Port Kernel Approaches: Pipeline, Tradeoffs, and Conservation

### 6.1 Common Solver Loop (all methods)
Every solver step executes the same broad-phase and collision stack; only the **port kernel** changes:

1. `update_bboxes_rigid` — reads `pos`, writes `bboxes_min/max`
2. `build_local_topology_rigid` — reads `pos`, `neighs`, `excl`, writes `neighs_local`, `excl_local`
3. `compute_collision_cluster_rigid` — reads `pos`, `radius`, `excl_local`, writes `dpos_coll`
4. **PORT KERNEL** (method-dependent) — writes `dpos_node`, `drot_node`, `dpos_neigh`
5. `apply_corrections_rigid_ports` — reads `pos`, `quat`, `dpos_node`, `drot_node`, `dpos_neigh`, `dpos_coll`, writes updated `pos`, `quat`, `dpos_mom`, `dquat_mom` (and `tips` for massless modes)

The remainder of this section describes each port kernel in detail, split by the fundamental distinction: **massfull** (physical rotational inertia) vs **massless** (pure geometric alignment).

---

### 6.2 Massfull Rotation (Physical Inertial DOF)
In this class, the quaternion is a true mechanical degree of freedom with inertia and angular momentum. Off-center port forces generate torque, which is absorbed by the rotational DOF, so total angular momentum `L_total = L_trans + L_spin` is conserved.

#### 6.2.1 `orig` / `current` — Full Rigid-Body XPBD
- **Rotation model**: Physical rotational inertia. Both linear and angular impulses are solved simultaneously.
- **Strategy**: Each node rotates its own ports into world space using its quaternion, then measures the distance from its tip to the **neighbor atom center**. The XPBD impulse is distributed into linear (`dpos_node`) and angular (`drot_node`) recoil using full rigid-body compliance (`w_ang` term).
- **Pipeline**: single kernel `compute_ports_cluster_rigid_orig`
  - **Consumes**: `pos`, `quat`, `port_local`, `stiffness`
  - **Produces**: `dpos_node`, `drot_node`, `dpos_neigh`
- **Conservation**: Exact `P` and `L_total = L_trans + L_spin`. The impulse is off-center (direction `n = (xj - tip_i)/|...|`), but the resulting torque is absorbed by the physical rotational DOF, so total angular momentum is conserved.
- **Tradeoff**: Most physically complete model; requires tuning the mass–inertia ratio and may need smaller time steps for stiff rotations.

---

### 6.3 Massless Rotation (Geometric Alignment, Zero Inertia)
In this class, the quaternion is **not** a mechanical DOF. It is purely a geometric orientation variable used to enforce correct bond/angle directionality. Because there is no rotational inertia to absorb torque, the linear impulse must be applied as a **central force** (along the atom–atom line) to keep translational angular momentum `L_trans = Σ r × p` conserved.

#### 6.3.1 Key Design Principle
When rotational inertia is zero, the system is effectively a **point-mass forcefield** (like UFF or similar valence forcefields). Off-center port forces would create unbalanced torques with no rotational channel to absorb them, violating translational angular momentum conservation. The fix is to **project the tip-atom constraint error onto the center–center line** and apply the impulse only along that line. The quaternion alignment still enforces correct bond directionality, but the linear dynamics behave like a standard central-force model.

#### 6.3.2 `substep_optimized` — Iterative Newton in Rotation Space
- **Strategy**: Inside each thread, perform several Newton–Raphson substeps in rotation-vector (`ω`) space. The local Hessian `H` is built from port lever arms, and `ω = H⁻¹ τ` is solved. After convergence, the optimized quaternion is used to compute linear recoil.
- **Pipeline**: single kernel `compute_ports_cluster_rigid_substep_optimized`
  - **Consumes**: `pos`, `quat`, `port_local`, `stiffness`
  - **Produces**: `dpos_node`, `drot_node`, `dpos_neigh`
- **Conservation**: Point-mass model. After the centerline projection fix, the linear impulse is applied **only along the center–center line** `n = (xj - xi)/|...|`. This guarantees exact conservation of `P` and translational `L = Σ r × p`.
- **Tradeoff**: Higher per-iteration cost because of repeated local 3×3 solves and quaternion updates, but often reaches geometric alignment faster than a single Jacobi sweep.

#### 6.3.3 `shapematch` — Kabsch / Polar Decomposition
- **Strategy**: Build the covariance matrix `A = Σ w · (xj - xi) ⊗ r_local` from neighbor centers and local port vectors, then extract the optimal rotation matrix `R` via polar decomposition (Newton iterations on `R ← ½ R (3I - RᵀR)`). Convert `R` to a quaternion and use it for the linear recoil pass.
- **Pipeline**: single kernel `compute_ports_cluster_rigid_shapematch`
  - **Consumes**: `pos`, `quat`, `port_local`, `stiffness`
  - **Produces**: `dpos_node`, `drot_node`, `dpos_neigh`
- **Conservation**: Same point-mass central-force guarantee as `substep_optimized`: exact `P` and `L_trans` via centerline projection.
- **Tradeoff**: Direct algebraic orientation solve (no iterative torque integration). Very fast for well-conditioned geometries, but can be sensitive to degenerate or near-planar configurations where the polar decomposition becomes ill-conditioned.

#### 6.3.4 `eigen` — Two-Pass Davenport q-Method with Tips Buffer
- **Strategy**: Explicitly decouple the orientation solve from the linear recoil by using a precomputed **tip-position buffer** (`tips`).
  - *Pass 1* computes the optimal quaternion via the **Davenport q-method** (eigenproblem on a 4×4 symmetric matrix built from weighted tip–neighbor-center vectors).
  - *Pass 2* reads the pre-rotated tip positions from the `tips` buffer and computes **only linear recoil**, with the impulse constrained to the center–center line.
- **Pipeline**:
  1. `compute_optimal_rotation_eigen`
     - **Consumes**: `pos`, `quat`, `port_local`, `stiffness`
     - **Produces**: `drot_node`, `quat_opt` (stored temporarily in `dquat_mom`)
  2. `compute_ports_cluster_rigid_eigen_tips`
     - **Consumes**: `pos`, `tips`, `stiffness`
     - **Produces**: `dpos_node`, `dpos_neigh`
  3. `apply_corrections_rigid_ports`
     - Updates `pos`, `quat`, and writes the new `tips` for the next iteration
- **Conservation**: Exact `P` and `L_trans`. The recoil uses a **bidirectional tip→atom** constraint (tip of A toward center of B, and vice-versa), but the impulse direction is the central line between the two atom centers, eliminating net torque on the point-mass system.
- **Tradeoff**: Requires an extra `tips` buffer (`nnode_tot × 4 × 16` bytes), but this is tiny compared to the atom buffers. The two-pass design is the cleanest for exact momentum conservation because it avoids any “double rotation” of neighbor tips in the same kernel and makes the symmetry explicit.

---

### 6.4 Comparison Summary

| Class | Kernel | Rot. inertia | Orientation solver | Recoil direction | Conserved quantities | Relative cost |
|-------|--------|--------------|-------------------|------------------|----------------------|---------------|
| **Massfull** | `orig` / `current` | Yes (physical) | XPBD joint solve | Off-center (tip→atom) | `P`, `L_total = L_trans + L_spin` | Baseline |
| **Massless** | `substep_optimized` | No (geometric) | Newton–Raphson in `ω` | Central (projected) | `P`, `L_trans` | Moderate |
| **Massless** | `shapematch` | No (geometric) | Kabsch / Polar Decomp. | Central (projected) | `P`, `L_trans` | Moderate |
| **Massless** | `eigen` | No (geometric) | Davenport q-method | Central (projected) | `P`, `L_trans` | 2 kernels |

