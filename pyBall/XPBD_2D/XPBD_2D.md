# XPBD_2D Design Documentation

## Overview

XPBD_2D is a 2D Position-Based Dynamics (PBD) simulator for rigid bodies with port-based constraints. Unlike its 3D counterpart that uses quaternions for rotation, XPBD_2D uses **complex numbers** to represent 2D rotations, simplifying the math while maintaining physical correctness.

## Key Design Decisions

### 1. Complex Numbers for 2D Rotation

In 2D, a rotation by angle θ can be represented as a complex number:
```
z = cos(θ) + i·sin(θ) = (cos(θ), sin(θ))
```

**Why complex numbers?**
- Simpler than quaternions (2 components vs 4)
- Direct representation of 2D rotation
- Efficient rotation via complex multiplication
- No gimbal lock issues (unlike Euler angles)

**Rotation formula:**
To rotate a vector `v = (vx, vy)` by complex `z = (zx, zy)`:
```
z * v = (zx·vx - zy·vy, zy·vx + zx·vy)
```

**Complex multiplication:**
```
z1 * z2 = (z1x·z2x - z1y·z2y, z1y·z2x + z1x·z2y)
```

### 2. Data Structures

| Variable | 3D (XPDB_new) | 2D (XPBD_2D) |
|----------|---------------|--------------|
| Position | float3 | float2 |
| Rotation | float4 (quaternion) | float2 (complex) |
| Velocity | float3 | float2 |
| Angular velocity | float3 | float (scalar) |
| Force | float3 | float2 |
| Torque | float3 | float (scalar) |

### 3. Port-Based Constraints

Each rigid node has up to 4 ports representing connection points:
- **Local port offset**: Vector from node center to port in local coordinates
- **Global port position**: `p_node + rotate(rotation, port_local)`
- **Constraint**: Connected ports should coincide

The constraint violation is:
```
C = |p_i + rotate(z_i, r_i) - p_j|
```

Where:
- `p_i`, `p_j`: positions of nodes i and j
- `z_i`: complex rotation of node i
- `r_i`: local port offset of node i

### 4. Solver Methods

#### 4.1 Explicit Force Integration

**Process:**
1. Clear force buffers
2. Gather spring forces from port constraints
3. Apply recoil forces to neighbors
4. Integrate linear and angular momentum

**Force calculation:**
```
f = K * (p_j - (p_i + rotate(z_i, r_i)))
```

**Torque calculation (2D cross product):**
```
τ = r.x * f.y - r.y * f.x
```

**Integration:**
```
v ← v·damp + (f/m)·dt
p ← p + v·dt
ω ← ω·damp + (τ/I)·dt
θ ← θ + ω·dt
z ← (cos(θ), sin(θ))
```

#### 4.2 PBD (Position-Based Dynamics)

Iterative constraint projection without velocity:

1. **Compute corrections**: For each constraint, calculate position and rotation corrections
2. **Apply corrections**: Update positions/rotations with relaxation factor

**Position correction:**
```
Δp_i = λ * w_i * n
Δp_j = -λ * w_j * n
```

**Rotation correction (2D):**
```
Δθ = (r × Δp) / |r|²
```

Where `×` is the 2D cross product (scalar result).

#### 4.3 XPBD (Extended PBD)

Adds compliance (inverse stiffness) and Lagrange multipliers:

```
α = 1 / (K · dt²)     # compliance
λ_new = λ_old + Δλ    # accumulate multiplier
Δλ = (-C - α·λ_old) / (w_total + α)
```

Benefits:
- Converges to exact constraint satisfaction
- Better energy conservation
- Tunable stiffness via compliance

### 5. Momentum Conservation

The solver uses a **gather-scatter** pattern to ensure momentum conservation:

1. **Compute phase**: Each node computes corrections for its ports
2. **Gather phase**: Corrections are accumulated symmetrically

This avoids race conditions and ensures:
- Linear momentum: `Σ m_i · v_i` conserved
- Angular momentum: `Σ (r_i × p_i + I_i · ω_i)` conserved

### 6. File Structure

```
XPBD_2D.cl     - OpenCL kernels
XPBD_2D.py     - Python wrapper class
XPBD_2D.md     - This documentation
test_XPBD_2D.py - Standalone test with animation
```

## Kernel Reference

### clear_2d_forces
Clears force buffer for all atoms.

### clear_2d_node_buffers
Clears port force and lever arm buffers for nodes.

### gather_port_forces_2d
Computes spring forces from port constraints.
- Input: positions, rotations, neighbors, stiffness
- Output: forces, recoil forces, lever arms

### integrate_2d_explicit
Explicit Euler integration with damping.
- Updates: positions, velocities, rotations, angular velocities
- Gathers recoil forces via bkSlots

### compute_corrections_2d
PBD: computes position and rotation corrections.
- Calculates constraint violations
- Computes weighted corrections
- Stores neighbor corrections for gather phase

### apply_corrections_2d
Applies PBD corrections with relaxation.
- Gathers corrections from connected nodes
- Updates positions and rotations

### compute_xpbd_corrections_2d
XPBD: computes corrections with compliance.
- Accumulates Lagrange multipliers
- Uses compliance α = 1/(K·dt²)

### gather_and_apply_xpbd_2d
Gathers and applies XPBD corrections.

### reset_lambda_2d
Resets Lagrange multiplier buffer.

## Usage Example

```python
from XPBD_2D import XPBD_2D, build_neighs_bk_from_bonds_2d, make_bk_slots_2d
import numpy as np

# Create simulator
sim = XPBD_2D(num_atoms=10)

# Define bonds (topology)
bonds = [(0, 1), (1, 2), (2, 3)]
neighs, bks = build_neighs_bk_from_bonds_2d(10, bonds)
bkSlots = make_bk_slots_2d(neighs, nnode=4, natoms=10)

# Upload topology
sim.upload_topology(neighs, bkSlots, stiffness)

# Define port geometry
port_local = np.zeros((4, 4, 2), dtype=np.float32)
port_local[0, 0] = [1.0, 0.0]   # Node 0, port 0: right
port_local[1, 0] = [-1.0, 0.0]  # Node 1, port 0: left
port_n = np.array([1, 1, 0, 0], dtype=np.uint8)

sim.upload_ports(port_local, port_n, nnode=4)

# Set initial state
pos = np.random.randn(10, 2).astype(np.float32)
sim.upload_state(pos)

# Run simulation
for _ in range(100):
    sim.step_xpbd(nnode=4, dt=0.01, iterations=10)
    pos, rot, vel, omega = sim.download_state()
```

## Differences from 3D (XPDB_new)

| Aspect | 3D | 2D |
|--------|-----|-----|
| Rotation | Quaternion (float4) | Complex (float2) |
| Vector ops | float3 | float2 |
| Angular vel | float3 (vector) | float (scalar) |
| Torque | float3 | float (scalar) |
| Inertia | Tensor (3x3) | Scalar |
| Cross product | 3D (returns vector) | 2D (returns scalar) |
| Memory | Higher | Lower |
| Speed | Slower | Faster |

## Future Extensions

Potential additions:
- Collision detection with 2D circles
- Polygon-based rigid bodies
- External forces (gravity, wind)
- Constraint types (distance, angle, pin)
- GUI for interactive molecule building
