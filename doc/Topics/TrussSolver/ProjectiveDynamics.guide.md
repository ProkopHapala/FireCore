# Projective Dynamics: A Fast and Stable Physics Solver

## Introduction

Projective Dynamics (PD) is a fast and robust numerical solver for physics-based animation, particularly useful for simulating deformable objects like cloth, soft bodies, and mass-spring systems. It was introduced by Bouaziz et al. in their 2014 paper "Projective Dynamics: Fusing Constraint Projections for Fast Simulation" [1].

The key advantages of PD include:
- Fast convergence compared to traditional explicit integration methods
- Unconditional stability (allows larger timesteps)
- Simple implementation
- Easy parallelization
- Natural handling of constraints

## Mathematical Foundation

### System Description

Consider a system of particles with positions $\mathbf{x}_i$ and masses $m_i$. These particles are connected by springs with rest lengths $l_0$ and stiffness coefficients $k$.

The system's energy consists of two main terms:
1. **Inertial energy**: keeps particles close to their predicted positions
2. **Constraint energy**: maintains spring rest lengths

### Energy Formulation

The total energy to be minimized at each timestep is:

$$E(\mathbf{x}) = \sum_i \frac{m_i}{2h^2} \|\mathbf{x}_i - \tilde{\mathbf{x}}_i\|^2 + \sum_j \frac{k_j}{2} \|\mathbf{d}_j - l_{0,j}\frac{\mathbf{d}_j}{\|\mathbf{d}_j\|}\|^2$$

where:
- $h$ is the timestep
- $\tilde{\mathbf{x}}_i$ is the predicted position
- $\mathbf{d}_j = \mathbf{x}_{j2} - \mathbf{x}_{j1}$ is the current spring vector
- $l_{0,j}$ is the rest length of spring $j$

## Implementation Details

### System Matrix Construction

The system matrix $A$ combines mass and spring contributions: $A = M/h^2 + L$ where:
- $M$ is the diagonal mass matrix
- $L$ is the Laplacian matrix of the spring network
- $h$ is the timestep

For each particle $i$:
1. Add mass term: $A_{ii} += m_i/h^2$
2. For each connected spring $j$:
   - Add stiffness to diagonal: $A_{ii} += k_j$
   - Add negative stiffness to off-diagonal: $A_{ij} = A_{ji} = -k_j$

```python
def make_pd_matrix(neighbs, bonds, masses, dt, ks):
    """Create the system matrix for projective dynamics"""
    np_total = len(masses)
    A = np.zeros((np_total, np_total))
    idt2 = 1.0 / (dt * dt)
    
    for i in range(np_total):
        # Mass term
        Aii = masses[i] * idt2
        # Spring terms
        for ib in neighbs[i]:
            k = ks[ib]
            i_, j_ = bonds[ib]
            Aii += k
            if j_ > i:
                A[i, j_] = -k
                A[j_, i] = -k
            elif i_ > i:
                A[i, i_] = -k
                A[i_, i] = -k
        A[i, i] = Aii
    
    Mt = masses * idt2
    return A, Mt
```

### Right-Hand Side Construction

The right-hand side vector $\mathbf{b}$ combines inertial and spring forces:

$$\mathbf{b}_i = \frac{m_i}{h^2}\tilde{\mathbf{x}}_i + \sum_j k_j l_{0,j} \frac{\mathbf{d}_j}{\|\mathbf{d}_j\|}$$

where:
1. First term: inertial force keeping particle near predicted position
2. Second term: spring forces trying to maintain rest lengths

```python
def make_pd_rhs(neighbs, bonds, masses, dt, ks, points, l0s, pnew):
    """Build the right-hand side of the system"""
    np_total = len(masses)
    b = np.zeros((np_total, 3))
    idt2 = 1.0 / (dt * dt)
    
    for i in range(np_total):
        # Mass term (inertial prediction)
        bi = pnew[i] * (masses[i] * idt2)
        
        # Spring terms
        for ib in neighbs[i]:
            k = ks[ib]
            i_, j_ = bonds[ib]
            j = j_ if i_ == i else i_
            
            # Using predicted positions for better propagation
            d = pnew[i] - pnew[j]
            d_norm = np.linalg.norm(d)
            if d_norm > 1e-10:  # Avoid division by zero
                d *= k * l0s[ib] / d_norm
                bi += d
        
        b[i] = bi
    
    return b
```

### Main Solver Loop

The solver implements a semi-implicit integration scheme:

1. **Prediction Step**:
   $$\tilde{\mathbf{x}} = \mathbf{x} + h\mathbf{v} + h^2\mathbf{g}$$
   where $\mathbf{g}$ is gravity

2. **Solve Linear System**:
   $$A\mathbf{x}_{t+1} = \mathbf{b}$$

3. **Update Velocity**:
   $$\mathbf{v}_{t+1} = (\mathbf{x}_{t+1} - \mathbf{x}_t)/h$$

```python
def solve_pd(points, velocity, bonds, masses, ks, dt=0.1, n_iter=100, 
             gravity=np.array([0, -9.81, 0]), fixed_points=None, 
             call_back=None, damping=0.01):
    """Solve the system using projective dynamics"""
    # Initialize
    n_points = len(points)
    neighbs = build_neighbor_list(bonds, n_points)
    A, Mt = make_pd_matrix(neighbs, bonds, masses, dt, ks)
    
    pos = points.copy()
    pos_cor = points.copy()
    pos_pred = points.copy()
    l0s = np.array([np.linalg.norm(points[j] - points[i]) 
                    for i, j in bonds])
    
    # Main simulation loop
    for itr in range(n_iter):
        # Apply damping and gravity
        velocity *= (1-damping)
        velocity[:,:] += gravity[None,:] * dt
        
        # Predict step
        pos_pred[:] = pos + velocity * dt
        
        # Apply constraints
        if fixed_points is not None:
            pos_pred[fixed_points] = points[fixed_points]
        
        # Build and solve system
        b = make_pd_rhs(neighbs, bonds, masses, dt, ks, pos, l0s, pos_pred)
        for i in range(3):
            pos_cor[:, i] = np.linalg.solve(A, b[:, i])
        
        # Update velocity and position
        velocity = (pos_cor - pos) / dt
        pos[:] = pos_cor[:]
        
        if call_back is not None:
            call_back(pos)
    
    return pos, velocity
```

### Efficient System Solver

For best performance, we use LDLT decomposition (root-free Cholesky) because:

1. The system matrix $A$ is symmetric positive definite
2. The decomposition $A = LDL^T$ where:
   - $L$ is lower triangular with ones on diagonal
   - $D$ is diagonal
   - Avoids square roots for better numerical stability

The system $A\mathbf{x} = \mathbf{b}$ is solved in three steps:
1. Forward substitution: $L\mathbf{z} = \mathbf{b}$
2. Diagonal scaling: $\mathbf{y} = D^{-1}\mathbf{z}$
3. Backward substitution: $L^T\mathbf{x} = \mathbf{y}$

```python
# Precompute decomposition
L, D, _ = scipy.linalg.ldl(A)

def solve_ldlt(L, D, b):
    # Solve Lz = b
    z = scipy.linalg.solve_triangular(L, b, lower=True)
    # Solve Dy = z
    y = z / np.diag(D)
    # Solve L^Tx = y
    x = scipy.linalg.solve_triangular(L.T, y, lower=False)
    return x
```

## References

[1] Bouaziz, S., Martin, S., Liu, T., Kavan, L., & Pauly, M. (2014). Projective dynamics: Fusing constraint projections for fast simulation. ACM Transactions on Graphics, 33(4), 1-11.

[2] Liu, T., Bouaziz, S., & Kavan, L. (2017). Quasi-newton methods for real-time simulation of hyperelastic materials. ACM Transactions on Graphics, 36(4), 1-16.

[3] Wang, H. (2015). A Chebyshev semi-iterative approach for accelerating projective and position-based dynamics. ACM Transactions on Graphics, 34(6), 1-9.
