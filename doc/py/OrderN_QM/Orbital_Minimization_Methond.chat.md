https://drive.google.com/file/d/10B4ojUt-DkFaRUYlPv91nnVFsZuO9N1h/view?usp=sharing, https://aistudio.google.com/app/prompts?state=%7B%22ids%22:%5B%221sRhqGDgxUKTm7Al8UO_WRhVAMuwD99Ux%22%5D,%22action%22:%22open%22,%22userId%22:%22100958146796876347936%22,%22resourceKeys%22:%7B%7D%7D&usp=sharing

## USER 1

OK, so Greens function solver works quite well, much better than FOE.

Now I would like to try my solver, based on localized orbitals. Basically I want to minimize energy of localized molecular orbitals (like Fragment molecular orbitals, Wannier functions, or natural bond orbitals).

The main computing cost in many methods (like Car-Parrignelo and others) is ensure Pauli exlusion principle which is equvalent to othrogonality between orbitals. We can do it by applying some lagrange multipliers and other orthogonalization steps.

How I want to do it is by "projective dynamics" simolar to what is used in game engines, where I do movement of orbitals down the hill (down the energy) and then solve system of equations (by fast method such as Cholsky on CPU and Jacobi or Gauss seidel on GPU) to ensure orthogonality constrains.

By construction my Molecular orbitals  have only finite support - they are centered around some site and cannot go more to certain distance due to additional confinement potential which is flat from the start and then grow steaply. Or some other confinement method, but I need to ensure that they are really zero beynd some cutof, so I can rely that overlap and hopping between orbitals futher appart are strictly zero. This will significantly reduce the computational cost. from O(n^3) to O(nk^2) where n is total number of basis functions and k is number of basis functions within cutof R_cut.

Are there any literature references for similar methods? Are there different variants of similar method? Do you have some ideas, or recomandations how to do it well?

---

This is a very solid approach. What you are describing aligns closely with the **Orbital Minimization Method (OMM)** and methods used in **Linear Scaling DFT** codes like **ONETEP** or **CONQUEST**, though your specific angle of "Projective Dynamics" (splitting the gradient step and the constraint step) is a very intuitive way to frame it, often called **Projected Gradient Descent (PGD)** on the Stiefel Manifold in mathematics.

Here is the breakdown of the literature, the theory, and a concrete algorithm recommendation.

### 1. Literature & Similar Methods

The family of methods you are looking for deals with **Non-orthogonal Generalized Wannier Functions (NGWFs)** or **Support Functions**.

*   **Mauri-Galli-Car (1993):** The foundational paper. They introduced a functional that, when minimized, naturally drives orbitals towards orthogonality without requiring an explicit $O(N^3)$ Gram-Schmidt step at every iteration.
*   **ONETEP (Skylaris et al.):** They use strictly localized orbitals (psinc functions inside a sphere). They optimize the coefficients of these orbitals. Their method is mathematically very similar to what you want: they minimize energy with respect to the orbital coefficients subject to a localization constraint.
*   **SIESTA (Ordejon et al.):** They use numerical atomic orbitals with a strict cutoff radius determined by an energy shift parameter (confinement potential), exactly as you proposed.
*   **Manifold Optimization:** In applied math, this is optimizing $F(X)$ subject to $X^T X = I$. The "Projective Dynamics" approach corresponds to a **Retraction** step on the Stiefel Manifold.

### 2. The Algorithm: Projected Gradient Descent with Localization

You want to solve:
$$ \min_{\{C\}} E[\{C\}] \quad \text{s.t.} \quad C^T S C = I \quad \text{and} \quad C_{i\alpha} = 0 \text{ if } r_\alpha > R_{cut} $$

Here is the "Projective Dynamics" recipe adapted for Quantum Chemistry:

#### Step A: Unconstrained Gradient Descent (The "Hill" Descent)
Calculate the gradient of the energy with respect to the orbital coefficients $C$.
$$ C^* = C_{old} - \eta \nabla E $$
*   **Localization:** Immediately apply your mask. If a coefficient corresponds to a basis function outside the cutoff radius of the orbital, force it to zero.
    $$ C^*_{i\alpha} = C^*_{i\alpha} \times \text{Mask}_{i\alpha} $$
*   *Note:* The gradient $\nabla E$ is roughly $H C$. Since $H$ is sparse and $C$ is sparse (due to localization), this matrix multiplication is $O(N)$.

#### Step B: The Constraint Projection (Orthogonalization)
Now $C^*$ is lower in energy but no longer orthogonal ($C^{*T} S C^* \neq I$). We must project it back to the "Orthogonality Manifold".

The "closest" valid set of orbitals (in a least-squares sense) is given by **Löwdin Symmetric Orthogonalization**:
$$ C_{new} = C^* (C^{*T} S C^*)^{-1/2} $$

This looks like it requires an inverse square root (expensive), **BUT** since your orbitals are localized and $C^*$ is already close to orthogonal (from the previous step), the overlap matrix $O = C^{*T} S C^*$ is **sparse** and **close to identity**.

#### Step C: Fast Solver for the Constraint (GPU Friendly)
You asked how to solve this efficiently without full diagonalization.

**1. Newton-Schulz Iteration (The GPU King)**
You can compute $X^{-1/2}$ purely using matrix multiplications (no Cholesky, no inverse). If $O \approx I$, this converges quadratically in very few steps.
Algorithm to find $Z = O^{-1/2}$:
1.  Initialize $Z_0 = I$ (or better, scaled $I$).
2.  Iterate: $Z_{k+1} = \frac{1}{2} Z_k (3I - O Z_k^2)$
3.  Apply: $C_{new} = C^* Z_{final}$

**2. Cholesky Projection (CPU/Sparse)**
If you treat the orthogonality as a linear system of Lagrange multipliers (Projective Dynamics style):
$$ (C^* + \delta C)^T S (C^* + \delta C) = I $$
Linearizing this leads to a system that looks like $A \lambda = B$.
You can solve this using Sparse Cholesky if the sparsity pattern of the orbitals is fixed.

### 3. Recommendation for your Implementation

Since you want **Linear Scaling ($O(N)$)** and **Sparsity**:

1.  **Strict Confinement:** Define a list of "neighboring basis functions" for each orbital center. The matrix $C$ becomes a "block-sparse" or "ragged" array. You only store and compute interactions within $R_{cut}$.
2.  **The Overlap Approximation:**
    Do not enforce perfect orthogonality ($10^{-14}$) at every single time step. In "Projective Dynamics" (and Car-Parrinello), it is sufficient that the orbitals stay *close* to the constraint surface.
    Use **iterative refinement**. Apply just 1 or 2 iterations of Newton-Schulz per physics step. The constraint error will decay alongside the energy minimization.
3.  **Visualizing the "Projection":**
    Think of the orthogonality constraint as a "spring" that pulls the orbitals apart if they overlap too much.

### 4. Python Prototype (Constraint Projection Solver)

Here is a script that demonstrates **Energy Minimization with Strict Localization and Iterative Orthogonalization**.

```python
import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt

def generate_1d_system(N=20, localized_dim=5):
    """
    Generates a sparse Hamiltonian H and Overlap S for a 1D chain.
    """
    # Simple tight binding H
    diag = np.ones(N) * 0.0
    off = np.ones(N-1) * -1.0
    H = np.diag(diag) + np.diag(off, k=1) + np.diag(off, k=-1)
    
    # Non-orthogonal S (tridiagonal)
    s_off = np.ones(N-1) * 0.2
    S = np.eye(N) + np.diag(s_off, k=1) + np.diag(s_off, k=-1)
    
    return H, S

def make_localization_mask(N, n_orb, radius=2):
    """
    Creates a boolean mask where Orb_i is allowed to be non-zero 
    only on sites [center-radius, center+radius].
    """
    mask = np.zeros((N, n_orb), dtype=bool)
    # Distribute orbital centers evenly
    centers = np.linspace(0, N-1, n_orb, dtype=int)
    
    for i, c in enumerate(centers):
        start = max(0, c - radius)
        end = min(N, c + radius + 1)
        mask[start:end, i] = True
    return mask

def solve_newton_schulz_inverse_sqrt(M, steps=5):
    """
    Computes M^(-1/2) using Newton-Schulz iteration.
    Requires M to be close to Identity (scaled).
    GPU friendly: only uses MatMul.
    """
    # 1. Scale M to ensure convergence (norm < 1)
    norm = np.linalg.norm(M, ord=np.inf)
    M_scaled = M / norm
    
    X = np.eye(M.shape[0])
    Y = M_scaled.copy()
    
    # Iteration: X_{k+1} = 0.5 * X_k * (3I - Y_k X_k)
    # Actually, efficient form for Inverse Sqrt Z = M^-0.5:
    # Z_0 = I
    # Z_{k+1} = 0.5 * Z_k * (3I - M_scaled * Z_k * Z_k)
    
    Z = np.eye(M.shape[0])
    I = np.eye(M.shape[0])
    
    for k in range(steps):
        # T = Z @ Z
        # diff = 3I - M_scaled @ T
        # Z = 0.5 * Z @ diff
        Z = 0.5 * Z @ (3.0 * I - M_scaled @ (Z @ Z))
        
    return Z / np.sqrt(norm)

def optimize_localized_orbitals(H, S, N, n_occ, steps=100, learning_rate=0.1):
    # 1. Define Localization Mask
    # Each orbital is strictly zero outside its radius
    R_cut = 3
    mask = make_localization_mask(N, n_occ, radius=R_cut)
    
    # 2. Random initialization within mask
    C = (np.random.rand(N, n_occ) - 0.5) * mask
    
    # Initial orthogonalization
    # C <- C * (C.T S C)^-1/2
    O = C.T @ S @ C
    inv_sqrt_O = solve_newton_schulz_inverse_sqrt(O, steps=15) # Rigid initial projection
    C = C @ inv_sqrt_O
    
    energies = []
    constraints_err = []
    
    print(f"#INFO Starting Optimization: N={N} Occ={n_occ} R_cut={R_cut}")
    
    for k in range(steps):
        # --- A. Unconstrained Gradient Descent ---
        # E = Tr(C.T H C)
        # dE/dC = 2 H C
        grad = 2.0 * H @ C
        
        # Update
        C_next = C - learning_rate * grad
        
        # --- B. Enforce Strict Localization (Masking) ---
        C_next = C_next * mask
        
        # --- C. Projective Dynamics (Orthogonalization) ---
        # We need to satisfy C.T S C = I
        # Calculate Overlap matrix of the updated orbitals
        M = C_next.T @ S @ C_next
        
        # Check error before projection
        err = np.max(np.abs(M - np.eye(n_occ)))
        constraints_err.append(err)
        
        # Compute correction matrix Z = M^-1/2
        # Since M is close to I, Newton-Schulz converges fast (2-3 steps is often enough)
        Z = solve_newton_schulz_inverse_sqrt(M, steps=3)
        
        # Apply Correction
        C_next = C_next @ Z
        
        # Update
        C = C_next
        
        # Calc Energy
        E = np.trace(C.T @ H @ C)
        energies.append(E)
        
        if k % 10 == 0:
            print(f"Iter {k:3d} E={E:.5f} MaxConstraintErr={err:.1e}")

    return C, energies, constraints_err

# --- Run ---
N = 50
n_occ = 10
H, S = generate_1d_system(N=N)
C_final, Es, errs = optimize_localized_orbitals(H, S, N, n_occ)

# Plot
fig, ax = plt.subplots(1, 2, figsize=(10, 4))
ax[0].plot(Es, 'b.-')
ax[0].set_title("Energy Minimization")
ax[0].set_xlabel("Step")
ax[0].set_ylabel("Energy [Ha]")

ax[1].plot(errs, 'r.-')
ax[1].set_yscale('log')
ax[1].set_title("Orthogonality Violation (Before Projection)")
ax[1].set_xlabel("Step")

plt.tight_layout()
plt.show()

# Visualize Matrix Sparsity
plt.figure(figsize=(5,5))
plt.imshow(np.abs(C_final) > 1e-4, aspect='auto', interpolation='nearest', cmap='Greys')
plt.title("Localized Orbital Coefficients (Sparsity)")
plt.xlabel("Orbital Index")
plt.ylabel("Basis Site")
plt.show()
```

### Key Takeaway for Linear Scaling

In the code above, `C` is stored as a dense $(N \times N_{occ})$ matrix for simplicity. To achieve true $O(N)$:

1.  **Storage:** `C` must be a sparse matrix or a list of small dense vectors (one per orbital).
2.  **Overlap $M$:** The matrix $M = C^T S C$ will also be sparse because orbitals far apart do not overlap (due to the mask).
3.  **Newton-Schulz:** The matrix multiplications $Z \times Z$ in the Newton-Schulz step become **Sparse-Sparse Matrix Multiplications (SpGEMM)**. This is the holy grail of GPU acceleration (libraries like cuSPARSE or optimized block-sparse routines).

This fits your goal perfectly: you project the dynamics back to the manifold using iterative sparse algebra.


---

## USER 2

From all the methods you mentioned and proposed, analyse each of them in term of computationa complezity O(n^k), considering the typical shape of the matrixes S and H (they are rather sparse, but when we try to multiply them or do LU / LDL factroization they quckly fill in). This is why I prefer linear solvers using Jacobi or Gauss seidel which was should to be quite efficient on GPU and correspond to very sparse matrix multiplication. Your secription of these mathods in matematical terms was abstrac and it was not obvious how costly it is but when I see Lowdin transform or Z=O^(-1/2) I smell it will be costly matrix inversion or diagonalization we I definitely want to avoid  


---

This is the most critical question for high-performance implementation. You are absolutely right: **Sparsity is fragile.**

In sparse linear algebra, **fill-in** is the enemy.
*   **Sparse $\times$ Dense vector** $\rightarrow$ **Dense** (Cheap, $O(N)$)
*   **Sparse $\times$ Sparse matrix** $\rightarrow$ **Less Sparse** (Expensive, fill-in grows)
*   **Sparse Inversion / Factorization** $\rightarrow$ **Dense** (Catastrophic, $O(N^2)$ or $O(N^3)$)

Here is the complexity analysis of the methods we discussed, specifically evaluating how they handle sparsity and GPU efficiency.

---

### 1. Green's Function (Contour Integration)
**Equation:** $P = \oint (zS - H)^{-1} dz$
**Core Operation:** Solving $Ax = b$ where $A = (zS - H)$.

*   **Variant A: Direct Solver (LU / Cholesky)**
    *   **1D Chains:** $O(N)$. Tridiagonal matrices do *not* suffer from fill-in during factorization.
    *   **2D/3D Grids:** **Fails linear scaling.** Factorization causes massive fill-in (the "skyline" of the matrix fills up). Complexity is $\approx O(N^{1.5})$ to $O(N^2)$.
    *   *Verdict:* Do not use for 3D on GPU.

*   **Variant B: Iterative Solver (Jacobi / CG / COCG)**
    *   **Operation:** Only requires matrix-vector multiplication ($A \cdot v$).
    *   **Sparsity:** strictly preserved.
    *   **Complexity:** $N_{poles} \times N_{iter} \times O(N)$.
    *   **GPU Suitability:** Excellent. SpMV (Sparse Matrix-Vector) is the benchmark for sparse GPU performance.
    *   *Bottleneck:* If the condition number is bad (small gap), $N_{iter}$ explodes. Preconditioning is required, which is hard to do linearly on GPU.

---

### 2. Orbital Minimization (Projective Dynamics)
**Equation:** Minimize $E(C)$ s.t. $C^T S C = I$.
**Core Operation:** Orthogonalization.

This is where your fear of $S^{-1/2}$ lies. Let's analyze the cost of the **Newton-Schulz** iteration ($Z_{k+1} = \frac{1}{2}Z_k(3I - M Z_k^2)$) vs **Lagrange Multipliers**.

#### The Problem: "Sparse $\times$ Sparse" (SpGEMM)
Even if your orbitals $C$ are localized (sparse), the overlap matrix $M = C^T S C$ is sparse, and $Z = M^{-1/2}$ is mathematically sparse (exponential decay), **multiplying them creates fill-in.**
*   If you multiply two sparse matrices with bandwidth $K$, the result has bandwidth $2K$.
*   After 5 Newton-Schulz steps, the bandwidth is $32K$. Your matrix effectively becomes dense locally.

#### Solution: The "Iterative Lagrange" Method (Pure SpMV)
Instead of computing $Z = M^{-1/2}$ (which requires SpGEMM), you can enforce orthogonality by solving for **Lagrange Multipliers** using only linear solvers (Jacobi/CG).

**The Math:**
We want $C_{new} = C_{old} + \Delta C$ such that constraints are satisfied.
Linearizing the constraint $(C+\Delta C)^T S (C+\Delta C) = I$ gives a linear system for the forces.

**Algorithm for $O(N)$ on GPU using only Jacobi:**

1.  **Gradient Step:**
    $C^* = C - \eta H C$
    *   *Cost:* Sparse Matrix $\times$ Sparse Matrix? NO.
    *   Since $C$ is a collection of localized vectors, you perform $N_{occ}$ parallel SpMV operations.
    *   **Complexity:** $O(N)$. **Sparsity Preserved.**

2.  **Overlap Check:**
    Compute $O_{ij} = \langle \psi_i | S | \psi_j \rangle$ only for neighboring orbitals.
    *   *Cost:* Local dot products.
    *   **Complexity:** $O(N)$.

3.  **Correction (The "Projective" Step):**
    We need to find a correction $\delta C$ to restore orthogonality.
    Instead of $C \leftarrow C (C^T S C)^{-1/2}$, we solve:
    $$ S \cdot \delta C_i = - \sum_{j \in neighbors} \lambda_{ij} \psi_j $$
    This looks like a linear system $S x = b$.
    *   **Solver:** Use **Jacobi** or **Chebyshev Semi-Iterative** on GPU.
    *   **Complexity:** $N_{iter} \times O(N)$.
    *   *Note:* You don't need to solve it perfectly. Just 2-3 Jacobi steps per dynamics step is enough to keep orbitals close to the manifold (like SHAKE/RATTLE in molecular dynamics).

### Comparison Summary

| Method | Core Math Operation | GPU Kernel | Complexity | Fill-in Risk |
| :--- | :--- | :--- | :--- | :--- |
| **FOE** | Polynomials of $H$ | **SpMV** | $O(N)$ | None |
| **Green's (Iterative)** | Linear Solve $(zS-H)x=b$ | **SpMV** | $O(N \times Iter)$ | None |
| **Green's (Direct)** | Factorize $(zS-H)$ | **Decomposition** | $O(N^2)$ (3D) | **High** |
| **OMM (Newton-Schulz)**| Matrix Powers $M^k$ | **SpGEMM** | $O(N \times Width)$ | Moderate |
| **OMM (Lagrange)** | Linear Solve $S \delta C = \text{Err}$ | **SpMV** | $O(N)$ | **None** |

### Recommendation

Since you essentially want to write a "Game Engine for Quantum Mechanics" (Projective Dynamics):

**Go with OMM (Orbital Minimization) using Iterative Lagrange Multipliers.**

1.  **Data Structure:** Store $C$ not as a matrix, but as a list of `{index, value}` pairs or a CSR matrix.
2.  **Strict Cutoff:** Hard-code the neighbor list. Atom $I$ only interacts with Atoms $\{J\}$ within $R_{cut}$.
3.  **The Loop:**
    *   **Step 1 (Physics):** Move coefficients $C_{i} \leftarrow C_{i} - \eta (H C)_i$. (Compute via SpMV).
    *   **Step 2 (Constraint):** Compute Overlap error $E_{ij} = \delta_{ij} - \psi_i^T S \psi_j$ for neighbors.
    *   **Step 3 (Projection):** Update $\psi_i \leftarrow \psi_i + \frac{1}{2} \sum_j E_{ij} \psi_j$.
        *   *Wait, isn't this too simple?* This is actually the first-order approximation of the Newton-Schulz loop ($Z \approx I + \frac{1}{2}(I - O)$).
        *   It involves **no matrix inversion** and **no linear solves**, just vector addition.
        *   If you do this small correction at every time step, the system stays orthogonal naturally.

This is exactly how game engines handle rigid body constraints (Projected Gauss-Seidel). It is $O(N)$, strictly sparse, and trivial to parallelize on GPU.


---


## USER 3

OK I still want to dig-in a bit more. It is not totally clear as it looks.

you proposed method which is  "first-order approximation of the Newton-Schulz loop"
what I do not like this this equation
$\psi_i \leftarrow \psi_i + \frac{1}{2} \sum_j E_{ij} \psi_j$
basicallu it is like 

v += 0.5<u|v>u
v +=c*u
the problem is that vector u has different support than vector v, and if I really do it like that, the orbtials will delocalize (fill-in) as vector u will 
introduce non zero term where the original orbital was zero.

What we should do instead (similarly as projective dynamics) is formulate the orthogonality as minimization problem.

In standard projective dynamics we have some inertial part p_pred scaled by M/dt^2 and then the constriains (often bond lenghts) and we basically minimize the error in constrains by the Jacobi or GS iteration.

I want to do something similar with my localized orbitals by with preserving their localization by construction. That is I variationally search for such updates of the allowed non zero components of each localized MO (keeping the empy zero by construction, not touching them) to minimize the error in orthogonality


For better discussion let me refresh how projective dynamics works

----


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

---

## Gemini

This is an excellent insight. You are absolutely right to draw the parallel with **Projective Dynamics (PD)**.

In classical PD, we define a "Local Step" (projection) and a "Global Step" (linear solve).
You want to apply the same logic:
1.  **Inertial Step:** Move orbitals "downhill" physically ($H$).
2.  **Constraint Step:** Variationally find the minimal update $\Delta C$ (within the sparse mask) that satisfies orthogonality.

Your concern about `v += c*u` causing fill-in is correct. If vector `u` (orbital $j$) extends beyond the cutoff of vector `v` (orbital $i$), adding them naively breaks the sparsity pattern.

Here is the formulation of **Quantum Projective Dynamics** that respects strict localization.

---

### The Mathematical Formulation

We define the total objective function exactly like PD: a weighted sum of an "Inertial Potential" and a "Constraint Potential".

$$ J(C) = \underbrace{\frac{1}{2\eta} \| C - \tilde{C} \|_F^2}_{\text{Inertial Term}} + \underbrace{\frac{k}{2} \sum_{i,j} \left( \psi_i^T S \psi_j - \delta_{ij} \right)^2}_{\text{Constraint Energy}} $$

Where:
*   $\tilde{C} = C_{old} - \eta H C_{old}$ (The predicted "Physical" state).
*   The norm $\| \cdot \|_F^2$ is restricted to the **allowed sparsity mask**.

We want to find $C$ that minimizes $J(C)$ subject to the mask constraints.
Taking the gradient with respect to the active coefficients of orbital $i$ ($\psi_i$):

$$ \nabla_{\psi_i} J = \frac{1}{\eta}(\psi_i - \tilde{\psi}_i) + k \sum_{j} 2 (\psi_i^T S \psi_j - \delta_{ij}) S \psi_j = 0 $$

To solve this efficiently (without a global Hessian), we use the **Block Coordinate Descent (Gauss-Seidel)** approach, similar to how games solve rigid body chains.

### The Algorithm: Iterative Constraint Minimization

Instead of a global linear solve (which would be $O(N^3)$), we solve for the update of orbital $i$ while holding neighbors fixed, then iterate. This is $O(N)$ and GPU-friendly.

#### Step 1: Prediction (Inertial Step)
Compute the unconstrained physical update.
$$ \tilde{\psi}_i = \psi_i^{old} - \eta (H \psi_i^{old}) $$
*   **Crucial:** Immediately apply the mask. If $(H \psi)_k$ is non-zero but $k$ is outside the mask of orbital $i$, drop it.

#### Step 2: Constraint Solver (The "Projective" Loop)
We perform $K$ iterations (e.g., 5-10) of minimizing the orthogonality error.
For each orbital $i$ (can be parallelized via graph coloring):

1.  **Calculate Errors:** Compute overlap with all neighbors $j$.
    $$ E_{ij} = \delta_{ij} - \psi_i^T S \psi_j $$
2.  **Compute Constraint Force:**
    The "force" pushing $\psi_i$ to restore orthogonality is a linear combination of neighbor orbitals:
    $$ F_i = \sum_{j \in \text{neigh}} E_{ij} (S \psi_j) $$
    *(Note: $S \psi_j$ is the gradient direction of the overlap constraint).*
3.  **Project Force to Mask:**
    This is the fix for your "fill-in" concern. We essentially truncate the force.
    $$ \tilde{F}_i = \text{Mask}_i( F_i ) $$
    *   This finds the "best possible direction" within the allowed support to fix the orthogonality.
4.  **Update Orbital:**
    $$ \psi_i \leftarrow \psi_i + \alpha \tilde{F}_i $$
    *(Where $\alpha$ is a step size related to the Hessian diagonal, usually $\approx 0.5$)*.

---

### Python Implementation

Here is the implementation of this specific logic. It uses sparse matrices but enforces strict dense-vector-like masks for the orbitals to keep computation linear.

```python
import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt

def make_1d_system_sparse(N=50):
    """Create sparse H and S for 1D chain."""
    # Hamiltonian: -2.0 onsite, -1.0 hopping
    diag = np.ones(N) * -2.0
    off = np.ones(N-1) * -1.0
    H = sp.diags([off, diag, off], [-1, 0, 1], format='csr')
    
    # Overlap: 1.0 onsite, 0.2 neighbor
    s_diag = np.ones(N) * 1.0
    s_off = np.ones(N-1) * 0.2
    S = sp.diags([s_off, s_diag, s_off], [-1, 0, 1], format='csr')
    return H, S

def get_orbital_mask_indices(N, center, radius):
    """Return indices of basis functions allowed for an orbital centered at 'center'."""
    start = max(0, center - radius)
    end = min(N, center + radius + 1)
    return np.arange(start, end)

class LocalizedOrbitalSolver:
    def __init__(self, N, n_orb, radius=3, dt=0.1, orth_stiffness=1.0):
        self.N = N
        self.n_orb = n_orb
        self.radius = radius
        self.dt = dt  # Time step (eta)
        self.k = orth_stiffness # Stiffness of orthogonality constraint
        
        # 1. Define support (Masks)
        # Store as list of index arrays for sparse access
        self.supports = []
        centers = np.linspace(0, N-1, n_orb, dtype=int)
        for c in centers:
            self.supports.append(get_orbital_mask_indices(N, c, radius))
            
        # 2. Initialize coefficients (Random within mask)
        # We store C as a list of small dense vectors, not a big matrix!
        # This guarantees O(N) storage and compute.
        self.C = []
        for i in range(n_orb):
            size = len(self.supports[i])
            vec = np.random.uniform(-1, 1, size)
            vec /= np.linalg.norm(vec) # Normalize initially
            self.C.append(vec)
            
    def get_full_matrix(self):
        """Helper to reconstruct full C matrix for visualization/checking."""
        C_mat = np.zeros((self.N, self.n_orb))
        for i in range(self.n_orb):
            C_mat[self.supports[i], i] = self.C[i]
        return C_mat

    def compute_overlap_matrix(self, S):
        """Compute M = C.T * S * C efficiently using sparsity."""
        # This is strictly restricted to neighbors.
        # In a real code, we would use a neighbor list. 
        # Here we do a naive double loop but checking bounds.
        M = np.zeros((self.n_orb, self.n_orb))
        
        # We only need to compute M_ij if orbitals i and j overlap physically
        # For 1D chain sorted by center, this is band-diagonal.
        
        # Convert S to CSR for fast multiplication
        # But actually, S * C_j is faster if C_j is dense-in-mask.
        
        for i in range(self.n_orb):
            # 1. Expand C_i to full vector (sparse)
            # In optimized C++/CUDA, we wouldn't expand, just gather.
            vec_i = np.zeros(self.N)
            vec_i[self.supports[i]] = self.C[i]
            
            # S_psi_i = S * vec_i
            S_psi_i = S.dot(vec_i)
            
            # Loop neighbors (optimization: only j near i)
            # In 1D, neighbor range is approx 2*radius / separation
            for j in range(self.n_orb):
                # Check bounding box overlap
                if abs(i - j) * (self.N/self.n_orb) > 2*self.radius + 2: 
                    continue
                
                vec_j = np.zeros(self.N)
                vec_j[self.supports[j]] = self.C[j]
                
                # Overlap <psi_i | S | psi_j>
                ovlp = np.dot(vec_i, S_psi_i) if i==j else np.dot(vec_j, S_psi_i)
                M[i, j] = ovlp
        return M

    def step_projective_dynamics(self, H, S, n_constraint_iter=5):
        """
        One step of 'Quantum Projective Dynamics'.
        1. Inertial Update (Gradient Descent on Energy)
        2. Constraint Solver (Projective dynamics on Orthogonality)
        """
        
        # --- PHASE 1: PREDICTION (Inertial Step) ---
        C_pred = []
        for i in range(self.n_orb):
            # Current orbital
            psi = np.zeros(self.N)
            psi[self.supports[i]] = self.C[i]
            
            # Gradient: H * psi
            # Sidenote: In variational optimization, grad is 2*H*psi. 
            # In dynamics, force is -grad.
            force = - H.dot(psi)
            
            # Predict: psi_new = psi + dt * force
            psi_pred = psi + self.dt * force
            
            # MASKING: Project prediction back to allowed support
            # This ensures we never fill-in.
            C_pred.append(psi_pred[self.supports[i]])

        # Update self.C to prediction
        self.C = C_pred

        # --- PHASE 2: CONSTRAINT SOLVER (Local Orthogonality Minimization) ---
        # We want to minimize sum((<i|S|j> - delta_ij)^2).
        # We use Block-Jacobi-ish iteration.
        
        for k in range(n_constraint_iter):
            # 1. Compute current Overlap Matrix M and S-products
            # Optimization: We can update these incrementally, but let's recompute for clarity.
            
            # Precompute S*psi for all orbitals (needed for gradient)
            S_psi_all = []
            for i in range(self.n_orb):
                vec = np.zeros(self.N)
                vec[self.supports[i]] = self.C[i]
                S_psi_all.append(S.dot(vec))
            
            # Compute overlaps M_ij
            M = np.zeros((self.n_orb, self.n_orb))
            for i in range(self.n_orb):
                for j in range(max(0, i-5), min(self.n_orb, i+6)): # Local neighbor loop
                    # Dot product of support vectors
                    # To be fast: use the precomputed S_psi
                    # M_ij = <psi_i | S | psi_j> = psi_i . (S_psi_j)
                    
                    # Extract psi_i full
                    vec_i = np.zeros(self.N)
                    vec_i[self.supports[i]] = self.C[i]
                    
                    val = np.dot(vec_i, S_psi_all[j])
                    M[i, j] = val
                    M[j, i] = val

            # 2. Compute Updates
            # Gradient of Constraint Energy w.r.t psi_i:
            # Grad_i = Sum_j 2 * (M_ij - delta_ij) * (S * psi_j)
            
            updates = []
            for i in range(self.n_orb):
                constraint_grad = np.zeros(self.N)
                
                # Iterate over neighbors involved in constraints
                for j in range(max(0, i-5), min(self.n_orb, i+6)):
                    error = M[i, j] - (1.0 if i==j else 0.0)
                    
                    # Force direction is S * psi_j
                    force_dir = S_psi_all[j]
                    
                    # Accumulate: error * direction
                    constraint_grad += 2.0 * error * force_dir
                
                # 3. PROJECT GRADIENT TO MASK
                # We only allow updates within the support of i
                local_grad = constraint_grad[self.supports[i]]
                
                # Store update (don't apply yet -> Jacobi style)
                # Apply update: psi = psi - learning_rate * grad
                # The learning rate for constraints in PD is usually alpha / (stiffness)
                # Here we use a fixed relaxation parameter
                updates.append(local_grad)
                
            # Apply updates
            alpha = 0.1 # Relaxation factor for constraints
            for i in range(self.n_orb):
                self.C[i] -= alpha * updates[i]

    def get_energy(self, H):
        C_mat = self.get_full_matrix()
        # E = Tr(C.T H C)
        return np.trace(C_mat.T @ H @ C_mat)
        
    def get_orth_error(self, S):
        C_mat = self.get_full_matrix()
        M = C_mat.T @ S @ C_mat
        I = np.eye(self.n_orb)
        return np.linalg.norm(M - I)

# --- Simulation ---
N = 100
n_orb = 20
H, S = make_1d_system_sparse(N)

solver = LocalizedOrbitalSolver(N, n_orb, radius=4, dt=0.05)

energies = []
errors = []

print("Starting Quantum Projective Dynamics...")
for step in range(100):
    solver.step_projective_dynamics(H, S, n_constraint_iter=5)
    
    E = solver.get_energy(H)
    err = solver.get_orth_error(S)
    
    energies.append(E)
    errors.append(err)
    
    if step % 10 == 0:
        print(f"Step {step:3d} Energy={E:.4f} OrthErr={err:.2e}")

# Plot
plt.figure(figsize=(10, 4))
plt.subplot(1,2,1)
plt.plot(energies)
plt.title("Energy Minimization")
plt.xlabel("Step"); plt.ylabel("Energy")

plt.subplot(1,2,2)
plt.plot(errors)
plt.title("Orthogonality Error")
plt.yscale('log')
plt.xlabel("Step")
plt.tight_layout()
plt.show()

# Visual check of localization
plt.figure(figsize=(6,6))
plt.imshow(np.abs(solver.get_full_matrix()) > 1e-4, aspect='auto', interpolation='nearest', cmap='Greys')
plt.title("Strictly Localized Orbitals (Masked)")
plt.xlabel("Orbital"); plt.ylabel("Basis Site")
plt.show()
```

### Why this works
1.  **Strict Sparsity:** The variable `self.C` is never a full matrix. It is a list of small vectors. We manually project any force or gradient onto the mask indices (`self.supports[i]`). This guarantees **zero fill-in**.
2.  **Constraint Minimization:** The inner loop (`n_constraint_iter`) calculates the gradient of the "Orthogonality Energy". By moving `psi_i` in the direction of the masked gradient of the constraint error, we find the best possible orthogonalization *that respects the user-defined localization*.
3.  **Linear Scaling:**
    *   Inertial step: $N_{orb} \times (\text{Support Size})$.
    *   Constraint step: $N_{orb} \times (\text{Neighbors}) \times (\text{Support Size})$.
    *   Total complexity: $O(N)$.

This is exactly the "Game Engine" approach: Move things physically, then iteratively crunch the positions to satisfy constraints, but strictly enforcing geometry (localization) at every sub-step.

---

## USER 4

OK there is article about SIESTA linear scaline, how this compare with our proposed method in terms of complexity / performance?


---

Based on the provided paper (Soler et al., 2002) and your proposed **Quantum Projective Dynamics (QPD)** method, here is a comparative analysis.

### Executive Summary

Both SIESTA and your proposed method fall under the category of **Orbital Minimization Methods (OMM)** using strictly localized support regions to achieve $O(N)$ scaling.

*   **SIESTA (Order-N)** uses the **Kim-Mauri-Galli (KMG) functional**. It treats orthogonality as an energy penalty. It minimizes a single scalar function $E[\{C\}]$ using Conjugate Gradients (CG).
*   **Your Method (QPD)** uses **Projective Dynamics**. It splits the problem into "Physics" (Hamiltonian descent) and "Geometry" (Orthogonality projection).

While mathematically they aim for the same ground state, **your method is likely more efficient on GPUs** and offers better stability control, whereas SIESTA's approach is mathematically more elegant (variational) but computationally heavier per step due to the "density matrix expansion."

---

### 1. The Mathematical Approach

#### SIESTA: The "Soft" Penalty (KMG Functional)
As described in **Section XIII (Order-N Functional)** of the PDF:
SIESTA minimizes a modified energy functional $E_{KMG}$ (Eq. 86).
$$ E_{KMG} = 2 \text{Tr}[(2I - S) (H - \eta S)] $$
This functional is designed such that its unconstrained minimum yields the orthogonal ground state.
*   **Mechanism:** It effectively performs a first-order purification (approximating $S^{-1} \approx 2I - S$) *inside* the energy expression.
*   **Optimization:** It uses Conjugate Gradients (CG). This requires calculating the gradient (Eq. 88), which involves terms like $S_{\mu\nu} c_{\nu j} c_{j\alpha} H_{\alpha\beta}$.

#### Your Method: The "Hard" Projection (Projective Dynamics)
You perform a gradient descent on the standard energy $E = \text{Tr}(C^T H C)$, followed by an explicit geometric correction step to enforce $C^T S C = I$.
*   **Mechanism:** You physically move the orbitals, then "snap" them back to the orthogonality manifold using a local constraint solver (iterative Lagrange/Jacobi).
*   **Optimization:** Semi-implicit integration (Projective Dynamics).

### 2. Computational Complexity & Matrix Operations

Here is where the "GPU-friendliness" difference becomes apparent.

#### SIESTA Complexity
To compute the gradient of the KMG functional (Eq. 88), SIESTA must compute the density matrix $\rho_{\mu\nu}$ (Eq. 87).
$$ \rho_{\mu\nu} = 4 \sum_i c_{\mu i} c_{i \nu} - 2 \sum_{i,j,\alpha,\beta} c_{\mu i} c_{i \alpha} S_{\alpha \beta} c_{\beta j} c_{j \nu} $$
*   **The Problem:** This contains a 4-term contraction. In matrix notation, this looks like $C C^T - C (C^T S C) C^T$.
*   **Fill-in Risk:** Even if $C$ and $S$ are sparse, calculating $\rho$ involves **Sparse Matrix-Matrix Multiplication (SpGEMM)**.
    *   Multiplying two sparse matrices increases the "bandwidth" (fill-in).
    *   SIESTA mitigates this by strictly truncating $\rho$ beyond a cutoff $R_{cut}$ (Section XIII, paragraph 2).
    *   **Cost:** SpGEMM is memory-bandwidth bound and difficult to optimize on GPUs because the output sparsity pattern is not known a priori.

#### Your Method Complexity
Your constraint solver (as defined in the previous Python code) calculates the error:
$$ E_{ij} = \delta_{ij} - \psi_i^T (S \psi_j) $$
and updates:
$$ \psi_i \leftarrow \psi_i + \sum_j E_{ij} (S \psi_j) $$
*   **The Advantage:** This is purely **Sparse Matrix-Vector (SpMV)** and Vector-Vector operations.
*   **No Fill-in:** You never explicitly form the matrix product $C^T S C$ or the density matrix $\rho$. You only compute the scalar overlaps $E_{ij}$ for neighbor pairs.
*   **Cost:** SpMV is the "Hello World" of GPU acceleration. It is strictly $O(N)$ with a very low prefactor and predictable memory access.

### 3. Stability and "Charge Sloshing"

**SIESTA (PDF Page 15):**
> "It is a common experience that the required number of selfconsistency iterations increases with the size of the system... mainly because of the 'charge sloshing' effect."

SIESTA solves this by using Pulay mixing (RMM-DIIS) on the density matrix. This introduces a "history" dependence and significant memory overhead.

**Your Method:**
Projective Dynamics is equivalent to an **Implicit Euler** integration. It is unconditionally stable for the constraints.
*   By treating orthogonality as a "stiff spring" and solving it with Jacobi iterations, you dampen high-frequency oscillations (charge sloshing) naturally.
*   You do not need Pulay mixing for the *orthogonality*; you might only need it for the Hartree potential updates.

### 4. The "Metal" Problem

**SIESTA (PDF Page 15):**
> "If the system is metallic... the KMG functional cannot be used in practice."

This limitation applies to **both** methods.
*   Both methods assume the density matrix decays exponentially (Nearsightedness of Electronic Matter). This is true for insulators/semiconductors but false for metals (algebraic decay).
*   In metals, the "strictly localized" orbitals cannot represent the delocalized electrons, leading to large errors.
*   **SIESTA's solution:** They revert to $O(N^3)$ diagonalization.
*   **Your solution:** You would likely need to implement the Green's Function probing method (discussed previously) for metals, as OMM fails here.

### Conclusion

| Feature | SIESTA (KMG) | Your Method (QPD/OMM) | Winner |
| :--- | :--- | :--- | :--- |
| **Mathematical Basis** | Variational Energy Penalty | Projective Dynamics | **SIESTA** (Theoretically cleaner) |
| **Matrix Operation** | SpGEMM (Matrix-Matrix) | SpMV (Matrix-Vector) | **Yours** (GPU Friendly) |
| **Memory Access** | Random (Fill-in management) | Predictable (Fixed mask) | **Yours** |
| **Pre-factor** | High (Density Matrix expansion) | Low (Direct Orbital Update) | **Yours** |
| **Stability** | Requires Pulay Mixing | Unconditionally Stable (Implicit) | **Yours** |
| **Accuracy** | Controlled by Basis + $\eta$ | Controlled by Basis + Cutoff | **Tie** |

**Verdict:**
Your method is a **modernized, GPU-centric evolution** of the ideas presented in SIESTA. While SIESTA was optimized for the CPUs of 2002 (where minimizing FLOPs was key), your method optimizes for **data locality and parallelism** (avoiding SpGEMM and fill-in), which is the correct strategy for 2025 hardware.

---

## USER 5

Back to Projective dynamics.
you egain expless it in term of force
"Compute Constraint Force" v_i += a*F_i  where F_i = sum_j( <v_i|v_j> ) v_j
Although I agree that it is kinda force I would prefer to formulate it in terms of Jacobi or Gauss-Seidel update.
Also what I do not right in your approach that you suggest to compute full vector  F_i and then truncate it my Mask only on support of v_i
However this is inefficient and it would in fact lead already to O(n^2*k) scaling as each vector F_i has full domenision n, we have n such vectos, and aach such vector F_i is obtained by sum of n k-dimensions vectors v_j 

Instead what I want is to update the components v_i variationally, that is ve compute variation of coefs v_i,mu only for mu within support of v_i with respect to the overlaps E_ij = <v_i,v_j> to minimize them, we will never even compute the elements outside of support of v_i

---

You are absolutely right. My previous explanation relied on a "lazy" mathematical shorthand (truncate after compute) which is indeed computationally wasteful.

To achieve strict $O(N)$ scaling and GPU efficiency, we must formulation the problem as **Projected Gradient Descent** (or Jacobi/Gauss-Seidel) directly on the **packed coefficients** $c_{i\mu}$.

Here is the strictly sparse, variational formulation.

### The Variational Formulation

We define the Constraint Energy (Orthogonality Error):
$$ E_{orth} = \frac{k}{2} \sum_{i} \sum_{j \in \text{neigh}(i)} \left( \Omega_{ij} - \delta_{ij} \right)^2 $$
where $\Omega_{ij} = \langle \psi_i | S | \psi_j \rangle$.

We want to find the update $\Delta c_{i\mu}$ for the coefficients of orbital $i$, **only for indices $\mu$ allowed by the mask of $i$**, to minimize this energy.

The gradient of this energy with respect to a specific allowed coefficient $c_{i\mu}$ is:
$$ \frac{\partial E_{orth}}{\partial c_{i\mu}} = \sum_{j \in \text{neigh}(i)} 2 (\Omega_{ij} - \delta_{ij}) \underbrace{\frac{\partial \Omega_{ij}}{\partial c_{i\mu}}}_{(S \psi_j)_\mu} $$

### The Efficient Algorithm (No Full Vectors)

To evaluate $(S \psi_j)_\mu$ efficiently without generating full vectors, we utilize the fact that $S$ is sparse and $\psi_j$ is sparse.

Let's define the **Interaction Vector** $\tilde{\psi}_j = S \psi_j$.
*   Since $\psi_j$ has finite support, $\tilde{\psi}_j$ also has finite support (slightly larger, expanded by the bandwidth of $S$).
*   We only store the non-zero elements of $\tilde{\psi}_j$.

#### The Algorithm Steps (Jacobi Style for GPU)

1.  **Precompute Interaction Vectors:**
    For each orbital $j$, compute $\tilde{\psi}_j = S \psi_j$.
    *   *Implementation:* Sparse Matrix-Vector multiplication restricted to the support.
    *   *Result:* A sparse vector (list of index/value pairs).

2.  **Compute Overlaps (Small Matrix):**
    For each pair of neighbors $(i, j)$:
    $$ \Omega_{ij} = \psi_i \cdot \tilde{\psi}_j $$
    *   *Implementation:* Dot product of two sparse vectors.
    *   *Cost:* Proportional to the intersection of supports.

3.  **Compute Variational Update (The Core):**
    We iterate over orbitals $i$. We want to update its coefficients $c_{i\mu}$.
    
    $$ \Delta c_{i\mu} = -\alpha \sum_{j \in \text{neigh}(i)} (\Omega_{ij} - \delta_{ij}) \times (\tilde{\psi}_j)_\mu $$

    *   **Crucial Efficiency Detail:** We do **not** loop over $\mu$ and then sum $j$. We loop over neighbors $j$, and "splat" (accumulate) their contribution into the coefficients of $i$.
    
    *   *Loop structure:*
        ```python
        Initialize delta_C[i] = 0 (size of mask i)
        For j in neighbors(i):
            scalar_residual = (delta_ij - Omega_ij)
            
            # We need to add: residual * S_psi_j
            # But ONLY at indices overlapping with mask(i)
            
            overlap_indices = intersect(mask(i), mask(j_expanded))
            
            for mu in overlap_indices:
                delta_C[i][mu] += scalar_residual * S_psi_j[mu]
        ```

4.  **Apply Update:**
    $c_{i\mu} \leftarrow c_{i\mu} + \Delta c_{i\mu}$

This approach guarantees we **never** touch a memory address corresponding to a zero coefficient.

### Python Prototype (Strictly Sparse)

Here is the implementation using explicit index mapping to simulate the "Packed Array" structure you would use on a GPU.

```python
import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt

def make_1d_system_sparse(N=50):
    # Standard 1D tight-binding system
    diag = np.ones(N) * -2.0
    off = np.ones(N-1) * -1.0
    H = sp.diags([off, diag, off], [-1, 0, 1], format='csr')
    
    s_diag = np.ones(N) * 1.0
    s_off = np.ones(N-1) * 0.2
    S = sp.diags([s_off, s_diag, s_off], [-1, 0, 1], format='csr')
    return H, S

class VariationalOrbitalSolver:
    def __init__(self, N, n_orb, radius=3, dt=0.1):
        self.N = N
        self.n_orb = n_orb
        self.dt = dt
        
        # 1. Define Geometry (The Masks)
        # indices[i] holds the global basis indices allowed for orbital i
        self.indices = []
        centers = np.linspace(0, N-1, n_orb, dtype=int)
        for c in centers:
            start = max(0, c - radius)
            end = min(N, c + radius + 1)
            self.indices.append(np.arange(start, end, dtype=np.int32))
            
        # 2. Packed Coefficients
        # coeffs[i] holds the values corresponding to indices[i]
        self.coeffs = []
        for i in range(n_orb):
            size = len(self.indices[i])
            vec = np.random.uniform(-1, 1, size)
            vec /= np.linalg.norm(vec)
            self.coeffs.append(vec)

    def compute_S_psi(self, S):
        """
        Compute S * psi_i for all i. 
        Result is stored as (indices, values) because S expands the support slightly.
        """
        S_psi_list = []
        for i in range(self.n_orb):
            # 1. Expand current packed psi to a temporary sparse vector
            # In C++/CUDA we would do SpMV directly. Here we use scipy.
            
            # Create sparse vector for psi_i
            rows = np.zeros(len(self.indices[i])) # flat
            cols = self.indices[i]
            data = self.coeffs[i]
            psi_sparse = sp.csr_matrix((data, (rows, cols)), shape=(1, self.N))
            
            # 2. Multiply: S_psi = psi * S (since S is symmetric)
            # Result is a sparse vector
            res = psi_sparse.dot(S)
            
            # Store simply as indices and data
            res_coo = res.tocoo()
            S_psi_list.append((res_coo.col, res_coo.data))
            
        return S_psi_list

    def step_variational_constraints(self, S, n_iter=5, alpha=0.5):
        """
        Jacobi-style update of coefficients to minimize orthogonality error.
        strictly respecting the mask.
        """
        for k in range(n_iter):
            # A. Precompute S*psi (The Interaction Vectors)
            # These are sparse vectors.
            S_psi_data = self.compute_S_psi(S)
            
            # B. Compute Overlap Matrix Omega (Small, N_orb x N_orb)
            # Only compute for neighbors
            Omega = np.zeros((self.n_orb, self.n_orb))
            
            for i in range(self.n_orb):
                # Unpack i (conceptually)
                idx_i = self.indices[i]
                val_i = self.coeffs[i]
                
                # Loop neighbors (naive loop for demo, use neighbor list in prod)
                for j in range(max(0, i-4), min(self.n_orb, i+5)):
                    # Dot product: psi_i . S_psi_j
                    # We utilize the precomputed S_psi_j
                    idx_j_S, val_j_S = S_psi_data[j]
                    
                    # Fast sparse dot: intersection of idx_i and idx_j_S
                    # Numpy's intersect1d is slow, but in C++ this is linear scan
                    # Here we map to global for the demo simplicity
                    
                    # Trick: use a dense accum for the dot product on the support of i
                    # (Since support is small, this is fast)
                    accum = 0.0
                    
                    # We need dot(psi_i, S_psi_j)
                    # We iterate the smaller set.
                    # Usually S_psi_j is larger than psi_i.
                    
                    # Create a lookup for i (or hash map)
                    # For demo: simple masked dot
                    # (In production: sorted arrays intersection)
                    
                    # Let's do it the "dense-in-mask" way for python speed
                    # Map S_psi_j onto the domain of i
                    mask_overlap = np.isin(idx_j_S, idx_i)
                    if not np.any(mask_overlap): continue
                    
                    # Global indices of overlap
                    overlap_glob = idx_j_S[mask_overlap]
                    overlap_vals = val_j_S[mask_overlap]
                    
                    # Map global to local index in i
                    # (This search is the cost we pay for sparsity in Python, 
                    # but on GPU with fixed stride it's O(1))
                    sorter = np.argsort(idx_i)
                    local_idx = sorter[np.searchsorted(idx_i, overlap_glob, sorter=sorter)]
                    
                    accum = np.dot(val_i[local_idx], overlap_vals)
                    Omega[i, j] = accum

            # C. Variational Update (Jacobi)
            # delta_c_i = sum_j (delta_ij - Omega_ij) * (S_psi_j)_masked_to_i
            
            updates = []
            
            for i in range(self.n_orb):
                update_accum = np.zeros_like(self.coeffs[i])
                idx_i = self.indices[i]
                
                for j in range(max(0, i-4), min(self.n_orb, i+5)):
                    residual = (1.0 if i==j else 0.0) - Omega[i,j]
                    
                    if abs(residual) < 1e-6: continue
                    
                    # Add residual * S_psi_j to update
                    # But ONLY components that fall into idx_i
                    
                    idx_j_S, val_j_S = S_psi_data[j]
                    
                    # Intersect again
                    mask_overlap = np.isin(idx_j_S, idx_i)
                    overlap_glob = idx_j_S[mask_overlap]
                    overlap_vals = val_j_S[mask_overlap]
                    
                    sorter = np.argsort(idx_i)
                    local_idx = sorter[np.searchsorted(idx_i, overlap_glob, sorter=sorter)]
                    
                    # Accumulate
                    update_accum[local_idx] += residual * overlap_vals
                
                updates.append(update_accum)
            
            # Apply updates
            for i in range(self.n_orb):
                self.coeffs[i] += alpha * updates[i]

    def step_physics(self, H, dt):
        """Standard inertial descent."""
        for i in range(self.n_orb):
            # Compute H * psi restricted to mask
            # For simplicity in python: expand -> mult -> restrict
            # In production: SpMV restricted
            rows = np.zeros(len(self.indices[i]))
            cols = self.indices[i]
            data = self.coeffs[i]
            psi_sparse = sp.csr_matrix((data, (rows, cols)), shape=(1, self.N))
            
            force_glob = psi_sparse.dot(-H).toarray().flatten()
            
            # Update only on mask
            self.coeffs[i] += dt * force_glob[self.indices[i]]

# --- Run ---
N = 100
n_orb = 20
H, S = make_1d_system_sparse(N)
solver = VariationalOrbitalSolver(N, n_orb, radius=4)

print("Running strictly sparse variational solver...")
for step in range(100):
    solver.step_physics(H, dt=0.05)
    solver.step_variational_constraints(S, n_iter=5, alpha=0.5)
    
    # Check orthogonality
    # (Just for visualization, this part is not sparse optimized)
    total_err = 0.0
    S_psi = solver.compute_S_psi(S)
    for i in range(n_orb):
        # dot with self
        idx_S, val_S = S_psi[i]
        # map to i
        mask = np.isin(idx_S, solver.indices[i])
        # ... calculation omitted for brevity ...
```

### Why this is better
1.  **Memory Access:** We iterate `j` (neighbors), calculate the scalar residual $r_{ij}$, and then "scatter" the vector $r_{ij} \times \tilde{\psi}_j$ onto the accumulator of $\psi_i$.
2.  **No Zeros:** We never compute, store, or load a zero value.
3.  **GPU Mapping:**
    *   Each orbital $i$ is a thread block.
    *   The coefficients $c_{i\mu}$ are in shared memory.
    *   We stream neighbors $j$.
    *   We stream the indices of $\tilde{\psi}_j$.
    *   If index matches, we atomicAdd to shared memory.
4.  **Scaling:** This scales as $O(N_{orb} \times N_{neigh} \times \text{SupportSize})$. All terms are constant relative to system size $N$.

---

## USER 6

when you say this
```
1.  **Precompute Interaction Vectors:**
    For each orbital $j$, compute $\tilde{\psi}_j = S \psi_j$.
    *   *Implementation:* Sparse Matrix-Vector multiplication restricted to the support.
```

I'm not sure this is correct approach, maybe it goes to far in atempt to prevent "fill-in"

I think to compute CSC matrix, it is efficient to compute SC (that is vectors v_j expanded to nearest neighbor sites by action of S), which lead to some fill in by one layer. and then search such coefs C (only for valied components within allowed support) which minimize the error

Still in your description I do not see exactly the formalism used in Jacobi and GS method ... I mean we 
have matrix A x = b
and update
(sum_(j!=i) A_ij)/A_ii
can you write it like that, what is A, what is b, what is x ?


---

You are absolutely right. Calculating $\tilde{\psi}_j = S \psi_j$ (which expands the support by one "shell" of neighbors) is the standard and correct way to handle this. It is efficient because you perform SpMV (Sparse Matrix-Vector) once per step.

Here is the formulation of the problem exactly as a linear system $Ax=b$ for the **unknown updates to the coefficients**, solved via Jacobi/Gauss-Seidel.

### The System: $A x = b$

We are looking for the update vectors $\Delta \mathbf{c}_i$ (strictly confined to the mask of orbital $i$) that minimize the orthogonality error.

Let's focus on a **single coefficient** $c_{i\mu}$ (the coefficient of basis function $\mu$ inside orbital $i$).
We treat this scalar as our unknown variable $x$.

The orthogonality constraints involving orbital $i$ are:
$$ \langle \psi_i | S | \psi_j \rangle = \delta_{ij} \quad \forall j \in \text{neigh}(i) $$

Linearizing this constraint with respect to the change $\Delta c_{i\mu}$:
$$ \text{Error}_{ij}^{new} \approx \text{Error}_{ij}^{old} + \frac{\partial \text{Error}_{ij}}{\partial c_{i\mu}} \Delta c_{i\mu} $$

The derivative term is simply the $\mu$-th component of the vector $S \psi_j$:
$$ \frac{\partial}{\partial c_{i\mu}} \langle \psi_i | S | \psi_j \rangle = (S \psi_j)_\mu $$

To find the optimal update (Gauss-Newton step), we minimize the sum of squared errors. This leads to the following linear equation for the update $\Delta c_{i\mu}$:

$$ \underbrace{\left( \sum_{j \in \text{neigh}(i)} [(S \psi_j)_\mu]^2 \right)}_{A_{\mu\mu}} \underbrace{\Delta c_{i\mu}}_{x} = \underbrace{-\sum_{j \in \text{neigh}(i)} (\Omega_{ij} - \delta_{ij}) (S \psi_j)_\mu}_{b} $$

---

### The Mapping

*   **The Unknown ($x$):**
    $\Delta c_{i\mu}$ — The change to the coefficient of basis atom $\mu$ in orbital $i$.
    *(Note: We only solve this for $\mu \in \text{Mask}_i$. If $\mu$ is outside, $x=0$).*

*   **The Matrix Diagonal ($A_{ii}$ or $A_{\mu\mu}$):**
    The "stiffness" or sensitivity of the constraints to this specific coefficient.
    $$ A_{\mu\mu} = \sum_{j} (\tilde{\psi}_{j})_\mu^2 $$
    *(Where $\tilde{\psi}_j = S \psi_j$ is the vector expanded by one shell).*

*   **The RHS ($b$):**
    The projection of the current constraint errors onto the direction of the gradient.
    $$ b = - \sum_{j} (\text{Overlap Error}_{ij}) \times (\tilde{\psi}_{j})_\mu $$

---

### The Algorithm (Jacobi / Gauss-Seidel)

Here is the strictly sparse algorithm. It requires **one layer of fill-in** for the auxiliary vectors $\tilde{\psi}$, but the coefficients $C$ remain strictly confined.

#### 1. Preparation (The "Fill-in" Step)
Compute $\tilde{\Psi} = S C$.
*   This creates vectors $\tilde{\psi}_j$ which are non-zero on `Support(j) + Neighbors`.
*   **Cost:** $N_{orb}$ SpMV operations.

#### 2. Compute Errors
Compute the small $N_{orb} \times N_{orb}$ error matrix:
$$ R_{ij} = (\psi_i \cdot \tilde{\psi}_j) - \delta_{ij} $$

#### 3. The Solver Loop (Iterating over Orbitals $i$)

For each orbital $i$, we update its active coefficients. We can do this in parallel (Jacobi) or sequentially (Gauss-Seidel).

**For a specific basis index $\mu$ inside `Mask(i)`:**

1.  **Gather Gradient ($b$):**
    Accumulate the error weighted by the interaction vector.
    $$ b_\mu = - \sum_{j \in \text{neigh}(i)} R_{ij} \times (\tilde{\psi}_j)_\mu $$
    *(Note: We only need $(\tilde{\psi}_j)_\mu$. Even though $\tilde{\psi}_j$ is expanded, we only access index $\mu$ if $\mu$ is in our mask).*

2.  **Gather Diagonal ($A$):**
    Accumulate the squared sensitivity.
    $$ A_{\mu\mu} = \sum_{j \in \text{neigh}(i)} [(\tilde{\psi}_j)_\mu]^2 $$

3.  **Update ($x$):**
    $$ c_{i\mu} \leftarrow c_{i\mu} + \frac{b_\mu}{A_{\mu\mu}} $$

---

### Python Prototype (Correct Linear Algebra Form)

This code implements the update exactly as $x = b/A$ per coefficient.

```python
import numpy as np
import scipy.sparse as sp

def solve_orthogonality_jacobi(C_list, S, masks, n_iter=5):
    """
    C_list: List of vectors (coefficients) for each orbital.
    S:      Sparse overlap matrix.
    masks:  List of indices (global) where C_list[i] is non-zero.
    """
    n_orb = len(C_list)
    
    for iteration in range(n_iter):
        
        # 1. PRE-COMPUTE INTERACTION VECTORS (S * C)
        # This includes the "fill-in" by 1 neighbor shell
        S_psi_list = [] # List of (indices, values) tuples
        
        for i in range(n_orb):
            # Construct sparse vector from current C
            # (In C++ this is just a pointer to the arrays)
            rows = np.zeros(len(masks[i]))
            cols = masks[i]
            data = C_list[i]
            psi = sp.csr_matrix((data, (rows, cols)), shape=(1, S.shape[0]))
            
            # SpMV: Expands support
            s_psi = psi.dot(S).tocoo()
            S_psi_list.append((s_psi.col, s_psi.data))

        # 2. COMPUTE ERROR MATRIX R_ij
        # R_ij = <psi_i | S | psi_j> - delta_ij
        R = np.zeros((n_orb, n_orb))
        
        # We need to access S_psi fast, let's map them for demo
        # (In production, use dense-in-mask arrays or hash maps)
        S_psi_maps = []
        for idx, val in S_psi_list:
            # Create a lookup table for the expanded vector
            # Only needed for the overlap calculation
            mapping = dict(zip(idx, val))
            S_psi_maps.append(mapping)

        for i in range(n_orb):
            for j in range(n_orb): # Neighbor loop optimization omitted
                # Dot product: sum( C_i[k] * S_psi_j[k] )
                # Iterate over mask of i (it's strictly smaller)
                dot = 0.0
                for k_local, k_global in enumerate(masks[i]):
                    if k_global in S_psi_maps[j]:
                        dot += C_list[i][k_local] * S_psi_maps[j][k_global]
                
                R[i,j] = dot - (1.0 if i==j else 0.0)

        # 3. JACOBI UPDATE COEFFICIENTS
        # Update C_i,mu based on the linear system A x = b
        
        new_C_list = []
        
        for i in range(n_orb):
            # We want to update every coefficient mu in mask[i]
            current_coeffs = C_list[i].copy()
            
            # For every basis function mu allowed for orbital i:
            for k_local, mu in enumerate(masks[i]):
                
                numerator_b = 0.0
                denominator_A = 0.0
                
                # Sum over all constraints j that affect this orbital
                for j in range(n_orb):
                    # We need the value (S*psi_j)_mu
                    # Check if S*psi_j is non-zero at mu
                    val_S_psi_j_mu = S_psi_maps[j].get(mu, 0.0)
                    
                    if val_S_psi_j_mu == 0.0: continue
                    
                    # Accumulate RHS (b)
                    # b = - sum_j (Error_ij * interaction_j_mu)
                    numerator_b -= R[i, j] * val_S_psi_j_mu
                    
                    # Accumulate Diagonal (A)
                    # A = sum_j (interaction_j_mu ^ 2)
                    denominator_A += val_S_psi_j_mu**2
                
                # Apply Jacobi Update: x = b / A
                # Damping factor 0.5 helps stability
                if denominator_A > 1e-12:
                    delta_x = numerator_b / denominator_A
                    current_coeffs[k_local] += 0.5 * delta_x
            
            new_C_list.append(current_coeffs)
            
        C_list = new_C_list
        
    return C_list
```

### Why this matches your intuition
1.  **$S \psi$ (Expansion):** We explicitly compute the expanded interaction vector. This correctly captures how a change in coefficient $c_{i\mu}$ affects the overlap with neighbors.
2.  **No coefficient fill-in:** We only loop over `mu` inside the original `masks[i]`. We read data from the expanded $S\psi$ vectors, but we only write to the restricted $C$ vectors.
3.  **$Ax=b$ Form:** The update `delta_x = numerator_b / denominator_A` is exactly the solution to the linearized least-squares problem for that single coefficient.

---

## USER 7

OK, in analogy with projective dynamics, I'm not sure about few things. 
1) I'm missing the interial term. In projective dynamics this terms help to stabilize the algorithm (makes it diagonally dominant), it gives some mixing / weighting of soft external forces with respect to hard constrains. Perhaps we can do something like this here as well to our benefit.
2) We should make clear if we minimize overlap or overlap square Sum_i Sum_j(<v_i|v_j>-delta_ij) or Sum_i Sum_j(<v_i|v_j>-delta_ij)^2 ?
the thing is that if we minimize only Sum_i Sum_j(<v_i|v_j>-delta_ij) it can lead to non-orthogonal vectors with opposite direction - i.e. large negative overlaps, rahter than zero overlap (orthogonality)
3) can you write the final script in terms of loops over indexes, rather then using some math libraries, this will allow us later to implement it in aby language (C/C++/OpenGL, Julia), I don't care the language, write it in python if you like, but write explicit loops and indexes, like old-fashioned, so it is clear. It will also help me relate better to my JAcobi / Projective dynamics algorithm

---

This is the correct direction. Let's formalize it exactly as you requested: **Quantum Projective Dynamics (QPD)** with an inertial term, minimizing the **squared error**, implemented using **low-level explicit loops**.

### 1. The Inertial Term
In Projective Dynamics, the system solves for the position $x$ that balances the "prediction" (inertia) against the "constraints".

$$ E_{total} = \underbrace{\frac{w}{2} \sum_{i,\mu} (c_{i\mu} - \tilde{c}_{i\mu})^2}_{\text{Inertial (Soft)}} + \underbrace{\frac{k}{2} \sum_{i \le j} (\Omega_{ij} - \delta_{ij})^2}_{\text{Constraint (Hard)}} $$

*   $\tilde{c}$: The predicted coefficient (from Hamiltonian descent).
*   $c$: The final corrected coefficient.
*   $w$: The inertial weight (equivalent to mass/$dt^2$).
*   $k$: Constraint stiffness.

The inertial term $w$ is added to the **diagonal** of your linear system. This guarantees the system is solvable (diagonally dominant) even if a coefficient doesn't participate in any constraints.

### 2. The Objective Function
We definitely minimize the **Squared Error**.
$$ \text{Minimize } \sum (\langle \psi_i | S | \psi_j \rangle - \delta_{ij})^2 $$
If we didn't square it, the derivative would be constant, and the optimizer would just push coefficients to $-\infty$ to make the overlap negative. Squaring creates a parabolic well around 0 (or 1).

### 3. The Low-Level Implementation (C-style Python)

Here is the complete logic written with explicit indices, accumulators, and arrays. This translates 1:1 to C++, CUDA, or GLSL Compute Shaders.

I define the data structures as flat arrays (buffers) to mimic how you would store them in memory.

```python
import numpy as np

# ==========================================
# C-STYLE DATA STRUCTURES (Simulating Raw Memory)
# ==========================================

# We assume Compressed Sparse Row (CSR) for Matrices H and S
# val: float array of non-zero values
# col: int array of column indices
# ptr: int array of row start offsets (size N+1)

# We assume Localized Orbitals (MOs) are stored as a Ragged Array
# C_val:     All active coefficients packed in one long array
# C_ind:     The global basis index for each coefficient
# C_ptr:     Start index of orbital i in C_val
# C_len:     Number of active coefficients for orbital i

def c_style_projective_dynamics(
    # System Size
    n_basis, 
    n_orbs,
    
    # Hamiltonian (Sparse CSR)
    H_val, H_col, H_ptr,
    
    # Overlap (Sparse CSR)
    S_val, S_col, S_ptr,
    
    # Orbitals (Ragged Array) - These are Input/Output
    C_val, C_ind, C_ptr, C_len,
    
    # Parameters
    dt,             # Time step
    w_inertial,     # Inertial weight (stabilizer)
    k_constraint,   # Constraint stiffness
    n_iter          # Jacobi iterations
):
    
    # ==========================================
    # STEP 1: PHYSICS PREDICTION (Inertial Step)
    # pred_C = C - dt * H * C
    # ==========================================
    
    # We store predicted values in a separate buffer (like p_pred in PD)
    pred_C_val = np.zeros_like(C_val)
    
    # Loop over all orbitals 'i'
    for i in range(n_orbs):
        
        # Loop over active coefficients 'mu' of orbital 'i'
        # Range in the packed array
        start_i = C_ptr[i]
        end_i   = start_i + C_len[i]
        
        for p in range(start_i, end_i):
            mu_global = C_ind[p] # Global basis index
            
            # We need to compute (H * C_i)[mu]
            # But H is sparse. H * C is hard to do efficiently this way.
            # Instead, we do Gradient Descent: C_new = C_old - dt * grad(E)
            # grad(E) w.r.t c_mu is roughly (H * C)_mu
            
            # Let's perform Sparse Matrix-Vector product explicitly: H[row_mu] . C[i]
            
            h_c_dot = 0.0
            
            # Row 'mu_global' in H matrix
            h_start = H_ptr[mu_global]
            h_end   = H_ptr[mu_global+1]
            
            for h_idx in range(h_start, h_end):
                col_h = H_col[h_idx] # The column index in H
                val_h = H_val[h_idx]
                
                # Now we need the coefficient C_i[col_h]
                # We have to search for 'col_h' inside the current orbital 'i'
                # In C/CUDA, this linear search is fast because support is small (e.g. 10-20 items)
                
                coef_val = 0.0
                for search_k in range(start_i, end_i):
                    if C_ind[search_k] == col_h:
                        coef_val = C_val[search_k]
                        break
                
                h_c_dot += val_h * coef_val
                
            # Apply Inertial update
            pred_C_val[p] = C_val[p] - dt * h_c_dot

    # Copy predicted back to C_val to serve as the target for constraints
    # In PD: x_curr is initialized to x_pred
    for k in range(len(C_val)):
        C_val[k] = pred_C_val[k]


    # ==========================================
    # STEP 2: PRECOMPUTE INTERACTION VECTORS (S * C)
    # We need these to evaluate constraints efficiently.
    # SC_val will hold expanded vectors.
    # ==========================================
    
    # Note: S*C usually has more non-zeros than C. 
    # For this demo, we assume we use a simplified dynamic list or pre-allocated max buffer.
    # Here, I simulate a ragged array construction.
    
    SC_vals = [] # List of lists (simulating ragged)
    SC_inds = [] 
    
    for i in range(n_orbs):
        # We want vector v = S * C_i
        # We iterate over rows of S that match active C_i indices? No.
        # We iterate over C_i non-zeros and add columns of S.
        
        # Dense accumulator for the support region (simulating local memory/cache)
        # We don't know the exact extent, so we use a map or large buffer.
        # For clarity: Dictionary = Hash Map
        accum_map = {} 
        
        start_i = C_ptr[i]
        end_i   = start_i + C_len[i]
        
        for p in range(start_i, end_i):
            mu_global = C_ind[p]
            c_val     = C_val[p]
            
            # Add column 'mu_global' of S scaled by c_val
            # But S is CSR (Row major). S is symmetric, so Col 'mu' == Row 'mu'
            
            s_row_start = S_ptr[mu_global]
            s_row_end   = S_ptr[mu_global+1]
            
            for s_idx in range(s_row_start, s_row_end):
                target_idx = S_col[s_idx]
                s_elem     = S_val[s_idx]
                
                if target_idx not in accum_map:
                    accum_map[target_idx] = 0.0
                accum_map[target_idx] += s_elem * c_val
        
        # Store result for orbital i
        # Convert map to arrays
        keys = sorted(accum_map.keys()) # Sorted for efficient intersection later
        vals = [accum_map[k] for k in keys]
        SC_inds.append(keys)
        SC_vals.append(vals)

    # ==========================================
    # STEP 3: JACOBI SOLVER FOR CONSTRAINTS
    # Minimize: w*(c - pred)^2 + k*(Overlap - delta)^2
    # ==========================================
    
    for iter_idx in range(n_iter):
        
        # A. Compute Overlap Matrix Omega_ij = C_i . (S*C_j)
        # -----------------------------------------------
        # We can treat this as small dense matrix
        Omega = np.zeros((n_orbs, n_orbs))
        
        for i in range(n_orbs):
            for j in range(i, n_orbs): # Symmetric loop
                
                # Check bounding boxes (omitted here, but essential in C++)
                
                dot_sum = 0.0
                
                # Intersection of C[i] and SC[j]
                # Since SC[j] is the "Interaction Vector", dot(C_i, SC_j) == <i|S|j>
                
                # Linear scan intersection (C-style)
                # Assumes C_ind and SC_inds are sorted
                ptr_a = C_ptr[i]
                end_a = C_ptr[i] + C_len[i]
                
                idx_b_list = SC_inds[j] # This is the S*C vector for j
                val_b_list = SC_vals[j]
                len_b = len(idx_b_list)
                ptr_b = 0
                
                while ptr_a < end_a and ptr_b < len_b:
                    idx_a = C_ind[ptr_a]
                    idx_b = idx_b_list[ptr_b]
                    
                    if idx_a < idx_b:
                        ptr_a += 1
                    elif idx_b < idx_a:
                        ptr_b += 1
                    else:
                        # Match
                        dot_sum += C_val[ptr_a] * val_b_list[ptr_b]
                        ptr_a += 1
                        ptr_b += 1
                
                Omega[i, j] = dot_sum
                Omega[j, i] = dot_sum

        # B. Jacobi Update
        # Solve A_diag * dx = b_residual
        # -----------------------------------------------
        
        # Buffer for updates (Jacobi requires not overwriting C immediately)
        C_updates = np.zeros_like(C_val)
        
        for i in range(n_orbs):
            
            # Loop over active coefficients 'mu' of orbital 'i'
            start_i = C_ptr[i]
            end_i   = start_i + C_len[i]
            
            for p in range(start_i, end_i):
                mu_global = C_ind[p]
                c_curr    = C_val[p]
                c_target  = pred_C_val[p]
                
                # --- Build Linear System (Diagonal Approximation) ---
                
                # 1. Inertial part
                # Energy = 0.5 * w * (c - pred)^2
                # dE/dc = w * (c - pred)
                # Diagonal Contribution = w
                # RHS Contribution      = w * (pred - c)
                
                A_diag = w_inertial
                b_rhs  = w_inertial * (c_target - c_curr)
                
                # 2. Constraint part
                # Energy = 0.5 * k * sum_j (<i|S|j> - delta)^2
                # Grad   = k * sum_j (<i|S|j> - delta) * (S*C_j)_mu
                
                # We assume Diagonal Hessian approx (Gauss-Newton):
                # A_diag += k * sum_j [ (S*C_j)_mu ]^2
                # RHS    -= Grad
                
                for j in range(n_orbs):
                    # Only if i and j are neighbors (Omega_ij non-zero)
                    # or strictly if bounding boxes overlap.
                    
                    # We need the value (S*C_j) at index 'mu'
                    # Linear Search in SC_inds[j]
                    sc_val_mu = 0.0
                    
                    # (Optimization: In C, use hash map or binary search)
                    idx_list = SC_inds[j]
                    val_list = SC_vals[j]
                    for k_scan in range(len(idx_list)):
                        if idx_list[k_scan] == mu_global:
                            sc_val_mu = val_list[k_scan]
                            break
                    
                    if sc_val_mu == 0.0: continue
                    
                    # Current Error
                    target_overlap = 1.0 if i == j else 0.0
                    error = Omega[i, j] - target_overlap
                    
                    # Add to System
                    # J^T J approximation for Hessian diagonal
                    A_diag += k_constraint * (sc_val_mu * sc_val_mu)
                    
                    # J^T * error for RHS
                    b_rhs -= k_constraint * error * sc_val_mu
                
                # 3. Solve dx
                # (A + epsilon) dx = b
                d_c = b_rhs / (A_diag + 1e-12)
                
                C_updates[p] = d_c
        
        # Apply updates
        for k in range(len(C_val)):
            C_val[k] += C_updates[k]

    return C_val
```

### Key Differences from Previous Code:

1.  **Inertial Term (`w_inertial`):** It is explicit now.
    *   `A_diag = w_inertial + ...`
    *   `b_rhs = w_inertial * (pred - curr) ...`
    This keeps the "soft" physics (Hamiltonian descent) balanced with the "hard" constraints (Orthogonality).
    *   If `w` is high, orbitals stick to the Physics solution.
    *   If `k` is high, orbitals stick to Orthogonality.

2.  **Squared Error Minimization:**
    The gradient term `b_rhs -= k * error * sc_val_mu` combined with the Hessian term `A_diag += k * sc_val_mu^2` is mathematically the **Gauss-Newton step** for minimizing the sum of squares.
    $$ \Delta x = -\frac{\nabla E}{\nabla^2 E} \approx -\frac{J^T r}{J^T J} $$

3.  **Strict Loops:**
    The intersection logic (`ptr_a`, `ptr_b`) and the linear searches mimic exactly how you would write this in C or a Compute Shader kernel where you don't have dictionaries.

This should now fully align with your mental model of Projective Dynamics: Predict the state, then iteratively relax it towards the constraint manifold using a localized linear solve.


--- 

USER 8

WHAT THE HELL IS THIS ?!?!?!

        # We can treat this as small dense matrix
        Omega = np.zeros((n_orbs, n_orbs))

we are aiming for order n method and you are computing dense matrix with n^2 elements? Why the fuck you think this matrix is true, it is neither small nor dense. 

You are like missing the key point, or you forgot, overlooked.
since the finite support of the moleculer orbitals only some orbitals overlap. We need to have some neighbor-list for the orbitals to know which volerpa which can be precalculated before hand just form geometry - position of site on which are centered and the cutoff R_cut

Also what I do not like is that your matrix operations are all in one block. I would prefere if you have function like "SparseMatrixVecttor" and "SparseJacobi" where you explicitly work with neighbolist of the sites. Since both H-matrix and S-matrix can be writen as some cooefs, and some indexis of neighbor sites.

Do not confuse neighbor list os sites (orbitals given by cutoff of basis functions) and neighbor list of molecular orbitals (given by finite cutoff/support of the molecular orbitals)

---

You are absolutely right to be angry. Using `np.zeros((n_orbs, n_orbs))` was a lazy shortcut in the python prototype that violates the fundamental premise of Order-N.

If we have 10,000 orbitals, that matrix is 100 million entries (800 MB) and computing it is $10^8$ operations. For linear scaling, we must **only** store and compute overlaps for physically overlapping orbitals.

Here is the strictly **Sparse, Order-N, Component-based** implementation.

I have separated the logic into low-level kernels as you requested:
1.  `sparse_mv_ragged`: Applies a sparse matrix (H or S) to the localized orbitals.
2.  `solve_constraints_sparse`: The Projective Dynamics solver using an **Orbital Neighbor List**.

### Data Structures Definition

To be crystal clear, here are the data structures we assume (C-style):

*   **Sparse Matrix (H, S)**: Standard CSR or Adjacency list.
    *   `val[k], col[k]`, `row_ptr[i]` (Standard CSR).
*   **Localized Orbitals (C)**: Ragged Array.
    *   `C_val`: Packed array of all active coefficients.
    *   `C_ptr[i]`: Start index in `C_val` for orbital `i`.
    *   `C_ind[k]`: The physical site index (basis function index) for coefficient `k`.
*   **Orbital Neighbor List**: The graph of interacting orbitals.
    *   `Orb_neigh_ptr[i]`: Start index for neighbors of orbital `i`.
    *   `Orb_neigh_col[k]`: The index of the neighbor orbital.

---

### The Code (Explicit Loops)

```python
import numpy as np

# =========================================================
# LOW LEVEL KERNEL 1: Sparse Matrix - Ragged Vector Product
# Computes Y = M * C
# This handles H * C (Physics) or S * C (Interaction Vectors)
# Complexity: O(N_orb * Support_Size * Avg_Connectivity_of_M)
# =========================================================
def sparse_mv_ragged(
    n_orbs,
    # Matrix M (CSR format)
    M_val, M_col, M_ptr,
    # Input Ragged Vectors C
    C_val, C_ind, C_ptr, C_len,
    # Output Ragged Vectors Y (Topology must match C, or be pre-expanded)
    # For H*C, topology is usually same as C if H is local.
    # For S*C, topology expands by 1 shell. We assume Y is pre-allocated.
    Y_val, Y_ind, Y_ptr 
):
    # Loop over every orbital 'i'
    for i in range(n_orbs):
        
        # Range of active coefficients for orbital i in INPUT
        c_start = C_ptr[i]
        c_end   = c_start + C_len[i]
        
        # Range of active coefficients for orbital i in OUTPUT
        y_start = Y_ptr[i]
        y_end   = y_start + C_len[i] # Assuming same sparsity for demo
        
        # For every component 'mu' in the output vector Y_i
        for p_out in range(y_start, y_end):
            row_idx = Y_ind[p_out] # Physical site index (row of M)
            
            # Dot product: Row(M, row_idx) . Vector(C_i)
            dot_val = 0.0
            
            # Iterate row of sparse matrix M
            m_start = M_ptr[row_idx]
            m_end   = M_ptr[row_idx+1]
            
            for m_k in range(m_start, m_end):
                col_idx = M_col[m_k] # Column index in M
                val_m   = M_val[m_k]
                
                # Find C_i[col_idx]
                # Linear search in C_i (fast because support is small ~20-50 atoms)
                # In GPU this is efficient coalesced read if sorted
                c_val_k = 0.0
                for p_in in range(c_start, c_end):
                    if C_ind[p_in] == col_idx:
                        c_val_k = C_val[p_in]
                        break
                
                if c_val_k != 0.0:
                    dot_val += val_m * c_val_k
            
            Y_val[p_out] = dot_val


# =========================================================
# LOW LEVEL KERNEL 2: Sparse Jacobi Constraint Solver
# Solves for dC to minimize: w|C-Pred|^2 + k|Overlap-I|^2
# Complexity: O(N_orb * N_orb_neighs * Support_Size)
# =========================================================
def solve_constraints_sparse(
    n_orbs,
    # Orbital Neighbor List (Who overlaps with whom)
    Orb_neigh_list, # List of lists for simplicity here, or CSR
    
    # Interaction Vectors (Precomputed SC = S * C)
    SC_val, SC_ind, SC_ptr, SC_len,
    
    # Current State (Input/Output)
    C_val, C_ind, C_ptr, C_len,
    
    # Target (Predicted) State
    Pred_val,
    
    # Parameters
    w_inertial, k_constraint
):
    
    # Loop over orbitals 'i' (Parallelizable on GPU)
    for i in range(n_orbs):
        
        # Identify active components of orbital i
        start_i = C_ptr[i]
        end_i   = start_i + C_len[i]
        
        # Iterate over every coefficient 'mu' of orbital 'i'
        for p in range(start_i, end_i):
            mu = C_ind[p]
            
            # --- BUILD LOCAL LINEAR SYSTEM (Diagonal Approx) ---
            # A * dx = b
            
            # 1. Inertial contribution
            # Energy = 0.5*w*(c - pred)^2
            A_diag = w_inertial
            b_rhs  = w_inertial * (Pred_val[p] - C_val[p])
            
            # 2. Constraint contribution
            # Loop ONLY over neighboring orbitals 'j'
            # (Strictly O(1) number of neighbors)
            neighbors = Orb_neigh_list[i] 
            
            for j in neighbors:
                
                # A. Calculate current Pairwise Overlap Omega_ij
                # Omega_ij = C_i . SC_j
                # (SC_j is the precomputed S*C_j vector)
                
                omega_ij = 0.0
                
                # Intersection loop (C_i AND SC_j)
                # Linear scan assuming sorted indices
                ptr_a = C_ptr[i]
                end_a = ptr_a + C_len[i]
                
                ptr_b = SC_ptr[j]
                end_b = ptr_b + SC_len[j]
                
                while ptr_a < end_a and ptr_b < end_b:
                    idx_a = C_ind[ptr_a]
                    idx_b = SC_ind[ptr_b]
                    
                    if idx_a < idx_b:
                        ptr_a += 1
                    elif idx_b < idx_a:
                        ptr_b += 1
                    else:
                        omega_ij += C_val[ptr_a] * SC_val[ptr_b]
                        ptr_a += 1
                        ptr_b += 1
                
                # B. Constraint Gradient Terms
                # Error = Omega_ij - delta_ij
                target = 1.0 if i == j else 0.0
                error  = omega_ij - target
                
                # We need the sensitivity of this error w.r.t our coefficient c_{i,mu}
                # d(Omega_ij)/d(c_{i,mu}) = (SC_j)_mu
                
                # Find (SC_j)_mu
                sc_j_mu = 0.0
                for search_k in range(SC_ptr[j], SC_ptr[j] + SC_len[j]):
                    if SC_ind[search_k] == mu:
                        sc_j_mu = SC_val[search_k]
                        break
                
                if sc_j_mu == 0.0: continue
                
                # Add to Linear System (Gauss-Newton)
                # A += k * J^T * J
                # b -= k * J^T * error
                
                term = k_constraint * sc_j_mu
                
                A_diag += term * sc_j_mu
                b_rhs  -= term * error
            
            # 3. Solve and Apply Update (Jacobi Step)
            if A_diag > 1e-12:
                d_c = b_rhs / A_diag
                C_val[p] += d_c # Can use a temporary buffer for pure Jacobi
                
```

### Explanation of the Sparse Logic

1.  **Complexity**:
    *   The loop structure is `For i in Orbitals` $\to$ `For mu in Support` $\to$ `For j in Neighbor_Orbitals`.
    *   Since `Support` is constant (volume of cutoff) and `Neighbor_Orbitals` is constant (packing density), the inner work is constant.
    *   Total Complexity: $O(N_{orbs})$.

2.  **No Dense Matrices**:
    *   We never form `Omega[N,N]`. We calculate `omega_ij` on the fly inside the neighbor loop.
    *   We never form full vectors. `sparse_mv_ragged` only writes to the pre-allocated sparse structure `Y`.

3.  **The "Projective" Part**:
    *   The term `b_rhs = w * (Pred - C)` pulls the orbital towards the Physics result.
    *   The term `b_rhs -= k * error * sc_j_mu` pulls the orbital towards the Orthogonality Manifold.
    *   `A_diag` weighs these two desires appropriately so the step is stable.

4.  **Implementation Notes**:
    *   **Pre-calculation**: You must pre-calculate `Orb_neigh_list`. This is easy: Iterate all orbitals, look at their site indices. If orbital $i$ and orbital $j$ share any site (or site-neighbor), add them to the list.
    *   **S*C expansion**: In `sparse_mv_ragged`, if you compute $S \times C$, the resulting vector has a slightly larger support than $C$. You should handle this by pre-allocating `Y` (the `SC` vector) with the expanded sparsity pattern (neighbors of neighbors of centers).


    ---

    ## USER 9

    Still I would prefer, resp. I in fact insist, if you extract function you would call "SparseJacobiStep" and "SparseDotMatrixVector" or somethig like that. Soemthing in alanlogy what I have in my projective dynamics

```C++
    void ProjectiveDynamics_d::updateJacobi_fly( Vec3d* ps_in, Vec3d* ps_out ){
        // This solver calculates bi and Aii on-the fly, preventing updatePD_RHS call and possibly improving numerical stability (?)
        // NOTE: this is goes beyond linar solution of Ap=b because we update b every iteration, and use non-linear operations like d_{ij}=(pi-pj)/|pi-pj| 
        for(int i=0; i<nPoint; i++){
            const Vec3d    pi = ps_in[i];
            const double idt2 = 1.0 / (dt * dt);
            double         Ii = points[i].w*idt2; 
            if(kFix)       Ii += kFix[i];
            Vec3d sum_j; sum_j.set_mul(pi, Ii );   // b_i    =  M_i/dt^2 p'_i   +  \sum_j ( K_{ij} d_{ij} )         
            double Aii = Ii;                       // A_{ii} =  M_i/dt^2        +  \sum_j   K_{ij} 
            const int j0 = i * nNeighMax;
            for(int jj=0; jj<nNeighMax; jj++){
                const int j  = j0 + jj;
                const int jG = neighs[j];
                if(jG == -1) break;
                
                //   update RHS bi
                const Quat4d par = params[j];
                const Vec3d  pj  = ps_in[jG];
                const Vec3d  dij = pi - pj;
                const double l   = dij.norm();
                const double k   = par.z;
                sum_j.add_mul( dij, k*par.x/l );   //  b_i  +=  \sum_j ( K_{ij} d_{ij} )   
                sum_j.add_mul( pj , k         );   //  s_j  +=  \sum_j ( K_{ij} p_j    )
                Aii  += k;

            }
            sum_j.mul( 1/Aii );
            ps_out[i]  = sum_j;
        }
    }

    void ProjectiveDynamics_d::updateGaussSeidel_fly( Vec3d* ps ){
        // This solver calculates bi and Aii on-the fly, preventing updatePD_RHS call and possibly improving numerical stability (?)
        // NOTE: this is goes beyond linar solution of Ap=b because we update b every iteration, and use non-linear operations like d_{ij}=(pi-pj)/|pi-pj| 
        for(int i=0; i<nPoint; i++){
            const Vec3d    pi = ps[i];
            const double idt2 = 1.0 / (dt * dt);
            double         Ii = points[i].w*idt2; 
            if(kFix)       Ii += kFix[i];
            Vec3d sum_j; sum_j.set_mul(pi, Ii );   // b_i    =  M_i/dt^2 p'_i   +  \sum_j ( K_{ij} d_{ij} )         
            double Aii = Ii;                       // A_{ii} =  M_i/dt^2        +  \sum_j   K_{ij} 
            const int j0 = i * nNeighMax;
            for(int jj=0; jj<nNeighMax; jj++){
                const int j  = j0 + jj;
                const int jG = neighs[j];
                if(jG == -1) break;
                
                //   update RHS bi
                const Quat4d par = params[j];
                const Vec3d  pj  = ps[jG];
                const Vec3d  dij = pi - pj;
                const double l   = dij.norm();
                const double k   = par.z;
                sum_j.add_mul( dij, k*par.x/l );   //  b_i  +=  \sum_j ( K_{ij} d_{ij} )   
                sum_j.add_mul( pj , k         );   //  s_j  +=  \sum_j ( K_{ij} p_j    )
                Aii  += k;

            }
            sum_j.mul( 1/Aii );
            ps[i]  = sum_j;
        }
    }

    void ProjectiveDynamics_d::updateIterativeMomentum( Vec3d* psa, Vec3d* psb ){
        //Vec3d* psa = ps_in;
        //Vec3d* psb = ps_out;
        //if(LinSolveMethod::JacobiMomentum == (LinSolveMethod)linSolveMethod  ){
        updatePD_RHS(psa, bvec ); // TODO: we need to calculate this only whem we use linear solvers, not if we use _fly solvers
        //}
        for (int i = 0; i < nSolverIters; i++) {  
            switch( (LinSolveMethod)linSolveMethod ){
                case LinSolveMethod::JacobiMomentum:    { updateJacobi_lin( psa, psb, bvec ); } break;
                case LinSolveMethod::JacobiFlyMomentum: { updateJacobi_fly( psa, psb );       } break;
                case LinSolveMethod::GSMomentum:     { 
                    for (int j=0; j<nPoint; j++){ psb[j]=psa[j]; }
                    updateGaussSeidel_lin( psb, bvec ); 
                } break; 
                case LinSolveMethod::GSFlyMomentum:     { 
                    for (int j=0; j<nPoint; j++){ psb[j]=psa[j]; }
                    updateGaussSeidel_fly( psb ); 
                } break; 
            }
            //updateJacobi_lin( psa, psb, bvec );      //    p'_{k+1}    = Solve(A,p_k,b)
            double bmix = mixer.get_bmix( i );       //    bmix_k = bmix(k)
            if( (i==0) || (i>=(nSolverIters-1) ) ) bmix = 0.0; // make sure we stop momentum mixing for the last iteration      
            for(int j=0; j<nPoint; j++){ 
                Vec3d p = psb[j]; 
                p.add_mul( linsolve_yy[j], bmix );   //    p_{k+1} = p'_k + bmix d_k
                linsolve_yy[j] = p - psa[j];         //    d_{k+1} = p_{k+1} - p_k
                psa[j]         = p;
            }   // we use linsolve_yy to store momentum
            //for (int i=0; i<nPoint; i++){ ps_pred[i]=ps_cor[i]; }
        }

        if( bApplyResudualForce ){
            double dt2 = dt * dt;
            evalTrussForce( psa, linsolve_b );
            //updateJacobi_lin( psa, 0, bvec, linsolve_b );
            // for (int i=0; i<nPoint; i++){ 
            //     if(kFix){ if(kFix[i]>1e-8) continue; }
            //     psa[i].add_mul( linsolve_b[i], residualForceFactor*dt2/points[i].w ); 
            // }
        }

        for (int i=0; i<nPoint; i++){ psb[i]=psa[i]; }
    }
```



---

## Gemini

This is exactly the right analogy. In your C++ code, `updateJacobi_fly` calculates the diagonal $A_{ii}$ and the force $b_i$ accumulation **on the fly** inside the loop to avoid storing a global matrix.

We will do exactly the same for the Quantum solver.

Here are the mappings:
*   **Particle $i$** $\rightarrow$ **Orbital $i$**
*   **Position $\vec{p}_i$** $\rightarrow$ **Coefficients vector $\mathbf{c}_i$**
*   **Spring Constraint** $\rightarrow$ **Orthogonality Constraint** $\langle \psi_i | \psi_j \rangle = 0$
*   **Bond Vector $\vec{d}_{ij}$** $\rightarrow$ **Interaction Vector** $\tilde{\psi}_j = S \psi_j$
*   **Mass $M_i$** $\rightarrow$ **Inertial Weight $w$**

Here is the implementation using strictly explicit loops and indices, designed to be ported 1:1 to your C++ engine.

### Data Structures (C-Style)

We assume **Ragged Arrays** for orbitals.
*   `C_val`, `C_ind` (basis index), `C_ptr` (start index per orbital).
*   `SC_val`, `SC_ind`, `SC_ptr` (Pre-calculated $S \times C$ vectors).
*   `Neigh_ind`, `Neigh_ptr` (List of neighboring orbitals).

---

### 1. `sparse_mv_ragged` (Precompute Interaction Vectors)
This corresponds to preparing the geometry. In PD this might be implicit, but here we need to compute $\tilde{\psi} = S \psi$ because the "direction" of the constraint force depends on the overlap matrix $S$.

```python
def sparse_mv_ragged(
    n_orbs,
    # Matrix S (CSR)
    S_val, S_col, S_ptr,
    # Input C (Ragged)
    C_val, C_ind, C_ptr, C_len,
    # Output SC (Ragged - pre-allocated with slightly wider mask)
    SC_val, SC_ind, SC_ptr
):
    """
    Computes SC = S * C using explicit sparse loops.
    Equivalent to updating bond vectors in PD, but linear.
    """
    for i in range(n_orbs):
        
        # Loop over active components 'mu' in the OUTPUT vector SC_i
        # (The mask of SC must be wide enough to cover S*C)
        sc_start = SC_ptr[i]
        sc_end   = SC_ptr[i+1] # Assuming packed format for simplicity
        
        for p_out in range(sc_start, sc_end):
            row_mu = SC_ind[p_out] # Physical basis site
            
            # Perform dot product: Row(S, row_mu) . Vector(C_i)
            sum_val = 0.0
            
            # Iterate row of S
            s_start = S_ptr[row_mu]
            s_end   = S_ptr[row_mu+1]
            
            for k in range(s_start, s_end):
                col_nu = S_col[k] # Neighbor basis site
                val_s  = S_val[k]
                
                # Find coefficient C_i[col_nu]
                # Linear scan in C_i (fast for small support)
                c_start = C_ptr[i]
                c_end   = c_start + C_len[i]
                
                c_val = 0.0
                for p_in in range(c_start, c_end):
                    if C_ind[p_in] == col_nu:
                        c_val = C_val[p_in]
                        break
                
                sum_val += val_s * c_val
            
            SC_val[p_out] = sum_val
```

---

### 2. `updateJacobi_fly`
This is the direct equivalent of your C++ function. It iterates over orbitals, and for each active coefficient (component), it builds the Diagonal ($A$) and RHS ($b$) by summing over constraints (neighbors).

```python
def update_jacobi_fly(
    n_orbs,
    # Topology
    Neigh_ind, Neigh_ptr, # Orbital Neighbor list
    
    # Data (Ragged)
    C_val, C_ind, C_ptr, C_len,        # Current State
    Pred_val,                          # Predicted State (Inertial target)
    SC_val, SC_ind, SC_ptr, SC_len,    # Interaction Vectors (S*C)
    
    # Output
    C_out_val, # Buffer for new coefficients
    
    # Physics Parameters
    w_inertial,    # Like Mass/dt^2
    k_constraint   # Stiffness
):
    
    # Loop over "Particles" (Orbitals)
    for i in range(n_orbs):
        
        # 1. Pre-calculate scalar Overlaps with neighbors for this orbital
        #    This is like calculating bond lengths |p_i - p_j| once per pair
        #    We store them in a temporary stack buffer
        
        neigh_start = Neigh_ptr[i]
        neigh_end   = Neigh_ptr[i+1]
        n_neigh     = neigh_end - neigh_start
        
        # Temp buffer for overlaps: omega[k] corresponds to neighbor k
        # In C++: double overlaps[64];
        overlaps = [0.0] * n_neigh 
        
        for k in range(n_neigh):
            j = Neigh_ind[neigh_start + k]
            
            # Compute Dot: <C_i | SC_j>
            # Linear intersection of sorted indices
            dot = 0.0
            p_a = C_ptr[i]
            end_a = p_a + C_len[i]
            p_b = SC_ptr[j]
            end_b = p_b + SC_len[j]
            
            while p_a < end_a and p_b < end_b:
                idx_a = C_ind[p_a]
                idx_b = SC_ind[p_b]
                
                if idx_a < idx_b: p_a += 1
                elif idx_b < idx_a: p_b += 1
                else: 
                    dot += C_val[p_a] * SC_val[p_b]
                    p_a += 1; p_b += 1
            
            # Constraint Error: (<v_i|v_j> - delta_ij)
            target = 1.0 if i == j else 0.0
            overlaps[k] = dot - target


        # 2. Update Components (The Fly Step)
        #    Loop over active components 'mu' of orbital 'i'
        
        c_start = C_ptr[i]
        c_end   = c_start + C_len[i]
        
        for p in range(c_start, c_end):
            mu_global = C_ind[p]
            
            # A. Inertial Term
            # A_ii = M / dt^2
            # b_i  = M / dt^2 * p_pred
            
            A_diag = w_inertial
            b_rhs  = w_inertial * Pred_val[p]
            
            # B. Constraint Terms (Sum over neighbors)
            # b_i += k * error * direction
            # A_ii += k * direction^2
            
            for k in range(n_neigh):
                j = Neigh_ind[neigh_start + k]
                
                # We need component 'mu' of vector SC[j]
                # "direction" = (S * psi_j)_mu
                
                direction = 0.0
                
                # Linear search in SC[j] for index 'mu'
                # (Fast because support is small)
                p_scan = SC_ptr[j]
                end_scan = p_scan + SC_len[j]
                
                for scan in range(p_scan, end_scan):
                    if SC_ind[scan] == mu_global:
                        direction = SC_val[scan]
                        break
                
                if direction == 0.0: continue
                
                # Accumulate Force (RHS) and Stiffness (LHS)
                # Note: overlaps[k] is the scalar error calculated above
                
                b_rhs  -= k_constraint * overlaps[k] * direction
                A_diag += k_constraint * direction * direction
            
            # C. Solve Linear System (Scalar division)
            # ps_out[i] = sum_j / Aii
            
            C_out_val[p] = b_rhs / (A_diag + 1e-12)
```

---

### 3. `updateGaussSeidel_fly`
Just like in your C++ code, GS is identical to Jacobi except it reads and writes to the same array (`ps` instead of `ps_in`/`ps_out`), allowing faster convergence.

```python
def update_gauss_seidel_fly(
    n_orbs,
    Neigh_ind, Neigh_ptr,
    C_val, C_ind, C_ptr, C_len,  # In-Place Update
    Pred_val,
    SC_val, SC_ind, SC_ptr, SC_len,
    w_inertial, k_constraint
):
    # Loop over Orbitals
    for i in range(n_orbs):
        
        # 1. Pre-calculate overlaps (Same as Jacobi)
        # BUT using current C_val which might have been updated by previous 'i'
        neigh_start = Neigh_ptr[i]
        n_neigh     = Neigh_ptr[i+1] - neigh_start
        overlaps    = [0.0] * n_neigh 
        
        for k in range(n_neigh):
            j = Neigh_ind[neigh_start + k]
            # ... Dot product logic (C_i . SC_j) ...
            # ... (Copied from Jacobi) ...
            dot = 0.0
            p_a = C_ptr[i]; end_a = p_a + C_len[i]
            p_b = SC_ptr[j]; end_b = p_b + SC_len[j]
            while p_a < end_a and p_b < end_b:
                # ... same intersection logic ...
                idx_a = C_ind[p_a]; idx_b = SC_ind[p_b]
                if idx_a < idx_b: p_a += 1
                elif idx_b < idx_a: p_b += 1
                else: 
                    dot += C_val[p_a] * SC_val[p_b]
                    p_a += 1; p_b += 1
            
            target = 1.0 if i == j else 0.0
            overlaps[k] = dot - target

        # 2. Update Components In-Place
        c_start = C_ptr[i]
        c_end   = c_start + C_len[i]
        
        for p in range(c_start, c_end):
            mu_global = C_ind[p]
            
            A_diag = w_inertial
            b_rhs  = w_inertial * Pred_val[p]
            
            for k in range(n_neigh):
                j = Neigh_ind[neigh_start + k]
                
                # Find (SC[j])_mu
                direction = 0.0
                p_scan = SC_ptr[j]; end_scan = p_scan + SC_len[j]
                for scan in range(p_scan, end_scan):
                    if SC_ind[scan] == mu_global:
                        direction = SC_val[scan]
                        break
                
                if direction == 0.0: continue
                
                b_rhs  -= k_constraint * overlaps[k] * direction
                A_diag += k_constraint * direction * direction
            
            # WRITE BACK DIRECTLY TO C_val
            C_val[p] = b_rhs / (A_diag + 1e-12)
```

### Why this structure matches your request
1.  **Analogy:**
    *   `A_diag` calculation $\leftrightarrow$ `Aii += k`
    *   `b_rhs` accumulation $\leftrightarrow$ `sum_j.add_mul(...)`
    *   Loop structure is identical.
2.  **Order N:** All inner loops (neighbors, linear search in support) are constant size with respect to total system size.
3.  **No Allocation:** It uses fixed buffers. The `overlaps` array is tiny (size of neighbors, ~10-20) and can be stack allocated in C++.


--- 

## USER 10


I'm not sure if I like your logic

# Linear intersection of sorted indices
 while p_a < end_a and p_b < end_b:
                # ... same intersection logic ...

it seems to me overly complex and inefficient. I looks like we are searching the index rather just using some precalculated index. Am I right? Cannt we replace this by some precalculated indexes? E.g. inverse mapping. We did not need this stuff in my jacobi_fly functions where we had the neighborlists 

How yould you implement this on GPU

notice we very much care about performance, the fact that we write in slow python maybe mask it but final target is C++ and GPU (OpenCL). Would it help you to write in C to be more performance aware?


---

You are absolutely right. In high-performance computing (and especially on GPUs), **memory lookups are cheap, but branching (searching/while-loops) is expensive.**

Your intuition is correct: We should not search for indices on the fly. We should **pre-calculate the pointers** exactly as you do with neighbor lists in MD/PD.

In MD, you have a list of pairs `(i,j)`.
Here, we need a list of **component pairs**. For every coefficient $\mu$ of orbital $i$, and for every neighbor orbital $j$, we need to know exactly **where** in memory the corresponding component of $j$ is stored.

### The Solution: The "Interaction Map" (Pre-calculated Indexing)

We will replace the "intersection search" with a direct lookup table.

#### Data Structure
We introduce a new integer array: `InterMap`.
*   **Dimensions:** `(Total_Number_Of_Active_Coeffs) x (Max_Neighbors)`
*   **Value:** `InterMap[coeff_index][n]` stores the memory address (integer index) in the `SC` array corresponding to the $n$-th neighbor. If they don't share that basis function, store `-1`.

This allows us to write the kernel with **zero searching** and **zero branching** (except checking for -1).

Here is the implementation in **C-style** (using Python syntax for the logic) to show exactly how this maps to GPU/OpenCL.

---

### 1. The Setup (CPU / Pre-computation)
This happens once (or only when neighbor lists change).

```python
def build_interaction_map(
    n_orbs,
    C_ptr, C_len, C_ind,      # Orbital Topology
    Neigh_ptr, Neigh_ind,     # Orbital-Orbital Neighbor List
    SC_ind, SC_ptr, SC_len    # Topology of interaction vectors
):
    # Allocate Map: Flat array
    # Size = (Total Active Coeffs) * (Max Neighbors per Orbital)
    # We use a fixed max_neigh for GPU friendliness (padding with -1)
    MAX_NEIGH = 20 # Example constant
    
    total_coeffs = C_ptr[n_orbs-1] + C_len[n_orbs-1]
    
    # Initialize with -1 (null pointer)
    # Dimensions: [Total_Coeffs][MAX_NEIGH]
    InterMap = np.full((total_coeffs, MAX_NEIGH), -1, dtype=np.int32)
    
    for i in range(n_orbs):
        
        # Get neighbors of i
        n_start = Neigh_ptr[i]
        n_cnt   = Neigh_ptr[i+1] - n_start
        
        # Loop over active coefficients (components) of i
        c_start = C_ptr[i]
        c_cnt   = C_len[i]
        
        for k in range(c_cnt):
            p_c = c_start + k    # Memory address of this coefficient
            mu  = C_ind[p_c]     # Global basis index
            
            # Loop over neighbors of i
            for n in range(min(n_cnt, MAX_NEIGH)):
                j = Neigh_ind[n_start + n]
                
                # We need to find where 'mu' is located in SC[j]
                # We do the search HERE (on CPU, once), so GPU doesn't have to.
                
                ptr_sc = -1
                
                # Linear scan of neighbor's support 
                # (This is the slow part we removed from the hot loop)
                sc_j_start = SC_ptr[j]
                sc_j_end   = sc_j_start + SC_len[j]
                
                for scan in range(sc_j_start, sc_j_end):
                    if SC_ind[scan] == mu:
                        ptr_sc = scan
                        break
                
                # Store the direct pointer
                InterMap[p_c, n] = ptr_sc
                
    return InterMap
```

---

### 2. The Solver Kernel (GPU / C-Style)
This is now a strictly linear stream of operations. 

*   **Thread Mapping:** 1 Thread per Orbital (Coarse) or 1 Thread per Coefficient (Fine).
*   **Logic:** `update_jacobi_fly_fast`

```python
def update_jacobi_fly_fast(
    n_orbs,
    # Topology
    Neigh_ptr, 
    InterMap,     # PRE-CALCULATED LOOKUP [Total_Coeffs x MAX_NEIGH]
    MAX_NEIGH,
    
    # Data
    C_val, C_ptr, C_len,
    Pred_val,
    SC_val,       # We only need values! Indices are baked into InterMap
    
    # Output
    C_out,
    
    # Params
    w_inertial, k_constraint
):
    # PARALLEL LOOP (In C++: #pragma omp parallel for or __kernel)
    for i in range(n_orbs):
        
        # ---------------------------------------------------
        # PHASE 1: Calculate Scalar Overlaps Omega_ij
        # ---------------------------------------------------
        
        # Local stack array for overlaps (like bond lengths)
        overlaps = np.zeros(MAX_NEIGH)
        
        c_start = C_ptr[i]
        c_cnt   = C_len[i]
        
        # Loop over my coefficients
        for k in range(c_cnt):
            p_my = c_start + k
            val_my = C_val[p_my]
            
            # Loop over neighbors to accumulate dot product
            # Unrolled loop on GPU
            for n in range(MAX_NEIGH):
                
                # DIRECT LOOKUP: Where is this component in neighbor n?
                p_neigh = InterMap[p_my, n]
                
                # Branchless check (in CUDA use predicate)
                if p_neigh != -1:
                    # Omega_ij += C_i[mu] * SC_j[mu]
                    overlaps[n] += val_my * SC_val[p_neigh]
        
        # Apply delta_ij (Target overlap)
        # We assume the first neighbor in list is always 'self' (i=i) for simplicity
        # or we check index. Here assuming overlaps[0] is self.
        overlaps[0] -= 1.0 
        
        # ---------------------------------------------------
        # PHASE 2: Update Coefficients
        # ---------------------------------------------------
        
        for k in range(c_cnt):
            p_my = c_start + k
            
            # Inertial Term
            A_diag = w_inertial
            b_rhs  = w_inertial * Pred_val[p_my]
            
            # Constraint Terms
            for n in range(MAX_NEIGH):
                p_neigh = InterMap[p_my, n]
                
                if p_neigh != -1:
                    # Access neighbor value directly
                    val_neigh = SC_val[p_neigh] # This is (S * psi_j)_mu
                    
                    # Accumulate
                    err = overlaps[n]
                    term = k_constraint * val_neigh
                    
                    A_diag += term * val_neigh
                    b_rhs  -= term * err
            
            # Write result
            C_out[p_my] = b_rhs / (A_diag + 1e-12)

```

### Why this is GPU optimal

1.  **Memory Coalescing:** `p_my` increments linearly. `InterMap` is accessed linearly (row-major). `SC_val` access is gathered, but that is unavoidable in sparse physics (texture cache handles this well).
2.  **No Searching:** The `while` loop intersection is gone. Replaced by `InterMap[p_my, n]`.
3.  **No Index Arrays in Kernel:** Notice we don't even pass `C_ind` or `SC_ind` or `Neigh_ind` to the kernel! The `InterMap` encodes the entire topology.
4.  **Register Usage:** We reuse `overlaps[]` which sits in registers.

### Is `InterMap` too big?

*   Let's check: 
    *   `N_atoms` = 1000
    *   `N_coeffs` (Total) $\approx 50 \times 1000 = 50,000$
    *   `MAX_NEIGH` $\approx 30$
    *   Total Map Size = $50,000 \times 30 = 1.5 \text{ Million ints}$.
    *   **Size = 6 MB**.

This is **tiny**. Even for 100,000 atoms, it is 600 MB, which fits easily on a modern 8GB+ GPU. This trades memory for raw compute speed, which is exactly the right trade-off for Projective Dynamics.

---

## USER 11

wait a moment, it does not seem to tiny. 
Let me understand why and if it is not better to store it differently.

in projective dynamics we have just one dense vector we solve for (particle positions), which has size n (natom,naparticles)
here we instead solve for n (norb=natom) of sparse vectors which are localized onky on k sites (k=50) right?
Then we have max_neigh (for sites - basis set cutoff) not max_neigh_MO (for molecular orbital, which is higher), right?
WE definitely mist distinguish these two

now it seems to me that this reindexing can be factorized into two subsequne indexing - 1) maping between index in dense/sparse vector v_i and indexing of neighbors.

I guess this factorization will save a lot of memory bandwidth, but it will introduce two jumps and wrek cache efficiency. But what we can probably do on GPU is pre-compute the final indexies for just small block inside local memory and reuse them, or somethig like this?

maybe it is tome how to impement this in OpenCL kernels, which expose us to thse memory access patern optimization decisions


---

## USER 12

OK, write the second kenrel

this approach of first pracalculation projection scalars and then variation derivs, remind me fititng be spline coefs to fit data on grid

we also first precompute spline (which is basically projection like SC) and then compute variation

```C++

__attribute__((hot)) 
double getVariations1D( const int n, const double* Gs, const double* Es, const double* Ws, double* fs, double* ps ){
    constexpr double B0=2.0/3.0;
    constexpr double B1=1.0/6.0;
    // --- evaluate current spline
    for(int i=0; i<n; i++){ fs[i]=0; ps[i]=0;  }
    // --- evaluate current spline (in order to evelauet approximation error)
    double err2sum=0;
    for(int i=1; i<n-1; i++){
        double val = Gs[i-1]*B1 + Gs[i]*B0 + Gs[i+1]*B1;
        double err = Es[i] - val;
        if(Ws){ err*=Ws[i]; }
        err2sum += err*err;
        ps[i] = err;
        //Ws[j] = err;
    }
    // --- distribute variational derivatives of approximation error
    for(int i=0; i<n; i++){
        if( (i>0)&&(i<(n-1)) ) [[likely]] {
            double val = ps[i-1]*B1 + ps[i]*B0 + ps[i+1]*B1;
            fs[i] = val;
        }else{
            fs[i] = 0;
        }
    }
    return err2sum;
}

__attribute__((hot)) 
int fit1D( const int n, double* Gs,  double* Es, double* Ws, double Ftol, int nmaxiter=100, double dt=0.1, double Kreg=-1.0, bool bPBC=false ){
    //if(verbosity>1)
    printf("Bspline::fit1D() bPBC=%i Kreg=%g Ftol=%g nmaxiter=%i dt=%g \n", bPBC, Kreg, Ftol, nmaxiter, dt );
    const double F2max = Ftol*Ftol;
    double* ps = new double[n];
    double* fs = new double[n];
    double* vs = new double[n];
    bool bRegForce = Kreg>0;
    Quat4i xqis[4];
    if(bRegForce){ make_inds_pbc( n, xqis ); }
    int itr=0;
    for(int i=0; i<n; i++){ vs[i]=0.0; } // clear velocity
    //if( ckeckRange(n, 1, Gs, -1.e+6, +1.e+6, "Gs", true ) && ckeckRange(n, 1, Gs, -1.e+6, +1.e+6, "Gs", true ) ){ printf("ERROR in fit1D()=> exit\n"); exit(0); };
    bool bErr = false;
    double err2sum=0;
    double f2reg=0;
    double f2sum=0;

    int nsub = 20;
    int nij = n*nsub;
    dgb_dBs.resize(nij);
    dgb_dFs.resize(nij);
    dgb_Bs.resize(nij);
    golbal_array_dict.insert( { "dgb_vG",   NDArray{ (double*)dgb_dBs.data(),   Quat4i{nij,4,-1,-1}} }  );
    golbal_array_dict.insert( { "dgb_dF",   NDArray{ (double*)dgb_dFs.data(),   Quat4i{nij,4,-1,-1}} }  );
    golbal_array_dict.insert( { "dgb_Bs",   NDArray{ (double*)dgb_Bs .data(),   Quat4i{nij,4,-1,-1}} }  );

    for(itr=0; itr<nmaxiter; itr++){        
        //if(bPBC){ err2sum = getVariations1D_pbc( n, Gs, Es, Ws, fs, ps ); }
        //else    { err2sum = getVariations1D    ( n, Gs, Es, Ws, fs, ps ); }
        if(bPBC){ err2sum = getVariations1D_pbc( n, Gs, Es, 0, fs, ps ); }
        else    { err2sum = getVariations1D    ( n, Gs, Es, 0, fs, ps ); }
        if( bRegForce ){
            //f2reg = regularizationForceMidpoint_1D( n, Gs, fs, Kreg, xqis );
            for(int i=0; i<n; i++){ Ws[i]=0.0; };
            //f2reg = regularizationForceMidpoint_1D( n, Gs, Ws, Kreg, xqis );
            f2reg = regularizationForceMidpoint_1D_sub( n, nsub, Gs, Ws, Kreg, xqis );
            //for(int i=0; i<n; i++){ fs[i]+=Ws[i]; };
        }
        Vec3d cfv = move(dt,n,Gs,fs, vs );
        //if(verbosity>2)printf( "fit1D |F[%i]|=%g |Err|=%g |Freg|=%g cos(f,v)=%g\n",itr,sqrt(f2sum), sqrt(err2sum), sqrt(f2reg), cfv.x/sqrt(cfv.y*cfv.z) );
        f2sum=cfv.y;
        if(f2sum<F2max){ break; };
    }
    //if(verbosity>1)
    //printf( "Bspline::fit1D() iter=%i err=%g \n", itr, sqrt(err2sum) );
    printf( "fit1D DONE |F[%i]|=%g |Err|=%g |Freg|=%g \n",itr,sqrt(f2sum), sqrt(err2sum), sqrt(f2reg) );
    delete [] ps;
    delete [] fs;
    delete [] vs;
    return itr;
}

```

---

Yes, the analogy is perfect.

In your B-spline code:
1.  `getVariations1D` calculates the **residuals** (errors) `ps[i] = Es[i] - val`.
2.  Then it distributes these residuals back to the control points via the basis functions (`B0`, `B1`).

In our Quantum Projective Dynamics:
1.  **Kernel 1 (This one):** Calculates the **residuals** (orthogonality errors) `Overlaps[i][n] = <psi_i | S psi_j> - delta_ij`.
2.  **Kernel 2 (Previous one):** Distributes these residuals back to the coefficients `C` weighted by the interaction vectors `SC`.

Here is the **Overlap Reduction Kernel**. It uses the same `Atom_Interaction_Map` to ensure memory access is perfectly aligned with the Update kernel.

### OpenCL Kernel: `compute_scalar_overlaps`

**Strategy:**
*   **1 Workgroup = 1 Orbital.**
*   We load the `Atom_Interaction_Map` into Local Memory (just like in the update kernel).
*   We compute the dot products for all neighbors.
*   We use a **Parallel Reduction** in shared memory to sum the partial dot products calculated by different threads.

```c
#define MAX_ATOMS_PER_ORB 64
#define MAX_NEIGHS 32
#define BASIS_SIZE 4
#define WORKGROUP_SIZE 64 // Adjust to hardware (e.g. 64 for AMD, 32 for NV)

__kernel void compute_scalar_overlaps(
    const int n_orbs,
    
    // Topology
    __global const int*   Neigh_Ptr,      
    __global const int*   Neigh_Ind,
    __global const short* Atom_Inter_Map, // The pre-calculated map
    
    // Data
    __global const int*   C_Ptr,
    __global const int*   Atom_Count,
    __global const float* C_Val,
    
    __global const float* SC_Val, // Precomputed S*C
    __global const int*   SC_Ptr,
    
    // Output: Flat array of scalar overlaps
    // Layout: [Orb0_Neigh0, Orb0_Neigh1, ... Orb1_Neigh0, ...]
    __global float*       Out_Overlaps
) {
    int orb_id  = get_group_id(0);
    int lid     = get_local_id(0);
    int lsize   = get_local_size(0);
    
    if (orb_id >= n_orbs) return;

    // ---------------------------------------------------------
    // 1. Load Map to Local Memory (Identical to Update Kernel)
    // ---------------------------------------------------------
    __local short local_map[MAX_ATOMS_PER_ORB * MAX_NEIGHS];
    
    int n_atoms     = Atom_Count[orb_id];
    int n_neighs    = Neigh_Ptr[orb_id+1] - Neigh_Ptr[orb_id];
    int map_offset  = (C_Ptr[orb_id] / BASIS_SIZE) * MAX_NEIGHS;
    
    // Cooperative load
    int total_map_entries = n_atoms * MAX_NEIGHS;
    for (int k = lid; k < total_map_entries; k += lsize) {
        local_map[k] = Atom_Inter_Map[map_offset + k];
    }
    
    barrier(CLK_LOCAL_MEM_FENCE);

    // ---------------------------------------------------------
    // 2. Compute Partial Dot Products
    // ---------------------------------------------------------
    // Each thread computes a partial dot product for ALL neighbors.
    // We store these in a local register array to avoid bank conflicts.
    
    float my_partials[MAX_NEIGHS];
    
    // Initialize registers
    for(int n=0; n<MAX_NEIGHS; n++) my_partials[n] = 0.0f;
    
    int c_start   = C_Ptr[orb_id];
    int n_coeffs  = n_atoms * BASIS_SIZE;

    // Loop over coefficients assigned to this thread
    for (int k = lid; k < n_coeffs; k += lsize) {
        
        int atom_local_idx = k / BASIS_SIZE;
        int basis_offset   = k % BASIS_SIZE;
        
        float val_c = C_Val[c_start + k];
        
        // Accumulate contribution to each neighbor
        for (int n = 0; n < n_neighs; n++) {
            
            // Where is this atom in neighbor n?
            short map_idx = local_map[atom_local_idx * MAX_NEIGHS + n];
            
            if (map_idx != -1) {
                // Get neighbor ID to find SC start
                int neigh_id = Neigh_Ind[Neigh_Ptr[orb_id] + n];
                int sc_start = SC_Ptr[neigh_id];
                
                // Read pre-computed interaction vector component
                int addr = sc_start + (map_idx * BASIS_SIZE) + basis_offset;
                float val_sc = SC_Val[addr];
                
                // Math: <psi_i | S psi_j> = sum( c_i[k] * (S*c_j)[k] )
                my_partials[n] += val_c * val_sc;
            }
        }
    }

    // ---------------------------------------------------------
    // 3. Parallel Reduction in Local Memory
    // ---------------------------------------------------------
    // We need to sum 'my_partials' across all threads in the workgroup
    
    // Re-use local_map memory or allocate new scratchpad? 
    // Let's allocate explicit reduction buffer.
    // Needs [Workgroup_Size * MAX_NEIGHS] floats.
    // 64 * 32 * 4 bytes = 8KB. Safe.
    __local float reduction_buffer[WORKGROUP_SIZE * MAX_NEIGHS];
    
    // Copy registers to shared mem
    for(int n=0; n<n_neighs; n++) {
        reduction_buffer[lid * MAX_NEIGHS + n] = my_partials[n];
    }
    
    barrier(CLK_LOCAL_MEM_FENCE);
    
    // Standard Tree Reduction
    for (int offset = lsize / 2; offset > 0; offset >>= 1) {
        if (lid < offset) {
            for(int n=0; n<n_neighs; n++) {
                reduction_buffer[lid * MAX_NEIGHS + n] += 
                    reduction_buffer[(lid + offset) * MAX_NEIGHS + n];
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    
    // ---------------------------------------------------------
    // 4. Write Result (Thread 0 only)
    // ---------------------------------------------------------
    if (lid == 0) {
        // We assume Out_Overlaps is packed: 
        // [N_Neig_Orb0, N_Neigh_Orb1, ...] 
        // We need an offset array for this output or assume fixed stride.
        // Assuming packed via Neigh_Ptr structure:
        int out_start = Neigh_Ptr[orb_id];
        
        for(int n=0; n<n_neighs; n++) {
            float total_overlap = reduction_buffer[0 * MAX_NEIGHS + n];
            
            // Subtract delta_ij (Orthogonality Constraint)
            int neigh_id = Neigh_Ind[out_start + n];
            if (neigh_id == orb_id) {
                total_overlap -= 1.0f;
            }
            
            // This is the "residual" (ps[i] in your spline code)
            Out_Overlaps[out_start + n] = total_overlap;
        }
    }
}
```

### Key Technical Details

1.  **Register Accumulation:**
    We calculate the dot products for *all 32 neighbors* simultaneously in registers (`my_partials`). This is extremely efficient because `val_c` is loaded once and reused 32 times (broadcast). This maximizes arithmetic density.

2.  **Shared Memory Reduction:**
    Since the dot product is split across threads (each thread handles a subset of coefficients), we must sum them up. The tree reduction (`offset >>= 1`) is the standard GPU pattern for this.

3.  **Memory Alignment:**
    The `local_map` logic is copy-pasted from the Update kernel. This guarantees that if we optimize the map layout for one, the other benefits too.

4.  **Output:**
    The kernel outputs the raw errors $\Omega_{ij} - \delta_{ij}$. The Update kernel will simply load these scalars:
    `float err = Out_Overlaps[neigh_idx];`
    instead of recomputing them or loading from a placeholder.