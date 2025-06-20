# Finding Localized Molecular Orbitals via Energy Minimization and Constraints

## Introduction

In standard Linear Combination of Atomic Orbitals (LCAO) methods, molecular orbitals (MOs) $\psi_i$ are typically found by solving the generalized eigenvalue problem $HC = SCE$, where $H$ is the Hamiltonian matrix, $S$ is the Overlap matrix, $C$ is the matrix of MO coefficients, and $E$ is the diagonal matrix of orbital energies. This approach yields delocalized canonical MOs that are mutually orthogonal ($\langle \psi_i | S | \psi_j \rangle = \delta_{ij}$) and minimize the total electronic energy for a given basis set.

However, for large systems or specific applications like Fragment Molecular Orbitals (FMO) or force fields based on localized orbitals, we might prefer to work with orbitals that are spatially confined. Solving the global eigenvalue problem can be computationally expensive (scaling as $O(N^3)$ with system size $N$), and the resulting orbitals are often delocalized across the entire system.

This document outlines an alternative approach: finding approximative, highly localized orbitals by directly minimizing their energies subject to explicit constraints on orthogonality and localization, without solving the full $HC=SCE$ problem. We assume a tight-binding context where the basis functions $\phi_\mu$ are localized on atomic sites $\mathbf{X}_\mu$, and the Overlap ($S_{\mu\nu} = \langle \phi_\mu | \phi_\nu \rangle$) and Hamiltonian ($H_{\mu\nu} = \langle \phi_\mu | \hat{H} | \phi_\nu \rangle$) matrix elements in this basis are pre-calculated or approximated.

## Orbital Representation

Each molecular orbital $\psi_i$ is expressed as a linear combination of the atomic-like basis functions $\phi_\mu$:
$$ \psi_i = \sum_\mu c_{i,\mu} \phi_\mu $$
where $c_{i,\mu}$ are the expansion coefficients for the $i$-th orbital on the $\mu$-th basis function.

## Energy Minimization Objective

For each occupied orbital $\psi_i$, we want to minimize its energy $E_i$. The energy of an orbital is given by:
$$ E_i = \frac{\langle \psi_i | \hat{H} | \psi_i \rangle}{\langle \psi_i | \psi_i \rangle} $$
In the non-orthogonal basis $\phi_\mu$, this becomes:
$$ E_i = \frac{\sum_\mu \sum_\nu c_{i,\mu}^* c_{i,\nu} H_{\mu\nu}}{\sum_\mu \sum_\nu c_{i,\mu}^* c_{i,\nu} S_{\mu\nu}} $$
If we enforce normalization $\langle \psi_i | \psi_i \rangle = 1$, the objective simplifies to minimizing:
$$ E_i = \sum_\mu \sum_\nu c_{i,\mu}^* c_{i,\nu} H_{\mu\nu} $$
subject to the constraint $\sum_\mu \sum_\nu c_{i,\mu}^* c_{i,\nu} S_{\mu\nu} = 1$.

## Constraints

We impose the following constraints on the orbitals:

### 1. Normalization

Each orbital $\psi_i$ must be normalized:
$$ \langle \psi_i | \psi_i \rangle = \sum_\mu \sum_\nu c_{i,\mu}^* c_{i,\nu} S_{\mu\nu} = 1 $$
This constraint can be enforced explicitly at each step of an iterative optimization procedure. After updating the coefficients $c_{i,\mu}$ for an orbital $\psi_i$, we calculate its current norm squared $S_{ii}^{current} = \sum_\mu \sum_\nu c_{i,\mu}^* c_{i,\nu} S_{\mu\nu}$ and rescale the coefficients:
$$ c_{i,\mu}^{new} = c_{i,\mu}^{old} / \sqrt{S_{ii}^{current}} $$

### 2. Orthogonality

The occupied orbitals $\psi_i$ must be mutually orthogonal with respect to the overlap metric:
$$ \langle \psi_i | \psi_j \rangle = \sum_\mu \sum_\nu c_{i,\mu}^* c_{j,\nu} S_{\mu\nu} = \delta_{ij} \quad \text{for all occupied } i, j $$
Enforcing this constraint while simultaneously minimizing energy and maintaining localization is challenging. Instead of Gram-Schmidt orthogonalization (which is order-dependent), we can use iterative symmetric orthogonalization methods. These methods typically involve calculating the off-diagonal overlaps $S_{ij}$ and applying updates to the coefficients of orbitals $i$ and $j$ to reduce $S_{ij}$. This can be viewed as applying "orthogonalization forces" or updates to the coefficient matrix $C$.

A common approach is based on the Löwdin symmetric orthogonalization idea, but applied iteratively to pairs or blocks of orbitals to drive off-diagonal overlaps towards zero. For a pair of orbitals $\psi_i$ and $\psi_j$, the goal is to reduce $S_{ij}$. The updates to $c_{i,\mu}$ and $c_{j,\mu}$ can be derived from minimizing a penalty term like $S_{ij}^2$ or $S_{ij}$ itself, subject to normalization.

### 3. Localization

To ensure the orbitals are highly localized, we can introduce a localization potential or explicitly constrain coefficients. A common approach is to add a penalty term to the energy functional that increases if an orbital spreads too far from a designated center $\mathbf{x}_i$. This center $\mathbf{x}_i$ could be the position of a nucleus, the centroid of the orbital's density, or a variational parameter itself.

The localization potential $V_{loc}(|\mathbf{r} - \mathbf{x}_i|)$ penalizes contributions from basis functions $\phi_\mu$ whose centers $\mathbf{X}_\mu$ are far from $\mathbf{x}_i$. The penalty term for orbital $\psi_i$ could be of the form:
$$ E_{i,loc} = \langle \psi_i | V_{loc}(|\mathbf{r} - \mathbf{x}_i|) | \psi_i \rangle = \sum_\mu \sum_\nu c_{i,\mu}^* c_{i,\nu} \langle \phi_\mu | V_{loc}(|\mathbf{r} - \mathbf{x}_i|) | \phi_\nu \rangle $$
A simpler approach, often used in tight-binding or localized orbital methods, is to directly penalize coefficients $c_{i,\mu}$ based on the distance $|\mathbf{X}_\mu - \mathbf{x}_i|$. For instance, we could add a term $\sum_\mu |c_{i,\mu}|^2 f(|\mathbf{X}_\mu - \mathbf{x}_i|)$ where $f$ is a function that increases with distance.

Alternatively, and perhaps more simply for finite-support basis functions, we can enforce localization by setting coefficients $c_{i,\mu}$ to zero if the basis function center $\mathbf{X}_\mu$ is beyond a certain cutoff distance $R_{loc}$ from the orbital center $\mathbf{x}_i$:
$$ c_{i,\mu} = 0 \quad \text{if} \quad |\mathbf{X}_\mu - \mathbf{x}_i| > R_{loc} $$
This explicitly limits the spatial extent of the orbital. The orbital center $\mathbf{x}_i$ could be fixed on a specific atom, or it could be a variational parameter optimized alongside the coefficients.

## Variational Optimization

The process of finding the localized orbitals involves an iterative optimization procedure. For each occupied orbital $\psi_i$, we want to update its coefficients $c_{i,\mu}$ to reduce its energy $E_i$ while satisfying the constraints.

1.  **Calculate Energy Gradient:** Compute the variational derivative of the energy $E_i$ with respect to the coefficients $c_{i,\mu}^*$:
    $$ \frac{\partial E_i}{\partial c_{i,\mu}^*} = \sum_\nu c_{i,\nu} H_{\mu\nu} $$
    This gives the direction of steepest descent for the energy in coefficient space (ignoring constraints for a moment).

2.  **Apply Localization Cutoff:** If using the cutoff method, set $c_{i,\mu} = 0$ for all $\mu$ where $|\mathbf{X}_\mu - \mathbf{x}_i| > R_{loc}$. This should ideally be done *before* calculating the gradient or applied as a projection step.

3.  **Incorporate Normalization:** The gradient calculated above is for the unconstrained energy. To account for normalization, we can use projected gradients or simply re-normalize the coefficients after each update step. The latter is simpler: update $c_{i,\mu}$ based on the gradient, then rescale by $1/\sqrt{\langle \psi_i | \psi_i \rangle}$.

4.  **Incorporate Orthogonality:** This is the most complex part. We need to calculate the off-diagonal overlaps $S_{ij}$ for all pairs of occupied orbitals. Based on these overlaps, we derive updates $\Delta c_{i,\mu}$ and $\Delta c_{j,\mu}$ that reduce $S_{ij}$. These updates are then applied simultaneously to all orbitals. This step can be seen as calculating "orthogonalization forces" on the coefficients. The total update for $c_{i,\mu}$ would be a combination of the energy gradient update and the orthogonalization update.

The iterative process would look something like this:
*   Initialize coefficients $c_{i,\mu}$ (e.g., from atomic orbitals or a previous step).
*   Iterate until convergence:
    *   For each occupied orbital $\psi_i$:
        *   Apply localization cutoff: $c_{i,\mu} \leftarrow 0$ if $|\mathbf{X}_\mu - \mathbf{x}_i| > R_{loc}$.
        *   Calculate energy gradient $\frac{\partial E_i}{\partial c_{i,\mu}^*}$.
    *   Calculate overlap matrix $S_{ij}$ for all occupied orbitals.
    *   Calculate orthogonalization updates based on $S_{ij}$ (for $i \neq j$).
    *   Combine energy gradient updates and orthogonalization updates for all $c_{i,\mu}$.
    *   Apply the combined updates to the coefficients.
    *   For each occupied orbital $\psi_i$:
        *   Renormalize coefficients: $c_{i,\mu} \leftarrow c_{i,\mu} / \sqrt{\langle \psi_i | \psi_i \rangle}$.
*   Optionally, optimize the orbital centers $\mathbf{x}_i$ alongside the coefficients.

This approach replaces the global matrix diagonalization with a series of local updates and constraint enforcements, potentially leading to better scaling and the desired localized orbitals.

## Implementation details


### 2.2. Iterative Symmetric Orthogonalization

While solving the generalized eigenvalue problem $HC=SCE$ directly yields a set of orthogonal molecular orbitals, there are scenarios where orbitals might be obtained or refined through other means, such as direct energy minimization under constraints, or when starting from a non-orthogonal set of guess orbitals. In such cases, or if one wishes to avoid the full diagonalization, iterative methods can be employed to enforce orthogonality. These methods must be order-independent, unlike Gram-Schmidt orthogonalization.

A common approach is a Jacobi-like iterative orthogonalization, which simultaneously updates all orbitals to reduce their pairwise overlaps. Let $\mathbf{c}_i$ be the column vector of coefficients for the $i$-th molecular orbital $\psi_i = \sum_\mu c_{i,\mu} \phi_\mu$. The overlap between two molecular orbitals $\psi_a$ and $\psi_b$ is $S_{ab}^{MO} = \mathbf{c}_a^\dagger S \mathbf{c}_b$, where $S$ is the basis function overlap matrix ($S_{\mu\nu}^{basis}$). We want $S_{ab}^{MO} = \delta_{ab}$.

The iterative update rule for the coefficient vector of orbital $\psi_i$ can be expressed as:
$$ \mathbf{c}_i^{(k+1)} = \mathbf{c}_i^{(k)} - c_{damp} \sum_{j \neq i} S_{ij}^{MO,(k)} \mathbf{c}_j^{(k)} $$
After this update, each orbital $\mathbf{c}_i^{(k+1)}$ must be re-normalized:
$$ \mathbf{c}_i^{(k+1)} \leftarrow \frac{\mathbf{c}_i^{(k+1)}}{\sqrt{(\mathbf{c}_i^{(k+1)})^\dagger S \mathbf{c}_i^{(k+1)}}} $$

Here, $S_{ij}^{MO,(k)} = (\mathbf{c}_i^{(k)})^\dagger S \mathbf{c}_j^{(k)}$ is the overlap between orbitals $\psi_i$ and $\psi_j$ at iteration $k$, and $c_{damp}$ is a damping factor, typically between 0 and 1.

**Pythonic Pseudocode:**

```python
import numpy as np

def iterative_orthogonalization(Cs, Suv, c_damp=-0.5, max_iter=2, tol=1e-7):
    norb = len(Cs)
    dC = np.zeros((norb, norb))
    for iteration in range(max_iter):
        # normalize orbitals
        for i in range(norb):
            Sii = sparse_dot(Cs[i],Cs[i], Suv)
            Cs[i] /= np.sqrt(Sii)
        # Calculate orthogonalization updates (like forces)
        Sij_max = 0
        for i in range(norb):
            for j in range(norb):
                if i == j:  continue
                Sij = sparse_dot(Cs[i], Cs[j], Suv)
                Sij_max = max(Sij_max, Sij)
                dC[i,:] = C[j,:]* Sij*c_damp
                dC[j,:] = C[i,:]* Sij*c_damp
        # Apply updates
        for i in range(norb):
            Cs[i] += dC[i,j]
        if Sij_max < tol:  break  # Check for convergence
    return current_coeffs
```

### 2.3. Exploiting Sparsity for Efficient Matrix Operations

In the context of localized orbitals and tight-binding Hamiltonians, both the molecular orbital coefficient vectors $\mathbf{c}_i$ and the basis Hamiltonian $H_{\mu\nu}$ (or overlap $S_{\mu\nu}$) matrices are often sparse.

*   **Localized Orbitals:** The coefficient $c_{i,\mu}$ for orbital $\psi_i$ is non-zero only for basis functions $\phi_\mu$ that are spatially close to the "center" of $\psi_i$. This means $\mathbf{c}_i$ has few non-zero entries.
*   **Sparse Basis Matrices:** In tight-binding, $H_{\mu\nu}$ and $S_{\mu\nu}$ are non-zero only if basis functions $\phi_\mu$ and $\phi_\nu$ are on the same or nearby atomic sites.

This dual sparsity can be leveraged to significantly accelerate calculations such as:
*   **MO Overlaps:** $S_{ij}^{MO} = \mathbf{c}_i^\dagger S \mathbf{c}_j = \sum_{\mu,\nu} c_{i,\mu}^* S_{\mu\nu} c_{j,\nu}$
*   **Orbital Energies:** $E_i = \mathbf{c}_i^\dagger H \mathbf{c}_i = \sum_{\mu,\nu} c_{i,\mu}^* H_{\mu\nu} c_{i,\nu}$
*   **Action of H on $\mathbf{c}_i$ (for gradients):** $(\mathbf{H}\mathbf{c}_i)_\mu = \sum_\nu H_{\mu\nu} c_{i,\nu}$

Instead of performing dense matrix-vector multiplications, one would iterate only over the non-zero elements. For example, to calculate $S_{ij}^{MO}$:
1. Identify the set of basis function indices $\mathcal{N}_i = \{ \mu \mid c_{i,\mu} \neq 0 \}$ for orbital $\psi_i$.
2. Identify the set of basis function indices $\mathcal{N}_j = \{ \nu \mid c_{j,\nu} \neq 0 \}$ for orbital $\psi_j$.
3. The sum becomes: $S_{ij}^{MO} = \sum_{\mu \in \mathcal{N}_i} \sum_{\nu \in \mathcal{N}_j, S_{\mu\nu} \neq 0} c_{i,\mu}^* S_{\mu\nu} c_{j,\nu}$.

This requires storing $\mathbf{c}_i$ in a sparse format or maintaining lists of their non-zero indices, and potentially using neighbor lists for the basis functions to quickly find non-zero $S_{\mu\nu}$ or $H_{\mu\nu}$.

