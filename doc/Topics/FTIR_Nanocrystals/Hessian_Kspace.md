# USER

I have classical forcefield like UFF, which I compute in periodic boundary conditions, now I want to fit the stiffness constants in that forcefield to experimental phonon spectra. How to actually do that. I can easily compute hessian and dynamical matrix by displacing atoms and computing forces (basically it is kinda like correlation matrix of forces by displacemen of two degrees of freedom, on on atom or the other, all possible pairs)
but I need to somehow implement bloch waves, the k-dependnece in that.
Do I need to do it from start of do I only tke one hessian and then apply the bloch k on that hessian ? You know what I mean?

I know how to compute vibrations for molecule, but I'm now a bit strugguling how to do it for crystal with k-space dependence.

---

# Gemini

You are on the exact right track. The short answer is **no, you do not need to recompute the Hessian from scratch for every $k$-point.** Instead, you compute the **real-space Hessian once** for a sufficiently large supercell (or using a cut-off radius appropriate for your UFF forcefield), and then you apply the Bloch phase factors to that static Hessian to construct the **Dynamical Matrix $D(k)$** for any $k$-point you want.

Here is exactly how the physics works and how to implement it step-by-step.

## 1. The Real-Space Hessian (What you already know)

For a molecule, you compute the mass-weighted Hessian. For a periodic crystal, you do almost the same thing, but you must track which periodic image of the unit cell the atoms belong to.

Let's define your reference unit cell as $L = 0$ (the origin cell). You need to compute the forces on atom $i$ in the origin cell ($L=0$) when you displace atom $j$ in cell $L$ (where $L$ is a lattice vector $\vec{R}_L$).

The real-space force constant matrix (Hessian component) is:

$$K_{i\alpha, j\beta}(0, L) = -\frac{\partial F_{i\alpha}(0)}{\partial u_{j\beta}(L)}$$

where:

* $i, j$ are atom indices within a single unit cell.
* $\alpha, \beta$ are Cartesian coordinates ($x, y, z$).
* $L$ represents the lattice vector $\vec{R}_L$ pointing to the cell containing the displaced atom.

> **Crucial Check for UFF:** Because UFF is a short-range forcefield (usually cutting off after a few Angstroms for non-bonded interactions), $K(0, L)$ will quickly become zero for cells $L$ that are far away. You only need to calculate displacements for atoms $j$ in neighboring cells $L$ that fall within your forcefield's interaction cutoff.

## 2. Applying the Bloch Condition (The $k$-dependence)

To get the phonon frequencies at a specific wavevector $\vec{k}$, you project this real-space Hessian into Fourier space using the Bloch theorem.

For a given $\vec{k}$, you build the **$k$-dependent Dynamical Matrix $D(\vec{k})$** by summing over all the lattice vectors $L$ where your force constants are non-zero:

$$D_{i\alpha, j\beta}(\vec{k}) = \frac{1}{\sqrt{m_i m_j}} \sum_{L} K_{i\alpha, j\beta}(0, L) e^{i \vec{k} \cdot \vec{R}_L}$$

where $m_i$ and $m_j$ are the masses of atoms $i$ and $j$.

### What this means visually/computationally:

1. **Size of the matrix:** If your unit cell has $N$ atoms, $D(\vec{k})$ is always a $3N \times 3N$ matrix, regardless of how many periodic images ($L$) you summed over.
2. **Complex Numbers:** Because of the phase factor $e^{i \vec{k} \cdot \vec{R}_L}$, $D(\vec{k})$ is a **complex Hermitian matrix** (except at specific high-symmetry points like $\Gamma$ where the imaginary part vanishes).

## 3. Getting the Phonon Frequencies

Once you have constructed the $3N \times 3N$ matrix $D(\vec{k})$ for your chosen $\vec{k}$:

1. Diagonalize $D(\vec{k})$.
2. The eigenvalues $\lambda$ are the squares of the angular frequencies: $\lambda = \omega^2$.
3. The frequencies are $\nu = \frac{\sqrt{\lambda}}{2\pi}$.

*(Note: If you get a tiny negative eigenvalue near the $\Gamma$ point, those are just your 3 acoustic translations. If you get large negative eigenvalues elsewhere, your structure isn't fully relaxed in the UFF forcefield).*

## 4. How to Fit the Stiffness Constants

Since your ultimate goal is to fit the UFF parameters to experimental phonon spectra (e.g., from inelastic neutron scattering or Raman/IR data), here is how the loop looks:

```
  [ Guess UFF Parameters ]
            │
            ▼
  [ Compute Real-Space Hessian K(0,L) ]  <── (Only needs finite differences once per iteration)
            │
            ▼
     For each experimental k-point:
        ├── Build D(k) using the Bloch sum
        └── Diagonalize to get calculated frequencies ω_calc(k)
            │
            ▼
  [ Compute Residual: Σ |ω_calc(k) - ω_exp(k)|^2 ]
            │
            ▼
  [ Optimizer updates UFF parameters ]

```

### Pro-Tip for Optimization Speed:

If your forcefield is linear with respect to the stiffness parameters you are fitting (like the spring constants $k_R$ in bond stretching or $k_\theta$ in angle bending), the total Hessian is just a linear combination of individual parameter contributions:

$$K(0,L) = \sum_{p} \alpha_p K^{(p)}(0,L)$$

If this holds true, you can compute the geometric part of the Hessian $K^{(p)}(0,L)$ for each parameter *exactly once* at the very beginning. Then, inside your optimization loop, you don't even need to rerun the forcefield code! You just scale the precomputed matrices by your changing stiffness coefficients $\alpha_p$, do the Bloch sum, and diagonalize. This makes the fitting routine incredibly fast.

Do you plan to code this directly in Python (using something like `ASE` or `SciPy`), or are you wrapping an existing molecular dynamics tool to extract the forces?

---

# DeepSeek-V4-pro

You’re exactly on the right track: **one set of force‑constant calculations (one “supercell Hessian”) is enough**. Once you have the real‑space interatomic force constants, you apply a Bloch‑wave phase factor to build the dynamical matrix at **any** wavevector **k**. You do **not** need to re‑compute forces for each **k**-point.

Below I outline the complete workflow, from the finite‑displacement forces in a supercell to the phonon dispersion ω(**k**).

## 1. What you need: the real‑space force‑constant matrix Φ

For a crystal, the potential energy expanded to second order is

\[
U = \frac12 \sum_{R,R'} \sum_{i,j} \sum_{\alpha,\beta} 
\Phi_{i\alpha,\,j\beta}(R'-R) \; u_{i\alpha}(R) \, u_{j\beta}(R'),
\]

where  
- \(R\) labels unit cells,  
- \(i,j\) label atoms inside a unit cell,  
- \(\alpha,\beta\) are Cartesian directions,  
- \(u_{i\alpha}(R)\) is the displacement of atom \(i\) in cell \(R\).  

The quantity \(\Phi_{i\alpha,\,j\beta}(R)\) is the force constant between atom \(i\) in a reference cell (cell \(0\)) and atom \(j\) in cell \(R\).  
It obeys \(\Phi_{i\alpha,\,j\beta}(R) = \Phi_{j\beta,\,i\alpha}(-R)\) and translational invariance requires \(\sum_{j,R} \Phi_{i\alpha,\,j\beta}(R) = 0\).

The **dynamical matrix** at wavevector **k** is then

\[
D_{i\alpha,\,j\beta}(\mathbf{k}) = \frac{1}{\sqrt{M_i M_j}} 
\sum_{R} \Phi_{i\alpha,\,j\beta}(R) \; e^{i\mathbf{k}\cdot \mathbf{R}},
\]

and the phonon frequencies ω(**k**) are the square roots of the eigenvalues of \(D(\mathbf{k})\).

## 2. Building Φ from finite‑displacement forces in a supercell

### 2.1 Choose a supercell
Pick a supercell large enough that force constants decay to near zero within half its dimensions. For molecular crystals or UFF, a \(2\times2\times2\) or \(3\times3\times3\) supercell is often sufficient, but you should check convergence.

### 2.2 Displace atoms and record forces
Do **finite differences** (or analytical derivatives if your code provides them):

1. Construct the supercell with atoms at their equilibrium positions.
2. Loop over **symmetry‑inequivalent atoms in the central unit cell** and over the three Cartesian directions.  
   For each \((i_0, \alpha)\):
   - Displace atom \(i_0\) (in the reference cell, \(R=0\)) by a small amount \(+\Delta\) along \(\alpha\).
   - Compute the forces on **every atom in the supercell**.
   - Displace by \(-\Delta\), compute forces again.
   - The force‑constant matrix element between the displaced atom \((i_0,\alpha)\) in cell \(0\) and any atom \((j,\beta)\) in cell \(R\) is  
     \[
     \Phi_{j\beta,\,i_0\alpha}(R) \approx -\frac{F_{j\beta}^{+\Delta} - F_{j\beta}^{-\Delta}}{2\Delta}
     \]
     (The negative sign appears because the force is the negative gradient of energy: \(F_j = -\partial U/\partial u_j\).)
   - Use symmetry to obtain remaining elements (e.g., \(\Phi_{i_0\alpha,\,j\beta}(R) = \Phi_{j\beta,\,i_0\alpha}(-R)\)).

After the loop you have a set of matrices \(\Phi(R)\) for all \(R\) vectors that connect atoms within the supercell.

### 2.3 Assign the correct \(R\) vector
When you displace atom \(i_0\) in the **central** cell, an atom \(j\) in the supercell belongs to some unit cell located at a lattice vector \(R\) relative to the central cell.  
If your supercell is built by replicating the unit cell with integer multiples of the lattice vectors \(a_1, a_2, a_3\), then \(R = n_1 a_1 + n_2 a_2 + n_3 a_3\) with \(n_i \in [-L_i/2, L_i/2)\) (wrap as needed, using the minimum image convention for the supercell).  

Make sure that **all \(R\) vectors are taken with respect to the same origin** (the central cell). The force‑constant you get from the displacement of \(i_0\) in cell 0 and the force on \(j\) in cell \(R\) gives \(\Phi_{j\beta,\,i_0\alpha}(R)\). By translation invariance you can relabel:  
\[
\Phi_{i\alpha,\,j\beta}(R) \; \text{(force on } j \text{ in cell } R \text{ due to displacement of } i \text{ in cell 0)}.
\]

## 3. Fourier transform: D(k) at any k

Once you have \(\Phi_{i\alpha,\,j\beta}(R)\) for a set of \(R\) vectors that covers all significant interactions, you can compute the dynamical matrix at **any** **k**-vector (not just commensurate with the supercell) by a simple phase sum:

\[
D_{i\alpha,\,j\beta}(\mathbf{k}) = \frac{1}{\sqrt{M_i M_j}} 
\sum_{R} \Phi_{i\alpha,\,j\beta}(R) \; e^{i\mathbf{k}\cdot \mathbf{R}}.
\]

The sum runs over all \(R\) for which you have force constants (including \(R=0\)).  
Because you only have force constants for a finite set of \(R\) (those within the supercell), this is an **interpolation**. It works well if the force constants are short‑ranged enough that the truncation error is negligible.

### Algorithm sketch
```python
# Precompute: for each pair of unit-cell atoms (i,j) and for each
# R vector, store the 3x3 block Phi[i,j,R] (alpha,beta)
# This block is Phi_{i alpha, j beta}(R)

def dynamical_matrix(k):
    D = zeros((3*N_atoms_cell, 3*N_atoms_cell), complex)
    for i in range(N_atoms_cell):
        for j in range(N_atoms_cell):
            block = zeros((3,3), complex)
            for R, Phi_R in Phi_blocks[(i,j)]:
                phase = exp(1j * dot(k, R))
                block += Phi_R * phase
            D[3*i:3*i+3, 3*j:3*j+3] = block / sqrt(M[i]*M[j])
    return D
```

Then for each **k** of interest:
- Build \(D(\mathbf{k})\),
- Diagonalize (it is Hermitian),
- The eigenvalues \(\omega^2\) give the phonon frequencies.

## 4. Important details and sanity checks

* **Acoustic sum rule (ASR)** – Because of numerical errors, your raw \(\Phi(R)\) will not satisfy \(\sum_{j,R} \Phi_{i\alpha,\,j\beta}(R) = 0\) exactly.  
  A simple correction: for each \(i,\alpha,\beta\), subtract a constant from the \(R=0\) block to enforce it. A common procedure is to add a correction \(\Delta_{i\alpha,i\alpha} = -\sum_{j\neq i, R} \Phi_{i\alpha,\,j\beta}(R)\) and symmetrize.

* **Check \(\Gamma\)-point (k=0)** – At **k**=0, \(D\) should reproduce the vibrational modes of the crystal unit cell with **periodic** boundary conditions. Compare with a separate molecule‑style calculation of the supercell’s Hessian (the eigenvalues of the supercell dynamical matrix at Γ should contain exact acoustic modes with ω=0). You can also check that you get three zero‑frequency acoustic modes.

* **Symmetry** – After fitting to experimental phonon spectra, you may need to refine the force constants. You already have all the ingredients to compute the theoretical spectrum for any **k**, so you can define a loss function comparing theory and experiment, and use it to adjust UFF parameters.

## 5. Direct answer to your question

> Do I need to do it from start or do I only take one Hessian and then apply the Bloch k on that Hessian?

You **only take one set of force constants** (from a single supercell calculation) and **then apply the Bloch phase factors** to build \(D(\mathbf{k})\).  
There is no need to repeat the displacement calculations for each **k**. This is the standard *real‑space force‑constant* (or *supercell*) method, and it works because the force constants are short‑ranged. The “Bloch k” simply appears in the Fourier sum.

If you are already computing the supercell Hessian (a giant 3N×3N matrix) you could, in principle, diagonalize it only at **k**=0. For arbitrary **k** you would have to extract the blocks Φ(R) as described above anyway, because the supercell Hessian is a collection of all Φ(R) arranged in a periodic way. The procedure above does exactly that in a clean, physically transparent manner.

Let me know if you’d like a more detailed explanation of any step or help with handling the periodic images and R‑vector assignment!