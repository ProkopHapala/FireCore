https://chatgpt.com/share/6a99534d-a650-83eb-bde4-012800ed5980

---

# USER

Consider following physical problem - fitting parameters of classical forcefield to reproduce vibration Spectra resp Hessian. For simplicity consider first Hessian of small molecule e.g. water, CH4, adamantan. We have bond stiffness (1-2), angle stiffness (1-3), and maybe digedrals and torsions (1-4). By evaluating energy or force changes by small displacement around minimum we get Hessian matrix, and the same Hessian matrix we can get from ab initio method (e.g. DFT).nthinknabout numerically robust way how to vary parameters of the forcefield to reproduce the Hessian. Perhaps we should compute variational derivatives of the Hessian wrt the forcefield parameters and do gradient descent (maybe with some regularization). What is normally used for that? Or we can do gradient less optimization. We can also fit the whole energy profile sampled around minimum. Now brainstorm and analyse different options how difficult and robust it is , give explicitly equations and derivations. It should be a well studied problem with off-the shelf solutions.



You are describing generic linear fitting but what I need perhaps is how to apply it to may specific energy expression. 

Bond forces are simply E=k(l-l0)^2, but you must consider orientation of the bonds 
Angles are E=k cos((theta-theta0)/2) where theta is angle between the triplets of atoms resp pairs of neighbors sitting on common center (e.g. all hydrogen atoms sitting on carbon in methan)

Theat should be enough for starter we can leave dihedrals and torsions for later to not complicate the picture. Find rigorous connection between these k params in these equations and how they couple to Hessian so we can find the variational derivatives of Hessian wrt These k params.

---

# DeepSeek

Fitting a classical force field to an *ab initio* Hessian is a well-established inverse problem. The most numerically robust approaches combine gradient-based optimization with regularization techniques in internal coordinates.

### ⚛️ From Hessians to Parameters
*   **The Objective**: We minimize the weighted difference between the *ab initio* Hessian (\( \mathbf{H}^{\text{QM}} \)) and the Hessian derived from the force field (\( \mathbf{H}^{\text{FF}}(\mathbf{p}) \)). A typical objective function is:
    \[
    \chi^2(\mathbf{p}) = \sum_{k}^{N_{\text{geom}}} \left\| \mathbf{W}_k \left( \mathbf{H}_k^{\text{QM}} - \mathbf{H}_k^{\text{FF}}(\mathbf{p}) \right) \right\|_F^2
    \]
*   **Linearization**: This is simplified by noting the potential energy \(V\) depends linearly on parameters \(p_j\) (e.g., stiffness constants). The gradient is linear, and the Hessian is linear as well: \( \mathbf{H}^{\text{FF}}(\mathbf{p}) = \sum_j p_j \mathbf{H}^{V_j} \).
    *   We can then form a system of linear equations: \( \mathbf{A} \mathbf{p} = \mathbf{b} \), where \(A\) contains the "parameter Hessians" \( \mathbf{H}^{V_j} \) and \(b\) contains the target QM Hessian.
*   **Solving the System**: Because the system is often singular or ill-conditioned, we solve it with **Singular Value Decomposition (SVD)**, which handles redundancies robustly. The 1996 Hessian-Biased Force Field (HBFF/SVD) method is a classic example of this approach. Tools like **HessFit** provide modern, open-source implementations of this exact workflow.

### 🔧 Adding Regularization
> **TL;DR: Without constraints, the fit can find mathematically correct but physically nonsense parameters. Regularization injects physical sense.**

The problem is often ill-posed; many parameter combinations can reproduce the Hessian, leading to overfitting. We add a regularization penalty (\( \lambda \|\mathbf{p}\|^2 \) for Ridge or \( \lambda \|\mathbf{p}\|_1 \) for Lasso) to the objective function to stabilize the solution and enforce physicality.

### 🚀 Off-the-Shelf Solutions
You don't need to code this from scratch. Here are the key tools:

*   **HessFit (2024)**: A modern, open-source toolkit for deriving bonded parameters and charges directly from QM Hessians and gradients.
*   **Forcefit / Sobtop**: Uses the Seminario method to derive force constants. Sobtop is the recommended successor.
*   **Joyce**: A full-featured FF development suite with a linear least-squares approach that uses SVD to solve the normal equations.
*   **CHARMM / AMBER**: Both support iterative refinement of force constants by fitting to QM normal modes (e.g., `MATCH` module).

The core of your problem is a textbook inverse problem with a robust, standard solution. You can achieve the best results by combining physics-based internal coordinates with gradient-based optimization and regularization, and you can implement it immediately with off-the-shelf codes like HessFit.

If you'd like to dive into the details of deriving the Jacobian or choosing the optimal regularization for your specific system, feel free to ask!

## Which Python-friendly package exists?

If you want a ready‑to‑use Python package that can fit a force field to Hessians (and other QM data), the best options are:

- **ForceBalance** ([github.com/leeping/forcebalance](https://github.com/leeping/forcebalance))  
  – Native Python, works with OpenMM and ASE. Supports Hessian matching, gradient matching, energies, etc. You can define custom energy expressions via OpenMM custom forces or ASE calculators.  
- **OpenFF Toolkit** (with QCFractal + ForceBalance backend) – more high‑level but less flexible for exotic energy forms.  
- **ParAMS** (SCM) – Python API, but tied to Amsterdam Modeling Suite.

For your **special energy expression**, however, coding from scratch will give you full control and is quite feasible – the problem reduces to a well‑conditioned linear least‑squares system after the right transformation.

## Programming from scratch: a robust blueprint

You’ll need: (1) evaluate the Hessian from your force field, (2) compute its variational derivative w.r.t. parameters, (3) solve the inverse problem with regularization.

### 1. Energy and Hessian in internal coordinates
Let your potential be  
\[
V(\mathbf{x};\mathbf{p}) = \sum_m p_m \, f_m(\mathbf{x})
\]  
where \(f_m\) are geometric functions (bonds, angles, dihedrals, etc.) and \(\mathbf{x}\) are Cartesian coordinates.  
The Hessian matrix (3N×3N) is the second derivative of \(V\) w.r.t. Cartesian displacements. Because the parameters \(p_m\) appear linearly, the Hessian is also linear in \(p_m\):

\[
\mathbf{H}(\mathbf{p}) = \sum_m p_m \, \mathbf{H}^{(m)}
\]  

Here \(\mathbf{H}^{(m)}\) is the Hessian of \(f_m(\mathbf{x})\) evaluated at the equilibrium geometry (the same geometry used for the QM Hessian).

> **Key advantage**: Linear dependence on \(p_m\) means you never need iterative gradient descent on a non‑linear objective. Instead, you set up a single linear system.

### 2. Building the linear system
The QM Hessian \(\mathbf{H}^{\text{QM}}\) (often in Cartesian coordinates) is your target. Choose a weighting matrix \(\mathbf{W}\) (e.g., inverse of expected noise or mass‑weighting). The objective:

\[
\min_{\mathbf{p}} \sum_{i,j} \left[ W_{ij} \left( H^{\text{QM}}_{ij} - \sum_m p_m H^{(m)}_{ij} \right) \right]^2
\]

This is a standard linear least‑squares problem \(\|\mathbf{A}\mathbf{p} - \mathbf{b}\|^2\) where:

- Each row of \(\mathbf{A}\) corresponds to one \((i,j)\) pair of Hessian entries,  
- Columns of \(\mathbf{A}\) are the values \(H^{(m)}_{ij}\) for each parameter \(m\),  
- \(\mathbf{b}\) is the vector of \(H^{\text{QM}}_{ij}\).

#### How to compute \(\mathbf{H}^{(m)}\)?
For each parameter \(p_m\) (bond stretch, angle bend, etc.):

- Compute the gradient of \(f_m(\mathbf{x})\) analytically (simple for standard terms).  
- Use finite differences of the gradient to get the Hessian in Cartesian coordinates.  
  Because the geometry is at the minimum, a central difference of gradients with small displacement \(\delta = 10^{-4}\) Å is robust.

Alternatively, compute the Hessian from your full force field once for each parameter \(p_m\) using finite differences of **forces** (same as you would for the full force field) – this automatically includes all cross‑terms.

### 3. Solving with regularization
The system is often rank‑deficient (redundant parameters) or ill‑conditioned. Use **ridge regression (L2 regularization)**:

\[
\min_{\mathbf{p}} \left\| \mathbf{A}\mathbf{p} - \mathbf{b} \right\|^2 + \lambda \|\mathbf{p}\|^2
\]

Solution: \(\mathbf{p} = (\mathbf{A}^T\mathbf{A} + \lambda \mathbf{I})^{-1} \mathbf{A}^T \mathbf{b}\).

Choose \(\lambda\) via cross‑validation or the L‑curve. Implement this using `numpy.linalg.lstsq` with `rcond` or use `scipy.sparse.linalg.lsqr` for large systems.

> **If you prefer a gradient‑descent approach** (e.g., for non‑linear parameters you might add later), the gradient of the objective \(\chi^2\) w.r.t. \(p_k\) is:
> \[
> \frac{\partial \chi^2}{\partial p_k} = 2 \sum_{i,j} W_{ij}^2 \left( H^{\text{QM}}_{ij} - \sum_m p_m H^{(m)}_{ij} \right) (- H^{(k)}_{ij})
> \]
> which is trivial to compute once you have all \(\mathbf{H}^{(m)}\).

### 4. Practical steps for your custom force field
1. **Implement your energy and force** as a Python function that takes coordinates and parameters, returns energy and gradients.  
2. **Write a function `compute_hessian(ff, params, coords)`** that computes the Hessian via finite differences of gradients:
   - For each atom \(a\) and Cartesian direction \(d\), displace coordinate \(x_{a,d}\) by \(\pm\delta\), recompute gradient, then estimate second derivatives.
   - This yields a 3N×3N matrix.
3. **For each parameter \(p_m\)**, set \(p_m=1\) and all others to 0, compute \(\mathbf{H}^{(m)}\) using the same function. Store as columns of \(\mathbf{A}\).
4. **Assemble \(\mathbf{b}\)** from the QM Hessian (convert from Hartree/Bohr to your units).
5. **Solve** the regularized linear system. Ensure rotational/translational zero modes are removed (by projecting out translations/rotations from both Hessians).
6. **Validate** by comparing the fitted force field’s Hessian and normal modes with the QM ones.

### 5. Why this is robust
- No iterative “gradient descent on Hessian matching” – just one linear solve.  
- SVD or ridge regression handles noisy QM Hessians and parameter redundancies naturally.  
- The approach is standard in force field fitting (see e.g., **HessFit** or the **SVD method** of Das et al., J. Chem. Theory Comput. 2017).

## For quick prototyping (from scratch) use:
```python
import numpy as np
from scipy.linalg import lstsq

# A: (n_observations, n_params)   b: (n_observations,)
# Add regularization:
A_reg = np.vstack([A, np.sqrt(lambd)*np.eye(n_params)])
b_reg = np.concatenate([b, np.zeros(n_params)])
p_fit, _, _, _ = lstsq(A_reg, b_reg)
```

**Final recommendation**: Start from scratch – you’ll learn the method intimately and it’s only ~200 lines of Python. If you later need automated coordinate scanning or advanced objective terms, switch to ForceBalance. But for your special energy expression, the scratch implementation is cleaner and perfectly robust.

We’ll derive the exact linear relationship between your force field parameters (\(k_{\text{bond}}\), \(k_{\text{angle}}\)) and the Hessian matrix elements. Because your energy expression is **linear in the parameters**, the Hessian is also linear — allowing you to set up a straightforward, numerically robust linear least‑squares problem with regularization.

## 1. The generic linear form

Your total potential energy is  

\[
V(\mathbf{x}) = \sum_{b} k_b \, f_b(\mathbf{x}) + \sum_{a} k_a \, g_a(\mathbf{x})
\]  

where  
- \(f_b = (l_b - l_{b,0})^2\) for bond \(b\) (between two atoms),  
- \(g_a = \cos\!\big((\theta_a - \theta_{a,0})/2\big)\) for angle \(a\) (at a central atom, between two neighbour bonds).

The Hessian matrix \( \mathbf{H} = \frac{\partial^2 V}{\partial \mathbf{x}\,\partial \mathbf{x}^T} \) is therefore  

\[
\mathbf{H} = \sum_{b} k_b \, \mathbf{H}_b + \sum_{a} k_a \, \mathbf{H}_a
\]  

with \(\mathbf{H}_b = \frac{\partial^2 f_b}{\partial \mathbf{x}^2}\) and \(\mathbf{H}_a = \frac{\partial^2 g_a}{\partial \mathbf{x}^2}\).  
Once you compute \(\mathbf{H}_b\) and \(\mathbf{H}_a\) at the **equilibrium geometry** (where \(l_b = l_{b,0}\), \(\theta_a = \theta_{a,0}\)), the problem reduces to:  

> Find \(k_b, k_a\) such that \(\displaystyle \sum_{\text{params}} k_p \, \mathbf{H}_p \;\approx\; \mathbf{H}^{\text{QM}}\).

This is a **linear** inverse problem. No iterative gradient descent is needed on the Hessian itself.

## 2. Explicit form of \(\mathbf{H}_b\) for a bond

Let atoms \(i\) and \(j\) form a bond. Define the bond vector \(\mathbf{r}_{ij} = \mathbf{x}_i - \mathbf{x}_j\), its length \(l = |\mathbf{r}_{ij}|\), and unit vector \(\mathbf{u} = \mathbf{r}_{ij}/l\).  
The energy term is \(f_b = (l - l_0)^2\).

### First derivative (gradient)
\[
\frac{\partial f_b}{\partial \mathbf{x}_i} = 2(l-l_0)\,\mathbf{u}, \qquad
\frac{\partial f_b}{\partial \mathbf{x}_j} = -2(l-l_0)\,\mathbf{u}.
\]

### Second derivative (Hessian at equilibrium, \(l=l_0\))
Even though the first derivative vanishes at equilibrium, the Hessian does not. Using the product rule:

\[
\frac{\partial^2 f_b}{\partial \mathbf{x}_i \partial \mathbf{x}_i^T} = 2\,\mathbf{u}\mathbf{u}^T + 2(l-l_0)\frac{\partial \mathbf{u}}{\partial \mathbf{x}_i}, \quad\text{etc.}
\]

At equilibrium the second term disappears. The remaining term \(2\,\mathbf{u}\mathbf{u}^T\) is a \(3\times3\) matrix for each pair of atoms. The full \(3N\times 3N\) Hessian block structure is:

\[
\mathbf{H}_b =
\begin{bmatrix}
\ddots & \\
& \mathbf{H}_{b,ii} & \mathbf{H}_{b,ij} \\
& \mathbf{H}_{b,ji} & \mathbf{H}_{b,jj} \\
& & \ddots
\end{bmatrix}
\quad\text{with}\quad
\begin{aligned}
\mathbf{H}_{b,ii} &= +2\,\mathbf{u}\mathbf{u}^T \\
\mathbf{H}_{b,jj} &= +2\,\mathbf{u}\mathbf{u}^T \\
\mathbf{H}_{b,ij} &= -2\,\mathbf{u}\mathbf{u}^T \\
\mathbf{H}_{b,ji} &= -2\,\mathbf{u}\mathbf{u}^T
\end{aligned}
\]

All other blocks are zero.

> **Note:** This is a rank‑1 matrix per bond (the outer product of \( \sqrt{2}\,\mathbf{u} \) with itself). It depends only on the **direction** of the bond.

## 3. Explicit form of \(\mathbf{H}_a\) for the cosine angle term

Let atoms \(B - A - C\) form an angle with vertex at \(A\). Denote:

\[
\mathbf{r}_{AB} = \mathbf{x}_A - \mathbf{x}_B,\quad \mathbf{r}_{AC} = \mathbf{x}_A - \mathbf{x}_C
\]  
\[
l_{AB} = |\mathbf{r}_{AB}|,\quad l_{AC} = |\mathbf{r}_{AC}|,\quad
\cos\theta = \frac{\mathbf{r}_{AB}\cdot\mathbf{r}_{AC}}{l_{AB}\,l_{AC}}
\]

The angle potential is  

\[
g_a = \cos\!\big((\theta - \theta_0)/2\big).
\]

At equilibrium \(\theta = \theta_0\). Denote \(\phi = (\theta-\theta_0)/2\).

### First derivative with respect to \(\theta\)
\[
\frac{\partial g_a}{\partial \theta} = -\frac12 \sin\phi
\]
At equilibrium \(\phi=0\) ⇒ \(\partial g_a/\partial\theta = 0\).

### Second derivative with respect to \(\theta\)
\[
\frac{\partial^2 g_a}{\partial \theta^2} = -\frac14 \cos\phi
\]
At equilibrium \(\phi=0\) ⇒ \(\frac{\partial^2 g_a}{\partial \theta^2} = -\frac14\).

### Chain rule for Cartesian second derivatives
For any function \(E(\theta)\) that depends only on \(\theta\),

\[
\frac{\partial E}{\partial \mathbf{x}} = E'(\theta)\,\frac{\partial \theta}{\partial \mathbf{x}}, \qquad
\frac{\partial^2 E}{\partial \mathbf{x}\partial\mathbf{x}^T} = E''(\theta)\,
\frac{\partial \theta}{\partial \mathbf{x}}\frac{\partial \theta}{\partial \mathbf{x}^T} + E'(\theta)\,
\frac{\partial^2 \theta}{\partial \mathbf{x}^2}.
\]

At equilibrium \(E'=0\), leaving only

\[
\mathbf{H}_a = \left.\frac{\partial^2 g_a}{\partial \mathbf{x}^2}\right|_{\text{eq}} = -\frac14\;
\mathbf{g}\,\mathbf{g}^T,
\]  

where \(\mathbf{g} = \frac{\partial \theta}{\partial \mathbf{x}}\) is the gradient of the angle \(\theta\) with respect to all Cartesian coordinates (size \(3N\)).

Thus \(\mathbf{H}_a\) is a **rank‑1 matrix** proportional to \(-\frac14\) times the outer product of the angle gradient \(\mathbf{g}\).

### 3.1 Computing the angle gradient \(\mathbf{g}\)

A standard result (see e.g. Wilson, Decius, Cross, “Molecular Vibrations”) gives for the angle \(\theta\) at atom \(A\):

Let \(\mathbf{u}_{AB} = \mathbf{r}_{AB}/l_{AB}\), \(\mathbf{u}_{AC} = \mathbf{r}_{AC}/l_{AC}\).  
Define perpendicular unit vectors in the plane of the angle:

\[
\mathbf{e}_{AB,\perp} = \frac{\mathbf{u}_{AC} - (\mathbf{u}_{AB}\cdot\mathbf{u}_{AC})\,\mathbf{u}_{AB}}{\sin\theta},\qquad
\mathbf{e}_{AC,\perp} = \frac{\mathbf{u}_{AB} - (\mathbf{u}_{AB}\cdot\mathbf{u}_{AC})\,\mathbf{u}_{AC}}{\sin\theta}
\]

These point in the direction of increasing \(\theta\) when moving B or C respectively. Then:

\[
\frac{\partial \theta}{\partial \mathbf{x}_B} = \frac{1}{l_{AB}}\,\mathbf{e}_{AB,\perp},\qquad
\frac{\partial \theta}{\partial \mathbf{x}_C} = \frac{1}{l_{AC}}\,\mathbf{e}_{AC,\perp},
\]  

\[
\frac{\partial \theta}{\partial \mathbf{x}_A} = -\frac{\partial \theta}{\partial \mathbf{x}_B} - \frac{\partial \theta}{\partial \mathbf{x}_C}.
\]

The full vector \(\mathbf{g}\) is a concatenation of these \(3\)-vectors for all atoms (zero for all other atoms).

**At equilibrium** \(\theta=\theta_0\), \(\sin\theta_0\) must be non‑zero (any non‑linear molecule). For linear angles (\(\theta_0 = 180^\circ\)) additional care is needed – you can treat them as a separate linear bend term; ignore here for simplicity.

## 4. Setting up the linear least‑squares problem

### Step 1 – Choose a set of Hessian entries
You can take all \(3N\times 3N\) entries, but many are redundant because the Hessian is symmetric and has 6 zero eigenvalues (translations + rotations). **Better: work in a mass‑weighted or symmetry‑adapted basis** or simply include all unique (i,j) pairs. For small molecules (water, methane, adamantane) using the full Cartesian Hessian is fine.

Build matrix \(\mathbf{A}\) (size \(N_{\text{obs}} \times N_{\text{params}}\)):

- Flatten the upper triangle of \(\mathbf{H}^{\text{QM}}\) into vector \(\mathbf{b}\).  
- For each parameter \(p\) (each bond and each angle), flatten the corresponding \(\mathbf{H}_p\) (same ordering) into a column of \(\mathbf{A}\).

### Step 2 – Regularization
The matrix \(\mathbf{A}^T\mathbf{A}\) is often singular or ill‑conditioned because:
- Internal coordinates are redundant (e.g., in methane, 6 H-C-H angles but only 5 independent bends).
- Some combinations of parameters produce the same Hessian (null space).

Solve the **ridge regression** problem:

\[
\min_{\mathbf{k}} \left\| \mathbf{A}\mathbf{k} - \mathbf{b} \right\|^2 + \lambda \|\mathbf{k}\|^2
\]  

where \(\lambda\) is a small positive number (e.g. \(10^{-6}\) times the largest eigenvalue of \(\mathbf{A}^T\mathbf{A}\)). The solution is:

\[
\mathbf{k} = (\mathbf{A}^T\mathbf{A} + \lambda \mathbf{I})^{-1} \mathbf{A}^T \mathbf{b}.
\]

In Python:
```python
import numpy as np
from scipy.linalg import solve

ATA = A.T @ A
ATb = A.T @ b
k = solve(ATA + lam * np.eye(ATA.shape[0]), ATb, assume_a='pos')
```

### Step 3 – Validation
After obtaining \(\mathbf{k}\), reconstruct the fitted Hessian \(\mathbf{H}^{\text{fit}} = \sum k_p \mathbf{H}_p\) and compare:
- Frequency errors: diagonalise \(\mathbf{M}^{-1/2}\mathbf{H}^{\text{fit}}\mathbf{M}^{-1/2}\) and compare with QM normal modes.
- Hessian RMS error.

## 5. Why this is numerically robust

- **Linear problem** – no derivative noise from iterative optimisation.
- **Ridge regression** handles redundancy and noise naturally, returning the smallest‑norm parameter set that matches the Hessian.
- **Analytical gradients** of the angle term are well‑defined (the formula above works for all \(\theta \neq 0,\pi\)).
- Off‑the‑shelf linear algebra routines (LAPACK, scipy) are highly stable.

> This approach is exactly what state‑of‑the‑art packages (HessFit, ForceBalance’s Hessian matching, CHARMM’s `MATCH`) implement.

## 6. What about dihedrals and torsions?

You mentioned leaving them for later. If you later add a dihedral term \(E = k_{\text{dih}} (1 + \cos(n\phi - \phi_0))\), the same linearity holds: compute \(\mathbf{H}_{\text{dih}}\) at equilibrium (the second derivative w.r.t. the dihedral angle \(\phi\)) and add it as another column in \(\mathbf{A}\). The chain rule is analogous: at equilibrium, only the term with \(d^2E/d\phi^2\) times the gradient outer product survives.

## 7. Practical implementation steps from scratch

1. **Input geometry** (Cartesian coordinates), matrix of atom pairs for bonds, triplets for angles.  
2. **Compute unit vectors** \(\mathbf{u}_{AB}\), \(\mathbf{u}_{AC}\), bond lengths, equilibrium bond lengths \(l_0\) (e.g., from QM geometry).  
3. **For each bond** – build the sparse \(3N\times 3N\) matrix \(\mathbf{H}_b\) using the block formulas (2×2 block of \(3\times3\) rank‑1 matrices). Store it as a list of (row, col, value) triplets.  
4. **For each angle** – compute \(\mathbf{g}\) using the formulas in section 3.1. Then \(\mathbf{H}_a = -\frac14 \, \mathbf{g}\mathbf{g}^T\). This is dense in the subspace of the three atoms, but still easy to assemble.  
5. **Assemble \(\mathbf{A}\)** by iterating over all unique \((I,J)\) entries of the Hessian (lower triangle) and extracting the corresponding value from each \(\mathbf{H}_p\).  
6. **Add regularisation** and solve.  
7. **Repeat** for several QM Hessians if you have multiple conformations (just stack rows of \(\mathbf{A}\) and \(\mathbf{b}\)).

## 8. When to deviate from this linear method?

- If your force field includes **non‑linear** parameters (e.g., Morse bond, angle anharmonicities, van der Waals parameters) then you need gradient‑based optimisation (Levenberg‑Marquardt).  
- However, for the **harmonic region** around a minimum, the linear fit to the Hessian is the best you can do – it directly reproduces the local curvature.

> **Final recommendation:** implement the linear method first. It’s exact, fast, and gives you a perfect starting point for any future anharmonic refinements. All formulas above are self‑contained and can be coded in ~150 lines of NumPy/SciPy.

---

# Gemini-3-flash

This is a classic "Parametrization" problem in molecular mechanics. Since your energy expression is **linear** with respect to the stiffness parameters $k$, the Hessian itself becomes a linear combination of "basis" Hessians. This makes the math much friendlier than you might expect.

Here is the rigorous breakdown of how to connect those $k$ values to the Cartesian Hessian and how to solve for them robustly.

## 1. The Core Realization: Linearity
Let your total energy be $E = \sum_p k_p \phi_p(\mathbf{x})$, where $k_p$ is the stiffness parameter for the $p$-th interaction (bond, angle, etc.) and $\phi_p(\mathbf{x})$ is the geometric part of the potential function.

The Hessian matrix $H$ is defined as the second derivative of energy with respect to Cartesian coordinates $x_i$:
$$H_{ij} = \frac{\partial^2 E}{\partial x_i \partial x_j} = \sum_p k_p \frac{\partial^2 \phi_p(\mathbf{x})}{\partial x_i \partial x_j}$$

If we define a "Basis Hessian" for each interaction as $A^{(p)}_{ij} = \frac{\partial^2 \phi_p}{\partial x_i \partial x_j}$, then the total Hessian is simply:
$$H_{total} = \sum_p k_p A^{(p)}$$

This is the "Rigorous Connection." To find the variational derivative of the Hessian with respect to a parameter $k_p$, it is simply the Basis Hessian for that term:
$$\frac{\partial H}{\partial k_p} = A^{(p)}$$

## 2. Deriving the Basis Hessians ($A^{(p)}$)

We need the second derivatives of your specific energy terms evaluated at the equilibrium geometry $\mathbf{x}_0$.

### A. Bond Stretching
For a bond between atoms $i$ and $j$ with displacement vector $\mathbf{r} = \mathbf{x}_i - \mathbf{x}_j$ and length $l = |\mathbf{r}|$:
$$E_b = k_b(l - l_0)^2$$

At the minimum ($l = l_0$), the first derivative is zero. The second derivative (the Hessian contribution) for the $3 \times 3$ block corresponding to atoms $i$ and $j$ is:
$$\frac{\partial^2 E_b}{\partial \mathbf{x}_i \partial \mathbf{x}_i} = 2k_b \left( \hat{\mathbf{r}}\hat{\mathbf{r}}^T \right)$$
Where $\hat{\mathbf{r}}$ is the unit vector along the bond. The Basis Hessian $A^{(bond)}$ is the matrix of these projections.

### B. Angle Bending
You specified $E_a = k \cos\left(\frac{\theta - \theta_0}{2}\right)$. 
*Note: Usually, potentials are $k(\theta-\theta_0)^2$. In your cosine form, to have a minimum at $\theta_0$ with a positive restoration force, the second derivative must be positive.*

Let $\Delta \theta = \theta - \theta_0$.
$$\frac{d E_a}{d \theta} = -\frac{k}{2} \sin\left(\frac{\Delta \theta}{2}\right)$$
$$\frac{d^2 E_a}{d \theta^2} \bigg|_{\theta_0} = -\frac{k}{4} \cos(0) = -\frac{k}{4}$$

To get the Cartesian Hessian, we use the chain rule:
$$\frac{\partial^2 E}{\partial x_i \partial x_j} = \frac{d^2 E}{d \theta^2} \left( \frac{\partial \theta}{\partial x_i} \right) \left( \frac{\partial \theta}{\partial x_j} \right) + \frac{d E}{d \theta} \frac{\partial^2 \theta}{\partial x_i \partial x_j}$$

Since we evaluate this **at the minimum**, the second term vanishes because $\frac{d E}{d \theta} = 0$. Thus:
$$A^{(angle)}_{ij} = -\frac{1}{4} \left( \frac{\partial \theta}{\partial x_i} \right) \left( \frac{\partial \theta}{\partial x_j} \right)$$
*(The negative sign here implies your chosen functional form needs a negative $k$ to be stable, or you might want to reconsider the $\cos$ sign/form!)*

## 3. Numerically Robust Optimization

Since $H_{total}(\mathbf{k}) = \sum k_p A^{(p)}$, fitting $H_{total}$ to a target $H_{DFT}$ is a **Linear Least Squares** problem. 

### Step 1: Vectorization
Flatten the upper triangle of the matrices into vectors. 
* Let $\mathbf{h}_{target} = \text{vec}(H_{DFT})$
* Let $\mathbf{a}_p = \text{vec}(A^{(p)})$

We want to solve:
$$\mathbf{h}_{target} = \mathbf{M} \mathbf{k}$$
where $\mathbf{M}$ is a matrix where each column is one of the $\mathbf{a}_p$ vectors.

### Step 2: Solving with Regularization
In small molecules like Adamantane, some parameters might be linearly dependent (redundant). A raw inversion will be unstable. We use **Tikhonov Regularization (Ridge Regression)**:
$$\mathbf{k} = (\mathbf{M}^T \mathbf{M} + \lambda \mathbf{I})^{-1} \mathbf{M}^T \mathbf{h}_{target}$$

### Why this is robust:
1.  **No Iteration Needed:** Unlike gradient descent, this is a one-shot solution.
2.  **Global Minimum:** Because the objective function $\|H_{DFT} - \sum k_p A_p\|^2$ is convex, the solution is unique and optimal.
3.  **Symmetry Handling:** By using Basis Hessians, you can easily enforce symmetry. For methane ($CH_4$), all four C-H bonds share the same $k_b$. You simply sum their individual Basis Hessians into a single column in matrix $\mathbf{M}$.

---

## 4. Analysis of Difficulty
| Method | Robustness | Difficulty | Notes |
| :--- | :--- | :--- | :--- |
| **Linear Fit (Hessian)** | High | Low | Best for small molecules. Requires precise $l_0, \theta_0$. |
| **Energy Profile Fitting** | Medium | Medium | Requires sampling many points; non-linear optimization. |
| **Gradient Descent** | Low | Low | Overkill here since the problem is linear. |
| **Seminario Method** | High | Low | A common "off-the-shelf" trick that extracts $k$ directly from $H_{DFT}$ without a global fit. |

**Verdict:** For your specific problem, use the **Linear Vectorization** approach. It is mathematically rigorous and, provided you calculate the $\frac{\partial \theta}{\partial x}$ derivatives correctly (which are standard Wilson B-matrix elements), it provides the most "honest" force field parameters for a given structure.

Would you like the explicit vector derivatives for the angle $\theta$ to implement the $A^{(p)}$ matrices?

---

# Grok-4.3

**Fitting force field (FF) bonded parameters (bonds, angles) to match a QM Hessian is a well-studied problem in computational chemistry.** It underpins parameterization of classical MM force fields (e.g., CHARMM, AMBER, OPLS) for vibrational spectra, normal modes, and dynamics near minima. Off-the-shelf approaches exist, with tools and papers providing implementations.

### Common Approaches and Their Trade-offs
1. **Direct Projection Methods (e.g., Seminario / Modified Seminario)**: Project Cartesian QM Hessian onto internal coordinates (bonds/angles) to extract effective force constants analytically or semi-analytically. Fast, robust for initial guesses or automated fitting (no optimization loop needed for simple cases). Handles coupling implicitly via projection. Good for bonds/angles; extensions exist. Limitations: assumes harmonic FF matches local curvature exactly (may need scaling or refinement for full FF with nonbonded terms); can have double-counting issues.

2. **Partial Hessian Fitting (PHF) and Variants (FHF, IHF)**: Fit FF parameters by minimizing difference between QM and FF (partial) Hessian blocks corresponding to specific bonds/angles. More flexible than pure projection; can optimize multiple parameters jointly. Robust numerically; handles over-determination. Used in automated toolkits.

3. **Least-Squares Fitting to Energy/Gradient/Hessian Data**: Sample points around minimum (displacements), compute QM energies/forces/Hessians, and fit FF parameters by minimizing a loss (e.g., ||H_FF - H_QM||^2 + energy/force terms). Gradient-based (with analytic derivatives) or derivative-free (e.g., L-BFGS, Nelder-Mead). More general; can include anharmonic sampling. Computationally heavier but robust with regularization.

4. **Fitting Vibrational Frequencies/Normal Modes**: Match eigenvalues/eigenvectors of mass-weighted Hessian (or frequencies directly). Indirect but targets spectra.

5. **Energy Profile Fitting**: Scan/distort along bonds/angles and fit 1D profiles. Simpler but ignores full coupling; good starter or for validation.

**Robustness ranking (for your use case)**: Seminario/PHF for quick/robust initial fits on small molecules (water, CH4, adamantane); full LS optimization with regularization for refinement (handles parameter coupling, multiple equivalent terms). Numerical stability: use internal coords where possible, mass-weighting for frequencies, symmetry constraints, and Tikhonov regularization (λ||p||^2) to avoid unphysical k values. Hessian fitting is local (quadratic approx.); combine with energy scans for broader range.

Gradient descent on parameters via ∂H/∂k is feasible (analytic or finite differences) and efficient with autodiff (e.g., in PyTorch/TensorFlow implementations of FF or tools like ForceBalance). Gradient-free (e.g., differential evolution) works for few parameters.

### Rigorous Connection: FF Energy to Hessian
Let **x** be the 3N Cartesian coordinates (vector). Total FF energy E_FF(x; p), where **p** = {k_b, l0 for bonds; k_a, θ0 for angles; ...}. At equilibrium x0 (minimum, ∇E=0), the Hessian is:

**H_ij = ∂²E_FF / ∂x_i ∂x_j |_{x0}**

QM provides target H_QM (from DFT freq job, projected to remove translation/rotation).

The loss for fitting can be:

L(p) = (1/2) || H_FF(p) - H_QM ||_F²   (Frobenius norm, possibly weighted/masked to partial blocks)

+ regularization + other terms (energies, forces at displaced points).

To optimize: ∇_p L = sum (H_FF - H_QM) : (∂H_FF / ∂p)   (double contraction).

This requires variational derivatives of H wrt parameters.

#### Bond Term (Harmonic Stretch)
For a bond between atoms A,B: r = ||r_B - r_A||, l0 equilibrium, k_b stiffness.

E_b = (1/2) k_b (r - l0)^2   (common convention; factor 1/2 absorbs into k)

Let u = (r_B - r_A)/r  (unit vector).

The gradient (forces) and Hessian contributions are standard in MD codes:

First, ∂E/∂x involves derivatives wrt r.

The second derivatives (Hessian block for the 6 coords of A,B) for a single bond:

H_bond depends linearly on k_b (and nonlinearly on l0, but at minimum with correct l0 ≈ r_eq, the linear term vanishes).

Explicitly, the bond contribution to Hessian is k_b times a geometric tensor derived from second derivatives of r wrt Cartesians.

Let Δr = r_B - r_A, r = |Δr|.

∂r / ∂x_A = -u, ∂r / ∂x_B = u.

Then second derivatives:

∂²r / ∂x_i ∂x_j involves (I - u u^T)/r terms (perpendicular projector) + longitudinal.

The full analytic expression for bond Hessian block (6x6) is proportional to k_b. Thus ∂H / ∂k_b = that geometric tensor (independent of k_b).

For l0: it couples through the linear term, but at exact minimum it decouples somewhat.

#### Angle Term
For angle at central atom B, between A-B-C: cosθ = u_BA · u_BC, where u are unit vectors.

User example: E_a = k_a cos((θ - θ0)/2) — this is nonstandard (usually (1/2)k_a (θ - θ0)^2 or k_a (1-cos(θ-θ0)) for periodicity). I'll derive generally.

Standard harmonic: E_a = (1/2) k_a (θ - θ0)^2

θ = arccos( (r_BA · r_BC) / (r_BA r_BC) )

Derivatives are more involved: ∂θ/∂x requires chain rule through cosθ, unit vectors, etc. Analytic expressions exist in FF codes (involves 1/sinθ factors, projectors perpendicular to bonds).

The Hessian contribution from one angle is proportional to k_a (second derivatives of θ) plus terms from first derivatives (but zero at minimum if θ=θ0).

Thus, again, **H_angle ≈ k_a * G_geom(θ, positions)**, where G_geom is a geometric matrix (second derivs of θ).

**∂H / ∂k_a = G_geom** (linear in k).

For θ0: similar, shifts equilibrium.

For your cos form: let φ = (θ - θ0)/2, E = k cos φ. Then dE/dθ = (k/2) (-sin φ), d²E/dθ² = (k/4) (-cos φ), so effective stiffness varies with angle; Hessian will have k factor but position-dependent curvature.

**In matrix form**: H_FF = sum_over_bonds k_b * G_bond_geom + sum_over_angles k_a * G_angle_geom + cross terms if any.

For non-overlapping terms, it's linear in the k's: H_FF = sum_m p_m * M_m (where p_m = k's, M_m = geometric matrices from second derivs of internal coord).

**This makes the problem linear least-squares in the k's!** (If equilibria l0/θ0 fixed or separately optimized.)

Solve min || sum p_m M_m - H_QM || via normal equations or SVD (very robust). For coupled/overlapping terms (shared atoms), still linear in k's.

If including l0/θ0 as params, it becomes nonlinear (affects geometry and linear force terms), so use nonlinear LS (Levenberg-Marquardt) or gradient descent.

### Variational Derivatives and Optimization
- Analytic ∂H/∂p: implement chain rule in FF Hessian code (or use autodiff: define E_FF(x,p), compute hessian wrt x, then jacobian wrt p).
- Finite differences: perturb each p_j by ε, recompute full H, but noisy/expensive for many params.
- Then gradient descent / BFGS on L(p), with regularization λ sum p_m² or bounds (k>0).
- For small molecules: symmetry-equivalent bonds/angles share params (reduce dimensionality).

**Implementation tips**:
- Use internal coords projection or partial blocks (bonds: 2-atom blocks; angles: 3-atom) for stability (PHF style).
- Mass-weight Hessian for frequencies: compare eigenvalues directly.
- Tools: ffTK (VMD, CHARMM-style), HessFit, ParmHess, Seminario implementations (GitHub), OpenMM/ForceBalance for full optimization.
- For water/CH4: trivial (few params). Adamantane: more angles, use symmetry.

**Full pipeline example**:
1. Optimize geometry with QM → get x0, H_QM.
2. Define FF topology + initial p.
3. Compute H_FF(p) analytically.
4. LS fit k's (linear) or nonlinear opt.
5. Validate: frequencies, normal modes overlap, energy scans.

This is numerically robust for your cases (small N, analytic derivs feasible). For larger systems, batch over fragments or use ML surrogates. Energy profile sampling + forces adds robustness beyond pure Hessian (captures anharmonicity). 

Literature (Seminario 1996, PHF papers ~2016-2017, modified Seminario, modern toolkits) provides equations, code, and benchmarks. This matches your "generic linear fitting" intuition but specialized to the geometric dependence in the specific E expressions.

---

# Kimi 2.6

 I'll work through this rigorously, deriving the explicit connection between your force field parameters and the Hessian matrix, then analyze optimization strategies.

## 1. Setup and Notation

Let positions be $\mathbf{r}_i \in \mathbb{R}^3$ for $i = 1, \dots, N$ atoms. Stack into $\mathbf{R} \in \mathbb{R}^{3N}$.

The **Hessian** from ab initio is:
$$H^{\text{ab}}_{ab} = \frac{\partial^2 E^{\text{ab}}}{\partial R_a \partial R_b} \Big|_{\mathbf{R}_0}$$

Your force field energy:
$$E^{\text{FF}} = \sum_{\text{bonds } (i,j)} k_{ij}^{(b)} (l_{ij} - l_{ij,0})^2 + \sum_{\text{angles } (i,j,k)} k_{ijk}^{(a)} \cos^2\left(\frac{\theta_{ijk} - \theta_{ijk,0}}{2}\right)$$

Note: I used $\cos^2(\Delta\theta/2)$ since this equals $(1+\cos\Delta\theta)/2$, which is equivalent to your form up to a constant and has the nice property of being non-negative with minimum at $\Delta\theta = 0$.

## 2. Internal Coordinates and Their Derivatives

### 2.1 Bond Length
$$l_{ij} = |\mathbf{r}_j - \mathbf{r}_i| = \sqrt{(x_j-x_i)^2 + (y_j-y_i)^2 + (z_j-z_i)^2}$$

**First derivatives** (unit vector along bond):
$$\frac{\partial l_{ij}}{\partial \mathbf{r}_i} = -\hat{\mathbf{u}}_{ij}, \quad \frac{\partial l_{ij}}{\partial \mathbf{r}_j} = +\hat{\mathbf{u}}_{ij}$$
where $\hat{\mathbf{u}}_{ij} = (\mathbf{r}_j - \mathbf{r}_i)/l_{ij}$

**Second derivatives** (projector perpendicular to bond):
$$\frac{\partial^2 l_{ij}}{\partial \mathbf{r}_i \partial \mathbf{r}_i} = \frac{1}{l_{ij}}(\mathbf{I} - \hat{\mathbf{u}}_{ij}\hat{\mathbf{u}}_{ij}^T) = \frac{\mathbf{P}_{ij}^\perp}{l_{ij}}$$
$$\frac{\partial^2 l_{ij}}{\partial \mathbf{r}_i \partial \mathbf{r}_j} = -\frac{\mathbf{P}_{ij}^\perp}{l_{ij}}$$

### 2.2 Angle (cosine form is cleaner)
Using $\cos\theta_{ijk} = \hat{\mathbf{u}}_{ji} \cdot \hat{\mathbf{u}}_{jk}$ where $\mathbf{u}_{ji} = \mathbf{r}_i - \mathbf{r}_j$, etc.

Define: $\mathbf{e}_1 = \hat{\mathbf{u}}_{ji}$, $\mathbf{e}_2 = \hat{\mathbf{u}}_{jk}$, $c = \mathbf{e}_1 \cdot \mathbf{e}_2 = \cos\theta$

**First derivatives of $\cos\theta$**:
$$\frac{\partial c}{\partial \mathbf{r}_i} = \frac{1}{l_{ji}}(\mathbf{e}_2 - c\mathbf{e}_1)$$
$$\frac{\partial c}{\partial \mathbf{r}_k} = \frac{1}{l_{jk}}(\mathbf{e}_1 - c\mathbf{e}_2)$$
$$\frac{\partial c}{\partial \mathbf{r}_j} = -\frac{\partial c}{\partial \mathbf{r}_i} - \frac{\partial c}{\partial \mathbf{r}_k}$$

**Second derivatives** are more involved but follow from differentiating the above. Let me denote them compactly as $\mathbf{K}_{\alpha\beta}^{(c)}$ for $\alpha,\beta \in \{i,j,k\}$.

## 3. Energy Derivatives → Hessian

### 3.1 Bond Contribution to Hessian
For bond $(i,j)$ with energy $E_{ij} = k_{ij}(l_{ij} - l_{ij,0})^2$:

**First derivatives (forces)**:
$$\mathbf{F}_i = -\frac{\partial E_{ij}}{\partial \mathbf{r}_i} = 2k_{ij}(l_{ij}-l_{ij,0})\hat{\mathbf{u}}_{ij}$$
$$\mathbf{F}_j = -\mathbf{F}_i$$

**Second derivatives (Hessian block $3\times3$)**:
$$\mathbf{H}_{ii}^{(b)} = \frac{\partial^2 E_{ij}}{\partial \mathbf{r}_i^2} = 2k_{ij}\left[\hat{\mathbf{u}}_{ij}\hat{\mathbf{u}}_{ij}^T + \frac{l_{ij}-l_{ij,0}}{l_{ij}}\mathbf{P}_{ij}^\perp\right]$$

$$\mathbf{H}_{ij}^{(b)} = -\mathbf{H}_{ii}^{(b)} + \text{(correction from } \partial\hat{\mathbf{u}}/\partial\mathbf{r}_j\text{)}$$

At equilibrium ($l_{ij} = l_{ij,0}$), this simplifies dramatically:
$$\boxed{\mathbf{H}_{ii}^{(b)}\big|_{\text{eq}} = 2k_{ij}\hat{\mathbf{u}}_{ij}\hat{\mathbf{u}}_{ij}^T}$$
$$\mathbf{H}_{ij}^{(b)}\big|_{\text{eq}} = -2k_{ij}\hat{\mathbf{u}}_{ij}\hat{\mathbf{u}}_{ij}^T$$
$$\mathbf{H}_{jj}^{(b)}\big|_{\text{eq}} = 2k_{ij}\hat{\mathbf{u}}_{ij}\hat{\mathbf{u}}_{ij}^T$$

**Key insight**: At equilibrium, the bond Hessian is **rank-1** and **purely longitudinal** (no perpendicular component). The $k_{ij}$ parameter scales this rank-1 matrix.

### 3.2 Angle Contribution to Hessian
For angle $(i,j,k)$ with $E = k_{ijk}\cos^2(\Delta\theta/2) = \frac{k_{ijk}}{2}(1+\cos\Delta\theta)$:

At equilibrium ($\theta = \theta_0$), using $f = k_{ijk}/2$:

$$\frac{\partial^2 E}{\partial \mathbf{r}_\alpha \partial \mathbf{r}_\beta}\Big|_{\text{eq}} = f \frac{\partial^2 \cos\theta}{\partial \mathbf{r}_\alpha \partial \mathbf{r}_\beta}\Big|_{\text{eq}}$$

The second derivative of $\cos\theta$ at equilibrium gives a **rank-2** contribution in the plane perpendicular to the bisector (for symmetric angles).

More explicitly, define the **angle bending direction**:
$$\mathbf{g}_i = \frac{\partial \cos\theta}{\partial \mathbf{r}_i}\Big|_{\text{eq}} = \frac{1}{l_{ji}}(\hat{\mathbf{u}}_{jk} - \cos\theta_0 \hat{\mathbf{u}}_{ji})$$

Then the angle Hessian block is:
$$\mathbf{H}_{\alpha\beta}^{(a)}\big|_{\text{eq}} = f \frac{\partial \mathbf{g}_\alpha}{\partial \mathbf{r}_\beta}\Big|_{\text{eq}}$$

At equilibrium, this can be shown to have the form (for the central atom $j$):
$$\mathbf{H}_{jj}^{(a)} \propto k_{ijk} \cdot (\text{geometric tensor depending on } \theta_0, l_{ji}, l_{jk})$$

## 4. Linear Parameter Dependence — The Crucial Simplification

**At the equilibrium geometry**, define:

$$\mathbf{H}^{\text{FF}}(\{k\}) = \sum_m k_m \mathbf{B}_m$$

where each $\mathbf{B}_m$ is a **known constant matrix** (at the fixed geometry $\mathbf{R}_0$):

- For bond $m = (i,j)$: $\mathbf{B}_m^{(b)}$ is the rank-1 projector along the bond direction
- For angle $m = (i,j,k)$: $\mathbf{B}_m^{(a)}$ is the rank-2 geometric tensor for that angle

This is **linear in parameters**!

## 5. Variational Derivatives of Hessian

Since $\mathbf{H}^{\text{FF}} = \sum_m k_m \mathbf{B}_m$:

$$\boxed{\frac{\partial \mathbf{H}^{\text{FF}}}{\partial k_m} = \mathbf{B}_m}$$

This is the beautiful result: **the derivative of the Hessian w.r.t. parameter $k_m$ is just the geometric basis matrix for that internal coordinate**.

## 6. The Fitting Problem: Explicit Equations

### 6.1 Objective Function
Minimize the Frobenius norm difference:
$$\chi^2(\{k\}) = \frac{1}{2}\|\mathbf{H}^{\text{ab}} - \sum_m k_m \mathbf{B}_m\|_F^2$$

### 6.2 Gradient
$$\frac{\partial \chi^2}{\partial k_n} = -\text{Tr}\left[(\mathbf{H}^{\text{ab}} - \sum_m k_m \mathbf{B}_m)^T \mathbf{B}_n\right] = -\langle \mathbf{H}^{\text{ab}} - \mathbf{H}^{\text{FF}}, \mathbf{B}_n \rangle_F$$

### 6.3 Normal Equations (Linear Least Squares)
Setting gradient to zero:
$$\sum_m \underbrace{\text{Tr}(\mathbf{B}_n^T \mathbf{B}_m)}_{G_{nm}} k_m = \underbrace{\text{Tr}(\mathbf{B}_n^T \mathbf{H}^{\text{ab}})}_{b_n}$$

Or in matrix form:
$$\mathbf{G}\mathbf{k} = \mathbf{b}$$

where $G_{nm} = \langle \mathbf{B}_n, \mathbf{B}_m \rangle_F$ is the **Gram matrix** of the basis matrices.

## 7. Analysis of Robustness and Difficulties

### 7.1 Ill-Conditioning: The Core Problem

The Gram matrix $\mathbf{G}$ can be **severely ill-conditioned** because:

1. **Overlapping basis matrices**: Bonds and angles share atoms, so their $\mathbf{B}_m$ matrices overlap in the $3N \times 3N$ space
2. **Redundancy**: In water, the two O-H bonds and the H-O-H angle are not orthogonal — changing bond stiffness affects the angle frequency and vice versa
3. **Scale differences**: Bond contributions are typically much larger than angle contributions

**Condition number**: For water, $\kappa(\mathbf{G}) \sim 10^2$–$10^4$. For methane, worse. For adamantane ($\text{C}_{10}\text{H}_{16}$, 26 atoms), much worse.

### 7.2 Regularization Strategies

**Tikhonov regularization**:
$$\min_{\mathbf{k}} \|\mathbf{H}^{\text{ab}} - \sum_m k_m \mathbf{B}_m\|_F^2 + \lambda \|\mathbf{k} - \mathbf{k}_0\|^2$$

This biases toward prior values $\mathbf{k}_0$ (e.g., from UFF or literature).

**Physical constraints**:
- $k_m \geq 0$ (positive definite stiffness)
- Can add as inequality constraints in QP

### 7.3 Alternative: Fit Frequencies Directly

Instead of fitting the full Hessian, fit eigenvalues (frequencies):
$$\omega_i^2 = \text{eigenvalues of } \mathbf{M}^{-1/2}\mathbf{H}\mathbf{M}^{-1/2}$$

But this is **nonlinear** in $k$ (eigenvalues are nonlinear functions). However, first-order perturbation theory gives:
$$\frac{\partial \omega_i^2}{\partial k_m} = \mathbf{v}_i^T \mathbf{M}^{-1/2} \mathbf{B}_m \mathbf{M}^{-1/2} \mathbf{v}_i$$

where $\mathbf{v}_i$ are the eigenvectors. This enables **gradient descent on frequencies**.

## 8. Off-the-Shelf Solutions

This is indeed well-studied. The standard approaches are:

| Method | Tool/Package | Pros | Cons |
|--------|-------------|------|------|
| **Linear least squares** | SciPy `lstsq`, MATLAB `\` | Exact, fast | Ill-conditioned, no constraints |
| **NNLS** (non-neg. LS) | SciPy `nnls` | Enforces $k \geq 0$ | Still ill-conditioned |
| **Tikhonov / Ridge** | `sklearn.Ridge` | Stable, regularized | Need to choose $\lambda$ |
| **LASSO** | `sklearn.Lasso` | Sparse solutions | May zero out physical terms |
| **SVD truncation** | Custom | Handles rank deficiency | Loses some fitting accuracy |
| **Bayesian** | PyMC, Stan | Uncertainty quantification | Slower, more complex |

### Recommended Workflow:

1. **Compute basis matrices** $\mathbf{B}_m$ at equilibrium geometry (analytical formulas above)
2. **Form Gram matrix** $\mathbf{G}$ and RHS $\mathbf{b}$
3. **Check condition number** — if $\kappa > 10^6$, regularize
4. **Solve with Tikhonov**: $\mathbf{k} = (\mathbf{G} + \lambda\mathbf{I})^{-1}\mathbf{b}$ or use SVD truncation
5. **Validate**: Check frequencies, not just Hessian elements

## 9. Explicit Algorithm for Your Case (Water Example)

For H₂O (3 atoms, 9 coordinates, 3 translational + 3 rotational = 3 vibrational modes):

**Internal coords**: 2 O-H bonds, 1 H-O-H angle → 3 parameters: $k_1, k_2, k_3$

**Basis matrices** (at equilibrium):
- $\mathbf{B}_1^{(b)} = 2\hat{\mathbf{u}}_{OH_1}\hat{\mathbf{u}}_{OH_1}^T$ (in $9\times9$, sparse)
- $\mathbf{B}_2^{(b)} = 2\hat{\mathbf{u}}_{OH_2}\hat{\mathbf{u}}_{OH_2}^T$  
- $\mathbf{B}_3^{(a)} = $ angle geometric tensor

**Project out translations/rotations** using the Eckart frame (or simply work in the 3 vibrational subspace).

**Solve**: $\mathbf{G}_{3\times3} \mathbf{k} = \mathbf{b}_3$

## 10. Extension: Non-Equilibrium Fitting

If you want to fit the **entire energy surface** around minimum (not just Hessian):

Sample $M$ displacements $\{\delta\mathbf{R}^{(s)}\}$, compute:
$$E^{\text{ab}}(\mathbf{R}_0 + \delta\mathbf{R}^{(s)}) \approx E_0 + \frac{1}{2}\delta\mathbf{R}^{(s)T} \mathbf{H}^{\text{ab}} \delta\mathbf{R}^{(s)}$$

Fit force field energy at these points:
$$\min_{\{k\}, \{l_0\}, \{\theta_0\}} \sum_s \left[E^{\text{FF}}(\mathbf{R}_0 + \delta\mathbf{R}^{(s)}; \{k,l_0,\theta_0\}) - E^{\text{ab},(s)}\right]^2$$

This is **nonlinear** in $l_0, \theta_0$ but still **linear in $k$** at fixed geometry. Use alternating optimization or full nonlinear least squares (Levenberg-Marquardt).

## Summary

The key equations are:

$$\mathbf{H}^{\text{FF}} = \sum_m k_m \mathbf{B}_m, \quad \frac{\partial \mathbf{H}^{\text{FF}}}{\partial k_m} = \mathbf{B}_m$$

$$\min_{\mathbf{k}} \|\mathbf{H}^{\text{ab}} - \sum_m k_m \mathbf{B}_m\|_F^2 \Rightarrow \mathbf{G}\mathbf{k} = \mathbf{b}$$

This is a **linear least squares problem** that becomes ill-conditioned for larger molecules. Use **Tikhonov regularization** or **SVD truncation**. For production, packages like **SciPy**, **scikit-learn**, or dedicated chemistry tools (OpenFF, ForceBalance) handle this.

Would you like me to derive the explicit $\mathbf{B}_m$ matrices for methane or work through a numerical example?

---

# Implementation Report: Hessian Parameter Fitting Refactor (May 2026)

## Overview

Refactored the Python-C++ interface for UFF Hessian basis matrices and stiffness parameters to consolidate many small getter/setter functions into two clean, rich functions. This reduces namespace pollution and makes the API easier to navigate.

## What Was Done

### 1. C++ Consolidation (`cpp/libs/Molecular/MMFF_lib.cpp`)

**Removed functions** (namespace pollution):
- `getNBondsHess()`, `getNAnglesHess()`, `getNAtomsHess()`
- `getBondHessianBases()`, `getAngleHessianBases()`
- `getBondAtomsHess()`, `getAngleAtomsHess()`
- `getBondParamsHess()`, `getAngleParamsHess()`
- `setBondParamsHess()`, `setAngleParamsHess()`

**Added consolidated functions**:
```cpp
// Get all Hessian fitting data in ONE call
void getHessianContext(
    int* n_bonds, int* n_angles, int* n_atoms,
    double* bond_bases,    // [nbonds*36]  (4,3,3) per bond: [Bii,Bij,Bji,Bjj]
    int* bond_atoms,        // [nbonds*2]   (2,) per bond: [i,j]
    double* bond_params,    // [nbonds]     k values
    double* angle_bases,   // [nangles*81] (9,9) per angle
    int* angle_atoms,      // [nangles*3]  (3,) per angle: [i,j,k]
    double* angle_params   // [nangles]     k values
);

// Set stiffness parameters in ONE call
void setParams(
    int n_bonds, const double* bond_params,
    int n_angles, const double* angle_params
);

// Kept for validation (assembles H from k vector)
void assembleHessianFromParams(int n_params, const double* params, double* out_hessian);
```

**Key design**: All output arrays are optional (pass `nullptr` to skip). First call gets counts, second call fills arrays.

### 2. Python Consolidation (`pyBall/MMFF.py`)

**Removed wrappers** (corresponding to C++ removals)

**Added consolidated wrappers**:
```python
def getHessianContext(get_bases=True, get_atoms=True, get_params=True):
    """
    Returns dict with:
        n_bonds, n_angles, n_atoms
        bond_bases (nbonds,4,3,3) or None
        bond_atoms (nbonds,2) or None
        bond_params (nbonds,) or None
        angle_bases (nangles,9,9) or None
        angle_atoms (nangles,3) or None
        angle_params (nangles,) or None
    """

def setParams(bond_params=None, angle_params=None):
    """Update stiffness parameters (optional per type)"""
```

**Also added**: `byref` import to fix `NameError`

### 3. Python Fitting Helper (`pyBall/FTIR.py`)

Added convenience wrapper:
```python
def fit_hessian_parameters_ctx(H_target, ctx, lam=1e-6, method='ridge'):
    """
    Fit stiffness parameters using context dict from getHessianContext.
    Returns dict with k_opt, k_bonds, k_angles, residual, r2, condition_number.
    """
```

### 4. Test Script (`tests/tMMFF/test_hessian_fitting.py`)

Complete workflow test:
1. **Init UFF** (bMMFF=True, bUFF=True) with proper data paths
2. **getBuffs_UFF()** - print topology/params
3. **setSwitchesUFF** - enable bonds+angles, disable non-bonded
4. **Relax geometry** with `run()` to ensure minimum
5. **getHessianContext()** - print all basis matrices, atom indices, params
6. **assembleHessianFromParams()** - print assembled Hessian
7. **Linear fit** - recover k values, validate (error ~1e-5, R²=1.0)

## Caveats and Problems Encountered

### 1. Build System Confusion
**Problem**: `/home/prokop/git/FireCore/cpp/Build` is a symlink to `Build-asan`. Using `make -C /home/prokop/git/FireCore/cpp/Build` builds the ASan version, not the optimized version. This caused ASan runtime errors when loading the library from Python.

**Solution**: Always use explicit paths:
- Build: `make -C /home/prokop/git/FireCore/cpp/Build MMFF_lib` (uses symlink → ASan)
- Or explicitly: `make -C /home/prokop/git/FireCore/cpp/Build-asan MMFF_lib`
- Run with ASan preload: `LD_PRELOAD=$(g++ -print-file-name=libasan.so):$(g++ -print-file-name=libfftw3.so) LSAN_OPTIONS=detect_leaks=0 python3 script.py`

### 2. UFF Initialization Requirements
**Problem**: Initial attempts used `bMMFF=False, bBUFF=True` which doesn't work. UFF requires `bMMFF=True` to activate the UFF subsystem.

**Solution**: Follow pattern from `tests/tUFF/test_UFF_multi.py`:
```python
MMFF.init(
    xyz_name=xyz_path,
    sElementTypes="data/ElementTypes.dat",
    sAtomTypes="data/AtomTypes.dat",
    sBondTypes="data/BondTypes.dat",
    sAngleTypes="data/AngleTypes.dat",
    sDihedralTypes="data/DihedralTypes.dat",
    bMMFF=True,   # MUST be True to use UFF
    bBUFF=True,
    bEpairs=False,
    nPBC=(1,1,1),
)
```

### 3. Buffer Access
**Problem**: `getBuffs()` is for MMFFsp3_loc, not UFF. UFF uses `getBuffs_UFF()`.

**Solution**: Use `MMFF.getBuffs_UFF()` which exposes UFF-specific buffers (bonAtoms, angAtoms, bonParams, angParams, neighs, neighBs, hneigh, etc.).

### 4. Switch Configuration
**Problem**: Need to enable UFF components and disable non-bonded forces.

**Solution**:
```python
MMFF.setSwitches(NonBonded=-1, SurfAtoms=-1, GridFF=-1)
MMFF.setSwitchesUFF(DoBond=1, DoAngle=1, DoDihedral=-1, DoInversion=-1, DoAssemble=1, SubtractBondNonBond=-1, ClampNonBonded=-1)
```

### 5. Pre-existing C++ Bug
**Problem**: Heap-buffer-overflow in `UFF::evalBonds()` at `UFF.h:725` when `idebug=1`. The debug print uses `bonParams[neighbor_atom_index]` instead of `bonParams[bond_index]`.

**Solution**: Keep `idebug=0` in test. This is a pre-existing bug in the C++ debug code path (not related to our changes). Should be fixed separately by changing line 725 to use bond index instead of neighbor atom index.

### 6. Angle Basis Computation
**Problem**: `getAllAngleHessianBases()` computed NaN values because `hneigh` (neighbor unit vectors) wasn't populated.

**Solution**: Added `bakeAngleNeighs()` call at start of `getAllAngleHessianBases()` to ensure `hneigh` is populated before computing angle gradients.

### 7. Missing Symbol
**Problem**: `scan_manipulation` was declared but not defined, causing undefined symbol error.

**Solution**: Added minimal stub in `MMFF_lib.cpp`:
```cpp
void scan_manipulation(...) {
    printf("MMFF_lib::scan_manipulation(...) stub\n");
}
```

## How to Use the Script

### Run the test:
```bash
cd /home/prokop/git/FireCore/tests/tMMFF
./run.sh
```

Or directly:
```bash
cd /home/prokop/git/FireCore/tests/tMMFF
LD_PRELOAD=$(g++ -print-file-name=libasan.so):$(g++ -print-file-name=libfftw3.so) LSAN_OPTIONS=detect_leaks=0 python3 test_hessian_fitting.py
```

### Use in your own code:
```python
from pyBall import MMFF, FTIR

# 1. Init UFF
MMFF.setVerbosity(verbosity=1, idebug=0)
MMFF.init(xyz_name="molecule.mol2", bMMFF=True, bBUFF=True, ...)
MMFF.getBuffs_UFF()

# 2. Configure switches
MMFF.setSwitches(NonBonded=-1, SurfAtoms=-1, GridFF=-1)
MMFF.setSwitchesUFF(DoBond=1, DoAngle=1, DoDihedral=-1, DoInversion=-1, DoAssemble=1, ...)

# 3. Relax to minimum
MMFF.setTrjName("trj.xyz", savePerNsteps=10)
MMFF.run(nstepMax=1000, dt=0.05, Fconv=1e-6, ialg=2, damping=0.1)

# 4. Get Hessian context
ctx = MMFF.getHessianContext(get_bases=True, get_atoms=True, get_params=True)

# 5. Assemble Hessian from current params
k_current = np.concatenate([ctx['bond_params'], ctx['angle_params']])
H = MMFF.assembleHessianFromParams(k_current, ctx['n_atoms'])

# 6. Fit parameters (e.g., to DFT Hessian)
result = FTIR.fit_hessian_parameters_ctx(H_DFT, ctx, lam=1e-6, method='ridge')

# 7. Set fitted parameters
MMFF.setParams(result['k_bonds'], result['k_angles'])
```

## Test Results

On H2O (3 atoms, 2 bonds, 1 angle):
- Relaxation converged in 130 steps (max|F| = 3.8e-7)
- Hessian assembled successfully (norm = 137.7)
- Linear fit recovered original k values:
  - Bond k error: 3.7e-6
  - Angle k error: 9.4e-6
  - R² = 1.000000
  - Condition number = 13.8 (well-conditioned)

## Files Modified

- `cpp/libs/Molecular/MMFF_lib.cpp` - Consolidated C++ API
- `cpp/common/molecular/UFF.h` - Added `bakeAngleNeighs()` call in `getAllAngleHessianBases()`
- `pyBall/MMFF.py` - Consolidated Python wrappers, added `byref` import
- `pyBall/FTIR.py` - Added `fit_hessian_parameters_ctx()` convenience wrapper
- `tests/tMMFF/test_hessian_fitting.py` - Complete workflow test with rich debug output
- `tests/tMMFF/run.sh` - Added test_hessian_fitting.py to run script





---

# USER

I'm thinking how to improve our forcefield. Currently we are starting from UFF and refit stiffness constants to fit vibration spectrum from DFT reference.

In fact I simplify UFF, I have 2 bond length stiffness konstants (Si-Si, Si-H) and 3 angles (Si-Si-Si, Si-Si-H, H-Si-H). Do you think I should ad something more? But I don't want to overcomplicate it as fitting is badly conditioned. I'm considering adding 3rd neighbors Si-Si-Si-Si or H-Si-Si-H but as distance rathe than dihedral. Another option is care about non covalent interactions mostly steric hindrance (Pauli repulsion) between surface hydrogen atoms which often bump to each other and which is perhaps critical for vibration spectrum of surface Si-H groups at 2000 cm-1.

Question is how to fit it well, I cannot do linear fit or even compute variational derivatives. I currently just randomly aneal the stiffens constants until I get close results.

Can you figure out better way how to do that fit?

I was thinking to fit directly local elements of the hessian rather than fit resulting vibration frequency. But the DFT hessian looks qualitatively very different it is much deneser (atoms far apart interact) while in my forcefield only the nearest and next nearest neighbors have nonzero elements.

Does it make sense if I compare just the nearest neighbors and neglected the rest? (This will also save me time computing the reference.

---

# Chat GPT 5.6 sol

Yes. I think you can improve both the **force field itself** and, more importantly, the **fitting procedure** quite substantially without adding many parameters.

The biggest point is this:

> **For your present harmonic bond/angle force field, fitting the Hessian can actually be made a linear least-squares problem.**

You do not need random annealing, and you do not need analytic derivatives with respect to the parameters.

Your current FF is already close to the ideal case for this.

## 1. Your current five parameters are a reasonable minimal basis

You have

$$
k_{\rm SiSi},\quad k_{\rm SiH}
$$

and

$$
k_{\rm SiSiSi},\quad
k_{\rm SiSiH},\quad
k_{\rm HSiH}.
$$

I would **not immediately add generic 1–4 springs**. First I would determine exactly what part of the DFT Hessian these five terms fail to reproduce.

In fact, your previous fitting work already suggests something interesting: environment typing improved the spectra substantially, while generic stretch–stretch terms in the Si network did little; the largest useful extra couplings were associated with H–Si–H stretch/bend physics. 

For the specific problem you care about now—**splitting and dispersion of the ~2000 cm\(^{-1}\) surface Si–H band**—I think there is one addition more natural than a generic 1–4 distance spring:

$$
\boxed{
E_{\rm SS}
=
\sum_{\langle a,b\rangle}
J_{ab}\,q_aq_b
}
$$

where \(q_a=\Delta r_{\rm SiH,a}\) is the stretch of one Si–H bond and \(q_b\) the stretch of a neighboring Si–H bond.

This is literally the mechanical analogue of a vibrational tight-binding hopping.

For a pair,

$$
E=
\frac12 k(q_1^2+q_2^2)+Jq_1q_2
$$

gives symmetric/antisymmetric eigenvalues

$$
\lambda_\pm=k\pm J.
$$

So **\(J\) controls the splitting directly**.

You might need only two such couplings:

* \(J_{\rm sameSi}\): between the two H bonds of SiH\(_2\);
* \(J_{\rm neighSi}\): between Si–H groups on neighboring surface Si atoms.

That is much more directly related to the band dispersion you want to describe than `Si-Si-Si-Si` 1–4 distance springs.

---

# 2. Separate two different pieces of physics

I think there are actually two effects hiding in your question.

### A. Harmonic electronic/mechanical coupling

One Si–H oscillator changes the force constant felt by the neighboring Si–H oscillator.

That is naturally described by

$$
J q_i q_j.
$$

It produces splitting even in a perfectly relaxed structure and costs only one or two parameters.

### B. Steric H···H interaction

Two surface hydrogens get close and repel each other.

That is different physics and should ideally be represented as an actual nonbonded potential,

$$
V_{\rm HH}(r).
$$

For example, even just

$$
V_{\rm HH}(r)=A e^{-br}
$$

could be enough.

This can be extremely relevant on SiH\(_2\)-rich {100}, edges, corners, etc. because an Si–H stretch changes H···H distances. Its Hessian therefore creates an **environment-dependent coupling automatically**.

And unlike a constant \(J\), it gets stronger when hydrogens are crowded.

I would probably fix \(b\) initially to a plausible value and fit only \(A\). Otherwise \(A,b\) will be strongly correlated.

So my minimal extension would be:

$$
\boxed{
5\ {\rm existing}
+
1\ {\rm SiH\!-\!SiH\ coupling}
+
1\ {\rm H-H\ repulsion\ amplitude}
}
$$

before adding generic 1–4 springs.

---

# 3. I would postpone the Si–Si–Si–Si 1–4 distance

It may eventually be useful, particularly because your simple bond+angle model is known to miss some soft/cage modes. But mechanically a 1–4 distance spring is mostly acting as an indirect torsion/longitudinal constraint.

It will mostly help:

* cage rigidity;
* torsional/shear modes;
* lower frequencies.

It seems less likely to be the missing physics behind the **Si–H stretching manifold**.

Likewise an H–Si–Si–H 1–4 spring might work empirically, but I would rather model the actual H···H distance directly if sterics is the intended physics.

---

# 4. The crucial fitting trick: your Hessian is linear in stiffness parameters

Suppose

$$
E(\mathbf R)
=
\sum_p k_p\,E_p(\mathbf R).
$$

Then at a fixed geometry,

$$
H^{FF}
=
\frac{\partial^2E}{\partial \mathbf R^2}
=
\sum_p k_p H_p
$$

where

$$
H_p =
\frac{\partial^2E_p}{\partial \mathbf R^2}.
$$

Therefore

$$
\boxed{
H^{FF}(\mathbf k)
=
\sum_p k_p H_p
}
$$

is linear in the parameters.

You don't even have to derive \(H_p\).

Just calculate the FF Hessian with:

```text
k1=1, all others=0 -> H1
k2=1, all others=0 -> H2
...
```

or finite-difference an existing Hessian implementation with respect to each parameter.

With only 5–10 parameters this cost is negligible.

Then the fit is simply

$$
\min_{\mathbf k}
\left\|
\sum_p k_p H_p-H^{DFT}
\right\|^2.
$$

This is ordinary linear least squares.

And you immediately get something random annealing cannot give you: the **singular values of the fitting matrix**, which tell you whether the parameters are actually identifiable.

---

# 5. But I would NOT fit the entire Cartesian DFT Hessian

Your observation is exactly correct.

The DFT Hessian is dense:

$$
H_{ij}^{DFT}\ne0
$$

even for relatively distant atoms, whereas your local FF has

$$
H_{ij}^{FF}=0
$$

outside its interaction graph.

If you minimize the Frobenius norm of the complete matrices,

$$
\|H^{FF}-H^{DFT}\|_F^2,
$$

you punish the FF for not representing physics it **by construction cannot represent**.

That can distort the local parameters.

But I also would **not simply delete all distant DFT matrix elements**. There is a better solution.

---

# 6. Fit a physically selected **projection of the Hessian**

This is, I think, the key idea.

Construct displacement vectors representing the coordinates your model is supposed to reproduce:

$$
Q=
[
q_{\rm SiH},
q_{\rm SiSi},
q_{\rm angle},
\ldots
].
$$

After suitable orthonormalization, compare

$$
\boxed{
K_Q = Q^T D Q
}
$$

where

$$
D=M^{-1/2}HM^{-1/2}.
$$

Do this for DFT and FF:

$$
K_Q^{DFT}=Q^T D^{DFT}Q,
$$

$$
K_Q^{FF}
=
\sum_p k_p Q^TD_pQ.
$$

Then solve

$$
\min_{\mathbf k}
\left\|
K_Q^{FF}-K_Q^{DFT}
\right\|_F^2.
$$

Still completely linear.

And now the dense DFT Hessian is **not a problem at all**: its long-range electronic response is integrated into the effective force constants within the coordinate subspace you care about.

This is very close conceptually to the Wilson-coordinate fitting machinery you already developed in the project.

---

# 7. For the Si–H problem specifically, use the **Si–H stretch subspace**

This gives a particularly elegant fit.

For every Si–H bond define a normalized local stretch coordinate

$$
g_a.
$$

Collect them into \(G\), orthogonalize them, and calculate

$$
K_{\rm H}^{DFT}=G^T D^{DFT}G.
$$

This matrix itself is scientifically interesting.

Its diagonal

$$
(K_H)_{aa}
$$

tells you the effective local Si–H stiffness.

Its off-diagonal terms

$$
(K_H)_{ab}
$$

measure coupling between Si–H oscillators.

So before fitting anything, plot

$$
(K_H^{DFT})_{ab}
$$

against:

* SiH–SiH graph distance;
* H···H distance;
* Si···Si distance;
* whether they share the same Si;
* same facet / across edge;
* relative bond orientation.

That tells you what additional FF term nature is asking for.

For example, if you discover

$$
K_{ab}^{DFT}\approx J_0
$$

for adjacent hydrides regardless of H···H distance, add \(Jq_aq_b\).

If instead you see something like

$$
K_{ab}^{DFT}
\sim e^{-b r_{HH}},
$$

then H···H Pauli repulsion is clearly the relevant model.

This is much better than guessing which 1–4 spring to add.

---

# 8. There is another extremely useful Hessian fitting objective: fit in the **DFT eigenmode basis**

Suppose \(U\) contains the DFT modes in the region you care about, say the Si–H modes around 1800–2200 cm\(^{-1}\).

DFT satisfies

$$
U^T D^{DFT} U
=
\mathrm{diag}(\omega_1^2,\ldots,\omega_m^2).
$$

Evaluate your FF in *the same fixed DFT basis*:

$$
K^{FF}_{mn}
=
u_m^TD^{FF}u_n.
$$

Then fit

$$
\boxed{
\min_{\mathbf k}
\left\|
U^TD^{FF}(\mathbf k)U
-
\mathrm{diag}(\omega_{DFT}^2)
\right\|^2
}
$$

Again, this is linear in \(\mathbf k\).

This has several advantages over fitting frequencies:

* no diagonalization inside the optimizer;
* no mode crossing problem;
* no Hungarian assignment every iteration;
* off-diagonal elements penalize wrong mode mixing;
* you can select exactly the spectral region important experimentally.

This may actually be the easiest first implementation.

---

# 9. I would combine several targets

For example,

$$
L =
w_H L_{\rm SiH}
+
w_{\rm loc}L_{\rm local}
+
w_{\rm mode}L_{\rm modes}
+
\lambda L_{\rm reg}.
$$

Where:

### Surface stretch matrix

$$
L_{\rm SiH}
=
\|G_H^T(D^{FF}-D^{DFT})G_H\|_F^2.
$$

This strongly constrains the experimental 2000 cm\(^{-1}\) problem.

### General local valence space

Include Si–Si stretches and angle bends:

$$
L_{\rm local}
=
\|Q^T(D^{FF}-D^{DFT})Q\|_F^2.
$$

### Selected DFT modes

$$
L_{\rm modes}
=
\|U^TD^{FF}U-\Omega_{DFT}^2\|_F^2.
$$

### Regularization

Keep parameters near sensible UFF/previous-fit values:

$$
L_{\rm reg}
=
\sum_p
\left(
\frac{k_p-k_p^0}{\sigma_p}
\right)^2.
$$

All of these remain a **quadratic linear least-squares problem**.

---

# 10. Use bounded / regularized least squares, not unconstrained inversion

Your conditioning concern is very real.

Build

$$
A\mathbf k\simeq\mathbf b
$$

and inspect the SVD:

$$
A=U\Sigma V^T.
$$

If two parameters give almost the same Hessian change, a small singular value will expose it immediately.

Then either:

* merge the parameters;
* fix one;
* or regularize them toward common values.

For ordinary positive stiffnesses use bounds

$$
k_p\ge0.
$$

Cross-couplings \(J\) can be signed.

In Python, `scipy.optimize.lsq_linear()` is already enough.

For ridge regularization just augment

$$
A'=
\begin{bmatrix}
A\\
\sqrt{\lambda}L
\end{bmatrix},
\qquad
b'=
\begin{bmatrix}
b\\
\sqrt{\lambda}Lk_0
\end{bmatrix}.
$$

No nonlinear optimizer is needed.

---

# 11. Should you compute only nearest-neighbor pieces of the DFT Hessian?

**Yes, potentially—but I would change exactly what you compute.**

Simply calculating only Cartesian blocks

$$
H_{ij}
$$

for bonded neighbors is useful, but it is not the most efficient target.

Instead compute the **directional Hessian matrix elements that correspond to your internal coordinates**:

$$
K_{ab}
=
g_a^T H g_b.
$$

You do **not need the full Hessian** to obtain these.

If \(g_a\) is one local displacement pattern, displace the molecule by

$$
\mathbf R_\pm
=
\mathbf R_0\pm\delta g_a
$$

and calculate the forces.

Then

$$
Hg_a
\simeq
-\frac{
F(\mathbf R_+)-F(\mathbf R_-)
}{2\delta}.
$$

From that one pair of DFT force calculations you obtain

$$
K_{ba}
=
g_b^T H g_a
$$

for **every** \(b\).

So if you have \(m\) selected internal coordinates, you need approximately

$$
2m
$$

DFT force calculations rather than

$$
6N
$$

Cartesian force calculations.

That is particularly attractive for larger PySCF/GPAW references.

For example, if your main question is the surface Si–H manifold and there are 20 surface Si–H bonds, 40 DFT force calculations give the complete

$$
20\times20
$$

coupling matrix between them.

That is much more useful than a sparse subset of arbitrary Cartesian Hessian entries.

---

# 12. There is one subtle point: local terms vs prestress

For the cleanest Hessian fit, I would initially evaluate all methods at **one common DFT geometry** and set the harmonic FF equilibrium values to that geometry:

$$
r_0=r^{DFT},\qquad
\theta_0=\theta^{DFT}.
$$

Then bond/angle terms have zero force individually and the Hessian reflects pure stiffness.

This makes the fit almost perfectly linear and interpretable.

Nonbonded H···H repulsion is different: it is nonzero-force at the DFT geometry, and the equilibrium comes from cancellation with the covalent forces.

So I would fit in two stages:

**Stage 1: harmonic model**

Fit the five valence parameters + perhaps \(J_{\rm SiH-SiH}\) to projected DFT Hessians.

**Stage 2: sterics/geometry**

Fit H···H repulsion using DFT forces/energies from geometries where neighboring hydrides are pushed toward one another.

Then recompute the Hessian and see how much the surface-band prediction improves.

That prevents the H···H repulsion from destabilizing the otherwise well-conditioned harmonic fit.

---

# 13. A very concrete development sequence

I would do the following before adding much more physics:

1. Take **SiH\(_4\), Si\(_2\)H\(_6\), Si\(_{10}\)H\(_{16}\)** and one small motif-complete NC.

2. Put all Hessians at the same DFT geometry.

3. Generate the five unit-parameter Hessians

$$
H_p.
$$

4. Fit the five existing parameters by projected linear least squares.

5. Examine residual

$$
R=D^{DFT}-D^{FF}.
$$

6. Project the residual onto the Si–H stretch basis:

$$
R_H=G_H^TRG_H.
$$

7. Plot off-diagonal residuals against H···H distance and surface topology.

This one plot will tell you what to add.

I predict there are three plausible outcomes:

* systematic nearest-neighbor Si–H coupling → add \(Jq_iq_j\);
* strong exponential dependence on H···H distance → add H–H repulsion;
* residual correlates mainly with stretch/bend coordinates on SiH\(_2\) → add one \(r\)-\(\theta\) cross term.

Only if the residual instead lives mostly in the Si backbone and soft modes would I add a Si–Si–Si–Si 1–4/torsion term.

That gives you a disciplined **one-extra-parameter-at-a-time model selection procedure**, rather than guessing more UFF terms.

And because every candidate term contributes another known basis Hessian \(H_p\), you can immediately see from the SVD whether the new parameter is independently identifiable before accepting it.


---

# USER

OK, I remember we were already doing some of these things, also I remember we made some localized hessian basis called 
 **Wilson B-matrix / GF method** 
can you review it ? Is it usefull? there were some problems, can we overcome it? Or should we pursure fundamentally different apporach like what you just proposed?

more is here
https://github.com/ProkopHapala/FireCore/tree/prokop

\---

From my coding agent

\---

# USER

OK, about the fitting - I remember we were doing some special ftiing of hessian in local frame (like bonds and angles, it was named after some guy, I do not remember can you find it)?

---

# Cursor Grok 4.6

Yes — that’s **Wilson** (E. Bright Wilson): the **Wilson B-matrix / GF method**. Internal coordinates are bonds and angles; each has a Wilson vector \(b = \partial q / \partial r\).

SSOT: [`doc/Topics/FFfit/HessianFitting_Theory.md`](doc/Topics/FFfit/HessianFitting_Theory.md)

Two different uses, don’t mix them:

1. **Least-norm Wilson \(F = C^{+T} D C^+\)** — diagnostic only. The diagonal is **not** a unique per-bond DFT stiffness (redundant bonds+angles). That was the “central Si–Si looks soft” artifact.
2. **Gauge-invariant Wilson row-space residual** — the actual fit, together with modes + graph-local Hessian. That is `fit_hybrid_hessian()` in `pyBall/FFfit.py` / `FFfit_utils.py`, CLI `tests/tSiNCs/test_FFfit.py`.

That is **not** the recent own-min stretch RMSE (`fit_mmff_kss_pyscf.py`). Wilson fitting is frozen-geometry Hessian matching (FFfit). Own-min stretch RMSE is FIRE then spectrum.

The recent C cube vs octahedron anneal never used this Wilson machinery.

---

# USER

OK, do we have some report of the results we get from it I remmber we had some rpoblems. Do we have some topical audit of this hessian fitting ?

---

# Cursor Grok 4.6

Yes on the **results**, no on a **dedicated FFfit topical audit**.

**Results report:** [`tests/tSiNCs/SiNCs_FFfit_summary.md`](tests/tSiNCs/SiNCs_FFfit_summary.md) (2026-07-10). Six PySCF systems: `SiH4` + `Si_R3p8`…`Si_R6p0`. Hybrid objective (modes + local Hessian + Wilson row-space).

| Model | Mean freq RMSE |
|-------|----------------|
| Elemental bond+angle | 41 cm⁻¹ |
| Hierarchical Si / SiH / SiH₂ / SiH₃ | **32 cm⁻¹** |
| + stretch–bend | Hessian relFrob 6.9% → 5.9%; spectrum barely better |

Fitted Si–Si \(k\) stayed ~**9.2–9.8 eV/Å²**. That is frozen-geometry Hessian matching, **not** an own-min FTIR spectrum.

**The problem you remember:** Wilson least-norm \(F\) made **central Si–Si look soft**. That diagonal is gauge-dependent in a redundant bond/angle basis. It is **diagnostic only**. Hierarchy regularization was the fix so bulk bonds are not dumped into acoustic modes. Theory: [`doc/Topics/FFfit/HessianFitting_Theory.md`](doc/Topics/FFfit/HessianFitting_Theory.md) §6.3–6.4. Chat log of the artifact: [`HessianFit.chat.md`](doc/Topics/FFfit/HessianFit.chat.md).

**Other recorded limits:** UFF torsions **exist but are not recommended** (indefinite prestress) — [`Dihedral_Torsion.md`](doc/Topics/FFfit/Dihedral_Torsion.md). Cross terms only with `--equilibrium local`. This benchmark is **small all-Si spheres**, not the L1 cube/octa PBE set.

**Audit status:** there is **no** `doc/topical_audit/FFfit.md`. FFfit is a subsection of [`topical_audit.md`](doc/topical_audit/topical_audit.md) §2d, plus [`SiNCs.md`](doc/topical_audit/SiNCs.md) §3.4 and graph/Wilson inventory in [`molecular_topology.md`](doc/topical_audit/molecular_topology.md). `Dihedral_Torsion.md` is a TopicalAudit, but only for dihedrals.

---

# Chat GPT 5.6 sol

Yes. I reviewed the current Wilson/GF machinery, the fitting code, the results report, and the current theory document in `prokop`. My conclusion is:

> **Do not abandon the Wilson machinery. But change how we use it.**
>
> The old attempt to extract a “DFT spring constant for each bond” from the pseudoinverse Wilson \(F\)-matrix was fundamentally ambiguous. The current code already fixes that mathematical mistake. However, the present **full Wilson row-space objective is more global than it looks**, and for your new Si–H surface problem I would replace it by **restricted, chemically chosen Wilson subspaces**—especially the Si–H stretch manifold.
>
> In fact, what I proposed in the previous message is not a fundamentally different approach. It is essentially the **localized/restricted version of the Wilson approach you already implemented**.

There is quite a lot worth reusing.

## What you already have is better than I realized

The underlying FFfit architecture is sound. At fixed geometry, it explicitly uses

$$
H^{FF}(\mathbf k)=\sum_p k_p A_p,
$$

with the exact sensitivity

$$
A_p=f''(q)\,bb^T+f'(q)\,C,
$$

where \(b=\partial q/\partial r\) is precisely the Wilson vector and the second term accounts for prestress. The code already solves the resulting problem with bounded linear least squares, column scaling and Tikhonov regularization. 

So for stiffness constants, **the random annealing should no longer be the primary fitter**. It is useful for testing a nonlinear end-to-end model after relaxation, but not for the frozen-geometry stiffness fit.

Your current production test actually combines three linear objectives: DFT normal modes, graph-local Cartesian Hessian blocks, and a Wilson internal-coordinate objective.  The graph-local part already implements exactly the idea we discussed before: ignore DFT Hessian blocks beyond the graph range that the local FF cannot represent. For bond+angle FF it uses graph distance \(\le2\); distance 3 is available when testing 1–4 interactions. 

And numerically the present fit isn't actually terribly ill-conditioned after scaling and regularization: the reported conditions are only about 1.3–1.7. The 5-parameter elemental model gives ~41 cm\(^{-1}\) mean frequency RMSE; hierarchical environment typing lowers it to ~32 cm\(^{-1}\), whereas adding many cross terms gives relatively modest further spectral improvement. 

That is a good starting point.

---

# The problem with the original Wilson \(F\)-matrix

Suppose

$$
q=B x,
$$

and in mass-weighted coordinates

$$
C=B M^{-1/2}.
$$

If the harmonic internal-coordinate energy is

$$
E=\frac12 q^TFq,
$$

then

$$
D=C^TFC,
\qquad
D=M^{-1/2}HM^{-1/2}.
$$

You tried to invert this as

$$
F=C^{+T}DC^+.
$$

The code still provides this as `internal_hessian_projection()`, but now correctly warns that when bonds + angles form a **redundant internal-coordinate set**, the resulting \(F\) is merely the minimum-norm representative. Its individual diagonal entries depend on the coordinate gauge and SVD cutoff; they are not uniquely defined bond stiffnesses. 

This explains the old phenomenon where a central Si–Si bond appeared strangely soft. The curvature wasn't necessarily physically soft; the pseudoinverse was free to distribute the same Cartesian curvature among many mutually redundant bond/angle coordinates.

So:

$$
\boxed{\text{Do not fit or interpret }F_{ii}^{DFT}\text{ as individual spring constants.}}
$$

That problem cannot really be “fixed” because it is not primarily numerical. It is non-uniqueness.

The current report correctly relegates this least-norm \(F\) to a diagnostic only. 

---

# The current row-space fix is mathematically correct — but there is another subtlety

The new code does something much better.

Take the SVD

$$
C=U\Sigma V^T,
$$

retain the independent directions, and define

$$
Q=V_r.
$$

Then fit

$$
L_W=
\left\|
Q^T(D^{FF}-D^{DFT})Q
\right\|_F^2.
$$

That is invariant to rescaling or nonsingular recombination of your redundant internal coordinates. This is the “gauge-invariant Wilson row-space” fit now used in `assemble_hybrid_hessian_system()`.  

That completely avoids the “central bond soft because of gauge” problem.

However, I see an important conceptual issue that isn't emphasized in the existing documentation.

For a normal connected nonlinear molecule, if the collection of all bonds + angles is sufficiently complete,

$$
\operatorname{rank}(B)=3N-6.
$$

Then \(Q\) spans essentially the **entire vibrational Cartesian space**.

Define

$$
P_{\rm vib}=QQ^T.
$$

Then

$$
\left\|Q^T\Delta DQ\right\|_F
=
\left\|P_{\rm vib}\Delta D P_{\rm vib}\right\|_F.
$$

In other words, once the Wilson bond+angle set spans all \(3N-6\) vibrational DOFs, the “internal-coordinate objective” has **lost most of its chemical locality**. It is basically the full vibrational Hessian norm with translations/rotations removed.

That is not bad—it is an excellent stable global Hessian objective—but it isn't really telling the fitter:

> reproduce these bonds and these angles.

It is saying:

> reproduce all curvature within the vibrational subspace.

And that makes it partly redundant with your **mode-space objective**.

Indeed, your mode objective evaluates

$$
K=V^TDV
$$

in the DFT eigenmode basis. If you used all modes, equal weights, and full off-diagonal mixing weight 1, then

$$
\|V^T\Delta DV\|_F
$$

would be mathematically the same global vibrational Frobenius norm, merely expressed in another orthonormal basis.

The current implementation differs because it deliberately frequency-weights the mode diagonals and uses only 0.1 by default for off-diagonal mode mixing. 

So the three current objectives are really:

$$
\begin{array}{ll}
\text{MODE} & \text{global vibrational Hessian, spectrally reweighted},\\
\text{WILSON row space} & \text{global vibrational Hessian, essentially unweighted},\\
\text{LOCAL} & \text{actually local Cartesian curvature}.
\end{array}
$$

This is why I would **not throw Wilson away**, but I also wouldn't automatically give all three weight 1 as the current CLI does. The current defaults indeed set mode/local/internal all to 1. 

For the new surface-FTIR fit, I think that is unnecessarily opaque.

---

# What I would do instead: **restricted Wilson spaces**

This is where Wilson becomes extremely useful for your current physics.

Instead of

$$
B=
B_{\rm all\ bonds+angles},
$$

construct several chemically meaningful \(B\)'s.

For example

$$
B_H =
\begin{bmatrix}
b_{\mathrm{SiH}_1}\\
b_{\mathrm{SiH}_2}\\
\vdots
\end{bmatrix},
$$

containing **only Si–H stretches**.

Or

$$
B_{\rm surf}
=
\begin{bmatrix}
B_{\rm SiH}\\
B_{\rm HSiH}\\
B_{\rm SiSiH}
\end{bmatrix}
$$

for the surface valence coordinates.

Or separately

$$
B_{\rm core}=
[B_{\rm SiSi},B_{\rm SiSiSi}].
$$

Now the row-space has dimension perhaps 20–100, rather than \(3N-6\), and

$$
Q_H^T D Q_H
$$

really does ask:

> What is the curvature inside the Si–H stretching manifold?

That is *exactly* what we need for Katerina's 2000–2200 cm\(^{-1}\) region.

I would therefore change the fitting objective conceptually from

$$
L_{\rm internal(all)}
$$

to something like

$$
L=
w_H L_{\rm SiH}
+
w_{\rm surf}L_{\rm surf}
+
w_{\rm local}L_{\rm Cartesian-local}
+
w_{\rm mode}L_{\rm selected\ modes}.
$$

And initially I would probably make \(w_H\) dominant.

---

# Even better: use a **localized orthogonal Wilson basis**

For fitting a subspace, arbitrary \(Q\) from the SVD is fine.

For interpreting coupling between individual Si–H groups, however, SVD vectors become delocalized combinations. We want a basis that remains as close as possible to the physical bonds.

For the Si–H Wilson matrix define

$$
C_H=B_HM^{-1/2}.
$$

Its row Gram matrix is

$$
G_H=C_HC_H^T.
$$

If the Si–H coordinates are linearly independent—which they almost certainly will be—perform symmetric/Löwdin orthogonalization:

$$
\boxed{
R_H=G_H^{-1/2}C_H
}
$$

so that

$$
R_HR_H^T=I.
$$

Unlike arbitrary SVD rotation, this produces orthogonal coordinates that remain maximally close to the original **individual local Si–H stretches**.

Then calculate

$$
\boxed{
K_H=R_HD R_H^T.
}
$$

This gives exactly the matrix I wanted in our previous discussion:

$$
K_H=
\begin{pmatrix}
\epsilon_1 & J_{12} & J_{13}&\cdots\\
J_{12}&\epsilon_2&J_{23}&\cdots\\
\vdots&&&
\end{pmatrix}.
$$

Now:

$$
\epsilon_i=(K_H)_{ii}
$$

is a meaningful effective local oscillator curvature, while

$$
J_{ij}=(K_H)_{ij}
$$

is the coupling between two local Si–H stretches.

This is basically your **Wilson machinery turned into a vibrational tight-binding model**.

And this is much closer to what we actually want to learn from DFT.

---

# This immediately helps decide what term your FF is missing

After fitting your present five parameters, calculate

$$
\Delta K_H
=
K_H^{DFT}-K_H^{FF}.
$$

Then inspect its off-diagonal terms.

Suppose two neighboring monohydride groups are

$$
H_i-\mathrm{Si}_i-\mathrm{Si}_j-H_j.
$$

Plot

$$
(\Delta K_H)_{ij}
$$

against H···H distance, Si···Si distance, bond orientations, facet identity and graph distance.

Then we can distinguish very cleanly:

$$
\Delta K_{ij}\approx\mathrm{constant}
$$

for nearest-neighbour SiH groups → explicit mechanical coupling

$$
E=Jq_iq_j
$$

is appropriate.

If instead

$$
\Delta K_{ij}
$$

strongly tracks \(r_{HH}\), that points toward the H···H Pauli/steric potential.

If the residual is mostly inside SiH\(_2\), correlated with H–Si–H motion, then stretch–bend coupling is the missing term.

This would make the **model design data-driven** instead of deciding beforehand between 1–4 springs, H–H repulsion and cross terms.

---

## One important detail about the cross term you already tested

Your existing “stretch–stretch” implementation is **not quite the coupling I was suggesting for facet dispersion**.

`compute_cross_sensitivity()` loops over angles \(i-j-k\) and couples the stretches of the two bonds sharing central atom \(j\):

$$
E_{rr}=K_{rr}\,q_{ij}q_{jk}.
$$

So for H–Si–H this couples the two hydrides on **one SiH\(_2\)** group. The implementation explicitly constructs the two bonds belonging to each angle and forms

$$
g_1g_2^T+g_2g_1^T.
$$

That is the cross term already tested. 

But the coupling responsible for finite {111} SiH surface bands is more naturally

$$
\boxed{
E_{\rm neigh}
=
J
\,q_{\mathrm{Si_i-H_i}}
\,q_{\mathrm{Si_j-H_j}}
}
$$

across a **Si–Si bond**:

$$
H_i-\mathrm{Si}_i-\mathrm{Si}_j-H_j.
$$

That is a bond-pair coupling over an H–Si–Si–H topology and is *not* the same term.

This matters because the current report says stretch–stretch coupling gave almost no improvement for the Si network, while the largest cross correction was H–Si–H stretch–bend.  But that result does **not** demonstrate that neighboring surface SiH oscillators are uncoupled; it tested a somewhat different coupling topology.

I would explicitly add this H–Si–Si–H bond-pair coupling as a **diagnostic parameter**, not necessarily immediately as a permanent FF term.

---

# Wilson can also save you DFT reference cost

And here the connection to our previous discussion becomes especially nice.

The current `internal_hessian_projection()` requires the complete \(H\).

But mathematically you do not need it.

Choose one relevant normalized displacement \(g_a\), for example a localized Si–H stretch. Compute DFT forces at

$$
x_\pm=x_0\pm\delta g_a.
$$

Then

$$
H^{DFT}g_a
\simeq
-\frac{F(x_+)-F(x_-)}{2\delta}.
$$

One pair of DFT calculations therefore gives an entire **Hessian-vector product**.

You can then calculate all selected couplings

$$
K_{ba}
=
g_b^T H g_a.
$$

If there are \(m\) Si–H stretch coordinates, the complete \(m\times m\) Si–H coupling matrix costs roughly

$$
2m
$$

DFT gradient calculations—not \(6N\).

This partial-reference idea is already captured in your uploaded discussion: 20 selected surface Si–H coordinates would need about 40 force calculations to obtain the complete 20×20 coupling matrix. 

There is an even stronger version for fitting.

You do not actually need to reduce the reference to \(g_b^THg_a\). Fit the whole Hessian-vector response:

$$
\boxed{
\sum_p k_p A_p g_a
\approx
H^{DFT}g_a.
}
$$

Since you already have the analytic sensitivity matrices \(A_p\), this is still linear in \(k_p\).

So each selected QM displacement gives \(3N\) fitting observations.

I really like this formulation because it avoids:

* computing the full QM Hessian;
* pseudoinverting Wilson coordinates;
* comparing irrelevant remote Hessian blocks;
* diagonalizing during fitting.

And the physical directions you probe are chosen deliberately.

---

# I would therefore simplify the existing hybrid fitter rather than replace it

For the next generation I would make the following concrete changes:

1. **Keep `FFfit` sensitivity matrices and bounded linear LSQ unchanged.** They are the valuable core.

2. Keep the least-norm Wilson \(F=C^{+T}DC^+\) **only for visualization/diagnostics**, with the current warning. Never derive individual FF constants from its diagonals.

3. Add `build_wilson_matrix(..., selection=...)`, or simpler specialist builders:
   `SiH_stretches`, `surface_bends`, `core_valence`.

4. Add a `lowdin_wilson_basis()` returning \(R=(CC^T)^{-1/2}C\), along with singular values/condition.

5. Add a `projected_coupling_matrix()` returning

   $$
   K=RDR^T.
   $$

   This should become the primary diagnostic for surface Si–H coupling.

6. Add a **directional/HVP fitting mode**:

   $$
   A_p g_a \rightarrow H_{\rm ref}g_a.
   $$

   It should accept either a stored complete Hessian or forces from ± directional displacements.

7. Change the default new Si-H fit from today's `mode=1, local=1, internal=1` to a chemically targeted objective—roughly **SiH-subspace + graph-local**, with the full mode objective used mainly as validation.

8. Add one candidate `H-Si-Si-H stretch–stretch` coupling parameter and separately an H···H repulsive potential; use \(\Delta K_H\) to decide which is justified.

I would also keep your existing 16-parameter Si/SiH/SiH₂/SiH₃ hierarchy, but **as a diagnostic oracle rather than necessarily the production FF**. It already showed that environment dependence is real, lowering the spectral RMSE from ~41 to ~32 cm\(^{-1}\).  If a simple H···H repulsion or one neighbouring-SiH coupling can reproduce those environment-dependent shifts, you can collapse many empirical subtype parameters back into a much smaller physically motivated model.

So conceptually I would aim for:

$$
\boxed{
\text{DFT directional response}
\;\rightarrow\;
K_{\rm SiH}^{DFT}
\;\rightarrow\;
\text{5-param FF}
\;\rightarrow\;
\Delta K_{\rm SiH}
\;\rightarrow\;
\text{one physically motivated new term}
}
$$

rather than continuing to make the global hybrid loss more sophisticated.

The **Wilson B matrix absolutely remains central**, but mainly as the language for saying *which local displacement subspace we care about*. That is, I think, the clean resolution of the problems you encountered before.








---

# USER

OK, thank you this is extremely valuable disicussionb and your undertaing and insight into the problem is obvious I learn a lot new. 

I completely agree that point of Willson and orironal motivation is to make localized set of chemically meaningfull vibrations ivolving mostly just individual bonds/angles (pairs, triples max quadrupoles of atoms) not global eigenmobeds of whole system! we definitely should restrict it and focus on these. I also very much the odea to restrict the cost of DFT reference so that I can do finite differences just fro pairs of nearby atoms which is ~O(N) elements for large system instead of O(n^2)

However some things are still not quite clear. You say I should simply compare just the nearest neighbor elements of the hessian between forcefield and DFT reference. But if other non-local element in DFT are large and non-zero and I will diagonalize that Hessian, I will get very different eigenfrequencies than form my forcefield even when all nearest neighbors ae exacly the same. Now my physical intuition a bit hit ist limit and I'm confused. If I imagine physically I'm bending the angle or streching the bond (basically holding the atoms) then I think the narest neighbor elements are all what matters. And the global deformation (like holed one side of the crystal or polymer) are just combination of these local which cna be obtained from the local by diagonalization/inversion (linear elestic response to external force RHS). This suggest intuitively I would expect that other than nearest neighbor elements would be near zero also in DFT (I should double check, but if I remember they were not), so why it is? Or how it should be? What should one expect? And how to reconsiliate these two views? Intuitivelly I would think that covalent crystal is relly just system of local springs (especially for large band gap like diamond, where the desity matrix decay very fast - exponentially with exponent being bandgap width). OK, there can be also some electrostatic and van-der-walls but those are like 2-3 orders of maginitude weaker than chemical bonds and we can neglect them.


---








# USER

Yes — this is exactly the subtle point. Your intuition is mostly right, but there are **three different objects that are easy to mentally mix together**:

$$
\boxed{\text{local force constants }H}
\qquad
\boxed{\text{relaxed response }H^{-1}}
\qquad
\boxed{\text{effective Hessian after eliminating hidden DOFs}}
$$

The first can be very local while the latter two are dense.

And I should slightly correct my previous wording: I would **not assume beforehand that matching only nearest-neighbor Cartesian Hessian blocks is sufficient**. I would first test whether the DFT Hessian is *effectively local* for the spectral subspace you care about. Your existing `local_hessian_mask()` already embodies this approximation: for bond+angle FF it keeps graph distance \(\le2\) because those are the blocks the model can actually generate, rather than penalizing the FF for missing longer-range QM blocks. 

## A purely local spring network already has completely global vibrations

Consider the simplest 1D chain with nearest-neighbor springs:

$$
E=\frac{k}{2}\sum_i(x_{i+1}-x_i)^2.
$$

Its Hessian is sparse,

$$
H=
k\begin{pmatrix}
1&-1&0&0&\cdots\\
-1&2&-1&0&\cdots\\
0&-1&2&-1&\cdots\\
\vdots
\end{pmatrix}.
$$

There is **zero direct force constant between distant atoms**.

Nevertheless, diagonalizing it gives perfectly global phonons,

$$
u_j(q)\sim e^{iqj},
$$

with dispersion

$$
\omega^2(q)=\frac{4k}{m}\sin^2\frac q2.
$$

So you absolutely do **not** need long-range Hessian elements to obtain global elastic modes, long-wavelength bending, sound waves, etc.

Likewise, applying a force at one end produces displacement everywhere because

$$
x=H^{-1}f,
$$

and although \(H\) is sparse,

$$
\boxed{H^{-1}\text{ is generally dense}.}
$$

That corresponds almost exactly to your intuition of holding/pushing one part of a crystal and seeing the deformation propagate through the local bonds.

So far, you are completely right.

---

# Then why is a DFT Hessian dense?

There are several reasons, and the first one is conceptually important.

DFT does **not** give you a nuclear Hessian before the electrons respond. It gives you the Born–Oppenheimer energy

$$
E_{\rm BO}(R)=\min_{\psi}E(R,\psi).
$$

The electronic degrees of freedom have already been eliminated.

Imagine schematically that near equilibrium the complete electron+nucleus problem had a Hessian

$$
K=
\begin{pmatrix}
K_{RR}&K_{Re}\\
K_{eR}&K_{ee}
\end{pmatrix}.
$$

When the nuclei move, the electrons relax:

$$
\delta e=-K_{ee}^{-1}K_{eR}\delta R.
$$

Substituting this back gives the nuclear Born–Oppenheimer Hessian

$$
\boxed{
H_{\rm BO}
=
K_{RR}
-
K_{Re}K_{ee}^{-1}K_{eR}.
}
$$

And \(K_{ee}^{-1}\), an electronic susceptibility/response function, need not be sparse.

So even if the **underlying electronic Hamiltonian is local**, integrating out the electrons generates effective interactions between nuclear displacements at different positions.

This is precisely analogous to eliminating intermediate atoms from a spring network.

If

$$
H=
\begin{pmatrix}
H_{AA}&H_{AB}\\
H_{BA}&H_{BB}
\end{pmatrix}
$$

and you allow \(B\) to relax, the effective stiffness seen by \(A\) becomes

$$
\boxed{
H_A^{\rm relaxed}
=
H_{AA}
-
H_{AB}H_{BB}^{-1}H_{BA}.
}
$$

That second term is generally dense.

So there are two distinct notions:

**Clamped-ion coupling:** move atom \(j\), hold every other nucleus fixed, measure force on atom \(i\). This is the Cartesian Hessian element \(H_{ij}\).

**Relaxed effective coupling:** move something and allow the rest of the crystal to respond. This involves inverses/Schur complements and is intrinsically nonlocal.

---

# But shouldn't the Born–Oppenheimer Hessian still become local in a wide-gap covalent material?

Yes. This is where your density-matrix intuition is valuable.

For a gapped insulator, electronic response is “nearsighted”: localized perturbations produce responses that decay rapidly with distance. Therefore the **short-range electronic contribution to interatomic force constants should also decay strongly**.

This is why real-space phonon force constants are routinely truncated after a finite number of neighbor shells.

For diamond in particular I would expect very strong locality.

For Si it should be somewhat less local because the gap is much smaller, but still reasonably short-ranged.

There is an important exception:

$$
\text{long-range electrostatics}.
$$

In polar crystals, displacement-induced dipoles give force constants behaving asymptotically like

$$
H_{ij}^{\rm dipole}\sim\frac{1}{r_{ij}^3}.
$$

This is the famous long-range dipole-dipole part of lattice dynamics.

Bulk elemental diamond/Si has no ordinary polar optical coupling of that kind. But **your H-terminated surface does** have Si–H bond dipoles and dynamic dipoles. So surface Si–H modes can have a weak but genuinely long-range coupling.

That is particularly relevant to exactly the feature you're interested in: **splitting of a nearly degenerate manifold of Si–H oscillators**.

---

# “Weak” interactions can matter surprisingly much for the 2000 cm⁻¹ fine structure

Suppose two identical Si–H oscillators have stiffness \(k\) and coupling \(J\):

$$
K=
\begin{pmatrix}
k&J\\
J&k
\end{pmatrix}.
$$

Then

$$
\omega_\pm^2=\frac{k\pm J}{m}.
$$

For \(J\ll k\),

$$
\frac{\Delta\omega}{\omega}
\sim
\frac{J}{k}.
$$

More precisely each branch moves by approximately

$$
\frac{\delta\omega}{\omega}\simeq \frac{J}{2k}.
$$

So if

$$
J/k=0.01,
$$

only a **1% force-constant coupling**, the shift is roughly

$$
0.5\%\times2000\ {\rm cm^{-1}}
\sim 10\ {\rm cm^{-1}}.
$$

That is already very visible in FTIR.

Even

$$
J/k\sim10^{-3}
$$

can give order-one cm\(^{-1}\) shifts.

Therefore something which is “100× weaker than a covalent bond” can safely be ignored for structural stability but **not necessarily for fine splitting of nearly resonant Si–H modes**.

That distinction is important.

---

# So the first thing we should determine is whether the DFT Hessian is actually strongly nonlocal

A color plot of a Hessian is misleading because analytic DFT will rarely give exact zeros.

What matters is the magnitude of \(3\times3\) atom blocks

$$
h_{ij}=\|H_{ij}\|_F
$$

versus graph distance and physical distance.

I would make one very simple diagnostic for Si\(_{10}\)H\(_{16}\) and the smallest NC:

$$
h(d)=
\left\{
\|H_{ij}\|_F:
d_{\rm graph}(i,j)=d
\right\}.
$$

Plot distributions for

$$
d=0,1,2,3,4,\ldots
$$

and separately classify

$$
\mathrm{Si-Si,\ Si-H,\ H-H,\ Si\cdots H}.
$$

That immediately tells us whether “dense” means

> many matrix elements of \(10^{-5}\) alongside local elements of \(1\)

or

> genuine long-range blocks of 5–10% of nearest-neighbor strength.

Those are completely different situations.

---

# Even more decisive: truncate the **DFT Hessian itself**

This is the experiment I would trust most.

Take the exact DFT Hessian and construct

$$
H^{(d)}
$$

by retaining only blocks with graph distance

$$
d_{\rm graph}(i,j)\le d.
$$

Test

$$
d=1,2,3,4,\ldots
$$

and diagonalize each one.

There is one important technical condition: after truncation, restore the acoustic sum rule,

$$
H_{ii}
=
-\sum_{j\ne i}H_{ij},
$$

otherwise simply deleting off-diagonal blocks artificially breaks translation invariance.

Then compare

$$
\omega^{(d)}
$$

with the complete DFT spectrum.

This cleanly answers your question without any force-field fitting.

If

$$
H^{(2)}
$$

already reproduces the relevant frequencies, then:

$$
\boxed{\text{a bond+angle local FF has enough spatial range.}}
$$

Its remaining errors come from the *form or parameterization* of those local interactions.

If you need

$$
d=3,
$$

that gives a strong argument for one 1–4 / neighboring-stretch term.

If convergence requires many shells, then a strictly local FF cannot reproduce DFT exactly and some explicit longer-range interaction is justified.

---

# And do this separately for the Si–H band and the complete spectrum

This matters enormously.

The low-frequency cage modes are likely to be much more sensitive to small errors spread through the entire Si network.

Your experimental priority, however, is approximately

$$
1800-2200\ {\rm cm^{-1}}.
$$

For every truncated DFT Hessian I would therefore calculate not only overall RMSD but also

$$
\mathrm{RMSE}_{\rm SiH},
$$

and mode-overlap/projected Si–H DOS.

It is entirely plausible that

$$
d=2
$$

is insufficient for 100 cm\(^{-1}\) cage modes but essentially exact for the Si–H stretches.

Then we shouldn't complicate the FF merely to fix irrelevant low-frequency modes.

---

# There is also a beautiful Wilson version of the same experiment

For the local Si–H stretch basis \(R_H\) we discussed, compute

$$
K_H^{DFT}
=
R_HD^{DFT}R_H^T.
$$

Now look directly at

$$
(K_H)_{ij}
$$

as a function of separation between hydrides.

This answers:

> Is the *Si–H oscillator Hamiltonian* local?

This is actually more directly relevant than asking whether the entire Cartesian nuclear Hessian is local.

You may find something like

$$
K_{ij}\approx
\begin{cases}
20 & \text{same Si}\\
5 & \text{adjacent surface Si}\\
0.8 & \text{second neighbor}\\
0.2 & \text{far}
\end{cases}
$$

in some consistent units.

Then the effective model is obvious.

The Wilson code you already have is well suited to constructing exactly these local stretch directions; the current generic implementation instead spans all bonds+angles and therefore tends toward the whole \(3N-6\) vibrational space. 

---

# One subtle correction about “only compute nearby Hessian elements”

There are two different shortcuts.

One is:

$$
\text{compute only }H_{ij}\text{ for nearby atom pairs}.
$$

The other is:

$$
\text{compute }Hg_a
$$

for selected **local physical displacement directions** \(g_a\).

I now strongly prefer the second.

Suppose \(g_a\) is a local Si–H stretch. Two DFT force calculations,

$$
x_\pm=x_0\pm\delta g_a,
$$

give

$$
Hg_a
\simeq
-\frac{F_+-F_-}{2\delta}.
$$

Although \(g_a\) itself involves only two atoms, the resulting force vector contains responses on **all atoms**.

So you get to directly measure:

> “If I stretch this particular Si–H bond, how local is the DFT response?”

That is almost exactly your physical thought experiment.

If the forces induced 10 Å away are negligible, wonderful: locality is directly demonstrated.

If there is a visible response far away, we quantify it rather than assuming it away.

And for \(m\) selected local coordinates it costs only \(2m\) DFT gradient evaluations. Your earlier discussion already identified this reduction from \(6N\) Cartesian displacements to \(2m\) chosen internal-coordinate displacements. 

---

## So I would run this locality experiment before changing the force field

Concretely:

1. On a small system where we already have the full DFT Hessian, plot \(\|H_{ij}\|\) versus graph/real-space distance.
2. Truncate the **DFT** Hessian at graph ranges \(d=1,2,3,\ldots\), repair the acoustic sum rule, and see how rapidly the full and Si–H-only spectra converge.
3. Build the localized Si–H Wilson matrix \(K_H\) and plot its off-diagonal couplings versus H···H distance and graph separation.
4. Compare those couplings with the same quantities from the present FF.
5. On larger crystals, stop calculating full DFT Hessians and instead perform ± local Si–H / angle displacements and save the full resulting force-response vectors.

My expectation—not something I would assume without testing—is that **diamond-like Si should be predominantly local in the covalent backbone**, and that most of the interesting residual nonlocality for your high-frequency surface band will turn out to be either nearest-neighbor SiH–SiH mechanical coupling or H···H/dynamic-dipole coupling along the surface.

If that is what the data show, then your original physical picture survives almost completely:

$$
\boxed{
\text{local valence springs}
+
\text{one weak surface coupling}
}
$$

can generate very complicated, crystal-shape-dependent global spectra by diagonalization. The dense appearance of the raw DFT Hessian would then mostly be a combination of small electronic-response tails and the fact that “exactly zero” is not the same thing as “physically important.”








---

# USER

On important question: you speak about 3N-6 modes after removing rigid translation and rotation. But I would like to ask about those, how to deal with them diring fitting. Thy can be removed there are rigorous wyas how to remove them from the spectrum. But when doing that we dipically densify the marixes and even originally localized hessian and Wilson localized modes become very non-local. I do not like that. Therefore I think it is better to fit in the basis of orginal cartesian hessian, respo transform it to lozalized willson, but not do this rigid-removal projection. We still want to make this rigid removal protiction when computing respose functon or maybe the spectrum (although having 6 modes at zero should not be problem if we know they come from the rigid translation/rotation we just ignore them), but for fitting we shoupd perhaps rather not do this removal. What you think ? Or are there any other issues we would face?

small notes - you say 1% error due to electrosatics is alrady visible in the FTIR spectrum - maybe depends, but for these nanocrystal the accuracy is no so great and our fit is way to far from 1% accuracy, even DFT is like 10% from experimental resultls. So I do not care so much about his. But I care moe about qualitative features and relative shift resp. relative ordering like bondin (in pahse) lower frequecy than antibonding (counter phase)... somthing like that. I want to fit qualitatively robuts model, not necessarily 1% accurate.


---

# Chat GPT 5.6 sol

Yes — I agree with your instinct here, and I would now make it an explicit principle of the fitting workflow:

> **Do not rigid-mode-project the Hessian used for local force-field fitting.**
>
> Keep the original Cartesian Hessian/local Wilson directions intact. Rigid translations and rotations should disappear because of the physics of the force field, not because we multiply the matrix by a dense global projector.

The current FireCore code does explicitly construct a mass-weighted rigid-body basis and uses \(PDP\) or a reduced basis \(Q^TDQ\) in parts of the spectral machinery.  The reference-mode fitter similarly diagonalizes \(D\) after restricting it to the complete non-rigid subspace.  I think that remains appropriate for **analysis**, but not for the new localized fitting objective.

## The nice thing: Wilson coordinates already remove rigid motion automatically

This is actually an important reason to like Wilson coordinates.

Take a bond coordinate

$$
q_{ij}=|\mathbf r_j-\mathbf r_i|.
$$

Its Wilson vector is

$$
b_{ij}=\frac{\partial q_{ij}}{\partial \mathbf r}.
$$

For a rigid translation \(t\),

$$
\delta\mathbf r_i=\mathbf t
$$

for every atom, and therefore

$$
b_{ij}\cdot t=0.
$$

Likewise for an infinitesimal rigid rotation,

$$
\delta\mathbf r_i
=
\boldsymbol\omega\times(\mathbf r_i-\mathbf R_0),
$$

the bond length does not change to first order:

$$
b_{ij}\cdot r_{\rm rot}=0.
$$

The same is true for an angle Wilson coordinate.

Therefore, ideally,

$$
\boxed{BR_{\rm rigid}=0.}
$$

So if we ask about local bond/angle response using \(B\), **the six rigid modes simply never enter the question**.

There is no need for

$$
P=I-RR^T
$$

and therefore no reason to turn a beautifully local representation into a dense one.

This is a major conceptual simplification.

---

## I would actually go one step further: don't orthogonalize the local Wilson vectors for fitting either

My previous suggestion of Löwdin orthogonalization,

$$
R=(CC^T)^{-1/2}C,
$$

is mathematically elegant for defining a nice orthonormal vibrational tight-binding basis.

But you have identified the same problem again:

$$
(CC^T)^{-1/2}
$$

is generally dense.

So even if every original Si–H stretch \(b_i\) lives on only two atoms, the orthogonalized \(R_i\) can contain tails over many bonds.

That is acceptable for **post-processing / visualization**, but it is contrary to the main philosophy of the fit.

For the fit I would keep the raw local coordinates.

---

# The clean local fitting object is therefore simply \(Hg_a\)

Let \(g_a\) be one local Cartesian displacement corresponding to:

* one Si–H stretch,
* one Si–Si stretch,
* one H–Si–H bend,
* one Si–Si–H bend, etc.

It involves only 2 or 3 atoms.

Then compare

$$
\boxed{
H^{FF}g_a
\quad\text{with}\quad
H^{DFT}g_a
}
$$

without any global projection.

This is extremely physical:

> “I impose this local bond stretch/bend. What force response does the system produce?”

And because

$$
H^{FF}=\sum_p k_p A_p,
$$

we have

$$
H^{FF}g_a
=
\sum_p k_p(A_pg_a),
$$

so the fitting problem is still completely linear.

Even better, you can choose how much of the resulting DFT force vector you want to fit.

For example,

$$
M_aH^{FF}g_a
\simeq
M_aH^{DFT}g_a
$$

where \(M_a\) selects atoms within graph distance 2 or 3 of the locally distorted bond.

That gives you a **strictly local-local fitting problem**:

$$
\boxed{\text{local displacement}\rightarrow\text{local force response}.}
$$

No eigenvectors.
No \(3N-6\) basis.
No dense projector.
No pseudoinverse.

I think this is probably the cleanest formulation of what you originally wanted from Wilson fitting.

---

# Alternatively use pairwise Wilson matrix elements directly

If \(g_a,g_b\) are two local Wilson directions, calculate

$$
K_{ab}=g_a^THg_b.
$$

Again, no projection is necessary.

This has a particularly intuitive interpretation:

* \(K_{aa}\): stiffness when coordinate \(a\) is distorted;
* \(K_{ab}\): direct harmonic coupling between local distortions \(a,b\).

And we only have to compare \(K_{ab}\) for chemically meaningful nearby pairs.

For example:

$$
K_{\rm SiH_i,SiH_i}
$$

for local stretch stiffness;

$$
K_{\rm SiH_1,SiH_2}
$$

for the two hydrides of one SiH\(_2\);

$$
K_{\rm SiH_i,SiH_j}
$$

for hydrides on adjacent surface Si atoms;

$$
K_{\rm SiH,angle}
$$

for stretch–bend coupling.

That becomes almost a **local force-constant table extracted from DFT**, but unlike the pseudoinverse Wilson \(F\), there is no claim that we have uniquely decomposed the entire DFT Hessian into independent springs.

We simply ask for the DFT curvature along two explicitly defined physical directions.

That quantity is unambiguous.

---

# A subtle but important distinction from the problematic old Wilson \(F\)

This is worth emphasizing.

The troublesome quantity was

$$
F=C^{+T}DC^+.
$$

There we asked:

> Find an internal-coordinate force matrix \(F\) whose transformation reproduces the complete DFT Hessian.

Because the internal coordinates are redundant, there are infinitely many answers.

Instead, now calculate

$$
\boxed{K_{ab}=g_a^TDg_b}
$$

directly.

We are asking:

> What curvature does the *actual DFT Hessian* have along these two particular displacement directions?

There is no inversion and therefore no gauge ambiguity.

So I would probably stop calling this quantity the Wilson \(F\) matrix at all.

Wilson gives us the **directions** \(g_a\), but we simply project DFT curvature onto them.

---

# What about normalization?

This is the only nuisance.

A bond coordinate naturally has units of Å, whereas an angle has radians. So if you want to compare absolute entries between stretch and bend coordinates, choose a convention.

For bonds, I like what the current code already does:

$$
q_r=\frac{\Delta r}{r_0}.
$$

Angles can use

$$
q_\theta=\Delta\theta.
$$

The existing code already scales bond Wilson rows by \(1/r_0\) for this reason. 

But for fitting you don't actually have to normalize everything identically. You can simply define a physically meaningful finite perturbation:

$$
\Delta r = 0.01\ {\rm \AA}
$$

or

$$
\Delta\theta=0.01\ {\rm rad},
$$

and normalize the residual by the response magnitude of each coordinate class.

That may actually be more transparent.

---

# Now the rigid modes themselves

For an **exact molecular potential at a stationary geometry**, the Cartesian Hessian should already satisfy six null modes:

$$
HT_x=HT_y=HT_z=0
$$

and

$$
HR_x=HR_y=HR_z=0.
$$

So fundamentally we should not need to remove them.

There are, however, two differences between translation and rotation.

### Translation is especially clean

For a potential depending only on relative coordinates,

$$
E(\{r_i+a\})=E(\{r_i\}),
$$

so

$$
H T=0
$$

exactly at **any geometry**, not just an equilibrium.

Equivalently, the Cartesian Hessian obeys the acoustic sum rule

$$
\boxed{
\sum_j H_{ij}=0.
}
$$

Your local valence FF should satisfy this basically to machine precision if implemented correctly.

This is also an excellent debugging invariant.

### Rotation has a slight subtlety

Rotational invariance means

$$
E(\{Rr_i\})=E(\{r_i\}).
$$

At a stationary geometry,

$$
\nabla E=0,
$$

this gives

$$
HR=0.
$$

But **away from a stationary point**, a rigid rotational path has a second-order geometrical contribution involving the existing forces. So simply applying the linear tangent vector \(\omega\times r_i\) to the Hessian does not necessarily give exactly zero unless the forces vanish.

This matters when fitting all methods on a common DFT geometry that is *not* the FF's equilibrium.

It is another reason I prefer fitting local Wilson responses rather than enforcing six zero eigenvalues during the fit.

---

# This connects to your prestress issue

Your FF sensitivity implementation already knows about this:

$$
A
=
f''(q)\,bb^T
+
f'(q)C.
$$

The second \(f'C\) term is exactly the geometric/prestress part. 

If we set every FF equilibrium bond and angle to the DFT geometry,

$$
r_0=r_{\rm DFT},\qquad
\theta_0=\theta_{\rm DFT},
$$

then

$$
f'=0
$$

term by term, and rotations become clean null directions.

But if we fit transferable \(r_0,\theta_0\), individual terms are prestressed even though total force can cancel.

So I would not use “how exactly zero are the six rigid eigenvalues?” as a primary fitting target.

Use it as a **consistency diagnostic**.

---

# How I would compute the final spectrum

Here too I agree with you: there is no fundamental need to project \(H\) first.

Simply mass-weight:

$$
D=M^{-1/2}HM^{-1/2}
$$

and diagonalize the original matrix.

If everything is good, the spectrum contains

$$
6
$$

very-low/zero eigenvalues followed by the physical modes.

Rather than saying

> “delete the lowest six”

I would identify them by overlap with the known analytic rigid-body vectors.

For each eigenvector \(v_n\),

$$
w_n^{rigid}
=
\sum_{\alpha=1}^{6}
|R_\alpha^Tv_n|^2.
$$

Then modes with \(w_n^{rigid}\approx1\) are rigid.

This is robust even if numerical noise gives them frequencies such as

$$
-3,\; 0.2,\; 1.7,\; 4\ {\rm cm^{-1}}
$$

rather than exact zeros.

And it avoids altering the Hessian.

Only if numerical contamination is so severe that rigid and genuine soft modes mix badly would I diagonalize in the explicitly projected \(3N-6\) subspace. That should be a **fallback analysis operation**, not something done to the stored/reference Hessian.

The current code presently projects \(D\) before diagonalization when `project_acoustic=True`.  I would add the alternative “classify/ignore rigid eigenvectors after ordinary diagonalization” and make that preferable for locality-sensitive work.

---

# Response functions are slightly different

For

$$
Hx=f
$$

there is a genuine mathematical problem: \(H\) has six zero modes, so it cannot be inverted.

There we **must choose a gauge**.

But again we don't have to densify \(H\).

Instead solve the constrained problem

$$
Hx=f,
$$

subject to

$$
R^TM^{1/2}x=0.
$$

Equivalently remove net translation and rotation from the solution.

A Krylov solver can implement projection matrix-free:

$$
v\rightarrow v-R(R^Tv)
$$

without ever constructing the dense

$$
P=I-RR^T.
$$

So even for response calculations, locality/sparsity of \(H\) itself can be preserved.

---

# And I agree with your accuracy correction

My previous “1% coupling can give 10 cm\(^{-1}\)” example was mathematically true but misleading relative to **your actual objective**.

You are not trying to build a spectroscopic force field accurate to 1%.

Your much more important questions are things like:

* does SiH lie below/above SiH\(_2\)?
* do {111}, {100}, edge and corner environments produce the correct ordering?
* does a patch of neighboring SiH groups form a band rather than independent identical oscillators?
* does the in-phase combination go to the correct side of the antisymmetric combination?
* does increasing facet size narrow/shift/converge the band sensibly?
* do defects locally perturb rather than completely rearrange the spectrum?

Those are **topological/qualitative properties of the effective dynamical matrix**.

For that purpose I would willingly tolerate a systematic overall error of 5–10% if the coupling structure is correct.

In fact, this suggests another useful fitting philosophy:

$$
D^{DFT}_{\rm local}
\approx
s\,D^{FF}_{\rm local}
$$

with one global frequency/stiffness scale \(s\), while giving more importance to the **dimensionless coupling ratios**

$$
\boxed{
\eta_{ij}
=
\frac{K_{ij}}
{\sqrt{K_{ii}K_{jj}}}.
}
$$

This is interesting for you because \(\eta_{ij}\) directly measures how strongly two local oscillators are coupled relative to their individual stiffnesses.

The **sign** of

$$
\eta_{ij}
$$

determines whether the symmetric/in-phase or antisymmetric/out-of-phase combination is softer.

The magnitude determines the splitting.

Thus, for your intended “qualitatively robust model,” I would inspect and perhaps partially fit

$$
K_{ii}
$$

for local chemistry,

but especially

$$
\boxed{\eta_{ij}}
$$

for local coupling topology.

That may ultimately be more informative than minimizing another global cm\(^{-1}\) RMSD.

So I think the architecture is becoming quite clean:

$$
\boxed{
\begin{array}{c}
\text{raw sparse Cartesian }H\\[2mm]
\downarrow\\
\text{local Wilson directions }g_a\\[2mm]
\downarrow\\
Hg_a\quad\text{or}\quad K_{ab}=g_a^THg_b\\[2mm]
\downarrow\\
\text{fit local stiffness + coupling ratios}\\[2mm]
\downarrow\\
\text{diagonalize original }H\\
\text{and merely identify/ignore 6 rigid modes}
\end{array}}
$$

No rigid-body projection is needed anywhere in the actual fitting path, and that is probably much closer to the original spirit of what you wanted Wilson coordinates to accomplish.





---

# USER

we are discussijg with coding agent the question of "Frozen Geometry Hessian", the thing is - we must be sure the Forcefield compute the frequencies in its stationdary point resp in its minimum so there are no imaginary frequencies. Can we ensure that in the fit - lake make such params which make sure the focefield have same minimum? But maybe not - especially for many different crystals if our forcefield is not able to reporcue equlibirum geometries from DFT exacly for all possible reference structures. Then it is still better to use relaxed geometry for given forcefield and just fit the elastic parameters around that minimum. This may mean that we canno direcly use the directions of eigenvectros from DFT or even cartesian Hessian, structli speaking if we want to be perfectly rigorous. But since the relaxatons are small (<0.1A) the direction of vectros are still more or lss same after this rectioxtions, so I would not be worid too much. What would you recomand? I would like to pretent the geometry is close engouh we can still use the same direction vectors (Wilson bonds and angles vectros, SVD eigenmodes ...) for both DFT-optimized and FF-optimized geometry and fit just the harmonic expansion around these to different minimum geometries pretentig like they are the same?

---


# CHAt GPT 5.6 sol

Yes. I think there is a clean way to reconcile these requirements, and I would **not require the FF and DFT to have exactly the same minimum**. With only a few transferable parameters that is probably impossible across many different nanocrystals.

The important distinction is:

> **Frozen-geometry Hessian fitting is a parameter-identification tool. The physically predicted spectrum must be evaluated at the FF's own minimum.**

Trying to force those two things to be identical creates unnecessary difficulty.

## 1. First, we can partially enforce the same minimum — and it is cheap

At the DFT minimum \(R_D\),

$$
F^{DFT}(R_D)=0.
$$

Ideally we would also like

$$
F^{FF}(R_D;k)=0.
$$

For fixed equilibrium parameters \(r_0,\theta_0\), the FF force is linear in the stiffnesses:

$$
F^{FF}(R_D;k)=\sum_p k_p f_p(R_D).
$$

So alongside the Hessian-fitting equations we can simply add **force-balance equations**

$$
\boxed{
\sum_p k_p f_p(R_D)\simeq0
}
$$

to the same linear least-squares problem.

Something like

$$
L =
L_{\rm local\,Hessian}
+
\lambda_F
\|F^{FF}(R_D)\|^2
+
L_{\rm reg}.
$$

No nonlinear optimizer is required.

This is very useful because it tells us whether our simple FF **can** have approximately the same equilibrium.

If after fitting we have very small residual forces at all DFT geometries, excellent. If not, this is not a numerical failure: it tells us that the chosen transferable \(r_0,\theta_0\) + 5 stiffnesses simply cannot reproduce all those equilibrium structures simultaneously.

The existing code already distinguishes `local` equilibrium coordinates from transferable `type-average` ones, and explicitly deals with prestress when the individual internal coordinates are not at their own minima.  The underlying theory likewise contains the \(f'(q)C\) prestress term rather than assuming that a relaxed *total* geometry means every individual bond/angle term is individually relaxed. 

I would definitely add the force residual to the fitter.

---

# 2. But I would not demand zero FF force at every DFT minimum

Suppose the DFT optimum is

$$
R_D
$$

and after fitting and relaxation the FF optimum is

$$
R_F=R_D+\Delta R.
$$

If

$$
|\Delta R|\lesssim 0.05-0.1\ \AA
$$

and there is no reconstruction/topology change, these are essentially two representations of the same chemical structure.

The quantities we really want to compare are then

$$
H_D(R_D)
$$

and

$$
H_F(R_F).
$$

Strictly speaking, they are Hessians at different points, so comparing their Cartesian entries element-by-element is not mathematically invariant.

But **local chemically defined coordinates give us a very natural correspondence**.

And I would improve slightly on your proposed “pretend the vectors are identical”:

> Use the **same chemical coordinate labels**, but recompute the Wilson vectors independently at the two minima.

For example,

$$
g_{\mathrm{Si_7-H_{32}}}^{DFT}
$$

is constructed from \(R_D\), while

$$
g_{\mathrm{Si_7-H_{32}}}^{FF}
$$

is constructed from \(R_F\).

Then compare

$$
\boxed{
(g_a^{DFT})^T H_D g_b^{DFT}
\quad\leftrightarrow\quad
(g_a^{FF})^T H_F g_b^{FF}.
}
$$

Same bond \(a\), same neighboring bond \(b\), but each expressed in its own slightly different geometry.

That is actually quite rigorous conceptually. You are comparing:

> curvature of the same chemically defined deformation at each model's own equilibrium.

No pretense that the Cartesian coordinates are exactly identical is needed.

---

# 3. For bonds/angles there is essentially no cost to doing this correctly

Take a Si–H bond whose direction changes by 2° after FF relaxation.

Why use the old DFT vector?

Just recompute

$$
\hat r_{SiH}^{FF}
$$

from the FF-relaxed geometry.

Likewise an H–Si–H angle Wilson vector is trivial to recompute.

Thus I would **not reuse exactly the same Wilson vectors**.

Instead preserve:

```text
coordinate ID:
    SiH stretch (Si=7,H=32)
    HSiH angle (H=31,Si=7,H=32)
    ...
```

and regenerate its geometric vector at each structure.

This completely removes one unnecessary approximation.

---

# 4. You can quantify whether the geometry difference matters

After a rigid Kabsch alignment of \(R_D\) and \(R_F\), calculate for every corresponding local coordinate

$$
s_a=
\left|
\hat g_a^{DFT}\cdot
\hat g_a^{FF}
\right|.
$$

If you find something like

$$
s_a=0.995-1.000
$$

for nearly everything, then you have numerical evidence that your original intuition—

> “the geometries are close enough that the local directions are essentially the same”

—is valid.

This would be a nice diagnostic to output with every fit:

```text
geometry RMSD       0.031 Å
Si-H Wilson overlap mean  0.9992
Si-H Wilson overlap min   0.9961
angle overlap mean        0.9978
```

Then there is no philosophical ambiguity about whether 0.1 Å matters.

---

# 5. Global eigenmodes are different

Here I would be much more cautious.

Suppose

$$
V_D
$$

are DFT normal modes and

$$
V_F
$$

are FF modes at its own minimum.

You can align the two geometries and calculate the overlap matrix

$$
J_{mn}
=
(v_m^{DFT})^T v_n^{FF}.
$$

This is essentially a **Duschinsky-type mode transformation**.

If

$$
J\approx I,
$$

then yes, modes correspond almost one-to-one.

But near degeneracies the individual modes can rotate arbitrarily:

$$
(v_1,v_2)\rightarrow
(v_1\cos\alpha+v_2\sin\alpha,\,
-v_1\sin\alpha+v_2\cos\alpha)
$$

with practically no physical change.

So I would **not make global DFT eigenvectors the primary fitting basis** once the minima differ.

Use them for validation:

* spectral ordering;
* mode-subspace overlap;
* SiH projected weight;
* localization;
* bands of nearly degenerate modes.

But use **local Wilson directions for fitting**.

This fits extremely well with the direction our discussion has taken.

---

# 6. There is an important mathematical price if we insist on fitting at the FF minimum

At fixed geometry,

$$
H^{FF}(R;k)=\sum_p k_p A_p(R)
$$

is linear in \(k\).

That is why our Hessian fitter is so nice.

But if we insist on

$$
R=R_F(k)
$$

where

$$
\nabla E(R_F(k);k)=0,
$$

then

$$
H^{FF}(R_F(k);k)
$$

is **no longer linear in \(k\)** because changing \(k\) changes the equilibrium geometry.

Formally,

$$
\frac{dH}{dk_p}
=
\frac{\partial H}{\partial k_p}
+
\frac{\partial H}{\partial R}
\frac{dR_F}{dk_p}.
$$

And from stationarity,

$$
H
\frac{dR_F}{dk_p}
=
-
\frac{\partial F}{\partial k_p}.
$$

So

$$
\frac{dR_F}{dk_p}
=
-H^+
\frac{\partial F}{\partial k_p}.
$$

Then evaluating \(dH/dk\) requires essentially **third derivatives** through \(\partial H/\partial R\).

That is exactly the complexity we do not want.

So I think it would be a mistake to throw away the linear Hessian fitting framework merely to be formally exact about a 0.03–0.1 Å geometry difference.

---

# 7. Therefore I would use a two-level / iterative scheme

This seems the best compromise.

### Inner fit — cheap and linear

At fixed geometries, solve

$$
k^*
=
\arg\min_k
\left[
L_{\rm local-curvature}
+
\lambda_F L_{\rm force}
+
L_{\rm reg}
\right].
$$

Use local Wilson/Hessian-response quantities rather than global modes.

### Outer validation

With the new \(k^*\):

1. Relax every FF structure:

   $$
   R_F=\arg\min E_{FF}(R;k^*).
   $$

2. Check:

   * RMS geometry difference;
   * residual DFT/FF local-coordinate overlap;
   * absence of non-rigid imaginary FF modes;
   * local curvature matrix at \(R_F\);
   * final vibration spectrum at \(R_F\).

3. If the geometry changed enough to matter, rebuild the FF local bases at \(R_F\) and refit.

Because your geometry changes are small, I would expect perhaps **one or two outer iterations**, not a long nonlinear optimization.

This is effectively self-consistent fitting:

$$
k^{(0)}
\rightarrow
R_F^{(0)}
\rightarrow
k^{(1)}
\rightarrow
R_F^{(1)}
\rightarrow\cdots
$$

while every \(k\)-update remains a simple linear least-squares problem.

---

# 8. There are actually two useful fits, and we should report both

I would deliberately distinguish them in the code/results.

### A. Frozen-reference fit

Evaluate FF and DFT response at

$$
R_D.
$$

Purpose:

> How good is the **force law itself**, with geometry held fixed?

This is extremely useful scientifically.

If it performs badly, the stiffness model is wrong.

### B. Self-relaxed prediction

Evaluate

$$
H_F(R_F).
$$

Purpose:

> How well does the **complete FF prediction** work, including its geometry error?

This is what should generate the reported FF spectrum.

So frozen geometry isn't some embarrassing approximation we have to eliminate. It answers a **different and very useful question**.

I would never interpret imaginary FF frequencies at \(R_D\) as a prediction, because \(R_D\) may not be an FF stationary point.

---

# 9. Ensuring the final FF calculation really is a minimum

This part should be strict.

After FF relaxation, calculate the Hessian at \(R_F\).

A stable isolated molecule/NC should have:

$$
\lambda_i>0
$$

for all internal modes and six near-zero rigid modes.

If a significant negative eigenvalue remains:

$$
\lambda_{\min}<-\epsilon,
$$

then we have found either:

* an insufficiently converged relaxation;
* a saddle;
* or an intrinsically unstable FF parameter set.

I would automate:

```text
relax
→ compute lowest Hessian eigenvalues
→ identify rigid modes
→ if real negative internal mode:
      displace slightly along it
      relax again
```

That way **the production spectrum is always calculated at an actual FF minimum**, regardless of what frozen-geometry fitting does.

---

# 10. The force penalty is particularly attractive for your minimal FF

I think this could become a very informative metric.

For each DFT reference structure report

$$
F_{\rm RMS}^{FF}(R_D).
$$

Maybe you find:

```text
SiH4      0.006 eV/Å
Si2H6     0.014 eV/Å
Si10H16   0.025 eV/Å
small NC  0.041 eV/Å
```

Then your claim

> DFT and FF minima are close enough that their local harmonic expansions can be meaningfully compared

has quantitative support.

Conversely, if one morphology gives

$$
F_{\rm RMS}=0.3\ {\rm eV/\AA},
$$

we immediately know that its geometry difference cannot simply be ignored.

And because these zero-force equations can be appended to the same linear fit, they are almost free computationally.

---

## So my concrete recommendation

I would **not** try to enforce “DFT minimum = FF minimum exactly.”

Instead:

$$
\boxed{
\begin{array}{l}
\textbf{Fit:}\\
\quad\text{DFT local Hessian/Wilson response}\\
\quad+\text{FF force residual at DFT geometry}\\
\quad+\text{regularization}\\[2mm]
\textbf{Predict:}\\
\quad\text{relax FF independently}\\
\quad\text{rebuild local Wilson vectors at FF minimum}\\
\quad\text{compute FF Hessian and spectrum there}\\[2mm]
\textbf{Validate:}\\
\quad\text{geometry RMSD}\\
\quad\text{local Wilson-vector overlap}\\
\quad\text{local curvature/coupling matrix}\\
\quad\text{qualitative spectrum}
\end{array}}
$$

And yes: for the accuracy level you care about, I would absolutely accept comparing the harmonic expansions around the **two nearby minima**. I just would not literally use the same Cartesian vectors when it costs essentially nothing to regenerate the corresponding Si–H/angle Wilson vectors at each geometry.

That gives us most of the rigor without sacrificing the simplicity and excellent conditioning of the linear fitting problem.
