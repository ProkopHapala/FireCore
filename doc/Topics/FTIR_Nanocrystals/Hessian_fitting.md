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