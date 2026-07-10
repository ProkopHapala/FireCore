# Hessian-Based Force-Field Fitting: Theory and Physics

This document provides a rigorous, self-contained derivation of the theory underlying the Hessian-based force-field fitting code in FireCore. The implementation spans:

- **C++ engine**: `cpp/common/molecular/FFfit.h` — analytic sensitivity matrices, linear least-squares, gradient descent
- **Python library**: `pyBall/FFfit.py` — hybrid objective, Wilson/GF projection, regularized bounded solver
- **Test driver**: `tests/tSiNCs/test_FFfit.py` — PySCF reference Hessians, topology, multi-system orchestration

---

## 1. Problem Statement

Given a reference (QM) Hessian $H_{\text{ref}}$ computed at a relaxed geometry, we seek classical force-field stiffness parameters $\{k_t\}$ such that the model Hessian $H_{\text{model}}$ reproduces $H_{\text{ref}}$ as closely as possible.

The force-field energy is a sum of bonded terms:

$$E(\mathbf{r}) = \sum_t k_t \, f_t\!\bigl(q_t(\mathbf{r})\bigr)$$

where:
- $t$ indexes interaction terms (bonds, angles, torsions, 1-4 springs)
- $k_t$ is the stiffness parameter for term $t$
- $f_t$ is the energy function per unit stiffness
- $q_t$ is the internal coordinate (bond length, angle, dihedral)
- $\mathbf{r} = (\mathbf{r}_1, \ldots, \mathbf{r}_N) \in \mathbb{R}^{3N}$ are Cartesian positions

**Key insight**: At fixed geometry, $H$ is *linear* in the stiffness parameters. This makes fitting a linear least-squares problem — no iterative diagonalization is needed for the core fit.

---

## 2. The Chain Rule: Sensitivity Matrices

### 2.1 General Formula

For a single term $E_t = k_t \, f_t(q_t)$, the Cartesian Hessian is:

$$H_t = k_t \, A_t, \qquad A_t \equiv \frac{\partial^2 f_t}{\partial \mathbf{r} \, \partial \mathbf{r}}$$

Applying the chain rule to $f_t(q_t(\mathbf{r}))$:

$$\boxed{A_{t,\alpha\beta} = f_t''(q_t) \, b_{t,\alpha} \, b_{t,\beta} + f_t'(q_t) \, C_{t,\alpha\beta}}$$

where:
- $b_{t,\alpha} = \partial q_t / \partial r_\alpha$ is the **Wilson vector** (gradient of the internal coordinate)
- $C_{t,\alpha\beta} = \partial^2 q_t / \partial r_\alpha \partial r_\beta$ is the **coordinate Hessian** of the internal coordinate
- $f_t'' \, \mathbf{b} \otimes \mathbf{b}^T$ is the **rank-one part** (dominant at equilibrium)
- $f_t' \, C$ is the **prestress part** (nonzero when geometry $\neq$ equilibrium)

### 2.2 Connection to the Wilson GF Method

The Wilson **B matrix** from molecular vibration theory ([Wilson, Decius, Cross, *Molecular Vibrations*, 1955](https://en.wikipedia.org/wiki/GF_method)) collects all Wilson vectors:

$$B_{t,\alpha} = \frac{\partial q_t}{\partial r_\alpha}$$

In the harmonic approximation near equilibrium, internal coordinates relate to Cartesian displacements by:

$$\Delta \mathbf{q} = B \, \Delta \mathbf{r}$$

The **GF method** expresses the dynamical problem in internal coordinates as:

$$G F \mathbf{P} = \lambda \mathbf{P}$$

where $G = B M^{-1} B^T$ is the kinematic (metric) matrix and $F$ is the valence-force constant matrix. Our fitting code uses the Wilson coordinates to define a physically relevant Cartesian subspace. The least-norm projected $F$ is retained for diagnosis, while the fit compares the gauge-invariant dynamical curvature in that subspace rather than fitting redundant entries of $F$ directly.

The implementation in `pyBall/FFfit.py`:

- `build_wilson_matrix(positions, bonds, angles)` — constructs $B$ analytically for bond lengths and angle coordinates
- `internal_hessian_projection(H, B, masses)` — computes $F = C^{+T} D C^+$ where $C = B M^{-1/2}$, giving the least-norm redundant valence force matrix
- `internal_coordinate_basis(B, masses)` — returns the orthonormal basis for the mass-weighted Cartesian space sampled by $B$, handling redundancy via SVD

The least-norm projection is a useful *diagnostic*, but it is not itself a unique
assignment of DFT curvature to individual bonds or angles. The gauge-invariant
row-space curvature used by the hybrid fit is derived in Section 6.4.

### 2.3 At Equilibrium: Rank-One Simplification

When the geometry is at the internal-coordinate minimum ($f_t'(q_0) = 0$) and the energy is normalized so $f_t''(q_0) = 1$:

$$A_t = \mathbf{b}_t \otimes \mathbf{b}_t^T$$

Each scalar internal coordinate contributes a **rank-one** matrix to the Cartesian Hessian. This is the foundation of the Seminario method and its variants.

**Important caveat**: A relaxed Cartesian geometry only implies $\sum_t f_t' \, \partial q_t / \partial \mathbf{r} = 0$ (force balance), *not* $f_t' = 0$ for every individual term. When shared (transferable) equilibrium parameters $r_0$, $\theta_0$ are used across multiple systems, real prestress remains, and the $f' \cdot C$ term must be retained.

---

## 3. Energy Forms and Their Sensitivities

### 3.1 Bond (Harmonic)

$$E = \frac{1}{2} k (r - r_0)^2$$

| Quantity | Expression |
|----------|-----------|
| $f(r)$ | $\frac{1}{2}(r - r_0)^2$ |
| $f'(r)$ | $r - r_0 = \Delta l$ |
| $f''(r)$ | $1$ |
| Wilson vector $\mathbf{b}$ | $\partial r / \partial \mathbf{r}_i = -\mathbf{e}$, $\partial r / \partial \mathbf{r}_j = +\mathbf{e}$ |
| Coord. Hessian $C$ | $C_{ii} = C_{jj} = \mathbf{P}/r$, $C_{ij} = -\mathbf{P}/r$ |

where $\mathbf{e} = (\mathbf{r}_j - \mathbf{r}_i)/r$ is the unit bond vector and $\mathbf{P} = \mathbf{I} - \mathbf{e}\mathbf{e}^T$ is the transverse projector.

**Sensitivity**:

$$A_{\text{bond}} = \mathbf{e} \otimes \mathbf{e}^T + \Delta l \, C_r$$

At equilibrium ($\Delta l = 0$): $A = \mathbf{e} \otimes \mathbf{e}^T$ (rank-one).

The coordinate Hessian $C_r$ has a clear geometric meaning: it accounts for the change in bond direction when atoms move transversely. The $\Delta l \cdot C_r$ prestress term captures the fact that a stretched bond has different transverse curvature than an unstretched one.

**Implementation**: `bond_wilson()`, `bond_coord_hessian()`, `bond_dHdk()` in `FFfit.h`.

### 3.2 Angle (UFF-Normalized Cosine Form)

$$E = \frac{1}{2} k_\theta \frac{(\cos\theta - \cos\theta_0)^2}{1 - \cos^2\theta_0}$$

This is the **primary angle form** used in the fitting code. The normalization by $s_0^2 = 1 - \cos^2\theta_0$ ensures that the fitted parameter $k_\theta$ has direct physical units of eV/rad² — it is the actual local angular curvature at $\theta_0$.

| Quantity | Expression |
|----------|-----------|
| $g(c)$ | $\frac{1}{2}(c - c_0)^2 / s_0^2$, where $c = \cos\theta$, $c_0 = \cos\theta_0$, $s_0^2 = 1 - c_0^2$ |
| $g'(c)$ | $(c - c_0) / s_0^2$ |
| $g''(c)$ | $1 / s_0^2$ |
| Coordinate gradient | $\partial c / \partial \mathbf{r}_i = (\mathbf{v} - c\mathbf{u})/|a|$, $\partial c / \partial \mathbf{r}_k = (\mathbf{u} - c\mathbf{v})/|b|$ |

where $\mathbf{u} = (\mathbf{r}_i - \mathbf{r}_j)/|\mathbf{r}_i - \mathbf{r}_j|$ and $\mathbf{v} = (\mathbf{r}_k - \mathbf{r}_j)/|\mathbf{r}_k - \mathbf{r}_j|$.

**Why $\cos\theta$ instead of $\theta$?** The direct gradient of $c = \cos\theta$ w.r.t. Cartesian positions avoids the artificial $1/\sin\theta$ singularity that appears in $\partial\theta/\partial\mathbf{r}$. For the energy written in $c$, the natural coordinate is $c$ itself, and the sensitivity is well-behaved even near $\theta = 0$ or $\theta = \pi$.

**Sensitivity**:

$$A_{\text{angle}} = \frac{1}{s_0^2} \left( \mathbf{g} \otimes \mathbf{g}^T + (c - c_0) \, C_{\cos} \right)$$

where $\mathbf{g} = \partial(\cos\theta)/\partial\mathbf{r}$ and $C_{\cos} = \partial^2(\cos\theta)/\partial\mathbf{r}^2$.

At equilibrium ($c = c_0$): $A = \mathbf{g} \otimes \mathbf{g} / s_0^2$. Since $\mathbf{g} = -\sin\theta_0 \, \mathbf{b}_\theta$ (where $\mathbf{b}_\theta = \partial\theta/\partial\mathbf{r}$ is the Wilson vector of the angle), we get $A = \mathbf{b}_\theta \otimes \mathbf{b}_\theta$, confirming that $k_\theta$ is the curvature in eV/rad².

**Coordinate Hessian of $\cos\theta$** (9 blocks for atoms $i,j,k$):

$$C_{ii} = \frac{-\mathbf{v}\otimes\mathbf{u} - \mathbf{u}\otimes\mathbf{v} + 3c\,\mathbf{u}\otimes\mathbf{u} - c\,\mathbf{I}}{|a|^2}$$

$$C_{kk} = \frac{-\mathbf{u}\otimes\mathbf{v} - \mathbf{v}\otimes\mathbf{u} + 3c\,\mathbf{v}\otimes\mathbf{v} - c\,\mathbf{I}}{|b|^2}$$

$$C_{ik} = \frac{\mathbf{I} - \mathbf{v}\otimes\mathbf{v} - \mathbf{u}\otimes\mathbf{u} + c\,\mathbf{u}\otimes\mathbf{v}}{|a|\,|b|}$$

The remaining blocks follow from translational invariance ($\partial/\partial\mathbf{r}_j = -\partial/\partial\mathbf{r}_i - \partial/\partial\mathbf{r}_k$):

$$C_{ij} = -C_{ii} - C_{ik}, \quad C_{jk} = -C_{ik} - C_{kk}, \quad C_{jj} = C_{ii} + C_{ik} + C_{ki} + C_{kk}$$

**Bug fix note**: $C_{jk}$ must use $C_{ik}$ (not $C_{ki}$), since $C_{ik} \neq C_{ki}$ in general for the mixed derivative of $\cos\theta$.

**Implementation**: `angle_grad_cos_direct()`, `angle_coord_hessian_cos()`, `angle_dHdk_cos()` in `FFfit.h`.

### 3.3 Angle (UFF Fourier Form)

$$E = K \sum_{n} C_n \cos(n\theta)$$

Unlike the cosine form, this energy is a function of $\theta$ *directly* (not $\cos\theta$). Therefore **no** $\sin^2\theta$ factor appears in the sensitivity:

$$A_{\text{UFF}} = f''(\theta) \, \mathbf{b}_\theta \otimes \mathbf{b}_\theta^T + f'(\theta) \, C_\theta$$

The Fourier coefficients are normalized so $f'(\theta_0) = 0$ and $f''(\theta_0) = 1$, giving $A = \mathbf{b}_\theta \otimes \mathbf{b}_\theta$ at equilibrium.

**Note**: The prestress term $f' \cdot C_\theta$ is currently skipped in the C++ implementation (marked TODO). $C_\theta$ can be derived from $C_{\cos}$ via:

$$C_\theta = -\frac{C_{\cos} + \cos\theta \, \mathbf{b}_\theta \otimes \mathbf{b}_\theta^T}{\sin\theta}$$

**Implementation**: `angle_dHdk_uff()`, `uff_angle_coeffs()` in `FFfit.h`.

### 3.4 Dihedral / Torsion (UFF/Prokop Form)

$$E = V \bigl(1 + d \cos(n\phi)\bigr)$$

The dihedral sensitivity $A_d = \partial H / \partial V$ is computed by finite differences of the analytical gradient in the Python test driver. The key finding (see Section 7) is that this term does not improve the fit for the Si/C nanocrystal systems studied.

---

## 4. Linear Least-Squares Fitting

### 4.1 Normal Equations

Since $H_{\text{model}} = \sum_p k_p A_p$ is linear in $\mathbf{k}$, the weighted least-squares problem:

$$\min_{\mathbf{k}} \left\| \sum_p k_p A_p - H_{\text{ref}} \right\|_W^2$$

has normal equations:

$$G \mathbf{k} = \mathbf{y}, \qquad G_{pq} = \langle A_p, A_q \rangle_W, \quad y_p = \langle A_p, H_{\text{ref}} \rangle_W$$

where $\langle A, B \rangle_W = \sum_\alpha w_\alpha^2 A_\alpha B_\alpha$.

**Implementation**: `FFfit::fit_hessian()` in `FFfit.h` (C++ Gaussian elimination with partial pivoting), `accumulate_normal_equations()` in `pyBall/FFfit.py`.

### 4.2 Gradient Descent (Alternative)

For very large systems where building the full $G$ matrix is expensive, a gradient descent with momentum is available. The gradient of the loss $L = \|H_{\text{model}} - H_{\text{ref}}\|_W^2$ w.r.t. $k_p$ is:

$$\frac{\partial L}{\partial k_p} = 2 \sum_{t \to p} \langle A_t, \Delta H \rangle_W$$

where $\Delta H = H_{\text{model}} - H_{\text{ref}}$ and the sum runs over all terms mapped to parameter $p$. This uses per-term 3×3 blocks only — no full 3N×3N sensitivity matrices are needed.

**Implementation**: `FFfit::fit_gradient_descent()` in `FFfit.h`, `compute_gradient_term_blocks()` in `test_FFfit.py`.

---

## 5. Multi-System Fitting and Parameter Sharing

### 5.1 Transferability Principle

Force-field parameters describe chemical bond *types*, not individual bonds. All Si–Si bonds share one stiffness $k_{\text{SiSi}}$, all H–Si–H angles share one $K_{\text{HSiH}}$. This is the **principle of transferability**.

Without sharing: 152 Si–H bonds → 152 parameters (underdetermined).
With sharing: 152 Si–H bonds → 1 parameter (well-constrained).

### 5.2 Parameter Map

A `ParamMap` assigns each bond/angle term to a free parameter index. Terms with the same index share one stiffness. The mapping is based on:

- **Element types** (default): Si–Si, Si–H, H–Si–H, Si–Si–H, etc.
- **Chemical environment** (Si subtyping): Si atoms are classified by H-coordination as SiH3, SiH2, SiH, or bulk Si. This splits a single "Si" type into four environment-resolved types, so bonds and angles involving surface Si atoms get distinct parameters from bulk Si atoms.
- **Manual assignment**: for custom symmetry groups

**Si type splitting** (`assign_si_environment_types()` in `test_FFfit.py`):
- `SiH3`: Si with ≥3 bonded H (including SiH4)
- `SiH2`: Si with 2 bonded H
- `SiH`: Si with 1 bonded H
- `Si`: bulk Si (no H neighbors)

This is physically motivated: surface Si atoms have different bonding geometry and stiffness from bulk Si atoms. The subtyping can be applied to the central atom only (`--si-subtype-scope central`) or to all three atoms in an angle (`--si-subtype-scope full`).

### 5.3 Multi-System Accumulation

When fitting across multiple systems (e.g. SiH4 + Si nanocrystals of different sizes), the same parameter appears in all systems. Normal equations accumulate:

$$G_{\text{total}} = \sum_{\text{sys}} G_{\text{sys}}, \qquad \mathbf{y}_{\text{total}} = \sum_{\text{sys}} \mathbf{y}_{\text{sys}}$$

This dramatically increases the data-to-parameter ratio and constrains transferable parameters.

### 5.4 Equilibrium Parameters

Two strategies for choosing $r_0$ and $\theta_0$:

- **Local** (`--equilibrium local`): Each term uses its actual relaxed DFT coordinate. Prestress is zero by construction. This is the Hessian-fitting default.
- **Type-average** (`--equilibrium type-average`): All terms of the same type share one mean $r_0 = \langle r \rangle$ and $c_0 = \langle \cos\theta \rangle$. This tests transferability but introduces prestress.

Averaging $c_0 = \langle \cos\theta \rangle$ (not $\cos\langle\theta\rangle$) is important because the angle energy is written in $c = \cos\theta$. The stored $\theta_0 = \arccos(c_0)$ is only for API compatibility.

---

## 6. Hybrid Objective: Mode + Local + Internal

The simple Hessian-element fit asks "which FF Hessian is closest to the QM Hessian element-by-element?" A better question for vibrational spectra is "which FF reproduces the action of the QM dynamical matrix on the vibrational modes that matter?"

The hybrid objective combines three terms, each normalized by its own reference norm:

### 6.1 Mode Objective

For each reference vibrational mode $\mathbf{v}_s$ with eigenvalue $\lambda_s = \omega_s^2$ (from the mass-weighted QM Hessian), the fit targets the Rayleigh quotient:

$$\mathbf{v}_s^T D(\mathbf{k}) \, \mathbf{v}_s = \lambda_s$$

where $D(\mathbf{k}) = M^{-1/2} H(\mathbf{k}) M^{-1/2}$ is the mass-weighted dynamical matrix. Off-diagonal mode mixing $\mathbf{v}_s^T D \, \mathbf{v}_t$ ($s \neq t$) is penalized separately with a mixing weight.

Mode balancing options:
- **equal**: $w_s = 1$ (all modes weighted equally)
- **frequency**: $w_s = 1/\max(|\lambda_s|, \lambda_{\text{floor}})$ (first-order weighting for squared frequency errors)
- **relative**: $w_s = 1/\max(|\lambda_s|, \lambda_{\text{floor}})^2$ (relative eigenvalue accuracy)

The frequency floor prevents low-frequency (translation/rotation-adjacent) modes from dominating.

**Rigid mode removal**: Translations and rotations are removed by SVD of the mass-weighted rigid-body basis (3 translations + 3 rotations for nonlinear molecules, 5 for linear). The remaining modes span the vibrational subspace.

**Implementation**: `reference_vibrational_modes()`, `rigid_and_vibrational_bases()` in `pyBall/FFfit.py`.

### 6.2 Local Hessian Objective

Only mass-weighted Cartesian atom blocks within a configurable bond-graph distance are fitted:

$$L_{\text{local}} = \sum_{\substack{i,j \\ d_{\text{graph}}(i,j) \leq d_{\max}}} \left\| D_{ij}^{\text{FF}} - D_{ij}^{\text{QM}} \right\|_F^2$$

This prevents absent long-range Hessian blocks from dominating the reduced local model. The graph distance is computed by BFS on the bond network (default $d_{\max} = 2$).

**Implementation**: `local_hessian_mask()` in `pyBall/FFfit.py`.

### 6.3 Wilson Decomposition: Diagnostic, Not a Unique Bond-Stiffness Extraction

For diagnostics, form scaled Wilson coordinates. For a bond $a$ of equilibrium
length $r_{0,a}$, use $q_a = \Delta r_a/r_{0,a}$; angles are already expressed
in radians. With the row-scaled Wilson matrix $B_s$ and

$$C = B_s M^{-1/2}, \qquad D=M^{-1/2}HM^{-1/2},$$

the SVD pseudoinverse gives the least-Frobenius-norm representative

$$F_{\rm min}=C^{+T}DC^+.$$

For a single isolated dimensionless bond coordinate,

$$F_{aa}=k_a r_{0,a}^2, \qquad k_a=F_{aa}/r_{0,a}^2.$$

Thus $F$ has units of energy, not dimensionless units. A previous plotting bug
used $r_0 B$ instead of $B/r_0$, then divided by $r_0^2$ again; this suppressed
the displayed bond indicator by $r_0^4$. The corrected implementation uses the
formula above.

#### Redundancy and Gauge Non-Uniqueness

Bond and angle coordinates in a connected tetrahedral network are redundant.
If $N$ satisfies

$$C^T N C = 0,$$

then $F_{\rm min}+N$ reconstructs exactly the same Cartesian dynamical matrix.
Therefore an individual diagonal element $F_{aa}$ depends on coordinate scaling,
the chosen coordinate list, the SVD threshold, and the least-norm gauge. It also
redistributes physical stretch-stretch and stretch-bend coupling into diagonal
entries. It must be called a **projected diagonal indicator**, not a uniquely
defined DFT stiffness of bond $a$.

This explains why raw projected Si-Si indicators can differ strongly between
bulk-like and surface-associated bonds even when a physically constrained fit
finds similar transferable Si-Si springs. The indicator remains useful for
locating missing couplings, but it must not be regularized merely to make its
histogram look uniform.

#### Rigid and Acoustic Modes

Exact translations and rotations satisfy $B R_{\rm rigid}=0$. The Wilson
row space is consequently orthogonal to the rigid-body subspace, so an explicit
projection $D \leftarrow P_{\rm vib}DP_{\rm vib}$ is a useful numerical
sum-rule diagnostic but cannot normally remove a large spread in $F_{aa}$.
Long-wavelength acoustic modes are physical: they are soft because neighboring
atoms move nearly together and $B\mathbf{u}$ is small, not because the local
bond springs must be soft. Removing low-frequency modes or damping $C^+$ to
force uniform diagonals would bias elastic and acoustic physics.

### 6.4 Gauge-Invariant Wilson Row-Space Objective

The fit avoids the gauge ambiguity of $F_{\rm min}$. Let the compact SVD be

$$C=U\Sigma V^T, \qquad Q=V_{[:,1:r]},$$

where $Q$ is an orthonormal basis of the mass-weighted Cartesian subspace
represented by the Wilson coordinates. The internal residual is

$$L_{\rm row}=\left\|Q^T(D_{\rm FF}-D_{\rm QM})Q\right\|_F^2.$$

This is invariant under any nonsingular rescaling or recombination of the
redundant Wilson rows. Only the upper triangle is stored, with off-diagonal
entries multiplied by $\sqrt2$, preserving the exact Frobenius norm while
halving residual storage. For the usual bond-plus-angle set of a nonlinear,
connected molecule, this row space commonly spans all $3N-6$ vibrational
degrees of freedom; it is then complementary to, but not statistically
independent from, the mode objective. The separate weights should therefore be
chosen by validation, not interpreted as independent experimental data.

**Implementation**: `internal_hessian_projection()` remains the diagnostic
projection; `internal_coordinate_basis()` and
`assemble_hybrid_hessian_system()` implement the fitted row-space residual in
`pyBall/FFfit.py`.

### 6.5 Regularization, Hierarchy, and Bounded Solving

All three terms are linear in $\mathbf{k}$, so the combined objective remains a single linear least-squares problem. Rather than forming normal equations (which squares the condition number), the code stacks the direct residual matrices and solves via:

1. **Column scaling**: Each parameter column is normalized by its norm
2. **Tikhonov regularization**: $\rho \sum_p (k_p - k_p^{(0)})^2 / \sigma_p^2$ pulls poorly constrained parameters toward physical priors
3. **Linear coupling rows**: arbitrary physically motivated constraints $R\mathbf{k}\approx\mathbf{t}$ can be appended without forming normal equations
4. **Bounds**: diagonal bond/angle stiffnesses obey $k_p \geq 0$; signed cross terms use symmetric finite bounds, enforced via `scipy.optimize.lsq_linear`
5. **Unobservable parameter detection**: Columns with near-zero norm raise explicit errors

#### Hierarchical Si Environment Regularization

Environment subtypes (`Si`, `SiH`, `SiH2`, `SiH3`) are useful because they can
improve spectra, but bulk `Si-Si` has far fewer observations than surface types.
Fitting all subtype parameters independently therefore risks treating sparse
data as a true chemical difference. For a family $g$ of subtype parameters,
define the fitted family mean

$$\bar{k}_g = \frac{1}{n_g}\sum_{i\in g}k_i.$$

The production hierarchy uses

$$L_{\rm hier}=
\lambda_{\rm dev}\sum_g\sum_{i\in g}\left(\frac{k_i-\bar{k}_g}{s_g}\right)^2
+
\lambda_{\rm mean}\sum_g\left(\frac{\bar{k}_g-k_g^{\rm elem}}{s_g}\right)^2.$$

$k_g^{\rm elem}$ is obtained from the elemental multi-system fit and is used
only as a weak prior on the **family mean**. The deviations remain data-driven.
The variance term is implemented by all pairwise subtype-difference rows; it
does not designate a privileged bulk subtype as the parent. Generic priors are
disabled for subtype members so an unrelated default (for example, 5 eV/A²)
cannot pull an entire family mean away from the DFT-supported value.

For the six-system all-Si fit, this hierarchy gave a stable Si-Si subtype range
of about 9.17--9.79 eV/A², despite least-norm Wilson diagonal indicators ranging
from roughly 2.78 to 7.93 eV/A². The mean frequency RMSE improved from
40.96 to 31.59 cm$^{-1}$ relative to elemental typing, while the mean Hessian
relative Frobenius error remained essentially unchanged (6.86%).

**Implementation**: `subtype_shrinkage_rows()` and
`family_mean_prior_rows()` in `pyBall/FFfit_utils.py`; the driver option is
`--subtype-shrinkage`.

**Implementation**: `solve_regularized_lsq()`, `fit_hybrid_hessian()` in `pyBall/FFfit.py`.

---

## 7. Finding: 1-4 Interactions and Torsions Do Not Help

### 7.1 1-4 Distance Springs

Third-neighbor (1-4) distance springs were added as a simplified dihedral/cross-term: a harmonic restraint on the distance between atoms separated by exactly 3 bonds. The test driver supports these via `--third-bonds` and `build_3rd_neighbor_bonds()`.

**Result**: The comparison framework (`--compare-models`) fits four progressive variants:
1. Bond + angle (baseline)
2. Bond + angle + 1-4
3. Bond + angle + UFF torsion
4. Bond + angle + UFF torsion + 1-4

For the Si and C diamond nanocrystal systems studied, adding 1-4 springs does not meaningfully improve the fit quality (measured by relative Frobenius norm, frequency RMSE, and frequency MAE). The 1-4 distance is too coarse a coordinate to capture the physics of torsional motion, and its sensitivity is highly correlated with the bond and angle sensitivities.

### 7.2 UFF/Prokop Torsions

Proper UFF torsion terms $E = V(1 + d\cos(n\phi))$ were implemented with analytical gradients and finite-difference Hessians (see `Dihedral_Torsion.md` for details).

**Result**: Adding torsions actively *degrades* the fit:

- **Indefinite sensitivity**: $A_d = \partial H / \partial V$ is not positive semidefinite because the $f'(\phi) \cdot C_\phi$ prestress term is non-zero when $\phi$ is not exactly at a potential minimum. For short Si–H/C–H bonds, $C_\phi$ has very large eigenvalues, producing strong negative curvature.
- **Angle collapse**: $A_d$ is highly correlated with the angle sensitivity $A_{\text{angle}}$, so the fit transfers stiffness from angles to torsions — $k_\theta$ collapses to near zero.
- **Phase problem**: For C diamond, many dihedrals have $\phi \approx 120°$, which is a *maximum* of $1 + \cos(3\phi)$ with $d=1$, so the torsion direction is destabilized.
- **Bounded hybrid fit drives $V \to 0$**: Even with correct signed angles and inverse bond lengths, the bounded solver sets torsion amplitudes to zero for Si_R3p8 at $d=1, n=3$.

**Physical interpretation**: For stiff tetrahedral networks (Si, C diamond), the bond and angle terms already capture the dominant Hessian curvature. Torsional motion in these systems is constrained by the geometry — the 3-body angle terms enforce near-tetrahedral angles, and the 4-body torsion is a soft perturbation that the Hessian fit cannot distinguish from noise. This is consistent with the Keating valence-force-field picture, where stretch and bend terms suffice for zone-center phonons, and couplings (stretch–stretch, stretch–bend) are the natural next order — not torsions.

### 7.3 Optional Stretch--Stretch and Stretch--Bend Couplings

The next minimal valence-force-field extension is not a 1-4 spring or a torsion,
but coupling between the two bonds and the angle meeting at a central atom. For
an angle $i$--$j$--$k$, define dimensionless stretches

$$q_1=\Delta r_{ij}/r_{0,ij}, \qquad q_2=\Delta r_{jk}/r_{0,jk},$$

and the angular displacement $\Delta\theta$. The optional symmetry-respecting
terms are

$$E_{rr}=K_{rr}q_1q_2,$$

$$E_{r\theta}=K_{r\theta}\frac{q_1+q_2}{\sqrt2}\Delta\theta.$$

The $1/\sqrt2$ normalizes the symmetric stretch combination. Both coefficients
have units of eV and are signed: a cross coefficient is not an independent
positive stiffness, so constraining it non-negative would be unphysical.

At a locally relaxed geometry ($q_1=q_2=\Delta\theta=0$), their Hessian
sensitivities are exact analytic outer products:

$$A_{rr}=\mathbf{g}_1\mathbf{g}_2^T+\mathbf{g}_2\mathbf{g}_1^T,$$

$$A_{r\theta}=\mathbf{g}_s\mathbf{g}_\theta^T+
\mathbf{g}_\theta\mathbf{g}_s^T, \qquad
\mathbf{g}_s=(\mathbf{g}_1+\mathbf{g}_2)/\sqrt2,$$

where $\mathbf{g}_{1,2}=\partial q_{1,2}/\partial\mathbf{r}$ and
$\mathbf{g}_\theta=\partial\theta/\partial\mathbf{r}$. These matrices are
symmetric and annihilate rigid Cartesian displacements. They are only enabled
for `--equilibrium local`: using type-averaged $r_0,\theta_0$ would require
the additional prestress derivatives of the cross energy, which are deliberately
not approximated silently.

Cross terms are typed by the same chemical key as their parent angle, receive
zero-centered regularization, and use finite symmetric bounds. They are disabled
unless `--stretch-stretch` and/or `--stretch-bend` are requested.

For the all-Si hierarchical fit, stretch--stretch alone changed little. Adding
both cross terms reduced the mean Hessian relative Frobenius error from 6.86%
to 5.86%, with condition number 1.67. The network Si--Si--Si coefficients were
small ($K_{rr}\approx0.012$ eV, $K_{r\theta}\approx4\times10^{-4}$ eV),
whereas the largest couplings were H--Si--H stretch--bend terms
(about 0.43--0.52 eV). Thus the principal extra curvature captured by this
minimal basis is Si-H valence coupling, not a hidden soft Si-Si torsion.

**Implementation**: `build_cross_param_maps()` and
`compute_cross_sensitivity()` in `pyBall/FFfit_utils.py`.

---

## 8. Mass Weighting and Dynamical Matrix

The mass-weighted dynamical matrix is:

$$D = M^{-1/2} H M^{-1/2}$$

where $M = \text{diag}(m_1, m_1, m_1, m_2, m_2, m_2, \ldots)$. Its eigenvalues are $\lambda_s = \omega_s^2$, with eigenvectors $\mathbf{v}_s$ in mass-weighted coordinates.

The conversion between Cartesian and mass-weighted representations is:

$$\mathbf{v}_s = M^{1/2} \mathbf{u}_s, \qquad \mathbf{u}_s = M^{-1/2} \mathbf{v}_s$$

where $\mathbf{u}_s$ are Cartesian displacement vectors.

Frequency in cm⁻¹: $\nu_s = \sqrt{|\lambda_s|} \times 521.5$ (conversion factor: $\sqrt{\text{eV}/(\text{Å}^2 \cdot \text{amu})} \to \text{cm}^{-1}$).

**Implementation**: `mass_weight_hessian()`, `get_frequencies_cm1()`, `get_reference_modes_and_freqs()` in `pyBall/FFfit.py` / `test_FFfit.py`.

---

## 9. Validation Diagnostics

The fitting code reports several complementary quality metrics:

- **Relative Frobenius norm**: $\|H_{\text{model}} - H_{\text{ref}}\| / \|H_{\text{ref}}\| \times 100\%$ — global Hessian accuracy
- **Frequency RMSE/MAE**: One-to-one minimum-cost assignment between reference and model frequencies (via Hungarian algorithm)
- **Negative modes count**: Number of model eigenvalues below $-(10/521.5)^2$ — indicates structural instability
- **Condition number**: $\sigma_{\max} / \sigma_{\min}$ of the stacked residual matrix — warns of ill-conditioning
- **Unobservable parameters**: Columns with near-zero norm in the residual matrix — parameters that the data cannot constrain
- **Rigid invariance**: every bond, angle, and analytic cross sensitivity must annihilate the six rigid Cartesian displacements
- **Wilson rank and scaling**: the SVD rank and bond-coordinate $1/r_0$ scaling are reported/checked; diagonal Wilson indicators are never treated as fitted parameters
- **Cross-term bounds**: signed cross coefficients are inspected together with the model negative-mode count; a lower residual must not be accepted as a physical improvement if it creates an unstable Hessian

---

## 10. Summary of Key Design Decisions

| Decision | Rationale |
|----------|-----------|
| $\cos\theta$ as angle coordinate | Avoids $1/\sin\theta$ singularity; direct Cartesian gradient |
| Normalized cosine form $E = \frac{1}{2}k_\theta (c-c_0)^2/s_0^2$ | $k_\theta$ has direct units of eV/rad² |
| Prestress term $f' \cdot C$ retained | Needed for transferable (type-averaged) equilibrium parameters |
| Hybrid mode + local + internal objective | Targets frequencies directly while preserving local Hessian structure |
| Direct stacked solve (not normal equations) | Avoids squaring the condition number |
| Non-negative bounds on stiffnesses | Physical constraint; prevents unphysical negative curvatures |
| Hierarchical Si environment subtyping | Allows surface shifts while shrinking sparse subtype data toward a data-supported family mean |
| Wilson least-norm $F$ used only diagnostically | Individual diagonal entries are gauge-dependent in redundant coordinates |
| Optional stretch--stretch/stretch--bend terms | Minimal signed valence couplings; local-equilibrium-only, zero-centered, bounded |
| 1-4 and torsions excluded | Do not improve fit for tetrahedral Si/C nanocrystals; torsion sensitivity is indefinite |

---

## References

1. E. B. Wilson, J. C. Decius, P. C. Cross, *Molecular Vibrations: The Theory of Infrared and Raman Vibrational Spectra* (McGraw-Hill, 1955). — Original GF method.
2. J. M. Seminario, "Calculation of intramolecular force fields from second-derivative matrices," *Int. J. Quantum Chem.* **60**, 1271 (1996). — Direct Hessian-block extraction.
3. R. M. Martin, "Elastic strain energy and invariance requirements," *Phys. Rev.* **145**, 637 (1966). — Keating VFF for diamond/zinc-blende.
4. [GF method — Wikipedia](https://en.wikipedia.org/wiki/GF_method)
