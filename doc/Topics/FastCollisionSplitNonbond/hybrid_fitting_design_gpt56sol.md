---
type: design
title: Hybrid B-Spline + Finite-Support Atomic Potential Fitting
tags: [potential-fitting, b-spline, atomic-basis, ewald, boltzmann-weighting, opencl]
description: Design and implementation of a linear hybrid potential fitting approach using tricubic B-spline + compact-support radial atomic basis with erf-damped Coulomb reference and Boltzmann-weighted least squares.
---

# Hybrid B-Spline + Finite-Support Atomic Potential Fitting

## 1. Goal

Approximate a periodic reference interaction potential

\[
V_{\mathrm{ref}}(\mathbf r)
\]

by a hybrid model

\[
V_{\mathrm{model}}(\mathbf r)
=
V_{\mathrm{BS}}(\mathbf r)
+
V_{\mathrm{atom}}(\mathbf r),
\]

where

- \(V_{\mathrm{BS}}\) is a smooth tricubic B-spline field on a coarse regular 3D grid;
- \(V_{\mathrm{atom}}\) is a sum of compact atom-centred radial functions evaluated only near selected substrate atoms.

The intended runtime decomposition is:

```text
smooth long-range / low-spatial-frequency field
    -> tricubic B-spline grid

localized short-range correction
    -> finite-support radial atomic functions
    -> accelerated by one-level AABB or uniform-cell lookup
```

All fitted parameters must enter **linearly**. There are no fitted Morse decay constants, equilibrium radii, cutoffs, or other nonlinear shape parameters.

The model must reproduce both the reference potential and its derivatives in the physically relevant region while remaining cheap and conservative at runtime.

---

## 2. Reference potential

The reference potential is computed using the existing tested periodic machinery with **erf-damped Coulomb** to avoid the 1/r singularity at short range.

For a given moving probe-atom type \(a\),

\[
V_{\mathrm{ref}}^{(a)}(\mathbf r)
=
V_{\mathrm{Morse}}^{(a)}(\mathbf r)
+
V_{\mathrm{Coul,erf}}^{(a)}(\mathbf r;\sigma),
\]

where the damped Coulomb is

\[
V_{\mathrm{Coul,erf}}(r;\sigma) = \frac{k_e\, q_j\, \mathrm{erf}(r/\sigma)}{r}.
\]

At \(r \gg \sigma\): \(\mathrm{erf}(r/\sigma) \to 1\), recovering full Coulomb.
At \(r \ll \sigma\): \(\mathrm{erf}(r/\sigma)/r \to 2/(\sigma\sqrt{\pi})\), finite and smooth.

This damping models electron screening at short range and prevents the Coulomb singularity from dominating the fit, which would otherwise bury the Morse repulsion (Morse \(D_e \sim 0.001\)–0.01 eV vs bare Coulomb \(\sim 10\) eV at 1.2 Å).

**Implementation**: The Ewald2D routine computes bare \(1/r\) Coulomb. The correction \(\Delta V = -k_e q\,\mathrm{erfc}(r/\sigma)/r\) is subtracted per-atom with analytic force:

\[
\frac{d}{dr}\left[-\frac{k_e q\,\mathrm{erfc}(r/\sigma)}{r}\right]
= k_e q\left[\frac{\mathrm{erfc}(r/\sigma)}{r^2} + \frac{2}{\sigma\sqrt{\pi}}\frac{e^{-r^2/\sigma^2}}{r}\right].
\]

The damping parameter \(\sigma\) (default 2.0 Å) is exposed as a fitting input. Smaller \(\sigma\) = stronger damping = Morse more visible.

Requirements:

1. The lateral directions \(x,y\) use the physical periodic unit cell.
2. Electrostatics use the existing 2D-periodic Ewald implementation, including the proper exponential dependence along \(z\).
3. Short-range terms include all required periodic images consistently with their physical cutoff.
4. The same electrostatic gauge is used for reference generation, fitting, and validation.
5. A convenient gauge is
   \[
   \langle V_{\mathrm{ref}}(x,y,z_{\mathrm{top}})\rangle_{x,y}=0.
   \]

The fitting code must not replace or approximate the already tested Ewald machinery. It consumes its output as the source of truth.

If several moving atom types are required, fit a separate hybrid field \(V_{\mathrm{model}}^{(a)}\) for each type, unless a later factorized representation is introduced explicitly.

---

## 3. Model definition

### 3.1 Smooth B-spline component

\[
V_{\mathrm{BS}}(\mathbf r)
=
\sum_k c_k B_k(\mathbf r),
\]

where

- \(B_k\) are tensor-product tricubic B-spline basis functions;
- \(c_k\) are fitted control coefficients;
- the control-grid spacing is fixed by design, typically
  \[
  \Delta_{\mathrm{BS}}\approx 1.0\ \text{\AA}.
  \]

The coarse spacing and cubic basis enforce smoothness by construction. No Laplacian smoothness penalty is required in the initial implementation.

Boundary conditions:

- periodic in \(x,y\);
- explicitly specified in \(z\), preferably clamped or fixed-gauge at the upper boundary;
- the coefficient indexing wraps periodically in \(x,y\).

Each sample point depends on at most \(4^3=64\) B-spline coefficients.

### 3.2 Finite-support atomic component

\[
V_{\mathrm{atom}}(\mathbf r)
=
\sum_t \sum_m p_{tm}
\sum_{j\in t}
\phi_m(|\mathbf r-\mathbf R_j|;R_{c,t}),
\]

where

- \(t\) labels substrate atom types or symmetry-equivalent site classes;
- \(p_{tm}\) are fitted linear coefficients shared by all equivalent atoms;
- \(R_{c,t}\) is fixed by design, not fitted;
- typical values are \(R_c=4\)--\(6\) Å, chosen together with the runtime AABB or cell size.

The initial radial basis is

\[
\phi_n(r;R_c)=
\begin{cases}
\left(1-\dfrac{r}{R_c}\right)^{2n}, & 0\le r<R_c,\\[2mm]
0, & r\ge R_c.
\end{cases}
\]

Recommended default:

\[
n=2,3,4,\ldots,N_n.
\]

Using \(n\ge2\) gives

\[
\phi_n(R_c)=0,\qquad
\phi_n'(R_c)=0,\qquad
\phi_n''(R_c)=0,
\]

so the atomic potential, force, and radial curvature join continuously to zero at the cutoff.

Its radial derivative is

\[
\frac{d\phi_n}{dr}
=
-\frac{2n}{R_c}
\left(1-\frac{r}{R_c}\right)^{2n-1},
\qquad r<R_c.
\]

The Cartesian gradient is

\[
\nabla\phi_n
=
\frac{d\phi_n}{dr}
\frac{\mathbf r-\mathbf R_j}{|\mathbf r-\mathbf R_j|}.
\]

The singular direction at \(r=0\) is irrelevant because a nuclear core region is excluded from fitting and handled separately at runtime.

If the simple sequence \((1-r/R_c)^{2n}\) is not expressive enough, enrich it without introducing nonlinear parameters:

\[
\phi_{mn}(r)
=
\left(\frac{r}{R_c}\right)^m
\left(1-\frac{r}{R_c}\right)^q,
\qquad q\ge3.
\]

This remains a fixed linear compact basis.

### 3.3 Periodic images of atomic functions

Atomic functions must obey the same lateral periodicity as the reference field:

\[
V_{\mathrm{atom}}(\mathbf r)
=
\sum_{j,\mathbf T_{xy}}
f_{t(j)}(|\mathbf r-\mathbf R_j-\mathbf T_{xy}|),
\]

including every lateral image whose support sphere intersects the evaluation point.

A minimum-image convention is insufficient when

\[
R_c > \frac12\min(L_x,L_y).
\]

For a small unit cell such as NaCl and \(R_c=4\)--\(6\) Å, several images of the same basis atom can contribute. Image translations must therefore be enumerated explicitly from \(R_c,L_x,L_y\).

---

## 4. Sampling domain

### 4.1 Oversampling relative to the B-spline grid

The fitting samples must not coincide only with B-spline nodes.

For

\[
\Delta_{\mathrm{BS}}\approx1.0\ \text{\AA},
\]

use a fitting spacing such as

\[
\Delta_{\mathrm{sample}}=0.5\ \text{\AA}
\]

or preferably

\[
\Delta_{\mathrm{sample}}=0.25\ \text{\AA}
\]

for the final fit.

This gives two to four samples per B-spline interval in each direction and constrains interpolation between control nodes.

The lateral sampling covers one complete periodic unit cell. The \(z\) interval covers the physically relevant region from the near-surface repulsive region to the upper vacuum boundary.

### 4.2 Excluded nuclear cores

Exclude every sample satisfying

\[
\min_j |\mathbf r_i-\mathbf R_j-\mathbf T_{xy}| < R_{\mathrm{core},t(j)}.
\]

Typical values are

\[
R_{\mathrm{core}}\approx1.0\text{--}1.5\ \text{\AA}.
\]

This avoids fitting:

- the Coulomb singularity;
- physically inaccessible configurations;
- regions intended for the separate hard-core collision continuation.

The fit is not required to define trustworthy behaviour inside the excluded cores.

At runtime, \(r<R_{\mathrm{core}}\) is handled by a separately constructed hard-core potential, preferably a quadratic continuation matched in value and derivative at \(R_{\mathrm{core}}\).

### 4.3 Training and validation grids

Use distinct interlaced grids:

- training grid: spacing \(0.5\) or \(0.25\) Å;
- validation grid: shifted by half a training-grid step, or a denser independent grid.

This detects errors hidden by sampling only at a regular set aligned with the spline lattice.

---

## 5. Fit potential values and forces

For fixed basis functions the model is linear:

\[
V_{\mathrm{model}}=A_V x,
\]

where

\[
x=
\begin{bmatrix}
c\\p
\end{bmatrix}.
\]

The force is also linear:

\[
\mathbf F_{\mathrm{model}}
=
-\nabla V_{\mathrm{model}}
=
A_F x.
\]

### 5.1 Boltzmann weighting

The physically relevant region for adsorption is the **attractive well** (\(V < 0\)), not the repulsive wall (\(V > 0\)) or the far field (\(V \approx 0\)). An unweighted LSQ is dominated by the large-magnitude repulsive region and the easy far field, leaving the attractive well poorly fit.

Apply **Boltzmann sample weights** to emphasize the attractive region:

\[
w_i = \min\!\left(\exp\!\left(-\beta\, V_{\mathrm{ref}}(\mathbf r_i)\right),\; w_{\max}\right),
\qquad
\beta = \frac{1}{k_B T_{\mathrm{eff}}},
\]

where \(T_{\mathrm{eff}}\) is an effective temperature (default 300 K, \(\beta \approx 38.9\) eV\(^{-1}\)) and \(w_{\max}=100\) prevents any single point from dominating.

This gives:
- **Attractive** (\(V < 0\)): \(w = e^{+\beta|V|} \gg 1\) → high weight (emphasized)
- **Far field** (\(V \approx 0\)): \(w \approx 1\) → moderate weight
- **Repulsive** (\(V > 0\)): \(w = e^{-\beta V} \ll 1\) → low weight (suppressed)

The weight is applied as \(\sqrt{w_i}\) to each row of the stacked LSQ system (both value and force rows).

### 5.2 Stacked least-squares problem

\[
\min_x
\left\|
\begin{bmatrix}
\sqrt{W}\, A_V\\
\sqrt{\eta_F}\,\sqrt{W}\, A_F
\end{bmatrix}
x
-
\begin{bmatrix}
\sqrt{W}\, V_{\mathrm{ref}}\\
\sqrt{\eta_F}\,\sqrt{W}\, \mathbf F_{\mathrm{ref}}
\end{bmatrix}
\right\|^2,
\]

where \(W = \mathrm{diag}(w_i)\) and \(\eta_F\) is the force weight (default 0.25).

The force block contains the three Cartesian components:

\[
A_F=
\begin{bmatrix}
-\partial_x A_V\\
-\partial_y A_V\\
-\partial_z A_V
\end{bmatrix}.
\]

All B-spline and atomic derivatives are evaluated analytically from the same basis used for the potential. Potential and force are therefore exactly conservative with respect to the fitted model.

A practical initial sequence is:

1. fit potential values only;
2. verify the basis construction;
3. add force rows;
4. add Boltzmann weighting;
5. tune \(\eta_F\) and \(T_{\mathrm{eff}}\) from physically relevant validation scans.

---

## 6. Defining the split between smooth and atomic components

A coarse B-spline basis strongly limits short-wavelength content, and finite support limits the atomic component spatially. However, broad atomic functions with \(R_c=4\)--\(6\) Å can still overlap with the low-frequency spline space.

The total model can fit well while the decomposition is poorly determined.

Rather than adding empirical smoothness penalties, define the split by construction using weighted orthogonalization.

Let

\[
B
\]

be the stacked value-and-force B-spline operator and

\[
\Phi
\]

the stacked atomic operator.

For every atomic basis column \(\phi_m\), solve

\[
q_m
=
\arg\min_q
\|W(Bq-\phi_m)\|^2.
\]

Define

\[
\widetilde\phi_m
=
\phi_m-Bq_m.
\]

Then

\[
\widetilde\Phi=\Phi-BQ
\]

contains only the part of each atomic basis column that cannot be represented by the chosen coarse spline basis under the fitting metric.

Fit

\[
Bx+\widetilde\Phi p.
\]

After fitting, convert back to the desired runtime representation:

\[
Bc+\widetilde\Phi p
=
B(c-Qp)+\Phi p.
\]

Therefore store

\[
c_{\mathrm{runtime}}=c-Qp
\]

together with the original compact radial coefficients \(p\).

Advantages:

- no nonlinear fitting;
- no smoothness penalty;
- a well-defined division of labour between spline and atomic basis;
- the runtime atomic functions remain exactly radial and compact;
- the spline absorbs the low-frequency component of the atomic basis.

This orthogonalization is recommended if simultaneous fitting shows large cancelling coefficients, poor conditioning, or unstable component separation. The first prototype may fit the raw combined basis and compare both approaches.

---

## 7. Symmetry and parameter sharing

Fit atomic coefficients by atom type or symmetry-equivalent site class, not independently for every atom:

\[
p_{jm}=p_{t(j)m}.
\]

This:

- preserves crystal symmetry;
- reduces parameter count;
- improves conditioning;
- avoids arbitrary compensation between neighbouring equivalent atoms;
- makes the local functions reusable.

For defects or deliberately inequivalent sites, introduce separate site classes explicitly.

The B-spline field is periodic by construction. Additional crystal symmetries can optionally be enforced by tying symmetry-related coefficients or by augmenting the sample set with symmetry-equivalent points.

---

## 8. Linear solver and matrix representation

### 8.1 Reference implementation

For one small unit cell, the number of parameters is likely modest. Begin with an explicit design matrix and a reliable dense or sparse least-squares solver:

```python
scipy.linalg.lstsq
```

for small dense systems, or

```python
scipy.sparse.linalg.lsqr
scipy.sparse.linalg.lsmr
```

for sparse systems.

This reference implementation should prioritize transparency and diagnostics over ultimate memory efficiency.

### 8.2 Column scaling

Before solving, scale each parameter column to comparable weighted norm:

\[
\widehat A_j=\frac{A_j}{s_j},
\qquad
s_j=\max(\|A_j\|,\epsilon).
\]

Solve for the scaled coefficients and transform back afterwards.

Column scaling is not regularization. It improves numerical conditioning when combining:

- potential and force rows;
- spline and atomic columns;
- radial basis functions with different powers.

### 8.3 Matrix-free extension

For larger grids, avoid explicitly storing all \(64N_{\mathrm{sample}}\) B-spline entries.

Use a `scipy.sparse.linalg.LinearOperator` implementing:

```python
matvec(x)   -> values and forces at all samples
rmatvec(y)  -> adjoint accumulation into spline and atomic coefficients
```

The spline operator can exploit its tensor-product structure; the atomic operator can exploit finite support.

This is an optimization phase, not a prerequisite for validating the method.

---

## 9. Runtime representation

The fitted model is evaluated as

\[
V_{\mathrm{model}}(\mathbf r)
=
V_{\mathrm{BS}}(\mathbf r)
+
\sum_{j\in\mathrm{near}}
f_{t(j)}(|\mathbf r-\mathbf R_j|).
\]

### B-spline part

- tricubic interpolation;
- analytic derivative of the same interpolant;
- periodic wrapping in \(x,y\);
- no independent force grid.

### Atomic part

- finite support \(R_c\);
- one-level AABB or uniform-cell/PIC lookup;
- no atom outside \(R_c\) is evaluated;
- fixed short loops over a small number of radial basis powers;
- analytic radial derivative.

Choose \(R_c\) jointly with the runtime spatial partition:

- \(R_c\approx4\) Å for smaller cells and stricter locality;
- \(R_c\approx6\) Å if a broader local correction materially improves accuracy;
- the spatial cell size should be \(R_c\) or a simple fraction of it.

The same basis and coefficients used in fitting are used at runtime. No refitting or parameter transformation occurs except the optional spline-coefficient adjustment after orthogonalization.

---

## 10. Validation

Report errors separately on training and independent validation grids.

### Field metrics

\[
\mathrm{RMSE}_V,\qquad
\max|V_{\mathrm{model}}-V_{\mathrm{ref}}|,
\]

\[
\mathrm{RMSE}_{F_x},
\quad
\mathrm{RMSE}_{F_y},
\quad
\mathrm{RMSE}_{F_z},
\]

\[
\max|\mathbf F_{\mathrm{model}}-\mathbf F_{\mathrm{ref}}|.
\]

### Physically relevant scans

Compare vertical and lateral scans through:

- top sites;
- bridge sites;
- hollow sites;
- positions close to the hard-core boundary;
- typical adsorption heights.

Measure errors in:

- equilibrium height;
- well depth;
- lateral corrugation/barriers;
- force at the hard-core matching radius.

### Decomposition diagnostics

Plot separately:

- \(V_{\mathrm{ref}}\);
- \(V_{\mathrm{model}}\);
- \(V_{\mathrm{BS}}\);
- \(V_{\mathrm{atom}}\);
- residual \(V_{\mathrm{ref}}-V_{\mathrm{model}}\);
- corresponding force residuals.

Also report:

- column norms;
- solver residual;
- estimated condition number when available;
- maximum fitted atomic coefficient;
- cancellation ratio
  \[
  \frac{|V_{\mathrm{BS}}|+|V_{\mathrm{atom}}|}
       {|V_{\mathrm{model}}|+\epsilon}.
  \]

Large cancellation indicates an ill-defined split and motivates atomic-column orthogonalization.

---

## 11. Implementation plan

### Phase 1: reference B-spline fit

1. Generate \(V_{\mathrm{ref}}\) and \(\mathbf F_{\mathrm{ref}}\) using the existing 2D-Ewald and Morse machinery.
2. Build a periodic tricubic B-spline basis with \(\Delta_{\mathrm{BS}}\approx1\) Å.
3. Sample at \(0.5\) Å or finer, not only at spline nodes.
4. Exclude all nuclear cores.
5. Fit B-spline coefficients only.
6. Validate on a shifted grid.

### Phase 2: finite-support atomic basis

1. Add shared atom-type basis functions
   \[
   \phi_n=(1-r/R_c)^{2n},
   \qquad n\ge2.
   \]
2. Explicitly sum all periodic images within \(R_c\).
3. Fit B-spline and atomic coefficients simultaneously.
4. Inspect conditioning and component cancellation.

### Phase 3: force-aware fit

1. Add analytic B-spline and radial derivative rows.
2. Fit values and forces jointly.
3. Tune the value/force weighting using physical scans.

### Phase 4: define the decomposition robustly

1. Orthogonalize atomic columns against the spline space.
2. Refit using \(\widetilde\Phi\).
3. Convert to runtime coefficients
   \[
   c_{\mathrm{runtime}}=c-Qp.
   \]
4. Compare accuracy and conditioning with the raw simultaneous fit.

### Phase 5: runtime implementation

1. Export spline coefficients and atomic type coefficients.
2. Evaluate the spline and its derivative on GPU.
3. Evaluate nearby finite-support atomic functions using one-level cell/AABB lookup.
4. Add the separate matched hard-core continuation for excluded nuclear regions.
5. Compare runtime output directly against the Python reference evaluator.

---

## 12. Implemented code structure

```text
pyBall/OCL/Surface_utils.py
    _hybrid_cubic(u, i, nx, periodic)        # tricubic B-spline basis eval at single point
    _hybrid_bspline(points, cell, zr, step)  # B-spline basis + gradients for all sample points
    _hybrid_atom_basis(points, apos, types, cell, Rc)  # compact radial basis + gradients
    _hybrid_reference(points, apos, qs, reqs, cell, sigma)  # erf-damped Ewald + Morse reference
    _hybrid_grid(cell, zr, h, shift)         # generate 3D sampling grid
    test_hybrid_potential_nacl(...)          # full fit + validation + plotting

    # Legacy (kept for reference, not used by current linear fit):
    eval_folded_morse(...)                   # windowed Morse with compact support
    eval_hardcore_quadratic_v2(...)          # quadratic + windowed Morse with derivative matching
    fit_bspline_3d(...)                      # Gaussian-smoothed B-spline fit (old approach)

tests/tMMFF/test_hybrid_potential.py
    # Thin caller: imports and invokes test_hybrid_potential_nacl()

cpp/common_resources/cl/Surface.cl
    eval_hardcore_grid  # OpenCL kernel: per-atom compact radial basis on a grid
    eval_hybrid_grid    # OpenCL kernel: tricubic B-spline + hardcore atomic basis
```

Actual main interface:

```python
test_hybrid_potential_nacl(
    save_dir='results_hybrid',
    Ng=80,             # grid resolution for 2D plots
    R_cut_hc=4.,       # atomic basis cutoff (Å)
    spline_step=1.,    # B-spline grid spacing (Å)
    sample_step=.5,    # training sample spacing (Å)
    force_weight=.25,  # η_F: force vs potential weight
    T_eff=300.,        # Boltzmann weighting temperature (K)
    sigma=2.0,         # erf-damping width for Coulomb (Å)
)
```

Returns: `{'coef', 'validation', 'plots', 'ewald'}`.

---

## 13. Design principles

1. **Purely linear fitting.**
2. **Smoothness and finite support enforced by basis construction.**
3. **Reference field generated by the existing proper 2D-periodic Ewald implementation with erf-damped Coulomb.**
4. **Oversampled fitting between B-spline nodes.**
5. **Potential and force fitted from the same conservative basis.**
6. **Near-nuclear singular regions excluded and handled by a separate hard core.**
7. **Equivalent substrate atoms share coefficients.**
8. **All periodic atomic images within support are included explicitly.**
9. **The smooth/local split is defined structurally, with optional exact orthogonalization rather than empirical penalties.**
10. **The final model is judged by force accuracy and physical scans, not only global potential RMSE.**
11. **Boltzmann weighting emphasizes the physically relevant attractive region over the repulsive wall.**

---

## 14. Validation results (NaCl 1x1 L3, σ=2.0 Å, T_eff=300 K)

| Metric | Value |
|--------|-------|
| n_train | 732 |
| n_valid | 680 |
| V_RMSE (unweighted) | 5.6e-04 eV |
| V_max | 5.4e-03 eV |
| F_RMSE | 2.2e-03 eV/Å |
| V_RMSE (Boltzmann-weighted) | 3.0e-04 eV |
| V_RMSE (attractive-only) | 2.3e-04 eV |

Plots generated:
- `linear_hybrid_xz_components.png`: 2×4 panel — reference, fit, error, B-spline, atom, Fz error, Boltzmann weight map (diagonal y=x cut)
- `linear_hybrid_zscan.png`: 1D vertical scan along diagonal (over Na at x=0, Cl at x=2)
