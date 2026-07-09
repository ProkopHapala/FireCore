# USER

I'm fitting parameters of classical potential to reproduce vibration spectrum. Perhaps it is better to fit directly from Hessian matrix elements.

In simple terms let's consider just bonding interactions (non-bonding are not present), so we have bonds, angles, maybe digedrals (torsions, improper).

Can you write explicit formuls for variational derivatives of Hessian matrix elements with respect to bond, angle, and dihedral stiffness constants. Assume typical formulation like in UFF. 

---

# ChatGPT 5.5

Yes. For a force field whose bonded terms are **linear in their stiffness parameters**, fitting Cartesian Hessian elements is especially attractive because the model Hessian is linear in those parameters at a fixed geometry. This is the principle behind several Hessian-fitting force-field procedures. ([American Chemical Society Publications][1])

UFF uses a harmonic bond term and cosine-Fourier forms for angle bending, proper torsions, and inversion terms. ([Materials Science Community Discourse][2])

# 1. General formula

Let

[
\mathbf x=(\mathbf r_1,\ldots,\mathbf r_N)\in\mathbb R^{3N}
]

and consider one bonded interaction (p),

[
E_p(\mathbf x)=k_p,F_p(q_p(\mathbf x)),
]

where (q_p) is a bond length, angle, dihedral, or improper coordinate.

Define its Wilson vector and coordinate Hessian

[
\mathbf b_p=\frac{\partial q_p}{\partial\mathbf x},
\qquad
\mathbf C_p=\frac{\partial^2q_p}{\partial\mathbf x,\partial\mathbf x}.
]

Then

[
\nabla_{\mathbf x}E_p
=====================

k_p F_p'(q_p)\mathbf b_p,
]

and the Cartesian Hessian contribution is

[
\boxed{
\mathbf H_p
===========

k_p
\left[
F_p''(q_p)\mathbf b_p\mathbf b_p^\mathrm T
+
F_p'(q_p)\mathbf C_p
\right].
}
]

Therefore, at fixed geometry,

[
\boxed{
\frac{\partial\mathbf H}{\partial k_p}
======================================

F_p''(q_p)\mathbf b_p\mathbf b_p^\mathrm T
+
F_p'(q_p)\mathbf C_p.
}
]

In element notation,

[
\boxed{
\frac{\partial H_{A\alpha,B\beta}}{\partial k_p}
================================================

F_p''(q_p)
\frac{\partial q_p}{\partial r_{A\alpha}}
\frac{\partial q_p}{\partial r_{B\beta}}
+
F_p'(q_p)
\frac{\partial^2q_p}
{\partial r_{A\alpha}\partial r_{B\beta}}.
}
]

The second term is sometimes called the **geometric-stiffness** or **prestress** contribution.

## At the internal-coordinate minimum

When

[
F_p'(q_p)=0,
]

the result simplifies to

[
\boxed{
\frac{\partial\mathbf H}{\partial k_p}
======================================

F_p''(q_p)\mathbf b_p\mathbf b_p^\mathrm T.
}
]

If (k_p) is defined to be the actual harmonic curvature,

[
\left.\frac{\partial^2 E_p}{\partial q_p^2}\right|_{q_0}=k_p,
]

then (F_p''(q_0)=1), and

[
\boxed{
\frac{\partial\mathbf H}{\partial k_p}
======================================

\mathbf b_p\mathbf b_p^\mathrm T.
}
]

Thus each scalar internal coordinate contributes a rank-one matrix to the Cartesian Hessian.

---

# 2. Bond stretching

For a bond (i-j), define

[
\mathbf d=\mathbf r_j-\mathbf r_i,
\qquad
r=|\mathbf d|,
\qquad
\mathbf e=\frac{\mathbf d}{r}.
]

The harmonic bond energy is

[
E_r=\frac12 k_r(r-r_0)^2.
]

The bond-coordinate derivatives are

[
\mathbf b_i=\frac{\partial r}{\partial\mathbf r_i}=-\mathbf e,
\qquad
\mathbf b_j=\frac{\partial r}{\partial\mathbf r_j}=+\mathbf e.
]

Define the transverse projector

[
\mathbf P=\mathbf I-\mathbf e\mathbf e^\mathrm T.
]

The second derivatives of the bond length are

[
\mathbf C_{ii}=\mathbf C_{jj}=\frac{\mathbf P}{r},
\qquad
\mathbf C_{ij}=\mathbf C_{ji}=-\frac{\mathbf P}{r}.
]

Because

[
F'(r)=r-r_0,
\qquad
F''(r)=1,
]

the exact derivative of the Hessian with respect to (k_r) is

[
\frac{\partial\mathbf H^{(r)}}{\partial k_r}
============================================

\begin{pmatrix}
\mathbf Q&-\mathbf Q\
-\mathbf Q&\mathbf Q
\end{pmatrix},
]

where

[
\boxed{
\mathbf Q
=========

\mathbf e\mathbf e^\mathrm T
+
\frac{r-r_0}{r}
\left(\mathbf I-\mathbf e\mathbf e^\mathrm T\right).
}
]

At (r=r_0),

[
\boxed{
\mathbf Q=\mathbf e\mathbf e^\mathrm T.
}
]

Therefore

[
\boxed{
\frac{\partial H_{i\alpha,i\beta}}{\partial k_r}
================================================

e_\alpha e_\beta,
\quad
\frac{\partial H_{i\alpha,j\beta}}{\partial k_r}
================================================

-e_\alpha e_\beta,
}
]

with analogous (jj) and (ji) blocks.

The transverse prestress term

[
\frac{r-r_0}{r}\mathbf P
]

is important when the bond is stretched or compressed at the reference geometry.

---

# 3. Angle bending

Consider the angle (i-j-k), with (j) central. Define

[
\mathbf a=\mathbf r_i-\mathbf r_j,
\qquad
\mathbf c=\mathbf r_k-\mathbf r_j,
]

[
a=|\mathbf a|,
\qquad
c_r=|\mathbf c|,
]

[
\mathbf u=\frac{\mathbf a}{a},
\qquad
\mathbf v=\frac{\mathbf c}{c_r},
]

[
\gamma=\mathbf u\cdot\mathbf v=\cos\theta,
\qquad
s=\sin\theta=\sqrt{1-\gamma^2}.
]

The angle-coordinate gradients are

[
\boxed{
\mathbf b_i
===========

\frac{\gamma\mathbf u-\mathbf v}{a,s},
}
]

[
\boxed{
\mathbf b_k
===========

\frac{\gamma\mathbf v-\mathbf u}{c_r,s},
}
]

and

[
\boxed{
\mathbf b_j=-\mathbf b_i-\mathbf b_k.
}
]

These three vectors automatically satisfy translational invariance.

## Harmonic angle

For

[
E_\theta=\frac12 k_\theta(\theta-\theta_0)^2,
]

the exact sensitivity is

[
\boxed{
\frac{\partial H^{(\theta)}*{A\alpha,B\beta}}
{\partial k*\theta}
===================

b_{A\alpha}b_{B\beta}
+
(\theta-\theta_0)
\frac{\partial^2\theta}
{\partial r_{A\alpha}\partial r_{B\beta}}.
}
]

At (\theta=\theta_0),

[
\boxed{
\frac{\partial\mathbf H^{(\theta)}}{\partial k_\theta}
======================================================

\begin{pmatrix}
\mathbf b_i\
\mathbf b_j\
\mathbf b_k
\end{pmatrix}
\begin{pmatrix}
\mathbf b_i^\mathrm T&
\mathbf b_j^\mathrm T&
\mathbf b_k^\mathrm T
\end{pmatrix}.
}
]

In (3\times3) block form,

[
\boxed{
\frac{\partial\mathbf H_{AB}}{\partial k_\theta}
================================================

\mathbf b_A\mathbf b_B^\mathrm T,
\qquad A,B\in{i,j,k}.
}
]

## UFF angle form

UFF uses

[
E_\theta
========

K_{ijk}
\left[
C_0+C_1\cos\theta+C_2\cos 2\theta
\right],
]

with

[
C_2=\frac{1}{4\sin^2\theta_0},
\qquad
C_1=-4C_2\cos\theta_0,
]

[
C_0=C_2(2\cos^2\theta_0+1).
]

Define

[
F(\theta)=C_0+C_1\cos\theta+C_2\cos2\theta.
]

Then

[
F'(\theta)
==========

-C_1\sin\theta-2C_2\sin2\theta,
]

[
F''(\theta)
===========

-C_1\cos\theta-4C_2\cos2\theta.
]

Therefore

[
\boxed{
\frac{\partial\mathbf H^{(\theta)}}{\partial K_{ijk}}
=====================================================

\left[-C_1\cos\theta-4C_2\cos2\theta\right]
\mathbf b_\theta\mathbf b_\theta^\mathrm T
+
\left[-C_1\sin\theta-2C_2\sin2\theta\right]
\mathbf C_\theta.
}
]

The UFF coefficients are normalized such that

[
F'(\theta_0)=0,
\qquad
F''(\theta_0)=1.
]

Consequently,

[
\boxed{
\left.
\frac{\partial\mathbf H^{(\theta)}}{\partial K_{ijk}}
\right|_{\theta_0}
==================

\mathbf b_\theta\mathbf b_\theta^\mathrm T.
}
]

## Explicit coordinate Hessian of an angle

For nonzero angular prestress, it is useful to calculate (\mathbf C_\theta) through (\gamma=\cos\theta):

[
\boxed{
\mathbf C_\theta
================

-\frac{1}{s}\mathbf C_\gamma
-\frac{\gamma}{s^3}
\mathbf g_\gamma\mathbf g_\gamma^\mathrm T,
}
]

where

[
\mathbf g_\gamma=\nabla\gamma.
]

In the six-dimensional coordinates ((\mathbf a,\mathbf c)), define

[
\mathbf w=\mathbf v-\gamma\mathbf u,
\qquad
\mathbf z=\mathbf u-\gamma\mathbf v,
]

[
\mathbf g_{\gamma,a}=\frac{\mathbf w}{a},
\qquad
\mathbf g_{\gamma,c}=\frac{\mathbf z}{c_r}.
]

The (3\times3) blocks of (\mathbf C_\gamma) are

[
\boxed{
\mathbf C_{\gamma,aa}
=====================

-\frac{
\mathbf w\mathbf u^\mathrm T+
\mathbf u\mathbf w^\mathrm T+
\gamma(\mathbf I-\mathbf u\mathbf u^\mathrm T)
}{a^2},
}
]

[
\boxed{
\mathbf C_{\gamma,cc}
=====================

-\frac{
\mathbf z\mathbf v^\mathrm T+
\mathbf v\mathbf z^\mathrm T+
\gamma(\mathbf I-\mathbf v\mathbf v^\mathrm T)
}{c_r^2},
}
]

[
\boxed{
\mathbf C_{\gamma,ac}
=====================

\frac{
\mathbf I-\mathbf v\mathbf v^\mathrm T
-\mathbf u\mathbf z^\mathrm T
}{a c_r},
\qquad
\mathbf C_{\gamma,ca}=\mathbf C_{\gamma,ac}^\mathrm T.
}
]

To transform these to atomic coordinates, use

[
\begin{pmatrix}
\mathbf a\
\mathbf c
\end{pmatrix}
=============

\underbrace{
\begin{pmatrix}
\mathbf I&-\mathbf I&0\
0&-\mathbf I&\mathbf I
\end{pmatrix}
}_{\mathbf T}
\begin{pmatrix}
\mathbf r_i\
\mathbf r_j\
\mathbf r_k
\end{pmatrix}.
]

Then

[
\mathbf C_\gamma^{\rm atoms}
============================

\mathbf T^\mathrm T
\mathbf C_\gamma^{(a,c)}
\mathbf T,
]

and similarly for (\mathbf C_\theta).

---

# 4. Proper dihedral or torsion

For atoms (i-j-k-l), define

[
\mathbf d_1=\mathbf r_j-\mathbf r_i,
\qquad
\mathbf d_2=\mathbf r_k-\mathbf r_j,
\qquad
\mathbf d_3=\mathbf r_l-\mathbf r_k,
]

[
\mathbf n_1=\mathbf d_1\times\mathbf d_2,
\qquad
\mathbf n_2=\mathbf d_2\times\mathbf d_3,
\qquad
L=|\mathbf d_2|.
]

One signed-dihedral convention is

[
\phi=
\operatorname{atan2}
\left(
L,\mathbf d_1\cdot\mathbf n_2,
\mathbf n_1\cdot\mathbf n_2
\right).
]

Define

[
a_1=\frac{\mathbf d_1\cdot\mathbf d_2}{L^2},
\qquad
a_3=\frac{\mathbf d_3\cdot\mathbf d_2}{L^2}.
]

The dihedral-coordinate gradients are

[
\boxed{
\mathbf b_i
===========

-\frac{L}{|\mathbf n_1|^2}\mathbf n_1,
}
]

[
\boxed{
\mathbf b_l
===========

+\frac{L}{|\mathbf n_2|^2}\mathbf n_2,
}
]

[
\boxed{
\mathbf b_j
===========

-(1+a_1)\mathbf b_i+a_3\mathbf b_l,
}
]

[
\boxed{
\mathbf b_k
===========

a_1\mathbf b_i-(1+a_3)\mathbf b_l.
}
]

A different sign convention for (\phi) changes all four (\mathbf b_A) signs simultaneously. The equilibrium Hessian contribution (\mathbf b\mathbf b^\mathrm T) is therefore unchanged.

## Harmonic torsion

For a locally harmonic torsion,

[
E_\phi=\frac12 k_\phi(\phi-\phi_0)^2,
]

[
\boxed{
\frac{\partial\mathbf H^{(\phi)}}{\partial k_\phi}
==================================================

\mathbf b_\phi\mathbf b_\phi^\mathrm T
+
(\phi-\phi_0)\mathbf C_\phi.
}
]

At the minimum,

[
\boxed{
\frac{\partial\mathbf H^{(\phi)}}{\partial k_\phi}
==================================================

\mathbf b_\phi\mathbf b_\phi^\mathrm T.
}
]

Thus every (3\times3) block is

[
\boxed{
\frac{\partial\mathbf H_{AB}}{\partial k_\phi}
==============================================

\mathbf b_A\mathbf b_B^\mathrm T,
\qquad
A,B\in{i,j,k,l}.
}
]

## UFF torsion form

UFF writes the torsional energy as

[
E_\phi
======

\frac{V_\phi}{2}
\left[
1-\cos(n\phi_0)\cos(n\phi)
\right].
]

Let

[
F(\phi)=\frac12
\left[
1-\cos(n\phi_0)\cos(n\phi)
\right].
]

Then

[
F'(\phi)
========

\frac n2\cos(n\phi_0)\sin(n\phi),
]

[
F''(\phi)
=========

\frac{n^2}{2}\cos(n\phi_0)\cos(n\phi).
]

Hence

[
\boxed{
\frac{\partial\mathbf H^{(\phi)}}{\partial V_\phi}
==================================================

\frac{n^2}{2}
\cos(n\phi_0)\cos(n\phi)
\mathbf b_\phi\mathbf b_\phi^\mathrm T
+
\frac n2
\cos(n\phi_0)\sin(n\phi)
\mathbf C_\phi.
}
]

At a UFF torsional minimum,

[
\cos(n\phi)=\cos(n\phi_0),
\qquad
\sin(n\phi)=0,
]

so

[
\boxed{
\left.
\frac{\partial\mathbf H^{(\phi)}}{\partial V_\phi}
\right|_{\phi_0}
================

\frac{n^2}{2}
\mathbf b_\phi\mathbf b_\phi^\mathrm T.
}
]

This distinction is important:

[
\boxed{
k_\phi^{\rm local}
==================

\left.
\frac{\partial^2E_\phi}{\partial\phi^2}
\right|_{\phi_0}
================

\frac{n^2V_\phi}{2}.
}
]

Thus (V_\phi) is a barrier parameter, not directly the local angular stiffness. If you fit the local stiffness (k_\phi^{\rm local}), then

[
\frac{\partial\mathbf H}{\partial k_\phi^{\rm local}}
=====================================================

\mathbf b_\phi\mathbf b_\phi^\mathrm T.
]

---

# 5. Improper torsion and UFF inversion

## Ordinary improper dihedral

If the improper is implemented simply as a signed dihedral,

[
E_\xi=\frac12k_\xi(\xi-\xi_0)^2,
]

use exactly the proper-dihedral formulas above, with the atom ordering appropriate to your improper convention:

[
\boxed{
\left.
\frac{\partial\mathbf H^{(\xi)}}{\partial k_\xi}
\right|_{\xi_0}
===============

\mathbf b_\xi\mathbf b_\xi^\mathrm T.
}
]

## Out-of-plane angle

For a central atom (i), plane atoms (j,k), and out-of-plane atom (l), define

[
\mathbf a=\mathbf r_j-\mathbf r_i,
\qquad
\mathbf c=\mathbf r_k-\mathbf r_i,
\qquad
\mathbf t=\mathbf r_l-\mathbf r_i,
]

[
\mathbf n=\mathbf a\times\mathbf c,
\qquad
\hat{\mathbf n}=\frac{\mathbf n}{|\mathbf n|},
\qquad
\hat{\mathbf t}=\frac{\mathbf t}{|\mathbf t|}.
]

Define the signed out-of-plane angle by

[
\sin\omega=\hat{\mathbf t}\cdot\hat{\mathbf n}.
]

Let

[
h=\hat{\mathbf t}\cdot\hat{\mathbf n},
\qquad
c_\omega=\sqrt{1-h^2},
]

[
\mathbf p_t=
\frac{\hat{\mathbf n}-h\hat{\mathbf t}}{|\mathbf t|},
\qquad
\mathbf p_n=
\frac{\hat{\mathbf t}-h\hat{\mathbf n}}{|\mathbf n|}.
]

Then

[
\boxed{
\mathbf b_l=\frac{\mathbf p_t}{c_\omega},
}
]

[
\boxed{
\mathbf b_j=\frac{\mathbf c\times\mathbf p_n}{c_\omega},
}
]

[
\boxed{
\mathbf b_k=-\frac{\mathbf a\times\mathbf p_n}{c_\omega},
}
]

[
\boxed{
\mathbf b_i=-\mathbf b_j-\mathbf b_k-\mathbf b_l.
}
]

For a harmonic out-of-plane term,

[
E_\omega=\frac12k_\omega(\omega-\omega_0)^2,
]

the equilibrium sensitivity is again

[
\boxed{
\frac{\partial\mathbf H^{(\omega)}}{\partial k_\omega}
======================================================

\mathbf b_\omega\mathbf b_\omega^\mathrm T.
}
]

UFF instead uses a one- or two-term cosine expansion,

[
E_\omega
========

K_\omega
\left(
C_0+C_1\cos\omega+C_2\cos2\omega
\right).
]

Thus

[
\boxed{
\frac{\partial\mathbf H^{(\omega)}}{\partial K_\omega}
======================================================

\left[-C_1\cos\omega-4C_2\cos2\omega\right]
\mathbf b_\omega\mathbf b_\omega^\mathrm T
+
\left[-C_1\sin\omega-2C_2\sin2\omega\right]
\mathbf C_\omega.
}
]

At a stationary inversion geometry, only the first term remains.

---

# 6. Shared atom-type parameters

Usually one parameter (k_t) is shared by many bonds or angles of type (t). Then

[
E_t=k_t\sum_{p\in t}F_p(q_p),
]

and

[
\boxed{
\frac{\partial\mathbf H}{\partial k_t}
======================================

\sum_{p\in t}
\frac{\partial\mathbf H_p}{\partial k_t}.
}
]

At termwise equilibrium,

[
\boxed{
\frac{\partial\mathbf H}{\partial k_t}
======================================

\sum_{p\in t}
\gamma_p\mathbf b_p\mathbf b_p^\mathrm T,
}
]

where

[
\gamma_p=
\begin{cases}
1, & \text{harmonic bond/angle/improper},[2mm]
n_p^2/2, & \text{UFF torsional barrier }V_p.
\end{cases}
]

---

# 7. Direct linear Hessian fitting

At a fixed geometry and with fixed equilibrium lengths, angles, periodicities, and phases,

[
\boxed{
\mathbf H_{\rm FF}
==================

\mathbf H_0+\sum_t k_t\mathbf A_t,
}
]

where

[
\mathbf A_t=\frac{\partial\mathbf H}{\partial k_t}.
]

Consequently, the stiffness fit is a linear least-squares problem:

[
\min_{\mathbf k}
\left|
\mathbf W^{1/2}
\left(
\mathbf H_0+\sum_t k_t\mathbf A_t-\mathbf H_{\rm ref}
\right)
\right|_F^2.
]

The normal equations are

[
\sum_q G_{pq}k_q=y_p,
]

with

[
G_{pq}
======

\left\langle
\mathbf A_p,\mathbf A_q
\right\rangle_W,
]

[
y_p
===

\left\langle
\mathbf A_p,
\mathbf H_{\rm ref}-\mathbf H_0
\right\rangle_W.
]

I would solve this by SVD or QR rather than forming the normal equations explicitly, because angle and torsion columns can be strongly correlated.

For nonnegative force constants, use nonnegative least squares or a bounded linear solver.

---

# 8. Mass-weighted Hessian

The vibrational dynamical matrix is

[
\mathbf D=\mathbf M^{-1/2}\mathbf H\mathbf M^{-1/2}.
]

Therefore,

[
\boxed{
\frac{\partial\mathbf D}{\partial k_t}
======================================

\mathbf M^{-1/2}
\frac{\partial\mathbf H}{\partial k_t}
\mathbf M^{-1/2}.
}
]

Elementwise,

[
\boxed{
\frac{\partial D_{A\alpha,B\beta}}{\partial k_t}
================================================

\frac{1}{\sqrt{m_A m_B}}
\frac{\partial H_{A\alpha,B\beta}}{\partial k_t}.
}
]

Fitting the mass-weighted Hessian may be more directly connected to frequencies, while fitting the Cartesian Hessian gives a more isotope-independent force field.

For a nondegenerate normal mode (s),

[
\mathbf D\mathbf v_s=\omega_s^2\mathbf v_s,
]

[
\boxed{
\frac{\partial\omega_s^2}{\partial k_t}
=======================================

\mathbf v_s^\mathrm T
\frac{\partial\mathbf D}{\partial k_t}
\mathbf v_s,
}
]

and

[
\boxed{
\frac{\partial\omega_s}{\partial k_t}
=====================================

\frac{1}{2\omega_s}
\mathbf v_s^\mathrm T
\frac{\partial\mathbf D}{\partial k_t}
\mathbf v_s.
}
]

Direct Hessian fitting avoids the mode-ordering and mode-crossing problems that appear when fitting frequencies individually.

---

# 9. Important prestress caveat

The familiar relation

[
\mathbf H=\mathbf B^\mathrm T\mathbf K\mathbf B
]

is exact only when the internal-coordinate forces vanish or when the geometric-stiffness terms are deliberately omitted.

The full relation is

[
\boxed{
\mathbf H
=========

\mathbf B^\mathrm T\mathbf K\mathbf B
+
\sum_p
\frac{\partial E}{\partial q_p}
\frac{\partial^2q_p}{\partial\mathbf x^2}.
}
]

Even if the **total Cartesian force is zero**, individual bonds and angles may be prestressed:

[
\sum_p
\frac{\partial E}{\partial q_p}\mathbf b_p=0
]

does not imply

[
\frac{\partial E}{\partial q_p}=0
]

for every interaction. This matters particularly in rings, strained molecules, and incompatible angle networks.

For practical fitting, I would use one of two consistent strategies:

1. Fit equilibrium bond lengths and angles so the reference geometry is nearly termwise force-free, then use the simple rank-one matrices.
2. Evaluate each unit-stiffness interaction Hessian at the actual reference geometry, including its (F'\mathbf C) term.

The second option remains completely linear in the stiffness constants and is the safer general formulation.

[1]: https://pubs.acs.org/doi/10.1021/acs.jcim.4c00540?utm_source=chatgpt.com "HessFit: A Toolkit to Derive Automated Force Fields from ..."
[2]: https://matsci.org/uploads/short-url/22ZVVhGYHxONHec5tvn4dv9mjfN.pdf "UFF, a full periodic table force field for molecular mechanics and molecular dynamics simulations"

---

# USER

OK we implemented it here, but it does not seem to work so well can you review it and eventually comment on it? Find bugs and tell if it is most effcijent way how to do it?

In genral I would prefer not to use sin and cosine function, I prefer to use dot product of normalized vectors ui, uj;   cij=dot(ui,uj); sij = sqrt(1-ci*cj)

maybe you can spot other things which are overcomplicated, inefficient or just plain wrong. I think we have problems 


```
Both tests pass (exit code 0). Results:

**Single-system (SiH4)**:
- LSQ vs GD max diff: 3.3e-04 (good agreement)
- At equilibrium, prestress = 0, only sin²θ scaling applies
- rel RMSD: 6.42%

**Multi-system (6 Si nanocrystals)**:
- LSQ vs GD max diff: 3.7e-04 (good agreement)
- Averaged equilibrium params computed with sensible values and small std deviations
- rel RMSD: 0.44–1.28% for larger NCs, 6.63% for SiH4 (expected — SiH4 is the outlier since its bonds/angles differ most from the multi-system average)

**Key observations**:
- `k_angle[Si-Si-H]` is negative (-0.125 eV/rad²), which is unphysical — this suggests the harmonic model can't fully capture the mixed angle contributions, or that more parameter types / cross-terms are needed
- The model Hessian norm is much smaller than the reference (e.g. ||H_model||=2.82 vs ||H_ref||=79.2 for SiH4), indicating the 2-body harmonic FF captures only a fraction of the full QM Hessian
- Frequency mismatch is large for SiH4 (it has high-frequency Si-H stretches ~2000-2500 cm⁻¹ that the fitted model underestimates), but the low-frequency modes of larger NCs are closerI think we have something wrong there, but write righ physical motivation/justifica/derivation including path formulas direcly in the code, atart in the header with general theory, and then in each function or key blokcs/lines put most relevant ideas and equation

Both tests pass (exit code 0). Results:

**Single-system (SiH4)**:
- LSQ vs GD max diff: 3.3e-04 (good agreement)
- At equilibrium, prestress = 0, only sin²θ scaling applies
- rel RMSD: 6.42%

**Multi-system (6 Si nanocrystals)**:
- LSQ vs GD max diff: 3.7e-04 (good agreement)
- Averaged equilibrium params computed with sensible values and small std deviations
- rel RMSD: 0.44–1.28% for larger NCs, 6.63% for SiH4 (expected — SiH4 is the outlier since its bonds/angles differ most from the multi-system average)

**Key observations**:
- `k_angle[Si-Si-H]` is negative (-0.125 eV/rad²), which is unphysical — this suggests the harmonic model can't fully capture the mixed angle contributions, or that more parameter types / cross-terms are needed
- The model Hessian norm is much smaller than the reference (e.g. ||H_model||=2.82 vs ||H_ref||=79.2 for SiH4), indicating the 2-body harmonic FF captures only a fraction of the full QM Hessian
- Frequency mismatch is large for SiH4 (it has high-frequency Si-H stretches ~2000-2500 cm⁻¹ that the fitted model underestimates), but the low-frequency modes of larger NCs are closer
```

---

# Chat GPT 5.5

## Verdict

Your suspicion is justified. The implementation contains one **definite mathematical bug**, several conceptual inconsistencies, and two likely sources of serious multi-system fitting errors. The agreement between LSQ and gradient descent does **not** validate the derivatives because both paths use essentially the same sensitivity formulas. 

I made a corrected version:

[Download the reviewed header](sandbox:/mnt/data/FFfit_reviewed.h)

[Download the unified diff](sandbox:/mnt/data/FFfit_review.patch)

I syntax-compiled it against small API-compatible `Vec3d`/`Mat3d` stubs and finite-difference tested the corrected angle-coordinate Hessian; the maximum error was (7.4\times10^{-7}). I could not compile it against your actual project because the real `Vec3.h`, `Mat3.h`, and tests were not included.

---

# 1. Definite bug in the angle coordinate Hessian

You currently construct

```cpp
// C_jk = -C_ki - C_kk
C[1][2].set(C[2][0]);
C[1][2].add(C[2][2]);
C[1][2].mul(-1.0);
```

But the correct relation is

[
\boxed{
C_{jk}=-(C_{ik}+C_{kk}),
}
]

not ( -(C_{ki}+C_{kk}) ).

The corrected code is:

```cpp
// C_jk = -(C_ik + C_kk)
C[1][2].set(C[0][2]);
C[1][2].add(C[2][2]);
C[1][2].mul(-1.0);

C[2][1].setT(C[1][2]);
```

Although (C_{ki}=C_{ik}^{T}), these matrices are generally not individually symmetric, so they cannot be exchanged.

This bug has an especially revealing consequence:

[
\sum_B C_{AB}\ne 0,
]

so the old prestress Hessian violates translational invariance and the acoustic sum rule. It can introduce spurious restoring forces for uniform translation.

It only contributes through

[
(c-c_0)C_c,
]

so it disappears when every angle is exactly at its individual equilibrium. That explains why the single-system equilibrium test passes while a multi-system fit using averaged equilibrium angles can behave badly.

---

# 2. Do not calculate (\partial\theta/\partial x) at all

Your current path is mathematically equivalent at ordinary angles, but unnecessarily complicated:

1. Compute

[
b_\theta=\frac{\partial\theta}{\partial x}
\propto \frac{1}{\sin\theta}.
]

2. Then multiply

[
b_\theta b_\theta^T
]

by (\sin^2\theta).

That introduces a large-small cancellation and a fake singularity near linear angles.

Since your energy is expressed through

[
c=\cos\theta=\mathbf u\cdot\mathbf v,
]

compute

[
g=\frac{\partial c}{\partial x}
]

directly.

For

[
\mathbf a=\mathbf r_i-\mathbf r_j,\qquad
\mathbf b=\mathbf r_k-\mathbf r_j,
]

[
\mathbf u=\frac{\mathbf a}{|\mathbf a|},\qquad
\mathbf v=\frac{\mathbf b}{|\mathbf b|},\qquad
c=\mathbf u\cdot\mathbf v,
]

the gradients are simply

[
\boxed{
\mathbf g_i=\frac{\mathbf v-c\mathbf u}{|\mathbf a|}
}
]

[
\boxed{
\mathbf g_k=\frac{\mathbf u-c\mathbf v}{|\mathbf b|}
}
]

[
\boxed{
\mathbf g_j=-\mathbf g_i-\mathbf g_k.
}
]

This uses only normalization and dot products. It needs neither `sin`, `cos`, `acos`, nor even

```cpp
sqrt(1.0 - c*c)
```

in the hot path.

---

# 3. Your present angle parameter is not an angular stiffness

The current energy is

[
E=K(c-c_0)^2.
]

At equilibrium,

[
\left.
\frac{\partial^2E}{\partial\theta^2}
\right|_{\theta_0}
==================

2K\sin^2\theta_0.
]

Therefore:

[
\boxed{
k_\theta^{\rm local}=2K(1-c_0^2).
}
]

So the fitted (K) has units of energy, not directly “eV/rad²”. Reporting the present `k_angle` as an angular stiffness is incorrect.

This also makes parameters for different equilibrium angles difficult to compare.

## Better UFF-normalized form

The standard non-linear UFF Fourier expression is exactly equivalent to

[
\boxed{
E_\theta=
\frac12 k_\theta
\frac{(c-c_0)^2}{1-c_0^2}.
}
]

This follows from

[
C_0+C_1\cos\theta+C_2\cos2\theta
================================

\frac{(c-c_0)^2}{2(1-c_0^2)}
]

with the usual UFF coefficients.

Now (k_\theta) is the actual local angular curvature:

[
\left.
\frac{\partial^2E}{\partial\theta^2}
\right|*{\theta_0}
=k*\theta.
]

The exact Hessian sensitivity is

[
\boxed{
\frac{\partial H}{\partial k_\theta}
====================================

\frac{
gg^T+(c-c_0)C_c
}{
1-c_0^2
}.
}
]

At equilibrium:

[
\boxed{
\frac{\partial H}{\partial k_\theta}
====================================

# \frac{gg^T}{1-c_0^2}

b_\theta b_\theta^T.
}
]

This is the form used in the reviewed header.

For a tetrahedral angle, (c_0=-1/3) and (1-c_0^2=8/9), so the numerical scaling difference is not enormous, but the physical meaning and conditioning are much cleaner.

---

# 4. The reported norms and RMSD are mutually inconsistent

You report

[
|H_{\rm model}|=2.82,
\qquad
|H_{\rm ref}|=79.2,
]

while claiming a relative RMSD of (6.42%).

For the same matrices, the triangle inequality alone gives

[
|H_{\rm model}-H_{\rm ref}|
\ge
79.2-2.82=76.38.
]

Therefore,

[
\frac{|H_{\rm model}-H_{\rm ref}|}{|H_{\rm ref}|}
\ge 96.4%.
]

For an unweighted LSQ projection with (H_0=0), the stronger orthogonal-projection identity holds:

[
|H_{\rm ref}|^2
===============

|H_{\rm model}|^2
+
|H_{\rm ref}-H_{\rm model}|^2,
]

which gives a relative residual of about

[
99.94%.
]

Thus at least one of the following is happening:

* the “relative RMSD” uses an unrelated normalization, such as matrix range or maximum element;
* the norms and RMSD use differently preprocessed matrices;
* one is Cartesian and the other mass-weighted;
* one includes (H_0) and the reported model norm does not;
* the matrix layouts or units differ;
* the reported model Hessian is not reconstructed with the fitted parameters.

This should be investigated before interpreting any fitted constants.

A useful invariant to print immediately after LSQ is

```cpp
dot(H_model, H_ref - H0 - H_model)
```

which should be approximately zero for an unweighted unconstrained fit.

Also print

```cpp
rel_frobenius =
    norm(H0 + H_model - H_ref) / norm(H_ref);
```

rather than a range-normalized element RMS.

---

# 5. Multi-system parameter indices may be mixed up

The current `auto_assign_types()` creates a new local type map based on the order in which types first appear.

That does not guarantee consistent indices across systems.

For example:

* SiH(_4): parameter 0 might be `Si-H`;
* a larger nanocrystal: parameter 0 might be `Si-Si`.

If their normal equations are then simply added by array index, physically unrelated parameters are mixed.

The class comments discuss global multi-system accumulation, but the implementation provides only local numbering.

The reviewed header introduces a shared `ParamRegistry`. For a multi-system fit, the registry should preferably be constructed in a first pass and then used to assign every system.

---

# 6. Angle types are not symmetric under (i\leftrightarrow k)

The original key is

```cpp
(symbol[i], symbol[j], symbol[k])
```

so

```text
H-Si-Si
```

and

```text
Si-Si-H
```

become different parameters even though they describe the same angle with reversed outer atoms.

The key must preserve the central atom but sort the two outer atoms:

```cpp
(min(symbol[i], symbol[k]),
 symbol[j],
 max(symbol[i], symbol[k]))
```

The reviewed header does this.

---

# 7. Averaging equilibrium angles is not rigorous

The header currently states that averaging (r_0) and (\theta_0) first and then fitting stiffnesses is rigorous. It is not.

First, for an energy expressed through (c=\cos\theta), the relevant coordinate average is

[
c_0\approx\langle \mathbf u\cdot\mathbf v\rangle,
]

not

[
\cos\langle\theta\rangle.
]

These differ because cosine is nonlinear.

More importantly, a relaxed structure only satisfies

[
\sum_t f_t'(q_t)g_t=0.
]

It does not imply

[
f_t'(q_t)=0
]

for each term. Bonds and angles can be mutually prestressed while their Cartesian forces cancel.

There are three defensible strategies:

1. **Pure local Hessian model:** set each term’s (q_0) to its current value and share only stiffnesses. This gives zero prestress but is not one globally transferable energy surface.

2. **Transferable energy model:** use shared (r_0,c_0) and fit forces together with Hessians.

3. **Linear reparameterization:** fit (k) and (kr_0), or (k) and (kc_0), because the Hessian is linear in these combinations. Then recover (r_0) or (c_0), with physical bounds and force constraints.

Simply averaging coordinates is a reasonable initial guess, but not a derivation.

---

# 8. The negative angle constant

A negative fitted coefficient can have several causes here:

1. The incorrect prestress block.
2. Mixed parameter indices across systems.
3. Reversed angle types being treated separately.
4. Wrong angle-constant units/convention.
5. An ill-conditioned normal equation.
6. Missing force-field terms.

After fixing the first five, a negative coefficient can indeed indicate insufficient model flexibility.

For SiH(_4), one shared Si–H stretch constant and one H–Si–H bend constant are generally too restrictive to reproduce the complete quantum Hessian. Symmetric and asymmetric modes also reflect stretch–stretch and bend–bend couplings. Useful additions would be:

[
(r_{ij}-r_0)(r_{ik}-r_0)
]

for stretch–stretch coupling,

[
(r_{ij}-r_0)(c-c_0)
]

for stretch–bend coupling, and possibly Urey–Bradley (i-k) terms.

I would not add these until the present bugs and unit checks are resolved.

---

# 9. Solver problems

The current singular-pivot handling is unsafe:

```cpp
if (maxval < 1e-15) {
    printf("WARNING ...");
    continue;
}
```

The back-substitution later still executes

```cpp
x[i] /= G[i*n+i];
```

and can divide by zero.

Other issues are:

* normal equations square the condition number;
* columns are not scaled;
* no condition or rank diagnostic is reported;
* unconstrained LSQ naturally permits negative stiffnesses.

The reviewed version:

* scales every sensitivity column by (\sqrt{G_{pp}});
* adds a small relative ridge;
* fails loudly on a singular pivot;
* never continues into invalid back-substitution.

For production-quality fitting, pivoted QR or SVD is still preferable. For only a few parameter types, scaled normal equations are usually adequate.

After validation, nonnegative least squares is appropriate if all fitted coefficients must be physical.

---

# 10. Efficiency

The most important efficiency changes are:

### Remove the angle trigonometric pathway

Replace

```cpp
sqrt(1-c*c)
1/s
b_theta
s*s*b_theta*b_theta
cos(theta0)
acos(c)
```

with direct (g=\partial c/\partial x). This is both faster and more stable.

### Simplify bond sensitivity

The four bond blocks are all generated by one matrix

[
Q=ee^T+\frac{r-r_0}{r}(I-ee^T).
]

So only one `Mat3d` needs to be constructed:

```cpp
Q, -Q, -Q, Q
```

### Do not use gradient descent for the fit

The problem is linear. Rebuilding an (N_3\times N_3) model Hessian and recomputing every term for hundreds of iterations is much more expensive than one linear solve.

LSQ–GD agreement is useful only as a temporary debugging check.

### Precompute sensitivities

At fixed geometry, all (A_p) matrices are constant. If gradient optimization is retained, precompute

[
G_{pq}=\langle A_p,A_q\rangle,
\qquad
y_p=\langle A_p,H_{\rm target}\rangle,
]

and use

[
\nabla L=2(Gk-y).
]

Then each iteration costs (O(n_p^2)), not repeated Hessian assembly.

### Use sparse (3\times3) blocks for large systems

Each bond touches four blocks and each angle touches nine blocks. Storing a full dense matrix per parameter costs

[
O(n_pN^2),
]

although almost all model blocks are zero.

For larger nanocrystals, aggregate sensitivity matrices as sparse atom-pair blocks. Then:

* (y_p) reads only reference blocks touched by parameter (p);
* (G_{pq}) is computed only where the two sparse block sets overlap;
* no dense (A_p) matrices are required.

For modest molecules, the dense implementation remains simpler and is probably fast enough.

---

# 11. Smaller code issues

The reviewed file also addresses these:

* The bond-coordinate Hessian comments have the wrong sign. The implementation is correct:

  [
  C_{ii}=C_{jj}=+\frac{I-ee^T}{r},
  \qquad
  C_{ij}=C_{ji}=-\frac{I-ee^T}{r}.
  ]

* `bond_energy_grad()` actually adds forces, not gradients.

* `set_geometry()` leaks memory when called more than once.

* Default copying of `FFfit` causes double deletion of `apos`.

* `<tuple>` is used but not explicitly included.

* `weight` is squared internally; it is really a `sqrt_weight`, not a direct statistical weight.

* `theta = acos(...)` in `angle_dHdk_uff()` is unused.

* `UFFAngleParams::K` and `C3` are unused.

* The supposedly exact UFF path skips its prestress contribution.

* No dihedral or improper term is implemented yet.

The first tests I would rerun are: finite-difference each individual bond and angle sensitivity away from equilibrium, verify every (3\times3) block-row sums to zero, verify Hessian symmetry, and verify the LSQ projection identity before looking at frequencies.

---

# USER




I’ll target the likely physical issues from the review: angle stiffness convention, averaging `cosθ` instead of `θ`, and diagnostics that can reveal whether the fit is physically meaningful. First I’ll inspect the current averaging and print paths.
I’ll correct the actual physical convention now: use an angle parameter whose fitted value is the local angular curvature (`eV/rad²`), not the coefficient of `(cosθ-cosθ0)²`; average equilibrium angles via `<cosθ>`; and update the Python/C++ mirrors consistently.







I found the main physical error: angle coefficients were still the coefficient of `(cosθ-cosθ0)^2`, while the code printed them as angular stiffness. I changed C++ to the normalized form where fitted `k_angle` is really `∂²E/∂θ²`. Now I’ll bring the Python reference and test plumbing into the same convention.




The Python plumbing now uses the same physical convention and `c0=<cosθ>` averaging. I’ll search for stale old-formula code, then rebuild/run.



Build passed. I’ll run single and multi-system tests with full output.

Single-system shows the units fix worked: angle stiffness is now `1.066 eV/rad²` instead of the old unnormalized coefficient. The Frobenius residual reveals the model captures only a small part of the QM Hessian (`96%` residual), so the earlier low relative error was misleading. Now I’ll run multi-system.



## What I corrected

I fixed the main physical convention errors in `@/home/prokop/git/FireCore/cpp/common/molecular/FFfit.h:60-68` and mirrored them in `@/home/prokop/git/FireCore/tests/tSiNCs/test_FFfit.py:305-306`.

- **Angle stiffness units:** changed angle fit from coefficient of `(cosθ-cosθ0)^2` to true angular stiffness:
  `E = 0.5*kθ*(cosθ-c0)^2/(1-c0^2)`

- **Direct angle gradient:** C++ now uses `g = ∂cosθ/∂r` directly, avoiding artificial `1/sinθ` singularity in `∂θ/∂r` at `@/home/prokop/git/FireCore/cpp/common/molecular/FFfit.h:339-383`.

- **Equilibrium angle averaging:** Python now averages `c0=<cosθ>`, not `cos(<θ>)`, at `@/home/prokop/git/FireCore/tests/tSiNCs/test_FFfit.py:496-535`.

- **Symmetric angle keys:** outer atoms are sorted consistently via [angle_type_key()](cci:1://file:///home/prokop/git/FireCore/tests/tSiNCs/test_FFfit.py:24:0-26:69) at `@/home/prokop/git/FireCore/tests/tSiNCs/test_FFfit.py:25-27`.

- **Physical residual reporting:** test now reports Frobenius relative residual, not misleading RMSD/norm combination, at `@/home/prokop/git/FireCore/tests/tSiNCs/test_FFfit.py:724-737`.

I also saved a project memory with this corrected FFfit convention.

## Test results

### Single-system `SiH4`

- **Build:** passed.
- **LSQ/GD agreement:** `5.46e-4`, good.
- **Fitted parameters:**
  - `k_Si-H = 4.3708 eV/A²`
  - `k_H-Si-H = 1.0662 eV/rad²`

This is more physically meaningful than before because `k_angle` is now truly angular curvature.

But the fit quality is poor:

- `||H_ref|| = 79.18`
- `||H_model|| = 21.55`
- `||diff|| = 76.19`
- `relFrob = 96.23%`

So for `SiH4`, the simple bond+angle model captures only a small subspace of the QM Hessian.

### Multi-system `all_Si`

- **Build/test:** passed.
- **LSQ/GD agreement:** `2.52e-4`, good.
- **Fitted parameters:**
  - `k_Si-H = 0.079 eV/A²`
  - `k_Si-Si = 0.233 eV/A²`
  - `k_H-Si-H = 0.447 eV/rad²`
  - `k_Si-Si-Si = 0.150 eV/rad²`
  - `k_H-Si-Si = -0.221 eV/rad²`

The negative `H-Si-Si` angle stiffness remains unphysical.

## Physical interpretation

The derivative/unit bugs are now corrected, but the model still does **not** make physical sense as a full QM Hessian fit:

- **Negative angle stiffness remains:** means the model basis is insufficient or the unconstrained LSQ is compensating for missing terms.
- **Frobenius residual is ~99-100% for multi-system:** the model Hessian norm is tiny compared to `H_ref`.
- **Frequencies remain far too low:** high-frequency modes are not captured.

## Conclusion

The implementation is now mathematically more physically consistent, but the current force-field form is too poor: only independent harmonic bonds + independent harmonic angles cannot reproduce the QM Hessian.

Next physical fixes should be model-level, not derivative-level:

- **Add coupling terms:** stretch-stretch, stretch-bend, bend-bend.
- **Try non-negative least squares:** to prevent unphysical negative stiffness.
- **Fit forces + Hessians:** if shared `l0/c0` should be a transferable energy model.
- **Use local per-system equilibria first:** to isolate whether bad results come from shared prestress or missing force-field terms.

---

# USER

there is some critique, but I do not think it is too justified. We should do what we can we should find effective nearest neighbor description using reduced model to fit vibration frequencies of the full model. Maybe fitting of hessian is not the bst way? aybe we should vary parametrs to get proper eigenvalues and eigenvectros? But I was thinking fitting hessian is more direcp. How it is normally done. I guess this is common problem to fit QM vibration spectra using just bond and angle spring (nearest-neighbor, nex nearest neighbor)

I’ll target the likely physical issues from the review: angle stiffness convention, averaging `cosθ` instead of `θ`, and diagnostics that can reveal whether the fit is physically meaningful. First I’ll inspect the current averaging and print paths. I’ll correct the actual physical convention now: use an angle parameter whose fitted value is the local angular curvature (`eV/rad²`), not the coefficient of `(cosθ-cosθ0)²`; average equilibrium angles via `<cosθ>`; and update the Python/C++ mirrors consistently. I found the main physical error: angle coefficients were still the coefficient of `(cosθ-cosθ0)^2`, while the code printed them as angular stiffness. I changed C++ to the normalized form where fitted `k_angle` is really `∂²E/∂θ²`. Now I’ll bring the Python reference and test plumbing into the same convention. The Python plumbing now uses the same physical convention and `c0=<cosθ>` averaging. I’ll search for stale old-formula code, then rebuild/run. Build passed. I’ll run single and multi-system tests with full output. Single-system shows the units fix worked: angle stiffness is now `1.066 eV/rad²` instead of the old unnormalized coefficient. The Frobenius residual reveals the model captures only a small part of the QM Hessian (`96%` residual), so the earlier low relative error was misleading. Now I’ll run multi-system. ## What I corrected I fixed the main physical convention errors in `@/home/prokop/git/FireCore/cpp/common/molecular/FFfit.h:60-68` and mirrored them in `@/home/prokop/git/FireCore/tests/tSiNCs/test_FFfit.py:305-306`. - **Angle stiffness units:** changed angle fit from coefficient of `(cosθ-cosθ0)^2` to true angular stiffness: `E = 0.5*kθ*(cosθ-c0)^2/(1-c0^2)` - **Direct angle gradient:** C++ now uses `g = ∂cosθ/∂r` directly, avoiding artificial `1/sinθ` singularity in `∂θ/∂r` at `@/home/prokop/git/FireCore/cpp/common/molecular/FFfit.h:339-383`. - **Equilibrium angle averaging:** Python now averages `c0=<cosθ>`, not `cos(<θ>)`, at `@/home/prokop/git/FireCore/tests/tSiNCs/test_FFfit.py:496-535`. - **Symmetric angle keys:** outer atoms are sorted consistently via [angle_type_key()](cci:1://file:///home/prokop/git/FireCore/tests/tSiNCs/test_FFfit.py:24:0-26:69) at `@/home/prokop/git/FireCore/tests/tSiNCs/test_FFfit.py:25-27`. - **Physical residual reporting:** test now reports Frobenius relative residual, not misleading RMSD/norm combination, at `@/home/prokop/git/FireCore/tests/tSiNCs/test_FFfit.py:724-737`. I also saved a project memory with this corrected FFfit convention. ## Test results ### Single-system `SiH4` - **Build:** passed. - **LSQ/GD agreement:** `5.46e-4`, good. - **Fitted parameters:** - `k_Si-H = 4.3708 eV/A²` - `k_H-Si-H = 1.0662 eV/rad²` This is more physically meaningful than before because `k_angle` is now truly angular curvature. But the fit quality is poor: - `||H_ref|| = 79.18` - `||H_model|| = 21.55` - `||diff|| = 76.19` - `relFrob = 96.23%` So for `SiH4`, the simple bond+angle model captures only a small subspace of the QM Hessian. ### Multi-system `all_Si` - **Build/test:** passed. - **LSQ/GD agreement:** `2.52e-4`, good. - **Fitted parameters:** - `k_Si-H = 0.079 eV/A²` - `k_Si-Si = 0.233 eV/A²` - `k_H-Si-H = 0.447 eV/rad²` - `k_Si-Si-Si = 0.150 eV/rad²` - `k_H-Si-Si = -0.221 eV/rad²` The negative `H-Si-Si` angle stiffness remains unphysical. ## Physical interpretation The derivative/unit bugs are now corrected, but the model still does **not** make physical sense as a full QM Hessian fit: - **Negative angle stiffness remains:** means the model basis is insufficient or the unconstrained LSQ is compensating for missing terms. - **Frobenius residual is ~99-100% for multi-system:** the model Hessian norm is tiny compared to `H_ref`. - **Frequencies remain far too low:** high-frequency modes are not captured. ## Conclusion The implementation is now mathematically more physically consistent, but the current force-field form is too poor: only independent harmonic bonds + independent harmonic angles cannot reproduce the QM Hessian. Next physical fixes should be model-level, not derivative-level: - **Add coupling terms:** stretch-stretch, stretch-bend, bend-bend. - **Try non-negative least squares:** to prevent unphysical negative stiffness. - **Fit forces + Hessians:** if shared `l0/c0` should be a transferable energy model. - **Use local per-system equilibria first:** to isolate whether bad results come from shared prestress or missing force-field terms.

---

# Chat GPT 5.5

You are right: the earlier conclusion “the model Hessian is very different, therefore the reduced force field is useless” is too strong.

A nearest-neighbor bond–angle model is **not intended to reproduce every element of the full QM Hessian**. It is intended to reproduce a selected part of the vibrational dynamics—typically acoustic modes, elastic response, important optical branches, or molecular normal modes. Hessian-based extraction and fitting are standard approaches: Seminario-type methods extract bonded constants from Cartesian QM Hessians, while full, partial, and internal Hessian fitting optimize a reduced valence-force model against selected Hessian information. ([ResearchGate][1])

For silicon specifically, this is exactly the motivation of the Keating valence-force model and its extensions. The simplest Keating model contains bond stretching and bond bending; modified VFF models add stretch–stretch, stretch–bend, and bend–bend/coplanar terms to improve phonon dispersion. ([APS Link][2])

Your current implementation performs an ordinary elementwise Cartesian-Hessian least-squares fit,

[
\min_{\mathbf k}
\left|
\sum_p k_p A_p-H_{\rm QM}
\right|_F^2,
]

with (A_p=\partial H/\partial k_p). That is mathematically valid, linear, and efficient. The problem is mainly that this particular norm may not correspond to what you actually care about. 

# 1. What full Cartesian Hessian fitting actually optimizes

Let

[
D_{\rm QM}=M^{-1/2}H_{\rm QM}M^{-1/2}
]

be the mass-weighted reference dynamical matrix and

[
D(\mathbf k)=\sum_p k_p D_p
]

your reduced model.

Elementwise fitting finds the projection of (D_{\rm QM}) onto the small linear space spanned by (D_p):

[
D_{\rm fit}
===========

\operatorname*{argmin}*{D\in{\rm span}{D_p}}
|D-D*{\rm QM}|_F^2.
]

That is a sensible definition of the “nearest reduced Hessian,” but it weights every matrix element equally.

It does **not** specifically minimize:

* low-frequency errors;
* relative frequency errors;
* errors in selected modes;
* acoustic velocities;
* optical band positions;
* errors in collective versus localized modes.

Therefore, a large global Frobenius residual is unsurprising when the true QM Hessian contains interactions absent from the reduced model.

However, the very small reported model norm still tells us something precise:

[
|H_{\rm model}|\ll |H_{\rm QM}|
]

means that the reduced model does not reproduce most of the full matrix in that norm. This is not necessarily fatal, but then the evaluation must explicitly show that the **selected modes** are reproduced.

# 2. Should you fit eigenvalues and eigenvectors instead?

Directly fitting recomputed eigenvalues is possible, but I would not use it as the first fitting method.

A naïve spectral objective is

[
L_\omega(\mathbf k)
===================

\sum_s w_s
\left[
\omega_s(\mathbf k)-\omega_s^{\rm QM}
\right]^2.
]

Its disadvantages are:

* nonlinear optimization;
* repeated diagonalization;
* mode reordering and crossings;
* ambiguity inside degenerate subspaces;
* potentially many parameter sets with similar eigenvalues but different force responses.

Fitting eigenvalues alone is especially underconstrained: distinct Hessians can have the same spectrum.

But there is a much better formulation that fits the QM eigenvalues **and eigenvectors while remaining a linear least-squares problem**.

# 3. Recommended method: fit the action of the dynamical matrix on QM modes

Let the reference normal modes satisfy

[
D_{\rm QM}v_s=\lambda_s v_s,
\qquad
\lambda_s=\omega_s^2.
]

Instead of fitting every element of (D), minimize

[
\boxed{
L_{\rm mode}(\mathbf k)
=======================

\sum_{s\in\mathcal T}
w_s
\left|
D(\mathbf k)v_s-\lambda_s v_s
\right|^2.
}
]

Here (\mathcal T) is the set of target modes.

Because

[
D(\mathbf k)=\sum_p k_pD_p,
]

we have

[
D(\mathbf k)v_s
===============

\sum_p k_p(D_pv_s).
]

Define

[
a_{sp}=D_pv_s,
\qquad
b_s=\lambda_s v_s.
]

Then the objective is simply

[
\boxed{
\min_{\mathbf k}
\sum_s w_s
\left|
\sum_p k_p a_{sp}-b_s
\right|^2,
}
]

which is again an ordinary linear least-squares problem.

This is probably the best objective for your stated purpose.

It simultaneously asks:

1. Is (v_s^{\rm QM}) approximately an eigenvector of the force-field matrix?
2. Does it have approximately the correct eigenvalue (\lambda_s)?
3. Does applying (D_{\rm FF}) to this mode generate spurious components in other directions?

The last point is important. Merely fitting

[
v_s^T D_{\rm FF}v_s=\lambda_s
]

fits the Rayleigh quotient but does not guarantee that (v_s) is an eigenvector. The vector-residual formulation does.

# 4. Equivalent formulation in the QM normal-mode basis

Collect target eigenvectors into

[
V=[v_1,v_2,\ldots,v_m],
\qquad
\Lambda=\operatorname{diag}(\lambda_1,\ldots,\lambda_m).
]

Then fit

[
\boxed{
L_{\rm mode}
============

\left|
\left[D(\mathbf k)V-V\Lambda\right]W
\right|_F^2.
}
]

This contains two kinds of errors.

Inside the target subspace,

[
V^T D(\mathbf k)V-\Lambda
]

measures frequency errors and unwanted coupling between reference modes.

Outside it,

[
(I-VV^T)D(\mathbf k)V
]

measures leakage of the target modes into the rest of configuration space.

So the full residual

[
D(\mathbf k)V-V\Lambda
]

is stronger and cleaner than comparing frequencies alone.

# 5. Degenerate and nearly degenerate modes

For degenerate modes, individual eigenvectors are not unique. Any rotation inside the degenerate eigenspace is equally valid.

Therefore, do not fit individual vectors inside, for example, a triply degenerate (T_2) manifold. Fit the whole subspace.

Let (V_g) contain all vectors in a degenerate group. Then use

[
\boxed{
L_g=
\left|
D(\mathbf k)V_g-V_g\Lambda_g
\right|_F^2.
}
]

If the target eigenvalue is exactly degenerate,

[
\Lambda_g=\lambda_g I.
]

This objective is invariant under

[
V_g\rightarrow V_gR
]

for any orthogonal rotation (R) inside the group.

That is particularly relevant for SiH(_4), where several modes belong to multidimensional irreducible representations.

# 6. Weighting is crucial

Since

[
\lambda_s=\omega_s^2,
]

an absolute error in (\lambda) overweights high-frequency modes.

For small errors,

[
\frac{\delta\omega_s}{\omega_s}
\approx
\frac{\delta\lambda_s}{2\lambda_s}.
]

Therefore, if you want approximately equal **relative frequency accuracy**, normalize each mode residual by its eigenvalue:

[
\boxed{
L_{\rm rel}
===========

\sum_s
w_s
\frac{
|D(\mathbf k)v_s-\lambda_sv_s|^2
}{
\max(\lambda_s,\lambda_{\rm floor})^2
}.
}
]

The floor prevents huge weights for translations, rotations, and very soft modes.

A practical categorization for your systems would be:

[
L =
w_{\rm acoustic}L_{\rm acoustic}
+
w_{\rm SiSi}L_{\rm SiSi}
+
w_{\rm surface}L_{\rm surface}
+
w_{\rm SiH}L_{\rm SiH}.
]

For example, the Si–H stretching modes around the top of the spectrum should not dominate all lower-frequency Si lattice modes merely because their (\omega^2) is much larger.

You should remove:

* three translations for a periodic cell;
* three translations plus three rotations for a nonlinear isolated cluster;
* five rigid modes for a linear molecule.

# 7. Hybrid Hessian–mode fitting

I would not entirely discard Hessian fitting. A hybrid objective is better:

[
\boxed{
L(\mathbf k)
============

\alpha L_{\rm mode}
+
\beta L_{\rm localH}
+
\rho L_{\rm reg}.
}
]

Here

[
L_{\rm mode}
============

\sum_{s\in\mathcal T}
w_s
\left|
D(\mathbf k)v_s-\lambda_sv_s
\right|^2
]

targets the vibration spectrum.

A local-Hessian objective can use only atom pairs that the reduced model can physically represent:

[
L_{\rm localH}
==============

\sum_{|R_i-R_j|<R_{\rm fit}}
w_{ij}
\left|
H_{ij}^{\rm FF}-H_{ij}^{\rm QM}
\right|_F^2.
]

This is more sensible than penalizing every long-range QM Hessian block that a nearest-neighbor model is structurally incapable of generating.

Finally,

[
L_{\rm reg}
===========

\sum_p
\left(
\frac{k_p-k_p^{(0)}}{\sigma_p}
\right)^2
]

prevents poorly constrained parameters from becoming unphysical.

All three terms remain quadratic in (\mathbf k), so the total fit remains a single linear solve.

# 8. Another useful diagnostic: transform the QM Hessian into your valence-coordinate basis

Your model has internal coordinates

[
q=
{
\Delta r_{ij},
\Delta\theta_{ijk}
}.
]

Near equilibrium,

[
q=B,x,
]

where (B) is the Wilson matrix.

A general harmonic valence-force field is

[
E=\frac12q^TFq,
]

giving

[
H=B^T F B.
]

Your current independent bond-and-angle model assumes that (F) is diagonal or block-diagonal by parameter type:

[
F_{pq}\approx k_p\delta_{pq}.
]

Instead of immediately adding arbitrary interaction terms, reconstruct or project the QM Hessian into this internal-coordinate representation and inspect (F).

If the largest omitted elements are:

* bond–bond between two bonds sharing an atom, add stretch–stretch;
* bond–angle sharing an atom, add stretch–bend;
* angle–angle around the same atom, add bend–bend;
* interactions across a bond, add torsion or coplanar bending.

Full, partial, and internal Hessian fitting methods were developed precisely because direct Cartesian fitting and local parameter extraction emphasize different pieces of the QM force field. ([ResearchGate][3])

# 9. For silicon, your model corresponds naturally to a Keating-like VFF

For a tetrahedral Si network, I would use dot-product coordinates throughout.

For a bond,

[
q_{r,ij}
========

\frac{|\mathbf d_{ij}|^2-r_{0,ij}^2}{2r_{0,ij}}.
]

For small displacement,

[
q_{r,ij}\approx r_{ij}-r_{0,ij}.
]

For an angle,

[
c_{ijk}=
\hat{\mathbf d}*{ji}\cdot\hat{\mathbf d}*{jk},
\qquad
c_0=\cos\theta_0,
]

[
q_{\theta,ijk}
==============

\frac{c_{ijk}-c_0}{\sqrt{1-c_0^2}}.
]

For small distortion,

[
q_{\theta,ijk}\approx-\Delta\theta_{ijk}.
]

Then

[
E_0=
\frac12\sum_{\langle ij\rangle}
k_r q_{r,ij}^2
+
\frac12\sum_{\langle ijk\rangle}
k_\theta q_{\theta,ijk}^2.
]

No `acos`, `sin`, or `cos` evaluation is required once (c_0) and

[
s_0=\sqrt{1-c_0^2}
]

are stored.

If two parameters are insufficient, add the standard local couplings:

[
E_{rr}
======

\sum_{j}\sum_{i<k}
k_{rr},
q_{r,ji}q_{r,jk},
]

[
E_{r\theta}
===========

\sum_{ijk}
k_{r\theta}
(q_{r,ji}+q_{r,jk})q_{\theta,ijk},
]

and optionally

[
E_{\theta\theta}
================

\sum
k_{\theta\theta}
q_{\theta,a}q_{\theta,b}.
]

Modified valence-force models for diamond/zinc-blende materials use this same hierarchy: the original stretch and bend terms are often insufficient for full phonon dispersion, while stretch coupling, stretch–bend coupling, and coplanar bend coupling substantially improve it. 

# 10. What is “normally done”?

For molecular force fields, a common workflow is:

1. optimize the QM geometry;
2. calculate the QM Hessian;
3. extract or fit bond and angle constants from the Hessian;
4. fit torsional profiles separately;
5. optionally refine parameters against vibrational frequencies, conformational energies, forces, or experimental observables.

Seminario and modified-Seminario methods directly derive bond and angle constants from Hessian blocks; more elaborate approaches use partial or internal Hessian fitting. Modern force-field optimizers commonly combine multiple QM targets rather than relying on only one observable. ([ResearchGate][1])

For crystalline silicon, parameters are more often judged against:

* elastic constants;
* sound velocities;
* phonon dispersion along the Brillouin zone;
* zone-center and zone-boundary optical frequencies;
* Grüneisen parameters or strain dependence.

The original two-parameter Keating model is elegant and local, but extended VFF models are used when accurate dispersion over the entire Brillouin zone is required. ([APS Link][2])

# 11. My concrete recommendation for your implementation

Keep your analytic Hessian sensitivities. They are useful and make the problem linear.

Replace the default global objective by:

[
\boxed{
L =
\sum_{\rm systems}
\frac{1}{N_{\rm target}}
\sum_{s\in\mathcal T}
\frac{
\left|
\sum_p k_pD_pv_s-\lambda_sv_s
\right|^2
}{
\max(\lambda_s,\lambda_{\rm floor})^2
}
+
\beta L_{\rm localH}
+
\rho L_{\rm reg}.
}
]

This gives you:

* direct targeting of frequencies;
* preservation of eigenvectors;
* no mode matching during optimization;
* no repeated diagonalizations;
* a linear least-squares solve;
* explicit control over which spectral region matters;
* compatibility with multiple nanocrystal sizes.

Then use actual eigenvalue/eigenvector comparison only for validation and perhaps a final nonlinear refinement.

The previous full-Hessian fit is not wrong. It is simply solving a different question:

> “Which bond–angle Hessian is closest to the full QM Hessian element by element?”

Your real question is:

> “Which bond–angle force field reproduces the action of the QM dynamical matrix on the vibrational modes that matter?”

The second objective is better aligned with your reduced-model philosophy.

[1]: https://www.researchgate.net/profile/Hamed-Haghshenas-3/post/How-to-determine-the-internal-strain-tensor/attachment/59d62b63c49f478072e9d7a2/AS%3A273532569751569%401442226654912/download/Seminario_Hess2FF.pdf?utm_source=chatgpt.com "Calculation of intramolecular force fields from second- ..."
[2]: https://link.aps.org/doi/10.1103/PhysRev.145.637?utm_source=chatgpt.com "Effect of Invariance Requirements on the Elastic Strain Energy ..."
[3]: https://www.researchgate.net/publication/305951874_Partial_hessian_fitting_for_determining_force_constant_parameters_in_molecular_mechanics?utm_source=chatgpt.com "Partial hessian fitting for determining force constant ..."

---

# SWE-1.7



**Status report**

**What is implemented**
- Added Python [dihedral_energy_gradient](cci:1://file:///home/prokop/git/FireCore/tests/tSiNCs/test_FFfit.py:99:0-152:18) and [dihedral_hessian](cci:1://file:///home/prokop/git/FireCore/tests/tSiNCs/test_FFfit.py:155:0-169:12) that mirror [evalDihedral_Prokop](cci:1://file:///home/prokop/git/FireCore/cpp/common/molecular/UFF.h:918:4-1000:5) in [cpp/common/molecular/UFF.h](cci:7://file:///home/prokop/git/FireCore/cpp/common/molecular/UFF.h:0:0-0:0) (no `bSubNonBond` non-bond subtraction).
- Added [build_dihedrals](cci:1://file:///home/prokop/git/FireCore/tests/tSiNCs/test_FFfit.py:178:0-205:20), [dihedral_type_key](cci:1://file:///home/prokop/git/FireCore/tests/tSiNCs/test_FFfit.py:172:0-175:52), [compute_dihedral_sensitivity](cci:1://file:///home/prokop/git/FireCore/tests/tSiNCs/test_FFfit.py:208:0-227:12), and [add_dihedral_hessian](cci:1://file:///home/prokop/git/FireCore/tests/tSiNCs/test_FFfit.py:230:0-233:21).
- Wired dihedral topology, parameter type map, and per-type sensitivity `A_p = dH/dV` into the mode-basis fit and `H_model` evaluation.
- Added CLI flags `--dihedrals`, `--dihedral-n`, `--dihedral-d`, `--dihedral-h`.
- Code compiles and runs (`--dihedrals` does not crash).

**What the tests show**
- Without dihedrals: mode fit `relFrob` ~5-6% for Si, ~10-12% for C diamond.
- With `--dihedrals` (default `d=1`, `n=3`):
  - `k_angle` collapses to ~0.001-0.07 eV/rad².
  - `k_dihedral` becomes ~0.004-0.035 eV.
  - Fit quality degrades: `relFrob` ~12% for Si, ~50-60% for C diamond.
  - Many model frequencies become too low or spurious.

**Root-cause problem**
- `A_d = dH/dV` for a single torsion is **not positive semidefinite** for real geometry, even when `phi` is close to the intended minimum.
- Reason: [evalDihedral_Prokop](cci:1://file:///home/prokop/git/FireCore/cpp/common/molecular/UFF.h:918:4-1000:5) Hessian is `H = V [f'' b⊗b + f' C]`, where `C = ∂²phi/∂r²`. The `f' C` term is non-zero whenever `phi` is not exactly a minimum of the `cos(n*phi)` potential.
- For the short Si-H/C-H bonds, `C` has very large eigenvalues, so `f' C` produces strong negative eigenvalues in `A_d`.
- Fitting a positive `V` then adds negative curvature to the model, and `A_d` is also highly correlated with `A_angle` (so `k_angle` collapses).
- For C diamond, many dihedrals have `phi` ≈ 120°, i.e. they are at a **maximum** of `1 + cos(3*phi)` with `d=1`, so `A_d` is negative in the torsion direction.

**Options where I need your guidance**
1. **Use the exact `A_d` as now** and stabilize the fit with strong regularization/NNLS around the C++ LSQ solution. This will keep dihedral `V` small and may not improve much.
2. **Switch the dihedral sensitivity to the pure curvature term `f'' b⊗b`** (positive semidefinite). This is the part of the Hessian that controls the torsion spring stiffness; it avoids the `f' C` prestress/negative-curvature terms. The model then behaves like a stable torsion spring around the reference geometry.
3. **Choose `d` per dihedral** (or per type) to set the nearest `1, 2, or 3-fold` minimum to the actual `phi` (e.g. `d = -sign(cos(n*phi))`). This removes the worst-case negative signs, but the `f' C` term remains.

Which way do you want to go? I can do (2) or (3) quickly.