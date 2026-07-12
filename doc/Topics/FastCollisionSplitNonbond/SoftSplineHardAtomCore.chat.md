https://chatgpt.com/share/6a5354d1-7710-83eb-a8bd-c4e48716a869

---

# USER

I want to make efficient evaluation of non-bodning interactions for molecules which can work both for rigid-body symilation (where molecule have just position and rotation) and for flexible body (where molecule have springs on bonds and angles), generally I want to use semiimpicit euler or other spliting scheme where I split hard DOFs and soft DOFs and  solve hard does suing linear solver or some sub steps (they should be fast short range O(n)) and there should be also slow non-colvalen long range interaction (mostly electrostatics O(n^2)) which I solve which as long time step as possible.

I was thinking that best way how to do this is to match value and derivative of my pairwise potential between two atoms (which is Morse+Coulomb) by quadratic polynominal at some cutoff distance derived from the two atomic Radii (if the parameters of the atoms like Qi,Qj and Ri, Ei, Rj, Ej for Morse are all static during simulation we can easily precalculate the fitting parameters of the quadratic polynominal before calculation producing zero run-time overhed), then we can use some fast spatical datastructure like AABB bounding boxes to solve this interaction.

I also want to use similar methods for interaction with rigid substrate. There I have two ways how to prepresent the potential. 1) smooth B-spline (e.g. tri-cubic B-spline with quite large boxes like 2-5A voxel size) which describe slowly varying potential (mostly electrostatics in space - that is good for long range), then I have some atoms with short range potential describing "hard" repulsive core, there are reperested by spherical function and bound within cells like PIC or AABB bounding boxes, the total potential is sum of the to smooth+hard, where the hard is sampled only sparsely if the atom of absorbate molecule (the moving part) is close to the substrate atoms. 

Think it throu and make comprehesive nots, scritinize the approach identify chanlanges, suggest how to over come them, and generally suggst idealt way hot o implement it in organized way. I want start with rigid body dynamics of this, so we do not have to bother with complexity of interal springs (bonds angles) and can fully focus on collision and non-covalent interactions, but in later stages we want to ad this interanl dofs so whole desing should be modular engouh that it cann be plugged easily. 

 I give you some exsistin document and code where we dicsuss compoentns relevant for this approach.

after reading it plese ask me for clarigication, and point out problems


```
Region          Potential                    Solver
r < R_cut       E = 0.5*k1*(r-R_min)^2 + E_min   → Linear constraint (PD implicit solver)
R_cut ≤ r < R_cut2  E = 0.5*k2*(r-R_cut2)^2     → External force (explicit)
r ≥ R_cut2      E = 0                          → No interaction

where:
  k1 = -2*E_min / (d1*(d1+d2))    (convex, positive)
  k2 = -k1 * (d1/d2)              (concave, negative)
  d1 = R_cut - R_min,  d2 = R_cut2 - R_cut
```

---

# GPT 5.6 sol

I read both uploaded files: [RRsp3(3).cl](sandbox:/mnt/data/RRsp3%283%29.cl) and [ToDo_FastCollision_3(1).md](sandbox:/mnt/data/ToDo_FastCollision_3%281%29.md). The latter describes the harmonic-plus-smoothing split and its proposed integration into the PD neighbor lists. 

# Overall assessment

The **high-level decomposition is good**:

[
\boxed{
\text{smooth long-range field}
+
\text{sparse short-range atom cores}
+
\text{fast implicit/contact solver}
}
]

It is particularly attractive for rigid molecules and rigid substrates because:

* molecular geometry is immutable;
* atomic coordinates in the body frame can be precomputed;
* AABBs, bounding spheres, local BVHs, multipoles and inertia tensors belong to a reusable rigid shape;
* substrate spatial indexing is entirely static;
* slow electrostatics and fast repulsion genuinely have different spatial and temporal scales.

However, I would change the central mathematical interpretation:

> Do not represent collision repulsion as a temporary bilateral harmonic bond whose equilibrium length is (R_{\min}).

Instead, use:

[
\boxed{
\text{one-sided convex contact}
+
\text{explicit soft residual}
+
\text{separate long-range electrostatics}
}
]

For rigid molecules, the contact solver must operate on **six rigid-body degrees of freedom**, including torque. The present collision kernel is still a particle/sphere collision solver, despite its name.

# 1. What the proposed piecewise quadratic actually does

You define

[
U_1(r)=\frac12 k_1(r-R_{\min})^2+E_{\min},
\qquad r<R_{\rm cut},
]

and

[
U_2(r)=\frac12 k_2(r-R_{\rm cut2})^2,
\qquad
R_{\rm cut}\le r<R_{\rm cut2},
]

where

[
d_1=R_{\rm cut}-R_{\min},\qquad
d_2=R_{\rm cut2}-R_{\rm cut},
]

[
k_1=-\frac{2E_{\min}}{d_1(d_1+d_2)},
\qquad
k_2=-k_1\frac{d_1}{d_2}.
]

For (E_{\min}<0), this gives (k_1>0) and (k_2<0).

The construction is mathematically consistent in a limited sense:

* (U) is continuous at (R_{\rm cut});
* (U'), and therefore force, is continuous there;
* (U(R_{\rm cut2})=0);
* (U'(R_{\rm cut2})=0).

So it is a (C^1) compactly supported artificial pair potential.

## But it is not a fit to Morse + Coulomb

The formulas construct a new potential from

[
E_{\min},R_{\min},R_{\rm cut},R_{\rm cut2}.
]

They do not use the value, force or curvature of the original Morse+Coulomb potential at (R_{\rm cut}). Thus they generally do not match that potential except at the selected minimum.

A quadratic outer segment

[
U_2(r)=\frac12 k_2(r-R_2)^2
]

has only one free coefficient after (R_2) is chosen. If the original potential at (R_c) has value (U_c) and force (F_c=-U'_c), exact matching requires

[
k_2=\frac{F_c}{R_2-R_c},
]

and consequently

[
U_c=\frac12F_c(R_2-R_c).
]

Therefore

[
R_2-R_c=\frac{2U_c}{F_c}.
]

You cannot independently prescribe (R_2) and simultaneously match arbitrary (U_c) and (F_c) with that quadratic.

A cubic Hermite tail is almost as cheap and removes this restriction. With

[
t=\frac{r-R_c}{R_2-R_c},
]

a tail matching (U_c,U'_c) at (R_c), and zero value and derivative at (R_2), is

[
U_{\rm tail}(t)=
(2t^3-3t^2+1)U_c
+
(t^3-2t^2+t)(R_2-R_c)U'_c.
]

That requires a few FMAs in Horner form. Compared with neighbor search, square root, memory access and rigid-body torque accumulation, the difference between quadratic and cubic evaluation is negligible.

---

# 2. The most important conceptual problem: the inner spring is bilateral

Inside (R_{\rm cut}), the force is

[
F=-k_1(r-R_{\min}).
]

Therefore:

* when (r<R_{\min}), it repels;
* when (R_{\min}<r<R_{\rm cut}), it attracts toward (R_{\min}).

That is a temporary **bond**, not a collision constraint.

Whenever a pair enters the active set, the implicit solver tries to bring it to (R_{\min}). This can be acceptable if you deliberately want a strongly associating pair potential, but it has several undesirable consequences for collision handling:

1. Contact activation can make particles stick together.
2. Pairs repeatedly entering and leaving the list produce active-set hysteresis.
3. A newly activated attractive spring can introduce a large positional correction.
4. The meaning of the solver changes from “prevent penetration” to “implicitly minimize the whole pair well”.
5. Like-charged or otherwise unbound pairs may have no meaningful (E_{\min}<0).

The appropriate hard term is one-sided:

[
U_{\rm hard}(r)
===============

\frac12 k_h,[R_h-r]_+^2,
]

where

[
[x]_+=\max(0,x).
]

It has force

[
F_{\rm hard}(r)=
\begin{cases}
k_h(R_h-r), & r<R_h,\
0, & r\ge R_h.
\end{cases}
]

This term is:

* convex where active;
* purely repulsive;
* zero outside the core;
* suitable for an inequality/contact solver;
* independent of whether the total physical pair potential has an attractive minimum.

The attractive well and the moderate outer repulsion should remain in the explicit residual.

# 3. A cleaner exact splitting of the target potential

Let (U^*(r)) be the finite-range short-range target after any desired cutoff treatment. Define

[
U^*(r)=U_{\rm hard}(r)+U_{\rm residual}(r),
]

with

[
U_{\rm hard}(r)
===============

\frac12 k_h[R_h-r]_+^2,
]

[
U_{\rm residual}(r)
===================

U^*(r)-U_{\rm hard}(r).
]

This has several advantages.

The stiff positive curvature is captured by the implicit term. The residual contains:

* attraction;
* the moderate part of repulsion;
* the difference between the target core and the quadratic contact approximation.

You can choose (k_h) so that the residual curvature is much smaller than the original curvature over the physically sampled range.

A useful criterion is not simply an atomic radius. Choose (R_h) from the timestep you want:

[
U^{*\prime\prime}(R_h)
\approx
k_{\max},
\qquad
k_{\max}\sim c,\frac{m_{\rm eff}}{\Delta t_{\rm inner}^2}.
]

Below that distance, the target potential is too stiff for explicit integration and is transferred to the contact solver. This gives (R_h) a clear numerical meaning.

For an approximate force field, (U_{\rm residual}) can then be replaced by a short piecewise cubic or quintic polynomial stored per pair type.

## Pair-type precomputation

You do not need parameters per atom pair. For (N_t) atom types, store a table of (N_t(N_t+1)/2) entries:

```cpp
struct PairType {
    float R_hard;
    float R_switch;
    float R_outer;
    float compliance;   // or k_hard
    float4 coeff0;      // residual polynomial
    float4 coeff1;
};
```

This is more compact and cache-friendly than constructing parameters from (R_i,E_i,R_j,E_j,Q_i,Q_j) repeatedly.

---

# 4. Keep Coulomb outside the collision potential

Combining Morse and Coulomb before constructing the hard quadratic causes difficult cases:

* like charges can remove the attractive minimum;
* Coulomb can shift (R_{\min}) strongly;
* different charge products change whether the pair is bound;
* Coulomb does not vanish at (R_{\rm cut2});
* a charge-dependent hard-core stiffness is not physically desirable.

A better decomposition is

[
U_{ij}
======

U^{\rm Pauli/Morse}*{ij,\rm SR}
+
U^{\rm electrostatic}*{ij,\rm SR}
+
U^{\rm electrostatic}_{ij,\rm LR}.
]

A particularly clean split is

[
\frac1r
=======

\frac{\operatorname{erfc}(r/\sigma)}r
+
\frac{\operatorname{erf}(r/\sigma)}r.
]

Then:

* the `erfc` part is short-ranged and belongs in the fast neighbor-list calculation;
* the `erf` part is smooth and can be evaluated less often, represented by a grid, or approximated by body multipoles.

This same split works for both molecule–molecule and molecule–substrate interactions.

It also resolves the contradiction in the document between “no interaction beyond (R_{\rm cut2})” and the intended presence of long-range electrostatics.

---

# 5. Rigid-body contact must include rotation

For rigid body (A), let an atomic site have local coordinate (s_a). Its world position is

[
p_a=x_A+R(q_A)s_a.
]

For a site pair (a,b),

[
d=p_a-p_b,\qquad
r=|d|,\qquad
n=d/r.
]

Let

[
r_A=p_a-x_A,\qquad
r_B=p_b-x_B.
]

For a contact constraint

[
C=r-R_h\ge 0,
]

the effective inverse mass is

[
w_{\rm eff}
===========

\frac1{M_A}
+
\frac1{M_B}
+
(r_A\times n)^T I_A^{-1}(r_A\times n)
+
(r_B\times n)^T I_B^{-1}(r_B\times n).
]

This is essential: an atom near the end of a long molecule is easier to move through rotation than an atom close to the centre of mass.

For an XPBD-style compliant contact,

[
\tilde\alpha=\frac{\alpha}{\Delta t^2},
]

[
\Delta\lambda
=============

\frac{-C-\tilde\alpha\lambda}
{w_{\rm eff}+\tilde\alpha},
]

followed by unilateral projection

[
\lambda_{\rm new}=\max(0,\lambda+\Delta\lambda).
]

The resulting corrections are

[
\Delta x_A= M_A^{-1}n,\Delta\lambda,
]

[
\Delta\theta_A=
I_A^{-1}(r_A\times n)\Delta\lambda,
]

[
\Delta x_B=-M_B^{-1}n,\Delta\lambda,
]

[
\Delta\theta_B=
-I_B^{-1}(r_B\times n)\Delta\lambda.
]

XPBD is better suited than ordinary PBD when you want a stiffness or compliance that is not strongly dependent on timestep and iteration count. Projective Dynamics itself is not restricted to globally linear constraints; it uses specially structured energies and alternating local projection/global minimization. ([projectivedynamics.org][1]) XPBD explicitly extends position-based methods to compliant elastic and dissipative potentials. ([Matthias Research][2])

---

# 6. The current collision kernel is not a rigid-body collision solver

In `RRsp3(3).cl`, `compute_collision_cluster_rigid()` only writes

```cpp
dpos_coll
```

and `apply_corrections_rigid_ports()` reads

```cpp
float3 dx_coll = dpos_coll[i].xyz;
```

The quaternion correction is obtained only from

```cpp
drot_node
```

associated with the port constraints. There is no `drot_coll`.

Thus a collision at any point on the atomic sphere passes through the position represented by `pos[i]` and produces no collision torque. It is a sphere-particle correction, not a general rigid-molecule contact.

For one rigid pose per molecule, collision output needs at least:

```cpp
dpos_body[ibody]
drot_body[ibody]
```

or generalized contact impulses gathered onto the body.

---

# 7. Specific problems in `RRsp3(3).cl`

## 7.1 Fixed substrate atoms cannot collide

In `update_bboxes_rigid()`, an object is included only if

```cpp
invM > 1e-12f
```

At lines 589–593 of the collision kernel, atoms with zero inverse mass are replaced by zero position and zero radius. At lines 631–633, a candidate neighbor with zero inverse mass is skipped.

Consequently an immobile substrate atom with (M^{-1}=0):

* is absent from the AABB;
* is absent from the local collision geometry;
* cannot repel a moving atom.

Validity and mobility must be separate concepts.

A static atom should be loaded normally with

[
w_j=0,
]

so the moving object receives the full correction. Only its own output correction should be zero.

Use something like:

```cpp
bool valid = (gid < num_atoms) && (radius[gid] > 0.0f);
float invM = valid ? pos[gid].w : 0.0f;
```

Do not use `invM` as a validity flag.

## 7.2 `k_coll` is unused

The kernel receives

```cpp
const float k_coll
```

but never uses it. Instead it contains

```cpp
float dl = (rsum - dist) / w_tot;
dl *= 0.5f;
```

The factor (0.5) is a hard-coded relaxation factor.

For a pair constraint, both threads evaluating the pair does not by itself justify halving the endpoint correction. Each thread writes only its own endpoint. For equal masses, the exact Jacobi correction already moves each endpoint by half the penetration.

The extra factor means that one isolated equal-mass contact removes only half the overlap per iteration before the outer relaxation factor is even applied.

If gradual correction is intentional, call it relaxation and make it explicit. Preferably replace it by XPBD compliance.

## 7.3 Coincident particles remain unresolved

The condition

```cpp
d2 > 1e-16f
```

skips exactly or almost coincident centres. That is precisely where a fallback is required.

Possible deterministic normals are:

* previous-step separation;
* relative velocity direction;
* a hash-derived unit vector from the pair IDs;
* a cached contact normal.

Never silently ignore a maximally penetrated pair.

## 7.4 Silent ghost-list truncation

The code performs

```cpp
int total_ghosts = min(l_ghost_counter, MAX_GHOSTS);
```

Any candidates beyond `MAX_GHOSTS` disappear without an error flag. Dense molecular packing or large overlapping group boxes can therefore miss collisions nondeterministically.

At minimum store:

```cpp
ghost_overflow[group] = (l_ghost_counter > MAX_GHOSTS);
```

and fail loudly during development. A production solver can fall back to a larger global list or split the work.

## 7.5 AABBs are not yet a spatial data structure

`build_local_topology_rigid()` loops over every other group:

```cpp
for (int other_g = 0; other_g < num_groups; other_g++)
```

That remains

[
O(N_{\rm groups}^2).
]

The AABB test reduces narrow-phase work, but not broad-phase asymptotic complexity.

For a small number of groups this may be perfectly acceptable. For scaling, use:

* a uniform grid or hashed grid;
* sweep-and-prune;
* Morton sorting plus LBVH;
* a static/dynamic BVH split.

The current scheme also depends strongly on storage locality. If consecutive groups are not spatially coherent, their AABBs become large and the ghost lists explode.

## 7.6 No swept broad phase

The pair list is based on instantaneous geometry. With a long timestep, two bodies can cross the entire contact region between updates.

Expand body AABBs by a displacement bound such as

[
\Delta_{\rm sweep}
==================

|v|\Delta t
+
2R_{\rm body}\sin\left(\frac{|\omega|\Delta t}{2}\right)
+
R_{\rm outer}
+
R_{\rm skin}.
]

Alternatively use continuous collision detection only for unusually fast pairs.

## 7.7 Heavy-ball momentum is dangerous for changing contacts

At lines 1692–1705, previous positional corrections are mixed into the new move through `dpos_mom` and `beta`.

The same momentum buffer combines port and collision corrections. When a contact disappears, its previous correction can continue moving the body. With changing unilateral constraints this can create:

* artificial bouncing;
* sticking;
* penetration on the opposite side;
* instability at contact activation/deactivation.

I would not use heavy-ball acceleration on collision constraints initially. Use it only for fixed-topology internal constraints, or keep separate momentum buffers and reset the collision component when the active contact set changes.

# 8. Problems in the Markdown implementation plan

## Mixing bonds and collisions in one neighbor array

The plan proposes inserting dynamic collision pairs into the same `pd.neighs` and `pd.params` used for bonds.

This is undesirable because bonds and collisions differ fundamentally:

| Bonds                        | Contacts                                |
| ---------------------------- | --------------------------------------- |
| Fixed topology               | Dynamic topology                        |
| Bilateral                    | Unilateral                              |
| Usually always active        | Active only under penetration/proximity |
| Small known degree           | Potentially large degree                |
| Persistent stiffness         | Pair-dependent compliance/contact state |
| No pair creation/destruction | Frequent activation/deactivation        |

Use separate arrays:

```cpp
BondConstraint bonds[];
ContactConstraint contacts[];
```

and solve them in separate passes.

A fixed `nNeighMax` chosen for covalent valence is not safe for nonbonded contacts.

## Updating only the RHS may be insufficient

The document assumes that changing `pd.neighs` and `pd.params`, followed by `updatePD_dRHS()`, is enough.

That depends on the unprovided `ProjectiveDynamics_d` implementation. In ordinary PD, changing constraint topology or stiffness changes the global operator

[
A=\frac{M}{\Delta t^2}+\sum_c k_c A_c^TA_c.
]

At least the Jacobi diagonal changes. If any part of the operator, preconditioner or normalization is cached by `prepare_LinearSystem()`, it must be updated.

I cannot verify this point without `ProjectiveDynamics_d.h/.cpp`.

## Wrong ordering risk

The external smoothing force must enter the velocity kick or predicted position before the constraint solve. Calling the neighbor evaluation from a generic `user_update` callback is safe only if its exact placement is known.

The correct order is conceptually:

1. evaluate/cache explicit forces;
2. kick velocities;
3. predict positions and orientations;
4. build contacts using predicted poses;
5. solve contacts;
6. reconstruct or correct velocities.

Building contacts at old positions and then predicting with a long timestep encourages tunnelling.

---

# 9. Recommended rigid-body architecture

## Static shape data

Each molecular species should have immutable shape data:

```cpp
struct RigidShape {
    int nSites;
    const float4* localPosType; // xyz + type
    float3 localAabbCenter;
    float3 localAabbHalf;
    float boundRadius;
    float invIbody[9];

    // Optional for larger molecules
    BVHNode* localBVH;
};
```

For repeated identical molecules, all instances reference the same shape. Only pose, velocity and body-specific state are duplicated.

## Dynamic body state

```cpp
struct RigidBody {
    float4 posInvMass;
    float4 quat;
    float4 velocity;
    float4 omega;
    int shapeId;
};
```

## Broad-phase body pairs

Construct candidate body pairs from body-level swept AABBs.

A rotated local AABB can be transformed cheaply:

[
h_{\rm world}=|R(q)|h_{\rm local}.
]

Then

[
AABB=[x+c_{\rm world}-h_{\rm world},,
x+c_{\rm world}+h_{\rm world}].
]

Do not put individual molecular atoms directly into the global broad phase unless the number of bodies is tiny.

## Narrow phase

For each candidate body pair:

* for molecules below perhaps a few dozen sites, brute-force site pairs may be fastest;
* for larger molecules, traverse precomputed local sphere trees/BVHs;
* emit atom-site contacts only for (r<R_{\rm hard}+{\rm skin});
* evaluate soft explicit pairs for (r<R_{\rm outer}).

After measurement, reduce many site contacts between the same body pair to a small contact manifold, perhaps the deepest 2–4 spatially separated contacts. Otherwise two large molecules can produce dozens of redundant atomic contacts and become artificially overconstrained.

## Force and torque accumulation

Every explicit atom-pair force must produce

[
F_A \mathrel{+}=f,
\qquad
\tau_A\mathrel{+}=r_A\times f,
]

[
F_B \mathrel{-}=f,
\qquad
\tau_B\mathrel{-}=r_B\times f.
]

The same rule applies to the substrate field:

[
F_{\rm body}=\sum_a f_a,
\qquad
\tau_{\rm body}=\sum_a r_a\times f_a.
]

# 10. Avoiding atomics on the GPU

A clean Jacobi implementation can use a contact-centric two-pass scheme.

### Pass 1: one thread per contact

Read both body poses and compute:

```cpp
deltaLambda[ic]
endpointA[ic] = {dxA, dthetaA}
endpointB[ic] = {dxB, dthetaB}
```

No body is written.

### Pass 2: one thread or workgroup per body

Use a CSR list of incident contact IDs and gather all corrections:

```cpp
for(contact : bodyContacts[ibody]){
    dx     += signedEndpointTranslation(contact, ibody);
    dtheta += signedEndpointRotation(contact, ibody);
}
```

Then write the body exactly once.

This avoids atomics and makes Jacobi semantics explicit. It requires two kernels per iteration, but the memory pattern is predictable.

A simpler first implementation can duplicate contact evaluation in each body’s adjacency list, as the current particle-centric kernel does, but then per-contact XPBD multipliers cannot be shared cleanly.

# 11. Rigid substrate representation

Your substrate idea is good, with one crucial refinement: the coarse grid must store only a genuinely smooth field.

Use

[
E_{\rm substrate}
=================

\sum_a q_a\phi_{\rm smooth}(p_a)
+
\sum_{a,J\ {\rm nearby}}
U^{\rm core}_{aJ}(|p_a-R_J|).
]

A strong version is

[
U^{\rm core}_{aJ}
=================

U^{\rm Pauli/Morse}_{aJ}
+
k_eq_aQ_J\frac{\operatorname{erfc}(r/\sigma)}r,
]

while the grid stores the complementary long-range electrostatic potential generated by

[
k_eQ_J\frac{\operatorname{erf}(r/\sigma)}r.
]

Then there is no ambiguous double counting.

## B-spline grid concerns

A voxel spacing of (2)–(5) Å is plausible only after the atomic high-frequency components have been removed. It is too coarse for raw ionic or point-charge potentials near substrate atoms.

Also, tricubic B-spline interpolation normally touches (4^3=64) coefficients per sample. On a GPU, a denser trilinear grid with eight reads may outperform a very coarse cubic grid. This should be benchmarked rather than assumed.

Whichever interpolation is chosen:

> Store/interpolate the potential and obtain force from the analytic derivative of the same interpolant.

Do not independently interpolate potential and force grids, because they will generally cease to be exactly conservative.

## Static core cell list

For an approximately planar substrate, a static 2D cell list is sufficient:

* cell size near the maximum core cutoff;
* each cell stores a compact contiguous list of substrate atom indices;
* each adsorbate atom checks nearby (x,y) cells;
* reject atoms by (z) and exact distance.

This is simpler than a dynamic AABB hierarchy and has essentially no update cost.

# 12. Multiple timestep integration

A sensible hierarchy is:

[
F=F_{\rm contact}+F_{\rm near}+F_{\rm far}.
]

* `contact`: hard one-sided constraints, several solver iterations;
* `near`: residual Morse/dispersion and short-range `erfc` electrostatics;
* `far`: smooth `erf` electrostatics or substrate grid.

For the first implementation I would use:

1. slow-force half kick;
2. several inner rigid-body steps;
3. recompute slow force;
4. slow-force half kick.

Inside each inner step:

1. evaluate near explicit forces and torques;
2. semi-implicit velocity kick;
3. predict pose;
4. build contacts from predicted pose;
5. perform XPBD contact iterations;
6. update velocities and angular velocities from corrected poses.

This is closely related to reversible multiple-timestep splitting such as r-RESPA, which was designed for systems with separated force timescales. ([Columbia University][3])

But the outer timestep cannot be made arbitrarily large. Multiple-timestep molecular dynamics can suffer resonances when the outer update couples to unresolved faster motion; this is a known limitation even when the slow force itself is smooth. ([arXiv][4])

Therefore add diagnostics based on:

* maximum displacement since the slow force was evaluated;
* maximum orientation change;
* relative change of slow force;
* energy drift in conservative tests.

# 13. Extension to flexible molecules

The modular boundary should not be “rigid solver versus flexible solver”. It should be the mapping from an interaction endpoint to generalized degrees of freedom.

For a flexible atom:

[
p_a=x_a,
]

and the Jacobian contains only translation.

For a rigid-body site:

[
p_a=x_A+R(q_A)s_a,
]

and the Jacobian contains translation and rotation.

Thus the same pair interaction can provide:

```cpp
PairGeometry {
    endpointA;
    endpointB;
    normal;
    distance;
    pairType;
}
```

while two solver variants map its force or constraint to:

* atom translations;
* rigid-body translation and rotation.

On the GPU I would not implement this with virtual classes. Use shared inline mathematical functions and separate specialized kernels.

Internal bonds and angles should remain separate from dynamic contacts:

```text
fixed bond constraints
fixed angle constraints
dynamic nonbonded contacts
explicit near forces
slow long-range forces
```

That makes it easy to introduce flexible molecules later without redesigning the nonbonded infrastructure.

# 14. Suggested implementation sequence

## Stage 1: correct rigid-body reference

Implement a double-precision CPU all-pairs reference with:

* one pose per molecule;
* force and torque from atom pairs;
* one-sided rigid-body contacts;
* no AABB;
* no multiple timesteps;
* no B-spline grid.

This establishes the equations.

## Stage 2: GPU rigid bodies

Implement:

* body state;
* transformed atomic sites;
* body-pair broad phase;
* atom-pair narrow phase;
* explicit force/torque accumulation;
* rigid XPBD contact.

Initially use a simple broad phase even if (O(N^2)).

## Stage 3: spatial acceleration

Add:

* swept body AABBs;
* grid/BVH broad phase;
* neighbor-list skin;
* overflow detection;
* contact manifold reduction.

## Stage 4: substrate decomposition

Add:

* smooth electrostatic grid;
* static substrate cell list;
* explicit atom-centred cores;
* direct-vs-split energy and force error tests.

## Stage 5: multiple timesteps

Only after the single-rate solver is stable, introduce:

* short-range inner steps;
* slow electrostatic outer steps;
* symmetric force kicks;
* drift and resonance diagnostics.

## Stage 6: flexible molecules

Replace rigid sites by atomic DOFs for selected bodies and add fixed-topology bond/angle constraint passes.

# 15. Essential tests

Before realistic molecules, I would require these tests to pass:

1. **Two equal spheres:** a single overlap is corrected symmetrically.
2. **Movable sphere against fixed sphere:** the static sphere participates but does not move.
3. **Off-centre rigid contact:** translation and rotation match the rigid effective-mass formula.
4. **Central rigid contact:** torque is numerically zero.
5. **Cutoff sweep:** energy and force are continuous across all piecewise boundaries.
6. **Pair without an attractive minimum:** parameters remain valid.
7. **Coincident centres:** deterministic separation without NaNs.
8. **Dense packing:** ghost/contact overflow is reported, never silently ignored.
9. **Angular neighbor-list skin:** rotation alone triggers rebuilding when necessary.
10. **Substrate decomposition:** grid plus core reproduces direct substrate energy and force.
11. **Conservative trajectory:** no damping, friction or restitution; quantify energy drift.
12. **MTS ratio scan:** compare outer ratios (1,2,4,8,\ldots) against the single-rate trajectory.

# Clarifications needed before fixing the design more concretely

1. By **rigid molecule**, do you mean exactly one position, quaternion, velocity and angular velocity for the whole molecule? The uploaded kernel appears closer to oriented atomic nodes connected through ports.

2. Is your primary objective physically meaningful molecular dynamics and approximate energy conservation, or rapid stable relaxation/collision avoidance? The former favours force-residual IMEX integration; the latter allows much harder XPBD projection.

3. Should nonbonded attraction actually bind molecules into minima, or is the immediate goal mostly steric repulsion plus electrostatics?

4. What are the typical numbers of molecules and atoms per molecule? A brute-force atom-pair narrow phase is often preferable for small rigid molecules, even when a sophisticated hierarchy sounds more elegant.

5. Are substrate atoms represented with zero inverse mass in the current arrays? If so, the existing collision kernel definitely excludes exactly the objects you need to collide with.

6. Can the electrostatics be split into smooth long-range and compact short-range components, or do you require the exact direct (q_iq_j/r) sum to remain one indivisible calculation?

7. Could you provide `ProjectiveDynamics_d.h/.cpp` and the host-side sequence that launches these kernels? The claim that only `updatePD_dRHS()` must change cannot be verified from the two current files.

[1]: https://www.projectivedynamics.org/?utm_source=chatgpt.com "Projective Dynamics: Fusing Constraint Projections for Fast ..."
[2]: https://matthias-research.github.io/pages/publications/XPBD.pdf?utm_source=chatgpt.com "XPBD: Position-Based Simulation of Compliant ..."
[3]: https://www.columbia.edu/cu/chemistry/groups/berne/papers/jcp_97_1990_1992.pdf?utm_source=chatgpt.com "Reversible multiple time scale molecular dynamics"
[4]: https://arxiv.org/abs/1307.1167?utm_source=chatgpt.com "Stochastic resonance-free multiple time-step algorithm for molecular dynamics with very large time steps"


---

# USER

OK, I should explain / clarify multiple thigns:
1) I do not thinkg my hard-core parabolic potential which "stick particles tograter" because by cutoff will be always deep in repulsion so the equlibirum distance R0 (minimum of the parabola) will be always outside the hard-core region and within the hardcore region it will be always repulsive
2) I do not worry abut the toques and rotational dynamics, this is different problem. When we solve interaction forces between the atoms we can easily compute teorque on the COM of the molecule by sanrad formula coss(r,f)
3) I do not like you idea of using erf(r)/r for electrostatics. It is unnecesarily compicated and we do not need, that matching using quadratic polynominal. If we take 1/r and at some Rc we cut it and replace by quadratic polynominal for r<Rc s solves the problem with divergence with minimal performence cost and no hessle, we do not need any costly erf(r)/r.
4) I do not understand why you do not like fittong total pairwise potential Morse+Coulomb with single quadratic hard-core polynominal. As I say if parameters of both Coulomb and Morse are fixed (Atom types and chanrges does not change during simulation) we can store the coeficients of theat polynominal precomputed wthout any performance cost, and it will be much faster than compute them separately.   

5) I think you missunderstood several things im my original desing: The switching distance Rc is not the Morse minimum R0. I think isntead we should put Rc alrady inside the Morse repulsion. The fast that for some harges there my be no minimum is fine. We should perhaps anyway use clamped electrostatics with some radional clamping function 

```
def soft_clamp(y, y1, y2):
    """Rational saturating clamp: y free below y1, asymptotes to y2 above y1.
    z=(y-y1)/(y2-y1), return y1+(y2-y1)*z/(1+z)"""
    y12 = y2 - y1
    z = (y - y1) / y12
    return np.where(y > y1, y1 + y12 * z / (1.0 + z), y)
```

5) asymptitic scaling of BVH is not so importaint. We are working with small system (<1000 atoms, often <100) so we do not care aboyt asymptotic but about pre factors. GPU prefer simple architectrue simple algorithm with predictable branching or branche less, therefore I do not want to use BVH with pointers and stuff. Our system simulate often the whole thing within one kernel, within ne workgroup (different workgorup simulate different replicas of the system, we have ofthen like 1000 independnet replicas of the system on single GPU in paralle, this is different use case than e.g. large-scale simulations of proteins with millions of atoms parllelized on multiple GPUs). Thereofre I prefer to have only single spatial-partitioning level for simplicity

add you questions:
1) yes by rigid molecule I mean y one position, quaternion, velocity and angular velocity for the whole molecule. But the code I geve you in RRsp3.cl is not the actually code which I want to make, it is just reference showing some aspects of it. Do not spend to much time criticising that particular piece of code.
2) My primary objective is rapid stable relaxation/collision avoidance, but in the end it would be good to derive both flavors, as this method can be used for both 
3) I think the non-bonding minim (at R0) will be in the "soft" region, so they will be solved with the long-time-step O(n^2) forcefield, but the short-time-step (O(n)) solver should deal only with hard-core collision for r<Rc<R0
4) My idea is to have dozen of molecule each with doezn of tom in one system ( e.g. 16 molecule each with 16 atom ), or it can be polymer units in one polymer. Sumbstrate can have hundrets or thousands of taoms but we will not explicitly represrnt all, only the top layer with some fitter/projected (we call it Folded atomic functions, but that is another storry)
5) Ithink you goo into to much details of the .cl code I gave you, I want to make new code from scrath. substrate atoms will be not presnt in the dyamics, they will be as external forcefield evauated separately, but using the same trick to split soft-long-range and hard-short-range non-boning interactions
6) no definitely I do not want shard 1/r divergence at zero, I will damp it anyway. But I also do not want the erf/refc sum, because I think the derivative mathcing polynominal is simple and faster solution
7) yes I agree that if we move atoms between the core and outer region, the collision neighbors change so the matrix change, this wuld be prblem for Cholesky LDL^T cactorization, but for Jacobi it is perhaps not problem. If we do position-based approach or projective dynamics we do JAcobi, that is good for relaxation (remove energy high stability) if we do proper energy-conserving molecular dynamics we use impulse based method (like sequantila impulses used in box2d/box3d resp. force time-steps)

What you think now? Can you update your critique and come with some new ideas after this clarification?

---

# GPT 5.6 sol

Yes. Your clarification removes two of my main objections. I had interpreted (R_0) as lying inside the active quadratic region, while you intend

[
r<R_c<R_0,
]

so the quadratic is sampled only on its repulsive left branch. It therefore does **not** bind or stick the particles. I also agree that an Ewald-like `erf/erfc` decomposition is unnecessary for your small-system, many-replica GPU use case.

The basic idea is sound:

[
\boxed{
\text{expensive soft total potential at long timestep}
+
\text{precomputed quadratic core at short timestep}
}
]

The main remaining issue is not whether the quadratic works, but how to split it mathematically between the slow force evaluation and the fast implicit/contact solver without double counting or timestep-dependent switching errors.

---

# 1. The core quadratic is valid

Let the physical outer pair potential be

[
U_{\rm out}(r)
==============

U_{\rm Morse}(r)+U_{\rm Coul}(r),
]

possibly with whatever damping you decide to use.

Choose (R_c) sufficiently far inside the repulsive branch. At (R_c), precompute

[
U_c=U_{\rm out}(R_c),
]

[
g_c=U'_{\rm out}(R_c)<0,
]

and choose a positive core curvature (k>0).

Then define, for (r<R_c),

[
U_{\rm core}(r)
===============

U_c+g_c(r-R_c)
+\frac12 k(r-R_c)^2.
]

This automatically matches value and derivative:

[
U_{\rm core}(R_c)=U_c,
\qquad
U'_{\rm core}(R_c)=g_c.
]

Its radial force magnitude is

[
F_{\rm core}(r)
===============

# -U'_{\rm core}(r)

F_c+k(R_c-r),
]

where

[
F_c=-g_c>0.
]

Therefore it is repulsive throughout (r<R_c), with increasing force under increasing penetration.

The extrapolated minimum of the parabola is

[
R_q=R_c+\frac{F_c}{k}>R_c.
]

That minimum is outside the region in which the parabola is used, exactly as you intended.

## Better parameterization

I would not store the parabola as ((R_{\min},E_{\min},k)), because this invites confusion with the physical Morse minimum. Store directly:

```cpp
struct PairCore {
    float Rc;
    float Rc2;
    float Fc;   // outward force at Rc
    float k;    // additional stiffness under penetration
    float Uc;   // needed only if energy is evaluated
};
```

Then the inner force is simply

```cpp
Fmag = Fc + k*(Rc-r);
```

The fictitious quadratic minimum (R_q) never needs to be calculated at runtime.

---

# 2. The unresolved choice is the curvature (k)

Matching value and derivative uses only two conditions. A general quadratic has three coefficients, so one parameter remains free.

There are three sensible ways to choose (k).

## A. Match the curvature at (R_c)

[
k=U''_{\rm out}(R_c).
]

Then the core is exactly the second-order Taylor approximation of the total Morse+Coulomb potential at (R_c).

This is the most direct mathematical fit, provided

[
U''_{\rm out}(R_c)>0.
]

It is probably the best first implementation.

## B. Fit an additional inner force value

Choose an inner reference distance

[
R_a<R_c
]

and a desired repulsive force (F_a). Then

[
k=\frac{F_a-F_c}{R_c-R_a}.
]

This gives direct control over how hard the core becomes under severe overlap.

## C. Choose stiffness from acceptable penetration

Suppose the largest expected compressive load is (F_{\max}), and the acceptable penetration is (\delta_{\max}). Then approximately

[
k\ge
\frac{F_{\max}-F_c}{\delta_{\max}}.
]

This may be more useful for a relaxation/collision solver than exact curvature matching.

I would support A and C:

```text
k = max( U''(Rc), k_required_for_penetration )
```

possibly with an upper clamp to avoid an unnecessarily ill-conditioned solver.

---

# 3. Fitting total Morse+Coulomb is completely reasonable

I agree with you here. If charges and atomic parameters do not change, fitting the **total pair potential** is simpler than separately evaluating Morse and regularized Coulomb in the core.

For each fixed pair:

1. evaluate total (U_c);
2. evaluate total (F_c=-U'(R_c));
3. choose or evaluate (k);
4. store the resulting core parameters.

Inside (R_c), no exponential and no division are needed apart from the distance normalization already required to obtain the force direction.

The only required condition is

[
F_c>0.
]

If an oppositely charged pair has such strong Coulomb attraction that the total force is still attractive at the proposed (R_c), move (R_c) deeper into the Morse repulsion until the total force is unambiguously outward.

Thus the condition for selecting (R_c) should be something like

[
F_{\rm total}(R_c)>F_{\rm safety}>0,
]

not merely a geometric combination of atomic radii.

The radii can provide the initial guess, but the total pair force should validate it.

---

# 4. You probably do not need a separate Coulomb clamp

If the raw Morse+Coulomb potential is evaluated only for

[
r\ge R_c,
]

and the total quadratic replaces it for

[
r<R_c,
]

then the (1/r) singularity is never evaluated. The quadratic itself is the Coulomb regularization.

That is simpler than applying `soft_clamp()` to the Coulomb part and then fitting the combined result.

Your proposed clamp also has one asymmetry:

```python
np.where(y > y1, ..., y)
```

It clamps only large positive values. For opposite charges,

[
y=\frac{q_iq_j}{r}\rightarrow-\infty,
]

so the attractive divergence remains unclamped.

A two-sided clamp would need to act on (|y|) and restore the sign. Also, if the energy is clamped, the force must include the chain-rule derivative

[
\frac{d\tilde y}{dy}
====================

\frac{1}{(1+z)^2}.
]

Otherwise energy and force become inconsistent.

My recommendation is therefore:

[
\boxed{
\text{Do not clamp Coulomb separately; replace the total potential below }R_c.
}
]

A separate clamp is useful only as a defensive fallback.

---

# 5. The most useful new decomposition

Your total inner quadratic can be separated exactly into a soft constant-slope component and a hard quadratic penetration component.

Let

[
x=r-R_c.
]

Inside the core,

[
U_{\rm core}
============

U_c+g_cx+\frac12kx^2.
]

Define the soft continuation

[
U_{\rm soft}(r)=
\begin{cases}
U_{\rm out}(r), & r\ge R_c,[2mm]
U_c+g_c(r-R_c), & r<R_c,
\end{cases}
]

and the hard component

[
U_{\rm hard}(r)
===============

\frac12 k[R_c-r]_+^2.
]

Then

[
U_{\rm total}=U_{\rm soft}+U_{\rm hard}
]

reproduces exactly:

* the original Morse+Coulomb potential outside (R_c);
* your fitted quadratic inside (R_c).

This is important because the two terms naturally belong to different solvers:

[
\boxed{
U_{\rm hard}
============

\frac12 k[R_c-r]_+^2
\quad\rightarrow\quad
\text{fast implicit/Jacobi/contact solver}
}
]

[
\boxed{
U_{\rm soft}
\quad\rightarrow\quad
\text{ordinary explicit force evaluation}
}
]

The implicit solver no longer tries to project the pair toward the extrapolated minimum (R_q>R_c). It only removes penetration toward the boundary (R_c). The constant repulsive force (F_c) is supplied by the soft component.

This is numerically cleaner than giving the PD solver a spring whose minimum lies outside its active domain.

---

# 6. A more rigorous multiple-timestep split without `erf`

There is still a subtle problem if the complete soft force is evaluated only once every several inner steps.

Suppose a pair is outside (R_c) during the slow-force evaluation and enters the core before the next slow update. Its cached force corresponds to its old distance. The contact solver fixes the severe penetration, but the soft contribution is temporarily stale.

For rapid relaxation this may be acceptable. For better dynamics, use a compact polynomial distance split.

Choose

[
R_c<R_1<R_2<R_0,
]

so the entire transition remains on the repulsive side of the physical minimum.

Define a smooth switching function (s(r)):

[
s=1,\quad r\le R_1,
]

[
s=0,\quad r\ge R_2.
]

For example, with

[
t=\frac{r-R_1}{R_2-R_1},
]

a cubic (C^1) switch is

[
s(t)=1-3t^2+2t^3.
]

Then split

[
U_{\rm near}=s(r)U_{\rm soft}(r),
]

[
U_{\rm far}=[1-s(r)]U_{\rm soft}(r).
]

The complete potential is

[
U=
U_{\rm hard}
+
U_{\rm near}
+
U_{\rm far}.
]

This gives:

* (U_{\rm hard}): implicit collision solver;
* (U_{\rm near}): compact short-range explicit force, evaluated at every inner step;
* (U_{\rm far}): slowly varying all-pairs force, evaluated at the long timestep.

Inside (R_1),

[
U_{\rm far}=0,
]

so no stale long-range force acts on a deeply colliding pair.

Outside (R_2),

[
U_{\rm near}=0,
\qquad
U_{\rm far}=U_{\rm out}.
]

Around the physical nonbonded minimum (R_0), the interaction therefore belongs entirely to the slow all-pairs calculation, as you wanted.

This is the polynomial equivalent of a range separation, but without expensive special functions.

## Performance refinement

The near component can itself be pre-fitted by a low-order polynomial over ([R_c,R_2]). Then the inner kernel needs only:

* distance;
* a few FMAs;
* multiplication by the normalized displacement.

No Morse exponential or Coulomb division is required in the inner loop.

---

# 7. Relaxation and energy-conserving dynamics should be separate modes

I agree that they should share the potential representation but use different solvers.

## Mode A: rapid relaxation / collision removal

Use soft forces to predict positions:

[
x^*
===

x_n+\Delta t,v_n
+\Delta t^2M^{-1}F_{\rm soft}.
]

Then approximately minimize

[
E(x)
====

\frac{1}{2\Delta t^2}
(x-x^*)^TM(x-x^*)
+
\sum_c
\frac12 k_c[R_{c}-r_c(x)]_+^2.
]

Solve using a few Jacobi/PBD/XPBD iterations.

This is deliberately dissipative and highly stable. It is appropriate for:

* structure relaxation;
* collision avoidance;
* basin exploration;
* preventing catastrophic overlaps during large optimization steps.

Dynamic topology is not a serious problem for Jacobi. You only need to recompute:

* active contacts;
* contact directions;
* effective inverse masses or diagonal weights.

No factorization and no persistent global matrix are needed.

## Mode B: approximately conservative molecular dynamics

PBD projection should not be used as the main conservative integrator because positional projection generally removes or adds energy.

For finite quadratic cores, use:

* velocity Verlet with hard-core substeps;
* implicit midpoint for the hard pair coordinate;
* or symmetric pairwise operator splitting.

For an ideal hard-wall approximation, a frictionless sequential impulse with restitution (e=1) conserves kinetic energy for an isolated collision. But that represents a hard contact, not the finite quadratic potential.

So I would distinguish:

```text
RELAX mode:
    Jacobi/XPBD positional correction

DYNAMICS mode:
    symmetric force splitting + hard-core substeps
    or elastic impulse limit
```

Both modes use the same (R_c,F_c,k) tables.

---

# 8. The quadratic is not fully linear in Cartesian coordinates

Precomputing the coefficients removes the expensive evaluation of Morse and Coulomb, but the pair energy

[
\frac12k(r-R_q)^2
]

is still nonlinear in Cartesian coordinates because

[
r=|\mathbf r_i-\mathbf r_j|
]

and the contact direction changes.

That is not a serious obstacle, but it means the method is not a single exact linear solve with a precomputed matrix.

For Jacobi/PBD:

1. recompute the current direction;
2. calculate scalar penetration or scalar force;
3. distribute the correction;
4. iterate.

This is exactly the kind of geometry-dependent local iteration that is cheap for a small system in one GPU workgroup.

---

# 9. Spatial organization for your actual system size

Given approximately

* 16 molecules;
* 16 atoms per molecule;
* one replica per workgroup;
* hundreds or thousands of replicas;

I agree that a BVH is probably the wrong architecture.

## Best first broad phase: molecule AABB plus dense atom pairs

For every rigid molecule:

1. transform its atoms into world coordinates;
2. compute one molecular AABB;
3. test all molecule pairs;
4. for overlapping molecular AABBs, evaluate all atom pairs.

With 16 molecules, there are only

[
\frac{16\cdot15}{2}=120
]

molecular pair tests.

An overlapping pair has

[
16\times16=256
]

atomic distance tests.

Even testing all molecular pairs fully is only about

[
120\times256\approx 3.1\times10^4
]

atom-pair candidate checks per replica. That is very manageable, especially because:

* data can remain in local memory;
* loops have fixed bounds;
* there are no pointers;
* branches are predictable;
* many replicas provide GPU occupancy.

A single molecular AABB layer is therefore a very reasonable design.

## For somewhat larger systems

For hundreds to roughly one thousand dynamic atoms, use a single uniform grid:

```text
cellCount[nCells]
cellStart[nCells]
cellAtoms[nAtoms]
```

Build it within the workgroup using:

1. count atoms per cell;
2. workgroup prefix sum;
3. scatter atom indices;
4. inspect neighboring cells.

No tree and no pointer chasing.

## Contact-list reuse

Build a contact candidate list using an inflated cutoff

[
R_{\rm list}=R_c+R_{\rm skin}.
]

Reuse it for several Jacobi iterations or inner substeps. Rebuild when translational plus rotational site displacement exceeds roughly

[
R_{\rm skin}/2.
]

For your molecule sizes, a fixed-capacity local-memory list with an explicit overflow flag should be adequate.

---

# 10. Force and torque accumulation

I agree that rotation is conceptually separate and straightforward.

For an atom-site force (\mathbf f_a) on a rigid molecule,

[
\mathbf F_{\rm body}
====================

\sum_a \mathbf f_a,
]

[
\boldsymbol\tau_{\rm body}
==========================

\sum_a
(\mathbf r_a-\mathbf R_{\rm COM})\times\mathbf f_a.
]

The same atom-pair kernel can therefore accumulate:

```cpp
bodyForce[ibody] += f;
bodyTorque[ibody] += cross(atomPos-COM, f);
```

For positional Jacobi corrections, the corresponding rigid-body rotational correction requires the inverse inertia tensor, but that does not affect the pair-potential representation.

---

# 11. Substrate decomposition

Your substrate representation is consistent with the same architecture:

[
V_{\rm substrate}
=================

V_{\rm smooth\ grid}
+
\sum_{J\ {\rm nearby}}V^{\rm hard}_J.
]

The most important implementation rule is:

> Fit the smooth grid to the residual after subtracting the atom-centred hard functions.

That is,

[
V_{\rm smooth}(\mathbf r)
=========================

## V_{\rm reference}(\mathbf r)

\sum_JV^{\rm hard}_J(|\mathbf r-\mathbf R_J|).
]

Otherwise the smooth grid and explicit cores will double count part of the repulsion.

For a planar substrate, use a static two-dimensional cell grid:

```text
cellStart[ix,iy]
cellCount[ix,iy]
substrateSiteIndices[]
```

Each moving atom checks only nearby substrate cells. The substrate sites are never part of the dynamic solver.

The smooth grid can be evaluated at the long timestep. The hard atomic functions use the same quadratic-core machinery as molecule–molecule interactions.

---

# 12. Proposed modular organization

```text
PairPotentialTable
    outer Morse+Coulomb parameters
    Rc
    Fc
    k
    Uc
    near/far switch parameters

RigidShape
    local atomic positions
    atom types
    local AABB
    mass and body inertia

RigidState
    position
    quaternion
    linear velocity
    angular velocity
    force
    torque

BroadPhase
    molecular AABBs
    optional uniform cell grid

NearInteraction
    compact explicit near polynomial
    active hard-core contacts

FarInteraction
    O(N²) slow Morse+Coulomb / long-range force

HardSolverRelax
    Jacobi / XPBD

HardSolverDynamics
    substepped force or impulse solver

SubstrateSmooth
    B-spline or other grid

SubstrateHard
    static cell list of atom-centred functions
```

When flexible molecules are introduced, only the endpoint accumulation changes:

* rigid site force → body force and torque;
* flexible atom force → atom force directly.

The pair tables, distance splitting, broad phase and substrate representation can remain unchanged.

---

# 13. What I would implement first

1. Precompute (U_c,F_c,k) for every pair type or fixed atomic pair.
2. Replace total Morse+Coulomb below (R_c) by the matched quadratic.
3. Implement the exact decomposition

[
U_{\rm core}
============

\underbrace{U_c+g_c(r-R_c)}*{\text{explicit}}
+
\underbrace{\frac12k(R_c-r)^2}*{\text{implicit}}.
]

4. Use molecular AABBs and dense atom-pair loops.
5. Implement a dissipative Jacobi/XPBD relaxation mode first.
6. Add the static substrate cell list using identical hard-core parameters.
7. Add the polynomial near/far split only after the single-rate version is stable.
8. Finally add the conservative/substepped dynamics mode.

The strongest conclusion after your clarification is:

[
\boxed{
\text{Your matched total quadratic is good, but the implicit solver should handle only its curvature term.}
}
]

The boundary force (F_c) and the attractive/soft region should remain explicit. This gives a clean contact constraint at (R_c), avoids projecting toward a fictitious distant minimum, and provides a natural route to both relaxation and conservative dynamics.
