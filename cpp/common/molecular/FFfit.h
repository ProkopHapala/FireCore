#pragma once
/// @file FFfit.h
/// @brief Force-field parameter fitting from reference Hessian data.
///
/// ============================================================================
/// GENERAL THEORY OF HESSIAN-BASED FORCE-FIELD FITTING
/// ============================================================================
///
/// Given a force field with energy:
///   E(r) = Σ_t k_t * f_t(q_t(r))
/// where:
///   t  = term index (bond, angle, ...)
///   k_t = stiffness parameter for term t
///   f_t = energy function for term t
///   q_t = internal coordinate (bond length r, angle θ, or cos θ)
///   r  = 3N Cartesian coordinates
///
/// The Hessian (second derivative of energy w.r.t. positions) is:
///   H_αβ = ∂²E/∂r_α ∂r_β = Σ_t k_t * ∂²f_t/∂r_α ∂r_β
///
/// At fixed geometry, H is LINEAR in stiffness parameters k_t:
///   H = Σ_t k_t * A_t,  where A_t = ∂²f_t/∂r_α ∂r_β  (sensitivity matrix)
///
/// Fitting is linear least-squares:
///   min_k ||Σ_t k_t A_t - H_ref||²_W
/// Normal equations: G k = y, where G_pq = <A_p, A_q>_W, y_p = <A_p, H_ref>_W
///
/// ============================================================================
/// CHAIN RULE: SENSITIVITY MATRIX FOR f(q(r))
/// ============================================================================
///
/// For a term E = k * f(q) where q = q(r) is an internal coordinate:
///
/// First derivative (force):
///   ∂f/∂r_α = f'(q) * ∂q/∂r_α = f'(q) * b_α
/// where b_α = ∂q/∂r_α is the WILSON VECTOR.
///
/// Second derivative (Hessian / sensitivity):
///   A_αβ = ∂²f/∂r_α ∂r_β = f''(q) * b_α b_β + f'(q) * C_αβ
///
/// where:
///   b_α b_β  = (∂q/∂r_α)(∂q/∂r_β)  — rank-one outer product (Wilson × Wilson)
///   C_αβ    = ∂²q/∂r_α ∂r_β         — COORDINATE HESSIAN of the internal coordinate
///   f''(q)  = second derivative of energy function w.r.t. internal coordinate
///   f'(q)   = first derivative (PRESTRESS term; vanishes at equilibrium)
///
/// Key insight: A = f'' * b⊗b^T + f' * C
///   - f'' * b⊗b^T is the RANK-ONE part (dominant at equilibrium)
///   - f' * C is the PRESTRESS part (nonzero when geometry ≠ equilibrium)
///
/// ============================================================================
/// ENERGY FORMS AND THEIR SENSITIVITIES
/// ============================================================================
///
/// Bond (harmonic):  E = (1/2) k (r - r0)²
///   f(r) = (1/2)(r - r0)²,  f' = (r - r0) = dl,  f'' = 1
///   A_bond = b⊗b^T + dl * C_r
///   At equilibrium (dl=0): A = b⊗b^T = e⊗e^T  (rank-one only)
///
/// Angle (UFF-normalized cosine form):  E = 0.5*kθ*(cos θ - cos θ0)²/(1-cos²θ0)
///   Let c = cos θ, c0 = cos θ0, s0² = 1-c0², g(c) = 0.5*(c-c0)²/s0²
///   g'(c) = (c-c0)/s0²,  g''(c) = 1/s0²
///   ∂c/∂r_α is computed directly from normalized bond vectors, avoiding ∂θ/∂r and 1/sinθ.
///   A_angle = ( ∂c/∂r ⊗ ∂c/∂r + (c-c0)*C_cos ) / s0²
///   At equilibrium (c=c0): A = (∂c/∂r⊗∂c/∂r)/s0² = bθ⊗bθ, so kθ is the actual local angular curvature in eV/rad².
///
///   This normalization fixes the physical units: the fitted angle parameter is kθ = ∂²E/∂θ² at θ0.
///   The older unnormalized form E=K(c-c0)² has local curvature 2*K*s0², so K is not directly eV/rad².
///
/// Angle (UFF Fourier):  E = K Σ Cn cos(nθ)
///   f is a function of θ directly (not cos θ), so NO sin²θ factor:
///   A_UFF = f''(θ) * b⊗b^T + f'(θ) * C_θ
///   At equilibrium: f'(θ0)=0, f''(θ0)=1 → A = b⊗b^T
///
/// ============================================================================
/// MULTI-SYSTEM FITTING STRATEGY
/// ============================================================================
///
/// For fitting across multiple systems (e.g. different nanocrystal sizes):
///
/// 1. EQUILIBRIUM PARAMETERS (l0, c0): Chosen before stiffness fitting.
///    - Single-system local Hessian model: l0 = r_eq and c0 = cos θ_eq, so each term has zero prestress.
///    - Multi-system transferable approximation: l0 = <r> and c0 = <cos θ> over all terms of the same type.
///    Averaging c0 rather than θ is important because the angle energy is written in c = cos θ.
///
/// 2. STIFFNESS PARAMETERS (k): Fitted by linear least-squares on the
///    combined normal equations: G_total = Σ_sys G_sys, y_total = Σ_sys y_sys
///
/// This two-stage approach is a practical approximation, not a proof of a globally transferable equilibrium surface:
///   - A relaxed Cartesian geometry only implies Σ_t f'_t(q_t)∂q_t/∂r = 0, not f'_t=0 for every term.
///   - Shared l0/c0 can leave real prestress, so the f'*C coordinate-Hessian term must be retained.
///   - A fully rigorous transferable energy fit should include force data or fit linear combinations such as k and k*c0.
///
/// ============================================================================
/// PARAMETER SHARING (SYMMETRY CONSTRAINTS)
/// ============================================================================
///
/// Force-field parameters describe chemical bond TYPES, not individual bonds.
/// All Si-Si bonds share one stiffness k, all H-Si-H angles share one K.
/// This is the PRINCIPLE OF TRANSFERABILITY — the same k applies wherever
/// that bond type appears, in any molecule or system.
///
/// Without sharing: 152 Si-H bonds → 152 parameters (underdetermined)
/// With sharing:    152 Si-H bonds → 1 parameter (well-constrained)
///
/// Multi-system: the same parameter appears in ALL systems, so constraints
/// accumulate: G_total = Σ_sys G_sys, y_total = Σ_sys y_sys.
/// See ParamMap below for implementation.

#include "Vec3.h"
#include "Mat3.h"
#include <vector>
#include <string>
#include <map>
#include <tuple>
#include <algorithm>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <stdexcept>

// ============================================================
// Bond and Angle structures
// ============================================================

struct BondDef { int i, j; double l0; };
struct AngleDef { int i, j, k; double theta0; }; // j is central atom

// ============================================================
// Bond: E = (1/2) k (r - r0)²
// ============================================================
// Internal coordinate: q = r = |r_j - r_i| (bond length)
// Energy: f(r) = (1/2)(r - r0)²,  f'(r) = (r - r0) = dl,  f''(r) = 1
// Sensitivity: A = f'' * b⊗b^T + f' * C = b⊗b^T + dl * C
// At equilibrium (dl=0): A = e⊗e^T  (rank-one, no prestress)
// ============================================================

/// Wilson vector of bond length: b_α = ∂r/∂r_α
/// For bond i-j: r = |r_j - r_i|, e = (r_j - r_i)/r
///   b_i = ∂r/∂r_i = -e    (moving atom i away from j shortens bond)
///   b_j = ∂r/∂r_j = +e    (moving atom j away from i lengthens bond)
inline void bond_wilson(const Vec3d* apos, const BondDef& b, Vec3d& bi, Vec3d& bj, double& r, Vec3d& e) {
    Vec3d d; d.set_sub(apos[b.j], apos[b.i]);
    r = d.normalize();
    e = d;
    bi.set_mul(e, -1.0);
    bj = e;
}

/// Coordinate Hessian of bond length: C_αβ = ∂²r/(∂r_α ∂r_β)
///
/// Derivation: r = |r_j - r_i|, ∂r/∂r_i = -e, ∂r/∂r_j = +e
///   ∂e/∂r_i = (e⊗e^T - I)/r = -P/r
///   ∂e/∂r_j = (I - e⊗e^T)/r = +P/r
/// Therefore:
///   ∂²r/∂r_i² = ∂(-e)/∂r_i = -(-P/r) = +P/r
///   ∂²r/∂r_j² = ∂(+e)/∂r_j = +P/r
///   ∂²r/∂r_i ∂r_j = ∂(-e)/∂r_j = -P/r
///   ∂²r/∂r_j ∂r_i = ∂(+e)/∂r_i = -P/r
///
/// where P = I - e⊗e^T is the projector perpendicular to bond direction.
/// The sensitivity A_bond = b⊗b^T + dl*C uses these C blocks directly.
inline void bond_coord_hessian(const Vec3d& e, double r, Mat3d& Cii, Mat3d& Cjj, Mat3d& Cij) {
    // P = I - e*e^T
    Mat3d P; P.setOne();
    Mat3d eeT; eeT.set_outer(e, e);
    P.sub(eeT);
    double inv_r = 1.0 / r;
    Cii.set(P); Cii.mul(inv_r);
    Cjj.set(P); Cjj.mul(inv_r);
    Cij.set(P); Cij.mul(-inv_r);
}

/// Energy and force of harmonic bond (f is negative gradient, in units force)
inline double bond_energy_force(const Vec3d* apos, const BondDef& b, double k, Vec3d* f) {
    Vec3d bi, bj, e; double r;
    bond_wilson(apos, b, bi, bj, r, e);
    double dl = r - b.l0;
    double E = 0.5 * k * dl * dl;
    Vec3d fb; fb.set_mul(e, k * dl);
    f[b.i].add(fb);
    f[b.j].sub(fb);
    return E;
}

/// Sensitivity matrix dH/dk for a bond term.
///
/// A_bond = f''(r) * b⊗b^T + f'(r) * C_r
///        = 1 * b⊗b^T + dl * C_r        (since f''=1, f'=dl)
///
/// 3x3 blocks (only 4 unique due to symmetry):
///   A_ii = e⊗e^T + dl * C_ii     A_ij = -e⊗e^T + dl * C_ij
///   A_ji = A_ij^T                 A_jj = e⊗e^T + dl * C_jj
///
/// where dl = r - r0, e = unit bond vector, C from bond_coord_hessian.
/// At equilibrium (dl=0): A = e⊗e^T (rank-one, no prestress).
/// H is (3N × 3N), but we only fill the 4 blocks touching atoms i,j.
inline void bond_dHdk(const Vec3d* apos, const BondDef& b, int natoms, double* dHdk_flat) {
    Vec3d bi, bj, e; double r;
    bond_wilson(apos, b, bi, bj, r, e);
    double dl = r - b.l0;
    // F'' = 1, F' = dl
    // b⊗b^T blocks:
    // (i,i): bi⊗bi = e⊗e
    // (i,j): bi⊗bj = -e⊗e
    // (j,i): bj⊗bi = -e⊗e
    // (j,j): bj⊗bj = e⊗e
    Mat3d eeT; eeT.set_outer(e, e);
    // Prestress: dl * C
    Mat3d Cii, Cjj, Cij;
    bond_coord_hessian(e, r, Cii, Cjj, Cij);
    // Add to dHdk_flat (3N x 3N, row-major)
    Mat3d block_ii, block_ij, block_ji, block_jj;
    block_ii.set(eeT);     block_ii.add_mul(Cii, dl);
    block_ij.set(eeT);     block_ij.mul(-1.0); block_ij.add_mul(Cij, dl);
    block_ji.set(block_ij); block_ji.makeT(); // C_ij^T = C_ji
    block_jj.set(eeT);     block_jj.add_mul(Cjj, dl);
    // Write into flat array
    auto write_block = [&](int a, int c, const Mat3d& blk) {
        for (int p = 0; p < 3; p++)
            for (int q = 0; q < 3; q++)
                dHdk_flat[(a * 3 + p) * (natoms * 3) + (c * 3 + q)] += blk.array[p * 3 + q];
    };
    write_block(b.i, b.i, block_ii);
    write_block(b.i, b.j, block_ij);
    write_block(b.j, b.i, block_ji);
    write_block(b.j, b.j, block_jj);
}

// ============================================================
// Angle (UFF-normalized cosine form): E = 0.5*kθ*(cos θ - cos θ0)²/(1-cos²θ0)
// ============================================================
// Internal coordinate: q = c = cos θ = u · v, where u = (r_i-r_j)/|r_i-r_j| and v = (r_k-r_j)/|r_k-r_j|
// Energy per unit kθ: g(c) = 0.5*(c-c0)²/s0², s0² = 1-c0²
//   g'(c) = (c-c0)/s0²,  g''(c) = 1/s0²
// Direct coordinate gradient of c:
//   g_i = ∂c/∂r_i = (v - c*u)/|a|,  g_k = (u - c*v)/|b|,  g_j = -g_i-g_k
// Sensitivity for local angular stiffness kθ:
//   A = ( g⊗g^T + (c-c0)*C_cos ) / s0²
// At equilibrium: g = -sinθ0*bθ, so A = g⊗g/s0² = bθ⊗bθ and kθ has units eV/rad².
// ============================================================

/// Wilson vector of angle θ: b_α = ∂θ/∂r_α  (NOT ∂(cos θ)/∂r_α)
///
/// For angle i-j-k (j central):
///   u = (r_i - r_j)/|r_i - r_j|,  v = (r_k - r_j)/|r_k - r_j|
///   cos θ = u · v = ct,  sin θ = s
///
/// RELATION to ∂(cos θ)/∂r:
///   By chain rule: ∂(cos θ)/∂r_α = -sin θ * ∂θ/∂r_α = -s * b_α
///   So: b_α = -∂(cos θ)/∂r_α / sin θ
///
///   ∂(cos θ)/∂r_i = (v - ct*u) / |a|     [derivative of u·v w.r.t. r_i]
///   ∂θ/∂r_i = -1/s * ∂(cos θ)/∂r_i = (ct*u - v) / (|a| * s)
///
/// So: b_i = (ct*u - v) / (a*s),  b_k = (ct*v - u) / (c*s),  b_j = -(b_i + b_k)
///
/// Note: ∂(cos θ)/∂r_α = -sin θ * b_α  (used in the sin²θ factor)
/// This Wilson vector is used by the UFF Fourier angle form (f is a function of θ).
/// The cosine angle form uses angle_grad_cos_direct instead (f is a function of cos θ).
inline void angle_wilson_cos(const Vec3d* apos, const AngleDef& ang,
                             Vec3d& bi, Vec3d& bj, Vec3d& bk,
                             double& cos_theta, double& sin_theta) {
    Vec3d a; a.set_sub(apos[ang.i], apos[ang.j]);
    Vec3d c; c.set_sub(apos[ang.k], apos[ang.j]);
    double al = a.normalize();
    double cl = c.normalize();
    cos_theta = a.dot(c);
    sin_theta = sqrt(std::max(0.0, 1.0 - cos_theta * cos_theta));
    double inv_s = 1.0 / (sin_theta + 1e-14);
    Vec3d w; w.set_lincomb(cos_theta, a, -1.0, c); w.mul(inv_s / al);
    Vec3d z; z.set_lincomb(cos_theta, c, -1.0, a); z.mul(inv_s / cl);
    bi = w;
    bk = z;
    bj.set_lincomb(-1.0, bi, -1.0, bk);
}

/// Coordinate Hessian of cos θ: C_cos[α][β] = ∂²(cos θ)/(∂r_α ∂r_β)
/// Fills 9 Mat3d blocks (3×3 grid for atoms i,j,k).
///
/// Derivation (index notation, p,q = Cartesian components):
///
/// cos θ = u·v,  u = a/|a|,  v = c/|c|,  a = r_i - r_j,  c = r_k - r_j
///
/// First derivatives:
///   ∂(cos θ)/∂r_i^p = (v_p - ct*u_p) / |a|     [since ∂u/∂r_i = (I-u⊗u)/|a|]
///   ∂(cos θ)/∂r_k^p = (u_p - ct*v_p) / |c|     [since ∂v/∂r_k = (I-v⊗v)/|c|]
///   ∂(cos θ)/∂r_j   = -(∂/∂r_i + ∂/∂r_k)       [translational invariance]
///
/// Second derivative C_ii (p,q components):
///   ∂²(cos θ)/∂r_i^p ∂r_i^q = (-v_q*u_p - v_p*u_q + 3*ct*u_p*u_q - ct*δ_pq) / |a|²
///
///   Derivation: differentiate ∂(cos θ)/∂r_i^p = (v_p - ct*u_p)/|a| w.r.t. r_i^q:
///     ∂(v_p)/∂r_i^q = 0  (v independent of r_i)
///     ∂(ct)/∂r_i^q = (v_q - ct*u_q)/|a|  (first derivative)
///     ∂(u_p)/∂r_i^q = (δ_pq - u_p*u_q)/|a|  (derivative of unit vector)
///     ∂|a|/∂r_i^q = u_q
///     Product rule → numerator: -(v_q-ct*u_q)*u_p - ct*(δ_pq-u_p*u_q) - (v_p-ct*u_p)*u_q
///     Simplify: -v_q*u_p + ct*u_q*u_p - ct*δ_pq + ct*u_p*u_q - v_p*u_q + ct*u_p*u_q
///             = -v_q*u_p - v_p*u_q + 3*ct*u_p*u_q - ct*δ_pq
///
/// In matrix notation:
///   C_ii = (-v⊗u - u⊗v + 3*ct*u⊗u - ct*I) / |a|²
///   C_kk = (-u⊗v - v⊗u + 3*ct*v⊗v - ct*I) / |c|²   [by symmetry i↔k, u↔v]
///   C_ik = (I - v⊗v - u⊗u + ct*u⊗v) / (|a|*|c|)   [mixed derivative]
///   C_ki = C_ik^T                                  [symmetry of mixed partials]
///
/// Translational invariance (∂/∂r_j = -∂/∂r_i - ∂/∂r_k):
///   C_ij = -C_ii - C_ik
///   C_jk = -C_ik - C_kk         [BUG FIX: use C_ik, not C_ki]
///   C_jj = C_ii + C_ik + C_ki + C_kk
inline void angle_coord_hessian_cos(const Vec3d* apos, const AngleDef& ang, Mat3d C[3][3]) {
    Vec3d a; a.set_sub(apos[ang.i], apos[ang.j]);
    Vec3d c; c.set_sub(apos[ang.k], apos[ang.j]);
    double al = a.normalize();  // a is now u
    double cl = c.normalize();  // c is now v
    const Vec3d& u = a;
    const Vec3d& v = c;
    double ct = u.dot(c);
    double al2 = al * al, cl2 = cl * cl, alcl = al * cl;
    Mat3d I; I.setOne();
    Mat3d uuT, vvT, vuT, uvT;
    uuT.set_outer(u, u); vvT.set_outer(v, v);
    vuT.set_outer(v, u); uvT.set_outer(u, v);
    // C_ii = (-v⊗u - u⊗v + 3*ct*u⊗u - ct*I) / |a|²
    C[0][0].set(I); C[0][0].mul(-ct); C[0][0].add_mul(uuT, 3.0*ct); C[0][0].sub(vuT); C[0][0].sub(uvT); C[0][0].mul(1.0/al2);
    // C_kk = (-u⊗v - v⊗u + 3*ct*v⊗v - ct*I) / |c|²
    C[2][2].set(I); C[2][2].mul(-ct); C[2][2].add_mul(vvT, 3.0*ct); C[2][2].sub(uvT); C[2][2].sub(vuT); C[2][2].mul(1.0/cl2);
    // C_ik = (I - v⊗v - u⊗u + ct*u⊗v) / (|a|*|c|)
    C[0][2].set(I); C[0][2].sub(vvT); C[0][2].sub(uuT); C[0][2].add_mul(uvT, ct); C[0][2].mul(1.0/alcl);
    // C_ki = C_ik^T
    C[2][0].setT(C[0][2]);
    // C_ij = -C_ii - C_ik
    C[0][1].set(C[0][0]); C[0][1].add(C[0][2]); C[0][1].mul(-1.0);
    // C_ji = C_ij^T
    C[1][0].setT(C[0][1]);
    // C_jk = -C_ik - C_kk   (use C_ik, not C_ki; C_ik ≠ C_ki in general)
    C[1][2].set(C[0][2]); C[1][2].add(C[2][2]); C[1][2].mul(-1.0);
    // C_kj = C_jk^T
    C[2][1].setT(C[1][2]);
    // C_jj = C_ii + C_ik + C_ki + C_kk
    C[1][1].set(C[0][0]); C[1][1].add(C[0][2]); C[1][1].add(C[2][0]); C[1][1].add(C[2][2]);
}

/// Direct gradient of c = cos θ with respect to Cartesian positions.
///
/// WHY cos θ INSTEAD OF θ?
/// The angle energy E = ½ kθ (cos θ - cos θ0)² / (1-cos²θ0) is naturally a
/// function of c = cos θ = u·v, NOT of θ = arccos(u·v).  Using c directly:
///   - Avoids the 1/sin θ singularity in ∂θ/∂r (which diverges at θ → 0, π)
///   - ∂c/∂r is always finite and well-conditioned (it's just a dot product derivative)
///   - The sensitivity A = (g⊗g + (c-c0)·C_cos) / s0² is smooth everywhere
///
/// DERIVATION of ∂(cos θ)/∂r:
///   c = u·v,  u = a/|a|,  v = c/|c|,  a = r_i - r_j,  c = r_k - r_j
///   ∂u/∂r_i = (I - u⊗u)/|a|  (derivative of unit vector: transverse projector / length)
///   ∂v/∂r_i = 0  (v depends only on r_j, r_k)
///   ∂c/∂r_i = ∂(u·v)/∂r_i = (∂u/∂r_i)·v = (v - (u·v)u)/|a| = (v - c*u)/|a|
///   By symmetry (i↔k, u↔v): ∂c/∂r_k = (u - c*v)/|c|
///   Translational invariance: ∂c/∂r_j = -(∂c/∂r_i + ∂c/∂r_k)
inline void angle_grad_cos_direct(const Vec3d* apos, const AngleDef& ang, Vec3d& gi, Vec3d& gj, Vec3d& gk, double& c) {
    Vec3d u; u.set_sub(apos[ang.i], apos[ang.j]);
    Vec3d v; v.set_sub(apos[ang.k], apos[ang.j]);
    double ul = u.normalize();
    double vl = v.normalize();
    c = u.dot(v);
    gi.set_lincomb(1.0, v, -c, u); gi.mul(1.0/ul); // ∂c/∂r_i = (v-c*u)/|a|
    gk.set_lincomb(1.0, u, -c, v); gk.mul(1.0/vl); // ∂c/∂r_k = (u-c*v)/|b|
    gj.set_lincomb(-1.0, gi, -1.0, gk);            // translational invariance
}

/// Sensitivity matrix dH/dkθ for normalized cosine angle: E = 0.5*kθ*(c-c0)²/(1-c0²)
///
/// A = ( g⊗g^T + (c-c0)*C_cos ) / (1-c0²), where g = ∂c/∂r.
/// At equilibrium: A = bθ⊗bθ, so kθ is directly the local angular stiffness in eV/rad².
inline void angle_dHdk_cos(const Vec3d* apos, const AngleDef& ang, int natoms, double* dHdk_flat) {
    Vec3d gi, gj, gk; double c;
    angle_grad_cos_direct(apos, ang, gi, gj, gk, c);
    double c0 = cos(ang.theta0);
    double s02 = 1.0 - c0*c0;
    if (s02 < 1e-14) throw std::runtime_error("angle_dHdk_cos: sin(theta0)^2 too small for normalized cosine angle");
    double dc = c - c0;
    double scale = 1.0 / s02;
    auto add_outer = [&](int a, int bidx, const Vec3d& va, const Vec3d& vb, double sc) {
        for (int p = 0; p < 3; p++)
            for (int q = 0; q < 3; q++)
                dHdk_flat[(a * 3 + p) * (natoms * 3) + (bidx * 3 + q)] += sc * va.array[p] * vb.array[q];
    };
    Vec3d g[3] = {gi, gj, gk};
    int atoms[3] = {ang.i, ang.j, ang.k};
    for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++)
            add_outer(atoms[a], atoms[b], g[a], g[b], scale);
    if (fabs(dc) > 1e-12) {
        Mat3d Ccos[3][3];
        angle_coord_hessian_cos(apos, ang, Ccos);
        for (int a = 0; a < 3; a++)
            for (int b = 0; b < 3; b++)
                for (int p = 0; p < 3; p++)
                    for (int q = 0; q < 3; q++)
                        dHdk_flat[(atoms[a] * 3 + p) * (natoms * 3) + (atoms[b] * 3 + q)] += scale * dc * Ccos[a][b].array[p * 3 + q];
    }
}

// ============================================================
// Angle (UFF Fourier form): E = K * (C0 + C1*cos(θ) + C2*cos(2θ))
// ============================================================
// Unlike the cosine form, UFF energy is a function of θ DIRECTLY (not cos θ).
// Therefore NO sin²θ factor in the sensitivity:
//   A_UFF = f''(θ) * b⊗b^T + f'(θ) * C_θ
// where C_θ = ∂²θ/∂r² (coordinate Hessian of θ, NOT of cos θ).
// At equilibrium: f'(θ0)=0, f''(θ0)=1 → A = b⊗b^T
// ============================================================

struct UFFAngleParams { double C0, C1, C2; };

/// Compute UFF Fourier coefficients from equilibrium angle.
///
/// The UFF angle energy is a Fourier series in θ:
///   E = K * f(θ),  f(θ) = C0 + C1*cos(θ) + C2*cos(2θ)
///
/// NORMALIZATION: We require f'(θ0) = 0 (equilibrium) and f''(θ0) = 1 (K is the curvature).
///
/// Derivation:
///   f'(θ)  = -C1*sin(θ) - 2*C2*sin(2θ) = -C1*s - 4*C2*s*c
///   f''(θ) = -C1*cos(θ) - 4*C2*cos(2θ) = -C1*c - 4*C2*(2c²-1)
///
/// At θ = θ0 (c = c0, s = s0):
///   f'(θ0) = -s0*(C1 + 4*C2*c0) = 0  →  C1 = -4*C2*c0
///   f''(θ0) = -C1*c0 - 4*C2*(2*c0²-1) = 4*C2*c0² - 4*C2*(2*c0²-1) = 4*C2*(1-c0²) = 4*C2*s0² = 1
///   →  C2 = 1/(4*s0²)
///   C0 = -C1*c0 - C2*(2*c0²-1) = 4*C2*c0² - C2*(2*c0²-1) = C2*(2*c0²+1)
///   (C0 ensures f(θ0) = 0 so that E = 0 at equilibrium)
inline UFFAngleParams uff_angle_coeffs(double theta0) {
    double c0 = cos(theta0), s0 = sin(theta0);
    double C2 = 1.0 / (4.0 * s0 * s0);
    double C1 = -4.0 * C2 * c0;
    double C0 = C2 * (2.0 * c0 * c0 + 1.0);
    return {C0, C1, C2}; // normalized UFF Fourier coefficients
}

/// Sensitivity dH/dK for UFF Fourier angle.
///
/// f(θ) = C0 + C1*cos(θ) + C2*cos(2θ)  (function of θ, NOT cos θ)
/// f'(θ)  = -C1*sin(θ) - 2*C2*sin(2θ)   = -C1*s - 4*C2*s*c
/// f''(θ) = -C1*cos(θ) - 4*C2*cos(2θ)   = -C1*c - 4*C2*(2c²-1)
///
/// A = f''(θ) * b⊗b^T + f'(θ) * C_θ
///
/// NO sin²θ factor (unlike cosine form) because f is directly a function of θ.
/// At equilibrium: f'(θ0)=0, f''(θ0)=1 → A = b⊗b^T
///
/// NOTE: Prestress term f'*C_θ is currently skipped (TODO).
///
/// RELATION between C_θ and C_cos:
///   c = cos θ, so ∂c/∂r = -sin θ * ∂θ/∂r = -s * b  (chain rule)
///   Differentiating again:
///     C_cos = ∂²c/∂r² = ∂(-s*b)/∂r = -c*(b⊗b^T)*s - s*C_θ
///   (using ∂s/∂r = ∂(sin θ)/∂r = cos θ * ∂θ/∂r = c*b, and ∂b/∂r = C_θ)
///   Solving for C_θ:
///     C_θ = -(C_cos + c * b⊗b^T) / s
///   This diverges at θ → 0, π (s → 0), which is why the cosine form is preferred
///   for the actual sensitivity computation.
inline void angle_dHdk_uff(const Vec3d* apos, const AngleDef& ang, const UFFAngleParams& par,
                            int natoms, double* dHdk_flat) {
    Vec3d bi, bj, bk; double c, s;
    angle_wilson_cos(apos, ang, bi, bj, bk, c, s);
    // F''(θ) = -C1*cos(θ) - 4*C2*cos(2θ)
    double Fpp = -par.C1 * c - 4.0 * par.C2 * (2.0 * c * c - 1.0);
    // F'(θ) = -C1*sin(θ) - 2*C2*sin(2θ)  (for prestress)
    double Fp = -par.C1 * s - 2.0 * par.C2 * 2.0 * s * c;
    // Rank-one part: F'' * b⊗b^T
    auto add_outer = [&](int a, int bidx, const Vec3d& va, const Vec3d& vb, double scale) {
        for (int p = 0; p < 3; p++)
            for (int q = 0; q < 3; q++)
                dHdk_flat[(a * 3 + p) * (natoms * 3) + (bidx * 3 + q)] += scale * va.array[p] * vb.array[q];
    };
    add_outer(ang.i, ang.i, bi, bi, Fpp);
    add_outer(ang.i, ang.j, bi, bj, Fpp);
    add_outer(ang.i, ang.k, bi, bk, Fpp);
    add_outer(ang.j, ang.i, bj, bi, Fpp);
    add_outer(ang.j, ang.j, bj, bj, Fpp);
    add_outer(ang.j, ang.k, bj, bk, Fpp);
    add_outer(ang.k, ang.i, bk, bi, Fpp);
    add_outer(ang.k, ang.j, bk, bj, Fpp);
    add_outer(ang.k, ang.k, bk, bk, Fpp);
    // Prestress: F' * C_θ — skipped (small at relaxed geometry)
    // TODO: add for full generality
}

// ============================================================
// Per-term sensitivity: compute only the 3x3 blocks for one bond/angle
// (for gradient descent — avoids building full 3N×3N matrices)
// ============================================================

/// Per-term bond sensitivity blocks: fills 4 Mat3d blocks (ii, ij, ji, jj)
inline void bond_dHdk_blocks(const Vec3d* apos, const BondDef& b,
                             Mat3d& Bii, Mat3d& Bij, Mat3d& Bji, Mat3d& Bjj) {
    Vec3d bi, bj, e; double r;
    bond_wilson(apos, b, bi, bj, r, e);
    double dl = r - b.l0;
    Mat3d eeT; eeT.set_outer(e, e);
    Mat3d Cii, Cjj, Cij;
    bond_coord_hessian(e, r, Cii, Cjj, Cij);
    Bii.set(eeT);     Bii.add_mul(Cii, dl);
    Bij.set(eeT);     Bij.mul(-1.0); Bij.add_mul(Cij, dl);
    Bji.set(Bij);     Bji.makeT();
    Bjj.set(eeT);     Bjj.add_mul(Cjj, dl);
}

/// Per-term angle sensitivity blocks (same as angle_dHdk_cos but block form).
/// A = ( g⊗g^T + (c-c0)*C_cos ) / (1-c0²), where g = ∂cosθ/∂r.
inline void angle_dHdk_cos_blocks(const Vec3d* apos, const AngleDef& ang,
                                  Mat3d B[3][3]) {
    Vec3d gi, gj, gk; double c;
    angle_grad_cos_direct(apos, ang, gi, gj, gk, c);
    double c0 = cos(ang.theta0);
    double s02 = 1.0 - c0*c0;
    if (s02 < 1e-14) throw std::runtime_error("angle_dHdk_cos_blocks: sin(theta0)^2 too small for normalized cosine angle");
    double dc = c - c0;
    double scale = 1.0 / s02;
    Vec3d g[3] = {gi, gj, gk};
    for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++) {
            B[a][b].set_outer(g[a], g[b]);
            B[a][b].mul(scale);
        }
    if (fabs(dc) > 1e-12) {
        Mat3d Ccos[3][3];
        angle_coord_hessian_cos(apos, ang, Ccos);
        for (int a = 0; a < 3; a++)
            for (int b = 0; b < 3; b++)
                B[a][b].add_mul(Ccos[a][b], scale*dc);
    }
}

/// Inner product <A_term, ΔH>_W for a bond term — only touches 4 blocks
inline double bond_term_inner_product(const Vec3d* apos, const BondDef& b,
                                       const double* dH, const double* weight, int natoms) {
    Mat3d Bii, Bij, Bji, Bjj;
    bond_dHdk_blocks(apos, b, Bii, Bij, Bji, Bjj);
    int n3 = natoms * 3;
    auto dot_blk = [&](int a, int c, const Mat3d& blk) -> double {
        double s = 0.0;
        for (int p = 0; p < 3; p++)
            for (int q = 0; q < 3; q++) {
                int idx = (a * 3 + p) * n3 + (c * 3 + q);
                double w = weight ? weight[idx] * weight[idx] : 1.0;
                s += w * blk.array[p * 3 + q] * dH[idx];
            }
        return s;
    };
    return dot_blk(b.i, b.i, Bii) + dot_blk(b.i, b.j, Bij) +
           dot_blk(b.j, b.i, Bji) + dot_blk(b.j, b.j, Bjj);
}

/// Inner product <A_term, ΔH>_W for an angle term — only touches 9 blocks
inline double angle_term_inner_product(const Vec3d* apos, const AngleDef& ang,
                                        const double* dH, const double* weight, int natoms) {
    Mat3d B[3][3];
    angle_dHdk_cos_blocks(apos, ang, B);
    int n3 = natoms * 3;
    int atoms[3] = {ang.i, ang.j, ang.k};
    double s = 0.0;
    for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++) {
            for (int p = 0; p < 3; p++)
                for (int q = 0; q < 3; q++) {
                    int idx = (atoms[a] * 3 + p) * n3 + (atoms[b] * 3 + q);
                    double w = weight ? weight[idx] * weight[idx] : 1.0;
                    s += w * B[a][b].array[p * 3 + q] * dH[idx];
                }
        }
    return s;
}

/// Accumulate term contribution into model Hessian (for gradient descent forward pass)
inline void bond_accumulate_H(Mat3d Bii, Mat3d Bij, Mat3d Bji, Mat3d Bjj,
                               int i, int j, double k_term, double* H, int natoms) {
    int n3 = natoms * 3;
    auto add_blk = [&](int a, int c, const Mat3d& blk) {
        for (int p = 0; p < 3; p++)
            for (int q = 0; q < 3; q++)
                H[(a * 3 + p) * n3 + (c * 3 + q)] += k_term * blk.array[p * 3 + q];
    };
    add_blk(i, i, Bii); add_blk(i, j, Bij);
    add_blk(j, i, Bji); add_blk(j, j, Bjj);
}

inline void angle_accumulate_H(Mat3d B[3][3], int atoms[3], double k_term, double* H, int natoms) {
    int n3 = natoms * 3;
    for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++) {
            for (int p = 0; p < 3; p++)
                for (int q = 0; q < 3; q++)
                    H[(atoms[a] * 3 + p) * n3 + (atoms[b] * 3 + q)] += k_term * B[a][b].array[p * 3 + q];
        }
}

// ============================================================
// ParamMap: symmetry-based parameter sharing
// ============================================================
//
// PHYSICAL MOTIVATION:
// --------------------
// In a typical nanocrystal or molecule, many bond/angle terms are
// CHEMICALLY EQUIVALENT — e.g. all Si-Si bonds in a diamond lattice
// have the same stiffness, all H-Si-H angles in SiH4 are identical.
//
// Rather than fitting one parameter per individual bond/angle (which
// would be underdetermined — more parameters than Hessian data), we
// exploit this equivalence:
//
//   - Bonds with the same element pair (e.g. Si-Si) share ONE stiffness k
//   - Angles with the same element triple (e.g. Si-Si-H) share ONE K
//
// This is the PRINCIPLE OF TRANSFERABILITY: force-field parameters
// describe chemical bond TYPES, not individual bonds.  The same k
// applies wherever that bond type appears — in different molecules,
// different positions, different systems.
//
// MULTI-SYSTEM IMPLICATION:
//   When fitting across multiple systems (e.g. SiH4 + Si nanocrystals),
//   the same parameter appears in ALL systems.  This dramatically
//   increases the data-to-parameter ratio:
//     - SiH4 alone: 4 Si-H bonds → 1 k_SiH parameter  (4 constraints)
//     - 6 systems:  152 Si-H bonds → 1 k_SiH parameter (152 constraints)
//   The normal equations accumulate: G_total = Σ_sys G_sys, y_total = Σ_sys y_sys
//
// IMPLEMENTATION:
//   ParamMap is a sparse index array: bond_param_idx[ib] → free param index.
//   Terms with the same index contribute to the same row/column of the
//   normal equation matrix G.  The mapping can be based on:
//     - Element types (current: Si-Si, Si-H, etc.)
//     - Chemical environment (future: Si-Si with 4 neighbors vs 3)
//     - Manual assignment (for custom symmetry groups)
//
// ============================================================

/// @brief Maps each bond/angle term to a free parameter index.
/// Implements symmetry constraints: chemically equivalent terms share
/// one parameter, reducing the number of fitted unknowns.
/// The mapping is a simple index array (sparse, one non-zero per row).
struct ParamMap {
    int n_free_params = 0;
    std::vector<int> bond_param_idx;   // [nbonds] → free param index
    std::vector<int> angle_param_idx;  // [nangles] → free param index

    void resize(int nbonds, int nangles) {
        bond_param_idx.resize(nbonds, -1);
        angle_param_idx.resize(nangles, -1);
    }

    /// Assign all bonds with same element pair to same parameter
    void assign_bond_types(const std::vector<BondDef>& bonds, const std::vector<std::string>& symbols) {
        std::map<std::pair<std::string,std::string>, int> type_map;
        for (size_t ib = 0; ib < bonds.size(); ib++) {
            auto key = std::minmax(symbols[bonds[ib].i], symbols[bonds[ib].j]);
            if (type_map.find(key) == type_map.end())
                type_map[key] = n_free_params++;
            bond_param_idx[ib] = type_map[key];
        }
    }

    /// Assign all angles with same element triple to same parameter.
    /// Key is symmetric under i↔k: outer atoms are sorted, central atom stays fixed.
    void assign_angle_types(const std::vector<AngleDef>& angles, const std::vector<std::string>& symbols) {
        int base = n_free_params;
        std::map<std::tuple<std::string,std::string,std::string>, int> type_map;
        for (size_t ia = 0; ia < angles.size(); ia++) {
            const std::string& si = symbols[angles[ia].i];
            const std::string& sk = symbols[angles[ia].k];
            const std::string& sj = symbols[angles[ia].j];
            auto key = (si < sk)
                ? std::make_tuple(si, sj, sk)
                : std::make_tuple(sk, sj, si);
            if (type_map.find(key) == type_map.end())
                type_map[key] = n_free_params++;
            angle_param_idx[ia] = type_map[key];
        }
    }
};

// ============================================================
// FFfit class: linear least-squares + gradient descent fitting
// ============================================================

class FFfit {
public:
    int natoms = 0;
    Vec3d* apos = nullptr;
    std::vector<BondDef> bonds;
    std::vector<AngleDef> angles;
    std::vector<std::string> symbols; // element symbols for type assignment

    // Parameter mapping (symmetry constraints)
    ParamMap param_map;

    // Results
    std::vector<double> k_free;  // fitted free parameters

    /// Set geometry from positions array.
    /// Deletes any existing geometry (safe to call multiple times).
    void set_geometry(int n, const Vec3d* positions) {
        if (apos) { delete[] apos; apos = nullptr; }
        natoms = n;
        apos = new Vec3d[n];
        for (int i = 0; i < n; i++) apos[i] = positions[i];
    }

    /// Explicit default constructor (required because copy/assignment are deleted).
    FFfit() = default;

    /// Disable copy and default assignment to avoid double-delete of apos.
    FFfit(const FFfit&) = delete;
    FFfit& operator=(const FFfit&) = delete;

    /// Move operations are also disabled to keep lifecycle simple.
    FFfit(FFfit&&) = delete;
    FFfit& operator=(FFfit&&) = delete;

    ~FFfit() { delete[] apos; }

    /// Add a bond (parameter index assigned later via param_map)
    void add_bond(int i, int j, double l0) {
        bonds.push_back({i, j, l0});
    }

    /// Add an angle (parameter index assigned later via param_map)
    void add_angle(int i, int j, int k, double theta0) {
        angles.push_back({i, j, k, theta0});
    }

    /// Set element symbols for type assignment
    void set_symbols(const std::vector<std::string>& syms) { symbols = syms; }

    /// Auto-assign parameter types based on element pairs/triples
    void auto_assign_types() {
        param_map.n_free_params = 0;
        param_map.resize(bonds.size(), angles.size());
        param_map.assign_bond_types(bonds, symbols);
        param_map.assign_angle_types(angles, symbols);
    }

    /// Manual: set bond parameter index
    void set_bond_param(int ibond, int param_idx) {
        if ((int)param_map.bond_param_idx.size() <= ibond) param_map.bond_param_idx.resize(bonds.size(), -1);
        param_map.bond_param_idx[ibond] = param_idx;
        param_map.n_free_params = std::max(param_map.n_free_params, param_idx + 1);
    }

    /// Manual: set angle parameter index
    void set_angle_param(int iangle, int param_idx) {
        if ((int)param_map.angle_param_idx.size() <= iangle) param_map.angle_param_idx.resize(angles.size(), -1);
        param_map.angle_param_idx[iangle] = param_idx;
        param_map.n_free_params = std::max(param_map.n_free_params, param_idx + 1);
    }

    // === Method 1: Linear least-squares (normal equations) ===
    //
    // DERIVATION of normal equations:
    //   Loss: L(k) = ||Σ_f k_f A_f - H_ref||²_W = Σ_α w_α² (Σ_f k_f A_f[α] - H_ref[α])²
    //   ∂L/∂k_p = 2 Σ_α w_α² A_p[α] (Σ_f k_f A_f[α] - H_ref[α])
    //            = 2 (Σ_f G_pf k_f - y_p) = 0
    //   → G k = y  (normal equations)
    //   where G_pf = <A_p, A_f>_W = Σ_α w_α² A_p[α] A_f[α]
    //         y_p  = <A_p, H_ref>_W = Σ_α w_α² A_p[α] H_ref[α]
    //
    //   The normal equations square the condition number (κ(G) = κ(A)²),
    //   so for ill-conditioned systems the direct stacked solve in Python
    //   (solve_regularized_lsq) is preferred.

    /// Build sensitivity matrices A_f = dH/dk_f for each free parameter
    /// Returns vector of (n_free_params) matrices, each (3N x 3N) flattened row-major
    std::vector<std::vector<double>> build_sensitivity_matrices() {
        int n3 = natoms * 3;
        int np = param_map.n_free_params;
        std::vector<std::vector<double>> A(np, std::vector<double>(n3 * n3, 0.0));
        for (size_t ib = 0; ib < bonds.size(); ib++) {
            int f = param_map.bond_param_idx[ib];
            if (f < 0) continue;
            bond_dHdk(apos, bonds[ib], natoms, A[f].data());
        }
        for (size_t ia = 0; ia < angles.size(); ia++) {
            int f = param_map.angle_param_idx[ia];
            if (f < 0) continue;
            angle_dHdk_cos(apos, angles[ia], natoms, A[f].data());
        }
        return A;
    }

    /// Fit via linear least-squares: min_k ||Σ k_f A_f - H_ref||^2_W
    std::vector<double> fit_hessian(const double* H_ref, const double* H_0 = nullptr,
                                     const double* weight = nullptr) {
        int n3 = natoms * 3;
        int np = param_map.n_free_params;
        auto A = build_sensitivity_matrices();
        // Normal equations: G[p][q] = <A_p, A_q>_W, y[p] = <A_p, H_ref-H_0>_W
        std::vector<std::vector<double>> G(np, std::vector<double>(np, 0.0));
        std::vector<double> y(np, 0.0);
        for (int p = 0; p < np; p++) {
            for (int q = p; q < np; q++) {
                double dot = 0.0;
                for (int i = 0; i < n3 * n3; i++) {
                    double w = weight ? weight[i] * weight[i] : 1.0;
                    dot += w * A[p][i] * A[q][i];
                }
                G[p][q] = dot; G[q][p] = dot;
            }
            for (int i = 0; i < n3 * n3; i++) {
                double w = weight ? weight[i] * weight[i] : 1.0;
                y[p] += w * A[p][i] * (H_ref[i] - (H_0 ? H_0[i] : 0.0));
            }
        }
        k_free.assign(np, 0.0);
        // Flatten G to 1D for solve_linear_system
        std::vector<double> G_flat(np * np, 0.0);
        for (int p = 0; p < np; p++)
            for (int q = 0; q < np; q++)
                G_flat[p * np + q] = G[p][q];
        solve_linear_system(G_flat.data(), y.data(), k_free.data(), np);
        return k_free;
    }

    // === Method 2: Gradient descent (variational derivatives only) ===
    //
    // DERIVATION of gradient:
    //   L(k) = ||H_model(k) - H_ref||²_W = Σ_α w_α² (H_model[α] - H_ref[α])²
    //   H_model = Σ_t k_{p(t)} A_t  (sum over terms, p(t) = param index of term t)
    //   ∂H_model[α]/∂k_p = Σ_{t→p} A_t[α]  (sum over terms mapped to param p)
    //   ∂L/∂k_p = 2 Σ_α w_α² (H_model[α] - H_ref[α]) * Σ_{t→p} A_t[α]
    //           = 2 Σ_{t→p} <A_t, ΔH>_W
    //   where ΔH = H_model - H_ref and <A, B>_W = Σ_α w_α² A[α] B[α].
    //
    //   KEY EFFICIENCY: Each term A_t only touches 4 (bond) or 9 (angle) 3×3 blocks,
    //   so <A_t, ΔH>_W can be computed without building the full 3N×3N A_t matrix.
    //   This makes gradient descent O(n_terms * 9) per iteration vs O(n_params * (3N)²)
    //   for the full sensitivity matrix approach.
    //
    //   Momentum update: v = μ*v - lr*grad; k += v
    //   This is heavy-ball gradient descent with momentum coefficient μ.
    //   Converges faster than plain GD for ill-conditioned G (different eigenvalues).

    /// Compute gradient of loss L = ||H_model - H_ref||^2_W w.r.t. free parameters
    /// Uses per-term sensitivity (no full 3N×3N matrices needed for gradient)
    /// Returns gradient[np] and optionally loss value
    std::vector<double> compute_gradient(const double* H_ref, const double* H_0,
                                          const double* weight, const std::vector<double>& k,
                                          double* loss_out = nullptr) {
        int n3 = natoms * 3;
        int np = param_map.n_free_params;
        // Forward pass: build H_model = Σ k_f * A_f
        std::vector<double> H_model(n3 * n3, 0.0);
        // Precompute per-term blocks and accumulate
        for (size_t ib = 0; ib < bonds.size(); ib++) {
            int f = param_map.bond_param_idx[ib];
            if (f < 0) continue;
            Mat3d Bii, Bij, Bji, Bjj;
            bond_dHdk_blocks(apos, bonds[ib], Bii, Bij, Bji, Bjj);
            bond_accumulate_H(Bii, Bij, Bji, Bjj, bonds[ib].i, bonds[ib].j, k[f], H_model.data(), natoms);
        }
        for (size_t ia = 0; ia < angles.size(); ia++) {
            int f = param_map.angle_param_idx[ia];
            if (f < 0) continue;
            Mat3d B[3][3];
            angle_dHdk_cos_blocks(apos, angles[ia], B);
            int atoms[3] = {angles[ia].i, angles[ia].j, angles[ia].k};
            angle_accumulate_H(B, atoms, k[f], H_model.data(), natoms);
        }
        // Residual: dH = H_model - H_ref + H_0
        std::vector<double> dH(n3 * n3);
        double loss = 0.0;
        for (int i = 0; i < n3 * n3; i++) {
            double h0 = H_0 ? H_0[i] : 0.0;
            dH[i] = H_model[i] - H_ref[i] + h0;
            double w = weight ? weight[i] * weight[i] : 1.0;
            loss += w * dH[i] * dH[i];
        }
        if (loss_out) *loss_out = loss;
        // Backward pass: gradient[f] = 2 * Σ_{term→f} <A_term, dH>_W
        std::vector<double> grad(np, 0.0);
        for (size_t ib = 0; ib < bonds.size(); ib++) {
            int f = param_map.bond_param_idx[ib];
            if (f < 0) continue;
            grad[f] += 2.0 * bond_term_inner_product(apos, bonds[ib], dH.data(), weight, natoms);
        }
        for (size_t ia = 0; ia < angles.size(); ia++) {
            int f = param_map.angle_param_idx[ia];
            if (f < 0) continue;
            grad[f] += 2.0 * angle_term_inner_product(apos, angles[ia], dH.data(), weight, natoms);
        }
        return grad;
    }

    /// Fit via gradient descent with momentum
    /// H_ref: (3N x 3N) reference Hessian
    /// H_0: constant part (can be nullptr)
    /// weight: per-element weight (can be nullptr)
    /// lr: learning rate, momentum: momentum coefficient, max_iter, tol: convergence tolerance
    std::vector<double> fit_gradient_descent(const double* H_ref, const double* H_0 = nullptr,
                                              const double* weight = nullptr,
                                              double lr = 1e-3, double momentum = 0.9,
                                              int max_iter = 1000, double tol = 1e-8,
                                              bool verbose = true) {
        int np = param_map.n_free_params;
        k_free.assign(np, 1.0); // start from uniform guess
        std::vector<double> velocity(np, 0.0);
        double prev_loss = 1e30;
        for (int iter = 0; iter < max_iter; iter++) {
            double loss = 0.0;
            auto grad = compute_gradient(H_ref, H_0, weight, k_free, &loss);
            // Check convergence
            double grad_norm = 0.0;
            for (int f = 0; f < np; f++) grad_norm += grad[f] * grad[f];
            grad_norm = sqrt(grad_norm);
            if (verbose && (iter % 100 == 0 || iter < 10))
                printf("  GD iter %4d: loss=%.6e grad_norm=%.6e\n", iter, loss, grad_norm);
            if (grad_norm < tol || (iter > 0 && fabs(prev_loss - loss) < tol * prev_loss)) {
                if (verbose) printf("  GD converged at iter %d: loss=%.6e\n", iter, loss);
                break;
            }
            prev_loss = loss;
            // Momentum update
            for (int f = 0; f < np; f++) {
                velocity[f] = momentum * velocity[f] - lr * grad[f];
                k_free[f] += velocity[f];
            }
        }
        return k_free;
    }

    /// Compute model Hessian from fitted parameters
    std::vector<double> compute_model_hessian() const {
        int n3 = natoms * 3;
        std::vector<double> H(n3 * n3, 0.0);
        for (size_t ib = 0; ib < bonds.size(); ib++) {
            int f = param_map.bond_param_idx[ib];
            if (f < 0) continue;
            Mat3d Bii, Bij, Bji, Bjj;
            bond_dHdk_blocks(apos, bonds[ib], Bii, Bij, Bji, Bjj);
            bond_accumulate_H(Bii, Bij, Bji, Bjj, bonds[ib].i, bonds[ib].j, k_free[f], H.data(), natoms);
        }
        for (size_t ia = 0; ia < angles.size(); ia++) {
            int f = param_map.angle_param_idx[ia];
            if (f < 0) continue;
            Mat3d B[3][3];
            angle_dHdk_cos_blocks(apos, angles[ia], B);
            int atoms[3] = {angles[ia].i, angles[ia].j, angles[ia].k};
            angle_accumulate_H(B, atoms, k_free[f], H.data(), natoms);
        }
        return H;
    }

    /// Print fitted parameters
    void print_params() const {
        printf("=== FFfit Parameters (%d free) ===\n", param_map.n_free_params);
        for (int f = 0; f < param_map.n_free_params; f++)
            printf("  k_free[%d] = %.6e\n", f, k_free[f]);
    }

private:
    void solve_linear_system(double* G, double* y, double* x, int n) {
        for (int i = 0; i < n; i++) {
            int p = i;
            double maxval = fabs(G[i * n + i]);
            for (int r = i + 1; r < n; r++) {
                double val = fabs(G[r * n + i]);
                if (val > maxval) { maxval = val; p = r; }
            }
            if (maxval < 1e-15) {
                char msg[256];
                snprintf(msg, sizeof(msg), "singular matrix in solve_linear_system (row %d, pivot=%.2e)", i, maxval);
                throw std::runtime_error(msg);
            }
            if (p != i) {
                for (int c = 0; c < n; c++) std::swap(G[i * n + c], G[p * n + c]);
                std::swap(y[i], y[p]);
            }
            for (int r = i + 1; r < n; r++) {
                double factor = G[r * n + i] / G[i * n + i];
                for (int c = i; c < n; c++) G[r * n + c] -= factor * G[i * n + c];
                y[r] -= factor * y[i];
            }
        }
        for (int i = n - 1; i >= 0; i--) {
            x[i] = y[i];
            for (int j = i + 1; j < n; j++) x[i] -= G[i * n + j] * x[j];
            x[i] /= G[i * n + i];
        }
    }
};
