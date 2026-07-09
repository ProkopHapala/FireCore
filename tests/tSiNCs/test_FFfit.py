#!/usr/bin/env python3
"""Fit bond/angle force-field parameters to PySCF reference Hessians.

Loads PySCF hessian .npy data from results directory, builds bond/angle topology
from relaxed .xyz geometry, and solves linear least-squares for stiffness parameters.

Usage:
    python3 tests/tSiNCs/test_FFfit.py
    python3 tests/tSiNCs/test_FFfit.py --case SiH4
    python3 tests/tSiNCs/test_FFfit.py --case Si_R5p0 --mass-weight
"""

import os, sys, json, argparse
from collections import deque
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pyBall import FFfit as FFfit_cpp

# --- Constants
BOHR_TO_ANG = 0.5291772109
HARTREE_TO_EV = 27.211386245988
BOHR_TO_ANG_INV = 1.0 / BOHR_TO_ANG
HARTREE_PER_BOHR2_TO_EV_PER_ANG2 = HARTREE_TO_EV * BOHR_TO_ANG_INV**2

RESULTS_DIR = "/home/prokop/SIMULATIONS/jobs_pyscf_vib_OUT_small_nc/results"

def angle_type_key(symbols, i, j, k):
    si, sk = symbols[i], symbols[k]
    return (si, symbols[j], sk) if si <= sk else (sk, symbols[j], si)

def load_pyscf_case(case_dir):
    """Load PySCF results for a case."""
    d = {}
    d['hessian'] = np.load(os.path.join(case_dir, 'hessian.npy'))      # (natoms,3,natoms,3) in Ha/Bohr^2
    d['modes']   = np.load(os.path.join(case_dir, 'modes.npy'))         # (nmodes, natoms, 3)
    d['freqs']   = np.load(os.path.join(case_dir, 'frequencies_cm1.npy'))
    d['masses']  = np.load(os.path.join(case_dir, 'masses.npy'))        # (natoms,)
    # Load relaxed geometry
    xyz_path = os.path.join(case_dir, 'relaxed.xyz')
    with open(xyz_path) as f:
        lines = f.readlines()
    natoms = int(lines[0].strip())
    symbols = []
    positions = np.zeros((natoms, 3))
    for i in range(natoms):
        parts = lines[2 + i].split()
        symbols.append(parts[0])
        positions[i] = [float(parts[1]), float(parts[2]), float(parts[3])]
    d['symbols'] = symbols
    d['positions'] = positions  # in Angstrom
    d['natoms'] = natoms
    # Load status
    with open(os.path.join(case_dir, 'status.json')) as f:
        d['status'] = json.load(f)
    return d

def hessian_ha_bohr_to_ev_ang2(H):
    """Convert Hessian from Hartree/Bohr^2 to eV/Angstrom^2."""
    return H * HARTREE_PER_BOHR2_TO_EV_PER_ANG2

def shortest_path_distances(bond_pairs, natoms):
    """BFS shortest-path distance between all pairs using the bond graph."""
    adj = [[] for _ in range(natoms)]
    for (i, j) in bond_pairs:
        adj[i].append(j)
        adj[j].append(i)
    dist = np.full((natoms, natoms), -1, dtype=int)
    for s in range(natoms):
        dist[s, s] = 0
        q = deque([s])
        while q:
            u = q.popleft()
            for v in adj[u]:
                if dist[s, v] == -1:
                    dist[s, v] = dist[s, u] + 1
                    q.append(v)
    return dist


def build_3rd_neighbor_bonds(symbols, positions, bond_pairs, max_dist=None):
    """Find atom pairs separated by exactly 3 bonds (1-4) and add a distance spring.

    These are not direct bonds (1-2) or angle pairs (1-3); they act like a
    simplified dihedral/cross-term by restraining the distance across four atoms.
    """
    natoms = len(symbols)
    dist = shortest_path_distances(bond_pairs, natoms)
    bonds3 = []
    for i in range(natoms):
        for j in range(i + 1, natoms):
            if dist[i, j] == 3:
                r = np.linalg.norm(positions[j] - positions[i])
                if max_dist is None or r < max_dist:
                    bonds3.append((i, j, r))
    return bonds3


def dihedral_energy_gradient(pos, d=1, n=3):
    """Energy and Cartesian gradient for a UFF/Prokop torsion term.

    Energy: E = V * (1 + d * cos(n * phi))   (V=1 here, so it returns dE/dV)
    Atoms: pos[0]=p1, pos[1]=p2, pos[2]=p3, pos[3]=p4  (i-j-k-l dihedral)
    Returns E, grad (4, 3) where grad = dE/dpos.
    This mirrors evalDihedral_Prokop() in cpp/common/molecular/UFF.h with
    bSubNonBond=False (no non-bonded subtraction).
    """
    p1, p2, p3, p4 = pos
    q12 = p1 - p2   # ji
    q32 = p3 - p2   # jk
    q43 = p4 - p3   # kl
    q12w = np.linalg.norm(q12)
    q32w = np.linalg.norm(q32)
    q43w = np.linalg.norm(q43)
    if q12w < 1e-12 or q32w < 1e-12 or q43w < 1e-12:
        return 0.0, np.zeros_like(pos)
    u12 = q12 / q12w
    u32 = q32 / q32w
    u43 = q43 / q43w

    n123 = np.cross(u12, u32)
    n234 = np.cross(u43, u32)
    nn123 = np.dot(n123, n123)
    nn234 = np.dot(n234, n234)
    if nn123 < 1e-14 or nn234 < 1e-14:
        return 0.0, np.zeros_like(pos)

    il2_123 = 1.0 / nn123
    il2_234 = 1.0 / nn234
    inv_n12 = np.sqrt(il2_123 * il2_234)

    csx = np.dot(n123, n234) * inv_n12
    csy = -np.dot(n123, u43) * inv_n12
    # cs = (cos phi, -sin phi)  => complex conjugate of e^{i phi}
    phi = np.arctan2(-csy, csx)
    nphi = n * phi
    csnx = np.cos(nphi)
    csny = -np.sin(nphi)

    E = 1.0 + d * csnx
    # f scalar used to build forces (f = -V*d*n*csn.y  with V=1)
    f = -d * n * csny

    fp1 = -f * n123 * il2_123 * q12w
    fp4 =  f * n234 * il2_234 * q43w
    c123 = np.dot(u32, u12) * (q32w / q12w)
    c432 = np.dot(u32, u43) * (q32w / q43w)
    fp3 = -c123 * fp1 - (c432 + 1.0) * fp4
    fp2 = (c123 - 1.0) * fp1 + c432 * fp4

    grad = -np.array([fp1, fp2, fp3, fp4])
    return E, grad


def dihedral_hessian(pos, d=1, n=3, h=1e-5):
    """Cartesian Hessian (12,12) for one torsion term via central finite differences.

    The Hessian is symmetric and is the derivative of dE/dpos with respect to pos.
    """
    pos = np.asarray(pos, dtype=float)
    n12 = pos.size
    gfun = lambda p: dihedral_energy_gradient(p, d, n)[1].ravel()
    H = np.zeros((n12, n12))
    for c in range(n12):
        posp = pos.copy(); posp.flat[c] += h
        posm = pos.copy(); posm.flat[c] -= h
        H[:, c] = (gfun(posp) - gfun(posm)) / (2.0 * h)
    H = (H + H.T) * 0.5
    return H


def dihedral_type_key(symbols, i, j, k, l):
    """Type key for a dihedral (i-j-k-l): central bond + sorted outer atoms."""
    return (tuple(sorted([symbols[j], symbols[k]])),
            tuple(sorted([symbols[i], symbols[l]])))


def build_dihedrals(symbols, positions, bonds, d=1, n=3, dihedral=False):
    """Build proper torsion list from bond topology (i-j-k-l).

    For each central bond j-k (taken once), enumerate all neighbours i of j and
    l of k to form the dihedral i-j-k-l.  d and n are the UFF phase/periodicity
    (default d=1, n=3 for sp3 tetrahedral torsions).
    """
    if not dihedral:
        return []
    natoms = len(symbols)
    neighs = [[] for _ in range(natoms)]
    bond_pairs = [(i, j) for (i, j, l0) in bonds]
    for (i, j) in bond_pairs:
        neighs[i].append(j)
        neighs[j].append(i)

    dihedrals = []
    for (j, k) in bond_pairs:
        if j > k:
            continue  # take each central bond only once
        for i in neighs[j]:
            if i == k:
                continue
            for l in neighs[k]:
                if l == j or l == i:
                    continue
                dihedrals.append((i, j, k, l, d, n))
    return dihedrals


def compute_dihedral_sensitivity(positions, symbols, dihedrals, global_dihedral_map, h=1e-5):
    """Compute full (3N,3N) sensitivity A_p = dH/dV for each dihedral type.

    Each A_p is the sum of the per-dihedral Hessians (with V=1) belonging to
    that type.  Returns a dict mapping parameter index to the sparse full matrix.
    """
    natoms = len(symbols)
    n3 = natoms * 3
    A = {p: np.zeros((n3, n3)) for p in global_dihedral_map.values()}
    for (i, j, k, l, d, n) in dihedrals:
        key = dihedral_type_key(symbols, i, j, k, l)
        p = global_dihedral_map[key]
        pos = positions[[i, j, k, l]]
        H = dihedral_hessian(pos, d, n, h=h)
        idx = np.array([i*3, i*3+1, i*3+2,
                        j*3, j*3+1, j*3+2,
                        k*3, k*3+1, k*3+2,
                        l*3, l*3+1, l*3+2])
        A[p][np.ix_(idx, idx)] += H
    return A


def add_dihedral_hessian(H, k, dihedral_A):
    """Add dihedral contribution to H in-place using precomputed A_p matrices."""
    for p, A in dihedral_A.items():
        H += k[p] * A


def build_topology(symbols, positions, bond_cutoff=1.85, third_bonds=False, third_bond_cutoff=None):
    """Build bond list and angle list from geometry using distance cutoff.
    
    For Si-H systems: Si-Si bonds ~2.35A, Si-H bonds ~1.48A
    Cutoff 1.85A captures both but not H-H (no bonds between H atoms).
    
    Returns:
        bonds: list of (i, j, l0)
        angles: list of (i, j, k, theta0) where j is central
        bonds3: list of (i, j, l0) for 3rd-neighbor 1-4 springs (empty if disabled)
    """
    natoms = len(symbols)
    # Find bonds by distance
    bonds = []
    bond_pairs = []
    for i in range(natoms):
        for j in range(i+1, natoms):
            if symbols[i] == 'H' and symbols[j] == 'H':
                continue  # no H-H bonds
            r = np.linalg.norm(positions[j] - positions[i])
            if r < bond_cutoff:
                bonds.append((i, j, r))  # l0 = current length (relaxed)
                bond_pairs.append((i, j))
    
    # Build neighbor list
    neighs = [[] for _ in range(natoms)]
    for (i, j) in bond_pairs:
        neighs[i].append(j)
        neighs[j].append(i)
    
    # Find angles: for each central atom j, all pairs (i, k) of neighbors
    angles = []
    for j in range(natoms):
        nn = neighs[j]
        for ii in range(len(nn)):
            for kk in range(ii+1, len(nn)):
                i = nn[ii]
                k = nn[kk]
                # Compute angle i-j-k
                ri = positions[i] - positions[j]
                rk = positions[k] - positions[j]
                ri /= np.linalg.norm(ri)
                rk /= np.linalg.norm(rk)
                cos_theta = np.dot(ri, rk)
                theta0 = np.arccos(np.clip(cos_theta, -1, 1))
                angles.append((i, j, k, theta0))

    bonds3 = []
    if third_bonds:
        bonds3 = build_3rd_neighbor_bonds(symbols, positions, bond_pairs, max_dist=third_bond_cutoff)
    
    return bonds, angles, bonds3

def assign_types(symbols, bonds, angles):
    """Assign bond/angle types based on element pairs/triples.
    
    For Si-H nanoparticles we expect:
      Bond types: Si-Si, Si-H
      Angle types: Si-Si-Si, Si-Si-H, H-Si-H, Si-Si-Si (same), etc.
    """
    # Bond types: sorted element pair -> type index
    bond_type_map = {}
    bond_types = []
    for (i, j, l0) in bonds:
        key = tuple(sorted([symbols[i], symbols[j]]))
        if key not in bond_type_map:
            bond_type_map[key] = len(bond_type_map)
        bond_types.append(bond_type_map[key])
    
    # Angle types: (element_i, element_j, element_k) where j is central.
    # Key is symmetric under i<->k: outer atoms are sorted.
    angle_type_map = {}
    angle_types = []
    for (i, j, k, theta0) in angles:
        key = angle_type_key(symbols, i, j, k)
        if key not in angle_type_map:
            angle_type_map[key] = len(angle_type_map)
        angle_types.append(angle_type_map[key])
    
    return bond_types, angle_types, bond_type_map, angle_type_map

class ParamMap:
    """Sparse mapping from interaction terms to free parameters (symmetry constraints).

    Mirrors C++ ParamMap struct in FFfit.h.

    PRINCIPLE OF TRANSFERABILITY:
        Force-field parameters describe chemical bond TYPES, not individual bonds.
        All Si-Si bonds share one stiffness k_SiSi, all H-Si-H angles share one K_HSiH.
        This reduces the number of fitted unknowns from O(n_bonds) to O(n_bond_types).

    WHY IT MATTERS:
        - Without sharing: 152 Si-H bonds → 152 parameters (underdetermined)
        - With sharing:    152 Si-H bonds → 1 parameter (well-constrained)
        - Multi-system: same k applies across ALL systems, so constraints accumulate:
          G_total = Σ_sys G_sys,  y_total = Σ_sys y_sys

    MAPPING CRITERIA (extensible):
        - Element types (current): Si-Si, Si-H, H-Si-H, Si-Si-H, etc.
        - Chemical environment (future): Si-Si with 4 vs 3 neighbors
        - Manual: custom symmetry groups

    Each bond/angle term maps to one free parameter index.
    Multiple terms can share the same parameter (e.g. all Si-H bonds → one k).
    """
    def __init__(self, nbonds, nangles):
        self.bond_param_idx = np.full(nbonds, -1, dtype=int)
        self.angle_param_idx = np.full(nangles, -1, dtype=int)
        self.n_free = 0
    
    @classmethod
    def from_element_types(cls, bonds, angles, symbols):
        """Auto-assign: same element pair/triple → same parameter."""
        pm = cls(len(bonds), len(angles))
        # Bonds: sorted element pair
        bmap = {}
        for ib, (i, j, l0) in enumerate(bonds):
            key = tuple(sorted([symbols[i], symbols[j]]))
            if key not in bmap:
                bmap[key] = pm.n_free; pm.n_free += 1
            pm.bond_param_idx[ib] = bmap[key]
        # Angles: element triple (i, j_central, k). Outer atoms sorted.
        amap = {}
        for ia, (i, j, k, theta0) in enumerate(angles):
            key = angle_type_key(symbols, i, j, k)
            if key not in amap:
                amap[key] = pm.n_free; pm.n_free += 1
            pm.angle_param_idx[ia] = amap[key]
        return pm, bmap, amap
    
    def set_bond_param(self, ibond, param_idx):
        if len(self.bond_param_idx) <= ibond:
            self.bond_param_idx = np.pad(self.bond_param_idx, ibond+1, constant_values=-1)
        self.bond_param_idx[ibond] = param_idx
        self.n_free = max(self.n_free, param_idx + 1)
    
    def set_angle_param(self, iangle, param_idx):
        if len(self.angle_param_idx) <= iangle:
            self.angle_param_idx = np.pad(self.angle_param_idx, iangle+1, constant_values=-1)
        self.angle_param_idx[iangle] = param_idx
        self.n_free = max(self.n_free, param_idx + 1)

def build_sensitivity_matrices(positions, bonds, angles, bond_types, angle_types, natoms):
    """Build dH/dk_t sensitivity matrices in Python (mirrors FFfit.h C++ logic).
    
    bond_types and angle_types are GLOBAL parameter indices (from ParamMap).
    Returns A: list of (3N x 3N) numpy arrays, one per free parameter.
    """
    n3 = natoms * 3
    nparams = max(max(bond_types) + 1 if bond_types else 0,
                  max(angle_types) + 1 if angle_types else 0)
    A = [np.zeros((n3, n3)) for _ in range(nparams)]
    
    # Bond contributions
    for ib, (i, j, l0) in enumerate(bonds):
        t = bond_types[ib]
        d = positions[j] - positions[i]
        r = np.linalg.norm(d)
        e = d / r
        dl = r - l0
        # F'' = 1, F' = dl
        # b_i = -e, b_j = +e
        # dH/dk = b⊗b^T + dl * C
        # C_ii = P/r, C_jj = P/r, C_ij = -P/r, P = I - e⊗e
        P = np.eye(3) - np.outer(e, e)
        inv_r = 1.0 / r
        Cii = P * inv_r
        Cjj = P * inv_r
        Cij = -P * inv_r
        # Blocks
        ii_idx = slice(i*3, i*3+3)
        jj_idx = slice(j*3, j*3+3)
        # Rank-one: b⊗b^T = e⊗e for (i,i) and (j,j), -e⊗e for (i,j) and (j,i)
        eeT = np.outer(e, e)
        A[t][ii_idx, ii_idx] += eeT + dl * Cii
        A[t][ii_idx, jj_idx] += -eeT + dl * Cij
        A[t][jj_idx, ii_idx] += -eeT + dl * Cij.T
        A[t][jj_idx, jj_idx] += eeT + dl * Cjj
    
    # Angle contributions: E = 0.5*kθ*(c-c0)^2/(1-c0^2), so kθ is angular curvature in eV/rad^2
    for ia, (i, j, k, theta0) in enumerate(angles):
        t = angle_types[ia]  # already global index from ParamMap
        a_vec = positions[i] - positions[j]
        c_vec = positions[k] - positions[j]
        al = np.linalg.norm(a_vec)
        cl = np.linalg.norm(c_vec)
        u = a_vec / al
        v = c_vec / cl
        cos_theta = np.dot(u, v)
        c0 = np.cos(theta0)
        s02 = 1.0 - c0*c0
        assert s02 > 1e-14
        gi = (v - cos_theta*u) / al
        gk = (u - cos_theta*v) / cl
        gj = -(gi + gk)
        scale = 1.0 / s02
        vecs = {i: gi, j: gj, k: gk}
        for a_atom in [i, j, k]:
            for b_atom in [i, j, k]:
                a_idx = slice(a_atom*3, a_atom*3+3)
                b_idx = slice(b_atom*3, b_atom*3+3)
                A[t][a_idx, b_idx] += scale * np.outer(vecs[a_atom], vecs[b_atom])
    
    return A

def accumulate_normal_equations(G, y, H_ref, A, weight=None):
    """Add one system's contribution to normal equations G, y (in-place)."""
    nparams = len(A)
    n3 = H_ref.shape[0]
    if weight is None:
        weight = np.ones(n3)
    for p in range(nparams):
        WAp = weight[:, None] * A[p]
        for q in range(p, nparams):
            WAq = weight[:, None] * A[q]
            val = np.sum(WAp * WAq)
            G[p, q] += val
            if q != p: G[q, p] += val
        y[p] += np.sum(WAp * H_ref)

def fit_hessian(H_ref, A, weight=None):
    """Solve linear least-squares for a single system."""
    nparams = len(A)
    G = np.zeros((nparams, nparams))
    y = np.zeros(nparams)
    accumulate_normal_equations(G, y, H_ref, A, weight)
    k = np.linalg.lstsq(G, y, rcond=None)[0]
    return k

def compute_model_hessian(A, k):
    """Compute model Hessian from fitted parameters."""
    H = np.zeros_like(A[0])
    for t in range(len(A)):
        H += k[t] * A[t]
    return H

# ============================================================================
# HESSIAN-BASED FF FITTING — Python reference implementation
# ============================================================================
# Theory (mirrors FFfit.h header):
#
# Energy: E(r) = Σ_t k_t * f_t(q_t(r))
# Hessian: H_αβ = Σ_t k_t * A_t,  where A_t = ∂²f_t/∂r_α ∂r_β
#
# Chain rule for f(q(r)):
#   A_αβ = f''(q) * b_α b_β + f'(q) * C_αβ
#   b_α = ∂q/∂r_α (Wilson vector),  C_αβ = ∂²q/∂r_α ∂r_β (coordinate Hessian)
#
# Bond: f(r) = ½(r-r0)², f'=dl, f''=1
#   A = b⊗b^T + dl * C_r,  at eq: A = e⊗e^T
#
# Angle (normalized cos): E = 0.5*kθ*(c-c0)²/(1-c0²), g=∂c/∂r
#   A = (g⊗g^T + (c-c0)*C_cos)/(1-c0²),  at eq: A = bθ⊗bθ and kθ is eV/rad²
# ============================================================================

def compute_gradient_term_blocks(positions, bonds, angles, param_map, k, H_ref, H0=None, weight=None):
    """Compute gradient of loss = ||H_model - H_ref||^2_W w.r.t. free parameters.

    Uses per-term sensitivity blocks (no full 3N×3N sensitivity matrices needed).
    Mirrors C++ FFfit::compute_gradient.
    """
    natoms = positions.shape[0]
    n3 = natoms * 3
    np_free = param_map.n_free
    # Forward pass: build H_model = Σ_t k_t * A_t
    H_model = np.zeros((n3, n3))
    # --- Bond sensitivity blocks: A = e⊗e^T + dl * C_r ---
    bond_blocks = []
    for ib, (i, j, l0) in enumerate(bonds):
        d = positions[j] - positions[i]
        r = np.linalg.norm(d)
        e = d / r
        dl = r - l0  # f'(r) = (r - r0)
        P = np.eye(3) - np.outer(e, e)  # projector ⊥ to bond
        inv_r = 1.0 / r
        eeT = np.outer(e, e)  # b⊗b^T = e⊗e^T (rank-one part, f''=1)
        # C_r: C_ii = P/r, C_jj = P/r, C_ij = -P/r
        # A_ii = e⊗e^T + dl*(P/r),  A_ij = -e⊗e^T + dl*(-P/r)
        Bii = eeT + dl * P * inv_r
        Bij = -eeT + dl * (-P * inv_r)
        Bji = Bij.T
        Bjj = eeT + dl * P * inv_r
        bond_blocks.append((Bii, Bij, Bji, Bjj, i, j))
        f = param_map.bond_param_idx[ib]
        if f >= 0:
            H_model[i*3:i*3+3, i*3:i*3+3] += k[f] * Bii
            H_model[i*3:i*3+3, j*3:j*3+3] += k[f] * Bij
            H_model[j*3:j*3+3, i*3:i*3+3] += k[f] * Bji
            H_model[j*3:j*3+3, j*3:j*3+3] += k[f] * Bjj
    # --- Angle sensitivity blocks (UFF-normalized cosine form) ---
    # E = 0.5*kθ*(c-c0)^2/(1-c0^2), g=∂c/∂r
    # A = (g⊗g^T + (c-c0)*C_cos)/(1-c0^2), so kθ is actual local angular curvature
    angle_blocks = []
    for ia, (i, j, k_atom, theta0) in enumerate(angles):
        a_vec = positions[i] - positions[j]
        c_vec = positions[k_atom] - positions[j]
        al = np.linalg.norm(a_vec)
        cl = np.linalg.norm(c_vec)
        u = a_vec / al; v = c_vec / cl
        cos_t = np.dot(u, v)  # c = cos θ
        gi = (v - cos_t*u) / al
        gk = (u - cos_t*v) / cl
        gj = -(gi + gk)
        c0 = np.cos(theta0)
        s02 = 1.0 - c0*c0
        assert s02 > 1e-14
        dc = cos_t - c0
        scale = 1.0 / s02
        vecs = [gi, gj, gk]
        atoms = [i, j, k_atom]
        # Rank-one part: B[a][b] = (g_a ⊗ g_b)/(1-c0²)
        B = [[scale * np.outer(vecs[a], vecs[b]) for b in range(3)] for a in range(3)]
        # Prestress: B[a][b] += Fp * C_cos[a][b]  (only if c ≠ c0)
        # C_cos = ∂²(cosθ)/∂r², analytical blocks (see FFfit.h for derivation):
        #   C_ii = (-v⊗u - u⊗v + 3*ct*u⊗u - ct*I) / |a|²
        #   C_kk = (-u⊗v - v⊗u + 3*ct*v⊗v - ct*I) / |c|²
        #   C_ik = (I - v⊗v - u⊗u + ct*u⊗v) / (|a|*|c|)
        #   C_ij = -C_ii - C_ik
        #   C_jk = -C_ik - C_kk   [BUG FIX: use C_ik, not C_ki]
        if abs(dc) > 1e-12:
            I3 = np.eye(3)
            uuT = np.outer(u, u); vvT = np.outer(v, v)
            vuT = np.outer(v, u); uvT = np.outer(u, v)
            C_ii = (-vuT - uvT + 3*cos_t*uuT - cos_t*I3) / (al*al)
            C_kk = (-uvT - vuT + 3*cos_t*vvT - cos_t*I3) / (cl*cl)
            C_ik = (I3 - vvT - uuT + cos_t*uvT) / (al*cl)
            C_ki = C_ik.T
            C_ij = -C_ii - C_ik; C_ji = C_ij.T
            C_jk = -C_ik - C_kk; C_kj = C_jk.T
            C_jj = C_ii + C_ik + C_ki + C_kk
            Ccos = [[C_ii, C_ij, C_ik], [C_ji, C_jj, C_jk], [C_ki, C_kj, C_kk]]
            for a in range(3):
                for b in range(3):
                    B[a][b] += scale * dc * Ccos[a][b]
        angle_blocks.append((B, atoms))
        f = param_map.angle_param_idx[ia]
        if f >= 0:
            for a in range(3):
                for b in range(3):
                    H_model[atoms[a]*3:atoms[a]*3+3, atoms[b]*3:atoms[b]*3+3] += k[f] * B[a][b]
    # Residual
    h0 = H0 if H0 is not None else 0.0
    dH = H_model - H_ref + h0
    if weight is not None:
        W2 = weight * weight  # weight is per-element (flattened)
        W2_mat = W2.reshape(n3, n3) if W2.size == n3*n3 else np.ones((n3, n3))
    else:
        W2_mat = np.ones((n3, n3))
    loss = np.sum(W2_mat * dH * dH)
    # Backward pass: gradient[f] = 2 * Σ_{term→f} <A_term, dH>_W
    grad = np.zeros(np_free)
    for ib, (Bii, Bij, Bji, Bjj, i, j) in enumerate(bond_blocks):
        f = param_map.bond_param_idx[ib]
        if f < 0: continue
        ii = slice(i*3, i*3+3); jj = slice(j*3, j*3+3)
        grad[f] += 2.0 * (np.sum(W2_mat[ii,ii] * Bii * dH[ii,ii]) +
                          np.sum(W2_mat[ii,jj] * Bij * dH[ii,jj]) +
                          np.sum(W2_mat[jj,ii] * Bji * dH[jj,ii]) +
                          np.sum(W2_mat[jj,jj] * Bjj * dH[jj,jj]))
    for ia, (B, atoms) in enumerate(angle_blocks):
        f = param_map.angle_param_idx[ia]
        if f < 0: continue
        for a in range(3):
            for b in range(3):
                ai = slice(atoms[a]*3, atoms[a]*3+3)
                bi_s = slice(atoms[b]*3, atoms[b]*3+3)
                grad[f] += 2.0 * np.sum(W2_mat[ai,bi_s] * B[a][b] * dH[ai,bi_s])
    return grad, loss

def fit_gradient_descent(positions, bonds, angles, param_map, H_ref, H0=None, weight=None,
                         lr=1e-3, momentum=0.9, max_iter=1000, tol=1e-8, verbose=True):
    """Fit via gradient descent with momentum (single system)."""
    np_free = param_map.n_free
    k = np.ones(np_free)
    velocity = np.zeros(np_free)
    prev_loss = 1e30
    for it in range(max_iter):
        grad, loss = compute_gradient_term_blocks(positions, bonds, angles, param_map, k, H_ref, H0, weight)
        grad_norm = np.linalg.norm(grad)
        if verbose and (it % 100 == 0 or it < 10):
            print(f"  GD iter {it:4d}: loss={loss:.6e} grad_norm={grad_norm:.6e}")
        if grad_norm < tol or (it > 0 and abs(prev_loss - loss) < tol * prev_loss):
            if verbose: print(f"  GD converged at iter {it}: loss={loss:.6e}")
            break
        prev_loss = loss
        velocity = momentum * velocity - lr * grad
        k += velocity
    return k

def fit_gradient_descent_multi(systems, lr=1e-4, momentum=0.9, max_iter=5000, tol=1e-10, verbose=True):
    """Fit via gradient descent across multiple systems simultaneously.
    
    systems: list of dicts with keys: positions, bonds, angles, param_map, H_ref, weight
    """
    np_free = systems[0]['param_map'].n_free
    k = np.ones(np_free)
    velocity = np.zeros(np_free)
    prev_loss = 1e30
    for it in range(max_iter):
        grad = np.zeros(np_free)
        total_loss = 0.0
        for sys in systems:
            g, l = compute_gradient_term_blocks(sys['positions'], sys['bonds'], sys['angles'],
                                                 sys['param_map'], k, sys['H_ref'], None, sys.get('weight'))
            grad += g
            total_loss += l
        grad_norm = np.linalg.norm(grad)
        if verbose and (it % 100 == 0 or it < 10):
            print(f"  GD iter {it:4d}: loss={total_loss:.6e} grad_norm={grad_norm:.6e}")
        if grad_norm < tol or (it > 0 and abs(prev_loss - total_loss) < tol * prev_loss):
            if verbose: print(f"  GD converged at iter {it}: loss={total_loss:.6e}")
            break
        prev_loss = total_loss
        velocity = momentum * velocity - lr * grad
        k += velocity
    return k

def build_global_param_map(all_bonds, all_angles, all_symbols, all_bonds3=None, all_dihedrals=None):
    """Build a global ParamMap from all systems' element types.
    
    all_bonds/all_angles: list of lists (one per system)
    all_symbols: list of symbol lists (one per system)
    all_bonds3: optional list of 3rd-neighbor (1-4) bond lists (one per system)
    all_dihedrals: optional list of dihedral tuples (one per system)
    Returns: (global_bond_type_map, global_bond3_type_map, global_angle_type_map, global_dihedral_map, n_free)
    """
    all_bonds3 = all_bonds3 or []
    all_dihedrals = all_dihedrals or []
    global_bond_map = {}
    global_bond3_map = {}
    global_angle_map = {}
    global_dihedral_map = {}
    # Pass 1: discover all 1-2 bond types
    for bonds, symbols in zip(all_bonds, all_symbols):
        for (i, j, l0) in bonds:
            key = tuple(sorted([symbols[i], symbols[j]]))
            if key not in global_bond_map:
                global_bond_map[key] = len(global_bond_map)
    # Pass 2: discover all 3rd-neighbor (1-4) bond types, offset after 1-2 bonds
    n_bond_types = len(global_bond_map)
    for bonds3, symbols in zip(all_bonds3, all_symbols):
        for (i, j, l0) in bonds3:
            key = tuple(sorted([symbols[i], symbols[j]]))
            if key not in global_bond3_map:
                global_bond3_map[key] = n_bond_types + len(global_bond3_map)
    # Pass 3: discover all angle types, offset after bond types
    n_bond3_types = len(global_bond3_map)
    for angles, symbols in zip(all_angles, all_symbols):
        for (i, j, k, theta0) in angles:
            key = angle_type_key(symbols, i, j, k)
            if key not in global_angle_map:
                global_angle_map[key] = n_bond_types + n_bond3_types + len(global_angle_map)
    # Pass 4: discover all dihedral types, offset after angle types
    n_angle_offset = n_bond_types + n_bond3_types
    for dihedrals, symbols in zip(all_dihedrals, all_symbols):
        for (i, j, k, l, d, n) in dihedrals:
            key = dihedral_type_key(symbols, i, j, k, l)
            if key not in global_dihedral_map:
                global_dihedral_map[key] = n_angle_offset + len(global_dihedral_map)
    n_free = n_angle_offset + len(global_dihedral_map)
    return global_bond_map, global_bond3_map, global_angle_map, global_dihedral_map, n_free

def compute_averaged_equilibrium(all_bonds, all_angles, all_symbols, all_positions,
                                 global_bond_map, global_angle_map, global_bond3_map=None, all_bonds3=None):
    """Compute averaged equilibrium bond lengths l0 and angle cosine c0 across all systems.

    The angle force is written in c=cos(theta), therefore the transferable equilibrium
    coordinate is c0=<cos(theta)>, not cos(<theta>). We store theta0=acos(c0) only
    because the C++ API currently stores AngleDef.theta0 in radians.
    """
    global_bond3_map = global_bond3_map or {}
    all_bonds3 = all_bonds3 or []
    bond_lengths = {}
    bond3_lengths = {}
    angle_cos_values = {}
    angle_theta_values = {}

    for bonds, angles, symbols, positions in zip(all_bonds, all_angles, all_symbols, all_positions):
        for (i, j, l0) in bonds:
            key = tuple(sorted([symbols[i], symbols[j]]))
            r = np.linalg.norm(positions[j] - positions[i])
            bond_lengths.setdefault(key, []).append(r)
        for (i, j, k, theta0) in angles:
            key = angle_type_key(symbols, i, j, k)
            ri = positions[i] - positions[j]
            rk = positions[k] - positions[j]
            ri /= np.linalg.norm(ri)
            rk /= np.linalg.norm(rk)
            cos_t = np.clip(np.dot(ri, rk), -1.0, 1.0)
            angle_cos_values.setdefault(key, []).append(cos_t)
            angle_theta_values.setdefault(key, []).append(np.arccos(cos_t))

    for bonds3, symbols, positions in zip(all_bonds3, all_symbols, all_positions):
        for (i, j, l0) in bonds3:
            key = tuple(sorted([symbols[i], symbols[j]]))
            r = np.linalg.norm(positions[j] - positions[i])
            bond3_lengths.setdefault(key, []).append(r)

    avg_l0 = {}
    for key, lengths in bond_lengths.items():
        avg_l0[key] = np.mean(lengths)
        print(f"  l0_avg[{key}] = {avg_l0[key]:.6f} A  (from {len(lengths)} bonds, std={np.std(lengths):.4f})")

    avg_l0_3 = {}
    for key, lengths in bond3_lengths.items():
        avg_l0_3[key] = np.mean(lengths)
        print(f"  l0_3_avg[{key}] = {avg_l0_3[key]:.6f} A  (from {len(lengths)} 3rd-neighbor bonds, std={np.std(lengths):.4f})")

    avg_theta0 = {}
    for key, cos_vals in angle_cos_values.items():
        c0 = np.clip(np.mean(cos_vals), -1.0, 1.0)
        avg_theta0[key] = np.arccos(c0)
        thetas = np.array(angle_theta_values[key])
        print(f"  c0_avg[{key}] = {c0:.8f}  theta0={np.degrees(avg_theta0[key]):.4f} deg  (from {len(cos_vals)} angles, theta_std={np.degrees(np.std(thetas)):.4f} deg)")

    return avg_l0, avg_l0_3, avg_theta0

def rebuild_topology_with_averaged(bonds, angles, bonds3, symbols, positions, avg_l0, avg_theta0, avg_l0_3=None):
    """Rebuild bond/angle lists with averaged equilibrium parameters.
    
    Returns new bonds/angles/bonds3 lists with l0 and theta0 replaced by averaged values.
    """
    new_bonds = []
    for (i, j, l0) in bonds:
        key = tuple(sorted([symbols[i], symbols[j]]))
        new_l0 = avg_l0[key]
        new_bonds.append((i, j, new_l0))

    new_angles = []
    for (i, j, k, theta0) in angles:
        key = angle_type_key(symbols, i, j, k)
        new_theta0 = avg_theta0[key]
        new_angles.append((i, j, k, new_theta0))

    new_bonds3 = []
    if avg_l0_3 is not None:
        for (i, j, l0) in bonds3:
            key = tuple(sorted([symbols[i], symbols[j]]))
            new_l0 = avg_l0_3[key]
            new_bonds3.append((i, j, new_l0))

    return new_bonds, new_angles, new_bonds3

def make_param_map_for_system(bonds, angles, symbols, global_bond_map, global_angle_map, bonds3=None, global_bond3_map=None):
    """Create a per-system ParamMap using global type indices."""
    bonds3 = bonds3 or []
    global_bond3_map = global_bond3_map or {}
    pm = ParamMap(len(bonds) + len(bonds3), len(angles))
    for ib, (i, j, l0) in enumerate(bonds):
        key = tuple(sorted([symbols[i], symbols[j]]))
        pm.bond_param_idx[ib] = global_bond_map[key]
    for ib3, (i, j, l0) in enumerate(bonds3):
        key = tuple(sorted([symbols[i], symbols[j]]))
        pm.bond_param_idx[len(bonds) + ib3] = global_bond3_map[key]
    for ia, (i, j, k, theta0) in enumerate(angles):
        key = angle_type_key(symbols, i, j, k)
        pm.angle_param_idx[ia] = global_angle_map[key]
    pm.n_free = len(global_bond_map) + len(global_bond3_map) + len(global_angle_map)
    return pm

def get_reference_modes_and_freqs(H_ref, masses, data=None, freq_floor_cm1=10.0):
    """Return mass-weighted eigenvectors (3N x n_modes) and eigenvalues (λ = ω²).

    H_ref: (3N, 3N) eV/Å² (fallback if PySCF modes not available)
    masses: (N,) amu
    data: dict with optional 'modes' (nmodes, natoms, 3) and 'freqs' (nmodes,) in cm^-1
    freq_floor_cm1: modes with |freq| < floor are skipped (translations/rotations)

    Returns:
        V: (3N, n_modes)  mass-weighted eigenvectors, normalized to 1
        lambdas: (n_modes,)  eV/(Å² amu)
        freqs_cm1: (n_modes,)
    """
    conv = 521.5  # eV/(A^2 amu) -> cm^-1
    if data is not None and 'modes' in data and 'freqs' in data:
        modes = data['modes']  # (nmodes, natoms, 3)  Cartesian displacement
        freqs = data['freqs']  # (nmodes,) cm^-1  (may be complex for imaginary modes)
        n_modes = modes.shape[0]
        n3 = len(masses) * 3
        sqrt_m = np.sqrt(np.repeat(masses, 3))
        V = (modes.reshape(n_modes, n3) * sqrt_m[None, :]).T  # (3N, n_modes)
        # Normalise (mass-weighted eigenvectors)
        norms = np.linalg.norm(V, axis=0)
        V = V / norms[None, :]
        # Use real part for real modes; pure imaginary modes are translations/rotations
        freqs_real = np.real(freqs)
        mask = freqs_real > freq_floor_cm1
        if not np.any(mask):
            mask = freqs_real > 0
        lam = (freqs_real[mask] / conv)**2
        return V[:, mask], lam, freqs_real[mask]
    # Fallback: derive from H_ref (may contain spurious translation/rotation modes)
    n3 = H_ref.shape[0]
    inv_sqrt_m = 1.0 / np.sqrt(np.repeat(masses, 3))
    D = (inv_sqrt_m[:, None] * inv_sqrt_m[None, :]) * H_ref
    lam, V = np.linalg.eigh(D)
    freqs = conv * np.sqrt(np.maximum(lam, 0))
    mask = freqs > freq_floor_cm1
    if not np.any(mask):
        mask = lam > 0
    return V[:, mask], lam[mask], freqs[mask]


def get_acoustic_projector(positions, masses):
    """Build projector P = I - V_ac V_ac^T onto vibrational subspace.

    V_ac is the orthonormal (3N, 6) matrix of mass-weighted translation and
    rotation vectors. Removing these eliminates the 6 rigid-body modes from
    the model Hessian so only internal vibrational frequencies are compared.
    """
    natoms = len(masses)
    n3 = natoms * 3
    sqrt_m = np.sqrt(np.repeat(masses, 3))
    # Center of mass
    cm = np.sum(positions * masses[:, None], axis=0) / np.sum(masses)
    r = positions - cm
    # Translation vectors (3 vectors)
    T = np.zeros((n3, 3))
    for a in range(3):
        T[a::3, a] = sqrt_m[a::3]
    # Rotation vectors (3 vectors)
    R = np.zeros((n3, 3))
    e = np.eye(3)
    for a in range(3):
        for i in range(natoms):
            v = np.cross(e[a], r[i])
            R[i*3:i*3+3, a] = v * np.sqrt(masses[i])
    # Orthonormalize
    V_ac = np.hstack([T, R])
    # Gram-Schmidt
    V_ac = np.linalg.qr(V_ac)[0]
    P = np.eye(n3) - V_ac @ V_ac.T
    return P


def get_frequencies_cm1(H, masses, data=None, freq_floor=10.0, positions=None, project_acoustic=False):
    """Return positive vibrational frequencies (cm^-1) from a Hessian.
    
    If data['freqs'] is available, use it as reference (real part). Otherwise
    compute from H via mass-weighting. If positions is provided and
    project_acoustic is True, the model Hessian is projected onto the
    complement of the rigid-body subspace before diagonalization.
    """
    conv = 521.5
    if data is not None and 'freqs' in data:
        freqs_cm1 = np.real(data['freqs'])
    else:
        inv_sqrt_m = 1.0 / np.sqrt(np.repeat(masses, 3))
        D = (inv_sqrt_m[:, None] * inv_sqrt_m[None, :]) * H
        if project_acoustic and positions is not None:
            P = get_acoustic_projector(positions, masses)
            D = P @ D @ P
        lam = np.linalg.eigvalsh(D)
        freqs_cm1 = conv * np.sqrt(np.maximum(0, lam))
    return np.sort(freqs_cm1[freqs_cm1 > freq_floor])


def compare_frequencies(H_ref, H_model, masses, data=None, label="", freq_floor=10.0, positions=None):
    """Compare vibrational frequencies from reference and model Hessians.

    Reference frequencies are taken from data['freqs'] if available; otherwise
    they are derived from H_ref.  Model frequencies are projected onto the
    internal vibrational subspace (removing translation/rotation) if positions
    are supplied, then matched to each reference frequency by nearest neighbour.
    """
    nz_ref = get_frequencies_cm1(H_ref, masses, data=data, freq_floor=freq_floor)
    nz_model = get_frequencies_cm1(H_model, masses, data=None, freq_floor=freq_floor,
                                   positions=positions, project_acoustic=True)
    print(f"  [{label}] Reference freqs (> {freq_floor} cm^-1): {len(nz_ref)}")
    print(f"  [{label}] Model freqs (> {freq_floor} cm^-1): {len(nz_model)}")
    if len(nz_ref) > 0 and len(nz_model) > 0:
        # Match each reference frequency to the nearest model frequency
        diffs = np.min(np.abs(nz_ref[:, None] - nz_model[None, :]), axis=1)
        max_diff = np.max(diffs)
        print(f"  [{label}] First 10 ref:  {nz_ref[:10]}")
        print(f"  [{label}] First 10 model:{nz_model[:10]}")
        print(f"  [{label}] Max nearest-neighbour freq diff: {max_diff:.2f} cm^-1")
        print(f"  [{label}] Mean nearest-neighbour freq diff: {np.mean(diffs):.2f} cm^-1")
    return nz_ref, nz_model

def broaden(freqs, x, width=20.0):
    """Sum of normalized Gaussians centered at freqs, evaluated on x grid."""
    spec = np.zeros_like(x)
    if len(freqs) == 0:
        return spec
    for f in freqs:
        spec += np.exp(-((x - f)**2) / (2 * width**2))
    return spec


def plot_spectrum(freqs_ref, freqs_model, label, outdir=None, exp_spectrum=None,
                  width=20.0, xmax=2500, noshow=True):
    """Plot reference vs model vibrational spectrum (stick + broadened).

    If exp_spectrum is a (2, n) array or a path to a CSV with two columns
    (frequency, intensity), it is overlaid after normalizing to the model peak.
    """
    x = np.arange(0, xmax + 1, 0.5)
    spec_ref = broaden(freqs_ref, x, width)
    spec_model = broaden(freqs_model, x, width)
    norm = spec_ref.max() if spec_ref.max() > 0 else 1.0
    if spec_model.max() > 0:
        spec_model /= spec_model.max()  # normalize to own peak
    if spec_ref.max() > 0:
        spec_ref /= spec_ref.max()

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
    fig.suptitle(f'Vibrational spectrum — {label}', fontsize=13)

    # Stick panel
    ax1.vlines(freqs_ref, 0, 1, colors='#2563eb', linewidth=0.6, alpha=0.75, label='reference (PySCF)')
    ax1.vlines(freqs_model, 0, 0.7, colors='#dc2626', linewidth=0.6, alpha=0.75, label='model (FF)')
    ax1.set_ylim(0, 1.1)
    ax1.set_ylabel('Stick (arb.)')
    ax1.legend(fontsize=9, loc='upper right')
    ax1.grid(True, alpha=0.25, ls=':')
    ax1.tick_params(labelbottom=False)

    # Broadened panel
    ax2.plot(x, spec_ref, color='#1d4ed8', lw=1.4, label=f'Reference (σ={width:.0f} cm⁻¹)')
    ax2.fill_between(x, 0, spec_ref, color='#3b82f6', alpha=0.15)
    ax2.plot(x, spec_model, color='#b91c1c', lw=1.4, label=f'Model (σ={width:.0f} cm⁻¹)')
    ax2.fill_between(x, 0, spec_model, color='#f87171', alpha=0.12)

    # Optional experimental spectrum overlay
    exp_x = exp_y = None
    if exp_spectrum is not None:
        try:
            if isinstance(exp_spectrum, str):
                exp = np.loadtxt(exp_spectrum)
            else:
                exp = np.asarray(exp_spectrum)
            exp_x, exp_y = exp[:, 0], exp[:, 1]
            if exp_y.max() > 0:
                exp_y = exp_y / exp_y.max()
            ax2.plot(exp_x, exp_y, 'k-', lw=1.0, alpha=0.7, label='Experiment')
        except Exception as e:
            print(f"  WARNING: could not load experimental spectrum: {e}")

    ax2.set_xlim(0, xmax)
    ax2.set_xlabel('Frequency (cm⁻¹)')
    ax2.set_ylabel('Intensity (normalized)')
    ax2.legend(fontsize=9, loc='upper right')
    ax2.grid(True, alpha=0.25, ls=':')

    plt.tight_layout()
    if outdir:
        os.makedirs(outdir, exist_ok=True)
        out = os.path.join(outdir, f'{label}_spectrum.png')
        fig.savefig(out, dpi=150)
        print(f"  Saved spectrum plot: {out}")
    if not noshow:
        plt.show()
    plt.close(fig)
    return out if outdir else None


def plot_spectra_overlay(system_data, outdir=None, width=20.0, xmax=2500, noshow=True):
    """Overlay broadened reference/model spectra across multiple systems.

    system_data: list of (label, freqs_ref, freqs_model) tuples.
    """
    x = np.arange(0, xmax + 1, 0.5)
    fig, ax = plt.subplots(figsize=(13, 7))
    colors = plt.cm.get_cmap('tab10')
    for i, (label, freqs_ref, freqs_model) in enumerate(system_data):
        spec_ref = broaden(freqs_ref, x, width)
        spec_model = broaden(freqs_model, x, width)
        if spec_ref.max() > 0:
            spec_ref /= spec_ref.max()
        if spec_model.max() > 0:
            spec_model /= spec_model.max()
        c = colors(i % 10)
        ax.plot(x, spec_ref + i, color=c, lw=1.2, alpha=0.6, linestyle='-', label=f'{label} ref')
        ax.plot(x, spec_model + i, color=c, lw=1.2, alpha=0.9, linestyle='--', label=f'{label} model')
    ax.set_xlim(0, xmax)
    ax.set_xlabel('Frequency (cm⁻¹)')
    ax.set_ylabel('Normalized intensity (offset by system)')
    ax.set_title('Overlay of reference vs model vibrational spectra')
    ax.legend(fontsize=8, loc='upper right', ncol=2)
    ax.grid(True, alpha=0.25, ls=':')
    if outdir:
        os.makedirs(outdir, exist_ok=True)
        out = os.path.join(outdir, 'spectrum_overlay.png')
        fig.savefig(out, dpi=150)
        print(f"  Saved overlay plot: {out}")
    if not noshow:
        plt.show()
    plt.close(fig)


def get_sensitivity_action(fitter, V, masses, dihedral_A=None):
    """Compute D_p v_s for each free parameter p and each reference mode s.

    D(k) = M^{-1/2} H(k) M^{-1/2}.  For a one-hot parameter vector e_p,
    H(e_p) = dH/dk_p = A_p, so D_p = M^{-1/2} A_p M^{-1/2}.
    The action on the mass-weighted eigenvector v_s is
        D_p v_s = M^{-1/2} [ A_p ( M^{-1/2} v_s ) ].

    dihedral_A: optional dict mapping parameter index -> full (3N,3N) sensitivity
    matrix for dihedral terms computed in Python.
    """
    n3 = V.shape[0]
    n_modes = V.shape[1]
    np_free = fitter.n_free
    inv_sqrt_m = 1.0 / np.sqrt(np.repeat(masses, 3))
    X = inv_sqrt_m[:, None] * V  # Cartesian coordinates = M^{-1/2} v_s
    actions = np.zeros((np_free, n3, n_modes))
    for p in range(np_free):
        k = np.zeros(np_free)
        k[p] = 1.0
        fitter.set_params(k)
        A_p = fitter.compute_model_hessian()  # (n3, n3)
        if dihedral_A and p in dihedral_A:
            A_p = A_p + dihedral_A[p]
        Y = A_p @ X  # (n3, n_modes)
        actions[p] = inv_sqrt_m[:, None] * Y
    return actions


def accumulate_mode_normal_equations(G, y, actions, V, lambdas,
                                     mode_weight='relative',
                                     lambda_floor=0.0,
                                     beta=0.0, H_ref=None, H0=None,
                                     reg=0.0, k0=None):
    """Accumulate one system's mode-basis fit into G, y (in-place).

    mode_weight: 'relative' -> w_s = 1/max(lambda_s, lambda_floor)^2;
                 'adaptive' -> w_s = 1/max(lambda_s, median(lambda))^2;
                 'equal'    -> w_s = 1;
                 'inverse'  -> w_s = 1/max(lambda_s, lambda_floor);
                 array      -> use provided weights.
    lambda_floor: floor for lambda used in 'relative' and 'inverse' weights.
    beta: optional local-Hessian weight (currently 0.0 means no local-H term).
    """
    np_free, n3, n_modes = actions.shape
    lambdas = np.asarray(lambdas, dtype=float)
    if mode_weight == 'relative':
        weights = 1.0 / np.maximum(lambdas, lambda_floor)**2
    elif mode_weight == 'adaptive':
        lam_floor = max(np.median(lambdas), lambda_floor)
        weights = 1.0 / np.maximum(lambdas, lam_floor)**2
    elif mode_weight == 'inverse':
        weights = 1.0 / np.maximum(lambdas, lambda_floor)
    elif mode_weight == 'equal':
        weights = np.ones(n_modes, dtype=float)
    else:
        weights = np.asarray(mode_weight, dtype=float).reshape(n_modes)
    B = V * lambdas[None, :]  # target b_s = lambda_s v_s  (n3, n_modes)
    # Scale each mode by sqrt(weight) for normal equations
    W = np.sqrt(weights)
    A_w = actions * W[None, None, :]
    B_w = B * W[None, :]
    A_mat = A_w.reshape(np_free, -1)
    B_vec = B_w.reshape(-1)
    G += A_mat @ A_mat.T
    y += A_mat @ B_vec
    if beta > 0.0 and H_ref is not None:
        # Optional local-Hessian objective: beta * ||H_model - H_ref||^2
        # This requires the same parameter-sensitivity matrices; we can add
        # later by reusing build_sensitivity_matrices if needed.
        pass
    if reg > 0.0:
        G += reg * np.eye(np_free)
        if k0 is not None:
            y += reg * np.asarray(k0)


def fit_modes_multi(fitters, systems, n_free, freq_floor_cm1=10.0,
                    mode_weight='relative', reg=0.0, k0=None,
                    use_nnls=False, verbose=True,
                    dihedral_A_per_system=None):
    """Multi-system force-field fit in the mass-weighted mode basis.

    For each system, compute the reference modes/eigenvalues from H_ref,
    compute the action of each parameter's sensitivity matrix on those modes,
    and solve the linear least-squares problem
        sum_systems sum_s w_s | D(k) v_s - lambda_s v_s |^2 -> min.
    """
    dihedral_A_per_system = dihedral_A_per_system or [None] * len(systems)
    lambda_floor = (freq_floor_cm1 / 521.5)**2
    G = np.zeros((n_free, n_free))
    y = np.zeros(n_free)
    for f, sys, A_d in zip(fitters, systems, dihedral_A_per_system):
        V, lam, freqs = get_reference_modes_and_freqs(sys['H_ref'], sys['data']['masses'], sys['data'], freq_floor_cm1)
        if V.shape[1] == 0:
            if verbose: print(f"  WARNING: {sys['name']} has no vibrational modes above floor")
            continue
        actions = get_sensitivity_action(f, V, sys['data']['masses'], dihedral_A=A_d)
        accumulate_mode_normal_equations(G, y, actions, V, lam,
                                         mode_weight=mode_weight, lambda_floor=lambda_floor,
                                         reg=0.0, k0=None)  # reg applied at the end
    if reg > 0.0:
        G += reg * np.eye(n_free)
        if k0 is not None:
            y += reg * np.asarray(k0)
    if use_nnls:
        try:
            from scipy.optimize import nnls
            k, _ = nnls(G, y)
            return k
        except ImportError:
            if verbose: print("  WARNING: scipy not available, using unconstrained LSQ")
    try:
        return np.linalg.solve(G, y)
    except np.linalg.LinAlgError:
        if verbose: print("  WARNING: G is singular, using least-squares")
        return np.linalg.lstsq(G, y, rcond=None)[0]


def fit_gradient_descent_multi_cpp(fitters, systems, n_free, lr=1e-4, momentum=0.9, max_iter=5000, tol=1e-10, verbose=True):
    """Multi-system gradient descent using C++ FFfit library."""
    k = np.ones(n_free)
    velocity = np.zeros(n_free)
    prev_loss = 1e30
    for it in range(max_iter):
        grad = np.zeros(n_free)
        total_loss = [0.0]
        for f, sys in zip(fitters, systems):
            f.compute_gradient_loss(sys['H_ref_w'].ravel(), k, weight=sys.get('weight'),
                                    grad_out=grad, loss_out=total_loss)
        grad_norm = np.linalg.norm(grad)
        if verbose and (it % 100 == 0 or it < 10):
            print(f"  GD iter {it:4d}: loss={total_loss[0]:.6e} grad_norm={grad_norm:.6e}")
        if grad_norm < tol or (it > 0 and abs(prev_loss - total_loss[0]) < tol * prev_loss):
            if verbose: print(f"  GD converged at iter {it}: loss={total_loss[0]:.6e}")
            break
        prev_loss = total_loss[0]
        velocity = momentum * velocity - lr * grad
        k += velocity
    return k

def main():
    parser = argparse.ArgumentParser(description='Fit FF parameters to PySCF Hessians via C++ FFfit library')
    parser.add_argument('--cases', default='SiH4', help='Comma-separated case names or "all_Si" or "all"')
    parser.add_argument('--mass-weight', action='store_true', help='Use mass-weighted Hessian')
    parser.add_argument('--bond-cutoff', type=float, default=2.5, help='Bond distance cutoff (Angstrom)')
    parser.add_argument('--third-bonds', action='store_true', help='Add 3rd-neighbor (1-4) distance springs to mimic dihedral/cross-terms')
    parser.add_argument('--third-bond-cutoff', type=float, default=4.0, help='Distance cutoff for 3rd-neighbor springs (Angstrom)')
    parser.add_argument('--dihedrals', action='store_true', help='Add UFF/Prokop torsion (dihedral) terms to the fit')
    parser.add_argument('--dihedral-n', type=int, default=3, help='Dihedral periodicity n (default 3 for sp3)')
    parser.add_argument('--dihedral-d', type=float, default=1.0, help='Dihedral phase sign d (default 1.0)')
    parser.add_argument('--dihedral-h', type=float, default=1e-5, help='Finite-difference step for dihedral Hessian (Angstrom)')
    parser.add_argument('--mode-weight', type=str, default='equal', help='Mode weight: equal (default), relative (1/lambda^2), adaptive (per-system median floor), inverse (1/lambda), or array string')
    parser.add_argument('--reg', type=float, default=0.0, help='Ridge regularization weight')
    parser.add_argument('--use-nnls', action='store_true', help='Use non-negative least squares for mode fit')
    parser.add_argument('--plot', action='store_true', help='Plot reference vs model vibrational spectra')
    parser.add_argument('--plot-dir', type=str, default='tests/tSiNCs/OUT_FFfit_plots', help='Output directory for spectrum plots')
    parser.add_argument('--exp-spectrum', type=str, default=None, help='Optional experimental IR/FTIR spectrum file (two columns: frequency cm^-1, intensity)')
    parser.add_argument('--spectrum-width', type=float, default=20.0, help='Gaussian broadening width (cm^-1)')
    parser.add_argument('--spectrum-xmax', type=float, default=2500.0, help='Maximum frequency in spectrum plots (cm^-1)')
    args = parser.parse_args()
    
    # Determine case list
    if args.cases == 'all_Si':
        case_names = ['SiH4', 'Si_R3p8', 'Si_R4p5', 'Si_R5p0', 'Si_R5p5', 'Si_R6p0']
    elif args.cases == 'all':
        case_names = ['SiH4', 'Si_R3p8', 'Si_R4p5', 'Si_R5p0', 'Si_R5p5', 'Si_R6p0',
                      'C_diamond_R2p8', 'C_diamond_R3p2', 'C_diamond_R3p8', 'C_diamond_R4p2']
    else:
        case_names = [c.strip() for c in args.cases.split(',')]
    
    # === Load all systems ===
    print(f"=== Loading {len(case_names)} systems: {case_names} ===")
    systems = []
    for name in case_names:
        case_dir = os.path.join(RESULTS_DIR, name)
        if not os.path.isdir(case_dir):
            print(f"  WARNING: skipping {name} (directory not found)")
            continue
        data = load_pyscf_case(case_dir)
        natoms = data['natoms']
        H_ref = hessian_ha_bohr_to_ev_ang2(
            data['hessian'].transpose(0, 2, 1, 3).reshape(natoms*3, natoms*3))
        bonds, angles, bonds3 = build_topology(data['symbols'], data['positions'], args.bond_cutoff,
                                                third_bonds=args.third_bonds, third_bond_cutoff=args.third_bond_cutoff)
        dihedrals = build_dihedrals(data['symbols'], data['positions'], bonds,
                                    d=args.dihedral_d, n=args.dihedral_n, dihedral=args.dihedrals)
        if args.mass_weight:
            masses = data['masses']
            inv_sqrt_m = 1.0 / np.sqrt(np.repeat(masses, 3))
            weight = (inv_sqrt_m[:, None] * inv_sqrt_m[None, :]).flatten()
            H_ref_w = (inv_sqrt_m[:, None] * inv_sqrt_m[None, :]) * H_ref
        else:
            weight = None
            H_ref_w = H_ref
        sys = {
            'name': name, 'data': data, 'natoms': natoms, 'H_ref': H_ref, 'H_ref_w': H_ref_w,
            'bonds': bonds, 'angles': angles, 'bonds3': bonds3, 'dihedrals': dihedrals, 'weight': weight,
            'symbols': data['symbols'], 'positions': data['positions'],
        }
        systems.append(sys)
        print(f"  {name}: natoms={natoms}, bonds={len(bonds)}, angles={len(angles)}, 3rd_bonds={len(bonds3)}, dihedrals={len(dihedrals)}, "
              f"H_ref [{H_ref.min():.2e}, {H_ref.max():.2e}]")
    
    if not systems:
        print("ERROR: no valid systems found"); sys.exit(1)
    
    # === Build global parameter mapping ===
    all_bonds = [s['bonds'] for s in systems]
    all_angles = [s['angles'] for s in systems]
    all_bonds3 = [s['bonds3'] for s in systems]
    all_dihedrals = [s['dihedrals'] for s in systems]
    all_symbols = [s['symbols'] for s in systems]
    all_positions = [s['positions'] for s in systems]
    global_bond_map, global_bond3_map, global_angle_map, global_dihedral_map, n_total = build_global_param_map(
        all_bonds, all_angles, all_symbols, all_bonds3=all_bonds3, all_dihedrals=all_dihedrals)
    n_cpp = len(global_bond_map) + len(global_bond3_map) + len(global_angle_map)
    print(f"\n=== Global Parameter Mapping ===")
    print(f"  bond types ({len(global_bond_map)}): {global_bond_map}")
    if global_bond3_map:
        print(f"  3rd-neighbor bond types ({len(global_bond3_map)}): {global_bond3_map}")
    print(f"  angle types ({len(global_angle_map)}): {global_angle_map}")
    if global_dihedral_map:
        print(f"  dihedral types ({len(global_dihedral_map)}): {global_dihedral_map}")
    print(f"  total free params: {n_total} (C++ params: {n_cpp}, dihedral params: {len(global_dihedral_map)})")
    
    # === Compute averaged equilibrium l0/theta0 across all systems ===
    print(f"\n=== Averaged Equilibrium Parameters (multi-system) ===")
    avg_l0, avg_l0_3, avg_theta0 = compute_averaged_equilibrium(
        all_bonds, all_angles, all_symbols, all_positions,
        global_bond_map, global_angle_map, global_bond3_map=global_bond3_map, all_bonds3=all_bonds3)
    
    # Rebuild each system's bonds/angles with averaged equilibrium parameters
    for sys in systems:
        sys['bonds'], sys['angles'], sys['bonds3'] = rebuild_topology_with_averaged(
            sys['bonds'], sys['angles'], sys['bonds3'], sys['symbols'], sys['positions'],
            avg_l0, avg_theta0, avg_l0_3=avg_l0_3 if args.third_bonds else None)

    # === Precompute per-system dihedral sensitivity matrices (Python side) ===
    if args.dihedrals:
        print(f"\n=== Precomputing Dihedral Sensitivities ===")
        for sys in systems:
            sys['dihedral_A'] = compute_dihedral_sensitivity(
                sys['positions'], sys['symbols'], sys['dihedrals'], global_dihedral_map, h=args.dihedral_h)
            print(f"  {sys['name']}: computed {len(sys['dihedral_A'])} dihedral sensitivity matrices")
    else:
        for sys in systems:
            sys['dihedral_A'] = {}

    # === Create C++ FFfit instances with global param indices ===
    fitters = []
    for sys in systems:
        f = FFfit_cpp.FFfit()
        f.set_geometry(sys['positions'])
        f.set_symbols(sys['symbols'])
        for (i, j, l0) in sys['bonds']:
            f.add_bond(i, j, l0)
        for (i, j, l0) in sys['bonds3']:
            f.add_bond(i, j, l0)
        for (i, j, k, t0) in sys['angles']:
            f.add_angle(i, j, k, t0)
        n_bonds = len(sys['bonds'])
        for ib, (i, j, l0) in enumerate(sys['bonds'] + sys['bonds3']):
            key = tuple(sorted([sys['symbols'][i], sys['symbols'][j]]))
            if ib < n_bonds:
                f.set_bond_param(ib, global_bond_map[key])
            else:
                f.set_bond_param(ib, global_bond3_map[key])
        for ia, (i, j, k, t0) in enumerate(sys['angles']):
            key = angle_type_key(sys['symbols'], i, j, k)
            f.set_angle_param(ia, global_angle_map[key])
        f.set_n_free(n_total)
        fitters.append(f)
    
    # C++ fitting knows only bond/angle/3rd-bond parameters; switch to n_cpp for Methods 1 and 2
    for f in fitters: f.set_n_free(n_cpp)

    # === Method 1: Multi-system linear least-squares (C++) ===
    print(f"\n=== Method 1: Multi-System Linear Least-Squares (C++) ===")
    G = np.zeros((n_cpp, n_cpp))
    y = np.zeros(n_cpp)
    for f, sys in zip(fitters, systems):
        f.accumulate_normal_equations(G, y, sys['H_ref_w'].ravel(), weight=sys.get('weight'))
    k_lsq = FFfit_cpp.FFfit.solve_normal_equations(G, y)
    for key, t in global_bond_map.items():
        print(f"  k_bond[{key}] = {k_lsq[t]:.6e} eV/A^2")
    for key, t in global_bond3_map.items():
        print(f"  k_bond3[{key}] = {k_lsq[t]:.6e} eV/A^2")
    for key, t in global_angle_map.items():
        print(f"  k_angle[{key}] = {k_lsq[t]:.6e} eV/rad^2")

    # === Method 2: Multi-system gradient descent (C++) ===
    print(f"\n=== Method 2: Multi-System Gradient Descent (C++) ===")
    k_gd = fit_gradient_descent_multi_cpp(fitters, systems, n_cpp, lr=1e-4, momentum=0.9, max_iter=5000, tol=1e-10, verbose=True)
    for key, t in global_bond_map.items():
        print(f"  k_bond[{key}] = {k_gd[t]:.6e} eV/A^2")
    for key, t in global_bond3_map.items():
        print(f"  k_bond3[{key}] = {k_gd[t]:.6e} eV/A^2")
    for key, t in global_angle_map.items():
        print(f"  k_angle[{key}] = {k_gd[t]:.6e} eV/rad^2")

    print(f"\n=== Method Comparison ===")
    print(f"  LSQ vs GD max |diff|: {np.max(np.abs(k_lsq - k_gd)):.6e}")
    
    # === Per-system Hessian comparison (C++) ===
    print(f"\n=== Per-System Hessian Comparison ===")
    k = k_lsq
    for f, sys in zip(fitters, systems):
        f.set_params(k)
        H_model = f.compute_model_hessian()
        H_diff = H_model - sys['H_ref']
        rmsd = np.sqrt(np.mean(H_diff**2))
        nref = np.linalg.norm(sys['H_ref'])
        ndiff = np.linalg.norm(H_diff)
        nmodel = np.linalg.norm(H_model)
        rel_frob = ndiff / nref * 100
        print(f"  {sys['name']:12s}: RMSD={rmsd:.4e} eV/A^2  relFrob={rel_frob:.2f}%  "
              f"||H_ref||={nref:.4e}  ||H_model||={nmodel:.4e}  ||diff||={ndiff:.4e}")
    
    # === Method 3: Mode-basis least-squares (Python) ===
    print(f"\n=== Method 3: Mode-Basis Fit (Python) ===")
    if args.mode_weight in ('relative', 'adaptive', 'equal', 'inverse'):
        mode_weight = args.mode_weight
    else:
        # parse comma-separated list of floats or JSON array
        mode_weight = args.mode_weight
    for f in fitters: f.set_n_free(n_total)

    k_mode = fit_modes_multi(fitters, systems, n_total,
                             freq_floor_cm1=10.0,
                             mode_weight=mode_weight,
                             reg=args.reg,
                             use_nnls=args.use_nnls,
                             verbose=True,
                             dihedral_A_per_system=[s.get('dihedral_A', {}) for s in systems])
    for key, t in global_bond_map.items():
        print(f"  k_bond[{key}] = {k_mode[t]:.6e} eV/A^2")
    for key, t in global_bond3_map.items():
        print(f"  k_bond3[{key}] = {k_mode[t]:.6e} eV/A^2")
    for key, t in global_angle_map.items():
        print(f"  k_angle[{key}] = {k_mode[t]:.6e} eV/rad^2")
    for key, t in global_dihedral_map.items():
        print(f"  k_dihedral[{key}] = {k_mode[t]:.6e} eV")

    # === Mode-basis fit diagnostics ===
    print(f"\n=== Mode-Basis Fit Diagnostics ===")
    k = k_mode
    for f, sys in zip(fitters, systems):
        f.set_params(k)
        H_model = f.compute_model_hessian()
        add_dihedral_hessian(H_model, k, sys.get('dihedral_A', {}))
        H_diff = H_model - sys['H_ref']
        rmsd = np.sqrt(np.mean(H_diff**2))
        nref = np.linalg.norm(sys['H_ref'])
        ndiff = np.linalg.norm(H_diff)
        nmodel = np.linalg.norm(H_model)
        rel_frob = ndiff / nref * 100
        print(f"  {sys['name']:12s}: RMSD={rmsd:.4e} eV/A^2  relFrob={rel_frob:.2f}%  "
              f"||H_ref||={nref:.4e}  ||H_model||={nmodel:.4e}  ||diff||={ndiff:.4e}")
    
    # === Per-system frequency comparison ===
    print(f"\n=== Per-System Frequency Comparison (Mode-Basis Fit) ===")
    system_spectra_data = []
    for f, sys in zip(fitters, systems):
        f.set_params(k_mode)
        H_model = f.compute_model_hessian()
        add_dihedral_hessian(H_model, k_mode, sys.get('dihedral_A', {}))
        nz_ref, nz_model = compare_frequencies(sys['H_ref'], H_model, sys['data']['masses'], data=sys['data'], label=sys['name'], positions=sys['positions'])
        if args.plot:
            plot_spectrum(nz_ref, nz_model, sys['name'], outdir=args.plot_dir,
                          exp_spectrum=args.exp_spectrum, width=args.spectrum_width, xmax=args.spectrum_xmax)
        system_spectra_data.append((sys['name'], nz_ref, nz_model))

    if args.plot:
        plot_spectra_overlay(system_spectra_data, outdir=args.plot_dir, width=args.spectrum_width, xmax=args.spectrum_xmax)

    print(f"\n=== Done ===")

if __name__ == '__main__':
    main()
