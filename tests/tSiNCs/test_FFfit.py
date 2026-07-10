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

def angle_type_key(symbols, i, j, k, elements=None, central_only=False):
    outer = elements if central_only else symbols
    si, sk = outer[i], outer[k]
    return (si, symbols[j], sk) if si <= sk else (sk, symbols[j], si)


def assign_si_environment_types(symbols, bonds, enabled=False):
    """Classify Si by bonded-H coordination: SiH3 (>=3), SiH2, SiH, or bulk Si."""
    if not enabled:
        return list(symbols)
    nH = np.zeros(len(symbols), dtype=int)
    for i, j, _ in bonds:
        if symbols[i] == 'Si' and symbols[j] == 'H': nH[i] += 1
        if symbols[j] == 'Si' and symbols[i] == 'H': nH[j] += 1
    labels = list(symbols)
    for i, symbol in enumerate(symbols):
        if symbol != 'Si': continue
        if nH[i] > 4:
            raise ValueError(f"Si atom {i} has impossible H coordination {nH[i]}")
        labels[i] = 'SiH3' if nH[i] >= 3 else ('SiH2' if nH[i] == 2 else ('SiH' if nH[i] == 1 else 'Si'))
    return labels


def interaction_type_counts(systems):
    """Count observations supporting each environment-resolved interaction type."""
    counts = {'bond': {}, 'angle': {}, '1-4': {}, 'torsion': {}}
    for sys in systems:
        labels = sys.get('atom_types', sys['symbols'])
        for i, j, _ in sys['bonds']:
            key = tuple(sorted((labels[i], labels[j])))
            counts['bond'][key] = counts['bond'].get(key, 0) + 1
        for i, j, k, _ in sys['angles']:
            key = angle_type_key(labels, i, j, k, elements=sys['symbols'], central_only=sys.get('angle_central_only', False))
            counts['angle'][key] = counts['angle'].get(key, 0) + 1
        for i, j, _ in sys.get('bonds3', []):
            key = tuple(sorted((labels[i], labels[j])))
            counts['1-4'][key] = counts['1-4'].get(key, 0) + 1
        for i, j, k, l, _, _ in sys.get('dihedrals', []):
            key = dihedral_type_key(labels, i, j, k, l)
            counts['torsion'][key] = counts['torsion'].get(key, 0) + 1
    return counts


def print_interaction_type_counts(systems, min_count=3):
    counts = interaction_type_counts(systems)
    for kind, table in counts.items():
        if not table: continue
        print(f"  {kind} type support:")
        for key, count in sorted(table.items(), key=lambda item: str(item[0])):
            marker = '  LOW SUPPORT' if count < min_count else ''
            print(f"    {key}: {count}{marker}")
    return counts

def load_pyscf_case(case_dir, geometry_only=False):
    """Load PySCF results for a case.

    If geometry_only=True, skip loading hessian.npy, modes.npy, and frequencies
    (the large/expensive files) — only load relaxed geometry and masses.
    """
    d = {}
    if not geometry_only:
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
    l12 = np.linalg.norm(q12)
    l32 = np.linalg.norm(q32)
    l43 = np.linalg.norm(q43)
    if l12 < 1e-12 or l32 < 1e-12 or l43 < 1e-12:
        return 0.0, np.zeros_like(pos)
    u12 = q12 / l12
    u32 = q32 / l32
    u43 = q43 / l43

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
    # This is the signed UFF/Prokop convention used by Vec2d complex powers.
    phi = np.arctan2(csy, csx)
    nphi = n * phi
    csnx = np.cos(nphi)
    csny = np.sin(nphi)

    E = 1.0 + d * csnx
    # f scalar used to build forces (f = -V*d*n*csn.y  with V=1)
    f = -d * n * csny

    # hneigh.w in UFF.h stores inverse bond length, not bond length.
    il12, il32, il43 = 1.0/l12, 1.0/l32, 1.0/l43
    fp1 = -f * n123 * il2_123 * il12
    fp4 =  f * n234 * il2_234 * il43
    c123 = np.dot(u32, u12) * (il32 / il12)
    c432 = np.dot(u32, u43) * (il32 / il43)
    fp3 = -c123 * fp1 - (c432 + 1.0) * fp4
    fp2 = (c123 - 1.0) * fp1 + c432 * fp4

    grad = -np.array([fp1, fp2, fp3, fp4])
    return E, grad


def dihedral_angle(pos):
    """Signed UFF/Prokop dihedral angle in radians for atoms i-j-k-l."""
    p1, p2, p3, p4 = np.asarray(pos, dtype=float)
    u12 = p1 - p2; u12 /= np.linalg.norm(u12)
    u32 = p3 - p2; u32 /= np.linalg.norm(u32)
    u43 = p4 - p3; u43 /= np.linalg.norm(u43)
    n123 = np.cross(u12, u32)
    n234 = np.cross(u43, u32)
    den = np.linalg.norm(n123) * np.linalg.norm(n234)
    if den < 1e-14:
        raise ValueError("dihedral angle is singular for collinear bond vectors")
    return np.arctan2(-np.dot(n123, u43) / den, np.dot(n123, n234) / den)


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
    nparams = max(int(bond_types.max()) + 1 if len(bond_types) else 0,
                  int(angle_types.max()) + 1 if len(angle_types) else 0)
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
    # Mirrors C++ angle_dHdk_cos: rank-one g⊗g + prestress (c-c0)*C_cos, all scaled by 1/s0²
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
        dc = cos_theta - c0
        scale = 1.0 / s02
        vecs = {i: gi, j: gj, k: gk}
        for a_atom in [i, j, k]:
            for b_atom in [i, j, k]:
                a_idx = slice(a_atom*3, a_atom*3+3)
                b_idx = slice(b_atom*3, b_atom*3+3)
                A[t][a_idx, b_idx] += scale * np.outer(vecs[a_atom], vecs[b_atom])
        # Prestress term: (c-c0)*C_cos / s0² — matches C++ angle_dHdk_cos
        if abs(dc) > 1e-12:
            I3 = np.eye(3)
            uuT = np.outer(u, u); vvT = np.outer(v, v)
            vuT = np.outer(v, u); uvT = np.outer(u, v)
            C_ii = (-vuT - uvT + 3*cos_theta*uuT - cos_theta*I3) / (al*al)
            C_kk = (-uvT - vuT + 3*cos_theta*vvT - cos_theta*I3) / (cl*cl)
            C_ik = (I3 - vvT - uuT + cos_theta*uvT) / (al*cl)
            C_ki = C_ik.T
            C_ij = -C_ii - C_ik; C_ji = C_ij.T
            C_jk = -C_ik - C_kk; C_kj = C_jk.T
            C_jj = C_ii + C_ik + C_ki + C_kk
            Ccos = {i: {i: C_ii, j: C_ij, k: C_ik},
                    j: {i: C_ji, j: C_jj, k: C_jk},
                    k: {i: C_ki, j: C_kj, k: C_kk}}
            for a_atom in [i, j, k]:
                for b_atom in [i, j, k]:
                    a_idx = slice(a_atom*3, a_atom*3+3)
                    b_idx = slice(b_atom*3, b_atom*3+3)
                    A[t][a_idx, b_idx] += scale * dc * Ccos[a_atom][b_atom]
    
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

def build_global_param_map(all_bonds, all_angles, all_symbols, all_bonds3=None, all_dihedrals=None, all_elements=None, angle_central_only=False):
    """Build a global ParamMap from all systems' element types.
    
    all_bonds/all_angles: list of lists (one per system)
    all_symbols: list of symbol lists (one per system)
    all_bonds3: optional list of 3rd-neighbor (1-4) bond lists (one per system)
    all_dihedrals: optional list of dihedral tuples (one per system)
    Returns: (global_bond_type_map, global_bond3_type_map, global_angle_type_map, global_dihedral_map, n_free)
    """
    all_bonds3 = all_bonds3 or []
    all_dihedrals = all_dihedrals or []
    all_elements = all_symbols if all_elements is None else all_elements
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
    for angles, symbols, elements in zip(all_angles, all_symbols, all_elements):
        for (i, j, k, theta0) in angles:
            key = angle_type_key(symbols, i, j, k, elements=elements, central_only=angle_central_only)
            if key not in global_angle_map:
                global_angle_map[key] = n_bond_types + n_bond3_types + len(global_angle_map)
    # Pass 4: discover all dihedral types, offset after angle types
    n_angle_offset = n_bond_types + n_bond3_types
    n_dihedral_offset = n_angle_offset + len(global_angle_map)
    for dihedrals, symbols in zip(all_dihedrals, all_symbols):
        for (i, j, k, l, d, n) in dihedrals:
            key = dihedral_type_key(symbols, i, j, k, l)
            if key not in global_dihedral_map:
                global_dihedral_map[key] = n_dihedral_offset + len(global_dihedral_map)
    n_free = n_dihedral_offset + len(global_dihedral_map)
    indices = sorted(list(global_bond_map.values()) + list(global_bond3_map.values()) + list(global_angle_map.values()) + list(global_dihedral_map.values()))
    if indices != list(range(n_free)):
        raise RuntimeError(f"global parameter map is not contiguous: {indices}, expected 0..{n_free-1}")
    return global_bond_map, global_bond3_map, global_angle_map, global_dihedral_map, n_free

def compute_averaged_equilibrium(all_bonds, all_angles, all_symbols, all_positions,
                                 global_bond_map, global_angle_map, global_bond3_map=None, all_bonds3=None,
                                 all_elements=None, angle_central_only=False):
    """Compute averaged equilibrium bond lengths l0 and angle cosine c0 across all systems.

    The angle force is written in c=cos(theta), therefore the transferable equilibrium
    coordinate is c0=<cos(theta)>, not cos(<theta>). We store theta0=acos(c0) only
    because the C++ API currently stores AngleDef.theta0 in radians.
    """
    global_bond3_map = global_bond3_map or {}
    all_bonds3 = all_bonds3 or []
    all_elements = all_symbols if all_elements is None else all_elements
    bond_lengths = {}
    bond3_lengths = {}
    angle_cos_values = {}
    angle_theta_values = {}

    for bonds, angles, symbols, elements, positions in zip(all_bonds, all_angles, all_symbols, all_elements, all_positions):
        for (i, j, l0) in bonds:
            key = tuple(sorted([symbols[i], symbols[j]]))
            r = np.linalg.norm(positions[j] - positions[i])
            bond_lengths.setdefault(key, []).append(r)
        for (i, j, k, theta0) in angles:
            key = angle_type_key(symbols, i, j, k, elements=elements, central_only=angle_central_only)
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

def rebuild_topology_with_averaged(bonds, angles, bonds3, symbols, positions, avg_l0, avg_theta0, avg_l0_3=None, elements=None, angle_central_only=False):
    """Rebuild bond/angle lists with averaged equilibrium parameters.
    
    Returns new bonds/angles/bonds3 lists with l0 and theta0 replaced by averaged values.
    """
    new_bonds = []
    for (i, j, l0) in bonds:
        key = tuple(sorted([symbols[i], symbols[j]]))
        new_l0 = avg_l0[key]
        new_bonds.append((i, j, new_l0))

    new_angles = []
    elements = symbols if elements is None else elements
    for (i, j, k, theta0) in angles:
        key = angle_type_key(symbols, i, j, k, elements=elements, central_only=angle_central_only)
        new_theta0 = avg_theta0[key]
        new_angles.append((i, j, k, new_theta0))

    new_bonds3 = []
    if avg_l0_3 is not None:
        for (i, j, l0) in bonds3:
            key = tuple(sorted([symbols[i], symbols[j]]))
            new_l0 = avg_l0_3[key]
            new_bonds3.append((i, j, new_l0))

    return new_bonds, new_angles, new_bonds3


def equilibrium_distributions(systems, include_third=True, include_dihedrals=True):
    """Collect relaxed internal-coordinate values by transferable interaction type."""
    values = {'bond': {}, 'angle': {}, '1-4': {}, 'dihedral': {}}
    for sys in systems:
        symbols, positions = sys.get('atom_types', sys['symbols']), sys['positions']
        for i, j, _ in sys['bonds']:
            key = '-'.join(sorted((symbols[i], symbols[j])))
            values['bond'].setdefault(key, []).append(np.linalg.norm(positions[j] - positions[i]))
        for i, j, k, _ in sys['angles']:
            key = '-'.join(angle_type_key(symbols, i, j, k, elements=sys['symbols'], central_only=sys.get('angle_central_only', False)))
            u = positions[i] - positions[j]; u /= np.linalg.norm(u)
            v = positions[k] - positions[j]; v /= np.linalg.norm(v)
            values['angle'].setdefault(key, []).append(np.degrees(np.arccos(np.clip(np.dot(u, v), -1.0, 1.0))))
        if include_third:
            for i, j, _ in sys.get('bonds3', []):
                key = '-'.join(sorted((symbols[i], symbols[j])))
                values['1-4'].setdefault(key, []).append(np.linalg.norm(positions[j] - positions[i]))
        if include_dihedrals:
            for i, j, k, l, _, _ in sys.get('dihedrals', []):
                key = '-'.join((symbols[i], symbols[j], symbols[k], symbols[l]))
                values['dihedral'].setdefault(key, []).append(np.degrees(dihedral_angle(positions[[i, j, k, l]])))
    return values


_SI_SUBTYPE_PREFIXES = ('SiH3', 'SiH2', 'SiH', 'Si')

def _element_family_key(subtype_key):
    """Map a subtype key string (e.g. 'H-SiH2', 'SiH-SiH2', 'H-SiH-Si') to element-only family (e.g. 'H-Si', 'Si-Si', 'H-Si-Si')."""
    parts = subtype_key.split('-')
    elem_parts = []
    for p in parts:
        matched = False
        for pref in _SI_SUBTYPE_PREFIXES:
            if p == pref:
                elem_parts.append('Si')
                matched = True
                break
        if not matched:
            elem_parts.append(p)
    return '-'.join(elem_parts)


def _plot_stacked_by_family(values, kind, title, unit, outdir, filename):
    """Plot one subplot per element-family, with subtypes as stacked histograms."""
    raw = values[kind]
    if not raw:
        fig, ax = plt.subplots(1, 1, figsize=(6, 4), constrained_layout=True)
        ax.text(0.5, 0.5, 'not enabled', ha='center', va='center', transform=ax.transAxes)
        ax.set_title(title)
        fig.savefig(os.path.join(outdir, filename), dpi=170)
        plt.close(fig)
        return
    # Group subtype keys by element family
    families = {}
    for sk in sorted(raw.keys()):
        fk = _element_family_key(sk)
        families.setdefault(fk, []).append(sk)
    fam_keys = sorted(families.keys())
    ncol = min(3, len(fam_keys))
    nrow = int(np.ceil(len(fam_keys) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.5*ncol, 4.0*nrow), constrained_layout=True, squeeze=False)
    rows = []
    for fi, fk in enumerate(fam_keys):
        ax = axes[fi // ncol, fi % ncol]
        subtypes = families[fk]
        all_vals = [np.asarray(raw[sk]) for sk in subtypes]
        vmin = min(v.min() for v in all_vals); vmax = max(v.max() for v in all_vals)
        span = vmax - vmin
        pad = 0.04*span if span > 1e-10 else max(abs(vmin)*1e-6, 1e-6)
        bins = np.linspace(vmin - pad, vmax + pad, 31)
        plot_vals = []
        labels = []
        for si, sk in enumerate(subtypes):
            v = np.asarray(raw[sk])
            plot_vals.append(v)
            labels.append(f'{sk} (n={v.size})')
            rows.append((kind, sk, v.size, np.mean(v), np.std(v), np.min(v), np.max(v), unit))
        ax.hist(plot_vals, bins=bins, stacked=True, histtype='stepfilled', alpha=0.7, linewidth=0.8, label=labels)
        ax.set_title(fk, fontsize=11)
        ax.set_xlabel(unit)
        ax.set_ylabel('count')
        ax.ticklabel_format(axis='x', style='plain', useOffset=False)
        ax.grid(True, alpha=0.2, ls=':')
        ax.legend(fontsize=7, loc='best')
    for fi in range(len(fam_keys), nrow*ncol):
        axes[fi // ncol, fi % ncol].axis('off')
    fig.suptitle(title, fontsize=13)
    out = os.path.join(outdir, filename)
    fig.savefig(out, dpi=170)
    plt.close(fig)
    return rows, out


def plot_equilibrium_distributions(systems, outdir):
    """Plot relaxed bond, angle, 1-4, and signed torsion distributions.

    Each element-family (e.g. Si-Si, Si-H for bonds; H-Si-H, Si-Si-H, Si-Si-Si
    for angles) gets its own subplot with an appropriate x-axis range.
    Subtypes (e.g. SiH, SiH2, SiH3) are shown as color-coded stacked histograms
    within each family panel.
    """
    values = equilibrium_distributions(systems)
    os.makedirs(outdir, exist_ok=True)
    all_rows = []
    specs = [
        ('bond',     'Bond length distributions by element family',   'Angstrom', 'equilibrium_bond_distributions.png'),
        ('angle',    'Bond angle distributions by element family',    'degree',   'equilibrium_angle_distributions.png'),
        ('1-4',      '1-4 endpoint distance by element family',       'Angstrom', 'equilibrium_1-4_distributions.png'),
        ('dihedral', 'Signed torsion distributions by element family','degree',   'equilibrium_dihedral_distributions.png'),
    ]
    for kind, title, unit, filename in specs:
        result = _plot_stacked_by_family(values, kind, title, unit, outdir, filename)
        if result is not None:
            all_rows.extend(result[0])
            print(f"  Saved {kind} distributions: {result[1]}")
    table = os.path.join(outdir, 'equilibrium_parameter_distributions.csv')
    with open(table, 'w') as f:
        f.write('coordinate,type,count,mean,std,min,max,unit\n')
        for row in all_rows: f.write(','.join(map(str, row)) + '\n')
    print(f"  Saved equilibrium statistics: {table}")
    return values, all_rows, outdir


def dft_stiffness_distributions(systems):
    """Extract DFT stiffnesses via Wilson GF projection (F = C^{+T} D C^+).

    For each system, project the DFT Hessian into internal coordinates.
    The diagonal of F gives the directly-extracted stiffness for each
    bond (eV/Å²) and angle (eV/rad²) — no fitting needed.

    Returns values dict with 'bond' and 'angle' keys, mapping type->list of stiffnesses.
    """
    values = {'bond': {}, 'angle': {}}
    for sys in systems:
        symbols = sys.get('atom_types', sys['symbols'])
        elements = sys['symbols']
        positions = sys['positions']
        masses = sys['data']['masses']
        H_ref = sys['H_ref']
        bonds = sys['bonds']
        angles = sys['angles']
        angle_central_only = sys.get('angle_central_only', False)
        B, labels = FFfit_cpp.build_wilson_matrix(positions, bonds, angles)
        # Coordinate scale: bonds scaled by r0 (dimensionless Δr/r0), angles in radians (scale=1)
        coord_scale = np.ones(B.shape[0])
        for ib, (i, j, _) in enumerate(bonds):
            coord_scale[ib] = np.linalg.norm(positions[j] - positions[i])
        F, info = FFfit_cpp.internal_hessian_projection(H_ref, B, masses, coordinate_scale=coord_scale)
        # F is (n_int, n_int); diagonal = directly-extracted stiffnesses
        Fdiag = np.diag(F)
        for ib, (i, j, _) in enumerate(bonds):
            key = '-'.join(sorted((symbols[i], symbols[j])))
            # F is in dimensionless coords (Δr/r0), so F_diag has units eV (energy per dimensionless displacement²)
            # Convert back to eV/Å²: k = F_diag / r0²
            r0 = coord_scale[ib]
            values['bond'].setdefault(key, []).append(Fdiag[ib] / (r0 * r0))
        for ia, (i, j, k, _) in enumerate(angles):
            key = '-'.join(FFfit_cpp.angle_type_key(symbols, i, j, k, elements=elements, central_only=angle_central_only)) if hasattr(FFfit_cpp, 'angle_type_key') else '-'.join(angle_type_key(symbols, i, j, k, elements=elements, central_only=angle_central_only))
            # Angles are already in radians (scale=1), so F_diag has units eV/rad²
            values['angle'].setdefault(key, []).append(Fdiag[len(bonds) + ia])
    # Convert lists to arrays
    for kind in values:
        for key in values[kind]:
            values[kind][key] = np.asarray(values[kind][key])
    return values


def plot_dft_stiffness_distributions(systems, outdir):
    """Plot DFT-extracted bond and angle stiffness histograms by element family.

    Uses the Wilson GF method (internal_hessian_projection) to directly extract
    valence force constants from the DFT Hessian — no fitting required.
    Stiffnesses are grouped by element family (Si-Si, Si-H for bonds;
    H-Si-H, H-Si-Si, Si-Si-Si for angles) with subtypes stacked.
    """
    values = dft_stiffness_distributions(systems)
    os.makedirs(outdir, exist_ok=True)
    all_rows = []
    specs = [
        ('bond',  'DFT bond stiffness (Wilson GF diagonal)',  'eV/Å²',  'dft_bond_stiffness_distributions.png'),
        ('angle', 'DFT angle stiffness (Wilson GF diagonal)', 'eV/rad²', 'dft_angle_stiffness_distributions.png'),
    ]
    for kind, title, unit, filename in specs:
        result = _plot_stacked_by_family(values, kind, title, unit, outdir, filename)
        if result is not None:
            all_rows.extend(result[0])
            print(f"  Saved {kind} stiffness distributions: {result[1]}")
    table = os.path.join(outdir, 'dft_stiffness_distributions.csv')
    with open(table, 'w') as f:
        f.write('coordinate,type,count,mean,std,min,max,unit\n')
        for row in all_rows: f.write(','.join(map(str, row)) + '\n')
    print(f"  Saved DFT stiffness statistics: {table}")
    return values, all_rows, outdir


# ========== Stiffness visualization: clustering + interactive HTML ==========

def cluster_1d(values, gap_threshold=2.0):
    """Identify disjunct clusters in 1D data via gap-based splitting.

    Sorts values, finds gaps between consecutive points where the gap exceeds
    gap_threshold * median_gap. Returns list of (center, std, indices) tuples.
    """
    v = np.sort(np.asarray(values, dtype=float))
    if len(v) < 2:
        return [(float(v.mean()), float(v.std()) if len(v) > 1 else 0.0, np.arange(len(v)))]
    gaps = np.diff(v)
    median_gap = np.median(gaps)
    if median_gap < 1e-12:
        return [(float(v.mean()), float(v.std()), np.arange(len(v)))]
    is_gap = gaps > gap_threshold * median_gap
    split_idx = np.where(is_gap)[0] + 1
    clusters = []
    prev = 0
    for si in split_idx:
        idx = np.arange(prev, si)
        clusters.append((float(v[idx].mean()), float(v[idx].std()), idx))
        prev = si
    idx = np.arange(prev, len(v))
    clusters.append((float(v[idx].mean()), float(v[idx].std()), idx))
    return clusters


def prepare_stiffness_viz_data(systems):
    """Compute per-bond/angle stiffnesses via Wilson GF, assign clusters, and compute family vmin/vmax.

    Returns list of dicts, one per system, each with:
      - positions, symbols, atom_types
      - bonds: list of (i, j, r0, stiffness, subtype_key, family_key, cluster_id)
      - angles: list of (i, j, k, theta_deg, stiffness, subtype_key, family_key, cluster_id)
      - family_ranges: {family_key: (vmin, vmax)} for bonds and angles separately
    """
    # First pass: collect all stiffnesses per subtype for clustering
    all_bond_stiff = {}  # subtype_key -> list of (sys_idx, bond_idx, stiffness)
    all_angle_stiff = {}
    sys_data = []
    for si, sys in enumerate(systems):
        symbols = sys.get('atom_types', sys['symbols'])
        elements = sys['symbols']
        positions = sys['positions']
        masses = sys['data']['masses']
        H_ref = sys['H_ref']
        bonds = sys['bonds']
        angles = sys['angles']
        angle_central_only = sys.get('angle_central_only', False)
        B, labels = FFfit_cpp.build_wilson_matrix(positions, bonds, angles)
        coord_scale = np.ones(B.shape[0])
        for ib, (i, j, _) in enumerate(bonds):
            coord_scale[ib] = np.linalg.norm(positions[j] - positions[i])
        F, info = FFfit_cpp.internal_hessian_projection(H_ref, B, masses, coordinate_scale=coord_scale)
        Fdiag = np.diag(F)
        bond_records = []
        for ib, (i, j, _) in enumerate(bonds):
            r0 = coord_scale[ib]
            k_val = Fdiag[ib] / (r0 * r0)
            sk = '-'.join(sorted((symbols[i], symbols[j])))
            fk = _element_family_key(sk)
            bond_records.append({'i': i, 'j': j, 'r0': r0, 'stiff': k_val, 'subtype': sk, 'family': fk})
            all_bond_stiff.setdefault(sk, []).append((si, ib, k_val))
        angle_records = []
        for ia, (i, j, k, _) in enumerate(angles):
            k_val = Fdiag[len(bonds) + ia]
            u = positions[i] - positions[j]; u /= np.linalg.norm(u)
            v = positions[k] - positions[j]; v /= np.linalg.norm(v)
            theta = np.degrees(np.arccos(np.clip(np.dot(u, v), -1.0, 1.0)))
            ak = angle_type_key(symbols, i, j, k, elements=elements, central_only=angle_central_only)
            sk = '-'.join(ak)
            fk = _element_family_key(sk)
            angle_records.append({'i': i, 'j': j, 'k': k, 'theta': theta, 'stiff': k_val, 'subtype': sk, 'family': fk})
            all_angle_stiff.setdefault(sk, []).append((si, ia, k_val))
        sys_data.append({
            'name': sys['name'], 'positions': positions, 'symbols': sys['symbols'],
            'atom_types': symbols, 'bonds': bond_records, 'angles': angle_records,
        })
    # Cluster per subtype
    bond_clusters = {}  # (sys_idx, bond_idx) -> cluster_id
    for sk, items in all_bond_stiff.items():
        vals = np.array([x[2] for x in items])
        clusters = cluster_1d(vals)
        for ci, (center, std, indices) in enumerate(clusters):
            for idx in indices:
                si, bi, _ = items[idx]
                bond_clusters[(si, bi)] = ci
    angle_clusters = {}
    for sk, items in all_angle_stiff.items():
        vals = np.array([x[2] for x in items])
        clusters = cluster_1d(vals)
        for ci, (center, std, indices) in enumerate(clusters):
            for idx in indices:
                si, ai, _ = items[idx]
                angle_clusters[(si, ai)] = ci
    # Assign cluster ids and compute family ranges
    bond_family_vals = {}
    angle_family_vals = {}
    for si, sd in enumerate(sys_data):
        for bi, br in enumerate(sd['bonds']):
            br['cluster'] = bond_clusters.get((si, bi), 0)
            bond_family_vals.setdefault(br['family'], []).append(br['stiff'])
        for ai, ar in enumerate(sd['angles']):
            ar['cluster'] = angle_clusters.get((si, ai), 0)
            angle_family_vals.setdefault(ar['family'], []).append(ar['stiff'])
    family_ranges = {'bond': {}, 'angle': {}}
    for fk, vals in bond_family_vals.items():
        family_ranges['bond'][fk] = (float(np.min(vals)), float(np.max(vals)))
    for fk, vals in angle_family_vals.items():
        family_ranges['angle'][fk] = (float(np.min(vals)), float(np.max(vals)))
    for sd in sys_data:
        sd['family_ranges'] = family_ranges
    return sys_data


def _json_safe(obj):
    """Recursively convert numpy types to native Python for JSON serialization."""
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    if isinstance(obj, dict):
        return {k: _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    return obj


def build_stiffness_html(sys_data, outdir):
    """Generate a self-contained interactive p5.js HTML for a single nanocrystal.

    Features:
      - 3D ball-and-stick rendering with mouse rotation
      - Bonds colored by stiffness (rainbow, vmin/vmax from family)
      - Angles shown as arcs colored by stiffness
      - Dropdown to select: all bonds, specific family, specific cluster
      - Toggle bonds/angles visibility
    """
    name = sys_data['name']
    positions = sys_data['positions']
    symbols = sys_data['symbols']
    atom_types = sys_data['atom_types']
    bonds = sys_data['bonds']
    angles = sys_data['angles']
    family_ranges = sys_data['family_ranges']
    # Build compact JSON data
    data = {
        'name': name,
        'positions': positions.tolist(),
        'symbols': symbols,
        'atom_types': atom_types,
        'bonds': [{'i': b['i'], 'j': b['j'], 'r0': b['r0'], 'stiff': b['stiff'], 'subtype': b['subtype'], 'family': b['family'], 'cluster': b['cluster']} for b in bonds],
        'angles': [{'i': a['i'], 'j': a['j'], 'k': a['k'], 'theta': a['theta'], 'stiff': a['stiff'], 'subtype': a['subtype'], 'family': a['family'], 'cluster': a['cluster']} for a in angles],
        'family_ranges': family_ranges,
    }
    data_json = json.dumps(_json_safe(data))
    # Collect unique families and subtypes for dropdowns
    bond_families = sorted(set(b['family'] for b in bonds))
    angle_families = sorted(set(a['family'] for a in angles))
    bond_subtypes = sorted(set(b['subtype'] for b in bonds))
    angle_subtypes = sorted(set(a['subtype'] for a in angles))
    html = f'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Stiffness Map — {name}</title>
<script src="https://cdnjs.cloudflare.com/ajax/libs/p5.js/1.9.0/p5.min.js"></script>
<style>
body {{ margin: 0; padding: 0; background: #1a1a2e; color: #eee; font-family: sans-serif; overflow: hidden; }}
#controls {{ position: absolute; top: 10px; left: 10px; z-index: 10; background: rgba(20,20,40,0.85); padding: 12px; border-radius: 8px; max-width: 320px; }}
#controls h3 {{ margin: 0 0 8px 0; font-size: 14px; }}
#controls label {{ font-size: 12px; display: block; margin: 6px 0 2px 0; }}
#controls select, #controls button {{ font-size: 12px; width: 100%; margin: 2px 0; }}
#info {{ position: absolute; bottom: 10px; left: 10px; z-index: 10; background: rgba(20,20,40,0.85); padding: 8px; border-radius: 8px; font-size: 11px; max-width: 400px; }}
canvas {{ display: block; }}
.legend-bar {{ width: 200px; height: 12px; border-radius: 4px; margin: 4px 0; }}
#histPanel {{ position: absolute; top: 10px; right: 10px; z-index: 10; background: rgba(20,20,40,0.9); padding: 10px; border-radius: 8px; width: 340px; }}
#histPanel h3 {{ margin: 0 0 6px 0; font-size: 13px; }}
#histCanvas {{ background: #111; border-radius: 4px; }}
#histPanel label {{ font-size: 11px; display: inline-block; margin: 4px 2px 0 0; }}
#histPanel input[type=range] {{ width: 120px; vertical-align: middle; }}
</style>
</head>
<body>
<div id="controls">
  <h3>{name} — Stiffness Map</h3>
  <label>Show:</label>
  <select id="modeSelect">
    <option value="bonds_all">All bonds</option>
    <option value="angles_all">All angles</option>
  </select>
  <label>Filter by family:</label>
  <select id="familySelect"><option value="all">All families</option></select>
  <label>Filter by cluster:</label>
  <select id="clusterSelect"><option value="all">All clusters</option></select>
  <label>Bond color mode:</label>
  <select id="colorMode">
    <option value="stiffness">Stiffness (rainbow)</option>
    <option value="length">Bond length (rainbow)</option>
    <option value="cluster">Cluster (discrete)</option>
    <option value="subtype">Subtype (discrete)</option>
  </select>
  <label>Atom size:</label>
  <input type="range" id="atomSize" min="5" max="30" value="14" style="width:100%">
  <div id="legendContainer"></div>
</div>
<div id="info"></div>
<div id="histPanel">
  <h3>Stiffness Histogram</h3>
  <canvas id="histCanvas" width="320" height="180"></canvas>
  <div id="histControls">
    <label>vmin: <input type="range" id="vminSlider" min="0" max="100" value="0" step="0.1"></label>
    <label>vmax: <input type="range" id="vmaxSlider" min="0" max="100" value="100" step="0.1"></label>
    <span id="rangeLabel" style="font-size:11px;color:#8af;"></span>
  </div>
</div>
<script>
const DATA = {data_json};
let rotX = -0.3, rotY = 0.5, camZoom = 1.0;
let dragging = false, lastMX = 0, lastMY = 0;
let mode = 'bonds_all', familyFilter = 'all', clusterFilter = 'all', colorMode = 'stiffness';
let atomSize = 14;
let center = [0, 0, 0];
let scaleF = 1.0;
let vminUser = null, vmaxUser = null; // user override for viz range

function setup() {{
  createCanvas(windowWidth, windowHeight, WEBGL);
  // Compute center and scale
  let pos = DATA.positions;
  for (let i = 0; i < pos.length; i++) {{
    center[0] += pos[i][0]; center[1] += pos[i][1]; center[2] += pos[i][2];
  }}
  center = [center[0]/pos.length, center[1]/pos.length, center[2]/pos.length];
  let maxR = 0;
  for (let i = 0; i < pos.length; i++) {{
    let d = dist(pos[i][0], pos[i][1], pos[i][2], center[0], center[1], center[2]);
    if (d > maxR) maxR = d;
  }}
  scaleF = 250 / maxR;
  // Populate dropdowns
  let famSel = document.getElementById('familySelect');
  let bondFams = {json.dumps(bond_families)};
  let angFams = {json.dumps(angle_families)};
  for (let f of bondFams) {{ let o = document.createElement('option'); o.value = 'bond_' + f; o.text = 'Bond: ' + f; famSel.appendChild(o); }}
  for (let f of angFams) {{ let o = document.createElement('option'); o.value = 'angle_' + f; o.text = 'Angle: ' + f; famSel.appendChild(o); }}
  let modeSel = document.getElementById('modeSelect');
  modeSel.addEventListener('change', e => {{ mode = e.target.value; updateClusterDropdown(); }});
  famSel.addEventListener('change', e => {{ familyFilter = e.target.value; updateClusterDropdown(); }});
  document.getElementById('clusterSelect').addEventListener('change', e => {{ clusterFilter = e.target.value; updateHistSliders(); }});
  document.getElementById('colorMode').addEventListener('change', e => {{ colorMode = e.target.value; updateLegend(); }});
  document.getElementById('atomSize').addEventListener('input', e => {{ atomSize = parseFloat(e.target.value); }});
  // Histogram range sliders
  let vminSlider = document.getElementById('vminSlider');
  let vmaxSlider = document.getElementById('vmaxSlider');
  vminSlider.addEventListener('input', e => {{
    let v = parseFloat(e.target.value);
    let range = getStiffRange();
    vminUser = range[0] + (v/100) * (range[1] - range[0]);
    updateRangeLabel();
  }});
  vmaxSlider.addEventListener('input', e => {{
    let v = parseFloat(e.target.value);
    let range = getStiffRange();
    vmaxUser = range[0] + (v/100) * (range[1] - range[0]);
    updateRangeLabel();
  }});
  updateLegend();
  updateHistSliders();
}}

function getStiffRange() {{
  let isAngle = mode.startsWith('angle');
  let items = isAngle ? DATA.angles : DATA.bonds;
  if (familyFilter !== 'all') {{
    let prefix = isAngle ? 'angle_' : 'bond_';
    if (familyFilter.startsWith(prefix)) {{
      let fam = familyFilter.substring(prefix.length);
      items = items.filter(it => it.family === fam);
    }} else return [0, 1];
  }}
  if (items.length === 0) return [0, 1];
  let vals = items.map(it => it.stiff);
  return [Math.min(...vals), Math.max(...vals)];
}}

function updateRangeLabel() {{
  let lbl = document.getElementById('rangeLabel');
  let isAngle = mode.startsWith('angle');
  let unit = isAngle ? 'eV/rad²' : 'eV/Å²';
  if (vminUser !== null && vmaxUser !== null)
    lbl.textContent = ' [' + vminUser.toFixed(3) + ' — ' + vmaxUser.toFixed(3) + ' ' + unit + ']';
  else
    lbl.textContent = ' [auto]';
}}

function updateHistSliders() {{
  let r = getStiffRange();
  vminUser = r[0]; vmaxUser = r[1];
  document.getElementById('vminSlider').value = 0;
  document.getElementById('vmaxSlider').value = 100;
  updateRangeLabel();
}}

function getFamilyFilteredItems() {{
  let isAngle = mode.startsWith('angle');
  let items = isAngle ? DATA.angles : DATA.bonds;
  if (familyFilter !== 'all') {{
    let prefix = isAngle ? 'angle_' : 'bond_';
    if (familyFilter.startsWith(prefix)) {{
      let fam = familyFilter.substring(prefix.length);
      items = items.filter(it => it.family === fam);
    }} else return [];
  }}
  return items;
}}

function drawHistogram() {{
  let canvas = document.getElementById('histCanvas');
  let ctx = canvas.getContext('2d');
  let W = canvas.width, H = canvas.height;
  ctx.clearRect(0, 0, W, H);
  ctx.fillStyle = '#111';
  ctx.fillRect(0, 0, W, H);
  let isAngle = mode.startsWith('angle');
  let items = getFamilyFilteredItems();
  if (items.length === 0) {{ ctx.fillStyle = '#666'; ctx.font = '12px sans-serif'; ctx.fillText('No data', W/2-30, H/2); return; }}
  let vals = items.map(it => it.stiff);
  let vMin = Math.min(...vals), vMax = Math.max(...vals);
  let pad = (vMax - vMin) * 0.05;
  vMin -= pad; vMax += pad;
  let nBins = 40;
  let binW = (vMax - vMin) / nBins;
  let bins = new Array(nBins).fill(0);
  let binClusters = Array.from({{length: nBins}}, () => new Set());
  for (let it of items) {{
    let bi = Math.min(nBins - 1, Math.max(0, Math.floor((it.stiff - vMin) / binW)));
    bins[bi]++;
    binClusters[bi].add(it.cluster);
  }}
  let maxCount = Math.max(...bins);
  if (maxCount === 0) maxCount = 1;
  let plotW = W - 50, plotH = H - 30, ox = 40, oy = 10;
  // Axes
  ctx.strokeStyle = '#444'; ctx.lineWidth = 1;
  ctx.beginPath(); ctx.moveTo(ox, oy); ctx.lineTo(ox, oy + plotH); ctx.lineTo(ox + plotW, oy + plotH); ctx.stroke();
  // Bars colored by cluster (if cluster filter active, highlight that cluster)
  for (let i = 0; i < nBins; i++) {{
    if (bins[i] === 0) continue;
    let x = ox + (i / nBins) * plotW;
    let bw = plotW / nBins - 1;
    let bh = (bins[i] / maxCount) * plotH;
    let clusters = [...binClusters[i]].sort((a,b)=>a-b);
    if (clusterFilter !== 'all') {{
      let c = parseInt(clusterFilter);
      if (clusters.includes(c)) {{
        let col = clusterColor(c);
        ctx.fillStyle = 'rgb(' + col.join(',') + ')';
      }} else {{
        ctx.fillStyle = 'rgba(80,80,80,0.3)';
      }}
    }} else {{
      // Color by stiffness rainbow
      let midVal = vMin + (i + 0.5) * binW;
      let fr = getFamilyRange(items[0].family, isAngle);
      let col = rainbow(midVal, fr[0], fr[1]);
      ctx.fillStyle = 'rgb(' + Math.round(col[0]) + ',' + Math.round(col[1]) + ',' + Math.round(col[2]) + ')';
    }}
    ctx.fillRect(x, oy + plotH - bh, bw, bh);
  }}
  // vmin/vmax markers
  if (vminUser !== null && vmaxUser !== null) {{
    let xvmin = ox + ((vminUser - vMin) / (vMax - vMin)) * plotW;
    let xvmax = ox + ((vmaxUser - vMin) / (vMax - vMin)) * plotW;
    ctx.strokeStyle = '#ff0'; ctx.lineWidth = 1.5;
    ctx.beginPath(); ctx.moveTo(xvmin, oy); ctx.lineTo(xvmin, oy + plotH); ctx.stroke();
    ctx.beginPath(); ctx.moveTo(xvmax, oy); ctx.lineTo(xvmax, oy + plotH); ctx.stroke();
    // Highlight range between
    ctx.fillStyle = 'rgba(255,255,0,0.08)';
    ctx.fillRect(xvmin, oy, xvmax - xvmin, plotH);
  }}
  // Axis labels
  ctx.fillStyle = '#aaa'; ctx.font = '10px sans-serif';
  ctx.fillText(vMin.toFixed(2), ox, H - 5);
  ctx.fillText(vMax.toFixed(2), ox + plotW - 30, H - 5);
  let unit = isAngle ? 'eV/rad²' : 'eV/Å²';
  ctx.fillText(unit, W - 50, H - 5);
  ctx.fillText('count', 5, oy + 10);
}}

function getFamilyRange(fk, isAngle) {{
  let ranges = isAngle ? DATA.family_ranges.angle : DATA.family_ranges.bond;
  return ranges[fk] ? ranges[fk] : [0, 1];
}}

function mouseWheel(event) {{
  camZoom *= event.delta > 0 ? 1.1 : 0.9;
  camZoom = Math.max(0.1, Math.min(10, camZoom));
  return false;
}}

function mousePressed() {{
  if (mouseButton === LEFT) {{ dragging = true; lastMX = mouseX; lastMY = mouseY; }}
}}
function mouseReleased() {{ dragging = false; }}
function mouseDragged() {{
  if (dragging) {{
    rotY += (mouseX - lastMX) * 0.01;
    rotX += (mouseY - lastMY) * 0.01;
    rotX = Math.max(-Math.PI/2, Math.min(Math.PI/2, rotX));
    lastMX = mouseX; lastMY = mouseY;
  }}
}}

function updateClusterDropdown() {{
  let sel = document.getElementById('clusterSelect');
  sel.innerHTML = '<option value="all">All clusters</option>';
  let items = getFamilyFilteredItems();
  let clusters = [...new Set(items.map(it => it.cluster))].sort((a,b)=>a-b);
  for (let c of clusters) {{
    let count = items.filter(it => it.cluster === c).length;
    let o = document.createElement('option'); o.value = c; o.text = 'Cluster ' + c + ' (' + count + ')'; sel.appendChild(o);
  }}
  clusterFilter = 'all';
  updateHistSliders();
}}

function getFilteredItems() {{
  let isAngle = mode.startsWith('angle');
  let items = isAngle ? DATA.angles : DATA.bonds;
  if (familyFilter !== 'all') {{
    let prefix = isAngle ? 'angle_' : 'bond_';
    if (familyFilter.startsWith(prefix)) {{
      let fam = familyFilter.substring(prefix.length);
      items = items.filter(it => it.family === fam);
    }} else {{
      return [];
    }}
  }}
  if (clusterFilter !== 'all') {{
    let c = parseInt(clusterFilter);
    items = items.filter(it => it.cluster === c);
  }}
  return items;
}}

function rainbow(t, vmin, vmax) {{
  let u = (t - vmin) / (vmax - vmin + 1e-12);
  u = Math.max(0, Math.min(1, u));
  let r, g, b;
  if (u < 0.25) {{ r = 0; g = 4*u; b = 1; }}
  else if (u < 0.5) {{ r = 0; g = 1; b = 1 - 4*(u-0.25); }}
  else if (u < 0.75) {{ r = 4*(u-0.5); g = 1; b = 0; }}
  else {{ r = 1; g = 1 - 4*(u-0.75); b = 0; }}
  return [r*255, g*255, b*255];
}}

function clusterColor(c) {{
  let colors = [[255,100,100],[100,255,100],[100,100,255],[255,255,100],[255,100,255],[100,255,255],[200,200,100],[200,100,200]];
  let col = colors[c % colors.length];
  return col;
}}

function subtypeColor(sk, allSubtypes) {{
  let idx = allSubtypes.indexOf(sk);
  let colors = [[255,100,100],[100,255,100],[100,100,255],[255,255,100],[255,100,255],[100,255,255],[200,200,100],[200,100,200],[150,200,50],[50,200,150]];
  return colors[idx % colors.length];
}}

function updateLegend() {{
  let container = document.getElementById('legendContainer');
  let isAngle = mode.startsWith('angle');
  let items = getFilteredItems();
  if (items.length === 0) {{ container.innerHTML = ''; return; }}
  if (colorMode === 'stiffness' || colorMode === 'length') {{
    let ranges = isAngle ? DATA.family_ranges.angle : DATA.family_ranges.bond;
    let prefix = isAngle ? 'angle_' : 'bond_';
    let unit = colorMode === 'stiffness' ? (isAngle ? 'eV/rad²' : 'eV/Å²') : 'Å';
    let label = colorMode === 'stiffness' ? 'Stiffness' : 'Bond length';
    let fams = [...new Set(items.map(it => it.family))].sort();
    let html = '';
    for (let fk of fams) {{
      let vmin = ranges[fk] ? ranges[fk][0] : 0;
      let vmax = ranges[fk] ? ranges[fk][1] : 1;
      html += '<label>' + label + ' [' + fk + ']: ' + vmin.toFixed(3) + ' — ' + vmax.toFixed(3) + ' ' + unit + '</label>' +
        '<div class="legend-bar" style="background: linear-gradient(to right, rgb(0,0,255), rgb(0,255,255), rgb(0,255,0), rgb(255,255,0), rgb(255,0,0));"></div>';
    }}
    container.innerHTML = html;
  }} else if (colorMode === 'cluster') {{
    let clusters = [...new Set(items.map(it => it.cluster))].sort((a,b)=>a-b);
    let html = '<label>Clusters:</label>';
    for (let c of clusters) {{
      let col = clusterColor(c);
      html += '<div style="display:inline-block;margin:2px 6px;"><span style="display:inline-block;width:10px;height:10px;background:rgb(' + col.join(',') + ');border-radius:50%;"></span> C' + c + '</div>';
    }}
    container.innerHTML = html;
  }} else {{
    let subs = isAngle ? {json.dumps(angle_subtypes)} : {json.dumps(bond_subtypes)};
    let html = '<label>Subtypes:</label>';
    for (let s of subs) {{
      let col = subtypeColor(s, subs);
      html += '<div style="display:inline-block;margin:2px 6px;"><span style="display:inline-block;width:10px;height:10px;background:rgb(' + col.join(',') + ');border-radius:50%;"></span> ' + s + '</div>';
    }}
    container.innerHTML = html;
  }}
}}

function draw() {{
  background(26, 26, 46);
  let hw = width / (2 * camZoom), hh = height / (2 * camZoom);
  ortho(-hw, hw, -hh, hh, -2000, 2000);
  rotateX(rotX);
  rotateY(rotY);
  // Lighting
  ambientLight(80);
  directionalLight(200, 200, 200, 0.5, 0.5, -1);
  // Atoms
  let pos = DATA.positions;
  let atomColors = {{'Si': [100,180,255], 'H': [255,255,255], 'SiH': [120,200,255], 'SiH2': [140,220,255], 'SiH3': [160,240,255], 'C': [100,100,100], 'O': [255,100,100], 'N': [100,100,255]}};
  for (let i = 0; i < pos.length; i++) {{
    push();
    let at = DATA.atom_types[i];
    let col = atomColors[at] || atomColors[DATA.symbols[i]] || [200,200,200];
    fill(col[0], col[1], col[2]);
    noStroke();
    translate((pos[i][0]-center[0])*scaleF, (pos[i][1]-center[1])*scaleF, (pos[i][2]-center[2])*scaleF);
    sphere(atomSize * (DATA.symbols[i] === 'Si' ? 1.0 : 0.6));
    pop();
  }}
  // Always draw all bonds as thin gray background
  let isAngle = mode.startsWith('angle');
  if (!isAngle) {{
    for (let b of DATA.bonds) {{
      let p1 = pos[b.i], p2 = pos[b.j];
      push();
      stroke(60, 60, 70);
      strokeWeight(1.5);
      let x1 = (p1[0]-center[0])*scaleF, y1 = (p1[1]-center[1])*scaleF, z1 = (p1[2]-center[2])*scaleF;
      let x2 = (p2[0]-center[0])*scaleF, y2 = (p2[1]-center[1])*scaleF, z2 = (p2[2]-center[2])*scaleF;
      line(x1, y1, z1, x2, y2, z2);
      pop();
    }}
  }}
  // Highlighted items (filtered)
  let items = getFilteredItems();
  let ranges = isAngle ? DATA.family_ranges.angle : DATA.family_ranges.bond;
  let allSubtypes = isAngle ? {json.dumps(angle_subtypes)} : {json.dumps(bond_subtypes)};
  for (let it of items) {{
    let col;
    let inRange = true;
    if (vminUser !== null && vmaxUser !== null) {{
      inRange = (it.stiff >= vminUser && it.stiff <= vmaxUser);
    }}
    if (colorMode === 'stiffness') {{
      let fr = ranges[it.family];
      let vmin = fr ? fr[0] : 0, vmax = fr ? fr[1] : 1;
      col = rainbow(it.stiff, vmin, vmax);
    }} else if (colorMode === 'length' && !isAngle) {{
      let vmin = Math.min(...items.map(x => x.r0)); vmax = Math.max(...items.map(x => x.r0));
      col = rainbow(it.r0, vmin, vmax);
    }} else if (colorMode === 'cluster') {{
      col = clusterColor(it.cluster);
    }} else {{
      col = subtypeColor(it.subtype, allSubtypes);
    }}
    if (!isAngle) {{
      let p1 = pos[it.i], p2 = pos[it.j];
      push();
      if (inRange) {{
        stroke(col[0], col[1], col[2]);
        strokeWeight(5);
      }} else {{
        stroke(col[0]*0.3, col[1]*0.3, col[2]*0.3);
        strokeWeight(2);
      }}
      let x1 = (p1[0]-center[0])*scaleF, y1 = (p1[1]-center[1])*scaleF, z1 = (p1[2]-center[2])*scaleF;
      let x2 = (p2[0]-center[0])*scaleF, y2 = (p2[1]-center[1])*scaleF, z2 = (p2[2]-center[2])*scaleF;
      line(x1, y1, z1, x2, y2, z2);
      pop();
    }} else {{
      let p1 = pos[it.i], p2 = pos[it.j], p3 = pos[it.k];
      push();
      if (inRange) {{
        stroke(col[0], col[1], col[2]);
        strokeWeight(5);
      }} else {{
        stroke(col[0]*0.3, col[1]*0.3, col[2]*0.3);
        strokeWeight(2);
      }}
      let x2 = (p2[0]-center[0])*scaleF, y2 = (p2[1]-center[1])*scaleF, z2 = (p2[2]-center[2])*scaleF;
      let x1 = (p1[0]-center[0])*scaleF, y1 = (p1[1]-center[1])*scaleF, z1 = (p1[2]-center[2])*scaleF;
      let x3 = (p3[0]-center[0])*scaleF, y3 = (p3[1]-center[1])*scaleF, z3 = (p3[2]-center[2])*scaleF;
      line(x2, y2, z2, x1, y1, z1);
      line(x2, y2, z2, x3, y3, z3);
      pop();
    }}
  }}
  // Draw histogram overlay
  drawHistogram();
  // Update info
  let info = document.getElementById('info');
  let inCount = items.filter(it => (vminUser === null || (it.stiff >= vminUser && it.stiff <= vmaxUser))).length;
  info.innerHTML = 'System: ' + DATA.name + ' | Atoms: ' + pos.length + ' | Showing: ' + inCount + '/' + items.length + ' ' + (isAngle ? 'angles' : 'bonds') + ' | Drag to rotate, scroll to zoom';
}}

function windowResized() {{
  resizeCanvas(windowWidth, windowHeight);
}}
</script>
</body>
</html>'''
    outpath = os.path.join(outdir, f'stiffness_map_{name}.html')
    with open(outpath, 'w') as f:
        f.write(html)
    return outpath


def generate_stiffness_html(systems, outdir):
    """Generate interactive stiffness HTML maps for each nanocrystal.

    Steps:
      1. Compute per-bond/angle DFT stiffnesses via Wilson GF projection
      2. Identify disjunct clusters per subtype
      3. Generate one self-contained p5.js HTML per system
    """
    os.makedirs(outdir, exist_ok=True)
    sys_data_list = prepare_stiffness_viz_data(systems)
    paths = []
    for sd in sys_data_list:
        p = build_stiffness_html(sd, outdir)
        paths.append(p)
        print(f"  Generated {sd['name']}: {p}")
    # Generate index page
    index_html = '<!DOCTYPE html>\n<html><head><meta charset="utf-8"><title>Stiffness Maps Index</title>\n'
    index_html += '<style>body{background:#1a1a2e;color:#eee;font-family:sans-serif;padding:20px;}a{color:#6cf;text-decoration:none;font-size:16px;display:block;margin:8px 0;}a:hover{text-decoration:underline;}</style>\n'
    index_html += '</head><body><h1>Stiffness Maps</h1>\n'
    for p in paths:
        fname = os.path.basename(p)
        name = fname.replace('stiffness_map_', '').replace('.html', '')
        index_html += f'<a href="{fname}">{name}</a>\n'
    index_html += '</body></html>'
    index_path = os.path.join(outdir, 'index.html')
    with open(index_path, 'w') as f:
        f.write(index_html)
    print(f"  Index: {index_path}")
    return paths

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


def plot_hessian_comparison(H_ref, H_model, label, outdir):
    """Plot signed values and log magnitudes of reference, model, and residual Hessians."""
    os.makedirs(outdir, exist_ok=True)
    matrices = [H_ref, H_model, H_model - H_ref]
    names = ['DFT reference', 'FF model', 'model - DFT']
    vmax = max(np.max(np.abs(H_ref)), np.max(np.abs(H_model)))
    floor = max(vmax * 1e-10, np.finfo(float).tiny)
    log_matrices = [np.log10(np.maximum(np.abs(H), floor)) for H in matrices]
    log_min = np.log10(floor)
    log_max = np.log10(vmax)
    fig, axes = plt.subplots(2, 3, figsize=(15, 9), constrained_layout=True)
    for i, (H, name) in enumerate(zip(matrices, names)):
        im = axes[0, i].imshow(H, cmap='RdBu_r', vmin=-vmax, vmax=vmax, interpolation='nearest')
        axes[0, i].set_title(name)
        axes[0, i].set_xlabel('Cartesian coordinate')
        axes[0, i].set_ylabel('Cartesian coordinate')
        fig.colorbar(im, ax=axes[0, i], shrink=0.8, label='eV/A^2')
        imlog = axes[1, i].imshow(log_matrices[i], cmap='magma', vmin=log_min, vmax=log_max, interpolation='nearest')
        axes[1, i].set_title(f'log10 |{name}|')
        axes[1, i].set_xlabel('Cartesian coordinate')
        axes[1, i].set_ylabel('Cartesian coordinate')
        fig.colorbar(imlog, ax=axes[1, i], shrink=0.8, label='log10(eV/A^2)')
    out = os.path.join(outdir, f'{label}_hessian.png')
    fig.savefig(out, dpi=160)
    plt.close(fig)
    print(f"  Saved Hessian comparison: {out}")
    return out


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
    # LEGACY: retained for numerical comparison. New fits use
    # pyBall.FFfit.fit_hybrid_hessian(), which solves the direct stacked system.
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


def make_cpp_fitter(sys, maps, n_total):
    """Build one C++ linear bond/angle/1-4 sensitivity model with global indices."""
    bond_map, bond3_map, angle_map = maps
    f = FFfit_cpp.FFfit()
    f.set_geometry(sys['positions'])
    labels = sys.get('atom_types', sys['symbols'])
    f.set_symbols(labels)
    for i, j, l0 in sys['bonds']: f.add_bond(i, j, l0)
    for i, j, l0 in sys['bonds3']: f.add_bond(i, j, l0)
    for i, j, k, t0 in sys['angles']: f.add_angle(i, j, k, t0)
    nb = len(sys['bonds'])
    for ib, (i, j, _) in enumerate(sys['bonds'] + sys['bonds3']):
        key = tuple(sorted((labels[i], labels[j])))
        f.set_bond_param(ib, bond_map[key] if ib < nb else bond3_map[key])
    for ia, (i, j, k, _) in enumerate(sys['angles']):
        key = angle_type_key(labels, i, j, k, elements=sys['symbols'], central_only=sys.get('angle_central_only', False))
        f.set_angle_param(ia, angle_map[key])
    f.set_n_free(n_total)
    return f


def parent_parameter_name(name):
    """Collapse SiH3/SiH2/SiH labels to elemental Si in a parameter name."""
    prefix, body = name.split(':', 1)
    fields = ['Si' if field.startswith('SiH') else field for field in body.split('-')]
    return prefix + ':' + '-'.join(fields)


def fit_comparison_variant(base_systems, use_third, use_dihedrals, args, parent_prior=None):
    """Fit one interaction-set variant with local DFT equilibrium coordinates."""
    systems = []
    for base in base_systems:
        sys = dict(base)
        sys['bonds'] = list(base['bonds'])
        sys['angles'] = list(base['angles'])
        sys['bonds3'] = list(base['bonds3']) if use_third else []
        sys['dihedrals'] = list(base['dihedrals']) if use_dihedrals else []
        systems.append(sys)
    all_bonds = [s['bonds'] for s in systems]; all_angles = [s['angles'] for s in systems]
    all_bonds3 = [s['bonds3'] for s in systems]; all_dihedrals = [s['dihedrals'] for s in systems]
    all_symbols = [s.get('atom_types', s['symbols']) for s in systems]
    all_elements = [s['symbols'] for s in systems]
    angle_central_only = any(s.get('angle_central_only', False) for s in systems)
    if angle_central_only and not all(s.get('angle_central_only', False) for s in systems):
        raise ValueError("all systems in one fit must use the same angle typing scope")
    bmap, b3map, amap, dmap, npar = build_global_param_map(all_bonds, all_angles, all_symbols, all_bonds3, all_dihedrals, all_elements, angle_central_only)
    print(f"  parameter types: bond={len(bmap)} 1-4={len(b3map)} angle={len(amap)} torsion={len(dmap)} total={npar}")
    if args.equilibrium == 'type-average':
        avg_l0, avg_l0_3, avg_t0 = compute_averaged_equilibrium(all_bonds, all_angles, all_symbols, [s['positions'] for s in systems], bmap, amap, b3map, all_bonds3, all_elements, angle_central_only)
        for sys in systems:
            sys['bonds'], sys['angles'], sys['bonds3'] = rebuild_topology_with_averaged(sys['bonds'], sys['angles'], sys['bonds3'], sys.get('atom_types', sys['symbols']), sys['positions'], avg_l0, avg_t0, avg_l0_3, sys['symbols'], angle_central_only)
    for sys in systems:
        sys['dihedral_A'] = compute_dihedral_sensitivity(sys['positions'], sys.get('atom_types', sys['symbols']), sys['dihedrals'], dmap, h=args.dihedral_h) if use_dihedrals else {}
    print("  sensitivity stage: torsions complete")
    fitters = [make_cpp_fitter(sys, (bmap, b3map, amap), npar) for sys in systems]
    hybrid_systems = []
    for f, sys in zip(fitters, systems):
        hybrid_systems.append({'A': FFfit_cpp.collect_sensitivity_matrices(f, extra=sys['dihedral_A']), 'H_ref': sys['H_ref'],
                               'positions': sys['positions'], 'masses': sys['data']['masses'], 'bonds': sys['bonds'], 'angles': sys['angles']})
    print("  sensitivity stage: stacked model complete")
    names = [''] * npar
    for key, t in bmap.items(): names[t] = 'bond:' + '-'.join(key)
    for key, t in b3map.items(): names[t] = '1-4:' + '-'.join(key)
    for key, t in amap.items(): names[t] = 'angle:' + '-'.join(key)
    for key, t in dmap.items(): names[t] = 'torsion:' + '/'.join(('-'.join(key[0]), '-'.join(key[1])))
    prior = np.ones(npar)
    for t in bmap.values(): prior[t] = 5.0
    for t in b3map.values(): prior[t] = 0.1
    for t in amap.values(): prior[t] = 1.0
    for t in dmap.values(): prior[t] = 0.1
    if parent_prior is not None:
        for i, name in enumerate(names):
            parent = parent_parameter_name(name)
            if parent not in parent_prior:
                raise KeyError(f"missing elemental parent prior {parent} for subtype parameter {name}")
            prior[i] = parent_prior[parent]
    scale = np.maximum(np.abs(prior), 0.1)
    k, diag = FFfit_cpp.fit_hybrid_hessian(
        hybrid_systems, mode_weight=args.hybrid_mode, local_weight=args.hybrid_local, internal_weight=args.hybrid_internal,
        mode_balance=args.mode_weight, mode_mixing=args.mode_mixing, frequency_floor_cm1=args.frequency_floor,
        local_graph_distance=args.local_graph_distance, prior=prior, regularization=args.reg, parameter_scale=scale, bounds=(0.0, np.inf))
    print(f"  solve complete: condition={diag['condition']:.3g} residual={diag['relative_residual']:.3g}")
    models = []
    for f, sys in zip(fitters, systems):
        f.set_params(k)
        H = f.compute_model_hessian()
        add_dihedral_hessian(H, k, sys['dihedral_A'])
        models.append(H)
    return systems, models, k, names, diag


def one_to_one_frequency_metrics(freq_ref, freq_model):
    """One-to-one minimum-cost frequency assignment; unmatched counts are reported."""
    from scipy.optimize import linear_sum_assignment
    if len(freq_ref) == 0 or len(freq_model) == 0:
        return np.nan, np.nan, abs(len(freq_ref) - len(freq_model))
    ir, im = linear_sum_assignment(np.abs(freq_ref[:, None] - freq_model[None, :]))
    diff = freq_model[im] - freq_ref[ir]
    return np.sqrt(np.mean(diff*diff)), np.mean(np.abs(diff)), abs(len(freq_ref) - len(freq_model))


def plot_comparison_spectra(case_name, freq_ref, model_freqs, outdir, width=20.0, xmax=2500.0):
    """Overlay normalized DFT and fitted-model vibrational spectra."""
    os.makedirs(outdir, exist_ok=True)
    x = np.arange(0.0, xmax + 0.5, 0.5)
    fig, ax = plt.subplots(figsize=(13, 7))
    curves = [('DFT reference', freq_ref, '#111827', 2.5, '-')]
    styles = [('#0072B2', 1.8, '--'), ('#009E73', 1.8, '-.'), ('#D55E00', 1.7, ':'), ('#CC79A7', 1.8, (0, (5, 2)))]
    for (label, freqs), style in zip(model_freqs.items(), styles): curves.append((label, freqs, *style))
    for label, freqs, color, lw, ls in curves:
        y = broaden(freqs, x, width)
        if np.max(y) > 0.0: y /= np.max(y)
        ax.plot(x, y, color=color, lw=lw, ls=ls, label=label)
    ax.set(xlim=(0, xmax), ylim=(0, 1.08), xlabel='Frequency (cm$^{-1}$)', ylabel='Normalized density', title=f'{case_name}: progressive force-field vibration fit')
    ax.grid(True, alpha=0.2, ls=':')
    ax.legend(frameon=False, fontsize=9, ncol=2)
    out = os.path.join(outdir, f'{case_name}_model_comparison_spectra.png')
    fig.savefig(out, dpi=180, bbox_inches='tight')
    plt.close(fig)
    return out


def run_model_comparison(base_systems, args):
    """Fit and compare bond/angle, 1-4, torsion, and combined models."""
    variants = [
        ('Bond + angle', False, False),
        ('Bond + angle + 1-4', True, False),
        ('Bond + angle + UFF torsion', False, True),
        ('Bond + angle + UFF torsion + 1-4', True, True),
    ]
    fitted = {}
    for label, use_third, use_dihedrals in variants:
        print(f"\n=== Comparison model: {label} ===")
        fitted[label] = fit_comparison_variant(base_systems, use_third, use_dihedrals, args)
    rows = []
    os.makedirs(args.plot_dir, exist_ok=True)
    for isys, base in enumerate(base_systems):
        freq_ref = get_frequencies_cm1(base['H_ref'], base['data']['masses'], data=base['data'], freq_floor=10.0)
        model_freqs = {}
        rows.append((base['name'], 'DFT reference', 0, 0.0, 0.0, 0.0, 0, 0, 1.0, ''))
        for label, _, _ in variants:
            systems, models, k, names, diag = fitted[label]
            sys, H = systems[isys], models[isys]
            freqs = get_frequencies_cm1(H, sys['data']['masses'], freq_floor=10.0, positions=sys['positions'], project_acoustic=True)
            model_freqs[label] = freqs
            rmse, mae, unmatched = one_to_one_frequency_metrics(freq_ref, freqs)
            lam = FFfit_cpp.reference_vibrational_modes(H, sys['positions'], sys['data']['masses'])[1]
            negative = int(np.sum(lam < -(10.0/521.5)**2))
            rel_frob = 100.0 * np.linalg.norm(H - sys['H_ref']) / np.linalg.norm(sys['H_ref'])
            ptxt = '; '.join(f'{name}={value:.4g}' for name, value in zip(names, k))
            rows.append((base['name'], label, len(k), rmse, mae, rel_frob, unmatched, negative, diag['condition'], ptxt))
        out = plot_comparison_spectra(base['name'], freq_ref, model_freqs, args.plot_dir, args.spectrum_width, args.spectrum_xmax)
        print(f"  Saved model spectrum comparison: {out}")
    header = ('case', 'model', 'npar', 'freq_RMSE_cm-1', 'freq_MAE_cm-1', 'relFrob_percent', 'unmatched_modes', 'negative_modes', 'condition', 'parameters')
    csv_path = os.path.join(args.plot_dir, 'model_comparison.csv')
    md_path = os.path.join(args.plot_dir, 'model_comparison.md')
    with open(csv_path, 'w') as f:
        f.write(','.join(header) + '\n')
        for row in rows: f.write(','.join(f'"{x}"' if isinstance(x, str) else str(x) for x in row) + '\n')
    with open(md_path, 'w') as f:
        f.write('| Case | Model | Npar | RMSE cm^-1 | MAE cm^-1 | relFrob % | Unmatched | Negative | Condition |\n')
        f.write('|---|---|---:|---:|---:|---:|---:|---:|---:|\n')
        for row in rows:
            f.write(f'| {row[0]} | {row[1]} | {row[2]} | {row[3]:.3f} | {row[4]:.3f} | {row[5]:.3f} | {row[6]} | {row[7]} | {row[8]:.3g} |\n')
    print(f"\n{'Model':43s} {'N':>3s} {'RMSE':>9s} {'MAE':>9s} {'Frob%':>8s} {'miss':>5s} {'neg':>4s} {'cond':>8s}")
    for row in rows:
        if row[1] == 'DFT reference': continue
        print(f"{row[1]:43s} {row[2]:3d} {row[3]:9.2f} {row[4]:9.2f} {row[5]:8.2f} {row[6]:5d} {row[7]:4d} {row[8]:8.2f}")
    print(f"  Saved comparison table: {md_path}")
    return rows


def run_typing_comparison(base_systems, args):
    """Compare elemental and H-coordination Si typing using bond+angle models only."""
    elemental = []
    subtyped = []
    for base in base_systems:
        se = dict(base); se['atom_types'] = list(base['symbols']); se['angle_central_only'] = False; elemental.append(se)
        ss = dict(base); ss['atom_types'] = assign_si_environment_types(base['symbols'], base['bonds'], enabled=True); ss['angle_central_only'] = args.si_subtype_scope == 'central'; subtyped.append(ss)
    print("\n=== Elemental typing support ===")
    print_interaction_type_counts(elemental, args.min_type_count)
    print("\n=== Si environment typing support ===")
    print_interaction_type_counts(subtyped, args.min_type_count)
    elemental_fit = fit_comparison_variant(elemental, False, False, args)
    elemental_prior = dict(zip(elemental_fit[3], elemental_fit[2]))
    fits = {
        'Elemental Si typing': elemental_fit,
        'SiH3/SiH2/SiH/Si typing': fit_comparison_variant(subtyped, False, False, args, parent_prior=elemental_prior),
    }
    rows = []
    for isys, base in enumerate(base_systems):
        freq_ref = get_frequencies_cm1(base['H_ref'], base['data']['masses'], data=base['data'], freq_floor=10.0)
        model_freqs = {}
        for label, (systems, models, k, names, diag) in fits.items():
            sys, H = systems[isys], models[isys]
            freqs = get_frequencies_cm1(H, sys['data']['masses'], freq_floor=10.0, positions=sys['positions'], project_acoustic=True)
            model_freqs[label] = freqs
            rmse, mae, unmatched = one_to_one_frequency_metrics(freq_ref, freqs)
            rel_frob = 100.0*np.linalg.norm(H - sys['H_ref'])/np.linalg.norm(sys['H_ref'])
            rows.append((base['name'], label, len(k), rmse, mae, rel_frob, unmatched, diag['condition']))
        out = plot_comparison_spectra(base['name'] + '_typing', freq_ref, model_freqs, args.plot_dir, args.spectrum_width, args.spectrum_xmax)
        print(f"  Saved typing spectrum comparison: {out}")
    csv_path = os.path.join(args.plot_dir, 'typing_comparison.csv')
    md_path = os.path.join(args.plot_dir, 'typing_comparison.md')
    param_path = os.path.join(args.plot_dir, 'typing_parameters.csv')
    with open(csv_path, 'w') as f:
        f.write('case,typing,npar,freq_RMSE_cm-1,freq_MAE_cm-1,relFrob_percent,unmatched_modes,condition\n')
        for row in rows: f.write(','.join(f'"{x}"' if isinstance(x, str) else str(x) for x in row) + '\n')
    with open(md_path, 'w') as f:
        f.write('| Case | Typing | Npar | RMSE cm^-1 | MAE cm^-1 | relFrob % | Unmatched | Condition |\n')
        f.write('|---|---|---:|---:|---:|---:|---:|---:|\n')
        for row in rows: f.write(f'| {row[0]} | {row[1]} | {row[2]} | {row[3]:.3f} | {row[4]:.3f} | {row[5]:.3f} | {row[6]} | {row[7]:.3g} |\n')
        f.write('| **Mean** | **Elemental Si typing** | | {:.3f} | {:.3f} | {:.3f} | | |\n'.format(*[np.mean([r[i] for r in rows if r[1] == 'Elemental Si typing']) for i in (3, 4, 5)]))
        f.write('| **Mean** | **SiH3/SiH2/SiH/Si typing** | | {:.3f} | {:.3f} | {:.3f} | | |\n'.format(*[np.mean([r[i] for r in rows if r[1].startswith('SiH3')]) for i in (3, 4, 5)]))
    with open(param_path, 'w') as f:
        f.write('typing,parameter,value\n')
        for label, (_, _, k, names, _) in fits.items():
            for name, value in zip(names, k): f.write(f'"{label}","{name}",{value}\n')
    print(f"  Saved typing comparison table: {md_path}")
    print(f"  Saved typing parameters: {param_path}")
    return rows

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
    parser.add_argument('--equilibrium', choices=['local', 'type-average'], default='local', help='Use each relaxed internal coordinate (default) or one transferable equilibrium value per type')
    parser.add_argument('--plot-distributions', action='store_true', help='Plot relaxed bond/angle/1-4/torsion distributions by type')
    parser.add_argument('--plot-only', action='store_true', help='Only plot equilibrium distributions (skip Hessian loading and fitting)')
    parser.add_argument('--plot-stiffness', action='store_true', help='Plot DFT-extracted bond/angle stiffnesses via Wilson GF projection (requires Hessian)')
    parser.add_argument('--stiffness-html', action='store_true', help='Generate interactive p5.js HTML stiffness maps per nanocrystal (requires Hessian)')
    parser.add_argument('--compare-models', action='store_true', help='Fit the four progressive interaction sets and write a comparison table/spectrum overlay')
    parser.add_argument('--si-subtypes', action='store_true', help='Type Si atoms as SiH3 (including SiH4), SiH2, SiH, or bulk Si by bonded-H count')
    parser.add_argument('--si-subtype-scope', choices=['central', 'full'], default='central', help='Subtype angle centers only (default) or all three angle atoms')
    parser.add_argument('--compare-typing', action='store_true', help='Compare elemental and Si H-coordination typing using bond+angle models')
    parser.add_argument('--min-type-count', type=int, default=3, help='Flag interaction types represented by fewer observations')
    parser.add_argument('--mode-weight', type=str, default='frequency', choices=['equal', 'frequency', 'relative'], help='All-mode balance: equal Hessian accuracy, frequency-error balance (default), or relative eigenvalue accuracy')
    parser.add_argument('--hybrid-mode', type=float, default=1.0, help='Weight of the complete vibrational-mode objective')
    parser.add_argument('--hybrid-local', type=float, default=1.0, help='Weight of graph-local mass-weighted Hessian blocks')
    parser.add_argument('--hybrid-internal', type=float, default=1.0, help='Weight of the Wilson internal-coordinate projection')
    parser.add_argument('--mode-mixing', type=float, default=0.1, help='Relative penalty for off-diagonal mixing of reference modes')
    parser.add_argument('--local-graph-distance', type=int, default=2, help='Maximum bond-graph distance included in local Hessian blocks')
    parser.add_argument('--frequency-floor', type=float, default=50.0, help='Frequency floor (cm^-1) for stable mode weighting')
    parser.add_argument('--reg', type=float, default=2e-3, help='Dimensionless per-parameter prior regularization weight')
    parser.add_argument('--use-nnls', action='store_true', help='Deprecated compatibility flag; hybrid fitting is non-negative by default')
    parser.add_argument('--allow-negative', action='store_true', help='Allow unphysical negative stiffnesses for diagnostic comparison only')
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
    geo_only = args.plot_only
    print(f"=== Loading {len(case_names)} systems: {case_names} ===" + (' (geometry only)' if geo_only else ''))
    systems = []
    for name in case_names:
        case_dir = os.path.join(RESULTS_DIR, name)
        if not os.path.isdir(case_dir):
            print(f"  WARNING: skipping {name} (directory not found)")
            continue
        data = load_pyscf_case(case_dir, geometry_only=geo_only)
        natoms = data['natoms']
        if geo_only:
            H_ref = np.zeros((natoms*3, natoms*3))  # placeholder, not used
            H_ref_w = H_ref
            weight = None
        else:
            H_ref = hessian_ha_bohr_to_ev_ang2(
                data['hessian'].transpose(0, 2, 1, 3).reshape(natoms*3, natoms*3))
            if args.mass_weight:
                masses = data['masses']
                inv_sqrt_m = 1.0 / np.sqrt(np.repeat(masses, 3))
                weight = (inv_sqrt_m[:, None] * inv_sqrt_m[None, :]).flatten()
                H_ref_w = (inv_sqrt_m[:, None] * inv_sqrt_m[None, :]) * H_ref
            else:
                weight = None
                H_ref_w = H_ref
        bonds, angles, bonds3 = build_topology(data['symbols'], data['positions'], args.bond_cutoff,
                                                third_bonds=args.third_bonds or args.compare_models, third_bond_cutoff=args.third_bond_cutoff)
        atom_types = assign_si_environment_types(data['symbols'], bonds, enabled=args.si_subtypes or args.compare_typing)
        dihedrals = build_dihedrals(data['symbols'], data['positions'], bonds,
                                    d=args.dihedral_d, n=args.dihedral_n, dihedral=args.dihedrals or args.compare_models)
        sys = {
            'name': name, 'data': data, 'natoms': natoms, 'H_ref': H_ref, 'H_ref_w': H_ref_w,
            'bonds': bonds, 'angles': angles, 'bonds3': bonds3, 'dihedrals': dihedrals, 'weight': weight,
            'symbols': data['symbols'], 'atom_types': atom_types, 'angle_central_only': (args.si_subtypes or args.compare_typing) and args.si_subtype_scope == 'central', 'positions': data['positions'],
        }
        systems.append(sys)
        print(f"  {name}: natoms={natoms}, bonds={len(bonds)}, angles={len(angles)}, 3rd_bonds={len(bonds3)}, dihedrals={len(dihedrals)}, "
              f"H_ref [{H_ref.min():.2e}, {H_ref.max():.2e}]")
        if args.si_subtypes or args.compare_typing:
            unique, count = np.unique(atom_types, return_counts=True)
            print(f"    Si environment populations: {dict(zip(unique.tolist(), count.tolist()))}")
    
    if not systems:
        print("ERROR: no valid systems found"); sys.exit(1)

    print_interaction_type_counts(systems, min_count=args.min_type_count)

    if args.plot_only:
        plot_equilibrium_distributions(systems, args.plot_dir)
        return

    if args.plot_stiffness:
        plot_dft_stiffness_distributions(systems, args.plot_dir)
        return

    if args.stiffness_html:
        generate_stiffness_html(systems, args.plot_dir)
        return

    if args.compare_typing:
        plot_equilibrium_distributions(systems, args.plot_dir)
        run_typing_comparison(systems, args)
        return

    if args.compare_models:
        plot_equilibrium_distributions(systems, args.plot_dir)
        run_model_comparison(systems, args)
        return
    
    # === Build global parameter mapping ===
    all_bonds = [s['bonds'] for s in systems]
    all_angles = [s['angles'] for s in systems]
    all_bonds3 = [s['bonds3'] for s in systems]
    all_dihedrals = [s['dihedrals'] for s in systems]
    all_symbols = [s.get('atom_types', s['symbols']) for s in systems]
    all_elements = [s['symbols'] for s in systems]
    all_positions = [s['positions'] for s in systems]
    global_bond_map, global_bond3_map, global_angle_map, global_dihedral_map, n_total = build_global_param_map(
        all_bonds, all_angles, all_symbols, all_bonds3=all_bonds3, all_dihedrals=all_dihedrals,
        all_elements=all_elements, angle_central_only=args.si_subtypes and args.si_subtype_scope == 'central')
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
        global_bond_map, global_angle_map, global_bond3_map=global_bond3_map, all_bonds3=all_bonds3,
        all_elements=all_elements, angle_central_only=args.si_subtypes and args.si_subtype_scope == 'central')
    
    if args.plot_distributions or args.plot:
        plot_equilibrium_distributions(systems, args.plot_dir)

    # Local Hessian fitting shares stiffnesses but expands each interaction at its
    # actual relaxed DFT coordinate. Type-averaged minima are a separate energy-model choice.
    if args.equilibrium == 'type-average':
        for sys in systems:
            sys['bonds'], sys['angles'], sys['bonds3'] = rebuild_topology_with_averaged(
                sys['bonds'], sys['angles'], sys['bonds3'], sys.get('atom_types', sys['symbols']), sys['positions'],
                avg_l0, avg_theta0, avg_l0_3=avg_l0_3 if args.third_bonds else None, elements=sys['symbols'],
                angle_central_only=args.si_subtypes and args.si_subtype_scope == 'central')

    # === Precompute per-system dihedral sensitivity matrices (Python side) ===
    if args.dihedrals:
        print(f"\n=== Precomputing Dihedral Sensitivities ===")
        for sys in systems:
            sys['dihedral_A'] = compute_dihedral_sensitivity(
                sys['positions'], sys.get('atom_types', sys['symbols']), sys['dihedrals'], global_dihedral_map, h=args.dihedral_h)
            print(f"  {sys['name']}: computed {len(sys['dihedral_A'])} dihedral sensitivity matrices")
    else:
        for sys in systems:
            sys['dihedral_A'] = {}

    # === Create C++ FFfit instances with global param indices ===
    fitters = [make_cpp_fitter(sys, (global_bond_map, global_bond3_map, global_angle_map), n_total) for sys in systems]
    
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
    
    # === Method 3: robust hybrid all-mode/local/internal fit (Python) ===
    print(f"\n=== Method 3: Hybrid All-Mode + Local + Internal-Coordinate Fit ===")
    for f in fitters: f.set_n_free(n_total)
    hybrid_systems = []
    for f, sys in zip(fitters, systems):
        A = FFfit_cpp.collect_sensitivity_matrices(f, extra=sys.get('dihedral_A', {}))
        hybrid_systems.append({'A': A, 'H_ref': sys['H_ref'], 'positions': sys['positions'], 'masses': sys['data']['masses'],
                               'bonds': sys['bonds'], 'angles': sys['angles']})
    k_prior = np.ones(n_total)
    for t in global_bond_map.values(): k_prior[t] = 5.0
    for t in global_bond3_map.values(): k_prior[t] = 0.1
    for t in global_angle_map.values(): k_prior[t] = 1.0
    for t in global_dihedral_map.values(): k_prior[t] = 0.1
    bounds = (-np.inf, np.inf) if args.allow_negative else (0.0, np.inf)
    k_mode, fit_diag = FFfit_cpp.fit_hybrid_hessian(
        hybrid_systems, mode_weight=args.hybrid_mode, local_weight=args.hybrid_local, internal_weight=args.hybrid_internal,
        mode_balance=args.mode_weight, mode_mixing=args.mode_mixing, frequency_floor_cm1=args.frequency_floor, local_graph_distance=args.local_graph_distance,
        prior=k_prior, regularization=args.reg, parameter_scale=k_prior, bounds=bounds)
    print(f"  solver={fit_diag['solver']} residual={fit_diag['relative_residual']:.6e} rank={fit_diag['rank']}/{n_total} condition={fit_diag['condition']:.6e}")
    for sys, info in zip(systems, fit_diag['systems']):
        print(f"  {sys['name']}: all_modes={info['n_modes']} internal_rank={info['internal_rank']}/{info['n_internal_coordinates']} rows={info['component_rows']}")
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
            plot_hessian_comparison(sys['H_ref'], H_model, sys['name'], args.plot_dir)
        system_spectra_data.append((sys['name'], nz_ref, nz_model))

    if args.plot:
        plot_spectra_overlay(system_spectra_data, outdir=args.plot_dir, width=args.spectrum_width, xmax=args.spectrum_xmax)

    print(f"\n=== Done ===")

if __name__ == '__main__':
    main()
