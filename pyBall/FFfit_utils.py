"""High-level utilities for force-field fitting: topology, type system, dihedral physics.

These functions build on top of pyBall.FFfit (the C++ wrapper) to provide
reusable Python-side logic for molecular topology construction, atom typing,
dihedral sensitivity computation, parameter mapping, frequency analysis,
and multi-system fitting orchestration.

Separation of concerns:
  - pyBall.FFfit      → C++ wrapper + core numerical kernels (Wilson, hybrid, LSQ)
  - pyBall.FFfit_utils → high-level Python utilities (this module)
  - pyBall.FFfit_plots → visualization (separate module)
"""

from collections import deque
import numpy as np

# === Unit conversion constants ===
BOHR_TO_ANG = 0.5291772109
HARTREE_TO_EV = 27.211386245988
BOHR_TO_ANG_INV = 1.0 / BOHR_TO_ANG
HARTREE_PER_BOHR2_TO_EV_PER_ANG2 = HARTREE_TO_EV * BOHR_TO_ANG_INV**2

def hessian_ha_bohr_to_ev_ang2(H):
    """Convert Hessian from Hartree/Bohr^2 to eV/Angstrom^2."""
    return H * HARTREE_PER_BOHR2_TO_EV_PER_ANG2

# === Type system ===

def angle_type_key(symbols, i, j, k, elements=None, central_only=False):
    outer = elements if central_only else symbols
    si, sk = outer[i], outer[k]
    return (si, symbols[j], sk) if si <= sk else (sk, symbols[j], si)

def dihedral_type_key(symbols, i, j, k, l):
    """Type key for a dihedral (i-j-k-l): central bond + sorted outer atoms."""
    return (tuple(sorted([symbols[j], symbols[k]])),
            tuple(sorted([symbols[i], symbols[l]])))

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

def parent_parameter_name(name):
    """Collapse SiH3/SiH2/SiH labels to elemental Si in a parameter name."""
    prefix, body = name.split(':', 1)
    fields = ['Si' if field.startswith('SiH') else field for field in body.split('-')]
    return prefix + ':' + '-'.join(fields)

# === Topology ===

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
    """Find atom pairs separated by exactly 3 bonds (1-4) and add a distance spring."""
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

def build_topology(symbols, positions, bond_cutoff=1.85, third_bonds=False, third_bond_cutoff=None):
    """Build bond list and angle list from geometry using distance cutoff.
    
    Returns:
        bonds: list of (i, j, l0)
        angles: list of (i, j, k, theta0) where j is central
        bonds3: list of (i, j, l0) for 3rd-neighbor 1-4 springs (empty if disabled)
    """
    natoms = len(symbols)
    bonds = []
    bond_pairs = []
    for i in range(natoms):
        for j in range(i+1, natoms):
            if symbols[i] == 'H' and symbols[j] == 'H':
                continue
            r = np.linalg.norm(positions[j] - positions[i])
            if r < bond_cutoff:
                bonds.append((i, j, r))
                bond_pairs.append((i, j))
    neighs = [[] for _ in range(natoms)]
    for (i, j) in bond_pairs:
        neighs[i].append(j)
        neighs[j].append(i)
    angles = []
    for j in range(natoms):
        nn = neighs[j]
        for ii in range(len(nn)):
            for kk in range(ii+1, len(nn)):
                i = nn[ii]
                k = nn[kk]
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

def build_dihedrals(symbols, positions, bonds, d=1, n=3, dihedral=False):
    """Build proper torsion list from bond topology (i-j-k-l)."""
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
            continue
        for i in neighs[j]:
            if i == k:
                continue
            for l in neighs[k]:
                if l == j or l == i:
                    continue
                dihedrals.append((i, j, k, l, d, n))
    return dihedrals

# === Dihedral physics ===

def dihedral_energy_gradient(pos, d=1, n=3):
    """Energy and Cartesian gradient for a UFF/Prokop torsion term.

    Energy: E = V * (1 + d * cos(n * phi))   (V=1 here, so it returns dE/dV)
    Atoms: pos[0]=p1, pos[1]=p2, pos[2]=p3, pos[3]=p4  (i-j-k-l dihedral)
    Returns E, grad (4, 3) where grad = dE/dpos.
    This mirrors evalDihedral_Prokop() in cpp/common/molecular/UFF.h with
    bSubNonBond=False (no non-bonded subtraction).
    """
    p1, p2, p3, p4 = pos
    q12 = p1 - p2; q32 = p3 - p2; q43 = p4 - p3
    l12 = np.linalg.norm(q12); l32 = np.linalg.norm(q32); l43 = np.linalg.norm(q43)
    if l12 < 1e-12 or l32 < 1e-12 or l43 < 1e-12:
        return 0.0, np.zeros_like(pos)
    u12 = q12 / l12; u32 = q32 / l32; u43 = q43 / l43
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
    phi = np.arctan2(csy, csx)
    nphi = n * phi
    csnx = np.cos(nphi); csny = np.sin(nphi)
    E = 1.0 + d * csnx
    f = -d * n * csny
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
    """Cartesian Hessian (12,12) for one torsion term via central finite differences."""
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

# === Parameter mapping ===

def assign_types(symbols, bonds, angles):
    """Assign bond/angle types based on element pairs/triples."""
    bond_type_map = {}
    bond_types = []
    for (i, j, l0) in bonds:
        key = tuple(sorted([symbols[i], symbols[j]]))
        if key not in bond_type_map:
            bond_type_map[key] = len(bond_type_map)
        bond_types.append(bond_type_map[key])
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
        bmap = {}
        for ib, (i, j, l0) in enumerate(bonds):
            key = tuple(sorted([symbols[i], symbols[j]]))
            if key not in bmap:
                bmap[key] = pm.n_free; pm.n_free += 1
            pm.bond_param_idx[ib] = bmap[key]
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
    for bonds, symbols in zip(all_bonds, all_symbols):
        for (i, j, l0) in bonds:
            key = tuple(sorted([symbols[i], symbols[j]]))
            if key not in global_bond_map:
                global_bond_map[key] = len(global_bond_map)
    n_bond_types = len(global_bond_map)
    for bonds3, symbols in zip(all_bonds3, all_symbols):
        for (i, j, l0) in bonds3:
            key = tuple(sorted([symbols[i], symbols[j]]))
            if key not in global_bond3_map:
                global_bond3_map[key] = n_bond_types + len(global_bond3_map)
    n_bond3_types = len(global_bond3_map)
    for angles, symbols, elements in zip(all_angles, all_symbols, all_elements):
        for (i, j, k, theta0) in angles:
            key = angle_type_key(symbols, i, j, k, elements=elements, central_only=angle_central_only)
            if key not in global_angle_map:
                global_angle_map[key] = n_bond_types + n_bond3_types + len(global_angle_map)
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

# === Sensitivity matrices (Python reference) ===

def build_sensitivity_matrices(positions, bonds, angles, bond_types, angle_types, natoms):
    """Build dH/dk_t sensitivity matrices in Python (mirrors FFfit.h C++ logic).
    
    bond_types and angle_types are GLOBAL parameter indices (from ParamMap).
    Returns A: list of (3N x 3N) numpy arrays, one per free parameter.
    """
    n3 = natoms * 3
    nparams = max(int(bond_types.max()) + 1 if len(bond_types) else 0,
                  int(angle_types.max()) + 1 if len(angle_types) else 0)
    A = [np.zeros((n3, n3)) for _ in range(nparams)]
    for ib, (i, j, l0) in enumerate(bonds):
        t = bond_types[ib]
        d = positions[j] - positions[i]
        r = np.linalg.norm(d)
        e = d / r
        dl = r - l0
        P = np.eye(3) - np.outer(e, e)
        inv_r = 1.0 / r
        Cii = P * inv_r; Cjj = P * inv_r; Cij = -P * inv_r
        ii_idx = slice(i*3, i*3+3)
        jj_idx = slice(j*3, j*3+3)
        eeT = np.outer(e, e)
        A[t][ii_idx, ii_idx] += eeT + dl * Cii
        A[t][ii_idx, jj_idx] += -eeT + dl * Cij
        A[t][jj_idx, ii_idx] += -eeT + dl * Cij.T
        A[t][jj_idx, jj_idx] += eeT + dl * Cjj
    for ia, (i, j, k, theta0) in enumerate(angles):
        t = angle_types[ia]
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

def compute_gradient_term_blocks(positions, bonds, angles, param_map, k, H_ref, H0=None, weight=None):
    """Compute gradient of loss = ||H_model - H_ref||^2_W w.r.t. free parameters.

    Uses per-term sensitivity blocks (no full 3N×3N sensitivity matrices needed).
    Mirrors C++ FFfit::compute_gradient.
    """
    natoms = positions.shape[0]
    n3 = natoms * 3
    np_free = param_map.n_free
    H_model = np.zeros((n3, n3))
    bond_blocks = []
    for ib, (i, j, l0) in enumerate(bonds):
        d = positions[j] - positions[i]
        r = np.linalg.norm(d)
        e = d / r
        dl = r - l0
        P = np.eye(3) - np.outer(e, e)
        inv_r = 1.0 / r
        eeT = np.outer(e, e)
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
    angle_blocks = []
    for ia, (i, j, k_atom, theta0) in enumerate(angles):
        a_vec = positions[i] - positions[j]
        c_vec = positions[k_atom] - positions[j]
        al = np.linalg.norm(a_vec)
        cl = np.linalg.norm(c_vec)
        u = a_vec / al; v = c_vec / cl
        cos_t = np.dot(u, v)
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
        B = [[scale * np.outer(vecs[a], vecs[b]) for b in range(3)] for a in range(3)]
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
    h0 = H0 if H0 is not None else 0.0
    dH = H_model - H_ref + h0
    if weight is not None:
        W2 = weight * weight
        W2_mat = W2.reshape(n3, n3) if W2.size == n3*n3 else np.ones((n3, n3))
    else:
        W2_mat = np.ones((n3, n3))
    loss = np.sum(W2_mat * dH * dH)
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

# === Equilibrium averaging ===

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

# === Frequency analysis ===

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
        norms = np.linalg.norm(V, axis=0)
        V = V / norms[None, :]
        freqs_real = np.real(freqs)
        mask = freqs_real > freq_floor_cm1
        if not np.any(mask):
            mask = freqs_real > 0
        lam = (freqs_real[mask] / conv)**2
        return V[:, mask], lam, freqs_real[mask]
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
    cm = np.sum(positions * masses[:, None], axis=0) / np.sum(masses)
    r = positions - cm
    T = np.zeros((n3, 3))
    for a in range(3):
        T[a::3, a] = sqrt_m[a::3]
    R = np.zeros((n3, 3))
    e = np.eye(3)
    for a in range(3):
        for i in range(natoms):
            v = np.cross(e[a], r[i])
            R[i*3:i*3+3, a] = v * np.sqrt(masses[i])
    V_ac = np.hstack([T, R])
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
        diffs = np.min(np.abs(nz_ref[:, None] - nz_model[None, :]), axis=1)
        max_diff = np.max(diffs)
        print(f"  [{label}] First 10 ref:  {nz_ref[:10]}")
        print(f"  [{label}] First 10 model:{nz_model[:10]}")
        print(f"  [{label}] Max nearest-neighbour freq diff: {max_diff:.2f} cm^-1")
        print(f"  [{label}] Mean nearest-neighbour freq diff: {np.mean(diffs):.2f} cm^-1")
    return nz_ref, nz_model

def one_to_one_frequency_metrics(freq_ref, freq_model):
    """One-to-one minimum-cost frequency assignment; unmatched counts are reported."""
    from scipy.optimize import linear_sum_assignment
    if len(freq_ref) == 0 or len(freq_model) == 0:
        return np.nan, np.nan, abs(len(freq_ref) - len(freq_model))
    ir, im = linear_sum_assignment(np.abs(freq_ref[:, None] - freq_model[None, :]))
    diff = freq_model[im] - freq_ref[ir]
    return np.sqrt(np.mean(diff*diff)), np.mean(np.abs(diff)), abs(len(freq_ref) - len(freq_model))

# === C++ bridge ===

def make_cpp_fitter(sys, maps, n_total):
    """Build one C++ linear bond/angle/1-4 sensitivity model with global indices."""
    from pyBall import FFfit as FFfit_cpp
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

# === Mode-basis fitting (legacy) ===

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
    W = np.sqrt(weights)
    A_w = actions * W[None, None, :]
    B_w = B * W[None, :]
    A_mat = A_w.reshape(np_free, -1)
    B_vec = B_w.reshape(-1)
    G += A_mat @ A_mat.T
    y += A_mat @ B_vec
    if beta > 0.0 and H_ref is not None:
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
                                         reg=0.0, k0=None)
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
