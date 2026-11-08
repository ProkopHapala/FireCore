"""Spherical diamond/sphalerite nanoparticle builder + tetrahedral H caps. Parity: scripts/gen_nanocrystals.mjs, EditableMolecule.js."""
import os
import re
import numpy as np

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
RES = os.path.join(REPO_ROOT, 'cpp', 'common_resources')


def load_primitive_xyz_with_lvs(fname):
    with open(fname, 'r') as f:
        n = int(f.readline().strip())
        comment = f.readline().strip().split()
        if len(comment) < 10 or comment[0].lower() != 'lvs':
            raise ValueError(f"Expected 'lvs' line in {fname}, got: {comment}")
        lv = np.array([float(x) for x in comment[1:10]], dtype=np.float64).reshape(3, 3)
        elems, pos = [], []
        for i in range(n):
            ws = f.readline().split()
            if len(ws) < 4:
                raise ValueError(f"Bad XYZ line {i} in {fname}: {ws}")
            elems.append(ws[0])
            pos.append([float(ws[1]), float(ws[2]), float(ws[3])])
    return np.array(elems), np.asarray(pos, dtype=np.float64), lv


def parse_bond_l0(bond_types_path, z_heavy, z_cap=1):
    z1, z2 = min(z_heavy, z_cap), max(z_heavy, z_cap)
    best = None
    with open(bond_types_path, 'r') as f:
        for line in f:
            line = line.split('#')[0].strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            a, b = parts[0], parts[1]
            from pyBall.MMParams import element_symbol_to_Z
            try:
                za = element_symbol_to_Z(a) if len(a) <= 2 else None
                zb = element_symbol_to_Z(b) if len(b) <= 2 else None
            except Exception:
                continue
            if za is None or zb is None:
                # fallback: parse from ElementTypes by symbol
                za = _symbol_z(a)
                zb = _symbol_z(b)
            if za is None or zb is None:
                continue
            lo, hi = min(za, zb), max(za, zb)
            if lo != z1 or hi != z2:
                continue
            order = int(parts[2])
            l0 = float(parts[3])
            if best is None or (order == 1 and best[0] != 1) or (order == best[0] and l0 < best[1]):
                best = (order, l0)
    if best is None:
        raise ValueError(f'parse_bond_l0: no BondTypes entry for Z=({z1},{z2}) in {bond_types_path}')
    return best[1]


def _symbol_z(sym):
    tbl = {'H': 1, 'C': 6, 'N': 7, 'O': 8, 'Si': 14}
    return tbl.get(sym.strip(), None)


def parse_bond_l0_simple(bond_types_path, z_heavy, z_cap=1):
    """Parse BondTypes.dat without MMParams import."""
    z1, z2 = min(z_heavy, z_cap), max(z_heavy, z_cap)
    sym = {1: 'H', 6: 'C', 14: 'Si'}
    sa, sb = sym.get(z1, '?'), sym.get(z2, '?')
    best = None
    with open(bond_types_path, 'r') as f:
        for line in f:
            line = line.split('#')[0].strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            a, b = parts[0], parts[1]
            lo, hi = sorted([a, b])
            if not ((a == sa and b == sb) or (a == sb and b == sa)):
                if not ((lo == sa and hi == sb) or (lo == sb and hi == sa)):
                    continue
            order = int(parts[2])
            l0 = float(parts[3])
            if best is None or (order == 1 and best[0] != 1) or (order == best[0] and l0 < best[1]):
                best = (order, l0)
    if best is None:
        raise ValueError(f'parse_bond_l0_simple: no entry for {sa}-{sb}')
    return best[1]


def _orthonormal_basis_from_dir(d):
    d = np.asarray(d, float)
    dn = np.linalg.norm(d)
    if dn < 1e-12:
        raise ValueError('_orthonormal_basis_from_dir: zero direction')
    d = d / dn
    a, b, c = abs(d[0]), abs(d[1]), abs(d[2])
    tmp = np.array([1.0, 0.0, 0.0]) if a < 0.9 else (np.array([0.0, 1.0, 0.0]) if b < 0.9 else np.array([0.0, 0.0, 1.0]))
    u = tmp - d * np.dot(tmp, d)
    un = np.linalg.norm(u)
    if un < 1e-12:
        raise ValueError('_orthonormal_basis_from_dir: failed')
    u /= un
    v = np.cross(d, u)
    vn = np.linalg.norm(v)
    if vn < 1e-12:
        raise ValueError('_orthonormal_basis_from_dir: failed (v)')
    return u, v / vn


def missing_dirs_vsepr_tetra(vs, n_missing):
    """Parity: web/molgui_webgpu/EditableMolecule.js missingDirsVSEPR (totalDomains=4)."""
    vs = [np.asarray(v, dtype=np.float64) for v in vs]
    nb = len(vs)
    if n_missing <= 0:
        return []
    if nb == 3 and n_missing == 1:
        m = vs[0] + vs[1] + vs[2]
        rn = np.linalg.norm(m)
        if rn < 1e-12:
            raise ValueError('missing_dirs_vsepr_tetra: nb=3 sum is zero')
        return [-m / rn]
    if nb == 2 and n_missing == 2:
        m_c = vs[0] + vs[1]
        rn = np.linalg.norm(m_c)
        if rn < 1e-12:
            raise ValueError('missing_dirs_vsepr_tetra: nb=2 v1+v2 is zero')
        m_c /= rn
        m_b = np.cross(vs[0], vs[1])
        rb = np.linalg.norm(m_b)
        if rb < 1e-12:
            u, _ = _orthonormal_basis_from_dir(m_c)
            m_b = u
        else:
            m_b /= rb
        cc, cb = 0.57735026919, 0.81649658092
        d1 = -cc * m_c + cb * m_b
        d2 = -cc * m_c - cb * m_b
        d1 /= max(np.linalg.norm(d1), 1e-12)
        d2 /= max(np.linalg.norm(d2), 1e-12)
        return [d1, d2]
    if nb == 1 and n_missing == 3:
        v1 = vs[0] / max(np.linalg.norm(vs[0]), 1e-12)
        u, w = _orthonormal_basis_from_dir(v1)
        a, b = -1.0 / 3.0, np.sqrt(8.0 / 9.0)
        c120, s120 = -0.5, 0.86602540378
        u2 = c120 * u + s120 * w
        u3 = c120 * u - s120 * w
        out = []
        for d in (a * v1 + b * u, a * v1 + b * u2, a * v1 + b * u3):
            dn = np.linalg.norm(d)
            if dn < 1e-12:
                raise ValueError('missing_dirs_vsepr_tetra: nb=1 normalize failed')
            out.append(d / dn)
        return out
    raise ValueError(f'missing_dirs_vsepr_tetra: unsupported nb={nb} n_missing={n_missing}')


def resolve_cap_h_clashes(pos, elems, heavy_sym='C', min_dist=1.8, max_iter=16):
    """Rotate H around parent heavy atom to separate clashing caps. Parity: EditableMolecule.resolveCapHClashes."""
    pos = np.array(pos, dtype=np.float64, copy=True)
    elems = list(elems)
    n = len(elems)
    heavy_sym = str(heavy_sym)

    def parent_heavy(ih):
        for j in range(n):
            if j == ih or elems[j] != heavy_sym:
                continue
            if np.linalg.norm(pos[j] - pos[ih]) < 1.6:
                return j
        return -1

    min_d2 = min_dist * min_dist
    for _ in range(max_iter):
        n_fix = 0
        for ih in range(n):
            if elems[ih] != 'H':
                continue
            for jh in range(ih + 1, n):
                if elems[jh] != 'H':
                    continue
                d2 = np.sum((pos[ih] - pos[jh]) ** 2)
                if d2 >= min_d2:
                    continue
                for hi in (ih, jh):
                    oj = jh if hi == ih else ih
                    si = parent_heavy(hi)
                    if si < 0:
                        raise ValueError(f'resolve_cap_h_clashes: H {hi} has no parent')
                    s, h = pos[si], pos[hi].copy()
                    dvec = h - s
                    r = np.linalg.norm(dvec)
                    if r < 1e-12:
                        raise ValueError('resolve_cap_h_clashes: zero bond')
                    dvec /= r
                    rep = h - pos[oj]
                    rn = np.linalg.norm(rep)
                    if rn < 1e-12:
                        continue
                    rep /= rn
                    tan = rep - rep.dot(dvec) * dvec
                    tn = np.linalg.norm(tan)
                    if tn < 1e-12:
                        continue
                    tan /= tn
                    dvec = dvec + 0.12 * tan
                    dn = np.linalg.norm(dvec)
                    if dn < 1e-12:
                        raise ValueError('resolve_cap_h_clashes: rotate failed')
                    dvec /= dn
                    pos[hi] = s + dvec * r
                n_fix += 1
        if n_fix == 0:
            break
    return pos


def find_cap_hh_pairs(pos, elems, max_dist=1.8):
    """0-based H-H index pairs closer than max_dist (parity: EditableMolecule.addCapHHBonds)."""
    pos = np.asarray(pos, dtype=np.float64)
    n = len(elems)
    max_d2 = max_dist * max_dist
    pairs = []
    for ih in range(n):
        if elems[ih] != 'H':
            continue
        for jh in range(ih + 1, n):
            if elems[jh] != 'H':
                continue
            if np.sum((pos[ih] - pos[jh]) ** 2) < max_d2:
                pairs.append((ih, jh))
    return pairs


def save_hh_bonds_json(fname, pairs):
    import json
    with open(fname, 'w') as f:
        json.dump([list(p) for p in pairs], f, indent=0)


def default_rcut_heavy(heavy_z):
    return 2.55 if (heavy_z | 0) == 14 else 1.75


def build_spherical_nanoparticle(prim_xyz=None, R=10.0, nrep=8, min_degree=2, heavy_z=6,
                                 rcut_heavy=None, bond_types_path=None, cap='H', outward_bias=0.0,
                                 resolve_clashes=True):
    """Build spherical nanoparticle from primitive XYZ (diamond/sphalerite)."""
    heavy_z = int(heavy_z)
    if rcut_heavy is None:
        rcut_heavy = default_rcut_heavy(heavy_z)
    if prim_xyz is None:
        name = 'Si_primitive.xyz' if heavy_z == 14 else 'diamond_primitive.xyz'
        prim_xyz = os.path.join(RES, 'crystals', name)
    elems0, pos0, lvs = load_primitive_xyz_with_lvs(prim_xyz)
    heavy_sym = 'Si' if heavy_z == 14 else 'C'
    if not np.all(elems0 == heavy_sym):
        raise ValueError(f'{prim_xyz}: expected all {heavy_sym}, got {set(elems0)}')

    a1, a2, a3 = lvs
    shifts = []
    for i in range(-nrep, nrep + 1):
        for j in range(-nrep, nrep + 1):
            for k in range(-nrep, nrep + 1):
                shifts.append(i * a1 + j * a2 + k * a3)
    shifts = np.asarray(shifts, dtype=np.float64)
    pos = (pos0[None, :, :] + shifts[:, None, :]).reshape(-1, 3)
    pos = pos - pos.mean(axis=0)
    r2 = np.sum(pos * pos, axis=1)
    pos = pos[r2 <= R * R]
    nat = pos.shape[0]
    if nat < 10:
        raise ValueError(f'Too few atoms after sphere cut (nat={nat}) R={R} nrep={nrep}')

    cell = rcut_heavy
    pmin = pos.min(axis=0) - cell
    inv = 1.0 / cell
    key2atoms = {}
    icell = np.floor((pos - pmin) * inv).astype(np.int32)
    for ia in range(nat):
        key = (int(icell[ia, 0]), int(icell[ia, 1]), int(icell[ia, 2]))
        key2atoms.setdefault(key, []).append(ia)

    neigh = np.full((nat, 4), -1, dtype=np.int32)
    deg = np.zeros(nat, dtype=np.int32)
    rcut2 = rcut_heavy * rcut_heavy
    for ia in range(nat):
        cx, cy, cz = icell[ia]
        pi = pos[ia]
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    lst = key2atoms.get((int(cx + dx), int(cy + dy), int(cz + dz)))
                    if not lst:
                        continue
                    for ja in lst:
                        if ja <= ia:
                            continue
                        d = pos[ja] - pi
                        if np.dot(d, d) > rcut2:
                            continue
                        di, dj = deg[ia], deg[ja]
                        if di < 4 and dj < 4:
                            neigh[ia, di] = ja
                            deg[ia] = di + 1
                            neigh[ja, dj] = ia
                            deg[ja] = dj + 1

    alive = np.ones(nat, dtype=np.int8)
    queue = [int(i) for i in np.where(deg < min_degree)[0]]
    while queue:
        ia = queue.pop()
        if alive[ia] == 0:
            continue
        alive[ia] = 0
        for t in range(4):
            ja = int(neigh[ia, t])
            if ja < 0 or alive[ja] == 0:
                continue
            for u in range(4):
                if neigh[ja, u] == ia:
                    neigh[ja, u] = -1
                    deg[ja] -= 1
                    break
            if deg[ja] < min_degree:
                queue.append(ja)
        neigh[ia, :] = -1
        deg[ia] = 0

    keep = np.where(alive != 0)[0]
    pos = pos[keep]
    neigh = neigh[keep]
    deg = deg[keep]
    old2new = -np.ones(nat, dtype=np.int32)
    old2new[keep] = np.arange(len(keep), dtype=np.int32)
    for ia in range(neigh.shape[0]):
        for t in range(4):
            ja = neigh[ia, t]
            if ja >= 0:
                neigh[ia, t] = old2new[ja]
    nat = pos.shape[0]

    if bond_types_path is None:
        bond_types_path = os.path.join(RES, 'BondTypes.dat')
    ch = parse_bond_l0_simple(bond_types_path, heavy_z, 1)

    si_cog = pos.mean(axis=0) if outward_bias > 0 else None
    Hpos = []
    for ia in range(nat):
        pi = pos[ia]
        nnb = int(deg[ia])
        n_missing = 4 - nnb
        if n_missing <= 0:
            continue
        vs = []
        for t in range(4):
            ja = int(neigh[ia, t])
            if ja < 0:
                continue
            d = pos[ja] - pi
            r = np.linalg.norm(d)
            if r > 1e-8:
                vs.append(d / r)
        if len(vs) != nnb:
            raise ValueError(f'atom {ia}: len(vs)={len(vs)} deg={nnb}')
        dirs = missing_dirs_vsepr_tetra(vs, n_missing)
        if outward_bias > 0 and n_missing != 2:
            outward = pi - si_cog
            on = np.linalg.norm(outward)
            if on < 1e-12:
                raise ValueError(f'outward zero ia={ia}')
            outward /= on
            new_dirs = []
            for d in dirs:
                nd = (1 - outward_bias) * d + outward_bias * outward
                nn = np.linalg.norm(nd)
                if nn < 1e-12:
                    raise ValueError('outward blend failed')
                new_dirs.append(nd / nn)
            dirs = new_dirs
        for d in dirs:
            Hpos.append(pi + d * ch)

    if len(Hpos) == 0:
        raise ValueError('No undercoordinated atoms => no H passivation')
    Hpos = np.asarray(Hpos, dtype=np.float64)
    elems = np.array([heavy_sym] * nat + ['H'] * len(Hpos))
    apos = np.vstack([pos, Hpos])
    if resolve_clashes:
        apos = resolve_cap_h_clashes(apos, elems, heavy_sym)
    return elems, apos


def save_xyz(fname, elems, apos, comment=''):
    with open(fname, 'w') as f:
        f.write(f'{len(elems)}\n')
        f.write((comment or '') + '\n')
        for e, p in zip(elems, apos):
            f.write(f'{e:2s} {p[0]:.8f} {p[1]:.8f} {p[2]:.8f}\n')


# backward compat alias
def build_diamond_nanoparticle(R=10.0, nrep=8, min_degree=2, ch=None, rcut_cc=1.75):
    kw = dict(R=R, nrep=nrep, min_degree=min_degree, heavy_z=6, rcut_heavy=rcut_cc)
    if ch is not None:
        raise ValueError('build_diamond_nanoparticle: ch= is deprecated; bond length from BondTypes.dat')
    return build_spherical_nanoparticle(**kw)
