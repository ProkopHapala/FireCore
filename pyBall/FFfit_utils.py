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
from pathlib import Path
import json
import numpy as np

# === Unit conversion constants ===
BOHR_TO_ANG = 0.5291772109
HARTREE_TO_EV = 27.211386245988
HARTREE_TO_CM1 = 219474.6313705  # ħω in Hartree → wavenumber
BOHR_TO_ANG_INV = 1.0 / BOHR_TO_ANG
HARTREE_PER_BOHR2_TO_EV_PER_ANG2 = HARTREE_TO_EV * BOHR_TO_ANG_INV**2
# DFTB+ L1 run.py MASS — use for PDOS of those jobs (same AMU as the modes bake)
DFTBPLUS_AMU = {1: 1.008, 6: 12.011, 14: 28.085}
SYMBOL_Z = {"H": 1, "C": 6, "Si": 14}
PYSCF_SMALL_NC = Path("/home/prokop/SIMULATIONS/jobs_pyscf_vib_OUT_small_nc/results")
CCU_CAGE_PYSCF = Path("/home/prokop/git/CompChemUtils/examples/tSiNCs/jobs_dftb_vib_cages")

def hessian_ha_bohr_to_ev_ang2(H):
    """Convert Hessian from Hartree/Bohr^2 to eV/Angstrom^2."""
    return H * HARTREE_PER_BOHR2_TO_EV_PER_ANG2

# Canonical L1 PBE/ccECP-cc-pVDZ bank (2026-09-03 tight Si). Prefixed dirs: luna_ / tight_ / modefollow_.
PYSCF_VIB_L1_ROOT = Path("/home/prokop/SIMULATIONS/SiNCs/pyscf_vib_results")
PYSCF_VIB_L1_LEGACY = Path("/home/prokop/SIMULATIONS/SiNCs/pySCF/jobs/results")  # Batch-1; Si off-min — superseded
_PYSCF_VIB_BATCH_RANK = {"modefollow": 3, "tight": 2, "luna": 1, "": 0}

def pyscf_vib_job_case_name(job_dir):
    """Canonical crystal id: status.json['case'], else directory name with luna_/tight_/modefollow_ stripped."""
    p = Path(job_dir)
    case = None
    st = p / "status.json"
    if st.is_file():
        case = json.loads(st.read_text()).get("case")
    name = p.name
    batch = ""
    for pref in ("modefollow_", "tight_", "luna_"):
        if name.startswith(pref):
            batch = pref[:-1]
            name = name[len(pref):]
            break
    return str(case) if case else name, batch

def n_imag_from_omega_cm(omega_cm, cut=10.0):
    w = as_signed_wavenumbers_cm1(omega_cm)
    return int(np.sum(w < -float(cut))), float(w.min()) if w.size else float("nan")

def resolve_pyscf_vib_case_dir(root, case):
    """Best job dir for a canonical name (cube_Si, …). Prefers fewer imaginaries, then later batch (modefollow>tight>luna)."""
    root = Path(root)
    if not root.is_dir():
        raise FileNotFoundError(root)
    case = str(case)
    cand = []
    for p in root.iterdir():
        if not p.is_dir() or not (p / "frequencies_cm1.npy").is_file():
            continue
        name, batch = pyscf_vib_job_case_name(p)
        if name != case:
            continue
        n_imag, wmin = n_imag_from_omega_cm(np.load(p / "frequencies_cm1.npy"))
        cand.append((n_imag, -_PYSCF_VIB_BATCH_RANK.get(batch, 0), p, batch, wmin))
    if not cand:
        raise FileNotFoundError(f"{root}: no PySCF vib job for case={case!r} (looked for {case}/ and luna_/tight_/modefollow_ prefixes)")
    cand.sort(key=lambda t: (t[0], t[1]))
    return cand[0][2]

def list_pyscf_vib_cases(root):
    """Unique canonical case names under a results root (best-of prefixes)."""
    root = Path(root)
    if not root.is_dir():
        raise FileNotFoundError(root)
    names = set()
    for p in root.iterdir():
        if p.is_dir() and (p / "frequencies_cm1.npy").is_file():
            names.add(pyscf_vib_job_case_name(p)[0])
    return sorted(names)

def list_pyscf_vib_case_dirs(root):
    """One Path per canonical case (best batch)."""
    return [resolve_pyscf_vib_case_dir(root, c) for c in list_pyscf_vib_cases(root)]

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

def hydride_counts_from_coords(pos, Z, r_SiH=1.70, r_CH=1.25):
    """Count XH/XH2/XH3 groups from coordinates (distance cutoff, no topology file)."""
    pos = np.asarray(pos, float); Z = np.asarray(Z, int)
    H = np.flatnonzero(Z == 1)
    out = {'natoms': int(len(Z)), 'nH': int(len(H)), 'SiH': 0, 'SiH2': 0, 'SiH3': 0, 'CH': 0, 'CH2': 0, 'CH3': 0}
    for z, rcut, keys in ((14, r_SiH, ('SiH', 'SiH2', 'SiH3')), (6, r_CH, ('CH', 'CH2', 'CH3'))):
        heavies = np.flatnonzero(Z == z)
        if len(heavies) == 0 or len(H) == 0: continue
        d = np.linalg.norm(pos[heavies][:, None, :] - pos[H][None, :, :], axis=2)
        nH = (d < rcut).sum(axis=1)
        out[keys[0]] = int((nH == 1).sum()); out[keys[1]] = int((nH == 2).sum()); out[keys[2]] = int((nH >= 3).sum())
        bad = nH[nH > 4]
        if bad.size: raise ValueError(f"hydride_counts_from_coords: Z={z} atom has nH={int(bad.max())} > 4")
    return out

def xh_bonds_from_topology(Z, bonds_ij):
    """Group explicit X–H bonds by hydride class of the heavy (XH / XH2 / XH3). Distance cutoffs are forbidden here — 5-ring closers are topology, not a cutoff."""
    Z = np.asarray(Z, dtype=int).reshape(-1)
    from pyBall.io.crystal_npz import as_bonds_ij
    bonds = as_bonds_ij(bonds_ij, n_atoms=len(Z), ctx='xh_bonds_from_topology')
    nH = np.zeros(len(Z), dtype=int)
    pairs = []
    for i, j in bonds:
        zi, zj = int(Z[i]), int(Z[j])
        if zi == 1 and zj == 1:
            raise ValueError("xh_bonds_from_topology: H–H covalent bond")
        if zi == 1:
            i, j, zi, zj = j, i, zj, zi
        if zj != 1:
            continue
        if zi not in (6, 14):
            raise ValueError(f"xh_bonds_from_topology: H bonded to Z={zi} (want C or Si)")
        nH[i] += 1
        pairs.append((i, j))
    groups = {'XH': [], 'XH2': [], 'XH3': []}
    for i, j in pairs:
        nh = int(nH[i])
        if nh <= 0 or nh > 4:
            raise ValueError(f"xh_bonds_from_topology: heavy {i} has nH={nh}")
        key = 'XH3' if nh >= 3 else ('XH2' if nh == 2 else 'XH')
        groups[key].append((i, j))
    out = {}
    for k, v in groups.items():
        out[k] = np.asarray(v, dtype=np.int32).reshape(-1, 2) if v else np.zeros((0, 2), dtype=np.int32)
    return out, nH

def stretch_mode_weights(pos, modes_cart, bond_groups):
    """Project Cartesian modes onto X–H stretch: q_b = rhat·(u_H−u_X); W[name,mode] = Σ_b q_b²."""
    pos = np.asarray(pos, dtype=np.float64)
    U = np.asarray(modes_cart, dtype=np.float64)
    natoms = pos.shape[0]
    if pos.ndim != 2 or pos.shape[1] != 3:
        raise ValueError(f"stretch_mode_weights: pos.shape={pos.shape}")
    if U.ndim != 2 or U.shape[0] != 3 * natoms:
        raise ValueError(f"stretch_mode_weights: modes_cart.shape={U.shape} expected ({3*natoms}, nmode)")
    U3 = U.reshape(natoms, 3, U.shape[1])
    out = {}
    for name, bij in bond_groups.items():
        bij = np.asarray(bij, dtype=np.int32)
        w = np.zeros(U.shape[1], dtype=np.float64)
        if bij.size == 0:
            out[name] = w
            continue
        if bij.ndim != 2 or bij.shape[1] != 2:
            raise ValueError(f"stretch_mode_weights: {name} bonds shape {bij.shape}")
        for a, b in bij:
            r = pos[b] - pos[a]
            nrm = float(np.linalg.norm(r))
            if nrm < 1e-12:
                raise ValueError(f"stretch_mode_weights: zero X–H length atoms {a},{b}")
            q = ((r / nrm)[:, None] * (U3[b] - U3[a])).sum(axis=0)
            w += q * q
        out[name] = w
    return out

def gaussian_spectrum(omega_cm, grid, sigma, weights=None):
    """Gaussian-broadened DOS on a shared grid. Vectorized; used for overlay and ΔS."""
    om = np.asarray(omega_cm, dtype=np.float64).reshape(-1)
    g = np.asarray(grid, dtype=np.float64).reshape(-1)
    w = np.ones(om.shape[0], dtype=np.float64) if weights is None else np.asarray(weights, dtype=np.float64).reshape(-1)
    if w.shape != om.shape:
        raise ValueError(f"gaussian_spectrum: weights {w.shape} != omega {om.shape}")
    if sigma <= 0:
        raise ValueError(f"gaussian_spectrum: sigma={sigma}")
    if om.size == 0:
        return np.zeros(g.shape[0], dtype=np.float64)
    d = (g[:, None] - om[None, :]) / float(sigma)
    return np.exp(-0.5 * d * d) @ w

def weighted_mean_cm1(omega_cm, weights, lo=None, hi=None):
    om = np.asarray(omega_cm, dtype=np.float64).reshape(-1)
    w = np.asarray(weights, dtype=np.float64).reshape(-1)
    if w.shape != om.shape:
        raise ValueError("weighted_mean_cm1: shape mismatch")
    m = np.ones(om.shape[0], dtype=bool)
    if lo is not None: m &= om >= lo
    if hi is not None: m &= om <= hi
    ww = w[m]
    s = float(ww.sum())
    if s < 1e-18:
        return float('nan'), 0.0
    return float((om[m] * ww).sum() / s), s

def modes_as_3N(modes, natoms):
    """Cartesian modes as (3N, nmode). Accepts (nmode,N,3) DFTB/PySCF or (3N, nmode) pipeline."""
    m = np.asarray(modes, dtype=np.float64)
    n = int(natoms)
    if m.ndim == 3 and m.shape[1] == n and m.shape[2] == 3:
        return np.transpose(m, (1, 2, 0)).reshape(3 * n, m.shape[0])
    if m.ndim == 2 and m.shape[0] == 3 * n:
        return m
    raise ValueError(f"modes_as_3N: shape {m.shape} incompatible with natoms={n}")

def masses_amu_from_Z(Z, table=None):
    """Atomic masses (amu) from Z. Default table = DFTB+ L1 run.py MASS."""
    tab = DFTBPLUS_AMU if table is None else table
    Z = np.asarray(Z, dtype=int).reshape(-1)
    missing = sorted({int(z) for z in Z if int(z) not in tab})
    if missing:
        raise KeyError(f"masses_amu_from_Z: no AMU for Z={missing}")
    return np.array([tab[int(z)] for z in Z], dtype=np.float64)

def masses_amu_from_symbols(symbols, table=None):
    """Atomic masses (amu) from element symbols. Fail if a symbol is not in SYMBOL_Z."""
    missing = sorted({s for s in symbols if s not in SYMBOL_Z})
    if missing:
        raise KeyError(f"masses_amu_from_symbols: no Z for {missing}")
    return masses_amu_from_Z([SYMBOL_Z[s] for s in symbols], table=table)

def flatten_cartesian_hessian(H):
    """(N,N,3,3) or (N,3,N,3) or (3N,3N) → (3N,3N) with atom-xyz blocking H[3i+α, 3j+β] = H_ijαβ."""
    H = np.asarray(H, dtype=np.float64)
    if H.ndim == 2:
        if H.shape[0] != H.shape[1] or H.shape[0] % 3 != 0:
            raise ValueError(f"2D Hessian shape {H.shape} is not 3N×3N")
        return H
    if H.ndim != 4:
        raise ValueError(f"Hessian ndim={H.ndim} want 2 or 4, got shape {H.shape}")
    n0, n1, n2, n3 = H.shape
    if n2 == 3 and n3 == 3 and n0 == n1:
        n = n0
        return np.transpose(H, (0, 2, 1, 3)).reshape(3 * n, 3 * n)
    if n1 == 3 and n3 == 3 and n0 == n2:
        return np.reshape(H, (3 * n0, 3 * n0))
    raise ValueError(f"unrecognized 4D Hessian shape {H.shape}")

def read_xyz_symbols_positions(path):
    """Standard XYZ → (symbols, positions Å)."""
    lines = Path(path).read_text().splitlines()
    n = int(lines[0].strip())
    if len(lines) < n + 2:
        raise ValueError(f"{path}: expected {n} atom lines, file too short")
    symbols, pos = [], np.zeros((n, 3), dtype=np.float64)
    for i in range(n):
        parts = lines[2 + i].split()
        symbols.append(parts[0])
        pos[i] = [float(parts[1]), float(parts[2]), float(parts[3])]
    return symbols, pos

def load_pyscf_hessian_case(case_dir, geometry_only=False):
    """Load a PySCF vib directory: relaxed.xyz + hessian.npy (Ha/Bohr²). masses.npy optional."""
    d = Path(case_dir)
    xyz = d / "relaxed.xyz"
    if not xyz.is_file():
        raise FileNotFoundError(xyz)
    symbols, positions = read_xyz_symbols_positions(xyz)
    natoms = len(symbols)
    out = {"dir": d, "symbols": symbols, "positions": positions, "natoms": natoms}
    mp = d / "masses.npy"
    out["masses"] = np.load(mp).astype(np.float64) if mp.is_file() else masses_amu_from_symbols(symbols)
    if out["masses"].shape != (natoms,):
        raise ValueError(f"{d}: masses {out['masses'].shape} vs natoms={natoms}")
    st = d / "status.json"
    out["status"] = json.loads(st.read_text()) if st.is_file() else {}
    if geometry_only:
        return out
    hp = d / "hessian.npy"
    if not hp.is_file():
        raise FileNotFoundError(hp)
    out["hessian"] = np.load(hp)
    H3 = flatten_cartesian_hessian(out["hessian"])
    if H3.shape != (3 * natoms, 3 * natoms):
        raise ValueError(f"{d}: flattened Hessian {H3.shape} vs 3N={3 * natoms}")
    out["H_ref"] = hessian_ha_bohr_to_ev_ang2(H3)
    fp, mpde = d / "frequencies_cm1.npy", d / "modes.npy"
    if fp.is_file():
        out["freqs"] = np.load(fp)
    if mpde.is_file():
        out["modes"] = np.load(mpde)
    return out

def parse_dftb_tagged(path):
    """DFTB+ tagged data file → dict of numpy arrays. Rank≥2 is Fortran column-major (`order='F'`)."""
    from pathlib import Path
    tokens = Path(path).read_text().split()
    i, out = 0, {}
    while i < len(tokens):
        if i + 1 >= len(tokens):
            break
        spec = tokens[i + 1]
        if not (spec.startswith(':') and spec.count(':') >= 3):
            i += 1
            continue
        name = tokens[i]
        parts = spec.split(':')  # ['', type, rank, dims]
        typ, rank, dims = parts[1], int(parts[2]), parts[3]
        shape = tuple(int(x) for x in dims.split(',') if x)
        n = 1
        for s in shape:
            n *= s
        i += 2
        raw = tokens[i:i + n]
        if len(raw) < n:
            raise ValueError(f"parse_dftb_tagged: {path} field {name!r} short {len(raw)}<{n}")
        i += n
        if typ == 'real':
            arr = np.array(raw, dtype=np.float64)
        elif typ == 'integer':
            arr = np.array(raw, dtype=np.int64)
        else:
            raise ValueError(f"parse_dftb_tagged: {path} unsupported type {typ!r} for {name!r}")
        if rank >= 2:
            arr = arr.reshape(shape, order='F')
        elif rank == 1:
            arr = arr.reshape(shape[0] if shape else n)
        else:
            arr = arr.reshape(())
        if name in out:
            raise ValueError(f"parse_dftb_tagged: {path} duplicate field {name!r}")
        out[name] = arr
    if not out:
        raise ValueError(f"parse_dftb_tagged: {path} no tagged records")
    return out

def load_dftb_vibrations_tag(path, natoms):
    """Cartesian modes + signed cm⁻¹ from DFTB+ `vibrations.tag`.

    Tag `frequencies` are ħω in Hartree (convert with HARTREE_TO_CM1).
    `eigenmodes_scaled` columns are Euclidean-unit Cartesian displacements — same convention as
    L2 DFTB `modes.npy` after `modes_as_3N` (not the mass-weighted `eigenmodes` block).
    """
    rec = parse_dftb_tagged(path)
    n3 = 3 * int(natoms)
    if 'frequencies' not in rec or 'eigenmodes_scaled' not in rec:
        raise KeyError(f"{path}: need frequencies + eigenmodes_scaled, have {sorted(rec)}")
    omega = np.asarray(rec['frequencies'], dtype=np.float64).reshape(-1) * HARTREE_TO_CM1
    U = np.asarray(rec['eigenmodes_scaled'], dtype=np.float64)
    if U.ndim != 2 or U.shape[0] != n3 or U.shape[1] != omega.size:
        raise ValueError(f"{path}: eigenmodes_scaled {U.shape} vs natoms={natoms} nfreq={omega.size} (want ({n3}, nmode))")
    if omega.size != n3:
        raise ValueError(f"{path}: nfreq={omega.size} != 3N={n3} (RemoveTranslation/Rotation should still write 3N, six ~0)")
    phys = np.abs(omega) > 10.0
    if not np.any(phys):
        raise ValueError(f"{path}: no |ν|>10 cm⁻¹ modes")
    norms = np.linalg.norm(U[:, phys], axis=0)
    err = float(np.max(np.abs(norms - 1.0)))
    if err > 1e-6:
        raise ValueError(f"{path}: eigenmodes_scaled physical columns not unit Euclidean (max |‖u‖−1|={err:.3e})")
    if not np.all(np.isfinite(omega)) or not np.all(np.isfinite(U)):
        raise ValueError(f"{path}: non-finite frequencies or modes")
    return omega, U

def xh_bonds_from_coords(pos, Z, r_SiH=1.70, r_CH=1.25):
    """XH pairs from distance (L2 DFTB xyz with no 5-ring). Prefer xh_bonds_from_topology when bonds exist."""
    pos = np.asarray(pos, dtype=np.float64); Z = np.asarray(Z, dtype=int).reshape(-1)
    H = np.flatnonzero(Z == 1)
    pairs = []; nH = np.zeros(len(Z), dtype=int)
    for z, rcut in ((14, r_SiH), (6, r_CH)):
        heavies = np.flatnonzero(Z == z)
        if heavies.size == 0 or H.size == 0: continue
        d = np.linalg.norm(pos[heavies][:, None, :] - pos[H][None, :, :], axis=2)
        for a, row in enumerate(d):
            js = H[row < rcut]
            if js.size > 4:
                raise ValueError(f"xh_bonds_from_coords: Z={z} atom {heavies[a]} nH={js.size}")
            nH[heavies[a]] = js.size
            for j in js:
                pairs.append((int(heavies[a]), int(j)))
    groups = {'XH': [], 'XH2': [], 'XH3': []}
    for i, j in pairs:
        nh = int(nH[i])
        key = 'XH3' if nh >= 3 else ('XH2' if nh == 2 else 'XH')
        groups[key].append((i, j))
    out = {k: (np.asarray(v, dtype=np.int32).reshape(-1, 2) if v else np.zeros((0, 2), dtype=np.int32)) for k, v in groups.items()}
    return out, nH

def facet_kind_from_vec(r, families=('100', '111', '110')):
    """Nearest Wulff family via support function. Restrict `families` to the crystal's actual faces (cube has no {111} facets)."""
    r = np.asarray(r, dtype=np.float64).reshape(3)
    if float(np.linalg.norm(r)) < 1e-9:
        return 'bulk'
    scores = {}
    if '100' in families:
        scores['100'] = float(np.max(np.abs(r)))
    if '111' in families:
        n111 = np.array([[1.0, 1.0, 1.0], [1.0, 1.0, -1.0], [1.0, -1.0, 1.0], [1.0, -1.0, -1.0]])
        scores['111'] = float(np.max(np.abs(n111 @ r))) / np.sqrt(3.0)
    if '110' in families:
        n110 = []
        for i, j in ((0, 1), (0, 2), (1, 2)):
            for s1 in (1.0, -1.0):
                for s2 in (1.0, -1.0):
                    v = np.zeros(3); v[i] = s1; v[j] = s2
                    n110.append(v)
        scores['110'] = float(np.max(np.asarray(n110) @ r)) / np.sqrt(2.0)
    if not scores:
        raise ValueError("facet_kind_from_vec: empty families")
    return max(scores, key=scores.get)

def miller_111_unit_normals():
    """Eight outward ⟨111⟩ unit vectors in the crystal/Cartesian frame. Unrotated crystals only."""
    n = np.array([[s0, s1, s2] for s0 in (-1.0, 1.0) for s1 in (-1.0, 1.0) for s2 in (-1.0, 1.0)], dtype=np.float64)
    return n / np.sqrt(3.0)

def miller_110_unit_normals():
    """Twelve signed unit ⟨110⟩ in the crystal/Cartesian frame (unrotated)."""
    n = []
    for i, j in ((0, 1), (0, 2), (1, 2)):
        for s1 in (1.0, -1.0):
            for s2 in (1.0, -1.0):
                v = np.zeros(3); v[i] = s1; v[j] = s2
                n.append(v)
    return np.asarray(n, dtype=np.float64) / np.sqrt(2.0)

def miller_110_unsigned_axes():
    """Six unsigned ⟨110⟩ unit axes (sign ignored). Same 12 ends as ``miller_110_unit_normals``."""
    n = np.array([[1.0, 1.0, 0.0], [1.0, -1.0, 0.0], [1.0, 0.0, 1.0], [1.0, 0.0, -1.0], [0.0, 1.0, 1.0], [0.0, 1.0, -1.0]], dtype=np.float64)
    return n / np.sqrt(2.0)

def heavies_near_110_extrema(pos, Z, below_A=0.5):
    """C/Si within ``below_A`` (Å) of the max or min projection on any of the 6 ⟨110⟩ axes.

    Those atoms are the ⟨110⟩ silhouette — intended as edges/vertices, not {111} terraces.
    On a flat {110} rhombic face the whole facet sits at the extremum, so this may paint the face as edge.
    """
    pos = np.asarray(pos, dtype=np.float64); Z = np.asarray(Z, dtype=int).reshape(-1)
    if pos.shape != (len(Z), 3):
        raise ValueError(f"heavies_near_110_extrema: pos {pos.shape} vs N={len(Z)}")
    below_A = float(below_A)
    if below_A < 0.0:
        raise ValueError(f"heavies_near_110_extrema: below_A={below_A} must be ≥ 0")
    heavies = np.flatnonzero(Z > 1)
    if heavies.size < 1:
        raise ValueError("heavies_near_110_extrema: no C/Si atoms")
    ph = pos[heavies]
    near = np.zeros(len(Z), dtype=bool)
    for n in miller_110_unsigned_axes():
        s = ph @ n
        near[heavies] |= (float(s.max()) - s) <= below_A
        near[heavies] |= (s - float(s.min())) <= below_A
    return near

def xh_miller_111_cosine(uh):
    """max n·û over the 8 ⟨111⟩. ``uh`` is C→H (need not be unit)."""
    uh = np.asarray(uh, dtype=np.float64).reshape(3)
    nrm = float(np.linalg.norm(uh))
    if nrm < 1e-12:
        raise ValueError("xh_miller_111_cosine: zero X–H")
    return float(np.max(miller_111_unit_normals() @ (uh / nrm)))

def is_xh_on_miller_111(uh, r_from_com, align_cos=0.90):
    """X–H sits on a {111} face iff it points along that face’s ⟨111⟩ (cos ≥ align_cos).

    Sitting face = argmax_n (r·n) among the 8 axis-aligned ⟨111⟩ from the heavy’s COM vector
    (the octant the atom sits in). Then require û_XH · n_sit ≥ align_cos.
    Default 0.90 = cosine within 10% of 1; 0.95 is 5%. Unrotated crystals only.
    This is not Wulff morphology and does not replace ``is_xh_align_terrace`` / ``facet_kind_from_vec``.
    """
    r = np.asarray(r_from_com, dtype=np.float64).reshape(3)
    if float(np.linalg.norm(r)) < 1e-12:
        raise ValueError("is_xh_on_miller_111: zero r_from_com")
    uh = np.asarray(uh, dtype=np.float64).reshape(3)
    nrm = float(np.linalg.norm(uh))
    if nrm < 1e-12:
        raise ValueError("is_xh_on_miller_111: zero X–H")
    u = uh / nrm
    N = miller_111_unit_normals()
    k = int(np.argmax(N @ r))
    return float(N[k] @ u) >= float(align_cos)


def heavy_neighbor_lists(pos, Z, bonds_ij=None, r_CC=1.85, r_SiSi=2.50):
    """Heavy–heavy adjacency. Prefer explicit bonds; else C–C / Si–Si distance (not for 5-ring closers). Fail if a heavy has >4 neighbors."""
    pos = np.asarray(pos, dtype=np.float64); Z = np.asarray(Z, dtype=int).reshape(-1)
    n = len(Z)
    adj = [[] for _ in range(n)]
    if bonds_ij is not None:
        from pyBall.io.crystal_npz import as_bonds_ij
        for i, j in as_bonds_ij(bonds_ij, n_atoms=n, ctx='heavy_neighbor_lists'):
            if int(Z[i]) > 1 and int(Z[j]) > 1:
                adj[int(i)].append(int(j)); adj[int(j)].append(int(i))
    else:
        heavies = np.flatnonzero(Z > 1)
        for a in range(len(heavies)):
            i = int(heavies[a]); zi = int(Z[i])
            rcut = r_SiSi if zi == 14 else r_CC
            for b in range(a + 1, len(heavies)):
                j = int(heavies[b])
                if int(Z[j]) != zi:
                    continue
                if float(np.linalg.norm(pos[i] - pos[j])) < rcut:
                    adj[i].append(j); adj[j].append(i)
    for i in range(n):
        if int(Z[i]) > 1 and len(adj[i]) > 4:
            raise ValueError(f"heavy_neighbor_lists: atom {i} Z={int(Z[i])} has {len(adj[i])} heavy neighbors")
        adj[i] = sorted(set(adj[i]))
    return adj


def xh_unit_dirs(pos, xh_groups):
    """Heavy → (n_XH, 3) unit X–H vectors. Fail on zero length."""
    pos = np.asarray(pos, dtype=np.float64)
    dirs = {}
    for bij in xh_groups.values():
        for i, j in np.asarray(bij, dtype=np.int32).reshape(-1, 2):
            r = pos[int(j)] - pos[int(i)]
            nrm = float(np.linalg.norm(r))
            if nrm < 1e-12:
                raise ValueError(f"xh_unit_dirs: zero X–H {i},{j}")
            dirs.setdefault(int(i), []).append(r / nrm)
    return {k: np.asarray(v, dtype=np.float64).reshape(-1, 3) for k, v in dirs.items()}


def _xh_dirs_aligned(A, B, align_cos):
    """Every bond in A has a partner in B with |n·n'| >= align_cos, and vice versa."""
    A = np.asarray(A, dtype=np.float64).reshape(-1, 3)
    B = np.asarray(B, dtype=np.float64).reshape(-1, 3)
    dots = np.abs(A @ B.T)
    return bool(np.all(dots.max(axis=1) >= align_cos) and np.all(dots.max(axis=0) >= align_cos))


def is_xh_align_terrace(i, nH, adj, dirs, align_cos=0.9, max_misaligned=1, xh2_rim_terrace=False):
    """Terrace vs ridge from hydride-neighbor chemistry (not Wulff).

    Diamond {111}: terrace CH bond only to bulk C, not to other hydrides — isolated => terrace.
    Diamond {100}: terrace CH₂ are not bonded to each other; one CH at the face rim is allowed (``xh2_rim_terrace``).
    Do not use that CH₂-rim rule on rhombic {110}: those CH₂ are vertices, not the face.
    Ridge: tetrahedral X–H (|n·n'|~1/3) or too many mixed XH/XH₂ contacts.
    """
    i = int(i)
    nh = int(nH[i])
    if nh <= 0:
        return False
    if i not in dirs:
        raise ValueError(f"is_xh_align_terrace: heavy {i} has nH={nh} but no X–H vectors")
    hyd_nbs = [int(j) for j in adj[i] if int(nH[j]) > 0]
    if not hyd_nbs:
        return True
    n_good = 0
    n_bad = 0
    n_same = 0
    for j in hyd_nbs:
        if j not in dirs:
            raise ValueError(f"is_xh_align_terrace: neighbor {j} has nH={int(nH[j])} but no X–H vectors")
        if int(nH[j]) != nh:
            n_bad += 1
            continue
        n_same += 1
        if _xh_dirs_aligned(dirs[i], dirs[j], align_cos):
            n_good += 1
        else:
            n_bad += 1
    if n_good >= 1 and n_bad <= int(max_misaligned):
        return True
    if xh2_rim_terrace and nh >= 2 and n_same == 0 and n_bad <= int(max_misaligned):
        return True
    return False


def wulff_families_from_name(name):
    s = str(name).lower()
    if 'trunc' in s:
        return ('100', '111', '110')
    if 'rhombic' in s:
        return ('110',)
    if 'cube' in s:
        return ('100', '110')
    if 'octa' in s:
        return ('111', '110')
    return ('100', '111', '110')

def heavy_cycles(Z, bonds_ij, length):
    """Simple cycles of given length on C/Si only (5-ring / 7-ring atoms)."""
    from pyBall.io.crystal_npz import as_bonds_ij
    Z = np.asarray(Z, dtype=int).reshape(-1)
    bonds = as_bonds_ij(bonds_ij, n_atoms=len(Z), ctx='heavy_cycles')
    adj = {i: [] for i in range(len(Z)) if int(Z[i]) > 1}
    for i, j in bonds:
        if i not in adj or j not in adj:
            continue
        adj[i].append(int(j)); adj[j].append(int(i))
    cycles = []
    for start in sorted(adj):
        stack = [(start, [start])]
        while stack:
            node, path = stack.pop()
            if len(path) == length:
                if start in adj[path[-1]] and path[1] < path[-1]:
                    cycles.append(tuple(path))
                continue
            for nb in adj[node]:
                if nb > start and nb not in path:
                    stack.append((nb, path + [nb]))
    return cycles

def neighborhood_xh_groups(pos, Z, xh_groups, nH, bonds_ij=None, ring_lengths=(), families=('100', '111', '110'),
                           facet_mode='wulff', face_families=None, align_cos=0.9, max_misaligned=1, ridge_below_A=0.5):
    """Split X–H bonds into XH@{111}, XH2@{100}, ring5, ring7, plus leftovers.

    ``facet_mode='wulff'``: Miller family from COM support (``facet_kind_from_vec``).
    ``facet_mode='xh_align'``: terrace vs edge from neighbor hydride class + X–H alignment;
    terrace Miller index still from Wulff among ``face_families`` (primary faces only).
    ``facet_mode='miller_111'``: sitting {111} octant from heavy−COM, then X–H vs that same ⟨111⟩ (``is_xh_on_miller_111``).
    ``facet_mode='ridge_110'``: C/Si within ridge_below_A of a ⟨110⟩ extremity are edge; leftover Miller from ``face_families``.
    Does not replace wulff / xh_align / miller_111.
    """
    if facet_mode not in ('wulff', 'xh_align', 'miller_111', 'ridge_110'):
        raise ValueError(f"neighborhood_xh_groups: facet_mode={facet_mode!r} want wulff|xh_align|miller_111|ridge_110")
    pos = np.asarray(pos, dtype=np.float64); Z = np.asarray(Z, dtype=int).reshape(-1)
    nH = np.asarray(nH, dtype=int).reshape(-1)
    com = pos[Z > 1].mean(axis=0) if np.any(Z > 1) else pos.mean(axis=0)
    ring5 = set(); ring7 = set()
    if bonds_ij is not None:
        if 5 in ring_lengths:
            for cyc in heavy_cycles(Z, bonds_ij, 5):
                ring5.update(cyc)
        if 7 in ring_lengths:
            for cyc in heavy_cycles(Z, bonds_ij, 7):
                ring7.update(cyc)
    adj = dirs = None
    ridge_mask = None
    ff = face_families if face_families is not None else tuple(f for f in families if f != '110') or families
    if facet_mode == 'xh_align':
        adj = heavy_neighbor_lists(pos, Z, bonds_ij=bonds_ij)
        dirs = xh_unit_dirs(pos, xh_groups)
    elif facet_mode == 'ridge_110':
        ridge_mask = heavies_near_110_extrema(pos, Z, below_A=ridge_below_A)
    out = {}
    def _add(name, pair):
        out.setdefault(name, []).append(pair)
    for cls, bij in xh_groups.items():
        for i, j in np.asarray(bij, dtype=np.int32).reshape(-1, 2):
            if facet_mode == 'wulff':
                fac = facet_kind_from_vec(pos[int(i)] - com, families=families)
            elif facet_mode == 'miller_111':
                fac = '111' if is_xh_on_miller_111(pos[int(j)] - pos[int(i)], pos[int(i)] - com, align_cos=align_cos) else 'edge'
            elif facet_mode == 'ridge_110':
                fac = 'edge' if ridge_mask[int(i)] else facet_kind_from_vec(pos[int(i)] - com, families=ff)
            elif is_xh_align_terrace(int(i), nH, adj, dirs, align_cos=align_cos, max_misaligned=max_misaligned, xh2_rim_terrace=('100' in ff)):
                fac = facet_kind_from_vec(pos[int(i)] - com, families=ff)
            else:
                fac = 'edge'
            _add(f"{cls}@{fac}", (int(i), int(j)))
            if int(i) in ring5:
                _add('ring5', (int(i), int(j)))
            if int(i) in ring7:
                _add('ring7', (int(i), int(j)))
    empty = np.zeros((0, 2), dtype=np.int32)
    return {k: (np.asarray(v, dtype=np.int32).reshape(-1, 2) if v else empty) for k, v in out.items()}, {'ring5': sorted(ring5), 'ring7': sorted(ring7)}

# Face ≠ edge hues. {111} CH = blue, {100} CH₂ = red, {110} face (rhombic) = purple, ridge = orange / brown.
NBHD_TAG_RGB = {
    "XH@111": np.array([0.13, 0.44, 0.71]), "XH2@100": np.array([0.80, 0.09, 0.11]),
    "XH@100": np.array([0.42, 0.68, 0.84]), "XH2@111": np.array([0.98, 0.42, 0.29]),
    "XH@110": np.array([0.42, 0.24, 0.60]), "XH2@110": np.array([0.68, 0.00, 0.40]),
    "XH@edge": np.array([0.90, 0.33, 0.05]), "XH2@edge": np.array([0.55, 0.32, 0.04]),
    "XH3@100": np.array([0.55, 0.35, 0.15]), "XH3@111": np.array([0.55, 0.35, 0.15]), "XH3@110": np.array([0.55, 0.35, 0.15]),
    "XH3@edge": np.array([0.55, 0.35, 0.15]),
    "bulk": np.array([0.78, 0.78, 0.78]),
    "ring5": np.array([0.14, 0.55, 0.27]), "ring7": np.array([0.85, 0.28, 0.00]),
}
NBHD_PLOT_STYLE = {
    "bulk": dict(color="#c6c6c6", ls="-"),
    "XH@111": dict(color="#2171b5", ls="-"), "XH2@100": dict(color="#cb181d", ls="-"),
    "XH@100": dict(color="#6baed6", ls="--"), "XH2@111": dict(color="#fb6a4a", ls="--"),
    "XH@110": dict(color="#6a51a3", ls="-"), "XH2@110": dict(color="#ae017e", ls="-"),
    "XH@edge": dict(color="#e6550d", ls=":"), "XH2@edge": dict(color="#8c510a", ls=":"),
    "ring5": dict(color="#238b45", ls="-"), "ring7": dict(color="#d94801", ls="-."),
}
_HYDRIDE_NAME = {"C": {"XH": "CH", "XH2": "CH₂", "XH3": "CH₃"}, "Si": {"XH": "SiH", "XH2": "SiH₂", "XH3": "SiH₃"}}

def primary_face_families_from_name(name):
    """Wulff *faces* (not edges). Cube {100}; octa {111}; trunc {100}+{111}; rhombic {110} is the face."""
    s = str(name).lower()
    if 'rhombic' in s:
        return ('110',)
    if 'trunc' in s:
        return ('100', '111')
    if 'cube' in s:
        return ('100',)
    if 'octa' in s:
        return ('111',)
    return ('100', '111', '110')

def nbhd_legend_label(tag, elem='C', face_families=('100', '111')):
    """e.g. CH@{111} face, CH@{110} edge. {110} is a face only on rhombic (pass face_families=('110',))."""
    if tag in ('ring5', 'ring7'):
        return '5-ring' if tag == 'ring5' else '7-ring'
    if tag == 'bulk':
        return 'bulk / core'
    if '@' not in tag:
        return tag
    cls, fac = tag.split('@', 1)
    hyd = _HYDRIDE_NAME.get(elem, _HYDRIDE_NAME['C']).get(cls, cls)
    if fac == 'edge':
        return f"{hyd} edge"
    kind = 'face' if fac in face_families else 'edge'
    return f"{hyd}@{{{fac}}} {kind}"

def halogen_symbol_for_nbhd_tag(tag, face_families):
    """Jmol: F = XH face, Cl = XH edge, Br = XH₂ face, I = XH₂ edge. Miller face vs edge from ``face_families``."""
    if tag in ('bulk', 'ring5', 'ring7') or '@' not in str(tag):
        raise ValueError(f"halogen_symbol_for_nbhd_tag: not a hydride tag {tag!r}")
    cls, fac = str(tag).split('@', 1)
    is_face = (fac != 'edge') and (fac in face_families)
    if cls == 'XH':
        return 'F' if is_face else 'Cl'
    if cls == 'XH2':
        return 'Br' if is_face else 'I'
    if cls == 'XH3':
        return 'At'
    raise ValueError(f"halogen_symbol_for_nbhd_tag: class {cls}")

def write_hydride_halogen_xyz(path, pos, Z, owner, face_families, comment=''):
    """Replace each H with a halogen so Jmol element colors show neighborhood groups. Heavies unchanged."""
    pos = np.asarray(pos, dtype=np.float64); Z = np.asarray(Z, dtype=int).reshape(-1)
    n = len(Z)
    if pos.shape != (n, 3):
        raise ValueError(f"write_hydride_halogen_xyz: pos {pos.shape} vs N={n}")
    heav = {6: 'C', 14: 'Si'}
    counts = {}
    body = []
    for i in range(n):
        zi = int(Z[i])
        if zi != 1:
            if zi not in heav:
                raise ValueError(f"write_hydride_halogen_xyz: unexpected Z={zi} atom {i}")
            body.append(f"{heav[zi]}  {pos[i,0]:.6f}  {pos[i,1]:.6f}  {pos[i,2]:.6f}")
            continue
        tag = owner.get(i)
        if tag is None:
            raise ValueError(f"write_hydride_halogen_xyz: H {i} has no neighborhood tag")
        sym = halogen_symbol_for_nbhd_tag(tag, face_families)
        key = f"{sym}:{tag}"
        counts[key] = counts.get(key, 0) + 1
        body.append(f"{sym}  {pos[i,0]:.6f}  {pos[i,1]:.6f}  {pos[i,2]:.6f}")
    bits = " ".join(f"{k}={v}" for k, v in sorted(counts.items()))
    hdr = "F=XH-face Cl=XH-edge Br=XH2-face I=XH2-edge. " + bits
    if comment:
        hdr = str(comment) + " | " + hdr
    from pathlib import Path
    p = Path(path); p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(f"{n}\n{hdr}\n" + "\n".join(body) + "\n")
    return p

def cartesian_modes_from_hessian(H, masses, freq_floor_cm1=10.0):
    """Mass-weighted eigh (eV/Å²). Returns vibrational ω(cm⁻¹) and Cartesian modes. Fails if the Hessian is indefinite (|ν|>floor)."""
    conv = 521.5
    masses = np.asarray(masses, dtype=np.float64).reshape(-1)
    inv_sqrt_m = 1.0 / np.sqrt(np.repeat(masses, 3))
    D = (inv_sqrt_m[:, None] * inv_sqrt_m[None, :]) * np.asarray(H, dtype=np.float64)
    lam, V = np.linalg.eigh(D)
    om_all = np.sign(lam) * conv * np.sqrt(np.abs(lam))
    assert_harmonic_spectrum_at_minimum(om_all, ctx="cartesian_modes_from_hessian: ")
    mask = om_all > float(freq_floor_cm1)
    if not np.any(mask):
        raise RuntimeError("cartesian_modes_from_hessian: no vibrational modes above floor — not a spectrum")
    return om_all[mask], (V * inv_sqrt_m[:, None])[:, mask]

def atom_tags_from_owner(n, owner):
    """(N,) object tags; untagged atoms are 'bulk'. Partition of {0..N-1}."""
    tags = np.empty(int(n), dtype=object)
    tags[:] = 'bulk'
    for i, t in owner.items():
        tags[int(i)] = t
    return tags

def apply_ring_tags(owner, nbhd, rings):
    """Exclusive PDOS partition: ring heavies (``heavy_cycles``) and their X–H hydrogens steal hydride tags.

    ``rings`` / ``nbhd['ring5'|'ring7']`` come from ``neighborhood_xh_groups`` (do not re-detect).
    """
    owner = {int(i): t for i, t in owner.items()}
    s5 = {int(i) for i in rings.get('ring5', ())}
    s7 = {int(i) for i in rings.get('ring7', ())}
    both = s5 & s7
    if both:
        raise ValueError(f"apply_ring_tags: atoms in both 5- and 7-ring: {sorted(both)}")
    for tag, atoms in (('ring5', s5), ('ring7', s7)):
        for i in atoms:
            owner[i] = tag
        bij = nbhd.get(tag)
        if bij is None or len(bij) == 0:
            continue
        for i, j in np.asarray(bij, dtype=np.int32).reshape(-1, 2):
            owner[int(i)] = tag
            owner[int(j)] = tag
    return owner

def atom_group_mode_weights(modes_cart, masses, atom_tags):
    """Per-mode mass-weighted |u|² partitioned by atom tag. Each returned array sums with the others to 1.

    p_i = m_i |u_i|² / Σ_j m_j |u_j|². Groups must partition the atoms. Fail-loud if the partition is incomplete.
    """
    U = np.asarray(modes_cart, dtype=np.float64)
    m = np.asarray(masses, dtype=np.float64).reshape(-1)
    tags = np.asarray(atom_tags).reshape(-1)
    n = m.size
    if tags.shape[0] != n:
        raise ValueError(f"atom_group_mode_weights: tags {tags.shape} vs N={n}")
    if U.ndim != 2 or U.shape[0] != 3 * n:
        raise ValueError(f"atom_group_mode_weights: modes_cart {U.shape} want ({3*n}, nmode)")
    if np.any(~np.isfinite(m)) or np.any(m <= 0):
        raise ValueError("atom_group_mode_weights: masses must be finite and positive")
    U3 = U.reshape(n, 3, U.shape[1])
    p = m[:, None] * np.sum(U3 * U3, axis=1)
    s = p.sum(axis=0)
    bad = s < 1e-18
    if np.any(bad):
        raise ValueError(f"atom_group_mode_weights: {int(bad.sum())} zero-norm modes")
    p = p / s[None, :]
    out = {}
    for t in tags:
        key = str(t)
        if key in out:
            continue
        out[key] = p[tags == t].sum(axis=0)
    tot = np.zeros(U.shape[1], dtype=np.float64)
    for v in out.values():
        tot += v
    err = float(np.max(np.abs(tot - 1.0)))
    if err > 1e-6:
        raise ValueError(f"atom_group_mode_weights: group weights do not sum to 1 (max |Δ|={err:.3e})")
    return out

def atom_stretch_in_window(pos, modes_cart, om_cm, bonds_xh, lo, hi):
    """Per-atom Σ q² of listed X–H bonds, summed over modes in [lo,hi]."""
    pos = np.asarray(pos, dtype=np.float64)
    U = np.asarray(modes_cart, dtype=np.float64)
    om = np.asarray(om_cm, dtype=np.float64).reshape(-1)
    natoms = pos.shape[0]
    if U.shape[0] != 3 * natoms or U.shape[1] != om.shape[0]:
        raise ValueError(f"atom_stretch_in_window: U {U.shape} om {om.shape} N={natoms}")
    w = np.zeros(natoms, dtype=np.float64)
    sel = np.flatnonzero((om >= lo) & (om <= hi))
    if sel.size == 0 or bonds_xh is None or len(bonds_xh) == 0:
        return w
    U3 = U.reshape(natoms, 3, U.shape[1])[:, :, sel]
    for a, b in np.asarray(bonds_xh, dtype=np.int32).reshape(-1, 2):
        r = pos[b] - pos[a]
        nrm = float(np.linalg.norm(r))
        if nrm < 1e-12:
            raise ValueError(f"atom_stretch_in_window: zero X–H {a},{b}")
        q = ((r / nrm)[:, None] * (U3[b] - U3[a])).sum(axis=0)
        s = float((q * q).sum())
        w[a] += s; w[b] += s
    return w

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


def subtype_shrinkage_rows(names, strength, parameter_scale):
    """Build reference-free rows penalizing variance among subtypes of one elemental interaction.

    For each elemental family, the rows implement strength*sum_i((k_i-mean(k)) / s)^2.
    The pairwise form leaves the data-determined family mean unconstrained, unlike a
    fixed prior centered on a separately fitted elemental parameter.
    """
    names = list(names)
    parameter_scale = np.asarray(parameter_scale, dtype=np.float64)
    if strength < 0.0:
        raise ValueError("subtype shrinkage strength must be non-negative")
    if parameter_scale.shape != (len(names),) or np.any(~np.isfinite(parameter_scale)) or np.any(parameter_scale <= 0.0):
        raise ValueError("parameter_scale must contain one finite positive scale per parameter")
    groups = {}
    for i, name in enumerate(names):
        groups.setdefault(parent_parameter_name(name), []).append(i)
    if strength == 0.0:
        return np.empty((0, len(names))), np.empty(0), groups
    rows = []
    for members in groups.values():
        if len(members) < 2:
            continue
        scale = np.mean(parameter_scale[members])
        fac = np.sqrt(strength / len(members)) / scale
        for ia, i in enumerate(members[:-1]):
            for j in members[ia+1:]:
                row = np.zeros(len(names))
                row[i] = fac
                row[j] = -fac
                rows.append(row)
    R = np.asarray(rows) if rows else np.empty((0, len(names)))
    return R, np.zeros(R.shape[0]), groups


def family_mean_prior_rows(names, parent_values, strength):
    """Constrain each subtype-family mean, not individual subtype values, to its elemental parent."""
    names = list(names)
    if strength < 0.0:
        raise ValueError("family mean prior strength must be non-negative")
    _, _, groups = subtype_shrinkage_rows(names, 0.0, np.ones(len(names)))
    if strength == 0.0:
        return np.empty((0, len(names))), np.empty(0)
    rows = []
    targets = []
    for parent, members in groups.items():
        if len(members) < 2:
            continue
        if parent not in parent_values:
            raise KeyError(f"missing elemental parent value {parent}")
        value = float(parent_values[parent])
        scale = max(abs(value), 0.1)
        fac = np.sqrt(strength) / scale
        row = np.zeros(len(names))
        row[members] = fac / len(members)
        rows.append(row)
        targets.append(fac * value)
    R = np.asarray(rows) if rows else np.empty((0, len(names)))
    return R, np.asarray(targets)

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
    for (j0, k0) in bond_pairs:
        j, k = sorted((j0, k0))
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

def add_linear_hessian(H, k, A_extra):
    """Add precomputed linear-in-parameter Hessian contributions in-place."""
    for p, A in A_extra.items():
        H += k[p] * A

def add_dihedral_hessian(H, k, dihedral_A):
    """Compatibility alias for precomputed torsion Hessian contributions."""
    add_linear_hessian(H, k, dihedral_A)

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


def build_cross_param_maps(all_angles, all_symbols, offset, all_elements=None, angle_central_only=False, stretch_stretch=False, stretch_bend=False):
    """Assign optional valence cross parameters by the same chemical type as their angle."""
    all_elements = all_symbols if all_elements is None else all_elements
    keys = []
    for angles, symbols, elements in zip(all_angles, all_symbols, all_elements):
        for i, j, k, _ in angles:
            key = angle_type_key(symbols, i, j, k, elements=elements, central_only=angle_central_only)
            if key not in keys:
                keys.append(key)
    ss_map = {key: offset + i for i, key in enumerate(keys)} if stretch_stretch else {}
    sb_offset = offset + len(ss_map)
    sb_map = {key: sb_offset + i for i, key in enumerate(keys)} if stretch_bend else {}
    return ss_map, sb_map, offset + len(ss_map) + len(sb_map)


def compute_cross_sensitivity(positions, bonds, angles, symbols, stretch_stretch_map=None, stretch_bend_map=None, elements=None, angle_central_only=False):
    """Analytic Hessian sensitivities for optional symmetric valence cross terms.

    At the relaxed geometry q_a=q_b=dtheta=0, the exact harmonic contributions are
    d2(q_a*q_b) = g_a g_b^T + g_b g_a^T and
    d2(((q_a+q_b)/sqrt(2))*dtheta) = g_s g_theta^T + g_theta g_s^T.
    Bond coordinates q=dr/r0 are dimensionless, hence both fitted coefficients are
    signed energies in eV.  This routine intentionally rejects no geometry: callers
    must use local equilibrium coordinates when these sensitivities are fitted.
    """
    stretch_stretch_map = {} if stretch_stretch_map is None else stretch_stretch_map
    stretch_bend_map = {} if stretch_bend_map is None else stretch_bend_map
    if not stretch_stretch_map and not stretch_bend_map:
        return {}
    positions = np.asarray(positions, dtype=np.float64)
    elements = symbols if elements is None else elements
    from pyBall.FFfit import build_wilson_matrix, dimensionless_wilson_scale
    B, labels = build_wilson_matrix(positions, bonds, angles)
    scale = dimensionless_wilson_scale(positions, labels)
    pair_to_bond = {tuple(sorted((i, j))): ib for ib, (i, j, _) in enumerate(bonds)}
    n3 = positions.size
    indices = set(stretch_stretch_map.values()) | set(stretch_bend_map.values())
    A = {p: np.zeros((n3, n3)) for p in indices}
    for ia, (i, j, k, _) in enumerate(angles):
        ib1 = pair_to_bond.get(tuple(sorted((i, j))))
        ib2 = pair_to_bond.get(tuple(sorted((j, k))))
        if ib1 is None or ib2 is None:
            raise RuntimeError(f"angle {(i, j, k)} is missing one of its two bond coordinates")
        key = angle_type_key(symbols, i, j, k, elements=elements, central_only=angle_central_only)
        g1 = scale[ib1] * B[ib1]
        g2 = scale[ib2] * B[ib2]
        gt = B[len(bonds) + ia]
        if key in stretch_stretch_map:
            A[stretch_stretch_map[key]] += np.outer(g1, g2) + np.outer(g2, g1)
        if key in stretch_bend_map:
            gs = (g1 + g2) / np.sqrt(2.0)
            A[stretch_bend_map[key]] += np.outer(gs, gt) + np.outer(gt, gs)
    return A

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
    freq_floor_cm1: modes with |freq| < floor are skipped (translations/rotations).
    FFfit-only: this path uses sqrt(max(λ,0)). Harmonic *spectra* must use signed_frequencies_cm1 + assert_harmonic_spectrum_at_minimum.

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

def signed_frequencies_cm1(H, masses):
    """All 3N frequencies (cm^-1); imaginary eigenvalues returned as negative. Do not hide them with sqrt(max(λ,0))."""
    conv = 521.5
    masses = np.asarray(masses, dtype=np.float64).reshape(-1)
    inv_sqrt_m = 1.0 / np.sqrt(np.repeat(masses, 3))
    D = (inv_sqrt_m[:, None] * inv_sqrt_m[None, :]) * np.asarray(H, dtype=np.float64)
    lam = np.linalg.eigvalsh(D)
    return np.sign(lam) * conv * np.sqrt(np.abs(lam))

def as_signed_wavenumbers_cm1(omega_cm):
    """Real signed cm⁻¹, or complex (PySCF/DFTB: unstable modes stored as i|ν| → negative)."""
    w = np.asarray(omega_cm)
    if np.iscomplexobj(w):
        re = np.real(w); im = np.imag(w)
        return np.where(np.abs(im) > 1e-8, -np.abs(im), re).astype(np.float64)
    return np.asarray(w, dtype=np.float64).reshape(-1)

def assert_harmonic_spectrum_at_minimum(omega_cm, ctx="", imag_cut_cm1=10.0, n_rigid_min=6):
    """Fail loud if this is not a harmonic spectrum at that energy's own minimum.

    |ν|>imag_cut imaginary = negative curvature (off-minimum). Do not clamp, drop, or plot.
    Fewer than n_rigid_min modes below imag_cut = leftover torque (typical: MMFF Hessian at DFTB q).
    """
    w = as_signed_wavenumbers_cm1(omega_cm)
    if w.size == 0:
        raise RuntimeError(f"{ctx}empty frequency array")
    if not np.all(np.isfinite(w)):
        raise RuntimeError(f"{ctx}non-finite frequencies ({int(np.sum(~np.isfinite(w)))} bad)")
    cut = float(imag_cut_cm1)
    n_imag = int(np.sum(w < -cut))
    n_rigid = int(np.sum(np.abs(w) < cut))
    if n_imag:
        worst = np.sort(w[w < -cut])[:8]
        raise RuntimeError(
            f"{ctx}{n_imag} imaginary modes with |ν|>{cut:.0f} cm⁻¹ (most negative {float(w.min()):.1f} cm⁻¹; "
            f"sample {np.array2string(worst, precision=1)}). "
            "Hessian is not at a stationary point of this forcefield. "
            "Do not sqrt(max(λ,0)), do not drop modes, do not plot a spectrum. "
            "Relax with the SAME energy (same scales, same switches), check |F|max, then recompute. "
            "See doc/Topics/FTIR_Nanocrystals/Hessian_at_own_minimum.md"
        )
    if n_rigid < int(n_rigid_min):
        raise RuntimeError(
            f"{ctx}only {n_rigid} modes with |ν|<{cut:.0f} cm⁻¹ (need ≥{int(n_rigid_min)} rigid-body). "
            "Missing rotations/translations mean leftover forces/torques — typically MMFF Hessian at a DFTB/DFT geometry. "
            "Not a spectrum. Relax this FF first. See doc/Topics/FTIR_Nanocrystals/Hessian_at_own_minimum.md"
        )
    return {"n_imag": 0, "n_rigid": n_rigid, "nu_min": float(w.min()), "nu_max": float(w.max()), "n_modes": int(w.size)}

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

def param_names_from_maps(bmap, amap, b3map=None, dmap=None):
    """Ordered parameter names matching global indices from build_global_param_map."""
    b3map = {} if b3map is None else b3map
    dmap = {} if dmap is None else dmap
    npar = len(bmap) + len(b3map) + len(amap) + len(dmap)
    names = [""] * npar
    for key, t in bmap.items():
        names[t] = "bond:" + "-".join(key)
    for key, t in b3map.items():
        names[t] = "1-4:" + "-".join(key)
    for key, t in amap.items():
        names[t] = "angle:" + "-".join(key)
    for key, t in dmap.items():
        names[t] = "torsion:" + "/".join(("-".join(key[0]), "-".join(key[1])))
    if any(n == "" for n in names):
        raise RuntimeError(f"param name map has holes: {names}")
    return names

def freeze_indices_from_names(names, pin):
    """Map {param_name: value} → {index: value}. Missing names fail loud."""
    idx = {n: i for i, n in enumerate(names)}
    freeze = {}
    missing = [n for n in pin if n not in idx]
    if missing:
        raise KeyError(f"pin names {missing} not in this system {names}")
    for n, val in pin.items():
        v = float(val)
        if not np.isfinite(v):
            raise ValueError(f"non-finite pin for {n}")
        freeze[idx[n]] = v
    return freeze

def assert_relative_param_shift(old, new, names, tol=0.05, ctx=""):
    """Fail if any overlapping named k moved by more than tol (relative)."""
    bad = []
    for n in names:
        if n not in old or n not in new:
            raise KeyError(f"{ctx}param {n} missing in old={n in old} new={n in new}")
        a, b = float(old[n]), float(new[n])
        if abs(a) < 1e-12:
            if abs(b) > 1e-12:
                bad.append(f"{n}: 0 → {b}")
            continue
        rel = abs(b - a) / abs(a)
        if rel > tol:
            bad.append(f"{n}: {a:.6g} → {b:.6g} ({100.0*rel:.1f}%)")
    if bad:
        raise RuntimeError(f"{ctx}parameters moved by more than {100.0*tol:.1f}%:\n  " + "\n  ".join(bad))

def stretch_rmse_nH(omega_a, omega_b, nH, floor=10.0):
    """RMSE of the nH highest frequencies above floor. Metric only — not an LS loss."""
    nH = int(nH)
    a = np.sort(np.asarray(omega_a, dtype=np.float64))
    b = np.sort(np.asarray(omega_b, dtype=np.float64))
    a = a[a > float(floor)]; b = b[b > float(floor)]
    if nH < 1 or a.size < nH or b.size < nH:
        raise ValueError(f"stretch_rmse_nH: nH={nH} a={a.size} b={b.size}")
    d = a[-nH:] - b[-nH:]
    return float(np.sqrt(np.mean(d * d))), float(a[-nH:].mean()), float(b[-nH:].mean())

def fit_local_kab_system(sys, freeze_names=None, prior_by_name=None, regularization=2e-3,
                         local_weight=1.0, wilson_kab_weight=1.0, wilson_rows="hydride", bond_cutoff=2.5):
    """One-system local Cartesian + raw Wilson K_ab fit. No mode projection, no full-Wilson Q."""
    from pyBall import FFfit as FFfit_cpp
    freeze_names = {} if freeze_names is None else dict(freeze_names)
    prior_by_name = {} if prior_by_name is None else dict(prior_by_name)
    symbols = list(sys["symbols"])
    positions = np.asarray(sys["positions"], dtype=np.float64)
    H_ref = np.asarray(sys["H_ref"], dtype=np.float64)
    masses = np.asarray(sys["masses"], dtype=np.float64)
    bonds = sys.get("bonds")
    angles = sys.get("angles")
    if bonds is None or angles is None:
        bonds, angles, _ = build_topology(symbols, positions, bond_cutoff)
    bonds3 = sys.get("bonds3") or []
    bmap, b3map, amap, dmap, npar = build_global_param_map([bonds], [angles], [symbols], [bonds3], [[]], [symbols])
    if dmap:
        raise RuntimeError("fit_local_kab_system: torsions are not on the pin-ladder path")
    names = param_names_from_maps(bmap, amap, b3map, dmap)
    rec = {"name": sys.get("name", ""), "positions": positions, "symbols": symbols, "bonds": bonds, "angles": angles,
           "bonds3": bonds3, "atom_types": symbols, "angle_central_only": False}
    fitter = make_cpp_fitter(rec, (bmap, b3map, amap), npar)
    A = FFfit_cpp.collect_sensitivity_matrices(fitter)
    hybrid = [{"A": A, "H_ref": H_ref, "positions": positions, "masses": masses, "bonds": bonds, "angles": angles, "symbols": symbols}]
    prior = np.ones(npar)
    for t in bmap.values():
        prior[t] = 5.0
    for t in amap.values():
        prior[t] = 1.0
    for t in b3map.values():
        prior[t] = 0.1
    for n, v in prior_by_name.items():
        if n in names:
            prior[names.index(n)] = float(v)
    freeze = freeze_indices_from_names(names, freeze_names) if freeze_names else None
    k, diag = FFfit_cpp.fit_hybrid_hessian(
        hybrid, mode_weight=0.0, local_weight=local_weight, internal_weight=0.0,
        wilson_kab_weight=wilson_kab_weight, wilson_rows=wilson_rows, freeze=freeze,
        prior=prior, regularization=regularization, parameter_scale=np.maximum(np.abs(prior), 0.1),
        bounds=(0.0, np.inf))
    fitter.set_params(k)
    H_model = fitter.compute_model_hessian()
    nH = int(sum(s == "H" for s in symbols))
    om_ref = signed_frequencies_cm1(H_ref, masses)
    om_model = signed_frequencies_cm1(H_model, masses)
    rmse_st, mean_st_m, mean_st_r = stretch_rmse_nH(om_model, om_ref, nH)
    n_imag_m = int(np.sum(om_model < -10.0))
    k_by_name = {n: float(v) for n, v in zip(names, k)}
    return {"k": k, "names": names, "k_by_name": k_by_name, "diag": diag, "H_model": H_model, "H_ref": H_ref,
            "masses": masses, "nH": nH, "om_ref": om_ref, "om_model": om_model, "rmse_stretch": rmse_st,
            "mean_stretch_model": mean_st_m, "mean_stretch_ref": mean_st_r, "n_imag_model": n_imag_m,
            "bonds": bonds, "angles": angles, "bmap": bmap, "amap": amap, "fitter": fitter}

def run_si_two_phase_fit(plot_dir, reg_phase1=0.0, reg_phase2=1e-3, bond_cutoff=2.5,
                         include_13=False, si_subtypes=False, l1_root=None):
    """Two-phase joint Si fit:
    Phase 1: SiH4 alone → k_H-Si, k_H-Si-H (+ optionally k_H..H 1-3 spring).
    Phase 2: all L1 Si crystals (0-imag), hydride k frozen → fit k_Si-Si, k_Si-Si-Si, k_H-Si-Si.

    include_13:  add H..H 1-3 neighbour spring to Phase 1 to split T2/A1 bends.
    si_subtypes: split H-Si bond type by Si coordination (H-SiH3 vs H-SiH vs H-SiH2).
                 Phase 1 pins k(H-SiH3) only; Phase 2 fits k(H-SiH) and k(H-SiH2) freely.
    """
    import os, json
    from pyBall import FFfit as FFfit_cpp
    from pyBall.FFfit_plots import plot_spectrum
    plot_dir = Path(plot_dir); plot_dir.mkdir(parents=True, exist_ok=True)
    l1_root = PYSCF_VIB_L1_ROOT if l1_root is None else Path(l1_root)

    # ── Phase 1: SiH4 alone ───────────────────────────────────────────────────
    print("\n=== Phase 1: SiH4 alone (pin hydride params) ===")
    d_sih4 = PYSCF_SMALL_NC / "SiH4"
    data1 = load_pyscf_hessian_case(d_sih4)
    bonds1, angles1, _ = build_topology(data1["symbols"], data1["positions"], bond_cutoff)
    if include_13:
        # Add geminal (1-3) H..H springs: pairs connected through common Si.
        from pyBall.FFfit_utils import shortest_path_distances
        bond_pairs1 = [(b[0], b[1]) for b in bonds1]
        dist1 = shortest_path_distances(bond_pairs1, len(data1["symbols"]))
        import numpy as _np
        bonds3_1 = [(i, j, float(_np.linalg.norm(data1["positions"][j] - data1["positions"][i])))
                    for i in range(len(data1["symbols"])) for j in range(i+1, len(data1["symbols"]))
                    if dist1[i, j] == 2 and data1["symbols"][i] == "H" and data1["symbols"][j] == "H"]
        print(f"  include_13: added {len(bonds3_1)} geminal H..H springs at ~{bonds3_1[0][2]:.3f} A" if bonds3_1 else "  include_13: no geminal H..H pairs found!")
    else:
        bonds3_1 = []
    atom_types1 = data1["symbols"]
    bmap1, b3map1, amap1, dmap1, npar1 = build_global_param_map(
        [bonds1], [angles1], [atom_types1], [bonds3_1], [[]], [atom_types1], False)
    print(f"  Phase1 params: bond={len(bmap1)} angle={len(amap1)} 1-3={len(b3map1)} total={npar1}")
    for k, t in {**bmap1, **amap1, **b3map1}.items(): print(f"    {'-'.join(k) if isinstance(k,tuple) else k} -> idx {t}")
    f1 = make_cpp_fitter({"bonds": bonds1, "angles": angles1, "bonds3": bonds3_1,
                           "atom_types": atom_types1, "symbols": atom_types1, "positions": data1["positions"]},
                          (bmap1, b3map1, amap1), npar1)
    A1 = FFfit_cpp.collect_sensitivity_matrices(f1, extra={})
    prior1 = np.ones(npar1)
    for t in bmap1.values(): prior1[t] = 5.0
    for t in amap1.values(): prior1[t] = 1.0
    for t in b3map1.values(): prior1[t] = 0.5
    hs1 = [{"A": A1, "H_ref": data1["H_ref"], "positions": data1["positions"],
             "masses": data1["masses"], "bonds": bonds1, "angles": angles1, "symbols": atom_types1}]
    k1, diag1 = FFfit_cpp.fit_hybrid_hessian(hs1, mode_weight=1.0, local_weight=0.0, internal_weight=0.0,
        prior=prior1, regularization=reg_phase1, parameter_scale=prior1, bounds=(0.0, np.inf))
    print(f"  solver={diag1['solver']} residual={diag1['relative_residual']:.4e} cond={diag1['condition']:.4e}")
    hydride_names = set()
    hydride_freeze = {}
    for key, t in bmap1.items():
        nm = "bond:" + "-".join(key)
        print(f"  {nm} = {k1[t]:.6f} eV/A^2")
        hydride_names.add(nm); hydride_freeze[nm] = float(k1[t])
    for key, t in amap1.items():
        nm = "angle:" + "-".join(key)
        print(f"  {nm} = {k1[t]:.6f} eV/rad^2")
        hydride_names.add(nm); hydride_freeze[nm] = float(k1[t])
    for key, t in b3map1.items():
        nm = "1-3:" + "-".join(key)
        print(f"  {nm} = {k1[t]:.6f} eV/A^2")
        hydride_names.add(nm); hydride_freeze[nm] = float(k1[t])
    f1.set_params(k1)
    om_ref1 = signed_frequencies_cm1(data1["H_ref"], data1["masses"])
    om_mod1 = signed_frequencies_cm1(f1.compute_model_hessian(), data1["masses"])
    rmse1, mae1, _ = one_to_one_frequency_metrics(om_ref1[om_ref1>10], om_mod1[om_mod1>10])
    print(f"  SiH4  phase1: freq RMSE={rmse1:.1f} cm-1  MAE={mae1:.1f} cm-1")
    out1 = str(plot_dir / "phase1_SiH4_spectrum.png")
    plot_spectrum(om_ref1[om_ref1>10], om_mod1[om_mod1>10], "Phase1_SiH4", outdir=str(plot_dir))

    # ── Phase 2: all L1 Si crystals, hydride frozen ────────────────────────────
    print("\n=== Phase 2: L1 Si crystals (hydride frozen) ===")
    crystal_dirs = []
    for case in sorted(list_pyscf_vib_cases(l1_root)):
        if "Si" not in case: continue
        try: d = resolve_pyscf_vib_case_dir(l1_root, case)
        except FileNotFoundError: continue
        nim, _ = n_imag_from_omega_cm(np.load(d / "frequencies_cm1.npy"))
        if nim == 0: crystal_dirs.append(d)
    print(f"  {len(crystal_dirs)} L1 systems: {[d.name for d in crystal_dirs]}")

    crystal_systems = []
    for d in crystal_dirs:
        data = load_pyscf_hessian_case(d)
        bonds, angles, bonds3 = build_topology(data["symbols"], data["positions"], bond_cutoff, third_bonds=include_13)
        atom_types = assign_si_environment_types(data["symbols"], bonds, enabled=si_subtypes)
        crystal_systems.append({
            "name": d.name, "symbols": data["symbols"], "atom_types": atom_types,
            "positions": data["positions"], "masses": data["masses"],
            "H_ref": data["H_ref"], "bonds": bonds, "angles": angles, "bonds3": bonds3,
            "dihedrals": [], "angle_central_only": False, "data": data,
        })
    all_bonds2  = [s["bonds"]      for s in crystal_systems]
    all_angles2 = [s["angles"]     for s in crystal_systems]
    all_bonds3_2= [s["bonds3"]     for s in crystal_systems]
    all_syms2   = [s["atom_types"] for s in crystal_systems]
    all_elems2  = [s["symbols"]    for s in crystal_systems]
    bmap2, b3map2, amap2, dmap2, npar2 = build_global_param_map(
        all_bonds2, all_angles2, all_syms2, all_bonds3_2, [[]]*len(crystal_systems), all_elems2, False)
    print(f"  Phase2 params: bond={len(bmap2)} angle={len(amap2)} total={npar2}")
    for k_, t in bmap2.items(): print(f"    bond:{'-'.join(k_)} -> idx {t}")
    for k_, t in amap2.items(): print(f"    angle:{'-'.join(k_)} -> idx {t}")

    # Build freeze dict: indices of hydride params that also appear in phase2 global map
    names2 = [""] * npar2
    for key, t in bmap2.items():  names2[t] = "bond:"  + "-".join(key)
    for key, t in amap2.items():  names2[t] = "angle:" + "-".join(key)
    for key, t in b3map2.items(): names2[t] = "1-3:"   + "-".join(key)
    freeze2 = {t: hydride_freeze[names2[t]] for t in range(npar2) if names2[t] in hydride_freeze}
    print(f"  Frozen: {[names2[t] for t in freeze2]}")

    fitters2 = [make_cpp_fitter(sys, (bmap2, b3map2, amap2), npar2) for sys in crystal_systems]
    prior2 = np.ones(npar2)
    for t in bmap2.values(): prior2[t] = 5.0
    for t in amap2.values(): prior2[t] = 1.0
    hs2 = []
    for f, sys in zip(fitters2, crystal_systems):
        A = FFfit_cpp.collect_sensitivity_matrices(f, extra={})
        hs2.append({"A": A, "H_ref": sys["H_ref"], "positions": sys["positions"],
                    "masses": sys["masses"], "bonds": sys["bonds"], "angles": sys["angles"],
                    "symbols": sys["symbols"]})
    k2, diag2 = FFfit_cpp.fit_hybrid_hessian(hs2, mode_weight=1.0, local_weight=0.0, internal_weight=0.0,
        prior=prior2, regularization=reg_phase2, parameter_scale=prior2, bounds=(0.0, np.inf), freeze=freeze2)
    print(f"  solver={diag2['solver']} residual={diag2['relative_residual']:.4e} cond={diag2['condition']:.4e}")
    print("\n=== Phase 2 fitted parameters ===")
    for key, t in bmap2.items():  print(f"  k_bond[{'-'.join(key)}] = {k2[t]:.6f} eV/A^2  {'(FROZEN)' if t in freeze2 else ''}")
    for key, t in amap2.items():  print(f"  k_angle[{'-'.join(key)}] = {k2[t]:.6f} eV/rad^2  {'(FROZEN)' if t in freeze2 else ''}")

    print("\n=== Phase 2 per-system frequency RMSE ===")
    results = []
    for f, sys in zip(fitters2, crystal_systems):
        f.set_params(k2)
        H_model = f.compute_model_hessian()
        om_ref = signed_frequencies_cm1(sys["H_ref"], sys["masses"])
        om_mod = signed_frequencies_cm1(H_model, sys["masses"])
        rmse, mae, _ = one_to_one_frequency_metrics(om_ref[om_ref>10], om_mod[om_mod>10])
        relF = 100.0 * np.linalg.norm(H_model - sys["H_ref"]) / np.linalg.norm(sys["H_ref"])
        print(f"  {sys['name']:35s}: RMSE={rmse:.1f} cm-1  MAE={mae:.1f} cm-1  relFrob={relF:.2f}%")
        results.append({"name": sys["name"], "rmse": rmse, "mae": mae, "relFrob": relF})
        plot_spectrum(om_ref[om_ref>10], om_mod[om_mod>10], sys["name"], outdir=str(plot_dir))
    # Save JSON summary
    summary = {"phase1": {nm: v for nm, v in hydride_freeze.items()},
               "phase2": {names2[t]: float(k2[t]) for t in range(npar2)},
               "per_system": results}
    out_json = plot_dir / "two_phase_si_summary.json"
    out_json.write_text(json.dumps(summary, indent=2))
    print(f"\n  Summary: {out_json}")
    return summary

def load_joint_si_systems(bond_cutoff=2.5, si_subtypes=False, l1_root=None):
    """Load SiH4 (PBE/def2svp) + all 0-imag tight-Si L1 (PBE/ccECP-cc-pVDZ) for joint fitting.

    Both use PBE; basis differs. Caveat printed. Returns list of system dicts compatible
    with build_global_param_map / make_cpp_fitter / fit_hybrid_hessian.
    """
    l1_root = PYSCF_VIB_L1_ROOT if l1_root is None else Path(l1_root)
    dirs = [PYSCF_SMALL_NC / "SiH4"]
    for case in sorted(list_pyscf_vib_cases(l1_root)):
        if "Si" not in case:
            continue
        try:
            d = resolve_pyscf_vib_case_dir(l1_root, case)
        except FileNotFoundError:
            continue
        nim, _ = n_imag_from_omega_cm(np.load(d / "frequencies_cm1.npy"))
        if nim == 0:
            dirs.append(d)
    print(f"  joint Si: {len(dirs)} systems: {[d.name for d in dirs]}")
    print("  CAVEAT: SiH4 is PBE/def2svp; L1 crystals are PBE/ccECP-cc-pVDZ. Separate basis — pin-transfer only, not a mixed-XC LS.")
    systems = []
    for d in dirs:
        data = load_pyscf_hessian_case(d)
        bonds, angles, _ = build_topology(data["symbols"], data["positions"], bond_cutoff)
        atom_types = assign_si_environment_types(data["symbols"], bonds, enabled=si_subtypes)
        systems.append({
            "name": d.name, "dir": d,
            "symbols": data["symbols"], "atom_types": atom_types,
            "positions": data["positions"], "masses": data["masses"],
            "H_ref": data["H_ref"], "bonds": bonds, "angles": angles, "bonds3": [],
            "dihedrals": [], "angle_central_only": si_subtypes, "data": data,
        })
        nH = int(sum(s == "H" for s in data["symbols"]))
        print(f"    {d.name}: N={data['natoms']} nH={nH} bonds={len(bonds)} angles={len(angles)}")
    return systems

def si_pbe_ladder_stages():
    """Si PBE pin-ladder. Stage B (Si2H6) is skipped — no PBE/def2svp Hessian in the bank."""
    return [
        {"stage": "A", "label": "SiH4", "xc": "PBE/def2svp",
         "dir": PYSCF_SMALL_NC / "SiH4", "note": "owns k_XH and k_HXH"},
        {"stage": "C", "label": "Si10H16", "xc": "PBE/def2svp",
         "dir": CCU_CAGE_PYSCF / "Si10H_pyscf_pbe_def2svp", "note": "unlock XX / XXX / XXH; pin hydride from A"},
        {"stage": "D", "label": "octahedron_Si", "xc": "PBE/ccECP-cc-pVDZ",
         "l1_case": "octahedron_Si", "note": "pin hydride from A; basis transfer from def2svp"},
        {"stage": "E", "label": "cube_Si", "xc": "PBE/ccECP-cc-pVDZ",
         "l1_case": "cube_Si", "note": "same pin as D; compare k_XX to D"},
    ]

def run_si_pbe_pin_ladder(plot_dir, regularization=2e-3, pin_tol=0.05, bond_cutoff=2.5, l1_root=None):
    """Si PBE A→C→D→E: freeze hydride k from SiH4; do not mix Hessians in one LS; no rigid-P."""
    from pyBall.FFfit_plots import plot_spectrum
    l1_root = PYSCF_VIB_L1_ROOT if l1_root is None else Path(l1_root)
    plot_dir = Path(plot_dir)
    plot_dir.mkdir(parents=True, exist_ok=True)
    pin = {}
    k_history = {}
    stages_out = []
    print("=== Si PBE pin-ladder (local Hessian + raw Wilson K_ab; mode_weight=0, internal_weight=0) ===")
    print("  Stage B (Si2H6) skipped: no PBE/def2svp Hessian in the bank.")
    print("  Do not cite these frozen-q frequencies as FTIR. Own-min FIRE is a later gate.")
    for spec in si_pbe_ladder_stages():
        if "l1_case" in spec:
            case_dir = resolve_pyscf_vib_case_dir(l1_root, spec["l1_case"])
        else:
            case_dir = Path(spec["dir"])
        data = load_pyscf_hessian_case(case_dir)
        if "freqs" in data:
            n_imag_ref, wmin = n_imag_from_omega_cm(data["freqs"])
        else:
            om = signed_frequencies_cm1(data["H_ref"], data["masses"])
            n_imag_ref, wmin = int(np.sum(om < -10.0)), float(om.min())
        print(f"\n--- Stage {spec['stage']} {spec['label']}  {spec['xc']}  N={data['natoms']}  n_imag_ref={n_imag_ref}  νmin={wmin:.2f} ---")
        print(f"  dir={case_dir}")
        print(f"  {spec['note']}")
        if n_imag_ref:
            raise RuntimeError(f"{spec['label']}: n_imag={n_imag_ref} — not a fit target (saddle).")
        if spec["stage"] in ("D", "E"):
            print("  CAVEAT: L1 basis is ccECP-cc-pVDZ; hydride k is a pin-transfer from PBE/def2svp, not a stacked LS mix.")
        rec = fit_local_kab_system({"name": spec["label"], "symbols": data["symbols"], "positions": data["positions"],
                                    "H_ref": data["H_ref"], "masses": data["masses"]},
                                   freeze_names=pin if pin else None, prior_by_name=k_history.get("C", {}).get("k_by_name") if spec["stage"] in ("D", "E") else None,
                                   regularization=regularization, bond_cutoff=bond_cutoff)
        for n, v in zip(rec["names"], rec["k"]):
            print(f"  {n} = {v:.6g}")
        print(f"  condition={rec['diag']['condition']:.3g} residual={rec['diag']['relative_residual']:.3g} n_imag_model(frozen-q)={rec['n_imag_model']}")
        print(f"  stretch RMSE (nH={rec['nH']} highest, metric): {rec['rmse_stretch']:.2f} cm-1  <st>_FF={rec['mean_stretch_model']:.1f} <st>_DFT={rec['mean_stretch_ref']:.1f}")
        if spec["stage"] == "A":
            pin = dict(rec["k_by_name"])
            print(f"  pin from A (frozen thereafter): {pin}")
        k_history[spec["stage"]] = rec
        om_ref = np.sort(rec["om_ref"][rec["om_ref"] > 10.0])
        om_model = np.sort(rec["om_model"][rec["om_model"] > 10.0])
        plot_spectrum(om_ref, om_model, f"ladder_{spec['stage']}_{spec['label']}", outdir=str(plot_dir), xmax=2500.0)
        stages_out.append({
            "stage": spec["stage"], "label": spec["label"], "xc": spec["xc"], "dir": str(case_dir),
            "natoms": int(data["natoms"]), "nH": rec["nH"], "n_imag_ref": int(n_imag_ref),
            "k": rec["k_by_name"], "condition": float(rec["diag"]["condition"]),
            "relative_residual": float(rec["diag"]["relative_residual"]),
            "rmse_stretch_nH": rec["rmse_stretch"], "mean_stretch_model": rec["mean_stretch_model"],
            "mean_stretch_ref": rec["mean_stretch_ref"], "n_imag_model_frozen_q": rec["n_imag_model"],
            "note": spec["note"],
        })
    d_k = k_history["D"]["k_by_name"]
    e_k = k_history["E"]["k_by_name"]
    overlap = sorted(set(d_k) & set(e_k) - set(pin))
    xx_names = [n for n in overlap if n.startswith("bond:") and "H" not in n]
    print(f"\n=== D vs E unfrozen overlap (k_XX fail-loud tol={100.0*pin_tol:.1f}%; angles diagnostic) ===")
    d_vs_e = {}
    for n in overlap:
        a, b = d_k[n], e_k[n]
        rel = abs(b - a) / abs(a) if abs(a) > 1e-12 else float("inf")
        d_vs_e[n] = {"D": a, "E": b, "rel": rel}
        print(f"  {n}: D={a:.6g}  E={b:.6g}  Δ={100.0*rel:.2f}%")
    angle_bad = [n for n in overlap if n.startswith("angle:") and d_vs_e[n]["rel"] > pin_tol]
    if angle_bad:
        print("  NOTE: unfrozen angles differ by more than tol (not a pin failure; cube vs octa see-saw): " + ", ".join(angle_bad))
    summary = {"protocol": "Si PBE pin-ladder A(SiH4)→C(Si10H16)→D(octa)→E(cube); hydride freeze; local+K_ab; no PDP",
               "pin": pin, "stages": stages_out, "pin_tol": pin_tol, "d_vs_e": d_vs_e}
    js = plot_dir / "ladder_si_pbe.json"
    js.write_text(json.dumps(summary, indent=2))
    print(f"  wrote {js}")
    if xx_names:
        assert_relative_param_shift(d_k, e_k, xx_names, tol=pin_tol, ctx="D vs E k_XX: ")
    return summary

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
