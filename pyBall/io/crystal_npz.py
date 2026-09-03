"""crystal_npz.py — mmap-friendly loaders for nanocrystal pipeline NPZ schema v1.2.

Contract: doc/Topics/FTIR_Nanocrystals/NPZ_Crystal_Schema.md. validate_topology_crystal_parity
checks natoms/Z/bond count only — topology pos may differ from relaxed coords by design.
is_viewer_listable_basename / find_nanocrystal_pipeline_stages implement viewer grid rules (§7).
"""
import os
import numpy as np
from pyBall import elements
from pyBall.MolGUI_common import colors_for_enames, sizes_for_enames

_Z_TO_SYM = {int(rec[elements.index_Z]): rec[elements.index_symbol] for rec in elements.ELEMENTS}

def enames_from_Z(Z):
    Z = np.asarray(Z, dtype=np.int32).reshape(-1)
    return [_Z_TO_SYM.get(int(z), 'X') for z in Z]

def colors_from_Z(Z, alpha=1.0):
    return colors_for_enames(enames_from_Z(Z), alpha=alpha)

def sizes_from_Z(Z, scale=1.0):
    return sizes_for_enames(enames_from_Z(Z), scale=scale)

def _np_load(path, mmap=True):
    return np.load(path, mmap_mode='r' if mmap else None)

def _require_keys(d, keys, ctx):
    missing = [k for k in keys if k not in d.files]
    if missing:
        raise KeyError(f"{ctx}: missing required keys {missing} in {d}")

def detect_npz_kind(path):
    with _np_load(path, mmap=True) as d:
        if 'group_bbox_min' in d.files and 'group_bbox_max' in d.files:
            return 'topology'
        if 'neigh_idx' in d.files:
            return 'topology'
        return 'crystal'

def topology_has_mmffl(topology):
    return topology is not None and 'neigh_idx' in topology and 'bond_k' in topology and 'bond_l0' in topology

def bonds_ij_to_segments(pos, bonds_ij):
    bonds_ij = np.asarray(bonds_ij, dtype=np.int32)
    if bonds_ij.size == 0:
        return np.zeros((0, 3), dtype=np.float32)
    if bonds_ij.ndim != 2 or bonds_ij.shape[1] != 2:
        raise ValueError(f"bonds_ij_to_segments: bonds_ij.shape={bonds_ij.shape} expected (Nb,2)")
    pos = np.asarray(pos, dtype=np.float32)
    segs = np.empty((bonds_ij.shape[0] * 2, 3), dtype=np.float32)
    segs[0::2] = pos[bonds_ij[:, 0]]
    segs[1::2] = pos[bonds_ij[:, 1]]
    return segs

def bonds_from_neigh_idx(neigh_idx, stick_class=None):
    neigh_idx = np.asarray(neigh_idx, dtype=np.int32)
    if neigh_idx.ndim != 2:
        raise ValueError(f"bonds_from_neigh_idx: neigh_idx.ndim={neigh_idx.ndim} expected 2")
    N, M = neigh_idx.shape
    i_idx = np.repeat(np.arange(N, dtype=np.int32), M)
    j_idx = neigh_idx.reshape(-1)
    valid = j_idx >= 0
    if stick_class is not None:
        sc = np.asarray(stick_class, dtype=np.int32).reshape(-1)
        valid &= (sc == 1)
    i_idx = i_idx[valid]; j_idx = j_idx[valid]
    keep = i_idx < j_idx
    if not np.any(keep):
        return np.zeros((0, 2), dtype=np.int32)
    return np.stack([i_idx[keep], j_idx[keep]], axis=1).astype(np.int32)

def bond_colors_by_stiffness(pos, bonds_ij, bond_k_flat, color_range=(0.0, 1.0)):
    """Map bond stiffness to blue (soft) → red (stiff) RGBA segment colors."""
    bonds_ij = np.asarray(bonds_ij, dtype=np.int32)
    k = np.asarray(bond_k_flat, dtype=np.float32).reshape(-1)
    if bonds_ij.shape[0] != k.shape[0]:
        raise ValueError(f"bond_colors_by_stiffness: bonds={bonds_ij.shape[0]} k={k.shape[0]}")
    vmin, vmax = float(np.min(k)), float(np.max(k))
    segs = bonds_ij_to_segments(pos, bonds_ij)
    colors = np.empty((bonds_ij.shape[0] * 2, 4), dtype=np.float32)
    for i, ki in enumerate(k):
        f = 0.5 if abs(vmax - vmin) < 1e-8 else (float(ki) - vmin) / (vmax - vmin)
        f = float(np.clip(f, color_range[0], color_range[1]))
        colors[2*i:2*i+2] = (f, 0.0, 1.0 - f, 0.85)
    return segs, colors

def extract_bond_k_for_pairs(neigh_idx, bond_k, bonds_ij, stick_class=None):
    """Return stiffness per bond pair from packed neigh_idx/bond_k slots."""
    neigh_idx = np.asarray(neigh_idx, dtype=np.int32)
    bond_k = np.asarray(bond_k, dtype=np.float32)
    bonds_ij = np.asarray(bonds_ij, dtype=np.int32)
    pair_to_k = {}
    N, M = neigh_idx.shape
    for i in range(N):
        for s in range(M):
            j = int(neigh_idx[i, s])
            if j < 0:
                continue
            if stick_class is not None and int(stick_class[i, s]) != 1:
                continue
            a, b = (i, j) if i < j else (j, i)
            pair_to_k[(a, b)] = float(bond_k[i, s])
    out = np.empty((bonds_ij.shape[0],), dtype=np.float32)
    for bi, (a, b) in enumerate(bonds_ij):
        out[bi] = pair_to_k.get((int(a), int(b)), np.nan)
    return out

def bboxes_to_atomscene(bmin, bmax):
    bmin = np.asarray(bmin, dtype=np.float32)
    bmax = np.asarray(bmax, dtype=np.float32)
    if bmin.shape != bmax.shape or bmin.ndim != 2 or bmin.shape[1] < 3:
        raise ValueError(f"bboxes_to_atomscene: bmin={bmin.shape} bmax={bmax.shape} expected (G,3+)")
    G = bmin.shape[0]
    bbmin = np.zeros((G, 4), dtype=np.float32); bbmax = np.zeros((G, 4), dtype=np.float32)
    bbmin[:, :3] = bmin[:, :3]; bbmax[:, :3] = bmax[:, :3]
    return bbmin, bbmax

def as_bonds_ij(bonds, n_atoms=None, ctx='bonds_ij'):
    """Require explicit (Nb,2) 0-based pairs. Distance-guessing drops 5-ring closers."""
    if bonds is None or len(bonds) == 0:
        raise ValueError(f"{ctx}: empty bonds — crystal NPZ/xyz must store explicit topology")
    arr = np.asarray(bonds, dtype=np.int32)
    if arr.ndim == 1:
        if arr.size % 2 != 0:
            raise ValueError(f"{ctx}: flat bonds length {arr.size} not even")
        arr = arr.reshape(-1, 2)
    if arr.ndim != 2 or arr.shape[1] != 2:
        raise ValueError(f"{ctx}: bonds.shape={arr.shape} expected (Nb,2)")
    if (arr < 0).any():
        raise ValueError(f"{ctx}: negative atom index in bonds")
    if n_atoms is not None and (arr >= int(n_atoms)).any():
        raise ValueError(f"{ctx}: bond index >= n_atoms={n_atoms}")
    return arr


def map_atoms_by_position(pos_src, pos_dst, tol=1e-3):
    """Index map src→dst by coordinates. MMFF packs heavies-then-H; mol2 may be interleaved."""
    pos_src = np.asarray(pos_src, dtype=np.float64)
    pos_dst = np.asarray(pos_dst, dtype=np.float64)
    if pos_src.shape != pos_dst.shape or pos_src.ndim != 2 or pos_src.shape[1] != 3:
        raise ValueError(f"map_atoms_by_position: shapes {pos_src.shape} vs {pos_dst.shape}")
    n = pos_src.shape[0]
    mapping = np.empty(n, dtype=np.int32)
    used = np.zeros(n, dtype=bool)
    for i in range(n):
        d2 = np.sum((pos_dst - pos_src[i]) ** 2, axis=1)
        j = int(np.argmin(d2))
        if d2[j] > tol * tol:
            raise ValueError(f"map_atoms_by_position: src[{i}] min distance {np.sqrt(d2[j]):.4g} Å > {tol}")
        if used[j]:
            raise ValueError(f"map_atoms_by_position: dst[{j}] already matched")
        used[j] = True
        mapping[i] = j
    return mapping


def remap_bonds_ij(bonds_ij, mapping, n_atoms=None):
    bonds_ij = as_bonds_ij(bonds_ij, n_atoms=len(mapping) if n_atoms is None else n_atoms, ctx='remap_bonds_ij')
    mapping = np.asarray(mapping, dtype=np.int32).reshape(-1)
    out = np.stack([mapping[bonds_ij[:, 0]], mapping[bonds_ij[:, 1]]], axis=1).astype(np.int32)
    swap = out[:, 0] > out[:, 1]
    out[swap] = out[swap][:, ::-1]
    return as_bonds_ij(out, n_atoms=int(mapping.size), ctx='remap_bonds_ij.out')


def load_bonds_ij_from_mol2(path):
    from pyBall.atomicUtils import loadMol2
    apos, _atypes, _enames, _qs, bonds, *_rest = loadMol2(path)
    return as_bonds_ij(bonds, n_atoms=len(apos), ctx=f'load_bonds_ij_from_mol2({path})')


def load_bonds_ij_from_xyz(path):
    with open(path, 'r') as f:
        lines = f.read().splitlines()
    if len(lines) < 2:
        raise ValueError(f"load_bonds_ij_from_xyz({path}): file too short")
    hdr = lines[0].split()
    n = int(hdr[0])
    nb = int(hdr[1]) if len(hdr) > 1 else 0
    if n <= 0 or nb <= 0:
        raise ValueError(f"load_bonds_ij_from_xyz({path}): need 'nAtoms nBonds' header, got {lines[0]!r}")
    need = 2 + n + nb
    if len(lines) < need:
        raise ValueError(f"load_bonds_ij_from_xyz({path}): expected {need} lines, got {len(lines)}")
    pairs = []
    for k in range(nb):
        w = lines[2 + n + k].split()
        if len(w) < 2:
            raise ValueError(f"load_bonds_ij_from_xyz({path}): bond line {k} parse failed")
        a1, a2 = int(w[0]), int(w[1])
        pairs.append((a1 - 1, a2 - 1))
    return as_bonds_ij(pairs, n_atoms=n, ctx=f'load_bonds_ij_from_xyz({path})')


def load_bonds_ij_from_init(path):
    suf = os.path.splitext(path)[1].lower()
    if suf == '.mol2':
        return load_bonds_ij_from_mol2(path)
    if suf == '.xyz':
        return load_bonds_ij_from_xyz(path)
    raise ValueError(f"load_bonds_ij_from_init: unsupported '{path}' (use .mol2 or .xyz with explicit bonds)")


def rewrite_crystal_npz_bonds(npz_path, bonds_ij):
    """Keep pos/Z/relax metadata; set required bonds_ij. Used to backfill atlas 02_relaxed.npz."""
    src = np.load(npz_path, allow_pickle=False)
    data = {k: np.asarray(src[k]) for k in src.files}
    n = int(np.asarray(data['Z']).reshape(-1).shape[0])
    data['bonds_ij'] = as_bonds_ij(bonds_ij, n_atoms=n, ctx=f'rewrite_crystal_npz_bonds({npz_path})')
    np.savez(npz_path, **data)


def load_crystal_npz(path, mmap=True):
    d = _np_load(path, mmap=mmap)
    _require_keys(d, ['pos', 'Z'], 'load_crystal_npz')
    pos = np.asarray(d['pos'], dtype=np.float64)
    Z = np.asarray(d['Z'], dtype=np.int32).reshape(-1)
    if pos.shape != (Z.shape[0], 3):
        raise ValueError(f"load_crystal_npz: pos.shape={pos.shape} Z.shape={Z.shape} mismatch in {path}")
    bonds_ij = None
    if 'bonds_ij' in d.files:
        bonds_ij = as_bonds_ij(d['bonds_ij'], n_atoms=Z.shape[0], ctx=f'load_crystal_npz({path})')
    return {'pos': pos, 'Z': Z, 'bonds_ij': bonds_ij, 'path': path}

_TOPOLOGY_MMFFL_KEYS = ('neigh_idx', 'bond_k', 'bond_l0', 'stick_class', 'neigh_count', 'KLs', 'max_neighbors', 'n_bond', 'n_angle', 'n_dihedral')
_TOPOLOGY_SURFACE_KEYS = ('radius', 'icol', 'icolGroup', 'group_atoms', 'group_nAtoms', 'n_groups', 'group_cap', 'excl_icol', 'excl_count')
_TOPOLOGY_META_KEYS = ('schema_version', 'source_mol2', 'defects_json')


def count_topology_sticks(topology):
    """Count active MMFFL sticks (bond/angle/dihedral slots) in packed topology."""
    if not topology_has_mmffl(topology):
        raise ValueError('count_topology_sticks: topology missing neigh_idx/bond_k/bond_l0')
    neigh_idx = np.asarray(topology['neigh_idx'], dtype=np.int32)
    sc = topology.get('stick_class')
    nc = topology.get('neigh_count')
    N, M = neigh_idx.shape
    seen = set()
    n_sticks = 0
    for i in range(N):
        m = int(nc[i]) if nc is not None else M
        for s in range(m):
            j = int(neigh_idx[i, s])
            if j < 0:
                continue
            if sc is not None and int(np.asarray(sc)[i, s]) == 0:
                continue
            a, b = (i, j) if i < j else (j, i)
            key = (a, b)
            if key in seen:
                continue
            seen.add(key)
            n_sticks += 1
    return n_sticks


def validate_topology_crystal_parity(crystal, topology):
    """Fail loudly if crystal and 03_topology.npz disagree on atom count, Z, or bond count.

    Topology connectivity and spring params (neigh_idx, bond_k, bond_l0, …) are fixed at export.
    Relaxation changes coordinates only — topology pos in the NPZ is an init snapshot, not required
  to match 02_relaxed.npz.
    """
    Z = np.asarray(crystal['Z'], dtype=np.int32).reshape(-1)
    natoms = Z.shape[0]
    pos = np.asarray(crystal['pos'], dtype=np.float64)
    if pos.shape != (natoms, 3):
        raise ValueError(f"validate_topology_crystal_parity: crystal pos.shape={pos.shape} Z={natoms}")
    if 'Z' in topology:
        tZ = np.asarray(topology['Z'], dtype=np.int32).reshape(-1)
        if tZ.shape[0] != natoms:
            raise ValueError(f"validate_topology_crystal_parity: topology natoms={tZ.shape[0]} crystal natoms={natoms}")
        if not np.array_equal(tZ, Z):
            raise ValueError("validate_topology_crystal_parity: topology Z != crystal Z")
    if 'pos' in topology:
        tpos = np.asarray(topology['pos'], dtype=np.float64)
        if tpos.shape != pos.shape:
            raise ValueError(f"validate_topology_crystal_parity: topology pos.shape={tpos.shape} crystal pos.shape={pos.shape}")
    bonds_ij = crystal.get('bonds_ij')
    if bonds_ij is not None and len(bonds_ij) > 0 and topology_has_mmffl(topology):
        nb = int(np.asarray(bonds_ij).shape[0])
        n_bond = int(topology['n_bond']) if 'n_bond' in topology else None
        if n_bond is not None and n_bond != nb:
            raise ValueError(f"validate_topology_crystal_parity: n_bond={n_bond} != bonds_ij={nb}")
    return {'natoms': natoms, 'n_sticks': count_topology_sticks(topology) if topology_has_mmffl(topology) else None}


def load_topology_npz(path, mmap=True, full=False):
    """Load topology NPZ: viewer AABB overlay required; MMFFL slots optional (parity: TopologyNpz.h)."""
    d = _np_load(path, mmap=mmap)
    _require_keys(d, ['group_bbox_min', 'group_bbox_max'], 'load_topology_npz')
    out = {
        'group_bbox_min': np.asarray(d['group_bbox_min'], dtype=np.float64),
        'group_bbox_max': np.asarray(d['group_bbox_max'], dtype=np.float64),
        'path': path,
    }
    if 'icolGroup' in d.files:
        out['icolGroup'] = np.asarray(d['icolGroup'], dtype=np.int32).reshape(-1)
    if 'pos' in d.files:
        out['pos'] = np.asarray(d['pos'], dtype=np.float64)
    if 'Z' in d.files:
        out['Z'] = np.asarray(d['Z'], dtype=np.int32).reshape(-1)
    key_lists = (_TOPOLOGY_MMFFL_KEYS + _TOPOLOGY_SURFACE_KEYS) if full else ('neigh_idx', 'bond_k', 'bond_l0', 'stick_class', 'neigh_count', 'radius')
    for k in key_lists:
        if k in d.files:
            out[k] = np.asarray(d[k])
    if full:
        for k in _TOPOLOGY_META_KEYS:
            if k in d.files:
                out[k] = d[k]
    return out

def load_npy_crystal(dir_path):
    pos_path = os.path.join(dir_path, 'pos.npy')
    Z_path = os.path.join(dir_path, 'Z.npy')
    if not os.path.isfile(pos_path):
        raise FileNotFoundError(f"load_npy_crystal: missing {pos_path}")
    if not os.path.isfile(Z_path):
        raise FileNotFoundError(f"load_npy_crystal: missing {Z_path}")
    pos = np.load(pos_path, mmap_mode='r')
    Z = np.load(Z_path, mmap_mode='r')
    pos = np.asarray(pos, dtype=np.float64)
    Z = np.asarray(Z, dtype=np.int32).reshape(-1)
    bonds_ij = None
    bonds_path = os.path.join(dir_path, 'bonds_ij.npy')
    if os.path.isfile(bonds_path):
        bonds_ij = np.asarray(np.load(bonds_path, mmap_mode='r'), dtype=np.int32)
    return {'pos': pos, 'Z': Z, 'bonds_ij': bonds_ij, 'path': dir_path}

def infer_bonds_if_missing(pos, Z, bonds_ij):
    if bonds_ij is not None and len(bonds_ij) > 0:
        return np.asarray(bonds_ij, dtype=np.int32)
    from pyBall.AtomicSystem import AtomicSystem
    enames = enames_from_Z(Z)
    sys = AtomicSystem(apos=np.asarray(pos, dtype=np.float64), atypes=np.asarray(Z, dtype=np.int32), enames=enames, bPreinit=False)
    sys.findBonds()
    if sys.bonds is None or len(sys.bonds) == 0:
        return np.zeros((0, 2), dtype=np.int32)
    return np.asarray(sys.bonds, dtype=np.int32)

_PIPELINE_STAGE_FILES = {'crystal': '01_crystal.npz', 'relaxed': '02_relaxed.npz', 'topology': '03_topology.npz', 'hessian': '04_hessian.npz', 'spectrum': '05_spectrum.npz'}

def is_viewer_listable_basename(name):
    """NPZ listing rule — mirror NPZ_Crystal_Schema §7 / C++ filterNonViewerNpzFiles."""
    low = os.path.basename(name).lower()
    if not low.endswith('.npz'):
        return True
    if low.startswith('04_') or low.startswith('05_'):
        return False
    if 'hessian' in low or 'spectrum' in low:
        return False
    return True

def find_nanocrystal_pipeline_stages(dir_path):
    """Return existing canonical stage NPZ paths in a crystal bundle directory."""
    d = os.path.abspath(dir_path)
    out = {}
    for stage, fn in _PIPELINE_STAGE_FILES.items():
        p = os.path.join(d, fn)
        if os.path.isfile(p):
            out[stage] = p
    return out

def pipeline_dir_for_molecule_path(mol_path):
    """Directory containing molecule NPZ (or npy bundle sentinel path)."""
    p = os.path.normpath(mol_path.rstrip(os.sep))
    return p if os.path.isdir(p) else os.path.dirname(p)
