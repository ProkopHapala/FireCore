import os
import sys
import re
import json
import numpy as np

_THIS_DIR = os.path.abspath(os.path.dirname(__file__))
_SHARED_DIR = os.path.abspath(os.path.join(_THIS_DIR, '..', 'shared'))
for _p in (_THIS_DIR, _SHARED_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from XPTB_utils import pack_molecules_contiguous, make_h2o_geometry, masses_from_elems, perturb_state

# ------------------------------------------------------------------
# Molecule Builders
# ------------------------------------------------------------------

def make_ch3oh_geometry():
    """Return (pos, bonds, nnode, elems) for methanol CH3OH.

    Atom order: 0=C, 1=O, 2=H(methyl), 3=H(methyl), 4=H(methyl), 5=H(hydroxyl)
    Nodes: C(0), O(1)  |  Caps: H(2,3,4,5)
    Bonds: C-O, C-Hx3, O-H
    """
    # C-O bond length ~1.43 A, C-H ~1.09 A, O-H ~0.96 A
    # Tetrahedral angle ~109.5 deg
    # Place C at origin, O along +x, methyl H's in tetrahedral arrangement
    c = np.array([0.0, 0.0, 0.0], dtype=np.float32)
    o = np.array([1.43, 0.0, 0.0], dtype=np.float32)

    tet = np.deg2rad(109.5)
    # One H on C opposite to O (along -x, slightly)
    h_c1 = np.array([-1.09, 0.0, 0.0], dtype=np.float32)
    # Two other H's on C in yz plane, tetrahedral angle from C-O bond
    # C-O is along +x. Tetrahedral from x-axis: spread into yz
    cos_tet = np.cos(tet)
    sin_tet = np.sin(tet)
    h_c2 = np.array([1.09 * cos_tet, 1.09 * sin_tet * np.cos(0.0), 1.09 * sin_tet * np.sin(0.0)], dtype=np.float32)
    h_c3 = np.array([1.09 * cos_tet, 1.09 * sin_tet * np.cos(np.deg2rad(120.0)), 1.09 * sin_tet * np.sin(np.deg2rad(120.0))], dtype=np.float32)

    # H on O: O-H bond, angle ~108.5 (bent), place in xy plane, away from C
    oh_ang = np.deg2rad(108.5)
    h_o = o + np.array([0.96 * np.cos(np.pi - oh_ang), 0.96 * np.sin(np.pi - oh_ang), 0.0], dtype=np.float32)

    pos = np.stack([c, o, h_c1, h_c2, h_c3, h_o], axis=0)
    bonds = [(0, 1), (0, 2), (0, 3), (0, 4), (1, 5)]
    nnode = 2  # C and O are nodes
    elems = ['C', 'O', 'H', 'H', 'H', 'H']
    return pos, bonds, nnode, elems


def make_tri3_geometry(a=1.5):
    """Return (pos, bonds, nnode, elems) for a 3-node triangle ring (no caps).

    This is a minimal nodes-only rigid body (no capping atoms) to debug momentum
    conservation without bkSlots recoil gathering.

    Atom order: 0,1,2 all nodes.
    Bonds: (0-1),(1-2),(2-0)
    """
    a = float(a)
    p0 = np.array([0.0, 0.0, 0.0], dtype=np.float32)
    p1 = np.array([a, 0.0, 0.0], dtype=np.float32)
    p2 = np.array([0.5 * a, 0.8660254 * a, 0.0], dtype=np.float32)
    pos = np.stack([p0, p1, p2], axis=0)
    bonds = [(0, 1), (1, 2), (2, 0)]
    nnode = 3
    elems = ['C', 'C', 'C']
    return pos, bonds, nnode, elems


def make_ch2nh_geometry():
    """Return (pos, bonds, nnode, elems) for methanimine CH2NH (formaldimine).

    Atom order: 0=C, 1=N, 2=H(onC), 3=H(onC), 4=H(onN)
    Nodes: C(0), N(1)  |  Caps: H(2,3,4)
    Bonds: C=N, C-Hx2, N-H
    """
    c = np.array([0.0, 0.0, 0.0], dtype=np.float32)
    n = np.array([1.27, 0.0, 0.0], dtype=np.float32)  # C=N ~1.27 A

    # H's on C: trigonal planar ~120 deg
    h_c1 = np.array([-1.09 * np.cos(np.deg2rad(60.0)), 1.09 * np.sin(np.deg2rad(60.0)), 0.0], dtype=np.float32)
    h_c2 = np.array([-1.09 * np.cos(np.deg2rad(60.0)), -1.09 * np.sin(np.deg2rad(60.0)), 0.0], dtype=np.float32)

    # H on N: ~107 deg from C=N bond
    nh_ang = np.deg2rad(107.0)
    h_n = n + np.array([1.01 * np.cos(np.pi - nh_ang), 1.01 * np.sin(np.pi - nh_ang), 0.0], dtype=np.float32)

    pos = np.stack([c, n, h_c1, h_c2, h_n], axis=0)
    bonds = [(0, 1), (0, 2), (0, 3), (1, 4)]
    nnode = 2
    elems = ['C', 'N', 'H', 'H', 'H']
    return pos, bonds, nnode, elems


def make_test_system(name):
    """Factory for test molecules. Returns dict with keys: elems, pos, bonds, nnode.

    Supported names:
        'h2o'      -> single water (3 atoms, 1 node)
        'ch3oh'    -> methanol (6 atoms, 2 nodes, asymmetric)
        'ch2nh'    -> methanimine (5 atoms, 2 nodes, asymmetric)
    """
    name = str(name).lower().strip()
    if name == 'h2o':
        pos, _ = make_h2o_geometry(add_angle=False)
        elems = ['O', 'H', 'H']
        bonds = [(0, 1), (0, 2)]
        nnode = 1
    elif name == 'ch3oh':
        pos, bonds, nnode, elems = make_ch3oh_geometry()
    elif name == 'ch2nh':
        pos, bonds, nnode, elems = make_ch2nh_geometry()
    elif name in ('tri3', 'tri3_nodes', 'tri3_ring'):
        pos, bonds, nnode, elems = make_tri3_geometry()
    else:
        raise ValueError(f"make_test_system: unknown molecule '{name}'")
    return {'elems': elems, 'pos': pos, 'bonds': bonds, 'nnode': nnode}


# ------------------------------------------------------------------
# Test System Assembly (packed for RRsp3)
# ------------------------------------------------------------------

def build_packed_system(molecule_names, shift=None, group_size=64, rng=None):
    """Build a packed multi-molecule system for RRsp3 testing.

    Args:
        molecule_names: list of molecule name strings (e.g., ['h2o','h2o'])
        shift: list of (x,y,z) shifts per molecule, or single shift for all
        group_size: must be 64 for RRsp3
        rng: optional np.random.Generator for perturbation

    Returns:
        dict with keys: sim, natoms, pos, quat, invm, elems_all, is_pad, real,
                        neighs, excl1, excl2, nnode_per_group, bkSlots,
                        port_local_atoms, K_atoms, rad, fixmask,
                        plus metadata: mols, packed, group_size
    """
    from RRsp3 import build_neighs_bk_from_bonds, make_bk_slots_clustered, make_exclusions_1st_2nd

    mols = [make_test_system(n) for n in molecule_names]
    nmol = len(mols)

    if shift is None:
        shift = [np.array([float(i) * 3.0, 0.0, 0.0], dtype=np.float32) for i in range(nmol)]
    elif isinstance(shift, (list, tuple)) and len(shift) == 3 and isinstance(shift[0], (int, float)):
        # single 3-vector, broadcast to all molecules with offset
        base = np.array(shift, dtype=np.float32)
        shift = [base * float(i) for i in range(nmol)]
    else:
        shift = [np.array(s, dtype=np.float32) for s in shift]

    for i, mol in enumerate(mols):
        mol['pos'] = mol['pos'] + shift[i]

    packed = pack_molecules_contiguous(mols, group_size=group_size, nodes_first=True, pad_to_group=True)
    natoms = int(packed['natoms_total'])
    pos = packed['pos']
    elems_all = packed['elems']
    is_pad = packed['is_padding']
    real = ~is_pad

    # masses
    elems_real = [e for e in elems_all if e != 'X']
    m_real = masses_from_elems(elems_real)
    m = np.zeros((natoms,), dtype=np.float32)
    m[real] = m_real
    invm = np.zeros_like(m)
    invm[real] = 1.0 / m[real]

    # topology
    neighs, _ = build_neighs_bk_from_bonds(natoms, packed['bonds'], max_deg=4)
    excl1, excl2 = make_exclusions_1st_2nd(neighs)

    nnode_per_group = int(packed['nnode_group'][0])
    if not np.all(packed['nnode_group'] == nnode_per_group):
        raise RuntimeError(f"build_packed_system: variable nnode_per_group not supported: {packed['nnode_group']}")

    bkSlots = make_bk_slots_clustered(neighs, group_size=group_size, nnode_per_group=nnode_per_group, natoms=natoms)

    # ports
    port_local_atoms, K_atoms = make_ports_from_neighs(pos, neighs, K=200.0)

    # radius
    rad = np.zeros((natoms,), dtype=np.float32)
    rad[real] = 1.0

    # fixmask (all free initially)
    fixmask = np.zeros((natoms,), dtype=np.int32)
    fixmask[is_pad] |= (1 | 2 | 4)

    # quat
    quat = np.zeros((natoms, 4), dtype=np.float32)
    quat[:, 3] = 1.0

    return {
        'natoms': natoms,
        'pos': pos,
        'quat': quat,
        'invm': invm,
        'm': m,
        'elems_all': elems_all,
        'is_pad': is_pad,
        'real': real,
        'neighs': neighs,
        'excl1': excl1,
        'excl2': excl2,
        'nnode_per_group': nnode_per_group,
        'bkSlots': bkSlots,
        'port_local_atoms': port_local_atoms,
        'K_atoms': K_atoms,
        'rad': rad,
        'fixmask': fixmask,
        'mols': mols,
        'packed': packed,
        'group_size': group_size,
    }


# ------------------------------------------------------------------
# Port / Topology Utilities
# ------------------------------------------------------------------

def make_ports_from_neighs(pos, neighs, K=200.0):
    pos = np.asarray(pos, dtype=np.float32)
    neighs = np.asarray(neighs, dtype=np.int32)
    n = int(pos.shape[0])
    if pos.shape != (n, 3):
        raise ValueError(f"make_ports_from_neighs: pos.shape={pos.shape} expected ({n},3)")
    if neighs.shape != (n, 4):
        raise ValueError(f"make_ports_from_neighs: neighs.shape={neighs.shape} expected ({n},4)")

    port_local = np.zeros((n, 4, 4), dtype=np.float32)
    Kflat = np.zeros((n, 4), dtype=np.float32)
    for i in range(n):
        for k in range(4):
            j = int(neighs[i, k])
            if j < 0:
                continue
            port_local[i, k, :3] = pos[j] - pos[i]
            Kflat[i, k] = float(K)
    return port_local, Kflat


def reorder_nodes_first(elems, xyz, bonds, *, is_node=None):
    """Reorder atoms: nodes first, caps second. Returns (elems2, xyz2, bonds2, nnode)."""
    elems = list(elems)
    xyz = np.asarray(xyz, dtype=np.float32)
    n = int(len(elems))
    if xyz.shape != (n, 3):
        raise ValueError(f"reorder_nodes_first: xyz.shape={xyz.shape} expected ({n},3)")
    if is_node is None:
        is_node = np.array([e != 'H' for e in elems], dtype=bool)
    is_node = np.asarray(is_node, dtype=bool)
    if is_node.shape != (n,):
        raise ValueError(f"reorder_nodes_first: is_node.shape={is_node.shape} expected ({n},)")
    if not np.any(is_node):
        is_node[:] = True
    order = np.concatenate([np.nonzero(is_node)[0], np.nonzero(~is_node)[0]]).astype(np.int32)
    perm_inv = np.empty((n,), dtype=np.int32)
    perm_inv[order] = np.arange(n, dtype=np.int32)
    elems2 = [elems[i] for i in order]
    xyz2 = xyz[order, :].copy()
    bonds2 = [(int(perm_inv[int(i)]), int(perm_inv[int(j)])) for (i, j) in bonds]
    nnode = int(np.sum(is_node))
    return elems2, xyz2, bonds2, nnode


def load_molecule_any_xyz(xyz_path):
    """Load XYZ, find bonds via AtomicSystem, reorder nodes first."""
    import os, sys
    _pyball_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
    if _pyball_dir not in sys.path:
        sys.path.insert(0, _pyball_dir)
    from pyBall.AtomicSystem import AtomicSystem
    mol = AtomicSystem(fname=xyz_path)
    if mol.bonds is None or len(mol.bonds) == 0:
        mol.findBonds()
    bonds = [(int(b[0]), int(b[1])) for b in mol.bonds]
    elems = list(mol.enames)
    xyz = np.asarray(mol.apos[:, :3], dtype=np.float32)
    elems, xyz, bonds, nnode = reorder_nodes_first(elems, xyz, bonds)
    return elems, xyz, bonds, nnode


def build_packed_system_from_xyz(xyz_path, group_size=64):
    """Load an XYZ file and build a packed RRsp3 sysdata dict.

    Returns the same dict shape as build_packed_system.
    """
    from RRsp3 import build_neighs_bk_from_bonds, make_bk_slots_clustered, make_exclusions_1st_2nd
    elems, xyz, bonds, nnode = load_molecule_any_xyz(xyz_path)
    mol = {'elems': elems, 'pos': xyz, 'bonds': bonds, 'nnode': nnode}
    packed = pack_molecules_contiguous([mol], group_size=group_size, nodes_first=True, pad_to_group=True)
    natoms = int(packed['natoms_total'])
    pos = packed['pos']
    elems_all = packed['elems']
    is_pad = packed['is_padding']
    real = ~is_pad

    elems_real = [e for e in elems_all if e != 'X']
    m_real = masses_from_elems(elems_real)
    m = np.zeros((natoms,), dtype=np.float32)
    m[real] = m_real
    invm = np.zeros_like(m)
    invm[real] = 1.0 / m[real]

    neighs, _ = build_neighs_bk_from_bonds(natoms, packed['bonds'], max_deg=4)
    excl1, excl2 = make_exclusions_1st_2nd(neighs)

    nnode_per_group = int(packed['nnode_group'][0])
    bkSlots = make_bk_slots_clustered(neighs, group_size=group_size, nnode_per_group=nnode_per_group, natoms=natoms)

    port_local_atoms, K_atoms = make_ports_from_neighs(pos, neighs, K=200.0)

    rad = np.zeros((natoms,), dtype=np.float32)
    rad[real] = 1.0

    fixmask = np.zeros((natoms,), dtype=np.int32)
    fixmask[is_pad] |= (1 | 2 | 4)

    quat = np.zeros((natoms, 4), dtype=np.float32)
    quat[:, 3] = 1.0

    return {
        'natoms': natoms,
        'pos': pos,
        'quat': quat,
        'invm': invm,
        'm': m,
        'elems_all': elems_all,
        'is_pad': is_pad,
        'real': real,
        'neighs': neighs,
        'excl1': excl1,
        'excl2': excl2,
        'nnode_per_group': nnode_per_group,
        'bkSlots': bkSlots,
        'port_local_atoms': port_local_atoms,
        'K_atoms': K_atoms,
        'rad': rad,
        'fixmask': fixmask,
        'group_size': group_size,
    }


def quat_rotate_vec(q, v):
    """Rotate vector(s) v by quaternion(s) q. q=[x,y,z,w]."""
    q = np.asarray(q, dtype=np.float32)
    v = np.asarray(v, dtype=np.float32)
    if q.shape[-1] != 4:
        raise ValueError("quat_rotate_vec: q must have last dim 4")
    if v.shape[-1] != 3:
        raise ValueError("quat_rotate_vec: v must have last dim 3")
    t = 2.0 * np.cross(q[..., :3], v)
    return v + q[..., 3:4] * t + np.cross(q[..., :3], t)


def compute_mean_constraint_error(pos, quat, neighs, port_local):
    """Mean port constraint violation: average ||pos[j] - tip_i|| over valid bonds."""
    pos = np.asarray(pos, dtype=np.float32)
    if pos.shape[1] == 4:
        pos = pos[:, :3]
    quat = np.asarray(quat, dtype=np.float32)
    neighs = np.asarray(neighs, dtype=np.int32)
    port_local = np.asarray(port_local, dtype=np.float32)
    natoms = int(pos.shape[0])
    errs = []
    for i in range(natoms):
        for k in range(4):
            j = int(neighs[i, k])
            if j < 0:
                continue
            tip = pos[i] + quat_rotate_vec(quat[i], port_local[i, k, :3])
            errs.append(np.linalg.norm(pos[j] - tip))
    return float(np.mean(errs)) if errs else 0.0


# ------------------------------------------------------------------
# Momentum Conservation Analysis
# ------------------------------------------------------------------

def gather_dx_port(natoms, nnode_per_group, group_size, dpos_node, dpos_neigh, bkSlots):
    """Gather per-atom port displacement from node + neigh buffers.

    This mirrors the kernel gather logic in apply_corrections_rigid_ports.
    dpos_neigh is (nnode_tot, 4, 4) where slot = inode*4 + k maps to (inode, k).
    """
    dx_port = np.zeros((natoms, 3), dtype=np.float32)
    ng = natoms // group_size
    for ig in range(ng):
        abase = ig * group_size
        inode_base = ig * nnode_per_group
        for il in range(nnode_per_group):
            ia = abase + il
            inode = inode_base + il
            dx_port[ia, :] += dpos_node[inode, :3]
    for ia in range(natoms):
        for k in range(4):
            slot = int(bkSlots[ia, k])
            if slot >= 0:
                inode = slot // 4
                k_neigh = slot % 4
                dx_port[ia, :] += dpos_neigh[inode, k_neigh, :3]
    return dx_port


def compute_momentum_deltas(pos, m, dx_coll, dx_port, dpos_node, drot_node, nnode_per_group, group_size, bkSlots, relax=1.0, massless_rot=False):
    """Compute linear and angular momentum change from correction deltas.

    Args:
        massless_rot: If True, rotation is massless geometric alignment (eigen/shapematch/substep).
                      In this case, do NOT include rotational angular momentum from dtheta in total dL.
                      Only check that linear displacement conserves momentum.
                      If False, rotation carries physical inertia (current/orig kernels).

    Returns (dP, dL) as float3 vectors.
    """
    natoms = int(pos.shape[0])
    dx = (dx_coll + dx_port) * float(relax)
    mdx = dx * m[:, None]
    dP = np.sum(mdx, axis=0)

    real = m > 0
    com = np.sum(pos[real] * m[real, None], axis=0) / (np.sum(m[real]) + 1e-12)
    r = pos[real] - com[None, :]
    dL_trans = np.sum(np.cross(r, mdx[real]), axis=0)

    if massless_rot:
        dL = dL_trans
    else:
        I = 0.4 * m
        dtheta = np.zeros((natoms, 3), dtype=np.float32)
        ng = natoms // group_size
        for ig in range(ng):
            abase = ig * group_size
            inode_base = ig * nnode_per_group
            for il in range(nnode_per_group):
                ia = abase + il
                inode = inode_base + il
                dtheta[ia, :] = drot_node[inode, :3] * float(relax)
        dL_rot = np.sum((I[:, None] * dtheta), axis=0)
        dL = dL_trans + dL_rot
    return dP.astype(np.float32), dL.astype(np.float32)


# ------------------------------------------------------------------
# Topology / Interaction Truth
# ------------------------------------------------------------------

def build_interaction_truth(molecule_names, packed, shift, collision_radius=1.0):
    """Auto-generate expected interaction table from molecular topology.

    Returns dict mapping (min_id, max_id) -> action string.
    Actions: 'EXCLUDE_BOND', 'EXCLUDE_ANGLE', 'COLLIDE', 'IGNORE_FAR'
    """
    mols = [make_test_system(n) for n in molecule_names]
    natoms_per_mol = [len(m['elems']) for m in mols]
    truth = {}

    def add(a, b, act):
        truth[(min(a, b), max(a, b))] = act

    # Intra-molecular exclusions
    offset = 0
    for mol in mols:
        n = len(mol['elems'])
        # bond exclusions (1-2)
        for (i, j) in mol['bonds']:
            add(offset + i, offset + j, 'EXCLUDE_BOND')
        # angle exclusions (1-3) from bond graph
        bond_dict = {}
        for (i, j) in mol['bonds']:
            bond_dict.setdefault(i, []).append(j)
            bond_dict.setdefault(j, []).append(i)
        for i in range(n):
            for j in bond_dict.get(i, []):
                for k in bond_dict.get(j, []):
                    if k != i:
                        add(offset + i, offset + k, 'EXCLUDE_ANGLE')
        offset += n

    # Inter-molecular collisions: compute distances using shifts
    if shift is None:
        shift = [np.zeros(3, dtype=np.float32) for _ in mols]
    elif isinstance(shift, (list, tuple)) and len(shift) == 3 and isinstance(shift[0], (int, float)):
        base = np.array(shift, dtype=np.float32)
        shift = [base * float(i) for i in range(len(mols))]
    else:
        shift = [np.array(s, dtype=np.float32) for s in shift]

    for i in range(len(mols)):
        for j in range(i + 1, len(mols)):
            for ai in range(len(mols[i]['elems'])):
                for aj in range(len(mols[j]['elems'])):
                    pi = mols[i]['pos'][ai] + shift[i]
                    pj = mols[j]['pos'][aj] + shift[j]
                    dist = float(np.linalg.norm(pi - pj))
                    gid_i = sum(natoms_per_mol[:i]) + ai
                    gid_j = sum(natoms_per_mol[:j]) + aj
                    if dist < 2.0 * collision_radius:
                        add(gid_i, gid_j, 'COLLIDE')
                    else:
                        add(gid_i, gid_j, 'IGNORE_FAR')

    return truth


def check_local_ranges(neighs_local, excl1_local, excl2_local, ghost_counts, *, group_size=64):
    """Validate that all local indices fit within group_size + ghost_count."""
    neighs_local = np.asarray(neighs_local, dtype=np.int32)
    excl1_local = np.asarray(excl1_local, dtype=np.int32)
    excl2_local = np.asarray(excl2_local, dtype=np.int32)
    ghost_counts = np.asarray(ghost_counts, dtype=np.int32)
    natoms = int(neighs_local.shape[0])
    if neighs_local.shape != (natoms, 4):
        raise ValueError(f"check_local_ranges: neighs_local.shape={neighs_local.shape} expected ({natoms},4)")

    ng = natoms // int(group_size)
    for ig in range(ng):
        g = int(ghost_counts[ig])
        hi = int(group_size) + g
        abase = ig * int(group_size)
        arrs = (neighs_local[abase:abase + group_size],
                excl1_local[abase:abase + group_size],
                excl2_local[abase:abase + group_size])
        for A in arrs:
            bad = (A >= hi) & (A != -1)
            if np.any(bad):
                w = np.where(bad)
                raise RuntimeError(f"local index out of range in group {ig}: hi={hi} example atom={w[0][0]} k={w[1][0]} val={A[w[0][0], w[1][0]]}")


# ------------------------------------------------------------------
# Debug Log Parsing
# ------------------------------------------------------------------

REGEX_COLL = re.compile(r"^COLL:\s*MeG=(\d+)\s+OtherL=(\d+)\s+OtherG=(\d+)\s+Dist=([0-9eE+\.-]+)\s+Action=(\w+)")
REGEX_TOPO = re.compile(r"^TOPOLOGY:\s*GID=(\d+)\s+LID=(\d+)\s+n_ghosts=(\d+)")
REGEX_MAP = re.compile(r"^MAP:\s*GID=(\d+)\s+L_IDX=(\d+)\s+G_IDX=(\d+)\s+TYPE=(\w+)")


def parse_kernel_debug_logs(stdout_text):
    """Parse COLL, TOPOLOGY, MAP lines from kernel stdout.

    Returns dict with keys: 'coll', 'topology', 'mapping'.
    Each is a list of parsed dicts.
    """
    out = {'coll': [], 'topology': [], 'mapping': []}
    for line in stdout_text.splitlines():
        line = line.strip()
        m = REGEX_COLL.match(line)
        if m:
            out['coll'].append({
                'me_g': int(m.group(1)),
                'other_l': int(m.group(2)),
                'other_g': int(m.group(3)),
                'dist': float(m.group(4)),
                'action': m.group(5),
            })
            continue
        m = REGEX_TOPO.match(line)
        if m:
            out['topology'].append({
                'gid': int(m.group(1)),
                'lid': int(m.group(2)),
                'n_ghosts': int(m.group(3)),
            })
            continue
        m = REGEX_MAP.match(line)
        if m:
            out['mapping'].append({
                'gid': int(m.group(1)),
                'l_idx': int(m.group(2)),
                'g_idx': int(m.group(3)),
                'type': m.group(4),
            })
            continue
    return out


def validate_interactions(parsed_logs, truth, group_size=64, dbg_gids=None):
    """Validate parsed collision logs against interaction truth table.

    Args:
        parsed_logs: output from parse_kernel_debug_logs()
        truth: dict from build_interaction_truth()
        group_size: 64 for RRsp3
        dbg_gids: set of global IDs that were inside the debug print window.
                  If None, requires logs for all truth entries.

    Returns:
        (pass_bool, report_list)
    """
    report = []
    passed = True
    seen = {}

    # Aggregate: keep strongest action per pair
    rank = {'COLLIDE': 3, 'SKIP_EXCL': 2, 'TOO_FAR': 1}
    for entry in parsed_logs.get('coll', []):
        a = int(entry['me_g'])
        b = int(entry['other_g'])
        k = (min(a, b), max(a, b))
        act = entry['action']
        prev = seen.get(k)
        if prev is None or rank.get(act, 0) > rank.get(prev, 0):
            seen[k] = act

    # Validate truth entries
    for (a, b), exp in truth.items():
        if exp.startswith('EXCLUDE'):
            # Only require if at least one atom was in debug window
            if dbg_gids is not None and not (a in dbg_gids or b in dbg_gids):
                continue
            got = seen.get((a, b))
            if got is None:
                report.append(f"FAIL: missing log for excluded pair {a}-{b} (expected {exp})")
                passed = False
            elif got != 'SKIP_EXCL':
                report.append(f"FAIL: excluded pair {a}-{b} expected SKIP_EXCL got {got}")
                passed = False
            else:
                report.append(f"PASS: excluded pair {a}-{b} correctly SKIP_EXCL")
        elif exp == 'COLLIDE':
            if dbg_gids is not None and not (a in dbg_gids or b in dbg_gids):
                continue
            got = seen.get((a, b))
            if got is None:
                report.append(f"INFO: no log for expected COLLIDE pair {a}-{b}")
            elif got != 'COLLIDE':
                report.append(f"FAIL: pair {a}-{b} expected COLLIDE got {got}")
                passed = False
            else:
                report.append(f"PASS: pair {a}-{b} correctly COLLIDE")

    # Check for unexpected inter-molecular collisions
    inter_coll = []
    for (a, b), act in seen.items():
        if (a < group_size and b >= group_size) and act == 'COLLIDE':
            inter_coll.append((a, b))
    if inter_coll:
        report.append(f"INFO: {len(inter_coll)} inter-molecular COLLIDE events logged")

    return passed, report


# ------------------------------------------------------------------
# Output / Reporting
# ------------------------------------------------------------------

def write_xyz_frame(f, elems, pos, comment=""):
    pos = np.asarray(pos, dtype=np.float32)
    f.write(f"{pos.shape[0]}\n")
    f.write(f"{comment}\n")
    for e, p in zip(elems, pos):
        f.write(f"{e} {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")


def write_test_report(path, data):
    """Write AI-evaluable JSON report."""
    with open(path, 'w') as f:
        json.dump(data, f, indent=2, default=str)


def print_report(data):
    """Print structured report to stdout (AI-parseable)."""
    print("=" * 50)
    print("RRSP3_TEST_REPORT")
    print("=" * 50)
    for k, v in data.items():
        if isinstance(v, (list, tuple)) and len(v) > 10:
            print(f"{k}: [{len(v)} items]")
        else:
            print(f"{k}: {v}")
    print("=" * 50)


# ------------------------------------------------------------------
# 2D Planar Constraint Helpers
# ------------------------------------------------------------------

def make_planar_fixmask(natoms, plane='xy'):
    """Return fixmask array that constrains atoms to a plane.

    plane='xy' -> fix z (bit 4)
    plane='xz' -> fix y (bit 2)
    plane='yz' -> fix x (bit 1)
    """
    fixmask = np.zeros((natoms,), dtype=np.int32)
    bit = {'xy': 4, 'xz': 2, 'yz': 1}.get(plane, 4)
    fixmask[:] |= bit
    return fixmask


def perturb_state_planar(pos, quat, pos_scale=0.05, rot_scale=0.0, rng=None, plane='xy'):
    """Perturb state but keep atoms in plane."""
    if rng is None:
        rng = np.random.default_rng(0)
    pos = np.asarray(pos, dtype=np.float32).copy()
    quat = np.asarray(quat, dtype=np.float32).copy()
    n = pos.shape[0]
    pos += rng.normal(0.0, float(pos_scale), size=(n, 3)).astype(np.float32)
    if plane == 'xy':
        pos[:, 2] = 0.0
    elif plane == 'xz':
        pos[:, 1] = 0.0
    elif plane == 'yz':
        pos[:, 0] = 0.0
    if rot_scale > 0.0:
        # small random rotation around plane-normal axis only
        axis_idx = {'xy': 2, 'xz': 1, 'yz': 0}.get(plane, 2)
        for i in range(n):
            angle = rng.normal(0.0, float(rot_scale))
            ax = np.zeros(3, dtype=np.float32)
            ax[axis_idx] = 1.0
            # quaternion from axis-angle
            half = angle * 0.5
            qrot = np.array([ax[0] * np.sin(half), ax[1] * np.sin(half), ax[2] * np.sin(half), np.cos(half)], dtype=np.float32)
            # multiply qrot * quat
            x1, y1, z1, w1 = qrot
            x2, y2, z2, w2 = quat[i]
            quat[i] = np.array([
                w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
                w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
                w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2,
                w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
            ], dtype=np.float32)
    return pos, quat


# ------------------------------------------------------------------
# Collision Reference & Reindexing Validation
# ------------------------------------------------------------------

def brute_force_collision_dpos(pos, radius, excl1, excl2, invm=None):
    """Fully vectorized O(n^2) brute-force collision displacements.

    Matches GPU kernel logic: dl = (rsum - dist) / (w_i + w_j) * 0.5.
    Uses numpy advanced indexing for exclusions — no Python atom loops.
    """
    pos = np.asarray(pos, dtype=np.float32)
    radius = np.asarray(radius, dtype=np.float32)
    excl1 = np.asarray(excl1, dtype=np.int32)
    excl2 = np.asarray(excl2, dtype=np.int32)
    n = pos.shape[0]

    # Pairwise displacement (n, n, 3) and distances
    d = pos[:, None, :] - pos[None, :, :]           # (n, n, 3)
    dists = np.linalg.norm(d, axis=2)               # (n, n)

    # Exclusion mask via advanced indexing (vectorized, no atom loops)
    excl = np.concatenate([excl1, excl2], axis=1)   # (n, 8)
    ii = np.broadcast_to(np.arange(n)[:, None], excl.shape)
    valid = excl >= 0
    excl_mask = np.zeros((n, n), dtype=bool)
    excl_mask[ii[valid], excl[valid]] = True

    # Mass weights
    if invm is None:
        invm = np.where(radius > 0, 1.0, 0.0).astype(np.float32)
    else:
        invm = np.asarray(invm, dtype=np.float32)

    # Collision mask
    rsum = radius[:, None] + radius[None, :]
    coll_mask = (dists < rsum) & (dists > 1e-6) & (invm[:, None] > 0) & (invm[None, :] > 0)
    np.fill_diagonal(coll_mask, False)
    coll_mask &= ~excl_mask

    # Impulses (same formula as GPU: halved symmetric repulsion)
    w_tot = invm[:, None] + invm[None, :] + 1e-12
    dl = np.where(coll_mask, (rsum - dists) / w_tot * 0.5, 0.0)

    # Normals + displacement for each atom
    n_vec = d / (dists[:, :, None] + 1e-12)
    dpos = np.sum(n_vec * dl[:, :, None] * invm[:, None, None], axis=1)
    return dpos.astype(np.float32)


def compute_wg_aabbs(pos, radius, group_size):
    """Compute per-workgroup AABBs from positions and radii."""
    pos = np.asarray(pos, dtype=np.float32)
    radius = np.asarray(radius, dtype=np.float32)
    n = pos.shape[0]
    ng = n // group_size
    pos_r = pos.reshape(ng, group_size, 3)
    rad_r = radius.reshape(ng, group_size, 1)
    real = rad_r[..., 0] > 0
    pos_masked = np.where(real[:, :, None], pos_r, np.nan)
    aabb_min = np.nanmin(pos_masked - rad_r, axis=1)
    aabb_max = np.nanmax(pos_masked + rad_r, axis=1)
    return aabb_min, aabb_max


def validate_ghost_list(pos, radius, ghost_indices_flat, ghost_counts, group_size, margin_sq):
    """Check that every ghost atom is within margin_sq of its workgroup's AABB.

    ghost_indices_flat: 1D array of size ng * max_ghosts.
    Returns (ok, report_list).
    """
    pos = np.asarray(pos, dtype=np.float32)
    radius = np.asarray(radius, dtype=np.float32)
    ghost_counts = np.asarray(ghost_counts, dtype=np.int32)
    ghost_indices_flat = np.asarray(ghost_indices_flat, dtype=np.int32)
    ng = len(ghost_counts)
    n = pos.shape[0]
    max_ghosts = ghost_indices_flat.shape[0] // ng
    ghost_2d = ghost_indices_flat.reshape(ng, max_ghosts)

    aabb_min, aabb_max = compute_wg_aabbs(pos, radius, group_size)
    margin = np.sqrt(float(margin_sq)) + 1e-6

    ok = True
    report = []
    for grp in range(ng):
        g_count = int(ghost_counts[grp])
        if g_count <= 0:
            continue
        ghosts = ghost_2d[grp, :g_count]
        # No duplicates
        if len(set(ghosts.tolist())) != g_count:
            report.append(f"FAIL group {grp}: duplicate ghost indices")
            ok = False
        # Each ghost within margin of workgroup AABB
        gmin = aabb_min[grp]
        gmax = aabb_max[grp]
        gpos = pos[ghosts]
        dx = np.maximum(0.0, np.maximum(gmin[0] - gpos[:, 0], gpos[:, 0] - gmax[0]))
        dy = np.maximum(0.0, np.maximum(gmin[1] - gpos[:, 1], gpos[:, 1] - gmax[1]))
        dz = np.maximum(0.0, np.maximum(gmin[2] - gpos[:, 2], gpos[:, 2] - gmax[2]))
        dists_sq = dx*dx + dy*dy + dz*dz
        bad = np.where(dists_sq > float(margin_sq) + 1e-6)[0]
        if len(bad) > 0:
            report.append(f"FAIL group {grp}: {len(bad)} ghosts beyond margin_sq (max dist={np.sqrt(dists_sq[bad].max()):.3f})")
            ok = False
        else:
            report.append(f"PASS group {grp}: {g_count} ghosts within margin")
    return ok, report


def validate_neighs_local(neighs_global, neighs_local, ghost_indices_flat, ghost_counts, group_size):
    """Round-trip check: global -> local -> global for neighbor and exclusion lists.

    Returns (ok, report_list).
    """
    neighs_global = np.asarray(neighs_global, dtype=np.int32)
    neighs_local = np.asarray(neighs_local, dtype=np.int32)
    ghost_counts = np.asarray(ghost_counts, dtype=np.int32)
    ghost_indices_flat = np.asarray(ghost_indices_flat, dtype=np.int32)
    n = neighs_global.shape[0]
    ng = n // group_size
    max_ghosts = ghost_indices_flat.shape[0] // ng
    ghost_2d = ghost_indices_flat.reshape(ng, max_ghosts)

    ok = True
    report = []
    for grp in range(ng):
        wg_start = grp * group_size
        wg_end = min(wg_start + group_size, n)
        g_count = int(ghost_counts[grp])
        errors = 0
        for i in range(wg_start, wg_end):
            for k in range(4):
                jglob = int(neighs_global[i, k])
                jloc = int(neighs_local[i, k])
                if jglob < 0:
                    if jloc >= 0:
                        report.append(f"FAIL atom {i} slot {k}: global=-1 but local={jloc}")
                        ok = False; errors += 1
                    continue
                if jloc < 0:
                    report.append(f"FAIL atom {i} slot {k}: global={jglob} mapped to local=-1")
                    ok = False; errors += 1
                    continue
                if jloc < group_size:
                    j_back = grp * group_size + jloc
                else:
                    gi = jloc - group_size
                    if gi >= g_count:
                        report.append(f"FAIL atom {i} slot {k}: ghost idx {gi} out of range (count={g_count})")
                        ok = False; errors += 1
                        continue
                    j_back = int(ghost_2d[grp, gi])
                if j_back != jglob:
                    report.append(f"FAIL atom {i} slot {k}: round-trip {jglob} -> {jloc} -> {j_back}")
                    ok = False; errors += 1
        if errors == 0:
            report.append(f"PASS group {grp}: neighbor/exclusion mapping ok")
    return ok, report


def plot_workgroups(pos, group_id, is_padding, outpath='workgroups.png'):
    """3D scatter plot of atoms colored by workgroup ID."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    pos = np.asarray(pos, dtype=np.float32)
    group_id = np.asarray(group_id, dtype=np.int32)
    is_padding = np.asarray(is_padding, dtype=bool)
    real = ~is_padding
    pos_r = pos[real]
    gid_r = group_id[real]

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    n_groups = int(gid_r.max()) + 1
    cmap = plt.cm.get_cmap('tab20', max(n_groups, 20))
    sc = ax.scatter(pos_r[:, 0], pos_r[:, 1], pos_r[:, 2], c=gid_r, cmap=cmap, s=50, alpha=0.8)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('Atoms colored by workgroup')
    plt.colorbar(sc, ax=ax, label='Workgroup ID')
    plt.tight_layout()
    plt.savefig(outpath, dpi=150)
    plt.close()
    print(f"Saved workgroup plot to {outpath}")
