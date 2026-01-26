import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.collections import LineCollection
import argparse
import sys
import os
import itertools
import pyopencl as cl
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from XPDB_new import XPDB_new

from pyBall.OCL.MMFFL import MMFFL
from pyBall import atomicUtils as au

# unlimited line length
np.set_printoptions(threshold=np.inf)
np.set_printoptions(linewidth=np.inf)


def build_bond_arrays_with_angles(
    apos,
    bonds_adj,
    n_max_bonded=16,
    default_L=None,
    default_K=200.0,
    alpha0_deg=120.0,
    k_angle=None,
    enable_angles=True
):
    """
    Build fixed-size bond index/length/stiffness arrays from first-order bonds
    and second-order (angle-derived) neighbors.

    Parameters
    ----------
    apos : np.ndarray (num_atoms, 3)
        Atom positions.
    bonds_adj : list of lists
        Each element is [(j, L0, K), ...] for atom i.
        If L0 or K missing, use default_L / default_K.
    n_max_bonded : int
        Maximum bonded neighbors per atom (e.g., 16).
    default_L : float or None
        Default bond rest length if not provided.
    default_K : float
        Default bond stiffness.
    alpha0_deg : float
        Default relaxed angle (degrees) for second-neighbor rest length.
    k_angle : float or None
        Stiffness for angle-derived bonds; if None, uses default_K.

    Returns
    -------
    bond_indices : np.ndarray (num_atoms, n_max_bonded)
        Neighbor indices (-1 for unused).
    bond_lengths : np.ndarray (num_atoms, n_max_bonded)
        Rest lengths.
    bond_stiffness : np.ndarray (num_atoms, n_max_bonded)
        Stiffness values.
    """
    num_atoms = len(apos)
    if default_L is None:
        default_L = 0.4
    if k_angle is None:
        k_angle = default_K

    alpha0_rad = math.radians(alpha0_deg)

    first_neighbors = [[] for _ in range(num_atoms)]
    for i, blist in enumerate(bonds_adj):
        for b in blist:
            j = b[0]
            L = b[1] if len(b) > 1 and b[1] is not None else default_L
            if L is None:
                L = float(np.linalg.norm(apos[i] - apos[j]))
            K = b[2] if len(b) > 2 and b[2] is not None else default_K
            first_neighbors[i].append([j, float(L), float(K)])

    second_neighbors = [[] for _ in range(num_atoms)]
    if enable_angles:
        for i, blist in enumerate(first_neighbors):
            n = len(blist)
            for a in range(n):
                for b in range(a + 1, n):
                    j, Lij, _ = blist[a]
                    k, Lik, _ = blist[b]
                    if Lij is None:
                        Lij = default_L if default_L is not None else float(np.linalg.norm(apos[i] - apos[j]))
                    if Lik is None:
                        Lik = default_L if default_L is not None else float(np.linalg.norm(apos[i] - apos[k]))
                    Ljk = Lij * Lij + Lik * Lik - 2 * Lij * Lik * math.cos(alpha0_rad)
                    if Ljk < 0:
                        Ljk = 0.0
                    Ljk = math.sqrt(Ljk)
                    second_neighbors[j].append([k, Ljk, k_angle])
                    second_neighbors[k].append([j, Ljk, k_angle])

    bond_indices = np.full((num_atoms, n_max_bonded), -1, dtype=np.int32)
    bond_lengths = np.zeros((num_atoms, n_max_bonded), dtype=np.float32)
    bond_stiffness = np.zeros((num_atoms, n_max_bonded), dtype=np.float32)
    bond_type_mask = np.zeros((num_atoms, n_max_bonded), dtype=np.int32)  # 1: primary, 2: secondary

    for i in range(num_atoms):
        merged = first_neighbors[i] + second_neighbors[i]
        seen = set()
        unique = []
        for entry in merged:
            idx = entry[0]
            if idx not in seen:
                seen.add(idx)
                unique.append(entry)
        if len(unique) > n_max_bonded:
            print(f"[WARN] Atom {i} has {len(unique)} bonded neighbors > n_max_bonded={n_max_bonded}; truncating.")
        for k, entry in enumerate(unique[:n_max_bonded]):
            bond_indices[i, k] = entry[0]
            bond_lengths[i, k] = entry[1]
            bond_stiffness[i, k] = entry[2]
            if entry in first_neighbors[i]:
                bond_type_mask[i, k] = 1
            else:
                bond_type_mask[i, k] = 2

    return bond_indices, bond_lengths, bond_stiffness, bond_type_mask


def load_molecule_bonds(fname, n_max_bonded=16, default_L=None, default_K=200.0, alpha0_deg=120.0, k_angle=None, enable_angles=True, print_topology=False):
    """
    Load a molecule file and build bond arrays with angle-derived second neighbors.

    Returns
    -------
    apos : np.ndarray (num_atoms, 3)
    bond_indices, bond_lengths, bond_stiffness : arrays ready for upload to GPU
    bonds_adj : list of lists
    bond_type_mask : np.ndarray (num_atoms, n_max_bonded)
    """
    raise RuntimeError("load_molecule_bonds() is deprecated; use load_molecule_topology_mmffl()")


def _bonds_to_adj(n, bonds, *, default_L=1.0, default_K=200.0, apos=None):
    bonds_adj = [[] for _ in range(n)]
    for (i, j) in bonds:
        i = int(i); j = int(j)
        if i == j:
            continue
        if apos is not None and default_L is None:
            L = float(np.linalg.norm(apos[i] - apos[j]))
        else:
            L = float(default_L) if default_L is not None else float(np.linalg.norm(apos[i] - apos[j]))
        bonds_adj[i].append([j, L, float(default_K)])
        bonds_adj[j].append([i, L, float(default_K)])
    return bonds_adj


def _pack_fixed_bonds(n, bonds_adj, *, n_max_bonded=16):
    bond_indices = np.full((n, n_max_bonded), -1, dtype=np.int32)
    bond_lengths = np.zeros((n, n_max_bonded), dtype=np.float32)
    bond_stiffness = np.zeros((n, n_max_bonded), dtype=np.float32)
    bond_type_mask = np.zeros((n, n_max_bonded), dtype=np.int32)
    for i in range(n):
        blist = bonds_adj[i]
        if len(blist) > n_max_bonded:
            print(f"[WARN] Atom {i} has {len(blist)} neighbors > n_max_bonded={n_max_bonded}; truncating")
        for k, (j, L, K) in enumerate(blist[:n_max_bonded]):
            bond_indices[i, k] = int(j)
            bond_lengths[i, k] = float(L)
            bond_stiffness[i, k] = float(K)
            bond_type_mask[i, k] = 1
    return bond_indices, bond_lengths, bond_stiffness, bond_type_mask


def load_molecule_topology_mmffl(
    fname,
    *,
    n_max_bonded=16,
    default_K=200.0,
    type_source="table",
    add_angle=True,
    add_pi=False,
    two_pi_dummies=False,
    add_pi_align=False,
    add_epair=False,
    add_epair_pairs=False,
    L_pi=1.0,
    L_epair=0.5,
    K_ang=0.0,
    Kpi_host=0.0,
    Kpi_orth=0.0,
    Kpi_align=0.0,
    Kep_host=0.0,
    Kep_orth=0.0,
    Kep_pair=0.0,
    print_topology=False,
    mol2_out=None,
):
    from pyBall import AtomicSystem
    mol = AtomicSystem.AtomicSystem(fname=fname)
    if mol.bonds is None or len(mol.bonds) == 0:
        mol.findBonds()
    mol.neighs()

    ff = MMFFL(
        L_pi=L_pi,
        two_pi_dummies=two_pi_dummies,
        Kang=K_ang,
        Kpi_host=Kpi_host,
        Kpi_orth=Kpi_orth,
        Kpi_align=Kpi_align,
        L_epair=L_epair,
        Kep_host=Kep_host,
        Kep_orth=Kep_orth,
        Kep_pair=Kep_pair,
        verbosity=1,
    )
    topo = ff.build_topology(
        mol,
        type_source=type_source,
        add_angle=bool(add_angle),
        add_pi=bool(add_pi),
        two_pi_dummies=bool(two_pi_dummies),
        add_pi_align=bool(add_pi_align),
        add_epair=bool(add_epair),
        add_epair_pairs=bool(add_epair_pairs),
        bUFF=False,
        report=True,
    )
    apos_all = topo["apos"].astype(np.float32)
    n_all = int(topo["n_all"])
    n_real = int(topo["n_real"])

    # Build adjacency lists for plotting and for XPDB buffers
    bonds_primary = topo["bonds_primary"]
    bonds_derived = []
    if add_angle:
        bonds_derived += topo["bonds_angle"]
    if add_pi:
        bonds_derived += topo["bonds_pi"]
    if add_pi and add_pi_align:
        bonds_derived += topo["bonds_pi_align"]
    if add_epair:
        bonds_derived += topo["bonds_epair"]
    if add_epair and add_epair_pairs:
        bonds_derived += topo["bonds_epair_pair"]

    bonds_all = sorted(set(bonds_primary + bonds_derived))
    bonds_adj = _bonds_to_adj(n_all, bonds_all, default_L=None, default_K=default_K, apos=apos_all)
    bond_indices, bond_lengths, bond_stiffness, bond_type_mask = _pack_fixed_bonds(n_all, bonds_adj, n_max_bonded=n_max_bonded)

    if print_topology:
        print("[TOPO] n_real=", n_real, " n_all=", n_all, " n_primary=", len(bonds_primary), " n_derived=", len(bonds_derived))
        print("[TOPO] type_names:")
        for ia, tn in enumerate(topo["type_names"]):
            print(f"  {ia:4d} {tn} {'DUMMY' if topo['is_dummy'][ia] else ''}")
        print("[TOPO] bonds_primary:", bonds_primary)
        print("[TOPO] bonds_angle:", topo["bonds_angle"])
        print("[TOPO] bonds_pi:", topo["bonds_pi"])
        if add_pi and add_pi_align:
            print("[TOPO] bonds_pi_align:", topo["bonds_pi_align"])
        print("[TOPO] bonds_epair:", topo["bonds_epair"])
        if add_epair and add_epair_pairs:
            print("[TOPO] bonds_epair_pair:", topo["bonds_epair_pair"])

    if mol2_out is not None:
        au.save_mol2(mol2_out, topo["type_names"], apos_all, bonds_all, qs=None, comment=f"n_real={n_real} n_all={n_all}")
        print("[TOPO] wrote", mol2_out)

    return topo, apos_all, bond_indices, bond_lengths, bond_stiffness, bonds_adj, bond_type_mask


def plot_topology_3d(apos, topo, *, title="Topology"):
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_title(title)
    p = np.asarray(apos, dtype=float)
    is_dummy = np.asarray(topo["is_dummy"], dtype=bool)
    type_names = np.asarray(topo["type_names"], dtype=object)
    is_pi_dummy = is_dummy & np.array([str(tn).startswith("Pi") for tn in type_names])
    is_ep_dummy = is_dummy & np.array([str(tn).startswith("E") for tn in type_names])
    is_other_dummy = is_dummy & ~(is_pi_dummy | is_ep_dummy)
    # atoms
    ax.scatter(p[~is_dummy, 0], p[~is_dummy, 1], p[~is_dummy, 2], s=80, c='k', depthshade=True, label='atoms')
    if np.any(is_pi_dummy):
        ax.scatter(p[is_pi_dummy, 0], p[is_pi_dummy, 1], p[is_pi_dummy, 2], s=50, c='magenta', depthshade=True, label='pi dummies')
    if np.any(is_ep_dummy):
        ax.scatter(p[is_ep_dummy, 0], p[is_ep_dummy, 1], p[is_ep_dummy, 2], s=50, c='cyan', depthshade=True, label='epair dummies')
    if np.any(is_other_dummy):
        ax.scatter(p[is_other_dummy, 0], p[is_other_dummy, 1], p[is_other_dummy, 2], s=40, c='r', depthshade=True, label='other dummies')

    def draw_bonds(bonds, color, lw, label):
        for (i, j) in bonds:
            i = int(i); j = int(j)
            ax.plot([p[i, 0], p[j, 0]], [p[i, 1], p[j, 1]], [p[i, 2], p[j, 2]], color=color, lw=lw)

    draw_bonds(topo["bonds_primary"], color='gray', lw=2.0, label='primary')
    if len(topo["bonds_angle"]) > 0:
        draw_bonds(topo["bonds_angle"], color='cyan', lw=1.0, label='angle')
    if len(topo["bonds_pi"]) > 0:
        draw_bonds(topo["bonds_pi"], color='magenta', lw=1.0, label='pi')
    if len(topo.get("bonds_pi_align", [])) > 0:
        draw_bonds(topo["bonds_pi_align"], color='purple', lw=1.0, label='pi-align')
    if len(topo["bonds_epair"]) > 0:
        draw_bonds(topo["bonds_epair"], color='orange', lw=1.0, label='epair')
    if len(topo.get("bonds_epair_pair", [])) > 0:
        draw_bonds(topo["bonds_epair_pair"], color='brown', lw=1.0, label='epair-pair')

    # labels
    for ia, tn in enumerate(topo["type_names"]):
        ax.text(p[ia, 0], p[ia, 1], p[ia, 2], f"{ia}:{tn}", fontsize=8)

    # aspect equal-ish
    mins = p.min(axis=0); maxs = p.max(axis=0)
    c = 0.5 * (mins + maxs)
    r = float(np.max(maxs - mins)) * 0.6 + 1e-6
    ax.set_xlim(c[0] - r, c[0] + r)
    ax.set_ylim(c[1] - r, c[1] + r)
    ax.set_zlim(c[2] - r, c[2] + r)
    ax.legend(loc='upper right')
    plt.show()


def make_dense_block(nx, ny, spacing, jitter, atom_rad):
    """Create a dense block of atoms"""
    pos = []
    for i in range(nx):
        for j in range(ny):
            x = (i - nx/2) * spacing + np.random.uniform(-jitter, jitter)
            y = (j - ny/2) * spacing + np.random.uniform(-jitter, jitter)
            pos.append([x, y, 0.0])
    pos = np.array(pos, dtype=np.float32)
    vel = np.zeros_like(pos)
    radius = np.full(len(pos), atom_rad, dtype=np.float32)
    mass = np.ones(len(pos), dtype=np.float32)
    return pos, vel, radius, mass


def make_random_bonds(pos, prob=0.25, rest_len=0.4, k_bond=200.0):
    """Create random bonds between nearby atoms"""
    num_atoms = len(pos)
    bonds_adj = [[] for _ in range(num_atoms)]
    
    for i in range(num_atoms):
        for j in range(i+1, num_atoms):
            dist = np.linalg.norm(pos[i] - pos[j])
            if dist < rest_len * 1.5 and np.random.random() < prob:
                bonds_adj[i].append([j, rest_len, k_bond])
                bonds_adj[j].append([i, rest_len, k_bond])
    
    return bonds_adj


def update_plot(scatter, scatter_sizes, bonds_primary, bonds_secondary, line_segments_primary, line_segments_secondary, force_lines, p, f_bond=None, f_coll=None):
    """Update plot for animation"""
    has_primary = bool(bonds_primary) and any(bonds_primary)
    has_secondary = bool(bonds_secondary) and any(bonds_secondary)
    scatter.set_offsets(p[:, :2])
    scatter.set_sizes(scatter_sizes)
    
    if has_primary:
        segs_primary = []
        for i, blist in enumerate(bonds_primary):
            for j, _, _ in blist:
                if j > i:
                    segs_primary.append([p[i, :2], p[j, :2]])
        line_segments_primary.set_segments(segs_primary)
    else:
        line_segments_primary.set_segments([])
    if has_secondary:
        segs_secondary = []
        for i, blist in enumerate(bonds_secondary):
            for j, _, _ in blist:
                if j > i:
                    segs_secondary.append([p[i, :2], p[j, :2]])
        line_segments_secondary.set_segments(segs_secondary)
    else:
        line_segments_secondary.set_segments([])

    # Plot forces as lines
    if f_bond is not None and f_coll is not None:
        f_segs = []
        f_colors = []
        scale = 0.1
        for i in range(len(p)):
            if np.linalg.norm(f_bond[i, :2]) > 1e-6:
                f_segs.append([p[i, :2], p[i, :2] + f_bond[i, :2] * scale])
                f_colors.append((0, 1, 0, 0.8))
            if np.linalg.norm(f_coll[i, :2]) > 1e-6:
                f_segs.append([p[i, :2], p[i, :2] + f_coll[i, :2] * scale])
                f_colors.append((1, 0, 0, 0.8))
        force_lines.set_segments(f_segs)
        force_lines.set_colors(f_colors)

    return (scatter, line_segments_primary, line_segments_secondary, force_lines)

'''

python -u test_TiledJacobi_molecules.py  --molecule ../../../cpp/common_resources/xyz/OCH2.xyz --use_mmffl 1 --enable_angles 1 --add_epair 1 --add_pi 1 --topology_only 1 --mol2_out /tmp/OCH2_topo.mmffl.mol2 --print_buffers 1

python -u test_TiledJacobi_molecules.py --molecule ../../cpp/common_resources/xyz/OCH2.xyz --use_mmffl 1 --type_source table --enable_angles 1 --add_pi 1 --two_pi 1 --add_pi_align 1 --Kpi_align 15 --add_epair 1 --add_epair_pairs 1 --Kep_pair 10 --topology_only 1 --print_buffers 1

python -u test_TiledJacobi_molecules.py  --molecule ../../cpp/common_resources/xyz/OCH2.xyz --topology_only 1 --mol2_out /tmp/OCH2_topo.mmffl.mol2 --print_buffers 1



node web/molgui_webgpu/dump_xpdb_topology.mjs --xyz          ./cpp/common_resources/xyz/guanine.xyz --mol2Out  OUT_js_guanine.mol2
python -u test_TiledJacobi_molecules.py       --molecule ../../cpp/common_resources/xyz/guanine.xyz --mol2_out OUT_py_guanine.mol2 --print_buffers 1 --topology_only 1 

node web/molgui_webgpu/dump_xpdb_topology.mjs --xyz          ./cpp/common_resources/xyz/pyridine.xyz --mol2Out  OUT_js_pyridine.mol2
python -u test_TiledJacobi_molecules.py       --molecule ../../cpp/common_resources/xyz/pyridine.xyz --mol2_out OUT_py_pyridine.mol2 --print_buffers 1 --topology_only 1 

node web/molgui_webgpu/dump_xpdb_topology.mjs --xyz          ./cpp/common_resources/xyz/pyrrol.xyz --mol2Out  OUT_js_pyrrol.mol2
python -u test_TiledJacobi_molecules.py       --molecule ../../cpp/common_resources/xyz/pyrrol.xyz --mol2_out OUT_py_pyrrol.mol2 --print_buffers 1 --topology_only 1 

'''

def run_test():
    parser = argparse.ArgumentParser(description="Tiled Jacobi Solver Test with Molecules (distances in Å, stiffness in eV/Å^2)")
    parser.add_argument("--molecule",        type=str,   default="../../cpp/common_resources/xyz/pentacene.xyz", help="Path to molecule file (.xyz, .mol, .mol2)")
    parser.add_argument("--atom_rad",        type=float, default=0.2, help="Collision sphere radius per atom [Å]")
    parser.add_argument("--omega",           type=float, default=0.7, help="Jacobi relaxation weight ω (dimensionless, 0<ω<1)")
    parser.add_argument("--momentum_beta",   type=float, default=0.0, help="Momentum/velocity blending factor β (dimensionless)")
    parser.add_argument("--dt",              type=float, default=100.0, help="Time step for XPDB updates (solver units, typically fs)")
    parser.add_argument("--k_coll",          type=float, default=100.0, help="Collision spring stiffness [eV/Å^2]")
    parser.add_argument("--bond_k",          type=float, default=200.0, help="Fallback bond stiffness when no MMFFL data [eV/Å^2]")
    parser.add_argument("--bond_len",        type=float, default=1.3, help="Fallback bond relaxed length [Å]")
    parser.add_argument("--solver_steps",    type=int,   default=10,  help="Jacobi iterations per MD step (inner solver loop)")
    parser.add_argument("--md_steps",        type=int,   default=-1, help="Number of PD time steps (outer loop); <0 to run until window close")
    parser.add_argument("--coll_scale",      type=float, default=2.0, help="Collision distance multiplier relative to atom radius (dimensionless)")
    parser.add_argument("--bbox_scale",      type=float, default=2.0, help="Simulation bounding box padding relative to molecule size (dimensionless)")
    parser.add_argument("--group_size",      type=int,   default=64, help="OpenCL work-group size for XPDB kernels")
    parser.add_argument("--max_ghosts",      type=int,   default=128, help="Maximum ghost interactions per atom for collision handling")
    parser.add_argument("--debug_viz",       type=int,   default=1, help="Enable matplotlib preview (1=yes,0=no)")
    parser.add_argument("--force_scale",     type=float, default=1.0, help="Scale factor for drawn force vectors in 2D debug view")
    parser.add_argument("--pick_radius",     type=float, default=0.4, help="Picking radius when interacting with the GUI [Å]")
    parser.add_argument("--alpha0_deg",      type=float, default=120.0, help="Default relaxed angle (degrees) used when no atom-specific data is available")
    parser.add_argument("--k_angle",         type=float, default=None, help="Override stiffness for geometric angle bonds [eV/Å^2]; None -> reuse --bond_k")
    parser.add_argument("--enable_angles",   type=int,   default=1, help="Enable MMFFL angle-derived second neighbors (1=yes,0=no)")
    parser.add_argument("--type_source",     type=str,   default="table", help="Atom type source: 'table' (AtomTypes.dat) or 'uff' (UFF inference)")
    parser.add_argument("--use_mmffl",       type=int,   default=1, help="Use MMFFL topology builder (1) or fall back to legacy geometric bonding (0)")
    parser.add_argument("--add_pi",          type=int,   default=1, help="Enable pi-orbital dummies and their bonds")
    parser.add_argument("--two_pi",          type=int,   default=1, help="If --add_pi, place dummies on both +/− sides of the pi plane")
    parser.add_argument("--add_pi_align",    type=int,   default=0, help="Add alignment bonds between pi dummies of bonded atoms")
    parser.add_argument("--add_epair",       type=int,   default=1, help="Enable electron-pair dummies (lone pairs) and their bonds")
    parser.add_argument("--add_epair_pairs", type=int,   default=0, help="Connect multiple electron-pair dummies on the same host to enforce spacing")
    parser.add_argument("--L_pi",            type=float, default=1.0, help="Distance from host atom to pi dummy center [Å]")
    parser.add_argument("--L_epair",         type=float, default=0.5, help="Distance from host atom to electron-pair dummy [Å]")
    parser.add_argument("--K_ang",           type=float, default=100.0, help="Stiffness for angle-derived bonds (k) between second neighbors; set 0 to disable custom value")
    parser.add_argument("--Kpi_host",        type=float, default=50.0, help="Stiffness for host↔pi-dummy bonds (controls sigma/pi coupling)")
    parser.add_argument("--Kpi_orth",        type=float, default=30.0, help="Stiffness for pi-dummy↔neighbor orthogonal bonds (keeps pi planes aligned)")
    parser.add_argument("--Kpi_align",       type=float, default=15.0, help="Stiffness for pi↔pi alignment bonds between bonded hosts (requires --add_pi_align)")
    parser.add_argument("--Kep_host",        type=float, default=40.0, help="Stiffness for host↔electron-pair dummy bonds")
    parser.add_argument("--Kep_orth",        type=float, default=25.0, help="Stiffness for electron-pair dummy↔neighbor orthogonal bonds")
    parser.add_argument("--Kep_pair",        type=float, default=10.0, help="Stiffness for electron-pair dummy↔dummy bonds on the same host (requires --add_epair_pairs)")
    parser.add_argument("--mol2_out",        type=str,   default=None, help="Optional MOL2 export path (writes Å coordinates + dummy atoms)")
    parser.add_argument("--topology_only",   type=int,   default=0, help="Stop after topology build/plot/export (1=yes,0=run XPDB)")
    parser.add_argument("--print_buffers",   type=int,   default=1, help="Dump XPDB input arrays (positions, bond idx/len/stiff) for inspection")
    parser.add_argument("--print_pos_every", type=int,   default=50, help="Print position bounding box every N frames (0=disable)")
    args = parser.parse_args()

    if args.molecule is None:
        np.random.seed(args.seed)
        pos0, vel0, radius, mass = make_dense_block(args.nx, args.ny, args.spacing, args.jitter, args.atom_rad)
        bonds_plot = make_random_bonds(pos0, prob=args.bond_prob, rest_len=args.bond_len, k_bond=args.bond_k)
        bond_indices, bond_lengths, bond_stiffness, bond_type_mask = build_bond_arrays_with_angles(
            pos0, bonds_plot, n_max_bonded=16, default_L=args.bond_len, default_K=args.bond_k,
            alpha0_deg=args.alpha0_deg, k_angle=args.k_angle, enable_angles=bool(args.enable_angles)
        )
    else:
        if args.use_mmffl:
            topo, pos0, bond_indices, bond_lengths, bond_stiffness, bonds_plot, bond_type_mask = load_molecule_topology_mmffl(
                args.molecule,
                n_max_bonded=16,
                default_K=args.bond_k,
                type_source=args.type_source,
                add_angle=bool(args.enable_angles),
                add_pi=bool(args.add_pi),
                two_pi_dummies=bool(args.two_pi),
                add_pi_align=bool(args.add_pi_align),
                add_epair=bool(args.add_epair),
                add_epair_pairs=bool(args.add_epair_pairs),
                L_pi=args.L_pi,
                L_epair=args.L_epair,
                K_ang=args.K_ang,
                Kpi_host=args.Kpi_host,
                Kpi_orth=args.Kpi_orth,
                Kpi_align=args.Kpi_align,
                Kep_host=args.Kep_host,
                Kep_orth=args.Kep_orth,
                Kep_pair=args.Kep_pair,
                print_topology=bool(args.print_buffers),
                mol2_out=args.mol2_out,
            )
            if args.topology_only:
                plot_topology_3d(pos0, topo, title=os.path.basename(args.molecule))
                return
        else:
            pos0, bond_indices, bond_lengths, bond_stiffness, bonds_plot, bond_type_mask = load_molecule_bonds(
                args.molecule, n_max_bonded=16, default_L=args.bond_len, default_K=args.bond_k,
                alpha0_deg=args.alpha0_deg, k_angle=args.k_angle, enable_angles=bool(args.enable_angles),
                print_topology=bool(args.print_buffers)
            )
        num_atoms = len(pos0)
        vel0 = np.zeros_like(pos0)
        radius = np.full(num_atoms, args.atom_rad, dtype=np.float32)
        mass = np.ones(num_atoms, dtype=np.float32)

    Rmax = float(np.max(radius))
    num_atoms = len(pos0)
    
    sim = XPDB_new(num_atoms, group_size=args.group_size)
    sim.upload_data(pos0, vel0, radius, mass)
    sim.upload_bonds_fixed_from_arrays(bond_indices, bond_lengths, bond_stiffness)

    if args.print_buffers:
        print("[DEBUG] pos0:\n", pos0)
        print("[DEBUG] bond_indices:\n", bond_indices)
        print("[DEBUG] bond_lengths:\n", bond_lengths)
        print("[DEBUG] bond_stiffness:\n", bond_stiffness)
        print("[DEBUG] bond_type_mask:\n", bond_type_mask)
    
    # Setup plot
    fig, ax = plt.subplots(figsize=(8, 8))
    margin = 1.0
    ax.set_xlim(np.min(pos0[:, 0]) - margin, np.max(pos0[:, 0]) + margin)
    ax.set_ylim(np.min(pos0[:, 1]) - margin, np.max(pos0[:, 1]) + margin)
    ax.set_aspect('equal')
    fig.canvas.draw()
    
    scatter_sizes = (radius * 1000).astype(int)
    cmap = plt.cm.get_cmap("tab20")
    cluster_colors = cmap((np.arange(sim.num_groups) % cmap.N) / cmap.N)
    atom_colors = cluster_colors[np.arange(num_atoms) // sim.group_size]
    scatter = ax.scatter(pos0[:, 0], pos0[:, 1], s=scatter_sizes, c=atom_colors, alpha=0.7, edgecolors='k', linewidths=0.3)
    line_segments_primary = LineCollection([], colors='gray', linewidths=1.5)
    line_segments_secondary = LineCollection([], colors='lightgray', linewidths=0.3)
    ax.add_collection(line_segments_secondary)
    ax.add_collection(line_segments_primary)
    
    force_lines = LineCollection([], linewidths=1.0)
    ax.add_collection(force_lines)

    host_pos = np.empty((num_atoms, 4), dtype=np.float32)
    host_pred = np.empty((num_atoms, 4), dtype=np.float32)
    # Seed host_pos with initial positions for picking
    sim.get_positions(out=host_pos)

    pick_state = {
        "idx": None,
        "active": False,
        "mouse": np.array([0.0, 0.0], dtype=np.float32),
    }

    def on_press(event):
        if event.inaxes != ax or event.xdata is None or event.ydata is None:
            return
        mouse_xy = np.array([event.xdata, event.ydata], dtype=np.float32)
        # Compute distance to all atoms (2D)
        d2 = np.sum((host_pos[:num_atoms, :2] - mouse_xy) ** 2, axis=1)
        i_min = int(np.argmin(d2))
        if d2[i_min] <= args.pick_radius * args.pick_radius:
            pick_state["idx"] = i_min
            pick_state["active"] = True
            pick_state["mouse"] = mouse_xy
            print(f"[DEBUG] on_press hit: idx={i_min} d={np.sqrt(d2[i_min]):.4f} radius={args.pick_radius} mouse={mouse_xy}")
        else:
            print(f"[DEBUG] on_press miss: closest idx={i_min} d={np.sqrt(d2[i_min]):.4f} radius={args.pick_radius} mouse={mouse_xy}")
            print(f"[DEBUG] pos closest={host_pos[i_min, :2]}")

    def on_release(event):
        if pick_state["active"]:
            print(f"[DEBUG] on_release: idx={pick_state['idx']}")
        pick_state["active"] = False
        pick_state["idx"] = None

    def on_motion(event):
        if not pick_state["active"]:
            return
        if event.xdata is None or event.ydata is None:
            return
        pick_state["mouse"] = np.array([event.xdata, event.ydata], dtype=np.float32)
        print(f"[DEBUG] on_motion: idx={pick_state['idx']} mouse={pick_state['mouse']}")

    cid_press = fig.canvas.mpl_connect('button_press_event', on_press)
    cid_release = fig.canvas.mpl_connect('button_release_event', on_release)
    cid_motion = fig.canvas.mpl_connect('motion_notify_event', on_motion)
    
    if bonds_plot is not None:
        bonds_plot_primary = []
        bonds_plot_secondary = []
        for i in range(num_atoms):
            primary_list = []
            secondary_list = []
            for k in range( bond_type_mask.shape[1] if bond_type_mask is not None else 0 ):
                j = bond_indices[i, k]
                if j == -1:
                    continue
                entry = [j, bond_lengths[i, k], bond_stiffness[i, k]]
                if bond_type_mask is not None and bond_type_mask[i, k] == 1:
                    primary_list.append(entry)
                elif bond_type_mask is not None and bond_type_mask[i, k] == 2:
                    secondary_list.append(entry)
            bonds_plot_primary.append(primary_list)
            bonds_plot_secondary.append(secondary_list)
    else:
        bonds_plot_primary = bonds_plot
        bonds_plot_secondary = None

    def update(md_step):
        # Build pred targets from current positions, override picked atom to mouse
        host_pred[:, :3] = host_pos[:, :3]
        host_pred[:, 3] = 0.0
        if pick_state["active"] and pick_state["idx"] is not None:
            idx = pick_state["idx"]
            host_pred[idx, 0] = pick_state["mouse"][0]
            host_pred[idx, 1] = pick_state["mouse"][1]
            host_pred[idx, 2] = 0.0
            # Force the current position to the mouse location for a strong drag
            host_pos[idx, 0] = pick_state["mouse"][0]
            host_pos[idx, 1] = pick_state["mouse"][1]
            host_pos[idx, 2] = 0.0
            cl.enqueue_copy(sim.queue, sim.cl_pos, host_pos)
        # Push pred targets to device
        cl.enqueue_copy(sim.queue, sim.cl_pred_pos, host_pred)

        # Build topology
        sim.build_local_topology(Rmax=Rmax, coll_scale=args.coll_scale, bbox_scale=args.bbox_scale)

        # Solve
        sim.solve_cluster_jacobi(
            dt=args.dt,
            iterations=args.solver_steps,
            k_coll=args.k_coll,
            omega=args.omega,
            momentum_beta=args.momentum_beta
        )
        
        # Download positions
        sim.get_positions(out=host_pos)
        if args.print_pos_every > 0 and md_step % args.print_pos_every == 0:
            bbox_min = host_pos[:, :3].min(axis=0)
            bbox_max = host_pos[:, :3].max(axis=0)
            print(f"[DEBUG] md_step {md_step}: bbox min {bbox_min} max {bbox_max}")
        
        # Download forces if debug viz enabled
        if args.debug_viz:
            f_bond, f_coll = sim.get_debug_forces()
            f_bond *= args.force_scale
            f_coll *= args.force_scale
            return update_plot(scatter, scatter_sizes, bonds_plot_primary, bonds_plot_secondary, line_segments_primary, line_segments_secondary, force_lines, host_pos, f_bond, f_coll)
        
        return update_plot(scatter, scatter_sizes, bonds_plot_primary, bonds_plot_secondary, line_segments_primary, line_segments_secondary, force_lines, host_pos)
    
    frames_iter = itertools.count() if args.md_steps < 0 else range(args.md_steps)
    ani = FuncAnimation(fig, update, frames=frames_iter, interval=30, blit=False, repeat=False, save_count=0, cache_frame_data=False)
    plt.show()


if __name__ == "__main__":
    run_test()
