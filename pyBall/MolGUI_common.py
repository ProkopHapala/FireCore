#!/usr/bin/env python3
"""
Shared utilities for molecular placement GUIs (both Matplotlib and VisPy versions).
Provides: element visuals, rotation matrices, xyz I/O, bond detection, sequence placement, hex lattice.
"""

import numpy as np
import os

# ======================== Element visuals ========================

ELEMENT_COLORS = {
    'H': '#C8C8C8', 'C': '#404040', 'N': '#0000FF', 'O': '#FF0000',
    'S': '#FFFF00', 'F': '#00FF00', 'Cl': '#00FF00', 'Br': '#8B0000',
    'I': '#9400D3', 'P': '#FFA500', 'B': '#FFB6C1', 'Si': '#DAA520',
    'Na': '#AB82FF', 'K': '#8B008B',
}
ELEMENT_SIZES = {
    'H': 15, 'C': 40, 'N': 40, 'O': 40, 'S': 50, 'F': 30, 'Cl': 45,
    'Br': 55, 'I': 70, 'P': 50, 'B': 35, 'Si': 50,
    'Na': 35, 'K': 45,
}

def hex_to_rgba(hex_str, alpha=1.0):
    """Convert '#RRGGBB' to (r,g,b,a) floats in [0,1]."""
    h = hex_str.lstrip('#')
    return (int(h[0:2],16)/255., int(h[2:4],16)/255., int(h[4:6],16)/255., alpha)

def colors_for_enames(enames, alpha=1.0):
    """Return (N,4) float32 RGBA array for element names."""
    out = np.empty((len(enames), 4), dtype=np.float32)
    for i, e in enumerate(enames):
        out[i] = hex_to_rgba(ELEMENT_COLORS.get(e, '#808080'), alpha)
    return out

def sizes_for_enames(enames, scale=1.0):
    """Return (N,) float32 array of marker sizes."""
    return np.array([ELEMENT_SIZES.get(e, 30)*scale for e in enames], dtype=np.float32)

# ======================== Rotation helpers ========================

def rotmat_x(a_deg):
    a = np.radians(a_deg); c, s = np.cos(a), np.sin(a)
    return np.array([[1,0,0],[0,c,-s],[0,s,c]])

def rotmat_y(a_deg):
    a = np.radians(a_deg); c, s = np.cos(a), np.sin(a)
    return np.array([[c,0,s],[0,1,0],[-s,0,c]])

def rotmat_z(a_deg):
    a = np.radians(a_deg); c, s = np.cos(a), np.sin(a)
    return np.array([[c,-s,0],[s,c,0],[0,0,1]])

def make_rotmat(rx, ry, rz):
    return rotmat_z(rz) @ rotmat_y(ry) @ rotmat_x(rx)

def flip_matrix_x():
    return np.array([[-1,0,0],[0,1,0],[0,0,1]], dtype=float)

def flip_matrix_y():
    return np.array([[1,0,0],[0,-1,0],[0,0,1]], dtype=float)

# ======================== XYZ I/O ========================

def load_molecule_xyz(fname):
    """Load molecule from xyz file, return (apos, enames). Minimal, no dependencies."""
    with open(fname, 'r') as f:
        lines = f.readlines()
    natoms = int(lines[0].strip())
    apos = np.empty((natoms, 3))
    enames = []
    for i in range(natoms):
        parts = lines[2 + i].split()
        enames.append(parts[0])
        apos[i] = [float(parts[1]), float(parts[2]), float(parts[3])]
    return apos, enames

def save_xyz(fname, apos, enames, lvec=None, comment=""):
    """Save .xyz with optional lvec in comment."""
    if lvec is not None:
        lv = lvec
        comment = ("lvs %g %g %g  %g %g %g  %g %g %g " % (
            lv[0,0],lv[0,1],lv[0,2], lv[1,0],lv[1,1],lv[1,2], lv[2,0],lv[2,1],lv[2,2])) + comment
    with open(fname, 'w') as f:
        f.write(f"{len(apos)}\n{comment}\n")
        for e, p in zip(enames, apos):
            f.write(f"{e} {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")
    print(f"Saved {fname} ({len(apos)} atoms)")

# ======================== Bond detection ========================

try:
    from pyBall import atomicUtils as au
    _HAS_AU = True
except Exception:
    try:
        import atomicUtils as au
        _HAS_AU = True
    except Exception:
        _HAS_AU = False
        au = None

def find_bonds(apos, enames=None, Rcut=1.8, RvdwCut=1.5, byRvdW=True):
    """Return (nbonds,2) int array of bond pairs, or None."""
    if not _HAS_AU or au is None:
        return None
    try:
        atypes = None if enames is None else [au.elements.ELEMENT_DICT[e][0] for e in enames]
        bonds, _ = au.findBondsNP(apos, atypes=atypes, Rcut=Rcut, RvdwCut=RvdwCut, byRvdW=byRvdW)
        return bonds
    except Exception as e:
        print(f"Warning: bond detection failed: {e}")
        return None

def replicate_bonds(mol_bonds_dict, sequence, mol_dict):
    """Replicate per-letter bonds for a placed sequence. Returns (nbonds,2) or None."""
    bonds_all = []
    offset = 0
    for ch in sequence:
        apos_ch, enames_ch = mol_dict[ch]
        nb = len(enames_ch)
        bonds = mol_bonds_dict.get(ch)
        if bonds is not None and len(bonds) > 0:
            bonds_all.append(bonds + offset)
        offset += nb
    if bonds_all:
        return np.vstack(bonds_all)
    return None

# ======================== Sequence placement ========================

def place_sequence(mol_dict, sequence, row_dir, row_spacing, mol_rotations,
                   origin=np.zeros(3), height=0.0):
    """Place a sequence of molecules along a row direction.
    Args:
        mol_dict:      dict letter -> (apos_centered, enames)
        sequence:      string of letters, e.g. "ABAB"
        row_dir:       2D unit vector [dx, dy]
        row_spacing:   distance between consecutive molecule centers
        mol_rotations: dict letter -> (rx, ry, rz) in degrees
        origin:        3D starting point
        height:        z-offset above substrate surface
    Returns:
        all_apos (N,3), all_enames (list)
    """
    assert len(sequence) > 0, "place_sequence: empty sequence"
    dir3 = np.array([row_dir[0], row_dir[1], 0.0])
    norm = np.linalg.norm(dir3[:2])
    if norm > 1e-12: dir3 /= norm
    all_apos   = []
    all_enames = []
    for idx, letter in enumerate(sequence):
        if letter not in mol_dict:
            raise ValueError(f"place_sequence: letter '{letter}' not in mol_dict (keys: {list(mol_dict.keys())})")
        apos_orig, enames = mol_dict[letter]
        rx, ry, rz = mol_rotations.get(letter, (0., 0., 0.))
        R = make_rotmat(rx, ry, rz)
        apos_rot = (R @ apos_orig.T).T
        shift = origin + dir3 * (idx * row_spacing) + np.array([0, 0, height])
        all_apos.append(apos_rot + shift[None, :])
        all_enames.extend(enames)
    return np.vstack(all_apos), all_enames

# ======================== Hex lattice ========================

class HexLattice:
    """Hexagonal/triangular lattice utilities."""
    def __init__(self, a=32.7):
        self.a = a
        self._update_vectors()

    def _update_vectors(self):
        a = self.a
        self.a1 = np.array([a, 0, 0])
        self.a2 = np.array([a * 0.5, a * np.sqrt(3) / 2, 0])
        self.tri_up   = (self.a1 + self.a2) / 3
        self.tri_down = 2 * (self.a1 + self.a2) / 3

    def set_a(self, a):
        self.a = a; self._update_vectors()

    def cell_origin(self, i, j):
        return i * self.a1 + j * self.a2

    def triangle_center(self, i, j, up=True):
        return self.cell_origin(i, j) + (self.tri_up if up else self.tri_down)

# ======================== Default molecules ========================

_THIS_DIR = os.path.abspath(os.path.dirname(__file__))
_ROOT_DIR = os.path.abspath(os.path.join(_THIS_DIR, os.pardir))

DEFAULT_MOLS = {
    'A': os.path.join(_ROOT_DIR, 'cpp', 'common_resources', 'xyz', 'PTCDA.xyz'),
    'B': os.path.join(_ROOT_DIR, 'cpp', 'common_resources', 'xyz', 'pentacene.xyz'),
}

def load_default_molecules():
    """Load default molecules, return (mol_dict, mol_bonds, mol_files) or empty dicts."""
    mol_dict  = {}
    mol_bonds = {}
    mol_files = {}
    for letter, path in DEFAULT_MOLS.items():
        if os.path.exists(path):
            apos, enames = load_molecule_xyz(path)
            apos -= apos.mean(axis=0)
            mol_dict[letter]  = (apos, enames)
            mol_files[letter] = path
            mol_bonds[letter] = find_bonds(apos)
    return mol_dict, mol_bonds, mol_files
