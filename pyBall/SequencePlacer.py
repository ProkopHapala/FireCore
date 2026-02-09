#!/usr/bin/env python3
"""
Interactive Matplotlib GUI for placing sequences of molecules on ionic substrate (NaCl step edges).

Features:
- NaCl simple-cubic step-edge substrate generator (parameterized)
- Load multiple molecule types, assign each to a letter
- Define a sequence string (e.g. "ABAB") to place molecules in a row
- Independent control of row direction (angle in surface plane) and molecule rotation (Rx,Ry,Rz)
- Row spacing along lattice direction
- Fast 3D scatter rendering with clip
- TextBox-based controls with mouse-wheel increment/decrement
- Export combined structure (substrate + molecules) to .xyz
- Uses AtomicSystem as underlying representation

Usage:
    python SequencePlacer.py                         # interactive GUI
    python SequencePlacer.py --headless ...           # headless export
"""

import numpy as np
import os, sys, json

# Ensure local pyBall package is importable when running script directly
_THIS_DIR = os.path.abspath(os.path.dirname(__file__))
_ROOT_DIR = os.path.abspath(os.path.join(_THIS_DIR, os.pardir))
for _p in (_ROOT_DIR, _THIS_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from pyBall import SubstrateBuilder as SB
try:
    from pyBall import atomicUtils as au
except Exception:
    au = None

# Default molecules (symbol -> relative path)
DEFAULT_MOLS = {
    'A': os.path.join(_ROOT_DIR, 'cpp', 'common_resources', 'xyz', 'PTCDA.xyz'),
    'B': os.path.join(_ROOT_DIR, 'cpp', 'common_resources', 'xyz', 'pentacene.xyz'),
}

# Element colors and sizes for visualization
ELEMENT_COLORS = {
    'H': '#FFFFFF', 'C': '#404040', 'N': '#0000FF', 'O': '#FF0000',
    'S': '#FFFF00', 'F': '#00FF00', 'Cl': '#00FF00', 'Br': '#8B0000',
    'I': '#9400D3', 'P': '#FFA500', 'B': '#FFB6C1', 'Si': '#DAA520',
    'Na': '#AB82FF', 'K': '#8B008B',
}
ELEMENT_SIZES = {
    'H': 15, 'C': 40, 'N': 40, 'O': 40, 'S': 50, 'F': 30, 'Cl': 45,
    'Br': 55, 'I': 70, 'P': 50, 'B': 35, 'Si': 50,
    'Na': 35, 'K': 45,
}

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

# ======================== Sequence placement logic ========================

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


def place_sequence(mol_dict, sequence, row_dir, row_spacing, mol_rotations,
                   origin=np.zeros(3), height=0.0):
    """Place a sequence of molecules along a row direction.
    
    Args:
        mol_dict:      dict letter -> (apos_centered, enames)  (centered at origin)
        sequence:      string of letters, e.g. "ABAB"
        row_dir:       2D unit vector [dx, dy] for row direction in surface plane
        row_spacing:   distance between consecutive molecule centers along row
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


# ======================== GUI ========================

class WheelTextBox:
    """Matplotlib TextBox with mouse-wheel increment/decrement."""
    def __init__(self, ax, label, initial, step=1.0, fmt="%.3f", on_change=None):
        from matplotlib.widgets import TextBox
        self.step = step
        self.fmt  = fmt
        self.val  = float(initial)
        self._on_change = on_change
        self.tb = TextBox(ax, label, initial=self.fmt % self.val)
        self.tb.on_submit(self._on_submit)
        ax.figure.canvas.mpl_connect('scroll_event', self._on_scroll)
        self._ax = ax

    def _on_submit(self, text):
        try:
            self.val = float(text)
        except ValueError:
            self.tb.set_val(self.fmt % self.val)
            return
        if self._on_change: self._on_change(self.val)

    def _on_scroll(self, event):
        if event.inaxes != self._ax: return
        if event.button == 'up':   self.val += self.step
        elif event.button == 'down': self.val -= self.step
        self.tb.set_val(self.fmt % self.val)
        if self._on_change: self._on_change(self.val)

    def set_val(self, v):
        self.val = float(v)
        self.tb.set_val(self.fmt % self.val)


class SequencePlacerGUI:
    """Interactive GUI for placing molecule sequences on NaCl step-edge substrate."""

    def __init__(self):
        import matplotlib
        import matplotlib.pyplot as plt
        from matplotlib.widgets import Button, TextBox
        from mpl_toolkits.mplot3d import Axes3D

        self.plt = plt
        self.Button = Button
        self.TextBox = TextBox

        # Data
        self.substrate = None        # AtomicSystem
        self.mol_dict  = {}          # letter -> (apos_centered, enames)
        self.mol_files = {}          # letter -> filename
        self.placed_apos   = None    # (N,3) placed molecule positions
        self.placed_enames = None

        # Substrate params
        self.sub_a  = 2.82065
        self.sub_nx = 15
        self.sub_ny = 8
        self.sub_nz = 3
        self.sub_nsteps = 1

        # Placement params
        self.row_angle   = 0.0     # degrees, angle of row in XY plane
        self.row_spacing = 5.0     # Angstrom between molecule centers
        self.mol_height  = 3.0     # z offset above substrate max z
        self.origin_x    = 0.0
        self.origin_y    = 0.0

        # Per-letter rotations (default 0,0,0)
        self.mol_rotations = {}

        # Sequence
        self.sequence = "A"

        # Molecule assignments string: "A=/path/to/mol.xyz,B=/path/to/other.xyz"
        self.mol_assign_str = ""

        # Rendering caches / artists
        self._sub_colors = None
        self._sub_sizes  = None
        self._mol_colors = None
        self._mol_sizes  = None
        self.sub_scatter3d = None
        self.sub_scatter2d = None
        self.mol_scatter3d = None
        self.mol_scatter2d = None
        self.bond_lines3d  = None
        self.bond_lines2d  = None
        self.mol_bonds = {}        # letter -> (nbonds,2) array
        self.placed_bonds = None

        self._setup_figure()
        # Auto-load defaults if available
        self._load_default_molecules()

    def _setup_figure(self):
        plt = self.plt
        self.fig = plt.figure(figsize=(18, 11))
        self.fig.suptitle('Sequence Placer - Molecules on NaCl Step Edge', fontsize=13)

        # ---- 3D view (main) ----
        self.ax3d = self.fig.add_axes([0.02, 0.38, 0.46, 0.57], projection='3d')
        self.ax3d.set_xlabel('X'); self.ax3d.set_ylabel('Y'); self.ax3d.set_zlabel('Z')

        # ---- Side view XZ ----
        self.ax_xz = self.fig.add_axes([0.52, 0.38, 0.46, 0.57])
        self.ax_xz.set_xlabel('X (A)'); self.ax_xz.set_ylabel('Z (A)')
        self.ax_xz.set_aspect('equal'); self.ax_xz.grid(True, alpha=0.3)

        # ---- Controls ----
        h  = 0.025   # row height
        g  = 0.005   # gap
        y0 = 0.30    # top of controls

        # Column 1: Substrate params
        c1 = 0.08; w1 = 0.14
        self.wb_a      = WheelTextBox(self.fig.add_axes([c1, y0,        w1, h]), 'a ',       self.sub_a,      step=0.1,  on_change=lambda v: self._set_sub('a', v))
        self.wb_nx     = WheelTextBox(self.fig.add_axes([c1, y0-1*(h+g), w1, h]), 'nx ',      self.sub_nx,     step=1, fmt="%.0f", on_change=lambda v: self._set_sub('nx', v))
        self.wb_ny     = WheelTextBox(self.fig.add_axes([c1, y0-2*(h+g), w1, h]), 'ny ',      self.sub_ny,     step=1, fmt="%.0f", on_change=lambda v: self._set_sub('ny', v))
        self.wb_nz     = WheelTextBox(self.fig.add_axes([c1, y0-3*(h+g), w1, h]), 'nz ',      self.sub_nz,     step=1, fmt="%.0f", on_change=lambda v: self._set_sub('nz', v))
        self.wb_nsteps = WheelTextBox(self.fig.add_axes([c1, y0-4*(h+g), w1, h]), 'steps ',   self.sub_nsteps,  step=1, fmt="%.0f", on_change=lambda v: self._set_sub('nsteps', v))

        # Column 2: Placement params
        c2 = 0.30; w2 = 0.14
        self.wb_row_ang  = WheelTextBox(self.fig.add_axes([c2, y0,        w2, h]), 'row ang ',  self.row_angle,   step=5.0,  on_change=lambda v: self._set_place('row_angle', v))
        self.wb_spacing  = WheelTextBox(self.fig.add_axes([c2, y0-1*(h+g), w2, h]), 'spacing ',  self.row_spacing,  step=0.5,  on_change=lambda v: self._set_place('row_spacing', v))
        self.wb_height   = WheelTextBox(self.fig.add_axes([c2, y0-2*(h+g), w2, h]), 'height ',   self.mol_height,   step=0.5,  on_change=lambda v: self._set_place('mol_height', v))
        self.wb_orig_x   = WheelTextBox(self.fig.add_axes([c2, y0-3*(h+g), w2, h]), 'orig X ',   self.origin_x,     step=1.0,  on_change=lambda v: self._set_place('origin_x', v))
        self.wb_orig_y   = WheelTextBox(self.fig.add_axes([c2, y0-4*(h+g), w2, h]), 'orig Y ',   self.origin_y,     step=1.0,  on_change=lambda v: self._set_place('origin_y', v))

        # Column 3: Molecule rotation (global, per-letter later)
        c3 = 0.52; w3 = 0.14
        self.wb_rx = WheelTextBox(self.fig.add_axes([c3, y0,        w3, h]), 'Rot X ',  0.0,  step=5.0,  on_change=lambda v: self._set_rot('rx', v))
        self.wb_ry = WheelTextBox(self.fig.add_axes([c3, y0-1*(h+g), w3, h]), 'Rot Y ',  0.0,  step=5.0,  on_change=lambda v: self._set_rot('ry', v))
        self.wb_rz = WheelTextBox(self.fig.add_axes([c3, y0-2*(h+g), w3, h]), 'Rot Z ',  0.0,  step=5.0,  on_change=lambda v: self._set_rot('rz', v))

        # Column 4: Molecule assignments and sequence (wide text boxes)
        c4 = 0.02; w4 = 0.50
        ax_mol = self.fig.add_axes([c4, 0.07, w4, h])
        self.tb_mol_assign = self.TextBox(ax_mol, 'Mols ', initial=self.mol_assign_str)
        self.tb_mol_assign.on_submit(self._on_mol_assign)

        ax_seq = self.fig.add_axes([c4, 0.035, w4, h])
        self.tb_sequence = self.TextBox(ax_seq, 'Seq ', initial=self.sequence)
        self.tb_sequence.on_submit(self._on_sequence)

        # Buttons
        bw = 0.09; bh = 0.035
        c5 = 0.75
        ax_build  = self.fig.add_axes([c5, y0,         bw, bh])
        ax_place  = self.fig.add_axes([c5, y0-1*(bh+g), bw, bh])
        ax_save   = self.fig.add_axes([c5, y0-2*(bh+g), bw, bh])
        ax_savem  = self.fig.add_axes([c5, y0-3*(bh+g), bw, bh])
        ax_reset  = self.fig.add_axes([c5, y0-4*(bh+g), bw, bh])

        self.btn_build = self.Button(ax_build, 'Build Sub')
        self.btn_place = self.Button(ax_place, 'Place Seq')
        self.btn_save  = self.Button(ax_save,  'Save XYZ')
        self.btn_savem = self.Button(ax_savem, 'Save Mols')
        self.btn_reset = self.Button(ax_reset, 'Reset View')

        self.btn_build.on_clicked(self._on_build_sub)
        self.btn_place.on_clicked(self._on_place_seq)
        self.btn_save.on_clicked(self._on_save_all)
        self.btn_savem.on_clicked(self._on_save_mols)
        self.btn_reset.on_clicked(self._on_reset_view)

        # Info / status
        self.status_text = self.fig.text(0.02, 0.005, 'Ready. Set molecule assignments (e.g. A=mol.xyz,B=other.xyz) then Build Sub / Place Seq.', fontsize=9)

    def _clear_mol_artists(self):
        if self.mol_scatter3d is not None:
            self.mol_scatter3d.remove(); self.mol_scatter3d = None
        if self.mol_scatter2d is not None:
            self.mol_scatter2d.remove(); self.mol_scatter2d = None
        if self.bond_lines3d is not None:
            self.bond_lines3d.remove(); self.bond_lines3d = None
        if self.bond_lines2d is not None:
            self.bond_lines2d.remove(); self.bond_lines2d = None

    def _load_default_molecules(self):
        """Load default molecules if files exist; populate assignment box."""
        available = {}
        parts = []
        for letter, path in DEFAULT_MOLS.items():
            if os.path.exists(path):
                available[letter] = path
                parts.append(f"{letter}={path}")
        if not available:
            return
        self.mol_assign_str = ",".join(parts)
        # Use programmatic loader to avoid widget submit jitter
        self.load_molecules_programmatic(available)
        self.tb_mol_assign.set_val(self.mol_assign_str)
        self.status_text.set_text(f"Loaded defaults: {list(available.keys())}")

    # ---- Param setters ----

    def _set_sub(self, key, val):
        if   key == 'a':      self.sub_a  = float(val)
        elif key == 'nx':     self.sub_nx = max(1, int(val))
        elif key == 'ny':     self.sub_ny = max(1, int(val))
        elif key == 'nz':     self.sub_nz = max(1, int(val))
        elif key == 'nsteps': self.sub_nsteps = max(1, int(val))

    def _set_place(self, key, val):
        if   key == 'row_angle':   self.row_angle   = float(val)
        elif key == 'row_spacing': self.row_spacing  = float(val)
        elif key == 'mol_height':  self.mol_height   = float(val)
        elif key == 'origin_x':    self.origin_x     = float(val)
        elif key == 'origin_y':    self.origin_y     = float(val)

    def _set_rot(self, key, val):
        """Set rotation globally for all loaded molecule types."""
        v = float(val)
        for letter in self.mol_dict:
            rx, ry, rz = self.mol_rotations.get(letter, (0., 0., 0.))
            if   key == 'rx': rx = v
            elif key == 'ry': ry = v
            elif key == 'rz': rz = v
            self.mol_rotations[letter] = (rx, ry, rz)

    # ---- Molecule assignment ----

    def _on_mol_assign(self, text):
        """Parse molecule assignment string like 'A=/path/mol.xyz,B=/path/other.xyz'"""
        self.mol_assign_str = text.strip()
        self.mol_dict.clear()
        self.mol_files.clear()
        self.mol_rotations.clear()
        self.mol_bonds.clear()
        self._clear_mol_artists()
        if not self.mol_assign_str:
            self.status_text.set_text('Molecule assignments cleared.')
            return
        parts = [p.strip() for p in self.mol_assign_str.split(',')]
        for part in parts:
            if '=' not in part:
                self.status_text.set_text(f'ERROR: bad format "{part}", expected LETTER=path')
                return
            letter, path = part.split('=', 1)
            letter = letter.strip()
            path   = path.strip()
            if len(letter) != 1:
                self.status_text.set_text(f'ERROR: key must be single letter, got "{letter}"')
                return
            if not os.path.exists(path):
                self.status_text.set_text(f'ERROR: file not found: {path}')
                return
            apos, enames = load_molecule_xyz(path)
            apos -= apos.mean(axis=0)  # center at origin
            self.mol_dict[letter]  = (apos, enames)
            self.mol_files[letter] = path
            self.mol_rotations[letter] = (self.wb_rx.val, self.wb_ry.val, self.wb_rz.val)
            if au is not None:
                bonds, _ = au.findBondsNP(apos, Rcut=1.8, byRvdW=False)
                self.mol_bonds[letter] = bonds
            else:
                self.mol_bonds[letter] = None
        self.status_text.set_text(f'Loaded {len(self.mol_dict)} molecule(s): {list(self.mol_dict.keys())}')
        self._update_mol_cache()

    def _on_sequence(self, text):
        self.sequence = text.strip()

    def load_molecules_programmatic(self, assignments):
        """Load molecules from a dict {letter: filepath}. For headless use."""
        for letter, path in assignments.items():
            assert os.path.exists(path), f"File not found: {path}"
            apos, enames = load_molecule_xyz(path)
            apos -= apos.mean(axis=0)
            self.mol_dict[letter]  = (apos, enames)
            self.mol_files[letter] = path
            self.mol_rotations[letter] = (0., 0., 0.)
            # detect bonds (distance-based) if atomicUtils available
            if au is not None:
                bonds, _ = au.findBondsNP(apos, Rcut=1.8, byRvdW=False)
                self.mol_bonds[letter] = bonds
            else:
                self.mol_bonds[letter] = None
        self._update_mol_cache()
        self._clear_mol_artists()

    # ---- Build / Place ----

    def _on_build_sub(self, event=None):
        self.build_substrate()
        self._update_plot()

    def build_substrate(self):
        self.substrate = SB.gen_nacl_step(
            a=self.sub_a, nx=self.sub_nx, ny=self.sub_ny, nz=self.sub_nz,
            nsteps=self.sub_nsteps)
        self._update_sub_cache()
        # Reset substrate artists so they are recreated
        if self.sub_scatter3d is not None:
            self.sub_scatter3d.remove(); self.sub_scatter3d = None
        if self.sub_scatter2d is not None:
            self.sub_scatter2d.remove(); self.sub_scatter2d = None
        self.status_text.set_text(f'Built substrate: {len(self.substrate.apos)} atoms, {self.sub_nsteps} step(s)')

    def _on_place_seq(self, event=None):
        self.place_molecules()
        self._update_plot()

    def place_molecules(self):
        if not self.mol_dict:
            raise ValueError("No molecules loaded. Set molecule assignments first.")
        for ch in self.sequence:
            if ch not in self.mol_dict:
                raise ValueError(f"Letter '{ch}' in sequence not loaded. Available: {list(self.mol_dict.keys())}")
        ang_rad = np.radians(self.row_angle)
        row_dir = np.array([np.cos(ang_rad), np.sin(ang_rad)])

        # Determine height: substrate max z + mol_height
        if self.substrate is not None:
            z_top = self.substrate.apos[:, 2].max()
        else:
            z_top = 0.0

        origin = np.array([self.origin_x, self.origin_y, 0.0])

        # regenerate placed atoms
        self.placed_apos, self.placed_enames = place_sequence(
            self.mol_dict, self.sequence, row_dir, self.row_spacing,
            self.mol_rotations, origin=origin, height=z_top + self.mol_height)

        # replicate bonds
        self.placed_bonds = None
        bonds_all = []
        offset = 0
        for ch in self.sequence:
            apos_ch, enames_ch = self.mol_dict[ch]
            nb = len(enames_ch)
            bonds = self.mol_bonds.get(ch)
            if bonds is not None and len(bonds) > 0:
                bonds_all.append(bonds + offset)
            offset += nb
        if bonds_all:
            import numpy as _np
            self.placed_bonds = _np.vstack(bonds_all)
        self._update_mol_cache()
        # Reset molecule artists to force refresh
        self._clear_mol_artists()
        self.status_text.set_text(f'Placed {len(self.sequence)} molecule(s), {len(self.placed_apos)} atoms')
        # In non-GUI usage, ensure plot updated
        self._update_plot()

    # ---- Rendering caches ----

    def _update_sub_cache(self):
        if self.substrate is None: return
        enames = self.substrate.enames
        self._sub_colors = [ELEMENT_COLORS.get(e, '#808080') for e in enames]
        self._sub_sizes  = np.array([ELEMENT_SIZES.get(e, 30) for e in enames])

    def _update_mol_cache(self):
        if self.placed_apos is None: return
        self._mol_colors = [ELEMENT_COLORS.get(e, '#808080') for e in self.placed_enames]
        self._mol_sizes  = np.array([ELEMENT_SIZES.get(e, 30) for e in self.placed_enames])

    # ---- Plotting ----

    def _update_plot(self):
        # Substrate artists (draw once)
        if self.substrate is not None and self._sub_colors is not None:
            sp = self.substrate.apos
            if self.sub_scatter3d is None:
                self.sub_scatter3d = self.ax3d.scatter(sp[:,0], sp[:,1], sp[:,2], c=self._sub_colors, s=self._sub_sizes, alpha=0.35, depthshade=True, edgecolors='none')
            if self.sub_scatter2d is None:
                self.sub_scatter2d = self.ax_xz.scatter(sp[:,0], sp[:,2], c=self._sub_colors, s=self._sub_sizes, alpha=0.35, edgecolors='none')

        # Molecule scatters (update offsets)
        if self.placed_apos is not None and self._mol_colors is not None:
            mp = self.placed_apos
            xs, ys, zs = mp[:,0], mp[:,1], mp[:,2]
            if self.mol_scatter3d is None:
                self.mol_scatter3d = self.ax3d.scatter(xs, ys, zs, c=self._mol_colors, s=self._mol_sizes, alpha=0.85, depthshade=True, edgecolors='none')
            else:
                self.mol_scatter3d._offsets3d = (xs, ys, zs)
            if self.mol_scatter2d is None:
                self.mol_scatter2d = self.ax_xz.scatter(xs, zs, c=self._mol_colors, s=self._mol_sizes, alpha=0.85, edgecolors='none')
            else:
                self.mol_scatter2d.set_offsets(mp[:,[0,2]])

            # Bonds
            if self.placed_bonds is not None and len(self.placed_bonds)>0:
                from mpl_toolkits.mplot3d.art3d import Line3DCollection
                from matplotlib.collections import LineCollection
                seg3d = [ [mp[i], mp[j]] for i,j in self.placed_bonds ]
                seg2d = [ [[mp[i,0], mp[i,2]],[mp[j,0], mp[j,2]]] for i,j in self.placed_bonds ]
                if self.bond_lines3d is None:
                    self.bond_lines3d = Line3DCollection(seg3d, colors='gray', linewidths=0.6, alpha=0.8)
                    self.ax3d.add_collection3d(self.bond_lines3d)
                else:
                    self.bond_lines3d.set_segments(seg3d)
                if self.bond_lines2d is None:
                    self.bond_lines2d = LineCollection(seg2d, colors='gray', linewidths=0.6, alpha=0.8)
                    self.ax_xz.add_collection(self.bond_lines2d)
                else:
                    self.bond_lines2d.set_segments(seg2d)

        # Axes labels and aspect once
        self.ax3d.set_xlabel('X'); self.ax3d.set_ylabel('Y'); self.ax3d.set_zlabel('Z')
        self.ax_xz.set_xlabel('X (A)'); self.ax_xz.set_ylabel('Z (A)')
        self.ax_xz.set_aspect('equal'); self.ax_xz.grid(True, alpha=0.3)

        # Auto-scale XZ
        all_x, all_z = [], []
        if self.substrate is not None:
            all_x.append(self.substrate.apos[:,0]); all_z.append(self.substrate.apos[:,2])
        if self.placed_apos is not None:
            all_x.append(self.placed_apos[:,0]); all_z.append(self.placed_apos[:,2])
        if all_x:
            ax = np.concatenate(all_x); az = np.concatenate(all_z)
            m = 3
            self.ax_xz.set_xlim(ax.min()-m, ax.max()+m)
            self.ax_xz.set_ylim(az.min()-m, az.max()+m)

        self.fig.canvas.draw_idle()

    def _on_save_all(self, event=None):
        """Save combined substrate + molecules to xyz."""
        all_apos = []; all_enames = []; lvec = None
        if self.substrate is not None:
            all_apos.append(self.substrate.apos)
            all_enames.extend(self.substrate.enames)
            lvec = self.substrate.lvec
        if self.placed_apos is not None:
            all_apos.append(self.placed_apos)
            all_enames.extend(self.placed_enames)
        if not all_apos:
            self.status_text.set_text('Nothing to save')
            return
        apos = np.vstack(all_apos)
        save_xyz('sequence_on_substrate.xyz', apos, all_enames, lvec=lvec, comment=" sequence=" + self.sequence)
        self.status_text.set_text(f'Saved sequence_on_substrate.xyz ({len(apos)} atoms)')

    def _on_save_mols(self, event=None):
        """Save only placed molecules to xyz."""
        if self.placed_apos is None:
            self.status_text.set_text('No molecules placed')
            return
        save_xyz('placed_molecules.xyz', self.placed_apos, self.placed_enames, comment=" sequence=" + self.sequence)
        self.status_text.set_text(f'Saved placed_molecules.xyz ({len(self.placed_apos)} atoms)')

    def _on_reset_view(self, event=None):
        self.ax3d.view_init(elev=25, azim=-60)
        self.fig.canvas.draw_idle()

    def save_combined(self, fname):
        """Programmatic save of combined structure."""
        all_apos = []; all_enames = []; lvec = None
        if self.substrate is not None:
            all_apos.append(self.substrate.apos)
            all_enames.extend(list(self.substrate.enames))
            lvec = self.substrate.lvec
        if self.placed_apos is not None:
            all_apos.append(self.placed_apos)
            all_enames.extend(self.placed_enames)
        assert len(all_apos) > 0, "Nothing to save"
        apos = np.vstack(all_apos)
        save_xyz(fname, apos, all_enames, lvec=lvec, comment=" sequence=" + self.sequence)

    def save_molecules_only(self, fname):
        """Programmatic save of placed molecules only."""
        assert self.placed_apos is not None, "No molecules placed"
        save_xyz(fname, self.placed_apos, self.placed_enames, comment=" sequence=" + self.sequence)

    def run(self):
        self.plt.show()


# ======================== Headless CLI ========================

def run_headless(mol_files, sequence="A", a=2.82065, nx=15, ny=8, nz=3, nsteps=1,
                 row_angle=0.0, row_spacing=5.0, mol_height=3.0,
                 origin_x=0.0, origin_y=0.0,
                 rx=0.0, ry=0.0, rz=0.0,
                 out_combined="sequence_on_substrate.xyz",
                 out_mols="placed_molecules.xyz"):
    """Run placement headlessly and save results.
    
    Args:
        mol_files: dict {letter: filepath}
        sequence:  e.g. "ABAB"
        (other args as in GUI)
    """
    # Build substrate
    substrate = SB.gen_nacl_step(a=a, nx=nx, ny=ny, nz=nz, nsteps=nsteps)

    # Load molecules
    mol_dict = {}
    for letter, path in mol_files.items():
        assert os.path.exists(path), f"File not found: {path}"
        apos, enames = load_molecule_xyz(path)
        apos -= apos.mean(axis=0)
        mol_dict[letter] = (apos, enames)

    # Rotations
    mol_rotations = {letter: (rx, ry, rz) for letter in mol_dict}

    # Place
    ang_rad = np.radians(row_angle)
    row_dir = np.array([np.cos(ang_rad), np.sin(ang_rad)])
    z_top   = substrate.apos[:, 2].max()
    origin  = np.array([origin_x, origin_y, 0.0])

    placed_apos, placed_enames = place_sequence(
        mol_dict, sequence, row_dir, row_spacing,
        mol_rotations, origin=origin, height=z_top + mol_height)

    # Save combined
    all_apos   = np.vstack([substrate.apos, placed_apos])
    all_enames = list(substrate.enames) + placed_enames
    save_xyz(out_combined, all_apos, all_enames, lvec=substrate.lvec, comment=" sequence=" + sequence)

    # Save molecules only
    save_xyz(out_mols, placed_apos, placed_enames, comment=" sequence=" + sequence)

    return substrate, placed_apos, placed_enames


# ======================== Main ========================

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Sequence Placer - molecules on NaCl step edge')
    parser.add_argument('--headless', action='store_true', help='Run without GUI')
    parser.add_argument('--mols', type=str, default='', help='Molecule assignments: A=file.xyz,B=file2.xyz')
    parser.add_argument('--seq',  type=str, default='A', help='Sequence string')
    parser.add_argument('--a',    type=float, default=2.82065)
    parser.add_argument('--nx',   type=int, default=15)
    parser.add_argument('--ny',   type=int, default=8)
    parser.add_argument('--nz',   type=int, default=3)
    parser.add_argument('--nsteps', type=int, default=1)
    parser.add_argument('--row_angle', type=float, default=0.0)
    parser.add_argument('--row_spacing', type=float, default=12.0)
    parser.add_argument('--mol_height', type=float, default=3.0)
    parser.add_argument('--rx', type=float, default=0.0)
    parser.add_argument('--ry', type=float, default=0.0)
    parser.add_argument('--rz', type=float, default=0.0)
    parser.add_argument('--origin_x', type=float, default=0.0)
    parser.add_argument('--origin_y', type=float, default=0.0)
    parser.add_argument('--out', type=str, default='sequence_on_substrate.xyz')
    parser.add_argument('--out_mols', type=str, default='placed_molecules.xyz')
    args = parser.parse_args()

    if args.headless:
        # Parse mol assignments
        mol_files = {}
        if args.mols:
            for part in args.mols.split(','):
                letter, path = part.strip().split('=', 1)
                mol_files[letter.strip()] = path.strip()
        assert mol_files, "No molecules specified (use --mols A=file.xyz,...)"
        run_headless(mol_files, sequence=args.seq, a=args.a, nx=args.nx, ny=args.ny, nz=args.nz,
                     nsteps=args.nsteps, row_angle=args.row_angle, row_spacing=args.row_spacing,
                     mol_height=args.mol_height, rx=args.rx, ry=args.ry, rz=args.rz,
                     origin_x=args.origin_x, origin_y=args.origin_y,
                     out_combined=args.out, out_mols=args.out_mols)
    else:
        gui = SequencePlacerGUI()
        gui.build_substrate()
        gui._update_plot()
        gui.run()
