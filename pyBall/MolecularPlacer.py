#!/usr/bin/env python3
"""
Interactive Matplotlib GUI for molecular placement with rotation and symmetry operations.
Designed for placing molecules into symmetric geometries (e.g., 3-fold, 6-fold symmetry).

Features:
- Load molecules from xyz files
- Rotate molecule with 3 rotation angles (Rx, Ry, Rz)
- N-fold rotational duplication around pivot
- Triangle-based positioning using lattice coordinates
- Dual-flower mode for unit cell with two triangles (inverted/rotated)
- Fast rendering with batched numpy operations
- Export final arrangement to xyz

Usage:
    python MolecularPlacer.py [xyz_file]
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, TextBox, RadioButtons
import os
import sys
import json

# Ensure local pyBall package is importable when running script directly
_THIS_DIR = os.path.abspath(os.path.dirname(__file__))
_ROOT_DIR = os.path.abspath(os.path.join(_THIS_DIR, os.pardir))
for _p in (_ROOT_DIR, _THIS_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

# Try to import atomicUtils for MOL2 export with bonds
try:
    from pyBall import atomicUtils as au
    HAS_ATOMICUTILS = True
except ImportError:
    try:
        import atomicUtils as au
        HAS_ATOMICUTILS = True
    except ImportError:
        HAS_ATOMICUTILS = False
        print("Warning: atomicUtils not found, MOL2 export disabled")

# Element colors for visualization
ELEMENT_COLORS = {
    'H': '#FFFFFF', 'C': '#404040', 'N': '#0000FF', 'O': '#FF0000',
    'S': '#FFFF00', 'F': '#00FF00', 'Cl': '#00FF00', 'Br': '#8B0000',
    'I': '#9400D3', 'P': '#FFA500', 'B': '#FFB6C1', 'Si': '#DAA520',
}
ELEMENT_SIZES = {
    'H': 20, 'C': 50, 'N': 50, 'O': 50, 'S': 60, 'F': 40, 'Cl': 55,
    'Br': 65, 'I': 80, 'P': 60, 'B': 45, 'Si': 60,
}

def rotation_matrix_x(angle_deg):
    """Rotation matrix around X axis."""
    a = np.radians(angle_deg)
    c, s = np.cos(a), np.sin(a)
    return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])

def rotation_matrix_y(angle_deg):
    """Rotation matrix around Y axis."""
    a = np.radians(angle_deg)
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])

def rotation_matrix_z(angle_deg):
    """Rotation matrix around Z axis."""
    a = np.radians(angle_deg)
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])

def flip_matrix_x():
    """Mirror/flip matrix across YZ plane (invert X)."""
    return np.array([[-1, 0, 0], [0, 1, 0], [0, 0, 1]])

def flip_matrix_y():
    """Mirror/flip matrix across XZ plane (invert Y)."""
    return np.array([[1, 0, 0], [0, -1, 0], [0, 0, 1]])


class HexLattice:
    """Hexagonal/triangular lattice utilities."""
    
    def __init__(self, a=32.7):
        self.a = a
        self._update_vectors()
    
    def _update_vectors(self):
        a = self.a
        self.a1 = np.array([a, 0, 0])                           # lattice vector 1
        self.a2 = np.array([a * 0.5, a * np.sqrt(3) / 2, 0])    # lattice vector 2
        # Triangle centers within unit cell (barycentric -> cartesian)
        self.tri_up   = (self.a1 + self.a2) / 3       # up-pointing triangle center
        self.tri_down = 2 * (self.a1 + self.a2) / 3   # down-pointing triangle center
    
    def set_a(self, a):
        self.a = a
        self._update_vectors()
    
    def cell_origin(self, i, j):
        """Get origin of cell (i,j)."""
        return i * self.a1 + j * self.a2
    
    def triangle_center(self, i, j, up=True):
        """Get center of up or down triangle in cell (i,j)."""
        origin = self.cell_origin(i, j)
        return origin + (self.tri_up if up else self.tri_down)
    
    def barycentric_to_cartesian(self, ba, bb, up=True):
        """Convert barycentric-like offsets (ba, bb) to cartesian position.
        ba, bb are fractions along a1, a2 directions from triangle center.
        """
        base = self.tri_up[:2] if up else self.tri_down[:2]
        offset = ba * self.a1[:2] + bb * self.a2[:2]
        return base + offset

def load_xyz_simple(fname):
    """Load xyz file without dependencies."""
    with open(fname, 'r') as f:
        lines = f.readlines()
    natoms = int(lines[0].strip())
    comment = lines[1].strip() if len(lines) > 1 else ""
    apos = []
    enames = []
    for i in range(2, 2 + natoms):
        parts = lines[i].split()
        enames.append(parts[0])
        apos.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return np.array(apos), enames, comment

def save_xyz_simple(fname, apos, enames, comment=""):
    """Save xyz file."""
    with open(fname, 'w') as f:
        f.write(f"{len(apos)}\n")
        f.write(f"{comment}\n")
        for i, (e, p) in enumerate(zip(enames, apos)):
            f.write(f"{e} {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")

def save_mol2_with_bonds(fname, apos, enames, bonds_per_mol, n_atoms_per_mol, comment=""):
    """Save MOL2 file with explicit bonding topology.
    bonds_per_mol: bonds for single molecule (0-indexed)
    n_atoms_per_mol: number of atoms in single molecule
    """
    n_mols = len(apos) // n_atoms_per_mol
    # Replicate bonds for each molecule copy
    all_bonds = []
    for imol in range(n_mols):
        offset = imol * n_atoms_per_mol
        for b in bonds_per_mol:
            all_bonds.append((b[0] + offset, b[1] + offset))
    all_bonds = np.array(all_bonds, dtype=np.int32) if all_bonds else None
    if HAS_ATOMICUTILS:
        au.save_mol2(fname, enames, apos, all_bonds, comment=comment)
    else:
        raise RuntimeError("atomicUtils not available for MOL2 export")

def save_poses(fname, poses, comment=""):
    """Save rigid-body poses to file.
    Each pose is (position[3], rot_matrix[3x3]) -> 12 floats per molecule.
    Format: pos_x pos_y pos_z  a_x a_y a_z  b_x b_y b_z  c_x c_y c_z
    where a,b,c are the columns of the rotation matrix.
    """
    with open(fname, 'w') as f:
        f.write(f"# {comment}\n")
        f.write(f"# pos_x pos_y pos_z  a_x a_y a_z  b_x b_y b_z  c_x c_y c_z\n")
        for pos, R in poses:
            f.write(f"{pos[0]:12.6f} {pos[1]:12.6f} {pos[2]:12.6f}  ")
            f.write(f"{R[0,0]:9.6f} {R[1,0]:9.6f} {R[2,0]:9.6f}  ")
            f.write(f"{R[0,1]:9.6f} {R[1,1]:9.6f} {R[2,1]:9.6f}  ")
            f.write(f"{R[0,2]:9.6f} {R[1,2]:9.6f} {R[2,2]:9.6f}\n")

class MolecularPlacer:
    """Interactive GUI for molecular placement and symmetry operations."""
    
    def __init__(self, xyz_file=None):
        self.apos_orig = None      # Original atomic positions
        self.enames = None         # Element names
        self.apos = None           # Current transformed positions
        
        # Rotation angles (degrees)
        self.rot_x = 0.0
        self.rot_y = 0.0
        self.rot_z = 0.0
        
        # Symmetry settings
        self.n_fold = 3            # n-fold symmetry
        self.radial_offset = 10.0  # distance from pivot to molecule center
        self.initial_angle = 0.0   # initial rotation angle for first copy
        
        # Lattice and triangle positioning
        self.lattice = HexLattice(a=32.7)
        self.tri_i = 0             # cell index i
        self.tri_j = 0             # cell index j  
        self.tri_ba = 0.0          # barycentric offset along a1
        self.tri_bb = 0.0          # barycentric offset along a2
        self.use_tri_up = True     # use up-triangle (True) or down-triangle (False)
        
        # Dual-flower mode (second flower in unit cell)
        self.dual_flower = False   # enable second flower
        self.dual_rot = 180.0      # rotation of second flower (degrees)
        self.dual_flip_x = False   # flip second flower in X
        self.dual_flip_y = False   # flip second flower in Y
        
        # View settings
        self.view_mode = 'xy'
        self.show_hex_grid = True
        
        # Precomputed for fast rendering
        self._colors_cache = None
        self._sizes_cache = None
        
        # Bond topology (computed on load)
        self.bonds = None
        
        self.setup_figure()
        
        if xyz_file and os.path.exists(xyz_file):
            self.load_molecule(xyz_file)
    
    def setup_figure(self):
        """Create the matplotlib figure and widgets."""
        self.fig = plt.figure(figsize=(16, 11))
        self.fig.suptitle('Molecular Placer - Interactive Symmetry Tool', fontsize=14)
        
        # Main plot area (top view XY)
        self.ax_main = self.fig.add_axes([0.03, 0.38, 0.50, 0.55])
        self.ax_main.set_aspect('equal')
        self.ax_main.set_xlabel('X (Å)')
        self.ax_main.set_ylabel('Y (Å)')
        self.ax_main.grid(True, alpha=0.3)
        
        # Side view (XZ)
        self.ax_side = self.fig.add_axes([0.56, 0.38, 0.42, 0.55])
        self.ax_side.set_aspect('equal')
        self.ax_side.set_xlabel('X (Å)')
        self.ax_side.set_ylabel('Z (Å)')
        self.ax_side.grid(True, alpha=0.3)
        
        # === LEFT COLUMN: Molecule rotation ===
        col1_x = 0.08
        ax_rx = self.fig.add_axes([col1_x, 0.30, 0.22, 0.025])
        ax_ry = self.fig.add_axes([col1_x, 0.265, 0.22, 0.025])
        ax_rz = self.fig.add_axes([col1_x, 0.23, 0.22, 0.025])
        
        self.slider_rx = Slider(ax_rx, 'Rot X', -180, 180, valinit=0)
        self.slider_ry = Slider(ax_ry, 'Rot Y', -180, 180, valinit=0)
        self.slider_rz = Slider(ax_rz, 'Rot Z', -180, 180, valinit=0)
        
        self.slider_rx.on_changed(self.on_rotation_change)
        self.slider_ry.on_changed(self.on_rotation_change)
        self.slider_rz.on_changed(self.on_rotation_change)
        
        # Flower symmetry controls
        ax_nfold   = self.fig.add_axes([col1_x, 0.18, 0.22, 0.025])
        ax_radius  = self.fig.add_axes([col1_x, 0.145, 0.22, 0.025])
        ax_initang = self.fig.add_axes([col1_x, 0.11, 0.22, 0.025])
        
        self.slider_nfold   = Slider(ax_nfold, 'N-fold', 1, 12, valinit=3, valstep=1)
        self.slider_radius  = Slider(ax_radius, 'Radius', 0, 50, valinit=10)
        self.slider_initang = Slider(ax_initang, 'Init φ', -180, 180, valinit=0)
        
        self.slider_nfold.on_changed(self.on_symmetry_change)
        self.slider_radius.on_changed(self.on_symmetry_change)
        self.slider_initang.on_changed(self.on_symmetry_change)
        
        # === MIDDLE COLUMN: Triangle/lattice positioning ===
        col2_x = 0.40
        ax_lat   = self.fig.add_axes([col2_x, 0.30, 0.18, 0.025])
        ax_tri_i = self.fig.add_axes([col2_x, 0.265, 0.18, 0.025])
        ax_tri_j = self.fig.add_axes([col2_x, 0.23, 0.18, 0.025])
        
        self.slider_lattice = Slider(ax_lat, 'Lattice a', 10, 60, valinit=32.7)
        self.slider_tri_i   = Slider(ax_tri_i, 'Cell i', -2, 2, valinit=0, valstep=1)
        self.slider_tri_j   = Slider(ax_tri_j, 'Cell j', -2, 2, valinit=0, valstep=1)
        
        self.slider_lattice.on_changed(self.on_lattice_change)
        self.slider_tri_i.on_changed(self.on_triangle_change)
        self.slider_tri_j.on_changed(self.on_triangle_change)
        
        # Fine position within cell (barycentric-like)
        ax_tri_ba = self.fig.add_axes([col2_x, 0.18, 0.18, 0.025])
        ax_tri_bb = self.fig.add_axes([col2_x, 0.145, 0.18, 0.025])
        
        self.slider_tri_ba = Slider(ax_tri_ba, 'Offset a', -0.5, 0.5, valinit=0)
        self.slider_tri_bb = Slider(ax_tri_bb, 'Offset b', -0.5, 0.5, valinit=0)
        
        self.slider_tri_ba.on_changed(self.on_triangle_change)
        self.slider_tri_bb.on_changed(self.on_triangle_change)
        
        # Triangle type toggle
        ax_tri_type = self.fig.add_axes([col2_x, 0.08, 0.08, 0.05])
        self.btn_tri_type = Button(ax_tri_type, '▲ Up')
        self.btn_tri_type.on_clicked(self.on_tri_type_click)
        
        # === RIGHT COLUMN: Dual flower controls ===
        col3_x = 0.68
        ax_dual_toggle = self.fig.add_axes([col3_x, 0.30, 0.10, 0.03])
        self.btn_dual = Button(ax_dual_toggle, 'Dual: OFF')
        self.btn_dual.on_clicked(self.on_dual_toggle)
        
        ax_dual_rot = self.fig.add_axes([col3_x, 0.255, 0.18, 0.025])
        self.slider_dual_rot = Slider(ax_dual_rot, 'Dual φ', -180, 180, valinit=180)
        self.slider_dual_rot.on_changed(self.on_dual_change)
        
        ax_dual_flip_x = self.fig.add_axes([col3_x, 0.20, 0.08, 0.03])
        ax_dual_flip_y = self.fig.add_axes([col3_x + 0.09, 0.20, 0.08, 0.03])
        self.btn_flip_x = Button(ax_dual_flip_x, 'Flip X')
        self.btn_flip_y = Button(ax_dual_flip_y, 'Flip Y')
        self.btn_flip_x.on_clicked(self.on_flip_x_click)
        self.btn_flip_y.on_clicked(self.on_flip_y_click)
        
        # PBC replication
        ax_pbc_n = self.fig.add_axes([col3_x, 0.145, 0.18, 0.025])
        self.slider_pbc_n = Slider(ax_pbc_n, 'PBC rep', 0, 3, valinit=0, valstep=1)
        self.slider_pbc_n.on_changed(self.on_pbc_change)
        
        # === BUTTONS ===
        btn_x = 0.88
        ax_load    = self.fig.add_axes([btn_x, 0.30, 0.10, 0.035])
        ax_reset   = self.fig.add_axes([btn_x, 0.255, 0.10, 0.035])
        ax_save    = self.fig.add_axes([btn_x, 0.21, 0.10, 0.035])
        ax_center  = self.fig.add_axes([btn_x, 0.165, 0.10, 0.035])
        ax_hexgrid = self.fig.add_axes([btn_x, 0.12, 0.10, 0.035])
        
        self.btn_load    = Button(ax_load, 'Load XYZ')
        self.btn_reset   = Button(ax_reset, 'Reset Rot')
        self.btn_save    = Button(ax_save, 'Save XYZ')
        self.btn_center  = Button(ax_center, 'Center Mol')
        self.btn_hexgrid = Button(ax_hexgrid, 'Toggle Grid')
        
        self.btn_load.on_clicked(self.on_load_click)
        self.btn_reset.on_clicked(self.on_reset_click)
        self.btn_save.on_clicked(self.on_save_click)
        self.btn_center.on_clicked(self.on_center_click)
        self.btn_hexgrid.on_clicked(self.on_hexgrid_click)
        
        # Additional export buttons (second row)
        btn_x2 = 0.88
        ax_save_mol2  = self.fig.add_axes([btn_x2, 0.075, 0.10, 0.035])
        ax_save_pose  = self.fig.add_axes([btn_x2, 0.030, 0.10, 0.035])
        ax_save_state = self.fig.add_axes([0.75, 0.030, 0.06, 0.035])
        ax_load_state = self.fig.add_axes([0.82, 0.030, 0.06, 0.035])
        
        self.btn_save_mol2  = Button(ax_save_mol2, 'Save MOL2')
        self.btn_save_pose  = Button(ax_save_pose, 'Save Poses')
        self.btn_save_state = Button(ax_save_state, 'Save St')
        self.btn_load_state = Button(ax_load_state, 'Load St')
        
        self.btn_save_mol2.on_clicked(self.on_save_mol2_click)
        self.btn_save_pose.on_clicked(self.on_save_pose_click)
        self.btn_save_state.on_clicked(self.on_save_state_click)
        self.btn_load_state.on_clicked(self.on_load_state_click)
        
        # Status text
        self.status_text = self.fig.text(0.03, 0.01, 'No molecule loaded', fontsize=10)
        self.info_text   = self.fig.text(0.03, 0.04, '', fontsize=9, family='monospace')
        
    def load_molecule(self, fname):
        """Load molecule from xyz file."""
        try:
            self.apos_orig, self.enames, comment = load_xyz_simple(fname)
            self.apos = self.apos_orig.copy()
            self.current_file = fname
            # Precompute colors and sizes for fast rendering
            self._colors_cache = [ELEMENT_COLORS.get(e, '#808080') for e in self.enames]
            self._sizes_cache  = np.array([ELEMENT_SIZES.get(e, 40) for e in self.enames])
            # Compute bonds for MOL2 export (use distance-based detection)
            if HAS_ATOMICUTILS:
                try:
                    # Use fixed cutoff instead of VdW radii (byRvdW=False)
                    self.bonds, _ = au.findBondsNP(self.apos_orig, Rcut=1.8, byRvdW=False)
                except Exception as e:
                    print(f"Warning: could not compute bonds: {e}")
                    self.bonds = None
            self.status_text.set_text(f'Loaded: {os.path.basename(fname)} ({len(self.apos)} atoms)')
            self.update_plot()
        except Exception as e:
            self.status_text.set_text(f'Error loading: {e}')
    
    def apply_rotation(self):
        """Apply current rotation to molecule."""
        if self.apos_orig is None:
            return
        # Center molecule first
        center = self.apos_orig.mean(axis=0)
        pos = self.apos_orig - center
        # Apply rotations in order X, Y, Z
        Rx = rotation_matrix_x(self.rot_x)
        Ry = rotation_matrix_y(self.rot_y)
        Rz = rotation_matrix_z(self.rot_z)
        R = Rz @ Ry @ Rx
        self.apos = (R @ pos.T).T + center
    
    def get_pivot_position(self, use_up=None):
        """Get pivot position based on triangle settings."""
        if use_up is None:
            use_up = self.use_tri_up
        i = int(self.slider_tri_i.val)
        j = int(self.slider_tri_j.val)
        ba = self.slider_tri_ba.val
        bb = self.slider_tri_bb.val
        
        # Base triangle center
        center = self.lattice.triangle_center(i, j, up=use_up)
        # Add fine offset
        offset = ba * self.lattice.a1 + bb * self.lattice.a2
        return center + offset
    
    def generate_flower(self, pivot, extra_rot=0.0, flip_x=False, flip_y=False):
        """Generate n-fold symmetric copies around given pivot.
        Returns list of (positions, flower_id) tuples.
        """
        if self.apos is None:
            return []
        
        copies = []
        n = int(self.slider_nfold.val)
        radius = self.slider_radius.val
        init_ang = self.slider_initang.val + extra_rot
        
        mol_center = self.apos.mean(axis=0)
        pos_centered = self.apos - mol_center  # center at origin
        
        # Pre-apply flip if needed
        if flip_x:
            pos_centered = (flip_matrix_x() @ pos_centered.T).T
        if flip_y:
            pos_centered = (flip_matrix_y() @ pos_centered.T).T
        
        for i in range(n):
            angle = init_ang + i * (360.0 / n)
            rad = np.radians(angle)
            offset = np.array([radius * np.cos(rad), radius * np.sin(rad), 0])
            Rz = rotation_matrix_z(angle)
            rotated = (Rz @ pos_centered.T).T
            final_pos = rotated + pivot + offset
            copies.append(final_pos)
        
        return copies
    
    def get_all_copies(self):
        """Generate all copies including dual flower and PBC replications."""
        if self.apos is None:
            return [], []
        
        all_copies = []
        all_flower_ids = []  # 0 = primary, 1 = dual
        
        # Primary flower at up-triangle
        pivot1 = self.get_pivot_position(use_up=self.use_tri_up)
        flower1 = self.generate_flower(pivot1)
        all_copies.extend(flower1)
        all_flower_ids.extend([0] * len(flower1))
        
        # Dual flower at opposite triangle
        if self.dual_flower:
            pivot2 = self.get_pivot_position(use_up=not self.use_tri_up)
            flower2 = self.generate_flower(pivot2, 
                                           extra_rot=self.slider_dual_rot.val,
                                           flip_x=self.dual_flip_x,
                                           flip_y=self.dual_flip_y)
            all_copies.extend(flower2)
            all_flower_ids.extend([1] * len(flower2))
        
        # PBC replication
        pbc_n = int(self.slider_pbc_n.val)
        if pbc_n > 0:
            base_copies = all_copies.copy()
            base_ids = all_flower_ids.copy()
            for di in range(-pbc_n, pbc_n + 1):
                for dj in range(-pbc_n, pbc_n + 1):
                    if di == 0 and dj == 0:
                        continue
                    shift = di * self.lattice.a1 + dj * self.lattice.a2
                    for pos, fid in zip(base_copies, base_ids):
                        all_copies.append(pos + shift)
                        all_flower_ids.append(fid + 2)  # mark as PBC copy
        
        return all_copies, all_flower_ids
    
    def draw_hexagonal_grid(self, ax):
        """Draw hexagonal lattice grid."""
        a1 = self.lattice.a1[:2]
        a2 = self.lattice.a2[:2]
        
        # Draw grid lines
        for i in range(-3, 4):
            for j in range(-3, 4):
                origin = i * a1 + j * a2
                corners = [origin, origin + a1, origin + a1 + a2, origin + a2, origin]
                xs = [c[0] for c in corners]
                ys = [c[1] for c in corners]
                ax.plot(xs, ys, 'b-', alpha=0.2, linewidth=0.5)
        
        # Mark triangle centers
        for i in range(-2, 3):
            for j in range(-2, 3):
                c1 = self.lattice.triangle_center(i, j, up=True)[:2]
                c2 = self.lattice.triangle_center(i, j, up=False)[:2]
                ax.plot(c1[0], c1[1], '^', color='red', markersize=6, alpha=0.4)
                ax.plot(c2[0], c2[1], 'v', color='blue', markersize=6, alpha=0.4)
        
        # Highlight current pivot(s)
        pivot1 = self.get_pivot_position(use_up=self.use_tri_up)[:2]
        ax.plot(pivot1[0], pivot1[1], 'o', color='red', markersize=12, mew=2, mfc='none')
        ax.plot(pivot1[0], pivot1[1], '+', color='red', markersize=15, mew=2)
        
        if self.dual_flower:
            pivot2 = self.get_pivot_position(use_up=not self.use_tri_up)[:2]
            ax.plot(pivot2[0], pivot2[1], 'o', color='blue', markersize=12, mew=2, mfc='none')
            ax.plot(pivot2[0], pivot2[1], '+', color='blue', markersize=15, mew=2)
    
    def update_plot(self):
        """Update the visualization - optimized for speed."""
        self.ax_main.clear()
        self.ax_side.clear()
        
        self.ax_main.set_xlabel('X (Å)')
        self.ax_main.set_ylabel('Y (Å)')
        self.ax_main.set_title('Top View (XY)')
        self.ax_main.grid(True, alpha=0.3)
        self.ax_main.set_aspect('equal')
        
        self.ax_side.set_xlabel('X (Å)')
        self.ax_side.set_ylabel('Z (Å)')
        self.ax_side.set_title('Side View (XZ)')
        self.ax_side.grid(True, alpha=0.3)
        self.ax_side.set_aspect('equal')
        
        if self.show_hex_grid:
            self.draw_hexagonal_grid(self.ax_main)
        
        if self.apos is None:
            self.fig.canvas.draw_idle()
            return
        
        # Get all copies (flowers + PBC)
        copies, flower_ids = self.get_all_copies()
        
        if not copies:
            self.fig.canvas.draw_idle()
            return
        
        # Batch rendering for speed: concatenate all positions
        n_atoms = len(self.apos)
        n_copies = len(copies)
        total_atoms = n_atoms * n_copies
        
        # Preallocate arrays
        all_x = np.empty(total_atoms)
        all_y = np.empty(total_atoms)
        all_z = np.empty(total_atoms)
        all_colors = []
        all_sizes = np.empty(total_atoms)
        all_edge_colors = []
        
        # Color scheme: primary=red edge, dual=blue edge, PBC=gray edge
        edge_color_map = {0: 'red', 1: 'blue', 2: 'gray', 3: 'lightgray'}
        
        for idx, (pos, fid) in enumerate(zip(copies, flower_ids)):
            i0, i1 = idx * n_atoms, (idx + 1) * n_atoms
            all_x[i0:i1] = pos[:, 0]
            all_y[i0:i1] = pos[:, 1]
            all_z[i0:i1] = pos[:, 2]
            all_colors.extend(self._colors_cache)
            all_sizes[i0:i1] = self._sizes_cache
            ec = edge_color_map.get(fid, 'gray')
            all_edge_colors.extend([ec] * n_atoms)
        
        # Single scatter call for top view - much faster
        alpha = 0.8 if n_copies < 20 else 0.5  # reduce alpha for many copies
        self.ax_main.scatter(all_x, all_y, c=all_colors, s=all_sizes,
                             alpha=alpha, edgecolors=all_edge_colors, linewidths=0.5)
        
        # Single scatter call for side view
        self.ax_side.scatter(all_x, all_z, c=all_colors, s=all_sizes,
                             alpha=alpha, edgecolors=all_edge_colors, linewidths=0.5)
        
        # Auto-scale
        margin = 5
        xmin, xmax = all_x.min() - margin, all_x.max() + margin
        ymin, ymax = all_y.min() - margin, all_y.max() + margin
        zmin, zmax = all_z.min() - margin, all_z.max() + margin
        
        if self.show_hex_grid:
            a = self.lattice.a
            xmin = min(xmin, -a)
            xmax = max(xmax, 2 * a)
            ymin = min(ymin, -a)
            ymax = max(ymax, 2 * a)
        
        self.ax_main.set_xlim(xmin, xmax)
        self.ax_main.set_ylim(ymin, ymax)
        self.ax_side.set_xlim(xmin, xmax)
        self.ax_side.set_ylim(zmin, zmax)
        
        # Update info text
        pivot1 = self.get_pivot_position()
        info = f"Atoms: {total_atoms} | Pivot: ({pivot1[0]:.1f}, {pivot1[1]:.1f})"
        if self.dual_flower:
            info += f" | Dual: ON (rot={self.slider_dual_rot.val:.0f}°)"
        self.info_text.set_text(info)
        
        self.fig.canvas.draw_idle()
    
    def on_rotation_change(self, val):
        self.rot_x = self.slider_rx.val
        self.rot_y = self.slider_ry.val
        self.rot_z = self.slider_rz.val
        self.apply_rotation()
        self.update_plot()
    
    def on_symmetry_change(self, val):
        self.update_plot()
    
    def on_triangle_change(self, val):
        self.update_plot()
    
    def on_lattice_change(self, val):
        self.lattice.set_a(self.slider_lattice.val)
        self.update_plot()
    
    def on_tri_type_click(self, event):
        self.use_tri_up = not self.use_tri_up
        self.btn_tri_type.label.set_text('▲ Up' if self.use_tri_up else '▼ Down')
        self.update_plot()
    
    def on_dual_toggle(self, event):
        self.dual_flower = not self.dual_flower
        self.btn_dual.label.set_text('Dual: ON' if self.dual_flower else 'Dual: OFF')
        self.update_plot()
    
    def on_dual_change(self, val):
        self.update_plot()
    
    def on_flip_x_click(self, event):
        self.dual_flip_x = not self.dual_flip_x
        self.btn_flip_x.label.set_text('Flip X ✓' if self.dual_flip_x else 'Flip X')
        self.update_plot()
    
    def on_flip_y_click(self, event):
        self.dual_flip_y = not self.dual_flip_y
        self.btn_flip_y.label.set_text('Flip Y ✓' if self.dual_flip_y else 'Flip Y')
        self.update_plot()
    
    def on_pbc_change(self, val):
        self.update_plot()
    
    def on_load_click(self, event):
        """Open file dialog to load xyz."""
        try:
            import tkinter as tk
            from tkinter import filedialog
            root = tk.Tk()
            root.withdraw()
            fname = filedialog.askopenfilename(
                title='Select XYZ file',
                filetypes=[('XYZ files', '*.xyz'), ('All files', '*.*')]
            )
            root.destroy()
            if fname:
                self.load_molecule(fname)
        except Exception as e:
            self.status_text.set_text(f'Error: {e}')
    
    def on_reset_click(self, event):
        """Reset rotation angles."""
        self.slider_rx.set_val(0)
        self.slider_ry.set_val(0)
        self.slider_rz.set_val(0)
    
    def on_center_click(self, event):
        """Center molecule at origin."""
        if self.apos_orig is not None:
            center = self.apos_orig.mean(axis=0)
            self.apos_orig = self.apos_orig - center
            self.apply_rotation()
            self.update_plot()
            self.status_text.set_text('Molecule centered at origin')
    
    def on_hexgrid_click(self, event):
        """Toggle hexagonal grid display."""
        self.show_hex_grid = not self.show_hex_grid
        self.update_plot()
    
    def on_save_click(self, event):
        """Save all copies to xyz file."""
        if self.apos is None:
            self.status_text.set_text('No molecule to save')
            return
        
        copies, _ = self.get_all_copies()
        all_pos = []
        all_names = []
        for pos in copies:
            all_pos.extend(pos.tolist())
            all_names.extend(self.enames)
        
        # Build comment with settings
        comment = f"N-fold={int(self.slider_nfold.val)} a={self.lattice.a:.2f}"
        if self.dual_flower:
            comment += f" dual_rot={self.slider_dual_rot.val:.0f}"
        
        try:
            import tkinter as tk
            from tkinter import filedialog
            root = tk.Tk()
            root.withdraw()
            fname = filedialog.asksaveasfilename(
                title='Save XYZ file',
                defaultextension='.xyz',
                filetypes=[('XYZ files', '*.xyz'), ('All files', '*.*')]
            )
            root.destroy()
            if fname:
                save_xyz_simple(fname, np.array(all_pos), all_names, comment)
                self.status_text.set_text(f'Saved: {fname} ({len(all_pos)} atoms)')
        except Exception as e:
            self.status_text.set_text(f'Error saving: {e}')
    
    def on_save_mol2_click(self, event):
        """Save all copies to MOL2 file with bond topology."""
        if self.apos is None:
            self.status_text.set_text('No molecule to save')
            return
        if not HAS_ATOMICUTILS:
            self.status_text.set_text('MOL2 export requires atomicUtils')
            return
        if self.bonds is None or len(self.bonds) == 0:
            self.status_text.set_text('No bonds detected in molecule')
            return
        
        copies, _ = self.get_all_copies()
        all_pos = []
        all_names = []
        for pos in copies:
            all_pos.extend(pos.tolist())
            all_names.extend(self.enames)
        
        comment = f"N-fold={int(self.slider_nfold.val)} a={self.lattice.a:.2f}"
        
        try:
            import tkinter as tk
            from tkinter import filedialog
            root = tk.Tk()
            root.withdraw()
            fname = filedialog.asksaveasfilename(
                title='Save MOL2 file',
                defaultextension='.mol2',
                filetypes=[('MOL2 files', '*.mol2'), ('All files', '*.*')]
            )
            root.destroy()
            if fname:
                save_mol2_with_bonds(fname, np.array(all_pos), all_names, 
                                     self.bonds, len(self.enames), comment)
                self.status_text.set_text(f'Saved MOL2: {fname} ({len(all_pos)} atoms, {len(self.bonds)*len(copies)} bonds)')
        except Exception as e:
            self.status_text.set_text(f'Error saving MOL2: {e}')
    
    def get_all_poses(self):
        """Get rigid-body poses (position, rotation_matrix) for all molecule copies."""
        if self.apos is None:
            return []
        
        poses = []
        n = int(self.slider_nfold.val)
        radius = self.slider_radius.val
        
        # Base rotation from sliders
        Rx = rotation_matrix_x(self.rot_x)
        Ry = rotation_matrix_y(self.rot_y)
        Rz_base = rotation_matrix_z(self.rot_z)
        R_mol = Rz_base @ Ry @ Rx
        
        def add_flower_poses(pivot, init_ang, flip_x, flip_y):
            R_flip = np.eye(3)
            if flip_x:
                R_flip = flip_matrix_x() @ R_flip
            if flip_y:
                R_flip = flip_matrix_y() @ R_flip
            
            for i in range(n):
                angle = init_ang + i * (360.0 / n)
                rad = np.radians(angle)
                offset = np.array([radius * np.cos(rad), radius * np.sin(rad), 0])
                Rz = rotation_matrix_z(angle)
                R_total = Rz @ R_flip @ R_mol
                pos = pivot + offset
                poses.append((pos, R_total))
        
        # Primary flower
        pivot1 = self.get_pivot_position(use_up=self.use_tri_up)
        init_ang1 = self.slider_initang.val
        add_flower_poses(pivot1, init_ang1, False, False)
        
        # Dual flower
        if self.dual_flower:
            pivot2 = self.get_pivot_position(use_up=not self.use_tri_up)
            init_ang2 = self.slider_initang.val + self.slider_dual_rot.val
            add_flower_poses(pivot2, init_ang2, self.dual_flip_x, self.dual_flip_y)
        
        # PBC replication
        pbc_n = int(self.slider_pbc_n.val)
        if pbc_n > 0:
            base_poses = poses.copy()
            for di in range(-pbc_n, pbc_n + 1):
                for dj in range(-pbc_n, pbc_n + 1):
                    if di == 0 and dj == 0:
                        continue
                    shift = di * self.lattice.a1 + dj * self.lattice.a2
                    for pos, R in base_poses:
                        poses.append((pos + shift, R))
        
        return poses
    
    def on_save_pose_click(self, event):
        """Save rigid-body poses to file."""
        if self.apos is None:
            self.status_text.set_text('No molecule to save')
            return
        
        poses = self.get_all_poses()
        comment = f"N-fold={int(self.slider_nfold.val)} a={self.lattice.a:.2f}"
        
        try:
            import tkinter as tk
            from tkinter import filedialog
            root = tk.Tk()
            root.withdraw()
            fname = filedialog.asksaveasfilename(
                title='Save Poses file',
                defaultextension='.poses',
                filetypes=[('Poses files', '*.poses'), ('Text files', '*.txt'), ('All files', '*.*')]
            )
            root.destroy()
            if fname:
                save_poses(fname, poses, comment)
                self.status_text.set_text(f'Saved poses: {fname} ({len(poses)} molecules)')
        except Exception as e:
            self.status_text.set_text(f'Error saving poses: {e}')
    
    def get_state_dict(self):
        """Get current slider/button state as dictionary."""
        return {
            'rot_x': self.slider_rx.val,
            'rot_y': self.slider_ry.val,
            'rot_z': self.slider_rz.val,
            'n_fold': int(self.slider_nfold.val),
            'radius': self.slider_radius.val,
            'init_angle': self.slider_initang.val,
            'lattice_a': self.slider_lattice.val,
            'cell_i': int(self.slider_tri_i.val),
            'cell_j': int(self.slider_tri_j.val),
            'offset_a': self.slider_tri_ba.val,
            'offset_b': self.slider_tri_bb.val,
            'use_tri_up': self.use_tri_up,
            'dual_flower': self.dual_flower,
            'dual_rot': self.slider_dual_rot.val,
            'dual_flip_x': self.dual_flip_x,
            'dual_flip_y': self.dual_flip_y,
            'pbc_n': int(self.slider_pbc_n.val),
            'molecule_file': getattr(self, 'current_file', None),
        }
    
    def set_state_dict(self, state):
        """Set slider/button state from dictionary."""
        if 'rot_x' in state: self.slider_rx.set_val(state['rot_x'])
        if 'rot_y' in state: self.slider_ry.set_val(state['rot_y'])
        if 'rot_z' in state: self.slider_rz.set_val(state['rot_z'])
        if 'n_fold' in state: self.slider_nfold.set_val(state['n_fold'])
        if 'radius' in state: self.slider_radius.set_val(state['radius'])
        if 'init_angle' in state: self.slider_initang.set_val(state['init_angle'])
        if 'lattice_a' in state: self.slider_lattice.set_val(state['lattice_a'])
        if 'cell_i' in state: self.slider_tri_i.set_val(state['cell_i'])
        if 'cell_j' in state: self.slider_tri_j.set_val(state['cell_j'])
        if 'offset_a' in state: self.slider_tri_ba.set_val(state['offset_a'])
        if 'offset_b' in state: self.slider_tri_bb.set_val(state['offset_b'])
        if 'use_tri_up' in state:
            self.use_tri_up = state['use_tri_up']
            self.btn_tri_type.label.set_text('▲ Up' if self.use_tri_up else '▼ Down')
        if 'dual_flower' in state:
            self.dual_flower = state['dual_flower']
            self.btn_dual.label.set_text('Dual: ON' if self.dual_flower else 'Dual: OFF')
        if 'dual_rot' in state: self.slider_dual_rot.set_val(state['dual_rot'])
        if 'dual_flip_x' in state:
            self.dual_flip_x = state['dual_flip_x']
            self.btn_flip_x.label.set_text('Flip X ✓' if self.dual_flip_x else 'Flip X')
        if 'dual_flip_y' in state:
            self.dual_flip_y = state['dual_flip_y']
            self.btn_flip_y.label.set_text('Flip Y ✓' if self.dual_flip_y else 'Flip Y')
        if 'pbc_n' in state: self.slider_pbc_n.set_val(state['pbc_n'])
        self.update_plot()
    
    def on_save_state_click(self, event):
        """Save current state to JSON file."""
        state = self.get_state_dict()
        try:
            import tkinter as tk
            from tkinter import filedialog
            root = tk.Tk()
            root.withdraw()
            fname = filedialog.asksaveasfilename(
                title='Save State JSON',
                defaultextension='.json',
                filetypes=[('JSON files', '*.json'), ('All files', '*.*')]
            )
            root.destroy()
            if fname:
                with open(fname, 'w') as f:
                    json.dump(state, f, indent=2)
                self.status_text.set_text(f'State saved: {fname}')
        except Exception as e:
            self.status_text.set_text(f'Error saving state: {e}')
    
    def on_load_state_click(self, event):
        """Load state from JSON file."""
        try:
            import tkinter as tk
            from tkinter import filedialog
            root = tk.Tk()
            root.withdraw()
            fname = filedialog.askopenfilename(
                title='Load State JSON',
                filetypes=[('JSON files', '*.json'), ('All files', '*.*')]
            )
            root.destroy()
            if fname:
                with open(fname, 'r') as f:
                    state = json.load(f)
                # Optionally load the molecule if path is in state and exists
                mol_file = state.get('molecule_file')
                if mol_file and os.path.exists(mol_file) and self.apos is None:
                    self.load_molecule(mol_file)
                self.set_state_dict(state)
                self.status_text.set_text(f'State loaded: {fname}')
        except Exception as e:
            self.status_text.set_text(f'Error loading state: {e}')
    
    def run(self):
        """Start the interactive GUI."""
        plt.show()


def main():
    """Main entry point."""
    xyz_file = None
    if len(sys.argv) > 1:
        xyz_file = sys.argv[1]
    
    # Default to one of the helicene molecules if available
    if xyz_file is None:
        default_paths = [
            'pyBall/XPDB_AVBD/pdb/DiTetraceno_helicene_1a.xyz',
            'XPDB_AVBD/pdb/DiTetraceno_helicene_1a.xyz',
        ]
        for p in default_paths:
            if os.path.exists(p):
                xyz_file = p
                break
    
    app = MolecularPlacer(xyz_file)
    app.run()


if __name__ == '__main__':
    main()
