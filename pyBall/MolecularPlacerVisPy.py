#!/usr/bin/env python3
"""
VisPy+PyQt5 GPU-accelerated version of MolecularPlacer.
Interactive molecular placement with rotation, N-fold symmetry, hex lattice, dual-flower, PBC.

Usage:
    python MolecularPlacerVisPy.py [xyz_file]
"""

import numpy as np
import os, sys, json

_THIS_DIR = os.path.abspath(os.path.dirname(__file__))
_ROOT_DIR = os.path.abspath(os.path.join(_THIS_DIR, os.pardir))
for _p in (_ROOT_DIR, _THIS_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from pyBall.MolGUI_common import (
    ELEMENT_COLORS, ELEMENT_SIZES,
    colors_for_enames, sizes_for_enames,
    rotmat_x, rotmat_y, rotmat_z, make_rotmat,
    flip_matrix_x, flip_matrix_y,
    load_molecule_xyz, save_xyz, find_bonds,
    HexLattice,
)

from PyQt5 import QtWidgets, QtCore
import vispy
vispy.use('pyqt5')
from vispy import scene
from vispy.scene import visuals


class MolecularPlacerVisPy(QtWidgets.QMainWindow):

    def __init__(self, xyz_file=None):
        super().__init__()
        self.setWindowTitle('Molecular Placer (VisPy) - Symmetry Tool')
        self.resize(1400, 900)

        # ---- Data ----
        self.apos_orig = None
        self.enames    = None
        self.apos      = None
        self.bonds     = None
        self.current_file = None
        self._colors_cache = None
        self._sizes_cache  = None

        # ---- Params ----
        self.rot_x = 0.0; self.rot_y = 0.0; self.rot_z = 0.0
        self.n_fold = 3; self.radial_offset = 10.0; self.initial_angle = 0.0
        self.lattice = HexLattice(a=32.7)
        self.tri_i = 0; self.tri_j = 0; self.tri_ba = 0.0; self.tri_bb = 0.0
        self.use_tri_up = True
        self.dual_flower = False; self.dual_rot = 180.0
        self.dual_flip_x = False; self.dual_flip_y = False
        self.pbc_n = 0

        self._build_ui()

        if xyz_file and os.path.exists(xyz_file):
            self.load_molecule(xyz_file)

    # ================================================================
    #  UI
    # ================================================================

    def _build_ui(self):
        central = QtWidgets.QWidget(); self.setCentralWidget(central)
        hbox = QtWidgets.QHBoxLayout(central)

        # ---- VisPy canvas ----
        self.canvas = scene.SceneCanvas(keys='interactive', bgcolor='white', show=False)
        self.view = self.canvas.central_widget.add_view()
        self.view.camera = scene.TurntableCamera(fov=0, distance=80, elevation=90, azimuth=0)  # top-down ortho

        self.mol_markers = visuals.Markers(parent=self.view.scene)
        self.bond_lines  = visuals.Line(parent=self.view.scene, color='gray', width=1.5, antialias=True, method='gl')
        self.grid_lines  = visuals.Line(parent=self.view.scene, color=(0,0,1,0.15), width=0.5, antialias=True, method='gl')
        self.pivot_markers = visuals.Markers(parent=self.view.scene)

        hbox.addWidget(self.canvas.native, stretch=3)

        # ---- Controls ----
        scroll = QtWidgets.QScrollArea()
        scroll.setWidgetResizable(True); scroll.setMaximumWidth(340)
        ctrl = QtWidgets.QWidget()
        vbox = QtWidgets.QVBoxLayout(ctrl)

        # Rotation
        grp_rot = QtWidgets.QGroupBox("Molecule Rotation")
        gl = QtWidgets.QGridLayout(grp_rot)
        self.sp_rx = self._dspin(gl, 0, "Rx (°)", 0, 5, -360, 360, 1)
        self.sp_ry = self._dspin(gl, 1, "Ry (°)", 0, 5, -360, 360, 1)
        self.sp_rz = self._dspin(gl, 2, "Rz (°)", 0, 5, -360, 360, 1)
        vbox.addWidget(grp_rot)

        # Symmetry
        grp_sym = QtWidgets.QGroupBox("Flower Symmetry")
        gl2 = QtWidgets.QGridLayout(grp_sym)
        self.sp_nfold  = self._ispin(gl2, 0, "N-fold", 3, 1, 1, 12)
        self.sp_radius = self._dspin(gl2, 1, "Radius", 10, 1, 0, 100, 1)
        self.sp_initang= self._dspin(gl2, 2, "Init φ (°)", 0, 5, -360, 360, 1)
        vbox.addWidget(grp_sym)

        # Lattice
        grp_lat = QtWidgets.QGroupBox("Lattice / Triangle")
        gl3 = QtWidgets.QGridLayout(grp_lat)
        self.sp_lat_a = self._dspin(gl3, 0, "Lattice a", 32.7, 1, 5, 100, 1)
        self.sp_tri_i = self._ispin(gl3, 1, "Cell i", 0, 1, -5, 5)
        self.sp_tri_j = self._ispin(gl3, 2, "Cell j", 0, 1, -5, 5)
        self.sp_tri_ba= self._dspin(gl3, 3, "Offset a", 0, 0.05, -1, 1, 3)
        self.sp_tri_bb= self._dspin(gl3, 4, "Offset b", 0, 0.05, -1, 1, 3)
        self.btn_tri_type = QtWidgets.QPushButton("▲ Up")
        self.btn_tri_type.clicked.connect(self._on_tri_type)
        gl3.addWidget(self.btn_tri_type, 5, 0, 1, 2)
        vbox.addWidget(grp_lat)

        # Dual flower
        grp_dual = QtWidgets.QGroupBox("Dual Flower")
        gl4 = QtWidgets.QGridLayout(grp_dual)
        self.btn_dual = QtWidgets.QPushButton("Dual: OFF")
        self.btn_dual.clicked.connect(self._on_dual_toggle)
        gl4.addWidget(self.btn_dual, 0, 0, 1, 2)
        self.sp_dual_rot = self._dspin(gl4, 1, "Dual φ (°)", 180, 5, -360, 360, 1)
        self.btn_flip_x = QtWidgets.QPushButton("Flip X")
        self.btn_flip_x.clicked.connect(self._on_flip_x)
        gl4.addWidget(self.btn_flip_x, 2, 0)
        self.btn_flip_y = QtWidgets.QPushButton("Flip Y")
        self.btn_flip_y.clicked.connect(self._on_flip_y)
        gl4.addWidget(self.btn_flip_y, 2, 1)
        self.sp_pbc = self._ispin(gl4, 3, "PBC rep", 0, 1, 0, 5)
        vbox.addWidget(grp_dual)

        # Buttons
        grp_btn = QtWidgets.QGroupBox("Actions")
        gl5 = QtWidgets.QVBoxLayout(grp_btn)
        btn_load = QtWidgets.QPushButton("Load XYZ"); btn_load.clicked.connect(self._on_load); gl5.addWidget(btn_load)
        btn_center = QtWidgets.QPushButton("Center Mol"); btn_center.clicked.connect(self._on_center); gl5.addWidget(btn_center)
        btn_reset = QtWidgets.QPushButton("Reset Rot"); btn_reset.clicked.connect(self._on_reset_rot); gl5.addWidget(btn_reset)
        btn_save = QtWidgets.QPushButton("Save XYZ"); btn_save.clicked.connect(self._on_save_xyz); gl5.addWidget(btn_save)
        bhbox = QtWidgets.QHBoxLayout()
        btn_save_st = QtWidgets.QPushButton("Save State"); btn_save_st.clicked.connect(self._on_save_state); bhbox.addWidget(btn_save_st)
        btn_load_st = QtWidgets.QPushButton("Load State"); btn_load_st.clicked.connect(self._on_load_state); bhbox.addWidget(btn_load_st)
        gl5.addLayout(bhbox)
        btn_view = QtWidgets.QPushButton("Reset View"); btn_view.clicked.connect(self._on_reset_view); gl5.addWidget(btn_view)
        vbox.addWidget(grp_btn)

        self.status = QtWidgets.QLabel("No molecule loaded")
        vbox.addWidget(self.status)
        vbox.addStretch()

        scroll.setWidget(ctrl)
        hbox.addWidget(scroll)

        # Connect all spinboxes for instant update
        for sp in (self.sp_rx, self.sp_ry, self.sp_rz,
                   self.sp_nfold, self.sp_radius, self.sp_initang,
                   self.sp_lat_a, self.sp_tri_i, self.sp_tri_j, self.sp_tri_ba, self.sp_tri_bb,
                   self.sp_dual_rot, self.sp_pbc):
            sp.valueChanged.connect(self._on_param_changed)

    def _dspin(self, grid, row, label, val, step, lo, hi, dec):
        grid.addWidget(QtWidgets.QLabel(label), row, 0)
        sp = QtWidgets.QDoubleSpinBox(); sp.setRange(lo, hi); sp.setSingleStep(step); sp.setDecimals(dec); sp.setValue(val)
        grid.addWidget(sp, row, 1); return sp

    def _ispin(self, grid, row, label, val, step, lo, hi):
        grid.addWidget(QtWidgets.QLabel(label), row, 0)
        sp = QtWidgets.QSpinBox(); sp.setRange(lo, hi); sp.setSingleStep(step); sp.setValue(val)
        grid.addWidget(sp, row, 1); return sp

    # ================================================================
    #  Molecule loading
    # ================================================================

    def load_molecule(self, fname):
        self.apos_orig, self.enames = load_molecule_xyz(fname)
        self.apos = self.apos_orig.copy()
        self.current_file = fname
        self._colors_cache = colors_for_enames(self.enames, alpha=0.9)
        self._sizes_cache  = sizes_for_enames(self.enames, scale=0.25)
        self.bonds = find_bonds(self.apos_orig)
        self.status.setText(f"Loaded: {os.path.basename(fname)} ({len(self.apos)} atoms)")
        self._apply_rotation()
        self._update_visuals()

    # ================================================================
    #  Rotation / symmetry
    # ================================================================

    def _apply_rotation(self):
        if self.apos_orig is None: return
        center = self.apos_orig.mean(axis=0)
        pos = self.apos_orig - center
        R = make_rotmat(self.sp_rx.value(), self.sp_ry.value(), self.sp_rz.value())
        self.apos = (R @ pos.T).T + center

    def _get_pivot(self, use_up=None):
        if use_up is None: use_up = self.use_tri_up
        i = self.sp_tri_i.value(); j = self.sp_tri_j.value()
        ba = self.sp_tri_ba.value(); bb = self.sp_tri_bb.value()
        center = self.lattice.triangle_center(i, j, up=use_up)
        offset = ba * self.lattice.a1 + bb * self.lattice.a2
        return center + offset

    def _generate_flower(self, pivot, extra_rot=0.0, flip_x=False, flip_y=False):
        if self.apos is None: return []
        n = self.sp_nfold.value()
        radius = self.sp_radius.value()
        init_ang = self.sp_initang.value() + extra_rot
        mol_center = self.apos.mean(axis=0)
        pos_centered = self.apos - mol_center
        if flip_x: pos_centered = (flip_matrix_x() @ pos_centered.T).T
        if flip_y: pos_centered = (flip_matrix_y() @ pos_centered.T).T
        copies = []
        for i in range(n):
            angle = init_ang + i * (360.0 / n)
            rad = np.radians(angle)
            offset = np.array([radius * np.cos(rad), radius * np.sin(rad), 0])
            Rz = rotmat_z(angle)
            rotated = (Rz @ pos_centered.T).T
            copies.append(rotated + pivot + offset)
        return copies

    def _get_all_copies(self):
        if self.apos is None: return [], []
        all_copies = []; all_ids = []
        pivot1 = self._get_pivot(self.use_tri_up)
        f1 = self._generate_flower(pivot1)
        all_copies.extend(f1); all_ids.extend([0]*len(f1))
        if self.dual_flower:
            pivot2 = self._get_pivot(not self.use_tri_up)
            f2 = self._generate_flower(pivot2, extra_rot=self.sp_dual_rot.value(),
                                       flip_x=self.dual_flip_x, flip_y=self.dual_flip_y)
            all_copies.extend(f2); all_ids.extend([1]*len(f2))
        pbc_n = self.sp_pbc.value()
        if pbc_n > 0:
            base_c = all_copies[:]; base_i = all_ids[:]
            for di in range(-pbc_n, pbc_n+1):
                for dj in range(-pbc_n, pbc_n+1):
                    if di==0 and dj==0: continue
                    shift = di*self.lattice.a1 + dj*self.lattice.a2
                    for p, fid in zip(base_c, base_i):
                        all_copies.append(p + shift)
                        all_ids.append(fid + 2)
        return all_copies, all_ids

    # ================================================================
    #  Visual updates
    # ================================================================

    def _on_param_changed(self, _=None):
        self.lattice.set_a(self.sp_lat_a.value())
        self._apply_rotation()
        self._update_visuals()

    def _update_visuals(self):
        copies, flower_ids = self._get_all_copies()
        if not copies:
            self.mol_markers.set_data(np.zeros((0,3), dtype=np.float32))
            self.bond_lines.set_data(np.zeros((0,3), dtype=np.float32))
            self._draw_hex_grid()
            self.canvas.update(); return

        n_atoms = len(self.apos)
        n_copies = len(copies)
        total = n_atoms * n_copies

        all_pos = np.empty((total, 3), dtype=np.float32)
        all_colors = np.empty((total, 4), dtype=np.float32)
        all_sizes = np.empty(total, dtype=np.float32)

        edge_tints = {0: np.array([1,0,0,1], dtype=np.float32),    # primary: red tint
                      1: np.array([0,0,1,1], dtype=np.float32),    # dual: blue
                      2: np.array([0.5,0.5,0.5,0.4], dtype=np.float32),
                      3: np.array([0.5,0.5,0.5,0.3], dtype=np.float32)}

        for idx, (pos, fid) in enumerate(zip(copies, flower_ids)):
            i0 = idx * n_atoms; i1 = i0 + n_atoms
            all_pos[i0:i1] = pos.astype(np.float32)
            c = self._colors_cache.copy()
            if fid >= 2:  # PBC: fade
                c[:, 3] *= 0.4
            all_colors[i0:i1] = c
            all_sizes[i0:i1] = self._sizes_cache

        self.mol_markers.set_data(all_pos, face_color=all_colors, size=all_sizes, edge_width=0.5, edge_color='black')

        # Bonds (vectorized interleave: start,end pairs)
        if self.bonds is not None and len(self.bonds) > 0:
            nb = len(self.bonds)
            segs = np.empty((n_copies * nb * 2, 3), dtype=np.float32)
            connect = np.zeros(n_copies * nb * 2, dtype=bool)
            for idx, pos in enumerate(copies):
                fp = pos.astype(np.float32)
                off = idx * nb * 2
                segs[off  :off+nb*2:2] = fp[self.bonds[:,0]]
                segs[off+1:off+nb*2:2] = fp[self.bonds[:,1]]
                connect[off:off+nb*2:2] = True
            self.bond_lines.set_data(segs, connect=connect, color=(0.3,0.3,0.3,0.7), width=1.2)
        else:
            self.bond_lines.set_data(np.zeros((0,3), dtype=np.float32))

        self._draw_hex_grid()
        self.status.setText(f"Atoms: {total} | Copies: {n_copies}")
        self.canvas.update()

    def _draw_hex_grid(self):
        """Draw hex lattice grid as lines."""
        a1 = self.lattice.a1; a2 = self.lattice.a2
        pts = []; conn = []
        idx = 0
        for i in range(-3, 4):
            for j in range(-3, 4):
                o = i*a1 + j*a2
                corners = [o, o+a1, o+a1+a2, o+a2, o]
                for k in range(4):
                    pts.append(corners[k]); pts.append(corners[k+1])
                    conn.append(True); conn.append(False)
        if pts:
            self.grid_lines.set_data(np.array(pts, dtype=np.float32), connect=np.array(conn, dtype=bool),
                                     color=(0,0,1,0.12), width=0.8)
        # pivot marker
        pivot1 = self._get_pivot(self.use_tri_up)
        pivots = [pivot1]
        colors = [np.array([1,0,0,1], dtype=np.float32)]
        if self.dual_flower:
            pivots.append(self._get_pivot(not self.use_tri_up))
            colors.append(np.array([0,0,1,1], dtype=np.float32))
        pp = np.array(pivots, dtype=np.float32)
        cc = np.array(colors, dtype=np.float32)
        self.pivot_markers.set_data(pp, face_color=cc, size=15, edge_width=2, edge_color='black', symbol='cross')

    # ================================================================
    #  Button handlers
    # ================================================================

    def _on_load(self):
        fname, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Load XYZ", "", "XYZ (*.xyz);;All (*)")
        if fname: self.load_molecule(fname)

    def _on_center(self):
        if self.apos_orig is not None:
            self.apos_orig -= self.apos_orig.mean(axis=0)
            self._apply_rotation(); self._update_visuals()

    def _on_reset_rot(self):
        self.sp_rx.setValue(0); self.sp_ry.setValue(0); self.sp_rz.setValue(0)

    def _on_tri_type(self):
        self.use_tri_up = not self.use_tri_up
        self.btn_tri_type.setText("▲ Up" if self.use_tri_up else "▼ Down")
        self._update_visuals()

    def _on_dual_toggle(self):
        self.dual_flower = not self.dual_flower
        self.btn_dual.setText("Dual: ON" if self.dual_flower else "Dual: OFF")
        self._update_visuals()

    def _on_flip_x(self):
        self.dual_flip_x = not self.dual_flip_x
        self.btn_flip_x.setText("Flip X ✓" if self.dual_flip_x else "Flip X")
        self._update_visuals()

    def _on_flip_y(self):
        self.dual_flip_y = not self.dual_flip_y
        self.btn_flip_y.setText("Flip Y ✓" if self.dual_flip_y else "Flip Y")
        self._update_visuals()

    def _on_reset_view(self):
        self.view.camera.set_state({'elevation': 90, 'azimuth': 0, 'distance': 80})
        self.canvas.update()

    def _on_save_xyz(self):
        if self.apos is None: return
        copies, _ = self._get_all_copies()
        all_pos = []; all_names = []
        for pos in copies:
            all_pos.append(pos)
            all_names.extend(self.enames)
        fname, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save XYZ", "placed.xyz", "XYZ (*.xyz)")
        if fname:
            save_xyz(fname, np.vstack(all_pos), all_names,
                     comment=f"N-fold={self.sp_nfold.value()} a={self.lattice.a:.2f}")
            self.status.setText(f"Saved {fname} ({len(all_names)} atoms)")

    def _on_save_state(self):
        state = {
            'rot_x': self.sp_rx.value(), 'rot_y': self.sp_ry.value(), 'rot_z': self.sp_rz.value(),
            'n_fold': self.sp_nfold.value(), 'radius': self.sp_radius.value(),
            'init_angle': self.sp_initang.value(), 'lattice_a': self.sp_lat_a.value(),
            'cell_i': self.sp_tri_i.value(), 'cell_j': self.sp_tri_j.value(),
            'offset_a': self.sp_tri_ba.value(), 'offset_b': self.sp_tri_bb.value(),
            'use_tri_up': self.use_tri_up, 'dual_flower': self.dual_flower,
            'dual_rot': self.sp_dual_rot.value(), 'dual_flip_x': self.dual_flip_x,
            'dual_flip_y': self.dual_flip_y, 'pbc_n': self.sp_pbc.value(),
            'molecule_file': self.current_file,
        }
        fname, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save State", "state.json", "JSON (*.json)")
        if fname:
            with open(fname, 'w') as f: json.dump(state, f, indent=2)
            self.status.setText(f"State saved: {fname}")

    def _on_load_state(self):
        fname, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Load State", "", "JSON (*.json);;All (*)")
        if not fname: return
        with open(fname, 'r') as f: state = json.load(f)
        mol_file = state.get('molecule_file')
        if mol_file and os.path.exists(mol_file) and self.apos is None:
            self.load_molecule(mol_file)
        # set spinboxes (will trigger _on_param_changed via signals)
        for key, sp in [('rot_x', self.sp_rx), ('rot_y', self.sp_ry), ('rot_z', self.sp_rz),
                        ('n_fold', self.sp_nfold), ('radius', self.sp_radius),
                        ('init_angle', self.sp_initang), ('lattice_a', self.sp_lat_a),
                        ('cell_i', self.sp_tri_i), ('cell_j', self.sp_tri_j),
                        ('offset_a', self.sp_tri_ba), ('offset_b', self.sp_tri_bb),
                        ('dual_rot', self.sp_dual_rot), ('pbc_n', self.sp_pbc)]:
            if key in state: sp.setValue(state[key])
        if 'use_tri_up' in state:
            self.use_tri_up = state['use_tri_up']
            self.btn_tri_type.setText("▲ Up" if self.use_tri_up else "▼ Down")
        if 'dual_flower' in state:
            self.dual_flower = state['dual_flower']
            self.btn_dual.setText("Dual: ON" if self.dual_flower else "Dual: OFF")
        if 'dual_flip_x' in state:
            self.dual_flip_x = state['dual_flip_x']
            self.btn_flip_x.setText("Flip X ✓" if self.dual_flip_x else "Flip X")
        if 'dual_flip_y' in state:
            self.dual_flip_y = state['dual_flip_y']
            self.btn_flip_y.setText("Flip Y ✓" if self.dual_flip_y else "Flip Y")
        self._update_visuals()
        self.status.setText(f"State loaded: {fname}")


# ======================== Main ========================

def main():
    xyz_file = None
    if len(sys.argv) > 1:
        xyz_file = sys.argv[1]

    # Default molecule
    if xyz_file is None:
        defaults = [
            os.path.join(_ROOT_DIR, 'cpp', 'common_resources', 'xyz', 'PTCDA.xyz'),
            'pyBall/XPDB_AVBD/pdb/DiTetraceno_helicene_1a.xyz',
            'XPDB_AVBD/pdb/DiTetraceno_helicene_1a.xyz',
        ]
        for p in defaults:
            if os.path.exists(p):
                xyz_file = p; break

    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
    win = MolecularPlacerVisPy(xyz_file)
    win.show()
    app.exec_()


if __name__ == '__main__':
    main()
