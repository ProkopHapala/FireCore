#!/usr/bin/env python3
"""
VisPy+PyQt5 GPU-accelerated version of SequencePlacer.
Places sequences of molecules on NaCl step-edge substrates with instant visual feedback.

Usage:
    python SequencePlacerVisPy.py                    # interactive GUI
    python SequencePlacerVisPy.py --headless ...     # headless (uses original SequencePlacer)
"""

import numpy as np
import os, sys

_THIS_DIR = os.path.abspath(os.path.dirname(__file__))
_ROOT_DIR = os.path.abspath(os.path.join(_THIS_DIR, os.pardir))
for _p in (_ROOT_DIR, _THIS_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from pyBall import SubstrateBuilder as SB
from pyBall.MolGUI_common import (
    ELEMENT_COLORS, ELEMENT_SIZES,
    colors_for_enames, sizes_for_enames, hex_to_rgba,
    rotmat_x, rotmat_y, rotmat_z, make_rotmat,
    load_molecule_xyz, save_xyz, find_bonds, replicate_bonds,
    place_sequence, load_default_molecules, DEFAULT_MOLS,
)

from PyQt5 import QtWidgets, QtCore
import vispy
vispy.use('pyqt5')
from vispy import scene
from vispy.scene import visuals


class SequencePlacerVisPy(QtWidgets.QMainWindow):

    def __init__(self):
        super().__init__()
        self.setWindowTitle('Sequence Placer (VisPy) - Molecules on NaCl Step Edge')
        self.resize(1400, 900)

        # ---- Data ----
        self.substrate = None
        self.mol_dict      = {}
        self.mol_files     = {}
        self.mol_rotations = {}
        self.mol_bonds     = {}
        self.placed_apos   = None
        self.placed_enames = None
        self.placed_bonds  = None

        # ---- Params ----
        self.sub_a = 2.82065; self.sub_nx = 15; self.sub_ny = 8; self.sub_nz = 3; self.sub_nsteps = 1
        self.row_angle = 0.0; self.row_spacing = 12.0; self.mol_height = 3.0
        self.origin_x = 0.0; self.origin_y = 0.0
        self.sequence = "ABAB"

        self._build_ui()
        self._load_defaults()

    # ================================================================
    #  UI setup
    # ================================================================

    def _build_ui(self):
        central = QtWidgets.QWidget(); self.setCentralWidget(central)
        hbox = QtWidgets.QHBoxLayout(central)

        # ---- Left: VisPy 3D canvas ----
        self.canvas = scene.SceneCanvas(keys='interactive', bgcolor='white', show=False)
        self.view = self.canvas.central_widget.add_view()
        self.view.camera = scene.TurntableCamera(fov=45, distance=80, elevation=30, azimuth=-60)

        # Visuals (created once, data updated)
        self.sub_markers = visuals.Markers(parent=self.view.scene)
        self.mol_markers = visuals.Markers(parent=self.view.scene)
        self.bond_lines  = visuals.Line(parent=self.view.scene, color='gray', width=1.5, antialias=True, method='gl')

        hbox.addWidget(self.canvas.native, stretch=3)

        # ---- Right: controls ----
        ctrl = QtWidgets.QWidget(); ctrl.setMaximumWidth(340)
        vbox = QtWidgets.QVBoxLayout(ctrl)

        # Substrate group
        grp_sub = QtWidgets.QGroupBox("Substrate")
        gl = QtWidgets.QGridLayout(grp_sub)
        self.sp_a  = self._add_double_spin(gl, 0, "a (Å)", self.sub_a, 0.1, 1.0, 10.0, 3)
        self.sp_nx = self._add_int_spin(gl, 1, "nx", self.sub_nx, 1, 1, 50)
        self.sp_ny = self._add_int_spin(gl, 2, "ny", self.sub_ny, 1, 1, 50)
        self.sp_nz = self._add_int_spin(gl, 3, "nz", self.sub_nz, 1, 1, 20)
        self.sp_nsteps = self._add_int_spin(gl, 4, "steps", self.sub_nsteps, 1, 1, 10)
        btn_build = QtWidgets.QPushButton("Build Substrate"); gl.addWidget(btn_build, 5, 0, 1, 2)
        btn_build.clicked.connect(self._on_build_sub)
        vbox.addWidget(grp_sub)

        # Placement group
        grp_place = QtWidgets.QGroupBox("Placement")
        gl2 = QtWidgets.QGridLayout(grp_place)
        self.sp_angle   = self._add_double_spin(gl2, 0, "Row angle (°)", self.row_angle, 5.0, -360, 360, 1)
        self.sp_spacing = self._add_double_spin(gl2, 1, "Spacing (Å)", self.row_spacing, 0.5, 0.1, 100, 1)
        self.sp_height  = self._add_double_spin(gl2, 2, "Height (Å)", self.mol_height, 0.5, -50, 50, 1)
        self.sp_ox      = self._add_double_spin(gl2, 3, "Origin X", self.origin_x, 1.0, -500, 500, 1)
        self.sp_oy      = self._add_double_spin(gl2, 4, "Origin Y", self.origin_y, 1.0, -500, 500, 1)
        vbox.addWidget(grp_place)

        # Rotation group
        grp_rot = QtWidgets.QGroupBox("Molecule Rotation")
        gl3 = QtWidgets.QGridLayout(grp_rot)
        self.sp_rx = self._add_double_spin(gl3, 0, "Rx (°)", 0, 5, -360, 360, 1)
        self.sp_ry = self._add_double_spin(gl3, 1, "Ry (°)", 0, 5, -360, 360, 1)
        self.sp_rz = self._add_double_spin(gl3, 2, "Rz (°)", 0, 5, -360, 360, 1)
        vbox.addWidget(grp_rot)

        # Sequence / Molecules
        grp_seq = QtWidgets.QGroupBox("Sequence")
        gl4 = QtWidgets.QFormLayout(grp_seq)
        self.le_seq = QtWidgets.QLineEdit(self.sequence)
        gl4.addRow("Sequence", self.le_seq)
        self.le_mols = QtWidgets.QLineEdit("")
        gl4.addRow("Mols", self.le_mols)
        vbox.addWidget(grp_seq)

        # Action buttons
        btn_place = QtWidgets.QPushButton("Place Sequence")
        btn_place.clicked.connect(self._on_place)
        vbox.addWidget(btn_place)

        bhbox = QtWidgets.QHBoxLayout()
        btn_save = QtWidgets.QPushButton("Save XYZ"); btn_save.clicked.connect(self._on_save_all)
        btn_savem = QtWidgets.QPushButton("Save Mols"); btn_savem.clicked.connect(self._on_save_mols)
        bhbox.addWidget(btn_save); bhbox.addWidget(btn_savem)
        vbox.addLayout(bhbox)

        btn_reset = QtWidgets.QPushButton("Reset View"); btn_reset.clicked.connect(self._on_reset_view)
        vbox.addWidget(btn_reset)

        # Connect spinboxes for instant update on placement params
        for sp in (self.sp_angle, self.sp_spacing, self.sp_height, self.sp_ox, self.sp_oy,
                   self.sp_rx, self.sp_ry, self.sp_rz):
            sp.valueChanged.connect(self._on_param_changed)
        self.le_seq.textChanged.connect(self._on_param_changed)

        self.status = QtWidgets.QLabel("Ready")
        vbox.addWidget(self.status)
        vbox.addStretch()

        hbox.addWidget(ctrl)

    # ---- Spin helpers ----

    def _add_double_spin(self, grid, row, label, val, step, lo, hi, decimals):
        grid.addWidget(QtWidgets.QLabel(label), row, 0)
        sp = QtWidgets.QDoubleSpinBox()
        sp.setRange(lo, hi); sp.setSingleStep(step); sp.setDecimals(decimals); sp.setValue(val)
        grid.addWidget(sp, row, 1)
        return sp

    def _add_int_spin(self, grid, row, label, val, step, lo, hi):
        grid.addWidget(QtWidgets.QLabel(label), row, 0)
        sp = QtWidgets.QSpinBox()
        sp.setRange(lo, hi); sp.setSingleStep(step); sp.setValue(val)
        grid.addWidget(sp, row, 1)
        return sp

    # ================================================================
    #  Data loading
    # ================================================================

    def _load_defaults(self):
        self.mol_dict, self.mol_bonds, self.mol_files = load_default_molecules()
        self.mol_rotations = {k: (0.,0.,0.) for k in self.mol_dict}
        parts = [f"{k}={v}" for k,v in self.mol_files.items()]
        self.le_mols.setText(",".join(parts))
        if self.mol_dict:
            self.status.setText(f"Loaded defaults: {list(self.mol_dict.keys())}")

    def _parse_mol_assign(self):
        """Parse mol assignment string from le_mols. Returns True on success."""
        text = self.le_mols.text().strip()
        if not text:
            return bool(self.mol_dict)
        new_dict = {}; new_files = {}; new_bonds = {}
        for part in text.split(','):
            part = part.strip()
            if '=' not in part: continue
            letter, path = part.split('=', 1)
            letter = letter.strip(); path = path.strip()
            if len(letter) != 1: continue
            if not os.path.exists(path):
                self.status.setText(f"ERROR: not found: {path}"); return False
            if letter in self.mol_dict and self.mol_files.get(letter) == path:
                new_dict[letter]  = self.mol_dict[letter]
                new_files[letter] = path
                new_bonds[letter] = self.mol_bonds.get(letter)
            else:
                apos, enames = load_molecule_xyz(path)
                apos -= apos.mean(axis=0)
                new_dict[letter]  = (apos, enames)
                new_files[letter] = path
                new_bonds[letter] = find_bonds(apos)
        self.mol_dict  = new_dict
        self.mol_files = new_files
        self.mol_bonds = new_bonds
        for k in self.mol_dict:
            if k not in self.mol_rotations:
                self.mol_rotations[k] = (0.,0.,0.)
        return True

    # ================================================================
    #  Build / Place
    # ================================================================

    def _on_build_sub(self):
        self.sub_a = self.sp_a.value()
        self.sub_nx = self.sp_nx.value(); self.sub_ny = self.sp_ny.value()
        self.sub_nz = self.sp_nz.value(); self.sub_nsteps = self.sp_nsteps.value()
        self.substrate = SB.gen_nacl_step(a=self.sub_a, nx=self.sub_nx, ny=self.sub_ny,
                                           nz=self.sub_nz, nsteps=self.sub_nsteps)
        self._update_sub_visual()
        self.status.setText(f"Substrate: {len(self.substrate.apos)} atoms")
        self._do_place()

    def _on_place(self):
        self._parse_mol_assign()
        self._do_place()

    def _on_param_changed(self, _=None):
        """Called on every spinbox/text change for instant feedback."""
        self._read_params()
        self._do_place()

    def _read_params(self):
        self.row_angle   = self.sp_angle.value()
        self.row_spacing = self.sp_spacing.value()
        self.mol_height  = self.sp_height.value()
        self.origin_x    = self.sp_ox.value()
        self.origin_y    = self.sp_oy.value()
        rx = self.sp_rx.value(); ry = self.sp_ry.value(); rz = self.sp_rz.value()
        for k in self.mol_dict:
            self.mol_rotations[k] = (rx, ry, rz)
        self.sequence = self.le_seq.text().strip()

    def _do_place(self):
        """Recompute placement and update visuals."""
        if not self.mol_dict or not self.sequence:
            self.mol_markers.set_data(np.zeros((0,3), dtype=np.float32))
            self.bond_lines.set_data(np.zeros((0,3), dtype=np.float32))
            return
        for ch in self.sequence:
            if ch not in self.mol_dict:
                self.status.setText(f"Letter '{ch}' not loaded"); return
        ang_rad = np.radians(self.row_angle)
        row_dir = np.array([np.cos(ang_rad), np.sin(ang_rad)])
        z_top = self.substrate.apos[:,2].max() if self.substrate is not None else 0.0
        origin = np.array([self.origin_x, self.origin_y, 0.0])
        self.placed_apos, self.placed_enames = place_sequence(
            self.mol_dict, self.sequence, row_dir, self.row_spacing,
            self.mol_rotations, origin=origin, height=z_top + self.mol_height)
        self.placed_bonds = replicate_bonds(self.mol_bonds, self.sequence, self.mol_dict)
        self._update_mol_visual()
        self.status.setText(f"Placed {len(self.sequence)} mol(s), {len(self.placed_apos)} atoms")

    # ================================================================
    #  VisPy visual updates
    # ================================================================

    def _update_sub_visual(self):
        if self.substrate is None:
            self.sub_markers.set_data(np.zeros((0,3), dtype=np.float32)); return
        sp = self.substrate.apos.astype(np.float32)
        colors = colors_for_enames(self.substrate.enames, alpha=0.45)
        sizes  = sizes_for_enames(self.substrate.enames, scale=0.15)
        self.sub_markers.set_data(sp, face_color=colors, size=sizes, edge_width=0, edge_color=None)

    def _update_mol_visual(self):
        if self.placed_apos is None or len(self.placed_apos) == 0:
            self.mol_markers.set_data(np.zeros((0,3), dtype=np.float32))
            self.bond_lines.set_data(np.zeros((0,3), dtype=np.float32))
            return
        mp = self.placed_apos.astype(np.float32)
        colors = colors_for_enames(self.placed_enames, alpha=0.95)
        sizes  = sizes_for_enames(self.placed_enames, scale=0.25)
        self.mol_markers.set_data(mp, face_color=colors, size=sizes, edge_width=0.5, edge_color='black')
        # Bonds
        if self.placed_bonds is not None and len(self.placed_bonds) > 0:
            segs = np.empty((len(self.placed_bonds)*2, 3), dtype=np.float32)
            segs[0::2] = mp[self.placed_bonds[:,0]]
            segs[1::2] = mp[self.placed_bonds[:,1]]
            connect = np.zeros(len(segs), dtype=bool)
            connect[0::2] = True  # connect each pair
            self.bond_lines.set_data(segs, connect=connect, color=(0.3,0.3,0.3,0.8), width=1.5)
        else:
            self.bond_lines.set_data(np.zeros((0,3), dtype=np.float32))

    # ================================================================
    #  Save / Reset
    # ================================================================

    def _on_save_all(self):
        all_apos = []; all_enames = []; lvec = None
        if self.substrate is not None:
            all_apos.append(self.substrate.apos); all_enames.extend(list(self.substrate.enames)); lvec = self.substrate.lvec
        if self.placed_apos is not None:
            all_apos.append(self.placed_apos); all_enames.extend(self.placed_enames)
        if not all_apos: self.status.setText("Nothing to save"); return
        fname, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save XYZ", "sequence_on_substrate.xyz", "XYZ (*.xyz)")
        if fname:
            save_xyz(fname, np.vstack(all_apos), all_enames, lvec=lvec, comment=f" sequence={self.sequence}")
            self.status.setText(f"Saved {fname}")

    def _on_save_mols(self):
        if self.placed_apos is None: self.status.setText("No molecules"); return
        fname, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save Mols XYZ", "placed_molecules.xyz", "XYZ (*.xyz)")
        if fname:
            save_xyz(fname, self.placed_apos, self.placed_enames, comment=f" sequence={self.sequence}")
            self.status.setText(f"Saved {fname}")

    def _on_reset_view(self):
        self.view.camera.set_state({'elevation': 30, 'azimuth': -60, 'distance': 80})
        self.canvas.update()


# ======================== Main ========================

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Sequence Placer (VisPy) - molecules on NaCl step edge')
    parser.add_argument('--headless', action='store_true', help='Run without GUI (uses original SequencePlacer)')
    args, remaining = parser.parse_known_args()

    if args.headless:
        from pyBall.SequencePlacer import run_headless
        # delegate to original headless
        import argparse as ap2
        p2 = ap2.ArgumentParser()
        p2.add_argument('--mols', type=str, default='')
        p2.add_argument('--seq',  type=str, default='A')
        p2.add_argument('--a',    type=float, default=2.82065)
        p2.add_argument('--nx',   type=int, default=15)
        p2.add_argument('--ny',   type=int, default=8)
        p2.add_argument('--nz',   type=int, default=3)
        p2.add_argument('--nsteps', type=int, default=1)
        p2.add_argument('--row_angle', type=float, default=0.0)
        p2.add_argument('--row_spacing', type=float, default=12.0)
        p2.add_argument('--mol_height', type=float, default=3.0)
        p2.add_argument('--rx', type=float, default=0.0)
        p2.add_argument('--ry', type=float, default=0.0)
        p2.add_argument('--rz', type=float, default=0.0)
        p2.add_argument('--origin_x', type=float, default=0.0)
        p2.add_argument('--origin_y', type=float, default=0.0)
        p2.add_argument('--out', type=str, default='sequence_on_substrate.xyz')
        p2.add_argument('--out_mols', type=str, default='placed_molecules.xyz')
        a2 = p2.parse_args(remaining)
        mol_files = {}
        if a2.mols:
            for part in a2.mols.split(','):
                letter, path = part.strip().split('=', 1)
                mol_files[letter.strip()] = path.strip()
        assert mol_files, "No molecules specified"
        run_headless(mol_files, sequence=a2.seq, a=a2.a, nx=a2.nx, ny=a2.ny, nz=a2.nz,
                     nsteps=a2.nsteps, row_angle=a2.row_angle, row_spacing=a2.row_spacing,
                     mol_height=a2.mol_height, rx=a2.rx, ry=a2.ry, rz=a2.rz,
                     origin_x=a2.origin_x, origin_y=a2.origin_y,
                     out_combined=a2.out, out_mols=a2.out_mols)
    else:
        app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
        win = SequencePlacerVisPy()
        win.show()
        win._on_build_sub()   # build substrate + initial placement on startup
        app.exec_()

if __name__ == '__main__':
    main()
