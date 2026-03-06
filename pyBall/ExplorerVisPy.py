#!/usr/bin/env python3
"""
ExplorerVisPy.py - VisPy+PyQt5 GUI for interactive molecule-substrate interaction scanning.
Based on SequencePlacerVisPy.py pattern. Uses InteractionScanner backend for GPU energy evaluation.

Usage:
    python ExplorerVisPy.py                              # default H2O on NaCl
    python ExplorerVisPy.py --mol PTCDA.xyz --sub NaCl_8x8_L3.xyz --type_map C=C_R
"""

import numpy as np
import os, sys, argparse

_THIS_DIR = os.path.abspath(os.path.dirname(__file__))
_ROOT_DIR = os.path.abspath(os.path.join(_THIS_DIR, os.pardir))
for _p in (_ROOT_DIR, _THIS_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from pyBall.MolGUI_common import (
    ELEMENT_COLORS, ELEMENT_SIZES,
    colors_for_enames, sizes_for_enames,
    find_bonds,
)

from PyQt5 import QtWidgets, QtCore
import vispy
vispy.use('pyqt5')
from vispy import scene
from vispy.scene import visuals

import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure


class ExplorerVisPy(QtWidgets.QMainWindow):
    """Interactive molecule-substrate interaction explorer with VisPy 3D view and scan controls."""

    def __init__(self, mol_file=None, sub_file=None, mol_type_map=None, sub_type_map=None):
        super().__init__()
        self.setWindowTitle('Molecule-Substrate Explorer (VisPy)')
        self.resize(1600, 950)

        # ---- Scanner backend (lazy init) ----
        self.scanner = None
        self._mol_file = mol_file
        self._sub_file = sub_file
        self._mol_type_map = mol_type_map or {}
        self._sub_type_map = sub_type_map or {}

        # ---- Molecule pose state ----
        self.mol_pos = np.array([0.0, 0.0, 3.0])   # translation (x, y, z)
        self.mol_rot = np.array([0.0, 0.0, 0.0])   # Euler angles (deg)

        # ---- Host data ----
        self.mol_apos   = None   # (nmol, 3) base positions
        self.mol_enames = None
        self.sub_apos   = None   # (nsub, 3)
        self.sub_enames = None
        self.mol_bonds  = None
        self.sub_bonds  = None

        # ---- Scan cache for visualization ----
        self.last_scan_positions = None  # (N,3) scan points in world coords

        # ---- Scan results ----
        self.last_results = None

        self._build_ui()
        self._init_scanner()

    # ================================================================
    #  UI setup
    # ================================================================

    def _build_ui(self):
        central = QtWidgets.QWidget(); self.setCentralWidget(central)
        hbox = QtWidgets.QHBoxLayout(central)

        # ---- Left: VisPy 3D canvas ----
        self.canvas = scene.SceneCanvas(keys='interactive', bgcolor='white', show=False)
        self.view = self.canvas.central_widget.add_view()
        # Orthographic camera by default (Turntable with fov=0 behaves orthographic)
        self.view.camera = scene.TurntableCamera(fov=0, distance=60, elevation=30, azimuth=-60)
        # Use 3D-looking atom markers with larger size and alpha to reduce z-fighting
        self.sub_markers = visuals.Markers(parent=self.view.scene)
        self.mol_markers = visuals.Markers(parent=self.view.scene)
        self.bond_lines  = visuals.Line(parent=self.view.scene, color=(0.5,0.5,0.5,0.4), width=1.5, antialias=True, method='gl')
        self.mol_bond_lines = visuals.Line(parent=self.view.scene, color=(0.8,0.8,0.2,0.8), width=2.0, antialias=True, method='gl')
        # Scan trajectory visuals
        self.scan_points = visuals.Markers(parent=self.view.scene)
        self.scan_path   = visuals.Line(parent=self.view.scene, color=(0.9,0.2,0.2,0.7), width=2.0, antialias=True, method='gl')

        # ---- Middle: Matplotlib plot ----
        self.fig = Figure(figsize=(4, 6), dpi=90, facecolor='white')
        self.fig_canvas = FigureCanvas(self.fig)
        self.fig_canvas.setMinimumWidth(200)
        self.ax = self.fig.add_subplot(111)
        self.ax.set_facecolor('white')
        self.ax.tick_params(colors='black')
        self.ax.xaxis.label.set_color('black')
        self.ax.yaxis.label.set_color('black')
        self.ax.title.set_color('black')
        for spine in self.ax.spines.values(): spine.set_edgecolor('#444')
        self._last_colorbar = None

        # ---- Right: controls ----
        scroll = QtWidgets.QScrollArea()
        scroll.setWidgetResizable(True)
        ctrl = QtWidgets.QWidget()
        vbox = QtWidgets.QVBoxLayout(ctrl)

        # --- Molecule Pose ---
        grp_pose = QtWidgets.QGroupBox("Molecule Pose")
        gl = QtWidgets.QGridLayout(grp_pose)
        self.sp_tx = self._add_dspin(gl, 0, "X (Å)", 0.0, 0.5, -100, 100, 2)
        self.sp_ty = self._add_dspin(gl, 1, "Y (Å)", 0.0, 0.5, -100, 100, 2)
        self.sp_tz = self._add_dspin(gl, 2, "Z (Å)", 3.0, 0.1, -200, 200, 2)
        self.sp_rx = self._add_dspin(gl, 3, "Rx (°)", 0, 5, -360, 360, 1)
        self.sp_ry = self._add_dspin(gl, 4, "Ry (°)", 0, 5, -360, 360, 1)
        self.sp_rz = self._add_dspin(gl, 5, "Rz (°)", 0, 5, -360, 360, 1)
        vbox.addWidget(grp_pose)

        # --- Physics ---
        grp_phys = QtWidgets.QGroupBox("Physics")
        pfl = QtWidgets.QVBoxLayout(grp_phys)
        self.cb_lj   = QtWidgets.QCheckBox("LJ (van der Waals)"); self.cb_lj.setChecked(True)
        self.cb_coul = QtWidgets.QCheckBox("Coulomb");            self.cb_coul.setChecked(True)
        self.cb_hb   = QtWidgets.QCheckBox("H-bond");             self.cb_hb.setChecked(False)
        self.cb_morse = QtWidgets.QCheckBox("Use Morse (instead of LJ)"); self.cb_morse.setChecked(False)
        pfl.addWidget(self.cb_lj); pfl.addWidget(self.cb_coul); pfl.addWidget(self.cb_hb); pfl.addWidget(self.cb_morse)
        vbox.addWidget(grp_phys)

        # --- Scan Controls ---
        grp_scan = QtWidgets.QGroupBox("Scan")
        sgl = QtWidgets.QGridLayout(grp_scan)
        sgl.addWidget(QtWidgets.QLabel("Type"), 0, 0)
        self.combo_scan = QtWidgets.QComboBox()
        self.combo_scan.addItems([ "Lateral XY", "Z approach","XZ slice", "Rotation", "Rot vs Z"])
        sgl.addWidget(self.combo_scan, 0, 1)
        self.sp_npts = self._add_ispin(sgl, 1, "N points", 60, 10, 10, 500)
        self.sp_zlo  = self._add_dspin(sgl, 2, "Z min (Å)", 2.0, 0.2, -200, 200, 1)
        self.sp_zhi  = self._add_dspin(sgl, 3, "Z max (Å)", 8.0, 0.5, -200, 200, 1)
        self.sp_xylo = self._add_dspin(sgl, 4, "XY min (Å)", -5.0, 1.0, -50, 50, 1)
        self.sp_xyhi = self._add_dspin(sgl, 5, "XY max (Å)", 25.0, 1.0, -50, 50, 1)
        self.sp_nxy  = self._add_ispin(sgl, 6, "N XY", 40, 5, 5, 200)
        self.cb_relax = QtWidgets.QCheckBox("Constrained relax"); self.cb_relax.setChecked(False)
        sgl.addWidget(self.cb_relax, 7, 0, 1, 2)
        btn_scan = QtWidgets.QPushButton("Run Scan"); btn_scan.clicked.connect(self._on_run_scan)
        sgl.addWidget(btn_scan, 8, 0, 1, 2)
        vbox.addWidget(grp_scan)

        # --- Relax params ---
        grp_relax = QtWidgets.QGroupBox("Relaxation")
        rgl = QtWidgets.QGridLayout(grp_relax)
        self.sp_springk = self._add_dspin(rgl, 0, "Spring K", 5.0, 1.0, 0.1, 100, 2)
        self.sp_rdt     = self._add_dspin(rgl, 1, "dt", 0.005, 0.001, 0.001, 0.1, 4)
        self.sp_rsteps  = self._add_ispin(rgl, 2, "N steps", 100, 10, 10, 1000)
        vbox.addWidget(grp_relax)

        # --- Plot Settings ---
        grp_plot = QtWidgets.QGroupBox("Plot Settings")
        pgl = QtWidgets.QGridLayout(grp_plot)
        self.cb_sym_cmap = QtWidgets.QCheckBox("Auto symmetric (vmax=-vmin)")
        self.cb_sym_cmap.setChecked(True)
        pgl.addWidget(self.cb_sym_cmap, 0, 0, 1, 2)
        self.cb_fix_clim = QtWidgets.QCheckBox("Fixed color limits")
        self.cb_fix_clim.setChecked(False)
        pgl.addWidget(self.cb_fix_clim, 1, 0, 1, 2)
        self.sp_vmin = self._add_dspin(pgl, 2, "vmin", -0.05, 0.01, -100, 100, 3)
        self.sp_vmax = self._add_dspin(pgl, 3, "vmax",  0.05, 0.01, -100, 100, 3)
        btn_replot = QtWidgets.QPushButton("Replot current data")
        btn_replot.clicked.connect(self._replot_current)
        pgl.addWidget(btn_replot, 4, 0, 1, 2)
        vbox.addWidget(grp_plot)

        # --- Energy display ---
        grp_e = QtWidgets.QGroupBox("Current Energy")
        efl = QtWidgets.QFormLayout(grp_e)
        self.lbl_etot = QtWidgets.QLabel("---"); efl.addRow("Total", self.lbl_etot)
        self.lbl_elj  = QtWidgets.QLabel("---"); efl.addRow("LJ/Morse", self.lbl_elj)
        self.lbl_ecoul = QtWidgets.QLabel("---"); efl.addRow("Coulomb", self.lbl_ecoul)
        self.lbl_ehb  = QtWidgets.QLabel("---"); efl.addRow("H-bond", self.lbl_ehb)
        vbox.addWidget(grp_e)

        # --- Scan path toggle ---
        self.cb_show_scan = QtWidgets.QCheckBox("Show scan trajectory")
        self.cb_show_scan.setChecked(True)
        self.cb_show_scan.stateChanged.connect(self._update_scan_visuals)
        vbox.addWidget(self.cb_show_scan)

        # --- Status ---
        self.status = QtWidgets.QLabel("Initializing..."); vbox.addWidget(self.status)
        vbox.addStretch()

        scroll.setWidget(ctrl)

        # Splitters: left (3D) | middle (plot) | right (controls)
        splitter_right = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
        splitter_right.addWidget(self.fig_canvas)
        splitter_right.addWidget(scroll)
        splitter_right.setStretchFactor(0, 2)
        splitter_right.setStretchFactor(1, 1)

        splitter_main = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
        splitter_main.addWidget(self.canvas.native)
        splitter_main.addWidget(splitter_right)
        splitter_main.setStretchFactor(0, 2)
        splitter_main.setStretchFactor(1, 3)

        hbox.addWidget(splitter_main)

        # Connect pose spinboxes for live update
        for sp in (self.sp_tx, self.sp_ty, self.sp_tz, self.sp_rx, self.sp_ry, self.sp_rz):
            sp.valueChanged.connect(self._on_pose_changed)
        for cb in (self.cb_lj, self.cb_coul, self.cb_hb, self.cb_morse):
            cb.stateChanged.connect(self._on_pose_changed)

        self.cb_sym_cmap.stateChanged.connect(self._replot_current)
        self.cb_fix_clim.stateChanged.connect(self._replot_current)
        self.sp_vmin.valueChanged.connect(self._replot_current)
        self.sp_vmax.valueChanged.connect(self._replot_current)

    # ---- Spin helpers ----
    def _add_dspin(self, grid, row, label, val, step, lo, hi, decimals):
        grid.addWidget(QtWidgets.QLabel(label), row, 0)
        sp = QtWidgets.QDoubleSpinBox()
        sp.setRange(lo, hi); sp.setSingleStep(step); sp.setDecimals(decimals); sp.setValue(val)
        grid.addWidget(sp, row, 1)
        return sp

    def _add_ispin(self, grid, row, label, val, step, lo, hi):
        grid.addWidget(QtWidgets.QLabel(label), row, 0)
        sp = QtWidgets.QSpinBox()
        sp.setRange(lo, hi); sp.setSingleStep(step); sp.setValue(val)
        grid.addWidget(sp, row, 1)
        return sp

    # ================================================================
    #  Scanner init
    # ================================================================

    def _init_scanner(self):
        from pyBall.OCL.InteractionEnergy import InteractionScanner
        self.scanner = InteractionScanner()
        if self._mol_file:
            if not os.path.exists(self._mol_file):
                raise FileNotFoundError(f"Molecule file not found: {self._mol_file}")
            apos, REQs, enames = self.scanner.load_molecule_xyz(self._mol_file, type_map=self._mol_type_map)
            self.mol_apos = apos; self.mol_enames = enames
            self.mol_bonds = find_bonds(apos)
        if self._sub_file:
            if not os.path.exists(self._sub_file):
                raise FileNotFoundError(f"Substrate file not found: {self._sub_file}")
            apos, REQs, enames = self.scanner.load_substrate_xyz(self._sub_file, type_map=self._sub_type_map)
            self.sub_apos = apos; self.sub_enames = enames
            self.sub_bonds = find_bonds(apos)
            # Match headless macro correction defaults used in tests
            self.scanner.nPBC[:] = (4, 4, 0)
            self.scanner.enable_macro = True
            self.scanner._update_macro_from_substrate()
        self._update_visuals()
        self._eval_current_pose()
        self.status.setText("Ready")

    # ================================================================
    #  Rotation matrix from Euler angles
    # ================================================================

    @staticmethod
    def _euler_to_rotmat(rx_deg, ry_deg, rz_deg):
        rx, ry, rz = np.radians(rx_deg), np.radians(ry_deg), np.radians(rz_deg)
        cx, sx = np.cos(rx), np.sin(rx)
        cy, sy = np.cos(ry), np.sin(ry)
        cz, sz = np.cos(rz), np.sin(rz)
        Rx = np.array([[1,0,0],[0,cx,-sx],[0,sx,cx]])
        Ry = np.array([[cy,0,sy],[0,1,0],[-sy,0,cy]])
        Rz = np.array([[cz,-sz,0],[sz,cz,0],[0,0,1]])
        return Rz @ Ry @ Rx

    def _current_rotmat(self):
        return self._euler_to_rotmat(self.sp_rx.value(), self.sp_ry.value(), self.sp_rz.value())

    def _current_translation(self):
        return np.array([self.sp_tx.value(), self.sp_ty.value(), self.sp_tz.value()])

    def _transform_mol(self):
        """Return transformed molecule positions (N,3)."""
        if self.mol_apos is None: return None
        R = self._current_rotmat()
        t = self._current_translation()
        return (self.mol_apos @ R.T) + t

    # ================================================================
    #  Visuals
    # ================================================================

    def _update_visuals(self):
        # Substrate
        if self.sub_apos is not None and len(self.sub_apos) > 0:
            cols = colors_for_enames(self.sub_enames, alpha=0.45)
            szs  = sizes_for_enames(self.sub_enames, scale=0.15)
            self.sub_markers.set_data(
                self.sub_apos.astype(np.float32), face_color=cols, size=szs,
                edge_width=0.0, symbol='o'
            )
            if self.sub_bonds is not None and len(self.sub_bonds) > 0:
                segs = []
                for i, j in self.sub_bonds:
                    segs.append(self.sub_apos[i]); segs.append(self.sub_apos[j])
                if segs:
                    self.bond_lines.set_data(np.array(segs, dtype=np.float32))
        # Molecule
        self._update_mol_visual()

    def _update_mol_visual(self):
        tpos = self._transform_mol()
        if tpos is None or len(tpos) == 0:
            self.mol_markers.set_data(np.zeros((0,3), dtype=np.float32))
            self.mol_bond_lines.set_data(np.zeros((0,3), dtype=np.float32))
            return
        cols = colors_for_enames(self.mol_enames, alpha=0.95)
        szs  = sizes_for_enames(self.mol_enames, scale=0.25)
        self.mol_markers.set_data(
            tpos.astype(np.float32), face_color=cols, size=szs,
            edge_width=0.5, edge_color='black', symbol='o'
        )
        if self.mol_bonds is not None and len(self.mol_bonds) > 0:
            segs = []
            for i, j in self.mol_bonds:
                segs.append(tpos[i]); segs.append(tpos[j])
            if segs:
                self.mol_bond_lines.set_data(np.array(segs, dtype=np.float32))

    # ================================================================
    #  Energy evaluation
    # ================================================================

    def _sync_physics(self):
        if self.scanner is None: return
        self.scanner.enable_LJ      = self.cb_lj.isChecked()
        self.scanner.enable_Coulomb  = self.cb_coul.isChecked()
        self.scanner.enable_HBond   = self.cb_hb.isChecked()
        self.scanner.enable_Morse   = self.cb_morse.isChecked()
        self.scanner.spring_k       = np.float32(self.sp_springk.value())
        self.scanner.relax_dt       = np.float32(self.sp_rdt.value())
        self.scanner.relax_nsteps   = self.sp_rsteps.value()

    def _eval_current_pose(self):
        if self.scanner is None or self.mol_apos is None or self.sub_apos is None: return
        self._sync_physics()
        R = self._current_rotmat()
        t = self._current_translation()
        res = self.scanner.evaluate_single(pos=t, R=R)
        self.lbl_etot.setText(f"{res['total']:.4f} eV")
        self.lbl_elj.setText(f"{res['LJ']:.4f} eV")
        self.lbl_ecoul.setText(f"{res['Coulomb']:.4f} eV")
        self.lbl_ehb.setText(f"{res['HBond']:.4f} eV")

    # ================================================================
    #  Event handlers
    # ================================================================

    def _on_pose_changed(self, _=None):
        self._update_mol_visual()
        self._eval_current_pose()
        self.canvas.update()

    def _on_run_scan(self):
        if self.scanner is None:
            self.status.setText("Scanner not initialized"); return
        self._sync_physics()
        scan_type = self.combo_scan.currentText()
        npts  = self.sp_npts.value()
        zlo   = self.sp_zlo.value()
        zhi   = self.sp_zhi.value()
        xylo  = self.sp_xylo.value()
        xyhi  = self.sp_xyhi.value()
        nxy   = self.sp_nxy.value()
        relax = self.cb_relax.isChecked()
        R = self._current_rotmat()
        pos_xy = (self.sp_tx.value(), self.sp_ty.value())
        self.status.setText(f"Running {scan_type} scan...")
        QtWidgets.QApplication.processEvents()
        try:
            if scan_type == "Z approach":
                results = self.scanner.scan_z(pos_xy=pos_xy, z_range=(zlo, zhi), nz=npts, R=R, relax=relax)
                self._store_scan_positions(results, mode='z')
                self._plot_1d(results, 'z')
            elif scan_type == "Lateral XY":
                z = self.sp_tz.value()
                # Use the main N points control for XY resolution (square grid npts x npts)
                results = self.scanner.scan_lateral(z=z, x_range=(xylo, xyhi), y_range=(xylo, xyhi), nx=npts, ny=npts, R=R, relax=relax)
                self._store_scan_positions(results, mode='xy')
                self._plot_2d(results, 'xy')
            elif scan_type == "XZ slice":
                y = self.sp_ty.value()
                results = self.scanner.scan_xz(y=y, x_range=(xylo, xyhi), z_range=(zlo, zhi), nx=nxy, nz=npts, R=R, relax=relax)
                self._store_scan_positions(results, mode='xz')
                self._plot_2d(results, 'xz')
            elif scan_type == "Rotation":
                pos = (self.sp_tx.value(), self.sp_ty.value(), self.sp_tz.value())
                results = self.scanner.scan_rotation(pos=pos, nrot=npts, relax=relax)
                self._store_scan_positions(results, mode='rot')
                self._plot_1d(results, 'angle')
            elif scan_type == "Rot vs Z":
                results = self.scanner.scan_rot_z(pos_xy=pos_xy, z_range=(zlo, zhi), nz=npts//2, nrot=npts//2, relax=relax)
                self._store_scan_positions(results, mode='rotz')
                self._plot_2d(results, 'rot_z')
            self.last_results = results
            self.status.setText(f"Scan done: {scan_type}")
        except Exception as e:
            self.status.setText(f"Scan error: {e}")
            import traceback; traceback.print_exc()

    # ================================================================
    #  Plotting
    # ================================================================

    def _replot_current(self):
        if self.last_results is None: return
        scan_type = self.combo_scan.currentText()
        if scan_type == "Z approach":
            self._plot_1d(self.last_results, 'z')
        elif scan_type == "Lateral XY":
            self._plot_2d(self.last_results, 'xy')
        elif scan_type == "XZ slice":
            self._plot_2d(self.last_results, 'xz')
        elif scan_type == "Rotation":
            self._plot_1d(self.last_results, 'angle')
        elif scan_type == "Rot vs Z":
            self._plot_2d(self.last_results, 'rot_z')

    def _plot_1d(self, results, xtype='z'):
        self.ax.clear()
        if self._last_colorbar is not None:
            self._last_colorbar.remove()
            self._last_colorbar = None
        info = results.get('scan_info', {})
        n = len(results['total'])
        if xtype == 'z':
            zr = info.get('z_range', (0, n))
            xs = np.linspace(zr[0], zr[1], n)
            xlabel = 'z (Å)'
        elif xtype == 'angle':
            ar = info.get('angle_range', (0, 2*np.pi))
            xs = np.degrees(np.linspace(ar[0], ar[1], n))
            xlabel = 'angle (°)'
        else:
            xs = np.arange(n)
            xlabel = 'step'
        self.ax.plot(xs, results['total'], color='black', lw=2, label='Total')
        self.ax.plot(xs, results['LJ'],      color='#e94560', lw=1, ls='--', label='LJ')
        self.ax.plot(xs, results['Coulomb'],  color='#0f3460', lw=1, ls='--', label='Coul')
        self.ax.plot(xs, results['HBond'],    color='#533483', lw=1, ls='--', label='HB')
        self.ax.set_xlabel(xlabel, color='black')
        self.ax.set_ylabel('Energy (eV)', color='black')
        self.ax.set_title('Interaction Energy', color='black')
        self.ax.legend(facecolor='white', edgecolor='#444', labelcolor='black', fontsize=8)
        self.ax.axhline(0, color='gray', lw=0.5, ls=':')
        self.fig_canvas.draw()

    # ================================================================
    #  Scan visualization (3D)
    # ================================================================

    def _store_scan_positions(self, results, mode='z'):
        """Compute scan point positions in world coordinates for visualization."""
        transforms = results.get('scan_info', {}).get('transforms')
        if transforms is None:
            self.last_scan_positions = None
            self._update_scan_visuals()
            return
        # transforms is (N,12) packed; unpack to positions of molecule origin (assume origin at (0,0,0))
        T = np.array(transforms, dtype=np.float32).reshape(-1, 3, 4)
        pts = T[:, :, 3]  # translations
        self.last_scan_positions = pts
        self._update_scan_visuals()

    def _update_scan_visuals(self):
        if not self.cb_show_scan.isChecked() or self.last_scan_positions is None:
            self.scan_points.set_data(np.zeros((0, 3), dtype=np.float32))
            self.scan_path.set_data(np.zeros((0, 3), dtype=np.float32))
            return
        pts = self.last_scan_positions.astype(np.float32)
        self.scan_points.set_data(pts, face_color=(0.9, 0.2, 0.2, 0.8), size=6.0, edge_width=0.5, edge_color='white')
        self.scan_path.set_data(pts, color=(0.9, 0.2, 0.2, 0.7), width=2.0)

    def _plot_2d(self, results, mode='xy'):
        self.ax.clear()
        if self._last_colorbar is not None:
            self._last_colorbar.remove()
            self._last_colorbar = None
        info = results.get('scan_info', {})
        E = results['total']
        if mode == 'xy':
            nx = info.get('nx', int(np.sqrt(len(E))))
            ny = info.get('ny', nx)
            xr = info.get('x_range', (0, nx))
            yr = info.get('y_range', (0, ny))
            E2d = E.reshape(nx, ny)
            im = self.ax.imshow(E2d.T, origin='lower', aspect='auto', extent=[xr[0], xr[1], yr[0], yr[1]], cmap='RdBu_r')
            self.ax.set_xlabel('x (Å)', color='black')
            self.ax.set_ylabel('y (Å)', color='black')
        elif mode == 'xz':
            shape = info.get('shape')
            if shape is None:
                nx = info.get('nx', int(np.sqrt(len(E))))
                nz = info.get('nz', nx)
            else:
                nx, nz = shape
            xr = info.get('x_range', (0, nx))
            zr = info.get('z_range', (0, nz))
            E2d = E.reshape(nx, nz)
            im = self.ax.imshow(E2d.T, origin='lower', aspect='auto', extent=[xr[0], xr[1], zr[0], zr[1]], cmap='RdBu_r')
            self.ax.set_xlabel('x (Å)', color='black')
            self.ax.set_ylabel('z (Å)', color='black')
        elif mode == 'rot_z':
            shape = info.get('shape')
            if shape is None:
                nz = info.get('nz', 15)
                nrot = info.get('nrot', 18)
            else:
                nz, nrot = shape
            zr = info.get('z_range', (2, 8))
            E2d = E.reshape(nz, nrot)
            im = self.ax.imshow(E2d.T, origin='lower', aspect='auto', extent=[zr[0], zr[1], 0, 360], cmap='RdBu_r')
            self.ax.set_xlabel('z (Å)', color='black')
            self.ax.set_ylabel('angle (°)', color='black')
        else:
            return
        
        valid_E = E[np.isfinite(E)]
        if self.cb_fix_clim.isChecked():
            vmin = self.sp_vmin.value()
            vmax = self.sp_vmax.value()
        elif self.cb_sym_cmap.isChecked():
            vmin = np.min(valid_E) if len(valid_E) > 0 else -1
            if vmin > 0: vmin = -0.01  # ensure some color range if all positive
            vmax = -vmin
        else:
            vmin, vmax = np.percentile(valid_E, [2, 98]) if len(valid_E) > 0 else (0, 1)
            
        im.set_clim(vmin, vmax)
        self.ax.set_title('Interaction Energy (eV)', color='black')
        self._last_colorbar = self.fig.colorbar(im, ax=self.ax, shrink=0.7)
        self.fig_canvas.draw()

    # ================================================================
    #  File loading
    # ================================================================

    def load_molecule(self, fname, type_map=None):
        if self.scanner is None: return
        apos, REQs, enames = self.scanner.load_molecule_xyz(fname, type_map=type_map)
        self.mol_apos = apos; self.mol_enames = enames
        self.mol_bonds = find_bonds(apos)
        self._update_visuals()
        self._eval_current_pose()

    def load_substrate(self, fname, type_map=None):
        if self.scanner is None: return
        apos, REQs, enames = self.scanner.load_substrate_xyz(fname, type_map=type_map)
        self.sub_apos = apos; self.sub_enames = enames
        self.sub_bonds = find_bonds(apos)
        self._update_visuals()
        self._eval_current_pose()


def parse_type_map(s):
    """Parse 'C=C_R,O=O_2' into dict."""
    if not s: return {}
    d = {}
    for part in s.split(','):
        part = part.strip()
        if '=' in part:
            k, v = part.split('=', 1)
            d[k.strip()] = v.strip()
    return d


def main():
    parser = argparse.ArgumentParser(description="Molecule-Substrate Interaction Explorer")
    parser.add_argument('--mol', type=str, default=None, help='Molecule XYZ file')
    parser.add_argument('--sub', type=str, default=None, help='Substrate XYZ file')
    parser.add_argument('--mol_types', type=str, default='', help='Molecule type map (e.g. C=C_R,O=O_2)')
    parser.add_argument('--sub_types', type=str, default='', help='Substrate type map')
    args = parser.parse_args()

    # Defaults
    data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'cpp', 'common_resources', 'xyz')
    mol_file = args.mol or os.path.join(data_dir, 'H2O_O.xyz')
    sub_file = args.sub or os.path.join(data_dir, 'NaCl_8x8_L3.xyz')

    app = QtWidgets.QApplication(sys.argv)
    win = ExplorerVisPy(
        mol_file=mol_file, sub_file=sub_file,
        mol_type_map=parse_type_map(args.mol_types),
        sub_type_map=parse_type_map(args.sub_types),
    )
    win.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
