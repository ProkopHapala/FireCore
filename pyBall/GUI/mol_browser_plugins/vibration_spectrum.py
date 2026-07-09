# === AUTO-DOC BEGIN ===
"""Vibration spectrum side panel for VispyMolBrowser.

Displays a zoomable FTIR histogram or mode-stick plot and overlays selected
vibrational eigenmodes on the 3D molecular structure (arrows + animation).

Two data source formats are supported (auto-detected, NPZ tried first):

  1. **Nanocrystal pipeline NPZ** (original path):
     - ``05_spectrum.npz``: keys ``omega_centers``, ``hist``, ``omegas_modes`` (internal units: sqrt(eV/amu)/Å)
     - ``04_hessian.npz``: keys ``K`` or ``K_projected``, ``M``, ``pos``, ``Z``, ``natoms``
     - Frequencies in internal units; converted to cm⁻¹ via ``omega_internal_to_cm1()``
     - Eigenvectors re-diagonalized from Hessian on load (not stored in 05 v1.2)

  2. **PySCF loose .npy** (added 2026-07):
     - ``frequencies_cm1.npy``: shape (3N,) float64 or complex128; already in cm⁻¹
     - ``modes.npy``: shape (n_vib, N, 3) float64; n_vib = 3N-6 (vib only, no trans/rot)
     - ``masses.npy``: shape (N,) — optional, not used by plugin
     - ``hessian.npy``: shape (N, N, 3, 3) float64 — optional, not used by plugin
     - ``relaxed.xyz``: standard XYZ geometry (loaded by browser, not plugin)
     - Frequencies already in cm⁻¹; ``_is_cm1`` flag skips unit conversion
     - Modes reshaped to (3N, n_vib) for compatibility with ``displacement_for_mode()``

Selection logic: ``on_molecule_selected`` tries NPZ stages first; if none found,
falls back to PySCF ``.npy`` files in the molecule directory.
"""
# === AUTO-DOC END ===
from __future__ import annotations

import os

import numpy as np
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from PyQt5 import QtCore, QtWidgets

from pyBall.GUI.mol_browser_plugins.base import MolBrowserContext, MolBrowserPlugin
from pyBall.io.crystal_npz import find_nanocrystal_pipeline_stages, pipeline_dir_for_molecule_path
from pyBall.nanocrystal_pipeline import load_spectrum_npz, omega_internal_to_cm1, solve_normal_modes_from_hessian_npz, vibrational_mode_mask

_VIEW_KINDS = (
    ('hist', 'Histogram'),
    ('sticks', 'Mode sticks (vib)'),
    ('probe_x', 'Probe |e·u|² (x)'),
    ('probe_y', 'Probe |e·u|² (y)'),
    ('probe_z', 'Probe |e·u|² (z)'),
)


class VibrationSpectrumPanel(QtWidgets.QWidget):
    mode_selected = QtCore.pyqtSignal(int, float)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._spectrum = None
        self._normal_modes = None
        self._view_kind = 'hist'
        self._mode_global_idx = None
        self._mode_marker = None
        v = QtWidgets.QVBoxLayout(self)
        v.setContentsMargins(4, 4, 4, 4)
        row = QtWidgets.QHBoxLayout()
        self.cmb_view = QtWidgets.QComboBox()
        for kid, label in _VIEW_KINDS:
            self.cmb_view.addItem(label, kid)
        self.cmb_view.currentIndexChanged.connect(self._on_view_changed)
        row.addWidget(QtWidgets.QLabel('Plot:'))
        row.addWidget(self.cmb_view, stretch=1)
        v.addLayout(row)
        self.fig = Figure(figsize=(4.2, 3.0), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvas(self.fig)
        self.toolbar = NavigationToolbar(self.canvas, self)
        v.addWidget(self.toolbar)
        v.addWidget(self.canvas, stretch=1)
        self.canvas.mpl_connect('button_press_event', self._on_canvas_click)
        ctrl = QtWidgets.QHBoxLayout()
        self.btn_prev = QtWidgets.QPushButton('◀')
        self.btn_next = QtWidgets.QPushButton('▶')
        self.spin_mode = QtWidgets.QSpinBox()
        self.spin_mode.setMinimum(0)
        self.spin_mode.setMaximum(0)
        self.lbl_mode = QtWidgets.QLabel('ω = —')
        self.btn_prev.clicked.connect(lambda: self._step_mode(-1))
        self.btn_next.clicked.connect(lambda: self._step_mode(+1))
        self.spin_mode.valueChanged.connect(self._on_spin_mode)
        ctrl.addWidget(self.btn_prev)
        ctrl.addWidget(self.spin_mode)
        ctrl.addWidget(self.btn_next)
        ctrl.addWidget(self.lbl_mode, stretch=1)
        v.addLayout(ctrl)
        vis = QtWidgets.QHBoxLayout()
        self.cb_arrows = QtWidgets.QCheckBox('Arrows')
        self.cb_arrows.setChecked(True)
        self.cb_animate = QtWidgets.QCheckBox('Animate')
        self.cb_animate.setChecked(False)
        self.sld_amp = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.sld_amp.setMinimum(1)
        self.sld_amp.setMaximum(100)
        self.sld_amp.setValue(25)
        vis.addWidget(self.cb_arrows)
        vis.addWidget(self.cb_animate)
        vis.addWidget(QtWidgets.QLabel('Amp'))
        vis.addWidget(self.sld_amp, stretch=1)
        v.addLayout(vis)
        self.cb_arrows.stateChanged.connect(self._emit_mode_vis)
        self.cb_animate.stateChanged.connect(self._emit_mode_vis)
        self.sld_amp.valueChanged.connect(self._emit_mode_vis)
        self.lbl_status = QtWidgets.QLabel('')
        self.lbl_status.setWordWrap(True)
        v.addWidget(self.lbl_status)
        self._clear_plot('Select a relaxed crystal with 05_spectrum.npz')

    def _clear_plot(self, msg=''):
        self.ax.clear()
        self.ax.set_title('Vibrational spectrum')
        self.ax.text(0.5, 0.5, msg or 'no data', ha='center', va='center', transform=self.ax.transAxes, fontsize=9)
        self.ax.set_xlabel('Wavenumber (cm$^{-1}$)')
        self.canvas.draw_idle()
        self.lbl_status.setText(msg)

    def set_bundle(self, stages: dict):
        self._spectrum = None
        self._normal_modes = None
        self._mode_global_idx = None
        spec_p = stages.get('spectrum')
        hess_p = stages.get('hessian')
        if spec_p is None:
            self._clear_plot('Missing 05_spectrum.npz in crystal directory')
            self.mode_selected.emit(-1, 0.0)
            return
        self._spectrum = load_spectrum_npz(spec_p)
        if hess_p is not None:
            self._normal_modes = solve_normal_modes_from_hessian_npz(hess_p)
        n_vib = int(self._normal_modes['vib_indices'].size) if self._normal_modes is not None else 0
        self.lbl_status.setText(f"spectrum: {os.path.basename(spec_p)}" + (f"  modes: {n_vib} vib" if n_vib else "  (no 04_hessian — pick disabled)"))
        self._redraw_plot()
        if n_vib > 0:
            self._select_vib_slot(0)

    def set_pyscf_bundle(self, dir_path: str):
        """Load PySCF loose .npy vibration files (frequencies_cm1.npy + modes.npy)."""
        self._spectrum = None
        self._normal_modes = None
        self._mode_global_idx = None
        d = os.path.abspath(dir_path)
        freq_p = os.path.join(d, 'frequencies_cm1.npy')
        modes_p = os.path.join(d, 'modes.npy')
        if not os.path.isfile(freq_p):
            self._clear_plot(f'No frequencies_cm1.npy in {os.path.basename(d)}')
            self.mode_selected.emit(-1, 0.0)
            return
        freqs = np.load(freq_p)
        # Handle complex frequencies (imaginary = unstable modes)
        if np.iscomplexobj(freqs):
            omegas_cm = np.where(freqs.imag > 1e-3, -freqs.imag, freqs.real)
        else:
            omegas_cm = freqs.real.astype(np.float64)
        # Vibrational mode mask: real, positive, above threshold (skip translation/rotation)
        vib_mask = omegas_cm > 10.0
        vib_indices = np.flatnonzero(vib_mask)
        omegas_cm_vib = omegas_cm[vib_mask]
        # Load modes if available — PySCF modes shape: (n_vib, N, 3), n_vib = 3N-6
        modes_cart = None
        if os.path.isfile(modes_p):
            modes = np.load(modes_p)  # (n_vib, N, 3)
            n_vib_m, N_atoms, _ = modes.shape
            modes_cart = modes.reshape(n_vib_m, N_atoms * 3).T  # (3N, n_vib): col i = mode i
        if modes_cart is None:
            self._normal_modes = None
        else:
            self._normal_modes = {
                'path': d,
                'omegas_cm': omegas_cm,
                'omegas_cm_vib': omegas_cm_vib,
                'vib_indices': vib_indices,
                'modes_cart': modes_cart,
            }
        # Build a simple histogram spectrum from frequencies
        xmax = max(omegas_cm_vib.max() * 1.1, 100.0) if len(omegas_cm_vib) else 100.0
        bins = np.arange(0, xmax + 10, 10.0)
        hist, edges = np.histogram(omegas_cm_vib, bins=bins)
        centers = 0.5 * (edges[:-1] + edges[1:])
        self._spectrum = {
            'path': d,
            'omega_centers': centers,
            'hist': hist.astype(np.float64),
            'omegas_modes': omegas_cm,
            'omegas_modes_vib': omegas_cm_vib,
            '_is_cm1': True,  # PySCF stores cm^-1 directly, no internal-unit conversion needed
        }
        n_vib = int(vib_indices.size)
        self.lbl_status.setText(f"PySCF: {os.path.basename(d)}  modes: {n_vib} vib" + (f"  ({len(omegas_cm) - n_vib} trans/rot/imag)" if len(omegas_cm) != n_vib else ""))
        self._redraw_plot()
        if n_vib > 0:
            self._select_vib_slot(0)

    def _vib_omegas_cm(self):
        if self._spectrum is None:
            return np.zeros(0, dtype=np.float64)
        if 'omegas_modes_vib' in self._spectrum:
            om = np.asarray(self._spectrum['omegas_modes_vib'], dtype=np.float64)
            if self._spectrum.get('_is_cm1'):
                return om  # PySCF path: already in cm^-1
            return omega_internal_to_cm1(om)  # NPZ path: internal units
        om = np.asarray(self._spectrum['omegas_modes'], dtype=np.float64)
        return omega_internal_to_cm1(om[vibrational_mode_mask(om)])

    def _redraw_plot(self):
        self.ax.clear()
        if self._spectrum is None:
            self._clear_plot()
            return
        kind = self._view_kind
        if kind == 'hist':
            x = self._spectrum['omega_centers']
            y = self._spectrum['hist']
            self.ax.plot(x, y, color='#1d4ed8', lw=1.4)
            self.ax.fill_between(x, 0, y, color='#3b82f6', alpha=0.18)
            self.ax.set_ylabel('Mode count')
        elif kind == 'sticks':
            om = self._vib_omegas_cm()
            self.ax.vlines(om, 0, 1, colors='#2563eb', linewidth=0.6, alpha=0.75)
            self.ax.set_ylim(0, 1.05)
            self.ax.set_yticks([])
            self.ax.set_ylabel('modes')
        else:
            pol = kind.split('_', 1)[1]
            key = f'probe_weight_{pol}'
            if key not in self._spectrum or self._normal_modes is None:
                self.ax.text(0.5, 0.5, f'missing {key}', ha='center', va='center', transform=self.ax.transAxes)
            else:
                om = np.asarray(self._normal_modes['omegas_cm_vib'], dtype=np.float64)
                w = np.asarray(self._spectrum[key], dtype=np.float64)[self._normal_modes['vib_indices']]
                idx = np.argsort(om)
                self.ax.vlines(om[idx], 0, w[idx], colors='#7c3aed', linewidth=0.8, alpha=0.85)
                self.ax.set_ylabel('|e·u|²')
        self.ax.set_xlabel('Wavenumber (cm$^{-1}$)')
        self.ax.set_title('Vibrational spectrum')
        self.ax.grid(True, alpha=0.25, ls=':')
        self._mode_marker = None
        if self._mode_global_idx is not None and self._normal_modes is not None:
            om_cm = float(self._normal_modes['omegas_cm'][self._mode_global_idx])
            self._mode_marker = self.ax.axvline(om_cm, color='#dc2626', lw=1.8, ls='--')
        self.canvas.draw_idle()

    def _on_view_changed(self, _idx):
        self._view_kind = str(self.cmb_view.currentData())
        self._redraw_plot()

    def _on_canvas_click(self, event):
        if event.inaxes != self.ax or event.button != 1 or event.xdata is None:
            return
        if self._normal_modes is None:
            return
        om_vib = self._normal_modes['omegas_cm_vib']
        if om_vib.size == 0:
            return
        slot = int(np.argmin(np.abs(om_vib - float(event.xdata))))
        self._select_vib_slot(slot)

    def _select_vib_slot(self, slot: int):
        if self._normal_modes is None:
            return
        vib_idx = np.asarray(self._normal_modes['vib_indices'], dtype=np.int32)
        if vib_idx.size == 0:
            return
        slot = int(np.clip(slot, 0, vib_idx.size - 1))
        gidx = int(vib_idx[slot])
        self._mode_global_idx = gidx
        self.spin_mode.blockSignals(True)
        self.spin_mode.setMaximum(max(0, vib_idx.size - 1))
        self.spin_mode.setValue(slot)
        self.spin_mode.blockSignals(False)
        om_cm = float(self._normal_modes['omegas_cm'][gidx])
        self.lbl_mode.setText(f'ω = {om_cm:.1f} cm⁻¹  (#{slot + 1}/{vib_idx.size})')
        self._redraw_plot()
        self.mode_selected.emit(gidx, om_cm)

    def _step_mode(self, delta: int):
        if self._normal_modes is None:
            return
        n = int(self._normal_modes['vib_indices'].size)
        if n == 0:
            return
        cur = self.spin_mode.value()
        self._select_vib_slot((cur + delta) % n)

    def _on_spin_mode(self, val: int):
        self._select_vib_slot(int(val))

    def _emit_mode_vis(self):
        if self._mode_global_idx is not None:
            om_cm = float(self._normal_modes['omegas_cm'][self._mode_global_idx]) if self._normal_modes is not None else 0.0
            self.mode_selected.emit(self._mode_global_idx, om_cm)

    def vis_options(self):
        return {
            'show_arrows': self.cb_arrows.isChecked(),
            'animate': self.cb_animate.isChecked(),
            'amplitude': float(self.sld_amp.value()) * 0.01,
        }

    def displacement_for_mode(self, global_idx: int):
        if self._normal_modes is None:
            return None
        vec = np.asarray(self._normal_modes['modes_cart'][:, int(global_idx)], dtype=np.float64)
        n3 = vec.reshape(-1, 3)
        norm = float(np.linalg.norm(n3))
        if norm > 1e-12:
            n3 = n3 / norm
        return n3


def _find_pyscf_vib_dir(path: str) -> str | None:
    """Return ``path`` if it contains ``frequencies_cm1.npy``, else None."""
    if not path or not os.path.isdir(path):
        return None
    freq_p = os.path.join(path, 'frequencies_cm1.npy')
    if os.path.isfile(freq_p):
        return path
    return None


class VibrationSpectrumPlugin(MolBrowserPlugin):
    plugin_id = 'vibration_spectrum'
    title = 'Vibration'
    priority = 10

    def __init__(self):
        self._panel = None
        self._ctx = None
        self._mode_idx = None

    def is_relevant(self, ctx: MolBrowserContext) -> bool:
        # Original NPZ pipeline path
        stages = ctx.pipeline_stages()
        if 'spectrum' in stages or 'hessian' in stages:
            return True
        # PySCF loose .npy path
        mol_dir = pipeline_dir_for_molecule_path(ctx.selected_path or ctx.directory)
        if _find_pyscf_vib_dir(mol_dir) is not None:
            return True
        if _find_pyscf_vib_dir(ctx.directory) is not None:
            return True
        return False

    def create_panel(self, parent: QtWidgets.QWidget) -> QtWidgets.QWidget:
        self._panel = VibrationSpectrumPanel(parent)
        self._panel.mode_selected.connect(self._on_mode_selected)
        self._panel.cb_arrows.stateChanged.connect(self._apply_viewer_mode)
        self._panel.cb_animate.stateChanged.connect(self._apply_viewer_mode)
        self._panel.sld_amp.valueChanged.connect(self._apply_viewer_mode)
        return self._panel

    def on_directory_changed(self, ctx: MolBrowserContext):
        self._ctx = ctx

    def on_molecule_selected(self, ctx: MolBrowserContext):
        self._ctx = ctx
        if self._panel is None:
            return
        mol_dir = pipeline_dir_for_molecule_path(ctx.selected_path or ctx.directory)
        # Try NPZ pipeline stages first (original path)
        stages = find_nanocrystal_pipeline_stages(mol_dir)
        if 'spectrum' in stages or 'hessian' in stages:
            self._panel.set_bundle(stages)
        else:
            # Try PySCF loose .npy format
            pyscf_dir = _find_pyscf_vib_dir(mol_dir) or _find_pyscf_vib_dir(ctx.directory)
            if pyscf_dir is not None:
                self._panel.set_pyscf_bundle(pyscf_dir)
            else:
                self._panel.set_bundle(stages)  # will show 'missing' message
        self._mode_idx = None
        self._clear_viewer_mode()

    def on_deactivate(self):
        self._clear_viewer_mode()

    def _on_mode_selected(self, global_idx: int, _omega_cm: float):
        self._mode_idx = None if global_idx < 0 else int(global_idx)
        self._apply_viewer_mode()

    def _apply_viewer_mode(self):
        if self._ctx is None or self._panel is None or self._mode_idx is None:
            self._clear_viewer_mode()
            return
        disp = self._panel.displacement_for_mode(self._mode_idx)
        if disp is None:
            self._clear_viewer_mode()
            return
        opts = self._panel.vis_options()
        self._ctx.viewer.set_vibration_mode(disp, amplitude=opts['amplitude'], show_arrows=opts['show_arrows'], animate=opts['animate'])

    def _clear_viewer_mode(self):
        if self._ctx is not None:
            self._ctx.viewer.clear_vibration_mode()
