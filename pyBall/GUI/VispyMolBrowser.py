#!/usr/bin/env python3
"""VispyMolBrowser.py — PyQt5 + Vispy ACDsee-style molecular browser with extensible analysis plugins.

Hosts thumbnail grid (01/02/03 NPZ only), interactive AtomScene 3D view, and MolBrowserPluginHost
east-side panels. VibrationSpectrumPlugin loads 05_spectrum + 04_hessian companions for FTIR plot
and eigenmode arrows/animation. CPU QPainter thumbnails avoid a second GL context.

Guide: doc/Topics/FTIR_Nanocrystals/Python_Vispy_MolBrowser_Plugins.md
Open issues: dense eigh on mode pick for large N; eigenvectors not in 05_spectrum v1.2.
"""
import argparse
import hashlib
import os
import sys

import numpy as np
from PyQt5 import QtCore, QtGui, QtWidgets

import vispy
vispy.use('pyqt5')

from pyBall.AtomicSystem import AtomicSystem
from pyBall.VispyUtils import AtomScene
from pyBall.MolGUI_common import ELEMENT_COLORS, ELEMENT_SIZES, rotmat_x, rotmat_y
from pyBall.io.crystal_npz import (
    bboxes_to_atomscene, bond_colors_by_stiffness, bonds_from_neigh_idx,
    colors_from_Z, detect_npz_kind, enames_from_Z, extract_bond_k_for_pairs,
    infer_bonds_if_missing, is_viewer_listable_basename, load_crystal_npz, load_npy_crystal, load_topology_npz,
    sizes_from_Z, topology_has_mmffl,
)
from pyBall.GUI.mol_browser_plugins import MolBrowserContext, default_plugin_registry
from pyBall.GUI.mol_browser_plugins.registry import MolBrowserPluginHost

_VIEW_EXTS = {'.xyz', '.mol', '.mol2', '.npz'}
_THUMB_IMG_W, _THUMB_IMG_H = 120, 90
_THUMB_CAPTION_H = 34
_THUMB_TILE_W = _THUMB_IMG_W
_THUMB_TILE_H = _THUMB_IMG_H + _THUMB_CAPTION_H
_THUMB_SIZE = (_THUMB_IMG_W, _THUMB_IMG_H)
_CACHE_DIRNAME = '.vispy_mol_browser_cache_v2'
_CAPTION_CSS = 'font-size:8px; color:#1a1a1a; background:transparent;'
_TILE_CSS = 'QFrame { border:1px solid #999; margin:0; padding:0; background:#fff; }'
_TILE_SEL_CSS = 'QFrame { border:2px solid #2060c0; margin:0; padding:0; background:#fff; }'
_GRID_GAP = 1
_THUMB_COLS = 2


def _repo_root():
    return os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))


def scan_molecule_entries(dir_path):
    """Directory scan: stat only, sorted molecule file paths (lazy load)."""
    entries = []
    if not os.path.isdir(dir_path):
        raise FileNotFoundError(f"scan_molecule_entries: not a directory: {dir_path}")
    for fname in sorted(os.listdir(dir_path)):
        fpath = os.path.join(dir_path, fname)
        ext = os.path.splitext(fname)[1].lower()
        if os.path.isfile(fpath) and ext in _VIEW_EXTS:
            if ext == '.npz' and not is_viewer_listable_basename(fname):
                continue
            entries.append(fpath)
    pos_npy = os.path.join(dir_path, 'pos.npy')
    Z_npy = os.path.join(dir_path, 'Z.npy')
    if os.path.isfile(pos_npy) and os.path.isfile(Z_npy):
        entries.append(dir_path + os.sep)  # sentinel: trailing sep marks npy dir
    return entries


def scan_subdir_entries(dir_path):
    """Immediate child directories (skip hidden + thumbnail cache)."""
    if not os.path.isdir(dir_path):
        return []
    out = []
    for name in sorted(os.listdir(dir_path)):
        if name.startswith('.') or name == _CACHE_DIRNAME:
            continue
        p = os.path.join(dir_path, name)
        if os.path.isdir(p):
            out.append(p)
    return out


def parent_directory(dir_path):
    dir_path = os.path.abspath(dir_path)
    parent = os.path.abspath(os.path.join(dir_path, os.pardir))
    if parent == dir_path:
        return None
    return parent


def _topology_companion(path):
    d = os.path.dirname(path)
    for name in ('03_topology.npz',):
        p = os.path.join(d, name)
        if os.path.isfile(p):
            return p
    return None


class MoleculeData:
    """Loaded geometry + optional topology overlay."""
    __slots__ = ('path', 'pos', 'Z', 'bonds_ij', 'topology', 'enames')

    def __init__(self, path, pos, Z, bonds_ij=None, topology=None):
        self.path = path
        self.pos = np.asarray(pos, dtype=np.float64)
        self.Z = np.asarray(Z, dtype=np.int32).reshape(-1)
        self.bonds_ij = None if bonds_ij is None else np.asarray(bonds_ij, dtype=np.int32)
        self.topology = topology
        self.enames = enames_from_Z(self.Z)

    @property
    def natoms(self):
        return int(self.Z.shape[0])


def load_molecule_entry(path):
    path = os.path.normpath(path)
    if path.endswith(os.sep):
        path = path.rstrip(os.sep)
    if os.path.isdir(path):
        raw = load_npy_crystal(path)
        bonds = infer_bonds_if_missing(raw['pos'], raw['Z'], raw.get('bonds_ij'))
        return MoleculeData(path, raw['pos'], raw['Z'], bonds)
    ext = os.path.splitext(path)[1].lower()
    topology = None
    if ext == '.npz':
        kind = detect_npz_kind(path)
        if kind == 'topology':
            topology = load_topology_npz(path)
            pos = topology.get('pos')
            Z = topology.get('Z')
            if pos is None or Z is None:
                raise KeyError(f"load_molecule_entry: topology NPZ missing pos/Z: {path}")
            if topology_has_mmffl(topology):
                bonds = bonds_from_neigh_idx(topology['neigh_idx'], topology.get('stick_class'))
            else:
                bonds = infer_bonds_if_missing(pos, Z, None)
        else:
            raw = load_crystal_npz(path)
            pos, Z = raw['pos'], raw['Z']
            bonds = infer_bonds_if_missing(pos, Z, raw.get('bonds_ij'))
            topo_path = _topology_companion(path)
            if topo_path:
                topology = load_topology_npz(topo_path)
    elif ext in ('.xyz', '.mol', '.mol2'):
        sys = AtomicSystem(fname=path, bPreinit=True)
        pos = np.asarray(sys.apos, dtype=np.float64)
        Z = np.asarray(sys.atypes, dtype=np.int32)
        bonds = None if sys.bonds is None else np.asarray(sys.bonds, dtype=np.int32)
        bonds = infer_bonds_if_missing(pos, Z, bonds)
    else:
        raise ValueError(f"load_molecule_entry: unsupported extension {ext} ({path})")
    return MoleculeData(path, pos, Z, bonds, topology)


def _fit_camera(scene, pos):
    pos = np.asarray(pos, dtype=np.float32)
    if pos.size == 0:
        return
    ctr = pos.mean(axis=0)
    ext = pos.max(axis=0) - pos.min(axis=0)
    rad = float(max(np.max(ext), 1.0))
    cam = scene.view.camera
    cam.center = tuple(ctr.tolist())
    cam.scale_factor = rad * 1.6
    cam.distance = max(rad * 3.0, 10.0)


def _group_palette(gid):
    h = float((int(gid) * 0.61803398875) % 1.0)
    r = abs(h * 6.0 - 3.0) - 1.0
    g = 2.0 - abs(h * 6.0 - 2.0)
    b = 2.0 - abs(h * 6.0 - 4.0)
    return tuple(np.clip([r, g, b], 0.0, 1.0).tolist()) + (0.95,)


def _thumb_cache_key(path):
    p = os.path.normpath(path.rstrip(os.sep))
    if os.path.isfile(p):
        st = os.stat(p)
        return (p, st.st_mtime_ns, st.st_size)
    if os.path.isdir(p):
        parts = []
        for name in ('pos.npy', 'Z.npy', 'bonds_ij.npy'):
            fp = os.path.join(p, name)
            if os.path.isfile(fp):
                st = os.stat(fp)
                parts.append((name, st.st_mtime_ns, st.st_size))
        return (p, tuple(parts))
    return (p, 0)


def _thumb_cache_file(cache_dir, key):
    h = hashlib.sha256(repr(key).encode('utf-8')).hexdigest()[:20]
    return os.path.join(cache_dir, f'{h}.png')


_THUMB_ELEVATION = 90.0
_THUMB_AZIMUTH = 0.0


def _format_caption_text(name, max_line=14):
    """Break long names for multiline captions (prefer underscores)."""
    name = name.strip()
    if len(name) <= max_line:
        return name
    parts = name.replace('-', '_').split('_')
    lines, cur = [], ''
    for p in parts:
        if not cur:
            cur = p
        elif len(cur) + 1 + len(p) <= max_line:
            cur = cur + '_' + p
        else:
            lines.append(cur)
            cur = p
    if cur:
        lines.append(cur)
    return '\n'.join(lines)


def _project_thumb_xy(pos, elevation=_THUMB_ELEVATION, azimuth=_THUMB_AZIMUTH):
    pos = np.asarray(pos, dtype=np.float64)
    p = pos - pos.mean(axis=0)
    if float(elevation) >= 89.0:
        xy = p[:, :2]
    else:
        R = rotmat_y(azimuth) @ rotmat_x(elevation)
        xy = (R @ p.T).T[:, :2]
    rad = float(max(np.max(np.linalg.norm(xy, axis=1)), 1e-3))
    return xy, rad


def render_molecule_thumbnail_image(path, size=_THUMB_SIZE):
    """Flat top-down CPU thumbnail — minimal margin, no faux-3D."""
    mol = load_molecule_entry(path)
    w, h = int(size[0]), int(size[1])
    xy, rad = _project_thumb_xy(mol.pos)
    margin = 2.0
    scale = (min(w, h) - 2.0 * margin) / (2.05 * rad)
    cx, cy = w * 0.5, h * 0.5
    px = cx + xy[:, 0] * scale
    py = cy - xy[:, 1] * scale
    img = QtGui.QImage(w, h, QtGui.QImage.Format_RGB32)
    img.fill(QtGui.QColor('#ffffff'))
    painter = QtGui.QPainter(img)
    painter.setRenderHint(QtGui.QPainter.Antialiasing, True)
    bonds = mol.bonds_ij
    if bonds is not None and bonds.size > 0:
        pen = QtGui.QPen(QtGui.QColor('#888888'))
        pen.setWidthF(0.8)
        painter.setPen(pen)
        for ia, ib in bonds:
            painter.drawLine(QtCore.QPointF(float(px[ia]), float(py[ia])), QtCore.QPointF(float(px[ib]), float(py[ib])))
    for i, e in enumerate(mol.enames):
        r_pt = max(1.2, float(ELEMENT_SIZES.get(e, 30)) * scale * 0.007)
        painter.setBrush(QtGui.QColor(ELEMENT_COLORS.get(e, '#808080')))
        painter.setPen(QtCore.Qt.NoPen)
        painter.drawEllipse(QtCore.QPointF(float(px[i]), float(py[i])), r_pt, r_pt)
    painter.end()
    return img


class ThumbnailCache:
    """Memory + on-disk PNG cache keyed by path mtime."""

    def __init__(self, browse_dir):
        self.cache_dir = os.path.join(os.path.abspath(browse_dir), _CACHE_DIRNAME)
        os.makedirs(self.cache_dir, exist_ok=True)
        self._mem = {}

    def try_disk_pixmap(self, path):
        key = _thumb_cache_key(path)
        if key in self._mem:
            return self._mem[key]
        png = _thumb_cache_file(self.cache_dir, key)
        if os.path.isfile(png) and os.path.getmtime(png) >= _thumb_source_mtime(path):
            pix = QtGui.QPixmap(png)
            if not pix.isNull():
                self._mem[key] = pix
                return pix
        return None

    def get_pixmap(self, path):
        key = _thumb_cache_key(path)
        if key in self._mem:
            return self._mem[key]
        png = _thumb_cache_file(self.cache_dir, key)
        if os.path.isfile(png) and os.path.getmtime(png) >= _thumb_source_mtime(path):
            pix = QtGui.QPixmap(png)
            if not pix.isNull():
                self._mem[key] = pix
                return pix
        img = render_molecule_thumbnail_image(path)
        pix = QtGui.QPixmap.fromImage(img)
        pix.save(png, 'PNG')
        self._mem[key] = pix
        return pix


def _thumb_source_mtime(path):
    p = os.path.normpath(path.rstrip(os.sep))
    if os.path.isfile(p):
        return os.path.getmtime(p)
    mt = 0.0
    if os.path.isdir(p):
        for name in ('pos.npy', 'Z.npy', 'bonds_ij.npy'):
            fp = os.path.join(p, name)
            if os.path.isfile(fp):
                mt = max(mt, os.path.getmtime(fp))
    return mt


class ThumbnailRenderer(QtCore.QThread):
    """Background CPU thumbnail render; emits QImage (safe across threads)."""
    thumb_ready = QtCore.pyqtSignal(str, object)

    def __init__(self, paths, cache, parent=None):
        super().__init__(parent)
        self.paths = list(paths)
        self.cache = cache

    def run(self):
        for path in self.paths:
            if self.isInterruptionRequested():
                return
            key = _thumb_cache_key(path)
            png = _thumb_cache_file(self.cache.cache_dir, key)
            if os.path.isfile(png) and os.path.getmtime(png) >= _thumb_source_mtime(path):
                img = QtGui.QImage(png)
            else:
                img = render_molecule_thumbnail_image(path)
                img.save(png, 'PNG')
            self.thumb_ready.emit(path, img)


class _CanvasWheelFilter(QtCore.QObject):
    def __init__(self, scene, parent=None):
        super().__init__(parent)
        self.scene = scene

    def eventFilter(self, obj, event):
        if event.type() == QtCore.QEvent.Wheel:
            delta = event.angleDelta().y()
            if delta:
                self.scene._cam_zoom(float(delta) / 120.0)
            return True
        return super().eventFilter(obj, event)


class VispyMolView(QtWidgets.QWidget):
    """Full interactive 3D viewer wrapping AtomScene."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._mol = None
        self._base_colors = None
        self._highlight_gid = None
        self._cache_key = None
        self._camera_initialized = False
        self._reset_camera_next = True
        self._mode_disp = None
        self._mode_amp = 0.25
        self._mode_arrows = False
        self._mode_animate = False
        self._anim_phase = 0.0
        self._anim_timer = QtCore.QTimer(self)
        self._anim_timer.setInterval(33)
        self._anim_timer.timeout.connect(self._anim_tick)
        self.scene = AtomScene(bgcolor='white')
        self.scene.set_lock_top_view(False)
        self.scene.set_view_navigation(True)
        self.scene.set_pick_mode('3d')
        self.scene.set_show_axes(True)
        self.scene.sig_atom_picked.connect(self._on_atom_picked)
        w = self.scene.widget
        w.setFocusPolicy(QtCore.Qt.StrongFocus)
        self._wheel_filter = _CanvasWheelFilter(self.scene, w)
        w.installEventFilter(self._wheel_filter)
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        toggles = QtWidgets.QHBoxLayout()
        self.cb_bonds = QtWidgets.QCheckBox('Bonds'); self.cb_bonds.setChecked(True)
        self.cb_aabb = QtWidgets.QCheckBox('AABB'); self.cb_aabb.setChecked(True)
        self.cb_radii = QtWidgets.QCheckBox('Radii'); self.cb_radii.setChecked(False)
        self.cb_labels = QtWidgets.QCheckBox('Labels'); self.cb_labels.setChecked(False)
        self.cb_bond_k = QtWidgets.QCheckBox('Bond k'); self.cb_bond_k.setChecked(False)
        for w in (self.cb_bonds, self.cb_aabb, self.cb_radii, self.cb_labels, self.cb_bond_k):
            toggles.addWidget(w); w.stateChanged.connect(self._refresh_display)
        toggles.addStretch(1)
        layout.addLayout(toggles)
        layout.addWidget(self.scene.widget, stretch=1)
        self.lbl_info = QtWidgets.QLabel('')
        layout.addWidget(self.lbl_info)

    def request_camera_reset(self):
        self._reset_camera_next = True

    def _save_camera(self):
        cam = self.scene.view.camera
        return {
            'center': tuple(cam.center),
            'azimuth': float(cam.azimuth),
            'elevation': float(cam.elevation),
            'scale_factor': float(getattr(cam, 'scale_factor', 1.0)),
            'distance': float(cam.distance),
        }

    def _restore_camera(self, state, mol_pos):
        cam = self.scene.view.camera
        ctr = np.asarray(mol_pos, dtype=np.float64).mean(axis=0)
        cam.center = tuple(ctr.tolist())
        cam.azimuth = state['azimuth']
        cam.elevation = state['elevation']
        cam.scale_factor = state['scale_factor']
        cam.distance = state['distance']

    def load_molecule(self, mol):
        self.clear_vibration_mode()
        key = (mol.path, os.path.getmtime(mol.path) if os.path.isfile(mol.path) else 0)
        if self._cache_key == key and self._mol is not None:
            return
        cam_state = None if self._reset_camera_next or not self._camera_initialized else self._save_camera()
        self._cache_key = key
        self._mol = mol
        self._highlight_gid = None
        self._base_colors = colors_from_Z(mol.Z)
        bonds = mol.bonds_ij
        if mol.topology is not None and topology_has_mmffl(mol.topology) and (bonds is None or bonds.shape[0] == 0):
            bonds = bonds_from_neigh_idx(mol.topology['neigh_idx'], mol.topology.get('stick_class'))
            mol.bonds_ij = bonds
        self.scene.set_data(mol.pos.astype(np.float32), colors=self._base_colors, sizes=sizes_from_Z(mol.Z, scale=0.35), bonds=bonds)
        if self._reset_camera_next or not self._camera_initialized:
            _fit_camera(self.scene, mol.pos)
            self.scene.view.camera.elevation = 25.0
            self.scene.view.camera.azimuth = -60.0
            self._camera_initialized = True
            self._reset_camera_next = False
        else:
            self._restore_camera(cam_state, mol.pos)
        self.lbl_info.setText(f"{os.path.basename(mol.path)}  N={mol.natoms}  bonds={0 if bonds is None else bonds.shape[0]}")
        self._refresh_display()

    def _refresh_display(self):
        mol = self._mol
        if mol is None:
            return
        colors = self._base_colors.copy() if self._base_colors is not None else colors_from_Z(mol.Z)
        if self._highlight_gid is not None and mol.topology is not None and 'icolGroup' in mol.topology:
            icol = np.asarray(mol.topology['icolGroup'], dtype=np.int32).reshape(-1)
            if icol.shape[0] != mol.natoms:
                icol = None
        else:
            icol = None
        if icol is not None and self._highlight_gid is not None:
            mask = icol == int(self._highlight_gid)
            pal = _group_palette(self._highlight_gid)
            colors[mask] = pal
            colors[~mask] = colors[~mask] * np.array([0.45, 0.45, 0.45, 1.0], dtype=np.float32)
        self.scene._colors = colors; self.scene._colors_base = colors.copy()
        show_bonds = self.cb_bonds.isChecked()
        use_bond_k = self.cb_bond_k.isChecked() and topology_has_mmffl(mol.topology) and mol.bonds_ij is not None and mol.bonds_ij.shape[0] > 0
        self.scene._bonds = None if (not show_bonds or use_bond_k) else mol.bonds_ij
        self.scene.set_show_bboxes(self.cb_aabb.isChecked())
        if mol.topology is not None and self.cb_aabb.isChecked():
            bbmin, bbmax = bboxes_to_atomscene(mol.topology['group_bbox_min'], mol.topology['group_bbox_max'])
            self.scene.set_bboxes(bbmin, bbmax)
        else:
            self.scene._bboxes_min = None; self.scene._bboxes_max = None
        if self.cb_radii.isChecked() and mol.topology is not None and 'radius' in mol.topology:
            self.scene.set_radius(np.asarray(mol.topology['radius'], dtype=np.float32))
            self.scene.set_show_radius(True)
        else:
            self.scene.set_show_radius(False)
        label_mode = 'global' if self.cb_labels.isChecked() else 'none'
        self.scene.set_label_mode(label_mode)
        if use_bond_k:
            k = extract_bond_k_for_pairs(mol.topology['neigh_idx'], mol.topology['bond_k'], mol.bonds_ij, mol.topology.get('stick_class'))
            valid = np.isfinite(k)
            if np.any(valid):
                segs, bcols = bond_colors_by_stiffness(mol.pos, mol.bonds_ij[valid], k[valid])
                self.scene.bond_lines.set_data(np.zeros((0, 3), dtype=np.float32))
                self.scene._line_set('bond_k', self.scene.bond_colored_lines, segs, color=bcols, width=2.5)
            else:
                self.scene.bond_colored_lines.set_data(np.zeros((0, 3), dtype=np.float32))
        else:
            self.scene.bond_colored_lines.set_data(np.zeros((0, 3), dtype=np.float32))
        self.scene._redraw()

    def clear_vibration_mode(self):
        self._anim_timer.stop()
        self._mode_disp = None
        self._mode_arrows = False
        self._mode_animate = False
        if self._mol is not None:
            self.scene._forces = None
            self.scene.set_data(self._mol.pos.astype(np.float32), colors=self._base_colors, sizes=sizes_from_Z(self._mol.Z, scale=0.35), bonds=self._mol.bonds_ij)
            self._refresh_display()

    def set_vibration_mode(self, displacement_N3, amplitude=0.25, show_arrows=True, animate=False):
        if self._mol is None:
            raise RuntimeError('VispyMolView.set_vibration_mode: no molecule loaded')
        disp = np.asarray(displacement_N3, dtype=np.float64).reshape(self._mol.natoms, 3)
        if disp.shape[0] != self._mol.natoms:
            raise ValueError(f'VispyMolView.set_vibration_mode: disp N={disp.shape[0]} != mol N={self._mol.natoms}')
        self._mode_disp = disp
        self._mode_amp = float(amplitude)
        self._mode_arrows = bool(show_arrows)
        self._mode_animate = bool(animate)
        self._anim_phase = 0.0
        if self._mode_animate:
            self._anim_timer.start()
        else:
            self._anim_timer.stop()
        self._apply_vibration_pose(phase=1.0 if not self._mode_animate else 0.0)

    def _apply_vibration_pose(self, phase=0.0):
        if self._mol is None or self._mode_disp is None:
            return
        base = self._mol.pos.astype(np.float32)
        s = float(np.sin(phase)) if self._mode_animate else 1.0
        pos = base + (self._mode_amp * s * self._mode_disp.astype(np.float32))
        forces = (self._mode_amp * self._mode_disp.astype(np.float32)) if self._mode_arrows else None
        self.scene.set_data(pos, colors=self._base_colors, sizes=sizes_from_Z(self._mol.Z, scale=0.35), bonds=self._mol.bonds_ij, forces=forces, force_scale=1.0)
        self._refresh_display()

    def _anim_tick(self):
        if self._mode_disp is None or not self._mode_animate:
            self._anim_timer.stop()
            return
        self._anim_phase += 0.12
        self._apply_vibration_pose(phase=self._anim_phase)

    def wheelEvent(self, ev):
        delta = ev.angleDelta().y()
        if delta:
            self.scene._cam_zoom(float(delta) / 120.0)
            ev.accept()
        else:
            super().wheelEvent(ev)

    def _on_atom_picked(self, idx):
        if self._mol is None or self._mol.topology is None or 'icolGroup' not in self._mol.topology:
            return
        icol = np.asarray(self._mol.topology['icolGroup'], dtype=np.int32).reshape(-1)
        if idx < 0 or idx >= icol.shape[0] or icol.shape[0] != self._mol.natoms:
            return
        gid = int(icol[idx])
        self._highlight_gid = None if gid < 0 else gid
        self._refresh_display()


class NavFolderWidget(QtWidgets.QFrame):
    """Folder tile: parent or subdir."""
    clicked = QtCore.pyqtSignal(str)

    def __init__(self, label, dir_path, parent=None):
        super().__init__(parent)
        self.dir_path = dir_path
        self.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.setFixedSize(_THUMB_TILE_W, _THUMB_TILE_H)
        self.setSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        self.setStyleSheet(_TILE_CSS)
        v = QtWidgets.QVBoxLayout(self); v.setContentsMargins(0, 0, 0, 0); v.setSpacing(0)
        icon = QtWidgets.QLabel('📁')
        icon.setFixedSize(_THUMB_IMG_W, _THUMB_IMG_H)
        icon.setAlignment(QtCore.Qt.AlignCenter)
        icon.setStyleSheet('font-size:28px; background:#eef2f8; border:none;')
        v.addWidget(icon)
        self.caption = QtWidgets.QLabel(_format_caption_text(label))
        self.caption.setAlignment(QtCore.Qt.AlignHCenter | QtCore.Qt.AlignTop)
        self.caption.setWordWrap(True)
        self.caption.setFixedHeight(_THUMB_CAPTION_H)
        self.caption.setStyleSheet(_CAPTION_CSS)
        v.addWidget(self.caption)

    def mousePressEvent(self, ev):
        self.clicked.emit(self.dir_path)
        super().mousePressEvent(ev)


class ThumbnailWidget(QtWidgets.QFrame):
    clicked = QtCore.pyqtSignal(str)

    def __init__(self, path, parent=None):
        super().__init__(parent)
        self.path = path
        self.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.setFixedSize(_THUMB_TILE_W, _THUMB_TILE_H)
        self.setSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        self.setStyleSheet(_TILE_CSS)
        v = QtWidgets.QVBoxLayout(self); v.setContentsMargins(0, 0, 0, 0); v.setSpacing(0)
        self.img = QtWidgets.QLabel()
        self.img.setFixedSize(_THUMB_IMG_W, _THUMB_IMG_H)
        self.img.setAlignment(QtCore.Qt.AlignCenter)
        self.img.setStyleSheet('background:#fff; border:none; margin:0; padding:0;')
        self.img.setScaledContents(True)
        self.img.setText('…')
        v.addWidget(self.img)
        label = os.path.basename(path.rstrip(os.sep)) or path
        self.caption = QtWidgets.QLabel(_format_caption_text(label))
        self.caption.setAlignment(QtCore.Qt.AlignHCenter | QtCore.Qt.AlignTop)
        self.caption.setWordWrap(True)
        self.caption.setFixedHeight(_THUMB_CAPTION_H)
        self.caption.setStyleSheet(_CAPTION_CSS)
        v.addWidget(self.caption)

    def set_selected(self, on):
        self.setStyleSheet(_TILE_SEL_CSS if on else _TILE_CSS)

    def set_pixmap(self, pix):
        self.img.setPixmap(pix)
        self.img.setText('')

    def mousePressEvent(self, ev):
        self.clicked.emit(self.path)
        super().mousePressEvent(ev)


class VispyMolBrowser(QtWidgets.QMainWindow):
    def __init__(self, dir_path=None):
        super().__init__()
        self.setWindowTitle('Vispy Molecular Browser')
        self.resize(1400, 900)
        self._entries = []
        self._thumb_widgets = {}
        self._thumb_cache = None
        self._thumb_worker = None
        self._browse_dir = None
        self._plugin_registry = default_plugin_registry()
        self._plugin_host = MolBrowserPluginHost(self._plugin_registry)
        central = QtWidgets.QWidget(); self.setCentralWidget(central)
        h = QtWidgets.QHBoxLayout(central)
        left = QtWidgets.QVBoxLayout()
        left.setContentsMargins(0, 0, 0, 0)
        left.setSpacing(2)
        nav_row = QtWidgets.QHBoxLayout()
        self.btn_up = QtWidgets.QPushButton('[..]')
        self.btn_up.setToolTip('Parent directory')
        self.btn_up.clicked.connect(self._go_parent)
        self.btn_open = QtWidgets.QPushButton('Open…')
        self.btn_open.setToolTip('Choose any directory')
        self.btn_open.clicked.connect(self._open_directory)
        self.path_edit = QtWidgets.QLineEdit()
        self.path_edit.setPlaceholderText('Directory path — Enter to open')
        self.path_edit.returnPressed.connect(self._open_path_from_edit)
        nav_row.addWidget(self.btn_up)
        nav_row.addWidget(self.btn_open)
        nav_row.addWidget(self.path_edit, stretch=1)
        left.addLayout(nav_row)
        self.scroll = QtWidgets.QScrollArea()
        self.scroll.setWidgetResizable(True)
        self.scroll.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.thumb_host = QtWidgets.QWidget()
        self.thumb_grid = QtWidgets.QGridLayout(self.thumb_host)
        self.thumb_grid.setSpacing(_GRID_GAP)
        self.thumb_grid.setContentsMargins(0, 0, 0, 0)
        self.thumb_grid.setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignTop)
        self.scroll.setWidget(self.thumb_host)
        self.scroll.verticalScrollBar().valueChanged.connect(self._on_scroll)
        left.addWidget(self.scroll, stretch=1)
        h.addLayout(left, stretch=0)
        self.viewer = VispyMolView()
        h.addWidget(self.viewer, stretch=1)
        h.addWidget(self._plugin_host, stretch=0)
        self._plugin_host.setMinimumWidth(360)
        if dir_path:
            self.load_directory(dir_path)

    def _open_directory(self):
        start = self._browse_dir if self._browse_dir else _repo_root()
        d = QtWidgets.QFileDialog.getExistingDirectory(self, 'Select molecule directory', start)
        if d:
            self.load_directory(d)

    def _open_path_from_edit(self):
        p = self.path_edit.text().strip()
        if not p:
            return
        p = os.path.abspath(os.path.expanduser(p))
        if not os.path.isdir(p):
            raise FileNotFoundError(f"VispyMolBrowser: not a directory: {p}")
        self.load_directory(p)

    def _go_parent(self):
        if not self._browse_dir:
            return
        parent = parent_directory(self._browse_dir)
        if parent is not None:
            self.load_directory(parent)

    def load_directory(self, dir_path):
        dir_path = os.path.abspath(dir_path)
        if not os.path.isdir(dir_path):
            raise FileNotFoundError(f"VispyMolBrowser.load_directory: not a directory: {dir_path}")
        self._browse_dir = dir_path
        self.viewer.request_camera_reset()
        if self._thumb_worker is not None and self._thumb_worker.isRunning():
            self._thumb_worker.requestInterruption()
            self._thumb_worker.wait(3000)
        self.setWindowTitle(f'Vispy Molecular Browser — {dir_path}')
        self.path_edit.setText(dir_path)
        parent = parent_directory(dir_path)
        self.btn_up.setEnabled(parent is not None)
        ctx = MolBrowserContext(directory=dir_path, selected_path=None, molecule=None, viewer=self.viewer)
        self._entries = self._plugin_registry.filter_directory_entries(scan_molecule_entries(dir_path), ctx)
        subdirs = scan_subdir_entries(dir_path)
        self._thumb_cache = ThumbnailCache(dir_path)
        while self.thumb_grid.count():
            item = self.thumb_grid.takeAt(0)
            if item.widget():
                item.widget().deleteLater()
        self._thumb_widgets = {}
        self._nav_widgets = []
        idx = 0
        cols = _THUMB_COLS
        if parent is not None:
            nw = NavFolderWidget('..', parent)
            nw.clicked.connect(self.load_directory)
            self.thumb_grid.addWidget(nw, idx // cols, idx % cols)
            self._nav_widgets.append(nw)
            idx += 1
        for sub in subdirs:
            nw = NavFolderWidget(os.path.basename(sub), sub)
            nw.clicked.connect(self.load_directory)
            self.thumb_grid.addWidget(nw, idx // cols, idx % cols)
            self._nav_widgets.append(nw)
            idx += 1
        for path in self._entries:
            tw = ThumbnailWidget(path)
            tw.clicked.connect(self._on_thumb_clicked)
            self.thumb_grid.addWidget(tw, idx // cols, idx % cols)
            self._thumb_widgets[path] = tw
            pix = self._thumb_cache.try_disk_pixmap(path)
            if pix is not None:
                tw.set_pixmap(pix)
            idx += 1
        self._tight_pack_thumb_grid(idx)
        ctx = MolBrowserContext(directory=dir_path, selected_path=None, molecule=None, viewer=self.viewer)
        self._plugin_registry.notify_directory_changed(ctx, self._plugin_host)
        if self._entries:
            self._select_path(self._entries[0])
        else:
            self.viewer.lbl_info.setText(f'(no molecules in {dir_path})')
        QtCore.QTimer.singleShot(100, lambda: self._start_thumb_render(self._entries))

    def _tight_pack_thumb_grid(self, n_items):
        """Prevent QGridLayout from stretching columns and wasting horizontal space."""
        cols = _THUMB_COLS
        nrows = max(1, (n_items + cols - 1) // cols)
        for c in range(cols):
            self.thumb_grid.setColumnStretch(c, 0)
            self.thumb_grid.setColumnMinimumWidth(c, _THUMB_TILE_W)
        for r in range(nrows):
            self.thumb_grid.setRowStretch(r, 0)
            self.thumb_grid.setRowMinimumHeight(r, _THUMB_TILE_H)
        host_w = cols * _THUMB_TILE_W + max(0, cols - 1) * _GRID_GAP
        host_h = nrows * _THUMB_TILE_H + max(0, nrows - 1) * _GRID_GAP
        self.thumb_host.setMinimumWidth(host_w)
        self.thumb_host.setMaximumWidth(host_w)
        self.thumb_host.setMinimumHeight(host_h)
        self.scroll.setMinimumWidth(host_w + 2)
        self.scroll.setMaximumWidth(host_w + 2)

    def _start_thumb_render(self, paths):
        pending = [p for p in paths if self._thumb_widgets[p].img.pixmap() is None or self._thumb_widgets[p].img.pixmap().isNull()]
        if not pending:
            return
        self._thumb_worker = ThumbnailRenderer(pending, self._thumb_cache, self)
        self._thumb_worker.thumb_ready.connect(self._on_thumb_ready)
        self._thumb_worker.start()

    def _on_thumb_ready(self, path, img):
        pix = QtGui.QPixmap.fromImage(img)
        if self._thumb_cache is not None:
            self._thumb_cache._mem[_thumb_cache_key(path)] = pix
        tw = self._thumb_widgets.get(path)
        if tw is not None:
            tw.set_pixmap(pix)

    def _on_scroll(self, _val):
        pass

    def _on_thumb_clicked(self, path):
        self._select_path(path)

    def _select_path(self, path):
        for p, tw in self._thumb_widgets.items():
            tw.set_selected(p == path)
        mol = load_molecule_entry(path)
        self.viewer.load_molecule(mol)
        ctx = MolBrowserContext(directory=self._browse_dir, selected_path=path, molecule=mol, viewer=self.viewer)
        self._plugin_registry.notify_molecule_selected(ctx, self._plugin_host)


def main(argv=None):
    parser = argparse.ArgumentParser(description='Vispy molecular browser (xyz/mol/mol2/npz)')
    parser.add_argument('--dir', default=None, help='Directory with molecule files')
    args = parser.parse_args(argv)
    app = QtWidgets.QApplication(sys.argv)
    default_dir = os.path.join(_repo_root(), 'tests', 'tSiNCs', 'fixtures', 'npz_viewer')
    if args.dir is None and os.path.isdir(default_dir):
        args.dir = default_dir
    win = VispyMolBrowser(dir_path=args.dir)
    win.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
