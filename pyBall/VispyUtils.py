import numpy as np

from PyQt5 import QtCore
import vispy
vispy.use('pyqt5')
from vispy import scene
from vispy.scene import visuals
from vispy.color import Colormap


def _as_f32(x):
    return np.asarray(x, dtype=np.float32)


def compute_bond_colors_by_length(bonds, pos, color_range=(0.0, 1.0)):
    """Compute bond colors based on bond length (blue=short, red=long).
    
    Args:
        bonds: List/array of (ia, ja) bond indices
        pos: (n,3) array of positions
        color_range: (min, max) color values for blue-red mapping
    
    Returns:
        bond_segs: (2*m, 3) array of segment endpoints
        bond_colors: (2*m, 4) array of RGBA colors
    """
    bond_lengths = []
    for b in bonds:
        ia, ja = b
        p1, p2 = pos[ia], pos[ja]
        d = np.linalg.norm(p1 - p2)
        bond_lengths.append(d)
    bond_lengths = np.array(bond_lengths)
    vmin, vmax = bond_lengths.min(), bond_lengths.max()
    
    bond_segs = []
    bond_colors = []
    for i, b in enumerate(bonds):
        ia, ja = b
        p1, p2 = pos[ia], pos[ja]
        bond_segs.append(p1)
        bond_segs.append(p2)
        
        if abs(vmax - vmin) < 1e-4:
            f = 0.5
        else:
            f = (bond_lengths[i] - vmin) / (vmax - vmin)
        color = (f, 0.0, 1.0 - f, 0.8)
        bond_colors.append(color)
        bond_colors.append(color)
    
    return np.array(bond_segs, dtype=np.float32), np.array(bond_colors, dtype=np.float32)


def generate_atom_labels(label_mode, pos, enames, atom_subtype=None, backend=None, bonds=None):
    """Generate text labels for atoms based on label_mode.
    
    Args:
        label_mode: String specifying label type
        pos: (n,3) array of positions
        enames: List/array of element names
        atom_subtype: Optional list of atom subtypes
        backend: Optional backend object for subtype queries
        bonds: Optional list of bonds for bond length labels
    
    Returns:
        lbl_pos: List of label positions
        lbl_texts: List of label text strings
    """
    lbl_pos = []
    lbl_texts = []
    
    if label_mode == 'Element+Index':
        for i, e in enumerate(enames):
            if e != 'H':
                lbl_pos.append(pos[i])
                lbl_texts.append(f"{e}{i}")
    elif label_mode == 'Atomic Type':
        for i, subtype in enumerate(atom_subtype or []):
            if enames[i] != 'H':
                lbl_pos.append(pos[i])
                if 'sp3' in subtype:
                    lbl_texts.append('sp3')
                elif 'sp2' in subtype:
                    lbl_texts.append('sp2')
                elif 'sp' in subtype:
                    lbl_texts.append('sp')
                else:
                    lbl_texts.append(subtype)
    elif label_mode == 'Pi Orbitals':
        for i in range(len(enames)):
            if i < len(atom_subtype or []):
                subtype = atom_subtype[i]
                if enames[i] != 'H':
                    lbl_pos.append(pos[i])
                    npi = backend._get_npi_from_subtype(subtype) if backend else 0
                    lbl_texts.append(str(npi))
    elif label_mode == 'Z-Height':
        for i, e in enumerate(enames):
            if e != 'H':
                lbl_pos.append(pos[i])
                lbl_texts.append(f"{pos[i, 2]:.2f}")
    elif label_mode == 'Charge':
        for i, e in enumerate(enames):
            if e != 'H':
                lbl_pos.append(pos[i])
                lbl_texts.append("0")
    elif label_mode == 'Bond Lengths':
        if bonds:
            for b in bonds:
                ia, ja = b
                p1, p2 = pos[ia], pos[ja]
                d = np.linalg.norm(p1 - p2)
                lbl_pos.append((p1 + p2) * 0.5)
                lbl_texts.append(f"{d:.3f}")
    
    return lbl_pos, lbl_texts


class AtomScene(QtCore.QObject):
    """Reusable Vispy widget for atoms (+ optional bonds) with orthographic top-down view.

    This is intended as a generic MD/FF viewer: operates on arrays (positions, colors, sizes, bonds)
    and optional per-atom vectors (forces) for visualization.

    Picking/dragging is implemented in a pseudo-2D mode:
    - camera is fixed to top-down view
    - drag moves atoms in XY plane (Z unchanged unless you set it externally)
    """

    sig_atom_picked = QtCore.pyqtSignal(int)
    sig_drag_state = QtCore.pyqtSignal(int, int, object)  # active(0/1), idx, pos3
    sig_atom_moved = QtCore.pyqtSignal(int, object)  # idx, pos3 - emitted during drag
    sig_rmb_remove = QtCore.pyqtSignal(int)  # idx to remove
    sig_selection_changed = QtCore.pyqtSignal(object)  # set of selected indices
    sig_camera_changed = QtCore.pyqtSignal()  # camera changed (zoom/pan/rotate)

    def __init__(self, *, bgcolor='white', backend=None):
        super().__init__(parent=None)
        self.backend = backend  # Reference to backend for authoritative geometry

        self.canvas = scene.SceneCanvas(keys='interactive', bgcolor=bgcolor, show=False)
        self.view = self.canvas.central_widget.add_view()
        # Ortho top-down like MolecularPlacerVisPy
        self.view.camera = scene.TurntableCamera(fov=0, distance=80, elevation=90, azimuth=0)
        self.view.camera.interactive = False  # disable default camera mouse handling; we handle RMB explicitly

        self._cam_debug = 0
        self._rmb_down = False
        self._rmb_last = None
        self._cam_rot_speed = 0.3  # deg per pixel
        self._cam_zoom_speed = 0.12
        self._cam_zoom_min = 1e-4
        self._cam_zoom_max = 1e+4
        self._cam_pan_speed = 2.0  # world units per key press

        # Draw ordering: radius behind everything, then bboxes/links/lines, then atom centers, then labels.
        self.radius_markers = visuals.Markers(parent=self.view.scene)
        self.bbox_lines = visuals.Line(parent=self.view.scene, color=(0.8, 0.0, 0.0, 0.9), width=2.0, antialias=True, method='gl')
        self.inbox_lines = visuals.Line(parent=self.view.scene, color=(0.0, 0.0, 0.0, 0.35), width=1.0, antialias=True, method='gl')
        self.halo_lines  = visuals.Line(parent=self.view.scene, color=(0.8, 0.1, 0.8, 0.35), width=1.0, antialias=True, method='gl')
        self.neigh_lines = visuals.Line(parent=self.view.scene, color=(0.2, 0.2, 0.2, 0.65), width=1.0, antialias=True, method='gl')
        self.port_lines = visuals.Line(parent=self.view.scene, color=(1.0, 0.55, 0.0, 0.55), width=1.0, antialias=True, method='gl')
        self.port_target_lines = visuals.Line(parent=self.view.scene, color=(0.0, 0.7, 0.9, 0.55), width=1.0, antialias=True, method='gl')
        self.dpos_lines = visuals.Line(parent=self.view.scene, color=(0.9, 0.0, 0.0, 0.75), width=1.6, antialias=True, method='gl')
        self.dpos_neigh_lines = visuals.Line(parent=self.view.scene, color=(0.2, 0.2, 1.0, 0.75), width=1.6, antialias=True, method='gl')
        self.bond_lines = visuals.Line(parent=self.view.scene, color='gray', width=1.5, antialias=True, method='gl')
        self.bond_colored_lines = visuals.Line(parent=self.view.scene, color='gray', width=3.0, antialias=True, method='gl')
        self.ch_bond_lines = visuals.Line(parent=self.view.scene, color=(0.4, 0.4, 0.4, 0.6), width=1.0, antialias=True, method='gl')
        self.hbond_lines = visuals.Line(parent=self.view.scene, color=(0.8, 0.2, 0.8, 0.5), width=1.5, antialias=True, method='gl')
        self.force_lines = visuals.Line(parent=self.view.scene, color=(1, 0, 0, 0.8), width=2.0, antialias=True, method='gl')
        self.atom_markers = visuals.Markers(parent=self.view.scene)
        self.axes = visuals.XYZAxis(parent=self.view.scene)
        self.text_labels = visuals.Text(parent=self.view.scene, color='black', font_size=10, anchor_x='left', anchor_y='bottom')
        # Hover visuals for debug visualization
        self.hover_bond_line = visuals.Line(parent=self.view.scene, color='lime', width=4.0, antialias=True, method='gl')
        self.hover_ring_lines = visuals.Line(parent=self.view.scene, color='cyan', width=2.0, antialias=True, method='gl')
        self.hover_ring_markers = visuals.Markers(parent=self.view.scene)
        self.hover_ring_text = visuals.Text(parent=self.view.scene, color='cyan', font_size=12, anchor_x='center', anchor_y='center')
        self.hover_atom_marker = visuals.Markers(parent=self.view.scene)
        # Selection rectangle (Line visual) - create lazily to avoid initialization issues
        self.selection_rect = None

        # Enforce z-order when supported
        for o, v in enumerate((self.radius_markers, self.bbox_lines, self.inbox_lines, self.halo_lines, self.neigh_lines, self.port_lines, self.port_target_lines, self.dpos_lines, self.dpos_neigh_lines, self.bond_lines, self.bond_colored_lines, self.ch_bond_lines, self.hbond_lines, self.force_lines, self.atom_markers, self.axes, self.text_labels, self.hover_bond_line, self.hover_ring_lines, self.hover_ring_markers, self.hover_ring_text, self.hover_atom_marker)):
            if hasattr(v, 'order'):
                v.order = int(o)

        # GL state: radius translucent and never blocks other overlays
        try:
            self.radius_markers.set_gl_state('translucent', depth_test=False)
            for v in (self.bbox_lines, self.inbox_lines, self.halo_lines, self.neigh_lines, self.port_lines, self.port_target_lines, self.dpos_lines, self.dpos_neigh_lines, self.bond_lines, self.bond_colored_lines, self.ch_bond_lines, self.hbond_lines, self.force_lines, self.hover_bond_line, self.hover_ring_lines):
                v.set_gl_state('translucent', depth_test=False)
            self.atom_markers.set_gl_state('translucent', depth_test=False)
            self.hover_ring_markers.set_gl_state('translucent', depth_test=False)
            self.hover_atom_marker.set_gl_state('translucent', depth_test=False)
            self.text_labels.set_gl_state('translucent', depth_test=False)
            self.hover_ring_text.set_gl_state('translucent', depth_test=False)
        except Exception:
            pass

        self._pos = np.zeros((0, 3), dtype=np.float32)
        self._colors = None
        self._sizes = None
        self._bonds = None
        self._forces = None
        self._radius = None

        self._render_mask = None
        self._group_size = 64
        self._color_by_group = False
        self._colors_base = None

        self._show_radius = False
        self._show_bboxes = False
        self._show_inbox_links = False
        self._show_halo_links = False
        self._show_neigh_bonds = False
        self._show_port_tips = False
        self._show_port_targets = False
        self._show_dpos = False
        self._show_dpos_neigh = False
        self._show_axes = True
        self._bboxes_min = None
        self._bboxes_max = None
        self._inbox_link_segs = None
        self._halo_link_segs = None
        self._inbox_link_gid = None
        self._halo_link_gid = None

        self._neigh_segs = None
        self._port_tip_segs = None
        self._port_target_segs = None
        self._dpos_segs = None
        self._dpos_neigh_segs = None

        self._label_mode = 'none'  # none|global|local|pair|radius
        self._labels_text = None

        self._marker_style = 'disc'      # vispy marker name
        self._radius_style = 'disc'

        self._pick_active = False
        self._pick_idx = -1
        self._pick_z = 0.0
        self.lock_drag = False       # Set True externally to suppress all atom drag

        self._pick_mode = '2d'   # '2d' or '3d'
        self._lock_top_view = True
        self._clamp_xy = False
        self._fixed = None

        # Selection state
        self._selection_mode = False  # If True, RMB drags to select atoms instead of camera rotation
        self._selected_indices = set()
        self._selection_start = None
        self._selection_end = None
        self._selected_colors_backup = None  # Store original colors for selected atoms

        self._drag_plane_p0 = None
        self._drag_plane_n = None

        self.canvas.events.mouse_press.connect(self._on_mouse_press)
        self.canvas.events.mouse_release.connect(self._on_mouse_release)
        self.canvas.events.mouse_move.connect(self._on_mouse_move)
        self.canvas.events.mouse_wheel.connect(self._on_mouse_wheel)
        self.canvas.events.key_press.connect(self._on_key_press)

        self._apply_camera_mode()

    @property
    def widget(self):
        return self.canvas.native

    def set_data(self, pos, *, colors=None, sizes=None, bonds=None, forces=None, force_scale=1.0):
        pos = _as_f32(pos)
        if pos.ndim != 2 or pos.shape[1] != 3:
            raise ValueError(f"AtomScene.set_data: pos.shape={pos.shape} expected (n,3)")
        # If backend is provided, use its authoritative positions as reference (not copy)
        # But keep a copy for rendering to avoid issues during camera changes
        if self.backend is not None:
            self._pos = self.backend.sys.apos.astype(np.float32).copy()
        else:
            self._pos = pos
        self._render_mask = None
        if (self._fixed is None) or (self._fixed.shape[0] != self._pos.shape[0]):
            self._fixed = np.zeros((self._pos.shape[0],), dtype=bool)
        self._colors = None if colors is None else _as_f32(colors)
        self._colors_base = None if colors is None else _as_f32(colors)
        self._sizes = None if sizes is None else _as_f32(sizes)
        self._bonds = None if bonds is None else np.asarray(bonds, dtype=np.int32)
        self._forces = None if forces is None else _as_f32(forces)
        self._force_scale = float(force_scale)
        self._redraw()

    def set_marker_style(self, style='disc'):
        style = str(style)
        self._marker_style = style
        self._redraw()

    def set_radius_style(self, style='disc'):
        style = str(style)
        self._radius_style = style
        self._redraw()

    def set_radius(self, radius):
        r = np.asarray(radius, dtype=np.float32)
        if r.shape != (self._pos.shape[0],):
            raise ValueError(f"AtomScene.set_radius: radius.shape={r.shape} expected ({self._pos.shape[0]},)")
        self._radius = r
        self._redraw()

    def set_render_mask(self, mask):
        mask = np.asarray(mask, dtype=bool)
        if mask.shape != (self._pos.shape[0],):
            raise ValueError(f"AtomScene.set_render_mask: mask.shape={mask.shape} expected ({self._pos.shape[0]},)")
        self._render_mask = mask.copy()
        self._redraw()

    def set_group_size(self, group_size):
        self._group_size = int(group_size)

    def set_color_by_group(self, enable):
        self._color_by_group = bool(enable)
        self._redraw()

    def set_show_radius(self, enable):
        self._show_radius = bool(enable)
        self._redraw()

    def set_show_bboxes(self, enable):
        self._show_bboxes = bool(enable)
        self._redraw()

    def set_show_inbox_links(self, show):
        self._show_inbox_links = bool(show)
        self._redraw()

    def set_show_halo_links(self, show):
        self._show_halo_links = bool(show)
        self._redraw()

    def set_show_neigh_bonds(self, show):
        self._show_neigh_bonds = bool(show)
        self._redraw()

    def set_show_port_tips(self, show):
        self._show_port_tips = bool(show)
        self._redraw()

    def set_show_port_targets(self, show):
        self._show_port_targets = bool(show)
        self._redraw()

    def set_show_dpos(self, show):
        self._show_dpos = bool(show)
        self._redraw()

    def set_show_dpos_neigh(self, show):
        self._show_dpos_neigh = bool(show)
        self._redraw()

    def set_show_axes(self, enable):
        self._show_axes = bool(enable)
        self.axes.visible = bool(enable)
        self.canvas.update()

    def set_inbox_links(self, segs, *, gids=None):
        if segs is None:
            self._inbox_link_segs = None
            self._inbox_link_gid = None
        else:
            s = _as_f32(segs)
            if s.ndim != 2 or s.shape[1] != 3 or (s.shape[0] % 2) != 0:
                raise ValueError(f"AtomScene.set_inbox_links: segs.shape={s.shape} expected (2*m,3)")
            self._inbox_link_segs = s
            if gids is None:
                self._inbox_link_gid = None
            else:
                g = np.asarray(gids, dtype=np.int32)
                if g.shape != (s.shape[0] // 2,):
                    raise ValueError(f"AtomScene.set_inbox_links: gids.shape={g.shape} expected ({s.shape[0]//2},)")
                self._inbox_link_gid = g
        self._redraw()

    def set_halo_links(self, segs, *, gids=None):
        if segs is None:
            self._halo_link_segs = None
            self._halo_link_gid = None
        else:
            s = _as_f32(segs)
            if s.ndim != 2 or s.shape[1] != 3 or (s.shape[0] % 2) != 0:
                raise ValueError(f"AtomScene.set_halo_links: segs.shape={s.shape} expected (2*m,3)")
            self._halo_link_segs = s
            if gids is None:
                self._halo_link_gid = None
            else:
                g = np.asarray(gids, dtype=np.int32)
                if g.shape != (s.shape[0] // 2,):
                    raise ValueError(f"AtomScene.set_halo_links: gids.shape={g.shape} expected ({s.shape[0]//2},)")
                self._halo_link_gid = g
        self._redraw()

    def set_neigh_bonds(self, segs):
        if segs is None:
            self._neigh_segs = None
        else:
            s = _as_f32(segs)
            if s.ndim != 2 or s.shape[1] != 3 or (s.shape[0] % 2) != 0:
                raise ValueError(f"AtomScene.set_neigh_bonds: segs.shape={s.shape} expected (2*m,3)")
            self._neigh_segs = s
        self._redraw()

    def set_port_tips(self, segs):
        if segs is None:
            self._port_tip_segs = None
        else:
            s = _as_f32(segs)
            if s.ndim != 2 or s.shape[1] != 3 or (s.shape[0] % 2) != 0:
                raise ValueError(f"AtomScene.set_port_tips: segs.shape={s.shape} expected (2*m,3)")
            self._port_tip_segs = s
        self._redraw()

    def set_port_targets(self, segs):
        if segs is None:
            self._port_target_segs = None
        else:
            s = _as_f32(segs)
            if s.ndim != 2 or s.shape[1] != 3 or (s.shape[0] % 2) != 0:
                raise ValueError(f"AtomScene.set_port_targets: segs.shape={s.shape} expected (2*m,3)")
            self._port_target_segs = s
        self._redraw()

    def set_dpos(self, segs):
        if segs is None:
            self._dpos_segs = None
        else:
            s = _as_f32(segs)
            if s.ndim != 2 or s.shape[1] != 3 or (s.shape[0] % 2) != 0:
                raise ValueError(f"AtomScene.set_dpos: segs.shape={s.shape} expected (2*m,3)")
            self._dpos_segs = s
        self._redraw()

    def set_dpos_neigh(self, segs):
        if segs is None:
            self._dpos_neigh_segs = None
        else:
            s = _as_f32(segs)
            if s.ndim != 2 or s.shape[1] != 3 or (s.shape[0] % 2) != 0:
                raise ValueError(f"AtomScene.set_dpos_neigh: segs.shape={s.shape} expected (2*m,3)")
            self._dpos_neigh_segs = s
        self._redraw()

    def set_bboxes(self, bmin, bmax):
        bmin = _as_f32(bmin)
        bmax = _as_f32(bmax)
        if bmin.shape != bmax.shape or bmin.ndim != 2 or bmin.shape[1] != 4:
            raise ValueError(f"AtomScene.set_bboxes: bmin.shape={bmin.shape} bmax.shape={bmax.shape} expected (ng,4)")
        self._bboxes_min = bmin
        self._bboxes_max = bmax
        self._redraw()

    def set_label_mode(self, mode):
        mode = str(mode).lower()
        if mode not in ('none', 'global', 'local', 'pair', 'radius'):
            raise ValueError(f"AtomScene.set_label_mode: mode={mode} expected none|global|local|pair|radius")
        self._label_mode = mode
        self._labels_text = None
        self._redraw()

    def _px_per_world_ortho(self):
        """Pixels per 1 world unit for TurntableCamera with fov=0.

        To keep glyph size independent of camera orientation, use only zoom (scale_factor) and viewport.
        """
        cam = self.view.camera
        if cam is None:
            return 1.0
        sf = float(getattr(cam, 'scale_factor', 1.0))
        if (not np.isfinite(sf)) or (sf <= 1e-12):
            return 1.0
        tr = self.view.scene.transform
        p0 = np.array(tr.imap((0.0, 0.0, 0.0)), dtype=np.float32)
        p1 = np.array(tr.imap((1.0, 0.0, 0.0)), dtype=np.float32)
        world_len = float(np.linalg.norm(p1[:2] - p0[:2]))
        if (not np.isfinite(world_len)) or (world_len <= 1e-12):
            return 1.0
        return 1.0 / world_len

    def get_zoom(self):
        cam = self.view.camera
        if cam is None:
            return 1.0
        sf = getattr(cam, 'scale_factor', 1.0)
        return float(sf)

    def set_zoom(self, zoom):
        cam = self.view.camera
        if cam is None:
            return
        z = float(zoom)
        if z < self._cam_zoom_min:
            z = self._cam_zoom_min
        if z > self._cam_zoom_max:
            z = self._cam_zoom_max
        cam.scale_factor = z
        self._redraw()
        self.sig_camera_changed.emit()

    def reset_view(self):
        cam = self.view.camera
        if cam is None:
            return
        cam.fov = 0
        cam.azimuth = 0
        cam.elevation = 90
        cam.roll = 0
        cam.scale_factor = 1.0
        self.canvas.update()

    def set_pick_mode(self, mode):
        mode = str(mode).lower()
        if mode not in ('2d', '3d'):
            raise ValueError(f"AtomScene.set_pick_mode: mode={mode} expected '2d'|'3d'")
        self._pick_mode = mode

    def set_lock_top_view(self, lock):
        self._lock_top_view = bool(lock)
        self._apply_camera_mode()

    def set_camera_debug(self, level=1):
        self._cam_debug = int(level)

    def set_clamp_xy(self, clamp):
        self._clamp_xy = bool(clamp)

    def set_fixed_mask(self, fixed):
        fixed = np.asarray(fixed, dtype=bool)
        if fixed.shape != (self._pos.shape[0],):
            raise ValueError(f"AtomScene.set_fixed_mask: fixed.shape={fixed.shape} expected ({self._pos.shape[0]},)")
        self._fixed = fixed.copy()

    def toggle_fixed(self, i):
        i = int(i)
        if i < 0 or i >= self._pos.shape[0]:
            raise ValueError(f"AtomScene.toggle_fixed: i={i} out of range")
        self._fixed[i] = ~self._fixed[i]
        return bool(self._fixed[i])

    def is_fixed(self, i):
        i = int(i)
        if self._fixed is None:
            return False
        return bool(self._fixed[i])

    def update_positions(self, pos):
        pos = _as_f32(pos)
        if pos.shape != self._pos.shape:
            raise ValueError(f"AtomScene.update_positions: pos.shape={pos.shape} != current {self._pos.shape}")
        self._pos = pos
        self._redraw()

    def _apply_camera_mode(self):
        cam = self.view.camera
        if cam is None:
            return
        cam.interactive = False
        if self._lock_top_view:
            cam.fov = 0
            cam.elevation = 90
            cam.azimuth = 0
            cam.roll = 0
        # leave distance/center as is

    def _cam_print(self, tag):
        if int(self._cam_debug) <= 0:
            return
        cam = self.view.camera
        if cam is None:
            return
        print(f"[CAM] {tag} az={float(cam.azimuth):.3f} el={float(cam.elevation):.3f} dist={float(cam.distance):.3f}")

    def _cam_rotate(self, dx_px, dy_px):
        if self._lock_top_view:
            return
        cam = self.view.camera
        if cam is None:
            return
        cam.azimuth = float(cam.azimuth) + float(dx_px) * float(self._cam_rot_speed)
        cam.elevation = float(cam.elevation) + float(dy_px) * float(self._cam_rot_speed)
        if cam.elevation > 89.0:
            cam.elevation = 89.0
        if cam.elevation < -89.0:
            cam.elevation = -89.0
        self._redraw()
        self.sig_camera_changed.emit()
        self._cam_print('rotate')

    def _cam_zoom(self, delta):
        cam = self.view.camera
        if cam is None:
            return
        # Orthographic zoom: change camera scale_factor (distance does not change zoom for fov=0).
        z0 = float(getattr(cam, 'scale_factor', 1.0))
        s = float(np.exp(-float(delta) * float(self._cam_zoom_speed)))
        z1 = z0 * s
        if z1 < self._cam_zoom_min:
            z1 = self._cam_zoom_min
        if z1 > self._cam_zoom_max:
            z1 = self._cam_zoom_max
        cam.scale_factor = z1
        self._redraw()
        if int(self._cam_debug) > 0:
            print(f"[CAM] zoom delta={float(delta):.6g} scale:{z0:.6g}->{z1:.6g}")

    def _cam_pan(self, dx, dy):
        cam = self.view.camera
        if cam is None:
            return
        # Pan camera center by dx, dy in world units
        center = np.array(cam.center)
        center[0] += float(dx) * float(self._cam_pan_speed)
        center[1] += float(dy) * float(self._cam_pan_speed)
        cam.center = tuple(center)
        self._redraw()
        self.sig_camera_changed.emit()
        if int(self._cam_debug) > 0:
            print(f"[CAM] pan dx={float(dx):.3f} dy={float(dy):.3f} center={tuple(center)}")

    def _ray_from_mouse(self, mouse_pos, z0=0.0, z1=1.0):
        # mouse_pos in canvas pixels
        # If the view is shifted (e.g. in a Grid), we must use local coordinates
        view_pos = np.array(mouse_pos) - self.view.pos[:2]
        tr = self.view.scene.transform
        p0 = np.array(tr.imap((view_pos[0], view_pos[1], float(z0)))[:3], dtype=np.float32)
        p1 = np.array(tr.imap((view_pos[0], view_pos[1], float(z1)))[:3], dtype=np.float32)
        d = p1 - p0
        dn = float(np.linalg.norm(d))
        if dn <= 1e-20:
            d = np.array([0.0, 0.0, 1.0], dtype=np.float32)
        else:
            d /= dn
        return p0, d

    def _pick_idx_from_ray(self, r0, rd):
        # closest point distance^2 to ray for each atom
        # d2 = |(p-r0) - rd*dot(rd,(p-r0))|^2
        valid = np.isfinite(self._pos).all(axis=1)
        if self._fixed is not None:
            valid &= (~self._fixed)
        if not np.any(valid):
            return -1, 1e30
        dp = self._pos - r0[None, :]
        t = np.sum(dp * rd[None, :], axis=1)
        q = dp - rd[None, :] * t[:, None]
        d2 = np.sum(q * q, axis=1)
        d2 = np.where(valid, d2, 1e30)
        i = int(np.argmin(d2))
        return i, float(d2[i])

    def _validated_segs(self, tag, segs):
        s = np.asarray(segs, dtype=np.float32)
        if s.ndim != 2 or s.shape[1] != 3 or (s.shape[0] % 2) != 0:
            raise ValueError(f"{tag}: segs.shape={s.shape} expected (2*m,3)")
        if not np.isfinite(s).all():
            bad = np.where(~np.isfinite(s))[0][:10]
            raise ValueError(f"{tag}: non-finite entries at idx={bad.tolist()}")
        return s

    def _line_set(self, tag, visual, segs, *, color=None, width=None, connect=None):
        if segs is None:
            visual.set_data(np.zeros((0, 3), dtype=np.float32))
            return
        try:
            s = self._validated_segs(tag, segs)
        except Exception as e:
            print(f"[VISPY-SEG-ERR] {tag} validation failed: {e}")
            visual.set_data(np.zeros((0, 3), dtype=np.float32))
            return
        if s.size == 0:
            visual.set_data(np.zeros((0, 3), dtype=np.float32))
            return
        conn = connect
        if conn is None:
            conn = np.zeros((s.shape[0],), dtype=bool); conn[0::2] = True
        # validate color length if array-like
        if hasattr(color, 'shape') and hasattr(color, '__len__'):
            clen = int(len(color))
            if clen not in (0, s.shape[0]):
                print(f"[VISPY-SEG-ERR] {tag} color length mismatch: len(color)={clen} verts={s.shape[0]}")
                color = None
        try:
            visual.set_data(s, connect=conn, color=color, width=width)
        except Exception as e:
            stats = {
                'tag': tag,
                'segs_shape': s.shape,
                'segs_min': float(np.min(s)) if s.size else None,
                'segs_max': float(np.max(s)) if s.size else None,
                'connect_shape': getattr(conn, 'shape', None),
                'color_shape': getattr(color, 'shape', None) if hasattr(color, 'shape') else None,
                'width': width,
            }
            print(f"[VISPY-LINE-ERR] {tag} set_data failed: {e} | stats={stats}")
            visual.set_data(np.zeros((0, 3), dtype=np.float32))

    def _intersect_ray_plane(self, r0, rd, p0, n):
        denom = float(np.dot(rd, n))
        if abs(denom) < 1e-12:
            return None
        t = float(np.dot((p0 - r0), n) / denom)
        return r0 + rd * t

    def _redraw(self):
        # Don't sync _pos from backend here - causes stale bond index issues during camera changes
        # The GUI should call refresh_view() when atoms actually change
        # Camera changes should only re-render, not re-sync data

        if self._pos.size == 0:
            idx = np.array([], dtype=int)
        else:
            m = self._render_mask
            if m is None:
                idx = np.arange(self._pos.shape[0], dtype=int)
            else:
                idx = np.where(m)[0].astype(int)

        if idx.size == 0:
            self.atom_markers.set_data(np.zeros((0, 3), dtype=np.float32))
            self.radius_markers.set_data(np.zeros((0, 3), dtype=np.float32))
            self.bbox_lines.set_data(np.zeros((0, 3), dtype=np.float32))
            self.inbox_lines.set_data(np.zeros((0, 3), dtype=np.float32))
            self.halo_lines.set_data(np.zeros((0, 3), dtype=np.float32))
            self.neigh_lines.set_data(np.zeros((0, 3), dtype=np.float32))
            self.port_lines.set_data(np.zeros((0, 3), dtype=np.float32))
            self.port_target_lines.set_data(np.zeros((0, 3), dtype=np.float32))
            self.dpos_lines.set_data(np.zeros((0, 3), dtype=np.float32))
            self.dpos_neigh_lines.set_data(np.zeros((0, 3), dtype=np.float32))
            self.bond_lines.set_data(np.zeros((0, 3), dtype=np.float32))
            self.bond_colored_lines.set_data(np.zeros((0, 3), dtype=np.float32))
            self.ch_bond_lines.set_data(np.zeros((0, 3), dtype=np.float32))
            self.hbond_lines.set_data(np.zeros((0, 3), dtype=np.float32))
            self.force_lines.set_data(np.zeros((0, 3), dtype=np.float32))
            return

        # Colors
        if self._colors is None:
            # Default coloring
            face_color = np.zeros((idx.size, 4), dtype=np.float32)
            face_color[:, 0] = 0.5
            face_color[:, 1] = 0.5
            face_color[:, 2] = 0.5
            face_color[:, 3] = 1.0

        if self._color_by_group:
            # deterministic HSV-like palette per group
            g = (idx // int(self._group_size)).astype(np.int32)
            c = np.empty((idx.size, 4), dtype=np.float32)
            for i, gi in enumerate(g):
                h = float((gi * 0.61803398875) % 1.0)
                r = abs(h * 6.0 - 3.0) - 1.0
                g1 = 2.0 - abs(h * 6.0 - 2.0)
                b = 2.0 - abs(h * 6.0 - 4.0)
                rgb = np.clip(np.array([r, g1, b], dtype=np.float32), 0.0, 1.0)
                c[i, :3] = rgb
                c[i, 3] = 0.9
            face_color = c
        else:
            if self._colors is None:
                face_color = (0.2, 0.2, 0.2, 1.0)
            else:
                face_color = self._colors[idx]
        if self._sizes is None:
            size = 8.0
            size = np.full((idx.size,), float(size), dtype=np.float32)
        else:
            size = _as_f32(self._sizes[idx])

        # marker style (disc/square, etc.)
        try:
            self.atom_markers.set_data(self._pos[idx], face_color=face_color, size=size, edge_width=0.5, edge_color='black', symbol=self._marker_style)
        except TypeError:
            # older vispy uses 'marker' kw
            self.atom_markers.set_data(self._pos[idx], face_color=face_color, size=size, edge_width=0.5, edge_color='black', marker=self._marker_style)

        if self._show_radius and (self._radius is not None):
            # world radius -> exact screen size scaling for orthographic camera (depends only on zoom)
            r = np.maximum(self._radius[idx], 0.0)
            px_per_world = float(self._px_per_world_ortho())
            sizeR = (2.0 * r * px_per_world).astype(np.float32)
            colR = np.zeros((idx.size, 4), dtype=np.float32)
            colR[:, :3] = face_color[:, :3] if isinstance(face_color, np.ndarray) else np.array(face_color[:3], dtype=np.float32)[None, :]
            colR[:, 3] = 0.10
            try:
                self.radius_markers.set_data(self._pos[idx], face_color=colR, size=sizeR, edge_width=0.0, symbol=self._radius_style)
            except TypeError:
                self.radius_markers.set_data(self._pos[idx], face_color=colR, size=sizeR, edge_width=0.0, marker=self._radius_style)
        else:
            self.radius_markers.set_data(np.zeros((0, 3), dtype=np.float32))

        # Debug link lines
        if self._show_inbox_links and (self._inbox_link_segs is not None) and (self._inbox_link_segs.size > 0):
            if self._color_by_group and (self._inbox_link_gid is not None):
                gid = self._inbox_link_gid
                col = np.empty((gid.size * 2, 4), dtype=np.float32)
                for i, gi in enumerate(gid):
                    h = float((int(gi) * 0.61803398875) % 1.0)
                    r = abs(h * 6.0 - 3.0) - 1.0
                    g1 = 2.0 - abs(h * 6.0 - 2.0)
                    b = 2.0 - abs(h * 6.0 - 4.0)
                    rgb = np.clip(np.array([r, g1, b], dtype=np.float32), 0.0, 1.0)
                    col[2*i+0, :3] = rgb; col[2*i+1, :3] = rgb
                    col[2*i+0, 3] = 0.35; col[2*i+1, 3] = 0.35
                self._line_set("inbox_links", self.inbox_lines, self._inbox_link_segs, color=col)
            else:
                self._line_set("inbox_links", self.inbox_lines, self._inbox_link_segs)
        else:
            self.inbox_lines.set_data(np.zeros((0, 3), dtype=np.float32))

        if self._show_halo_links and (self._halo_link_segs is not None) and (self._halo_link_segs.size > 0):
            if self._color_by_group and (self._halo_link_gid is not None):
                gid = self._halo_link_gid
                col = np.empty((gid.size * 2, 4), dtype=np.float32)
                for i, gi in enumerate(gid):
                    h = float((int(gi) * 0.61803398875) % 1.0)
                    r = abs(h * 6.0 - 3.0) - 1.0
                    g1 = 2.0 - abs(h * 6.0 - 2.0)
                    b = 2.0 - abs(h * 6.0 - 4.0)
                    rgb = np.clip(np.array([r, g1, b], dtype=np.float32), 0.0, 1.0)
                    col[2*i+0, :3] = rgb; col[2*i+1, :3] = rgb
                    col[2*i+0, 3] = 0.35; col[2*i+1, 3] = 0.35
                self._line_set("halo_links", self.halo_lines, self._halo_link_segs, color=col)
            else:
                self._line_set("halo_links", self.halo_lines, self._halo_link_segs)
        else:
            self.halo_lines.set_data(np.zeros((0, 3), dtype=np.float32))

        if self._show_neigh_bonds and (self._neigh_segs is not None) and (self._neigh_segs.size > 0):
            self._line_set("neigh_bonds", self.neigh_lines, self._neigh_segs)
        else:
            self.neigh_lines.set_data(np.zeros((0, 3), dtype=np.float32))

        if self._show_port_tips and (self._port_tip_segs is not None) and (self._port_tip_segs.size > 0):
            self._line_set("port_tips", self.port_lines, self._port_tip_segs)
        else:
            self.port_lines.set_data(np.zeros((0, 3), dtype=np.float32))

        if self._show_port_targets and (self._port_target_segs is not None) and (self._port_target_segs.size > 0):
            self._line_set("port_targets", self.port_target_lines, self._port_target_segs)
        else:
            self.port_target_lines.set_data(np.zeros((0, 3), dtype=np.float32))

        if self._show_dpos and (self._dpos_segs is not None) and (self._dpos_segs.size > 0):
            self._line_set("dpos", self.dpos_lines, self._dpos_segs)
        else:
            self.dpos_lines.set_data(np.zeros((0, 3), dtype=np.float32))

        if self._show_dpos_neigh and (self._dpos_neigh_segs is not None) and (self._dpos_neigh_segs.size > 0):
            self._line_set("dpos_neigh", self.dpos_neigh_lines, self._dpos_neigh_segs)
        else:
            self.dpos_neigh_lines.set_data(np.zeros((0, 3), dtype=np.float32))

        # Bonds: draw segment pairs (normal bonds - GUI handles CH/H-bonds separately)
        if (self._bonds is not None) and (self._bonds.size > 0):
            b = self._bonds
            if m is not None:
                mb = m[b[:, 0]] & m[b[:, 1]]
                b = b[mb]
            if b.size > 0:
                segs = np.empty((b.shape[0] * 2, 3), dtype=np.float32)
                segs[0::2] = self._pos[b[:, 0]]
                segs[1::2] = self._pos[b[:, 1]]
                self._line_set("bonds", self.bond_lines, segs, color=(0.3, 0.3, 0.3, 0.8), width=1.5)
            else:
                self.bond_lines.set_data(np.zeros((0, 3), dtype=np.float32))
        else:
            self.bond_lines.set_data(np.zeros((0, 3), dtype=np.float32))

        # Forces: per-atom line from pos to pos+f*scale
        if self._forces is not None:
            f = self._forces
            if f.shape != self._pos.shape:
                raise ValueError(f"AtomScene._redraw: forces.shape={f.shape} expected {self._pos.shape}")
            segs = np.empty((idx.size * 2, 3), dtype=np.float32)
            segs[0::2] = self._pos[idx]
            segs[1::2] = self._pos[idx] + f[idx] * self._force_scale
            self._line_set("forces", self.force_lines, segs)
        else:
            self.force_lines.set_data(np.zeros((0, 3), dtype=np.float32))

        # Bounding boxes (clusters)
        if self._show_bboxes and (self._bboxes_min is not None) and (self._bboxes_max is not None):
            bmin = self._bboxes_min
            bmax = self._bboxes_max
            ng = int(bmin.shape[0])
            segs = []
            connect = []
            for ig in range(ng):
                mn = bmin[ig, :3]; mx = bmax[ig, :3]
                v = np.array([
                    [mn[0], mn[1], mn[2]], [mx[0], mn[1], mn[2]],
                    [mx[0], mx[1], mn[2]], [mn[0], mx[1], mn[2]],
                    [mn[0], mn[1], mx[2]], [mx[0], mn[1], mx[2]],
                    [mx[0], mx[1], mx[2]], [mn[0], mx[1], mx[2]],
                ], dtype=np.float32)
                e = [(0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),(0,4),(1,5),(2,6),(3,7)]
                for (a,b) in e:
                    segs.append(v[a]); segs.append(v[b])
                    connect.append(True); connect.append(False)
            if len(segs) > 0:
                segs = np.asarray(segs, dtype=np.float32)
                connect = np.asarray(connect, dtype=bool)
                self.bbox_lines.set_data(segs, connect=connect)
            else:
                self.bbox_lines.set_data(np.zeros((0, 3), dtype=np.float32))
        else:
            self.bbox_lines.set_data(np.zeros((0, 3), dtype=np.float32))

        # Labels - only manage internally if no backend (GUI manages labels externally)
        if self.backend is None:
            if self._label_mode == 'none':
                self.text_labels.text = ['']
                self.text_labels.pos = np.zeros((1, 3), dtype=np.float32)
                self.text_labels.visible = False
            else:
                if self._labels_text is None:
                    txt = []
                    for ii in idx:
                        if self._label_mode == 'global':
                            txt.append(str(int(ii)))
                        elif self._label_mode == 'local':
                            txt.append(str(int(ii % int(self._group_size))))
                        elif self._label_mode == 'radius':
                            if self._radius is None:
                                txt.append('nan')
                            else:
                                txt.append(f"{float(self._radius[ii]):.2f}")
                        else:
                            txt.append(f"{int(ii//int(self._group_size))},{int(ii%int(self._group_size))}")
                    self._labels_text = txt
                self.text_labels.text = self._labels_text
                self.text_labels.pos = (self._pos[idx] + np.array([0.02, 0.02, 0.02], dtype=np.float32)[None, :]).astype(np.float32)
                self.text_labels.visible = True

        self.canvas.update()

    def _mouse_to_world_xy(self, mouse_pos, z=0.0):
        # Works best with top-down orthographic camera.
        # mouse_pos is in canvas pixels (x,y).
        tr = self.view.scene.transform
        p = tr.imap((mouse_pos[0], mouse_pos[1], float(z)))
        return np.array([p[0], p[1]], dtype=np.float32)

    def _pick_idx_from_mouse(self, pos):
        if self._pos.shape[0] == 0:
            return -1
        if self._pick_mode == '2d':
            xy = self._mouse_to_world_xy(pos, z=0.0)
            valid = np.isfinite(self._pos).all(axis=1)
            if self._fixed is not None:
                valid &= (~self._fixed)
            if not np.any(valid):
                return -1
            d2 = np.sum((self._pos[:, :2] - xy[None, :]) ** 2, axis=1)
            d2 = np.where(valid, d2, 1e30)
            return int(np.argmin(d2))
        r0, rd = self._ray_from_mouse(pos)
        i, _ = self._pick_idx_from_ray(r0, rd)
        return i

    def _on_mouse_press(self, ev):
        if ev.button in (2, 3):
            if self._selection_mode:
                # Create selection rectangle (Line visual) on first use
                if self.selection_rect is None:
                    self.selection_rect = visuals.Line(parent=self.view.scene, color=(0.2, 0.6, 1.0, 1.0), width=2.0)
                    self.selection_rect.visible = False
                # Start selection rectangle - use same method as picking for consistency
                r0, rd = self._ray_from_mouse(ev.pos)
                p = self._intersect_ray_plane(r0, rd, np.zeros(3), np.array([0,0,1]))
                if p is not None:
                    self._selection_start = np.array([p[0], p[1]], dtype=np.float32)
                else:
                    self._selection_start = self._mouse_to_world_xy(ev.pos, z=0.0)
                self._selection_end = self._selection_start.copy()
                self.selection_rect.visible = True
                self._update_selection_rect()
                ev.handled = True
                return
            i = self._pick_idx_from_mouse(ev.pos)
            if i >= 0:
                self.sig_rmb_remove.emit(i)
                ev.handled = True
                return
            self._rmb_down = True
            self._rmb_last = np.array(ev.pos, dtype=np.float32)
            self._cam_print('rmb_down')
            ev.handled = True
            return
        if ev.button != 1:
            return

        # External lock: suppress all drag (e.g. Ring mode)
        if self.lock_drag:
            i = self._pick_idx_from_mouse(ev.pos)
            self._pick_idx = i  # still allow pick detection
            self._pick_active = False
            return

        # In selection mode, always drag selected atoms regardless of where you click
        if self._selection_mode and self._selected_indices:
            self._pick_active = True
            self._pick_idx = -1  # No specific atom clicked
            self._pick_z = 0.0
            # Store initial positions of all selected atoms
            self._selected_initial_pos = {idx: self._pos[idx].copy() for idx in self._selected_indices}
            # Store initial mouse position for delta calculation
            r0, rd = self._ray_from_mouse(ev.pos)
            self._drag_start_mouse = np.array(ev.pos, dtype=np.float32)
            if int(self._cam_debug) > 0:
                print(f"[DRAG] down selected={len(self._selected_indices)} atoms (selection mode)")
            ev.handled = True
            return

        i = self._pick_idx_from_mouse(ev.pos)
        if i < 0:
            return

        # If atoms are selected and we click on one of them, drag all selected
        if self._selected_indices and i in self._selected_indices:
            self._pick_active = True
            self._pick_idx = i  # Track which atom was clicked for delta calculation
            self._pick_z = 0.0 if self._clamp_xy else float(self._pos[i, 2])
            # Store initial positions of all selected atoms
            self._selected_initial_pos = {idx: self._pos[idx].copy() for idx in self._selected_indices}
            self.sig_atom_picked.emit(i)
            self.sig_drag_state.emit(1, i, self._pos[i].copy())
            if int(self._cam_debug) > 0:
                print(f"[DRAG] down selected={len(self._selected_indices)} atoms, anchor idx={int(i)}")
            if self._pick_mode == '3d':
                r0, rd = self._ray_from_mouse(ev.pos)
                self._drag_plane_p0 = self._pos[i].copy()
                self._drag_plane_n = rd.copy()
            ev.handled = True
            return

        if self.is_fixed(i):
            # still allow pick, but not drag
            self._pick_active = False
            self._pick_idx = i
            self.sig_atom_picked.emit(i)
            ev.handled = True
            return

        self._pick_active = True
        self._pick_idx = i
        self._pick_z = 0.0 if self._clamp_xy else float(self._pos[i, 2])
        self.sig_atom_picked.emit(i)
        self.sig_drag_state.emit(1, i, self._pos[i].copy())
        if int(self._cam_debug) > 0:
            print(f"[DRAG] down idx={int(i)} pos=({self._pos[i,0]:.3f},{self._pos[i,1]:.3f},{self._pos[i,2]:.3f}) mode={self._pick_mode}")

        if self._pick_mode == '3d':
            r0, rd = self._ray_from_mouse(ev.pos)
            self._drag_plane_p0 = self._pos[i].copy()
            self._drag_plane_n = rd.copy()
        ev.handled = True

    def _on_mouse_release(self, ev):
        if ev.button in (2, 3):
            if self._selection_mode and self.selection_rect is not None and self.selection_rect.visible:
                # Finalize selection
                self.selection_rect.visible = False
                self._finalize_selection()
                ev.handled = True
                return
            self._rmb_down = False
            self._rmb_last = None
            self._cam_print('rmb_up')
            ev.handled = True
            return
        self._pick_active = False
        if self._pick_idx >= 0:
            self.sig_drag_state.emit(0, int(self._pick_idx), self._pos[int(self._pick_idx)].copy())
            if int(self._cam_debug) > 0:
                i = int(self._pick_idx)
                print(f"[DRAG] up idx={i} pos=({self._pos[i,0]:.3f},{self._pos[i,1]:.3f},{self._pos[i,2]:.3f})")
        self._pick_idx = -1
        self._drag_plane_p0 = None
        self._drag_plane_n = None
        # Clean up selected initial positions
        if hasattr(self, '_selected_initial_pos'):
            del self._selected_initial_pos

    def _on_mouse_move(self, ev):
        if self._selection_mode and self.selection_rect is not None and self.selection_rect.visible:
            # Update selection rectangle - use same method as picking for consistency
            r0, rd = self._ray_from_mouse(ev.pos)
            p = self._intersect_ray_plane(r0, rd, np.zeros(3), np.array([0,0,1]))
            if p is not None:
                self._selection_end = np.array([p[0], p[1]], dtype=np.float32)
            else:
                self._selection_end = self._mouse_to_world_xy(ev.pos, z=0.0)
            self._update_selection_rect()
            ev.handled = True
            return
        if self._rmb_down:
            if self._rmb_last is not None:
                cur = np.array(ev.pos, dtype=np.float32)
                d = cur - self._rmb_last
                self._rmb_last = cur
                self._cam_rotate(d[0], d[1])
            ev.handled = True
            return
        if not self._pick_active:
            return
        i = self._pick_idx

        # Use authoritative geometry directly if backend is available
        if self.backend is not None:
            p = self.backend.sys.apos
        else:
            p = self._pos.copy()

        if self._pick_mode == '2d':
            # Use ray casting for consistent coordinate handling with axis widgets
            r0, rd = self._ray_from_mouse(ev.pos)
            new_xy = self._intersect_ray_plane(r0, rd, np.zeros(3), np.array([0,0,1]))
            if new_xy is not None:
                # If dragging selected atoms, move all of them
                if self._selected_indices and hasattr(self, '_selected_initial_pos'):
                    if i >= 0:
                        # Clicked on an atom - use its position as reference
                        delta = new_xy[:2] - self._selected_initial_pos[i][:2]
                    else:
                        # Selection mode - calculate delta from mouse movement
                        r0_start, rd_start = self._ray_from_mouse(self._drag_start_mouse)
                        start_xy = self._intersect_ray_plane(r0_start, rd_start, np.zeros(3), np.array([0,0,1]))
                        if start_xy is not None:
                            delta = new_xy[:2] - start_xy[:2]
                        else:
                            delta = np.array([0.0, 0.0])
                    for idx in self._selected_indices:
                        p[idx, 0] = self._selected_initial_pos[idx][0] + delta[0]
                        p[idx, 1] = self._selected_initial_pos[idx][1] + delta[1]
                        p[idx, 2] = self._selected_initial_pos[idx][2]
                else:
                    p[i, 0] = new_xy[0]
                    p[i, 1] = new_xy[1]
                    p[i, 2] = self._pick_z
        else:
            if (self._drag_plane_p0 is None) or (self._drag_plane_n is None):
                return
            r0, rd = self._ray_from_mouse(ev.pos)
            x = self._intersect_ray_plane(r0, rd, self._drag_plane_p0, self._drag_plane_n)
            if x is None:
                return
            if self._clamp_xy:
                x[2] = 0.0
            # If dragging selected atoms, move all of them
            if self._selected_indices and hasattr(self, '_selected_initial_pos'):
                delta = x - self._selected_initial_pos[i]
                for idx in self._selected_indices:
                    p[idx] = self._selected_initial_pos[idx] + delta
            else:
                p[i, :] = x

        # If using backend proxy, update _pos for rendering (but backend is authoritative)
        if self.backend is not None:
            self._pos = self.backend.sys.apos.astype(np.float32)
        else:
            self._pos = p
        # Emit signal for parent to track drag position
        self.sig_atom_moved.emit(int(i), self._pos[i].copy())
        self._redraw()
        ev.handled = True

    def _on_mouse_wheel(self, ev):
        # Manual zoom (do not rely on camera.interactive)
        delta = None
        raw = {}
        if hasattr(ev, 'delta') and (ev.delta is not None):
            raw['delta'] = ev.delta
            d = ev.delta
            try:
                delta = float(d[1])
            except Exception:
                try:
                    delta = float(d)
                except Exception:
                    delta = None
        elif hasattr(ev, 'delta_y'):
            raw['delta_y'] = ev.delta_y
            delta = float(ev.delta_y)
        elif hasattr(ev, 'dy'):
            raw['dy'] = ev.dy
            delta = float(ev.dy)
        elif hasattr(ev, 'step'):
            raw['step'] = ev.step
            delta = float(ev.step)

        # fallback: if tuple and y is 0, try x
        if (delta is not None) and (abs(delta) < 1e-12) and isinstance(raw.get('delta', None), (tuple, list)):
            try:
                delta = float(raw['delta'][0])
            except Exception:
                pass

        # normalize common wheel conventions (some give +-120 per notch)
        if delta is None:
            if int(self._cam_debug) > 0:
                print(f"[WHEEL] no-delta fields={list(raw.keys())}")
            ev.handled = True
            return

        if abs(delta) > 50.0:
            delta /= 120.0
        if int(self._cam_debug) > 0:
            print(f"[WHEEL] delta={float(delta):.6g} raw={raw}")
        if abs(delta) < 1e-12:
            ev.handled = True
            return
        self._cam_zoom(delta)
        ev.handled = True

    def _on_key_press(self, ev):
        """Handle keyboard events for camera panning with arrow keys."""
        if int(self._cam_debug) > 0:
            print(f"[KEY] key={ev.key} text={ev.text}")
        if ev.key == 'ArrowUp':
            self._cam_pan(0, 1)
            ev.handled = True
        elif ev.key == 'ArrowDown':
            self._cam_pan(0, -1)
            ev.handled = True
        elif ev.key == 'ArrowLeft':
            self._cam_pan(-1, 0)
            ev.handled = True
        elif ev.key == 'ArrowRight':
            self._cam_pan(1, 0)
            ev.handled = True
        elif ev.key in ('Up', 'Down', 'Left', 'Right'):
            # Try alternative key names
            if ev.key == 'Up':
                self._cam_pan(0, 1)
            elif ev.key == 'Down':
                self._cam_pan(0, -1)
            elif ev.key == 'Left':
                self._cam_pan(-1, 0)
            elif ev.key == 'Right':
                self._cam_pan(1, 0)
            ev.handled = True

    def _update_selection_rect(self):
        """Update selection rectangle visualization (Line visual)."""
        if self._selection_start is None or self._selection_end is None:
            return
        x0, y0 = self._selection_start
        x1, y1 = self._selection_end
        # Ensure proper ordering
        x_min, x_max = min(x0, x1), max(x0, x1)
        y_min, y_max = min(y0, y1), max(y0, y1)
        # Create rectangle vertices (5 points to close the loop)
        vertices = np.array([
            [x_min, y_min, 0],
            [x_max, y_min, 0],
            [x_max, y_max, 0],
            [x_min, y_max, 0],
            [x_min, y_min, 0]
        ], dtype=np.float32)
        self.selection_rect.set_data(pos=vertices)
        self.canvas.update()

    def _finalize_selection(self):
        """Finalize selection and select atoms within rectangle."""
        if self._selection_start is None or self._selection_end is None:
            return
        x0, y0 = self._selection_start
        x1, y1 = self._selection_end
        x_min, x_max = min(x0, x1), max(x0, x1)
        y_min, y_max = min(y0, y1), max(y0, y1)

        # Find atoms within rectangle
        selected = set()
        for i in range(len(self._pos)):
            x, y = self._pos[i, 0], self._pos[i, 1]
            if x_min <= x <= x_max and y_min <= y <= y_max:
                selected.add(i)

        # Update selection
        self._selected_indices = selected
        self._highlight_selected()
        self.sig_selection_changed.emit(selected)

        self._selection_start = None
        self._selection_end = None

    def _highlight_selected(self):
        """Highlight selected atoms by changing their color."""
        if self._colors is None:
            return
        # Restore original colors (only if sizes match)
        if self._selected_colors_backup is not None and self._selected_colors_backup.shape == self._colors.shape:
            self._colors[:] = self._selected_colors_backup
        # Store original colors if first selection or if sizes don't match
        else:
            self._selected_colors_backup = self._colors.copy()
        # Highlight selected atoms
        for i in self._selected_indices:
            if i < len(self._colors):
                self._colors[i] = (1.0, 0.5, 0.0, 1.0)  # Orange highlight
        self.atom_markers.set_data(self._pos, edge_color=None, face_color=self._colors, size=self._sizes)
        self.canvas.update()

    def set_selection_mode(self, enabled):
        """Enable or disable selection mode."""
        self._selection_mode = enabled
        # Clear selection when exiting selection mode
        if not enabled:
            self.clear_selection()

    def get_selected_indices(self):
        """Return set of selected atom indices."""
        return self._selected_indices.copy()

    def set_selected_indices(self, indices):
        """Set selected atom indices."""
        self._selected_indices = set(indices)
        self._highlight_selected()
        self.sig_selection_changed.emit(self._selected_indices)

    def clear_selection(self):
        """Clear selection and restore original colors."""
        self._selected_indices.clear()
        if self._selected_colors_backup is not None:
            # Check shape consistency before restoring (atoms may have been added/removed)
            if self._selected_colors_backup.shape == self._colors.shape:
                self._colors[:] = self._selected_colors_backup
                self.atom_markers.set_data(self._pos, edge_color=None, face_color=self._colors, size=self._sizes)
            else:
                # Shape mismatch - refresh colors from backend
                self._selected_colors_backup = None
        self.canvas.update()
        self.sig_selection_changed.emit(set())


def normalize_scalar_field(vals, vmin=None, vmax=None, symmetric=False):
    a = np.asarray(vals, dtype=np.float64)
    finite = np.isfinite(a)
    if not np.any(finite):
        raise ValueError("normalize_scalar_field(): field has no finite values")
    if symmetric:
        if vmin is None and vmax is None:
            vmax = float(np.max(np.abs(a[finite])))
            vmin = -vmax
        elif vmin is None:
            vmin = -float(vmax)
        elif vmax is None:
            vmax = -float(vmin)
    else:
        if vmin is None:
            vmin = float(np.min(a[finite]))
        if vmax is None:
            vmax = float(np.max(a[finite]))
    if not np.isfinite(vmin) or not np.isfinite(vmax):
        raise ValueError(f"normalize_scalar_field(): invalid limits vmin={vmin} vmax={vmax}")
    if vmax <= vmin:
        if vmax == vmin:
            out = np.zeros_like(a, dtype=np.float32)
            out[~finite] = 0.0
            return out, float(vmin), float(vmax)
        raise ValueError(f"normalize_scalar_field(): vmax={vmax} must be > vmin={vmin}")
    out = np.zeros_like(a, dtype=np.float32)
    out[finite] = np.clip((a[finite] - vmin) / (vmax - vmin), 0.0, 1.0)
    out[~finite] = 0.0
    return out, float(vmin), float(vmax)


def make_grid_mesh_data(xs, ys, zs, colors=None, mask=None):
    xs = np.asarray(xs, dtype=np.float32)
    ys = np.asarray(ys, dtype=np.float32)
    zs = np.asarray(zs, dtype=np.float32)
    if zs.shape != (len(xs), len(ys)):
        raise ValueError(f"make_grid_mesh_data(): zs.shape={zs.shape} expected ({len(xs)},{len(ys)})")
    if mask is None:
        mask = np.isfinite(zs)
    else:
        mask = np.asarray(mask, dtype=bool)
        if mask.shape != zs.shape:
            raise ValueError(f"make_grid_mesh_data(): mask.shape={mask.shape} expected {zs.shape}")
    X, Y = np.meshgrid(xs, ys, indexing='ij')
    verts = np.stack([X, Y, zs], axis=2).reshape(-1, 3).astype(np.float32)
    if colors is not None:
        cols = np.asarray(colors, dtype=np.float32)
        if cols.shape[:2] != zs.shape:
            raise ValueError(f"make_grid_mesh_data(): colors.shape[:2]={cols.shape[:2]} expected {zs.shape}")
        cols = cols.reshape(-1, cols.shape[-1]).astype(np.float32)
    else:
        cols = None
    faces = []
    for ix in range(len(xs) - 1):
        for iy in range(len(ys) - 1):
            i00 = ix * len(ys) + iy
            i10 = (ix + 1) * len(ys) + iy
            i01 = ix * len(ys) + (iy + 1)
            i11 = (ix + 1) * len(ys) + (iy + 1)
            if mask[ix, iy] and mask[ix + 1, iy] and mask[ix + 1, iy + 1]:
                faces.append((i00, i10, i11))
            if mask[ix, iy] and mask[ix + 1, iy + 1] and mask[ix, iy + 1]:
                faces.append((i00, i11, i01))
    faces = np.asarray(faces, dtype=np.uint32)
    return verts, faces, cols


def colormap_rgba(vals, cmap='coolwarm', vmin=None, vmax=None, symmetric=False, alpha=1.0):
    t, vmin, vmax = normalize_scalar_field(vals, vmin=vmin, vmax=vmax, symmetric=symmetric)
    cm = Colormap(cmap) if isinstance(cmap, (list, tuple)) else vispy.color.get_colormap(cmap)
    mapped = cm.map(t.ravel())
    rgba = np.asarray(mapped, dtype=np.float32).reshape(t.shape + (4,))
    rgba[..., 3] = float(alpha)
    return rgba, vmin, vmax


def make_surface_mesh(xs, ys, zs, scalar=None, cmap='coolwarm', vmin=None, vmax=None, symmetric=False, alpha=1.0, mask=None):
    if scalar is None:
        rgba = None
        clim = None
    else:
        rgba, vmin, vmax = colormap_rgba(scalar, cmap=cmap, vmin=vmin, vmax=vmax, symmetric=symmetric, alpha=alpha)
        clim = (vmin, vmax)
    verts, faces, cols = make_grid_mesh_data(xs, ys, zs, colors=rgba, mask=mask)
    return {'vertices': verts, 'faces': faces, 'vertex_colors': cols, 'clim': clim}


def create_surface_visual(parent, mesh_data, shading='smooth'):
    v = np.asarray(mesh_data['vertices'], dtype=np.float32)
    f = np.asarray(mesh_data['faces'], dtype=np.uint32)
    if len(v) == 0 or len(f) == 0:
        raise ValueError(f"create_surface_visual(): empty mesh vertices={v.shape} faces={f.shape}")
    vc = mesh_data.get('vertex_colors', None)
    if vc is not None:
        return visuals.Mesh(vertices=v, faces=f, vertex_colors=np.asarray(vc, dtype=np.float32), shading=shading, parent=parent)
    return visuals.Mesh(vertices=v, faces=f, color=(0.7, 0.7, 0.9, 1.0), shading=shading, parent=parent)


def render_surface_png(out_path, mesh_data, atom_points=None, atom_colors=None, atom_sizes=None, title=None, bgcolor='white', azimuth=-60.0, elevation=35.0, scale=1.2):
    canvas = scene.SceneCanvas(keys=None, bgcolor=bgcolor, show=False, size=(1200, 900))
    view = canvas.central_widget.add_view()
    view.camera = scene.TurntableCamera(fov=0.0, elevation=float(elevation), azimuth=float(azimuth))
    mesh = create_surface_visual(view.scene, mesh_data, shading='smooth')
    mesh.set_gl_state('translucent', depth_test=True)
    if atom_points is not None:
        pts = np.asarray(atom_points, dtype=np.float32)
        if pts.ndim != 2 or pts.shape[1] != 3:
            raise ValueError(f"render_surface_png(): atom_points.shape={pts.shape} expected (n,3)")
        mk = visuals.Markers(parent=view.scene)
        mk.set_data(pts, face_color=atom_colors if atom_colors is not None else (0.1, 0.1, 0.1, 0.5), size=atom_sizes if atom_sizes is not None else 8.0, edge_width=0.0)
    if title:
        visuals.Text(text=str(title), pos=np.array([[0.0, 0.0, 0.0]], dtype=np.float32), color='black', font_size=12, parent=view.scene)
    vv = np.asarray(mesh_data['vertices'], dtype=np.float32)
    ctr = vv.mean(axis=0)
    ext = np.max(vv, axis=0) - np.min(vv, axis=0)
    rad = float(np.max(ext[:2]))
    if rad <= 1e-6:
        rad = 1.0
    view.camera.center = tuple(ctr.tolist())
    view.camera.scale_factor = float(rad * scale)
    img = canvas.render(alpha=False)
    from vispy.io import write_png
    write_png(out_path, img)
    return img


def create_heatmap_window(data_2d, extent, title="Heatmap", cmap='bwr', symmetric=True, atom_pos=None, atom_types=None):
    """Create a VisPy window to display 2D heatmap (orbital/density) with optional atom overlay.

    Args:
        data_2d: 2D numpy array (ny, nx) of scalar values where data[i,j] corresponds to (x[j], y[i])
        extent: [xmin, xmax, ymin, ymax] in world coordinates
        title: Window title
        cmap: Colormap name ('bwr', 'hot', 'viridis', etc.)
        symmetric: If True, colormap is symmetric around zero
        atom_pos: Optional (n,3) array of atom positions
        atom_types: Optional array of atom types for coloring

    Returns:
        (canvas, view) tuple for further manipulation if needed
    """
    from PyQt5 import QtWidgets
    data = np.asarray(data_2d, dtype=np.float32)
    ny, nx = data.shape

    # Create colormap
    if symmetric:
        vmax = max(abs(np.min(data)), abs(np.max(data)))
        if vmax < 1e-30:
            vmax = 1.0
        vmin = -vmax
    else:
        vmin, vmax = np.min(data), np.max(data)
        if vmax - vmin < 1e-30:
            vmax = vmin + 1.0

    # Generate RGBA colors
    rgba, _, _ = colormap_rgba(data.ravel(), cmap=cmap, vmin=vmin, vmax=vmax, symmetric=symmetric, alpha=1.0)
    rgba = rgba.reshape(ny, nx, 4).astype(np.float32)

    # Create mesh vertices for image (quad grid)
    xmin, xmax, ymin, ymax = extent
    xs = np.linspace(xmin, xmax, nx)
    ys = np.linspace(ymin, ymax, ny)
    X, Y = np.meshgrid(xs, ys)

    # Create vertices (nx*ny grid points)
    verts = np.stack([X.ravel(), Y.ravel(), np.zeros_like(X.ravel())], axis=1).astype(np.float32)

    # Create faces (two triangles per grid cell)
    n_cells = (nx - 1) * (ny - 1)
    faces = np.zeros((n_cells * 2, 3), dtype=np.uint32)
    for i in range(nx - 1):
        for j in range(ny - 1):
            cell_idx = i * (ny - 1) + j
            v00 = i * ny + j
            v01 = i * ny + (j + 1)
            v10 = (i + 1) * ny + j
            v11 = (i + 1) * ny + (j + 1)
            faces[cell_idx * 2 + 0] = [v00, v01, v10]
            faces[cell_idx * 2 + 1] = [v01, v11, v10]

    # Map colors to vertices (use corner colors)
    vertex_colors = rgba.reshape(nx * ny, 4)

    # Create window
    canvas = scene.SceneCanvas(keys=None, bgcolor='white', show=True, size=(800, 600))
    canvas.title = str(title)
    view = canvas.central_widget.add_view()
    view.camera = scene.PanZoomCamera(aspect=1.0)

    # Create heatmap mesh
    mesh = visuals.Mesh(vertices=verts, faces=faces, vertex_colors=vertex_colors, shading='flat', parent=view.scene)
    mesh.set_gl_state('translucent', depth_test=False)

    # Add atoms if provided
    if atom_pos is not None:
        pos = np.asarray(atom_pos, dtype=np.float32)
        print(f"DEBUG create_heatmap_window atom_pos: pos.shape={pos.shape}, pos.ndim={pos.ndim}")
        if pos.ndim != 2 or pos.shape[1] != 3:
            raise ValueError(f"atom_pos.shape={pos.shape} expected (n,3)")
        # Project to 2D (use x,y, set z=0 for visibility)
        pos_2d = pos.copy()
        pos_2d[:, 2] = 0.0

        # Color atoms by type if types provided, else green
        if atom_types is not None:
            from pyBall import elements
            colors = []
            for atype in atom_types:
                c = elements.getColor(atype)
                colors.append((c[0], c[1], c[2], 1.0))
            colors = np.array(colors, dtype=np.float32)
            print(f"DEBUG create_heatmap_window atom_types: atom_types={atom_types}, colors.shape={colors.shape}")
        else:
            colors = (0.0, 0.5, 0.0, 1.0)
            print(f"DEBUG create_heatmap_window atom_types: using default colors")

        print(f"DEBUG create_heatmap_window: calling atom_markers.set_data with pos_2d.shape={pos_2d.shape}, face_color type={type(colors)}")
        atom_markers = visuals.Markers(parent=view.scene)
        atom_markers.set_data(pos_2d, face_color=colors, size=5.0, edge_width=0.5, edge_color='black')

    # Center camera on extent
    cx = (xmin + xmax) / 2.0
    cy = (ymin + ymax) / 2.0
    w = xmax - xmin
    h = ymax - ymin
    view.camera.center = (cx, cy, 0.0)
    view.camera.rect = (-w/2, -h/2, w, h)

    return canvas, view
