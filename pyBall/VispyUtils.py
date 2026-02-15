import numpy as np

from PyQt5 import QtCore
import vispy
vispy.use('pyqt5')
from vispy import scene
from vispy.scene import visuals


def _as_f32(x):
    return np.asarray(x, dtype=np.float32)


class AtomScene(QtCore.QObject):
    """Reusable Vispy widget for atoms (+ optional bonds) with orthographic top-down view.

    This is intended as a generic MD/FF viewer: operates on arrays (positions, colors, sizes, bonds)
    and optional per-atom vectors (forces) for visualization.

    Picking/dragging is implemented in a pseudo-2D mode:
    - camera is fixed to top-down view
    - drag moves atoms in XY plane (Z unchanged unless you set it externally)
    """

    sig_atom_picked = QtCore.pyqtSignal(int)
    sig_atom_dragged = QtCore.pyqtSignal(int, object)  # idx, new_pos3
    sig_drag_state = QtCore.pyqtSignal(int, int, object)  # active(0/1), idx, pos3

    def __init__(self, *, bgcolor='white'):
        super().__init__(parent=None)

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

        # Draw ordering: radius behind everything, then bboxes/links/lines, then atom centers, then labels.
        self.radius_markers = visuals.Markers(parent=self.view.scene)
        self.bbox_lines = visuals.Line(parent=self.view.scene, color=(0.2, 0.6, 1.0, 0.6), width=1.2, antialias=True, method='gl')
        self.inbox_lines = visuals.Line(parent=self.view.scene, color=(0.0, 0.0, 0.0, 0.35), width=1.0, antialias=True, method='gl')
        self.halo_lines  = visuals.Line(parent=self.view.scene, color=(0.8, 0.1, 0.8, 0.35), width=1.0, antialias=True, method='gl')
        self.bond_lines = visuals.Line(parent=self.view.scene, color='gray', width=1.5, antialias=True, method='gl')
        self.force_lines = visuals.Line(parent=self.view.scene, color=(1, 0, 0, 0.8), width=2.0, antialias=True, method='gl')
        self.atom_markers = visuals.Markers(parent=self.view.scene)
        self.axes = visuals.XYZAxis(parent=self.view.scene)
        self.text_labels = visuals.Text(parent=self.view.scene, color='black', font_size=10, anchor_x='left', anchor_y='bottom')

        # Enforce z-order when supported
        for o, v in enumerate((self.radius_markers, self.bbox_lines, self.inbox_lines, self.halo_lines, self.bond_lines, self.force_lines, self.atom_markers, self.axes, self.text_labels)):
            if hasattr(v, 'order'):
                v.order = int(o)

        # GL state: radius translucent and never blocks other overlays
        try:
            self.radius_markers.set_gl_state('translucent', depth_test=False)
            for v in (self.bbox_lines, self.inbox_lines, self.halo_lines, self.bond_lines, self.force_lines):
                v.set_gl_state('translucent', depth_test=False)
            self.atom_markers.set_gl_state('translucent', depth_test=False)
            self.text_labels.set_gl_state('translucent', depth_test=False)
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
        self._show_axes = True
        self._bboxes_min = None
        self._bboxes_max = None
        self._inbox_link_segs = None
        self._halo_link_segs = None

        self._label_mode = 'none'  # none|global|local|pair|radius
        self._labels_text = None

        self._marker_style = 'disc'      # vispy marker name
        self._radius_style = 'disc'

        self._pick_active = False
        self._pick_idx = -1
        self._pick_z = 0.0

        self._pick_mode = '2d'   # '2d' or '3d'
        self._lock_top_view = True
        self._clamp_xy = False
        self._fixed = None

        self._drag_plane_p0 = None
        self._drag_plane_n = None

        self.canvas.events.mouse_press.connect(self._on_mouse_press)
        self.canvas.events.mouse_release.connect(self._on_mouse_release)
        self.canvas.events.mouse_move.connect(self._on_mouse_move)
        self.canvas.events.mouse_wheel.connect(self._on_mouse_wheel)

        self._apply_camera_mode()

    @property
    def widget(self):
        return self.canvas.native

    def set_data(self, pos, *, colors=None, sizes=None, bonds=None, forces=None, force_scale=1.0):
        pos = _as_f32(pos)
        if pos.ndim != 2 or pos.shape[1] != 3:
            raise ValueError(f"AtomScene.set_data: pos.shape={pos.shape} expected (n,3)")
        self._pos = pos
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

    def set_show_inbox_links(self, enable):
        self._show_inbox_links = bool(enable)
        self._redraw()

    def set_show_halo_links(self, enable):
        self._show_halo_links = bool(enable)
        self._redraw()

    def set_show_axes(self, enable):
        self._show_axes = bool(enable)
        self.axes.visible = bool(enable)
        self.canvas.update()

    def set_inbox_links(self, segs):
        if segs is None:
            self._inbox_link_segs = None
        else:
            s = _as_f32(segs)
            if s.ndim != 2 or s.shape[1] != 3 or (s.shape[0] % 2) != 0:
                raise ValueError(f"AtomScene.set_inbox_links: segs.shape={s.shape} expected (2*m,3)")
            self._inbox_link_segs = s
        self._redraw()

    def set_halo_links(self, segs):
        if segs is None:
            self._halo_link_segs = None
        else:
            s = _as_f32(segs)
            if s.ndim != 2 or s.shape[1] != 3 or (s.shape[0] % 2) != 0:
                raise ValueError(f"AtomScene.set_halo_links: segs.shape={s.shape} expected (2*m,3)")
            self._halo_link_segs = s
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

    def _ray_from_mouse(self, mouse_pos, z0=0.0, z1=1.0):
        # mouse_pos in canvas pixels
        tr = self.view.scene.transform
        p0 = np.array(tr.imap((mouse_pos[0], mouse_pos[1], float(z0)))[:3], dtype=np.float32)
        p1 = np.array(tr.imap((mouse_pos[0], mouse_pos[1], float(z1)))[:3], dtype=np.float32)
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

    def _intersect_ray_plane(self, r0, rd, p0, n):
        denom = float(np.dot(rd, n))
        if abs(denom) < 1e-12:
            return None
        t = float(np.dot((p0 - r0), n) / denom)
        return r0 + rd * t

    def _redraw(self):
        n = self._pos.shape[0]
        if n == 0:
            self.atom_markers.set_data(np.zeros((0, 3), dtype=np.float32))
            self.radius_markers.set_data(np.zeros((0, 3), dtype=np.float32))
            self.bond_lines.set_data(np.zeros((0, 3), dtype=np.float32))
            self.force_lines.set_data(np.zeros((0, 3), dtype=np.float32))
            self.bbox_lines.set_data(np.zeros((0, 3), dtype=np.float32))
            # vispy Text requires at least one position; keep a hidden dummy.
            self.text_labels.text = ['']
            self.text_labels.pos = np.zeros((1, 3), dtype=np.float32)
            self.text_labels.visible = False
            self.canvas.update();
            return

        if self._render_mask is None:
            m = np.isfinite(self._pos).all(axis=1)
            if not np.all(m):
                raise RuntimeError("AtomScene._redraw: found NaN/Inf in positions but no render_mask set")
            self._render_mask = np.ones((n,), dtype=bool)

        # Strict rule: never render invalid atoms; if any non-masked atom is invalid -> crash loudly.
        finite = np.isfinite(self._pos).all(axis=1)
        if np.any((~finite) & (self._render_mask)):
            bad = np.where((~finite) & (self._render_mask))[0]
            raise RuntimeError(f"AtomScene._redraw: invalid pos for render_mask atoms; bad indices (first 10)={bad[:10]}")
        m = self._render_mask & finite
        idx = np.where(m)[0]

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
            connect = np.zeros((self._inbox_link_segs.shape[0],), dtype=bool); connect[0::2] = True
            self.inbox_lines.set_data(self._inbox_link_segs, connect=connect)
        else:
            self.inbox_lines.set_data(np.zeros((0, 3), dtype=np.float32))

        if self._show_halo_links and (self._halo_link_segs is not None) and (self._halo_link_segs.size > 0):
            connect = np.zeros((self._halo_link_segs.shape[0],), dtype=bool); connect[0::2] = True
            self.halo_lines.set_data(self._halo_link_segs, connect=connect)
        else:
            self.halo_lines.set_data(np.zeros((0, 3), dtype=np.float32))

        # Bonds: draw segment pairs
        if (self._bonds is not None) and (self._bonds.size > 0):
            b = self._bonds
            mb = m[b[:, 0]] & m[b[:, 1]]
            b = b[mb]
            if b.size > 0:
                segs = np.empty((b.shape[0] * 2, 3), dtype=np.float32)
                segs[0::2] = self._pos[b[:, 0]]
                segs[1::2] = self._pos[b[:, 1]]
                connect = np.zeros(segs.shape[0], dtype=bool)
                connect[0::2] = True
                self.bond_lines.set_data(segs, connect=connect, color=(0.3, 0.3, 0.3, 0.8), width=1.5)
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
            connect = np.zeros(segs.shape[0], dtype=bool)
            connect[0::2] = True
            self.force_lines.set_data(segs, connect=connect)
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

        # Labels
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
                            txt.append(f"{float(self._radius[int(ii)]):.3f}")
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
            self._rmb_down = True
            self._rmb_last = np.array(ev.pos, dtype=np.float32)
            self._cam_print('rmb_down')
            ev.handled = True
            return
        if ev.button != 1:
            return
        i = self._pick_idx_from_mouse(ev.pos)
        if i < 0:
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

    def _on_mouse_move(self, ev):
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
        p = self._pos.copy()

        if self._pick_mode == '2d':
            new_xy = self._mouse_to_world_xy(ev.pos, z=0.0)
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
            p[i, :] = x

        self._pos = p
        self._redraw()
        self.sig_atom_dragged.emit(i, p[i].copy())
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
