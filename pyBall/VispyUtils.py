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

        self.atom_markers = visuals.Markers(parent=self.view.scene)
        self.bond_lines = visuals.Line(parent=self.view.scene, color='gray', width=1.5, antialias=True, method='gl')
        self.force_lines = visuals.Line(parent=self.view.scene, color=(1, 0, 0, 0.8), width=2.0, antialias=True, method='gl')

        self._pos = np.zeros((0, 3), dtype=np.float32)
        self._colors = None
        self._sizes = None
        self._bonds = None
        self._forces = None

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
        self._sizes = None if sizes is None else _as_f32(sizes)
        self._bonds = None if bonds is None else np.asarray(bonds, dtype=np.int32)
        self._forces = None if forces is None else _as_f32(forces)
        self._force_scale = float(force_scale)
        self._redraw()

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
        self.canvas.update()
        self._cam_print('rotate')

    def _cam_zoom(self, delta):
        cam = self.view.camera
        if cam is None:
            return
        d = float(cam.distance)
        s = float(np.exp(-float(delta) * float(self._cam_zoom_speed)))
        d2 = d * s
        if d2 < 0.1:
            d2 = 0.1
        cam.distance = d2
        self.canvas.update()
        if int(self._cam_debug) > 0:
            print(f"[CAM] zoom delta={float(delta):.6g} dist:{d:.6g}->{d2:.6g}")

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
            self.bond_lines.set_data(np.zeros((0, 3), dtype=np.float32))
            self.force_lines.set_data(np.zeros((0, 3), dtype=np.float32))
            self.canvas.update();
            return

        if self._colors is None:
            face_color = (0.2, 0.2, 0.2, 1.0)
        else:
            face_color = self._colors
        if self._sizes is None:
            size = 8.0
        else:
            size = self._sizes

        self.atom_markers.set_data(self._pos, face_color=face_color, size=size, edge_width=0.5, edge_color='black')

        # Bonds: draw segment pairs
        if (self._bonds is not None) and (self._bonds.size > 0):
            b = self._bonds
            segs = np.empty((b.shape[0] * 2, 3), dtype=np.float32)
            segs[0::2] = self._pos[b[:, 0]]
            segs[1::2] = self._pos[b[:, 1]]
            connect = np.zeros(segs.shape[0], dtype=bool)
            connect[0::2] = True
            self.bond_lines.set_data(segs, connect=connect, color=(0.3, 0.3, 0.3, 0.8), width=1.5)
        else:
            self.bond_lines.set_data(np.zeros((0, 3), dtype=np.float32))

        # Forces: per-atom line from pos to pos+f*scale
        if self._forces is not None:
            f = self._forces
            if f.shape != self._pos.shape:
                raise ValueError(f"AtomScene._redraw: forces.shape={f.shape} expected {self._pos.shape}")
            segs = np.empty((n * 2, 3), dtype=np.float32)
            segs[0::2] = self._pos
            segs[1::2] = self._pos + f * self._force_scale
            connect = np.zeros(segs.shape[0], dtype=bool)
            connect[0::2] = True
            self.force_lines.set_data(segs, connect=connect)
        else:
            self.force_lines.set_data(np.zeros((0, 3), dtype=np.float32))

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
