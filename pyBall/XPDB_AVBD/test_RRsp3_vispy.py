import os, sys
import numpy as np

from PyQt5 import QtWidgets, QtCore

_THIS_DIR = os.path.abspath(os.path.dirname(__file__))
_REPO_ROOT = os.path.abspath(os.path.join(_THIS_DIR, os.pardir, os.pardir))
_PYBALL_DIR = os.path.abspath(os.path.join(_REPO_ROOT, 'pyBall'))
for _p in (_REPO_ROOT, _PYBALL_DIR, _THIS_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from MolGUI_common import colors_for_enames, sizes_for_enames
from VispyUtils import AtomScene

from RRsp3 import RRsp3, build_neighs_bk_from_bonds, make_bk_slots_clustered, make_exclusions_1st_2nd
from XPTB_utils import pack_molecules_contiguous, make_h2o_geometry, masses_from_elems, load_xyz, perturb_state


def reorder_nodes_first(elems, xyz, bonds, *, is_node=None):
    elems = list(elems)
    xyz = np.asarray(xyz, dtype=np.float32)
    n = int(len(elems))
    if xyz.shape != (n, 3):
        raise ValueError(f"reorder_nodes_first: xyz.shape={xyz.shape} expected ({n},3)")
    if is_node is None:
        is_node = np.array([e != 'H' for e in elems], dtype=bool)
    is_node = np.asarray(is_node, dtype=bool)
    if is_node.shape != (n,):
        raise ValueError(f"reorder_nodes_first: is_node.shape={is_node.shape} expected ({n},)")
    if not np.any(is_node):
        is_node[:] = True
    order = np.concatenate([np.nonzero(is_node)[0], np.nonzero(~is_node)[0]]).astype(np.int32)
    perm_inv = np.empty((n,), dtype=np.int32)
    perm_inv[order] = np.arange(n, dtype=np.int32)
    elems2 = [elems[i] for i in order]
    xyz2 = xyz[order, :].copy()
    bonds2 = [(int(perm_inv[int(i)]), int(perm_inv[int(j)])) for (i, j) in bonds]
    nnode = int(np.sum(is_node))
    return elems2, xyz2, bonds2, nnode


def load_molecule_any_xyz(xyz_path):
    # Use AtomicSystem bond finding for general molecules
    from pyBall.AtomicSystem import AtomicSystem
    mol = AtomicSystem(fname=xyz_path)
    if mol.bonds is None or len(mol.bonds) == 0:
        mol.findBonds()
    bonds = [(int(b[0]), int(b[1])) for b in mol.bonds]
    elems = list(mol.enames)
    xyz = np.asarray(mol.apos[:, :3], dtype=np.float32)
    elems, xyz, bonds, nnode = reorder_nodes_first(elems, xyz, bonds)
    return elems, xyz, bonds, nnode


def make_ports_from_neighs(pos, neighs, K=200.0):
    pos = np.asarray(pos, dtype=np.float32)
    neighs = np.asarray(neighs, dtype=np.int32)
    n = int(pos.shape[0])
    port_local = np.zeros((n, 4, 4), dtype=np.float32)
    Kflat = np.zeros((n, 4), dtype=np.float32)
    for i in range(n):
        for k in range(4):
            j = int(neighs[i, k])
            if j < 0:
                continue
            port_local[i, k, :3] = pos[j] - pos[i]
            Kflat[i, k] = float(K)
    return port_local, Kflat


class RRsp3VisDebug(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle('RRsp3 Vis Debug (Vispy)')
        self.resize(1400, 900)

        self.scene = AtomScene(bgcolor='white')

        central = QtWidgets.QWidget(); self.setCentralWidget(central)
        hbox = QtWidgets.QHBoxLayout(central)
        hbox.addWidget(self.scene.widget, stretch=3)

        # ---- Controls ----
        ctrl = QtWidgets.QWidget(); ctrl.setMaximumWidth(360)
        vbox = QtWidgets.QVBoxLayout(ctrl)

        self.btn_step = QtWidgets.QPushButton('Step (1 iter)')
        self.btn_step.clicked.connect(self.on_step)
        vbox.addWidget(self.btn_step)

        self.btn_run = QtWidgets.QPushButton('Run')
        self.btn_run.setCheckable(True)
        self.btn_run.clicked.connect(self.on_run_toggled)
        vbox.addWidget(self.btn_run)

        self.sp_dtick = QtWidgets.QSpinBox(); self.sp_dtick.setRange(5, 1000); self.sp_dtick.setSingleStep(5); self.sp_dtick.setValue(30)
        vbox.addWidget(QtWidgets.QLabel('run dt [ms]'))
        vbox.addWidget(self.sp_dtick)

        self.cb_coll = QtWidgets.QCheckBox('Enable collisions'); self.cb_coll.setChecked(True)
        self.cb_ports = QtWidgets.QCheckBox('Enable ports'); self.cb_ports.setChecked(True)
        vbox.addWidget(self.cb_coll)
        vbox.addWidget(self.cb_ports)

        # ---- Molecule loading ----
        self.xyz_dir = '/home/prokop/git/FireCore/web/common_resources/xyz'
        self.cb_preset = QtWidgets.QComboBox(); self.cb_preset.addItems(['h2o', 'ch2nh', 'hcooh', 'pyrrole', 'guanine', 'pentacene'])
        self.cb_preset.setCurrentText('guanine')
        vbox.addWidget(QtWidgets.QLabel('molecule preset'))
        vbox.addWidget(self.cb_preset)

        self.sp_nmol = QtWidgets.QSpinBox(); self.sp_nmol.setRange(1, 200); self.sp_nmol.setValue(2)
        vbox.addWidget(QtWidgets.QLabel('num molecules (=num clusters)'))
        vbox.addWidget(self.sp_nmol)

        self.sp_shift = QtWidgets.QDoubleSpinBox(); self.sp_shift.setRange(0.0, 100.0); self.sp_shift.setSingleStep(0.5); self.sp_shift.setValue(8.0)
        vbox.addWidget(QtWidgets.QLabel('cluster spacing'))
        vbox.addWidget(self.sp_shift)

        self.sp_pert_pos = QtWidgets.QDoubleSpinBox(); self.sp_pert_pos.setRange(0.0, 5.0); self.sp_pert_pos.setSingleStep(0.01); self.sp_pert_pos.setDecimals(4); self.sp_pert_pos.setValue(0.0)
        self.sp_pert_rot = QtWidgets.QDoubleSpinBox(); self.sp_pert_rot.setRange(0.0, 5.0); self.sp_pert_rot.setSingleStep(0.01); self.sp_pert_rot.setDecimals(4); self.sp_pert_rot.setValue(0.0)
        vbox.addWidget(QtWidgets.QLabel('perturb pos'))
        vbox.addWidget(self.sp_pert_pos)
        vbox.addWidget(QtWidgets.QLabel('perturb rot'))
        vbox.addWidget(self.sp_pert_rot)

        self.sp_seed = QtWidgets.QSpinBox(); self.sp_seed.setRange(0, 10**9); self.sp_seed.setValue(0)
        vbox.addWidget(QtWidgets.QLabel('perturb seed'))
        vbox.addWidget(self.sp_seed)

        self.btn_reload = QtWidgets.QPushButton('Reload molecule(s)')
        self.btn_reload.clicked.connect(self.on_reload)
        vbox.addWidget(self.btn_reload)

        self.cb_lock = QtWidgets.QCheckBox('Lock top view (2D pick)'); self.cb_lock.setChecked(False)
        self.cb_clamp = QtWidgets.QCheckBox('Clamp dragged atom to XY'); self.cb_clamp.setChecked(False)
        self.cb_pick3d = QtWidgets.QCheckBox('3D picking/drag (ortho)'); self.cb_pick3d.setChecked(True)
        vbox.addWidget(self.cb_lock)
        vbox.addWidget(self.cb_clamp)
        vbox.addWidget(self.cb_pick3d)

        self.cb_color_group = QtWidgets.QCheckBox('Color by groups'); self.cb_color_group.setChecked(False)
        vbox.addWidget(self.cb_color_group)

        self.cb_show_rad = QtWidgets.QCheckBox('Show collision radius'); self.cb_show_rad.setChecked(False)
        vbox.addWidget(self.cb_show_rad)

        self.sp_rcoll = QtWidgets.QDoubleSpinBox(); self.sp_rcoll.setRange(0.0, 100.0); self.sp_rcoll.setSingleStep(0.1); self.sp_rcoll.setDecimals(4); self.sp_rcoll.setValue(1.0)
        vbox.addWidget(QtWidgets.QLabel('collision radius R_coll'))
        vbox.addWidget(self.sp_rcoll)

        self.cb_show_bbox = QtWidgets.QCheckBox('Show group bboxes'); self.cb_show_bbox.setChecked(False)
        vbox.addWidget(self.cb_show_bbox)

        self.cb_bbox_mode = QtWidgets.QComboBox(); self.cb_bbox_mode.addItems(['tight', 'overlap', 'halo'])
        vbox.addWidget(QtWidgets.QLabel('bbox mode'))
        vbox.addWidget(self.cb_bbox_mode)

        self.cb_axes = QtWidgets.QCheckBox('Show axes'); self.cb_axes.setChecked(True)
        vbox.addWidget(self.cb_axes)

        self.cb_inbox_links = QtWidgets.QCheckBox('Links: bbox center -> in-box atoms'); self.cb_inbox_links.setChecked(False)
        self.cb_halo_links = QtWidgets.QCheckBox('Links: bbox center -> halo atoms'); self.cb_halo_links.setChecked(False)
        vbox.addWidget(self.cb_inbox_links)
        vbox.addWidget(self.cb_halo_links)

        self.cb_labels = QtWidgets.QComboBox(); self.cb_labels.addItems(['none', 'global', 'local', 'pair', 'radius'])
        vbox.addWidget(QtWidgets.QLabel('labels'))
        vbox.addWidget(self.cb_labels)

        self.sp_zoom = QtWidgets.QDoubleSpinBox(); self.sp_zoom.setRange(1e-4, 1e4); self.sp_zoom.setDecimals(6); self.sp_zoom.setSingleStep(0.05); self.sp_zoom.setValue(1.0)
        vbox.addWidget(QtWidgets.QLabel('zoom (ortho scale)'))
        vbox.addWidget(self.sp_zoom)

        self.btn_reset_view = QtWidgets.QPushButton('Reset view')
        self.btn_reset_view.clicked.connect(self.on_reset_view)
        vbox.addWidget(self.btn_reset_view)

        self.btn_pin = QtWidgets.QPushButton('Pin/Unpin picked')
        self.btn_pin.clicked.connect(self.on_pin_toggle)
        vbox.addWidget(self.btn_pin)

        self.sp_relax = QtWidgets.QDoubleSpinBox(); self.sp_relax.setRange(0.0, 2.0); self.sp_relax.setSingleStep(0.05); self.sp_relax.setValue(0.5)
        vbox.addWidget(QtWidgets.QLabel('relaxation'))
        vbox.addWidget(self.sp_relax)

        self.sp_kcoll = QtWidgets.QDoubleSpinBox(); self.sp_kcoll.setRange(0.0, 500.0); self.sp_kcoll.setSingleStep(5.0); self.sp_kcoll.setValue(50.0)
        vbox.addWidget(QtWidgets.QLabel('k_coll'))
        vbox.addWidget(self.sp_kcoll)

        self.sp_dt = QtWidgets.QDoubleSpinBox(); self.sp_dt.setRange(1e-6, 10.0); self.sp_dt.setSingleStep(0.01); self.sp_dt.setDecimals(4); self.sp_dt.setValue(0.1)
        vbox.addWidget(QtWidgets.QLabel('dt'))
        vbox.addWidget(self.sp_dt)

        self.status = QtWidgets.QLabel('')
        vbox.addWidget(self.status)
        vbox.addStretch()

        hbox.addWidget(ctrl)

        self.scene.sig_atom_picked.connect(self.on_pick)
        self.scene.sig_atom_dragged.connect(self.on_drag)
        self.scene.sig_drag_state.connect(self.on_drag_state)

        self.scene.set_camera_debug(2)

        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self.on_timer)

        self._picked = -1
        self._fixed = None
        self._drag_active = False
        self._drag_idx = -1
        self._drag_pos = None

        self._init_system()

        self.cb_lock.stateChanged.connect(self.on_view_mode_changed)
        self.cb_clamp.stateChanged.connect(self.on_view_mode_changed)
        self.cb_pick3d.stateChanged.connect(self.on_view_mode_changed)
        self.cb_color_group.stateChanged.connect(self.on_viz_changed)
        self.cb_show_rad.stateChanged.connect(self.on_viz_changed)
        self.cb_show_bbox.stateChanged.connect(self.on_viz_changed)
        self.cb_bbox_mode.currentIndexChanged.connect(self.on_viz_changed)
        self.cb_axes.stateChanged.connect(self.on_viz_changed)
        self.cb_inbox_links.stateChanged.connect(self.on_viz_changed)
        self.cb_halo_links.stateChanged.connect(self.on_viz_changed)
        self.cb_labels.currentIndexChanged.connect(self.on_viz_changed)
        self.sp_rcoll.valueChanged.connect(self.on_radius_changed)
        self.sp_zoom.valueChanged.connect(self.on_zoom_changed)
        self.on_view_mode_changed()
        self.on_viz_changed()

    def on_run_toggled(self, checked):
        if bool(checked):
            self.btn_run.setText('Stop')
            self.timer.start(int(self.sp_dtick.value()))
        else:
            self.btn_run.setText('Run')
            self.timer.stop()

    def on_timer(self):
        # Keep GUI responsive; one step per tick.
        # If dragging is active, re-apply dragged atom position before and after stepping
        if self._drag_active and (self._drag_idx >= 0) and (self._drag_pos is not None):
            i = int(self._drag_idx)
            if (i >= 0) and (i < self.natoms) and (not bool(self.is_pad[i])) and (not bool(self._fixed[i])):
                print(f"[DRAG-HOLD] pre-step idx={i} pos=({self._drag_pos[0]:.3f},{self._drag_pos[1]:.3f},{self._drag_pos[2]:.3f})")
                self.pos[i, :] = self._drag_pos
                self.sim.upload_state(self.pos, self.invm, quat=self.quat)

        self.on_step()

        if self._drag_active and (self._drag_idx >= 0) and (self._drag_pos is not None):
            i = int(self._drag_idx)
            if (i >= 0) and (i < self.natoms) and (not bool(self.is_pad[i])) and (not bool(self._fixed[i])):
                print(f"[DRAG-HOLD] post-step idx={i} pos=({self._drag_pos[0]:.3f},{self._drag_pos[1]:.3f},{self._drag_pos[2]:.3f})")
                self.pos[i, :] = self._drag_pos
                self.sim.upload_state(self.pos, self.invm, quat=self.quat)
                self.scene.update_positions(self.pos)

    def _init_system(self):
        group_size = 64
        name = str(self.cb_preset.currentText()).strip().lower()
        xyz_path = os.path.join(self.xyz_dir, f"{name}.xyz")
        if os.path.isfile(xyz_path):
            try:
                elems0, xyz0, bonds0, _nnode_loaded = load_molecule_any_xyz(xyz_path)
            except Exception:
                elems0, xyz0, _q = load_xyz(xyz_path)
                bonds0 = []
        else:
            pos_h2o, _ = make_h2o_geometry(add_angle=False)
            elems0 = ['O', 'H', 'H']
            xyz0 = pos_h2o
            bonds0 = [(0, 1), (0, 2)]

        # Infer nnode the same way pack_molecules_contiguous does (deg>1 => node). Avoid mismatches.
        deg = np.zeros((len(elems0),), dtype=np.int32)
        for (i, j) in bonds0:
            deg[int(i)] += 1
            deg[int(j)] += 1
        nnode = int(np.sum(deg > 1))
        if nnode <= 0:
            nnode = len(elems0)  # fallback: all atoms as nodes if no bonds/isolated

        nmol = int(self.sp_nmol.value())
        shift = float(self.sp_shift.value())
        mols = []
        for i in range(nmol):
            sh = np.array([shift * (float(i) - 0.5 * float(nmol - 1)), 0.0, 0.0], dtype=np.float32)
            mols.append({'elems': list(elems0), 'pos': xyz0 + sh[None, :], 'bonds': list(bonds0), 'nnode': int(nnode)})

        packed = pack_molecules_contiguous(mols, group_size=group_size, nodes_first=True, pad_to_group=True)
        self.group_size = group_size
        self.natoms = int(packed['natoms_total'])
        self.pos = packed['pos'].copy()
        self.elems = packed['elems']
        self.is_pad = packed['is_padding']
        self.real = ~self.is_pad

        # masses / inv masses
        elems_real = [e for e in self.elems if e != 'X']
        m_real = masses_from_elems(elems_real)
        self.m = np.zeros((self.natoms,), dtype=np.float32)
        self.m[self.real] = m_real
        self.invm = np.zeros_like(self.m)
        self.invm[self.real] = 1.0 / self.m[self.real]
        self.invm0 = self.invm.copy()
        self._fixed = np.zeros((self.natoms,), dtype=bool)
        # padding atoms always fixed (and should not be pickable)
        self._fixed[self.is_pad] = True

        # quat (+ optional perturb)
        self.quat = np.zeros((self.natoms, 4), dtype=np.float32); self.quat[:, 3] = 1.0
        try:
            from numpy.random import default_rng
            rng = default_rng(int(self.sp_seed.value()))
            pos_init, quat_init = perturb_state(self.pos, self.quat, float(self.sp_pert_pos.value()), float(self.sp_pert_rot.value()), rng)
            # keep padding inert and non-confusing
            gs = int(self.group_size)
            ng = int(self.natoms // gs)
            for ig in range(ng):
                a0 = ig * gs
                sl = slice(a0, a0 + gs)
                pad_g = np.asarray(self.is_pad[sl], dtype=bool)
                if not np.any(pad_g):
                    continue
                real_idx = np.nonzero(~pad_g)[0]
                if real_idx.size == 0:
                    continue
                iref = a0 + int(real_idx[0])
                idx_pad = a0 + np.nonzero(pad_g)[0]
                pos_init[idx_pad, :] = pos_init[iref, :][None, :]
                quat_init[idx_pad, :] = (0.0, 0.0, 0.0, 1.0)
            self.pos = pos_init
            self.quat = quat_init
        except Exception:
            pass

        # neighs/excl
        self.neighs, _ = build_neighs_bk_from_bonds(self.natoms, packed['bonds'], max_deg=4)
        self.excl1, self.excl2 = make_exclusions_1st_2nd(self.neighs)

        self.nnode_per_group = int(packed['nnode_group'][0])
        self.bkSlots = make_bk_slots_clustered(self.neighs, group_size=group_size, nnode_per_group=self.nnode_per_group, natoms=self.natoms)

        # ports
        self.port_local_atoms, self.K_atoms = make_ports_from_neighs(self.pos, self.neighs, K=200.0)

        # radii (constant for now; will later be per-element)
        self.rad = np.zeros((self.natoms,), dtype=np.float32)
        self.rad[self.real] = float(self.sp_rcoll.value())

        self.sim = RRsp3(self.natoms, group_size=group_size, prefer_gpu=True)
        self.sim.upload_state(self.pos, self.invm, quat=self.quat)
        self.sim.upload_radius(self.rad)
        self.sim.upload_neighs_and_exclusions(self.neighs, self.excl1, self.excl2)
        self.sim.upload_cluster_ports(self.port_local_atoms, self.K_atoms, nnode_per_group=self.nnode_per_group)
        self.sim.upload_bkSlots(self.bkSlots)
        self._upload_fixmask()
        self.sim.run_bboxes_and_topology(bbox_margin=0.5)

        colors = colors_for_enames(self.elems, alpha=0.9)
        sizes = sizes_for_enames(self.elems, scale=0.25)
        bonds_np = np.array([(0, 1), (0, 2), (64, 65), (64, 66)], dtype=np.int32)

        self.scene.set_data(self.pos, colors=colors, sizes=sizes, bonds=bonds_np)
        # strict: never render padding atoms (they are invalid and should carry NaNs)
        self.scene.set_render_mask(self.real)
        self.scene.set_group_size(self.group_size)
        self.scene.set_radius(self.rad)
        self.scene.set_fixed_mask(self._fixed)
        bmin, bmax = self.sim.download_bboxes()
        self._set_scene_bboxes_from_device(bmin, bmax)
        self._update_debug_links_from_device()
        self.status.setText(f'natoms={self.natoms} groups={self.natoms//group_size} nnode_per_group={self.nnode_per_group}')

    def on_reload(self):
        if self.btn_run.isChecked():
            self.btn_run.setChecked(False)
            self.on_run_toggled(False)
        self._init_system()

    def on_pick(self, idx):
        self._picked = int(idx)
        self.status.setText(f'picked idx={idx} elem={self.elems[idx]} pad={bool(self.is_pad[idx])} fixed={bool(self._fixed[idx])}')

    def on_drag(self, idx, new_pos3):
        idx = int(idx)
        if bool(self.is_pad[idx]):
            return
        if bool(self._fixed[idx]):
            return
        self.pos[idx, :] = np.asarray(new_pos3, dtype=np.float32)
        self._drag_pos = self.pos[idx, :].copy()
        self.sim.upload_state(self.pos, self.invm, quat=self.quat)

    def on_drag_state(self, active, idx, pos3):
        self._drag_active = bool(active)
        self._drag_idx = int(idx)
        self._drag_pos = None if pos3 is None else np.asarray(pos3, dtype=np.float32).copy()

    def on_view_mode_changed(self):
        lock = bool(self.cb_lock.isChecked())
        clamp = bool(self.cb_clamp.isChecked())
        pick3d = bool(self.cb_pick3d.isChecked())
        self.scene.set_lock_top_view(lock)
        self.scene.set_clamp_xy(clamp)
        self.scene.set_pick_mode('3d' if pick3d else '2d')
        self._upload_fixmask()

    def on_viz_changed(self):
        self.scene.set_color_by_group(bool(self.cb_color_group.isChecked()))
        self.scene.set_show_radius(bool(self.cb_show_rad.isChecked()))
        self.scene.set_show_bboxes(bool(self.cb_show_bbox.isChecked()))
        self.scene.set_radius_style('disc')
        self.scene.set_label_mode(str(self.cb_labels.currentText()))
        self.scene.set_show_axes(bool(self.cb_axes.isChecked()))
        self.scene.set_show_inbox_links(bool(self.cb_inbox_links.isChecked()))
        self.scene.set_show_halo_links(bool(self.cb_halo_links.isChecked()))
        # update bbox mode immediately from last downloaded device bboxes
        if hasattr(self, '_bmin_dev') and hasattr(self, '_bmax_dev'):
            self._set_scene_bboxes_from_device(self._bmin_dev, self._bmax_dev)

    def on_zoom_changed(self, v):
        self.scene.set_zoom(float(v))

    def on_reset_view(self):
        self.scene.reset_view()
        self.sp_zoom.blockSignals(True)
        self.sp_zoom.setValue(self.scene.get_zoom())
        self.sp_zoom.blockSignals(False)

    def _set_scene_bboxes_from_device(self, bmin, bmax):
        self._bmin_dev = np.asarray(bmin, dtype=np.float32)
        self._bmax_dev = np.asarray(bmax, dtype=np.float32)
        mode = str(self.cb_bbox_mode.currentText())
        bbmin = self._bmin_dev.copy(); bbmax = self._bmax_dev.copy()
        if mode != 'tight':
            rad = self.sim.download_radius()
            rmax = float(np.max(rad[np.isfinite(rad)])) if rad.size else 0.0
            bbox_margin = 0.5
            if mode == 'overlap':
                ext = float(bbox_margin)
            else:
                ext = float(2.0 * rmax + bbox_margin)
            bbmin[:, :3] -= ext
            bbmax[:, :3] += ext
        self.scene.set_bboxes(bbmin, bbmax)

    def _update_debug_links_from_device(self):
        bmin, bmax = self._bmin_dev, self._bmax_dev
        ng = int(bmin.shape[0])
        c = 0.5 * (bmin[:, :3] + bmax[:, :3])

        seg_in = []
        if bool(self.cb_inbox_links.isChecked()):
            for ig in range(ng):
                a0 = ig * self.group_size
                a1 = a0 + self.group_size
                p = self.pos[a0:a1, :]
                m = self.real[a0:a1] & np.isfinite(p).all(axis=1)
                for k, ok in enumerate(m):
                    if not bool(ok):
                        continue
                    seg_in.append(c[ig]); seg_in.append(p[k])

        seg_h = []
        if bool(self.cb_halo_links.isChecked()):
            gi, gc = self.sim.download_ghosts()
            for ig in range(ng):
                n = int(gc[ig])
                if n <= 0:
                    continue
                for j in gi[ig, :min(n, gi.shape[1])]:
                    jj = int(j)
                    if jj < 0 or jj >= self.natoms:
                        continue
                    if bool(self.is_pad[jj]):
                        continue
                    if not np.isfinite(self.pos[jj]).all():
                        continue
                    seg_h.append(c[ig]); seg_h.append(self.pos[jj])

        self.scene.set_inbox_links(None if len(seg_in) == 0 else np.asarray(seg_in, dtype=np.float32))
        self.scene.set_halo_links(None if len(seg_h) == 0 else np.asarray(seg_h, dtype=np.float32))

    def on_radius_changed(self, v):
        self.rad[:] = 0.0
        self.rad[self.real] = float(v)
        self.sim.upload_radius(self.rad)
        self.scene.set_radius(self.rad)

    def keyPressEvent(self, ev):
        k = ev.key()
        if k == QtCore.Qt.Key_Space:
            self.btn_run.setChecked(not self.btn_run.isChecked())
            self.on_run_toggled(self.btn_run.isChecked())
            ev.accept();
            return
        if k == QtCore.Qt.Key_P:
            self.on_pin_toggle()
            ev.accept();
            return
        super().keyPressEvent(ev)

    def _upload_fixmask(self):
        if self._fixed is None:
            return
        m = np.zeros((self.natoms,), dtype=np.int32)
        # pinned atoms
        m[self._fixed] |= (1 | 2 | 4)
        # padding always pinned
        m[self.is_pad] |= (1 | 2 | 4)
        # clamp z for all atoms if enabled
        if bool(self.cb_clamp.isChecked()):
            m |= 8
        self.sim.upload_fixmask(m)

    def on_pin_toggle(self):
        i = int(self._picked)
        if i < 0:
            return
        if bool(self.is_pad[i]):
            return
        self._fixed[i] = ~self._fixed[i]
        if self._fixed[i]:
            self.invm[i] = 0.0
        else:
            self.invm[i] = self.invm0[i]
        self._fixed[self.is_pad] = True
        self.scene.set_fixed_mask(self._fixed)
        self.sim.upload_state(self.pos, self.invm, quat=self.quat)
        self._upload_fixmask()
        self.status.setText(f'picked idx={i} elem={self.elems[i]} pad={bool(self.is_pad[i])} fixed={bool(self._fixed[i])}')

    def on_step(self):
        # rebuild program with toggles (simple, safe)
        build_opts = []
        if not self.cb_coll.isChecked():
            build_opts.append('-DENABLE_COLL=0')
        if not self.cb_ports.isChecked():
            build_opts.append('-DENABLE_PORT=0')

        if build_opts:
            # recreate sim to force rebuild options
            self.sim = RRsp3(self.natoms, group_size=self.group_size, prefer_gpu=True, build_options=build_opts)
            self.sim.upload_state(self.pos, self.invm, quat=self.quat)
            self.sim.upload_radius(self.rad)
            self.sim.upload_neighs_and_exclusions(self.neighs, self.excl1, self.excl2)
            self.sim.upload_cluster_ports(self.port_local_atoms, self.K_atoms, nnode_per_group=self.nnode_per_group)
            self.sim.upload_bkSlots(self.bkSlots)
            self._upload_fixmask()
            self.sim.run_bboxes_and_topology(bbox_margin=0.5)

        self.sim.step_cluster(nnode_per_group=self.nnode_per_group, dt=float(self.sp_dt.value()), k_coll=float(self.sp_kcoll.value()), relaxation=float(self.sp_relax.value()), bbox_margin=0.5)
        pos4, quat4 = self.sim.download_pos_quat()
        # strict: padding atoms remain NaN (invalid) always
        self.pos = pos4[:, :3].copy()
        self.quat = quat4.copy()
        self.scene.update_positions(self.pos)
        bmin, bmax = self.sim.download_bboxes()
        self._set_scene_bboxes_from_device(bmin, bmax)
        self._update_debug_links_from_device()


def main():
    app = QtWidgets.QApplication(sys.argv)
    w = RRsp3VisDebug()
    w.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
