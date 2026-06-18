import os, sys
import numpy as np

from PyQt5 import QtWidgets, QtCore

_THIS_DIR = os.path.abspath(os.path.dirname(__file__))
_SHARED_DIR = os.path.abspath(os.path.join(_THIS_DIR, '..', 'shared'))
_REPO_ROOT = os.path.abspath(os.path.join(_THIS_DIR, os.pardir))
_PYBALL_DIR = os.path.abspath(os.path.join(_THIS_DIR, '..', '..'))
for _p in (_PYBALL_DIR, _REPO_ROOT, _THIS_DIR, _SHARED_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from MolGUI_common import colors_for_enames, sizes_for_enames
from VispyUtils import AtomScene

from RRsp3 import RRsp3, build_neighs_bk_from_bonds, make_bk_slots_clustered, make_exclusions_1st_2nd
from XPTB_utils import pack_molecules_contiguous, make_h2o_geometry, masses_from_elems, load_xyz, perturb_state
from RRsp3_utils import load_molecule_any_xyz, make_ports_from_neighs, quat_rotate_vec, reorder_nodes_first


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

        self.sp_debug = QtWidgets.QSpinBox(); self.sp_debug.setRange(0, 3); self.sp_debug.setValue(1)
        vbox.addWidget(QtWidgets.QLabel('debug verbosity'))
        vbox.addWidget(self.sp_debug)

        self.cb_kprints = QtWidgets.QCheckBox('Kernel debug prints (OpenCL printf)'); self.cb_kprints.setChecked(False)
        vbox.addWidget(self.cb_kprints)

        self.cb_mode = QtWidgets.QComboBox(); self.cb_mode.addItems(['Relaxation', 'Dynamics'])
        vbox.addWidget(QtWidgets.QLabel('mode'))
        vbox.addWidget(self.cb_mode)

        self.cb_coll = QtWidgets.QCheckBox('Enable collisions'); self.cb_coll.setChecked(True)
        self.cb_ports = QtWidgets.QCheckBox('Enable ports'); self.cb_ports.setChecked(True)
        vbox.addWidget(self.cb_coll)
        vbox.addWidget(self.cb_ports)

        # ---- Molecule loading ----
        self.xyz_dir = '/home/prokophapala/git/FireCore/cpp/common_resources/xyz'
        self.xyz_dir_alt = '/home/prokophapala/git/FireCore/tests/tKekuleExplorer/out_topology'
        presets = ['h2o', 'ch2nh', 'hcooh', 'pyrrole', 'guanine', 'pentacene', 'benzene', 'naphthalene', 'anthracene', 'pyridine']
        self.cb_preset = QtWidgets.QComboBox(); self.cb_preset.addItems(presets)
        self.cb_preset.setCurrentText('anthracene')
        vbox.addWidget(QtWidgets.QLabel('molecule preset'))
        vbox.addWidget(self.cb_preset)

        self.le_custom_xyz = QtWidgets.QLineEdit()
        self.le_custom_xyz.setPlaceholderText('Or type full path to .xyz file')
        vbox.addWidget(QtWidgets.QLabel('custom xyz path'))
        vbox.addWidget(self.le_custom_xyz)

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
        self.cb_clamp = QtWidgets.QCheckBox('Constrain all atoms to Z=0 (2D mode)'); self.cb_clamp.setChecked(False)
        self.cb_pick3d = QtWidgets.QCheckBox('3D picking/drag (ortho)'); self.cb_pick3d.setChecked(True)
        vbox.addWidget(self.cb_lock)
        vbox.addWidget(self.cb_clamp)
        vbox.addWidget(self.cb_pick3d)

        self.sp_drag_max = QtWidgets.QDoubleSpinBox(); self.sp_drag_max.setRange(0.0, 1e+6); self.sp_drag_max.setSingleStep(0.5); self.sp_drag_max.setDecimals(4); self.sp_drag_max.setValue(5.0)
        vbox.addWidget(QtWidgets.QLabel('drag max step (3D)'))
        vbox.addWidget(self.sp_drag_max)

        self.cb_color_group = QtWidgets.QCheckBox('Color by groups'); self.cb_color_group.setChecked(False)
        vbox.addWidget(self.cb_color_group)

        self.cb_show_rad = QtWidgets.QCheckBox('Show collision radius'); self.cb_show_rad.setChecked(False)
        vbox.addWidget(self.cb_show_rad)

        self.sp_rcoll = QtWidgets.QDoubleSpinBox(); self.sp_rcoll.setRange(0.0, 100.0); self.sp_rcoll.setSingleStep(0.1); self.sp_rcoll.setDecimals(4); self.sp_rcoll.setValue(1.0)
        vbox.addWidget(QtWidgets.QLabel('collision radius R_coll'))
        vbox.addWidget(self.sp_rcoll)

        self.cb_show_bbox = QtWidgets.QCheckBox('Show group bboxes'); self.cb_show_bbox.setChecked(True)
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

        self.cb_neigh_bonds = QtWidgets.QCheckBox('Bonds: neighs (global list)'); self.cb_neigh_bonds.setChecked(False)
        self.cb_port_tips = QtWidgets.QCheckBox('Bonds: node -> port tip'); self.cb_port_tips.setChecked(False)
        self.cb_port_targets = QtWidgets.QCheckBox('Bonds: port tip -> target atom'); self.cb_port_targets.setChecked(False)
        self.cb_dpos = QtWidgets.QCheckBox('Debug: dpos (total)'); self.cb_dpos.setChecked(False)
        self.cb_dpos_neigh = QtWidgets.QCheckBox('Debug: dpos_neigh (slot recoils)'); self.cb_dpos_neigh.setChecked(False)
        vbox.addWidget(self.cb_neigh_bonds)
        vbox.addWidget(self.cb_port_tips)
        vbox.addWidget(self.cb_port_targets)
        vbox.addWidget(self.cb_dpos)
        vbox.addWidget(self.cb_dpos_neigh)

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
        vbox.addWidget(QtWidgets.QLabel('spring scale (dt)'))
        vbox.addWidget(self.sp_dt)

        self.sp_damp = QtWidgets.QDoubleSpinBox(); self.sp_damp.setRange(0.0, 1.0); self.sp_damp.setSingleStep(0.01); self.sp_damp.setDecimals(4); self.sp_damp.setValue(1.0)
        vbox.addWidget(QtWidgets.QLabel('damping (dynamics)'))
        vbox.addWidget(self.sp_damp)

        self.status = QtWidgets.QLabel('')
        vbox.addWidget(self.status)
        vbox.addStretch()

        hbox.addWidget(ctrl)

        self.scene.sig_atom_picked.connect(self.on_pick)
        self.scene.sig_drag_state.connect(self.on_drag_state)
        self.scene.sig_atom_moved.connect(self.on_atom_moved)

        self.scene.set_camera_debug(2)

        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self.on_timer)

        self._picked = -1
        self._fixed = None
        self._drag_active = False
        self._drag_idx = -1
        self._drag_pos = None
        self._drag_invm_backup = 0.0
        self._drag_fixed_backup = False
        self._last_build_opts = ()  # cache to avoid recompiling on every step

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
        self.cb_neigh_bonds.stateChanged.connect(self.on_viz_changed)
        self.cb_port_tips.stateChanged.connect(self.on_viz_changed)
        self.cb_port_targets.stateChanged.connect(self.on_viz_changed)
        self.cb_dpos.stateChanged.connect(self.on_viz_changed)
        self.cb_dpos_neigh.stateChanged.connect(self.on_viz_changed)
        self.cb_labels.currentIndexChanged.connect(self.on_viz_changed)
        self.sp_rcoll.valueChanged.connect(self.on_radius_changed)
        self.sp_zoom.valueChanged.connect(self.on_zoom_changed)
        self.on_view_mode_changed()
        self.on_viz_changed()

    def _debug_level(self):
        try:
            return int(self.sp_debug.value())
        except Exception:
            return 0

    def _cluster_stats(self):
        """Return per-cluster stats for real atoms only.

        Returns:
            cog (ng,3), mn (ng,3), mx (ng,3), nreal (ng,)
        """
        gs = int(self.group_size)
        ng = int(self.natoms // gs)
        cog = np.zeros((ng, 3), dtype=np.float32)
        mn = np.zeros((ng, 3), dtype=np.float32)
        mx = np.zeros((ng, 3), dtype=np.float32)
        nreal = np.zeros((ng,), dtype=np.int32)
        for ig in range(ng):
            a0 = ig * gs
            sl = slice(a0, a0 + gs)
            m = np.asarray(self.real[sl], dtype=bool)
            if not np.any(m):
                cog[ig, :] = np.nan; mn[ig, :] = np.nan; mx[ig, :] = np.nan
                continue
            p = np.asarray(self.pos[sl, :], dtype=np.float32)[m, :]
            cog[ig, :] = np.mean(p, axis=0)
            mn[ig, :] = np.min(p, axis=0)
            mx[ig, :] = np.max(p, axis=0)
            nreal[ig] = int(p.shape[0])
        return cog, mn, mx, nreal

    def _dump_cluster_stats(self, tag, *, max_groups=16):
        lvl = self._debug_level()
        if lvl <= 0:
            return
        cog, mn, mx, nreal = self._cluster_stats()
        ng = int(cog.shape[0])
        print(f"[{tag}] clusters ng={ng} (showing first {min(ng, max_groups)})")
        for ig in range(min(ng, int(max_groups))):
            c = cog[ig]; a = mn[ig]; b = mx[ig]
            print(f"  ig={ig:3d} nreal={int(nreal[ig]):3d} cog=({c[0]: .3f},{c[1]: .3f},{c[2]: .3f}) min=({a[0]: .3f},{a[1]: .3f},{a[2]: .3f}) max=({b[0]: .3f},{b[1]: .3f},{b[2]: .3f})")

    def _assert_real_finite(self, tag):
        """Fail loudly if any REAL atom has non-finite position."""
        p = np.asarray(self.pos, dtype=np.float32)
        bad = ~np.isfinite(p).all(axis=1)
        bad_real = np.where(bad & self.real)[0]
        if bad_real.size:
            print(f"[FATAL] {tag}: non-finite REAL atom positions at idx={bad_real[:32].tolist()} (showing up to 32)")
            self._dump_cluster_stats(f"{tag}/cluster_stats")
            raise RuntimeError(f"{tag}: non-finite REAL atom positions")

    def on_run_toggled(self, checked):
        if bool(checked):
            self.btn_run.setText('Stop')
            self.timer.start(int(self.sp_dtick.value()))
        else:
            self.btn_run.setText('Run')
            self.timer.stop()

    def on_timer(self):
        # Keep GUI responsive; one step per tick.
        # Atom dragging: on_drag_state pins the atom by fixmask during drag,
        # on_atom_moved updates its position. No need for drag-hold logic here.
        self.on_step()

    def _init_system(self):
        group_size = 64
        custom = self.le_custom_xyz.text().strip()
        if custom and os.path.isfile(custom):
            xyz_path = custom
            name = os.path.splitext(os.path.basename(custom))[0].lower()
        else:
            name = str(self.cb_preset.currentText()).strip().lower()
            xyz_path = os.path.join(self.xyz_dir, f"{name}.xyz")
            if not os.path.isfile(xyz_path):
                xyz_path = os.path.join(self.xyz_dir_alt, f"{name}.xyz")
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

        bonds_np = np.asarray(packed['bonds'], dtype=np.int32)
        self.scene.set_data(self.pos, colors=colors, sizes=sizes, bonds=bonds_np)
        # strict: never render padding atoms (they are invalid and should carry NaNs)
        self.scene.set_render_mask(self.real)
        self.scene.set_group_size(self.group_size)
        self.scene.set_radius(self.rad)
        self.scene.set_fixed_mask(self._fixed)
        bmin, bmax = self.sim.download_bboxes()
        self._set_scene_bboxes_from_device(bmin, bmax)
        self._update_debug_links_from_device()
        self._update_debug_overlays_from_device()
        self.status.setText(f'natoms={self.natoms} groups={self.natoms//group_size} nnode_per_group={self.nnode_per_group}')

    def on_reload(self):
        if self.btn_run.isChecked():
            self.btn_run.setChecked(False)
            self.on_run_toggled(False)
        self._init_system()

    def on_pick(self, idx):
        self._picked = int(idx)
        self.status.setText(f'picked idx={idx} elem={self.elems[idx]} pad={bool(self.is_pad[idx])} fixed={bool(self._fixed[idx])}')

    def on_drag_state(self, active, idx, pos3):
        """Handle drag start/stop.

        IMPORTANT: do NOT set invm=0 for real atoms.
        RRsp3.upload_state() uses invm==0 as a diagnostic marker and overwrites pos->NaN.
        Pinned/dragged atoms must be constrained by fixmask only.
        """
        was_active = self._drag_active
        self._drag_active = bool(active)
        self._drag_idx = int(idx)
        self._drag_pos = None if pos3 is None else np.asarray(pos3, dtype=np.float32).copy()
        
        i = int(idx)
        if i < 0 or i >= self.natoms or bool(self.is_pad[i]):
            return
            
        if active and not was_active:
            # Drag start: constrain atom by fixmask (NOT invm=0) and reset solver momentum
            self._drag_invm_backup = float(self.invm[i])
            self._drag_fixed_backup = bool(self._fixed[i])
            self._fixed[i] = True
            self.sim.reset_momentum()
            self.sim.upload_state(self.pos, self.invm, quat=self.quat)
            self._upload_fixmask()
            print(f"[DRAG-START] idx={i} fixed atom (fixmask), invm kept={self.invm[i]:.4f}")
            if self._debug_level() >= 2:
                self._dump_cluster_stats('DRAG-START')
        elif not active and was_active:
            # Drag end: restore original mobility
            self._fixed[i] = self._drag_fixed_backup
            self.sim.reset_momentum()
            self.sim.upload_state(self.pos, self.invm, quat=self.quat)
            self._upload_fixmask()
            print(f"[DRAG-END] idx={i} restored fixed={self._fixed[i]}")
            if self._debug_level() >= 2:
                self._dump_cluster_stats('DRAG-END')
            
    def on_atom_moved(self, idx, pos3):
        """Handle atom position updates during drag (from mouse movement)."""
        if not self._drag_active or idx != self._drag_idx:
            return
        i = int(idx)
        if i < 0 or i >= self.natoms or bool(self.is_pad[i]):
            return
        # Update position during drag
        new_pos = np.asarray(pos3, dtype=np.float32)
        if bool(self.cb_pick3d.isChecked()):
            max_step = float(self.sp_drag_max.value())
            if max_step > 0.0 and self._drag_pos is not None:
                dp = new_pos - self._drag_pos
                d = float(np.linalg.norm(dp))
                if d > max_step:
                    print(f"[DRAG-CLAMP] idx={i} step={d:.3f} > {max_step:.3f}; clamping")
                    new_pos = self._drag_pos + dp * (max_step / (d + 1e-16))
        self.pos[i, :] = new_pos
        self._drag_pos = new_pos.copy()
        # Upload to GPU so solver uses this as the fixed position.
        # Reset momentum because this is a discontinuous constraint target change.
        self.sim.reset_momentum()
        self.sim.upload_state(self.pos, self.invm, quat=self.quat)

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
        self.scene.set_show_neigh_bonds(bool(self.cb_neigh_bonds.isChecked()))
        self.scene.set_show_port_tips(bool(self.cb_port_tips.isChecked()))
        self.scene.set_show_port_targets(bool(self.cb_port_targets.isChecked()))
        self.scene.set_show_dpos(bool(self.cb_dpos.isChecked()))
        self.scene.set_show_dpos_neigh(bool(self.cb_dpos_neigh.isChecked()))
        # update bbox mode immediately from last downloaded device bboxes
        if hasattr(self, '_bmin_dev') and hasattr(self, '_bmax_dev'):
            self._set_scene_bboxes_from_device(self._bmin_dev, self._bmax_dev)
        self._update_debug_overlays_from_device()

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
        gid_in = []
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
                    gid_in.append(int(ig))

        seg_h = []
        gid_h = []
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

                    gid_h.append(int(ig))

        if len(seg_in) == 0:
            self.scene.set_inbox_links(None)
        else:
            self.scene.set_inbox_links(np.asarray(seg_in, dtype=np.float32), gids=np.asarray(gid_in, dtype=np.int32))
        if len(seg_h) == 0:
            self.scene.set_halo_links(None)
        else:
            self.scene.set_halo_links(np.asarray(seg_h, dtype=np.float32), gids=np.asarray(gid_h, dtype=np.int32))

    def _update_debug_overlays_from_device(self):
        if not hasattr(self, 'sim'):
            return
        if not hasattr(self, 'pos'):
            return
        if not hasattr(self, 'quat'):
            return
        pos = np.asarray(self.pos, dtype=np.float32)
        quat = np.asarray(self.quat, dtype=np.float32)
        gs = int(self.group_size)
        ng = int(self.natoms // gs)

        # neigh bonds (global)
        if bool(self.cb_neigh_bonds.isChecked()):
            neighs = np.asarray(self.neighs, dtype=np.int32)
            segs = []
            for i in range(self.natoms):
                if bool(self.is_pad[i]):
                    continue
                for k in range(4):
                    j = int(neighs[i, k])
                    if j < 0:
                        continue
                    if j <= i:
                        continue
                    if j >= self.natoms:
                        raise RuntimeError(f"neighs out of range: i={i} k={k} j={j} natoms={self.natoms}")
                    if bool(self.is_pad[j]):
                        continue
                    segs.append(pos[i]); segs.append(pos[j])
            self.scene.set_neigh_bonds(None if len(segs) == 0 else np.asarray(segs, dtype=np.float32))
        else:
            self.scene.set_neigh_bonds(None)

        # port data (only meaningful when ports enabled and buffers allocated)
        ports_ok = bool(self.cb_ports.isChecked())
        if ports_ok and (self.sim.cl_port_local is None):
            ports_ok = False

        if ports_ok and (bool(self.cb_port_tips.isChecked()) or bool(self.cb_port_targets.isChecked()) or bool(self.cb_dpos.isChecked()) or bool(self.cb_dpos_neigh.isChecked())):
            port_local = self.sim.download_port_local(nnode_per_group=self.nnode_per_group)[:, :, :3].copy()  # (nnode_tot,4,3)
            bk = self.sim.download_bkSlots()  # (natoms,4)
        else:
            port_local = None
            bk = None

        # node -> tip
        if ports_ok and bool(self.cb_port_tips.isChecked()):
            segs = []
            nnode_pg = int(self.nnode_per_group)
            for ig in range(ng):
                a0 = ig * gs
                for il in range(nnode_pg):
                    ia = a0 + il
                    if bool(self.is_pad[ia]):
                        continue
                    inode = ig * nnode_pg + il
                    qi = quat[ia]
                    xi = pos[ia]
                    rloc = port_local[inode]
                    rarm = quat_rotate_vec(qi[None, :], rloc).reshape(4, 3)
                    tips = xi[None, :] + rarm
                    for k in range(4):
                        segs.append(xi); segs.append(tips[k])
            self.scene.set_port_tips(None if len(segs) == 0 else np.asarray(segs, dtype=np.float32))
        else:
            self.scene.set_port_tips(None)

        # tip -> target atom via bkSlots
        if ports_ok and bool(self.cb_port_targets.isChecked()):
            segs = []
            nnode_pg = int(self.nnode_per_group)
            for ia in range(self.natoms):
                if bool(self.is_pad[ia]):
                    continue
                for k in range(4):
                    slot = int(bk[ia, k])
                    if slot < 0:
                        continue
                    inode = slot >> 2
                    pk = slot & 3
                    ig = inode // nnode_pg
                    il = inode - ig * nnode_pg
                    if ig < 0 or ig >= ng:
                        raise RuntimeError(f"bkSlots bad inode->ig: ia={ia} slot={slot} inode={inode} ig={ig} ng={ng}")
                    node_atom = ig * gs + il
                    if node_atom < 0 or node_atom >= self.natoms:
                        raise RuntimeError(f"bkSlots bad node_atom: ia={ia} slot={slot} node_atom={node_atom}")
                    if bool(self.is_pad[node_atom]):
                        continue
                    tip = pos[node_atom] + quat_rotate_vec(quat[node_atom], port_local[inode, pk])
                    segs.append(tip); segs.append(pos[ia])
            self.scene.set_port_targets(None if len(segs) == 0 else np.asarray(segs, dtype=np.float32))
        else:
            self.scene.set_port_targets(None)

        # dpos (total)
        if ports_ok and bool(self.cb_dpos.isChecked()):
            dpos_coll = self.sim.download_dpos_coll()[:, :3]
            dpos_node = self.sim.download_dpos_node(nnode_per_group=self.nnode_per_group)[:, :3]
            dpos_neigh = self.sim.download_dpos_neigh(nnode_per_group=self.nnode_per_group)[:, :, :3]
            dx = np.zeros_like(pos)
            nnode_pg = int(self.nnode_per_group)
            for ig in range(ng):
                a0 = ig * gs
                inode_base = ig * nnode_pg
                for il in range(nnode_pg):
                    ia = a0 + il
                    if bool(self.is_pad[ia]):
                        continue
                    inode = inode_base + il
                    dx[ia, :] += dpos_node[inode]
            for ia in range(self.natoms):
                if bool(self.is_pad[ia]):
                    continue
                for k in range(4):
                    slot = int(bk[ia, k])
                    if slot < 0:
                        continue
                    inode = slot >> 2
                    pk = slot & 3
                    dx[ia, :] += dpos_neigh[inode, pk]
            dx += dpos_coll
            segs = []
            s = 1.0
            for ia in range(self.natoms):
                if bool(self.is_pad[ia]):
                    continue
                segs.append(pos[ia]); segs.append(pos[ia] + dx[ia] * s)
            self.scene.set_dpos(None if len(segs) == 0 else np.asarray(segs, dtype=np.float32))
        else:
            self.scene.set_dpos(None)

        # dpos_neigh vectors (slot recoils)
        if ports_ok and bool(self.cb_dpos_neigh.isChecked()):
            dpos_neigh = self.sim.download_dpos_neigh(nnode_per_group=self.nnode_per_group)[:, :, :3]
            segs = []
            nnode_pg = int(self.nnode_per_group)
            for ia in range(self.natoms):
                if bool(self.is_pad[ia]):
                    continue
                for k in range(4):
                    slot = int(bk[ia, k])
                    if slot < 0:
                        continue
                    inode = slot >> 2
                    pk = slot & 3
                    ig = inode // nnode_pg
                    il = inode - ig * nnode_pg
                    node_atom = ig * gs + il
                    if bool(self.is_pad[node_atom]):
                        continue
                    segs.append(pos[ia]); segs.append(pos[ia] + dpos_neigh[inode, pk])
            self.scene.set_dpos_neigh(None if len(segs) == 0 else np.asarray(segs, dtype=np.float32))
        else:
            self.scene.set_dpos_neigh(None)

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
        self._fixed[self.is_pad] = True
        self.scene.set_fixed_mask(self._fixed)
        self.sim.reset_momentum()
        self.sim.upload_state(self.pos, self.invm, quat=self.quat)
        self._upload_fixmask()
        self.status.setText(f'picked idx={i} elem={self.elems[i]} pad={bool(self.is_pad[i])} fixed={bool(self._fixed[i])}')

    def on_step(self):
        # rebuild program with toggles only when changed
        build_opts = []
        if not self.cb_coll.isChecked():
            build_opts.append('-DENABLE_COLL=0')
        if not self.cb_ports.isChecked():
            build_opts.append('-DENABLE_PORT=0')

        if bool(self.cb_kprints.isChecked()):
            lvl = int(self.sp_debug.value())
            if lvl < 1: lvl = 1
            gid0 = 0
            wg = -1
            if int(self._picked) >= 0:
                gid0 = (int(self._picked) // int(self.group_size)) * int(self.group_size)
                wg = int(self._picked) // int(self.group_size)
            build_opts.append('-DENABLE_DEBUG_PRINTS')
            build_opts.append(f'-DDEBUG_VERBOSITY={lvl}')
            build_opts.append(f'-DDEBUG_COMPONENTS={1|2|4|8}')
            build_opts.append(f'-DDEBUG_TARGET_WG={wg}')
            build_opts.append(f'-DDEBUG_GID_START={gid0}')
            build_opts.append(f'-DDEBUG_GID_END={gid0 + int(self.group_size)}')
        build_opts_tuple = tuple(build_opts)

        if build_opts_tuple != self._last_build_opts:
            # recreate sim only when build options change
            self.sim = RRsp3(self.natoms, group_size=self.group_size, prefer_gpu=True, build_options=build_opts)
            self.sim.upload_state(self.pos, self.invm, quat=self.quat)
            self.sim.upload_radius(self.rad)
            self.sim.upload_neighs_and_exclusions(self.neighs, self.excl1, self.excl2)
            self.sim.upload_cluster_ports(self.port_local_atoms, self.K_atoms, nnode_per_group=self.nnode_per_group)
            self.sim.upload_bkSlots(self.bkSlots)
            self._upload_fixmask()
            self.sim.run_bboxes_and_topology(bbox_margin=0.5)
            self._last_build_opts = build_opts_tuple

        if str(self.cb_mode.currentText()).strip().lower().startswith('dyn'):
            self.sim.step_dynamics(nnode_per_group=self.nnode_per_group, dt=float(self.sp_dt.value()), k_coll=float(self.sp_kcoll.value()), relaxation=float(self.sp_relax.value()), bbox_margin=0.5, damp=float(self.sp_damp.value()))
        else:
            self.sim.step_cluster(nnode_per_group=self.nnode_per_group, dt=float(self.sp_dt.value()), k_coll=float(self.sp_kcoll.value()), relaxation=float(self.sp_relax.value()), bbox_margin=0.5)
        pos4, quat4 = self.sim.download_pos_quat()
        # strict: padding atoms remain NaN (invalid) always
        self.pos = pos4[:, :3].copy()
        self.quat = quat4.copy()
        self._assert_real_finite('on_step/post_download')
        if self._debug_level() >= 3:
            self._dump_cluster_stats('STEP')
        self.scene.update_positions(self.pos)
        bmin, bmax = self.sim.download_bboxes()
        self._set_scene_bboxes_from_device(bmin, bmax)
        self._update_debug_links_from_device()
        self._update_debug_overlays_from_device()


def main():
    app = QtWidgets.QApplication(sys.argv)
    w = RRsp3VisDebug()
    w.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
