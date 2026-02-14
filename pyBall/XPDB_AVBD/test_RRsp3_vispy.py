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
from XPTB_utils import pack_molecules_contiguous, make_h2o_geometry, masses_from_elems


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

        self.cb_lock = QtWidgets.QCheckBox('Lock top view (2D pick)'); self.cb_lock.setChecked(False)
        self.cb_clamp = QtWidgets.QCheckBox('Clamp dragged atom to XY'); self.cb_clamp.setChecked(False)
        self.cb_pick3d = QtWidgets.QCheckBox('3D picking/drag (ortho)'); self.cb_pick3d.setChecked(True)
        vbox.addWidget(self.cb_lock)
        vbox.addWidget(self.cb_clamp)
        vbox.addWidget(self.cb_pick3d)

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
        self.on_view_mode_changed()

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
        pos_h2o, _ = make_h2o_geometry(add_angle=False)
        elems = ['O', 'H', 'H']
        bonds = [(0, 1), (0, 2)]
        nnode = 1

        mols = [
            {'elems': elems, 'pos': pos_h2o, 'bonds': bonds, 'nnode': nnode},
            {'elems': elems, 'pos': pos_h2o + np.array([2.2, 0.0, 0.0], dtype=np.float32), 'bonds': bonds, 'nnode': nnode},
        ]
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

        # quat
        self.quat = np.zeros((self.natoms, 4), dtype=np.float32); self.quat[:, 3] = 1.0

        # neighs/excl
        self.neighs, _ = build_neighs_bk_from_bonds(self.natoms, packed['bonds'], max_deg=4)
        self.excl1, self.excl2 = make_exclusions_1st_2nd(self.neighs)

        self.nnode_per_group = int(packed['nnode_group'][0])
        self.bkSlots = make_bk_slots_clustered(self.neighs, group_size=group_size, nnode_per_group=self.nnode_per_group, natoms=self.natoms)

        # ports
        self.port_local_atoms, self.K_atoms = make_ports_from_neighs(self.pos, self.neighs, K=200.0)

        # radii
        self.rad = np.zeros((self.natoms,), dtype=np.float32)
        self.rad[self.real] = 1.0

        self.sim = RRsp3(self.natoms, group_size=group_size, prefer_gpu=True)
        self.sim.upload_state(self.pos, self.invm, quat=self.quat)
        self.sim.upload_radius(self.rad)
        self.sim.upload_neighs_and_exclusions(self.neighs, self.excl1, self.excl2)
        self.sim.upload_cluster_ports(self.port_local_atoms, self.K_atoms, nnode_per_group=self.nnode_per_group)
        self.sim.upload_bkSlots(self.bkSlots)

        colors = colors_for_enames(self.elems, alpha=0.9)
        sizes = sizes_for_enames(self.elems, scale=0.25)
        bonds_np = np.array([(0, 1), (0, 2), (64, 65), (64, 66)], dtype=np.int32)

        self.scene.set_data(self.pos, colors=colors, sizes=sizes, bonds=bonds_np)
        self.scene.set_fixed_mask(self._fixed)
        self.status.setText(f'natoms={self.natoms} groups={self.natoms//group_size} nnode_per_group={self.nnode_per_group}')

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

        self.sim.step_cluster(nnode_per_group=self.nnode_per_group, dt=float(self.sp_dt.value()), k_coll=float(self.sp_kcoll.value()), relaxation=float(self.sp_relax.value()), bbox_margin=0.5)
        pos4, quat4 = self.sim.download_pos_quat()
        pos_new = pos4[:, :3].copy()
        # GPU padding slots may be NaN by design; keep viewer positions stable for padding
        mnan = ~np.isfinite(pos_new).all(axis=1)
        if np.any(mnan):
            pos_new[mnan, :] = self.pos[mnan, :]
        self.pos = pos_new
        self.quat = quat4.copy()
        self.scene.update_positions(self.pos)


def main():
    app = QtWidgets.QApplication(sys.argv)
    w = RRsp3VisDebug()
    w.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
