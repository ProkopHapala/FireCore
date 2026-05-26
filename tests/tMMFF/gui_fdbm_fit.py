#!/usr/bin/env python3
"""Interactive GUI for FDBM GridFF fitting — tune P_H, P_O, q_H, dz with live 2x2 plot.

Uses shared backend from pyBall/OCL/Surface_utils.py (same as gen_gridff_nacl_gpu.py).
Run from tests/tMMFF/:
    python gui_fdbm_fit.py [--orient Na|Cl|hollow]
"""

import os, sys, json, argparse
import numpy as np

def _repo_root():
    here = os.path.dirname(os.path.abspath(__file__))
    return os.path.realpath(os.path.join(here, '..', '..'))

sys.path.insert(0, _repo_root())

from PyQt5 import QtWidgets, QtCore
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from pyBall.OCL import Surface_utils as su
from pyBall.OCL.Surface_utils import load_substrate_xyz_with_lvec

bohr2ang = 0.5291772109217
hartree2ev = 27.211396132


def load_gridff_and_data(orient='Na', Eint_max=0.5, kT_weight=0.2):
    """Load GridFF, DFT data, compute weights and pre-sample panels. Returns everything needed for GUI."""
    ROOT = _repo_root()
    dft_data_dir = os.path.join(ROOT, 'tests', 'tSurf', 'small_mols_NaCl_New')

    # Load GridFF + metadata
    gridff_path = os.path.join(dft_data_dir, 'Bspline_PLQd.npy')
    meta_path = os.path.join(dft_data_dir, 'Bspline_PLQd_meta.json')
    if not os.path.exists(gridff_path):
        raise RuntimeError(f'GridFF not found: {gridff_path}\nRun gen_gridff_nacl_gpu.py first.')
    gridff = np.load(gridff_path)
    with open(meta_path) as f:
        meta = json.load(f)
    g0 = np.array(meta['g0'], dtype=np.float32)
    dg = tuple(meta['dg'])

    # Load substrate for z_top
    sub_path = os.path.join(dft_data_dir, '0-geoms', 'NaCl.xyz')
    sub_data = load_substrate_xyz_with_lvec(sub_path)
    z_top_sub = sub_data['apos'][:, 2].max()

    # Load DFT data
    dft_paths = [
        os.path.join(ROOT, 'tests', 'tMMFF', 'H2O-H_matched.npz'),
        os.path.join(ROOT, 'tests', 'tMMFF', 'H2O-O_matched.npz'),
    ]
    data = su.load_dft_scan_data(dft_paths, z_top_sub)

    # Weights
    fit_scan_configs = [('H2O-H', 0, 0), ('H2O-O', 6, 0)]
    weights = su.compute_fit_weights(data['E_int'], data['mol_tag'], data['ix'], data['iy'],
                                      data['orientz'], fit_scan_configs, scan_tilt=0,
                                      Eint_max=Eint_max, kT_weight=kT_weight)
    data['weights'] = weights

    # Load fitted coefficients as initial values
    fit_json = os.path.join(ROOT, 'tests', 'tMMFF', 'out_fdbm_dft_gridff', 'fitted_PLQ_coeffs.json')
    if os.path.exists(fit_json):
        with open(fit_json) as f:
            fc = json.load(f)
        P_H0 = fc.get('diagnostic_P_H', 7.877)
        P_O0 = fc.get('diagnostic_P_O', 27.778)
        q_H0 = fc.get('resp_q_H', 0.2)
        dz0 = fc.get('fit_dz_shift', 0.0)
    else:
        P_H0 = 7.877; P_O0 = 27.778; q_H0 = 0.2; dz0 = 0.0

    # Panel config
    panel_mols  = ['H2O-H', 'H2O-H', 'H2O-O', 'H2O-O']
    panel_sites = ['Na', 'Cl', 'Na', 'Cl']
    panel_ix    = [6, 0, 6, 0]
    panel_iy    = [0, 0, 0, 0]
    scan_tilt = 0

    print(f'Pre-sampling GridFF for {orient} orientation (this takes a moment)...')
    panels = su.prepare_scan_panel_data(data, gridff, g0, dg, panel_mols, panel_ix, panel_iy, panel_sites,
                                         orient, scan_tilt, dz0)
    print('Done. Starting GUI.')

    return dict(panels=panels, data=data, gridff=gridff, g0=g0, dg=dg,
                panel_mols=panel_mols, panel_sites=panel_sites, panel_ix=panel_ix, panel_iy=panel_iy,
                orient=orient, scan_tilt=scan_tilt,
                P_H0=P_H0, P_O0=P_O0, q_H0=q_H0, dz0=dz0)


class FDBMFitGUI(QtWidgets.QMainWindow):
    def __init__(self, ctx):
        super().__init__()
        self.ctx = ctx
        self.panels = ctx['panels']
        self.setWindowTitle('FDBM GridFF Fit — Interactive P_H / P_O tuning')
        self.resize(1200, 900)

        central = QtWidgets.QWidget(); self.setCentralWidget(central)
        layout = QtWidgets.QVBoxLayout(central)

        # Matplotlib canvas
        self.fig = Figure(figsize=(12, 9))
        self.canvas = FigureCanvas(self.fig)
        layout.addWidget(self.canvas)

        # Spinboxes row (mouse wheel support)
        spinbox_layout = QtWidgets.QGridLayout()
        layout.addLayout(spinbox_layout)

        self.spinboxes = {}
        params = [
            ('P_H',  0.0, 1000.0, ctx['P_H0'], 0.5),
            ('P_O',  0.0, 1000.0, ctx['P_O0'], 0.5),
            ('H_H',  -10000.0, 10000.0, 0.0, 0.5),
            ('H_O',  -10000.0, 10000.0, 0.0, 0.5),
            ('q_H',  -2.0, 2.0, ctx['q_H0'], 0.01),
            ('beta', 0.5, 4.0, 1.0, 0.05),
            ('dz',   -2.0, 2.0, ctx['dz0'],  0.01),
        ]
        for col, (name, vmin, vmax, v0, step) in enumerate(params):
            lbl = QtWidgets.QLabel(f'{name}:')
            spinbox = QtWidgets.QDoubleSpinBox()
            spinbox.setRange(vmin, vmax)
            spinbox.setSingleStep(step)
            spinbox.setValue(v0)
            spinbox.setDecimals(3 if step < 0.1 else 2)
            spinbox.valueChanged.connect(self._update_plot)
            self.spinboxes[name] = spinbox
            spinbox_layout.addWidget(lbl, 0, col)
            spinbox_layout.addWidget(spinbox, 1, col)

        # View toggle (total vs per-atom)
        view_layout = QtWidgets.QHBoxLayout()
        self.view_total = QtWidgets.QRadioButton('Total (summed)')
        self.view_per_atom = QtWidgets.QRadioButton('Per-atom (H1, H2, O)')
        self.view_total.setChecked(True)
        self.view_total.toggled.connect(self._resample)
        self.view_per_atom.toggled.connect(self._resample)
        view_layout.addWidget(self.view_total)
        view_layout.addWidget(self.view_per_atom)
        layout.addLayout(view_layout)

        # Resample button (for dz changes)
        btn_resample = QtWidgets.QPushButton('Resample GridFF (after dz change)')
        btn_resample.clicked.connect(self._resample)
        layout.addWidget(btn_resample)

        # Save/Load buttons
        btn_layout = QtWidgets.QHBoxLayout()
        btn_save = QtWidgets.QPushButton('Save Params to JSON')
        btn_save.clicked.connect(self._save_params)
        btn_load = QtWidgets.QPushButton('Load Params from JSON')
        btn_load.clicked.connect(self._load_params)
        btn_layout.addWidget(btn_save)
        btn_layout.addWidget(btn_load)
        layout.addLayout(btn_layout)

        # RMS label
        self.rms_label = QtWidgets.QLabel('Total weighted RMS: ---')
        self.rms_label.setStyleSheet('font-size: 14px; font-weight: bold;')
        layout.addWidget(self.rms_label)

        self._update_plot()

    def _get_val(self, name):
        return self.spinboxes[name].value()

    def _save_params(self):
        """Save current parameter values to JSON file."""
        params = {name: self._get_val(name) for name in self.spinboxes.keys()}
        filepath, _ = QtWidgets.QFileDialog.getSaveFileName(self, 'Save Parameters', '', 'JSON files (*.json)')
        if filepath:
            with open(filepath, 'w') as f:
                json.dump(params, f, indent=2)
            print(f'Saved parameters to {filepath}')

    def _load_params(self):
        """Load parameter values from JSON file."""
        filepath, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Load Parameters', '', 'JSON files (*.json)')
        if filepath:
            with open(filepath, 'r') as f:
                params = json.load(f)
            for name, value in params.items():
                if name in self.spinboxes:
                    self.spinboxes[name].setValue(value)
            print(f'Loaded parameters from {filepath}')
            self._update_plot()

    def _resample(self):
        """Re-sample GridFF panels with new dz or view mode."""
        dz = self._get_val('dz')
        per_atom = self.view_per_atom.isChecked()
        ctx = self.ctx
        print(f'Re-sampling with dz={dz:.3f}, per_atom={per_atom}...')
        self.panels = su.prepare_scan_panel_data(
            ctx['data'], ctx['gridff'], ctx['g0'], ctx['dg'],
            ctx['panel_mols'], ctx['panel_ix'], ctx['panel_iy'], ctx['panel_sites'],
            ctx['orient'], ctx['scan_tilt'], dz, per_atom=per_atom)
        print('Done.')
        self._update_plot()

    def _update_plot(self):
        P_H = self._get_val('P_H'); P_O = self._get_val('P_O'); H_H = self._get_val('H_H'); H_O = self._get_val('H_O')
        q_H = self._get_val('q_H'); beta = self._get_val('beta')
        per_atom = self.view_per_atom.isChecked()
        self.fig.clear()
        # 2x3 layout: energy plots (2x2) + geometry (2x1)
        gs = self.fig.add_gridspec(2, 3)
        axs_energy = [self.fig.add_subplot(gs[i, j]) for i in range(2) for j in range(2)]
        axs_geom = [self.fig.add_subplot(gs[i, 2]) for i in range(2)]
        total_rms_sum = 0.0; total_n = 0
        for pidx, panel in enumerate(self.panels):
            ax = axs_energy[pidx]; axr = ax.twinx()
            if panel is None: continue
            Edft_s = panel['E_int_dft']; z_s = panel['z_s']; w_s = panel['w_s']
            if per_atom:
                Em_dict = su.compute_model_Eint_per_atom(panel, P_H, P_O, q_H, beta, H_H, H_O)
                Em_s = Em_dict['H1'] + Em_dict['H2'] + Em_dict['O']
                ax.plot(z_s, Edft_s, ls=':', lw=1.5, color='k', label='DFT')
                ax.plot(z_s, Em_s,   ls='-', lw=0.8, color='k', label='Total Model')
                ax.plot(z_s, Em_dict['H1'], ls='--', lw=0.5, color='b', alpha=0.6, label='H1')
                ax.plot(z_s, Em_dict['H2'], ls='--', lw=0.5, color='c', alpha=0.6, label='H2')
                ax.plot(z_s, Em_dict['O'],  ls='--', lw=0.5, color='m', alpha=0.6, label='O')
            else:
                Em_s = su.compute_model_Eint(panel, P_H, P_O, q_H, beta, H_H, H_O)
                ax.plot(z_s, Edft_s, ls=':', lw=1.5, color='k', label='DFT')
                ax.plot(z_s, Em_s,   ls='-', lw=0.8, color='k', label='Model')
            resid_w = (Em_s - Edft_s) * w_s
            rms = np.sqrt(np.mean(resid_w**2))
            total_rms_sum += np.sum(resid_w**2); total_n += len(resid_w)
            axr.plot(z_s, resid_w, ls='-', lw=0.6, color='r', alpha=0.5, label=f'w*(Em-Ed) RMS={rms:.4f}')
            ax.axhline(0, color='k', lw=0.3, alpha=0.3)
            axr.axhline(0, color='k', lw=0.3, alpha=0.1)
            ax.set_title(panel['title'], fontsize=9)
            ax.set_xlabel('z [A]'); ax.set_ylabel('E_int [eV]')
            axr.set_ylabel('w*resid [eV]', color='r'); axr.set_ylim(-0.5, 0.5); ax.set_ylim(-0.5, 0.5)
            ax.legend(loc='upper left', fontsize=6); axr.legend(loc='upper right', fontsize=7)
            ax.grid(True, alpha=0.2)
            # Plot geometry (map 4 panels to 2 geometry subplots)
            geom_idx = 0 if pidx < 2 else 1
            axg = axs_geom[geom_idx]
            if 'apos' in panel:
                apos = panel['apos']  # (nconf, natoms, 3)
                names = panel.get('names', ['O', 'H', 'H'])
                # Plot the first configuration as reference
                for iatom, name in enumerate(names):
                    color = 'r' if name == 'O' else 'b'
                    axg.scatter(apos[0, iatom, 0], apos[0, iatom, 2], c=color, s=100, label=name)
                    axg.text(apos[0, iatom, 0], apos[0, iatom, 2], name, fontsize=8)
                axg.set_xlabel('x [A]'); axg.set_ylabel('z [A]')
                axg.set_title('Geometry (xz plane)', fontsize=8)
                axg.grid(True, alpha=0.3)
                axg.set_aspect('equal')
                axg.legend(fontsize=6)
        dz = self._get_val('dz')
        total_rms = np.sqrt(total_rms_sum / total_n) if total_n > 0 else 0.0
        view_tag = 'per-atom' if per_atom else 'total'
        self.fig.suptitle(f'P_H={P_H:.3f} P_O={P_O:.3f} H_H={H_H:.3f} H_O={H_O:.3f} q_H={q_H:.3f} beta={beta:.2f} dz={dz:.3f} view={view_tag} | wRMS={total_rms:.4f}', fontsize=11)
        self.fig.tight_layout(rect=(0, 0, 1, 0.96))
        self.canvas.draw()
        self.rms_label.setText(f'Total weighted RMS: {total_rms:.6f}')


def main():
    parser = argparse.ArgumentParser(description='Interactive FDBM GridFF fitting GUI')
    parser.add_argument('--orient', type=str, default='Na', choices=['Na', 'Cl', 'hollow'])
    parser.add_argument('--Eint_max', type=float, default=0.5)
    parser.add_argument('--kT', type=float, default=0.2)
    args = parser.parse_args()

    ctx = load_gridff_and_data(orient=args.orient, Eint_max=args.Eint_max, kT_weight=args.kT)

    app = QtWidgets.QApplication(sys.argv)
    win = FDBMFitGUI(ctx)
    win.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
