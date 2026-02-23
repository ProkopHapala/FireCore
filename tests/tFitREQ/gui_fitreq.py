#!/usr/bin/python3 -u
import sys
import os
import numpy as np
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import QKeySequence
import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

sys.path.append("../../")
from pyBall import FitREQ as fit

def parse_dof_file(fname):
    dofs = []
    with open(fname, 'r') as f:
        for line in f:
            if line.strip().startswith('#') or not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 10:
                dofs.append({
                    'typename': parts[0],
                    'comp': int(parts[1]),
                    'Min': float(parts[2]),
                    'Max': float(parts[3]),
                    'xlo': float(parts[4]),
                    'xhi': float(parts[5]),
                    'Klo': float(parts[6]),
                    'Khi': float(parts[7]),
                    'K0': float(parts[8]),
                    'xstart': float(parts[9]),
                    'invMass': float(parts[10]) if len(parts)>10 else 1.0
                })
    return dofs

def write_dof_file(fname, dofs):
    with open(fname, 'w') as f:
        f.write("# typename comp Min Max xlo xhi Klo Khi K0 xstart invMass\n")
        for d in dofs:
            f.write(f"{d['typename']} {d['comp']} {d['Min']} {d['Max']} {d['xlo']} {d['xhi']} {d['Klo']} {d['Khi']} {d['K0']} {d['xstart']} {d['invMass']}\n")

class FitReqSession:
    def __init__(self, xyz_path, dof_file):
        self.xyz_path = xyz_path
        self.dof_file = dof_file
        
        self.dof_data = parse_dof_file(self.dof_file)
        
        fit.setVerbosity(0)
        # Regularize=1 is needed to apply K-constraints during relaxation
        fit.setup(imodel=7, EvalJ=1, WriteJ=1, Regularize=1, SaveJustElementXYZ=-1)
        fit.setGlobalParams(kMorse=1.7, Lepairs=0.9)
        fit.loadTypes()
        fit.loadDOFSelection(fname=self.dof_file)
        fit.loadXYZ(self.xyz_path, bAddEpairs=True)
        fit.getBuffs()
        
        self.dof_names = fit.loadDOFnames(self.dof_file)
        self.nDOFs = len(self.dof_names)
        
        self.Gref, self.seq, self.axis, self.distances, self.angles = fit.parse_xyz_mapping(self.xyz_path)
        self.Epanel_ref = None
        
        GRS, _, _ = fit.shift_grid(self.Gref)
        if GRS is not None:
            self.Epanel_ref = GRS.T
        self.Xpanel = np.tile(np.asarray(self.distances, dtype=float), (len(self.angles), 1))
        
    def evaluate(self, dof_values, kMorse, Lepairs):
        fit.setGlobalParams(kMorse=kMorse, Lepairs=Lepairs)
        for i in range(self.nDOFs):
            fit.DOFs[i] = dof_values[i]
        
        _, Es, _ = fit.getEs(bOmp=False, bDOFtoTypes=True, bEs=True, bFs=False)
        
        Gm = np.empty_like(self.Gref)
        Gm[:] = np.nan
        nmap = min(len(Es), len(self.seq))
        for i in range(nmap):
            idist, iang = self.seq[i]
            Gm[idist, iang] = Es[i]
            
        return Gm

class FitReqGUI(QMainWindow):
    def __init__(self, session):
        super().__init__()
        self.session = session
        self.setWindowTitle("FitREQ Interactive GUI")
        self.resize(1100, 700)
        
        # To fix plot normalization
        self.ref_e_min, self.ref_e_max = None, None
        self.ref_r_min, self.ref_r_max = None, None
        
        self.main_widget = QWidget()
        self.setCentralWidget(self.main_widget)
        self.layout = QHBoxLayout(self.main_widget)
        
        self.create_controls()
        self.create_plot()
        
        self.update_timer = QTimer()
        self.update_timer.setSingleShot(True)
        self.update_timer.timeout.connect(self.do_evaluate_and_plot)
        
        # Setup hotkey for relax
        self.shortcut_relax = QShortcut(QKeySequence("R"), self)
        self.shortcut_relax.activated.connect(self.run_relax)
        
        print("Calling evaluate and plot"); sys.stdout.flush(); self.do_evaluate_and_plot()
        
    def create_controls(self):
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFixedWidth(450)
        
        content = QWidget()
        self.ctrl_layout = QFormLayout(content)
        
        # File operations
        btn_row = QHBoxLayout()
        btn_save = QPushButton("Save DOFs")
        btn_save.clicked.connect(self.save_dofs)
        btn_load = QPushButton("Load DOFs")
        btn_load.clicked.connect(self.load_dofs)
        btn_row.addWidget(btn_save)
        btn_row.addWidget(btn_load)
        self.ctrl_layout.addRow(btn_row)
        
        # Global params
        self.spin_kMorse = QDoubleSpinBox(); self.spin_kMorse.setRange(0, 10); self.spin_kMorse.setSingleStep(0.05); self.spin_kMorse.setValue(1.7)
        self.spin_Lepairs = QDoubleSpinBox(); self.spin_Lepairs.setRange(0, 5); self.spin_Lepairs.setSingleStep(0.05); self.spin_Lepairs.setValue(0.9)
        self.chk_kcal = QCheckBox("Use kcal/mol")
        self.chk_kcal.setChecked(True)
        
        self.ctrl_layout.addRow("<b>Globals</b>", QWidget())
        self.ctrl_layout.addRow("kMorse:", self.spin_kMorse)
        self.ctrl_layout.addRow("Lepairs:", self.spin_Lepairs)
        self.ctrl_layout.addRow("", self.chk_kcal)
        
        # Relaxation settings
        self.ctrl_layout.addRow("<b>Relaxation Params</b>", QWidget())
        
        self.chk_soft_clamp = QCheckBox("Enable Soft Clamp")
        self.chk_soft_clamp.setChecked(True)
        self.spin_sc_start = QDoubleSpinBox(); self.spin_sc_start.setRange(0, 10); self.spin_sc_start.setValue(4.0)
        self.spin_sc_max = QDoubleSpinBox(); self.spin_sc_max.setRange(0, 10); self.spin_sc_max.setValue(6.0)
        
        self.chk_user_weights = QCheckBox("Enable User Weights")
        self.chk_user_weights.setChecked(True)
        self.spin_w_a = QDoubleSpinBox(); self.spin_w_a.setRange(0, 10); self.spin_w_a.setValue(1.0)
        self.spin_w_alpha = QDoubleSpinBox(); self.spin_w_alpha.setRange(0, 10); self.spin_w_alpha.setValue(4.0)
        self.spin_n_before = QSpinBox(); self.spin_n_before.setRange(1, 1000); self.spin_n_before.setValue(100)
        self.spin_emin_min = QDoubleSpinBox(); self.spin_emin_min.setRange(-1.0, 1.0); self.spin_emin_min.setSingleStep(0.01); self.spin_emin_min.setValue(-0.02)
        self.spin_nstep = QSpinBox(); self.spin_nstep.setRange(1, 10000); self.spin_nstep.setValue(100)
        
        self.ctrl_layout.addRow(self.chk_soft_clamp)
        self.ctrl_layout.addRow("SoftClamp Start:", self.spin_sc_start)
        self.ctrl_layout.addRow("SoftClamp Max:", self.spin_sc_max)
        
        self.ctrl_layout.addRow(self.chk_user_weights)
        self.ctrl_layout.addRow("Weight a:", self.spin_w_a)
        self.ctrl_layout.addRow("Weight alpha:", self.spin_w_alpha)
        self.ctrl_layout.addRow("n_before:", self.spin_n_before)
        self.ctrl_layout.addRow("Emin min:", self.spin_emin_min)
        self.ctrl_layout.addRow("n_steps:", self.spin_nstep)
        
        btn_relax = QPushButton("Relax (Hotkey: R)")
        btn_relax.setStyleSheet("background-color: #4CAF50; color: white; font-weight: bold;")
        btn_relax.clicked.connect(self.run_relax)
        self.ctrl_layout.addRow(btn_relax)
        
        # DOFs
        self.dof_spins = []
        self.k_spins = []
        self.ctrl_layout.addRow(QLabel("<b>DOF Parameters</b>"), QWidget())
        for i, d in enumerate(self.session.dof_data):
            row = QHBoxLayout()
            
            val_spin = QDoubleSpinBox()
            val_spin.setRange(-10.0, 10.0)
            val_spin.setDecimals(5)
            val_spin.setSingleStep(0.01)
            val_spin.setValue(d['xstart'])
            val_spin.valueChanged.connect(self.trigger_update)
            
            k_spin = QDoubleSpinBox()
            k_spin.setRange(0.0, 1000.0)
            k_spin.setDecimals(3)
            k_spin.setSingleStep(0.1)
            k_spin.setValue(d['K0'])
            # K constraint doesn't need to trigger immediate plot update, only used on relax
            
            row.addWidget(val_spin)
            row.addWidget(QLabel(" K:"))
            row.addWidget(k_spin)
            
            self.dof_spins.append(val_spin)
            self.k_spins.append(k_spin)
            self.ctrl_layout.addRow(self.session.dof_names[i] + ":", row)
            
        self.spin_kMorse.valueChanged.connect(self.trigger_update)
        self.spin_Lepairs.valueChanged.connect(self.trigger_update)
        self.chk_kcal.stateChanged.connect(self.trigger_update)
        
        scroll.setWidget(content)
        self.layout.addWidget(scroll)
        
    def save_dofs(self):
        fname, _ = QFileDialog.getSaveFileName(self, "Save DOF Selection", "dofSelection_saved.dat", "Data Files (*.dat);;All Files (*)")
        if fname:
            for i, d in enumerate(self.session.dof_data):
                d['xstart'] = self.dof_spins[i].value()
                d['K0'] = self.k_spins[i].value()
            write_dof_file(fname, self.session.dof_data)
            QMessageBox.information(self, "Saved", f"Successfully saved to {fname}")
            
    def load_dofs(self):
        fname, _ = QFileDialog.getOpenFileName(self, "Load DOF Selection", "", "Data Files (*.dat);;All Files (*)")
        if fname:
            new_data = parse_dof_file(fname)
            if len(new_data) == self.session.nDOFs:
                self.session.dof_data = new_data
                for i, d in enumerate(self.session.dof_data):
                    self.dof_spins[i].blockSignals(True)
                    self.k_spins[i].blockSignals(True)
                    
                    self.dof_spins[i].setValue(d['xstart'])
                    self.k_spins[i].setValue(d['K0'])
                    
                    self.dof_spins[i].blockSignals(False)
                    self.k_spins[i].blockSignals(False)
                self.trigger_update()
                QMessageBox.information(self, "Loaded", f"Successfully loaded {fname}")
            else:
                QMessageBox.warning(self, "Error", "Mismatch in number of DOFs in file.")
                
    def run_relax(self):
        print("Starting relaxation...")
        # 1. Update dof_data
        for i, d in enumerate(self.session.dof_data):
            d['xstart'] = self.dof_spins[i].value()
            d['K0'] = self.k_spins[i].value()
            
        # 2. Write tmp and reload
        tmp_dof = "_tmp_dofs.dat"
        write_dof_file(tmp_dof, self.session.dof_data)
        fit.loadDOFSelection(tmp_dof)
        fit.getBuffs()
        
        # 3. Setup global params & Soft clamp
        kM = self.spin_kMorse.value()
        Lep = self.spin_Lepairs.value()
        if self.chk_soft_clamp.isChecked():
            fit.setGlobalParams(kMorse=kM, Lepairs=Lep, 
                                softClamp_start=self.spin_sc_start.value(), 
                                softClamp_max=self.spin_sc_max.value())
        else:
            fit.setGlobalParams(kMorse=kM, Lepairs=Lep)
            
        # 4. Setup user weights
        if self.chk_user_weights.isChecked():
            Erefs, x0s = fit.read_xyz_data(self.session.xyz_path)
            weight_func = lambda E: fit.exp_weight_func(E, a=self.spin_w_a.value(), alpha=self.spin_w_alpha.value())
            weights0, lens = fit.split_and_weight_curves(Erefs, x0s, 
                                                         n_before_min=self.spin_n_before.value(), 
                                                         weight_func=weight_func, 
                                                         EminMin=self.spin_emin_min.value())
            fit.setWeights(weights0)
            
        # 5. Run relaxation
        nstep = self.spin_nstep.value()
        self.trjs = fit.setTrjBuffs(niter=nstep)
        fit.run(iparallel=0, ialg=1, nstep=nstep, Fmax=1e-8, dt=0.05, damping=0.0, max_step=0.1, bClamp=True)
        
        fit.getBuffs()
        
        # 6. Read back DOFs into GUI
        vals = np.array(fit.DOFs)
        for i in range(self.session.nDOFs):
            self.dof_spins[i].blockSignals(True)
            self.dof_spins[i].setValue(vals[i])
            self.dof_spins[i].blockSignals(False)
            
        print("Relaxation complete."); sys.stdout.flush()
        print("Calling evaluate and plot"); sys.stdout.flush(); self.do_evaluate_and_plot()
        
    def create_plot(self):
        self.fig = Figure(figsize=(6, 4))
        self.canvas = FigureCanvas(self.fig)
        self.layout.addWidget(self.canvas)
        
        self.axR = self.fig.add_subplot(121)
        self.axE = self.fig.add_subplot(122)
        
    def trigger_update(self):
        self.update_timer.start(150) # debounce 150ms
        
    def do_evaluate_and_plot(self):
        dof_values = [s.value() for s in self.dof_spins]
        Gm = self.session.evaluate(dof_values, self.spin_kMorse.value(), self.spin_Lepairs.value())
        
        GMS, _, _ = fit.shift_grid(Gm)
        Epanel_mod = GMS.T if GMS is not None else None
        
        to_kcal = self.chk_kcal.isChecked()
        ev2kcal = 23.0609
        
        rR, eR = fit.compute_min_lines_from_panel(self.session.Epanel_ref, self.session.Xpanel, self.session.angles)
        rM = eM = None
        if Epanel_mod is not None:
            rM, eM = fit.compute_min_lines_from_panel(Epanel_mod, self.session.Xpanel, self.session.angles)
            
        if to_kcal:
            eR = [e * ev2kcal if e is not None else None for e in eR]
            if eM is not None:
                eM = [e * ev2kcal if e is not None else None for e in eM]
        
        # Determine fixed limits from reference data
        if self.ref_e_min is None:
            eR_clean = [e for e in eR if e is not None]
            if eR_clean:
                span = max(eR_clean) - min(eR_clean)
                self.ref_e_min = min(eR_clean) - 0.2 * span
                self.ref_e_max = max(eR_clean) + 0.2 * span
            rR_clean = [r for r in rR if r is not None]
            if rR_clean:
                span = max(rR_clean) - min(rR_clean)
                self.ref_r_min = min(rR_clean) - 0.2 * span
                self.ref_r_max = max(rR_clean) + 0.2 * span
                
        self.axR.clear()
        self.axE.clear()
        
        self.axR.plot(self.session.angles, rR, '.-k', lw=1, ms=4, label='Ref')
        if rM is not None:
            self.axR.plot(self.session.angles, rM, '.-r', lw=1, ms=4, label='Model')
        self.axR.set_title('r_min(angle)')
        self.axR.set_xlabel('Angle (deg)')
        self.axR.set_ylabel('Distance x0 [Å]')
        if self.ref_r_min is not None:
            self.axR.set_ylim(self.ref_r_min, self.ref_r_max)
        self.axR.legend()
        
        self.axE.plot(self.session.angles, eR, '.-k', lw=1, ms=4, label='Ref')
        if eM is not None:
            self.axE.plot(self.session.angles, eM, '.-r', lw=1, ms=4, label='Model')
        self.axE.set_title('E_min(angle)')
        self.axE.set_xlabel('Angle (deg)')
        self.axE.set_ylabel('Energy [kcal/mol]' if to_kcal else 'Energy [eV]')
        if self.ref_e_min is not None:
            self.axE.set_ylim(self.ref_e_min, self.ref_e_max)
        self.axE.legend()
        
        self.fig.tight_layout()
        self.canvas.draw()

if __name__ == "__main__":
    xyz_path = "../tFitREQ_PN/wb97m-split/H2O-A1_H2O-D1-y.xyz"
    dof_file = "dofSelection_MorseSR_H2O.dat"
    
    app = QApplication(sys.argv)
    session = FitReqSession(xyz_path, dof_file)
    gui = FitReqGUI(session)
    gui.show()
    sys.exit(app.exec_())
