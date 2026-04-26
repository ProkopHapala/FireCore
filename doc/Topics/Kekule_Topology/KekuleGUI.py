"""
KekuleGUI.py - Interactive GUI for Kekule Topology
==================================================
Usage: python KekuleGUI.py
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from PyQt5 import QtWidgets, QtCore
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

# Add FireCore path for pyBall imports
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../")))
from pyBall.OGL.BaseGUI import BaseGUI
from KekuleSolver import KekuleSolver

class KekuleGUI(BaseGUI):
    def __init__(self, solver):
        super().__init__(title="Kekule Topology Explorer (Enhanced)")
        self.solver = solver
        self.init_ui()
        self.regen_lattice()

    def init_ui(self):
        l_main = QtWidgets.QHBoxLayout(self.main_widget)
        l_ctrl = QtWidgets.QVBoxLayout()
        l_main.addLayout(l_ctrl, 1)
        
        self.fig = Figure(figsize=(14, 10))
        self.canvas = FigureCanvas(self.fig)
        l_main.addWidget(self.canvas, 5)
        
        # --- Control Panel ---
        self.ch_auto = self.checkBox("Auto-Refresh", checked=True, layout=l_ctrl)
        self.button("Calculate", callback=self.update_physics, layout=l_ctrl)
        
        form = self.add_form_layout(l_ctrl)
        
        # Physics Params
        self.sp_delta = self.spinBox(self.solver.delta, 0.05, label="Delta 1", layout=form, callback=self.on_param_change)
        self.sp_delta2 = self.spinBox(self.solver.delta2, 0.05, label="Delta 2", layout=form, callback=self.on_param_change)
        self.sp_phi1 = self.spinBox(self.solver.phi1, 0.1, label="Phase 1", layout=form, callback=self.on_param_change)
        self.sp_phi2 = self.spinBox(self.solver.phi2, 0.1, label="Phase 2", layout=form, callback=self.on_param_change)
        self.sp_onsite = self.spinBox(0.0, 0.1, label="Boundary Onsite (N)", layout=form, callback=self.on_param_change)
        
        # Energy Window for LDOS
        self.sp_emin = self.spinBox(-0.1, 0.05, label="LDOS E-Min", layout=form, callback=self.draw)
        self.sp_emax = self.spinBox(0.1, 0.05, label="LDOS E-Max", layout=form, callback=self.draw)
        
        # Lattice Params
        self.sp_Lx = self.spinBox(20.0, 1.0, label="Target Lx", layout=form, callback=self.regen_lattice)
        self.sp_Ly = self.spinBox(20.0, 1.0, label="Target Ly", layout=form, callback=self.regen_lattice)
        
        # Patch Params
        # For Step mode, Lx/2 is a good default (center)
        self.sp_Pw = self.spinBox(10.0, 1.0, label="Patch Width / X-Bound", layout=form, callback=self.on_param_change)
        self.sp_Ph = self.spinBox(5.0, 1.0, label="Patch Height", layout=form, callback=self.on_param_change)
        
        # Actual size labels
        self.lbl_actual = QtWidgets.QLabel("Actual Size: -")
        l_ctrl.addWidget(self.lbl_actual)
        
        # Periodicity
        l_pbc = QtWidgets.QHBoxLayout()
        self.ch_pbc_x = self.checkBox("PBC X", callback=self.regen_lattice, layout=l_pbc)
        self.ch_pbc_y = self.checkBox("PBC Y", callback=self.regen_lattice, layout=l_pbc)
        l_ctrl.addLayout(l_pbc)
        
        self.cb_domain = self.comboBox(["Rect", "Circle", "Step"], callback=self.on_param_change, layout=l_ctrl)
        
        self.button("Save Image", callback=self.save_fig, layout=l_ctrl)

    def on_param_change(self):
        if self.ch_auto.isChecked():
            self.update_physics()

    def regen_lattice(self):
        self.solver.Lx = self.sp_Lx.value()
        self.solver.Ly = self.sp_Ly.value()
        self.solver.pbc = (self.ch_pbc_x.isChecked(), self.ch_pbc_y.isChecked())
        self.solver.generate_lattice()
        
        # Update labels with snapped size
        self.lbl_actual.setText(f"Actual Size: {self.solver.Lx:.2f} x {self.solver.Ly:.2f}")
        
        if self.ch_auto.isChecked():
            self.update_physics()

    def update_physics(self):
        self.solver.delta = self.sp_delta.value()
        self.solver.delta2 = self.sp_delta2.value()
        self.solver.phi1 = self.sp_phi1.value()
        self.solver.phi2 = self.sp_phi2.value()
        self.solver.domain_type = self.cb_domain.currentText()
        self.solver.patch_size = [self.sp_Pw.value(), self.sp_Ph.value()]
        
        self.solver.set_onsite_at_boundary(width=1.0, energy=self.sp_onsite.value())
        self.solver.build_hamiltonian()
        self.solver.solve()
        self.draw()

    def draw(self):
        self.fig.clear()
        # 2x2 Grid
        ax_lat = self.fig.add_subplot(221)
        ax_phi = self.fig.add_subplot(222)
        ax_spec = self.fig.add_subplot(223)
        ax_dos = self.fig.add_subplot(224)
        
        # Color Scheme
        from matplotlib.colors import LinearSegmentedColormap
        # Bone-like: light yellow/gray for 0, Deep blue for high
        custom_cmap = LinearSegmentedColormap.from_list("kekule", ["#f0f0d0", "#000080"])
        
        # 1. Bonds & LDOS
        ax_lat.set_aspect('equal')
        ax_lat.set_facecolor('#ffffff')
        t_vals = [b[2] for b in self.solver.bonds]
        if t_vals:
            t_min, t_max = min(t_vals), max(t_vals)
            for i, j, t, r_mid, k, s in self.solver.bonds:
                p1 = self.solver.pos[i]
                p2 = self.solver.pos[j] + s
                lw = 0.5 + 4.0 * (t - t_min) / (t_max - t_min + 1e-9)
                ax_lat.plot([p1[0], p2[0]], [p1[1], p2[1]], color='gray', linewidth=lw, alpha=0.3, zorder=1)
                if np.any(s):
                    p1_alt = p1 - s
                    p2_alt = self.solver.pos[j]
                    ax_lat.plot([p1_alt[0], p2_alt[0]], [p1_alt[1], p2_alt[1]], color='gray', linewidth=lw, alpha=0.3, zorder=1)
        
        ax_lat.set_xlim(-self.solver.Lx/2, self.solver.Lx/2)
        ax_lat.set_ylim(-self.solver.Ly/2, self.solver.Ly/2)
        
        ldos = self.solver.get_ldos(self.sp_emin.value(), self.sp_emax.value())
        if np.any(ldos):
            # Use vmax to ensure even small densities are visible if they are the max in the window
            vmax = np.max(ldos) if np.max(ldos) > 1e-6 else 1.0
            sc = ax_lat.scatter(self.solver.pos[:,0], self.solver.pos[:,1], c=ldos, cmap=custom_cmap, s=100, 
                                edgecolors='k', linewidth=0.3, zorder=2, vmin=0, vmax=vmax)
            self.fig.colorbar(sc, ax=ax_lat, label="LDOS")
        ax_lat.set_title(f"Localized States (LDOS)")
        
        # 2. Phase Map
        ax_phi.set_aspect('equal')
        gx, gy = np.meshgrid(np.linspace(-self.solver.Lx/2, self.solver.Lx/2, 25), 
                             np.linspace(-self.solver.Ly/2, self.solver.Ly/2, 25))
        gphi = np.zeros_like(gx)
        for i in range(gx.shape[0]):
            for j in range(gx.shape[1]):
                gphi[i,j], _ = self.solver.get_phi_delta([gx[i,j], gy[i,j]])
        ax_phi.contourf(gx, gy, gphi % (2*np.pi), levels=50, cmap='hsv', alpha=0.3)
        ax_phi.quiver(gx, gy, np.cos(gphi), np.sin(gphi), color='black', alpha=0.3)
        ax_phi.set_title("Kekule Phase Map")
        
        # 3. Spectrum
        ax_spec.plot(self.solver.evals, 'o', markersize=3, color='#000080', alpha=0.6)
        ax_spec.axhline(0, color='red', linestyle='--', alpha=0.3)
        ax_spec.axhspan(self.sp_emin.value(), self.sp_emax.value(), color='yellow', alpha=0.3)
        ax_spec.set_ylabel("Energy (eV)")
        ax_spec.set_xlabel("State Index")
        ax_spec.set_title("Energy Spectrum")
        
        # 4. DOS
        # Use a dynamic range for DOS focused on the area of interest
        e_min_dos = min(-1.0, self.sp_emin.value() - 0.5)
        e_max_dos = max( 1.0, self.sp_emax.value() + 0.5)
        e_range = np.linspace(e_min_dos, e_max_dos, 500)
        dos = self.solver.get_dos(e_range, sigma=0.02)
        ax_dos.fill_between(e_range, dos, color='#000080', alpha=0.2)
        ax_dos.plot(e_range, dos, color='#000080', linewidth=1.5)
        ax_dos.axvline(0, color='red', linestyle='--', alpha=0.3)
        ax_dos.axvspan(self.sp_emin.value(), self.sp_emax.value(), color='yellow', alpha=0.2)
        ax_dos.set_xlabel("Energy (eV)")
        ax_dos.set_title("Density of States (DOS)")
        ax_dos.set_ylim(0, None)
        ax_dos.set_xlim(e_min_dos, e_max_dos)
        
        self.canvas.draw()

    def save_fig(self):
        self.fig.savefig("Kekule_Enhanced_Output.png")
        print("Saved to Kekule_Enhanced_Output.png")

if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    solver = KekuleSolver()
    gui = KekuleGUI(solver)
    gui.show()
    sys.exit(app.exec_())
