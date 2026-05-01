"""
GrapheneRibbonGUI.py - GUI for building graphene ribbons with zigzag edges
========================================================================
Interactive GUI to generate and visualize graphene nanoribbons with various passivations.
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from PyQt5 import QtWidgets, QtCore, QtGui
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

# Add FireCore path for pyBall imports
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../")))
from pyBall.OGL.BaseGUI import BaseGUI
from GrapheneRibbonBuilder import GrapheneRibbonBuilder

class GrapheneRibbonGUI(BaseGUI):
    def __init__(self):
        super().__init__(title="Graphene Ribbon Builder")
        self.builder = GrapheneRibbonBuilder()
        self.init_ui()
        self.build_and_draw()

    def init_ui(self):
        l_main = QtWidgets.QHBoxLayout(self.main_widget)
        l_ctrl = QtWidgets.QVBoxLayout()
        l_main.addLayout(l_ctrl, 1)
        
        self.fig = Figure(figsize=(12, 8))
        self.canvas = FigureCanvas(self.fig)
        l_main.addWidget(self.canvas, 5)
        
        # --- Control Panel ---
        self.ch_auto = self.checkBox("Auto-Refresh", checked=True, layout=l_ctrl)
        self.button("Build Ribbon", callback=self.build_and_draw, layout=l_ctrl)
        self.button("Save XYZ", callback=self.save_xyz, layout=l_ctrl)
        
        form = self.add_form_layout(l_ctrl)
        
        # Ribbon Parameters
        self.sp_width = self.spinBox(4, 1, 80, label="Width (chains)", layout=form, callback=self.on_param_change)
        self.sp_length = self.spinBox(8, 1, 80, label="Length (cells)", layout=form, callback=self.on_param_change)
        
        # Bond Length
        self.sp_aCC = self.spinBox(1.42, 0.1, 80, vmin=0.1, vmax=3.0, decimals=3, label="C-C Bond (A)", layout=form, callback=self.on_param_change)
        
        # Scaling Factors (anisotropic strain) - in percent
        self.sp_scale_x = self.spinBox(100.0, 80, vmin=50.0, vmax=300.0, decimals=1, label="Scale X (%)", layout=form, callback=self.on_param_change)
        self.sp_scale_x.setSingleStep(0.1)
        self.sp_scale_y = self.spinBox(100.0, 80, vmin=50.0, vmax=300.0, decimals=1, label="Scale Y (%)", layout=form, callback=self.on_param_change)
        self.sp_scale_y.setSingleStep(0.1)
        
        # Passivation
        passivation_options = ["N", "NH", "CH", "H", "O", "C=O", "C-OH"]
        self.cb_passivation = self.comboBox(passivation_options, callback=self.on_param_change, layout=l_ctrl)
        
        # Additional Options
        self.chk_start_A = self.checkBox("Start with A strip", checked=True, layout=l_ctrl, callback=self.on_param_change)
        
        # Y-offsets for passivation
        self.sp_y_top = self.spinBox(0.0, 80, vmin=-2.0, vmax=2.0, decimals=3, label="Y-Top Offset (A)", layout=form, callback=self.on_param_change)
        self.sp_y_bot = self.spinBox(0.0, 80, vmin=-2.0, vmax=2.0, decimals=3, label="Y-Bot Offset (A)", layout=form, callback=self.on_param_change)
        
        # Info labels
        self.lbl_info = QtWidgets.QLabel("Atoms: -")
        l_ctrl.addWidget(self.lbl_info)
        
        l_ctrl.addStretch()

    def on_param_change(self):
        if self.ch_auto.isChecked():
            self.build_and_draw()

    def build_and_draw(self):
        # Update builder parameters
        self.builder.a_CC = self.sp_aCC.value()
        
        # Build ribbon
        width = int(self.sp_width.value())
        length = int(self.sp_length.value())
        passivation = self.cb_passivation.currentText()
        start_with_A = self.chk_start_A.isChecked()
        y_top_offset = self.sp_y_top.value() if self.sp_y_top.value() != 0.0 else None
        y_bottom_offset = self.sp_y_bot.value() if self.sp_y_bot.value() != 0.0 else None
        scale_x = self.sp_scale_x.value() / 100.0  # Convert percent to scaling factor
        scale_y = self.sp_scale_y.value() / 100.0  # Convert percent to scaling factor
        
        pos, elems, bonds = self.builder.build_zigzag_ribbon(
            width, length, passivation,
            start_with_A=start_with_A,
            y_top_offset=y_top_offset,
            y_bottom_offset=y_bottom_offset,
            scale_x=scale_x,
            scale_y=scale_y
        )
        
        # Update info
        unique_elems, counts = np.unique(elems, return_counts=True)
        elem_str = ", ".join([f"{e}:{c}" for e, c in zip(unique_elems, counts)])
        self.lbl_info.setText(f"Atoms: {len(elems)} ({elem_str})")
        
        # Draw
        self.draw(pos, elems, bonds)

    def draw(self, pos, elems, bonds):
        self.fig.clear()
        ax = self.fig.add_subplot(111)
        
        # Color scheme for elements
        colors = {
            'C': '#333333',
            'N': '#3050F8',
            'O': '#FF0D0D',
            'H': '#FFFFFF'
        }
        
        # Plot atoms
        for elem in np.unique(elems):
            mask = elems == elem
            ax.scatter(pos[mask, 0], pos[mask, 1], 
                      c=colors.get(elem, 'gray'), 
                      s=200 if elem == 'C' else 100,
                      edgecolors='black' if elem == 'H' else None,
                      linewidth=0.5 if elem == 'H' else 0,
                      label=elem,
                      zorder=2)
        
        # Draw bonds using explicit bond list
        for bond in bonds:
            i, j = bond
            ax.plot([pos[i, 0], pos[j, 0]], 
                   [pos[i, 1], pos[j, 1]], 
                   'k-', linewidth=1.5, alpha=0.5, zorder=1)
        
        ax.set_aspect('equal')
        ax.set_xlabel('x (Å)')
        ax.set_ylabel('y (Å)')
        ax.set_title(f'Zigzag Graphene Ribbon - {self.cb_passivation.currentText()} Passivation')
        ax.legend(loc='upper right')
        ax.grid(True, alpha=0.3)
        
        self.canvas.draw()

    def save_xyz(self):
        filename, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Save XYZ", "", "XYZ Files (*.xyz);;All Files (*)"
        )
        if filename:
            self.builder.save_xyz(filename)
            print(f"Saved to {filename}")

if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    gui = GrapheneRibbonGUI()
    gui.show()
    sys.exit(app.exec_())
