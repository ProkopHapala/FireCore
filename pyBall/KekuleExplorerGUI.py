#!/usr/bin/env python3
"""
CODE STYLE POLICIES:
- Strive for concise, general, and reusable code
- Modularity and composability over duplication
- Minimize code duplication
- Prefer single-line function calls and messages
- Extract repeated logic to shared utilities (e.g., VispyUtils.py)
- Use BaseGUI helper methods for widget creation to reduce boilerplate
- Use polymorphic functions with default arguments instead of specialized variants
- Consolidate similar functions (e.g., spinBox/spinBoxInt → spinBox(int_mode=True))
- Refactor if-else labyrinths into general functions with callbacks
"""

import sys
import os
import numpy as np
from PyQt5 import QtWidgets, QtCore, QtGui
from vispy import scene

# Add FireCore to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from pyBall.KekuleBackend import KekuleBackend, parse_passivation_string
from pyBall.AtomicSystem import AtomicSystem
from pyBall import VispyUtils as vu
from pyBall import atomicUtils as au
from pyBall import elements
from pyBall.GUI.BaseGUI import BaseGUI
from pyBall.VispyUtils import compute_bond_colors_by_length, generate_atom_labels

try:
    from pyBall import FireCore as fc
except Exception as e:
    print(f"ERROR: FireCore not available: {e}")
    raise

class KekuleExplorerWindow(BaseGUI):
    def __init__(self):
        super().__init__("Kekule Structure Explorer")
        self.resize(1024, 768)

        self.backend = KekuleBackend()
        self.cur_atom_type = 'C'
        self.edit_mode = 'Hex1'  # 'Hex1' (paint), 'Hex2' (toggle), 'Atom', 'pi', 'Select'
        self.last_clicked_node = None
        self.label_mode = 'Element+Index'

        # Load settings
        self.settings = QtCore.QSettings("FireCore", "KekuleExplorer")
        self.fdata_path = self.settings.value("fdata_path", "/home/prokop/Fireball/Fdata_HCNOS")
        
        self.initUI()
        # Scene now operates directly on backend.sys.apos, no drag signal needed
        self.scene.sig_rmb_remove.connect(self.on_atom_remove)
        self.scene.sig_camera_changed.connect(self.refresh_view)
        self.refresh_view()

    def initUI(self):
        # --- Central Widget (Vispy Scene) ---
        self.scene = vu.AtomScene(bgcolor=(0.95, 0.95, 0.95), backend=self.backend)
        
        # Link axes to view
        self.scene.view.parent = None # Re-parent from central_widget to grid
        grid = self.scene.canvas.central_widget.add_grid(spacing=0, margin=10)
        
        self.axis_x = scene.AxisWidget(orientation='bottom', axis_label='x (A)', font_size=8)
        self.axis_y = scene.AxisWidget(orientation='left', axis_label='y (A)', font_size=8)
        
        self.axis_x.height_max = 30
        self.axis_y.width_max = 40

        grid.add_widget(self.axis_y, row=0, col=0)
        grid.add_widget(self.scene.view, row=0, col=1)
        grid.add_widget(self.axis_x, row=1, col=1)
        
        self.scene.view.stretch = (1, 1)
        
        self.axis_x.link_view(self.scene.view)
        self.axis_y.link_view(self.scene.view)
        
        # Configure axis after linking to view
        self.axis_x.axis.text_color = 'black'
        self.axis_y.axis.text_color = 'black'
        self.axis_x.axis.tick_color = 'black'
        self.axis_y.axis.tick_color = 'black'

        # --- Main Layout with Side Panel ---
        main_widget = QtWidgets.QWidget()
        main_layout = QtWidgets.QHBoxLayout(main_widget)
        
        # Create side panel
        side_panel = QtWidgets.QFrame()
        side_panel.setFrameStyle(QtWidgets.QFrame.StyledPanel)
        side_panel.setFixedWidth(280)
        side_layout = QtWidgets.QVBoxLayout(side_panel)
        
        # Add sections
        side_layout.addWidget(self.create_ribbon_section())
        side_layout.addWidget(self.create_builder_section())
        side_layout.addWidget(self.create_editor_section())
        side_layout.addWidget(self.create_dftb_section())
        side_layout.addWidget(self.create_fireball_section())
        side_layout.addStretch()
        
        # Add to main layout
        main_layout.addWidget(side_panel)
        main_layout.addWidget(self.scene.canvas.native)
        
        self.setCentralWidget(main_widget)

        # Add grid guide markers
        self.grid_markers = scene.visuals.Markers(parent=self.scene.view.scene)
        self.grid_markers.set_gl_state('translucent', depth_test=False)
        self.grid_markers.order = -1  # Behind everything

        # Add mouse cursor (cross) for debugging
        self.cursor_markers = scene.visuals.Markers(parent=self.scene.view.scene)
        self.cursor_markers.set_gl_state('translucent', depth_test=False)
        self.cursor_markers.order = 10  # On top
        self.cursor_markers.set_data(
            pos=np.zeros((1, 3)),
            symbol='cross',
            edge_width=2,
            edge_color='red',
            face_color='transparent',
            size=10
        )

        # Add debug markers for grid node -> atom mappings (cyan)
        self.debug_markers = scene.visuals.Markers(parent=self.scene.view.scene)
        self.debug_markers.set_gl_state('translucent', depth_test=False)
        self.debug_markers.order = 5  # Behind atoms but in front of grid

        # Add hover markers for highlighting hexagon under mouse
        self.hover_markers = scene.visuals.Markers(parent=self.scene.view.scene)
        self.hover_markers.set_gl_state('translucent', depth_test=False)
        self.hover_markers.order = 6  # On top of debug markers

        # Add debug lines for node -> atom connections
        self.debug_lines = scene.visuals.Line(parent=self.scene.view.scene)
        self.debug_lines.set_gl_state('translucent', depth_test=False)
        self.debug_lines.order = 4  # Behind debug markers

        # Help / Status
        self.statusBar().showMessage("LMB: Add/Toggle | RMB: Remove | Middle-Click: Toggle H | Scroll: Zoom")
        self.scene.lock_drag = True   # Default mode is Ring, no dragging
        self.scene.canvas.events.mouse_press.connect(self.on_mouse_press)
        self.scene.canvas.events.mouse_move.connect(self.on_mouse_move)
        self.scene.sig_selection_changed.connect(self.on_selection_changed)
        self.copied_atoms = None  # (enames, apos) tuple
        self.scene.canvas.events.key_press.connect(self.on_key_press)
        self.create_menus()
        self.error_print = True      # Print to stdout
        self.error_raise = True      # Raise exception
        self.error_dialog = True     # Show QMessageBox
        self.error_statusbar = True  # Update status bar

    def _raise(self, msg, title="Error", dialog_type="critical"):
        """Reusable error handling function.
        
        Args:
            msg: Error message
            title: Dialog title
            dialog_type: 'critical', 'warning', or 'information'
        """
        if self.error_print:
            print(msg)
        if self.error_statusbar:
            self.statusBar().showMessage(msg)
        if self.error_dialog:
            if dialog_type == "critical":
                QtWidgets.QMessageBox.critical(self, title, str(msg))
            elif dialog_type == "warning":
                QtWidgets.QMessageBox.warning(self, title, str(msg))
            elif dialog_type == "information":
                QtWidgets.QMessageBox.information(self, title, str(msg))
        if self.error_raise:
            raise RuntimeError(msg)

    def create_builder_section(self):
        layout = QtWidgets.QVBoxLayout()
        self.label("Edit Mode:", layout=layout)
        self.mode_combo = self.comboBox(["Hex1", "Hex2", "Atom", "pi", "Select"], self.set_edit_mode, layout=layout)
        self.label("Atom Type:", layout=layout)
        self.atom_combo = self.comboBox(["C", "N", "O"], self.set_atom_type, layout=layout)
        return self.group("Builder", layout)

    def create_editor_section(self):
        layout = QtWidgets.QVBoxLayout()
        self.button("Snap to Grid", self.reset_offsets, layout=layout)
        self.button("Adjust H", self.adjust_h, layout=layout)
        self.button("Recalc Bonds", self.recalc_bonds, layout=layout)
        layout.addSpacing(10)
        self.label("Labels:", layout=layout)
        self.label_combo = self.comboBox(["Element+Index", "Atomic Type", "Pi Orbitals", "Z-Height", "Charge", "Bond Lengths"], self.set_label_mode, layout=layout)
        self.bond_viz_mode = False
        self.button("Bond Colors", self.toggle_bond_viz, layout=layout).setCheckable(True)
        layout.addSpacing(10)
        self.debug_view_mode = True
        self.debug_btn = self.button("Debug View", self.toggle_debug_view, layout=layout)
        self.debug_btn.setCheckable(True)
        self.debug_btn.setChecked(True)
        layout.addSpacing(10)
        self.button("Show XYZ", self.show_xyz, layout=layout)
        self.button("Export XYZ", self.export_xyz, layout=layout)
        return self.group("Editor", layout)

    def create_dftb_section(self):
        layout = QtWidgets.QVBoxLayout()
        self.button("Relax (DFTB+)", self.run_relaxation, layout=layout)
        self.dftb_status_label = self.label("Status: Ready", layout=layout, word_wrap=True)
        return self.group("DFTB+", layout)

    def create_fireball_section(self):
        layout = QtWidgets.QVBoxLayout()
        self.button("Compute SCF", self.compute_orbitals, layout=layout)
        self.orbital_info_label = self.label("Orbitals: Not computed", layout=layout, word_wrap=True)
        layout.addSpacing(10)
        self.label("Z-height (Å):", layout=layout)
        self.z_height_spinbox = self.spinBox(2.0, 0.5, layout=layout, vmin=-10.0, vmax=20.0)
        layout.addSpacing(10)
        self.label("Orbital Index:", layout=layout)
        self.orbital_spinbox = self.spinBox(0, vmin=0, vmax=999, enabled=False, callback=self.update_orbital_energy_label, layout=layout, int_mode=True)
        layout.addSpacing(10)
        self.plot_orb_btn = self.button("Plot Orbital", self.plot_orbital_from_spinbox, layout=layout)
        self.plot_orb_btn.setEnabled(False)
        layout.addSpacing(10)
        self.plot_density_btn = self.button("Plot Density", self.plot_density, layout=layout)
        self.plot_density_btn.setEnabled(False)
        layout.addSpacing(10)
        self.plot_delta_btn = self.button("Plot Delta-Rho", self.plot_delta_rho, layout=layout)
        self.plot_delta_btn.setEnabled(False)
        layout.addSpacing(10)
        self.button("Set Fdata Path", self.set_fdata_path, layout=layout)
        return self.group("Fireball", layout)

    def create_ribbon_section(self):
        layout = QtWidgets.QVBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(2)
        
        # Shared inputs (used by both single and two-ribbon)
        shared_layout = QtWidgets.QHBoxLayout()
        shared_layout.setContentsMargins(0, 0, 0, 0)
        shared_layout.setSpacing(2)
        
        rows_label = QtWidgets.QLabel("Rows:")
        rows_label.setFixedWidth(30)
        self.ribbon_rows_spinbox = self.spinBox(4, 1.0, max_width=50, vmin=1, vmax=20, int_mode=True)
        shared_layout.addWidget(rows_label)
        shared_layout.addWidget(self.ribbon_rows_spinbox)
        
        bottom_label = QtWidgets.QLabel("Bot:")
        bottom_label.setFixedWidth(25)
        self.ribbon_bottom_edit = QtWidgets.QLineEdit()
        self.ribbon_bottom_edit.setPlaceholderText("n/N/o/O/H/h")
        self.ribbon_bottom_edit.setMaximumWidth(80)
        shared_layout.addWidget(bottom_label)
        shared_layout.addWidget(self.ribbon_bottom_edit)
        
        top_label = QtWidgets.QLabel("Top:")
        top_label.setFixedWidth(25)
        self.ribbon_top_edit = QtWidgets.QLineEdit()
        self.ribbon_top_edit.setPlaceholderText("n/N/o/O/H/h")
        self.ribbon_top_edit.setMaximumWidth(80)
        shared_layout.addWidget(top_label)
        shared_layout.addWidget(self.ribbon_top_edit)
        
        layout.addLayout(shared_layout)
        
        # Generate buttons side by side
        btn_layout = QtWidgets.QHBoxLayout()
        btn_layout.setContentsMargins(0, 0, 0, 0)
        btn_layout.setSpacing(2)
        self.button("Single", self.generate_single_ribbon, layout=btn_layout)
        self.button("Two", self.generate_two_ribbons, layout=btn_layout)
        layout.addLayout(btn_layout)
        
        # Two-ribbon specific inputs (collapsible)
        two_ribbon_group = QtWidgets.QGroupBox("Two-Ribbon Options")
        two_ribbon_group.setCheckable(True)
        two_ribbon_group.setChecked(False)
        two_ribbon_group.toggled.connect(lambda checked: two_ribbon_group.setVisible(checked))
        two_ribbon_layout = QtWidgets.QVBoxLayout()
        two_ribbon_layout.setContentsMargins(0, 0, 0, 0)
        two_ribbon_layout.setSpacing(2)
        
        # Ribbon 2 inputs
        r2_layout = QtWidgets.QHBoxLayout()
        r2_layout.setContentsMargins(0, 0, 0, 0)
        r2_layout.setSpacing(2)
        
        r2_rows_label = QtWidgets.QLabel("R2:")
        r2_rows_label.setFixedWidth(20)
        self.ribbon2_rows_spinbox = self.spinBox(4, 1.0, max_width=50, vmin=1, vmax=20, int_mode=True)
        r2_layout.addWidget(r2_rows_label)
        r2_layout.addWidget(self.ribbon2_rows_spinbox)
        
        r2_bottom_label = QtWidgets.QLabel("Bot:")
        r2_bottom_label.setFixedWidth(25)
        self.ribbon2_bottom_edit = QtWidgets.QLineEdit()
        self.ribbon2_bottom_edit.setPlaceholderText("n/N/o/O/H/h")
        self.ribbon2_bottom_edit.setMaximumWidth(80)
        r2_layout.addWidget(r2_bottom_label)
        r2_layout.addWidget(self.ribbon2_bottom_edit)
        
        r2_top_label = QtWidgets.QLabel("Top:")
        r2_top_label.setFixedWidth(25)
        self.ribbon2_top_edit = QtWidgets.QLineEdit()
        self.ribbon2_top_edit.setPlaceholderText("n/N/o/O/H/h")
        self.ribbon2_top_edit.setMaximumWidth(80)
        r2_layout.addWidget(r2_top_label)
        r2_layout.addWidget(self.ribbon2_top_edit)
        
        two_ribbon_layout.addLayout(r2_layout)
        
        # H-bond spacing
        hb_layout = QtWidgets.QHBoxLayout()
        hb_layout.setContentsMargins(0, 0, 0, 0)
        hb_layout.setSpacing(2)
        hb_label = QtWidgets.QLabel("H-bond:")
        hb_label.setFixedWidth(50)
        self.ribbon_L_Hb_spinbox = self.spinBox(3.0, 0.1, max_width=60, vmin=2.0, vmax=10.0)
        hb_layout.addWidget(hb_label)
        hb_layout.addWidget(self.ribbon_L_Hb_spinbox)
        two_ribbon_layout.addLayout(hb_layout)
        
        two_ribbon_group.setLayout(two_ribbon_layout)
        layout.addWidget(two_ribbon_group)
        
        return self.group("Ribbon", layout)

    def generate_single_ribbon(self):
        """Generate single periodic ribbon from passivation strings."""
        bottom_str = self.ribbon_bottom_edit.text().strip()
        top_str = self.ribbon_top_edit.text().strip()
        
        if not bottom_str or not top_str:
            QtWidgets.QMessageBox.warning(self, "Warning", "Please provide passivation strings for both bottom and top edges.")
            return
        
        try:
            bottom_passivation = parse_passivation_string(bottom_str)
            top_passivation = parse_passivation_string(top_str)
            length_cells = len(bottom_passivation)
            width_chains = self.ribbon_rows_spinbox.value()
            Lx = 2.4
            
            self.backend = KekuleBackend()
            self.backend.build_zigzag_ribbon(width_chains=width_chains, length_cells=length_cells, passivation_bottom=bottom_passivation, passivation_top=top_passivation, scale_x=Lx / (2.0 * 1.42 * np.cos(np.pi / 6)), bPeriodicX=True)
            
            self.backend.recalc_bonds()
            self.scene.backend = self.backend
            
            n_C = sum(1 for e in self.backend.sys.enames if e == 'C')
            n_N = sum(1 for e in self.backend.sys.enames if e == 'N')
            n_O = sum(1 for e in self.backend.sys.enames if e == 'O')
            n_H = sum(1 for e in self.backend.sys.enames if e == 'H')
            
            msg = f"Generated single ribbon: C={n_C}, N={n_N}, O={n_O}, H={n_H}"
            print(msg)
            self.statusBar().showMessage(msg)
            self.refresh_view()
            
        except Exception as e:
            self._raise(f"Ribbon generation FAILED: {e}", title="Ribbon Error")

    def generate_two_ribbons(self):
        """Generate two-ribbon system from passivation strings."""
        bottom1_str = self.ribbon_bottom_edit.text().strip()
        top1_str = self.ribbon_top_edit.text().strip()
        bottom2_str = self.ribbon2_bottom_edit.text().strip()
        top2_str = self.ribbon2_top_edit.text().strip()
        
        if not all([bottom1_str, top1_str, bottom2_str, top2_str]):
            QtWidgets.QMessageBox.warning(self, "Warning", "Please provide passivation strings for all four edges.")
            return
        
        try:
            bottom1_passivation = parse_passivation_string(bottom1_str)
            top1_passivation = parse_passivation_string(top1_str)
            bottom2_passivation = parse_passivation_string(bottom2_str)
            top2_passivation = parse_passivation_string(top2_str)
            length_cells = len(bottom1_passivation)
            width_chains1 = self.ribbon_rows_spinbox.value()
            width_chains2 = self.ribbon2_rows_spinbox.value()
            Lx = 2.4
            L_Hb = self.ribbon_L_Hb_spinbox.value()
            
            # Build bottom ribbon
            bottom_ribbon = KekuleBackend()
            bottom_ribbon.build_zigzag_ribbon(width_chains=width_chains1, length_cells=length_cells, passivation_bottom=bottom1_passivation, passivation_top=top1_passivation, scale_x=Lx / (2.0 * 1.42 * np.cos(np.pi / 6)), bPeriodicX=True)
            
            # Build top ribbon
            top_ribbon = KekuleBackend()
            top_ribbon.build_zigzag_ribbon(width_chains=width_chains2, length_cells=length_cells, passivation_bottom=bottom2_passivation, passivation_top=top2_passivation,  scale_x=Lx / (2.0 * 1.42 * np.cos(np.pi / 6)), bPeriodicX=True)
            
            # Combine ribbons
            self.backend = KekuleBackend()
            self.backend.combine_ribbons(bottom_ribbon, top_ribbon, L_Hb=L_Hb, shift_x=0.0)
            
            self.backend.recalc_bonds()
            self.scene.backend = self.backend
            
            n_C = sum(1 for e in self.backend.sys.enames if e == 'C')
            n_N = sum(1 for e in self.backend.sys.enames if e == 'N')
            n_O = sum(1 for e in self.backend.sys.enames if e == 'O')
            n_H = sum(1 for e in self.backend.sys.enames if e == 'H')
            
            msg = f"Generated two-ribbon system: C={n_C}, N={n_N}, O={n_O}, H={n_H}"
            print(msg)
            self.statusBar().showMessage(msg)
            self.refresh_view()
            
        except Exception as e:
            self._raise(f"Two-ribbon generation FAILED: {e}", title="Ribbon Error")

    def create_menus(self):
        # Settings Menu
        self.settings_menu = self.menuBar().addMenu("Settings")
        # Fdata path is now in Fireball section, but keep menu for convenience
        self.settings_menu.addAction("Set Fdata Path", self.set_fdata_path)

    def update_orbital_energy_label(self, value):
        """Update orbital energy label when spinbox value changes."""
        if hasattr(self, '_eigen'):
            idx = int(value)
            if 0 <= idx < len(self._eigen):
                self.orbital_info_label.setText(f"Orbital {idx} E: {self._eigen[idx]:.3f} eV")
            else:
                self.orbital_info_label.setText("Orbital: Invalid index")

    def plot_orbital_from_spinbox(self):
        """Plot orbital at the index selected in spinbox."""
        if not hasattr(self, '_eigen'):
            msg = "Please run Compute SCF first."
            print(f"INFO: {msg}")
            QtWidgets.QMessageBox.information(self, "Info", msg)
            return
        mo_idx = self.orbital_spinbox.value()
        if mo_idx < 0 or mo_idx >= len(self._eigen):
            msg = f"Invalid orbital index: {mo_idx}"
            print(f"WARNING: {msg}")
            QtWidgets.QMessageBox.warning(self, "Warning", msg)
            return
        z_height = self.z_height_spinbox.value()
        self.statusBar().showMessage(f"Projecting MO {mo_idx + 1}...")
        QtWidgets.QApplication.processEvents()
        try:
            apos = self.backend.sys.apos
            grid_origin, size, center_z = self._compute_extent_from_geometry(apos)
            points, extent, n = self._make_2d_grid(grid_origin, size, center_z, z_height)
            flat_data = self._evaluate_on_grid(points, 'orbital', orb_index=mo_idx)
            data_2d = np.asarray(flat_data, dtype=np.float64).reshape(n, n)
            E = self._eigen[mo_idx]
            pos = apos.astype(np.float32)
            enames = self.backend.sys.enames
            self._plot_2d_projection(data_2d, extent, title=f"MO {mo_idx + 1} E={E:+.3f} eV  z={z_height:.1f}Å", cmap='bwr', symmetric=True, atom_pos=pos, atom_types=enames)
            msg = f"Plotted MO {mo_idx + 1}"
            print(msg)
            self.statusBar().showMessage(msg)
        except Exception as e:
            self._raise(f"Plot FAILED: {e}", title="Plot Error")

    def set_edit_mode(self, mode):
        self.edit_mode = mode
        print(f"Edit Mode: {mode}")
        # Sync backend hex_mode when switching Hex1/Hex2
        if mode == 'Hex1':
            self.backend.hex_mode = 'Hex1'
        elif mode == 'Hex2':
            self.backend.hex_mode = 'Hex2'
        # Auto-switch label mode to Pi Orbitals when in pi mode
        if mode == 'pi':
            self.set_label_mode('Pi Orbitals')
        # Scene and UI settings
        if mode == 'Select':
            self.scene.set_selection_mode(True)
            self.scene.lock_drag = False
            self.statusBar().showMessage("Selection Mode: RMB drag to select | Delete: Remove | Ctrl-C: Copy | Ctrl-V: Paste | LMB: Drag selected")
        elif mode in ('Hex1', 'Hex2'):
            self.scene.set_selection_mode(False)
            self.scene.lock_drag = True   # No atom dragging in hex mode
            mode_str = "Hex1 (paint: force add/remove)" if mode == 'Hex1' else "Hex2 (toggle: preserve shared)"
            self.statusBar().showMessage(f"{mode_str}: LMB: Add | RMB: Remove")
        else:
            self.scene.set_selection_mode(False)
            self.scene.lock_drag = False
            self.statusBar().showMessage("LMB: Add/Toggle | RMB: Remove | Middle-Click: Toggle H | Scroll: Zoom")

    def set_atom_type(self, atype):
        self.cur_atom_type = atype
        print(f"Atom Type: {atype}")

    def set_label_mode(self, mode):
        """Set label display mode."""
        self.label_mode = mode
        # Update combo box to reflect current mode
        index = self.label_combo.findText(mode)
        if index >= 0:
            self.label_combo.blockSignals(True)
            self.label_combo.setCurrentIndex(index)
            self.label_combo.blockSignals(False)
        self.refresh_view()

    def toggle_bond_viz(self):
        """Toggle bond color visualization mode."""
        self.bond_viz_mode = not self.bond_viz_mode
        self.refresh_view()

    def toggle_debug_view(self):
        """Toggle debug visualization mode."""
        self.debug_view_mode = not self.debug_view_mode
        self.refresh_view()

    def on_selection_changed(self, selected_indices):
        """Handle selection change from Vispy scene."""
        n_selected = len(selected_indices)
        if n_selected > 0:
            self.statusBar().showMessage(f"Selected {n_selected} atoms | Delete: Remove | Ctrl-C: Copy | Ctrl-V: Paste | LMB: Drag selected")
        elif self.edit_mode == 'Select':
            self.statusBar().showMessage("Selection Mode: RMB drag to select | Delete: Remove | Ctrl-C: Copy | Ctrl-V: Paste | LMB: Drag selected")

    def on_key_press(self, event):
        """Handle keyboard shortcuts."""
        selected = self.scene.get_selected_indices()
        if not selected:
            return

        # Check for Control modifier (Vispy uses tuple of strings)
        ctrl_pressed = 'Control' in event.modifiers if isinstance(event.modifiers, (tuple, list)) else False

        # Delete key - remove selected atoms
        if event.key == 'Delete':
            self.delete_selected_atoms()
        # Ctrl-C - copy selected atoms
        elif event.key == 'C' and ctrl_pressed:
            self.copy_selected_atoms()
        # Ctrl-V - paste copied atoms
        elif event.key == 'V' and ctrl_pressed:
            self.paste_copied_atoms()

    def delete_selected_atoms(self):
        """Delete currently selected atoms."""
        selected = list(self.scene.get_selected_indices())
        if not selected:
            return
        # Sort in descending order to avoid index issues
        selected.sort(reverse=True)
        for idx in selected:
            self.backend._rebuild_after_delete([idx])
        self.backend.recalc_bonds()
        self.scene.clear_selection()
        self.refresh_view()
        print(f"Deleted {len(selected)} atoms")

    def copy_selected_atoms(self):
        """Copy currently selected atoms to clipboard."""
        selected = list(self.scene.get_selected_indices())
        if not selected:
            return
        enames = [self.backend.sys.enames[i] for i in selected]
        apos = [self.backend.sys.apos[i].copy() for i in selected]
        self.copied_atoms = (enames, apos)
        print(f"Copied {len(selected)} atoms")

    def paste_copied_atoms(self):
        """Paste copied atoms at original position (duplicate in place)."""
        if self.copied_atoms is None:
            print("No atoms copied")
            return
        enames, apos_orig = self.copied_atoms
        # Track indices of newly added atoms
        new_indices = []
        # Add atoms at original positions using _append_atom to avoid adjust_h()
        for ename, pos in zip(enames, apos_orig):
            idx = self.backend._append_atom(pos=list(pos.copy()), ename=ename, pin=None,parent=None,subtype=f"{ename}_sp2"  )
            new_indices.append(idx)
        self.backend.recalc_bonds()
        # Don't call adjust_h() - it adds H atoms and shifts indices
        # Refresh view first to update scene arrays
        self.refresh_view()
        # Select the newly pasted atoms
        self.scene.set_selected_indices(new_indices)
        print(f"Pasted {len(enames)} atoms at original positions")

    def on_mouse_move(self, event):
        """Update cursor cross position on mouse move and highlight hexagon under cursor."""
        r0, rd = self.scene._ray_from_mouse(event.pos)
        p_world = self.scene._intersect_ray_plane(r0, rd, np.zeros(3), np.array([0,0,1]))
        if p_world is not None:
            self.cursor_markers.set_data(pos=np.array([p_world]),symbol='cross',edge_width=2,edge_color='red',face_color='transparent',size=10 )

            # Highlight hexagon under mouse if in hex mode
            if self.edit_mode in ('Hex1', 'Hex2') and hasattr(self.backend, 'snap_to_ring'):
                from pyBall.KekuleBackend import honeycomb_ring_nodes, snap_to_grid
                q, r = self.backend.snap_to_ring(p_world[0], p_world[1])
                ring_nodes = honeycomb_ring_nodes(q, r, self.backend.a_CC)
                hover_pos = []
                for node in ring_nodes:
                    nk = snap_to_grid(node, self.backend.a_CC)
                    hover_pos.append([nk[0], nk[1], -0.08])
                if hover_pos:
                    self.hover_markers.set_data(
                        pos=np.array(hover_pos, dtype=np.float32),
                        symbol='disc', edge_width=2, edge_color='orange', face_color='transparent', size=12
                    )
                    self.hover_markers.visible = True
                else:
                    self.hover_markers.visible = False
            else:
                self.hover_markers.visible = False

    def on_atom_remove(self, idx):
        """Remove atom at index and refresh view. Only in non-Hex modes."""
        if self.edit_mode in ('Hex1', 'Hex2'):
            return   # In Hex modes, RMB on atoms is ignored (hex removal handled by handle_click)
        # Find which node_key corresponds to this atom index
        node_to_remove = None
        if idx < len(self.backend.atom_pin):
            node_to_remove = self.backend.atom_pin[idx]
        if node_to_remove:
            self.backend.remove_atom(node_to_remove)
        else:
            # Atom not pinned to grid (e.g. H) - remove by object reference
            atom_list, *_ = self.backend.graph.to_arrays()
            if 0 <= idx < len(atom_list):
                atom = atom_list[idx]
                self.backend.graph.remove_atom(atom)
                self.backend.recalc_bonds()
        self.refresh_view()

    def on_mouse_press(self, event):
        # In Select mode, let Vispy handle everything (RMB selection, LMB drag)
        if self.edit_mode == 'Select':
            return

        # In Hex modes, prevent dragging - we only want add/remove hex operations
        if self.edit_mode in ('Hex1', 'Hex2'):
            # Continue to handle_click for hex operations, but don't let Vispy handle drag
            pass
        else:
            # For non-ring modes, if atoms are selected and LMB, let Vispy handle dragging
            selected = self.scene.get_selected_indices()
            if selected and event.button == 1:
                return

        # If atom picked and LMB in atom/pi mode, handle atom change instead of drag
        picked = self.scene._pick_idx
        if picked >= 0 and event.button == 1:
            if self.edit_mode == 'Atom':
                # Change atom type at picked position
                self.backend.sys.enames[picked] = self.cur_atom_type
                self.backend.sys.atypes[picked] = elements.ELEMENT_DICT[self.cur_atom_type][0]
                self.backend.atom_subtype[picked] = self.backend._get_element_default_subtype(self.cur_atom_type)
                self.backend.recalc_bonds()
                self.backend.adjust_h()
                self.refresh_view()
                return
            elif self.edit_mode == 'pi':
                # Cycle pi orbitals on picked atom
                subtype = self.backend.atom_subtype[picked]
                current_npi = self.backend._get_npi_from_subtype(subtype)
                new_npi = (current_npi + 1) % 3
                e = self.backend.sys.enames[picked]
                sp_map = {0: 'sp3', 1: 'sp2', 2: 'sp'}
                self.backend.atom_subtype[picked] = f"{e}_{sp_map.get(new_npi, 'sp2')}"
                self.backend.adjust_h()
                self.refresh_view()
                return
            elif self.edit_mode == 'Ring':
                # In ring mode, don't allow dragging atoms - ignore atom picks
                return

        if event.button == 1: # LMB
            self.handle_click(event.pos, action='add')
        elif event.button == 2: # RMB
            self.handle_click(event.pos, action='remove')
        elif event.button == 3: # Middle / Scroll click
            self.handle_click(event.pos, action='toggle_h')

    def reset_offsets(self):
        self.backend.snap_atoms_to_grid()
        self.refresh_view()

    def show_xyz(self):
        """Show current structure in a text dialog."""
        xyz_str = self.backend.get_xyz_string()
        dialog = QtWidgets.QDialog(self)
        dialog.setWindowTitle("Current XYZ Structure")
        layout = QtWidgets.QVBoxLayout(dialog)
        self.textEdit(xyz_str, read_only=True, min_size=(400, 500), layout=layout, plain=True)
        self.button("Close", dialog.accept, layout=layout)
        dialog.exec_()

    def export_xyz(self):
        """Export current structure to an XYZ file."""
        fname = self.fileDialog(mode="save", title="Export XYZ", filter_str="XYZ Files (*.xyz)")
        if fname:
            self.backend.save_xyz(fname)
            self.statusBar().showMessage(f"Exported to {fname}")

    def adjust_h(self):
        """Manually trigger H passivation."""
        self.backend.adjust_h()
        self.refresh_view()

    def recalc_bonds(self):
        """Manually trigger bond recalculation and refresh view."""
        self.backend.recalc_bonds()
        self.refresh_view()

    def handle_click(self, mouse_pos, action='add'):
        # 1. Get world coordinates on z=0 plane
        # Vispy mouse pos is (x, y) from top-left
        # We need to use Vispy's internal ray casting
        # ray = self.scene._ray_from_mouse(mouse_pos)
        # intersect = self.scene._intersect_ray_plane(ray[0], ray[1], np.zeros(3), np.array([0,0,1]))
        
        # Helper to get world pos
        r0, rd = self.scene._ray_from_mouse(mouse_pos)
        p_world = self.scene._intersect_ray_plane(r0, rd, np.zeros(3), np.array([0,0,1]))
        
        if p_world is None: return
        x, y = p_world[0], p_world[1]
        
        # 2. Snap to grid
        if self.edit_mode in ('Hex1', 'Hex2'):
            q, r = self.backend.snap_to_ring(x, y)
            node_key = None
        else:
            node_key = self.backend.snap_to_node(x, y)
            q, r = (None, None)
        
        print(f"Click at ({x:.2f}, {y:.2f}) -> Mode={self.edit_mode} Ring={(q,r)} Node={node_key} | Action: {action}")

        # Track last clicked node for pi mode
        if node_key:
            self.last_clicked_node = node_key

        # 3. Modify backend
        if action == 'add':
            if self.edit_mode in ('Hex1', 'Hex2'):
                if q is not None and r is not None:
                    self.backend.add_ring(q, r)
            elif self.edit_mode == 'pi':
                # Cycle pi orbitals: 0 -> 1 -> 2 -> 0
                if node_key:
                    ia = self.backend._build_node_to_atom().get(node_key)
                    if ia is not None:
                        subtype = self.backend.atom_subtype[ia]
                        current_npi = self.backend._get_npi_from_subtype(subtype)
                        new_npi = (current_npi + 1) % 3
                        self.backend.set_atom_valency(node_key, new_npi)
                        print(f"Set atom {node_key} to npi={new_npi}")
            elif node_key:
                print(f"DEBUG: Setting node {node_key} to {self.cur_atom_type}")
                self.backend.set_atom_type(node_key, self.cur_atom_type)
        elif action == 'remove':
            if self.edit_mode in ('Hex1', 'Hex2'):
                if q is not None and r is not None:
                    self.backend.remove_ring(q, r)
            elif node_key:
                self.backend.remove_atom(node_key)
        elif action == 'toggle_h':
            nk = self.backend.snap_to_node(x, y)
            if nk:
                self.backend.toggle_h_state(nk)
        
        self.refresh_view()

    def refresh_view(self):
        # 0. Update Guide Grid
        guides = self.backend.get_guide_points()
        self.grid_markers.set_data(
            pos=np.column_stack([guides, np.full(len(guides), -0.1)]).astype(np.float32),
            symbol='disc', edge_width=0, size=2,
            face_color=(0.3, 0.3, 0.3, 0.3)
        )

        # 0.5. Debug view: for each atom, draw its pin node (cyan disc) and a line atom->pin
        if self.debug_view_mode and hasattr(self.backend, 'atom_pin'):
            pin_pos = []
            line_segs = []
            for ia, pin in enumerate(self.backend.atom_pin):
                if pin is not None and ia < len(self.backend.sys.apos):
                    atom_pos = self.backend.sys.apos[ia]
                    pin_pos.append([pin[0], pin[1], 0.05])          # pin node in z=0 plane
                    line_segs.append([pin[0], pin[1], 0.05])        # line: pin -> atom
                    line_segs.append([atom_pos[0], atom_pos[1], atom_pos[2]])
            if pin_pos:
                self.debug_markers.set_data(
                    pos=np.array(pin_pos, dtype=np.float32),
                    symbol='disc', edge_width=0, face_color=(0.0, 1.0, 1.0, 0.7), size=6
                )
                self.debug_markers.visible = True
                segs = np.array(line_segs, dtype=np.float32)
                conn = np.zeros(len(segs), dtype=bool); conn[0::2] = True  # isolated pairs
                self.debug_lines.set_data(pos=segs, connect=conn, color=(0.0, 1.0, 1.0, 0.6), width=1.5)
                self.debug_lines.visible = True
            else:
                self.debug_markers.visible = False
                self.debug_lines.visible = False
        else:
            self.debug_markers.visible = False
            self.debug_lines.visible = False

        # 1. Use persistent sys directly
        sys = self.backend.sys
        pos = sys.apos.astype(np.float32)

        if pos.size == 0:
            self.scene.set_data(np.zeros((0,3)))
            return

        # Colors based on elements
        colors = []
        sizes = []
        for e in sys.enames:
            c = elements.getColor(e)
            if e == 'H':
                colors.append((0.4, 0.4, 0.4, 1.0))
                sizes.append(8.0)
            else:
                colors.append((c[0], c[1], c[2], 1.0))
                sizes.append(15.0)
        
        colors = np.array(colors, dtype=np.float32)
        sizes = np.array(sizes, dtype=np.float32)
        
        # Bonds
        if sys.bonds is not None:
            is_heavy = np.array([sys.enames[i] != 'H' for i in range(len(sys.enames))])
            bonds_heavy = []
            bonds_h = []
            
            for b in sys.bonds:
                if is_heavy[b[0]] and is_heavy[b[1]]:
                    bonds_heavy.append(b)
                else:
                    bonds_h.append(b)
            
            # Bond color visualization mode
            if self.bond_viz_mode and bonds_heavy:
                bond_segs, bond_colors = compute_bond_colors_by_length(bonds_heavy, pos)
                self.scene._line_set("bonds-colored", self.scene.bond_colored_lines, bond_segs, color=bond_colors, width=5.0)
                self.scene.bond_colored_lines.visible = True
                self.scene.bond_lines.visible = False
                self.scene.set_data(pos, colors=colors, sizes=sizes, bonds=None)
            elif self.bond_viz_mode:
                # Bond viz mode but no bonds - hide both bond visuals
                self.scene.bond_colored_lines.visible = False
                self.scene.bond_lines.visible = False
                self.scene.set_data(pos, colors=colors, sizes=sizes, bonds=None)
            else:
                # Normal bond rendering - use standard bond_lines
                self.scene.bond_colored_lines.visible = False
                self.scene.bond_lines.visible = True
                self.scene.set_data(pos, colors=colors, sizes=sizes, bonds=bonds_heavy)
            
            if bonds_h:
                h_segs = pos[np.array(bonds_h)].reshape(-1, 3)
                self.scene._line_set("CH-bonds", self.scene.ch_bond_lines, h_segs, color=(0.4, 0.4, 0.4, 0.6), width=1.0)
            else:
                self.scene.ch_bond_lines.set_data(np.zeros((0, 3), dtype=np.float32))

            hbonds = sys.find_hbonds(bPrint=False)
            if hbonds:
                hb_segs = []
                for d, h, a, dist, ang in hbonds:
                    hb_segs.append(pos[h])
                    hb_segs.append(pos[a])
                hb_segs = np.array(hb_segs, dtype=np.float32)
                self.scene._line_set("H-bonds", self.scene.hbond_lines, hb_segs, color=(0.8, 0.2, 0.8, 0.5), width=1.5)
            else:
                self.scene.hbond_lines.set_data(np.zeros((0, 3), dtype=np.float32))

            # Labels based on label_mode
            lbl_pos, lbl_texts = generate_atom_labels(self.label_mode, pos, sys.enames, self.backend.atom_subtype, self.backend, bonds_heavy)
            if lbl_pos:
                self.scene.text_labels.text = lbl_texts
                self.scene.text_labels.pos = np.array(lbl_pos, dtype=np.float32)
                self.scene.text_labels.color = np.array([(0, 0, 0, 1)] * len(lbl_texts), dtype=np.float32)
                self.scene.text_labels.visible = True
            else:
                self.scene.text_labels.visible = False

    def run_relaxation(self):
        self.dftb_status_label.setText("Status: Relaxing...")
        self.statusBar().showMessage("Relaxing... please wait")
        QtWidgets.QApplication.processEvents()

        try:
            E, forces, lvs = self.backend.run_relaxation(workdir='gui_relax')
            msg = f"Relaxation done. E = {E:.4f} eV"
            self.statusBar().showMessage(msg)
            self.dftb_status_label.setText(f"Status: Done\nE = {E:.4f} eV")
            self.refresh_view()
        except Exception as e:
            msg = f"Relaxation FAILED: {e}"
            self.statusBar().showMessage(msg)
            self.dftb_status_label.setText(f"Status: FAILED\n{e}")
            self._raise(msg, title="Relaxation Error")

    def compute_orbitals(self):
        if len(self.backend.sys.apos) == 0:
            msg = "No atoms to compute orbitals for."
            print(f"WARNING: {msg}")
            QtWidgets.QMessageBox.warning(self, "Warning", msg)
            return

        # Setup Fdata symlink (needed by FireCore)
        _THIS_DIR = os.path.dirname(os.path.abspath(__file__))
        FDATA_DIR = os.path.join(_THIS_DIR, "Fdata")
        FDATA_TARGET = self.fdata_path

        if not os.path.exists(FDATA_DIR):
            if os.path.exists(FDATA_TARGET):
                os.symlink(FDATA_TARGET, FDATA_DIR)
                print(f"Created symlink: {FDATA_DIR} -> {FDATA_TARGET}")
            else:
                msg = f"Neither {FDATA_DIR} nor {FDATA_TARGET} exists. Please download Fdata_HC_minimal from fireball-qmd.github.io"
                print(f"ERROR: {msg}")
                QtWidgets.QMessageBox.critical(self, "Fdata Error", msg)
                return

        self.statusBar().showMessage("Running FireCore SCF...")
        QtWidgets.QApplication.processEvents()
        try:
            atypes = np.array([elements.ELEMENT_DICT[e][0] for e in self.backend.sys.enames], dtype=np.int32)
            apos = np.array(self.backend.sys.apos, dtype=np.float64)
            fc.setVerbosity(0)
            fc.initialize(atomType=atypes, atomPos=apos)
            fc.evalForce(apos, nmax_scf=200)
            dims = fc.get_HS_dims()
            norb = int(dims.norbitals)
            self._eigen = fc.get_eigen(ikp=1, norb=norb)
            self._wfcoef = fc.get_wfcoef(norb=norb)
            self._norb = norb
            occ = np.where(self._eigen < 0.0)[0]
            self._homo = int(occ[-1]) if len(occ) > 0 else len(self._eigen) // 2 - 1
            self._lumo = self._homo + 1

            # Calculate total valence electrons
            valence_dict = {'H': 1, 'C': 4, 'N': 5, 'O': 6}
            total_electrons = sum([valence_dict.get(e, 0) for e in self.backend.sys.enames])
            occupied_orbitals = total_electrons // 2

            # Update orbital info label
            info_text = (f"Total Orbitals: {norb}\n"
                        f"HOMO: {self._homo + 1} (E={self._eigen[self._homo]:.3f} eV)\n"
                        f"LUMO: {self._lumo + 1} (E={self._eigen[self._lumo]:.3f} eV)\n"
                        f"Occupied: {occupied_orbitals} (e-/2)")
            self.orbital_info_label.setText(info_text)

            # Enable orbital controls
            self.orbital_spinbox.setEnabled(True)
            self.setSpinBox(self.orbital_spinbox, vmin=0, vmax=norb-1, value=self._homo)
            self.plot_orb_btn.setEnabled(True)
            self.plot_density_btn.setEnabled(True)
            self.plot_delta_btn.setEnabled(True)

            # Update orbital energy label for current selection
            self.update_orbital_energy_label(self._homo)

            msg = f"SCF done. HOMO={self._homo + 1} E={self._eigen[self._homo]:.3f} eV  LUMO={self._lumo + 1} E={self._eigen[self._lumo]:.3f} eV"
            print(msg)
            self.statusBar().showMessage(msg)
            QtWidgets.QMessageBox.information(self, "SCF Done", f"HOMO={self._homo + 1} E={self._eigen[self._homo]:.3f} eV\nLUMO={self._lumo + 1} E={self._eigen[self._lumo]:.3f} eV")
        except Exception as e:
            self._raise(f"SCF FAILED: {e}", title="SCF Error")

    def _compute_extent_from_geometry(self, apos, padding_factor=0.1, default_size=14.0):
        """Compute grid extent and origin from atomic positions."""
        if len(apos) > 0:
            apos_2d = apos[:, :2]  # Only x,y
            min_pos = apos_2d.min(axis=0)
            max_pos = apos_2d.max(axis=0)
            center_z = apos[:, 2].mean()
            # Add padding
            padding = (max_pos - min_pos) * padding_factor
            grid_origin = min_pos - padding
            size = (max_pos - min_pos + 2 * padding).max()
        else:
            # No atoms - cannot compute extent from geometry
            # This should not happen in normal usage (plotting requires atoms)
            raise ValueError("Cannot compute grid extent: no atoms in system. Add atoms first.")
        return grid_origin, size, center_z

    def _make_2d_grid(self, grid_origin, size, center_z, z_height, n=100):
        """Generate 2D grid points for projection."""
        xs = np.linspace(grid_origin[0], grid_origin[0] + size, n)
        ys = np.linspace(grid_origin[1], grid_origin[1] + size, n)
        X, Y = np.meshgrid(xs, ys)
        Z = np.zeros_like(X) + (center_z + z_height)
        points = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)
        extent = [grid_origin[0], grid_origin[0] + size, grid_origin[1], grid_origin[1] + size]
        return points, extent, n

    def _evaluate_on_grid(self, points, what, orb_index=None):
        """Call Fortran evaluator on grid points.
        
        Args:
            points: (npoints, 3) array of grid points
            what: 'orbital', 'density', or 'delta_rho'
            orb_index: orbital index (1-based, required for 'orbital')
        
        Returns:
            flat array of evaluated values
        """
        points_f64 = points.astype(np.float64)
        if what == 'orbital':
            if orb_index is None:
                raise ValueError("orb_index required for orbital evaluation")
            return fc.orb2points(points_f64, iMO=int(orb_index + 1), ikpoint=1)
        elif what == 'density':
            return fc.dens2points(points_f64, f_den=1.0, f_den0=0.0)
        elif what == 'delta_rho':
            return fc.dens2points(points_f64, f_den=1.0, f_den0=-1.0)
        else:
            raise ValueError(f"Unknown what={what}, expected 'orbital', 'density', or 'delta_rho'")

    def _plot_2d_projection(self, data_2d, extent, title, cmap, symmetric, atom_pos, atom_types):
        """Plot 2D projection using VisPy heatmap."""
        canvas, view = vu.create_heatmap_window(data_2d, extent, title=title,  cmap=cmap,  symmetric=symmetric, atom_pos=atom_pos, atom_types=atom_types )
        return canvas, view

    def plot_density(self):
        if not hasattr(self, '_eigen'):
            msg = "Please run Compute SCF first."
            print(f"INFO: {msg}")
            QtWidgets.QMessageBox.information(self, "Info", msg)
            return
        self.statusBar().showMessage("Computing electron density...")
        QtWidgets.QApplication.processEvents()
        try:
            apos = self.backend.sys.apos
            z_height = self.z_height_spinbox.value()
            grid_origin, size, center_z = self._compute_extent_from_geometry(apos)
            points, extent, n = self._make_2d_grid(grid_origin, size, center_z, z_height)
            flat_data = self._evaluate_on_grid(points, 'density')
            data_2d = np.asarray(flat_data, dtype=np.float64).reshape(n, n)
            pos = apos.astype(np.float32)
            enames = list(self.backend.sys.enames)
            self._plot_2d_projection( data_2d, extent,title=f"Electron Density (z={z_height:.1f}Å)",cmap='bwr',   symmetric=False,  atom_pos=pos, atom_types=enames  )
            msg = "Density plotted"
            print(msg)
            self.statusBar().showMessage(msg)
        except Exception as e:
            self._raise(f"Density plot FAILED: {e}", title="Plot Error")

    def plot_delta_rho(self):
        if not hasattr(self, '_eigen'):
            msg = "Please run Compute SCF first."
            print(f"INFO: {msg}")
            QtWidgets.QMessageBox.information(self, "Info", msg)
            return
        self.statusBar().showMessage("Computing delta-rho (rho_SCF - rho_NA)...")
        QtWidgets.QApplication.processEvents()
        try:
            apos = self.backend.sys.apos
            z_height = self.z_height_spinbox.value()
            grid_origin, size, center_z = self._compute_extent_from_geometry(apos)
            points, extent, n = self._make_2d_grid(grid_origin, size, center_z, z_height)
            flat_data = self._evaluate_on_grid(points, 'delta_rho')
            data_2d = np.asarray(flat_data, dtype=np.float64).reshape(n, n)
            pos = apos.astype(np.float32)
            enames = list(self.backend.sys.enames)
            self._plot_2d_projection(data_2d, extent, title=f"Delta-Rho (z={z_height:.1f}Å)", cmap='bwr', symmetric=True, atom_pos=pos, atom_types=enames)
            msg = "Delta-rho plotted"
            print(msg)
            self.statusBar().showMessage(msg)
        except Exception as e:
            self._raise(f"Delta-rho plot FAILED: {e}", title="Plot Error")

    def set_fdata_path(self):
        """Open dialog to set Fdata path and save to settings."""
        selected = self.fileDialog(mode="directory", title="Select Fdata Directory", start_dir=self.fdata_path)
        if selected:
            self.fdata_path = selected
            self.settings.setValue("fdata_path", selected)
            print(f"Set Fdata path to: {selected}")
            self.statusBar().showMessage(f"Fdata path set to: {selected}")
            QtWidgets.QMessageBox.information(self, "Settings Saved", f"Fdata path set to:\n{selected}")

if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    window = KekuleExplorerWindow()
    window.show()
    sys.exit(app.exec_())
