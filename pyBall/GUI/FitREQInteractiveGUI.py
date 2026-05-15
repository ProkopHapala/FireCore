#!/usr/bin/env python3
"""
FitREQInteractiveGUI.py - Interactive parameter fitting GUI for pyOpenCL FitREQ

Uses matplotlib for proper axes and PyQt5 for GUI controls.
"""

import sys
import os
import numpy as np

# Remove custom pyopencl from path to use system version
sys.path = [p for p in sys.path if 'SW/pyopencl' not in p]

# Add FireCore to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from PyQt5 import QtWidgets, QtCore, QtGui
from pyBall.GUI.BaseGUI import BaseGUI
from pyBall.OCL.FittingDriver import FittingDriver
from pyBall.OCL.NonBondFitting import extract_macro_block

# Import proper XYZ parsing for grid reshaping
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))
from tests.tFitREQ.split_scan_imshow_new import (
    parse_headers_ra, read_scan_atomicutils, parse_xyz_blocks,
    compute_ra_vec, detect_rows_by_r, reshape_to_grid as reshape_to_grid_proper
)

# Use matplotlib for proper axes
import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

# VisPy for 3D molecular visualization
import vispy
vispy.use('pyqt5')
from vispy import scene
from vispy.scene import visuals as vispy_visuals


class SimpleMoleculeViewer(QtWidgets.QWidget):
    """Simple 3D molecular viewer with camera rotation."""
    
    def __init__(self):
        super().__init__()
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        
        # Create VisPy canvas
        self.canvas = scene.SceneCanvas(keys='interactive', bgcolor='white')
        self.view = self.canvas.central_widget.add_view()
        self.view.camera = scene.TurntableCamera(elevation=30, azimuth=45, distance=20)
        
        # Create visual elements
        self.markers = vispy_visuals.Markers(parent=self.view.scene)
        self.lines = vispy_visuals.Line(parent=self.view.scene, connect='segments', 
                                         color='gray', width=1.5, antialias=True)
        self.text = vispy_visuals.Text(parent=self.view.scene, color='black', 
                                       font_size=10, anchor_x='center', anchor_y='center')
        
        # Add axes
        vispy_visuals.XYZAxis(parent=self.view.scene)
        
        layout.addWidget(self.canvas.native)
    
    def set_molecule(self, positions, enames, bonds=None):
        """Set molecule data for visualization."""
        # Element colors (CPK-like)
        element_colors = {
            'H': (1.0, 1.0, 1.0, 1.0),    # White
            'C': (0.5, 0.5, 0.5, 1.0),    # Gray
            'N': (0.0, 0.0, 1.0, 1.0),    # Blue
            'O': (1.0, 0.0, 0.0, 1.0),    # Red
            'E': (0.0, 1.0, 0.0, 0.5),    # Green (E-pairs)
        }
        
        # Element sizes (van der Waals radii scaled)
        element_sizes = {
            'H': 10,
            'C': 20,
            'N': 20,
            'O': 20,
            'E': 12,
        }
        
        # Get colors and sizes
        colors = np.array([element_colors.get(e, (0.5, 0.5, 0.5, 1.0)) for e in enames], dtype=np.float32)
        sizes = np.array([element_sizes.get(e, 10) for e in enames], dtype=np.float32)
        
        # Set atom markers
        self.markers.set_data(pos=positions, face_color=colors, size=sizes, edge_width=0)
        
        # Set bonds
        if bonds is not None and len(bonds) > 0:
            bond_segments = []
            for i, j in bonds:
                bond_segments.append(positions[i])
                bond_segments.append(positions[j])
            self.lines.set_data(pos=np.array(bond_segments, dtype=np.float32))
        else:
            self.lines.set_data(pos=np.zeros((0, 3), dtype=np.float32))
        
        # Set atom labels
        lbl_texts = [f"{ename}{i}" for i, ename in enumerate(enames)]
        self.text.text = lbl_texts
        self.text.pos = positions.astype(np.float32)
        
        # Center camera on molecule
        center = np.mean(positions, axis=0)
        self.view.camera.center = center
        self.view.camera.distance = np.max(np.ptp(positions, axis=0)) * 3
        
        self.canvas.update()


# Default parameter values - matching atom types in XYZ file (O_3 and H_O)
# E values are stored as sqrt(EvdW) in GPU, but we input absolute values here
DEFAULT_SIMPLE_PARAMS = """# Oxygen O_3 parameters
O_3.R   1.75
O_3.E   0.0026
O_3.Q   0.0
O_3.H   -0.95

# Hydrogen H_O parameters  
H_O.R   1.487
H_O.E   0.00068
H_O.Q   0.0
H_O.H   0.93"""

DEFAULT_DOF_PARAMS = """# typename comp   Min Max     xlo  xhi        Klo   Khi   K0       xstart
E_HO  2     0.0     1.0          0.00   1.00     0.0   0.0   0.0    +0.05
E_O3  2    -1.0    -0.0         -1.00   0.00     0.0   0.0   0.0    -0.05
O_3   3    -1.0     0.0         -0.90  -0.50     1.0   1.0   0.0    -0.5"""


class FitREQInteractiveWindow(BaseGUI):
    """Main window for interactive FitREQ parameter fitting."""
    
    def __init__(self, xyz_file=None, atom_types_file=None, model_name='ENERGY_MorseQ_PAIR'):
        super().__init__("FitREQ Interactive Fitting")
        self.resize(1400, 900)
        
        # File paths
        self.xyz_file = xyz_file or 'wb97m_input/H2O-A1_H2O-D1.xyz'
        self.atom_types_file = atom_types_file or '../../cpp/common_resources/AtomTypes.dat'
        self.model_name = model_name
        
        # Initialize FittingDriver
        self.drv = None
        self.Vref = None
        self.Erefs = None
        self.x0s = None
        
        # Parameter storage
        self.param_dict = {}  # name -> value
        self.param_widgets = {}
        
        self.initUI()
        self.init_fitting_driver()
        
    def initUI(self):
        """Initialize the GUI layout."""
        # Main widget with horizontal layout
        main_widget = QtWidgets.QWidget()
        main_layout = QtWidgets.QHBoxLayout(main_widget)
        
        # --- Left panel: Controls ---
        left_panel = QtWidgets.QFrame()
        left_panel.setFrameStyle(QtWidgets.QFrame.StyledPanel)
        left_panel.setFixedWidth(400)
        left_layout = QtWidgets.QVBoxLayout(left_panel)
        
        # Simple parameters text area
        left_layout.addWidget(QtWidgets.QLabel("Simple Parameters (type.comp value):"))
        self.simple_params_edit = self.textEdit(
            text=DEFAULT_SIMPLE_PARAMS,
            plain=True,
            min_size=(380, 200)
        )
        left_layout.addWidget(self.simple_params_edit)
        
        # Apply simple params button
        self.button("Apply Simple Params", self.apply_simple_params, layout=left_layout)
        
        # Full DOF parameters text area
        left_layout.addWidget(QtWidgets.QLabel("Full DOF Parameters:"))
        self.dof_params_edit = self.textEdit(
            text=DEFAULT_DOF_PARAMS,
            plain=True,
            min_size=(380, 150)
        )
        left_layout.addWidget(self.dof_params_edit)
        
        # Apply DOF params button
        self.button("Apply DOF Params", self.apply_dof_params, layout=left_layout)
        
        # --- Parameter selector ---
        left_layout.addWidget(QtWidgets.QLabel("Parameter Selector:"))
        selector_layout = QtWidgets.QHBoxLayout()
        
        self.param_combo = QtWidgets.QComboBox()
        self.param_combo.setEditable(True)
        self.populate_param_combo()
        self.param_combo.currentTextChanged.connect(self.on_param_selected)
        selector_layout.addWidget(self.param_combo)
        
        self.param_spin = QtWidgets.QDoubleSpinBox()
        self.param_spin.setRange(-100, 100)
        self.param_spin.setDecimals(6)
        self.param_spin.setSingleStep(0.01)
        # Connect valueChanged to apply params and recompute
        self.param_spin.valueChanged.connect(self.on_param_value_changed)
        selector_layout.addWidget(self.param_spin)
        
        left_layout.addLayout(selector_layout)
        
        # Mousewheel sensitivity
        left_layout.addWidget(QtWidgets.QLabel("Mousewheel Step:"))
        self.wheel_step_spin = QtWidgets.QDoubleSpinBox()
        self.wheel_step_spin.setRange(0.0001, 1.0)
        self.wheel_step_spin.setDecimals(4)
        self.wheel_step_spin.setValue(0.01)
        self.wheel_step_spin.setSingleStep(0.001)
        left_layout.addWidget(self.wheel_step_spin)
        
        # Recompute button
        self.button("Recompute Energy", self.recompute_energy, layout=left_layout)
        
        # Status label
        self.status_label = QtWidgets.QLabel("Status: Ready")
        left_layout.addWidget(self.status_label)
        
        # --- Pixel Selection Controls ---
        left_layout.addWidget(QtWidgets.QLabel("Pixel Selection (row, col):"))
        pixel_layout = QtWidgets.QHBoxLayout()
        
        self.pixel_row_spin = QtWidgets.QSpinBox()
        self.pixel_row_spin.setRange(0, 9999)
        self.pixel_row_spin.setValue(0)
        self.pixel_row_spin.valueChanged.connect(self.on_pixel_changed)
        pixel_layout.addWidget(QtWidgets.QLabel("Row:"))
        pixel_layout.addWidget(self.pixel_row_spin)
        
        self.pixel_col_spin = QtWidgets.QSpinBox()
        self.pixel_col_spin.setRange(0, 9999)
        self.pixel_col_spin.setValue(0)
        self.pixel_col_spin.valueChanged.connect(self.on_pixel_changed)
        pixel_layout.addWidget(QtWidgets.QLabel("Col:"))
        pixel_layout.addWidget(self.pixel_col_spin)
        
        left_layout.addLayout(pixel_layout)
        
        # Show 3D and Decomposition buttons
        self.button("Show 3D Geometry", self.show_3d_geometry, layout=left_layout)
        self.button("Show Interaction Matrix", self.show_interaction_matrix, layout=left_layout)
        
        left_layout.addStretch()
        
        # --- Right panel: Tab widget with multiple views ---
        self.tab_widget = QtWidgets.QTabWidget()
        
        # Tab 1: Energy Maps (matplotlib)
        self.energy_tab = QtWidgets.QWidget()
        energy_layout = QtWidgets.QVBoxLayout(self.energy_tab)
        
        self.fig = Figure(figsize=(14, 5), dpi=100)
        self.fig.patch.set_facecolor('white')
        self.canvas = FigureCanvas(self.fig)
        
        # Create 3 subplots with shared axes
        self.ax_ref = self.fig.add_subplot(131)
        self.ax_mod = self.fig.add_subplot(132, sharex=self.ax_ref, sharey=self.ax_ref)
        self.ax_diff = self.fig.add_subplot(133, sharex=self.ax_ref, sharey=self.ax_ref)
        
        # Initialize empty plots (will be updated with real data)
        self.im_ref = self.ax_ref.imshow([[0]], aspect='auto', cmap='viridis', origin='lower')
        self.im_mod = self.ax_mod.imshow([[0]], aspect='auto', cmap='viridis', origin='lower')
        self.im_diff = self.ax_diff.imshow([[0]], aspect='auto', cmap='bwr', origin='lower')
        
        # Add colorbars once
        self.cbar_ref = self.fig.colorbar(self.im_ref, ax=self.ax_ref, fraction=0.046, pad=0.04)
        self.cbar_mod = self.fig.colorbar(self.im_mod, ax=self.ax_mod, fraction=0.046, pad=0.04)
        self.cbar_diff = self.fig.colorbar(self.im_diff, ax=self.ax_diff, fraction=0.046, pad=0.04)
        
        # Set labels once
        for ax in [self.ax_ref, self.ax_mod, self.ax_diff]:
            ax.set_xlabel('Distance (A)')
            ax.set_ylabel('Angle (deg)')
        
        # Tight layout once
        self.fig.tight_layout()
        
        # Connect mousewheel to parameter adjustment (on canvas)
        self.canvas.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.canvas.wheelEvent = self.on_matplotlib_wheel
        
        # Connect click event for pixel selection
        self.canvas.mpl_connect('button_press_event', self.on_plot_click)
        
        energy_layout.addWidget(self.canvas)
        self.tab_widget.addTab(self.energy_tab, "Energy Maps")
        
        # Tab 2: 3D Geometry (VisPy)
        self.geometry_tab = QtWidgets.QWidget()
        geometry_layout = QtWidgets.QVBoxLayout(self.geometry_tab)
        
        # Create SimpleMoleculeViewer for 3D visualization
        self.molecule_viewer = SimpleMoleculeViewer()
        geometry_layout.addWidget(self.molecule_viewer)
        
        self.tab_widget.addTab(self.geometry_tab, "3D Geometry")
        
        # Tab 3: Interaction Matrix (matplotlib)
        self.matrix_tab = QtWidgets.QWidget()
        matrix_layout = QtWidgets.QVBoxLayout(self.matrix_tab)
        
        self.fig_matrix = Figure(figsize=(14, 5), dpi=100)
        self.fig_matrix.patch.set_facecolor('white')
        self.canvas_matrix = FigureCanvas(self.fig_matrix)
        
        # Create 4 subplots for each interaction component
        self.ax_pauli = self.fig_matrix.add_subplot(141)
        self.ax_london = self.fig_matrix.add_subplot(142)
        self.ax_elec = self.fig_matrix.add_subplot(143)
        self.ax_hbond = self.fig_matrix.add_subplot(144)
        
        self.ax_pauli.set_title('Pauli Repulsion')
        self.ax_london.set_title('London Dispersion')
        self.ax_elec.set_title('Electrostatic')
        self.ax_hbond.set_title('H-Bond')
        
        self.fig_matrix.tight_layout()
        matrix_layout.addWidget(self.canvas_matrix)
        
        self.tab_widget.addTab(self.matrix_tab, "Interaction Matrix")
        
        # Add panels to main layout
        main_layout.addWidget(left_panel)
        main_layout.addWidget(self.tab_widget)
        
        self.setCentralWidget(main_widget)
        
        # Status bar
        self.statusBar().showMessage("Ready - Load XYZ and click Recompute")
        
    def populate_param_combo(self):
        """Populate parameter combo box with available parameters."""
        # Parse simple params to get names
        params = self.parse_simple_params()
        self.param_combo.clear()
        for name in sorted(params.keys()):
            self.param_combo.addItem(name)
        
    def init_fitting_driver(self):
        """Initialize the pyOpenCL FittingDriver."""
        try:
            self.status_label.setText("Initializing FittingDriver...")
            QtWidgets.QApplication.processEvents()
            
            # Resolve paths
            if not os.path.isabs(self.xyz_file):
                self.xyz_file = os.path.join(os.path.dirname(__file__), '../../tests/tFitREQ', self.xyz_file)
            if not os.path.isabs(self.atom_types_file):
                self.atom_types_file = os.path.join(os.path.dirname(__file__), self.atom_types_file)
            
            # Create driver (will select NVIDIA GPU automatically)
            self.drv = FittingDriver(verbose=1)
            self.drv.load_atom_types(self.atom_types_file)
            
            # Compile energy kernel template (required before evaluate_energies)
            # Go up 3 levels: GUI -> pyBall -> FireCore root
            base_path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
            forces_path = os.path.join(base_path, 'cpp/common_resources/cl/Forces.cl')
            macro_code = extract_macro_block(forces_path, self.model_name)
            macros = {'MODEL_PAIR_ENERGY': macro_code}
            self.drv.compile_energy_with_model(macros=macros, bPrint=False)
            
            self.status_label.setText("FittingDriver initialized")
            
        except Exception as e:
            self.status_label.setText(f"Error: {e}")
            print(f"Error initializing FittingDriver: {e}")
            
    def parse_simple_params(self):
        """Parse simple parameter format: 'E_HO.R  1.15' -> {E_HO.R: 1.15}"""
        params = {}
        text = self.simple_params_edit.toPlainText()
        for line in text.strip().split('\n'):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 2:
                name = parts[0]
                try:
                    value = float(parts[1])
                    params[name] = value
                except ValueError:
                    continue
        return params
        
    def parse_dof_params(self):
        """Parse full DOF parameter format."""
        dofs = []
        text = self.dof_params_edit.toPlainText()
        for line in text.strip().split('\n'):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 10:
                dofs.append({
                    'typename': parts[0],
                    'comp': int(parts[1]),
                    'min': float(parts[2]),
                    'max': float(parts[3]),
                    'xlo': float(parts[4]),
                    'xhi': float(parts[5]),
                    'Klo': float(parts[6]),
                    'Khi': float(parts[7]),
                    'K0': float(parts[8]),
                    'xstart': float(parts[9])
                })
        return dofs
        
    def apply_simple_params(self):
        """Apply simple parameters to FittingDriver."""
        if self.drv is None:
            self.status_label.setText("Error: FittingDriver not initialized")
            return
            
        params = self.parse_simple_params()
        
        # Build type overrides
        type_overrides = {}
        for name, value in params.items():
            if '.' in name:
                type_name, comp = name.split('.')
                comp_idx = {'R': 0, 'E': 1, 'Q': 2, 'H': 3}.get(comp, -1)
                if comp_idx >= 0:
                    if type_name not in type_overrides:
                        type_overrides[type_name] = [None, None, None, None]
                    type_overrides[type_name][comp_idx] = value
        
        # Apply overrides using FittingDriver method
        if hasattr(self.drv, 'apply_type_overrides'):
            # Convert to type_set format: [(type_name, comp, value), ...]
            type_set = []
            comp_names = ['R', 'E', 'Q', 'H']
            for type_name, values in type_overrides.items():
                for i, val in enumerate(values):
                    if val is not None:
                        type_set.append((type_name, comp_names[i], val))
            self.drv.type_set = type_set
            self.drv.apply_type_overrides()
            self.status_label.setText(f"Applied {len(type_set)} type overrides")
        else:
            # Manual override via tREQHs_base
            for type_name, reqh in type_overrides.items():
                if type_name in self.drv.base_params:
                    base = self.drv.base_params[type_name]
                    if reqh[0] is not None:
                        base['R'] = reqh[0]
                    if reqh[1] is not None:
                        base['E'] = reqh[1]
                    if reqh[2] is not None:
                        base['Q'] = reqh[2]
                    if reqh[3] is not None:
                        base['H'] = reqh[3]
            self.status_label.setText(f"Applied {len(type_overrides)} type overrides (manual)")
        
        # Update parameter combo
        self.populate_param_combo()
        
        # Recompute
        self.recompute_energy()
        
    def apply_dof_params(self):
        """Apply full DOF parameters."""
        if self.drv is None:
            self.status_label.setText("Error: FittingDriver not initialized")
            return
            
        dofs = self.parse_dof_params()
        self.drv.dof_definitions = dofs
        self.status_label.setText(f"Applied {len(dofs)} DOF definitions")
        
    def on_param_value_changed(self, value):
        """Handle spin box value change - update params and recompute."""
        param_name = self.param_combo.currentText()
        if not param_name or '.' not in param_name:
            return
        
        print(f"\n=== SPIN BOX VALUE CHANGED ===")
        print(f"  param_name: {param_name}")
        print(f"  new_value: {value}")
        
        # Update simple params text
        self.update_simple_param_text(param_name, value)
        
        # Apply params and recompute
        print(f"  calling recompute_energy()...")
        self.recompute_energy()
        print(f"=== END SPIN BOX VALUE CHANGED ===\n")
        
    def on_param_selected(self, param_name):
        """Handle parameter selection change - update spin box to show current value."""
        if not param_name or '.' not in param_name:
            return
        
        # Parse current simple params to get value
        params = self.parse_simple_params()
        current_value = params.get(param_name, 0.0)
        
        # Update spin box without triggering callbacks
        self.param_spin.blockSignals(True)
        self.param_spin.setValue(current_value)
        self.param_spin.blockSignals(False)
        
        print(f"Selected parameter: {param_name} = {current_value}")
            
    def on_matplotlib_wheel(self, event):
        """Handle mousewheel on matplotlib canvas to adjust selected parameter."""
        if not self.param_combo.currentText():
            return
            
        param_name = self.param_combo.currentText()
        step = self.wheel_step_spin.value()
        current_value = self.param_spin.value()
        
        # Adjust based on wheel direction
        if event.angleDelta().y() > 0:
            new_value = current_value + step
        else:
            new_value = current_value - step
            
        # Update spin box
        self.param_spin.blockSignals(True)
        self.param_spin.setValue(new_value)
        self.param_spin.blockSignals(False)
        
        # Update simple params text and recompute
        self.update_simple_param_text(param_name, new_value)
        self.recompute_energy()
        
        # Accept event
        event.accept()
        
    def update_simple_param_text(self, param_name, value):
        """Update the simple params text area with new value."""
        print(f"    [update_simple_param_text] param={param_name}, value={value}")
        text = self.simple_params_edit.toPlainText()
        lines = text.split('\n')
        new_lines = []
        found = False
        for line in lines:
            if line.strip().startswith(param_name):
                new_line = f"{param_name}  {value:.6g}"
                new_lines.append(new_line)
                found = True
                print(f"      REPLACED: '{line}' -> '{new_line}'")
            else:
                new_lines.append(line)
        if not found:
            new_line = f"{param_name}  {value:.6g}"
            new_lines.append(new_line)
            print(f"      ADDED new line: '{new_line}'")
        new_text = '\n'.join(new_lines)
        self.simple_params_edit.setPlainText(new_text)
        print(f"    [update_simple_param_text] text area updated")
        
    def apply_current_params(self):
        """Apply current simple parameters to FittingDriver tREQHs."""
        print(f"  [apply_current_params] START")
        
        if self.drv is None:
            print(f"    ERROR: drv is None")
            return
        if not hasattr(self.drv, 'tREQHs_base'):
            print(f"    ERROR: tREQHs_base not found")
            return
            
        # Parse current simple params
        print(f"    parsing simple params...")
        params = self.parse_simple_params()
        print(f"    parsed {len(params)} params: {list(params.items())[:4]}...")
        
        # Build type_set list for FittingDriver
        type_set = []
        for name, value in params.items():
            if '.' in name:
                type_name, comp = name.split('.')
                if comp in ['R', 'E', 'Q', 'H']:
                    type_set.append((type_name, comp, value))
        
        print(f"    built type_set with {len(type_set)} entries")
        
        # Apply to FittingDriver
        if type_set:
            print(f"\n    --- Applying Parameter Overrides ---")
            for type_name, comp, value in type_set:
                if type_name in self.drv.atom_type_map:
                    i = self.drv.atom_type_map[type_name]
                    ci = {'R':0, 'E':1, 'Q':2, 'H':3}[comp]
                    base_val = self.drv.tREQHs_base[i, ci]
                    if comp == 'E':
                        base_val = base_val**2
                    print(f"      {type_name}.{comp}: {base_val:.6f} -> {value:.6f}")
                else:
                    print(f"      {type_name}.{comp}: {value:.6f} (TYPE NOT IN MAP!)")
            
            print(f"    calling drv.apply_type_overrides()...")
            self.drv.type_set = type_set
            self.drv.apply_type_overrides()
            
            print(f"    uploading tREQHs to GPU...")
            import pyopencl as cl
            cl.enqueue_copy(self.drv.queue, self.drv.tREQHs_buff, self.drv.tREQHs_base)
            self.drv.queue.finish()
            print(f"    GPU upload complete")
            
            print(f"    Updated tREQHs:")
            for i, (name, params_arr) in enumerate(zip(self.drv.atom_type_names, self.drv.tREQHs_base)):
                R, E_sqrt, Q, H = params_arr
                E = E_sqrt**2
                print(f"      {name}: R={R:.4f} E={E:.6f} Q={Q:.4f} H={H:.4f}")
            print(f"    [apply_current_params] DONE - Applied {len(type_set)} overrides\n")
        else:
            print(f"    WARNING: type_set is empty, nothing to apply")
            
    def recompute_energy(self):
        """Recompute energy using pyOpenCL and update plots."""
        if self.drv is None:
            self.status_label.setText("Error: FittingDriver not initialized")
            return
            
        try:
            self.status_label.setText("Recomputing...")
            self.statusBar().showMessage("Recomputing energy...")
            QtWidgets.QApplication.processEvents()
            
            # Load data if not already loaded
            if not hasattr(self.drv, 'host_ranges'):
                self.drv.load_data(self.xyz_file)
                self.drv.init_and_upload_energy_only()
                self.drv.setup_energy_kernel()
                
                # Extract reference energies
                self.Erefs = self.drv.host_ErefW if hasattr(self.drv, 'host_ErefW') else None
                if self.Erefs is not None:
                    # Handle 2D case (n_samples, 2) - take first column
                    if self.Erefs.ndim > 1:
                        self.Erefs = self.Erefs[:, 0]
                
                # Parse XYZ geometry (r=radius/distance, a=angle) from headers
                # Try atomicUtils first, fallback to local parser
                Es_geo, Ps = read_scan_atomicutils(self.xyz_file)
                if Es_geo.size == 0:
                    # Fallback to local parser
                    _, _, Ps = parse_xyz_blocks(self.xyz_file, natoms=None)
                
                # Compute r, a vectors from positions
                self.r, self.a = compute_ra_vec(Ps, signed=True)
                
                # Parse headers for r, a overrides (returns Eh, Rh, Ah)
                Eh_header, Rh, Ah = parse_headers_ra(self.xyz_file)
                if Rh.size == self.r.size and np.any(np.isfinite(Rh)):
                    self.r = np.where(np.isfinite(Rh), Rh, self.r)
                if Ah.size == self.a.size and np.any(np.isfinite(Ah)):
                    self.a = np.where(np.isfinite(Ah), Ah, self.a)
                
                # Detect rows by r values
                self.rows, _ = detect_rows_by_r(self.r)
                
                # Reshape reference energies to grid
                if self.Erefs is not None:
                    self.Vref, self.Rg, self.Arow, self.rv = reshape_to_grid_proper(
                        self.Erefs, self.r, self.a, self.rows
                    )
                    print(f"Grid shape: Vref={self.Vref.shape}, rows={len(self.rows)}")
                    print(f"  Distance range: [{self.rv.min():.2f}, {self.rv.max():.2f}] A")
                    print(f"  Angle range: [{self.Arow.min():.1f}, {self.Arow.max():.1f}] deg")
                    
                    # Update pixel spin box ranges
                    ny, nx = self.Vref.shape
                    self.pixel_row_spin.setRange(0, ny - 1)
                    self.pixel_col_spin.setRange(0, nx - 1)
                    
                    # Load and store all XYZ configurations
                    self._load_xyz_configurations()
                    
                    # Create pixel to index mapping
                    self._create_pixel_mapping()
                    
                    self.update_plots(self.Vref, self.Vref)
            
            # Apply current parameters to FittingDriver
            self.apply_current_params()
            
            # Evaluate model energies
            print(f"    calling drv.evaluate_energies()...")
            Em = self.drv.evaluate_energies()
            print(f"    evaluate_energies returned: shape={Em.shape}, min={np.nanmin(Em):.6f}, max={np.nanmax(Em):.6f}")
            
            # Reshape model energies using same geometry
            print(f"    reshaping model energies...")
            if hasattr(self, 'r') and hasattr(self, 'rows'):
                Vmod, _, _, _ = reshape_to_grid_proper(Em, self.r, self.a, self.rows)
                print(f"    reshaped to: Vmod shape={Vmod.shape}")
            else:
                print(f"    WARNING: using fallback reshape (no geometry)")
                Vmod = self.reshape_to_grid_simple(Em)
            
            # Update plots
            print(f"    updating plots...")
            if self.Vref is not None:
                self.update_plots(self.Vref, Vmod)
                
                # Compute and display error
                diff = Vmod - self.Vref
                rmse = np.sqrt(np.nanmean(diff**2))
                print(f"    RMSE: {rmse:.6f} eV")
                self.status_label.setText(f"RMSE: {rmse:.6f} eV")
                self.statusBar().showMessage(f"Recompute done - RMSE: {rmse:.6f} eV")
            else:
                print(f"    no reference data, plotting model only")
                self.status_label.setText("No reference data")
                self.update_plots(Vmod, Vmod)
            print(f"  [recompute_energy] DONE\n")
                
        except Exception as e:
            self.status_label.setText(f"Error: {e}")
            self.statusBar().showMessage(f"Error: {e}")
            print(f"Error recomputing energy: {e}")
            import traceback
            traceback.print_exc()
            
    def reshape_to_grid_simple(self, E):
        """Simple factor-based reshaping (fallback only - use reshape_to_grid_proper for correct scan geometry)."""
        # Handle 2D case - take first column if needed
        if E.ndim > 1:
            E = E[:, 0] if E.shape[1] > 0 else E.flatten()
        
        # Try to detect grid shape from data
        n = len(E)
        
        # Try to find square-ish factors
        nx = int(np.sqrt(n))
        while nx > 1 and n % nx != 0:
            nx -= 1
        ny = n // nx if nx > 0 else n
        
        # Reshape
        try:
            return E.reshape((ny, nx))
        except:
            # Fallback: pad with NaN
            padded = np.full(nx * ny, np.nan)
            padded[:n] = E
            return padded.reshape((ny, nx))
            
    def update_plots(self, Vref, Vmod):
        """Update matplotlib plots with new data using blitting."""
        Vdif = Vmod - Vref
        
        # Get extent for proper axis scaling
        if hasattr(self, 'rv') and hasattr(self, 'Arow'):
            r_min, r_max = np.nanmin(self.rv), np.nanmax(self.rv)
            a_min, a_max = np.nanmin(self.Arow), np.nanmax(self.Arow)
            extent = [r_min, r_max, a_min, a_max]
        else:
            extent = None
        
        # Update data only (fast blitting)
        self.im_ref.set_data(Vref)
        self.im_mod.set_data(Vmod)
        self.im_diff.set_data(Vdif)
        
        # Update extents
        if extent:
            self.im_ref.set_extent(extent)
            self.im_mod.set_extent(extent)
            self.im_diff.set_extent(extent)
        
        # Update clim (color limits)
        vref_min, vref_max = np.nanmin(Vref), np.nanmax(Vref)
        vmod_min, vmod_max = np.nanmin(Vmod), np.nanmax(Vmod)
        vdif_min, vdif_max = np.nanmin(Vdif), np.nanmax(Vdif)
        
        self.im_ref.set_clim(vref_min, vref_max)
        self.im_mod.set_clim(vmod_min, vmod_max)
        self.im_diff.set_clim(vdif_min, vdif_max)
        
        # Update colorbar limits (mappable.set_clim updates the colorbar)
        self.im_ref.set_clim(vref_min, vref_max)
        self.im_mod.set_clim(vmod_min, vmod_max)
        self.im_diff.set_clim(vdif_min, vdif_max)
        
        # Update titles
        self.ax_ref.set_title(f'Ref [{vref_min:.3f}, {vref_max:.3f}]')
        self.ax_mod.set_title(f'Mod [{vmod_min:.3f}, {vmod_max:.3f}]')
        self.ax_diff.set_title(f'Diff [{vdif_min:.3f}, {vdif_max:.3f}]')
        
        # Update canvas (fast - only redraws changed artists)
        self.im_ref.axes.figure.canvas.draw_idle()
        self.canvas.draw()
    
    def _load_xyz_configurations(self):
        """Load all XYZ configurations from file."""
        self.xyz_data = []
        
        try:
            with open(self.xyz_file, 'r') as f:
                content = f.read()
            
            # Split into blocks (each configuration is a block)
            lines = content.strip().split('\n')
            i = 0
            config_idx = 0
            
            while i < len(lines):
                # Read number of atoms
                if not lines[i].strip():
                    i += 1
                    continue
                    
                try:
                    natoms = int(lines[i].strip())
                except ValueError:
                    i += 1
                    continue
                
                # Read comment line
                i += 1
                comment = lines[i].strip() if i < len(lines) else ""
                
                # Read atom lines
                i += 1
                atoms = []
                for j in range(natoms):
                    if i >= len(lines):
                        break
                    parts = lines[i].split()
                    if len(parts) >= 4:
                        atom_name = parts[0]
                        x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                        atoms.append({
                            'name': atom_name,
                            'pos': np.array([x, y, z], dtype=np.float32)
                        })
                    i += 1
                
                if atoms:
                    self.xyz_data.append({
                        'idx': config_idx,
                        'natoms': natoms,
                        'comment': comment,
                        'atoms': atoms
                    })
                    config_idx += 1
            
            print(f"Loaded {len(self.xyz_data)} XYZ configurations")
            
        except Exception as e:
            print(f"Error loading XYZ configurations: {e}")
            self.xyz_data = []
    
    def _create_pixel_mapping(self):
        """Create mapping from (row, col) pixel to configuration index."""
        self.pixel_to_idx = {}
        
        if not hasattr(self, 'rows') or self.Vref is None:
            return
        
        ny, nx = self.Vref.shape
        idx = 0
        
        # Map each pixel to its corresponding configuration index
        for row in range(ny):
            for col in range(nx):
                if idx < len(self.xyz_data):
                    self.pixel_to_idx[(row, col)] = idx
                    idx += 1
                else:
                    break
        
        print(f"Created pixel mapping for {len(self.pixel_to_idx)} pixels")
    
    def on_plot_click(self, event):
        """Handle click on matplotlib plot to select pixel."""
        if event.inaxes in [self.ax_ref, self.ax_mod, self.ax_diff]:
            # Get pixel coordinates from data coordinates
            x, y = event.xdata, event.ydata
            if x is None or y is None:
                return
            
            # Convert to pixel indices using extent
            if hasattr(self, 'rv') and hasattr(self, 'Arow') and self.Vref is not None:
                ny, nx = self.Vref.shape
                r_min, r_max = np.nanmin(self.rv), np.nanmax(self.rv)
                a_min, a_max = np.nanmin(self.Arow), np.nanmax(self.Arow)
                
                # Map data coordinates to pixel indices
                col = int((x - r_min) / (r_max - r_min) * (nx - 1))
                row = int((y - a_min) / (a_max - a_min) * (ny - 1))
                
                # Clamp to valid range
                col = max(0, min(col, nx - 1))
                row = max(0, min(row, ny - 1))
                
                # Update spin boxes
                self.pixel_row_spin.blockSignals(True)
                self.pixel_col_spin.blockSignals(True)
                self.pixel_row_spin.setValue(row)
                self.pixel_col_spin.setValue(col)
                self.pixel_row_spin.blockSignals(False)
                self.pixel_col_spin.blockSignals(False)
                
                # Trigger pixel change
                self.on_pixel_changed()
                
                print(f"Selected pixel: row={row}, col={col}")
    
    def on_pixel_changed(self):
        """Handle pixel selection change."""
        row = self.pixel_row_spin.value()
        col = self.pixel_col_spin.value()
        
        # Update status
        self.status_label.setText(f"Selected pixel: ({row}, {col})")
    
    def show_3d_geometry(self):
        """Display 3D geometry for selected pixel."""
        if not hasattr(self, 'xyz_data') or not self.xyz_data:
            print("No XYZ data loaded")
            return
        
        row = self.pixel_row_spin.value()
        col = self.pixel_col_spin.value()
        
        # Map pixel to configuration index
        if hasattr(self, 'pixel_to_idx'):
            idx = self.pixel_to_idx.get((row, col), 0)
        else:
            # Fallback: assume regular grid
            if self.Vref is not None:
                ny, nx = self.Vref.shape
                idx = row * nx + col
            else:
                idx = 0
        
        # Get configuration data
        if idx >= len(self.xyz_data):
            print(f"Invalid config index: {idx}")
            return
        
        config = self.xyz_data[idx]
        atoms = config['atoms']
        
        # Extract positions and element names
        positions = []
        enames = []
        
        for i, atom in enumerate(atoms):
            ename = atom['name']
            pos = atom['pos']
            positions.append(pos)
            enames.append(ename)
        
        positions = np.array(positions, dtype=np.float32)
        
        # Compute bonds (simple distance-based)
        bonds = self._compute_bonds(positions, enames)
        
        # Update molecule viewer
        self.molecule_viewer.set_molecule(positions, enames, bonds)
        
        # Switch to geometry tab
        self.tab_widget.setCurrentIndex(1)
        
        print(f"Showing 3D geometry for config {idx} (pixel {row},{col}) with {len(positions)} atoms")
    
    def _compute_bonds(self, positions, enames, max_bond_length=2.0):
        """Compute bonds based on distance."""
        bonds = []
        n = len(positions)
        for i in range(n):
            for j in range(i + 1, n):
                dist = np.linalg.norm(positions[i] - positions[j])
                # Simple bond detection: H-X bonds shorter, others longer
                if enames[i] == 'H' or enames[j] == 'H':
                    if dist < 1.2:
                        bonds.append((i, j))
                else:
                    if dist < max_bond_length:
                        bonds.append((i, j))
        return np.array(bonds, dtype=np.int32) if bonds else np.zeros((0, 2), dtype=np.int32)
    
    def _get_atom_styles(self, enames):
        """Get atom colors and sizes based on element names."""
        colors = []
        sizes = []
        
        # Element colors (CPK-like)
        element_colors = {
            'H': (1.0, 1.0, 1.0, 1.0),    # White
            'C': (0.5, 0.5, 0.5, 1.0),    # Gray
            'N': (0.0, 0.0, 1.0, 1.0),    # Blue
            'O': (1.0, 0.0, 0.0, 1.0),    # Red
            'E': (0.0, 1.0, 0.0, 0.5),    # Green (E-pairs, transparent)
        }
        
        # Element sizes (van der Waals radii scaled)
        element_sizes = {
            'H': 5,
            'C': 15,
            'N': 15,
            'O': 15,
            'E': 8,  # E-pairs smaller
        }
        
        for ename in enames:
            color = element_colors.get(ename, (0.5, 0.5, 0.5, 1.0))
            size = element_sizes.get(ename, 10)
            colors.append(color)
            sizes.append(size)
        
        return np.array(colors, dtype=np.float32), np.array(sizes, dtype=np.float32)
    
    def _generate_labels(self, positions, enames):
        """Generate atom labels with indices."""
        lbl_pos = []
        lbl_texts = []
        
        for i, ename in enumerate(enames):
            lbl_pos.append(positions[i])
            lbl_texts.append(f"{ename}{i}")
        
        return lbl_pos, lbl_texts
    
    def show_interaction_matrix(self):
        """Display interaction matrix for selected pixel."""
        if not hasattr(self, 'drv') or self.drv is None:
            print("FittingDriver not initialized")
            return
        
        if not hasattr(self, 'pixel_to_idx'):
            print("No pixel mapping available. Recompute energy first.")
            return
        
        row = self.pixel_row_spin.value()
        col = self.pixel_col_spin.value()
        
        # Map pixel to configuration index
        idx = self.pixel_to_idx.get((row, col), 0)
        
        try:
            print(f"Evaluating interaction matrix for config {idx} (pixel {row},{col})...")
            
            # Call the decomposition kernel
            interactions = self.drv.evaluate_interaction_matrix(idx)
            
            print(f"Got interaction matrix: shape={interactions.shape}")
            print(f"  Pauli: [{interactions[:,:,0].min():.6f}, {interactions[:,:,0].max():.6f}]")
            print(f"  London: [{interactions[:,:,1].min():.6f}, {interactions[:,:,1].max():.6f}]")
            print(f"  Electro: [{interactions[:,:,2].min():.6f}, {interactions[:,:,2].max():.6f}]")
            print(f"  H-bond: [{interactions[:,:,3].min():.6f}, {interactions[:,:,3].max():.6f}]")
            
            # Update plots
            self._update_interaction_plots(interactions)
            
            # Switch to matrix tab
            self.tab_widget.setCurrentIndex(2)
            
            self.status_label.setText(f"Interaction matrix for pixel ({row}, {col})")
            
        except Exception as e:
            print(f"Error evaluating interaction matrix: {e}")
            import traceback
            traceback.print_exc()
    
    def _update_interaction_plots(self, interactions):
        """Update interaction matrix plots."""
        ni, nj, _ = interactions.shape
        
        # Clear previous plots
        self.ax_pauli.clear()
        self.ax_london.clear()
        self.ax_elec.clear()
        self.ax_hbond.clear()
        
        # Plot each component
        im1 = self.ax_pauli.imshow(interactions[:,:,0], aspect='auto', cmap='hot', origin='lower')
        im2 = self.ax_london.imshow(interactions[:,:,1], aspect='auto', cmap='cool', origin='lower')
        im3 = self.ax_elec.imshow(interactions[:,:,2], aspect='auto', cmap='seismic', origin='lower')
        im4 = self.ax_hbond.imshow(interactions[:,:,3], aspect='auto', cmap='viridis', origin='lower')
        
        # Add colorbars
        self.fig_matrix.colorbar(im1, ax=self.ax_pauli, fraction=0.046, pad=0.04)
        self.fig_matrix.colorbar(im2, ax=self.ax_london, fraction=0.046, pad=0.04)
        self.fig_matrix.colorbar(im3, ax=self.ax_elec, fraction=0.046, pad=0.04)
        self.fig_matrix.colorbar(im4, ax=self.ax_hbond, fraction=0.046, pad=0.04)
        
        # Set titles
        self.ax_pauli.set_title(f'Pauli [{interactions[:,:,0].min():.3f}, {interactions[:,:,0].max():.3f}]')
        self.ax_london.set_title(f'London [{interactions[:,:,1].min():.3f}, {interactions[:,:,1].max():.3f}]')
        self.ax_elec.set_title(f'Electro [{interactions[:,:,2].min():.3f}, {interactions[:,:,2].max():.3f}]')
        self.ax_hbond.set_title(f'H-bond [{interactions[:,:,3].min():.3f}, {interactions[:,:,3].max():.3f}]')
        
        # Set labels
        for ax in [self.ax_pauli, self.ax_london, self.ax_elec, self.ax_hbond]:
            ax.set_xlabel('Fragment 2 Atom')
            ax.set_ylabel('Fragment 1 Atom')
        
        # Update canvas
        self.canvas_matrix.draw()


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Interactive FitREQ parameter fitting')
    parser.add_argument('--xyz', type=str, default=None, help='XYZ file')
    parser.add_argument('--atom-types', type=str, default=None, help='Atom types file')
    parser.add_argument('--model', type=str, default='ENERGY_MorseQ_PAIR', help='Energy model')
    args = parser.parse_args()
    
    app = QtWidgets.QApplication(sys.argv)
    app.setFont(QtGui.QFont("Sans", 8))
    
    window = FitREQInteractiveWindow(
        xyz_file=args.xyz,
        atom_types_file=args.atom_types,
        model_name=args.model
    )
    window.show()
    
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
