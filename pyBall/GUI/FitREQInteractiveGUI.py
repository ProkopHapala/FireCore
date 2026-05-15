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
        
        left_layout.addStretch()
        
        # --- Right panel: matplotlib canvas (with blitting for speed) ---
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
        
        # Add panels to main layout
        main_layout.addWidget(left_panel)
        main_layout.addWidget(self.canvas)
        
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
