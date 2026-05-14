import os
import sys
import copy
import json
import subprocess
from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                             QHBoxLayout, QLabel, QPushButton, QSpinBox, 
                             QComboBox, QScrollArea, QGridLayout, QMessageBox, 
                             QDoubleSpinBox, QMenu, QCheckBox)
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QColor, QPalette, QEnterEvent

import numpy as np
import pyqtgraph.opengl as gl

def parse_xyz(filepath):
    """Loads an XYZ file and extracts atoms, charges, and parameters from the comment line."""
    with open(filepath, 'r') as f:
        lines = [line.strip() for line in f.readlines() if line.strip()]
    
    if not lines:
        return None
        
    natoms = int(lines[0])
    comment = lines[1]
    
    # Default parameters
    attach_idx = -1
    shift_x = 0.0
    shift_y = 0.0
    shift_z = 0.0
    
    for token in comment.split():
        if token.startswith("attach_idx="):
            # Assuming 1-based index in XYZ file, converting to 0-based for Python
            attach_idx = int(token.split("=")[1]) - 1
        elif token.startswith("shift_x="):
            shift_x = float(token.split("=")[1])
        elif token.startswith("shift_y="):
            shift_y = float(token.split("=")[1])
        elif token.startswith("shift_z="):
            shift_z = float(token.split("=")[1])
            
    atoms = []
    for line in lines[2:2+natoms]:
        parts = line.split()
        if len(parts) >= 4:
            q = float(parts[4]) if len(parts) > 4 else 0.0
            atoms.append({
                "elem": parts[0],
                "x": float(parts[1]),
                "y": float(parts[2]),
                "z": float(parts[3]),
                "q": q
            })
            
    return {
        "name": os.path.basename(filepath).replace('.xyz', ''),
        "filepath": filepath,
        "natoms": natoms,
        "comment": comment,
        "attach_idx": attach_idx,
        "shift_x": shift_x,
        "shift_y": shift_y,
        "shift_z": shift_z,
        "atoms": atoms
    }

class ItemButton(QPushButton):
    """Button representing a grid cell in the schema."""
    def __init__(self, item_type, callback, hover_callback=None):
        super().__init__("Empty")
        self.item_type = item_type # 'monomer' or 'molecule'
        self.callback = callback
        self.hover_callback = hover_callback
        self.selected_item = None
        self.setFixedSize(80, 80)
        self.clicked.connect(self.on_click)
        self.update_style()

    def set_item(self, item):
        self.selected_item = item
        if item:
            self.setText(item["name"])
        else:
            self.setText("Empty")
        self.update_style()

    def on_click(self):
        self.callback(self)

    def enterEvent(self, event: QEnterEvent):
        if self.hover_callback:
            self.hover_callback(self.selected_item)
        super().enterEvent(event)

    def update_style(self):
        if self.selected_item is None:
            self.setStyleSheet("background-color: #e0e0e0; border: 1px dashed #aaa; border-radius: 5px;")
        else:
            # Simple color hashing based on the item name
            h = hash(self.selected_item["name"])
            r = (h & 0xFF0000) >> 16
            g = (h & 0x00FF00) >> 8
            b = (h & 0x0000FF)
            # Lighten the color so black text remains readable
            r = int(r * 0.5 + 127)
            g = int(g * 0.5 + 127)
            b = int(b * 0.5 + 127)
            
            if self.item_type == 'monomer':
                self.setStyleSheet(f"background-color: rgb({r},{g},{b}); border: 2px solid #333; font-weight: bold;")
            else:
                self.setStyleSheet(f"background-color: rgb({r},{g},{b}); border: 2px solid #333; border-radius: 40px;")

class BuilderWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Polymer System Builder")
        self.resize(1000, 700)
        
        # Library Data
        self.monomers = []
        self.molecules = []
        
        # Grid structure: [row][col] where col = [0: BB_L, 1: Mol_L, 2: Gap, 3: Mol_R, 4: BB_R]
        self.grid_items = []
        
        self.init_ui()
        self.load_libraries()
        self.load_session()

    def init_ui(self):
        main_widget = QWidget()
        main_layout = QHBoxLayout(main_widget)
        
        # --- Left Panel: Controls ---
        control_panel = QWidget()
        control_layout = QVBoxLayout(control_panel)
        control_panel.setFixedWidth(250)
        
        # Polymer Length
        length_layout = QHBoxLayout()
        length_layout.addWidget(QLabel("Polymer Length:"))
        self.length_spin = QSpinBox()
        self.length_spin.setRange(1, 100)
        self.length_spin.setValue(5)
        self.length_spin.valueChanged.connect(self.rebuild_grid)
        length_layout.addWidget(self.length_spin)
        control_layout.addLayout(length_layout)
        
        # Orientation
        control_layout.addWidget(QLabel("Complement Orientation:"))
        self.orient_combo = QComboBox()
        self.orient_combo.addItems(["Parallel", "Antiparallel"])
        control_layout.addWidget(self.orient_combo)
        
        # X Distance (Padding)
        dist_layout = QHBoxLayout()
        dist_layout.addWidget(QLabel("Strand Distance (Padding) [Å]:"))
        self.dist_spin = QDoubleSpinBox()
        self.dist_spin.setRange(0.0, 20.0)
        self.dist_spin.setValue(2.0)
        self.dist_spin.setSingleStep(0.5)
        dist_layout.addWidget(self.dist_spin)
        control_layout.addLayout(dist_layout)
        
        # Right Strand Y Offset
        y_offset_layout = QHBoxLayout()
        y_offset_layout.addWidget(QLabel("Right Strand Y Offset [Å]:"))
        self.y_offset_spin = QDoubleSpinBox()
        self.y_offset_spin.setRange(-20.0, 20.0)
        self.y_offset_spin.setValue(0.0)
        self.y_offset_spin.setSingleStep(0.5)
        y_offset_layout.addWidget(self.y_offset_spin)
        control_layout.addLayout(y_offset_layout)
        
        # Refresh Library Button
        btn_refresh = QPushButton("Reload Library")
        btn_refresh.clicked.connect(self.load_libraries)
        control_layout.addWidget(btn_refresh)
        
        # Batch Fill
        batch_layout = QVBoxLayout()
        batch_layout.addWidget(QLabel("Batch Fill Column:"))
        
        self.batch_col_combo = QComboBox()
        self.batch_col_combo.addItems(["Left Backbone", "Left Molecule", "Right Molecule", "Right Backbone"])
        self.batch_col_combo.currentIndexChanged.connect(self.update_batch_items)
        batch_layout.addWidget(self.batch_col_combo)
        
        self.batch_item_combo = QComboBox()
        batch_layout.addWidget(self.batch_item_combo)
        
        btn_batch_fill = QPushButton("Fill Column")
        btn_batch_fill.clicked.connect(self.batch_fill)
        batch_layout.addWidget(btn_batch_fill)
        
        control_layout.addLayout(batch_layout)
        
        control_layout.addStretch()
        
        # Export Button & Jmol Option
        btn_export = QPushButton("Export XYZ")
        btn_export.setStyleSheet("background-color: #4CAF50; color: white; font-weight: bold; padding: 10px;")
        btn_export.clicked.connect(self.export_system)
        control_layout.addWidget(btn_export)
        
        self.chk_jmol = QCheckBox("Show structure in Jmol")
        self.chk_jmol.setChecked(True)
        control_layout.addWidget(self.chk_jmol)
        
        # --- 3D Preview Widget ---
        control_layout.addWidget(QLabel("3D Preview:"))
        self.preview_canvas = gl.GLViewWidget()
        self.preview_canvas.setMinimumHeight(250)
        self.preview_canvas.opts['distance'] = 15
        control_layout.addWidget(self.preview_canvas)
        
        # --- Right Panel: Canvas & PBC ---
        right_panel = QWidget()
        right_layout = QVBoxLayout(right_panel)
        
        self.canvas_area = QScrollArea()
        self.canvas_area.setWidgetResizable(True)
        self.canvas_widget = QWidget()
        self.canvas_layout = QGridLayout(self.canvas_widget)
        self.canvas_layout.setAlignment(Qt.AlignmentFlag.AlignTop | Qt.AlignmentFlag.AlignHCenter)
        self.canvas_area.setWidget(self.canvas_widget)
        right_layout.addWidget(self.canvas_area)
        
        # PBC Controls (Under the canvas)
        pbc_widget = QWidget()
        pbc_layout = QHBoxLayout(pbc_widget)
        pbc_layout.addWidget(QLabel("PBC Mode:"))
        
        self.chk_pbc_1d = QCheckBox("1D (Y)")
        self.chk_pbc_3d = QCheckBox("3D (X, Y, Z)")
        self.chk_pbc_1d.toggled.connect(lambda checked: self.chk_pbc_3d.setChecked(False) if checked else None)
        self.chk_pbc_3d.toggled.connect(lambda checked: self.chk_pbc_1d.setChecked(False) if checked else None)
        pbc_layout.addWidget(self.chk_pbc_1d)
        pbc_layout.addWidget(self.chk_pbc_3d)
        
        pbc_layout.addWidget(QLabel(" Lx [Å]:"))
        self.lx_spin = QDoubleSpinBox()
        self.lx_spin.setRange(0.0, 500.0)
        self.lx_spin.setValue(15.0)
        pbc_layout.addWidget(self.lx_spin)
        
        pbc_layout.addWidget(QLabel(" Lz [Å]:"))
        self.lz_spin = QDoubleSpinBox()
        self.lz_spin.setRange(0.0, 500.0)
        self.lz_spin.setValue(15.0)
        pbc_layout.addWidget(self.lz_spin)
        
        pbc_layout.addStretch()
        right_layout.addWidget(pbc_widget)
        
        main_layout.addWidget(control_panel)
        main_layout.addWidget(right_panel)
        self.setCentralWidget(main_widget)
        
        self.rebuild_grid()

    def closeEvent(self, event):
        """Prompt to save session state when the window is closed."""
        reply = QMessageBox.question(
            self, 'Save Session?', 
            'Do you want to save the current work for the next session?',
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No | QMessageBox.StandardButton.Cancel
        )
        
        if reply == QMessageBox.StandardButton.Yes:
            self.save_session()
            event.accept()
        elif reply == QMessageBox.StandardButton.No:
            # Remove the session file so it starts fresh next time
            session_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'session.json')
            if os.path.exists(session_file):
                try: os.remove(session_file)
                except Exception: pass
            event.accept()
        else:
            event.ignore()

    def save_session(self):
        session_data = {
            "length": self.length_spin.value(),
            "orient": self.orient_combo.currentText(),
            "dist": self.dist_spin.value(),
            "y_offset": self.y_offset_spin.value(),
            "pbc_1d": self.chk_pbc_1d.isChecked(),
            "pbc_3d": self.chk_pbc_3d.isChecked(),
            "lx": self.lx_spin.value(),
            "lz": self.lz_spin.value(),
            "grid": []
        }
        for row in self.grid_items:
            row_data = []
            for btn in row:
                if btn is None:
                    row_data.append(None)
                elif btn.selected_item:
                    row_data.append({"name": btn.selected_item["name"], "type": btn.item_type})
                else:
                    row_data.append(None)
            session_data["grid"].append(row_data)
        
        try:
            session_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'session.json')
            with open(session_file, 'w') as f:
                json.dump(session_data, f)
        except Exception as e:
            print("Failed to save session:", e)

    def load_session(self):
        session_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'session.json')
        if not os.path.exists(session_file):
            return
        try:
            with open(session_file, 'r') as f:
                session_data = json.load(f)
            
            # Block signals to prevent the grid from redrawing empty after the first change
            self.length_spin.blockSignals(True)
            self.length_spin.setValue(session_data.get("length", 5))
            self.length_spin.blockSignals(False)
            
            self.orient_combo.setCurrentText(session_data.get("orient", "Parallel"))
            self.dist_spin.setValue(session_data.get("dist", 2.0))
            self.y_offset_spin.setValue(session_data.get("y_offset", 0.0))
            
            self.chk_pbc_1d.setChecked(session_data.get("pbc_1d", False))
            self.chk_pbc_3d.setChecked(session_data.get("pbc_3d", False))
            self.lx_spin.setValue(session_data.get("lx", 15.0))
            self.lz_spin.setValue(session_data.get("lz", 15.0))
            
            self.rebuild_grid()
            
            grid_data = session_data.get("grid", [])
            for r, row_data in enumerate(grid_data):
                if r >= len(self.grid_items): break
                for c, cell_data in enumerate(row_data):
                    if cell_data and self.grid_items[r][c] is not None:
                        name = cell_data["name"]
                        item_type = cell_data["type"]
                        collection = self.monomers if item_type == "monomer" else self.molecules
                        found = next((item for item in collection if item["name"] == name), None)
                        if found:
                            self.grid_items[r][c].set_item(found)
        except Exception as e:
            print("Failed to load session:", e)

    def load_libraries(self):
        """Loads monomers and molecules from their respective directories."""
        self.monomers = []
        self.molecules = []
        
        base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        monomer_dir = os.path.join(base_dir, "monomers")
        molecule_dir = os.path.join(base_dir, "molecules")
        
        os.makedirs(monomer_dir, exist_ok=True)
        os.makedirs(molecule_dir, exist_ok=True)
        
        for f in os.listdir(monomer_dir):
            if f.endswith('.xyz'):
                data = parse_xyz(os.path.join(monomer_dir, f))
                if data: self.monomers.append(data)
                
        for f in os.listdir(molecule_dir):
            if f.endswith('.xyz'):
                data = parse_xyz(os.path.join(molecule_dir, f))
                if data: self.molecules.append(data)
                
        print(f"Loaded monomers: {len(self.monomers)}")
        print(f"Loaded molecules: {len(self.molecules)}")
        
        self.update_batch_items()

    def update_batch_items(self):
        if not hasattr(self, 'batch_item_combo'): return
        self.batch_item_combo.clear()
        col_name = self.batch_col_combo.currentText()
        items = self.monomers if "Backbone" in col_name else self.molecules
        self.batch_item_combo.addItem("Remove (Empty)", None)
        for item in items:
            self.batch_item_combo.addItem(item["name"], item)

    def batch_fill(self):
        col_name = self.batch_col_combo.currentText()
        item = self.batch_item_combo.currentData()
        
        col_idx = 0
        if col_name == "Left Backbone": col_idx = 0
        elif col_name == "Left Molecule": col_idx = 1
        elif col_name == "Right Molecule": col_idx = 3
        elif col_name == "Right Backbone": col_idx = 4
        
        for row in self.grid_items:
            btn = row[col_idx]
            if btn:
                btn.set_item(item)

    def rebuild_grid(self):
        """Redraws the grid on the canvas according to the requested length."""
        # Clear layout
        for i in reversed(range(self.canvas_layout.count())): 
            self.canvas_layout.itemAt(i).widget().setParent(None)
            
        self.grid_items = []
        n_rows = self.length_spin.value()
        
        # Headers
        self.canvas_layout.addWidget(QLabel("Left Backbone"), 0, 0, alignment=Qt.AlignmentFlag.AlignCenter)
        self.canvas_layout.addWidget(QLabel("Left Molecule"), 0, 1, alignment=Qt.AlignmentFlag.AlignCenter)
        self.canvas_layout.addWidget(QLabel("Gap"), 0, 2, alignment=Qt.AlignmentFlag.AlignCenter)
        self.canvas_layout.addWidget(QLabel("Right Molecule"), 0, 3, alignment=Qt.AlignmentFlag.AlignCenter)
        self.canvas_layout.addWidget(QLabel("Right Backbone"), 0, 4, alignment=Qt.AlignmentFlag.AlignCenter)
        
        for r in range(n_rows):
            row_buttons = []
            
            # L_BB (Col 0)
            btn_l_bb = ItemButton('monomer', self.show_selection_menu, self.update_preview)
            self.canvas_layout.addWidget(btn_l_bb, r+1, 0)
            row_buttons.append(btn_l_bb)
            
            # L_Mol (Col 1)
            btn_l_mol = ItemButton('molecule', self.show_selection_menu, self.update_preview)
            self.canvas_layout.addWidget(btn_l_mol, r+1, 1)
            row_buttons.append(btn_l_mol)
            
            # Gap (Col 2)
            line = QWidget()
            line.setFixedSize(50, 2)
            line.setStyleSheet("background-color: #ccc;")
            self.canvas_layout.addWidget(line, r+1, 2, alignment=Qt.AlignmentFlag.AlignCenter)
            row_buttons.append(None) # Placeholder to keep column indices correct
            
            # R_Mol (Col 3)
            btn_r_mol = ItemButton('molecule', self.show_selection_menu, self.update_preview)
            self.canvas_layout.addWidget(btn_r_mol, r+1, 3)
            row_buttons.append(btn_r_mol)
            
            # R_BB (Col 4)
            btn_r_bb = ItemButton('monomer', self.show_selection_menu, self.update_preview)
            self.canvas_layout.addWidget(btn_r_bb, r+1, 4)
            row_buttons.append(btn_r_bb)
            
            self.grid_items.append(row_buttons)

    def show_selection_menu(self, btn):
        """Shows a popup menu to select a structure from the library."""
        menu = QMenu(self)
        
        # Option to clear the cell
        action_clear = menu.addAction("Remove (Empty)")
        action_clear.triggered.connect(lambda: btn.set_item(None))
        menu.addSeparator()
        
        items = self.monomers if btn.item_type == 'monomer' else self.molecules
        
        if not items:
            menu.addAction("Library is empty!").setEnabled(False)
        else:
            for item in items:
                action = menu.addAction(item["name"])
                action.setData(item) # Store data for hover event
                # Closure trick to capture iterator
                action.triggered.connect(lambda checked, i=item: btn.set_item(i))
                
        menu.hovered.connect(self.on_menu_hover)
        menu.exec(btn.mapToGlobal(btn.rect().bottomLeft()))

    def on_menu_hover(self, action):
        item = action.data()
        if item:
            self.update_preview(item)

    def update_preview(self, item):
        """Updates the pyqtgraph 3D widget with the selected molecule."""
        self.preview_canvas.clear()
        
        if not item:
            return
            
        pos = []
        colors = []
        sizes = []
        
        elems = [a['elem'] for a in item['atoms']]
        for a in item['atoms']:
            pos.append([a['x'], a['y'], a['z']])
            e_up = a['elem'].upper()
            if e_up == 'C':   c, s = (0.5, 0.5, 0.5, 1.0), 20
            elif e_up == 'O': c, s = (0.9, 0.1, 0.1, 1.0), 20
            elif e_up == 'N': c, s = (0.1, 0.2, 0.9, 1.0), 20
            elif e_up == 'H': c, s = (0.8, 0.8, 0.8, 1.0), 10
            elif e_up == 'P': c, s = (1.0, 0.6, 0.1, 1.0), 25
            elif e_up == 'S': c, s = (0.8, 0.8, 0.1, 1.0), 25
            elif e_up == 'SI':c, s = (0.6, 0.8, 0.6, 1.0), 25
            else:             c, s = (0.1, 0.9, 0.1, 1.0), 20
            colors.append(c)
            sizes.append(s)
            
        if not pos:
            return
            
        pos = np.array(pos)
        colors = np.array(colors)
        sizes = np.array(sizes)
        
        scatter = gl.GLScatterPlotItem(pos=pos, color=colors, size=sizes, pxMode=True)
        self.preview_canvas.addItem(scatter)
        
        # Lines for bonds using covalent radii heuristic
        radii = {'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'S': 1.05, 'P': 1.07, 'SI': 1.11}
        bond_lines = []
        for i in range(len(pos)):
            for j in range(i+1, len(pos)):
                e1 = elems[i].upper()
                e2 = elems[j].upper()
                if e1 == 'H' and e2 == 'H': 
                    continue # No H-H bonds
                    
                dist = np.linalg.norm(pos[i] - pos[j])
                r1 = radii.get(e1, 0.7)
                r2 = radii.get(e2, 0.7)
                
                if dist < (r1 + r2 + 0.3): # covalent bond threshold
                    bond_lines.append(pos[i])
                    bond_lines.append(pos[j])
        
        if bond_lines:
            line_pos = np.array(bond_lines)
            lines = gl.GLLinePlotItem(pos=line_pos, color=(0.5, 0.5, 0.5, 1.0), width=2.0, mode='lines')
            self.preview_canvas.addItem(lines)

    def _build_strand(self, bb_col, mol_col, invert_x=False):
        """Internal function to build a single strand along the Y axis, oriented in X."""
        strand_atoms = []
        current_y = 0.0
        
        for r, row in enumerate(self.grid_items):
            bb_item = row[bb_col].selected_item
            mol_item = row[mol_col].selected_item
            
            if bb_item is None:
                continue # Empty spot in backbone
                
            bb_atoms = copy.deepcopy(bb_item["atoms"])
            bb_att_idx = bb_item["attach_idx"]
            
            # 1. Add backbone with Y shift
            for a in bb_atoms:
                a['y'] += current_y
                if invert_x: a['x'] = -a['x']
                
            current_y += bb_item.get("shift_y", 0.0)
            
            # 2. Attach molecule, if it exists
            if mol_item:
                mol_atoms = copy.deepcopy(mol_item["atoms"])
                mol_att_idx = mol_item.get("attach_idx", -1)
                
                if bb_att_idx < 0 or bb_att_idx >= len(bb_atoms):
                    self.warnings.append(f"Monomer '{bb_item['name']}' has no valid attach_idx!")
                    strand_atoms.extend(bb_atoms)
                    strand_atoms.extend(mol_atoms)
                elif mol_att_idx < 0 or mol_att_idx >= len(mol_atoms):
                    self.warnings.append(f"Molecule '{mol_item['name']}' has no valid attach_idx!")
                    strand_atoms.extend(bb_atoms)
                    strand_atoms.extend(mol_atoms)
                else:
                    bb_att_atom = bb_atoms[bb_att_idx]
                    
                    # Mirror molecule FIRST if invert_x is True to calculate dx correctly
                    if invert_x:
                        for a in mol_atoms:
                            a['x'] = -a['x']
                            
                    mol_att_atom = mol_atoms[mol_att_idx]
                    
                    # Overlap attach atoms exactly (no additional shifts)
                    dx = bb_att_atom['x'] - mol_att_atom['x']
                    dy = bb_att_atom['y'] - mol_att_atom['y']
                    dz = bb_att_atom['z'] - mol_att_atom['z']
                    
                    for a in mol_atoms:
                        a['x'] += dx
                        a['y'] += dy
                        a['z'] += dz
                        
                    # Determine which attach atom to keep (prefer heavy atom)
                    bb_e = bb_att_atom['elem'].upper()
                    mol_e = mol_att_atom['elem'].upper()
                    
                    if mol_e != 'H' and bb_e == 'H':
                        kept_atom = mol_att_atom
                        drop_bb_idx = bb_att_idx
                        drop_mol_idx = -1 # keep mol_att_atom
                    else:
                        kept_atom = bb_att_atom
                        drop_bb_idx = -1 
                        drop_mol_idx = mol_att_idx
                        
                    merged_atoms = []
                    for i, a in enumerate(mol_atoms):
                        if i != drop_mol_idx:
                            merged_atoms.append(a)
                    for i, a in enumerate(bb_atoms):
                        if i != drop_bb_idx:
                            merged_atoms.append(a)
                        
                    strand_atoms.extend(merged_atoms)
            else:
                # No molecule attached - add full backbone (including its free hydrogen)
                strand_atoms.extend(bb_atoms)
                
        return strand_atoms

    def export_system(self):
        """Builds the system and exports it to a file."""
        self.warnings = []
        # 1. Build left strand (Molecules face +X)
        left_strand = self._build_strand(bb_col=0, mol_col=1, invert_x=False)
        
        # 2. Build right strand (Molecules originally face +X, so after invert_x=True they face -X)
        right_strand = self._build_strand(bb_col=4, mol_col=3, invert_x=True)
        
        if not left_strand and not right_strand:
            QMessageBox.warning(self, "Error", "System is empty, nothing to export.")
            return
            
        # 3. Orient right strand (Z-rotation for Antiparallel configuration)
        if self.orient_combo.currentText() == "Antiparallel":
            # Rotate around Z: x -> x (already inverted in build), y -> -y
            for a in right_strand:
                a['y'] = -a['y']
                
            # Must move the right strand up in Y axis so it aligns
            if right_strand and left_strand:
                min_y_right = min([a['y'] for a in right_strand])
                min_y_left = min([a['y'] for a in left_strand])
                offset_y = min_y_left - min_y_right
                for a in right_strand:
                    a['y'] += offset_y

        # 4. Calculate extreme points and shift in X axis
        if left_strand and right_strand:
            max_x_left = max([a['x'] for a in left_strand])
            min_x_right = min([a['x'] for a in right_strand])
            
            # Shift right_strand so its minimum X matches left strand's max X + padding
            padding = self.dist_spin.value()
            offset_x = max_x_left - min_x_right + padding
            
            for a in right_strand:
                a['x'] += offset_x
                
        # Apply manual Y offset to the right strand if defined
        manual_y_offset = self.y_offset_spin.value()
        if manual_y_offset != 0.0:
            for a in right_strand:
                a['y'] += manual_y_offset
                
        if self.warnings:
            QMessageBox.warning(self, "Missing Attachment Indices", 
                "Some components could not be connected due to missing or invalid `attach_idx` "
                "in their comment line:\n\n" + "\n".join(set(self.warnings)) + 
                "\n\nMolecule atoms were still copied to the system (without alignment).")

        # 5. Global valency check and removal of excess hydrogens across the entire system
        all_atoms = left_strand + right_strand
        
        atoms_to_delete = set()
        valencies = {'C': 4, 'N': 3, 'O': 2, 'S': 2, 'P': 5, 'SI': 4}
        heavy_valency_warnings = []
        
        # Calculate cell size for periodic boundaries
        base_shift_y = 0.0
        for row in self.grid_items:
            if row[0].selected_item:
                base_shift_y = row[0].selected_item.get("shift_y", 0.0)
                break
            elif row[4].selected_item:
                base_shift_y = row[4].selected_item.get("shift_y", 0.0)
                break
                
        ly = abs(self.length_spin.value() * base_shift_y)
        pbc_enabled = self.chk_pbc_1d.isChecked() or self.chk_pbc_3d.isChecked()
        
        for i, a in enumerate(all_atoms):
            if i in atoms_to_delete: continue
            if a['elem'].upper() == 'H': continue
            
            # Find all bonds to this heavy atom (considering Y-axis PBC)
            h_bonds = []
            heavy_bonds = 0
            for j, b in enumerate(all_atoms):
                if i == j: continue
                if j in atoms_to_delete: continue
                
                dx = a['x'] - b['x']
                dy = a['y'] - b['y']
                dz = a['z'] - b['z']
                
                if pbc_enabled and ly > 0:
                    dy = dy - ly * round(dy / ly)
                    
                dist = (dx**2 + dy**2 + dz**2)**0.5
                if dist < 1.8: # Tolerance for covalent bond
                    if b['elem'].upper() == 'H': h_bonds.append(j)
                    else: heavy_bonds += 1
                        
            max_bonds = valencies.get(a['elem'].upper(), 4)
            if heavy_bonds > max_bonds:
                heavy_valency_warnings.append(f"Atom {a['elem']} has {heavy_bonds} heavy bonds (exceeds valency {max_bonds})!")
                
            allowed_h = max(0, max_bonds - heavy_bonds)
            excess_h = len(h_bonds) - allowed_h
            
            if excess_h > 0:
                # Remove those hydrogens that collide (overlap) most with other atoms
                h_clash_scores = []
                for h_idx in h_bonds:
                    hx, hy, hz = all_atoms[h_idx]['x'], all_atoms[h_idx]['y'], all_atoms[h_idx]['z']
                    min_dist_to_other = float('inf')
                    for k, c in enumerate(all_atoms):
                        if k == i or k == h_idx or c['elem'].upper() == 'H': continue
                        
                        dx_c, dy_c, dz_c = hx - c['x'], hy - c['y'], hz - c['z']
                        if pbc_enabled and ly > 0:
                            dy_c = dy_c - ly * round(dy_c / ly)
                            
                        dist_other = (dx_c**2 + dy_c**2 + dz_c**2)**0.5
                        if dist_other < min_dist_to_other:
                            min_dist_to_other = dist_other
                            
                    h_clash_scores.append((min_dist_to_other, h_idx))
                
                h_clash_scores.sort() # The shortest distances are the colliding ones (boundary crossing)
                for k in range(min(excess_h, len(h_clash_scores))):
                    atoms_to_delete.add(h_clash_scores[k][1])
                    
        all_atoms = [a for i, a in enumerate(all_atoms) if i not in atoms_to_delete]

        if heavy_valency_warnings:
            QMessageBox.warning(self, "Critical Valency Error!", 
                "Some heavy atoms overlap after connection and have more bonds than the C++ backend (MMFF) can handle. "
                "The calculation will likely crash (NaNs).\n\n" + "\n".join(set(heavy_valency_warnings)))

        # 6. Final output to file
        n_hydrogens = sum([1 for a in all_atoms if a['elem'].upper() == 'H'])
        
        if pbc_enabled:
            lx = self.lx_spin.value()
            lz = self.lz_spin.value()
            
            # Shift atoms to the center of the defined cell (Lx, Lz)
            min_x = min([a['x'] for a in all_atoms])
            max_x = max([a['x'] for a in all_atoms])
            min_z = min([a['z'] for a in all_atoms])
            max_z = max([a['z'] for a in all_atoms])
            
            shift_to_center_x = (lx / 2.0) - ((min_x + max_x) / 2.0)
            shift_to_center_z = (lz / 2.0) - ((min_z + max_z) / 2.0)
            
            for a in all_atoms:
                a['x'] += shift_to_center_x
                a['z'] += shift_to_center_z
                
            mode_str = "1D PBC" if self.chk_pbc_1d.isChecked() else "3D PBC"
            comment_line = f"lvs {lx:.5f} 0.0 0.0    0.0 {ly:.5f} 0.0    0.0 0.0 {lz:.5f} | Generated by GUI Builder | {mode_str}"
        else:
            comment_line = f"Generated by GUI Builder | Total Hydrogens: {n_hydrogens}"
        
        output_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "output")
        os.makedirs(output_dir, exist_ok=True)
        out_file = os.path.join(output_dir, "generated_system.xyz")
        
        with open(out_file, 'w') as f:
            f.write(f"{len(all_atoms)}\n")
            f.write(f"{comment_line}\n")
            for a in all_atoms:
                f.write(f"{a['elem']:<4} {a['x']:10.5f} {a['y']:10.5f} {a['z']:10.5f} {a['q']:10.5f}\n")
                
        QMessageBox.information(self, "Success", 
            f"System saved to:\n{out_file}\n\nTotal atoms: {len(all_atoms)}\nTotal hydrogens (H): {n_hydrogens}")
            
        if self.chk_jmol.isChecked():
            try:
                # Run Jmol in the background so it doesn't block the GUI
                subprocess.Popen(["jmol.sh", out_file])
            except Exception as e:
                QMessageBox.warning(self, "Jmol Error", f"Failed to launch Jmol (jmol.sh):\n{e}")

if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = BuilderWindow()
    window.show()
    sys.exit(app.exec())