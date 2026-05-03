#!/usr/bin/env python3
import sys
import os
import numpy as np
from PyQt5 import QtWidgets, QtCore, QtGui
from vispy import scene

# Add FireCore to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from pyBall.KekuleBackend import KekuleBackend
from pyBall.AtomicSystem import AtomicSystem
from pyBall import VispyUtils as vu
from pyBall import atomicUtils as au
from pyBall import elements

class KekuleExplorerWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Kekule Structure Explorer")
        self.resize(1024, 768)

        self.backend = KekuleBackend()
        self.cur_atom_type = 'C'
        self.edit_mode = 'Ring' # 'Ring' or 'Atom'
        self.bAutoPassivate = False
        self.bAutoRecalcBonds = False
        self.cached_sys = None
        
        self.initUI()
        self.scene.sig_atom_dragged.connect(self.on_atom_dragged)
        self.idx_to_node = {}
        self.refresh_view()

    def initUI(self):
        # --- Central Widget (Vispy Scene) ---
        self.scene = vu.AtomScene(bgcolor=(0.95, 0.95, 0.95))
        
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

        self.setCentralWidget(self.scene.widget)
        
        # Add grid guide markers
        self.grid_markers = scene.visuals.Markers(parent=self.scene.view.scene)
        self.grid_markers.set_gl_state('translucent', depth_test=False)
        self.grid_markers.order = -1  # Behind everything

        # --- Sidebar / ToolBar ---
        toolbar = self.addToolBar("Controls")
        
        # Mode Selection
        self.mode_combo = QtWidgets.QComboBox()
        self.mode_combo.addItems(["Ring", "Atom"])
        self.mode_combo.currentTextChanged.connect(self.set_edit_mode)
        toolbar.addWidget(QtWidgets.QLabel(" Mode: "))
        toolbar.addWidget(self.mode_combo)

        # Atom Selection
        self.atom_combo = QtWidgets.QComboBox()
        self.atom_combo.addItems(["C", "N", "O"])
        self.atom_combo.currentTextChanged.connect(self.set_atom_type)
        toolbar.addWidget(QtWidgets.QLabel(" Atom: "))
        toolbar.addWidget(self.atom_combo)

        toolbar.addSeparator()

        # Relaxation Button
        relax_btn = QtWidgets.QPushButton("Relax (DFTB+)")
        relax_btn.clicked.connect(self.run_relaxation)
        toolbar.addWidget(relax_btn)

        # Auto-passivate checkbox
        self.passivate_cb = QtWidgets.QCheckBox("Auto-Passivate")
        self.passivate_cb.stateChanged.connect(self.toggle_passivate)
        self.passivate_cb.setChecked(self.bAutoPassivate)
        toolbar.addWidget(self.passivate_cb)

        toolbar.addSeparator()

        # Snap to grid button
        snap_btn = QtWidgets.QPushButton("Snap to Grid")
        snap_btn.clicked.connect(self.reset_offsets)
        toolbar.addWidget(snap_btn)

        toolbar.addSeparator()

        # Manual control buttons
        h_btn = QtWidgets.QPushButton("Adjust H")
        h_btn.clicked.connect(self.adjust_h)
        toolbar.addWidget(h_btn)

        bonds_btn = QtWidgets.QPushButton("Recalc Bonds")
        bonds_btn.clicked.connect(self.recalc_bonds)
        toolbar.addWidget(bonds_btn)

        toolbar.addSeparator()

        # XYZ Buttons
        show_xyz_btn = QtWidgets.QPushButton("Show XYZ")
        show_xyz_btn.clicked.connect(self.show_xyz)
        toolbar.addWidget(show_xyz_btn)

        export_xyz_btn = QtWidgets.QPushButton("Export XYZ")
        export_xyz_btn.clicked.connect(self.export_xyz)
        toolbar.addWidget(export_xyz_btn)

        # Help / Status
        self.statusBar().showMessage("LMB: Add/Toggle | RMB: Remove | Middle-Click: Toggle H | Scroll: Zoom")

        # --- Connect Mouse Events ---
        self.scene.canvas.events.mouse_press.connect(self.on_mouse_press)

    def set_edit_mode(self, mode):
        self.edit_mode = mode
        print(f"Edit Mode: {mode}")

    def set_atom_type(self, atype):
        self.cur_atom_type = atype
        print(f"Atom Type: {atype}")

    def toggle_passivate(self, state):
        self.bAutoPassivate = (state == QtCore.Qt.Checked)
        self.refresh_view()

    def on_mouse_press(self, event):
        # Guard: If Vispy's AtomScene already picked an atom for dragging, skip grid editing
        if self.scene._pick_idx >= 0:
            print(f"DEBUG: Mouse press ignored (picked atom {self.scene._pick_idx})")
            return

        if event.button == 1: # LMB
            self.handle_click(event.pos, action='add')
        elif event.button == 2: # RMB
            self.handle_click(event.pos, action='remove')
        elif event.button == 3: # Middle / Scroll click
            self.handle_click(event.pos, action='toggle_h')

    def on_atom_dragged(self, idx, pos):
        """Update backend offset when an atom is dragged in the UI."""
        if idx in self.idx_to_node:
            node_key = self.idx_to_node[idx]
            offset = pos[:2] - np.array(node_key)
            self.backend.update_node_offset(node_key, offset)
            # We don't call refresh_view here to avoid flickering during drag
            # but we update the bond labels if needed
            self.update_bond_labels()

    def reset_offsets(self):
        self.backend.reset_offsets()
        self.refresh_view()

    def show_xyz(self):
        """Show current structure in a text dialog."""
        # Use cached_sys if available and matches passivate setting
        sys = self.cached_sys if self.cached_sys else self.backend.build_system(bPassivate=self.bAutoPassivate)
        xyz_str = self.backend.get_xyz_string(sys=sys)
        
        dialog = QtWidgets.QDialog(self)
        dialog.setWindowTitle("Current XYZ Structure")
        layout = QtWidgets.QVBoxLayout(dialog)
        text_edit = QtWidgets.QPlainTextEdit()
        text_edit.setPlainText(xyz_str)
        text_edit.setReadOnly(True)
        text_edit.setMinimumSize(400, 500)
        layout.addWidget(text_edit)
        close_btn = QtWidgets.QPushButton("Close")
        close_btn.clicked.connect(dialog.accept)
        layout.addWidget(close_btn)
        dialog.exec_()

    def export_xyz(self):
        """Export current structure to an XYZ file."""
        fname, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Export XYZ", "", "XYZ Files (*.xyz)")
        if fname:
            sys = self.cached_sys if self.cached_sys else self.backend.build_system(bPassivate=self.bAutoPassivate)
            self.backend.save_xyz(fname, sys=sys)
            self.statusBar().showMessage(f"Exported to {fname}")

    def adjust_h(self):
        """Manually trigger H passivation and make it persistent."""
        self.passivate_cb.setChecked(True) # This will trigger refresh_view(bPassivate=True)

    def recalc_bonds(self):
        """Manually trigger bond recalculation and refresh view."""
        self.refresh_view(bRecalcBonds=True)

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
        if self.edit_mode == 'Ring':
            q, r = self.backend.snap_to_ring(x, y)
            node_key = None
        else:
            node_key = self.backend.snap_to_node(x, y)
            q, r = (None, None)
        
        print(f"Click at ({x:.2f}, {y:.2f}) -> Mode={self.edit_mode} Ring={(q,r)} Node={node_key} | Action: {action}")

        # 3. Modify backend
        if action == 'add':
            if self.edit_mode == 'Ring':
                self.backend.toggle_ring(q, r)
            elif node_key:
                print(f"DEBUG: Setting node {node_key} to {self.cur_atom_type}")
                self.backend.set_atom(node_key, self.cur_atom_type)
        elif action == 'remove':
            if self.edit_mode == 'Ring':
                self.backend.remove_ring(q, r)
            elif node_key:
                self.backend.remove_atom(node_key)
        elif action == 'toggle_h':
            # For toggle_h, we always snap to the nearest node regardless of mode
            nk = self.backend.snap_to_node(x, y)
            if nk:
                self.backend.toggle_h_state(nk)
        
        self.refresh_view()

    def refresh_view(self, bPassivate=None, bRecalcBonds=None):
        # Use defaults if not specified
        if bPassivate is None: bPassivate = self.bAutoPassivate
        if bRecalcBonds is None: bRecalcBonds = self.bAutoRecalcBonds

        # 0. Update Guide Grid
        guides = self.backend.get_guide_points()
        print(f"DEBUG: refresh_view guides={len(guides)} bPassivate={bPassivate}", flush=True)
        self.grid_markers.set_data(
            pos=np.column_stack([guides, np.full(len(guides), -0.1)]).astype(np.float32),
            symbol='disc', edge_width=0, size=5,
            face_color=(0.3, 0.3, 0.3, 0.8)
        )

        # 1. Build system
        sys = self.backend.build_system(bPassivate=bPassivate)
        self.cached_sys = sys # store for manual bond reuse if we want
        
        # 2. Prepare visual data
        pos = sys.apos.astype(np.float32)
        
        # Build index-to-node mapping for heavy atoms
        self.idx_to_node = {}
        node_keys = sorted(self.backend.nodes.keys())
        print(f"DEBUG: refresh_view backend_nodes={len(node_keys)} sys_atoms={len(sys.apos)} sys_enames[:10]={sys.enames[:10]}")
        for i, nk in enumerate(node_keys):
            self.idx_to_node[i] = nk
            if i < len(sys.enames) and sys.enames[i] != self.backend.nodes[nk]:
                print(f"WARNING: Element mismatch at node {nk}: Backend={self.backend.nodes[nk]} System={sys.enames[i]}")

        if pos.size == 0:
            self.scene.set_data(np.zeros((0,3)))
            return

        # Colors based on elements
        colors = []
        sizes = []
        for e in sys.enames:
            c = elements.getColor(e)
            if e == 'H':
                # Make Hydrogen darker gray for visibility on light background
                colors.append((0.4, 0.4, 0.4, 1.0))
                sizes.append(8.0)
            else:
                colors.append((c[0], c[1], c[2], 1.0))
                sizes.append(15.0)
        
        colors = np.array(colors, dtype=np.float32)
        sizes = np.array(sizes, dtype=np.float32)
        
        # Bonds
        if sys.bonds is not None:
            # 1. Heavy-Heavy Bonds (Thick)
            # 2. Heavy-H Bonds (Thin)
            # 3. H-H Bonds (Thin)
            
            # Identify bond types
            is_heavy = np.array([sys.enames[i] != 'H' for i in range(len(sys.enames))])
            bonds_heavy = []
            bonds_h = []
            
            for b in sys.bonds:
                if is_heavy[b[0]] and is_heavy[b[1]]:
                    bonds_heavy.append(b)
                else:
                    bonds_h.append(b)
            
            # We use 'bond_lines' for heavy-heavy and 'neigh_lines' for H-involved
            # AtomScene.set_data doesn't expose neigh_lines easily, so we use _line_set
            
            self.scene.set_data(pos, colors=colors, sizes=sizes, bonds=bonds_heavy)
            
            # Add H-bonds/C-H bonds using neigh_lines
            if bonds_h:
                h_segs = pos[np.array(bonds_h)].reshape(-1, 3)
                self.scene._line_set("CH-bonds", self.scene.neigh_lines, h_segs, color=(0.4, 0.4, 0.4, 0.6), width=1.0)
            else:
                self.scene.neigh_lines.set_data(np.zeros((0, 3), dtype=np.float32))

            # 3. H-bonds (Halo lines)
            hbonds = sys.find_hbonds(bPrint=False)
            if hbonds:
                hb_segs = []
                for d, h, a, dist, ang in hbonds:
                    hb_segs.append(pos[h])
                    hb_segs.append(pos[a])
                hb_segs = np.array(hb_segs, dtype=np.float32)
                self.scene._line_set("H-bonds", self.scene.halo_lines, hb_segs, color=(0.8, 0.2, 0.8, 0.5), width=1.5)
            else:
                self.scene.halo_lines.set_data(np.zeros((0, 3), dtype=np.float32))

            # 4. Bond Labels (Bond lengths)
            if len(bonds_heavy) > 0:
                bls = []
                lbl_pos = []
                lbl_texts = []
                lbl_colors = []
                for b in bonds_heavy:
                    ia, ja = b
                    p1, p2 = pos[ia], pos[ja]
                    d = np.linalg.norm(p1 - p2)
                    bls.append(d)
                    lbl_pos.append((p1 + p2) * 0.5)
                    lbl_texts.append(f"{d:.3f}")
                
                bls = np.array(bls)
                vmin, vmax = bls.min(), bls.max()
                if abs(vmax - vmin) < 1e-4:
                    f = np.zeros_like(bls)
                else:
                    f = (bls - vmin) / (vmax - vmin)
                
                # Color labels by length (Blue to Red)
                for fi in f:
                    lbl_colors.append((fi, 0.5 * (1-fi), 1.0 - fi, 1.0))
                
                self.scene.text_labels.text = lbl_texts
                self.scene.text_labels.pos = np.array(lbl_pos, dtype=np.float32)
                self.scene.text_labels.color = np.array(lbl_colors, dtype=np.float32)
                self.scene.text_labels.visible = True
            else:
                self.scene.text_labels.visible = False
        else:
            self.scene.set_data(pos, colors=colors, sizes=sizes)
            self.scene.neigh_lines.set_data(np.zeros((0, 3), dtype=np.float32))
            self.scene.halo_lines.set_data(np.zeros((0, 3), dtype=np.float32))
            self.scene.text_labels.visible = False

    def update_bond_labels(self):
        """Update bond labels (distances) based on current scene positions (e.g., during dragging)."""
        if self.scene.text_labels.visible:
            pos = self.scene._pos
            sys = self.backend.build_system(bPassivate=self.bAutoPassivate)
            if sys.bonds is not None:
                lbl_pos = []
                lbl_texts = []
                for b in sys.bonds:
                    is_heavy = (sys.enames[b[0]] != 'H' and sys.enames[b[1]] != 'H')
                    if is_heavy:
                        p1, p2 = pos[b[0]], pos[b[1]]
                        d = np.linalg.norm(p1 - p2)
                        lbl_pos.append((p1 + p2) * 0.5)
                        lbl_texts.append(f"{d:.3f}")
                
                if lbl_texts:
                    self.scene.text_labels.text = lbl_texts
                    self.scene.text_labels.pos = np.array(lbl_pos, dtype=np.float32)

    def run_relaxation(self):
        self.statusBar().showMessage("Relaxing... please wait")
        QtWidgets.QApplication.processEvents()
        
        try:
            # Build system and track original CoG to undo shift
            sys_init = self.backend.build_system(bPassivate=self.bAutoPassivate)
            cog0 = sys_init.apos.mean(axis=0) if len(sys_init.apos)>0 else np.zeros(3)
            
            E, apos_new, forces, lvs = self.backend.run_relaxation(sys_init, workdir='gui_relax')
            self.statusBar().showMessage(f"Relaxation done. E = {E:.4f} eV")
            
            # Undo centering shift if backend did it
            # run_relaxation centers atoms in cell. We want to preserve absolute coordinates.
            cog1 = apos_new.mean(axis=0) if len(apos_new)>0 else np.zeros(3)
            # Actually, the backend centers based on LVs. 
            # We just want to ensure that the difference relative to the grid nodes is correct.
            
            # Save relaxed positions back to grid offsets for persistence
            for i, nk in self.idx_to_node.items():
                if i < len(apos_new):
                    # We need the absolute position relative to the grid
                    # If the whole molecule shifted, the offset should NOT include the shift
                    # but should represent the distortion.
                    # Wait, if we undo the shift on apos_new, we are back to original coordinate system.
                    
                    # Better: calculate distortion from sys_init.apos
                    # distortion = apos_new[i] - sys_init.apos[i]
                    # This way it doesn't matter where it is in space.
                    
                    p_new = apos_new[i]
                    p_old = sys_init.apos[i]
                    distortion = p_new - p_old
                    
                    # Current offset in backend
                    off_old = self.backend.node_offsets.get(nk, np.zeros(2))
                    off_new = off_old + distortion[:2]
                    self.backend.update_node_offset(nk, off_new)
            
            self.refresh_view()
            # We'll just update the view with relaxed coords.
            
            # To update the view while keeping the grid logic, we might need a "relaxed" mode.
            # For now, let's just refresh with the new apos.
            self.scene.set_data(apos_new.astype(np.float32), bonds=sys_init.bonds)
            
        except Exception as e:
            self.statusBar().showMessage(f"Relaxation FAILED: {e}")
            QtWidgets.QMessageBox.critical(self, "Relaxation Error", str(e))
            print(f"Relaxation Error: {e}")

if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    window = KekuleExplorerWindow()
    window.show()
    sys.exit(app.exec_())
