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
        self.edit_mode = 'Ring' # 'Ring', 'Atom', or 'pi'
        self.last_clicked_node = None
        self.label_mode = 'Element+Index'
        
        self.initUI()
        self.scene.sig_atom_dragged.connect(self.on_atom_dragged)
        self.scene.sig_rmb_remove.connect(self.on_atom_remove)
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
        
        # Configure axis after linking to view
        self.axis_x.axis.text_color = 'black'
        self.axis_y.axis.text_color = 'black'
        self.axis_x.axis.tick_color = 'black'
        self.axis_y.axis.tick_color = 'black'

        self.setCentralWidget(self.scene.widget)
        
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

        # --- Sidebar / ToolBar ---
        toolbar = self.addToolBar("Controls")
        
        # Mode Selection
        self.mode_combo = QtWidgets.QComboBox()
        self.mode_combo.addItems(["Ring", "Atom", "pi", "Select"])
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

        # Label Mode Selection
        self.label_combo = QtWidgets.QComboBox()
        self.label_combo.addItems(["Element+Index", "Atomic Type", "Pi Orbitals", "Z-Height", "Charge"])
        self.label_combo.currentTextChanged.connect(self.set_label_mode)
        toolbar.addWidget(QtWidgets.QLabel(" Labels: "))
        toolbar.addWidget(self.label_combo)

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
        self.scene.canvas.events.mouse_move.connect(self.on_mouse_move)
        self.scene.sig_selection_changed.connect(self.on_selection_changed)

        # Copy/paste buffer
        self.copied_atoms = None  # (enames, apos) tuple

        # Keyboard shortcuts
        self.scene.canvas.events.key_press.connect(self.on_key_press)

    def set_edit_mode(self, mode):
        self.edit_mode = mode
        print(f"Edit Mode: {mode}")
        # Enable/disable selection mode in Vispy scene
        if mode == 'Select':
            self.scene.set_selection_mode(True)
            self.statusBar().showMessage("Selection Mode: RMB drag to select | Delete: Remove | Ctrl-C: Copy | Ctrl-V: Paste | LMB: Drag selected")
        else:
            self.scene.set_selection_mode(False)
            self.statusBar().showMessage("LMB: Add/Toggle | RMB: Remove | Middle-Click: Toggle H | Scroll: Zoom")

    def set_atom_type(self, atype):
        self.cur_atom_type = atype
        print(f"Atom Type: {atype}")

    def set_label_mode(self, mode):
        self.label_mode = mode
        print(f"Label Mode: {mode}")
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
            idx = self.backend._append_atom(
                pos=list(pos.copy()),
                ename=ename,
                pin=None,
                parent=None,
                subtype=f"{ename}_sp2"
            )
            new_indices.append(idx)
        self.backend.recalc_bonds()
        # Don't call adjust_h() - it adds H atoms and shifts indices
        # Refresh view first to update scene arrays
        self.refresh_view()
        # Select the newly pasted atoms
        self.scene.set_selected_indices(new_indices)
        print(f"Pasted {len(enames)} atoms at original positions")

    def on_mouse_move(self, event):
        """Update cursor cross position on mouse move."""
        r0, rd = self.scene._ray_from_mouse(event.pos)
        p_world = self.scene._intersect_ray_plane(r0, rd, np.zeros(3), np.array([0,0,1]))
        if p_world is not None:
            self.cursor_markers.set_data(
                pos=np.array([p_world]),
                symbol='cross',
                edge_width=2,
                edge_color='red',
                face_color='transparent',
                size=10
            )

    def on_atom_remove(self, idx):
        """Remove atom at index and refresh view."""
        # Find which node_key corresponds to this atom index
        node_to_remove = None
        for nk, ia in list(self.backend.node_to_atom.items()):
            if ia == idx:
                node_to_remove = nk
                break
        if node_to_remove:
            self.backend.remove_atom(node_to_remove)
        else:
            # Atom not pinned to grid (e.g. H) - remove by index
            # Use backend's _rebuild_after_delete to keep tracking arrays consistent
            self.backend._rebuild_after_delete([idx])
            self.backend.recalc_bonds()
        self.refresh_view()

    def on_mouse_press(self, event):
        # In Select mode, let Vispy handle everything (RMB selection, LMB drag)
        if self.edit_mode == 'Select':
            return

        # If selection mode and atoms selected, skip normal LMB handling
        selected = self.scene.get_selected_indices()
        if self.edit_mode == 'Select' and selected and event.button == 1:
            # Let Vispy handle dragging of selected atoms
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

        if event.button == 1: # LMB
            self.handle_click(event.pos, action='add')
        elif event.button == 2: # RMB
            self.handle_click(event.pos, action='remove')
        elif event.button == 3: # Middle / Scroll click
            self.handle_click(event.pos, action='toggle_h')

    def on_atom_dragged(self, idx, pos):
        """Update backend atom position when dragged in the UI."""
        # Update ALL atoms, not just pinned ones
        self.backend.sys.apos[idx, 0] = pos[0]
        self.backend.sys.apos[idx, 1] = pos[1]
        self.backend.sys.apos[idx, 2] = pos[2]
        self.update_bond_labels()

    def reset_offsets(self):
        self.backend.snap_atoms_to_grid()
        self.refresh_view()

    def show_xyz(self):
        """Show current structure in a text dialog."""
        xyz_str = self.backend.get_xyz_string()
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
        if self.edit_mode == 'Ring':
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
            if self.edit_mode == 'Ring':
                self.backend.toggle_ring(q, r)
            elif self.edit_mode == 'pi':
                # Cycle pi orbitals: 0 -> 1 -> 2 -> 0
                if node_key:
                    ia = self.backend.node_to_atom.get(node_key)
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
            if self.edit_mode == 'Ring':
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
            symbol='disc', edge_width=0, size=5,
            face_color=(0.3, 0.3, 0.3, 0.8)
        )

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
            
            self.scene.set_data(pos, colors=colors, sizes=sizes, bonds=bonds_heavy)
            
            if bonds_h:
                h_segs = pos[np.array(bonds_h)].reshape(-1, 3)
                self.scene._line_set("CH-bonds", self.scene.neigh_lines, h_segs, color=(0.4, 0.4, 0.4, 0.6), width=1.0)
            else:
                self.scene.neigh_lines.set_data(np.zeros((0, 3), dtype=np.float32))

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

            # Labels based on label_mode
            if self.label_mode == 'Element+Index':
                # Show element + index (e.g. N1, C12)
                lbl_pos = []
                lbl_texts = []
                for i, e in enumerate(sys.enames):
                    if e != 'H':  # Only show for heavy atoms
                        lbl_pos.append(pos[i])
                        lbl_texts.append(f"{e}{i}")
                if lbl_pos:
                    self.scene.text_labels.text = lbl_texts
                    self.scene.text_labels.pos = np.array(lbl_pos, dtype=np.float32)
                    self.scene.text_labels.color = np.array([(0, 0, 0, 1)] * len(lbl_texts), dtype=np.float32)
                    self.scene.text_labels.visible = True
                else:
                    self.scene.text_labels.visible = False
            elif self.label_mode == 'Atomic Type':
                # Show atomic type (sp1, sp2, sp3)
                lbl_pos = []
                lbl_texts = []
                for i, subtype in enumerate(self.backend.atom_subtype):
                    if sys.enames[i] != 'H':
                        lbl_pos.append(pos[i])
                        if 'sp3' in subtype:
                            lbl_texts.append('sp3')
                        elif 'sp2' in subtype:
                            lbl_texts.append('sp2')
                        elif 'sp' in subtype:
                            lbl_texts.append('sp')
                        else:
                            lbl_texts.append(subtype)
                if lbl_pos:
                    self.scene.text_labels.text = lbl_texts
                    self.scene.text_labels.pos = np.array(lbl_pos, dtype=np.float32)
                    self.scene.text_labels.color = np.array([(0, 0, 0, 1)] * len(lbl_texts), dtype=np.float32)
                    self.scene.text_labels.visible = True
                else:
                    self.scene.text_labels.visible = False
            elif self.label_mode == 'Pi Orbitals':
                # Show number of pi orbitals
                lbl_pos = []
                lbl_texts = []
                for i in range(len(sys.enames)):
                    if i < len(self.backend.atom_subtype):
                        subtype = self.backend.atom_subtype[i]
                        if sys.enames[i] != 'H':
                            lbl_pos.append(pos[i])
                            npi = self.backend._get_npi_from_subtype(subtype)
                            lbl_texts.append(str(npi))
                if lbl_pos:
                    self.scene.text_labels.text = lbl_texts
                    self.scene.text_labels.pos = np.array(lbl_pos, dtype=np.float32)
                    self.scene.text_labels.color = np.array([(0, 0, 0, 1)] * len(lbl_texts), dtype=np.float32)
                    self.scene.text_labels.visible = True
                else:
                    self.scene.text_labels.visible = False
            elif self.label_mode == 'Z-Height':
                # Show Z coordinate
                lbl_pos = []
                lbl_texts = []
                for i, e in enumerate(sys.enames):
                    if e != 'H':
                        lbl_pos.append(pos[i])
                        lbl_texts.append(f"{pos[i, 2]:.2f}")
                if lbl_pos:
                    self.scene.text_labels.text = lbl_texts
                    self.scene.text_labels.pos = np.array(lbl_pos, dtype=np.float32)
                    self.scene.text_labels.color = np.array([(0, 0, 0, 1)] * len(lbl_texts), dtype=np.float32)
                    self.scene.text_labels.visible = True
                else:
                    self.scene.text_labels.visible = False
            elif self.label_mode == 'Charge':
                # Show charge (placeholder - not implemented yet)
                lbl_pos = []
                lbl_texts = []
                for i, e in enumerate(sys.enames):
                    if e != 'H':
                        lbl_pos.append(pos[i])
                        lbl_texts.append("0")  # Placeholder
                if lbl_pos:
                    self.scene.text_labels.text = lbl_texts
                    self.scene.text_labels.pos = np.array(lbl_pos, dtype=np.float32)
                    self.scene.text_labels.color = np.array([(0, 0, 0, 1)] * len(lbl_texts), dtype=np.float32)
                    self.scene.text_labels.visible = True
                else:
                    self.scene.text_labels.visible = False
            else:
                # Default: show bond distances
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
            sys = self.backend.sys
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
            E, forces, lvs = self.backend.run_relaxation(workdir='gui_relax')
            self.statusBar().showMessage(f"Relaxation done. E = {E:.4f} eV")
            self.refresh_view()
        except Exception as e:
            self.statusBar().showMessage(f"Relaxation FAILED: {e}")
            QtWidgets.QMessageBox.critical(self, "Relaxation Error", str(e))
            print(f"Relaxation Error: {e}")

if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    window = KekuleExplorerWindow()
    window.show()
    sys.exit(app.exec_())
