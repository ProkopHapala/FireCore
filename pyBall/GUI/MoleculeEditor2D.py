#!/usr/bin/env python3
"""Planar molecule editor for FireCore.

Usage examples:

```
python -m pyBall.GUI.MoleculeEditor2D
python -m pyBall.GUI.MoleculeEditor2D tests/tAttach/backbones/porphirin.mol2
```

The editor embeds a Matplotlib canvas inside a PyQt5 UI and operates directly on
`AtomicSystem` instances.
"""

import sys
from math import sqrt
from pathlib import Path

import numpy as np
from PyQt5.QtCore import QSize
from PyQt5.QtWidgets import (
    QApplication,
    QMainWindow,
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QComboBox,
    QFileDialog,
    QListWidget,
    QListWidgetItem,
    QToolButton,
    QButtonGroup,
    QFrame,
)
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.collections import LineCollection

if __package__ in {None, ""}:
    package_root = Path(__file__).resolve().parents[1]
    sys.path.append(str(package_root))
    from AtomicSystem import AtomicSystem
    import atomicUtils as au
    import elements
else:
    from ..AtomicSystem import AtomicSystem
    from .. import atomicUtils as au
    from .. import elements

ELEMENT_DICT = elements.ELEMENT_DICT
INDEX_Z = elements.index_Z
INDEX_Rvdw = elements.index_Rvdw
INDEX_VALENCE = elements.index_val_elec

DEFAULT_ELEMENTS = ["H", "C", "N", "O", "F"]
MODE_DRAW = "draw"
MODE_MOVE = "move"


def canonical_element(name):
    if not name:
        return "C"
    if name in ELEMENT_DICT:
        return name
    cleaned = name.replace(".", "_")
    for token in cleaned.split("_"):
        if not token:
            continue
        token_cap = token[0].upper() + token[1:].lower()
        if token_cap in ELEMENT_DICT:
            return token_cap
    letters = ''.join(ch for ch in name if ch.isalpha())
    if letters:
        token_cap = letters[0].upper() + letters[1:].lower()
        if token_cap in ELEMENT_DICT:
            return token_cap
        upper = letters.upper()
        for width in (1, 2):
            sym = upper[:width]
            if not sym:
                continue
            formatted = sym[0].upper() + sym[1:].lower()
            if formatted in ELEMENT_DICT:
                return formatted
    return "C"


class MoleculeDocument:
    def __init__(self, system=None):
        if system is None:
            system = AtomicSystem(
                apos=np.zeros((0, 3), dtype=np.float64),
                enames=[],
                atypes=np.zeros(0, dtype=np.int32),
                qs=np.zeros(0, dtype=np.float64),
                Rs=np.zeros(0, dtype=np.float64),
                bonds=np.zeros((0, 2), dtype=np.int32),
                bPreinit=False,
            )
        if system.apos is None:
            system.apos = np.zeros((0, 3), dtype=np.float64)
        if system.enames is None:
            system.enames = []
        system.enames = [canonical_element(str(e)) for e in system.enames]
        if system.atypes is None:
            system.atypes = np.zeros(len(system.enames), dtype=np.int32)
        if system.qs is None:
            system.qs = np.zeros(len(system.enames), dtype=np.float64)
        if system.Rs is None:
            system.Rs = np.zeros(len(system.enames), dtype=np.float64)
        if system.bonds is None:
            system.bonds = np.zeros((0, 2), dtype=np.int32)
        self.system = system
        self._ensure_mutable_lists()
        self._ensure_numeric_arrays()
        self._ensure_aux_labels()
        if system.bonds.shape[0] == 0:
            self.bond_orders = np.zeros(0, dtype=np.int8)
        else:
            self.bond_orders = np.ones(system.bonds.shape[0], dtype=np.int8)
        self.selected_atoms = set()
        self.selected_bonds = set()

    @property
    def positions(self):
        return self.system.apos

    def set_positions(self, apos):
        self.system.apos = apos

    def atom_count(self):
        return self.system.apos.shape[0]

    def bond_count(self):
        return self.system.bonds.shape[0]

    def ensure_bond_orders(self):
        nb = self.bond_count()
        if self.bond_orders.shape[0] != nb:
            if nb == 0:
                self.bond_orders = np.zeros(0, dtype=np.int8)
            else:
                self.bond_orders = np.ones(nb, dtype=np.int8)

    def _ensure_mutable_lists(self):
        enames = self.system.enames
        if enames is None:
            self.system.enames = []
            return
        if isinstance(enames, np.ndarray):
            self.system.enames = [canonical_element(str(e)) for e in enames]
        else:
            self.system.enames = [canonical_element(str(e)) for e in list(enames)]

    def _ensure_numeric_arrays(self):
        if self.system.apos is None:
            self.system.apos = np.zeros((0, 3), dtype=np.float64)
        else:
            self.system.apos = np.asarray(self.system.apos, dtype=np.float64)
        if self.system.atypes is None:
            self.system.atypes = np.zeros(len(self.system.enames), dtype=np.int32)
        else:
            self.system.atypes = np.asarray(self.system.atypes, dtype=np.int32)
        if self.system.qs is None:
            self.system.qs = np.zeros(len(self.system.enames), dtype=np.float64)
        else:
            self.system.qs = np.asarray(self.system.qs, dtype=np.float64)
        if self.system.Rs is None:
            self.system.Rs = np.zeros(len(self.system.enames), dtype=np.float64)
        else:
            self.system.Rs = np.asarray(self.system.Rs, dtype=np.float64)
        if self.system.bonds is None:
            self.system.bonds = np.zeros((0, 2), dtype=np.int32)
        else:
            self.system.bonds = np.asarray(self.system.bonds, dtype=np.int32)

    def _ensure_aux_labels(self):
        aux = getattr(self.system, "aux_labels", None)
        n = self.atom_count()
        if aux is None:
            self.system.aux_labels = [str(i) for i in range(n)]
            return
        if isinstance(aux, np.ndarray):
            aux = aux.tolist()
        else:
            aux = list(aux)
        if len(aux) < n:
            start = len(aux)
            aux.extend(str(i) for i in range(start, n))
        elif len(aux) > n:
            aux = aux[:n]
        self.system.aux_labels = aux

    def auto_bonds(self, Rcut=3.0, RvdwCut=0.6):
        atypes = self.system.atypes
        natoms = self.atom_count()
        if atypes is None or len(atypes) != natoms:
            atypes = np.array(
                [ELEMENT_DICT[canonical_element(e)][INDEX_Z] for e in self.system.enames],
                dtype=np.int32,
            )
            self.system.atypes = atypes
        if natoms == 0:
            self.system.bonds = np.zeros((0, 2), dtype=np.int32)
            self.bond_orders = np.zeros(0, dtype=np.int8)
            return
        bonds, _ = au.findBondsNP(self.system.apos, atypes=atypes, Rcut=Rcut, RvdwCut=RvdwCut)
        self.system.bonds = np.asarray(bonds, dtype=np.int32)
        self.bond_orders = np.ones(len(bonds), dtype=np.int8)
        self.system.neighs()
        self._ensure_aux_labels()

    def add_atom(self, xy, element, z_height=0.0):
        self._ensure_mutable_lists()
        self._ensure_numeric_arrays()
        self._ensure_aux_labels()
        pos = np.array([xy[0], xy[1], z_height], dtype=np.float64)
        apos = self.system.apos
        if apos.size == 0:
            self.system.apos = pos.reshape(1, 3)
        else:
            self.system.apos = np.vstack([apos, pos])
        element = canonical_element(element)
        self.system.enames.append(element)
        atomic_number = ELEMENT_DICT[element][INDEX_Z]
        if self.system.atypes.size == 0:
            self.system.atypes = np.array([atomic_number], dtype=np.int32)
        else:
            self.system.atypes = np.append(self.system.atypes, atomic_number)
        valence = ELEMENT_DICT[element][INDEX_VALENCE]
        radius = ELEMENT_DICT[element][INDEX_Rvdw]
        if self.system.qs.size == 0:
            self.system.qs = np.array([valence], dtype=np.float64)
        else:
            self.system.qs = np.append(self.system.qs, valence)
        if self.system.Rs.size == 0:
            self.system.Rs = np.array([radius], dtype=np.float64)
        else:
            self.system.Rs = np.append(self.system.Rs, radius)
        self._ensure_aux_labels()
        self.system.aux_labels[-1] = str(len(self.system.aux_labels) - 1)
        self.ensure_bond_orders()
        return self.atom_count() - 1

    def delete_atom(self, index):
        if index < 0 or index >= self.atom_count():
            return
        mask = (self.system.bonds[:, 0] != index) & (self.system.bonds[:, 1] != index)
        kept_bonds = self.system.bonds[mask]
        kept_orders = self.bond_orders[mask] if self.bond_orders.shape[0] == mask.shape[0] else np.ones(len(kept_bonds), dtype=np.int8)
        # reindex bonds after deletion
        def remap(i):
            if i < index:
                return i
            return i - 1
        if kept_bonds.size:
            reindexed = np.vectorize(remap)(kept_bonds)
            self.system.bonds = reindexed.astype(np.int32)
        else:
            self.system.bonds = kept_bonds
        self.bond_orders = kept_orders.astype(np.int8) if kept_orders.size else np.zeros(0, dtype=np.int8)
        self.system.delete_atoms([index])
        self.system.neighs()
        self._ensure_mutable_lists()
        self._ensure_numeric_arrays()
        self._ensure_aux_labels()
        self.selected_atoms.clear()
        self.selected_bonds.clear()

    def find_bond_index(self, i, j):
        if i > j:
            i, j = j, i
        bonds = self.system.bonds
        if bonds.shape[0] == 0:
            return -1
        hits = np.where((bonds[:, 0] == i) & (bonds[:, 1] == j))[0]
        if hits.size == 0:
            return -1
        return int(hits[0])

    def set_bond(self, i, j, order):
        if i == j:
            return
        if i > j:
            i, j = j, i
        idx = self.find_bond_index(i, j)
        if idx >= 0:
            self.bond_orders[idx] = order
            return idx
        bond = np.array([[i, j]], dtype=np.int32)
        if self.system.bonds.size == 0:
            self.system.bonds = bond
            self.bond_orders = np.array([order], dtype=np.int8)
        else:
            self.system.bonds = np.vstack([self.system.bonds, bond])
            self.bond_orders = np.append(self.bond_orders, order).astype(np.int8)
        self.system.neighs()
        self._ensure_numeric_arrays()
        self._ensure_aux_labels()
        return self.bond_count() - 1

    def delete_bond(self, i, j):
        idx = self.find_bond_index(i, j)
        if idx < 0:
            return
        keep = np.ones(self.bond_count(), dtype=bool)
        keep[idx] = False
        self.system.bonds = self.system.bonds[keep]
        self.bond_orders = self.bond_orders[keep]
        self.system.neighs()
        self._ensure_numeric_arrays()

    def move_atom(self, index, xy):
        self.system.apos[index, 0] = xy[0]
        self.system.apos[index, 1] = xy[1]

    def change_element(self, index, element):
        self._ensure_mutable_lists()
        element = canonical_element(element)
        self.system.enames[index] = element
        atomic_number = ELEMENT_DICT[element][INDEX_Z]
        self.system.atypes[index] = atomic_number
        self.system.qs[index] = ELEMENT_DICT[element][INDEX_VALENCE]
        self.system.Rs[index] = ELEMENT_DICT[element][INDEX_Rvdw]

    def load(self, path):
        ext = Path(path).suffix.lower()
        if ext not in {".xyz", ".mol", ".mol2"}:
            raise ValueError(f"Unsupported file format: {ext}")
        system = AtomicSystem(fname=str(path))
        if system.enames is None:
            system.enames = []
        system.enames = [canonical_element(str(e)) for e in system.enames]
        if system.atypes is None:
            system.atypes = np.array([ELEMENT_DICT[e][INDEX_Z] for e in system.enames], dtype=np.int32)
        if system.qs is None:
            system.qs = np.zeros(len(system.enames), dtype=np.float64)
        if system.Rs is None:
            system.Rs = np.zeros(len(system.enames), dtype=np.float64)
        system.findBonds()
        bond_orders = np.ones(system.bonds.shape[0], dtype=np.int8)
        system.neighs()
        self.system = system
        self._ensure_mutable_lists()
        self._ensure_numeric_arrays()
        self._ensure_aux_labels()
        self.bond_orders = bond_orders
        self.selected_atoms.clear()
        self.selected_bonds.clear()

    def save(self, path):
        ext = Path(path).suffix.lower()
        if ext == ".xyz":
            self.system.saveXYZ(path)
            return
        if ext == ".mol":
            self.system.save_mol(path)
            return
        if ext == ".mol2":
            self.system.save_mol2(path)
            return
        raise ValueError(f"Unsupported file format: {ext}")


class MoleculeCanvas(FigureCanvas):
    def __init__(self, document, parent=None):
        self.document = document
        self.figure = Figure(figsize=(5, 5))
        super().__init__(self.figure)
        self.setParent(parent)
        self.ax = self.figure.add_subplot(111)
        self.ax.set_aspect("equal")
        self.ax.set_axis_off()
        self.atom_artist = None
        self.selection_artist = None
        self.bond_collection = None
        self.mode = MODE_DRAW
        self.current_element = DEFAULT_ELEMENTS[1]
        self.current_bond_order = 1
        self.pending_atom = None
        self.dragging_atom = None
        self.last_mouse = None
        self._connect_events()
        self.refresh()

    def _connect_events(self):
        self.mpl_connect("button_press_event", self.on_mouse_press)
        self.mpl_connect("button_release_event", self.on_mouse_release)
        self.mpl_connect("motion_notify_event", self.on_mouse_move)

    def set_mode(self, mode):
        self.mode = mode
        self.pending_atom = None
        self.dragging_atom = None

    def set_element(self, element):
        self.current_element = element

    def set_bond_order(self, order):
        self.current_bond_order = order

    def refresh(self):
        self.ax.clear()
        self.ax.set_aspect("equal")
        self.ax.set_axis_off()
        doc = self.document
        apos = doc.positions
        if apos.shape[0] > 0:
            xy = apos[:, :2]
            enames = list(doc.system.enames)
            n_atoms = xy.shape[0]
            if len(enames) != n_atoms:
                enames = [canonical_element(enames[i]) if i < len(enames) else "C" for i in range(n_atoms)]
                doc.system.enames = enames
            colors = [elements.getColor(e) for e in enames]
            radii = np.array([ELEMENT_DICT[e][INDEX_Rvdw] for e in enames], dtype=np.float64) * 600.0
            self.atom_artist = self.ax.scatter(xy[:, 0], xy[:, 1], c=colors, s=radii, zorder=3)
            if doc.selected_atoms:
                sel_idx = [i for i in doc.selected_atoms if 0 <= i < n_atoms]
                doc.selected_atoms = set(sel_idx)
                if sel_idx:
                    sel_xy = xy[sel_idx]
                    self.selection_artist = self.ax.scatter(sel_xy[:, 0], sel_xy[:, 1], facecolors="none", edgecolors="yellow", s=np.ones(len(sel_idx)) * 700, linewidths=2, zorder=4)
        bonds = doc.system.bonds
        if bonds.shape[0] > 0 and apos.shape[0] > 0:
            segments = []
            widths = []
            xy = doc.positions[:, :2]
            for ib, (i, j) in enumerate(bonds):
                segments.append([xy[i], xy[j]])
                widths.append(1.0 + 2.0 * doc.bond_orders[ib])
            self.bond_collection = LineCollection(segments, colors="black", linewidths=widths, zorder=2)
            self.ax.add_collection(self.bond_collection)
            if doc.selected_bonds:
                doc.selected_bonds = {ib for ib in doc.selected_bonds if 0 <= ib < len(segments)}
        self.figure.canvas.draw_idle()

    def hit_atom(self, pos, pixel_tol=10):
        if pos is None:
            return None
        apos = self.document.positions
        if apos.shape[0] == 0:
            return None
        xy = apos[:, :2]
        dx = xy[:, 0] - pos[0]
        dy = xy[:, 1] - pos[1]
        d2 = dx * dx + dy * dy
        if d2.size == 0:
            return None
        tol = (pixel_tol / 100.0) ** 2
        idx = int(np.argmin(d2))
        if d2[idx] < tol:
            return idx
        return None

    @staticmethod
    def point_segment_distance(point, a, b):
        ap = point - a
        ab = b - a
        ab2 = np.dot(ab, ab)
        if ab2 == 0:
            return sqrt(np.dot(ap, ap))
        t = np.clip(np.dot(ap, ab) / ab2, 0.0, 1.0)
        closest = a + ab * t
        return sqrt(np.sum((point - closest) ** 2))

    def hit_bond(self, pos, tol=0.05):
        bonds = self.document.system.bonds
        xy = self.document.positions[:, :2]
        for ib, (i, j) in enumerate(bonds):
            a = xy[i]
            b = xy[j]
            if self.point_segment_distance(np.array(pos), a, b) < tol:
                return ib
        return None

    def select_atom(self, index, additive=False):
        if not additive:
            self.document.selected_atoms.clear()
        if index is not None:
            self.document.selected_atoms.add(index)
        self.document.selected_bonds.clear()
        self.refresh()

    def select_bond(self, index):
        self.document.selected_atoms.clear()
        if index is not None:
            self.document.selected_bonds = {index}
        else:
            self.document.selected_bonds.clear()
        self.refresh()

    def on_mouse_press(self, event):
        if event.inaxes != self.ax:
            return
        pos = (event.xdata, event.ydata)
        doc = self.document
        if event.button == 1:
            if self.mode == MODE_DRAW:
                idx = self.hit_atom(pos)
                if idx is None:
                    new_idx = doc.add_atom(pos, self.current_element)
                    doc.selected_atoms = {new_idx}
                    doc.selected_bonds.clear()
                    self.pending_atom = new_idx
                    self.refresh()
                else:
                    if self.pending_atom is None:
                        self.pending_atom = idx
                        self.select_atom(idx)
                    else:
                        other = idx
                        doc.set_bond(self.pending_atom, other, self.current_bond_order)
                        doc.selected_bonds = {doc.find_bond_index(self.pending_atom, other)}
                        doc.selected_atoms.clear()
                        self.pending_atom = None
                        self.refresh()
            else:  # move mode
                idx = self.hit_atom(pos)
                if idx is not None:
                    self.dragging_atom = idx
                    self.select_atom(idx)
                    self.last_mouse = pos
        elif event.button == 3:
            idx = self.hit_atom(pos)
            if idx is not None:
                doc.delete_atom(idx)
                self.pending_atom = None
                self.dragging_atom = None
                self.refresh()
                return
            bond_idx = self.hit_bond(pos)
            if bond_idx is not None:
                i, j = doc.system.bonds[bond_idx]
                doc.delete_bond(int(i), int(j))
                self.refresh()

    def on_mouse_release(self, event):
        if event.button != 1:
            return
        self.dragging_atom = None
        self.last_mouse = None

    def on_mouse_move(self, event):
        if event.inaxes != self.ax:
            return
        if self.dragging_atom is None:
            return
        pos = (event.xdata, event.ydata)
        if pos[0] is None or pos[1] is None:
            return
        self.document.move_atom(self.dragging_atom, pos)
        self.refresh()


class ToolPanel(QWidget):
    def __init__(self, canvas):
        super().__init__(canvas)
        self.canvas = canvas
        layout = QVBoxLayout()
        layout.setContentsMargins(4, 4, 4, 4)
        layout.setSpacing(6)

        layout.addWidget(QLabel("Mode"))
        mode_row = QHBoxLayout()
        draw_button = QToolButton()
        draw_button.setText("Draw")
        draw_button.setCheckable(True)
        move_button = QToolButton()
        move_button.setText("Move")
        move_button.setCheckable(True)
        self.mode_group = QButtonGroup(self)
        self.mode_group.setExclusive(True)
        self.mode_group.addButton(draw_button)
        self.mode_group.addButton(move_button)
        draw_button.setChecked(True)
        self.mode_group.buttonClicked.connect(self._mode_changed)
        mode_row.addWidget(draw_button)
        mode_row.addWidget(move_button)
        layout.addLayout(mode_row)

        layout.addWidget(QLabel("Elements"))
        self.element_list = QListWidget()
        for sym in DEFAULT_ELEMENTS:
            item = QListWidgetItem(sym)
            item.setSizeHint(QSize(40, 24))
            self.element_list.addItem(item)
        self.element_list.setCurrentRow(1)
        self.element_list.currentTextChanged.connect(canvas.set_element)
        layout.addWidget(self.element_list)

        layout.addWidget(QLabel("Bond order"))
        self.bond_combo = QComboBox()
        self.bond_combo.addItems(["1", "2", "3"])
        self.bond_combo.currentTextChanged.connect(self._bond_changed)
        layout.addWidget(self.bond_combo)

        auto_button = QPushButton("Auto bonds")
        auto_button.clicked.connect(self._auto_bonds)
        layout.addWidget(auto_button)

        layout.addStretch(1)
        self.setLayout(layout)

    def _mode_changed(self, button):
        mode = MODE_DRAW if button.text().lower() == "draw" else MODE_MOVE
        self.canvas.set_mode(mode)

    def _bond_changed(self, text):
        try:
            order = int(text)
        except ValueError:
            order = 1
        self.canvas.set_bond_order(order)

    def _auto_bonds(self):
        self.canvas.document.auto_bonds()
        self.canvas.refresh()


class MoleculeEditorWindow(QMainWindow):
    def __init__(self, document):
        super().__init__()
        self.setWindowTitle("Molecule Editor 2D")
        self.document = document
        self.canvas = MoleculeCanvas(document, self)
        self.tool_panel = ToolPanel(self.canvas)

        central = QWidget()
        layout = QHBoxLayout(central)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        layout.addWidget(self.canvas, 1)

        side_frame = QFrame()
        side_layout = QVBoxLayout(side_frame)
        side_layout.setContentsMargins(4, 4, 4, 4)
        side_layout.addWidget(self.tool_panel)
        side_layout.addStretch(1)
        layout.addWidget(side_frame)

        self.setCentralWidget(central)
        self._build_menu()

    def _build_menu(self):
        menubar = self.menuBar()
        file_menu = menubar.addMenu("File")

        open_act = file_menu.addAction("Open...")
        open_act.triggered.connect(self.open_file)

        save_act = file_menu.addAction("Save As...")
        save_act.triggered.connect(self.save_file)

        file_menu.addSeparator()
        exit_act = file_menu.addAction("Exit")
        exit_act.triggered.connect(self.close)

    def open_file(self):
        path, _ = QFileDialog.getOpenFileName(self, "Open molecule", "", "Molecules (*.xyz *.mol *.mol2)")
        if not path:
            return
        self.document.load(path)
        self.canvas.refresh()
        self.setWindowTitle(f"Molecule Editor 2D - {Path(path).name}")

    def save_file(self):
        path, _ = QFileDialog.getSaveFileName(self, "Save molecule", "", "XYZ (*.xyz);;MOL (*.mol);;MOL2 (*.mol2)")
        if not path:
            return
        self.document.save(path)

def launch_editor(path=None):
    doc = MoleculeDocument()
    if path:
        doc.load(path)
    app = QApplication.instance() or QApplication(sys.argv)
    win = MoleculeEditorWindow(doc)
    win.resize(1000, 700)
    win.show()
    return app.exec_()

if __name__ == "__main__":
    argv = sys.argv[1:]
    path = None
    if len(argv) == 1:
        path = argv[0]
    sys.exit(launch_editor(path))