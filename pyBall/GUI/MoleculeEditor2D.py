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
from PyQt5.QtCore import QSize, Qt
from PyQt5.QtWidgets import (
    QApplication,
    QMainWindow,
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QFormLayout,
    QGroupBox,
    QLabel,
    QPushButton,
    QComboBox,
    QFileDialog,
    QListWidget,
    QListWidgetItem,
    QDoubleSpinBox,
    QToolButton,
    QButtonGroup,
    QFrame,
)
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from matplotlib.collections import LineCollection
from matplotlib.patches import Rectangle

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
MODE_BOND = "bond"


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

    def auto_bonds(self, Rcut=None, RvdwCut=None):
        if Rcut is None:
            Rcut = self.autobond_rcut
        if RvdwCut is None:
            RvdwCut = self.autobond_rvdw
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

    def move_atoms(self, indices, delta):
        if not indices:
            return
        arr = self.system.apos
        if arr.shape[0] == 0:
            return
        idx = np.array(list(indices), dtype=int)
        if idx.size == 0:
            return
        arr[idx, 0] += delta[0]
        arr[idx, 1] += delta[1]

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


class MoleculeCanvas(FigureCanvasQTAgg):
    def __init__(self, document, parent=None):
        self.document = document
        self.figure = Figure(figsize=(5, 5))
        super().__init__(self.figure)
        self.setParent(parent)
        self.ax = self.figure.add_subplot(111)
        self.figure.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=1.0)
        self.ax.set_position([0.0, 0.0, 1.0, 1.0])
        self.ax.set_aspect("equal")
        self.ax.set_axis_off()
        for spine in self.ax.spines.values():
            spine.set_visible(True)
            spine.set_color("#999999")
            spine.set_linewidth(1.0)
        self.atom_artist = None
        self.selection_artist = None
        self.bond_collection = None
        self.mode = MODE_DRAW
        self.current_element = DEFAULT_ELEMENTS[1]
        self.current_bond_order = 1
        self.pending_atom = None
        self.dragging_atom = None
        self.last_mouse = None
        self.drag_select_start = None
        self.selection_rect = None
        self.selection_mode = "replace"
        self.arrow_step = 0.1
        self.pick_radius_scale = 0.6
        self.autobond_rcut = 3.0
        self.autobond_rvdw = 1.2
        self.view_center = np.zeros(2, dtype=np.float64)
        self.view_zoom = 1.0
        self.tool_panel = None
        self._connect_events()
        self.autoscale_view(refresh=False)
        self.refresh()

    def _connect_events(self):
        self.mpl_connect("button_press_event", self.on_mouse_press)
        self.mpl_connect("button_release_event", self.on_mouse_release)
        self.mpl_connect("motion_notify_event", self.on_mouse_move)

    def set_mode(self, mode):
        self.mode = mode
        self.pending_atom = None
        self.dragging_atom = None
        self.last_mouse = None

    def set_element(self, element):
        self.current_element = element

    def set_bond_order(self, order):
        self.current_bond_order = order

    def refresh(self):
        self.ax.clear()
        self.ax.set_position([0.0, 0.0, 1.0, 1.0])
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
                    sel_sizes = np.ones(len(sel_idx)) * 900.0
                    self.selection_artist = self.ax.scatter(
                        sel_xy[:, 0],
                        sel_xy[:, 1],
                        s=sel_sizes,
                        c="yellow",
                        edgecolors="black",
                        linewidths=1.5,
                        zorder=4,
                    )
                    self.selection_artist.set_alpha(0.6)
                else:
                    self.selection_artist = None
            else:
                self.selection_artist = None
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
        scale = max(0.1, float(self.view_zoom))
        cx, cy = self.view_center
        half_size = 5.0 / scale
        self.ax.set_xlim(cx - half_size, cx + half_size)
        self.ax.set_ylim(cy - half_size, cy + half_size)
        if self.selection_rect is not None:
            if self.selection_rect not in self.ax.patches:
                self.ax.add_patch(self.selection_rect)
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
        idx = int(np.argmin(d2))
        doc = self.document
        radii = getattr(doc.system, "Rs", None)
        if radii is None or len(radii) != xy.shape[0]:
            radii = np.array([
                ELEMENT_DICT[canonical_element(e)][INDEX_Rvdw]
                for e in doc.system.enames[: xy.shape[0]]
            ], dtype=np.float64)
        else:
            radii = np.asarray(radii, dtype=np.float64)
        if radii.size == 0:
            limit = self.pick_radius_scale
        else:
            limit = max(0.01, self.pick_radius_scale * radii[idx])
        if d2[idx] <= limit * limit:
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

    def _cycle_bond_order(self, bond_idx):
        if bond_idx is None or bond_idx < 0:
            return
        doc = self.document
        if bond_idx >= len(doc.bond_orders):
            return
        current = int(doc.bond_orders[bond_idx]) if len(doc.bond_orders) else 1
        new_order = 1 if current >= 3 else current + 1
        doc.bond_orders[bond_idx] = new_order
        self.refresh()

    def select_atom(self, index, additive=False, subtract=False):
        if not additive and not subtract:
            self.document.selected_atoms.clear()
        if index is not None:
            if subtract:
                self.document.selected_atoms.discard(index)
            else:
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

    def _in_axes(self, pos):
        if pos[0] is None or pos[1] is None:
            return False
        xlim = self.ax.get_xlim()
        ylim = self.ax.get_ylim()
        return (xlim[0] <= pos[0] <= xlim[1]) and (ylim[0] <= pos[1] <= ylim[1])

    def on_mouse_press(self, event):
        pos = (event.xdata, event.ydata)
        inside = (event.inaxes == self.ax) and self._in_axes(pos)
        add_mode, subtract_mode = self._modifier_state(event)
        doc = self.document
        if event.button == 1:
            if not inside:
                return
            if self.mode == MODE_MOVE:
                idx = self.hit_atom(pos)
                if idx is not None:
                    self.dragging_atom = idx
                    self.last_mouse = pos
                return
            bond_idx = self.hit_bond(pos)
            if bond_idx is not None:
                self._cycle_bond_order(bond_idx)
                return
            if self.mode == MODE_DRAW:
                idx = self.hit_atom(pos)
                if idx is None:
                    new_idx = doc.add_atom(pos, self.current_element)
                    doc.selected_atoms = {new_idx}
                    doc.selected_bonds.clear()
                    self.refresh()
                else:
                    # allow tapping atom to prepare for bond placement
                    self.pending_atom = idx
                    self.select_atom(idx, additive=False, subtract=False)
            elif self.mode == MODE_BOND:
                idx = self.hit_atom(pos)
                if idx is None:
                    return
                if self.pending_atom is None:
                    self.pending_atom = idx
                    self.select_atom(idx)
                else:
                    other = idx
                    if other != self.pending_atom:
                        ib = doc.set_bond(self.pending_atom, other, self.current_bond_order)
                        if ib >= 0:
                            doc.selected_bonds = {ib}
                        doc.selected_atoms.clear()
                    self.pending_atom = None
                    self.refresh()
        elif event.button == 3:
            if not inside:
                return
            if self.mode == MODE_MOVE:
                idx = self.hit_atom(pos)
                if idx is not None:
                    self.select_atom(idx, additive=add_mode, subtract=subtract_mode)
                    return
                bond_idx = self.hit_bond(pos)
                if bond_idx is not None:
                    self.select_bond(bond_idx)
                    return
                self._begin_rect_selection(pos, add_mode, subtract_mode)
                return
            bond_idx = self.hit_bond(pos)
            if bond_idx is not None:
                i, j = doc.system.bonds[bond_idx]
                doc.delete_bond(int(i), int(j))
                doc.selected_bonds.clear()
                doc.selected_atoms.clear()
                self.refresh()
                return
            idx = self.hit_atom(pos)
            if idx is not None:
                doc.delete_atom(idx)
                doc.selected_atoms.clear()
                doc.selected_bonds.clear()
                self.pending_atom = None
                self.dragging_atom = None
                self.refresh()
                return
            self._begin_rect_selection(pos, add_mode, subtract_mode)

    def on_mouse_release(self, event):
        if event.button == 3:
            if self.drag_select_start is not None:
                pos = (event.xdata, event.ydata)
                self._finalize_rect_selection(pos)
            return
        if event.button != 1:
            return
        if self.drag_select_start is not None:
            pos = (event.xdata, event.ydata)
            self._finalize_rect_selection(pos)
            return
        self.dragging_atom = None
        self.last_mouse = None

    def on_mouse_move(self, event):
        if event.inaxes != self.ax:
            return
        pos = (event.xdata, event.ydata)
        if self.drag_select_start is not None:
            if pos[0] is None or pos[1] is None:
                return
            self._update_rect_selection(pos)
            return
        if self.dragging_atom is None:
            return
        anchor = self.last_mouse
        if anchor is None or pos[0] is None or pos[1] is None:
            return
        dx = pos[0] - anchor[0]
        dy = pos[1] - anchor[1]
        if dx == 0 and dy == 0:
            return
        target_indices = list(self.document.selected_atoms)
        if target_indices and self.dragging_atom in target_indices:
            self.document.move_atoms(target_indices, (dx, dy))
        else:
            self.document.move_atom(self.dragging_atom, pos)
        self.last_mouse = pos
        self.refresh()

    def keyPressEvent(self, event):
        key = event.key()
        doc = self.document
        if key in (Qt.Key_Delete, Qt.Key_Backspace):
            removed = sorted(doc.selected_atoms)
            for idx in reversed(removed):
                doc.delete_atom(idx)
            if doc.selected_bonds:
                for ib in sorted(doc.selected_bonds, reverse=True):
                    i, j = doc.system.bonds[ib]
                    doc.delete_bond(int(i), int(j))
            doc.selected_atoms.clear()
            doc.selected_bonds.clear()
            self.refresh()
            return
        step = float(self.arrow_step)
        delta = None
        if key == Qt.Key_Left:
            delta = (-step, 0.0)
        elif key == Qt.Key_Right:
            delta = (step, 0.0)
        elif key == Qt.Key_Up:
            delta = (0.0, step)
        elif key == Qt.Key_Down:
            delta = (0.0, -step)
        if delta is not None:
            if doc.selected_atoms:
                doc.move_atoms(doc.selected_atoms, delta)
                self.refresh()
            return
        if key in (Qt.Key_Plus, Qt.Key_Equal, Qt.Key_Plus + Qt.KeypadModifier):
            self.set_view_zoom(self.view_zoom * 1.1)
            return
        if key in (Qt.Key_Minus, Qt.Key_Minus + Qt.KeypadModifier):
            self.set_view_zoom(self.view_zoom / 1.1)
            return

    def set_pick_radius_scale(self, scale):
        value = float(scale)
        if value <= 0.0:
            raise ValueError("pick radius scale must be positive")
        self.pick_radius_scale = value

    def set_autobond_params(self, Rcut=None, Rvdw=None):
        if Rcut is not None:
            value = float(Rcut)
            if value <= 0.0:
                raise ValueError("Rcut must be positive")
            self.autobond_rcut = value
        if Rvdw is not None:
            value = float(Rvdw)
            if value <= 0.0:
                raise ValueError("Rvdw multiplier must be positive")
            self.autobond_rvdw = value
        if self.tool_panel is not None:
            self.tool_panel.sync_autobond_controls()

    def set_view_center(self, x=None, y=None, update_controls=True):
        changed = False
        if x is not None:
            self.view_center[0] = float(x)
            changed = True
        if y is not None:
            self.view_center[1] = float(y)
            changed = True
        if changed:
            if update_controls and self.tool_panel is not None:
                self.tool_panel.sync_view_controls()
            self.refresh()

    def set_view_zoom(self, zoom, update_controls=True):
        value = float(zoom)
        if value <= 0.0:
            raise ValueError("zoom must be positive")
        self.view_zoom = np.clip(value, 0.1, 10.0)
        if update_controls and self.tool_panel is not None:
            self.tool_panel.sync_view_controls()
        self.refresh()

    def autoscale_view(self, refresh=True, margin=0.2):
        doc = self.document
        apos = doc.positions
        if apos.shape[0] == 0:
            self.view_center[:] = 0.0
            self.view_zoom = 1.0
        else:
            xy = apos[:, :2]
            min_xy = xy.min(axis=0)
            max_xy = xy.max(axis=0)
            span = max_xy - min_xy
            extent = float(np.max(span))
            if extent < 1e-3:
                extent = 1.0
            width = extent * (1.0 + margin)
            zoom = 10.0 / width
            self.view_center[:] = 0.5 * (min_xy + max_xy)
            self.view_zoom = np.clip(zoom, 0.1, 10.0)
        if self.tool_panel is not None:
            self.tool_panel.sync_view_controls()
        if refresh:
            self.refresh()

    def _modifier_state(self, event):
        qevent = getattr(event, "guiEvent", None)
        if qevent is None:
            return False, False
        modifiers = qevent.modifiers()
        add_mode = bool(modifiers & Qt.ShiftModifier)
        subtract_mode = bool(modifiers & Qt.ControlModifier)
        return add_mode, subtract_mode

    def _begin_rect_selection(self, pos, add_mode, subtract_mode):
        if pos[0] is None or pos[1] is None:
            return
        self.drag_select_start = pos
        if self.selection_rect is None:
            self.selection_rect = Rectangle(
                (pos[0], pos[1]),
                0.0,
                0.0,
                facecolor="yellow",
                edgecolor="orange",
                alpha=0.2,
                linewidth=1.0,
                zorder=1,
            )
        else:
            self.selection_rect.set_xy((pos[0], pos[1]))
            self.selection_rect.set_width(0.0)
            self.selection_rect.set_height(0.0)
        self.selection_mode = "replace"
        if add_mode:
            self.selection_mode = "add"
        elif subtract_mode:
            self.selection_mode = "subtract"
        if self.selection_rect not in self.ax.patches:
            self.ax.add_patch(self.selection_rect)
        self.figure.canvas.draw_idle()

    def _update_rect_selection(self, pos):
        if self.selection_rect is None or self.drag_select_start is None:
            return
        x0, y0 = self.drag_select_start
        if pos[0] is None or pos[1] is None:
            return
        min_x = min(x0, pos[0])
        min_y = min(y0, pos[1])
        width = abs(pos[0] - x0)
        height = abs(pos[1] - y0)
        self.selection_rect.set_xy((min_x, min_y))
        self.selection_rect.set_width(width)
        self.selection_rect.set_height(height)
        self.figure.canvas.draw_idle()

    def _finalize_rect_selection(self, pos):
        if self.drag_select_start is None:
            return
        if pos[0] is None or pos[1] is None:
            pos = self.drag_select_start
        x0, y0 = self.drag_select_start
        min_x = min(x0, pos[0])
        max_x = max(x0, pos[0])
        min_y = min(y0, pos[1])
        max_y = max(y0, pos[1])
        indices = []
        if self.document.atom_count() > 0:
            xy = self.document.positions[:, :2]
            mask = (
                (xy[:, 0] >= min_x)
                & (xy[:, 0] <= max_x)
                & (xy[:, 1] >= min_y)
                & (xy[:, 1] <= max_y)
            )
            indices = np.where(mask)[0].tolist()
        if self.selection_mode == "replace":
            self.document.selected_atoms = set(indices)
        elif self.selection_mode == "add":
            self.document.selected_atoms.update(indices)
        elif self.selection_mode == "subtract":
            for idx in indices:
                self.document.selected_atoms.discard(idx)
        self.document.selected_bonds.clear()
        self.drag_select_start = None
        if self.selection_rect is not None:
            self.selection_rect.set_width(0.0)
            self.selection_rect.set_height(0.0)
            self.selection_rect.set_xy((0.0, 0.0))
            if self.selection_rect in self.ax.patches:
                self.selection_rect.remove()
            self.selection_rect = None
        self.figure.canvas.draw_idle()
        self.refresh()


class ToolPanel(QWidget):
    def __init__(self, canvas):
        super().__init__(canvas)
        self.canvas = canvas
        main_layout = QVBoxLayout()
        main_layout.setContentsMargins(4, 4, 4, 4)
        main_layout.setSpacing(8)

        mode_box = QGroupBox("Mode")
        mode_layout = QHBoxLayout()
        draw_button = QToolButton()
        draw_button.setText("Draw")
        draw_button.setCheckable(True)
        draw_button.setProperty("mode", MODE_DRAW)
        bond_button = QToolButton()
        bond_button.setText("Bond")
        bond_button.setCheckable(True)
        bond_button.setProperty("mode", MODE_BOND)
        move_button = QToolButton()
        move_button.setText("Move")
        move_button.setCheckable(True)
        move_button.setProperty("mode", MODE_MOVE)
        self.mode_group = QButtonGroup(self)
        self.mode_group.setExclusive(True)
        self.mode_group.addButton(draw_button)
        self.mode_group.addButton(bond_button)
        self.mode_group.addButton(move_button)
        draw_button.setChecked(True)
        self.mode_group.buttonClicked.connect(self._mode_changed)
        mode_layout.addWidget(draw_button)
        mode_layout.addWidget(bond_button)
        mode_layout.addWidget(move_button)
        mode_box.setLayout(mode_layout)
        main_layout.addWidget(mode_box)

        elements_box = QGroupBox("Elements")
        elements_layout = QVBoxLayout()
        self.element_list = QListWidget()
        for sym in DEFAULT_ELEMENTS:
            item = QListWidgetItem(sym)
            item.setSizeHint(QSize(40, 24))
            self.element_list.addItem(item)
        self.element_list.setCurrentRow(1)
        self.element_list.itemSelectionChanged.connect(self._element_changed)
        elements_layout.addWidget(self.element_list)
        elements_box.setLayout(elements_layout)
        main_layout.addWidget(elements_box)

        bond_box = QGroupBox("Bond settings")
        bond_form = QFormLayout()
        self.bond_combo = QComboBox()
        self.bond_combo.addItems(["1", "2", "3"])
        self.bond_combo.currentTextChanged.connect(self._bond_changed)
        bond_form.addRow("Bond order", self.bond_combo)

        self.autobond_rcut_spin = QDoubleSpinBox()
        self.autobond_rcut_spin.setDecimals(2)
        self.autobond_rcut_spin.setRange(0.5, 10.0)
        self.autobond_rcut_spin.setSingleStep(0.1)
        self.autobond_rcut_spin.setValue(self.canvas.autobond_rcut)
        self.autobond_rcut_spin.valueChanged.connect(self._autobond_rcut_changed)
        bond_form.addRow("Rcut (Å)", self.autobond_rcut_spin)

        self.autobond_rvdw_spin = QDoubleSpinBox()
        self.autobond_rvdw_spin.setDecimals(2)
        self.autobond_rvdw_spin.setRange(0.1, 2.0)
        self.autobond_rvdw_spin.setSingleStep(0.05)
        self.autobond_rvdw_spin.setValue(self.canvas.autobond_rvdw)
        self.autobond_rvdw_spin.valueChanged.connect(self._autobond_rvdw_changed)
        bond_form.addRow("Rvdw ×", self.autobond_rvdw_spin)

        auto_button = QPushButton("Auto bonds")
        auto_button.clicked.connect(self._auto_bonds)
        bond_form.addRow(auto_button)
        bond_box.setLayout(bond_form)
        main_layout.addWidget(bond_box)

        selection_box = QGroupBox("Selection & move")
        selection_form = QFormLayout()
        self.move_step = QDoubleSpinBox()
        self.move_step.setDecimals(2)
        self.move_step.setRange(0.01, 10.0)
        self.move_step.setSingleStep(0.05)
        self.move_step.setValue(self.canvas.arrow_step)
        self.move_step.valueChanged.connect(self._move_step_changed)
        selection_form.addRow("Arrow step", self.move_step)

        self.pick_radius_spin = QDoubleSpinBox()
        self.pick_radius_spin.setDecimals(2)
        self.pick_radius_spin.setRange(0.05, 5.0)
        self.pick_radius_spin.setSingleStep(0.05)
        self.pick_radius_spin.setValue(self.canvas.pick_radius_scale)
        self.pick_radius_spin.valueChanged.connect(self._pick_radius_changed)
        selection_form.addRow("Pick radius", self.pick_radius_spin)
        selection_box.setLayout(selection_form)
        main_layout.addWidget(selection_box)

        view_box = QGroupBox("View")
        view_form = QFormLayout()
        self.view_zoom_spin = QDoubleSpinBox()
        self.view_zoom_spin.setDecimals(2)
        self.view_zoom_spin.setRange(0.1, 10.0)
        self.view_zoom_spin.setSingleStep(0.1)
        self.view_zoom_spin.setValue(self.canvas.view_zoom)
        self.view_zoom_spin.valueChanged.connect(self._view_zoom_changed)
        view_form.addRow("Zoom", self.view_zoom_spin)

        self.view_center_x_spin = QDoubleSpinBox()
        self.view_center_x_spin.setDecimals(2)
        self.view_center_x_spin.setRange(-1000.0, 1000.0)
        self.view_center_x_spin.setSingleStep(0.1)
        self.view_center_x_spin.setValue(self.canvas.view_center[0])
        self.view_center_x_spin.valueChanged.connect(self._view_center_x_changed)
        view_form.addRow("Center X", self.view_center_x_spin)

        self.view_center_y_spin = QDoubleSpinBox()
        self.view_center_y_spin.setDecimals(2)
        self.view_center_y_spin.setRange(-1000.0, 1000.0)
        self.view_center_y_spin.setSingleStep(0.1)
        self.view_center_y_spin.setValue(self.canvas.view_center[1])
        self.view_center_y_spin.valueChanged.connect(self._view_center_y_changed)
        view_form.addRow("Center Y", self.view_center_y_spin)

        autoscale_button = QPushButton("Autoscale view")
        autoscale_button.clicked.connect(self._autoscale_view)
        view_form.addRow(autoscale_button)
        view_box.setLayout(view_form)
        main_layout.addWidget(view_box)

        main_layout.addStretch(1)
        self.setLayout(main_layout)
        self.canvas.tool_panel = self
        self.sync_autobond_controls()
        self.sync_view_controls()

    def _mode_changed(self, button):
        mode = button.property("mode") or MODE_DRAW
        self.canvas.set_mode(mode)

    def _element_changed(self):
        items = self.element_list.selectedItems()
        if not items:
            return
        self.canvas.set_element(items[0].text())

    def _bond_changed(self, text):
        try:
            order = int(text)
        except ValueError:
            order = 1
        self.canvas.set_bond_order(order)

    def _move_step_changed(self, value):
        self.canvas.arrow_step = float(value)

    def _pick_radius_changed(self, value):
        self.canvas.set_pick_radius_scale(value)

    def _autobond_rcut_changed(self, value):
        self.canvas.set_autobond_params(Rcut=value)

    def _autobond_rvdw_changed(self, value):
        self.canvas.set_autobond_params(Rvdw=value)

    def _view_zoom_changed(self, value):
        self.canvas.set_view_zoom(value, update_controls=False)

    def _view_center_x_changed(self, value):
        self.canvas.set_view_center(x=value, update_controls=False)

    def _view_center_y_changed(self, value):
        self.canvas.set_view_center(y=value, update_controls=False)

    def _autoscale_view(self):
        self.canvas.autoscale_view()

    def _auto_bonds(self):
        self.canvas.document.auto_bonds(
            Rcut=self.canvas.autobond_rcut,
            RvdwCut=self.canvas.autobond_rvdw,
        )
        self.canvas.refresh()

    def sync_autobond_controls(self):
        self.autobond_rcut_spin.blockSignals(True)
        self.autobond_rcut_spin.setValue(self.canvas.autobond_rcut)
        self.autobond_rcut_spin.blockSignals(False)
        self.autobond_rvdw_spin.blockSignals(True)
        self.autobond_rvdw_spin.setValue(self.canvas.autobond_rvdw)
        self.autobond_rvdw_spin.blockSignals(False)

    def sync_view_controls(self):
        self.view_zoom_spin.blockSignals(True)
        self.view_zoom_spin.setValue(self.canvas.view_zoom)
        self.view_zoom_spin.blockSignals(False)
        self.view_center_x_spin.blockSignals(True)
        self.view_center_x_spin.setValue(self.canvas.view_center[0])
        self.view_center_x_spin.blockSignals(False)
        self.view_center_y_spin.blockSignals(True)
        self.view_center_y_spin.setValue(self.canvas.view_center[1])
        self.view_center_y_spin.blockSignals(False)


class MoleculeEditor(QMainWindow):
    def __init__(self, document=None, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Molecule Editor 2D")
        self.resize(1000, 800)
        self.document = document
        self.canvas = MoleculeCanvas(self.document)
        self.canvas.setFocusPolicy(Qt.StrongFocus)
        self.canvas.setFocus()
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
        central.setLayout(layout)
        self.setCentralWidget(central)
        self.setFocusPolicy(Qt.StrongFocus)
        self.canvas.setFocusProxy(self)
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

    def keyPressEvent(self, event):
        if self.canvas is not None:
            self.canvas.keyPressEvent(event)

def launch_editor(path=None):
    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)
    doc = MoleculeDocument()
    if path is not None:
        doc.load(path)
    win = MoleculeEditor(doc)
    win.resize(1000, 700)
    win.show()
    return app.exec_()

if __name__ == "__main__":
    argv = sys.argv[1:]
    path = None
    if len(argv) == 1:
        path = argv[0]
    sys.exit(launch_editor(path))