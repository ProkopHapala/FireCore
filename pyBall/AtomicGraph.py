"""
AtomicGraph.py — Object-graph representation of an atomic structure.

Design principles:
  - Atom, Bond, Ring are plain Python objects with stable identity (not integer indices).
  - Deletion of any object does NOT renumber or invalidate any other object.
  - Integer indices for interop with numpy/vispy are generated on demand via to_arrays().
  - No parallel arrays that must be kept in sync; all per-atom data lives on the Atom object.

Public API:
  graph = AtomicGraph()
  a = graph.add_atom(pos, ename, pin=None, parent=None, subtype='C_sp2')
  graph.remove_atom(a)          # removes a and all its bonds; caller handles rings
  b = graph.add_bond(a1, a2)
  graph.remove_bond(b)
  r = graph.add_ring(q, r_coord, atoms)
  graph.remove_ring(r)
  atoms, enames, apos, bonds = graph.to_arrays()   # for numpy/vispy rendering
"""

import numpy as np

# ─── Atom ───────────────────────────────────────────────────────────────────

class Atom:
    __slots__ = ('pos', 'ename', 'atype', 'pin', 'parent', 'subtype', 'bonds', '_id')
    _counter = 0

    def __init__(self, pos, ename, atype, pin=None, parent=None, subtype=''):
        Atom._counter += 1
        self._id = Atom._counter
        self.pos     = np.asarray(pos, dtype=np.float64)   # (3,)
        self.ename   = ename        # element symbol string
        self.atype   = atype        # integer Z
        self.pin     = pin          # (rx, ry) grid node key or None
        self.parent  = parent       # Atom object (heavy atom this H belongs to) or None
        self.subtype = subtype      # 'C_sp2', 'N_sp3', 'H_cap', etc.
        self.bonds   = []           # list of Bond objects involving this atom

    def __repr__(self):
        return f"Atom({self._id} {self.ename} pin={self.pin} pos={self.pos[:2]})"


# ─── Bond ───────────────────────────────────────────────────────────────────

class Bond:
    __slots__ = ('a', 'b', 'order', '_id')
    _counter = 0

    def __init__(self, a: Atom, b: Atom, order=1):
        Bond._counter += 1
        self._id  = Bond._counter
        self.a    = a
        self.b    = b
        self.order = order

    def other(self, atom: Atom) -> Atom:
        return self.b if atom is self.a else self.a

    def __repr__(self):
        return f"Bond({self.a._id}-{self.b._id} o={self.order})"


# ─── Ring ───────────────────────────────────────────────────────────────────

class Ring:
    __slots__ = ('q', 'r', 'atoms', '_id')
    _counter = 0

    def __init__(self, q: int, r: int, atoms):
        Ring._counter += 1
        self._id  = Ring._counter
        self.q    = q
        self.r    = r
        self.atoms = list(atoms)    # [Atom, ...] — object refs, not indices

    def __repr__(self):
        return f"Ring({self._id} ({self.q},{self.r}) natoms={len(self.atoms)})"


# ─── AtomicGraph ────────────────────────────────────────────────────────────

class AtomicGraph:
    """Mutable graph of atoms, bonds, and rings.
    All collections are dicts keyed by object id for O(1) lookup and deletion.
    """

    def __init__(self):
        self.atoms  = {}    # id -> Atom
        self.bonds  = {}    # id -> Bond
        self.rings  = {}    # (q,r) -> Ring   (axial key is the natural unique key)
        self._pin_to_atom = {}   # (rx,ry) -> Atom  — kept in sync with every add/remove

    # ── Atom operations ──────────────────────────────────────────────────────

    def add_atom(self, pos, ename, atype, pin=None, parent=None, subtype='') -> Atom:
        a = Atom(pos, ename, atype, pin=pin, parent=parent, subtype=subtype)
        self.atoms[a._id] = a
        if pin is not None:
            assert pin not in self._pin_to_atom, f"Duplicate pin {pin} (existing={self._pin_to_atom[pin]}, new={a})"
            self._pin_to_atom[pin] = a
        return a

    def remove_atom(self, atom: Atom):
        """Remove atom, its bonds, and its pin mapping. Does NOT touch rings."""
        if atom._id not in self.atoms:
            return
        for b in list(atom.bonds):
            self._remove_bond_internal(b)
        if atom.pin is not None:
            self._pin_to_atom.pop(atom.pin, None)
        del self.atoms[atom._id]

    def atom_at_pin(self, pin) -> 'Atom | None':
        return self._pin_to_atom.get(pin)

    # ── Bond operations ──────────────────────────────────────────────────────

    def add_bond(self, a: Atom, b: Atom, order=1) -> Bond:
        for bond in a.bonds:
            if bond.other(a) is b:
                return bond   # already exists
        bond = Bond(a, b, order)
        self.bonds[bond._id] = bond
        a.bonds.append(bond)
        b.bonds.append(bond)
        return bond

    def remove_bond(self, bond: Bond):
        self._remove_bond_internal(bond)

    def _remove_bond_internal(self, bond: Bond):
        if bond._id not in self.bonds:
            return
        bond.a.bonds = [b for b in bond.a.bonds if b is not bond]
        bond.b.bonds = [b for b in bond.b.bonds if b is not bond]
        del self.bonds[bond._id]

    def get_bond(self, a: Atom, b: Atom) -> 'Bond | None':
        for bond in a.bonds:
            if bond.other(a) is b:
                return bond
        return None

    # ── Ring operations ──────────────────────────────────────────────────────

    def add_ring(self, q: int, r: int, atoms) -> Ring:
        key = (q, r)
        assert key not in self.rings, f"Ring {key} already exists"
        ring = Ring(q, r, atoms)
        self.rings[key] = ring
        return ring

    def remove_ring(self, q: int, r: int):
        self.rings.pop((q, r), None)

    def ring_at(self, q: int, r: int) -> 'Ring | None':
        return self.rings.get((q, r))

    # ── Bulk bond rebuild ─────────────────────────────────────────────────────

    def recalc_bonds(self, bond_length=1.42, tol_factor=0.35):
        """Remove all bonds and recompute from distance threshold."""
        for bond in list(self.bonds.values()):
            self._remove_bond_internal(bond)
        atoms = list(self.atoms.values())
        threshold = bond_length * (1.0 + tol_factor)
        threshold_sq = threshold ** 2
        for i, a in enumerate(atoms):
            for j in range(i + 1, len(atoms)):
                b = atoms[j]
                d2 = float(np.sum((a.pos - b.pos) ** 2))
                if d2 < threshold_sq:
                    self.add_bond(a, b)

    # ── Export for numpy/vispy ─────────────────────────────────────────────────

    def to_arrays(self):
        """Return (atom_list, enames, apos, bonds_idx) for rendering.
        atom_list[i] is the Atom object at index i.
        bonds_idx is (N,2) int array of indices into atom_list.
        Index assignment is stable within one call; call again after mutations.
        """
        atom_list = list(self.atoms.values())
        idx = {a._id: i for i, a in enumerate(atom_list)}
        enames = np.array([a.ename for a in atom_list], dtype=object)
        apos   = np.array([a.pos   for a in atom_list], dtype=np.float64)
        atypes = np.array([a.atype for a in atom_list], dtype=np.int32)
        bond_pairs = []
        for bond in self.bonds.values():
            ia = idx.get(bond.a._id)
            ib = idx.get(bond.b._id)
            if ia is not None and ib is not None:
                bond_pairs.append((ia, ib))
        bonds = np.array(bond_pairs, dtype=np.int32).reshape(-1, 2) if bond_pairs else np.zeros((0, 2), dtype=np.int32)
        return atom_list, enames, apos, atypes, bonds

    # ── Convenience queries ───────────────────────────────────────────────────

    def heavy_atoms(self):
        return [a for a in self.atoms.values() if a.ename not in ('H', 'E')]

    def h_children(self, atom: Atom):
        return [a for a in self.atoms.values() if a.parent is atom]

    def neighbors(self, atom: Atom):
        return [b.other(atom) for b in atom.bonds]

    def __repr__(self):
        return f"AtomicGraph(atoms={len(self.atoms)}, bonds={len(self.bonds)}, rings={len(self.rings)})"
