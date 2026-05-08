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
    __slots__ = ('atoms', 'bonds', 'cog', '_id')
    _counter = 0

    def __init__(self, atoms, bonds):
        """Ring as real geometry cycle (n-gon).
        Args:
            atoms: list[Atom] - ordered list of atoms in the cycle
            bonds: list[Bond] - ordered list of bonds in the cycle
        """
        Ring._counter += 1
        self._id   = Ring._counter
        self.atoms = list(atoms)    # [Atom, ...] — ordered cyclically
        self.bonds = list(bonds)    # [Bond, ...] — ordered cyclically
        self.cog   = self._compute_cog()

    def _compute_cog(self):
        """Compute center of geometry as average of atom positions."""
        if not self.atoms:
            return np.zeros(3)
        positions = np.array([a.pos for a in self.atoms])
        return np.mean(positions, axis=0)

    def __repr__(self):
        return f"Ring({self._id} natoms={len(self.atoms)})"


# ─── AtomicGraph ────────────────────────────────────────────────────────────

class AtomicGraph:
    """Mutable graph of atoms, bonds, and rings.
    All collections are dicts keyed by object id for O(1) lookup and deletion.
    """

    def __init__(self):
        self.atoms  = {}    # id -> Atom
        self.bonds  = {}    # id -> Bond
        self.rings  = {}    # id -> Ring
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

    def add_ring(self, atoms, bonds) -> Ring:
        """Add a geometry-based ring (n-gon cycle).
        Args:
            atoms: list[Atom] - ordered list of atoms in the cycle
            bonds: list[Bond] - ordered list of bonds in the cycle
        """
        ring = Ring(atoms, bonds)
        self.rings[ring._id] = ring
        return ring

    def remove_ring(self, ring: Ring):
        self.rings.pop(ring._id, None)

    def detect_rings(self, max_ring_size=8):
        """Detect all rings (cycles) in the bond graph using DFS.
        Returns list of Ring objects.
        """
        # Build adjacency list
        adj = {a._id: [b.other(a) for b in a.bonds] for a in self.atoms.values()}
        visited = set()
        rings = []

        def dfs(start, current, path_atoms, path_bonds, visited_edges):
            if len(path_atoms) > max_ring_size:
                return
            if current._id in visited:
                return
            for neighbor in adj[current._id]:
                edge = self.get_bond(current, neighbor)
                if edge is None:
                    continue
                edge_key = frozenset((current._id, neighbor._id))
                if edge_key in visited_edges:
                    continue
                if neighbor._id == start._id and len(path_atoms) >= 3:
                    # Found a cycle
                    rings.append(self.add_ring(path_atoms + [neighbor], path_bonds + [edge]))
                    continue
                if neighbor._id not in [a._id for a in path_atoms]:
                    dfs(start, neighbor, path_atoms + [neighbor], path_bonds + [edge],
                        visited_edges | {edge_key})

        for atom in self.atoms.values():
            if atom._id in visited:
                continue
            dfs(atom, atom, [], [], set())

        return rings

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
        """Return (atom_list, enames, apos, atypes, bonds_idx, bond_list, ring_list) for rendering.
        atom_list[i] is the Atom object at index i.
        bonds_idx is (N,2) int array of indices into atom_list.
        bond_list[i] is the Bond object at index i (parallel to bonds_idx).
        ring_list is list of Ring objects.
        Index assignment is stable within one call; call again after mutations.
        """
        atom_list = list(self.atoms.values())
        idx = {a._id: i for i, a in enumerate(atom_list)}
        enames = np.array([a.ename for a in atom_list], dtype=object)
        apos   = np.array([a.pos   for a in atom_list], dtype=np.float64)
        atypes = np.array([a.atype for a in atom_list], dtype=np.int32)
        bond_pairs = []
        bond_list = []
        for bond in self.bonds.values():
            ia = idx.get(bond.a._id)
            ib = idx.get(bond.b._id)
            if ia is not None and ib is not None:
                bond_pairs.append((ia, ib))
                bond_list.append(bond)
        bonds = np.array(bond_pairs, dtype=np.int32).reshape(-1, 2) if bond_pairs else np.zeros((0, 2), dtype=np.int32)
        ring_list = list(self.rings.values())
        return atom_list, enames, apos, atypes, bonds, bond_list, ring_list

    # ── Convenience queries ───────────────────────────────────────────────────

    def heavy_atoms(self):
        return [a for a in self.atoms.values() if a.ename not in ('H', 'E')]

    def h_children(self, atom: Atom):
        return [a for a in self.atoms.values() if a.parent is atom]

    def neighbors(self, atom: Atom):
        return [b.other(atom) for b in atom.bonds]

    # ── Picking helpers ────────────────────────────────────────────────────────

    def pick_atom(self, pos, radius=0.5):
        """Find atom within radius of position. Returns Atom or None."""
        for atom in self.atoms.values():
            if np.linalg.norm(atom.pos - pos) < radius:
                return atom
        return None

    def pick_bond(self, pos, radius=0.5):
        """Find bond whose center is within radius of position. Returns Bond or None."""
        for bond in self.bonds.values():
            center = (bond.a.pos + bond.b.pos) / 2
            if np.linalg.norm(center - pos) < radius:
                return bond
        return None

    def pick_ring(self, pos, radius=1.0):
        """Find ring whose COG is within radius of position. Returns Ring or None."""
        for ring in self.rings.values():
            if np.linalg.norm(ring.cog - pos) < radius:
                return ring
        return None

    def __repr__(self):
        return f"AtomicGraph(atoms={len(self.atoms)}, bonds={len(self.bonds)}, rings={len(self.rings)})"
