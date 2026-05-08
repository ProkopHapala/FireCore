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
    __slots__ = ('pos', 'ename', 'atype', 'pin', 'parent', 'subtype', 'bonds', 'neighbors', '_id', 'alive')
    _counter = 0

    def __init__(self, pos, ename, atype, pin=None, parent=None, subtype=''):
        Atom._counter += 1
        self._id = Atom._counter
        self.alive   = True         # False = marked for deletion, will be cleaned up
        self.pos     = np.asarray(pos, dtype=np.float64)   # (3,)
        self.ename   = ename        # element symbol string
        self.atype   = atype        # integer Z
        self.pin     = pin          # (rx, ry) grid node key or None
        self.parent  = parent       # Atom object (heavy atom this H belongs to) or None
        self.subtype = subtype      # 'C_sp2', 'N_sp3', 'H_cap', etc.
        self.bonds   = []           # list of Bond objects involving this atom
        self.neighbors = []         # list of neighboring Atoms (derived from bonds)

    def __repr__(self):
        status = "" if self.alive else "[DEAD]"
        return f"Atom({self._id}{status} {self.ename} pin={self.pin} pos={self.pos[:2]})"


# ─── Bond ───────────────────────────────────────────────────────────────────

class Bond:
    __slots__ = ('a', 'b', 'order', '_id', 'alive')
    _counter = 0

    def __init__(self, a: Atom, b: Atom, order=1):
        Bond._counter += 1
        self._id   = Bond._counter
        self.alive = True           # False = marked for deletion, will be cleaned up
        self.a     = a
        self.b     = b
        self.order = order

    def other(self, atom: Atom) -> Atom:
        return self.b if atom is self.a else self.a

    def __repr__(self):
        status = "" if self.alive else "[DEAD]"
        return f"Bond({self._id}{status} {self.a._id}-{self.b._id} o={self.order})"


# ─── Ring ───────────────────────────────────────────────────────────────────

class Ring:
    __slots__ = ('atoms', 'bonds', 'cog', '_id', 'alive')
    _counter = 0

    def __init__(self, atoms, bonds):
        """Ring as real geometry cycle (n-gon).
        Args:
            atoms: list[Atom] - ordered list of atoms in the cycle
            bonds: list[Bond] - ordered list of bonds in the cycle
        """
        Ring._counter += 1
        self._id   = Ring._counter
        self.alive = True           # False = marked for deletion, will be cleaned up
        self.atoms = list(atoms)    # [Atom, ...] — ordered cyclically
        self.bonds = list(bonds)    # [Bond, ...] — ordered cyclically
        self.cog   = self._compute_cog()

    def _compute_cog(self):
        """Compute center of geometry as average of atom positions."""
        # Only count alive atoms
        alive_atoms = [a for a in self.atoms if a.alive]
        if not alive_atoms:
            return np.zeros(3)
        positions = np.array([a.pos for a in alive_atoms])
        return np.mean(positions, axis=0)

    def __repr__(self):
        status = "" if self.alive else "[DEAD]"
        return f"Ring({self._id}{status} natoms={len(self.atoms)})"


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

    def remove_atom(self, atom: Atom, soft=True):
        """Remove atom. If soft=True, mark as dead and cleanup later.
        If soft=False, immediate hard removal."""
        if atom._id not in self.atoms:
            return
        if soft:
            # Soft deletion: mark as dead, will be cleaned up later
            atom.alive = False
            # Also mark all its bonds as dead
            for b in atom.bonds:
                b.alive = False
        else:
            # Hard deletion: immediate removal
            for b in list(atom.bonds):
                self._remove_bond_internal(b, hard=True)
            if atom.pin is not None:
                self._pin_to_atom.pop(atom.pin, None)
            del self.atoms[atom._id]

    def cleanup_invalid(self):
        """Remove all dead (alive=False) atoms, bonds, and rings.
        Clean up references from other objects."""
        # First, remove dead bonds from atom bond lists
        for atom in self.atoms.values():
            if atom.alive:
                atom.bonds = [b for b in atom.bonds if b.alive]
        
        # Remove dead rings
        dead_ring_ids = [rid for rid, r in self.rings.items() if not r.alive]
        for rid in dead_ring_ids:
            del self.rings[rid]
        
        # Remove dead bonds from main bonds dict
        dead_bond_ids = [bid for bid, b in self.bonds.items() if not b.alive]
        for bid in dead_bond_ids:
            del self.bonds[bid]
        
        # Remove dead atoms and update pin mapping
        dead_atom_ids = [aid for aid, a in self.atoms.items() if not a.alive]
        for aid in dead_atom_ids:
            atom = self.atoms[aid]
            if atom.pin is not None:
                self._pin_to_atom.pop(atom.pin, None)
            del self.atoms[aid]
        
        return len(dead_atom_ids), len(dead_bond_ids), len(dead_ring_ids)

    def sync_neighbor_lists(self):
        """Rebuild neighbor lists from alive bonds.
        Call this after any bond topology change."""
        # Clear all neighbor lists
        for atom in self.atoms.values():
            if atom.alive:
                atom.neighbors = []
        # Rebuild from alive bonds
        for bond in self.bonds.values():
            if bond.alive:
                bond.a.neighbors.append(bond.b)
                bond.b.neighbors.append(bond.a)

    def h_children(self, heavy_atom: Atom) -> list:
        """Return list of H atoms (subtype='H_cap') that have parent=heavy_atom."""
        return [a for a in self.atoms.values() 
                if a.alive and a.subtype == 'H_cap' and a.parent is heavy_atom]

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

    def _remove_bond_internal(self, bond: Bond, hard=False):
        if bond._id not in self.bonds:
            return
        if hard:
            bond.a.bonds = [b for b in bond.a.bonds if b is not bond]
            bond.b.bonds = [b for b in bond.b.bonds if b is not bond]
            del self.bonds[bond._id]
        else:
            bond.alive = False

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
        # Build adjacency list (only for alive atoms with alive bonds to alive atoms)
        adj = {a._id: [b.other(a) for b in a.bonds if b.alive and b.other(a).alive] for a in self.atoms.values() if a.alive}
        visited = set()
        rings = []

        def dfs(start, current, path_atoms, path_bonds, visited_edges):
            if len(path_atoms) > max_ring_size:
                return
            if current._id in visited:
                return
            for neighbor in adj.get(current._id, []):
                if neighbor._id not in adj:
                    continue  # Skip dead/removed neighbors
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
            if atom._id in visited or not atom.alive:
                continue
            dfs(atom, atom, [], [], set())

        return rings

    # ── Bulk bond rebuild ─────────────────────────────────────────────────────

    def recalc_bonds(self, bond_length=1.42, tol_factor=0.35):
        """Remove all bonds and recompute from distance threshold."""
        # Use hard delete to immediately remove bonds from atom bond lists
        for bond in list(self.bonds.values()):
            self._remove_bond_internal(bond, hard=True)
        # Only consider alive atoms
        atoms = [a for a in self.atoms.values() if a.alive]
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
        Only alive objects are included.
        """
        # Only include alive atoms
        atom_list = [a for a in self.atoms.values() if a.alive]
        idx = {a._id: i for i, a in enumerate(atom_list)}
        enames = np.array([a.ename for a in atom_list], dtype=object)
        apos   = np.array([a.pos   for a in atom_list], dtype=np.float64)
        atypes = np.array([a.atype for a in atom_list], dtype=np.int32)
        
        # Only include alive bonds between alive atoms
        bond_pairs = []
        bond_list = []
        for bond in self.bonds.values():
            if not bond.alive:
                continue
            # Both atoms must be alive
            if not bond.a.alive or not bond.b.alive:
                continue
            ia = idx.get(bond.a._id)
            ib = idx.get(bond.b._id)
            if ia is not None and ib is not None:
                bond_pairs.append((ia, ib))
                bond_list.append(bond)
        bonds = np.array(bond_pairs, dtype=np.int32).reshape(-1, 2) if bond_pairs else np.zeros((0, 2), dtype=np.int32)
        
        # Only include alive rings
        ring_list = [r for r in self.rings.values() if r.alive]
        return atom_list, enames, apos, atypes, bonds, bond_list, ring_list

    # ── Position update ───────────────────────────────────────────────────────

    def update_positions_from_array(self, apos):
        """Update atom positions from array (same order as to_arrays()).
        
        Args:
            apos: (N,3) array of positions, where N matches len(atoms) and order matches to_arrays()
        
        This updates geometry only (atom positions), not topology (bonds, rings).
        Used after external geometry relaxation (e.g., DFTB) to sync relaxed positions back to graph.
        """
        atom_list = [a for a in self.atoms.values() if a.alive]
        if len(atom_list) != len(apos):
            raise ValueError(f"Position array length {len(apos)} does not match number of alive atoms {len(atom_list)}")
        for i, atom in enumerate(atom_list):
            atom.pos[:] = apos[i]

    # ── Convenience queries ───────────────────────────────────────────────────

    def heavy_atoms(self):
        return [a for a in self.atoms.values() if a.alive and a.ename not in ('H', 'E')]

    def h_children(self, atom: Atom):
        return [a for a in self.atoms.values() if a.alive and a.parent is atom]

    def neighbors(self, atom: Atom):
        return [b.other(atom) for b in atom.bonds if b.alive]

    # ── Picking helpers ────────────────────────────────────────────────────────

    def pick_atom(self, pos, radius=0.5):
        """Find atom within radius of position. Returns Atom or None."""
        for atom in self.atoms.values():
            if atom.alive and np.linalg.norm(atom.pos - pos) < radius:
                return atom
        return None

    def pick_bond(self, pos, radius=0.5):
        """Find bond whose center is within radius of position. Returns Bond or None."""
        for bond in self.bonds.values():
            if not bond.alive:
                continue
            # Both atoms must be alive
            if not bond.a.alive or not bond.b.alive:
                continue
            center = (bond.a.pos + bond.b.pos) / 2
            if np.linalg.norm(center - pos) < radius:
                return bond
        return None

    def pick_ring(self, pos, radius=1.0):
        """Find ring whose COG is within radius of position. Returns Ring or None."""
        for ring in self.rings.values():
            if ring.alive and np.linalg.norm(ring.cog - pos) < radius:
                return ring
        return None

    def __repr__(self):
        return f"AtomicGraph(atoms={len(self.atoms)}, bonds={len(self.bonds)}, rings={len(self.rings)})"
