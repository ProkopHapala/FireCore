"""
KekuleBackend.py - Backend logic for the Kekule Structure Explorer
===================================================================
Manages persistent AtomicSystem state with hexagonal grid metadata.
AtomicSystem is the authoritative store; grid mappings are auxiliary.
Each UI action triggers exactly one mutation of the backend state.
"""

import sys
import os
import numpy as np
from pyBall.AtomicSystem import AtomicSystem
from pyBall.AtomicGraph import AtomicGraph
from pyBall import elements
from pyBall import atomicUtils as au

# Global verbosity level for debug prints (sync with KekuleExplorerGUI)
# 0: Only exceptions and explicit prints
# 1: Warnings and complex operation reports
# 2: Click and action prints (default)
# 3: Hovered prints (most verbose)
VERBOSITY_LEVEL = 2

def debug_print(level, message):
    """Print message if verbosity level is >= specified level."""
    global VERBOSITY_LEVEL
    if VERBOSITY_LEVEL >= level:
        print(message)

# ============ Honeycomb geometry helpers ============

def honeycomb_ring_nodes(q, r, a_CC=1.42):
    """Return the 6 node positions (in Cartesian) of a hexagonal ring at axial coords (q,r).
    
    Uses pointy-top hexagon orientation. The ring center is at:
        cx = a_CC * sqrt(3) * (q + r/2)
        cy = a_CC * 1.5 * r
    
    Args:
        q, r: axial hex coordinates of the ring center
        a_CC: C-C bond length (Angstrom)
    Returns:
        (6,2) array of node positions in xy
    """
    s3 = np.sqrt(3.0)
    cx = a_CC * s3 * (q + r * 0.5)
    cy = a_CC * 1.5 * r
    # 6 vertices of a pointy-top hexagon with circumradius = a_CC
    angles = np.arange(6) * (np.pi / 3.0) + np.pi / 6.0  # start at 30 degrees
    nodes = np.column_stack([cx + a_CC * np.cos(angles), cy + a_CC * np.sin(angles)])
    return nodes

def snap_to_grid(pos_xy, a_CC=1.42, tol=0.15):
    """Snap a Cartesian position to the nearest honeycomb node. Returns rounded tuple as key."""
    decimals = 4
    return (round(float(pos_xy[0]), decimals), round(float(pos_xy[1]), decimals))

# ============ KekuleBackend class ============

# Passivation group definitions: list of (element, x, y, z) coordinates relative to C atom
# First atom at (0,0,0) replaces the C atom; others are added at C_pos + coords
# For top/bottom edges, y coordinates are scaled by direction (+1 or -1)
PASSIVATION_GROUPS = {
    'N': [('N', 0.0, 0.0, 0.0)],
    'NH': [('N', 0.0, 0.0, 0.0), ('H', 0.0, 1.01, 0.0)],
    'CH': [('H', 0.0, 1.09, 0.0)],
    'H': [('H', 0.0, 1.09, 0.0)],
    'O': [('O', 0.0, 0.0, 0.0)],
    'C=O': [('O', 0.0, 1.23, 0.0)],
    'C-OH': [('O', 0.0, 1.43, 0.0), ('H', 0.31, 2.34, 0.0)],  # H at 109.5° from y-axis
}

# Passivation string encoding mapping for CLI
# Each character in the string represents one passivation group at one site along the edge
PASSIVATION_ENCODING = {
    'n': 'NH',
    'N': 'N',
    'o': 'C=O',
    'O': 'O',
    'H': 'CH',
    'h': 'C-OH'
}

def parse_passivation_string(s):
    """Convert passivation string to list of passivation group names.
    
    Encoding:
    - n -> NH
    - N -> N
    - o -> C=O
    - O -> O
    - H -> CH
    - h -> C-OH
    
    Each character in the string represents one passivation group at one site.
    """
    if s is None:
        return None
    result = []
    for char in s:
        if char not in PASSIVATION_ENCODING:
            raise ValueError(f"Unknown passivation character: '{char}'. Valid: {list(PASSIVATION_ENCODING.keys())}")
        result.append(PASSIVATION_ENCODING[char])
    return result

bond_length_cutoff = 2.0  # Slightly larger than typical C-C bond (1.42)
class KekuleBackend:
    """Manages molecular editing state on a hexagonal grid.

    Authoritative state: self.graph (AtomicGraph) — object graph, no integer indices.
    Render state:        self.sys  (AtomicSystem) — synced on demand via _sync_sys().

    AtomicGraph contains Atom objects.  Each Atom carries:
        .pin     : (rx,ry) grid node key, or None for H caps
        .parent  : Atom object (heavy atom this H belongs to), or None
        .subtype : 'C_sp2', 'H_cap', etc.
    Rings are stored in graph.rings as {(q,r): Ring(atoms=[Atom,...])}
    """

    def __init__(self, a_CC=1.42):
        self.a_CC = a_CC
        self.graph = AtomicGraph()
        self.sys   = AtomicSystem(apos=np.zeros((0,3)), atypes=np.array([],dtype=np.int32), enames=np.array([],dtype=object))
        self.pbc_x = False
        self.pbc_y = False
        self.vacuum_gap = 15.0
        self.hex_mode = 'Hex1'  # 'Hex1' = paint mode (force add/remove), 'Hex2' = toggle mode (preserve shared)
        # Separate storage for hex grid tiles (for Hex1/Hex2 editing)
        self.hex_tiles = set()  # set of (q,r) tuples
        self._rings_dirty = True  # flag to re-detect geometry rings
        self.auto_h_cap = True  # auto-adjust hydrogen caps when structure changes
        self.auto_recalc_bonds = False  # auto-recalculate bonds based on distance (DANGEROUS - can create spurious bonds)

    # ── Properties for GUI compatibility ─────────────────────────────────────

    @property
    def rings(self):
        """Return set of hex tile coordinates for Hex1/Hex2 mode (legacy)."""
        return self.hex_tiles

    @property
    def geometry_rings(self):
        """Return list of geometry-based Ring objects from bond graph."""
        return list(self.graph.rings.values())

    @property
    def atom_pin(self):
        """List of pin values in same order as to_arrays() atom_list."""
        atom_list, *_ = self.graph.to_arrays()
        return [a.pin for a in atom_list]

    @property
    def atom_parent(self):
        atom_list, *_ = self.graph.to_arrays()
        idx = {a._id: i for i, a in enumerate(atom_list)}
        return [idx.get(a.parent._id) if a.parent is not None else None for a in atom_list]

    @property
    def atom_subtype(self):
        atom_list, *_ = self.graph.to_arrays()
        return [a.subtype for a in atom_list]

    @atom_subtype.setter
    def atom_subtype(self, value):
        """Allow list assignment (used by adjust_h legacy path)."""
        atom_list, *_ = self.graph.to_arrays()
        for i, a in enumerate(atom_list):
            if i < len(value):
                a.subtype = value[i]

    def set_atom_subtype_by_index(self, atom_idx, subtype):
        """Set subtype of atom at specific index directly."""
        atom_list, *_ = self.graph.to_arrays()
        if 0 <= atom_idx < len(atom_list):
            atom_list[atom_idx].subtype = subtype

    @property
    def ring_atoms(self):
        """Dict {ring_id: [int indices]} for geometry rings."""
        atom_list, *_ = self.graph.to_arrays()
        idx = {a._id: i for i, a in enumerate(atom_list)}
        result = {}
        for ring_id, ring in self.graph.rings.items():
            result[ring_id] = [idx[a._id] for a in ring.atoms if a._id in idx]
        return result

    # ── Internal helpers ──────────────────────────────────────────────────────

    def _build_node_to_atom(self):
        """Return {pin_key: int_index} built fresh from graph."""
        atom_list, *_ = self.graph.to_arrays()
        return {a.pin: i for i, a in enumerate(atom_list) if a.pin is not None}

    def _build_node_to_atom_obj(self):
        """Return {pin_key: Atom} — the fast O(1) version backed by graph._pin_to_atom."""
        return self.graph._pin_to_atom  # already maintained by AtomicGraph

    def _atom_to_index(self, atom):
        """Convert Atom object to int index in sys arrays."""
        atom_list, *_ = self.graph.to_arrays()
        idx_map = {a._id: i for i, a in enumerate(atom_list)}
        return idx_map.get(atom._id)

    def _rebuild_after_delete(self, indices_to_remove):
        """Remove atoms by index and rebuild sys arrays from graph."""
        # Remove atoms from graph
        atom_list, *_ = self.graph.to_arrays()
        for idx in sorted(indices_to_remove, reverse=True):
            if idx < len(atom_list):
                a = atom_list[idx]
                # Remove H children first
                for h in self.graph.h_children(a):
                    self.graph.remove_atom(h)
                self.graph.remove_atom(a)
        # Rebuild sys from graph
        self._sync_sys()

    def _append_atom(self, pos, ename, pin=None, parent=None, subtype=''):
        """Add atom to graph. parent may be int index (legacy) or Atom object."""
        atype = elements.ELEMENT_DICT[ename][0]
        if isinstance(parent, int):
            atom_list, *_ = self.graph.to_arrays()
            parent = atom_list[parent] if 0 <= parent < len(atom_list) else None
        a = self.graph.add_atom(np.array([pos[0], pos[1], pos[2] if len(pos) > 2 else 0.0]),
                                ename, atype, pin=pin, parent=parent, subtype=subtype)
        return a  # return Atom object; callers that need int index use graph.to_arrays()

    def _sync_sys(self):
        """Rebuild self.sys arrays from graph for rendering/export. Called before any sys read."""
        atom_list, enames, apos, atypes, bonds, bond_list, ring_list = self.graph.to_arrays()
        self.sys.apos    = apos.astype(np.float64)
        self.sys.enames  = enames
        self.sys.atypes  = atypes
        self.sys.natoms  = len(atom_list)
        self.sys.bonds   = bonds if len(bonds) else None
        self.sys.ngs     = None  # invalidate neighbor cache
        # Store bond_list and ring_list for picking/visualization
        self._bond_list = bond_list
        self._ring_list = ring_list

    def detect_geometry_rings(self, max_ring_size=8):
        """Detect geometry rings from bond graph and store them in graph.rings.
        Only re-detects if _rings_dirty flag is set (structure changed).
        """
        if not self._rings_dirty:
            return list(self.graph.rings.values())
        # Clear existing geometry rings (but keep hex tiles)
        self.graph.rings.clear()
        # Detect rings from bond graph
        rings = self.graph.detect_rings(max_ring_size=max_ring_size)
        self._rings_dirty = False
        return rings

    # ── Picking helpers ────────────────────────────────────────────────────────

    def pick_atom(self, pos, radius=0.5):
        """Find atom within radius of position. Returns Atom or None."""
        return self.graph.pick_atom(pos, radius)

    def pick_bond(self, pos, radius=0.5):
        """Find bond whose center is within radius of position. Returns Bond or None."""
        return self.graph.pick_bond(pos, radius)

    def pick_ring(self, pos, radius=1.0):
        """Find ring whose COG is within radius of position. Returns Ring or None."""
        return self.graph.pick_ring(pos, radius)
    
    def _get_npi_from_subtype(self, subtype):
        """Extract npi (number of pi bonds) from subtype string."""
        if 'sp3' in subtype: return 0
        elif 'sp2' in subtype: return 1
        elif 'sp' in subtype: return 2
        else: return 1  # Default to sp2

    def _get_element_default_subtype(self, element):
        """Return default subtype for a newly added element."""
        if element == 'C':
            return 'C_sp2'
        elif element == 'N':
            return 'N_sp2'
        elif element == 'O':
            return 'O_sp2'
        elif element == 'H':
            return 'H_cap'
        return f'{element}_sp2'
    
    def _target_sigma(self, element, npi):
        """Target sigma bonds: nsigma = nval - npi - nepair.

        All 2nd-period atoms (C,N,O) want octet = 4 electron pairs.
        nval = 4 for C,N,O; nval = 1 for H (duet = 1 pair).
        nepair = 0 (C,H), 1 (N), 2 (O).
        """
        nval = 1 if element == 'H' else 4
        nepair = {'C': 0, 'N': 1, 'O': 2, 'H': 0}.get(element, 0)
        return nval - npi - nepair
    
    # --- Public mutation methods (each = exactly one mutation) ---

    def add_ring(self, q, r):
        """Add a benzene ring at axial position (q, r).

        Behavior depends on self.hex_mode:
        - Hex1 (paint): Add atoms at all 6 nodes if not present.
        - Hex2 (toggle): Add atoms only at empty nodes. Idempotent if already present.
        """
        key = (q, r)
        n2a = self.graph._pin_to_atom
        new_atoms = []
        for node in honeycomb_ring_nodes(q, r, self.a_CC):
            nk = snap_to_grid(node, self.a_CC)
            if nk not in n2a:
                a = self.graph.add_atom(np.array([nk[0], nk[1], 0.0]),
                                        'C', elements.ELEMENT_DICT['C'][0],
                                        pin=nk, parent=None, subtype='C_sp2')
                new_atoms.append(a)
        self.hex_tiles.add(key)
        if new_atoms:
            # Create bonds between all ring atoms at proper C-C distance
            # This bonds new atoms to each other AND to existing ring atoms
            all_ring_atoms = [n2a[snap_to_grid(node, self.a_CC)] 
                             for node in honeycomb_ring_nodes(q, r, self.a_CC)
                             if snap_to_grid(node, self.a_CC) in n2a]
            self._create_bonds_for_ring_atoms(all_ring_atoms)
            self.graph.sync_neighbor_lists()  # Sync neighbors after bond creation
        self._rings_dirty = True
        if self.auto_h_cap:
            self.adjust_h()
    
    def remove_ring(self, q, r):
        """Remove atoms at the 6 node positions of hexagon at axial (q,r).

        Behavior depends on self.hex_mode:
        - Hex1 (paint): Remove all atoms at 6 nodes (no sharing check).
        - Hex2 (toggle): Remove only atoms NOT shared with other hex tiles.

        In both modes, H atoms attached to removed heavy atoms are also removed.
        """
        key = (q, r)
        n2a = self.graph._pin_to_atom
        to_remove = []   # Atom objects
        for node in honeycomb_ring_nodes(q, r, self.a_CC):
            nk = snap_to_grid(node, self.a_CC)
            atom = n2a.get(nk)
            if atom is not None and atom not in to_remove:
                if self.hex_mode == 'Hex2':
                    # Toggle mode: preserve atoms shared with other hex tiles
                    shared = False
                    for other_q, other_r in self.hex_tiles:
                        if (other_q, other_r) == key:
                            continue
                        for other_node in honeycomb_ring_nodes(other_q, other_r, self.a_CC):
                            other_nk = snap_to_grid(other_node, self.a_CC)
                            if other_nk == nk:
                                shared = True
                                break
                        if shared:
                            break
                    if not shared:
                        to_remove.append(atom)
                else:
                    # Paint mode: remove all atoms
                    to_remove.append(atom)
        self.hex_tiles.discard(key)
        # Remove atoms and their H children
        for atom in to_remove:
            for h in self.graph.h_children(atom):
                self.graph.remove_atom(h)
            self.graph.remove_atom(atom)
        # Clean up dead bonds and sync (no recalc_bonds!)
        self.graph.cleanup_invalid()
        self.graph.sync_neighbor_lists()
        self._rings_dirty = True
        if self.auto_h_cap:
            self.adjust_h()

    def toggle_ring(self, q, r):
        """Toggle a benzene ring at axial position (q, r)."""
        if (q, r) in self.rings:
            self.remove_ring(q, r)
        else:
            self.add_ring(q, r)
    
    def snap_to_ring(self, x, y):
        """Find the axial coordinates (q, r) of the hexagon whose center is closest to (x, y)."""
        s3 = np.sqrt(3.0)
        r_exact = y / (1.5 * self.a_CC)
        q_exact = x / (s3 * self.a_CC) - r_exact * 0.5
        return int(round(q_exact)), int(round(r_exact))

    def snap_to_node(self, x, y, tol=0.2):
        """Find the rounded Cartesian coordinates (rx, ry) of the honeycomb node closest to (x, y)."""
        # Find the nearest hexagon first
        q, r = self.snap_to_ring(x, y)
        
        best_node = None
        min_dist = float('inf')
        # Check nodes of this hexagon and its 6 neighbors
        for dq, dr in [(0,0), (1,0), (-1,0), (0,1), (0,-1), (1,-1), (-1,1)]:
            nodes = honeycomb_ring_nodes(q+dq, r+dr, self.a_CC)
            for node in nodes:
                d = np.sqrt((node[0] - x)**2 + (node[1] - y)**2)
                if d < min_dist:
                    min_dist = d
                    best_node = node
        
        if min_dist < tol:
            return snap_to_grid(best_node, self.a_CC)
        return None

    def get_guide_points(self, qrange=(-10, 10), rrange=(-10, 10)):
        """Return all node positions in the specified range for UI guide dots."""
        nodes = set()
        for q in range(qrange[0], qrange[1] + 1):
            for r in range(rrange[0], rrange[1] + 1):
                ring_nodes = honeycomb_ring_nodes(q, r, self.a_CC)
                for node in ring_nodes:
                    nodes.add(snap_to_grid(node, self.a_CC))
        return np.array(list(nodes))

    def set_atom_type(self, node_key, element):
        """Set or change the element at a pinned grid node. Adds new atom if node is empty."""
        a = self.graph._pin_to_atom.get(node_key)
        if a is not None:
            a.ename   = element
            a.atype   = elements.ELEMENT_DICT[element][0]
            a.subtype = self._get_element_default_subtype(element)
        else:
            self.graph.add_atom(np.array([node_key[0], node_key[1], 0.0]),
                                element, elements.ELEMENT_DICT[element][0],
                                pin=node_key, parent=None,
                                subtype=self._get_element_default_subtype(element))
            # Create bond to nearest heavy atom (not H cap) within cutoff
            a = self.graph._pin_to_atom[node_key]
            self._create_bond_to_nearest_heavy(a)
            # Sync neighbor lists after adding bonds
            self.graph.sync_neighbor_lists()
        
        if self.auto_h_cap:
            self.adjust_h()

    def set_atom_type_by_index(self, atom_idx, element):
        """Set or change the element at a specific atom index (independent of grid)."""
        atom_list, *_ = self.graph.to_arrays()
        if 0 <= atom_idx < len(atom_list):
            a = atom_list[atom_idx]
            a.ename = element
            a.atype = elements.ELEMENT_DICT[element][0]
            a.subtype = self._get_element_default_subtype(element)
            if self.auto_h_cap:
                self.adjust_h()
            self._sync_sys()

    def _create_bond_to_nearest_heavy(self, atom):
        """Create bond to nearest heavy atom within bond length cutoff.
        
        This is used when adding atoms to ensure bonds are created immediately
        without the brute-force recalc_bonds that removes all bonds.
        """
        nearest_atom = None
        min_dist = float('inf')
        
        # Find nearest heavy atom (not H cap)
        for a in self.graph.atoms.values():
            if not a.alive or a == atom or a.subtype == 'H_cap':
                continue
            dist = np.linalg.norm(atom.pos - a.pos)
            if dist < min_dist:
                min_dist = dist
                nearest_atom = a
        
        # Create bond if within cutoff
        if nearest_atom and min_dist < bond_length_cutoff:
            # Check if bond already exists
            bond_exists = False
            for bond in atom.bonds:
                if bond.alive and (bond.a == nearest_atom or bond.b == nearest_atom):
                    bond_exists = True
                    break
            if not bond_exists:
                self.graph.add_bond(atom, nearest_atom)

    def _create_bonds_for_ring_atoms(self, ring_atoms):
        """Create bonds between all ring atoms that are at proper C-C bond distance.
        
        For a hexagon ring, each atom should be bonded to its two neighbors
        at ~1.42 Å distance (C-C bond length).
        """
        cc_bond_sq = (self.a_CC * 1.1) ** 2  # Slightly larger than a_CC for tolerance
        
        for i, a in enumerate(ring_atoms):
            if not a.alive:
                continue
            for j, b in enumerate(ring_atoms):
                if i >= j or not b.alive:
                    continue
                dist_sq = np.sum((a.pos - b.pos) ** 2)
                if dist_sq < cc_bond_sq:
                    # Check if bond already exists
                    bond_exists = False
                    for bond in a.bonds:
                        if bond.alive and (bond.a == b or bond.b == b):
                            bond_exists = True
                            break
                    if not bond_exists:
                        self.graph.add_bond(a, b)

    def set_atom_valency(self, node_key, npi):
        """Set npi (pi bond count) for atom at node_key. npi in {0,1,2}."""
        a = self.graph._pin_to_atom.get(node_key)
        if a is None: return
        sp_map = {0: 'sp3', 1: 'sp2', 2: 'sp'}
        a.subtype = f"{a.ename}_{sp_map.get(npi, 'sp2')}"
        if self.auto_h_cap:
            self.adjust_h()

    def add_atom_at_position(self, pos, element, npi=1):
        """Add an atom at arbitrary position (not on grid).

        Args:
            pos: (x, y, z) position
            element: element symbol (C, N, O, H)
            npi: number of pi bonds (default 1 for sp2)
        """
        ia = self._append_atom(
            pos=list(pos),
            ename=element,
            pin=None,  # Not pinned to grid
            parent=None,
            subtype=f"{element}_sp2" if npi == 1 else f"{element}_sp3"
        )
        self.recalc_bonds(skip_sync=True)
        if self.auto_h_cap:
            self.adjust_h()
        return ia

    def remove_atom(self, node_key):
        """Remove atom at snapped grid position node_key.
        
        Also removes any H atoms attached to the removed heavy atom.
        Does NOT call recalc_bonds() - just removes from graph and syncs.
        """
        n2a = self.graph._pin_to_atom
        a = n2a.get(node_key)
        if a is None:
            return
        ia = self._atom_to_index(a)
        self._rebuild_after_delete([ia])
        # Clean up dead bonds and sync neighbor lists (no recalc_bonds!)
        self.graph.cleanup_invalid()
        self.graph.sync_neighbor_lists()
        if self.auto_h_cap:
            self.adjust_h()

    def remove_atom_by_index(self, atom_idx):
        """Remove atom at specific index (independent of grid).
        
        Also removes any H atoms attached to the removed heavy atom.
        Does NOT call recalc_bonds() - just removes from graph and syncs.
        """
        self._rebuild_after_delete([atom_idx])
        # Clean up dead bonds and sync neighbor lists (no recalc_bonds!)
        self.graph.cleanup_invalid()
        self.graph.sync_neighbor_lists()
        if self.auto_h_cap:
            self.adjust_h()

    def insert_atom_into_bond(self, bond, new_ename='C'):
        """Insert a new atom into the middle of a bond, pushing original atoms aside.
        
        Original bond A-B becomes A-C-B where:
        - C is at the center of original A-B bond
        - A and B are pushed aside to satisfy A-C and B-C bond lengths
        
        Args:
            bond: Bond object (from graph.pick_bond())
            new_ename: Element name for new atom (default 'C')
        
        Returns:
            The newly created Atom object
        """
        if bond is None or not bond.alive:
            raise ValueError("Invalid bond (None or not alive)")
        
        # Get atom objects from bond (using object references, not indices!)
        atom_a = bond.a
        atom_b = bond.b
        
        debug_print(1, f"Inserting atom into bond {bond._id} between Atom({atom_a._id}) and Atom({atom_b._id})")
        
        pos_a = atom_a.pos
        pos_b = atom_b.pos
        
        # Center of original bond
        center = (pos_a + pos_b) / 2.0
        
        # Bond length for new element
        bond_lengths = {'C': 1.42, 'N': 1.33, 'O': 1.36, 'H': 1.09}
        bl = bond_lengths.get(new_ename, 1.42)
        
        # Direction from center to A
        dir_a = au.normalize(pos_a - center)
        # Direction from center to B (opposite)
        dir_b = -dir_a
        
        # Push A and B aside to satisfy bond lengths
        new_pos_a = center + dir_a * bl
        new_pos_b = center + dir_b * bl
        
        # Update positions of A and B (in the Atom objects!)
        atom_a.pos = new_pos_a
        atom_b.pos = new_pos_b
        
        debug_print(1, f"  Moved Atom({atom_a._id}) to {new_pos_a[:2]}")
        debug_print(1, f"  Moved Atom({atom_b._id}) to {new_pos_b[:2]}")
        
        # Mark original bond as dead (soft delete)
        bond.alive = False
        debug_print(1, f"  Marked bond {bond._id} as dead")
        
        # Add new atom C at center (using graph!)
        new_atom = self.graph.add_atom(center, new_ename, elements.ELEMENT_DICT[new_ename][0],
                                       pin=None, parent=None, subtype=f'{new_ename}_sp2')
        debug_print(1, f"  Created new Atom({new_atom._id}) at center {center[:2]}")
        
        # Add new bonds A-C and B-C (using graph!)
        bond_ac = self.graph.add_bond(atom_a, new_atom)
        bond_bc = self.graph.add_bond(atom_b, new_atom)
        debug_print(1, f"  Created new bonds: {bond_ac._id} (A-C), {bond_bc._id} (B-C)")
        
        # Cleanup dead objects
        n_atoms, n_bonds, n_rings = self.graph.cleanup_invalid()
        debug_print(1, f"  Cleanup removed: {n_atoms} atoms, {n_bonds} bonds, {n_rings} rings")
        
        # Sync neighbor lists from bonds (derived data)
        self.graph.sync_neighbor_lists()
        
        # Sync sys from graph
        self._sync_sys()
        self._rings_dirty = True
        
        # Auto-adjust H if enabled (after topology operation completes)
        if self.auto_h_cap:
            self.adjust_h()
        
        return new_atom

    def collapse_bond(self, bond, mouse_pos):
        """Collapse a bond by removing one atom and transferring its bonds to the other.
        
        Original bond A-B becomes A (or B) where:
        - The atom farther from mouse position survives
        - The surviving atom is moved to the center of original A-B bond
        - All neighbors of the removed atom are transferred to the survivor
        
        Args:
            bond: Bond object (from graph.pick_bond())
            mouse_pos: Mouse position (x, y) to determine which atom survives
        
        Returns:
            The surviving Atom object
        """
        if bond is None or not bond.alive:
            raise ValueError("Invalid bond (None or not alive)")
        
        # Get atom objects from bond (using object references, not indices!)
        atom_a = bond.a
        atom_b = bond.b
        
        debug_print(1, f"Collapsing bond {bond._id} between Atom({atom_a._id}) and Atom({atom_b._id})")
        
        pos_a = atom_a.pos
        pos_b = atom_b.pos
        center = (pos_a + pos_b) / 2.0
        
        # Determine which atom survives (farther from mouse)
        dist_a = np.linalg.norm(pos_a[:2] - mouse_pos)
        dist_b = np.linalg.norm(pos_b[:2] - mouse_pos)
        
        if dist_a > dist_b:
            # A survives, B is removed
            survivor = atom_a
            to_remove = atom_b
            debug_print(1, f"  Atom({atom_a._id}) survives (farther from mouse), Atom({atom_b._id}) removed")
        else:
            # B survives, A is removed
            survivor = atom_b
            to_remove = atom_a
            debug_print(1, f"  Atom({atom_b._id}) survives (farther from mouse), Atom({atom_a._id}) removed")
        
        # Move survivor to center (in the Atom object!)
        survivor.pos = center
        debug_print(1, f"  Moved survivor Atom({survivor._id}) to center {center[:2]}")
        
        # Mark original bond as dead
        bond.alive = False
        debug_print(1, f"  Marked bond {bond._id} as dead")
        
        # Transfer bonds from removed atom to survivor
        # Find all bonds involving the removed atom (excluding the bond we're collapsing)
        bonds_to_transfer = []
        for b in to_remove.bonds:
            if b is bond:
                debug_print(1, f"    Skipping bond {b._id} (is the bond being collapsed)")
                continue  # Skip the bond being collapsed
            if not b.alive:
                debug_print(1, f"    Skipping bond {b._id} (already dead)")
                continue  # Skip already dead bonds
            # Get the other atom in this bond
            other = b.other(to_remove)
            if other is not survivor:
                bonds_to_transfer.append((b, other))
                debug_print(1, f"  Will transfer bond {b._id} (to Atom({other._id})) to survivor")
            else:
                debug_print(1, f"    Skipping bond {b._id} (other atom is survivor)")
        
        # Mark bonds to removed atom as dead (they'll be recreated from survivor)
        for b, _ in bonds_to_transfer:
            b.alive = False
            debug_print(1, f"  Marked bond {b._id} as dead (will recreate)")
        
        # Create new bonds from survivor to other atoms
        for _, other in bonds_to_transfer:
            new_bond = self.graph.add_bond(survivor, other)
            debug_print(1, f"  Created new bond {new_bond._id} between survivor({survivor._id}) and Atom({other._id})")
        
        # Note: survivor's existing bonds (to atoms other than to_remove) are preserved
        # Only the bond between survivor and to_remove is marked dead
        
        # Remove H children of the atom being removed (they become orphans)
        h_children = self.graph.h_children(to_remove)
        for h in h_children:
            h.alive = False
            debug_print(1, f"  Marked H child Atom({h._id}) as dead (orphan)")
        
        # Mark removed atom as dead (soft delete)
        to_remove.alive = False
        debug_print(1, f"  Marked Atom({to_remove._id}) as dead")
        
        # Cleanup dead objects
        n_atoms, n_bonds, n_rings = self.graph.cleanup_invalid()
        debug_print(1, f"  Cleanup removed: {n_atoms} atoms, {n_bonds} bonds, {n_rings} rings")
        
        # Sync neighbor lists from bonds (derived data)
        self.graph.sync_neighbor_lists()
        
        # Sync sys from graph
        self._sync_sys()
        self._rings_dirty = True
        
        # Auto-adjust H if enabled (after topology operation completes)
        if self.auto_h_cap:
            self.adjust_h()
        
        return survivor

    def remove_h_caps(self):
        """Remove all H cap atoms from the graph (soft delete)."""
        for h in [a for a in list(self.graph.atoms.values()) if a.subtype == 'H_cap']:
            h.alive = False
        self.graph.cleanup_invalid()

    def add_h_caps(self):
        """Add H caps to undercoordinated heavy atoms based on current topology."""
        # Bond lengths for H caps
        bond_lengths = {'C': 1.09, 'N': 1.01, 'O': 0.97}
        added_h_atoms = []
        
        # Find undercoordinated heavy atoms and add H caps
        # Collect atoms first to avoid dict-changed-during-iteration error
        heavy_atoms = [a for a in self.graph.atoms.values() if a.alive and a.ename not in {'H', 'E'}]
        for a in heavy_atoms:
            if not a.alive or a.ename in {'H', 'E'}:
                continue
            
            npi = self._get_npi_from_subtype(a.subtype)
            target = self._target_sigma(a.ename, npi)
            
            # Count heavy atom neighbors (excluding H caps) using graph.neighbors
            heavy_neighbors = [n for n in a.neighbors if n.subtype != 'H_cap']
            current = len(heavy_neighbors)
            n_missing = target - current
            
            if n_missing <= 0:
                continue
            
            ename = a.ename
            bl = bond_lengths.get(ename, 1.09)
            pos_a = a.pos
            nb = len(heavy_neighbors)
            
            # Calculate directions using make_epair_geom logic (now using atom objects)
            directions = self._calc_h_directions_atom(a, npi, nb, heavy_neighbors)
            
            for direction in directions[:n_missing]:
                h_pos = pos_a + direction * bl
                # Add H cap to graph with explicit bond
                h_atom = self.graph.add_atom(h_pos, 'H', elements.ELEMENT_DICT['H'][0],
                                            pin=None, parent=a, subtype='H_cap')
                self.graph.add_bond(a, h_atom)  # Explicit bond between heavy atom and H
                added_h_atoms.append(h_atom)
                debug_print(1, f"  Added H cap to Atom({a._id}) at {h_pos[:2]}")
        
        # Sync neighbor lists to include new H atoms
        self.graph.sync_neighbor_lists()
        self._sync_sys()
        return added_h_atoms

    def adjust_h(self):
        """Add/remove H caps based on electron counting: nsigma = nvalence - npi - nepair.
        
        Uses graph.neighbors (derived from bond topology) for neighbor counting.
        - sp (npi=2, acetylene-like): 180° linear
        - sp3 (npi=0): 109.5° tetrahedral
        - sp2 (npi=1): 120° trigonal planar (O uses 109° water-like)
        
        NOTE: This uses the persistent bond topology from the graph, NOT distance-based detection.
        """
        self.remove_h_caps()
        self.add_h_caps()
        # Note: add_h_caps() already calls _sync_sys() and _rings_dirty = True

    def _calc_h_directions_atom(self, atom, npi, nb, heavy_neighbors):
        """Calculate H placement directions for 2D systems using atom objects.
        
        Args:
            atom: The heavy atom to place H caps on
            npi: Number of pi orbitals (0=sp3, 1=sp2, 2=sp)
            nb: Number of heavy atom neighbors
            heavy_neighbors: List of neighboring Atom objects (not H caps)
        
        Returns:
            List of direction vectors for H placement
        """
        pos = atom.pos
        if nb == 0:
            # No neighbors: use default directions based on hybridization
            if npi == 2:  # sp (linear in-plane)
                return [np.array([1.0, 0.0, 0.0]), np.array([-1.0, 0.0, 0.0])]
            elif npi == 0:  # sp3 (tetrahedral, but 2D: 3 in-plane directions)
                return [
                    np.array([1.0, 0.0, 0.0]),
                    np.array([-0.5, np.sqrt(3)/2, 0.0]),
                    np.array([-0.5, -np.sqrt(3)/2, 0.0])
                ]
            else:  # sp2 (trigonal planar in-plane)
                angles = [0, 2*np.pi/3, 4*np.pi/3]
                return [np.array([np.cos(a), np.sin(a), 0.0]) for a in angles]
        
        v1 = au.normalize(heavy_neighbors[0].pos - pos)
        if nb > 1: v2 = au.normalize(heavy_neighbors[1].pos - pos)
        if nb > 2: v3 = au.normalize(heavy_neighbors[2].pos - pos)
        
        z_hat = np.array([0.0, 0.0, 1.0])
        if npi == 0:  # sp3 (tetrahedral, 109.5° bond angles)
            if nb == 3:  # NH3-like - place H opposite to centroid of neighbors
                centroid = au.normalize(v1 + v2 + v3)
                return [-centroid]
            elif nb == 2:  # sp3 with 2 heavy neighbors: 2 Hs symmetric above/below plane
                bisect = au.normalize(v1 + v2)  # points toward neighbors avg
                # The 2 Hs lie in the plane containing z_hat and -bisect
                # component along bisect: -1/sqrt(3) (from tetrahedral constraint)
                # component along z_hat:  +/- sqrt(2/3) (above/below plane)
                c_along = -1.0/np.sqrt(3.0)  # ≈ -0.577
                c_z     =  np.sqrt(2.0/3.0)  # ≈ 0.816
                return [au.normalize(bisect*c_along + z_hat*c_z),
                        au.normalize(bisect*c_along - z_hat*c_z)]
            elif nb == 1:  # CH3-like: 3 Hs in tetrahedral cone around -v1
                # In-plane perp to v1
                perp_ip = np.array([-v1[1], v1[0], 0.0])
                if np.linalg.norm(perp_ip) < 1e-8: perp_ip = np.array([0.0, 1.0, 0.0])
                perp_ip = au.normalize(perp_ip)
                # Tetrahedral cone: dot(h, v1) = -1/3, transverse magnitude = 2√2/3
                cc = -1.0/3.0           # component along v1
                ct = 2.0*np.sqrt(2.0)/3.0  # transverse magnitude = sqrt(8/9) ≈ 0.9428
                # 3 Hs at azimuthal 0°,120°,240° around v1: (perp_ip, z_hat) basis
                return [
                    au.normalize(v1*cc + perp_ip*ct),                                        # in-plane
                    au.normalize(v1*cc + perp_ip*(ct*(-0.5)) + z_hat*(ct*(np.sqrt(3)/2))),  # above
                    au.normalize(v1*cc + perp_ip*(ct*(-0.5)) - z_hat*(ct*(np.sqrt(3)/2))),  # below
                ]
        elif npi == 1:  # sp2 (2D: H in-plane)
            if nb == 2:  # =N- like - H opposite to bisector in-plane
                bisect = au.normalize(v1 + v2)
                return [-bisect]
            elif nb == 1:  # =O like - 2 H at 120° in-plane
                perp = np.array([-v1[1], v1[0], 0.0])
                if np.linalg.norm(perp) < 1e-8: perp = np.array([0.0, 1.0, 0.0])
                perp = au.normalize(perp)
                # 2 H at +/-120° from -v1 in-plane
                return [
                    au.normalize(-v1 * 0.5 + perp * (np.sqrt(3)/2)),
                    au.normalize(-v1 * 0.5 - perp * (np.sqrt(3)/2))
                ]
        elif npi == 2:  # sp (linear in-plane)
            if nb == 1:
                return [au.normalize(pos - heavy_neighbors[0].pos)]
            else:
                return [au.normalize(-v1)]
        
        return []

    def _calc_h_directions(self, i, npi, nb, neighbors):
        """Calculate H placement directions for 2D systems (xy-plane, pi along z).
        
        Legacy function - now wraps _calc_h_directions_atom for compatibility.
        """
        # Create dummy Atom object with position from sys
        class DummyAtom:
            def __init__(self, pos, _id):
                self.pos = pos
                self._id = _id
        
        atom = DummyAtom(self.sys.apos[i], i)
        heavy_neighbors = [DummyAtom(self.sys.apos[j], j) for j in neighbors]
        return self._calc_h_directions_atom(atom, npi, nb, heavy_neighbors)

    def recalc_bonds(self, skip_sync=False):
        """Recompute bonds from distance threshold.
        
        Args:
            skip_sync: If True, skip _sync_sys() call (use when adjust_h will be called after)
        """
        self.graph.recalc_bonds(self.a_CC)
        self._rings_dirty = True
        if not skip_sync:
            self._sync_sys()
        if self.sys.apos is not None and len(self.sys.apos) > 0:
            self.sys.findBonds(Rcut=3.0, RvdwCut=1.2)
            self.sys.neighs()
        else:
            self.sys.bonds = None; self.sys.ngs = None

    def update_node_offset(self, node_key, offset):
        """Set absolute position of a pinned atom relative to its grid node."""
        a = self.graph._pin_to_atom.get(node_key)
        if a is not None and a.pin is not None:
            a.pos[0] = a.pin[0] + offset[0]
            a.pos[1] = a.pin[1] + offset[1]

    def snap_atoms_to_grid(self):
        """Snap all pinned atoms back to their grid node positions."""
        for a in self.graph.atoms.values():
            if a.pin is not None:
                a.pos[0] = a.pin[0]; a.pos[1] = a.pin[1]; a.pos[2] = 0.0
    
    def build_system(self):
        """Return the persistent AtomicSystem (sys is now the authoritative state)."""
        return self.sys
    
    def build_lattice_vectors(self):
        """Compute lattice vectors based on self.sys bounding box and PBC settings."""
        if self.sys.apos is None or len(self.sys.apos) == 0:
            return np.eye(3) * 20.0
        apos = self.sys.apos
        xmin, xmax = apos[:,0].min(), apos[:,0].max()
        ymin, ymax = apos[:,1].min(), apos[:,1].max()
        if self.pbc_x:
            s3 = np.sqrt(3.0)
            period = self.a_CC * s3
            span = xmax - xmin
            ncells = max(1, round(span / period))
            Lx = ncells * period
        else:
            Lx = (xmax - xmin) + self.vacuum_gap
        if self.pbc_y:
            span = ymax - ymin
            Ly = span + 2.0
        else:
            Ly = (ymax - ymin) + self.vacuum_gap
        Lz = self.vacuum_gap
        return np.array([[Lx, 0.0, 0.0], [0.0, Ly, 0.0], [0.0, 0.0, Lz]])
    
    def run_relaxation(self, workdir='kekule_relax', **kwargs):
        """Run DFTB+ geometry optimization on self.sys."""
        from . import dftb_utils
        lvs = self.build_lattice_vectors()
        if len(self.sys.apos) == 0:
            raise ValueError("Cannot relax an empty system.")
        # Store original CoG to undo centering shift after relaxation (only needed for non-PBC)
        cog_orig = self.sys.apos.mean(axis=0)
        if not (self.pbc_x or self.pbc_y):
            center = 0.5 * (lvs[0] + lvs[1] + lvs[2])
            self.sys.apos += (center - cog_orig)[None, :]
        enames = list(self.sys.enames)
        debug_print(1, f"DEBUG: run_relaxation starting with {len(self.sys.apos)} atoms: {dict(zip(*np.unique(enames, return_counts=True)))}")
        nk_x = max(1, int(8 / max(1, lvs[0,0] / 2.46))) if self.pbc_x else 1
        nk_y = max(1, int(8 / max(1, lvs[1,1] / 2.46))) if self.pbc_y else 1
        E, apos_out, forces = dftb_utils.run_pbc(
            self.sys.apos, enames, lvs,
            do_relax=True,
            nk=(nk_x, nk_y, 1),
            k_shift=(0.5 if self.pbc_x else 0.0, 0.5 if self.pbc_y else 0.0, 0.0),
            workdir=workdir,
            **kwargs
        )
        # Undo centering shift if we applied it
        if not (self.pbc_x or self.pbc_y):
            cog_relaxed = apos_out.mean(axis=0)
            self.sys.apos[:] = apos_out - (cog_relaxed - cog_orig)[None, :]
        else:
            self.sys.apos[:] = apos_out
        # Sync relaxed positions back to AtomicGraph (authoritative source)
        # This ensures topological operations don't lose the relaxed geometry
        self.graph.update_positions_from_array(self.sys.apos)
        return E, forces, lvs

    def get_xyz_string(self):
        """Return the current system as an XYZ formatted string."""
        s = f"{len(self.sys.apos)}\n"
        s += "Kekule Structure Explorer Export\n"
        for i in range(len(self.sys.apos)):
            p = self.sys.apos[i]
            s += f"{self.sys.enames[i]} {p[0]:10.5f} {p[1]:10.5f} {p[2]:10.5f}\n"
        return s

    def save_xyz(self, fname, comment=""):
        """Save the current system to an XYZ file with optional comment."""
        with open(fname, 'w') as f:
            f.write(f"{len(self.sys.apos)}\n")
            f.write(f"{comment}\n")
            for i, (pos, ename) in enumerate(zip(self.sys.apos, self.sys.enames)):
                f.write(f"{ename} {pos[0]:.6f} {pos[1]:.6f} {pos[2]:.6f}\n")

    def load_xyz(self, fname):
        """Load a system from an XYZ file and map it back to the grid."""
        self.graph = AtomicGraph()   # full reset
        from . import atomicUtils as au
        apos, Zs, es, qs, comment = au.load_xyz(fname)
        # Add heavy atoms snapped to grid
        for i, e in enumerate(es):
            if e in ('H', 'E'): continue
            nk = snap_to_grid(apos[i], self.a_CC)
            if nk not in self.graph._pin_to_atom:
                self.graph.add_atom(np.array([nk[0], nk[1], 0.0]),
                                    e, elements.ELEMENT_DICT[e][0],
                                    pin=nk, parent=None,
                                    subtype=self._get_element_default_subtype(e))
        # Create bonds for all heavy atoms (no recalc_bonds!)
        heavy_atoms = [a for a in self.graph.atoms.values() if a.alive and a.ename not in ('H', 'E')]
        for a in heavy_atoms:
            self._create_bond_to_nearest_heavy(a)
        self.graph.sync_neighbor_lists()
        self._guess_rings()
        # Add H atoms, finding nearest heavy atom by object ref
        heavy = [a for a in self.graph.atoms.values() if a.ename not in ('H', 'E')]
        for i, e in enumerate(es):
            if e != 'H': continue
            p = apos[i]
            best_d, best_a = float('inf'), None
            for a in heavy:
                d = float(np.linalg.norm(p - a.pos))
                if d < best_d: best_d = d; best_a = a
            if best_a is not None and best_d < 1.5:
                self.graph.add_atom(np.array(p, dtype=np.float64), 'H',
                                    elements.ELEMENT_DICT['H'][0],
                                    pin=None, parent=best_a, subtype='H_cap')
        # Sync after loading H atoms (bonds already created above)
        self.graph.sync_neighbor_lists()

    def _guess_rings(self):
        """Heuristic: infer ring axial coords from heavy atom grid positions."""
        n2a = self.graph._pin_to_atom
        for q in range(-20, 21):
            for r in range(-20, 21):
                if (q, r) in self.graph.rings: continue
                ring_nodes = honeycomb_ring_nodes(q, r, self.a_CC)
                ring_atom_objs = [n2a[snap_to_grid(nd, self.a_CC)]
                                  for nd in ring_nodes
                                  if snap_to_grid(nd, self.a_CC) in n2a]
                if len(ring_atom_objs) == 6:
                    self.graph.add_ring(q, r, ring_atom_objs)

    # ============ Ribbon construction (replaces GrapheneRibbonBuilder) ============

    def _set_atom_element(self, ia, element):
        """Set element of atom ia (int index into to_arrays() order)."""
        atom_list, *_ = self.graph.to_arrays()
        a = atom_list[ia]
        a.ename   = element
        a.atype   = elements.ELEMENT_DICT[element][0]
        a.subtype = self._get_element_default_subtype(element)

    def _add_passivation_h(self, ia, bond_length=1.09):
        """Add a single H atom in the missing bond direction of atom ia."""
        if self.sys.ngs is None:
            self.sys.neighs()
        pos = self.sys.apos[ia]
        heavy_neighs = [j for j in self.sys.ngs[ia] if self.sys.enames[j] not in ('H', 'E')]
        nb = len(heavy_neighs)
        direction = self.sys._missing_sp2_direction(ia, heavy_neighs, nb, 0)
        h_pos = pos + direction * bond_length
        self._append_atom(h_pos, 'H', parent=ia, subtype='H_cap')

    def _identify_edge_atoms(self):
        """Identify edge atoms of a ribbon patch.

        Returns
        -------
        zigzag_edge : list of int
            Atom indices on zigzag edges (top/bottom).
        armchair_edge : list of int
            Atom indices on armchair edges (left/right ends).
        """
        if self.sys.ngs is None:
            self.sys.neighs()

        edge_atoms = []
        for i in range(len(self.sys.apos)):
            e = self.sys.enames[i]
            if e in ('H', 'E'):
                continue
            heavy_neighs = [j for j in self.sys.ngs[i] if self.sys.enames[j] not in ('H', 'E')]
            if len(heavy_neighs) < 3:
                edge_atoms.append(i)

        if not edge_atoms:
            return [], []

        ys = self.sys.apos[edge_atoms, 1]
        xs = self.sys.apos[edge_atoms, 0]
        y_margin = self.a_CC * 0.8
        x_margin = self.a_CC * 1.0

        top_edge    = [edge_atoms[i] for i, y in enumerate(ys) if y > ys.max() - y_margin]
        bottom_edge = [edge_atoms[i] for i, y in enumerate(ys) if y < ys.min() + y_margin]
        left_edge   = [edge_atoms[i] for i, x in enumerate(xs) if x < xs.min() + x_margin]
        right_edge  = [edge_atoms[i] for i, x in enumerate(xs) if x > xs.max() - x_margin]

        # For periodic x, armchair edges are left/right (these wrap around)
        # Zigzag edges are top/bottom (these don't wrap)
        # Corner atoms (top-left, top-right, bottom-left, bottom-right) are in both
        # For periodic x: exclude ONLY the pure armchair edge atoms (not corners)
        zigzag_edge   = list(set(top_edge + bottom_edge))
        armchair_edge = list(set(left_edge + right_edge))
        return zigzag_edge, armchair_edge

    def _passivate_edges(self, passivation, bPeriodicX=False):
        """Passivate edge atoms of the current structure.

        Args
        ----
        passivation : str or None
            Passivation type ('N', 'NH', 'CH', 'H', 'O', 'C=O', 'C-OH') or None for no passivation.
        bPeriodicX : bool
            If True, exclude atoms at extreme x positions (armchair edges that wrap).
        """
        if passivation is None:
            return  # No passivation requested

        if passivation not in PASSIVATION_GROUPS:
            raise ValueError(f"Unknown passivation type: {passivation}")

        edge_atoms = []
        for i in range(len(self.sys.apos)):
            e = self.sys.enames[i]
            if e in ('H', 'E'):
                continue
            if self.sys.ngs is None:
                self.sys.neighs()
            heavy_neighs = [j for j in self.sys.ngs[i] if self.sys.enames[j] not in ('H', 'E')]
            if len(heavy_neighs) < 3:
                edge_atoms.append(i)

        if not edge_atoms:
            return

        # For periodic x, exclude atoms at extreme x positions (armchair edges)
        if bPeriodicX:
            xs = self.sys.apos[edge_atoms, 0]
            x_min, x_max = xs.min(), xs.max()
            x_margin = self.a_CC * 1.2
            edge_atoms = [edge_atoms[i] for i, x in enumerate(xs)
                          if x > x_min + x_margin and x < x_max - x_margin]

        if not edge_atoms:
            return

        y_center = self.sys.apos[:, 1].mean() if len(self.sys.apos) > 0 else 0.0

        # Remove any existing H atoms first (clean slate)
        h_indices = [i for i, e in enumerate(self.sys.enames) if e == 'H']
        if h_indices:
            self._rebuild_after_delete(h_indices)
            # Clean up and re-find edge atoms after H removal (no recalc_bonds!)
            self.graph.cleanup_invalid()
            self.graph.sync_neighbor_lists()
            edge_atoms = [i for i in range(len(self.sys.apos))
                          if self.sys.enames[i] not in ('H', 'E') and
                          len([j for j in self.sys.ngs[i] if self.sys.enames[j] not in ('H', 'E')]) < 3]
            if bPeriodicX:
                xs = self.sys.apos[edge_atoms, 0]
                x_min, x_max = xs.min(), xs.max()
                x_margin = self.a_CC * 0.5
                edge_atoms = [edge_atoms[i] for i, x in enumerate(xs)
                              if x > x_min + x_margin and x < x_max - x_margin]
            if not edge_atoms:
                return

        # Use data-driven PASSIVATION_GROUPS
        group = PASSIVATION_GROUPS[passivation]
        for ia in edge_atoms:
            if self.sys.enames[ia] != 'C':
                continue
            pos_C = self.sys.apos[ia]
            is_top = pos_C[1] > y_center
            direction = 1.0 if is_top else -1.0
            self._apply_passivation_group(ia, passivation, is_top)
        # Sync after passivation (no recalc_bonds!)
        self.graph.sync_neighbor_lists()

    def _apply_y_offsets(self, y_bottom_offset, y_top_offset):
        """Apply y offsets to top/bottom edge atoms."""
        if len(self.sys.apos) == 0:
            return
        ys = self.sys.apos[:, 1]
        y_margin = self.a_CC * 0.6
        if y_bottom_offset is not None:
            bottom_mask = ys < ys.min() + y_margin
            self.sys.apos[bottom_mask, 1] -= y_bottom_offset
        if y_top_offset is not None:
            top_mask = ys > ys.max() - y_margin
            self.sys.apos[top_mask, 1] += y_top_offset

    def build_zigzag_ribbon(self, width_chains, length_cells, passivation='N', passivation_bottom=None, passivation_top=None, start_with_A=True, y_top_offset=None, y_bottom_offset=None, scale_x=1.0, scale_y=1.0, bPeriodicX=False, side_passivation='CH'):
        """Build a zigzag graphene ribbon.

        For periodic x (bPeriodicX=True), uses strip-based construction to ensure
        armchair edges wrap correctly without passivation.
        For non-periodic x (bPeriodicX=False), uses ring-based construction with
        edge-specific passivation: top/bottom edges use 'passivation', side edges use 'side_passivation'.

        Parameters
        ----------
        width_chains : int
            Number of atom rows across the ribbon width.
        length_cells : int
            Number of unit cells along the ribbon length.
        passivation : str
            Passivation type for both top and bottom edges (used if passivation_bottom/passivation_top not specified).
        passivation_bottom : str or list of str
            Passivation type for bottom edge only.
        passivation_top : str or list of str
            Passivation type for top edge only.
        bPeriodicX : bool
            If True, ribbon is periodic along x (armchair edges wrap, no side passivation).
        side_passivation : str
            Passivation type for side edges when bPeriodicX=False (default: 'CH').
        """
        # Handle passivation parameters
        if passivation_bottom is None and passivation_top is None:
            passivation_bottom = passivation
            passivation_top = passivation
        elif passivation_bottom is None:
            passivation_bottom = passivation
        elif passivation_top is None:
            passivation_top = passivation

        if bPeriodicX:
            # PBC mode: use strip-based construction
            self._build_strip_ribbon(width_chains, length_cells, passivation_bottom, passivation_top,
                                     start_with_A, y_top_offset, y_bottom_offset, scale_x, scale_y)
        else:
            # Non-periodic mode: build using rings, then apply selective passivation
            self._build_ribbon_from_rings(width_chains, length_cells, start_with_A)
            # Apply top/bottom passivation only (not to side edges)
            self._passivate_edges_top_bottom_only_separate(passivation_bottom, passivation_top)
            # Apply side passivation if specified
            if side_passivation:
                self._apply_side_passivation_single_ribbon(side_passivation)

        # Apply anisotropic scaling if requested
        if scale_x != 1.0 or scale_y != 1.0:
            self.sys.apos[:, 0] *= scale_x
            self.sys.apos[:, 1] *= scale_y

        return self

    def _apply_side_passivation_single_ribbon(self, passivation):
        """Apply passivation to side edges (left/right) of a single ribbon.

        Side edges are identified as atoms with < 3 heavy neighbors at extreme x positions.
        Passivation groups are applied with x-direction orientation (H points left/right).
        """
        if self.sys.ngs is None:
            self.sys.neighs()

        # Find edge atoms (less than 3 heavy neighbors)
        edge_atoms = []
        for i in range(len(self.sys.apos)):
            if self.sys.enames[i] in ('H', 'E'):
                continue
            heavy_neighs = [j for j in self.sys.ngs[i] if self.sys.enames[j] not in ('H', 'E')]
            if len(heavy_neighs) < 3:
                edge_atoms.append(i)

        if not edge_atoms:
            return

        # Identify side edges (extreme x positions)
        xs = self.sys.apos[edge_atoms, 0]
        x_min, x_max = xs.min(), xs.max()
        x_margin = self.a_CC * 0.5

        # Filter to side edges only (left/right) and determine direction
        side_atoms = []
        for ia, x in zip(edge_atoms, xs):
            if x < x_min + x_margin:
                side_atoms.append((ia, -1.0))  # Left edge, point left
            elif x > x_max - x_margin:
                side_atoms.append((ia, 1.0))   # Right edge, point right

        if not side_atoms:
            return

        # Apply passivation to side atoms using PASSIVATION_GROUPS
        if passivation not in PASSIVATION_GROUPS:
            raise ValueError(f"Unknown passivation type: {passivation}")

        group = PASSIVATION_GROUPS[passivation]
        for ia, direction in side_atoms:
            pos_C = self.sys.apos[ia]
            # Process each atom in the group (swap x and y for side edges)
            for i, (elem, x, y, z) in enumerate(group):
                if i == 0 and x == 0.0 and y == 0.0 and z == 0.0:
                    # First atom at origin replaces the C atom
                    self._set_atom_element(ia, elem)
                else:
                    # Add atom at C_pos + coords (swap x and y for side edges)
                    # y becomes x-direction for side edges
                    pos_new = [pos_C[0] + y * direction, pos_C[1] + x, pos_C[2] + z]
                    a = self._append_atom(pos_new, elem, subtype='H_cap')
                    # Create bond to nearest heavy atom
                    self._create_bond_to_nearest_heavy(a)

        self.graph.sync_neighbor_lists()

    def _passivate_edges_top_bottom_only_separate(self, passivation_bottom, passivation_top):
        """Passivate only top/bottom edge atoms (zigzag edges), not side edges (armchair edges).
        
        Allows different passivation for top vs bottom edges.
        Supports list-based passivation where each element is applied to successive edge sites.
        """
        if self.sys.ngs is None:
            self.sys.neighs()

        # Convert strings to single-element lists for uniform handling
        if isinstance(passivation_bottom, str):
            passivation_bottom = [passivation_bottom]
        if isinstance(passivation_top, str):
            passivation_top = [passivation_top]

        # Find edge atoms (less than 3 heavy neighbors)
        edge_atoms = []
        for i in range(len(self.sys.apos)):
            if self.sys.enames[i] in ('H', 'E'):
                continue
            heavy_neighs = [j for j in self.sys.ngs[i] if self.sys.enames[j] not in ('H', 'E')]
            if len(heavy_neighs) < 3:
                edge_atoms.append(i)

        if not edge_atoms:
            return

        # Separate edge atoms by position
        y_center = self.sys.apos[:, 1].mean()
        xs = self.sys.apos[edge_atoms, 0]
        x_min, x_max = xs.min(), xs.max()
        x_margin = self.a_CC * 0.5

        bottom_edge_atoms = []
        top_edge_atoms = []
        
        for ia in edge_atoms:
            y = self.sys.apos[ia, 1]
            x = self.sys.apos[ia, 0]
            
            # Skip side edge atoms (extreme x positions)
            if x < x_min + x_margin or x > x_max - x_margin:
                continue
            
            # Classify as top or bottom edge
            if y < y_center:
                bottom_edge_atoms.append(ia)
            else:
                top_edge_atoms.append(ia)

        # Sort by x position to apply passivation in order along the edge
        bottom_edge_atoms.sort(key=lambda ia: self.sys.apos[ia, 0])
        top_edge_atoms.sort(key=lambda ia: self.sys.apos[ia, 0])

        # Apply passivation using PASSIVATION_GROUPS
        for idx, ia in enumerate(bottom_edge_atoms):
            p = passivation_bottom[idx % len(passivation_bottom)] if passivation_bottom else None
            if p and p in PASSIVATION_GROUPS:
                self._apply_passivation_group(ia, p, is_top=False)
        
        for idx, ia in enumerate(top_edge_atoms):
            p = passivation_top[idx % len(passivation_top)] if passivation_top else None
            if p and p in PASSIVATION_GROUPS:
                self._apply_passivation_group(ia, p, is_top=True)
        # Sync after passivation (no recalc_bonds!)
        self.graph.sync_neighbor_lists()

    def _passivate_edges_top_bottom_only(self, passivation):
        """Passivate only top/bottom edge atoms (zigzag edges), not side edges (armchair edges)."""
        if passivation is None:
            return

        if passivation not in PASSIVATION_GROUPS:
            raise ValueError(f"Unknown passivation type: {passivation}")

        if self.sys.ngs is None:
            self.sys.neighs()

        # Find edge atoms (less than 3 heavy neighbors)
        edge_atoms = []
        for i in range(len(self.sys.apos)):
            if self.sys.enames[i] in ('H', 'E'):
                continue
            heavy_neighs = [j for j in self.sys.ngs[i] if self.sys.enames[j] not in ('H', 'E')]
            if len(heavy_neighs) < 3:
                edge_atoms.append(i)

        if not edge_atoms:
            return

        # Filter to top/bottom edges only (exclude extreme x positions = side edges)
        xs = self.sys.apos[edge_atoms, 0]
        x_min, x_max = xs.min(), xs.max()
        x_margin = self.a_CC * 0.5

        top_bottom_edge_atoms = [edge_atoms[i] for i, x in enumerate(xs)
                                if x > x_min + x_margin and x < x_max - x_margin]

        if not top_bottom_edge_atoms:
            return

        y_center = self.sys.apos[:, 1].mean() if len(self.sys.apos) > 0 else 0.0

        # Remove any existing H atoms first (clean slate)
        h_indices = [i for i, e in enumerate(self.sys.enames) if e == 'H']
        if h_indices:
            self._rebuild_after_delete(h_indices)
            # Clean up and re-find edge atoms after H removal (no recalc_bonds!)
            self.graph.cleanup_invalid()
            self.graph.sync_neighbor_lists()
            # Re-find edge atoms after deletion
            edge_atoms = []
            for i in range(len(self.sys.apos)):
                if self.sys.enames[i] in ('H', 'E'):
                    continue
                heavy_neighs = [j for j in self.sys.ngs[i] if self.sys.enames[j] not in ('H', 'E')]
                if len(heavy_neighs) < 3:
                    edge_atoms.append(i)
            xs = self.sys.apos[edge_atoms, 0]
            top_bottom_edge_atoms = [edge_atoms[i] for i, x in enumerate(xs)
                                    if x > x_min + x_margin and x < x_max - x_margin]

        # Apply passivation using PASSIVATION_GROUPS
        for ia in top_bottom_edge_atoms:
            if self.sys.enames[ia] != 'C':
                continue
            pos_C = self.sys.apos[ia]
            is_top = pos_C[1] > y_center
            self._apply_passivation_group(ia, passivation, is_top)
        # Sync after passivation (no recalc_bonds!)
        self.graph.sync_neighbor_lists()

    def _build_ribbon_from_rings(self, width_chains, length_cells, start_with_A):
        """Build ribbon using ring-based construction (for non-PBC mode)."""
        for r in range(width_chains):
            for q in range(length_cells):
                self.add_ring(q, r)

    def _build_strip_ribbon(self, width_chains, length_cells, passivation_bottom, passivation_top, start_with_A, y_top_offset, y_bottom_offset, scale_x, scale_y):
        """Strip-based construction for periodic ribbons (mirrors GrapheneRibbonBuilder logic).
        
        For PBC along x, only creates atoms within one periodicity (x < x_periodicity).
        Supports list-based passivation for variable passivation per site along edges.
        """
        L = self.a_CC
        xa = L * np.cos(np.pi / 6)
        ya = L * np.sin(np.pi / 6)
        yb = L

        xa *= scale_x
        ya *= scale_y
        yb *= scale_y

        x_periodicity = 2 * xa

        strip_types = []
        for row in range(width_chains):
            if start_with_A:
                row_mod = row % 4
                is_A_strip = (row_mod == 0) or (row_mod == 3)
            else:
                row_mod = row % 4
                is_A_strip = (row_mod == 2) or (row_mod == 3)
            strip_types.append(is_A_strip)

        y_positions = [0.0]
        for r in range(1, width_chains):
            prev_is_A = strip_types[r-1]
            curr_is_A = strip_types[r]
            if prev_is_A and not curr_is_A:
                y_positions.append(y_positions[-1] + ya)
            elif not prev_is_A and not curr_is_A:
                y_positions.append(y_positions[-1] + yb)
            elif not prev_is_A and curr_is_A:
                y_positions.append(y_positions[-1] + ya)
            else:
                y_positions.append(y_positions[-1] + yb)

        # Build atoms - create all atoms (DFT code handles PBC wrapping)
        self.sys.bonds = np.empty((0, 2), dtype=np.int32)
        
        for row in range(width_chains):
            is_A_strip = strip_types[row]
            y = y_positions[row]
            x_shift = 0.0 if is_A_strip else xa

            for i in range(length_cells):
                x = i * x_periodicity + x_shift
                self._append_atom([x, y, 0.0], 'C', subtype='C_sp2')

                # Add bonds to previous row
                if row > 0:
                    prev_row_start = (row - 1) * length_cells
                    prev_is_A = strip_types[row - 1]
                    atom_idx = len(self.sys.apos) - 1
                    prev_idx = prev_row_start + i
                    self.sys.bonds = np.append(self.sys.bonds, [[prev_idx, atom_idx]], axis=0)

                    if is_A_strip and not prev_is_A:
                        if i > 0:
                            prev_idx = prev_row_start + (i - 1)
                            self.sys.bonds = np.append(self.sys.bonds, [[prev_idx, atom_idx]], axis=0)
                    elif not is_A_strip and prev_is_A:
                        if i < length_cells - 1:
                            prev_idx = prev_row_start + (i + 1)
                            self.sys.bonds = np.append(self.sys.bonds, [[prev_idx, atom_idx]], axis=0)

        # Apply y offsets
        if y_bottom_offset is not None:
            for i in range(length_cells):
                self.sys.apos[i, 1] -= y_bottom_offset

        if y_top_offset is not None:
            start_idx = (width_chains - 1) * length_cells
            for i in range(length_cells):
                self.sys.apos[start_idx + i, 1] += y_top_offset

        # Passivate top and bottom rows only (zigzag edges)
        self._strip_passivate(width_chains, length_cells, passivation_bottom, passivation_top)

    def _strip_passivate(self, width_chains, length_cells, passivation_bottom, passivation_top):
        """Passivate only top and bottom rows for strip-based construction (zigzag edges).

        For periodic x, passivate ALL atoms in top/bottom rows (zigzag edges).
        The armchair edges (left/right) wrap in PBC, so they are not separate edges.
        
        Supports list-based passivation where each element is applied to successive sites.
        """
        # Convert strings to single-element lists for uniform handling
        if isinstance(passivation_bottom, str):
            passivation_bottom = [passivation_bottom]
        if isinstance(passivation_top, str):
            passivation_top = [passivation_top]

        # Bottom row (row 0)
        for i in range(length_cells):
            pb = passivation_bottom[i % len(passivation_bottom)] if passivation_bottom else None
            if pb is not None:
                self._strip_passivate_atom(i, pb, is_top=False)

        # Top row (last row)
        start_idx = (width_chains - 1) * length_cells
        for i in range(length_cells):
            pt = passivation_top[i % len(passivation_top)] if passivation_top else None
            if pt is not None:
                self._strip_passivate_atom(start_idx + i, pt, is_top=True)

        self.sys.neighs()

    def _strip_passivate_atom(self, ia, passivation, is_top):
        """Passivate a single edge atom using data-driven PASSIVATION_GROUPS.

        Args:
            ia: atom index to passivate
            passivation: passivation type (string)
            is_top: True for top edge, False for bottom edge
        """
        if passivation not in PASSIVATION_GROUPS:
            raise ValueError(f"Unknown passivation type: {passivation}")

        group = PASSIVATION_GROUPS[passivation]
        pos_C = self.sys.apos[ia]
        direction = 1.0 if is_top else -1.0

        # Process each atom in the group
        for i, (elem, x, y, z) in enumerate(group):
            if i == 0 and x == 0.0 and y == 0.0 and z == 0.0:
                # First atom at origin replaces the C atom
                self._set_atom_element(ia, elem)
            else:
                # Add atom at C_pos + scaled coords
                pos_new = [pos_C[0] + x, pos_C[1] + y * direction, pos_C[2] + z]
                self._append_atom(pos_new, elem, subtype='H_cap')

    def combine_ribbons(self, backend1, backend2, L_Hb=2.0, shift_x=0.0):
        """Combine two single ribbons into a two-ribbon system with hydrogen-bond gap.

        Parameters
        ----------
        backend1 : KekuleBackend
            Bottom ribbon.
        backend2 : KekuleBackend
            Top ribbon.
        L_Hb : float
            Hydrogen-bond separation between ribbons (Angstrom).
        shift_x : float
            Relative shift along x (fraction of Lx).

        Returns
        -------
        self (for chaining)
        """
        apos_N = backend1.sys.apos.copy()
        apos_NH = backend2.sys.apos.copy()
        enames_N = backend1.sys.enames
        enames_NH = backend2.sys.enames

        # Center y positions
        apos_N[:, 1] -= apos_N[:, 1].mean()
        apos_NH[:, 1] -= apos_NH[:, 1].mean()

        # Apply x shift
        Lx = apos_N[:, 0].max() - apos_N[:, 0].min()
        apos_NH[:, 0] += shift_x * Lx

        # Position top ribbon above bottom ribbon with hydrogen-bond gap
        y_max_N = np.max(apos_N[:, 1])
        y_min_NH = np.min(apos_NH[:, 1])
        apos_NH[:, 1] += (y_max_N + L_Hb) - y_min_NH

        # Combine into this backend
        self.__init__(a_CC=self.a_CC)
        added_atoms = []
        for pos, ename in zip(apos_N, enames_N):
            a = self._append_atom(pos, ename)
            added_atoms.append(a)
        for pos, ename in zip(apos_NH, enames_NH):
            a = self._append_atom(pos, ename)
            added_atoms.append(a)
        # Create bonds for all added atoms (no recalc_bonds!)
        for a in added_atoms:
            self._create_bond_to_nearest_heavy(a)
        self.graph.sync_neighbor_lists()
        return self

    def build_two_ribbon_cell(self, width_chains=4, length_cells=1, Lx=2.4, L_Hb=2.0, shift_x=0.0, 
                           bottom_passivation='N', top_passivation='NH', bPeriodicX=True, side_passivation='CH'):
        """Build a cell with two ribbons separated by a hydrogen-bond gap.

        Composable approach: generates two single ribbons separately, then combines them.

        Parameters
        ----------
        width_chains : int
            Number of atom rows per ribbon.
        length_cells : int
            Number of unit cells along ribbon length.
        Lx : float
            Target periodic length along x (Angstrom).
        L_Hb : float
            Hydrogen-bond separation between ribbons (Angstrom).
        shift_x : float
            Relative shift along x (fraction of Lx).
        bottom_passivation : str or list of str
            Passivation type(s) for bottom ribbon (both edges for single ribbon).
        top_passivation : str or list of str
            Passivation type(s) for top ribbon (both edges for single ribbon).
        bPeriodicX : bool
            If True, ribbon is periodic along x (armchair edges wrap, no side passivation).
        side_passivation : str
            Passivation type for side edges when bPeriodicX=False (default: 'CH').

        Returns
        -------
        self (for chaining)
        """
        # Convert single passivation to list for uniform handling
        if isinstance(bottom_passivation, str):
            bottom_passivation = [bottom_passivation]
        if isinstance(top_passivation, str):
            top_passivation = [top_passivation]
        
        # Ensure both lists have the same length
        if len(bottom_passivation) != len(top_passivation):
            raise ValueError(f"bottom_passivation and top_passivation lists must have the same length: {len(bottom_passivation)} vs {len(top_passivation)}")
        
        # If lists are provided, use length to determine length_cells
        if len(bottom_passivation) > 1:
            length_cells = len(bottom_passivation)

        # Build bottom ribbon
        backend1 = KekuleBackend(a_CC=self.a_CC)
        backend1.build_zigzag_ribbon(width_chains, length_cells, passivation_bottom=bottom_passivation, passivation_top=bottom_passivation,
                                     scale_x=Lx / (2.0 * self.a_CC * np.cos(np.pi / 6)), bPeriodicX=bPeriodicX, side_passivation=side_passivation if not bPeriodicX else None)

        # Build top ribbon
        backend2 = KekuleBackend(a_CC=self.a_CC)
        backend2.build_zigzag_ribbon(width_chains, length_cells, passivation_bottom=top_passivation, passivation_top=top_passivation,
                                     scale_x=Lx / (2.0 * self.a_CC * np.cos(np.pi / 6)), bPeriodicX=bPeriodicX, side_passivation=side_passivation if not bPeriodicX else None)

        # Combine ribbons
        self.combine_ribbons(backend1, backend2, L_Hb=L_Hb, shift_x=shift_x)

        return self

    def _apply_edge_passivation(self, bottom_passivation, top_passivation, L_Hb):
        """Apply passivation to top/bottom edges of two-ribbon system.

        Bottom ribbon gets bottom_passivation on its bottom edge only.
        Top ribbon gets top_passivation on its top edge only.
        Side edges are NOT passivated here (handled by _apply_side_passivation).
        
        If passivation is a list, each element is applied to successive edge sites along x.
        """
        if self.sys.ngs is None:
            self.sys.neighs()

        # Convert to lists for uniform handling
        if isinstance(bottom_passivation, str):
            bottom_passivation = [bottom_passivation]
        if isinstance(top_passivation, str):
            top_passivation = [top_passivation]

        # Find edge atoms (less than 3 heavy neighbors)
        edge_atoms = []
        for i in range(len(self.sys.apos)):
            if self.sys.enames[i] in ('H', 'E'):
                continue
            heavy_neighs = [j for j in self.sys.ngs[i] if self.sys.enames[j] not in ('H', 'E')]
            if len(heavy_neighs) < 3:
                edge_atoms.append(i)

        if not edge_atoms:
            return

        # Separate edge atoms by position
        y_center = self.sys.apos[:, 1].mean()
        xs = self.sys.apos[edge_atoms, 0]
        x_min, x_max = xs.min(), xs.max()
        x_margin = self.a_CC * 0.5

        bottom_edge_atoms = []
        top_edge_atoms = []
        
        for ia in edge_atoms:
            y = self.sys.apos[ia, 1]
            x = self.sys.apos[ia, 0]
            
            # Skip side edge atoms (extreme x positions)
            if x < x_min + x_margin or x > x_max - x_margin:
                continue
            
            # Classify as top or bottom edge
            if y < y_center:
                bottom_edge_atoms.append(ia)
            else:
                top_edge_atoms.append(ia)

        # Sort edge atoms by x position to apply passivation in order
        bottom_edge_atoms.sort(key=lambda ia: self.sys.apos[ia, 0])
        top_edge_atoms.sort(key=lambda ia: self.sys.apos[ia, 0])

        # Apply passivation using PASSIVATION_GROUPS
        for idx, ia in enumerate(bottom_edge_atoms):
            passivation = bottom_passivation[idx % len(bottom_passivation)]
            if passivation and passivation in PASSIVATION_GROUPS:
                self._apply_passivation_group(ia, passivation, is_top=False)
        
        for idx, ia in enumerate(top_edge_atoms):
            passivation = top_passivation[idx % len(top_passivation)]
            if passivation and passivation in PASSIVATION_GROUPS:
                self._apply_passivation_group(ia, passivation, is_top=True)

    def _apply_passivation_group(self, ia, passivation, is_top):
        """Apply passivation group to a single edge atom."""
        if passivation not in PASSIVATION_GROUPS:
            return
        
        group = PASSIVATION_GROUPS[passivation]
        pos_C = self.sys.apos[ia]
        direction = 1.0 if is_top else -1.0

        for i, (elem, x, y, z) in enumerate(group):
            if i == 0 and x == 0.0 and y == 0.0 and z == 0.0:
                self._set_atom_element(ia, elem)
            else:
                pos_new = [pos_C[0] + x, pos_C[1] + y * direction, pos_C[2] + z]
                self._append_atom(pos_new, elem, subtype='H_cap')

    def _apply_side_passivation(self, passivation):
        """Apply passivation to side edges (left/right) of the combined two-ribbon system.

        Side edges are identified as atoms with < 3 heavy neighbors at extreme x positions.
        Passivation groups are applied with x-direction orientation (H points left/right).
        """
        if self.sys.ngs is None:
            self.sys.neighs()

        # Find edge atoms (less than 3 heavy neighbors)
        edge_atoms = []
        for i in range(len(self.sys.apos)):
            if self.sys.enames[i] in ('H', 'E'):
                continue
            heavy_neighs = [j for j in self.sys.ngs[i] if self.sys.enames[j] not in ('H', 'E')]
            if len(heavy_neighs) < 3:
                edge_atoms.append(i)

        if not edge_atoms:
            return

        # Identify side edges (extreme x positions)
        xs = self.sys.apos[edge_atoms, 0]
        x_min, x_max = xs.min(), xs.max()
        x_margin = self.a_CC * 0.5

        # Filter to side edges only (left/right) and determine direction
        side_atoms = []
        for ia, x in zip(edge_atoms, xs):
            if x < x_min + x_margin:
                side_atoms.append((ia, -1.0))  # Left edge, point left
            elif x > x_max - x_margin:
                side_atoms.append((ia, 1.0))   # Right edge, point right

        if not side_atoms:
            return

        # Apply passivation to side atoms using PASSIVATION_GROUPS
        if passivation not in PASSIVATION_GROUPS:
            raise ValueError(f"Unknown passivation type: {passivation}")

        group = PASSIVATION_GROUPS[passivation]
        for ia, direction in side_atoms:
            pos_C = self.sys.apos[ia]
            # Process each atom in the group (swap x and y for side edges)
            for i, (elem, x, y, z) in enumerate(group):
                if i == 0 and x == 0.0 and y == 0.0 and z == 0.0:
                    # First atom at origin replaces the C atom
                    self._set_atom_element(ia, elem)
                else:
                    # Add atom at C_pos + coords (swap x and y for side edges)
                    # y becomes x-direction for side edges
                    pos_new = [pos_C[0] + y * direction, pos_C[1] + x, pos_C[2] + z]
                    self._append_atom(pos_new, elem, subtype='H_cap')

    def report_state(self):
        """Print summary of the backend state for debugging."""
        debug_print(2, "=== KekuleBackend State ===")
        debug_print(2, f"  rings={self.rings}")
        debug_print(2, f"  natoms={len(self.sys.apos)}")
        if len(self.sys.apos) > 0:
            elems = dict(zip(*np.unique(self.sys.enames, return_counts=True)))
            debug_print(2, f"  elements={elems}")
        debug_print(2, f"  n_pinned={sum(1 for p in self.atom_pin if p is not None)}")
        debug_print(2, f"  n_hydrogens={sum(1 for e in self.sys.enames if e == 'H')}")
        debug_print(2, f"  n_bonds={len(self.sys.bonds) if self.sys.bonds is not None else 0}")


# ============ Module-level convenience functions ============

def build_ribbon(passivation, width_chains, length_cells, Lx, a_CC=1.42):
    """Build a ribbon and return arrays (mirrors deprecated GrapheneRibbonBuilder.build_ribbon API).

    Parameters
    ----------
    passivation : str
        Edge passivation type ('N', 'NH', 'CH', 'H', 'O', 'C=O', 'C-OH').
    width_chains : int
        Number of atom rows across the ribbon width.
    length_cells : int
        Number of unit cells along the ribbon length.
    Lx : float
        Target periodic length along x (Angstrom).
    a_CC : float
        C-C bond length (Angstrom).

    Returns
    -------
    pos2d : np.ndarray, shape (n_atoms, 2)
        2D atom positions.
    atypes : np.ndarray, dtype int32
        Atomic numbers.
    elems : list of str
        Element symbols.
    """
    backend = KekuleBackend(a_CC=a_CC)
    xa_nom = a_CC * np.cos(np.pi / 6)
    scale_x = Lx / (2.0 * xa_nom)
    backend.build_zigzag_ribbon(width_chains=width_chains, length_cells=length_cells,
                                  passivation=passivation, scale_x=scale_x, bPeriodicX=False)
    elems = list(backend.sys.enames)
    atypes = backend.sys.atypes
    pos2d = backend.sys.apos[:, :2].copy()
    return pos2d, atypes, elems


def build_two_ribbon_cell(width_chains=4, length_cells=1, Lx=2.4, a_CC=1.42, L_Hb=2.0, shift_x=0.0):
    """Build a cell with two ribbons separated by a hydrogen-bond gap.

    Mirrors the deprecated GrapheneRibbonBuilder.build_two_ribbon_cell API.

    Parameters
    ----------
    width_chains : int
        Number of atom rows per ribbon.
    length_cells : int
        Number of unit cells along ribbon length.
    Lx : float
        Target periodic length along x (Angstrom).
    a_CC : float
        C-C bond length (Angstrom).
    L_Hb : float
        Hydrogen-bond separation between ribbons (Angstrom).
    shift_x : float
        Relative shift along x (fraction of Lx).

    Returns
    -------
    apos : np.ndarray, shape (n_atoms, 3)
        3D atom positions.
    atypes : np.ndarray, dtype int32
        Atomic numbers.
    elems : list of str
        Element symbols.
    lvs : np.ndarray, shape (3, 3)
        Lattice vectors.
    """
    backend1 = KekuleBackend(a_CC=a_CC)
    backend2 = KekuleBackend(a_CC=a_CC)
    xa_nom = a_CC * np.cos(np.pi / 6)
    scale_x = Lx / (2.0 * xa_nom)
    backend1.build_zigzag_ribbon(width_chains, length_cells, passivation='N', scale_x=scale_x, bPeriodicX=False)
    backend2.build_zigzag_ribbon(width_chains, length_cells, passivation='NH', scale_x=scale_x, bPeriodicX=False)

    apos_N = backend1.sys.apos.copy()
    apos_NH = backend2.sys.apos.copy()
    apos_N[:, 1]  -= apos_N[:, 1].mean()
    apos_NH[:, 1] -= apos_NH[:, 1].mean()
    apos_NH[:, 0] += shift_x * Lx

    y_max_N  = np.max(apos_N[:, 1])
    y_min_NH = np.min(apos_NH[:, 1])
    apos_NH[:, 1] += (y_max_N + L_Hb) - y_min_NH

    y_span_N  = np.max(apos_N[:, 1])  - np.min(apos_N[:, 1])
    y_span_NH = np.max(apos_NH[:, 1]) - np.min(apos_NH[:, 1])
    Ly = y_span_N + y_span_NH + 2 * L_Hb

    apos   = np.vstack([apos_N, apos_NH])
    atypes = np.concatenate([backend1.sys.atypes, backend2.sys.atypes])
    elems  = list(backend1.sys.enames) + list(backend2.sys.enames)

    apos[:, 2] = 0.0
    apos[:, 1] -= apos[:, 1].mean()
    Lz = 20.0
    apos[:, 2] += 0.5 * Lz
    lvs = np.array([[Lx, 0.0, 0.0], [0.0, Ly, 0.0], [0.0, 0.0, Lz]])
    return apos, atypes, elems, lvs
