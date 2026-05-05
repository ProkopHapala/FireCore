"""
KekuleBackend.py - Backend logic for the Kekule Structure Explorer
===================================================================
Manages persistent AtomicSystem state with hexagonal grid metadata.
AtomicSystem is the authoritative store; grid mappings are auxiliary.
Each UI action triggers exactly one mutation of the backend state.
"""

import numpy as np
from . import atomicUtils as au
from . import elements
from .AtomicSystem import AtomicSystem

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

class KekuleBackend:
    """Manages persistent AtomicSystem with hexagonal grid metadata.
    
    Persistent State:
        self.sys          = AtomicSystem()  # Authoritative: apos, enames, atypes, bonds, ngs
        self.atom_pin     = []            # For each atom i: grid node (rx,ry) or None
        self.atom_parent  = []            # For each atom i: parent heavy atom index or None
        self.atom_subtype = []            # For each atom i: e.g. 'C_sp2', 'N_pyridinic', 'H_cap'
        self.rings        = set()         # Set of (q,r) axial coords
        self.ring_atoms   = {}            # Dict {(q,r): [atom_indices]}
        self.node_to_atom = {}            # Dict {(rx,ry): atom_index}
    """
    
    def __init__(self, a_CC=1.42):
        self.a_CC = a_CC
        self.sys = AtomicSystem(apos=np.zeros((0,3)), atypes=np.array([],dtype=np.int32), enames=np.array([],dtype=object))
        self.atom_pin = []
        self.atom_parent = []
        self.atom_subtype = []
        self.rings = set()
        self.ring_atoms = {}
        self.node_to_atom = {}
        self.pbc_x = False
        self.pbc_y = False
        self.vacuum_gap = 15.0
    
    # --- Internal helpers ---
    
    def _append_atom(self, pos, ename, pin=None, parent=None, subtype=''):
        """Append a new atom to self.sys and tracking arrays."""
        atype = elements.ELEMENT_DICT[ename][0]
        if self.sys.apos is None or len(self.sys.apos) == 0:
            self.sys.apos = np.array([pos])
            self.sys.atypes = np.array([atype], dtype=np.int32)
            self.sys.enames = np.array([ename], dtype=object)
            self.sys.natoms = 1
        else:
            self.sys.apos = np.append(self.sys.apos, [pos], axis=0)
            self.sys.atypes = np.append(self.sys.atypes, atype)
            self.sys.enames = np.append(self.sys.enames, [ename])
            self.sys.natoms += 1
        self.atom_pin.append(pin)
        self.atom_parent.append(parent)
        self.atom_subtype.append(subtype)
        return len(self.sys.apos) - 1
    
    def _rebuild_ring_atoms(self):
        """Rebuild node_to_atom and ring_atoms from atom_pin and rings."""
        self.node_to_atom = {}
        for i, pin in enumerate(self.atom_pin):
            if pin is not None:
                self.node_to_atom[pin] = i
        
        self.ring_atoms = {}
        for q, r in self.rings:
            ring_nodes = honeycomb_ring_nodes(q, r, self.a_CC)
            atom_set = set()
            for node in ring_nodes:
                nk = snap_to_grid(node, self.a_CC)
                if nk in self.node_to_atom:
                    atom_set.add(self.node_to_atom[nk])
            self.ring_atoms[(q, r)] = sorted(list(atom_set))
    
    def _rebuild_after_delete(self, to_remove):
        """Remove atoms from sys and rebuild all mapping arrays."""
        if not to_remove:
            return
        
        to_remove = sorted(list(set(to_remove)))
        
        # Build old->new index mapping
        old_to_new = {}
        shift = 0
        for i in range(len(self.sys.apos)):
            if i in to_remove:
                shift += 1
                old_to_new[i] = None
            else:
                old_to_new[i] = i - shift
        
        # Delete from AtomicSystem
        self.sys.delete_atoms(to_remove)
        
        # Rebuild tracking arrays
        new_atom_pin = []
        new_atom_parent = []
        new_atom_subtype = []
        for i in range(len(self.atom_pin)):
            if i not in to_remove:
                new_atom_pin.append(self.atom_pin[i])
                parent = self.atom_parent[i]
                if parent is not None:
                    parent = old_to_new.get(parent)
                new_atom_parent.append(parent)
                new_atom_subtype.append(self.atom_subtype[i])
        
        self.atom_pin = new_atom_pin
        self.atom_parent = new_atom_parent
        self.atom_subtype = new_atom_subtype
        
        # Rebuild derived mappings
        self._rebuild_ring_atoms()
    
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
        """Add a benzene ring at axial position (q, r). Idempotent if already present."""
        key = (q, r)
        if key in self.rings:
            return
        
        ring_nodes = honeycomb_ring_nodes(q, r, self.a_CC)
        added_indices = []
        
        for node in ring_nodes:
            nk = snap_to_grid(node, self.a_CC)
            if nk not in self.node_to_atom:
                ia = self._append_atom(
                    pos=[nk[0], nk[1], 0.0],
                    ename='C',
                    pin=nk,
                    parent=None,
                    subtype='C_sp2'
                )
                added_indices.append(ia)
                self.node_to_atom[nk] = ia
        
        self.rings.add(key)
        self._rebuild_ring_atoms()
        
        if added_indices:
            self.recalc_bonds()
    
    def remove_ring(self, q, r):
        """Remove a benzene ring at axial position (q, r). Shared atoms are preserved."""
        key = (q, r)
        if key not in self.rings:
            return
        
        ring_atoms = self.ring_atoms.get(key, [])
        to_remove = set()
        
        for ia in ring_atoms:
            shared = False
            for other_key, other_atoms in self.ring_atoms.items():
                if other_key == key:
                    continue
                if ia in other_atoms:
                    shared = True
                    break
            if not shared:
                to_remove.add(ia)
                for j, parent in enumerate(self.atom_parent):
                    if parent == ia:
                        to_remove.add(j)
        
        self.rings.discard(key)
        self._rebuild_after_delete(to_remove)
        self.recalc_bonds()
    
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

    def snap_to_node(self, x, y, tol=0.5):
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
        if node_key in self.node_to_atom:
            ia = self.node_to_atom[node_key]
            self.sys.enames[ia] = element
            self.sys.atypes[ia] = elements.ELEMENT_DICT[element][0]
            self.atom_subtype[ia] = self._get_element_default_subtype(element)
        else:
            self._append_atom(
                pos=[node_key[0], node_key[1], 0.0],
                ename=element,
                pin=node_key,
                parent=None,
                subtype=self._get_element_default_subtype(element)
            )
            self._rebuild_ring_atoms()
            self.recalc_bonds()
            self.adjust_h()

    def set_atom_valency(self, node_key, npi):
        """Set npi (pi bond count) for atom at node_key. npi in {0,1,2}."""
        ia = self.node_to_atom.get(node_key)
        if ia is None: return
        e = self.sys.enames[ia]
        sp_map = {0: 'sp3', 1: 'sp2', 2: 'sp'}
        self.atom_subtype[ia] = f"{e}_{sp_map.get(npi, 'sp2')}"
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
        self.recalc_bonds()
        self.adjust_h()
        return ia

    def remove_atom(self, node_key):
        """Remove the atom at a grid node and any attached H atoms."""
        if node_key not in self.node_to_atom:
            return
        ia = self.node_to_atom[node_key]
        to_remove = {ia}
        for j, parent in enumerate(self.atom_parent):
            if parent == ia:
                to_remove.add(j)
        # Remove rings that would become empty or invalid
        rings_to_remove = []
        for key, atoms in self.ring_atoms.items():
            remaining = [a for a in atoms if a not in to_remove]
            if len(remaining) < len(atoms):
                ring_nodes = honeycomb_ring_nodes(key[0], key[1], self.a_CC)
                has_all_nodes = True
                for node in ring_nodes:
                    nk = snap_to_grid(node, self.a_CC)
                    if nk not in self.node_to_atom or self.node_to_atom[nk] in to_remove:
                        has_all_nodes = False
                        break
                if not has_all_nodes:
                    rings_to_remove.append(key)
        for key in rings_to_remove:
            self.rings.discard(key)
        self._rebuild_after_delete(to_remove)
        self.recalc_bonds()

    def remove_ring(self, q, r):
        """Remove a benzene ring at axial position (q, r). Shared atoms are preserved."""
        key = (q, r)
        if key not in self.rings:
            return
        ring_atoms = self.ring_atoms.get(key, [])
        to_remove = set()
        for ia in ring_atoms:
            shared = False
            for other_key, other_atoms in self.ring_atoms.items():
                if other_key == key:
                    continue
                if ia in other_atoms:
                    shared = True
                    break
            if not shared:
                to_remove.add(ia)
                for j, parent in enumerate(self.atom_parent):
                    if parent == ia:
                        to_remove.add(j)
        self.rings.discard(key)
        self._rebuild_after_delete(to_remove)
        self.recalc_bonds()

    def adjust_h(self):
        """Add/remove H caps based on electron counting: nsigma = nvalence - npi - nepair."""
        if self.sys.ngs is None: self.sys.neighs()
        # Remove existing H caps
        h_indices = [i for i, st in enumerate(self.atom_subtype) if st == 'H_cap']
        if h_indices:
            self._rebuild_after_delete(h_indices)
            self.recalc_bonds()
        # Build target sigma count per atom
        tv = {}
        for i, subtype in enumerate(self.atom_subtype):
            e = self.sys.enames[i]
            if e in {'H', 'E'}: continue
            npi = self._get_npi_from_subtype(subtype)
            tv[i] = self._target_sigma(e, npi)
        # Add H where needed
        added = self.sys.add_capping_h_sp2(target_valence=tv)
        for h_idx in added:
            parent = None
            if self.sys.bonds is not None:
                for b in self.sys.bonds:
                    if b[0] == h_idx: parent = b[1]; break
                    elif b[1] == h_idx: parent = b[0]; break
            self.atom_pin.append(None); self.atom_parent.append(parent); self.atom_subtype.append('H_cap')
        self.sys.neighs()
        return added

    def recalc_bonds(self):
        """Recompute bonds and neighbor lists from current positions."""
        if len(self.sys.apos) > 0:
            self.sys.findBonds(Rcut=3.0, RvdwCut=1.2)
            self.sys.neighs()
        else:
            self.sys.bonds = None; self.sys.ngs = None

    def update_node_offset(self, node_key, offset):
        """Set absolute position of a pinned atom relative to its grid node."""
        if node_key in self.node_to_atom:
            ia = self.node_to_atom[node_key]
            pin = self.atom_pin[ia]
            if pin is not None:
                self.sys.apos[ia, 0] = pin[0] + offset[0]
                self.sys.apos[ia, 1] = pin[1] + offset[1]

    def snap_atoms_to_grid(self):
        """Snap all pinned atoms back to their grid node positions."""
        for i, pin in enumerate(self.atom_pin):
            if pin is not None:
                self.sys.apos[i, 0] = pin[0]
                self.sys.apos[i, 1] = pin[1]
                self.sys.apos[i, 2] = 0.0
    
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
        print(f"DEBUG: run_relaxation starting with {len(self.sys.apos)} atoms: {dict(zip(*np.unique(enames, return_counts=True)))}")
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
        self.rings = set(); self.ring_atoms = {}; self.node_to_atom = {}
        self.atom_pin = []; self.atom_parent = []; self.atom_subtype = []
        from . import atomicUtils as au
        apos, Zs, es, qs, comment = au.load_xyz(fname)
        heavy_atoms = {}
        for i, e in enumerate(es):
            if e in ('H', 'E'): continue
            p = apos[i]
            nk = snap_to_grid(p, self.a_CC)
            heavy_atoms[nk] = (i, e)
        for nk, (orig_idx, e) in heavy_atoms.items():
            ia = self._append_atom(
                pos=[nk[0], nk[1], 0.0],
                ename=e,
                pin=nk,
                parent=None,
                subtype=self._get_element_default_subtype(e)
            )
            self.node_to_atom[nk] = ia
        self.recalc_bonds()
        # Guess rings from topology
        self._guess_rings()
        # Map H atoms
        for i, e in enumerate(es):
            if e != 'H': continue
            p = apos[i]
            # Find nearest heavy atom
            best_d, best_ia = float('inf'), None
            for j in range(len(self.sys.apos)):
                if self.sys.enames[j] in ('H', 'E'): continue
                d = np.linalg.norm(apos[i] - self.sys.apos[j])
                if d < best_d:
                    best_d = d; best_ia = j
            if best_ia is not None and best_d < 1.5:
                self._append_atom(
                    pos=p,
                    ename='H',
                    pin=None,
                    parent=best_ia,
                    subtype='H_cap'
                )
        self.recalc_bonds()

    def _guess_rings(self):
        """Heuristic: infer ring axial coords from heavy atom grid positions."""
        seen_nodes = set()
        for nk in self.node_to_atom:
            seen_nodes.add(nk)
        s3 = np.sqrt(3.0)
        for q in range(-20, 21):
            for r in range(-20, 21):
                ring_nodes = honeycomb_ring_nodes(q, r, self.a_CC)
                nodes_ok = True
                for node in ring_nodes:
                    nk = snap_to_grid(node, self.a_CC)
                    if nk not in seen_nodes:
                        nodes_ok = False; break
                if nodes_ok:
                    self.rings.add((q, r))
        self._rebuild_ring_atoms()

    # ============ Ribbon construction (replaces GrapheneRibbonBuilder) ============

    def _set_atom_element(self, ia, element):
        """Set element of atom ia, updating sys and atom_subtype tracking."""
        self.sys.enames[ia] = element
        self.sys.atypes[ia] = elements.ELEMENT_DICT[element][0]
        self.atom_subtype[ia] = self._get_element_default_subtype(element)

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
            self.recalc_bonds()
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
        self.recalc_bonds()

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
                    self._append_atom(pos_new, elem, subtype='H_cap')

        self.recalc_bonds()

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
        self.recalc_bonds()

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
            self.recalc_bonds()
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
        self.recalc_bonds()

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
        for pos, ename in zip(apos_N, enames_N):
            self._append_atom(pos, ename)
        for pos, ename in zip(apos_NH, enames_NH):
            self._append_atom(pos, ename)

        self.recalc_bonds()
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
        print("=== KekuleBackend State ===")
        print(f"  rings={self.rings}")
        print(f"  natoms={len(self.sys.apos)}")
        if len(self.sys.apos) > 0:
            elems = dict(zip(*np.unique(self.sys.enames, return_counts=True)))
            print(f"  elements={elems}")
        print(f"  n_pinned={sum(1 for p in self.atom_pin if p is not None)}")
        print(f"  n_hydrogens={sum(1 for e in self.sys.enames if e == 'H')}")
        print(f"  n_bonds={len(self.sys.bonds) if self.sys.bonds is not None else 0}")


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
