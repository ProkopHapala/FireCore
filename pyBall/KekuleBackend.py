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
    
    def _element_target_valence(self, element, npi):
        """Calculate target sigma bonds using electron counting: nvalence = nsigma + npi + nepair.
        
        For aromatic sp2 systems:
        - C: needs 3 sigma bonds (nvalence=4, nepair=0, npi=1) -> nsigma=3
        - N: pyrrolic needs 3 sigma (nvalence=5, nepair=1 in pi, npi=1) -> nsigma=3
             pyridinic needs 2 sigma (nvalence=5, nepair=2, npi=1) -> nsigma=2
        - O: ketone needs 1 sigma (nvalence=6, nepair=2, npi=1) -> nsigma=3? No, ketone=1
             hydroxyl needs 2 sigma (nvalence=6, nepair=2, npi=1) -> nsigma=3? No, hydroxyl=2
        
        Returns: target number of sigma bonds (nsigma)
        """
        # Use subtype-based logic for aromatic systems (more robust)
        subtype = f"{element}_sp2" if npi == 1 else f"{element}_sp3"
        if 'pyridinic' in subtype:
            return 2
        elif 'pyrrolic' in subtype or 'N_sp2' in subtype:
            return 3
        elif 'ketone' in subtype or 'O_sp2' in subtype:
            return 1
        elif 'hydroxyl' in subtype:
            return 2
        elif 'C_sp2' in subtype:
            return 3
        elif 'C_sp3' in subtype:
            return 4
        elif 'N_sp3' in subtype:
            return 3
        elif 'O_sp3' in subtype:
            return 2
        return 3
    
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
            # Try to find H atom near this position
            nk = snap_to_grid(node_key, self.a_CC)
            nk_3d = np.array([nk[0], nk[1], 0.0])
            found_h = False
            for i, (pin, parent) in enumerate(zip(self.atom_pin, self.atom_parent)):
                if pin is None and parent is not None and self.sys.enames[i] == 'H':
                    p = self.sys.apos[i]
                    if np.linalg.norm(p - nk_3d) < 0.5:
                        # Replace H with new element at grid position
                        self.sys.enames[i] = element
                        self.sys.atypes[i] = elements.ELEMENT_DICT[element][0]
                        self.atom_pin[i] = nk
                        self.atom_parent[i] = None
                        self.atom_subtype[i] = self._get_element_default_subtype(element)
                        self._rebuild_ring_atoms()
                        found_h = True
                        break
            # If no H found and node is empty, add new atom
            if not found_h:
                self._append_atom(
                    pos=[nk[0], nk[1], 0.0],
                    ename=element,
                    pin=nk,
                    parent=None,
                    subtype=self._get_element_default_subtype(element)
                )
                self._rebuild_ring_atoms()
                self.recalc_bonds()
                self.adjust_h()

    def set_atom_valency(self, node_key, npi):
        """Set the number of pi bonds for an atom (npi).
        
        Electron counting: nvalence = nsigma + npi + nepair
        - nvalence: 8 for N,O (octet), 4 for C
        - nsigma: number of sigma bonds (neighbors)
        - nepair: number of lone pairs (1 for N, 2 for O, 0 for C)
        - npi: number of pi bonds (0=sp3, 1=sp2, 2=sp)
        
        Setting npi changes the target valence and triggers H passivation.
        """
        if node_key not in self.node_to_atom:
            return
        ia = self.node_to_atom[node_key]
        e = self.sys.enames[ia]
        if e == 'H':
            return
        # Store npi in subtype as f"{element}_sp{npi+3}" for compatibility
        # but the actual logic should use the electron counting
        sp_map = {0: 'sp3', 1: 'sp2', 2: 'sp'}
        self.atom_subtype[ia] = f"{e}_{sp_map.get(npi, 'sp2')}"
        # Trigger H passivation to adjust to new valency
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

    def toggle_h_state(self, node_key):
        """Toggle H-passivation state for an atom at a node (N or O) using electron counting."""
        if node_key not in self.node_to_atom:
            return
        ia = self.node_to_atom[node_key]
        e = self.sys.enames[ia]
        subtype = self.atom_subtype[ia]
        npi = self._get_npi_from_subtype(subtype)
        
        if e == 'N':
            # Toggle between pyridinic (npi=1, target sigma=2) and pyrrolic (npi=1, target sigma=3)
            if 'pyridinic' in subtype:
                self.atom_subtype[ia] = 'N_pyrrolic'
            else:
                self.atom_subtype[ia] = 'N_pyridinic'
        elif e == 'O':
            # Toggle between ketone (npi=1, target sigma=1) and hydroxyl (npi=1, target sigma=2)
            if 'ketone' in subtype:
                self.atom_subtype[ia] = 'O_hydroxyl'
            else:
                self.atom_subtype[ia] = 'O_ketone'
        
        # Rebuild H passivation with new subtype
        self.adjust_h()

    def adjust_h(self):
        """Rebuild H passivation from scratch using electron counting."""
        if self.sys.ngs is None: self.sys.neighs()
        
        # Remove all existing H_cap atoms
        h_indices = [i for i, st in enumerate(self.atom_subtype) if st == 'H_cap']
        if h_indices:
            self._rebuild_after_delete(h_indices)
            self.recalc_bonds()
        
        # Build target valence dict using electron counting
        tv = {}
        for i, subtype in enumerate(self.atom_subtype):
            if self.sys.enames[i] in {'H', 'E'}: continue
            element = self.sys.enames[i]
            npi = self._get_npi_from_subtype(subtype)
            tv[i] = self._element_target_valence(element, npi)
        
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

    def save_xyz(self, fname):
        """Save the current system to an XYZ file."""
        with open(fname, 'w') as f:
            f.write(self.get_xyz_string())

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
