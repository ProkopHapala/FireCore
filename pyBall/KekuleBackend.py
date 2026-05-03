"""
KekuleBackend.py - Backend logic for the Kekule Structure Explorer
===================================================================
Manages hexagonal grid state, converts to AtomicSystem, handles passivation
and DFTB+ relaxation. Reuses GrapheneRibbonBuilder geometry and AtomicSystem topology.
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
    """Manages the hexagonal grid state and generates AtomicSystem instances.
    
    State is stored as:
        self.rings: set of (q, r) tuples - which hex cells have benzene rings
        self.nodes: dict mapping (rx, ry) rounded Cartesian -> element name ('C','N','O')
        self.a_CC: bond length
    """
    
    def __init__(self, a_CC=1.42):
        self.a_CC = a_CC
        self.rings = set()       # set of (q, r) axial coords
        self.nodes = {}          # {(rx,ry): 'C'|'N'|'O'}
        self.node_valences = {}  # {(rx,ry): target_valence}
        self.node_offsets = {}   # {(rx,ry): np.array([dx, dy])}
        self.pbc_x = False
        self.pbc_y = False
        self.vacuum_gap = 15.0   # Angstrom vacuum gap for non-periodic directions
    
    def toggle_ring(self, q, r):
        """Add or remove a benzene ring at axial position (q, r)."""
        key = (q, r)
        if key in self.rings:
            self.rings.discard(key)
            # Remove nodes that belong ONLY to this ring (not shared with other rings)
            ring_nodes = honeycomb_ring_nodes(q, r, self.a_CC)
            for node in ring_nodes:
                nk = snap_to_grid(node, self.a_CC)
                # Check if any other ring shares this node
                shared = False
                for oq, oR in self.rings:
                    other_nodes = honeycomb_ring_nodes(oq, oR, self.a_CC)
                    for on in other_nodes:
                        if snap_to_grid(on, self.a_CC) == nk:
                            shared = True; break
                    if shared: break
                if not shared:
                    self.nodes.pop(nk, None)
                    self.node_valences.pop(nk, None)
        else:
            self.rings.add(key)
            ring_nodes = honeycomb_ring_nodes(q, r, self.a_CC)
            for node in ring_nodes:
                nk = snap_to_grid(node, self.a_CC)
                if nk not in self.nodes:
                    self.nodes[nk] = 'C'  # default element
    
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

    def set_atom(self, node_key, element, target_valence=None):
        """Set or change the element at a node position. element='C','N','O' or None to remove."""
        if element is None:
            self.nodes.pop(node_key, None)
            self.node_valences.pop(node_key, None)
        else:
            self.nodes[node_key] = element
            if target_valence is not None:
                self.node_valences[node_key] = target_valence
            elif element == 'N':
                self.node_valences[node_key] = 3
            elif element == 'O':
                self.node_valences[node_key] = 2

    def remove_atom(self, node_key):
        self.set_atom(node_key, None)

    def remove_ring(self, q, r):
        key = (q, r)
        if key in self.rings:
            self.toggle_ring(q, r)

    def toggle_h_state(self, key):
        """Toggle H state based on valence."""
        if key not in self.nodes: return
        e = self.nodes[key]
        if e == 'N':
            v = self.node_valences.get(key, 3)
            self.node_valences[key] = 3 if v == 2 else 2
        elif e == 'O':
            v = self.node_valences.get(key, 2)
            self.node_valences[key] = 2 if v == 1 else 1

    def update_node_offset(self, node_key, offset):
        """Store a manual displacement for a node."""
        self.node_offsets[node_key] = np.asarray(offset)[:2]

    def reset_offsets(self):
        """Clear all manual displacements."""
        self.node_offsets = {}
    
    def build_system(self, bPassivate=True):
        """Convert current grid state to an AtomicSystem with bonds and optional H passivation.
        
        Returns:
            AtomicSystem instance with 3D positions (z=0), bonds, and optional H atoms
        """
        if not self.nodes:
            return AtomicSystem(apos=np.zeros((0,3)), atypes=np.array([],dtype=np.int32), enames=np.array([],dtype=object))
        
        # Build position and element arrays from node dict
        node_keys = sorted(self.nodes.keys())  # deterministic order
        natoms = len(node_keys)
        apos = np.zeros((natoms, 3))
        enames = []
        atypes = []
        key_to_idx = {}
        for i, nk in enumerate(node_keys):
            apos[i, 0] = nk[0]
            apos[i, 1] = nk[1]
            # Apply offset if present
            if nk in self.node_offsets:
                off = self.node_offsets[nk]
                apos[i, 0] += off[0]
                apos[i, 1] += off[1]
            
            enames.append(self.nodes[nk])
            atypes.append(elements.ELEMENT_DICT[self.nodes[nk]][0])
            key_to_idx[nk] = i
        
        enames = np.array(enames, dtype=object)
        atypes = np.array(atypes, dtype=np.int32)
        
        sys = AtomicSystem(apos=apos, atypes=atypes, enames=enames)
        sys.findBonds(Rcut=3.0, RvdwCut=1.2)
        sys.neighs()
        
        if bPassivate:
            # Build target_valence dict for AtomicSystem
            # Default to sp2 valence: C=3, N=3, O=2
            # If PBC is active, C might be 4
            v_def = {'C': 3, 'N': 3, 'O': 2}
            if self.pbc_x or self.pbc_y:
                v_def['C'] = 4
                
            tv = {}
            for i, k in enumerate(node_keys):
                el = self.nodes[k]
                # Use per-node valence if available, else element default
                tv[i] = self.node_valences.get(k, v_def.get(el, 3))
            
            sys.add_capping_h_sp2(target_valence=tv)
            sys.neighs()
        
        print(f"DEBUG: build_system nodes={len(self.nodes)} total_atoms={len(sys.apos)} elements={dict(zip(*np.unique(sys.enames, return_counts=True)))}")
        return sys
    
    def build_lattice_vectors(self, sys):
        """Compute lattice vectors based on system bounding box and PBC settings.
        
        Returns:
            (3,3) numpy array of lattice vectors
        """
        if sys.apos is None or len(sys.apos) == 0:
            return np.eye(3) * 20.0
        
        xmin, xmax = sys.apos[:,0].min(), sys.apos[:,0].max()
        ymin, ymax = sys.apos[:,1].min(), sys.apos[:,1].max()
        
        if self.pbc_x:
            # For periodic ribbons, Lx should be a multiple of the zigzag period
            s3 = np.sqrt(3.0)
            period = self.a_CC * s3
            span = xmax - xmin
            ncells = max(1, round(span / period))
            Lx = ncells * period
        else:
            Lx = (xmax - xmin) + self.vacuum_gap
        
        if self.pbc_y:
            span = ymax - ymin
            Ly = span + 2.0  # small buffer for periodic y
        else:
            Ly = (ymax - ymin) + self.vacuum_gap
        
        Lz = self.vacuum_gap
        return np.array([[Lx, 0.0, 0.0], [0.0, Ly, 0.0], [0.0, 0.0, Lz]])
    
    def run_relaxation(self, sys=None, workdir='kekule_relax', **kwargs):
        """Run DFTB+ geometry optimization on the current system.
        
        Args:
            sys: AtomicSystem to relax (if None, builds from current grid state)
            workdir: working directory for DFTB+ files
            **kwargs: passed to dftb_utils.run_pbc
        Returns:
            (E, apos_relaxed, forces)
        """
        from . import dftb_utils
        if sys is None:
            sys = self.build_system(bPassivate=True)
        
        lvs = self.build_lattice_vectors(sys)
        
        # Center atoms in the cell
        if len(sys.apos) > 0:
            center = 0.5 * (lvs[0] + lvs[1] + lvs[2])
            cog = sys.apos.mean(axis=0)
            sys.apos += (center - cog)[None, :]
        else:
            raise ValueError("Cannot relax an empty system.")
        
        enames = list(sys.enames)
        print(f"DEBUG: run_relaxation starting with {len(sys.apos)} atoms: {dict(zip(*np.unique(enames, return_counts=True)))}")
        nk_x = max(1, int(8 / max(1, lvs[0,0] / 2.46))) if self.pbc_x else 1
        nk_y = max(1, int(8 / max(1, lvs[1,1] / 2.46))) if self.pbc_y else 1
        
        E, apos_out, forces = dftb_utils.run_pbc(
            sys.apos, enames, lvs,
            do_relax=True,
            nk=(nk_x, nk_y, 1),
            k_shift=(0.5 if self.pbc_x else 0.0, 0.5 if self.pbc_y else 0.0, 0.0),
            workdir=workdir,
            **kwargs
        )
        return E, apos_out, forces, lvs

    def get_xyz_string(self, sys=None, bPassivate=True):
        """Return the current system as an XYZ formatted string."""
        if sys is None:
            sys = self.build_system(bPassivate=bPassivate)
        s = f"{len(sys.apos)}\n"
        s += "Kekule Structure Explorer Export\n"
        for i in range(len(sys.apos)):
            p = sys.apos[i]
            s += f"{sys.enames[i]} {p[0]:10.5f} {p[1]:10.5f} {p[2]:10.5f}\n"
        return s

    def save_xyz(self, fname, sys=None, bPassivate=True):
        """Save the current system to an XYZ file."""
        with open(fname, 'w') as f:
            f.write(self.get_xyz_string(sys, bPassivate))

    def load_xyz(self, fname):
        """Load a system from an XYZ file and try to map it back to the grid."""
        self.nodes = {}
        self.node_offsets = {}
        self.rings = set()
        
        from . import atomicUtils as au
        apos, es, _, _ = au.loadAtomsXYZ(fname)
        for i, e in enumerate(es):
            if e in ('H', 'E'): continue
            p = apos[i]
            nk = snap_to_grid(p, self.a_CC)
            self.nodes[nk] = e
            off = p[:2] - np.array(nk)
            if np.linalg.norm(off) > 1e-4:
                self.node_offsets[nk] = off
