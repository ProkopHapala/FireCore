"""
GrapheneRibbonBuilder.py - Build periodic graphene ribbons with zigzag edges
===========================================================================
Generates graphene nanoribbons with zigzag edges and various passivation options.
"""

import numpy as np
import sys
import os

class GrapheneRibbonBuilder:
    """Builds graphene ribbons with zigzag edges and passivation."""
    
    def __init__(self, a_CC=1.42, lattice_constant=None):
        """
        Initialize ribbon builder.
        
        Args:
            a_CC: C-C bond length in Angstroms (default 1.42)
            lattice_constant: Lattice constant along ribbon direction (None = auto from a_CC)
        """
        self.a_CC = a_CC
        if lattice_constant is None:
            # For zigzag ribbon, lattice constant along ribbon is sqrt(3) * a_CC
            self.lattice_constant = np.sqrt(3) * a_CC
        else:
            self.lattice_constant = lattice_constant
        
        # Atomic data
        self.elements = []
        self.positions = []
        self.bonds = []  # List of (i, j) atom index pairs
        
    def build_zigzag_ribbon(self, width_chains, length_cells, passivation='N', 
                           start_with_A=True, y_top_offset=None, y_bottom_offset=None,
                           scale_x=1.0, scale_y=1.0):
        """
        Build a zigzag graphene ribbon.
        Systematic approach: periodic A/B strips along x-axis.
        
        Args:
            width_chains: Number of rows (ribbon width)
            length_cells: Number of atoms per strip (ribbon length)
            passivation: Passivation type ('N', 'NH', 'CH', 'H', 'O', 'C=O', 'C-OH')
            start_with_A: If True, start with A strip, else start with B strip
            y_top_offset: Override y-distance for top passivation strip (None = use default)
            y_bottom_offset: Override y-distance for bottom passivation strip (None = use default)
            scale_x: Scaling factor for x-axis (anisotropic strain)
            scale_y: Scaling factor for y-axis (anisotropic strain)
        """
        self.elements = []
        self.positions = []
        self.bonds = []  # List of (i, j) atom index pairs
        
        # Step 1: Simple hexagonal lattice with A/B strips
        # Hexagon oriented with tips up/down along y-axis
        # Pattern: A (tip), B (2 atoms), B (2 atoms), A (tip)
        # B-B strips: bonded vertically, y-distance = L
        # A-B strips: bonded by tilted bond, decomposed into:
        #   ya = L * sin(30°) (vertical component)
        #   xa = L * cos(30°) (horizontal component, also phase shift)
        # Periodicity along x in each strip: 2*xa
        # NO bonds within x-strips, all bonds are between strips
        L = self.a_CC  # bond length
        xa = L * np.cos(np.pi/6)  # L * cos(30°) = horizontal component
        ya = L * np.sin(np.pi/6)  # L * sin(30°) = vertical component
        yb = L  # vertical distance between B-B strips
        
        # Apply scaling factors
        xa *= scale_x
        ya *= scale_y
        yb *= scale_y
        
        x_periodicity = 2 * xa  # atoms in same strip are 2*xa apart
        
        # Pattern: A, B, B, A, A, B, B, A, ...
        # If start_with_A=False, shift the pattern: B, B, A, A, B, B, A, A, ...
        
        # First determine strip type for each row
        strip_types = []  # True for A, False for B
        for row in range(width_chains):
            if start_with_A:
                # Pattern: A(0), B(1), B(2), A(3), A(4), B(5), B(6), A(7), ...
                row_mod = row % 4
                is_A_strip = (row_mod == 0) or (row_mod == 3)
            else:
                # Pattern: B(0), B(1), A(2), A(3), B(4), B(5), A(6), A(7), ...
                row_mod = row % 4
                is_A_strip = (row_mod == 2) or (row_mod == 3)
            strip_types.append(is_A_strip)
        
        # Now calculate y positions based on strip type transitions
        y_positions = [0.0]  # First row at y=0
        for r in range(1, width_chains):
            prev_is_A = strip_types[r-1]
            curr_is_A = strip_types[r]
            
            if prev_is_A and not curr_is_A:
                # A to B: tilted bond
                y_positions.append(y_positions[-1] + ya)
            elif not prev_is_A and not curr_is_A:
                # B to B: vertical bond
                y_positions.append(y_positions[-1] + yb)
            elif not prev_is_A and curr_is_A:
                # B to A: tilted bond
                y_positions.append(y_positions[-1] + ya)
            else:  # prev_is_A and curr_is_A
                # A to A: vertical bond
                y_positions.append(y_positions[-1] + yb)
        
        # Generate atoms
        for row in range(width_chains):
            is_A_strip = strip_types[row]
            y = y_positions[row]
            
            # x shift: A and B strips are shifted by xa
            x_shift = 0.0 if is_A_strip else xa
            
            # Add atoms along x-axis
            # Periodicity is 2*xa (NO bonds within strip)
            for i in range(length_cells):
                x = i * x_periodicity + x_shift
                atom_idx = len(self.elements)
                self.elements.append('C')
                self.positions.append([x, y])
                
                # Add bonds to atoms in previous row based on strip pattern
                if row > 0:
                    prev_row_start = (row - 1) * length_cells
                    prev_is_A = strip_types[row - 1]
                    
                    # Base bond: index i to index i
                    prev_idx = prev_row_start + i
                    self.bonds.append([prev_idx, atom_idx])
                    
                    # Additional bond at AB junctions
                    if is_A_strip and not prev_is_A:
                        # A→B junction: add A[i] to B[i-1]
                        if i > 0:
                            prev_idx = prev_row_start + (i - 1)
                            self.bonds.append([prev_idx, atom_idx])
                    elif not is_A_strip and prev_is_A:
                        # B→A junction: add B[i] to A[i+1]
                        if i < length_cells - 1:
                            prev_idx = prev_row_start + (i + 1)
                            self.bonds.append([prev_idx, atom_idx])
        
        # Step 2: Replace top and bottom rows with passivation element
        # Bond length table for chemical groups (in Angstroms)
        bond_lengths = {
            'C-H': 1.09,    # C-H bond
            'C-N': 1.47,    # C-N bond
            'N-H': 1.01,    # N-H bond
            'C-O': 1.43,    # C-O single bond
            'C=O': 1.23,    # C=O double bond
            'O-H': 0.97,    # O-H bond
        }
        
        # Determine passivation element and add additional atoms if needed
        if passivation == 'N':
            pass_elem = 'N'
        elif passivation == 'NH':
            pass_elem = 'N'  # NH uses N atom, will add H later
        elif passivation in ['CH', 'H']:
            pass_elem = 'C'  # CH uses C atom, will add H later
        elif passivation == 'O':
            pass_elem = 'O'
        elif passivation == 'C=O':
            pass_elem = 'C'  # C=O uses C atom, will add O later
        elif passivation == 'C-OH':
            pass_elem = 'C'  # C-OH uses C atom, will add O and H later
        else:
            pass_elem = 'N'  # default
        
        # Bottom row (row 0) - apply y_bottom_offset if specified
        if len(self.elements) >= length_cells:
            for i in range(length_cells):
                idx = i
                self.elements[idx] = pass_elem
            if y_bottom_offset is not None:
                # Adjust y positions of bottom row
                for i in range(length_cells):
                    self.positions[i][1] -= y_bottom_offset
            
            # Add additional atoms for chemical groups (bottom row)
            self._add_chemical_group_atoms(length_cells, passivation, bond_lengths, 
                                           is_top=False, y_offset=y_bottom_offset)
        
        # Top row (last row) - apply y_top_offset if specified
        if len(self.elements) >= width_chains * length_cells:
            start_idx = (width_chains - 1) * length_cells
            for i in range(length_cells):
                self.elements[start_idx + i] = pass_elem
            if y_top_offset is not None:
                # Adjust y positions of top row
                for i in range(length_cells):
                    self.positions[start_idx + i][1] += y_top_offset
            
            # Add additional atoms for chemical groups (top row)
            self._add_chemical_group_atoms(length_cells, passivation, bond_lengths,
                                           is_top=True, y_offset=y_top_offset, 
                                           row_start_idx=start_idx)
        
        # Convert to numpy arrays
        self.positions = np.array(self.positions)
        self.elements = np.array(self.elements)
        self.bonds = np.array(self.bonds)
        
        return self.positions, self.elements, self.bonds
    
    def _add_chemical_group_atoms(self, length_cells, passivation, bond_lengths, 
                                   is_top=False, y_offset=None, row_start_idx=0):
        """Add additional atoms for chemical groups (H in CH, NH, C-OH; O in C=O, C-OH)."""
        
        if passivation == 'NH':
            # Add H atom bonded to N
            # N-H bond length
            b_len = bond_lengths['N-H']
            for i in range(length_cells):
                idx = row_start_idx + i
                pos_N = np.array(self.positions[idx])
                # Add H atom below N (for bottom) or above N (for top)
                if is_top:
                    pos_H = pos_N + np.array([0.0, b_len])
                else:
                    pos_H = pos_N + np.array([0.0, -b_len])
                h_idx = len(self.elements)
                self.elements.append('H')
                self.positions.append(pos_H)
                # Add explicit N-H bond
                self.bonds.append([idx, h_idx])
        
        elif passivation == 'CH':
            # Add H atom bonded to C
            # C-H bond length
            b_len = bond_lengths['C-H']
            for i in range(length_cells):
                idx = row_start_idx + i
                pos_C = np.array(self.positions[idx])
                # Add H atom below C (for bottom) or above C (for top)
                if is_top:
                    pos_H = pos_C + np.array([0.0, b_len])
                else:
                    pos_H = pos_C + np.array([0.0, -b_len])
                h_idx = len(self.elements)
                self.elements.append('H')
                self.positions.append(pos_H)
                # Add explicit C-H bond
                self.bonds.append([idx, h_idx])
        
        elif passivation == 'C=O':
            # Add O atom bonded to C (double bond)
            # C=O bond length
            b_len = bond_lengths['C=O']
            for i in range(length_cells):
                idx = row_start_idx + i
                pos_C = np.array(self.positions[idx])
                # Add O atom below C (for bottom) or above C (for top)
                if is_top:
                    pos_O = pos_C + np.array([0.0, b_len])
                else:
                    pos_O = pos_C + np.array([0.0, -b_len])
                o_idx = len(self.elements)
                self.elements.append('O')
                self.positions.append(pos_O)
                # Add explicit C-O bond
                self.bonds.append([idx, o_idx])
        
        elif passivation == 'C-OH':
            # Add O and H atoms for C-OH group
            # C-O bond length and O-H bond length
            b_CO = bond_lengths['C-O']
            b_OH = bond_lengths['O-H']
            angle = 109.5 * np.pi / 180  # tetrahedral angle in radians
            
            for i in range(length_cells):
                idx = row_start_idx + i
                pos_C = np.array(self.positions[idx])
                
                # Add O atom
                if is_top:
                    pos_O = pos_C + np.array([0.0, b_CO])
                else:
                    pos_O = pos_C + np.array([0.0, -b_CO])
                o_idx = len(self.elements)
                self.elements.append('O')
                self.positions.append(pos_O)
                # Add explicit C-O bond
                self.bonds.append([idx, o_idx])
                
                # Add H atom at tetrahedral angle from C-O bond
                # Position H relative to O to ensure O-H bond length is correct
                if is_top:
                    # C-O points up, O-H at 109.5° from upward
                    # H is to the right and down from O
                    pos_H = pos_O + np.array([b_OH * np.sin(angle), -b_OH * np.cos(angle)])
                else:
                    # C-O points down, O-H at 109.5° from downward
                    # H is to the left and up from O
                    pos_H = pos_O + np.array([-b_OH * np.sin(angle), b_OH * np.cos(angle)])
                h_idx = len(self.elements)
                self.elements.append('H')
                self.positions.append(pos_H)
                # Add explicit O-H bond
                self.bonds.append([o_idx, h_idx])
    
    def _add_zigzag_passivation(self, width_chains, length_cells, passivation):
        """Add passivation atoms at zigzag edges."""
        
        # Bond length for passivation (C-X bond)
        bond_lengths = {
            'N': 1.47,   # C-N
            'NH': 1.47,  # C-N (for NH group)
            'CH': 1.09,  # C-H
            'H': 1.09,   # C-H
            'C=O': 1.43, # C=O
            'C-OH': 1.43 # C-O (for OH)
        }
        
        # Get bond length
        if passivation in bond_lengths:
            b_len = bond_lengths[passivation]
        else:
            b_len = 1.47  # Default to C-N
        
        # Find edge atoms - atoms with only 2 or fewer C neighbors within bond distance
        pos_array = np.array(self.positions)
        elem_array = np.array(self.elements)
        c_mask = elem_array == 'C'
        c_pos = pos_array[c_mask]
        
        # Count neighbors for each C atom
        neighbor_counts = []
        for i, pos in enumerate(c_pos):
            count = 0
            for j, other_pos in enumerate(c_pos):
                if i == j:
                    continue
                dist = np.linalg.norm(pos - other_pos)
                if abs(dist - self.a_CC) < 0.2 * self.a_CC:
                    count += 1
            neighbor_counts.append(count)
        
        # Edge atoms have < 3 neighbors
        edge_indices = [i for i, count in enumerate(neighbor_counts) if count < 3]
        
        # Add passivation to edge atoms
        for idx in edge_indices:
            pos_edge = c_pos[idx]
            
            # Determine passivation direction (away from center of ribbon)
            y_center = np.mean(c_pos[:, 1])
            if pos_edge[1] < y_center:
                # Bottom edge - passivate downward
                pos_pass = pos_edge + np.array([0.0, -b_len])
            else:
                # Top edge - passivate upward
                pos_pass = pos_edge + np.array([0.0, b_len])
            
            if passivation == 'N':
                self.elements.append('N')
            elif passivation == 'NH':
                self.elements.append('N')
                pos_H = pos_pass + np.array([0.0, 1.01 if pos_edge[1] >= y_center else -1.01])
                self.elements.append('H')
                self.positions.append(pos_H)
            elif passivation in ['CH', 'H']:
                self.elements.append('H')
            elif passivation == 'C=O':
                self.elements.append('O')
                pos_O2 = pos_edge + np.array([b_len * 0.5, b_len * 0.866 if pos_edge[1] >= y_center else -b_len * 0.866])
                self.elements.append('O')
                self.positions.append(pos_O2)
            elif passivation == 'C-OH':
                self.elements.append('O')
                pos_H = pos_pass + np.array([0.8, 0.6] if pos_edge[1] >= y_center else [-0.8, -0.6])
                self.elements.append('H')
                self.positions.append(pos_H)
            else:
                self.elements.append('N')
            
            self.positions.append(pos_pass)
    
    def get_structure(self):
        """Return the atomic structure."""
        return self.positions, self.elements
    
    def save_xyz(self, filename):
        """Save structure to XYZ file."""
        with open(filename, 'w') as f:
            f.write(f"{len(self.elements)}\n")
            f.write(f"Graphene zigzag ribbon - a_CC={self.a_CC:.3f} A\n")
            for elem, pos in zip(self.elements, self.positions):
                # Add z-coordinate (0 for 2D)
                z = pos[2] if len(pos) > 2 else 0.0
                f.write(f"{elem} {pos[0]:.6f} {pos[1]:.6f} {z:.6f}\n")
        print(f"Saved to {filename}")
    
    def get_lattice_vectors(self):
        """Return lattice vectors for periodic system."""
        a1 = np.array([self.lattice_constant, 0.0, 0.0])
        a2 = np.array([0.0, 0.0, 20.0])  # Vacuum in y
        a3 = np.array([0.0, 0.0, 20.0])  # Vacuum in z
        return a1, a2, a3


ELEM_MAP = {'H': 1, 'C': 6, 'N': 7, 'O': 8}
ELEM_MAP_INV = {v: k for k, v in ELEM_MAP.items()}

def build_ribbon(passivation, width_chains, length_cells, Lx, a_CC=1.42):
    from doc.Topics.Kekule_Topology.GrapheneRibbonBuilder import GrapheneRibbonBuilder
    xa_nom = a_CC * np.cos(np.pi / 6)
    b = GrapheneRibbonBuilder(a_CC=a_CC)
    scale_x = Lx / (2.0 * xa_nom)
    pos2d, elems, bonds = b.build_zigzag_ribbon(width_chains=width_chains, length_cells=length_cells, passivation=passivation, scale_x=scale_x)
    atypes = np.array([ELEM_MAP[e] for e in elems], dtype=np.int32)
    return np.array(pos2d), atypes, elems

def build_two_ribbon_cell(width_chains=4, length_cells=1, Lx=2.4, a_CC=1.42, L_Hb=2.0, shift_x=0.0):
    pos2d_N,  atypes_N,  elems_N  = build_ribbon('N',  width_chains, length_cells, Lx, a_CC)
    pos2d_NH, atypes_NH, elems_NH = build_ribbon('NH', width_chains, length_cells, Lx, a_CC)
    apos_N  = np.zeros((len(atypes_N),  3));  apos_N[:,  0:2] = pos2d_N
    apos_NH = np.zeros((len(atypes_NH), 3));  apos_NH[:, 0:2] = pos2d_NH
    apos_N[:,  1] -= np.mean(apos_N[:,  1])
    apos_NH[:, 1] -= np.mean(apos_NH[:, 1])
    apos_NH[:, 0] += shift_x * Lx
    y_max_N  = np.max(apos_N[:,  1])
    y_min_NH = np.min(apos_NH[:, 1])
    apos_NH[:, 1] += (y_max_N + L_Hb) - y_min_NH
    y_span_N  = np.max(apos_N[:,  1]) - np.min(apos_N[:,  1])
    y_span_NH = np.max(apos_NH[:, 1]) - np.min(apos_NH[:, 1])
    Ly = y_span_N + y_span_NH + 2 * L_Hb
    apos  = np.vstack([apos_N, apos_NH])
    atypes = np.concatenate([atypes_N, atypes_NH])
    elems  = list(elems_N) + list(elems_NH)
    apos[:, 2]  = 0.0
    apos[:, 1] -= np.mean(apos[:, 1])
    Lz = 20.0
    apos[:, 2] += 0.5 * Lz
    lvs = np.array([[Lx, 0.0, 0.0], [0.0, Ly, 0.0], [0.0, 0.0, Lz]])
    return apos, atypes, elems, lvs



if __name__ == "__main__":
    # Test the builder
    builder = GrapheneRibbonBuilder(a_CC=1.42)
    
    # Build a ribbon with 4 zigzag chains, 8 unit cells, N passivation
    pos, elems = builder.build_zigzag_ribbon(width_chains=4, length_cells=8, passivation='N')
    
    print(f"Built ribbon with {len(elems)} atoms")
    print(f"Elements: {np.unique(elems, return_counts=True)}")
    
    # Save to XYZ
    builder.save_xyz("zigzag_ribbon_N.xyz")
