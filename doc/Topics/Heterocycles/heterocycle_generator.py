#!/usr/bin/env python3
"""
Heterocycle generator from a sparse rectangular-grid description of a hexagonal lattice.

Rows are either E (edge) layers or D (dimer) layers:
    E layers: first and last row, single atom per unit cell
    D layers: all internal rows, one or two atoms per unit cell

Coordinate system:
    x : zig-zag direction
    y : armchair direction

D-layer mode markers (inside a unit spec):
    default / '|' : two atoms form a vertical dimer along y (armchair)
    '.'           : one atom at the nominal D-row position
    '-'           : two atoms form a horizontal dimer along x (zig-zag),
                    creates pentagons/heptagons

Input format: a list of row specifications. Each row is one of:
    'CCNCn'                              E-layer, starts at x=0, atoms from string
    (start_x, 'CCNCn')                   E-layer, explicit start index
    (start_x, length)                    D-layer, default vertical dimer, all carbon
    (start_x, length, '1N,2.,3-OC')      D-layer sparse compact spec
    (start_x, length, [(i,'N'), ...])    D-layer sparse tuple spec
    (start_x, length, {i:'N', ...})      D-layer sparse dict spec

D-layer spec rules:
    'N' or '|N'  -> vertical dimer, both atoms N
    'NC' or '|NC'-> vertical dimer, atoms N and C
    '.'          -> single C atom
    '.N'         -> single N atom
    '-'          -> horizontal C dimer
    '-OC'        -> horizontal dimer, atoms O and C

Legacy formats are still accepted for backward compatibility:
    (start_x, length, 'S'), (start_x, length, 'R'), (start_x, length, ['D','S','R','D'])

Example:
    system = [
        (1, 2),   # row 0 E layer
        (1, 3),   # row 1 D layer, 3 vertical C dimers
        (1, 2),   # row 2 E layer
    ]

Outputs:
    - SVG plot of the generated geometry (via pyBall.plotSystem)
    - Optional XYZ file (via AtomicSystem.saveXYZ)

Depends on matplotlib and pyBall (AtomicSystem, plotUtils, elements).
"""

import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt

# Use pyBall's AtomicSystem and plotSystem
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..'))
from pyBall.AtomicSystem import AtomicSystem
from pyBall import plotUtils as _plot_utils
from pyBall import elements as _elements
from pyBall.KekulePure import KekulePure, make_n_pi

# ---------------------------------------------------------------------------
# Default geometry parameters (graphene-like units)
# ---------------------------------------------------------------------------
A_CC = 1.42
DX = np.sqrt(3.0) * A_CC          # horizontal spacing between E columns
DY = 2.0 * A_CC                    # vertical spacing between E layers (E-row to E-row)
DIMER_OFFSET = A_CC / 2.0          # half of the D-layer dimer bond length along y (armchair)

# ---------------------------------------------------------------------------
# Parser
# ---------------------------------------------------------------------------
def parse_row(row_spec, row_idx, n_rows):
    """Return (start_x, atoms) for one row.

    Layer type is determined by the spec.  Explicit string D-layer specs
    (e.g. '|', '-', '0NC,1.') are always parsed as D layers.  Otherwise, the
    default convention is: first and last rows are E (edge, single atom per
    unit cell), all internal rows are D (dimer, two atoms per unit cell).
    Each atom tuple is (x_index, element_name, subtype).
    """
    is_edge = (row_idx == 0) or (row_idx == n_rows - 1)
    if isinstance(row_spec, str):
        start_x = 0
        atom_str = row_spec
    elif isinstance(row_spec, (list, tuple)):
        if len(row_spec) < 2:
            raise ValueError(f"Row spec must have at least 2 elements: {row_spec}")
        start_x = int(row_spec[0])
        second = row_spec[1]
        if isinstance(second, str):
            # E-layer explicit string
            atom_str = second
            if len(row_spec) > 2:
                raise ValueError(f"E-layer string row accepts only (start, str): {row_spec}")
        elif isinstance(second, int) or isinstance(second, float):
            # Sparse length-based row
            length = int(second)
            has_d_spec = (len(row_spec) >= 3 and isinstance(row_spec[2], str))
            if is_edge and not has_d_spec:
                # E layer: length single atoms, optional substitutions
                substitutions = _parse_substitutions(row_spec[2]) if len(row_spec) >= 3 else {}
                return _build_e_row(start_x, length, substitutions)
            else:
                # D layer: length units, optional per-unit specs
                unit_specs = _parse_d_row_args(row_spec, length)
                return _build_d_row(start_x, length, unit_specs)
        else:
            raise ValueError(f"Second row-spec item must be int (length) or str (atom string): {second}")
    else:
        raise ValueError(f"Unsupported row spec type: {type(row_spec)}")

    # Build E-layer atoms from string (preserves case for sp2/sp3 distinction)
    atoms = []
    for i, ch in enumerate(atom_str):
        atoms.append((start_x + i, ch, 'E'))
    return start_x, atoms


def _parse_substitutions(subs):
    if subs is None:
        return {}
    result = {}
    for item in subs:
        if isinstance(item, (list, tuple)) and len(item) == 2:
            pos = int(item[0])
            spec = item[1]
            if isinstance(spec, str):
                spec = spec  # preserve case for sp2/sp3 distinction
            result[pos] = spec
        else:
            raise ValueError(f"Substitution item must be (position, spec): {item}")
    return result


def _parse_d_spec(spec):
    """Parse a D-layer unit spec string into (mode, element_spec).

    Mode markers (returned directly as the mode string):
        '.' -> single atom
        '-' -> horizontal/rotated dimer
        '|' -> vertical dimer (explicit)
        otherwise -> vertical dimer (implicit)

    Returns:
        mode (str): '|', '.', or '-'
        element_spec (str): element characters for the atom(s)
    """
    spec = str(spec).strip()
    if not spec:
        return '|', ''
    marker = spec[0]
    if marker == '.':
        return '.', spec[1:]
    if marker == '-':
        return '-', spec[1:]
    if marker == '|':
        return '|', spec[1:]
    return '|', spec


def _parse_d_modifications(third, length):
    """Parse a D-layer modification argument into a dict position -> spec.

    Accepts:
        'D' / 'S' / 'R'             legacy uniform mode (mapped to '|', '.', '-')
        '1NC,2.,3-OC'               compact sparse string
        [(1,'N'), (2,'.')]          sparse tuple list
        {1:'N', 2:'.'}              sparse dict
        ['D','S','R','D']           legacy full mode list
    """
    if third is None:
        return {}
    legacy_map = {'D': '', 'S': '.', 'R': '-'}
    if isinstance(third, str):
        third = third.upper().strip()
        if third in legacy_map:
            # Legacy uniform mode
            return {i: legacy_map[third] for i in range(length)}
        if third and third[0].isdigit():
            # Compact sparse string (e.g. '1.', '1-', '1NC,2.,3-OC')
            if ',' in third:
                result = {}
                for item in third.split(','):
                    item = item.strip()
                    if not item:
                        continue
                    pos_str = ''
                    for j, ch in enumerate(item):
                        if ch.isdigit():
                            pos_str += ch
                        else:
                            break
                    if not pos_str:
                        raise ValueError(f"Compact spec item must start with position: {item}")
                    pos = int(pos_str)
                    spec = item[len(pos_str):]
                    result[pos] = spec
                return result
            # Single compact spec without commas: '1.', '2-OC'
            pos_str = ''
            for j, ch in enumerate(third):
                if ch.isdigit():
                    pos_str += ch
                else:
                    break
            spec = third[len(pos_str):]
            if pos_str and spec:
                return {int(pos_str): spec}
        # Uniform spec for all positions (e.g. 'NC', '.N', '-OC')
        return {i: third for i in range(length)}
    if isinstance(third, dict):
        return {int(pos): str(spec) for pos, spec in third.items()}
    if isinstance(third, (list, tuple)):
        if len(third) == length and all(isinstance(x, str) for x in third):
            # Legacy full mode list
            return {i: legacy_map[x.upper()] for i, x in enumerate(third)}
        # Sparse tuple list
        return _parse_substitutions(third)
    raise ValueError(f"Unsupported D-layer modification argument: {third}")


def _parse_d_row_args(row_spec, length):
    """Parse optional D-layer arguments into a list of unit specs.

    Returns list of length `length` with one spec string per unit.
    """
    unit_specs = [''] * length
    if len(row_spec) >= 3:
        mods = _parse_d_modifications(row_spec[2], length)
        for pos, spec in mods.items():
            if 0 <= pos < length:
                unit_specs[pos] = spec
            else:
                raise ValueError(f"D-layer unit position {pos} out of range [0, {length-1}]")
    if len(row_spec) >= 4:
        mods = _parse_d_modifications(row_spec[3], length)
        for pos, spec in mods.items():
            if 0 <= pos < length:
                unit_specs[pos] = spec
            else:
                raise ValueError(f"D-layer unit position {pos} out of range [0, {length-1}]")
    if len(row_spec) >= 5:
        mods = _parse_d_modifications(row_spec[4], length)
        for pos, spec in mods.items():
            if 0 <= pos < length:
                unit_specs[pos] = spec
            else:
                raise ValueError(f"D-layer unit position {pos} out of range [0, {length-1}]")
    return unit_specs


def _build_e_row(start_x, length, substitutions):
    """Build E-layer atom list. Each unit is one atom."""
    atoms = []
    for i in range(length):
        x_idx = start_x + i
        ename = substitutions.get(i, 'C')
        if not isinstance(ename, str):
            ename = str(ename)
        atoms.append((x_idx, ename, 'E'))
    return start_x, atoms


def _build_d_row(start_x, length, unit_specs):
    """Build D-layer atom list. Each unit is one or two atoms."""
    atoms = []
    for i in range(length):
        x_idx = start_x + i
        mode, elem_spec = _parse_d_spec(unit_specs[i] if i < len(unit_specs) else '')
        if mode == '.':
            ename = elem_spec[0] if elem_spec else 'C'
            atoms.append((x_idx, ename, '.'))
        elif mode in ('|', '-'):
            if len(elem_spec) == 0:
                e1 = e2 = 'C'
            elif len(elem_spec) == 1:
                e1 = e2 = elem_spec[0]
            elif len(elem_spec) == 2:
                e1, e2 = elem_spec[0], elem_spec[1]
            else:
                raise ValueError(f"D-layer dimer spec must be 0, 1 or 2 chars: {elem_spec}")
            atoms.append((x_idx, e1, f'{mode}1'))
            atoms.append((x_idx, e2, f'{mode}2'))
        else:
            raise ValueError(f"Unknown D-layer mode: {mode}")
    return start_x, atoms


# ---------------------------------------------------------------------------
# Geometry generation
# ---------------------------------------------------------------------------
def build_geometry(system, dx=DX, dy=DY, dimer_offset=DIMER_OFFSET):
    """Generate (positions, elements, subtypes, original_case) from the system description."""
    pos = []
    enames = []
    enames_original = []
    subtypes = []
    row_info = []  # (row_idx, unit_idx, subtype)
    n_rows = len(system)

    # Pre-parse all rows and determine layer types from atom subtypes.
    parsed = []
    for row_idx, row_spec in enumerate(system):
        start_x, atoms = parse_row(row_spec, row_idx, n_rows)
        layer = 'E' if all(st == 'E' for _, _, st in atoms) else 'D'
        parsed.append((row_idx, start_x, atoms, layer))

    # Cumulative y: D-D transition needs 3*dy/4 (3*a_CC/2) to form hexagons;
    # E-D and D-E use dy/2 (a_CC) as normal.
    y_pos = 0.0
    y_list = []
    for k, (row_idx, _, _, layer) in enumerate(parsed):
        if k > 0:
            prev_layer = parsed[k - 1][3]
            y_pos += 3 * dy / 4.0 if (prev_layer == 'D' and layer == 'D') else dy / 2.0
        y_list.append(y_pos)

    for row_idx, start_x, atoms, layer in parsed:
        y = y_list[row_idx]
        shift = 0.0 if row_idx % 2 == 0 else -0.5 * dx
        for atom in atoms:
            x_idx, ename, st = atom
            cx = x_idx * dx + shift
            if st == 'E':
                pos.append([cx, y, 0.0])
                enames.append(ename.upper())
                enames_original.append(ename)
                subtypes.append('E')
                row_info.append((row_idx, x_idx, 'E'))
            elif st == '.':
                pos.append([cx, y, 0.0])
                enames.append(ename.upper())
                enames_original.append(ename)
                subtypes.append('.')
                row_info.append((row_idx, x_idx, '.'))
            elif st in ('|1', '|2'):
                offset = -dimer_offset if st == '|1' else dimer_offset
                pos.append([cx, y + offset, 0.0])
                enames.append(ename.upper())
                enames_original.append(ename)
                subtypes.append(st)
                row_info.append((row_idx, x_idx, st))
            elif st in ('-1', '-2'):
                offset = -dimer_offset if st == '-1' else dimer_offset
                pos.append([cx + offset, y, 0.0])
                enames.append(ename.upper())
                enames_original.append(ename)
                subtypes.append(st)
                row_info.append((row_idx, x_idx, st))
            else:
                raise ValueError(f"Unknown atom subtype: {st}")

    return np.array(pos), enames, enames_original, subtypes, row_info


def build_atomic_system(system, dx=DX, dy=DY, dimer_offset=DIMER_OFFSET):
    """Build an AtomicSystem from the sparse grid description.

    Returns an AtomicSystem with apos, enames, bonds and aux_labels set.
    Bonds are derived from the grid topology for each D-layer mode.
    """
    pos, enames, enames_original, subtypes, row_info = build_geometry(system, dx=dx, dy=dy, dimer_offset=dimer_offset)
    atoms = AtomicSystem(apos=np.array(pos), enames=enames)
    atoms.atypes = [_elements.ELEMENT_DICT[e][0] - 1 for e in enames]
    atoms.aux_labels = [f"{r}:{x}:{st}" for r, x, st in row_info]
    atoms._enames_original = enames_original

    # Build bonds explicitly from the lattice topology.
    # Each D-layer atom connects to the nearest atoms in the adjacent rows at the
    # relevant x positions. This handles both E-E (original) and D-D adjacency.
    bonds = []
    idx = {(r, x, st): i for i, (r, x, st) in enumerate(row_info)}

    # Collect atoms by row and x_idx (multiple atoms per x_idx are possible in D rows)
    row_x = {}
    for i, (r, x, st) in enumerate(row_info):
        row_x.setdefault(r, {}).setdefault(x, []).append((i, st))

    def _nearest(targets, pi):
        if not targets:
            return None
        return min(targets, key=lambda j: np.linalg.norm(pos[j] - pi))

    def _candidates(r, x_target):
        return [j for j, st in row_x.get(r, {}).get(x_target, [])]

    def _adj_x(r_cur, r_adj, x):
        # Even rows have shift=0, odd rows have shift=-dx/2.
        # Even->odd: nearest atoms are at x and x+1 in the odd row.
        # Odd->even (or same parity): nearest atoms are at x-1 and x.
        if r_cur % 2 == 0 and r_adj % 2 == 1:
            return (x, x + 1)
        return (x - 1, x)

    for i, (r, x, st) in enumerate(row_info):
        if st == 'E':
            continue
        if st == '.':
            for dr in (-1, 1):
                for ax in _adj_x(r, r + dr, x):
                    cands = _candidates(r + dr, ax)
                    j = _nearest(cands, pos[i])
                    if j is not None and np.linalg.norm(pos[i] - pos[j]) < 1.5 * A_CC:
                        bonds.append((min(i, j), max(i, j)))
        elif st in ('|1', '|2'):
            dr = -1 if st == '|1' else 1
            for ax in _adj_x(r, r + dr, x):
                cands = _candidates(r + dr, ax)
                j = _nearest(cands, pos[i])
                if j is not None and np.linalg.norm(pos[i] - pos[j]) < 1.5 * A_CC:
                    bonds.append((min(i, j), max(i, j)))
            partner = '|2' if st == '|1' else '|1'
            j = idx.get((r, x, partner))
            if j is not None and i < j:
                bonds.append((i, j))
        elif st in ('-1', '-2'):
            for dr in (-1, 1):
                for ax in _adj_x(r, r + dr, x):
                    cands = _candidates(r + dr, ax)
                    j = _nearest(cands, pos[i])
                    if j is not None and np.linalg.norm(pos[i] - pos[j]) < 1.5 * A_CC:
                        bonds.append((min(i, j), max(i, j)))
            partner = '-2' if st == '-1' else '-1'
            j = idx.get((r, x, partner))
            if j is not None and i < j:
                bonds.append((i, j))

    atoms.bonds = np.array(sorted(set(bonds)), dtype=np.int32)
    return atoms


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def plot_system(atoms, title=None, fname='heterocycle.svg',
                figsize=(8, 8), dpi=150, show=False, sz=50., bond_orders=None,
                n_pi=None):
    """Plot the AtomicSystem using plotSystem from plotUtils and save as SVG.

    If bond_orders is provided, it must be a 1-D array of total bond orders
    (one per bond in atoms.bonds).  Double bonds are drawn thicker and
    aromatic bonds are drawn green.
    """
    fig, ax = plt.subplots(figsize=figsize)
    _plot_utils.plotSystem(atoms, axes=(0, 1), bBonds=True, bLabels=True, sz=sz)

    if bond_orders is not None and atoms.bonds is not None:
        bo = np.asarray(bond_orders)
        if len(bo) != len(atoms.bonds):
            raise ValueError(f"bond_orders length {len(bo)} != number of bonds {len(atoms.bonds)}")
        lws = np.ones(len(bo), dtype=float)
        colors = np.array(['k'] * len(bo), dtype=object)
        # aromatic ~ 1.5, double ~ 2.0, single ~ 1.0
        aromatic = np.abs(bo - 1.5) < 0.15
        double = bo > 1.7
        lws[aromatic] = 2.5
        colors[aromatic] = 'green'
        lws[double] = 3.5
        _plot_utils.plotBonds(links=atoms.bonds, ps=atoms.apos, lws=lws,
                              colors=colors, axes=(0, 1))

    if n_pi is not None:
        n_pi = np.asarray(n_pi)
        for i, (x, y) in enumerate(atoms.apos[:, [0, 1]]):
            ax.text(x, y, f" {int(n_pi[i])}", fontsize=8, ha='left', va='top',
                    color='red', zorder=6,
                    bbox=dict(boxstyle='round,pad=0.1', facecolor='yellow',
                              edgecolor='none', alpha=0.6))

    if title:
        ax.set_title(title)
    ax.text(0.02, 0.98, f"N = {atoms.natoms}", transform=ax.transAxes,
            fontsize=12, verticalalignment='top', color='black')
    ax.set_aspect('equal')
    ax.axis('off')
    fig.tight_layout()
    fig.savefig(fname, format='svg', dpi=dpi, bbox_inches='tight')
    if show:
        plt.show()
    plt.close(fig)
    print(f"Saved SVG: {fname}")


def _bond_midpoints(atoms):
    """Return (x, y) midpoints of each bond."""
    bonds = np.asarray(atoms.bonds)
    p0 = atoms.apos[bonds[:, 0]]
    p1 = atoms.apos[bonds[:, 1]]
    return 0.5 * (p0 + p1)


def _plot_bond_phase(ax, atoms, pi_bos, sz, title, n_pi=None, bLabels=True):
    """Plot one phase: atoms + bonds colored by bond order + pi labels."""
    plt.sca(ax)
    _plot_utils.plotSystem(atoms, axes=(0, 1), bBonds=True, bLabels=True, sz=sz)

    bo = np.asarray(pi_bos)
    total = 1.0 + bo
    lws = np.ones(len(bo), dtype=float)
    colors = np.array(['k'] * len(bo), dtype=object)
    aromatic = np.abs(total - 1.5) < 0.15
    double = total > 1.7
    lws[aromatic] = 2.5
    colors[aromatic] = 'green'
    lws[double] = 3.5
    plt.sca(ax)
    _plot_utils.plotBonds(links=atoms.bonds, ps=atoms.apos, lws=lws,
                          colors=colors, axes=(0, 1))

    if n_pi is not None:
        n_pi = np.asarray(n_pi)
        for i, (x, y) in enumerate(atoms.apos[:, [0, 1]]):
            ax.text(x, y, f" {int(n_pi[i])}", fontsize=8, ha='left', va='top',
                    color='red', zorder=6,
                    bbox=dict(boxstyle='round,pad=0.1', facecolor='yellow',
                              edgecolor='none', alpha=0.6))

    if bLabels:
        mids = _bond_midpoints(atoms)[:, [0, 1]]
        for (x, y), b in zip(mids, bo):
            ax.text(x, y, f"{b:.2f}", fontsize=8, ha='center', va='center',
                    color='darkblue', zorder=5,
                    bbox=dict(boxstyle='round,pad=0.15', facecolor='white',
                              edgecolor='none', alpha=0.7))

    if title:
        ax.set_title(title)
    ax.set_aspect('equal')
    ax.axis('off')


def plot_kekule_phases(atoms, k, bo_raw=None, bo_snap=None, title=None,
                       fname='kekule_phases.svg', figsize=(14, 7), dpi=150,
                       show=False, sz=50.):
    """Plot raw and snapped Kekule bond orders side by side with labels."""
    if bo_raw is None:
        bo_raw = k.pi_bond_orders()
    if bo_snap is None:
        bo_snap = k.snap().copy()
    n_pi = getattr(k, 'n_pi', None)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
    _plot_bond_phase(ax1, atoms, bo_raw, sz, title='phase 1: raw',
                     n_pi=n_pi, bLabels=True)
    _plot_bond_phase(ax2, atoms, bo_snap, sz, title='phase 2: snapped',
                     n_pi=n_pi, bLabels=True)

    if title:
        fig.suptitle(title, fontsize=14)
    fig.tight_layout()
    fig.savefig(fname, format='svg', dpi=dpi, bbox_inches='tight')
    if show:
        plt.show()
    plt.close(fig)
    print(f"Saved SVG: {fname}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
EXAMPLES = {
    'naphthalene': [
        (1, 2),   # E:   . .                           
        (1, 3),   # D:  | | |
        (1, 2),   # E:   . .
    ],
    'fulvalene': [
        (1, 2),         # E:    . .
        (1, 3, '1-'),   # D:   | - |
        (1, 2),         # E:    . .
    ],
    'pentacross': [
        (1, 2),         # E:   . .
        (1, 3, '1.'),   # D:  | . |   
        (1, 2),         # E:   . .
    ],
    'purin': [
        (1, 'CN'),         # E:    . N
        (1, 3, '0NC,2.'),  # D:   | | .    : 1NC
        (1, 'Nn'),         # E:    N n
    ],
    'uracil': [
        (0, 'OnO'),      # E:    O n O  
        (1, 2, '0Cn'),   # D:     | |       : 0Cn
        (1,1),           # E:      .
    ],
    'cytosin': [
        (0, 'ONn'),     # E:    O N n  
        (1, 2, '0Cn'),  # D:     | |        : 0Cn
        (1,1),          # E:      .
    ],
    'guanin': [
        (0,'nnO'),      # E:    n  n O    
        (1, 2, '0CN'),  # D:     | |        : 0CN
        (1,2,'0Cn,1.N'),        # D:      | .       : 0Cn,1.N
        (2,1),          # E:       .       
    ],
    'NTCDA': [
        (0,'OoO'),     # E:  O o O
        (1, 2),        # D:   | |
        (0, 3),        # D:  | | |
        (1, 2),        # D:   | | 
        (0,'OoO'),     # E:  O o O
    ],
}


def main():
    parser = argparse.ArgumentParser(description='Generate heterocycle geometry SVG from sparse grid description')
    parser.add_argument('--out', '-o', default='/tmp/kekule/heterocycle.svg', help='Output SVG file')
    parser.add_argument('--xyz', default=None, help='Optional output XYZ file')
    parser.add_argument('--input', '-i', default=None, help='Python file containing a "system" variable')
    parser.add_argument('--example', '-e', default='naphthalene', choices=list(EXAMPLES.keys()), help='Built-in example system to use when --input is not given')
    parser.add_argument('--title', '-t', default=None, help='Plot title')
    parser.add_argument('--size', type=float, default=120, help='Atom marker size (scales van der Waals radii)')
    parser.add_argument('--aCC', type=float, default=A_CC, help='C-C bond length [Angstrom]')
    parser.add_argument('--plot', type=int, default=1, help='Show matplotlib figure (1) or suppress it (0)')
    parser.add_argument('--kekule', type=int, default=0, help='Run Kekule pi-bond optimization and draw bond orders (1) or not (0)')
    parser.add_argument('--Kval', type=float, default=50.0, help='Atom-sum stiffness (K_atom_sum)')
    parser.add_argument('--Kloc', type=float, default=5.0, help='Snap/localization stiffness (K_snap)')
    parser.add_argument('--Karo', type=float, default=0.5, help='Aromatization stiffness (K_arom)')
    parser.add_argument('--aromatic', type=int, default=1, help='Allow aromatic (0.5) bond orders (1) or force integer (0/1) localization (0)')
    parser.add_argument('--solver', default='linsolve', choices=['linsolve', 'gd'], help='Solver: linsolve (recommended, quadratic) or gradient descent (fallback)')
    parser.add_argument('--kkt', type=int, default=0, help='Save KKT matrix plot (1) or not (0)')
    parser.add_argument('--kkt_mode', default='signed', choices=['signed', 'logabs'], help='KKT plot mode')
    parser.add_argument('--kkt_cmap', default=None, help='KKT colormap (signed: seismic; logabs: magma/inferno)')
    args = parser.parse_args()

    out_dir = os.path.dirname(os.path.abspath(args.out))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    if args.input:
        import importlib.util
        spec = importlib.util.spec_from_file_location("user_system", args.input)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        system = mod.system
    else:
        system = EXAMPLES[args.example]

    dx = np.sqrt(3.0) * args.aCC
    dy = 2.0 * args.aCC
    dimer = args.aCC / 2.0

    atoms = build_atomic_system(system, dx=dx, dy=dy, dimer_offset=dimer)
    atoms.neighs()

    if args.kekule:
        n_pi = make_n_pi(atoms)
        allow_aromatic = (args.aromatic != 0)
        k = KekulePure(atoms, n_pi=n_pi, Kval=args.Kval, Kloc=0.0, Karo=args.Karo,
                       Kbound=1.0, allow_aromatic=allow_aromatic)
        if args.kkt:
            if args.kkt_cmap is None:
                kkt_cmap = 'seismic' if args.kkt_mode == 'signed' else 'magma'
            else:
                kkt_cmap = args.kkt_cmap
            root, ext = os.path.splitext(os.path.abspath(args.out))
            kkt_out = root + f"_kkt_{args.kkt_mode}.png"
            k.plot_kkt_matrix(mode=args.kkt_mode, cmap=kkt_cmap, fname=kkt_out, show=False)
            print(f"Saved KKT matrix: {kkt_out}")
        if args.solver == 'linsolve':
            k.solve_quadratic(Kloc=0.0)
            F2 = 0.0
        else:
            F2 = k.relax(dt=0.05, nmax=5000, tol=1e-6)
        bo_raw = k.pi_bond_orders()
        print(f"  phase 1 F2={F2:.3e}  pi_BO={np.round(bo_raw,2)}")
        print(f"  n_pi={np.asarray(n_pi, dtype=int)}")
        # stage 2: localization
        k.Kloc = args.Kloc
        if args.solver == 'linsolve':
            k.solve_snap(niter=50)
            F2 = 0.0
        else:
            F2 = k.relax(dt=0.05, nmax=5000, tol=1e-6)
        cls = k.classify()
        print(f"  phase 2 F2={F2:.3e}  pi_BO={np.round(k.pi_bond_orders(),2)}")
        print(f"  single={np.sum(cls==0):2d} aromatic={np.sum(cls==1):2d} double={np.sum(cls==2):2d}")
        bo_snap = k.snap().copy()
        plot_kekule_phases(atoms, k, bo_raw=bo_raw, bo_snap=bo_snap,
                           title=args.title, fname=args.out,
                           show=args.plot != 0, sz=args.size)
    else:
        plot_system(atoms, title=args.title, fname=args.out, show=args.plot != 0,
                    sz=args.size)
    if args.xyz:
        atoms.saveXYZ(args.xyz)


if __name__ == '__main__':
    main()
