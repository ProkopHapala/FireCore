#!/usr/bin/env python3
"""
ASCII art heterocycle builder.

Supports two ASCII input formats:

1. Single-atom format: every row is an atom layer, characters are element
   symbols, spaces are empty.  Bonds are inferred from the zig-zag lattice
   topology.

       O o O
        c c
        c c
       c c c
       c c c
        c c
       O o O

2. Dimer format: rows alternate between atom rows and bond rows.
   Atom rows contain atom symbols (always converted to carbon in this format).
   Bond rows contain '|' (vertical dimer) and '-' (horizontal dimer).

       O o O
        | |
       | | |
        | |
       O o O

Each character column is 0.5 of the zig-zag lattice constant (dx/2).
Each character row is 0.5 of the armchair lattice constant (dy/2).
"""

import sys
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, '../../..'))

from pyBall.AtomicSystem import AtomicSystem
from pyBall import elements as _elements
from pyBall.KekulePure import KekulePure, make_n_pi
from heterocycle_generator import plot_system, plot_kekule_phases

A_CC = 1.42


# ---------------------------------------------------------------------------
# Bond-length relaxation
# ---------------------------------------------------------------------------
def jacobi_relax_bond_lengths(atoms, L0=1.42, n_iters=1, bmix=0.0):
    """Multi-step Jacobi bond-length relaxation with momentum acceleration.

    Each step is order-independent: a frozen copy of the current positions is
    used to compute all bond corrections, the per-atom displacements are
    accumulated, and only then the positions are updated.

    For each bond (i,j) we compute rij = ri - rj, the current distance d,
    and the vector delta = rij * (L0/d - 1) which would bring the bond exactly
    to length L0 if applied alone.  Half of delta is added to atom i and
    subtracted from atom j.

    Momentum: on the first step the full Jacobi displacement is applied.  From
    the second step onward the displacement is ``v = bmix * v + acc``, where
    ``v`` is the previous step's displacement.
    """
    pos = atoms.apos.copy()
    v = np.zeros_like(pos)
    for step in range(n_iters):
        acc = np.zeros_like(pos)
        for i, j in atoms.bonds:
            rij = pos[i] - pos[j]
            d = np.linalg.norm(rij)
            if d == 0.0:
                continue
            delta = rij * (L0 / d - 1.0)  # vector that takes this bond to length L0
            acc[i] += 0.5 * delta
            acc[j] -= 0.5 * delta
        if step == 0:
            v = acc
        else:
            v = bmix * v + acc
        pos += v
    atoms.apos = pos


def _lines(text):
    return [line.rstrip('\n') for line in text.strip('\n').splitlines()]


def _is_bond_row(line):
    return '|' in line or '-' in line


def _atom_tokens(line):
    """Return list of (col, char) for non-space characters."""
    return [(c, ch) for c, ch in enumerate(line) if ch != ' ']


def _col_to_xidx(r, c):
    """Map ASCII column to lattice x-index (one x-unit = 2 columns)."""
    return c // 2


# ---------------------------------------------------------------------------
# Dimer format -> direct AtomicSystem (x-shift is explicit via spaces)
# ---------------------------------------------------------------------------
def _build_dimer(lines, aCC=A_CC):
    dx = np.sqrt(3.0) * aCC
    pos = []
    enames = []
    enames_original = []  # keep case to distinguish sp2 (upper) vs sp3 (lower)
    rows = []  # per-atom row index
    xidxs = []  # per-atom x index (integer column//2)
    bonds = set()

    row_kind = {}
    row_parity = {}
    for r, line in enumerate(lines):
        tokens = _atom_tokens(line)
        if not tokens:
            continue
        row_kind[r] = 'E' if any(ch.isalpha() for _, ch in tokens) else 'D'
        if row_kind[r] == 'E':
            for c, ch in tokens:
                if ch.isalpha():
                    row_parity[r] = c % 2
                    break

    y_list = {}
    y_pos = 0.0
    for r in range(len(lines)):
        if r == 0:
            y_list[r] = y_pos
            continue
        prev = row_kind.get(r - 1, 'E')
        cur = row_kind.get(r, 'E')
        if prev == 'E' and cur == 'E':
            p0 = row_parity.get(r - 1)
            p1 = row_parity.get(r)
            if (p0 is not None) and (p1 is not None) and (p0 != p1):
                y_pos += aCC / 2.0
            else:
                y_pos += aCC
        else:
            y_pos += (1.5 * aCC) if (prev == 'D' and cur == 'D') else aCC
        y_list[r] = y_pos

    def _add_atom(x, y, e, r, xi):
        i = len(pos)
        pos.append([x, y, 0.0])
        enames.append(e.upper())
        enames_original.append(e)
        rows.append(r)
        xidxs.append(xi)
        return i

    row_atoms = {}
    for r, line in enumerate(lines):
        y = -y_list.get(r, r * aCC)
        tokens = _atom_tokens(line)
        if not tokens:
            continue
        if any(ch.isalpha() for _, ch in tokens):
            for c, ch in tokens:
                xi = c // 2
                i = _add_atom(c * dx / 2.0, y, ch, r, xi)
                row_atoms.setdefault(r, []).append(i)
        else:
            for c, ch in tokens:
                if ch not in ('|', '-', '.'):
                    continue
                xi = c // 2
                x = c * dx / 2.0
                if ch == '.':
                    i = _add_atom(x, y, 'C', r, xi)
                    row_atoms.setdefault(r, []).append(i)
                elif ch == '|':
                    i1 = _add_atom(x, y - aCC / 2.0, 'C', r, xi)
                    i2 = _add_atom(x, y + aCC / 2.0, 'C', r, xi)
                    bonds.add((min(i1, i2), max(i1, i2)))
                    row_atoms.setdefault(r, []).extend([i1, i2])
                elif ch == '-':
                    i1 = _add_atom(x - aCC / 2.0, y, 'C', r, xi)
                    i2 = _add_atom(x + aCC / 2.0, y, 'C', r, xi)
                    bonds.add((min(i1, i2), max(i1, i2)))
                    row_atoms.setdefault(r, []).extend([i1, i2])

    row_x = {}
    for i, (r, xi) in enumerate(zip(rows, xidxs)):
        row_x.setdefault(r, {}).setdefault(xi, []).append(i)

    for r in range(len(lines) - 1):
        rs = row_atoms.get(r, [])
        rt = row_atoms.get(r + 1, [])
        if (not rs) or (not rt):
            continue
        for i in rs:
            xi = xidxs[i]
            p = np.array(pos[i])
            for dx_i in (-1, 0, 1):
                for j in row_x.get(r + 1, {}).get(xi + dx_i, []):
                    if np.linalg.norm(p - np.array(pos[j])) < 1.6 * aCC:
                        bonds.add((min(i, j), max(i, j)))

    atoms = AtomicSystem(apos=np.array(pos), enames=enames)
    atoms.atypes = [_elements.ELEMENT_DICT[e][0] - 1 for e in enames]
    atoms.bonds = sorted(bonds)
    atoms._enames_original = enames_original
    return atoms


# ---------------------------------------------------------------------------
# Single-atom format -> direct AtomicSystem
# ---------------------------------------------------------------------------
def _build_single(lines, aCC=A_CC):
    """Build an AtomicSystem from single-atom ASCII using honeycomb topology."""
    dx = np.sqrt(3.0) * aCC

    atoms = {}  # (r, c) -> element (preserves case)
    atoms_original = {}
    row_parity = {}  # r -> c % 2 of first atom in row
    bond_marks = []  # (r, c, ch) where ch in {'_', '/', '\\'}
    for r, line in enumerate(lines):
        tokens = _atom_tokens(line)
        if not tokens:
            continue
        atom_tokens = [(c, ch) for c, ch in tokens if ch.isalpha()]
        if atom_tokens:
            row_parity[r] = atom_tokens[0][0] % 2
        for c, ch in tokens:
            if ch in ('_', '/', '\\'):
                bond_marks.append((r, c, ch))
                continue
            if not ch.isalpha():
                continue
            atoms[(r, c)] = ch.upper()
            atoms_original[(r, c)] = ch

    # Cumulative y: adjacent rows of the same parity are dimer-bonded (dy=A_CC),
    # adjacent rows of different parity are diagonal-bonded (dy=A_CC/2).
    y = [0.0]
    for r in range(1, len(lines)):
        same_parity = row_parity.get(r) == row_parity.get(r - 1)
        y.append(y[-1] - (aCC if same_parity else aCC / 2.0))

    bonds = []
    for (r, c) in atoms:
        for nr, nc in [(r + 1, c), (r + 1, c - 1), (r + 1, c + 1)]:
            if (nr, nc) in atoms:
                bonds.append(tuple(sorted(((r, c), (nr, nc)))))

    row_atoms_cols = {}
    for (r, c) in atoms:
        row_atoms_cols.setdefault(r, []).append(c)
    for r in row_atoms_cols:
        row_atoms_cols[r].sort()

    def _nearest_left(r, c):
        cols = row_atoms_cols.get(r)
        if not cols:
            return None
        for cc in reversed(cols):
            if cc < c:
                return (r, cc)
        return None

    def _nearest_right(r, c):
        cols = row_atoms_cols.get(r)
        if not cols:
            return None
        for cc in cols:
            if cc > c:
                return (r, cc)
        return None

    def _closest_in_row(r, c):
        cols = row_atoms_cols.get(r)
        if not cols:
            return None
        cc = min(cols, key=lambda x: (abs(x - c), x))
        return (r, cc)

    for r, c, ch in bond_marks:
        if ch == '_':
            a = _nearest_left(r, c)
            b = _nearest_right(r, c)
            if (a is not None) and (b is not None):
                bonds.append(tuple(sorted((a, b))))
        elif ch == '/':
            # shortcut bond through an omitted atom at (r,c): (r-1,c) -- (r+1,c-1)
            a = _closest_in_row(r - 1, c)
            b = _closest_in_row(r + 1, c - 1)
            if (a is not None) and (b is not None):
                bonds.append(tuple(sorted((a, b))))
            else:
                # fallback: local diagonal mark
                a = _nearest_left(r, c)
                b = _nearest_left(r + 1, c)
                if (a is not None) and (b is not None):
                    bonds.append(tuple(sorted((a, b))))
        elif ch == '\\':
            # shortcut bond through an omitted atom at (r,c): (r-1,c) -- (r+1,c+1)
            a = _closest_in_row(r - 1, c)
            b = _closest_in_row(r + 1, c + 1)
            if (a is not None) and (b is not None):
                bonds.append(tuple(sorted((a, b))))
            else:
                # fallback: local diagonal mark
                a = _nearest_right(r, c)
                b = _nearest_right(r + 1, c)
                if (a is not None) and (b is not None):
                    bonds.append(tuple(sorted((a, b))))

    bonds = list(set(bonds))

    idx_map = {}
    pos = []
    enames = []
    for i, ((r, c), el) in enumerate(atoms.items()):
        idx_map[(r, c)] = i
        pos.append([c * dx / 2.0, y[r], 0.0])
        enames.append(el)

    atoms = AtomicSystem(apos=np.array(pos), enames=enames)
    atoms.bonds = [(idx_map[a], idx_map[b]) for a, b in bonds]
    atoms.atypes = [_elements.ELEMENT_DICT[e][0] - 1 for e in enames]
    atoms._enames_original = [atoms_original[k] for k in atoms_original.keys()]
    return atoms


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------
def parse_ascii_art(text):
    """Return an AtomicSystem built from ASCII art."""
    lines = _lines(text)
    if any(_is_bond_row(line) for line in lines):
        return _build_dimer(lines)
    return _build_single(lines)


# ---------------------------------------------------------------------------
# ASCII equivalents of the built-in examples from heterocycle_generator.py.
# In dimer format atom symbols are only visual; the parser converts them to C.
# ---------------------------------------------------------------------------
ASCII_EXAMPLES = {
    'naphthalene': """
  C C
 | | |
  C C
""",
    'naphthalene2': """
  C C
 C C C
 C C C
  C C
""",
    'fulvalene': """
  C C
 | - |
  C C
""",
    'pentacross': """
  C C
 | . |
  C C
""",
    'biphenyl': """
  C 
 | | 
  |
 | | 
  C
""", 
  'biphenyl2': """
  C 
 C C
 C C 
  C
  C
 C C
 C C
  C 
""",
    'phenanthrene': """
  C 
 | | 
  | |
 | | 
  C 
""",
    'phenanthrene2': """
  C
 C C
 C C
  C C
  C C
 C C
 C C
  C
""",
    'perylene': """
  C C
 | | |
  | |
 | | |
  C C
""",
    'perylene2': """
  C C
 C C C
 C C C
  C C
  C C
 C C C
 C C C
  C C  
""",
    'purin': """
  C N
 N C C
 C C /
  N n
""",

    'purin_x': """
  C C
 N C N
  \\C C
  n N
""",
    'purin_y': """
  N n
 C C\\
 N C C
  C N
""",
    '7azaindol': """
  N n
 | | .
  C C
# """,
    'karbazol': """
  C n C
 C C C C
 C C_C C
  C   C
""",
    'biphenylene': """
  C   C
 C C_C C
 C C_C C
  C   C
""",
    'uracil': """
O n O
 C C
 n C
  C
""",
    'cytosin': """
O N n
 C C
 n C
  C
""",
    'guanin': """
n n O
 C C
 N C
  C n
   -
""",
    'NTCDA': """
O o O
 | |
| | |
 | |
O o O
""",
    'TAP': """
  C
 N N 
 C C
| | |
 C C
 N N
  C 
""",
}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description='Generate heterocycle geometry from ASCII art')
    parser.add_argument('--out', '-o', default='/tmp/kekule/heterocycle.svg', help='Output SVG file')
    parser.add_argument('--xyz', default=None, help='Optional output XYZ file')
    parser.add_argument('--example', '-e', default='naphthalene', choices=list(ASCII_EXAMPLES.keys()),  help='Built-in ASCII example to use')
    parser.add_argument('--title', '-t', default=None, help='Plot title')
    parser.add_argument('--size', type=float, default=120, help='Atom marker size')
    parser.add_argument('--aCC', type=float, default=A_CC, help='C-C bond length [Angstrom]')
    parser.add_argument('--plot', type=int, default=1, help='Show matplotlib figure (1) or suppress it (0)')
    parser.add_argument('--dump_bonds', type=int, default=0, help='Print explicit bond list (1) or not (0)')
    parser.add_argument('--kekule', type=int, default=1, help='Run Kekule pi-bond optimization and draw bond orders (1) or not (0)')
    parser.add_argument('--kekule_single', type=int, default=0, help='When --kekule 1, draw a single-panel system plot instead of two-phase raw+snapped')
    parser.add_argument('--relax_bonds', type=int, default=0, help='Number of Jacobi bond-length relaxation steps (0=off)')
    parser.add_argument('--relax_bmix', type=float, default=0.5, help='Momentum mixing factor from 2nd relaxation step onward')
    parser.add_argument('--sym_break', type=int, default=0.0, help='Symmetry-breaking noise added to bond orders before phase-2 localization (0=off)')
    parser.add_argument('--seed', type=int, default=0, help='Random seed used for --sym_break (0 means do not set seed)')
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

    art = ASCII_EXAMPLES[args.example]
    atoms = parse_ascii_art(art)
    atoms.neighs()

    if args.relax_bonds:
        jacobi_relax_bond_lengths(atoms, L0=args.aCC, n_iters=args.relax_bonds, bmix=args.relax_bmix)

    lengths = [np.linalg.norm(atoms.apos[i] - atoms.apos[j]) for i, j in atoms.bonds]
    uniq = sorted({round(d, 3) for d in lengths})
    print(f"{args.example:12s}: atoms={atoms.natoms:2d} bonds={len(atoms.bonds):2d} lengths={uniq}")

    if args.dump_bonds:
        enames_original = getattr(atoms, '_enames_original', None)
        print("ATOMS:")
        for i, p in enumerate(atoms.apos):
            e = atoms.enames[i]
            if enames_original is not None:
                e = enames_original[i]
            print(f"  {i+1:2d} {e:2s}  ({p[0]:7.3f},{p[1]:7.3f},{p[2]:7.3f})")
        print("BONDS (1-indexed):")
        for i, j in atoms.bonds:
            print(f"  {i+1:2d} - {j+1:2d}")
        if atoms.natoms >= 9:
            i, j = 5, 9
            has = (min(i-1, j-1), max(i-1, j-1)) in set((min(a,b), max(a,b)) for a, b in atoms.bonds)
            print(f"CHECK bond {i}-{j}: {has}")

    bond_orders = None
    if args.kekule:
        bo_raw = np.zeros(len(atoms.bonds))
        bo_snap = np.zeros(len(atoms.bonds))
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
        err = None
        try:
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
            if args.sym_break:
                if args.seed:
                    np.random.seed(args.seed)
                k.bo = np.clip(k.bo + args.sym_break * (np.random.rand(k.nbond) - 0.5), 0.0, 1.0)
            if args.solver == 'linsolve':
                k.solve_snap(niter=50)
                F2 = 0.0
            else:
                F2 = k.relax(dt=0.05, nmax=5000, tol=1e-6)
            cls = k.classify()
            print(f"  phase 2 F2={F2:.3e}  pi_BO={np.round(k.pi_bond_orders(),2)}")
            print(f"  single={np.sum(cls==0):2d} aromatic={np.sum(cls==1):2d} double={np.sum(cls==2):2d}")
            err2 = k._A @ k.bo - n_pi
            max_err2 = float(np.max(np.abs(err2))) if err2.size else 0.0
            print(f"  phase 2 constraint max|A@bo-n_pi|={max_err2:.3e}")
            # IMPORTANT: do NOT do naive rounding snap here. Rounding breaks A@bo=n_pi.
            # Phase-2 result is already a localized, constraint-preserving solution.
            bo_snap = k.pi_bond_orders().copy()
        except Exception as e:
            err = e
            print(f"ERROR: {err}")
        title = args.title
        if err is not None:
            title = (title + '\n' if title else '') + f"ERROR: {err}"
        if args.kekule_single:
            # single-panel system plot with snapped Kekule bond orders and n_pi labels
            bo_total = 1.0 + bo_snap
            plot_system(atoms, title=title, fname=args.out, show=args.plot != 0,
                        sz=args.size, n_pi=n_pi, bond_orders=bo_total)
        else:
            plot_kekule_phases(atoms, k, bo_raw=bo_raw, bo_snap=bo_snap,
                               title=title, fname=args.out,
                               show=args.plot != 0, sz=args.size)
    else:
        n_pi = make_n_pi(atoms)
        plot_system(atoms, title=args.title, fname=args.out, show=args.plot != 0,
                    sz=args.size, n_pi=n_pi)
    if args.xyz:
        atoms.saveXYZ(args.xyz)


if __name__ == '__main__':
    main()
