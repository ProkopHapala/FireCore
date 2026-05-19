#!/usr/bin/env python3
"""
Replicate a unit cell defined by lattice vectors in .xyz comment line.
Format: first comment line contains 'lvs a1 a2 a3 b1 b2 b3 c1 c2 c3'
Usage: python replicate_cell.py input.xyz na nb nc [output.xyz]
"""

import sys
import argparse

def parse_lattice_vectors(comment_line):
    """Parse lattice vectors from comment line starting with 'lvs'."""
    parts = comment_line.strip().split()
    if parts[0].lower() != 'lvs':
        raise ValueError(f"Comment line must start with 'lvs', got: {parts[0]}")
    if len(parts) != 10:
        raise ValueError(f"Expected 9 numbers after 'lvs', got {len(parts)-1}")
    
    nums = list(map(float, parts[1:10]))
    a = nums[0:3]
    b = nums[3:6]
    c = nums[6:9]
    return a, b, c

def read_xyz(filename):
    """Read .xyz file with lattice vectors in comment line."""
    with open(filename, 'r') as f:
        natoms = int(f.readline().strip())
        comment = f.readline().strip()
        a, b, c = parse_lattice_vectors(comment)
        
        atoms = []
        for _ in range(natoms):
            line = f.readline().strip()
            parts = line.split()
            elem = parts[0]
            coords = list(map(float, parts[1:4]))
            extra = parts[4:] if len(parts) > 4 else []
            atoms.append((elem, coords, extra))
    
    return natoms, a, b, c, atoms, comment

def replicate_cell(natoms, a, b, c, atoms, na, nb, nc, original_comment):
    """Replicate the unit cell na x nb x nc times."""
    new_atoms = []
    
    for ia in range(na):
        for ib in range(nb):
            for ic in range(nc):
                offset = [ia*a[i] + ib*b[i] + ic*c[i] for i in range(3)]
                for elem, coords, extra in atoms:
                    new_coords = [coords[i] + offset[i] for i in range(3)]
                    new_atoms.append((elem, new_coords, extra))
    
    # Scale lattice vectors
    new_a = [na * a[i] for i in range(3)]
    new_b = [nb * b[i] for i in range(3)]
    new_c = [nc * c[i] for i in range(3)]
    
    # Build new comment line
    new_comment = f"lvs  {new_a[0]:.3f} {new_a[1]:.3f} {new_a[2]:.3f}   {new_b[0]:.3f} {new_b[1]:.3f} {new_b[2]:.3f}  {new_c[0]:.3f} {new_c[1]:.3f} {new_c[2]:.3f}"
    
    return len(new_atoms), new_a, new_b, new_c, new_atoms, new_comment

def write_xyz(filename, natoms, a, b, c, atoms, comment):
    """Write .xyz file with lattice vectors in comment line."""
    with open(filename, 'w') as f:
        f.write(f"{natoms}\n")
        f.write(f"{comment}\n")
        for elem, coords, extra in atoms:
            line = f"  {elem}  {coords[0]:.5f} {coords[1]:.5f} {coords[2]:.5f}"
            if extra:
                line += " " + " ".join(extra)
            f.write(line + "\n")

def main():
    parser = argparse.ArgumentParser(description='Replicate unit cell from .xyz file with lattice vectors')
    parser.add_argument('input', help='Input .xyz file')
    parser.add_argument('na', type=int, help='Replications along a vector')
    parser.add_argument('nb', type=int, help='Replications along b vector')
    parser.add_argument('nc', type=int, help='Replications along c vector')
    parser.add_argument('output', nargs='?', help='Output .xyz file (default: input_replicated.xyz)')
    
    args = parser.parse_args()
    
    if args.output is None:
        base = args.input.rsplit('.', 1)[0]
        args.output = f"{base}_replicated.xyz"
    
    # Read input
    natoms, a, b, c, atoms, comment = read_xyz(args.input)
    print(f"Read {natoms} atoms from {args.input}")
    print(f"Lattice vectors: a={a}, b={b}, c={c}")
    
    # Replicate
    new_natoms, new_a, new_b, new_c, new_atoms, new_comment = replicate_cell(
        natoms, a, b, c, atoms, args.na, args.nb, args.nc, comment
    )
    print(f"Replicated {args.na}x{args.nb}x{args.nc}: {new_natoms} atoms")
    print(f"New lattice vectors: a={new_a}, b={new_b}, c={new_c}")
    
    # Write output
    write_xyz(args.output, new_natoms, new_a, new_b, new_c, new_atoms, new_comment)
    print(f"Wrote {args.output}")

if __name__ == '__main__':
    main()
