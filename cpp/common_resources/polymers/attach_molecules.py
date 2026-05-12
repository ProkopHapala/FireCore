#!/usr/bin/env python3
import argparse
import math
import os
import glob

# --- Helper math functions for 3D geometry ---

def norm(v):
    return math.sqrt(sum(x*x for x in v))

def normalize(v):
    l = norm(v)
    return [x/l for x in v] if l > 0 else v

def sub(v1, v2):
    return [x - y for x, y in zip(v1, v2)]

def add(v1, v2):
    return [x + y for x, y in zip(v1, v2)]

def cross(a, b):
    return [a[1]*b[2] - a[2]*b[1],
            a[2]*b[0] - a[0]*b[2],
            a[0]*b[1] - a[1]*b[0]]

def dot(a, b):
    return sum(x*y for x, y in zip(a, b))

def rotate_point(p, axis, theta):
    # Rodrigues rotation formula
    c = math.cos(theta)
    s = math.sin(theta)
    cp = cross(axis, p)
    dp = dot(axis, p)
    return [p[i]*c + cp[i]*s + axis[i]*dp*(1-c) for i in range(3)]

def align_vectors(v_from, v_to):
    # Returns a function that rotates a point to align v_from to v_to
    a = normalize(v_from)
    b = normalize(v_to)
    v = cross(a, b)
    s = norm(v)
    c = dot(a, b)
    if s < 1e-5:
        if c > 0:
            return lambda p: p
        else:
            ortho = cross(a, [1,0,0])
            if norm(ortho) < 1e-5:
                ortho = cross(a, [0,1,0])
            ortho = normalize(ortho)
            return lambda p: rotate_point(p, ortho, math.pi)
    axis = normalize(v)
    theta = math.atan2(s, c)
    return lambda p: rotate_point(p, axis, theta)

def dist(p1, p2):
    return norm(sub(p1, p2))

# --- Functions for reading and writing XYZ ---

def read_xyz(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()
    if not lines: return [], ""
    n_atoms = int(lines[0].strip())
    comment = lines[1].strip()
    atoms = []
    for line in lines[2:2+n_atoms]:
        parts = line.split()
        if len(parts) >= 4:
            atoms.append((parts[0], float(parts[1]), float(parts[2]), float(parts[3])))
    return atoms, comment

def write_xyz(filepath, atoms, comment):
    with open(filepath, 'w') as f:
        f.write(f"{len(atoms)}\n{comment}\n")
        for el, x, y, z in atoms:
            f.write(f"{el:<4} {x:10.5f} {y:10.5f} {z:10.5f}\n")

def process_anchor(atoms):
    si_idx = -1
    for i, (el, x, y, z) in enumerate(atoms):
        if el == 'Si':
            si_idx = i
            break
    if si_idx == -1:
        return None

    si_pos = atoms[si_idx][1:]
    neighbors = []
    for i, (el, x, y, z) in enumerate(atoms):
        if i == si_idx: continue
        d = dist(si_pos, (x, y, z))
        if el == 'H' and d < 1.3:
            neighbors.append((i, el, (x, y, z), d))
        elif el != 'H' and d < 2.0:
            neighbors.append((i, el, (x, y, z), d))
            
    heavy_neighbors = [n for n in neighbors if n[1] != 'H']
    h_neighbors = [n for n in neighbors if n[1] == 'H']
    
    # Carbon should have exactly 4 bonds (tetravalence).
    # One bond will be reserved for connecting the two structures.
    needed_h = 4 - 1 - len(heavy_neighbors)
    
    to_remove = []
    bond_dir = None
    
    # If we have more hydrogens than needed, delete the excess
    if len(h_neighbors) > needed_h:
        h_to_remove = len(h_neighbors) - needed_h
        for i in range(h_to_remove):
            to_remove.append(h_neighbors[i][0])
            if i == 0:
                # The direction of the new bond will be where the first deleted hydrogen was
                bond_dir = normalize(sub(h_neighbors[i][2], si_pos))
                
    if bond_dir is None:
        # If there was nothing to delete, find empty space as the 
        # opposite direction to the sum of existing bond vectors
        v_sum = [0.0, 0.0, 0.0]
        for n in neighbors:
            v = normalize(sub(n[2], si_pos))
            v_sum = add(v_sum, v)
        if norm(v_sum) > 1e-5:
            bond_dir = normalize([-x for x in v_sum])
        else:
            bond_dir = [1.0, 0.0, 0.0]
            
    return {
        'si_idx': si_idx,
        'anchor_pos': si_pos,
        'bond_dir': bond_dir,
        'to_remove': to_remove
    }

# --- Main processing ---

def main():
    parser = argparse.ArgumentParser(description="Attaches molecules from a directory to a polymer (Silicon to Silicon).")
    parser.add_argument("-p", "--polymer", required=True, help="XYZ file of the polymer (containing Si anchor).")
    parser.add_argument("-m", "--molecules-dir", required=True, help="Directory with molecules (containing Si anchors).")
    parser.add_argument("-o", "--output-dir", required=True, help="Target directory for the resulting merged molecules.")
    parser.add_argument("-b", "--bond-length", type=float, default=1.54, help="Target length of the new bond (C-C is approx 1.54A).")
    
    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    poly_atoms, poly_comment = read_xyz(args.polymer)
    p_info = process_anchor(poly_atoms)
    if not p_info:
        print("Error: Polymer does not contain a silicon (Si) anchor atom!")
        return

    # Rename Si to C and remove excess hydrogens
    poly_atoms_clean = []
    for i, a in enumerate(poly_atoms):
        if i in p_info['to_remove']: continue
        if i == p_info['si_idx']:
            poly_atoms_clean.append(('C', *a[1:]))
        else:
            poly_atoms_clean.append(a)

    mol_files = glob.glob(os.path.join(args.molecules_dir, '*.xyz'))
    if not mol_files:
        print(f"Warning: No .xyz files found in the directory '{args.molecules_dir}'.")
        return

    for mol_file in mol_files:
        mol_name = os.path.basename(mol_file)
        mol_atoms, _ = read_xyz(mol_file)
        m_info = process_anchor(mol_atoms)
        if not m_info:
            print(f"Skipping {mol_name} - no Si anchor atom found.")
            continue

        m_atoms_clean = []
        for i, a in enumerate(mol_atoms):
            if i in m_info['to_remove']: continue
            if i == m_info['si_idx']:
                m_atoms_clean.append(('C', *a[1:]))
            else:
                m_atoms_clean.append(a)

        rotate = align_vectors(m_info['bond_dir'], [-x for x in p_info['bond_dir']])
        target_pos = add(p_info['anchor_pos'], [x * args.bond_length for x in p_info['bond_dir']])
        
        mol_final = [ (a[0], *add(rotate(sub(a[1:], m_info['anchor_pos'])), target_pos)) for a in m_atoms_clean ]
        
        out_file = os.path.join(args.output_dir, f"attached_{mol_name}")
        write_xyz(out_file, poly_atoms_clean + mol_final, poly_comment)
        print(f"✅ Created: {out_file} (Aligned {mol_name} to polymer)")

if __name__ == "__main__":
    main()