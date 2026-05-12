#!/usr/bin/env python3
import argparse
import os
import math

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

def dist_pbc(p1, p2, ly):
    dx = p1[0] - p2[0]
    dy = p1[1] - p2[1]
    dz = p1[2] - p2[2]
    # Minimum Image Convention along Y axis
    if ly > 0:
        dy = dy - ly * round(dy / ly)
    return math.sqrt(dx*dx + dy*dy + dz*dz)

def fix_hydrogens(atoms, ly, n_units=1, monomer_len=0, backbone=None):
    # Target valencies for standard elements
    target_valency = {'C': 4, 'Si': 4, 'N': 3, 'O': 2, 'S': 2, 'P': 3}
    to_remove = set()
    
    check_indices = set(range(len(atoms)))
    if backbone:
        check_indices = set()
        h_idx, t_idx = backbone
        for i in range(n_units):
            check_indices.add(h_idx + i * monomer_len)
            check_indices.add(t_idx + i * monomer_len)
    
    for i in check_indices:
        el, x, y, z = atoms[i]
        if el == 'H' or el not in target_valency:
            continue
        
        # Find neighbors considering Y-axis periodicity
        heavy_neighbors = []
        h_neighbors = []
        for j, (el2, x2, y2, z2) in enumerate(atoms):
            if i == j: continue
            d = dist_pbc((x, y, z), (x2, y2, z2), ly)
            if el2 == 'H' and d < 1.3:
                h_neighbors.append((j, d))
            elif el2 != 'H' and d < 2.0:
                heavy_neighbors.append((j, d))
        
        # Check for over-valency (e.g. from capping hydrogens at monomer boundaries)
        current_bonds = len(heavy_neighbors) + len(h_neighbors)
        excess_h = current_bonds - target_valency[el]
        
        if excess_h > 0 and len(h_neighbors) > 0:
            # Score hydrogens by their proximity to OTHER heavy atoms.
            # The ones closest to other heavy atoms are the overlapping caps.
            h_clash_scores = []
            for h_idx, d_hi in h_neighbors:
                hx, hy, hz = atoms[h_idx][1:]
                min_dist_to_other = float('inf')
                for k, (el3, x3, y3, z3) in enumerate(atoms):
                    if k == i or k == h_idx or el3 == 'H': continue
                    d_other = dist_pbc((hx, hy, hz), (x3, y3, z3), ly)
                    if d_other < min_dist_to_other:
                        min_dist_to_other = d_other
                h_clash_scores.append((min_dist_to_other, h_idx))
            
            # Sort ascending (smallest distance -> most clashing -> remove first)
            h_clash_scores.sort()
            for k in range(min(excess_h, len(h_clash_scores))):
                to_remove.add(h_clash_scores[k][1])
    
    return [a for idx, a in enumerate(atoms) if idx not in to_remove]

def main():
    parser = argparse.ArgumentParser(description="Builds a 1D periodic polymer from a given monomer XYZ file.")
    parser.add_argument("-m", "--monomer", required=True, help="Path to the monomer XYZ file.")
    parser.add_argument("-n", "--n-units", type=int, default=10, help="Number of monomer units in the polymer.")
    parser.add_argument("-dy", "--shift-y", type=float, required=True, help="Translation along Y-axis per monomer unit (in Angstroms).")
    parser.add_argument("-a", "--anchor", type=int, required=True, help="1-based index of the anchor atom in the monomer to be replaced by Si (e.g., 20 for the 20th atom).")
    parser.add_argument("-b", "--backbone", type=int, nargs=2, help="1-based indices of the head and tail atoms of the monomer's backbone.")
    parser.add_argument("--lx", type=float, default=15.0, help="Cell size in X direction.")
    parser.add_argument("--lz", type=float, default=15.0, help="Cell size in Z direction.")
    parser.add_argument("-o", "--output", required=True, help="Output XYZ file path.")

    args = parser.parse_args()

    monomer_atoms, _ = read_xyz(args.monomer)
    if not monomer_atoms:
        print(f"Error: Failed to load monomer from {args.monomer}.")
        return

    anchor_idx = args.anchor - 1 # Convert 1-based to 0-based
    if anchor_idx < 0 or anchor_idx >= len(monomer_atoms):
        print(f"Error: Anchor index (--anchor {args.anchor}) is out of range (1 to {len(monomer_atoms)}).")
        return

    poly_atoms = []
    
    # Identify the middle monomer unit to place the anchor (Si)
    middle_unit = args.n_units // 2

    for i in range(args.n_units):
        sy = i * args.shift_y
        for j, (el, x, y, z) in enumerate(monomer_atoms):
            # Replace the designated atom with Si in the middle monomer
            if i == middle_unit and j == anchor_idx:
                current_el = 'Si'
            else:
                current_el = el
            
            # Append atom with Y-axis translation, centering X and Z
            poly_atoms.append((current_el, x + args.lx/2.0, y + sy, z + args.lz/2.0))

    # Calculate exact periodic cell dimension in Y axis
    ly = args.n_units * args.shift_y
    
    # Setup backbone indices if provided (convert to 0-based)
    backbone_idx = None
    if args.backbone:
        h_idx = args.backbone[0] - 1
        t_idx = args.backbone[1] - 1
        backbone_idx = (h_idx, t_idx)
        
        if args.n_units > 1:
            tail_pos = poly_atoms[t_idx][1:]
            head_next_pos = poly_atoms[h_idx + len(monomer_atoms)][1:]
            dist = math.sqrt(sum((a - b)**2 for a, b in zip(tail_pos, head_next_pos)))
            print(f"ℹ️  Vzdálenost páteřních atomů ({args.backbone[0]} a {args.backbone[1]}) mezi monomery: {dist:.3f} Å (uprav -dy, pokud je to daleko od ~1.5 Å)")

    # Remove overlapping capping hydrogens at monomer junctions
    poly_atoms = fix_hydrogens(poly_atoms, ly, args.n_units, len(monomer_atoms), backbone_idx)

    lvs_line = f"lvs {args.lx:.5f} 0.0 0.0    0.0 {ly:.5f} 0.0    0.0 0.0 {args.lz:.5f}"

    os.makedirs(os.path.dirname(os.path.abspath(args.output)) or ".", exist_ok=True)
    write_xyz(args.output, poly_atoms, lvs_line)
    
    print(f"✅ Polymer created and saved to: {args.output}")

if __name__ == "__main__":
    main()