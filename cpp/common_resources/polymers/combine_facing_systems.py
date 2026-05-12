#!/usr/bin/env python3
import argparse
import os

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

def parse_lvs(comment):
    parts = comment.split()
    if len(parts) >= 10 and parts[0] == 'lvs':
        # Format: lvs ax ay az bx by bz cx cy cz
        return float(parts[1]), float(parts[5]), float(parts[9])
    return None, None, None

def write_xyz(filepath, atoms, comment):
    with open(filepath, 'w') as f:
        f.write(f"{len(atoms)}\n{comment}\n")
        for el, x, y, z in atoms:
            f.write(f"{el:<4} {x:10.5f} {y:10.5f} {z:10.5f}\n")

def main():
    parser = argparse.ArgumentParser(description="Combines two systems so that the attached molecules face each other.")
    parser.add_argument("-s1", "--system1", required=True, help="XYZ file of the first system.")
    parser.add_argument("-s2", "--system2", required=True, help="XYZ file of the second system.")
    parser.add_argument("-o", "--output", required=True, help="Output XYZ file.")
    parser.add_argument("-dx", "--distance-x", type=float, default=15.0, help="Distance between the axes of both polymers in the X direction.")
    parser.add_argument("-dy", "--shift-y", type=float, default=0.0, help="Mutual shift in the Y direction (allows interleaving of molecules).")
    parser.add_argument("-dz", "--shift-z", type=float, default=0.0, help="Mutual shift in the Z direction.")

    args = parser.parse_args()

    atoms1, comm1 = read_xyz(args.system1)
    atoms2, comm2 = read_xyz(args.system2)

    lx1, ly1, lz1 = parse_lvs(comm1)
    lx2, ly2, lz2 = parse_lvs(comm2)

    if lx1 is None or lx2 is None:
        print("Error: Failed to load 'lvs' lattice from the comment of one of the files.")
        return

    # Rotate the second system by 180° around the center (lx2/2, lz2/2) and shift by dx, dy, dz
    atoms2_transformed = []
    for el, x, y, z in atoms2:
        atoms2_transformed.append((el, (lx2 - x) + args.distance_x, y + args.shift_y, (lz2 - z) + args.shift_z))

    new_lx = lx1 + args.distance_x
    new_comm = f"lvs {new_lx:.5f} 0.0 0.0    0.0 {ly1:.5f} 0.0    0.0 0.0 {lz1:.5f}"

    os.makedirs(os.path.dirname(os.path.abspath(args.output)) or ".", exist_ok=True)
    write_xyz(args.output, atoms1 + atoms2_transformed, new_comm)
    print(f"✅ Combined system created and saved to: {args.output}")
    print(f"📐 New cell (lvs): {new_comm}")

if __name__ == "__main__":
    main()