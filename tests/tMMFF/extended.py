import numpy as np
import sys

def read_xyz(filename):
    """Reads an XYZ file and extracts the number of atoms, lattice vectors, and atomic positions."""
    with open(filename, "r") as f:
        lines = f.readlines()

    num_atoms = int(lines[0].strip())
    lattice_vector_line = lines[1].strip().split()
    
    # Extract lattice vectors
    lvs = np.array([
        [float(lattice_vector_line[1]), float(lattice_vector_line[2]), float(lattice_vector_line[3])],
        [float(lattice_vector_line[4]), float(lattice_vector_line[5]), float(lattice_vector_line[6])],
        [float(lattice_vector_line[7]), float(lattice_vector_line[8]), float(lattice_vector_line[9])]
    ])
    
    # Extract atomic positions
    atoms = []
    for line in lines[2:num_atoms+2]:
        parts = line.split()
        element = parts[0]
        x, y, z = map(float, parts[1:4])
        charge = float(parts[4]) if len(parts) > 4 else 0.0  # Handle missing charge
        atoms.append((element, x, y, z, charge))

    return num_atoms, lvs, atoms

def generate_supercell(input_file, Nx, Ny, Nz):
    """Generates a supercell by replicating the original structure Nx × Ny × Nz times."""
    num_atoms, lattice_vectors, atoms = read_xyz(input_file)

    # Compute new lattice vectors
    new_lattice_vectors = lattice_vectors.copy()
    new_lattice_vectors[0] *= Nx  # Scale x-axis
    new_lattice_vectors[1] *= Ny  # Scale y-axis
    new_lattice_vectors[2] *= Nz  # Scale z-axis

    # Generate new atomic positions
    new_atoms = []
    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz):
                dx, dy, dz = (i * lattice_vectors[0][0], j * lattice_vectors[1][1], k * lattice_vectors[2][2])
                for atom in atoms:
                    element, x, y, z, charge = atom
                    new_atoms.append((element, x + dx, y + dy, z + dz, charge))

    return new_lattice_vectors, new_atoms

def write_xyz(output_file, new_lattice_vectors, new_atoms):
    """Writes the new supercell structure to an XYZ file."""
    with open(output_file, "w") as f:
        f.write(f"{len(new_atoms)}\n")
        f.write(f"lvs {new_lattice_vectors[0][0]} {new_lattice_vectors[0][1]} {new_lattice_vectors[0][2]} "
                f"{new_lattice_vectors[1][0]} {new_lattice_vectors[1][1]} {new_lattice_vectors[1][2]} "
                f"{new_lattice_vectors[2][0]} {new_lattice_vectors[2][1]} {new_lattice_vectors[2][2]}\n")
        for atom in new_atoms:
            f.write(f"{atom[0]} {atom[1]:.6f} {atom[2]:.6f} {atom[3]:.6f} {atom[4]:.1f}\n")

def main():
    if len(sys.argv) != 5:
        print("Usage: python extend_xyz.py <input_file> <Nx> <Ny> <Nz>")
        sys.exit(1)

    input_file = sys.argv[1]
    Nx, Ny, Nz = map(int, sys.argv[2:5])

    output_file = f"{input_file.split('.')[0]}_{Nx}x{Ny}x{Nz}.xyz"

    new_lattice_vectors, new_atoms = generate_supercell(input_file, Nx, Ny, Nz)
    write_xyz(output_file, new_lattice_vectors, new_atoms)

    print(f"Supercell written to {output_file}")

if __name__ == "__main__":
    main()

