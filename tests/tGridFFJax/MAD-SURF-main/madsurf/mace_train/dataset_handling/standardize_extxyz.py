import sys
from ase.io import read, write

if len(sys.argv) != 3:
    print("Usage: python standardize_extxyz.py input.xyz output.xyz")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

# Step 1: Read and write only necessary columns
atoms_list = read(input_file, index=':')
write(output_file, atoms_list, format='extxyz', columns=['symbols', 'positions', 'forces'])

# Step 2: Replace metadata field names in the file
with open(output_file, 'r') as f:
    lines = f.readlines()

with open(output_file, 'w') as f:
    for line in lines:
        line = line.replace(' energy=', ' REF_energy=')
        line = line.replace(':forces:', ':REF_forces:')  # No leading space here to catch 'forces=' inside arrays
        f.write(line)

print(f"Standardized and renamed file written to: {output_file}")

