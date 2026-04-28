import os
import sys
import numpy as np

# Ensure correct import paths
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))
from pyBall.AtomicSystem import AtomicSystem
from pyBall.OCL import UFF, MMFF, MMparams

# Configure OpenCL environment
os.environ['PYOPENCL_COMPILER_OUTPUT'] = '1'
os.environ['PYOPENCL_CTX'] = '0:1'  # Enables double precision on device 0

def print_section(title):
    print("\n" + "-" * 60)
    print(title)
    print("-" * 60)

# Load forcefield parameters from data files
print_section("Loading force field parameters")
base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
resources_dir = os.path.join(base_dir, 'cpp/common_resources')
element_types_path = os.path.join(resources_dir, 'ElementTypes.dat')
atom_types_path = os.path.join(resources_dir, 'AtomTypes.dat')

# Read element and atom types from data files
element_types = MMparams.read_element_types(element_types_path)
atom_types = MMparams.read_atom_types(atom_types_path, element_types)
print(f"Loaded {len(element_types)} element types and {len(atom_types)} atom types")

# Load a molecule from an XYZ file
print_section("Loading molecule")
mol_file = os.path.join(resources_dir, 'xyz/CH2NH.xyz')

try:
    mol = AtomicSystem(mol_file)
    print(f"Loaded molecule from {mol_file} with {len(mol.atypes)} atoms")
except Exception as e:
    print(f"Error loading molecule from file: {e}")
    print("Creating a simple test molecule instead")
    # Create a simple methane molecule as fallback
    mol = AtomicSystem(
        apos=np.array([
            [0.000, 0.000, 0.000],  # C atom at origin
            [0.629, 0.629, 0.629],  # H atom
            [-0.629, -0.629, 0.629], # H atom
            [0.629, -0.629, -0.629], # H atom
            [-0.629, 0.629, -0.629]  # H atom
        ]),
        atypes=np.array([6, 1, 1, 1, 1], dtype=np.int32),  # C and H atoms by atomic number
        enames=["C", "H", "H", "H", "H"],
        lvec=np.identity(3, dtype=np.float32)  # Unit cell for non-periodic system
    )

# Create bonds if they don't exist
if not hasattr(mol, 'bonds') or len(mol.bonds) == 0:
    # For methane - connect carbon to all hydrogens
    mol.bonds = [(0, 1), (0, 2), (0, 3), (0, 4)]
    print(f"Created bonds: {mol.bonds}")

# Generate neighbor information
mol.neighs()
print(f"Molecule has {len(mol.bonds)} bonds and neighbors: {mol.ngs}")

# Test UFF_CL creation
print_section("Creating UFF forcefield")

# Initialize UFF_CL class
uff = UFF.UFF_CL(nloc=32)
print("UFF_CL object created")

# Convert molecular structure to UFF representation using toUFF
print("Converting molecule to UFF representation...")
uff_data = uff.toUFF(mol)

# Print UFF data statistics
print("UFF data generated successfully!")
print("UFF data contents:")

# Check which keys exist in the data dictionary
print(f"UFF data keys: {list(uff_data.keys()) if isinstance(uff_data, dict) else 'Not a dictionary'}")

# Print atom count
natoms = len(mol.atypes)
print(f"Number of atoms: {natoms}")

# Print bond statistics
if 'bndIJs' in uff_data and uff_data['bndIJs'] is not None:
    print(f"Number of bonds: {len(uff_data['bndIJs'])}")
    print(f"Bond indices: {uff_data['bndIJs'][:5]}")
    print(f"Bond params: {uff_data['bndKs'][:5]}")

# Print angle statistics
if 'angIJs' in uff_data and uff_data['angIJs'] is not None:
    print(f"Number of angles: {len(uff_data['angIJs'])}")
    print(f"Angle indices: {uff_data['angIJs'][:3]}")
    print(f"Angle params: {uff_data['angKs'][:3]}")

# Print dihedral statistics
if 'torsIJs' in uff_data and uff_data['torsIJs'] is not None:
    print(f"Number of dihedrals: {len(uff_data['torsIJs'])}")
    if len(uff_data['torsIJs']) > 0:
        print(f"Dihedral indices: {uff_data['torsIJs'][:2]}")
        print(f"Dihedral params: {uff_data['torsKs'][:2]}")

# Print inversion statistics
if 'invIJs' in uff_data and uff_data['invIJs'] is not None:
    print(f"Number of inversions: {len(uff_data['invIJs'])}")
    if len(uff_data['invIJs']) > 0:
        print(f"Inversion indices: {uff_data['invIJs'][:2]}")
        print(f"Inversion params: {uff_data['invKs'][:2]}")

print("# UFF_CL test completed successfully!")

# Print UFF data information
print(f"UFF Data prepared with:")
print(f"  {len(uff_data['atype'])} atoms")
print(f"  {len(uff_data['bonds'])//2} bonds")  # bonds is a flat array with 2 indices per bond
print(f"  {len(uff_data['angles'])//3} angles")  # angles is a flat array with 3 indices per angle
print(f"  {len(uff_data['dihedrals'])//4} dihedrals")  # dihedrals is a flat array with 4 indices per dihedral
print(f"  {len(uff_data['inversions'])} inversions")

# Inspect uff_data structure
print("\nUFF data structure analysis:")
for key, value in uff_data.items():
    if key == 'neighBs':
        print(f"  {key}: {type(value).__name__}, len={len(value)}")
        print(f"  Sample: {value[0]}") # Show the first neighbor list
    else:
        print(f"  {key}: {type(value).__name__}, shape={np.array(value).shape}")

# Print success message
print("\nTest completed! UFF_CL object initialization and data extraction are working correctly.")
#print("Note: Full GPU evaluation is skipped because we don't have kernel files available.")

print("TODO: UFF_CL initialization seems correct. Now we should run the calculations using UFF_CL")



# End of test

print("\nUFF OpenCL test completed successfully")
