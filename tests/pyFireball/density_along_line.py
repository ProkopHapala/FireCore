import sys
import os
import numpy as np
import matplotlib.pyplot as plt

# Adjust the path to your FireCore/pyBall directory if necessary
sys.path.append("../../") # Assuming this script is in a subdirectory of FireCore/tests/pyFireball or similar
from pyBall.AtomicSystem import AtomicSystem
from pyBall import FireCore as fc

# --- Configuration ---
xyz_file_path = "data/xyz/H2O.xyz"  # Path to your H2O.xyz file
num_points_on_line = 100
output_plot_file = "h2o_density_profile.png"

print("#=========== DEBUG 1 ")
mol = AtomicSystem(fname=xyz_file_path)
print("Molecule loaded:")
mol.print()

def find_line_atoms(mol, e1='O', e2='H'):
    # Identify Oxygen and one Hydrogen atom (assuming standard H2O order: O, H1, H2)
    # You might need to adjust indices based on your specific XYZ file
    idx1 = -1
    idx2 = -1
    for i, ename in enumerate(mol.enames):
        if ename == e1 and idx1 == -1:
            idx1 = i
        elif ename == e2 and idx2 == -1:
            idx2 = i
    if idx1 == -1 or idx2 == -1:
        print("Error: Could not find atom %s and atom %s in the XYZ file." % (e1, e2))
        sys.exit(1)
    print(f"Atom {e1} index: {idx1}, Atom {e2} atom index: {idx2}")
    pos1 = mol.apos[idx1]
    pos2 = mol.apos[idx2]
    return pos1, pos2

e1 = 'O'
e2 = 'H'
print("#=========== DEBUG 1 ")
pos1, pos2 = find_line_atoms(mol)

print("#=========== DEBUG 1 ")

# 2. Initialize FireCore and perform SCF
print("Initializing FireCore...")
fc.setVerbosity(1) # Set verbosity (0 for less output, 1 or more for more)
print( "mol.atypes ", mol.atypes )
print( "mol.apos   ", mol.apos )
fc.initialize(atomType=mol.atypes, atomPos=mol.apos)
print("Performing SCF calculation (via evalForce)...")
# evalForce performs an SCF cycle and calculates energy/forces
# We are interested in the density matrix (rho) that is populated after SCF
forces, energies = fc.evalForce(mol.apos, nmax_scf=100) # Use default 100 SCF iterations
print(f"SCF converged. Total Energy: {energies[0]:.6f} eV")

# 3. Define Points along the O-H bond
print(f"Defining {num_points_on_line} points along the O-H bond...")
points_on_line = np.zeros((num_points_on_line, 3))
distances = np.linspace(0, 1, num_points_on_line) # Normalized distance parameter
print("#=========== DEBUG 2 ")

for i, t in enumerate(distances):
    points_on_line[i, :] = (1 - t) * pos1 + t * pos2

# The `dens2points` function expects points in shape (N, 3)
# The Fortran backend `firecore_dens2points` expects (3, N),
# but your Python wrapper `dens2points` should handle the len(points) and pass N.
# If your wrapper expects (3,N) directly, you'd do points_on_line.T

# 4. Calculate Density
print("Calculating density at specified points using dens2points...")
# We want rho_SCF, so f_den=1.0 and f_den0=0.0
# The `points` argument for your simplified `dens2points` should be (N,3)
density_values = fc.dens2points(points_on_line, f_den=1.0, f_den0=0.0)
print(f"Calculated {len(density_values)} density values.")
print("#=========== DEBUG 3 ")

# 5. Plot
print(f"Plotting density profile to {output_plot_file}...")
# Calculate actual distances from Oxygen for the x-axis
actual_distances = np.linalg.norm(points_on_line - pos1, axis=1)

plt.figure(figsize=(10, 6))
plt.plot(actual_distances, density_values, 'b.-', label='SCF Density (rho_SCF)')
plt.xlabel(f'along {e2}-{e1} bond (Å)')
plt.ylabel('Electron Density (arbitrary units)')
plt.title(f'Electron Density Profile along {e2}-{e1} bond in {xyz_file_path}')
plt.axvline(0, color='gray', linestyle='--', label=f'{e1}')
plt.axvline(np.linalg.norm(pos2 - pos1), color='lightcoral', linestyle='--', label=f'{e2}')
plt.legend()
plt.grid(True)
plt.tight_layout()
#plt.savefig(output_plot_file)
print(f"Plot saved to {output_plot_file}")
plt.show() # Uncomment to display the plot interactively

print("Script finished.")
