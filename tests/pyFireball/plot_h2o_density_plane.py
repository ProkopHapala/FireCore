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
grid_span       = 4.0  # Grid will go from -grid_span to +grid_span Angstroms
num_grid_points = 400  # Number of points along each axis of the 2D grid
output_plot_file = "h2o_density_plane.png"

# 1. Load Molecule
print(f"Loading molecule from: {xyz_file_path}")
if not os.path.exists(xyz_file_path):
    print(f"Error: {xyz_file_path} not found.")
    sys.exit(1)

mol = AtomicSystem(fname=xyz_file_path)
print("Molecule loaded:")
mol.print()

# Identify Oxygen and Hydrogen atoms
# Assumes O is the first atom, H1 is the second, H2 is the third.
# Adjust if your XYZ file has a different order.
try:
    idx_O  = -1
    idx_H1 = -1
    idx_H2 = -1
    h_indices = []
    for i, ename in enumerate(mol.enames):
        if ename == 'O' and idx_O == -1:
            idx_O = i
        elif ename == 'H':
            h_indices.append(i)
    
    if idx_O == -1: raise ValueError("Oxygen atom not found.")
    if len(h_indices) < 2: raise ValueError("At least two Hydrogen atoms not found.")
    idx_H1 = h_indices[0]
    idx_H2 = h_indices[1]

except (ValueError, IndexError) as e:
    print(f"Error identifying atoms: {e}")
    print("Please ensure H2O.xyz contains O, H, H in a recognizable order.")
    sys.exit(1)

pos_O  = mol.apos[idx_O]
pos_H1 = mol.apos[idx_H1]
pos_H2 = mol.apos[idx_H2]

print(f"Oxygen atom (idx {idx_O}) at {pos_O}")
print(f"Hydrogen 1 (idx {idx_H1}) at {pos_H1}")
print(f"Hydrogen 2 (idx {idx_H2}) at {pos_H2}")

# 2. Initialize FireCore and perform SCF
print("Initializing FireCore...")
fc.setVerbosity(0) # Set verbosity (0 for less output, 1 or more for more)
fc.initialize(atomType=mol.atypes, atomPos=mol.apos)
print("Performing SCF calculation (via evalForce)...")
forces, energies = fc.evalForce(mol.apos, nmax_scf=100)
print(f"SCF converged. Total Energy: {energies[0]:.6f} eV")

# 3. Define the 2D grid plane and points
print("Defining 2D grid in the molecular plane...")
v_O_H1 = pos_H1 - pos_O
v_O_H2 = pos_H2 - pos_O

# Define grid axes based on your specification
# y_axis_dir along (O-H1) + (O-H2)
y_axis_unnormalized = v_O_H1 + v_O_H2
if np.linalg.norm(y_axis_unnormalized) < 1e-6:
    print("Error: O-H1 + O-H2 vector is near zero. Cannot define y-axis.")
    # Fallback: use one of the O-H vectors if the sum is zero (e.g. linear molecule H-O-H)
    # This case shouldn't happen for a typical H2O molecule.
    # If it does, a different strategy for defining the plane might be needed.
    y_axis_unnormalized = v_O_H1
    if np.linalg.norm(y_axis_unnormalized) < 1e-6:
         print("Error: O-H1 vector is also near zero. Cannot define y-axis.")
         sys.exit(1)

y_axis = y_axis_unnormalized / np.linalg.norm(y_axis_unnormalized)

# x_axis_dir perpendicular to y_axis_dir, in the plane of O, H1, H2
# We can get this by taking the cross product of y_axis with the normal to the molecular plane.
# The normal to the molecular plane can be v_O_H1 x v_O_H2
plane_normal_unnormalized = np.cross(v_O_H1, v_O_H2)
if np.linalg.norm(plane_normal_unnormalized) < 1e-6:
    print("Error: Atoms are collinear. Cannot define a unique molecular plane with v_O_H1 x v_O_H2.")
    # Fallback for collinear case: find an arbitrary vector perpendicular to y_axis
    # This is a simplified fallback. For robust handling of linear molecules,
    # one might need to pick an arbitrary perpendicular or allow user input.
    if abs(y_axis[0]) > 1e-6 or abs(y_axis[1]) > 1e-6 :
        x_axis_unnormalized = np.array([-y_axis[1], y_axis[0], 0.0])
    else: # y_axis is along z
        x_axis_unnormalized = np.array([1.0, 0.0, 0.0])
else:
    plane_normal = plane_normal_unnormalized / np.linalg.norm(plane_normal_unnormalized)
    x_axis_unnormalized = np.cross(y_axis, plane_normal)

if np.linalg.norm(x_axis_unnormalized) < 1e-6:
    print("Error: Could not define x-axis perpendicular to y-axis in the plane.")
    sys.exit(1)
x_axis = x_axis_unnormalized / np.linalg.norm(x_axis_unnormalized)

print(f"Grid origin (Oxygen position): {pos_O}")
print(f"Grid x-axis (normalized): {x_axis}")
print(f"Grid y-axis (normalized): {y_axis}")

# Create 1D coordinates for the local grid
local_coords_1d = np.linspace(-grid_span, grid_span, num_grid_points)
X_local, Y_local = np.meshgrid(local_coords_1d, local_coords_1d)

# Transform local grid coordinates to 3D Cartesian coordinates
points_on_grid_3d = np.zeros((num_grid_points * num_grid_points, 3))
density_on_grid_2d = np.zeros((num_grid_points, num_grid_points)) # To store results for imshow

idx = 0
for i in range(num_grid_points):      # Iterate over y_local
    for j in range(num_grid_points):  # Iterate over x_local
        points_on_grid_3d[idx, :] = pos_O + X_local[i, j] * x_axis + Y_local[i, j] * y_axis
        idx += 1

# 4. Calculate Density
print(f"Calculating density at {num_grid_points*num_grid_points} grid points...")
density_values_flat = fc.dens2points(points_on_grid_3d, f_den=1.0, f_den0=0.0)

# Reshape flat density values to 2D grid for plotting
density_on_grid_2d = density_values_flat.reshape((num_grid_points, num_grid_points))
print(f"Calculated density values. Shape: {density_on_grid_2d.shape}")

# 5. Plot
print(f"Plotting 2D density map to {output_plot_file}...")
plt.figure(figsize=(10, 8))
plt.imshow(density_on_grid_2d, origin='lower', cmap='viridis',
           extent=[-grid_span, grid_span, -grid_span, grid_span])
plt.colorbar(label='Electron Density (arbitrary units)')

# Project atom positions onto the 2D grid plane for plotting
def project_to_plane(pos_3d, origin_3d, x_axis_plane, y_axis_plane):
    vec = pos_3d - origin_3d
    x_coord = np.dot(vec, x_axis_plane)
    y_coord = np.dot(vec, y_axis_plane)
    return x_coord, y_coord

atom_colors = {'O': 'red', 'H': 'white'}
for i, ename in enumerate(mol.enames):
    px, py = project_to_plane(mol.apos[i], pos_O, x_axis, y_axis)
    plt.scatter(px, py, s=100, c=atom_colors.get(ename, 'gray'), edgecolors='black', label=f'{ename}{i-idx_O if ename=="O" else i-idx_H1 if i==idx_H1 else i-idx_H2}')

# Create a legend that doesn't duplicate labels
handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys())

plt.xlabel(f'X\' (along defined x-axis, Å)')
plt.ylabel(f'Y\' (along defined y-axis, Å)')
plt.title(f'Electron Density in the H2O molecular plane (centered on O)')
plt.axhline(0, color='gray', linestyle='--', lw=0.5)
plt.axvline(0, color='gray', linestyle='--', lw=0.5)
plt.axis('equal') # Ensure aspect ratio is equal
plt.grid(True, alpha=0.3)
plt.tight_layout()
#plt.savefig(output_plot_file)
print(f"Plot saved to {output_plot_file}")
plt.show() # Uncomment to display the plot interactively

print("Script finished.")

