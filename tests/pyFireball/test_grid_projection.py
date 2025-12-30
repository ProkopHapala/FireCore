import numpy as np
import os
import sys
import matplotlib.pyplot as plt

# Adjust path to FireCore/pyBall directory
sys.path.append("../../")
from pyBall import FireCore as fc
from pyBall.FireballOCL import Grid as ocl_grid

def plot_density_slices(data, title="Density Slices"):
    nx, ny, nz = data.shape
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle(title)
    
    axes[0].imshow(data[nx//2, :, :].T, origin='lower')
    axes[0].set_title("X slice")
    
    axes[1].imshow(data[:, ny//2, :].T, origin='lower')
    axes[1].set_title("Y slice")
    
    axes[2].imshow(data[:, :, nz//2].T, origin='lower')
    axes[2].set_title("Z slice")
    
    plt.show()

def test_pentacene_projection():
    # 1. Load Pentacene
    xyz_path = "../../tests/tUFF/data/xyz/pentacene.xyz"
    if not os.path.exists(xyz_path):
        print(f"Error: XYZ file not found at {xyz_path}")
        return

    with open(xyz_path, 'r') as f:
        lines = f.readlines()
    natoms = int(lines[0])
    atomTypes = []
    atomPos = []
    for line in lines[2:2+natoms]:
        parts = line.split()
        sym = parts[0]
        z = 6 if sym == 'C' else 1
        atomTypes.append(z)
        atomPos.append([float(parts[1]), float(parts[2]), float(parts[3])])
    
    atomTypes = np.array(atomTypes, dtype=np.int32)
    atomPos = np.array(atomPos, dtype=np.float64)

    # 2. Initialize FireCore and run 1 SCF step
    print("Initializing FireCore...")
    fc.initialize(atomType=atomTypes, atomPos=atomPos)
    fc.setVerbosity(1)
    
    print("Running 1 SCF iteration...")
    fc.SCF(atomPos, nmax_scf=1)
    
    # 3. Get Sparse Data
    print("Fetching sparse data...")
    dims = fc.get_HS_dims()
    neighs = fc.get_HS_neighs(dims)
    rho = fc.get_rho_sparse(dims, data=neighs).rho
    
    # 4. Setup Grid Projector
    # Use local test Fdata (contains basis/*.wf1 for H, C)
    fdata_dir = os.path.join(os.path.dirname(__file__), "Fdata", "basis")
    projector = ocl_grid.GridProjector(fdata_dir)
    
    print("Loading basis functions...")
    species_nz = sorted(list(set(atomTypes)))
    projector.load_basis(species_nz)
    
    # 5. Define Grid
    # centered around benzene/pentacene
    pos_min = np.min(atomPos, axis=0) - 4.0
    pos_max = np.max(atomPos, axis=0) + 4.0
    ngrid = [64, 64, 64]
    dCell = (pos_max - pos_min) / ngrid
    grid_spec = {
        'origin': pos_min,
        'dA': [dCell[0], 0, 0],
        'dB': [0, dCell[1], 0],
        'dC': [0, 0, dCell[2]],
        'ngrid': ngrid
    }
    
    # 6. Project Density
    print("Projecting density to grid on GPU...")
    # Map atom data for projector
    atoms_dict = {
        'pos': atomPos,
        'Rcut': np.array([4.5 if z == 6 else 3.5 for z in atomTypes]), # Dummy Rcuts
        'type': atomTypes
    }
    
    dens = projector.project(rho, neighs, atoms_dict, grid_spec)
    
    if dens is not None:
        print(f"Density range: {dens.min()} to {dens.max()}")
        print(f"Total charge (integrated): {np.sum(dens) * np.prod(dCell)}")
        plot_density_slices(dens, title="Pentacene Density (GPU)")
    else:
        print("Projection failed or returned None.")

if __name__ == "__main__":
    test_pentacene_projection()
