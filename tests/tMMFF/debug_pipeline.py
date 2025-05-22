#!/usr/bin/env python3
import numpy as np
import os
import sys
import matplotlib.pyplot as plt

# Add parent directory to path so we can import pyBall
sys.path.append('/home/indranil/git/FireCore')
from pyBall.OCL.GridFF import GridFF_cl

def save_grid(grid, name, data_path="./data/"):
    """Save grid data and create visualization"""
    os.makedirs(data_path, exist_ok=True)
    np.save(f"{data_path}/{name}.npy", grid)
    
    # Create visualization of middle slices
    plt.figure(figsize=(15, 5))
    mid_z = grid.shape[0]//2
    mid_y = grid.shape[1]//2
    mid_x = grid.shape[2]//2
    
    plt.subplot(131)
    plt.imshow(grid[mid_z,:,:], origin='lower')
    plt.colorbar()
    plt.title(f'XY plane (z={mid_z})')
    
    plt.subplot(132)
    plt.imshow(grid[:,mid_y,:], origin='lower')
    plt.colorbar()
    plt.title(f'XZ plane (y={mid_y})')
    
    plt.subplot(133)
    plt.imshow(grid[:,:,mid_x], origin='lower')
    plt.colorbar()
    plt.title(f'YZ plane (x={mid_x})')
    
    plt.suptitle(f'{name}')
    plt.savefig(f"{data_path}/{name}.png")
    plt.close()

def find_max_position(grid):
    """Find position of maximum absolute value in grid"""
    abs_grid = np.abs(grid)
    max_idx = np.unravel_index(np.argmax(abs_grid), abs_grid.shape)
    return max_idx, grid[max_idx]

def debug_pipeline():
    """Debug the GPU pipeline step by step"""
    print("=== GPU Pipeline Debug ===")
    
    # Initialize GridFF_cl
    grid = GridFF_cl()
    
    # Load test data - adjust paths as needed
    data_path = "./data/NaCl_1x1_L1/"
    
    # Step 1: Load CPU Qgrid for reference
    cpu_qgrid = np.load(f"{data_path}/Qgrid.npy")
    print(f"CPU Qgrid shape: {cpu_qgrid.shape}")
    cpu_max_pos, cpu_max_val = find_max_position(cpu_qgrid)
    print(f"CPU Qgrid max at {cpu_max_pos}, value: {cpu_max_val:.6f}")
    
    # Step 2: Run GPU pipeline with instrumentation
    # Set up grid with same parameters as the test
    nx, ny, nz = 40, 40, 400  # Adjust based on your data
    grid.set_grid((nx, ny, nz), spacing=(0.1, 0.1, 0.1), origin=(-2.0, -2.0, -20.0))
    
    # Load atoms from file (adjust as needed)
    atoms = []
    try:
        with open(f"{data_path}/atoms.xyz", 'r') as f:
            lines = f.readlines()
            for line in lines[2:]:  # Skip header
                parts = line.split()
                if len(parts) >= 4:
                    x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                    q = 1.0  # Default charge, adjust if needed
                    atoms.append([x, y, z, q])
    except:
        print("Could not load atoms.xyz, using default test atoms")
        # Create simple test atoms if file not available
        atoms = [
            [0.0, 0.0, 0.0, 1.0],
            [0.0, 0.0, 2.0, -1.0]
        ]
    
    atoms_np = np.array(atoms, dtype=np.float32)
    print(f"Using {len(atoms_np)} atoms for testing")
    
    # Checkpoint 1: Project atoms to grid
    print("\n=== Checkpoint 1: Atom Projection ===")
    grid.makeCoulombEwald_slab(atoms_np, bOld=False, bSaveQgrid=True)
    
    # Load the GPU Qgrid that was saved
    gpu_qgrid = np.load(f"{data_path}/Qgrid_gpu.npy")
    print(f"GPU Qgrid shape: {gpu_qgrid.shape}")
    gpu_max_pos, gpu_max_val = find_max_position(gpu_qgrid)
    print(f"GPU Qgrid max at {gpu_max_pos}, value: {gpu_max_val:.6f}")
    
    # Checkpoint 2: After FFT (Vgrid)
    print("\n=== Checkpoint 2: After Poisson/FFT ===")
    # Run poisson separately to capture intermediate state
    vgrid = grid.poisson(bReturn=True)
    save_grid(vgrid, "Vgrid_gpu_debug", data_path)
    vgrid_real = vgrid[..., 0]  # Extract real part
    vgrid_max_pos, vgrid_max_val = find_max_position(vgrid_real)
    print(f"Vgrid shape: {vgrid.shape}")
    print(f"Vgrid (real) max at {vgrid_max_pos}, value: {vgrid_max_val:.6f}")
    
    # Checkpoint 3: After Laplace
    print("\n=== Checkpoint 3: After Laplace ===")
    vin = grid.laplace_real_loop_inert(niter=16, bReturn=True)
    save_grid(vin, "Vin_gpu_debug", data_path)
    vin_max_pos, vin_max_val = find_max_position(vin)
    print(f"Vin shape: {vin.shape}")
    print(f"Vin max at {vin_max_pos}, value: {vin_max_val:.6f}")
    
    # Checkpoint 4: After slabPotential
    print("\n=== Checkpoint 4: After slabPotential ===")
    # Create a buffer from Vin
    import pyopencl as cl
    mf = cl.mem_flags
    vin_buff = cl.Buffer(grid.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=vin)
    vcoul = grid.slabPotential(vin_buff, nz_slab=nz, dipol=0.0, bDownload=True)
    save_grid(vcoul, "VCoul_gpu_debug", data_path)
    vcoul_max_pos, vcoul_max_val = find_max_position(vcoul)
    print(f"VCoul shape: {vcoul.shape}")
    print(f"VCoul max at {vcoul_max_pos}, value: {vcoul_max_val:.6f}")
    
    # Summary
    print("\n=== Pipeline Summary ===")
    print(f"CPU Qgrid max at {cpu_max_pos}, value: {cpu_max_val:.6f}")
    print(f"GPU Qgrid max at {gpu_max_pos}, value: {gpu_max_val:.6f}")
    print(f"Vgrid max at {vgrid_max_pos}, value: {vgrid_max_val:.6f}")
    print(f"Vin max at {vin_max_pos}, value: {vin_max_val:.6f}")
    print(f"VCoul max at {vcoul_max_pos}, value: {vcoul_max_val:.6f}")
    
    # Calculate shifts between stages
    qgrid_to_vgrid = tuple(v - q for v, q in zip(vgrid_max_pos, gpu_max_pos))
    vgrid_to_vin = tuple(i - v for i, v in zip(vin_max_pos, vgrid_max_pos))
    vin_to_vcoul = tuple(c - i for c, i in zip(vcoul_max_pos, vin_max_pos))
    
    print("\n=== Shifts Between Stages ===")
    print(f"Qgrid → Vgrid shift: {qgrid_to_vgrid}")
    print(f"Vgrid → Vin shift: {vgrid_to_vin}")
    print(f"Vin → VCoul shift: {vin_to_vcoul}")
    print(f"Total shift: {tuple(c - q for c, q in zip(vcoul_max_pos, gpu_max_pos))}")

if __name__ == "__main__":
    debug_pipeline()
