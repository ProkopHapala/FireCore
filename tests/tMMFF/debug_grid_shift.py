#!/usr/bin/env python3

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import time
import pyopencl as cl

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from pyBall.OCL.GridFF import GridFF_cl
import pyBall.OCL.clUtils as clu

print("=== Grid Shift Debugging Script ===")

def init_gpu_grid(name):
    """Initialize a GPU grid with the given dataset name"""
    
    data_path = f"./data/{name}/"
    print(f"Loading data from {data_path}")
    
    grid_ff = GridFF_cl()
    
    print("Loading grid parameters...")
    grid_ff.loadGrid(name, data_path=data_path)
    
    print("Grid parameters:")
    print(f"  origin (g0): {grid_ff.gcl.g0}")
    print(f"  spacing (dg): {grid_ff.gcl.dg}")
    print(f"  dimensions (ns): {grid_ff.gsh.ns}")
    
    return grid_ff

def load_grid_data(filepath):
    """Load grid data from a numpy file"""
    if os.path.exists(filepath):
        data = np.load(filepath)
        print(f"Loaded {filepath}, shape: {data.shape}")
        return data
    else:
        print(f"Error: File {filepath} not found")
        return None

def compare_grid_positions(cpu_grid, gpu_grid, name="comparison"):
    """Compare positions of features in CPU and GPU grids"""
    
    # Ensure grids have the same shape
    if cpu_grid.shape != gpu_grid.shape:
        print(f"Error: CPU grid shape {cpu_grid.shape} != GPU grid shape {gpu_grid.shape}")
        return
        
    print(f"\n=== {name} Grid Position Analysis ===")
    
    # Basic stats
    print("CPU grid range:", cpu_grid.min(), cpu_grid.max())
    print("GPU grid range:", gpu_grid.min(), gpu_grid.max())
    
    # Find positions of significant values (non-zero elements)
    cpu_nonzero = np.where(np.abs(cpu_grid) > 1e-5)
    gpu_nonzero = np.where(np.abs(gpu_grid) > 1e-5)
    
    print(f"CPU grid has {len(cpu_nonzero[0])} significant points")
    print(f"GPU grid has {len(gpu_nonzero[0])} significant points")
    
    # Find positions of maximum absolute values
    cpu_max_idx = np.unravel_index(np.argmax(np.abs(cpu_grid)), cpu_grid.shape)
    gpu_max_idx = np.unravel_index(np.argmax(np.abs(gpu_grid)), gpu_grid.shape)
    
    print(f"CPU grid max at (z,y,x): {cpu_max_idx}, value: {cpu_grid[cpu_max_idx]:.6f}")
    print(f"GPU grid max at (z,y,x): {gpu_max_idx}, value: {gpu_grid[gpu_max_idx]:.6f}")
    
    # Check for shift
    z_diff = gpu_max_idx[0] - cpu_max_idx[0]
    y_diff = gpu_max_idx[1] - cpu_max_idx[1]
    x_diff = gpu_max_idx[2] - cpu_max_idx[2]
    
    print(f"Position difference (GPU - CPU): z={z_diff}, y={y_diff}, x={x_diff}")
    
    if z_diff != 0 or y_diff != 0 or x_diff != 0:
        print(f"SHIFT DETECTED in {name}!")
    else:
        print(f"No shift detected in {name}.")
        
    # Plot comparison between CPU and GPU grids
    plot_grid_comparison(cpu_grid, gpu_grid, z_diff, y_diff, x_diff, name)
    
    return z_diff, y_diff, x_diff

def plot_grid_comparison(cpu_grid, gpu_grid, z_shift=0, y_shift=0, x_shift=0, name="comparison"):
    """Create plots comparing CPU and GPU grids"""
    
    # Get grid dimensions
    nz, ny, nx = cpu_grid.shape
    
    # Get middle slice indices
    mid_z = nz // 2
    mid_y = ny // 2
    mid_x = nx // 2
    
    # Create figure
    fig, axs = plt.subplots(3, 3, figsize=(15, 12))
    fig.suptitle(f"CPU vs GPU Grid Comparison: {name}\nShift detected: z={z_shift}, y={y_shift}, x={x_shift}", fontsize=16)
    
    # Define color normalization for better visualization
    vmin = min(cpu_grid.min(), gpu_grid.min())
    vmax = max(cpu_grid.max(), gpu_grid.max())
    norm = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
    
    # Plot XY plane (z fixed)
    axs[0,0].set_title(f"CPU: XY plane (z={mid_z})")
    im1 = axs[0,0].imshow(cpu_grid[mid_z,:,:], origin='lower', cmap='viridis', norm=norm)
    divider = make_axes_locatable(axs[0,0])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im1, cax=cax)
    
    axs[0,1].set_title(f"GPU: XY plane (z={mid_z})")
    im2 = axs[0,1].imshow(gpu_grid[mid_z,:,:], origin='lower', cmap='viridis', norm=norm)
    divider = make_axes_locatable(axs[0,1])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im2, cax=cax)
    
    axs[0,2].set_title(f"Difference: XY plane")
    im3 = axs[0,2].imshow(gpu_grid[mid_z,:,:] - cpu_grid[mid_z,:,:], origin='lower', cmap='coolwarm')
    divider = make_axes_locatable(axs[0,2])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im3, cax=cax)
    
    # Plot XZ plane (y fixed)
    axs[1,0].set_title(f"CPU: XZ plane (y={mid_y})")
    im4 = axs[1,0].imshow(cpu_grid[:,mid_y,:], origin='lower', cmap='viridis', norm=norm)
    divider = make_axes_locatable(axs[1,0])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im4, cax=cax)
    
    axs[1,1].set_title(f"GPU: XZ plane (y={mid_y})")
    im5 = axs[1,1].imshow(gpu_grid[:,mid_y,:], origin='lower', cmap='viridis', norm=norm)
    divider = make_axes_locatable(axs[1,1])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im5, cax=cax)
    
    axs[1,2].set_title(f"Difference: XZ plane")
    im6 = axs[1,2].imshow(gpu_grid[:,mid_y,:] - cpu_grid[:,mid_y,:], origin='lower', cmap='coolwarm')
    divider = make_axes_locatable(axs[1,2])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im6, cax=cax)
    
    # Plot YZ plane (x fixed)
    axs[2,0].set_title(f"CPU: YZ plane (x={mid_x})")
    im7 = axs[2,0].imshow(cpu_grid[:,:,mid_x], origin='lower', cmap='viridis', norm=norm)
    divider = make_axes_locatable(axs[2,0])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im7, cax=cax)
    
    axs[2,1].set_title(f"GPU: YZ plane (x={mid_x})")
    im8 = axs[2,1].imshow(gpu_grid[:,:,mid_x], origin='lower', cmap='viridis', norm=norm)
    divider = make_axes_locatable(axs[2,1])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im8, cax=cax)
    
    axs[2,2].set_title(f"Difference: YZ plane")
    im9 = axs[2,2].imshow(gpu_grid[:,:,mid_x] - cpu_grid[:,:,mid_x], origin='lower', cmap='coolwarm')
    divider = make_axes_locatable(axs[2,2])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im9, cax=cax)
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(f"./data/{name}_comparison.png", dpi=150)
    print(f"Saved plot to ./data/{name}_comparison.png")

def compare_stages(name="NaCl_1x1_L1"):
    """Main function to compare different stages between CPU and GPU"""
    
    data_path = f"./data/{name}/"
    
    # 1. First compare Qgrid (charge density)
    print("\n=== Stage 1: Comparing Qgrid (Charge Density) ===")
    cpu_qgrid = load_grid_data(f"{data_path}Qgrid.npy")
    gpu_qgrid = load_grid_data(f"{data_path}Qgrid_gpu.npy")
    
    if cpu_qgrid is not None and gpu_qgrid is not None:
        qgrid_shift = compare_grid_positions(cpu_qgrid, gpu_qgrid, "Qgrid")
    else:
        print("Skipping Qgrid comparison due to missing data")
        qgrid_shift = (0, 0, 0)
    
    # 2. Next compare Vin_buff (after FFT/Poisson)
    print("\n=== Stage 2: Comparing Potential (VCoul) ===")
    cpu_vcoul = load_grid_data(f"{data_path}VCoul.npy")
    gpu_vcoul = load_grid_data(f"{data_path}V_Coul_gpu_before.npy")
    
    if cpu_vcoul is not None and gpu_vcoul is not None:
        # Handle dimension mismatch - transpose GPU data if needed
        print(f"CPU VCoul shape: {cpu_vcoul.shape}")
        print(f"GPU VCoul shape: {gpu_vcoul.shape}")
        
        if cpu_vcoul.shape != gpu_vcoul.shape:
            print("Detected dimension mismatch! Transposing GPU data for proper comparison...")
            if gpu_vcoul.shape == (40, 40, 200) and cpu_vcoul.shape == (200, 40, 40):
                # Transpose GPU data from XYZ to ZYX order
                gpu_vcoul = np.transpose(gpu_vcoul, (2, 1, 0))
                print(f"After transposing: GPU VCoul shape: {gpu_vcoul.shape}")
            else:
                print("Warning: Unexpected dimension mismatch, attempting general transpose")
                # Try to match dimensions by permuting axes
                for perm in [(0,1,2), (0,2,1), (1,0,2), (1,2,0), (2,0,1), (2,1,0)]:
                    transposed = np.transpose(gpu_vcoul, perm)
                    if transposed.shape == cpu_vcoul.shape:
                        gpu_vcoul = transposed
                        print(f"Found matching transpose with permutation {perm}")
                        print(f"After transposing: GPU VCoul shape: {gpu_vcoul.shape}")
                        break
        
        vcoul_shift = compare_grid_positions(cpu_vcoul, gpu_vcoul, "VCoul")
    else:
        print("Skipping VCoul comparison due to missing data")
        vcoul_shift = (0, 0, 0)
    
    # 3. Compare Poisson output (Vgrid)
    print("\n=== Stage 3: Comparing Poisson output (Vgrid) ===")
    cpu_vgrid = load_grid_data(f"{data_path}Vgrid.npy")
    gpu_vgrid = load_grid_data("./data/Vgrid_gpu.npy")
    if cpu_vgrid is not None and gpu_vgrid is not None:
        # No dimension mismatch expected
        vgrid_shift = compare_grid_positions(cpu_vgrid, gpu_vgrid, "Vgrid")
    else:
        print("Skipping Vgrid comparison due to missing data")
        vgrid_shift = (0, 0, 0)
    
    # 4. Compare Laplace output (Vin)
    print("\n=== Stage 4: Comparing Laplace output (Vin) ===")
    cpu_vin = load_grid_data(f"{data_path}Vin.npy")
    gpu_vin = load_grid_data("./data/Vin_gpu.npy")
    if cpu_vin is not None and gpu_vin is not None:
        # Ensure proper shape
        if cpu_vin.shape != gpu_vin.shape:
            print("Detected shape mismatch for Vin, attempting transpose if needed...")
            for perm in [(0,1,2),(2,1,0),(1,0,2)]:
                trans = np.transpose(gpu_vin, perm)
                if trans.shape == cpu_vin.shape:
                    gpu_vin = trans
                    print(f"Transposed gpu_vin with perm {perm}, new shape: {gpu_vin.shape}")
                    break
        vin_shift = compare_grid_positions(cpu_vin, gpu_vin, "Vin")
    else:
        print("Skipping Vin comparison due to missing data")
        vin_shift = (0, 0, 0)
    
    # 5. Display combined results
    print("\n=== Detailed Shift Analysis ===")
    print(f"Qgrid shift: {qgrid_shift}")
    print(f"VCoul shift: {vcoul_shift}")
    print(f"Vgrid shift: {vgrid_shift}")
    print(f"Vin shift: {vin_shift}")

if __name__ == "__main__":
    # Parse command line arguments
    import argparse
    parser = argparse.ArgumentParser(description="Debug grid shifts between CPU and GPU implementations")
    parser.add_argument("--name", default="NaCl_1x1_L1", help="Dataset name (default: NaCl_1x1_L1)")
    args = parser.parse_args()
    
    # Run the comparison
    compare_stages(args.name)
