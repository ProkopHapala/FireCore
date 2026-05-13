#!/usr/bin/env python3
"""
Test GPU gradient computation with visualization.
Compares 4 methods: analytical, numpy, clFFT, cl central differences.
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

# Add FireCore to path
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.realpath(os.path.join(_THIS_DIR, '..', '..'))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from pyBall.OCL.AFM import AFMulator

# ═══════════════════════════════════════════════════════════════════════════════
# Analytical Potentials and Their Gradients
# ═══════════════════════════════════════════════════════════════════════════════

def gaussian_potential(x, y, z, sigma=1.0):
    """E = exp(-(x²+y²+z²)/(2σ²))"""
    r2 = x*x + y*y + z*z
    return np.exp(-r2 / (2*sigma*sigma))

def gaussian_gradient(x, y, z, sigma=1.0):
    """∇E = -(x,y,z)/σ² * exp(-r²/(2σ²))"""
    r2 = x*x + y*y + z*z
    exp_val = np.exp(-r2 / (2*sigma*sigma))
    factor = -exp_val / (sigma*sigma)
    return factor*x, factor*y, factor*z

def periodic_potential(x, y, z, L=6.4):
    """E = cos(2π*x/L) + cos(2π*y/L) + cos(2π*z/L) - truly periodic with period L"""
    return np.cos(2*np.pi*x/L) + np.cos(2*np.pi*y/L) + np.cos(2*np.pi*z/L)

def periodic_gradient(x, y, z, L=6.4):
    """∇E = -(2π/L)sin(2π*x/L) + ..."""
    factor = -2*np.pi/L
    return factor*np.sin(2*np.pi*x/L), factor*np.sin(2*np.pi*y/L), factor*np.sin(2*np.pi*z/L)

# ═══════════════════════════════════════════════════════════════════════════════
# Gradient Computation Methods
# ═══════════════════════════════════════════════════════════════════════════════

def compute_analytical_gradient(grad_func, X, Y, Z, **kwargs):
    """Compute gradient using analytical function."""
    gx, gy, gz = grad_func(X, Y, Z, **kwargs)
    return np.stack([gx, gy, gz], axis=-1)

def compute_numpy_gradient(field, step):
    """Compute gradient using numpy with periodic boundary conditions."""
    # np.gradient doesn't support true periodic BC, so we implement manually
    # Pad with periodic wrapping, compute gradient, then unpad
    pad_width = 1
    field_padded = np.pad(field, pad_width, mode='wrap')
    
    # Central differences on padded array
    grad_padded = np.stack([np.gradient(field_padded, step, axis=i) for i in range(3)], axis=-1)
    
    # Remove padding
    grad = grad_padded[pad_width:-pad_width, pad_width:-pad_width, pad_width:-pad_width, :]
    return grad

def compute_cl_central_gradient(afmulator, E_field, step):
    """Compute gradient using OpenCL central difference kernel."""
    grads_gpu = afmulator.compute_gradient_cl(E_field, step, bAlloc=True)
    # grads_gpu is (nx, ny, nz, 4) with (Fx, Fy, Fz, E)
    # Force = -gradient, so gradients are negative of force components
    grads_gpu_f = -grads_gpu[:, :, :, :3]
    return grads_gpu_f

def compute_cl_fft_gradient(afmulator, E_field, step):
    """Compute gradient using OpenCL FFT spectral differentiation."""
    grads_gpu = afmulator.compute_gradient_fft_cl(E_field, step, bAlloc=True)
    # grads_gpu is (nx, ny, nz, 4) with (Fx, Fy, Fz, E)
    # Force = -gradient, so gradients are negative of force components
    grads_gpu_f = -grads_gpu[:, :, :, :3]
    return grads_gpu_f

# ═══════════════════════════════════════════════════════════════════════════════
# Visualization
# ═══════════════════════════════════════════════════════════════════════════════

def plot_gradient_comparison(grads_dict, axes, z_slice_idx=None, output_dir=None):
    """
    Plot gradient comparison for all methods.
    
    Args:
        grads_dict: dict with keys like 'analytical', 'numpy', 'cl_central', 'cl_fft'
                   Each value is (nx, ny, nz, 3) array
        axes: tuple of (x, y, z) coordinate arrays
        z_slice_idx: index of z-slice to plot (default: middle)
        output_dir: directory to save plots
    """
    if z_slice_idx is None:
        nz = axes[2].shape[0]
        z_slice_idx = nz // 2
    
    methods = [m for m in grads_dict.keys() if grads_dict[m] is not None]
    n_methods = len(methods)
    
    if n_methods == 0:
        print("No valid gradients to plot")
        return None
    
    # Create figure with n_methods columns
    fig, axes_plot = plt.subplots(3, n_methods, figsize=(4*n_methods, 12))
    if n_methods == 1:
        axes_plot = axes_plot.reshape(3, 1)
    fig.suptitle('Gradient Comparison (z-slice)', fontsize=14)
    
    # Plot x-component for each method
    for i, method in enumerate(methods):
        ax = axes_plot[0, i]
        data = grads_dict[method][:, :, z_slice_idx, 0]
        
        vabs = max(abs(data.min()), abs(data.max()), 1e-10)
        im = ax.imshow(data.T, origin='lower', cmap='seismic', vmin=-vabs, vmax=vabs)
        ax.set_title(f'{method}\nFx (z={z_slice_idx})')
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    
    # Plot y-component for each method
    for i, method in enumerate(methods):
        ax = axes_plot[1, i]
        data = grads_dict[method][:, :, z_slice_idx, 1]
        
        vabs = max(abs(data.min()), abs(data.max()), 1e-10)
        im = ax.imshow(data.T, origin='lower', cmap='seismic', vmin=-vabs, vmax=vabs)
        ax.set_title(f'{method}\nFy')
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    
    # Plot z-component for each method
    for i, method in enumerate(methods):
        ax = axes_plot[2, i]
        data = grads_dict[method][:, :, z_slice_idx, 2]
        
        vabs = max(abs(data.min()), abs(data.max()), 1e-10)
        im = ax.imshow(data.T, origin='lower', cmap='seismic', vmin=-vabs, vmax=vabs)
        ax.set_title(f'{method}\nFz')
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    
    plt.tight_layout()
    
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        fname = os.path.join(output_dir, 'gradient_comparison.png')
        plt.savefig(fname, dpi=150, bbox_inches='tight')
        print(f"Saved plot to {fname}")
    
    # Create separate figure for differences
    if 'analytical' in methods:
        ref = grads_dict['analytical']
        diff_methods = [m for m in methods if m != 'analytical']
        
        if len(diff_methods) > 0:
            fig_diff, axes_diff = plt.subplots(3, len(diff_methods), figsize=(4*len(diff_methods), 12))
            if len(diff_methods) == 1:
                axes_diff = axes_diff.reshape(3, 1)
            fig_diff.suptitle('Gradient Differences vs Analytical', fontsize=14)
            
            for i, method in enumerate(diff_methods):
                # x-component difference
                ax = axes_diff[0, i]
                diff = grads_dict[method][:, :, z_slice_idx, 0] - ref[:, :, z_slice_idx, 0]
                vabs = max(abs(diff.min()), abs(diff.max()), 1e-10)
                im = ax.imshow(diff.T, origin='lower', cmap='seismic', vmin=-vabs, vmax=vabs)
                ax.set_title(f'{method}\nFx diff')
                plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
                
                # y-component difference
                ax = axes_diff[1, i]
                diff = grads_dict[method][:, :, z_slice_idx, 1] - ref[:, :, z_slice_idx, 1]
                vabs = max(abs(diff.min()), abs(diff.max()), 1e-10)
                im = ax.imshow(diff.T, origin='lower', cmap='seismic', vmin=-vabs, vmax=vabs)
                ax.set_title(f'{method}\nFy diff')
                plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
                
                # z-component difference
                ax = axes_diff[2, i]
                diff = grads_dict[method][:, :, z_slice_idx, 2] - ref[:, :, z_slice_idx, 2]
                vabs = max(abs(diff.min()), abs(diff.max()), 1e-10)
                im = ax.imshow(diff.T, origin='lower', cmap='seismic', vmin=-vabs, vmax=vabs)
                ax.set_title(f'{method}\nFz diff')
                plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
            
            plt.tight_layout()
            fname_diff = os.path.join(output_dir, 'gradient_differences.png')
            plt.savefig(fname_diff, dpi=150, bbox_inches='tight')
            print(f"Saved difference plot to {fname_diff}")
    
    return fig

def main():
    print("="*60)
    print("GPU Gradient Kernel Visualization Test")
    print("="*60)
    
    # Initialize AFMulator
    print("\nInitializing AFMulator...")
    afmulator = AFMulator(nloc=8)
    
    # Test configurations - multiple grid sizes including non-cubic
    test_configs = [
        (32, 32, 32, 0.2),   # Small cubic grid
        (64, 64, 64, 0.1),   # Medium cubic grid
        (128, 64, 32, 0.05), # Non-cubic: x > y > z
        (64, 128, 32, 0.05), # Non-cubic: y > x > z
        (32, 64, 128, 0.05), # Non-cubic: z > y > x
        (128, 128, 64, 0.05), # Larger non-cubic
    ]
    
    for nx, ny, nz, step in test_configs:
        origin = (0, 0, 0)
        print(f"\n{'='*60}")
        print(f"Testing grid: {nx}x{ny}x{nz}, step={step}, origin={origin}")
        print(f"{'='*60}")
        
        # Create grid
        x = np.arange(nx) * step + origin[0]
        y = np.arange(ny) * step + origin[1]
        z = np.arange(nz) * step + origin[2]
        X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
        
        # Center coordinates for Gaussian - make it smaller so it fits in grid
        cx, cy, cz = (nx*step/2, ny*step/2, nz*step/2)
        Xc, Yc, Zc = X - cx, Y - cy, Z - cz
        
        # Test with Gaussian potential - use smaller sigma so it decays to zero before boundaries
        print("\nComputing Gaussian potential...")
        sigma = 0.3  # Even smaller to ensure true zero at boundaries
        E_field = gaussian_potential(Xc, Yc, Zc, sigma=sigma)
        print(f"  E_field range: [{E_field.min():.6f}, {E_field.max():.6f}]")
        print(f"  E_field at boundaries: min={E_field.min():.6e}")
        
        # Compute gradients with different methods
        print("\nComputing gradients with different methods...")
        
        # 1. Analytical
        print("  1. Analytical gradient...")
        grads_analytical = compute_analytical_gradient(gaussian_gradient, Xc, Yc, Zc, sigma=sigma)
        
        # 2. NumPy central differences
        print("  2. NumPy gradient (central differences)...")
        grads_numpy = compute_numpy_gradient(E_field, step)
        
        # 3. OpenCL central differences
        print("  3. OpenCL gradient (central differences)...")
        grads_cl_central = compute_cl_central_gradient(afmulator, E_field, step)
        
        # 4. OpenCL FFT spectral differentiation
        print("  4. OpenCL gradient (FFT)...")
        try:
            grads_cl_fft = compute_cl_fft_gradient(afmulator, E_field, step)
        except Exception as e:
            print(f"    FFT gradient failed: {e}")
            grads_cl_fft = None
        
        # Store results
        grads_dict = {
            'analytical': grads_analytical,
            'numpy': grads_numpy,
            'cl_central': grads_cl_central,
            'cl_fft': grads_cl_fft,
        }
        
        # Print error metrics
        print("\n" + "="*60)
        print("Error Metrics")
        print("="*60)
        
        def print_error(name, gpu, ref):
            diff = gpu - ref
            abs_error = np.abs(diff)
            interior_slice = slice(2, -2)
            interior_error = abs_error[interior_slice, interior_slice, interior_slice, :]
            
            print(f"\n{name}:")
            print(f"  Max error (all):     {np.max(abs_error):.6e}")
            print(f"  Mean error (all):    {np.mean(abs_error):.6e}")
            print(f"  Max error (interior): {np.max(interior_error):.6e}")
            print(f"  Mean error (interior): {np.mean(interior_error):.6e}")
        
        print_error('NumPy vs Analytical', grads_numpy, grads_analytical)
        print_error('CL Central vs Analytical', grads_cl_central, grads_analytical)
        print_error('CL Central vs NumPy', grads_cl_central, grads_numpy)
        
        # Create visualization
        print("\nCreating visualization...")
        output_dir = '/home/prokop/git/FireCore/tests/tAFM/gradient_plots'
        fig = plot_gradient_comparison(grads_dict, (x, y, z), z_slice_idx=nz//2, output_dir=output_dir)
        
        # Save with gaussian-specific name
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
            fname_gauss = os.path.join(output_dir, f'gradient_comparison_gaussian_{nx}x{ny}x{nz}.png')
            plt.savefig(fname_gauss, dpi=150, bbox_inches='tight')
            print(f"Saved Gaussian plot to {fname_gauss}")
    
        # Also test with periodic cosine potential
        print("\n" + "="*60)
        print("Testing with Periodic Cosine Potential")
        print("="*60)
        
        # Use grid from 0 to L where L is exactly the period
        # For true periodicity, we need the grid to cover exactly one period
        # so that first and last points connect smoothly via periodic BC
        L = nx * step  # Grid length = one period
        origin_per = (0, 0, 0)
        x_per = np.arange(nx) * step + origin_per[0]
        y_per = np.arange(ny) * step + origin_per[1]
        z_per = np.arange(nz) * step + origin_per[2]
        X_per, Y_per, Z_per = np.meshgrid(x_per, y_per, z_per, indexing='ij')
        
        print(f"  Periodic grid: {nx}x{ny}x{nz}, step={step:.6f}, L={L:.6f}")
        print(f"  Periodicity check: cos(0)={np.cos(0):.6f}, cos(L)={np.cos(2*np.pi):.6f}")
        
        E_per = periodic_potential(X_per, Y_per, Z_per, L=L)
        print(f"  E range: [{E_per.min():.6f}, {E_per.max():.6f}]")
        
        grads_per_analytical = compute_analytical_gradient(periodic_gradient, X_per, Y_per, Z_per, L=L)
        grads_per_numpy = compute_numpy_gradient(E_per, step)
        grads_per_cl_central = compute_cl_central_gradient(afmulator, E_per, step)
        
        # FFT gradient
        print("  Computing FFT gradient...")
        try:
            grads_per_fft = compute_cl_fft_gradient(afmulator, E_per, step)
        except Exception as e:
            print(f"    FFT gradient failed: {e}")
            grads_per_fft = None
        
        # Debug boundary values
        print(f"  E[0,0,0] = {E_per[0,0,0]:.6f}, E[-1,0,0] = {E_per[-1,0,0]:.6f}")
        print(f"  grad_numpy[0,0,0,0] = {grads_per_numpy[0,0,0,0]:.6f}, grad_numpy[-1,0,0,0] = {grads_per_numpy[-1,0,0,0]:.6f}")
        print(f"  grad_cl[0,0,0,0] = {grads_per_cl_central[0,0,0,0]:.6f}, grad_cl[-1,0,0,0] = {grads_per_cl_central[-1,0,0,0]:.6f}")
        print(f"  grad_ana[0,0,0,0] = {grads_per_analytical[0,0,0,0]:.6f}, grad_ana[-1,0,0,0] = {grads_per_analytical[-1,0,0,0]:.6f}")
        
        grads_per_dict = {
            'analytical': grads_per_analytical,
            'numpy': grads_per_numpy,
            'cl_central': grads_per_cl_central,
            'cl_fft': grads_per_fft,
        }
        
        fig_per = plot_gradient_comparison(grads_per_dict, (x_per, y_per, z_per), 
                                            z_slice_idx=nz//2, output_dir=output_dir)
        
        # Save with periodic-specific name
        if output_dir:
            fname_per = os.path.join(output_dir, f'gradient_comparison_periodic_{nx}x{ny}x{nz}.png')
            plt.savefig(fname_per, dpi=150, bbox_inches='tight')
            print(f"Saved Periodic plot to {fname_per}")
    
    print("\n" + "="*60)
    print("Test completed!")
    print("="*60)

if __name__ == "__main__":
    main()
