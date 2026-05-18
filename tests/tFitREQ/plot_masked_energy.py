#!/usr/bin/env python3
"""
Plot masked energy decomposition using reference data and OpenCL model computation.
Properly loads reference energies, computes model with decomposition, and plots both.
"""
import argparse
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from pyBall.OCL.FittingDriver import FittingDriver
from pyBall.FitREQutils import plot_system_panel, _build_frame_grid, parse_xyz_with_headers

def main():
    parser = argparse.ArgumentParser(description='Plot masked energy visualization')
    #parser.add_argument('--xyz', default='/home/prokop/git/FireCore-fitREQH/tests/tFitREQ/add_epairs_full/H2O-A1_HF-D1-y.xyz',help='XYZ file to process')
    #parser.add_argument('--xyz', default='/home/prokop/git/FireCore-fitREQH/tests/tFitREQ/add_epairs_full/H2O-A1_H2O-D1-y.xyz',help='XYZ file to process')
    parser.add_argument('--xyz', default='/home/prokop/git/FireCore-fitREQH/tests/tFitREQ/add_epairs_full/HCN-A1_HF-D1-y.xyz',help='XYZ file to process')

    parser.add_argument('--atypes', default='/home/prokop/git/FireCore-fitREQH/tests/tFitREQ/data/AtomTypes.dat',help='Atom types file')
    parser.add_argument('--output', help='Output image file')
    parser.add_argument('--kcal', type=int, default=1, help='Convert to kcal/mol (default: 0.0 = no conversion)')
    parser.add_argument('--sym',  type=int, default=1, help='Symmetric color scale (like plot_ref.py --sym)')
    parser.add_argument('--same-scale', type=int, default=1, help='Force same color scale for Reference vs Model Etot panels')
    args = parser.parse_args()
    
    print(f"Loading data from {args.xyz}")
    
    # Step 1: Load reference energy grid from XYZ file using _build_frame_grid
    print("Step 1: Building reference energy grid from XYZ file...")
    V_ref, rv, A, shift, frame_idx, Ps_raw, Ts_raw = _build_frame_grid(args.xyz)
    
    print(f"Reference energy grid shape: {V_ref.shape}")
    print(f"Reference energy range: {np.min(V_ref):.6f} to {np.max(V_ref):.6f}")
    
    # Check for valid data
    if not np.any(np.isfinite(V_ref)):
        print("ERROR: All-NaN reference energy grid!")
        return
    
    # Step 2: Initialize FittingDriver and load atom data
    print("Step 2: Initializing FittingDriver...")
    drv = FittingDriver()
    drv.load_atom_types(args.atypes)
    drv.load_data(args.xyz)
    
    # Step 3: Compile with MODEL_MorseQ_PAIR_DECOMP macro
    print("Step 3: Compiling with Morse decomposition macro...")
    
    # Read Forces.cl content to extract macro
    forces_path = "/home/prokop/git/FireCore-fitREQH/cpp/common_resources/cl/Forces.cl"
    with open(forces_path, 'r') as f:
        forces_content = f.read()
    
    # Extract the actual macro content for MODEL_MorseQ_PAIR_DECOMP
    start_idx = forces_content.find('//>>>macro MODEL_MorseQ_PAIR_DECOMP')
    if start_idx == -1:
        raise ValueError("Could not find MODEL_MorseQ_PAIR_DECOMP macro start")
    
    open_brace = forces_content.find('{', start_idx)
    if open_brace == -1:
        raise ValueError("Could not find opening brace for MODEL_MorseQ_PAIR_DECOMP macro")
    
    # Find the matching closing brace
    brace_count = 0
    pos = open_brace
    while pos < len(forces_content):
        if forces_content[pos] == '{':
            brace_count += 1
        elif forces_content[pos] == '}':
            brace_count -= 1
            if brace_count == 0:
                close_brace = pos
                break
        pos += 1
    
    if brace_count != 0:
        raise ValueError("Could not find matching closing brace for MODEL_MorseQ_PAIR_DECOMP macro")
    
    morse_decomp = forces_content[open_brace + 1:close_brace].strip()
    
    # Use the actual macro content
    drv.compile_with_model(macros={
        'MODEL_PAIR_DECOMP': morse_decomp
    })
    drv.init_and_upload_energy_only()
    
    # Get dimensions
    ni = drv.host_ranges[0][2]  # atoms in fragment 1
    nj = drv.host_ranges[0][3]  # atoms in fragment 2
    print(f"Interaction matrix dimensions: {ni} x {nj}")

    # Step 4: Compute model energies for all samples on GPU (parallel)
    print("Step 4: Computing model energies for all samples (GPU parallel)...")

    # Number of samples is defined by loaded dataset
    n_samples = int(drv.n_samples)
    print(f"n_samples = {n_samples}")

    # Create masks (component weights)
    mask_all      = drv.create_mask(ni, nj, pauli=1.0, london=1.0, electro=1.0, hbond=1.0)
    mask_hbond    = drv.create_mask(ni, nj, pauli=0.0, london=0.0, electro=0.0, hbond=1.0)
    mask_baseline = drv.create_mask(ni, nj, pauli=1.0, london=1.0, electro=1.0, hbond=0.0)

    # IMPORTANT: evaluate_energies_masked returns (Emols_total, Emols_masked)
    # We must use the *masked* result for decomposition plots.
    _, Em_all      = drv.evaluate_energies_masked(mask_all)
    _, Em_hbond    = drv.evaluate_energies_masked(mask_hbond)
    _, Em_baseline = drv.evaluate_energies_masked(mask_baseline)

    assert Em_all.shape[0] == n_samples
    assert Em_hbond.shape[0] == n_samples
    assert Em_baseline.shape[0] == n_samples

    # Step 5: Energy decomposition completed
    print("Step 5: Energy decomposition completed")
    
    # Step 6: Reshape model energies to match reference grid shape
    print("Step 6: Reshaping model energies to match grid...")
    ny, nx = V_ref.shape

    # The frame_idx from _build_frame_grid tells us how to map samples -> grid positions
    V_model_total    = np.full((ny, nx), np.nan, dtype=np.float32)
    V_model_hbond    = np.full((ny, nx), np.nan, dtype=np.float32)
    V_model_baseline = np.full((ny, nx), np.nan, dtype=np.float32)

    for iy in range(ny):
        for ix in range(nx):
            isamp = int(frame_idx[iy, ix])
            if isamp >= 0 and isamp < n_samples:
                V_model_total[iy, ix]    = Em_all[isamp]
                V_model_hbond[iy, ix]    = Em_hbond[isamp]
                V_model_baseline[iy, ix] = Em_baseline[isamp]
    
    # Compute model Eout = total - hbond
    V_model_eout = V_model_total - V_model_hbond
    
    # Step 7: Plot reference and model energies with proper vmin/vmax
    print("Step 7: Plotting energy components...")
    
    # Create figure with 4 columns
    fig, axes = plt.subplots(1, 4, figsize=(16, 4))
    
    # Set symmetric color limits: vmin = 1.2 * min(E_ref), vmax = -vmin
    # Apply to ALL panels for consistent comparison
    # Note: plot_imshow handles eV->kcal conversion (23.060548) when kcal=True
    vmin_ref = np.nanmin(V_ref)
    if args.kcal:
        vmin_ref *= 23.060548  # Convert eV to kcal/mol (to match plot_imshow's conversion)
    vmin = 1.2 * vmin_ref
    vmax = -vmin
    
    print(f"Symmetric color limits: vmin={vmin:.6f}, vmax={vmax:.6f}")
    
    # Plot each energy component using plot_system_panel with proper vmin/vmax
    plot_system_panel(V_ref, rv, A, axes[0], 'Reference\n(Total)', args.kcal, sym=False, overlay_rmin=False)
    plot_system_panel(V_model_total, rv, A, axes[1], 'Model Etot\n(Total)', args.kcal, sym=False, overlay_rmin=False)
    plot_system_panel(V_model_hbond, rv, A, axes[2], 'Model Ein\n(H-bond)', args.kcal, sym=False, overlay_rmin=False)
    plot_system_panel(V_model_eout, rv, A, axes[3], 'Model Eout\n(Baseline)', args.kcal, sym=False, overlay_rmin=False)

    # Apply symmetric color limits to ALL panels
    for ax in axes:
        if ax.images:
            ax.images[0].set_clim(vmin, vmax)
    
    plt.suptitle(f'Energy Decomposition: {os.path.basename(args.xyz)}', fontsize=14)
    plt.tight_layout()
    
    if args.output:
        plt.savefig(args.output, dpi=200, bbox_inches='tight')
        print(f"Saved: {args.output}")
    else:
        plt.show()
    
    # Print statistics
    print("\nEnergy Statistics:")
    print(f"  Reference: min={np.nanmin(V_ref):.6f}, max={np.nanmax(V_ref):.6f}, mean={np.nanmean(V_ref):.6f}")
    print(f"  Model Total: min={np.nanmin(V_model_total):.6f}, max={np.nanmax(V_model_total):.6f}, mean={np.nanmean(V_model_total):.6f}")
    print(f"  Model H-bond: min={np.nanmin(V_model_hbond):.6f}, max={np.nanmax(V_model_hbond):.6f}, mean={np.nanmean(V_model_hbond):.6f}")
    print(f"  Model Eout: min={np.nanmin(V_model_eout):.6f}, max={np.nanmax(V_model_eout):.6f}, mean={np.nanmean(V_model_eout):.6f}")
    
    # Compute correlation between reference and model
    valid_mask = np.isfinite(V_ref) & np.isfinite(V_model_total)
    if np.any(valid_mask):
        correlation = np.corrcoef(V_ref[valid_mask], V_model_total[valid_mask])[0, 1]
        print(f"\nReference-Model correlation: {correlation:.6f}")

if __name__ == "__main__":
    main()
