#!/usr/bin/env python3
import sys
import os
import numpy as np
import matplotlib.pyplot as plt

# Add FireCore to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from pyBall.OCL.FittingDriver import FittingDriver
from pyBall.GUI.FitREQUtil import plot_stacked_1d, reshape_to_grid_proper, EV_TO_KCAL

def run_cli(xyz_file, atom_types_file, model_name='ENERGY_MorseQ_PAIR', group1="0,1,2", group2="0,1,2", out_dir="plots_headless"):
    """Headless CLI for FitREQ diagnostic mapping."""
    if not os.path.exists(out_dir): os.makedirs(out_dir)
    
    print(f"--- FitREQ Headless Diagnostics ---")
    print(f"  XYZ: {xyz_file}")
    print(f"  AtomTypes: {atom_types_file}")
    print(f"  Model: {model_name}")
    print(f"  Groups: G1=[{group1}], G2=[{group2}]")
    
    # 1. Initialize Driver
    drv = FittingDriver(verbose=1)
    drv.load_atom_types(atom_types_file)
    drv.compile_all_with_model(model_name, bPrint=False)
    drv.load_data(xyz_file)
    drv.init_and_upload_energy_only()
    drv.setup_energy_kernel()
    
    # 2. Load Geometry and Reshape
    r, a, rows, Ps = drv.load_scan_geometry(xyz_file)
    Erefs = drv.host_ErefW if hasattr(drv, 'host_ErefW') else None
    if Erefs is not None:
        if Erefs.ndim > 1: Erefs = Erefs[:, 0]
        Erefs = Erefs * EV_TO_KCAL # Convert to kcal/mol
    
    Vref, Rg, Arow, rv, rows = reshape_to_grid_proper(Erefs, r, a, rows)
    
    drv.apply_type_overrides() # Ensure defaults are applied
    Em = drv.evaluate_energies() * EV_TO_KCAL # Convert to kcal/mol
    Vmod, _, _, _, _ = reshape_to_grid_proper(Em, r, a, rows)
    
    # 3. Find Global Minimum of Reference
    if Vref is not None and np.any(np.isfinite(Vref)):
        min_idx = np.nanargmin(Vref)
        row, col = np.unravel_index(min_idx, Vref.shape)
    else:
        row, col = 0, 0
    
    print(f"Global minimum of reference at: row={row}, col={col} (r={rv[col]:.3f} A, a={Arow[row]:.1f} deg)")
    
    # 4. Plot 2D Maps
    fig, axs = plt.subplots(1, 3, figsize=(18, 5))
    
    # Calculate symmetric limits based on reference minimum
    vmin_ref = np.nanmin(Vref) if Vref is not None else np.nanmin(Vmod)
    vmin_val = 1.2 * vmin_ref if vmin_ref < 0 else -10.0 # fallback if min is positive or 0
    vmax_val = -vmin_val
    vdif_max = 0.10 * vmax_val
    
    ny, nx = Vref.shape
    
    # Use imshow to plot pixel grid and avoid non-monotonic grid issues
    im1 = axs[0].imshow(Vref, origin='lower', aspect='auto', cmap='seismic', vmin=vmin_val, vmax=vmax_val)
    axs[0].plot(col, row, 'wo', ms=8, mew=1.5, mfc='none') # Hollow circle for visibility
    axs[0].set_title('Reference Energy (kcal/mol)')
    plt.colorbar(im1, ax=axs[0])
    
    im2 = axs[1].imshow(Vmod, origin='lower', aspect='auto', cmap='seismic', vmin=vmin_val, vmax=vmax_val)
    axs[1].plot(col, row, 'wo', ms=8, mew=1.5, mfc='none')
    axs[1].set_title('Model (kcal/mol)')
    plt.colorbar(im2, ax=axs[1])
    
    Vdiff = Vmod - Vref
    im3 = axs[2].imshow(Vdiff, origin='lower', aspect='auto', cmap='bwr', vmin=-vdif_max, vmax=vdif_max)
    axs[2].plot(col, row, 'ko', ms=8, mew=1.5, mfc='none')
    axs[2].set_title('Difference (kcal/mol)')
    plt.colorbar(im3, ax=axs[2])
    
    for ax in axs:
        ax.set_xlabel('Distance (A)')
        ax.set_ylabel('Angle (deg)')
        
        # Ticks every 5th pixel for distance
        tick_indices = np.arange(0, nx, 5)
        ax.set_xticks(tick_indices)
        ax.set_xticklabels([f"{rv[i]:.2f}" for i in tick_indices])
        
        # Ticks for angle
        y_tick_indices = np.arange(0, ny, max(1, ny // 10))
        ax.set_yticks(y_tick_indices)
        ax.set_yticklabels([f"{Arow[i]:.1f}" for i in y_tick_indices])
    
    plt.tight_layout()
    plt.savefig(f"{out_dir}/energy_maps.png")
    plt.close()
    
    # 5. Plot 1D Cuts with Decomposition
    mask1 = [int(i.strip()) for i in group1.split(',') if i.strip()]
    mask2 = [int(i.strip()) for i in group2.split(',') if i.strip()]
    
    # Radial Cut
    s, e = rows[row]
    indices = list(range(s, e))
    Ms_rad = drv.evaluate_interaction_matrices_along_cut(indices)
    E_total_rad = np.array([M.sum() for M in Ms_rad]) * EV_TO_KCAL
    E_group_rad = np.array([drv.sum_interaction_matrix(M, mask1, mask2).sum() for M in Ms_rad]) * EV_TO_KCAL
    E_other_rad = E_total_rad - E_group_rad
    rv_cut = rv[:len(indices)]
    Vref_cut = Vref[row, :len(indices)] if Vref is not None else None
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    plot_stacked_1d(ax1, rv_cut, E_total_rad, [E_group_rad, E_other_rad], 
                   ['Group Contribution', 'Other'], ['green', 'blue'], ref=Vref_cut,
                   title=f'Radial Cut @ Angle={Arow[row]:.1f} deg',
                   xlabel='Distance (A)', ylabel='Energy (kcal/mol)')
    
    # Angular Cut
    ang_indices = []
    for r_idx, (s, e) in enumerate(rows):
        if col < (e - s): ang_indices.append(s + col)
        else: ang_indices.append(-1)
    valid_rows = [i for i, idx in enumerate(ang_indices) if idx >= 0]
    valid_indices = [ang_indices[i] for i in valid_rows]
    
    if valid_indices:
        Ms_ang = drv.evaluate_interaction_matrices_along_cut(valid_indices)
        E_total_ang = np.array([M.sum() for M in Ms_ang]) * EV_TO_KCAL
        E_group_ang = np.array([drv.sum_interaction_matrix(M, mask1, mask2).sum() for M in Ms_ang]) * EV_TO_KCAL
        E_other_ang = E_total_ang - E_group_ang
        A_cut = Arow[valid_rows]
        Vref_ang = Vref[valid_rows, col] if Vref is not None else None
        
        plot_stacked_1d(ax2, A_cut, E_total_ang, [E_group_ang, E_other_ang], 
                       ['Group Contribution', 'Other'], ['green', 'blue'], ref=Vref_ang,
                       title=f'Angular Cut @ Dist={rv[col]:.2f} A',
                       xlabel='Angle (deg)', ylabel='Energy (kcal/mol)')
    
    plt.tight_layout()
    plt.savefig(f"{out_dir}/1d_cuts.png")
    plt.close()
    
    # 6. Interaction Matrix at Minimum
    # Map pixel (row, col) to sample index
    idx_map = {}
    curr = 0
    for r_i in range(Vref.shape[0]):
        for c_i in range(Vref.shape[1]):
            if curr < len(Ps):
                idx_map[(r_i, c_i)] = curr
                curr += 1
    
    idx_min = idx_map.get((row, col), 0)
    interactions = drv.evaluate_interaction_matrix(idx_min) * EV_TO_KCAL # Convert to kcal/mol
    
    fig, axs = plt.subplots(2, 2, figsize=(10, 10))
    im_p = axs[0,0].imshow(interactions[:,:,0], cmap='hot', origin='lower'); axs[0,0].set_title('Pauli Repulsion (kcal/mol)'); plt.colorbar(im_p, ax=axs[0,0])
    im_l = axs[0,1].imshow(interactions[:,:,1], cmap='cool', origin='lower'); axs[0,1].set_title('London Dispersion (kcal/mol)'); plt.colorbar(im_l, ax=axs[0,1])
    im_e = axs[1,0].imshow(interactions[:,:,2], cmap='seismic', origin='lower'); axs[1,0].set_title('Electrostatic (kcal/mol)'); plt.colorbar(im_e, ax=axs[1,0])
    im_h = axs[1,1].imshow(interactions[:,:,3], cmap='viridis', origin='lower'); axs[1,1].set_title('H-Bond (kcal/mol)'); plt.colorbar(im_h, ax=axs[1,1])
    
    plt.suptitle(f"Interaction Matrix @ r={rv[col]:.2f}, a={Arow[row]:.1f}")
    plt.tight_layout()
    plt.savefig(f"{out_dir}/interaction_matrix.png")
    plt.close()
    
    print(f"--- SUCCESS: Plots saved to '{out_dir}/' ---")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="FitREQ Headless Diagnostic Tool")
    parser.add_argument("--xyz", type=str, default="../../tests/tFitREQ_PN/wb97m-split/H2O-A1_H2O-D1-y.xyz", help="Input XYZ scan file")
    parser.add_argument("--atypes", type=str, default="../../tests/tFitREQ_PN/data/AtomTypes.dat", help="Atom types parameter file")
    parser.add_argument("--model", type=str, default="ENERGY_MorseQ_PAIR", help="Energy model macro name")
    parser.add_argument("--g1", type=str, default="0,1,2", help="Atom indices in group 1 (comma separated)")
    parser.add_argument("--g2", type=str, default="0,1,2", help="Atom indices in group 2 (comma separated)")
    parser.add_argument("--out", type=str, default="plots_headless", help="Output directory")
    args = parser.parse_args()
    
    # Resolve relative paths relative to script location if needed
    base_dir = os.path.dirname(os.path.abspath(__file__))
    xyz = os.path.abspath(os.path.join(base_dir, args.xyz)) if not os.path.isabs(args.xyz) else args.xyz
    atypes = os.path.abspath(os.path.join(base_dir, args.atypes)) if not os.path.isabs(args.atypes) else args.atypes
    
    run_cli(xyz, atypes, args.model, args.g1, args.g2, args.out)
