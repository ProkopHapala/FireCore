#!/usr/bin/env python3
"""
Plot masked energy decomposition using reference data and OpenCL model computation.
Properly loads reference energies, computes model with decomposition, and plots both.
"""
import argparse
import json
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from pyBall.OCL.FittingDriver import FittingDriver
from pyBall.FitREQutils import (plot_system_panel, plot_polar_symmetric, _build_frame_grid, 
                                 parse_xyz_with_headers, plot_energy_panel, plot_profiles, plot_molecule,
                                 extract_min_curves)

def plot_profile_row(fig, axes, V_ref, V_model_total, V_model_hbond, V_model_eout, rv, A, 
                    frame_idx, Ps_raw, Ts_raw, n0_first, kcal, Rmax1D):
    """Plot second row with radial slice, angular slice, E_min(angle), and geometry."""
    # Use existing axes from row 2
    ax_radial = axes[1, 0]
    ax_angular = axes[1, 1]
    ax_emin = axes[1, 2]
    ax_geom = axes[1, 3]
    
    # Find global minimum indices from Reference
    iy_g, ix_g = np.unravel_index(np.nanargmin(V_ref), V_ref.shape)
    angle_min = A[iy_g]
    r_min = rv[ix_g] if ix_g < len(rv) else np.nan
    
    # 1) Radial slice at minimum angle (Eref, Etot, Ein, Eout)
    fac = 23.060548 if kcal else 1.0
    ax_radial.plot(rv, V_ref[iy_g, :] * fac, 'k:', lw=1.5, label='Eref')
    ax_radial.plot(rv, V_model_total[iy_g, :] * fac, 'r-', lw=0.5, label='Etot')
    ax_radial.plot(rv, V_model_hbond[iy_g, :] * fac, 'b-', lw=0.5, label='Ein')
    ax_radial.plot(rv, V_model_eout[iy_g, :] * fac, 'g-', lw=0.5, label='Eout')
    ax_radial.axvline(r_min, color='gray', linestyle='--', alpha=0.5, label='r_min')
    ax_radial.set_xlabel('Distance (Å)')
    ax_radial.set_ylabel('E (kcal/mol)' if kcal else 'E (eV)')
    ax_radial.set_title(f'Radial slice @ {angle_min:.1f}°')
    ax_radial.set_xlim(1.4, Rmax1D)
    ax_radial.legend(fontsize=8)
    ax_radial.grid(alpha=0.3)
    
    # 2) Angular slice at minimum radius (Eref, Etot, Ein, Eout)
    o = np.argsort(A)
    ax_angular.plot(A[o], V_ref[:, ix_g][o] * fac, 'k:', lw=1.5, label='Eref')
    ax_angular.plot(A[o], V_model_total[:, ix_g][o] * fac, 'r-', lw=0.5, label='Etot')
    ax_angular.plot(A[o], V_model_hbond[:, ix_g][o] * fac, 'b-', lw=0.5, label='Ein')
    ax_angular.plot(A[o], V_model_eout[:, ix_g][o] * fac, 'g-', lw=0.5, label='Eout')
    ax_angular.axvline(angle_min, color='gray', linestyle='--', alpha=0.5, label='angle_min')
    ax_angular.set_xlabel('Angle (deg)')
    ax_angular.set_ylabel('E (kcal/mol)' if kcal else 'E (eV)')
    ax_angular.set_title(f'Angular slice @ r={r_min:.2f}Å')
    ax_angular.legend(fontsize=8)
    ax_angular.grid(alpha=0.3)
    
    # 3) E_min(angle) using plot_energy_panel for reference
    plot_energy_panel(V_ref, rv, A, ax_emin, kcal)
    # Overlay model curves on top
    _, emin_tot = extract_min_curves(A, rv, V_model_total.T)
    _, emin_hbond = extract_min_curves(A, rv, V_model_hbond.T)
    _, emin_eout = extract_min_curves(A, rv, V_model_eout.T)
    ax_emin.plot(A[o], emin_tot[o] * fac, 'r-', lw=0.5, label='Etot')
    ax_emin.plot(A[o], emin_hbond[o] * fac, 'b-', lw=0.5, label='Ein')
    ax_emin.plot(A[o], emin_eout[o] * fac, 'g-', lw=0.5, label='Eout')
    ax_emin.legend(fontsize=8)
    
    # 4) Geometry at global minimum using plot_molecule (2D)
    iframe = frame_idx[iy_g, ix_g]
    if iframe >= 0 and iframe < len(Ps_raw):
        apos_mol = Ps_raw[iframe]
        enames_mol = (Ts_raw[iframe] if Ts_raw is not None and iframe < len(Ts_raw) else ["?"] * len(apos_mol))
        n0_mol = n0_first
        plot_molecule(ax_geom, enames_mol, apos_mol, n0=n0_mol, title=f"Geometry @ min\nFrame {iframe}")
    else:
        ax_geom.text(0.5, 0.5, 'No geometry', transform=ax_geom.transAxes, ha='center', va='center')

def main():
    parser = argparse.ArgumentParser(description='Plot masked energy visualization')
    #parser.add_argument('--xyz', default='/home/prokop/git/FireCore-fitREQH/tests/tFitREQ/add_epairs_full/H2O-A1_HF-D1-y.xyz',help='XYZ file to process')
    #parser.add_argument('--xyz', default='/home/prokop/git/FireCore-fitREQH/tests/tFitREQ/add_epairs_full/H2O-A1_H2O-D1-y.xyz',help='XYZ file to process')
    parser.add_argument('--xyz', default='/home/prokop/git/FireCore-fitREQH/tests/tFitREQ/add_epairs_full/HCN-A1_HF-D1-y.xyz',help='XYZ file to process')
    parser.add_argument('--atypes',      default='/home/prokop/git/FireCore-fitREQH/tests/tFitREQ/data/AtomTypes.dat',help='Atom types file')
    parser.add_argument('--params',      help='JSON file with parameter overrides (overrides AtomTypes.dat and global epair settings)')
    parser.add_argument('--output',      help='Output image file')
    parser.add_argument('--kcal',        type=int, default=1, help='Convert to kcal/mol (default: 0.0 = no conversion)')
    parser.add_argument('--sym',         type=int, default=1, help='Symmetric color scale (like plot_ref.py --sym)')
    parser.add_argument('--same-scale',  type=int, default=1, help='Force same color scale for Reference vs Model Etot panels')
    parser.add_argument('--polar',       type=int, default=1, help='Use polar plotting with symmetric color scale and Rcut=6.0A')
    parser.add_argument('--rcut',        type=float, default=5.0, help='Radial cutoff for polar plots in Angstrom (default: 6.0)')
    parser.add_argument('--Rmax1D',       type=float, default=10.0, help='Maximum distance for 1D profile plots (default: 10.0)')
    parser.add_argument('--Erange',      type=float, default=None, help='Energy range for color scale (if None, use symmetric vmax=-vmin)')
    parser.add_argument('--Emin_factor', type=float, default=1.0, help='Factor to multiply min(E_ref) for vmin (default: 1.0)')
    parser.add_argument('--cmap',         default='seismic', help='Colormap for plots (default: seismic)')
    args = parser.parse_args()
    
    print(f"Loading data from {args.xyz}")
    
    # Step 1: Load reference energy grid from XYZ file using _build_frame_grid
    print("Step 1: Building reference energy grid from XYZ file...")
    V_ref, rv, A, shift, frame_idx, Ps_raw, Ts_raw = _build_frame_grid(args.xyz)
    
    # Extract geometry from first frame for polar plotting
    _, _, _, N0s = parse_xyz_with_headers(args.xyz)
    n0_first = int(N0s[0]) if N0s is not None and len(N0s) > 0 else None
    geometry_first = (Ts_raw[0], Ps_raw[0]) if Ts_raw is not None and len(Ts_raw) > 0 else None
    
    print(f"Reference energy grid shape: {V_ref.shape}")
    print(f"Reference energy range: {np.nanmin(V_ref):.6f} to {np.nanmax(V_ref):.6f}")
    
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
    
    # Extract macro using the new FittingDriver method
    morse_decomp = drv.extract_macro_from_forces('MODEL_MorseQ_PAIR_DECOMP')
    
    # Compile with the extracted macro
    drv.compile_with_model(macros={ 'MODEL_PAIR_DECOMP': morse_decomp})
    
    # Apply parameter overrides from JSON if provided
    if args.params:
        # Resolve path relative to script directory if not absolute
        params_path = args.params
        if not os.path.isabs(params_path):
            script_dir = os.path.dirname(os.path.abspath(__file__))
            params_path = os.path.join(script_dir, params_path)
        drv.apply_parameter_overrides(params_path)
    
    drv.init_and_upload_energy_only()
    
    # Print parameters used for epair calculation
    print("\n" + "="*70)
    print("PARAMETERS USED FOR EPAIR CALCULATION")
    print("="*70)
    print(f"Global epair settings:")
    print(f"  iEpairs = {drv.iEpairs} (1=exp, 2=Gaussian, 3=sigmoid)")
    print(f"  sEpairs = {drv.sEpairs} (global scaling factor)")
    print(f"\nAtom type parameters (from AtomTypes.dat):")
    for type_name in drv.atom_type_names:
        params = drv.base_params.get(type_name)
        if params:
            print(f"  {type_name:8s}: R={params['R']:.4f}, E={params['E']:.6f}, Q={params['Q']:.4f}, H={params['H']:.4f}")
            if isinstance(type_name, str) and type_name.startswith('E'):
                print(f"             ^^^ Epair type: H column used as R0ep, Q column used as Qep")
    print("="*70 + "\n")
    
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
    
    # Detect scan plane from filename (y or z)
    # y-scan varies angle in y direction → project onto xz plane to see z variation
    # z-scan varies angle in z direction → project onto xy plane to see y variation
    scan_plane = 'xz'  # default for y-scan
    if '-z.xyz' in args.xyz or '_z.xyz' in args.xyz:
        scan_plane = 'xy'
    elif '-y.xyz' in args.xyz or '_y.xyz' in args.xyz:
        scan_plane = 'xz'
    print(f"Detected scan plane: {scan_plane} (geometry projection)")
    
    if args.polar:
        # Use polar plotting with symmetric color scale and Rcut=6.0A
        # Row 1: 4 polar plots (Reference, Etot, Ein, Eout)
        # Row 2: 4 profile plots (radial slice, angular slice, E_min(angle), geometry)
        fig, axes = plt.subplots(2, 4, figsize=(16, 8), subplot_kw={'projection': 'polar'})
        
        print("Using polar plotting with symmetric color scale and Rcut=6.0A")
        
        # Calculate vmin from Reference
        vmin_ref = np.nanmin(V_ref)
        if args.kcal:
            vmin_ref *= 23.060548  # Convert eV to kcal/mol
        vmin = args.Emin_factor * vmin_ref
        
        # Calculate vmax: symmetric if Erange is None, else vmin+Erange
        if args.Erange is None:
            vmax = -vmin
            print(f"Symmetric color limits from Reference: vmin={vmin:.6f}, vmax={vmax:.6f}")
        else:
            vmax = vmin + args.Erange
            print(f"Asymmetric color limits: vmin={vmin:.6f}, vmax={vmax:.6f}")
        
        # Plot each energy component using plot_polar_symmetric with same vmin/vmax and geometry (Row 1)
        plot_polar_symmetric(V_ref, rv, A, title='Reference\n(Total)', cmap=args.cmap, kcal=args.kcal, ax=axes[0,0], bColorbar=True, rmax=args.rcut, vmin=vmin, vmax=vmax, geometry=geometry_first, n0=n0_first, plane=scan_plane)
        plot_polar_symmetric(V_model_total, rv, A, title='Model Etot\n(Total)', cmap=args.cmap, kcal=args.kcal, ax=axes[0,1], bColorbar=False, rmax=args.rcut, vmin=vmin, vmax=vmax, geometry=geometry_first, n0=n0_first, plane=scan_plane)
        plot_polar_symmetric(V_model_hbond, rv, A, title='Model Ein\n(H-bond)', cmap=args.cmap, kcal=args.kcal, ax=axes[0,2], bColorbar=False, rmax=args.rcut, vmin=vmin, vmax=vmax, geometry=geometry_first, n0=n0_first, plane=scan_plane)
        plot_polar_symmetric(V_model_eout, rv, A, title='Model Eout\n(Baseline)', cmap=args.cmap, kcal=args.kcal, ax=axes[0,3], bColorbar=False, rmax=args.rcut, vmin=vmin, vmax=vmax, geometry=geometry_first, n0=n0_first, plane=scan_plane)
        
        # Row 2: Convert to regular axes for profile plots (not polar)
        for i in range(4):
            axes[1,i].axis('off')  # Turn off polar projection
            fig.delaxes(axes[1,i])
        # Create new regular axes for row 2 and update axes array
        gs = fig.add_gridspec(2, 4)
        for i in range(4):
            axes[1, i] = fig.add_subplot(gs[1, i])
        
        # Plot second row using reusable function
        plot_profile_row(fig, axes, V_ref, V_model_total, V_model_hbond, V_model_eout, rv, A,
                        frame_idx, Ps_raw, Ts_raw, n0_first, args.kcal, args.Rmax1D)
    else:
        # Create figure with 2 rows x 4 columns for regular plotting
        fig, axes = plt.subplots(2, 4, figsize=(16, 8))
        
        # Calculate vmin from Reference
        vmin_ref = np.nanmin(V_ref)
        if args.kcal:
            vmin_ref *= 23.060548  # Convert eV to kcal/mol (to match plot_imshow's conversion)
        vmin = args.Emin_factor * vmin_ref
        
        # Calculate vmax: symmetric if Erange is None, else vmin+Erange
        if args.Erange is None:
            vmax = -vmin
            print(f"Symmetric color limits: vmin={vmin:.6f}, vmax={vmax:.6f}")
        else:
            vmax = vmin + args.Erange
            print(f"Asymmetric color limits: vmin={vmin:.6f}, vmax={vmax:.6f}")
        
        # Plot each energy component using plot_system_panel with proper vmin/vmax (Row 1)
        plot_system_panel(V_ref, rv, A, axes[0,0], 'Reference\n(Total)', args.kcal, sym=False, overlay_rmin=False, cmap=args.cmap)
        plot_system_panel(V_model_total, rv, A, axes[0,1], 'Model Etot\n(Total)', args.kcal, sym=False, overlay_rmin=False, cmap=args.cmap)
        plot_system_panel(V_model_hbond, rv, A, axes[0,2], 'Model Ein\n(H-bond)', args.kcal, sym=False, overlay_rmin=False, cmap=args.cmap)
        plot_system_panel(V_model_eout, rv, A, axes[0,3], 'Model Eout\n(Baseline)', args.kcal, sym=False, overlay_rmin=False, cmap=args.cmap)

        # Apply symmetric color limits to ALL panels
        for ax in axes[0]:
            if ax.images:
                ax.images[0].set_clim(vmin, vmax)
        
        # Plot second row using reusable function (axes[1,:] already exist)
        plot_profile_row(fig, axes, V_ref, V_model_total, V_model_hbond, V_model_eout, rv, A,
                        frame_idx, Ps_raw, Ts_raw, n0_first, args.kcal, args.Rmax1D)
    
    # Build parameter caption
    caption = f"Epairs: iEpairs={drv.iEpairs}, sEpairs={drv.sEpairs}"
    epair_type = drv.base_params.get('E')
    if epair_type:
        caption += f", E: R0ep={epair_type['H']:.3f}, Qep={epair_type['Q']:.3f}"
    else:
        caption += ", E: (not found)"
    
    plt.suptitle(f'Energy Decomposition: {os.path.basename(args.xyz)}', fontsize=14)
    plt.figtext(0.5, 0.02, caption, ha='center', fontsize=10,  bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
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
