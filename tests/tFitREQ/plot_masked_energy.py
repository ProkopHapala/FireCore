#!/usr/bin/env python3
"""
Plot masked energy decomposition using reference data and OpenCL model computation.
Properly loads reference energies, computes model with decomposition, and plots both.
"""
import argparse
import json
import os
import sys
import tempfile

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from pyBall.OCL.FittingDriver import FittingDriver
from pyBall.FitREQutils import (plot_system_panel, plot_polar_symmetric, _build_frame_grid, parse_xyz_with_headers, plot_profile_row)

def generate_dof_file(drv, output_path):
    """
    Auto-generate DOF file from current driver parameters.
    Creates reasonable bounds based on current parameter values.
    """
    comp_labels = {0:'R', 1:'E', 2:'Q', 3:'H'}
    
    # Get unique atom types present in the loaded data
    unique_types = set(drv.host_atypes)
    type_names_present = [drv.atom_type_names[i] for i in unique_types if i < len(drv.atom_type_names)]
    
    lines = []
    lines.append("# Auto-generated DOF file for MC optimization")
    lines.append("# typename component   Min         Max         xlo    xhi    Klo  Khi  K0  xstart  InvMass")
    
    for typename in type_names_present:
        if not isinstance(typename, str):
            continue
        
        # Get current parameters for this type
        params = drv.base_params.get(typename)
        if params is None:
            continue
        
        # Generate DOF entries for each component
        for comp in range(4):
            comp_name = comp_labels[comp]
            current_val = params[comp_name]
            
            # Skip if value is zero (likely not used)
            if abs(current_val) < 1e-10 and comp != 2:  # Q can be zero
                continue
            
            # Set reasonable bounds (20% around current value)
            if comp == 1:  # E parameter (always positive)
                min_val = max(1e-10, current_val * 0.5)
                max_val = current_val * 2.0
            elif comp == 2:  # Q parameter (can be negative)
                min_val = current_val - 2.0
                max_val = current_val + 2.0
            else:  # R and H
                min_val = current_val - 1.0
                max_val = current_val + 1.0
            
            line = f"{typename:12s} {comp:1d}           {min_val:10.6f}  {max_val:10.6f}  0.0    0.0    0.0  0.0  0.0  {current_val:10.6f}  1.0"
            lines.append(line)
    
    with open(output_path, 'w') as f:
        f.write('\n'.join(lines) + '\n')
    
    print(f"Generated DOF file: {output_path}")
    print(f"  Found {len([l for l in lines if not l.startswith('#')])} DOF entries")

def mc_history_to_json(history, atom_type_names):
    """
    Convert MC history (best DOFs) to JSON format compatible with epair_defaults.json.
    Maps DOF definitions back to type names and components.
    """
    json_out = {}
    comp_labels = {0:'R', 1:'E', 2:'Q', 3:'H'}
    
    for i, dof_def in enumerate(history['dof_defs']):
        typename = dof_def['typename']
        comp = dof_def['comp']
        value = history['best_dofs'][i]
        
        if typename not in json_out:
            json_out[typename] = {}
        
        json_out[typename][comp_labels[comp]] = value
    
    return json_out

def run_mc_optimization(drv, dof_file,
                        max_steps=500, step_size=0.1, temperature=0.0,
                        out_dir="mc_output"):
    """
    Run Monte Carlo optimization directly using optimizer_montecarlo.
    Returns history dict and optimized parameters.
    """
    # Import here to avoid circular dependency
    from pyBall.OCL.NonBondFitting import optimizer_montecarlo, plot_mc_convergence
    
    print(f"\n{'='*70}")
    print(f"RUNNING MC OPTIMIZATION")
    print(f"  DOFs: {dof_file}")
    print(f"  max_steps={max_steps}, step_size={step_size}, temperature={temperature}")
    print(f"{'='*70}")
    
    # Load DOF definitions
    drv.load_dofs(dof_file)
    
    # Re-initialize driver for MC optimization (energy-only path with DOFs loaded)
    drv.init_and_upload_energy_only()
    drv.setup_energy_kernel()
    
    # Get initial DOFs
    initial_dofs = np.array([d['xstart'] for d in drv.dof_definitions], dtype=np.float64)
    
    # Run Monte Carlo optimization
    history = optimizer_montecarlo(
        drv, initial_dofs,
        max_steps=max_steps, step_size=step_size,
        temperature=temperature,
        verbose=max(1, max_steps // 20)
    )
    
    # Generate convergence plot
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    conv_path = os.path.join(out_dir, "mc_convergence.png")
    fig_conv = plot_mc_convergence(history, out_path=conv_path)
    plt.close(fig_conv)
    print(f"  Saved convergence plot: {conv_path}")
    
    return history

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
    parser.add_argument('--Rmax1D',      type=float, default=10.0, help='Maximum distance for 1D profile plots (default: 10.0)')
    parser.add_argument('--Erange',      type=float, default=None, help='Energy range for color scale (if None, use symmetric vmax=-vmin)')
    parser.add_argument('--Emin_factor', type=float, default=1.0, help='Factor to multiply min(E_ref) for vmin (default: 1.0)')
    parser.add_argument('--cmap',        default='seismic', help='Colormap for plots (default: seismic)')
    parser.add_argument('--mc_before',   action='store_true', help='Run MC optimization before plotting')
    parser.add_argument('--mc_after',    action='store_true', help='Run MC optimization after plotting')
    parser.add_argument('--dofs',        help='DOF selection file for MC optimization (auto-generated if not provided)')
    parser.add_argument('--max_steps',   type=int, default=500, help='MC max steps (default: 500)')
    parser.add_argument('--step_size',   type=float, default=0.1, help='MC step size (default: 0.1)')
    parser.add_argument('--temperature', type=float, default=0.0, help='MC temperature for simulated annealing (default: 0.0 = greedy)')
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
    
    # MC Optimization: Before (if requested)
    history_before = None
    if args.mc_before:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        
        # Generate DOF file if not provided
        if args.dofs is None:
            dof_file = os.path.join(script_dir, "auto_dof.dat")
            generate_dof_file(drv, dof_file)
        else:
            dof_file = args.dofs
            if not os.path.isabs(dof_file):
                dof_file = os.path.join(script_dir, dof_file)
        
        # Run MC optimization
        history_before = run_mc_optimization(
            drv, dof_file,
            max_steps=args.max_steps, step_size=args.step_size,
            temperature=args.temperature, out_dir="mc_before_output"
        )
        
        # Convert MC results to JSON and apply
        optimized_json = mc_history_to_json(history_before, drv.atom_type_names)
        
        # Save optimized parameters to temporary JSON file
        temp_params_path = os.path.join(script_dir, "mc_optimized_before.json")
        with open(temp_params_path, 'w') as f:
            json.dump(optimized_json, f, indent=2)
        
        # Apply optimized parameters
        drv.apply_parameter_overrides(temp_params_path)
        drv.init_and_upload_energy_only()  # Re-upload with new parameters
        
        print(f"\nMC Before Results:")
        print(f"  Initial error: {history_before['initial_error']:.6e}")
        print(f"  Best error:    {history_before['best_error']:.6e}")
        print(f"  Steps accepted: {history_before['n_accepted']} / {history_before['n_steps']}")
    
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
    
    # MC Optimization: After (if requested)
    if args.mc_after:
        print(f"\n{'='*70}")
        print(f"RUNNING MC OPTIMIZATION AFTER PLOTTING")
        print(f"{'='*70}")
        
        script_dir = os.path.dirname(os.path.abspath(__file__))
        
        # Generate DOF file if not provided
        if args.dofs is None:
            dof_file = os.path.join(script_dir, "auto_dof.dat")
            generate_dof_file(drv, dof_file)
        else:
            dof_file = args.dofs
            if not os.path.isabs(dof_file):
                dof_file = os.path.join(script_dir, dof_file)
        
        # Run MC optimization
        history_after = run_mc_optimization(
            drv, dof_file,
            max_steps=args.max_steps, step_size=args.step_size,
            temperature=args.temperature, out_dir="mc_after_output"
        )
        
        # Convert MC results to JSON and apply
        optimized_json = mc_history_to_json(history_after, drv.atom_type_names)
        
        # Save optimized parameters to temporary JSON file
        temp_params_path = os.path.join(script_dir, "mc_optimized_after.json")
        with open(temp_params_path, 'w') as f:
            json.dump(optimized_json, f, indent=2)
        
        # Apply optimized parameters
        drv.apply_parameter_overrides(temp_params_path)
        drv.init_and_upload_energy_only()
        
        # Re-compute energies with optimized parameters
        _, Em_all_opt      = drv.evaluate_energies_masked(mask_all)
        _, Em_hbond_opt    = drv.evaluate_energies_masked(mask_hbond)
        _, Em_baseline_opt = drv.evaluate_energies_masked(mask_baseline)
        
        # Reshape optimized energies to match grid
        V_model_total_opt    = np.full((ny, nx), np.nan, dtype=np.float32)
        V_model_hbond_opt    = np.full((ny, nx), np.nan, dtype=np.float32)
        V_model_baseline_opt = np.full((ny, nx), np.nan, dtype=np.float32)
        
        for iy in range(ny):
            for ix in range(nx):
                isamp = int(frame_idx[iy, ix])
                if isamp >= 0 and isamp < n_samples:
                    V_model_total_opt[iy, ix]    = Em_all_opt[isamp]
                    V_model_hbond_opt[iy, ix]    = Em_hbond_opt[isamp]
                    V_model_baseline_opt[iy, ix] = Em_baseline_opt[isamp]
        
        V_model_eout_opt = V_model_total_opt - V_model_hbond_opt
        
        print(f"\nMC After Results:")
        print(f"  Initial error: {history_after['initial_error']:.6e}")
        print(f"  Best error:    {history_after['best_error']:.6e}")
        print(f"  Steps accepted: {history_after['n_accepted']} / {history_after['n_steps']}")
        
        # Create before/after comparison figure
        print("\nCreating before/after comparison plot...")
        fig_comp, axes_comp = plt.subplots(2, 4, figsize=(16, 8))
        
        # Row 1: Before optimization
        plot_system_panel(V_ref, rv, A, axes_comp[0,0], 'Reference\n(Total)', args.kcal, sym=False, overlay_rmin=False, cmap=args.cmap)
        plot_system_panel(V_model_total, rv, A, axes_comp[0,1], 'Model Etot\n(Before)', args.kcal, sym=False, overlay_rmin=False, cmap=args.cmap)
        plot_system_panel(V_model_hbond, rv, A, axes_comp[0,2], 'Model Ein\n(Before)', args.kcal, sym=False, overlay_rmin=False, cmap=args.cmap)
        plot_system_panel(V_model_eout, rv, A, axes_comp[0,3], 'Model Eout\n(Before)', args.kcal, sym=False, overlay_rmin=False, cmap=args.cmap)
        
        # Row 2: After optimization
        plot_system_panel(V_ref, rv, A, axes_comp[1,0], 'Reference\n(Total)', args.kcal, sym=False, overlay_rmin=False, cmap=args.cmap)
        plot_system_panel(V_model_total_opt, rv, A, axes_comp[1,1], 'Model Etot\n(After)', args.kcal, sym=False, overlay_rmin=False, cmap=args.cmap)
        plot_system_panel(V_model_hbond_opt, rv, A, axes_comp[1,2], 'Model Ein\n(After)', args.kcal, sym=False, overlay_rmin=False, cmap=args.cmap)
        plot_system_panel(V_model_eout_opt, rv, A, axes_comp[1,3], 'Model Eout\n(After)', args.kcal, sym=False, overlay_rmin=False, cmap=args.cmap)
        
        # Apply symmetric color limits to all panels
        for ax_row in axes_comp:
            for ax in ax_row:
                if ax.images:
                    ax.images[0].set_clim(vmin, vmax)
        
        # Build parameter caption
        caption_comp = f"MC Optimization: {args.max_steps} steps, step_size={args.step_size}, T={args.temperature}"
        plt.suptitle(f'Before/After Comparison: {os.path.basename(args.xyz)}', fontsize=14)
        plt.figtext(0.5, 0.02, caption_comp, ha='center', fontsize=10, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        plt.tight_layout()
        
        if args.output:
            comp_output = args.output.replace('.png', '_comparison.png')
            plt.savefig(comp_output, dpi=200, bbox_inches='tight')
            print(f"Saved comparison plot: {comp_output}")
        else:
            plt.show()
        
        # Print statistics for optimized parameters
        print("\nOptimized Energy Statistics:")
        print(f"  Model Total (After): min={np.nanmin(V_model_total_opt):.6f}, max={np.nanmax(V_model_total_opt):.6f}, mean={np.nanmean(V_model_total_opt):.6f}")
        print(f"  Model H-bond (After): min={np.nanmin(V_model_hbond_opt):.6f}, max={np.nanmax(V_model_hbond_opt):.6f}, mean={np.nanmean(V_model_hbond_opt):.6f}")
        print(f"  Model Eout (After): min={np.nanmin(V_model_eout_opt):.6f}, max={np.nanmax(V_model_eout_opt):.6f}, mean={np.nanmean(V_model_eout_opt):.6f}")
        
        # Compute correlation for optimized parameters
        valid_mask_opt = np.isfinite(V_ref) & np.isfinite(V_model_total_opt)
        if np.any(valid_mask_opt):
            correlation_opt = np.corrcoef(V_ref[valid_mask_opt], V_model_total_opt[valid_mask_opt])[0, 1]
            print(f"\nReference-Model correlation (After): {correlation_opt:.6f}")
            print(f"Correlation improvement: {correlation_opt - correlation:.6f}")

if __name__ == "__main__":
    main()
