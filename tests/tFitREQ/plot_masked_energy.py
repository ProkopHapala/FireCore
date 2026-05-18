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

def create_dof_definitions_from_current(drv, config=None):
    """
    Create DOF definitions programmatically from current driver state (tREQHs_base).
    No DOF file needed - uses current parameters (including JSON overrides) as starting point.
    If config is provided, only optimize parameters specified in the config.
    Config file format: {"typename": ["R", "E", "Q", "H"]} or {"typename": {"R": true, "E": false, ...}}
    """
    comp_labels = {0:'R', 1:'E', 2:'Q', 3:'H'}
    comp_indices_map = {'R':0, 'E':1, 'Q':2, 'H':3}
    
    # Get unique atom types present in the loaded data
    unique_types = set(drv.host_atypes)
    type_names_present = [drv.atom_type_names[i] for i in unique_types if i < len(drv.atom_type_names)]
    
    dof_definitions = []
    
    for typename in type_names_present:
        if not isinstance(typename, str):
            continue
        
        # Check if this typename is in config (if config provided)
        if config and typename not in config:
            continue
        
        # Get which components to optimize for this typename
        if config and typename in config:
            opt_comps = config[typename]
            if isinstance(opt_comps, list):
                # Convert component names to indices (reverse lookup)
                comp_labels_reverse = {v: k for k, v in comp_labels.items()}
                comp_indices = [comp_labels_reverse.get(c, -1) for c in opt_comps]
                comp_indices = [c for c in comp_indices if c >= 0]
            elif isinstance(opt_comps, dict):
                # Dict with component names as keys and bool values
                comp_labels_reverse = {v: k for k, v in comp_labels.items()}
                comp_indices = [comp_labels_reverse.get(c, -1) for c, val in opt_comps.items() if val]
                comp_indices = [c for c in comp_indices if c >= 0]
            else:
                comp_indices = list(range(4))  # All components
        else:
            comp_indices = list(range(4))  # All components
        
        # Get type index
        try:
            type_idx = drv.atom_type_names.index(typename)
        except ValueError:
            continue
        
        # Get current parameters from tREQHs_base (includes JSON overrides)
        # tREQHs_base is [n_types, 4] with columns: R, sqrt(E), Q, H
        if not hasattr(drv, 'tREQHs_base') or type_idx >= drv.tREQHs_base.shape[0]:
            continue
        
        treqh = drv.tREQHs_base[type_idx]  # [R, sqrt(E), Q, H]
        
        # For epair types, E is stored directly (not sqrt)
        is_epair = isinstance(typename, str) and typename.startswith('E')
        
        # Generate DOF entries for each component
        for comp in comp_indices:
            comp_name = comp_labels[comp]
            
            if comp == 0:  # R
                current_val = treqh[0]
            elif comp == 1:  # E (need to handle sqrt)
                current_val = treqh[1] if is_epair else (treqh[1] ** 2)
            elif comp == 2:  # Q
                current_val = treqh[2]
            else:  # H
                current_val = treqh[3]
            
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
            elif comp == 0:  # R parameter (always positive)
                min_val = max(0.1, current_val * 0.5)
                max_val = current_val * 2.0
            else:  # H parameter (homogeneous correction - can be positive or negative)
                min_val = current_val - 2.0
                max_val = current_val + 2.0
            
            dof_definitions.append({
                'typename': typename,
                'comp': comp,
                'min': min_val,
                'max': max_val,
                'xlo': 0.0,
                'xhi': 0.0,
                'Klo': 0.0,
                'Khi': 0.0,
                'K0': 0.0,
                'xstart': current_val
            })
            print(f"  DOF: {typename}.{comp_labels[comp]} = {current_val:.6f}  range=[{min_val:.6f}, {max_val:.6f}]")
    
    print(f"Created {len(dof_definitions)} DOF definitions from current parameters")
    return dof_definitions

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

def load_mc_config_with_softclamp(config_file):
    """
    Load MC config file and extract soft-clamp parameters.
    Returns (config_dict, soft_clamp, clamp_start, clamp_max)
    """
    config = None
    soft_clamp = False
    clamp_start = 4.0
    clamp_max = 6.0
    
    if config_file:
        if not os.path.isabs(config_file):
            script_dir = os.path.dirname(os.path.abspath(__file__))
            config_file = os.path.join(script_dir, config_file)
        if os.path.exists(config_file):
            with open(config_file, 'r') as f:
                config = json.load(f)
            print(f"Loaded MC config from {config_file}")
            # Extract soft-clamp parameters if present
            soft_clamp = config.pop('soft_clamp', False)
            clamp_start = config.pop('clamp_start', 4.0)
            clamp_max = config.pop('clamp_max', 6.0)
        else:
            print(f"Warning: MC config file not found: {config_file}, optimizing all parameters")
    
    return config, soft_clamp, clamp_start, clamp_max

def _map_samples_to_grid(frame_idx, V_ref_shape, values_1d, fill=np.nan):
    ny, nx = V_ref_shape
    V = np.full((ny, nx), fill, dtype=np.float32)
    n = len(values_1d)
    for iy in range(ny):
        for ix in range(nx):
            isamp = int(frame_idx[iy, ix])
            if isamp >= 0 and isamp < n:
                V[iy, ix] = values_1d[isamp]
    return V

def _set_im_clim(im, v0, v1, kcal=False):
    fac = 23.060548 if kcal else 1.0
    im.set_clim(v0 * fac, v1 * fac)

def run_mc_optimization(drv,
                        max_steps=500, step_size=0.1, temperature=0.0,
                        soft_clamp=False, clamp_start=4.0, clamp_max=6.0,
                        kcal_objective=False,
                        out_dir="mc_output", config_file=None):
    """
    Run Monte Carlo optimization directly using optimizer_montecarlo.
    Creates DOF definitions programmatically from current driver state (no DOF file needed).
    If config_file is provided, only optimizes parameters specified in the config.
    Returns history dict with optimized parameters.
    """
    # Import here to avoid circular dependency
    from pyBall.OCL.NonBondFitting import optimizer_montecarlo, plot_mc_convergence
    
    print(f"\n{'='*70}")
    print(f"RUNNING MC OPTIMIZATION")
    print(f"  max_steps={max_steps}, step_size={step_size}, temperature={temperature}")
    if soft_clamp:
        print(f"  soft_clamp: start={clamp_start}, max={clamp_max}")
    print(f"{'='*70}")
    
    # Create DOF definitions programmatically from current parameters
    dof_definitions = create_dof_definitions_from_current(drv, config=config_file)
    
    # Set DOF definitions directly on driver (no file needed)
    drv.dof_definitions = dof_definitions
    drv.n_dofs = len(dof_definitions)
    
    # Setup energy kernel (don't re-init, just setup kernel for objective evaluation)
    drv.setup_energy_kernel()
    
    # Get initial DOFs from current driver state (tREQHs_base includes JSON overrides)
    initial_dofs = np.array([d['xstart'] for d in dof_definitions], dtype=np.float64)
    
    print(f"  Initial DOFs: {initial_dofs}")
    
    # Run Monte Carlo optimization
    history = optimizer_montecarlo(
        drv, initial_dofs,
        max_steps=max_steps, step_size=step_size,
        temperature=temperature,
        soft_clamp=soft_clamp, clamp_start=clamp_start, clamp_max=clamp_max,
        kcal_objective=kcal_objective,
        verbose=10
    )
    
    # Store DOF definitions in history for later conversion
    history['dof_defs'] = dof_definitions
    
    # Generate convergence plot
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    conv_path = os.path.join(out_dir, "mc_convergence.png")
    fig_conv = plot_mc_convergence(history, out_path=conv_path)
    plt.close(fig_conv)
    print(f"  Convergence plot: {os.path.abspath(conv_path)}")
    
    return history


'''
HOW TO RUN (Examples):

python3 plot_masked_energy.py --xyz add_epairs_full/H2O-A1_H2O-D1-y.xyz --params epair_defaults.json --polar 0 --output mc_dbg22.png --mc_after --mc_config mc_config.json --max_steps 20 --step_size 0.2

'''


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
    parser.add_argument('--mc_config',   help='JSON config file specifying which parameters to optimize (e.g., mc_config.json)')
    parser.add_argument('--max_steps',   type=int, default=500, help='MC max steps (default: 500)')
    parser.add_argument('--step_size',   type=float, default=0.1, help='MC step size (default: 0.1)')
    parser.add_argument('--temperature', type=float, default=0.0, help='MC temperature for simulated annealing (default: 0.0 = greedy)')
    parser.add_argument('--soft_clamp',  action='store_true', help='Enable soft-clamping of energy differences before squaring')
    parser.add_argument('--clamp_start', type=float, default=4.0, help='Soft-clamp start threshold (default: 4.0)')
    parser.add_argument('--clamp_max',   type=float, default=6.0, help='Soft-clamp maximum threshold (default: 6.0)')
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
    
    # Apply parameter overrides from JSON BEFORE compilation
    # This ensures tREQHs_base is updated before DOF definitions are created
    if args.params:
        # Resolve path relative to script directory if not absolute
        params_path = args.params
        if not os.path.isabs(params_path):
            script_dir = os.path.dirname(os.path.abspath(__file__))
            params_path = os.path.join(script_dir, params_path)
        drv.apply_parameter_overrides(params_path)
    
    # Step 3: Compile with MODEL_MorseQ_PAIR_DECOMP macro
    print("Step 3: Compiling with Morse decomposition macro...")
    
    # Extract macro using the new FittingDriver method
    morse_decomp = drv.extract_macro_from_forces('MODEL_MorseQ_PAIR_DECOMP')
    
    # CRITICAL FIX: Disable HBOND_GATE to allow homogeneous correction (H-component) to work
    # HBOND_GATE gates H to zero when positive, which breaks optimization of O_3.H (homogeneous correction)
    drv.compile_with_model(macros={
        'MODEL_PAIR_DECOMP': morse_decomp,
        'HBOND_GATE_DEFINE': '#define HBOND_GATE 0'
    })
    
    drv.init_and_upload_energy_only()

    # IMPORTANT: make the parameter state explicit for subsequent energy evaluations.
    # init_and_upload_energy_only() does not guarantee that the current `tREQHs_base` is uploaded
    # into `tREQHs_buff`. If we don't upload it here, the first plot can use stale GPU parameters
    # from a previous call (which is exactly the kind of inconsistency we're hunting).
    assert hasattr(drv, 'tREQHs_base')
    drv.toGPU_(drv.tREQHs_buff, np.array(drv.tREQHs_base, dtype=np.float32, copy=False))
    
    # Load MC config file (if provided) to extract soft-clamp parameters
    mc_config_dict, mc_soft_clamp, mc_clamp_start, mc_clamp_max = load_mc_config_with_softclamp(args.mc_config)
    
    # MC Optimization: Before (if requested)
    history_before = None
    if args.mc_before:
        # Use CLI args if provided, otherwise use config values
        soft_clamp_use = args.soft_clamp if args.soft_clamp else mc_soft_clamp
        clamp_start_use = args.clamp_start if args.clamp_start != 4.0 else mc_clamp_start
        clamp_max_use = args.clamp_max if args.clamp_max != 6.0 else mc_clamp_max
        
        # Run MC optimization (DOF definitions created programmatically from current parameters)
        history_before = run_mc_optimization(
            drv,
            max_steps=args.max_steps, step_size=args.step_size,
            temperature=args.temperature,
            soft_clamp=soft_clamp_use, clamp_start=clamp_start_use, clamp_max=clamp_max_use,
            out_dir="mc_before_output",
            config_file=mc_config_dict
        )
        
        # Convert MC results to JSON and apply
        optimized_json = mc_history_to_json(history_before, drv.atom_type_names)
        
        # Save optimized parameters to temporary JSON file
        temp_params_path = os.path.join(script_dir, "mc_optimized_before.json")
        temp_params_path_abs = os.path.abspath(temp_params_path)
        with open(temp_params_path_abs, 'w') as f:
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

    # If we are going to run MC after plotting, we must ensure the *initial* plot uses
    # exactly the same parameter upload mechanism as MC (DOF-derived tREQHs).
    # Otherwise `tREQHs_base` vs `tREQHs_from_dofs(initial_dofs)` can disagree and the plots diverge.
    if args.mc_after:
        dof_definitions_plot = create_dof_definitions_from_current(drv, config=mc_config_dict)
        drv.dof_definitions = dof_definitions_plot
        drv.n_dofs = len(dof_definitions_plot)
        initial_dofs_plot = np.array([d['xstart'] for d in dof_definitions_plot], dtype=np.float64)
        T_plot = drv.tREQHs_from_dofs(initial_dofs_plot)
        drv.toGPU_(drv.tREQHs_buff, T_plot)

    # Step 4: Compute model energies for all samples on GPU (parallel)
    print("Step 4: Computing model energies for all samples (GPU parallel)...")

    # Number of samples is defined by loaded dataset
    n_samples = int(drv.n_samples)
    print(f"n_samples = {n_samples}")

    # Create masks (component weights)
    mask_hbond    = drv.create_mask(ni, nj, pauli=0.0, london=0.0, electro=0.0, hbond=1.0)

    # IMPORTANT: evaluate_energies_masked returns (Etot, Emasked).
    # For decomposition we must take Etot from the FIRST return value,
    # and Ein from the masked value for the selected mask.
    # Compute energies ONCE and reuse for all plots (avoid redundant kernel calls).
    Em_tot, Em_hbond = drv.evaluate_energies_masked(mask_hbond)

    assert Em_tot.shape[0] == n_samples
    assert Em_hbond.shape[0] == n_samples

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
                V_model_total[iy, ix]    = Em_tot[isamp]
                V_model_hbond[iy, ix]    = Em_hbond[isamp]
                V_model_baseline[iy, ix] = Em_tot[isamp] - Em_hbond[isamp]
    
    # Compute model Eout = total - hbond
    V_model_eout = V_model_total - V_model_hbond
    
    # Step 7: Plot reference and model energies with proper vmin/vmax
    print("Step 7: Plotting energy components...")

    def _save_main_plot(V_ref_, V_model_total_, V_model_hbond_, V_model_eout_, output_path_abs):
        # Detect scan plane from filename (y or z)
        scan_plane = 'xz'  # default for y-scan
        if '-z.xyz' in args.xyz or '_z.xyz' in args.xyz:
            scan_plane = 'xy'
        elif '-y.xyz' in args.xyz or '_y.xyz' in args.xyz:
            scan_plane = 'xz'
        print(f"Detected scan plane: {scan_plane} (geometry projection)")

        if args.polar:
            fig, axes = plt.subplots(2, 4, figsize=(16, 8), subplot_kw={'projection': 'polar'})

            print("Using polar plotting with symmetric color scale and Rcut=6.0A")

            vmin_ref = np.nanmin(V_ref_)
            if args.kcal:
                vmin_ref *= 23.060548
            vmin = args.Emin_factor * vmin_ref

            if args.Erange is None:
                vmax = -vmin
                print(f"Symmetric color limits from Reference: vmin={vmin:.6f}, vmax={vmax:.6f}")
            else:
                vmax = vmin + args.Erange
                print(f"Asymmetric color limits: vmin={vmin:.6f}, vmax={vmax:.6f}")

            plot_polar_symmetric(V_ref_, rv, A, title='Reference\n(Total)', cmap=args.cmap, kcal=args.kcal, ax=axes[0,0], bColorbar=True, rmax=args.rcut, vmin=vmin, vmax=vmax, geometry=geometry_first, n0=n0_first, plane=scan_plane)
            plot_polar_symmetric(V_model_total_, rv, A, title='Model Etot\n(Total)', cmap=args.cmap, kcal=args.kcal, ax=axes[0,1], bColorbar=False, rmax=args.rcut, vmin=vmin, vmax=vmax, geometry=geometry_first, n0=n0_first, plane=scan_plane)
            plot_polar_symmetric(V_model_hbond_, rv, A, title='Model Ein\n(H-bond)', cmap=args.cmap, kcal=args.kcal, ax=axes[0,2], bColorbar=False, rmax=args.rcut, vmin=vmin, vmax=vmax, geometry=geometry_first, n0=n0_first, plane=scan_plane)
            plot_polar_symmetric(V_model_eout_, rv, A, title='Model Eout\n(Baseline)', cmap=args.cmap, kcal=args.kcal, ax=axes[0,3], bColorbar=False, rmax=args.rcut, vmin=vmin, vmax=vmax, geometry=geometry_first, n0=n0_first, plane=scan_plane)

            for i in range(4):
                axes[1,i].axis('off')
                fig.delaxes(axes[1,i])
            gs = fig.add_gridspec(2, 4)
            for i in range(4):
                axes[1, i] = fig.add_subplot(gs[1, i])

            plot_profile_row(fig, axes, V_ref_, V_model_total_, V_model_hbond_, V_model_eout_, rv, A, frame_idx, Ps_raw, Ts_raw, n0_first, args.kcal, args.Rmax1D)
        else:
            fig, axes = plt.subplots(2, 4, figsize=(16, 8))

            vmin_ref = np.nanmin(V_ref_)
            if args.kcal:
                vmin_ref *= 23.060548
            vmin = args.Emin_factor * vmin_ref

            if args.Erange is None:
                vmax = -vmin
                print(f"Symmetric color limits: vmin={vmin:.6f}, vmax={vmax:.6f}")
            else:
                vmax = vmin + args.Erange
                print(f"Asymmetric color limits: vmin={vmin:.6f}, vmax={vmax:.6f}")

            plot_system_panel(V_ref_, rv, A, axes[0,0], 'Reference\n(Total)', args.kcal, sym=False, overlay_rmin=False, cmap=args.cmap)
            plot_system_panel(V_model_total_, rv, A, axes[0,1], 'Model Etot\n(Total)', args.kcal, sym=False, overlay_rmin=False, cmap=args.cmap)
            plot_system_panel(V_model_hbond_, rv, A, axes[0,2], 'Model Ein\n(H-bond)', args.kcal, sym=False, overlay_rmin=False, cmap=args.cmap)
            plot_system_panel(V_model_eout_, rv, A, axes[0,3], 'Model Eout\n(Baseline)', args.kcal, sym=False, overlay_rmin=False, cmap=args.cmap)

            for ax in axes[0]:
                if ax.images:
                    ax.images[0].set_clim(vmin, vmax)

            plot_profile_row(fig, axes, V_ref_, V_model_total_, V_model_hbond_, V_model_eout_, rv, A, frame_idx, Ps_raw, Ts_raw, n0_first, args.kcal, args.Rmax1D)

        caption = f"Epairs: iEpairs={drv.iEpairs}, sEpairs={drv.sEpairs}"
        epair_type = drv.base_params.get('E')
        if epair_type:
            caption += f", E: R0ep={epair_type['H']:.3f}, Qep={epair_type['Q']:.3f}"
        else:
            caption += ", E: (not found)"

        plt.suptitle(f'Energy Decomposition: {os.path.basename(args.xyz)}', fontsize=14)
        plt.figtext(0.5, 0.02, caption, ha='center', fontsize=10,  bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        plt.tight_layout()
        plt.savefig(output_path_abs, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f"Output file: {output_path_abs}")
    
    output_path = os.path.abspath(args.output)
    _save_main_plot(V_ref, V_model_total, V_model_hbond, V_model_eout, output_path)
    
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
        
        # Use CLI args if provided, otherwise use config values
        soft_clamp_use = args.soft_clamp if args.soft_clamp else mc_soft_clamp
        clamp_start_use = args.clamp_start if args.clamp_start != 4.0 else mc_clamp_start
        clamp_max_use = args.clamp_max if args.clamp_max != 6.0 else mc_clamp_max
        
        # Run MC optimization (DOF definitions created programmatically from current parameters)
        history_after = run_mc_optimization(drv,
                                             max_steps=args.max_steps,
                                             step_size=args.step_size,
                                             temperature=args.temperature,
                                             soft_clamp=soft_clamp_use,
                                             clamp_start=clamp_start_use,
                                             clamp_max=clamp_max_use,
                                             kcal_objective=args.kcal,
                                             config_file=mc_config_dict)
        
        # Convert MC results to JSON and apply
        optimized_json = mc_history_to_json(history_after, drv.atom_type_names)
        
        # Save optimized parameters to temporary JSON file
        temp_params_path = os.path.join(script_dir, "mc_optimized_after.json")
        temp_params_path_abs = os.path.abspath(temp_params_path)
        with open(temp_params_path_abs, 'w') as f:
            json.dump(optimized_json, f, indent=2)
        
        # RIGOROUS decomposition from a SINGLE kernel call per state:
        # evaluate_energies_masked(mask_hbond) returns (Etot, Emasked) where Emasked is our Ein correction.
        # Therefore Eout must be Etot - Ein and must not change when only epair params (inside mask) change.
        params_initial_path = args.params
        assert params_initial_path is not None

        # We must use exactly the same tREQHs used by the MC optimizer (derived from DOFs),
        # otherwise we can end up comparing different parameterizations even if the JSON looks the same.
        dofs0 = history_after.get('initial_dofs', None)
        dofs1 = history_after.get('best_dofs', None)
        assert (dofs0 is not None) and (dofs1 is not None)

        # Set base parameters (tREQHs_base) from initial JSON once
        drv.apply_parameter_overrides(params_initial_path)
        drv.init_and_upload_energy_only()

        # --- Initial (Before): upload tREQHs from DOFs and evaluate (single kernel call) ---
        T0 = drv.tREQHs_from_dofs(dofs0)
        drv.toGPU_(drv.tREQHs_buff, T0)
        Etot0, Ein0 = drv.evaluate_energies_masked(mask_hbond)
        Eout0 = Etot0 - Ein0

        # --- Optimized (After): upload tREQHs from best DOFs and evaluate ---
        T1 = drv.tREQHs_from_dofs(dofs1)
        drv.toGPU_(drv.tREQHs_buff, T1)
        Etot1, Ein1 = drv.evaluate_energies_masked(mask_hbond)
        Eout1 = Etot1 - Ein1

        # RIGOROUS CHECK: the main plot energies must match the MC-initial (T0) energies.
        # If they don't, it means mc_dbgXX.png was computed from a different kernel state.
        V_model_total_mc0 = _map_samples_to_grid(frame_idx, V_ref.shape, Etot0)
        V_model_hbond_mc0 = _map_samples_to_grid(frame_idx, V_ref.shape, Ein0)
        V_model_eout_mc0  = V_model_total_mc0 - V_model_hbond_mc0

        V_model_total_mc1 = _map_samples_to_grid(frame_idx, V_ref.shape, Etot1)
        V_model_hbond_mc1 = _map_samples_to_grid(frame_idx, V_ref.shape, Ein1)
        V_model_eout_mc1  = V_model_total_mc1 - V_model_hbond_mc1
        m = np.isfinite(V_model_total) & np.isfinite(V_model_total_mc0)
        if np.any(m):
            d_tot = (V_model_total[m] - V_model_total_mc0[m]).astype(np.float64)
            d_in  = (V_model_hbond[m] - V_model_hbond_mc0[m]).astype(np.float64)
            d_out = (V_model_eout[m]  - V_model_eout_mc0[m]).astype(np.float64)
            print("MC DIAG: RIGOROUS CHECK - main plot (single) vs MC-initial (T0) grids")
            print(f"  max|Etot_plot - Etot_T0| = {float(np.max(np.abs(d_tot))):.6e}   mean={float(np.mean(d_tot)):.6e}")
            print(f"  max|Ein_plot  - Ein_T0 | = {float(np.max(np.abs(d_in ))):.6e}   mean={float(np.mean(d_in )):.6e}")
            print(f"  max|Eout_plot - Eout_T0| = {float(np.max(np.abs(d_out))):.6e}   mean={float(np.mean(d_out)):.6e}")
            if (np.max(np.abs(d_tot)) > 1e-6) or (np.max(np.abs(d_in)) > 1e-6) or (np.max(np.abs(d_out)) > 1e-6):
                raise RuntimeError("BUG: mc_dbgXX.png energies differ from MC-initial (T0) energies; plot was computed from different state")

        # Save before/after main plots with 1D cuts + geometry for easier inspection.
        # Keep the original output name as an alias for "before" (historical behavior).
        out_root, out_ext = os.path.splitext(output_path)
        out_before = out_root + "_before" + out_ext
        out_after  = out_root + "_after"  + out_ext

        print(f"Saving main plot (MC-initial/T0): {out_before}")
        _save_main_plot(V_ref, V_model_total_mc0, V_model_hbond_mc0, V_model_eout_mc0, out_before)

        print(f"Saving main plot (MC-best/T1): {out_after}")
        _save_main_plot(V_ref, V_model_total_mc1, V_model_hbond_mc1, V_model_eout_mc1, out_after)

        print("Overwriting main output plot using authoritative MC-initial energies (T0)...")
        _save_main_plot(V_ref, V_model_total_mc0, V_model_hbond_mc0, V_model_eout_mc0, output_path)

        # Rigorous check: Eout must be invariant (this is your key test)
        dEout = Eout1 - Eout0
        dEout_max = float(np.max(np.abs(dEout)))
        dEout_mean = float(np.mean(dEout))
        print("MC DIAG: RIGOROUS CHECK - Eout invariance (Eout = Etot - Ein(masked))")
        print(f"  max|Eout_after - Eout_before| = {dEout_max:.6e}   mean={dEout_mean:.6e}")
        # float32 tolerance: this should be ~1e-6..1e-5 if truly invariant
        if dEout_max > 1e-4:
            print("  WARNING: Eout changed => either mask does NOT capture all epair effects or epair params leak into baseline")
        else:
            print("  OK: Eout invariant within tolerance")

        # Build panels for comparison plot from these authoritative arrays
        V_model_total = np.full((ny, nx), np.nan, dtype=np.float32)
        V_model_hbond = np.full((ny, nx), np.nan, dtype=np.float32)   # Ein
        V_model_eout  = np.full((ny, nx), np.nan, dtype=np.float32)   # baseline remainder

        V_model_total_opt = np.full((ny, nx), np.nan, dtype=np.float32)
        V_model_hbond_opt = np.full((ny, nx), np.nan, dtype=np.float32)
        V_model_eout_opt  = np.full((ny, nx), np.nan, dtype=np.float32)

        for iy in range(ny):
            for ix in range(nx):
                isamp = int(frame_idx[iy, ix])
                if isamp >= 0 and isamp < n_samples:
                    V_model_total[iy, ix]     = Etot0[isamp]
                    V_model_hbond[iy, ix]     = Ein0[isamp]
                    V_model_eout[iy, ix]      = Eout0[isamp]
                    V_model_total_opt[iy, ix] = Etot1[isamp]
                    V_model_hbond_opt[iy, ix] = Ein1[isamp]
                    V_model_eout_opt[iy, ix]  = Eout1[isamp]

        print(f"\nMC After Results:")
        print(f"  Initial error: {history_after['initial_error']:.6e}")
        print(f"  Best error:    {history_after['best_error']:.6e}")
        print(f"  Steps accepted: {history_after['n_accepted']} / {history_after['n_steps']}")

        # --- RIGOROUS DIAGNOSTICS: verify what optimizer actually used ---
        W = history_after.get('W', None)
        Eref = history_after.get('Eref', None)
        dE0 = history_after.get('dE_initial', None)
        dE1 = history_after.get('dE_best', None)
        Jps0 = history_after.get('J_per_sample_initial', None)
        Jps1 = history_after.get('J_per_sample_best', None)
        assert (W is not None) and (Eref is not None) and (dE0 is not None) and (dE1 is not None) and (Jps0 is not None) and (Jps1 is not None)

        print("\nMC DIAG: weights (W) used by objective")
        print(f"  W: dtype={W.dtype} shape={W.shape} min={float(np.min(W)):.6e} max={float(np.max(W)):.6e} unique~={len(np.unique(W))}")

        # Consistency: total J must equal sum of per-sample contributions
        sJ0 = float(np.sum(Jps0))
        sJ1 = float(np.sum(Jps1))
        print("MC DIAG: objective consistency")
        print(f"  sum(J_per_sample_initial)={sJ0:.9e}  history.initial_error={float(history_after['initial_error']):.9e}  diff={sJ0-float(history_after['initial_error']):.3e}")
        print(f"  sum(J_per_sample_best)   ={sJ1:.9e}  history.best_error   ={float(history_after['best_error']):.9e}  diff={sJ1-float(history_after['best_error']):.3e}")

        # Get soft-clamp and kcal_objective flags for formula check and diagnostics
        kcal_obj = history_after.get('kcal_objective', False)
        sc_enabled = history_after.get('soft_clamp', False)
        sc_start = float(history_after.get('clamp_start', 4.0))
        sc_max = float(history_after.get('clamp_max', 6.0))

        # Consistency: per-sample formula check
        # NOTE: dE0 and dE1 are raw (unclamped) eV values, but J is computed with soft-clamp in kcal
        # Apply kcal conversion, then soft-clamp to dE for formula check to match what objective actually used
        if sc_enabled:
            # Convert to kcal first (clamp parameters are in kcal when kcal_objective=True)
            dE0c = dE0 * (23.060548 if kcal_obj else 1.0)
            dE1c = dE1 * (23.060548 if kcal_obj else 1.0)
            # Then apply soft-clamp
            dE0c = drv._soft_clamp(dE0c, sc_start, sc_max)
            dE1c = drv._soft_clamp(dE1c, sc_start, sc_max)
        else:
            dE0c = dE0 * (23.060548 if kcal_obj else 1.0)
            dE1c = dE1 * (23.060548 if kcal_obj else 1.0)

        Jps0_chk = 0.5 * W * dE0c * dE0c
        Jps1_chk = 0.5 * W * dE1c * dE1c
        print("MC DIAG: per-sample formula check")
        print(f"  max|Jps0 - 0.5*W*dE_clamped^2| = {float(np.max(np.abs(Jps0 - Jps0_chk))):.3e}")
        print(f"  max|Jps1 - 0.5*W*dE_clamped^2| = {float(np.max(np.abs(Jps1 - Jps1_chk))):.3e}")
        print(f"  dE_initial (raw, eV): min={float(np.min(dE0)):.6e} max={float(np.max(dE0)):.6e}")
        print(f"  dE_best    (raw, eV): min={float(np.min(dE1)):.6e} max={float(np.max(dE1)):.6e}")

        # Softclamp parameters and units
        unit_str = "kcal" if kcal_obj else "eV"
        j_unit_str = "kcal^2" if kcal_obj else "eV^2"
        print(f"MC DIAG: softclamp settings (applied to dE in {unit_str})")
        print(f"  soft_clamp={sc_enabled}")
        if sc_enabled:
            print(f"  clamp_start={sc_start:.3f} {unit_str}")
            print(f"  clamp_max={sc_max:.3f} {unit_str}")
        print(f"  Note: J = 0.5 * W * softclamp(dE)^2, so J units are {j_unit_str}")

        # Hard consistency check: plotted Etot arrays are those returned from the same kernel call as Ein.
        # Therefore internal consistency here is guaranteed by construction.
        print("MC DIAG: energy kernel consistency (Etot/Ei_masked)")
        print(f"  Before: Etot0 dtype={Etot0.dtype} shape={Etot0.shape} min={float(np.min(Etot0)):.6e} max={float(np.max(Etot0)):.6e}")
        print(f"  After : Etot1 dtype={Etot1.dtype} shape={Etot1.shape} min={float(np.min(Etot1)):.6e} max={float(np.max(Etot1)):.6e}")

        # RIGOROUS CHECK: optimizer dE must match (Emols - Eref) using energies stored by optimizer.
        # Note: Etot0/Etot1 here come from evaluate_energies_masked(mask_hbond) and are used for Ein/Eout split.
        # The objective, however, uses unmasked energies from driver.evaluate_energies() which are stored in history.
        Emols0_obj = history_after.get('Emols_initial', None)
        Emols1_obj = history_after.get('Emols_best', None)
        assert (Emols0_obj is not None) and (Emols1_obj is not None)
        dE_from_obj_before = Emols0_obj - Eref
        dE_from_obj_after  = Emols1_obj - Eref
        
        # Compare
        dE_diff_before = dE_from_obj_before - dE0
        dE_diff_after = dE_from_obj_after - dE1
        print("MC DIAG: RIGOROUS CHECK - objective dE (Emols-Eref) vs stored MC dE")
        print(f"  Before: max|dE_plot - dE_MC| = {float(np.max(np.abs(dE_diff_before))):.6e}   mean={float(np.mean(dE_diff_before)):.6e}")
        print(f"  After : max|dE_plot - dE_MC| = {float(np.max(np.abs(dE_diff_after))):.6e}   mean={float(np.mean(dE_diff_after)):.6e}")
        # Tolerance: 1e-3 eV (~0.023 kcal/mol) accounts for float32 precision differences
        if np.max(np.abs(dE_diff_before)) > 1e-3 or np.max(np.abs(dE_diff_after)) > 1e-3:
            print("  WARNING: PLOTTED Etot and MC dE ARE INCONSISTENT - different parameters or kernel!")
            print(f"  This is a BUG - the plotted energies don't match what the optimizer used!")
        else:
            print("  OK: plotted Etot and MC dE are consistent (within numerical tolerance)")
        
        # Create before/after comparison figure
        print("\nCreating before/after comparison plot...")
        fig_comp, axes_comp = plt.subplots(2, 6, figsize=(24, 8))

        # --- RIGOROUS PANEL CONSISTENCY CHECKS (grid-level) ---
        # The user-visible invariant must hold on the plotted grids:
        #   dE_grid_shown == Etot_grid_shown - Eref_grid_shown
        # Any discrepancy here means we plotted inconsistent sources/units.
        tol_grid = 1e-6  # eV-level; should be exactly identical up to float32 mapping noise

        # Per-sample check: Etot from masked kernel vs objective Emols stored in history
        # Also recompute unmasked energies *right now* to disambiguate bMask-kernel issues vs history/state issues.

        # Evaluate unmasked energies with EXACTLY the same DOF-derived tREQHs as Etot0/Etot1
        drv.apply_parameter_overrides(params_initial_path)
        drv.init_and_upload_energy_only()
        drv.toGPU_(drv.tREQHs_buff, T0)
        Emols0_now = drv.evaluate_energies()

        drv.toGPU_(drv.tREQHs_buff, T1)
        Emols1_now = drv.evaluate_energies()

        dEtot0_hist = Etot0 - Emols0_obj
        dEtot1_hist = Etot1 - Emols1_obj
        dEtot0_now  = Etot0 - Emols0_now
        dEtot1_now  = Etot1 - Emols1_now
        dEobj0      = Emols0_now - Emols0_obj
        dEobj1      = Emols1_now - Emols1_obj

        print("MC DIAG: RIGOROUS CHECK - Etot consistency across masked/unmasked/objective")
        print(f"  Before: max|Etot0(masked)-Emols0_obj(history)| = {float(np.max(np.abs(dEtot0_hist))):.6e}   mean={float(np.mean(dEtot0_hist)):.6e}")
        print(f"  Before: max|Etot0(masked)-Emols0_now(unmasked)| = {float(np.max(np.abs(dEtot0_now ))):.6e}   mean={float(np.mean(dEtot0_now )):.6e}")
        print(f"  Before: max|Emols0_now - Emols0_obj|           = {float(np.max(np.abs(dEobj0     ))):.6e}   mean={float(np.mean(dEobj0     )):.6e}")
        print(f"  After : max|Etot1(masked)-Emols1_obj(history)| = {float(np.max(np.abs(dEtot1_hist))):.6e}   mean={float(np.mean(dEtot1_hist)):.6e}")
        print(f"  After : max|Etot1(masked)-Emols1_now(unmasked)| = {float(np.max(np.abs(dEtot1_now ))):.6e}   mean={float(np.mean(dEtot1_now )):.6e}")
        print(f"  After : max|Emols1_now - Emols1_obj|           = {float(np.max(np.abs(dEobj1     ))):.6e}   mean={float(np.mean(dEobj1     )):.6e}")

        # Fail loudly if masked vs unmasked mismatch exists NOW under identical tREQHs
        if (np.max(np.abs(dEtot0_now)) > 1e-5) or (np.max(np.abs(dEtot1_now)) > 1e-5):
            raise RuntimeError("BUG: evaluate_energies_masked() Etot differs from evaluate_energies() under the same DOF-derived tREQHs")
        
        # Compute difference and error maps FROM MC OPTIMIZER DATA (authoritative)
        # Also map Eref used by optimizer for consistent Etot-Ref visualization
        V_Eref = _map_samples_to_grid(frame_idx, V_ref.shape, Eref)
        V_diff_before = _map_samples_to_grid(frame_idx, V_ref.shape, dE0)
        V_diff_after  = _map_samples_to_grid(frame_idx, V_ref.shape, dE1)

        # Grid check: Etot_grid - Eref_grid must equal dE_grid
        V_dE_from_panels_before = V_model_total - V_Eref
        V_dE_from_panels_after  = V_model_total_opt - V_Eref
        m0 = np.isfinite(V_dE_from_panels_before) & np.isfinite(V_diff_before)
        m1 = np.isfinite(V_dE_from_panels_after) & np.isfinite(V_diff_after)
        dd0 = (V_dE_from_panels_before[m0] - V_diff_before[m0]).astype(np.float64)
        dd1 = (V_dE_from_panels_after[m1]  - V_diff_after[m1]).astype(np.float64)
        print("MC DIAG: RIGOROUS CHECK - grid identity dE == Etot_panel - Eref_panel")
        print(f"  Before: max|Δ|={float(np.max(np.abs(dd0))):.6e} meanΔ={float(np.mean(dd0)):.6e}  frac(dE>0)={(float(np.mean(V_diff_before[m0]>0.0))):.3f}  frac(Etot-Eref>0)={(float(np.mean(V_dE_from_panels_before[m0]>0.0))):.3f}")
        print(f"  After : max|Δ|={float(np.max(np.abs(dd1))):.6e} meanΔ={float(np.mean(dd1)):.6e}  frac(dE>0)={(float(np.mean(V_diff_after[m1]>0.0))):.3f}  frac(Etot-Eref>0)={(float(np.mean(V_dE_from_panels_after[m1]>0.0))):.3f}")
        if (np.max(np.abs(dd0)) > tol_grid) or (np.max(np.abs(dd1)) > tol_grid):
            raise RuntimeError("BUG: dE panel does not equal Etot_panel - Eref_panel on the grid (units/source mismatch)")
        
        # Compute weighted error using the same formula as optimizer: J = 0.5 * sum(W * dE^2)
        # Use the error values already computed during MC optimization (from history)
        total_error_before = history_after['initial_error']
        total_error_after = history_after['best_error']
        
        # Use the per-sample error data stored during MC optimization (actual data used by optimizer)
        J_per_sample_before = Jps0
        J_per_sample_after  = Jps1

        # Map per-sample errors to grid for visualization (NO reshape assumptions)
        V_weighted_error_before = _map_samples_to_grid(frame_idx, V_ref.shape, J_per_sample_before)
        V_weighted_error_after  = _map_samples_to_grid(frame_idx, V_ref.shape, J_per_sample_after)

        print("MC DIAG: plotted grid ranges (what you actually see in col#4/#5)")
        print(f"  dE grid before: min={float(np.nanmin(V_diff_before)):.6e} max={float(np.nanmax(V_diff_before)):.6e}")
        print(f"  dE grid after : min={float(np.nanmin(V_diff_after )):.6e} max={float(np.nanmax(V_diff_after )):.6e}")
        print(f"  Jps grid before: min={float(np.nanmin(V_weighted_error_before)):.6e} max={float(np.nanmax(V_weighted_error_before)):.6e}")
        print(f"  Jps grid after : min={float(np.nanmin(V_weighted_error_after )):.6e} max={float(np.nanmax(V_weighted_error_after )):.6e}")
        
        # Row 1: Before optimization
        plot_system_panel(V_Eref, rv, A, axes_comp[0,0], 'Reference\n(Total)', args.kcal, sym=False, overlay_rmin=False, cmap=args.cmap)
        plot_system_panel(V_model_total, rv, A, axes_comp[0,1], 'Model Etot\n(Before)', args.kcal, sym=False, overlay_rmin=False, cmap=args.cmap)
        plot_system_panel(V_model_hbond, rv, A, axes_comp[0,2], 'Model Ein\n(Before)', args.kcal, sym=False, overlay_rmin=False, cmap=args.cmap)
        plot_system_panel(V_model_eout, rv, A, axes_comp[0,3], 'Model Eout\n(Before)', args.kcal, sym=False, overlay_rmin=False, cmap=args.cmap)
        plot_system_panel(V_diff_before, rv, A, axes_comp[0,4], f'Etot-Ref\n(Before)\nJ={total_error_before:.2f}', args.kcal, sym=False, overlay_rmin=False, cmap=args.cmap)
        # Weighted error: units are energy^2, label based on kcal_objective
        j_unit_label = "kcal²" if kcal_obj else "eV²"
        plot_system_panel(V_weighted_error_before, rv, A, axes_comp[0,5], f'Weighted Error\n(Before)\nJ={total_error_before:.2f}', False, sym=False, overlay_rmin=False, cmap='inferno', unit_label=j_unit_label)
        
        # Row 2: After optimization
        plot_system_panel(V_Eref, rv, A, axes_comp[1,0], 'Reference\n(Total)', args.kcal, sym=False, overlay_rmin=False, cmap=args.cmap)
        plot_system_panel(V_model_total_opt, rv, A, axes_comp[1,1], 'Model Etot\n(After)', args.kcal, sym=False, overlay_rmin=False, cmap=args.cmap)
        plot_system_panel(V_model_hbond_opt, rv, A, axes_comp[1,2], 'Model Ein\n(After)', args.kcal, sym=False, overlay_rmin=False, cmap=args.cmap)
        plot_system_panel(V_model_eout_opt, rv, A, axes_comp[1,3], 'Model Eout\n(After)', args.kcal, sym=False, overlay_rmin=False, cmap=args.cmap)
        plot_system_panel(V_diff_after, rv, A, axes_comp[1,4], f'Etot-Ref\n(After)\nJ={total_error_after:.2f}', args.kcal, sym=False, overlay_rmin=False, cmap=args.cmap)
        # Weighted error: units are energy^2, label based on kcal_objective
        plot_system_panel(V_weighted_error_after, rv, A, axes_comp[1,5], f'Weighted Error\n(After)\nJ={total_error_after:.2f}', False, sym=False, overlay_rmin=False, cmap='inferno', unit_label=j_unit_label)
        
        # Apply color limits:
        # - energy-like panels (0-3): use symmetric vmin/vmax from reference (same as main plot)
        # - dE panel (4): symmetric around 0 with its own range
        # - weighted error panel (5): non-negative with its own range (no kcal conversion)

        # Recompute vmin/vmax in this scope (we no longer keep these as outer-scope variables).
        vmin_ref = float(np.nanmin(V_Eref))
        if args.kcal:
            vmin_ref *= 23.060548
        vmin = args.Emin_factor * vmin_ref
        if args.Erange is None:
            vmax = -vmin
        else:
            vmax = vmin + args.Erange

        for irow in range(2):
            for icol in range(6):
                ax = axes_comp[irow, icol]
                if not ax.images: continue
                im = ax.images[0]
                if icol in (0, 1, 2, 3):
                    # Use same symmetric limits as main plot (reference-derived)
                    im.set_clim(vmin, vmax)
                elif icol == 4:
                    # dE map: symmetric around 0 with its own range
                    dEmax = float(np.nanmax(np.abs(V_diff_before)))
                    dEmax = max(dEmax, float(np.nanmax(np.abs(V_diff_after))))
                    if dEmax > 0:
                        _set_im_clim(im, -dEmax, dEmax, kcal=args.kcal)
                elif icol == 5:
                    # J_per_sample map: non-negative
                    Jmax = float(np.nanmax(V_weighted_error_before))
                    Jmax = max(Jmax, float(np.nanmax(V_weighted_error_after)))
                    if Jmax > 0:
                        # IMPORTANT: weighted error is energy^2; do NOT apply kcal conversion
                        _set_im_clim(im, 0.0, Jmax, kcal=False)

                # Diagnostics: print what is actually plotted vs clim (for error panels)
                if icol in (4, 5):
                    try:
                        arr = im.get_array()
                        amax = float(np.nanmax(arr))
                        amin = float(np.nanmin(arr))
                        c0, c1 = im.get_clim()
                        print(f"MC DIAG: im[{irow},{icol}] data[min,max]=({amin:.6e},{amax:.6e}) clim=({float(c0):.6e},{float(c1):.6e})")
                    except Exception:
                        pass
        
        # Compute correlation after optimization
        valid_mask_opt = np.isfinite(V_ref) & np.isfinite(V_model_total_opt)
        correlation_opt = 0.0
        if np.any(valid_mask_opt):
            correlation_opt = np.corrcoef(V_ref[valid_mask_opt], V_model_total_opt[valid_mask_opt])[0, 1]
        
        # Build parameter caption with error metrics and optimized parameters
        # Format optimized parameters for display
        opt_params_flat = []
        for typename, params in optimized_json.items():
            if isinstance(params, dict):
                for comp, val in params.items():
                    opt_params_flat.append(f"{typename}.{comp}={val:.3f}")
        opt_params_str = ", ".join(opt_params_flat[:6])  # Show first 6 parameters
        if len(opt_params_flat) > 6:
            opt_params_str += "..."
        
        soft_clamp_info = f", soft_clamp: {soft_clamp_use}" if soft_clamp_use else ""
        caption_comp = (f"MC Optimization: {args.max_steps} steps, step_size={args.step_size}, T={args.temperature}{soft_clamp_info}\n"
                       f"Total Error (J): Before={total_error_before:.4f}, After={total_error_after:.4f}, "
                       f"Correlation: Before={correlation:.4f}, After={correlation_opt:.4f}\n"
                       f"Optimized: {opt_params_str}")
        plt.suptitle(f'Before/After Comparison: {os.path.basename(args.xyz)}', fontsize=14)
        fig_comp.text(0.5, 0.02, caption_comp, ha='center', fontsize=9, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        plt.tight_layout(rect=[0, 0.08, 1, 0.95])  # Make room for caption at bottom
        comp_path = os.path.splitext(args.output)[0] + '_comparison.png'
        comp_path_abs = os.path.abspath(comp_path)
        plt.savefig(comp_path_abs, dpi=150, bbox_inches='tight')
        plt.close(fig_comp)
        print(f"Comparison plot: {comp_path_abs}")
        
        # Print optimized statistics
        print(f"\nOptimized Energy Statistics:")
        print(f"  Model Total (After): min={V_model_total_opt.min():.6f}, max={V_model_total_opt.max():.6f}, mean={V_model_total_opt.mean():.6f}")
        print(f"  Model H-bond (After): min={V_model_hbond_opt.min():.6f}, max={V_model_hbond_opt.max():.6f}, mean={V_model_hbond_opt.mean():.6f}")
        print(f"  Model Eout (After): min={V_model_eout_opt.min():.6f}, max={V_model_eout_opt.max():.6f}, mean={V_model_eout_opt.mean():.6f}")
        
        if np.any(valid_mask_opt):
            print(f"Reference-Model correlation (After): {correlation_opt:.6f}")
            print(f"Correlation improvement: {correlation_opt - correlation:.6f}")
            print(f"Total fitting error J: Before={total_error_before:.4f}, After={total_error_after:.4f}")

if __name__ == "__main__":
    main()
