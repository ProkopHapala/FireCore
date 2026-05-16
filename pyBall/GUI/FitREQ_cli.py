#!/usr/bin/env python3
import sys
import os
import numpy as np
import matplotlib.pyplot as plt

# Add FireCore to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from pyBall.OCL.FittingDriver import FittingDriver
from pyBall.OCL.NonBondFitting import optimizer_montecarlo, plot_mc_convergence, plot_mc_energy_maps, plot_mc_error_maps
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

def run_mc_cli(xyz_file, atom_types_file, dof_file, model_name='ENERGY_MorseQ_PAIR',
               max_steps=500, step_size=0.1, temperature=0.0,
               soft_clamp=False, clamp_start=4.0, clamp_max=6.0,
               out_dir="plots_mc"):
    """Run Monte Carlo optimizer and save convergence + energy-map plots."""
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    print(f"--- FitREQ Monte Carlo Optimizer ---")
    print(f"  XYZ:       {xyz_file}")
    print(f"  AtomTypes: {atom_types_file}")
    print(f"  DOFs:      {dof_file}")
    print(f"  Model:     {model_name}")
    print(f"  max_steps={max_steps}  step_size={step_size}  temperature={temperature}")
    if soft_clamp:
        print(f"  Soft clamp: start={clamp_start}, max={clamp_max}")

    # 1. Setup driver (energy-only path)
    drv = FittingDriver(verbose=0)
    drv.load_atom_types(atom_types_file)
    drv.compile_all_with_model(model_name, bPrint=False)
    drv.load_data(xyz_file)
    drv.load_dofs(dof_file)
    drv.init_and_upload_energy_only()
    drv.setup_energy_kernel()

    # 2. Load geometry for 2D energy maps
    r, a, rows, Ps = drv.load_scan_geometry(xyz_file)
    Erefs = drv.host_ErefW[:, 0] * EV_TO_KCAL
    Vref, _, Arow, rv, rows = reshape_to_grid_proper(Erefs, r, a, rows)

    # 3. Evaluate initial model energies
    initial_dofs = np.array([d['xstart'] for d in drv.dof_definitions], dtype=np.float64)
    Em_initial = drv.evaluate_energies() * EV_TO_KCAL
    Vmod_initial, _, _, _, _ = reshape_to_grid_proper(Em_initial, r, a, rows)

    # 4. Run Monte Carlo
    history = optimizer_montecarlo(drv, initial_dofs,
                                   max_steps=max_steps, step_size=step_size,
                                   temperature=temperature,
                                   soft_clamp=soft_clamp, clamp_start=clamp_start, clamp_max=clamp_max,
                                   verbose=max(1, max_steps // 20))

    # 5. Evaluate final model energies with best DOFs
    Em_final = drv.evaluate_energies() * EV_TO_KCAL   # tREQHs already updated by last evaluate_objective call
    # Re-evaluate explicitly with best DOFs to be sure
    drv.evaluate_objective(history['best_dofs'])
    Em_final = drv.evaluate_energies() * EV_TO_KCAL
    Vmod_final, _, _, _, _ = reshape_to_grid_proper(Em_final, r, a, rows)

    # 6. Compute per-sample errors for error maps
    J_initial_per_sample = drv.evaluate_objective_per_sample(initial_dofs, soft_clamp=False)
    J_final_per_sample = drv.evaluate_objective_per_sample(history['best_dofs'], soft_clamp=False)
    if soft_clamp:
        J_softclamp_per_sample = drv.evaluate_objective_per_sample(history['best_dofs'],
                                                                    soft_clamp=True, clamp_start=clamp_start, clamp_max=clamp_max)
    else:
        J_softclamp_per_sample = J_final_per_sample

    # Reshape to 2D grid
    J_initial_grid, _, _, _, _ = reshape_to_grid_proper(J_initial_per_sample, r, a, rows)
    J_final_grid, _, _, _, _ = reshape_to_grid_proper(J_final_per_sample, r, a, rows)
    J_softclamp_grid, _, _, _, _ = reshape_to_grid_proper(J_softclamp_per_sample, r, a, rows)

    # 7. Save plots
    conv_path = os.path.join(out_dir, "mc_convergence.png")
    maps_path = os.path.join(out_dir, "mc_energy_maps.png")
    err_path = os.path.join(out_dir, "mc_error_maps.png")

    fig_conv = plot_mc_convergence(history, out_path=conv_path)
    plt.close(fig_conv)

    fig_maps = plot_mc_energy_maps(Vref, Vmod_initial, Vmod_final, rv, Arow, out_path=maps_path)
    plt.close(fig_maps)

    fig_err = plot_mc_error_maps(J_initial_grid, J_final_grid, J_softclamp_grid,
                                  rv, Arow, r, a,
                                  soft_clamp=soft_clamp, clamp_start=clamp_start, clamp_max=clamp_max,
                                  out_path=err_path)
    plt.close(fig_err)

    print(f"\n--- MC Results ---")
    print(f"  Initial error : {history['initial_error']:.6e}")
    print(f"  Best error    : {history['best_error']:.6e}")
    print(f"  Steps accepted: {history['n_accepted']} / {history['n_steps']}")
    comp_labels = {0:'R', 1:'E', 2:'Q', 3:'H'}
    print(f"  Best DOFs:")
    for i, d in enumerate(drv.dof_definitions):
        print(f"    {d['typename']:8s} {comp_labels.get(d['comp'],'?')} : {initial_dofs[i]:.6f} -> {history['best_dofs'][i]:.6f}  (range [{d['min']:.4f}, {d['max']:.4f}])")
    print(f"\n  Plots saved to: {os.path.abspath(out_dir)}/")
    print(f"    {conv_path}")
    print(f"    {maps_path}")
    print(f"    {err_path}")
    return history


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="FitREQ Headless Diagnostic / MC Optimizer Tool")
    parser.add_argument("--mode", type=str, default="diag", choices=["diag", "mc"],
                        help="Run mode: 'diag' = energy diagnostics, 'mc' = Monte Carlo optimizer")
    parser.add_argument("--xyz", type=str, default="../../tests/tFitREQ_PN/wb97m-split/H2O-A1_H2O-D1-y.xyz", help="Input XYZ scan file")
    parser.add_argument("--atypes", type=str, default="../../tests/tFitREQ_PN/data/AtomTypes.dat", help="Atom types parameter file")
    parser.add_argument("--dofs", type=str, default="../../tests/tFitREQ_PN/data/dofSelection_H2O.dat", help="DOF selection file (for MC mode)")
    parser.add_argument("--model", type=str, default="ENERGY_MorseQ_PAIR", help="Energy model macro name")
    parser.add_argument("--g1", type=str, default="0,1,2", help="Atom indices in group 1 (comma separated)")
    parser.add_argument("--g2", type=str, default="0,1,2", help="Atom indices in group 2 (comma separated)")
    parser.add_argument("--out", type=str, default="plots_headless", help="Output directory")
    parser.add_argument("--max_steps", type=int, default=500, help="MC: max number of trial steps")
    parser.add_argument("--step_size", type=float, default=0.1, help="MC: step size as fraction of parameter range")
    parser.add_argument("--temperature", type=float, default=0.0, help="MC: temperature for simulated annealing (0=greedy)")
    parser.add_argument("--soft_clamp", action="store_true", help="Enable soft clamp on energy differences")
    parser.add_argument("--clamp_start", type=float, default=4.0, help="Soft clamp: start clamping when |dE| > start")
    parser.add_argument("--clamp_max", type=float, default=6.0, help="Soft clamp: saturate toward max value")
    args = parser.parse_args()

    base_dir = os.path.dirname(os.path.abspath(__file__))
    def _abs(p): return os.path.abspath(os.path.join(base_dir, p)) if not os.path.isabs(p) else p
    xyz    = _abs(args.xyz)
    atypes = _abs(args.atypes)
    dofs   = _abs(args.dofs)

    if args.mode == "mc":
        run_mc_cli(xyz, atypes, dofs, args.model,
                   max_steps=args.max_steps, step_size=args.step_size,
                   temperature=args.temperature,
                   soft_clamp=args.soft_clamp, clamp_start=args.clamp_start, clamp_max=args.clamp_max,
                   out_dir=args.out)
    else:
        run_cli(xyz, atypes, args.model, args.g1, args.g2, args.out)
