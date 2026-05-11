#!/usr/bin/env python3
import os, sys, argparse
import numpy as np

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.realpath(os.path.join(_THIS_DIR, '..', '..', '..'))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from pyBall.DFTB import TestUtils as tu
from pyBall import plotUtils as pu
from pyBall.OCL import AFM as afm
from pyBall import atomicUtils as au
from scipy.optimize import curve_fit

def run_comparison(basis='mio-1-1', output_dir=None, xyz_path=None, target_idx=0, fb_dir=None, db_dir=None):
    FB_DIR = fb_dir if fb_dir else 'debug'
    DB_DIR = db_dir if db_dir else (f'debug_dftb_pentacene' if basis == 'mio-1-1' else f'debug_dftb_pentacene_{basis.replace("-", "_")}')
    OUT_DIR = output_dir if output_dir else f'debug_comparison_{basis.replace("-", "_")}'
    os.makedirs(OUT_DIR, exist_ok=True)

    # Read grid specs from logs (molecule-agnostic)
    orig_fb, ngrid_fb, step_fb = afm.read_grid_spec_from_log(os.path.join(FB_DIR, 'step1_density', 'log.txt'))
    if orig_fb is None:
        orig_fb = np.array([-11.16, -6.8, -4.2]); step_fb = 0.10
        print(f"WARNING: Could not read Fireball grid spec, using fallback")
    else:
        step_fb = step_fb if step_fb is not None else 0.10

    orig_db, ngrid_db, step_db = afm.read_grid_spec_from_log(os.path.join(DB_DIR, 'step1_density', 'log.txt'))
    if orig_db is None:
        orig_db = np.array([-11.36, -6.6, -4.2]); step_db = 0.15
        print(f"WARNING: Could not read DFTB grid spec, using fallback")
    else:
        step_db = step_db if step_db is not None else 0.15

    # Target atom position from XYZ
    if xyz_path and os.path.exists(xyz_path):
        pos, _, names, _, _ = au.load_xyz(xyz_path)
        pos = np.array(pos, dtype=np.float64)
        pos_atom = pos[target_idx]
        print(f"Target atom {target_idx} ({names[target_idx]}) at {pos_atom}")
    else:
        pos_atom = np.array([-6.0056, -0.8000, 0.0000])
        print(f"Using fallback target position: {pos_atom}")

    print("--- Loading Densities ---")
    rho_fb = np.load(os.path.join(FB_DIR, 'step1_density/rho_grid.npy'))
    rho_db = np.load(os.path.join(DB_DIR, 'step1_density/rho_grid.npy'))
    
    tu.check_density_integration(rho_fb, step_fb, expected_n=102.0, label="Fireball")
    tu.check_density_integration(rho_db, step_db, expected_n=102.0, label="DFTB")

    print("\n--- Extracting Z-Profiles ---")
    z_fb, r_fb = tu.get_z_profile(rho_fb, pos_atom, orig_fb, step_fb)
    z_db, r_db = tu.get_z_profile(rho_db, pos_atom, orig_db, step_db)
    
    path_dens = os.path.join(OUT_DIR, 'density_comparison.png')
    pu.plot_density_comparison(z_fb, r_fb, z_db, r_db, fname=path_dens)
    print(f"  Produced: {path_dens}")

    # Numerical summary: Density tail ratio
    z_contact = 2.1
    val_fb = np.interp(z_contact, z_fb, r_fb)
    val_db = np.interp(z_contact, z_db, r_db)
    print(f"\nNumerical Summary (at z={z_contact} A):")
    print(f"  Density FB: {val_fb:12.6f} e/A^3")
    print(f"  Density DB: {val_db:12.6f} e/A^3")
    print(f"  Ratio FB/DB: {val_fb/val_db:10.2f} x")

    print("\n--- Loading Pauli Fields ---")
    E_pauli_fb = np.load(os.path.join(FB_DIR, 'step3_pauli/E_Pauli_field.npy'))
    E_pauli_db = np.load(os.path.join(DB_DIR, 'step3_pauli/E_Pauli_field.npy'))
    
    _, ep_fb = tu.get_z_profile(E_pauli_fb, pos_atom, orig_fb, step_fb)
    _, ep_db = tu.get_z_profile(E_pauli_db, pos_atom, orig_db, step_db)
    
    path_pauli = os.path.join(OUT_DIR, 'pauli_comparison.png')
    pu.plot_pauli_comparison(z_fb, ep_fb, z_db, ep_db, fname=path_pauli)
    print(f"  Produced: {path_pauli}")
    
    val_p_fb = np.interp(z_contact, z_fb, ep_fb)
    val_p_db = np.interp(z_contact, z_db, ep_db)
    print(f"  Pauli FB: {val_p_fb:12.6f} eV")
    print(f"  Pauli DB: {val_p_db:12.6f} eV")
    print(f"  Ratio FB/DB: {val_p_fb/val_p_db:10.2f} x")

    print("\n--- Loading DFTB Reference Energy ---")
    z_ref, e_ref = tu.get_dftb_zscan_energies(os.path.join('debug_dftb_comparison', 'zscan_results.txt'))
    if z_ref is not None:
        # Relative energy (binding energy)
        e_ref = e_ref - e_ref[-1]
        
        # --- Formal Fitting ---
        print("\n--- Fitting Pauli Parameters (z=2.1 - 3.0 A) ---")
        # Mask for fitting range
        mask_ref = (z_ref >= 2.0) & (z_ref <= 3.0)
        z_fit = z_ref[mask_ref]
        e_fit = e_ref[mask_ref]
        
        # We assume E_pauli = A * overlap^beta
        # ln(E) = ln(A) + beta * ln(overlap)
        # However, it's easier to fit: ln(E) = -alpha_ref * z + C_ref
        # And ln(overlap) = -alpha_ovl * z + C_ovl
        
        # Log of reference energy
        log_e_ref = np.log(e_fit)
        # Interpolate overlap to z_fit
        # ep_db is A_old * overlap. Let's get raw overlap by dividing by 16.0
        ovl_fit = np.interp(z_fit, z_db, ep_db) / 16.0
        log_ovl = np.log(ovl_fit)
        
        # Slopes
        slope_ref, intercept_ref = np.polyfit(z_fit, log_e_ref, 1)
        slope_ovl, intercept_ovl = np.polyfit(z_fit, log_ovl, 1)
        
        beta_fitted = slope_ref / slope_ovl
        # ln(A_fitted) = intercept_ref - beta_fitted * intercept_ovl
        A_fitted = np.exp(intercept_ref - beta_fitted * intercept_ovl)
        
        print(f"  Fitted beta_pauli: {beta_fitted:10.4f}")
        print(f"  Fitted A_pauli:    {A_fitted:10.2f}")
        
        # Generate fitted curve
        e_fitted = A_fitted * (ovl_fit**beta_fitted)
        
        path_final = os.path.join(OUT_DIR, 'final_fit_comparison.png')
        pu.plot_compare_1d(z_fit, e_fit, z_fit, e_fitted, 
                           labels=["DFTB Ref", f"FDBM (A={A_fitted:.0f}, b={beta_fitted:.2f})"],
                           title="Final Pauli Fit (Log Scale)",
                           ylabel="Energy [eV]",
                           fname=path_final)
        print(f"  Produced: {path_final}")
        
        # Linear scale plot for verification
        path_final_lin = os.path.join(OUT_DIR, 'final_fit_comparison_lin.png')
        pu.plot_compare_1d(z_fit, e_fit, z_fit, e_fitted, 
                           labels=["DFTB Ref", "FDBM Fitted"],
                           title="Final Pauli Fit (Linear Scale)",
                           ylabel="Energy [eV]",
                           log=False,
                           fname=path_final_lin)
        print(f"  Produced: {path_final_lin}")

    print(f"\nAll plots and summary files are saved in: {os.path.abspath(OUT_DIR)}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--basis', type=str, default='mio-1-1', help='Basis set: mio-1-1 or 3ob-3-1')
    parser.add_argument('--xyz', type=str, default='pentacene.xyz', help='Molecule XYZ file')
    parser.add_argument('--target_idx', type=int, default=0, help='Target atom index')
    parser.add_argument('--fb_dir', type=str, default=None, help='Fireball output directory')
    parser.add_argument('--db_dir', type=str, default=None, help='DFTB output directory')
    parser.add_argument('--output_dir', type=str, default=None, help='Custom output directory')
    args = parser.parse_args()
    xyz_path = os.path.join(_THIS_DIR, args.xyz) if args.xyz else None
    run_comparison(basis=args.basis, output_dir=args.output_dir, xyz_path=xyz_path,
                   target_idx=args.target_idx, fb_dir=args.fb_dir, db_dir=args.db_dir)
