#!/usr/bin/env python3
import os
import sys
import numpy as np

# Add project root to path
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.realpath(os.path.join(_THIS_DIR, '..', '..', '..'))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from pyBall.DFTB import TestUtils as tu
from pyBall import plotUtils as pu
from pyBall.OCL import AFM as afm

def run_comparison():
    # --- Config ---
    FB_DIR = 'debug'
    DB_DIR = 'debug_dftb_pentacene'
    OUT_DIR = 'debug_comparison_v2'
    os.makedirs(OUT_DIR, exist_ok=True)

    # Fireball grid config
    orig_fb = np.array([-11.16, -6.8, -4.2])
    step_fb = 0.10
    # DFTB grid config
    orig_db = np.array([-11.36, -6.6, -4.2])
    step_db = 0.15
    
    # Target atom (Carbon in Pentacene)
    pos_atom = np.array([-6.0056, -0.8000, 0.0000])

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
        
        # Compare DFTB Ref vs FDBM Pauli (scaled)
        # Note: ep_db is already a z-slice. We need to align it with z_ref.
        ep_interp = np.interp(z_ref, z_db, ep_db)
        
        path_fit = os.path.join(OUT_DIR, 'fitted_comparison.png')
        pu.plot_compare_1d(z_ref, e_ref, z_ref, ep_interp * (1500/16.0), 
                           labels=["DFTB Total", "FDBM Pauli (A=1500)"],
                           title="Fitted Pauli vs DFTB Reference",
                           ylabel="Energy [eV]",
                           fname=path_fit)
        print(f"  Produced: {path_fit}")
        
        val_ref = np.interp(z_contact, z_ref, e_ref)
        val_fit = np.interp(z_contact, z_ref, ep_interp * (1500/16.0))
        print(f"  DFTB Total Ref: {val_ref:12.6f} eV")
        print(f"  FDBM Pauli (A=1500): {val_fit:12.6f} eV")
        print(f"  Error at contact: {val_fit - val_ref:10.4f} eV")

    print(f"\nAll plots and summary files are saved in: {os.path.abspath(OUT_DIR)}")

if __name__ == "__main__":
    run_comparison()
