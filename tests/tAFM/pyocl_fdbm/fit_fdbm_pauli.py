#!/usr/bin/env python3
"""
fit_fdbm_pauli.py - Fit FDBM Pauli parameters against DFTB SCF reference.

For each basis set (mio-1-1, 3ob-3-1) and each target atom:
  1. Load FDBM forcefield grids (Pauli, ES, vdW)
  2. Load DFTB rigid z-scan reference energies
  3. Extract profiles at matching z-positions
  4. Fit: E_DFTB(z) = A * overlap(z)^beta
  5. Save compact plots and parameters to dedicated output dir.

Usage:
    cd tests/tAFM/pyocl_fdbm
    PYTHONPATH=/path/to/FireCore:$PYTHONPATH python fit_fdbm_pauli.py --basis mio-1-1 --target_indices 0,1,20,21
    PYTHONPATH=/path/to/FireCore:$PYTHONPATH python fit_fdbm_pauli.py --basis all --target_indices 0,1,20,21 --xyz mymol.xyz
"""

import os, sys, argparse, json, time
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.realpath(os.path.join(_THIS_DIR, '..', '..', '..'))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from pyBall import atomicUtils as au
from pyBall.DFTB import TestUtils as tu
from pyBall.OCL import AFM as afm

A_PAULI_DEFAULT = 16.0

BASIS_CONFIG = {
    'mio-1-1': {'fdbm_dir': 'debug_dftb_mio_1_1'},
    '3ob-3-1': {'fdbm_dir': 'debug_dftb_3ob_3_1'},
}

# ── Helpers ──────────────────────────────────────────────────────────────────
def log(msg, log_path):
    print(msg)
    with open(log_path, 'a') as f:
        f.write(msg + '\n')

def load_fdbm_grids(fdbm_dir):
    paths = {
        'pauli': os.path.join(fdbm_dir, 'step3_pauli', 'E_Pauli_field.npy'),
        'es':    os.path.join(fdbm_dir, 'step4_electrostatics_conv', 'E_ES_field.npy'),
        'vdw':   os.path.join(fdbm_dir, 'step5_dispersion', 'E_vdw_field.npy'),
    }
    grids = {}
    for key, path in paths.items():
        grids[key] = np.load(path) if os.path.exists(path) else None
    return grids

def load_dftb_zscan(zscan_dir):
    z_path = os.path.join(zscan_dir, 'zscan_z.npy')
    e_path = os.path.join(zscan_dir, 'zscan_energy_eV.npy')
    if not (os.path.exists(z_path) and os.path.exists(e_path)):
        return None, None
    z = np.load(z_path)
    e = np.load(e_path)
    return z, e - e[-1]

# ── Fitting ──────────────────────────────────────────────────────────────────
def fit_pauli_powerlaw(z, overlap_raw, e_ref, z_min=2.0, z_max=3.5):
    mask = (z >= z_min) & (z <= z_max)
    if mask.sum() < 3:
        raise ValueError(f"Need >=3 points in fit range [{z_min},{z_max}]")
    z_fit = z[mask]
    o_fit = overlap_raw[mask]
    e_fit = e_ref[mask]
    pos_mask = (o_fit > 1e-15) & (e_fit > 1e-15)
    if pos_mask.sum() < 3:
        raise ValueError("Not enough positive points")
    log_o = np.log(o_fit[pos_mask])
    log_e = np.log(e_fit[pos_mask])
    beta_ll, lnA_ll = np.polyfit(log_o, log_e, 1)
    A_ll = np.exp(lnA_ll)
    def model(overlap, A, beta):
        return A * overlap**beta
    try:
        popt, _ = curve_fit(model, o_fit, e_fit, p0=[A_ll, beta_ll],
                            bounds=([0.0, 0.0], [1e6, 5.0]))
        A_nls, beta_nls = popt
        e_pred = model(o_fit, A_nls, beta_nls)
        ss_res = np.sum((e_fit - e_pred)**2)
        ss_tot = np.sum((e_fit - np.mean(e_fit))**2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
    except Exception as e:
        print(f"  WARNING: Nonlinear fit failed ({e}), using log-linear")
        A_nls, beta_nls = A_ll, beta_ll
        e_pred = model(o_fit, A_nls, beta_nls)
        r2 = 0.0
    return A_nls, beta_nls, r2, e_pred

# ── Plotting ───────────────────────────────────────────────────────────────
def plot_pauli_fit(z, e_ref, e_fitted, A, beta, fname, title, z_min=2.0, z_max=3.5):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    mask = (z >= z_min) & (z <= z_max)
    # Linear
    ax = axes[0]
    ax.plot(z, e_ref, 'o-', color='tab:blue', markersize=3, label='DFTB Ref', zorder=3)
    ax.plot(z[mask], e_fitted[mask], 's--', color='tab:red', markersize=3, label=f'Fit A={A:.0f} b={beta:.3f}', zorder=2)
    ax.axvspan(z_min, z_max, alpha=0.08, color='gray', label='Fit range')
    ax.set_xlabel('z [Å]')
    ax.set_ylabel('Energy [eV]')
    ax.set_title(f'{title} (Linear)')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    # Log
    ax = axes[1]
    pos = e_ref > 1e-12
    ax.semilogy(z[pos], e_ref[pos], 'o-', color='tab:blue', markersize=3, label='DFTB Ref')
    ax.semilogy(z[mask], e_fitted[mask], 's--', color='tab:red', markersize=3, label='Fit')
    ax.axvspan(z_min, z_max, alpha=0.08, color='gray')
    ax.set_xlabel('z [Å]')
    ax.set_ylabel('Energy [eV]')
    ax.set_title(f'{title} (Log)')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(fname, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {fname}")

def plot_summary(all_results, fname, basis, z_min, z_max):
    """Plot summary comparing all atoms."""
    n_atoms = len(all_results)
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    # Panel 1: A_pauli per atom
    ax = axes[0, 0]
    idxs = [r['idx'] for r in all_results]
    As = [r['A'] for r in all_results]
    ax.bar(idxs, As, color='tab:orange')
    ax.set_xlabel('Atom index')
    ax.set_ylabel('A_pauli')
    ax.set_title(f'A_pauli per atom ({basis})')
    ax.grid(True, alpha=0.3, axis='y')
    # Panel 2: beta per atom
    ax = axes[0, 1]
    betas = [r['beta'] for r in all_results]
    ax.bar(idxs, betas, color='tab:green')
    ax.set_xlabel('Atom index')
    ax.set_ylabel('beta_pauli')
    ax.set_title(f'beta_pauli per atom ({basis})')
    ax.grid(True, alpha=0.3, axis='y')
    # Panel 3: RMSE per atom
    ax = axes[1, 0]
    rmses = [r['rmse_fit'] for r in all_results]
    ax.bar(idxs, rmses, color='tab:red')
    ax.set_xlabel('Atom index')
    ax.set_ylabel('RMSE fit [eV]')
    ax.set_title(f'RMSE(fit range) per atom ({basis})')
    ax.grid(True, alpha=0.3, axis='y')
    # Panel 4: All fitted curves overlaid
    ax = axes[1, 1]
    for r in all_results:
        z = r['z']
        e_fit = r['e_fitted']
        ax.plot(z, e_fit, '-', lw=1.0, label=f"atom {r['idx']}")
    ax.set_xlabel('z [Å]')
    ax.set_ylabel('Fitted Pauli [eV]')
    ax.set_title(f'Fitted Pauli curves ({basis})')
    ax.set_yscale('log')
    ax.legend(fontsize=7, ncol=2)
    ax.grid(True, alpha=0.3)
    ax.axvspan(z_min, z_max, alpha=0.08, color='gray')
    plt.suptitle(f'Multi-Atom Summary: {basis}', fontsize=12)
    plt.tight_layout()
    plt.savefig(fname, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved summary: {fname}")

# ── Main workflow ──────────────────────────────────────────────────────────
def process_atom(basis, target_idx, target_name, target_pos, grids, zscan_dir, atom_out_dir, origin, step, ngrid, args):
    """Fit one target atom."""
    os.makedirs(atom_out_dir, exist_ok=True)
    log_path = os.path.join(atom_out_dir, 'fit.log')
    if os.path.exists(log_path):
        os.remove(log_path)

    log(f"Atom {target_idx} ({target_name}) at [{target_pos[0]:.4f}, {target_pos[1]:.4f}, {target_pos[2]:.4f}]", log_path)

    # Load DFTB z-scan
    z_ref, e_ref = load_dftb_zscan(zscan_dir)
    if z_ref is None:
        log(f"  ERROR: No z-scan found in {zscan_dir}", log_path)
        return None

    # Extract profiles using shared helper
    e_pauli = tu.extract_z_profile(grids['pauli'], target_pos, origin, step, z_distances=z_ref - target_pos[2])
    if e_pauli is None:
        log("  ERROR: Pauli field missing", log_path)
        return None

    # Overlap and fit
    overlap_raw = e_pauli / A_PAULI_DEFAULT
    overlap_safe = np.clip(overlap_raw, 1e-30, None)

    try:
        A_fit, beta_fit, r2, e_fitted_range = fit_pauli_powerlaw(
            z_ref, overlap_raw, e_ref, z_min=args.z_min, z_max=args.z_max
        )
    except Exception as e:
        log(f"  Fit failed: {e}", log_path)
        return None

    e_fitted_all = A_fit * overlap_safe**beta_fit
    mask_fit = (z_ref >= args.z_min) & (z_ref <= args.z_max)
    rmse_fit = np.sqrt(np.mean((e_ref[mask_fit] - e_fitted_range)**2))
    rmse_all = np.sqrt(np.mean((e_ref - e_fitted_all)**2))

    log(f"  A={A_fit:.2f}, beta={beta_fit:.4f}, R2={r2:.6f}, RMSE(fit)={rmse_fit:.4f}, RMSE(all)={rmse_all:.4f}", log_path)

    # Save
    params = {
        'basis': basis, 'atom_idx': target_idx, 'atom_name': target_name,
        'A_pauli': float(A_fit), 'beta_pauli': float(beta_fit),
        'R2_fit': float(r2), 'RMSE_fit': float(rmse_fit), 'RMSE_all': float(rmse_all),
        'fit_z_min': args.z_min, 'fit_z_max': args.z_max,
    }
    with open(os.path.join(atom_out_dir, 'params.json'), 'w') as f:
        json.dump(params, f, indent=2)

    np.save(os.path.join(atom_out_dir, 'z_ref.npy'), z_ref)
    np.save(os.path.join(atom_out_dir, 'e_ref.npy'), e_ref)
    np.save(os.path.join(atom_out_dir, 'e_pauli.npy'), e_pauli)
    np.save(os.path.join(atom_out_dir, 'overlap_raw.npy'), overlap_raw)

    # Plot
    plot_pauli_fit(
        z_ref, e_ref, e_fitted_all, A_fit, beta_fit,
        fname=os.path.join(atom_out_dir, 'fit_pauli.png'),
        title=f'{target_name}{target_idx} ({basis})',
        z_min=args.z_min, z_max=args.z_max
    )

    return {
        'idx': target_idx, 'name': target_name, 'pos': target_pos,
        'A': A_fit, 'beta': beta_fit, 'r2': r2,
        'rmse_fit': rmse_fit, 'rmse_all': rmse_all,
        'z': z_ref, 'e_fitted': e_fitted_all, 'e_ref': e_ref,
    }

def process_basis(basis, out_dir, args):
    t0 = time.time()
    os.makedirs(out_dir, exist_ok=True)

    xyz_path = os.path.join(_THIS_DIR, args.xyz)
    mol_pos, _, mol_names, _, _ = au.load_xyz(xyz_path)
    mol_pos = np.array(mol_pos, dtype=np.float64)
    target_indices = [int(x.strip()) for x in args.target_indices.split(',')]
    for idx in target_indices:
        if idx < 0 or idx >= len(mol_names):
            raise ValueError(f"Target index {idx} out of range (0-{len(mol_names)-1})")

    # Use explicit fdbm_dir if given, otherwise from BASIS_CONFIG
    fdbm_dir = os.path.join(_THIS_DIR, args.fdbm_dir) if args.fdbm_dir else os.path.join(_THIS_DIR, BASIS_CONFIG[basis]['fdbm_dir'])
    if not os.path.exists(fdbm_dir):
        raise FileNotFoundError(f"FDBM directory not found: {fdbm_dir}")

    grids = load_fdbm_grids(fdbm_dir)
    if grids['pauli'] is None:
        raise FileNotFoundError(f"Pauli field not found in {fdbm_dir}")

    # Read grid spec from FDBM step1 log (molecule-agnostic)
    log_path_grid = os.path.join(fdbm_dir, 'step1_density', 'log.txt')
    origin, ngrid, step = afm.read_grid_spec_from_log(log_path_grid)
    if origin is None:
        # Fallback: compute from molecule positions
        grid_spec, origin, ngrid = afm.setup_density_grid(mol_pos, step=0.15)
        step = 0.15
        print(f"  WARNING: Could not read grid from {log_path_grid}, computed from molecule")
    print(f"  Grid: origin={origin.round(2)} ngrid={ngrid} step={step}")

    zscan_dir_base = os.path.join(_THIS_DIR, args.zscan_dir) if args.zscan_dir else os.path.join(_THIS_DIR, f'zscan_{basis.replace("-", "_")}')

    log(f"{'='*70}", os.path.join(out_dir, 'fit.log'))
    log(f"FDBM Pauli Fitting: basis={basis}, atoms={target_indices}", os.path.join(out_dir, 'fit.log'))
    log(f"Molecule: {args.xyz} ({len(mol_names)} atoms)", os.path.join(out_dir, 'fit.log'))
    log(f"FDBM dir: {fdbm_dir}", os.path.join(out_dir, 'fit.log'))
    log(f"Output dir: {out_dir}", os.path.join(out_dir, 'fit.log'))

    all_results = []
    for target_idx in target_indices:
        target_name = mol_names[target_idx]
        target_pos = mol_pos[target_idx]
        atom_out_dir = os.path.join(out_dir, f'atom_{target_idx}')
        zscan_dir = os.path.join(zscan_dir_base, f'atom_{target_idx}')

        result = process_atom(basis, target_idx, target_name, target_pos, grids, zscan_dir, atom_out_dir, origin, step, ngrid, args)
        if result:
            all_results.append(result)

    if len(all_results) > 1:
        plot_summary(all_results, os.path.join(out_dir, 'summary_all_atoms.png'), basis, args.z_min, args.z_max)

    # Write summary table
    with open(os.path.join(out_dir, 'summary.txt'), 'w') as f:
        f.write("FDBM Pauli Fitting Summary\n")
        f.write("="*70 + "\n")
        f.write(f"Basis: {basis}\n")
        f.write(f"Atoms: {[r['idx'] for r in all_results]}\n")
        f.write(f"Fit range: z=[{args.z_min}, {args.z_max}] Å\n\n")
        f.write(f"{'Atom':>6} {'Name':>4} {'A_pauli':>10} {'beta':>8} {'R2':>10} {'RMSE_fit':>10} {'RMSE_all':>10}\n")
        f.write("-"*70 + "\n")
        for r in all_results:
            f.write(f"{r['idx']:6d} {r['name']:>4} {r['A']:10.2f} {r['beta']:8.4f} {r['r2']:10.6f} {r['rmse_fit']:10.4f} {r['rmse_all']:10.4f}\n")
        f.write(f"\nTime: {time.time()-t0:.1f}s\n")

    print(f"\nAll results saved to: {out_dir}/")
    return True

# ── CLI ────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description='Fit FDBM Pauli to DFTB reference')
    parser.add_argument('--basis', type=str, default='mio-1-1',
                        choices=['mio-1-1', '3ob-3-1', 'all'],
                        help='Basis set to fit')
    parser.add_argument('--xyz', type=str, default='pentacene.xyz',
                        help='Molecule XYZ file')
    parser.add_argument('--fdbm_dir', type=str, default=None,
                        help='FDBM pipeline output directory (default: debug_dftb_pentacene per basis)')
    parser.add_argument('--output_prefix', type=str, default='fit_pauli',
                        help='Output directory prefix')
    parser.add_argument('--zscan_dir', type=str, default=None,
                        help='Z-scan base directory (contains atom_X/ subdirs)')
    parser.add_argument('--target_indices', type=str, default='0',
                        help='Comma-separated target atom indices (e.g., "0,1,20,21")')
    parser.add_argument('--z_min', type=float, default=2.0,
                        help='Minimum z for Pauli fit range [Å]')
    parser.add_argument('--z_max', type=float, default=3.0,
                        help='Maximum z for Pauli fit range [Å]')
    args = parser.parse_args()

    basis_list = ['mio-1-1', '3ob-3-1'] if args.basis == 'all' else [args.basis]

    for basis in basis_list:
        out_dir = os.path.join(_THIS_DIR, f"{args.output_prefix}_{basis.replace('-', '_')}")
        try:
            process_basis(basis, out_dir, args)
        except Exception as e:
            print(f"\nERROR processing {basis}: {e}")
            import traceback
            traceback.print_exc()

if __name__ == "__main__":
    main()
