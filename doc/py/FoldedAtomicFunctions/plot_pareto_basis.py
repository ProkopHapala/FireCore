"""Utilities to visualise Pareto front of basis-set optimisation results.

Outputs generated in current directory:
  1. pareto_scatter.png           – All DB points (grey) + Pareto front (red).
  2. pf_<idx>_basis.png           – Basis-function shapes for Pareto point *idx*.
  3. pf_<idx>_fits.png            – Example fits for same basis.

Usage:
  python plot_pareto_basis.py --db results.jsonl --samples Morse1.json

The sample set JSON is the same format used by `Optimize_basis_batch.py`.
"""
from __future__ import annotations
import argparse, pathlib, sys, ast, json, random, datetime
from typing import List

import numpy as np
import matplotlib.pyplot as plt

# Local project modules
sys.path.append(str(pathlib.Path(__file__).parent))
from results_db import ResultsDB
from plot_db_utils import basis_key_to_struct
from basis_utils import (get_basis_metrics, construct_composite_cutoff_basis, weighted_rmse, load_morse_samples_json)
from Optimize_basis_batch import load_samples_json
from optimize import fit_coefficients
import plot_utils as pu

###############################################################################
# Helper functions
###############################################################################

bView = False
bSave = True

def metric_nb_rmse(key: str, rec: dict) -> tuple[int, float]:
    """Return (nBasis, RMSE) for given DB record."""
    _z0, bd = basis_key_to_struct(key)
    nb, *_ = get_basis_metrics(bd)
    return nb, rec["rmse"]

def scatter_all_and_pareto(db: ResultsDB, pf_keys: List[str], out_png: str):
    xs_all, ys_all, xs_pf, ys_pf = [], [], [], []
    for k, rec in db.items():
        x, y = metric_nb_rmse(k, rec)
        xs_all.append(x); ys_all.append(y)
        if k in pf_keys:
            xs_pf.append(x); ys_pf.append(y)
    fig, ax = plt.subplots(figsize=(6,4))
    ax.scatter(xs_all, ys_all, c='lightgray', s=20, marker='.')
    ax.scatter(xs_pf, ys_pf, c='tab:red', s=50, marker='o', label='Pareto')
    ax.set_xlabel('nBasis'); ax.set_ylabel('RMSE'); ax.set_yscale('log'); ax.grid(True, linestyle=':')
    ax.legend(); 
    fig.tight_layout(); 
    if bSave:     fig.savefig(out_png); 
    if not bView: plt.close(fig)
    print(f'Saved scatter to {out_png}')


def plot_basis_functions(idx: int, bd, rmse: float, z_grid: np.ndarray, z0: float):
    phi, _ = construct_composite_cutoff_basis(z_grid, bd, z0=z0)
    if phi.shape[0] == 0: return
    fig, ax = plt.subplots(figsize=(5,3))
    for row in phi:
        ax.plot(z_grid, row)
    ax.set_title(f'PF #{idx}  nBasis={phi.shape[0]}  RMSE={rmse:.2e}')
    ax.set_xlabel('z'); ax.set_ylabel('basis'); ax.set_ylim(-0.1,1.1)
    fig.tight_layout(); fname = f'pf_{idx}_nbas{phi.shape[0]}_basis.png'; 
    if bSave:     fig.savefig(fname); 
    if not bView: plt.close(fig)
    print('  basis ->', fname)


def plot_fits(idx: int, bd, z_grid: np.ndarray, z0: float,
              Y_samples: np.ndarray, weights: np.ndarray, prms_list=None):
    phi, _ = construct_composite_cutoff_basis(z_grid, bd, z0=z0)
    if phi.shape[0] == 0:
        return
    coeffs = fit_coefficients(Y_samples, phi, weights=weights)  # (P, M)
    Y_recon = coeffs.T @ phi  # (M, Nz)

    rng = np.random.default_rng(123 + idx)
    n_plot = min(5, Y_samples.shape[0])
    idx_sel = rng.choice(Y_samples.shape[0], size=n_plot, replace=False)
    pairs = []
    label_pairs = []
    for j in idx_sel:
        pairs.append((Y_samples[j], Y_recon[j]))
        if prms_list and j < len(prms_list):
            prm = prms_list[j]
            ref_lab = f"Morse[r0={prm['r0']:.2f},a={prm['a']:.2f}]"
        else:
            ref_lab = f"Sample {j}"
        coeff_vec = coeffs[:, j]
        coeff_label = '[' + ', '.join(f"{c:.3f}" for c in coeff_vec) + ']'
        app_lab = f"Coeff {coeff_label}"
        label_pairs.append((ref_lab, app_lab))

    title = f'PF #{idx} fits (nBasis={phi.shape[0]})'
    fig, _axes = pu.plotMultiFunctionApprox(z_grid, pairs, bError=True, errMax=0.05, scMin=1.2, title=title, label_pairs=label_pairs)
    fname = f'pf_{idx}_nbas{phi.shape[0]}_fits.png'
    if bSave:
        fig.savefig(fname)
    if not bView:
        plt.close(fig)
    print('  fits  ->', fname)


def plot_basis_and_fits(idx: int, bd, z_grid: np.ndarray, z0: float,
                        Y_samples: np.ndarray, weights: np.ndarray,
                        prms_list=None, rmse: float=None):
    """Plot basis functions (with labels) and corresponding fits in one figure."""
    phi, labels = construct_composite_cutoff_basis(z_grid, bd, z0=z0)
    if phi.shape[0] == 0:
        return
    coeffs = fit_coefficients(Y_samples, phi, weights=weights)
    Y_recon = coeffs.T @ phi

    fig, (ax_basis, ax_fit) = plt.subplots(2, 1, figsize=(10,14), sharex=False)
    # --- top: basis -------------------
    pu.plot1D(z_grid, phi, title=f'PF #{idx} basis (nBasis={phi.shape[0]})',
              labels=labels, ylims=(-0.1, 1.1), ax=ax_basis)
    ax_basis.set_xlabel('z')

    # --- bottom: sample fits -----------
    rng = np.random.default_rng(123 + idx)
    n_plot = min(5, Y_samples.shape[0])
    idx_sel = rng.choice(Y_samples.shape[0], size=n_plot, replace=False)
    pairs, label_pairs = [], []
    for j in idx_sel:
        pairs.append((Y_samples[j], Y_recon[j]))
        if prms_list and j < len(prms_list):
            prm = prms_list[j]
            ref_lab = f"Morse[r0={prm['r0']:.2f},a={prm['a']:.2f}]"
        else:
            ref_lab = f"Sample {j}"
        coeff_vec = coeffs[:, j]
        coeff_label = '[' + ', '.join(f"{c:.3f}" for c in coeff_vec) + ']'
        app_lab = f"Coeff {coeff_label}"
        label_pairs.append((ref_lab, app_lab))
    pu.plotMultiFunctionApprox(z_grid, pairs, bError=True, errMax=0.05, scMin=1.2,
                               title=f'PF #{idx} fits', label_pairs=label_pairs, ax1=ax_fit)
    fig.tight_layout()
    fname = f'pf_{idx}_nbas{phi.shape[0]}_basis_fits.png'
    if bSave:
        fig.savefig(fname)
    if not bView:
        plt.close(fig)
    print('  combined ->', fname)

###############################################################################
# Main
###############################################################################

if __name__ == '__main__':
    p = argparse.ArgumentParser(description='Plot Pareto front and basis details')
    p.add_argument('--db', default='results.jsonl')
    p.add_argument('--samples', default='Morse1.json')
    p.add_argument('--view', action='store_true')
    p.add_argument('--save', action='store_false')
    args = p.parse_args()

    bView = args.view
    bSave = args.save

    db = ResultsDB(args.db)
    if not db.keys():
        print('Database empty');
        sys.exit()

    # Load sample data (for fits)
    sample_fname = args.samples
    z0, z_grid, Y_samples, weights_mask, prms_list = load_morse_samples_json(sample_fname)
    # Morse parameters are already properly formatted from the shared function

    # Pareto with unique nBasis buckets handled by ResultsDB
    pf_keys = db.pareto_front(metric_nb_rmse, unique_x=True)
    print('Pareto (unique nBasis) size:', len(pf_keys))

    # Scatter plot
    scatter_all_and_pareto(db, pf_keys, 'pareto_scatter.png')

    # Individual basis & fits
    for idx, key in enumerate(pf_keys):
        rmse = db.get(key)['rmse']
        _z0, bd = basis_key_to_struct(key)
        print(f'PF #{idx}: key={key}')
        plot_basis_and_fits(idx, bd, z_grid, z0, Y_samples, weights_mask, prms_list, rmse)
    if bView: plt.show()
