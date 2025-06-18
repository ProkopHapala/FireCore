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
from basis_utils import (get_basis_metrics, construct_composite_cutoff_basis,
                         weighted_rmse)
from Optimize_basis_batch import load_samples_json
from optimize import fit_coefficients

###############################################################################
# Helper functions
###############################################################################

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
    ax.legend(); fig.tight_layout(); fig.savefig(out_png); plt.close(fig)
    print(f'Saved scatter to {out_png}')


def plot_basis_functions(idx: int, bd, rmse: float, z_grid: np.ndarray, z0: float):
    phi, _ = construct_composite_cutoff_basis(z_grid, bd, z0=z0)
    if phi.shape[0] == 0: return
    fig, ax = plt.subplots(figsize=(5,3))
    for row in phi:
        ax.plot(z_grid, row)
    ax.set_title(f'PF #{idx}  nBasis={phi.shape[0]}  RMSE={rmse:.2e}')
    ax.set_xlabel('z'); ax.set_ylabel('basis'); ax.set_ylim(-0.1,1.1)
    fig.tight_layout(); fname = f'pf_{idx}_basis.png'; fig.savefig(fname); plt.close(fig)
    print('  basis ->', fname)


def plot_fits(idx: int, bd, z_grid: np.ndarray, z0: float,
              Y_samples: np.ndarray, weights: np.ndarray):
    phi, _ = construct_composite_cutoff_basis(z_grid, bd, z0=z0)
    if phi.shape[0] == 0: return
    coeffs = fit_coefficients(Y_samples, phi, weights=weights)
    Y_fit = coeffs.T @ phi
    rng = np.random.default_rng(123)
    idxs = rng.choice(Y_samples.shape[0], size=min(4, Y_samples.shape[0]), replace=False)
    fig, ax = plt.subplots(figsize=(5,3))
    for j in idxs:
        ax.plot(z_grid, Y_samples[j], 'k-', lw=1)
        ax.plot(z_grid, Y_fit[j], 'r--', lw=1)
    ax.set_title(f'PF #{idx} example fits'); ax.set_xlabel('z'); ax.set_ylabel('V(z)')
    fig.tight_layout(); fname = f'pf_{idx}_fits.png'; fig.savefig(fname); plt.close(fig)
    print('  fits  ->', fname)

###############################################################################
# Main
###############################################################################

def main():
    p = argparse.ArgumentParser(description='Plot Pareto front and basis details')
    p.add_argument('--db', default='results.jsonl')
    p.add_argument('--samples', default='Morse1.json')
    args = p.parse_args()

    db = ResultsDB(args.db)
    if not db.keys():
        print('Database empty'); return

    # Load sample data (for fits)
    z0, z_grid, Y_samples, weights = load_samples_json(args.samples)

    # Compute Pareto
    pf_keys = db.pareto_front(metric_nb_rmse)
    print('Pareto size:', len(pf_keys))

    # Scatter plot
    scatter_all_and_pareto(db, pf_keys, 'pareto_scatter.png')

    # Individual basis & fits
    for idx, key in enumerate(pf_keys):
        rmse = db.get(key)['rmse']
        _z0, bd = basis_key_to_struct(key)
        print(f'PF #{idx}: key={key}')
        plot_basis_functions(idx, bd, rmse, z_grid, z0)
        plot_fits(idx, bd, z_grid, z0, Y_samples, weights)

if __name__ == '__main__':
    main()
