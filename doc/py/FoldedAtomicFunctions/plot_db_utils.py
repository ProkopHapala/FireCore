"""Plotting helpers for results database."""
from __future__ import annotations
import matplotlib.pyplot as plt
from typing import List, Tuple, Dict
import json, ast
from basis_utils import get_basis_metrics


def basis_key_to_struct(key: str):
    """Parse key of form '<z0>,<basis_str>' back to (z0, basis_def)."""
    z0_str, basis_str = key.split(',', 1)
    z0 = float(z0_str)
    # ast.literal_eval is safe for lists/tuples of literals
    basis_def = ast.literal_eval(basis_str)
    return z0, basis_def


def scatter_rmse_vs_nbasis(db_items: Dict[str, dict], *, show=True, ax=None):
    """Plot RMSE against nBasis for records in *db_items*."""
    n_basis_vals: List[int] = []
    rmses: List[float] = []

    for key, rec in db_items.items():
        rmse = rec.get("rmse")
        if rmse is None:
            continue
        _z0, bd = basis_key_to_struct(key)
        nb, *_ = get_basis_metrics(bd)
        n_basis_vals.append(nb)
        rmses.append(rmse)

    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 4))
    sc = ax.scatter(n_basis_vals, rmses, c='tab:blue')
    ax.set_xlabel('nBasis')
    ax.set_ylabel('RMSE')
    ax.set_yscale('log')
    ax.grid(True, linestyle=':')
    if show:
        plt.tight_layout()
        plt.show()
    return ax
