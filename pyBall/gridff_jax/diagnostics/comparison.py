"""Multi-method comparison framework for MLIP benchmarking.

Generic comparison of any pair of evaluation methods (e.g., MAD-SURF vs DFT,
GridFF vs MAD-SURF, GridFF vs DFT) with parity plots, error maps, and
z-resolved / site-resolved statistics.

Usage
-----
    from pyBall.gridff_jax.diagnostics.comparison import (
        MethodResult, compare_methods, parity_statistics,
    )
    results = {
        "MAD-SURF": MethodResult(energies=..., forces=...),
        "GridFF":   MethodResult(energies=..., forces=...),
    }
    stats = compare_methods(results, reference="MAD-SURF")
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np


@dataclass
class MethodResult:
    """Container for evaluation results from any method."""
    energies: np.ndarray
    forces: np.ndarray | None = None
    pose_params: np.ndarray | None = None
    site_labels: list[str] | None = None
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass
class ComparisonStatistics:
    """Pairwise comparison statistics between a method and a reference."""
    method_name: str = ""
    reference_name: str = ""
    n_poses: int = 0
    energy_rmse_meV: float = 0.0
    energy_mae_meV: float = 0.0
    energy_max_error_meV: float = 0.0
    energy_bias_meV: float = 0.0
    force_rmse_eV_A: float = 0.0
    force_mae_eV_A: float = 0.0
    force_max_error_eV_A: float = 0.0
    energy_r2: float = 0.0
    per_site_energy_rmse: dict[str, float] = field(default_factory=dict)
    per_site_force_rmse: dict[str, float] = field(default_factory=dict)
    z_binned_energy_rmse: dict[str, float] = field(default_factory=dict)
    corrugation_ref_meV: float = 0.0
    corrugation_method_meV: float = 0.0
    metadata: dict[str, Any] = field(default_factory=dict)


def parity_statistics(
    energies_ref: np.ndarray,
    energies_pred: np.ndarray,
    forces_ref: np.ndarray | None = None,
    forces_pred: np.ndarray | None = None,
) -> dict[str, float]:
    """Compute basic parity statistics between reference and prediction."""
    e_ref = np.asarray(energies_ref, dtype=float).ravel()
    e_pred = np.asarray(energies_pred, dtype=float).ravel()
    e_err = e_pred - e_ref

    stats = {
        "energy_rmse_meV": float(np.sqrt(np.mean(e_err ** 2))) * 1000.0,
        "energy_mae_meV": float(np.mean(np.abs(e_err))) * 1000.0,
        "energy_max_error_meV": float(np.max(np.abs(e_err))) * 1000.0,
        "energy_bias_meV": float(np.mean(e_err)) * 1000.0,
        "n_poses": len(e_ref),
    }

    ss_res = np.sum(e_err ** 2)
    ss_tot = np.sum((e_ref - np.mean(e_ref)) ** 2)
    stats["energy_r2"] = float(1.0 - ss_res / max(ss_tot, 1.0e-30))

    if forces_ref is not None and forces_pred is not None:
        f_ref = np.asarray(forces_ref, dtype=float)
        f_pred = np.asarray(forces_pred, dtype=float)
        f_err = f_pred - f_ref
        flat_err = f_err.reshape(-1)
        stats["force_rmse_eV_A"] = float(np.sqrt(np.mean(flat_err ** 2)))
        stats["force_mae_eV_A"] = float(np.mean(np.abs(flat_err)))
        stats["force_max_error_eV_A"] = float(np.max(np.abs(flat_err)))

    return stats


def compare_methods(
    results: dict[str, MethodResult],
    reference: str,
    z_bins: list[tuple[float, float]] | None = None,
) -> dict[str, ComparisonStatistics]:
    """Compare multiple methods against a reference.

    Parameters
    ----------
    results : dict[str, MethodResult]
        Mapping of method name to evaluation results.
    reference : str
        Name of the reference method (must be a key in ``results``).
    z_bins : list[tuple[float, float]], optional
        Z-height bins for z-resolved statistics.  Each entry is (z_lo, z_hi).
        Default: [(1.5, 2.5), (2.5, 3.5), (3.5, 5.0), (5.0, 8.0)].

    Returns
    -------
    dict[str, ComparisonStatistics]
        One entry per non-reference method.
    """
    if reference not in results:
        raise KeyError(f"Reference method '{reference}' not found in results")

    ref = results[reference]
    e_ref = np.asarray(ref.energies, dtype=float).ravel()

    if z_bins is None:
        z_bins = [(1.5, 2.5), (2.5, 3.5), (3.5, 5.0), (5.0, 8.0)]

    comparisons = {}
    for name, method_result in results.items():
        if name == reference:
            continue

        e_pred = np.asarray(method_result.energies, dtype=float).ravel()
        e_err = e_pred - e_ref

        f_ref = np.asarray(ref.forces, dtype=float) if ref.forces is not None else None
        f_pred = np.asarray(method_result.forces, dtype=float) if method_result.forces is not None else None

        stats = ComparisonStatistics(
            method_name=name,
            reference_name=reference,
            n_poses=len(e_ref),
            energy_rmse_meV=float(np.sqrt(np.mean(e_err ** 2))) * 1000.0,
            energy_mae_meV=float(np.mean(np.abs(e_err))) * 1000.0,
            energy_max_error_meV=float(np.max(np.abs(e_err))) * 1000.0,
            energy_bias_meV=float(np.mean(e_err)) * 1000.0,
            corrugation_ref_meV=float(e_ref.max() - e_ref.min()) * 1000.0,
            corrugation_method_meV=float(e_pred.max() - e_pred.min()) * 1000.0,
        )

        ss_res = np.sum(e_err ** 2)
        ss_tot = np.sum((e_ref - np.mean(e_ref)) ** 2)
        stats.energy_r2 = float(1.0 - ss_res / max(ss_tot, 1.0e-30))

        if f_ref is not None and f_pred is not None:
            f_err = f_pred - f_ref
            flat_err = f_err.reshape(-1)
            stats.force_rmse_eV_A = float(np.sqrt(np.mean(flat_err ** 2)))
            stats.force_mae_eV_A = float(np.mean(np.abs(flat_err)))
            stats.force_max_error_eV_A = float(np.max(np.abs(flat_err)))

        # Per-site statistics
        labels = ref.site_labels or method_result.site_labels
        if labels is not None and len(labels) == len(e_ref):
            unique_labels = sorted(set(labels))
            for lbl in unique_labels:
                mask = np.array([lab == lbl for lab in labels])
                if mask.sum() == 0:
                    continue
                site_err = e_err[mask]
                stats.per_site_energy_rmse[lbl] = float(np.sqrt(np.mean(site_err ** 2))) * 1000.0
                if f_ref is not None and f_pred is not None:
                    site_f_err = (f_pred[mask] - f_ref[mask]).reshape(-1)
                    stats.per_site_force_rmse[lbl] = float(np.sqrt(np.mean(site_f_err ** 2)))

        # Z-binned statistics
        pose_params = ref.pose_params
        if pose_params is None:
            pose_params = method_result.pose_params
        if pose_params is not None and pose_params.shape[0] == len(e_ref):
            z_vals = pose_params[:, 2]
            for z_lo, z_hi in z_bins:
                mask = (z_vals >= z_lo) & (z_vals < z_hi)
                if mask.sum() == 0:
                    continue
                bin_err = e_err[mask]
                label = f"{z_lo:.1f}-{z_hi:.1f}"
                stats.z_binned_energy_rmse[label] = float(np.sqrt(np.mean(bin_err ** 2))) * 1000.0

        comparisons[name] = stats

    return comparisons
