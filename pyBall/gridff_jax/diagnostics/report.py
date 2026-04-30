"""Automated diagnostic report generation.

Aggregates results from smoothness, force consistency, corrugation, and
comparison analyses into a single structured report (JSON + optional plots).

Usage
-----
    from pyBall.gridff_jax.diagnostics.report import generate_diagnostic_report
    report = generate_diagnostic_report(
        smoothness_result=...,
        consistency_result=...,
        corrugation_result=...,
        comparison_stats=...,
        out_dir="/path/to/output",
    )
"""

from __future__ import annotations

import json
from dataclasses import asdict, is_dataclass
from pathlib import Path
from typing import Any

import numpy as np

from .comparison import ComparisonStatistics
from .corrugation import CorrugationResult
from .force_consistency import ConsistencyResult
from .smoothness import SmoothnessResult


def _to_serializable(obj: Any) -> Any:
    """Recursively convert dataclasses and numpy types to JSON-safe form."""
    if is_dataclass(obj) and not isinstance(obj, type):
        return {k: _to_serializable(v) for k, v in asdict(obj).items()
                if not isinstance(v, np.ndarray)}
    if isinstance(obj, dict):
        return {str(k): _to_serializable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_to_serializable(item) for item in obj]
    if isinstance(obj, np.ndarray):
        return None
    if isinstance(obj, (np.floating, np.integer)):
        return float(obj)
    return obj


def generate_diagnostic_report(
    smoothness_result: SmoothnessResult | None = None,
    consistency_result: ConsistencyResult | None = None,
    corrugation_result: CorrugationResult | None = None,
    comparison_stats: dict[str, ComparisonStatistics] | None = None,
    out_dir: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Generate a comprehensive diagnostic report.

    Parameters
    ----------
    smoothness_result : SmoothnessResult, optional
        From ``smoothness_analysis()``.
    consistency_result : ConsistencyResult, optional
        From ``force_energy_consistency()``.
    corrugation_result : CorrugationResult, optional
        From ``corrugation_analysis()``.
    comparison_stats : dict[str, ComparisonStatistics], optional
        From ``compare_methods()``.
    out_dir : str or Path, optional
        If provided, save JSON report and any plots here.
    verbose : bool
        Print summary to stdout.

    Returns
    -------
    dict
        Structured report dictionary.
    """
    report: dict[str, Any] = {"summary": {}, "sections": {}}

    if consistency_result is not None:
        section = {
            "n_poses": consistency_result.n_poses,
            "fd_step_angstrom": consistency_result.fd_step,
            "max_inconsistency_eV_A": consistency_result.max_inconsistency,
            "rms_inconsistency_eV_A": consistency_result.rms_inconsistency,
            "mae_inconsistency_eV_A": consistency_result.mae_inconsistency,
            "worst_pose_index": consistency_result.worst_pose_index,
            "worst_atom_index": consistency_result.worst_atom_index,
            "worst_component": ["x", "y", "z"][consistency_result.worst_component],
        }
        if consistency_result.per_component_rms is not None:
            section["per_component_rms"] = {
                c: float(consistency_result.per_component_rms[i])
                for i, c in enumerate(["x", "y", "z"])
            }
        report["sections"]["force_consistency"] = section
        report["summary"]["force_consistency_rms_eV_A"] = consistency_result.rms_inconsistency

    if smoothness_result is not None:
        section = {
            "n_grid": smoothness_result.n_grid,
            "roughness_per_z": {
                f"{z:.2f}": v for z, v in smoothness_result.roughness_per_z.items()
            },
            "max_gradient_per_z": {
                f"{z:.2f}": v for z, v in smoothness_result.max_gradient_per_z.items()
            },
        }
        report["sections"]["smoothness"] = section
        if smoothness_result.roughness_per_z:
            report["summary"]["avg_roughness"] = float(
                np.mean(list(smoothness_result.roughness_per_z.values()))
            )

    if corrugation_result is not None:
        section: dict[str, Any] = {"n_grid": corrugation_result.n_grid}
        for method_name in corrugation_result.corrugation_amplitude:
            method_sec = {}
            for z_h in corrugation_result.z_heights:
                amp = corrugation_result.corrugation_amplitude.get(method_name, {}).get(z_h, 0.0)
                rms = corrugation_result.corrugation_rms.get(method_name, {}).get(z_h, 0.0)
                method_sec[f"z={z_h:.2f}"] = {
                    "amplitude_meV": amp,
                    "rms_meV": rms,
                }
            section[method_name] = method_sec
        report["sections"]["corrugation"] = section

    if comparison_stats is not None:
        section = {}
        for method_name, stats in comparison_stats.items():
            section[method_name] = {
                "reference": stats.reference_name,
                "energy_rmse_meV": stats.energy_rmse_meV,
                "energy_mae_meV": stats.energy_mae_meV,
                "energy_max_error_meV": stats.energy_max_error_meV,
                "energy_bias_meV": stats.energy_bias_meV,
                "energy_r2": stats.energy_r2,
                "force_rmse_eV_A": stats.force_rmse_eV_A,
                "force_mae_eV_A": stats.force_mae_eV_A,
                "corrugation_ref_meV": stats.corrugation_ref_meV,
                "corrugation_method_meV": stats.corrugation_method_meV,
                "per_site_energy_rmse": stats.per_site_energy_rmse,
                "z_binned_energy_rmse": stats.z_binned_energy_rmse,
            }
        report["sections"]["comparison"] = section

    if verbose:
        print("\n" + "=" * 60)
        print("  DIAGNOSTIC REPORT SUMMARY")
        print("=" * 60)
        for key, val in report["summary"].items():
            print(f"  {key}: {val:.6f}" if isinstance(val, float) else f"  {key}: {val}")
        if comparison_stats:
            for method_name, stats in comparison_stats.items():
                print(f"\n  {method_name} vs {stats.reference_name}:")
                print(f"    Energy RMSE = {stats.energy_rmse_meV:.1f} meV")
                print(f"    Force  RMSE = {stats.force_rmse_eV_A:.4f} eV/A")
                print(f"    R^2         = {stats.energy_r2:.4f}")
        print("=" * 60)

    if out_dir is not None:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        report_path = out_dir / "diagnostic_report.json"
        serializable = _to_serializable(report)
        with report_path.open("w", encoding="utf-8") as fh:
            json.dump(serializable, fh, indent=2, default=str)
        if verbose:
            print(f"\n[report] Saved to {report_path}")

    return report
