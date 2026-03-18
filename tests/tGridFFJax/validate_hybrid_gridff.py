#!/usr/bin/env python3
"""Validate fitted hybrid GridFF artifacts against teacher datasets."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.append(str(ROOT))

from pyBall.gridff_jax import backend_summary, load_config
from pyBall.gridff_jax.config import FeatureToggles
from pyBall.gridff_jax.density_backends import make_density_backend
from pyBall.gridff_jax.fit import load_pose_dataset
from pyBall.gridff_jax.hybrid_energy import HybridGridFFModel, HybridParameters
from pyBall.gridff_jax.utils import ensure_dir, save_json
from pyBall.gridff_jax.validation import (
    compute_error_metrics,
    plot_convergence,
    plot_error_histogram,
    plot_grid_slices,
    plot_parity,
    plot_z_profile,
)


def parse_args():
    parser = argparse.ArgumentParser(description="Validate the fitted hybrid GridFF student model")
    parser.add_argument("--config", default=str(ROOT / "tests/tGridFFJax/ag_fit/config_used.json"), help="Path to config_used.json")
    parser.add_argument("--dataset-dir", default=str(ROOT / "tests/tGridFFJax/ag_dataset/datasets"), help="Directory containing *.npz dataset files")
    parser.add_argument("--fit-dir", default=str(ROOT / "tests/tGridFFJax/ag_fit"), help="Directory with fit outputs")
    parser.add_argument("--out-dir", default=str(ROOT / "tests/tGridFFJax/ag_fit/validation"), help="Validation output directory")
    parser.add_argument("--mode", default="", choices=["", "plq", "plq_reactive", "plq_ct", "plq_ct_image", "full"], help="Optional override of active physics terms during validation")
    parser.add_argument("--adsorbates", default="", help="Optional comma-separated adsorbate filter for focused validation")
    parser.add_argument("--alpha-morse", type=float, default=None, help="Optional override of substrate Morse alpha")
    parser.add_argument("--ag-radius-scale", type=float, default=None, help="Optional override of Ag substrate vdW radius scale")
    parser.add_argument("--ag-energy-scale", type=float, default=None, help="Optional override of Ag substrate vdW energy scale")
    parser.add_argument("--prefer-jax", dest="prefer_jax", action="store_true", help="Use the JAX student backend when available")
    parser.add_argument("--no-prefer-jax", dest="prefer_jax", action="store_false", help="Force the NumPy/SciPy fallback backend")
    parser.set_defaults(prefer_jax=True)
    return parser.parse_args()


def _load_params(path):
    with Path(path).open("r", encoding="utf8") as fin:
        payload = json.load(fin)
    return HybridParameters(
        pauli=payload["pauli"],
        london=payload["london"],
        reactive=payload["reactive"],
        static_charge=payload["static_charge"],
        req_radius_offset=payload.get("req_radius_offset", {}),
        req_energy_scale=payload.get("req_energy_scale", {}),
        chi=payload["chi"],
        hardness=payload["hardness"],
        image_scale=payload["image_scale"],
        image_plane=payload["image_plane"],
        image_damping=payload["image_damping"],
        sample_shift_z=payload.get("sample_shift_z", 0.0),
        coulomb_shift_z=payload.get("coulomb_shift_z", 0.0),
        reservoir_chi=payload["reservoir_chi"],
        reservoir_hardness=payload["reservoir_hardness"],
        total_charge=payload["total_charge"],
    )


def _mode_toggles(mode: str):
    if mode == "plq":
        return FeatureToggles(True, True, False, False, False, True)
    if mode == "plq_reactive":
        return FeatureToggles(True, True, False, False, True, True)
    if mode == "plq_ct":
        return FeatureToggles(True, True, True, False, False, True)
    if mode == "plq_ct_image":
        return FeatureToggles(True, True, True, True, False, True)
    return FeatureToggles(True, True, True, True, True, True)


def _z_bin_metrics(z_values, energy_errors, force_errors, edges):
    z_values = np.asarray(z_values, dtype=float)
    payload = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        if hi == edges[-1]:
            mask = (z_values >= lo) & (z_values <= hi)
        else:
            mask = (z_values >= lo) & (z_values < hi)
        if np.count_nonzero(mask) == 0:
            continue
        payload.append(
            {
                "z_min": float(lo),
                "z_max": float(hi),
                "n_samples": int(np.count_nonzero(mask)),
                "energy": compute_error_metrics(np.zeros(np.count_nonzero(mask), dtype=float), energy_errors[mask]),
                "force": compute_error_metrics(np.zeros(np.count_nonzero(mask) * force_errors.shape[1] * force_errors.shape[2], dtype=float), force_errors[mask].reshape(-1)),
            }
        )
    return payload


def main():
    args = parse_args()
    out_dir = ensure_dir(args.out_dir)
    config = load_config(args.config)
    config.grid.builder_mode = "metal_density_plq"
    config.grid.interpolation_order = 3
    if args.alpha_morse is not None:
        config.grid.alpha_morse = args.alpha_morse
    if args.ag_radius_scale is not None:
        config.grid.req_scale_radius["Ag"] = float(args.ag_radius_scale)
    if args.ag_energy_scale is not None:
        config.grid.req_scale_energy["Ag"] = float(args.ag_energy_scale)
    if args.mode:
        config.toggles = _mode_toggles(args.mode)
        config.hybrid_model.use_qeq = config.toggles.use_ct_qeq
        config.hybrid_model.use_image = config.toggles.use_image_charge
        config.hybrid_model.use_reactive_grid = config.toggles.use_reactive_grid
    params = _load_params(Path(args.fit_dir) / "fit_params.json")
    with (Path(args.fit_dir) / "fit_summary.json").open("r", encoding="utf8") as fin:
        fit_summary = json.load(fin)
    density = make_density_backend(config.density_backend, surface=config.surface, grid=config.grid).load()

    dataset_dir = Path(args.dataset_dir)
    selected_adsorbates = {name.strip() for name in args.adsorbates.split(",") if name.strip()}
    datasets = [load_pose_dataset(path) for path in sorted(dataset_dir.glob("*.npz")) if not selected_adsorbates or path.stem in selected_adsorbates]
    if not datasets:
        raise ValueError("No datasets selected for validation")
    reactive_elements = tuple(sorted({symbol for poses, _, _ in datasets for symbol in poses.adsorbate.symbols}))
    model = HybridGridFFModel(
        density=density,
        reactive_elements=reactive_elements,
        toggles=config.toggles,
        grid_config=config.grid,
        hybrid_config=config.hybrid_model,
        prefer_jax=args.prefer_jax,
    )

    split_indices = fit_summary["metrics"]["split_indices"]
    summary = {
        "backend": {
            "prefer_jax": bool(args.prefer_jax),
            "model_use_jax": bool(model.use_jax),
            "student_backend": backend_summary(),
        },
        "mode": args.mode or fit_summary.get("mode", "full"),
        "adsorbates_selected": sorted(selected_adsorbates) if selected_adsorbates else "all",
        "adsorbates": {},
    }
    aggregate = {name: {"energy": [], "force": []} for name in ("train", "val", "test", "all")}
    for poses, teacher, _ in datasets:
        eval_t0 = time.perf_counter()
        pred = model.evaluate_batch(poses, params=params)
        eval_elapsed = time.perf_counter() - eval_t0
        energy_metrics = compute_error_metrics(teacher.energies, pred.energies)
        force_metrics = compute_error_metrics(teacher.forces.reshape(-1), pred.forces.reshape(-1))
        adsorbate_summary = {
            "n_samples": int(len(poses.positions)),
            "seconds_total": eval_elapsed,
            "seconds_per_pose": eval_elapsed / max(len(poses.positions), 1),
            "energy": energy_metrics,
            "force": force_metrics,
            "splits": {},
        }
        energy_errors = pred.energies - teacher.energies
        force_errors_tensor = pred.forces - teacher.forces
        force_errors = force_errors_tensor.reshape(-1)
        adsorbate_summary["z_bins"] = _z_bin_metrics(
            poses.pose_params[:, 2],
            energy_errors,
            force_errors_tensor,
            edges=np.array([1.0, 2.0, 2.5, 3.0, 4.0, 6.0], dtype=float),
        )
        aggregate["all"]["energy"].append(energy_errors)
        aggregate["all"]["force"].append(force_errors)
        for split_name in ("train", "val", "test"):
            indices = np.asarray(split_indices[poses.adsorbate.name][split_name], dtype=int)
            if len(indices) == 0:
                continue
            split_energy_errors = energy_errors[indices]
            split_force_errors = (pred.forces[indices] - teacher.forces[indices]).reshape(-1)
            adsorbate_summary["splits"][split_name] = {
                "n_samples": int(len(indices)),
                "energy": compute_error_metrics(np.zeros_like(split_energy_errors), split_energy_errors),
                "force": compute_error_metrics(np.zeros_like(split_force_errors), split_force_errors),
            }
            aggregate[split_name]["energy"].append(split_energy_errors)
            aggregate[split_name]["force"].append(split_force_errors)
        summary["adsorbates"][poses.adsorbate.name] = adsorbate_summary

        plot_parity(teacher.energies, pred.energies, out_dir, f"{poses.adsorbate.name}_energy_parity.png")
        plot_parity(teacher.forces.reshape(-1), pred.forces.reshape(-1), out_dir, f"{poses.adsorbate.name}_force_parity.png")
        plot_error_histogram(pred.energies - teacher.energies, out_dir, f"{poses.adsorbate.name}_energy_error.png")
        plot_error_histogram(pred.energies - teacher.energies, out_dir, f"{poses.adsorbate.name}_energy_error_log.png", log_scale=True)
        plot_error_histogram(pred.forces - teacher.forces, out_dir, f"{poses.adsorbate.name}_force_error.png")
        plot_error_histogram(pred.forces - teacher.forces, out_dir, f"{poses.adsorbate.name}_force_error_log.png", log_scale=True)
        plot_z_profile(poses.pose_params[:, 2], teacher.energies, pred.energies, out_dir, f"{poses.adsorbate.name}_z_profile.png")
        for site in ("top", "bridge", "hollow"):
            mask = np.array([label == site for label in poses.site_labels], dtype=bool)
            if np.count_nonzero(mask) < 2:
                continue
            plot_z_profile(
                poses.pose_params[mask, 2],
                teacher.energies[mask],
                pred.energies[mask],
                out_dir,
                f"{poses.adsorbate.name}_{site}_z_profile.png",
            )

    summary["aggregate"] = {}
    for split_name, payload in aggregate.items():
        energy_concat = np.concatenate(payload["energy"])
        force_concat = np.concatenate(payload["force"])
        summary["aggregate"][split_name] = {
            "n_energy": int(energy_concat.size),
            "n_force": int(force_concat.size),
            "energy": compute_error_metrics(np.zeros_like(energy_concat), energy_concat),
            "force": compute_error_metrics(np.zeros_like(force_concat), force_concat),
        }
    plot_grid_slices(model.grids, out_dir, prefix="ag_grid")

    plot_convergence(fit_summary.get("history", []), out_dir, name_prefix="training_loss")

    save_json(summary, out_dir / "validation_metrics.json")
    print(f"validation complete -> {out_dir}")


if __name__ == "__main__":
    main()
