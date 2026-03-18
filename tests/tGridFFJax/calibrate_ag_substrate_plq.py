#!/usr/bin/env python3
"""Calibrate Ag strict-PLQ substrate builder parameters on atomic probe datasets."""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.append(str(ROOT))

from pyBall.gridff_jax import load_config, save_config
from pyBall.gridff_jax.density_backends import make_density_backend
from pyBall.gridff_jax.fit import fit_hybrid_parameters, load_pose_dataset
from pyBall.gridff_jax.hybrid_energy import HybridGridFFModel
from pyBall.gridff_jax.utils import ensure_dir, save_json
from pyBall.gridff_jax.validation.metrics import compute_error_metrics


def parse_args():
    parser = argparse.ArgumentParser(description="Calibrate Ag substrate PLQ builder parameters")
    parser.add_argument("--config", default=str(ROOT / "tests/tGridFFJax/ag_dataset/config.json"))
    parser.add_argument("--dataset-dir", default=str(ROOT / "tests/tGridFFJax/ag_dataset/datasets"))
    parser.add_argument("--out-dir", default=str(ROOT / "tests/tGridFFJax/ag_substrate_calibration"))
    parser.add_argument("--adsorbates", default="H,C,O")
    parser.add_argument("--alpha-values", default="1.2,1.5,1.8")
    parser.add_argument("--radius-scales", default="0.9,1.0,1.1")
    parser.add_argument("--energy-scales", default="0.5,1.0,2.0")
    parser.add_argument("--max-steps", type=int, default=180)
    parser.add_argument("--force-weight", type=float, default=20.0)
    parser.add_argument("--learning-rate", type=float, default=1.0e-2)
    parser.add_argument("--prefer-jax", dest="prefer_jax", action="store_true")
    parser.add_argument("--no-prefer-jax", dest="prefer_jax", action="store_false")
    parser.set_defaults(prefer_jax=True)
    return parser.parse_args()


def _parse_values(text: str):
    return [float(item.strip()) for item in text.split(",") if item.strip()]


def _evaluate_all(model, datasets, params):
    energy_errors = []
    force_errors = []
    per_adsorbate = {}
    for poses, teacher in datasets:
        pred = model.evaluate_batch(poses, params=params)
        e_err = pred.energies - teacher.energies
        f_err = pred.forces - teacher.forces
        energy_errors.append(e_err)
        force_errors.append(f_err.reshape(-1))
        per_adsorbate[poses.adsorbate.name] = {
            "energy": compute_error_metrics(np.zeros_like(e_err), e_err),
            "force": compute_error_metrics(np.zeros(f_err.size, dtype=float), f_err.reshape(-1)),
            "n_samples": int(len(poses.positions)),
        }
    energy_concat = np.concatenate(energy_errors)
    force_concat = np.concatenate(force_errors)
    return {
        "aggregate": {
            "energy": compute_error_metrics(np.zeros_like(energy_concat), energy_concat),
            "force": compute_error_metrics(np.zeros_like(force_concat), force_concat),
        },
        "adsorbates": per_adsorbate,
    }


def main():
    args = parse_args()
    out_dir = ensure_dir(args.out_dir)
    config = load_config(args.config)
    config.grid.builder_mode = "metal_density_plq"
    config.grid.interpolation_order = 3
    config.training.max_steps = int(args.max_steps)
    config.training.force_weight = float(args.force_weight)
    config.training.learning_rate = float(args.learning_rate)
    config.training.fit_static_charge = False
    config.training.fit_sample_shift_z = False
    config.training.fit_coulomb_shift_z = False

    selected_adsorbates = {name.strip() for name in args.adsorbates.split(",") if name.strip()}
    datasets = []
    for path in sorted(Path(args.dataset_dir).glob("*.npz")):
        if selected_adsorbates and path.stem not in selected_adsorbates:
            continue
        poses, teacher, _ = load_pose_dataset(path)
        datasets.append((poses, teacher))
    if not datasets:
        raise ValueError("No datasets selected for substrate calibration")

    reactive_elements = tuple(sorted({symbol for poses, _ in datasets for symbol in poses.adsorbate.symbols}))
    alpha_values = _parse_values(args.alpha_values)
    radius_scales = _parse_values(args.radius_scales)
    energy_scales = _parse_values(args.energy_scales)

    trials = []
    best_record = None
    for alpha, radius_scale, energy_scale in itertools.product(alpha_values, radius_scales, energy_scales):
        trial_config = load_config(args.config)
        trial_config.grid.builder_mode = "metal_density_plq"
        trial_config.grid.interpolation_order = 3
        trial_config.toggles = config.toggles
        trial_config.hybrid_model = config.hybrid_model
        trial_config.training = config.training
        trial_config.grid.alpha_morse = float(alpha)
        trial_config.grid.req_scale_radius["Ag"] = float(radius_scale)
        trial_config.grid.req_scale_energy["Ag"] = float(energy_scale)
        density = make_density_backend(trial_config.density_backend, surface=trial_config.surface, grid=trial_config.grid).load()
        model = HybridGridFFModel(
            density=density,
            reactive_elements=reactive_elements,
            toggles=trial_config.toggles,
            grid_config=trial_config.grid,
            hybrid_config=trial_config.hybrid_model,
            prefer_jax=args.prefer_jax,
        )
        fit_result = fit_hybrid_parameters(
            density=density,
            datasets=datasets,
            model=model,
            training=trial_config.training,
        )
        metrics = _evaluate_all(model, datasets, fit_result.params)
        objective = float(metrics["aggregate"]["force"]["rmse"] + 0.25 * metrics["aggregate"]["energy"]["rmse"])
        record = {
            "alpha_morse": float(alpha),
            "ag_radius_scale": float(radius_scale),
            "ag_energy_scale": float(energy_scale),
            "objective": objective,
            "metrics": metrics,
            "fit_params": {
                "pauli": dict(fit_result.params.pauli),
                "london": dict(fit_result.params.london),
                "static_charge": dict(fit_result.params.static_charge),
            },
        }
        trials.append(record)
        if (best_record is None) or (objective < best_record["objective"]):
            best_record = record
            save_config(trial_config, out_dir / "best_config.json")
            save_json(best_record, out_dir / "best_result.json")

    trials.sort(key=lambda item: item["objective"])
    save_json(
        {
            "trials": trials,
            "best": best_record,
            "selected_adsorbates": sorted(selected_adsorbates),
        },
        out_dir / "calibration_results.json",
    )
    print(json.dumps(best_record, indent=2))


if __name__ == "__main__":
    main()
