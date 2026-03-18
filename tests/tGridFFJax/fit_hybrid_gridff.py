#!/usr/bin/env python3
"""Fit the hybrid GridFF student model and export FireCore artifacts."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.append(str(ROOT))

from pyBall.gridff_jax import backend_summary, load_config, save_config
from pyBall.gridff_jax.config import FeatureToggles
from pyBall.gridff_jax.density_backends import make_density_backend
from pyBall.gridff_jax.export import export_firecore_artifacts
from pyBall.gridff_jax.fit import fit_hybrid_parameters, load_pose_dataset
from pyBall.gridff_jax.hybrid_energy import HybridGridFFModel, HybridParameters
from pyBall.gridff_jax.utils import ensure_dir, save_json


def parse_args():
    parser = argparse.ArgumentParser(description="Fit the hybrid GridFF student model")
    parser.add_argument("--config", default=str(ROOT / "tests/tGridFFJax/ag_dataset/config.json"), help="Path to dataset config.json")
    parser.add_argument("--dataset-dir", default=str(ROOT / "tests/tGridFFJax/ag_dataset/datasets"), help="Directory containing *.npz dataset files")
    parser.add_argument("--out-dir", default=str(ROOT / "tests/tGridFFJax/ag_fit"), help="Output directory for fitted artifacts")
    parser.add_argument("--max-steps", type=int, default=None, help="Override optimizer max iterations")
    parser.add_argument("--force-weight", type=float, default=None, help="Override force loss weight")
    parser.add_argument("--learning-rate", type=float, default=None, help="Override optimizer learning rate")
    parser.add_argument("--mode", default="plq", choices=["plq", "plq_reactive", "plq_ct", "plq_ct_image", "full"], help="Select which physics terms are active during fitting")
    parser.add_argument("--adsorbates", default="", help="Optional comma-separated adsorbate filter for focused fits")
    parser.add_argument("--init-fit-dir", default="", help="Optional directory containing fit_params.json used to initialize the optimizer")
    parser.add_argument("--alpha-morse", type=float, default=None, help="Optional override of substrate Morse alpha")
    parser.add_argument("--ag-radius-scale", type=float, default=None, help="Optional override of Ag substrate vdW radius scale")
    parser.add_argument("--ag-energy-scale", type=float, default=None, help="Optional override of Ag substrate vdW energy scale")
    parser.add_argument("--prefer-jax", dest="prefer_jax", action="store_true", help="Use the JAX student backend when available")
    parser.add_argument("--no-prefer-jax", dest="prefer_jax", action="store_false", help="Force the NumPy/SciPy fallback backend")
    parser.set_defaults(prefer_jax=True)
    return parser.parse_args()


def _params_to_dict(params):
    return {
        "pauli": dict(params.pauli),
        "london": dict(params.london),
        "reactive": dict(params.reactive),
        "static_charge": dict(params.static_charge),
        "req_radius_offset": dict(params.req_radius_offset),
        "req_energy_scale": dict(params.req_energy_scale),
        "chi": dict(params.chi),
        "hardness": dict(params.hardness),
        "image_scale": params.image_scale,
        "image_plane": params.image_plane,
        "image_damping": params.image_damping,
        "sample_shift_z": params.sample_shift_z,
        "coulomb_shift_z": params.coulomb_shift_z,
        "reservoir_chi": params.reservoir_chi,
        "reservoir_hardness": params.reservoir_hardness,
        "total_charge": params.total_charge,
    }


def _load_params(path):
    payload = Path(path)
    with payload.open("r", encoding="utf8") as fin:
        data = json.load(fin)
    return HybridParameters(
        pauli=data["pauli"],
        london=data["london"],
        reactive=data["reactive"],
        static_charge=data["static_charge"],
        req_radius_offset=data.get("req_radius_offset", {}),
        req_energy_scale=data.get("req_energy_scale", {}),
        chi=data["chi"],
        hardness=data["hardness"],
        image_scale=data["image_scale"],
        image_plane=data["image_plane"],
        image_damping=data["image_damping"],
        sample_shift_z=data.get("sample_shift_z", 0.0),
        coulomb_shift_z=data.get("coulomb_shift_z", 0.0),
        reservoir_chi=data["reservoir_chi"],
        reservoir_hardness=data["reservoir_hardness"],
        total_charge=data["total_charge"],
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


def _apply_mode(config, mode: str):
    config.grid.builder_mode = "metal_density_plq"
    config.grid.interpolation_order = 3
    config.toggles = _mode_toggles(mode)
    config.hybrid_model.use_qeq = config.toggles.use_ct_qeq
    config.hybrid_model.use_image = config.toggles.use_image_charge
    config.hybrid_model.use_reactive_grid = config.toggles.use_reactive_grid
    config.training.fit_chi = config.toggles.use_ct_qeq
    config.training.fit_hardness = config.toggles.use_ct_qeq and config.training.fit_hardness
    config.training.fit_image_plane = config.toggles.use_image_charge
    config.training.fit_reactive = config.toggles.use_reactive_grid


def main():
    args = parse_args()
    config = load_config(args.config)
    _apply_mode(config, args.mode)
    if args.max_steps is not None:
        config.training.max_steps = args.max_steps
    if args.force_weight is not None:
        config.training.force_weight = args.force_weight
    if args.learning_rate is not None:
        config.training.learning_rate = args.learning_rate
    if args.alpha_morse is not None:
        config.grid.alpha_morse = args.alpha_morse
    if args.ag_radius_scale is not None:
        config.grid.req_scale_radius["Ag"] = float(args.ag_radius_scale)
    if args.ag_energy_scale is not None:
        config.grid.req_scale_energy["Ag"] = float(args.ag_energy_scale)

    density = make_density_backend(config.density_backend, surface=config.surface, grid=config.grid).load()
    dataset_dir = Path(args.dataset_dir)
    selected_adsorbates = {name.strip() for name in args.adsorbates.split(",") if name.strip()}
    datasets = []
    for path in sorted(dataset_dir.glob("*.npz")):
        if selected_adsorbates and path.stem not in selected_adsorbates:
            continue
        poses, teacher, _ = load_pose_dataset(path)
        datasets.append((poses, teacher))
    if not datasets:
        raise ValueError("No datasets selected for fitting")

    reactive_elements = tuple(sorted({symbol for poses, _ in datasets for symbol in poses.adsorbate.symbols}))
    model = HybridGridFFModel(
        density=density,
        reactive_elements=reactive_elements,
        toggles=config.toggles,
        grid_config=config.grid,
        hybrid_config=config.hybrid_model,
        prefer_jax=args.prefer_jax,
    )
    initial_params = None
    if args.init_fit_dir:
        initial_params = _load_params(Path(args.init_fit_dir) / "fit_params.json")
    fit_t0 = time.perf_counter()
    fit_result = fit_hybrid_parameters(
        density=density,
        datasets=datasets,
        model=model,
        training=config.training,
        initial_params=initial_params,
    )
    fit_elapsed = time.perf_counter() - fit_t0

    out_dir = ensure_dir(args.out_dir)
    export_paths = export_firecore_artifacts(
        out_dir=out_dir,
        density=density,
        model=model,
        fit_result=fit_result,
        toggles=config.toggles,
        teacher_backend_id=config.teacher_backend.kind,
    )
    save_config(config, out_dir / "config_used.json")
    save_json(_params_to_dict(fit_result.params), out_dir / "fit_params.json")
    save_json(
        {
            "metrics": fit_result.metrics,
            "optimizer": fit_result.optimizer_result,
            "history": fit_result.history,
            "timing": {
                "fit_seconds_total": fit_elapsed,
            },
            "backend": {
                "prefer_jax": bool(args.prefer_jax),
                "model_use_jax": bool(model.use_jax),
                "student_backend": backend_summary(),
            },
            "mode": args.mode,
            "adsorbates": sorted(selected_adsorbates) if selected_adsorbates else "all",
            "initial_params_source": str(Path(args.init_fit_dir) / "fit_params.json") if args.init_fit_dir else "",
            "export_paths": export_paths,
        },
        out_dir / "fit_summary.json",
    )
    print(f"fit complete -> {out_dir}")


if __name__ == "__main__":
    main()
