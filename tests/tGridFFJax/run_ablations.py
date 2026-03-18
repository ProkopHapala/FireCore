#!/usr/bin/env python3
"""Fit and compare GridFF component ablations against the real teacher datasets."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.append(str(ROOT))

from pyBall.gridff_jax import load_config
from pyBall.gridff_jax.benchmarking import (
    ABLATION_CASES,
    evaluate_postfit_ablation_suite,
    evaluate_model_on_datasets,
    load_fit_params,
    plot_ablation_summary,
)
from pyBall.gridff_jax.density_backends import make_density_backend
from pyBall.gridff_jax.fit import load_pose_dataset
from pyBall.gridff_jax.hybrid_energy import HybridGridFFModel
from pyBall.gridff_jax.utils import ensure_dir, save_json


def parse_args():
    base = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description="Run GridFF ablation fits and compare their errors")
    parser.add_argument("--config", default=str(base / "ag_dataset" / "config.json"))
    parser.add_argument("--dataset-dir", default=str(base / "ag_dataset" / "datasets"))
    parser.add_argument("--out-dir", default=str(base / "ablations"))
    parser.add_argument("--fit-params", default=str(base / "ag_fit" / "fit_params.json"))
    parser.add_argument("--cases", default="plq,plq_reactive,plq_ct,plq_ct_image,full")
    parser.add_argument("--prefer-jax", dest="prefer_jax", action="store_true", help="Use the JAX student backend when available")
    parser.add_argument("--no-prefer-jax", dest="prefer_jax", action="store_false", help="Force the NumPy/SciPy fallback backend")
    parser.set_defaults(prefer_jax=True)
    return parser.parse_args()


def main():
    args = parse_args()
    config = load_config(args.config)
    density = make_density_backend(config.density_backend, surface=config.surface, grid=config.grid).load()
    datasets = []
    for path in sorted(Path(args.dataset_dir).glob("*.npz")):
        poses, teacher, _ = load_pose_dataset(path)
        datasets.append((poses, teacher))
    case_names = [item.strip() for item in args.cases.split(",") if item.strip()]
    out_dir = ensure_dir(args.out_dir)
    params = load_fit_params(args.fit_params)
    summary = evaluate_postfit_ablation_suite(
        density=density,
        datasets=datasets,
        config=config,
        params=params,
        prefer_jax=args.prefer_jax,
        case_names=case_names,
    )

    for case_name in case_names:
        model = HybridGridFFModel(
            density=density,
            reactive_elements=tuple(sorted({symbol for poses, _ in datasets for symbol in poses.adsorbate.symbols})),
            toggles=ABLATION_CASES[case_name],
            grid_config=config.grid,
            hybrid_config=config.hybrid_model,
            prefer_jax=args.prefer_jax,
        )
        case_metrics = evaluate_model_on_datasets(model, params, datasets)
        save_json(case_metrics, ensure_dir(out_dir / case_name) / "evaluation.json")
    save_json(summary, out_dir / "ablation_evaluations.json")
    plot_ablation_summary(summary, out_dir)
    print(f"ablations -> {out_dir}")


if __name__ == "__main__":
    main()
