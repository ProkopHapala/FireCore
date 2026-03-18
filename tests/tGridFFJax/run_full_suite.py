#!/usr/bin/env python3
"""Orchestrate the full real-teacher benchmark workflow for tGridFFJax."""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]


def parse_args():
    base = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description="Run the full MAD-SURF -> GridFF benchmark suite")
    parser.add_argument("--work-dir", default=str(base / "suite_run"))
    parser.add_argument("--teacher-device", default="cpu")
    parser.add_argument("--samples-per-site", type=int, default=8)
    parser.add_argument("--random-samples", type=int, default=64)
    parser.add_argument("--representative-sites-per-label", type=int, default=1)
    parser.add_argument("--max-steps", type=int, default=400)
    parser.add_argument("--force-weight", type=float, default=5.0)
    parser.add_argument("--full-sample-size", type=int, default=4096)
    parser.add_argument("--prefer-jax", dest="prefer_jax", action="store_true")
    parser.add_argument("--no-prefer-jax", dest="prefer_jax", action="store_false")
    parser.set_defaults(prefer_jax=True)
    return parser.parse_args()


def _run(cmd):
    print("RUN:", " ".join(cmd))
    subprocess.run(cmd, check=True)


def main():
    args = parse_args()
    base = Path(__file__).resolve().parent
    work = Path(args.work_dir)
    dataset_dir = work / "ag_dataset"
    fit_dir = work / "ag_fit"
    ablation_dir = work / "ablations"
    path_dir = work / "path_benchmarks"
    sanity_dir = work / "teacher_sanity"

    env_prefix = ["env", "MPLCONFIGDIR=/tmp/mpl_gridff"]
    _run(env_prefix + [
        "python3", str(base / "run_madsurf_teacher_sanity.py"),
        "--out-dir", str(sanity_dir),
        "--device", args.teacher_device,
        "--full-sample-size", str(args.full_sample_size),
    ])
    _run(env_prefix + [
        "python3", str(base / "build_ag_dataset.py"),
        "--teacher", "madsurf",
        "--device", args.teacher_device,
        "--out-dir", str(dataset_dir),
        "--samples-per-site", str(args.samples_per_site),
        "--random-samples", str(args.random_samples),
        "--representative-sites-per-label", str(args.representative_sites_per_label),
    ])
    _run(env_prefix + [
        "python3", str(base / "fit_hybrid_gridff.py"),
        "--config", str(dataset_dir / "config.json"),
        "--dataset-dir", str(dataset_dir / "datasets"),
        "--out-dir", str(fit_dir),
        "--mode", "plq",
        "--max-steps", str(args.max_steps),
        "--force-weight", str(args.force_weight),
        "--prefer-jax" if args.prefer_jax else "--no-prefer-jax",
    ])
    _run(env_prefix + [
        "python3", str(base / "validate_hybrid_gridff.py"),
        "--config", str(fit_dir / "config_used.json"),
        "--dataset-dir", str(dataset_dir / "datasets"),
        "--fit-dir", str(fit_dir),
        "--out-dir", str(fit_dir / "validation"),
        "--mode", "plq",
        "--prefer-jax" if args.prefer_jax else "--no-prefer-jax",
    ])
    _run(env_prefix + [
        "python3", str(base / "run_ablations.py"),
        "--config", str(dataset_dir / "config.json"),
        "--dataset-dir", str(dataset_dir / "datasets"),
        "--out-dir", str(ablation_dir),
        "--fit-params", str(fit_dir / "fit_params.json"),
        "--prefer-jax" if args.prefer_jax else "--no-prefer-jax",
    ])
    _run(env_prefix + [
        "python3", str(base / "benchmark_pose_paths.py"),
        "--config", str(fit_dir / "config_used.json"),
        "--fit-params", str(fit_dir / "fit_params.json"),
        "--out-dir", str(path_dir),
        "--prefer-jax" if args.prefer_jax else "--no-prefer-jax",
    ])
    print(f"full suite -> {work}")


if __name__ == "__main__":
    main()
