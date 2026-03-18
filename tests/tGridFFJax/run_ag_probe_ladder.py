#!/usr/bin/env python3
"""Run the Ag strict-PLQ probe ladder from single atoms up to small molecules."""

from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.append(str(ROOT))


def parse_args():
    base = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description="Run the Ag strict-PLQ probe ladder")
    parser.add_argument("--work-dir", default=str(base / "ag_probe_ladder"))
    parser.add_argument("--teacher-device", default="cpu")
    parser.add_argument("--adsorbates", default="H,C,O,CO,H2O")
    parser.add_argument("--grid-shape", default="", help="Optional density grid override as nx,ny,nz")
    parser.add_argument("--samples-per-site", type=int, default=8)
    parser.add_argument("--systematic-orientations", type=int, default=1)
    parser.add_argument("--random-samples", type=int, default=32)
    parser.add_argument("--representative-sites-per-label", type=int, default=1)
    parser.add_argument("--z-bias-power", type=float, default=2.0)
    parser.add_argument("--jitter-uv", type=float, default=0.03)
    parser.add_argument("--max-steps", type=int, default=200)
    parser.add_argument("--force-weight", type=float, default=5.0)
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
    dataset_dir = work / "dataset"
    fit_dir = work / "fit"
    validation_dir = fit_dir / "validation"

    env_prefix = ["env", "MPLCONFIGDIR=/tmp/mpl_gridff"]
    _run(env_prefix + [
        "python3", str(base / "build_ag_dataset.py"),
        "--teacher", "madsurf",
        "--device", args.teacher_device,
        "--out-dir", str(dataset_dir),
        "--adsorbates", args.adsorbates,
        "--grid-shape", args.grid_shape,
        "--samples-per-site", str(args.samples_per_site),
        "--systematic-orientations", str(args.systematic_orientations),
        "--random-samples", str(args.random_samples),
        "--representative-sites-per-label", str(args.representative_sites_per_label),
        "--z-bias-power", str(args.z_bias_power),
        "--jitter-uv", str(args.jitter_uv),
    ])
    _run(env_prefix + [
        "python3", str(base / "fit_hybrid_gridff.py"),
        "--config", str(dataset_dir / "config.json"),
        "--dataset-dir", str(dataset_dir / "datasets"),
        "--out-dir", str(fit_dir),
        "--mode", "plq",
        "--max-steps", str(args.max_steps),
        "--force-weight", str(args.force_weight),
        "--adsorbates", args.adsorbates,
        "--prefer-jax" if args.prefer_jax else "--no-prefer-jax",
    ])
    _run(env_prefix + [
        "python3", str(base / "validate_hybrid_gridff.py"),
        "--config", str(fit_dir / "config_used.json"),
        "--dataset-dir", str(dataset_dir / "datasets"),
        "--fit-dir", str(fit_dir),
        "--out-dir", str(validation_dir),
        "--mode", "plq",
        "--adsorbates", args.adsorbates,
        "--prefer-jax" if args.prefer_jax else "--no-prefer-jax",
    ])
    print(f"ag probe ladder -> {work}")


if __name__ == "__main__":
    main()
