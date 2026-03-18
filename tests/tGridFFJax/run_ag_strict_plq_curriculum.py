#!/usr/bin/env python3
"""Run a strict-PLQ Ag curriculum: atomic probes first, then molecules."""

from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.append(str(ROOT))


def parse_args():
    base = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description="Run the strict-PLQ Ag curriculum benchmark")
    parser.add_argument("--work-dir", default=str(base / "ag_strict_plq_curriculum"))
    parser.add_argument("--teacher-device", default="cpu")
    parser.add_argument("--grid-shape", default="", help="Optional density grid override as nx,ny,nz")
    parser.add_argument("--probe-samples-per-site", type=int, default=12)
    parser.add_argument("--probe-random-samples", type=int, default=48)
    parser.add_argument("--probe-max-steps", type=int, default=250)
    parser.add_argument("--probe-force-weight", type=float, default=20.0)
    parser.add_argument("--molecule-samples-per-site", type=int, default=10)
    parser.add_argument("--molecule-systematic-orientations", type=int, default=4)
    parser.add_argument("--molecule-random-samples", type=int, default=128)
    parser.add_argument("--molecule-max-steps", type=int, default=400)
    parser.add_argument("--molecule-force-weight", type=float, default=25.0)
    parser.add_argument("--learning-rate", type=float, default=1.0e-2)
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

    probe_dataset_dir = work / "probe_dataset"
    probe_fit_dir = work / "probe_fit"
    probe_validation_dir = probe_fit_dir / "validation"

    molecule_dataset_dir = work / "molecule_dataset"
    molecule_fit_dir = work / "molecule_fit"
    molecule_validation_dir = molecule_fit_dir / "validation"

    env_prefix = ["env", "MPLCONFIGDIR=/tmp/mpl_gridff"]
    jax_flag = "--prefer-jax" if args.prefer_jax else "--no-prefer-jax"

    _run(env_prefix + [
        "python3", str(base / "build_ag_dataset.py"),
        "--teacher", "madsurf",
        "--device", args.teacher_device,
        "--out-dir", str(probe_dataset_dir),
        "--adsorbates", "H,C,O",
        "--grid-shape", args.grid_shape,
        "--samples-per-site", str(args.probe_samples_per_site),
        "--systematic-orientations", "1",
        "--random-samples", str(args.probe_random_samples),
        "--representative-sites-per-label", "1",
        "--z-bias-power", "2.5",
        "--jitter-uv", "0.02",
    ])
    _run(env_prefix + [
        "python3", str(base / "fit_hybrid_gridff.py"),
        "--config", str(probe_dataset_dir / "config.json"),
        "--dataset-dir", str(probe_dataset_dir / "datasets"),
        "--out-dir", str(probe_fit_dir),
        "--mode", "plq",
        "--max-steps", str(args.probe_max_steps),
        "--force-weight", str(args.probe_force_weight),
        "--learning-rate", str(args.learning_rate),
        "--adsorbates", "H,C,O",
        jax_flag,
    ])
    _run(env_prefix + [
        "python3", str(base / "validate_hybrid_gridff.py"),
        "--config", str(probe_fit_dir / "config_used.json"),
        "--dataset-dir", str(probe_dataset_dir / "datasets"),
        "--fit-dir", str(probe_fit_dir),
        "--out-dir", str(probe_validation_dir),
        "--mode", "plq",
        "--adsorbates", "H,C,O",
        jax_flag,
    ])

    _run(env_prefix + [
        "python3", str(base / "build_ag_dataset.py"),
        "--teacher", "madsurf",
        "--device", args.teacher_device,
        "--out-dir", str(molecule_dataset_dir),
        "--adsorbates", "H,CO,H2O",
        "--grid-shape", args.grid_shape,
        "--samples-per-site", str(args.molecule_samples_per_site),
        "--systematic-orientations", str(args.molecule_systematic_orientations),
        "--random-samples", str(args.molecule_random_samples),
        "--representative-sites-per-label", "1",
        "--z-bias-power", "2.5",
        "--jitter-uv", "0.02",
    ])
    _run(env_prefix + [
        "python3", str(base / "fit_hybrid_gridff.py"),
        "--config", str(molecule_dataset_dir / "config.json"),
        "--dataset-dir", str(molecule_dataset_dir / "datasets"),
        "--out-dir", str(molecule_fit_dir),
        "--mode", "plq",
        "--max-steps", str(args.molecule_max_steps),
        "--force-weight", str(args.molecule_force_weight),
        "--learning-rate", str(args.learning_rate),
        "--adsorbates", "H,CO,H2O",
        "--init-fit-dir", str(probe_fit_dir),
        jax_flag,
    ])
    _run(env_prefix + [
        "python3", str(base / "validate_hybrid_gridff.py"),
        "--config", str(molecule_fit_dir / "config_used.json"),
        "--dataset-dir", str(molecule_dataset_dir / "datasets"),
        "--fit-dir", str(molecule_fit_dir),
        "--out-dir", str(molecule_validation_dir),
        "--mode", "plq",
        "--adsorbates", "H,CO,H2O",
        jax_flag,
    ])

    print(f"strict PLQ curriculum -> {work}")


if __name__ == "__main__":
    main()
