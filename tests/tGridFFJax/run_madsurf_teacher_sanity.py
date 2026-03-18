#!/usr/bin/env python3
"""Evaluate the local MAD-SURF model against the supplied reference extxyz dataset."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.append(str(ROOT))

from pyBall.gridff_jax.benchmarking import run_teacher_sanity_suite
from pyBall.gridff_jax.utils import ensure_dir, save_json


def parse_args():
    base = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description="Run MAD-SURF sanity metrics on the provided extxyz datasets")
    parser.add_argument("--model-path", default=str(base / "mad-surf_data" / "models" / "full_dataset_config_weights" / "MACE_model.model"))
    parser.add_argument("--small-dataset-path", default=str(base / "mad-surf_data" / "dataset" / "test_small_dataset_std.extxyz"))
    parser.add_argument("--full-dataset-path", default=str(base / "mad-surf_data" / "full_train_test_std_config_types.extxyz"))
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--full-sample-size", type=int, default=4096)
    parser.add_argument("--out-dir", default=str(base / "teacher_sanity"))
    return parser.parse_args()


def main():
    args = parse_args()
    out_dir = ensure_dir(args.out_dir)
    metrics = run_teacher_sanity_suite(
        model_path=args.model_path,
        small_extxyz=args.small_dataset_path,
        full_extxyz=args.full_dataset_path,
        full_sample_size=args.full_sample_size,
        device=args.device,
    )
    save_json(metrics, out_dir / "teacher_sanity.json")
    print(f"teacher sanity -> {out_dir}")


if __name__ == "__main__":
    main()
