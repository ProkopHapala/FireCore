#!/usr/bin/env python3
"""Run single-atom H/C/N/O Ag z-scan diagnostics with the strict-PLQ benchmark."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[2]
BENCHMARK = Path(__file__).resolve().with_name("benchmark_ag_zscan.py")
DEFAULT_MODEL = Path(__file__).resolve().parent / "mad-surf_data" / "models" / "full_dataset_config_weights" / "MACE_model.model"
DEFAULT_CHGCAR = "/home/niel/git/ORR_HER_Ag_Colab/results/Ag_ORR_HER/slab_clean/final_scf_12x12x1/CHGCAR"
DEFAULT_LOCPOT = "/home/niel/git/ORR_HER_Ag_Colab/results/Ag_ORR_HER/slab_clean/workfunc_12x12x1/LOCPOT"


def parse_args():
    parser = argparse.ArgumentParser(description="Run strict-PLQ Ag atomic anchor diagnostics")
    parser.add_argument("--out-dir", default=str(ROOT / "tests/tGridFFJax/ag_anchor_diagnostics"))
    parser.add_argument("--adsorbates", default="H,C,N,O")
    parser.add_argument("--chgcar", default=DEFAULT_CHGCAR)
    parser.add_argument("--locpot", default=DEFAULT_LOCPOT)
    parser.add_argument("--model-path", default=str(DEFAULT_MODEL))
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--z-min", type=float, default=2.0)
    parser.add_argument("--stress-z-min", type=float, default=1.8)
    parser.add_argument("--z-max", type=float, default=5.6)
    parser.add_argument("--eval-z-points", type=int, default=121)
    parser.add_argument("--train-low-count", type=int, default=18)
    parser.add_argument("--train-high-count", type=int, default=12)
    parser.add_argument("--z-dense-cutoff", type=float, default=2.8)
    parser.add_argument("--max-steps", type=int, default=400)
    parser.add_argument("--force-weight", type=float, default=10.0)
    parser.add_argument("--learning-rate", type=float, default=1.0e-2)
    parser.add_argument("--teacher-chunk-size", type=int, default=64)
    parser.add_argument("--student-chunk-size", type=int, default=64)
    parser.add_argument("--prefer-jax", dest="prefer_jax", action="store_true")
    parser.add_argument("--no-prefer-jax", dest="prefer_jax", action="store_false")
    parser.set_defaults(prefer_jax=True)
    return parser.parse_args()


def _run_one(args, adsorbate: str, out_dir: Path):
    summary_path = out_dir / "zscan_summary.json"
    if summary_path.exists():
        with summary_path.open("r", encoding="utf8") as fin:
            return json.load(fin)
    cmd = [
        sys.executable,
        str(BENCHMARK),
        "--adsorbate-name",
        adsorbate,
        "--out-dir",
        str(out_dir),
        "--chgcar",
        args.chgcar,
        "--locpot",
        args.locpot,
        "--model-path",
        args.model_path,
        "--device",
        args.device,
        "--z-min",
        str(args.z_min),
        "--stress-z-min",
        str(args.stress_z_min),
        "--z-max",
        str(args.z_max),
        "--eval-z-points",
        str(args.eval_z_points),
        "--train-low-count",
        str(args.train_low_count),
        "--train-high-count",
        str(args.train_high_count),
        "--z-dense-cutoff",
        str(args.z_dense_cutoff),
        "--max-steps",
        str(args.max_steps),
        "--force-weight",
        str(args.force_weight),
        "--learning-rate",
        str(args.learning_rate),
        "--teacher-chunk-size",
        str(args.teacher_chunk_size),
        "--student-chunk-size",
        str(args.student_chunk_size),
        "--no-fit-static-charge",
    ]
    if args.prefer_jax:
        cmd.append("--prefer-jax")
    else:
        cmd.append("--no-prefer-jax")
    env = dict(os.environ)
    env["MPLBACKEND"] = "Agg"
    subprocess.run(cmd, check=True, env=env)
    with summary_path.open("r", encoding="utf8") as fin:
        return json.load(fin)


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    adsorbates = [item.strip() for item in args.adsorbates.split(",") if item.strip()]
    summaries = {}
    for adsorbate in adsorbates:
        ads_dir = out_dir / adsorbate
        ads_dir.mkdir(parents=True, exist_ok=True)
        summaries[adsorbate] = _run_one(args, adsorbate, ads_dir)
    aggregate = {}
    for adsorbate, summary in summaries.items():
        aggregate[adsorbate] = {
            "primary_energy_rmse": summary["windows"]["primary"]["all"]["energy"]["rmse"],
            "primary_force_rmse": summary["windows"]["primary"]["all"]["force"]["rmse"],
            "stress_energy_rmse": summary["windows"]["stress"]["all"]["energy"]["rmse"],
            "stress_force_rmse": summary["windows"]["stress"]["all"]["force"]["rmse"],
            "constraint_limited": summary["acceptance"]["constraint_limited"],
        }
    result = {
        "adsorbates": adsorbates,
        "aggregate": aggregate,
        "runs": {name: str((out_dir / name).resolve()) for name in adsorbates},
    }
    with (out_dir / "anchor_diagnostics_summary.json").open("w", encoding="utf8") as fout:
        json.dump(result, fout, indent=2)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
