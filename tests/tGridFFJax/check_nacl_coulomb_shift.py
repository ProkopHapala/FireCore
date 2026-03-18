#!/usr/bin/env python3
"""Check whether a one-pixel Coulomb-grid shift explains the remaining NaCl Q mismatch."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.append(str(ROOT))

from pyBall.gridff_jax.utils import ensure_dir, save_json


GRIDFF_TEST_DIR = ROOT / "tests" / "tGridFF"
COMMON_RESOURCES = ROOT / "cpp" / "common_resources"
XYZ_DIR = COMMON_RESOURCES / "xyz"
RUN_HEADLESS_SCAN = GRIDFF_TEST_DIR / "run_headless_scan.py"
DEFAULT_CANDIDATE_GRID = Path("/tmp/nacl_parity_benchmark_real2/jax_grid/Bspline_PLQd.npy")
SCAN_PREFIX = {"coulomb": "coul", "total": "total"}
SHIFT_CASES = {
    "baseline": (0, 0),
    "shift_x_plus": (0, 1),
    "shift_x_minus": (0, -1),
    "shift_y_plus": (1, 1),
    "shift_y_minus": (1, -1),
    "shift_z_plus": (2, 1),
    "shift_z_minus": (2, -1),
}


def parse_args():
    parser = argparse.ArgumentParser(description="Diagnose the NaCl Coulomb mismatch by shifting Q by one pixel")
    parser.add_argument("--surface", default="NaCl_1x1_L1")
    parser.add_argument("--molecule", default="H2O")
    parser.add_argument("--candidate-grid", default=str(DEFAULT_CANDIDATE_GRID))
    parser.add_argument("--candidate-surface-prefix", default="NaCl_1x1_L1_qshift")
    parser.add_argument("--scan-types", default="coulomb,total")
    parser.add_argument("--cases", default="baseline,shift_z_plus,shift_z_minus", help="Comma-separated shift cases to test")
    parser.add_argument("--nscan", type=int, default=200)
    parser.add_argument("--span-min", type=float, default=2.5)
    parser.add_argument("--span-max", type=float, default=15.0)
    parser.add_argument("--ld-preload", default="/usr/lib/gcc/x86_64-linux-gnu/13/libasan.so")
    parser.add_argument("--out-dir", default=str(ROOT / "tests" / "tGridFFJax" / "nacl_qshift_check"))
    return parser.parse_args()


def run_command(cmd: list[str], cwd: Path, env: dict[str, str], log_path: Path):
    t0 = time.perf_counter()
    result = subprocess.run(cmd, cwd=str(cwd), env=env, capture_output=True, text=True)
    elapsed = time.perf_counter() - t0
    log_path.write_text(result.stdout + "\n[stderr]\n" + result.stderr, encoding="utf8")
    if result.returncode != 0:
        raise subprocess.CalledProcessError(result.returncode, cmd, output=result.stdout, stderr=result.stderr)
    return elapsed


def copy_surface_grid(grid: np.ndarray, source_xyz: Path, surface_name: str):
    surface_dir = COMMON_RESOURCES / surface_name
    surface_xyz = XYZ_DIR / f"{surface_name}.xyz"
    ensure_dir(surface_dir)
    np.save(surface_dir / "Bspline_PLQd.npy", grid)
    shutil.copy2(source_xyz, surface_xyz)
    return surface_xyz


def load_scan_payload(path: Path):
    payload = np.load(path)
    return {
        "ts": np.asarray(payload["ts"], dtype=np.float64),
        "energies_norm": np.asarray(payload["energies_norm"], dtype=np.float64),
        "forces": np.asarray(payload["forces"], dtype=np.float64),
    }


def error_metrics(left, right):
    delta = np.asarray(left, dtype=np.float64) - np.asarray(right, dtype=np.float64)
    return {
        "rmse": float(np.sqrt(np.mean(delta * delta))),
        "mae": float(np.mean(np.abs(delta))),
        "max_abs": float(np.max(np.abs(delta))),
    }


def force_z_trace(forces):
    return np.sum(np.asarray(forces, dtype=np.float64)[:, :, 2], axis=1)


def plot_shift_errors(summary: dict, scan_type: str, out_path: Path):
    cases = list(summary["cases"].keys())
    energy_rmse = [summary["cases"][case][scan_type]["energy"]["rmse"] for case in cases]
    force_rmse = [summary["cases"][case][scan_type]["force_trace_fz_sum"]["rmse"] for case in cases]

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    xpos = np.arange(len(cases))

    axes[0].bar(xpos, energy_rmse, color="tab:red")
    axes[0].set_xticks(xpos, cases, rotation=45, ha="right")
    axes[0].set_title(f"{scan_type} energy RMSE")
    axes[0].set_ylabel("RMSE (eV)")
    axes[0].grid(True, axis="y", alpha=0.3)

    axes[1].bar(xpos, np.maximum(force_rmse, 1.0e-16), color="tab:orange")
    axes[1].set_xticks(xpos, cases, rotation=45, ha="right")
    axes[1].set_title(f"{scan_type} force-trace RMSE")
    axes[1].set_ylabel("RMSE (eV/A)")
    axes[1].set_yscale("log")
    axes[1].grid(True, axis="y", which="both", alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main():
    args = parse_args()
    out_dir = ensure_dir(args.out_dir)
    logs_dir = ensure_dir(out_dir / "logs")
    scans_dir = ensure_dir(out_dir / "scans")
    scan_types = tuple(item.strip() for item in args.scan_types.split(",") if item.strip())
    cases = tuple(item.strip() for item in args.cases.split(",") if item.strip())

    base_env = dict(os.environ)
    base_env["GRIDFF_SKIP_PLOTS"] = "1"
    base_env["ASAN_OPTIONS"] = "detect_leaks=0"
    base_env["LD_PRELOAD"] = args.ld_preload

    source_xyz = XYZ_DIR / f"{args.surface}.xyz"
    molecule_xyz = XYZ_DIR / f"{args.molecule}.xyz"
    candidate_grid = np.load(args.candidate_grid)

    reference_results = {}
    for scan_type in scan_types:
        ref_out = ensure_dir(scans_dir / "reference" / scan_type)
        run_command(
            [
                "python3",
                str(RUN_HEADLESS_SCAN),
                "--molecule",
                str(molecule_xyz),
                "--substrate",
                str(source_xyz),
                "--output-dir",
                str(ref_out),
                "--scan-type",
                scan_type,
                "--nscan",
                str(args.nscan),
                "--span-min",
                str(args.span_min),
                "--span-max",
                str(args.span_max),
            ],
            cwd=GRIDFF_TEST_DIR,
            env=base_env,
            log_path=logs_dir / f"reference_{scan_type}.log",
        )
        reference_results[scan_type] = load_scan_payload(ref_out / f"{args.molecule}_{SCAN_PREFIX[scan_type]}.npz")

    summary = {
        "surface": args.surface,
        "molecule": args.molecule,
        "candidate_grid": str(args.candidate_grid),
        "scan_types": list(scan_types),
        "cases": {},
    }

    for case_name in cases:
        axis, shift = SHIFT_CASES[case_name]
        shifted = np.array(candidate_grid, copy=True)
        if shift != 0:
            shifted[..., 2] = np.roll(shifted[..., 2], shift=shift, axis=axis)
        surface_name = f"{args.candidate_surface_prefix}_{case_name}"
        surface_xyz = copy_surface_grid(shifted, source_xyz, surface_name)
        case_summary = {}
        for scan_type in scan_types:
            scan_out = ensure_dir(scans_dir / case_name / scan_type)
            elapsed = run_command(
                [
                    "python3",
                    str(RUN_HEADLESS_SCAN),
                    "--molecule",
                    str(molecule_xyz),
                    "--substrate",
                    str(surface_xyz),
                    "--output-dir",
                    str(scan_out),
                    "--scan-type",
                    scan_type,
                    "--nscan",
                    str(args.nscan),
                    "--span-min",
                    str(args.span_min),
                    "--span-max",
                    str(args.span_max),
                ],
                cwd=GRIDFF_TEST_DIR,
                env=base_env,
                log_path=logs_dir / f"{case_name}_{scan_type}.log",
            )
            payload = load_scan_payload(scan_out / f"{args.molecule}_{SCAN_PREFIX[scan_type]}.npz")
            reference = reference_results[scan_type]
            case_summary[scan_type] = {
                "energy": error_metrics(reference["energies_norm"], payload["energies_norm"]),
                "force_tensor": error_metrics(reference["forces"], payload["forces"]),
                "force_trace_fz_sum": error_metrics(force_z_trace(reference["forces"]), force_z_trace(payload["forces"])),
                "scan_seconds": float(elapsed),
            }
        summary["cases"][case_name] = case_summary

    for scan_type in scan_types:
        plot_shift_errors(summary, scan_type, out_dir / f"{scan_type}_shift_errors.png")

    save_json(summary, out_dir / "summary.json")
    print(f"NaCl Coulomb shift check written to {out_dir}")


if __name__ == "__main__":
    main()
