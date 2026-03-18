#!/usr/bin/env python3
"""Real one-to-one NaCl parity benchmark between OpenCL GridFF and JAX GridFF."""

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

import jax


GRIDFF_TEST_DIR = ROOT / "tests" / "tGridFF"
COMMON_RESOURCES = ROOT / "cpp" / "common_resources"
XYZ_DIR = COMMON_RESOURCES / "xyz"
RUN_HEADLESS_SCAN = GRIDFF_TEST_DIR / "run_headless_scan.py"
RUN_NACL_PARITY = ROOT / "tests" / "tGridFFJax" / "run_nacl_parity.py"
DEFAULT_SURFACE = "NaCl_1x1_L1"
DEFAULT_MOLECULE = "H2O"
DEFAULT_SCAN_TYPES = ("total", "morse", "pauli", "london", "coulomb")
SCAN_PREFIX = {
    "total": "total",
    "morse": "morse",
    "pauli": "pauli",
    "london": "london",
    "coulomb": "coul",
}


def parse_args():
    parser = argparse.ArgumentParser(description="Benchmark OpenCL vs JAX NaCl GridFF parity")
    parser.add_argument("--surface", default=DEFAULT_SURFACE, help="Surface base name without .xyz")
    parser.add_argument("--molecule", default=DEFAULT_MOLECULE, help="Molecule base name without .xyz")
    parser.add_argument("--desired-voxel", type=float, default=0.1, help="Grid voxel size in Angstrom")
    parser.add_argument("--nscan", type=int, default=200, help="Number of scan points")
    parser.add_argument("--span-min", type=float, default=2.5, help="Start of z-scan")
    parser.add_argument("--span-max", type=float, default=15.0, help="End of z-scan")
    parser.add_argument("--candidate-surface-name", default="NaCl_1x1_L1_jax_bench", help="Temporary surface name for the JAX-generated grid")
    parser.add_argument("--out-dir", default=str(ROOT / "tests" / "tGridFFJax" / "nacl_parity_benchmark"), help="Output directory")
    parser.add_argument("--ld-preload", default="/usr/lib/gcc/x86_64-linux-gnu/13/libasan.so", help="ASan preload used by MMFF scan binaries")
    parser.add_argument("--scan-types", default=",".join(DEFAULT_SCAN_TYPES), help="Comma-separated scan types")
    return parser.parse_args()


def run_command(cmd: list[str], cwd: Path, env: dict[str, str], log_path: Path):
    t0 = time.perf_counter()
    result = subprocess.run(
        cmd,
        cwd=str(cwd),
        env=env,
        capture_output=True,
        text=True,
    )
    elapsed = time.perf_counter() - t0
    log_path.write_text(result.stdout + "\n[stderr]\n" + result.stderr, encoding="utf8")
    if result.returncode != 0:
        raise subprocess.CalledProcessError(result.returncode, cmd, output=result.stdout, stderr=result.stderr)
    return elapsed


def copy_candidate_surface(source_grid: Path, source_xyz: Path, candidate_surface_name: str):
    surface_dir = COMMON_RESOURCES / candidate_surface_name
    xyz_path = XYZ_DIR / f"{candidate_surface_name}.xyz"
    ensure_dir(surface_dir)
    shutil.copy2(source_grid, surface_dir / "Bspline_PLQd.npy")
    shutil.copy2(source_xyz, xyz_path)
    return xyz_path


def load_scan_payload(path: Path):
    payload = np.load(path)
    return {
        "ts": np.asarray(payload["ts"], dtype=np.float64),
        "energies": np.asarray(payload["energies"], dtype=np.float64),
        "energies_norm": np.asarray(payload["energies_norm"], dtype=np.float64),
        "forces": np.asarray(payload["forces"], dtype=np.float64),
        "positions": np.asarray(payload["positions"], dtype=np.float64),
        "scan_positions": np.asarray(payload["scan_positions"], dtype=np.float64),
    }


def error_metrics(left, right):
    delta = np.asarray(left, dtype=np.float64) - np.asarray(right, dtype=np.float64)
    return {
        "rmse": float(np.sqrt(np.mean(delta * delta))),
        "mae": float(np.mean(np.abs(delta))),
        "max_abs": float(np.max(np.abs(delta))),
    }


def force_z_trace(forces):
    force_array = np.asarray(forces, dtype=np.float64)
    if force_array.ndim != 3:
        raise ValueError(f"Expected forces with shape (nscan, natom, 3), got {force_array.shape}")
    return np.sum(force_array[:, :, 2], axis=1)


def plot_scan_comparison(scan_type: str, ref_payload, cand_payload, ref_timing: dict, cand_timing: dict, out_path: Path):
    ts = ref_payload["ts"]
    energy_ref = ref_payload["energies_norm"]
    energy_cand = cand_payload["energies_norm"]
    energy_abs = np.abs(energy_cand - energy_ref)

    force_ref = force_z_trace(ref_payload["forces"])
    force_cand = force_z_trace(cand_payload["forces"])
    force_abs = np.abs(force_cand - force_ref)

    eps = 1.0e-12
    fig, axes = plt.subplots(2, 3, figsize=(18, 8), sharex="col")

    axes[0, 0].plot(ts, energy_ref, label="OpenCL", lw=2.0)
    axes[0, 0].plot(ts, energy_cand, label="JAX", lw=1.4)
    axes[0, 0].set_ylabel("Energy (eV)")
    axes[0, 0].set_title(f"{scan_type} energy")
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)

    axes[0, 1].plot(ts, energy_cand - energy_ref, color="tab:red", lw=1.5)
    axes[0, 1].set_title("Energy difference (JAX - OpenCL)")
    axes[0, 1].grid(True, alpha=0.3)

    axes[0, 2].semilogy(ts, np.maximum(energy_abs, eps), color="tab:red", lw=1.5)
    axes[0, 2].set_title("Energy |difference| log scale")
    axes[0, 2].grid(True, which="both", alpha=0.3)

    axes[1, 0].plot(ts, force_ref, label="OpenCL", lw=2.0)
    axes[1, 0].plot(ts, force_cand, label="JAX", lw=1.4)
    axes[1, 0].set_xlabel("z scan coordinate (A)")
    axes[1, 0].set_ylabel("Sum Fz (eV/A)")
    axes[1, 0].set_title("Force trace")
    axes[1, 0].grid(True, alpha=0.3)

    axes[1, 1].plot(ts, force_cand - force_ref, color="tab:orange", lw=1.5)
    axes[1, 1].set_xlabel("z scan coordinate (A)")
    axes[1, 1].set_title("Force difference (JAX - OpenCL)")
    axes[1, 1].grid(True, alpha=0.3)

    axes[1, 2].semilogy(ts, np.maximum(force_abs, eps), color="tab:orange", lw=1.5)
    axes[1, 2].set_xlabel("z scan coordinate (A)")
    axes[1, 2].set_title("Force |difference| log scale")
    axes[1, 2].grid(True, which="both", alpha=0.3)

    fig.suptitle(
        (
            f"NaCl parity benchmark: {scan_type}\n"
            f"OpenCL scan {ref_timing['scan']:.3f}s, JAX scan {cand_timing['scan']:.3f}s"
        ),
        fontsize=13,
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_timing(summary: dict, out_path: Path):
    labels = ["grid_opencl", "grid_jax"]
    values = [
        summary["timings_s"]["grid_generation_opencl"],
        summary["timings_s"]["grid_generation_jax"],
    ]
    for scan_type in summary["scan_types"]:
        labels.append(f"{scan_type}_scan_opencl")
        values.append(summary["scan_results"][scan_type]["timing_s"]["opencl_scan"])
        labels.append(f"{scan_type}_scan_jax")
        values.append(summary["scan_results"][scan_type]["timing_s"]["jax_scan"])

    fig, ax = plt.subplots(figsize=(10, 5))
    ypos = np.arange(len(labels))
    ax.barh(ypos, values, color=["tab:blue", "tab:green"] + ["tab:gray", "tab:purple"] * len(summary["scan_types"]))
    ax.set_yticks(ypos, labels)
    ax.set_xlabel("Wall time (s)")
    ax.set_title("NaCl parity benchmark timings")
    ax.grid(True, axis="x", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main():
    args = parse_args()
    out_dir = ensure_dir(args.out_dir)
    logs_dir = ensure_dir(out_dir / "logs")
    scans_dir = ensure_dir(out_dir / "scans")
    scan_types = tuple(item.strip() for item in args.scan_types.split(",") if item.strip())

    surface_xyz = XYZ_DIR / f"{args.surface}.xyz"
    molecule_xyz = XYZ_DIR / f"{args.molecule}.xyz"
    reference_grid = COMMON_RESOURCES / args.surface / "Bspline_PLQd.npy"
    candidate_grid_dir = out_dir / "jax_grid"
    candidate_grid = candidate_grid_dir / "Bspline_PLQd.npy"

    base_env = dict(os.environ)
    base_env["JAX_ENABLE_X64"] = "1"
    base_env["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"

    opencl_time = run_command(
        [
            "python3",
            "generate_grid.py",
            "--fname",
            f"{args.surface}.xyz",
            "--desired_voxel",
            str(args.desired_voxel),
        ],
        cwd=GRIDFF_TEST_DIR,
        env=base_env,
        log_path=logs_dir / "opencl_generate.log",
    )

    jax_time = run_command(
        [
            "python3",
            str(RUN_NACL_PARITY),
            "--xyz",
            str(surface_xyz),
            "--desired-voxel",
            str(args.desired_voxel),
            "--reference-grid",
            str(reference_grid),
            "--out-dir",
            str(candidate_grid_dir),
            "--prefer-jax",
        ],
        cwd=ROOT,
        env=base_env,
        log_path=logs_dir / "jax_generate.log",
    )

    candidate_surface_xyz = copy_candidate_surface(candidate_grid, surface_xyz, args.candidate_surface_name)

    scan_env = dict(base_env)
    scan_env["GRIDFF_SKIP_PLOTS"] = "1"
    scan_env["ASAN_OPTIONS"] = "detect_leaks=0"
    scan_env["LD_PRELOAD"] = args.ld_preload

    summary = {
        "surface": args.surface,
        "molecule": args.molecule,
        "scan_types": list(scan_types),
        "grid_generation": {
            "reference_grid": str(reference_grid),
            "candidate_grid": str(candidate_grid),
            "candidate_surface_xyz": str(candidate_surface_xyz),
        },
        "jax_runtime": {
            "default_backend": jax.default_backend(),
            "devices": [str(device) for device in jax.devices()],
        },
        "timings_s": {
            "grid_generation_opencl": float(opencl_time),
            "grid_generation_jax": float(jax_time),
        },
        "scan_results": {},
    }

    parity_summary = json.loads((candidate_grid_dir / "parity_summary.json").read_text(encoding="utf8"))
    summary["grid_channel_errors"] = parity_summary.get("channel_errors", {})
    summary["grid_builder_metadata"] = parity_summary.get("builder_metadata", {})

    for scan_type in scan_types:
        ref_out = ensure_dir(scans_dir / "opencl" / scan_type)
        cand_out = ensure_dir(scans_dir / "jax" / scan_type)

        ref_time = run_command(
            [
                "python3",
                str(RUN_HEADLESS_SCAN),
                "--molecule",
                str(molecule_xyz),
                "--substrate",
                str(surface_xyz),
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
            env=scan_env,
            log_path=logs_dir / f"scan_opencl_{scan_type}.log",
        )
        cand_time = run_command(
            [
                "python3",
                str(RUN_HEADLESS_SCAN),
                "--molecule",
                str(molecule_xyz),
                "--substrate",
                str(candidate_surface_xyz),
                "--output-dir",
                str(cand_out),
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
            env=scan_env,
            log_path=logs_dir / f"scan_jax_{scan_type}.log",
        )

        prefix = SCAN_PREFIX[scan_type]
        ref_npz = ref_out / f"{args.molecule}_{prefix}.npz"
        cand_npz = cand_out / f"{args.molecule}_{prefix}.npz"
        ref_metrics = json.loads((ref_out / f"{args.molecule}_{prefix}_metrics.json").read_text(encoding="utf8"))
        cand_metrics = json.loads((cand_out / f"{args.molecule}_{prefix}_metrics.json").read_text(encoding="utf8"))

        ref_payload = load_scan_payload(ref_npz)
        cand_payload = load_scan_payload(cand_npz)

        energy_metrics = error_metrics(ref_payload["energies_norm"], cand_payload["energies_norm"])
        force_metrics = error_metrics(ref_payload["forces"], cand_payload["forces"])
        force_trace_metrics = error_metrics(force_z_trace(ref_payload["forces"]), force_z_trace(cand_payload["forces"]))

        plot_scan_comparison(
            scan_type=scan_type,
            ref_payload=ref_payload,
            cand_payload=cand_payload,
            ref_timing=ref_metrics["timing_s"],
            cand_timing=cand_metrics["timing_s"],
            out_path=out_dir / f"{scan_type}_comparison.png",
        )

        summary["scan_results"][scan_type] = {
            "energy_error": energy_metrics,
            "force_error_tensor": force_metrics,
            "force_error_trace_fz_sum": force_trace_metrics,
            "timing_s": {
                "opencl_scan": float(ref_time),
                "jax_scan": float(cand_time),
                "opencl_scan_internal": float(ref_metrics["timing_s"]["scan"]),
                "jax_scan_internal": float(cand_metrics["timing_s"]["scan"]),
                "opencl_init_internal": float(ref_metrics["timing_s"]["init"]),
                "jax_init_internal": float(cand_metrics["timing_s"]["init"]),
            },
        }

    plot_timing(summary, out_dir / "timings.png")
    save_json(summary, out_dir / "summary.json")
    print(f"NaCl parity benchmark written to {out_dir}")


if __name__ == "__main__":
    main()
