#!/usr/bin/env python3
"""Sweep London damping parameters efficiently: run MACE teacher ONCE, then
rebuild grids and refit for each (d0, width) combination.

This avoids the bottleneck of running MACE N times when only the grid
builder changes.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

# Limit thread contention
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("MKL_NUM_THREADS", "4")
os.environ.setdefault("JAX_PLATFORMS", "cpu")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.append(str(ROOT))

from pyBall.gridff_jax import RunConfig, save_config
from pyBall.gridff_jax.density_backends import make_density_backend
from pyBall.gridff_jax.fit import fit_hybrid_parameters
from pyBall.gridff_jax.hybrid_energy import HybridGridFFModel
from pyBall.gridff_jax.interfaces import PoseBatch, TeacherResult
from pyBall.gridff_jax.pose_sampling import (
    get_adsorbate, infer_reference_planes, infer_surface_sites, transform_adsorbate,
)
from pyBall.gridff_jax.teacher_backends import make_teacher_backend
from pyBall.gridff_jax.utils import ensure_dir, save_json
from pyBall.gridff_jax.validation import compute_error_metrics

DEFAULT_CHGCAR = "/home/niel/git/ORR_HER_Ag_Colab/results/Ag_ORR_HER/slab_clean/final_scf_12x12x1/CHGCAR"
DEFAULT_LOCPOT = "/home/niel/git/ORR_HER_Ag_Colab/results/Ag_ORR_HER/slab_clean/workfunc_12x12x1/LOCPOT"
DEFAULT_MODEL  = str(Path(__file__).resolve().parent / "mad-surf_data" / "models" / "full_dataset_config_weights" / "MACE_model.model")


def parse_args():
    p = argparse.ArgumentParser(description="Sweep London damping d0/width on Ag(111) z-scan")
    p.add_argument("--chgcar", default=DEFAULT_CHGCAR)
    p.add_argument("--locpot", default=DEFAULT_LOCPOT)
    p.add_argument("--model-path", default=DEFAULT_MODEL)
    p.add_argument("--device", default="cpu")
    p.add_argument("--out-dir", default=str(Path(__file__).resolve().parent / "damping_sweep"))
    p.add_argument("--sites", default="top,bridge,hollow")
    p.add_argument("--z-min", type=float, default=2.0)
    p.add_argument("--z-max", type=float, default=8.0)
    p.add_argument("--stress-z-min", type=float, default=1.8)
    p.add_argument("--eval-z-points", type=int, default=63)
    p.add_argument("--train-low-count", type=int, default=12)
    p.add_argument("--train-high-count", type=int, default=8)
    p.add_argument("--z-dense-cutoff", type=float, default=3.5)
    p.add_argument("--max-steps", type=int, default=600)
    p.add_argument("--force-weight", type=float, default=10.0)
    p.add_argument("--learning-rate", type=float, default=0.01)
    # Damping sweep ranges
    p.add_argument("--d0-values", default="0.0,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0",
                   help="Comma-separated d0 values to sweep (0.0 = no damping = baseline)")
    p.add_argument("--width-values", default="0.5",
                   help="Comma-separated width values to sweep")
    return p.parse_args()


def _rotation_quaternion(axis, angle_rad):
    axis = np.asarray(axis, dtype=float)
    norm = np.linalg.norm(axis)
    if norm < 1e-12:
        return np.array([1.0, 0.0, 0.0, 0.0])
    axis = axis / norm
    half = 0.5 * float(angle_rad)
    return np.array([np.cos(half), axis[0]*np.sin(half), axis[1]*np.sin(half), axis[2]*np.sin(half)])


def _site_uv_lookup(density, requested_labels):
    requested = [l.strip() for l in requested_labels.split(",") if l.strip()]
    lookup = {}
    for site in infer_surface_sites(density):
        if site.label in requested and site.label not in lookup:
            lookup[site.label] = np.asarray(site.uv, dtype=float)
    missing = [l for l in requested if l not in lookup]
    if missing:
        raise ValueError(f"Sites not found: {missing}")
    return {l: lookup[l] for l in requested}


def _build_zscan_poses(density, adsorbate, site_uvs, z_values, quaternion):
    planes = infer_reference_planes(density)
    z_ref = float(planes["z_ref"])
    z_image = float(planes["z_image"])
    positions, pose_params, site_labels = [], [], []
    for label, uv in site_uvs.items():
        for z_h in z_values:
            positions.append(transform_adsorbate(adsorbate, density, uv, float(z_h), quaternion, z_ref))
            pose_params.append(np.concatenate([uv, [float(z_h)], quaternion]))
            site_labels.append(label)
    return PoseBatch(
        adsorbate=adsorbate,
        positions=np.asarray(positions),
        pose_params=np.asarray(pose_params),
        site_labels=site_labels,
        metadata={"z_ref": z_ref, "z_image": z_image, "scan_z_values": z_values.tolist()},
    )


def _slice_pose_batch(poses, indices):
    idx = np.asarray(indices, dtype=int)
    return PoseBatch(
        adsorbate=poses.adsorbate,
        positions=poses.positions[idx].copy(),
        pose_params=poses.pose_params[idx].copy(),
        site_labels=[poses.site_labels[int(i)] for i in idx],
        metadata=dict(poses.metadata),
    )


def _slice_teacher(teacher, indices):
    idx = np.asarray(indices, dtype=int)
    return TeacherResult(
        energies=teacher.energies[idx].copy(),
        forces=teacher.forces[idx].copy(),
        metadata=dict(teacher.metadata),
    )


def _select_even(indices, count):
    indices = np.asarray(indices, dtype=int)
    if len(indices) <= count:
        return indices
    coords = np.linspace(0, len(indices)-1, count)
    return np.unique(indices[np.round(coords).astype(int)])


def _train_indices(n_sites, nz, z_values, z_min, z_dense_cutoff, z_max, n_low, n_high):
    z = np.asarray(z_values)
    low = np.where((z >= z_min) & (z < z_dense_cutoff))[0]
    high = np.where((z >= z_dense_cutoff) & (z <= z_max))[0]
    local = np.unique(np.concatenate([_select_even(low, n_low), _select_even(high, n_high)]))
    indices = []
    for i in range(n_sites):
        indices.extend((i * nz + local).tolist())
    return np.asarray(indices, dtype=int)


def _evaluate_teacher_chunked(teacher_backend, density, poses, chunk=64):
    t0 = time.perf_counter()
    if len(poses.positions) <= chunk:
        result = teacher_backend.evaluate_batch(density, poses)
        return result, time.perf_counter() - t0
    energies, forces = [], []
    for start in range(0, len(poses.positions), chunk):
        stop = min(start + chunk, len(poses.positions))
        partial = teacher_backend.evaluate_batch(density, _slice_pose_batch(poses, np.arange(start, stop)))
        energies.append(np.asarray(partial.energies))
        forces.append(np.asarray(partial.forces))
    elapsed = time.perf_counter() - t0
    return TeacherResult(energies=np.concatenate(energies), forces=np.concatenate(forces), metadata={}), elapsed


def _metrics(teacher_e, student_e, teacher_f, student_f, mask):
    mask = np.asarray(mask, dtype=bool)
    if mask.sum() == 0:
        return {"n": 0, "e_rmse": 0, "e_mae": 0, "f_rmse": 0}
    e_err = (student_e[mask] - teacher_e[mask])
    f_err = (student_f[mask] - teacher_f[mask]).reshape(-1)
    return {
        "n": int(mask.sum()),
        "e_rmse": float(np.sqrt(np.mean(e_err**2))),
        "e_mae": float(np.mean(np.abs(e_err))),
        "e_max": float(np.max(np.abs(e_err))),
        "e_bias": float(np.mean(e_err)),
        "f_rmse": float(np.sqrt(np.mean(f_err**2))),
    }


def _make_config(args, d0, width):
    """Create RunConfig with the specified damping parameters."""
    config = RunConfig()
    config.density_backend.chgcar_path = args.chgcar
    config.density_backend.locpot_path = args.locpot
    config.teacher_backend.kind = "madsurf"
    config.teacher_backend.model_path = args.model_path
    config.teacher_backend.device = args.device
    config.grid.builder_mode = "metal_density_plq"
    config.grid.interpolation_order = 3
    config.grid.metal_density_pauli_power = 1.0
    config.grid.metal_density_london_power = 0.5
    if d0 > 0.0:
        config.grid.london_damping_d0 = float(d0)
        config.grid.london_damping_width = float(width)
    config.toggles.use_ct_qeq = False
    config.toggles.use_image_charge = False
    config.toggles.use_reactive_grid = False
    config.hybrid_model.use_qeq = False
    config.hybrid_model.use_image = False
    config.hybrid_model.use_reactive_grid = False
    config.hybrid_model.use_req_plq = True
    config.training.fit_req_radius_offset = True
    config.training.fit_req_energy_scale = True
    config.training.fit_chi = False
    config.training.fit_hardness = False
    config.training.fit_image_plane = False
    config.training.fit_reactive = False
    config.training.fit_static_charge = False
    config.training.fit_sample_shift_z = False
    config.training.fit_coulomb_shift_z = False
    config.training.req_radius_regularization = 5e-2
    config.training.req_energy_regularization = 5e-3
    config.training.constraint_bound_fraction = 5e-2
    config.training.max_steps = int(args.max_steps)
    config.training.learning_rate = float(args.learning_rate)
    config.training.force_weight = float(args.force_weight)
    config.training.validation_split = 0.0
    config.training.test_split = 0.0
    return config


def _run_single_damping(args, density_loader, adsorbate, site_uvs, z_values,
                        poses, teacher, train_indices, d0, width, out_dir):
    """Build grid with damping, fit, evaluate, return metrics."""
    config = _make_config(args, d0, width)
    label = f"d0={d0:.1f}_w={width:.1f}"
    print(f"\n{'='*60}")
    print(f"  Damping: d0={d0:.2f} A, width={width:.2f} A")
    print(f"{'='*60}", flush=True)

    run_dir = ensure_dir(out_dir / f"d0_{d0:.1f}_w_{width:.1f}".replace(".", "p"))
    save_json({"d0": d0, "width": width}, run_dir / "damping_params.json")

    # Reload density with new damping config — grid builder uses config.grid
    t0 = time.perf_counter()
    density = density_loader(config)
    grid_time = time.perf_counter() - t0
    print(f"  Grid built in {grid_time:.1f}s", flush=True)

    reactive_elements = tuple(sorted(set(adsorbate.symbols)))
    model = HybridGridFFModel(
        density=density,
        reactive_elements=reactive_elements,
        toggles=config.toggles,
        grid_config=config.grid,
        hybrid_config=config.hybrid_model,
        prefer_jax=True,
    )

    train_poses = _slice_pose_batch(poses, train_indices)
    train_teacher = _slice_teacher(teacher, train_indices)

    t0 = time.perf_counter()
    fit_result = fit_hybrid_parameters(
        density=density,
        datasets=[(train_poses, train_teacher)],
        model=model,
        training=config.training,
        initial_params=None,
    )
    fit_time = time.perf_counter() - t0
    print(f"  Fit done in {fit_time:.1f}s", flush=True)

    # Evaluate on ALL poses (dense grid)
    t0 = time.perf_counter()
    student = model.evaluate_batch(poses, params=fit_result.params, compute_forces=True)
    eval_time = time.perf_counter() - t0

    teacher_e = np.asarray(teacher.energies)
    student_e = np.asarray(student.energies)
    teacher_f = np.asarray(teacher.forces)
    student_f = np.asarray(student.forces)
    nz = len(z_values)

    # Primary window mask (z >= z_min)
    z_vec = poses.pose_params[:, 2]
    primary_mask = (z_vec >= args.z_min) & (z_vec <= args.z_max)

    overall = _metrics(teacher_e, student_e, teacher_f, student_f, primary_mask)

    # Per-site metrics
    site_metrics = {}
    site_names = list(site_uvs.keys())
    for si, sname in enumerate(site_names):
        lo, hi = si * nz, (si + 1) * nz
        sm = np.zeros(len(z_vec), dtype=bool)
        sm[lo:hi] = True
        sm &= primary_mask
        site_metrics[sname] = _metrics(teacher_e, student_e, teacher_f, student_f, sm)

    # Z-region breakdown
    z_regions = {}
    for zlo, zhi, rname in [(2.0, 2.5, "2.0-2.5"), (2.5, 3.5, "2.5-3.5"),
                             (3.5, 5.0, "3.5-5.0"), (5.0, 8.0, "5.0-8.0")]:
        rm = (z_vec >= zlo) & (z_vec < zhi)
        z_regions[rname] = _metrics(teacher_e, student_e, teacher_f, student_f, rm)

    # Fitted parameters
    params = fit_result.params
    fit_params = {
        "req_radius_offset": dict(params.req_radius_offset),
        "req_energy_scale": dict(params.req_energy_scale),
    }

    # Save per-site traces
    traces_dir = ensure_dir(run_dir / "traces")
    for si, sname in enumerate(site_names):
        lo, hi = si * nz, (si + 1) * nz
        np.savez(
            traces_dir / f"{sname}_zscan_trace.npz",
            z_values=z_values,
            teacher_energy=teacher_e[lo:hi],
            student_energy=student_e[lo:hi],
            teacher_force=teacher_f[lo:hi],
            student_force=student_f[lo:hi],
            pauli=np.asarray(student.components["pauli"][lo:hi]),
            london=np.asarray(student.components["london"][lo:hi]),
            coulomb=np.asarray(student.components["coulomb"][lo:hi]),
        )

    result = {
        "d0": d0,
        "width": width,
        "overall": overall,
        "sites": site_metrics,
        "z_regions": z_regions,
        "fit_params": fit_params,
        "timing": {"grid_s": grid_time, "fit_s": fit_time, "eval_s": eval_time},
    }
    save_json(result, run_dir / "result.json")

    e_rmse_mev = overall["e_rmse"] * 1000
    print(f"  E_RMSE = {e_rmse_mev:.1f} meV  (bias={overall['e_bias']*1000:.1f} meV)")
    for sname, sm in site_metrics.items():
        print(f"    {sname:8s}: E_RMSE={sm['e_rmse']*1000:.1f} meV, bias={sm['e_bias']*1000:.1f} meV")
    for rname, rm in z_regions.items():
        print(f"    z={rname}: E_RMSE={rm['e_rmse']*1000:.1f} meV, bias={rm['e_bias']*1000:.1f} meV")
    print(f"  req_radius_offset: {fit_params['req_radius_offset']}")
    print(f"  req_energy_scale:  {fit_params['req_energy_scale']}")
    sys.stdout.flush()
    return result


def _plot_sweep_summary(results, out_dir):
    """Plot RMSE vs d0 for each width, plus z-region breakdown."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # Group by width
    by_width = {}
    for r in results:
        w = r["width"]
        by_width.setdefault(w, []).append(r)

    # Plot 1: Overall E_RMSE vs d0
    ax = axes[0]
    for w, runs in sorted(by_width.items()):
        runs_sorted = sorted(runs, key=lambda x: x["d0"])
        d0s = [r["d0"] for r in runs_sorted]
        rmses = [r["overall"]["e_rmse"] * 1000 for r in runs_sorted]
        ax.plot(d0s, rmses, "o-", label=f"w={w:.1f} A", markersize=6)
    ax.set_xlabel("London damping d0 [A]")
    ax.set_ylabel("E_RMSE [meV]")
    ax.set_title("Overall Energy RMSE vs Damping Distance")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 2: Per-site RMSE for the best width
    ax = axes[1]
    best_width = min(by_width.keys(), key=lambda w: min(r["overall"]["e_rmse"] for r in by_width[w]))
    runs_sorted = sorted(by_width[best_width], key=lambda x: x["d0"])
    d0s = [r["d0"] for r in runs_sorted]
    for site in runs_sorted[0]["sites"]:
        rmses = [r["sites"][site]["e_rmse"] * 1000 for r in runs_sorted]
        ax.plot(d0s, rmses, "o-", label=site, markersize=5)
    ax.set_xlabel("London damping d0 [A]")
    ax.set_ylabel("E_RMSE [meV]")
    ax.set_title(f"Per-Site RMSE (width={best_width:.1f} A)")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 3: Z-region breakdown for the best width
    ax = axes[2]
    for region in ["2.0-2.5", "2.5-3.5", "3.5-5.0", "5.0-8.0"]:
        rmses = [r["z_regions"][region]["e_rmse"] * 1000 for r in runs_sorted]
        ax.plot(d0s, rmses, "o-", label=f"z={region} A", markersize=5)
    ax.set_xlabel("London damping d0 [A]")
    ax.set_ylabel("E_RMSE [meV]")
    ax.set_title(f"Z-Region RMSE (width={best_width:.1f} A)")
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(out_dir / "damping_sweep_summary.png", dpi=170)
    plt.close()
    print(f"\nSummary plot saved: {out_dir / 'damping_sweep_summary.png'}")

    # Plot 4: Best z-scan curves comparison (baseline vs best damped)
    fig, axes2 = plt.subplots(1, 3, figsize=(16, 5))
    # Find baseline (d0=0) and best
    baseline = next((r for r in results if r["d0"] == 0.0), None)
    best_all = min(results, key=lambda r: r["overall"]["e_rmse"])
    if baseline and best_all["d0"] != 0.0:
        sites = list(runs_sorted[0]["sites"].keys())
        for si, site in enumerate(sites):
            ax = axes2[si]
            # Load traces
            base_dir = out_dir / "d0_0p0_w_0p5" / "traces"
            best_dir = out_dir / f"d0_{best_all['d0']:.1f}_w_{best_all['width']:.1f}".replace(".", "p") / "traces"
            try:
                bt = np.load(base_dir / f"{site}_zscan_trace.npz")
                dt = np.load(best_dir / f"{site}_zscan_trace.npz")
                z = bt["z_values"]
                ax.plot(z, bt["teacher_energy"] * 1000, "k-", lw=2, label="MACE teacher")
                ax.plot(z, bt["student_energy"] * 1000, "r--", lw=1.5, label=f"No damp (d0=0)")
                ax.plot(z, dt["student_energy"] * 1000, "b--", lw=1.5, label=f"Best (d0={best_all['d0']:.1f})")
                ax.set_xlabel("z [A]")
                ax.set_ylabel("E_int [meV]")
                ax.set_title(f"{site}")
                ax.legend(fontsize=8)
                ax.grid(True, alpha=0.3)
            except Exception as e:
                ax.text(0.5, 0.5, f"No traces:\n{e}", transform=ax.transAxes, ha="center")
        plt.suptitle(f"Z-scan: baseline vs best damping (d0={best_all['d0']:.1f})")
        plt.tight_layout()
        plt.savefig(out_dir / "damping_zscan_comparison.png", dpi=170)
        plt.close()
        print(f"Comparison plot saved: {out_dir / 'damping_zscan_comparison.png'}")


def main():
    args = parse_args()
    out_dir = ensure_dir(args.out_dir)

    d0_values = [float(x) for x in args.d0_values.split(",")]
    width_values = [float(x) for x in args.width_values.split(",")]
    n_combos = len(d0_values) * len(width_values)

    print(f"London Damping Sweep")
    print(f"  d0 values:    {d0_values}")
    print(f"  width values: {width_values}")
    print(f"  Total combos: {n_combos}")
    print(f"  Output:       {out_dir}")
    print(flush=True)

    # Step 1: Load density (without damping — just for pose building)
    print("\n[1/4] Loading density backend for pose construction...", flush=True)
    base_config = _make_config(args, d0=0.0, width=0.5)

    def density_loader(config):
        """Reload density with specific grid config (damping params change the grid)."""
        return make_density_backend(config.density_backend, surface=config.surface, grid=config.grid).load()

    density_base = density_loader(base_config)

    # Step 2: Build poses (same for all damping values)
    print("[2/4] Building scan poses...", flush=True)
    adsorbate = get_adsorbate(name="CHONH2")
    site_uvs = _site_uv_lookup(density_base, args.sites)
    z_values = np.linspace(args.stress_z_min, args.z_max, args.eval_z_points)
    quaternion = _rotation_quaternion([0, 0, 1], 0.0)
    poses = _build_zscan_poses(density_base, adsorbate, site_uvs, z_values, quaternion)
    print(f"  {len(poses.positions)} poses ({len(site_uvs)} sites x {len(z_values)} z-points)")

    train_idx = _train_indices(
        len(site_uvs), len(z_values), z_values,
        args.z_min, args.z_dense_cutoff, args.z_max,
        args.train_low_count, args.train_high_count,
    )
    print(f"  {len(train_idx)} training points", flush=True)

    # Step 3: Run MACE teacher ONCE
    print("\n[3/4] Running MACE teacher evaluation (ONE TIME)...", flush=True)
    teacher_backend = make_teacher_backend(base_config.teacher_backend)
    teacher, teacher_time = _evaluate_teacher_chunked(teacher_backend, density_base, poses)
    print(f"  MACE done in {teacher_time:.1f}s ({teacher_time/len(poses.positions):.3f}s/pose)")

    # Cache teacher results
    np.savez(
        out_dir / "teacher_cache.npz",
        energies=np.asarray(teacher.energies),
        forces=np.asarray(teacher.forces),
    )
    print(f"  Teacher cache saved: {out_dir / 'teacher_cache.npz'}", flush=True)

    # Free MACE model memory
    del teacher_backend
    import gc; gc.collect()

    # Step 4: Sweep damping parameters
    print(f"\n[4/4] Sweeping {n_combos} damping combinations...", flush=True)
    results = []
    for d0 in d0_values:
        for width in width_values:
            result = _run_single_damping(
                args, density_loader, adsorbate, site_uvs, z_values,
                poses, teacher, train_idx, d0, width, out_dir,
            )
            results.append(result)

    # Summary table
    print(f"\n{'='*80}")
    print(f"  DAMPING SWEEP SUMMARY")
    print(f"{'='*80}")
    print(f"{'d0':>6s} {'width':>6s} | {'E_RMSE':>10s} {'E_bias':>10s} {'F_RMSE':>10s} | {'z<2.5':>8s} {'2.5-3.5':>8s} {'3.5-5':>8s} {'5-8':>8s}")
    print("-" * 90)
    for r in sorted(results, key=lambda x: (x["width"], x["d0"])):
        z25 = r["z_regions"]["2.0-2.5"]["e_rmse"] * 1000
        z35 = r["z_regions"]["2.5-3.5"]["e_rmse"] * 1000
        z50 = r["z_regions"]["3.5-5.0"]["e_rmse"] * 1000
        z80 = r["z_regions"]["5.0-8.0"]["e_rmse"] * 1000
        print(f"{r['d0']:6.1f} {r['width']:6.1f} | "
              f"{r['overall']['e_rmse']*1000:10.1f} {r['overall']['e_bias']*1000:10.1f} {r['overall']['f_rmse']*1000:10.1f} | "
              f"{z25:8.1f} {z35:8.1f} {z50:8.1f} {z80:8.1f}")

    best = min(results, key=lambda r: r["overall"]["e_rmse"])
    print(f"\nBest: d0={best['d0']:.1f}, width={best['width']:.1f} → E_RMSE = {best['overall']['e_rmse']*1000:.1f} meV")
    print(f"  Fit params: {best['fit_params']}")

    # Save summary
    save_json(
        {"results": results, "best": {"d0": best["d0"], "width": best["width"],
         "e_rmse_mev": best["overall"]["e_rmse"] * 1000}},
        out_dir / "sweep_summary.json",
    )

    # Generate plots
    _plot_sweep_summary(results, out_dir)

    print(f"\nAll results saved in: {out_dir}")


if __name__ == "__main__":
    main()
