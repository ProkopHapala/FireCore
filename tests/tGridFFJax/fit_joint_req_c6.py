#!/usr/bin/env python3
"""Joint REQ + C₆ optimization for metal GridFF.

Compares three configurations:
  A) PLQ-only (REQ Morse, no C₆) — baseline
  B) PLQ + C₆ grid (joint REQ + C₆ optimization)
  C) PLQ + C₆ grid (C₆ only, frozen REQ from baseline A)

Uses the same teacher data and fitting pipeline as sweep_london_damping.py.
"""

from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("MKL_NUM_THREADS", "4")
os.environ.setdefault("JAX_PLATFORMS", "cpu")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.append(str(ROOT))

from pyBall.gridff_jax import RunConfig
from pyBall.gridff_jax.density_backends import make_density_backend
from pyBall.gridff_jax.fit import fit_hybrid_parameters
from pyBall.gridff_jax.hybrid_energy import HybridGridFFModel, HybridParameters
from pyBall.gridff_jax.interfaces import PoseBatch, TeacherResult
from pyBall.gridff_jax.pose_sampling import (
    get_adsorbate, infer_reference_planes, infer_surface_sites, transform_adsorbate,
)
from pyBall.gridff_jax.utils import ensure_dir, save_json

TEACHER_CACHE = Path(__file__).resolve().parent / "damping_sweep" / "teacher_cache.npz"
OUT_DIR = ensure_dir(Path(__file__).resolve().parent / "joint_req_c6_fit")

DEFAULT_CHGCAR = "/home/niel/git/ORR_HER_Ag_Colab/results/Ag_ORR_HER/slab_clean/final_scf_12x12x1/CHGCAR"
DEFAULT_LOCPOT = "/home/niel/git/ORR_HER_Ag_Colab/results/Ag_ORR_HER/slab_clean/workfunc_12x12x1/LOCPOT"


def _make_config(use_c6=False, fit_c6=False, d0=3.70, width=0.35, max_steps=600):
    config = RunConfig()
    config.density_backend.chgcar_path = DEFAULT_CHGCAR
    config.density_backend.locpot_path = DEFAULT_LOCPOT
    config.grid.builder_mode = "metal_density_plq"
    config.grid.interpolation_order = 3
    config.grid.metal_density_pauli_power = 1.0
    config.grid.metal_density_london_power = 0.5
    config.grid.london_damping_d0 = float(d0)
    config.grid.london_damping_width = float(width)
    config.grid.use_pairwise_c6 = use_c6
    config.toggles.use_ct_qeq = False
    config.toggles.use_image_charge = False
    config.toggles.use_reactive_grid = False
    config.hybrid_model.use_qeq = False
    config.hybrid_model.use_image = False
    config.hybrid_model.use_reactive_grid = False
    config.hybrid_model.use_req_plq = True
    config.training.fit_req_radius_offset = True
    config.training.fit_req_energy_scale = True
    config.training.fit_c6_coeff = fit_c6
    config.training.fit_chi = False
    config.training.fit_hardness = False
    config.training.fit_image_plane = False
    config.training.fit_reactive = False
    config.training.fit_static_charge = False
    config.training.fit_sample_shift_z = False
    config.training.fit_coulomb_shift_z = False
    config.training.req_radius_regularization = 5e-2
    config.training.req_energy_regularization = 5e-3
    config.training.max_steps = max_steps
    config.training.learning_rate = 0.01
    config.training.force_weight = 10.0
    config.training.validation_split = 0.0
    config.training.test_split = 0.0
    return config


def _build_poses(density, adsorbate, site_uvs, z_values):
    planes = infer_reference_planes(density)
    z_ref = float(planes["z_ref"])
    z_image = float(planes["z_image"])
    quat = np.array([1, 0, 0, 0], dtype=float)
    positions, pose_params, labels = [], [], []
    for label, uv in site_uvs.items():
        for zh in z_values:
            positions.append(transform_adsorbate(adsorbate, density, uv, float(zh), quat, z_ref))
            pose_params.append(np.concatenate([uv, [float(zh)], quat]))
            labels.append(label)
    return PoseBatch(
        adsorbate=adsorbate,
        positions=np.asarray(positions),
        pose_params=np.asarray(pose_params),
        site_labels=labels,
        metadata={"z_ref": z_ref, "z_image": z_image},
    )


def _train_indices(n_sites, nz, z_values, z_min=2.0, z_dense_cutoff=3.5, z_max=8.0, n_low=12, n_high=8):
    """Select training subset: denser at short range, sparser at long range."""
    z = np.asarray(z_values)
    low = np.where((z >= z_min) & (z < z_dense_cutoff))[0]
    high = np.where((z >= z_dense_cutoff) & (z <= z_max))[0]
    local = np.unique(np.concatenate([
        np.linspace(0, len(low) - 1, min(n_low, len(low))).astype(int) if len(low) > 0 else np.array([], dtype=int),
        np.linspace(0, len(high) - 1, min(n_high, len(high))).astype(int) + len(low) if len(high) > 0 else np.array([], dtype=int),
    ]).astype(int))
    # Map back to full z indices
    all_valid = np.where((z >= z_min) & (z <= z_max))[0]
    selected = all_valid[np.clip(local, 0, len(all_valid) - 1)]
    indices = []
    for i in range(n_sites):
        indices.extend((i * nz + selected).tolist())
    return np.asarray(indices, dtype=int)


def _slice(poses, teacher, indices):
    idx = np.asarray(indices, dtype=int)
    return (
        PoseBatch(
            adsorbate=poses.adsorbate,
            positions=poses.positions[idx],
            pose_params=poses.pose_params[idx],
            site_labels=[poses.site_labels[int(i)] for i in idx],
            metadata=dict(poses.metadata),
        ),
        TeacherResult(
            energies=teacher.energies[idx],
            forces=teacher.forces[idx],
            metadata=dict(teacher.metadata),
        ),
    )


def _per_site_metrics(teacher_e, student_e, labels, z_params, nz, z_values, sites):
    """Compute per-site RMSE and per-z-region errors."""
    results = {}
    for site in sites:
        mask = np.array([l == site for l in labels])
        if mask.sum() == 0:
            continue
        err = (student_e[mask] - teacher_e[mask]) * 1000  # meV
        z_site = np.array([z_params[k] for k in range(len(labels)) if labels[k] == site])
        results[site] = {
            "rmse_meV": float(np.sqrt(np.mean(err ** 2))),
            "mae_meV": float(np.mean(np.abs(err))),
            "max_meV": float(np.max(np.abs(err))),
            "bias_meV": float(np.mean(err)),
            "well_depth_mace": float(np.min(teacher_e[mask]) * 1000),
            "well_depth_gridff": float(np.min(student_e[mask]) * 1000),
        }
        # Short range (z < 3.0) vs long range
        z_short = z_site < 3.0
        z_long = z_site >= 3.0
        if z_short.sum() > 0:
            results[site]["rmse_short_meV"] = float(np.sqrt(np.mean(err[z_short] ** 2)))
        if z_long.sum() > 0:
            results[site]["rmse_long_meV"] = float(np.sqrt(np.mean(err[z_long] ** 2)))
    return results


def _run_fit(label, config, density, poses, teacher, train_indices, max_steps, initial_params=None):
    """Run a single fitting configuration and evaluate on ALL poses."""
    print(f"\n{'='*60}")
    print(f"  {label}")
    print(f"{'='*60}", flush=True)

    reactive_elements = tuple(sorted(set(poses.adsorbate.symbols)))
    t0 = time.perf_counter()
    model = HybridGridFFModel(
        density=density,
        reactive_elements=reactive_elements,
        toggles=config.toggles,
        grid_config=config.grid,
        hybrid_config=config.hybrid_model,
        prefer_jax=True,
    )
    model_time = time.perf_counter() - t0
    print(f"  Model built in {model_time:.1f}s", flush=True)

    # Check dispersion grid
    if config.grid.use_pairwise_c6:
        disp = model.grids["dispersion_zyxc"]
        print(f"  Dispersion grid: shape={disp.shape}, min={disp.min():.2f}, max={disp.max():.6f}")

    train_poses, train_teacher = _slice(poses, teacher, train_indices)
    print(f"  Training on {len(train_indices)} / {len(poses.positions)} poses")

    config.training.max_steps = max_steps
    t0 = time.perf_counter()
    fit_result = fit_hybrid_parameters(
        density=density,
        datasets=[(train_poses, train_teacher)],
        model=model,
        training=config.training,
        initial_params=initial_params,
    )
    fit_time = time.perf_counter() - t0
    nit = fit_result.optimizer_result.get("nit", "?") if fit_result.optimizer_result else "?"
    print(f"  Fit done in {fit_time:.1f}s ({nit} iterations)", flush=True)

    # Evaluate on ALL poses
    student = model.evaluate_batch(poses, params=fit_result.params, compute_forces=True)
    teacher_e = np.asarray(teacher.energies)
    student_e = np.asarray(student.energies)
    teacher_f = np.asarray(teacher.forces)
    student_f = np.asarray(student.forces)

    # Overall metrics
    rmse_e = float(np.sqrt(np.mean((student_e - teacher_e) ** 2))) * 1000
    mae_e = float(np.mean(np.abs(student_e - teacher_e))) * 1000
    rmse_f = float(np.sqrt(np.mean((student_f - teacher_f) ** 2)))

    print(f"\n  Results ({label}):")
    print(f"    Energy RMSE: {rmse_e:.1f} meV")
    print(f"    Energy MAE:  {mae_e:.1f} meV")
    print(f"    Force RMSE:  {rmse_f:.4f} eV/A")

    # Print fitted parameters
    print(f"\n  Fitted parameters:")
    for el in reactive_elements:
        ro = fit_result.params.req_radius_offset.get(el, 0.0)
        es = fit_result.params.req_energy_scale.get(el, 1.0)
        c6 = fit_result.params.c6_coeff.get(el, 1.0)
        print(f"    {el}: req_radius_offset={ro:.4f}, req_energy_scale={es:.4f}, c6_coeff={c6:.4f}")

    return {
        "label": label,
        "model": model,
        "fit_result": fit_result,
        "student": student,
        "rmse_e_meV": rmse_e,
        "mae_e_meV": mae_e,
        "rmse_f": rmse_f,
        "teacher_e": teacher_e,
        "student_e": student_e,
        "teacher_f": teacher_f,
        "student_f": student_f,
        "fit_time": fit_time,
        "model_time": model_time,
    }


def main():
    print("=" * 70)
    print("  JOINT REQ + C6 OPTIMIZATION")
    print("  CHONH2 / Ag(111) — vdW-surf screened parameters")
    print("=" * 70)
    print(f"Output: {OUT_DIR}")
    t_start = time.perf_counter()

    # Load teacher data
    if not TEACHER_CACHE.exists():
        print(f"ERROR: Teacher cache not found: {TEACHER_CACHE}")
        print("Run sweep_london_damping.py first to generate teacher data.")
        return
    tc = np.load(TEACHER_CACHE)
    teacher = TeacherResult(energies=tc["energies"], forces=tc["forces"], metadata={})
    print(f"Teacher: {len(teacher.energies)} poses loaded from cache")

    # Build poses (same as sweep script)
    config_base = _make_config(use_c6=False)
    print("Loading density for pose construction...", flush=True)
    density_base = make_density_backend(
        config_base.density_backend, surface=config_base.surface, grid=config_base.grid
    ).load()

    adsorbate = get_adsorbate(name="CHONH2")
    site_uvs = {}
    for site in infer_surface_sites(density_base):
        if site.label in ["top", "bridge", "hollow"]:
            site_uvs[site.label] = np.asarray(site.uv)
    sites = list(site_uvs.keys())

    z_values = np.linspace(1.8, 8.0, 63)
    nz = len(z_values)
    poses = _build_poses(density_base, adsorbate, site_uvs, z_values)

    # Training subset
    train_idx = _train_indices(len(sites), nz, z_values)
    print(f"Training indices: {len(train_idx)} / {len(poses.positions)}")

    # ============================================================
    # CONFIG A: PLQ-only (baseline)
    # ============================================================
    config_a = _make_config(use_c6=False, fit_c6=False, max_steps=600)
    result_a = _run_fit(
        "A: PLQ-only (baseline)",
        config_a, density_base, poses, teacher, train_idx, max_steps=600,
    )

    # ============================================================
    # CONFIG B: Joint REQ + C₆ (full optimization)
    # ============================================================
    config_b = _make_config(use_c6=True, fit_c6=True, max_steps=800)
    print("\nBuilding density + C₆ grid for Config B...", flush=True)
    density_c6 = make_density_backend(
        config_b.density_backend, surface=config_b.surface, grid=config_b.grid
    ).load()
    result_b = _run_fit(
        "B: Joint REQ + C6 optimization",
        config_b, density_c6, poses, teacher, train_idx, max_steps=800,
    )

    # ============================================================
    # CONFIG C: C₆ only (frozen REQ from baseline A)
    # ============================================================
    config_c = _make_config(use_c6=True, fit_c6=True, max_steps=800)
    config_c.training.fit_req_radius_offset = False
    config_c.training.fit_req_energy_scale = False
    result_c = _run_fit(
        "C: C6-only (frozen REQ from A)",
        config_c, density_c6, poses, teacher, train_idx, max_steps=800,
        initial_params=result_a["fit_result"].params,
    )

    # ============================================================
    # PER-SITE ANALYSIS
    # ============================================================
    print("\n" + "=" * 70)
    print("  PER-SITE ANALYSIS")
    print("=" * 70)

    z_params = [poses.pose_params[k, 2] for k in range(len(poses.site_labels))]
    for result in [result_a, result_b, result_c]:
        print(f"\n  {result['label']}:")
        site_metrics = _per_site_metrics(
            result["teacher_e"], result["student_e"],
            poses.site_labels, z_params, nz, z_values, sites,
        )
        for site, m in site_metrics.items():
            print(f"    {site:8s}: RMSE={m['rmse_meV']:.1f} meV, "
                  f"MAE={m['mae_meV']:.1f}, max={m['max_meV']:.1f}, "
                  f"bias={m['bias_meV']:+.1f}")
            if "rmse_short_meV" in m:
                print(f"             short(z<3)={m['rmse_short_meV']:.1f} meV, "
                      f"long(z>=3)={m.get('rmse_long_meV', 0):.1f} meV")
            print(f"             well: MACE={m['well_depth_mace']:.1f} meV, "
                  f"GridFF={m['well_depth_gridff']:.1f} meV")

    # ============================================================
    # FORCE ACCURACY
    # ============================================================
    print("\n" + "=" * 70)
    print("  FORCE ACCURACY (z-component on molecular atoms)")
    print("=" * 70)
    for result in [result_a, result_b, result_c]:
        fz_err = (result["student_f"][:, :, 2] - result["teacher_f"][:, :, 2]).ravel()
        print(f"  {result['label']}: Fz RMSE = {np.sqrt(np.mean(fz_err**2)):.4f} eV/A")

    # ============================================================
    # PLOTS
    # ============================================================
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))

    colors = {"A: PLQ-only (baseline)": "blue",
              "B: Joint REQ + C6 optimization": "red",
              "C: C6-only (frozen REQ from A)": "green"}

    for idx, site in enumerate(sites):
        ax_e = axes[0, idx]
        ax_r = axes[1, idx]
        mask = np.array([l == site for l in poses.site_labels])
        z_site = np.array([z_params[k] for k in range(len(poses.site_labels)) if poses.site_labels[k] == site])
        e_mace = teacher.energies[mask] * 1000

        ax_e.plot(z_site, e_mace, 'ko-', ms=3, lw=1, label='MACE', zorder=10)
        ax_r.axhline(0, color='grey', lw=0.5)

        for result in [result_a, result_b, result_c]:
            e_student = result["student_e"][mask] * 1000
            short_label = result["label"].split(":")[0] + f" ({result['rmse_e_meV']:.0f}meV)"
            color = list(colors.values())[list(colors.keys()).index(result["label"])]
            ax_e.plot(z_site, e_student, '-', lw=1.5, color=color, label=short_label)
            ax_r.plot(z_site, e_student - e_mace, '-', lw=1.5, color=color, label=short_label)

        ax_e.set_title(f'{site} site', fontsize=12)
        ax_e.set_ylabel('Energy (meV)')
        ax_e.legend(fontsize=7)
        ax_e.grid(True, alpha=0.3)

        ax_r.set_xlabel('z (A above surface)')
        ax_r.set_ylabel('Error (meV)')
        ax_r.legend(fontsize=7)
        ax_r.grid(True, alpha=0.3)

    plt.suptitle('Joint REQ + C₆ Optimization: CHONH2 / Ag(111)\n'
                 f'A={result_a["rmse_e_meV"]:.1f} meV → B={result_b["rmse_e_meV"]:.1f} meV '
                 f'({(1 - result_b["rmse_e_meV"]/result_a["rmse_e_meV"])*100:.0f}% improvement)',
                 fontsize=13)
    plt.tight_layout()
    plt.savefig(OUT_DIR / "joint_req_c6_comparison.png", dpi=150)
    print(f"\nPlot saved: {OUT_DIR / 'joint_req_c6_comparison.png'}")

    # ============================================================
    # SAVE RESULTS
    # ============================================================
    summary = {
        "configs": {},
        "best_config": None,
        "best_rmse_meV": 999.0,
    }
    for result in [result_a, result_b, result_c]:
        key = result["label"]
        params = result["fit_result"].params
        summary["configs"][key] = {
            "rmse_e_meV": result["rmse_e_meV"],
            "mae_e_meV": result["mae_e_meV"],
            "rmse_f": result["rmse_f"],
            "fit_time_s": result["fit_time"],
            "model_time_s": result["model_time"],
            "params": {
                "req_radius_offset": dict(params.req_radius_offset),
                "req_energy_scale": dict(params.req_energy_scale),
                "c6_coeff": dict(params.c6_coeff),
                "sample_shift_z": params.sample_shift_z,
                "coulomb_shift_z": params.coulomb_shift_z,
            },
        }
        if result["rmse_e_meV"] < summary["best_rmse_meV"]:
            summary["best_rmse_meV"] = result["rmse_e_meV"]
            summary["best_config"] = key

    save_json(summary, OUT_DIR / "joint_fit_summary.json")

    # Save best params for downstream use
    best = [r for r in [result_a, result_b, result_c] if r["label"] == summary["best_config"]][0]
    best_params = best["fit_result"].params
    save_json({
        "pauli": dict(best_params.pauli),
        "london": dict(best_params.london),
        "reactive": dict(best_params.reactive),
        "static_charge": dict(best_params.static_charge),
        "c6_coeff": dict(best_params.c6_coeff),
        "req_radius_offset": dict(best_params.req_radius_offset),
        "req_energy_scale": dict(best_params.req_energy_scale),
        "chi": dict(best_params.chi),
        "hardness": dict(best_params.hardness),
        "image_scale": best_params.image_scale,
        "image_plane": best_params.image_plane,
        "image_damping": best_params.image_damping,
        "sample_shift_z": best_params.sample_shift_z,
        "coulomb_shift_z": best_params.coulomb_shift_z,
    }, OUT_DIR / "best_fit_params.json")

    # ============================================================
    # FINAL SUMMARY
    # ============================================================
    elapsed = time.perf_counter() - t_start
    print("\n" + "=" * 70)
    print("  FINAL RESULTS")
    print("=" * 70)
    for result in [result_a, result_b, result_c]:
        marker = " *** BEST ***" if result["label"] == summary["best_config"] else ""
        print(f"  {result['label']}:")
        print(f"    Energy RMSE: {result['rmse_e_meV']:.1f} meV  |  Force RMSE: {result['rmse_f']:.4f} eV/A{marker}")
    if result_b["rmse_e_meV"] < result_a["rmse_e_meV"]:
        improvement = (1 - result_b["rmse_e_meV"] / result_a["rmse_e_meV"]) * 100
        print(f"\n  Joint optimization improvement: {improvement:.1f}%")
    print(f"\n  Total elapsed: {elapsed:.0f}s")
    print(f"  Results saved: {OUT_DIR}")
    print("=" * 70)


if __name__ == "__main__":
    main()
