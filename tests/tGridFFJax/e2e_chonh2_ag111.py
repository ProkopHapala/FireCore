#!/usr/bin/env python3
"""End-to-end CHONH2/Ag(111) GridFF benchmark.

Three modes:
    --mode generate   Generate ultra-fine 6D reference dataset (~2M poses)
    --mode fit        Load dataset, fit GridFF params (energy + force)
    --mode scan       Run 4 validation scans (z, x, y, xy)
    --mode all        Run all phases sequentially

Best two-stage result (target ≤ 16 meV RMSE):
    python e2e_chonh2_ag111.py --mode all --two-stage-c6 --raw-r6 --prefer-jax
"""

from __future__ import annotations

import argparse
import json
import os
import time
from copy import deepcopy
from pathlib import Path
import sys

# Force JAX to CPU to avoid CuDNN version mismatch on GPU
os.environ.setdefault("JAX_PLATFORMS", "cpu")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from pyBall.gridff_jax import RunConfig, backend_summary, save_config
from pyBall.gridff_jax.density_backends import make_density_backend
from pyBall.gridff_jax.export import export_firecore_artifacts
from pyBall.gridff_jax.fit import fit_hybrid_parameters, fit_two_stage_c6
from pyBall.gridff_jax.fit.fit_6d import (
    Dataset6D,
    SplitDataset6D,
    generate_6d_reference_data,
    load_6d_dataset,
    save_6d_dataset,
    split_6d_dataset,
)
from pyBall.gridff_jax.hybrid_energy import HybridGridFFModel, HybridParameters
from pyBall.gridff_jax.interfaces import ModelEvaluation, PoseBatch, TeacherResult
from pyBall.gridff_jax.pose_sampling import (
    get_adsorbate,
    infer_reference_planes,
    infer_surface_sites,
    transform_adsorbate,
)
from pyBall.gridff_jax.pose_sampling.sampler_6d import Sampler6DConfig
from pyBall.gridff_jax.teacher_backends import make_teacher_backend
from pyBall.gridff_jax.utils import (
    ensure_dir,
    fractional_to_cartesian,
    quaternion_to_matrix,
    save_json,
)
from pyBall.gridff_jax.validation import (
    compute_error_metrics,
    plot_convergence,
    plot_error_histogram,
    plot_parity,
)

# ─── Default paths ────────────────────────────────────────────────────
_BASE = Path(__file__).resolve().parent
_DFT_ROOT = Path("/home/indranil/git/ORR_HER_Ag_Colab/results/slab_clean")
DEFAULT_CHGCAR = str(_DFT_ROOT / "final_scf_12x12x1" / "CHGCAR")
DEFAULT_LOCPOT = str(_DFT_ROOT / "workfunc_12x12x1" / "LOCPOT")
DEFAULT_MODEL = str(_BASE / "MAD-SURF-main" / "models" / "full_dataset_config_weights" / "MACE_model.model")

# ─── Plot colours ─────────────────────────────────────────────────────
CLR_TEACHER = "tab:blue"
CLR_STUDENT = "tab:orange"
CLR_ERROR = "#cc3333"


# =====================================================================
#  CLI
# =====================================================================

def parse_args():
    p = argparse.ArgumentParser(
        description="End-to-end CHONH2 / Ag(111) GridFF benchmark",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--mode", default="all", choices=["generate", "fit", "scan", "all"])
    p.add_argument("--out-dir", default=str(_BASE / "e2e_chonh2_ag111_out"))

    # Data paths
    p.add_argument("--chgcar", default=DEFAULT_CHGCAR)
    p.add_argument("--locpot", default=DEFAULT_LOCPOT)
    p.add_argument("--model-path", default=DEFAULT_MODEL)
    p.add_argument("--device", default="cuda")

    # Adsorbate
    p.add_argument("--adsorbate-name", default="CHONH2")
    p.add_argument("--xyz-path", default=None)
    p.add_argument("--anchor-index", type=int, default=0)
    p.add_argument("--use-input-charges", action="store_true")

    # 6D sampling (Phase 1)
    p.add_argument("--n-u", type=int, default=30)
    p.add_argument("--n-v", type=int, default=30)
    p.add_argument("--n-z", type=int, default=80)
    p.add_argument("--z-min", type=float, default=2.0)
    p.add_argument("--z-max", type=float, default=20.0)
    p.add_argument("--z-bias-power", type=float, default=2.0)
    p.add_argument("--n-orient", type=int, default=24)
    p.add_argument("--tilt-max-deg", type=float, default=60.0)
    p.add_argument("--n-yaw", type=int, default=6)
    p.add_argument("--random-fraction", type=float, default=0.15)
    p.add_argument("--teacher-chunk-size", type=int, default=32)

    # Fitting (Phase 2)
    p.add_argument("--max-steps", type=int, default=800)
    p.add_argument("--force-weight", type=float, default=10.0)
    p.add_argument("--learning-rate", type=float, default=1.0e-2)
    p.add_argument("--london-damp-d0", type=float, default=3.70)
    p.add_argument("--london-damp-width", type=float, default=0.35)
    p.add_argument("--two-stage-c6", action="store_true")
    p.add_argument("--raw-r6", action="store_true")
    p.add_argument("--energy-cutoff", type=float, default=2.0,
                   help="Discard poses with E > this value (eV) before fitting. "
                        "Collision poses from tilted molecules at low z can reach "
                        "tens of eV and destroy the loss landscape.")
    p.add_argument("--max-train-poses", type=int, default=1500,
                   help="Subsample training set for Stage 1 gradient descent. "
                        "L-BFGS-B with finite differences is O(N*P) per iteration. "
                        "Stage 2 lstsq always uses ALL filtered data.")
    p.add_argument("--teacher-tile", type=str, default="auto")
    p.add_argument("--prefer-jax", dest="prefer_jax", action="store_true", default=True)
    p.add_argument("--no-prefer-jax", dest="prefer_jax", action="store_false")

    # Scan (Phase 3)
    p.add_argument("--scan-z-min", type=float, default=2.0)
    p.add_argument("--scan-z-max", type=float, default=20.0)
    p.add_argument("--scan-z-points", type=int, default=300)
    p.add_argument("--scan-lateral-max", type=float, default=6.0,
                   help="Cartesian range for x/y/xy scans (Angstrom)")
    p.add_argument("--scan-lateral-points", type=int, default=100)
    p.add_argument("--scan-xy-points", type=int, default=40)
    p.add_argument("--scan-z-height", type=float, default=None,
                   help="Fixed z for lateral scans; auto-detected from z-scan minimum if None")

    p.add_argument("--student-chunk-size", type=int, default=64)

    return p.parse_args()


# =====================================================================
#  Helpers
# =====================================================================

def _setup_config(args) -> RunConfig:
    """Build a RunConfig from CLI arguments."""
    config = RunConfig()
    config.density_backend.chgcar_path = args.chgcar
    config.density_backend.locpot_path = args.locpot
    config.teacher_backend.kind = "madsurf"
    config.teacher_backend.model_path = args.model_path
    config.teacher_backend.device = args.device

    tile_str = args.teacher_tile
    if tile_str.strip().lower() == "auto":
        config.teacher_backend.teacher_tile = (0, 0)
    else:
        parts = tile_str.split(",")
        config.teacher_backend.teacher_tile = (int(parts[0]), int(parts[1]))

    config.grid.builder_mode = "metal_density_plq"
    config.grid.interpolation_order = 3
    if args.london_damp_d0 > 0.0:
        config.grid.london_damping_d0 = float(args.london_damp_d0)
        config.grid.london_damping_width = float(args.london_damp_width)
    if args.two_stage_c6:
        config.grid.use_pairwise_c6 = True
        config.grid.c6_rcut = 15.0
        config.training.two_stage_c6 = True

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
    config.training.req_radius_regularization = 5.0e-2
    config.training.req_energy_regularization = 5.0e-3
    config.training.constraint_bound_fraction = 5.0e-2
    config.training.max_steps = int(args.max_steps)
    config.training.learning_rate = float(args.learning_rate)
    config.training.force_weight = float(args.force_weight)
    config.training.validation_split = 0.0
    config.training.test_split = 0.0
    config.export.out_dir = str(args.out_dir)
    return config


def _build_model(config, density, adsorbate, prefer_jax=True, install_raw_r6=False):
    reactive_elements = tuple(sorted(set(adsorbate.symbols)))
    model = HybridGridFFModel(
        density=density,
        reactive_elements=reactive_elements,
        toggles=config.toggles,
        grid_config=config.grid,
        hybrid_config=config.hybrid_model,
        prefer_jax=prefer_jax,
    )
    if install_raw_r6:
        print("[e2e] Installing raw 1/r^6 dispersion grid ...", flush=True)
        model.install_raw_r6_grid()
    return model


def _slice_pose_batch(poses: PoseBatch, indices):
    indices = np.asarray(indices, dtype=int)
    return PoseBatch(
        adsorbate=poses.adsorbate,
        positions=np.asarray(poses.positions[indices], dtype=float),
        pose_params=np.asarray(poses.pose_params[indices], dtype=float),
        site_labels=[poses.site_labels[int(i)] for i in indices],
        metadata=dict(poses.metadata),
    )


def _evaluate_teacher_chunked(teacher_backend, density, poses, chunk_size=32):
    n = len(poses.positions)
    energies, forces = [], []
    t0 = time.perf_counter()
    for start in range(0, n, chunk_size):
        stop = min(start + chunk_size, n)
        idx = np.arange(start, stop, dtype=int)
        partial = teacher_backend.evaluate_batch(density, _slice_pose_batch(poses, idx))
        energies.append(np.asarray(partial.energies, dtype=float))
        forces.append(np.asarray(partial.forces, dtype=float))
        if (stop % (chunk_size * 20) == 0) or stop == n:
            print(f"  MAD-SURF: {stop}/{n} ({100.0 * stop / n:.0f}%)", flush=True)
    elapsed = time.perf_counter() - t0
    return TeacherResult(
        energies=np.concatenate(energies, axis=0),
        forces=np.concatenate(forces, axis=0),
        metadata={"teacher_eval_seconds": elapsed},
    ), elapsed


def _evaluate_student_chunked(model, poses, params, chunk_size=64):
    n = len(poses.positions)
    all_e, all_f, all_comp = [], [], {}
    t0 = time.perf_counter()
    for start in range(0, n, chunk_size):
        stop = min(start + chunk_size, n)
        idx = np.arange(start, stop, dtype=int)
        partial = model.evaluate_batch(
            _slice_pose_batch(poses, idx), params=params, compute_forces=True
        )
        all_e.append(np.asarray(partial.energies, dtype=float))
        all_f.append(np.asarray(partial.forces, dtype=float))
        for key, val in partial.components.items():
            all_comp.setdefault(key, []).append(np.asarray(val, dtype=float))
    elapsed = time.perf_counter() - t0
    return ModelEvaluation(
        energies=np.concatenate(all_e, axis=0),
        forces=np.concatenate(all_f, axis=0),
        components={k: np.concatenate(v, axis=0) for k, v in all_comp.items()},
        charges=None,
    ), elapsed


def _force_norms(forces_array):
    f = np.asarray(forces_array, dtype=float)
    return np.linalg.norm(f.reshape(f.shape[0], -1), axis=1)


def _metrics(teacher_e, student_e, teacher_f, student_f):
    e_err = np.asarray(student_e, dtype=float) - np.asarray(teacher_e, dtype=float)
    f_err = np.asarray(student_f, dtype=float).ravel() - np.asarray(teacher_f, dtype=float).ravel()
    return {
        "energy": compute_error_metrics(
            np.asarray(teacher_e, dtype=float),
            np.asarray(student_e, dtype=float),
        ),
        "force": compute_error_metrics(
            np.asarray(teacher_f, dtype=float).ravel(),
            np.asarray(student_f, dtype=float).ravel(),
        ),
    }


def _fit_param_payload(params: HybridParameters):
    payload = {
        "pauli": dict(params.pauli),
        "london": dict(params.london),
        "reactive": dict(params.reactive),
        "static_charge": dict(params.static_charge),
        "req_radius_offset": dict(params.req_radius_offset),
        "req_energy_scale": dict(params.req_energy_scale),
        "chi": dict(params.chi),
        "hardness": dict(params.hardness),
        "image_scale": float(params.image_scale),
        "image_plane": float(params.image_plane),
        "image_damping": float(params.image_damping),
        "sample_shift_z": float(params.sample_shift_z),
        "coulomb_shift_z": float(params.coulomb_shift_z),
        "reservoir_chi": float(params.reservoir_chi),
        "reservoir_hardness": float(params.reservoir_hardness),
        "total_charge": float(params.total_charge),
        "energy_offset": float(params.energy_offset),
    }
    if params.c6_coeff:
        payload["c6_coeff"] = dict(params.c6_coeff)
    return payload


def _load_fitted_params(path) -> HybridParameters:
    with open(path) as fh:
        raw = json.load(fh)
    return HybridParameters(
        pauli=dict(raw.get("pauli", {})),
        london=dict(raw.get("london", {})),
        reactive=dict(raw.get("reactive", {})),
        static_charge=dict(raw.get("static_charge", {})),
        req_radius_offset=dict(raw.get("req_radius_offset", {})),
        req_energy_scale=dict(raw.get("req_energy_scale", {})),
        chi=dict(raw.get("chi", {})),
        hardness=dict(raw.get("hardness", {})),
        c6_coeff=dict(raw.get("c6_coeff", {})),
        image_scale=float(raw.get("image_scale", 1.0)),
        image_plane=float(raw.get("image_plane", 0.0)),
        image_damping=float(raw.get("image_damping", 0.5)),
        sample_shift_z=float(raw.get("sample_shift_z", 0.0)),
        coulomb_shift_z=float(raw.get("coulomb_shift_z", 0.0)),
        reservoir_chi=float(raw.get("reservoir_chi", -4.5)),
        reservoir_hardness=float(raw.get("reservoir_hardness", 7.0)),
        total_charge=float(raw.get("total_charge", 0.0)),
        energy_offset=float(raw.get("energy_offset", 0.0)),
    )


def _build_scan_poses(density, adsorbate, uv_array, z_height, quaternion):
    """Build PoseBatch for a set of (u,v) positions at fixed z and orientation."""
    planes = infer_reference_planes(density)
    z_ref = float(planes["z_ref"])
    positions, pose_params = [], []
    for uv in uv_array:
        uv = np.asarray(uv, dtype=float)
        pos = transform_adsorbate(adsorbate, density, uv, float(z_height), quaternion, z_ref)
        positions.append(pos)
        pose_params.append(np.concatenate([uv, [float(z_height)], quaternion]))
    return PoseBatch(
        adsorbate=adsorbate,
        positions=np.asarray(positions, dtype=float),
        pose_params=np.asarray(pose_params, dtype=float),
        site_labels=["scan"] * len(uv_array),
        metadata={"z_ref": z_ref, "z_image": float(planes["z_image"]),
                  "z_height": float(z_height), "scan_type": "lateral"},
    )


def _build_zscan_poses(density, adsorbate, site_uvs, z_values, quaternion):
    """Build PoseBatch for z-scan at multiple high-symmetry sites."""
    planes = infer_reference_planes(density)
    z_ref = float(planes["z_ref"])
    positions, pose_params, labels = [], [], []
    for label, uv in site_uvs.items():
        for z_h in z_values:
            positions.append(
                transform_adsorbate(adsorbate, density, uv, float(z_h), quaternion, z_ref)
            )
            pose_params.append(np.concatenate([uv, [float(z_h)], quaternion]))
            labels.append(label)
    return PoseBatch(
        adsorbate=adsorbate,
        positions=np.asarray(positions, dtype=float),
        pose_params=np.asarray(pose_params, dtype=float),
        site_labels=labels,
        metadata={"z_ref": z_ref, "z_image": float(planes["z_image"]),
                  "scan_type": "z_scan"},
    )


def _site_uv_lookup(density, labels):
    lookup = {}
    for site in infer_surface_sites(density):
        if site.label in labels and site.label not in lookup:
            lookup[site.label] = np.asarray(site.uv, dtype=float)
    missing = [lab for lab in labels if lab not in lookup]
    if missing:
        raise ValueError(f"Site labels not found: {missing}")
    return {lab: lookup[lab] for lab in labels}


def _cartesian_to_fractional_uv(xy_cart, cell):
    """Convert Cartesian (x, y) to fractional (u, v) using the 2D cell."""
    cell2 = np.asarray(cell[:2, :2], dtype=float)
    inv_cell2 = np.linalg.inv(cell2)
    return np.asarray(xy_cart, dtype=float) @ inv_cell2.T


# =====================================================================
#  Phase 1: Generate reference dataset
# =====================================================================

def phase_generate(args, config, density, adsorbate, teacher_backend, out_dir):
    print("\n" + "=" * 60)
    print("  PHASE 1: Generate ultra-fine 6D reference dataset")
    print("=" * 60 + "\n")

    sampler_config = Sampler6DConfig(
        n_u=args.n_u,
        n_v=args.n_v,
        n_z=args.n_z,
        z_range=(args.z_min, args.z_max),
        z_bias_power=args.z_bias_power,
        n_orient=args.n_orient,
        tilt_max_deg=args.tilt_max_deg,
        n_yaw=args.n_yaw,
        include_high_symmetry_sites=True,
        random_fraction=args.random_fraction,
        seed=42,
    )

    print(f"[generate] Sampling config:")
    print(f"  Lateral: {args.n_u}x{args.n_v} = {args.n_u * args.n_v} points")
    print(f"  Height:  {args.n_z} points in [{args.z_min}, {args.z_max}] A (bias={args.z_bias_power})")
    print(f"  Orient:  {args.n_orient} (tilt_max={args.tilt_max_deg}°, n_yaw={args.n_yaw})")
    systematic = args.n_u * args.n_v * args.n_z * args.n_orient
    total_est = int(systematic * (1.0 + args.random_fraction))
    print(f"  Systematic: {systematic:,}  Random: {int(systematic * args.random_fraction):,}")
    print(f"  Total estimated: {total_est:,} poses")
    print(f"  Estimated GPU time: {total_est / 100.0 / 3600.0:.1f}–{total_est / 80.0 / 3600.0:.1f} hours")
    print(flush=True)

    t0 = time.perf_counter()
    dataset = generate_6d_reference_data(
        density=density,
        adsorbate=adsorbate,
        teacher_backend=teacher_backend,
        sampler_config=sampler_config,
        chunk_size=args.teacher_chunk_size,
        verbose=True,
    )
    elapsed = time.perf_counter() - t0

    print(f"\n[generate] Done! {dataset.n_poses:,} poses in {elapsed:.1f}s "
          f"({elapsed / 3600.0:.2f}h)")
    print(f"[generate] Energy range: [{dataset.energies.min():.4f}, {dataset.energies.max():.4f}] eV")

    save_6d_dataset(dataset, out_dir)
    dataset_path = out_dir / f"{adsorbate.name}_6d.npz"
    print(f"[generate] Saved to {dataset_path}")

    save_json({
        "phase": "generate",
        "n_poses": dataset.n_poses,
        "elapsed_seconds": elapsed,
        "elapsed_hours": elapsed / 3600.0,
        "energy_range_eV": [float(dataset.energies.min()), float(dataset.energies.max())],
        "sampler_config": sampler_config.__dict__,
    }, out_dir / "generate_summary.json")

    return dataset


# =====================================================================
#  Phase 2: Fit parameters
# =====================================================================

def _filter_by_energy(dataset: Dataset6D, e_max: float) -> Dataset6D:
    """Remove poses with teacher energy above e_max (eV).

    Collision poses from tilted molecules at low z can reach tens of eV and
    completely dominate the MSE loss, preventing the optimizer from fitting
    the physically relevant energy region.
    """
    energies = np.asarray(dataset.teacher.energies, dtype=float)
    mask = energies <= e_max
    n_keep = int(np.count_nonzero(mask))
    n_drop = len(energies) - n_keep
    idx = np.where(mask)[0]

    if n_drop > 0:
        print(f"[filter] Energy cutoff {e_max:.2f} eV: keeping {n_keep:,} / "
              f"{len(energies):,} poses ({n_drop:,} dropped, "
              f"max dropped E = {energies[~mask].max():.2f} eV)", flush=True)
        e_kept = energies[mask]
        print(f"[filter] Kept energy range: [{e_kept.min():.4f}, {e_kept.max():.4f}] eV", flush=True)

    return Dataset6D(
        adsorbate=dataset.adsorbate,
        poses=PoseBatch(
            adsorbate=dataset.adsorbate,
            positions=dataset.poses.positions[idx],
            pose_params=dataset.poses.pose_params[idx],
            site_labels=[dataset.poses.site_labels[int(i)] for i in idx],
            metadata=dict(dataset.poses.metadata),
        ),
        teacher=TeacherResult(
            energies=dataset.teacher.energies[idx],
            forces=dataset.teacher.forces[idx],
            metadata=dict(dataset.teacher.metadata),
        ),
        metadata={"energy_cutoff_eV": e_max, "n_dropped": n_drop,
                  "n_kept": n_keep, "original_n": len(energies)},
    )


def phase_fit(args, config, density, adsorbate, out_dir, dataset=None):
    print("\n" + "=" * 60)
    print("  PHASE 2: Fit GridFF parameters (energy + force)")
    print("=" * 60 + "\n")

    # Load dataset if not provided
    if dataset is None:
        dataset_path = out_dir / f"{adsorbate.name}_6d.npz"
        print(f"[fit] Loading dataset from {dataset_path} ...", flush=True)
        dataset = load_6d_dataset(dataset_path)
        print(f"[fit] Loaded {dataset.n_poses:,} poses", flush=True)

    # FIX #1: Filter out high-energy collision poses before fitting
    e_cutoff = float(args.energy_cutoff)
    raw_energies = np.asarray(dataset.teacher.energies, dtype=float)
    print(f"[fit] Raw energy range: [{raw_energies.min():.4f}, {raw_energies.max():.4f}] eV "
          f"({dataset.n_poses:,} poses)", flush=True)

    n_above = int(np.count_nonzero(raw_energies > e_cutoff))
    if n_above > 0:
        print(f"[fit] {n_above:,} poses ({100.0 * n_above / len(raw_energies):.1f}%) "
              f"above {e_cutoff:.1f} eV cutoff — filtering ...", flush=True)
        dataset = _filter_by_energy(dataset, e_cutoff)
    else:
        print(f"[fit] All poses within {e_cutoff:.1f} eV cutoff — no filtering needed", flush=True)

    # Split filtered dataset
    print("[fit] Splitting dataset (70/15/15) ...", flush=True)
    splits = split_6d_dataset(
        dataset,
        val_fraction=0.15,
        test_fraction=0.15,
        stratify_by_z=True,
        stratify_by_site=True,
        z_bins=8,
        seed=42,
    )
    print(f"[fit] Train: {splits.train.n_poses:,}  Val: {splits.val.n_poses:,}  "
          f"Test: {splits.test.n_poses:,}", flush=True)

    # Build model
    install_raw = bool(args.raw_r6 and args.two_stage_c6)
    model = _build_model(config, density, adsorbate,
                         prefer_jax=args.prefer_jax,
                         install_raw_r6=install_raw)

    # Subsample training data for Stage 1 gradient descent if needed.
    # L-BFGS-B uses finite differences → O(N*P) per iteration.
    # 14k poses × 8 params × 2 = ~230k model evaluations per L-BFGS step.
    # Subsampling to ~1000-1500 makes this tractable in minutes.
    max_train = int(args.max_train_poses)
    train_poses_s1 = splits.train.poses
    train_teacher_s1 = splits.train.teacher
    if splits.train.n_poses > max_train:
        rng = np.random.default_rng(42)
        idx = rng.choice(splits.train.n_poses, size=max_train, replace=False)
        idx = np.sort(idx)
        train_poses_s1 = _slice_pose_batch(splits.train.poses, idx)
        train_teacher_s1 = TeacherResult(
            energies=splits.train.teacher.energies[idx],
            forces=splits.train.teacher.forces[idx],
            metadata=dict(splits.train.teacher.metadata),
        )
        print(f"[fit] Subsampled Stage 1 training: {max_train:,} / "
              f"{splits.train.n_poses:,} poses", flush=True)

    train_data_s1 = [(train_poses_s1, train_teacher_s1)]
    # FIX #3: For lstsq Stage 2, use ALL filtered data (train+val+test) since
    # C₆ lstsq is a 5-parameter linear fit with negligible overfitting risk.
    # This matches the benchmark_ag_zscan.py approach.
    all_filtered_data = [(dataset.poses, dataset.teacher)]

    # Fit
    print("[fit] Starting parameter optimization ...", flush=True)
    t0 = time.perf_counter()

    if args.two_stage_c6:
        fit_result = fit_two_stage_c6(
            density=density,
            datasets=train_data_s1,
            model=model,
            training=config.training,
            initial_params=None,
            lstsq_datasets=all_filtered_data,
            use_raw_r6=args.raw_r6,
        )
    else:
        fit_result = fit_hybrid_parameters(
            density=density,
            datasets=train_data_s1,
            model=model,
            training=config.training,
            initial_params=None,
        )

    fit_elapsed = time.perf_counter() - t0
    print(f"[fit] Optimization done in {fit_elapsed:.1f}s", flush=True)

    # Save fitted params
    save_json(_fit_param_payload(fit_result.params), out_dir / "fit_params.json")
    save_json(fit_result.metrics.get("constraint_report", {}), out_dir / "fit_constraints.json")

    # Evaluate on all splits
    plots_dir = ensure_dir(out_dir / "plots")
    print("[fit] Evaluating student on all splits ...", flush=True)

    for split_name, split_ds in [("train", splits.train), ("val", splits.val), ("test", splits.test)]:
        student_eval, _ = _evaluate_student_chunked(
            model, split_ds.poses, fit_result.params, chunk_size=args.student_chunk_size
        )
        m = _metrics(split_ds.teacher.energies, student_eval.energies,
                     split_ds.teacher.forces, student_eval.forces)
        e_rmse = m["energy"]["rmse"] * 1000
        f_rmse = m["force"]["rmse"] * 1000
        print(f"  {split_name:5s}: E_RMSE = {e_rmse:.1f} meV,  F_RMSE = {f_rmse:.1f} meV/A", flush=True)

    # Parity plots on filtered dataset
    print("[fit] Generating parity plots ...", flush=True)
    full_student, _ = _evaluate_student_chunked(
        model, dataset.poses, fit_result.params, chunk_size=args.student_chunk_size
    )
    plot_parity(dataset.teacher.energies, full_student.energies, plots_dir, "energy_parity.png")
    plot_parity(dataset.teacher.forces.ravel(), full_student.forces.ravel(), plots_dir, "force_parity.png")
    plot_error_histogram(full_student.energies - dataset.teacher.energies, plots_dir, "energy_hist.png")
    plot_error_histogram(full_student.forces.ravel() - dataset.teacher.forces.ravel(), plots_dir, "force_hist.png")
    plot_convergence(fit_result.history, plots_dir, name_prefix="training_loss")

    # Export grid artifacts
    export_paths = export_firecore_artifacts(
        out_dir=out_dir,
        density=density,
        model=model,
        fit_result=fit_result,
        toggles=config.toggles,
        teacher_backend_id=config.teacher_backend.kind,
    )
    print(f"[fit] Grid exported: {export_paths['plq_path']}", flush=True)

    full_metrics = _metrics(dataset.teacher.energies, full_student.energies,
                            dataset.teacher.forces, full_student.forces)
    save_json({
        "phase": "fit",
        "fit_seconds": fit_elapsed,
        "energy_cutoff_eV": e_cutoff,
        "n_filtered": dataset.n_poses,
        "n_train": splits.train.n_poses,
        "n_val": splits.val.n_poses,
        "n_test": splits.test.n_poses,
        "full_dataset_metrics": full_metrics,
        "optimizer_metrics": fit_result.metrics,
    }, out_dir / "fit_summary.json")

    return fit_result, model


# =====================================================================
#  Phase 3: Validation scans
# =====================================================================

def phase_scan(args, config, density, adsorbate, teacher_backend, out_dir,
               fit_result=None, model=None):
    print("\n" + "=" * 60)
    print("  PHASE 3: Validation scans (z, x, y, xy)")
    print("=" * 60 + "\n")

    plots_dir = ensure_dir(out_dir / "plots")
    traces_dir = ensure_dir(out_dir / "traces")

    # Load or reuse params and model
    if fit_result is None:
        params = _load_fitted_params(out_dir / "fit_params.json")
        print("[scan] Loaded fitted params from disk", flush=True)
    else:
        params = fit_result.params

    if model is None:
        install_raw = bool(args.raw_r6 and args.two_stage_c6)
        model = _build_model(config, density, adsorbate,
                             prefer_jax=args.prefer_jax,
                             install_raw_r6=install_raw)

    cell = np.asarray(density.cell, dtype=float)
    quaternion = np.array([1.0, 0.0, 0.0, 0.0], dtype=float)  # upright
    site_uvs = _site_uv_lookup(density, ["top", "bridge", "hollow"])
    print(f"[scan] Sites: {  {k: v.tolist() for k, v in site_uvs.items()}  }", flush=True)

    scan_summary = {}

    # ── 3a. Z-scan ────────────────────────────────────────────────
    print("\n--- Z-scan (2–20 Å) at top/bridge/hollow ---", flush=True)
    # Biased z-sampling: more points near surface
    t_lin = np.linspace(0.0, 1.0, args.scan_z_points, dtype=float)
    z_values = args.scan_z_min + (args.scan_z_max - args.scan_z_min) * np.power(t_lin, 1.5)

    zscan_poses = _build_zscan_poses(density, adsorbate, site_uvs, z_values, quaternion)
    print(f"[scan] Z-scan: {len(zscan_poses.positions)} poses "
          f"({len(site_uvs)} sites × {len(z_values)} z-points)", flush=True)

    zscan_teacher, zt_elapsed = _evaluate_teacher_chunked(
        teacher_backend, density, zscan_poses, chunk_size=args.teacher_chunk_size
    )
    print(f"[scan] Z-scan teacher: {zt_elapsed:.1f}s", flush=True)

    zscan_student, zs_elapsed = _evaluate_student_chunked(
        model, zscan_poses, params, chunk_size=args.student_chunk_size
    )
    print(f"[scan] Z-scan student: {zs_elapsed:.3f}s", flush=True)

    nz = len(z_values)
    eq_z_height = None  # will auto-detect

    for i_site, site_label in enumerate(site_uvs.keys()):
        sl = slice(i_site * nz, (i_site + 1) * nz)
        te = np.asarray(zscan_teacher.energies[sl], dtype=float)
        se = np.asarray(zscan_student.energies[sl], dtype=float)
        tf = np.asarray(zscan_teacher.forces[sl], dtype=float)
        sf = np.asarray(zscan_student.forces[sl], dtype=float)

        residual = se - te
        e_rmse = float(np.sqrt(np.mean(residual ** 2))) * 1000
        e_mae = float(np.mean(np.abs(residual))) * 1000
        f_rmse_val = float(np.sqrt(np.mean((sf - tf) ** 2))) * 1000
        print(f"  {site_label}: E_RMSE={e_rmse:.1f} meV, E_MAE={e_mae:.1f} meV, F_RMSE={f_rmse_val:.1f} meV/A")

        # Auto-detect equilibrium z from top site
        if site_label == "top" and eq_z_height is None:
            eq_idx = int(np.argmin(te))
            eq_z_height = float(z_values[eq_idx])
            print(f"  [auto] Equilibrium z = {eq_z_height:.2f} A (from {site_label} minimum)")

        scan_summary[f"zscan_{site_label}"] = {
            "energy_rmse_meV": e_rmse,
            "energy_mae_meV": e_mae,
            "force_rmse_meV_per_A": f_rmse_val,
        }

        # Save trace
        np.savez_compressed(
            traces_dir / f"{site_label}_zscan_trace.npz",
            z_values=z_values, teacher_energy=te, student_energy=se,
            teacher_force=tf, student_force=sf,
        )

        # Plot energy trace
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), sharex=True)
        ax1.plot(z_values, te * 1000, label="MAD-SURF", color=CLR_TEACHER, lw=2.5)
        ax1.plot(z_values, se * 1000, label="GridFF", color=CLR_STUDENT, lw=2.0, ls="--")
        ax1.set_ylabel("Energy [meV]")
        ax1.set_title(f"{site_label} z-scan  |  E_RMSE = {e_rmse:.1f} meV")
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        ax2.axhline(0, color="black", lw=0.5, alpha=0.5)
        ax2.plot(z_values, residual * 1000, color=CLR_ERROR, lw=1.5)
        ax2.set_xlabel("z [A]")
        ax2.set_ylabel("Residual [meV]")
        ax2.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(plots_dir / f"{site_label}_zscan.png", dpi=170)
        plt.close(fig)

        # Plot force trace
        tf_norm = _force_norms(tf)
        sf_norm = _force_norms(sf)
        fig, ax = plt.subplots(figsize=(8, 4.5))
        ax.plot(z_values, tf_norm, label="MAD-SURF", color=CLR_TEACHER, lw=2.5)
        ax.plot(z_values, sf_norm, label="GridFF", color=CLR_STUDENT, lw=2.0, ls="--")
        ax.set_xlabel("z [A]")
        ax.set_ylabel("|Force| [eV/A]")
        ax.set_title(f"{site_label} z-scan force norm")
        ax.legend()
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(plots_dir / f"{site_label}_zscan_force.png", dpi=170)
        plt.close(fig)

    # Use auto-detected or user-specified z for lateral scans
    if args.scan_z_height is not None:
        eq_z_height = args.scan_z_height
    if eq_z_height is None:
        eq_z_height = 3.0  # fallback
    print(f"\n[scan] Lateral scan z-height = {eq_z_height:.2f} A", flush=True)

    # ── 3b. X-scan (Cartesian x = 0..6 Å) ────────────────────────
    print(f"\n--- X-scan (0–{args.scan_lateral_max} A along x) ---", flush=True)
    x_vals = np.linspace(0.0, args.scan_lateral_max, args.scan_lateral_points, dtype=float)
    xy_cart_x = np.column_stack([x_vals, np.zeros_like(x_vals)])
    uv_x = _cartesian_to_fractional_uv(xy_cart_x, cell)

    xscan_poses = _build_scan_poses(density, adsorbate, uv_x, eq_z_height, quaternion)
    xscan_teacher, _ = _evaluate_teacher_chunked(
        teacher_backend, density, xscan_poses, chunk_size=args.teacher_chunk_size
    )
    xscan_student, _ = _evaluate_student_chunked(
        model, xscan_poses, params, chunk_size=args.student_chunk_size
    )
    xscan_resid = np.asarray(xscan_student.energies, dtype=float) - np.asarray(xscan_teacher.energies, dtype=float)
    x_rmse = float(np.sqrt(np.mean(xscan_resid ** 2))) * 1000
    x_mae = float(np.mean(np.abs(xscan_resid))) * 1000
    print(f"  X-scan: E_RMSE = {x_rmse:.1f} meV, E_MAE = {x_mae:.1f} meV")
    scan_summary["xscan"] = {"energy_rmse_meV": x_rmse, "energy_mae_meV": x_mae}

    np.savez_compressed(traces_dir / "x_scan.npz",
                        x_vals=x_vals,
                        teacher_energy=xscan_teacher.energies,
                        student_energy=xscan_student.energies,
                        teacher_force=xscan_teacher.forces,
                        student_force=xscan_student.forces)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
    ax1.plot(x_vals, np.asarray(xscan_teacher.energies) * 1000, label="MAD-SURF", color=CLR_TEACHER, lw=2.5)
    ax1.plot(x_vals, np.asarray(xscan_student.energies) * 1000, label="GridFF", color=CLR_STUDENT, lw=2.0, ls="--")
    ax1.set_ylabel("Energy [meV]")
    ax1.set_title(f"X-scan at z={eq_z_height:.1f} A  |  RMSE={x_rmse:.1f} meV")
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax2.axhline(0, color="black", lw=0.5)
    ax2.plot(x_vals, xscan_resid * 1000, color=CLR_ERROR, lw=1.5)
    ax2.set_xlabel("x [A]")
    ax2.set_ylabel("Residual [meV]")
    ax2.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(plots_dir / "x_scan_energy.png", dpi=170)
    plt.close(fig)

    # X-scan force
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(x_vals, _force_norms(xscan_teacher.forces), label="MAD-SURF", color=CLR_TEACHER, lw=2.5)
    ax.plot(x_vals, _force_norms(xscan_student.forces), label="GridFF", color=CLR_STUDENT, lw=2.0, ls="--")
    ax.set_xlabel("x [A]")
    ax.set_ylabel("|Force| [eV/A]")
    ax.set_title(f"X-scan force at z={eq_z_height:.1f} A")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(plots_dir / "x_scan_force.png", dpi=170)
    plt.close(fig)

    # ── 3c. Y-scan (Cartesian y = 0..6 Å) ────────────────────────
    print(f"\n--- Y-scan (0–{args.scan_lateral_max} A along y) ---", flush=True)
    y_vals = np.linspace(0.0, args.scan_lateral_max, args.scan_lateral_points, dtype=float)
    xy_cart_y = np.column_stack([np.zeros_like(y_vals), y_vals])
    uv_y = _cartesian_to_fractional_uv(xy_cart_y, cell)

    yscan_poses = _build_scan_poses(density, adsorbate, uv_y, eq_z_height, quaternion)
    yscan_teacher, _ = _evaluate_teacher_chunked(
        teacher_backend, density, yscan_poses, chunk_size=args.teacher_chunk_size
    )
    yscan_student, _ = _evaluate_student_chunked(
        model, yscan_poses, params, chunk_size=args.student_chunk_size
    )
    yscan_resid = np.asarray(yscan_student.energies, dtype=float) - np.asarray(yscan_teacher.energies, dtype=float)
    y_rmse = float(np.sqrt(np.mean(yscan_resid ** 2))) * 1000
    y_mae = float(np.mean(np.abs(yscan_resid))) * 1000
    print(f"  Y-scan: E_RMSE = {y_rmse:.1f} meV, E_MAE = {y_mae:.1f} meV")
    scan_summary["yscan"] = {"energy_rmse_meV": y_rmse, "energy_mae_meV": y_mae}

    np.savez_compressed(traces_dir / "y_scan.npz",
                        y_vals=y_vals,
                        teacher_energy=yscan_teacher.energies,
                        student_energy=yscan_student.energies,
                        teacher_force=yscan_teacher.forces,
                        student_force=yscan_student.forces)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
    ax1.plot(y_vals, np.asarray(yscan_teacher.energies) * 1000, label="MAD-SURF", color=CLR_TEACHER, lw=2.5)
    ax1.plot(y_vals, np.asarray(yscan_student.energies) * 1000, label="GridFF", color=CLR_STUDENT, lw=2.0, ls="--")
    ax1.set_ylabel("Energy [meV]")
    ax1.set_title(f"Y-scan at z={eq_z_height:.1f} A  |  RMSE={y_rmse:.1f} meV")
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax2.axhline(0, color="black", lw=0.5)
    ax2.plot(y_vals, yscan_resid * 1000, color=CLR_ERROR, lw=1.5)
    ax2.set_xlabel("y [A]")
    ax2.set_ylabel("Residual [meV]")
    ax2.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(plots_dir / "y_scan_energy.png", dpi=170)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(y_vals, _force_norms(yscan_teacher.forces), label="MAD-SURF", color=CLR_TEACHER, lw=2.5)
    ax.plot(y_vals, _force_norms(yscan_student.forces), label="GridFF", color=CLR_STUDENT, lw=2.0, ls="--")
    ax.set_xlabel("y [A]")
    ax.set_ylabel("|Force| [eV/A]")
    ax.set_title(f"Y-scan force at z={eq_z_height:.1f} A")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(plots_dir / "y_scan_force.png", dpi=170)
    plt.close(fig)

    # ── 3d. XY-scan (2D heatmap, 0..6 x 0..6 Å) ─────────────────
    print(f"\n--- XY-scan ({args.scan_xy_points}x{args.scan_xy_points}, "
          f"0–{args.scan_lateral_max} x 0–{args.scan_lateral_max} A) ---", flush=True)
    nx = ny = args.scan_xy_points
    x_grid = np.linspace(0.0, args.scan_lateral_max, nx, dtype=float)
    y_grid = np.linspace(0.0, args.scan_lateral_max, ny, dtype=float)
    xx, yy = np.meshgrid(x_grid, y_grid, indexing="ij")
    xy_flat = np.column_stack([xx.ravel(), yy.ravel()])
    uv_flat = _cartesian_to_fractional_uv(xy_flat, cell)

    xyscan_poses = _build_scan_poses(density, adsorbate, uv_flat, eq_z_height, quaternion)
    print(f"[scan] XY-scan: {len(xyscan_poses.positions)} poses", flush=True)

    xyscan_teacher, xyt_elapsed = _evaluate_teacher_chunked(
        teacher_backend, density, xyscan_poses, chunk_size=args.teacher_chunk_size
    )
    print(f"[scan] XY teacher: {xyt_elapsed:.1f}s", flush=True)

    xyscan_student, xys_elapsed = _evaluate_student_chunked(
        model, xyscan_poses, params, chunk_size=args.student_chunk_size
    )
    print(f"[scan] XY student: {xys_elapsed:.3f}s", flush=True)

    te_2d = np.asarray(xyscan_teacher.energies, dtype=float).reshape(nx, ny)
    se_2d = np.asarray(xyscan_student.energies, dtype=float).reshape(nx, ny)
    err_2d = (se_2d - te_2d)
    xy_rmse = float(np.sqrt(np.mean(err_2d ** 2))) * 1000
    xy_mae = float(np.mean(np.abs(err_2d))) * 1000
    print(f"  XY-scan: E_RMSE = {xy_rmse:.1f} meV, E_MAE = {xy_mae:.1f} meV")
    scan_summary["xyscan"] = {"energy_rmse_meV": xy_rmse, "energy_mae_meV": xy_mae}

    np.savez_compressed(traces_dir / "xy_2d_grid.npz",
                        x_grid=x_grid, y_grid=y_grid,
                        teacher_energy_2d=te_2d, student_energy_2d=se_2d)

    # 2D heatmap plots
    vmin = min(te_2d.min(), se_2d.min()) * 1000
    vmax = max(te_2d.max(), se_2d.max()) * 1000
    err_abs_max = max(abs(err_2d.min()), abs(err_2d.max())) * 1000

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), dpi=170)
    titles = ["MAD-SURF", "GridFF", "Error (GridFF - MAD-SURF)"]
    data = [te_2d * 1000, se_2d * 1000, err_2d * 1000]
    cmaps = ["viridis", "viridis", "RdBu_r"]
    vmins_list = [vmin, vmin, -err_abs_max]
    vmaxs_list = [vmax, vmax, err_abs_max]

    for i, ax in enumerate(axes):
        im = ax.pcolormesh(x_grid, y_grid, data[i].T, cmap=cmaps[i],
                          vmin=vmins_list[i], vmax=vmaxs_list[i], shading="auto")
        cb = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.03)
        cb.set_label("meV", fontsize=9)
        ax.set_xlabel("x [A]")
        if i == 0:
            ax.set_ylabel("y [A]")
        ax.set_title(titles[i], fontsize=10)
        ax.set_aspect("equal")

    fig.suptitle(
        f"CHONH2 / Ag(111) XY energy heatmap  z={eq_z_height:.1f} A"
        f"  |  RMSE={xy_rmse:.1f} meV",
        fontsize=11, y=1.02,
    )
    fig.tight_layout()
    fig.savefig(plots_dir / "xy_heatmap_energy.png", dpi=170, bbox_inches="tight")
    plt.close(fig)

    # XY force heatmap
    tf_norms_2d = _force_norms(xyscan_teacher.forces).reshape(nx, ny)
    sf_norms_2d = _force_norms(xyscan_student.forces).reshape(nx, ny)
    ferr_2d = sf_norms_2d - tf_norms_2d
    fvmin = min(tf_norms_2d.min(), sf_norms_2d.min())
    fvmax = max(tf_norms_2d.max(), sf_norms_2d.max())
    ferr_max = max(abs(ferr_2d.min()), abs(ferr_2d.max()))

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), dpi=170)
    for i, (d, cm, lo, hi, ttl) in enumerate([
        (tf_norms_2d, "viridis", fvmin, fvmax, "MAD-SURF |F|"),
        (sf_norms_2d, "viridis", fvmin, fvmax, "GridFF |F|"),
        (ferr_2d, "RdBu_r", -ferr_max, ferr_max, "Error"),
    ]):
        im = axes[i].pcolormesh(x_grid, y_grid, d.T, cmap=cm, vmin=lo, vmax=hi, shading="auto")
        cb = fig.colorbar(im, ax=axes[i], shrink=0.85, pad=0.03)
        cb.set_label("eV/A", fontsize=9)
        axes[i].set_xlabel("x [A]")
        if i == 0:
            axes[i].set_ylabel("y [A]")
        axes[i].set_title(ttl, fontsize=10)
        axes[i].set_aspect("equal")

    fig.suptitle(
        f"CHONH2 / Ag(111) XY force heatmap  z={eq_z_height:.1f} A",
        fontsize=11, y=1.02,
    )
    fig.tight_layout()
    fig.savefig(plots_dir / "xy_heatmap_force.png", dpi=170, bbox_inches="tight")
    plt.close(fig)

    # Save scan summary
    scan_summary["equilibrium_z_A"] = eq_z_height
    save_json(scan_summary, out_dir / "scan_summary.json")

    print("\n" + "=" * 60)
    print("  SCAN RESULTS SUMMARY")
    print("=" * 60)
    for key, val in scan_summary.items():
        if isinstance(val, dict):
            rmse = val.get("energy_rmse_meV", "?")
            mae = val.get("energy_mae_meV", "?")
            print(f"  {key:20s}: E_RMSE = {rmse:.1f} meV, E_MAE = {mae:.1f} meV")
    print("=" * 60)

    return scan_summary


# =====================================================================
#  Main
# =====================================================================

def main():
    args = parse_args()
    out_dir = ensure_dir(args.out_dir)

    print("=" * 60)
    print("  CHONH2 / Ag(111) End-to-End GridFF Benchmark")
    print("=" * 60)
    print(f"  Mode:       {args.mode}")
    print(f"  Output:     {out_dir}")
    print(f"  CHGCAR:     {args.chgcar}")
    print(f"  LOCPOT:     {args.locpot}")
    print(f"  MACE model: {args.model_path}")
    print(f"  Device:     {args.device}")
    print(f"  Two-stage:  {args.two_stage_c6}")
    print(f"  Raw r6:     {args.raw_r6}")
    print("=" * 60)
    print(flush=True)

    # Build config
    config = _setup_config(args)
    save_config(config, out_dir / "config.json")

    # Load density
    print("[e2e] Loading density (CHGCAR + LOCPOT) ...", flush=True)
    density = make_density_backend(
        config.density_backend, surface=config.surface, grid=config.grid
    ).load()
    print(f"[e2e] Cell: {density.cell[0,0]:.3f} x {density.cell[1,1]:.3f} x {density.cell[2,2]:.3f} A")
    print(f"[e2e] Grid: {density.grid_shape_xyz_hint}", flush=True)

    # Load adsorbate
    adsorbate = get_adsorbate(
        name=args.adsorbate_name,
        xyz_path=args.xyz_path,
        anchor_index=args.anchor_index,
        use_input_charges=args.use_input_charges,
    )
    print(f"[e2e] Adsorbate: {adsorbate.name} ({adsorbate.natoms} atoms: {adsorbate.symbols})")

    # Load teacher
    need_teacher = args.mode in ("generate", "scan", "all")
    teacher_backend = None
    if need_teacher:
        print("[e2e] Loading MAD-SURF teacher backend ...", flush=True)
        teacher_backend = make_teacher_backend(config.teacher_backend)
        print("[e2e] Teacher ready", flush=True)

    # Dispatch
    dataset = None
    fit_result = None
    model = None

    if args.mode in ("generate", "all"):
        dataset = phase_generate(args, config, density, adsorbate, teacher_backend, out_dir)

    if args.mode in ("fit", "all"):
        fit_result, model = phase_fit(args, config, density, adsorbate, out_dir, dataset=dataset)

    if args.mode in ("scan", "all"):
        phase_scan(args, config, density, adsorbate, teacher_backend, out_dir,
                   fit_result=fit_result, model=model)

    print(f"\n[e2e] All done! Results in: {out_dir}")


if __name__ == "__main__":
    main()
