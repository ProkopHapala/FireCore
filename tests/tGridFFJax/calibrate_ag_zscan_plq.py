#!/usr/bin/env python3
"""Calibrate Ag strict-PLQ substrate builder parameters against a focused CHON z-scan benchmark."""

from __future__ import annotations

import argparse
import itertools
from pathlib import Path
import sys
import time

import matplotlib
matplotlib.use("Agg")
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.append(str(ROOT))

from pyBall.gridff_jax import RunConfig, backend_summary, save_config
from pyBall.gridff_jax.density_backends import make_density_backend
from pyBall.gridff_jax.fit import fit_hybrid_parameters
from pyBall.gridff_jax.hybrid_energy import HybridGridFFModel
from pyBall.gridff_jax.interfaces import ModelEvaluation, PoseBatch, TeacherResult
from pyBall.gridff_jax.pose_sampling import get_adsorbate, infer_reference_planes, infer_surface_sites, transform_adsorbate
from pyBall.gridff_jax.teacher_backends import make_teacher_backend
from pyBall.gridff_jax.utils import ensure_dir, save_json
from pyBall.gridff_jax.validation import compute_error_metrics


DEFAULT_CHGCAR = "/home/niel/git/ORR_HER_Ag_Colab/results/Ag_ORR_HER/slab_clean/final_scf_12x12x1/CHGCAR"
DEFAULT_LOCPOT = "/home/niel/git/ORR_HER_Ag_Colab/results/Ag_ORR_HER/slab_clean/workfunc_12x12x1/LOCPOT"
DEFAULT_MADSURF_MODEL = str((Path(__file__).resolve().parent / "mad-surf_data" / "models" / "full_dataset_config_weights" / "MACE_model.model"))
def parse_args():
    parser = argparse.ArgumentParser(description="Calibrate Ag substrate PLQ builder parameters on a focused CHON z scan")
    parser.add_argument("--chgcar", default=DEFAULT_CHGCAR)
    parser.add_argument("--locpot", default=DEFAULT_LOCPOT)
    parser.add_argument("--model-path", default=DEFAULT_MADSURF_MODEL)
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--out-dir", default=str(ROOT / "tests/tGridFFJax/ag_zscan_calibration"))
    parser.add_argument("--xyz-path", default=None)
    parser.add_argument("--adsorbate-name", default="CHONH2")
    parser.add_argument("--anchor-index", type=int, default=0)
    parser.add_argument("--use-input-charges", action="store_true")
    parser.add_argument("--sites", default="top,bridge,hollow")
    parser.add_argument("--z-min", type=float, default=2.0)
    parser.add_argument("--stress-z-min", type=float, default=1.8)
    parser.add_argument("--z-max", type=float, default=5.60)
    parser.add_argument("--eval-z-points", type=int, default=121)
    parser.add_argument("--train-low-count", type=int, default=18)
    parser.add_argument("--train-high-count", type=int, default=12)
    parser.add_argument("--yaw-deg", type=float, default=0.0)
    parser.add_argument("--tilt-x-deg", type=float, default=0.0)
    parser.add_argument("--tilt-y-deg", type=float, default=0.0)
    parser.add_argument("--fit-static-charge", dest="fit_static_charge", action="store_true")
    parser.add_argument("--no-fit-static-charge", dest="fit_static_charge", action="store_false")
    parser.set_defaults(fit_static_charge=False)
    parser.add_argument("--calibration-mode", choices=("req_only", "power_sweep"), default="req_only")
    parser.add_argument("--alpha-values", default="1.5")
    parser.add_argument("--radius-scales", default="1.0")
    parser.add_argument("--energy-scales", default="1.0")
    parser.add_argument("--pauli-powers", default="1.0", help="Comma-separated metal_density_pauli_power values to sweep")
    parser.add_argument("--london-powers", default="0.5", help="Comma-separated metal_density_london_power values to sweep")
    parser.add_argument("--z-dense-cutoff", type=float, default=2.8, help="Upper edge of the dense chemisorption window")
    parser.add_argument("--max-steps", type=int, default=300)
    parser.add_argument("--force-weight", type=float, default=10.0)
    parser.add_argument("--learning-rate", type=float, default=1.0e-2)
    parser.add_argument("--teacher-chunk-size", type=int, default=64)
    parser.add_argument("--student-chunk-size", type=int, default=64)
    parser.add_argument("--w-holdout-force", type=float, default=1.0)
    parser.add_argument("--w-holdout-energy", type=float, default=0.35)
    parser.add_argument("--w-gap-energy", type=float, default=0.10)
    parser.add_argument("--prefer-jax", dest="prefer_jax", action="store_true")
    parser.add_argument("--no-prefer-jax", dest="prefer_jax", action="store_false")
    parser.set_defaults(prefer_jax=True)
    return parser.parse_args()


def _parse_values(text: str):
    return [float(item.strip()) for item in text.split(",") if item.strip()]


def _rotation_quaternion(axis, angle_rad: float):
    axis = np.asarray(axis, dtype=float)
    norm = np.linalg.norm(axis)
    if norm < 1.0e-12:
        return np.array([1.0, 0.0, 0.0, 0.0], dtype=float)
    axis = axis / norm
    half = 0.5 * float(angle_rad)
    return np.array(
        [
            np.cos(half),
            axis[0] * np.sin(half),
            axis[1] * np.sin(half),
            axis[2] * np.sin(half),
        ],
        dtype=float,
    )


def _quaternion_multiply(q2, q1):
    w2, x2, y2, z2 = np.asarray(q2, dtype=float)
    w1, x1, y1, z1 = np.asarray(q1, dtype=float)
    return np.array(
        [
            w2 * w1 - x2 * x1 - y2 * y1 - z2 * z1,
            w2 * x1 + x2 * w1 + y2 * z1 - z2 * y1,
            w2 * y1 - x2 * z1 + y2 * w1 + z2 * x1,
            w2 * z1 + x2 * y1 - y2 * x1 + z2 * w1,
        ],
        dtype=float,
    )


def _compose_fixed_quaternion(yaw_deg: float, tilt_x_deg: float, tilt_y_deg: float):
    q_yaw = _rotation_quaternion([0.0, 0.0, 1.0], np.deg2rad(float(yaw_deg)))
    q_tilt_x = _rotation_quaternion([1.0, 0.0, 0.0], np.deg2rad(float(tilt_x_deg)))
    q_tilt_y = _rotation_quaternion([0.0, 1.0, 0.0], np.deg2rad(float(tilt_y_deg)))
    return _quaternion_multiply(q_yaw, _quaternion_multiply(q_tilt_x, q_tilt_y))


def _site_uv_lookup(density, requested_labels):
    requested = [label.strip() for label in requested_labels.split(",") if label.strip()]
    lookup = {}
    for site in infer_surface_sites(density):
        if site.label in requested and site.label not in lookup:
            lookup[site.label] = np.asarray(site.uv, dtype=float)
    missing = [label for label in requested if label not in lookup]
    if missing:
        raise ValueError(f"Requested site labels not found on surface: {missing}")
    return {label: lookup[label] for label in requested}


def _build_zscan_pose_batch(density, adsorbate, site_uvs, z_values, quaternion):
    planes = infer_reference_planes(density)
    z_ref = float(planes["z_ref"])
    positions = []
    pose_params = []
    site_labels = []
    for label, uv in site_uvs.items():
        for z_height in z_values:
            positions.append(transform_adsorbate(adsorbate, density, uv, float(z_height), quaternion, z_ref))
            pose_params.append(np.concatenate([uv, np.array([float(z_height)], dtype=float), quaternion]))
            site_labels.append(label)
    return PoseBatch(
        adsorbate=adsorbate,
        positions=np.asarray(positions, dtype=float),
        pose_params=np.asarray(pose_params, dtype=float),
        site_labels=site_labels,
        metadata={
            "z_ref": z_ref,
            "z_image": float(planes["z_image"]),
            "surface_sites": [{"label": label, "uv": uv.tolist()} for label, uv in site_uvs.items()],
            "scan_type": "z_scan",
            "scan_z_values": np.asarray(z_values, dtype=float).tolist(),
            "orientation_quaternion": np.asarray(quaternion, dtype=float).tolist(),
        },
    )


def _slice_pose_batch(poses: PoseBatch, indices):
    indices = np.asarray(indices, dtype=int)
    return PoseBatch(
        adsorbate=poses.adsorbate,
        positions=np.asarray(poses.positions[indices], dtype=float),
        pose_params=np.asarray(poses.pose_params[indices], dtype=float),
        site_labels=[poses.site_labels[int(i)] for i in indices],
        metadata=dict(poses.metadata),
    )


def _slice_teacher_result(teacher: TeacherResult, indices):
    indices = np.asarray(indices, dtype=int)
    return TeacherResult(
        energies=np.asarray(teacher.energies[indices], dtype=float),
        forces=np.asarray(teacher.forces[indices], dtype=float),
        metadata=dict(teacher.metadata),
    )


def _evaluate_teacher_chunked(teacher_backend, density, poses, chunk_size: int):
    if len(poses.positions) <= chunk_size:
        t0 = time.perf_counter()
        result = teacher_backend.evaluate_batch(density, poses)
        elapsed = time.perf_counter() - t0
        return result, elapsed
    energies = []
    forces = []
    t0 = time.perf_counter()
    for start in range(0, len(poses.positions), int(chunk_size)):
        stop = min(start + int(chunk_size), len(poses.positions))
        partial = teacher_backend.evaluate_batch(density, _slice_pose_batch(poses, np.arange(start, stop, dtype=int)))
        energies.append(np.asarray(partial.energies, dtype=float))
        forces.append(np.asarray(partial.forces, dtype=float))
    elapsed = time.perf_counter() - t0
    return TeacherResult(
        energies=np.concatenate(energies, axis=0),
        forces=np.concatenate(forces, axis=0),
        metadata={
            "teacher_eval_seconds": elapsed,
            "seconds_per_pose": elapsed / max(len(poses.positions), 1),
        },
    ), elapsed


def _evaluate_student_chunked(model, poses, params, chunk_size: int):
    if len(poses.positions) <= chunk_size:
        return model.evaluate_batch(poses, params=params, compute_forces=True)
    energies = []
    forces = []
    components = {}
    charges = []
    for start in range(0, len(poses.positions), int(chunk_size)):
        stop = min(start + int(chunk_size), len(poses.positions))
        partial = model.evaluate_batch(_slice_pose_batch(poses, np.arange(start, stop, dtype=int)), params=params, compute_forces=True)
        energies.append(np.asarray(partial.energies, dtype=float))
        forces.append(np.asarray(partial.forces, dtype=float))
        for key, value in partial.components.items():
            components.setdefault(key, []).append(np.asarray(value, dtype=float))
        if partial.charges is not None:
            charges.append(np.asarray(partial.charges, dtype=float))
    return ModelEvaluation(
        energies=np.concatenate(energies, axis=0),
        forces=np.concatenate(forces, axis=0),
        components={key: np.concatenate(value, axis=0) for key, value in components.items()},
        charges=np.concatenate(charges, axis=0) if charges else None,
    )


def _set_strict_plq_config(config: RunConfig, args):
    config.density_backend.chgcar_path = args.chgcar
    config.density_backend.locpot_path = args.locpot
    config.teacher_backend.kind = "madsurf"
    config.teacher_backend.model_path = args.model_path
    config.teacher_backend.device = args.device
    config.grid.builder_mode = "metal_density_plq"
    config.grid.interpolation_order = 3
    config.toggles.use_ct_qeq = False
    config.toggles.use_image_charge = False
    config.toggles.use_reactive_grid = False
    config.hybrid_model.use_qeq = False
    config.hybrid_model.use_image = False
    config.hybrid_model.use_reactive_grid = False
    config.hybrid_model.use_req_plq = True
    config.training.fit_chi = False
    config.training.fit_hardness = False
    config.training.fit_image_plane = False
    config.training.fit_reactive = False
    config.training.fit_static_charge = bool(args.fit_static_charge)
    config.training.fit_req_radius_offset = True
    config.training.fit_req_energy_scale = True
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


def _select_even_indices(indices, count: int):
    indices = np.asarray(indices, dtype=int)
    if len(indices) == 0 or count <= 0:
        return np.zeros(0, dtype=int)
    if len(indices) <= count:
        return indices
    coords = np.linspace(0, len(indices) - 1, int(count), dtype=float)
    return np.unique(indices[np.round(coords).astype(int)])


def _train_indices_per_site(n_sites: int, nz: int, z_values, primary_z_min: float, z_dense_cutoff: float, z_max: float, n_train_low: int, n_train_high: int):
    z_values = np.asarray(z_values, dtype=float)
    low_positions = np.where((z_values >= float(primary_z_min)) & (z_values < float(z_dense_cutoff)))[0]
    high_positions = np.where((z_values >= float(z_dense_cutoff)) & (z_values <= float(z_max)))[0]
    local = np.unique(
        np.concatenate(
            [
                _select_even_indices(low_positions, int(n_train_low)),
                _select_even_indices(high_positions, int(n_train_high)),
            ]
        )
    )
    indices = []
    for isite in range(n_sites):
        offset = isite * nz
        indices.extend((offset + local).tolist())
    return np.asarray(indices, dtype=int)


def _window_masks(z_values, primary_z_min: float, stress_z_min: float, z_max: float):
    z_values = np.asarray(z_values, dtype=float)
    full = (z_values >= float(stress_z_min)) & (z_values <= float(z_max))
    primary = (z_values >= float(primary_z_min)) & (z_values <= float(z_max))
    stress = (z_values >= float(stress_z_min)) & (z_values < float(primary_z_min))
    return {"primary": primary, "stress": stress, "full": full}


def _metrics_for_mask(teacher_energy, student_energy, teacher_force, student_force, mask):
    mask = np.asarray(mask, dtype=bool)
    n_samples = int(np.count_nonzero(mask))
    if n_samples == 0:
        return {
            "n_samples": 0,
            "energy": {"rmse": 0.0, "mae": 0.0, "max_abs": 0.0, "mean_signed": 0.0},
            "force_components": {"rmse": 0.0, "mae": 0.0, "max_abs": 0.0, "mean_signed": 0.0},
            "force_norm": {"rmse": 0.0, "mae": 0.0, "max_abs": 0.0, "mean_signed": 0.0},
        }
    e_true = np.asarray(teacher_energy, dtype=float)[mask]
    e_pred = np.asarray(student_energy, dtype=float)[mask]
    f_true = np.asarray(teacher_force, dtype=float)[mask]
    f_pred = np.asarray(student_force, dtype=float)[mask]
    f_true_flat = f_true.reshape(-1)
    f_pred_flat = f_pred.reshape(-1)
    f_err_norm = np.linalg.norm((f_pred - f_true).reshape(len(f_true), -1), axis=1) if len(f_true) > 0 else np.zeros(0, dtype=float)
    return {
        "n_samples": n_samples,
        "energy": compute_error_metrics(e_true, e_pred),
        "force_components": compute_error_metrics(f_true_flat, f_pred_flat),
        "force_norm": compute_error_metrics(np.zeros_like(f_err_norm), f_err_norm),
    }


def _fit_param_payload(params):
    return {
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
    }


def _window_metrics(z_values, teacher_energy, student_energy, teacher_force, student_force, train_mask, primary_z_min: float, stress_z_min: float, z_max: float):
    train_mask = np.asarray(train_mask, dtype=bool)
    windows = _window_masks(z_values, primary_z_min=primary_z_min, stress_z_min=stress_z_min, z_max=z_max)
    return {
        "primary": {
            "all": _metrics_for_mask(teacher_energy, student_energy, teacher_force, student_force, windows["primary"]),
            "train": _metrics_for_mask(teacher_energy, student_energy, teacher_force, student_force, windows["primary"] & train_mask),
            "holdout": _metrics_for_mask(teacher_energy, student_energy, teacher_force, student_force, windows["primary"] & (~train_mask)),
        },
        "stress": {
            "all": _metrics_for_mask(teacher_energy, student_energy, teacher_force, student_force, windows["stress"]),
        },
        "full": {
            "all": _metrics_for_mask(teacher_energy, student_energy, teacher_force, student_force, windows["full"]),
        },
    }


def _trial_objective(window_metrics, args):
    primary_holdout = window_metrics["primary"]["holdout"]
    primary_train = window_metrics["primary"]["train"]
    return float(
        args.w_holdout_force * primary_holdout["force_norm"]["rmse"]
        + args.w_holdout_energy * primary_holdout["energy"]["rmse"]
        + args.w_gap_energy * abs(primary_train["energy"]["rmse"] - primary_holdout["energy"]["rmse"])
    )


def main():
    args = parse_args()
    out_dir = ensure_dir(args.out_dir)

    base_config = RunConfig()
    _set_strict_plq_config(base_config, args)
    save_config(base_config, out_dir / "config.json")

    base_density = make_density_backend(base_config.density_backend, surface=base_config.surface, grid=base_config.grid).load()
    teacher_backend = make_teacher_backend(base_config.teacher_backend)
    adsorbate = get_adsorbate(
        name=args.adsorbate_name,
        xyz_path=args.xyz_path,
        anchor_index=args.anchor_index,
        use_input_charges=args.use_input_charges,
    )
    site_uvs = _site_uv_lookup(base_density, args.sites)
    z_values = np.linspace(float(args.stress_z_min), float(args.z_max), int(args.eval_z_points), dtype=float)
    quaternion = _compose_fixed_quaternion(args.yaw_deg, args.tilt_x_deg, args.tilt_y_deg)
    dense_poses = _build_zscan_pose_batch(base_density, adsorbate, site_uvs, z_values, quaternion)

    teacher_dense, teacher_elapsed = _evaluate_teacher_chunked(
        teacher_backend=teacher_backend,
        density=base_density,
        poses=dense_poses,
        chunk_size=args.teacher_chunk_size,
    )
    save_json(
        {
            "teacher_seconds_total": float(teacher_elapsed),
            "teacher_seconds_per_pose": float(teacher_elapsed / max(len(dense_poses.positions), 1)),
            "backend": backend_summary(),
            "n_total": int(len(dense_poses.positions)),
        },
        out_dir / "teacher_cache_summary.json",
    )

    train_indices = _train_indices_per_site(
        len(site_uvs),
        len(z_values),
        z_values=z_values,
        primary_z_min=float(args.z_min),
        z_dense_cutoff=float(args.z_dense_cutoff),
        z_max=float(args.z_max),
        n_train_low=int(args.train_low_count),
        n_train_high=int(args.train_high_count),
    )
    train_poses = _slice_pose_batch(dense_poses, train_indices)
    train_teacher = _slice_teacher_result(teacher_dense, train_indices)
    train_mask = np.zeros(len(dense_poses.positions), dtype=bool)
    train_mask[train_indices] = True

    alpha_values = _parse_values(args.alpha_values)
    radius_scales = _parse_values(args.radius_scales)
    energy_scales = _parse_values(args.energy_scales)
    pauli_powers = _parse_values(args.pauli_powers)
    london_powers = _parse_values(args.london_powers)
    if args.calibration_mode == "req_only":
        pauli_powers = [float(pauli_powers[0])]
        london_powers = [float(london_powers[0])]
        alpha_values = [float(alpha_values[0])]
        radius_scales = [float(radius_scales[0])]
        energy_scales = [float(energy_scales[0])]
    trials = []
    best_record = None
    best_config = None
    reactive_elements = tuple(sorted(set(adsorbate.symbols)))

    for alpha_morse, radius_scale, energy_scale, pauli_power, london_power in itertools.product(
        alpha_values, radius_scales, energy_scales, pauli_powers, london_powers
    ):
        trial_config = RunConfig()
        _set_strict_plq_config(trial_config, args)
        trial_config.grid.alpha_morse = float(alpha_morse)
        trial_config.grid.req_scale_radius["Ag"] = float(radius_scale)
        trial_config.grid.req_scale_energy["Ag"] = float(energy_scale)
        trial_config.grid.metal_density_pauli_power = float(pauli_power)
        trial_config.grid.metal_density_london_power = float(london_power)

        density = make_density_backend(trial_config.density_backend, surface=trial_config.surface, grid=trial_config.grid).load()
        model = HybridGridFFModel(
            density=density,
            reactive_elements=reactive_elements,
            toggles=trial_config.toggles,
            grid_config=trial_config.grid,
            hybrid_config=trial_config.hybrid_model,
            prefer_jax=args.prefer_jax,
        )
        fit_t0 = time.perf_counter()
        fit_result = fit_hybrid_parameters(
            density=density,
            datasets=[(train_poses, train_teacher)],
            model=model,
            training=trial_config.training,
            initial_params=None,
        )
        fit_elapsed = time.perf_counter() - fit_t0
        pred_dense = _evaluate_student_chunked(
            model=model,
            poses=dense_poses,
            params=fit_result.params,
            chunk_size=args.student_chunk_size,
        )

        teacher_energy = np.asarray(teacher_dense.energies, dtype=float)
        student_energy = np.asarray(pred_dense.energies, dtype=float)
        teacher_force = np.asarray(teacher_dense.forces, dtype=float)
        student_force = np.asarray(pred_dense.forces, dtype=float)
        split_metrics = {
            "train": _metrics_for_mask(teacher_energy, student_energy, teacher_force, student_force, train_mask),
            "holdout": _metrics_for_mask(teacher_energy, student_energy, teacher_force, student_force, ~train_mask),
            "all": _metrics_for_mask(teacher_energy, student_energy, teacher_force, student_force, np.ones_like(train_mask, dtype=bool)),
        }
        windows = _window_metrics(
            z_values=dense_poses.pose_params[:, 2],
            teacher_energy=teacher_energy,
            student_energy=student_energy,
            teacher_force=teacher_force,
            student_force=student_force,
            train_mask=train_mask,
            primary_z_min=float(args.z_min),
            stress_z_min=float(args.stress_z_min),
            z_max=float(args.z_max),
        )
        objective = _trial_objective(windows, args)
        record = {
            "calibration_mode": args.calibration_mode,
            "alpha_morse": float(alpha_morse),
            "ag_radius_scale": float(radius_scale),
            "ag_energy_scale": float(energy_scale),
            "pauli_power": float(pauli_power),
            "london_power": float(london_power),
            "objective": objective,
            "timing": {
                "fit_seconds_total": float(fit_elapsed),
            },
            "split_metrics": split_metrics,
            "window_metrics": windows,
            "fit_params": _fit_param_payload(fit_result.params),
            "constraint_report": fit_result.metrics.get("constraint_report", {}),
            "optimizer_result": fit_result.optimizer_result,
        }
        trials.append(record)
        if (best_record is None) or (record["objective"] < best_record["objective"]):
            best_record = record
            best_config = trial_config
            save_config(best_config, out_dir / "best_config.json")
            save_json(best_record, out_dir / "best_result.json")

    trials.sort(key=lambda item: item["objective"])
    save_json(
        {
            "teacher_cache_summary": {
                "teacher_seconds_total": float(teacher_elapsed),
                "teacher_seconds_per_pose": float(teacher_elapsed / max(len(dense_poses.positions), 1)),
            },
            "objective_weights": {
                "w_holdout_force": float(args.w_holdout_force),
                "w_holdout_energy": float(args.w_holdout_energy),
                "w_gap_energy": float(args.w_gap_energy),
            },
            "scan": {
                "sites": list(site_uvs.keys()),
                "eval_z_points": int(args.eval_z_points),
                "train_low_count": int(args.train_low_count),
                "train_high_count": int(args.train_high_count),
                "stress_z_min": float(args.stress_z_min),
                "z_min": float(args.z_min),
                "z_max": float(args.z_max),
            },
            "trials": trials,
            "best": best_record,
        },
        out_dir / "calibration_results.json",
    )
    print(f"best objective={best_record['objective']:.6f}")
    print(best_record)


if __name__ == "__main__":
    main()
