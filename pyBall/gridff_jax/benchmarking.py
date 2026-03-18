"""Benchmark helpers for teacher sanity checks, ablations, and pose-path smoothness."""

from __future__ import annotations

import json
import time
from pathlib import Path

from ase.io import iread
import matplotlib.pyplot as plt
from mace.calculators import MACECalculator
import numpy as np

from .config import FeatureToggles
from .hybrid_energy import HybridGridFFModel, HybridParameters
from .interfaces import ModelEvaluation, PoseBatch, TeacherResult
from .utils import ensure_dir, fractional_to_cartesian, quaternion_to_matrix, save_json
from .validation import compute_error_metrics


ABLATION_CASES = {
    "plq": FeatureToggles(
        use_density_charge=True,
        use_locpot=True,
        use_ct_qeq=False,
        use_image_charge=False,
        use_reactive_grid=False,
        use_teacher_residual=True,
    ),
    "plq_reactive": FeatureToggles(
        use_density_charge=True,
        use_locpot=True,
        use_ct_qeq=False,
        use_image_charge=False,
        use_reactive_grid=True,
        use_teacher_residual=True,
    ),
    "plq_ct": FeatureToggles(
        use_density_charge=True,
        use_locpot=True,
        use_ct_qeq=True,
        use_image_charge=False,
        use_reactive_grid=False,
        use_teacher_residual=True,
    ),
    "plq_ct_image": FeatureToggles(
        use_density_charge=True,
        use_locpot=True,
        use_ct_qeq=True,
        use_image_charge=True,
        use_reactive_grid=False,
        use_teacher_residual=True,
    ),
    "full": FeatureToggles(),
}


ABLATION_CASES["full"] = FeatureToggles(
    use_density_charge=True,
    use_locpot=True,
    use_ct_qeq=True,
    use_image_charge=True,
    use_reactive_grid=True,
    use_teacher_residual=True,
)


def load_fit_params(path: str | Path):
    with Path(path).open("r", encoding="utf8") as fin:
        payload = json.load(fin)
    return HybridParameters(
        pauli=payload["pauli"],
        london=payload["london"],
        reactive=payload["reactive"],
        static_charge=payload["static_charge"],
        req_radius_offset=payload.get("req_radius_offset", {}),
        req_energy_scale=payload.get("req_energy_scale", {}),
        chi=payload["chi"],
        hardness=payload["hardness"],
        image_scale=payload["image_scale"],
        image_plane=payload["image_plane"],
        image_damping=payload["image_damping"],
        sample_shift_z=payload.get("sample_shift_z", 0.0),
        coulomb_shift_z=payload.get("coulomb_shift_z", 0.0),
        reservoir_chi=payload["reservoir_chi"],
        reservoir_hardness=payload["reservoir_hardness"],
        total_charge=payload["total_charge"],
    )


def count_extxyz_structures(extxyz_path: str | Path):
    count = 0
    for _ in iread(str(extxyz_path), index=":"):
        count += 1
    return count


def teacher_sanity_on_extxyz(
    model_path: str,
    extxyz_path: str,
    limit: int | None = None,
    device: str = "cpu",
    stride: int = 1,
    total_count: int | None = None,
):
    calc = MACECalculator(model_paths=model_path, device=device, default_dtype="float64")
    energy_errors = []
    force_maes = []
    n_selected = 0
    t0 = time.perf_counter()
    for index, atoms in enumerate(iread(str(extxyz_path), index=":")):
        if index % max(stride, 1) != 0:
            continue
        if limit is not None and n_selected >= limit:
            break
        ref_energy = atoms.info.get("REF_energy", atoms.info.get("energy"))
        ref_forces = atoms.arrays.get("REF_forces", atoms.arrays.get("forces"))
        test_atoms = atoms.copy()
        test_atoms.calc = calc
        pred_energy = float(test_atoms.get_potential_energy())
        pred_forces = np.asarray(test_atoms.get_forces(), dtype=float)
        energy_errors.append(pred_energy - ref_energy)
        force_maes.append(np.mean(np.linalg.norm(pred_forces - ref_forces, axis=1)))
        n_selected += 1
    elapsed = time.perf_counter() - t0
    energy_errors = np.asarray(energy_errors, dtype=float)
    force_maes = np.asarray(force_maes, dtype=float)
    if energy_errors.size == 0:
        energy_errors = np.zeros(1, dtype=float)
        force_maes = np.zeros(1, dtype=float)
    return {
        "source_path": str(extxyz_path),
        "n_source": int(total_count if total_count is not None else n_selected),
        "n_selected": int(n_selected),
        "stride": int(stride),
        "seconds_total": elapsed,
        "seconds_per_structure": elapsed / max(n_selected, 1),
        "energy": compute_error_metrics(np.zeros_like(energy_errors), energy_errors),
        "force_mae": {
            "mean": float(np.mean(force_maes)),
            "median": float(np.median(force_maes)),
            "max": float(np.max(force_maes)),
        },
    }


def run_teacher_sanity_suite(model_path: str, small_extxyz: str, full_extxyz: str, full_sample_size: int, device: str):
    small_total = count_extxyz_structures(small_extxyz)
    full_total = count_extxyz_structures(full_extxyz)
    full_stride = max(full_total // max(full_sample_size, 1), 1)
    return {
        "model_path": model_path,
        "device": device,
        "small_full": teacher_sanity_on_extxyz(
            model_path=model_path,
            extxyz_path=small_extxyz,
            limit=None,
            device=device,
            stride=1,
            total_count=small_total,
        ),
        "full_uniform_sample": teacher_sanity_on_extxyz(
            model_path=model_path,
            extxyz_path=full_extxyz,
            limit=full_sample_size,
            device=device,
            stride=full_stride,
            total_count=full_total,
        ),
    }


def evaluate_model_on_datasets(model, params, datasets):
    summary = {}
    total_energy = []
    total_force = []
    for poses, teacher in datasets:
        t0 = time.perf_counter()
        pred = model.evaluate_batch(poses, params=params, compute_forces=True)
        elapsed = time.perf_counter() - t0
        eerr = pred.energies - teacher.energies
        ferr = (pred.forces - teacher.forces).reshape(-1)
        summary[poses.adsorbate.name] = {
            "n_samples": int(len(poses.positions)),
            "energy": compute_error_metrics(np.zeros_like(eerr), eerr),
            "force": compute_error_metrics(np.zeros_like(ferr), ferr),
            "student_seconds_total": elapsed,
            "student_seconds_per_pose": elapsed / max(len(poses.positions), 1),
        }
        total_energy.append(eerr)
        total_force.append(ferr)
    energy_concat = np.concatenate(total_energy)
    force_concat = np.concatenate(total_force)
    summary["all_energy"] = compute_error_metrics(np.zeros_like(energy_concat), energy_concat)
    summary["all_force"] = compute_error_metrics(np.zeros_like(force_concat), force_concat)
    return summary


def evaluate_postfit_ablation_suite(density, datasets, config, params, prefer_jax: bool = True, case_names=None):
    case_names = case_names or list(ABLATION_CASES.keys())
    reactive_elements = tuple(sorted({symbol for poses, _ in datasets for symbol in poses.adsorbate.symbols}))
    summary = {}
    for case_name in case_names:
        model = HybridGridFFModel(
            density=density,
            reactive_elements=reactive_elements,
            toggles=ABLATION_CASES[case_name],
            grid_config=config.grid,
            hybrid_config=config.hybrid_model,
            prefer_jax=prefer_jax,
        )
        summary[case_name] = evaluate_model_on_datasets(model, params, datasets)
    full_energy = summary["full"]["all_energy"]["rmse"]
    full_force = summary["full"]["all_force"]["rmse"]
    for case_name in case_names:
        summary[case_name]["delta_vs_full"] = {
            "energy_rmse": summary[case_name]["all_energy"]["rmse"] - full_energy,
            "force_rmse": summary[case_name]["all_force"]["rmse"] - full_force,
        }
    return summary


def _rotation_quaternion(axis, angle):
    axis = np.asarray(axis, dtype=float)
    axis = axis / max(np.linalg.norm(axis), 1.0e-12)
    half = 0.5 * angle
    return np.array([np.cos(half), *(np.sin(half) * axis)], dtype=float)


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


def _transform_local(local, quaternion, anchor):
    rot = quaternion_to_matrix(quaternion)
    return local @ rot.T + anchor[None, :]


def _site_uv_lookup(surface_sites):
    lookup = {}
    for site in surface_sites:
        label = site["label"]
        if label not in lookup:
            lookup[label] = np.asarray(site["uv"], dtype=float)
    return lookup


def _anchor_from_uv_z(density, uv, z_ref, z_height):
    anchor = fractional_to_cartesian(np.array([uv[0], uv[1], 0.0], dtype=float), density.cell)
    anchor[2] = z_ref + z_height
    return anchor


def _compose_pose(local, density, uv, z_ref, z_height, quaternion):
    anchor = _anchor_from_uv_z(density, uv, z_ref, z_height)
    positions = _transform_local(local, quaternion, anchor)
    pose = np.concatenate([np.asarray(uv, dtype=float), np.array([z_height], dtype=float), np.asarray(quaternion, dtype=float)])
    return positions, pose


def _path_error_plots(step_axis, teacher_energy, student_energy, teacher_force, student_force, out_dir, name):
    step_axis = np.asarray(step_axis, dtype=float)
    teacher_energy = np.asarray(teacher_energy, dtype=float)
    student_energy = np.asarray(student_energy, dtype=float)
    teacher_force = np.asarray(teacher_force, dtype=float).reshape(len(step_axis), -1)
    student_force = np.asarray(student_force, dtype=float).reshape(len(step_axis), -1)
    energy_abs = np.abs(student_energy - teacher_energy)
    force_abs = np.linalg.norm(student_force - teacher_force, axis=1)
    teacher_force_norm = np.linalg.norm(teacher_force, axis=1)
    student_force_norm = np.linalg.norm(student_force, axis=1)
    markevery = max(len(step_axis) // 12, 1)

    plt.figure(figsize=(6, 4))
    plt.plot(step_axis, teacher_energy, label="MAD-SURF", color="tab:blue", lw=2.4, alpha=0.95, zorder=3)
    plt.plot(
        step_axis,
        student_energy,
        label="GridFF",
        color="tab:orange",
        lw=2.0,
        ls="--",
        marker="o",
        ms=2.6,
        markevery=markevery,
        alpha=0.9,
        zorder=4,
    )
    plt.xlabel("Path coordinate")
    plt.ylabel("Energy [eV]")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_dir / f"{name}_energy_path.png", dpi=160)
    plt.close()

    plt.figure(figsize=(6, 4))
    plt.plot(step_axis, teacher_force_norm, label="MAD-SURF", color="tab:blue", lw=2.4, alpha=0.95, zorder=3)
    plt.plot(
        step_axis,
        student_force_norm,
        label="GridFF",
        color="tab:orange",
        lw=2.0,
        ls="--",
        marker="o",
        ms=2.6,
        markevery=markevery,
        alpha=0.9,
        zorder=4,
    )
    plt.xlabel("Path coordinate")
    plt.ylabel("|Force| [eV/A]")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_dir / f"{name}_force_norm_path.png", dpi=160)
    plt.close()

    for logscale, suffix in ((False, "png"), (True, "log.png")):
        plt.figure(figsize=(6, 4))
        plt.plot(step_axis, energy_abs + 1.0e-16)
        if logscale:
            plt.yscale("log")
        plt.xlabel("Path coordinate")
        plt.ylabel("|Energy error| [eV]")
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(out_dir / f"{name}_energy_abs_error.{suffix}", dpi=160)
        plt.close()

        plt.figure(figsize=(6, 4))
        plt.plot(step_axis, force_abs + 1.0e-16)
        if logscale:
            plt.yscale("log")
        plt.xlabel("Path coordinate")
        plt.ylabel("|Force error| [eV/A]")
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(out_dir / f"{name}_force_abs_error.{suffix}", dpi=160)
        plt.close()


def _save_path_trace(step_axis, teacher_energy, student_energy, teacher_force, student_force, out_dir, name):
    np.savez_compressed(
        out_dir / f"{name}_trace.npz",
        step_axis=np.asarray(step_axis, dtype=float),
        teacher_energy=np.asarray(teacher_energy, dtype=float),
        student_energy=np.asarray(student_energy, dtype=float),
        teacher_force=np.asarray(teacher_force, dtype=float),
        student_force=np.asarray(student_force, dtype=float),
    )


def _plane_error_plots(step_x, step_y, teacher_energy, student_energy, out_dir, name, xlabel, ylabel):
    teacher_energy = np.asarray(teacher_energy, dtype=float)
    student_energy = np.asarray(student_energy, dtype=float)
    error = student_energy - teacher_energy
    abs_error = np.abs(error)
    extent = [np.min(step_x), np.max(step_x), np.min(step_y), np.max(step_y)]
    plots = [
        ("teacher", teacher_energy),
        ("student", student_energy),
        ("error", error),
        ("abs_error", abs_error),
        ("abs_error_log10", np.log10(abs_error + 1.0e-12)),
    ]
    for label, values in plots:
        plt.figure(figsize=(6, 5))
        plt.imshow(values, origin="lower", extent=extent, aspect="auto")
        plt.colorbar()
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.tight_layout()
        plt.savefig(out_dir / f"{name}_{label}.png", dpi=160)
        plt.close()


def _save_plane_trace(step_x, step_y, teacher_energy, student_energy, teacher_force, student_force, out_dir, name):
    np.savez_compressed(
        out_dir / f"{name}_trace.npz",
        step_x=np.asarray(step_x, dtype=float),
        step_y=np.asarray(step_y, dtype=float),
        teacher_energy=np.asarray(teacher_energy, dtype=float),
        student_energy=np.asarray(student_energy, dtype=float),
        teacher_force=np.asarray(teacher_force, dtype=float),
        student_force=np.asarray(student_force, dtype=float),
    )


def _slice_pose_batch(poses, start, stop):
    return PoseBatch(
        adsorbate=poses.adsorbate,
        positions=np.asarray(poses.positions[start:stop], dtype=float),
        pose_params=np.asarray(poses.pose_params[start:stop], dtype=float),
        site_labels=list(poses.site_labels[start:stop]),
        metadata=poses.metadata,
    )


def _concat_model_evaluations(evaluations):
    return ModelEvaluation(
        energies=np.concatenate([item.energies for item in evaluations], axis=0),
        forces=np.concatenate([item.forces for item in evaluations], axis=0),
        components={
            key: np.concatenate([item.components[key] for item in evaluations], axis=0)
            for key in evaluations[0].components.keys()
        }
        if evaluations[0].components
        else {},
        charges=np.concatenate([item.charges for item in evaluations], axis=0) if evaluations[0].charges is not None else None,
    )


def _trim_model_evaluation(evaluation, size):
    size = int(size)
    return ModelEvaluation(
        energies=np.asarray(evaluation.energies[:size], dtype=float),
        forces=np.asarray(evaluation.forces[:size], dtype=float),
        components={key: np.asarray(value[:size], dtype=float) for key, value in evaluation.components.items()},
        charges=np.asarray(evaluation.charges[:size], dtype=float) if evaluation.charges is not None else None,
    )


def _pad_pose_batch(poses, target_size):
    target_size = int(target_size)
    current_size = int(len(poses.positions))
    if current_size >= target_size:
        return poses, current_size
    pad_count = target_size - current_size
    positions_tail = np.repeat(np.asarray(poses.positions[-1:], dtype=float), pad_count, axis=0)
    pose_tail = np.repeat(np.asarray(poses.pose_params[-1:], dtype=float), pad_count, axis=0)
    site_tail = [poses.site_labels[-1]] * pad_count
    padded = PoseBatch(
        adsorbate=poses.adsorbate,
        positions=np.concatenate([np.asarray(poses.positions, dtype=float), positions_tail], axis=0),
        pose_params=np.concatenate([np.asarray(poses.pose_params, dtype=float), pose_tail], axis=0),
        site_labels=list(poses.site_labels) + site_tail,
        metadata=poses.metadata,
    )
    return padded, current_size


def _evaluate_teacher_chunked(teacher_backend, density, poses, chunk_size):
    if len(poses.positions) <= chunk_size:
        t0 = time.perf_counter()
        result = teacher_backend.evaluate_batch(density, poses)
        elapsed = time.perf_counter() - t0
        return result, elapsed
    results = []
    t0 = time.perf_counter()
    for start in range(0, len(poses.positions), chunk_size):
        stop = min(start + chunk_size, len(poses.positions))
        results.append(teacher_backend.evaluate_batch(density, _slice_pose_batch(poses, start, stop)))
    elapsed = time.perf_counter() - t0
    metadata = dict(results[-1].metadata) if results and results[-1].metadata else {}
    metadata["teacher_eval_seconds"] = elapsed
    metadata["seconds_per_pose"] = elapsed / max(len(poses.positions), 1)
    return TeacherResult(
        energies=np.concatenate([item.energies for item in results], axis=0),
        forces=np.concatenate([item.forces for item in results], axis=0),
        metadata=metadata,
    ), elapsed


def _evaluate_model_chunked(model, poses, params, compute_forces, chunk_size, warmup=False):
    if model.use_jax:
        warm_stop = min(chunk_size, len(poses.positions))
        warm_chunk = _slice_pose_batch(poses, 0, warm_stop)
        warm_chunk, _ = _pad_pose_batch(warm_chunk, chunk_size)
        if warmup:
            model.evaluate_batch(warm_chunk, params=params, compute_forces=compute_forces)
    elif len(poses.positions) <= chunk_size:
        t0 = time.perf_counter()
        result = model.evaluate_batch(poses, params=params, compute_forces=compute_forces)
        elapsed = time.perf_counter() - t0
        return result, elapsed
    results = []
    t0 = time.perf_counter()
    for start in range(0, len(poses.positions), chunk_size):
        stop = min(start + chunk_size, len(poses.positions))
        chunk = _slice_pose_batch(poses, start, stop)
        if model.use_jax:
            padded_chunk, real_size = _pad_pose_batch(chunk, chunk_size)
            result = model.evaluate_batch(padded_chunk, params=params, compute_forces=compute_forces)
            result = _trim_model_evaluation(result, real_size)
        else:
            result = model.evaluate_batch(chunk, params=params, compute_forces=compute_forces)
        results.append(result)
    elapsed = time.perf_counter() - t0
    return _concat_model_evaluations(results), elapsed


def benchmark_pose_paths(
    teacher_backend,
    model,
    params,
    density,
    adsorbate,
    pose_metadata,
    out_dir,
    teacher_chunk_size: int = 128,
    student_chunk_size: int = 64,
):
    out_dir = ensure_dir(out_dir)
    z_ref = pose_metadata["z_ref"]
    results = {}
    plane_results = {}

    def save_partial():
        save_json(
            {
                "paths": results,
                "planes": plane_results,
            },
            out_dir / "path_benchmark.json",
        )

    site_lookup = _site_uv_lookup(pose_metadata["surface_sites"])
    if "top" in site_lookup:
        base_uv = site_lookup["top"]
    else:
        first_key = sorted(site_lookup.keys())[0]
        base_uv = site_lookup[first_key]

    local = adsorbate.positions - adsorbate.positions[adsorbate.anchor_index]
    q_identity = np.array([1.0, 0.0, 0.0, 0.0], dtype=float)
    q_yaw_30 = _rotation_quaternion(np.array([0.0, 0.0, 1.0], dtype=float), np.deg2rad(30.0))
    q_tilt_x_25 = _rotation_quaternion(np.array([1.0, 0.0, 0.0], dtype=float), np.deg2rad(25.0))
    q_tilt_y_m20 = _rotation_quaternion(np.array([0.0, 1.0, 0.0], dtype=float), np.deg2rad(-20.0))
    q_skew = _quaternion_multiply(q_yaw_30, _quaternion_multiply(q_tilt_x_25, q_tilt_y_m20))

    def evaluate_path(name, positions_batch, pose_params, path_axis=None):
        poses = PoseBatch(
            adsorbate=adsorbate,
            positions=np.asarray(positions_batch, dtype=float),
            pose_params=np.asarray(pose_params, dtype=float),
            site_labels=[name] * len(positions_batch),
            metadata=pose_metadata,
        )
        teacher, teacher_elapsed = _evaluate_teacher_chunked(teacher_backend, density, poses, teacher_chunk_size)
        pred, student_elapsed = _evaluate_model_chunked(
            model,
            poses,
            params=params,
            compute_forces=True,
            chunk_size=student_chunk_size,
            warmup=True,
        )

        eerr = pred.energies - teacher.energies
        ferr = (pred.forces - teacher.forces).reshape(-1)
        teacher_energy = teacher.energies
        student_energy = pred.energies
        teacher_force_series = teacher.forces.reshape(len(positions_batch), -1)
        student_force_series = pred.forces.reshape(len(positions_batch), -1)
        teacher_roughness = float(np.mean(np.abs(np.diff(teacher_energy, n=2)))) if len(teacher_energy) > 2 else 0.0
        student_roughness = float(np.mean(np.abs(np.diff(student_energy, n=2)))) if len(student_energy) > 2 else 0.0
        teacher_force_delta = float(np.max(np.linalg.norm(np.diff(teacher_force_series, axis=0), axis=1))) if len(positions_batch) > 1 else 0.0
        student_force_delta = float(np.max(np.linalg.norm(np.diff(student_force_series, axis=0), axis=1))) if len(positions_batch) > 1 else 0.0
        if path_axis is None:
            path_axis = np.arange(len(positions_batch), dtype=float)
        path_axis = np.asarray(path_axis, dtype=float)
        _path_error_plots(path_axis, teacher_energy, student_energy, teacher.forces, pred.forces, out_dir, name)
        _save_path_trace(path_axis, teacher_energy, student_energy, teacher.forces, pred.forces, out_dir, name)

        results[name] = {
            "n_points": int(len(positions_batch)),
            "energy_error": compute_error_metrics(np.zeros_like(eerr), eerr),
            "force_error": compute_error_metrics(np.zeros_like(ferr), ferr),
            "teacher_seconds_total": teacher_elapsed,
            "teacher_seconds_per_pose": teacher_elapsed / max(len(positions_batch), 1),
            "student_seconds_total": student_elapsed,
            "student_seconds_per_pose": student_elapsed / max(len(positions_batch), 1),
            "speedup_teacher_over_student": teacher_elapsed / max(student_elapsed, 1.0e-12),
            "teacher_energy_roughness": teacher_roughness,
            "student_energy_roughness": student_roughness,
            "teacher_force_delta_max": teacher_force_delta,
            "student_force_delta_max": student_force_delta,
        }
        save_partial()

    def evaluate_parameterized_path(name, coordinates, builder, xlabel):
        positions_batch = []
        pose_params = []
        for value in coordinates:
            positions, pose = builder(float(value))
            positions_batch.append(positions)
            pose_params.append(pose)
        evaluate_path(name, positions_batch, pose_params, path_axis=np.asarray(coordinates, dtype=float))
        results[name]["path_coordinate_min"] = float(np.min(coordinates))
        results[name]["path_coordinate_max"] = float(np.max(coordinates))
        results[name]["path_label"] = xlabel

    def evaluate_plane(name, xs, ys, builder, xlabel, ylabel):
        positions_batch = []
        pose_params = []
        for yval in ys:
            for xval in xs:
                positions, pose = builder(float(xval), float(yval))
                positions_batch.append(positions)
                pose_params.append(pose)
        poses = PoseBatch(
            adsorbate=adsorbate,
            positions=np.asarray(positions_batch, dtype=float),
            pose_params=np.asarray(pose_params, dtype=float),
            site_labels=[name] * len(positions_batch),
            metadata=pose_metadata,
        )
        teacher, teacher_elapsed = _evaluate_teacher_chunked(teacher_backend, density, poses, teacher_chunk_size)
        pred, student_elapsed = _evaluate_model_chunked(
            model,
            poses,
            params=params,
            compute_forces=True,
            chunk_size=student_chunk_size,
            warmup=True,
        )
        eerr = pred.energies - teacher.energies
        ferr = (pred.forces - teacher.forces).reshape(-1)
        nx = len(xs)
        ny = len(ys)
        teacher_energy = teacher.energies.reshape(ny, nx)
        student_energy = pred.energies.reshape(ny, nx)
        _plane_error_plots(np.asarray(xs, dtype=float), np.asarray(ys, dtype=float), teacher_energy, student_energy, out_dir, name, xlabel, ylabel)
        _save_plane_trace(
            np.asarray(xs, dtype=float),
            np.asarray(ys, dtype=float),
            teacher_energy,
            student_energy,
            teacher.forces.reshape(ny, nx, -1),
            pred.forces.reshape(ny, nx, -1),
            out_dir,
            name,
        )
        plane_results[name] = {
            "shape": [int(ny), int(nx)],
            "energy_error": compute_error_metrics(np.zeros_like(eerr), eerr),
            "force_error": compute_error_metrics(np.zeros_like(ferr), ferr),
            "teacher_seconds_total": teacher_elapsed,
            "teacher_seconds_per_pose": teacher_elapsed / max(len(positions_batch), 1),
            "student_seconds_total": student_elapsed,
            "student_seconds_per_pose": student_elapsed / max(len(positions_batch), 1),
            "speedup_teacher_over_student": teacher_elapsed / max(student_elapsed, 1.0e-12),
            "x_min": float(np.min(xs)),
            "x_max": float(np.max(xs)),
            "y_min": float(np.min(ys)),
            "y_max": float(np.max(ys)),
            "xlabel": xlabel,
            "ylabel": ylabel,
        }
        save_partial()

    z_values = np.linspace(1.4, 5.5, 80)
    for label, uv in site_lookup.items():
        evaluate_parameterized_path(
            f"{label}_z_line",
            z_values,
            lambda z_height, uv=uv: _compose_pose(local, density, uv, z_ref, z_height, q_identity),
            "z [A]",
        )

    uv_span = np.linspace(-0.35, 0.35, 81)
    evaluate_parameterized_path(
        "top_x_line",
        uv_span,
        lambda du: _compose_pose(local, density, (base_uv + np.array([du, 0.0])) % 1.0, z_ref, 2.8, q_identity),
        "du",
    )
    evaluate_parameterized_path(
        "top_y_line",
        uv_span,
        lambda dv: _compose_pose(local, density, (base_uv + np.array([0.0, dv])) % 1.0, z_ref, 2.8, q_identity),
        "dv",
    )
    evaluate_parameterized_path(
        "top_xy_diag_line",
        uv_span,
        lambda t: _compose_pose(local, density, (base_uv + np.array([t, t])) % 1.0, z_ref, 2.8, q_identity),
        "d(u=v)",
    )
    evaluate_parameterized_path(
        "top_xy_skew_line",
        uv_span,
        lambda t: _compose_pose(local, density, (base_uv + np.array([t, 0.37 * t])) % 1.0, z_ref, 2.8, q_identity),
        "d(u,0.37u)",
    )

    t_values = np.linspace(0.0, 1.0, 81)
    evaluate_parameterized_path(
        "top_xyz_skew_line",
        t_values,
        lambda t: _compose_pose(local, density, (base_uv + np.array([0.28 * t, 0.17 * t])) % 1.0, z_ref, 1.6 + 3.4 * t, q_identity),
        "t",
    )

    if adsorbate.natoms > 1:
        yaw_angles = np.linspace(0.0, 2.0 * np.pi, 96, endpoint=False)
        evaluate_parameterized_path(
            "top_yaw_rotation",
            yaw_angles,
            lambda angle: _compose_pose(local, density, base_uv, z_ref, 2.8, _rotation_quaternion(np.array([0.0, 0.0, 1.0], dtype=float), angle)),
            "yaw [rad]",
        )
        tilt_angles = np.linspace(-np.pi / 3.0, np.pi / 3.0, 64)
        evaluate_parameterized_path(
            "top_tilt_x_rotation",
            tilt_angles,
            lambda angle: _compose_pose(local, density, base_uv, z_ref, 2.8, _rotation_quaternion(np.array([1.0, 0.0, 0.0], dtype=float), angle)),
            "tilt_x [rad]",
        )
        evaluate_parameterized_path(
            "top_tilt_y_rotation",
            tilt_angles,
            lambda angle: _compose_pose(local, density, base_uv, z_ref, 2.8, _rotation_quaternion(np.array([0.0, 1.0, 0.0], dtype=float), angle)),
            "tilt_y [rad]",
        )
        evaluate_parameterized_path(
            "top_skew_orientation_line",
            uv_span,
            lambda t: _compose_pose(local, density, (base_uv + np.array([t, 0.21 * t])) % 1.0, z_ref, 2.4 + 0.8 * t, q_skew),
            "skew t",
        )

    plane_axis = np.linspace(-0.35, 0.35, 41)
    evaluate_plane(
        "top_xy_plane_z2p2",
        plane_axis,
        plane_axis,
        lambda du, dv: _compose_pose(local, density, (base_uv + np.array([du, dv])) % 1.0, z_ref, 2.2, q_identity),
        "du",
        "dv",
    )
    evaluate_plane(
        "top_xy_plane_z3p0",
        plane_axis,
        plane_axis,
        lambda du, dv: _compose_pose(local, density, (base_uv + np.array([du, dv])) % 1.0, z_ref, 3.0, q_identity),
        "du",
        "dv",
    )
    evaluate_plane(
        "top_xz_plane",
        np.linspace(-0.35, 0.35, 41),
        np.linspace(1.4, 5.5, 41),
        lambda du, z_height: _compose_pose(local, density, (base_uv + np.array([du, 0.0])) % 1.0, z_ref, z_height, q_identity),
        "du",
        "z [A]",
    )

    summary = {
        "paths": results,
        "planes": plane_results,
    }
    save_partial()
    return summary


def plot_ablation_summary(summary, out_dir):
    out_dir = ensure_dir(out_dir)
    cases = list(summary.keys())
    energy_rmse = [summary[case]["all_energy"]["rmse"] for case in cases]
    force_rmse = [summary[case]["all_force"]["rmse"] for case in cases]

    plt.figure(figsize=(8, 4))
    x = np.arange(len(cases))
    plt.bar(x, energy_rmse, width=0.6)
    plt.xticks(x, cases, rotation=30, ha="right")
    plt.ylabel("Energy RMSE [eV]")
    plt.grid(True, axis="y", alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_dir / "ablation_energy_rmse.png", dpi=160)
    plt.close()

    for scale, suffix in (("linear", "png"), ("log", "log.png")):
        plt.figure(figsize=(8, 4))
        plt.bar(x, force_rmse, width=0.6)
        plt.xticks(x, cases, rotation=30, ha="right")
        plt.ylabel("Force RMSE [eV/A]")
        if scale == "log":
            plt.yscale("log")
        plt.grid(True, axis="y", alpha=0.3)
        plt.tight_layout()
        plt.savefig(out_dir / f"ablation_force_rmse.{suffix}", dpi=160)
        plt.close()
