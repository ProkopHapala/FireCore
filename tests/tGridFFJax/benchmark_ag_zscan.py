#!/usr/bin/env python3
"""Focused Ag(111) z-scan benchmark for one fixed adsorbate orientation.

Typical CLI usage
-----------------
# PLQ-only baseline with London Fermi damping (37.1 meV holdout):
  python benchmark_ag_zscan.py --london-damp-d0 3.70 --london-damp-width 0.35 --prefer-jax

# Two-stage PLQ + raw 1/r⁶ lstsq correction (16.9 meV holdout):
  python benchmark_ag_zscan.py --london-damp-d0 3.70 --london-damp-width 0.35 \\
      --two-stage-c6 --raw-r6 --prefer-jax

Two-stage workflow (--two-stage-c6):
  Stage 1: Fit REQ parameters with c6_coeff=0 (PLQ-only, gradient descent)
  Stage 2: Freeze REQ, fit C₆ coefficients + energy_offset via np.linalg.lstsq
  Stage 3 (optional, --stage3-refine): Gentle joint REQ+C₆ gradient refinement

The --raw-r6 flag replaces the TS-damped dispersion grid with raw -Σ 1/|r-r_j|⁶
(no Fermi damping, c6_pair=1.0).  This gives a basis function with different
lateral variation than exponential density decay, enabling site differentiation
that PLQ alone cannot achieve.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.append(str(ROOT))

from pyBall.gridff_jax import RunConfig, backend_summary, save_config
from pyBall.gridff_jax.density_backends import make_density_backend
from pyBall.gridff_jax.export import export_firecore_artifacts
from pyBall.gridff_jax.fit import fit_hybrid_parameters, fit_two_stage_c6
from pyBall.gridff_jax.hybrid_energy import HybridGridFFModel, HybridParameters
from pyBall.gridff_jax.interfaces import ModelEvaluation, PoseBatch, TeacherResult
from pyBall.gridff_jax.pose_sampling import get_adsorbate, infer_reference_planes, infer_surface_sites, transform_adsorbate
from pyBall.gridff_jax.teacher_backends import make_teacher_backend
from pyBall.gridff_jax.utils import ensure_dir, quaternion_to_matrix, save_json
from pyBall.gridff_jax.validation import compute_error_metrics, plot_convergence, plot_error_histogram, plot_parity


DEFAULT_CHGCAR = "/home/niel/git/ORR_HER_Ag_Colab/results/Ag_ORR_HER/slab_clean/final_scf_12x12x1/CHGCAR"
DEFAULT_LOCPOT = "/home/niel/git/ORR_HER_Ag_Colab/results/Ag_ORR_HER/slab_clean/workfunc_12x12x1/LOCPOT"
DEFAULT_MADSURF_MODEL = str((Path(__file__).resolve().parent / "mad-surf_data" / "models" / "full_dataset_config_weights" / "MACE_model.model"))
def parse_args():
    base = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description="Run a dense top/bridge/hollow Ag z-scan benchmark for one fixed CHON adsorbate")
    parser.add_argument("--chgcar", default=DEFAULT_CHGCAR, help="Path to Ag CHGCAR")
    parser.add_argument("--locpot", default=DEFAULT_LOCPOT, help="Path to Ag LOCPOT")
    parser.add_argument("--model-path", default=DEFAULT_MADSURF_MODEL, help="Path to MAD-SURF/MACE model")
    parser.add_argument("--device", default="cuda", help="Teacher device")
    parser.add_argument("--out-dir", default=str(base / "ag_zscan_chon"), help="Output directory")
    parser.add_argument("--xyz-path", default=None, help="Optional adsorbate xyz file; when omitted, use the built-in adsorbate")
    parser.add_argument("--adsorbate-name", default="CHONH2", help="Built-in adsorbate name or display name for --xyz-path")
    parser.add_argument("--anchor-index", type=int, default=0, help="Anchor atom index in the xyz file")
    parser.add_argument("--use-input-charges", action="store_true", help="Use charges from the xyz file instead of fitting free static charges")
    parser.add_argument("--sites", default="top,bridge,hollow", help="Comma-separated site labels to scan")
    parser.add_argument("--z-min", type=float, default=2.0, help="Primary fit/eval lower bound above the Ag reference plane")
    parser.add_argument("--stress-z-min", type=float, default=1.8, help="Stress-test lower bound used only for reporting")
    parser.add_argument("--z-max", type=float, default=5.60, help="Maximum z height above the Ag reference plane")
    parser.add_argument("--eval-z-points", type=int, default=181, help="Dense evaluation points per site")
    parser.add_argument("--train-low-count", type=int, default=18, help="Training points per site in the dense chemisorption window")
    parser.add_argument("--train-high-count", type=int, default=12, help="Training points per site in the outer window")
    parser.add_argument("--z-dense-cutoff", type=float, default=2.8, help="Upper edge of the dense chemisorption window")
    parser.add_argument("--tilt-x-deg", type=float, default=0.0, help="Fixed tilt around x for the whole scan")
    parser.add_argument("--tilt-y-deg", type=float, default=0.0, help="Fixed tilt around y for the whole scan")
    parser.add_argument("--yaw-deg", type=float, default=0.0, help="Fixed yaw around z for the whole scan")
    parser.add_argument("--heldout-tilt-x-deg", type=float, default=None, help="Optional held-out evaluation tilt around x")
    parser.add_argument("--heldout-tilt-y-deg", type=float, default=None, help="Optional held-out evaluation tilt around y")
    parser.add_argument("--heldout-yaw-deg", type=float, default=None, help="Optional held-out evaluation yaw around z")
    parser.add_argument("--fit-static-charge", dest="fit_static_charge", action="store_true", help="Fit molecule-side static charges in the strict PLQ stage")
    parser.add_argument("--no-fit-static-charge", dest="fit_static_charge", action="store_false", help="Disable fitting of molecule-side static charges")
    parser.set_defaults(fit_static_charge=False)
    parser.add_argument("--max-steps", type=int, default=600, help="Optimizer steps")
    parser.add_argument("--force-weight", type=float, default=10.0, help="Force loss weight")
    parser.add_argument("--learning-rate", type=float, default=1.0e-2, help="Optimizer learning rate")
    parser.add_argument("--alpha-morse", type=float, default=None, help="Optional override of substrate Morse alpha")
    parser.add_argument("--ag-radius-scale", type=float, default=None, help="Optional override of Ag vdW radius scale")
    parser.add_argument("--ag-energy-scale", type=float, default=None, help="Optional override of Ag vdW energy scale")
    parser.add_argument("--pauli-power", type=float, default=None, help="metal_density_pauli_power for grid building")
    parser.add_argument("--london-power", type=float, default=None, help="metal_density_london_power for grid building")
    parser.add_argument("--direct-plq", action="store_true", help="Use direct pauli/london scalars instead of REQ Morse coupling")
    parser.add_argument("--fit-z-shift", action="store_true", help="Fit sample_shift_z (global z offset for grid sampling)")
    parser.add_argument("--london-damp-d0", type=float, default=0.0, help="London damping d0 (Fermi midpoint above surface, Angstrom)")
    parser.add_argument("--london-damp-width", type=float, default=0.5, help="London damping Fermi width (Angstrom)")
    parser.add_argument("--use-c6", action="store_true", help="Enable pairwise C6/R6 dispersion grid channel (vdW-surf screened)")
    parser.add_argument("--fit-c6", action="store_true", help="Fit per-element c6_coeff scaling factors (requires --use-c6)")
    parser.add_argument("--c6-rcut", type=float, default=15.0, help="Cutoff for pairwise C6 sum (Angstrom)")
    parser.add_argument("--two-stage-c6", action="store_true", help="Enable two-stage REQ→C6 optimization (Stage1: PLQ-only, Stage2: C6-only)")
    parser.add_argument("--stage2-max-steps", type=int, default=400, help="Steps for C6-only stage (default 400)")
    parser.add_argument("--stage2-lr", type=float, default=5.0e-3, help="Learning rate for C6-only stage (default 5e-3)")
    parser.add_argument("--stage2-force-weight", type=float, default=10.0, help="Force loss weight for C6-only stage (default 10.0)")
    parser.add_argument("--raw-r6", action="store_true", help="Use raw 1/r^6 basis (no TS damping) for Stage 2 lstsq — matches proof-of-concept")
    parser.add_argument("--stage3-refine", action="store_true", help="Enable optional joint REQ+C6 refinement stage")
    parser.add_argument("--stage3-max-steps", type=int, default=200, help="Steps for joint refinement (default 200)")
    parser.add_argument("--stage3-lr", type=float, default=1.0e-4, help="Learning rate for joint refinement (default 1e-4)")
    # DISABLED(2026-03): density-derived channels confirmed ineffective in ablation scan.
    parser.add_argument("--use-density-lap", action="store_true", help="[DISABLED] Density Laplacian — confirmed ineffective")
    parser.add_argument("--use-density-grad", action="store_true", help="[DISABLED] Density gradient — confirmed ineffective")
    parser.add_argument("--use-london-alt", action="store_true", help="[DISABLED] Alt London power — confirmed ineffective")
    parser.add_argument("--london-alt-power", type=float, default=1.5, help="[DISABLED] Power for alt London channel")
    parser.add_argument("--builder-mode", type=str, default="metal_density_plq", help="Grid builder mode (metal_density_plq, metal_dft_plq, parity_core)")
    parser.add_argument("--teacher-chunk-size", type=int, default=64, help="Teacher batch chunk size")
    parser.add_argument("--student-chunk-size", type=int, default=64, help="Student batch chunk size")
    parser.add_argument("--teacher-tile", type=str, default="1,1",
                        help="Slab tiling for teacher evaluation: 'NX,NY' (e.g. '2,2') or 'auto' to compute from molecule extent")
    parser.add_argument("--prefer-jax", dest="prefer_jax", action="store_true", help="Use the JAX student backend when available")
    parser.add_argument("--no-prefer-jax", dest="prefer_jax", action="store_false", help="Force the NumPy fallback backend")
    parser.set_defaults(prefer_jax=True)
    return parser.parse_args()


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
        t0 = time.perf_counter()
        result = model.evaluate_batch(poses, params=params, compute_forces=True)
        elapsed = time.perf_counter() - t0
        return result, elapsed
    energies = []
    forces = []
    components = {}
    charges = []
    t0 = time.perf_counter()
    for start in range(0, len(poses.positions), int(chunk_size)):
        stop = min(start + int(chunk_size), len(poses.positions))
        partial = model.evaluate_batch(_slice_pose_batch(poses, np.arange(start, stop, dtype=int)), params=params, compute_forces=True)
        energies.append(np.asarray(partial.energies, dtype=float))
        forces.append(np.asarray(partial.forces, dtype=float))
        for key, value in partial.components.items():
            components.setdefault(key, []).append(np.asarray(value, dtype=float))
        if partial.charges is not None:
            charges.append(np.asarray(partial.charges, dtype=float))
    elapsed = time.perf_counter() - t0
    return (
        ModelEvaluation(
            energies=np.concatenate(energies, axis=0),
            forces=np.concatenate(forces, axis=0),
            components={key: np.concatenate(value, axis=0) for key, value in components.items()},
            charges=np.concatenate(charges, axis=0) if charges else None,
        ),
        elapsed,
    )


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
    }
    if params.c6_coeff:
        payload["c6_coeff"] = dict(params.c6_coeff)
    if params.energy_offset != 0.0:
        payload["energy_offset"] = float(params.energy_offset)
    return payload


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


def _safe_corr(x, y):
    x = np.asarray(x, dtype=float).reshape(-1)
    y = np.asarray(y, dtype=float).reshape(-1)
    if x.size < 2 or y.size < 2:
        return 0.0
    if np.std(x) < 1.0e-12 or np.std(y) < 1.0e-12:
        return 0.0
    return float(np.corrcoef(x, y)[0, 1])


def _roughness(values):
    values = np.asarray(values, dtype=float).reshape(-1)
    if values.size <= 2:
        return 0.0
    return float(np.mean(np.abs(np.diff(values, n=2))))


def _force_delta_max(forces):
    forces = np.asarray(forces, dtype=float)
    if len(forces) <= 1:
        return 0.0
    return float(np.max(np.linalg.norm(np.diff(forces.reshape(len(forces), -1), axis=0), axis=1)))


def _set_strict_plq_config(config: RunConfig, args):
    config.density_backend.chgcar_path = args.chgcar
    config.density_backend.locpot_path = args.locpot
    config.teacher_backend.kind = "madsurf"
    config.teacher_backend.model_path = args.model_path
    config.teacher_backend.device = args.device
    # Parse teacher tiling: "auto" or "NX,NY"
    tile_str = getattr(args, "teacher_tile", "1,1")
    if tile_str.strip().lower() == "auto":
        config.teacher_backend.teacher_tile = (0, 0)  # 0,0 triggers auto in madsurf
    else:
        parts = tile_str.split(",")
        config.teacher_backend.teacher_tile = (int(parts[0]), int(parts[1]))
    config.grid.builder_mode = str(getattr(args, 'builder_mode', 'metal_density_plq'))
    config.grid.interpolation_order = 3
    if args.alpha_morse is not None:
        config.grid.alpha_morse = float(args.alpha_morse)
    if args.pauli_power is not None:
        config.grid.metal_density_pauli_power = float(args.pauli_power)
    if args.london_power is not None:
        config.grid.metal_density_london_power = float(args.london_power)
    if args.london_damp_d0 > 0.0:
        config.grid.london_damping_d0 = float(args.london_damp_d0)
        config.grid.london_damping_width = float(args.london_damp_width)
    if args.use_c6 or args.two_stage_c6:
        config.grid.use_pairwise_c6 = True
        config.grid.c6_rcut = float(args.c6_rcut)
        config.training.fit_c6_coeff = bool(args.fit_c6)
    if args.two_stage_c6:
        config.training.two_stage_c6 = True
        config.training.stage2_max_steps = int(args.stage2_max_steps)
        config.training.stage2_learning_rate = float(args.stage2_lr)
        config.training.stage2_force_weight = float(args.stage2_force_weight)
        config.training.stage3_refine = bool(args.stage3_refine)
        config.training.stage3_max_steps = int(args.stage3_max_steps)
        config.training.stage3_learning_rate = float(args.stage3_lr)
    # DISABLED(2026-03): density-derived channels confirmed ineffective in ablation scan.
    # if args.use_density_lap:
    #     config.grid.use_density_laplacian = True
    #     config.training.fit_density_lap = True
    # if args.use_density_grad:
    #     config.grid.use_density_gradient = True
    #     config.training.fit_density_grad = True
    # if args.use_london_alt:
    #     config.grid.use_london_alt = True
    #     config.grid.london_alt_power = float(args.london_alt_power)
    #     config.training.fit_london_alt = True
    if args.ag_radius_scale is not None:
        config.grid.req_scale_radius["Ag"] = float(args.ag_radius_scale)
    if args.ag_energy_scale is not None:
        config.grid.req_scale_energy["Ag"] = float(args.ag_energy_scale)
    config.toggles.use_ct_qeq = False
    config.toggles.use_image_charge = False
    config.toggles.use_reactive_grid = False
    config.hybrid_model.use_qeq = False
    config.hybrid_model.use_image = False
    config.hybrid_model.use_reactive_grid = False
    if args.direct_plq:
        config.hybrid_model.use_req_plq = False
        config.training.fit_req_radius_offset = False
        config.training.fit_req_energy_scale = False
    else:
        config.hybrid_model.use_req_plq = True
        config.training.fit_req_radius_offset = True
        config.training.fit_req_energy_scale = True
    config.training.fit_chi = False
    config.training.fit_hardness = False
    config.training.fit_image_plane = False
    config.training.fit_reactive = False
    config.training.fit_static_charge = bool(args.fit_static_charge)
    config.training.fit_sample_shift_z = bool(getattr(args, 'fit_z_shift', False))
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


def _metrics_for_mask(teacher_energy, student_energy, teacher_force, student_force, mask):
    mask = np.asarray(mask, dtype=bool)
    n_samples = int(np.count_nonzero(mask))
    if n_samples == 0:
        return {
            "n_samples": 0,
            "energy": {"rmse": 0.0, "mae": 0.0, "max_abs": 0.0, "mean_signed": 0.0},
            "force": {"rmse": 0.0, "mae": 0.0, "max_abs": 0.0, "mean_signed": 0.0},
        }
    e_true = np.asarray(teacher_energy, dtype=float)[mask]
    e_pred = np.asarray(student_energy, dtype=float)[mask]
    f_true = np.asarray(teacher_force, dtype=float)[mask]
    f_pred = np.asarray(student_force, dtype=float)[mask]
    return {
        "n_samples": n_samples,
        "energy": compute_error_metrics(e_true, e_pred),
        "force": compute_error_metrics(f_true.reshape(-1), f_pred.reshape(-1)),
    }


def _write_stage(out_dir, stage: str, payload=None):
    info = {"stage": stage}
    if payload is not None:
        info.update(payload)
    save_json(info, Path(out_dir) / "progress.json")


def _plot_site_traces(out_dir, site_label, z_values, teacher_energy, student_energy, teacher_force, student_force, components, train_mask, primary_z_min: float, stress_z_min: float):
    z_values = np.asarray(z_values, dtype=float)
    teacher_energy = np.asarray(teacher_energy, dtype=float)
    student_energy = np.asarray(student_energy, dtype=float)
    teacher_force_norm = np.linalg.norm(np.asarray(teacher_force, dtype=float).reshape(len(z_values), -1), axis=1)
    student_force_norm = np.linalg.norm(np.asarray(student_force, dtype=float).reshape(len(z_values), -1), axis=1)
    train_mask = np.asarray(train_mask, dtype=bool)
    energy_error = np.abs(student_energy - teacher_energy)
    force_error = np.linalg.norm(
        np.asarray(student_force, dtype=float).reshape(len(z_values), -1)
        - np.asarray(teacher_force, dtype=float).reshape(len(z_values), -1),
        axis=1,
    )
    force_norm_residual = student_force_norm - teacher_force_norm

    plt.figure(figsize=(7, 4.8))
    plt.plot(z_values, teacher_energy, label="MAD-SURF", color="tab:blue", lw=2.5)
    plt.plot(z_values, student_energy, label="GridFF", color="tab:orange", lw=2.0, ls="--")
    plt.xlabel("z [A]")
    plt.ylabel("Energy [eV]")
    plt.title(f"{site_label} z scan")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_dir / f"{site_label}_energy_trace.png", dpi=170)
    plt.close()

    plt.figure(figsize=(7, 4.8))
    plt.plot(z_values, teacher_force_norm, label="MAD-SURF", color="tab:blue", lw=2.5)
    plt.plot(z_values, student_force_norm, label="GridFF", color="tab:orange", lw=2.0, ls="--")
    plt.xlabel("z [A]")
    plt.ylabel("|Force| [eV/A]")
    plt.title(f"{site_label} force norm")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_dir / f"{site_label}_force_trace.png", dpi=170)
    plt.close()

    plt.figure(figsize=(7, 5.4))
    plt.plot(z_values, components["pauli"], label="Pauli", color="tab:red", lw=2.0)
    plt.plot(z_values, components["london"], label="London", color="tab:green", lw=2.0)
    plt.plot(z_values, components["coulomb"], label="Coulomb", color="tab:purple", lw=2.0)
    plt.xlabel("z [A]")
    plt.ylabel("Predicted component [eV]")
    plt.title(f"{site_label} GridFF components")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_dir / f"{site_label}_components.png", dpi=170)
    plt.close()

    plt.figure(figsize=(7, 4.8))
    plt.axhline(0.0, color="black", lw=1.0, alpha=0.5)
    plt.axvspan(float(stress_z_min), float(primary_z_min), color="tab:red", alpha=0.08)
    plt.plot(z_values, student_energy - teacher_energy, color="tab:red", lw=2.0)
    plt.xlabel("z [A]")
    plt.ylabel("E(GridFF) - E(MAD-SURF) [eV]")
    plt.title(f"{site_label} signed energy residual")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_dir / f"{site_label}_energy_residual_signed.png", dpi=170)
    plt.close()

    plt.figure(figsize=(7, 4.8))
    plt.axhline(0.0, color="black", lw=1.0, alpha=0.5)
    plt.axvspan(float(stress_z_min), float(primary_z_min), color="tab:red", alpha=0.08)
    plt.plot(z_values, force_norm_residual, color="tab:orange", lw=2.0)
    plt.xlabel("z [A]")
    plt.ylabel("|F|_GridFF - |F|_MAD-SURF [eV/A]")
    plt.title(f"{site_label} signed force-norm residual")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_dir / f"{site_label}_force_residual_signed.png", dpi=170)
    plt.close()

    for log_scale, suffix in ((False, "png"), (True, "log.png")):
        plt.figure(figsize=(7, 4.8))
        plt.plot(z_values, energy_error + 1.0e-16, color="tab:red", lw=2.0)
        if log_scale:
            plt.yscale("log")
        plt.xlabel("z [A]")
        plt.ylabel("|Energy error| [eV]")
        plt.title(f"{site_label} energy error")
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(out_dir / f"{site_label}_energy_error.{suffix}", dpi=170)
        plt.close()

        plt.figure(figsize=(7, 4.8))
        plt.plot(z_values, force_error + 1.0e-16, color="tab:orange", lw=2.0)
        if log_scale:
            plt.yscale("log")
        plt.xlabel("z [A]")
        plt.ylabel("|Force error| [eV/A]")
        plt.title(f"{site_label} force error")
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(out_dir / f"{site_label}_force_error.{suffix}", dpi=170)
        plt.close()


def main():
    args = parse_args()
    out_dir = ensure_dir(args.out_dir)
    plots_dir = ensure_dir(out_dir / "plots")
    traces_dir = ensure_dir(out_dir / "traces")

    config = RunConfig()
    _set_strict_plq_config(config, args)
    save_config(config, out_dir / "config.json")
    _write_stage(out_dir, "config_saved")

    density = make_density_backend(config.density_backend, surface=config.surface, grid=config.grid).load()
    teacher_backend = make_teacher_backend(config.teacher_backend)
    adsorbate = get_adsorbate(
        name=args.adsorbate_name,
        xyz_path=args.xyz_path,
        anchor_index=args.anchor_index,
        use_input_charges=args.use_input_charges,
    )
    site_uvs = _site_uv_lookup(density, args.sites)
    z_values = np.linspace(float(args.stress_z_min), float(args.z_max), int(args.eval_z_points), dtype=float)
    quaternion = _compose_fixed_quaternion(args.yaw_deg, args.tilt_x_deg, args.tilt_y_deg)
    heldout_requested = any(value is not None for value in (args.heldout_yaw_deg, args.heldout_tilt_x_deg, args.heldout_tilt_y_deg))
    heldout_quaternion = None
    if heldout_requested:
        heldout_quaternion = _compose_fixed_quaternion(
            args.yaw_deg if args.heldout_yaw_deg is None else args.heldout_yaw_deg,
            args.tilt_x_deg if args.heldout_tilt_x_deg is None else args.heldout_tilt_x_deg,
            args.tilt_y_deg if args.heldout_tilt_y_deg is None else args.heldout_tilt_y_deg,
        )
    dense_poses = _build_zscan_pose_batch(density, adsorbate, site_uvs, z_values, quaternion)
    _write_stage(
        out_dir,
        "poses_ready",
        {
            "n_total": int(len(dense_poses.positions)),
            "sites": list(site_uvs.keys()),
            "eval_z_points": int(len(z_values)),
            "primary_z_min": float(args.z_min),
            "stress_z_min": float(args.stress_z_min),
        },
    )
    print(f"[zscan] poses ready: n={len(dense_poses.positions)}", flush=True)

    teacher_dense, teacher_elapsed = _evaluate_teacher_chunked(
        teacher_backend=teacher_backend,
        density=density,
        poses=dense_poses,
        chunk_size=args.teacher_chunk_size,
    )
    _write_stage(
        out_dir,
        "teacher_done",
        {
            "teacher_seconds_total": float(teacher_elapsed),
            "teacher_seconds_per_pose": float(teacher_elapsed / max(len(dense_poses.positions), 1)),
        },
    )
    print(f"[zscan] teacher done in {teacher_elapsed:.3f}s", flush=True)

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
    holdout_mask = np.ones(len(dense_poses.positions), dtype=bool)
    holdout_mask[train_indices] = False
    train_poses = _slice_pose_batch(dense_poses, train_indices)
    train_teacher = _slice_teacher_result(teacher_dense, train_indices)

    reactive_elements = tuple(sorted(set(adsorbate.symbols)))
    model = HybridGridFFModel(
        density=density,
        reactive_elements=reactive_elements,
        toggles=config.toggles,
        grid_config=config.grid,
        hybrid_config=config.hybrid_model,
        prefer_jax=args.prefer_jax,
    )
    if getattr(args, "raw_r6", False):
        print("[zscan] Installing raw 1/r^6 dispersion grid (no TS damping)...", flush=True)
        model.install_raw_r6_grid()

    fit_t0 = time.perf_counter()
    if config.training.two_stage_c6:
        # Build primary-window dataset for lstsq Stage 2 (all z >= z_min points).
        # REQ Stage 1 fits on training data; lstsq Stage 2 uses all primary-window
        # data because it's a 5-parameter linear fit with negligible overfitting risk.
        primary_mask = np.zeros(len(dense_poses.positions), dtype=bool)
        n_z = len(z_values)
        for si in range(len(site_uvs)):
            lo, hi = si * n_z, (si + 1) * n_z
            primary_mask[lo:hi] = z_values >= float(args.z_min)
        primary_indices = np.where(primary_mask)[0]
        primary_poses = _slice_pose_batch(dense_poses, primary_indices)
        primary_teacher = _slice_teacher_result(teacher_dense, primary_indices)
        fit_result = fit_two_stage_c6(
            density=density,
            datasets=[(train_poses, train_teacher)],
            model=model,
            training=config.training,
            initial_params=None,
            lstsq_datasets=[(primary_poses, primary_teacher)],
            use_raw_r6=getattr(args, "raw_r6", False),
        )
    else:
        fit_result = fit_hybrid_parameters(
            density=density,
            datasets=[(train_poses, train_teacher)],
            model=model,
            training=config.training,
            initial_params=None,
        )
    fit_elapsed = time.perf_counter() - fit_t0
    save_json(_fit_param_payload(fit_result.params), out_dir / "fit_params.json")
    save_json(fit_result.metrics.get("constraint_report", {}), out_dir / "fit_constraints.json")
    _write_stage(
        out_dir,
        "fit_done",
        {
            "fit_seconds_total": float(fit_elapsed),
            "history_points": int(len(fit_result.history)),
        },
    )
    print(f"[zscan] fit done in {fit_elapsed:.3f}s", flush=True)

    # Export FireCore-compatible grid artifacts (Bspline_PLQd.npy, etc.)
    export_paths = export_firecore_artifacts(
        out_dir=out_dir,
        density=density,
        model=model,
        fit_result=fit_result,
        toggles=config.toggles,
        teacher_backend_id=config.teacher_backend.kind,
    )
    print(f"[zscan] grid exported: {export_paths['plq_path']}", flush=True)

    student_dense, student_elapsed = _evaluate_student_chunked(
        model=model,
        poses=dense_poses,
        params=fit_result.params,
        chunk_size=args.student_chunk_size,
    )
    _write_stage(
        out_dir,
        "student_done",
        {
            "student_seconds_total": float(student_elapsed),
            "student_seconds_per_pose": float(student_elapsed / max(len(dense_poses.positions), 1)),
        },
    )
    print(f"[zscan] student done in {student_elapsed:.3f}s", flush=True)

    teacher_energy = np.asarray(teacher_dense.energies, dtype=float)
    student_energy = np.asarray(student_dense.energies, dtype=float)
    teacher_force = np.asarray(teacher_dense.forces, dtype=float)
    student_force = np.asarray(student_dense.forces, dtype=float)
    train_mask = np.zeros(len(dense_poses.positions), dtype=bool)
    train_mask[train_indices] = True

    summary = {
        "adsorbate": {
            "name": adsorbate.name,
            "natoms": int(adsorbate.natoms),
            "symbols": list(adsorbate.symbols),
            "source_file": str(args.xyz_path) if args.xyz_path else None,
            "used_input_charges": bool(args.use_input_charges),
        },
        "orientation": {
            "yaw_deg": float(args.yaw_deg),
            "tilt_x_deg": float(args.tilt_x_deg),
            "tilt_y_deg": float(args.tilt_y_deg),
            "quaternion": quaternion.tolist(),
            "rotation_matrix": quaternion_to_matrix(quaternion).tolist(),
        },
        "sampling": {
            "sites": list(site_uvs.keys()),
            "eval_z_points": int(len(z_values)),
            "train_low_count_per_site": int(args.train_low_count),
            "train_high_count_per_site": int(args.train_high_count),
            "train_z_points_per_site": int(len(np.unique(train_indices % len(z_values)))),
            "primary_z_min": float(args.z_min),
            "stress_z_min": float(args.stress_z_min),
            "z_dense_cutoff": float(args.z_dense_cutoff),
            "z_min": float(z_values.min()),
            "z_max": float(z_values.max()),
            "n_total": int(len(dense_poses.positions)),
            "n_train": int(np.count_nonzero(train_mask)),
            "n_holdout": int(np.count_nonzero(holdout_mask)),
        },
        "timing": {
            "teacher_seconds_total": float(teacher_elapsed),
            "teacher_seconds_per_pose": float(teacher_elapsed / max(len(dense_poses.positions), 1)),
            "fit_seconds_total": float(fit_elapsed),
            "student_seconds_total": float(student_elapsed),
            "student_seconds_per_pose": float(student_elapsed / max(len(dense_poses.positions), 1)),
            "speedup_teacher_over_student_eval": float(teacher_elapsed / max(student_elapsed, 1.0e-12)),
        },
        "backend": {
            "teacher_device": args.device,
            "student_backend": backend_summary(),
            "student_use_jax": bool(model.use_jax),
        },
        "fit": {
            "optimizer_metrics": fit_result.metrics,
            "optimizer_result": fit_result.optimizer_result,
        },
        "splits": {
            "train": _metrics_for_mask(teacher_energy, student_energy, teacher_force, student_force, train_mask),
            "holdout": _metrics_for_mask(teacher_energy, student_energy, teacher_force, student_force, holdout_mask),
            "all": _metrics_for_mask(teacher_energy, student_energy, teacher_force, student_force, np.ones_like(train_mask, dtype=bool)),
        },
        "sites": {},
    }
    if heldout_requested:
        summary["heldout_orientation"] = {
            "yaw_deg": float(args.yaw_deg if args.heldout_yaw_deg is None else args.heldout_yaw_deg),
            "tilt_x_deg": float(args.tilt_x_deg if args.heldout_tilt_x_deg is None else args.heldout_tilt_x_deg),
            "tilt_y_deg": float(args.tilt_y_deg if args.heldout_tilt_y_deg is None else args.heldout_tilt_y_deg),
            "quaternion": heldout_quaternion.tolist(),
            "rotation_matrix": quaternion_to_matrix(heldout_quaternion).tolist(),
        }
    z_vector = dense_poses.pose_params[:, 2]
    window_masks = _window_masks(z_vector, primary_z_min=float(args.z_min), stress_z_min=float(args.stress_z_min), z_max=float(args.z_max))
    summary["windows"] = {
        "primary": {
            "all": _metrics_for_mask(teacher_energy, student_energy, teacher_force, student_force, window_masks["primary"]),
            "train": _metrics_for_mask(teacher_energy, student_energy, teacher_force, student_force, window_masks["primary"] & train_mask),
            "holdout": _metrics_for_mask(teacher_energy, student_energy, teacher_force, student_force, window_masks["primary"] & holdout_mask),
        },
        "stress": {
            "all": _metrics_for_mask(teacher_energy, student_energy, teacher_force, student_force, window_masks["stress"]),
        },
        "full": {
            "all": _metrics_for_mask(teacher_energy, student_energy, teacher_force, student_force, window_masks["full"]),
        },
    }
    primary_train_rmse = summary["windows"]["primary"]["train"]["energy"]["rmse"]
    primary_holdout_rmse = summary["windows"]["primary"]["holdout"]["energy"]["rmse"]
    summary["acceptance"] = {
        "primary_energy_rmse": float(summary["windows"]["primary"]["all"]["energy"]["rmse"]),
        "primary_force_rmse": float(summary["windows"]["primary"]["all"]["force"]["rmse"]),
        "primary_train_holdout_gap": float(abs(primary_train_rmse - primary_holdout_rmse)),
        "stress_energy_rmse": float(summary["windows"]["stress"]["all"]["energy"]["rmse"]),
        "stress_force_rmse": float(summary["windows"]["stress"]["all"]["force"]["rmse"]),
        "constraint_limited": bool(fit_result.metrics.get("constraint_report", {}).get("constraint_limited", False)),
    }

    residual = student_energy - teacher_energy
    plot_parity(teacher_energy, student_energy, plots_dir, "energy_parity.png")
    plot_parity(teacher_force.reshape(-1), student_force.reshape(-1), plots_dir, "force_parity.png")
    plot_error_histogram(student_energy - teacher_energy, plots_dir, "energy_abs_hist.png", log_scale=False)
    plot_error_histogram(student_energy - teacher_energy, plots_dir, "energy_abs_hist_log.png", log_scale=True)
    plot_error_histogram(student_force - teacher_force, plots_dir, "force_abs_hist.png", log_scale=False)
    plot_error_histogram(student_force - teacher_force, plots_dir, "force_abs_hist_log.png", log_scale=True)
    plot_convergence(fit_result.history, plots_dir, name_prefix="training_loss")

    for site_index, site_label in enumerate(site_uvs.keys()):
        lo = site_index * len(z_values)
        hi = (site_index + 1) * len(z_values)
        sl = slice(lo, hi)
        site_train_mask = train_mask[sl]
        site_residual = residual[sl]
        site_force_error = np.linalg.norm((student_force[sl] - teacher_force[sl]).reshape(len(z_values), -1), axis=1)
        site_components = {
            "pauli": np.asarray(student_dense.components["pauli"][sl], dtype=float),
            "london": np.asarray(student_dense.components["london"][sl], dtype=float),
            "coulomb": np.asarray(student_dense.components["coulomb"][sl], dtype=float),
        }
        site_summary = {
            "metrics_train": _metrics_for_mask(
                teacher_energy[sl],
                student_energy[sl],
                teacher_force[sl],
                student_force[sl],
                site_train_mask,
            ),
            "metrics_holdout": _metrics_for_mask(
                teacher_energy[sl],
                student_energy[sl],
                teacher_force[sl],
                student_force[sl],
                ~site_train_mask,
            ),
            "metrics_all": _metrics_for_mask(
                teacher_energy[sl],
                student_energy[sl],
                teacher_force[sl],
                student_force[sl],
                np.ones(len(z_values), dtype=bool),
            ),
            "energy_residual_mean_signed": float(np.mean(site_residual)),
            "energy_residual_max_abs": float(np.max(np.abs(site_residual))),
            "force_error_norm_mean": float(np.mean(site_force_error)),
            "force_error_norm_max": float(np.max(site_force_error)),
            "teacher_energy_roughness": _roughness(teacher_energy[sl]),
            "student_energy_roughness": _roughness(student_energy[sl]),
            "teacher_force_delta_max": _force_delta_max(teacher_force[sl]),
            "student_force_delta_max": _force_delta_max(student_force[sl]),
            "corr_energy_residual_z": _safe_corr(site_residual, z_values),
            "corr_force_error_z": _safe_corr(site_force_error, z_values),
            "corr_energy_residual_pauli": _safe_corr(site_residual, site_components["pauli"]),
            "corr_energy_residual_london": _safe_corr(site_residual, site_components["london"]),
            "corr_energy_residual_coulomb": _safe_corr(site_residual, site_components["coulomb"]),
            "corr_force_error_pauli": _safe_corr(site_force_error, site_components["pauli"]),
            "corr_force_error_london": _safe_corr(site_force_error, site_components["london"]),
            "corr_force_error_coulomb": _safe_corr(site_force_error, site_components["coulomb"]),
            "worst_energy_index": int(np.argmax(np.abs(site_residual))),
            "worst_force_index": int(np.argmax(site_force_error)),
        }
        worst_energy_index = int(np.argmax(np.abs(site_residual)))
        worst_force_index = int(np.argmax(site_force_error))
        site_summary["worst_energy_point"] = {
            "z": float(z_values[worst_energy_index]),
            "teacher_energy": float(teacher_energy[sl][worst_energy_index]),
            "student_energy": float(student_energy[sl][worst_energy_index]),
            "energy_residual": float(site_residual[worst_energy_index]),
        }
        site_summary["worst_force_point"] = {
            "z": float(z_values[worst_force_index]),
            "teacher_force_norm": float(np.linalg.norm(teacher_force[sl][worst_force_index].reshape(-1))),
            "student_force_norm": float(np.linalg.norm(student_force[sl][worst_force_index].reshape(-1))),
            "force_error_norm": float(site_force_error[worst_force_index]),
        }
        site_window_masks = _window_masks(z_values, primary_z_min=float(args.z_min), stress_z_min=float(args.stress_z_min), z_max=float(args.z_max))
        site_summary["windows"] = {
            "primary": {
                "all": _metrics_for_mask(teacher_energy[sl], student_energy[sl], teacher_force[sl], student_force[sl], site_window_masks["primary"]),
                "train": _metrics_for_mask(teacher_energy[sl], student_energy[sl], teacher_force[sl], student_force[sl], site_window_masks["primary"] & site_train_mask),
                "holdout": _metrics_for_mask(teacher_energy[sl], student_energy[sl], teacher_force[sl], student_force[sl], site_window_masks["primary"] & (~site_train_mask)),
            },
            "stress": {
                "all": _metrics_for_mask(teacher_energy[sl], student_energy[sl], teacher_force[sl], student_force[sl], site_window_masks["stress"]),
            },
            "full": {
                "all": _metrics_for_mask(teacher_energy[sl], student_energy[sl], teacher_force[sl], student_force[sl], site_window_masks["full"]),
            },
        }
        summary["sites"][site_label] = site_summary

        np.savez_compressed(
            traces_dir / f"{site_label}_zscan_trace.npz",
            z_values=z_values,
            teacher_energy=np.asarray(teacher_energy[sl], dtype=float),
            student_energy=np.asarray(student_energy[sl], dtype=float),
            teacher_force=np.asarray(teacher_force[sl], dtype=float),
            student_force=np.asarray(student_force[sl], dtype=float),
            pauli=np.asarray(site_components["pauli"], dtype=float),
            london=np.asarray(site_components["london"], dtype=float),
            coulomb=np.asarray(site_components["coulomb"], dtype=float),
            train_mask=np.asarray(site_train_mask, dtype=bool),
        )
        _plot_site_traces(
            out_dir=plots_dir,
            site_label=site_label,
            z_values=z_values,
            teacher_energy=teacher_energy[sl],
            student_energy=student_energy[sl],
            teacher_force=teacher_force[sl],
            student_force=student_force[sl],
            components=site_components,
            train_mask=site_train_mask,
            primary_z_min=float(args.z_min),
            stress_z_min=float(args.stress_z_min),
        )

    if heldout_requested:
        heldout_poses = _build_zscan_pose_batch(density, adsorbate, site_uvs, z_values, heldout_quaternion)
        heldout_teacher, heldout_teacher_elapsed = _evaluate_teacher_chunked(
            teacher_backend=teacher_backend,
            density=density,
            poses=heldout_poses,
            chunk_size=args.teacher_chunk_size,
        )
        heldout_student, heldout_student_elapsed = _evaluate_student_chunked(
            model=model,
            poses=heldout_poses,
            params=fit_result.params,
            chunk_size=args.student_chunk_size,
        )
        heldout_masks = _window_masks(
            heldout_poses.pose_params[:, 2],
            primary_z_min=float(args.z_min),
            stress_z_min=float(args.stress_z_min),
            z_max=float(args.z_max),
        )
        summary["heldout_orientation"]["timing"] = {
            "teacher_seconds_total": float(heldout_teacher_elapsed),
            "student_seconds_total": float(heldout_student_elapsed),
        }
        summary["heldout_orientation"]["windows"] = {
            "primary": {
                "all": _metrics_for_mask(heldout_teacher.energies, heldout_student.energies, heldout_teacher.forces, heldout_student.forces, heldout_masks["primary"]),
            },
            "stress": {
                "all": _metrics_for_mask(heldout_teacher.energies, heldout_student.energies, heldout_teacher.forces, heldout_student.forces, heldout_masks["stress"]),
            },
            "full": {
                "all": _metrics_for_mask(heldout_teacher.energies, heldout_student.energies, heldout_teacher.forces, heldout_student.forces, heldout_masks["full"]),
            },
        }

    save_json(summary, out_dir / "zscan_summary.json")
    _write_stage(out_dir, "complete")
    print(f"Ag z scan benchmark -> {out_dir}")


if __name__ == "__main__":
    main()
