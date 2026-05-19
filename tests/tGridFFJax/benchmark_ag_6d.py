#!/usr/bin/env python3
"""6D apples-to-apples benchmark: MAD-SURF vs fitted GridFF on the SAME poses.

Pipeline
--------
1. Build a 6D rigid pose batch (u, v, z, quaternion) over the surface unit cell.
2. Evaluate the MAD-SURF (MACE) teacher at every pose -> E_teacher, F_teacher.
3. Stratified train/val/test split (by z-bin and by site).
4. Fit GridFF parameters on the train split (Stage 1 PLQ + Stage 2 raw 1/r^6 lstsq).
5. Re-evaluate the FITTED GridFF at the SAME poses -> E_student, F_student.
6. Diagnostics (per split, per z-bin, per site, per channel) + plots.

Output (under --out-dir):
    apples_dataset.npz   all positions, pose_params, energies, forces, split idx
    fit_params.json      fitted GridFF parameters
    metrics.json         per-split comparison stats
    report.txt           human-readable summary
    figures/             parity plots, per-z RMSE, 6D slice heatmaps, channels

The teacher run is cacheable: if apples_dataset.npz exists with a matching pose
hash, the teacher step is skipped (use --rebuild-teacher to force).

Typical CLI
-----------
  # Smoke run (~250 poses, ~1-2 min on CUDA):
  python benchmark_ag_6d.py --smoke --adsorbate CHONH2 \
      --out-dir tests/tGridFFJax/runs/ag6d_chonh2_smoke

  # Production CHONH2 (~1900 poses, ~3-5 min on CUDA):
  python benchmark_ag_6d.py --adsorbate CHONH2 \
      --n-u 6 --n-v 6 --n-z 12 --n-orient 4 --n-yaw 2 \
      --out-dir tests/tGridFFJax/runs/ag6d_chonh2

  # PTCDA (auto teacher tiling for the supercell MACE call):
  python benchmark_ag_6d.py --adsorbate cpp/common_resources/xyz/PTCDA.xyz \
      --teacher-tile auto --n-u 4 --n-v 4 --n-z 10 --n-orient 2 --n-yaw 1 \
      --out-dir tests/tGridFFJax/runs/ag6d_ptcda
"""

from __future__ import annotations

import argparse
import hashlib
import sys
import time
from dataclasses import asdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.append(str(ROOT))

from pyBall.gridff_jax import RunConfig
from pyBall.gridff_jax.density_backends import make_density_backend
from pyBall.gridff_jax.diagnostics import MethodResult, compare_methods
from pyBall.gridff_jax.fit import fit_two_stage_c6
from pyBall.gridff_jax.hybrid_energy import HybridGridFFModel
from pyBall.gridff_jax.interfaces import PoseBatch, TeacherResult
from pyBall.gridff_jax.pose_sampling.molecules import (
    benchmark_adsorbates,
    load_adsorbate_from_xyz,
)
from pyBall.gridff_jax.pose_sampling.sampler_6d import (
    Sampler6DConfig,
    _transform_adsorbate,
    build_6d_pose_batch,
)
from pyBall.gridff_jax.pose_sampling.rigid import infer_reference_planes
from pyBall.gridff_jax.teacher_backends import make_teacher_backend
from pyBall.gridff_jax.utils import ensure_dir, save_json


DEFAULT_CHGCAR = "/home/niel/git/ORR_HER_Ag_Colab/results/Ag_ORR_HER/slab_clean/final_scf_12x12x1/CHGCAR"
DEFAULT_LOCPOT = "/home/niel/git/ORR_HER_Ag_Colab/results/Ag_ORR_HER/slab_clean/workfunc_12x12x1/LOCPOT"
DEFAULT_MODEL = str(
    Path(__file__).resolve().parent
    / "mad-surf_data" / "models" / "full_dataset_config_weights" / "MACE_model.model"
)


# ---------------------------------------------------------------------------
#  CLI
# ---------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    # Data paths
    p.add_argument("--chgcar", default=DEFAULT_CHGCAR)
    p.add_argument("--locpot", default=DEFAULT_LOCPOT)
    p.add_argument("--model-path", default=DEFAULT_MODEL)
    p.add_argument("--adsorbate", default="CHONH2",
                   help="Built-in name (CHONH2, H2O, ...) or path to XYZ file")
    p.add_argument("--out-dir", required=True)

    # Pose grid
    p.add_argument("--n-u", type=int, default=6)
    p.add_argument("--n-v", type=int, default=6)
    p.add_argument("--n-z", type=int, default=12)
    p.add_argument("--z-min", type=float, default=2.0)
    p.add_argument("--z-max", type=float, default=5.5)
    p.add_argument("--z-bias", type=float, default=1.5)
    p.add_argument("--n-orient", type=int, default=4)
    p.add_argument("--n-yaw", type=int, default=2)
    p.add_argument("--tilt-max-deg", type=float, default=60.0)
    p.add_argument("--random-fraction", type=float, default=0.1)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--smoke", action="store_true",
                   help="Override grid to small values for fast end-to-end test")

    # Splits
    p.add_argument("--val-fraction", type=float, default=0.1)
    p.add_argument("--test-fraction", type=float, default=0.1)
    p.add_argument("--z-bins-stratify", type=int, default=5)

    # Teacher
    p.add_argument("--device", default="cuda",
                   help="MAD-SURF MACE device: cuda / cpu")
    p.add_argument("--teacher-tile", default="auto",
                   help='Tiling for teacher MACE eval: "auto" (default) or "NX,NY"')
    p.add_argument("--max-abs-teacher-eV", type=float, default=3.0,
                   help="Drop poses where |E_teacher| > this many eV (data hygiene). "
                        "Catastrophic overlaps at low z + high tilt produce >>1 eV "
                        "interaction energies that pull the fit to bad minima.")
    p.add_argument("--teacher-chunk", type=int, default=64)
    p.add_argument("--rebuild-teacher", action="store_true",
                   help="Force re-evaluation of teacher even if cached dataset exists")
    p.add_argument("--load-apples", default=None,
                   help="Skip pose generation, teacher eval, and outlier filter. "
                        "Load poses + teacher + splits from this apples_dataset.npz file. "
                        "Used to refit on existing data with a different --fit-mode.")

    # Fit
    p.add_argument("--london-damp-d0", type=float, default=3.70)
    p.add_argument("--london-damp-width", type=float, default=0.35)
    p.add_argument("--alpha-morse", type=float, default=1.5)
    p.add_argument("--c6-rcut", type=float, default=15.0)
    p.add_argument("--max-steps", type=int, default=400)
    p.add_argument("--learning-rate", type=float, default=1.0e-2)
    p.add_argument("--fit-mode", choices=["energy", "joint"], default="joint",
                   help="Stage 2 lstsq objective: 'energy' (E only) or 'joint' (E + F).")
    p.add_argument("--force-weight", type=float, default=1.0,
                   help="Weight on force residuals in joint lstsq. Force has 3N "
                        "components per pose vs 1 energy. w=1.0 lets each force "
                        "component contribute equally to one energy data point.")

    # Student
    p.add_argument("--student-chunk", type=int, default=64)
    p.add_argument("--prefer-jax", action="store_true", default=True)

    return p.parse_args()


def apply_smoke_overrides(args):
    args.n_u = 4
    args.n_v = 4
    args.n_z = 8
    args.n_orient = 2
    args.n_yaw = 1
    args.random_fraction = 0.0


# ---------------------------------------------------------------------------
#  Adsorbate loader
# ---------------------------------------------------------------------------
def load_adsorbate(spec: str):
    builtins = benchmark_adsorbates()
    if spec in builtins:
        return builtins[spec]
    p = Path(spec)
    if p.suffix.lower() in (".xyz", ".mol2"):
        return load_adsorbate_from_xyz(xyz_path=str(p), name=p.stem,
                                       anchor_index=0, use_input_charges=False)
    raise ValueError(f"Adsorbate '{spec}' is neither a built-in nor an XYZ path")


# ---------------------------------------------------------------------------
#  Config builder (mirrors the 16 meV recipe)
# ---------------------------------------------------------------------------
def build_run_config(args) -> RunConfig:
    cfg = RunConfig()
    cfg.density_backend.kind = "vasp_volumetric"
    cfg.density_backend.chgcar_path = args.chgcar
    cfg.density_backend.locpot_path = args.locpot

    cfg.teacher_backend.kind = "madsurf"
    cfg.teacher_backend.model_path = args.model_path
    cfg.teacher_backend.device = args.device
    cfg.teacher_backend.default_dtype = "float64"
    cfg.teacher_backend.interaction_energy = True
    tile = args.teacher_tile.strip().lower()
    if tile == "auto":
        cfg.teacher_backend.teacher_tile = (0, 0)
    else:
        parts = [int(x) for x in tile.split(",")]
        if len(parts) != 2:
            raise ValueError(f"--teacher-tile must be 'auto' or 'NX,NY', got {args.teacher_tile!r}")
        cfg.teacher_backend.teacher_tile = (parts[0], parts[1])

    cfg.grid.builder_mode = "metal_density_plq"
    cfg.grid.alpha_morse = float(args.alpha_morse)
    cfg.grid.london_damping_d0 = float(args.london_damp_d0)
    cfg.grid.london_damping_width = float(args.london_damp_width)
    cfg.grid.use_pairwise_c6 = True
    cfg.grid.c6_rcut = float(args.c6_rcut)

    cfg.toggles.use_ct_qeq = False
    cfg.toggles.use_image_charge = False
    cfg.toggles.use_reactive_grid = False

    cfg.hybrid_model.use_qeq = False
    cfg.hybrid_model.use_image = False
    cfg.hybrid_model.use_reactive_grid = False
    cfg.hybrid_model.use_req_plq = True

    cfg.training.fit_req_radius_offset = True
    cfg.training.fit_req_energy_scale = True
    cfg.training.fit_static_charge = False
    cfg.training.fit_chi = False
    cfg.training.fit_hardness = False
    cfg.training.fit_image_plane = False
    cfg.training.fit_reactive = False
    cfg.training.fit_c6_coeff = True
    cfg.training.two_stage_c6 = True
    cfg.training.max_steps = int(args.max_steps)
    cfg.training.learning_rate = float(args.learning_rate)
    cfg.training.validation_split = 0.0
    cfg.training.test_split = 0.0
    cfg.training.req_radius_regularization = 5.0e-2
    cfg.training.req_energy_regularization = 5.0e-3
    cfg.training.constraint_bound_fraction = 5.0e-2
    return cfg


# ---------------------------------------------------------------------------
#  Pose helpers
# ---------------------------------------------------------------------------
def pose_signature(poses: PoseBatch, adsorbate_name: str) -> str:
    """Stable hash of the pose grid + adsorbate for cache validation."""
    h = hashlib.sha256()
    h.update(adsorbate_name.encode())
    h.update(np.ascontiguousarray(poses.pose_params, dtype=np.float64).tobytes())
    h.update(np.ascontiguousarray(poses.positions, dtype=np.float64).tobytes())
    return h.hexdigest()[:16]


def slice_poses(poses: PoseBatch, idx: np.ndarray) -> PoseBatch:
    return PoseBatch(
        adsorbate=poses.adsorbate,
        positions=poses.positions[idx],
        pose_params=poses.pose_params[idx],
        site_labels=[poses.site_labels[int(i)] for i in idx],
        metadata=dict(poses.metadata),
    )


def slice_teacher(t: TeacherResult, idx: np.ndarray) -> TeacherResult:
    return TeacherResult(
        energies=np.asarray(t.energies, dtype=float)[idx],
        forces=np.asarray(t.forces, dtype=float)[idx],
        metadata=dict(t.metadata),
    )


# ---------------------------------------------------------------------------
#  Stratified split: by (site, z-bin)
# ---------------------------------------------------------------------------
def stratified_split_indices(poses: PoseBatch, val_frac: float, test_frac: float,
                             z_bins: int, seed: int):
    rng = np.random.default_rng(seed)
    n = len(poses.positions)
    idx_all = np.arange(n)
    z_vals = poses.pose_params[:, 2]
    labels = np.asarray(poses.site_labels)

    z_edges = np.linspace(z_vals.min(), z_vals.max(), z_bins + 1)
    z_bin = np.clip(np.digitize(z_vals, z_edges[1:-1]), 0, z_bins - 1)
    label_to_id = {lab: i for i, lab in enumerate(sorted(set(labels)))}
    label_id = np.array([label_to_id[lab] for lab in labels])
    stratum = z_bin * len(label_to_id) + label_id

    train_idx, val_idx, test_idx = [], [], []
    for s in np.unique(stratum):
        members = idx_all[stratum == s]
        rng.shuffle(members)
        n_s = len(members)
        n_test = max(1, int(round(n_s * test_frac)))
        n_val = max(1, int(round(n_s * val_frac)))
        n_train = max(1, n_s - n_test - n_val)
        if n_train + n_val + n_test > n_s:
            n_test = max(0, n_s - n_train - n_val)
        test_idx.extend(members[:n_test])
        val_idx.extend(members[n_test:n_test + n_val])
        train_idx.extend(members[n_test + n_val:])
    return (np.asarray(train_idx, dtype=int),
            np.asarray(val_idx, dtype=int),
            np.asarray(test_idx, dtype=int))


# ---------------------------------------------------------------------------
#  Teacher / student evaluation (chunked)
# ---------------------------------------------------------------------------
def evaluate_teacher(teacher_backend, density, poses, chunk: int):
    n = len(poses.positions)
    if n <= chunk:
        t0 = time.perf_counter()
        r = teacher_backend.evaluate_batch(density, poses)
        return r, time.perf_counter() - t0
    es, fs = [], []
    t0 = time.perf_counter()
    for s in range(0, n, chunk):
        e = min(s + chunk, n)
        idx = np.arange(s, e, dtype=int)
        r = teacher_backend.evaluate_batch(density, slice_poses(poses, idx))
        es.append(np.asarray(r.energies, dtype=float))
        fs.append(np.asarray(r.forces, dtype=float))
        print(f"  [teacher] {e}/{n} ({100.0 * e / n:.0f}%)", flush=True)
    elapsed = time.perf_counter() - t0
    return TeacherResult(
        energies=np.concatenate(es), forces=np.concatenate(fs),
        metadata={"teacher_eval_seconds": elapsed,
                  "seconds_per_pose": elapsed / max(n, 1)},
    ), elapsed


def evaluate_student(model, poses, params, chunk: int):
    n = len(poses.positions)
    es, fs, comps = [], [], {}
    t0 = time.perf_counter()
    for s in range(0, n, chunk):
        e = min(s + chunk, n)
        idx = np.arange(s, e, dtype=int)
        r = model.evaluate_batch(slice_poses(poses, idx), params=params, compute_forces=True)
        es.append(np.asarray(r.energies, dtype=float))
        fs.append(np.asarray(r.forces, dtype=float))
        for k, v in r.components.items():
            comps.setdefault(k, []).append(np.asarray(v, dtype=float))
    elapsed = time.perf_counter() - t0
    return {
        "energies": np.concatenate(es),
        "forces": np.concatenate(fs),
        "components": {k: np.concatenate(v) for k, v in comps.items()},
        "elapsed": elapsed,
    }


# ---------------------------------------------------------------------------
#  Cache load / save (apples_dataset.npz)
# ---------------------------------------------------------------------------
def save_apples(out_dir: Path, poses, teacher, splits, sig, density,
                student=None, fit_params_payload=None):
    np.savez_compressed(
        out_dir / "apples_dataset.npz",
        positions=poses.positions,
        pose_params=poses.pose_params,
        site_labels=np.asarray(poses.site_labels),
        E_teacher=teacher.energies,
        F_teacher=teacher.forces,
        E_student=(student["energies"] if student is not None else np.zeros(0)),
        F_student=(student["forces"] if student is not None else np.zeros(0)),
        train_idx=splits[0], val_idx=splits[1], test_idx=splits[2],
        adsorbate_symbols=np.asarray(poses.adsorbate.symbols),
        adsorbate_positions=poses.adsorbate.positions,
        adsorbate_anchor=np.int32(poses.adsorbate.anchor_index),
        cell=density.cell, origin=density.origin,
        signature=sig,
    )
    if student is not None and student["components"]:
        np.savez_compressed(out_dir / "student_components.npz", **student["components"])
    if fit_params_payload is not None:
        save_json(fit_params_payload, out_dir / "fit_params.json")


def load_cached_teacher(out_dir: Path, sig: str):
    p = out_dir / "apples_dataset.npz"
    if not p.exists():
        return None
    d = np.load(p, allow_pickle=True)
    if str(d["signature"]) != sig:
        return None
    return TeacherResult(
        energies=np.asarray(d["E_teacher"], dtype=float),
        forces=np.asarray(d["F_teacher"], dtype=float),
        metadata={"loaded_from_cache": True},
    )


# ---------------------------------------------------------------------------
#  Diagnostics
# ---------------------------------------------------------------------------
def split_stats(poses, teacher, student, idx):
    if len(idx) == 0:
        return {}
    ref = MethodResult(
        energies=teacher.energies[idx],
        forces=teacher.forces[idx],
        pose_params=poses.pose_params[idx],
        site_labels=[poses.site_labels[int(i)] for i in idx],
    )
    pred = MethodResult(
        energies=student["energies"][idx],
        forces=student["forces"][idx],
        pose_params=poses.pose_params[idx],
        site_labels=[poses.site_labels[int(i)] for i in idx],
    )
    cmp_ = compare_methods({"MAD-SURF": ref, "GridFF": pred}, reference="MAD-SURF")
    return asdict(cmp_["GridFF"])


# ---------------------------------------------------------------------------
#  Plots
# ---------------------------------------------------------------------------
def plot_parity_energy(_poses, teacher, student, splits, out_path):
    fig, ax = plt.subplots(figsize=(5.4, 5.4))
    colors = {"train": "#888888", "val": "#ff7f0e", "test": "#d62728"}
    for name, idx, marker in [("train", splits[0], "."),
                              ("val", splits[1], "s"),
                              ("test", splits[2], "o")]:
        if len(idx) == 0:
            continue
        x = teacher.energies[idx]
        y = student["energies"][idx]
        rmse = np.sqrt(np.mean((y - x) ** 2)) * 1000.0
        ax.scatter(x, y, c=colors[name], s=10, marker=marker, alpha=0.55,
                   label=f"{name} (n={len(idx)}, RMSE={rmse:.1f} meV)")
    lo = min(teacher.energies.min(), student["energies"].min())
    hi = max(teacher.energies.max(), student["energies"].max())
    ax.plot([lo, hi], [lo, hi], "k--", lw=0.8)
    ax.set_xlabel("MAD-SURF E [eV]")
    ax.set_ylabel("GridFF E [eV]")
    ax.set_title("Energy parity (apples-to-apples)")
    ax.legend(loc="upper left", fontsize=8, framealpha=0.85)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_parity_force(_poses, teacher, student, splits, out_path):
    fig, ax = plt.subplots(figsize=(5.4, 5.4))
    colors = {"train": "#888888", "val": "#ff7f0e", "test": "#d62728"}
    f_ref_norm = np.linalg.norm(teacher.forces, axis=-1).reshape(len(teacher.forces), -1)
    f_pred_norm = np.linalg.norm(student["forces"], axis=-1).reshape(len(student["forces"]), -1)
    for name, idx, marker in [("train", splits[0], "."),
                              ("val", splits[1], "s"),
                              ("test", splits[2], "o")]:
        if len(idx) == 0:
            continue
        x = f_ref_norm[idx].ravel()
        y = f_pred_norm[idx].ravel()
        rmse = np.sqrt(np.mean(
            (student["forces"][idx] - teacher.forces[idx]) ** 2
        ))
        ax.scatter(x, y, c=colors[name], s=6, marker=marker, alpha=0.45,
                   label=f"{name} (per-atom |F| pairs; vec RMSE={rmse:.3f} eV/Å)")
    lo = 0.0
    hi = max(f_ref_norm.max(), f_pred_norm.max())
    ax.plot([lo, hi], [lo, hi], "k--", lw=0.8)
    ax.set_xlabel("MAD-SURF |F_atom| [eV/Å]")
    ax.set_ylabel("GridFF |F_atom| [eV/Å]")
    ax.set_title("Force-magnitude parity (per atom)")
    ax.legend(loc="upper left", fontsize=8, framealpha=0.85)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_per_z_rmse(poses, teacher, student, splits, out_path, n_bins=8):
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    z_all = poses.pose_params[:, 2]
    edges = np.linspace(z_all.min(), z_all.max(), n_bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    width = (edges[1] - edges[0]) / 4.0
    colors = {"train": "#888888", "val": "#ff7f0e", "test": "#d62728"}
    offsets = {"train": -width, "val": 0.0, "test": +width}
    for name, idx in zip(["train", "val", "test"], splits):
        if len(idx) == 0:
            continue
        z = poses.pose_params[idx, 2]
        err = (student["energies"][idx] - teacher.energies[idx]) * 1000.0
        b = np.digitize(z, edges[1:-1])
        rmse = np.array([
            np.sqrt(np.mean(err[b == i] ** 2)) if np.any(b == i) else 0.0
            for i in range(n_bins)
        ])
        ax.bar(centers + offsets[name], rmse, width, color=colors[name],
               label=name, alpha=0.85)
    ax.set_xlabel("z height [Å]")
    ax.set_ylabel("Energy RMSE [meV]")
    ax.set_title("Apples-to-apples RMSE vs z (per split)")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_per_site_rmse(poses, teacher, student, out_path):
    labels = np.asarray(poses.site_labels)
    unique = sorted(set(labels.tolist()))
    err = (student["energies"] - teacher.energies) * 1000.0
    rmse = [np.sqrt(np.mean(err[labels == lab] ** 2)) for lab in unique]
    counts = [int(np.sum(labels == lab)) for lab in unique]
    fig, ax = plt.subplots(figsize=(max(4.0, 0.5 * len(unique) + 2), 3.6))
    bars = ax.bar(unique, rmse, color="#1f77b4")
    for b, c in zip(bars, counts):
        ax.text(b.get_x() + b.get_width() / 2, b.get_height(), f"n={c}",
                ha="center", va="bottom", fontsize=8)
    ax.set_ylabel("Energy RMSE [meV]")
    ax.set_title("Per-site energy RMSE (all poses)")
    ax.grid(True, axis="y", alpha=0.3)
    fig.autofmt_xdate(rotation=30)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_uv_slice(density, adsorbate, model, fit_params, teacher_backend,
                  z_height: float, n_grid: int, out_path):
    """Apples-to-apples 2D slice: teacher E(u,v) vs student E(u,v) at fixed z, identity quat."""
    z_ref = float(infer_reference_planes(density)["z_ref"])
    quat = np.array([1.0, 0.0, 0.0, 0.0])
    u = np.linspace(0.0, 1.0, n_grid, endpoint=False)
    v = np.linspace(0.0, 1.0, n_grid, endpoint=False)
    UU, VV = np.meshgrid(u, v, indexing="ij")
    pose_rows = []
    cart = []
    for ui, vi in zip(UU.ravel(), VV.ravel()):
        pose_rows.append(np.concatenate([[ui, vi, z_height], quat]))
        cart.append(_transform_adsorbate(adsorbate, density,
                                         np.array([ui, vi]), z_height, quat, z_ref))
    poses = PoseBatch(
        adsorbate=adsorbate,
        positions=np.asarray(cart, dtype=float),
        pose_params=np.asarray(pose_rows, dtype=float),
        site_labels=["slice"] * len(pose_rows),
        metadata={"z_ref": z_ref, "z_image": float(infer_reference_planes(density)["z_image"])},
    )
    print(f"  [slice] teacher eval @ z={z_height:.2f} Å ({len(poses.positions)} poses)", flush=True)
    t_res, _ = evaluate_teacher(teacher_backend, density, poses, chunk=64)
    s_res = evaluate_student(model, poses, fit_params, chunk=64)
    Et = t_res.energies.reshape(n_grid, n_grid) * 1000.0
    Es = s_res["energies"].reshape(n_grid, n_grid) * 1000.0
    Ed = Es - Et

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.0), constrained_layout=True)
    vmin = min(Et.min(), Es.min())
    vmax = max(Et.max(), Es.max())
    extent = [0, 1, 0, 1]
    for ax, M, title in [(axes[0], Et, f"MAD-SURF E(u,v)  z={z_height:.2f}Å"),
                         (axes[1], Es, "GridFF E(u,v)")]:
        im = ax.imshow(M.T, origin="lower", extent=extent, aspect="equal",
                       vmin=vmin, vmax=vmax, cmap="viridis")
        ax.set_xlabel("u (frac)"); ax.set_ylabel("v (frac)")
        ax.set_title(f"{title}\n[corrugation {M.max() - M.min():.1f} meV]", fontsize=10)
        plt.colorbar(im, ax=ax, label="E [meV]")
    em = max(abs(Ed.min()), abs(Ed.max()))
    im = axes[2].imshow(Ed.T, origin="lower", extent=extent, aspect="equal",
                        vmin=-em, vmax=em, cmap="RdBu_r")
    axes[2].set_xlabel("u (frac)"); axes[2].set_ylabel("v (frac)")
    rmse = np.sqrt(np.mean(Ed ** 2))
    axes[2].set_title(f"GridFF − MAD-SURF\n[RMSE {rmse:.1f} meV]", fontsize=10)
    plt.colorbar(im, ax=axes[2], label="ΔE [meV]")
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return {"corrugation_teacher_meV": float(Et.max() - Et.min()),
            "corrugation_student_meV": float(Es.max() - Es.min()),
            "slice_rmse_meV": float(rmse),
            "z_height": float(z_height)}


def plot_components_vs_z(density, adsorbate, model, fit_params, teacher_backend,
                         uv: tuple, z_values: np.ndarray, out_path):
    """Channel-by-channel decomposition vs z at one fixed (u,v) and identity orientation."""
    z_ref = float(infer_reference_planes(density)["z_ref"])
    quat = np.array([1.0, 0.0, 0.0, 0.0])
    pose_rows, cart = [], []
    for z in z_values:
        pose_rows.append(np.concatenate([[uv[0], uv[1], z], quat]))
        cart.append(_transform_adsorbate(adsorbate, density,
                                         np.array(uv), float(z), quat, z_ref))
    poses = PoseBatch(
        adsorbate=adsorbate,
        positions=np.asarray(cart, dtype=float),
        pose_params=np.asarray(pose_rows, dtype=float),
        site_labels=["zscan"] * len(z_values),
        metadata={"z_ref": z_ref, "z_image": float(infer_reference_planes(density)["z_image"])},
    )
    t_res, _ = evaluate_teacher(teacher_backend, density, poses, chunk=64)
    s_res = evaluate_student(model, poses, fit_params, chunk=64)
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    ax.plot(z_values, t_res.energies * 1000, "k-", lw=2, label="MAD-SURF total")
    ax.plot(z_values, s_res["energies"] * 1000, "C0--", lw=2, label="GridFF total")
    palette = ["#9467bd", "#2ca02c", "#ff7f0e", "#8c564b", "#17becf", "#e377c2"]
    for k, (name, vals) in enumerate(s_res["components"].items()):
        ax.plot(z_values, np.asarray(vals).ravel() * 1000,
                color=palette[k % len(palette)], lw=1.0, alpha=0.8,
                label=f"GridFF {name}")
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xlabel("z [Å]")
    ax.set_ylabel("E [meV]")
    ax.set_title(f"Channel decomposition at (u,v)={uv}, identity quaternion")
    ax.legend(fontsize=8, loc="best")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
#  Force-augmented Stage 2 lstsq:  E_disp and F_disp are both linear in C₆.
#  Build the design matrix once via per-element probe evaluations, then solve
#  the joint system  [E_residual ; w·F_residual]  for [c6_1..c6_n , offset].
# ---------------------------------------------------------------------------
def fit_c6_joint_lstsq(model, poses, teacher, baseline_params, force_weight=1.0,
                       chunk: int = 64, verbose: bool = True):
    """Fit C₆ + offset to BOTH energy and force residuals simultaneously.

    Returns (fitted_params, diagnostics_dict). The dict contains the baseline
    eval (c6=0) and the per-element probe deltas, which are reused for the
    fit-quality plots downstream.
    """
    from copy import deepcopy

    elements = sorted(baseline_params.c6_coeff.keys())
    n_el = len(elements)

    # ---- 1. Baseline eval (c6 = 0) ----
    base_p = deepcopy(baseline_params)
    for el in base_p.c6_coeff:
        base_p.c6_coeff[el] = 0.0
    base_p.energy_offset = 0.0
    base = evaluate_student(model, poses, base_p, chunk)
    E_base = base["energies"]
    F_base = base["forces"]
    n_pose, n_atom, _ = F_base.shape

    # ---- 2. Per-element probe evaluations: c6_el = 1, others 0 ----
    X_E = np.zeros((n_pose, n_el), dtype=float)
    X_F = np.zeros((n_pose, n_atom, 3, n_el), dtype=float)
    for j, el in enumerate(elements):
        probe_p = deepcopy(base_p)
        probe_p.c6_coeff[el] = 1.0
        probe = evaluate_student(model, poses, probe_p, chunk)
        X_E[:, j] = probe["energies"] - E_base
        X_F[:, :, :, j] = probe["forces"] - F_base
        if verbose:
            print(f"  [joint lstsq] probe element {el}: "
                  f"|E_probe-E_base| max = {np.abs(X_E[:, j]).max() * 1000:.1f} meV, "
                  f"|F_probe-F_base| max = {np.abs(X_F[:, :, :, j]).max():.4f} eV/Å",
                  flush=True)

    # ---- 3. Targets ----
    E_target = np.asarray(teacher.energies, dtype=float)
    F_target = np.asarray(teacher.forces, dtype=float)
    E_residual = E_target - E_base
    F_residual = F_target - F_base

    # ---- 4. Build joint linear system ----
    # parameters d = [c6_1, ..., c6_n, offset]
    A_E = np.column_stack([X_E, np.ones(n_pose)])           # (n_pose, n_el+1)
    b_E = E_residual                                          # (n_pose,)
    A_F = np.zeros((n_pose * n_atom * 3, n_el + 1))
    A_F[:, :n_el] = X_F.reshape(-1, n_el)                    # offset col stays 0
    b_F = F_residual.reshape(-1)                             # (n_pose*n_atom*3,)

    w = float(force_weight)
    A = np.vstack([A_E, w * A_F])
    b = np.concatenate([b_E, w * b_F])

    # ---- 5. Solve ----
    D, _, rank, _ = np.linalg.lstsq(A, b, rcond=None)
    c6_values = D[:n_el]
    offset = float(D[n_el])

    # ---- 6. Reconstruct predictions and report ----
    E_pred = E_base + X_E @ c6_values + offset
    F_pred = F_base + np.tensordot(X_F, c6_values, axes=([3], [0]))
    e_rmse_before = np.sqrt(np.mean(E_residual ** 2)) * 1000
    e_rmse_after = np.sqrt(np.mean((E_target - E_pred) ** 2)) * 1000
    f_rmse_before = np.sqrt(np.mean(F_residual ** 2))
    f_rmse_after = np.sqrt(np.mean((F_target - F_pred) ** 2))
    if verbose:
        print(f"  [joint lstsq] system rank={rank}/{n_el + 1}", flush=True)
        print(f"  [joint lstsq] C₆: {dict(zip(elements, c6_values.tolist()))}", flush=True)
        print(f"  [joint lstsq] offset: {offset * 1000:+.2f} meV", flush=True)
        print(f"  [joint lstsq] train E RMSE: {e_rmse_before:.1f} → {e_rmse_after:.1f} meV", flush=True)
        print(f"  [joint lstsq] train F RMSE: {f_rmse_before:.4f} → {f_rmse_after:.4f} eV/Å", flush=True)

    fitted = deepcopy(baseline_params)
    for j, el in enumerate(elements):
        fitted.c6_coeff[el] = float(c6_values[j])
    fitted.energy_offset = offset

    diagnostics = {
        "elements": elements,
        "E_baseline": E_base,
        "F_baseline": F_base,
        "X_E": X_E,
        "X_F": X_F,
        "c6": c6_values,
        "offset": offset,
        "rank": int(rank),
        "force_weight": w,
        "metrics": {
            "train_E_RMSE_before_meV": float(e_rmse_before),
            "train_E_RMSE_after_meV": float(e_rmse_after),
            "train_F_RMSE_before_eV_A": float(f_rmse_before),
            "train_F_RMSE_after_eV_A": float(f_rmse_after),
        },
    }
    return fitted, diagnostics


# ---------------------------------------------------------------------------
#  Stage 1 PLQ-only fit (used as the baseline for both energy-only and joint
#  Stage 2). Re-uses the framework's gradient-descent fitter, but never
#  invokes the framework's own Stage 2 lstsq (we do that ourselves above).
# ---------------------------------------------------------------------------
def fit_stage1_plq_only(density, train_poses, train_teacher, model, training_cfg):
    """Stage 1: REQ + PLQ gradient descent with c6_coeff frozen at 0."""
    from copy import deepcopy
    from pyBall.gridff_jax.fit import fit_hybrid_parameters
    from pyBall.gridff_jax.hybrid_energy import default_hybrid_parameters

    elements = tuple(sorted(set(train_poses.adsorbate.symbols)))
    z_image = train_poses.metadata["z_image"]
    init = default_hybrid_parameters(elements, z_image=z_image)
    init = deepcopy(init)
    for el in init.c6_coeff:
        init.c6_coeff[el] = 0.0

    s1_cfg = deepcopy(training_cfg)
    s1_cfg.fit_c6_coeff = False
    s1_cfg.two_stage_c6 = False  # we control staging here
    return fit_hybrid_parameters(
        density=density,
        datasets=[(train_poses, train_teacher)],
        model=model,
        training=s1_cfg,
        initial_params=init,
    )


# ---------------------------------------------------------------------------
#  Fit-quality plots
# ---------------------------------------------------------------------------
def plot_residual_histograms(teacher, student, splits, out_path):
    """Histogram of energy residuals (train/val/test) and force component residuals."""
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.0))
    colors = {"train": "#888888", "val": "#ff7f0e", "test": "#d62728"}
    for name, idx in zip(["train", "val", "test"], splits):
        if len(idx) == 0:
            continue
        e_err = (student["energies"][idx] - teacher.energies[idx]) * 1000.0  # meV
        rmse_e = float(np.sqrt(np.mean(e_err ** 2)))
        bias_e = float(np.mean(e_err))
        axes[0].hist(e_err, bins=60, alpha=0.55, color=colors[name],
                     label=f"{name} bias={bias_e:+.1f} σ={rmse_e:.1f} meV")
        f_err = (student["forces"][idx] - teacher.forces[idx]).ravel()
        rmse_f = float(np.sqrt(np.mean(f_err ** 2)))
        bias_f = float(np.mean(f_err))
        axes[1].hist(f_err, bins=80, alpha=0.55, color=colors[name],
                     label=f"{name} bias={bias_f:+.3f} σ={rmse_f:.3f} eV/Å")
    axes[0].axvline(0, color="k", lw=0.5)
    axes[1].axvline(0, color="k", lw=0.5)
    axes[0].set_xlabel("E_GridFF − E_MAD-SURF [meV]")
    axes[0].set_ylabel("count"); axes[0].set_title("Energy fit residuals")
    axes[0].legend(fontsize=8); axes[0].grid(True, alpha=0.3)
    axes[1].set_xlabel("F_GridFF − F_MAD-SURF [eV/Å, per component]")
    axes[1].set_ylabel("count"); axes[1].set_title("Force-component fit residuals")
    axes[1].legend(fontsize=8); axes[1].grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_baseline_vs_post(E_baseline, E_post, E_teacher, splits, out_path):
    """Two-panel: PLQ-only baseline vs MAD-SURF, and post-fit (PLQ+C₆+offset) vs MAD-SURF."""
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 5.2), constrained_layout=True)
    colors = {"train": "#888888", "val": "#ff7f0e", "test": "#d62728"}
    for ax, y_arr, title in [
        (axes[0], E_baseline, "Stage 1 baseline (PLQ only, c₆=0, offset=0)"),
        (axes[1], E_post,     "After Stage 2 (PLQ + C₆ + offset)"),
    ]:
        for name, idx in zip(["train", "val", "test"], splits):
            if len(idx) == 0:
                continue
            x = E_teacher[idx]
            y = y_arr[idx]
            rmse = np.sqrt(np.mean((y - x) ** 2)) * 1000.0
            ax.scatter(x, y, c=colors[name], s=10, alpha=0.55,
                       label=f"{name} RMSE={rmse:.1f} meV")
        lo = min(E_teacher.min(), E_baseline.min(), E_post.min())
        hi = max(E_teacher.max(), E_baseline.max(), E_post.max())
        ax.plot([lo, hi], [lo, hi], "k--", lw=0.7)
        ax.set_xlabel("MAD-SURF E [eV]")
        ax.set_ylabel("GridFF E [eV]")
        ax.set_title(title, fontsize=10)
        ax.legend(fontsize=8); ax.grid(True, alpha=0.3)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_force_per_element(symbols, F_teacher, F_student, out_path):
    """Bar chart of force RMSE per element symbol (averaged over atoms+poses+xyz)."""
    symbols = list(symbols)
    unique = sorted(set(symbols))
    F_err = F_student - F_teacher  # (n_pose, n_atom, 3)
    rmse_per_el = []
    counts = []
    for el in unique:
        atom_mask = np.array([s == el for s in symbols])
        if atom_mask.sum() == 0:
            rmse_per_el.append(0.0); counts.append(0)
            continue
        sub = F_err[:, atom_mask, :].ravel()
        rmse_per_el.append(float(np.sqrt(np.mean(sub ** 2))))
        counts.append(int(atom_mask.sum()))
    fig, ax = plt.subplots(figsize=(max(4.5, 0.7 * len(unique) + 2.0), 3.6))
    bars = ax.bar(unique, rmse_per_el, color="#1f77b4")
    for b, c in zip(bars, counts):
        ax.text(b.get_x() + b.get_width() / 2, b.get_height(),
                f"n_atoms={c}", ha="center", va="bottom", fontsize=8)
    ax.set_ylabel("Force RMSE [eV/Å]")
    ax.set_title("Per-element force-fit quality (all poses, all xyz components)")
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_c6_contribution_vs_z(poses, E_baseline, E_post, out_path):
    """Per-pose C₆ + offset correction (E_post − E_baseline) vs z. Shows where the
    Stage 2 correction is doing its work."""
    z_vals = poses.pose_params[:, 2]
    correction_meV = (E_post - E_baseline) * 1000.0
    fig, ax = plt.subplots(figsize=(7.5, 4.4))
    labels = np.asarray(poses.site_labels)
    palette = {"top": "#d62728", "bridge": "#2ca02c", "hollow": "#1f77b4",
               "grid": "#888888", "random": "#ff7f0e"}
    for lab in sorted(set(labels.tolist())):
        m = labels == lab
        ax.scatter(z_vals[m], correction_meV[m], s=10, alpha=0.6,
                   color=palette.get(lab, "#9467bd"), label=lab)
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xlabel("z height [Å]")
    ax.set_ylabel("ΔE = E_post − E_baseline [meV]")
    ax.set_title("Stage 2 correction (C₆ + offset) per pose")
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_fit_loss_curve(history, out_path):
    """Stage 1 PLQ gradient-descent loss curve. history is a list of dicts with
    keys 'step', 'train_loss', and optionally 'val_loss'."""
    if not history:
        return
    steps = [h.get("step", i) for i, h in enumerate(history)]
    train = [h.get("train_loss", np.nan) for h in history]
    val = [h.get("val_loss", np.nan) for h in history]
    fig, ax = plt.subplots(figsize=(7.5, 4.0))
    ax.semilogy(steps, train, "C0-", lw=1.6, label="train loss")
    if any(not np.isnan(v) for v in val):
        ax.semilogy(steps, val, "C1--", lw=1.4, label="val loss")
    ax.set_xlabel("optimization step")
    ax.set_ylabel("loss (energy + force_weight·force MSE)")
    ax.set_title("Stage 1 PLQ fit convergence")
    ax.legend(fontsize=9); ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
#  Main
# ---------------------------------------------------------------------------
def main():
    args = parse_args()
    if args.smoke:
        apply_smoke_overrides(args)

    out_dir = ensure_dir(args.out_dir)
    fig_dir = ensure_dir(out_dir / "figures")

    print("=" * 72)
    print(f"  6D apples-to-apples benchmark  →  {out_dir}")
    print("=" * 72)
    cfg = build_run_config(args)
    save_json(asdict(cfg), out_dir / "run_config.json")

    print("[1/6] Loading density (CHGCAR/LOCPOT)...", flush=True)
    density = make_density_backend(cfg.density_backend).load()
    print(f"   slab atoms={density.natoms}  cell={np.diag(density.cell).round(3)}",
          flush=True)

    print(f"[2/6] Loading adsorbate '{args.adsorbate}'...", flush=True)
    adsorbate = load_adsorbate(args.adsorbate)
    print(f"   atoms={adsorbate.natoms}  symbols={adsorbate.symbols}", flush=True)

    # ---- Fast path: skip teacher and load everything from existing npz ----
    if args.load_apples:
        print(f"[3/6 → 6/6] Loading existing apples dataset from {args.load_apples} ...",
              flush=True)
        d = np.load(args.load_apples, allow_pickle=True)
        from pyBall.gridff_jax.interfaces import AdsorbateDefinition
        ads_loaded = AdsorbateDefinition(
            name=adsorbate.name,
            symbols=[str(s) for s in d["adsorbate_symbols"].tolist()],
            positions=np.asarray(d["adsorbate_positions"], dtype=float),
            anchor_index=int(d["adsorbate_anchor"]),
        )
        adsorbate = ads_loaded
        z_ref_loaded = float(infer_reference_planes(density)["z_ref"])
        z_image_loaded = float(infer_reference_planes(density)["z_image"])
        poses = PoseBatch(
            adsorbate=adsorbate,
            positions=np.asarray(d["positions"], dtype=float),
            pose_params=np.asarray(d["pose_params"], dtype=float),
            site_labels=[str(s) for s in d["site_labels"].tolist()],
            metadata={"z_ref": z_ref_loaded, "z_image": z_image_loaded,
                      "loaded_from": args.load_apples},
        )
        teacher = TeacherResult(
            energies=np.asarray(d["E_teacher"], dtype=float),
            forces=np.asarray(d["F_teacher"], dtype=float),
            metadata={"loaded_from": args.load_apples},
        )
        splits = (np.asarray(d["train_idx"], dtype=int),
                  np.asarray(d["val_idx"], dtype=int),
                  np.asarray(d["test_idx"], dtype=int))
        sig = str(d["signature"])
        n = len(poses.positions)
        teacher_elapsed = 0.0
        teacher_backend = make_teacher_backend(cfg.teacher_backend)
        e_min = float(teacher.energies.min()) * 1000
        e_max = float(teacher.energies.max()) * 1000
        print(f"   {n} poses loaded, signature={sig}", flush=True)
        print(f"   E range = [{e_min:.1f}, {e_max:.1f}] meV  "
              f"(train={len(splits[0])}, val={len(splits[1])}, test={len(splits[2])})",
              flush=True)
    else:
        # ---- Slow path: build poses + teacher + filter ----
        print(f"[3/6] Building 6D pose batch "
              f"(n_u={args.n_u}, n_v={args.n_v}, n_z={args.n_z}, "
              f"n_orient={args.n_orient}, n_yaw={args.n_yaw})...", flush=True)
        sampler_cfg = Sampler6DConfig(
            n_u=args.n_u, n_v=args.n_v, n_z=args.n_z,
            z_range=(args.z_min, args.z_max), z_bias_power=args.z_bias,
            n_orient=args.n_orient, tilt_max_deg=args.tilt_max_deg,
            n_yaw=args.n_yaw, include_high_symmetry_sites=True,
            random_fraction=args.random_fraction, seed=args.seed,
        )
        poses = build_6d_pose_batch(density, adsorbate, sampler_cfg)
        n = len(poses.positions)
        sig = pose_signature(poses, adsorbate.name)
        print(f"   {n} poses generated, signature={sig}", flush=True)

        splits = stratified_split_indices(
            poses, args.val_fraction, args.test_fraction, args.z_bins_stratify, args.seed)
        print(f"   split: train={len(splits[0])}, val={len(splits[1])}, test={len(splits[2])}",
              flush=True)

        print("[4/6] Teacher (MAD-SURF/MACE) evaluation...", flush=True)
        cached = None if args.rebuild_teacher else load_cached_teacher(out_dir, sig)
        if cached is not None:
            teacher = cached
            teacher_elapsed = float(teacher.metadata.get("teacher_eval_seconds", 0.0))
            print(f"   loaded from cache ({len(teacher.energies)} poses)", flush=True)
            teacher_backend = make_teacher_backend(cfg.teacher_backend)
        else:
            teacher_backend = make_teacher_backend(cfg.teacher_backend)
            teacher, teacher_elapsed = evaluate_teacher(
                teacher_backend, density, poses, args.teacher_chunk)
            print(f"   teacher done in {teacher_elapsed:.1f}s "
                  f"({teacher_elapsed / max(n, 1):.3f}s/pose)", flush=True)
            save_apples(out_dir, poses, teacher, splits, sig, density)

        e_min = float(teacher.energies.min()) * 1000
        e_max = float(teacher.energies.max()) * 1000
        print(f"   teacher E range = [{e_min:.1f}, {e_max:.1f}] meV", flush=True)

        # ---- Outlier filter: drop poses with unphysical |E_teacher| ----
        keep_mask = np.abs(teacher.energies) <= args.max_abs_teacher_eV
        n_dropped = int(np.sum(~keep_mask))
        if n_dropped:
            print(f"   [filter] dropping {n_dropped}/{len(teacher.energies)} poses "
                  f"with |E| > {args.max_abs_teacher_eV} eV (overlap artefacts)", flush=True)
            keep_idx = np.where(keep_mask)[0]
            old_to_new = -np.ones(len(teacher.energies), dtype=int)
            old_to_new[keep_idx] = np.arange(len(keep_idx))
            poses = slice_poses(poses, keep_idx)
            teacher = slice_teacher(teacher, keep_idx)
            splits = tuple(
                old_to_new[s[np.isin(s, keep_idx)]] for s in splits
            )
            sig = pose_signature(poses, adsorbate.name)
            e_min = float(teacher.energies.min()) * 1000
            e_max = float(teacher.energies.max()) * 1000
            print(f"   post-filter E range = [{e_min:.1f}, {e_max:.1f}] meV  "
                  f"(n_train={len(splits[0])}, n_val={len(splits[1])}, n_test={len(splits[2])})",
                  flush=True)
            save_apples(out_dir, poses, teacher, splits, sig, density)
        n = len(teacher.energies)

    print(f"[5/6] Fitting GridFF (mode={args.fit_mode!r}, "
          f"force_weight={args.force_weight}) — "
          f"Stage 1 PLQ + Stage 2 raw 1/r⁶ lstsq...", flush=True)
    train_poses = slice_poses(poses, splits[0])
    train_teacher = slice_teacher(teacher, splits[0])
    elements = tuple(sorted(set(adsorbate.symbols)))
    model = HybridGridFFModel(
        density=density, reactive_elements=elements,
        toggles=cfg.toggles, grid_config=cfg.grid,
        hybrid_config=cfg.hybrid_model, prefer_jax=args.prefer_jax,
    )
    model.install_raw_r6_grid()  # raw 1/r^6 basis (the 16 meV recipe)

    fit_t0 = time.perf_counter()
    if args.fit_mode == "joint":
        # Stage 1: PLQ-only gradient descent (uses force_weight from cfg.training)
        stage1 = fit_stage1_plq_only(density, train_poses, train_teacher, model, cfg.training)
        # Stage 2: joint E+F lstsq for C₆ + offset
        print("[two-stage-joint] Stage 2: fitting C₆ to E AND F residuals...", flush=True)
        fitted_params, joint_diag = fit_c6_joint_lstsq(
            model=model, poses=train_poses, teacher=train_teacher,
            baseline_params=stage1.params, force_weight=args.force_weight,
            chunk=args.student_chunk, verbose=True,
        )
        stage1_history = stage1.history
        stage2_metrics = joint_diag["metrics"]
    else:
        # energy-only path: existing fit_two_stage_c6 (Stage 2 = energy-only lstsq)
        fit_result = fit_two_stage_c6(
            density=density,
            datasets=[(train_poses, train_teacher)],
            model=model,
            training=cfg.training,
            initial_params=None,
            lstsq_datasets=[(train_poses, train_teacher)],
            use_raw_r6=True,
        )
        fitted_params = fit_result.params
        stage1_history = fit_result.history
        stage2_metrics = fit_result.metrics
        joint_diag = None

    fit_elapsed = time.perf_counter() - fit_t0
    print(f"   fit done in {fit_elapsed:.1f}s", flush=True)
    fit_payload = {
        "fit_mode": args.fit_mode,
        "force_weight": float(args.force_weight),
        "elements": list(elements),
        "pauli": dict(fitted_params.pauli),
        "london": dict(fitted_params.london),
        "c6_coeff": dict(fitted_params.c6_coeff),
        "energy_offset_meV": float(fitted_params.energy_offset) * 1000,
        "req_radius_offset": dict(fitted_params.req_radius_offset),
        "req_energy_scale": dict(fitted_params.req_energy_scale),
        "metrics": stage2_metrics,
    }
    save_json(fit_payload, out_dir / "fit_params.json")
    fit_result_params = fitted_params  # for use below

    print("[6/6] Re-evaluating fitted GridFF at the SAME 6D poses...", flush=True)
    student = evaluate_student(model, poses, fit_result_params, args.student_chunk)
    print(f"   student done in {student['elapsed']:.1f}s", flush=True)

    # Baseline (PLQ only, c6=0, offset=0) on ALL poses — needed for the
    # baseline-vs-post fit-quality plot. Reuse the joint_diag baseline if
    # available (it was computed on train poses only); otherwise compute fresh.
    from copy import deepcopy as _dc
    baseline_only_params = _dc(fit_result_params)
    for el in baseline_only_params.c6_coeff:
        baseline_only_params.c6_coeff[el] = 0.0
    baseline_only_params.energy_offset = 0.0
    baseline_eval = evaluate_student(model, poses, baseline_only_params, args.student_chunk)
    E_baseline_all = baseline_eval["energies"]

    save_apples(out_dir, poses, teacher, splits, sig, density,
                student=student, fit_params_payload=fit_payload)

    # ---- Diagnostics ----
    print("\n--- Per-split apples-to-apples comparison ---", flush=True)
    metrics = {}
    for name, idx in zip(["train", "val", "test"], splits):
        s = split_stats(poses, teacher, student, idx)
        metrics[name] = s
        if not s:
            print(f"  {name}: (empty split)")
            continue
        print(f"  {name:5s}  n={s['n_poses']:5d}  "
              f"E RMSE={s['energy_rmse_meV']:7.2f} meV  "
              f"E MAE={s['energy_mae_meV']:7.2f} meV  "
              f"E max={s['energy_max_error_meV']:7.2f} meV  "
              f"F RMSE={s['force_rmse_eV_A']:.4f} eV/Å  "
              f"R²={s['energy_r2']:.4f}", flush=True)

    print("\n--- Per-site (all poses) ---", flush=True)
    full = split_stats(poses, teacher, student, np.arange(len(poses.positions)))
    for k, v in (full.get("per_site_energy_rmse") or {}).items():
        print(f"  site {k:8s}  E RMSE = {v:7.2f} meV  (F RMSE = "
              f"{(full.get('per_site_force_rmse') or {}).get(k, 0.0):.4f} eV/Å)",
              flush=True)

    print("\n--- Per-z-bin (all poses) ---", flush=True)
    for k, v in (full.get("z_binned_energy_rmse") or {}).items():
        print(f"  z {k:>10s} Å   E RMSE = {v:7.2f} meV", flush=True)

    print(f"\nCorrugation: teacher = {full.get('corrugation_ref_meV', 0):.1f} meV, "
          f"GridFF = {full.get('corrugation_method_meV', 0):.1f} meV", flush=True)

    metrics["all_poses"] = full
    save_json(metrics, out_dir / "metrics.json")

    # ---- Plots ----
    print("\n[plots] generating figures...", flush=True)
    # Apples-to-apples (prediction quality on the same poses)
    plot_parity_energy(poses, teacher, student, splits, fig_dir / "parity_energy.png")
    plot_parity_force(poses, teacher, student, splits, fig_dir / "parity_force.png")
    plot_per_z_rmse(poses, teacher, student, splits, fig_dir / "per_z_rmse.png")
    plot_per_site_rmse(poses, teacher, student, fig_dir / "per_site_rmse.png")
    # Fit-quality (what the FIT did — residuals, baseline vs post, per element)
    plot_residual_histograms(teacher, student, splits,
                             fig_dir / "fit_residual_histograms.png")
    plot_baseline_vs_post(E_baseline_all, student["energies"], teacher.energies,
                          splits, fig_dir / "fit_baseline_vs_post.png")
    plot_force_per_element(adsorbate.symbols, teacher.forces, student["forces"],
                           fig_dir / "fit_force_per_element.png")
    plot_c6_contribution_vs_z(poses, E_baseline_all, student["energies"],
                              fig_dir / "fit_c6_contribution_vs_z.png")
    plot_fit_loss_curve(stage1_history, fig_dir / "fit_stage1_loss_curve.png")

    # 6D slice + channel decomposition (these re-call the teacher on small batches)
    z_eq = float(np.median(poses.pose_params[poses.pose_params[:, 2] < 3.5, 2]) if
                 np.any(poses.pose_params[:, 2] < 3.5) else 2.5)
    n_slice = 16 if args.smoke else 24
    slice_info = plot_uv_slice(density, adsorbate, model, fit_result_params,
                               teacher_backend, z_height=z_eq, n_grid=n_slice,
                               out_path=fig_dir / "slice_uv.png")
    print(f"   slice corrugation: teacher={slice_info['corrugation_teacher_meV']:.1f} meV, "
          f"GridFF={slice_info['corrugation_student_meV']:.1f} meV, "
          f"slice RMSE={slice_info['slice_rmse_meV']:.1f} meV", flush=True)

    z_ch = np.linspace(args.z_min, args.z_max, 30)
    plot_components_vs_z(density, adsorbate, model, fit_result_params, teacher_backend,
                         uv=(0.0, 0.0), z_values=z_ch,
                         out_path=fig_dir / "channels_vs_z_uv00.png")

    # ---- Human-readable report ----
    with (out_dir / "report.txt").open("w") as fh:
        fh.write("6D apples-to-apples benchmark report\n")
        fh.write(f"  adsorbate: {adsorbate.name} ({adsorbate.natoms} atoms, {elements})\n")
        fh.write(f"  poses:     {n}  signature={sig}\n")
        fh.write(f"  teacher:   {teacher_elapsed:.1f}s ({teacher_elapsed / max(n, 1):.3f}s/pose)\n")
        fh.write(f"  fit:       {fit_elapsed:.1f}s\n")
        fh.write(f"  student:   {student['elapsed']:.1f}s\n\n")
        for name in ("train", "val", "test", "all_poses"):
            s = metrics.get(name, {})
            if not s:
                continue
            fh.write(f"[{name}] n={s.get('n_poses', 0)}  "
                     f"E RMSE={s.get('energy_rmse_meV', 0):.2f} meV  "
                     f"F RMSE={s.get('force_rmse_eV_A', 0):.4f} eV/Å  "
                     f"R²={s.get('energy_r2', 0):.4f}\n")
        fh.write("\nFitted C₆:\n")
        for el, c in fit_payload["c6_coeff"].items():
            fh.write(f"  {el}: {c:+.4f}\n")
        fh.write(f"  energy_offset = {fit_payload['energy_offset_meV']:+.2f} meV\n")
        fh.write(f"\nCorrugation (slice u,v at z={slice_info['z_height']:.2f} Å): "
                 f"teacher={slice_info['corrugation_teacher_meV']:.1f} meV, "
                 f"GridFF={slice_info['corrugation_student_meV']:.1f} meV\n")

    print(f"\n DONE. Outputs in {out_dir}", flush=True)
    print(f"   - apples_dataset.npz     (poses + E/F teacher and student)", flush=True)
    print(f"   - fit_params.json        (fitted GridFF params)", flush=True)
    print(f"   - metrics.json           (per-split, per-site, per-z RMSE)", flush=True)
    print(f"   - report.txt             (human-readable summary)", flush=True)
    print(f"   - figures/               (parity, per-z, per-site, slice, channels)", flush=True)


if __name__ == "__main__":
    main()
