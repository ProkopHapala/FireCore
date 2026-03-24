#!/usr/bin/env python3
"""Lateral scan benchmark: 1D path and 2D heatmap (energy + force).

Loads the fitted raw-r6 model and evaluates GridFF vs MAD-SURF at fixed z
while varying the lateral (u,v) position across one surface unit cell.

Outputs
-------
1D path scan:  top -> bridge -> hollow -> top  (~60 points)
2D grid scan:  25x25 fractional (u,v) grid over [0,1] x [0,1]
Plots:         lateral_1d_{energy,force}.pdf, lateral_2d_{energy,force}_heatmaps.pdf
Data:          lateral_1d_path.npz/.csv, lateral_2d_grid.npz/.csv

Usage
-----
    python lateral_scan.py [--result-dir /tmp/ag_raw_r6_final] [--adsorbate CHONH2]
                           [--z-height 2.7] [--n-path 60] [--n-grid 25]
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("JAX_PLATFORMS", "cpu")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

ROOT = Path("/home/niel/git/FireCore")
sys.path.insert(0, str(ROOT))

from pyBall.gridff_jax import RunConfig
from pyBall.gridff_jax.density_backends import make_density_backend
from pyBall.gridff_jax.hybrid_energy import HybridGridFFModel, HybridParameters
from pyBall.gridff_jax.interfaces import PoseBatch, TeacherResult
from pyBall.gridff_jax.pose_sampling import (
    get_adsorbate,
    infer_reference_planes,
    infer_surface_sites,
    transform_adsorbate,
)
from pyBall.gridff_jax.teacher_backends import make_teacher_backend
from pyBall.gridff_jax.utils import fractional_to_cartesian

# ── Plot style ────────────────────────────────────────────────────────
CLR_MADSURF = "orange"
CLR_GRIDFF  = "#1f77b4"
CLR_ERROR   = "#cc3333"


# ═════════════════════════════════════════════════════════════════════
#  Helpers
# ═════════════════════════════════════════════════════════════════════

def load_config(config_path) -> RunConfig:
    with open(config_path) as fh:
        raw = json.load(fh)
    cfg = RunConfig()
    for section in ("surface", "density_backend", "teacher_backend", "toggles",
                    "grid", "hybrid_model", "training", "export"):
        for key, val in raw.get(section, {}).items():
            obj = getattr(cfg, section)
            if hasattr(obj, key):
                setattr(obj, key, val)
    return cfg


def load_fitted_params(params_path) -> HybridParameters:
    with open(params_path) as fh:
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


def build_model(config, density, adsorbate, prefer_jax=True):
    reactive_elements = tuple(sorted(set(adsorbate.symbols)))
    model = HybridGridFFModel(
        density=density,
        reactive_elements=reactive_elements,
        toggles=config.toggles,
        grid_config=config.grid,
        hybrid_config=config.hybrid_model,
        prefer_jax=prefer_jax,
    )
    print("[lateral] Installing raw 1/r^6 dispersion grid ...", flush=True)
    model.install_raw_r6_grid()
    return model


def _site_uv_lookup(density, labels):
    lookup = {}
    for site in infer_surface_sites(density):
        if site.label in labels and site.label not in lookup:
            lookup[site.label] = np.asarray(site.uv, dtype=float)
    missing = [lab for lab in labels if lab not in lookup]
    if missing:
        raise ValueError(f"Site labels not found: {missing}")
    return lookup


def build_lateral_poses(density, adsorbate, uv_array, z_height, quaternion):
    planes = infer_reference_planes(density)
    z_ref = float(planes["z_ref"])
    positions = []
    pose_params = []
    for uv in uv_array:
        uv = np.asarray(uv, dtype=float)
        pos = transform_adsorbate(adsorbate, density, uv, float(z_height), quaternion, z_ref)
        positions.append(pos)
        pose_params.append(np.concatenate([uv, np.array([float(z_height)], dtype=float), quaternion]))
    return PoseBatch(
        adsorbate=adsorbate,
        positions=np.asarray(positions, dtype=float),
        pose_params=np.asarray(pose_params, dtype=float),
        site_labels=["scan"] * len(uv_array),
        metadata={"z_ref": z_ref, "z_height": float(z_height), "scan_type": "lateral"},
    )


def evaluate_teacher_chunked(teacher_backend, density, poses, chunk_size=32):
    n = len(poses.positions)
    if n <= chunk_size:
        t0 = time.perf_counter()
        result = teacher_backend.evaluate_batch(density, poses)
        elapsed = time.perf_counter() - t0
        return result, elapsed

    energies, forces = [], []
    t0 = time.perf_counter()
    for start in range(0, n, chunk_size):
        stop = min(start + chunk_size, n)
        idx = np.arange(start, stop, dtype=int)
        partial_poses = PoseBatch(
            adsorbate=poses.adsorbate,
            positions=np.asarray(poses.positions[idx], dtype=float),
            pose_params=np.asarray(poses.pose_params[idx], dtype=float),
            site_labels=[poses.site_labels[int(i)] for i in idx],
            metadata=dict(poses.metadata),
        )
        partial = teacher_backend.evaluate_batch(density, partial_poses)
        energies.append(np.asarray(partial.energies, dtype=float))
        forces.append(np.asarray(partial.forces, dtype=float))
        print(f"  MAD-SURF: {stop}/{n} ({100.0*stop/n:.0f}%)", flush=True)
    elapsed = time.perf_counter() - t0
    return TeacherResult(
        energies=np.concatenate(energies, axis=0),
        forces=np.concatenate(forces, axis=0),
        metadata={"teacher_eval_seconds": elapsed},
    ), elapsed


def evaluate_student_chunked(model, poses, params, chunk_size=64, compute_forces=True):
    """Evaluate GridFF student in chunks, with optional force computation."""
    n = len(poses.positions)
    all_energies = []
    all_forces = []
    t0 = time.perf_counter()
    for start in range(0, n, chunk_size):
        stop = min(start + chunk_size, n)
        idx = np.arange(start, stop, dtype=int)
        partial_poses = PoseBatch(
            adsorbate=poses.adsorbate,
            positions=np.asarray(poses.positions[idx], dtype=float),
            pose_params=np.asarray(poses.pose_params[idx], dtype=float),
            site_labels=[poses.site_labels[int(i)] for i in idx],
            metadata=dict(poses.metadata),
        )
        partial = model.evaluate_batch(partial_poses, params=params, compute_forces=compute_forces)
        all_energies.append(np.asarray(partial.energies, dtype=float))
        if compute_forces:
            all_forces.append(np.asarray(partial.forces, dtype=float))
    elapsed = time.perf_counter() - t0
    result_energies = np.concatenate(all_energies, axis=0)
    result_forces = np.concatenate(all_forces, axis=0) if compute_forces else None
    return type("StudentResult", (), {
        "energies": result_energies,
        "forces": result_forces,
    })(), elapsed


def force_norms(forces_array):
    """Compute |F| per pose from (n_poses, n_atoms, 3) array."""
    f = np.asarray(forces_array, dtype=float)
    return np.linalg.norm(f.reshape(f.shape[0], -1), axis=1)


# ═════════════════════════════════════════════════════════════════════
#  1D path:  top -> bridge -> hollow -> top
# ═════════════════════════════════════════════════════════════════════

def make_1d_path(site_uvs, n_per_segment=20):
    waypoints = [
        ("top",    site_uvs["top"]),
        ("bridge", site_uvs["bridge"]),
        ("hollow", site_uvs["hollow"]),
        ("top",    site_uvs["top"]),
    ]
    uv_list = []
    for i in range(len(waypoints) - 1):
        uv_start = waypoints[i][1]
        uv_end = waypoints[i + 1][1]
        duv = uv_end - uv_start
        duv = duv - np.round(duv)
        for t in np.linspace(0.0, 1.0, n_per_segment, endpoint=(i == len(waypoints) - 2)):
            uv_list.append((uv_start + t * duv).copy())
    return np.asarray(uv_list, dtype=float), [w[0] for w in waypoints]


def compute_path_distance(uv_array, cell):
    cell2 = cell[:2, :2]
    dists = [0.0]
    for i in range(1, len(uv_array)):
        duv = uv_array[i] - uv_array[i - 1]
        duv = duv - np.round(duv)
        dxy = duv @ cell2
        dists.append(dists[-1] + np.linalg.norm(dxy))
    return np.asarray(dists, dtype=float)


def find_segment_boundaries(uv_array, cell, n_per_segment=20):
    dist = compute_path_distance(uv_array, cell)
    n = n_per_segment
    boundaries = []
    for i in range(4):
        idx = min(i * n, len(dist) - 1) if i < 3 else len(dist) - 1
        boundaries.append(dist[idx])
    return boundaries


def _plot_1d_generic(dist, Y_madsurf, Y_gridff, boundaries, z_height,
                     ylabel, unit, title_suffix, out_path, adsorbate_name):
    """Generic 1D lateral scan plot for energy or force."""
    error = Y_gridff - Y_madsurf

    fig, ax1 = plt.subplots(figsize=(6.0, 3.8), dpi=300)

    ax1.plot(dist, Y_madsurf, color=CLR_MADSURF, lw=2.2, label="MAD-SURF", zorder=3)
    ax1.plot(dist, Y_gridff, color=CLR_GRIDFF, lw=1.8, ls=":", label="GridFF", zorder=3)

    ax1.set_xlabel(r"Path distance ($\AA$)", fontsize=11)
    ax1.set_ylabel(ylabel, fontsize=11)

    site_labels = ["Top", "Bridge", "Hollow", "Top"]
    for bd, lab in zip(boundaries, site_labels):
        ax1.axvline(bd, color="grey", ls="--", lw=0.7, alpha=0.5)
    ax1.set_xlim(dist[0], dist[-1])

    # Hide main right spine before twinx
    ax1.spines["right"].set_visible(False)

    # Error on twin axis
    ax2 = ax1.twinx()
    ax2.fill_between(dist, 0, np.abs(error), color=CLR_ERROR, alpha=0.15, zorder=1)
    ax2.plot(dist, error, color=CLR_ERROR, lw=1.0, alpha=0.6, label="Error", zorder=2)
    ax2.set_ylabel(f"Error ({unit})", fontsize=10, color=CLR_ERROR)
    ax2.tick_params(axis="y", colors=CLR_ERROR)
    ax2.spines["right"].set_color(CLR_ERROR)
    ax2.spines["right"].set_alpha(0.5)
    ax2.spines["left"].set_visible(False)
    ax2.spines["top"].set_visible(False)

    # Site labels at top
    ymin, ymax = ax1.get_ylim()
    for bd, lab in zip(boundaries, site_labels):
        ax1.text(bd, ymax + 0.02 * (ymax - ymin), lab,
                 ha="center", va="bottom", fontsize=9, fontweight="bold", color="dimgrey")

    # Metrics
    mae = np.mean(np.abs(error))
    rmse = float(np.sqrt(np.mean(error ** 2)))
    corrugation_ms = float(Y_madsurf.max() - Y_madsurf.min())
    corrugation_gf = float(Y_gridff.max() - Y_gridff.min())
    ax1.text(
        0.98, 0.04,
        f"MAE = {mae:.2f} {unit}\n"
        f"RMSE = {rmse:.2f} {unit}\n"
        f"corrugation: {corrugation_ms:.1f} (MAD-SURF) / {corrugation_gf:.1f} (GridFF)",
        transform=ax1.transAxes, fontsize=7,
        ha="right", va="bottom",
        bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.85, ec="lightgrey"),
    )

    ax1.legend(loc="upper left", fontsize=8.5, framealpha=0.9)
    ax1.set_title(f"{adsorbate_name} / Ag(111)  lateral {title_suffix}  (z = {z_height:.1f} " + r"$\AA$)",
                  fontsize=11, pad=16)
    ax1.grid(True, alpha=0.2)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[lateral] 1D plot saved: {out_path}", flush=True)


def plot_1d_energy(dist, E_ms, E_gf, boundaries, z_height, out_path, adsorbate_name):
    _plot_1d_generic(dist, E_ms * 1000, E_gf * 1000, boundaries, z_height,
                     "Interaction energy (meV)", "meV", "energy scan",
                     out_path, adsorbate_name)

def plot_1d_force(dist, F_ms, F_gf, boundaries, z_height, out_path, adsorbate_name):
    _plot_1d_generic(dist, F_ms, F_gf, boundaries, z_height,
                     r"Force norm (eV/$\AA$)", r"eV/$\AA$", "force scan",
                     out_path, adsorbate_name)


# ═════════════════════════════════════════════════════════════════════
#  2D grid scan
# ═════════════════════════════════════════════════════════════════════

def make_2d_grid(n_grid=25):
    u_vals = np.linspace(0.0, 1.0, n_grid, endpoint=False)
    v_vals = np.linspace(0.0, 1.0, n_grid, endpoint=False)
    uu, vv = np.meshgrid(u_vals, v_vals, indexing="ij")
    uv_flat = np.column_stack([uu.ravel(), vv.ravel()])
    return u_vals, v_vals, uv_flat


def _plot_2d_generic(u_vals, v_vals, Z_ms_2d, Z_gf_2d, site_uvs, z_height,
                     cbar_label, title_suffix, out_path, adsorbate_name):
    """Generic 2D heatmap: MAD-SURF, GridFF, Error."""
    error_2d = Z_gf_2d - Z_ms_2d

    vmin = min(Z_ms_2d.min(), Z_gf_2d.min())
    vmax = max(Z_ms_2d.max(), Z_gf_2d.max())
    err_abs_max = max(abs(error_2d.min()), abs(error_2d.max()))

    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.8), dpi=300)

    titles = ["MAD-SURF", "GridFF", "Error (GridFF - MAD-SURF)"]
    data = [Z_ms_2d, Z_gf_2d, error_2d]
    cmaps = ["viridis", "viridis", "RdBu_r"]
    vmins = [vmin, vmin, -err_abs_max]
    vmaxs = [vmax, vmax, err_abs_max]

    for i, ax in enumerate(axes):
        im = ax.pcolormesh(
            u_vals, v_vals, data[i].T,
            cmap=cmaps[i], vmin=vmins[i], vmax=vmaxs[i], shading="auto",
        )
        cb = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.03)
        cb.set_label(cbar_label, fontsize=9)
        ax.set_xlabel("u (frac.)", fontsize=10)
        if i == 0:
            ax.set_ylabel("v (frac.)", fontsize=10)
        ax.set_title(titles[i], fontsize=10)
        ax.set_aspect("equal")

        for label, uv in site_uvs.items():
            marker_color = "white" if i < 2 else "black"
            ax.plot(uv[0], uv[1], "o", color=marker_color, ms=6, mew=1.5, mfc="none", zorder=5)
            ax.text(uv[0] + 0.03, uv[1] + 0.03, label.capitalize(),
                    fontsize=7, color=marker_color, fontweight="bold", zorder=5)

    mae = np.mean(np.abs(error_2d))
    rmse = np.sqrt(np.mean(error_2d ** 2))
    fig.suptitle(
        f"{adsorbate_name} / Ag(111)  2D lateral {title_suffix}  "
        f"(z = {z_height:.1f} " + r"$\AA$" + f")    MAE = {mae:.2f},  RMSE = {rmse:.2f}",
        fontsize=10.5, y=1.02,
    )

    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[lateral] 2D plot saved: {out_path}", flush=True)


def plot_2d_energy(u_vals, v_vals, E_ms_2d, E_gf_2d, site_uvs, z_height, out_path, adsorbate_name):
    _plot_2d_generic(u_vals, v_vals, E_ms_2d * 1000, E_gf_2d * 1000, site_uvs, z_height,
                     "meV", "energy scan", out_path, adsorbate_name)

def plot_2d_force(u_vals, v_vals, F_ms_2d, F_gf_2d, site_uvs, z_height, out_path, adsorbate_name):
    _plot_2d_generic(u_vals, v_vals, F_ms_2d, F_gf_2d, site_uvs, z_height,
                     r"eV/$\AA$", "force scan", out_path, adsorbate_name)


# ═════════════════════════════════════════════════════════════════════
#  Main
# ═════════════════════════════════════════════════════════════════════

def parse_args():
    p = argparse.ArgumentParser(description="Lateral scan benchmark")
    p.add_argument("--result-dir", type=str, default="/tmp/ag_raw_r6_final",
                   help="Directory with config.json and fit_params.json from z-scan")
    p.add_argument("--adsorbate", type=str, default="CHONH2",
                   help="Adsorbate name (CHONH2, PTCDA, etc.)")
    p.add_argument("--xyz-path", type=str, default=None,
                   help="Path to adsorbate xyz file (needed for non-builtin molecules)")
    p.add_argument("--z-height", type=float, default=2.7,
                   help="Height above surface reference plane (A)")
    p.add_argument("--n-path", type=int, default=60,
                   help="Total points along 1D path (split across 3 segments)")
    p.add_argument("--n-grid", type=int, default=25,
                   help="Number of grid points per dimension for 2D scan")
    p.add_argument("--teacher-chunk", type=int, default=32,
                   help="MAD-SURF batch chunk size")
    p.add_argument("--skip-teacher-2d", action="store_true",
                   help="Skip 2D MAD-SURF eval (load from cache if available)")
    p.add_argument("--prefer-jax", dest="prefer_jax", action="store_true", default=True)
    p.add_argument("--no-prefer-jax", dest="prefer_jax", action="store_false")
    p.add_argument("--device", default="cuda", help="MAD-SURF MACE device")
    return p.parse_args()


def main():
    args = parse_args()

    result_dir = Path(args.result_dir)
    config_path = result_dir / "config.json"
    params_path = result_dir / "fit_params.json"
    out_dir = result_dir / "presentation_data"
    out_dir.mkdir(parents=True, exist_ok=True)
    adsorbate_name = args.adsorbate

    print("=" * 60)
    print(f"  {adsorbate_name} / Ag(111) lateral scan benchmark")
    print("=" * 60)

    # ── Load config and build infrastructure ──────────────────────
    config = load_config(config_path)
    config.teacher_backend.device = args.device

    print("[lateral] Loading density ...", flush=True)
    density = make_density_backend(
        config.density_backend, surface=config.surface, grid=config.grid
    ).load()

    print("[lateral] Loading adsorbate ...", flush=True)
    adsorbate = get_adsorbate(name=adsorbate_name, xyz_path=args.xyz_path)

    print("[lateral] Building model ...", flush=True)
    model = build_model(config, density, adsorbate, prefer_jax=args.prefer_jax)
    params = load_fitted_params(params_path)

    print("[lateral] Loading MAD-SURF backend ...", flush=True)
    teacher_backend = make_teacher_backend(config.teacher_backend)

    # ── Get high-symmetry sites ──────────────────────────────────
    site_uvs = _site_uv_lookup(density, ["top", "bridge", "hollow"])
    print(f"[lateral] Sites: { {k: v.tolist() for k, v in site_uvs.items()} }", flush=True)

    z_height = float(args.z_height)
    quaternion = np.array([1.0, 0.0, 0.0, 0.0], dtype=float)
    print(f"[lateral] z_height = {z_height:.2f} A above surface", flush=True)

    cell = np.asarray(density.cell, dtype=float)

    # ==============================================================
    #  1D path scan:  top -> bridge -> hollow -> top
    # ==============================================================
    print("\n--- 1D path scan ---", flush=True)
    n_per_seg = max(args.n_path // 3, 5)
    uv_path, _ = make_1d_path(site_uvs, n_per_segment=n_per_seg)
    print(f"[lateral] 1D path: {len(uv_path)} points, {n_per_seg} per segment", flush=True)

    poses_1d = build_lateral_poses(density, adsorbate, uv_path, z_height, quaternion)

    # MAD-SURF
    print("[lateral] Evaluating MAD-SURF on 1D path ...", flush=True)
    teacher_1d, t_teacher_1d = evaluate_teacher_chunked(
        teacher_backend, density, poses_1d, chunk_size=args.teacher_chunk
    )
    print(f"[lateral] MAD-SURF 1D: {t_teacher_1d:.1f}s ({t_teacher_1d/len(uv_path):.2f}s/pose)", flush=True)

    # GridFF (with forces)
    print("[lateral] Evaluating GridFF on 1D path ...", flush=True)
    student_1d, t_student_1d = evaluate_student_chunked(
        model, poses_1d, params, chunk_size=64, compute_forces=True
    )
    print(f"[lateral] GridFF 1D: {t_student_1d:.3f}s", flush=True)

    # Distances along path
    path_dist = compute_path_distance(uv_path, cell)
    boundaries = find_segment_boundaries(uv_path, cell, n_per_segment=n_per_seg)

    E_ms_1d = np.asarray(teacher_1d.energies, dtype=float)
    E_gf_1d = np.asarray(student_1d.energies, dtype=float)

    # Force norms
    F_ms_1d = force_norms(teacher_1d.forces)
    F_gf_1d = force_norms(student_1d.forces)

    # Cartesian path coordinates for CSV
    cart_path = np.array([uv @ cell[:2, :2] for uv in uv_path])

    # Save 1D data
    np.savez(
        out_dir / "lateral_1d_path.npz",
        uv=uv_path,
        path_distance=path_dist,
        madsurf_energy=E_ms_1d,
        gridff_energy=E_gf_1d,
        madsurf_force_norm=F_ms_1d,
        gridff_force_norm=F_gf_1d,
        madsurf_forces=np.asarray(teacher_1d.forces, dtype=float),
        gridff_forces=np.asarray(student_1d.forces, dtype=float),
        z_height=z_height,
        boundaries=boundaries,
    )
    # CSV
    with open(out_dir / "lateral_1d_path.csv", "w") as fh:
        fh.write("path_dist_A,x_A,y_A,u,v,madsurf_energy_eV,gridff_energy_eV,energy_error_eV,"
                 "madsurf_force_norm_eV_A,gridff_force_norm_eV_A,force_error_eV_A\n")
        for i in range(len(uv_path)):
            fh.write(f"{path_dist[i]:.6f},{cart_path[i,0]:.6f},{cart_path[i,1]:.6f},"
                     f"{uv_path[i,0]:.6f},{uv_path[i,1]:.6f},"
                     f"{E_ms_1d[i]:.8f},{E_gf_1d[i]:.8f},{E_gf_1d[i]-E_ms_1d[i]:.8f},"
                     f"{F_ms_1d[i]:.8f},{F_gf_1d[i]:.8f},{F_gf_1d[i]-F_ms_1d[i]:.8f}\n")
    print(f"[lateral] 1D data saved", flush=True)

    # 1D Metrics
    err_E_1d = (E_gf_1d - E_ms_1d) * 1000.0
    err_F_1d = F_gf_1d - F_ms_1d
    print(f"[lateral] 1D energy  MAE = {np.mean(np.abs(err_E_1d)):.1f} meV  "
          f"RMSE = {np.sqrt(np.mean(err_E_1d**2)):.1f} meV")
    print(f"[lateral] 1D force   MAE = {np.mean(np.abs(err_F_1d)):.3f} eV/A  "
          f"RMSE = {np.sqrt(np.mean(err_F_1d**2)):.3f} eV/A")
    print(f"[lateral] 1D corrugation: MAD-SURF = {(E_ms_1d.max()-E_ms_1d.min())*1000:.0f} meV  "
          f"GridFF = {(E_gf_1d.max()-E_gf_1d.min())*1000:.0f} meV")

    # Plot 1D energy
    plot_1d_energy(path_dist, E_ms_1d, E_gf_1d, boundaries, z_height,
                   out_dir / "lateral_1d_energy.pdf", adsorbate_name)
    # Plot 1D force
    plot_1d_force(path_dist, F_ms_1d, F_gf_1d, boundaries, z_height,
                  out_dir / "lateral_1d_force.pdf", adsorbate_name)

    # ==============================================================
    #  2D grid scan
    # ==============================================================
    print("\n--- 2D grid scan ---", flush=True)
    n_grid = int(args.n_grid)
    u_vals, v_vals, uv_flat = make_2d_grid(n_grid)
    n_2d = len(uv_flat)
    print(f"[lateral] 2D grid: {n_grid}x{n_grid} = {n_2d} points", flush=True)

    cache_2d = out_dir / "lateral_2d_grid.npz"

    # MAD-SURF 2D
    if args.skip_teacher_2d and cache_2d.exists():
        print("[lateral] Loading cached 2D MAD-SURF data ...", flush=True)
        cached = np.load(cache_2d)
        E_ms_2d_flat = cached["madsurf_energy"] if "madsurf_energy" in cached else cached["teacher_energy"]
        F_ms_2d_flat = cached["madsurf_force_norm"] if "madsurf_force_norm" in cached else None
        if F_ms_2d_flat is None:
            print("[lateral] Warning: no cached force data, re-running MAD-SURF 2D", flush=True)
            args.skip_teacher_2d = False  # force re-evaluation

    if not args.skip_teacher_2d or not cache_2d.exists():
        poses_2d = build_lateral_poses(density, adsorbate, uv_flat, z_height, quaternion)
        print(f"[lateral] Evaluating MAD-SURF on 2D grid ({n_2d} poses) ...", flush=True)
        teacher_2d, t_teacher_2d = evaluate_teacher_chunked(
            teacher_backend, density, poses_2d, chunk_size=args.teacher_chunk
        )
        E_ms_2d_flat = np.asarray(teacher_2d.energies, dtype=float)
        F_ms_2d_flat = force_norms(teacher_2d.forces)
        print(f"[lateral] MAD-SURF 2D: {t_teacher_2d:.1f}s ({t_teacher_2d/n_2d:.2f}s/pose)", flush=True)

    # GridFF 2D (with forces)
    poses_2d = build_lateral_poses(density, adsorbate, uv_flat, z_height, quaternion)
    print("[lateral] Evaluating GridFF on 2D grid ...", flush=True)
    student_2d, t_student_2d = evaluate_student_chunked(
        model, poses_2d, params, chunk_size=64, compute_forces=True
    )
    E_gf_2d_flat = np.asarray(student_2d.energies, dtype=float)
    F_gf_2d_flat = force_norms(student_2d.forces)
    print(f"[lateral] GridFF 2D: {t_student_2d:.3f}s", flush=True)

    # Reshape to 2D
    E_ms_2d = E_ms_2d_flat.reshape(n_grid, n_grid)
    E_gf_2d = E_gf_2d_flat.reshape(n_grid, n_grid)
    F_ms_2d = F_ms_2d_flat.reshape(n_grid, n_grid)
    F_gf_2d = F_gf_2d_flat.reshape(n_grid, n_grid)

    # Cartesian grid for CSV
    cart_flat = np.array([uv @ cell[:2, :2] for uv in uv_flat])

    # Save 2D data
    np.savez(
        cache_2d,
        u_vals=u_vals, v_vals=v_vals,
        madsurf_energy=E_ms_2d_flat,
        gridff_energy=E_gf_2d_flat,
        madsurf_force_norm=F_ms_2d_flat,
        gridff_force_norm=F_gf_2d_flat,
        madsurf_energy_2d=E_ms_2d,
        gridff_energy_2d=E_gf_2d,
        madsurf_force_2d=F_ms_2d,
        gridff_force_2d=F_gf_2d,
        z_height=z_height,
    )
    # CSV
    with open(out_dir / "lateral_2d_grid.csv", "w") as fh:
        fh.write("x_A,y_A,u,v,madsurf_energy_eV,gridff_energy_eV,energy_error_eV,"
                 "madsurf_force_norm_eV_A,gridff_force_norm_eV_A,force_error_eV_A\n")
        for i in range(n_2d):
            fh.write(f"{cart_flat[i,0]:.6f},{cart_flat[i,1]:.6f},"
                     f"{uv_flat[i,0]:.6f},{uv_flat[i,1]:.6f},"
                     f"{E_ms_2d_flat[i]:.8f},{E_gf_2d_flat[i]:.8f},"
                     f"{E_gf_2d_flat[i]-E_ms_2d_flat[i]:.8f},"
                     f"{F_ms_2d_flat[i]:.8f},{F_gf_2d_flat[i]:.8f},"
                     f"{F_gf_2d_flat[i]-F_ms_2d_flat[i]:.8f}\n")
    print(f"[lateral] 2D data saved", flush=True)

    # 2D Metrics
    err_E_2d = (E_gf_2d_flat - E_ms_2d_flat) * 1000.0
    err_F_2d = F_gf_2d_flat - F_ms_2d_flat
    print(f"[lateral] 2D energy  MAE = {np.mean(np.abs(err_E_2d)):.1f} meV  "
          f"RMSE = {np.sqrt(np.mean(err_E_2d**2)):.1f} meV")
    print(f"[lateral] 2D force   MAE = {np.mean(np.abs(err_F_2d)):.3f} eV/A  "
          f"RMSE = {np.sqrt(np.mean(err_F_2d**2)):.3f} eV/A")
    print(f"[lateral] 2D corrugation: MAD-SURF = {(E_ms_2d_flat.max()-E_ms_2d_flat.min())*1000:.0f} meV  "
          f"GridFF = {(E_gf_2d_flat.max()-E_gf_2d_flat.min())*1000:.0f} meV")

    # Plot 2D energy
    plot_2d_energy(u_vals, v_vals, E_ms_2d, E_gf_2d, site_uvs, z_height,
                   out_dir / "lateral_2d_energy_heatmaps.pdf", adsorbate_name)
    # Plot 2D force
    plot_2d_force(u_vals, v_vals, F_ms_2d, F_gf_2d, site_uvs, z_height,
                  out_dir / "lateral_2d_force_heatmaps.pdf", adsorbate_name)

    # ── Summary ──────────────────────────────────────────────────
    summary = {
        "adsorbate": adsorbate_name,
        "z_height_A": z_height,
        "scan_1d": {
            "n_points": len(uv_path),
            "energy_mae_meV": float(np.mean(np.abs(err_E_1d))),
            "energy_rmse_meV": float(np.sqrt(np.mean(err_E_1d ** 2))),
            "force_mae_eV_A": float(np.mean(np.abs(err_F_1d))),
            "force_rmse_eV_A": float(np.sqrt(np.mean(err_F_1d ** 2))),
            "corrugation_madsurf_meV": float((E_ms_1d.max() - E_ms_1d.min()) * 1000),
            "corrugation_gridff_meV": float((E_gf_1d.max() - E_gf_1d.min()) * 1000),
        },
        "scan_2d": {
            "n_grid": n_grid,
            "n_points": n_2d,
            "energy_mae_meV": float(np.mean(np.abs(err_E_2d))),
            "energy_rmse_meV": float(np.sqrt(np.mean(err_E_2d ** 2))),
            "force_mae_eV_A": float(np.mean(np.abs(err_F_2d))),
            "force_rmse_eV_A": float(np.sqrt(np.mean(err_F_2d ** 2))),
            "corrugation_madsurf_meV": float((E_ms_2d_flat.max() - E_ms_2d_flat.min()) * 1000),
            "corrugation_gridff_meV": float((E_gf_2d_flat.max() - E_gf_2d_flat.min()) * 1000),
        },
    }
    with open(out_dir / "lateral_scan_metrics.json", "w") as fh:
        json.dump(summary, fh, indent=2)

    print(f"\n{'='*60}")
    print(f"  {adsorbate_name} lateral scan complete")
    print(f"  1D energy MAE = {summary['scan_1d']['energy_mae_meV']:.1f} meV   "
          f"2D energy MAE = {summary['scan_2d']['energy_mae_meV']:.1f} meV")
    print(f"  1D force  MAE = {summary['scan_1d']['force_mae_eV_A']:.3f} eV/A  "
          f"2D force  MAE = {summary['scan_2d']['force_mae_eV_A']:.3f} eV/A")
    print("=" * 60)


if __name__ == "__main__":
    main()
