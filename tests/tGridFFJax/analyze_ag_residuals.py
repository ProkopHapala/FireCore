#!/usr/bin/env python3
"""Residual analysis for strict-PLQ Ag GridFF fits against MAD-SURF labels."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.append(str(ROOT))

from pyBall.gridff_jax import load_config
from pyBall.gridff_jax.config import FeatureToggles
from pyBall.gridff_jax.density_backends import make_density_backend
from pyBall.gridff_jax.fit import load_pose_dataset
from pyBall.gridff_jax.hybrid_energy import HybridGridFFModel, HybridParameters
from pyBall.gridff_jax.utils import ensure_dir, save_json
from pyBall.gridff_jax.validation import compute_error_metrics


def parse_args():
    parser = argparse.ArgumentParser(description="Analyze residuals of a fitted Ag GridFF model against teacher labels")
    parser.add_argument("--config", default=str(ROOT / "tests/tGridFFJax/ag_fit/config_used.json"))
    parser.add_argument("--dataset-dir", default=str(ROOT / "tests/tGridFFJax/ag_dataset/datasets"))
    parser.add_argument("--fit-dir", default=str(ROOT / "tests/tGridFFJax/ag_fit"))
    parser.add_argument("--out-dir", default=str(ROOT / "tests/tGridFFJax/ag_fit/residual_analysis"))
    parser.add_argument("--mode", default="", choices=["", "plq", "plq_reactive", "plq_ct", "plq_ct_image", "full"])
    parser.add_argument("--adsorbates", default="", help="Optional comma-separated adsorbate filter")
    parser.add_argument("--path-benchmark-roots", default="", help="Optional comma-separated directories containing <adsorbate>/path_benchmark.json")
    parser.add_argument("--prefer-jax", dest="prefer_jax", action="store_true")
    parser.add_argument("--no-prefer-jax", dest="prefer_jax", action="store_false")
    parser.set_defaults(prefer_jax=True)
    return parser.parse_args()


def _load_params(path):
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


def _mode_toggles(mode: str):
    if mode == "plq":
        return FeatureToggles(True, True, False, False, False, True)
    if mode == "plq_reactive":
        return FeatureToggles(True, True, False, False, True, True)
    if mode == "plq_ct":
        return FeatureToggles(True, True, True, False, False, True)
    if mode == "plq_ct_image":
        return FeatureToggles(True, True, True, True, False, True)
    return FeatureToggles(True, True, True, True, True, True)


def _component_scalar(values):
    array = np.asarray(values, dtype=float)
    if array.ndim <= 1:
        return array.reshape(-1)
    return array.reshape(array.shape[0], -1).sum(axis=1)


def _safe_corr(x, y):
    x = np.asarray(x, dtype=float).reshape(-1)
    y = np.asarray(y, dtype=float).reshape(-1)
    if x.size == 0 or y.size == 0:
        return 0.0
    if np.std(x) < 1.0e-12 or np.std(y) < 1.0e-12:
        return 0.0
    return float(np.corrcoef(x, y)[0, 1])


def _z_edges(z_values):
    z_values = np.asarray(z_values, dtype=float)
    zmin = float(np.floor(z_values.min() * 2.0) / 2.0)
    zmax = float(np.ceil(z_values.max() * 2.0) / 2.0)
    edges = np.arange(zmin, zmax + 0.51, 0.5, dtype=float)
    if edges.size < 3:
        edges = np.array([zmin, 0.5 * (zmin + zmax), zmax], dtype=float)
    return edges


def _site_color_map(labels):
    unique = sorted(set(labels))
    colors = plt.cm.tab10(np.linspace(0.0, 1.0, max(len(unique), 1)))
    return {label: colors[index] for index, label in enumerate(unique)}


def _plot_residual_vs_z(z_values, residuals, site_labels, ylabel, out_path):
    z_values = np.asarray(z_values, dtype=float)
    residuals = np.asarray(residuals, dtype=float)
    color_map = _site_color_map(site_labels)
    plt.figure(figsize=(6.5, 4.5))
    for label in sorted(color_map.keys()):
        mask = np.array([item == label for item in site_labels], dtype=bool)
        if not np.any(mask):
            continue
        plt.scatter(
            z_values[mask],
            residuals[mask],
            s=16,
            alpha=0.45,
            color=color_map[label],
            label=label,
        )
    edges = _z_edges(z_values)
    centers = 0.5 * (edges[:-1] + edges[1:])
    means = np.full(edges.size - 1, np.nan, dtype=float)
    for iedge, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
        mask = (z_values >= lo) & (z_values < hi if iedge < (edges.size - 2) else z_values <= hi)
        if np.any(mask):
            means[iedge] = np.mean(residuals[mask])
    valid = np.isfinite(means)
    if np.any(valid):
        plt.plot(centers[valid], means[valid], color="black", lw=2.4, label="binned mean", zorder=5)
    plt.axhline(0.0, color="black", lw=1.0, ls=":")
    plt.xlabel("z [A]")
    plt.ylabel(ylabel)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path, dpi=160)
    plt.close()


def _plot_site_z_heatmap(z_values, site_labels, values, out_path, title, cmap="coolwarm"):
    z_values = np.asarray(z_values, dtype=float)
    values = np.asarray(values, dtype=float)
    sites = sorted(set(site_labels))
    edges = _z_edges(z_values)
    grid = np.full((len(sites), edges.size - 1), np.nan, dtype=float)
    for isite, site in enumerate(sites):
        site_mask = np.array([label == site for label in site_labels], dtype=bool)
        for ibin, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
            z_mask = (z_values >= lo) & (z_values < hi if ibin < (edges.size - 2) else z_values <= hi)
            mask = site_mask & z_mask
            if np.any(mask):
                grid[isite, ibin] = np.mean(values[mask])
    plt.figure(figsize=(8, 3.6))
    plt.imshow(grid, origin="lower", aspect="auto", cmap=cmap)
    plt.colorbar()
    plt.yticks(np.arange(len(sites)), sites)
    plt.xticks(np.arange(edges.size - 1), [f"{lo:.1f}-{hi:.1f}" for lo, hi in zip(edges[:-1], edges[1:])], rotation=45, ha="right")
    plt.title(title)
    plt.xlabel("z bin [A]")
    plt.ylabel("site")
    plt.tight_layout()
    plt.savefig(out_path, dpi=160)
    plt.close()


def _plot_component_scatter(component, residuals, out_path, xlabel):
    component = np.asarray(component, dtype=float)
    residuals = np.asarray(residuals, dtype=float)
    plt.figure(figsize=(5.4, 4.4))
    plt.scatter(component, residuals, s=16, alpha=0.45)
    plt.axhline(0.0, color="black", lw=1.0, ls=":")
    plt.xlabel(xlabel)
    plt.ylabel("Energy residual [GridFF - MAD-SURF] [eV]")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_path, dpi=160)
    plt.close()


def _top_outliers(z_values, site_labels, teacher_energy, pred_energy, force_error_norm, components, pose_params, count=12):
    residual = pred_energy - teacher_energy
    order_e = np.argsort(np.abs(residual))[::-1][:count]
    order_f = np.argsort(force_error_norm)[::-1][:count]

    def pack(indices):
        payload = []
        for index in indices:
            payload.append(
                {
                    "index": int(index),
                    "site": str(site_labels[index]),
                    "z": float(z_values[index]),
                    "teacher_energy": float(teacher_energy[index]),
                    "pred_energy": float(pred_energy[index]),
                    "energy_residual": float(residual[index]),
                    "force_error_norm": float(force_error_norm[index]),
                    "components": {key: float(values[index]) for key, values in components.items()},
                    "pose_params": [float(v) for v in pose_params[index]],
                }
            )
        return payload

    return {
        "energy": pack(order_e),
        "force": pack(order_f),
    }


def _path_summaries(adsorbate, roots):
    summaries = []
    for root in roots:
        candidate = Path(root) / adsorbate / "path_benchmark.json"
        if candidate.is_file():
            with candidate.open("r", encoding="utf8") as fin:
                payload = json.load(fin)
            worst_path_energy = max(payload["paths"].items(), key=lambda item: item[1]["energy_error"]["rmse"])
            worst_path_force = max(payload["paths"].items(), key=lambda item: item[1]["force_error"]["rmse"])
            worst_plane_energy = max(payload["planes"].items(), key=lambda item: item[1]["energy_error"]["rmse"])
            worst_plane_force = max(payload["planes"].items(), key=lambda item: item[1]["force_error"]["rmse"])
            summaries.append(
                {
                    "root": str(root),
                    "path_file": str(candidate),
                    "worst_path_energy": {
                        "name": worst_path_energy[0],
                        "rmse": float(worst_path_energy[1]["energy_error"]["rmse"]),
                    },
                    "worst_path_force": {
                        "name": worst_path_force[0],
                        "rmse": float(worst_path_force[1]["force_error"]["rmse"]),
                    },
                    "worst_plane_energy": {
                        "name": worst_plane_energy[0],
                        "rmse": float(worst_plane_energy[1]["energy_error"]["rmse"]),
                    },
                    "worst_plane_force": {
                        "name": worst_plane_force[0],
                        "rmse": float(worst_plane_force[1]["force_error"]["rmse"]),
                    },
                }
            )
    return summaries


def main():
    args = parse_args()
    out_dir = ensure_dir(args.out_dir)
    config = load_config(args.config)
    if args.mode:
        config.toggles = _mode_toggles(args.mode)
        config.hybrid_model.use_qeq = config.toggles.use_ct_qeq
        config.hybrid_model.use_image = config.toggles.use_image_charge
        config.hybrid_model.use_reactive_grid = config.toggles.use_reactive_grid
    params = _load_params(Path(args.fit_dir) / "fit_params.json")
    density = make_density_backend(config.density_backend, surface=config.surface, grid=config.grid).load()

    selected_adsorbates = {name.strip() for name in args.adsorbates.split(",") if name.strip()}
    path_roots = [item.strip() for item in args.path_benchmark_roots.split(",") if item.strip()]
    dataset_dir = Path(args.dataset_dir)
    datasets = [load_pose_dataset(path) for path in sorted(dataset_dir.glob("*.npz")) if not selected_adsorbates or path.stem in selected_adsorbates]
    if not datasets:
        raise ValueError("No datasets selected for residual analysis")

    reactive_elements = tuple(sorted({symbol for poses, _, _ in datasets for symbol in poses.adsorbate.symbols}))
    model = HybridGridFFModel(
        density=density,
        reactive_elements=reactive_elements,
        toggles=config.toggles,
        grid_config=config.grid,
        hybrid_config=config.hybrid_model,
        prefer_jax=args.prefer_jax,
    )

    summary = {
        "mode": args.mode or "fit_config",
        "adsorbates_selected": sorted(selected_adsorbates) if selected_adsorbates else "all",
        "adsorbates": {},
    }

    for poses, teacher, _ in datasets:
        ads_dir = ensure_dir(out_dir / poses.adsorbate.name)
        pred = model.evaluate_batch(poses, params=params)
        z_values = np.asarray(poses.pose_params[:, 2], dtype=float)
        site_labels = list(poses.site_labels)
        teacher_energy = np.asarray(teacher.energies, dtype=float)
        pred_energy = np.asarray(pred.energies, dtype=float)
        energy_residual = pred_energy - teacher_energy
        force_error = np.asarray(pred.forces - teacher.forces, dtype=float)
        force_error_norm = np.linalg.norm(force_error.reshape(len(force_error), -1), axis=1)
        force_teacher_norm = np.linalg.norm(np.asarray(teacher.forces, dtype=float).reshape(len(force_error), -1), axis=1)
        components = {
            key: _component_scalar(values)
            for key, values in pred.components.items()
            if key in ("pauli", "london", "coulomb", "qeq", "image", "reactive")
        }

        np.savez_compressed(
            ads_dir / "residual_data.npz",
            z_values=z_values,
            teacher_energy=teacher_energy,
            pred_energy=pred_energy,
            energy_residual=energy_residual,
            teacher_forces=np.asarray(teacher.forces, dtype=float),
            pred_forces=np.asarray(pred.forces, dtype=float),
            force_error=force_error,
            force_error_norm=force_error_norm,
            teacher_force_norm=force_teacher_norm,
            pose_params=np.asarray(poses.pose_params, dtype=float),
            site_labels=np.asarray(site_labels),
            **{f"component_{key}": values for key, values in components.items()},
        )

        _plot_residual_vs_z(z_values, energy_residual, site_labels, "Energy residual [eV]", ads_dir / "energy_residual_vs_z.png")
        _plot_residual_vs_z(z_values, force_error_norm, site_labels, "Force error norm [eV/A]", ads_dir / "force_error_norm_vs_z.png")
        _plot_site_z_heatmap(z_values, site_labels, energy_residual, ads_dir / "energy_residual_site_z.png", "Mean energy residual by site/z", cmap="coolwarm")
        _plot_site_z_heatmap(z_values, site_labels, force_error_norm, ads_dir / "force_error_site_z.png", "Mean force error norm by site/z", cmap="viridis")
        for key, values in components.items():
            _plot_component_scatter(values, energy_residual, ads_dir / f"energy_residual_vs_{key}.png", f"{key} component [eV]")

        per_site = {}
        for site in sorted(set(site_labels)):
            mask = np.array([label == site for label in site_labels], dtype=bool)
            per_site[site] = {
                "n_samples": int(np.count_nonzero(mask)),
                "energy": compute_error_metrics(np.zeros(np.count_nonzero(mask), dtype=float), energy_residual[mask]),
                "force_norm": compute_error_metrics(np.zeros(np.count_nonzero(mask), dtype=float), force_error_norm[mask]),
            }

        component_correlations = {
            key: {
                "energy_residual": _safe_corr(values, energy_residual),
                "force_error_norm": _safe_corr(values, force_error_norm),
            }
            for key, values in components.items()
        }

        ads_summary = {
            "n_samples": int(len(poses.positions)),
            "energy": compute_error_metrics(np.zeros_like(energy_residual), energy_residual),
            "force_norm": compute_error_metrics(np.zeros_like(force_error_norm), force_error_norm),
            "force_components": compute_error_metrics(np.zeros(force_error.size, dtype=float), force_error.reshape(-1)),
            "per_site": per_site,
            "component_correlations": component_correlations,
            "diagnostics": {
                "energy_residual_vs_z_corr": _safe_corr(z_values, energy_residual),
                "force_error_norm_vs_z_corr": _safe_corr(z_values, force_error_norm),
                "energy_residual_vs_teacher_energy_corr": _safe_corr(teacher_energy, energy_residual),
                "force_error_norm_vs_teacher_force_norm_corr": _safe_corr(force_teacher_norm, force_error_norm),
            },
            "outliers": _top_outliers(
                z_values=z_values,
                site_labels=site_labels,
                teacher_energy=teacher_energy,
                pred_energy=pred_energy,
                force_error_norm=force_error_norm,
                components=components,
                pose_params=np.asarray(poses.pose_params, dtype=float),
            ),
            "path_benchmark": _path_summaries(poses.adsorbate.name, path_roots),
        }
        summary["adsorbates"][poses.adsorbate.name] = ads_summary
        save_json(ads_summary, ads_dir / "residual_summary.json")

    save_json(summary, out_dir / "residual_analysis.json")
    print(f"residual analysis -> {out_dir}")


if __name__ == "__main__":
    main()
