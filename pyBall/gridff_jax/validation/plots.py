"""Validation plots matching the current GridFF diagnostics style."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from ..utils import ensure_dir


def _figure_path(out_dir, name):
    out_dir = ensure_dir(out_dir)
    return Path(out_dir) / name


def plot_parity(y_true, y_pred, out_dir, name):
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    lo = min(y_true.min(), y_pred.min())
    hi = max(y_true.max(), y_pred.max())
    plt.figure(figsize=(5, 5))
    plt.scatter(y_true, y_pred, s=14, alpha=0.55, color="tab:blue", edgecolors="none")
    plt.plot([lo, hi], [lo, hi], "--k", lw=1)
    plt.xlabel("Teacher")
    plt.ylabel("Student")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(_figure_path(out_dir, name), dpi=160)
    plt.close()


def plot_error_histogram(errors, out_dir, name, log_scale: bool = False):
    errors = np.asarray(errors, dtype=float).reshape(-1)
    plt.figure(figsize=(6, 4))
    plt.hist(np.abs(errors), bins=40, alpha=0.8)
    if log_scale:
        plt.yscale("log")
    plt.xlabel("|Error|")
    plt.ylabel("Count")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(_figure_path(out_dir, name), dpi=160)
    plt.close()


def plot_convergence(history, out_dir, name_prefix="loss"):
    if not history:
        return
    steps = np.arange(len(history))
    train = np.array([item["train_loss"] for item in history], dtype=float)
    val = np.array([item["val_loss"] for item in history], dtype=float)
    for scale, suffix in [("linear", "linear.png"), ("log", "log.png")]:
        plt.figure(figsize=(6, 4))
        plt.plot(steps, train, label="train")
        plt.plot(steps, val, label="val")
        if scale == "log":
            plt.yscale("log")
        plt.xlabel("Iteration")
        plt.ylabel("Loss")
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.tight_layout()
        plt.savefig(_figure_path(out_dir, f"{name_prefix}_{suffix}"), dpi=160)
        plt.close()


def plot_z_profile(z_values, teacher_energy, pred_energy, out_dir, name):
    z_values = np.asarray(z_values, dtype=float)
    teacher_energy = np.asarray(teacher_energy, dtype=float)
    pred_energy = np.asarray(pred_energy, dtype=float)
    order = np.argsort(z_values)
    z_sorted = z_values[order]
    teacher_sorted = teacher_energy[order]
    pred_sorted = pred_energy[order]
    nbins = int(min(24, max(6, np.sqrt(len(z_sorted)))))
    edges = np.linspace(z_sorted.min(), z_sorted.max(), nbins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    teacher_mean = np.full(nbins, np.nan, dtype=float)
    pred_mean = np.full(nbins, np.nan, dtype=float)
    for ibin, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
        if ibin == (nbins - 1):
            mask = (z_sorted >= lo) & (z_sorted <= hi)
        else:
            mask = (z_sorted >= lo) & (z_sorted < hi)
        if np.count_nonzero(mask) == 0:
            continue
        teacher_mean[ibin] = np.mean(teacher_sorted[mask])
        pred_mean[ibin] = np.mean(pred_sorted[mask])
    plt.figure(figsize=(6, 4))
    plt.scatter(z_sorted, teacher_sorted, s=10, alpha=0.25, color="black", label="teacher samples")
    plt.scatter(z_sorted, pred_sorted, s=10, alpha=0.25, color="tab:red", label="student samples")
    valid = np.isfinite(teacher_mean) & np.isfinite(pred_mean)
    if np.any(valid):
        plt.plot(centers[valid], teacher_mean[valid], "-k", lw=2, label="teacher mean")
        plt.plot(centers[valid], pred_mean[valid], "--r", lw=2, label="student mean")
    plt.xlabel("z [A]")
    plt.ylabel("Energy [eV]")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(_figure_path(out_dir, name), dpi=160)
    plt.close()


def plot_grid_slices(grids, out_dir, prefix="grid"):
    scalar_grids = {
        "pauli_zyx": np.asarray(grids["pauli_zyx"], dtype=float),
        "london_zyx": np.asarray(grids["london_zyx"], dtype=float),
        "coulomb_zyx": np.asarray(grids["coulomb_zyx"], dtype=float),
    }
    reactive = np.asarray(grids["reactive_zyxc"], dtype=float)
    for channel in range(reactive.shape[-1]):
        scalar_grids[f"reactive_zyx_c{channel}"] = reactive[..., channel]
    for key, values in scalar_grids.items():
        iz = values.shape[0] // 2
        iy = values.shape[1] // 2
        ix = values.shape[2] // 2
        fig, axes = plt.subplots(1, 3, figsize=(14, 4))
        axes[0].imshow(values[iz], origin="lower", aspect="auto")
        axes[0].set_title(f"{key} XY")
        axes[1].imshow(values[:, iy, :], origin="lower", aspect="auto")
        axes[1].set_title(f"{key} XZ")
        axes[2].imshow(values[:, :, ix], origin="lower", aspect="auto")
        axes[2].set_title(f"{key} YZ")
        for ax in axes:
            ax.set_xlabel("grid")
            ax.set_ylabel("grid")
        plt.tight_layout()
        plt.savefig(_figure_path(out_dir, f"{prefix}_{key}.png"), dpi=160)
        plt.close(fig)
