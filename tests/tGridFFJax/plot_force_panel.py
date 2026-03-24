#!/usr/bin/env python3
"""Publication-quality 1x3 force z-scan panel for CHONH2/Ag(111) GridFF benchmark.

Each panel shows |Force| magnitude vs z for MAD-SURF (teacher) and GridFF (student),
plus force error on a right y-axis.

Style matches plot_publication_panel.py (energy panel).

Usage:
    python plot_force_panel.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------- configuration ----------
SITES = ["top", "bridge", "hollow"]
SITE_LABELS = {"top": "Top", "bridge": "Bridge", "hollow": "Hollow"}
TRACE_DIR = Path(__file__).parent.parent / "traces"
OUT_DIR = Path(__file__).parent
Z_MIN = 2.0
Z_MAX = 5.6

# Colours — same as energy panel
CLR_TEACHER = "orange"
CLR_STUDENT = "#1f77b4"
CLR_ERROR   = "#cc3333"


# ---------- helpers ----------
def load_site_force(site: str):
    """Load trace and compute per-pose force magnitude."""
    d = np.load(TRACE_DIR / f"{site}_zscan_trace.npz")
    z = d["z_values"]
    nz = len(z)

    # teacher_force / student_force: (nz, n_atoms, 3)
    F_teacher = d["teacher_force"]
    F_student = d["student_force"]

    # |F| per pose = norm of the full 3*n_atoms force vector
    F_teacher_mag = np.linalg.norm(F_teacher.reshape(nz, -1), axis=1)
    F_student_mag = np.linalg.norm(F_student.reshape(nz, -1), axis=1)

    mask = z >= Z_MIN
    return {
        "z": z[mask],
        "F_teacher": F_teacher_mag[mask],
        "F_student": F_student_mag[mask],
        "train": d["train_mask"][mask],
    }


def export_csv(site: str, data: dict):
    """Export force z-scan to CSV."""
    csv_path = OUT_DIR / f"{site}_force_zscan.csv"
    header = "z_A,F_teacher_eV_A,F_student_eV_A,error_eV_A"
    err = data["F_student"] - data["F_teacher"]
    arr = np.column_stack([data["z"], data["F_teacher"], data["F_student"], err])
    np.savetxt(csv_path, arr, header=header, delimiter=",",
               fmt="%.6f", comments="")
    print(f"Saved CSV: {csv_path}")


def make_figure():
    # Pre-compute global limits
    all_F = []
    all_err = []
    site_data = {}
    for site in SITES:
        data = load_site_force(site)
        site_data[site] = data
        all_F.extend(data["F_teacher"])
        all_F.extend(data["F_student"])
        all_err.extend(data["F_student"] - data["F_teacher"])

    F_hi = max(all_F) * 1.10
    err_lim = max(abs(min(all_err)), abs(max(all_err))) * 1.2

    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.6), dpi=300,
                             gridspec_kw=dict(wspace=0.30,
                                              left=0.08, right=0.92,
                                              top=0.88, bottom=0.17))

    for col, site in enumerate(SITES):
        data = site_data[site]
        z = data["z"]
        F_t = data["F_teacher"]
        F_s = data["F_student"]
        err = F_s - F_t

        ax_e = axes[col]

        # Force curves
        ax_e.plot(z, F_t, color=CLR_TEACHER, linewidth=2, linestyle="-",
                  label="MAD-SURF ", zorder=3)
        ax_e.plot(z, F_s, color=CLR_STUDENT, linewidth=2, linestyle=":",
                  label="GridFF ", zorder=3)

        ax_e.set_xlim(Z_MIN, Z_MAX)
        ax_e.set_ylim(0, F_hi)

        ax_e.set_xlabel("z (\u00c5)", fontsize=10)
        if col == 0:
            ax_e.set_ylabel("|Force| (eV/\u00c5)", fontsize=10)
        else:
            ax_e.tick_params(axis="y", labelleft=False)
        ax_e.tick_params(axis="both", labelsize=8.5, direction="in", top=True)

        # CRITICAL: hide right spine before twinx
        ax_e.spines["right"].set_visible(False)

        # Error on right y-axis
        ax_err = ax_e.twinx()
        ax_err.plot(z, err, color=CLR_ERROR, linewidth=0.8, alpha=0.50, zorder=2)
        ax_err.set_ylim(-err_lim, err_lim)

        ax_err.spines["left"].set_visible(False)
        ax_err.spines["top"].set_visible(False)

        if col == 2:
            ax_err.set_ylabel("Force error (eV/\u00c5)", fontsize=10,
                              color=CLR_ERROR, alpha=0.8)
            ax_err.tick_params(axis="y", labelsize=8, colors=CLR_ERROR,
                               direction="in")
        else:
            ax_err.tick_params(axis="y", labelright=False)

        ax_err.spines["right"].set_color(CLR_ERROR)
        ax_err.spines["right"].set_alpha(0.5)

        # Panel label
        panel_letter = chr(ord("a") + col)
        ax_e.text(0.04, 0.96, f"({panel_letter})", transform=ax_e.transAxes,
                  fontsize=10, fontweight="bold", va="top", ha="left")

        # RMSE annotation
        rmse = np.sqrt(np.mean(err ** 2))
        ax_e.text(0.96, 0.96,
                  f"RMSE {rmse:.3f} eV/\u00c5",
                  transform=ax_e.transAxes, fontsize=7, va="top", ha="right",
                  color="black",
                  bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#dddddd",
                            alpha=0.90, lw=0.5))

        ax_e.set_title(SITE_LABELS[site], fontsize=8, fontweight="bold", pad=4)

        # Legend only on first panel
        if col == 0:
            ax_e.legend(loc="lower right", fontsize=7.5, frameon=True,
                        fancybox=True, framealpha=0.9, edgecolor="#cccccc")

        # Export CSV
        export_csv(site, data)

    # Save
    png_path = OUT_DIR / "force_panel.png"
    pdf_path = OUT_DIR / "force_panel.pdf"
    fig.savefig(png_path, dpi=300, bbox_inches="tight", pad_inches=0.05)
    print(f"Saved: {png_path}")
    fig.savefig(pdf_path, bbox_inches="tight", pad_inches=0.05)
    print(f"Saved: {pdf_path}")
    plt.close(fig)


# ---------- main ----------
if __name__ == "__main__":
    make_figure()
