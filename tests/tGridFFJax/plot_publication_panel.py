#!/usr/bin/env python3
"""Publication-quality 2×3 panel figure for CHONH2/Ag(111) GridFF benchmark.

Top row:    Molecular structure snapshots at top / bridge / hollow sites
            (loaded from pre-rendered PNG images — see pymol_render_sites.py)
Bottom row: Energy z-scan curves + error on dual y-axis

Style matches the FireCore/LAMMPS benchmark figure (Morse+Coulomb on NaCl).

Usage:
    python plot_publication_panel.py [--structure-dir DIR] [--out figure_panel.pdf]

If structure images are not yet rendered, the top row shows schematic placeholders.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.image as mpimg
import numpy as np

# ---------- configuration ----------
SITES = ["top", "bridge", "hollow"]
SITE_LABELS = {"top": "Top", "bridge": "Bridge", "hollow": "Hollow"}
TRACE_DIR = Path(__file__).parent.parent / "traces"
Z_MIN = 2.0      # only plot z >= 2.0 A
Z_MAX = 5.6

# Colours — inspired by the reference image
CLR_TEACHER  = 'orange'    # near-black solid line
CLR_STUDENT  = "#1f77b4"    # dark green dashed
CLR_ERROR    = "#cc3333"    # red, faded for error axis

# ---------- helpers ----------
def load_site(site: str):
    d = np.load(TRACE_DIR / f"{site}_zscan_trace.npz")
    z = d["z_values"]
    mask = z >= Z_MIN
    # return {
    #     "z": z[mask],
    #     "E_teacher": d["teacher_energy"][mask] * 1000,   # meV
    #     "E_student": d["student_energy"][mask] * 1000,
    #     "train": d["train_mask"][mask],
    # }
    # Calculate force norm (eV/A) from (N, Atoms, 3) arrays
    ft = d["teacher_force"][mask]
    fs = d["student_force"][mask]
    return {
        "z": z[mask],
        "E_teacher": np.linalg.norm(ft.reshape(ft.shape[0], -1), axis=1),
        "E_student": np.linalg.norm(fs.reshape(fs.shape[0], -1), axis=1),
        "train": d["train_mask"][mask],
    }

def make_figure(structure_dir: Path | None, out_path: Path):
    fig = plt.figure(figsize=(7.2, 5.4), dpi=300)

    # GridSpec: 2 rows, 3 cols; top row shorter (structure images)
    gs = gridspec.GridSpec(
        2, 3,
        height_ratios=[1.0, 1.4],
        hspace=0.10, wspace=0.30,
        left=0.08, right=0.92, top=0.95, bottom=0.10,
    )

    # ===== TOP ROW: structure images =====
    # for col, site in enumerate(SITES):
    #     ax_img = fig.add_subplot(gs[0, col])
    #     img_path = None
    #     if structure_dir:
    #         for ext in (".png", ".jpg", ".tiff"):
    #             candidate = structure_dir / f"{site}_site{ext}"
    #             if candidate.exists():
    #                 img_path = candidate
    #                 break
    #     if img_path:
    #         img = mpimg.imread(str(img_path))
    #         ax_img.imshow(img, aspect="equal")
    #     else:
    #         # Placeholder — draw a schematic
    #         _draw_placeholder(ax_img, site)
    #     ax_img.set_title(SITE_LABELS[site], fontsize=11, fontweight="bold", pad=4)
    #     ax_img.axis("off")

    # ===== BOTTOM ROW: energy + error panels =====
    # Pre-compute global y-limits for consistent axes
    all_E = []
    all_err = []
    for site in SITES:
        data = load_site(site)
        all_E.extend(data["E_teacher"])
        all_E.extend(data["E_student"])
        all_err.extend(data["E_student"] - data["E_teacher"])
    
    # Auto-scale y-limits
    E_lo, E_hi = min(all_E) * 0.95, max(all_E) * 1.05
    err_lim = max(abs(min(all_err)), abs(max(all_err))) * 1.1

    for col, site in enumerate(SITES):
        data = load_site(site)
        z = data["z"]
        E_t = data["E_teacher"]
        E_s = data["E_student"]
        err = E_s - E_t

        ax_e = fig.add_subplot(gs[1, col])

        # Energy curves
        ax_e.plot(z, E_t, color=CLR_TEACHER, linewidth=2, linestyle="-",
                  label="MAD-SURF ", zorder=3)
        ax_e.plot(z, E_s, color=CLR_STUDENT, linewidth=2, linestyle=":",
                  label="GridFF ", zorder=3)

        # ax_e.axhline(0, color="#aaaaaa", linewidth=0.4, zorder=1)
        ax_e.set_xlim(Z_MIN, Z_MAX)
        ax_e.set_ylim(E_lo, E_hi)

        ax_e.set_xlabel("z (Å)", fontsize=10)
        if col == 0:
            ax_e.set_ylabel("Force Norm (eV/Å)", fontsize=10)
        else:
            ax_e.tick_params(axis="y", labelleft=False)
        ax_e.tick_params(axis="both", labelsize=8.5, direction="in", top=True)
        ax_e.spines["right"].set_visible(False)   # hide so twinx red spine isn't blocked

        # Error on right y-axis (kcal/mol, faded red)
        ax_err = ax_e.twinx()
        # ax_err.fill_between(z, 0, err_kcal, color=CLR_ERROR, alpha=0.10, zorder=0)
        ax_err.plot(z, err, color=CLR_ERROR, linewidth=0.8, alpha=0.50, zorder=2)
        # ax_err.axhline(0, color=CLR_ERROR, linewidth=0.3, alpha=0.25, zorder=0)
        ax_err.set_ylim(-err_lim, err_lim)

        # Hide left spine of error axis (inherited from main axis)
        ax_err.spines["left"].set_visible(False)
        ax_err.spines["top"].set_visible(False)

        if col == 2:
            ax_err.set_ylabel("Error (eV/Å)", fontsize=10, color=CLR_ERROR, alpha=0.8)
            ax_err.tick_params(axis="y", labelsize=8, colors=CLR_ERROR, direction="in")
        else:
            ax_err.tick_params(axis="y", labelright=False)

        ax_err.spines["right"].set_color(CLR_ERROR)
        ax_err.spines["right"].set_alpha(0.5)

        # Panel label
        panel_letter = chr(ord("a") + col)
        ax_e.text(0.04, 0.96, f"({panel_letter})", transform=ax_e.transAxes,
                  fontsize=10, fontweight="bold", va="top", ha="left")

        # RMSE annotation — compact box in upper-right
        rmse = np.sqrt(np.mean((E_s - E_t) ** 2))
        ax_e.text(0.96, 0.96,
                  f"RMSE {rmse:.3f} eV/Å",
                  transform=ax_e.transAxes, fontsize=7, va="top", ha="right",
                  color='black',
                  bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#dddddd",
                            alpha=0.90, lw=0.5))

        ax_e.set_title(SITE_LABELS[site], fontsize=8, fontweight="bold", pad=4)
        # Legend only on first panel
        if col == 0:
            ax_e.legend(loc="lower right", fontsize=7.5, frameon=True,
                        fancybox=True, framealpha=0.9, edgecolor="#cccccc")
        

    # Save
    fig.savefig(out_path, dpi=300, bbox_inches="tight", pad_inches=0.05)
    print(f"Saved: {out_path}")
    # Also save PDF for publication
    pdf_path = out_path.with_suffix(".pdf")
    fig.savefig(pdf_path, bbox_inches="tight", pad_inches=0.05)
    print(f"Saved: {pdf_path}")
    plt.close(fig)


def _draw_placeholder(ax, site):
    """Draw a minimal schematic placeholder when structure images are unavailable."""
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-0.5, 2.0)
    ax.set_aspect("equal")

    # Draw Ag surface atoms as grey circles
    offsets = {"top": 0.0, "bridge": 0.5, "hollow": 0.333}
    y_base = 0.0
    for ix in range(-2, 3):
        x = ix * 1.0
        ax.add_patch(plt.Circle((x, y_base), 0.35, fc="#c0c0c0", ec="#888888", lw=0.5, zorder=1))

    # Molecule marker
    mol_x = offsets.get(site, 0.0)
    ax.plot(mol_x, 1.3, "v", color=CLR_STUDENT, markersize=10, zorder=5)
    ax.text(mol_x, 1.7, "CHONH₂", ha="center", fontsize=7, color=CLR_STUDENT)

    # Arrow showing z-scan
    ax.annotate("", xy=(mol_x + 0.7, 0.5), xytext=(mol_x + 0.7, 1.5),
                arrowprops=dict(arrowstyle="<->", color="#666666", lw=0.8))
    ax.text(mol_x + 0.85, 1.0, "z", fontsize=7, color="#666666")
    ax.set_facecolor("#f8f8f8")


# ---------- main ----------
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--structure-dir", type=str, default=None,
                        help="Directory with top_site.png, bridge_site.png, hollow_site.png")
    parser.add_argument("--out", type=str, default=None,
                        help="Output filename (default: figure_panel.png in presentation_data/)")
    args = parser.parse_args()

    out_path = Path(args.out) if args.out else Path(__file__).parent / "figure_panel.png"
    struct_dir = Path(args.structure_dir) if args.structure_dir else None

    make_figure(struct_dir, out_path)
