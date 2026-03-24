#!/usr/bin/env python3
"""Render CHONH2 at top/bridge/hollow sites on Ag(111) using ASE + matplotlib.

Produces three transparent PNG images that slot into the top row of the
publication panel figure.

Usage:
    python render_sites_ase.py [--out-dir .]

Then:
    python plot_publication_panel.py --structure-dir .
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
import numpy as np

from ase import Atoms
from ase.build import fcc111
from ase.visualize.plot import plot_atoms

# ---------- Configuration ----------
# Molecule: CHONH2 (formamide)
MOL_SYMBOLS = ["C", "O", "N", "H", "H", "H"]
MOL_POSITIONS = np.array([
    [ 0.000,  0.000,  0.000],  # C (anchor)
    [ 1.200,  0.000,  0.000],  # O
    [-0.689,  1.134,  0.000],  # N
    [-0.548, -0.948,  0.000],  # H on C
    [-1.693,  1.101,  0.000],  # H on N
    [-0.220,  2.038,  0.000],  # H on N
])

# CPK colors
ELEMENT_COLORS = {
    "Ag": (0.75, 0.75, 0.78),
    "C":  (0.20, 0.20, 0.20),
    "O":  (0.85, 0.15, 0.15),
    "N":  (0.15, 0.25, 0.85),
    "H":  (0.95, 0.95, 0.95),
}

# Atom display radii (for visual, not physical)
DISPLAY_RADII = {"Ag": 0.8, "C": 0.45, "O": 0.42, "N": 0.44, "H": 0.28}


def _site_uv(site: str) -> tuple[float, float]:
    """Fractional uv coordinates for each high-symmetry site on fcc(111)."""
    if site == "top":
        return (0.0, 0.0)
    elif site == "bridge":
        return (0.5, 0.0)
    elif site == "hollow":
        return (1.0 / 3.0, 1.0 / 3.0)
    raise ValueError(site)


def build_scene(site: str, z_height: float = 2.7) -> Atoms:
    """Build an ASE Atoms object: Ag(111) slab + CHONH2 at the given site."""
    slab = fcc111("Ag", size=(4, 4, 3), a=4.071, vacuum=6.0, periodic=True)
    cell = slab.cell

    # Place molecule at site position
    u, v = _site_uv(site)
    anchor_xy = u * cell[0, :2] + v * cell[1, :2]
    z_surface = slab.positions[:, 2].max()

    mol_pos = MOL_POSITIONS.copy()
    # Tilt molecule 15° from vertical for visual interest (rotate around x-axis)
    theta = np.radians(15)
    rot = np.array([[1, 0, 0], [0, np.cos(theta), -np.sin(theta)], [0, np.sin(theta), np.cos(theta)]])
    mol_pos = mol_pos @ rot.T
    mol_pos[:, 0] += anchor_xy[0]
    mol_pos[:, 1] += anchor_xy[1]
    mol_pos[:, 2] += z_surface + z_height

    mol = Atoms(symbols=MOL_SYMBOLS, positions=mol_pos)
    scene = slab + mol
    return scene


def render_site(site: str, out_path: Path, figsize=(3.0, 3.0)):
    """Render one site view from above (top-down) with slight perspective."""
    scene = build_scene(site, z_height=2.7)

    fig, ax = plt.subplots(1, 1, figsize=figsize, dpi=200)

    # ASE plot_atoms: top-down view (rotation='0x,0y,0z' = looking down z)
    plot_atoms(scene, ax, radii=0.5, rotation=("-75x,0y,0z"), show_unit_cell=0)

    # Add site label
    ax.set_title(f"{site.capitalize()} site", fontsize=11, fontweight="bold", pad=8)

    ax.axis("off")
    ax.set_aspect("equal")

    fig.savefig(out_path, dpi=200, bbox_inches="tight", pad_inches=0.02,
                transparent=True)
    plt.close(fig)
    print(f"  {out_path.name}")


def render_site_matplotlib(site: str, out_path: Path, figsize=(2.8, 2.4)):
    """Render one site view using custom matplotlib drawing (side view for z-scan context)."""
    scene = build_scene(site, z_height=2.7)
    symbols = scene.get_chemical_symbols()
    pos = scene.positions.copy()

    # Centre the view on the molecule anchor
    mol_mask = np.array([s != "Ag" for s in symbols])
    mol_centre_x = pos[mol_mask, 0].mean()

    fig, ax = plt.subplots(1, 1, figsize=figsize, dpi=250)

    # Side view: x horizontal, z vertical
    # Only show Ag atoms within ±5 A of molecule centre (avoid clutter)
    view_half = 5.5

    # Sort by y (depth) for painter's algorithm
    order = np.argsort(-pos[:, 1])

    for i in order:
        sym = symbols[i]
        x, y, z = pos[i]
        # Skip Ag atoms far from molecule
        if sym == "Ag" and abs(x - mol_centre_x) > view_half + 1.0:
            continue
        r = DISPLAY_RADII.get(sym, 0.4)
        color = ELEMENT_COLORS.get(sym, (0.5, 0.5, 0.5))
        # Depth cue
        y_range = pos[:, 1].max() - pos[:, 1].min()
        depth_frac = (y - pos[:, 1].min()) / max(y_range, 1e-6)
        brightness = 0.65 + 0.35 * depth_frac
        display_color = tuple(min(1.0, c * brightness) for c in color)
        display_r = r * (0.85 + 0.25 * depth_frac)

        # Molecule atoms: larger radii + thicker edge for emphasis
        if sym != "Ag":
            display_r *= 1.5
            lw = 0.6
        else:
            lw = 0.3

        circle = plt.Circle((x, z), display_r, fc=display_color,
                             ec=(0.25, 0.25, 0.25), lw=lw,
                             zorder=int(depth_frac * 100) + (200 if sym != "Ag" else 0))
        ax.add_patch(circle)

        # Element labels on molecule atoms
        if sym != "Ag":
            label_color = "white" if sym in ("C", "N") else "black"
            ax.text(x, z, sym, ha="center", va="center", fontsize=6,
                    fontweight="bold", color=label_color, zorder=500)

    # z-scan arrow
    z_surf = max(p[2] for p, s in zip(pos, symbols) if s == "Ag")
    z_mol_top = max(p[2] for p, s in zip(pos, symbols) if s != "Ag")
    x_arrow = mol_centre_x + 3.5

    ax.annotate("", xy=(x_arrow, z_surf + 0.2), xytext=(x_arrow, float(z_mol_top) + 0.8),
                arrowprops=dict(arrowstyle="<->", color="#555555", lw=1.0))
    ax.text(x_arrow + 0.4, float((z_surf + z_mol_top) / 2) + 0.4, "z",
            fontsize=10, color="#555555", fontstyle="italic", va="center")

    # Surface dashed line
    ax.axhline(z_surf, color="#bbbbbb", linewidth=0.5, linestyle="--", zorder=0)

    # Mark the site position with a small triangle on the surface
    site_x = mol_centre_x
    ax.plot(site_x, z_surf + 0.15, "^", color="#1a6e1a", markersize=5, zorder=150, alpha=0.7)

    ax.set_xlim(mol_centre_x - view_half, mol_centre_x + view_half)
    ax.set_ylim(pos[~mol_mask, 2].min() - 0.8, z_mol_top + 2.5)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_facecolor("white")

    fig.savefig(out_path, dpi=250, bbox_inches="tight", pad_inches=0.02,
                transparent=True)
    plt.close(fig)
    print(f"  {out_path.name}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=str, default=".")
    parser.add_argument("--mode", choices=["ase", "custom"], default="custom",
                        help="'ase' uses ASE plot_atoms; 'custom' uses hand-drawn side view")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(exist_ok=True)

    print("Rendering site structures:")
    for site in ["top", "bridge", "hollow"]:
        out_path = out_dir / f"{site}_site.png"
        if args.mode == "ase":
            render_site(site, out_path)
        else:
            render_site_matplotlib(site, out_path)

    print(f"\nDone! Now run:")
    print(f"  python plot_publication_panel.py --structure-dir {out_dir}")
