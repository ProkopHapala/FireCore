#!/usr/bin/env python3
"""visualize_top_layer_xy.py

Utility script to produce an illustrative XY-plane plot that combines

1.  The *topmost* layer of a (rigid) crystalline substrate
2.  The 2-D trajectories of two selected atoms of an adsorbed molecule –
    typically a *fixed* reference atom and the atom diagonally opposite to
    it – recorded during a scan/relax trajectory.
3.  A handful of whole-molecule snapshots along the trajectory, drawn in
    faint colours to give a feeling for the molecular conformation changes
    without obscuring the trajectories.

The script is intentionally self-contained and can be executed directly
from the command line, e.g.::

    python visualize_top_layer_xy.py \
        --traj tests/tMMFF/relax_perfect_line/dir_1.0_1.0_0.0/cons_26/\
                PTCDA_12x12_26_total_trajectory.xyz \
        --fixed 26 --opposite 29

Interactivity is not required – the result is saved as PNG (and optionally
shown).
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable, Sequence

import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl

# Re-use the generic visualiser that already exists in FireCore
from plot_utils import MoleculeTrajectoryVisualizer


# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

def _detect_substrate_indices(
    atom_types: Sequence[str],
    allowed_types: Iterable[str] = ("Na", "Cl", "K", "Br", "F", "I"),
) -> np.ndarray:
    """Return indices of atoms whose element *symbol* is in *allowed_types*."""

    allowed_set = set(allowed_types)
    return np.array([i for i, t in enumerate(atom_types) if t in allowed_set])


def _top_layer_indices(positions: np.ndarray, indices: np.ndarray, tol: float = 0.35) -> np.ndarray:
    """Return subset of *indices* that belong to the *uppermost* z-layer.

    Parameters
    ----------
    positions
        ``(N, 3)`` array for *one* frame.
    indices
        1-D array with indices of atoms to consider (typically substrate).
    tol
        Thickness (Å) regarded as *top* layer.
    """
    zs = positions[indices, 2]
    z_max = zs.max()
    return indices[zs >= z_max - tol]


# -----------------------------------------------------------------------------
# Core plotting routine
# -----------------------------------------------------------------------------

def _plot_single_projection(
    ax: "plt.Axes",
    vis: MoleculeTrajectoryVisualizer,
    ix: int,
    iy: int,
    coord_labels: tuple[str, str],
    top_layer_sub: np.ndarray,
    molecule_indices: np.ndarray,
    fixed_atom_idx: int,
    opposite_atom_idx: int,
    frame_sel: np.ndarray,
    n_sample_structures: int,
    bond_length_thresh: float,
):
    """Low-level helper that draws one 2-D projection (XY, XZ, or YZ)."""

    # --- substrate top layer (colour-coded) ---------------------------------
    element_colors = {"Na": "orange", "Cl": "cyan"}
    default_sub_col = "lightgray"
    for elem in np.unique([vis.atom_types[i] for i in top_layer_sub]):
        elem_indices = top_layer_sub[[vis.atom_types[i] == elem for i in top_layer_sub]]
        xy = vis.atom_positions[0, elem_indices][:, [ix, iy]]
        ax.scatter(
            xy[:, 0],
            xy[:, 1],
            s=60,
            color=element_colors.get(elem, default_sub_col),
            edgecolors="k",
            linewidths=0.5,
            label=f"{elem} (top layer)" if elem in element_colors else "Substrate atom",
            zorder=0,
            alpha=0.3,
        )

    # --- molecule snapshots --------------------------------------------------
    cmap = mpl.colormaps.get_cmap("rainbow")
    snap_colors = [cmap(i) for i in np.linspace(0, 1, len(frame_sel))]
    for c, fr in zip(snap_colors, frame_sel):
        pos = vis.atom_positions[fr, molecule_indices][:, [ix, iy]]
        ax.scatter(pos[:, 0], pos[:, 1], s=10, color=c, alpha=0.4, zorder=1)
        # bonds in 2-D projection
        for i in range(pos.shape[0]):
            d2 = np.sum((pos - pos[i]) ** 2, axis=1)
            neigh = np.where((d2 < bond_length_thresh**2) & (d2 > 0))[0]
            for j in neigh:
                if i < j:
                    ax.plot(
                        (pos[i, 0], pos[j, 0]),
                        (pos[i, 1], pos[j, 1]),
                        color=c,
                        alpha=0.8,
                        linewidth=0.8,
                        zorder=1,
                    )

    # --- two special atom trajectories --------------------------------------
    for idx, col, lab in (
        (fixed_atom_idx, "red", f"Fixed atom {fixed_atom_idx}"),
        (opposite_atom_idx, "blue", f"Opposite atom {opposite_atom_idx}"),
    ):
        traj = vis.atom_positions[:, idx][:, [ix, iy]]
        ax.plot(traj[:, 0], traj[:, 1], color=col, lw=2, label=lab, zorder=3)
        ax.scatter(traj[0, 0], traj[0, 1], color=col, marker="o", s=60, zorder=4)
        ax.scatter(traj[-1, 0], traj[-1, 1], color=col, marker="s", s=60, zorder=4)

    ax.set_xlabel(f"{coord_labels[0]} (Å)")
    ax.set_ylabel(f"{coord_labels[1]} (Å)")
    ax.set_aspect("equal", adjustable="box")
    # ax.grid(True)


def plot_top_layer_projections(
    *,
    trajectory_file: Path | str,
    fixed_atom_idx: int,
    opposite_atom_idx: int,
    projections: Iterable[str] = ("xy",),
    n_sample_structures: int = 6,
    substrate_types: Iterable[str] | None = None,
    top_layer_tol: float = 0.35,
    bond_length_thresh: float = 2.0,
    figsize_per_plot: int = 6,
    out_png: Path | None = None,
    show: bool = True,
):
    """Generate 2-D projection plots (xy, xz, yz) of the trajectory."""

    traj_path = Path(trajectory_file)
    if not traj_path.is_file():
        raise FileNotFoundError(traj_path)

    # ---------------------------------------------------------------------
    # Load trajectory
    # ---------------------------------------------------------------------
    vis = MoleculeTrajectoryVisualizer()
    vis.read_xyz_trajectory(str(traj_path))

    # Identify substrate atoms via their element symbols and grab the top layer.
    if substrate_types is None:
        substrate_indices = _detect_substrate_indices(vis.atom_types)
    else:
        substrate_indices = _detect_substrate_indices(vis.atom_types, substrate_types)

    if len(substrate_indices) == 0:
        raise RuntimeError(
            "Could not detect any substrate atoms – please specify --substrate-types."
        )

    top_layer_sub = _top_layer_indices(
        vis.atom_positions[0], substrate_indices, tol=top_layer_tol
    )

    # Molecule atoms = everything that is *not* substrate
    molecule_indices = np.setdiff1d(np.arange(vis.num_atoms), substrate_indices)

    # Frames to draw full molecule structures
    n_sample_structures = max(2, n_sample_structures)
    frame_sel = np.linspace(0, vis.num_frames - 1, n_sample_structures, dtype=int)

    # ---------------------------------------------------------------------
    # Plotting for requested projections
    # ---------------------------------------------------------------------
    n_plots = len(projections)
    fig, axes = plt.subplots(
        1,
        n_plots,
        figsize=(figsize_per_plot * n_plots, figsize_per_plot),
        squeeze=False,
    )
    axes = axes[0]  # 1D list

    # mapping of projection string to coordinate indices and labels
    proj_map = {
        "xy": (0, 1, ("X", "Y")),
        "xz": (0, 2, ("X", "Z")),
        "yz": (1, 2, ("Y", "Z")),
    }

    for ax, proj in zip(axes, projections):
        if proj not in proj_map:
            raise ValueError(f"Unknown projection '{proj}'. Choose from xy,xz,yz.")
        ix, iy, labels = proj_map[proj]
        _plot_single_projection(
            ax,
            vis,
            ix,
            iy,
            labels,
            top_layer_sub,
            molecule_indices,
            fixed_atom_idx,
            opposite_atom_idx,
            frame_sel,
            n_sample_structures,
            bond_length_thresh,
        )
        ax.set_title(proj.upper())

    # common legend (take from first axis)
    handles, labels = axes[0].get_legend_handles_labels()
    uniq = {}
    for h, l in zip(handles, labels):
        if l not in uniq:
            uniq[l] = h
    fig.legend(uniq.values(), uniq.keys(), loc="upper left")

    fig.tight_layout()

    # save if requested
    if out_png is not None:
        fig.savefig(out_png, dpi=300, bbox_inches='tight')
        print(f"Saved figure to {out_png}")

    if show:
        plt.show()
    plt.close(fig)
# end plot_top_layer_projections
    element_colors = {
        "Cl": "cyan",
        "Na": "orange",
    }
    # fall-back colour
    default_sub_col = "lightgray"

    for elem in np.unique([vis.atom_types[i] for i in top_layer_sub]):
        elem_indices = top_layer_sub[[vis.atom_types[i] == elem for i in top_layer_sub]]
        xy = vis.atom_positions[0, elem_indices, :2]
        ax.scatter(
            xy[:, 0],
            xy[:, 1],
            s=60,
            color=element_colors.get(elem, default_sub_col),
            edgecolors="k",
            linewidths=0.5,
            label=f"{elem} (top layer)" if elem in element_colors else "Substrate atom",
            zorder=0,
            alpha=0.3,
        )

    # --- full-molecule snapshots
    cmap = mpl.colormaps.get_cmap("rainbow")
    snap_colors = [cmap(i) for i in np.linspace(0, 1, len(frame_sel))]

    for c, fr in zip(snap_colors, frame_sel):
        pos = vis.atom_positions[fr, molecule_indices]
        # scatter
        ax.scatter(pos[:, 0], pos[:, 1], s=10, color=c, alpha=0.4, zorder=1)
        # simple bonding within threshold (2-D projection)
        for i in range(pos.shape[0]):
            d2 = np.sum((pos - pos[i]) ** 2, axis=1)
            neigh = np.where((d2 < bond_length_thresh**2) & (d2 > 0))[0]
            for j in neigh:
                if i < j:
                    ax.plot(
                        (pos[i, 0], pos[j, 0]),
                        (pos[i, 1], pos[j, 1]),
                        color=c,
                        alpha=0.4,
                        linewidth=0.8,
                        zorder=1,
                    )

    # --- trajectories of the two special atoms
    for idx, col, lab in (
        (fixed_atom_idx, "red", f"Fixed atom {fixed_atom_idx}"),
        (opposite_atom_idx, "blue", f"Opposite atom {opposite_atom_idx}"),
    ):
        traj_xy = vis.atom_positions[:, idx, :2]
        ax.plot(traj_xy[:, 0], traj_xy[:, 1], color=col, lw=2, label=lab, zorder=3)
        ax.scatter(traj_xy[0, 0], traj_xy[0, 1], color=col, marker="o", s=60, zorder=4)
        ax.scatter(traj_xy[-1, 0], traj_xy[-1, 1], color=col, marker="s", s=60, zorder=4)

    ax.set_xlabel("X (Å)")
    ax.set_ylabel("Y (Å)")
    ax.set_aspect("equal", adjustable="box")

    title = f"XY trajectory (fixed atom {fixed_atom_idx})"
    ax.set_title(title)
    fig.tight_layout()

    if out_png is None:
        out_png = traj_path.with_name(traj_path.stem + "_plot_xy.png")
    fig.savefig(out_png, dpi=300)
    print(f"Saved figure to {out_png}")

    if show:
        plt.show()
    plt.close(fig)


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def _parse_arguments() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Plot XY trajectory with substrate top layer.")
    p.add_argument("--traj", required=True, help="Path to XYZ trajectory file")
    p.add_argument("--fixed", type=int, required=True, help="Index of fixed atom (0-based)")
    p.add_argument(
        "--opposite", type=int, required=True, help="Index of opposite/diagonal atom (0-based)"
    )
    p.add_argument(
        "--substrate-types",
        default="Na,Cl",
        help="Comma-separated list of element symbols considered substrate (default: Na,Cl)",
    )
    p.add_argument(
        "--samples",
        type=int,
        default=6,
        help="Number of molecule snapshots to display along trajectory (default: 6)",
    )
    p.add_argument(
        "--projections",
        default="xy",
        help="Comma-separated list of 2D projections to show: xy,xz,yz (default: xy)",
    )
    # Mutually exclusive show / no-show flags
    p.add_argument(
        "--show",
        dest="show",
        action="store_true",
        default=True,
        help="Show matplotlib GUI window (default)",
    )
    p.add_argument(
        "--no-show",
        dest="show",
        action="store_false",
        help="Do not display the matplotlib window",
    )
    p.add_argument(
        "--out",
        default=None,
        help="Optional output PNG filename (defaults to <traj>_xy.png)",
    )
    return p.parse_args()


def main() -> None:
    args = _parse_arguments()
    projections = [p.strip() for p in args.projections.split(',') if p.strip()]
    plot_top_layer_projections(
        trajectory_file=args.traj,
        fixed_atom_idx=args.fixed,
        opposite_atom_idx=args.opposite,
        projections=projections,
        n_sample_structures=args.samples,
        substrate_types=[s.strip() for s in args.substrate_types.split(',') if s.strip()],
        show=args.show,
        out_png=args.out,
    )


if __name__ == "__main__":
    main()


'''
 python /home/indranil/git/FireCore/doc/py/visualize_top_layer_xy.py  --traj dir_2.0_1.0_0.0/cons_26/PTCDA_20x20_26_total_trajectory.xyz --fixed 26 --opposite 29  

###For all the projections plot (xy, xz, yz) can also be done for any individual and the pairs
 python /home/indranil/git/FireCore/doc/py/visualize_top_layer_xy.py  --traj dir_2.0_1.0_0.0/cons_26/PTCDA_20x20_26_total_trajectory.xyz --fixed 26 --opposite 29 --projections xy,xz,yz

python /home/indranil/git/FireCore/doc/py/visualize_top_layer_xy.py  --traj dir_2.0_1.0_0.0/cons_26/PTCDA_20x20_26_total_trajectory.xyz --fixed 26 --opposite 29 --projections xy,xz,yz



'''
