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
from scipy.spatial import ConvexHull
from matplotlib import font_manager
import os

# 2. Register all Times New Roman variants explicitly
font_files = [
    '/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman.ttf',
    '/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman_Bold.ttf',
    '/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman_Italic.ttf',
    '/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman_Bold_Italic.ttf'
]

for font_file in font_files:
    if os.path.exists(font_file):
        font_manager.fontManager.addfont(font_file)
    else:
        print(f"Warning: Font file not found - {font_file}")

mpl.rcParams['font.family'] = 'Times New Roman'

l_w=2.5
f_s=25
# Adjusting Matplotlib's default settings
mpl.rcParams.update({
    'axes.linewidth': l_w,          # Thickness of the axix lines
    'xtick.major.width': l_w,       # Thickness of major ticks on x-axis
    'ytick.major.width': l_w,       # Thickness of major ticks on y-axis
    'xtick.minor.width': 1.5,       # Thickness of minor ticks on x-axis
    'ytick.minor.width': 1.5,       # Thickness of minor ticks on y-axis
    'xtick.major.size': 6,          # Length of major ticks on x-axis
    'ytick.major.size': 6,          # Length of major ticks on y-axis
    'xtick.minor.size': 4,          # Length of minor ticks on x-axis
    'ytick.minor.size': 4,          # Length of minor ticks on y-axis
    'xtick.labelsize': f_s,         # Font Size of x-axis tick labels
    'ytick.labelsize': f_s,         # Font Size of y-axis tick labels
    'axes.labelsize': f_s,          # Font Size of axis labels
    'legend.fontsize': f_s,         # Font Size of legend
    'axes.titlesize': f_s,          # Font Size of axis titles
    'figure.titlesize': f_s,        # Font Size of figure titles
    'xtick.direction': 'in',        # X-axis ticks point inward
    'ytick.direction': 'in',        # Y-axis ticks point inward
    'xtick.top': True,              # Show ticks on top axis
    'xtick.bottom': True,           # Show ticks on bottom axis (default)
    'ytick.left': True,             # Show ticks on left axis (default)
    'ytick.right': True,            # Show ticks on right axis
    'axes.grid.which': 'both',      # Grid lines at both major and minor ticks
    'xtick.major.size': 4,          # # Optional: control tick length for inward ticks Slightly longer major ticks
    'xtick.minor.size': 2,          # Minor ticks
    # 'ytick.major.size': 8,
    # 'ytick.minor.size': 4
    # 'axes.titlepad': 20,
    # 'axes.labelpad': 20,


})


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


def _multi_layer_indices(positions: np.ndarray, indices: np.ndarray, num_layers: int = 3, tol: float = 0.35) -> np.ndarray:
    """Return subset of *indices* that belong to the top *num_layers* z-layers.

    Parameters
    ----------
    positions
        ``(N, 3)`` array for *one* frame.
    indices
        1-D array with indices of atoms to consider (typically substrate).
    num_layers
        Number of top layers to include.
    tol
        Thickness (Å) regarded as one layer.
    """
    zs = positions[indices, 2]

    # Find distinct z-layers by clustering z-coordinates
    z_sorted = np.sort(zs)[::-1]  # Sort in descending order
    layer_z_values = [z_sorted[0]]  # Start with the highest z

    for z in z_sorted[1:]:
        # If this z is significantly different from the last layer, it's a new layer
        if layer_z_values[-1] - z > tol:
            layer_z_values.append(z)
            if len(layer_z_values) >= num_layers:
                break

    # Include atoms from the identified layers
    selected_indices = []
    for layer_z in layer_z_values:
        layer_mask = np.abs(zs - layer_z) <= tol/2  # Half tolerance for layer membership
        selected_indices.extend(indices[layer_mask])

    return np.array(selected_indices, dtype=int)


def _calculate_common_limits(vis: MoleculeTrajectoryVisualizer, substrate_indices: np.ndarray, molecule_indices: np.ndarray) -> dict:
    """Calculate common axis limits for consistent subplot dimensions.

    Parameters
    ----------
    vis
        MoleculeTrajectoryVisualizer instance
    substrate_indices
        Indices of substrate atoms
    molecule_indices
        Indices of molecule atoms

    Returns
    -------
    dict
        Dictionary with 'x', 'y', 'z' keys containing (min, max) tuples
    """
    # Get all positions (substrate + molecule) across all frames
    all_indices = np.concatenate([substrate_indices, molecule_indices])
    all_positions = vis.atom_positions[:, all_indices, :]

    # Calculate limits for each dimension with some padding
    padding = 2.0  # Angstroms
    limits = {}
    ranges = {}

    for i, coord in enumerate(['x', 'y', 'z']):
        coord_min = np.min(all_positions[:, :, i]) - padding
        coord_max = np.max(all_positions[:, :, i]) + padding
        limits[coord] = (coord_min, coord_max)
        ranges[coord] = coord_max - coord_min

    # Make all ranges equal to the maximum range for square subplots
    max_range = max(ranges.values())

    for coord in ['x', 'y', 'z']:
        current_center = (limits[coord][0] + limits[coord][1]) / 2
        half_max_range = max_range / 2
        limits[coord] = (current_center - half_max_range, current_center + half_max_range)

    return limits


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
    num_mol_snapshots: int = 2,
    custom_mol_frames: list[int] | None = None,
    mol_colors: list[str] | None = None,
    custom_sub: bool = False,
    xlim: tuple[float, float] | None = None,
    ylim: tuple[float, float] | None = None,
):
    """Low-level helper that draws one 2-D projection (XY, XZ, or YZ)."""

    # --- substrate layer(s) (colour-coded) ---------------------------------
    element_colors = {"Na": "orange", "Cl": "cyan"}
    default_sub_col = "lightgray"

    # Create different alpha values for different layers to show depth
    unique_elements = np.unique([vis.atom_types[i] for i in top_layer_sub])

    for elem in unique_elements:
        elem_indices = top_layer_sub[[vis.atom_types[i] == elem for i in top_layer_sub]]
        if len(elem_indices) == 0:
            continue

        # Get z-coordinates to distinguish layers
        elem_positions = vis.atom_positions[0, elem_indices]
        z_coords = elem_positions[:, 2]

        # Group atoms by z-layers for different alpha values
        z_unique = np.unique(z_coords)
        layer_alphas = np.linspace(0.6, 0.2, len(z_unique))  # Higher layers more opaque

        for i, z_val in enumerate(z_unique):
            layer_mask = np.abs(z_coords - z_val) < 0.1  # Tolerance for same layer
            layer_indices = elem_indices[layer_mask]

            if len(layer_indices) > 0:
                xy = vis.atom_positions[0, layer_indices][:, [ix, iy]]

                if custom_sub and (xlim is not None or ylim is not None):
                    # Extend substrate across the entire axis range(s)
                    original_xy = xy.copy()

                    # Find the substrate pattern periods in both directions
                    x_coords = original_xy[:, 0]
                    y_coords = original_xy[:, 1]
                    x_min_orig = np.min(x_coords)
                    x_max_orig = np.max(x_coords)
                    y_min_orig = np.min(y_coords)
                    y_max_orig = np.max(y_coords)
                    x_period = x_max_orig - x_min_orig
                    y_period = y_max_orig - y_min_orig

                    # Create extended substrate positions
                    extended_xy = []

                    # Determine extension ranges
                    if xlim is not None:
                        x_target_min, x_target_max = xlim
                        if x_period > 0:
                            n_periods_x_left = int(np.ceil((x_min_orig - x_target_min) / x_period))
                            n_periods_x_right = int(np.ceil((x_target_max - x_max_orig) / x_period))
                        else:
                            n_periods_x_left = n_periods_x_right = 0
                    else:
                        x_target_min, x_target_max = x_min_orig, x_max_orig
                        n_periods_x_left = n_periods_x_right = 0

                    if ylim is not None:
                        y_target_min, y_target_max = ylim
                        if y_period > 0:
                            n_periods_y_left = int(np.ceil((y_min_orig - y_target_min) / y_period))
                            n_periods_y_right = int(np.ceil((y_target_max - y_max_orig) / y_period))
                        else:
                            n_periods_y_left = n_periods_y_right = 0
                    else:
                        y_target_min, y_target_max = y_min_orig, y_max_orig
                        n_periods_y_left = n_periods_y_right = 0

                    # Create periodic copies in both x and y directions
                    for nx in range(-n_periods_x_left, n_periods_x_right + 1):
                        for ny in range(-n_periods_y_left, n_periods_y_right + 1):
                            shifted_xy = original_xy.copy()
                            if x_period > 0:
                                shifted_xy[:, 0] += nx * x_period
                            if y_period > 0:
                                shifted_xy[:, 1] += ny * y_period
                            extended_xy.append(shifted_xy)

                    # Combine all positions
                    xy = np.vstack(extended_xy)

                    # Filter to only show atoms within the target ranges
                    mask = ((xy[:, 0] >= x_target_min) & (xy[:, 0] <= x_target_max) &
                            (xy[:, 1] >= y_target_min) & (xy[:, 1] <= y_target_max))
                    xy = xy[mask]

                ax.scatter(
                    xy[:, 0],
                    xy[:, 1],
                    s=60,
                    color=element_colors.get(elem, default_sub_col),
                    edgecolors="none",
                    linewidths=0.5,
                    label=f"{elem} (layer {i+1})" if i == 0 and elem in element_colors else None,
                    zorder=0,
                    alpha=layer_alphas[i],
                )
    

    # --- molecule convex hulls --------------------------------------------------
    # Only draw polygons if n_sample_structures > 0
    if n_sample_structures > 0:
        for frame_idx, fr in enumerate(frame_sel):
            pos = vis.atom_positions[fr, molecule_indices][:, [ix, iy]]
            hull = ConvexHull(pos)

            fixed_atom_idx_in_mol = molecule_indices[fixed_atom_idx]
            fixed_atom_pos = vis.atom_positions[fr, fixed_atom_idx_in_mol, [ix, iy]]
            opposite_atom_idx_in_mol = molecule_indices[opposite_atom_idx]
            opposite_atom_pos = vis.atom_positions[fr, opposite_atom_idx_in_mol, [ix, iy]]

            # Determine color for this frame to match molecule colors
            if mol_colors is None or len(mol_colors) == 0:
                anchor_color = 'k'
            elif len(mol_colors) == 1:
                anchor_color = mol_colors[0]
            else:
                anchor_color = mol_colors[frame_idx % len(mol_colors)]

            # Draw anchor points with matching colors and labels
            ax.scatter(fixed_atom_pos[0], fixed_atom_pos[1], color='r', marker="o", s=10, zorder=4)
            ax.scatter(opposite_atom_pos[0], opposite_atom_pos[1], color='b', marker="o", s=10, zorder=4)

            # Add numbered labels only on fixed atoms
            label_text = str(frame_idx + 1)  # 1-based numbering
            ax.annotate(label_text, (fixed_atom_pos[0], fixed_atom_pos[1]),
                       xytext=(5, 5), textcoords='offset points',
                       fontsize=12, fontweight='bold', color=anchor_color, zorder=5)

            # Draw convex hull polygon with same color scheme as molecules
            hull_poly = plt.Polygon(
                pos[hull.vertices],
                closed=True,
                linewidth=0.8,
                edgecolor=anchor_color,
                facecolor=anchor_color,
                alpha=0.2,
                linestyle='-',
                zorder=1
            )
            ax.add_patch(hull_poly)
    else:
        # When --samples 0, still show anchor points on custom molecule frames
        if custom_mol_frames is not None:
            snapshot_frames = np.array(custom_mol_frames, dtype=int)
            for frame_idx, fr in enumerate(snapshot_frames):
                fixed_atom_idx_in_mol = molecule_indices[fixed_atom_idx]
                fixed_atom_pos = vis.atom_positions[fr, fixed_atom_idx_in_mol, [ix, iy]]
                opposite_atom_idx_in_mol = molecule_indices[opposite_atom_idx]
                opposite_atom_pos = vis.atom_positions[fr, opposite_atom_idx_in_mol, [ix, iy]]

                # Determine color for this frame to match molecule colors
                if mol_colors is None or len(mol_colors) == 0:
                    anchor_color = 'k'
                elif len(mol_colors) == 1:
                    anchor_color = mol_colors[0]
                else:
                    anchor_color = mol_colors[frame_idx % len(mol_colors)]

                # Draw anchor points with matching colors and labels
                ax.scatter(fixed_atom_pos[0], fixed_atom_pos[1], color='r', marker="o", s=10, zorder=4)
                ax.scatter(opposite_atom_pos[0], opposite_atom_pos[1], color='b', marker="o", s=10, zorder=4)

                # Add numbered labels only on fixed atoms
                label_text = str(frame_idx + 1)  # 1-based numbering
                ax.annotate(label_text, (fixed_atom_pos[0], fixed_atom_pos[1]),
                           xytext=(5, 5), textcoords='offset points',
                           fontsize=12, fontweight='bold', color=anchor_color, zorder=5)

    ##--- molecule snapshots --------------------------------------------------
    if custom_mol_frames is not None:
        # Use custom molecule snapshot frames if provided
        snapshot_frames = np.array(custom_mol_frames, dtype=int)
        # Validate frame indices
        if np.any(snapshot_frames >= vis.num_frames) or np.any(snapshot_frames < 0):
            raise ValueError(f"Custom molecule frame indices must be between 0 and {vis.num_frames-1}")
    else:
        # Use default evenly spaced frames
        snapshot_frames = np.linspace(0, vis.num_frames-1, num_mol_snapshots, dtype=int)

    for frame_idx, fr in enumerate(snapshot_frames):
        pos_3d = vis.atom_positions[fr, molecule_indices]  # Keep 3D positions for bond calculation
        pos_2d = pos_3d[:, [ix, iy]]  # 2D projection for plotting

        # Determine color for this frame
        if mol_colors is None or len(mol_colors) == 0:
            # Default: all black
            mol_color = 'k'
        elif len(mol_colors) == 1:
            # Single color for all frames
            mol_color = mol_colors[0]
        else:
            # Cycle through provided colors
            mol_color = mol_colors[frame_idx % len(mol_colors)]

        # ax.scatter(pos_2d[:, 0], pos_2d[:, 1], s=10, color=mol_color, alpha=0.4, zorder=1)

        # bonds in 2-D projection - use 3D distances but plot in 2D
        # This prevents false bonds that appear close in 2D but are far in 3D
        for i in range(pos_3d.shape[0]):
            # Calculate 3D distances to avoid false bonds in projections
            d2_3d = np.sum((pos_3d - pos_3d[i]) ** 2, axis=1)
            neigh = np.where((d2_3d < bond_length_thresh**2) & (d2_3d > 0))[0]
            for j in neigh:
                if i < j:
                    ax.plot((pos_2d[i, 0], pos_2d[j, 0]),(pos_2d[i, 1], pos_2d[j, 1]),color=mol_color,alpha=0.8,linewidth=1.5, zorder=1,)

    # --- two special atom trajectories --------------------------------------
    for idx, col, lab in (
        (fixed_atom_idx, "red", f"Fixed atom {fixed_atom_idx}"),
        (opposite_atom_idx, "blue", f"Opposite atom {opposite_atom_idx}"),
    ):
        traj = vis.atom_positions[:, idx][:, [ix, iy]]
        ax.plot(traj[:, 0], traj[:, 1], color=col, marker="o", ms=0.5,lw=0.2, label=lab, zorder=3)
        ax.scatter(traj[0, 0], traj[0, 1], color=col, marker="o", s=10, zorder=4)
        ax.scatter(traj[-1, 0], traj[-1, 1], color=col, marker="s", s=10, zorder=4)

    ax.set_xlabel(f"{coord_labels[0]} (Å)", fontsize=f_s)
    ax.set_ylabel(f"{coord_labels[1]} (Å)", fontsize=f_s)
    ax.set_aspect("equal", adjustable="box")


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
    num_mol_snapshots: int = 2,
    ax=None,  # New parameter for external axis
    custom_frames: list[int] | None = None,  # Custom frame indices for polygons
    custom_mol_frames: list[int] | None = None,  # Custom frame indices for molecule snapshots
    mol_colors: list[str] | None = None,  # Colors for molecule snapshots
    xlim: tuple[float, float] | None = None,  # X-axis limits
    ylim: tuple[float, float] | None = None,  # Y-axis limits
    custom_sub: bool = False,  # Enable extended substrate display
):
    """Generate 2-D projection plots (xy, xz, yz) of the trajectory."""
    traj_path = Path(trajectory_file)
    if not traj_path.is_file():
        raise FileNotFoundError(traj_path)

    # Load trajectory
    vis = MoleculeTrajectoryVisualizer()
    vis.read_xyz_trajectory(str(traj_path))

    # Identify substrate and molecule atoms
    if substrate_types is None:
        substrate_indices = _detect_substrate_indices(vis.atom_types)
    else:
        substrate_indices = _detect_substrate_indices(vis.atom_types, substrate_types)

    if len(substrate_indices) == 0:
        raise RuntimeError(
            "Could not detect any substrate atoms – please specify --substrate-types."
        )

    molecule_indices = np.setdiff1d(np.arange(vis.num_atoms), substrate_indices)

    # Calculate common limits for consistent subplot dimensions
    common_limits = _calculate_common_limits(vis, substrate_indices, molecule_indices)

    # Frames to draw full molecule structures (polygons)
    if custom_frames is not None:
        # Use custom frame indices if provided
        frame_sel = np.array(custom_frames, dtype=int)
        # Validate frame indices
        if np.any(frame_sel >= vis.num_frames) or np.any(frame_sel < 0):
            raise ValueError(f"Custom frame indices must be between 0 and {vis.num_frames-1}")
        n_sample_structures = len(frame_sel)
    elif n_sample_structures > 0:
        n_sample_structures = max(2, n_sample_structures)
        frame_sel = np.linspace(0, vis.num_frames - 1, n_sample_structures, dtype=int)
    else:
        frame_sel = np.array([], dtype=int)
        n_sample_structures = 0

    # Handle plotting - either create new figure or use provided axis
    if ax is None:
        # Original behavior - create new figure
        n_plots = len(projections)
        fig, axes = plt.subplots(
            1,
            n_plots,
            figsize=(figsize_per_plot * n_plots, figsize_per_plot),
            squeeze=False,
        )
        axes = axes[0]
    else:
        # Use provided axis for single projection
        if len(projections) > 1:
            raise ValueError("Cannot plot multiple projections when ax is provided")
        axes = [ax]
        fig = ax.figure

    # Plot each projection
    proj_map = {
        "xy": (0, 1, ("X", "Y")),
        "xz": (0, 2, ("X", "Z")),
        "yz": (1, 2, ("Y", "Z")),
    }

    for ax, proj in zip(axes, projections):
        if proj not in proj_map:
            raise ValueError(f"Unknown projection '{proj}'. Choose from xy,xz,yz.")
        ix, iy, labels = proj_map[proj]

        # Use appropriate substrate layers based on projection
        if proj == "xy":
            # For XY projection, use only top layer
            substrate_layer = _top_layer_indices(
                vis.atom_positions[0], substrate_indices, tol=top_layer_tol
            )
        else:
            # For XZ and YZ projections, use multiple layers to show 3D structure
            substrate_layer = _multi_layer_indices(
                vis.atom_positions[0], substrate_indices, num_layers=3, tol=top_layer_tol
            )

        _plot_single_projection(
            ax,
            vis,
            ix,
            iy,
            labels,
            substrate_layer,
            molecule_indices,
            fixed_atom_idx,
            opposite_atom_idx,
            frame_sel,
            n_sample_structures,
            bond_length_thresh,
            num_mol_snapshots,
            custom_mol_frames,
            mol_colors,
            custom_sub,
            xlim,
            ylim,
        )

        # Apply axis limits - use custom limits if provided, otherwise use common limits
        coord_names = ['x', 'y', 'z']
        x_coord = coord_names[ix]
        y_coord = coord_names[iy]

        if xlim is not None:
            ax.set_xlim(xlim)
        else:
            ax.set_xlim(common_limits[x_coord])

        if ylim is not None:
            ax.set_ylim(ylim)
        else:
            ax.set_ylim(common_limits[y_coord])

    # Only add legend if creating new figure
    if ax is None:
        handles, labels = axes[0].get_legend_handles_labels()
        uniq = {}
        for h, l in zip(handles, labels):
            if l not in uniq:
                uniq[l] = h
        fig.tight_layout()

    if out_png is not None:
        fig.savefig(out_png, dpi=600, bbox_inches='tight')
        print(f"Saved figure to {out_png}")

    if show and ax is None:  # Only show if we created the figure
        plt.show()
    if ax is None:  # Only close if we created the figure
        plt.close(fig)
# end plot_top_layer_projections


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
        default=8,
        help="Number of molecule snapshots to display along trajectory (default: 8)",
    )
    p.add_argument(
        "--projections",
        default="xy",
        help="Comma-separated list of 2D projections to show: xy,xz,yz (default: xy)",
    )
    p.add_argument(
        "--num-mol-snapshots",
        type=int,
        default=2,
        help="Number of molecule snapshots to show. Default: 2",
    )
    p.add_argument(
        "--custom-frames",
        default=None,
        help="Comma-separated list of specific frame indices for polygons/snapshots (e.g., '0,10,20,30'). Overrides --samples if provided.",
    )
    p.add_argument(
        "--custom-mol-frames",
        default=None,
        help="Comma-separated list of specific frame indices for molecule snapshots (e.g., '5,15'). Overrides --num-mol-snapshots if provided.",
    )
    p.add_argument(
        "--mol-colors",
        default="auto",
        help="Color scheme for molecule snapshots: 'auto' (different colors per frame), 'same' (all black), or comma-separated color list (e.g., 'red,blue,green')",
    )
    p.add_argument(
        "--xlim",
        default=None,
        help="X-axis limits as 'min,max' (e.g., '-10,10'). If not provided, automatic limits are used.",
    )
    p.add_argument(
        "--ylim",
        default=None,
        help="Y-axis limits as 'min,max' (e.g., '-10,10'). If not provided, automatic limits are used.",
    )
    p.add_argument(
        "--custom-sub",
        action="store_true",
        help="Extend substrate atoms across the entire x-axis range specified by --xlim. Creates periodic substrate pattern.",
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

    # Parse custom frame indices if provided
    custom_frames = None
    if args.custom_frames:
        custom_frames = [int(f.strip()) for f in args.custom_frames.split(',') if f.strip()]

    custom_mol_frames = None
    if args.custom_mol_frames:
        custom_mol_frames = [int(f.strip()) for f in args.custom_mol_frames.split(',') if f.strip()]

    # Parse molecule colors
    mol_colors = None
    if args.mol_colors == "auto":
        # Generate different colors automatically
        import matplotlib.pyplot as plt
        cmap = plt.cm.tab10
        mol_colors = [cmap(i) for i in range(10)]  # Use first 10 colors from tab10
    elif args.mol_colors == "same":
        mol_colors = ['k']  # All black
    elif args.mol_colors != "auto":
        # Custom color list
        mol_colors = [c.strip() for c in args.mol_colors.split(',') if c.strip()]

    # Parse axis limits
    xlim = None
    if args.xlim:
        xlim_parts = [float(x.strip()) for x in args.xlim.split(',')]
        if len(xlim_parts) != 2:
            raise ValueError("--xlim must be in format 'min,max'")
        xlim = tuple(xlim_parts)

    ylim = None
    if args.ylim:
        ylim_parts = [float(y.strip()) for y in args.ylim.split(',')]
        if len(ylim_parts) != 2:
            raise ValueError("--ylim must be in format 'min,max'")
        ylim = tuple(ylim_parts)

    plot_top_layer_projections(
        trajectory_file=args.traj,
        fixed_atom_idx=args.fixed,
        opposite_atom_idx=args.opposite,
        projections=projections,
        n_sample_structures=args.samples,
        substrate_types=[s.strip() for s in args.substrate_types.split(',') if s.strip()],
        show=args.show,
        out_png=args.out,
        num_mol_snapshots=args.num_mol_snapshots,
        custom_frames=custom_frames,
        custom_mol_frames=custom_mol_frames,
        mol_colors=mol_colors,
        xlim=xlim,
        ylim=ylim,
        custom_sub=args.custom_sub,
    )


if __name__ == "__main__":
    main()


'''


###For all the projections plot (xy, xz, yz) can also be done for any individual and the pairs
 python /home/indranil/git/FireCore/doc/py/visualize_top_layer_xy.py  --traj dir_2.0_1.0_0.0/cons_26/PTCDA_20x20_26_total_trajectory.xyz --fixed 26 --opposite 29 --projections xy,xz,yz


## For one projection
python /home/indranil/git/FireCore/doc/py/visualize_top_layer_xy.py  --traj relax_defect_aligned_line/dir_1.0_1.0_0.0/cons_26/PTCDA_20x20_26_total_trajectory.xyz --fixed 26 --opposite 29 --samples 6 --num-mol-snapshots 6 --out PTCDA_defect_aligned_20x20_26_trajectory_xy.png


## For custom frames
python /home/indranil/git/FireCore/doc/py/visualize_top_layer_xy.py  --traj PTCDA_data_trial_1d_relax_z/old_mol_old_sub_PTCDA_total_trajectory.xyz --fixed 26 --opposite 29 --custom-frames "0,5,10,15" --custom-mol-frames "2,8" --out test_custom_frames.png --projections xy,xz,yz


## For custom colors
python /home/indranil/git/FireCore/doc/py/visualize_top_layer_xy.py  --traj PTCDA_data_trial_1d_relax_z/old_mol_old_sub_PTCDA_total_trajectory.xyz --fixed 26 --opposite 29 --custom-frames "0,5,10,15" --custom-mol-frames "2,8" --out test_custom_colors.png --projections xy,xz,yz --mol-colors "auto"


## For same color
python /home/indranil/git/FireCore/doc/py/visualize_top_layer_xy.py  --traj PTCDA_data_trial_1d_relax_z/old_mol_old_sub_PTCDA_total_trajectory.xyz --fixed 26 --opposite 29 --custom-frames "0,5,10,15" --custom-mol-frames "2,8" --out test_same_color.png --projections xy,xz,yz --mol-colors "same"



python /home/indranil/git/FireCore/doc/py/visualise_trajectory_atoms.py  --traj PTCDA_data_trial_1d_relax_z/old_mol_old_sub_PTCDA_total_trajectory.xyz --fixed 26 --opposite 29 --custom-frames "28,84,101,127,200" --custom-mol-frames "28,84,101,127,200" --out test_custom_colors.png --projections xz --mol-colors "auto" --xlim="-15,10" --ylim="-0.5,30"



python /home/indranil/git/FireCore/doc/py/visualise_trajectory_atoms.py  --traj PTCDA_data_trial_1d_relax_z/old_mol_old_sub_PTCDA_total_trajectory.xyz --fixed 26 --opposite 29 --custom-frames "28,84,101,127,200" --custom-mol-frames "28,84,101,127,200" --out test_custom_colors.png --projections xz --mol-colors "auto" --xlim="-12,10" --ylim="-0.6,30" --custom-sub



'''
