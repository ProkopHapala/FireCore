#!/usr/bin/env python3
"""Milestone A: extract Fireball SCF data for PyOpenCL workflows.

This script drives the Fortran Fireball library via :mod:`pyBall.FireCore`
so that downstream PyOpenCL steps can operate on NumPy artefacts instead of
C/C++ bindings. It follows the test harness conventions requested for the AFM
migration: a CLI built with :mod:`argparse`, structured outputs, and diagnostic
plots created with Matplotlib.

Outputs generated per run:
- ``firecore_density_matrix.npy`` – molecular-orbital coefficient matrix
  (shape ``norbitals × norbitals``) saved as double precision.
- ``firecore_density_grid.npy`` – real-space density on the Fireball grid
  (3D NumPy array, single precision).
- ``firecore_grid_meta.json`` – grid metadata (``ngrid``, ``dcell``, ``lvs``).
- ``firecore_density_slice.png`` – visual slice through the density field.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Sequence

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 - activates 3D projection
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

import sys
sys.path.append("../../")

from pyBall import FireCore as fc
from pyBall import atomicUtils as au
from pyBall.DFT import utils as dft_utils


@dataclass
class GridMetadata:
    ngrid: tuple[int, int, int]
    dcell: list[list[float]]
    lvs: list[list[float]]
    def to_json(self) -> str:
        return json.dumps(asdict(self), indent=2)

def ensure_output_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path

def run_fireball(atom_types: Sequence[int], atom_pos: np.ndarray, *, verbosity: int, nmax_scf: int, ecut: float):
    atom_types_np = np.asarray(atom_types, dtype=np.int32)
    atom_pos_np = np.asarray(atom_pos, dtype=np.float64)
    print("py.run_fireball(): verbosity", verbosity)
    fc.setVerbosity(verbosity)
    fc.initialize(atomType=atom_types_np, atomPos=atom_pos_np, verbosity=verbosity)
    fc.SCF(atom_pos_np, nmax_scf=nmax_scf)
    ngrid, dcell, lvs = fc.setupGrid(Ecut=ecut)
    density = fc.getGridDens(ngrid=tuple(int(x) for x in ngrid))
    i0orb = dft_utils.countOrbs(atom_types_np)
    wfcoef = fc.get_wfcoef(norb=i0orb[-1])
    return wfcoef, density, GridMetadata(ngrid=tuple(int(x) for x in ngrid), dcell=dcell.tolist(), lvs=lvs.tolist())

def plot_density(
    out_path: Path,
    density: np.ndarray,
    *,
    mode: str = "slice",
    slice_axis: int = 2,
    slice_index: int | None = None,
    iso_level: float = 0.5,
    voxel_step: int = 1,
    show: bool = False,
) -> None:
    if density.ndim != 3:
        raise ValueError("density grid must be 3D")

    if mode == "slice":
        axis = slice_axis % 3
        if slice_index is None:
            slice_index = density.shape[axis] // 2
        slice_index = max(0, min(int(slice_index), density.shape[axis] - 1))

        if axis == 0:
            plane = density[slice_index, :, :]
            xlabel, ylabel, title_axis = "y index", "z index", "x"
        elif axis == 1:
            plane = density[:, slice_index, :]
            xlabel, ylabel, title_axis = "x index", "z index", "y"
        else:
            plane = density[:, :, slice_index]
            xlabel, ylabel, title_axis = "x index", "y index", "z"

        fig, ax = plt.subplots(figsize=(5, 4))
        im = ax.imshow(plane.T, origin="lower", cmap="viridis")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Density")
        ax.set_title(f"Fireball density slice {title_axis}={slice_index}")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        fig.tight_layout()

    elif mode == "voxels":
        if voxel_step < 1:
            raise ValueError("voxel_step must be >= 1")
        sub_density = density[::voxel_step, ::voxel_step, ::voxel_step]
        level = float(iso_level)
        if 0.0 < level <= 1.0:
            threshold = level * float(sub_density.max(initial=0.0))
        else:
            threshold = level
        mask = sub_density >= threshold
        if not np.any(mask):
            raise ValueError("iso-surface threshold produces empty mask; adjust iso_level")

        from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 - activates 3D projection
        #from mpl_toolkits.mplot3d.art3d import Poly3DCollection

        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(projection="3d")
        colors = np.zeros(mask.shape + (4,), dtype=np.float32)
        colors[mask] = (0.1, 0.4, 0.9, 0.3)
        ax.voxels(mask, facecolors=colors, edgecolors=None)
        ax.set_xlabel("x index")
        ax.set_ylabel("y index")
        ax.set_zlabel("z index")
        ax.set_title("Fireball density voxels")
        fig.tight_layout()

    elif mode == "isosurface":
        if voxel_step < 1:
            raise ValueError("voxel_step must be >= 1")
        try:
            from skimage import measure
        except ImportError as exc:
            raise RuntimeError(
                "Isosurface mode requires scikit-image (install via `pip install scikit-image`)."
            ) from exc

        sub_density = density[::voxel_step, ::voxel_step, ::voxel_step]
        level = float(iso_level)
        if 0.0 < level <= 1.0:
            threshold = level * float(sub_density.max(initial=0.0))
        else:
            threshold = level
        if threshold <= 0.0 or not np.isfinite(threshold):
            raise ValueError("iso-level must be positive for isosurface mode")

        verts, faces, normals, values = measure.marching_cubes(sub_density, level=threshold)
        verts *= voxel_step

        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(projection="3d")
        mesh = Poly3DCollection(verts[faces], alpha=0.35)
        mesh.set_facecolor((0.1, 0.4, 0.9, 0.35))
        mesh.set_edgecolor((0.05, 0.2, 0.6, 0.3))
        ax.add_collection3d(mesh)

        ax.set_xlim(0, sub_density.shape[0])
        ax.set_ylim(0, sub_density.shape[1])
        ax.set_zlim(0, sub_density.shape[2])
        ax.set_xlabel("x index")
        ax.set_ylabel("y index")
        ax.set_zlabel("z index")
        ax.set_title("Fireball density isosurface")
        fig.tight_layout()

    else:
        raise ValueError(f"Unknown plot mode '{mode}'")

    fig.savefig(out_path, dpi=150)
    if show:
        plt.show()
    else:
        plt.close(fig)


def plot_density_matrix(out_path: Path, matrix: np.ndarray, show: bool = False) -> None:
    fig, ax = plt.subplots(figsize=(5, 4))
    data = np.abs(matrix)
    im = ax.imshow(data, origin="lower", cmap="viridis")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="|density matrix|")
    ax.set_title("Fireball density matrix magnitude")
    ax.set_xlabel("orbital index")
    ax.set_ylabel("orbital index")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    if show:
        plt.show()
    else:
        plt.close(fig)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fireball SCF data extractor")
    parser.add_argument("xyz",           type=Path, nargs="?", default=Path("../../cpp/common_resources/xyz/CH4.xyz"), help="Path to XYZ file (default: CH4 example)")
    parser.add_argument("--output-dir",  type=Path,  default=Path("firecore_outputs"), help="Directory where NumPy arrays and plots will be stored")
    parser.add_argument("--verbosity",   type=int,   default=3,                        help="Firecore verbosity flag (passed to FireCore.initialize)")
    parser.add_argument("--nmax-scf",    type=int,   default=150,                      help="Maximum SCF iterations (passed to FireCore.SCF)")
    parser.add_argument("--ecut",        type=float, default=100.0,                    help="Plane-wave cutoff passed to FireCore.setupGrid")
    parser.add_argument("--save-prefix", type=str,   default="firecore",               help="Prefix for generated output files")
    parser.add_argument("--plot-mode",   type=str,   default="slice", choices=("slice", "voxels", "isosurface"), help="Visualization mode: 2D slice, 3D voxels, or smooth isosurface")
    parser.add_argument("--plot-slice",  type=int,   default=None,                     help="Optional index along slice axis (default: middle of grid)")
    parser.add_argument("--slice-axis",  type=int,   default=2,                       help="Axis for slice mode (0=x, 1=y, 2=z)")
    parser.add_argument("--iso-level",   type=float, default=0.5,                     help="Iso threshold for voxel mode; <=1 treated as fraction of max density")
    parser.add_argument("--voxel-step",  type=int,   default=1,                       help="Stride for voxel subsampling to keep plots manageable")
    parser.add_argument("--show",        type=int,   default=1,                   help="Display plots interactively instead of closing after save (0/1)")
    args = parser.parse_args()

    out_dir = ensure_output_dir(args.output_dir)

    # Load molecule
    apos, atom_types, _, _ = au.loadAtomsNP(args.xyz)

    # Run Fireball SCF and extract artefacts
    wfcoef, density, grid_meta = run_fireball(
        atom_types,
        apos,
        verbosity=args.verbosity,
        nmax_scf=args.nmax_scf,
        ecut=args.ecut,
    )

    # Persist arrays
    np.save(out_dir / f"{args.save_prefix}_density_matrix.npy", wfcoef)
    np.save(out_dir / f"{args.save_prefix}_density_grid.npy", density.astype(np.float32))
    (out_dir / f"{args.save_prefix}_grid_meta.json").write_text(grid_meta.to_json(), encoding="utf8")

    # Plot diagnostic slice
    show_flag = bool(args.show)

    plot_density(
        out_dir / f"{args.save_prefix}_density_{args.plot_mode}.png",
        density,
        mode=args.plot_mode,
        slice_axis=args.slice_axis,
        slice_index=args.plot_slice,
        iso_level=args.iso_level,
        voxel_step=args.voxel_step,
        show=show_flag,
    )

    plot_density_matrix(
        out_dir / f"{args.save_prefix}_density_matrix.png",
        wfcoef,
        show=show_flag,
    )

    plt.show()

    print(f"Saved Fireball artefacts to {out_dir.resolve()}")
