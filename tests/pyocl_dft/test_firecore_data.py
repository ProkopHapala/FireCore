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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Fireball SCF data extractor")
    parser.add_argument("xyz",           type=Path, nargs="?", default=Path("tests/pyFireball/molecules/CH4.xyz"), help="Path to XYZ file (default: CH4 example)")
    parser.add_argument("--output-dir",  type=Path,  default=Path("firecore_outputs"), help="Directory where NumPy arrays and plots will be stored")
    parser.add_argument("--verbosity",   type=int,   default=1,                        help="Firecore verbosity flag (passed to FireCore.initialize)")
    parser.add_argument("--nmax-scf",    type=int,   default=150,                      help="Maximum SCF iterations (passed to FireCore.SCF)")
    parser.add_argument("--ecut",        type=float, default=100.0,                    help="Plane-wave cutoff passed to FireCore.setupGrid")
    parser.add_argument("--save-prefix", type=str,   default="firecore",               help="Prefix for generated output files")
    parser.add_argument("--plot-slice",  type=int,   default=None,                     help="Optional index along Z for density slice (default: middle plane)")
    parser.add_argument("--show",        type=int,   default=0, action="store_true",   help="Display plots interactively instead of closing after save")
    return parser.parse_args()


def ensure_output_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def run_fireball(atom_types: Sequence[int], atom_pos: np.ndarray, *, verbosity: int, nmax_scf: int, ecut: float):
    atom_types_np = np.asarray(atom_types, dtype=np.int32)
    atom_pos_np = np.asarray(atom_pos, dtype=np.float64)
    fc.initialize(atomType=atom_types_np, atomPos=atom_pos_np, verbosity=verbosity)
    fc.SCF(atom_pos_np, nmax_scf=nmax_scf)
    ngrid, dcell, lvs = fc.setupGrid(Ecut=ecut)
    density = fc.getGridDens(ngrid=tuple(int(x) for x in ngrid))
    i0orb = dft_utils.countOrbs(atom_types_np)
    wfcoef = fc.get_wfcoef(norb=i0orb[-1])
    return wfcoef, density, GridMetadata(ngrid=tuple(int(x) for x in ngrid), dcell=dcell.tolist(), lvs=lvs.tolist())


def save_numpy(path: Path, array: np.ndarray) -> None:
    np.save(path, array)


def plot_density_slice(out_path: Path, density: np.ndarray, iz: int | None, show: bool = False) -> None:
    if density.ndim != 3:
        raise ValueError("density grid must be 3D")
    if iz is None:
        iz = density.shape[2] // 2
    iz = max(0, min(int(iz), density.shape[2] - 1))
    fig, ax = plt.subplots(figsize=(5, 4))
    im = ax.imshow(density[:, :, iz].T, origin="lower", cmap="viridis")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Density")
    ax.set_title(f"Fireball density slice z={iz}")
    ax.set_xlabel("x index")
    ax.set_ylabel("y index")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    if show:
        plt.show()
    else:
        plt.close(fig)


def main() -> None:
    args = parse_args()
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
    save_numpy(out_dir / f"{args.save_prefix}_density_matrix.npy", wfcoef)
    save_numpy(out_dir / f"{args.save_prefix}_density_grid.npy", density.astype(np.float32))
    (out_dir / f"{args.save_prefix}_grid_meta.json").write_text(grid_meta.to_json(), encoding="utf8")

    # Plot diagnostic slice
    plot_density_slice(out_dir / f"{args.save_prefix}_density_slice.png", density, args.plot_slice, show=args.show)

    print(f"Saved Fireball artefacts to {out_dir.resolve()}")


if __name__ == "__main__":
    main()
