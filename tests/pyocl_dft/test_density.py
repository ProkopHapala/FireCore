#!/usr/bin/env python3
"""Milestone B: GPU density projection validation using PyOpenCL.

This script bridges Fireball SCF data with the experimental PyOpenCL host.
It loads precomputed molecular-orbital coefficients and density grids created
by :mod:`tests.pyocl_dft.test_firecore_data`, runs the OpenCL texture-based
projection, and compares the GPU output against the reference density.

Outputs per run:
- ``*_gpu_density.npy`` – projected complex grid (float32 complex, stored as
  two channels) produced by the PyOpenCL pipeline.
- ``*_comparison_metrics.json`` – RMS / L∞ error metrics for diagnostics.
- ``*_comparison_plot.png`` – overlays of GPU vs reference slices.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Sequence

import matplotlib.pyplot as plt
import numpy as np
import pyopencl as cl

from pyBall import atomicUtils as au
from pyBall.DFT import utils as dft_utils
from pyBall.pyocl_dft import assets
from pyBall.pyocl_dft.density import OCLDensityProjector, projector_singleton

DEFAULT_PREFIX = "firecore"
DEFAULT_OUTPUT = Path("pyocl_dft_density")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="PyOpenCL density projection benchmark")
    parser.add_argument("--xyz",          type=Path, default=Path("tests/pyFireball/molecules/CH4.xyz"), help="XYZ file for atom ordering metadata (default matches milestone A)")
    parser.add_argument("--coeffs",       type=Path, default=DEFAULT_OUTPUT / f"{DEFAULT_PREFIX}_density_matrix.npy", help="Path to MO coefficient matrix saved by test_firecore_data.py")
    parser.add_argument("--density",      type=Path, default=DEFAULT_OUTPUT / f"{DEFAULT_PREFIX}_density_grid.npy", help="Reference real-space density grid (float32) for comparison")
    parser.add_argument("--grid-meta",    type=Path, default=DEFAULT_OUTPUT / f"{DEFAULT_PREFIX}_grid_meta.json", help="Grid metadata JSON produced by test_firecore_data.py")
    parser.add_argument("--basis-dir",    type=Path, default=None, help="Optional override for Fireball basis directory (Fdata/basis)")
    parser.add_argument("--output-dir",   type=Path, default=DEFAULT_OUTPUT, help="Directory to store projected grids and comparison plots")
    parser.add_argument("--prefix",       type=str, default=DEFAULT_PREFIX, help="Prefix for output artefacts")
    parser.add_argument("--plot-slice",   type=int, default=None, help="Optional index along Z for density slice (default: middle plane)")
    parser.add_argument("--show",         type=int, default=0, action="store_true", help="Display plots interactively instead of closing after save")
    parser.add_argument("--device-index", type=int, default=0, help="Device index forwarded to OpenCLBase (default: auto-select first NVIDIA)")
    parser.add_argument("--dry-run",      type=int, default=0, help="Skip kernel launch; useful for verifying data loading")
    return parser.parse_args()

def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def load_grid_meta(path: Path) -> dict:
    data = json.loads(path.read_text(encoding="utf8"))
    required = {"ngrid", "dcell", "lvs"}
    missing = required - data.keys()
    if missing:
        raise KeyError(f"Grid metadata file {path} missing keys: {missing}")
    return data


def as_float4(vec: Sequence[float]) -> np.ndarray:
    if len(vec) != 3:
        raise ValueError("Expected 3-vector for grid basis")
    return np.asarray(list(vec) + [0.0], dtype=np.float32)


def build_grid_vectors(meta: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    dcell = np.asarray(meta["dcell"], dtype=np.float32)
    if dcell.shape != (3, 3):
        raise ValueError(f"dcell must be 3x3, got {dcell.shape}")
    return tuple(as_float4(row) for row in dcell)


def make_basis_image(ctx: cl.Context, basis_dir: Path, elements: Sequence[int]) -> cl.Image:
    """Upload basis textures as a 2D image using legacy wf files.

    The existing C++ path loads radial basis functions from ``*.wf1``/``*.wf2``
    pairs. For parity we reuse the same naming convention here. The image stores
    ``float2`` values with slots for ``s`` and ``p`` harmonics.
    """

    # Discover basis files
    symbols = [f"{z:03d}" for z in sorted(set(int(z) for z in elements))]
    assets.verify_basis(basis_dir, symbols)  # validates files exist

    nsamp = 512
    nslots = len(symbols) * 2
    host = np.zeros((nsamp, nslots, 2), dtype=np.float32)

    for slot, sym in enumerate(symbols):
        for wf_idx, suffix in enumerate(("wf1", "wf2")):
            fname = basis_dir / f"{sym}_450.{suffix}"  # TODO: expose cutoffs in metadata
            data = np.fromfile(fname, dtype=np.float32, count=nsamp * 2)
            if data.size != nsamp * 2:
                raise RuntimeError(f"Unexpected samples in {fname}: {data.size}")
            host[:, slot * 2 + wf_idx, :] = data.reshape(nsamp, 2)

    fmt = cl.ImageFormat(cl.channel_order.RG, cl.channel_type.FLOAT)
    image = cl.Image(ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, fmt, host.shape[:2], hostbuf=host)
    return image


def load_atoms(xyz_path: Path) -> tuple[np.ndarray, np.ndarray]:
    apos, atom_types, _, _ = au.loadAtomsNP(xyz_path)
    return apos.astype(np.float64), np.asarray(atom_types, dtype=np.int32)


def prepare_coefficients(wfcoef: np.ndarray, atom_types: np.ndarray) -> np.ndarray:
    i0orb = dft_utils.countOrbs(atom_types)
    if wfcoef.shape[0] != wfcoef.shape[1]:
        raise ValueError("Wavefunction coefficient matrix must be square")
    if i0orb[-1] != wfcoef.shape[0]:
        raise ValueError("Counted orbitals do not match coefficient dimension")
    # Collapse into float4 blocks per atom as expected by the kernel
    coefs = []
    start = 0
    for stop in i0orb:
        block = wfcoef[start:stop, :].copy()
        block4 = np.zeros((block.shape[1], 4), dtype=np.float32)
        block4[:, : min(block.shape[0], 4)] = block[:4, :].T
        coefs.append(block4)
        start = stop
    return np.vstack(coefs)


def format_atoms(apos: np.ndarray, atom_types: np.ndarray) -> np.ndarray:
    atoms4 = np.zeros((len(atom_types), 4), dtype=np.float32)
    atoms4[:, :3] = apos.astype(np.float32)
    atoms4[:, 3] = atom_types.astype(np.float32)
    return atoms4


def compute_metrics(reference: np.ndarray, gpu: np.ndarray) -> dict:
    diff = gpu[..., 0] - reference
    rms = float(np.sqrt(np.mean(diff**2)))
    linf = float(np.max(np.abs(diff)))
    return {"rms": rms, "linf": linf}


def save_metrics(path: Path, metrics: dict) -> None:
    path.write_text(json.dumps(metrics, indent=2), encoding="utf8")


def plot_slice(path: Path, reference: np.ndarray, gpu: np.ndarray, iz: int | None) -> None:
    if iz is None:
        iz = reference.shape[2] // 2
    iz = max(0, min(int(iz), reference.shape[2] - 1))

    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    im0 = axes[0].imshow(reference[:, :, iz].T, origin="lower", cmap="viridis")
    axes[0].set_title("Reference density")
    fig.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04)

    im1 = axes[1].imshow(gpu[:, :, iz, 0].T, origin="lower", cmap="viridis")
    axes[1].set_title("GPU density (real)")
    fig.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04)

    diff = gpu[:, :, iz, 0] - reference[:, :, iz]
    im2 = axes[2].imshow(diff.T, origin="lower", cmap="coolwarm")
    axes[2].set_title("Difference (real)")
    fig.colorbar(im2, ax=axes[2], fraction=0.046, pad=0.04)

    for ax in axes:
        ax.set_xlabel("x index")
        ax.set_ylabel("y index")

    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    out_dir = ensure_dir(args.output_dir)

    apos, atom_types = load_atoms(args.xyz)
    wfcoef = np.load(args.coeffs)
    density_ref = np.load(args.density)
    grid_meta = load_grid_meta(args.grid_meta)

    grid_shape = tuple(int(v) for v in grid_meta["ngrid"])
    if density_ref.shape != tuple(grid_shape):
        raise ValueError(f"Reference density shape {density_ref.shape} != grid shape {grid_shape}")

    basis_dir = args.basis_dir or (Path(__file__).resolve().parents[2] / "Fdata" / "basis")
    if not basis_dir.is_dir():
        raise FileNotFoundError(f"Basis directory not found: {basis_dir}")

    projector = projector_singleton()
    atoms4 = format_atoms(apos, atom_types)
    coefs4 = prepare_coefficients(wfcoef, atom_types)
    bufs = projector.ensure_buffers(atoms=atoms4, coefs=coefs4, grid_shape=grid_shape)

    if args.dry_run:
        print("Dry run completed: buffers prepared, skipping kernel launch")
        return

    basis_image = make_basis_image(projector.ctx, basis_dir, atom_types)
    grid_origin = np.asarray(grid_meta.get("pos0", [0.0, 0.0, 0.0]), dtype=np.float32)
    grid_vectors = build_grid_vectors(grid_meta)
    acum_coef = assets.ensure_acum_coef()

    gpu_density = projector.project_density(
        bufs,
        basis_image=basis_image,
        grid_origin=grid_origin,
        grid_vectors=grid_vectors,
        acum_coef=acum_coef,
    )

    metrics = compute_metrics(density_ref, gpu_density)
    metrics_path = out_dir / f"{args.prefix}_comparison_metrics.json"
    save_metrics(metrics_path, metrics)

    gpu_path = out_dir / f"{args.prefix}_gpu_density.npy"
    np.save(gpu_path, gpu_density)

    plot_path = out_dir / f"{args.prefix}_comparison_plot.png"
    plot_slice(plot_path, density_ref, gpu_density, args.plot_slice)

    print(f"Stored metrics in {metrics_path}")
    print(f"Stored GPU density in {gpu_path}")
    print(f"Stored comparison plot in {plot_path}")


if __name__ == "__main__":
    main()
