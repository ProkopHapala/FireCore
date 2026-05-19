#!/usr/bin/env python3
"""Generate a NaCl parity-core PLQ grid and optionally compare it to a reference export."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.append(str(ROOT))

from pyBall.gridff_jax.config import DensityBackendConfig, FeatureToggles, GridConfig
from pyBall.gridff_jax.density_backends import make_density_backend
from pyBall.gridff_jax.hybrid_energy import HybridGridFFModel
from pyBall.gridff_jax.utils import ensure_dir, save_json


DEFAULT_XYZ = str(ROOT / "cpp/common_resources/xyz/NaCl_1x1_L3.xyz")


def parse_args():
    parser = argparse.ArgumentParser(description="Build a NaCl parity-core PLQ grid")
    parser.add_argument("--xyz", default=DEFAULT_XYZ, help="Path to substrate xyz")
    parser.add_argument("--out-dir", default=str(ROOT / "tests/tGridFFJax/nacl_parity"), help="Output directory")
    parser.add_argument("--desired-voxel", type=float, default=0.1, help="Grid voxel size in Angstrom")
    parser.add_argument("--reference-grid", default="", help="Optional reference Bspline_PLQd.npy for comparison")
    parser.add_argument("--prefer-jax", dest="prefer_jax", action="store_true")
    parser.add_argument("--no-prefer-jax", dest="prefer_jax", action="store_false")
    parser.set_defaults(prefer_jax=True)
    return parser.parse_args()


def _error_metrics(left, right):
    delta = np.asarray(left, dtype=float) - np.asarray(right, dtype=float)
    return {
        "rmse": float(np.sqrt(np.mean(delta * delta))),
        "mae": float(np.mean(np.abs(delta))),
        "max_abs": float(np.max(np.abs(delta))),
    }


def main():
    args = parse_args()
    out_dir = ensure_dir(args.out_dir)
    density_cfg = DensityBackendConfig(kind="surface_xyz", xyz_path=args.xyz)
    grid_cfg = GridConfig(
        spacing=args.desired_voxel,
        builder_mode="parity_core",
        interpolation_order=3,
    )
    # Keyword construction only — positional FeatureToggles silently maps the
    # 6th positional arg onto use_reactive_grid (dataclass field #6) which
    # ENABLES the reactive grid instead of the intended use_teacher_residual.
    toggles = FeatureToggles(
        use_density_charge=True, use_locpot=True,
        use_ct_qeq=False, use_image_charge=False,
        use_image_charge_fixed=False, use_reactive_grid=False,
        use_teacher_residual=True,
    )
    density = make_density_backend(density_cfg, grid=grid_cfg).load()
    model = HybridGridFFModel(
        density=density,
        reactive_elements=tuple(sorted(set(density.symbols))),
        toggles=toggles,
        grid_config=grid_cfg,
        prefer_jax=args.prefer_jax,
    )

    plq = np.stack(
        [
            np.asarray(model.grids["pauli_coeff_zyx"], dtype=float),
            np.asarray(model.grids["london_coeff_zyx"], dtype=float),
            np.asarray(model.grids["coulomb_coeff_zyx"], dtype=float),
        ],
        axis=-1,
    ).transpose((2, 1, 0, 3))
    np.save(out_dir / "Bspline_PLQd.npy", plq)

    summary = {
        "xyz_path": str(args.xyz),
        "grid_shape_xyz": list(density.grid_shape_xyz),
        "origin": density.origin.tolist(),
        "voxel": density.voxel.tolist(),
        "builder_metadata": model.grids["metadata"],
        "model_use_jax": bool(model.use_jax),
    }
    if args.reference_grid:
        reference = np.load(args.reference_grid)
        if reference.shape != plq.shape:
            raise ValueError(f"Reference grid shape {reference.shape} does not match {plq.shape}")
        summary["reference_grid"] = str(args.reference_grid)
        summary["channel_errors"] = {
            "pauli": _error_metrics(plq[..., 0], reference[..., 0]),
            "london": _error_metrics(plq[..., 1], reference[..., 1]),
            "coulomb": _error_metrics(plq[..., 2], reference[..., 2]),
        }
    save_json(summary, out_dir / "parity_summary.json")
    print(f"nacl parity grid -> {out_dir}")


if __name__ == "__main__":
    main()
