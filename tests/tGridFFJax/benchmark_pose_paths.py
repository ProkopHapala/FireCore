#!/usr/bin/env python3
"""Benchmark smoothness and throughput of MAD-SURF vs distilled GridFF along continuous pose paths."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.append(str(ROOT))

from pyBall.gridff_jax import load_config
from pyBall.gridff_jax.benchmarking import benchmark_pose_paths, load_fit_params
from pyBall.gridff_jax.density_backends import make_density_backend
from pyBall.gridff_jax.hybrid_energy import HybridGridFFModel
from pyBall.gridff_jax.pose_sampling import get_adsorbate
from pyBall.gridff_jax.teacher_backends import make_teacher_backend
from pyBall.gridff_jax.utils import ensure_dir


def parse_args():
    base = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description="Run path smoothness and timing benchmarks for a fitted GridFF model")
    parser.add_argument("--config", default=str(base / "ag_fit" / "config_used.json"))
    parser.add_argument("--fit-params", default=str(base / "ag_fit" / "fit_params.json"))
    parser.add_argument("--adsorbates", default="H,CO,H2O")
    parser.add_argument("--out-dir", default=str(base / "path_benchmarks"))
    parser.add_argument("--prefer-jax", dest="prefer_jax", action="store_true", help="Use the JAX student backend when available")
    parser.add_argument("--no-prefer-jax", dest="prefer_jax", action="store_false", help="Force the NumPy/SciPy fallback backend")
    parser.set_defaults(prefer_jax=True)
    return parser.parse_args()


def main():
    args = parse_args()
    config = load_config(args.config)
    density = make_density_backend(config.density_backend, surface=config.surface, grid=config.grid).load()
    teacher = make_teacher_backend(config.teacher_backend)
    params = load_fit_params(args.fit_params)
    model = HybridGridFFModel(
        density=density,
        reactive_elements=tuple(sorted(config.hybrid_model.reactive_channels)),
        toggles=config.toggles,
        grid_config=config.grid,
        hybrid_config=config.hybrid_model,
        prefer_jax=args.prefer_jax,
    )
    out_dir = ensure_dir(args.out_dir)
    from pyBall.gridff_jax.pose_sampling.sites import infer_surface_sites
    from pyBall.gridff_jax.pose_sampling.rigid import infer_reference_planes

    pose_meta = infer_reference_planes(density)
    site_lookup = {}
    for site in infer_surface_sites(density):
        if site.label not in site_lookup:
            site_lookup[site.label] = site.uv.tolist()
    pose_meta["surface_sites"] = [{"label": label, "uv": uv} for label, uv in sorted(site_lookup.items())]
    for name in [item.strip() for item in args.adsorbates.split(",") if item.strip()]:
        adsorbate = get_adsorbate(name)
        benchmark_pose_paths(
            teacher_backend=teacher,
            model=model,
            params=params,
            density=density,
            adsorbate=adsorbate,
            pose_metadata=pose_meta,
            out_dir=ensure_dir(out_dir / name),
        )
    print(f"path benchmarks -> {out_dir}")


if __name__ == "__main__":
    main()
