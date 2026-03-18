#!/usr/bin/env python3
"""Build teacher datasets for rigid adsorbates on Ag(111)."""

from __future__ import annotations

import argparse
from copy import deepcopy
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.append(str(ROOT))

from pyBall.gridff_jax import RunConfig, backend_summary, save_config
from pyBall.gridff_jax.config import AdsorbateConfig
from pyBall.gridff_jax.density_backends import make_density_backend
from pyBall.gridff_jax.fit import build_teacher_dataset, save_pose_dataset
from pyBall.gridff_jax.teacher_backends import make_teacher_backend
from pyBall.gridff_jax.utils import ensure_dir, save_json
from pyBall.gridff_jax.validation.volumetric import validate_vasp_volumetrics


DEFAULT_CHGCAR = "/home/niel/git/ORR_HER_Ag_Colab/results/Ag_ORR_HER/slab_clean/final_scf_12x12x1/CHGCAR"
DEFAULT_LOCPOT = "/home/niel/git/ORR_HER_Ag_Colab/results/Ag_ORR_HER/slab_clean/workfunc_12x12x1/LOCPOT"
DEFAULT_MADSURF_MODEL = str((Path(__file__).resolve().parent / "mad-surf_data" / "models" / "full_dataset_config_weights" / "MACE_model.model"))


def parse_args():
    parser = argparse.ArgumentParser(description="Build Ag(111) rigid-pose datasets for hybrid GridFF fitting")
    parser.add_argument("--chgcar", default=DEFAULT_CHGCAR, help="Path to Ag CHGCAR")
    parser.add_argument("--locpot", default=DEFAULT_LOCPOT, help="Path to Ag LOCPOT")
    parser.add_argument("--out-dir", default=str(ROOT / "tests/tGridFFJax/ag_dataset"), help="Dataset output directory")
    parser.add_argument("--teacher", default="madsurf", choices=["synthetic", "madsurf"], help="Teacher backend")
    parser.add_argument("--model-path", default=DEFAULT_MADSURF_MODEL, help="Path to MAD-SURF/MACE model")
    parser.add_argument("--device", default="cpu", help="Teacher device")
    parser.add_argument("--adsorbates", default="H,CO,H2O", help="Comma-separated adsorbate list")
    parser.add_argument("--grid-shape", default="", help="Optional grid shape override as nx,ny,nz")
    parser.add_argument("--alpha-morse", type=float, default=None, help="Override substrate Morse alpha used for Ag P/L fields")
    parser.add_argument("--ag-radius-scale", type=float, default=None, help="Scale Ag substrate vdW radius in the PL builder")
    parser.add_argument("--ag-energy-scale", type=float, default=None, help="Scale Ag substrate vdW energy in the PL builder")
    parser.add_argument("--samples-per-site", type=int, default=8, help="Systematic z samples per surface site")
    parser.add_argument("--systematic-orientations", type=int, default=4, help="Deterministic orientations per systematic site/z ladder")
    parser.add_argument("--random-samples", type=int, default=96, help="Random SE(3) samples per adsorbate")
    parser.add_argument("--representative-sites-per-label", type=int, default=1, help="Representative top/bridge/hollow sites per label")
    parser.add_argument("--z-bias-power", type=float, default=2.0, help="Bias systematic/random z sampling toward low heights when >1")
    parser.add_argument("--jitter-uv", type=float, default=0.03, help="Gaussian jitter in fractional uv for systematic site samples")
    parser.add_argument("--fit-static-charge", action="store_true", help="Allow fitting molecule-side static charges in the strict PLQ stage")
    parser.add_argument("--seed", type=int, default=12345, help="Random seed")
    return parser.parse_args()


def main():
    args = parse_args()
    out_dir = ensure_dir(args.out_dir)
    dataset_dir = ensure_dir(out_dir / "datasets")

    config = RunConfig()
    config.density_backend.chgcar_path = args.chgcar
    config.density_backend.locpot_path = args.locpot
    if args.grid_shape:
        config.density_backend.grid_shape = tuple(int(item) for item in args.grid_shape.split(","))
    if args.alpha_morse is not None:
        config.grid.alpha_morse = args.alpha_morse
    if args.ag_radius_scale is not None:
        config.grid.req_scale_radius["Ag"] = float(args.ag_radius_scale)
    if args.ag_energy_scale is not None:
        config.grid.req_scale_energy["Ag"] = float(args.ag_energy_scale)
    config.teacher_backend.kind = args.teacher
    config.teacher_backend.model_path = args.model_path
    config.teacher_backend.device = args.device
    config.grid.builder_mode = "metal_density_plq"
    config.toggles.use_ct_qeq = False
    config.toggles.use_image_charge = False
    config.toggles.use_reactive_grid = False
    config.hybrid_model.use_qeq = False
    config.hybrid_model.use_image = False
    config.hybrid_model.use_reactive_grid = False
    config.training.fit_static_charge = bool(args.fit_static_charge)
    config.training.fit_sample_shift_z = False
    config.training.fit_coulomb_shift_z = False
    config.export.out_dir = str(out_dir)
    selected = [name.strip() for name in args.adsorbates.split(",") if name.strip()]
    existing = {item.name: item for item in config.adsorbates}
    config.adsorbates = [existing.get(name, AdsorbateConfig(name=name, z_range=(1.2, 5.0))) for name in selected]
    for item in config.adsorbates:
        item.samples_per_site = args.samples_per_site
        item.systematic_orientations = args.systematic_orientations
        item.random_samples = args.random_samples
        item.representative_sites_per_label = args.representative_sites_per_label
        item.z_bias_power = args.z_bias_power
        item.jitter_uv = args.jitter_uv

    density_backend = make_density_backend(config.density_backend, surface=config.surface, grid=config.grid)
    density = density_backend.load()
    teacher_backend = make_teacher_backend(config.teacher_backend)

    save_config(config, out_dir / "config.json")
    save_json(
        {
            "teacher_backend": config.teacher_backend.kind,
            "teacher_device": config.teacher_backend.device,
            "student_backend": backend_summary(),
        },
        out_dir / "teacher.json",
    )
    save_json(validate_vasp_volumetrics(density), out_dir / "volumetric_validation.json")

    dataset_summary = {
        "teacher_backend": config.teacher_backend.kind,
        "teacher_device": config.teacher_backend.device,
        "adsorbates": {},
    }
    for offset, adsorbate_cfg in enumerate(config.adsorbates):
        local_cfg = deepcopy(adsorbate_cfg)
        poses, teacher = build_teacher_dataset(
            density=density,
            adsorbate_config=local_cfg,
            teacher_backend=teacher_backend,
            seed=args.seed + offset,
        )
        save_pose_dataset(dataset_dir, poses, teacher, density_metadata=density.metadata)
        dataset_summary["adsorbates"][adsorbate_cfg.name] = {
            "n_samples": int(len(poses.positions)),
            "n_sites_total": int(poses.metadata["n_sites_total"]),
            "n_sites_used": int(poses.metadata["n_sites_used"]),
            "teacher_seconds_total": float(teacher.metadata.get("teacher_eval_seconds", 0.0)),
            "teacher_seconds_per_pose": float(teacher.metadata.get("seconds_per_pose", 0.0)),
        }
        print(f"saved {adsorbate_cfg.name}: n={len(poses.positions)} -> {dataset_dir / (adsorbate_cfg.name + '.npz')}")
    dataset_summary["n_total"] = int(sum(item["n_samples"] for item in dataset_summary["adsorbates"].values()))
    save_json(dataset_summary, out_dir / "dataset_summary.json")


if __name__ == "__main__":
    main()
