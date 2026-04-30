#!/usr/bin/env python3
"""Generate 6D reference datasets for one or more molecules on a surface.

Reads molecules from a text file (or uses built-in names) and generates
systematic 6D (u,v,z,alpha,beta,gamma) pose batches, evaluates them with
a teacher backend (MAD-SURF by default), and saves the resulting datasets.

Usage
-----
    # Single molecule (built-in):
    python run_6d_sampling.py --molecules CHONH2

    # Multiple molecules from a file:
    python run_6d_sampling.py --molecule-file molecules.txt

    # Custom grid density:
    python run_6d_sampling.py --molecules H,CO,H2O --n-u 12 --n-v 12 --n-z 25

    # With specific substrate data:
    python run_6d_sampling.py --molecules CHONH2 \
        --chgcar /path/to/CHGCAR --locpot /path/to/LOCPOT \
        --model-path /path/to/MACE_model.model
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.append(str(ROOT))

from pyBall.gridff_jax import RunConfig, backend_summary, save_config
from pyBall.gridff_jax.pbc import PBCInfo, check_molecule_cell_compatibility, compute_minimum_teacher_tiling
from pyBall.gridff_jax.density_backends import make_density_backend
from pyBall.gridff_jax.fit.fit_6d import (
    Dataset6D,
    generate_6d_reference_data,
    save_6d_dataset,
)
from pyBall.gridff_jax.pose_sampling import get_adsorbate, load_molecule_list
from pyBall.gridff_jax.pose_sampling.sampler_6d import Sampler6DConfig
from pyBall.gridff_jax.teacher_backends import make_teacher_backend
from pyBall.gridff_jax.utils import ensure_dir


DEFAULT_CHGCAR = "/home/niel/git/ORR_HER_Ag_Colab/results/Ag_ORR_HER/slab_clean/final_scf_12x12x1/CHGCAR"
DEFAULT_LOCPOT = "/home/niel/git/ORR_HER_Ag_Colab/results/Ag_ORR_HER/slab_clean/workfunc_12x12x1/LOCPOT"
DEFAULT_MADSURF_MODEL = str(
    Path(__file__).resolve().parent
    / "mad-surf_data"
    / "models"
    / "full_dataset_config_weights"
    / "MACE_model.model"
)


def parse_args():
    p = argparse.ArgumentParser(
        description="Generate 6D reference datasets for molecule-on-surface systems",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Substrate
    p.add_argument("--chgcar", default=DEFAULT_CHGCAR, help="Path to substrate CHGCAR")
    p.add_argument("--locpot", default=DEFAULT_LOCPOT, help="Path to substrate LOCPOT")
    p.add_argument("--model-path", default=DEFAULT_MADSURF_MODEL, help="Path to teacher MLIP model")
    p.add_argument("--device", default="cuda", help="Teacher device (cuda/cpu)")
    p.add_argument("--substrate", default="Ag", help="Substrate element")
    p.add_argument("--facet", default="111", help="Surface facet")

    # Molecules
    p.add_argument("--molecules", default="CHONH2",
                   help="Comma-separated built-in molecule names (H, CO, H2O, CHONH2, ...)")
    p.add_argument("--molecule-file", default=None,
                   help="Path to a text file listing molecules (one per line)")

    # 6D sampler configuration
    p.add_argument("--n-u", type=int, default=10, help="Lateral grid points along u")
    p.add_argument("--n-v", type=int, default=10, help="Lateral grid points along v")
    p.add_argument("--n-z", type=int, default=20, help="Height points per lateral position")
    p.add_argument("--z-min", type=float, default=1.5, help="Minimum z above surface (A)")
    p.add_argument("--z-max", type=float, default=5.5, help="Maximum z above surface (A)")
    p.add_argument("--z-bias", type=float, default=1.5, help="Z sampling bias power")
    p.add_argument("--n-orient", type=int, default=8, help="Number of orientations")
    p.add_argument("--tilt-max", type=float, default=60.0, help="Max tilt from surface normal (deg)")
    p.add_argument("--n-yaw", type=int, default=4, help="Yaw angles per tilt direction")
    p.add_argument("--random-fraction", type=float, default=0.1, help="Fraction of random poses")
    p.add_argument("--seed", type=int, default=42, help="Random seed")

    # Teacher configuration
    p.add_argument("--teacher-chunk", type=int, default=64, help="Teacher batch chunk size")
    p.add_argument("--teacher-tile", type=str, default="1,1",
                   help="Slab tiling for teacher: 'NX,NY'")

    # Grid building (for density backend)
    p.add_argument("--london-damp-d0", type=float, default=3.70, help="London damping d0")
    p.add_argument("--london-damp-width", type=float, default=0.35, help="London damping width")

    # Output
    p.add_argument("--out-dir", default=None,
                   help="Output directory (default: tests/tGridFFJax/datasets_6d/)")

    return p.parse_args()


def main():
    args = parse_args()

    out_dir = Path(args.out_dir) if args.out_dir else (
        Path(__file__).resolve().parent / "datasets_6d"
    )
    ensure_dir(out_dir)

    print("=" * 60)
    print("  6D Reference Dataset Generation")
    print("=" * 60)

    # --- Build sampler config ---
    sampler_config = Sampler6DConfig(
        n_u=args.n_u,
        n_v=args.n_v,
        n_z=args.n_z,
        z_range=(args.z_min, args.z_max),
        z_bias_power=args.z_bias,
        n_orient=args.n_orient,
        tilt_max_deg=args.tilt_max,
        n_yaw=args.n_yaw,
        random_fraction=args.random_fraction,
        seed=args.seed,
    )

    n_est = args.n_u * args.n_v * args.n_z * args.n_orient
    n_random_est = int(n_est * args.random_fraction)
    print(f"  Grid: {args.n_u}x{args.n_v} lateral, {args.n_z} z-points, {args.n_orient} orientations")
    print(f"  Estimated poses per molecule: ~{n_est + n_random_est}")

    # --- Load molecules ---
    if args.molecule_file:
        print(f"  Loading molecules from: {args.molecule_file}")
        molecules = load_molecule_list(args.molecule_file)
    else:
        names = [n.strip() for n in args.molecules.split(",") if n.strip()]
        print(f"  Molecules: {names}")
        molecules = [get_adsorbate(name=n) for n in names]

    print(f"  Total molecules: {len(molecules)}")
    for mol in molecules:
        print(f"    - {mol.name}: {mol.natoms} atoms ({', '.join(mol.symbols)})")

    # --- Build RunConfig for density/teacher ---
    config = RunConfig()
    config.density_backend.chgcar_path = args.chgcar
    config.density_backend.locpot_path = args.locpot
    config.surface.metal = args.substrate
    config.surface.facet = args.facet
    config.grid.london_damping_d0 = args.london_damp_d0
    config.grid.london_damping_width = args.london_damp_width

    tile_parts = args.teacher_tile.split(",")
    config.teacher_backend.teacher_tile = (int(tile_parts[0]), int(tile_parts[1]))
    config.teacher_backend.model_path = args.model_path
    config.teacher_backend.device = args.device
    config.teacher_backend.kind = "madsurf"

    # --- Load density ---
    print("\n[sampling] Loading substrate density ...", flush=True)
    density = make_density_backend(
        config.density_backend, surface=config.surface, grid=config.grid
    ).load()
    print(f"  Cell: {density.cell[0,0]:.2f} x {density.cell[1,1]:.2f} A")
    print(f"  Atoms: {density.natoms}")

    # --- PBC diagnostics ---
    pbc_info = PBCInfo.from_density(density)
    print(f"\n[PBC] Cell: {pbc_info.cell_lengths[0]:.2f} x {pbc_info.cell_lengths[1]:.2f} x {pbc_info.cell_lengths[2]:.2f} A")
    print(f"[PBC] Lateral area: {pbc_info.lateral_area:.1f} A^2")
    print(f"[PBC] Supercell multiplicity: {pbc_info.supercell_multiplicity[0]}x{pbc_info.supercell_multiplicity[1]}")

    max_tile = (1, 1)
    for mol in molecules:
        ok, msg = check_molecule_cell_compatibility(mol, density)
        tile = compute_minimum_teacher_tiling(mol, density.cell)
        max_tile = (max(max_tile[0], tile[0]), max(max_tile[1], tile[1]))
        status = "OK" if ok else "WARNING"
        print(f"[PBC] {mol.name}: {status} — {msg}")
        if tile != (1, 1):
            print(f"       Min teacher tiling: {tile[0]}x{tile[1]}")

    # Auto-set teacher tiling if user specified "auto" (0,0) or if needed
    if config.teacher_backend.teacher_tile == (0, 0) or config.teacher_backend.teacher_tile == (1, 1):
        if max_tile != (1, 1):
            config.teacher_backend.teacher_tile = max_tile
            print(f"[PBC] Auto-setting teacher tile to {max_tile[0]}x{max_tile[1]}")

    # --- Load teacher ---
    print("\n[sampling] Loading teacher backend ...", flush=True)
    teacher_backend = make_teacher_backend(config.teacher_backend)

    # --- Generate datasets ---
    t_total_start = time.perf_counter()
    for mol in molecules:
        print(f"\n{'='*40}")
        print(f"  Molecule: {mol.name}")
        print(f"{'='*40}")

        t0 = time.perf_counter()
        dataset = generate_6d_reference_data(
            density=density,
            adsorbate=mol,
            teacher_backend=teacher_backend,
            sampler_config=sampler_config,
            chunk_size=args.teacher_chunk,
            verbose=True,
        )
        elapsed = time.perf_counter() - t0

        save_6d_dataset(dataset, out_dir)

        e = dataset.energies
        print(f"\n  {mol.name}: {dataset.n_poses} poses in {elapsed:.1f}s")
        print(f"  E range: [{e.min()*1000:.1f}, {e.max()*1000:.1f}] meV")
        print(f"  E mean: {e.mean()*1000:.1f} meV, std: {e.std()*1000:.1f} meV")
        print(f"  Saved to: {out_dir / (mol.name + '_6d.npz')}")

    t_total = time.perf_counter() - t_total_start

    # --- Save summary ---
    summary = {
        "substrate": args.substrate,
        "facet": args.facet,
        "molecules": [mol.name for mol in molecules],
        "sampler_config": sampler_config.__dict__,
        "total_time_seconds": t_total,
    }
    summary_path = out_dir / "sampling_summary.json"
    with summary_path.open("w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, default=str)

    # --- Save config ---
    save_config(config, out_dir / "config.json")

    print(f"\n{'='*60}")
    print(f"  6D sampling complete for {len(molecules)} molecule(s)")
    print(f"  Total time: {t_total:.1f}s")
    print(f"  Output: {out_dir}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
