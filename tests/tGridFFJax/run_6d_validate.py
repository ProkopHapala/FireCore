#!/usr/bin/env python3
"""Validate 6D-fitted parameters with lateral + z scans and held-out data.

Loads fitted parameters from ``run_6d_fit.py`` and evaluates them on:
  1. The held-out test set from the 6D split
  2. Dense lateral scans (1D path + 2D heatmap) at multiple z-heights
  3. Dense z-scans at high-symmetry sites
  4. Optionally, transferability test on molecules NOT in the training set

Usage
-----
    python run_6d_validate.py --fit-dir datasets_6d/fit_results/
    python run_6d_validate.py --fit-dir datasets_6d/fit_results/ \
        --test-molecules CH2O,C6H6 --z-heights 2.0,2.5,3.0
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

from pyBall.gridff_jax import RunConfig, load_config
from pyBall.gridff_jax.density_backends import make_density_backend
from pyBall.gridff_jax.diagnostics.comparison import MethodResult, compare_methods, parity_statistics
from pyBall.gridff_jax.hybrid_energy import HybridGridFFModel
from pyBall.gridff_jax.interfaces import PoseBatch, TeacherResult
from pyBall.gridff_jax.pose_sampling import get_adsorbate, infer_reference_planes, infer_surface_sites, transform_adsorbate
from pyBall.gridff_jax.teacher_backends import make_teacher_backend
from pyBall.gridff_jax.utils import ensure_dir, save_json


def parse_args():
    p = argparse.ArgumentParser(
        description="Validate 6D-fitted parameters with lateral and z-scans",
    )
    p.add_argument("--fit-dir", required=True,
                   help="Directory with fit_params.json and config.json")
    p.add_argument("--dataset-dir", default=None,
                   help="Directory with *_6d.npz (default: parent of fit-dir)")
    p.add_argument("--molecules", default="CHONH2",
                   help="Molecules to validate (comma-separated)")
    p.add_argument("--test-molecules", default=None,
                   help="Additional held-out molecules for transferability test")
    p.add_argument("--z-heights", default="2.0,2.5,3.0,3.5",
                   help="Z-heights for lateral scans")
    p.add_argument("--n-lateral", type=int, default=25,
                   help="Grid points per dimension for 2D lateral scan")
    p.add_argument("--n-zscan", type=int, default=60,
                   help="Points per z-scan")
    p.add_argument("--z-scan-range", default="1.8,5.5",
                   help="Z-scan range: z_min,z_max")
    p.add_argument("--device", default="cuda")
    p.add_argument("--out-dir", default=None)
    p.add_argument("--prefer-jax", dest="prefer_jax", action="store_true", default=True)
    p.add_argument("--no-prefer-jax", dest="prefer_jax", action="store_false")
    return p.parse_args()


def _build_zscan_poses(density, adsorbate, site_uvs, z_values, quaternion):
    planes = infer_reference_planes(density)
    z_ref = float(planes["z_ref"])
    positions, params, labels = [], [], []
    for label, uv in site_uvs.items():
        for z_h in z_values:
            pos = transform_adsorbate(adsorbate, density, uv, float(z_h), quaternion, z_ref)
            positions.append(pos)
            params.append(np.concatenate([uv, [float(z_h)], quaternion]))
            labels.append(label)
    return PoseBatch(
        adsorbate=adsorbate,
        positions=np.asarray(positions, dtype=float),
        pose_params=np.asarray(params, dtype=float),
        site_labels=labels,
        metadata={"z_ref": z_ref},
    )


def _build_lateral_poses(density, adsorbate, n_grid, z_height, quaternion):
    planes = infer_reference_planes(density)
    z_ref = float(planes["z_ref"])
    u_vals = np.linspace(0.0, 1.0, n_grid, endpoint=False)
    v_vals = np.linspace(0.0, 1.0, n_grid, endpoint=False)
    uu, vv = np.meshgrid(u_vals, v_vals, indexing="ij")
    uv_flat = np.column_stack([uu.ravel(), vv.ravel()])
    positions, params = [], []
    for uv in uv_flat:
        pos = transform_adsorbate(adsorbate, density, uv, z_height, quaternion, z_ref)
        positions.append(pos)
        params.append(np.concatenate([uv, [z_height], quaternion]))
    return PoseBatch(
        adsorbate=adsorbate,
        positions=np.asarray(positions, dtype=float),
        pose_params=np.asarray(params, dtype=float),
        site_labels=["grid"] * len(uv_flat),
        metadata={"z_ref": z_ref, "z_height": z_height},
    ), u_vals, v_vals


def _eval_chunked(fn, poses, chunk=64):
    n = len(poses.positions)
    all_e, all_f = [], []
    for s in range(0, n, chunk):
        e = min(s + chunk, n)
        idx = np.arange(s, e)
        partial = PoseBatch(
            adsorbate=poses.adsorbate,
            positions=poses.positions[idx],
            pose_params=poses.pose_params[idx],
            site_labels=[poses.site_labels[int(i)] for i in idx],
            metadata=dict(poses.metadata),
        )
        r = fn(partial)
        all_e.append(np.asarray(r.energies, dtype=float))
        if r.forces is not None:
            all_f.append(np.asarray(r.forces, dtype=float))
    energies = np.concatenate(all_e, axis=0)
    forces = np.concatenate(all_f, axis=0) if all_f else None
    return type("R", (), {"energies": energies, "forces": forces})()


def main():
    args = parse_args()
    fit_dir = Path(args.fit_dir)
    dataset_dir = Path(args.dataset_dir) if args.dataset_dir else fit_dir.parent
    out_dir = Path(args.out_dir) if args.out_dir else fit_dir / "validation"
    ensure_dir(out_dir)

    mol_names = [n.strip() for n in args.molecules.split(",") if n.strip()]
    z_heights = [float(z) for z in args.z_heights.split(",")]
    z_lo, z_hi = [float(x) for x in args.z_scan_range.split(",")]

    print("=" * 60)
    print("  6D Validation Suite")
    print("=" * 60)

    # Load config + params
    config = load_config(fit_dir / "config.json")
    config.teacher_backend.device = args.device

    with (fit_dir / "fit_params.json").open("r") as fh:
        from pyBall.gridff_jax.hybrid_energy.model import HybridParameters
        params_raw = json.load(fh)
        params = HybridParameters(**{
            k: v for k, v in params_raw.items()
            if k in HybridParameters.__dataclass_fields__
        })

    print("[val] Loading density ...", flush=True)
    density = make_density_backend(
        config.density_backend, surface=config.surface, grid=config.grid
    ).load()

    print("[val] Loading teacher ...", flush=True)
    teacher_backend = make_teacher_backend(config.teacher_backend)

    # Site lookup
    sites = infer_surface_sites(density)
    site_uvs = {}
    for s in sites:
        if s.label not in site_uvs:
            site_uvs[s.label] = np.asarray(s.uv, dtype=float)

    quaternion = np.array([1.0, 0.0, 0.0, 0.0], dtype=float)

    all_metrics = {}

    for mol_name in mol_names:
        mol_out = ensure_dir(out_dir / mol_name)
        print(f"\n{'='*50}")
        print(f"  Validating: {mol_name}")
        print(f"{'='*50}")

        adsorbate = get_adsorbate(name=mol_name)
        all_elements = tuple(sorted(set(adsorbate.symbols)))

        model = HybridGridFFModel(
            density=density,
            reactive_elements=all_elements,
            toggles=config.toggles,
            grid_config=config.grid,
            hybrid_config=config.hybrid_model,
            prefer_jax=args.prefer_jax,
        )
        model.install_raw_r6_grid()

        def gridff_fn(poses):
            return model.evaluate_batch(poses, params=params, compute_forces=True)

        def teacher_fn(poses):
            return teacher_backend.evaluate_batch(density, poses)

        mol_metrics = {"molecule": mol_name}

        # --- Z-scans at high-symmetry sites ---
        print("\n  --- Z-scans ---")
        z_vals = np.linspace(z_lo, z_hi, args.n_zscan)
        zscan_poses = _build_zscan_poses(density, adsorbate, site_uvs, z_vals, quaternion)

        teacher_z = _eval_chunked(teacher_fn, zscan_poses)
        gridff_z = _eval_chunked(gridff_fn, zscan_poses)

        z_stats = parity_statistics(teacher_z.energies, gridff_z.energies,
                                    teacher_z.forces, gridff_z.forces)
        mol_metrics["zscan"] = z_stats
        print(f"  Z-scan: E_RMSE={z_stats['energy_rmse_meV']:.1f} meV, "
              f"F_RMSE={z_stats.get('force_rmse_eV_A', 0):.4f} eV/A")

        np.savez_compressed(
            mol_out / "zscan_data.npz",
            z_values=z_vals,
            teacher_energies=teacher_z.energies,
            gridff_energies=gridff_z.energies,
            site_labels=zscan_poses.site_labels,
        )

        # --- 2D lateral scans at each z-height ---
        print("\n  --- Lateral scans ---")
        for z_h in z_heights:
            print(f"  z = {z_h:.2f} A ...", end=" ", flush=True)
            lat_poses, u_vals, v_vals = _build_lateral_poses(
                density, adsorbate, args.n_lateral, z_h, quaternion,
            )

            teacher_lat = _eval_chunked(teacher_fn, lat_poses)
            gridff_lat = _eval_chunked(gridff_fn, lat_poses)

            lat_stats = parity_statistics(teacher_lat.energies, gridff_lat.energies)
            print(f"E_RMSE={lat_stats['energy_rmse_meV']:.1f} meV, "
                  f"corrugation: teacher={float(teacher_lat.energies.max()-teacher_lat.energies.min())*1000:.0f} meV, "
                  f"GridFF={float(gridff_lat.energies.max()-gridff_lat.energies.min())*1000:.0f} meV")

            mol_metrics[f"lateral_z{z_h:.1f}"] = lat_stats

            n_g = args.n_lateral
            np.savez_compressed(
                mol_out / f"lateral_z{z_h:.1f}.npz",
                u_vals=u_vals, v_vals=v_vals,
                teacher_energy_2d=teacher_lat.energies.reshape(n_g, n_g),
                gridff_energy_2d=gridff_lat.energies.reshape(n_g, n_g),
                z_height=z_h,
            )

        all_metrics[mol_name] = mol_metrics

    save_json(all_metrics, out_dir / "validation_metrics.json")

    print(f"\n{'='*60}")
    print(f"  Validation complete for {len(mol_names)} molecule(s)")
    print(f"  Output: {out_dir}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
