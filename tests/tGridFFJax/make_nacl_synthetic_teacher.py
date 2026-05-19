#!/usr/bin/env python3
"""Generate a synthetic 'precomputed teacher' npz for NaCl + H2O.

Why this exists
---------------
For the NaCl/ionic substrate the project's MAD-SURF MACE model is
out-of-distribution (it was trained on coin metals + organic adsorbates,
not bulk ionic surfaces), and the user's real DFT teacher data (from Psi4
or VASP) is not yet on disk. To validate the framework end-to-end on the
ionic path we need *some* physically sensible (E, F) per pose. We compute
that here analytically: pure pairwise Coulomb between the H2O point
charges (TIP3P-like) and the Na/Cl substrate point charges, with PBC
images out to a configurable cutoff. This is dominantly the right
physics for ionic substrates.

This script also serves as a template — when real Psi4 output is
available, replace the ``coulomb_energy_and_forces`` function with a call
to a Psi4 wrapper that produces the same (energies, forces) for each
pose. The npz schema is fixed and consumed by
``pyBall/gridff_jax/teacher_backends/precomputed.py``.

NPZ output schema (consumed by PrecomputedTeacherBackend)
---------------------------------------------------------
  positions          (n_pose, n_atom, 3)  float64  Å
  pose_params        (n_pose, 7)          float64  [u, v, z, qw, qx, qy, qz]
  energies           (n_pose,)            float64  eV
  forces             (n_pose, n_atom, 3)  float64  eV/Å
  adsorbate_symbols  (n_atom,)            str

Usage
-----
  python make_nacl_synthetic_teacher.py \\
      --density-xyz cpp/common_resources/xyz/NaCl_1x1_L1_jax.xyz \\
      --out runs/nacl_psi4_dataset_synthetic.npz \\
      --n-u 4 --n-v 4 --n-z 8 --n-orient 1 --n-yaw 1
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.append(str(ROOT))

from pyBall.gridff_jax.config import DensityBackendConfig, GridConfig
from pyBall.gridff_jax.density_backends import make_density_backend
from pyBall.gridff_jax.pose_sampling.molecules import benchmark_adsorbates
from pyBall.gridff_jax.pose_sampling.sampler_6d import (
    Sampler6DConfig, build_6d_pose_batch,
)


# Coulomb constant in eV·Å / e²
COULOMB_K = 14.399645


def _h2o_with_tip3p_charges():
    """Override H2O charges with TIP3P (O=-0.834, H=+0.417)."""
    h2o = benchmark_adsorbates()["H2O"]
    h2o.charges = np.array([-0.834, +0.417, +0.417], dtype=float)
    return h2o


def coulomb_energy_and_forces(mol_positions, mol_charges,
                              sub_positions, sub_charges,
                              cell, image_range=2):
    """Pairwise Coulomb E and per-atom F with PBC images out to ±image_range.

    Parameters
    ----------
    mol_positions : (n_atom, 3)
    mol_charges   : (n_atom,)
    sub_positions : (n_sub, 3)
    sub_charges   : (n_sub,)
    cell          : (3, 3)   lattice vectors as rows
    image_range   : int      include images i,j ∈ [-image_range, +image_range]

    Returns (E_eV, F_per_atom (n_atom, 3) eV/Å).
    """
    a, b = cell[0], cell[1]
    image_shifts = []
    for i in range(-image_range, image_range + 1):
        for j in range(-image_range, image_range + 1):
            image_shifts.append(i * a + j * b)
    image_shifts = np.array(image_shifts)  # (n_image, 3)

    E = 0.0
    F = np.zeros_like(mol_positions)
    for shift in image_shifts:
        sub_pos_shifted = sub_positions + shift
        diff = mol_positions[:, None, :] - sub_pos_shifted[None, :, :]  # (n_mol, n_sub, 3)
        r = np.linalg.norm(diff, axis=-1)                                # (n_mol, n_sub)
        r_safe = np.maximum(r, 1.0e-6)
        qq = mol_charges[:, None] * sub_charges[None, :]                 # (n_mol, n_sub)
        E += float((COULOMB_K * qq / r_safe).sum())
        # F_i = sum_j qq * (r_i - r_j) / |r|^3 * COULOMB_K
        F += (COULOMB_K * qq[:, :, None] * diff / (r_safe[:, :, None] ** 3)).sum(axis=1)
    return E, F


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--density-xyz", required=True)
    p.add_argument("--out", required=True, help="Output npz path")
    p.add_argument("--n-u", type=int, default=4)
    p.add_argument("--n-v", type=int, default=4)
    p.add_argument("--n-z", type=int, default=8)
    p.add_argument("--z-min", type=float, default=2.0)
    p.add_argument("--z-max", type=float, default=5.5)
    p.add_argument("--n-orient", type=int, default=1)
    p.add_argument("--n-yaw", type=int, default=1)
    p.add_argument("--random-fraction", type=float, default=0.0)
    p.add_argument("--image-range", type=int, default=2,
                   help="Include PBC images i,j in [-N, N] for the Coulomb sum.")
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args()


def main():
    args = parse_args()

    print(f"[1/3] Loading substrate from {args.density_xyz} ...", flush=True)
    cfg = DensityBackendConfig(kind="surface_xyz", xyz_path=args.density_xyz)
    density = make_density_backend(cfg, grid=GridConfig()).load()
    print(f"      {density.natoms} atoms  cell={np.diag(density.cell).round(3)}",
          flush=True)

    adsorbate = _h2o_with_tip3p_charges()
    print(f"[2/3] Building 6D poses (n_u={args.n_u}, n_v={args.n_v}, "
          f"n_z={args.n_z}, n_orient={args.n_orient}, n_yaw={args.n_yaw}) ...",
          flush=True)
    sampler_cfg = Sampler6DConfig(
        n_u=args.n_u, n_v=args.n_v, n_z=args.n_z,
        z_range=(args.z_min, args.z_max), z_bias_power=1.5,
        n_orient=args.n_orient, n_yaw=args.n_yaw,
        tilt_max_deg=60.0,
        include_high_symmetry_sites=True,
        random_fraction=args.random_fraction, seed=args.seed,
    )
    poses = build_6d_pose_batch(density, adsorbate, sampler_cfg)
    n = len(poses.positions)
    print(f"      {n} poses generated", flush=True)

    print(f"[3/3] Computing analytical Coulomb (PBC images ±{args.image_range}) ...",
          flush=True)
    energies = np.zeros(n, dtype=float)
    forces = np.zeros((n, adsorbate.natoms, 3), dtype=float)
    sub_charges = np.asarray(density.charges, dtype=float)
    for i in range(n):
        mol_pos = np.asarray(poses.positions[i], dtype=float)
        # H2O charges are absolute; the rigid-body rotation does NOT permute them
        # because they are atom-indexed.
        e, f = coulomb_energy_and_forces(
            mol_pos, np.asarray(adsorbate.charges, dtype=float),
            density.positions, sub_charges,
            density.cell, image_range=args.image_range,
        )
        energies[i] = e
        forces[i] = f
    print(f"      E range = [{energies.min() * 1000:.1f}, {energies.max() * 1000:.1f}] meV",
          flush=True)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        out_path,
        positions=poses.positions,
        pose_params=poses.pose_params,
        energies=energies,
        forces=forces,
        adsorbate_symbols=np.asarray(adsorbate.symbols, dtype=object),
        metadata=np.array(
            {"source": "synthetic_coulomb", "image_range": args.image_range,
             "physics": "pairwise_coulomb_with_PBC_images_only",
             "warning": "no Pauli, no London, no dispersion — placeholder until real DFT"},
            dtype=object,
        ),
    )
    print(f"\n DONE. Wrote {out_path}", flush=True)
    print(f"   {n} poses, schema matches PrecomputedTeacherBackend", flush=True)


if __name__ == "__main__":
    main()
