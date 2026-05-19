#!/usr/bin/env python3
"""Convert DFT (or any QM/MLIP) output → precomputed-teacher npz.

This is a TEMPLATE script. You will need to write 1 function — the one
that extracts (energy, forces) from your DFT/QM/MLIP output files — and
this driver will:

  1. Build the same 6D pose batch the GridFF-JAX driver builds,
  2. For each pose: generate the molecule positions in lab frame,
     write whatever input files your code needs, call your code, and
     parse the result via your extractor function,
  3. Save everything in the npz schema that
     PrecomputedTeacherBackend understands.

The result drops directly into the GridFF-JAX pipeline:

    python benchmark_substrate_6d.py \\
        --substrate-class ionic \\
        --density-xyz NaCl_1x1_L1_jax.xyz \\
        --adsorbate H2O \\
        --teacher-kind precomputed \\
        --teacher-npz path/to/your/output.npz

What you provide
----------------
Two things:

  (a) A function ``extract_energy_force(out_dir, mol_symbols, mol_positions,
                                        substrate_symbols, substrate_positions,
                                        cell) -> (E, F)``
      where E is a float (eV, interaction energy preferred) and F is an
      (n_atom, 3) ndarray (eV/Å, forces on the MOLECULE atoms only).
      If your method only gives energies, return (E, None) and the
      backend will fill F with zeros.

  (b) The CLI flags (substrate xyz, adsorbate, pose-grid parameters).
      Pose generation is identical to the benchmark driver's, so the
      pose_params will match exactly.

Three example extractors are provided below — pick the one that matches
your DFT/MLIP code (or write your own). They are placeholders, you have
to fill in the actual code-invocation logic.

Usage (with the synthetic Coulomb example shipped here):
  python dft_to_npz_template.py --substrate-xyz <surface_xyz> \\
      --adsorbate H2O --extractor coulomb_only --out path/to/out.npz \\
      --n-u 4 --n-v 4 --n-z 8

For real DFT (Psi4), set ``--extractor psi4_template`` and edit the
``run_psi4`` function below to match your group's standard input.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path
from typing import Callable

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.append(str(ROOT))

from pyBall.gridff_jax.config import DensityBackendConfig, GridConfig
from pyBall.gridff_jax.density_backends import make_density_backend
from pyBall.gridff_jax.pose_sampling.molecules import (
    benchmark_adsorbates, load_adsorbate_from_xyz,
)
from pyBall.gridff_jax.pose_sampling.sampler_6d import (
    Sampler6DConfig, build_6d_pose_batch,
)


# ===========================================================================
#  EXTRACTOR (a) — example: pure pairwise Coulomb (no DFT call).
#  Useful as a sanity check; matches what make_nacl_synthetic_teacher.py does.
# ===========================================================================
COULOMB_K = 14.399645   # eV·Å / e²


def coulomb_only(work_dir, mol_symbols, mol_positions, sub_symbols,
                 sub_positions, sub_charges, mol_charges, cell, image_range=2):
    a, b = cell[0], cell[1]
    shifts = [i * a + j * b
              for i in range(-image_range, image_range + 1)
              for j in range(-image_range, image_range + 1)]
    E = 0.0
    F = np.zeros_like(mol_positions)
    for shift in shifts:
        sub_pos = sub_positions + shift
        diff = mol_positions[:, None, :] - sub_pos[None, :, :]
        r = np.linalg.norm(diff, axis=-1)
        r_safe = np.maximum(r, 1.0e-6)
        qq = mol_charges[:, None] * sub_charges[None, :]
        E += float((COULOMB_K * qq / r_safe).sum())
        F += (COULOMB_K * qq[:, :, None] * diff / (r_safe[:, :, None] ** 3)).sum(axis=1)
    return E, F


# ===========================================================================
#  EXTRACTOR (b) — Psi4 template.  Edit run_psi4 to match your inputs.
# ===========================================================================
def psi4_template(work_dir, mol_symbols, mol_positions, sub_symbols,
                  sub_positions, sub_charges, mol_charges, cell,
                  method="b3lyp", basis="def2-svp", do_forces=False):
    """Example Psi4 wrapper. Writes an input file, calls psi4 via subprocess,
    parses the output for energy and (optionally) forces.

    Replace the placeholder commented-out subprocess call with your actual
    psi4 invocation. The energy returned should be the INTERACTION energy:
        E_interaction = E_complex − E_substrate − E_molecule

    The substrate is the 2-atom (or N-atom) point-charge xyz; do the slab
    + molecule full DFT, then subtract the slab-only and isolated-molecule
    references that you compute separately (and cache on first call)."""
    work_dir = Path(work_dir); work_dir.mkdir(parents=True, exist_ok=True)
    inp = work_dir / "input.dat"
    out = work_dir / "output.dat"

    # ---- write Psi4 input ----
    lines = [f"memory 4 gb", f"set basis {basis}", f"set scf_type df", ""]
    lines.append("molecule complex {")
    lines.append("0 1   # charge, spin multiplicity")
    for s, p in zip(sub_symbols, sub_positions):
        lines.append(f"  {s}  {p[0]:.6f}  {p[1]:.6f}  {p[2]:.6f}")
    for s, p in zip(mol_symbols, mol_positions):
        lines.append(f"  {s}  {p[0]:.6f}  {p[1]:.6f}  {p[2]:.6f}")
    lines.append("}")
    lines.append("")
    if do_forces:
        lines.append(f"E, wfn = gradient('{method}', return_wfn=True)")
        lines.append("forces = wfn.gradient().to_array()")
    else:
        lines.append(f"E = energy('{method}')")
    lines.append("print(f'E_TOTAL_eV {E * 27.211385:.10f}')")
    inp.write_text("\n".join(lines) + "\n")

    # ---- call Psi4 (PLACEHOLDER — uncomment & adjust for your setup) ----
    # subprocess.run(["psi4", str(inp), str(out)], check=True)

    # ---- parse output (PLACEHOLDER — adapt to your psi4 output format) ----
    # E_total = ...
    # F = ...   # forces on the molecule atoms, eV/Å
    # Subtract slab + molecule references (compute once and cache):
    # E_interaction = E_total − E_slab_ref − E_mol_ref
    # return E_interaction, F  (or (E_interaction, None) if no forces)

    raise NotImplementedError(
        "psi4_template is a stub. Fill in the subprocess call + output "
        "parser for your Psi4 setup, then return (E_interaction, F_or_None)."
    )


# ===========================================================================
#  EXTRACTOR (c) — VASP template.  Edit run_vasp to match your inputs.
# ===========================================================================
def vasp_template(work_dir, mol_symbols, mol_positions, sub_symbols,
                  sub_positions, sub_charges, mol_charges, cell):
    """Example VASP wrapper. Writes POSCAR + INCAR + KPOINTS + POTCAR-symlink,
    runs vasp_std, parses OUTCAR for energy and forces.

    Forces parsing is standard ASE; use ``ase.io.vasp.read_vasp_out`` or
    ``ase.calculators.vasp.Vasp.read_forces()``. Interaction energy =
    E_total − E_slab_ref − E_mol_ref (compute and cache the two references
    once at the top of the script)."""
    raise NotImplementedError(
        "vasp_template is a stub. Fill in POSCAR generation + vasp_std + "
        "OUTCAR parsing, then return (E_interaction, F).")


# ===========================================================================
#  Driver
# ===========================================================================
EXTRACTORS = {
    "coulomb_only": coulomb_only,
    "psi4_template": psi4_template,
    "vasp_template": vasp_template,
}


def _load_adsorbate(spec: str):
    builtins = benchmark_adsorbates()
    if spec in builtins:
        ads = builtins[spec]
    elif Path(spec).suffix.lower() in (".xyz", ".mol2"):
        ads = load_adsorbate_from_xyz(xyz_path=spec, name=Path(spec).stem,
                                      anchor_index=0, use_input_charges=False)
    else:
        raise ValueError(f"unknown adsorbate {spec!r}")
    # H2O TIP3P default charges if absent (common case)
    if ads.name == "H2O" and ads.charges is None:
        ads.charges = np.array([-0.834, +0.417, +0.417], dtype=float)
    return ads


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--substrate-xyz", required=True)
    p.add_argument("--adsorbate", required=True, help="Built-in name or XYZ path")
    p.add_argument("--out", required=True, help="Output npz path")
    p.add_argument("--extractor", required=True, choices=list(EXTRACTORS.keys()),
                   help="Which (E, F) extractor to use per pose")
    p.add_argument("--work-dir", default="/tmp/dft_to_npz_work",
                   help="Per-pose subdirs go here (DFT scratch)")
    p.add_argument("--n-u", type=int, default=4)
    p.add_argument("--n-v", type=int, default=4)
    p.add_argument("--n-z", type=int, default=8)
    p.add_argument("--z-min", type=float, default=2.0)
    p.add_argument("--z-max", type=float, default=5.5)
    p.add_argument("--n-orient", type=int, default=1)
    p.add_argument("--n-yaw", type=int, default=1)
    p.add_argument("--random-fraction", type=float, default=0.0)
    p.add_argument("--energies-only", action="store_true",
                   help="Skip force computation/return; npz will have no 'forces' array.")
    p.add_argument("--max-poses", type=int, default=None,
                   help="Cap the number of poses (useful for smoke tests).")
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args()


def main():
    args = parse_args()
    extractor: Callable = EXTRACTORS[args.extractor]

    print(f"[1/3] Loading substrate from {args.substrate_xyz} ...", flush=True)
    dcfg = DensityBackendConfig(kind="surface_xyz", xyz_path=args.substrate_xyz)
    density = make_density_backend(dcfg, grid=GridConfig()).load()
    print(f"      {density.natoms} atoms, cell={np.diag(density.cell).round(3)}",
          flush=True)

    adsorbate = _load_adsorbate(args.adsorbate)
    print(f"[2/3] Building 6D poses for {adsorbate.name} "
          f"(n_u={args.n_u}, n_v={args.n_v}, n_z={args.n_z}, "
          f"n_orient={args.n_orient}, n_yaw={args.n_yaw}) ...", flush=True)
    sampler = Sampler6DConfig(
        n_u=args.n_u, n_v=args.n_v, n_z=args.n_z,
        z_range=(args.z_min, args.z_max), z_bias_power=1.5,
        n_orient=args.n_orient, n_yaw=args.n_yaw, tilt_max_deg=60.0,
        include_high_symmetry_sites=True,
        random_fraction=args.random_fraction, seed=args.seed,
    )
    poses = build_6d_pose_batch(density, adsorbate, sampler)
    n = len(poses.positions)
    if args.max_poses and args.max_poses < n:
        n = args.max_poses
        print(f"      capping to {n} poses (--max-poses)", flush=True)

    print(f"[3/3] Running extractor '{args.extractor}' on {n} poses ...", flush=True)
    energies = np.zeros(n, dtype=float)
    forces = np.zeros((n, adsorbate.natoms, 3), dtype=float)
    forces_available = not args.energies_only

    work_root = Path(args.work_dir)
    work_root.mkdir(parents=True, exist_ok=True)

    for i in range(n):
        mol_pos = np.asarray(poses.positions[i], dtype=float)
        sub_charges = np.asarray(density.charges, dtype=float) if density.charges is not None \
                      else np.zeros(density.natoms)
        mol_charges = np.asarray(adsorbate.charges, dtype=float) if adsorbate.charges is not None \
                      else np.zeros(adsorbate.natoms)
        E, F = extractor(
            work_root / f"pose_{i:05d}",
            adsorbate.symbols, mol_pos,
            density.symbols, density.positions, sub_charges, mol_charges,
            density.cell,
        )
        energies[i] = E
        if F is None or args.energies_only:
            forces_available = False
            forces[i] = 0.0
        else:
            forces[i] = F
        if i % max(1, n // 20) == 0:
            print(f"      pose {i}/{n}  E={E * 1000:+.1f} meV", flush=True)

    # ---- Write npz ----
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    payload = dict(
        pose_params=poses.pose_params[:n],
        energies=energies,
        positions=poses.positions[:n],
        adsorbate_symbols=np.asarray(adsorbate.symbols, dtype=object),
        metadata=np.array(
            {"extractor": args.extractor,
             "substrate_xyz": args.substrate_xyz,
             "adsorbate": adsorbate.name,
             "n_poses": int(n)},
            dtype=object,
        ),
    )
    if forces_available:
        payload["forces"] = forces
        print(f"      Wrote {n} poses with energies AND forces.", flush=True)
    else:
        print(f"      Wrote {n} poses with energies ONLY (driver will auto-switch "
              f"--fit-mode energy).", flush=True)
    np.savez(out_path, **payload)
    print(f"\n DONE. Wrote {out_path}", flush=True)
    print(f"   Feed this to:\n     "
          f"python tests/tGridFFJax/benchmark_substrate_6d.py \\\n"
          f"       --substrate-class ionic --density-xyz {args.substrate_xyz} \\\n"
          f"       --adsorbate {adsorbate.name} \\\n"
          f"       --teacher-kind precomputed --teacher-npz {out_path} \\\n"
          f"       --out-dir tests/tGridFFJax/runs/<your_run_name>", flush=True)


if __name__ == "__main__":
    main()
