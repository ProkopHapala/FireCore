#!/usr/bin/env python3
"""Generate documented sparse Hessian benchmark .npz files for external iterative-solver repos.

Run from repo root:
  python3 tests/tSiNCs/export_vibration_benchmarks.py
  python3 tests/tSiNCs/export_vibration_benchmarks.py --calibrate-only

Outputs (gitignored):
  tests/tSiNCs/fixtures/vibration_benchmarks/
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
TMMFF = REPO / "tests/tMMFF"
DEFAULT_OUT = REPO / "tests/tSiNCs/fixtures/vibration_benchmarks"

from pyBall.nanocrystal_gen import build_spherical_nanoparticle, save_xyz
from pyBall.vibration_benchmark import BenchmarkMeta, compute_mmff_hessian_bundle, export_benchmark_npz


def _ensure_tmmff_data():
    data = TMMFF / "data"
    if not data.exists():
        os.symlink(REPO / "cpp/common_resources", data)


def _bundle_from_xyz(name: str, xyz: Path, *, relax: bool, description: str, notes: list[str]) -> tuple[dict, BenchmarkMeta]:
    import numpy as np
    bundle = compute_mmff_hessian_bundle(str(xyz.resolve()), relax=relax, cwd=str(TMMFF))
    nnz = int(bundle["H_sparse_projected"].nnz)
    if bundle["n_negative_projected"] > 0:
        raise ValueError(f"{name}: projected Hessian has {bundle['n_negative_projected']} negative eigenvalues")
    if bundle["recon_vs_full_max_abs"] > 1e-2:
        notes.append(f"WARNING: shell blocks vs full H max|diff|={bundle['recon_vs_full_max_abs']:.3e}")
    meta = BenchmarkMeta(
        name=name,
        description=description,
        source_xyz=str(xyz.resolve()),
        natoms=int(bundle["natoms"]),
        ndof=int(bundle["ndof"]),
        dx=float(bundle["dx"]),
        n_shells=int(bundle["n_shells"]),
        max_neigh=int(bundle["neigh_idx"].shape[1]),
        rigid_shift=float(bundle["rigid_shift"]),
        relax_steps=int(bundle["relax_steps"]),
        nnz=nnz,
        n_negative_eigvals_raw=int(bundle["n_negative_raw"]),
        n_vibrational_modes=int(np.sum(bundle["omegas_modes_projected"] > 1e-4)),
        omega_min_vib=float(bundle["omega_min_vib"]),
        omega_max=float(bundle["omegas_modes_projected"].max()),
        eigh_seconds=float(bundle["eigh_seconds"]),
        notes=notes,
    )
    return bundle, meta


def calibrate_eigh_sizes(out_dir: Path):
    """Find diamond NC radius where dense eigh ~ 1s (target ndof for reference suite)."""
    _ensure_tmmff_data()
    work = out_dir / "calibration"
    work.mkdir(parents=True, exist_ok=True)
    print("=== EIGH size calibration (C diamond NC, relaxed, projected H) ===")
    rows = []
    for R in (3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0):
        nrep = max(4, int(R) + 2)
        elems, apos = build_spherical_nanoparticle(R=R, nrep=nrep, heavy_z=6)
        xyz = work / f"calib_R{R:.1f}.xyz"
        save_xyz(str(xyz), elems, apos, comment=f"calibration R={R}")
        try:
            b = compute_mmff_hessian_bundle(str(xyz), relax=True, cwd=str(TMMFF))
            rows.append((R, b["natoms"], b["ndof"], b["eigh_seconds"], b["n_negative_projected"]))
            print(f"  R={R:4.1f}  N={b['natoms']:4d}  DOF={b['ndof']:5d}  eigh={b['eigh_seconds']:.3f}s  neg_proj={b['n_negative_projected']}")
        except Exception as e:
            print(f"  R={R:4.1f}  FAILED: {e}")
    if not rows:
        return
    ok = [r for r in rows if r[3] < 1.5 and r[4] == 0]
    if ok:
        best = max(ok, key=lambda r: r[2])
        print(f"=== Recommended dense-reference max: DOF={best[2]} (R={best[0]}, N={best[1]}, eigh={best[3]:.2f}s) ===")
    else:
        print("=== No configuration met eigh<1.5s; use smallest systems only for dense reference ===")


def export_suite(out_dir: Path):
    _ensure_tmmff_data()
    struct = out_dir / "structures"
    struct.mkdir(parents=True, exist_ok=True)
    suite = []

    mol2 = REPO / "cpp/common_resources/mol/adamantane.mol2"
    suite.append((
        "adamantane",
        mol2,
        False,
        "C10H16 adamantane cage; MMFF minimum geometry (no relax). Dense eigh reference.",
        ["Internal vibrational modes only after rigid_shift=1e6."],
    ))

    for R, tag in ((4.0, "nc_C_R4"), (5.0, "nc_C_R5"), (6.0, "nc_C_R6"), (7.0, "nc_C_R7"), (8.0, "nc_C_R8")):
        nrep = max(4, int(R) + 2)
        elems, apos = build_spherical_nanoparticle(R=R, nrep=nrep, heavy_z=6, resolve_clashes=True)
        xyz = struct / f"{tag}_relaxed.xyz"
        save_xyz(str(xyz), elems, apos, comment=f"C diamond NC R={R} caps VSEPR+clash resolve")
        suite.append((
            tag,
            xyz,
            True,
            f"H-passivated carbon diamond nanocrystal, spherical R={R} Ang, MMFF-relaxed before Hessian.",
            ["H caps: pyBall.nanocrystal_gen VSEPR tetrahedral (not legacy -v bisector placement)."],
        ))

    fix_si = REPO / "tests/tSiNCs/fixtures/vibration_parallel/structures/si_G1_caps_only.xyz"
    if fix_si.is_file():
        suite.append((
            "si_G1_caps_only",
            fix_si,
            True,
            "Si nanocrystal from gen_nanocrystals.mjs (fixture G1), MMFF-relaxed.",
            ["From JS generator; skipped if MMFF topology build fails."],
        ))

    manifest_lines = ["# Vibration benchmark suite\n"]
    for name, xyz, relax, desc, notes in suite:
        print(f"\n=== Export {name} ===")
        try:
            bundle, meta = _bundle_from_xyz(name, Path(xyz), relax=relax, description=desc, notes=notes)
        except Exception as e:
            print(f"  SKIP {name}: {e}")
            manifest_lines.append(f"- `{name}.npz`: **SKIPPED** ({e})")
            continue
        out_npz = out_dir / f"{name}.npz"
        include_dense = bundle["ndof"] <= 3000
        export_benchmark_npz(str(out_npz), bundle, meta, include_dense=include_dense)
        manifest_lines.append(f"- `{out_npz.name}`: N={meta.natoms} DOF={meta.ndof} nnz={meta.nnz} eigh={meta.eigh_seconds:.3f}s dense={'yes' if include_dense else 'no'}")

    readme = out_dir / "README.md"
    readme.write_text(
        "\n".join(manifest_lines)
        + "\n\nRegenerate: `python3 tests/tSiNCs/export_vibration_benchmarks.py`\n"
        + "Format: see `doc/Topics/FTIR_Nanocrystals/Vibration_benchmark_npz.guide.md`\n"
        + "**Gitignored** — copy to external benchmark repo manually.\n"
    )
    print(f"\nWrote {readme}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--calibrate-only", action="store_true")
    ap.add_argument("--skip-calibrate", action="store_true", help="Skip eigh size sweep (use when suite only)")
    ap.add_argument("--outdir", type=str, default=str(DEFAULT_OUT))
    args = ap.parse_args()
    out_dir = Path(args.outdir)
    if args.calibrate_only:
        calibrate_eigh_sizes(out_dir)
    else:
        if not args.skip_calibrate:
            calibrate_eigh_sizes(out_dir)
        export_suite(out_dir)


if __name__ == "__main__":
    main()
