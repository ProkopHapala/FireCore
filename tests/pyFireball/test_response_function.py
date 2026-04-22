#!/usr/bin/env python3
"""STM Response Function Parity Test

Compares three methods for computing the STM response amplitude:
  1. Full matrix: explicit inversion of augmented (tip+sample) A_full = (E+iη)S-H
  2. Block elimination (G0): precomputed G0 = A_ss^{-1}, per-point O(n_s^2)
  3. LU factorization: precomputed LU of A_ss, per-point forward/back substitution

Also compares against OpenCL GPU kernel (response_amplitude_exp).
Outputs 2x3 panel plots per MO: ψ(r), LDOS, and all three response methods.

Run from:
  cd tests/pyFireball
  python test_response_function.py
  python test_response_function.py --xyz ../../cpp/common_resources/xyz/PTCDA.xyz --mo 74
"""

import os
import sys
import time
import argparse
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.path.dirname(__file__), "..", ".."))

from pyBall.AtomicSystem import AtomicSystem
from pyBall import FireCore as fc
from pyBall.FireballOCL.STM_utils import (
    set_export_dir, save_plot, plot_atoms,
    get_orbital_layout, build_plane_grid,
    project_orbital_to_points
)
from pyBall.FireballOCL.STM import (
    build_inter_system_blocks_exp_sk, _atom_orb_starts,
    response_amplitude_map, response_amplitude_full_matrix, response_amplitude_lu
)

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.normpath(os.path.join(_THIS_DIR, "..", ".."))


def _ensure_fdata():
    src = os.path.join(_REPO_ROOT, "tests", "Fireball", "Fdata_HCNOS")
    dst = os.path.join(_THIS_DIR, "Fdata")
    if os.path.realpath(dst) if os.path.exists(dst) else "" != os.path.realpath(src):
        if os.path.lexists(dst): os.unlink(dst)
        os.symlink(src, dst)


def _homo_lumo(eigen):
    occ = np.where(np.asarray(eigen) < 0.0)[0]
    homo = int(occ[-1]) if len(occ) > 0 else len(eigen)//2 - 1
    return homo, homo + 1


def main():
    parser = argparse.ArgumentParser(description="STM Response Function Parity Test")
    parser.add_argument("--xyz", default="../../cpp/common_resources/xyz/PTCDA.xyz")
    parser.add_argument("--mo", type=int, default=None, help="MO index (1-based). Default: HOMO")
    parser.add_argument("--z_orbital", type=float, default=1.0, help="Orbital projection height (Å)")
    parser.add_argument("--z_tip", type=float, default=3.0, help="Tip height for STM response (Å)")
    parser.add_argument("--size", type=float, default=20.0, help="Grid extent (Å)")
    parser.add_argument("--n", type=int, default=40, help="Grid resolution (small for CPU speed)")
    parser.add_argument("--rcut", type=float, default=8.0, help="Tip-sample coupling cutoff (Å)")
    parser.add_argument("--beta", type=float, default=1.0, help="Exponential decay β (Å⁻¹)")
    parser.add_argument("--r0", type=float, default=3.0, help="Reference distance r0 (Å)")
    parser.add_argument("--eta", type=float, default=1e-2, help="Broadening η (eV)")
    parser.add_argument("--E_tip", type=float, default=0.0, help="Tip onsite energy (eV)")
    parser.add_argument("--outdir", default="export/response_function")
    parser.add_argument("--skip_full", action="store_true", help="Skip full matrix method (slow for large molecules)")
    parser.add_argument("--fortran_mode", type=int, default=None, help="Fortran tipResponse2points mode: 0=expSK, 1=doscentros. If set, adds Fortran response panel.")
    parser.add_argument("--tipZ", type=int, default=2, help="Tip species index for Fortran (for H2O: H=1, O=2)")
    parser.add_argument("--fortran_simple", type=int, default=None, help="Fortran tipResponseSimple2points mode: 0=expSK, 1=doscentros. If set, adds simplified Fortran response panel.")
    parser.add_argument("--tip_src", type=str, default="s", help="Tip source coefficients. 's' activates only s. 'all' activates all orbitals equally. Or comma list of floats length ntip.")
    parser.add_argument("--python_simple", action="store_true", help="Add Python SIMPLE response panel (rhs=A_ts^H tip_src, solve A_ss y=rhs, resp=|c^H y|^2) for debugging vs Fortran SIMPLE.")
    args = parser.parse_args()

    os.chdir(_THIS_DIR)
    _ensure_fdata()

    xyz_path = args.xyz if os.path.isabs(args.xyz) else os.path.join(_THIS_DIR, args.xyz)
    mol = AtomicSystem(fname=xyz_path)
    mol_name = os.path.splitext(os.path.basename(xyz_path))[0]
    export_dir = set_export_dir(os.path.join(_THIS_DIR, args.outdir, mol_name))

    fc.setVerbosity(0)
    fc.initialize(atomType=mol.atypes, atomPos=mol.apos, verbosity=0)
    dims = fc.get_HS_dims()
    norb = int(dims.norbitals)
    natoms = int(dims.natoms)

    print(f"Running SCF for {mol_name} (natoms={natoms}, norb={norb})...")
    fc.evalForce(mol.apos, nmax_scf=200)
    eigen = fc.get_eigen(ikp=1, norb=norb)

    homo, lumo = _homo_lumo(eigen)
    mo0 = (args.mo - 1) if args.mo is not None else homo
    E = float(eigen[mo0])
    label = "HOMO" if mo0 == homo else ("LUMO" if mo0 == lumo else f"MO{mo0+1}")
    print(f"MO {mo0+1} ({label}) E={E:+.4f} eV")

    # Get H, S, eigenvectors
    C_fc = fc.get_wfcoef(norb=norb)

    # Dense H, S from Fortran ktransform (Gamma point) — same as test_mo_vs_ldos.py
    kpt = np.array((0.0, 0.0, 0.0), dtype=np.float64)
    Hk_fc, Sk_fc = fc.get_HS_k(kpt, norb)
    H_s = np.asarray(Hk_fc, dtype=np.complex128)
    S_s = np.asarray(Sk_fc, dtype=np.complex128)

    # neighbor metadata for orbital layout (coupling builder only)
    data = fc.get_HS_neighs(dims)
    data = fc.get_HS_sparse(dims, data)
    norb_per, _ = get_orbital_layout(data, natoms)
    print(f"  H_s shape: {H_s.shape}  (norb={norb})")

    # Build orbital projection grid (z_orbital for short-range orbitals)
    X_orb, Y_orb, pts_orb, extent_orb = build_plane_grid(mol.apos, z=args.z_orbital, size=args.size, n=args.n)
    n = args.n
    pts_orb_f64 = pts_orb.astype(np.float64)
    print(f"  Orbital grid: {n}x{n}={len(pts_orb)} points, z={args.z_orbital} Å, size={args.size} Å")

    # Build tip scan grid (z_tip for longer-range STM response)
    X_tip, Y_tip, pts_tip, extent_tip = build_plane_grid(mol.apos, z=args.z_tip, size=args.size, n=args.n)
    pts_tip_f64 = pts_tip.astype(np.float64)
    print(f"  Tip scan grid: {n}x{n}={len(pts_tip)} points, z={args.z_tip} Å, size={args.size} Å")

    # Tip metadata (single H atom, 1s orbital)
    tip_norb_per = np.array([1], dtype=np.int32)
    smp_norb_per = norb_per.astype(np.int32)

    def _parse_tip_src(ntip):
        s = str(args.tip_src).strip().lower()
        if s == 's':
            v = np.zeros(ntip, dtype=np.float64); v[0] = 1.0; return v
        if s == 'all':
            v = np.ones(ntip, dtype=np.float64); return v
        parts = [p.strip() for p in str(args.tip_src).split(',') if p.strip()!='']
        v = np.array([float(p) for p in parts], dtype=np.float64)
        if v.size != ntip:
            raise ValueError(f"--tip_src expects ntip={ntip} values, got {v.size}")
        return v

    # Coupling builder (exponential SK)
    def coupling_exp(tip_pos):
        return build_inter_system_blocks_exp_sk(
            np.array([1], dtype=np.int32), tip_pos, tip_norb_per,
            mol.atypes, mol.apos, smp_norb_per,
            rcut=args.rcut, beta=args.beta, r0=args.r0,
            A_ss=-1.0, A_sp=-1.0, A_pp_sig=-1.0, A_pp_pi=1.0, overlap_scale=0.0
        )

    # ── (1) Fortran ψ(r) reference
    print("  Computing Fortran ψ(r)...")
    psi_f = np.array(fc.orb2points(pts_orb_f64, iMO=mo0+1, ikpoint=1), dtype=np.float64).reshape(n, n).T

    # ── (2) OpenCL ψ(r) true basis
    print("  Computing OpenCL ψ(r)...")
    fdata_basis = os.path.join(_THIS_DIR, "Fdata", "basis")
    orb2atom = np.array([ia for ia in range(natoms) for _ in range(int(norb_per[ia]))], dtype=np.int32)
    psi_cl = np.array(project_orbital_to_points(
        C_fc, mo0, mol.atypes, mol.apos, orb2atom, smp_norb_per, fdata_basis, points=pts_orb
    ), dtype=np.float64).reshape(n, n).T

    # ── (3) Method 1: Block elimination with G0 (existing reference)
    print("  Computing response: block elimination (G0)...")
    t0 = time.time()
    resp_g0 = response_amplitude_map(
        pts_tip_f64, H_s, S_s, C_fc, mo0, E,
        np.array([1], dtype=np.int32), np.zeros((1,3)), tip_norb_per,
        mol.atypes, mol.apos, smp_norb_per,
        coupling_exp, eta=args.eta, E_tip=args.E_tip
    ).reshape(n, n).T
    t_g0 = time.time() - t0
    print(f"    done in {t_g0:.1f}s")

    # ── (4) Method 2: LU factorization
    print("  Computing response: LU factorization...")
    t0 = time.time()
    resp_lu = response_amplitude_lu(
        pts_tip_f64, H_s, S_s, C_fc, mo0, E,
        tip_norb_per, smp_norb_per,
        coupling_exp, eta=args.eta, E_tip=args.E_tip
    ).reshape(n, n).T
    t_lu = time.time() - t0
    print(f"    done in {t_lu:.1f}s")

    d_lu = float(np.max(np.abs(resp_lu - resp_g0)))
    print(f"  [PARITY] max|resp_LU - resp_G0| = {d_lu:.3e}")

    # ── (4b) Python SIMPLE (matches Fortran SIMPLE algebra)
    resp_py_simple = None
    if args.python_simple:
        print("  Computing response: Python SIMPLE (LU)...")
        t0 = time.time()
        ntip = 1
        tip_src = np.zeros(ntip, dtype=np.complex128); tip_src[0] = 1.0
        tip_norb_per1 = np.array([1], dtype=np.int32)
        from pyBall.FireballOCL.STM import response_amplitude_simple_lu
        resp_py_simple = response_amplitude_simple_lu(
            pts_tip_f64, H_s, S_s, C_fc, mo0, E,
            tip_norb_per1, smp_norb_per,
            coupling_exp, tip_src=tip_src, eta=args.eta
        ).reshape(n, n).T
        t_py_simple = time.time() - t0
        L = resp_py_simple[:, :n//2]; R = resp_py_simple[:, n//2:][:, ::-1]
        d_sym_py_simple = float(np.max(np.abs(L - R)))
        print(f"    done in {t_py_simple:.1f}s")
        print(f"  [SYMMETRY PY SIMPLE] max|L-R| = {d_sym_py_simple:.3e}  (mean={np.mean(np.abs(L-R)):.3e})")

    # ── (5) Method 3: Full matrix inversion (skip for large systems)
    resp_full = None
    if not args.skip_full and norb <= 200:
        print(f"  Computing response: full matrix (norb={norb}, {len(pts_tip)} points)...")
        t0 = time.time()
        resp_full = response_amplitude_full_matrix(
            pts_tip_f64, H_s, S_s, C_fc, mo0, E,
            tip_norb_per, smp_norb_per,
            coupling_exp, eta=args.eta, E_tip=args.E_tip
        ).reshape(n, n).T
        t_full = time.time() - t0
        print(f"    done in {t_full:.1f}s")
        d_full = float(np.max(np.abs(resp_full - resp_g0)))
        print(f"  [PARITY] max|resp_Full - resp_G0| = {d_full:.3e}")
    elif not args.skip_full:
        print(f"  Skipping full matrix (norb={norb} > 200). Use --skip_full to suppress this check.")

    # ── Fortran tipResponse2points (if requested)
    resp_fortran = None
    if args.fortran_mode is not None:
        mode_name = "expSK" if args.fortran_mode == 0 else "doscentros"
        print(f"  Computing response: Fortran mode {args.fortran_mode} ({mode_name})...")
        t0 = time.time()
        resp_fortran = fc.tipResponse2points(
            pts_tip_f64,
            E=E, eta=args.eta, mode=args.fortran_mode, iMO=mo0+1, ikpoint=1, tipZ=args.tipZ,
            rcut=args.rcut, beta=args.beta, r0=args.r0,
            A_ss=-1.0, A_sp=-1.0, A_pp_sig=-1.0, A_pp_pi=+1.0, A_ps=-1.0,
            overlap_scale=0.0, E_tip=args.E_tip
        ).reshape(n, n).T
        t_fortran = time.time() - t0
        print(f"    done in {t_fortran:.1f}s")
        # Symmetry check
        L = resp_fortran[:, :n//2]; R = resp_fortran[:, n//2:][:, ::-1]
        d_sym = float(np.max(np.abs(L - R)))
        print(f"  [SYMMETRY] max|L-R| = {d_sym:.3e}  (mean={np.mean(np.abs(L-R)):.3e})")

    # ── Fortran simplified response (if requested)
    resp_fortran_simple = None
    if args.fortran_simple is not None:
        mode_name = "expSK" if args.fortran_simple == 0 else "doscentros"
        ntip = int(fc.get_num_orb(args.tipZ))
        tip_src = _parse_tip_src(ntip)
        print(f"  Computing response: Fortran SIMPLE mode {args.fortran_simple} ({mode_name}), ntip={len(tip_src)}...")
        t0 = time.time()
        resp_fortran_simple = fc.tipResponseSimple2points(
            pts_tip_f64, tip_src,
            E=E, eta=args.eta, mode=args.fortran_simple, iMO=mo0+1, ikpoint=1, tipZ=args.tipZ,
            rcut=args.rcut, beta=args.beta, r0=args.r0,
            A_ss=-1.0, A_sp=-1.0, A_pp_sig=-1.0, A_pp_pi=+1.0, A_ps=-1.0,
            overlap_scale=0.0
        ).reshape(n, n).T
        t_fortran_simple = time.time() - t0
        print(f"    done in {t_fortran_simple:.1f}s")
        L = resp_fortran_simple[:, :n//2]; R = resp_fortran_simple[:, n//2:][:, ::-1]
        d_sym_simple = float(np.max(np.abs(L - R)))
        print(f"  [SYMMETRY SIMPLE] max|L-R| = {d_sym_simple:.3e}  (mean={np.mean(np.abs(L-R)):.3e})")

    # ── Plot
    def _panel_mo(ax, Z, title):
        vmax = max(abs(np.min(Z)), abs(np.max(Z)))
        if vmax < 1e-30: vmax = 1.0
        im = ax.imshow(Z, origin='lower', extent=extent_orb, cmap='bwr', vmin=-vmax, vmax=vmax, aspect='equal')
        ax.set_title(title); ax.set_xlabel('x (Å)'); ax.set_ylabel('y (Å)')
        plot_atoms(ax, mol.apos, mol.atypes, color='green', ms=4)
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    def _panel_pos(ax, Z, title, cmap='hot'):
        vmax = np.max(Z)
        if vmax < 1e-30: vmax = 1.0
        im = ax.imshow(Z, origin='lower', extent=extent_tip, cmap=cmap, vmin=0, vmax=vmax, aspect='equal')
        ax.set_title(title); ax.set_xlabel('x (Å)'); ax.set_ylabel('y (Å)')
        plot_atoms(ax, mol.apos, mol.atypes, color='cyan', ms=4)
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    # Determine rows: 2 base rows + optional full matrix + optional Fortran
    base_rows = 2
    extra_rows = 0
    if resp_full is not None: extra_rows += 1
    if resp_py_simple is not None: extra_rows += 1
    if resp_fortran is not None: extra_rows += 1
    if resp_fortran_simple is not None: extra_rows += 1
    nrows = base_rows + extra_rows

    fig, axes = plt.subplots(nrows, 3, figsize=(18, 6*nrows))
    fig.suptitle(f"{mol_name} MO {mo0+1} ({label}) E={E:+.4f} eV  z_orb={args.z_orbital:.1f} Å  z_tip={args.z_tip:.1f} Å  η={args.eta:.0e}")

    _panel_mo(axes[0, 0], psi_f,  "ψ(r) Fortran orb2points")
    _panel_mo(axes[0, 1], psi_cl, "ψ(r) OpenCL true basis")
    _panel_pos(axes[0, 2], resp_g0, f"Response G0 ({t_g0:.1f}s)", cmap='hot')

    _panel_pos(axes[1, 0], resp_lu, f"Response LU ({t_lu:.1f}s) Δ={d_lu:.1e}", cmap='hot')
    _panel_pos(axes[1, 1], np.abs(resp_lu - resp_g0), "Diff |LU - G0|", cmap='viridis')

    row_idx = 2
    if resp_full is not None:
        _panel_pos(axes[row_idx, 0], resp_full, f"Response Full ({t_full:.1f}s)", cmap='hot')
        _panel_pos(axes[row_idx, 1], np.abs(resp_full - resp_g0), "Diff |Full - G0|", cmap='viridis')
        _panel_pos(axes[row_idx, 2], np.abs(resp_full - resp_lu), "Diff |Full - LU|", cmap='viridis')
        row_idx += 1

    if resp_py_simple is not None:
        _panel_pos(axes[row_idx, 0], resp_py_simple, f"Python SIMPLE (LU, {t_py_simple:.1f}s)\nmax|L-R|={d_sym_py_simple:.1e}", cmap='hot')
        _panel_pos(axes[row_idx, 1], np.abs(resp_py_simple - resp_g0), "Diff |PySimple - G0|", cmap='viridis')
        axes[row_idx, 2].axis('off')
        axes[row_idx, 2].text(0.5, 0.5,
            f"max|PySimple-G0| = {np.max(np.abs(resp_py_simple - resp_g0)):.3e}\nmax|L-R| = {d_sym_py_simple:.3e}",
            ha='center', va='center', fontsize=12, transform=axes[row_idx, 2].transAxes)
        row_idx += 1
    elif resp_fortran is None:
        # If no full matrix and no Fortran, use placeholder
        axes[1, 2].axis('off')
        axes[1, 2].text(0.5, 0.5,
            f"max|LU-G0| = {d_lu:.3e}\n\nFull matrix skipped\n(norb={norb} > 200)",
            ha='center', va='center', fontsize=12, transform=axes[1, 2].transAxes)

    if resp_fortran is not None:
        mode_name = "expSK" if args.fortran_mode == 0 else "doscentros"
        _panel_pos(axes[row_idx, 0], resp_fortran, f"Fortran mode {args.fortran_mode} ({mode_name}, {t_fortran:.1f}s)\nmax|L-R|={d_sym:.1e}", cmap='hot')
        _panel_pos(axes[row_idx, 1], np.abs(resp_fortran - resp_g0), "Diff |Fortran - G0|", cmap='viridis')
        axes[row_idx, 2].axis('off')
        axes[row_idx, 2].text(0.5, 0.5,
            f"max|Fortran-G0| = {np.max(np.abs(resp_fortran - resp_g0)):.3e}\nmax|L-R| = {d_sym:.3e}",
            ha='center', va='center', fontsize=12, transform=axes[row_idx, 2].transAxes)
        row_idx += 1

    if resp_fortran_simple is not None:
        mode_name = "expSK" if args.fortran_simple == 0 else "doscentros"
        _panel_pos(axes[row_idx, 0], resp_fortran_simple, f"Fortran SIMPLE mode {args.fortran_simple} ({mode_name}, {t_fortran_simple:.1f}s)\nmax|L-R|={d_sym_simple:.1e}", cmap='hot')
        _panel_pos(axes[row_idx, 1], np.abs(resp_fortran_simple - resp_g0), "Diff |FortranSimple - G0|", cmap='viridis')
        axes[row_idx, 2].axis('off')
        axes[row_idx, 2].text(0.5, 0.5,
            f"max|FortranSimple-G0| = {np.max(np.abs(resp_fortran_simple - resp_g0)):.3e}\nmax|L-R| = {d_sym_simple:.3e}",
            ha='center', va='center', fontsize=12, transform=axes[row_idx, 2].transAxes)

    plt.tight_layout()
    save_plot(fig, f"mo{mo0+1:04d}_{label}_response_parity", export_dir)
    plt.close(fig)

    print(f"\nDone. Output in {export_dir}")


if __name__ == "__main__":
    main()
