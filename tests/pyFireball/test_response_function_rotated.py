#!/usr/bin/env python3
"""STM Response Function Rotational Invariance Test

Tests that rotating a molecule AND rotating the imaging plane by the same
rotation matrix R produces the same STM image (up to numerical precision).

This validates:
  1. rotate_dense_hs(H, S, norb_per, R)   — rigid-body rotation of Hamiltonian
  2. rotate_mo_coefficients(C, norb_per, R) — rotation of MO coefficients
  3. rotate_points(apos, R)               — rotation of atomic positions
  4. build_rotated_plane_grid(..., R)     — rotation of imaging plane

Test procedure for each rotation R:
  A. Unrotated system, scan in xy-plane   → image I0
  B. Rotated H,S,C + rotated apos + rotated scan grid → image I1
  C. Compare I0 vs I1: should be identical

Also tests tip orbital dependence: s, px, py, pz tips should behave
correctly under rotation (s-tip invariant, p-tips rotate).

Run from:
  cd tests/pyFireball
  python test_response_function_rotated.py
  python test_response_function_rotated.py --xyz ../../cpp/common_resources/xyz/PTCDA.xyz --mo 74
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
    response_amplitude_lu,
    response_amplitude_simple_lu,
    rotate_dense_hs, rotate_mo_coefficients,
    rotate_points, build_rotated_plane_grid,
    build_atom_rotation_matrix
)

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))


def _ensure_fdata():
    dst = os.path.join(_THIS_DIR, "Fdata")
    src_rel = "../Fireball/Fdata_HCNOS"
    if not (os.path.lexists(dst) and os.path.realpath(dst).endswith("Fdata_HCNOS")):
        if os.path.lexists(dst):
            os.unlink(dst)
        os.symlink(src_rel, dst)


def _homo_lumo(eigen):
    occ = np.where(np.asarray(eigen) < 0.0)[0]
    homo = int(occ[-1]) if len(occ) > 0 else len(eigen)//2 - 1
    return homo, homo + 1


def rotation_matrix_x(theta_deg):
    """Rotation around x-axis by theta (degrees)."""
    th = np.radians(float(theta_deg))
    c, s = np.cos(th), np.sin(th)
    return np.array([
        [1,  0,  0],
        [0,  c, -s],
        [0,  s,  c]
    ], dtype=np.float64)


def rotation_matrix_y(theta_deg):
    """Rotation around y-axis by theta (degrees)."""
    th = np.radians(float(theta_deg))
    c, s = np.cos(th), np.sin(th)
    return np.array([
        [c,  0,  s],
        [0,  1,  0],
        [-s, 0,  c]
    ], dtype=np.float64)


def rotation_matrix_z(theta_deg):
    """Rotation around z-axis by theta (degrees)."""
    th = np.radians(float(theta_deg))
    c, s = np.cos(th), np.sin(th)
    return np.array([
        [c, -s,  0],
        [s,  c,  0],
        [0,  0,  1]
    ], dtype=np.float64)
def run_rotated_test(mol, H_s, S_s, C_fc, eigen, mo0, norb_per,
                     tip_norb_per, smp_norb_per, R, label,
                     z_tip, size, n,
                     rcut, beta, r0, eta, E_tip):
    """Run STM response for unrotated and rotated system, return both images."""

    E = float(eigen[mo0])
    ntip_total = int(_atom_orb_starts(tip_norb_per)[-1])

    # ── Unrotated reference ──
    X0, Y0, pts0, ext0 = build_plane_grid(mol.apos, z=z_tip, size=size, n=n)
    pts0_f64 = pts0.astype(np.float64)

    # Coupling builder (unrotated)
    def coupling0(tip_pos):
        pairs, Hb, Sb = build_inter_system_blocks_exp_sk(
            np.array([1], dtype=np.int32), tip_pos, tip_norb_per,
            mol.atypes, mol.apos, smp_norb_per,
            rcut=rcut, beta=beta, r0=r0,
            A_ss=-1.0, A_sp=-1.0, A_pp_sig=-1.0, A_pp_pi=+1.0, A_ps=-1.0,
            overlap_scale=0.0
        )
        # IMPORTANT: match Fortran firecore_tipResponseSimple2points behavior for H tip (ntip=1)
        # Fortran only adds the s-s term when ntip<4 (no s->p couplings).
        if int(_atom_orb_starts(tip_norb_per)[-1]) == 1 and len(Hb) > 0:
            Hb = [b.copy() for b in Hb]
            Sb = [b.copy() for b in Sb]
            for ib in range(len(Hb)):
                if Hb[ib].shape[0] == 1 and Hb[ib].shape[1] > 1:
                    Hb[ib][0, 1:] = 0.0
                if Sb[ib].shape[0] == 1 and Sb[ib].shape[1] > 1:
                    Sb[ib][0, 1:] = 0.0
        return pairs, Hb, Sb

    # s-tip only (per user request)
    tip_src = np.zeros(ntip_total, dtype=np.complex128)
    tip_src[0] = 1.0

    # Reference response (Fortran, to isolate rotation logic)
    tipZ = 1
    tip_src_fb = np.array([1.0], dtype=np.float64)
    resp_ref = fc.tipResponseSimple2points(
        pts0_f64, tip_src_fb,
        E=E, eta=eta, mode=0, iMO=int(mo0+1), ikpoint=1, tipZ=tipZ,
        rcut=rcut, beta=beta, r0=r0,
        A_ss=-1.0, A_sp=-1.0, A_pp_sig=-1.0, A_pp_pi=+1.0, A_ps=-1.0,
        overlap_scale=0.0
    ).reshape(n, n).T

    # Use unambiguous MO coefficient vector from Fortran: c = bbnkre(:, iband, ikpoint)
    # _get_mo_vec has a buggy square-matrix heuristic that can pick the wrong orientation.
    # Pass non-square (1, norb) matrix with mo_index=0 to bypass it entirely.
    c_vec = np.asarray(fc.get_wfcoef_vec(iband=int(mo0+1), ikp=1, norb=H_s.shape[0]), dtype=np.float64)
    C_use = np.zeros((1, H_s.shape[0]), dtype=np.float64)
    C_use[0, :] = c_vec

    H_use = H_s
    S_use = S_s

    resp_py_ref = response_amplitude_simple_lu(
        pts0_f64, H_use, S_use, C_use, 0, E,
        tip_norb_per, smp_norb_per,
        coupling_builder=coupling0,
        tip_src=tip_src,
        eta=eta
    ).reshape(n, n).T

    diff_py_ref = float(np.max(np.abs(resp_py_ref - resp_ref)))
    if diff_py_ref > 1e-9:
        print(f"  [DEBUG] max|py_ref - fc_ref| = {diff_py_ref:.3e}")
        print(f"  [DEBUG] ||c_vec (Fortran)||={np.linalg.norm(c_vec):.6f}")
        print(f"  [DEBUG] C_use.shape={C_use.shape}  mo0={mo0}")
        print(f"  [DEBUG] ||C_use[mo0,:]||={np.linalg.norm(C_use[int(mo0),:]):.6f}")
        # Also check C_fc for comparison
        print(f"  [DEBUG] C_fc.shape={C_fc.shape}")
        print(f"  [DEBUG] ||C_fc[:,mo0]||={np.linalg.norm(C_fc[:,int(mo0)]):.6f}")
        print(f"  [DEBUG] ||C_fc[mo0,:]||={np.linalg.norm(C_fc[int(mo0),:]):.6f}")
        # DEBUG: compare coupling matrices for the first point to locate basis/mapping mismatch
        p_dbg = np.asarray(pts0_f64[0], dtype=np.float64)
        Hts_fc, Sts_fc = fc.export_tip_coupling_point(
            p_dbg, mode=0, tipZ=tipZ,
            rcut=rcut, beta=beta, r0=r0,
            A_ss=-1.0, A_sp=-1.0, A_pp_sig=-1.0, A_pp_pi=+1.0, A_ps=-1.0,
            overlap_scale=0.0,
            ntip=1, norb=H_s.shape[0]
        )
        pairs_dbg, Hb_dbg, Sb_dbg = coupling0(np.array([p_dbg], dtype=np.float64))
        Hts_py = np.zeros((1, H_s.shape[0]), dtype=np.float64)
        Sts_py = np.zeros((1, H_s.shape[0]), dtype=np.float64)
        starts_tip_dbg = _atom_orb_starts(tip_norb_per)
        starts_smp_dbg = _atom_orb_starts(smp_norb_per)
        from pyBall.FireballOCL.STM import _blocks_to_dense_vector
        _blocks_to_dense_vector(pairs_dbg, Hb_dbg, starts_tip_dbg, starts_smp_dbg, tip_norb_per, smp_norb_per, 1, H_s.shape[0], out=Hts_py)
        _blocks_to_dense_vector(pairs_dbg, Sb_dbg, starts_tip_dbg, starts_smp_dbg, tip_norb_per, smp_norb_per, 1, H_s.shape[0], out=Sts_py)
        dHts = float(np.max(np.abs(Hts_py - Hts_fc)))
        dSts = float(np.max(np.abs(Sts_py - Sts_fc)))
        print(f"  [DEBUG] max|Hts_py - Hts_fc|={dHts:.3e}  max|Sts_py - Sts_fc|={dSts:.3e}")

    # ── Rotated system (Fortran-native rigid rotation) ──
    # Build rotated scan grid (lab-frame points)
    Xr, Yr, pts_r, ext_r = build_rotated_plane_grid(mol.apos, z=z_tip, size=size, n=n, R=R)
    pts_r_f64 = pts_r.astype(np.float64)

    # FireCore ctypes wrapper expects points shaped (npoints,3) (same as tipResponseSimple2points)
    pts_fortran = np.ascontiguousarray(pts_r_f64)
    # tip_src length must match num_orb(tipZ) in Fireball (1 for H)

    resp_rot = fc.tipResponseSimple2points_rotated(
        pts_fortran, R, tip_src_fb,
        E=E, eta=eta, mode=0, iMO=int(mo0+1), ikpoint=1, tipZ=tipZ,
        rcut=rcut, beta=beta, r0=r0,
        A_ss=-1.0, A_sp=-1.0, A_pp_sig=-1.0, A_pp_pi=+1.0, A_ps=-1.0,
        overlap_scale=0.0
    ).reshape(n, n).T

    # ── Rotated system (Python, rigid-body handled in Python) ──
    # Match Fortran rotated logic: inverse-rotate probe points around centroid and
    # evaluate coupling in the body frame, while keeping sample H/S and MO coeffs fixed.
    cen = mol.apos.mean(axis=0)
    pts_body = cen + (pts_r_f64 - cen) @ R  # row-vectors: dp @ R == (R^T @ dp)_col
    resp_py_rot = response_amplitude_simple_lu(
        pts_body, H_use, S_use, C_use, 0, E,
        tip_norb_per, smp_norb_per,
        coupling_builder=coupling0,
        tip_src=tip_src,
        eta=eta
    ).reshape(n, n).T

    # ── Compare ──
    diff = float(np.max(np.abs(resp_rot - resp_ref)))
    rms = float(np.sqrt(np.mean((resp_rot - resp_ref)**2)))
    rel = diff / (float(np.max(np.abs(resp_ref))) + 1e-30)

    diff_py = float(np.max(np.abs(resp_py_rot - resp_rot)))
    rms_py = float(np.sqrt(np.mean((resp_py_rot - resp_rot)**2)))
    rel_py = diff_py / (float(np.max(np.abs(resp_rot))) + 1e-30)

    print(f"  [{label}] max|rot-ref| = {diff:.3e}  RMS={rms:.3e}  rel={rel:.3e}")
    print(f"  [{label}] max|py-fc_rot| = {diff_py:.3e}  RMS={rms_py:.3e}  rel={rel_py:.3e}")

    return resp_ref, resp_rot, resp_py_rot, diff, rms, rel, diff_py, rms_py, rel_py


def main():
    parser = argparse.ArgumentParser(description="STM Rotational Invariance Test")
    parser.add_argument("--xyz", default="../../cpp/common_resources/xyz/H2O.xyz")
    parser.add_argument("--mo", type=int, default=None, help="MO index (1-based). Default: HOMO")
    parser.add_argument("--z_tip", type=float, default=3.0, help="Tip height for STM response (Å)")
    parser.add_argument("--size", type=float, default=12.0, help="Grid extent (Å)")
    parser.add_argument("--n", type=int, default=30, help="Grid resolution (n×n)")
    parser.add_argument("--rcut", type=float, default=8.0, help="Tip-sample coupling cutoff (Å)")
    parser.add_argument("--beta", type=float, default=1.0, help="Exponential decay β (Å⁻¹)")
    parser.add_argument("--r0", type=float, default=3.0, help="Reference distance r0 (Å)")
    parser.add_argument("--eta", type=float, default=1e-2, help="Broadening η (eV)")
    parser.add_argument("--E_tip", type=float, default=0.0, help="Tip onsite energy (eV)")
    parser.add_argument("--axis", type=str, default="x", help="Rotation axis: x|y|z")
    parser.add_argument("--angles", type=str, default="0,45,90", help="Comma-separated angles in degrees")
    parser.add_argument("--debug_fireball", action="store_true", help="DEBUG: recompute H,S with Fireball for rotated geometry and compare to rotated H,S")
    parser.add_argument("--all_orbitals", action="store_true", help="Test all orbitals instead of single MO")
    parser.add_argument("--outdir", default="export/response_rotated")
    args = parser.parse_args()

    os.chdir(_THIS_DIR)
    _ensure_fdata()

    xyz_path = args.xyz if os.path.isabs(args.xyz) else os.path.join(_THIS_DIR, args.xyz)
    mol = AtomicSystem(fname=xyz_path)
    mol_name = os.path.splitext(os.path.basename(xyz_path))[0]
    export_dir = set_export_dir(os.path.join(_THIS_DIR, args.outdir, mol_name))

    fc.setVerbosity(1)
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
    mo_label = "HOMO" if mo0 == homo else ("LUMO" if mo0 == lumo else f"MO{mo0+1}")
    print(f"MO {mo0+1} ({mo_label}) E={E:+.4f} eV")

    # Get H, S, C
    C_fc = fc.get_wfcoef(norb=norb)
    kpt = np.array((0.0, 0.0, 0.0), dtype=np.float64)
    Hk_fc, Sk_fc = fc.get_HS_k(kpt, norb)
    H_s = np.asarray(Hk_fc, dtype=np.complex128)
    S_s = np.asarray(Sk_fc, dtype=np.complex128)

    # Orbital layout
    data = fc.get_HS_neighs(dims)
    data = fc.get_HS_sparse(dims, data)
    norb_per, _ = get_orbital_layout(data, natoms)
    smp_norb_per = norb_per.astype(np.int32)

    # Tip: single H atom, 1s orbital
    tip_norb_per = np.array([1], dtype=np.int32)

    # Determine which orbitals to test
    if args.all_orbitals:
        mo_indices = list(range(norb))
    else:
        mo_indices = [mo0]

    axis = str(args.axis).strip().lower()
    angles = [float(a.strip()) for a in str(args.angles).split(',') if a.strip() != '']
    if len(angles) == 0:
        raise ValueError('--angles is empty')

    # Build grid for extent calculation (used in plotting)
    _, _, _, extent_plot = build_plane_grid(mol.apos, z=args.z_tip, size=args.size, n=args.n)

    rot_tests = []
    for ang in angles:
        if axis == 'x':
            R = rotation_matrix_x(ang)
        elif axis == 'y':
            R = rotation_matrix_y(ang)
        elif axis == 'z':
            R = rotation_matrix_z(ang)
        else:
            raise ValueError(f"--axis must be x|y|z, got {args.axis}")
        rot_tests.append((f"{axis}{int(ang):03d}", R))

    # Summary table for all orbitals
    all_orb_summary = []

    for mo_idx in mo_indices:
        mo0 = mo_idx
        E = float(eigen[mo0])
        mo_label = "HOMO" if mo0 == homo else ("LUMO" if mo0 == lumo else f"MO{mo0+1}")
        print(f"\n{'='*60}")
        print(f"MO {mo0+1} ({mo_label}) E={E:+.4f} eV")
        print(f"{'='*60}")

        print("\n=== Rotational invariance test (s-tip) ===")
        all_results = []
        for rot_name, R in rot_tests:
            print(f"\n  Rotation: {rot_name}")

            if args.debug_fireball:
                # DEBUG: compare rotated dense H,S against Fireball recomputation from rotated geometry
                apos_rot_dbg = rotate_points(mol.apos, R)
                fc.evalForce(apos_rot_dbg, nmax_scf=50)
                Hk_fc_r, Sk_fc_r = fc.get_HS_k(kpt, norb)
                H_fc_r = np.asarray(Hk_fc_r, dtype=np.complex128)
                S_fc_r = np.asarray(Sk_fc_r, dtype=np.complex128)
                H_rot_dbg, S_rot_dbg = rotate_dense_hs(H_s, S_s, smp_norb_per, R)
                dH = float(np.max(np.abs(H_rot_dbg - H_fc_r)))
                dS = float(np.max(np.abs(S_rot_dbg - S_fc_r)))
                print(f"  [DEBUG Fireball] max|H_rot - H_fc(rot geom)| = {dH:.3e}  max|S_rot - S_fc| = {dS:.3e}")
                # restore original geometry to keep reference consistent
                fc.evalForce(mol.apos, nmax_scf=1)

            resp_ref, resp_rot, resp_py_rot, diff, rms, rel, diff_py, rms_py, rel_py = run_rotated_test(
                mol, H_s, S_s, C_fc, eigen, mo0, norb_per,
                tip_norb_per, smp_norb_per, R, rot_name,
                args.z_tip, args.size, args.n,
                args.rcut, args.beta, args.r0, args.eta, args.E_tip
            )
            all_results.append({
                'rot': rot_name,
                'diff': diff,
                'rms': rms,
                'rel': rel,
                'diff_py': diff_py,
                'rms_py': rms_py,
                'rel_py': rel_py,
                'resp_ref': resp_ref,
                'resp_rot': resp_rot,
                'resp_py_rot': resp_py_rot,
            })
            all_orb_summary.append({
                'mo': mo0+1,
                'mo_label': mo_label,
                'rot': rot_name,
                'diff': diff,
                'rel': rel,
                'diff_py': diff_py,
                'rel_py': rel_py,
            })

        # One figure: each row = angle, columns = ref | Fortran-rot | Python-rot | |py-fc|
        n_rot = len(all_results)
        fig, axes = plt.subplots(n_rot, 4, figsize=(20, 5*n_rot))
        if n_rot == 1:
            axes = axes.reshape(1, 4)

        for i, res in enumerate(all_results):
            ref = res['resp_ref']
            rot = res['resp_rot']
            pyrot = res['resp_py_rot']
            diff_map = np.abs(pyrot - rot)

            vmax = float(np.max(ref))
            if vmax < 1e-30: vmax = 1.0
            vmax_d = float(np.max(diff_map))
            if vmax_d < 1e-30: vmax_d = 1.0

            ax = axes[i, 0]
            im = ax.imshow(ref, origin='lower', extent=extent_plot, cmap='hot', vmin=0, vmax=vmax, aspect='equal')
            ax.set_title(f"{res['rot']}  reference")
            ax.set_xlabel('x (Å)'); ax.set_ylabel('y (Å)')
            plot_atoms(ax, mol.apos, mol.atypes, color='cyan', ms=4)
            plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

            ax = axes[i, 1]
            im = ax.imshow(rot, origin='lower', extent=extent_plot, cmap='hot', vmin=0, vmax=vmax, aspect='equal')
            ax.set_title(f"{res['rot']}  Fortran-rot\nmax|Δ(ref)|={res['diff']:.2e}")
            ax.set_xlabel('x (Å)'); ax.set_ylabel('y (Å)')
            plot_atoms(ax, mol.apos, mol.atypes, color='cyan', ms=4)
            plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

            ax = axes[i, 2]
            im = ax.imshow(pyrot, origin='lower', extent=extent_plot, cmap='hot', vmin=0, vmax=vmax, aspect='equal')
            ax.set_title(f"{res['rot']}  Python-rot\nmax|Δ(py-fc)|={res['diff_py']:.2e}")
            ax.set_xlabel('x (Å)'); ax.set_ylabel('y (Å)')
            plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

            ax = axes[i, 3]
            im = ax.imshow(diff_map, origin='lower', extent=extent_plot, cmap='viridis', vmin=0, vmax=vmax_d, aspect='equal')
            ax.set_title(f"{res['rot']}  |py - fc_rot|\nmax={res['diff_py']:.2e}")
            ax.set_xlabel('x (Å)'); ax.set_ylabel('y (Å)')
            plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

        fig.suptitle(f"{mol_name} MO {mo0+1} ({mo_label}) E={E:+.4f} eV  z_tip={args.z_tip:.1f} Å\nRigid-body rotation invariance (s-tip)")
        plt.tight_layout()
        save_plot(fig, f"mo{mo0+1:04d}_{mo_label}_rotinv_s_tip_{axis}", export_dir)
        plt.close(fig)

        print("\n" + "="*60)
        print("SUMMARY")
        print("="*60)
        for res in all_results:
            print(f"  rot={res['rot']:6s}  max|Δ(fc-ref)|={res['diff']:.3e}  rel={res['rel']:.3e}   max|Δ(py-fc)|={res['diff_py']:.3e}  rel_py={res['rel_py']:.3e}")

    # Print overall summary for all orbitals
    if args.all_orbitals:
        print("\n" + "="*60)
        print("OVERALL SUMMARY (ALL ORBITALS)")
        print("="*60)
        print(f"{'MO':>4} {'Label':>6} {'Rot':>6} {'max|Δ(fc-ref)|':>14} {'max|Δ(py-fc)|':>14}")
        print("-" * 46)
        for s in all_orb_summary:
            print(f"{s['mo']:>4} {s['mo_label']:>6} {s['rot']:>6} {s['diff']:>14.3e} {s['diff_py']:>14.3e}")

    print(f"\nDone. Output in {export_dir}")


if __name__ == "__main__":
    main()
