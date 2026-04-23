#!/usr/bin/env python3

import os, sys, argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.path.dirname(__file__), "..", ".."))

from pyBall.AtomicSystem import AtomicSystem
from pyBall import FireCore as fc

from pyBall.FireballOCL.STM_utils import (
    set_export_dir, save_plot, plot_atoms,
    project_orbital_to_points,
    get_orbital_layout, build_plane_grid,
)

from pyBall.FireballOCL.STM import (
    rotate_mo_coefficients, rotate_points,
    response_amplitude_simple_lu, make_coupling_builder_molecular_tip_exp_sk,
)
from pyBall.FireballOCL.Grid import GridProjector
from pyBall.FireballOCL.STM_utils import remap_coeffs_fortran_to_grid

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))


def _ensure_fdata():
    dst = os.path.join(_THIS_DIR, "Fdata")
    src_rel = "../Fireball/Fdata_HCNOS"
    if not (os.path.lexists(dst) and os.path.realpath(dst).endswith("Fdata_HCNOS")):
        if os.path.lexists(dst):
            os.unlink(dst)
        os.symlink(src_rel, dst)


def _homo_index(eigen):
    occ = np.where(np.asarray(eigen) < 0.0)[0]
    if len(occ) == 0:
        return max(0, len(eigen)//2 - 1)
    return int(occ[-1])


def rotation_matrix(axis, angle_deg):
    th = np.radians(float(angle_deg))
    c, s = np.cos(th), np.sin(th)
    if axis == 'x':
        return np.array([[1,0,0],[0,c,-s],[0,s,c]], dtype=np.float64)
    if axis == 'y':
        return np.array([[c,0,s],[0,1,0],[-s,0,c]], dtype=np.float64)
    if axis == 'z':
        return np.array([[c,-s,0],[s,c,0],[0,0,1]], dtype=np.float64)
    raise ValueError(f"axis must be x|y|z, got {axis}")


def _ao_local_index_fortran(ao):
    # Fireball/Fortran order on sp atoms: [s, py, pz, px]
    if ao == 's':  return 0
    if ao == 'py': return 1
    if ao == 'pz': return 2
    if ao == 'px': return 3
    raise ValueError(f"ao must be s|px|py|pz, got {ao}")


def build_mock_C(norb, starts, norb_per, atom, ao):
    C = np.zeros((1, norb), dtype=np.float64)
    ia = int(atom)
    no = int(norb_per[ia])
    if no == 1:
        if ao != 's':
            raise ValueError(f"mock orbital on atom {ia}: only s is available (norb_per=1), got ao={ao}")
        C[0, int(starts[ia]) + 0] = 1.0
        return C
    if no == 4:
        il = _ao_local_index_fortran(ao)
        C[0, int(starts[ia]) + il] = 1.0
        return C
    raise ValueError(f"mock orbital: unsupported norb_per[{ia}]={no} (expected 1 or 4)")


def quat_from_axis_angle(axis, angle_deg):
    th = np.radians(float(angle_deg))
    s = np.sin(0.5*th)
    c = np.cos(0.5*th)
    if axis == 'x':
        return np.array([s, 0.0, 0.0, c], dtype=np.float32)
    if axis == 'y':
        return np.array([0.0, s, 0.0, c], dtype=np.float32)
    if axis == 'z':
        return np.array([0.0, 0.0, s, c], dtype=np.float32)
    raise ValueError(f"axis must be x|y|z, got {axis}")


def get_orbital_layout_from_fireball(dims):
    data = fc.get_HS_neighs(dims)
    data = fc.get_HS_sparse(dims, data)
    n_orb_atom, starts = get_orbital_layout(data, dims.natoms)
    orb2atom = np.array([ia for ia in range(dims.natoms) for _ in range(int(n_orb_atom[ia]))], dtype=np.int32)
    return n_orb_atom.astype(np.int32), starts, orb2atom

'''
EXAMPLE OF USAGE:

# Mock AO mode with different tip/sample orbitals and rotations:
python3 test_stm_orbital_rotated.py --do_overlap --mock_ao --mock_tip_atom 1 --mock_smp_atom 1 --mock_tip_ao px --mock_smp_ao py --tip_axis z --tip_angles 0,45,90 --smp_axis z --smp_angles 0 --n 80 --size 10.0 --ztip 3.0 --zmid 1.5 --beta 1.0 --r0 3.0 --rcut 8.0 --nmax_scf 30

# True orbitals with different tip and sample MOs:
python3 test_stm_orbital_rotated.py --do_overlap --tip_mo 1 --smp_mo 2 --tip_axis z --tip_angles 0,45,90 --smp_axis z --smp_angles 0 --n 80 --size 10.0 --ztip 3.0 --zmid 1.5 --beta 1.0 --r0 3.0 --rcut 8.0 --nmax_scf 30

# True orbitals with same MO for both (backward compatible):
python3 test_stm_orbital_rotated.py --do_overlap --mo_list 1 --tip_axis z --tip_angles 0,45,90 --smp_axis z --smp_angles 0 --n 80 --size 10.0 --ztip 3.0 --zmid 1.5 --beta 1.0 --r0 3.0 --rcut 8.0 --nmax_scf 30
'''


def main():
    ap = argparse.ArgumentParser(description='CH2O rotated-orbital projection test')
    ap.add_argument('--xyz', default='../../cpp/common_resources/xyz/CH2O.xyz')
    ap.add_argument('--mo', type=int, default=None, help='MO index 1-based for both tip and sample; default HOMO-1')
    ap.add_argument('--mo_list', default='', help='Comma-separated list of MO indices (1-based) for both tip and sample; overrides --mo')
    ap.add_argument('--tip_mo', type=int, default=None, help='Tip MO index 1-based; overrides --mo if set')
    ap.add_argument('--tip_mo_list', default='', help='Comma-separated tip MO indices (1-based); overrides --tip_mo')
    ap.add_argument('--smp_mo', type=int, default=None, help='Sample MO index 1-based; overrides --mo if set')
    ap.add_argument('--smp_mo_list', default='', help='Comma-separated sample MO indices (1-based); overrides --smp_mo')
    ap.add_argument('--z', type=float, default=2.0)
    ap.add_argument('--size', type=float, default=10.0)
    ap.add_argument('--n', type=int, default=120)
    ap.add_argument('--axis', default='x', help='Rotation axis (deprecated; use --tip_axis)')
    ap.add_argument('--angle', type=float, default=90.0, help='Rotation angle (deprecated; use --tip_angles)')
    ap.add_argument('--angles', default='', help='Comma-separated rotation angles (deg) to run; overrides --angle (deprecated)')
    ap.add_argument('--tip_axis', default='x', help='Tip rotation axis')
    ap.add_argument('--tip_angle', type=float, default=90.0, help='Tip rotation angle (deg)')
    ap.add_argument('--tip_angles', default='', help='Comma-separated tip rotation angles (deg) to run; overrides --tip_angle')
    ap.add_argument('--smp_axis', default='x', help='Sample rotation axis (default identity: no rotation)')
    ap.add_argument('--smp_angle', type=float, default=0.0, help='Sample rotation angle (deg)')
    ap.add_argument('--smp_angles', default='', help='Comma-separated sample rotation angles (deg) to run; overrides --smp_angle')
    ap.add_argument('--nmax_scf', type=int, default=200)
    ap.add_argument('--outdir', default='export/stm_orbital_rotated')
    ap.add_argument('--do_proj3', action='store_true', help='DEBUG: plot 3-panel unrot/rot/diff orbital projection')
    ap.add_argument('--do_response', action='store_true', help='Compute STM response map with H-tip')
    ap.add_argument('--do_overlap', action='store_true', help='Compute two-molecule overlap/correlation + convolution debug maps')
    ap.add_argument('--mock_ao', action='store_true', help='Use simple mock orbitals: single AO on chosen atom(s) to debug rotation')
    ap.add_argument('--mock_tip_atom', type=int, default=0, help='Tip atom index (0-based) for --mock_ao')
    ap.add_argument('--mock_smp_atom', type=int, default=0, help='Sample atom index (0-based) for --mock_ao')
    ap.add_argument('--mock_tip_ao', default='px', help='Tip AO for --mock_ao: s|px|py|pz')
    ap.add_argument('--mock_smp_ao', default='px', help='Sample AO for --mock_ao: s|px|py|pz')
    ap.add_argument('--ztip', type=float, default=3.0, help='Tip molecule center shift in z (Å), sample stays at z=0')
    ap.add_argument('--zmid', type=float, default=None, help='Mid imaging plane z offset from sample center (Å); default ztip/2')
    ap.add_argument('--eta', type=float, default=1e-2)
    ap.add_argument('--rcut', type=float, default=8.0)
    ap.add_argument('--beta', type=float, default=1.0)
    ap.add_argument('--r0', type=float, default=3.0)
    args = ap.parse_args()

    os.chdir(_THIS_DIR)
    _ensure_fdata()

    xyz_path = args.xyz if os.path.isabs(args.xyz) else os.path.join(_THIS_DIR, args.xyz)
    mol = AtomicSystem(fname=xyz_path)
    mol_name = os.path.splitext(os.path.basename(xyz_path))[0]

    export_dir = set_export_dir(os.path.join(_THIS_DIR, args.outdir, mol_name.lower()))
    fdata_basis_dir = os.path.join(_THIS_DIR, 'Fdata', 'basis')
    fdata_dir = os.path.join(_THIS_DIR, 'Fdata')

    fc.setVerbosity(0)
    fc.initialize(atomType=mol.atypes, atomPos=mol.apos, verbosity=0)
    dims = fc.get_HS_dims()
    norb = int(dims.norbitals)
    fc.evalForce(mol.apos, nmax_scf=args.nmax_scf)
    eigen = fc.get_eigen(ikp=1, norb=norb)
    C_fc = fc.get_wfcoef(norb=norb)

    norb_per, starts, orb2atom = get_orbital_layout_from_fireball(dims)

    homo = _homo_index(eigen)

    # Parse tip and sample MO lists (support both old and new CLI)
    if args.tip_mo_list.strip():
        tip_mo_list = [int(s)-1 for s in args.tip_mo_list.split(',') if s.strip()]
    elif args.tip_mo is not None:
        tip_mo_list = [args.tip_mo-1]
    elif args.mo_list.strip():
        tip_mo_list = [int(s)-1 for s in args.mo_list.split(',') if s.strip()]
    elif args.mo is not None:
        tip_mo_list = [args.mo-1]
    else:
        mo0 = max(0, homo-1)
        tip_mo_list = [int(mo0)]

    if args.smp_mo_list.strip():
        smp_mo_list = [int(s)-1 for s in args.smp_mo_list.split(',') if s.strip()]
    elif args.smp_mo is not None:
        smp_mo_list = [args.smp_mo-1]
    elif args.mo_list.strip():
        smp_mo_list = [int(s)-1 for s in args.mo_list.split(',') if s.strip()]
    elif args.mo is not None:
        smp_mo_list = [args.mo-1]
    else:
        mo0 = max(0, homo-1)
        smp_mo_list = [int(mo0)]

    tip_mo_list = [int(np.clip(m, 0, norb-1)) for m in tip_mo_list]
    smp_mo_list = [int(np.clip(m, 0, norb-1)) for m in smp_mo_list]

    # Parse tip and sample rotation parameters (support both old and new CLI)
    tip_axis = args.tip_axis
    if args.tip_angles.strip():
        tip_ang_list = [float(s) for s in args.tip_angles.split(',') if s.strip()]
    else:
        tip_ang_list = [float(args.tip_angle)]
    smp_axis = args.smp_axis
    if args.smp_angles.strip():
        smp_ang_list = [float(s) for s in args.smp_angles.split(',') if s.strip()]
    else:
        smp_ang_list = [float(args.smp_angle)]
    # Backward compatibility: if old --axis/--angles used, apply to tip
    if args.angles.strip():
        tip_ang_list = [float(s) for s in args.angles.split(',') if s.strip()]
        tip_axis = args.axis
    elif args.angle != 90.0:
        tip_ang_list = [float(args.angle)]
        tip_axis = args.axis

    if args.do_proj3:
        X, Y, pts, extent = build_plane_grid(mol.apos, z=args.z, size=args.size, n=args.n)
        for ang in tip_ang_list:
            R = rotation_matrix(tip_axis, ang)
            apos_rot = rotate_points(mol.apos, R)
            if args.mock_ao:
                C_mock_tip = build_mock_C(norb, starts, norb_per, args.mock_tip_atom, args.mock_tip_ao)
                C_mock_smp = build_mock_C(norb, starts, norb_per, args.mock_smp_atom, args.mock_smp_ao)
                C_tip0 = C_mock_tip
                C_smp0 = C_mock_smp
            else:
                C_tip0 = C_fc
                C_smp0 = C_fc
            C_tip_rot = rotate_mo_coefficients(C_tip0.T, norb_per, R).T
            for mo0 in tip_mo_list:
                E = float(eigen[mo0])

                psi0 = project_orbital_to_points(C_smp0,    mo0, mol.atypes, mol.apos,   orb2atom, norb_per, fdata_basis_dir, pts.astype(np.float32))
                psi1 = project_orbital_to_points(C_tip_rot, mo0, mol.atypes, apos_rot,   orb2atom, norb_per, fdata_basis_dir, pts.astype(np.float32))

                psi0 = np.asarray(psi0, dtype=np.float64).reshape(args.n, args.n)
                psi1 = np.asarray(psi1, dtype=np.float64).reshape(args.n, args.n)

                vmax = max(np.max(np.abs(psi0)), np.max(np.abs(psi1)))
                if vmax < 1e-30: vmax = 1.0

                fig, ax = plt.subplots(1, 3, figsize=(18, 6))
                if args.mock_ao:
                    fig.suptitle(f"{mol_name} MO{mo0+1} E={E:+.4f} eV  z={args.z}  R={tip_axis}{ang}  (mock AO: tip={args.mock_tip_atom}:{args.mock_tip_ao} smp={args.mock_smp_atom}:{args.mock_smp_ao})")
                else:
                    fig.suptitle(f"{mol_name} MO{mo0+1} E={E:+.4f} eV  z={args.z}  R={tip_axis}{ang}")

                im0 = ax[0].imshow(psi0.T, origin='lower', extent=extent, cmap='bwr', vmin=-vmax, vmax=vmax, aspect='equal')
                ax[0].set_title('unrotated ψ'); plot_atoms(ax[0], mol.apos, mol.atypes)
                fig.colorbar(im0, ax=ax[0], fraction=0.046, pad=0.04)

                im1 = ax[1].imshow(psi1.T, origin='lower', extent=extent, cmap='bwr', vmin=-vmax, vmax=vmax, aspect='equal')
                ax[1].set_title('rotated ψ'); plot_atoms(ax[1], apos_rot, mol.atypes)
                fig.colorbar(im1, ax=ax[1], fraction=0.046, pad=0.04)

                im2 = ax[2].imshow((psi1-psi0).T, origin='lower', extent=extent, cmap='bwr', aspect='equal')
                ax[2].set_title('difference (rot-unrot)'); plot_atoms(ax[2], mol.apos, mol.atypes)
                fig.colorbar(im2, ax=ax[2], fraction=0.046, pad=0.04)

                save_plot(fig, f"orbital_projection_mo{mo0+1:03d}_z{args.z:.2f}_{tip_axis}{int(ang):03d}", export_dir=export_dir)
                plt.close(fig)

    if args.do_overlap:
        zmid = (0.5*args.ztip) if (args.zmid is None) else float(args.zmid)
        cen0 = np.asarray(mol.apos.mean(axis=0), dtype=np.float64)
        smp_pos = np.asarray(mol.apos, dtype=np.float64) + np.array([0.0, 0.0, -cen0[2]], dtype=np.float64)  # center at z=0

        # Build mid-plane grid centered on sample in xy at absolute z=zmid
        xs = np.linspace(cen0[0] - args.size*0.5, cen0[0] + args.size*0.5, args.n)
        ys = np.linspace(cen0[1] - args.size*0.5, cen0[1] + args.size*0.5, args.n)
        X2, Y2 = np.meshgrid(xs, ys, indexing='ij')
        pts_mid = np.zeros((args.n*args.n, 3), dtype=np.float64)
        pts_mid[:, 0] = X2.ravel(); pts_mid[:, 1] = Y2.ravel(); pts_mid[:, 2] = zmid
        extent2 = [xs[0], xs[-1], ys[0], ys[-1]]

        # Axes in shift coordinates (Å)
        dx = (xs[1] - xs[0])
        dy = (ys[1] - ys[0])
        sx = (np.arange(args.n) - args.n//2) * dx
        sy = (np.arange(args.n) - args.n//2) * dy
        extent_shift = [sx[0], sx[-1], sy[0], sy[-1]]

        gp = GridProjector(fdata_dir=fdata_dir)

        # Nested loops: tip rotations x sample rotations
        for tip_ang in tip_ang_list:
            R_tip = rotation_matrix(tip_axis, tip_ang)
            apos_tip_rot = rotate_points(mol.apos, R_tip)
            if args.mock_ao:
                C_mock_tip = build_mock_C(norb, starts, norb_per, args.mock_tip_atom, args.mock_tip_ao)
                C_tip0 = C_mock_tip
            else:
                C_tip0 = C_fc
            C_tip_rot = rotate_mo_coefficients(C_tip0.T, norb_per, R_tip).T

            # Tip: rotated geometry, then shifted to z=ztip
            tip_pos_rot = np.asarray(apos_tip_rot, dtype=np.float64) + np.array([0.0, 0.0, -cen0[2] + args.ztip], dtype=np.float64)

            # Define rigid tip geometry in *local unrotated* coordinates (relative to tip center)
            tip_pos_rel0 = np.asarray(mol.apos, dtype=np.float32)
            tip_pos_rel0 = tip_pos_rel0 - tip_pos_rel0.mean(axis=0)

            # Tip-center scan points: keep z fixed at ztip, scan over the same xy grid
            tip_centers = pts_mid.copy().astype(np.float32)
            tip_centers[:, 2] = float(args.ztip)

            # Per-pixel rotation quaternions (tip rotation only)
            q_tip = quat_from_axis_angle(tip_axis, tip_ang)
            tip_quat = np.tile(q_tip[None, :], (tip_centers.shape[0], 1)).astype(np.float32)

            for smp_ang in smp_ang_list:
                R_smp = rotation_matrix(smp_axis, smp_ang)
                apos_smp_rot = rotate_points(mol.apos, R_smp)
                smp_pos = np.asarray(apos_smp_rot, dtype=np.float64) + np.array([0.0, 0.0, -cen0[2]], dtype=np.float64)

                if args.mock_ao:
                    C_mock_smp = build_mock_C(norb, starts, norb_per, args.mock_smp_atom, args.mock_smp_ao)
                    C_smp0 = C_mock_smp
                else:
                    C_smp0 = C_fc
                C_smp_rot = rotate_mo_coefficients(C_smp0.T, norb_per, R_smp).T

                for tip_mo in tip_mo_list:
                    for smp_mo in smp_mo_list:
                        psi_s = project_orbital_to_points(C_smp_rot, smp_mo, mol.atypes, smp_pos, orb2atom, norb_per, fdata_basis_dir, pts_mid.astype(np.float32))
                        psi_t = project_orbital_to_points(C_tip_rot, tip_mo, mol.atypes, tip_pos_rot, orb2atom, norb_per, fdata_basis_dir, pts_mid.astype(np.float32))
                        psi_s = np.asarray(psi_s, dtype=np.float64).reshape(args.n, args.n)
                        psi_t = np.asarray(psi_t, dtype=np.float64).reshape(args.n, args.n)

                        Fs = np.fft.fft2(psi_s)
                        Ft = np.fft.fft2(psi_t)
                        corr = np.fft.ifft2(Fs * np.conj(Ft)).real
                        corr = np.fft.fftshift(corr)
                        corr2 = corr * corr

                        # IMPORTANT: kernel applies rotation per pixel using tip_quat, so we pass unrotated tip coefficients here
                        coeff_tip0 = remap_coeffs_fortran_to_grid(C_tip0[tip_mo, :], norb_per)
                        coeff_smp0 = remap_coeffs_fortran_to_grid(C_smp0[smp_mo, :], norb_per)
                        t_amp, t_I = gp.mo_overlap_points_exp_sk(
                            tip_centers=tip_centers,
                            tip_quat=tip_quat,
                            tip_pos_rel=tip_pos_rel0,
                            smp_pos=smp_pos.astype(np.float32),
                            coeffs_tip=coeff_tip0.astype(np.float32),
                            coeffs_smp=coeff_smp0.astype(np.float32),
                            beta=float(args.beta), r0=float(args.r0), rcut=float(args.rcut),
                        )
                        t_I = np.asarray(t_I, dtype=np.float64).reshape(args.n, args.n)

                        fig, ax = plt.subplots(2, 2, figsize=(14, 12))
                        if args.mock_ao:
                            fig.suptitle(f"Overlap-STM vs xcorr (mock AO) sample={args.mock_smp_atom}:{args.mock_smp_ao} tip={args.mock_tip_atom}:{args.mock_tip_ao} zmid={zmid:.2f} ztip={args.ztip:.2f} R_tip={tip_axis}{int(tip_ang)} R_smp={smp_axis}{int(smp_ang)}")
                        else:
                            if tip_mo == smp_mo:
                                fig.suptitle(f"Overlap-STM vs xcorr (mid-plane) MO{tip_mo+1} zmid={zmid:.2f} ztip={args.ztip:.2f} R_tip={tip_axis}{int(tip_ang)} R_smp={smp_axis}{int(smp_ang)}")
                            else:
                                fig.suptitle(f"Overlap-STM vs xcorr (mid-plane) tip MO{tip_mo+1} sample MO{smp_mo+1} zmid={zmid:.2f} ztip={args.ztip:.2f} R_tip={tip_axis}{int(tip_ang)} R_smp={smp_axis}{int(smp_ang)}")

                        vmax_s = float(np.max(np.abs(psi_s)))
                        if vmax_s < 1e-30: vmax_s = 1.0
                        im = ax[0,0].imshow(psi_s.T, origin='lower', extent=extent2, cmap='bwr', vmin=-vmax_s, vmax=vmax_s, aspect='equal')
                        ax[0,0].set_title('1) ψ_sample(mid)'); ax[0,0].set_xlabel('x (Å)'); ax[0,0].set_ylabel('y (Å)')
                        plot_atoms(ax[0,0], smp_pos, mol.atypes)
                        fig.colorbar(im, ax=ax[0,0], fraction=0.046, pad=0.04)

                        vmax_t = float(np.max(np.abs(psi_t)))
                        if vmax_t < 1e-30: vmax_t = 1.0
                        im = ax[0,1].imshow(psi_t.T, origin='lower', extent=extent2, cmap='bwr', vmin=-vmax_t, vmax=vmax_t, aspect='equal')
                        ax[0,1].set_title('2) ψ_tip(mid)'); ax[0,1].set_xlabel('x (Å)'); ax[0,1].set_ylabel('y (Å)')
                        plot_atoms(ax[0,1], tip_pos_rot, mol.atypes)
                        fig.colorbar(im, ax=ax[0,1], fraction=0.046, pad=0.04)

                        vI = float(np.max(t_I))
                        if vI <= 0: vI = 1.0
                        im = ax[1,0].imshow(t_I.T, origin='lower', extent=extent2, cmap='viridis', vmin=0.0, vmax=vI, aspect='equal')
                        ax[1,0].set_title('3) STM overlap |t|^2 (LCAO shift tip)'); ax[1,0].set_xlabel('x (Å)'); ax[1,0].set_ylabel('y (Å)')
                        plot_atoms(ax[1,0], smp_pos, mol.atypes)
                        fig.colorbar(im, ax=ax[1,0], fraction=0.046, pad=0.04)

                        vC = float(np.max(corr2))
                        if vC <= 0: vC = 1.0
                        im = ax[1,1].imshow(corr2.T, origin='lower', extent=extent_shift, cmap='viridis', vmin=0.0, vmax=vC, aspect='equal')
                        ax[1,1].set_title('4) xcorr(ψ_sample,ψ_tip)^2'); ax[1,1].set_xlabel('dx (Å)'); ax[1,1].set_ylabel('dy (Å)')
                        fig.colorbar(im, ax=ax[1,1], fraction=0.046, pad=0.04)

                        save_plot(fig, f"overlap_stm_vs_xcorr_tipMO{tip_mo+1:03d}_smpMO{smp_mo+1:03d}_zmid{zmid:.2f}_ztip{args.ztip:.2f}_tip{tip_axis}{int(tip_ang):03d}_smp{smp_axis}{int(smp_ang):03d}", export_dir=export_dir)
                        plt.close(fig)

    if args.do_response:
        # Dense H,S for sample response (wrapper already handles transpose)
        kpt = np.array((0.0, 0.0, 0.0), dtype=np.float64)
        H_s, S_s = fc.get_HS_k(kpt, norb)
        H_s = np.asarray(H_s, dtype=np.complex128)
        S_s = np.asarray(S_s, dtype=np.complex128)

        # Sample MO coefficients must be (norb, nmo)
        C_s = np.asarray(C_fc.T, dtype=np.float64)
        C_s_rot = np.asarray(C_rot.T, dtype=np.float64)

        # H-tip model (single s orbital)
        tip_types = np.array([1], dtype=np.int32)
        tip_pos_rel = np.array([[0.0, 0.0, 0.0]], dtype=np.float64)
        tip_norb_per = np.array([1], dtype=np.int32)

        smp_types = np.asarray(mol.atypes, dtype=np.int32)
        smp_pos0 = np.asarray(mol.apos, dtype=np.float64)
        smp_pos1 = np.asarray(apos_rot, dtype=np.float64)
        smp_norb_per = np.asarray(norb_per, dtype=np.int32)

        coupling0 = make_coupling_builder_molecular_tip_exp_sk(
            tip_types, tip_pos_rel, tip_norb_per,
            smp_types, smp_pos0, smp_norb_per,
            rcut=args.rcut, beta=args.beta, r0=args.r0,
            overlap_scale=0.0
        )
        coupling1 = make_coupling_builder_molecular_tip_exp_sk(
            tip_types, tip_pos_rel, tip_norb_per,
            smp_types, smp_pos1, smp_norb_per,
            rcut=args.rcut, beta=args.beta, r0=args.r0,
            overlap_scale=0.0
        )

        tip_src = np.array([1.0+0j], dtype=np.complex128)
        resp0 = response_amplitude_simple_lu(pts, H_s, S_s, C_s, mo0, E, tip_norb_per, smp_norb_per, coupling0, tip_src=tip_src, eta=args.eta)
        resp1 = response_amplitude_simple_lu(pts, H_s, S_s, C_s_rot, mo0, E, tip_norb_per, smp_norb_per, coupling1, tip_src=tip_src, eta=args.eta)
        resp0 = np.asarray(resp0, dtype=np.float64).reshape(args.n, args.n)
        resp1 = np.asarray(resp1, dtype=np.float64).reshape(args.n, args.n)

        fig, ax = plt.subplots(1, 3, figsize=(18, 6))
        fig.suptitle(f"STM response (H-tip) {mol_name} MO{mo0+1} E={E:+.4f} eV  eta={args.eta:g}  z={args.z}  R={args.axis}{args.angle}")

        v0 = float(np.max(resp0))
        v1 = float(np.max(resp1))
        vmax = max(v0, v1)
        if vmax <= 0: vmax = 1.0

        im0 = ax[0].imshow(resp0.T, origin='lower', extent=extent, cmap='viridis', vmin=0.0, vmax=vmax, aspect='equal')
        ax[0].set_title('unrotated resp'); plot_atoms(ax[0], mol.apos, mol.atypes)
        fig.colorbar(im0, ax=ax[0], fraction=0.046, pad=0.04)

        im1 = ax[1].imshow(resp1.T, origin='lower', extent=extent, cmap='viridis', vmin=0.0, vmax=vmax, aspect='equal')
        ax[1].set_title('rotated resp'); plot_atoms(ax[1], apos_rot, mol.atypes)
        fig.colorbar(im1, ax=ax[1], fraction=0.046, pad=0.04)

        im2 = ax[2].imshow((resp1-resp0).T, origin='lower', extent=extent, cmap='bwr', aspect='equal')
        ax[2].set_title('difference (rot-unrot)'); plot_atoms(ax[2], mol.apos, mol.atypes)
        fig.colorbar(im2, ax=ax[2], fraction=0.046, pad=0.04)

        save_plot(fig, f"response_map_mo{mo0+1:03d}_z{args.z:.2f}_{args.axis}{int(args.angle):03d}", export_dir=export_dir)
        plt.close(fig)

    print('DONE')
    print('export_dir=', export_dir)


if __name__ == '__main__':
    main()
