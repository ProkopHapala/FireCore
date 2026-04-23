#!/usr/bin/env python3

import os, sys, argparse, subprocess
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
    get_orbital_layout,
)

from pyBall.FireballOCL.STM import (
    rotate_mo_coefficients, rotate_points,
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


def run_scf_for_mol(mol, nmax_scf=200):
    fc.setVerbosity(0)
    fc.initialize(atomType=mol.atypes, atomPos=mol.apos, verbosity=0)
    dims = fc.get_HS_dims()
    norb = int(dims.norbitals)
    fc.evalForce(mol.apos, nmax_scf=nmax_scf)
    eigen = fc.get_eigen(ikp=1, norb=norb)
    C_fc = fc.get_wfcoef(norb=norb)
    norb_per, starts, orb2atom = get_orbital_layout_from_fireball(dims)
    return {
        'dims': dims,
        'norb': norb,
        'eigen': np.array(eigen, dtype=np.float64),
        'C': np.array(C_fc, dtype=np.float64),
        'norb_per': np.array(norb_per, dtype=np.int32),
        'starts': np.array(starts, dtype=np.int32),
        'orb2atom': np.array(orb2atom, dtype=np.int32),
    }


def save_scf_cache(fname, mol, scf):
    np.savez(
        fname,
        atypes=np.asarray(mol.atypes, dtype=np.int32),
        apos=np.asarray(mol.apos, dtype=np.float64),
        eigen=np.asarray(scf['eigen'], dtype=np.float64),
        C=np.asarray(scf['C'], dtype=np.float64),
        norb_per=np.asarray(scf['norb_per'], dtype=np.int32),
        starts=np.asarray(scf['starts'], dtype=np.int32),
        orb2atom=np.asarray(scf['orb2atom'], dtype=np.int32),
        norb=np.int32(scf['norb']),
        natoms=np.int32(int(len(mol.atypes))),
    )


def load_scf_cache(fname):
    d = np.load(fname)
    out = {
        'norb': int(d['norb']),
        'eigen': np.asarray(d['eigen'], dtype=np.float64),
        'C': np.asarray(d['C'], dtype=np.float64),
        'norb_per': np.asarray(d['norb_per'], dtype=np.int32),
        'starts': np.asarray(d['starts'], dtype=np.int32),
        'orb2atom': np.asarray(d['orb2atom'], dtype=np.int32),
        'atypes': np.asarray(d['atypes'], dtype=np.int32),
        'apos': np.asarray(d['apos'], dtype=np.float64),
    }
    return out


def ensure_scf_cache(xyz_path, cache_path, nmax_scf):
    if os.path.isfile(cache_path):
        return
    cmd = [
        sys.executable, os.path.abspath(__file__),
        '--dump_scf',
        '--xyz_one', os.path.abspath(xyz_path),
        '--cache_one', os.path.abspath(cache_path),
        '--nmax_scf', str(int(nmax_scf)),
    ]
    res = subprocess.run(cmd)
    if res.returncode != 0:
        raise RuntimeError(f"SCF cache subprocess failed (returncode={res.returncode}): {' '.join(cmd)}")


def _parse_mo_list(s, one_based=True):
    if not s.strip():
        return []
    out = [int(x) for x in s.split(',') if x.strip()]
    if one_based:
        out = [x-1 for x in out]
    return out

'''
python3 tests/pyFireball/test_stm_orbital_rotated_2mol.py \
  --xyz_smp ../../cpp/common_resources/xyz/PTCDA.xyz \
  --xyz_tip ../../cpp/common_resources/xyz/CH2O.xyz \
  --smp_mo 70 \
  --tip_axis z --tip_angles 0,30,45,60,90 \
  --smp_axis z --smp_angles 0 \
  --n 80 --size 20.0 --ztip 3.0 --zmid 1.5 \
  --beta 1.0 --r0 3.0 --rcut 8.0 \
  --nmax_scf 30


python3 tests/pyFireball/test_stm_orbital_rotated_2mol.py --xyz_smp ../../cpp/common_resources/xyz/PTCDA.xyz --xyz_tip ../../cpp/common_resources/xyz/CH2O.xyz --smp_mo 70 --tip_mo 2 --tip_axis x --tip_angles 0,30,45,60,90 --smp_axis z --smp_angles 0 --n 80 --size 20.0 --ztip 3.0 --zmid 1.5 --beta 1.0 --r0 3.0 --rcut 8.0 --nmax_scf 30

'''
def main():
    ap = argparse.ArgumentParser(description='Two-molecule STM overlap test (tip molecule != sample molecule)')
    ap.add_argument('--dump_scf', action='store_true', help='Internal: run SCF for one molecule and dump cache')
    ap.add_argument('--xyz_one', default='', help='Internal: XYZ for --dump_scf')
    ap.add_argument('--cache_one', default='', help='Internal: output .npz for --dump_scf')
    ap.add_argument('--xyz_tip', default='../../cpp/common_resources/xyz/CH2O.xyz')
    ap.add_argument('--xyz_smp', default='../../cpp/common_resources/xyz/PTCDA.xyz')

    ap.add_argument('--tip_mo', type=int, default=None, help='Tip MO index 1-based; default HOMO-1')
    ap.add_argument('--tip_mo_list', default='', help='Comma-separated tip MO indices (1-based); overrides --tip_mo')
    ap.add_argument('--smp_mo', type=int, default=None, help='Sample MO index 1-based; default HOMO-1')
    ap.add_argument('--smp_mo_list', default='', help='Comma-separated sample MO indices (1-based); overrides --smp_mo')

    ap.add_argument('--tip_axis', default='z')
    ap.add_argument('--tip_angle', type=float, default=0.0)
    ap.add_argument('--tip_angles', default='')
    ap.add_argument('--smp_axis', default='z')
    ap.add_argument('--smp_angle', type=float, default=0.0)
    ap.add_argument('--smp_angles', default='')

    ap.add_argument('--nmax_scf', type=int, default=200)
    ap.add_argument('--outdir', default='export/stm_orbital_rotated_2mol')

    ap.add_argument('--n', type=int, default=80)
    ap.add_argument('--size', type=float, default=20.0)
    ap.add_argument('--ztip', type=float, default=3.0)
    ap.add_argument('--zmid', type=float, default=None)
    ap.add_argument('--rcut', type=float, default=8.0)
    ap.add_argument('--beta', type=float, default=1.0)
    ap.add_argument('--r0', type=float, default=3.0)

    args = ap.parse_args()

    os.chdir(_THIS_DIR)
    _ensure_fdata()

    fdata_basis_dir = os.path.join(_THIS_DIR, 'Fdata', 'basis')
    fdata_dir = os.path.join(_THIS_DIR, 'Fdata')

    if args.dump_scf:
        if not args.xyz_one:
            raise ValueError('--dump_scf requires --xyz_one')
        if not args.cache_one:
            raise ValueError('--dump_scf requires --cache_one')
        mol = AtomicSystem(fname=args.xyz_one)
        scf = run_scf_for_mol(mol, nmax_scf=args.nmax_scf)
        save_scf_cache(args.cache_one, mol, scf)
        return

    xyz_tip = args.xyz_tip if os.path.isabs(args.xyz_tip) else os.path.join(_THIS_DIR, args.xyz_tip)
    xyz_smp = args.xyz_smp if os.path.isabs(args.xyz_smp) else os.path.join(_THIS_DIR, args.xyz_smp)

    tip_name = os.path.splitext(os.path.basename(xyz_tip))[0]
    smp_name = os.path.splitext(os.path.basename(xyz_smp))[0]

    export_dir = set_export_dir(os.path.join(_THIS_DIR, args.outdir, f"tip_{tip_name.lower()}__smp_{smp_name.lower()}"))

    cache_tip = os.path.join(export_dir, f"scf_tip_{tip_name}.npz")
    cache_smp = os.path.join(export_dir, f"scf_smp_{smp_name}.npz")
    ensure_scf_cache(xyz_tip, cache_tip, args.nmax_scf)
    ensure_scf_cache(xyz_smp, cache_smp, args.nmax_scf)

    tip_scf = load_scf_cache(cache_tip)
    smp_scf = load_scf_cache(cache_smp)

    mol_tip = AtomicSystem()
    mol_tip.atypes = tip_scf['atypes']
    mol_tip.apos   = tip_scf['apos']
    mol_smp = AtomicSystem()
    mol_smp.atypes = smp_scf['atypes']
    mol_smp.apos   = smp_scf['apos']

    homo_tip = _homo_index(tip_scf['eigen'])
    homo_smp = _homo_index(smp_scf['eigen'])

    tip_mo_list = _parse_mo_list(args.tip_mo_list)
    smp_mo_list = _parse_mo_list(args.smp_mo_list)
    if not tip_mo_list:
        tip_mo_list = [args.tip_mo-1] if (args.tip_mo is not None) else [max(0, homo_tip-1)]
    if not smp_mo_list:
        smp_mo_list = [args.smp_mo-1] if (args.smp_mo is not None) else [max(0, homo_smp-1)]

    tip_mo_list = [int(np.clip(m, 0, tip_scf['norb']-1)) for m in tip_mo_list]
    smp_mo_list = [int(np.clip(m, 0, smp_scf['norb']-1)) for m in smp_mo_list]

    tip_ang_list = [float(s) for s in args.tip_angles.split(',') if s.strip()] if args.tip_angles.strip() else [float(args.tip_angle)]
    smp_ang_list = [float(s) for s in args.smp_angles.split(',') if s.strip()] if args.smp_angles.strip() else [float(args.smp_angle)]

    zmid = (0.5*args.ztip) if (args.zmid is None) else float(args.zmid)

    # Center sample at z=0 (only z shift), keep xy as-is
    cen_smp = np.asarray(mol_smp.apos.mean(axis=0), dtype=np.float64)

    xs = np.linspace(cen_smp[0] - args.size*0.5, cen_smp[0] + args.size*0.5, args.n)
    ys = np.linspace(cen_smp[1] - args.size*0.5, cen_smp[1] + args.size*0.5, args.n)
    X2, Y2 = np.meshgrid(xs, ys, indexing='ij')
    pts_mid = np.zeros((args.n*args.n, 3), dtype=np.float64)
    pts_mid[:, 0] = X2.ravel(); pts_mid[:, 1] = Y2.ravel(); pts_mid[:, 2] = zmid
    extent2 = [xs[0], xs[-1], ys[0], ys[-1]]

    dx = (xs[1] - xs[0])
    dy = (ys[1] - ys[0])
    sx = (np.arange(args.n) - args.n//2) * dx
    sy = (np.arange(args.n) - args.n//2) * dy
    extent_shift = [sx[0], sx[-1], sy[0], sy[-1]]

    # Tip centers follow the same xy scan, fixed z=ztip
    tip_centers = pts_mid.copy().astype(np.float32)
    tip_centers[:, 2] = float(args.ztip)

    # Tip local geometry in unrotated coordinates (relative to tip center)
    tip_pos_rel0 = np.asarray(mol_tip.apos, dtype=np.float32)
    tip_pos_rel0 = tip_pos_rel0 - tip_pos_rel0.mean(axis=0)

    gp = GridProjector(fdata_dir=fdata_dir)

    for tip_ang in tip_ang_list:
        R_tip = rotation_matrix(args.tip_axis, tip_ang)
        q_tip = quat_from_axis_angle(args.tip_axis, tip_ang)
        tip_quat = np.tile(q_tip[None, :], (tip_centers.shape[0], 1)).astype(np.float32)

        # For plotting only: build rotated absolute positions of tip (shifted to z=ztip)
        tip_pos_rot = rotate_points(mol_tip.apos, R_tip)
        tip_pos_rot = np.asarray(tip_pos_rot, dtype=np.float64)
        tip_pos_rot = tip_pos_rot - tip_pos_rot.mean(axis=0)
        tip_pos_rot = tip_pos_rot + np.array([cen_smp[0], cen_smp[1], float(args.ztip)], dtype=np.float64)

        for smp_ang in smp_ang_list:
            R_smp = rotation_matrix(args.smp_axis, smp_ang)
            smp_pos_rot = rotate_points(mol_smp.apos, R_smp)
            smp_pos = np.asarray(smp_pos_rot, dtype=np.float64) + np.array([0.0, 0.0, -cen_smp[2]], dtype=np.float64)

            # For visualization of MOs at mid-plane
            C_tip0 = tip_scf['C']
            C_smp0 = smp_scf['C']
            C_tip_rot = rotate_mo_coefficients(C_tip0.T, tip_scf['norb_per'], R_tip).T
            C_smp_rot = rotate_mo_coefficients(C_smp0.T, smp_scf['norb_per'], R_smp).T

            for tip_mo in tip_mo_list:
                for smp_mo in smp_mo_list:
                    psi_s = project_orbital_to_points(C_smp_rot, smp_mo, mol_smp.atypes, smp_pos, smp_scf['orb2atom'], smp_scf['norb_per'], fdata_basis_dir, pts_mid.astype(np.float32))
                    psi_t = project_orbital_to_points(C_tip_rot, tip_mo, mol_tip.atypes, tip_pos_rot, tip_scf['orb2atom'], tip_scf['norb_per'], fdata_basis_dir, pts_mid.astype(np.float32))
                    psi_s = np.asarray(psi_s, dtype=np.float64).reshape(args.n, args.n)
                    psi_t = np.asarray(psi_t, dtype=np.float64).reshape(args.n, args.n)

                    Fs = np.fft.fft2(psi_s)
                    Ft = np.fft.fft2(psi_t)
                    corr = np.fft.ifft2(Fs * np.conj(Ft)).real
                    corr = np.fft.fftshift(corr)
                    corr2 = corr * corr

                    coeff_tip0 = remap_coeffs_fortran_to_grid(C_tip0[tip_mo, :], tip_scf['norb_per'])
                    coeff_smp0 = remap_coeffs_fortran_to_grid(C_smp0[smp_mo, :], smp_scf['norb_per'])

                    t_amp, t_I = gp.mo_overlap_points_exp_sk_2mol(
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
                    fig.suptitle(f"2mol overlap-STM vs xcorr: tip {tip_name} MO{tip_mo+1} smp {smp_name} MO{smp_mo+1} zmid={zmid:.2f} ztip={args.ztip:.2f} R_tip={args.tip_axis}{int(tip_ang)} R_smp={args.smp_axis}{int(smp_ang)}")

                    vmax_s = float(np.max(np.abs(psi_s)))
                    if vmax_s < 1e-30: vmax_s = 1.0
                    im = ax[0,0].imshow(psi_s.T, origin='lower', extent=extent2, cmap='bwr', vmin=-vmax_s, vmax=vmax_s, aspect='equal')
                    ax[0,0].set_title('1) ψ_sample(mid)'); ax[0,0].set_xlabel('x (Å)'); ax[0,0].set_ylabel('y (Å)')
                    plot_atoms(ax[0,0], smp_pos, mol_smp.atypes)
                    fig.colorbar(im, ax=ax[0,0], fraction=0.046, pad=0.04)

                    vmax_t = float(np.max(np.abs(psi_t)))
                    if vmax_t < 1e-30: vmax_t = 1.0
                    im = ax[0,1].imshow(psi_t.T, origin='lower', extent=extent2, cmap='bwr', vmin=-vmax_t, vmax=vmax_t, aspect='equal')
                    ax[0,1].set_title('2) ψ_tip(mid)'); ax[0,1].set_xlabel('x (Å)'); ax[0,1].set_ylabel('y (Å)')
                    plot_atoms(ax[0,1], tip_pos_rot, mol_tip.atypes)
                    fig.colorbar(im, ax=ax[0,1], fraction=0.046, pad=0.04)

                    vI = float(np.max(t_I))
                    if vI <= 0: vI = 1.0
                    im = ax[1,0].imshow(t_I.T, origin='lower', extent=extent2, cmap='viridis', vmin=0.0, vmax=vI, aspect='equal')
                    ax[1,0].set_title('3) STM overlap |t|^2 (2mol)'); ax[1,0].set_xlabel('x (Å)'); ax[1,0].set_ylabel('y (Å)')
                    plot_atoms(ax[1,0], smp_pos, mol_smp.atypes)
                    fig.colorbar(im, ax=ax[1,0], fraction=0.046, pad=0.04)

                    vC = float(np.max(corr2))
                    if vC <= 0: vC = 1.0
                    im = ax[1,1].imshow(corr2.T, origin='lower', extent=extent_shift, cmap='viridis', vmin=0.0, vmax=vC, aspect='equal')
                    ax[1,1].set_title('4) xcorr(ψ_sample,ψ_tip)^2'); ax[1,1].set_xlabel('dx (Å)'); ax[1,1].set_ylabel('dy (Å)')
                    fig.colorbar(im, ax=ax[1,1], fraction=0.046, pad=0.04)

                    save_plot(fig, f"overlap2mol_tip{tip_name}_MO{tip_mo+1:03d}_smp{smp_name}_MO{smp_mo+1:03d}_zmid{zmid:.2f}_ztip{args.ztip:.2f}_tip{args.tip_axis}{int(tip_ang):03d}_smp{args.smp_axis}{int(smp_ang):03d}", export_dir=export_dir)
                    plt.close(fig)


if __name__ == "__main__":
    main()
