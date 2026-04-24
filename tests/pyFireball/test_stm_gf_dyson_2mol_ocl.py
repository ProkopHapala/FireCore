#!/usr/bin/env python3

import os, sys, argparse, subprocess
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import time

sys.path.append(os.path.join(os.path.dirname(__file__), "..", ".."))

from pyBall.AtomicSystem import AtomicSystem
from pyBall import FireCore as fc

from pyBall.FireballOCL.STM_utils import (
    set_export_dir, save_plot, plot_atoms,
    project_orbital_to_points,
    get_orbital_layout,
    remap_coeffs_fortran_to_grid,
)

from pyBall.FireballOCL.STM import (
    rotate_mo_coefficients, rotate_points,
)

from pyBall.FireballOCL.Grid import GridProjector

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
    Hk, Sk = fc.get_HS_k((0.0, 0.0, 0.0), norb)
    norb_per, starts, orb2atom = get_orbital_layout_from_fireball(dims)
    return {
        'dims': dims,
        'norb': norb,
        'eigen': np.array(eigen, dtype=np.float64),
        'C': np.array(C_fc, dtype=np.float64),
        'H': np.array(Hk, dtype=np.complex128),
        'S': np.array(Sk, dtype=np.complex128),
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
        H=np.asarray(scf['H'], dtype=np.complex128),
        S=np.asarray(scf['S'], dtype=np.complex128),
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
        'H': np.asarray(d['H'], dtype=np.complex128),
        'S': np.asarray(d['S'], dtype=np.complex128),
        'norb_per': np.asarray(d['norb_per'], dtype=np.int32),
        'starts': np.asarray(d['starts'], dtype=np.int32),
        'orb2atom': np.asarray(d['orb2atom'], dtype=np.int32),
        'atypes': np.asarray(d['atypes'], dtype=np.int32),
        'apos': np.asarray(d['apos'], dtype=np.float64),
    }
    return out


def get_green_from_cache(scf, E, eta):
    z = float(E) + 1j*float(eta)
    H = np.asarray(scf['H'], dtype=np.complex128)
    S = np.asarray(scf['S'], dtype=np.complex128)
    A = z*S - H
    return np.linalg.inv(A)


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


def main():
    ap = argparse.ArgumentParser(description='Two-molecule STM GF(Dyson) GPU kernel test (work-item per pixel)')
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

    ap.add_argument('--nmax_scf', type=int, default=50)
    ap.add_argument('--outdir', default='export/stm_gf_dyson_2mol_ocl')

    ap.add_argument('--n', type=int, default=80)
    ap.add_argument('--size', type=float, default=20.0)
    ap.add_argument('--ztip', type=float, default=3.0)
    ap.add_argument('--zmid', type=float, default=None)
    ap.add_argument('--rcut', type=float, default=8.0)
    ap.add_argument('--beta', type=float, default=1.0)
    ap.add_argument('--r0', type=float, default=3.0)
    ap.add_argument('--eta', type=float, default=1e-2)

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

    mol_tip = AtomicSystem(); mol_tip.atypes = tip_scf['atypes']; mol_tip.apos = tip_scf['apos']
    mol_smp = AtomicSystem(); mol_smp.atypes = smp_scf['atypes']; mol_smp.apos = smp_scf['apos']

    fc.setVerbosity(0)
    init_atypes = np.concatenate([np.asarray(mol_tip.atypes, dtype=np.int32), np.asarray(mol_smp.atypes, dtype=np.int32)])
    init_apos = np.vstack([
        np.asarray(mol_tip.apos, dtype=np.float64),
        np.asarray(mol_smp.apos, dtype=np.float64) + np.array([0.0, 0.0, 20.0], dtype=np.float64),
    ])
    fc.initialize(atomType=init_atypes, atomPos=init_apos, verbosity=0)
    species_map = {}
    next_ispec = 1
    for z in init_atypes:
        z = int(z)
        if z not in species_map:
            species_map[z] = next_ispec
            next_ispec += 1
    tip_species = np.array([species_map[int(z)] for z in mol_tip.atypes], dtype=np.int32)
    smp_species = np.array([species_map[int(z)] for z in mol_smp.atypes], dtype=np.int32)

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

    cen_smp = np.asarray(mol_smp.apos.mean(axis=0), dtype=np.float64)

    xs = np.linspace(cen_smp[0] - args.size*0.5, cen_smp[0] + args.size*0.5, args.n)
    ys = np.linspace(cen_smp[1] - args.size*0.5, cen_smp[1] + args.size*0.5, args.n)
    X2, Y2 = np.meshgrid(xs, ys, indexing='ij')
    pts_mid = np.zeros((args.n*args.n, 3), dtype=np.float64)
    pts_mid[:, 0] = X2.ravel(); pts_mid[:, 1] = Y2.ravel(); pts_mid[:, 2] = zmid
    extent2 = [xs[0], xs[-1], ys[0], ys[-1]]

    tip_centers = pts_mid.copy().astype(np.float32)
    tip_centers[:, 2] = float(args.ztip)

    tip_pos_rel0 = np.asarray(mol_tip.apos, dtype=np.float32)
    tip_pos_rel0 = tip_pos_rel0 - tip_pos_rel0.mean(axis=0)

    gp = GridProjector(fdata_dir=fdata_dir)

    for tip_ang in tip_ang_list:
        R_tip = rotation_matrix(args.tip_axis, tip_ang)

        tip_pos_rot = rotate_points(mol_tip.apos, R_tip)
        tip_pos_rot = np.asarray(tip_pos_rot, dtype=np.float64)
        tip_pos_rot = tip_pos_rot - tip_pos_rot.mean(axis=0)
        tip_pos_rot = tip_pos_rot + np.array([cen_smp[0], cen_smp[1], float(args.ztip)], dtype=np.float64)

        for smp_ang in smp_ang_list:
            R_smp = rotation_matrix(args.smp_axis, smp_ang)
            smp_pos_rot = rotate_points(mol_smp.apos, R_smp)
            smp_pos = np.asarray(smp_pos_rot, dtype=np.float64) + np.array([0.0, 0.0, -cen_smp[2]], dtype=np.float64)

            C_tip0 = tip_scf['C']
            C_smp0 = smp_scf['C']
            C_tip_rot = rotate_mo_coefficients(C_tip0.T, tip_scf['norb_per'], R_tip).T
            C_smp_rot = rotate_mo_coefficients(C_smp0.T, smp_scf['norb_per'], R_smp).T

            for tip_mo in tip_mo_list:
                for smp_mo in smp_mo_list:
                    E_tip = float(tip_scf['eigen'][tip_mo])
                    E_smp = float(smp_scf['eigen'][smp_mo])
                    E_transport = 0.5*(E_tip + E_smp)
                    G_tip = get_green_from_cache(tip_scf, E=E_transport, eta=args.eta)
                    G_smp = get_green_from_cache(smp_scf, E=E_transport, eta=args.eta)

                    nt = int(G_tip.shape[0])
                    ns = int(G_smp.shape[0])
                    nt_expected = int(len(tip_scf['orb2atom']))
                    ns_expected = int(len(smp_scf['orb2atom']))
                    if nt != nt_expected:
                        raise ValueError(f"GT size mismatch: got {nt}, expected {nt_expected}")
                    if ns != ns_expected:
                        raise ValueError(f"GS size mismatch: got {ns}, expected {ns_expected}")
                    c_tip = np.asarray(C_tip_rot[tip_mo, :], dtype=np.complex128)
                    c_smp = np.asarray(C_smp_rot[smp_mo, :], dtype=np.complex128)

                    psi_s = project_orbital_to_points(C_smp_rot, smp_mo, mol_smp.atypes, smp_pos, smp_scf['orb2atom'], smp_scf['norb_per'], fdata_basis_dir, pts_mid.astype(np.float32))
                    psi_t = project_orbital_to_points(C_tip_rot, tip_mo, mol_tip.atypes, tip_pos_rot, tip_scf['orb2atom'], tip_scf['norb_per'], fdata_basis_dir, pts_mid.astype(np.float32))
                    psi_s = np.asarray(psi_s, dtype=np.float64).reshape(args.n, args.n)
                    psi_t = np.asarray(psi_t, dtype=np.float64).reshape(args.n, args.n)

                    # --- Fortran reference (full GF chain) ---
                    print(f"Fortran GF tipMO={tip_mo+1} smpMO={smp_mo+1} E_tip={E_tip:+.4f} E_smp={E_smp:+.4f} E={E_transport:+.4f}")
                    t0 = time.time()
                    gf_I_fort = fc.stm_gf_2mol_mo_2points(
                        points=np.asarray(tip_centers, dtype=np.float64),
                        tip_pos_rel=np.asarray(tip_pos_rel0, dtype=np.float64),
                        tip_atypes=tip_species,
                        tip_orb2atom=np.asarray(tip_scf['orb2atom'], dtype=np.int32),
                        smp_pos=np.asarray(smp_pos, dtype=np.float64),
                        smp_atypes=smp_species,
                        smp_orb2atom=np.asarray(smp_scf['orb2atom'], dtype=np.int32),
                        GT_global=G_tip,
                        GS_global=G_smp,
                        c_tip=c_tip,
                        c_smp=c_smp,
                        beta=float(args.beta), r0=float(args.r0), rcut=float(args.rcut), overlap_scale=0.0,
                    )
                    t_fort = time.time() - t0
                    print(f"  Fortran time = {t_fort:.3f} s")
                    gf_I_fort = np.asarray(gf_I_fort, dtype=np.float64).reshape(args.n, args.n)

                    # --- GPU Green's Function (simplified exponential SK) ---
                    print(f"GPU GF tipMO={tip_mo+1} smpMO={smp_mo+1}")
                    t0 = time.time()
                    gf_I_gpu = gp.stm_gf_dyson_2mol_mo_scan(
                        tip_centers=tip_centers,
                        tip_pos_rel=tip_pos_rel0,
                        smp_pos=smp_pos.astype(np.float32),
                        GT_global=G_tip,
                        GS_global=G_smp,
                        c_tip=c_tip,
                        c_smp=c_smp,
                        tip_norb_per=tip_scf['norb_per'],
                        smp_norb_per=smp_scf['norb_per'],
                        beta=float(args.beta),
                        r0=float(args.r0),
                        rcut=float(args.rcut),
                    )
                    t_gpu = time.time() - t0
                    print(f"  GPU time = {t_gpu*1e3:.2f} ms  (×{t_fort/max(t_gpu,1e-9):.0f} faster)")
                    gf_I_gpu = np.asarray(gf_I_gpu, dtype=np.float64).reshape(args.n, args.n)

                    # --- Overlap-STM (existing code path) for comparison ---
                    coeff_tip0 = remap_coeffs_fortran_to_grid(C_tip0[tip_mo, :], tip_scf['norb_per'])
                    coeff_smp0 = remap_coeffs_fortran_to_grid(C_smp0[smp_mo, :], smp_scf['norb_per'])
                    _, ov_I = gp.mo_overlap_points_exp_sk_2mol(
                        tip_centers=tip_centers,
                        tip_pos_rel=tip_pos_rel0,
                        smp_pos=smp_pos.astype(np.float32),
                        coeffs_tip=coeff_tip0.astype(np.float32),
                        coeffs_smp=coeff_smp0.astype(np.float32),
                        beta=float(args.beta), r0=float(args.r0), rcut=float(args.rcut),
                    )
                    ov_I = np.asarray(ov_I, dtype=np.float64).reshape(args.n, args.n)

                    # Compute statistics for scatter panel
                    gf_flat_fort = gf_I_fort.ravel()
                    gf_flat_gpu  = gf_I_gpu.ravel()
                    try:
                        corr = np.corrcoef(gf_flat_fort, gf_flat_gpu)[0, 1]
                    except Exception:
                        corr = np.nan
                    try:
                        # log-space correlation — more sensitive to pattern than scale
                        gf_log_fort = np.log10(gf_flat_fort + 1e-30)
                        gf_log_gpu  = np.log10(gf_flat_gpu  + 1e-30)
                        corr_log = np.corrcoef(gf_log_fort, gf_log_gpu)[0, 1]
                    except Exception:
                        corr_log = np.nan

                    # Compute overlap-STM correlation with GF results
                    ov_flat = ov_I.ravel()
                    try:
                        corr_ov_fort = np.corrcoef(ov_flat, gf_flat_fort)[0, 1]
                        corr_ov_gpu  = np.corrcoef(ov_flat, gf_flat_gpu)[0, 1]
                    except Exception:
                        corr_ov_fort = corr_ov_gpu = np.nan

                    speedup = t_fort / max(t_gpu, 1e-9)

                    # === Plot: 2 rows × 3 columns (6 panels, clean linear scales) ===
                    fig, ax = plt.subplots(2, 3, figsize=(20, 13))

                    suptitle = (
                        f"STM GF 2mol: tip={tip_name} MO{tip_mo+1}  smp={smp_name} MO{smp_mo+1}  "
                        f"ztip={args.ztip:.1f}A  zmid={zmid:.1f}A  beta={args.beta}  r0={args.r0}"
                    )
                    fig.suptitle(suptitle, fontsize=11, fontweight='bold')

                    # Row 1: orbitals and overlap STM (each with own symmetric scale)
                    # 1) psi_sample at mid-plane
                    vmax_s = float(np.max(np.abs(psi_s)))
                    if vmax_s < 1e-30: vmax_s = 1.0
                    im = ax[0, 0].imshow(psi_s.T, origin='lower', extent=extent2,
                                           cmap='bwr', vmin=-vmax_s, vmax=vmax_s, aspect='equal')
                    ax[0, 0].set_title(f'sample MO{smp_mo+1}  $\\psi_s$  at z={zmid:.1f}A', fontsize=10)
                    plot_atoms(ax[0, 0], smp_pos, mol_smp.atypes)
                    fig.colorbar(im, ax=ax[0, 0], fraction=0.046, pad=0.04)

                    # 2) psi_tip at mid-plane
                    vmax_t = float(np.max(np.abs(psi_t)))
                    if vmax_t < 1e-30: vmax_t = 1.0
                    im = ax[0, 1].imshow(psi_t.T, origin='lower', extent=extent2,
                                           cmap='bwr', vmin=-vmax_t, vmax=vmax_t, aspect='equal')
                    ax[0, 1].set_title(f'tip MO{tip_mo+1}  $\\psi_t$  at z={zmid:.1f}A', fontsize=10)
                    plot_atoms(ax[0, 1], tip_pos_rot, mol_tip.atypes)
                    fig.colorbar(im, ax=ax[0, 1], fraction=0.046, pad=0.04)

                    # 3) Overlap STM (intensity, own scale)
                    vmax_ov = float(np.max(ov_I)) if np.max(ov_I) > 0 else 1.0
                    im = ax[0, 2].imshow(ov_I.T, origin='lower', extent=extent2,
                                           cmap='viridis', vmin=0, vmax=vmax_ov, aspect='equal')
                    ax[0, 2].set_title(f'overlap-STM  $|t|^2$  max={vmax_ov:.2e}', fontsize=10)
                    plot_atoms(ax[0, 2], smp_pos, mol_smp.atypes)
                    fig.colorbar(im, ax=ax[0, 2], fraction=0.046, pad=0.04)

                    # Row 2: Fortran GF, GPU GF, scatter plot (each with OWN scale)
                    # 4) Fortran GF
                    vmax_f = float(np.max(gf_I_fort)) if np.max(gf_I_fort) > 0 else 1.0
                    im = ax[1, 0].imshow(gf_I_fort.T, origin='lower', extent=extent2,
                                           cmap='viridis', vmin=0, vmax=vmax_f, aspect='equal')
                    ax[1, 0].set_title(
                        f'Fortran GF  $|G_t M_{{ts}} G_s|^2$  '
                        f'max={vmax_f:.2e}  {t_fort:.1f}s',
                        fontsize=10)
                    plot_atoms(ax[1, 0], smp_pos, mol_smp.atypes)
                    fig.colorbar(im, ax=ax[1, 0], fraction=0.046, pad=0.04)

                    # 5) GPU GF
                    vmax_g = float(np.max(gf_I_gpu)) if np.max(gf_I_gpu) > 0 else 1.0
                    im = ax[1, 1].imshow(gf_I_gpu.T, origin='lower', extent=extent2,
                                           cmap='viridis', vmin=0, vmax=vmax_g, aspect='equal')
                    ax[1, 1].set_title(
                        f'GPU GF (OCL)  $|G_t M_{{ts}} G_s|^2$  '
                        f'max={vmax_g:.2e}  {t_gpu*1e3:.2f}ms  ×{speedup:.0f}',
                        fontsize=10)
                    plot_atoms(ax[1, 1], smp_pos, mol_smp.atypes)
                    fig.colorbar(im, ax=ax[1, 1], fraction=0.046, pad=0.04)

                    # 6) Scatter: GPU GF vs Fortran GF (pixel-by-pixel linear correlation)
                    mask = (gf_flat_fort > 0) & (gf_flat_gpu > 0)
                    xv = gf_flat_fort[mask]; yv = gf_flat_gpu[mask]
                    if len(xv) > 2:
                        ax[1, 2].loglog(xv, yv, '.', ms=2, alpha=0.5, color='navy')
                        # diagonal reference line
                        allmax = max(xv.max(), yv.max(), 1e-30)
                        allmin = max(xv.min(), yv.min(), 1e-30)
                        ax[1, 2].loglog([allmin*0.5, allmax*2], [allmin*0.5, allmax*2],
                                         'r--', lw=1.0, alpha=0.6, label='y=x')
                        # linear fit in log-log space
                        lx = np.log10(xv + 1e-30); ly = np.log10(yv + 1e-30)
                        slope, intercept = np.polyfit(lx, ly, 1)
                        fit_x = np.array([lx.min(), lx.max()])
                        ax[1, 2].plot(10**fit_x, 10**(slope*fit_x + intercept),
                                       'g-', lw=1.5, label=f'slope={slope:.2f}')
                        ax[1, 2].legend(fontsize=7, loc='upper left')
                        ax[1, 2].set_xlabel('Fortran GF intensity'); ax[1, 2].set_ylabel('GPU GF intensity')
                        ax[1, 2].set_title(
                            f'pixel scatter  corr={corr:.4f}  log-corr={corr_log:.4f}  '
                            f'n={len(xv)}\n'
                            f'ov<->F={corr_ov_fort:.3f}  ov<->G={corr_ov_gpu:.3f}',
                            fontsize=9)
                        ax[1, 2].grid(True, alpha=0.3)
                    else:
                        ax[1, 2].text(0.5, 0.5, 'no data', ha='center', va='center',
                                       transform=ax[1, 2].transAxes)
                        ax[1, 2].set_title('scatter (no signal)')

                    save_plot(fig,
                        f"ocl_gf_dyson2mol_tip{tip_name}_MO{tip_mo+1:03d}_smp{smp_name}_MO{smp_mo+1:03d}"
                        f"_zmid{zmid:.2f}_ztip{args.ztip:.2f}"
                        f"_tip{args.tip_axis}{int(tip_ang):03d}_smp{args.smp_axis}{int(smp_ang):03d}",
                        export_dir=export_dir)
                    plt.close(fig)


if __name__ == "__main__":
    main()
