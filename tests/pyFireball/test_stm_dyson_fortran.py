#!/usr/bin/env python3

import os, sys, argparse, time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.path.dirname(__file__), "..", ".."))

from pyBall.AtomicSystem import AtomicSystem
from pyBall import FireCore as fc
from pyBall.FireballOCL.STM_utils import set_export_dir, save_plot, plot_atoms, build_plane_grid, project_orbital_to_points, get_orbital_layout

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
    return homo, min(homo+1, len(eigen)-1)


def main():
    ap = argparse.ArgumentParser(description='Fortran-native Dyson STM image scan')
    ap.add_argument('--xyz', default='../../cpp/common_resources/xyz/PTCDA.xyz')
    ap.add_argument('--mo', type=int, default=None, help='Sample MO index 1-based; default HOMO')
    ap.add_argument('--nmax_scf', type=int, default=100)
    ap.add_argument('--n', type=int, default=80)
    ap.add_argument('--size', type=float, default=20.0)
    ap.add_argument('--z', type=float, default=3.0)
    ap.add_argument('--eta', type=float, default=1e-2)
    ap.add_argument('--tipZ', type=int, default=1)
    ap.add_argument('--mode', type=int, default=0)
    ap.add_argument('--rcut', type=float, default=8.0)
    ap.add_argument('--beta', type=float, default=1.0)
    ap.add_argument('--r0', type=float, default=3.0)
    ap.add_argument('--A_ss', type=float, default=-1.0)
    ap.add_argument('--A_sp', type=float, default=-1.0)
    ap.add_argument('--A_pp_sig', type=float, default=-1.0)
    ap.add_argument('--A_pp_pi', type=float, default=1.0)
    ap.add_argument('--A_ps', type=float, default=None)
    ap.add_argument('--overlap_scale', type=float, default=0.0)
    ap.add_argument('--E_tip', type=float, default=0.0)
    ap.add_argument('--outdir', default='export/stm_dyson_fortran')
    args = ap.parse_args()

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
    print(f'DEBUG SCF start mol={mol_name} natoms={natoms} norb={norb}')
    fc.evalForce(mol.apos, nmax_scf=args.nmax_scf)
    eigen = fc.get_eigen(ikp=1, norb=norb)
    C = fc.get_wfcoef(norb=norb)

    homo, lumo = _homo_lumo(eigen)
    mo0 = (args.mo-1) if (args.mo is not None) else homo
    E = float(eigen[mo0])
    label = 'HOMO' if mo0 == homo else ('LUMO' if mo0 == lumo else f'MO{mo0+1}')
    print(f'DEBUG using MO={mo0+1} label={label} E={E:+.8f}')

    data = fc.get_HS_neighs(dims)
    data = fc.get_HS_sparse(dims, data)
    norb_per, starts = get_orbital_layout(data, natoms)
    orb2atom = np.array([ia for ia in range(natoms) for _ in range(int(norb_per[ia]))], dtype=np.int32)

    X, Y, pts, extent = build_plane_grid(mol.apos, z=args.z, size=args.size, n=args.n)
    fdata_basis_dir = os.path.join(_THIS_DIR, 'Fdata', 'basis')
    psi = project_orbital_to_points(np.asarray(C, dtype=np.float64), mo0, mol.atypes, mol.apos, orb2atom, norb_per, fdata_basis_dir, pts.astype(np.float32))
    psi = np.asarray(psi, dtype=np.float64).reshape(args.n, args.n)

    A_ps = args.A_sp if (args.A_ps is None) else args.A_ps
    print('DEBUG running Fortran Dyson point scan')
    t0 = time.time()
    I = fc.stm_dyson_pointscan(
        points=np.asarray(pts, dtype=np.float64),
        E=E,
        eta=args.eta,
        mode=args.mode,
        ikpoint=1,
        tipZ=args.tipZ,
        rcut=args.rcut,
        beta=args.beta,
        r0=args.r0,
        A_ss=args.A_ss,
        A_sp=args.A_sp,
        A_pp_sig=args.A_pp_sig,
        A_pp_pi=args.A_pp_pi,
        A_ps=A_ps,
        overlap_scale=args.overlap_scale,
        E_tip=args.E_tip,
    )
    dt = time.time() - t0
    I = np.asarray(I, dtype=np.float64).reshape(args.n, args.n)
    Ilog = np.log10(I + 1e-30)
    print(f'DEBUG Fortran Dyson scan done in {dt:.3f} s  Imax={np.max(I):.6e}')

    fig, ax = plt.subplots(1, 3, figsize=(18, 6))
    fig.suptitle(f'Fortran Dyson STM {mol_name} MO {mo0+1} ({label}) E={E:+.4f} z={args.z:.2f} eta={args.eta:.2e}')

    vmax = float(np.max(np.abs(psi)))
    if vmax < 1e-30:
        vmax = 1.0
    im = ax[0].imshow(psi.T, origin='lower', extent=extent, cmap='bwr', vmin=-vmax, vmax=vmax, aspect='equal')
    ax[0].set_title('1) sample MO')
    ax[0].set_xlabel('x (Å)')
    ax[0].set_ylabel('y (Å)')
    plot_atoms(ax[0], mol.apos, mol.atypes)
    fig.colorbar(im, ax=ax[0], fraction=0.046, pad=0.04)

    im = ax[1].imshow(Ilog.T, origin='lower', extent=extent, cmap='viridis', aspect='equal')
    ax[1].set_title('2) Fortran Dyson log10(current)')
    ax[1].set_xlabel('x (Å)')
    ax[1].set_ylabel('y (Å)')
    plot_atoms(ax[1], mol.apos, mol.atypes)
    fig.colorbar(im, ax=ax[1], fraction=0.046, pad=0.04)

    im = ax[2].imshow((Ilog - np.log10(psi*psi + 1e-30)).T, origin='lower', extent=extent, cmap='bwr', aspect='equal')
    ax[2].set_title('3) log10(Dyson) - log10(psi^2)')
    ax[2].set_xlabel('x (Å)')
    ax[2].set_ylabel('y (Å)')
    plot_atoms(ax[2], mol.apos, mol.atypes)
    fig.colorbar(im, ax=ax[2], fraction=0.046, pad=0.04)

    plt.tight_layout()
    save_plot(fig, f'fortran_dyson_mo{mo0+1:04d}_{label.lower()}', export_dir)
    plt.close(fig)
    print(f'DEBUG output in {export_dir}')


if __name__ == '__main__':
    main()
