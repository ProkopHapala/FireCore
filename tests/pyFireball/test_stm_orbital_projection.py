#!/usr/bin/env python3
"""STM Orbital-Projection Debug Test (PTCDA)

For each MO around HOMO/LUMO (+/-4), and for z in {2.0,10.0} Å, make a 1x3 panel:
1) Orbital projection ψ(r) on xy grid (OpenCL exact orbital projection kernel)
2) STM overlap amplitude t(r) using TRUE Fireball overlap S_TS from Fdata (finite cutoff)
3) STM overlap amplitude t(r) using EXP-extended overlap via exp(-beta*(r-r0)) * SK angular (long-range)

Tip: single H atom with only s orbital.

Run from:
  cd tests/pyFireball
  python test_stm_orbital_projection.py

Outputs:
  export/stm_orbital_projection/ptcda_z{z}/moXXXX_LABEL_panels.png
"""

import os
import sys
import warnings
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.path.dirname(__file__), "..", ".."))

from pyBall.AtomicSystem import AtomicSystem
from pyBall import FireCore as fc
from pyBall.FireballOCL.STM_utils import (
    set_export_dir, save_plot, plot_atoms, project_orbital_to_points, project_orbital_to_points_exp,
    get_orbital_layout, sparse_to_dense
)
from pyBall.FireballOCL.STM import (
    build_inter_system_blocks_fdata, build_inter_system_blocks_exp_sk,
    _atom_orb_starts, response_amplitude_map
)
from pyBall.FireballOCL.Grid import GridProjector


def _get_orbital_mapping_from_fireball(dims):
    data = fc.get_HS_neighs(dims)
    data = fc.get_HS_sparse(dims, data)
    norb_per, starts = get_orbital_layout(data, dims.natoms)
    if int(starts[-1]) != int(dims.norbitals):
        raise RuntimeError(f"Orbital count mismatch: mapped {int(starts[-1])} vs dims.norbitals={int(dims.norbitals)}")
    orb2atom = np.array([ia for ia in range(dims.natoms) for _ in range(int(norb_per[ia]))], dtype=np.int32)
    return orb2atom, norb_per


def _homo_lumo_from_eigen(eigen):
    occ = np.where(np.asarray(eigen) < 0.0)[0]
    if len(occ) == 0:
        homo = len(eigen)//2 - 1
    else:
        homo = int(occ[-1])
    lumo = homo + 1
    return homo, lumo


def _mo_label(i, homo, lumo):
    if i == homo: return "HOMO"
    if i == lumo: return "LUMO"
    if i < homo:  return f"HOMO{i-homo:+d}"  # negative
    return f"LUMO{i-lumo:+d}"                # positive


def _build_xy_grid(mol_pos, z, size=20.0, n=80):
    origin = mol_pos.mean(axis=0)
    step = float(size) / int(n)
    grid_origin = origin - np.array([size/2, size/2, 0.0])
    xs = np.linspace(grid_origin[0] + step*0.5, grid_origin[0] + size - step*0.5, int(n))
    ys = np.linspace(grid_origin[1] + step*0.5, grid_origin[1] + size - step*0.5, int(n))
    X, Y = np.meshgrid(xs, ys, indexing='ij')
    Z = np.zeros_like(X) + (origin[2] + float(z))
    points = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)
    points = np.ascontiguousarray(points, dtype=np.float32)
    extent_abs = [grid_origin[0], grid_origin[0] + size, grid_origin[1], grid_origin[1] + size]
    return X, Y, points, extent_abs


def main():
    os.chdir(os.path.dirname(__file__))

    # ensure Fdata symlink exists (same convention as other scripts)
    _THIS_DIR  = os.path.dirname(os.path.abspath(__file__))
    _REPO_ROOT = os.path.normpath(os.path.join(_THIS_DIR, "..", ".."))
    _FDATA_HCNOS = os.path.join(_REPO_ROOT, "tests", "Fireball", "Fdata_HCNOS")
    _FDATA_LOCAL = os.path.join(_THIS_DIR, "Fdata")
    if os.path.realpath(_FDATA_LOCAL) if os.path.exists(_FDATA_LOCAL) else "" != os.path.realpath(_FDATA_HCNOS):
        if os.path.lexists(_FDATA_LOCAL):
            os.unlink(_FDATA_LOCAL)
        os.symlink(_FDATA_HCNOS, _FDATA_LOCAL)

    xyz = "../../cpp/common_resources/xyz/PTCDA.xyz"
    mol = AtomicSystem(fname=xyz)

    fc.setVerbosity(0)
    fc.initialize(atomType=mol.atypes, atomPos=mol.apos, verbosity=0)
    dims = fc.get_HS_dims()
    norb = int(dims.norbitals)

    print(f"Running SCF for PTCDA (norb={norb})...")
    fc.evalForce(mol.apos, nmax_scf=200)
    eigen = fc.get_eigen(ikp=1, norb=norb)

    homo, lumo = _homo_lumo_from_eigen(eigen)
    print(f"HOMO={homo+1} E={eigen[homo]:.4f} eV")
    print(f"LUMO={lumo+1} E={eigen[lumo]:.4f} eV")

    mo_list = list(range(homo-4, homo+5))
    mo_list = [i for i in mo_list if 0 <= i < norb]

    # coefficients matrix for orbital projection and overlap
    C_fc = fc.get_wfcoef(norb=norb)

    # authoritative mapping
    orb2atom, norb_per_smp = _get_orbital_mapping_from_fireball(dims)

    # extract sample H, S as dense matrices (once)
    natoms = dims.natoms
    data_hs = fc.get_HS_neighs(dims)
    data_hs = fc.get_HS_sparse(dims, data_hs)
    H_s = _blocked_to_dense(data_hs, data_hs.h_mat, natoms)
    S_s = _blocked_to_dense(data_hs, data_hs.s_mat, natoms)
    print(f"Sample H,S dense shape: {H_s.shape}")

    # basis for orbital projection kernel
    fdata_basis = os.path.join(_THIS_DIR, "Fdata", "basis")
    fdata_dir = os.path.join(_THIS_DIR, "Fdata")

    # scan setup
    size = 20.0
    n = 80
    z_list = [2.0, 10.0]

    # exponential vacuum radial params (used only for exp-kernel)
    rcut_exp  = 20.0
    beta = 1.0
    r0 = 3.0

    # tip metadata (single H atom, 1 s-orbital)
    tip_types = np.array([1], dtype=np.int32)
    tip_norb_per = np.array([1], dtype=np.int32)
    tip_pos_base = np.zeros((1, 3), dtype=np.float64)

    for z in z_list:
        outdir = set_export_dir(os.path.join(_THIS_DIR, "export", "stm_orbital_projection", f"ptcda_z{z:.1f}"))
        print(f"\n# z={z}Å  outdir={outdir}")

        X, Y, points_plane, extent_abs = _build_xy_grid(mol.apos, z=z, size=size, n=n)

        # tip positions are the same as points_plane
        points_tip = points_plane.copy().astype(np.float64)

        for mo0 in mo_list:
            label = _mo_label(mo0, homo, lumo)
            E = float(eigen[mo0])
            print(f"  MO {mo0+1:4d} {label:8s} E={E:+.4f} eV")

            # (1) Fortran reference ψ (orb2points)
            psi_fortran_flat = fc.orb2points(points_plane.astype(np.float64), iMO=int(mo0+1), ikpoint=1)
            psi_fortran = np.asarray(psi_fortran_flat, dtype=np.float64).reshape(n, n).T

            # (2) OpenCL true basis ψ (finite cutoff from basis)
            psi_opencl_flat = project_orbital_to_points(
                C_fc, int(mo0),
                mol.atypes, mol.apos,
                orb2atom, norb_per_smp,
                fdata_basis, points=points_plane.astype(np.float32)
            )
            psi_opencl = np.asarray(psi_opencl_flat, dtype=np.float64).reshape(n, n).T

            # (3) OpenCL exponential extended ψ (vacuum tail)
            psi_exp_flat = project_orbital_to_points_exp(
                C_fc, int(mo0),
                mol.atypes, mol.apos,
                orb2atom, norb_per_smp,
                fdata_basis, points=points_plane.astype(np.float32),
                beta=beta, r0=r0, rcut=rcut_exp
            )
            psi_exp = np.asarray(psi_exp_flat, dtype=np.float64).reshape(n, n).T

            # (4) Response amplitude — GPU OpenCL with exponential coupling
            print(f"    computing response (GPU OpenCL exp)...")
            # Precompute on CPU (once per MO)
            z = E + 1j * 1e-6
            A_ss = z * S_s.astype(np.complex64) - H_s.astype(np.complex64)
            G0 = np.linalg.inv(A_ss).astype(np.complex64)
            v = (C_fc[:, mo0].astype(np.complex64) @ G0)

            starts_s = _atom_orb_starts(norb_per_smp)
            atoms_dict_s = {
                'pos': mol.apos.astype(np.float32),
                'Rcut': np.full(len(mol.apos), rcut_exp, dtype=np.float32),
                'type': mol.atypes.astype(np.int32),
            }

            gp = GridProjector(fdata_dir=fdata_dir)
            gp.load_basis(species_nz=sorted(set(mol.atypes.astype(int))))
            resp_exp = gp.response_amplitude_exp(
                points_plane.astype(np.float32),
                atoms_dict_s, norb_per_smp, starts_s,
                v, G0, E=E, eta=1e-6, E_tip=0.0,
                beta=beta, r0=r0, A_ss=-1.0, A_sp=-1.0, rcut=rcut_exp
            )
            resp_exp = resp_exp.reshape(n, n).T

            # (5) True-coupling response placeholder — CPU path with Fdata is too slow
            # TODO: implement OpenCL kernel with Fdata radial lookups for true coupling
            resp_true = np.zeros_like(resp_exp)

            fig, axes = plt.subplots(2, 3, figsize=(18, 10))
            fig.suptitle(f"PTCDA MO {mo0+1} ({label}) E={E:+.4f} eV  z={z:.1f}Å")

            # Row 1: signed orbital projection (bwr)
            def _panel(ax, Z, title):
                vmax = max(abs(np.min(Z)), abs(np.max(Z)))
                if vmax < 1e-30:
                    vmax = 1.0
                im = ax.imshow(Z, origin='lower', extent=extent_abs, cmap='bwr', vmin=-vmax, vmax=vmax, aspect='equal')
                ax.set_title(title)
                ax.set_xlabel('x (Å)'); ax.set_ylabel('y (Å)')
                plot_atoms(ax, mol.apos, mol.atypes, color='green', ms=3)
                plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

            _panel(axes[0, 0], psi_fortran, "ψ(r) Fortran orb2points")
            _panel(axes[0, 1], psi_opencl,  "ψ(r) OpenCL true basis")
            _panel(axes[0, 2], psi_exp,     f"ψ(r) OpenCL exp tail (β={beta})")

            # Row 2: response amplitude (hot, positive only)
            def _panel_resp(ax, Z, title):
                vmax = np.max(Z)
                if vmax < 1e-30:
                    vmax = 1.0
                im = ax.imshow(Z, origin='lower', extent=extent_abs, cmap='hot', vmin=0, vmax=vmax, aspect='equal')
                ax.set_title(title)
                ax.set_xlabel('x (Å)'); ax.set_ylabel('y (Å)')
                plot_atoms(ax, mol.apos, mol.atypes, color='green', ms=3)
                plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

            _panel_resp(axes[1, 0], np.zeros_like(resp_true), "[reserved]")
            _panel_resp(axes[1, 1], resp_true, "Response |C^T x_s|^2 (true coupling)")
            _panel_resp(axes[1, 2], resp_exp,  "Response |C^T x_s|^2 (exp coupling)")

            plt.tight_layout()
            save_plot(fig, f"mo{mo0+1:04d}_{label}_panels", outdir, dpi=140)

    print("Done.")


if __name__ == "__main__":
    main()
