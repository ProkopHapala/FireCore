#!/usr/bin/env python3
r"""H2O MO vs Green's-function LDOS projection test.

Goal:
- Verify Fireball sparse H,S mapping to dense AO matrix (6x6 for H2O)
- Compare orbital wavefunction \psi_n(r) against LDOS(r;E=eps_n)

Outputs:
- export/h2o_mo_vs_ldos/moXXXX_panels.png  (2-panel: MO and LDOS)

Run from:
  cd tests/pyFireball
  python test_h2o_mo_vs_ldos.py
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.path.dirname(__file__), "..", ".."))

from pyBall import FireCore as fc
from pyBall.atomicUtils import load_xyz
from pyBall.FireballOCL.Grid import GridProjector
from pyBall.FireballOCL.STM_utils import set_export_dir, save_plot, plot_atoms, project_orbital_to_points

np.set_printoptions(precision=6, suppress=True, linewidth=np.inf)

_THIS_DIR = os.path.dirname(__file__)


def _orbital_layout(sparse_data, natoms):
    nzx = np.array(sparse_data.nzx, dtype=np.int32)
    iatyp = np.array(sparse_data.iatyp, dtype=np.int32)
    num_orb = np.array(sparse_data.num_orb, dtype=np.int32)
    n_orb_atom = np.zeros(natoms, dtype=np.int32)
    for ia in range(natoms):
        Z = int(iatyp[ia])
        w = np.where(nzx == Z)[0]
        if w.size == 0:
            raise RuntimeError(f"Cannot map atom Z={Z} to nzx={nzx}")
        n_orb_atom[ia] = int(num_orb[int(w[0])])
    offs = np.zeros(natoms + 1, dtype=np.int32)
    offs[1:] = np.cumsum(n_orb_atom)
    return n_orb_atom, offs


def _blocked_to_dense(sparse_data, H_blocks, natoms):
    """Map FireCore blocked sparse -> dense AO matrix.

    FireCore Python wrapper keeps the *linear* Fortran layout by reversing axes,
    so the natural view is:
      H_blocks[iatom, ineigh, imu, inu]
    where imu runs orbitals on iatom (row) and inu runs orbitals on jatom (col).
    This is exactly what the OpenCL projector kernels expect as well (row-major
    4x4 blocks when cast to float4).
    """
    n_orb_atom, offs = _orbital_layout(sparse_data, natoms)
    norb = int(offs[-1])
    M = np.zeros((norb, norb), dtype=np.float64)
    neighn = np.array(sparse_data.neighn, dtype=np.int32)
    neigh_j = np.array(sparse_data.neigh_j, dtype=np.int32)
    neigh_b = np.array(sparse_data.neigh_b, dtype=np.int32)

    neigh_self = np.full(natoms, -1, dtype=np.int32)
    for i in range(natoms):
        ii = i + 1
        for ineigh in range(int(neighn[i])):
            if int(neigh_j[i, ineigh]) == ii and int(neigh_b[i, ineigh]) == 0:
                neigh_self[i] = ineigh
                break

    for i in range(natoms):
        ni = int(n_orb_atom[i]); i0 = int(offs[i])
        for ineigh in range(int(neighn[i])):
            if ineigh == int(neigh_self[i]):
                j = i
            else:
                j = int(neigh_j[i, ineigh]) - 1
            if j < 0 or j >= natoms:
                continue
            nj = int(n_orb_atom[j]); j0 = int(offs[j])
            # H_blocks last axes are (imu, inu) => (ni, nj)
            blk = H_blocks[i, ineigh, :ni, :nj]
            M[i0:i0+ni, j0:j0+nj] += blk
    return M


def _blocked_to_dense_T(sparse_data, H_blocks, natoms):
    """Same as _blocked_to_dense but assumes blocks are stored as (inu, imu) and transposes last 2 axes."""
    n_orb_atom, offs = _orbital_layout(sparse_data, natoms)
    norb = int(offs[-1])
    M = np.zeros((norb, norb), dtype=np.float64)
    neighn = np.array(sparse_data.neighn, dtype=np.int32)
    neigh_j = np.array(sparse_data.neigh_j, dtype=np.int32)
    neigh_b = np.array(sparse_data.neigh_b, dtype=np.int32)

    neigh_self = np.full(natoms, -1, dtype=np.int32)
    for i in range(natoms):
        ii = i + 1
        for ineigh in range(int(neighn[i])):
            if int(neigh_j[i, ineigh]) == ii and int(neigh_b[i, ineigh]) == 0:
                neigh_self[i] = ineigh
                break

    for i in range(natoms):
        ni = int(n_orb_atom[i]); i0 = int(offs[i])
        for ineigh in range(int(neighn[i])):
            if ineigh == int(neigh_self[i]):
                j = i
            else:
                j = int(neigh_j[i, ineigh]) - 1
            if j < 0 or j >= natoms:
                continue
            nj = int(n_orb_atom[j]); j0 = int(offs[j])
            blk = H_blocks[i, ineigh, :nj, :ni].T
            M[i0:i0+ni, j0:j0+nj] += blk
    return M


def _dense_to_sparse_blocks(sparse_data, M_dense, natoms, numorb_max):
    """Pack dense AO matrix into sparse blocks M_blocks[iatom,ineigh,nu,mu] matching Grid.cl kernel expectations."""
    n_orb_atom, offs = _orbital_layout(sparse_data, natoms)
    neighn = np.array(sparse_data.neighn, dtype=np.int32)
    neigh_j = np.array(sparse_data.neigh_j, dtype=np.int32)
    neigh_b = np.array(sparse_data.neigh_b, dtype=np.int32)

    neigh_self = np.full(natoms, -1, dtype=np.int32)
    for i in range(natoms):
        ii = i + 1
        for ineigh in range(int(neighn[i])):
            if int(neigh_j[i, ineigh]) == ii and int(neigh_b[i, ineigh]) == 0:
                neigh_self[i] = ineigh
                break

    neigh_max = int(M_dense.shape[0])  # not used
    neigh_max = int(sparse_data.neigh_j.shape[1])

    blocks = np.zeros((natoms, neigh_max, numorb_max, numorb_max), dtype=np.float32)

    for i in range(natoms):
        ni = int(n_orb_atom[i]); i0 = int(offs[i])
        for ineigh in range(int(neighn[i])):
            if ineigh == int(neigh_self[i]):
                j = i
            else:
                j = int(neigh_j[i, ineigh]) - 1
            if j < 0 or j >= natoms:
                continue
            nj = int(n_orb_atom[j]); j0 = int(offs[j])
            # Store as (imu, inu) so the kernel reads float4 rows correctly
            blk = M_dense[i0:i0+ni, j0:j0+nj]  # (imu_i, inu_j)
            blocks[i, ineigh, :ni, :nj] = blk.astype(np.float32)

    return blocks


def _dense_vec_to_atom4(Ccol, n_orb_atom, offs):
    """Pack AO coefficients (global AO order) into per-atom (natoms,4) coeffs for OpenCL orbital projector.

    Assumptions for H2O:
    - O has 4 orbitals ordered [s,px,py,pz]
    - H has 1 orbital [s]
    OpenCL expects [px,py,pz,s].
    """
    natoms = len(n_orb_atom)
    coeffs = np.zeros((natoms, 4), dtype=np.float32)
    for ia in range(natoms):
        no = int(n_orb_atom[ia])
        i0 = int(offs[ia])
        if no == 4:
            s  = float(Ccol[i0+0])
            px = float(Ccol[i0+1])
            py = float(Ccol[i0+2])
            pz = float(Ccol[i0+3])
            coeffs[ia, 0] = px
            coeffs[ia, 1] = py
            coeffs[ia, 2] = pz
            coeffs[ia, 3] = s
        elif no == 1:
            coeffs[ia, 3] = float(Ccol[i0])
        else:
            raise RuntimeError(f"_dense_vec_to_atom4: unsupported norb={no} on atom {ia} (this test assumes H2O)")
    return coeffs


def _build_plane_grid(atomPos, z=2.0, size=8.0, n=128):
    cen = atomPos.mean(axis=0)
    xs = np.linspace(cen[0] - size * 0.5, cen[0] + size * 0.5, n)
    ys = np.linspace(cen[1] - size * 0.5, cen[1] + size * 0.5, n)
    X, Y = np.meshgrid(xs, ys, indexing='ij')
    pts = np.zeros((n * n, 3), dtype=np.float64)
    pts[:, 0] = X.ravel()
    pts[:, 1] = Y.ravel()
    pts[:, 2] = cen[2] + float(z)
    extent = [xs[0], xs[-1], ys[0], ys[-1]]
    return X, Y, pts, extent


def main():
    xyz = os.path.join(_THIS_DIR, "..", "..", "cpp", "common_resources", "xyz", "H2O.xyz")
    atomPos, atomTypes, enames, _, _comment = load_xyz(xyz)
    atomTypes = np.array(atomTypes, dtype=np.int32)
    atomPos = np.array(atomPos, dtype=np.float64)

    outdir = set_export_dir(os.path.join(_THIS_DIR, "export", "h2o_mo_vs_ldos"))

    fc.setVerbosity(0)
    fc.initialize(atomType=atomTypes, atomPos=atomPos, verbosity=0)

    # Do one SCF (small system)
    fc.evalForce(atomPos, nmax_scf=50)

    # Make sparse HS export include VNL, to match ktransform() used by Fortran Green's function.
    # export_mode: 0/1 without VNL, 2 with VNL added into h_mat_out where possible.
    fc.set_export_mode(2)

    dims = fc.get_HS_dims()
    data = fc.get_HS_neighs(dims)
    data = fc.get_HS_sparse(dims, data)

    natoms = int(dims.natoms)
    norb = int(dims.norbitals)
    assert natoms == 3
    assert norb == 6

    # Use mapping that best matches Fortran ktransform(Gamma) for diagnostic printing
    Hk_fc, Sk_fc = fc.get_HS_k(np.array((0.0, 0.0, 0.0), dtype=np.float64), norb)
    Hk_fc = np.asarray(Hk_fc, dtype=np.complex128)
    Sk_fc = np.asarray(Sk_fc, dtype=np.complex128)
    H0 = _blocked_to_dense(data, data.h_mat, natoms)
    S0 = _blocked_to_dense(data, data.s_mat, natoms)
    H1 = _blocked_to_dense_T(data, data.h_mat, natoms)
    S1 = _blocked_to_dense_T(data, data.s_mat, natoms)
    dH0 = float(np.max(np.abs(H0 - np.real(Hk_fc)))); dS0 = float(np.max(np.abs(S0 - np.real(Sk_fc))))
    dH1 = float(np.max(np.abs(H1 - np.real(Hk_fc)))); dS1 = float(np.max(np.abs(S1 - np.real(Sk_fc))))
    use_T = (dH1 + dS1) < (dH0 + dS0)
    H = H1 if use_T else H0
    S = S1 if use_T else S0

    print("H (6x6) =\n", H)
    print("S (6x6) =\n", S)
    eig = fc.get_eigen(ikp=1, norb=norb)
    print("eigen (eV)=", np.round(eig, 6))

    C = fc.get_wfcoef(norb=norb)

    # mapping for OpenCL basis/projectors
    n_orb_atom, offs = _orbital_layout(data, natoms)
    norb_per = n_orb_atom.astype(np.int32)

    fdata_dir = os.path.join(_THIS_DIR, "Fdata")

    # grid plane
    z_plane = 2.0
    nxy = 160
    size = 8.0
    X, Y, pts, extent = _build_plane_grid(atomPos, z=z_plane, size=size, n=nxy)

    gp = GridProjector(fdata_dir=fdata_dir)
    gp.verbosity = 1
    gp.load_basis(species_nz=sorted(set(atomTypes.tolist())))

    atoms_dict = {
        'pos': atomPos.astype(np.float32),
        'Rcut': np.full(natoms, 6.0, dtype=np.float32),
        'type': atomTypes.astype(np.int32),
    }

    xs = np.linspace(extent[0], extent[1], nxy)
    ys = np.linspace(extent[2], extent[3], nxy)
    dx = float(xs[1] - xs[0])
    dy = float(ys[1] - ys[0])
    origin = np.array([xs[0], ys[0], atomPos[:, 2].mean() + z_plane], dtype=np.float32)
    grid_spec = {
        'origin': origin,
        'dA': np.array([dx, 0.0, 0.0], dtype=np.float32),
        'dB': np.array([0.0, dy, 0.0], dtype=np.float32),
        'dC': np.array([0.0, 0.0, 1.0], dtype=np.float32),
        'ngrid': np.array([nxy, nxy, 1], dtype=np.int32),
    }

    # Choose orbitals to plot (all 6)
    mo_list = list(range(norb))

    eta = 1e-2  # eV

    n_orb_atom, offs = _orbital_layout(data, natoms)

    # --- HS mapping parity check vs Fortran ktransform(Gamma)
    Hk_fc, Sk_fc = fc.get_HS_k(np.array((0.0, 0.0, 0.0), dtype=np.float64), norb)
    Hk_fc = np.asarray(Hk_fc, dtype=np.complex128)
    Sk_fc = np.asarray(Sk_fc, dtype=np.complex128)
    Hs0 = _blocked_to_dense(data, data.h_mat, natoms)
    Ss0 = _blocked_to_dense(data, data.s_mat, natoms)
    Hs1 = _blocked_to_dense_T(data, data.h_mat, natoms)
    Ss1 = _blocked_to_dense_T(data, data.s_mat, natoms)
    dH0 = float(np.max(np.abs(Hs0 - np.real(Hk_fc))))
    dS0 = float(np.max(np.abs(Ss0 - np.real(Sk_fc))))
    dH1 = float(np.max(np.abs(Hs1 - np.real(Hk_fc))))
    dS1 = float(np.max(np.abs(Ss1 - np.real(Sk_fc))))
    print(f"[HS PARITY] direct : max|H-Hk|={dH0:.3e}  max|S-Sk|={dS0:.3e}")
    print(f"[HS PARITY] transT : max|H-Hk|={dH1:.3e}  max|S-Sk|={dS1:.3e}")
    use_T = (dH1 + dS1) < (dH0 + dS0)
    print(f"[HS PARITY] selected mapping = {'transpose_last2' if use_T else 'direct'}")

    for mo in mo_list:
        E = float(eig[mo])

        # --- Fortran reference panels
        psiF = np.asarray(fc.orb2points(pts.astype(np.float64), iMO=int(mo + 1), ikpoint=1), dtype=np.float64).reshape(nxy, nxy).T
        ldosF = np.asarray(fc.ldos2points(pts.astype(np.float64), E=E, eta=eta, mode=1, iMO=int(mo + 1), ikpoint=1), dtype=np.float64).reshape(nxy, nxy).T

        # --- Python Green's function from sparse->dense mapping (checkpoint) + OpenCL projection
        if use_T:
            Hs = _blocked_to_dense_T(data, data.h_mat, natoms)
            Ss = _blocked_to_dense_T(data, data.s_mat, natoms)
        else:
            Hs = _blocked_to_dense(data, data.h_mat, natoms)
            Ss = _blocked_to_dense(data, data.s_mat, natoms)
        z = E + 1j * eta
        Gnp = np.linalg.inv(z * Ss.astype(np.complex128) - Hs.astype(np.complex128))

        # Compare against Fortran G(k) (Gamma)
        Gfc = fc.get_G_k(kpoint=(0.0, 0.0, 0.0), E=E, eta=eta, norb=norb)
        dG = np.max(np.abs(Gnp - Gfc))
        print(f"MO {mo+1}  E={E:+.6f}  max|Gnp-Gfc|={dG:.3e}")

        spec = -(1.0 / np.pi) * np.imag(Gnp)
        rho_blocks = _dense_to_sparse_blocks(data, spec, natoms, int(dims.numorb_max))
        tasks_np, task_atoms_np = gp.build_tasks(atoms_dict, grid_spec, nMaxAtom=64)
        ldosCL3d = gp.project(rho_blocks, data, atoms_dict, grid_spec, tasks=(tasks_np, task_atoms_np), nMaxAtom=64, use_gpu_tasks=False, use_tiled=True)
        ldosCL = np.asarray(ldosCL3d[:, :, 0], dtype=np.float64).T

        # OpenCL orbital projection: reuse the already-validated packing/mapping helper
        # (see tests/pyFireball/test_stm_orbital_projection.py)
        orb2atom = np.array([ia for ia in range(natoms) for _ in range(int(n_orb_atom[ia]))], dtype=np.int32)
        psiCL_flat = project_orbital_to_points(C, int(mo), atomTypes, atomPos, orb2atom, n_orb_atom, os.path.join(fdata_dir, 'basis'), points=pts.astype(np.float32))
        psiCL = np.asarray(psiCL_flat, dtype=np.float64).reshape(nxy, nxy).T

        fig, axes = plt.subplots(2, 2, figsize=(11, 8))
        fig.suptitle(f"H2O MO {mo+1}  E={E:+.4f} eV  z={z_plane:.1f} Å")

        vmax = float(max(abs(np.min(psiF)), abs(np.max(psiF)), abs(np.min(psiCL)), abs(np.max(psiCL))))
        if vmax < 1e-12: vmax = 1.0
        im00 = axes[0,0].imshow(psiF, origin='lower', extent=extent, cmap='bwr', vmin=-vmax, vmax=vmax, aspect='equal')
        axes[0,0].set_title("Fortran MO ψ(r)")
        plot_atoms(axes[0,0], atomPos, atomTypes, color='green', ms=6)
        plt.colorbar(im00, ax=axes[0,0], fraction=0.046, pad=0.04)

        vmaxF = float(np.max(ldosF));  vmaxF = max(vmaxF, 1e-16)
        im01 = axes[0,1].imshow(ldosF, origin='lower', extent=extent, cmap='hot', vmin=0.0, vmax=vmaxF, aspect='equal')
        axes[0,1].set_title("Fortran LDOS(r;E)")
        plot_atoms(axes[0,1], atomPos, atomTypes, color='green', ms=6)
        plt.colorbar(im01, ax=axes[0,1], fraction=0.046, pad=0.04)

        im10 = axes[1,0].imshow(psiCL, origin='lower', extent=extent, cmap='bwr', vmin=-vmax, vmax=vmax, aspect='equal')
        axes[1,0].set_title("OpenCL MO ψ(r)")
        plot_atoms(axes[1,0], atomPos, atomTypes, color='green', ms=6)
        plt.colorbar(im10, ax=axes[1,0], fraction=0.046, pad=0.04)

        vmaxCL = float(np.max(ldosCL)); vmaxCL = max(vmaxCL, 1e-16)
        im11 = axes[1,1].imshow(ldosCL, origin='lower', extent=extent, cmap='hot', vmin=0.0, vmax=vmaxCL, aspect='equal')
        axes[1,1].set_title("OpenCL LDOS(r;E)")
        plot_atoms(axes[1,1], atomPos, atomTypes, color='green', ms=6)
        plt.colorbar(im11, ax=axes[1,1], fraction=0.046, pad=0.04)

        for ax in axes.ravel():
            ax.set_xlabel('x (Å)'); ax.set_ylabel('y (Å)')

        plt.tight_layout()
        save_plot(fig, f"mo{mo+1:04d}_parity", outdir, dpi=160)

    print("Done.")


if __name__ == "__main__":
    main()
