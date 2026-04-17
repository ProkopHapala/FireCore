#!/usr/bin/env python3
"""
Compare Fortran orb2points vs OpenCL orbital projection for PTCDA HOMO.

This script:
1. Uses Fortran orb2points as reference (correct implementation)
2. Uses OpenCL GridProjector with remapped coefficients
3. Plots signed wavefunctions side by side (not squared)
4. Uses plane at z=1Å above molecular plane (pi orbital is zero at atom plane)
"""

import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.path.dirname(__file__), "..", ".."))

from pyBall.AtomicSystem import AtomicSystem
from pyBall import FireCore as fc
from pyBall.FireballOCL.STM_utils import set_export_dir, save_plot, plot_atoms, project_orbital_to_grid_v2


def main():
    parser = argparse.ArgumentParser(description="Compare Fortran vs OpenCL orbital projection")
    parser.add_argument("--xyz", default="../../cpp/common_resources/xyz/PTCDA.xyz",
                       help="Path to XYZ file")
    parser.add_argument("--mo", type=int, default=None,
                       help="MO index (1-based). If None, use HOMO")
    parser.add_argument("--z", type=float, default=1.0,
                       help="Height above molecular plane for projection (Å)")
    parser.add_argument("--size", type=float, default=20.0,
                       help="Grid size in Å")
    parser.add_argument("--n", type=int, default=80,
                       help="Grid resolution (nx=ny=n)")
    parser.add_argument("--nmax_scf", type=int, default=200,
                       help="Max SCF iterations")
    parser.add_argument("--verbosity", type=int, default=0,
                       help="Fireball verbosity")
    parser.add_argument("--outdir", default="export/orbital_projection_compare",
                       help="Output directory")
    args = parser.parse_args()

    export_dir = set_export_dir(args.outdir)

    # Load molecule
    print(f"Loading molecule from {args.xyz}")
    mol = AtomicSystem(fname=args.xyz)
    print(f"Loaded {len(mol.apos)} atoms")

    # Initialize Fireball
    print("Initializing Fireball...")
    fc.setVerbosity(args.verbosity)
    fc.initialize(atomType=mol.atypes, atomPos=mol.apos, verbosity=args.verbosity)

    # Get dimensions
    dims = fc.get_HS_dims()
    norb = dims.norbitals
    nelec = dims.nelec
    homo_idx = nelec // 2
    print(f"norb={norb}, nelec={nelec}, HOMO index (1-based)={homo_idx}")

    # Run SCF
    print(f"Running SCF (nmax_scf={args.nmax_scf})...")
    forces, energies = fc.evalForce(mol.apos, nmax_scf=args.nmax_scf)
    print(f"SCF done. Etot={energies[0]:.6f} eV")

    # Get eigenvalues
    eigen = fc.get_eigen(ikp=1, norb=norb)
    print(f"Got {len(eigen)} eigenvalues")

    # Select MO
    if args.mo is None:
        mo_idx = homo_idx
    else:
        mo_idx = args.mo
    E_mo = eigen[mo_idx - 1]
    print(f"Selected MO {mo_idx} (HOMO{mo_idx - homo_idx:+d}) with E={E_mo:.4f} eV")

    # Build grid in molecular plane
    print("Building projection grid...")
    origin = mol.apos.mean(axis=0)
    normal = np.array([0.0, 0.0, 1.0])  # Assume z-normal for planar molecules

    # Create orthonormal in-plane basis
    trial = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(trial, normal)) > 0.9:
        trial = np.array([0.0, 1.0, 0.0])
    x_axis = trial - np.dot(trial, normal) * normal
    x_axis /= np.linalg.norm(x_axis)
    y_axis = np.cross(normal, x_axis)
    y_axis /= np.linalg.norm(y_axis)

    # Grid parameters - use same for both Fortran and OpenCL
    step = args.size / args.n
    xs = np.linspace(-args.size/2, args.size/2, args.n)
    ys = np.linspace(-args.size/2, args.size/2, args.n)
    X, Y = np.meshgrid(xs, ys)

    # Generate 3D points at height z above plane
    plane_origin = origin + args.z * normal
    points = plane_origin[None, :] + X[..., None] * x_axis[None, None, :] + Y[..., None] * y_axis[None, None, :]
    points = points.reshape(-1, 3)

    print(f"Grid: {args.n}x{args.n} points, step={step:.3f} Å, z={args.z} Å above plane")

    # ═══ Fortran orb2points (reference) ═══════════════════════════════════════
    print(f"\nProjecting MO {mo_idx} using Fortran orb2points...")
    # Fortran expects C-contiguous array, but with swapped strides (transpose the grid)
    # Create 2D grid, transpose it, then flatten to get correct stride order for raw pointer
    points_2d = plane_origin[None, None, :] + X[:, :, None] * x_axis[None, None, :] + Y[:, :, None] * y_axis[None, None, :]
    points_2d_T = points_2d.transpose(1, 0, 2)  # Swap x,y strides
    points_c = np.ascontiguousarray(points_2d_T.reshape(-1, 3), dtype=np.float64)
    mo_fortran = fc.orb2points(points_c, iMO=mo_idx, ikpoint=1)
    mo_fortran_2d = mo_fortran.reshape(args.n, args.n).T  # Transpose back for plotting

    print(f"Fortran range: [{mo_fortran_2d.min():.3e}, {mo_fortran_2d.max():.3e}]")

    # ═══ OpenCL GridProjector with remapped coefficients ═════════════════════════
    print(f"\nProjecting MO {mo_idx} using OpenCL GridProjector...")

    # Get MO coefficients directly from Fireball (avoid running SCF again)
    print("Getting MO coefficients from Fireball...")
    C_fc = fc.get_wfcoef(norb=norb)
    print(f"Got MO coefficients shape: {C_fc.shape}")

    # Build authoritative orbital mapping from Fireball sparse H/S metadata
    data = fc.get_HS_neighs(dims)
    data = fc.get_HS_sparse(dims, data)
    nzx = np.array(data.nzx, dtype=np.int32)
    iatyp = np.array(data.iatyp, dtype=np.int32)
    num_orb = np.array(data.num_orb, dtype=np.int32)
    norb_per = np.zeros(dims.natoms, dtype=np.int32)
    for ia in range(dims.natoms):
        Z = int(iatyp[ia])
        w = np.where(nzx == Z)[0]
        if w.size == 0:
            raise RuntimeError(f"Cannot map atom Z={Z} to nzx species list {nzx}")
        norb_per[ia] = int(num_orb[w[0]])
    starts = np.zeros(dims.natoms + 1, dtype=np.int32)
    starts[1:] = np.cumsum(norb_per)
    if int(starts[-1]) != norb:
        raise RuntimeError(f"Orbital count mismatch: mapped {int(starts[-1])} vs dims.norbitals={norb}")
    orb2atom = np.array([ia for ia in range(dims.natoms) for _ in range(int(norb_per[ia]))], dtype=np.int32)
    print(f"Orbital mapping: norb_per={norb_per[:5]}..., total={sum(norb_per)}")

    # Extract MO coefficients and project - use same grid as Fortran
    fdata_basis = os.path.join(os.path.dirname(__file__), "Fdata", "basis")

    # Build grid spec that matches Fortran grid exactly
    # Fortran grid: centered at origin, size x size, n x n points at z=args.z
    # Use ngrid=[80,80,8] for task building (build_tasks needs minimum z thickness)
    # Grid should span from z=args.z to z=args.z+step*8, take middle slice at z=args.z+step*4
    grid_origin = plane_origin - np.array([args.size/2, args.size/2, 0.0])
    grid_spec = {
        'origin': grid_origin,
        'dA': [step, 0, 0],
        'dB': [0, step, 0],
        'dC': [0, 0, step],
        'ngrid': [args.n, args.n, 8]  # z=8 for task building
    }

    # Use new orbital projection kernel (not density kernel)
    psi_opencl = project_orbital_to_grid_v2(
        C_fc, mo_idx - 1,  # 0-based index
        mol.atypes, mol.apos,
        orb2atom, norb_per,
        fdata_basis, grid_spec=grid_spec
    )

    print(f"OpenCL output shape: {psi_opencl.shape}")

    # Take z-slice at z=args.z (slice index = args.z/step = 1.0/0.25 = 4)
    z_slice_idx = int(args.z / 0.25)  # Fixed step for comparison
    psi_opencl_2d = psi_opencl[:, :, z_slice_idx].T

    print(f"OpenCL signed range: [{psi_opencl_2d.min():.3e}, {psi_opencl_2d.max():.3e}]")
    print(f"OpenCL non-zero count: {np.count_nonzero(psi_opencl_2d)} / {psi_opencl_2d.size}")

    # ═══ Comparison plots - SIGNED WAVEFUNCTIONS SIDE BY SIDE ═══════════════════
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle(f"MO {mo_idx} (HOMO{mo_idx - homo_idx:+d}) E={E_mo:.4f} eV at z={args.z}Å - Fortran vs OpenCL (SIGNED)")

    vmax_f = max(abs(mo_fortran_2d.min()), abs(mo_fortran_2d.max()))
    vmax_o = max(abs(psi_opencl_2d.min()), abs(psi_opencl_2d.max()))

    # Fortran
    im0 = axes[0].imshow(mo_fortran_2d, origin='lower',
                       extent=[xs.min(), xs.max(), ys.min(), ys.max()],
                       cmap='bwr', vmin=-vmax_f, vmax=vmax_f, aspect='equal')
    plot_atoms(axes[0], mol.apos, mol.atypes, color='green', ms=4)
    axes[0].set_title("Fortran: ψ(r)")
    axes[0].set_xlabel("x (Å)")
    axes[0].set_ylabel("y (Å)")
    plt.colorbar(im0, ax=axes[0])

    # OpenCL
    im1 = axes[1].imshow(psi_opencl_2d, origin='lower',
                       extent=[xs.min(), xs.max(), ys.min(), ys.max()],
                       cmap='bwr', vmin=-vmax_o, vmax=vmax_o, aspect='equal')
    plot_atoms(axes[1], mol.apos, mol.atypes, color='green', ms=4)
    axes[1].set_title("OpenCL: ψ(r) (remapped)")
    axes[1].set_xlabel("x (Å)")
    axes[1].set_ylabel("y (Å)")
    plt.colorbar(im1, ax=axes[1])

    # Difference
    diff = mo_fortran_2d - psi_opencl_2d
    vmax_d = max(abs(diff.min()), abs(diff.max()))
    im2 = axes[2].imshow(diff, origin='lower',
                       extent=[xs.min(), xs.max(), ys.min(), ys.max()],
                       cmap='bwr', vmin=-vmax_d, vmax=vmax_d, aspect='equal')
    axes[2].set_title(f"Difference (max|diff|={vmax_d:.3e})")
    axes[2].set_xlabel("x (Å)")
    axes[2].set_ylabel("y (Å)")
    plt.colorbar(im2, ax=axes[2])

    plt.tight_layout()
    save_plot(fig, f"mo{mo_idx:04d}_compare_signed", export_dir, dpi=140)

    # Correlation plot
    fig2, ax2 = plt.subplots(1, 1, figsize=(6, 6))
    ax2.scatter(mo_fortran_2d.ravel(), psi_opencl_2d.ravel(), alpha=0.3, s=1)
    ax2.plot([mo_fortran_2d.min(), mo_fortran_2d.max()],
             [mo_fortran_2d.min(), mo_fortran_2d.max()], 'r--')
    ax2.set_xlabel("Fortran ψ")
    ax2.set_ylabel("OpenCL ψ")
    ax2.set_title(f"Correlation (r={np.corrcoef(mo_fortran_2d.ravel(), psi_opencl_2d.ravel())[0,1]:.6f})")
    ax2.axis('equal')
    plt.tight_layout()
    save_plot(fig2, f"mo{mo_idx:04d}_correlation", export_dir, dpi=140)

    # Save numerical data
    np.savez(os.path.join(export_dir, f"mo{mo_idx:04d}_compare_signed.npz"),
             X=X, Y=Y, Z_fortran=mo_fortran_2d, Z_opencl=psi_opencl_2d,
             mo_idx=mo_idx, E=E_mo, z=args.z)

    print(f"\nCorrelation coefficient: {np.corrcoef(mo_fortran_2d.ravel(), psi_opencl_2d.ravel())[0,1]:.6f}")
    print(f"RMS difference: {np.sqrt(np.mean((mo_fortran_2d - psi_opencl_2d)**2)):.3e}")
    print(f"Saved to {export_dir}")
    print("Done.")


if __name__ == "__main__":
    main()
