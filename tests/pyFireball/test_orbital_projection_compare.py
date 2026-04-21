#!/usr/bin/env python3
"""
Compare Fortran orb2points vs OpenCL orbital projection for any molecule.

This script:
1. Uses Fortran orb2points as reference (correct implementation)
2. Uses OpenCL project_orbital_points kernel with correct coefficient remapping
3. Plots signed wavefunctions side by side (not squared)
4. Computes correlation statistics between methods
5. Uses plane at z-height above molecular plane

Usage:
    python test_orbital_projection_compare.py --xyz PTCDA.xyz
    python test_orbital_projection_compare.py --xyz H2O.xyz --mo 5 --z 1.0
    python test_orbital_projection_compare.py --xyz CH4.xyz --orbitals 0-7
"""

import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.path.dirname(__file__), "..", ".."))

from pyBall.AtomicSystem import AtomicSystem
from pyBall import FireCore as fc
from pyBall.FireballOCL.STM_utils import (
    set_export_dir, save_plot, plot_atoms, compute_correlation_stats,
    write_scaling_summary, plot_orbital_comparison, parse_orbital_indices,
    generate_mo_labels, project_orbital_to_points, build_grid_spec
)


def main():
    parser = argparse.ArgumentParser(description="Compare Fortran vs OpenCL orbital projection")
    parser.add_argument("--xyz", default="../../cpp/common_resources/xyz/PTCDA.xyz",
                       help="Path to XYZ file")
    parser.add_argument("--mo", type=int, default=None,
                       help="MO index (1-based). If None, use HOMO")
    parser.add_argument("--orbitals", type=str, default=None,
                       help='MO indices to compare: "0,1,2" or "0-5". If None, uses --mo or HOMO')
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
    parser.add_argument("--save-npz", action="store_true",
                       help="Save numerical data to .npz file")
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
    print(f"norb={norb}")

    # Run SCF
    print(f"Running SCF (nmax_scf={args.nmax_scf})...")
    forces, energies = fc.evalForce(mol.apos, nmax_scf=args.nmax_scf)
    print(f"SCF done. Etot={energies[0]:.6f} eV")

    # Get eigenvalues
    eigen = fc.get_eigen(ikp=1, norb=norb)
    print(f"Got {len(eigen)} eigenvalues")

    # Determine HOMO from eigenvalues (last occupied = last negative energy)
    # HOMO is the highest energy orbital that is still occupied (negative or closest to zero)
    occupied = np.where(eigen < 0)[0]
    if len(occupied) > 0:
        homo_idx = occupied[-1] + 1  # Convert to 1-based
    else:
        # If no negative eigenvalues, use the lowest positive as LUMO, HOMO would be last "filled"
        # This is a fallback for unusual cases
        homo_idx = len(eigen) // 2
    
    nelec = homo_idx * 2  # Estimate nelec from HOMO position
    print(f"HOMO detected: MO {homo_idx} at {eigen[homo_idx-1]:.4f} eV")
    print(f"LUMO detected: MO {homo_idx+1} at {eigen[homo_idx]:.4f} eV")

    # Determine which MOs to compare
    if args.orbitals is not None:
        # Multiple orbitals specified
        mo_indices = parse_orbital_indices(args.orbitals, norb)
        mo_indices = [i + 1 for i in mo_indices]  # Convert to 1-based
        print(f"Selected {len(mo_indices)} orbitals: {mo_indices}")
    elif args.mo is not None:
        # Single MO specified
        mo_indices = [args.mo]
    else:
        # Default to HOMO
        mo_indices = [homo_idx]
    
    # Generate labels for all selected MOs
    mo_labels = generate_mo_labels([i - 1 for i in mo_indices], eigen, nelec, 
                                   mol_name=os.path.basename(args.xyz).replace('.xyz', ''))
    print(f"MO labels: {mo_labels}")

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
    # IMPORTANT: Use absolute coordinates for extent to match atom positions
    step = args.size / args.n
    # Grid bounds in absolute coordinates (centered on molecular origin)
    # Match test_h2o_orbital_comparison.py: grid_origin = atom_center - [Lx/2, Ly/2, Lz/2]
    grid_origin = origin - np.array([args.size/2, args.size/2, 0.0])
    
    xs = np.linspace(grid_origin[0] + step*0.5, grid_origin[0] + args.size - step*0.5, args.n)  # Voxel centers
    ys = np.linspace(grid_origin[1] + step*0.5, grid_origin[1] + args.size - step*0.5, args.n)
    X, Y = np.meshgrid(xs, ys, indexing='ij')  # IMPORTANT: indexing='ij' for correct orientation
    
    # Extent for plotting (absolute coordinates, matching atom positions)
    extent_abs = [grid_origin[0], grid_origin[0] + args.size, grid_origin[1], grid_origin[1] + args.size]

    # Generate 3D points at height z above plane - match test_h2o_orbital_comparison.py exactly
    z_height = args.z
    Z_plane = np.zeros_like(X) + origin[2] + z_height
    points_c = np.stack([X.ravel(), Y.ravel(), Z_plane.ravel()], axis=1)
    points_c = np.ascontiguousarray(points_c, dtype=np.float64)

    print(f"Grid: {args.n}x{args.n} points, step={step:.3f} Å, z={z_height} Å")
    print(f"  Extent: x=[{grid_origin[0]:.2f}, {grid_origin[0] + args.size:.2f}], y=[{grid_origin[1]:.2f}, {grid_origin[1] + args.size:.2f}]")

    # ═══ Get MO coefficients and orbital mapping (common for all MOs) ════════════
    print("\nGetting MO coefficients from Fireball...")
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

    fdata_basis = os.path.join(os.path.dirname(__file__), "Fdata", "basis")

    # ═══ Process each MO ═══════════════════════════════════════════════════════
    all_results = []
    
    for mo_idx, mo_label in zip(mo_indices, mo_labels):
        E_mo = eigen[mo_idx - 1]
        print(f"\n{'='*60}")
        print(f"Processing MO {mo_idx} ({mo_label})  E={E_mo:.4f} eV")
        print(f"{'='*60}")

        # ═══ Fortran orb2points (reference) ═══════════════════════════════════
        print(f"  Fortran orb2points...")
        mo_fortran_flat = fc.orb2points(points_c, iMO=mo_idx, ikpoint=1)
        mo_fortran_2d = mo_fortran_flat.reshape(args.n, args.n)  # NO transpose (matches test_h2o)
        print(f"  Fortran range: [{mo_fortran_2d.min():.3e}, {mo_fortran_2d.max():.3e}]")

        # ═══ OpenCL project_orbital_points (exact points kernel) ══════════════
        print(f"  OpenCL project_orbital_points (exact kernel)...")
        
        psi_opencl_flat = project_orbital_to_points(
            C_fc, mo_idx - 1,  # 0-based index
            mol.atypes, mol.apos,
            orb2atom, norb_per,
            fdata_basis, points=points_c.astype(np.float32)  # Same points as Fortran
        )
        psi_opencl_2d = psi_opencl_flat.reshape(args.n, args.n)  # NO transpose (matches test_h2o)
        
        print(f"  OpenCL range: [{psi_opencl_2d.min():.3e}, {psi_opencl_2d.max():.3e}]")

        # ═══ Comparison statistics ════════════════════════════════════════════
        stats = compute_correlation_stats(mo_fortran_2d, psi_opencl_2d)
        stats['mo_idx'] = mo_idx
        stats['mo_label'] = mo_label
        all_results.append(stats)
        
        print(f"  Correlation: {stats['corr']:.6f}")
        print(f"  Scale (linreg): {stats['scale_linreg']:.4f}")
        print(f"  RMS difference: {np.sqrt(np.mean((mo_fortran_2d - psi_opencl_2d)**2)):.3e}")

        # ═══ Comparison plots ═════════════════════════════════════════════════
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))
        fig.suptitle(f"MO {mo_idx} ({mo_label}) E={E_mo:.4f} eV at z={args.z}Å - Fortran vs OpenCL (EXACT POINTS)")

        vmax_f = max(abs(mo_fortran_2d.min()), abs(mo_fortran_2d.max()))
        vmax_o = max(abs(psi_opencl_2d.min()), abs(psi_opencl_2d.max()))

        # Fortran - use absolute coordinates for both extent and atom positions
        plot_orbital_comparison(axes[0], mo_fortran_2d, mol.apos, mol.atypes,
                              extent_abs,
                              f"Fortran: ψ(r)", enames=mol.enames)

        # OpenCL
        plot_orbital_comparison(axes[1], psi_opencl_2d, mol.apos, mol.atypes,
                              extent_abs,
                              f"OpenCL: ψ(r)", enames=mol.enames)

        # Difference
        diff = mo_fortran_2d - psi_opencl_2d
        vmax_d = max(abs(diff.min()), abs(diff.max()))
        im2 = axes[2].imshow(diff, origin='lower',
                           extent=extent_abs,
                           cmap='bwr', vmin=-vmax_d, vmax=vmax_d, aspect='equal')
        axes[2].set_title(f"Difference (max|diff|={vmax_d:.3e})")
        axes[2].set_xlabel("x (Å)")
        axes[2].set_ylabel("y (Å)")
        plt.colorbar(im2, ax=axes[2])
        plot_atoms(axes[2], mol.apos, mol.atypes, color='green', ms=4)

        plt.tight_layout()
        save_plot(fig, f"mo{mo_idx:04d}_{mo_label}_compare", export_dir, dpi=140)

        # Correlation plot
        fig2, ax2 = plt.subplots(1, 1, figsize=(6, 6))
        ax2.scatter(mo_fortran_2d.ravel(), psi_opencl_2d.ravel(), alpha=0.3, s=1)
        ax2.plot([mo_fortran_2d.min(), mo_fortran_2d.max()],
                 [mo_fortran_2d.min(), mo_fortran_2d.max()], 'r--')
        ax2.set_xlabel("Fortran ψ")
        ax2.set_ylabel("OpenCL ψ")
        ax2.set_title(f"Correlation (r={stats['corr']:.6f})")
        ax2.axis('equal')
        plt.tight_layout()
        save_plot(fig2, f"mo{mo_idx:04d}_{mo_label}_correlation", export_dir, dpi=140)

        # Save numerical data if requested
        if args.save_npz:
            np.savez(os.path.join(export_dir, f"mo{mo_idx:04d}_{mo_label}_data.npz"),
                     X=X, Y=Y, Z_fortran=mo_fortran_2d, Z_opencl=psi_opencl_2d,
                     mo_idx=mo_idx, mo_label=mo_label, E=E_mo, z=args.z, stats=stats)

    # ═══ Summary statistics ═══════════════════════════════════════════════════
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    
    for r in all_results:
        print(f"MO {r['mo_idx']:3d} ({r['mo_label']:8s}): Corr={r['corr']:8.6f}, Scale={r['scale_linreg']:7.4f}")
    
    avg_corr = np.mean([r['corr'] for r in all_results])
    avg_scale = np.mean([r['scale_linreg'] for r in all_results])
    print(f"\nAverage correlation: {avg_corr:.6f}")
    print(f"Average scale: {avg_scale:.4f}")

    # Write scaling summary file
    summary_path = write_scaling_summary(all_results, export_dir, 'scaling_summary.txt')
    print(f"\nSummary written to: {summary_path}")
    print(f"Plots saved to: {export_dir}")
    print("Done.")

    return all_results


if __name__ == "__main__":
    main()
