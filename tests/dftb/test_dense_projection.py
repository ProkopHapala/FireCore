#!/usr/bin/env python3
"""
test_dense_projection.py

Test dense-matrix projection methods vs old sparse methods for sp-only systems.
Uses point evaluation (orb2points/density2points style) for easy comparison.
Generates plots comparing dense vs sparse methods.

This validates that the new dense kernels produce identical results to the old
sparse kernels for systems without d-orbitals.

Usage:
    python test_dense_projection.py --dftb-dir tests/tAFM/pyocl_fdbm/test_pentacene_mio/dftb_work --plot
"""

import sys, os, argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

REPO_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

from pyBall.DFTB.DFTBplusParser import (parse_basis_hsd_ang, parse_detailed_xml_custom, 
                                         parse_eigenvec_bin_custom, evec_to_kernel_coeffs,
                                         parse_wfc_hsd, convert_wfc_to_species_list_ang)
from pyBall.DFTB.Grid_dftb import GridProjector, setup_gridprojector_from_dftb
from pyBall.DFTB.TestUtils import generate_2d_point_grid
from pyBall.plotUtils import plot_comparison_2d
from pyBall.FireballOCL.STM_utils import plot_orbital_comparison

BOHR2ANG = 0.5291772109

def parse_args():
    p = argparse.ArgumentParser(description='Test dense vs sparse projection',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('--dftb-dir', type=str, required=True,help='DFTB+ run dir containing detailed.xml, eigenvec.bin')
    p.add_argument('--basis', type=str, default='mio-1-1', choices=['mio-1-1', '3ob-3-1'],help='Basis set to use (mio-1-1 or 3ob-3-1)')
    p.add_argument('--basis-path', type=str, default=None, help='Path to basis files (default: pyBall/DFTB/data/)')
    p.add_argument('--tol', type=float, default=1e-6,  help='Tolerance for comparing dense vs sparse results')
    p.add_argument('--step', type=float, default=0.1,  help='Grid step in Angstrom for 2D plots (default: 0.1)')
    p.add_argument('--z-offset', type=float, default=2.0, help='Z offset in Angstrom for XY plane (default: 2.0)')
    p.add_argument('--plot', action='store_true', help='Generate comparison plots')
    p.add_argument('--output-dir', type=str, default=None,   help='Output directory for plots (default: same as dftb-dir)')
    p.add_argument('--dpi', type=int, default=150,help='DPI for output images (default: 150)')
    return p.parse_args()


def main():
    args = parse_args()
    dftb_dir = Path(args.dftb_dir)
    output_dir = Path(args.output_dir) if args.output_dir else dftb_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 70)
    print("Dense vs Sparse Projection Test")
    print("=" * 70)
    
    # Check required files
    detailed_xml = dftb_dir / 'detailed.xml'
    eigenvec_bin = dftb_dir / 'eigenvec.bin'
    waveplot_in = dftb_dir / 'waveplot_in.hsd'
    
    if not detailed_xml.exists():
        raise RuntimeError(f"detailed.xml not found in {dftb_dir}")
    if not eigenvec_bin.exists():
        raise RuntimeError(f"eigenvec.bin not found in {dftb_dir}")
    if not waveplot_in.exists():
        raise RuntimeError(f"waveplot_in.hsd not found in {dftb_dir}")
    
    print(f"\nLoading data from {dftb_dir}")
    
    # Parse detailed.xml
    print("  Parsing detailed.xml...")
    detailed = parse_detailed_xml_custom(detailed_xml)
    coords_bohr = detailed['coords_bohr']
    species_per_atom = detailed['species_per_atom']
    species_names = detailed['species_names']
    natoms = len(coords_bohr)
    nstates = detailed['nstates']
    norb_total = detailed['norb']
    
    # Parse eigenvectors
    print("  Parsing eigenvec.bin...")
    evecs = parse_eigenvec_bin_custom(eigenvec_bin, nstates, norb_total)  # (nStates, nOrb)
    print(f"    nStates={nstates}, nOrb={norb_total}")
    
    # Load real basis from wfc file
    if args.basis_path is None:
        basis_path = REPO_ROOT / 'pyBall/DFTB/data' / f'wfc.{args.basis}.hsd'
    else:
        basis_path = Path(args.basis_path) / f'wfc.{args.basis}.hsd'
    
    print(f"  Loading basis from {basis_path}")
    basis_data = parse_wfc_hsd(str(basis_path))
    basis = convert_wfc_to_species_list_ang(basis_data, resolution_bohr=0.04)
    print(f"  Loaded {len(basis)} species: {[sp['name'] for sp in basis]}")
    
    # Build projector with max_shells=2 (sp only, for old sparse method)
    print("\nSetting up projector (sparse, max_shells=2)...")
    dftb_data_sparse = {
        'coords_bohr': coords_bohr,
        'species_per_atom': species_per_atom,
        'species_names': species_names,
    }
    projector_sparse, atoms_dict = setup_gridprojector_from_dftb( dftb_data_sparse, basis, verbosity=1, max_shells=2 )
    
    # Build projector with max_shells=2 (for dense method on sp system)
    print("Setting up projector (dense, max_shells=2)...")
    projector_dense, _ = setup_gridprojector_from_dftb(  dftb_data_sparse, basis, verbosity=1, max_shells=2 )
    
    # Compute norb_per_atom and orb_offsets for dense method
    sp_by_name = {sp['name']: sp for sp in basis}
    norb_per_atom = []
    max_l = 0
    for ia in range(natoms):
        sp_name = species_names[species_per_atom[ia]]
        sp_info = sp_by_name[sp_name]
        norb = sum(2*orb['l']+1 for orb in sp_info['orbitals'])
        for orb in sp_info['orbitals']:
            max_l = max(max_l, orb['l'])
        norb_per_atom.append(norb)
    norb_per_atom = np.array(norb_per_atom, dtype=np.int32)
    orb_offsets = np.zeros(natoms + 1, dtype=np.int32)
    orb_offsets[1:] = np.cumsum(norb_per_atom)
    
    print(f"  norb_per_atom: {norb_per_atom}")
    print(f"  orb_offsets: {orb_offsets}")
    print(f"  norb_total: {orb_offsets[-1]}")
    
    # Generate 2D grid for visualization (XY plane at fixed Z offset)
    coords_ang = coords_bohr * BOHR2ANG
    rmin = coords_ang.min() - 2.0
    rmax = coords_ang.max() + 2.0
    range_size = rmax - rmin
    ngrid = int(np.ceil(range_size / args.step))
    
    print(f"\nGenerating 2D XY grid (step={args.step} Å, z-offset={args.z_offset} Å)")
    print(f"  Range: [{rmin:.2f}, {rmax:.2f}] Å")
    print(f"  Grid size: {ngrid}x{ngrid}")
    
    # Use library function for grid generation
    points, extent = generate_2d_point_grid(plane='xy', npoints=ngrid, z_offset=args.z_offset, xy_range=(rmin, rmax))
    points = points.astype(np.float32)
    
    # Test orbital projection for a few MOs (HOMO, HOMO-1, LUMO)
    # For H2O, HOMO is around MO 5 (1-based), so index 4 (0-based)
    mo_indices = [3, 4, 5]  # Test a few occupied MOs
    
    print(f"\nTesting orbital projection for MOs {mo_indices}")
    print("-" * 70)
    
    max_orb_diff = 0.0
    orb_results = {}  # Store results for plotting
    for imo in mo_indices:
        if imo >= nstates:
            print(f"  Skipping MO{imo+1} (index out of range)")
            continue
        
        # OLD: sparse method
        coeffs_sparse = evec_to_kernel_coeffs(  evecs[imo], natoms, species_per_atom, species_names, basis  )
        psi_sparse = projector_sparse.project_orbital_points(  points, coeffs_sparse, norb_per_atom, atoms_dict    )
        
        # NEW: dense method
        coeffs_dense = evecs[imo].astype(np.float32)
        psi_dense = projector_dense.project_orbital_dense_points(  points, coeffs_dense, norb_per_atom, orb_offsets, atoms_dict )
        
        # Compare
        diff = np.max(np.abs(psi_sparse - psi_dense))
        max_orb_diff = max(max_orb_diff, diff)
        print(f"  MO{imo+1}: max|diff| = {diff:.2e}")
        
        if diff > args.tol:
            print(f"    WARNING: difference exceeds tolerance {args.tol}")
            # Show some statistics
            print(f"    sparse: mean={psi_sparse.mean():.2e}, std={psi_sparse.std():.2e}")
            print(f"    dense:  mean={psi_dense.mean():.2e}, std={psi_dense.std():.2e}")
        
        # Store for plotting
        orb_results[imo] = {
            'sparse': psi_sparse,
            'dense': psi_dense,
            'diff': psi_sparse - psi_dense
        }
    
    # Test density projection
    print(f"\nTesting density projection")
    print("-" * 70)
    
    # Build density matrix from occupied orbitals
    # For pentacene, assume first 51 orbitals are occupied (102 electrons)
    nocc = 51
    dm_dense = np.zeros((norb_total, norb_total), dtype=np.float32)
    for iocc in range(nocc):
        c = evecs[iocc].astype(np.float32)
        dm_dense += np.outer(c, c)
    
    # OLD: sparse method uses sum of orbitals (no direct density projection in old code)
    # We'll just compute sum of occupied orbitals as reference
    rho_sparse = np.zeros(len(points), dtype=np.float32)
    for iocc in range(nocc):
        coeffs_sparse = evec_to_kernel_coeffs(evecs[iocc], natoms, species_per_atom, species_names, basis )
        psi = projector_sparse.project_orbital_points( points, coeffs_sparse, norb_per_atom, atoms_dict )
        rho_sparse += psi ** 2
    
    # NEW: dense method using density matrix
    rho_dense_dm = projector_dense.project_density_dense_points( points, dm_dense, norb_per_atom, orb_offsets, atoms_dict )
    
    # Also compute sum of orbitals with dense method for comparison
    rho_dense_sum = np.zeros(len(points), dtype=np.float32)
    for iocc in range(nocc):
        coeffs_dense = evecs[iocc].astype(np.float32)
        psi = projector_dense.project_orbital_dense_points(points, coeffs_dense, norb_per_atom, orb_offsets, atoms_dict  )
        rho_dense_sum += psi ** 2
    
    # Compare
    diff_sparse_dense = np.max(np.abs(rho_sparse - rho_dense_sum))
    diff_sum_dm = np.max(np.abs(rho_dense_sum - rho_dense_dm))
    
    print(f"  sparse vs dense sum: max|diff| = {diff_sparse_dense:.2e}")
    print(f"  dense sum vs dense DM: max|diff| = {diff_sum_dm:.2e}")
    
    # Store density results for plotting
    density_results = {
        'sparse': rho_sparse,
        'dense_sum': rho_dense_sum,
        'dense_dm': rho_dense_dm,
        'diff_sparse_dense': rho_sparse - rho_dense_sum,
        'diff_sum_dm': rho_dense_sum - rho_dense_dm
    }
    
    # Generate plots if requested
    if args.plot:
        print("\nGenerating comparison plots...")
        
        # Reshape 1D point data to 2D grid for imshow
        def to_grid(data_1d, ngrid):
            return data_1d.reshape(ngrid, ngrid)
        
        # Prepare data for plot_comparison_2d
        nstates = len(mo_indices)
        sparse_vals = np.zeros((nstates, ngrid, ngrid))
        dense_vals = np.zeros((nstates, ngrid, ngrid))
        diff_vals = np.zeros((nstates, ngrid, ngrid))
        energies = []
        for i, imo in enumerate(mo_indices):
            if imo not in orb_results:
                continue
            res = orb_results[imo]
            sparse_vals[i] = to_grid(res['sparse'], ngrid)
            dense_vals[i] = to_grid(res['dense'], ngrid)
            diff_vals[i] = to_grid(res['diff'], ngrid)
            energies.append(0.0)  # Placeholder - we don't have energies from eigenvec.bin
        
        homo = mo_indices[-1]  # Assume last tested orbital is HOMO
        plane_desc = f'XY plane  z={args.z_offset:.3f} Å'
        title_prefix = f'{args.basis} basis'
        
        # Use library function for orbital comparison
        plot_path = output_dir / 'orbital_comparison.png'
        plot_comparison_2d(sparse_vals, dense_vals, diff_vals, extent, title_prefix, plane_desc, 'orb2points', mo_indices, energies, homo, plot_path, dpi=args.dpi, atom_coords=coords_ang)
        print(f"  Saved: {plot_path}")
        
        # Plot density comparisons using library function
        # Prepare data for plot_comparison_2d (sparse vs dense sum)
        rho_sparse_grid = to_grid(density_results['sparse'], ngrid)
        rho_dense_sum_grid = to_grid(density_results['dense_sum'], ngrid)
        rho_diff_grid = to_grid(density_results['diff_sparse_dense'], ngrid)
        
        # Use plot_comparison_2d for density comparison
        density_sparse_vals = rho_sparse_grid[np.newaxis, :, :]  # Add state dimension
        density_dense_vals = rho_dense_sum_grid[np.newaxis, :, :]
        density_diff_vals = rho_diff_grid[np.newaxis, :, :]
        
        plot_path = output_dir / 'density_comparison.png'
        plot_comparison_2d(density_sparse_vals, density_dense_vals, density_diff_vals, extent, title_prefix, plane_desc, 'density2points', [0], [0.0], 0,   plot_path, dpi=args.dpi, atom_coords=coords_ang)
        print(f"  Saved: {plot_path}")
    
    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"  Max orbital difference: {max_orb_diff:.2e}")
    print(f"  Max density difference (sparse vs dense sum): {diff_sparse_dense:.2e}")
    print(f"  Max density difference (dense sum vs dense DM): {diff_sum_dm:.2e}")
    
    if max_orb_diff < args.tol and diff_sparse_dense < args.tol and diff_sum_dm < args.tol:
        print("\n✓ All tests PASSED within tolerance")
        return 0
    else:
        print(f"\n✗ Some tests FAILED (tolerance = {args.tol})")
        return 1


if __name__ == '__main__':
    sys.exit(main())
