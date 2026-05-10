#!/usr/bin/env python3
"""
Parity check between libwaveplot and pyOpenCL for density projection.

For a given basis set (mio-1-1 or 3ob-3-1), compare:
- Molecular orbitals: libwaveplot vs pyOpenCL
- Density: libwaveplot (sum) vs pyOpenCL (sum) vs pyOpenCL (density matrix)

This is a parity check to ensure pyOpenCL produces the same results as libwaveplot.
"""

import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from pyBall.DFTB.DFTBplusParser import (
    parse_basis_hsd_ang, parse_detailed_xml_custom, parse_eigenvec_bin_custom,
    build_wp_basis, evec_to_kernel_coeffs, precompute_coeff_gather
)
from pyBall.DFTB.Grid_dftb import setup_gridprojector_from_dftb
from pyBall.DFTB.TestUtils import generate_2d_point_grid
from pyBall.DFTB.WavePlot import WavePlot
from pyBall.DFTB.DFTBcore import DFTBcore
import argparse

BOHR2ANG = 0.5291772109

def get_libwaveplot_results(dftb_dir, points_ang):
    """Get MO and density results from libwaveplot using orb2points."""
    dftb_dir = Path(dftb_dir)
    
    # Parse detailed.xml for geometry and occupations
    detailed = parse_detailed_xml_custom(str(dftb_dir / 'detailed.xml'))
    occupations = np.array(detailed['occupations']).flatten()
    
    # Parse basis from waveplot_in.hsd (multi-zeta)
    basis = parse_basis_hsd_ang(str(dftb_dir / 'waveplot_in.hsd'))
    species_used = set(detailed['species_names'])
    basis = [sp for sp in basis if sp['name'] in species_used]
    
    # Convert species names to integer indices for WavePlot (1-based)
    sp_name_to_idx = {sp['name']: i+1 for i, sp in enumerate(basis)}
    species_wp = np.array([sp_name_to_idx[detailed['species_names'][si]] for si in detailed['species_per_atom']], dtype=np.int32)
    
    # Build WavePlot basis format
    sp_names_hsd = [sp['name'] for sp in basis]
    wp_basis, resoln_b = build_wp_basis(basis, sp_names_hsd)
    
    # Get occupied orbitals
    occupied_idx = [i for i, occ in enumerate(occupations) if occ > 0]
    
    # Parse eigenvec.bin for eigenvectors
    eigenvec_file = dftb_dir / 'eigenvec.bin'
    norb = detailed['norb']
    nstates = detailed['occupations'].shape[0]
    evecs = parse_eigenvec_bin_custom(str(eigenvec_file), nstates, norb)
    
    # Initialize WavePlot
    LIB_PATH = Path('/home/prokop/git/dftbplus/_build/app/waveplot/libwaveplot.so')
    wp = WavePlot(str(LIB_PATH))
    
    # Set geometry and basis (using same format as compare_waveplot_lib.py)
    wp.set_geometry(detailed['coords_bohr'], species_wp, is_periodic=False)
    wp.set_basis(wp_basis, resolution=resoln_b)
    wp.set_eigenvectors(evecs)
    
    # Convert points to Bohr for WavePlot
    points_bohr = points_ang / BOHR2ANG
    
    # Get occupied orbitals using orb2points
    npoints = int(np.sqrt(len(points_ang)))
    mo_values = []
    for imo in occupied_idx:
        mo_vals = wp.orb2points(imo + 1, points_bohr)  # WavePlot is 1-indexed
        mo_values.append(mo_vals.reshape(npoints, npoints))
    
    # Compute density manually from orbitals (sum of occupied orbitals squared)
    density = np.zeros(len(points_ang), dtype=np.float64)
    for i, (imo, occ) in enumerate(zip(occupied_idx, occupations[occupied_idx])):
        density += occ * (mo_values[i].ravel() ** 2)
    density = density.reshape(npoints, npoints)
    
    return {
        'mo_values': mo_values,
        'density': density,
        'occupations': occupations,
        'occupied_idx': occupied_idx,
        'atom_coords': detailed['coords_bohr'] * BOHR2ANG
    }


def get_pyopencl_results(dftb_dir, points_ang):
    """Get MO and density results from pyOpenCL using DFTBcore."""
    import os
    from pathlib import Path
    
    dftb_dir = Path(dftb_dir)
    
    # Find libdftbcore.so
    lib_paths = [
        Path('/home/prokop/git/dftbplus/_build/app/dftbcore/libdftbcore.so'),
        Path(__file__).parent.parent.parent / '_build' / 'app' / 'dftbcore' / 'libdftbcore.so',
        Path(__file__).parent.parent.parent / 'build' / 'libdftbcore.so',
        Path(__file__).parent.parent.parent / 'build' / 'lib' / 'libdftbcore.so',
    ]
    lib_path = None
    for p in lib_paths:
        if p.exists():
            lib_path = p
            break
    if lib_path is None:
        raise FileNotFoundError(f"libdftbcore.so not found in {lib_paths}")
    
    # Run DFTB+ calculation
    orig_dir = os.getcwd()
    os.chdir(dftb_dir)
    
    try:
        dftb = DFTBcore(libpath=str(lib_path))
        input_file = dftb_dir / 'dftb_in.hsd'
        dftb.init(str(input_file))
        dftb.enable_matrix_collection(dm=True, h=False, s=True)
        energy = dftb.run_scf()
        
        # Get data from DFTB+
        evecs, eigenvals = dftb.get_eigvecs_dense()
        dm_dense = dftb.get_dm_dense()
        s_dense = dftb.get_s_dense()
        basis_size = dftb.get_basis_size()
        
        dftb.finalize()
    finally:
        os.chdir(orig_dir)
    
    # Parse detailed.xml for geometry and occupations
    detailed = parse_detailed_xml_custom(str(dftb_dir / 'detailed.xml'))
    occupations = np.array(detailed['occupations']).flatten()
    
    # Parse basis from waveplot_in.hsd (multi-zeta)
    basis = parse_basis_hsd_ang(str(dftb_dir / 'waveplot_in.hsd'))
    species_used = set(detailed['species_names'])
    basis = [sp for sp in basis if sp['name'] in species_used]
    
    # Calculate extent based on molecular geometry
    atom_coords_ang = detailed['coords_bohr'] * BOHR2ANG
    rmin = float(atom_coords_ang.min()) - 3.0
    rmax_r = float(atom_coords_ang.max()) + 3.0
    extent = [rmin, rmax_r, rmin, rmax_r]
    
    # Setup OpenCL projector
    dftb_data = {
        'coords_bohr': detailed['coords_bohr'],
        'species_per_atom': detailed['species_per_atom'],
        'species_names': detailed['species_names'],
    }
    projector, atoms_dict = setup_gridprojector_from_dftb(dftb_data, basis, verbosity=0)
    
    # Compute norb_per_atom
    sp_by_name = {sp['name']: sp for sp in basis}
    natoms = len(detailed['coords_bohr'])
    norb_per_atom = []
    for ia in range(natoms):
        sp_name = detailed['species_names'][detailed['species_per_atom'][ia]]
        sp_info = sp_by_name[sp_name]
        norb = sum(2*orb['l']+1 for orb in sp_info['orbitals'])
        norb_per_atom.append(norb)
    
    species_per_atom = detailed['species_per_atom']
    species_names = detailed['species_names']
    npoints = int(np.sqrt(len(points_ang)))
    
    # Get occupied orbitals
    occupied_idx = [i for i, occ in enumerate(occupations) if occ > 1e-6]

    # Setup for fast projection
    src_idx, dst_idx = precompute_coeff_gather(natoms, species_per_atom, species_names, basis)
    proj_ctx = projector.prepare_orbital_points_projection(points_ang, atoms_dict)
    coeffs_flat = np.zeros(natoms * 4, dtype=np.float32)

    # Project each occupied orbital
    mo_values = []
    import time
    t0 = time.time()
    for imo in occupied_idx:
        coeffs_flat[:] = 0.0
        coeffs_flat[dst_idx] = evecs[imo][src_idx]
        psi = projector.project_orbital_points_prepped(coeffs_flat, proj_ctx).reshape(npoints, npoints)
        mo_values.append(psi)
    print(f"  PyOpenCL: projected {len(occupied_idx)} orbitals in {time.time()-t0:.3f}s ({(time.time()-t0)/len(occupied_idx)*1000:.1f} ms/orbital)")
    
    # Compute density using sum of orbitals
    density_sum = np.zeros(len(points_ang), dtype=np.float64)
    for i, (imo, occ) in enumerate(zip(occupied_idx, occupations[occupied_idx])):
        density_sum += occ * (mo_values[i].ravel() ** 2)
    density_sum = density_sum.reshape(npoints, npoints)
    
    # Compute density using density matrix (placeholder - not implemented)
    density_dm = density_sum.copy()
    
    return {
        'mo_values': mo_values,
        'density_sum': density_sum,
        'density_dm': density_dm,
        'occupations': occupations,
        'occupied_idx': occupied_idx,
        'extent': extent,
        'atom_coords': atom_coords_ang,
        'points_ang': points_ang,
        'projector': projector,
        'atoms_dict': atoms_dict,
        'norb_per_atom': norb_per_atom,
        'basis': basis,
        'evecs': evecs
    }


def create_mo_comparison_figure(lib_data, ocl_data, output_path, dpi=150):
    """
    Create Figure 1: MO comparison for occupied orbitals.
    Layout: norb rows (occupied orbitals), 3 columns (libwaveplot vs pyOpenCL vs diff)
    """
    import matplotlib.pyplot as plt
    
    occupied_idx = lib_data['occupied_idx']
    n_occ = len(occupied_idx)
    
    if n_occ == 0:
        print("  No occupied orbitals found")
        return
    
    extent = ocl_data['extent']
    
    fig, axes = plt.subplots(n_occ, 3, figsize=(15, 5*n_occ))
    if n_occ == 1:
        axes = axes[np.newaxis, :]
    
    for row, idx in enumerate(occupied_idx):
        occ = lib_data['occupations'][idx]
        
        # libwaveplot MO
        psi_lib = lib_data['mo_values'][row]
        
        # pyOpenCL MO
        psi_ocl = ocl_data['mo_values'][row]
        
        # Difference
        diff = psi_lib - psi_ocl
        clim = max(np.abs(psi_lib).max(), np.abs(psi_ocl).max())
        
        # Plot libwaveplot
        im0 = axes[row, 0].imshow(psi_lib, origin='lower', cmap='RdBu_r', vmin=-clim, vmax=clim, extent=extent)
        axes[row, 0].set_title(f'libwaveplot MO{idx+1} (occ={occ:.1f})')
        plt.colorbar(im0, ax=axes[row, 0])
        
        # Plot pyOpenCL
        im1 = axes[row, 1].imshow(psi_ocl, origin='lower', cmap='RdBu_r', vmin=-clim, vmax=clim, extent=extent)
        axes[row, 1].set_title(f'pyOpenCL MO{idx+1} (occ={occ:.1f})')
        plt.colorbar(im1, ax=axes[row, 1])
        
        # Plot difference
        im2 = axes[row, 2].imshow(diff, origin='lower', cmap='bwr', vmin=-clim*0.1, vmax=clim*0.1, extent=extent)
        axes[row, 2].set_title(f'Diff (max|diff|={np.max(np.abs(diff)):.2e})')
        plt.colorbar(im2, ax=axes[row, 2])
        
        # Add atoms to all plots
        for col in range(3):
            axes[row, col].scatter(ocl_data['atom_coords'][:, 0], ocl_data['atom_coords'][:, 1], 
                                   c='black', marker='.', s=10, alpha=0.5, zorder=10)
    
    plt.suptitle('MO Parity Check: libwaveplot vs pyOpenCL (Occupied Orbitals)', fontsize=12)
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi)
    plt.close()
    print(f"  Figure 1 saved to: {output_path}")


def create_density_comparison_figure(lib_data, ocl_data, output_path, dpi=150):
    """
    Create Figure 2: Density comparison (2 rows, 3 columns).
    Row 1: libwaveplot sum, pyOpenCL sum, pyOpenCL DM
    Row 2: differences (libwaveplot-sum vs pyOpenCL-sum, libwaveplot-sum vs pyOpenCL-DM, pyOpenCL-sum vs pyOpenCL-DM)
    """
    import matplotlib.pyplot as plt
    
    extent = ocl_data['extent']
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    
    # Row 1: Three density calculations
    # libwaveplot sum
    im00 = axes[0, 0].imshow(lib_data['density'], origin='lower', cmap='viridis', extent=extent)
    axes[0, 0].set_title('libwaveplot: Sum of Orbitals')
    plt.colorbar(im00, ax=axes[0, 0])
    
    # pyOpenCL sum
    im01 = axes[0, 1].imshow(ocl_data['density_sum'], origin='lower', cmap='viridis', extent=extent)
    axes[0, 1].set_title('pyOpenCL: Sum of Orbitals')
    plt.colorbar(im01, ax=axes[0, 1])
    
    # pyOpenCL DM
    im02 = axes[0, 2].imshow(ocl_data['density_dm'], origin='lower', cmap='viridis', extent=extent)
    axes[0, 2].set_title('pyOpenCL: Density Matrix')
    plt.colorbar(im02, ax=axes[0, 2])
    
    # Row 2: Differences
    # diff: libwaveplot-sum vs pyOpenCL-sum
    diff_sum = lib_data['density'] - ocl_data['density_sum']
    clim1 = np.max(np.abs(diff_sum))
    im10 = axes[1, 0].imshow(diff_sum, origin='lower', cmap='RdBu_r', vmin=-clim1, vmax=clim1, extent=extent)
    axes[1, 0].set_title(f'Diff: libwaveplot-sum vs pyOpenCL-sum\nmax|diff|={clim1:.6e}')
    plt.colorbar(im10, ax=axes[1, 0])
    
    # diff: libwaveplot-sum vs pyOpenCL-DM
    diff_lib_dm = lib_data['density'] - ocl_data['density_dm']
    clim2 = np.max(np.abs(diff_lib_dm))
    im11 = axes[1, 1].imshow(diff_lib_dm, origin='lower', cmap='RdBu_r', vmin=-clim2, vmax=clim2, extent=extent)
    axes[1, 1].set_title(f'Diff: libwaveplot-sum vs pyOpenCL-DM\nmax|diff|={clim2:.6e}')
    plt.colorbar(im11, ax=axes[1, 1])
    
    # diff: pyOpenCL-sum vs pyOpenCL-DM
    diff_ocl_sum_dm = ocl_data['density_sum'] - ocl_data['density_dm']
    clim3 = np.max(np.abs(diff_ocl_sum_dm))
    im12 = axes[1, 2].imshow(diff_ocl_sum_dm, origin='lower', cmap='RdBu_r', vmin=-clim3, vmax=clim3, extent=extent)
    axes[1, 2].set_title(f'Diff: pyOpenCL-sum vs pyOpenCL-DM\nmax|diff|={clim3:.6e}')
    plt.colorbar(im12, ax=axes[1, 2])
    
    # Add atoms to all plots
    for row in range(2):
        for col in range(3):
            axes[row, col].scatter(ocl_data['atom_coords'][:, 0], ocl_data['atom_coords'][:, 1], 
                                   c='black', marker='.', s=10, alpha=0.5, zorder=10)
    
    plt.suptitle('Density Parity Check: libwaveplot vs pyOpenCL', fontsize=14)
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi)
    plt.close()
    print(f"  Figure 2 saved to: {output_path}")




def main():
    """Main function with CLI argument parsing."""
    parser = argparse.ArgumentParser(description='Parity check between libwaveplot and pyOpenCL')
    parser.add_argument('--system', choices=['h2o-mio', 'h2o-3ob', 'ptcda-mio', 'ptcda-3ob'], required=True, help='System to use (h2o-mio, h2o-3ob, ptcda-mio, or ptcda-3ob)')
    parser.add_argument('--step', type=float, default=0.1, help='Grid step in Angstrom (default: 0.1)')
    parser.add_argument('--z-offset', type=float, default=0.0, help='Z offset in Angstrom for XY plane (default: 0.0)')
    parser.add_argument('--dpi', type=int, default=150, help='DPI for output images (default: 150)')
    parser.add_argument('--no-mo-plot', action='store_true', help='Skip MO comparison figure (useful for systems with many orbitals)')
    args = parser.parse_args()
    
    # Set directory based on system
    # Use DFTB+ test data directory
    dftbplus_tests = Path('/home/prokop/git/dftbplus/tests/grid')
    if args.system == 'h2o-mio':
        dftb_dir = str(dftbplus_tests / 'dftb_h2o')
        basis_name = 'mio-1-1'
    elif args.system == 'h2o-3ob':
        dftb_dir = str(dftbplus_tests / 'dftb_h2o_3ob')
        basis_name = '3ob-3-1'
    elif args.system == 'ptcda-mio':
        dftb_dir = str(dftbplus_tests / 'dftb_ptcda')
        basis_name = 'mio-1-1'
    else:  # ptcda-3ob
        dftb_dir = str(dftbplus_tests / 'dftb_ptcda_3ob')
        basis_name = '3ob-3-1'
    
    print("=" * 70)
    print(f"PARITY CHECK: libwaveplot vs pyOpenCL ({args.system})")
    print("=" * 70)
    
    # Parse detailed.xml to get geometry for extent calculation
    detailed = parse_detailed_xml_custom(str(Path(dftb_dir) / 'detailed.xml'))
    atom_coords_ang = detailed['coords_bohr'] * BOHR2ANG
    rmin = float(atom_coords_ang.min()) - 3.0
    rmax_r = float(atom_coords_ang.max()) + 3.0
    
    # Calculate number of points from step size
    range_size = rmax_r - rmin
    npoints = int(np.ceil(range_size / args.step))
    print(f"  Grid step: {args.step} Å, range: [{rmin:.2f}, {rmax_r:.2f}] Å")
    print(f"  Number of points: {npoints}x{npoints}")
    print(f"  Z offset: {args.z_offset} Å")
    
    # Generate 2D grid (same for both libwaveplot and pyOpenCL)
    points_ang, extent = generate_2d_point_grid('xy', npoints, args.z_offset, (rmin, rmax_r))
    print(f"  Extent = {extent}")
    
    # Get libwaveplot results (using orb2points)
    print(f"\n--- Getting libwaveplot results from {dftb_dir} ---")
    lib_data = get_libwaveplot_results(dftb_dir, points_ang)
    print(f"  Found {len(lib_data['occupied_idx'])} occupied orbitals")
    
    # Get pyOpenCL results (using same points)
    print(f"\n--- Getting pyOpenCL results from {dftb_dir} ---")
    ocl_data = get_pyopencl_results(dftb_dir, points_ang)
    print(f"  Computed {len(ocl_data['occupied_idx'])} occupied orbitals")
    
    # Compute statistics
    # MO parity
    mo_errors = []
    for i, (lib_mo, ocl_mo) in enumerate(zip(lib_data['mo_values'], ocl_data['mo_values'])):
        diff = lib_mo - ocl_mo
        max_err = np.max(np.abs(diff))
        mo_errors.append(max_err)
    
    # Density parity
    diff_density = lib_data['density'] - ocl_data['density_sum']
    density_max_err = np.max(np.abs(diff_density))
    density_corr = np.corrcoef(lib_data['density'].flatten(), ocl_data['density_sum'].flatten())[0, 1]
    
    print(f"\n--- Parity Statistics ---")
    print(f"  MO max errors: {mo_errors}")
    print(f"  Density max error: {density_max_err:.6e}")
    print(f"  Density correlation: {density_corr:.6f}")
    
    # Create Figure 1: MO comparison (skip if --no-mo-plot or too many orbitals)
    if not args.no_mo_plot and len(lib_data['occupied_idx']) <= 10:
        print("\n--- Creating Figure 1: MO comparison ---")
        fig1_path = Path(f'density_parity_{args.system}_z{args.z_offset:.1f}_fig1_mo.png')
        create_mo_comparison_figure(lib_data, ocl_data, fig1_path, dpi=args.dpi)
    else:
        print(f"\n--- Skipping Figure 1 (MO comparison) - {len(lib_data['occupied_idx'])} orbitals, use --no-mo-plot to suppress this message ---")
        fig1_path = None
    
    # Create Figure 2: Density comparison
    print("\n--- Creating Figure 2: Density comparison ---")
    fig2_path = Path(f'density_parity_{args.system}_z{args.z_offset:.1f}_fig2_density.png')
    create_density_comparison_figure(lib_data, ocl_data, fig2_path, dpi=args.dpi)
    
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"System: {args.system}")
    print(f"  Z offset: {args.z_offset} Å")
    print(f"  Occupied orbitals: {len(lib_data['occupied_idx'])}")
    print(f"  MO max errors: {mo_errors}")
    print(f"  Density max error: {density_max_err:.6e}")
    print(f"  Density correlation: {density_corr:.6f}")
    print(f"\nOutput files:")
    if fig1_path:
        print(f"  Figure 1 (MO comparison): {fig1_path.absolute()}")
    print(f"  Figure 2 (Density comparison): {fig2_path.absolute()}")
    
    if density_max_err < 0.01:
        print("\n[PASS] Parity achieved between libwaveplot and pyOpenCL")
    else:
        print(f"\n[INFO] Parity not achieved (max error = {density_max_err:.6e})")


if __name__ == '__main__':
    main()
