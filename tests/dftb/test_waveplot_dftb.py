#!/usr/bin/env python3
"""
DFTB+ OpenCL Waveplot Test Script

Tests orbital projection and density computation using OpenCL.
Supports H2O (cluster) and PTCDA (periodic) test systems.

Usage:
    python test_waveplot_dftb.py --system H2O --step 0.2
    python test_waveplot_dftb.py --system PTCDA --z-height 2.0 --step 0.15
"""

import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Add repo root to path
_REPO_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(_REPO_ROOT))

from pyBall.DFTB.DFTBplusParser import (
    DFTBplusParser, compute_sto_radial, parse_basis_hsd_ang,
    parse_detailed_xml_custom, parse_eigenvec_bin_custom,
    build_wp_basis, evec_to_kernel_coeffs
)
from pyBall.DFTB.Grid_dftb import GridProjector, setup_gridprojector_from_dftb
from pyBall import plotUtils as pu
from pyBall.DFTB.TestUtils import print_eigenvecs


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='DFTB+ OpenCL Waveplot Test',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    # DFTB+ directory (instead of system selection)
    parser.add_argument('--dftb-dir', type=str, default=None, help='DFTB+ run directory containing detailed.xml, eigenvec.bin, waveplot_in.hsd')
    
    # Grid parameters
    parser.add_argument('--step', type=float, default=0.15,  help='Grid spacing in Angstrom')
    parser.add_argument('--margin', type=float, default=3.0,  help='Margin around molecule in Angstrom')
    parser.add_argument('--z-height', type=float, default=None,  help='Height above molecular plane for 2D slice (PTCDA only)')
    parser.add_argument('--ngrid', type=int, nargs=3, default=None, metavar=('NX', 'NY', 'NZ'),  help='Explicit grid dimensions (overrides auto)')
    
    # Point evaluation (Method 2)
    parser.add_argument('--points', action='store_true',  help='Method 2: evaluate at points along z-axis instead of 3D grid')
    parser.add_argument('--z-range', type=float, nargs=2, default=[-3.0, 3.0],  metavar=('ZMIN', 'ZMAX'), help='Z-range for point evaluation (Å)')
    parser.add_argument('--npoints', type=int, default=301, help='Number of points for evaluation')
    
    # Output options
    parser.add_argument('--output-dir', type=str, default='waveplot_output',  help='Output directory for plots')
    parser.add_argument('--prefix', type=str, default='dftb_waveplot',  help='Filename prefix for output files')
    parser.add_argument('--dpi', type=int, default=150,  help='DPI for saved figures')
    parser.add_argument('--no-show', action='store_true',  help='Do not display plots (only save)')
    
    # Plotting options
    parser.add_argument('--plot-mos', type=int, nargs='+', default=None,  help='MO indices to plot (default: HOMO-2 to LUMO+2)')
    parser.add_argument('--plot-density', action='store_true',  help='Plot total electron density')
    parser.add_argument('--plot-slices', action='store_true',  help='Plot 3 orthogonal slices through 3D grid')
    parser.add_argument('--cmap-orbital', type=str, default='RdBu_r',  help='Colormap for orbital plots (signed)')
    parser.add_argument('--cmap-density', type=str, default='hot',  help='Colormap for density plots (positive)')
    
    # OpenCL options
    parser.add_argument('--nmax-atom', type=int, default=64,  help='Max atoms per task block')
    parser.add_argument('--gpu-tasks', action='store_true',  help='Use GPU-based task builder')
    
    # Validation options
    parser.add_argument('--vmin-min', type=float, default=1e-6,  help='Minimum expected absolute value')
    parser.add_argument('--vmax-max', type=float, default=1e6,  help='Maximum expected absolute value')
    parser.add_argument('--check-nan', action='store_true', default=True,  help='Check for NaN values')
    
    # Verbosity
    parser.add_argument('-v', '--verbose', action='count', default=0,  help='Increase verbosity (use -v, -vv, or -vvv)')
    parser.add_argument('--debug', action='store_true',  help='Enable debug mode')
    parser.add_argument('--print-eigenvec', action='store_true', help='Print eigenvectors from eigenvec.bin and exit')
    
    return parser.parse_args()


def main():
    """Main entry point - generic for any DFTB+ directory."""
    args = parse_args()
    
    # Determine DFTB+ directory
    dftb_dir = Path(args.dftb_dir) if args.dftb_dir else Path(__file__).parent / 'dftb_h2o'
    assert dftb_dir.exists(), f"DFTB+ directory not found: {dftb_dir}"
    system_name = dftb_dir.name
    
    # Print eigenvectors if requested
    if args.print_eigenvec:
        print_eigenvecs(dftb_dir / 'eigenvec.bin', dftb_dir / 'detailed.xml', dftb_dir / 'waveplot_in.hsd', max_orbitals=6)
        return
    output_dir = Path(args.output_dir) / system_name
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("DFTB+ OpenCL Waveplot Test")
    print("=" * 60)
    print(f"Directory: {dftb_dir}")
    print(f"Step: {args.step} Å")
    print(f"Margin: {args.margin} Å")
    
    # Parse DFTB+ data using shared functions
    geo = parse_detailed_xml_custom(dftb_dir / 'detailed.xml')
    atom_coords_b = geo['coords_bohr']
    species_names = geo['species_names']
    species_per_atom = geo['species_per_atom']
    natoms = len(species_per_atom)
    nstates = geo['nstates']
    norb = geo['norb']
    occupations = geo['occupations']
    
    # Convert to Angstrom for OpenCL
    BOHR2ANG = 0.5291772109
    coords = atom_coords_b * BOHR2ANG
    atom_symbols = [species_names[si] for si in species_per_atom]
    
    print(f"  natoms={natoms}, nstates={nstates}, norb={norb}")
    print(f"  species: {set(atom_symbols)}")
    
    # Parse eigenvectors
    evecs_full = parse_eigenvec_bin_custom(dftb_dir / 'eigenvec.bin', nstates, norb)
    
    # Parse STO basis
    hsd_path = dftb_dir / 'waveplot_in.hsd'
    species_list_sto = parse_basis_hsd_ang(hsd_path)
    
    # Setup projector using shared function
    dftb_data_ocl = {
        'coords_bohr': atom_coords_b,
        'species_per_atom': species_per_atom,
        'species_names': species_names
    }
    projector, atoms_dict = setup_gridprojector_from_dftb(dftb_data_ocl, species_list_sto, verbosity=args.verbose)
    
    # Compute norb_per_atom
    sp_by_name = {sp['name']: sp for sp in species_list_sto}
    norb_per = np.array([sum(2*o['l']+1 for o in sp_by_name[species_names[si]]['orbitals'])
                          for si in species_per_atom], dtype=np.int32)
    numorb_max = norb_per.max()
    
    # Read energies from band.out if available
    energies_ev = np.zeros(nstates)
    band_path = dftb_dir / 'band.out'
    if band_path.exists():
        with open(band_path) as f:
            for line in f:
                parts = line.split()
                if len(parts) == 3:
                    try:
                        idx = int(parts[0])
                        if 1 <= idx <= nstates:
                            energies_ev[idx-1] = float(parts[1])
                    except ValueError:
                        pass
    
    # MO selection
    if args.plot_mos:
        mo_indices = [i-1 for i in args.plot_mos]
    else:
        # Default: HOMO-2 to LUMO+2
        homo_idx = np.where(occupations > 0.5)[0][-1] if len(occupations) > 0 else nstates // 2
        mo_start = max(0, homo_idx - 2)
        mo_end = min(nstates - 1, homo_idx + 2)
        mo_indices = list(range(mo_start, mo_end + 1))
    
    # Method 2: Point evaluation
    if args.points:
        print(f"\n[Method 2: Point evaluation along z-axis]")
        z_vals = np.linspace(args.z_range[0], args.z_range[1], args.npoints)
        center = coords.mean(axis=0)
        points_ang = np.column_stack([np.full_like(z_vals, center[0]), np.full_like(z_vals, center[1]),z_vals])
        print(f"  Evaluating at {args.npoints} points along z-axis")
        
        ocl_vals_list = []
        for imo in mo_indices:
            coeffs_k = evec_to_kernel_coeffs(evecs_full[imo], natoms, species_per_atom, species_names, species_list_sto)
            psi = projector.project_orbital_points( points_ang.astype(np.float32), coeffs_k, norb_per, atoms_dict  )
            ocl_vals_list.append(psi.astype(np.float64))
            print(f"  MO{imo+1} max = {np.abs(psi).max():.6f}")
        
        # Plot
        fig, axes = plt.subplots(len(mo_indices), 2, figsize=(14, 3*len(mo_indices)))
        if len(mo_indices) == 1:
            axes = np.array([axes])
        axes = axes.flatten()
        
        for idx, imo in enumerate(mo_indices):
            ax_lin = axes[2*idx]
            ax_log = axes[2*idx+1]
            
            ax_lin.plot(z_vals, ocl_vals_list[idx], 'r-', lw=2, label='OpenCL')
            ax_lin.axhline(0, c='gray', lw=0.5)
            ax_lin.set_xlabel('z (Å)')
            ax_lin.set_ylabel('ψ')
            ax_lin.set_title(f'MO{imo+1} E={energies_ev[imo]:.2f}eV')
            ax_lin.legend(fontsize=8)
            
            mask = z_vals >= 0
            eps = 1e-12
            y_oc = np.abs(ocl_vals_list[idx][mask]).clip(eps)
            ax_log.semilogy(z_vals[mask], y_oc, 'r-', lw=2, label='OpenCL')
            ax_log.set_xlabel('z (Å)')
            ax_log.set_ylabel('|ψ| (log)')
            ax_log.set_title(f'MO{imo+1} — log scale')
            ax_log.legend(fontsize=8)
        
        fig.tight_layout()
        out_file = output_dir / f'points_n{args.npoints}.png'
        fig.savefig(str(out_file), dpi=args.dpi)
        print(f"\nSaved: {out_file}")
        if not args.no_show:
            plt.show()
        return 0
    
    # Method 1: 3D grid projection
    print(f"\n[Method 1: 3D grid projection]")
    
    # Auto grid setup
    pos_min = coords.min(axis=0) - args.margin
    pos_max = coords.max(axis=0) + args.margin
    if args.ngrid:
        ngrid = np.array(args.ngrid, dtype=np.int32)
    else:
        span = pos_max - pos_min
        ngrid = np.ceil(span / args.step).astype(int)
        # Round to block size
        block = 8
        ngrid = ((ngrid + block - 1) // block) * block
    
    origin = pos_min
    dA = np.array([args.step, 0.0, 0.0])
    dB = np.array([0.0, args.step, 0.0])
    dC = np.array([0.0, 0.0, args.step])
    grid_spec = {'origin': origin, 'dA': dA, 'dB': dB, 'dC': dC, 'ngrid': ngrid}
    
    print(f"  Grid: {ngrid}, origin: {origin}, step: {args.step}")
    
    # Project MOs
    print(f"\n[Projecting {len(mo_indices)} MOs]")
    grids = []
    for imo in mo_indices:
        coeffs_k = evec_to_kernel_coeffs(evecs_full[imo], natoms, species_per_atom,  species_names, species_list_sto)
        grid_3d = projector.project_orbital(coeffs_k, norb_per, atoms_dict, grid_spec,  nMaxAtom=args.nmax_atom)
        grids.append(grid_3d)
        gmax = np.abs(grid_3d).max()
        print(f"  MO{imo+1} E={energies_ev[imo]:.2f}eV  |ψ|max={gmax:.4e}")
    
    # Plot combined figure
    n_cols = min(3, len(mo_indices))
    n_rows = (len(mo_indices) + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows))
    if len(mo_indices) == 1:
        axes = np.array([axes])
    axes = np.array(axes).flatten()
    
    for idx, imo in enumerate(mo_indices):
        grid_3d = grids[idx]
        # Extract slice at specified z-height (or middle if not specified)
        if args.z_height is not None:
            iz = int((args.z_height - origin[2]) / args.step)
            iz = max(0, min(iz, grid_3d.shape[2] - 1))
        else:
            iz = grid_3d.shape[2] // 2
        grid_2d = grid_3d[:, :, iz]
        z_val = origin[2] + iz * args.step
        
        extent = [origin[0], origin[0] + ngrid[0]*args.step,origin[1], origin[1] + ngrid[1]*args.step]
        
        ax = axes[idx]
        clim = max(abs(grid_2d.min()), abs(grid_2d.max())) or 1e-10
        im = ax.imshow(grid_2d.T, origin='lower', extent=extent,   cmap='RdBu_r', vmin=-clim, vmax=clim)
        pu.plotAtoms(apos=coords, es=None, sizes=100, colors='gray',marker='o', axes=(0, 1))
        ax.set_title(f"MO{imo+1} E={energies_ev[imo]:.2f}eV\nz={z_val:.2f}Å", fontsize=8)
        plt.colorbar(im, ax=ax, fraction=0.046)
    
    for ax in axes[len(mo_indices):]:
        ax.set_visible(False)
    
    fig.suptitle(f"{system_name} Molecular Orbitals", fontsize=12)
    fig.tight_layout()
    combined_path = output_dir / f"{args.prefix}_all_MOs.png"
    fig.savefig(combined_path, dpi=args.dpi)
    print(f"\nSaved: {combined_path}")
    if not args.no_show:
        plt.show()
    
    print(f"\n[Outputs saved to: {output_dir}]")
    return 0


if __name__ == '__main__':
    sys.exit(main())
