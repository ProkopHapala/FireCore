#!/usr/bin/env python3
"""
Compare OpenCL vs Fortran orbital projection for molecules.

Usage:
    python test_h2o_orbital_comparison.py [options]
    
Options:
    -m, --mol       Molecule name (default: H2O)
    -o, --orbitals  MO indices to plot (default: all occupied), e.g. "0,1,2,3,4" or "1-5"
    -z, --zplane    z-height for imaging plane in Å (default: 1.0)
    -n, --grid      Grid resolution as "nx,ny,nz" (default: "50,50,20")
    -e, --extent    Grid extent as "Lx,Ly" in Å (default: "8.0,8.0")
    -v, --verbose   Verbose output
    -h, --help      Show this help

Examples:
    # Plot all occupied orbitals of H2O at z=1.0Å
    python test_h2o_orbital_comparison.py
    
    # Plot only HOMO (MO 5) and HOMO-1 (MO 4) of H2O at z=0.5Å
    python test_h2o_orbital_comparison.py -o "3,4" -z 0.5
    
    # Plot CH4 orbitals with custom grid
    python test_h2o_orbital_comparison.py -m CH4 -o "0-7" -n "60,60,30" -z 1.5
"""

import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.path.dirname(__file__), "..", ".."))

from pyBall import FireCore as fc
from pyBall.atomicUtils import load_xyz
from pyBall.FireballOCL import Grid as ocl_grid
from pyBall.FireballOCL.STM_utils import set_export_dir

def plot_orbital_comparison(ax, data_xy, atomPos, atomTypes, extent, title, coeffs_text, enames):
    """Plot orbital with atoms overlaid and symmetric color range."""
    vmax = np.max(np.abs(data_xy))
    if vmax < 1e-10:
        vmax = 1e-3
    vmin = -vmax
    
    im = ax.imshow(data_xy.T, origin='lower', extent=extent, cmap='bwr', aspect='equal', vmin=vmin, vmax=vmax)
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    
    # Overlay atoms
    ax.scatter(atomPos[:, 0], atomPos[:, 1], c='green', s=100,  edgecolors='black', linewidths=1, zorder=10, marker='o')
    
    # Add atom labels
    for i, (pos, name) in enumerate(zip(atomPos, enames)):
        ax.annotate(f'{name}{i}', (pos[0], pos[1]), 
                   textcoords="offset points", xytext=(5, 5), fontsize=8)
    
    # Add coefficients text
    ax.text(0.02, 0.98, coeffs_text, transform=ax.transAxes, fontsize=8,
            verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    ax.set_xlabel('x (Å)')
    ax.set_ylabel('y (Å)')
    ax.set_title(title, fontsize=10)

def parse_orbital_indices(arg_str, max_mo):
    """Parse orbital indices from string like '0,1,2,3,4' or '1-5' or 'all'"""
    if arg_str.lower() == 'all':
        return list(range(max_mo))
    
    indices = []
    parts = arg_str.split(',')
    for part in parts:
        if '-' in part:
            start, end = part.split('-')
            indices.extend(range(int(start), int(end) + 1))
        else:
            indices.append(int(part))
    return sorted(set(indices))  # Remove duplicates and sort


def main():
    # Parse CLI arguments
    parser = argparse.ArgumentParser(
        description='Compare OpenCL vs Fortran orbital projection',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python test_h2o_orbital_comparison.py -m H2O -o 0-4 -z 1.0
  python test_h2o_orbital_comparison.py -m CH4 -o 0,1,2,3,4 -z 1.5
        """)
    
    parser.add_argument('-m', '--mol',      type=str,   default='H2O',  help='Molecule name (default: H2O)')
    parser.add_argument('-o', '--orbitals', type=str,   default='0-5', help='MO indices: "0,1,2" or "0-5" or "all" (default: 0-5 = all MOs)')
    parser.add_argument('-z', '--zplane',   type=float, default=1.0, help='z-height for imaging plane in Å (default: 1.0)')
    parser.add_argument('-n', '--grid',     type=str,   default='50,50,20', help='Grid resolution as "nx,ny,nz" (default: 50,50,20)')
    parser.add_argument('-e', '--extent',   type=str,   default='8.0,8.0',  help='Grid extent as "Lx,Ly" in Å (default: 8.0,8.0)')
    parser.add_argument('-v', '--verbose', action='store_true',  help='Verbose output')
    parser.add_argument('--debug-atom', type=int, default=None, help='Debug mode: atom index for single-basis test')
    parser.add_argument('--debug-orb', type=str, default=None, help='Debug mode: orbital type (s, px, py, pz)')
    parser.add_argument('--debug-val', type=float, default=1.0, help='Debug mode: coefficient value (default: 1.0)')
    args = parser.parse_args()
    
    # Parse grid parameters
    nx, ny, nz = [int(x) for x in args.grid.split(',')]
    Lx, Ly = [float(x) for x in args.extent.split(',')]
    z_height = args.zplane
    verbose = args.verbose
    
    mol_name = args.mol
    export_dir = set_export_dir(f"export/{mol_name.lower()}_orbital_comparison")
    
    # Load molecule
    xyz_path = f"../../cpp/common_resources/xyz/{mol_name}.xyz"
    if not os.path.exists(xyz_path):
        raise FileNotFoundError(f"XYZ file not found: {xyz_path}")
    
    atomPos, atomTypes, enames, qs, comment = load_xyz(xyz_path)
    natoms = len(atomPos)
    
    print(f"{mol_name}: {natoms} atoms")
    print(f"Atom positions:")
    for i, (pos, name) in enumerate(zip(atomPos, enames)):
        print(f"  {i}: {name} at ({pos[0]:.3f}, {pos[1]:.3f}, {pos[2]:.3f})")
    
    # Initialize Fireball
    print("\nInitializing Fireball...")
    fc.setVerbosity(0)
    fc.initialize(atomType=atomTypes, atomPos=atomPos, verbosity=0)
    fc.SCF(positions=atomPos, nmax_scf=200)
    
    # Get orbital info
    dims = fc.get_HS_dims()
    sparse = fc.get_HS_neighs(dims)
    norb_total = dims.norbitals
    
    # Get coefficients and eigenvalues
    wfcoef = np.zeros((norb_total, norb_total), dtype=np.float64)
    fc.get_wfcoef(wfcoef=wfcoef, ikp=1)
    eigenvalues = fc.get_eigen(ikp=1, norb=norb_total)  # Orbital energies in eV
    
    # Per-atom orbital info
    norb_per = np.array([sparse.num_orb[np.where(sparse.nzx[:dims.nspecies] == Z)[0][0]] 
                         for Z in atomTypes])
    atom_coeff_offs = sparse.degelec[:natoms]
    
    # Setup OpenCL
    fdata_dir = os.path.join(os.path.dirname(__file__), "Fdata", "basis")
    projector = ocl_grid.GridProjector(fdata_dir, verbosity=0)
    projector.load_basis(species_nz=sorted(set(atomTypes)))
    
    # Grid setup from CLI args
    Lz = 6.0  # Larger z extent to ensure z_height plane is inside grid
    atom_center = atomPos.mean(axis=0)
    grid_origin = atom_center - np.array([Lx/2.0, Ly/2.0, Lz/2.0])  # Center grid on molecule
    
    grid_spec = {
        'origin': grid_origin,
        'dA': np.array([Lx/nx, 0.0, 0.0]),
        'dB': np.array([0.0, Ly/ny, 0.0]),
        'dC': np.array([0.0, 0.0, Lz/nz]),
        'ngrid': np.array([nx, ny, nz], dtype=np.int32)
    }
    
    atoms_dict = {
        'pos': atomPos.astype(np.float32),
        'type': atomTypes.astype(np.int32),
        'Rcut': np.array([6.0 if z == 8 else 4.0 for z in atomTypes], dtype=np.float32)
    }
    
    extent = [grid_origin[0], grid_origin[0] + Lx, grid_origin[1], grid_origin[1] + Ly]
    
    # Determine which MOs to plot
    nmo_total = wfcoef.shape[0]
    
    if args.orbitals == 'auto':
        # Auto = occupied orbitals only (electrons/2)
        # For H2O: 10 electrons = 5 occupied
        # For CH4: 8 electrons = 4 occupied  
        # Need to determine from system
        nelec = nmo_total * 2  # rough estimate
        nocc = nmo_total // 2 if nmo_total <= 10 else 5
        mo_indices = list(range(nocc))
        print(f"\nAuto-selected {nocc} occupied orbitals: {[i+1 for i in mo_indices]}")
    else:
        mo_indices = parse_orbital_indices(args.orbitals, nmo_total)
        print(f"\nUser-selected orbitals: {[i+1 for i in mo_indices]}")
    
    # Default labels
    mo_labels = []
    for idx in mo_indices:
        if idx == 0:
            mo_labels.append('O1s' if mol_name == 'H2O' else f'MO{idx+1}')
        elif idx == 1:
            mo_labels.append('O2s' if mol_name == 'H2O' else f'MO{idx+1}')
        elif nmo_total <= 10 and idx == nmo_total - 1:
            mo_labels.append('LUMO')
        elif nmo_total <= 10 and idx == nmo_total - 2:
            mo_labels.append('HOMO')
        elif nmo_total <= 10 and idx == nmo_total - 3:
            mo_labels.append('HOMO-1')
        elif nmo_total <= 10 and idx == nmo_total - 4:
            mo_labels.append('HOMO-2')
        else:
            mo_labels.append(f'MO{idx+1}')
    
    # Create text file for coefficient comparison
    coef_comp_file = os.path.join(export_dir, 'fortran_vs_python_coefficients.txt')
    with open(coef_comp_file, 'w') as fcomp:
        fcomp.write("=" * 70 + "\n")
        fcomp.write("FORT vs PYTH: COEFFICIENT EXTRACTION VERIFICATION\n")
        fcomp.write("=" * 70 + "\n\n")
        
        # Loop over ALL MOs to compare extraction
        for mo_idx in range(wfcoef.shape[0]):
            mo_fortran = mo_idx + 1
            
            fcomp.write(f"\n{'=' * 70}\n")
            fcomp.write(f"MO {mo_fortran} (Python index {mo_idx})\n")
            fcomp.write(f"{'=' * 70}\n")
            
            # Get Fortran output
            import io
            from contextlib import redirect_stdout
            fort_capture = io.StringIO()
            with redirect_stdout(fort_capture):
                fc.print_orb_coefs(iMO=mo_fortran, ikpoint=1)
            fcomp.write("FORT:\n" + fort_capture.getvalue())
            
            # Get Python extraction (row indexing)
            fcomp.write("PYTH (wfcoef[mo_idx, :] row):\n")
            for ia in range(natoms):
                i0 = atom_coeff_offs[ia]
                no = norb_per[ia]
                py_vals = wfcoef[mo_idx, i0:i0+no]
                fcomp.write(f"  Atom {ia} ({enames[ia]}): [{i0}:{i0+no}] = ")
                fcomp.write(", ".join([f"{v:+.6f}" for v in py_vals]) + "\n")
            
            # Also try column indexing for comparison
            if wfcoef.ndim == 2 and mo_idx < wfcoef.shape[1]:
                fcomp.write("PYTH_ALT (wfcoef[:, mo_idx] column):\n")
                for ia in range(natoms):
                    i0 = atom_coeff_offs[ia]
                    no = norb_per[ia]
                    py_vals_col = wfcoef[i0:i0+no, mo_idx]
                    fcomp.write(f"  Atom {ia} ({enames[ia]}): [{i0}:{i0+no}, {mo_idx}] = ")
                    fcomp.write(", ".join([f"{v:+.6f}" for v in py_vals_col]) + "\n")
        
        fcomp.write("\n" + "=" * 70 + "\n")
        fcomp.write("END OF COMPARISON\n")
    
    print(f"\nCoefficient comparison written to: {coef_comp_file}")
    
    for mo_idx, mo_label in zip(mo_indices, mo_labels):
        mo_fortran = mo_idx + 1  # Fortran uses 1-based
        
        # Extract coefficients for this MO
        numorb_max = max(norb_per)
        coeffs_fortran = np.zeros((natoms, numorb_max), dtype=np.float64)
        for ia in range(natoms):
            i0 = atom_coeff_offs[ia]
            no = norb_per[ia]
            coeffs_fortran[ia, :no] = wfcoef[mo_idx, i0:i0+no]
        
        # DEBUG MODE: Override coefficients for single-basis testing
        if args.debug_atom is not None and args.debug_orb is not None:
            debug_ia = args.debug_atom
            debug_orb = args.debug_orb
            debug_val = args.debug_val
            
            print(f"\n{'='*60}")
            print(f"DEBUG MODE: Testing single basis function")
            print(f"  Atom {debug_ia} ({enames[debug_ia]}), orbital={debug_orb}, value={debug_val}")
            print(f"{'='*60}")
            
            # Create artificial wfcoef: all zeros except one coefficient
            wfcoef_debug = np.zeros((norb_total, norb_total), dtype=np.float64)
            
            # Map orbital type to Fortran coefficient index
            # Fortran order: [s, py, pz, px] for O, [s] for H
            orb_idx_map = {'s': 0, 'py': 1, 'pz': 2, 'px': 3}
            if debug_orb not in orb_idx_map:
                raise ValueError(f"Unknown orbital type: {debug_orb}")
            
            orb_idx = orb_idx_map[debug_orb]
            coeff_global_idx = atom_coeff_offs[debug_ia] + orb_idx
            
            # Set the single coefficient
            wfcoef_debug[0, coeff_global_idx] = debug_val  # Put in MO index 0
            
            # Set in Fortran
            fc.set_wfcoef(wfcoef_debug[0, :], iMO=1, ikp=1)
            
            # Verify with print
            print("\nVerifying set coefficients (Fortran print_orb_coefs):")
            fc.print_orb_coefs(iMO=1, ikpoint=1)
            
            # Override coeffs_fortran for this debug run
            coeffs_fortran = np.zeros((natoms, numorb_max), dtype=np.float64)
            coeffs_fortran[debug_ia, orb_idx] = debug_val
            
            # Update mo_fortran to use MO 1 (our artificial MO)
            mo_fortran = 1
            mo_label = f"DEBUG_{enames[debug_ia]}_{debug_orb}"
        
        # Print Fortran coefficients for reference
        print(f"\n{'='*60}")
        print(f"MO {mo_fortran} ({mo_label})  E = {eigenvalues[mo_idx]:.4f} eV")
        fc.print_orb_coefs(iMO=mo_fortran, ikpoint=1)
        
        # Print Python extraction
        print(f"Python extraction (row {mo_idx}):")
        for ia in range(natoms):
            i0 = atom_coeff_offs[ia]
            no = norb_per[ia]
            py_vals = wfcoef[mo_idx, i0:i0+no]
            print(f"  Atom {ia} ({enames[ia]}): {np.round(py_vals, 6)}")
        
        # Remap to OpenCL order [px, py, pz, s] matching DFT/utils.py convCoefs()
        # Fortran [s, py, pz, px] -> OpenCL [px, py, pz, s]
        coeffs_opencl = np.zeros((natoms, 4), dtype=np.float64)
        for ia in range(natoms):
            no = norb_per[ia]
            if no == 1:  # H atom: just s
                coeffs_opencl[ia, 3] = coeffs_fortran[ia, 0]  # s -> index 3
            elif no == 4:  # O atom: s, py, pz, px -> px, py, pz, s
                coeffs_opencl[ia, 0] = -coeffs_fortran[ia, 3]  # px -> index 0 (FLIPPED)
                coeffs_opencl[ia, 1] = coeffs_fortran[ia, 1]  # py -> index 1
                coeffs_opencl[ia, 2] = coeffs_fortran[ia, 2]  # pz -> index 2
                coeffs_opencl[ia, 3] = coeffs_fortran[ia, 0]  # s -> index 3
        
        # Build coefficient text
        coeffs_text_f = []
        coeffs_text_o = []
        for ia in range(natoms):
            co = norb_per[ia]
            cf_str = f"{enames[ia]}: " + ",".join([f"{coeffs_fortran[ia,i]:+.3f}" for i in range(co)])
            # Show all 4 OpenCL values [px, py, pz, s]
            co_str = f"{enames[ia]}: {coeffs_opencl[ia,0]:+.3f},{coeffs_opencl[ia,1]:+.3f},{coeffs_opencl[ia,2]:+.3f},{coeffs_opencl[ia,3]:+.3f}"
            coeffs_text_f.append(cf_str)
            coeffs_text_o.append(co_str)
        coeffs_text_f = "Fortran [s,py,pz,px]:\n" + "\n".join(coeffs_text_f)
        coeffs_text_o = "OpenCL [px,py,pz,s]:\n" + "\n".join(coeffs_text_o)
        
        # Fortran projection at z = z_height (1.0 Å above molecular plane)
        X = np.linspace(grid_origin[0] + 8.0/nx * 0.5, grid_origin[0] + 8.0/nx * (nx-0.5), nx)
        Y = np.linspace(grid_origin[1] + 8.0/ny * 0.5, grid_origin[1] + 8.0/ny * (ny-0.5), ny)
        XX, YY = np.meshgrid(X, Y, indexing='ij')
        Z_plane = np.zeros_like(XX) + atom_center[2] + z_height
        points = np.stack([XX.ravel(), YY.ravel(), Z_plane.ravel()], axis=1)
        
        # Fortran uses 1-based MO indexing
        psi_fort_flat = fc.orb2points(points, iMO=mo_fortran, ikpoint=1)
        psi_fort = psi_fort_flat.reshape(nx, ny)
        
        # OpenCL projection - get 3D then extract z=z_height slice
        psi_ocl_3d = projector.project_orbital(
            coeffs=coeffs_opencl.astype(np.float32),
            norb_per=norb_per,
            atoms_dict=atoms_dict,
            grid_spec=grid_spec,
            nMaxAtom=64
        )
        # Extract z=z_height slice
        iz_target = int((atom_center[2] + z_height - grid_origin[2]) / grid_spec['dC'][2])
        iz_target = max(0, min(iz_target, nz - 1))
        psi_ocl = psi_ocl_3d[:, :, iz_target]
        
        # Plot side-by-side
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        
        E_str = f"E={eigenvalues[mo_idx]:.2f}eV"
        plot_orbital_comparison(axes[0], psi_fort, atomPos, atomTypes, extent, f"Fortran MO {mo_fortran} ({mo_label}) {E_str}", coeffs_text_f, enames)
        plot_orbital_comparison(axes[1], psi_ocl, atomPos, atomTypes, extent, f"OpenCL MO {mo_fortran} ({mo_label}) {E_str}", coeffs_text_o, enames)
        
        # Compute correlation
        corr = np.corrcoef(psi_fort.flatten(), psi_ocl.flatten())[0, 1]
        
        # Add stats
        fig.text(0.5, 0.02, f"z={z_height}Å | Fortran: [{psi_fort.min():.3f}, {psi_fort.max():.3f}] | "
                           f"OpenCL: [{psi_ocl.min():.3f}, {psi_ocl.max():.3f}] | Corr: {corr:.3f}",
                 ha='center', fontsize=10)
        
        plt.tight_layout(rect=[0, 0.03, 1, 0.97])
        
        filename = f"h2o_MO{mo_fortran}_{mo_label}_z{z_height}_fortran_vs_opencl.png"
        filepath = os.path.join(export_dir, filename)
        plt.savefig(filepath, dpi=150, bbox_inches='tight')
        print(f"\nSaved: {filepath}")
        print(f"  Fortran: [{psi_fort.min():.4f}, {psi_fort.max():.4f}]")
        print(f"  OpenCL:  [{psi_ocl.min():.4f}, {psi_ocl.max():.4f}]")
        print(f"  Correlation: {corr:.4f}")
    
    plt.close('all')
    print(f"\nAll plots saved to: {export_dir}")

if __name__ == "__main__":
    main()
