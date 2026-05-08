"""
Test utilities for waveplot comparison and debugging.

Provides high-level functions for comparing libwaveplot and OpenCL evaluations,
including plotting utilities and RMS error calculations.
"""

import numpy as np
from pathlib import Path


def compute_rms_error(vals1, vals2):
    """
    Compute RMS error between two arrays.
    
    Args:
        vals1: array of values
        vals2: array of values (same shape)
    
    Returns:
        rms: float, root-mean-square error
        max_err: float, maximum absolute error
    """
    diff = vals1 - vals2
    rms = np.sqrt(np.mean(diff**2))
    max_err = np.max(np.abs(diff))
    return rms, max_err


def compare_point_evaluations(wp_vals, ocl_vals, mo_indices, energies, homo):
    """
    Compare libwaveplot and OpenCL point evaluations.
    
    Args:
        wp_vals: list of arrays, each [npoints] from libwaveplot
        ocl_vals: list of arrays, each [npoints] from OpenCL
        mo_indices: list of MO indices
        energies: list of energies (eV)
        homo: HOMO index
    
    Returns:
        results: list of dicts with rms, max_err for each MO
    """
    results = []
    for i, (wp, ocl) in enumerate(zip(wp_vals, ocl_vals)):
        rms, max_err = compute_rms_error(wp, ocl)
        mo_idx = mo_indices[i]
        tag = " [HOMO]" if mo_idx == homo else (" [LUMO]" if mo_idx == homo+1 else "")
        results.append({
            'mo_index': mo_idx,
            'tag': tag,
            'energy': energies[i],
            'rms': rms,
            'max_err': max_err
        })
    return results


def compare_grid_evaluations(grid_ref, grid_lib, mo_indices, energies, homo):
    """
    Compare reference cube and libwaveplot grid evaluations.
    
    Args:
        grid_ref: list of arrays, each (nx, ny, nz) from reference cube
        grid_lib: list of arrays, each (nx, ny, nz) from libwaveplot
        mo_indices: list of MO indices
        energies: list of energies (eV)
        homo: HOMO index
    
    Returns:
        results: list of dicts with rms, max_err for each MO
    """
    results = []
    for i, (ref, lib) in enumerate(zip(grid_ref, grid_lib)):
        rms, max_err = compute_rms_error(ref, lib)
        mo_idx = mo_indices[i]
        tag = " [HOMO]" if mo_idx == homo else (" [LUMO]" if mo_idx == homo+1 else "")
        results.append({
            'mo_index': mo_idx,
            'tag': tag,
            'energy': energies[i],
            'rms': rms,
            'max_err': max_err
        })
    return results


def print_comparison_results(results, method_name=""):
    """
    Print comparison results in a formatted table.
    
    Args:
        results: list of dicts from compare_*_evaluations
        method_name: string to identify the method (e.g., "orb2points")
    """
    print(f"\n  RMS ({method_name}):")
    for r in results:
        print(f"    MO {r['mo_index']:3d}{r['tag']:8s} E={r['energy']:8.3f}eV  RMS={r['rms']:.3e}  max={r['max_err']:.3e}")


def generate_2d_point_grid(plane='xy', npoints=64, z_offset=0.0, xy_range=None):
    """
    Generate 2D grid of points for orbital evaluation.
    
    Args:
        plane: 'xy', 'xz', or 'yz'
        npoints: number of points per axis
        z_offset: fixed coordinate for out-of-plane axis (Å)
        xy_range: optional (min, max) for in-plane axes (Å)
    
    Returns:
        points: (npoints*npoints, 3) array in Angstrom
        extent: [xmin, xmax, ymin, ymax] for plotting
    """
    if xy_range is None:
        xy_range = (-3.0, 8.0)
    x = np.linspace(xy_range[0], xy_range[1], npoints)
    y = np.linspace(xy_range[0], xy_range[1], npoints)
    X, Y = np.meshgrid(y, x)    
    if plane == 'xy':
        points = np.column_stack([X.ravel(), Y.ravel(), np.full(X.size, z_offset)])
        extent = [xy_range[0], xy_range[1], xy_range[0], xy_range[1]]
    elif plane == 'xz':
        points = np.column_stack([X.ravel(), np.full(X.size, z_offset), Y.ravel()])
        extent = [xy_range[0], xy_range[1], xy_range[0], xy_range[1]]
    elif plane == 'yz':
        points = np.column_stack([np.full(X.size, z_offset), X.ravel(), Y.ravel()])
        extent = [xy_range[0], xy_range[1], xy_range[0], xy_range[1]]
    else:
        raise ValueError(f"Unknown plane: {plane}")
    return points, extent


def generate_1d_z_scan(npoints=100, z_range=(-3.0, 3.0)):
    """
    Generate 1D scan along z-axis.
    
    Args:
        npoints: number of points
        z_range: (zmin, zmax) in Angstrom
    
    Returns:
        points: (npoints, 3) array with x=y=0
        z_vals: (npoints,) z coordinates
    """
    z = np.linspace(z_range[0], z_range[1], npoints)
    points = np.column_stack([np.zeros_like(z), np.zeros_like(z), z])
    return points, z


def print_eigenvecs(eigenvec_path, detailed_xml_path=None, waveplot_in_path=None, max_orbitals=None):
    """
    Pretty print eigenvectors from eigenvec.bin with atom/orbital labels.
    
    Args:
        eigenvec_path: path to eigenvec.bin file
        detailed_xml_path: optional path to detailed.xml for atom indices
        waveplot_in_path: optional path to waveplot_in.hsd for orbital info
        max_orbitals: maximum number of orbitals to print (None for all)
    """
    import sys
    sys.path.insert(0, str(Path(__file__).parent.parent.parent))
    from pyBall.DFTB.DFTBplusParser import parse_eigenvec_bin_custom, parse_detailed_xml_custom, parse_basis_hsd_ang
    
    # Parse detailed.xml for atom indices
    species_per_atom = None
    species_names = None
    if detailed_xml_path:
        data = parse_detailed_xml_custom(detailed_xml_path)
        nstates = data['nstates']
        norb = data['norb']
        species_per_atom = data['species_per_atom']
        species_names = data['species_names']
    else:
        # Read eigenvec.bin to get dimensions
        with open(eigenvec_path, 'rb') as f:
            raw = f.read()
        nstates = int.from_bytes(raw[0:4], byteorder='little', signed=False)
        nvals = (len(raw) - 4) // 8
        norb = nvals // nstates
    
    # Parse waveplot_in.hsd for orbital info
    species_orbitals = {}  # species -> list of orbital labels
    if waveplot_in_path:
        species_list = parse_basis_hsd_ang(waveplot_in_path)
        for sp in species_list:
            orb_labels = []
            for orb in sp['orbitals']:
                l = orb['l']
                if l == 0:
                    orb_labels.append('s')
                elif l == 1:
                    orb_labels.extend(['px', 'py', 'pz'])
                elif l == 2:
                    orb_labels.extend(['dxy', 'dyz', 'dz2', 'dxz', 'dx2'])
                else:
                    orb_labels.append(f'l{l}')
            species_orbitals[sp['name']] = orb_labels
    
    # Build atom/orbital map
    atom_orbital_map = []
    if species_per_atom is not None and species_names is not None and species_orbitals:
        orb_idx = 0
        for iatom, species_idx in enumerate(species_per_atom):
            species_name = species_names[species_idx]
            orb_names = species_orbitals.get(species_name, ['s'])
            
            for orb_name in orb_names:
                if orb_idx >= norb:
                    break
                label = f"{species_name}{iatom}{orb_name}"
                atom_orbital_map.append(label)
                orb_idx += 1
    elif species_per_atom is not None and species_names is not None:
        # Fallback: H has s, others have s,px,py,pz
        orb_idx = 0
        for iatom, species_idx in enumerate(species_per_atom):
            species_name = species_names[species_idx]
            if species_name == 'H':
                orb_names = ['s']
            else:
                orb_names = ['s', 'px', 'py', 'pz']
            
            for orb_name in orb_names:
                if orb_idx >= norb:
                    break
                label = f"{species_name}{iatom}{orb_name}"
                atom_orbital_map.append(label)
                orb_idx += 1
    else:
        atom_orbital_map = [f"AO{i:03d}" for i in range(norb)]
    
    # Pad to match norb if needed
    while len(atom_orbital_map) < norb:
        atom_orbital_map.append(f"AO{len(atom_orbital_map)}")
    
    # Parse eigenvectors
    evecs = parse_eigenvec_bin_custom(eigenvec_path, nstates, norb)
    
    if max_orbitals is not None:
        nstates = min(nstates, max_orbitals)
    
    print(f"\n{'='*80}")
    print(f"EIGENVECTORS from {eigenvec_path}")
    print(f"nStates={nstates}, nOrb={norb}")
    print(f"{'='*80}\n")
    
    # Build column headers
    col_headers = atom_orbital_map[:norb]
    
    # Calculate column width
    col_width = max(12, max(len(h) for h in col_headers))
    
    # Print header row
    header_str = f"{'MOs':<6}  {'Coefficients':<12}  |  " + "  ".join(f"{h:>{col_width}s}" for h in col_headers)
    print(header_str)
    print("-" * len(header_str))
    
    # Print each MO row
    for istate in range(nstates):
        mo_label = f"MO{istate:03d}"
        coeffs = [f"{coeff:{col_width}.6f}" for coeff in evecs[istate, :]]
        row_str = f"{mo_label:<6}  {'':<12}  |  " + "  ".join(coeffs)
        print(row_str)
