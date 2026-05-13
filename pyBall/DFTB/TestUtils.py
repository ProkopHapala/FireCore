"""
Test utilities for waveplot comparison and debugging.

Provides high-level functions for comparing libwaveplot and OpenCL evaluations,
including plotting utilities and RMS error calculations.
"""

import numpy as np
import os
import json
import hashlib
from pathlib import Path
from typing import Dict, Any, Optional, Union, Tuple


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



def load_xyz(fname):
    """Simple XYZ parser for atomic types and positions."""
    ELEM_Z = {'H':1, 'C':6, 'N':7, 'O':8, 'P':15, 'S':16}
    with open(fname, 'r') as f:
        lines = f.readlines()
    natoms = int(lines[0])
    atomTypes = []; atomPos = []
    for line in lines[2:2+natoms]:
        p = line.split()
        sym = p[0]
        z = ELEM_Z.get(sym, 0)
        atomTypes.append(z)
        atomPos.append([float(p[1]), float(p[2]), float(p[3])])
    return np.array(atomTypes, dtype=np.int32), np.array(atomPos, dtype=np.float64)

def check_density_integration(rho, step, expected_n=None, label=""):
    """Check and print the integrated electron density."""
    dV = step**3
    total_q = np.sum(rho) * dV
    msg = f"[{label}] Integrated electrons: {total_q:.4f}"
    if expected_n is not None:
        msg += f" (expected: {expected_n:.2f})"
    print(msg)
    return total_q

def get_z_profile(rho, pos, origin, step):
    """Extract 1D z-profile from 3D grid above a specific XY position."""
    nx, ny, nz = rho.shape
    ix = int(round((pos[0] - origin[0]) / step))
    iy = int(round((pos[1] - origin[1]) / step))
    ix = np.clip(ix, 0, nx-1)
    iy = np.clip(iy, 0, ny-1)
    
    profile = rho[ix, iy, :]
    z = origin[2] + np.arange(nz) * step
    return z, profile

def atom_to_grid_idx(atom_pos, origin, step, ngrid):
    """Convert atom position to nearest grid indices.

    Args:
        atom_pos: (3,) array of position in Angstrom
        origin: (3,) array of grid origin
        step: grid spacing in Angstrom
        ngrid: (3,) array of grid dimensions or int for 1D

    Returns:
        (ix, iy, iz) clipped to [0, ngrid-1]
    """
    frac = (np.asarray(atom_pos) - np.asarray(origin)) / step
    idx = np.round(frac).astype(np.int32)
    ngrid_arr = np.array(ngrid, dtype=np.int32)
    if ngrid_arr.size == 1:
        ngrid_arr = np.array([ngrid_arr, ngrid_arr, ngrid_arr], dtype=np.int32)
    idx = np.clip(idx, [0, 0, 0], ngrid_arr - 1)
    return idx[0], idx[1], idx[2]

def extract_z_profile(field, atom_pos, origin, step, z_distances=None):
    """Extract z-profile from 3D field at atom XY position.

    Args:
        field: (nx, ny, nz) 3D array
        atom_pos: (3,) target position in Angstrom
        origin: (3,) grid origin
        step: grid spacing
        z_distances: optional 1D array of z offsets from atom_pos[2];
                     if None returns full grid column (z_grid, values)

    Returns:
        If z_distances is None: (z_grid, column)
        If z_distances given: interpolated values at z = atom_pos[2] + z_distances
    """
    tix, tiy, _ = atom_to_grid_idx(atom_pos, origin, step, field.shape)
    z_grid_vals = origin[2] + np.arange(field.shape[2]) * step
    col = field[tix, tiy, :]
    if z_distances is None:
        return z_grid_vals, col
    z_abs = atom_pos[2] + z_distances
    return np.interp(z_abs, z_grid_vals, col, left=col[0], right=col[-1])

def get_dftb_zscan_energies(log_path):
    """Parse DFTB total energies from a z-scan log file."""
    zvals = []
    energies = []
    if not os.path.exists(log_path):
        return None, None
    with open(log_path, 'r') as f:
        lines = f.readlines()
    
    # Simple parser for the table format in zscan_results.txt
    start_line = -1
    for i, line in enumerate(lines):
        if ('z [' in line and 'E [' in line) or ('-------' in line and i > 5):
            if 'z [' in line:
                start_line = i + 2
            else:
                start_line = i + 1
            break
    
    if start_line != -1:
        for line in lines[start_line:]:
            if not line.strip() or '---' in line: break
            p = line.split()
            if len(p) >= 3:
                zvals.append(float(p[0]))
                energies.append(float(p[2]))
    
    return np.array(zvals), np.array(energies)


def write_cube(filename, data, atoms, origin, vectors, title="Electron density"):
    """
    Write Gaussian cube file.

    Args:
        filename: output file path
        data: (nx, ny, nz) array (any memory layout; written z-fastest)
        atoms: list of tuples (atomic_number, x, y, z) in same units as origin/vectors
        origin: (3,) array
        vectors: (3, 3) array; each row is a grid step vector
        title: first comment line
    """
    nx, ny, nz = data.shape
    with open(filename, 'w') as f:
        f.write(f"{title}\n")
        f.write("Generated by FireCore test suite\n")
        f.write(f"{len(atoms):5d} {origin[0]:12.6f} {origin[1]:12.6f} {origin[2]:12.6f}\n")
        f.write(f"{nx:5d} {vectors[0,0]:12.6f} {vectors[0,1]:12.6f} {vectors[0,2]:12.6f}\n")
        f.write(f"{ny:5d} {vectors[1,0]:12.6f} {vectors[1,1]:12.6f} {vectors[1,2]:12.6f}\n")
        f.write(f"{nz:5d} {vectors[2,0]:12.6f} {vectors[2,1]:12.6f} {vectors[2,2]:12.6f}\n")
        for atom in atoms:
            z, x, y, zc = atom
            f.write(f"{z:5d} {float(z):12.6f} {x:12.6f} {y:12.6f} {zc:12.6f}\n")
        # Write data z-fastest (cube convention)
        count = 0
        for ix in range(nx):
            for iy in range(ny):
                for iz in range(nz):
                    f.write(f"{data[ix, iy, iz]:14.5e}")
                    count += 1
                    if count % 6 == 0:
                        f.write("\n")
                # End of z-slice for fixed iy
        if count % 6 != 0:
            f.write("\n")


# ═══════════════════════════════════════════════════════════════════════════════
# Parity Checkpoint System (Generic, Reusable)
# ═══════════════════════════════════════════════════════════════════════════════

def _hash_config(config: Dict[str, Any]) -> str:
    """
    Compute a stable hash of configuration dictionary.

    IMPORTANT: This hash is ONLY for configuration parameters (grid size,
    physics parameters, basis set, etc.) to detect when the computational
    setup changes. Data arrays are NOT hashed - they are compared using
    numerical tolerances (RMS/max error) to handle float32/float64 precision
    differences and operation ordering variations.

    Args:
        config: Dictionary of configuration parameters

    Returns:
        Hexadecimal hash string
    """
    # Convert to JSON with sorted keys for stable hashing
    # Convert numpy arrays to lists for JSON serialization
    def sanitize(obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, dict):
            return {k: sanitize(v) for k, v in sorted(obj.items())}
        elif isinstance(obj, (list, tuple)):
            return [sanitize(v) for v in obj]
        else:
            return obj

    config_str = json.dumps(sanitize(config), sort_keys=True)
    return hashlib.sha256(config_str.encode()).hexdigest()


def save_checkpoint(
    checkpoint_dir: Union[str, Path],
    checkpoint_name: str,
    data: Dict[str, Any],
    config: Optional[Dict[str, Any]] = None,
    overwrite: bool = False
) -> Path:
    """
    Save checkpoint data and configuration to disk.

    This is a generic function for saving computational checkpoints
    that can be loaded later for parity testing.

    Args:
        checkpoint_dir: Directory to store checkpoints
        checkpoint_name: Name of the checkpoint (e.g., 'density', 'pauli_field')
        data: Dictionary of data to save (arrays, scalars, etc.)
        config: Dictionary of configuration parameters that affect the computation
                (e.g., grid parameters, physics parameters, basis set info)
        overwrite: If False, raise error if checkpoint already exists

    Returns:
        Path to the saved checkpoint file

    Example:
        >>> data = {'rho': rho_grid, 'origin': origin, 'step': step}
        >>> config = {'step': 0.2, 'margin': 2.0, 'basis': 'mio-1-1'}
        >>> save_checkpoint('/tmp/checkpoints', 'density', data, config)
    """
    checkpoint_dir = Path(checkpoint_dir)
    checkpoint_dir.mkdir(parents=True, exist_ok=True)

    checkpoint_file = checkpoint_dir / f"{checkpoint_name}.npz"
    meta_file = checkpoint_dir / f"{checkpoint_name}.json"

    if not overwrite and checkpoint_file.exists():
        raise FileExistsError(f"Checkpoint already exists: {checkpoint_file}")

    # Compute config hash if provided
    config_hash = None
    if config is not None:
        config_hash = _hash_config(config)

    # Save data arrays
    np.savez_compressed(checkpoint_file, **data)

    # Save metadata
    metadata = {
        'checkpoint_name': checkpoint_name,
        'config_hash': config_hash,
        'config': config,
        'data_keys': list(data.keys()),
        'data_shapes': {k: (v.shape if hasattr(v, 'shape') else str(type(v)))
                        for k, v in data.items()}
    }
    with open(meta_file, 'w') as f:
        json.dump(metadata, f, indent=2, default=str)

    print(f"[Checkpoint] Saved: {checkpoint_file}")
    if config_hash:
        print(f"[Checkpoint] Config hash: {config_hash[:16]}...")

    return checkpoint_file


def load_checkpoint(
    checkpoint_dir: Union[str, Path],
    checkpoint_name: str,
    config: Optional[Dict[str, Any]] = None,
    strict_config: bool = True
) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """
    Load checkpoint data and metadata from disk.

    Args:
        checkpoint_dir: Directory containing checkpoints
        checkpoint_name: Name of the checkpoint to load
        config: Configuration dictionary to verify against saved config
        strict_config: If True, verify config hash matches (raises error if mismatch)

    Returns:
        (data, metadata) tuple where:
            data: Dictionary of loaded data
            metadata: Dictionary of metadata (config, shapes, etc.)

    Raises:
        FileNotFoundError: If checkpoint file doesn't exist
        ValueError: If config hash mismatch and strict_config=True

    Example:
        >>> data, meta = load_checkpoint('/tmp/checkpoints', 'density')
        >>> rho = data['rho']
        >>> origin = data['origin']
    """
    checkpoint_dir = Path(checkpoint_dir)
    checkpoint_file = checkpoint_dir / f"{checkpoint_name}.npz"
    meta_file = checkpoint_dir / f"{checkpoint_name}.json"

    if not checkpoint_file.exists():
        raise FileNotFoundError(f"Checkpoint not found: {checkpoint_file}")

    # Load data
    data = dict(np.load(checkpoint_file))

    # Load metadata
    with open(meta_file, 'r') as f:
        metadata = json.load(f)

    # Verify config hash if provided
    if config is not None and strict_config:
        saved_hash = metadata.get('config_hash')
        if saved_hash:
            current_hash = _hash_config(config)
            if saved_hash != current_hash:
                raise ValueError(
                    f"Config hash mismatch for checkpoint '{checkpoint_name}':\n"
                    f"  Saved: {saved_hash[:16]}...\n"
                    f"  Current: {current_hash[:16]}...\n"
                    f"  Set strict_config=False to bypass verification."
                )

    print(f"[Checkpoint] Loaded: {checkpoint_file}")
    return data, metadata


def compare_checkpoint(
    checkpoint_dir: Union[str, Path],
    checkpoint_name: str,
    new_data: Dict[str, Any],
    config: Optional[Dict[str, Any]] = None,
    tolerances: Optional[Dict[str, float]] = None,
    strict_config: bool = True,
    verbose: bool = True
) -> Dict[str, Any]:
    """
    Load reference checkpoint and compare with new data.

    Args:
        checkpoint_dir: Directory containing checkpoints
        checkpoint_name: Name of the checkpoint to compare
        new_data: Dictionary of new data to compare against reference
        config: Configuration dictionary for hash verification
        tolerances: Dictionary of tolerances per key (e.g., {'rho': 1e-6})
                    If not provided, uses default tolerances
        strict_config: If True, verify config hash matches
        verbose: If True, print comparison results

    Returns:
        Dictionary with comparison results:
            - 'match': bool, overall match status
            - 'details': dict with per-key comparison results
            - 'rms': dict with RMS errors per array key
            - 'max_err': dict with max absolute errors per array key

    Raises:
        FileNotFoundError: If reference checkpoint doesn't exist
        ValueError: If config hash mismatch or shape mismatch

    Example:
        >>> result = compare_checkpoint('/tmp/checkpoints', 'density',
        ...                             {'rho': new_rho, 'origin': origin},
        ...                             tolerances={'rho': 1e-6})
        >>> if result['match']:
        ...     print("Parity achieved!")
    """
    # Load reference
    ref_data, metadata = load_checkpoint(
        checkpoint_dir, checkpoint_name, config, strict_config
    )

    # Default tolerances
    if tolerances is None:
        tolerances = {
            'default': 1e-6,  # For float arrays
            'position': 1e-8,  # For positions (Å)
            'energy': 1e-4,    # For energies (eV)
            'force': 1e-4,    # For forces (eV/Å)
            'density': 1e-5,  # For densities (e/Å³)
            'gradient': 1e-5  # For gradients
        }

    # Compare each key
    details = {}
    rms_errors = {}
    max_errors = {}
    overall_match = True

    for key in ref_data.keys():
        if key not in new_data:
            if verbose:
                print(f"[Compare] WARNING: Key '{key}' in reference but not in new data")
            overall_match = False
            continue

        ref_val = ref_data[key]
        new_val = new_data[key]

        # Get tolerance for this key
        tol = tolerances.get(key, tolerances.get('default', 1e-6))

        # Handle arrays
        if isinstance(ref_val, np.ndarray):
            if not isinstance(new_val, np.ndarray):
                error_msg = f"Type mismatch for '{key}': ref=ndarray, new={type(new_val)}"
                details[key] = {'match': False, 'error': error_msg}
                overall_match = False
                if verbose:
                    print(f"[Compare] {error_msg}")
                continue

            if ref_val.shape != new_val.shape:
                error_msg = f"Shape mismatch for '{key}': ref={ref_val.shape}, new={new_val.shape}"
                details[key] = {'match': False, 'error': error_msg}
                overall_match = False
                if verbose:
                    print(f"[Compare] {error_msg}")
                continue

            # Compute RMS and max error
            rms, max_err = compute_rms_error(ref_val, new_val)
            rms_errors[key] = rms
            max_errors[key] = max_err

            # Check tolerance
            match = max_err < tol
            details[key] = {
                'match': match,
                'rms': rms,
                'max_err': max_err,
                'tolerance': tol
            }

            if not match:
                overall_match = False

            if verbose:
                status = "PASS" if match else "FAIL"
                print(f"[Compare] {key:20s}: {status}  RMS={rms:.3e}  max={max_err:.3e}  tol={tol:.1e}")

        # Handle scalars (including floats, ints, etc.)
        else:
            try:
                diff = abs(float(ref_val) - float(new_val))
            except (TypeError, ValueError):
                error_msg = f"Cannot compare '{key}': ref={type(ref_val)}, new={type(new_val)}"
                details[key] = {'match': False, 'error': error_msg}
                overall_match = False
                if verbose:
                    print(f"[Compare] {error_msg}")
                continue

            match = diff < tol
            details[key] = {
                'match': match,
                'diff': diff,
                'tolerance': tol,
                'ref': ref_val,
                'new': new_val
            }

            if not match:
                overall_match = False

            if verbose:
                status = "PASS" if match else "FAIL"
                print(f"[Compare] {key:20s}: {status}  diff={diff:.3e}  tol={tol:.1e}")

    result = {
        'match': overall_match,
        'details': details,
        'rms': rms_errors,
        'max_err': max_errors
    }

    if verbose:
        print(f"[Compare] Overall: {'PASS' if overall_match else 'FAIL'}")

    return result


class CheckpointManager:
    """
    High-level manager for parity checkpoints.

    Provides a convenient interface for saving/loading checkpoints
    with automatic directory management and configuration tracking.

    Example:
        >>> mgr = CheckpointManager('/tmp/checkpoints', 'afm_pipeline')
        >>> mgr.save('density', {'rho': rho}, config={'step': 0.2})
        >>> data = mgr.load('density')
        >>> result = mgr.compare('density', {'rho': new_rho})
    """

    def __init__(self, base_dir: Union[str, Path], suite_name: str):
        """
        Initialize checkpoint manager.

        Args:
            base_dir: Base directory for all checkpoints
            suite_name: Name of the test suite (creates subdirectory)
        """
        self.base_dir = Path(base_dir)
        self.suite_name = suite_name
        self.checkpoint_dir = self.base_dir / suite_name
        self.checkpoint_dir.mkdir(parents=True, exist_ok=True)

    def get_path(self, checkpoint_name: str) -> Path:
        """Get full path to checkpoint file."""
        return self.checkpoint_dir / f"{checkpoint_name}.npz"

    def save(self, checkpoint_name: str, data: Dict[str, Any],
             config: Optional[Dict[str, Any]] = None, overwrite: bool = False) -> Path:
        """Save checkpoint."""
        return save_checkpoint(
            self.checkpoint_dir, checkpoint_name, data, config, overwrite
        )

    def load(self, checkpoint_name: str,
             config: Optional[Dict[str, Any]] = None, strict_config: bool = True) -> Dict[str, Any]:
        """Load checkpoint (returns data only)."""
        data, _ = load_checkpoint(
            self.checkpoint_dir, checkpoint_name, config, strict_config
        )
        return data

    def compare(self, checkpoint_name: str, new_data: Dict[str, Any],
                config: Optional[Dict[str, Any]] = None,
                tolerances: Optional[Dict[str, float]] = None,
                strict_config: bool = True, verbose: bool = True) -> Dict[str, Any]:
        """Compare new data against checkpoint."""
        return compare_checkpoint(
            self.checkpoint_dir, checkpoint_name, new_data, config,
            tolerances, strict_config, verbose
        )

    def exists(self, checkpoint_name: str) -> bool:
        """Check if checkpoint exists."""
        return self.get_path(checkpoint_name).exists()
