"""
Surface_utils.py - GridFF alignment verification and visualization utilities

High-level utilities for:
- Loading GridFF .npy files with metadata tracking
- Visualizing grids with atom overlay
- Sampling grids at atom positions
- Verifying grid-atom alignment
- Detecting proper shift conventions

Design principle: Glue layer that imports/reuses existing modules
(GridFF.py, RigidBodyAFM.py, InteractionEnergy.py) with minimal new code.
"""

import os
import sys
import numpy as np
import json

# matplotlib with non-interactive backend for headless use
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Import from existing modules
from .InteractionEnergy import load_xyz_with_REQs
from .RigidBodyAFM import sample_gridff_single_atom


# =============================================================================
# Section 1: GridFF I/O Utilities (copied/adapted from existing modules)
# =============================================================================

def load_gridff_metadata(grid_path):
    """
    Load GridFF metadata from JSON file if available.
    
    Args:
        grid_path: Path to Bspline_PLQd.npy or similar
        
    Returns:
        dict with keys: g0, dg, ns, lvec, z0, grid_type, generation_script
        or None if metadata file not found
    """
    # Derive metadata path from grid path
    base_dir = os.path.dirname(grid_path)
    base_name = os.path.splitext(os.path.basename(grid_path))[0]
    meta_path = os.path.join(base_dir, f"{base_name}_meta.json")
    
    if os.path.exists(meta_path):
        with open(meta_path, 'r') as f:
            metadata = json.load(f)
        print(f"Loaded GridFF metadata from: {meta_path}")
        return metadata
    else:
        print(f"Warning: GridFF metadata not found at {meta_path}")
        return None




def load_gridff_array(path):
    """
    Load GridFF from .npy file with validation.
    Copied from GridFFRelaxedScan.py
    
    Args:
        path: Path to Bspline_PLQd.npy or similar
        
    Returns:
        grid: (nx, ny, nz, nch) float32 array, nch=3 or 4
    """
    arr = np.load(path)
    if arr.ndim != 4:
        raise ValueError(f"GridFF must be 4D, got {arr.shape} from {path}")
    if arr.shape[3] == 3:
        arr4 = np.zeros(arr.shape[:3] + (4,), dtype=np.float32)
        arr4[:, :, :, :3] = arr.astype(np.float32)
        return arr4
    if arr.shape[3] == 4:
        return np.ascontiguousarray(arr.astype(np.float32))
    raise ValueError(f"GridFF channels must be 3 or 4, got {arr.shape} from {path}")


def find_generated_gridff(workdir, src_xyz):
    """
    Find generated GridFF file in workdir.
    Copied from GridFFRelaxedScan.py
    """
    bname = os.path.splitext(os.path.basename(src_xyz))[0]
    cands = [
        os.path.join(workdir, 'data', bname, 'Bspline_PLQd.npy'),
        os.path.join(workdir, 'data', bname, 'Bspline_PLQd_ocl.npy'),
    ]
    for p in cands:
        if os.path.exists(p):
            return p
    raise FileNotFoundError(f"Generated GridFF not found, tried {cands}")


def load_substrate_xyz_with_lvec(path):
    """
    Load substrate XYZ and extract lattice vectors from comment line.
    
    Args:
        path: Path to .xyz file
        
    Returns:
        dict with:
            'apos': (n,3) atom positions
            'enames': element names
            'lvec': (3,3) lattice vectors from comment line
            'z_top': float, max z coordinate
    """
    apos, REQs, enames, Zs, lvec = load_xyz_with_REQs(path)
    
    # Extract lattice vectors from first line if available
    with open(path, 'r') as f:
        n = int(f.readline())
        comment = f.readline().strip()
        # Try to parse lvec from comment
        lvec_from_comment = None
        for prefix in ["lvec:", "lvs"]:
            if prefix in comment:
                idx = comment.find(prefix) + len(prefix)
                parts = comment[idx:].split()
                try:
                    vals = [float(v) for v in parts if v.strip()]
                    if len(vals) >= 9:
                        lvec_from_comment = np.array(vals[:9]).reshape(3, 3).astype(np.float32)
                        break
                except ValueError:
                    pass
    
    # Use lvec from load_xyz_with_REQs if not found in comment
    if lvec_from_comment is not None:
        lvec = lvec_from_comment
    
    z_top = float(np.max(apos[:, 2]))
    
    return {
        'apos': apos,
        'enames': enames,
        'lvec': lvec,
        'z_top': z_top,
        'n_atoms': len(apos)
    }


def infer_grid_metadata(grid_path, substrate_info):
    """
    Infer grid origin (g0) and spacing (dg) from grid shape and substrate lattice.
    
    Args:
        grid_path: Path to GridFF .npy file
        substrate_info: Dict from load_substrate_xyz_with_lvec()
        
    Returns:
        dict with:
            'ns': (nx, ny, nz) grid shape
            'dg': (dx, dy, dz) grid spacing
            'g0': (x0, y0, z0) grid origin
            'Ls': (Lx, Ly, Lz) lattice dimensions
            'convention': str describing detected convention
    """
    grid = np.load(grid_path)
    ns = grid.shape[:3]  # (nx, ny, nz) - note: grid is stored as (nx, ny, nz, nch)
    
    lvec = substrate_info['lvec']
    Lx = float(np.linalg.norm(lvec[0]))
    Ly = float(np.linalg.norm(lvec[1]))
    Lz = float(np.linalg.norm(lvec[2]))
    
    # Grid spacing
    dg = (Lx / ns[0], Ly / ns[1], Lz / ns[2])
    
    # Try different origin conventions
    z_top = substrate_info['z_top']
    z_min_atom = np.min(substrate_info['apos'][:, 2])
    
    conventions = {
        'centered_xy_bottom_z': (-Lx/2, -Ly/2, z_min_atom - 2.0),  # Centered XY, z below atoms
        'centered_xy_zero_z': (-Lx/2, -Ly/2, 0.0),  # Centered XY, z at 0
        'corner_xy_bottom_z': (0.0, 0.0, z_min_atom - 2.0),  # Corner XY, z below atoms
        'corner_xy_zero_z': (0.0, 0.0, 0.0),  # Corner XY, z at 0
    }
    
    return {
        'ns': ns,
        'dg': dg,
        'Ls': (Lx, Ly, Lz),
        'lvec': lvec,
        'z_top': z_top,
        'conventions': conventions,
        'grid_shape_full': grid.shape
    }


# =============================================================================
# Section 2: Visualization Utilities (adapted from RigidBodyAFM.py)
# =============================================================================

def plot_gridff_diagnostics(grid_data, sub_apos, sub_enames, lvec, iz_slices=None, iy_slice=None, save_path='grid_diagnostics.png', g0=None, dg=None, z_marks=None, channel_name=None, z_atom_range=5.0):
    """
    Diagnostic tool to plot GridFF channels and overlay substrate atoms.
    Adapted from RigidBodyAFM.py with added g0/dg support.
    
    Args:
        grid_data: (nx, ny, nz, nch) GridFF array
        sub_apos: (n, 3) substrate atom positions
        sub_enames: List of element names
        lvec: (3, 3) lattice vectors
        iz_slices: List of z-indices to plot for XY slices (default: [nz//4, nz//2, 3*nz//4])
        iy_slice: Y-index for XZ slice (default: ny//2)
        save_path: Output PNG path
        g0: Grid origin (x0, y0, z0) - if None, assumes (0,0,0)
        dg: Grid spacing (dx, dy, dz) - if None, infers from substrate lattice
        z_marks: List of z-values to mark with dashed lines on XZ subplot
        channel_name: Override channel name for single-channel plots (e.g., 'total', 'pauli_vdw', 'electrostatic')
    """
    nx, ny, nz, nch = grid_data.shape
    
    # Compute grid parameters
    if dg is None:
        ax = float(np.linalg.norm(lvec[0]))
        ay = float(np.linalg.norm(lvec[1]))
        az = float(np.linalg.norm(lvec[2]))
        dg = (ax/nx, ay/ny, az/nz)
    else:
        ax, ay, az = nx*dg[0], ny*dg[1], nz*dg[2]
    
    if g0 is None:
        g0 = (0.0, 0.0, 0.0)
    
    dx, dy, dz = dg
    
    # Color mapping for atoms
    colors = ['purple' if e in ['Ca', 'Na'] else 'green' if e in ['F', 'Cl'] else 'gray' 
              for e in sub_enames]
    
    # Default slice indices
    if iz_slices is None:
        iz_slices = [nz//4, nz//2, 3*nz//4]
    if iy_slice is None:
        iy_slice = ny // 2
    
    n_slices = len(iz_slices)
    fig, axs = plt.subplots(nch, n_slices + 1, figsize=(5*(n_slices+1), 4*nch))
    if nch == 1: 
        axs = axs[None, :]
    if n_slices == 0:
        axs = axs[:, None]

    if channel_name is not None:
        names = [channel_name]
    else:
        names = ['Pauli (P)', 'London (L)', 'Electrostatic (Q)', 'Hydrogen (H)'][:nch]

    # Show atoms within z_atom_range below the top atom on XY plots
    z_top_atom = np.max(sub_apos[:, 2])
    xy_atom_mask = sub_apos[:, 2] >= (z_top_atom - z_atom_range)
    
    for i in range(nch):
        # XY Slices
        for j, iz in enumerate(iz_slices):
            z_val = g0[2] + iz * dz
            extent = [g0[0], g0[0] + ax, g0[1], g0[1] + ay]
            im = axs[i, j].imshow(grid_data[:, :, iz, i].T, extent=extent, 
                                   origin='lower', cmap='bwr', aspect='equal')
            # Overlay all surface atoms (within z_atom_range below top)
            if np.any(xy_atom_mask):
                axs[i, j].scatter(sub_apos[xy_atom_mask, 0], sub_apos[xy_atom_mask, 1], c=np.array(colors)[xy_atom_mask], s=40, alpha=0.7, edgecolors='black')
            axs[i, j].set_title(f"{names[i]} XY at z={z_val:.2f} A (iz={iz})")
            axs[i, j].set_xlabel('x [A]')
            axs[i, j].set_ylabel('y [A]')
            plt.colorbar(im, ax=axs[i, j])

        # XZ Slice
        y_val = g0[1] + iy_slice * dy
        extent_xz = [g0[0], g0[0] + ax, g0[2], g0[2] + az]
        xz_data = grid_data[:, iy_slice, :, i].T
        # Compute symmetric vmin/vmax from data above first XY slice to avoid saturation
        if iz_slices is not None and len(iz_slices) > 0:
            iz_min = iz_slices[0]  # First XY slice index
            # Only consider data for iz > iz_min (exclude the slice itself)
            xz_data_above = xz_data[iz_min+1:, :]  # xz_data is (nz, nx) after .T — z is axis 0
            max_abs = max(abs(xz_data_above.min()), abs(xz_data_above.max()))
            if max_abs > 0:
                vmin, vmax = -max_abs, max_abs
            else:
                vmin, vmax = None, None
        else:
            max_abs = max(abs(xz_data.min()), abs(xz_data.max()))
            if max_abs > 0:
                vmin, vmax = -max_abs, max_abs
            else:
                vmin, vmax = None, None
        im_xz = axs[i, -1].imshow(xz_data, extent=extent_xz, origin='lower', cmap='bwr', aspect='equal', vmin=vmin, vmax=vmax)
        mask_y = np.abs(sub_apos[:, 1] - y_val) < 2.0
        if np.any(mask_y):
            axs[i, -1].scatter(sub_apos[mask_y, 0], sub_apos[mask_y, 2],  c=np.array(colors)[mask_y], s=20, alpha=0.5, edgecolors='black')
        if z_marks is not None:
            for zm in z_marks:
                axs[i, -1].axhline(zm, color='k', linestyle='--', linewidth=1.0, alpha=0.7, label=f'z={zm:.2f}')
        axs[i, -1].set_title(f"{names[i]} XZ at y={y_val:.2f} A (iy={iy_slice})")
        axs[i, -1].set_xlabel('x [A]')
        axs[i, -1].set_ylabel('z [A]')
        axs[i, -1].set_ylim(g0[2], g0[2] + min(az, 20))  # Limit z range for visibility
        plt.colorbar(im_xz, ax=axs[i, -1])

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved GridFF diagnostics to {save_path}")
    return save_path


def plot_alignment_summary(grid_data, g0, dg, atoms_xyz, atoms_enames, save_path, iz_top=None, iy_center=None):
    """
    Generate comprehensive alignment diagnostic figure.
    
    Args:
        grid_data: (nx, ny, nz, nch) GridFF array
        g0: Grid origin (x0, y0, z0)
        dg: Grid spacing (dx, dy, dz)
        atoms_xyz: (n, 3) atom positions
        atoms_enames: Element names for color coding
        save_path: Output PNG path
        iz_top: Z-index for top layer (default: auto-detect from atoms)
        iy_center: Y-index for center slice (default: ny//2)
    """
    nx, ny, nz, nch = grid_data.shape
    dx, dy, dz = dg
    
    # Auto-detect top layer if not specified
    if iz_top is None:
        z_top = np.max(atoms_xyz[:, 2])
        iz_top = int((z_top - g0[2]) / dz)
        iz_top = max(0, min(nz-1, iz_top))
    
    if iy_center is None:
        iy_center = ny // 2
    
    # Color mapping
    colors = ['purple' if e in ['Ca', 'Na'] else 'green' if e in ['F', 'Cl'] else 'gray' 
              for e in atoms_enames]
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # 1. XY slice at top layer - Pauli
    ax = axes[0, 0]
    z_val = g0[2] + iz_top * dz
    extent = [g0[0], g0[0] + nx*dx, g0[1], g0[1] + ny*dy]
    im = ax.imshow(grid_data[:, :, iz_top, 0].T, extent=extent, 
                   origin='lower', cmap='bwr', aspect='equal')
    ax.scatter(atoms_xyz[:, 0], atoms_xyz[:, 1], c=colors, s=30, alpha=0.6, 
               edgecolors='black', linewidths=1)
    ax.set_title(f'Pauli Potential XY at z={z_val:.2f}A (iz={iz_top})')
    ax.set_xlabel('x [A]')
    ax.set_ylabel('y [A]')
    plt.colorbar(im, ax=ax)
    
    # 2. XZ slice at center Y - Pauli
    ax = axes[0, 1]
    y_val = g0[1] + iy_center * dy
    extent_xz = [g0[0], g0[0] + nx*dx, g0[2], g0[2] + nz*dz]
    im = ax.imshow(grid_data[:, iy_center, :, 0].T, extent=extent_xz,
                   origin='lower', cmap='bwr', aspect='auto')
    mask_y = np.abs(atoms_xyz[:, 1] - y_val) < 2.0
    if np.any(mask_y):
        ax.scatter(atoms_xyz[mask_y, 0], atoms_xyz[mask_y, 2], 
                   c=np.array(colors)[mask_y], s=30, alpha=0.6, edgecolors='black')
    ax.set_title(f'Pauli Potential XZ at y={y_val:.2f}A (iy={iy_center})')
    ax.set_xlabel('x [A]')
    ax.set_ylabel('z [A]')
    plt.colorbar(im, ax=ax)
    
    # 3. 1D profiles through atom centers
    ax = axes[1, 0]
    # Find center atom
    cx = nx // 2
    cy = ny // 2
    
    # X-profile at center Y, top Z
    x_coords = g0[0] + np.arange(nx) * dx
    ax.plot(x_coords, grid_data[:, cy, iz_top, 0], 'b-', label='Pauli X-profile')
    ax.plot(x_coords, grid_data[:, cy, iz_top, 1], 'r-', label='London X-profile')
    ax.axhline(0, color='k', linestyle='--', alpha=0.3)
    ax.set_xlabel('x [A]')
    ax.set_ylabel('Energy [eV]')
    ax.set_title(f'X-Profile at y={g0[1]+cy*dy:.2f}A, z={z_val:.2f}A')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 4. Z-profile at center
    ax = axes[1, 1]
    z_coords = g0[2] + np.arange(nz) * dz
    ax.plot(z_coords, grid_data[cx, cy, :, 0], 'b-', label='Pauli')
    ax.plot(z_coords, grid_data[cx, cy, :, 1], 'r-', label='London')
    ax.axhline(0, color='k', linestyle='--', alpha=0.3)
    ax.axvline(z_val, color='g', linestyle='--', alpha=0.5, label=f'z_top={z_val:.2f}A')
    ax.set_xlabel('z [A]')
    ax.set_ylabel('Energy [eV]')
    ax.set_title(f'Z-Profile at x={g0[0]+cx*dx:.2f}A, y={g0[1]+cy*dy:.2f}A')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_ylim(-5, 10)  # Focus on relevant range
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved alignment summary to {save_path}")
    return save_path




def sample_gridff_opencl(gridff_path, sub_xyz, positions, grid_p0=(0.0, 0.0, 0.0), grid_step=None,
                         atom_req=(1.487, 0.0006808, 0.0, 0.0), 
                         atom_mass=1.008, alpha_morse=1.5, debug=False):
    """
    Sample GridFF at arbitrary positions using OpenCL B-spline interpolation.
    This uses the same implementation as molecular simulations.
    
    Wraps RigidBodyAFM.sample_gridff_single_atom().
    
    Args:
        gridff_path: Path to GridFF .npy file
        sub_xyz: Path to substrate .xyz file (for lattice vectors)
        positions: (n, 3) positions to sample
        grid_p0: Grid origin (x0, y0, z0)
        grid_step: Grid spacing (dx, dy, dz)
        atom_req: Tuple (R, E, Q, H) for test atom
        atom_mass: Mass of test atom
        alpha_morse: Alpha parameter for REQ->PLQ conversion
        debug: Enable debug output
        
    Returns:
        dict with 'forces' (n,3), 'energies' (n,)
    """
    forces, energies = sample_gridff_single_atom(
        scan_positions=positions,
        gridff_path=gridff_path,
        sub_xyz=sub_xyz,
        atom_req=atom_req,
        atom_mass=atom_mass,
        alpha_morse=alpha_morse,
        debug=debug,
        grid_p0=grid_p0,
        grid_step=grid_step
    )
    
    return {'forces': forces, 'energies': energies}


def compare_sampling_methods(grid_data, g0, dg, gridff_path, sub_xyz, positions, 
                             atom_req, alpha_morse, verbose=True):
    """
    Compare simple coefficient summation vs OpenCL B-spline sampling for 3 components:
    1. Total potential (all PLQ)
    2. Pauli+vdW (Q=0)
    3. Electrostatics (E=0)
    
    This verifies that the grid transformation is consistent between methods.
    They will not give exactly the same values (different interpolation),
    but should show similar trends and spatial patterns.
    
    Args:
        grid_data: (nx, ny, nz, nch) array
        g0: Grid origin
        dg: Grid spacing
        gridff_path: Path to GridFF .npy file
        sub_xyz: Path to substrate .xyz file
        positions: (n, 3) positions to sample
        atom_req: Tuple (R, E, Q, H) for test atom
        alpha_morse: Alpha parameter
        verbose: Print comparison statistics
        
    Returns:
        dict with comparison statistics for each component
    """
    R, E, Q, H = atom_req
    
    # Define 3 components
    # Note: Grid has 3 channels (0=Pauli, 1=London, 2=Coulomb)
    # Channel 3 (Hydrogen) is not present in this grid
    # For fair comparison, use atom_req that gives unit coefficients for the channels we want
    components = {
        'total': {
            'name': 'Total (PLQ)',
            'channels': [0, 1, 2],  # Pauli, London, Coulomb
            'atom_req': (R, E, Q, H)
        },
        'pauli_vdw': {
            'name': 'Pauli+vdW (Q=0)',
            'channels': [0, 1],  # Pauli, London only (no Hydrogen channel)
            'atom_req': (R, E, 0.0, H)  # Use same REQ as original
        },
        'electrostatic': {
            'name': 'Electrostatic (E=0, Q=1)',
            'channels': [2],  # Coulomb only
            'atom_req': (R, 0.0, 1.0, 0.0)  # Set Q=1 to sample electrostatics
        }
    }
    
    results = {}
    for comp_key, comp in components.items():
        # Simple: sum numpy arrays at nearest grid points (no interpolation)
        # Convert positions to grid indices
        nx, ny, nz, nch = grid_data.shape
        dx, dy, dz = dg
        x0, y0, z0 = g0
        
        positions_arr = np.asarray(positions, dtype=np.float32)
        ix = ((positions_arr[:, 0] - x0) / dx).astype(np.int32)
        iy = ((positions_arr[:, 1] - y0) / dy).astype(np.int32)
        iz = ((positions_arr[:, 2] - z0) / dz).astype(np.int32)
        
        # Clamp to valid range
        ix = np.clip(ix, 0, nx-1)
        iy = np.clip(iy, 0, ny-1)
        iz = np.clip(iz, 0, nz-1)
        
        # Compute PLQ coefficients for this component's atom_req
        from .RigidBodyDynamics import _reqs_to_plq
        comp_req = np.array([comp['atom_req']], dtype=np.float32)
        comp_plq = _reqs_to_plq(comp_req, alpha=alpha_morse)[0]
        
        # Sum channels at nearest grid points, weighted by PLQ coefficients
        # For electrostatic (channel 2 only), use Q coefficient directly
        simple_vals = np.zeros(len(positions), dtype=np.float32)
        if comp_key == 'electrostatic':
            # For electrostatic, use Q coefficient (index 2 in PLQ)
            simple_vals += grid_data[ix, iy, iz, 2] * comp_plq[2]
        else:
            for i, ch in enumerate(comp['channels']):
                simple_vals += grid_data[ix, iy, iz, ch] * comp_plq[i]
        
        # OpenCL: sample with B-spline interpolation
        opencl_result = sample_gridff_opencl(gridff_path, sub_xyz, positions, grid_p0=g0, grid_step=dg,
                                              atom_req=comp['atom_req'], alpha_morse=alpha_morse, debug=False)
        opencl_vals = opencl_result['energies']
        
        # Compute statistics
        corr = np.corrcoef(simple_vals, opencl_vals)[0, 1] if len(simple_vals) > 1 else np.nan
        rmse = np.sqrt(np.mean((simple_vals - opencl_vals)**2))
        
        results[comp_key] = {
            'name': comp['name'],
            'simple': simple_vals,
            'opencl': opencl_vals,
            'correlation': corr,
            'rmse': rmse
        }
        
        if verbose:
            print(f"  {comp['name']}:")
            print(f"    Simple: mean={np.mean(simple_vals):.4f}, std={np.std(simple_vals):.4f}")
            print(f"    OpenCL: mean={np.mean(opencl_vals):.4f}, std={np.std(opencl_vals):.4f}")
            print(f"    Correlation: {corr:.4f}, RMSE: {rmse:.4f}")
    
    return results


def sample_grid_at_atoms_opencl(gridff_path, sub_xyz, atoms_xyz, atom_req, alpha_morse, 
                                  grid_p0=(0.0, 0.0, 0.0), verbose=True):
    """
    Sample GridFF at atom positions using OpenCL B-spline interpolation.
    
    Args:
        gridff_path: Path to GridFF .npy file
        sub_xyz: Path to substrate .xyz file
        atoms_xyz: (n, 3) atom positions
        grid_p0: Grid origin (x0, y0, z0)
        atom_req: Tuple (R, E, Q, H) for test atom
        alpha_morse: Alpha parameter
        verbose: Print statistics
        
    Returns:
        dict with 'forces', 'energies', 'stats'
    """
    result = sample_gridff_opencl(gridff_path, sub_xyz, atoms_xyz, grid_p0=grid_p0, grid_step=None,
                                   atom_req=atom_req, alpha_morse=alpha_morse, debug=False)
    
    energies = result['energies']
    stats = {
        'mean': float(np.mean(energies)),
        'min': float(np.min(energies)),
        'max': float(np.max(energies)),
        'std': float(np.std(energies))
    }
    
    if verbose:
        print(f"Sampling GridFF (OpenCL) at {len(atoms_xyz)} atom positions:")
        print(f"  Mean: {stats['mean']:.4f} eV")
        print(f"  Min:  {stats['min']:.4f} eV")
        print(f"  Max:  {stats['max']:.4f} eV")
        print(f"  Std:  {stats['std']:.4f} eV")
    
    return {'forces': result['forces'], 'energies': energies, 'stats': stats}


# =============================================================================
# Section 4: Alignment Verification
# =============================================================================

def find_grid_minima(grid_data, g0, dg, component='london', threshold=-0.5):
    """
    Find local minima in grid using scipy.ndimage.minimum_filter.
    
    Args:
        grid_data: GridFF array
        g0: Grid origin
        dg: Grid spacing
        component: 'pauli' (0), 'london' (1), or 'coulomb' (2)
        threshold: Only return minima below this value
        
    Returns:
        minima_xyz: (n, 3) positions of minima
    """
    from scipy.ndimage import minimum_filter
    
    ch_map = {'pauli': 0, 'london': 1, 'coulomb': 2, 'total': 3}
    ch = ch_map.get(component, 1)
    
    if ch >= grid_data.shape[3]:
        print(f"Warning: Channel {component} not available, using channel 0")
        ch = 0
    
    grid_ch = grid_data[:, :, :, ch]
    
    # Find local minima
    footprint = np.ones((3, 3, 3))
    local_min = minimum_filter(grid_ch, footprint=footprint, mode='constant') == grid_ch
    
    # Filter by threshold
    mask = (grid_ch < threshold) & local_min
    
    if not np.any(mask):
        print(f"Warning: No minima found below threshold {threshold}")
        # Return empty array
        return np.zeros((0, 3))
    
    # Get indices
    indices = np.argwhere(mask)
    
    # Convert to physical coordinates
    dx, dy, dz = dg
    x0, y0, z0 = g0
    
    minima_xyz = np.zeros_like(indices, dtype=np.float32)
    minima_xyz[:, 0] = x0 + indices[:, 0] * dx
    minima_xyz[:, 1] = y0 + indices[:, 1] * dy
    minima_xyz[:, 2] = z0 + indices[:, 2] * dz
    
    return minima_xyz


def verify_atom_grid_alignment(atoms_xyz, grid_minima_xyz, threshold=0.2, verbose=True):
    """
    Verify alignment between atoms and grid minima.
    
    Args:
        atoms_xyz: (n, 3) atom positions
        grid_minima_xyz: (m, 3) detected grid minima
        threshold: Maximum allowed distance error in Angstroms
        verbose: Print statistics
        
    Returns:
        dict with alignment statistics
    """
    if len(grid_minima_xyz) == 0:
        return {
            'n_atoms': len(atoms_xyz),
            'n_minima': 0,
            'mean_error': float('inf'),
            'max_error': float('inf'),
            'aligned': False,
            'error': 'No grid minima found'
        }
    
    # Compute distances from each atom to nearest minimum
    discrepancies = []
    for atom in atoms_xyz:
        distances = np.linalg.norm(grid_minima_xyz - atom, axis=1)
        discrepancies.append(np.min(distances))
    
    discrepancies = np.array(discrepancies)
    
    stats = {
        'n_atoms': len(atoms_xyz),
        'n_minima': len(grid_minima_xyz),
        'mean_error': float(np.mean(discrepancies)),
        'max_error': float(np.max(discrepancies)),
        'min_error': float(np.min(discrepancies)),
        'std_error': float(np.std(discrepancies)),
        'aligned': float(np.max(discrepancies)) < threshold
    }
    
    if verbose:
        print(f"\nAlignment Verification:")
        print(f"  Atoms: {stats['n_atoms']}, Grid minima: {stats['n_minima']}")
        print(f"  Mean distance: {stats['mean_error']:.4f} Å")
        print(f"  Max distance:  {stats['max_error']:.4f} Å")
        print(f"  Std distance:  {stats['std_error']:.4f} Å")
        print(f"  Threshold:     {threshold:.4f} Å")
        print(f"  Aligned:       {stats['aligned']}")
    
    return stats


def test_shift_convention(grid_data, atoms_xyz, g0, dg, convention_name):
    """
    Test a specific shift convention by sampling at atom positions.
    
    Returns error metric (lower is better alignment).
    """
    # For Pauli potential, values at atom positions should be strongly negative (repulsive core)
    # Sample Pauli channel at nearest grid points (no interpolation)
    nx, ny, nz, nch = grid_data.shape
    dx, dy, dz = dg
    x0, y0, z0 = g0
    
    positions_arr = np.asarray(atoms_xyz, dtype=np.float32)
    ix = ((positions_arr[:, 0] - x0) / dx).astype(np.int32)
    iy = ((positions_arr[:, 1] - y0) / dy).astype(np.int32)
    iz = ((positions_arr[:, 2] - z0) / dz).astype(np.int32)
    
    # Clamp to valid range
    ix = np.clip(ix, 0, nx-1)
    iy = np.clip(iy, 0, ny-1)
    iz = np.clip(iz, 0, nz-1)
    
    # Get Pauli values at nearest grid points
    values = grid_data[ix, iy, iz, 0]
    
    # We want values to be minimum (most negative) at atom centers
    # Metric: mean value (more negative = better alignment)
    metric = -np.mean(values)  # Negative because we want to maximize negativity
    
    return metric


def auto_detect_shift(grid_data, atoms_xyz, substrate_lvec, verbose=True):
    """
    Try different shift conventions and find best alignment.
    
    Returns:
        dict with best convention and error metrics
    """
    nx, ny, nz, nch = grid_data.shape
    
    # Compute lattice dimensions
    Lx = float(np.linalg.norm(substrate_lvec[0]))
    Ly = float(np.linalg.norm(substrate_lvec[1]))
    Lz = float(np.linalg.norm(substrate_lvec[2]))
    
    dg = (Lx/nx, Ly/ny, Lz/nz)
    
    # Try different conventions
    z_atoms = atoms_xyz[:, 2]
    z_min_atom = np.min(z_atoms)
    z_max_atom = np.max(z_atoms)
    
    conventions = {
        'centered_xy_bottom_z': (-Lx/2, -Ly/2, z_min_atom - 2.0),  # Centered XY, z below atoms
        'centered_xy_zero_z': (-Lx/2, -Ly/2, 0.0),  # Centered XY, z at 0
        'corner_xy_bottom_z': (0.0, 0.0, z_min_atom - 2.0),  # Corner XY, z below atoms
        'corner_xy_zero_z': (0.0, 0.0, 0.0),  # Corner XY, z at 0
    }
    
    results = {}
    for name, g0_test in conventions.items():
        metric = test_shift_convention(grid_data, atoms_xyz, g0_test, dg, name)
        results[name] = {'g0': g0_test, 'metric': metric}
    
    # Find best (highest metric = most negative mean = best alignment)
    best_name = max(results.keys(), key=lambda k: results[k]['metric'])
    best = results[best_name]
    
    if verbose:
        print(f"\nShift Convention Detection:")
        for name, res in results.items():
            marker = " <-- BEST" if name == best_name else ""
            print(f"  {name:25s}: metric={res['metric']:.4f}{marker}")
        print(f"\nBest convention: {best_name}")
        print(f"  g0 = {best['g0']}")
    
    return {
        'best_convention': best_name,
        'best_g0': best['g0'],
        'dg': dg,
        'all_results': results
    }


# =============================================================================
# Section 5: Main Orchestration Function
# =============================================================================

def run_alignment_verification(grid_path, substrate_path, save_dir, 
                               test_conventions=None, atom_req=(1.487, 0.0006808, 0.0, 0.0), 
                               alpha_morse=1.5, verbose=True):
    """
    Run complete alignment verification workflow.
    
    Args:
        grid_path: Path to GridFF .npy file
        substrate_path: Path to substrate .xyz file
        save_dir: Directory for output plots and reports
        test_conventions: List of convention names to test (None = test all)
        atom_req: Tuple (R, E, Q, H) for test atom
        alpha_morse: Alpha parameter for REQ to PLQ conversion
        verbose: Print progress
        
    Returns:
        dict with verification results
    """
    os.makedirs(save_dir, exist_ok=True)
    
    base_name = os.path.splitext(os.path.basename(grid_path))[0]
    
    if verbose:
        print(f"\n{'='*60}")
        print(f"GridFF Alignment Verification")
        print(f"{'='*60}")
        print(f"Grid:    {grid_path}")
        print(f"Substrate: {substrate_path}")
        print(f"Output:  {save_dir}")
    
    # Load data
    if verbose:
        print(f"\n[1/5] Loading GridFF and substrate...")
    grid_data = load_gridff_array(grid_path)
    sub_info = load_substrate_xyz_with_lvec(substrate_path)
    
    if verbose:
        print(f"  Grid shape: {grid_data.shape}")
        print(f"  Substrate atoms: {sub_info['n_atoms']}")
        print(f"  Lattice vectors:")
        print(f"    a = {sub_info['lvec'][0]}")
        print(f"    b = {sub_info['lvec'][1]}")
        print(f"    c = {sub_info['lvec'][2]}")
    
    # Try to load metadata from JSON first, fall back to inference
    if verbose:
        print(f"\n[2/5] Loading grid metadata...")
    metadata_json = load_gridff_metadata(grid_path)
    
    if metadata_json is not None:
        # Use metadata from JSON
        metadata = {
            'ns': tuple(metadata_json['ns']),
            'dg': tuple(metadata_json['dg']),
            'g0': tuple(metadata_json['g0']),
            'Ls': tuple([np.linalg.norm(sub_info['lvec'][i]) for i in range(3)]),
            'z_top': float(metadata_json['z0']),
            'conventions': {
                'from_json': tuple(metadata_json['g0'])
            }
        }
        if verbose:
            print(f"  Using metadata from JSON file")
            print(f"  Grid shape (nx, ny, nz): {metadata['ns']}")
            print(f"  Grid spacing (dx, dy, dz): {metadata['dg']}")
            print(f"  Grid origin (g0): {metadata['g0']}")
            print(f"  Lattice (Lx, Ly, Lz): {metadata['Ls']}")
            print(f"  z_top (substrate): {metadata['z_top']:.2f} Å")
    else:
        # Fall back to inference
        if verbose:
            print(f"  JSON metadata not found, inferring from grid and substrate...")
        metadata = infer_grid_metadata(grid_path, sub_info)
        if verbose:
            print(f"  Grid shape (nx, ny, nz): {metadata['ns']}")
            print(f"  Grid spacing (dx, dy, dz): {metadata['dg']}")
            print(f"  Lattice (Lx, Ly, Lz): {metadata['Ls']}")
            print(f"  z_top (substrate): {metadata['z_top']:.2f} Å")
    
    # Generate diagnostic plots with different conventions
    if verbose:
        print(f"\n[3/5] Generating diagnostic plots...")
    
    plot_results = []
    conventions_to_test = test_conventions or list(metadata['conventions'].keys())
    
    for conv_name in conventions_to_test:
        g0 = metadata['conventions'][conv_name]
        dg = metadata['dg']
        
        save_path = os.path.join(save_dir, f"{base_name}_diagnostic_{conv_name}.png")
        
        try:
            plot_gridff_diagnostics(
                grid_data, sub_info['apos'], sub_info['enames'], sub_info['lvec'],
                iz_slices=[metadata['ns'][2]//4, metadata['ns'][2]//2, 3*metadata['ns'][2]//4],
                iy_slice=metadata['ns'][1]//2,
                save_path=save_path,
                g0=g0, dg=dg
            )
            plot_results.append({'convention': conv_name, 'path': save_path, 'status': 'ok'})
        except Exception as e:
            plot_results.append({'convention': conv_name, 'path': save_path, 'status': 'error', 'error': str(e)})
            if verbose:
                print(f"  Warning: Failed to plot {conv_name}: {e}")
    
    # Generate alignment summary for best convention
    if verbose:
        print(f"\n[4/5] Testing shift conventions...")
    
    # If metadata from JSON is available, use it directly
    if metadata_json is not None:
        best_g0 = metadata['g0']
        best_dg = metadata['dg']
        shift_detection = {
            'best_convention': 'from_json',
            'best_g0': best_g0,
            'dg': best_dg,
            'all_results': {}
        }
        if verbose:
            print(f"  Using grid origin from metadata: g0 = {best_g0}")
    else:
        # Fall back to auto-detection
        shift_detection = auto_detect_shift(
            grid_data, sub_info['apos'], sub_info['lvec'], verbose=verbose
        )
        best_g0 = shift_detection['best_g0']
        best_dg = shift_detection['dg']
    
    summary_path = os.path.join(save_dir, f"{base_name}_alignment_summary.png")
    plot_alignment_summary(
        grid_data, best_g0, best_dg,
        sub_info['apos'], sub_info['enames'],
        summary_path
    )
    
    # Sample at positions above substrate (z0+2.0, z0+6.0) with best convention
    if verbose:
        print(f"\n[5/5] Sampling GridFF at positions above substrate (z0+2.0, z0+6.0)...")
    
    # Generate test positions at z0+2.0 and z0+6.0 above substrate
    # Use same (x,y) as substrate atoms, but different z heights relative to top atom
    z0 = metadata['z_top']  # z of the topmost substrate atom
    test_positions = []
    for z_offset in [2.0, 6.0]:
        z_height = z0 + z_offset
        for pos in sub_info['apos']:
            test_positions.append([pos[0], pos[1], z_height])
    test_positions = np.array(test_positions, dtype=np.float32)
    
    # Use OpenCL sampling (same as molecular simulations)
    # Pass the correct grid origin from metadata or auto-detection
    sampling_results = sample_grid_at_atoms_opencl(
        grid_path, substrate_path, test_positions,
        grid_p0=best_g0,  # Use grid origin from metadata or auto-detection
        atom_req=atom_req,
        alpha_morse=alpha_morse,
        verbose=verbose
    )
    
    # Compare with simple coefficient summation for verification
    if verbose:
        print(f"\n[5/6] Comparing sampling methods (OpenCL vs simple)...")
        comparison = compare_sampling_methods(
            grid_data, best_g0, best_dg, grid_path, substrate_path,
            test_positions,  # Use positions above substrate
            atom_req=atom_req,
            alpha_morse=alpha_morse,
            verbose=True
        )
        
        # Generate comparison plots for 3 components
        print(f"\n[6/6] Generating sampling comparison plots...")
        
        R, E, Q, H = atom_req
        components = {
            'total': {
                'name': 'total',
                'channels': [0, 1, 2],  # Grid has 3 channels
                'atom_req': (R, E, Q, H)
            },
            'pauli_vdw': {
                'name': 'pauli_vdw',
                'channels': [0, 1],  # Pauli, London only
                'atom_req': (R, E, 0.0, H)
            },
            'electrostatic': {
                'name': 'electrostatic',
                'channels': [2],  # Coulomb only
                'atom_req': (R, 0.0, 1.0, 0.0)  # Set Q=1 to sample electrostatics
            }
        }
        
        nx, ny, nz, nch = grid_data.shape
        ax, ay, az = nx*best_dg[0], ny*best_dg[1], nz*best_dg[2]
        x_sample = np.linspace(best_g0[0], best_g0[0] + ax, nx)
        y_sample = np.linspace(best_g0[1], best_g0[1] + ay, ny)
        z_sample = np.linspace(best_g0[2], best_g0[2] + az, nz)
        # Compute iz_slices at z0+2.0 and z0+6.0 above top surface atom
        z0_surf = metadata['z_top']  # z of the topmost substrate atom
        iz_z2 = int((z0_surf + 2.0 - best_g0[2]) / best_dg[2])
        iz_z6 = int((z0_surf + 6.0 - best_g0[2]) / best_dg[2])
        iz_slices = [iz_z2, iz_z6]
        iy_slice = metadata['ns'][1]//2
        if verbose:
            print(f"  XY slice heights: iz={iz_z2} (z={z0_surf+2.0:.2f}A), iz={iz_z6} (z={z0_surf+6.0:.2f}A)")
        
        # z-values of the XY slices for marking on XZ subplot
        z_marks = [best_g0[2] + iz * best_dg[2] for iz in iz_slices]
        
        # Compute PLQ coefficients for each component (same as GPU)
        from .RigidBodyDynamics import _reqs_to_plq
        
        component_plots = {}
        for comp_key, comp in components.items():
            # Simple plot: weight grid channels by PLQ coefficients (same as GPU kernel)
            comp_req = np.array([comp['atom_req']], dtype=np.float32)
            comp_plq = _reqs_to_plq(comp_req, alpha=alpha_morse)[0]  # [cP, cL, Q, cH]
            plq_weights = [comp_plq[ch] for ch in comp['channels']]  # weight for each channel
            simple_grid = np.zeros((nx, ny, nz, 1), dtype=np.float32)
            for ch, w in zip(comp['channels'], plq_weights):
                simple_grid[:, :, :, 0] += grid_data[:, :, :, ch] * w
            
            simple_plot_path = os.path.join(save_dir, f"{base_name}_diagnostic_{comp['name']}_simple.png")
            plot_gridff_diagnostics(
                simple_grid, sub_info['apos'], sub_info['enames'], sub_info['lvec'],
                iz_slices=iz_slices,
                iy_slice=iy_slice,
                save_path=simple_plot_path,
                g0=best_g0, dg=best_dg, z_marks=z_marks, channel_name=comp['name']
            )
            
            # OpenCL plot: sample with B-spline
            opencl_xy_slices = {}
            for iz in iz_slices:
                z_val = best_g0[2] + iz * best_dg[2]
                positions = []
                for ix in range(nx):
                    for iy in range(ny):
                        positions.append([x_sample[ix], y_sample[iy], z_val])
                positions = np.array(positions, dtype=np.float32)
                
                result = sample_gridff_opencl(grid_path, substrate_path, positions, grid_p0=best_g0, grid_step=best_dg,
                                              atom_req=comp['atom_req'], alpha_morse=alpha_morse, debug=False)
                energies = result['energies'].reshape(nx, ny)
                opencl_xy_slices[iz] = energies
            
            y_val = best_g0[1] + iy_slice * best_dg[1]
            positions_xz = []
            # Sample only valid B-spline region (iz from 2 to nz-2) to avoid zeros at boundaries
            iz_min = 2
            iz_max = nz - 2
            for ix in range(nx):
                for iz in range(iz_min, iz_max):
                    positions_xz.append([x_sample[ix], y_val, z_sample[iz]])
            positions_xz = np.array(positions_xz, dtype=np.float32)
            
            result_xz = sample_gridff_opencl(grid_path, substrate_path, positions_xz, grid_p0=best_g0, grid_step=best_dg,
                                               atom_req=comp['atom_req'], alpha_morse=alpha_morse, debug=False)
            # Reshape to (nx, nz_valid)
            nz_valid = iz_max - iz_min
            opencl_xz = result_xz['energies'].reshape(nx, nz_valid)
            
            # Pad with zeros to match full grid size for plotting
            opencl_xz_full = np.zeros((nx, nz), dtype=np.float32)
            opencl_xz_full[:, iz_min:iz_max] = opencl_xz
            opencl_xz = opencl_xz_full
            
            opencl_grid = np.zeros((nx, ny, nz, 1), dtype=np.float32)
            for iz in iz_slices:
                opencl_grid[:, :, iz, 0] = opencl_xy_slices[iz]
            opencl_grid[:, iy_slice, :, 0] = opencl_xz
            
            opencl_plot_path = os.path.join(save_dir, f"{base_name}_diagnostic_{comp['name']}_opencl.png")
            plot_gridff_diagnostics(
                opencl_grid, sub_info['apos'], sub_info['enames'], sub_info['lvec'],
                iz_slices=iz_slices,
                iy_slice=iy_slice,
                save_path=opencl_plot_path,
                g0=best_g0, dg=best_dg, z_marks=z_marks, channel_name=comp['name']
            )
            
            component_plots[comp_key] = {
                'simple': simple_plot_path,
                'opencl': opencl_plot_path
            }
    
    # Verify alignment with grid minima
    minima_xyz = find_grid_minima(grid_data, best_g0, best_dg, component='london', threshold=-0.5)
    alignment_stats = verify_atom_grid_alignment(
        sub_info['apos'], minima_xyz, threshold=0.5, verbose=verbose
    )
    
    # Save JSON report
    report = {
        'grid_path': grid_path,
        'substrate_path': substrate_path,
        'grid_shape': grid_data.shape,
        'n_atoms': sub_info['n_atoms'],
        'lattice_vectors': sub_info['lvec'].tolist(),
        'best_convention': shift_detection['best_convention'],
        'best_g0': [float(x) for x in shift_detection['best_g0']],
        'dg': list(shift_detection['dg']),
        'alignment_stats': alignment_stats,
        'sampling_results': sampling_results['stats'],
        'sampling_comparison': {k: {kk: float(vv) if isinstance(vv, (np.float32, np.float64)) else 
                                   (vv.tolist() if isinstance(vv, np.ndarray) else vv) 
                                   for kk, vv in v.items()} 
                                for k, v in comparison.items()} if 'comparison' in locals() else None,
        'component_plots': component_plots if 'component_plots' in locals() else None,
        'plot_results': plot_results,
        'summary_plot': summary_path
    }
    
    report_path = os.path.join(save_dir, f"{base_name}_alignment_report.json")
    with open(report_path, 'w') as f:
        json.dump(report, f, indent=2)
    
    if verbose:
        print(f"\n{'='*60}")
        print(f"Verification Complete")
        print(f"{'='*60}")
        print(f"Report saved: {report_path}")
        print(f"Summary plot: {summary_path}")
        print(f"Diagnostic plots ({len(plot_results)} conventions):")
        for pr in plot_results:
            if pr['status'] == 'ok':
                print(f"  - {pr['path']}")
    
    return report
