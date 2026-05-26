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
import time

# matplotlib with non-interactive backend for headless use
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Import from existing modules
from .InteractionEnergy import load_xyz_with_REQs
from .RigidBodyAFM import sample_gridff_single_atom
from .RigidBodyDynamics import RigidBodyDynamics, _reqs_to_plq

# CRITICAL: _reqs_to_plq expects REQ.y to be sqrt(EvdW), NOT raw EvdW.
# When reading from ElementTypes.dat, you MUST sqrt the E value before calling _reqs_to_plq.
# This matches the GridFF generation convention in ocl_GridFF_new.py::make_atoms_arrays(bSqrtEvdw=True).

from .MolecularDynamics import MolecularDynamics

# Import Ewald2D for electrostatics comparison
from pyBall.Ewald2D import Ewald2D


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


def load_bspline_gridff(grid_path):
    """Load Bspline GridFF and its JSON metadata."""
    meta = load_gridff_metadata(grid_path)
    if meta is None:
        raise FileNotFoundError(f"load_bspline_gridff(): missing metadata JSON for grid '{grid_path}'")
    grid = load_gridff_array(grid_path)
    ns = tuple(int(x) for x in meta['ns'])
    if tuple(grid.shape[:3]) != ns:
        raise ValueError(f"load_bspline_gridff(): grid.shape[:3]={grid.shape[:3]} != meta.ns={ns} for '{grid_path}'")
    return grid, meta


def init_gridff_sampler_md(grid_path, apos0, nSystems, use_texture=False):
    """Initialize MolecularDynamics for fast GridFF sampling on many rigid transforms.

    - Reuses OpenCL kernel `sampleGridFF_Bspline_points`.
    - Does NOT recompute GridFF; it only uploads existing Bspline grid once.
    """
    grid, meta = load_bspline_gridff(grid_path)
    g0 = np.array(meta['g0'], dtype=np.float32)
    dg = np.array(meta['dg'], dtype=np.float32)
    ns = tuple(int(x) for x in meta['ns'])
    print(f"init_gridff_sampler_md(): grid_path='{grid_path}' ns={ns} g0={g0.tolist()} dg={dg.tolist()} use_texture={use_texture}")
    apos0 = np.asarray(apos0, dtype=np.float32)
    if apos0.ndim != 2 or apos0.shape[1] < 3:
        raise ValueError(f"init_gridff_sampler_md(): apos0 must have shape (natoms,3+) got {apos0.shape}")
    natoms = int(apos0.shape[0])
    REQs0 = np.zeros((natoms, 4), dtype=np.float32)
    md = MolecularDynamics(nloc=32, debug_build_options='-DDBG_UFF=0')
    md.init_rigid_molecule_batch(apos0[:, :3].copy(), REQs0, nSystems=int(nSystems))
    md.initGridFF(grid_shape=ns, bspline_data=grid, grid_p0=g0, grid_step=dg, use_texture=bool(use_texture), r_damp=0.0, alpha_morse=0.0, bKernels=True)
    return md, meta


def sample_gridff_channels_rigid(md, transforms, PLQH_channels=((1.0, 0.0, 0.0, 0.0), (0.0, 1.0, 0.0, 0.0), (0.0, 0.0, 1.0, 0.0))):
    """Sample GridFF channels (Pauli/London/Coulomb) at atom positions for many rigid transforms.

    Returns:
        Es: (nconf,natoms,nch) float32 energies per atom per channel
    """
    T = np.asarray(transforms, dtype=np.float32).reshape(-1, 3, 4)
    nconf = int(T.shape[0])
    natoms = int(md.natoms)
    nch = int(len(PLQH_channels))
    Es = np.empty((nconf, natoms, nch), dtype=np.float32)
    print(f"sample_gridff_channels_rigid(): nconf={nconf} natoms={natoms} nch={nch}")
    chunk = int(md.nSystems)
    if chunk <= 0:
        raise ValueError(f"sample_gridff_channels_rigid(): invalid md.nSystems={md.nSystems}")
    for i0 in range(0, nconf, chunk):
        nchunk = min(chunk, nconf - i0)
        md.upload_rigid_transforms(T[i0:i0+nchunk], iSys0=0)
        for ich, plqh in enumerate(PLQH_channels):
            md.run_sampleGridFF_Bspline_points(nSystems=nchunk, PLQH=plqh)
            aforce = np.empty((nchunk, md.nvecs, 4), dtype=np.float32)
            md.fromGPU('aforce', aforce)
            Es[i0:i0+nchunk, :, ich] = -aforce[:, :natoms, 3]
    return Es


def fdbm_build_feature_matrix(Es_PLQ, type_ids, ntypes, use_london=True):
    """Build feature matrix for linear fitting: columns are per-type sums of sampled channels."""
    Es_PLQ = np.asarray(Es_PLQ)
    if Es_PLQ.ndim != 3 or Es_PLQ.shape[2] < 3:
        raise ValueError(f"fdbm_build_feature_matrix(): Es_PLQ must have shape (nconf,natoms,>=3) got {Es_PLQ.shape}")
    type_ids = np.asarray(type_ids, dtype=np.int32).reshape(-1)
    natoms = int(Es_PLQ.shape[1])
    if len(type_ids) != natoms:
        raise ValueError(f"fdbm_build_feature_matrix(): len(type_ids)={len(type_ids)} != natoms={natoms}")
    ntypes = int(ntypes)
    G = np.zeros((natoms, ntypes), dtype=np.float64)
    G[np.arange(natoms), type_ids] = 1.0
    Psum = Es_PLQ[:, :, 0].astype(np.float64) @ G
    if use_london:
        Lsum = Es_PLQ[:, :, 1].astype(np.float64) @ G
        return np.concatenate([Psum, Lsum], axis=1)
    return Psum


def fdbm_make_mock_reference(Es_PLQ, charges, type_ids, P_true, L_true=None, use_london=True, pauli_alpha=1.0, pauli_rescale=1.0, london_rescale=1.0, coulomb_rescale=1.0, noise_sigma=0.0, rng=None):
    """Generate mock reference energies from GridFF samples, with controlled perturbations."""
    Es = np.asarray(Es_PLQ, dtype=np.float64)
    if Es.ndim != 3 or Es.shape[2] < 3:
        raise ValueError(f"fdbm_make_mock_reference(): Es_PLQ must have shape (nconf,natoms,>=3) got {Es.shape}")
    nconf, natoms = int(Es.shape[0]), int(Es.shape[1])
    charges = np.asarray(charges, dtype=np.float64).reshape(-1)
    type_ids = np.asarray(type_ids, dtype=np.int32).reshape(-1)
    if len(charges) != natoms:
        raise ValueError(f"fdbm_make_mock_reference(): len(charges)={len(charges)} != natoms={natoms}")
    if len(type_ids) != natoms:
        raise ValueError(f"fdbm_make_mock_reference(): len(type_ids)={len(type_ids)} != natoms={natoms}")
    P_true = np.asarray(P_true, dtype=np.float64).reshape(-1)
    ntypes = int(P_true.shape[0])
    if use_london:
        if L_true is None:
            raise ValueError("fdbm_make_mock_reference(): use_london=True requires L_true")
        L_true = np.asarray(L_true, dtype=np.float64).reshape(-1)
        if len(L_true) != ntypes:
            raise ValueError(f"fdbm_make_mock_reference(): len(L_true)={len(L_true)} != ntypes={ntypes}")
    P_i = P_true[type_ids]
    L_i = np.zeros(natoms, dtype=np.float64) if (not use_london) else L_true[type_ids]
    Vp = pauli_rescale * Es[:, :, 0]
    Vl = london_rescale * Es[:, :, 1]
    Vq = coulomb_rescale * Es[:, :, 2]
    if pauli_alpha != 1.0:
        Vp_eff = np.power(np.clip(Vp, 0.0, None), float(pauli_alpha))
    else:
        Vp_eff = Vp
    E_pauli = (Vp_eff * P_i[None, :]).sum(axis=1)
    E_london = (Vl * L_i[None, :]).sum(axis=1) if use_london else np.zeros(nconf, dtype=np.float64)
    E_coul = (Vq * charges[None, :]).sum(axis=1)
    if rng is None:
        rng = np.random.default_rng(0)
    noise = rng.normal(0.0, float(noise_sigma), size=nconf) if (noise_sigma != 0.0) else np.zeros(nconf, dtype=np.float64)
    E_ref = E_pauli + E_london + E_coul + noise
    parts = {'pauli': E_pauli, 'london': E_london, 'coulomb': E_coul, 'noise': noise}
    print(f"fdbm_make_mock_reference(): E_ref range=[{E_ref.min():.6f},{E_ref.max():.6f}] eV noise_sigma={noise_sigma} pauli_alpha={pauli_alpha} use_london={use_london}")
    return E_ref, parts


def save_xyz_movie_with_energies(fname, enames, apos_list, energies, qs=None):
    """Save multi-frame XYZ with total energy in comment line.

    Args:
        fname: output .xyz path
        enames: list of element symbols (natoms,)
        apos_list: list of (natoms,3) arrays or (nframes,natoms,3) array
        energies: (nframes,) array of total energies [eV]
        qs: optional (natoms,) charges to include as 5th column
    """
    from ..atomicUtils import writeToXYZ
    apos_list = np.asarray(apos_list)
    if apos_list.ndim == 2:
        apos_list = apos_list[None, :]
    nframes = int(apos_list.shape[0])
    energies = np.asarray(energies).reshape(-1)
    if len(energies) != nframes:
        raise ValueError(f"save_xyz_movie_with_energies(): len(energies)={len(energies)} != nframes={nframes}")
    with open(fname, 'w') as f:
        for i in range(nframes):
            comment = f"E={energies[i]:.10f} eV"
            writeToXYZ(f, enames, apos_list[i], qs=qs, comment=comment, bHeader=True)


def load_xyz_movie_with_energies(fname):
    """Load multi-frame XYZ and parse energies from comment lines.

    Returns:
        enames: list of element symbols (natoms,)
        apos: (nframes,natoms,3) positions
        energies: (nframes,) energies parsed from comments
        qs: (nframes,natoms) charges if present, else None
    """
    from ..atomicUtils import load_xyz_movie
    trj = load_xyz_movie(fname)
    nframes = len(trj)
    enames = trj[0][0]
    natoms = len(enames)
    apos = np.zeros((nframes, natoms, 3), dtype=np.float64)
    qs = np.zeros((nframes, natoms), dtype=np.float64)
    energies = np.zeros(nframes, dtype=np.float64)
    has_charges = False
    for i, (es_i, apos_i, qs_i, rs_i, comment_i) in enumerate(trj):
        if es_i != enames:
            raise ValueError(f"load_xyz_movie_with_energies(): element symbols differ in frame {i}")
        apos[i] = apos_i
        qs[i] = qs_i
        if np.any(qs_i != 0.0):
            has_charges = True
        if comment_i:
            comment_i = comment_i.strip()
            for token in comment_i.split():
                if token.startswith('E='):
                    try:
                        energies[i] = float(token[2:].split()[0])
                    except (ValueError, IndexError):
                        pass
                    break
    return enames, apos, energies, qs if has_charges else None


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
        'from_json': (-Lx/2, -Ly/2, z_top),  # GridFF generation convention: centered XY, z at top atom
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

def plot_atoms_overlay(ax, atoms_xyz, atoms_enames, z_atom_range=None, g0=None, dg=None, 
                      marker='.', s=10, alpha=0.7, coords='xy', extra_mask=None):
    """
    Plot atoms as overlay with consistent styling and optional z-filtering.
    
    Args:
        ax: Matplotlib axis
        atoms_xyz: (n, 3) atom positions
        atoms_enames: Element names for color coding
        z_atom_range: Z-range below top atom to show (None = show all)
        g0: Grid origin (x0, y0, z0) for z-filtering
        dg: Grid spacing (dx, dy, dz) for z-filtering
        marker: Marker style (default '.')
        s: Marker size (default 10)
        alpha: Transparency (default 0.7)
        coords: Which coordinates to plot ('xy' or 'xz')
        extra_mask: Additional boolean mask to apply (e.g., y-filtering for XZ)
    """
    # Color mapping
    colors = ['purple' if e in ['Ca', 'Na'] else 'green' if e in ['F', 'Cl'] else 'gray' 
              for e in atoms_enames]
    
    # Apply z-filtering if requested
    mask = np.ones(len(atoms_xyz), dtype=bool)
    if z_atom_range is not None and g0 is not None and dg is not None:
        z_top_atom = np.max(atoms_xyz[:, 2])
        mask &= atoms_xyz[:, 2] >= (z_top_atom - z_atom_range)
    
    # Apply extra mask if provided
    if extra_mask is not None:
        mask &= extra_mask
    
    atoms_xyz = atoms_xyz[mask]
    colors = np.array(colors)[mask]
    
    if len(atoms_xyz) > 0:
        if coords == 'xy':
            ax.scatter(atoms_xyz[:, 0], atoms_xyz[:, 1], c=colors, s=s, alpha=alpha, marker=marker)
        elif coords == 'xz':
            ax.scatter(atoms_xyz[:, 0], atoms_xyz[:, 2], c=colors, s=s, alpha=alpha, marker=marker)




def plot_gridff_diagnostics(grid_data, sub_apos, sub_enames, lvec, iz_slices=None, iy_slice=None, save_path='grid_diagnostics.png', g0=None, dg=None, z_marks=None, channel_name=None, z_atom_range=5.0, mol_apos=None, mol_enames=None):
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
        z_atom_range: Z-range below top atom to show in XY plots (default 5.0 A)
        mol_apos: (nmol, 3) or (nframes, nmol, 3) molecule sample positions to overlay (optional)
        mol_enames: List of element names for molecule atoms (optional)
    """
    t0 = time.perf_counter()
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

    # Handle molecule samples
    mol_apos_flat = None
    mol_colors = None
    if mol_apos is not None:
        mol_apos = np.asarray(mol_apos)
        if mol_apos.ndim == 3:
            mol_apos_flat = mol_apos.reshape(-1, mol_apos.shape[1])
        else:
            mol_apos_flat = mol_apos
        if mol_enames is not None:
            # Repeat colors for each atom in each sample
            natoms_per_sample = len(mol_enames)
            mol_colors = ['red' if e == 'O' else 'orange' if e == 'H' else 'blue' for e in mol_enames]
            mol_colors = mol_colors * (len(mol_apos_flat) // natoms_per_sample)
        else:
            mol_colors = ['red'] * len(mol_apos_flat)

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

    for i in range(nch):
        # XY Slices
        for j, iz in enumerate(iz_slices):
            z_val = g0[2] + iz * dz
            extent = [g0[0], g0[0] + ax, g0[1], g0[1] + ay]
            im = axs[i, j].imshow(grid_data[:, :, iz, i].T, extent=extent,
                                   origin='lower', cmap='bwr', aspect='equal')
            # Overlay substrate atoms with z-filtering
            plot_atoms_overlay(axs[i, j], sub_apos, sub_enames, z_atom_range=z_atom_range, g0=g0, dg=dg, coords='xy')
            # Overlay molecule samples (no z-filter, show all samples)
            if mol_apos_flat is not None:
                axs[i, j].scatter(mol_apos_flat[:, 0], mol_apos_flat[:, 1], c=mol_colors, s=15, alpha=0.6, marker='o', edgecolors='black', linewidth=0.5, label='mol samples')
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
        # Overlay substrate atoms with y-filtering and z-filtering
        mask_y = np.abs(sub_apos[:, 1] - y_val) < 2.0
        plot_atoms_overlay(axs[i, -1], sub_apos, sub_enames, z_atom_range=z_atom_range, g0=g0, dg=dg, coords='xz', extra_mask=mask_y)
        # Overlay molecule samples (y-filtered to slice)
        if mol_apos_flat is not None:
            mask_mol_y = np.abs(mol_apos_flat[:, 1] - y_val) < 2.0
            if np.any(mask_mol_y):
                axs[i, -1].scatter(mol_apos_flat[mask_mol_y, 0], mol_apos_flat[mask_mol_y, 2],
                                   c=np.array(mol_colors)[mask_mol_y], s=15, alpha=0.6, marker='o', edgecolors='black', linewidth=0.5, label='mol samples')
        if z_marks is not None:
            for zm in z_marks:
                axs[i, -1].axhline(zm, color='k', linestyle='--', linewidth=1.0, alpha=0.7, label=f'z={zm:.2f}')
        axs[i, -1].set_title(f"{names[i]} XZ at y={y_val:.2f} A (iy={iy_slice})")
        axs[i, -1].set_xlabel('x [A]')
        axs[i, -1].set_ylabel('z [A]')
        axs[i, -1].set_ylim(g0[2], g0[2] + min(az, 20))  # Limit z range for visibility
        plt.colorbar(im_xz, ax=axs[i, -1])

    plt.tight_layout()
    t1 = time.perf_counter()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    t2 = time.perf_counter()
    plt.close(fig)
    print(f"Saved GridFF diagnostics to {save_path} (render: {t1-t0:.3f}s, save: {t2-t1:.3f}s, total: {t2-t0:.3f}s)")
    return save_path


def plot_alignment_summary(grid_data, g0, dg, atoms_xyz, atoms_enames, save_path, iz_top=None, iy_center=None, z_atom_range=2.0, mol_apos=None, mol_enames=None, plq_coeffs=None, plq_coeffs2=None, zmin_offset=2.0, z_ylim=None, vmin_vmax_xz=None, plot_diagnostics=True):
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
        z_atom_range: Z-range below top atom to show (default 2.0)
        mol_apos: (nmol, 3) or (nframes, nmol, 3) molecule sample positions to overlay (optional)
        mol_enames: List of element names for molecule atoms (optional)
        plq_coeffs: Tuple (P, L, Q) coefficients for total potential (default: (1.0, 1.0, 1.0))
        plq_coeffs2: Optional second set of (P, L, Q) coefficients (e.g. for second atom type) plotted alongside first in Z-profiles
        zmin_offset: Offset above g0[2] for XZ color scale normalization (default 2.0 A)
        z_ylim: Y-axis limits for Z-profile (default: +/-0.5 eV)
        vmin_vmax_xz: Optional tuple (vmin, vmax) for XZ slice color scale (default: auto from data above zmin_offset)
        plot_diagnostics: If False, skip plotting and return None (default: True)
    """
    if not plot_diagnostics:
        return None

    nx, ny, nz, nch = grid_data.shape
    dx, dy, dz = dg

    # Default PLQ coefficients (P, L, Q)
    if plq_coeffs is None:
        plq_coeffs = (1.0, 1.0, 1.0)
    P, L, Q = plq_coeffs

    # Compute total potential: E = P*Pauli + L*London + Q*Coulomb
    total_potential = P * grid_data[..., 0:1] + L * grid_data[..., 1:2] + Q * grid_data[..., 2:3]

    # Auto-detect top layer if not specified
    if iz_top is None:
        z_top = np.max(atoms_xyz[:, 2])
        iz_top = int((z_top - g0[2]) / dz)
        iz_top = max(0, min(nz-1, iz_top))

    if iy_center is None:
        iy_center = ny // 2

    # Handle molecule samples
    mol_apos_flat = None
    mol_colors = None
    if mol_apos is not None:
        mol_apos = np.asarray(mol_apos)
        if mol_apos.ndim == 3:
            mol_apos_flat = mol_apos.reshape(-1, mol_apos.shape[1])
        else:
            mol_apos_flat = mol_apos
        if mol_enames is not None:
            natoms_per_sample = len(mol_enames)
            mol_colors = ['red' if e == 'O' else 'orange' if e == 'H' else 'blue' for e in mol_enames]
            mol_colors = mol_colors * (len(mol_apos_flat) // natoms_per_sample)
        else:
            mol_colors = ['red'] * len(mol_apos_flat)

    # Color scale: XY uses symmetric range from that slice only, XZ uses data above zmin
    # XY slice - symmetric range from this slice only
    xy_slice_data = total_potential[:, :, iz_top, 0]
    vmax_xy = max(abs(xy_slice_data.min()), abs(xy_slice_data.max()))
    vmin_xy, vmax_xy = -vmax_xy, vmax_xy if vmax_xy > 0 else (None, None)

    # XZ slice - data above zmin (compute from the XZ slice only, not whole grid)
    # zmin_offset is absolute z (e.g., 2.0 means 2A above surface at z=0)
    zmin = zmin_offset
    iz_min = int((zmin - g0[2]) / dz)
    iz_min = max(0, min(nz-1, iz_min))
    # For display: start from z=0 (surface)
    iz_surface = int((0.0 - g0[2]) / dz)
    iz_surface = max(0, min(nz-1, iz_surface))
    xz_slice_data = total_potential[:, iy_center, :, 0].T  # Get XZ slice first
    xz_data_above = xz_slice_data[iz_min:, :]  # Only data above zmin for vmin/vmax
    vmax_xz = max(abs(xz_data_above.min()), abs(xz_data_above.max()))
    vmin_xz, vmax_xz = -vmax_xz, vmax_xz if vmax_xz > 0 else (None, None)

    # Default Z-profile y-limits
    if z_ylim is None:
        z_ylim = (-0.5, 0.5)

    # Precompute atom colors once (no on-the-fly filtering)
    atoms_colors = ['purple' if e in ['Ca', 'Na'] else 'green' if e in ['F', 'Cl'] else 'gray'
                    for e in atoms_enames]

    fig, axes = plt.subplots(2, 3, figsize=(18, 12))

    # Helper function for XZ slice plotting
    def plot_xz_slice(ax, iy, title_suffix):
        """Plot XZ slice at given iy index with consistent styling."""
        y_val = g0[1] + iy * dy
        # Display from z=0 (surface) up, but vmin/vmax from z > zmin_offset
        xz_slice_full = total_potential[:, iy, :, 0].T
        xz_slice_data = xz_slice_full[iz_surface:, :]  # Data from z=0 (surface)
        zmin_display = g0[2] + iz_surface * dz  # z=0 (surface)
        extent_xz = [g0[0], g0[0] + nx*dx, zmin_display, g0[2] + nz*dz]
        # Use manual vmin/vmax if provided, otherwise auto from z > zmin_offset
        if vmin_vmax_xz is not None:
            vmin_xz_use, vmax_xz_use = vmin_vmax_xz
        else:
            vmin_xz_use, vmax_xz_use = vmin_xz, vmax_xz
        im = ax.imshow(xz_slice_data, extent=extent_xz,
                       origin='lower', cmap='bwr', aspect='equal', vmin=vmin_xz_use, vmax=vmax_xz_use)
        # Plot substrate atoms (same array for all panels)
        ax.scatter(atoms_xyz[:, 0], atoms_xyz[:, 2],
                   c=atoms_colors, s=10, alpha=0.7, marker='.')
        # Plot molecule samples (same array for all panels)
        if mol_apos_flat is not None:
            ax.scatter(mol_apos_flat[:, 0], mol_apos_flat[:, 2],
                       c=mol_colors, s=15, alpha=0.6, marker='.', edgecolors='black', linewidth=0.5, label='mol samples')
        ax.set_title(f'Total Potential (P={P:.3f},L={L:.3f},Q={Q:.3f}) XZ at y={y_val:.3f}A (iy={iy}) {title_suffix}')
        ax.set_xlabel('x [A]')
        ax.set_ylabel('z [A]')
        ax.set_ylim(zmin_display, g0[2] + min(nz*dz, 20))  # Show from z=0 (surface)
        plt.colorbar(im, ax=ax)

    # 1. XY slice at top layer - Total potential
    ax = axes[0, 0]
    z_val = g0[2] + iz_top * dz
    extent = [g0[0], g0[0] + nx*dx, g0[1], g0[1] + ny*dy]
    im = ax.imshow(total_potential[:, :, iz_top, 0].T, extent=extent,
                   origin='lower', cmap='bwr', aspect='equal', vmin=vmin_xy, vmax=vmax_xy)
    plot_atoms_overlay(ax, atoms_xyz, atoms_enames, z_atom_range=z_atom_range, g0=g0, dg=dg)
    # Overlay molecule samples
    if mol_apos_flat is not None:
        ax.scatter(mol_apos_flat[:, 0], mol_apos_flat[:, 1], c=mol_colors, s=15, alpha=0.6, marker='.', edgecolors='black', linewidth=0.5, label='mol samples')
    ax.set_title(f'Total Potential (P={P:.3f},L={L:.3f},Q={Q:.3f}) XY at z={z_val:.3f}A (iz={iz_top})')
    ax.set_xlabel('x [A]')
    ax.set_ylabel('y [A]')
    plt.colorbar(im, ax=ax)

    # 2. XZ slice at center Y
    plot_xz_slice(axes[0, 1], iy_center, '')

    # 3. XZ slice at 1/4 of cell size (iy=32, ~2.82 Å)
    plot_xz_slice(axes[0, 2], ny//4, '(1/4 cell)')
    
    # 4. 1D profiles through atom centers
    # Helper function for 1D profile plotting (X or Z)
    def plot_1d_profile(ax, axis, fixed_idx, title_suffix):
        """Plot 1D profile along specified axis with consistent styling and atom markers.
        
        Args:
            axis: 'x' for X-profile (vary x, fixed y,z), 'z' for Z-profile (vary z, fixed x,y)
            fixed_idx: tuple of (fixed1, fixed2) indices for the other dimensions
            title_suffix: suffix for title
        """
        if axis == 'x':
            cy, cz = fixed_idx
            x_coords = g0[0] + np.arange(nx) * dx
            ax.plot(x_coords, total_potential[:, cy, cz, 0], 'k-', linewidth=2, label=f'Total (P={P:.3f},L={L:.3f},Q={Q:.3f})')
            ax.plot(x_coords, grid_data[:, cy, cz, 0], 'b--', alpha=0.5, label='Pauli')
            ax.plot(x_coords, grid_data[:, cy, cz, 1], 'r--', alpha=0.5, label='London')
            ax.plot(x_coords, grid_data[:, cy, cz, 2], 'g--', alpha=0.5, label='Coulomb')
            ax.axhline(0, color='k', linestyle=':', alpha=0.3)
            # Mark all atom positions along x (use precomputed colors)
            for x_atom, color in zip(atoms_xyz[:, 0], atoms_colors):
                ax.axvline(x_atom, color=color, linestyle=':', alpha=0.5, linewidth=1.5)
            ax.set_xlabel('x [A]')
            ax.set_ylabel('Energy [eV]')
            ax.set_title(f'X-Profile at y={g0[1]+cy*dy:.3f}A, z={g0[2]+cz*dz:.3f}A {title_suffix}')
        elif axis == 'z':
            cx, cy = fixed_idx
            z_coords = g0[2] + np.arange(nz) * dz
            ax.plot(z_coords, total_potential[cx, cy, :, 0], 'k-', linewidth=2, label=f'H-Total (P={P:.3f},L={L:.3f},Q={Q:.3f})')
            if plq_coeffs2 is not None:
                P2, L2, Q2 = plq_coeffs2[0], plq_coeffs2[1], plq_coeffs2[2]
                total2 = P2*grid_data[...,0] + L2*grid_data[...,1] + Q2*grid_data[...,2]
                ax.plot(z_coords, total2[cx, cy, :], 'm-', linewidth=2, label=f'O-Total (P={P2:.3f},L={L2:.3f},Q={Q2:.3f})')
            ax.plot(z_coords, grid_data[cx, cy, :, 0], 'b--', alpha=0.5, label='Pauli')
            ax.plot(z_coords, grid_data[cx, cy, :, 1], 'r--', alpha=0.5, label='London')
            ax.plot(z_coords, grid_data[cx, cy, :, 2], 'g--', alpha=0.5, label='Coulomb')
            ax.axhline(0, color='k', linestyle=':', alpha=0.3)
            ax.axvline(z_val, color='purple', linestyle='--', alpha=0.5, label=f'z_top={z_val:.3f}A')
            # Mark all atom positions along z (use precomputed colors)
            for z_atom, color in zip(atoms_xyz[:, 2], atoms_colors):
                ax.axvline(z_atom, color=color, linestyle=':', alpha=0.5, linewidth=1.5)
            ax.set_xlabel('z [A]')
            ax.set_ylabel('Energy [eV]')
            ax.set_title(f'Z-Profile at x={g0[0]+cx*dx:.3f}A, y={g0[1]+cy*dy:.3f}A {title_suffix}')
            ax.set_xlim(g0[2], g0[2] + min(nz*dz, 20))  # Consistent z-range limit
            ax.set_ylim(z_ylim)
        
        ax.legend()
        ax.grid(True, alpha=0.3)

    # 4. X-profile at center Y, top Z
    cx = nx // 2
    cy = ny // 2
    plot_1d_profile(axes[1, 0], 'x', (cy, iz_top), '')

    # 5. Z-profile above Cl (at ix=0, iy=0)
    plot_1d_profile(axes[1, 1], 'z', (0, 0), 'above Cl-')

    # 6. Z-profile above Na (at ix=32, iy=0 = 2.82 A along x)
    ix_Na = int(round(2.821 / dx))
    plot_1d_profile(axes[1, 2], 'z', (ix_Na, 0), 'above Na+')
    
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


def compare_sampling_methods(grid_data, g0, dg, rbd, grid, positions, 
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
        rbd: Pre-initialized RigidBodyDynamics instance
        grid: Grid data array
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
        
        # OpenCL: sample with B-spline interpolation using reusable RigidBodyDynamics
        positions_np = np.asarray(positions, dtype=np.float32)
        n_bodies = len(positions_np)
        
        rbd.realloc(n_bodies=n_bodies, num_atoms=1)
        rbd.enames = ['TestAtom']
        rbd.atom_types_assigned = ['TestAtom']
        
        reqs = np.array([comp['atom_req']], dtype=np.float32)
        rbd.atom_REQ = reqs
        rbd.atom_masses = np.array([1.008], dtype=np.float32)
        rbd.mass_physical = 1.008
        rbd.mass_trans = 1.008
        rbd.mass_rot = 1.008
        
        atom_plq_single = _reqs_to_plq(reqs, alpha=alpha_morse)
        atom_plq = np.repeat(atom_plq_single[None, :, :], n_bodies, axis=0).reshape(n_bodies, 4)
        rbd.atom_PLQ = atom_plq.copy()
        
        pos4 = np.zeros((n_bodies, 4), dtype=np.float32)
        pos4[:, :3] = positions_np
        pos4[:, 3] = 1.008
        
        quat4 = np.zeros((n_bodies, 4), dtype=np.float32)
        quat4[:, 3] = 1.0
        zero4 = np.zeros((n_bodies, 4), dtype=np.float32)
        
        Iinv_relax = np.eye(3, dtype=np.float32)
        atom_body = np.zeros((n_bodies, 1, 3), dtype=np.float32)
        
        rbd.upload_state(pos4, quat4, zero4, zero4, rbd.mass_trans, 1.0 / rbd.mass_trans, np.repeat(Iinv_relax[None, :, :], n_bodies, axis=0), atom_body, atom_PLQ=atom_plq)
        
        # Re-init grid after realloc
        rbd.init_gridff(grid, grid_p0=g0, grid_step=dg)
        
        # Run 1 step with dt=0.0 to evaluate forces without moving
        rbd.run_gridff(num_steps=1, dt=0.0)
        
        outputs = rbd.download_selected(('atom_force',))
        opencl_vals = outputs['atom_force'][:, 0, 3]
        
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


def sample_grid_at_atoms_opencl(rbd, grid, atoms_xyz, atom_req, alpha_morse, 
                                  grid_p0=(0.0, 0.0, 0.0), grid_step=None, verbose=True):
    """
    Sample GridFF at atom positions using OpenCL B-spline interpolation.
    
    Args:
        rbd: Pre-initialized RigidBodyDynamics instance
        grid: Grid data array
        atoms_xyz: (n, 3) atom positions
        atom_req: Tuple (R, E, Q, H) for test atom
        alpha_morse: Alpha parameter
        grid_p0: Grid origin (x0, y0, z0)
        grid_step: Grid spacing (dx, dy, dz)
        verbose: Print statistics
        
    Returns:
        dict with 'forces', 'energies', 'stats'
    """
    positions_np = np.asarray(atoms_xyz, dtype=np.float32)
    n_bodies = len(positions_np)
    
    rbd.realloc(n_bodies=n_bodies, num_atoms=1)
    rbd.enames = ['TestAtom']
    rbd.atom_types_assigned = ['TestAtom']
    
    reqs = np.array([atom_req], dtype=np.float32)
    rbd.atom_REQ = reqs
    rbd.atom_masses = np.array([1.008], dtype=np.float32)
    rbd.mass_physical = 1.008
    rbd.mass_trans = 1.008
    rbd.mass_rot = 1.008
    
    atom_plq_single = _reqs_to_plq(reqs, alpha=alpha_morse)
    atom_plq = np.repeat(atom_plq_single[None, :, :], n_bodies, axis=0).reshape(n_bodies, 4)
    rbd.atom_PLQ = atom_plq.copy()
    
    pos4 = np.zeros((n_bodies, 4), dtype=np.float32)
    pos4[:, :3] = positions_np
    pos4[:, 3] = 1.008
    
    quat4 = np.zeros((n_bodies, 4), dtype=np.float32)
    quat4[:, 3] = 1.0
    zero4 = np.zeros((n_bodies, 4), dtype=np.float32)
    
    Iinv_relax = np.eye(3, dtype=np.float32)
    atom_body = np.zeros((n_bodies, 1, 3), dtype=np.float32)
    
    rbd.upload_state(pos4, quat4, zero4, zero4, rbd.mass_trans, 1.0 / rbd.mass_trans, 
                    np.repeat(Iinv_relax[None, :, :], n_bodies, axis=0), atom_body, atom_PLQ=atom_plq)
    
    # Re-init grid after realloc
    rbd.init_gridff(grid, grid_p0=grid_p0, grid_step=grid_step)
    
    # Run 1 step with dt=0.0 to evaluate forces without moving
    rbd.run_gridff(num_steps=1, dt=0.0)
    
    outputs = rbd.download_selected(('atom_force',))
    forces = outputs['atom_force'][:, 0, :3]
    energies = outputs['atom_force'][:, 0, 3]
    
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
    
    return {'forces': forces, 'energies': energies, 'stats': stats}


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
                               alpha_morse=1.5, z_atom_range=2.0, verbose=True):
    """
    Run complete alignment verification workflow.
    
    Args:
        grid_path: Path to GridFF .npy file
        substrate_path: Path to substrate .xyz file
        save_dir: Directory for output plots and reports
        test_conventions: List of convention names to test (None = test all)
        atom_req: Tuple (R, E, Q, H) for test atom
        alpha_morse: Alpha parameter for REQ to PLQ conversion
        z_atom_range: Z-range below top atom to show in XY plots (default 2.0 A)
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
    
    # If specific conventions are requested, use the first one directly
    if test_conventions is not None and len(test_conventions) > 0:
        conv_name = test_conventions[0]
        if conv_name in metadata['conventions']:
            best_g0 = metadata['conventions'][conv_name]
            best_dg = metadata['dg']
            shift_detection = {
                'best_convention': conv_name,
                'best_g0': best_g0,
                'dg': best_dg,
                'all_results': {}
            }
            if verbose:
                print(f"Using specified convention: {conv_name}")
                print(f"  g0 = {best_g0}")
        else:
            if verbose:
                print(f"Warning: Convention '{conv_name}' not found, running auto-detection...")
            shift_detection = auto_detect_shift(
                grid_data, sub_info['apos'], sub_info['lvec'], verbose=verbose
            )
            best_g0 = shift_detection['best_g0']
            best_dg = shift_detection['dg']
    # If metadata from JSON is available, use it directly
    elif metadata_json is not None:
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
        summary_path,
        z_atom_range=z_atom_range
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
    
    # Compare with simple coefficient summation for verification
    if verbose:
        print(f"\n[5/6] Comparing sampling methods (OpenCL vs simple)...")
        
        # Initialize RigidBodyDynamics once for all sampling
        rbd = RigidBodyDynamics(debug=False)
        rbd.realloc(n_bodies=1, num_atoms=1)
        rbd.init_gridff(grid_data, grid_p0=best_g0, grid_step=best_dg)
        
        comparison = compare_sampling_methods(
            grid_data, best_g0, best_dg, rbd, grid_data,
            test_positions,  # Use positions above substrate
            atom_req=atom_req,
            alpha_morse=alpha_morse,
            verbose=True
        )
        
        # Use the same rbd for the earlier sampling call
        sampling_results = sample_grid_at_atoms_opencl(
            rbd, grid_data, test_positions,
            atom_req=atom_req,
            alpha_morse=alpha_morse,
            grid_p0=best_g0,
            grid_step=best_dg,
            verbose=verbose
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
            # Vectorized: sum weighted channels in single operation
            simple_grid = np.sum(grid_data[:, :, :, comp['channels']] * np.array(plq_weights), axis=3, keepdims=True)
            
            simple_plot_path = os.path.join(save_dir, f"{base_name}_diagnostic_{comp['name']}_simple.png")
            plot_gridff_diagnostics(
                simple_grid, sub_info['apos'], sub_info['enames'], sub_info['lvec'],
                iz_slices=iz_slices,
                iy_slice=iy_slice,
                save_path=simple_plot_path,
                g0=best_g0, dg=best_dg, z_marks=z_marks, channel_name=comp['name'],
                z_atom_range=z_atom_range
            )
            
            # OpenCL plot: sample with B-spline using reusable rbd
            # Fully vectorized - NO Python loops over grid
            # Create 3D grid directly with meshgrid, single allocation
            z_vals = np.array([best_g0[2] + iz * best_dg[2] for iz in iz_slices], dtype=np.float32)
            xx, yy, zz = np.meshgrid(x_sample, y_sample, z_vals, indexing='ij')
            # Stack to (nx, ny, nz, 3), transpose to (nz, nx, ny, 3), reshape to (n_points, 3)
            all_positions = np.stack([xx, yy, zz], axis=-1).transpose(2, 0, 1, 3).reshape(-1, 3).astype(np.float32)
            n_bodies = len(all_positions)
            
            # Single realloc/init for all XY slices
            rbd.realloc(n_bodies=n_bodies, num_atoms=1)
            rbd.enames = ['TestAtom']
            rbd.atom_types_assigned = ['TestAtom']
            
            reqs = np.array([comp['atom_req']], dtype=np.float32)
            rbd.atom_REQ = reqs
            rbd.atom_masses = np.array([1.008], dtype=np.float32)
            rbd.mass_physical = 1.008
            rbd.mass_trans = 1.008
            rbd.mass_rot = 1.008
            
            atom_plq_single = _reqs_to_plq(reqs, alpha=alpha_morse)
            atom_plq = np.repeat(atom_plq_single[None, :, :], n_bodies, axis=0).reshape(n_bodies, 4)
            rbd.atom_PLQ = atom_plq.copy()
            
            pos4 = np.zeros((n_bodies, 4), dtype=np.float32)
            pos4[:, :3] = all_positions
            pos4[:, 3] = 1.008
            
            quat4 = np.zeros((n_bodies, 4), dtype=np.float32)
            quat4[:, 3] = 1.0
            zero4 = np.zeros((n_bodies, 4), dtype=np.float32)
            
            Iinv_relax = np.eye(3, dtype=np.float32)
            atom_body = np.zeros((n_bodies, 1, 3), dtype=np.float32)
            
            rbd.upload_state(pos4, quat4, zero4, zero4, rbd.mass_trans, 1.0 / rbd.mass_trans, 
                            np.repeat(Iinv_relax[None, :, :], n_bodies, axis=0), atom_body, atom_PLQ=atom_plq)
            
            # Re-init grid after realloc
            rbd.init_gridff(grid_data, grid_p0=best_g0, grid_step=best_dg)
            
            # Run 1 step with dt=0.0 to evaluate forces without moving
            rbd.run_gridff(num_steps=1, dt=0.0)
            
            outputs = rbd.download_selected(('atom_force',))
            all_energies = outputs['atom_force'][:, 0, 3]
            
            # Split results back into slices (vectorized)
            all_energies = all_energies.reshape(len(iz_slices), nx, ny)
            opencl_xy_slices = {}
            for i, iz in enumerate(iz_slices):
                opencl_xy_slices[iz] = all_energies[i]
            
            y_val = best_g0[1] + iy_slice * best_dg[1]
            # Sample only valid B-spline region (iz from 2 to nz-2) to avoid zeros at boundaries
            iz_min = 2
            iz_max = nz - 2
            nz_valid = iz_max - iz_min
            
            # Create 2D grid directly with meshgrid, single allocation
            xx, zz = np.meshgrid(x_sample, z_sample[iz_min:iz_max], indexing='ij')
            yy = np.full_like(xx, y_val)
            # Stack to (nx, nz, 3) then reshape to (n_points, 3)
            positions_xz = np.stack([xx, yy, zz], axis=-1).reshape(-1, 3).astype(np.float32)
            
            n_bodies = len(positions_xz)
            rbd.realloc(n_bodies=n_bodies, num_atoms=1)
            rbd.enames = ['TestAtom']
            rbd.atom_types_assigned = ['TestAtom']
            
            reqs = np.array([comp['atom_req']], dtype=np.float32)
            rbd.atom_REQ = reqs
            rbd.atom_masses = np.array([1.008], dtype=np.float32)
            rbd.mass_physical = 1.008
            rbd.mass_trans = 1.008
            rbd.mass_rot = 1.008
            
            atom_plq_single = _reqs_to_plq(reqs, alpha=alpha_morse)
            atom_plq = np.repeat(atom_plq_single[None, :, :], n_bodies, axis=0).reshape(n_bodies, 4)
            rbd.atom_PLQ = atom_plq.copy()
            
            pos4 = np.zeros((n_bodies, 4), dtype=np.float32)
            pos4[:, :3] = positions_xz
            pos4[:, 3] = 1.008
            
            quat4 = np.zeros((n_bodies, 4), dtype=np.float32)
            quat4[:, 3] = 1.0
            zero4 = np.zeros((n_bodies, 4), dtype=np.float32)
            
            Iinv_relax = np.eye(3, dtype=np.float32)
            atom_body = np.zeros((n_bodies, 1, 3), dtype=np.float32)
            
            rbd.upload_state(pos4, quat4, zero4, zero4, rbd.mass_trans, 1.0 / rbd.mass_trans, 
                            np.repeat(Iinv_relax[None, :, :], n_bodies, axis=0), atom_body, atom_PLQ=atom_plq)
            
            # Re-init grid after realloc
            rbd.init_gridff(grid_data, grid_p0=best_g0, grid_step=best_dg)
            
            # Run 1 step with dt=0.0 to evaluate forces without moving
            rbd.run_gridff(num_steps=1, dt=0.0)
            
            outputs = rbd.download_selected(('atom_force',))
            opencl_xz = outputs['atom_force'][:, 0, 3]
            # Reshape to (nx, nz_valid)
            opencl_xz = opencl_xz.reshape(nx, nz_valid)
            
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
                g0=best_g0, dg=best_dg, z_marks=z_marks, channel_name=comp['name'],
                z_atom_range=z_atom_range
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


# =============================================================================
# Section 6: Electrostatics Comparison (GridFF vs Ewald2D vs Brute Force)
# =============================================================================

def compare_electrostatics_methods(sys_at, gridff_path, sub_xyz_path, ew, save_dir='.', prefix='electrostatics',
                                  N_rep=20, coarse_spacing=1.0, grid_res=200,
                                  z_scan_heights=None, z_min_for_error=2.0, verbose=True):
    """
    Compare electrostatic potential from three methods:
    1. GridFF (OpenCL B-spline sampling, channel 2 - electrostatic)
    2. Ewald2D (2D Fourier representation)
    3. Brute force (direct Coulomb sum over periodic replicas)

    Evaluates all three methods at the same 1D coordinate arrays for fair comparison.
    Generates 2D slices for GridFF vs Ewald2D comparison.

    Parameters
    ----------
    sys_at : AtomicSystem
        Substrate with lattice vectors and charges
    gridff_path : str
        Path to GridFF Bspline_PLQd.npy file
    sub_xyz_path : str
        Path to substrate .xyz file (for grid loading)
    ew : Ewald2D
        Initialized Ewald2D instance
    save_dir : str
        Output directory for plots
    prefix : str
        Prefix for output filenames
    N_rep : int
        Number of periodic replicas for brute force (default: 20)
    coarse_spacing : float
        Point spacing (Å) for 1D brute force scans (default: 1.0)
    grid_res : int
        Resolution for 2D slices (default: 200)
    z_scan_heights : list
        Heights for XY slices, default [z_max+0.5, z_max+1.0, z_max+2.0]
    verbose : bool
        Print progress

    Returns
    -------
    dict with error statistics and paths to output figures
    """
    os.makedirs(save_dir, exist_ok=True)

    grid_data = np.load(gridff_path)
    nx, ny, nz = grid_data.shape[:3]
    meta = load_gridff_metadata(gridff_path)
    if meta is not None and ('g0' in meta) and ('dg' in meta):
        g0 = tuple(meta['g0'])
        dg = tuple(meta['dg'])
        if verbose:
            print(f"Loaded GridFF metadata g0={g0} dg={dg}")
    else:
        lvec = sys_at.lvec
        assert lvec is not None, "compare_electrostatics_methods(): sys_at.lvec is None"
        Lx_ = float(np.linalg.norm(lvec[0]))
        Ly_ = float(np.linalg.norm(lvec[1]))
        Lz_ = float(np.linalg.norm(lvec[2]))
        dg = (Lx_/nx, Ly_/ny, Lz_/nz)
        g0 = (-Lx_/2, -Ly_/2, np.max(sys_at.apos[:, 2]))
        if verbose:
            print(f"WARNING: GridFF metadata missing; inferred g0={g0} dg={dg}")

    apos = sys_at.apos
    assert apos is not None, "compare_electrostatics_methods(): sys_at.apos is None"
    lvec = sys_at.lvec
    assert lvec is not None, "compare_electrostatics_methods(): sys_at.lvec is None"
    Lx_cell = float(np.linalg.norm(lvec[0]))
    Ly_cell = float(np.linalg.norm(lvec[1]))
    z_top = float(np.max(apos[:, 2]))
    x_min = float(np.min(apos[:, 0]))
    y_min = float(np.min(apos[:, 1]))
    x_max = float(np.max(apos[:, 0]))
    y_max = float(np.max(apos[:, 1]))
    if verbose:
        print(f"AtomicSystem apos x=[{x_min:.3f},{x_max:.3f}] y=[{y_min:.3f},{y_max:.3f}] z_top={z_top:.3f}")
        print(f"Cell lengths from lvec |a|={Lx_cell:.3f} |b|={Ly_cell:.3f}")

    Lx = ew.Lx
    Ly = ew.Ly
    z_min = ew.z_min
    z_max = ew.z_max

    # We only compare electrostatics in vacuum above surface
    z_lo = z_top + 0.2
    z_hi = z_top + 5.0
    if z_scan_heights is None:
        z_scan_heights = [z_top + 3.0, z_top + 4.0, z_top + 5.0]

    # Define 1D coordinate arrays
    # Brute force: coarse spacing (expensive)
    z_scan_brute = np.arange(z_lo, z_hi + 1e-9, coarse_spacing)
    x_scan_brute = np.arange(x_min, x_min + Lx_cell + 1e-9, coarse_spacing)
    
    # GridFF/Ewald2D: fine spacing (fast) - use 0.1A or finer
    fine_spacing = min(0.1, coarse_spacing)
    z_scan_fine = np.arange(z_lo, z_hi + 1e-9, fine_spacing)
    x_scan_fine = np.arange(x_min, x_min + Lx_cell + 1e-9, fine_spacing)
    
    if verbose:
        print(f"  Brute force sampling: {coarse_spacing}Å ({len(z_scan_brute)} z-points, {len(x_scan_brute)} x-points)")
        print(f"  GridFF/Ewald2D sampling: {fine_spacing}Å ({len(z_scan_fine)} z-points, {len(x_scan_fine)} x-points)")

    # Avoid singularities in brute-force Coulomb sum: do not evaluate exactly on ions
    # (division-by-zero). We choose scan points away from known ionic sites.
    eps_xy = 1e-3
    x_hollow = x_min + 0.25 * Lx_cell
    y_hollow = y_min + 0.25 * Ly_cell
    
    # Debug: Print grid alignment info
    if verbose:
        print(f"\n  Grid alignment check:")
        print(f"    Atomic z_top: {z_top:.3f} Å")
        print(f"    GridFF g0[2] (z-origin): {g0[2]:.3f} Å")
        print(f"    GridFF z-extent: {g0[2]} to {g0[2] + nz*dg[2]:.3f} Å")
        print(f"    Sampling z-range (vacuum only): {z_lo:.3f} to {z_hi:.3f} Å")
        print(f"    Sampling x-range: {x_scan_fine[0]:.3f} to {x_scan_fine[-1]:.3f} Å")
        if (z_lo < g0[2]) or (z_hi > (g0[2] + nz*dg[2])):
            print(f"WARNING: z sampling range not fully inside GridFF z-extent")

    # Scan locations - each entry: (name, x0, y0, z_brute_arr, z_fine_arr, z0_for_xscan)
    # For z-scans: use arrays; for x-scan: z0 is scalar height
    # Based on NaCl_1x1_L3.xyz geometry:
    #   Na at (0,0,-3.25), (0,0,-8.91), (2,2,-6.08)
    #   Cl at (0,0,-6.08), (2,2,-3.25)
    # Use small offset to avoid exact singularity in brute force
    eps = 0.05
    scan_configs = [
        ('z_on_Na', 0.0 + eps, 0.0 + eps, z_scan_brute, z_scan_fine, None),       # Z-scan above Na at (0,0)
        ('z_on_Cl', 0.0 + eps, 0.0 + eps, z_scan_brute, z_scan_fine, None),     # Same XY but Cl is below
        ('z_midpoint', 1.0, 0.0, z_scan_brute, z_scan_fine, None),                # Z-scan at (1.0, 0.0) as suggested
    ]
    # Note: (0,0) has both Na and Cl at different z, so we scan same XY for both
    # The difference is just which atom is directly below at that XY

    # ==========================================================================
    # PART 1: Brute Force 1D Lines (Expensive - coarse sampling)
    # ==========================================================================
    if verbose:
        print(f"\n{'='*60}")
        print(f"Part 1: Brute Force 1D Lines (N_rep={N_rep})")
        print(f"{'='*60}")

    brute_results = {}
    ewald_1d_results = {}
    gridff_1d_results = {}

    for config in scan_configs:
        name = config[0]
        
        if name.startswith('z_'):
            # Z-scan: (name, x0, y0, z_brute_arr, z_fine_arr, None)
            _, x0, y0, z_brute_arr, z_fine_arr, _ = config
            
            # Brute force on coarse grid
            phi_brute = ew.phi_brute_1d(x0, y0, z_brute_arr, N_rep=N_rep)
            brute_results[name] = {'z': z_brute_arr, 'phi': phi_brute}
            
            # Ewald2D on fine grid
            phi_ewald = ew.phi_full_1d(x0, y0, z_fine_arr)
            ewald_1d_results[name] = {'z': z_fine_arr, 'phi': phi_ewald}
            
            # GridFF on fine grid
            x_arr = np.full_like(z_fine_arr, x0)
            y_arr = np.full_like(z_fine_arr, y0)
            positions = np.column_stack([x_arr, y_arr, z_fine_arr])
            result = sample_gridff_opencl(gridff_path, sub_xyz_path, positions, 
                                          grid_p0=g0, grid_step=dg,
                                          atom_req=(0.0, 0.0, 1.0, 0.0))
            phi_gridff = result['energies']
            gridff_1d_results[name] = {'z': z_fine_arr, 'phi': phi_gridff}
            
            if verbose:
                print(f"  {name}: ({x0:.2f}, {y0:.2f}) z=[{z_brute_arr[0]:.2f}, {z_brute_arr[-1]:.2f}], "
                      f"brute_n={len(z_brute_arr)}, fine_n={len(z_fine_arr)}")
        
        elif name == 'x_scan':
            # X-scan: (name, x_brute_arr, x_fine_arr, None, None, z0)
            _, x_brute_arr, x_fine_arr, _, _, z0 = config
            
            # Brute force on coarse grid
            phi_brute = np.zeros(len(x_brute_arr))
            for i, x in enumerate(x_brute_arr):
                phi_brute[i] = ew.phi_brute_1d(x, y_hollow, np.array([z0]), N_rep=N_rep)[0]
            brute_results[name] = {'x': x_brute_arr, 'phi': phi_brute}
            
            # Ewald2D on fine grid
            phi_ewald = np.zeros(len(x_fine_arr))
            for i, x in enumerate(x_fine_arr):
                phi_ewald[i] = ew.phi_full_1d(x, y_hollow, np.array([z0]))[0]
            ewald_1d_results[name] = {'x': x_fine_arr, 'phi': phi_ewald}
            
            # GridFF on fine grid
            y_arr = np.full_like(x_fine_arr, y_hollow)
            z_arr = np.full_like(x_fine_arr, z0)
            positions = np.column_stack([x_fine_arr, y_arr, z_arr])
            result = sample_gridff_opencl(gridff_path, sub_xyz_path, positions,
                                          grid_p0=g0, grid_step=dg,
                                          atom_req=(0.0, 0.0, 1.0, 0.0))
            phi_gridff = result['energies']
            gridff_1d_results[name] = {'x': x_fine_arr, 'phi': phi_gridff}
            
            if verbose:
                print(f"  {name}: y={y_hollow:.2f}, z={z0:.2f}, x=[{x_brute_arr[0]:.2f}, {x_brute_arr[-1]:.2f}], "
                      f"brute_n={len(x_brute_arr)}, fine_n={len(x_fine_arr)}")

    # ==========================================================================
    # PART 2: GridFF vs Ewald2D 2D Slices (Full Resolution)
    # ==========================================================================
    if verbose:
        print(f"\n{'='*60}")
        print(f"Part 2: GridFF vs Ewald2D 2D Slices")
        print(f"{'='*60}")
        print(f"  Resolution: {grid_res}x{grid_res}")
        print(f"  XY heights: {[f'{z:.2f}' for z in z_scan_heights]}")

    # XY slices: use GridFF g0 as origin so both methods sample the same spatial region
    xv = np.linspace(g0[0], g0[0] + Lx_cell, grid_res)
    yv = np.linspace(g0[1], g0[1] + Ly_cell, grid_res)
    X_xy, Y_xy = np.meshgrid(xv, yv)
    gridff_xy_slices = []
    ewald_xy_slices = []

    for zh in z_scan_heights:
        # Ewald2D vacuum XY (correct formula for z > all ions)
        phi_ewald_xy = ew.phi_vacuum_xy(X_xy, Y_xy, zh)
        ewald_xy_slices.append(phi_ewald_xy)

        # GridFF: sample at all grid points (electrostatic: Q=1, others=0)
        positions_xy = np.column_stack([X_xy.ravel(), Y_xy.ravel(), np.full(X_xy.size, zh)])
        result_xy = sample_gridff_opencl(gridff_path, sub_xyz_path, positions_xy,
                                         grid_p0=g0, grid_step=dg,
                                         atom_req=(0.0, 0.0, 1.0, 0.0))
        phi_gridff_xy = result_xy['energies'].reshape(X_xy.shape)
        gridff_xy_slices.append(phi_gridff_xy)

    # XZ slice at y = g0[1] + 0.5*Ly_cell (consistent with GridFF origin)
    y_fixed = g0[1] + 0.5 * Ly_cell
    xv = np.linspace(g0[0], g0[0] + Lx_cell, grid_res)
    zv = np.linspace(z_lo, z_hi, grid_res)
    X_xz, Z_xz = np.meshgrid(xv, zv)
    Y_xz = np.full_like(X_xz, y_fixed)

    # Ewald2D: compute row by row (each row = same z, use phi_vacuum_xy for consistency)
    phi_ewald_xz = np.zeros_like(X_xz)
    for iz, z_row in enumerate(zv):
        phi_ewald_xz[iz, :] = ew.phi_full_1d(xv, np.full_like(xv, y_fixed), np.array([z_row]))[0] if False else ew.phi_vacuum_xy(X_xz[iz:iz+1,:], Y_xz[iz:iz+1,:], z_row)[0, :]

    # GridFF: sample at all grid points (electrostatic: Q=1, others=0)
    positions_xz = np.column_stack([X_xz.ravel(), Y_xz.ravel(), Z_xz.ravel()])
    result_xz = sample_gridff_opencl(gridff_path, sub_xyz_path, positions_xz,
                                     grid_p0=g0, grid_step=dg,
                                     atom_req=(0.0, 0.0, 1.0, 0.0))
    phi_gridff_xz = result_xz['energies'].reshape(X_xz.shape)

    if verbose:
        print(f"  XY slices computed: {len(z_scan_heights)}")
        print(f"  XZ slice computed at y={Ly/2:.2f}")

    # ==========================================================================
    # PART 3: Error Statistics
    # ==========================================================================
    if verbose:
        print(f"\n{'='*60}")
        print(f"Part 3: Error Statistics")
        print(f"{'='*60}")

    stats = {}

    # 1D line comparisons
    for name in brute_results.keys():
        phi_brute = brute_results[name]['phi']
        coord_brute = brute_results[name].get('z', brute_results[name].get('x'))
        
        phi_ewald_fine = ewald_1d_results[name]['phi']
        coord_ewald_fine = ewald_1d_results[name].get('z', ewald_1d_results[name].get('x'))
        
        phi_gridff_fine = gridff_1d_results[name]['phi']
        
        # Interpolate fine-grid results to coarse grid for comparison with brute force
        phi_ewald_interp = np.interp(coord_brute, coord_ewald_fine, phi_ewald_fine)
        phi_gridff_interp = np.interp(coord_brute, coord_ewald_fine, phi_gridff_fine)
        
        # Filter: only compute errors for points sufficiently far from surface
        # (methods don't converge very close to ions)
        if 'z' in brute_results[name]:
            mask = coord_brute >= (z_top + z_min_for_error)
            coord_brute_filtered = coord_brute[mask]
            phi_brute_filtered = phi_brute[mask]
            phi_ewald_interp_filtered = phi_ewald_interp[mask]
            phi_gridff_interp_filtered = phi_gridff_interp[mask]
            coord_ewald_fine_filtered = coord_ewald_fine[coord_ewald_fine >= (z_top + z_min_for_error)]
            phi_ewald_fine_filtered = phi_ewald_fine[coord_ewald_fine >= (z_top + z_min_for_error)]
            phi_gridff_fine_filtered = phi_gridff_fine[coord_ewald_fine >= (z_top + z_min_for_error)]
        else:
            # For x_scan, use all points
            coord_brute_filtered = coord_brute
            phi_brute_filtered = phi_brute
            phi_ewald_interp_filtered = phi_ewald_interp
            phi_gridff_interp_filtered = phi_gridff_interp
            coord_ewald_fine_filtered = coord_ewald_fine
            phi_ewald_fine_filtered = phi_ewald_fine
            phi_gridff_fine_filtered = phi_gridff_fine
        
        err_ewald = phi_ewald_interp_filtered - phi_brute_filtered
        err_gridff = phi_gridff_interp_filtered - phi_brute_filtered

        stats[name] = {
            'ewald_vs_brute': {'rmse': float(np.sqrt(np.mean(err_ewald**2))),
                              'max_err': float(np.max(np.abs(err_ewald)))},
            'gridff_vs_brute': {'rmse': float(np.sqrt(np.mean(err_gridff**2))),
                               'max_err': float(np.max(np.abs(err_gridff)))},
            'gridff_vs_ewald': {'rmse': float(np.sqrt(np.mean((phi_gridff_fine_filtered - phi_ewald_fine_filtered)**2))),
                               'max_err': float(np.max(np.abs(phi_gridff_fine_filtered - phi_ewald_fine_filtered)))}
        }

        if verbose:
            print(f"\n  {name}:")
            print(f"    Ewald2D vs brute:    RMSE={stats[name]['ewald_vs_brute']['rmse']:.4e}, max_err={stats[name]['ewald_vs_brute']['max_err']:.4e}")
            print(f"    GridFF vs brute:     RMSE={stats[name]['gridff_vs_brute']['rmse']:.4e}, max_err={stats[name]['gridff_vs_brute']['max_err']:.4e}")
            print(f"    GridFF vs Ewald2D:   RMSE={stats[name]['gridff_vs_ewald']['rmse']:.4e}, max_err={stats[name]['gridff_vs_ewald']['max_err']:.4e}")

    # 2D slice comparisons
    for i, zh in enumerate(z_scan_heights):
        phi_g = gridff_xy_slices[i]
        phi_e = ewald_xy_slices[i]
        err = phi_g - phi_e
        stats[f'xy_z{zh:.1f}'] = {
            'gridff_vs_ewald': {'rmse': float(np.sqrt(np.mean(err**2))),
                               'max_err': float(np.max(np.abs(err)))}
        }

    err_xz = phi_gridff_xz - phi_ewald_xz
    stats['xz'] = {
        'gridff_vs_ewald': {'rmse': float(np.sqrt(np.mean(err_xz**2))),
                           'max_err': float(np.max(np.abs(err_xz)))}
    }

    # ==========================================================================
    # PART 4: Plotting
    # ==========================================================================
    if verbose:
        print(f"\n{'='*60}")
        print(f"Part 4: Generating Plots")
        print(f"{'='*60}")

    # Figure 1: Brute Force 1D Lines with all three methods
    fig1, axes1 = plt.subplots(1, 3, figsize=(15, 4))
    fig1.suptitle(f'Electrostatics 1D Line Scans (N_rep={N_rep})', fontsize=12)

    for ax, config in zip(axes1, scan_configs):
        name = config[0]
        _, x0, y0, _, _, _ = config
        z_brute = brute_results[name]['z']
        z_fine = ewald_1d_results[name]['z']
        ax.plot(z_brute, brute_results[name]['phi'], 'ko', ms=4, label='Brute force')
        ax.plot(z_fine, ewald_1d_results[name]['phi'], 'r-', lw=1, label='Ewald2D')
        ax.plot(z_fine, gridff_1d_results[name]['phi'], 'b--', lw=1, label='GridFF')
        ax.set_xlabel('z (Å)')
        ax.set_ylabel('φ (eV)')
        ax.set_title(f'{name}\n({x0:.2f}, {y0:.2f})')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    fig1.tight_layout()
    fig1_path = os.path.join(save_dir, f'{prefix}_fig1_1d_lines.png')
    fig1.savefig(fig1_path, dpi=150)
    plt.close(fig1)
    if verbose:
        print(f"  Saved {fig1_path}")

    # Figure 2: GridFF vs Ewald2D XY Slices
    n_heights = len(z_scan_heights)
    fig2, axes2 = plt.subplots(2, n_heights, figsize=(5*n_heights, 8))
    if n_heights == 1:
        axes2 = axes2[:, None]

    for i, zh in enumerate(z_scan_heights):
        # GridFF
        vmax = np.max(np.abs([gridff_xy_slices[i], ewald_xy_slices[i]]))
        im0 = axes2[0, i].imshow(gridff_xy_slices[i], extent=[0, Lx, 0, Ly],
                                  origin='lower', cmap='RdBu_r', vmin=-vmax, vmax=vmax)
        axes2[0, i].set_title(f'GridFF z={zh:.1f}Å')
        axes2[0, i].set_xlabel('x (Å)')
        axes2[0, i].set_ylabel('y (Å)')
        plt.colorbar(im0, ax=axes2[0, i], label='φ (eV)')

        # Ewald2D
        im1 = axes2[1, i].imshow(ewald_xy_slices[i], extent=[0, Lx, 0, Ly],
                                  origin='lower', cmap='RdBu_r', vmin=-vmax, vmax=vmax)
        axes2[1, i].set_title(f'Ewald2D z={zh:.1f}Å')
        axes2[1, i].set_xlabel('x (Å)')
        axes2[1, i].set_ylabel('y (Å)')
        plt.colorbar(im1, ax=axes2[1, i], label='φ (eV)')

    fig2.tight_layout()
    fig2_path = os.path.join(save_dir, f'{prefix}_fig2_xy_slices.png')
    fig2.savefig(fig2_path, dpi=150)
    plt.close(fig2)
    if verbose:
        print(f"  Saved {fig2_path}")

    # Figure 3: GridFF vs Ewald2D XZ Slice with difference
    fig3, axes3 = plt.subplots(1, 3, figsize=(15, 5))

    vmax = 1.0  # Fixed symmetric color limit
    im0 = axes3[0].pcolormesh(X_xz, Z_xz, phi_gridff_xz, shading='auto',
                               cmap='RdBu_r', vmin=-vmax, vmax=vmax)
    axes3[0].set_title('GridFF XZ')
    axes3[0].set_xlabel('x (Å)')
    axes3[0].set_ylabel('z (Å)')
    plt.colorbar(im0, ax=axes3[0], label='φ (eV)')

    im1 = axes3[1].pcolormesh(X_xz, Z_xz, phi_ewald_xz, shading='auto',
                               cmap='RdBu_r', vmin=-vmax, vmax=vmax)
    axes3[1].set_title('Ewald2D XZ')
    axes3[1].set_xlabel('x (Å)')
    axes3[1].set_ylabel('z (Å)')
    plt.colorbar(im1, ax=axes3[1], label='φ (eV)')

    vmax_err = 1.0  # Fixed symmetric error color limit
    im2 = axes3[2].pcolormesh(X_xz, Z_xz, err_xz, shading='auto',
                               cmap='RdBu_r', vmin=-vmax_err, vmax=vmax_err)
    axes3[2].set_title(f'Difference (GridFF - Ewald2D)')
    axes3[2].set_xlabel('x (Å)')
    axes3[2].set_ylabel('z (Å)')
    plt.colorbar(im2, ax=axes3[2], label='Δφ (eV)')

    fig3.tight_layout()
    fig3_path = os.path.join(save_dir, f'{prefix}_fig3_xz_slice.png')
    fig3.savefig(fig3_path, dpi=150)
    plt.close(fig3)
    if verbose:
        print(f"  Saved {fig3_path}")

    # Figure 4: Error plots for 1D lines
    fig4, axes4 = plt.subplots(2, 3, figsize=(15, 8))
    fig4.suptitle('1D Line Errors vs Brute Force Reference (interpolated)', fontsize=12)

    for col, config in enumerate(scan_configs):
        name = config[0]
        _, x0, y0, _, _, _ = config
        z_brute = brute_results[name]['z']
        z_fine = ewald_1d_results[name]['z']
        # Interpolate to brute grid for comparison
        phi_ewald_interp = np.interp(z_brute, z_fine, ewald_1d_results[name]['phi'])
        phi_gridff_interp = np.interp(z_brute, z_fine, gridff_1d_results[name]['phi'])
        # Top row: Ewald2D - brute
        axes4[0, col].plot(z_brute, phi_ewald_interp - brute_results[name]['phi'], 'r-', lw=1)
        axes4[0, col].set_title(f'{name}: Ewald2D error')
        axes4[0, col].set_xlabel('z (Å)')
        axes4[0, col].set_ylabel('Δφ (e/Å)')
        axes4[0, col].axhline(0, color='k', lw=0.5)
        axes4[0, col].grid(True, alpha=0.3)

        # Bottom row: GridFF - brute
        axes4[1, col].plot(z_brute, phi_gridff_interp - brute_results[name]['phi'], 'b-', lw=1)
        axes4[1, col].set_title(f'{name}: GridFF error')
        axes4[1, col].set_xlabel('z (Å)')
        axes4[1, col].set_ylabel('Δφ (e/Å)')
        axes4[1, col].axhline(0, color='k', lw=0.5)
        axes4[1, col].grid(True, alpha=0.3)

    fig4.tight_layout()
    fig4_path = os.path.join(save_dir, f'{prefix}_fig4_1d_errors.png')
    fig4.savefig(fig4_path, dpi=150)
    plt.close(fig4)
    if verbose:
        print(f"  Saved {fig4_path}")

    # Save JSON report
    report = {
        'prefix': prefix,
        'N_rep': N_rep,
        'coarse_spacing': coarse_spacing,
        'grid_res': grid_res,
        'z_scan_heights': [float(z) for z in z_scan_heights],
        'stats': stats,
        'figures': {
            'fig1_1d_lines': fig1_path,
            'fig2_xy_slices': fig2_path,
            'fig3_xz_slice': fig3_path,
            'fig4_1d_errors': fig4_path
        }
    }
    report_path = os.path.join(save_dir, f'{prefix}_report.json')
    with open(report_path, 'w') as f:
        json.dump(report, f, indent=2)
    if verbose:
        print(f"  Saved {report_path}")

    if verbose:
        print(f"\n{'='*60}")
        print(f"Comparison Complete")
        print(f"{'='*60}")

    return report


# ======================================================================
# OpenCL Ewald Comparison Functions
# ======================================================================

def compare_ewald_opencl_python(sys_at, n_harm=4, verbose=True):
    """
    Compare OpenCL Ewald implementation against Python reference.
    
    Parameters:
        sys_at: AtomicSystem with positions, charges, and lattice vectors
        n_harm: Ewald harmonic truncation
        verbose: print progress
        
    Returns:
        dict with comparison results
    """
    import time
    import numpy as np
    
    # Import OpenCL Ewald
    try:
        from pyBall.OCL.SurfaceEwald import SurfaceEwaldCL
        from pyBall.Ewald2D import Ewald2D
    except ImportError as e:
        if verbose:
            print(f"Cannot run OpenCL comparison: {e}")
        return {'error': str(e)}
    
    if verbose:
        print(f"\n{'='*60}")
        print(f"OpenCL vs Python Ewald Comparison")
        print(f"{'='*60}")
    
    # Initialize Python Ewald
    t0 = time.time()
    ew_py = Ewald2D.from_AtomicSystem(sys_at, n_harm=n_harm)
    t1 = time.time()
    if verbose:
        print(f"Python Ewald init: {t1-t0:.3f} s")
    
    # Initialize OpenCL Ewald
    t0 = time.time()
    ew_cl = SurfaceEwaldCL()
    ion_data = np.column_stack([
        sys_at.apos[:, 0],
        sys_at.apos[:, 1],
        sys_at.apos[:, 2],
        sys_at.qs
    ])
    a_vec = sys_at.lvec[0, :2]
    b_vec = sys_at.lvec[1, :2]
    ew_cl.prepare_system(ion_data, a_vec, b_vec, n_harm=n_harm)
    t1 = time.time()
    if verbose:
        print(f"OpenCL Ewald init: {t1-t0:.3f} s")
    
    results = {}
    
    # Test 1: Vacuum evaluation
    if verbose:
        print(f"\n{'='*60}")
        print(f"Test 1: Vacuum Evaluation")
        print(f"{'='*60}")
    
    xv = np.linspace(0, a_vec[0], 50)
    yv = np.linspace(0, b_vec[1], 50)
    X, Y = np.meshgrid(xv, yv)
    z_test = 2.0
    
    # OpenCL
    t0 = time.time()
    phi_cl = ew_cl.eval_vacuum(X, Y, z_test)
    t1 = time.time()
    time_cl = t1 - t0
    
    # Python
    t0 = time.time()
    phi_py = ew_py.phi_vacuum_xy(X, Y, z_test)
    t1 = time.time()
    time_py = t1 - t0
    
    # Compare
    diff = phi_cl - phi_py
    rmse = float(np.sqrt(np.mean(diff**2)))
    max_err = float(np.max(np.abs(diff)))
    
    results['vacuum'] = {
        'rmse': rmse,
        'max_err': max_err,
        'time_cl': time_cl,
        'time_py': time_py,
        'speedup': time_py / time_cl if time_cl > 0 else float('inf'),
        'N_points': X.size
    }
    
    if verbose:
        print(f"  Grid: {X.shape[0]}x{X.shape[1]} = {X.size} points")
        print(f"  OpenCL time: {time_cl:.3f} s")
        print(f"  Python time: {time_py:.3f} s")
        print(f"  Speedup: {results['vacuum']['speedup']:.1f}x")
        print(f"  RMSE: {rmse:.6e} eV")
        print(f"  Max error: {max_err:.6e} eV")
        if rmse < 1e-5:
            print(f"  ✓ PASS")
        else:
            print(f"  ✗ FAIL")
    
    # Test 2: Full evaluation
    if verbose:
        print(f"\n{'='*60}")
        print(f"Test 2: Full Evaluation (1D line)")
        print(f"{'='*60}")
    
    x0, y0 = 0.5, 0.5
    z_arr = np.linspace(-0.5, 5.0, 100)
    
    X_line = np.full((1, len(z_arr)), x0, dtype=np.float32)
    Y_line = np.full((1, len(z_arr)), y0, dtype=np.float32)
    Z_line = z_arr.reshape(1, -1).astype(np.float32)
    
    # OpenCL
    t0 = time.time()
    phi_cl_line = ew_cl.eval_full(X_line, Y_line, Z_line)[0, :]
    t1 = time.time()
    time_cl = t1 - t0
    
    # Python
    t0 = time.time()
    phi_py_line = ew_py.phi_full_1d(x0, y0, z_arr)
    t1 = time.time()
    time_py = t1 - t0
    
    # Compare
    diff_line = phi_cl_line - phi_py_line
    rmse_line = float(np.sqrt(np.mean(diff_line**2)))
    max_err_line = float(np.max(np.abs(diff_line)))
    
    results['full'] = {
        'rmse': rmse_line,
        'max_err': max_err_line,
        'time_cl': time_cl,
        'time_py': time_py,
        'speedup': time_py / time_cl if time_cl > 0 else float('inf'),
        'N_points': len(z_arr)
    }
    
    if verbose:
        print(f"  Points: {len(z_arr)}")
        print(f"  OpenCL time: {time_cl:.3f} s")
        print(f"  Python time: {time_py:.3f} s")
        print(f"  Speedup: {results['full']['speedup']:.1f}x")
        print(f"  RMSE: {rmse_line:.6e} eV")
        print(f"  Max error: {max_err_line:.6e} eV")
        if rmse_line < 1e-5:
            print(f"  ✓ PASS")
        else:
            print(f"  ✗ FAIL")
    
    # Overall result
    results['pass'] = (results['vacuum']['rmse'] < 1e-5 and 
                       results['full']['rmse'] < 1e-5)
    
    if verbose:
        print(f"\n{'='*60}")
        print(f"Overall: {'PASS' if results['pass'] else 'FAIL'}")
        print(f"{'='*60}")
    
    return results


def compare_all_methods(sys_at, gridff_path, sub_xyz_path, n_harm=4, N_rep=20, verbose=True):
    """
    Comprehensive comparison of all electrostatics methods:
    - GridFF (OpenCL B-spline)
    - Ewald2D Python
    - Ewald2D OpenCL
    - Brute Force (reference)
    
    Parameters:
        sys_at: AtomicSystem
        gridff_path: path to GridFF .npy file
        sub_xyz_path: path to substrate .xyz file
        n_harm: Ewald harmonic truncation
        N_rep: Brute force PBC shells
        verbose: print progress
        
    Returns:
        dict with all comparison results
    """
    all_results = {}
    
    # Run GridFF vs Ewald2D comparison
    from pyBall.Ewald2D import Ewald2D
    ew = Ewald2D.from_AtomicSystem(sys_at, n_harm=n_harm)
    
    report = compare_electrostatics_methods(
        sys_at=sys_at,
        gridff_path=gridff_path,
        sub_xyz_path=sub_xyz_path,
        ew=ew,
        save_dir='results_electrostatics',
        prefix='all_methods',
        N_rep=N_rep,
        verbose=verbose
    )
    all_results['gridff_vs_ewald'] = report
    
    # Run OpenCL vs Python comparison
    opencl_results = compare_ewald_opencl_python(sys_at, n_harm=n_harm, verbose=verbose)
    all_results['opencl_vs_python'] = opencl_results
    
    # Summary
    if verbose:
        print(f"\n{'='*60}")
        print(f"Comprehensive Comparison Summary")
        print(f"{'='*60}")
        print(f"\n1. GridFF vs Python Ewald:")
        if 'gridff_vs_ewald' in report['stats']:
            stats = report['stats']['gridff_vs_ewald']
            print(f"   GridFF vs Ewald: RMSE={stats['rmse']:.4e}, max={stats['max_err']:.4e}")
        
        print(f"\n2. OpenCL Ewald vs Python Ewald:")
        if 'vacuum' in opencl_results:
            print(f"   Vacuum: RMSE={opencl_results['vacuum']['rmse']:.4e}")
            print(f"   Full:   RMSE={opencl_results['full']['rmse']:.4e}")
            print(f"   Speedup: {opencl_results['vacuum']['speedup']:.1f}x (vacuum)")
        
        all_pass = (
            opencl_results.get('pass', False) and
            all(s['gridff_vs_ewald']['rmse'] < 1e-5 
                for key, s in report['stats'].items() 
                if isinstance(s, dict) and 'gridff_vs_ewald' in s)
        )
        print(f"\nOverall: {'PASS' if all_pass else 'FAIL'}")
        print(f"{'='*60}")
    
    return all_results
