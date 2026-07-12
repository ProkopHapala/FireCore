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


def plot_alignment_summary(grid_data, g0, dg, atoms_xyz, atoms_enames, save_path, iz_top=None, iy_center=None, z_atom_range=2.0, mol_apos=None, mol_enames=None, plq_coeffs=None, plq_coeffs2=None, zmin_offset=2.0, z_ylim=None, vmin_vmax_xz=None, plot_diagnostics=True, plot_mol_samples=True, xy_same_scale_as_xz=False, z_profile_range=None, elem_name='H', elem_name2='O', H_H=0.0, H_O=0.0):
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
        plot_mol_samples: If False, skip plotting molecule sample atoms in 2D imshow plots (default: True)
        xy_same_scale_as_xz: If True, use same color scale for XY plot as XZ plot (default: False)
        z_profile_range: Optional tuple (zmin, zmax) for Z-profile x-axis range in Angstrom (default: 0 to 6 A)
        elem_name: Element name for first plq_coeffs (default: 'H')
        elem_name2: Element name for second plq_coeffs2 (default: 'O')
        H_H: H-bond correction coefficient for H (default: 0.0)
        H_O: H-bond correction coefficient for O (default: 0.0)
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
    # XY slice - symmetric range from this slice only (or use XZ scale if requested)
    xy_slice_data = total_potential[:, :, iz_top, 0]
    if xy_same_scale_as_xz:
        # Use XZ scale for XY (will be computed below)
        vmin_xy, vmax_xy = None, None
    else:
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
        if plot_mol_samples and mol_apos_flat is not None:
            ax.scatter(mol_apos_flat[:, 0], mol_apos_flat[:, 2],
                       c=mol_colors, s=15, alpha=0.6, marker='.', edgecolors='black', linewidth=0.5, label='mol samples')
        ax.set_title(f'Total Potential\n(P={P:.3f},L={L:.3f},Q={Q:.3f})\nXZ at y={y_val:.3f}A (iy={iy}) {title_suffix}')
        ax.set_xlabel('x [A]')
        ax.set_ylabel('z [A]')
        ax.set_ylim(zmin_display, g0[2] + min(nz*dz, 20))  # Show from z=0 (surface)
        plt.colorbar(im, ax=ax, shrink=0.8)

    # 1. XY slice at top layer - Total potential
    ax = axes[0, 0]
    z_val = g0[2] + iz_top * dz
    extent = [g0[0], g0[0] + nx*dx, g0[1], g0[1] + ny*dy]
    # Use XZ scale for XY if requested
    vmin_xy_use, vmax_xy_use = (vmin_xz, vmax_xz) if xy_same_scale_as_xz else (vmin_xy, vmax_xy)
    im = ax.imshow(total_potential[:, :, iz_top, 0].T, extent=extent,
                   origin='lower', cmap='bwr', aspect='equal', vmin=vmin_xy_use, vmax=vmax_xy_use)
    plot_atoms_overlay(ax, atoms_xyz, atoms_enames, z_atom_range=z_atom_range, g0=g0, dg=dg)
    # Overlay molecule samples
    if plot_mol_samples and mol_apos_flat is not None:
        ax.scatter(mol_apos_flat[:, 0], mol_apos_flat[:, 1], c=mol_colors, s=15, alpha=0.6, marker='.', edgecolors='black', linewidth=0.5, label='mol samples')
    ax.set_title(f'Total Potential\n(P={P:.3f},L={L:.3f},Q={Q:.3f})\nXY at z={z_val:.3f}A (iz={iz_top})')
    ax.set_xlabel('x [A]')
    ax.set_ylabel('y [A]')
    plt.colorbar(im, ax=ax, shrink=0.8)

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
            # Get H-bond coefficient for this element
            H_val = H_H if elem_name == 'H' else H_O
            ax.plot(z_coords, total_potential[cx, cy, :, 0], 'k-', linewidth=2, label=f'{elem_name}-Total (P={P:.3f},L={L:.3f},Q={Q:.3f},H={H_val:.1f})')
            if plq_coeffs2 is not None:
                P2, L2, Q2 = plq_coeffs2[0], plq_coeffs2[1], plq_coeffs2[2]
                total2 = P2*grid_data[...,0] + L2*grid_data[...,1] + Q2*grid_data[...,2]
                H_val2 = H_H if elem_name2 == 'H' else H_O
                ax.plot(z_coords, total2[cx, cy, :], 'm-', linewidth=2, label=f'{elem_name2}-Total (P={P2:.3f},L={L2:.3f},Q={Q2:.3f},H={H_val2:.1f})')
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
            if z_profile_range is not None:
                ax.set_xlim(z_profile_range[0], z_profile_range[1])
            else:
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


# ============================================================
# FDBM Fitting Utilities — shared between batch and interactive tools
# ============================================================

def sample_gridff_trilinear(pos, gridff, g0, dg):
    """Sample GridFF numpy array at position (x, y, z) using trilinear interpolation.
    - x/y: periodic wrap
    - z: clamp
    Returns: vals (nch,) float32, clamped_z bool
    """
    nx, ny, nz, nch = gridff.shape
    fx = (pos[0] - g0[0]) / dg[0]
    fy = (pos[1] - g0[1]) / dg[1]
    fz = (pos[2] - g0[2]) / dg[2]
    ix0 = int(np.floor(fx)); iy0 = int(np.floor(fy)); iz0 = int(np.floor(fz))
    tx = fx - ix0; ty = fy - iy0; tz = fz - iz0
    ixx = ix0 % nx;  iyy = iy0 % ny
    clamped_z = (iz0 < 0) or (iz0 >= nz - 1)
    iz = max(0, min(iz0, nz - 2));  iz1 = iz + 1
    ix1 = (ixx + 1) % nx;  iy1 = (iyy + 1) % ny
    c000 = gridff[ixx, iyy, iz , :]; c100 = gridff[ix1, iyy, iz , :]
    c010 = gridff[ixx, iy1, iz , :]; c110 = gridff[ix1, iy1, iz , :]
    c001 = gridff[ixx, iyy, iz1, :]; c101 = gridff[ix1, iyy, iz1, :]
    c011 = gridff[ixx, iy1, iz1, :]; c111 = gridff[ix1, iy1, iz1, :]
    c00 = c000*(1.0-tx) + c100*tx; c10 = c010*(1.0-tx) + c110*tx
    c01 = c001*(1.0-tx) + c101*tx; c11 = c011*(1.0-tx) + c111*tx
    c0  = c00 *(1.0-ty) + c10 *ty; c1  = c01 *(1.0-ty) + c11 *ty
    vals = c0  *(1.0-tz) + c1  *tz
    return vals, clamped_z


def load_dft_scan_data(dft_paths, z_top_sub):
    """Load and concatenate DFT matched NPZ datasets, compute interaction energies.
    Args:
        dft_paths: list of paths to *_matched.npz files
        z_top_sub: z coordinate of top substrate atom (for coordinate frame shift)
    Returns: dict with keys: coords, energies, enames, ix, iy, orientxy, orientz, zdist, mol_tag, E_int, E0, weights (all zero initially)
    """
    dsets = []
    for p in dft_paths:
        if os.path.exists(p):
            dsets.append((os.path.basename(p).replace('_matched.npz', ''), np.load(p)))
    if len(dsets) == 0:
        raise RuntimeError(f'No DFT matched datasets found at {dft_paths}')
    coords   = np.concatenate([ds['coords']   for _, ds in dsets], axis=0).copy()
    energies = np.concatenate([ds['energies'] for _, ds in dsets], axis=0)
    enames   = dsets[0][1]['enames']
    ix_arr       = np.concatenate([ds['ix']       for _, ds in dsets], axis=0)
    iy_arr       = np.concatenate([ds['iy']       for _, ds in dsets], axis=0)
    orientxy = np.concatenate([ds['orientxy'] for _, ds in dsets], axis=0)
    orientz  = np.concatenate([ds['orientz']  for _, ds in dsets], axis=0)
    zdist    = np.concatenate([ds['zdist']    for _, ds in dsets], axis=0)
    mol_tag  = np.concatenate([np.full(len(ds['energies']), name, dtype=object) for name, ds in dsets], axis=0)
    coords[:, :, 2] -= z_top_sub
    # Baseline subtraction per group
    E0 = np.zeros_like(energies, dtype=np.float64)
    uniq_keys = {}
    for i in range(len(energies)):
        k = (str(mol_tag[i]), int(ix_arr[i]), int(iy_arr[i]), str(orientxy[i]), int(orientz[i]))
        uniq_keys.setdefault(k, []).append(i)
    for k, idxs in uniq_keys.items():
        zz = zdist[idxs]
        zmax = float(np.max(zz))
        mask = np.isclose(zz, zmax, rtol=0.0, atol=1e-6)
        imax = np.array(idxs, dtype=int)[mask]
        E0[idxs] = float(np.mean(energies[imax]))
    E_int = energies - E0
    return dict(coords=coords, energies=energies, enames=enames, ix=ix_arr, iy=iy_arr,
                orientxy=orientxy, orientz=orientz, zdist=zdist, mol_tag=mol_tag,
                E_int=E_int, E0=E0, dset_names=[n for n, _ in dsets])


def compute_fit_weights(E_int, mol_tag, ix, iy, orientz, fit_scan_configs, scan_tilt=0, Eint_max=0.5, kT_weight=0.2):
    """Compute exponential fit weights. Returns weights array (same length as E_int), zero for excluded points."""
    mask_scan = np.zeros(len(E_int), dtype=bool)
    for mol_f, ix_f, iy_f in fit_scan_configs:
        mask_scan |= (mol_tag == mol_f) & (ix == ix_f) & (iy == iy_f) & (orientz == scan_tilt)
    mask_fit = mask_scan & (E_int <= Eint_max)
    weights = np.exp(-E_int / kT_weight)
    weights[~mask_fit] = 0.0
    if np.any(mask_fit):
        weights[mask_fit] /= weights[mask_fit].max()
    return weights


def prepare_scan_panel_data(data, gridff, g0, dg, panel_mols, panel_ix, panel_iy, panel_sites,
                             orient_name, scan_tilt, dz_shift, idx_pauli=0, idx_coulomb=2, idx_polar=3, per_atom=False):
    """Pre-sample GridFF channels for each of the 4 scan panels.
    Args:
        per_atom: if True, return per-atom rho/phi/tau (rhoH1, rhoH2, rhoO, phiH1, phiH2, phiO, tauH1, tauH2, tauO)
                  if False, return summed (rhoH, rhoO, phiH, phiO, tauH, tauO)
    Returns list of dicts, one per panel.
    """
    panels = []
    for pidx in range(len(panel_mols)):
        mol_p = panel_mols[pidx]; ix_p = panel_ix[pidx]; iy_p = panel_iy[pidx]; site_label = panel_sites[pidx]
        mask = ((data['mol_tag'] == mol_p) & (data['ix'] == ix_p) & (data['iy'] == iy_p) &
                (data['orientxy'] == orient_name) & (data['orientz'] == scan_tilt))
        if np.count_nonzero(mask) == 0:
            panels.append(None); continue
        z_sel = data['zdist'][mask]
        E_sel = data['energies'][mask]
        c_sel = data['coords'][mask]
        w_sel = data['weights'][mask]
        zmax = float(np.max(z_sel))
        mask_zmax = np.isclose(z_sel, zmax, rtol=0.0, atol=1e-6)
        E0 = float(np.mean(E_sel[mask_zmax]))
        E_int_dft = E_sel - E0
        enames = data['enames']
        n_frames = len(z_sel)
        if per_atom:
            rhoH1 = np.zeros(n_frames); rhoH2 = np.zeros(n_frames); rhoO = np.zeros(n_frames)
            phiH1 = np.zeros(n_frames); phiH2 = np.zeros(n_frames); phiO = np.zeros(n_frames)
            tauH1 = np.zeros(n_frames); tauH2 = np.zeros(n_frames); tauO = np.zeros(n_frames)
            for i in range(n_frames):
                for j, ename in enumerate(enames):
                    pos = c_sel[i, j].copy(); pos[2] += dz_shift
                    vals, _ = sample_gridff_trilinear(pos, gridff, g0, dg)
                    if ename == 'H':
                        if j == 0:
                            rhoH1[i] += vals[idx_pauli]; phiH1[i] += vals[idx_coulomb]; tauH1[i] += vals[idx_polar]
                        else:
                            rhoH2[i] += vals[idx_pauli]; phiH2[i] += vals[idx_coulomb]; tauH2[i] += vals[idx_polar]
                    elif ename == 'O':
                        rhoO[i] += vals[idx_pauli]; phiO[i] += vals[idx_coulomb]; tauO[i] += vals[idx_polar]
            sort_idx = np.argsort(z_sel)
            panels.append(dict(
                z_s=z_sel[sort_idx], E_int_dft=E_int_dft[sort_idx], w_s=w_sel[sort_idx],
                rhoH1=rhoH1[sort_idx], rhoH2=rhoH2[sort_idx], rhoO=rhoO[sort_idx],
                phiH1=phiH1[sort_idx], phiH2=phiH2[sort_idx], phiO=phiO[sort_idx],
                tauH1=tauH1[sort_idx], tauH2=tauH2[sort_idx], tauO=tauO[sort_idx],
                mask_zmax=mask_zmax[sort_idx],
                apos=c_sel[sort_idx], names=list(enames),
                title=f'{mol_p} over {site_label} (ix={ix_p} iy={iy_p}, tilt={scan_tilt})'
            ))
        else:
            rhoH = np.zeros(n_frames); rhoO = np.zeros(n_frames)
            phiH = np.zeros(n_frames); phiO = np.zeros(n_frames)
            tauH = np.zeros(n_frames); tauO = np.zeros(n_frames)
            for i in range(n_frames):
                for j, ename in enumerate(enames):
                    pos = c_sel[i, j].copy(); pos[2] += dz_shift
                    vals, _ = sample_gridff_trilinear(pos, gridff, g0, dg)
                    if ename == 'H':
                        rhoH[i] += vals[idx_pauli]; phiH[i] += vals[idx_coulomb]; tauH[i] += vals[idx_polar]
                    elif ename == 'O':
                        rhoO[i] += vals[idx_pauli]; phiO[i] += vals[idx_coulomb]; tauO[i] += vals[idx_polar]
            sort_idx = np.argsort(z_sel)
            panels.append(dict(
                z_s=z_sel[sort_idx], E_int_dft=E_int_dft[sort_idx], w_s=w_sel[sort_idx],
                rhoH=rhoH[sort_idx], rhoO=rhoO[sort_idx], phiH=phiH[sort_idx], phiO=phiO[sort_idx],
                tauH=tauH[sort_idx], tauO=tauO[sort_idx],
                mask_zmax=mask_zmax[sort_idx],
                apos=c_sel[sort_idx], names=list(enames),
                title=f'{mol_p} over {site_label} (ix={ix_p} iy={iy_p}, tilt={scan_tilt})'
            ))
    return panels


def compute_model_Eint(panel, P_H, P_O, q_H, beta=1.0, H_H=0.0, H_O=0.0):
    """Compute model interaction energy for a pre-sampled panel (summed mode).
    E = P_H * rho_H^beta + P_O * rho_O^beta + H_H * tau_H + H_O * tau_O + E_coulomb
    Returns: E_int_model array (same length as panel['z_s'])
    """
    q_O = -2.0 * q_H
    Ec = q_H * panel['phiH'] + q_O * panel['phiO']
    rhoH_b = np.abs(panel['rhoH'])**beta
    rhoO_b = np.abs(panel['rhoO'])**beta
    tauH = panel.get('tauH', np.zeros_like(panel['rhoH']))
    tauO = panel.get('tauO', np.zeros_like(panel['rhoO']))
    Em = P_H * rhoH_b + P_O * rhoO_b + H_H * tauH + H_O * tauO + Ec
    Em0 = float(np.mean(Em[panel['mask_zmax']]))
    return Em - Em0


def compute_model_Eint_per_atom(panel, P_H, P_O, q_H, beta=1.0, H_H=0.0, H_O=0.0):
    """Compute per-atom model interaction energies (for decomposition view).
    Returns dict with keys: H1, H2, O (each an array)
    """
    q_O = -2.0 * q_H
    rhoH1_b = np.abs(panel['rhoH1'])**beta
    rhoH2_b = np.abs(panel['rhoH2'])**beta
    rhoO_b  = np.abs(panel['rhoO'])**beta
    tauH1 = panel.get('tauH1', np.zeros_like(panel['rhoH1']))
    tauH2 = panel.get('tauH2', np.zeros_like(panel['rhoH2']))
    tauO  = panel.get('tauO',  np.zeros_like(panel['rhoO']))
    EcH1 = q_H * panel['phiH1']
    EcH2 = q_H * panel['phiH2']
    EcO  = q_O * panel['phiO']
    EmH1 = P_H * rhoH1_b + H_H * tauH1 + EcH1
    EmH2 = P_H * rhoH2_b + H_H * tauH2 + EcH2
    EmO  = P_O * rhoO_b  + H_O * tauO  + EcO
    EmH1_0 = float(np.mean(EmH1[panel['mask_zmax']]))
    EmH2_0 = float(np.mean(EmH2[panel['mask_zmax']]))
    EmO_0  = float(np.mean(EmO[panel['mask_zmax']]))
    return dict(H1=EmH1-EmH1_0, H2=EmH2-EmH2_0, O=EmO-EmO_0)


def compute_cl_repulsive_on_grid(cl_apos, g0, dg, ns, lvec, nPBC=(4,4,0),
                                  R_eq=1.80, D_e=0.0116, alphaMorse=1.5, R_damp=0.1):
    """Compute Morse repulsive (Pauli) potential on 3D grid from Cl atoms only, using GPU via GridFF_cl.make_MorseFF.

    Reuses GridFF_cl.make_MorseFF (cpp/common_resources/cl/GridFF.cl) which evaluates:
        E_Paul(r) = sum_i D_e_i * exp(-2*alpha*(|r-r_i| - R_eq_i))
    over all atoms (Cl only here) with PBC. Only V_Paul is returned (London discarded).

    Args:
        cl_apos: (n_cl, 3) array of Cl atom positions in Angstrom
        g0: (3,) grid origin in Angstrom
        dg: (dx, dy, dz) grid spacing in Angstrom
        ns: (nx, ny, nz) grid shape
        lvec: (3, 3) lattice vectors [[ax,ay,az],[bx,by,bz],[cx,cy,cz]] in Angstrom
        nPBC: (npx, npy, npz) PBC replication counts (default (4,4,0) for 2D surface)
        R_eq: Morse equilibrium distance for Cl in Angstrom
        D_e: Morse well depth for Cl in eV
        alphaMorse: Morse exponent in 1/Angstrom (GFFParams.y)
        R_damp: damping radius for numerical stability (GFFParams.x)

    Returns:
        V_Paul: (nx, ny, nz) float32 array of Cl repulsive potential in eV
    """
    from pyBall.OCL.GridFF import GridFF_cl
    nx, ny, nz = ns
    n_cl = len(cl_apos)
    # atoms: (n, 4) float32 positions with w=0
    atoms = np.zeros((n_cl, 4), dtype=np.float32)
    atoms[:, :3] = cl_apos
    # REQs: (n, 4) float32 (R_eq, D_e, q=0, w=0) -- no Coulomb, pure Pauli
    REQs = np.zeros((n_cl, 4), dtype=np.float32)
    REQs[:, 0] = R_eq
    REQs[:, 1] = D_e
    # lvec as numpy array for GridShape (GridFF_cl will handle conversion internally)
    lvec_arr = np.asarray(lvec, dtype=np.float64) if not isinstance(lvec, np.ndarray) else lvec
    GFFParams = (R_damp, alphaMorse, 0.0, 0.0)
    gff = GridFF_cl()
    # make_MorseFF returns V_Paul, V_Lond both in shape ns[::-1] = (nz, ny, nx)
    V_Paul_zyx, _ = gff.make_MorseFF(
        atoms, REQs, nPBC=nPBC, dg=dg, ng=ns, lvec=lvec_arr, g0=tuple(g0), GFFParams=GFFParams
    )
    # Transpose from (nz, ny, nx) -> (nx, ny, nz) to match GridFF C-order convention
    V_Paul = np.ascontiguousarray(V_Paul_zyx.transpose(2, 1, 0))
    return V_Paul


# ============================================================
# Hybrid Potential: B-spline + Short-Range HardCore
# ============================================================
# Combines a low-resolution tricubic B-spline potential (smooth
# long-range electrostatics via erf split) with per-atom radial
# functions (Morse + erfc Coulomb core) for short-range.
#
# Decomposition (from SoftSplineHardAtomCore.chat.md §4):
#   V_total = V_smooth + V_core
#   V_smooth = Σ_j k_e * Q_j * erf(r_j/σ) / r_j   (long-range, on B-spline grid)
#   V_core   = Σ_j [Morse(r_j) + k_e * Q_j * erfc(r_j/σ) / r_j]  (short-range, per-atom)
#
# This is an EXACT decomposition — no approximation needed.
# σ controls the split length scale (typically ~2-4 Å).
#
# Reference: cpp/common_resources/cl/Surface.cl::eval_hybrid_grid()
# Parity: OpenCL kernel eval_hardcore_grid, eval_hybrid_grid

def eval_morse_coulomb_ref(atoms_pos, atoms_q, atoms_REQ, eval_points, nPBC=(4,4,0), lvec=None, R_damp=0.1):
    """Compute reference Morse+Coulomb potential at eval_points.
    
    Args:
        atoms_pos: (N,3) substrate atom positions
        atoms_q: (N,) charges
        atoms_REQ: (N,4) Morse params (R_eq, D_e, q, 0) per atom
        eval_points: (M,3) evaluation positions
        nPBC: (nx, ny, nz) PBC replica counts
        lvec: (3,3) lattice vectors
        R_damp: damping radius for Coulomb stability
    
    Returns:
        V_ref: (M,) reference potential in eV
        V_morse: (M,) Morse part
        V_coul: (M,) Coulomb part
    """
    COULOMB_CONST = 14.3996448915
    atoms_pos = np.asarray(atoms_pos, dtype=np.float64)
    atoms_q = np.asarray(atoms_q, dtype=np.float64)
    atoms_REQ = np.asarray(atoms_REQ, dtype=np.float64)
    eval_points = np.asarray(eval_points, dtype=np.float64)
    M = eval_points.shape[0]
    N = atoms_pos.shape[0]
    
    V_morse = np.zeros(M, dtype=np.float64)
    V_coul = np.zeros(M, dtype=np.float64)
    
    if lvec is None:
        lvec = np.eye(3, dtype=np.float64) * 20.0
    lvec = np.asarray(lvec, dtype=np.float64)
    
    npx, npy, npz = nPBC
    for ix in range(-npx, npx+1):
        for iy in range(-npy, npy+1):
            for iz in range(-npz, npz+1):
                shift = ix*lvec[0] + iy*lvec[1] + iz*lvec[2]
                for i in range(N):
                    dp = eval_points - (atoms_pos[i] + shift)
                    r = np.sqrt(np.sum(dp**2, axis=1))
                    r_safe = np.maximum(r, 1e-10)
                    # Morse: D_e * (exp(-2*alpha*(r-R_eq)) - 2*exp(-alpha*(r-R_eq)))
                    R_eq = atoms_REQ[i, 0]
                    D_e = atoms_REQ[i, 1]
                    alpha = atoms_REQ[i, 3] if atoms_REQ.shape[1] > 3 else 1.5
                    if D_e != 0.0:
                        e_arg = alpha * (r_safe - R_eq)
                        V_morse += D_e * (np.exp(-2*e_arg) - 2*np.exp(-e_arg))
                    # Damped Coulomb
                    q_i = atoms_q[i]
                    if q_i != 0.0:
                        V_coul += COULOMB_CONST * q_i / np.sqrt(r**2 + R_damp**2)
    
    V_ref = V_morse + V_coul
    return V_ref, V_morse, V_coul


def eval_hybrid_split(atoms_pos, atoms_q, atoms_REQ, eval_points, sigma=3.0, nPBC=(4,4,0), lvec=None, R_damp=0.1):
    """Compute erfc/erf split of Morse+Coulomb potential.
    
    V_smooth = Σ_j k_e * Q_j * erf(r_j/σ) / r_j   → B-spline grid
    V_core   = Σ_j [Morse(r_j) + k_e * Q_j * erfc(r_j/σ) / r_j]  → per-atom hardcore
    V_total  = V_smooth + V_core  (exact, no approximation)
    
    Args:
        sigma: split length scale (Å). Controls how much Coulomb goes to smooth vs core.
        R_damp: damping for r→0 stability in the core part
    
    Returns:
        V_total, V_smooth, V_core, V_morse, V_coul_erf, V_coul_erfc
    """
    from scipy.special import erf, erfc
    COULOMB_CONST = 14.3996448915
    atoms_pos = np.asarray(atoms_pos, dtype=np.float64)
    atoms_q = np.asarray(atoms_q, dtype=np.float64)
    atoms_REQ = np.asarray(atoms_REQ, dtype=np.float64)
    eval_points = np.asarray(eval_points, dtype=np.float64)
    M = eval_points.shape[0]
    N = atoms_pos.shape[0]
    
    V_morse = np.zeros(M, dtype=np.float64)
    V_coul_erf = np.zeros(M, dtype=np.float64)   # smooth part
    V_coul_erfc = np.zeros(M, dtype=np.float64)  # core part
    
    if lvec is None:
        lvec = np.eye(3, dtype=np.float64) * 20.0
    lvec = np.asarray(lvec, dtype=np.float64)
    
    npx, npy, npz = nPBC
    for ix in range(-npx, npx+1):
        for iy in range(-npy, npy+1):
            for iz in range(-npz, npz+1):
                shift = ix*lvec[0] + iy*lvec[1] + iz*lvec[2]
                for i in range(N):
                    dp = eval_points - (atoms_pos[i] + shift)
                    r = np.sqrt(np.sum(dp**2, axis=1))
                    r_safe = np.maximum(r, 1e-10)
                    # Morse
                    R_eq = atoms_REQ[i, 0]
                    D_e = atoms_REQ[i, 1]
                    alpha = atoms_REQ[i, 3] if atoms_REQ.shape[1] > 3 else 1.5
                    if D_e != 0.0:
                        e_arg = alpha * (r_safe - R_eq)
                        V_morse += D_e * (np.exp(-2*e_arg) - 2*np.exp(-e_arg))
                    # Coulomb erfc/erf split
                    q_i = atoms_q[i]
                    if q_i != 0.0:
                        r_damped = np.sqrt(r**2 + R_damp**2)
                        # erf part (smooth, long-range)
                        V_coul_erf += COULOMB_CONST * q_i * erf(r_safe / sigma) / r_safe
                        # erfc part (short-range, compact)
                        V_coul_erfc += COULOMB_CONST * q_i * erfc(r_safe / sigma) / r_damped
    
    V_smooth = V_coul_erf  # smooth part goes on B-spline grid
    V_core = V_morse + V_coul_erfc  # short-range part evaluated per-atom
    V_total = V_smooth + V_core
    
    return V_total, V_smooth, V_core, V_morse, V_coul_erf, V_coul_erfc


def fit_bspline_3d(V_ref_3d, ns_grid, lvec, z_range=None, n_smooth=1, step=None):
    """Fit tricubic B-spline to 3D reference potential.
    
    Uses scipy B-spline fitting. Returns fitted grid and residual.
    
    Args:
        V_ref_3d: (nx, ny, nz) reference potential on a regular grid
        ns_grid: (nx, ny, nz) grid shape
        lvec: (3,3) lattice vectors
        z_range: (z_min, z_max) for fitting range (None = all)
        n_smooth: smoothing factor
        step: B-spline node spacing in Å (None = use ns//16 heuristic)
    
    Returns:
        V_bspline: (nx, ny, nz) B-spline fitted potential
        residual: (nx, ny, nz) V_ref - V_bspline
    """
    from scipy.ndimage import gaussian_filter
    
    if step is not None:
        nx, ny, nz = ns_grid
        dx = np.linalg.norm(lvec[0]) / nx
        dy = np.linalg.norm(lvec[1]) / ny
        dz = np.linalg.norm(lvec[2]) / nz
        sigma = tuple(max(1, int(round(s / d))) for s, d in zip([step, step, step], [dx, dy, dz]))
    else:
        sigma = tuple(max(1, s // 16) for s in ns_grid)
    
    V_bspline = gaussian_filter(V_ref_3d.astype(np.float64), sigma=sigma, mode=('wrap', 'wrap', 'nearest'))
    
    residual = V_ref_3d - V_bspline
    return V_bspline, residual


def fit_hardcore_residual(residual_3d, grid_coords, hc_atoms_pos, hc_types, R_cuts, n_basis=4):
    """Fit per-atom radial functions to the residual potential.
    
    Uses compact polynomial basis: (1 - r/R_cut)^(2n) for n=0..n_basis-1.
    Least-squares fit per atom type.
    
    Args:
        residual_3d: (nx, ny, nz) residual potential
        grid_coords: (nx, ny, nz, 3) grid point coordinates
        hc_atoms_pos: (N_hc, 3) hardcore atom positions (top layer)
        hc_types: (N_hc,) atom type indices
        R_cuts: (n_types,) cutoff radii per type
        n_basis: number of basis functions per atom
    
    Returns:
        coeffs: (n_types, n_basis) fitted coefficients
        V_hardcore: (nx, ny, nz) reconstructed hardcore potential
    """
    n_types = len(R_cuts)
    HC_MAX_BASIS = 8
    n_basis = min(n_basis, HC_MAX_BASIS)
    
    nx, ny, nz = residual_3d.shape
    M = nx * ny * nz  # total grid points
    N_hc = hc_atoms_pos.shape[0]
    
    # Flatten grid
    resid_flat = residual_3d.ravel()
    coords_flat = grid_coords.reshape(-1, 3)
    
    # Build design matrix: for each grid point, sum over hardcore atoms
    # A[i, n*ntype + t] = sum over atoms of type t contributing basis_n
    # This is a per-type fitting: group all atoms of same type
    
    # Simpler approach: fit per-atom independently, then sum
    # But per-atom fitting is underdetermined if atoms overlap.
    # Instead: build a combined design matrix with per-type coefficients.
    
    n_cols = n_types * n_basis
    A = np.zeros((M, n_cols), dtype=np.float64)
    
    for ia in range(N_hc):
        ityp = hc_types[ia]
        R_cut = R_cuts[ityp]
        dp = coords_flat - hc_atoms_pos[ia]
        r = np.sqrt(np.sum(dp**2, axis=1))
        mask = r < R_cut
        x = np.zeros(M, dtype=np.float64)
        x[mask] = 1.0 - r[mask] / R_cut
        x2 = x * x
        xp = np.ones(M, dtype=np.float64)
        for n in range(n_basis):
            col_idx = ityp * n_basis + n
            A[:, col_idx] += xp
            xp *= x2
    
    # Regularized least squares
    lam = 1e-6 * np.max(np.abs(A))**2  # small regularization
    AtA = A.T @ A + lam * np.eye(n_cols)
    Atb = A.T @ resid_flat
    coeffs_flat = np.linalg.solve(AtA, Atb)
    
    coeffs = coeffs_flat.reshape(n_types, n_basis)
    
    # Reconstruct hardcore potential
    V_hardcore = (A @ coeffs_flat).reshape(nx, ny, nz)
    
    return coeffs, V_hardcore


def eval_folded_morse(atoms_pos, atoms_REQ, eval_points, atom_mask=None,
                      R_fold=4.0, lvec=None, nPBC=(1,1,0), n_win=2):
    """Evaluate folded exponential (Morse) basis with compact-support window.
    
    V_folded = Σ_{j in mask} Morse_j(r_j) * w(r_j)  where w(r) = (1-(r/R_fold))^(2n) for r < R_fold, else 0
    The window ensures V_folded = 0 and C^(2n-1) continuous at r = R_fold (no discontinuity).
    
    Args:
        atoms_pos: (N,3) all atom positions
        atoms_REQ: (N,4) Morse params (R_eq, D_e, q, alpha)
        eval_points: (M,3) evaluation positions
        atom_mask: (N,) bool mask of atoms to include (e.g. top layer). None = all.
        R_fold: cutoff radius for folded basis
        lvec: (3,3) lattice vectors for PBC
        nPBC: (nx, ny, nz) periodic replicas
        n_win: window exponent parameter (w = (1-r/R_fold)^(2*n_win), default n_win=2 → C3 at R_fold)
    
    Returns:
        V_folded: (M,) folded Morse potential
    """
    atoms_pos = np.asarray(atoms_pos, dtype=np.float64)
    atoms_REQ = np.asarray(atoms_REQ, dtype=np.float64)
    eval_points = np.asarray(eval_points, dtype=np.float64)
    M = eval_points.shape[0]
    N = atoms_pos.shape[0]
    
    if atom_mask is None:
        atom_mask = np.ones(N, dtype=bool)
    
    if lvec is None:
        lvec = np.eye(3, dtype=np.float64) * 20.0
    lvec = np.asarray(lvec, dtype=np.float64)
    npx, npy, npz = nPBC
    
    V_folded = np.zeros(M, dtype=np.float64)
    
    for j in range(N):
        if not atom_mask[j]:
            continue
        R_eq = atoms_REQ[j, 0]
        D_e = atoms_REQ[j, 1]
        alpha = atoms_REQ[j, 3] if atoms_REQ.shape[1] > 3 else 1.5
        
        for ix in range(-npx, npx+1):
            for iy in range(-npy, npy+1):
                for iz in range(-npz, npz+1):
                    shift = ix*lvec[0] + iy*lvec[1] + iz*lvec[2]
                    atom_pos = atoms_pos[j] + shift
                    dp = eval_points - atom_pos
                    r = np.sqrt(np.sum(dp**2, axis=1))
                    r_safe = np.maximum(r, 1e-10)
                    
                    mask_fold = r < R_fold
                    e_arg = alpha * (r_safe - R_eq)
                    V_morse = D_e * (np.exp(-2*e_arg) - 2*np.exp(-e_arg)) if D_e != 0 else np.zeros_like(r_safe)
                    # Compact-support window: w(r) = (1-(r/R_fold))^(2n) for r < R_fold, else 0
                    w = np.where(mask_fold, (1.0 - r / R_fold)**(2*n_win), 0.0)
                    V_folded += V_morse * w
    
    return V_folded


def eval_hardcore_quadratic_v2(atoms_pos, atoms_REQ, eval_points, atom_mask=None,
                               R_c=2.0, R_fold=4.0, lvec=None, nPBC=(1,1,0), n_win=2):
    """Evaluate hardcore potential with compact-support window.
    
    For r < R_c: quadratic p(r) = a + b*r + c*r² matching V_w(R_c), V_w'(R_c), V_w''(R_c)
                 where V_w(r) = Morse(r) * w(r) is the windowed Morse
    For R_c ≤ r < R_fold: Morse(r) * w(r) where w(r) = (1-(r/R_fold))^(2n)
    For r ≥ R_fold: 0 (w=0, no discontinuity)
    
    C2 continuous at R_c (quadratic↔windowed-Morse switch, matches value+1st+2nd deriv).
    C(2n-1) continuous at R_fold (window goes to 0 with 2n-1 zero derivatives).
    
    Args:
        atoms_pos: (N,3) all atom positions
        atoms_REQ: (N,4) Morse params
        eval_points: (M,3) evaluation positions
        atom_mask: (N,) bool mask of atoms to include. None = all.
        R_c: quadratic↔Morse switch radius
        R_fold: cutoff radius (window goes to 0)
        lvec: (3,3) lattice vectors for PBC
        nPBC: (nx, ny, nz) periodic replicas
        n_win: window exponent (w = (1-r/R_fold)^(2*n_win), default 2)
    
    Returns:
        V_hardcore: (M,) hardcore potential
        match_params: list of dicts with per-atom a, b, c, R_c
    """
    atoms_pos = np.asarray(atoms_pos, dtype=np.float64)
    atoms_REQ = np.asarray(atoms_REQ, dtype=np.float64)
    eval_points = np.asarray(eval_points, dtype=np.float64)
    M = eval_points.shape[0]
    N = atoms_pos.shape[0]
    
    if atom_mask is None:
        atom_mask = np.ones(N, dtype=bool)
    
    if lvec is None:
        lvec = np.eye(3, dtype=np.float64) * 20.0
    lvec = np.asarray(lvec, dtype=np.float64)
    npx, npy, npz = nPBC
    
    V_hardcore = np.zeros(M, dtype=np.float64)
    match_params = []
    
    for j in range(N):
        if not atom_mask[j]:
            match_params.append(None)
            continue
        R_eq = atoms_REQ[j, 0]
        D_e = atoms_REQ[j, 1]
        alpha = atoms_REQ[j, 3] if atoms_REQ.shape[1] > 3 else 1.5
        
        # --- Raw Morse and its derivatives at R_c ---
        e_c = alpha * (R_c - R_eq)
        exp1 = np.exp(-e_c); exp2 = np.exp(-2*e_c)
        M_rc  = D_e * (exp2 - 2*exp1) if D_e != 0 else 0.0      # Morse(R_c)
        dM_rc = 2*alpha*D_e*(exp1 - exp2) if D_e != 0 else 0.0   # Morse'(R_c)
        ddM_rc = 2*alpha**2*D_e*(2*exp2 - exp1) if D_e != 0 else 0.0  # Morse''(R_c)
        
        # --- Window w(r) = (1 - r/R_fold)^(2n) and its derivatives at R_c ---
        u_c = 1.0 - R_c / R_fold
        w_rc   = u_c**(2*n_win)
        dw_rc  = -2*n_win / R_fold * u_c**(2*n_win - 1)
        ddw_rc = 2*n_win*(2*n_win - 1) / R_fold**2 * u_c**(2*n_win - 2)
        
        # --- Windowed Morse V_w = Morse * w and its derivatives at R_c (product rule) ---
        E_rc  = M_rc * w_rc
        dE_rc = dM_rc * w_rc + M_rc * dw_rc
        ddE_rc = ddM_rc * w_rc + 2*dM_rc * dw_rc + M_rc * ddw_rc
        
        # --- Fit quadratic p(r) = a + b*r + c*r² to windowed Morse at R_c ---
        c_coef = ddE_rc / 2.0
        b_coef = dE_rc - ddE_rc * R_c
        a_coef = E_rc - dE_rc * R_c + ddE_rc * R_c**2 / 2.0
        
        match_params.append({'R_c': R_c, 'R_fold': R_fold, 'a': a_coef, 'b': b_coef, 'c': c_coef,
                             'E_rc': E_rc, 'dE_rc': dE_rc, 'ddE_rc': ddE_rc,
                             'M_rc': M_rc, 'w_rc': w_rc})
        
        for ix in range(-npx, npx+1):
            for iy in range(-npy, npy+1):
                for iz in range(-npz, npz+1):
                    shift = ix*lvec[0] + iy*lvec[1] + iz*lvec[2]
                    atom_pos = atoms_pos[j] + shift
                    dp = eval_points - atom_pos
                    r = np.sqrt(np.sum(dp**2, axis=1))
                    r_safe = np.maximum(r, 1e-10)
                    
                    e_arg = alpha * (r_safe - R_eq)
                    V_morse = D_e * (np.exp(-2*e_arg) - 2*np.exp(-e_arg)) if D_e != 0 else np.zeros_like(r_safe)
                    
                    # quadratic for r < R_c, Morse*window for R_c ≤ r < R_fold, 0 for r ≥ R_fold
                    mask_quad = r < R_c
                    mask_fold = r < R_fold
                    w = np.where(mask_fold, (1.0 - r / R_fold)**(2*n_win), 0.0)
                    V_hc_atom = np.where(mask_quad, a_coef + b_coef * r + c_coef * r**2, V_morse * w)
                    V_hardcore += V_hc_atom
    
    return V_hardcore, match_params


def eval_hardcore_potential(eval_points, hc_atoms_pos, hc_types, coeffs, R_cuts, n_basis=4):
    """Evaluate hardcore potential at arbitrary points (Python reference for parity).
    
    Args:
        eval_points: (M, 3) evaluation positions
        hc_atoms_pos: (N_hc, 3) hardcore atom positions
        hc_types: (N_hc,) type indices
        coeffs: (n_types, n_basis) fitted coefficients
        R_cuts: (n_types,) cutoff radii
        n_basis: number of basis functions
    
    Returns:
        V_hardcore: (M,) hardcore potential
    """
    eval_points = np.asarray(eval_points, dtype=np.float64)
    hc_atoms_pos = np.asarray(hc_atoms_pos, dtype=np.float64)
    M = eval_points.shape[0]
    N_hc = hc_atoms_pos.shape[0]
    V = np.zeros(M, dtype=np.float64)
    
    for ia in range(N_hc):
        ityp = hc_types[ia]
        R_cut = R_cuts[ityp]
        dp = eval_points - hc_atoms_pos[ia]
        r = np.sqrt(np.sum(dp**2, axis=1))
        mask = r < R_cut
        x = np.zeros(M, dtype=np.float64)
        x[mask] = 1.0 - r[mask] / R_cut
        x2 = x * x
        xp = np.ones(M, dtype=np.float64)
        for n in range(n_basis):
            V += coeffs[ityp, n] * xp
            xp *= x2
    
    return V


def make_nacl_1x1_l3_system():
    """Create NaCl 1x1 3-layer test system.
    
    Returns:
        dict with apos, qs, enames, lvec, REQs, z_top, top_layer_mask
    """
    # From NaCl_1x1_L3.xyz
    apos = np.array([
        [0.0, 0.0, -3.25],   # Na
        [0.0, 0.0, -6.078427], # Cl
        [0.0, 0.0, -8.906854], # Na
        [2.0, 2.0, -3.25],   # Cl
        [2.0, 2.0, -6.078427], # Na
        [2.0, 2.0, -8.906854], # Cl
    ], dtype=np.float64)
    qs = np.array([+0.7, -0.7, +0.7, -0.7, +0.7, -0.7], dtype=np.float64)
    enames = ['Na', 'Cl', 'Na', 'Cl', 'Na', 'Cl']
    lvec = np.array([
        [4.0, 0.0, 0.0],
        [0.0, 4.0, 0.0],
        [0.0, 0.0, 20.0],
    ], dtype=np.float64)
    
    # Morse parameters from ElementTypes.dat
    # Na: RvdW=1.4915, EvdW=0.00130, alpha~1.5
    # Cl: RvdW=1.9735, EvdW=0.00984, alpha~1.5
    REQs = np.zeros((6, 4), dtype=np.float64)
    for i, e in enumerate(enames):
        if e == 'Na':
            REQs[i] = [1.4915, 0.00130, 0.7, 1.5]
        elif e == 'Cl':
            REQs[i] = [1.9735, 0.00984, -0.7, 1.5]
    
    z_top = float(np.max(apos[:, 2]))  # -3.25
    # Top layer: atoms at z = z_top
    top_layer_mask = np.abs(apos[:, 2] - z_top) < 0.01
    
    return {
        'apos': apos, 'qs': qs, 'enames': enames, 'lvec': lvec,
        'REQs': REQs, 'z_top': z_top, 'top_layer_mask': top_layer_mask,
    }


def test_hybrid_potential_nacl(save_dir='results_hybrid', Ng=80, n_basis=4, R_cut_hc=5.0, nPBC_ref=(6,6,0), sigma_split=3.0, R_match=None,
                               R_fold=4.0, R_c_quad=2.0, spline_step=1.0):
    """Test hybrid potential on NaCl 1x1 L3 — two-step approach.
    
    Step 1: Fit V_total (Morse+Coulomb, all atoms, 6×6 supercell) using
            V_model1 = B-spline(1.0 Å step) + folded Morse (top-layer atoms, cutoff R_fold)
            Error = V_total - V_model1  →  shows quality of B-spline + folded atomic basis
    
    Step 2: Replace folded Morse with quadratic derivative matching for r < R_c:
            V_hardcore = p(r) for r < R_c, Morse for R_c ≤ r < R_fold, 0 for r ≥ R_fold
            V_model2 = B-spline(1.0 Å step) + V_hardcore
            Error = V_total - V_model2  →  shows effect of quadratic replacement
    
    Args:
        save_dir: output directory
        Ng: grid resolution for 1D scans
        R_fold: cutoff for folded Morse basis (default 4.0 Å)
        R_c_quad: quadratic↔Morse switch radius (default 2.0 Å)
        spline_step: B-spline node spacing in Å (default 1.0)
        nPBC_ref: PBC replicas for reference V_total
    """
    os.makedirs(save_dir, exist_ok=True)
    
    print("="*60)
    print("Hybrid Potential Test: NaCl 1x1 L3 (two-step: B-spline+folded Morse, then quadratic)")
    print("="*60)
    
    # --- Setup system ---
    sys = make_nacl_1x1_l3_system()
    apos = sys['apos']; qs = sys['qs']; enames = sys['enames']
    lvec = sys['lvec']; REQs = sys['REQs']; z_top = sys['z_top']
    top_mask = sys['top_layer_mask']
    
    print(f"System: {len(apos)} atoms, z_top={z_top:.3f}")
    print(f"Top layer atoms: {np.sum(top_mask)} ({[enames[i] for i in range(len(enames)) if top_mask[i]]})")
    print(f"Lattice: {lvec[0,:2]}, {lvec[1,:2]}")
    print(f"Parameters: R_fold={R_fold} Å, R_c={R_c_quad} Å, spline_step={spline_step} Å")
    
    # --- Build evaluation grid ---
    x_min, x_max = 0.0, 4.0
    z_min = z_top + 0.5
    z_max = z_top + 6.0
    y_fixed = 0.0
    
    xv = np.linspace(x_min, x_max, Ng)
    zv = np.linspace(z_min, z_max, Ng)
    X, Z = np.meshgrid(xv, zv, indexing='ij')
    Y = np.full_like(X, y_fixed)
    
    # 3D grid for B-spline fitting
    print(f"\nBuilding 3D grid for B-spline fitting...")
    nx_3d = 48; ny_3d = 48; nz_3d = 64
    x_3d = np.linspace(0, 4.0, nx_3d)
    y_3d = np.linspace(0, 4.0, ny_3d)
    z_3d = np.linspace(z_min, z_max, nz_3d)
    X3d, Y3d, Z3d = np.meshgrid(x_3d, y_3d, z_3d, indexing='ij')
    coords_3d = np.stack([X3d, Y3d, Z3d], axis=-1)
    eval_3d = coords_3d.reshape(-1, 3)
    
    # --- Compute V_total reference (all atoms, large supercell) ---
    print(f"\nComputing V_total reference (Morse+Coulomb, {nPBC_ref} PBC)...")
    V_total_flat, _, _, _, _, _ = eval_hybrid_split(
        apos, qs, REQs, eval_3d, sigma=1e6, nPBC=nPBC_ref, lvec=lvec, R_damp=0.1
    )
    V_total_3d = V_total_flat.reshape(nx_3d, ny_3d, nz_3d)
    print(f"  V_total range: [{V_total_3d.min():.4f}, {V_total_3d.max():.4f}] eV")
    
    # === STEP 1: B-spline + folded Morse ===
    print(f"\n{'='*40}")
    print(f"Step 1: B-spline({spline_step} Å) + folded Morse (top layer, R_fold={R_fold} Å)")
    print(f"{'='*40}")
    
    V_folded_flat = eval_folded_morse(apos, REQs, eval_3d, atom_mask=top_mask,
                                      R_fold=R_fold, lvec=lvec, nPBC=(1,1,0))
    V_folded_3d = V_folded_flat.reshape(nx_3d, ny_3d, nz_3d)
    
    # Residual for B-spline: V_total - V_folded
    V_residual_3d = V_total_3d - V_folded_3d
    print(f"  V_folded range:  [{V_folded_3d.min():.6f}, {V_folded_3d.max():.6f}] eV")
    print(f"  V_residual range: [{V_residual_3d.min():.4f}, {V_residual_3d.max():.4f}] eV")
    
    V_bspline_3d, bspline_err_3d = fit_bspline_3d(V_residual_3d, (nx_3d, ny_3d, nz_3d), lvec, step=spline_step)
    print(f"  V_bspline range: [{V_bspline_3d.min():.4f}, {V_bspline_3d.max():.4f}] eV")
    
    V_model1_3d = V_bspline_3d + V_folded_3d
    err1_3d = V_total_3d - V_model1_3d
    print(f"  Step 1 error RMSE: {np.sqrt(np.mean(err1_3d**2)):.6f} eV")
    print(f"  Step 1 error max:  {np.max(np.abs(err1_3d)):.6f} eV")
    
    # === STEP 2: Replace folded Morse with quadratic derivative matching ===
    print(f"\n{'='*40}")
    print(f"Step 2: Quadratic p(r) for r<{R_c_quad} Å, Morse for {R_c_quad}≤r<{R_fold} Å")
    print(f"{'='*40}")
    
    V_hardcore_flat, match_params = eval_hardcore_quadratic_v2(
        apos, REQs, eval_3d, atom_mask=top_mask,
        R_c=R_c_quad, R_fold=R_fold, lvec=lvec, nPBC=(1,1,0)
    )
    V_hardcore_3d = V_hardcore_flat.reshape(nx_3d, ny_3d, nz_3d)
    
    # Residual for B-spline: V_total - V_hardcore
    V_residual2_3d = V_total_3d - V_hardcore_3d
    print(f"  V_hardcore range:  [{V_hardcore_3d.min():.6f}, {V_hardcore_3d.max():.6f}] eV")
    print(f"  V_residual2 range: [{V_residual2_3d.min():.4f}, {V_residual2_3d.max():.4f}] eV")
    
    V_bspline2_3d, bspline_err2_3d = fit_bspline_3d(V_residual2_3d, (nx_3d, ny_3d, nz_3d), lvec, step=spline_step)
    print(f"  V_bspline2 range: [{V_bspline2_3d.min():.4f}, {V_bspline2_3d.max():.4f}] eV")
    
    V_model2_3d = V_bspline2_3d + V_hardcore_3d
    err2_3d = V_total_3d - V_model2_3d
    print(f"  Step 2 error RMSE: {np.sqrt(np.mean(err2_3d**2)):.6f} eV")
    print(f"  Step 2 error max:  {np.max(np.abs(err2_3d)):.6f} eV")
    
    for j, mp in enumerate(match_params):
        if mp is None: continue
        print(f"  Atom {j} ({enames[j]}): R_c={mp['R_c']:.3f}, a={mp['a']:.6f}, b={mp['b']:.6f}, c={mp['c']:.6f}, Morse(Rc)={mp['E_rc']:.6f}")
    
    # --- min_dist for masking ---
    coords_flat = coords_3d.reshape(-1, 3)
    min_dist = np.full(len(coords_flat), 1e10)
    npx2, npy2, npz2 = 2, 2, 0
    for j in range(len(apos)):
        for ix in range(-npx2, npx2+1):
            for iy in range(-npy2, npy2+1):
                for iz in range(-npz2, npz2+1):
                    shift = ix*lvec[0] + iy*lvec[1] + iz*lvec[2]
                    d = np.sqrt(np.sum((coords_flat - (apos[j] + shift))**2, axis=1))
                    min_dist = np.minimum(min_dist, d)
    min_dist_3d = min_dist.reshape(nx_3d, ny_3d, nz_3d)
    mask_phys = min_dist_3d >= R_c_quad
    
    # --- Extract XZ slices (at y=0) ---
    iy = np.argmin(np.abs(y_3d - y_fixed))
    
    V_total_xz = V_total_3d[:, iy, :]
    V_folded_xz = V_folded_3d[:, iy, :]
    V_bspline_xz = V_bspline_3d[:, iy, :]
    V_model1_xz = V_model1_3d[:, iy, :]
    err1_xz = err1_3d[:, iy, :]
    
    V_hardcore_xz = V_hardcore_3d[:, iy, :]
    V_bspline2_xz = V_bspline2_3d[:, iy, :]
    V_model2_xz = V_model2_3d[:, iy, :]
    err2_xz = err2_3d[:, iy, :]
    min_dist_xz = min_dist_3d[:, iy, :]
    
    # --- Plotting ---
    print(f"\nGenerating plots...")
    
    def overlay_atoms(ax, ms=8):
        for i in range(len(apos)):
            if abs(apos[i, 1] - y_fixed) < 0.5:
                color = 'purple' if enames[i] == 'Na' else 'green'
                ax.plot(apos[i, 0], apos[i, 2], 'o', color=color, ms=ms, alpha=0.8, zorder=5)
    
    def draw_Rc_circles(ax, R_val, alpha=0.5):
        for i in range(len(apos)):
            if abs(apos[i, 1] - y_fixed) < 0.5 and top_mask[i]:
                c = 'purple' if enames[i] == 'Na' else 'green'
                circle = plt.Circle((apos[i, 0], apos[i, 2]), R_val, fill=False, linestyle='--', color=c, alpha=alpha, lw=1.5)
                ax.add_patch(circle)
    
    # Spline node spacing
    dx = x_3d[1] - x_3d[0]; dz = z_3d[1] - z_3d[0]
    sigma_x = max(1, int(round(spline_step / dx))); sigma_z = max(1, int(round(spline_step / dz)))
    node_dx = sigma_x * dx; node_dz = sigma_z * dz
    n_nodes_x = int((x_3d[-1] - x_3d[0]) / node_dx) + 1
    n_nodes_z = int((z_3d[-1] - z_3d[0]) / node_dz) + 1
    print(f"  Spline nodes: {n_nodes_x} x {n_nodes_z} (spacing {node_dx:.2f} x {node_dz:.2f} Å)")
    
    def draw_spline_nodes(ax, alpha=0.15):
        xs = np.arange(x_3d[0], x_3d[-1] + node_dx*0.5, node_dx)
        zs = np.arange(z_3d[0], z_3d[-1] + node_dz*0.5, node_dz)
        for xi in xs: ax.axvline(xi, color='gray', alpha=alpha, lw=0.3)
        for zi in zs: ax.axhline(zi, color='gray', alpha=alpha, lw=0.3)
        ZX, XX = np.meshgrid(zs, xs)
        ax.plot(XX.ravel(), ZX.ravel(), '.', color='gray', alpha=alpha*2, ms=2, zorder=1)
    
    # ========== Figure 1: Step 1 — B-spline + folded Morse ==========
    fig1, axes1 = plt.subplots(2, 3, figsize=(18, 12))
    fig1.suptitle(f'Step 1: B-spline({spline_step} Å) + Folded Morse (top L1, R_fold={R_fold} Å)', fontsize=14)
    
    vmax = max(np.abs(V_total_xz).max(), 1e-10)
    
    # (1) Reference
    ax = axes1[0, 0]
    im = ax.pcolormesh(x_3d, z_3d, V_total_xz.T, shading='auto', cmap='RdBu_r', vmin=-vmax, vmax=vmax)
    ax.set_title('(1) Reference: V_total (Morse+Coulomb, 6×6 PBC)', fontsize=10)
    ax.set_xlabel('x (Å)'); ax.set_ylabel('z (Å)'); ax.set_aspect('equal')
    plt.colorbar(im, ax=ax, label='eV', fraction=0.046); overlay_atoms(ax)
    
    # (2) Model1 = B-spline + folded Morse
    ax = axes1[0, 1]
    im = ax.pcolormesh(x_3d, z_3d, V_model1_xz.T, shading='auto', cmap='RdBu_r', vmin=-vmax, vmax=vmax)
    ax.set_title('(2) Model1: B-spline + folded Morse', fontsize=10)
    ax.set_xlabel('x (Å)'); ax.set_ylabel('z (Å)'); ax.set_aspect('equal')
    plt.colorbar(im, ax=ax, label='eV', fraction=0.046); overlay_atoms(ax); draw_spline_nodes(ax, alpha=0.08)
    
    # (3) Error1
    ax = axes1[0, 2]
    vmax_e1 = max(np.abs(err1_xz).max(), 1e-10)
    im = ax.pcolormesh(x_3d, z_3d, err1_xz.T, shading='auto', cmap='RdBu_r', vmin=-vmax_e1, vmax=vmax_e1)
    ax.set_title(f'(3) Error1: V_total - V_model1 (RMSE={np.sqrt(np.mean(err1_3d**2)):.4f} eV)', fontsize=10)
    ax.set_xlabel('x (Å)'); ax.set_ylabel('z (Å)'); ax.set_aspect('equal')
    plt.colorbar(im, ax=ax, label='eV', fraction=0.046); overlay_atoms(ax, ms=6); draw_Rc_circles(ax, R_fold)
    
    # (4) B-spline component
    ax = axes1[1, 0]
    vmax_bs = max(np.abs(V_bspline_xz).max(), 1e-10)
    im = ax.pcolormesh(x_3d, z_3d, V_bspline_xz.T, shading='auto', cmap='RdBu_r', vmin=-vmax_bs, vmax=vmax_bs)
    ax.set_title(f'(4) B-spline({n_nodes_x}×{n_nodes_z} nodes, σ={node_dx:.1f}×{node_dz:.1f}Å)', fontsize=10)
    ax.set_xlabel('x (Å)'); ax.set_ylabel('z (Å)'); ax.set_aspect('equal')
    plt.colorbar(im, ax=ax, label='eV', fraction=0.046); overlay_atoms(ax, ms=6); draw_spline_nodes(ax)
    
    # (5) Folded Morse component
    ax = axes1[1, 1]
    vmax_f = max(np.abs(V_folded_xz).max(), 1e-10)
    im = ax.pcolormesh(x_3d, z_3d, V_folded_xz.T, shading='auto', cmap='RdBu_r', vmin=-vmax_f, vmax=vmax_f)
    ax.set_title(f'(5) Folded Morse (top L1, {int(np.sum(top_mask))} atoms, R_fold={R_fold} Å)', fontsize=10)
    ax.set_xlabel('x (Å)'); ax.set_ylabel('z (Å)'); ax.set_aspect('equal')
    plt.colorbar(im, ax=ax, label='eV', fraction=0.046); overlay_atoms(ax, ms=6); draw_Rc_circles(ax, R_fold)
    
    # (6) 1D scans
    ax = axes1[1, 2]
    ix_na = np.argmin(np.abs(x_3d - 0.0))
    ax.plot(z_3d, V_total_3d[ix_na, iy, :], 'k-', lw=2, label='V_total (ref)')
    ax.plot(z_3d, V_model1_3d[ix_na, iy, :], 'r--', lw=1.5, label='V_model1')
    ax.plot(z_3d, V_bspline_3d[ix_na, iy, :], 'b:', lw=1, alpha=0.7, label='B-spline')
    ax.plot(z_3d, V_folded_3d[ix_na, iy, :], 'g-.', lw=1, alpha=0.7, label='folded Morse')
    ax2 = ax.twinx(); ax2.plot(z_3d, err1_3d[ix_na, iy, :], 'm-', lw=1, alpha=0.5)
    ax2.set_ylabel('Error (eV)', color='magenta')
    ax.set_xlabel('z (Å)'); ax.set_ylabel('V (eV)')
    ax.set_title('1D scan above Na (x=0, y=0)', fontsize=10)
    ax.legend(fontsize=7, loc='upper left'); ax.grid(True, alpha=0.3); ax.axhline(0, color='gray', lw=0.5)
    
    fig1.tight_layout()
    fig1_path = os.path.join(save_dir, 'step1_bspline_folded_morse.png')
    fig1.savefig(fig1_path, dpi=150); plt.close(fig1)
    print(f"  Saved {fig1_path}")
    
    # ========== Figure 2: Step 2 — Quadratic replacement ==========
    fig2, axes2 = plt.subplots(2, 3, figsize=(18, 12))
    fig2.suptitle(f'Step 2: Quadratic p(r) for r<{R_c_quad} Å, Morse for {R_c_quad}≤r<{R_fold} Å (same B-spline {spline_step} Å)', fontsize=14)
    
    # (1) Reference
    ax = axes2[0, 0]
    im = ax.pcolormesh(x_3d, z_3d, V_total_xz.T, shading='auto', cmap='RdBu_r', vmin=-vmax, vmax=vmax)
    ax.set_title('(1) Reference: V_total', fontsize=10)
    ax.set_xlabel('x (Å)'); ax.set_ylabel('z (Å)'); ax.set_aspect('equal')
    plt.colorbar(im, ax=ax, label='eV', fraction=0.046); overlay_atoms(ax)
    
    # (2) Model2 = B-spline + hardcore(quadratic+Morse)
    ax = axes2[0, 1]
    im = ax.pcolormesh(x_3d, z_3d, V_model2_xz.T, shading='auto', cmap='RdBu_r', vmin=-vmax, vmax=vmax)
    ax.set_title('(2) Model2: B-spline + V_hardcore(quad+Morse)', fontsize=10)
    ax.set_xlabel('x (Å)'); ax.set_ylabel('z (Å)'); ax.set_aspect('equal')
    plt.colorbar(im, ax=ax, label='eV', fraction=0.046); overlay_atoms(ax); draw_spline_nodes(ax, alpha=0.08)
    
    # (3) Error2
    ax = axes2[0, 2]
    vmax_e2 = max(np.abs(err2_xz).max(), 1e-10)
    im = ax.pcolormesh(x_3d, z_3d, err2_xz.T, shading='auto', cmap='RdBu_r', vmin=-vmax_e2, vmax=vmax_e2)
    ax.set_title(f'(3) Error2: V_total - V_model2 (RMSE={np.sqrt(np.mean(err2_3d**2)):.4f} eV)', fontsize=10)
    ax.set_xlabel('x (Å)'); ax.set_ylabel('z (Å)'); ax.set_aspect('equal')
    plt.colorbar(im, ax=ax, label='eV', fraction=0.046); overlay_atoms(ax, ms=6)
    draw_Rc_circles(ax, R_c_quad, alpha=0.7); draw_Rc_circles(ax, R_fold, alpha=0.3)
    
    # (4) B-spline2 component
    ax = axes2[1, 0]
    vmax_bs2 = max(np.abs(V_bspline2_xz).max(), 1e-10)
    im = ax.pcolormesh(x_3d, z_3d, V_bspline2_xz.T, shading='auto', cmap='RdBu_r', vmin=-vmax_bs2, vmax=vmax_bs2)
    ax.set_title(f'(4) B-spline2 (residual after quadratic)', fontsize=10)
    ax.set_xlabel('x (Å)'); ax.set_ylabel('z (Å)'); ax.set_aspect('equal')
    plt.colorbar(im, ax=ax, label='eV', fraction=0.046); overlay_atoms(ax, ms=6); draw_spline_nodes(ax)
    
    # (5) HardCore component (quadratic + Morse)
    ax = axes2[1, 1]
    vmax_hc = max(np.abs(V_hardcore_xz).max(), 1e-10)
    im = ax.pcolormesh(x_3d, z_3d, V_hardcore_xz.T, shading='auto', cmap='RdBu_r', vmin=-vmax_hc, vmax=vmax_hc)
    ax.set_title(f'(5) V_hardcore: p(r) for r<{R_c_quad}, Morse for {R_c_quad}≤r<{R_fold}', fontsize=9)
    ax.set_xlabel('x (Å)'); ax.set_ylabel('z (Å)'); ax.set_aspect('equal')
    plt.colorbar(im, ax=ax, label='eV', fraction=0.046); overlay_atoms(ax, ms=6)
    draw_Rc_circles(ax, R_c_quad, alpha=0.7); draw_Rc_circles(ax, R_fold, alpha=0.3)
    
    # (6) 1D scans
    ax = axes2[1, 2]
    ax.plot(z_3d, V_total_3d[ix_na, iy, :], 'k-', lw=2, label='V_total (ref)')
    ax.plot(z_3d, V_model2_3d[ix_na, iy, :], 'r--', lw=1.5, label='V_model2')
    ax.plot(z_3d, V_bspline2_3d[ix_na, iy, :], 'b:', lw=1, alpha=0.7, label='B-spline2')
    ax.plot(z_3d, V_hardcore_3d[ix_na, iy, :], 'g-.', lw=1, alpha=0.7, label='V_hardcore')
    ax2 = ax.twinx(); ax2.plot(z_3d, err2_3d[ix_na, iy, :], 'm-', lw=1, alpha=0.5)
    ax2.set_ylabel('Error (eV)', color='magenta')
    ax.set_xlabel('z (Å)'); ax.set_ylabel('V (eV)')
    ax.set_title('1D scan above Na (x=0, y=0)', fontsize=10)
    ax.legend(fontsize=7, loc='upper left'); ax.grid(True, alpha=0.3); ax.axhline(0, color='gray', lw=0.5)
    
    fig2.tight_layout()
    fig2_path = os.path.join(save_dir, 'step2_quadratic_replacement.png')
    fig2.savefig(fig2_path, dpi=150); plt.close(fig2)
    print(f"  Saved {fig2_path}")
    
    # ========== Figure 3: Error comparison ==========
    fig3, axes3 = plt.subplots(1, 3, figsize=(18, 6))
    fig3.suptitle('Error Comparison: Step 1 (B-spline+folded Morse) vs Step 2 (quadratic replacement)', fontsize=13)
    
    # Error1 map
    ax = axes3[0]
    im = ax.pcolormesh(x_3d, z_3d, err1_xz.T, shading='auto', cmap='RdBu_r', vmin=-vmax_e1, vmax=vmax_e1)
    ax.set_title(f'Step 1 error (RMSE={np.sqrt(np.mean(err1_3d**2)):.4f} eV)', fontsize=10)
    ax.set_xlabel('x (Å)'); ax.set_ylabel('z (Å)'); ax.set_aspect('equal')
    plt.colorbar(im, ax=ax, label='eV', fraction=0.046); overlay_atoms(ax, ms=6); draw_Rc_circles(ax, R_fold)
    
    # Error2 map
    ax = axes3[1]
    im = ax.pcolormesh(x_3d, z_3d, err2_xz.T, shading='auto', cmap='RdBu_r', vmin=-vmax_e2, vmax=vmax_e2)
    ax.set_title(f'Step 2 error (RMSE={np.sqrt(np.mean(err2_3d**2)):.4f} eV)', fontsize=10)
    ax.set_xlabel('x (Å)'); ax.set_ylabel('z (Å)'); ax.set_aspect('equal')
    plt.colorbar(im, ax=ax, label='eV', fraction=0.046); overlay_atoms(ax, ms=6)
    draw_Rc_circles(ax, R_c_quad, alpha=0.7); draw_Rc_circles(ax, R_fold, alpha=0.3)
    
    # 1D error comparison
    ax = axes3[2]
    ax.plot(z_3d, err1_3d[ix_na, iy, :], 'b-', lw=1.5, label='Step 1 error (above Na)')
    ax.plot(z_3d, err2_3d[ix_na, iy, :], 'r-', lw=1.5, label='Step 2 error (above Na)')
    ix_cl = np.argmin(np.abs(x_3d - 2.0))
    ax.plot(z_3d, err1_3d[ix_cl, iy, :], 'b--', lw=1, alpha=0.5, label='Step 1 error (above Cl)')
    ax.plot(z_3d, err2_3d[ix_cl, iy, :], 'r--', lw=1, alpha=0.5, label='Step 2 error (above Cl)')
    ax.set_xlabel('z (Å)'); ax.set_ylabel('Error (eV)')
    ax.set_title('1D error comparison'); ax.legend(fontsize=8); ax.grid(True, alpha=0.3); ax.axhline(0, color='gray', lw=0.5)
    
    fig3.tight_layout()
    fig3_path = os.path.join(save_dir, 'error_comparison.png')
    fig3.savefig(fig3_path, dpi=150); plt.close(fig3)
    print(f"  Saved {fig3_path}")
    
    # ========== Figure 4: Basis functions ==========
    n_win_plot = 2
    fig4, axes4 = plt.subplots(1, 3, figsize=(18, 5))
    fig4.suptitle(f'Basis Functions: {int(np.sum(top_mask))} folded Morse×window (top L1, R_fold={R_fold} Å, n_win={n_win_plot}) + quadratic p(r) (R_c={R_c_quad} Å)', fontsize=13)
    
    r_plot = np.linspace(0.3, R_fold + 1, 300)
    colors_type = {'Na': 'purple', 'Cl': 'green'}
    
    # Panel 1: Morse*window + V_hardcore (quadratic + Morse*window)
    ax = axes4[0]
    plotted_types = set()
    for j in range(len(apos)):
        if not top_mask[j]: continue
        name = enames[j]
        if name in plotted_types: continue
        plotted_types.add(name)
        R_eq = REQs[j, 0]; D_e = REQs[j, 1]; alpha = REQs[j, 3]
        col = colors_type.get(name, 'black')
        e_arg = alpha * (r_plot - R_eq)
        V_morse_r = D_e * (np.exp(-2*e_arg) - 2*np.exp(-e_arg))
        w_r = np.where(r_plot < R_fold, (1.0 - r_plot / R_fold)**(2*n_win_plot), 0.0)
        V_morse_win = V_morse_r * w_r
        
        mp = match_params[j]
        a_c = mp['a']; b_c = mp['b']; c_c = mp['c']; Rc = mp['R_c']
        V_poly_r = a_c + b_c * r_plot + c_c * r_plot**2
        V_hc_r = np.where(r_plot < Rc, V_poly_r, V_morse_win)
        
        ax.plot(r_plot, V_morse_r, '-', color=col, lw=1, alpha=0.4, label=f'Morse_{name}(r) [raw]')
        ax.plot(r_plot, V_morse_win, '-', color=col, lw=1.5, alpha=0.7, label=f'Morse_{name}×w(r) [windowed]')
        ax.plot(r_plot, V_hc_r, '--', color=col, lw=2, label=f'V_hardcore_{name}: p(r) for r<R_c, Morse×w for R_c≤r<R_fold')
        ax.axvline(Rc, color=col, ls=':', lw=1, alpha=0.5)
    
    ax.axvline(R_c_quad, color='red', ls=':', lw=1.5, label=f'R_c={R_c_quad}')
    ax.axvline(R_fold, color='blue', ls=':', lw=1.5, label=f'R_fold={R_fold}')
    ax.axvspan(0.3, R_c_quad, alpha=0.05, color='orange', label='quadratic region')
    ax.axvspan(R_c_quad, R_fold, alpha=0.03, color='cyan', label='Morse×w region')
    ax.set_xlabel('r (Å)'); ax.set_ylabel('V (eV)')
    ax.set_title(f'Morse×window + quadratic replacement')
    ax.legend(fontsize=6); ax.grid(True, alpha=0.3); ax.axhline(0, color='gray', lw=0.5)
    
    # Panel 2: Window function w(r) = (1-(r/R_fold))^(2n)
    ax = axes4[1]
    for nw in [1, 2, 3, 4]:
        w_demo = np.where(r_plot < R_fold, (1.0 - r_plot / R_fold)**(2*nw), 0.0)
        ax.plot(r_plot, w_demo, lw=1.5, label=f'n_win={nw}: (1-r/R_fold)^{2*nw}')
    ax.axvline(R_fold, color='blue', ls=':', lw=1.5, label=f'R_fold={R_fold}')
    ax.set_xlabel('r (Å)'); ax.set_ylabel('w(r)')
    ax.set_title(f'Compact-support window w(r) = (1-r/R_fold)^(2n)')
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3); ax.axhline(0, color='gray', lw=0.5)
    
    # Panel 3: Quadratic polynomials
    ax = axes4[2]
    plotted_types = set()
    for j in range(len(apos)):
        if not top_mask[j]: continue
        name = enames[j]
        if name in plotted_types: continue
        plotted_types.add(name)
        col = colors_type.get(name, 'black')
        mp = match_params[j]
        a_c = mp['a']; b_c = mp['b']; c_c = mp['c']; Rc = mp['R_c']
        V_poly_r = a_c + b_c * r_plot + c_c * r_plot**2
        mask = r_plot < Rc
        ax.plot(r_plot[mask], V_poly_r[mask], '-', color=col, lw=2, label=f'p_{name}(r) = {a_c:.4f} + {b_c:.4f}*r + {c_c:.6f}*r²')
        ax.plot(r_plot[~mask], V_poly_r[~mask], ':', color=col, lw=0.5, alpha=0.3)
    
    ax.axvline(R_c_quad, color='red', ls=':', lw=1.5, label=f'R_c={R_c_quad}')
    ax.axvspan(0.3, R_c_quad, alpha=0.05, color='orange', label='valid region (r < R_c)')
    ax.set_xlabel('r (Å)'); ax.set_ylabel('p(r) (eV)')
    ax.set_title('Quadratic polynomials p(r) = a + b·r + c·r² (derivative matching at R_c)')
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3); ax.axhline(0, color='gray', lw=0.5)
    
    fig4.tight_layout()
    fig4_path = os.path.join(save_dir, 'basis_functions.png')
    fig4.savefig(fig4_path, dpi=150); plt.close(fig4)
    print(f"  Saved {fig4_path}")
    
    # --- Summary ---
    print(f"\n{'='*60}")
    print(f"Summary")
    print(f"{'='*60}")
    print(f"  V_total range:      [{V_total_3d.min():.4f}, {V_total_3d.max():.4f}] eV")
    print(f"  Step 1 (B-spline+folded Morse):")
    print(f"    V_folded range:   [{V_folded_3d.min():.6f}, {V_folded_3d.max():.6f}] eV")
    print(f"    V_bspline range:  [{V_bspline_3d.min():.4f}, {V_bspline_3d.max():.4f}] eV")
    print(f"    Error RMSE:       {np.sqrt(np.mean(err1_3d**2)):.6f} eV")
    print(f"    Error max:        {np.max(np.abs(err1_3d)):.6f} eV")
    print(f"  Step 2 (quadratic replacement):")
    print(f"    V_hardcore range: [{V_hardcore_3d.min():.6f}, {V_hardcore_3d.max():.6f}] eV")
    print(f"    V_bspline2 range: [{V_bspline2_3d.min():.4f}, {V_bspline2_3d.max():.4f}] eV")
    print(f"    Error RMSE:       {np.sqrt(np.mean(err2_3d**2)):.6f} eV")
    print(f"    Error max:        {np.max(np.abs(err2_3d)):.6f} eV")
    print(f"  Spline: {nx_3d}x{ny_3d}x{nz_3d} grid, step={spline_step} Å → {n_nodes_x}x{n_nodes_z} nodes in XZ")
    print(f"  Folded Morse: {int(np.sum(top_mask))} atoms (top L1), R_fold={R_fold} Å, {len(set(enames[i] for i in range(len(enames)) if top_mask[i]))} types")
    print(f"  Quadratic: {int(np.sum(top_mask))} polynomials, R_c={R_c_quad} Å")
    print(f"  Plots saved to: {save_dir}/")
    
    return {
        'V_total_3d': V_total_3d, 'V_folded_3d': V_folded_3d, 'V_hardcore_3d': V_hardcore_3d,
        'V_bspline_3d': V_bspline_3d, 'V_bspline2_3d': V_bspline2_3d,
        'V_model1_3d': V_model1_3d, 'V_model2_3d': V_model2_3d,
        'err1_3d': err1_3d, 'err2_3d': err2_3d, 'match_params': match_params,
        'grid_coords': coords_3d, 'mask_phys': mask_phys,
    }


# === AUTO-DOC BEGIN ===
# Linear hybrid potential fitting: tricubic B-spline + compact-support radial atomic basis.
# Reference: erf-damped Ewald2D Coulomb + Morse. Boltzmann-weighted least squares.
# See: doc/Topics/FastCollisionSplitNonbond/hybrid_fitting_design_gpt56sol.md
# === AUTO-DOC END ===

test_hybrid_potential_nacl_legacy = test_hybrid_potential_nacl

def _hybrid_cubic(t):
    t=np.asarray(t); return np.stack(((1-t)**3,3*t**3-6*t*t+4,-3*t**3+3*t*t+3*t+1,t**3),-1)/6, np.stack((-3*(1-t)**2,9*t*t-12*t,-9*t*t+6*t+3,3*t*t),-1)/6

def _hybrid_bspline(points, cell, zr, step):
    """Periodic-XY tricubic B-spline value and analytic-gradient design matrix."""
    Lx,Ly=np.linalg.norm(cell[0]),np.linalg.norm(cell[1]); nx=max(4,round(Lx/step)); ny=max(4,round(Ly/step)); dz=float(step); nz=int(np.ceil((zr[1]-zr[0])/dz))+3
    p=np.asarray(points); us=((p[:,0]%Lx)/(Lx/nx),(p[:,1]%Ly)/(Ly/ny),(p[:,2]-zr[0])/dz+1); ii=[np.floor(u).astype(int) for u in us]; ws=[]; ds=[]
    for u,i,h in zip(us,ii,(Lx/nx,Ly/ny,dz)): w,d=_hybrid_cubic(u-i); ws.append(w); ds.append(d/h)
    A=np.zeros((len(p),nx*ny*nz)); G=np.zeros((3,len(p),nx*ny*nz)); rows=np.arange(len(p))
    for a in range(4):
      for b in range(4):
       for c in range(4):
        col=(((ii[0]+a-1)%nx)*ny+((ii[1]+b-1)%ny))*nz+np.clip(ii[2]+c-1,0,nz-1)
        A[rows,col]+=ws[0][:,a]*ws[1][:,b]*ws[2][:,c]
        G[0,rows,col]+=ds[0][:,a]*ws[1][:,b]*ws[2][:,c]; G[1,rows,col]+=ws[0][:,a]*ds[1][:,b]*ws[2][:,c]; G[2,rows,col]+=ws[0][:,a]*ws[1][:,b]*ds[2][:,c]
    return A,G

def _hybrid_atom_basis(points, apos, type_ids, cell, Rc, powers=(2,3,4,5)):
    """Shared finite-support radial basis and gradients; sums every required XY image."""
    p=np.asarray(points); A=np.zeros((len(p),(int(np.max(type_ids))+1)*len(powers))); G=np.zeros((3,*A.shape)); nx=int(np.ceil(Rc/np.linalg.norm(cell[0])))+1; ny=int(np.ceil(Rc/np.linalg.norm(cell[1])))+1
    for ia,pos in enumerate(apos):
     for tx in range(-nx,nx+1):
      for ty in range(-ny,ny+1):
       d=p-(pos+tx*cell[0]+ty*cell[1]); r=np.linalg.norm(d,axis=1); m=r<Rc
       if np.any(m):
        u=1-r[m]/Rc
        for ib,n in enumerate(powers):
         q=2*n; col=type_ids[ia]*len(powers)+ib; A[m,col]+=u**q; dv=-q*u**(q-1)/Rc
         for k in range(3): G[k,m,col]+=dv*d[m,k]/np.maximum(r[m],1e-12)
    return A,G

def _hybrid_reference(points, apos, qs, reqs, cell, n_harm=12, sigma=2.0):
    """NaCl unit-charge probe reference: erf-damped Coulomb + Morse images.
    
    V_ref = Σ_j [ COUL*q_j*erf(r_j/σ)/r_j + Morse_j(r_j) ]
    
    The erf(r/σ)/r damping replaces bare 1/r: at r>>σ it → 1/r (full Coulomb),
    at r<<σ it → 2/(σ√π) (finite, no singularity). This models electron screening
    and makes the Morse repulsion visible at short range.
    
    Forces: analytic Morse + finite-diff on Ewald erf-part + analytic erfc correction.
    """
    from scipy.special import erf, erfc
    COUL=14.3996448915
    ew=Ewald2D(cell[0,:2],cell[1,:2],apos[:,0],apos[:,1],apos[:,2],qs,n_harm=n_harm)
    def e(p): return ew.phi_full_2d(p[:,0],p[:,1],p[:,2])
    V=e(points); F=np.empty((len(points),3)); h=1e-4
    for k in range(3): d=np.zeros(3); d[k]=h; F[:,k]=-(e(points+d)-e(points-d))/(2*h)
    for ia,pos in enumerate(apos):
     R,D,_,alpha=reqs[ia]; q=qs[ia]
     for tx in range(-3,4):
      for ty in range(-3,4):
       dp=points-(pos+tx*cell[0]+ty*cell[1]); r=np.maximum(np.linalg.norm(dp,axis=1),1e-12)
       # Morse
       e1=np.exp(-alpha*(r-R)); e2=e1*e1; V+=D*(e2-2*e1); F-=2*alpha*D*(e1-e2)[:,None]*dp/r[:,None]
       # Replace bare Coulomb 1/r with erf(r/σ)/r: ΔV = -COUL*q*erfc(r/σ)/r
       erfc_r=erfc(r/sigma); V+=COUL*q*(erf(r/sigma)-1.0)/r
       # Analytic force from erfc correction: d/dr[-COUL*q*erfc(r/σ)/r]
       dVdr=COUL*q*(erfc_r/(r*r) + 2.0/(sigma*np.sqrt(np.pi))*np.exp(-(r/sigma)**2)/r)
       F-=dVdr[:,None]*dp/r[:,None]
    return V,F,ew

def _hybrid_grid(cell,zr,h,shift=0.):
    x=(np.arange(0,np.linalg.norm(cell[0]),h)+shift)%np.linalg.norm(cell[0]); y=(np.arange(0,np.linalg.norm(cell[1]),h)+shift)%np.linalg.norm(cell[1]); z=np.arange(zr[0]+shift,zr[1]+1e-12,h); X,Y,Z=np.meshgrid(x,y,z,indexing='ij'); return np.c_[X.ravel(),Y.ravel(),Z.ravel()]

def test_hybrid_potential_nacl(save_dir='results_hybrid', Ng=80, R_cut_hc=4., spline_step=1., sample_step=.5, force_weight=.25, T_eff=300., sigma=2.0, **_ignored):
    """Fit and validate linear B-spline + folded-atom potential against Ewald2D on NaCl."""
    os.makedirs(save_dir,exist_ok=True); s=make_nacl_1x1_l3_system(); apos,qs,reqs,cell=s['apos'],s['qs'],s['REQs'],s['lvec']; types=np.array([0 if e=='Na' else 1 for e in s['enames']]); zr=(s['z_top']+.5,s['z_top']+6.)
    def mask_core(p):
      d=np.full(len(p),np.inf)
      for tx in range(-2,3):
       for ty in range(-2,3): d=np.minimum(d,np.min(np.linalg.norm(p[:,None]-apos[None]-tx*cell[0]-ty*cell[1],axis=2),axis=1))
      return d>=1.2
    train0=_hybrid_grid(cell,zr,sample_step); train=train0[mask_core(train0)]; V,F,ew=_hybrid_reference(train,apos,qs,reqs,cell,sigma=sigma)
    beta=1./(8.617e-5*T_eff); wb=np.minimum(np.exp(-beta*V),100.); swb=np.sqrt(wb)
    B,dB=_hybrid_bspline(train,cell,zr,spline_step); P,dP=_hybrid_atom_basis(train,apos,types,cell,R_cut_hc); A=np.c_[B,P]; dA=np.concatenate((dB,dP),axis=2)
    wf=swb*np.sqrt(force_weight); M=np.vstack((swb[:,None]*A,-wf[:,None]*dA[0],-wf[:,None]*dA[1],-wf[:,None]*dA[2])); rhs=np.r_[swb*V,wf*F[:,0],wf*F[:,1],wf*F[:,2]]
    scale=np.maximum(np.linalg.norm(M,axis=0),1e-14); coef=np.linalg.lstsq(M/scale,rhs,rcond=None)[0]/scale
    valid0=_hybrid_grid(cell,zr,sample_step,.25); valid=valid0[mask_core(valid0)]; Vr,Fr,_=_hybrid_reference(valid,apos,qs,reqs,cell,sigma=sigma); Bv,dBv=_hybrid_bspline(valid,cell,zr,spline_step); Pv,dPv=_hybrid_atom_basis(valid,apos,types,cell,R_cut_hc); Vv=np.c_[Bv,Pv]@coef; Fv=-np.einsum('knp,p->nk',np.concatenate((dBv,dPv),axis=2),coef); E=Vv-Vr; FE=Fv-Fr
    wb_valid=np.minimum(np.exp(-beta*Vr),100.)
    print(f'linear hybrid (T_eff={T_eff}K, beta={beta:.2f} eV^-1): ntrain={len(train)} nvalid={len(valid)} V_RMSE={np.sqrt(np.mean(E*E)):.3e} eV Vmax={np.max(abs(E)):.3e} eV F_RMSE={np.sqrt(np.mean(FE*FE)):.3e} eV/A')
    print(f'  weighted V_RMSE (Boltzmann): {np.sqrt(np.mean(wb_valid*E*E)):.3e} eV  (attractive-only V_RMSE: {np.sqrt(np.mean(E[Vr<0]**2)):.3e} eV)')
    x=np.linspace(0,np.linalg.norm(cell[0]),Ng,endpoint=False); z=np.linspace(*zr,Ng); X,Z=np.meshgrid(x,z,indexing='ij'); cut=np.c_[X.ravel(),X.ravel()%np.linalg.norm(cell[1]),Z.ravel()]; Vref,Fref,_=_hybrid_reference(cut,apos,qs,reqs,cell,sigma=sigma); Bc,dBc=_hybrid_bspline(cut,cell,zr,spline_step); Pc,dPc=_hybrid_atom_basis(cut,apos,types,cell,R_cut_hc); nB=Bc.shape[1]; Vbs=Bc@coef[:nB]; Vat=Pc@coef[nB:]; Vmod=Vbs+Vat; Verr=Vref-Vmod; Fzerr=Fref[:,2]+np.einsum('np,p->n',dBc[2],coef[:nB])+np.einsum('np,p->n',dPc[2],coef[nB:])
    fit_mask=mask_core(cut); Verr_plot=np.where(fit_mask,Verr,np.nan); Fzerr_plot=np.where(fit_mask,Fzerr,np.nan)
    wb_cut=np.minimum(np.exp(-beta*Vref),100.); wb_plot=np.where(fit_mask,wb_cut,np.nan)
    maps=(Vref,Vmod,Verr_plot,Vbs,Vat,Fzerr_plot,wb_plot); titles=(f'Ewald2D(erf σ={sigma}) + Morse ref','linear hybrid fit','potential error (ref-fit)','B-spline component','folded-atom component','Fz error',f'Boltzmann weight (T={T_eff}K)')
    fig,axs=plt.subplots(2,4,figsize=(20,9)); limits=(max(abs(Vref).max(),abs(Vmod).max()),max(abs(Vref).max(),abs(Vmod).max()),max(abs(Verr_plot).max(),1e-12),max(abs(Vbs).max(),1e-12),max(abs(Vat).max(),1e-12),max(abs(Fzerr_plot).max(),1e-12),max(np.nanmax(wb_plot),1e-12))
    for ax,a,title,lim in zip(axs.flat,maps,titles,limits):
      is_wt='Boltzmann' in title
      im=ax.pcolormesh(x,z,a.reshape(Ng,Ng).T,shading='auto',cmap='hot' if is_wt else 'RdBu_r',vmin=0 if is_wt else -lim,vmax=lim); ax.set(title=title,xlabel='x=y diag (A)',ylabel='z (A)'); ax.set_aspect('equal'); plt.colorbar(im,ax=ax,label='weight' if is_wt else ('eV/A' if title=='Fz error' else 'eV'))
    fig.tight_layout(); pmap=os.path.join(save_dir,'linear_hybrid_xz_components.png'); fig.savefig(pmap,dpi=160); plt.close(fig)
    fig,ax=plt.subplots(figsize=(10,5)); sl=slice(0,Ng); ax.plot(z,Vref[sl],'k:',lw=1.5,label='reference'); ax.plot(z,Vmod[sl],'r-',lw=.8,label='fit'); ax.plot(z,Vbs[sl],'b-',lw=.8,label='B-spline'); ax.plot(z,Vat[sl],'g-',lw=.8,label='folded atom'); ax2=ax.twinx(); ax2.plot(z,Verr[sl],'m-',lw=.8); ax2.set_ylabel('error (eV)',color='m'); ax.set(xlabel='z (A)',ylabel='V (eV)',title='Diagonal cut y=x (over Na at x=0, Cl at x=2)'); ax.legend(); ax.grid(alpha=.3); fig.tight_layout(); pscan=os.path.join(save_dir,'linear_hybrid_zscan.png'); fig.savefig(pscan,dpi=160); plt.close(fig)
    return {'coef':coef,'validation':{'V_RMSE':np.sqrt(np.mean(E*E)),'F_RMSE':np.sqrt(np.mean(FE*FE))},'plots':(pmap,pscan),'ewald':ew}
