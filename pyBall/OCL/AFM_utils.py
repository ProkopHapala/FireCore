#!/usr/bin/env python3
"""
AFM_utils.py - High-level AFM utilities, plotting, and orchestration.

This module provides:
- Plotting utilities for AFM fields and images
- High-level orchestration functions for AFM simulation
- Integration with QM backends for density providers
- I/O utilities for saving/loading AFM data

Design principle: AFM.py contains pure physics (no matplotlib, no QM).
This module depends on AFM.py and adds plotting, I/O, and orchestration.
"""

import numpy as np
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

# Import core AFM physics
from . import AFM as afm

# ═══════════════════════════════════════════════════════════════════════════════
# Plotting Utilities (moved from AFM.py)
# ═══════════════════════════════════════════════════════════════════════════════

def safe_norm(data_2d, pct=99):
    """Symmetric ±vabs TwoSlopeNorm for diverging colormaps."""
    vabs = max(float(np.percentile(np.abs(data_2d), pct)), 1e-6)
    return TwoSlopeNorm(vmin=-vabs, vcenter=0, vmax=vabs)


def plot_xy_slice(data, origin, step, iz, title, fname, save_dir, sym=False, cmap='magma'):
    """Plot xy slice at given z-index."""
    import matplotlib.pyplot as plt
    slice_data = data[:, :, iz].T
    nx, ny = data.shape[0], data.shape[1]
    x_min = origin[0]
    x_max = origin[0] + nx * step
    y_min = origin[1]
    y_max = origin[1] + ny * step
    
    fig, ax = plt.subplots(figsize=(6, 5))
    norm = None
    if sym:
        vmax = np.max(np.abs(slice_data))
        norm = plt.Normalize(-vmax, vmax)
    
    im = ax.imshow(slice_data, origin='lower', extent=[x_min, x_max, y_min, y_max], cmap=cmap, norm=norm, aspect='equal')
    ax.set_title(title)
    ax.set_xlabel('x [A]')
    ax.set_ylabel('y [A]')
    plt.colorbar(im, ax=ax)
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, fname), dpi=120, bbox_inches='tight')
    plt.close()
    print(f"  Saved {fname}")


def save_afm_images(df, scan_xs, scan_ys, heights, out_dir, prefix='df'):
    """Save AFM frequency-shift images at all heights.

    Args:
        df: (nx, ny, nz) frequency-shift array
        scan_xs, scan_ys: 1D scan coordinate arrays
        heights: 1D probe height array
        out_dir: directory for PNG output
        prefix: filename prefix (e.g. 'df' -> df_h3.0.png)
    """
    for i in range(len(heights)):
        h = heights[i]
        fig, ax = plt.subplots(figsize=(5,4))
        im = ax.imshow(df[:,:,i].T, origin='lower', extent=[scan_xs[0], scan_xs[-1], scan_ys[0], scan_ys[-1]], cmap='afmhot')
        ax.set_title(f"{prefix} at h={h:.1f} A")
        plt.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
        plt.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.9)
        fname = os.path.join(out_dir, f"{prefix}_h{h:.1f}.png")
        plt.savefig(fname, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  Saved: {fname}")


def plot_slices(data, title, fname, sym=False, cmap='magma', save_dir='.'):
    """Plot central XY/XZ/YZ slices + 1D profiles of a 3D field."""
    nx, ny, nz = data.shape
    cx, cy, cz = nx//2, ny//2, nz//2
    if sym: cmap = 'bwr'
    fig, axes = plt.subplots(2, 3, figsize=(16, 8)); fig.suptitle(title)
    norm = safe_norm(data) if sym else None
    kw = dict(origin='lower', cmap=cmap, aspect='auto', norm=norm)
    for ax, sl, tl in zip(axes[0],
        [data[cx,:,:].T, data[:,cy,:].T, data[:,:,cz].T],
        [f'ix={cx} (YZ)', f'iy={cy} (XZ)', f'iz={cz} (XY)']):
        im = ax.imshow(sl, **kw); ax.set_title(tl); plt.colorbar(im, ax=ax, shrink=0.8)
    axes[1,0].plot(data[cx,cy,:]); axes[1,0].set_xlabel('iz'); axes[1,0].set_title('z-profile center')
    axes[1,1].plot(data[:,cy,cz]); axes[1,1].set_xlabel('ix'); axes[1,1].set_title('x-profile center')
    axes[1,2].plot(data[cx,:,cz]); axes[1,2].set_xlabel('iy'); axes[1,2].set_title('y-profile center')
    for ax in axes[1]: ax.axhline(0, color='k', lw=0.5)
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, fname), dpi=90, bbox_inches='tight'); plt.close()
    print(f"Saved {fname}")


def plot_grid_Fz(Fz, heights, label, fname, x_ext=None, y_ext=None, ncols=7, save_dir='.'):
    """Plot grid of 2D Fz images at all heights with per-slice colorbars."""
    nz_p = len(heights)
    nrows = int(np.ceil(nz_p / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(2.5*ncols, 2.8*nrows))
    axes = np.array(axes).reshape(nrows, ncols)
    fig.suptitle(f"{label} (eV/Å) [per-slice]", fontsize=10)
    ext = [x_ext[0], x_ext[1], y_ext[0], y_ext[1]] if x_ext is not None and y_ext is not None else None
    kw = dict(origin='lower', cmap='bwr', aspect='auto')
    if ext: kw['extent'] = ext
    for k in range(nz_p):
        r, c = divmod(k, ncols); ax = axes[r, c]
        vabs = max(float(np.percentile(np.abs(Fz[:,:,k]), 99)), 1e-6)
        norm = safe_norm(Fz[:,:,k])
        im = ax.imshow(Fz[:,:,k].T, norm=norm, **kw)
        ax.set_title(f"h={heights[k]:.1f}Å ±{vabs:.2g}", fontsize=6); ax.tick_params(labelsize=4)
        plt.colorbar(im, ax=ax, shrink=0.8)
    for k in range(nz_p, nrows*ncols):
        r, c = divmod(k, ncols); axes[r, c].set_visible(False)
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, fname), dpi=90, bbox_inches='tight'); plt.close()
    print(f"Saved {fname}")


# ═══════════════════════════════════════════════════════════════════════════════
# Density Provider Adapters (standard interface)
# ═══════════════════════════════════════════════════════════════════════════════

def _make_grid_spec(atomPos, step, margin, z_extra):
    """Build grid_spec + return (grid_spec, origin, ngrid, step). Wraps afm.setup_density_grid."""
    grid_spec, origin, ngrid = afm.setup_density_grid(atomPos, step=step, margin=margin, z_extra=z_extra)
    return grid_spec, origin, ngrid, step


def _project_densities(geo, evecs, basis, grid_spec, verbosity=0):
    """Shared projection logic: returns (rho_scf, rho_na, rho_diff) using Grid_dftb backends."""
    from pyBall.DFTB import Grid_dftb as dg
    dftb_data = {k: geo[k] for k in ('coords_bohr', 'species_per_atom', 'species_names')}
    projector, atoms_dict = dg.setup_gridprojector_from_dftb(dftb_data, basis, verbosity=verbosity)
    rho_scf = dg.project_dftb_density(geo, evecs, projector, atoms_dict, grid_spec, basis)
    rho_na  = dg.project_neutral_density(geo, projector, atoms_dict, grid_spec, basis)
    return rho_scf, rho_na, (rho_scf - rho_na).astype(np.float32)


def get_density_from_dftb_plus(atomPos, atomTypes, basis, slako_prefix, work_dir,
                                grid_spec=None, step=0.1, margin=4.0, z_extra=6.0, verbosity=0):
    """
    Run DFTB+ SCF for density projection and return density grids.

    Returns dict with 'rho_scf', 'rho_na', 'rho_diff', 'V_ES', 'origin', 'ngrid', 'grid_spec'.
    """
    from pyBall import dftb_utils as du
    ELEM_Z = {'H':1,'C':6,'N':7,'O':8,'P':15,'S':16}
    inv_z  = {v: k for k, v in ELEM_Z.items()}
    enames = [inv_z[int(z)] for z in atomTypes]

    if grid_spec is None:
        grid_spec, origin, ngrid, step = _make_grid_spec(atomPos, step, margin, z_extra)
    else:
        origin, ngrid, step = grid_spec['origin'], grid_spec['ngrid'], grid_spec['dA'][0]

    geo, evecs = du.run_dftb_for_density(work_dir, enames, atomPos, slako_prefix)
    
    # Parse basis HSD file for density projection
    # Use waveplot_in.hsd from DFTB output if it exists (matches actual calculation)
    # Otherwise fall back to pre-defined basis file
    from pyBall.DFTB.DFTBplusParser import parse_basis_hsd_ang
    waveplot_hsd = os.path.join(work_dir, 'waveplot_in.hsd')
    if os.path.exists(waveplot_hsd):
        species_list_ang = parse_basis_hsd_ang(waveplot_hsd)
    else:
        basis_hsd_path = du.WFC_HSD_PATHS.get(basis)
        if basis_hsd_path is None:
            raise ValueError(f"No basis HSD file defined for basis '{basis}'. Available: {list(du.WFC_HSD_PATHS.keys())}")
        species_list_ang = parse_basis_hsd_ang(basis_hsd_path)
    
    # Validate that all atoms in the molecule are present in the basis
    basis_species = set(sp['name'] for sp in species_list_ang)
    molecule_species = set(geo['species_names'])
    missing_species = molecule_species - basis_species
    if missing_species:
        raise ValueError( f"Atoms in molecule not supported by basis '{basis}': {missing_species}.  Basis contains: {sorted(basis_species)}")
    
    rho_scf, rho_na, rho_diff = _project_densities(geo, evecs, species_list_ang, grid_spec, verbosity)
    V_ES = afm.fft_poisson(rho_diff, step)
    return {'rho_scf': rho_scf, 'rho_na': rho_na, 'rho_diff': rho_diff, 'V_ES': V_ES,
            'origin': origin, 'ngrid': ngrid, 'grid_spec': grid_spec}


def get_density_from_dftb(atomPos, atomTypes, dftb_dir, basis=None,
                           grid_spec=None, step=0.15, margin=4.0, z_extra=6.0, verbosity=0):
    """
    Get density grids from pre-computed DFTB+ output files (detailed.xml + eigenvec.bin).

    Returns dict with 'rho_scf', 'rho_na', 'rho_diff', 'V_ES', 'origin', 'ngrid', 'grid_spec'.
    """
    from pyBall.DFTB.DFTBplusParser import parse_detailed_xml_custom, parse_eigenvec_bin_custom, parse_basis_hsd_ang

    if grid_spec is None:
        grid_spec, origin, ngrid, step = _make_grid_spec(atomPos, step, margin, z_extra)
    else:
        origin, ngrid, step = grid_spec['origin'], grid_spec['ngrid'], grid_spec['dA'][0]

    geo   = parse_detailed_xml_custom(os.path.join(dftb_dir, 'detailed.xml'))
    evecs = parse_eigenvec_bin_custom(os.path.join(dftb_dir, 'eigenvec.bin'), geo['nstates'], geo['norb'])

    if basis is None:
        hsd = os.path.join(dftb_dir, 'waveplot_in.hsd')
        if not os.path.exists(hsd):
            raise FileNotFoundError(f"waveplot_in.hsd not found in {dftb_dir}")
        basis = parse_basis_hsd_ang(hsd)

    rho_scf, rho_na, rho_diff = _project_densities(geo, evecs, basis, grid_spec, verbosity)
    V_ES = afm.fft_poisson(rho_diff, step)
    return {'rho_scf': rho_scf, 'rho_na': rho_na, 'rho_diff': rho_diff, 'V_ES': V_ES,
            'origin': origin, 'ngrid': ngrid, 'grid_spec': grid_spec}


def get_density_from_fireball(atomPos, atomTypes, grid_spec, fdata_dir, fc_instance=None, step=0.15, margin=4.0, z_extra=6.0, verbosity=0):
    """
    Get electron density from Fireball SCF.
    
    Args:
        atomPos: (natoms, 3) positions in Angstrom
        atomTypes: (natoms,) atomic numbers
        grid_spec: dict with origin, dA, dB, dC, ngrid (optional, will auto-generate if None)
        fdata_dir: directory with Fireball basis files
        fc_instance: optional FireCore instance (will create if None)
        step: grid spacing in Angstrom (if grid_spec not provided)
        margin: margin around molecule for grid
        z_extra: extra margin in z direction
        verbosity: logging level
        
    Returns:
        dict with 'rho_scf', 'rho_na', 'rho_diff', 'V_ES', 'origin', 'ngrid', 'grid_spec'
    """
    from pyBall.FireballOCL import Grid as ocl_grid
    
    # Auto-generate grid spec if not provided
    if grid_spec is None:
        origin, ngrid, step = afm.setup_density_grid(atomPos, step=step, margin=margin, z_extra=z_extra)
        grid_spec = {
            'origin': origin,
            'dA': [step, 0., 0.], 'dB': [0., step, 0.], 'dC': [0., 0., step],
            'ngrid': ngrid.astype(int),
        }
    else:
        origin = grid_spec['origin']
        ngrid = grid_spec['ngrid']
        step = grid_spec['dA'][0]
    
    if fc_instance is None:
        raise NotImplementedError("Fireball density provider needs FireCore instance to compute SCF and get density matrices")
    
    # Get density from FireCore and project using Grid projector
    # This would require:
    # 1. Get sparse density matrices from FireCore
    # 2. Convert to format expected by Grid projector
    # 3. Project to grid
    
    raise NotImplementedError("Fireball density provider needs implementation with density matrix extraction and projection")


# ═══════════════════════════════════════════════════════════════════════════════
# Orchestration Functions (glue AFM.py physics with I/O and plotting)
# ═══════════════════════════════════════════════════════════════════════════════

def compose_and_relax(grads_pauli, grads_es, grads_vdw, scan_xs, scan_ys, heights,
                     origin, step, atomPos, K_LAT=0.5):
    """
    Compose force fields and run probe particle relaxation to get AFM frequency shift.
    
    This is orchestration - combines AFM.py physics functions with force interpolation.
    
    Args:
        grads_pauli: (nx, ny, nz, 3) Pauli gradients from afm.compute_pauli_field
        grads_es: (nx, ny, nz, 3) Electrostatic gradients from afm.compute_es_conv_field
        grads_vdw: (nx, ny, nz, 3) Dispersion gradients from afm.compute_dispersion_grid
        scan_xs: (nx_s,) scan x coordinates
        scan_ys: (ny_s,) scan y coordinates
        heights: (nz_s,) probe heights
        origin: (3,) grid origin
        step: grid spacing
        atomPos: (natoms, 3) atom positions (for mol_z)
        K_LAT: lateral stiffness
        
    Returns:
        df: (nx_s, ny_s, nz_s) frequency shift array
    """
    from scipy.ndimage import map_coordinates
    
    # F = -grad E
    F_total = -(grads_pauli + grads_es + grads_vdw)
    
    def force_func(positions):
        """Interpolate forces at arbitrary positions from scan-grid force field."""
        ix = (positions[:, 0] - origin[0]) / step
        iy = (positions[:, 1] - origin[1]) / step
        iz = (positions[:, 2] - origin[2]) / step
        coords = np.vstack([ix, iy, iz])
        fx = map_coordinates(F_total[..., 0], coords, order=1)
        fy = map_coordinates(F_total[..., 1], coords, order=1)
        fz = map_coordinates(F_total[..., 2], coords, order=1)
        return np.stack([fx, fy, fz], axis=-1)
    
    mol_z = atomPos[:,2].max()
    FEs_relax = afm.pp_relax_2d(force_func, scan_xs, scan_ys, heights, mol_z=mol_z, K_LAT=K_LAT, N_RELAX=50, step=step)
    df = afm.compute_df(FEs_relax[:,:,:,2], heights[1]-heights[0])
    return df



# ═══════════════════════════════════════════════════════════════════════════════
# Step Plotting Functions (separate from computation)
# ═══════════════════════════════════════════════════════════════════════════════

def plot_step1_outputs(rho_grid, rho_na_grid, rho_diff, step_dir, origin, step):
    """Plot step 1 density outputs."""
    from pyBall import plotUtils as pu
    z_slice = 2.0
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle('Step 1: Density Projection')
    
    pu.plot_field_slice(axes[0], rho_grid, origin, step, z_slice, cmap='magma', 
                       title='SCF Density [e/Å³]')
    pu.plot_field_slice(axes[1], rho_na_grid, origin, step, z_slice, cmap='magma',
                       title='Neutral Atom Density [e/Å³]')
    pu.plot_field_slice(axes[2], rho_diff, origin, step, z_slice, cmap='bwr', sym=True,
                       title='Delta Density [e/Å³]')
    
    plt.tight_layout()
    plt.savefig(os.path.join(step_dir, 'step1_rho_slices.png'), dpi=120, bbox_inches='tight')
    plt.close()
    
    print(f"  Saved step1 density plots")


def plot_step2_outputs(V_ES, step_dir, origin, step):
    """Plot step 2 electrostatics outputs."""
    from pyBall import plotUtils as pu
    z_slice = 2.0
    
    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    pu.plot_field_slice(ax, V_ES, origin, step, z_slice, cmap='bwr', sym=True,
                       title='Electrostatic Potential [eV]')
    plt.tight_layout()
    plt.savefig(os.path.join(step_dir, 'step2_VES_slices.png'), dpi=120, bbox_inches='tight')
    plt.close()
    print(f"  Saved step2 electrostatics plot")


def plot_step3_outputs(E_pauli_field, grads_pauli, step_dir, origin, step, A_pauli, beta_pauli):
    """Plot step 3 Pauli outputs."""
    from pyBall import plotUtils as pu
    z_slice = 2.0
    
    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    pu.plot_field_slice(ax, E_pauli_field, origin, step, z_slice, cmap='magma',
                       title=f'Pauli Energy [eV] (A={A_pauli:.1f}, b={beta_pauli:.3f})')
    plt.tight_layout()
    plt.savefig(os.path.join(step_dir, 'step3_Epauli_slices.png'), dpi=120, bbox_inches='tight')
    plt.close()
    print(f"  Saved step3 Pauli plot")


def plot_step4_outputs(E_ES_field, grads_ES, step_dir, origin, step):
    """Plot step 4 electrostatics convolution outputs."""
    from pyBall import plotUtils as pu
    z_slice = 2.0
    
    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    pu.plot_field_slice(ax, E_ES_field, origin, step, z_slice, cmap='bwr', sym=True,
                       title='ES Energy [eV]')
    plt.tight_layout()
    plt.savefig(os.path.join(step_dir, 'step4_EES_slices.png'), dpi=120, bbox_inches='tight')
    plt.close()
    print(f"  Saved step4 ES convolution plot")


def plot_step5_outputs(E_vdw, grads_vdw, step_dir, origin, step):
    """Plot step 5 dispersion outputs."""
    from pyBall import plotUtils as pu
    z_slice = 2.0
    
    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    pu.plot_field_slice(ax, E_vdw, origin, step, z_slice, cmap='magma',
                       title='Dispersion Energy [eV]')
    plt.tight_layout()
    plt.savefig(os.path.join(step_dir, 'step5_Evdw_slices.png'), dpi=120, bbox_inches='tight')
    plt.close()
    print(f"  Saved step5 dispersion plot")


def plot_step6_outputs(df, scan_xs, scan_ys, heights, step_dir):
    """Plot step 6 final AFM images."""
    save_afm_images(df, scan_xs, scan_ys, heights, step_dir, prefix='df')
    print(f"  Saved step6 AFM images")


# ═══════════════════════════════════════════════════════════════════════════════
# I/O Utilities
# ═══════════════════════════════════════════════════════════════════════════════

def save_grid_spec(grid_spec, step_dir):
    """Save grid specification to file."""
    np.save(os.path.join(step_dir, 'origin.npy'), grid_spec['origin'])
    np.save(os.path.join(step_dir, 'ngrid.npy'), grid_spec['ngrid'])
    with open(os.path.join(step_dir, 'step.txt'), 'w') as f:
        f.write(str(grid_spec['dA'][0]))  # step is same for all axes


def load_grid_spec(step_dir):
    """Load grid specification from file."""
    origin = np.load(os.path.join(step_dir, 'origin.npy'))
    ngrid = np.load(os.path.join(step_dir, 'ngrid.npy'))
    with open(os.path.join(step_dir, 'step.txt'), 'r') as f:
        step = float(f.read().strip())
    grid_spec = {
        'origin': origin,
        'dA': [step, 0., 0.], 'dB': [0., step, 0.], 'dC': [0., 0., step],
        'ngrid': ngrid.astype(int),
    }
    return grid_spec, origin, step


# ═══════════════════════════════════════════════════════════════════════════════
# High-Level Orchestration Function
# ═══════════════════════════════════════════════════════════════════════════════

def run_afm_pipeline(
    rho_grid, rho_na_grid, rho_diff, V_ES,
    rho_tip_total, rho_tip_delta,
    atomPos, atomTypes,
    origin, step, ngrid,
    scan_xs, scan_ys, heights,
    output_dir,
    pauli_params={'A': None, 'beta': None},
    pauli_fit_params=None,
    fit_pauli=False,
    fit_pauli_params=None,  # Dict with zscan_dir, target_indices, z_min, z_max, basis
    vdw_params={'C6_CO': 30.0},
    relax_params={'K_LAT': 0.5},
    plot_steps=True
):
    """
    High-level AFM simulation pipeline using pre-computed densities.
    
    This function runs steps 2-6 of the AFM simulation, assuming step 1
    (density projection) has already been done separately.
    
    Args:
        rho_grid: (nx, ny, nz) sample SCF density
        rho_na_grid: (nx, ny, nz) neutral atom density
        rho_diff: (nx, ny, nz) delta density
        V_ES: (nx, ny, nz) electrostatic potential (optional, can compute)
        rho_tip_total: (nx, ny, nz) CO tip total density
        rho_tip_delta: (nx, ny, nz) CO tip delta density
        atomPos: (natoms, 3) atom positions
        atomTypes: (natoms,) atomic numbers
        origin: (3,) grid origin
        step: grid spacing
        ngrid: (3,) grid dimensions
        scan_xs: (nx_s,) scan x coordinates
        scan_ys: (ny_s,) scan y coordinates
        heights: (nz_s,) probe heights
        output_dir: directory for outputs
        pauli_params: dict with 'A', 'beta'
        pauli_fit_params: dict with fitted 'A', 'beta' (from external fit)
        fit_pauli: if True, fit Pauli parameters internally after computing raw overlap
        fit_pauli_params: dict with 'zscan_dir', 'target_indices', 'z_min', 'z_max', 'basis'
        vdw_params: dict with 'C6_CO'
        relax_params: dict with 'K_LAT'
        plot_steps: whether to generate plots
        
    Returns:
        dict with 'df', 'intermediates', 'grid_spec'
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Step 2: Electrostatics (if V_ES not provided)
    if V_ES is None:
        print("\nStep 2: Computing electrostatics...")
        V_ES = afm.fft_poisson(rho_diff, step)
        if plot_steps:
            plot_step2_outputs(V_ES, output_dir, origin, step)
        np.save(os.path.join(output_dir, 'V_ES.npy'), V_ES)
    else:
        print("\nStep 2: Using provided V_ES")
        if plot_steps:
            plot_step2_outputs(V_ES, output_dir, origin, step)
    
    # Step 3a: Compute raw Pauli overlap (A=1, beta=1 — pure density convolution)
    print("\nStep 3a: Computing raw Pauli overlap (A=1, beta=1)...")
    overlap_raw = afm.compute_pauli_overlap(rho_grid, rho_tip_total, step, tip_rolled=True)
    np.save(os.path.join(output_dir, 'overlap_raw.npy'), overlap_raw)
    print(f"  overlap_raw: shape={overlap_raw.shape}  range=[{overlap_raw.min():.4e}, {overlap_raw.max():.4e}]")
    
    # Step 3b: Fit Pauli parameters or use provided / default values
    A_pauli = pauli_params.get('A') if pauli_params else None
    beta_pauli = pauli_params.get('beta') if pauli_params else None
    
    if fit_pauli and fit_pauli_params is not None:
        # Fit Pauli parameters internally using DFTB reference
        print("\nStep 3b: Fitting Pauli parameters using DFTB z-scan reference...")
        zscan_dir = fit_pauli_params['zscan_dir']
        target_indices = fit_pauli_params['target_indices']
        z_min = fit_pauli_params.get('z_min', 2.0)
        z_max = fit_pauli_params.get('z_max', 3.0)
        
        from pyBall.DFTB import TestUtils as tu
        
        # Load DFTB reference for each target atom
        all_A, all_beta = [], []
        for idx in target_indices:
            atom_dir = os.path.join(zscan_dir, f'atom_{idx}')
            z_ref = np.load(os.path.join(atom_dir, 'zscan_z.npy'))
            e_ref_abs = np.load(os.path.join(atom_dir, 'zscan_energy_eV.npy'))
            e_ref = e_ref_abs - e_ref_abs[-1]  # Reference to far distance
            
            target_pos = atomPos[idx]
            overlap_profile = tu.extract_z_profile(overlap_raw, target_pos, origin, step, z_distances=z_ref)
            
            # Fit power law (returns tuple: A, beta, r2, e_pred)
            A_fit, beta_fit, r2_fit, _ = _fit_pauli_powerlaw(z_ref, overlap_profile, e_ref, z_min, z_max)
            all_A.append(A_fit)
            all_beta.append(beta_fit)
            
            print(f"  Atom {idx}: A={A_fit:.2f}, beta={beta_fit:.4f}, R2={r2_fit:.4f}")
        
        A_pauli = np.mean(all_A)
        beta_pauli = np.mean(all_beta)
        print(f"\nStep 3b: Fitted Pauli params: A={A_pauli:.4f}, beta={beta_pauli:.4f}")
    elif pauli_fit_params is not None:
        # Fit was done externally; use those results
        A_pauli   = pauli_fit_params['A']
        beta_pauli = pauli_fit_params['beta']
        print(f"\nStep 3b: Using externally fitted Pauli params: A={A_pauli:.4f}, beta={beta_pauli:.4f}")
    elif A_pauli is None or beta_pauli is None:
        raise ValueError(
            "Pauli parameters A and beta must be provided via pauli_params or pauli_fit_params, "
            "or set fit_pauli=True with fit_pauli_params."
        )
    else:
        print(f"\nStep 3b: Using provided Pauli params: A={A_pauli:.4f}, beta={beta_pauli:.4f}")
    
    # Step 3c: Scale overlap into energy field
    print(f"\nStep 3c: Scaling E_pauli = {A_pauli:.4f} * overlap^{beta_pauli:.4f}")
    E_pauli_field, grads_pauli = afm.scale_pauli_field(overlap_raw, step, A_pauli, beta_pauli)
    
    # Consistency diagnostics
    print(f"  overlap_raw at max: {overlap_raw.max():.4e}")
    print(f"  E_pauli_field: range=[{E_pauli_field.min():.4e}, {E_pauli_field.max():.4e}]")
    print(f"  Check: A*overlap_max^beta = {A_pauli:.4f}*{overlap_raw.max():.4e}^{beta_pauli:.4f} = {A_pauli * float(overlap_raw.max())**beta_pauli:.4e}")
    
    if plot_steps:
        plot_step3_outputs(E_pauli_field, grads_pauli, output_dir, origin, step, A_pauli, beta_pauli)
    np.save(os.path.join(output_dir, 'E_Pauli_field.npy'), E_pauli_field)
    np.save(os.path.join(output_dir, 'grads_E_Pauli.npy'), grads_pauli)
    
    # Step 4: Electrostatic convolution
    print("\nStep 4: Computing electrostatic convolution...")
    E_ES_field, grads_ES = afm.compute_es_conv_field(V_ES, rho_tip_delta, step, tip_rolled=True)
    if plot_steps:
        plot_step4_outputs(E_ES_field, grads_ES, output_dir, origin, step)
    np.save(os.path.join(output_dir, 'E_ES_field.npy'), E_ES_field)
    np.save(os.path.join(output_dir, 'grads_E_ES.npy'), grads_ES)
    
    # Step 5: Dispersion
    print("\nStep 5: Computing dispersion...")
    E_vdw, grads_vdw = afm.compute_dispersion_grid(
        atomPos, atomTypes, origin, step, ngrid,
        C6_CO=vdw_params['C6_CO']
    )
    if plot_steps:
        plot_step5_outputs(E_vdw, grads_vdw, output_dir, origin, step)
    np.save(os.path.join(output_dir, 'E_vdw_field.npy'), E_vdw)
    np.save(os.path.join(output_dir, 'grads_E_vdw.npy'), grads_vdw)
    
    # Step 6: Compose and relax
    print("\nStep 6: Composing force fields and running probe relaxation...")
    df = compose_and_relax(
        grads_pauli, grads_ES, grads_vdw,
        scan_xs, scan_ys, heights,
        origin, step, atomPos, K_LAT=relax_params['K_LAT']
    )
    np.save(os.path.join(output_dir, 'df.npy'), df)
    if plot_steps:
        plot_step6_outputs(df, scan_xs, scan_ys, heights, output_dir)
    
    # Return results
    grid_spec_out = {
        'origin': origin,
        'dA': [step, 0., 0.], 'dB': [0., step, 0.], 'dC': [0., 0., step],
        'ngrid': ngrid.astype(int),
    }
    
    return {
        'df': df,
        'intermediates': {
            'V_ES': V_ES,
            'E_pauli_field': E_pauli_field,
            'grads_pauli': grads_pauli,
            'E_ES_field': E_ES_field,
            'grads_ES': grads_ES,
            'E_vdw': E_vdw,
            'grads_vdw': grads_vdw,
        },
        'grid_spec': grid_spec_out,
    }


def _compute_co_tip_grid(step=0.1, margin=4.0):
    """Return grid_spec for CO tip computation with O at grid center."""
    co_span = np.array([0.0, 0.0, 1.13])  # C is at z=1.13 relative to O
    ngrid = np.ceil((2 * margin + co_span) / step).astype(np.int32)
    # Round up to nearest multiple of 8 for GPU
    ngrid = ((ngrid + 7) // 8) * 8
    origin = np.array([-margin, -margin, -margin], dtype=np.float32)
    grid_spec = {
        'origin': origin,
        'dA': np.array([step, 0.0, 0.0], dtype=np.float32),
        'dB': np.array([0.0, step, 0.0], dtype=np.float32),
        'dC': np.array([0.0, 0.0, step], dtype=np.float32),
        'ngrid': ngrid,
    }
    return grid_spec, ngrid, origin


def _call_compute_co_tip_script(out_dir, grid_spec, step, nscf, fdata_dir, fdata_basis):
    """Call compute_co_tip.py as subprocess."""
    import json, subprocess, sys
    _THIS_FILE = os.path.abspath(__file__)
    # __file__ is pyBall/OCL/AFM_utils.py; repo root is 2 levels up
    repo_root = os.path.normpath(os.path.join(os.path.dirname(_THIS_FILE), '..', '..'))
    script = os.path.join(repo_root, 'tests', 'tAFM', 'pyocl_fdbm', 'compute_co_tip.py')

    # Convert numpy arrays to lists for JSON serialization
    grid_spec_json = {k: (v.tolist() if hasattr(v, 'tolist') else v) for k, v in grid_spec.items()}
    grid_spec_str = json.dumps(grid_spec_json)

    cmd = [sys.executable, script, out_dir, grid_spec_str, str(step), str(nscf), fdata_dir, fdata_basis]
    print(f"  Running CO tip computation: {' '.join(cmd[:5])} ...")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"compute_co_tip.py failed:\n{result.stderr}\n{result.stdout}")
    print(result.stdout)
    return out_dir


def _pad_and_roll_co_tip(co_rho, target_shape):
    """Pad CO density with zeros to target grid and roll so O atom is at index 0.

    CO is centered in the target grid, then rolled by n//2 so O (at center)
    moves to index 0, matching FFT convolution convention.
    """
    nx_t, ny_t, nz_t = target_shape
    nx_c, ny_c, nz_c = co_rho.shape

    # Center CO in target grid
    ox = (nx_t - nx_c) // 2
    oy = (ny_t - ny_c) // 2
    oz = (nz_t - nz_c) // 2

    padded = np.zeros(target_shape, dtype=np.float32)
    padded[ox:ox+nx_c, oy:oy+ny_c, oz:oz+nz_c] = co_rho

    # Roll so O atom (at center) moves to index 0
    padded = np.roll(padded, -(nx_t // 2), axis=0)
    padded = np.roll(padded, -(ny_t // 2), axis=1)
    padded = np.roll(padded, -(nz_t // 2), axis=2)

    return padded


def _plot_co_tip_diagnostics(co_rho_total, co_rho_delta, output_dir, origin, step, title_suffix=""):
    """Plot diagnostic slices of padded+rolled CO tip density."""
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle(f'CO Tip Diagnostics {title_suffix}')

    def _plot_slice(ax, data, axis, idx, title):
        if axis == 0:
            sl = data[idx, :, :].T
            exts = [origin[1], origin[1] + data.shape[1]*step, origin[2], origin[2] + data.shape[2]*step]
            xl, yl = 'y [A]', 'z [A]'
        elif axis == 1:
            sl = data[:, idx, :].T
            exts = [origin[0], origin[0] + data.shape[0]*step, origin[2], origin[2] + data.shape[2]*step]
            xl, yl = 'x [A]', 'z [A]'
        else:
            sl = data[:, :, idx].T
            exts = [origin[0], origin[0] + data.shape[0]*step, origin[1], origin[1] + data.shape[1]*step]
            xl, yl = 'x [A]', 'y [A]'
        im = ax.imshow(sl, origin='lower', extent=exts, cmap='magma')
        ax.set_title(title)
        ax.set_xlabel(xl)
        ax.set_ylabel(yl)
        plt.colorbar(im, ax=ax, fraction=0.03)

    nx, ny, nz = co_rho_total.shape
    # Slice through origin (0,0,0) where oxygen should be after roll
    ix, iy, iz = 0, 0, 0

    _plot_slice(axes[0, 0], co_rho_total, 0, ix, f'Total YZ (ix={ix} - through origin)')
    _plot_slice(axes[0, 1], co_rho_total, 1, iy, f'Total XZ (iy={iy} - through origin)')
    _plot_slice(axes[0, 2], co_rho_total, 2, iz, f'Total XY (iz={iz} - through origin)')

    _plot_slice(axes[1, 0], co_rho_delta, 0, ix, f'Delta YZ (ix={ix} - through origin)')
    _plot_slice(axes[1, 1], co_rho_delta, 1, iy, f'Delta XZ (iy={iy} - through origin)')
    _plot_slice(axes[1, 2], co_rho_delta, 2, iz, f'Delta XY (iz={iz} - through origin)')

    fname = os.path.join(output_dir, f'co_tip_diagnostics{title_suffix.replace(" ", "_")}.png')
    plt.tight_layout()
    plt.savefig(fname, dpi=120, bbox_inches='tight')
    plt.close()
    print(f"  CO tip diagnostic plot: {fname}")


def run_afm_from_xyz(
    xyz_file,
    output_dir,
    basis,
    slako_prefix='mio-1-1',
    co_tip_dir=None,
    fdata_dir=None,
    fdata_basis=None,
    work_dir=None,
    step=0.1, margin=4.0, z_extra=6.0,
    scan_range=3.0, scan_step=0.1,
    height_range=(3.0, 6.5), height_step=0.1,
    pauli_params=None,
    pauli_fit_params=None,
    fit_pauli=False,
    fit_pauli_params=None,
    vdw_params={'C6_CO': 30.0},
    relax_params={'K_LAT': 0.5},
    plot_steps=True
):
    """
    Full AFM simulation pipeline from .xyz to AFM images via DFTB+ density.

    Args:
        xyz_file: path to .xyz file
        output_dir: all outputs go here
        basis: basis list from parse_basis_hsd_ang (required)
        slako_prefix: Slater-Koster prefix for DFTB+
        co_tip_dir: directory with co_rho_total.npy + co_rho_delta.npy (optional;
                    if not provided or missing, CO is computed on-the-fly)
        fdata_dir: Fireball Fdata directory (required if co_tip_dir not provided)
        fdata_basis: OpenCL basis directory (required if co_tip_dir not provided)
        work_dir: DFTB+ scratch dir (default: output_dir/dftb_work)
        step/margin/z_extra: grid parameters
        scan_range/scan_points/height_range/height_step: scan parameters
        pauli_params/vdw_params/relax_params: physics parameters
        plot_steps: save intermediate plots

    Returns:
        dict with 'df', 'intermediates', 'grid_spec'
    """
    import pyBall.atomicUtils as au
    ELEM_Z = {'H':1,'C':6,'N':7,'O':8,'P':15,'S':16}

    os.makedirs(output_dir, exist_ok=True)
    if work_dir is None:
        work_dir = os.path.join(output_dir, 'dftb_work')

    # Load molecule
    print(f"\nLoading molecule from {xyz_file}")
    pos, _, names, _, _ = au.load_xyz(xyz_file)
    atomPos  = np.array(pos, dtype=np.float64)
    atomTypes = np.array([ELEM_Z.get(e, 6) for e in names], dtype=np.int32)
    print(f"  {len(atomPos)} atoms")

    # Scan grid (compute points from step size)
    x_min, x_max = atomPos[:,0].min()-scan_range, atomPos[:,0].max()+scan_range
    y_min, y_max = atomPos[:,1].min()-scan_range, atomPos[:,1].max()+scan_range
    scan_points_x = int(np.ceil((x_max - x_min) / scan_step))
    scan_points_y = int(np.ceil((y_max - y_min) / scan_step))
    scan_xs = np.linspace(x_min, x_max, scan_points_x)
    scan_ys = np.linspace(y_min, y_max, scan_points_y)
    heights  = np.arange(height_range[0], height_range[1], height_step)

    # Set up Slater-Koster path
    from pyBall import dftb_utils as du
    if slako_prefix == 'mio-1-1':
        slako_prefix = du.SK_PATHS.get('mio-1-1', slako_prefix)
    elif slako_prefix == '3ob-3-1':
        slako_prefix = du.SK_PATHS.get('3ob-3-1', slako_prefix)
    
    # Get densities from DFTB+
    if work_dir is None:
        work_dir = os.path.join(output_dir, 'dftb_work')
    d = get_density_from_dftb_plus(atomPos, atomTypes, basis, slako_prefix, work_dir,
                                    step=step, margin=margin, z_extra=z_extra)

    # Plot density slices to check anisotropy
    if plot_steps:
        z_heights = [0.0, 2.0, 2.5]
        for z in z_heights:
            iz = int(np.clip(np.round((z - d['origin'][2]) / step), 0, d['rho_scf'].shape[2]-1))
            plot_xy_slice(d['rho_scf'], d['origin'], step, iz, f'SCF Density z={z}A', f'step1_rho_scf_z{z:.1f}.png', output_dir)
            plot_xy_slice(d['rho_na'], d['origin'], step, iz, f'Neutral Atom Density z={z}A', f'step1_rho_na_z{z:.1f}.png', output_dir)
            plot_xy_slice(d['rho_diff'], d['origin'], step, iz, f'Delta Density z={z}A', f'step1_rho_diff_z{z:.1f}.png', output_dir, sym=True, cmap='bwr')

    # Save grid spec for later fitting
    grid_spec_path = os.path.join(output_dir, 'grid_spec.txt')
    with open(grid_spec_path, 'w') as f:
        f.write(f"origin = {d['origin'].tolist()}\n")
        f.write(f"ngrid = {d['ngrid'].tolist()}\n")
        f.write(f"step = {step}\n")
    print(f"  Saved grid spec to {grid_spec_path}")

    # CO tip: load precomputed or compute on-the-fly
    target_shape = tuple(d['ngrid'])
    co_origin = None
    if co_tip_dir is not None and os.path.isdir(co_tip_dir):
        print(f"\nLoading precomputed CO tip from {co_tip_dir}...")
        co_rho_total_raw = np.load(os.path.join(co_tip_dir, 'co_rho_total.npy'))
        co_rho_delta_raw = np.load(os.path.join(co_tip_dir, 'co_rho_delta.npy'))
        print(f"  Raw CO tip shape: {co_rho_total_raw.shape}")
    else:
        print(f"\nComputing CO tip on-the-fly (step={step})...")
        if fdata_dir is None or fdata_basis is None:
            _ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
            fdata_dir = fdata_dir or os.path.join(_ROOT, 'tests', 'pyFireball', 'Fdata')
            fdata_basis = fdata_basis or os.path.join(fdata_dir, 'basis')
        co_tip_work = os.path.join(output_dir, 'co_tip_work')
        os.makedirs(co_tip_work, exist_ok=True)
        co_grid_spec, co_ngrid, co_origin = _compute_co_tip_grid(step=step, margin=margin)
        print(f"  CO grid: ngrid={co_ngrid}, origin={co_origin}")
        _call_compute_co_tip_script(co_tip_work, co_grid_spec, step, 100, fdata_dir, fdata_basis)
        co_rho_total_raw = np.load(os.path.join(co_tip_work, 'co_rho_total.npy'))
        co_rho_delta_raw = np.load(os.path.join(co_tip_work, 'co_rho_delta.npy'))
        print(f"  Raw CO tip shape: {co_rho_total_raw.shape}")

    # Pad with zeros and roll so O atom is at index 0
    print(f"  Padding CO tip to target shape {target_shape}...")
    co_rho_total = _pad_and_roll_co_tip(co_rho_total_raw, target_shape)
    co_rho_delta = _pad_and_roll_co_tip(co_rho_delta_raw, target_shape)
    print(f"  Padded+rolled CO tip shape: {co_rho_total.shape}")

    # Diagnostic plots
    if plot_steps:
        co_diag_dir = os.path.join(output_dir, 'co_tip_diagnostics')
        os.makedirs(co_diag_dir, exist_ok=True)
        if co_origin is not None:
            _plot_co_tip_diagnostics(co_rho_total_raw, co_rho_delta_raw, co_diag_dir, co_origin, step, title_suffix="_raw")
        _plot_co_tip_diagnostics(co_rho_total, co_rho_delta, co_diag_dir, d['origin'], step, title_suffix="_padded_rolled")
        # Also save central profiles to verify symmetry
        nx, ny, nz = co_rho_total.shape
        cx, cy, cz = nx // 2, ny // 2, nz // 2
        fig, axes = plt.subplots(1, 3, figsize=(15, 4))
        xs = np.arange(nx) * step + d['origin'][0]
        axes[0].plot(xs, co_rho_total[:, cy, cz], 'b-', label='x profile')
        axes[0].axvline(xs[cx], color='r', ls='--', label='center')
        axes[0].set_title('X profile through center')
        axes[0].set_xlabel('x [A]')
        axes[0].legend()
        ys = np.arange(ny) * step + d['origin'][1]
        axes[1].plot(ys, co_rho_total[cx, :, cz], 'g-', label='y profile')
        axes[1].axvline(ys[cy], color='r', ls='--', label='center')
        axes[1].set_title('Y profile through center')
        axes[1].set_xlabel('y [A]')
        axes[1].legend()
        zs = np.arange(nz) * step + d['origin'][2]
        axes[2].plot(zs, co_rho_total[cx, cy, :], 'm-', label='z profile')
        axes[2].axvline(zs[cz], color='r', ls='--', label='center')
        axes[2].set_title('Z profile through center')
        axes[2].set_xlabel('z [A]')
        axes[2].legend()
        plt.tight_layout()
        prof_path = os.path.join(co_diag_dir, 'co_tip_center_profiles.png')
        plt.savefig(prof_path, dpi=120, bbox_inches='tight')
        plt.close()
        print(f"  CO tip profiles: {prof_path}")

    return run_afm_pipeline(
        d['rho_scf'], d['rho_na'], d['rho_diff'], d['V_ES'],
        co_rho_total, co_rho_delta,
        atomPos, atomTypes,
        d['origin'], step, d['ngrid'],
        scan_xs, scan_ys, heights,
        output_dir,
        pauli_params=pauli_params, pauli_fit_params=pauli_fit_params,
        fit_pauli=fit_pauli, fit_pauli_params=fit_pauli_params,
        vdw_params=vdw_params, relax_params=relax_params, plot_steps=plot_steps
    )


def plot_diagnostic_panel(E_pauli, E_es, E_vdw, E_total, origin, step, heights, output_dir):
    """Plot diagnostic panel with 4 columns (Total, Pauli, Electrostatics, vdW) and n-rows for heights.

    Each subplot has symmetric vmin/vmax zero-centered with its own colorbar (seismic colormap).
    Shows field slices at z=0 (molecular plane) for all heights to show field structure.
    """
    import matplotlib.pyplot as plt
    from scipy.ndimage import map_coordinates
    n_heights = len(heights)
    fig, axes = plt.subplots(n_heights, 4, figsize=(16, 3*n_heights))
    if n_heights == 1:
        axes = axes.reshape(1, -1)

    # Use z=0 slice (molecular plane) for all heights to show field structure
    iz_0 = int(np.clip(np.round((0.0 - origin[2]) / step), 0, E_total.shape[2]-1))
    
    for iz, z in enumerate(heights):
        # Compute actual z-index from physical z coordinate
        iz_grid = int(np.clip(np.round((z - origin[2]) / step), 0, E_total.shape[2]-1))
        for icol, (field, title) in enumerate([
            (E_total, 'Total'),
            (E_pauli, 'Pauli'),
            (E_es, 'Electrostatics'),
            (E_vdw, 'vdW'),
        ]):
            ax = axes[iz, icol]
            slice_data = field[:, :, iz_grid]
            vmax = np.max(np.abs(slice_data))
            im = ax.imshow(slice_data.T, origin='lower', cmap='seismic', vmin=-vmax, vmax=vmax)
            ax.set_title(f'{title} z={z:.1f}Å iz={iz_grid}')
            plt.colorbar(im, ax=ax, fraction=0.03, pad=0.02)

    plt.subplots_adjust(left=0.02, right=0.98, bottom=0.02, top=0.95, wspace=0.25, hspace=0.2)
    plt.savefig(os.path.join(output_dir, 'diagnostic_panel.png'), dpi=120, bbox_inches='tight')
    plt.close()
    print(f"Saved diagnostic panel: {os.path.join(output_dir, 'diagnostic_panel.png')}")


def plot_diagnostic_slices(E_pauli, E_es, E_vdw, origin, step, heights, output_dir):
    """Plot 3x3 diagnostic: Pauli, ES, vdW with XY, XZ, YZ slices through origin.

    All slices pass through origin (0,0,0) to show field structure.
    Probe heights are marked with gray dotted lines on XZ and YZ slices.
    """
    import matplotlib.pyplot as plt
    
    fig, axes = plt.subplots(3, 3, figsize=(15, 15))
    fig.suptitle('Energy Field Slices Through Origin (0,0,0)')
    
    # Slice indices through origin
    ix = 0
    iy = 0
    iz = int(np.clip(np.round((0.0 - origin[2]) / step), 0, E_pauli.shape[2]-1))
    
    fields = [(E_pauli, 'Pauli'), (E_es, 'Electrostatics'), (E_vdw, 'vdW')]
    
    for row, (field, title) in enumerate(fields):
        # XY slice
        xy_slice = field[:, :, iz].T
        vmax = np.max(np.abs(xy_slice))
        x_min = origin[0]
        x_max = origin[0] + field.shape[0] * step
        y_min = origin[1]
        y_max = origin[1] + field.shape[1] * step
        im = axes[row, 0].imshow(xy_slice, origin='lower', extent=[x_min, x_max, y_min, y_max], 
                                 cmap='seismic', vmin=-vmax, vmax=vmax, aspect='equal')
        axes[row, 0].set_title(f'{title} XY (iz={iz})')
        axes[row, 0].set_xlabel('x [Å]')
        axes[row, 0].set_ylabel('y [Å]')
        plt.colorbar(im, ax=axes[row, 0], fraction=0.03, pad=0.02)
        
        # XZ slice
        xz_slice = field[ix, :, :].T
        vmax = np.max(np.abs(xz_slice))
        y_min = origin[1]
        y_max = origin[1] + field.shape[1] * step
        z_min = origin[2]
        z_max = origin[2] + field.shape[2] * step
        im = axes[row, 1].imshow(xz_slice, origin='lower', extent=[y_min, y_max, z_min, z_max], 
                                 cmap='seismic', vmin=-vmax, vmax=vmax, aspect='equal')
        axes[row, 1].set_title(f'{title} XZ (ix={ix})')
        axes[row, 1].set_xlabel('y [Å]')
        axes[row, 1].set_ylabel('z [Å]')
        plt.colorbar(im, ax=axes[row, 1], fraction=0.03, pad=0.02)
        # Mark probe heights with gray dotted lines
        for h in heights:
            axes[row, 1].axhline(y=h, color='gray', linestyle=':', alpha=0.7, linewidth=1)
        
        # YZ slice
        yz_slice = field[:, iy, :].T
        vmax = np.max(np.abs(yz_slice))
        x_min = origin[0]
        x_max = origin[0] + field.shape[0] * step
        im = axes[row, 2].imshow(yz_slice, origin='lower', extent=[x_min, x_max, z_min, z_max], 
                                 cmap='seismic', vmin=-vmax, vmax=vmax, aspect='equal')
        axes[row, 2].set_title(f'{title} YZ (iy={iy})')
        axes[row, 2].set_xlabel('x [Å]')
        axes[row, 2].set_ylabel('z [Å]')
        plt.colorbar(im, ax=axes[row, 2], fraction=0.03, pad=0.02)
        # Mark probe heights with gray dotted lines
        for h in heights:
            axes[row, 2].axhline(y=h, color='gray', linestyle=':', alpha=0.7, linewidth=1)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'diagnostic_panel_slices.png'), dpi=120, bbox_inches='tight')
    plt.close()
    print(f"Saved diagnostic panel slices: {os.path.join(output_dir, 'diagnostic_panel_slices.png')}")


# ═══════════════════════════════════════════════════════════════════════════════
# Pauli Parameter Fitting (modular, reusable, portable)
# ═══════════════════════════════════════════════════════════════════════════════

def _fit_pauli_powerlaw(z, overlap_raw, e_ref, z_min=2.0, z_max=3.5):
    """Fit Pauli power-law model: E_DFTB(z) = A * overlap(z)^beta."""
    from scipy.optimize import curve_fit
    
    mask = (z >= z_min) & (z <= z_max)
    if mask.sum() < 3:
        raise ValueError(f"Need >=3 points in fit range [{z_min},{z_max}]")
    z_fit = z[mask]
    o_fit = overlap_raw[mask]
    e_fit = e_ref[mask]
    pos_mask = (o_fit > 1e-15) & (e_fit > 1e-15)
    if pos_mask.sum() < 3:
        raise ValueError("Not enough positive points")
    log_o = np.log(o_fit[pos_mask])
    log_e = np.log(e_fit[pos_mask])
    beta_ll, lnA_ll = np.polyfit(log_o, log_e, 1)
    A_ll = np.exp(lnA_ll)
    def model(overlap, A, beta):
        return A * overlap**beta
    try:
        popt, _ = curve_fit(model, o_fit, e_fit, p0=[A_ll, beta_ll],
                            bounds=([0.0, 0.0], [1e6, 5.0]))
        A_nls, beta_nls = popt
        e_pred = model(o_fit, A_nls, beta_nls)
        ss_res = np.sum((e_fit - e_pred)**2)
        ss_tot = np.sum((e_fit - np.mean(e_fit))**2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
    except Exception as e:
        print(f"  WARNING: Nonlinear fit failed ({e}), using log-linear")
        A_nls, beta_nls = A_ll, beta_ll
        e_pred = model(o_fit, A_nls, beta_nls)
        r2 = 0.0
    return A_nls, beta_nls, r2, e_pred


def _load_fdbm_grids(fdbm_dir):
    """Load FDBM forcefield grids from directory.
    
    Loads overlap_raw (raw Pauli overlap, A=1 beta=1) for fitting,
    plus E_Pauli_field, E_ES_field, E_vdw_field for diagnostics.
    """
    paths = {
        'overlap_raw': os.path.join(fdbm_dir, 'overlap_raw.npy'),
        'pauli': os.path.join(fdbm_dir, 'E_Pauli_field.npy'),
        'es':    os.path.join(fdbm_dir, 'E_ES_field.npy'),
        'vdw':   os.path.join(fdbm_dir, 'E_vdw_field.npy'),
    }
    grids = {}
    for key, path in paths.items():
        grids[key] = np.load(path) if os.path.exists(path) else None
    return grids


def _load_dftb_zscan(zscan_dir):
    """Load DFTB z-scan reference data."""
    z_path = os.path.join(zscan_dir, 'zscan_z.npy')
    e_path = os.path.join(zscan_dir, 'zscan_energy_eV.npy')
    if not (os.path.exists(z_path) and os.path.exists(e_path)):
        return None, None
    z = np.load(z_path)
    e = np.load(e_path)
    return z, e - e[-1]  # Relative energy


def _plot_pauli_fit(z, e_ref, e_fitted, A, beta, fname, title, z_min=2.0, z_max=3.5):
    """Plot per-atom Pauli fit (linear + log)."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    mask = (z >= z_min) & (z <= z_max)
    # Linear
    ax = axes[0]
    ax.plot(z, e_ref, 'o-', color='tab:blue', markersize=3, label='DFTB Ref', zorder=3)
    ax.plot(z[mask], e_fitted[mask], 's--', color='tab:red', markersize=3, label=f'Fit A={A:.0f} b={beta:.3f}', zorder=2)
    ax.axvspan(z_min, z_max, alpha=0.08, color='gray', label='Fit range')
    ax.set_xlabel('z [Å]')
    ax.set_ylabel('Energy [eV]')
    ax.set_title(f'{title} (Linear)')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    # Log
    ax = axes[1]
    pos = e_ref > 1e-12
    ax.semilogy(z[pos], e_ref[pos], 'o-', color='tab:blue', markersize=3, label='DFTB Ref')
    ax.semilogy(z[mask], e_fitted[mask], 's--', color='tab:red', markersize=3, label='Fit')
    ax.axvspan(z_min, z_max, alpha=0.08, color='gray')
    ax.set_xlabel('z [Å]')
    ax.set_ylabel('Energy [eV]')
    ax.set_title(f'{title} (Log)')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(fname, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {fname}")


def _plot_fitting_summary(all_results, fname, basis, z_min, z_max):
    """Plot summary comparing all atoms."""
    import matplotlib.pyplot as plt
    n_atoms = len(all_results)
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    # Panel 1: A_pauli per atom
    ax = axes[0, 0]
    idxs = [r['idx'] for r in all_results]
    As = [r['A'] for r in all_results]
    ax.bar(idxs, As, color='tab:orange')
    ax.set_xlabel('Atom index')
    ax.set_ylabel('A_pauli')
    ax.set_title(f'A_pauli per atom ({basis})')
    ax.grid(True, alpha=0.3, axis='y')
    # Panel 2: beta per atom
    ax = axes[0, 1]
    betas = [r['beta'] for r in all_results]
    ax.bar(idxs, betas, color='tab:green')
    ax.set_xlabel('Atom index')
    ax.set_ylabel('beta_pauli')
    ax.set_title(f'beta_pauli per atom ({basis})')
    ax.grid(True, alpha=0.3, axis='y')
    # Panel 3: RMSE per atom
    ax = axes[1, 0]
    rmses = [r['rmse_fit'] for r in all_results]
    ax.bar(idxs, rmses, color='tab:red')
    ax.set_xlabel('Atom index')
    ax.set_ylabel('RMSE fit [eV]')
    ax.set_title(f'RMSE(fit range) per atom ({basis})')
    ax.grid(True, alpha=0.3, axis='y')
    # Panel 4: All fitted curves overlaid
    ax = axes[1, 1]
    for r in all_results:
        z = r['z']
        e_fit = r['e_fitted']
        ax.plot(z, e_fit, '-', lw=1.0, label=f"atom {r['idx']}")
    ax.set_xlabel('z [Å]')
    ax.set_ylabel('Fitted Pauli [eV]')
    ax.set_title(f'Fitted Pauli curves ({basis})')
    ax.set_yscale('log')
    ax.legend(fontsize=7, ncol=2)
    ax.grid(True, alpha=0.3)
    ax.axvspan(z_min, z_max, alpha=0.08, color='gray')
    plt.suptitle(f'Multi-Atom Summary: {basis}', fontsize=12)
    plt.tight_layout()
    plt.savefig(fname, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved summary: {fname}")


def _run_dftb_zscan_for_atom(target_idx, mol_names, mol_pos, tip_names, sk_prefix, 
                             z_distances, out_dir, xyz_path, tip_path):
    """Run DFTB z-scan for a single target atom. Returns z_vals, e_vals."""
    import time
    from pyBall import dftb_utils as du
    from pyBall import atomicUtils as au
    
    HAU2EV = 27.211386245988
    target_name = mol_names[target_idx]
    target_pos = mol_pos[target_idx]
    atom_dir = os.path.join(out_dir, f'atom_{target_idx}')
    os.makedirs(atom_dir, exist_ok=True)

    print(f"\n{'='*60}")
    print(f"Target atom {target_idx}: {target_name} at [{target_pos[0]:.4f}, {target_pos[1]:.4f}, {target_pos[2]:.4f}]")
    print(f"Output: {atom_dir}")
    print(f"{'='*60}")

    cache_path = os.path.join(atom_dir, 'zscan_results_cache.npz')
    results = []
    if os.path.exists(cache_path):
        print(f"Loading cache: {cache_path}")
        cache = np.load(cache_path, allow_pickle=True)
        cached_z = cache['z_distances']
        if np.allclose(cached_z, z_distances):
            results = cache['results'].tolist()
            print(f"  Using {len(results)} cached results")
        else:
            print("  Cache z-range mismatch, recomputing")

    if len(results) != len(z_distances):
        combined_names = list(mol_names) + list(tip_names)
        for iz, z in enumerate(z_distances):
            print(f"\n[z-scan {iz+1}/{len(z_distances)}] z = {z:.2f} Å")
            o_pos = np.array([target_pos[0], target_pos[1], target_pos[2] + z])
            c_pos = np.array([target_pos[0], target_pos[1], target_pos[2] + z + 1.13])
            co_pos_shifted = np.array([o_pos, c_pos])
            combined_pos = np.vstack([mol_pos, co_pos_shifted])

            work_dir = os.path.join(atom_dir, f'zscan_z{z:.2f}')
            t_start = time.time()
            try:
                energy_ha = du.run_dftb_sp(work_dir, combined_names, combined_pos, sk_prefix)
                energy_ev = energy_ha * HAU2EV
                t_elapsed = time.time() - t_start
                print(f"  Energy: {energy_ha:.8f} Ha = {energy_ev:.6f} eV  ({t_elapsed:.1f}s)")
                results.append({'z': float(z), 'energy_Ha': float(energy_ha), 'energy_eV': float(energy_ev)})
            except Exception as e:
                print(f"  ERROR: {e}")
                raise

        np.savez(cache_path, results=results, z_distances=z_distances)
        print(f"\nSaved cache to {cache_path}")

    z_vals = np.array([r['z'] for r in results])
    e_vals = np.array([r['energy_eV'] for r in results])
    np.save(os.path.join(atom_dir, 'zscan_z.npy'), z_vals)
    np.save(os.path.join(atom_dir, 'zscan_energy_eV.npy'), e_vals)

    with open(os.path.join(atom_dir, 'zscan_results.txt'), 'w') as f:
        f.write("DFTB Z-Scan Results\n")
        f.write("="*70 + "\n")
        f.write(f"Target: {target_name} [{target_pos[0]:.4f}, {target_pos[1]:.4f}, {target_pos[2]:.4f}]\n")
        f.write(f"CO bond: 1.13 Å\n\n")
        f.write(f"{'z [Å]':>10}  {'E [Ha]':>16}  {'E [eV]':>16}\n")
        f.write("-"*70 + "\n")
        for r in results:
            f.write(f"{r['z']:10.2f}  {r['energy_Ha']:16.8f}  {r['energy_eV']:16.6f}\n")

    e_rel = e_vals - e_vals[-1]
    print(f"\nAtom {target_idx} complete: {len(z_vals)} points")
    print(f"  Rel energy at contact: {e_rel[0]:.4f} eV")
    return z_vals, e_vals


def fit_pauli_parameters(xyz_file, basis='mio-1-1', target_indices=[0], 
                         fdbm_dir=None, zscan_dir=None, output_dir='fit_pauli',
                         z_min=2.0, z_max=3.0, generate_ref=False,
                         step=0.1, margin=4.0, z_extra=6.0,
                         sk_prefix=None, tip_xyz='CO.xyz',
                         scan_range=3.0, scan_step=0.1, height_range=[2.8, 3.6], height_step=0.1):
    """High-level modular function to fit Pauli parameters against DFTB reference.
    
    This function integrates the full fitting workflow:
    1. If fdbm_dir is None/missing: run FDBM pipeline with new CO tip handling
    2. If zscan_dir is None/missing and generate_ref=True: run DFTB z-scan
    3. Load FDBM grids and DFTB z-scan data
    4. For each target atom: extract profile, fit power-law, save results
    5. Generate summary plots and table
    
    Args:
        xyz_file: Path to molecule XYZ file
        basis: DFTB+ basis set (mio-1-1, 3ob-3-1, etc.)
        target_indices: List of atom indices to fit (e.g., [0, 1, 20, 21])
        fdbm_dir: Pre-computed FDBM grid directory (if None, generates on-the-fly)
        zscan_dir: Pre-computed DFTB z-scan directory (if None and generate_ref=True, generates)
        output_dir: Output directory for fitting results
        z_min, z_max: Fit range in Å (contact region)
        generate_ref: Whether to generate DFTB z-scan if missing
        step, margin, z_extra: Grid parameters for FDBM generation
        sk_prefix: DFTB+ Slater-Koster path (if None, uses default from dftb_utils)
        tip_xyz: Tip molecule XYZ file (default: CO.xyz)
        scan_range, scan_step: AFM scan parameters (for FDBM grid generation)
        height_range, height_step: AFM height parameters (for FDBM grid generation)
    
    Returns:
        dict: Fitting results with keys:
            - 'basis': basis set name
            - 'atoms': list of per-atom results (dict with A, beta, rmse, etc.)
            - 'A_mean', 'beta_mean': mean values across atoms
            - 'A_std', 'beta_std': standard deviations
    """
    import json
    import time
    from pyBall import atomicUtils as au
    from pyBall.DFTB import TestUtils as tu
    from pyBall import dftb_utils as du
    
    A_PAULI_DEFAULT = 16.0
    t0 = time.time()
    os.makedirs(output_dir, exist_ok=True)
    
    # Load molecule
    mol_pos, _, mol_names, _, _ = au.load_xyz(xyz_file)
    mol_pos = np.array(mol_pos, dtype=np.float64)
    for idx in target_indices:
        if idx < 0 or idx >= len(mol_names):
            raise ValueError(f"Target index {idx} out of range (0-{len(mol_names)-1})")
    
    # Set up SK path
    if sk_prefix is None:
        if basis in du.SK_PATHS:
            sk_prefix = du.SK_PATHS[basis]
        else:
            raise ValueError(f"Basis '{basis}' not found in dftb_utils.SK_PATHS; provide sk_prefix explicitly")
    
    # Step 1: Generate FDBM grids if needed
    if fdbm_dir is None or not os.path.isdir(fdbm_dir):
        print(f"\nGenerating FDBM grids for {basis}...")
        fdbm_dir = os.path.join(output_dir, f'fdbm_grids_{basis.replace("-", "_")}')
        os.makedirs(fdbm_dir, exist_ok=True)
        
        run_afm_from_xyz(
            xyz_file, output_dir=fdbm_dir, basis=basis,
            step=step, margin=margin, z_extra=z_extra,
            scan_range=scan_range, scan_step=scan_step,
            height_range=height_range, height_step=height_step,
            co_tip_dir=None,  # Force on-the-fly CO computation with new padding/rolling
            plot_steps=False
        )
        print(f"  FDBM grids saved to: {fdbm_dir}")
    
    # Load FDBM grids
    grids = _load_fdbm_grids(fdbm_dir)
    if grids['overlap_raw'] is None:
        raise FileNotFoundError(f"overlap_raw.npy not found in {fdbm_dir}. Run pipeline first (it saves raw overlap).")
    
    # Read grid spec
    log_path_grid = os.path.join(fdbm_dir, 'step1_density', 'log.txt')
    grid_spec_path = os.path.join(fdbm_dir, 'grid_spec.txt')
    origin, ngrid, step_grid = afm.read_grid_spec_from_log(log_path_grid)
    
    # Try grid_spec.txt (new format)
    if origin is None and os.path.exists(grid_spec_path):
        with open(grid_spec_path, 'r') as f:
            lines = f.readlines()
            for line in lines:
                if line.startswith('origin ='):
                    origin = np.array(eval(line.split('=')[1].strip()))
                elif line.startswith('ngrid ='):
                    ngrid = np.array(eval(line.split('=')[1].strip()))
                elif line.startswith('step ='):
                    step_grid = float(line.split('=')[1].strip())
    
    if origin is None:
        # Fallback: compute from grid shape and molecule
        pauli_shape = grids['pauli'].shape
        ngrid = pauli_shape
        step_grid = step  # Use the step parameter
        # Better estimate: center grid on molecule
        mol_center = mol_pos.mean(axis=0)
        grid_size = np.array(ngrid) * step_grid
        origin = mol_center - 0.5 * grid_size
        print(f"  WARNING: Could not read grid spec, estimated from molecule center")
    print(f"  Grid: origin={origin.round(2)} ngrid={ngrid} step={step_grid}")
    
    # Step 2: Generate DFTB z-scan if needed
    zscan_dir_base = zscan_dir if zscan_dir else os.path.join(output_dir, f'zscan_{basis.replace("-", "_")}')
    
    if generate_ref:
        print(f"\nGenerating DFTB z-scan reference for {basis}...")
        tip_path = os.path.join(os.path.dirname(xyz_file), tip_xyz)
        tip_pos, _, tip_names, _, _ = au.load_xyz(tip_path)
        
        z_distances = np.arange(2.0, 10.0 + 0.15*0.5, 0.15)
        print(f"  Z-scan: {len(z_distances)} points from {z_distances.min():.2f} to {z_distances.max():.2f} Å")
        
        for target_idx in target_indices:
            _run_dftb_zscan_for_atom(
                target_idx, mol_names, mol_pos, tip_names, sk_prefix,
                z_distances, zscan_dir_base, xyz_file, tip_path
            )
        print(f"  DFTB z-scan saved to: {zscan_dir_base}")
    
    # Step 3: Fit each atom
    all_results = []
    for target_idx in target_indices:
        target_name = mol_names[target_idx]
        target_pos = mol_pos[target_idx]
        atom_out_dir = os.path.join(output_dir, f'atom_{target_idx}')
        zscan_atom_dir = os.path.join(zscan_dir_base, f'atom_{target_idx}')
        os.makedirs(atom_out_dir, exist_ok=True)
        
        print(f"\nAtom {target_idx} ({target_name}) at [{target_pos[0]:.4f}, {target_pos[1]:.4f}, {target_pos[2]:.4f}]")
        
        # Load DFTB z-scan
        z_ref, e_ref = _load_dftb_zscan(zscan_atom_dir)
        if z_ref is None:
            raise FileNotFoundError(f"No z-scan found in {zscan_atom_dir}. Generate with --fit_generate_ref or provide valid --fit_zscan_dir")
        
        # Extract raw overlap profile at atom XY, at absolute z = target_pos[2] + z_ref
        # z_ref are distances ABOVE the atom (e.g. 2.0..10.0 Å)
        # extract_z_profile uses: z_abs = atom_pos[2] + z_distances
        overlap_col = tu.extract_z_profile(grids['overlap_raw'], target_pos, origin, step_grid, z_distances=z_ref)
        if overlap_col is None:
            raise ValueError(f"overlap_raw extraction failed for atom {target_idx}")
        overlap_safe = np.clip(overlap_col, 1e-30, None)
        
        # Diagnostics: verify z-grid alignment
        grid_z_range = [origin[2], origin[2] + grids['overlap_raw'].shape[2] * step_grid]
        print(f"  Grid z: [{grid_z_range[0]:.2f}, {grid_z_range[1]:.2f}] Å")
        print(f"  z_ref: [{z_ref.min():.2f}, {z_ref.max():.2f}] Å (above atom at z={target_pos[2]:.2f})")
        print(f"  z_abs at contact: {target_pos[2] + z_ref.min():.2f} Å")
        print(f"  overlap at z_ref[0]={z_ref[0]:.2f}: {overlap_col[0]:.4e}")
        idx_z3 = np.argmin(np.abs(z_ref - 3.0))
        print(f"  overlap at z_ref~3.0 (={z_ref[idx_z3]:.2f}): {overlap_col[idx_z3]:.4e}")
        print(f"  e_ref at z_ref~3.0: {e_ref[idx_z3]:.4e} eV")
        
        # Fit: E_DFTB = A_fit * overlap^beta_fit
        A_fit, beta_fit, r2, e_fitted_range = _fit_pauli_powerlaw(
            z_ref, overlap_safe, e_ref, z_min=z_min, z_max=z_max
        )
        
        e_fitted_all = A_fit * overlap_safe**beta_fit
        mask_fit = (z_ref >= z_min) & (z_ref <= z_max)
        rmse_fit = np.sqrt(np.mean((e_ref[mask_fit] - e_fitted_range)**2))
        rmse_all = np.sqrt(np.mean((e_ref - e_fitted_all)**2))
        
        print(f"  Fit: A={A_fit:.4f}, beta={beta_fit:.4f}, R2={r2:.6f}, RMSE(fit)={rmse_fit:.4e} eV")
        print(f"  Consistency check at z_ref~3.0: E_fit={e_fitted_all[idx_z3]:.4e} vs E_DFTB={e_ref[idx_z3]:.4e} eV")
        
        # Save
        params = {
            'basis': basis, 'atom_idx': target_idx, 'atom_name': target_name,
            'A_pauli': float(A_fit), 'beta_pauli': float(beta_fit),
            'R2_fit': float(r2), 'RMSE_fit': float(rmse_fit), 'RMSE_all': float(rmse_all),
            'fit_z_min': z_min, 'fit_z_max': z_max,
        }
        with open(os.path.join(atom_out_dir, 'params.json'), 'w') as f:
            json.dump(params, f, indent=2)
        
        np.save(os.path.join(atom_out_dir, 'z_ref.npy'), z_ref)
        np.save(os.path.join(atom_out_dir, 'e_ref.npy'), e_ref)
        np.save(os.path.join(atom_out_dir, 'overlap_col.npy'), overlap_col)
        np.save(os.path.join(atom_out_dir, 'e_fitted.npy'), e_fitted_all)
        
        # Plot
        _plot_pauli_fit(
            z_ref, e_ref, e_fitted_all, A_fit, beta_fit,
            fname=os.path.join(atom_out_dir, 'fit_pauli.png'),
            title=f'{target_name}{target_idx} ({basis})',
            z_min=z_min, z_max=z_max
        )
        
        all_results.append({
            'idx': target_idx, 'name': target_name, 'pos': target_pos,
            'A': A_fit, 'beta': beta_fit, 'r2': r2,
            'rmse_fit': rmse_fit, 'rmse_all': rmse_all,
            'z': z_ref, 'e_fitted': e_fitted_all, 'e_ref': e_ref,
        })
    
    # Step 4: Summary
    if len(all_results) > 1:
        _plot_fitting_summary(all_results, os.path.join(output_dir, 'summary_all_atoms.png'), basis, z_min, z_max)
    
    # Write summary table
    with open(os.path.join(output_dir, 'summary.txt'), 'w') as f:
        f.write("FDBM Pauli Fitting Summary\n")
        f.write("="*70 + "\n")
        f.write(f"Basis: {basis}\n")
        f.write(f"Atoms: {[r['idx'] for r in all_results]}\n")
        f.write(f"Fit range: z=[{z_min}, {z_max}] Å\n\n")
        f.write(f"{'Atom':>6} {'Name':>4} {'A_pauli':>10} {'beta':>8} {'R2':>10} {'RMSE_fit':>10} {'RMSE_all':>10}\n")
        f.write("-"*70 + "\n")
        for r in all_results:
            f.write(f"{r['idx']:6d} {r['name']:>4} {r['A']:10.2f} {r['beta']:8.4f} {r['r2']:10.6f} {r['rmse_fit']:10.4f} {r['rmse_all']:10.4f}\n")
        
        if len(all_results) > 1:
            As = [r['A'] for r in all_results]
            betas = [r['beta'] for r in all_results]
            f.write(f"\nMean ± std:\n")
            f.write(f"  A_pauli: {np.mean(As):.2f} ± {np.std(As):.2f}\n")
            f.write(f"  beta:    {np.mean(betas):.4f} ± {np.std(betas):.4f}\n")
        f.write(f"\nTime: {time.time()-t0:.1f}s\n")
    
    print(f"\nAll results saved to: {output_dir}/")
    
    # Return structured results
    result_dict = {
        'basis': basis,
        'atoms': all_results,
        'A_mean': np.mean([r['A'] for r in all_results]) if all_results else None,
        'beta_mean': np.mean([r['beta'] for r in all_results]) if all_results else None,
        'A_std': np.std([r['A'] for r in all_results]) if all_results else None,
        'beta_std': np.std([r['beta'] for r in all_results]) if all_results else None,
    }
    return result_dict
