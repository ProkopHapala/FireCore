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
                                grid_spec=None, step=0.15, margin=4.0, z_extra=6.0, verbosity=0):
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
    rho_scf, rho_na, rho_diff = _project_densities(geo, evecs, basis, grid_spec, verbosity)
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
    pauli_params={'A': 16.0, 'beta': 1.0},
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
    
    # Step 3: Pauli repulsion
    print("\nStep 3: Computing Pauli repulsion...")
    E_pauli_field, grads_pauli = afm.compute_pauli_field(
        rho_grid, rho_tip_total, step,
        A_pauli=pauli_params['A'], beta_pauli=pauli_params['beta']
    )
    if plot_steps:
        plot_step3_outputs(E_pauli_field, grads_pauli, output_dir, origin, step,
                          pauli_params['A'], pauli_params['beta'])
    np.save(os.path.join(output_dir, 'E_Pauli_field.npy'), E_pauli_field)
    np.save(os.path.join(output_dir, 'grads_E_Pauli.npy'), grads_pauli)
    
    # Step 4: Electrostatic convolution
    print("\nStep 4: Computing electrostatic convolution...")
    E_ES_field, grads_ES = afm.compute_es_conv_field(V_ES, rho_tip_delta, step)
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


def run_afm_from_xyz(
    xyz_file,
    output_dir,
    basis,
    slako_prefix='mio-1-1',
    co_tip_dir=None,
    work_dir=None,
    step=0.15, margin=4.0, z_extra=6.0,
    scan_range=3.0, scan_step=0.1,
    height_range=(3.0, 5.5), height_step=0.1,
    pauli_params={'A': 16.0, 'beta': 1.0},
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
        co_tip_dir: directory with co_rho_total.npy + co_rho_delta.npy (required)
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

    # Step 1: density from DFTB+
    print(f"\nStep 1: DFTB+ SCF + density projection...")
    d = get_density_from_dftb_plus(atomPos, atomTypes, basis, slako_prefix, work_dir,
                                    step=step, margin=margin, z_extra=z_extra)

    # Load CO tip densities (must be pre-computed)
    if co_tip_dir is None or not os.path.isdir(co_tip_dir):
        raise FileNotFoundError(f"CO tip directory not found: {co_tip_dir}\n"
                                 "Run compute_co_tip.py first to generate co_rho_total.npy + co_rho_delta.npy")
    co_rho_total = np.load(os.path.join(co_tip_dir, 'co_rho_total.npy'))
    co_rho_delta = np.load(os.path.join(co_tip_dir, 'co_rho_delta.npy'))
    print(f"  CO tip loaded: {co_rho_total.shape}")

    # Interpolate CO tip to match density grid if shapes differ
    from scipy.ndimage import zoom
    target_shape = d['ngrid']
    if co_rho_total.shape != tuple(target_shape):
        zoom_factors = [target_shape[i] / co_rho_total.shape[i] for i in range(3)]
        co_rho_total = zoom(co_rho_total, zoom_factors, order=1).astype(np.float32)
        co_rho_delta = zoom(co_rho_delta, zoom_factors, order=1).astype(np.float32)
        print(f"  CO tip interpolated to: {co_rho_total.shape}")

    return run_afm_pipeline(
        d['rho_scf'], d['rho_na'], d['rho_diff'], d['V_ES'],
        co_rho_total, co_rho_delta,
        atomPos, atomTypes,
        d['origin'], step, d['ngrid'],
        scan_xs, scan_ys, heights,
        output_dir,
        pauli_params=pauli_params, vdw_params=vdw_params,
        relax_params=relax_params, plot_steps=plot_steps
    )


def plot_diagnostic_panel(E_pauli, E_es, E_vdw, E_total, origin, step, heights, output_dir):
    """Plot diagnostic panel with 4 columns (Total, Pauli, Electrostatics, vdW) and n-rows for heights.

    Each subplot has symmetric vmin/vmax zero-centered with its own colorbar (seismic colormap).
    """
    import matplotlib.pyplot as plt
    n_heights = len(heights)
    fig, axes = plt.subplots(n_heights, 4, figsize=(16, 3*n_heights))
    if n_heights == 1:
        axes = axes.reshape(1, -1)

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
            ax.set_title(f'{title} z={z:.1f}Å')
            plt.colorbar(im, ax=ax, fraction=0.03, pad=0.02)

    plt.subplots_adjust(left=0.02, right=0.98, bottom=0.02, top=0.95, wspace=0.25, hspace=0.2)
    fname = os.path.join(output_dir, 'diagnostic_panel.png')
    plt.savefig(fname, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {fname}")
