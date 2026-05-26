"""Test FDBM linear fitting using DFT reference data from small_mols_NaCl_New.

This script adapts test_fdbm_fit_gridff_mock.py to use actual DFT data:
- Cube files: charge_density.cube, electrostatic_potential.cube
- DFT energies: 3-results/H2O-H/*.dat, 3-results/H2O-O/*.dat
- Molecule configs: from 1-inputs/confs/ or regenerate

Key changes from mock version:
1. Load cube files using ppafm.io.loadCUBE (handles Bohr→Å, Hartree→eV conversion)
2. Sample cube data at atom positions via trilinear interpolation
3. Load DFT reference energies from .dat files
4. Fit P_i coefficients (Pauli) with fixed Q_i (charges)

Outputs:
- tests/tMMFF/out_fdbm_dft/*png  (diagnostic plots)
- tests/tMMFF/out_fdbm_dft/*.npz (data + fitted coeffs)
"""

import os
import sys
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def _repo_root():
    here = os.path.dirname(os.path.abspath(__file__))
    return os.path.realpath(os.path.join(here, '..', '..'))


# Add repository root to path for pyBall imports
sys.path.insert(0, _repo_root())

# Unit conversion constants
bohr2ang = 0.5291772109217
hartree2ev = 27.211396132


def _load_xyz_with_charges(xyz_path):
    """Load XYZ file with optional charges in 5th column."""
    with open(xyz_path, 'r') as f:
        lines = f.readlines()
    n = int(lines[0].strip())
    names = []
    pos = []
    qs = []
    for i in range(2, 2 + n):
        w = lines[i].split()
        names.append(w[0])
        pos.append((float(w[1]), float(w[2]), float(w[3])))
        qs.append(float(w[4]) if len(w) > 4 else 0.0)
    return np.array(pos, dtype=np.float32), np.array(names), np.array(qs, dtype=np.float32)


def _rot_z(phi):
    c = float(np.cos(phi)); s = float(np.sin(phi))
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]], dtype=np.float32)


def _make_transforms_zscan(xy=(1.0, 1.0), z_vals=None, phis=None):
    """Generate z-scan transforms (for diagnostic comparison with DFT)."""
    if z_vals is None:
        z_vals = np.linspace(0.5, 6.0, 56)
    if phis is None:
        phis = np.array([0.0], dtype=np.float32)
    z_vals = np.asarray(z_vals, dtype=np.float32)
    phis = np.asarray(phis, dtype=np.float32)
    Ts = []
    labels = []
    for phi in phis:
        R = _rot_z(phi)
        for z in z_vals:
            T = np.zeros((3, 4), dtype=np.float32)
            T[:, :3] = R
            T[:, 3] = np.array([xy[0], xy[1], float(z)], dtype=np.float32)
            Ts.append(T)
            labels.append((float(z), float(phi)))
    return np.array(Ts, dtype=np.float32), np.array(labels, dtype=np.float32)


def _load_cube_to_gridff(density_cube, potential_cube):
    """Load DFT cube files and convert to FireCore GridFF format.
    
    Cube file format:
    - Line 1: comment
    - Line 2: comment
    - Line 3: n_atoms origin_x origin_y origin_z (Bohr)
    - Line 4: n_x spacing_x spacing_y spacing_z (Bohr)
    - Line 5: n_y spacing_x spacing_y spacing_z (Bohr)
    - Line 6: n_z spacing_x spacing_y spacing_z (Bohr)
    - Lines 7-6+n_atoms: atom lines
    - Remaining: volumetric data (n_x * n_y * n_z values)
    
    Handles:
    - Bohr → Angstrom conversion for grid spacing and origin
    - Hartree → eV conversion for data values (for potential)
    - Transposition to (nx, ny, nz) order
    
    Returns:
        gridff: (nx, ny, nz, 4) array
        g0: (3,) grid origin in Angstrom
        dg: (3,) grid spacing in Angstrom
        lvec: (4,3) lattice vectors
    """
    
    def load_single_cube(fname, convert_to_ev=False):
        """Load a single cube file."""
        with open(fname, 'r') as f:
            # Skip first 2 comment lines
            f.readline()
            f.readline()
            
            # Line 3: n_atoms origin
            line3 = f.readline().split()
            n_atoms = int(line3[0])
            origin = np.array([float(line3[1]), float(line3[2]), float(line3[3])]) * bohr2ang
            
            # Lines 4-6: grid dimensions and spacing
            line4 = f.readline().split()
            nx = int(line4[0])
            dx = float(line4[1]) * bohr2ang
            
            line5 = f.readline().split()
            ny = int(line5[0])
            dy = float(line5[2]) * bohr2ang
            
            line6 = f.readline().split()
            nz = int(line6[0])
            dz = float(line6[3]) * bohr2ang
            
            # Skip atom lines
            for _ in range(n_atoms):
                f.readline()
            
            # Read volumetric data
            # CRITICAL: Cube file stores data with Z varying fastest (inner loop).
            # C-order reshape to (nx, ny, nz) makes last index (Z) vary fastest.
            n_vals = nx * ny * nz
            data = np.fromfile(f, sep=' ', count=n_vals)
            data = data.reshape((nx, ny, nz))  # C-order: z fastest
            
            # Convert units:
            # - Potential: Hartree/e → eV/e (multiply by Hartree2eV)
            # - Density: e/Bohr³ → e/Å³ (1 Bohr = bohr2ang Å, so 1/Bohr³ = 1/bohr2ang³ /Å³)
            if convert_to_ev:
                data *= hartree2ev  # Potential: Hartree → eV
            else:
                data /= (bohr2ang ** 3)  # Density: e/Bohr³ → e/Å³
            
            return data, origin, (nx, ny, nz), (dx, dy, dz)
    
    # Load density: convert e/Bohr³ → e/Å³
    rho, origin_rho, nDim_rho, dg_rho = load_single_cube(density_cube, convert_to_ev=False)
    
    # Load potential: convert Hartree/e → eV/e
    phi, origin_phi, nDim_phi, dg_phi = load_single_cube(potential_cube, convert_to_ev=True)
    
    # Both should have same shape
    nx, ny, nz = nDim_rho
    
    # Create 4-channel array [nx, ny, nz, 4]
    # GridFF format: (nx, ny, nz, nch) with array[x, y, z, ch]
    gridff = np.zeros((nx, ny, nz, 4), dtype=np.float32)
    
    # Channel 0: Pauli proxy = electron density (rho)
    # rho is already (nx, ny, nz) from C-order reshape
    gridff[:,:,:,0] = rho
    
    # Channel 1: London (zero for now)
    gridff[:,:,:,1] = 0.0
    
    # Channel 2: Coulomb = electrostatic potential (in eV)
    # phi is already (nx, ny, nz) from C-order reshape
    gridff[:,:,:,2] = phi
    
    # Channel 3: Hbond (zero for NaCl)
    gridff[:,:,:,3] = 0.0
    
    # Grid parameters
    g0 = origin_rho  # Origin in Angstrom
    dg = dg_rho  # Spacing in Angstrom
    
    # Build lattice vectors
    lvec = np.zeros((4, 3), dtype=np.float64)
    lvec[0] = g0
    lvec[1] = [nx * dg[0], 0.0, 0.0]
    lvec[2] = [0.0, ny * dg[1], 0.0]
    lvec[3] = [0.0, 0.0, nz * dg[2]]
    
    print(f'Cube grid: shape=({nx}, {ny}, {nz}) g0={g0} dg={dg}')
    
    return gridff, g0, dg, lvec


def _sample_cube_at_positions(grid_data, positions, g0, dg):
    """Sample cube data at arbitrary positions using trilinear interpolation.
    
    Args:
        grid_data: (nx, ny, nz, nch) array
        positions: (n, 3) positions in Angstrom
        g0: (3,) grid origin
        dg: (3,) grid spacing
    
    Returns:
        sampled: (n, nch) array of sampled values
    """
    from scipy.ndimage import map_coordinates
    
    nx, ny, nz, nch = grid_data.shape
    dx, dy, dz = dg
    
    # Convert positions to grid indices (fractional)
    # idx = (pos - g0) / dg
    idx = (positions - g0) / np.array(dg)
    
    # Apply PBC: wrap to [0, nx), [0, ny), [0, nz)
    idx[:, 0] = idx[:, 0] % nx
    idx[:, 1] = idx[:, 1] % ny
    idx[:, 2] = idx[:, 2] % nz
    
    sampled = np.zeros((len(positions), nch), dtype=np.float32)
    
    # Sample each channel
    for ch in range(nch):
        # grid_data is (nx, ny, nz) in C-order, map_coordinates expects coords matching array dims
        channel_data = grid_data[:, :, :, ch]  # (nx, ny, nz)
        
        # coords correspond to (x, y, z) axes, shape (npoints, 3)
        coords = idx[:, [0, 1, 2]]
        
        # Sample with order=1 (linear interpolation)
        sampled[:, ch] = map_coordinates(channel_data, coords, order=1, mode='wrap')
    
    return sampled


def _sample_cube_rigid_batch(grid_data, transforms, apos0, g0, dg):
    """Sample cube data for all rigid transforms.
    
    Args:
        grid_data: (nx, ny, nz, nch) array
        transforms: (nconf, 3, 4) rigid transforms
        apos0: (natoms, 3) reference atom positions
        g0: (3,) grid origin
        dg: (3,) grid spacing
    
    Returns:
        Es: (nconf, natoms, nch) sampled energies per atom per channel
    """
    nconf = transforms.shape[0]
    natoms = apos0.shape[0]
    nch = grid_data.shape[3]
    
    Es = np.zeros((nconf, natoms, nch), dtype=np.float32)
    
    for i, T in enumerate(transforms):
        R = T[:, :3]
        t = T[:, 3]
        # Transform atom positions: pos = apos0 @ R.T + t
        pos = (apos0 @ R.T) + t
        # Sample at these positions
        Es[i] = _sample_cube_at_positions(grid_data, pos, g0, dg)
    
    return Es


def _load_dft_results(dft_results_dir, molecule_site='H2O-H'):
    """Load DFT energy results from 3-results directory.
    
    Args:
        dft_results_dir: Path to 3-results directory
        molecule_site: e.g., 'H2O-H' or 'H2O-O'
    
    Returns:
        results: dict mapping config keys to (z_vals, E_vals) arrays
    """
    site_dir = os.path.join(dft_results_dir, molecule_site)
    if not os.path.exists(site_dir):
        raise FileNotFoundError(f'DFT results directory not found: {site_dir}')
    
    results = {}
    
    # Parse all .dat files
    for fname in os.listdir(site_dir):
        if not fname.endswith('.dat'):
            continue
        
        # Parse filename: ixX_iyY-orientxy_XY-orientz_Z.dat
        key = fname.replace('.dat', '')
        fpath = os.path.join(site_dir, fname)
        
        # Load data (2 columns: z [Å], E [eV])
        data = np.loadtxt(fpath)
        if data.ndim == 1:
            data = data.reshape(-1, 1)
        
        # The files contain 3 repeats of the same z-scan (33 points × 3 = 99 lines)
        # Average the repeats
        n_unique = 33
        if len(data) == n_unique * 3:
            z_vals = data[:n_unique, 0]
            E_vals = data[:n_unique, 1]  # Take first repeat (or could average)
        else:
            z_vals = data[:, 0]
            E_vals = data[:, 1]
        
        results[key] = (z_vals, E_vals)
    
    print(f'Loaded {len(results)} DFT result files from {molecule_site}')
    return results


def _fit_linear(A, b):
    """Linear least squares fit."""
    x, *_ = np.linalg.lstsq(A, b, rcond=None)
    return x


def main():
    import argparse
    parser = argparse.ArgumentParser(description='FDBM DFT-fit test with cube files')
    parser.add_argument('--no-diagnostics', action='store_true', help='Skip diagnostic plots')
    parser.add_argument('--zmin-offset', type=float, default=2.0, help='Z offset for XZ color scale')
    parser.add_argument('--z-ylim', type=float, nargs=2, default=(-0.5, 0.5), help='Z-profile y-limits')
    parser.add_argument('--xy-range', type=float, nargs=2, default=(-2.0, 2.0), help='XY range for configs')
    parser.add_argument('--z-range', type=float, nargs=2, default=(1.4, 5.0), help='Z range (matches DFT)')
    parser.add_argument('--molecule-site', type=str, default='H2O-H', help='Molecule-site (H2O-H or H2O-O)')
    args = parser.parse_args()

    ROOT = _repo_root()
    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out_fdbm_dft')
    os.makedirs(out_dir, exist_ok=True)

    from pyBall.OCL import Surface_utils as su

    # Paths
    dft_data_dir = os.path.join(ROOT, 'tests', 'tSurf', 'small_mols_NaCl_New')
    density_cube = os.path.join(dft_data_dir, '4-elec_pot_chr_den', 'charge_density.cube')
    potential_cube = os.path.join(dft_data_dir, '4-elec_pot_chr_den', 'electrostatic_potential.cube')
    mol_path = os.path.join(ROOT, 'cpp', 'common_resources', 'xyz', 'H2O.xyz')
    sub_path = os.path.join(dft_data_dir, '0-geoms', 'NaCl.xyz')
    dft_results_dir = os.path.join(dft_data_dir, '3-results')

    print('=== FDBM DFT-fit test (cube files) ===')
    print('density_cube =', density_cube)
    print('potential_cube =', potential_cube)
    print('mol_path =', mol_path)
    print('sub_path =', sub_path)
    print('dft_results_dir =', dft_results_dir)
    print('molecule_site =', args.molecule_site)
    print('out_dir =', out_dir)

    # Load molecule
    apos, names, charges = _load_xyz_with_charges(mol_path)
    # Rigid body frame: center molecule
    apos0 = (apos - apos.mean(axis=0)).astype(np.float32)

    # Type mapping: per-element type (H,O)
    uniq = []
    type_ids = np.zeros(len(names), dtype=np.int32)
    for i, e in enumerate(names.tolist()):
        if e not in uniq:
            uniq.append(e)
        type_ids[i] = uniq.index(e)
    ntypes = len(uniq)
    print('types=', uniq, 'type_ids=', type_ids.tolist())
    print('charges=', charges.tolist())

    # Load GridFF from file (generated by gen_gridff_nacl_gpu.py)
    from pyBall.OCL.Surface_utils import load_bspline_gridff
    gridff_path = os.path.join(os.path.dirname(os.path.dirname(sub_path)), 'Bspline_PLQd.npy')
    gridff, meta = load_bspline_gridff(gridff_path)
    g0 = np.array(meta['g0'])
    dg = np.array(meta['dg'])
    ns = meta['ns']
    nx, ny, nz = ns
    # Build lvec from dg and ns if not in metadata
    if 'lvec' in meta:
        lvec = np.array(meta['lvec'])
    else:
        lvec = np.array([[nx * dg[0], 0.0, 0.0], [0.0, ny * dg[1], 0.0], [0.0, 0.0, nz * dg[2]]])
    print(f'Loaded GridFF from: {gridff_path}')
    print(f'  Grid shape: ({nx}, {ny}, {nz})')
    print(f'  g0: {g0}')
    print(f'  dg: {dg}')

    # Load substrate for visualization
    sub_data = su.load_substrate_xyz_with_lvec(sub_path)
    sub_apos = sub_data['apos']
    sub_enames = sub_data['enames']
    lvec_sub = sub_data['lvec']

    # Generate transforms (use z-scan matching DFT range)
    z_vals = np.linspace(args.z_range[0], args.z_range[1], 33)  # 33 points like DFT
    T_z, zphi = _make_transforms_zscan(xy=(0.0, 0.0), z_vals=z_vals, phis=[0.0])
    print(f'Generated {len(T_z)} z-scan transforms')

    # Sample GridFF at atom positions for all configs
    print('Sampling GridFF data at atom positions...')
    Es = _sample_cube_rigid_batch(gridff, T_z, apos0, g0, dg)
    print(f'Sampled Es shape: {Es.shape} (nconf={Es.shape[0]}, natoms={Es.shape[1]}, nch={Es.shape[2]})')

    # Load DFT reference energies
    print('Loading DFT reference energies...')
    dft_results = _load_dft_results(dft_results_dir, args.molecule_site)
    
    # Use config matching our xy=(0,0) scan: ix0_iy0-orientxy_Na-orientz_0
    sample_key = 'ix0_iy0-orientxy_Na-orientz_0'
    if sample_key not in dft_results:
        # Fallback to first available
        sample_key = list(dft_results.keys())[0]
        print(f'Warning: {sample_key} not found, using {sample_key}')
    
    z_dft, E_dft = dft_results[sample_key]
    print(f'Using DFT config: {sample_key}')
    print(f'  z range: [{z_dft.min():.2f}, {z_dft.max():.2f}] Å')
    print(f'  E range: [{E_dft.min():.2f}, {E_dft.max():.2f}] eV')
    
    # DFT energies are absolute total energies. Convert to interaction energies.
    # Reference: asymptotic value at large z (slab + molecule)
    E_asymptote = E_dft[-1]  # Last point (largest z)
    E_int = E_dft - E_asymptote
    print(f'  Asymptotic energy (z={z_dft[-1]:.2f} Å): {E_asymptote:.3f} eV')
    print(f'  Interaction energy range: [{E_int.min():.3f}, {E_int.max():.3f}] eV')

    # Build feature matrix for GridFF fitting
    # GridFF formula: E_total = sum_i (P_i * Pauli_i + L_i * London_i + Q_i * Coulomb_i)
    # where Pauli_i, London_i, Coulomb_i are sampled from grid channels at atom position
    
    # Sampled Es shape: (nconf, natoms, nch=4)
    # ch0: Pauli, ch1: London, ch2: Coulomb, ch3: Hbond
    
    # Approach: Fit only P (Pauli), with Coulomb fixed by known charges
    # London channel may not match DFT physics well, so skip it for now
    # E_int = sum_i (P_i * Pauli_i) + sum_i (Q_i * Coulomb_i)
    
    # Per-atom sampled values
    Pauli_per_atom = Es[:, :, 0]  # (nconf, natoms)
    London_per_atom = Es[:, :, 1]  # (nconf, natoms)
    Coulomb_per_atom = Es[:, :, 2]  # (nconf, natoms)
    
    # Diagnostics: print ranges of sampled values
    print(f'Diagnostics - sampled GridFF values:')
    print(f'  Pauli_per_atom range: [{Pauli_per_atom.min():.6f}, {Pauli_per_atom.max():.6f}]')
    print(f'  London_per_atom range: [{London_per_atom.min():.6f}, {London_per_atom.max():.6f}]')
    print(f'  Coulomb_per_atom range: [{Coulomb_per_atom.min():.6f}, {Coulomb_per_atom.max():.6f}]')
    
    # Compute Coulomb contribution (fixed by known charges)
    E_coulomb = np.sum(charges * Coulomb_per_atom, axis=1)  # (nconf,)
    
    # Subtract Coulomb to get Pauli reference
    E_P_ref = E_int - E_coulomb
    
    # Build feature matrix for P fitting only
    # Each row: [Pauli_O, Pauli_H] for that config
    A = np.zeros((len(T_z), ntypes), dtype=np.float64)
    for itype in range(ntypes):
        mask = (type_ids == itype)
        A[:, itype] = Pauli_per_atom[:, mask].sum(axis=1)  # Sum Pauli for this type
    
    # Fit P coefficients
    P_fit = _fit_linear(A, E_P_ref)
    
    print('--- Fit results (Pauli only) ---')
    for it, e in enumerate(uniq):
        print(f'  type {e}: P_fit={P_fit[it]:.6f}')
    
    # Predict and evaluate
    E_P_pred = A @ P_fit
    E_pred = E_P_pred + E_coulomb
    resid = (E_int - E_pred)
    rmse = float(np.sqrt(np.mean(resid**2)))
    mae = float(np.mean(np.abs(resid)))
    print(f'RMSE={rmse:.6e} eV  MAE={mae:.6e} eV')

    # Convert transforms to atomic positions for visualization
    apos_all = np.zeros((len(T_z), len(apos0), 3), dtype=np.float32)
    for i, T in enumerate(T_z):
        R = T[:, :3]
        t = T[:, 3]
        apos_all[i] = (apos0 @ R.T) + t

    # Diagnostic plots (reuse existing visualization)
    # Z-scan comparison
    fig, ax = plt.subplots(figsize=(7.0, 4.0))
    ax.plot(z_dft, E_int, 'o-', ms=4, lw=1.0, label='DFT (interaction)')
    ax.plot(z_vals, E_pred, 's--', ms=4, lw=1.0, label='Fit')
    ax.set_xlabel('z [Å]')
    ax.set_ylabel('E [eV]')
    ax.set_title(f'DFT vs Fit ({args.molecule_site}, RMSE={rmse:.3f} eV)')
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, 'zscan_comparison.png'), dpi=160)
    plt.close(fig)
    print('Saved z-scan comparison:', os.path.join(out_dir, 'zscan_comparison.png'))
    
    # Grid slice visualization with atom overlays
    if not args.no_diagnostics:
        print('Generating grid slice visualization...')
        
        # Plot density slice (channel 0) - substrate only
        su.plot_gridff_diagnostics(
            gridff[:,:,:,0:1], sub_apos, sub_enames, lvec_sub,
            iz_slices=[nz//4, nz//2, 3*nz//4],
            iy_slice=ny//2,
            save_path=os.path.join(out_dir, 'grid_density_slices.png'),
            g0=g0, dg=dg,
            channel_name='Pauli (Morse)'
        )
        print('Saved density slices:', os.path.join(out_dir, 'grid_density_slices.png'))
        
        # Plot potential slice (channel 2) - substrate only
        su.plot_gridff_diagnostics(
            gridff[:,:,:,2:3], sub_apos, sub_enames, lvec_sub,
            iz_slices=[nz//4, nz//2, 3*nz//4],
            iy_slice=ny//2,
            save_path=os.path.join(out_dir, 'grid_potential_slices.png'),
            g0=g0, dg=dg,
            channel_name='Coulomb (Ewald)'
        )
        print('Saved potential slices:', os.path.join(out_dir, 'grid_potential_slices.png'))
        
        # Alignment summary plot (matching test_gridff_alignment.py and test_fdbm_fit_gridff_mock.py)
        print('Generating alignment summary plot...')
        
        # Get PLQ coefficients for H and O from fitted P values
        # Use fitted P_fit as Pauli coefficient, L=0 for now (not fitted)
        # For alignment visualization, use unit coefficients to show raw grid
        plq_H = (1.0, 1.0, charges[np.where(np.array(names) == 'H')[0][0]])
        plq_O = (1.0, 1.0, charges[np.where(np.array(names) == 'O')[0][0]])
        
        # Use subset of molecule positions for visualization
        n_samples_plot = min(100, len(apos_all))
        mol_apos_subset = apos_all[:n_samples_plot]
        
        # H atom alignment plot
        su.plot_alignment_summary(
            gridff, g0, dg, sub_apos, sub_enames,
            save_path=os.path.join(out_dir, 'alignment_H.png'),
            iz_top=145, iy_center=ny//2,
            z_atom_range=5.0, mol_apos=mol_apos_subset, mol_enames=names,
            plq_coeffs=plq_H,
            zmin_offset=2.0, z_ylim=(-2.0, 2.0)
        )
        print('Saved H alignment:', os.path.join(out_dir, 'alignment_H.png'))
        
        # O atom alignment plot
        su.plot_alignment_summary(
            gridff, g0, dg, sub_apos, sub_enames,
            save_path=os.path.join(out_dir, 'alignment_O.png'),
            iz_top=145, iy_center=ny//2,
            z_atom_range=5.0, mol_apos=mol_apos_subset, mol_enames=names,
            plq_coeffs=plq_O,
            zmin_offset=2.0, z_ylim=(-2.0, 2.0)
        )
        print('Saved O alignment:', os.path.join(out_dir, 'alignment_O.png'))

    # Save data
    np.savez(os.path.join(out_dir, 'fdbm_dft_fit_results.npz'),
             transforms=T_z, names=names, charges=charges, type_ids=type_ids, uniq_types=np.array(uniq),
             Es=Es, E_dft=E_dft, E_pred=E_pred, E_coulomb=E_coulomb, resid=resid,
             P_fit=P_fit, z_vals=z_vals, dft_key=sample_key)

    print('Saved results:')
    print(' ', os.path.join(out_dir, 'fdbm_dft_fit_results.npz'))
    print(' ', os.path.join(out_dir, 'zscan_comparison.png'))


if __name__ == '__main__':
    main()
