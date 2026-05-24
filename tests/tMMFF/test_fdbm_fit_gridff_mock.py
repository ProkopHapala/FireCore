"""Test FDBM linear fitting using existing NaCl GridFF Bspline_PLQd as mock substrate.

This test intentionally reuses existing GPU sampling infrastructure:
- GridFF Bspline grid: cpp/common_resources/NaCl_1x1_L3/Bspline_PLQd.npy (+ meta JSON)
- OpenCL kernel: sampleGridFF_Bspline_points via pyBall.OCL.MolecularDynamics

We generate mock reference energies by composing the same sampled channels,
optionally adding controlled perturbations (noise, rescale, pauli_alpha != 1).

Outputs:
- tests/tMMFF/out_fdbm_mock/*png  (parity + z-scan)
- tests/tMMFF/out_fdbm_mock/*.npz (data + fitted coeffs)
"""

import os
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def _repo_root():
    here = os.path.dirname(os.path.abspath(__file__))
    return os.path.realpath(os.path.join(here, '..', '..'))


def _load_xyz_with_charges(xyz_path):
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


def _get_element_REQ(element_name, element_types_path):
    """Extract REQ parameters from ElementTypes.dat for a given element.
    
    CRITICAL: ElementTypes.dat stores EvdW (energy), but GridFF generation uses sqrt(EvdW) in REQ.y.
    When converting to PLQ for GridFF sampling, you MUST sqrt the E value.
    
    Returns: (R, EvdW, Q, H) where EvdW is the raw energy value from ElementTypes.dat.
    """
    with open(element_types_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 10 and parts[0] == element_name:
                R = float(parts[7])   # RvdW
                E = float(parts[8])   # EvdW (raw energy, NOT sqrt)
                Q = float(parts[9])   # Quff
                H = 0.0  # Not in ElementTypes.dat, set to 0
                return (R, E, Q, H)
    raise ValueError(f"Element {element_name} not found in {element_types_path}")


def _req_to_plq(req, alpha=1.5):
    """Convert REQ to PLQ coefficients using the same formula as RigidBodyDynamics._reqs_to_plq.
    
    CRITICAL: This function expects E to be sqrt(EvdW), NOT raw EvdW.
    If reading from ElementTypes.dat, you MUST sqrt the E value before calling this.
    
    Formula (matching C++ REQ2PLQ in Forces.h):
        e  = exp(alpha * R)
        cL = e * E              # London coefficient
        cP = e * cL = e^2 * E   # Pauli coefficient
        cH = e^2 * H           # H-bond coefficient (usually 0)
    
    The sqrt(E) convention ensures proper mixed interaction:
        Eij = sqrt(Eii * Ejj) when GridFF channels contain substrate sqrt(Ejj)
    
    Args:
        req: (R, sqrt(EvdW), Q, H) - E MUST be sqrt(EvdW) !!!
        alpha: alphaMorse parameter (default 1.8, must match GridFF generation)
    
    Returns:
        (cP, cL, Q, cH) - PLQ coefficients for GridFF sampling
    """
    R, E, Q, H = req
    e = np.exp(alpha * R)
    cL = e * E
    cP = e * cL
    cH = e * e * H
    return (cP, cL, Q, cH)


def _rot_z(phi):
    c = float(np.cos(phi)); s = float(np.sin(phi))
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]], dtype=np.float32)


def _make_transforms_zscan(xy=(1.0, 1.0), z_vals=None, phis=None):
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


def _make_transforms_random(n, xy_box=((0.0, 4.0), (0.0, 4.0)), z_range=(0.5, 6.0), seed=1):
    from pyBall.atomicUtils import rotation_matrix
    rng = np.random.default_rng(int(seed))
    xs = rng.uniform(float(xy_box[0][0]), float(xy_box[0][1]), size=n)
    ys = rng.uniform(float(xy_box[1][0]), float(xy_box[1][1]), size=n)
    zs = rng.uniform(float(z_range[0]), float(z_range[1]), size=n)
    # Random 3D rotations: random axis + random angle
    axes = rng.normal(size=(n, 3))
    axes = axes / np.linalg.norm(axes, axis=1, keepdims=True)
    angles = rng.uniform(0.0, 2.0*np.pi, size=n)
    Ts = np.zeros((n, 3, 4), dtype=np.float32)
    for i in range(n):
        Ts[i, :, :3] = rotation_matrix(axes[i], angles[i]).astype(np.float32)
        Ts[i, :, 3] = (xs[i], ys[i], zs[i])
    return Ts


def _fit_linear(A, b):
    x, *_ = np.linalg.lstsq(A, b, rcond=None)
    return x


def main():
    import argparse
    parser = argparse.ArgumentParser(description='FDBM mock-fit test with GridFF')
    parser.add_argument('--no-diagnostics', action='store_true', help='Skip diagnostic plots')
    parser.add_argument('--zmin-offset', type=float, default=2.0, help='Z offset for XZ color scale (default 2.0 A)')
    parser.add_argument('--z-ylim', type=float, nargs=2, default=(-0.5, 0.5), help='Z-profile y-limits (default -0.5 0.5 eV)')
    parser.add_argument('--xy-range', type=float, nargs=2, default=(-2.0, 2.0), help='XY range for random configs (default -2.0 2.0 A)')
    parser.add_argument('--z-range', type=float, nargs=2, default=(-1.0, 5.0), help='Z range for random configs (default -1.0 5.0 A)')
    args = parser.parse_args()

    ROOT = _repo_root()
    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out_fdbm_mock')
    os.makedirs(out_dir, exist_ok=True)

    from pyBall.OCL import Surface_utils as su

    grid_path = os.path.join(ROOT, 'cpp', 'common_resources', 'NaCl_1x1_L3', 'Bspline_PLQd.npy')
    mol_path = os.path.join(ROOT, 'cpp', 'common_resources', 'xyz', 'H2O.xyz')

    print('=== FDBM mock-fit test (GridFF Bspline_PLQd) ===')
    print('grid_path=', grid_path)
    print('mol_path =', mol_path)
    print('out_dir  =', out_dir)

    apos, names, charges = _load_xyz_with_charges(mol_path)
    # Rigid body frame: center molecule (prevents translation mixing into sampling)
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

    # True parameters (we will try to recover)
    use_london_fit = True
    P_true = np.array([1.30, 0.60], dtype=np.float64)  # (O,H)
    L_true = np.array([0.80, 0.20], dtype=np.float64)  # (O,H)

    # Configs: random + z-scan diagnostic
    xy_box = ((args.xy_range[0], args.xy_range[1]), (args.xy_range[0], args.xy_range[1]))
    T_rand = _make_transforms_random(800, xy_box=xy_box, z_range=args.z_range, seed=3)
    z_vals = np.linspace(args.z_range[0], args.z_range[1], 66)
    T_z, zphi = _make_transforms_zscan(xy=(0.0, 0.0), z_vals=z_vals, phis=[0.0])
    transforms = np.concatenate([T_rand, T_z], axis=0)
    print('nconf total=', len(transforms), ' (rand=', len(T_rand), ' zscan=', len(T_z), ')')

    # Initialize sampler ONCE (uploads grid once)
    md, meta = su.init_gridff_sampler_md(grid_path, apos0=apos0, nSystems=1024, use_texture=False)

    # Load substrate for alignment visualization
    sub_path = os.path.join(ROOT, 'cpp', 'common_resources', 'xyz', 'NaCl_1x1_L3.xyz')
    sub_data = su.load_substrate_xyz_with_lvec(sub_path)
    sub_apos = sub_data['apos']
    sub_enames = sub_data['enames']
    lvec = sub_data['lvec']

    # Sample channels (P,L,Q) for all configs
    Es = su.sample_gridff_channels_rigid(md, transforms)
    if Es.shape[2] < 3:
        raise RuntimeError(f'Expected >=3 channels (P,L,Q), got Es.shape={Es.shape}')

    # Mock reference with perturbations
    mock_opts = {
        'use_london': True,
        'pauli_alpha': 1.00,
        'pauli_rescale': 1.00,
        'london_rescale': 1.00,
        'coulomb_rescale': 1.00,
        'noise_sigma': 0.00,
    }
    E_ref, parts = su.fdbm_make_mock_reference(Es, charges=charges, type_ids=type_ids, P_true=P_true, L_true=L_true, **mock_opts)

    # Build linear system for fit of P,L (Coulomb fixed)
    # E_ref = A*[P_types, L_types] + E_coul
    E_coul = parts['coulomb']
    b = (E_ref - E_coul).astype(np.float64)
    A = su.fdbm_build_feature_matrix(Es, type_ids=type_ids, ntypes=ntypes, use_london=use_london_fit)

    x = _fit_linear(A, b)
    if use_london_fit:
        P_fit = x[:ntypes]
        L_fit = x[ntypes:]
    else:
        P_fit = x
        L_fit = np.zeros(ntypes)

    print('--- Fit results ---')
    for it, e in enumerate(uniq):
        print(f'  type {e}: P_fit={P_fit[it]: .6f}  P_true={P_true[it]: .6f}   L_fit={L_fit[it]: .6f}  L_true={L_true[it]: .6f}')

    # Predict energies and stats
    E_pred = (A @ x) + E_coul
    resid = (E_ref - E_pred)
    rmse = float(np.sqrt(np.mean(resid**2)))
    mae = float(np.mean(np.abs(resid)))
    print(f'RMSE={rmse:.6e} eV  MAE={mae:.6e} eV')

    # Convert transforms to actual atomic positions for XYZ save
    # T[:, :3] = rotation, T[:, 3] = translation
    # apos = apos0 @ R.T + t
    apos_all = np.zeros((len(transforms), len(apos0), 3), dtype=np.float32)
    for i, T in enumerate(transforms):
        R = T[:, :3]
        t = T[:, 3]
        apos_all[i] = (apos0 @ R.T) + t

    # Load GridFF for diagnostic visualization
    grid, meta_grid = su.load_bspline_gridff(grid_path)
    g0 = np.array(meta_grid['g0'], dtype=np.float32)
    dg = np.array(meta_grid['dg'], dtype=np.float32)

    # Get proper REQ parameters from ElementTypes.dat and convert to PLQ
    # ElementTypes.dat stores EvdW, but GridFF generation uses sqrt(EvdW) in REQ.y
    element_types_path = os.path.join(ROOT, 'cpp', 'common_resources', 'ElementTypes.dat')
    req_H = _get_element_REQ('H', element_types_path)
    req_O = _get_element_REQ('O', element_types_path)

    # Replace Q with actual partial charges from H2O.xyz
    # charges array corresponds to names order: ['O', 'H', 'H']
    charge_O = charges[np.where(np.array(names) == 'O')[0][0]]
    charge_H = charges[np.where(np.array(names) == 'H')[0][0]]
    
    # IMPORTANT: ElementTypes.dat stores EvdW, but GridFF uses sqrt(EvdW) in REQ.y
    # We must sqrt the E value to match the GridFF generation convention
    req_H = (req_H[0], np.sqrt(req_H[1]), charge_H, req_H[3])
    req_O = (req_O[0], np.sqrt(req_O[1]), charge_O, req_O[3])

    # Use alpha_morse=1.8 to match the regenerated GridFF
    plq_H = _req_to_plq(req_H, alpha=1.8)  # (cP, cL, Q, cH)
    plq_O = _req_to_plq(req_O, alpha=1.8)  # (cP, cL, Q, cH)
    print(f'H REQ: {req_H} -> PLQ (alpha=1.8): P={plq_H[0]:.6f}, L={plq_H[1]:.6f}, Q={plq_H[2]:.6f}')
    print(f'O REQ: {req_O} -> PLQ (alpha=1.8): P={plq_O[0]:.6f}, L={plq_O[1]:.6f}, Q={plq_O[2]:.6f}')

    # Plot total potential summary with molecule samples overlay (use subset to avoid clutter)
    n_samples_plot = min(100, len(apos_all))
    mol_apos_subset = apos_all[:n_samples_plot]

    # H atom plot
    diag_path_h = os.path.join(out_dir, 'total_potential_H_atom.png')
    su.plot_alignment_summary(grid, g0, dg, sub_apos, sub_enames,
                              save_path=diag_path_h, iz_top=33, iy_center=20,
                              z_atom_range=3.0, mol_apos=mol_apos_subset, mol_enames=names,
                              plq_coeffs=(plq_H[0], plq_H[1], plq_H[2]),
                              zmin_offset=args.zmin_offset, z_ylim=args.z_ylim,
                              plot_diagnostics=not args.no_diagnostics)
    if not args.no_diagnostics:
        print('Saved H atom diagnostic plot:', diag_path_h)

    # O atom plot
    diag_path_o = os.path.join(out_dir, 'total_potential_O_atom.png')
    su.plot_alignment_summary(grid, g0, dg, sub_apos, sub_enames,
                              save_path=diag_path_o, iz_top=33, iy_center=20,
                              z_atom_range=3.0, mol_apos=mol_apos_subset, mol_enames=names,
                              plq_coeffs=(plq_O[0], plq_O[1], plq_O[2]),
                              zmin_offset=args.zmin_offset, z_ylim=args.z_ylim,
                              plot_diagnostics=not args.no_diagnostics)
    if not args.no_diagnostics:
        print('Saved O atom diagnostic plot:', diag_path_o)

    # Save XYZ movie with energies in comment lines
    xyz_path = os.path.join(out_dir, 'mock_configs.xyz')
    su.save_xyz_movie_with_energies(xyz_path, names, apos_all, E_ref, qs=charges)
    print('Saved XYZ movie:', xyz_path)

    # Verify by loading back
    enames_back, apos_back, energies_back, qs_back = su.load_xyz_movie_with_energies(xyz_path)
    print('Loaded back XYZ movie: nframes=', len(energies_back), ' natoms=', len(enames_back))
    assert np.allclose(energies_back, E_ref, rtol=1e-8), 'Energies mismatch after round-trip'
    max_pos_diff = np.max(np.abs(apos_back - apos_all))
    print(f'Max position difference after round-trip: {max_pos_diff:.2e}')
    assert np.allclose(apos_back, apos_all, rtol=1e-5, atol=1e-6), 'Positions mismatch after round-trip'
    print('Round-trip save/load verified')

    # Save data (NPZ still useful for full feature matrix etc.)
    np.savez(os.path.join(out_dir, 'fdbm_mock_fit_results.npz'),
             transforms=transforms, names=names, charges=charges, type_ids=type_ids, uniq_types=np.array(uniq),
             Es=Es, E_ref=E_ref, E_pred=E_pred, E_coul=E_coul, resid=resid,
             P_true=P_true, L_true=L_true, P_fit=P_fit, L_fit=L_fit,
             mock_opts=mock_opts, grid_meta=meta)

    # Plot parity
    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    ax.plot(E_ref, E_pred, 'k.', ms=3, alpha=0.6)
    v0 = float(min(E_ref.min(), E_pred.min()))
    v1 = float(max(E_ref.max(), E_pred.max()))
    ax.plot([v0, v1], [v0, v1], 'r--', lw=1)
    ax.set_xlabel('E_ref [eV]')
    ax.set_ylabel('E_pred [eV]')
    ax.set_title(f'FDBM mock fit parity (RMSE={rmse:.2e} eV)')
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, 'parity.png'), dpi=160)
    plt.close(fig)

    # Plot z-scan
    iz0 = len(T_rand)
    zscan_ref = E_ref[iz0:]
    zscan_pred = E_pred[iz0:]
    fig, ax = plt.subplots(figsize=(7.0, 4.0))
    ax.plot(z_vals, zscan_ref, 'o-', ms=3, lw=1.0, label='ref')
    ax.plot(z_vals, zscan_pred, 's--', ms=3, lw=1.0, label='fit')
    ax.set_xlabel('z [Å]')
    ax.set_ylabel('E [eV]')
    ax.set_title('Z-scan (x=y=1Å, phi=0)')
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, 'zscan.png'), dpi=160)
    plt.close(fig)

    print('Saved results:')
    print(' ', os.path.join(out_dir, 'fdbm_mock_fit_results.npz'))
    print(' ', os.path.join(out_dir, 'parity.png'))
    print(' ', os.path.join(out_dir, 'zscan.png'))

    # Reasonableness check
    # In the unperturbed case (alpha=1, noise=0) this should fit almost exactly.
    if (mock_opts['pauli_alpha'] == 1.0) and (mock_opts['noise_sigma'] == 0.0):
        if rmse > 1e-5:
            raise RuntimeError(f'Unexpectedly large RMSE={rmse} for self-fit with no perturbations')


if __name__ == '__main__':
    main()
