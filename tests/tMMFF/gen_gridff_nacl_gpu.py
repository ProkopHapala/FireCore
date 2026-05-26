"""Generate Bspline_PLQd.npy for NaCl using GPU Morse + DFT Coulomb.

Uses existing GridFF_cl GPU infrastructure:
- Pauli/London: computed via make_MorseFF kernel (GPU)
- Coulomb: loaded from DFT electrostatic_potential.cube
- Hbond: 0
"""

import os
import sys
import json
import numpy as np

import matplotlib
matplotlib.use('Agg')


def _repo_root():
    here = os.path.dirname(os.path.abspath(__file__))
    return os.path.realpath(os.path.join(here, '..', '..'))


sys.path.insert(0, _repo_root())

bohr2ang = 0.5291772109217
hartree2ev = 27.211396132


def _load_cube(cube_path, is_potential):
    """Load cube file and convert units.

    - potential: Hartree -> eV
    - density:   e/Bohr^3 -> e/Angstrom^3

    Returns array in (nx, ny, nz) C-order (z fastest).
    """
    with open(cube_path, 'r') as f:
        f.readline(); f.readline()
        line3 = f.readline().split()
        n_atoms = int(line3[0])
        origin = np.array([float(line3[1]), float(line3[2]), float(line3[3])]) * bohr2ang
        
        line4 = f.readline().split()
        nx = int(line4[0])
        dx = float(line4[1]) * bohr2ang
        
        line5 = f.readline().split()
        ny = int(line5[0])
        dy = float(line5[2]) * bohr2ang
        
        line6 = f.readline().split()
        nz = int(line6[0])
        dz = float(line6[3]) * bohr2ang
        
        for _ in range(n_atoms):
            f.readline()
        
        n_vals = nx * ny * nz
        data = np.fromfile(f, sep=' ', count=n_vals)
        data = data.reshape((nx, ny, nz))  # C-order: z fastest
        if is_potential:
            # DFT electrostatic potential is already in Volts (or eV), NOT Hartree
            # No conversion needed
            pass
        else:
            # Density: convert e/Bohr^3 to e/Angstrom^3
            data = data / (bohr2ang ** 3)
        return data, origin, (nx, ny, nz), (dx, dy, dz)


def main():
    ROOT = _repo_root()
    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out_fdbm_dft_gridff')
    os.makedirs(out_dir, exist_ok=True)
    
    dft_data_dir = os.path.join(ROOT, 'tests', 'tSurf', 'small_mols_NaCl_New')
    density_cube   = os.path.join(dft_data_dir, '4-elec_pot_chr_den', 'charge_density.cube')
    potential_cube = os.path.join(dft_data_dir, '4-elec_pot_chr_den', 'electrostatic_potential.cube')
    sub_path = os.path.join(dft_data_dir, '0-geoms', 'NaCl.xyz')
    
    print('=== Generating GridFF for NaCl (DFT cube density+potential) ===')
    
    # Load DFT density & potential (authoritative)
    rho, g0_rho, (nx, ny, nz), dg_rho = _load_cube(density_cube, is_potential=False)
    phi, g0_phi, n_phi, dg_phi = _load_cube(potential_cube, is_potential=True)
    if tuple(n_phi) != (nx, ny, nz):
        raise RuntimeError(f'Cube shape mismatch: density {nx,ny,nz} vs potential {tuple(n_phi)}')
    if not (np.allclose(g0_rho, g0_phi) and np.allclose(dg_rho, dg_phi)):
        raise RuntimeError(f'Cube grid mismatch: g0_rho={g0_rho} g0_phi={g0_phi} dg_rho={dg_rho} dg_phi={dg_phi}')
    g0 = g0_rho
    dg = dg_rho
    # Correct electrostatic potential: subtract V0 (average of top slice at z = end of cell)
    # This removes the arbitrary constant from FFT Poisson solver
    V0 = phi[:, :, -1].mean()
    phi = phi - V0
    print(f'Corrected electrostatics by subtracting V0={V0:.6f} eV (top slice average)')
    # Invert sign: DFT Hartree potential is for electrons (opposite sign to physical electrostatic)
    phi = -phi
    print(f'Inverted sign of electrostatic potential (Hartree for electrons -> physical)')
    
    print(f'DFT grid: ({nx}, {ny}, {nz}), g0={g0}, dg={dg}')
    print(f'DFT density range:   [{rho.min():.6e}, {rho.max():.6e}] e/A^3')
    print(f'DFT potential range: [{phi.min():.6f}, {phi.max():.6f}] eV (after V0 correction)')
    
    # Load substrate atoms
    from pyBall.OCL.Surface_utils import load_substrate_xyz_with_lvec
    sub_data = load_substrate_xyz_with_lvec(sub_path)
    sub_apos = sub_data['apos']
    sub_enames = sub_data['enames']
    # lvec from substrate might be None, use DFT grid lattice instead
    lvec = sub_data['lvec']
    if lvec is None:
        # Build lattice from DFT grid parameters
        lvec = np.zeros((4, 3), dtype=np.float64)
        lvec[0] = g0
        lvec[1] = [nx * dg[0], 0.0, 0.0]
        lvec[2] = [0.0, ny * dg[1], 0.0]
        lvec[3] = [0.0, 0.0, nz * dg[2]]

    print(f'Substrate: {len(sub_apos)} atoms')

    # Channels from DFT cubes (authoritative)
    pauli   = np.ascontiguousarray(rho.astype(np.float32))
    london  = np.zeros((nx, ny, nz), dtype=np.float32)
    coulomb = np.ascontiguousarray(phi.astype(np.float32))

    # Shift g0[2] so that z=0 is at the top substrate atom
    z_top_sub = sub_apos[:, 2].max()
    g0 = np.array([g0[0], g0[1], g0[2] - z_top_sub], dtype=np.float32)
    print(f'Shifted g0 so z=0 at top substrate atom (z_top={z_top_sub:.3f} A)')
    
    # Assemble GridFF: (nx, ny, nz, 4) as float32 contiguous
    gridff = np.ascontiguousarray(np.stack([pauli, london, coulomb, np.zeros_like(pauli)], axis=3).astype(np.float32))
    
    # Save
    gridff_path = os.path.join(dft_data_dir, 'Bspline_PLQd.npy')
    np.save(gridff_path, gridff)
    print(f'\nSaved GridFF to: {gridff_path}')
    
    # Save metadata
    meta = {
        'ns': [nx, ny, nz],
        'g0': g0.tolist(),
        'dg': list(dg),
        'grid_type': 'Bspline_PLQd',
        'generation_script': 'gen_gridff_nacl_gpu.py',
        'channels': ['Pauli (charge_density.cube)', 'London (0)', 'Coulomb (electrostatic_potential.cube)', 'Hbond'],
        'note': 'Pauli channel loaded from DFT charge_density.cube (converted e/Bohr^3 -> e/A^3). Coulomb channel loaded from DFT electrostatic_potential.cube (Hartree -> eV).'
    }
    meta_path = os.path.join(dft_data_dir, 'Bspline_PLQd_meta.json')
    with open(meta_path, 'w') as f:
        json.dump(meta, f, indent=2)
    print(f'Saved metadata to: {meta_path}')
    
    # Generate diagnostic plots
    from pyBall.OCL import Surface_utils as su
    
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
    
    mol_path = os.path.join(ROOT, 'cpp', 'common_resources', 'xyz', 'H2O.xyz')
    apos, names, charges = _load_xyz_with_charges(mol_path)
    apos0 = (apos - apos.mean(axis=0)).astype(np.float32)
    
    z_vals = np.linspace(1.4, 10.0, 33)
    def _make_transforms_zscan(xy=(0.0, 0.0), z_vals=None, phis=None):
        if z_vals is None:
            z_vals = np.linspace(1.4, 10.0, 33)
        if phis is None:
            phis = np.array([0.0], dtype=np.float32)
        z_vals = np.asarray(z_vals, dtype=np.float32)
        phis = np.asarray(phis, dtype=np.float32)
        Ts = []
        for phi in phis:
            c = float(np.cos(phi)); s = float(np.sin(phi))
            R = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]], dtype=np.float32)
            for z in z_vals:
                T = np.zeros((3, 4), dtype=np.float32)
                T[:, :3] = R
                T[:, 3] = np.array([xy[0], xy[1], float(z)], dtype=np.float32)
                Ts.append(T)
        return np.array(Ts, dtype=np.float32)
    
    T_z = _make_transforms_zscan(xy=(0.0, 0.0), z_vals=z_vals, phis=[0.0])
    apos_all = np.zeros((len(T_z), len(apos0), 3), dtype=np.float32)
    for i, T in enumerate(T_z):
        R = T[:, :3]
        t = T[:, 3]
        apos_all[i] = (apos0 @ R.T) + t
    
    n_samples_plot = min(100, len(apos_all))
    mol_apos_subset = apos_all[:n_samples_plot]

    from pyBall.OCL.MMparams import read_element_types
    element_types = read_element_types(os.path.join(ROOT, 'cpp', 'common_resources', 'ElementTypes.dat'))
    alpha_morse = 1.8
    
    def _req_to_plq(req, alpha=1.8):
        R, E, Q, H = req
        e = np.exp(alpha * R)
        cL = e * E
        cP = e * cL
        return (cP, cL, Q, 0.0)
    
    req_H = (element_types['H'].RvdW, np.sqrt(element_types['H'].EvdW), 
             charges[np.where(np.array(names) == 'H')[0][0]], 0.0)
    req_O = (element_types['O'].RvdW, np.sqrt(element_types['O'].EvdW), 
             charges[np.where(np.array(names) == 'O')[0][0]], 0.0)
    
    plq_H = _req_to_plq(req_H, alpha=alpha_morse)
    plq_O = _req_to_plq(req_O, alpha=alpha_morse)
    
    # iz_top at z=2 A above surface: iz = (2.0 - g0[2]) / dg[2]
    iz_top = int((2.0 - g0[2]) / dg[2])
    print(f'Using iz_top={iz_top} (z={g0[2] + iz_top*dg[2]:.3f} A)')

    # Single combined plot with H as primary (XY/XZ) and both H+O in Z-profiles
    diag_path_h = os.path.join(out_dir, 'total_potential_H_atom.png')
    su.plot_alignment_summary(gridff, g0, dg, sub_apos, sub_enames,
                              save_path=diag_path_h, iz_top=iz_top, iy_center=64,
                              z_atom_range=5.0, mol_apos=mol_apos_subset, mol_enames=names,
                              plq_coeffs=(plq_H[0], plq_H[1], plq_H[2]),
                              plq_coeffs2=(plq_O[0], plq_O[1], plq_O[2]),
                              zmin_offset=2.0, z_ylim=(-2.0, 2.0), vmin_vmax_xz=(-0.5, 0.5))
    print('Saved H atom plot:', diag_path_h)

    # O as primary (XY/XZ) and both H+O in Z-profiles
    diag_path_o = os.path.join(out_dir, 'total_potential_O_atom.png')
    su.plot_alignment_summary(gridff, g0, dg, sub_apos, sub_enames,
                              save_path=diag_path_o, iz_top=iz_top, iy_center=64,
                              z_atom_range=5.0, mol_apos=mol_apos_subset, mol_enames=names,
                              plq_coeffs=(plq_O[0], plq_O[1], plq_O[2]),
                              plq_coeffs2=(plq_H[0], plq_H[1], plq_H[2]),
                              zmin_offset=2.0, z_ylim=(-2.0, 2.0), vmin_vmax_xz=(-0.5, 0.5))
    print('Saved O atom plot:', diag_path_o)

    # === Fit P, L and charges from DFT data ===
    print('\n=== Fitting P, L coefficients from DFT data ===')

    # Load consolidated DFT data (use both H2O-H and H2O-O if available)
    dft_paths = [
        os.path.join(ROOT, 'tests', 'tMMFF', 'H2O-H_matched.npz'),
        os.path.join(ROOT, 'tests', 'tMMFF', 'H2O-O_matched.npz'),
    ]
    dsets = []
    for p in dft_paths:
        if os.path.exists(p):
            dsets.append((os.path.basename(p).replace('_matched.npz', ''), np.load(p)))
        else:
            print(f'Warning: DFT data not found at {p}')
    if len(dsets) == 0:
        print('No DFT matched datasets found. Run: python consolidate_dft_data.py match --molecule H2O-H (and H2O-O)')
    else:
        # concatenate datasets
        coords = np.concatenate([ds['coords'] for _, ds in dsets], axis=0).copy()
        energies = np.concatenate([ds['energies'] for _, ds in dsets], axis=0)
        enames = dsets[0][1]['enames']
        ix = np.concatenate([ds['ix'] for _, ds in dsets], axis=0)
        iy = np.concatenate([ds['iy'] for _, ds in dsets], axis=0)
        orientxy = np.concatenate([ds['orientxy'] for _, ds in dsets], axis=0)
        orientz = np.concatenate([ds['orientz'] for _, ds in dsets], axis=0)
        zdist = np.concatenate([ds['zdist'] for _, ds in dsets], axis=0)
        mol_tag = np.concatenate([np.full(len(ds['energies']), name, dtype=object) for name, ds in dsets], axis=0)

        # CRITICAL: grid origin g0 was shifted so that z=0 is at top substrate atom.
        # DFT xyz coordinates are still in the original frame (same as substrate xyz).
        # Therefore shift molecule coordinates by the same z_top_sub before sampling.
        coords[:, :, 2] -= z_top_sub

        print(f'Loaded {len(coords)} configurations from DFT data')
        print(f'Element names: {list(enames)}')
        print(f'Energy range: [{energies.min():.6f}, {energies.max():.6f}] eV')

        # Compute interaction energy by subtracting baseline at farthest z within each group.
        # Group key: (mol, ix, iy, orientxy, orientz)
        E0 = np.zeros_like(energies, dtype=np.float64)
        uniq_keys = {}
        for i in range(len(energies)):
            k = (str(mol_tag[i]), int(ix[i]), int(iy[i]), str(orientxy[i]), int(orientz[i]))
            uniq_keys.setdefault(k, []).append(i)
        for k, idxs in uniq_keys.items():
            zz = zdist[idxs]
            zmax = float(np.max(zz))
            mask = np.isclose(zz, zmax, rtol=0.0, atol=1e-6)
            imax = np.array(idxs, dtype=int)[mask]
            if len(imax) == 0: raise RuntimeError(f'No max-z samples for group {k} zmax={zmax}')
            E0_val = float(np.mean(energies[imax]))
            E0[idxs] = E0_val
        E_int = energies - E0
        print(f'Interaction energy range (E-E0): [{E_int.min():.6f}, {E_int.max():.6f}] eV')

        # Restrict fit: only tilt=0, specific scan configs, E_int <= Eint_max
        # Also apply exponential weighting to prioritize near-minimum energies
        Eint_max = 0.5
        scan_tilt_fit = 0
        # scan configs: (mol, ix, iy) that we want to fit
        # Focus only on the two key bonding profiles:
        # - H2O-H over Cl (H-down hydrogen bond)
        # - H2O-O over Na (electrostatic minimum)
        fit_scan_configs = [
            ('H2O-H', 0, 0),  # Cl site
            ('H2O-O', 6, 0),  # Na site
        ]
        mask_scan = np.zeros(len(E_int), dtype=bool)
        for mol_f, ix_f, iy_f in fit_scan_configs:
            mask_scan |= (mol_tag == mol_f) & (ix == ix_f) & (iy == iy_f) & (orientz == scan_tilt_fit)
        mask_fit = mask_scan & (E_int <= Eint_max)
        n_keep = int(np.count_nonzero(mask_fit))
        if n_keep == 0: raise RuntimeError(f'No samples satisfy fit criteria')
        print(f'Fit filter: tilt={scan_tilt_fit}, scan configs only, E_int <= {Eint_max:.3f} eV  ->  {n_keep} / {len(E_int)} configs')
        # Exponential weighting: w = exp(-E_int / kT), prioritize near-minimum
        kT_weight = 0.2  # eV
        weights = np.exp(-E_int / kT_weight)
        weights[~mask_fit] = 0.0
        weights[mask_fit] /= weights[mask_fit].max()  # normalize so max weight = 1
        print(f'Weighting: exp(-E_int/{kT_weight:.2f}), effective weight range [{weights[mask_fit].min():.4f}, {weights[mask_fit].max():.4f}]')
        # crude per-dataset stats
        for name, _ in dsets:
            m = (mol_tag == name) & mask_fit
            nk = int(np.count_nonzero(m))
            nt = int(np.count_nonzero(mol_tag == name))
            print(f'  {name}: {nk} / {nt} kept')

        # Sample GridFF channels at atom positions for all configurations
        n_configs = len(coords)
        n_atoms = coords.shape[1]
        n_types = 2

        idx_pauli = 0
        idx_coulomb = 2

        print(f'Sampling GridFF channels at {n_configs * n_atoms} atom positions...')

        n_clamp = 0
        n_wrapx = 0
        n_wrapy = 0
        ix_min =  10**9; iy_min =  10**9; iz_min =  10**9
        ix_max = -10**9; iy_max = -10**9; iz_max = -10**9

        def sample_gridff(pos, gridff, g0, dg):
            """Sample GridFF numpy array at position (x, y, z) using trilinear interpolation.

            - x/y: periodic wrap
            - z: clamp

            Returns:
                vals: (nch,) float32
                clamped_z: bool
                (ix0,iy0,iz0): int base grid indices before wrap/clamp (floor)
                (ix,iy,iz): int wrapped/clamped base indices used on grid
            """
            nx, ny, nz, nch = gridff.shape
            fx = (pos[0] - g0[0]) / dg[0]
            fy = (pos[1] - g0[1]) / dg[1]
            fz = (pos[2] - g0[2]) / dg[2]
            ix0 = int(np.floor(fx)); iy0 = int(np.floor(fy)); iz0 = int(np.floor(fz))
            tx = fx - ix0; ty = fy - iy0; tz = fz - iz0

            ix = ix0 % nx
            iy = iy0 % ny
            clamped_z = (iz0 < 0) or (iz0 >= nz - 1)
            iz = max(0, min(iz0, nz - 2))
            iz1 = iz + 1

            ix1 = (ix + 1) % nx
            iy1 = (iy + 1) % ny

            # gather 8 corners
            c000 = gridff[ix , iy , iz , :]
            c100 = gridff[ix1, iy , iz , :]
            c010 = gridff[ix , iy1, iz , :]
            c110 = gridff[ix1, iy1, iz , :]
            c001 = gridff[ix , iy , iz1, :]
            c101 = gridff[ix1, iy , iz1, :]
            c011 = gridff[ix , iy1, iz1, :]
            c111 = gridff[ix1, iy1, iz1, :]

            # trilinear interpolation
            c00 = c000*(1.0-tx) + c100*tx
            c10 = c010*(1.0-tx) + c110*tx
            c01 = c001*(1.0-tx) + c101*tx
            c11 = c011*(1.0-tx) + c111*tx
            c0  = c00 *(1.0-ty) + c10 *ty
            c1  = c01 *(1.0-ty) + c11 *ty
            vals = c0  *(1.0-tz) + c1  *tz
            return vals, clamped_z, (ix0, iy0, iz0), (ix, iy, iz)

        pauli_samp = []
        coul_samp  = []
        b = E_int
        idx_fit = np.where(mask_fit)[0]
        if len(idx_fit) == 0: raise RuntimeError('mask_fit is empty')
        w_fit = weights[idx_fit]
        b_fit0 = b[idx_fit]
        sqw = np.sqrt(w_fit)

        dz_grid = np.linspace(-1.0, 1.0, 81)
        best = None
        for idz, dz_shift_fit in enumerate(dz_grid):
            # accumulate only for fit configurations (massive speedup)
            nfit = len(idx_fit)
            rhoH = np.zeros(nfit, dtype=np.float64)
            rhoO = np.zeros(nfit, dtype=np.float64)
            phiH = np.zeros(nfit, dtype=np.float64)
            phiO = np.zeros(nfit, dtype=np.float64)
            for ii, i in enumerate(idx_fit):
                for j, ename in enumerate(enames):
                    pos = coords[i, j].copy(); pos[2] += dz_shift_fit
                    vals, clamped_z, (ix0, iy0, iz0), (ixg, iyg, izg) = sample_gridff(pos, gridff, g0, dg)
                    if idz == 0:
                        if clamped_z: n_clamp += 1
                        if ix0 != ixg: n_wrapx += 1
                        if iy0 != iyg: n_wrapy += 1
                        if ixg < ix_min: ix_min = ixg
                        if iyg < iy_min: iy_min = iyg
                        if izg < iz_min: iz_min = izg
                        if ixg > ix_max: ix_max = ixg
                        if iyg > iy_max: iy_max = iyg
                        if izg > iz_max: iz_max = izg
                        pauli_samp.append(float(vals[idx_pauli])); coul_samp.append(float(vals[idx_coulomb]))
                    if ename == 'H':
                        rhoH[ii] += float(vals[idx_pauli])
                        phiH[ii] += float(vals[idx_coulomb])
                    elif ename == 'O':
                        rhoO[ii] += float(vals[idx_pauli])
                        phiO[ii] += float(vals[idx_coulomb])

            aH = rhoH * sqw
            aO = rhoO * sqw
            aC = (phiH - 2.0 * phiO) * sqw
            b_fit = b_fit0 * sqw

            qH_fit = float(plq_H[2]); qO_fit = float(-2.0 * qH_fit)
            bW = b_fit - qH_fit * aC

            denomO = float(np.dot(aO, aO))
            if denomO <= 0.0: continue

            P_H_grid = np.linspace(0.0, 60.0, 601)
            for PH in P_H_grid:
                rhs = bW - PH * aH
                PO = float(np.dot(aO, rhs) / denomO)
                if PO < 0.0: PO = 0.0
                d = aH*PH + aO*PO - bW
                r = float(np.dot(d, d))
                if (best is None) or (r < best[0]):
                    best = (r, float(dz_shift_fit), float(PH), float(PO))

        if best is None: raise RuntimeError('Constrained fit failed')
        residuals, dz_shift_fit, PH_fit, PO_fit = best
        qH_fit = float(plq_H[2]); qO_fit = float(-2.0 * qH_fit)
        c = np.array([PH_fit, PO_fit, qH_fit], dtype=np.float64)

        print(f'Constrained Pauli fit (q fixed, P_H>=0,P_O>=0): dz in [{dz_grid.min():.2f},{dz_grid.max():.2f}] found dz={dz_shift_fit:.3f}')
        rank = 2

        print(f'\n=== Fitted coefficients ===')
        print(f'P_H = {c[0]:.6f} (diagnostic from atom params: {plq_H[0]:.6f})')
        print(f'P_O = {c[1]:.6f} (diagnostic from atom params: {plq_O[0]:.6f})')
        print(f'q_H = {qH_fit:.6f} (RESP ref: {plq_H[2]:.6f})')
        print(f'q_O = {qO_fit:.6f} (RESP ref: {plq_O[2]:.6f})')
        print(f'residuals: {residuals:.6f}')

        # Save fitted coefficients
        fit_coeffs = {
            'P_H': float(c[0]), 'P_O': float(c[1]),
            'q_H': qH_fit, 'q_O': qO_fit,
            'diagnostic_P_H': float(plq_H[0]), 'diagnostic_P_O': float(plq_O[0]),
            'resp_q_H': float(plq_H[2]), 'resp_q_O': float(plq_O[2]),
            'residual': float(residuals),
            'rank': int(rank),
            'fit_dz_shift': float(dz_shift_fit),
            'interaction_energy_range': [float(E_int.min()), float(E_int.max())],
            'fit_Eint_max': float(Eint_max),
            'fit_n_keep': int(n_keep),
            'z_clamp_fraction': float(n_clamp / (n_configs * n_atoms)),
            'xy_wrap_counts': [int(n_wrapx), int(n_wrapy)],
            'datasets': [name for name, _ in dsets]
        }
        fit_path = os.path.join(out_dir, 'fitted_PLQ_coeffs.json')
        with open(fit_path, 'w') as f:
            json.dump(fit_coeffs, f, indent=2)
        print(f'Saved fitted coefficients to: {fit_path}')

        # Generate comparison plots with fitted coefficients
        print('\nGenerating comparison plot with fitted coefficients...')
        diag_path_fit = os.path.join(out_dir, 'total_potential_H_atom_fitted.png')
        su.plot_alignment_summary(gridff, g0, dg, sub_apos, sub_enames,
                                  save_path=diag_path_fit, iz_top=iz_top, iy_center=64,
                                  z_atom_range=5.0, mol_apos=mol_apos_subset, mol_enames=names,
                                  plq_coeffs=(c[0], 0.0, qH_fit),  # fitted P_H, q_H; London=0
                                  plq_coeffs2=(c[1], 0.0, qO_fit),  # fitted P_O, q_O; London=0
                                  zmin_offset=2.0, z_ylim=(-2.0, 2.0), vmin_vmax_xz=(-0.5, 0.5))
        print('Saved fitted H atom plot:', diag_path_fit)

        diag_path_fit_o = os.path.join(out_dir, 'total_potential_O_atom_fitted.png')
        su.plot_alignment_summary(gridff, g0, dg, sub_apos, sub_enames,
                                  save_path=diag_path_fit_o, iz_top=iz_top, iy_center=64,
                                  z_atom_range=5.0, mol_apos=mol_apos_subset, mol_enames=names,
                                  plq_coeffs=(c[1], 0.0, qO_fit),
                                  plq_coeffs2=(c[0], 0.0, qH_fit),
                                  zmin_offset=2.0, z_ylim=(-2.0, 2.0), vmin_vmax_xz=(-0.5, 0.5))
        print('Saved fitted O atom plot:', diag_path_fit_o)

        # === 1D scan diagnostic: DFT vs fitted interaction energies ===
        # Panels: H2O-H/Na, H2O-H/Cl, H2O-O/Na, H2O-O/Cl
        # Lines within each panel: one per orientxy (Na, Cl, hollow), tilt=0 only
        print('\n=== 1D scan diagnostic: DFT vs fitted ===')
        import matplotlib.pyplot as plt
        panel_mols = ['H2O-H', 'H2O-H', 'H2O-O', 'H2O-O']
        panel_sites = ['Na', 'Cl', 'Na', 'Cl']   # substrate site (ix=0,iy=0 -> Cl; ix=6,iy=0 -> Na)
        panel_ix    = [6, 0, 6, 0]
        panel_iy    = [0, 0, 0, 0]
        orient_colors = {'Na': 'r', 'Cl': 'g', 'hollow': 'b'}
        scan_tilt = 0

        P_H0 = float(plq_H[0]); P_O0 = float(plq_O[0])
        def _plot_scan(fig_path, P_H_use, P_O_use, title_tag):
            fig, axs = plt.subplots(2, 2, figsize=(12, 10))
            axs = axs.flatten()
            for pidx in range(4):
                ax = axs[pidx]
                axr = ax.twinx()
                mol_p = panel_mols[pidx]
                ix_p = panel_ix[pidx]
                iy_p = panel_iy[pidx]
                site_label = panel_sites[pidx]
                for orient_name, col in orient_colors.items():
                    mask = (mol_tag == mol_p) & (ix == ix_p) & (iy == iy_p) & (orientxy == orient_name) & (orientz == scan_tilt)
                    if np.count_nonzero(mask) == 0: continue
                    z_sel = zdist[mask]
                    E_sel = energies[mask]
                    c_sel = coords[mask]
                    zmax = float(np.max(z_sel))
                    mask_zmax = np.isclose(z_sel, zmax, rtol=0.0, atol=1e-6)
                    E0 = float(np.mean(E_sel[mask_zmax]))
                    E_int_dft = E_sel - E0

                    n_frames = len(z_sel)
                    Em = np.zeros(n_frames, dtype=np.float64)
                    for i in range(n_frames):
                        rhoH_i = 0.0; rhoO_i = 0.0; phiH_i = 0.0; phiO_i = 0.0
                        for j, ename in enumerate(enames):
                            pos = c_sel[i, j].copy(); pos[2] += dz_shift_fit
                            vals, _, _, _ = sample_gridff(pos, gridff, g0, dg)
                            if ename == 'H':
                                rhoH_i += vals[idx_pauli]
                                phiH_i += vals[idx_coulomb]
                            elif ename == 'O':
                                rhoO_i += vals[idx_pauli]
                                phiO_i += vals[idx_coulomb]
                        Ec = qH_fit*(phiH_i - 2.0*phiO_i)
                        Em[i] = P_H_use*rhoH_i + P_O_use*rhoO_i + Ec

                    Em0 = float(np.mean(Em[mask_zmax]))
                    E_int_m = Em - Em0

                    sort_idx = np.argsort(z_sel)
                    z_s = z_sel[sort_idx]
                    Edft_s = E_int_dft[sort_idx]
                    Em_s = E_int_m[sort_idx]
                    ax.plot(z_s, Edft_s, ls=':',  lw=1.5, color=col, label=f'{orient_name} (DFT)')
                    ax.plot(z_s, Em_s,  ls='-',  lw=0.6, color=col, label=f'{orient_name} ({title_tag})')
                    axr.plot(z_s, Em_s-Edft_s, ls='-', lw=0.6, color=col, alpha=0.35)

                ax.axhline(0, color='k', lw=0.5, alpha=0.5)
                axr.axhline(0, color='k', lw=0.5, alpha=0.15)
                ax.set_title(f'{mol_p} over {site_label} (ix={ix_p} iy={iy_p}, tilt={scan_tilt})')
                ax.set_xlabel('zdist [A]')
                ax.set_ylabel('E_int [eV]')
                axr.set_ylabel('residual [eV]')
                axr.set_ylim(-0.5, 0.5)
                ax.set_ylim(-0.5, 0.5)
                ax.legend()
                ax.grid(True, alpha=0.3)

            cap = f"{title_tag}: P_H={P_H_use:.3f} P_O={P_O_use:.3f} q_H={qH_fit:.3f} q_O={qO_fit:.3f} dz={dz_shift_fit:.3f} | Eint_max={Eint_max:.2f} kT={kT_weight:.2f} tilt={scan_tilt}"
            fig.suptitle(cap, fontsize=10)
            plt.tight_layout(rect=(0, 0, 1, 0.96))
            plt.savefig(fig_path, dpi=150)
            plt.close(fig)

        scan_path_def = os.path.join(out_dir, 'scan_DFT_vs_default.png')
        _plot_scan(scan_path_def, P_H0, P_O0, 'default')
        print('Saved 1D scan diagnostic (default):', scan_path_def)

        scan_path_fit = os.path.join(out_dir, 'scan_DFT_vs_fitted.png')
        _plot_scan(scan_path_fit, float(c[0]), float(c[1]), 'fitted')
        print('Saved 1D scan diagnostic (fitted):', scan_path_fit)

    print('\nDone. GridFF saved to:')
    print(f'  {gridff_path}')
    print(f'  {meta_path}')
    print('Diagnostic plots saved to:')
    print(f'  {out_dir}/')


if __name__ == '__main__':
    main()
