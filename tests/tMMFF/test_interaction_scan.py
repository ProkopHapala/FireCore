#!/usr/bin/env python3
"""
test_interaction_scan.py - Test molecule-substrate interaction scanning
Test cases: H2O on NaCl and PTCDA on NaCl

Uses MMparams/AtomTypes.dat for proper MMFF REQ parameter initialization.

Demonstrates:
1. Z-approach curve (1D) with H2O
2. Lateral (x,y) energy map (2D)
3. Rotation scan (1D)
4. SLERP path scan
5. Constrained relaxation scan
6. PTCDA on NaCl (larger molecule with proper atom typing)
"""

import sys, os
import numpy as np

# Ensure pyBall is importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from pyBall.OCL.InteractionEnergy import InteractionScanner, load_xyz_with_REQs
from pyBall.OCL import ScanUtils

# ======== Paths ========
XYZ_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'cpp', 'common_resources', 'xyz')
OUT_DIR = os.path.join(os.path.dirname(__file__), 'output_interaction_scan')


def run_macro_reference_scan():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    os.makedirs(OUT_DIR, exist_ok=True)
    mol_file = os.path.join(XYZ_DIR, 'H2O_O.xyz')
    sub_file = os.path.join(XYZ_DIR, 'NaCl_8x8_L3.xyz')
    scanner = InteractionScanner(nloc=32)
    scanner.load_molecule_xyz(mol_file)
    scanner.load_substrate_xyz(sub_file)
    scanner.enable_LJ = True
    scanner.enable_Coulomb = True
    scanner.enable_HBond = False
    scanner.enable_Morse = True
    scanner.enable_macro = True
    scanner.nPBC[:] = (4, 4, 0)
    scanner._update_macro_from_substrate()
    p0 = np.array([1.0, 1.0, 1.0], dtype=float)
    step_x = np.array([4.0, 0.0, 0.0], dtype=float)
    step_y = np.array([0.0, 4.0, 0.0], dtype=float)
    ns = 9
    tx = [p0 + i * step_x for i in range(ns)]
    ty = [p0 + i * step_y for i in range(ns)]
    transforms_x = ScanUtils.pack_transforms([np.eye(3)] * len(tx), tx)
    transforms_y = ScanUtils.pack_transforms([np.eye(3)] * len(ty), ty)
    Ex = np.array(scanner.evaluate(transforms_x)['total'], dtype=float)
    Ey = np.array(scanner.evaluate(transforms_y)['total'], dtype=float)
    dEx = Ex - Ex[0]
    dEy = Ey - Ey[0]
    labs = np.arange(ns)
    print('Equivalent-site line errors along x [eV]:')
    for i, de in zip(labs, dEx):
        print(f'  ix={i:2d}  dE={de:+.6e}')
    print('Equivalent-site line errors along y [eV]:')
    for i, de in zip(labs, dEy):
        print(f'  iy={i:2d}  dE={de:+.6e}')
    fig, ax = plt.subplots(1, 1, figsize=(8, 4.5))
    ax.plot(labs, dEx, 'o-', lw=1.5, label='step (4,0) Å')
    ax.plot(labs, dEy, 's-', lw=1.5, label='step (0,4) Å')
    ax.axhline(0.0, c='k', lw=0.8)
    ax.set_xlabel('Equivalent-site index along a+b')
    ax.set_ylabel('ΔE [eV]')
    ax.set_title('Equivalent-site energy error, H2O_O on NaCl8x8, z=1.0 Å')
    ax.grid(True, alpha=0.3)
    ax.legend(frameon=False)
    line_png = os.path.join(OUT_DIR, 'H2O_NaCl8x8_equiv_line_macro.png')
    fig.tight_layout()
    fig.savefig(line_png, dpi=180)
    plt.close(fig)

    res_xy = scanner.scan_lateral(z=1.0, x_range=(-5.0, 35.0), y_range=(-5.0, 35.0), nx=161, ny=161)
    E2d = np.array(res_xy['total'], dtype=float).reshape(res_xy['scan_info']['shape'])
    xs = res_xy['scan_info']['x']
    ys = res_xy['scan_info']['y']
    vabs = np.nanmax(np.abs(E2d - np.nanmean(E2d)))
    if not np.isfinite(vabs) or (vabs <= 0.0):
        raise ValueError(f'Invalid XY energy span {vabs}')
    fig, ax = plt.subplots(1, 1, figsize=(8, 7))
    im = ax.imshow((E2d - np.mean(E2d)).T, extent=[xs[0], xs[-1], ys[0], ys[-1]], origin='lower', aspect='equal', cmap='bwr', vmin=-vabs, vmax=vabs)
    ax.set_xlabel('x [Å]')
    ax.set_ylabel('y [Å]')
    ax.set_title('H2O_O on NaCl8x8 XY scan at z=1.0 Å (macro corrected, mean-shifted)')
    plt.colorbar(im, ax=ax, label='E - <E> [eV]')
    xy_png = os.path.join(OUT_DIR, 'H2O_NaCl8x8_XYscan_-5_35_npbc4_macro.png')
    fig.tight_layout()
    fig.savefig(xy_png, dpi=180)
    plt.close(fig)
    return {'line_png': line_png, 'xy_png': xy_png, 'dEx': dEx, 'dEy': dEy, 'scanner': scanner}


def test_h2o_on_nacl():
    """Test H2O on NaCl - basic small system."""
    print("=" * 60)
    print("Test 1: H2O on NaCl (MMFF REQ params from AtomTypes.dat)")
    print("=" * 60)

    mol_file = os.path.join(XYZ_DIR, 'H2O_O.xyz')
    sub_file = os.path.join(XYZ_DIR, 'NaCl_1x1_L3.xyz')

    scanner = InteractionScanner(nloc=32)
    mol_pos, mol_REQs, mol_names = scanner.load_molecule_xyz(mol_file)
    sub_pos, sub_REQs, sub_names = scanner.load_substrate_xyz(sub_file)

    print(f"Molecule: {len(mol_names)} atoms: {mol_names}")
    print(f"  positions:\n{mol_pos}")
    print(f"  REQs (R, sqrt(E), Q, H):\n{mol_REQs}")
    print(f"Substrate: {len(sub_names)} atoms: {sub_names}")
    print(f"  REQs (R, sqrt(E), Q, H):\n{sub_REQs}")

    scanner.enable_LJ      = True
    scanner.enable_Coulomb  = True
    scanner.enable_HBond    = True
    scanner.enable_Morse    = False

    # -------- PBC invariance sanity check --------
    # With finite image sums (nPBC), invariance under translation by lattice vectors is only guaranteed
    # if we wrap the molecule translation into the primary cell (InteractionScanner.wrap_PBC=True).
    print("\n--- 0. PBC invariance check (E(t) == E(t+a) == E(t+b)) ---")
    assert scanner.lvec is not None, "NaCl_1x1_L3.xyz must provide lattice vectors in comment line"
    a = np.array(scanner.lvec[0], dtype=float)
    b = np.array(scanner.lvec[1], dtype=float)
    pos = np.array([0.1, 0.2, 3.0], dtype=float)
    E0 = scanner.evaluate_single(pos=pos)['total']
    Ea = scanner.evaluate_single(pos=pos + a)['total']
    Eb = scanner.evaluate_single(pos=pos + b)['total']
    dEa = float(Ea - E0)
    dEb = float(Eb - E0)
    print(f"  E0={E0:+.6e}  Ea={Ea:+.6e}  dEa={dEa:+.3e}  Eb={Eb:+.6e}  dEb={dEb:+.3e}")
    assert abs(dEa) < 1e-5 and abs(dEb) < 1e-5, "PBC invariance failed; check wrap_PBC or kernel PBC loops"

    # -------- 1. Z-approach curve --------
    print("\n--- 1. Z-approach scan ---")
    res_z = scanner.scan_z(pos_xy=(0.0, 0.0), z_range=(1.5, 8.0), nz=60)
    zs = res_z['scan_info']['z']
    imin = np.argmin(res_z['total'])
    print(f"  Min total energy: {res_z['total'][imin]:.4f} eV at z = {zs[imin]:.3f} Ang")
    print(f"    LJ:      {res_z['LJ'][imin]:.4f}")
    print(f"    Coulomb: {res_z['Coulomb'][imin]:.4f}")
    print(f"    HBond:   {res_z['HBond'][imin]:.4f}")

    # -------- 2. Lateral 2D scan --------
    print("\n--- 2. Lateral (x,y) scan at z=3.0 ---")
    res_xy = scanner.scan_lateral(z=3.0, x_range=(0, 4), y_range=(0, 4), nx=30, ny=30)
    E2d = res_xy['total'].reshape(res_xy['scan_info']['shape'])
    imin2d = np.unravel_index(np.argmin(E2d), E2d.shape)
    xs, ys = res_xy['scan_info']['x'], res_xy['scan_info']['y']
    print(f"  Min energy: {E2d[imin2d]:.4f} eV at x={xs[imin2d[0]]:.2f}, y={ys[imin2d[1]]:.2f}")

    # -------- 3. Rotation scan --------
    print("\n--- 3. Rotation scan (z-axis) at pos=(0,0,3) ---")
    res_rot = scanner.scan_rotation(pos=(0.0, 0.0, 3.0), axis=(0,0,1), nrot=36)
    angles = res_rot['scan_info']['angles']
    imin_rot = np.argmin(res_rot['total'])
    print(f"  Min energy: {res_rot['total'][imin_rot]:.4f} eV at angle={angles[imin_rot]:.1f} deg")

    # -------- 4. SLERP path --------
    print("\n--- 4. SLERP path scan ---")
    q0 = ScanUtils.quat_from_axis_angle([0,0,1], 0)
    q1 = ScanUtils.quat_from_axis_angle([0,0,1], np.pi)
    res_sl = scanner.scan_slerp(q0, q1, t0=[0,0,3], t1=[4,0,3], npts=40)
    imin_sl = np.argmin(res_sl['total'])
    print(f"  Min energy along SLERP path: {res_sl['total'][imin_sl]:.4f} eV at t={res_sl['scan_info']['t'][imin_sl]:.3f}")

    # -------- 5. Constrained relaxation --------
    print("\n--- 5. Z-approach with constrained relaxation ---")
    scanner.spring_k    = np.float32(5.0)
    scanner.relax_dt    = np.float32(0.005)
    scanner.relax_nsteps = 100
    res_relax = scanner.scan_z(pos_xy=(0.0, 0.0), z_range=(2.0, 6.0), nz=40, relax=True)
    zs_r = res_relax['scan_info']['z']
    imin_r = np.argmin(res_relax['total'])
    print(f"  Min relaxed energy: {res_relax['total'][imin_r]:.4f} eV at z = {zs_r[imin_r]:.3f}")

    return scanner, res_z, res_xy, res_rot, res_sl, res_relax


def test_ptcda_on_nacl():
    """Test PTCDA on NaCl - larger molecule with proper atom typing."""
    print("\n" + "=" * 60)
    print("Test 2: PTCDA on NaCl (38 atoms, proper atom typing)")
    print("=" * 60)

    mol_file = os.path.join(XYZ_DIR, 'PTCDA.xyz')
    sub_file = os.path.join(XYZ_DIR, 'NaCl_8x8_L3.xyz')  # 384 atoms

    # PTCDA atom type mapping:
    # - All C are aromatic/conjugated -> C_R
    # - Bridge O (ether in anhydride) -> O_3
    # - Carbonyl O -> O_2
    # - H bonded to aromatic C -> H
    # But since we load from xyz (only element names), we need a type_map
    # The xyz has 24 C, 6 O, 8 H
    # We map C->C_R for all carbons (simplification; some are C_2 carbonyl)
    ptcda_type_map = {'C': 'C_R', 'O': 'O_2'}  # O_2 works for both ether and carbonyl vdW

    scanner = InteractionScanner(nloc=32)
    mol_pos, mol_REQs, mol_names = scanner.load_molecule_xyz(mol_file, type_map=ptcda_type_map)
    sub_pos, sub_REQs, sub_names = scanner.load_substrate_xyz(sub_file)

    print(f"Molecule: {len(mol_names)} atoms")
    print(f"  Elements: {set(mol_names)}")
    print(f"  REQs sample (first 5):\n{mol_REQs[:5]}")
    print(f"Substrate: {len(sub_names)} atoms")
    print(f"  Elements: {set(sub_names)}")

    scanner.enable_LJ      = True
    scanner.enable_Coulomb  = True
    scanner.enable_HBond    = False
    scanner.enable_Morse    = False

    # -------- Z-approach --------
    print("\n--- PTCDA Z-approach scan ---")
    res_z = scanner.scan_z(pos_xy=(0.0, 0.0), z_range=(2.5, 10.0), nz=80)
    zs = res_z['scan_info']['z']
    imin = np.argmin(res_z['total'])
    print(f"  Min total energy: {res_z['total'][imin]:.4f} eV at z = {zs[imin]:.3f} Ang")
    print(f"    LJ:      {res_z['LJ'][imin]:.4f}")
    print(f"    Coulomb: {res_z['Coulomb'][imin]:.4f}")

    # -------- Lateral 2D scan --------
    print("\n--- PTCDA Lateral scan at z=3.5 ---")
    res_xy = scanner.scan_lateral(z=3.5, x_range=(-5, 5), y_range=(-5, 5), nx=40, ny=40)
    E2d = res_xy['total'].reshape(res_xy['scan_info']['shape'])
    imin2d = np.unravel_index(np.argmin(E2d), E2d.shape)
    xs, ys = res_xy['scan_info']['x'], res_xy['scan_info']['y']
    print(f"  Min energy: {E2d[imin2d]:.4f} eV at x={xs[imin2d[0]]:.2f}, y={ys[imin2d[1]]:.2f}")

    # -------- Rotation scan --------
    print("\n--- PTCDA Rotation scan at z=3.5 ---")
    res_rot = scanner.scan_rotation(pos=(0, 0, 3.5), axis=(0,0,1), nrot=72)
    angles = res_rot['scan_info']['angles']
    imin_rot = np.argmin(res_rot['total'])
    print(f"  Min energy: {res_rot['total'][imin_rot]:.4f} eV at angle={angles[imin_rot]:.1f} deg")
    print(f"  Energy range: [{res_rot['total'].min():.4f}, {res_rot['total'].max():.4f}] eV")

    # -------- Single pose evaluation --------
    print("\n--- PTCDA single pose ---")
    e_single = scanner.evaluate_single(pos=(0, 0, 3.5))
    print(f"  E_total={e_single['total']:.4f}  E_LJ={e_single['LJ']:.4f}  E_Coul={e_single['Coulomb']:.4f}")

    return scanner, res_z, res_xy, res_rot


def main():
    macro_info = run_macro_reference_scan()
    print(f"Macro 1D equivalent-site plot saved to: {macro_info['line_png']}")
    print(f"Macro 2D XY scan saved to: {macro_info['xy_png']}")

    # Run H2O test
    scanner_h2o, res_z_h2o, res_xy_h2o, res_rot_h2o, res_sl_h2o, res_relax_h2o = test_h2o_on_nacl()

    # Run PTCDA test
    scanner_ptcda, res_z_ptcda, res_xy_ptcda, res_rot_ptcda = test_ptcda_on_nacl()

    # -------- Plotting --------
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(2, 4, figsize=(20, 10))
        fig.suptitle('Molecule-Substrate Interaction Scans (MMFF REQ params)', fontsize=14)

        # --- H2O plots ---
        zs = res_z_h2o['scan_info']['z']

        ax = axes[0, 0]
        ax.plot(zs, res_z_h2o['total'], 'k-', lw=2, label='Total')
        ax.plot(zs, res_z_h2o['LJ'], 'b--', label='LJ')
        ax.plot(zs, res_z_h2o['Coulomb'], 'r--', label='Coulomb')
        ax.plot(zs, res_z_h2o['HBond'], 'g--', label='HBond')
        ax.set_xlabel('z [Å]'); ax.set_ylabel('E [eV]')
        ax.set_title('H2O Z-approach')
        ax.legend(fontsize=7); ax.grid(True, alpha=0.3)
        ax.set_ylim(max(-5, res_z_h2o['total'].min()-0.5), min(5, res_z_h2o['total'].max()+0.5))

        ax = axes[0, 1]
        E2d = res_xy_h2o['total'].reshape(res_xy_h2o['scan_info']['shape'])
        xs, ys = res_xy_h2o['scan_info']['x'], res_xy_h2o['scan_info']['y']
        im = ax.imshow(E2d.T, extent=[xs[0], xs[-1], ys[0], ys[-1]], origin='lower', aspect='equal', cmap='RdBu_r')
        ax.set_xlabel('x [Å]'); ax.set_ylabel('y [Å]')
        ax.set_title('H2O Lateral (z=3.0)')
        plt.colorbar(im, ax=ax, label='E [eV]')

        ax = axes[0, 2]
        angles = res_rot_h2o['scan_info']['angles']
        ax.plot(angles, res_rot_h2o['total'], 'k-', lw=2)
        ax.set_xlabel('angle [°]'); ax.set_ylabel('E [eV]')
        ax.set_title('H2O Rotation')
        ax.grid(True, alpha=0.3)

        ax = axes[0, 3]
        zs_r = res_relax_h2o['scan_info']['z']
        zs_rigid = res_z_h2o['scan_info']['z']
        ax.plot(zs_rigid, res_z_h2o['total'], 'b-', lw=2, label='Rigid')
        ax.plot(zs_r, res_relax_h2o['total'], 'r-', lw=2, label='Relaxed')
        ax.set_xlabel('z [Å]'); ax.set_ylabel('E [eV]')
        ax.set_title('H2O Rigid vs Relaxed')
        ax.legend(fontsize=8); ax.grid(True, alpha=0.3)
        ylo = max(-5, min(res_z_h2o['total'].min(), res_relax_h2o['total'].min()) - 0.5)
        yhi = min(5, max(res_z_h2o['total'].max(), res_relax_h2o['total'].max()) + 0.5)
        ax.set_ylim(ylo, yhi)

        # --- PTCDA plots ---
        zs_p = res_z_ptcda['scan_info']['z']

        ax = axes[1, 0]
        ax.plot(zs_p, res_z_ptcda['total'], 'k-', lw=2, label='Total')
        ax.plot(zs_p, res_z_ptcda['LJ'], 'b--', label='LJ')
        ax.plot(zs_p, res_z_ptcda['Coulomb'], 'r--', label='Coulomb')
        ax.set_xlabel('z [Å]'); ax.set_ylabel('E [eV]')
        ax.set_title('PTCDA Z-approach')
        ax.legend(fontsize=7); ax.grid(True, alpha=0.3)
        ax.set_ylim(max(-10, res_z_ptcda['total'].min()-1), min(10, res_z_ptcda['total'].max()+1))

        ax = axes[1, 1]
        E2d_p = res_xy_ptcda['total'].reshape(res_xy_ptcda['scan_info']['shape'])
        xs_p, ys_p = res_xy_ptcda['scan_info']['x'], res_xy_ptcda['scan_info']['y']
        vmin, vmax = np.percentile(E2d_p[np.isfinite(E2d_p)], [2, 98])
        im = ax.imshow(E2d_p.T, extent=[xs_p[0], xs_p[-1], ys_p[0], ys_p[-1]], origin='lower', aspect='equal', cmap='RdBu_r', vmin=vmin, vmax=vmax)
        ax.set_xlabel('x [Å]'); ax.set_ylabel('y [Å]')
        ax.set_title('PTCDA Lateral (z=3.5)')
        plt.colorbar(im, ax=ax, label='E [eV]')

        ax = axes[1, 2]
        angles_p = res_rot_ptcda['scan_info']['angles']
        ax.plot(angles_p, res_rot_ptcda['total'], 'k-', lw=2)
        ax.set_xlabel('angle [°]'); ax.set_ylabel('E [eV]')
        ax.set_title('PTCDA Rotation (z=3.5)')
        ax.grid(True, alpha=0.3)

        ax = axes[1, 3]
        ax.plot(res_sl_h2o['scan_info']['t'], res_sl_h2o['total'], 'k-', lw=2, label='Total')
        ax.plot(res_sl_h2o['scan_info']['t'], res_sl_h2o['LJ'], 'b--', label='LJ')
        ax.plot(res_sl_h2o['scan_info']['t'], res_sl_h2o['Coulomb'], 'r--', label='Coulomb')
        ax.set_xlabel('Path t'); ax.set_ylabel('E [eV]')
        ax.set_title('H2O SLERP path')
        ax.legend(fontsize=7); ax.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig('interaction_scan_results.png', dpi=150)
        print(f"\nPlot saved to: interaction_scan_results.png")
    except ImportError:
        print("\nMatplotlib not available, skipping plots.")

    print("\n" + "=" * 60)
    print("All scans completed successfully!")
    print("=" * 60)


if __name__ == '__main__':
    main()
