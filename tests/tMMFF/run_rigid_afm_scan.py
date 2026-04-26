#!/usr/bin/env python3

import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(BASE)

from pyBall.OCL.RigidBodyAFM import RigidBodyAFM, sample_gridff_single_atom, plot_gridff_diagnostics

def read_xyz(fname):
    with open(fname, 'r') as f:
        lines = f.readlines()
    natoms = int(lines[0].strip())
    enames = []
    apos = []
    for line in lines[2:2+natoms]:
        parts = line.split()
        if len(parts) >= 4:
            enames.append(parts[0])
            apos.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return np.array(apos, dtype=np.float32), enames

def write_xyz_frame(fout, enames, pos3, comment=""):
    fout.write(f"{len(enames)}\n")
    fout.write(f"{comment}\n")
    for e, p in zip(enames, pos3):
        fout.write(f"{e:3s} {p[0]:12.6f} {p[1]:12.6f} {p[2]:12.6f}\n")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mol', default='xyz/PTCDA.xyz')
    ap.add_argument('--sub', default='CaF2_6L_Ni3_rect_nx2_nz1_L2_top')
    ap.add_argument('--grid', default='Bspline_PLQd.npy')
    ap.add_argument('--outdir', default='afm_results', help='Directory to store results')
    ap.add_argument('--anchor', type=int, default=28, help='Anchor atom index')
    ap.add_argument('--kanchor', type=float, default=20.0, help='Anchor stiffness')
    ap.add_argument('--nscan', type=int, default=100)
    ap.add_argument('--x0', type=float, default=0.0)
    ap.add_argument('--x1', type=float, default=50.0)
    ap.add_argument('--y0', type=float, default=2.0)
    ap.add_argument('--z0', type=float, default=6.0)
    # 2D scan parameters
    ap.add_argument('--nx', type=int, default=40, help='Number of X pixels for 2D scan')
    ap.add_argument('--ny', type=int, default=40, help='Number of Y pixels for 2D scan')
    ap.add_argument('--x0_2d', type=float, default=0.0, help='X start for 2D scan')
    ap.add_argument('--x1_2d', type=float, default=None, help='X end for 2D scan (default: lattice vector ax)')
    ap.add_argument('--y0_2d', type=float, default=0.0, help='Y start for 2D scan')
    ap.add_argument('--y1_2d', type=float, default=None, help='Y end for 2D scan (default: lattice vector ay)')
    ap.add_argument('--dt', type=float, default=0.05)
    ap.add_argument('--nsteps', type=int, default=1000)
    ap.add_argument('--fconv', type=float, default=1e-3)
    ap.add_argument('--tconv', type=float, default=1e-3)
    ap.add_argument('--test_z', type=float, default=5.0, help='Z height for the test grid scan')
    ap.add_argument('--diag', action='store_true', help='Plot GridFF diagnostics (raw grids)')
    ap.add_argument('--test_grid', action='store_true', help='Run 2D single-atom alignment tests')
    ap.add_argument('--scan1d', action='store_true', help='Run 1D relaxed AFM scan (whole molecule)')
    ap.add_argument('--scan2d', action='store_true', help='Run 2D relaxed AFM scan (parallel molecule relaxation)')
    ap.add_argument('--all', action='store_true', help='Run everything')
    ap.add_argument('--debug', action='store_true', help='Enable diagnostic plots (diag and test_grid) with scan2d')
    args = ap.parse_args()

    # Determine tasks
    do_diag = args.diag or args.all or (args.debug and args.scan2d)
    do_test_grid = args.test_grid or args.all or (args.debug and args.scan2d)
    do_scan1d = args.scan1d or args.all or (not (args.diag or args.test_grid or args.scan2d))
    do_scan2d = args.scan2d or args.all

    # Setup paths
    data_root = os.path.join(BASE, 'cpp', 'common_resources')
    mol_path = os.path.join(data_root, args.mol)
    sub_dir = os.path.join(data_root, args.sub)
    if not os.path.exists(sub_dir): sub_dir = data_root
    sub_xyz = os.path.join(sub_dir, args.sub + '.xyz')
    if not os.path.exists(sub_xyz):
        sub_xyz = os.path.join(data_root, 'Substrates', 'generated_rect', args.sub + '.xyz')
        if not os.path.exists(sub_xyz): sub_xyz = os.path.join(data_root, 'xyz', args.sub + '.xyz')
    gridff_path = os.path.join(sub_dir, args.grid)

    if not os.path.exists(args.outdir): os.makedirs(args.outdir)

    print(f"--- AFM Run Configuration ---")
    print(f"Molecule  : {mol_path}")
    print(f"Substrate : {sub_xyz}")
    print(f"Grid      : {gridff_path}")
    print(f"OutDir    : {args.outdir}")

    # Load substrate atoms and lattice
    sub_apos, sub_enames = read_xyz(sub_xyz)
    with open(sub_xyz, 'r') as f:
        lines = f.readlines()
        comment = lines[1]
        lvec = None
        for pref in ["lvec:", "lvs"]:
            if pref in comment:
                vals = [float(x) for x in comment.split(pref)[1].split()[:9]]
                lvec = np.array(vals).reshape(3,3)
                break
    if lvec is None: raise ValueError("Lattice vectors not found in substrate XYZ comment line")
    ax, ay, az = np.linalg.norm(lvec, axis=1)
    print(f"Lattice detected: {ax:.3f}, {ay:.3f}, {az:.3f}")

    # Task 1: Raw Grid Diagnostics
    if do_diag:
        print("\n[1/4] Plotting GridFF Diagnostics (Raw Grids)...")
        plot_gridff_diagnostics(gridff_path, sub_xyz, ax, ay, az, save_path=os.path.join(args.outdir, 'diag_raw_grids.png'))

    # Task 2: Single Atom 2D Alignment Scan
    if do_test_grid:
        print(f"\n[2/4] Running Single-Atom 2D Alignment Test at Z={args.test_z}...")
        nx, ny = 100, 100
        xs = np.linspace(0, ax, nx)
        ys = np.linspace(0, ay, ny)
        X, Y = np.meshgrid(xs, ys)
        scan_pos = np.zeros((nx * ny, 3), dtype=np.float32)
        scan_pos[:, 0] = X.flatten(); scan_pos[:, 1] = Y.flatten(); scan_pos[:, 2] = args.test_z
        
        req_neutral = (1.6612, 0.0091, 0.0, 0.0)
        req_charged = (1.6612, 0.0091, -0.2, 0.0)
        
        _, E_neu = sample_gridff_single_atom(scan_pos, gridff_path, sub_xyz, atom_req=req_neutral)
        _, E_char = sample_gridff_single_atom(scan_pos, gridff_path, sub_xyz, atom_req=req_charged)
        
        plt.figure(figsize=(12, 5))
        plt.subplot(1, 2, 1); plt.imshow(E_neu.reshape(ny, nx), extent=[0, ax, 0, ay], origin='lower', cmap='bwr')
        plt.scatter(sub_apos[:,0], sub_apos[:,1], c='k', s=5, alpha=0.3); plt.title("Neutral O"); plt.colorbar()
        plt.subplot(1, 2, 2); plt.imshow(E_char.reshape(ny, nx), extent=[0, ax, 0, ay], origin='lower', cmap='bwr')
        plt.scatter(sub_apos[:,0], sub_apos[:,1], c='k', s=5, alpha=0.3); plt.title("Charged O (Q=-0.2)"); plt.colorbar()
        plt.savefig(os.path.join(args.outdir, 'test_single_atom_2d.png'))

    # Task 3: 1D Relaxed AFM Scan (Trajectory)
    if do_scan1d:
        print("\n[3/4] Running 1D Relaxed AFM Scan (Whole Molecule)...")
        afm = RigidBodyAFM(mol_path, gridff_path, sub_xyz, anchor_idx=args.anchor, anchor_k=args.kanchor)
        mol_apos, mol_enames = read_xyz(mol_path)
        r_anchor = mol_apos[args.anchor]
        tip_pos = np.array([args.x0, args.y0, args.z0], dtype=np.float32)
        initial_com = tip_pos - r_anchor
        afm.prepare(n_bodies=1, initial_positions=np.array([initial_com], dtype=np.float32))
        
        out_xyz = os.path.join(args.outdir, 'scan_molecule.xyz')
        dx = (args.x1 - args.x0) / args.nscan
        
        trajectory = []
        with open(out_xyz, 'w') as fout:
            for i in range(args.nscan):
                afm.set_anchor_positions(tip_pos)
                outputs, converged = afm.relax_to_constraint(nsteps=args.nsteps, dt=args.dt, fconv=args.fconv, tconv=args.tconv)
                
                apos = outputs['atom_positions'][0]
                aforce = outputs['atom_force'][0]
                E = aforce[args.anchor, 3]
                
                print(f"Step {i:3d}/{args.nscan}: tip=({tip_pos[0]:.2f}, {tip_pos[1]:.2f}, {tip_pos[2]:.2f}) E={E:9.4f} conv={converged[0]}")
                
                trajectory.append({
                    'tip': tip_pos.copy(),
                    'apos': apos.copy(),
                    'E': E,
                    'F': aforce[args.anchor, :3].copy(),
                    'pos': outputs['pos'][0, :3].copy()
                })
                
                write_xyz_frame(fout, sub_enames + mol_enames, np.vstack([sub_apos, apos[:, :3]]), comment=f"frame {i} tip {tip_pos} E {E:.4f}")
                tip_pos[0] += dx

        # Plotting 1D results
        traj_data = {k: np.array([t[k] for t in trajectory]) for k in trajectory[0].keys()}
        
        plt.figure(figsize=(10, 6))
        plt.plot(traj_data['tip'][:, 0], traj_data['E'], label='Energy (eV)')
        plt.plot(traj_data['tip'][:, 0], traj_data['F'][:, 2] * 10, label='Force Z (eV/A x10)')
        plt.xlabel('Tip X (A)'); plt.ylabel('Energy / Force'); plt.legend(); plt.title('AFM 1D Scan Profile')
        plt.savefig(os.path.join(args.outdir, 'afm_scan_1d_profile.png'))

        plt.figure(figsize=(12, 5))
        plt.subplot(1, 2, 1); plt.plot(traj_data['pos'][:, 0], traj_data['pos'][:, 1], 'k-', label='COM')
        plt.plot(traj_data['apos'][:, args.anchor, 0], traj_data['apos'][:, args.anchor, 1], 'r--', label='Anchor')
        plt.scatter(sub_apos[:, 0], sub_apos[:, 1], c='gray', s=10, alpha=0.2); plt.title('XY Projection'); plt.legend(); plt.axis('equal')
        plt.subplot(1, 2, 2); plt.plot(traj_data['pos'][:, 0], traj_data['pos'][:, 2], 'k-', label='COM')
        plt.plot(traj_data['apos'][:, args.anchor, 0], traj_data['apos'][:, args.anchor, 2], 'r--', label='Anchor')
        plt.scatter(sub_apos[:, 0], sub_apos[:, 2], c='gray', s=10, alpha=0.2); plt.title('XZ Projection'); plt.legend(); plt.axis('equal')
        plt.savefig(os.path.join(args.outdir, 'afm_scan_1d_projections.png'))

    # Task 4: 2D Parallel AFM Scan (Images)
    if do_scan2d:
        print("\n[4/4] Running 2D Parallel AFM Scan (Full Relaxation per Pixel)...")
        nx, ny = args.nx, args.ny
        x0_2d, x1_2d = args.x0_2d, (args.x1_2d if args.x1_2d is not None else ax)
        y0_2d, y1_2d = args.y0_2d, (args.y1_2d if args.y1_2d is not None else ay)
        xs = np.linspace(x0_2d, x1_2d, nx); ys = np.linspace(y0_2d, y1_2d, ny)
        X, Y = np.meshgrid(xs, ys)
        n_bodies = nx * ny
        print(f"2D Scan: {nx}×{ny} = {n_bodies} pixels over X=[{x0_2d:.2f}, {x1_2d:.2f}] Y=[{y0_2d:.2f}, {y1_2d:.2f}] Å")
        
        mol_apos, mol_enames = read_xyz(mol_path)
        r_anchor = mol_apos[args.anchor]
        initial_coms = np.zeros((n_bodies, 3), dtype=np.float32)
        initial_coms[:, 0] = X.flatten() - r_anchor[0]
        initial_coms[:, 1] = Y.flatten() - r_anchor[1]
        initial_coms[:, 2] = args.z0 - r_anchor[2]
        
        afm = RigidBodyAFM(mol_path, gridff_path, sub_xyz, anchor_idx=args.anchor, anchor_k=args.kanchor)
        afm.prepare(n_bodies=n_bodies, initial_positions=initial_coms)
        afm.set_anchor_positions(np.vstack([X.flatten(), Y.flatten(), np.full(n_bodies, args.z0)]).T)
        
        print(f"Relaxing {n_bodies} molecules in parallel...")
        outputs, converged = afm.relax_to_constraint(nsteps=args.nsteps, dt=args.dt, fconv=args.fconv, tconv=args.tconv)
        
        # Diagnostics: check anchor displacement after relaxation
        apos_final = outputs['atom_positions']
        anchor_final = apos_final[:, args.anchor, :3]
        anchor_target = np.vstack([X.flatten(), Y.flatten(), np.full(n_bodies, args.z0)]).T
        anchor_disp = np.linalg.norm(anchor_final - anchor_target, axis=1)
        print(f"Anchor displacement after relaxation: mean={anchor_disp.mean():.4f} Å, max={anchor_disp.max():.4f} Å")
        print(f"Converged replicas: {np.sum(converged)}/{n_bodies}")
        
        force_z = outputs['atom_force'][:, args.anchor, 2].reshape(ny, nx)
        energy = outputs['atom_force'][:, args.anchor, 3].reshape(ny, nx)
        
        # Extract COM positions and quaternions
        pos = outputs['pos'][:, :3]
        quat = outputs['quats']
        
        # COM shift from initial positions
        com_shift_x = (pos[:, 0] - initial_coms[:, 0]).reshape(ny, nx)
        com_shift_y = (pos[:, 1] - initial_coms[:, 1]).reshape(ny, nx)
        com_shift_z = (pos[:, 2] - initial_coms[:, 2]).reshape(ny, nx)
        
        # Convert quaternion to rotation angle (in degrees)
        # quat = [qx, qy, qz, qw], angle = 2 * acos(|qw|)
        quat_w = np.abs(quat[:, 3])
        quat_w = np.clip(quat_w, -1.0, 1.0)
        rot_angle = 2.0 * np.arccos(quat_w) * 180.0 / np.pi
        rot_angle = rot_angle.reshape(ny, nx)
        
        # Rotation axis (normalized) for tilt direction
        quat_xyz = quat[:, :3]
        axis_norm = np.linalg.norm(quat_xyz, axis=1)
        axis_norm[axis_norm < 1e-10] = 1.0  # Avoid division by zero
        axis_x = (quat_xyz[:, 0] / axis_norm).reshape(ny, nx)
        axis_y = (quat_xyz[:, 1] / axis_norm).reshape(ny, nx)
        axis_z = (quat_xyz[:, 2] / axis_norm).reshape(ny, nx)
        
        # Create comprehensive figure
        fig = plt.figure(figsize=(18, 12))
        
        # Row 1: Force and Energy
        ax1 = plt.subplot(3, 4, 1)
        im1 = ax1.imshow(force_z, extent=[x0_2d, x1_2d, y0_2d, y1_2d], origin='lower', cmap='afmhot')
        ax1.set_title("Force Z (eV/Å)"); plt.colorbar(im1, ax=ax1)
        
        ax2 = plt.subplot(3, 4, 2)
        im2 = ax2.imshow(energy, extent=[x0_2d, x1_2d, y0_2d, y1_2d], origin='lower', cmap='viridis')
        ax2.set_title("Energy (eV)"); plt.colorbar(im2, ax=ax2)
        
        # Row 1: Convergence
        conv_map = converged.reshape(ny, nx)
        ax3 = plt.subplot(3, 4, 3)
        im3 = ax3.imshow(conv_map, extent=[x0_2d, x1_2d, y0_2d, y1_2d], origin='lower', cmap='RdYlGn')
        ax3.set_title("Converged"); plt.colorbar(im3, ax=ax3)
        
        # Row 2: COM shifts
        ax4 = plt.subplot(3, 4, 4)
        im4 = ax4.imshow(com_shift_z, extent=[x0_2d, x1_2d, y0_2d, y1_2d], origin='lower', cmap='coolwarm')
        ax4.set_title("COM Z Shift (Å)"); plt.colorbar(im4, ax=ax4)
        
        # XY COM shift as complex (HSV) plot
        ax5 = plt.subplot(3, 4, 5)
        # Combine X and Y shifts as complex number for HSV color
        com_shift_complex = com_shift_x + 1j * com_shift_y
        magnitude = np.abs(com_shift_complex)
        angle = np.angle(com_shift_complex)
        # Use HSV: Hue = angle, Saturation = 1, Value = magnitude (normalized)
        from matplotlib.colors import hsv_to_rgb
        hue = (angle + np.pi) / (2 * np.pi)  # Normalize to [0, 1]
        saturation = np.ones_like(hue)
        value = magnitude / (magnitude.max() + 1e-10)
        hsv = np.stack([hue, saturation, value], axis=-1)
        rgb = hsv_to_rgb(hsv)
        im5 = ax5.imshow(rgb, extent=[x0_2d, x1_2d, y0_2d, y1_2d], origin='lower')
        ax5.set_title("COM XY Shift (HSV: angle=magnitude)")
        # Add colorbar for magnitude
        sm = plt.cm.ScalarMappable(cmap='hsv', norm=plt.Normalize(vmin=0, vmax=magnitude.max()))
        sm.set_array([])
        plt.colorbar(sm, ax=ax5, label='Magnitude (Å)')
        
        ax6 = plt.subplot(3, 4, 6)
        im6 = ax6.imshow(com_shift_x, extent=[x0_2d, x1_2d, y0_2d, y1_2d], origin='lower', cmap='RdBu')
        ax6.set_title("COM X Shift (Å)"); plt.colorbar(im6, ax=ax6)
        
        ax7 = plt.subplot(3, 4, 7)
        im7 = ax7.imshow(com_shift_y, extent=[x0_2d, x1_2d, y0_2d, y1_2d], origin='lower', cmap='RdBu')
        ax7.set_title("COM Y Shift (Å)"); plt.colorbar(im7, ax=ax7)
        
        # Row 3: Quaternion tilt
        ax8 = plt.subplot(3, 4, 8)
        im8 = ax8.imshow(rot_angle, extent=[x0_2d, x1_2d, y0_2d, y1_2d], origin='lower', cmap='plasma')
        ax8.set_title("Rotation Angle (deg)"); plt.colorbar(im8, ax=ax8)
        
        ax9 = plt.subplot(3, 4, 9)
        im9 = ax9.imshow(axis_x, extent=[x0_2d, x1_2d, y0_2d, y1_2d], origin='lower', cmap='RdBu', vmin=-1, vmax=1)
        ax9.set_title("Tilt Axis X"); plt.colorbar(im9, ax=ax9)
        
        ax10 = plt.subplot(3, 4, 10)
        im10 = ax10.imshow(axis_y, extent=[x0_2d, x1_2d, y0_2d, y1_2d], origin='lower', cmap='RdBu', vmin=-1, vmax=1)
        ax10.set_title("Tilt Axis Y"); plt.colorbar(im10, ax=ax10)
        
        ax11 = plt.subplot(3, 4, 11)
        im11 = ax11.imshow(axis_z, extent=[x0_2d, x1_2d, y0_2d, y1_2d], origin='lower', cmap='RdBu', vmin=-1, vmax=1)
        ax11.set_title("Tilt Axis Z"); plt.colorbar(im11, ax=ax11)
        
        plt.tight_layout()
        plt.savefig(os.path.join(args.outdir, 'scan_molecule_2d_full.png'), dpi=150)
        
        # Also save the simple 2-panel plot for quick reference
        plt.figure(figsize=(12, 5))
        plt.subplot(1, 2, 1); plt.imshow(force_z, extent=[x0_2d, x1_2d, y0_2d, y1_2d], origin='lower', cmap='afmhot'); plt.title("Force Z"); plt.colorbar()
        plt.subplot(1, 2, 2); plt.imshow(energy, extent=[x0_2d, x1_2d, y0_2d, y1_2d], origin='lower', cmap='viridis'); plt.title("Energy"); plt.colorbar()
        plt.savefig(os.path.join(args.outdir, 'scan_molecule_2d.png'))

    print(f"\nDone! All results saved to: {args.outdir}")

if __name__ == "__main__":
    main()
