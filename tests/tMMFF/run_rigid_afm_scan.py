#!/usr/bin/env python3

import os
import sys
import argparse
import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(BASE)

from pyBall.OCL.RigidBodyAFM import RigidBodyAFM
import matplotlib.pyplot as plt

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
    ap.add_argument('--anchor', type=int, default=28, help='Anchor atom index')
    ap.add_argument('--kanchor', type=float, default=20.0, help='Anchor stiffness')
    ap.add_argument('--nscan', type=int, default=250)
    ap.add_argument('--x0', type=float, default=0.0)
    ap.add_argument('--x1', type=float, default=50.0)
    ap.add_argument('--y0', type=float, default=2.0)
    ap.add_argument('--z0', type=float, default=8.0)
    ap.add_argument('--dt', type=float, default=0.05)
    ap.add_argument('--nsteps', type=int, default=2000)
    ap.add_argument('--fconv', type=float, default=1e-3)
    ap.add_argument('--tconv', type=float, default=1e-3)
    args = ap.parse_args()

    data_root = os.path.join(BASE, 'cpp', 'common_resources')
    mol_path = os.path.join(data_root, args.mol)
    
    sub_dir = os.path.join(data_root, args.sub)
    if not os.path.exists(sub_dir):
        sub_dir = data_root
        sub_xyz = os.path.join(data_root, args.sub)
    else:
        sub_xyz = os.path.join(sub_dir, args.sub + '.xyz')
        if not os.path.exists(sub_xyz):
            sub_xyz = os.path.join(data_root, 'Substrates', 'generated_rect', args.sub + '.xyz')
            
    gridff_path = os.path.join(sub_dir, args.grid)

    print(f"molecule  : {mol_path}")
    print(f"substrate : {sub_xyz}")
    print(f"grid      : {gridff_path}")

    # Load substrate
    sub_apos, sub_enames = read_xyz(sub_xyz)

    afm = RigidBodyAFM(
        mol_path, gridff_path, sub_xyz,
        type_map={'C': 'C_R', 'O': 'O_2', 'H': 'H'},
        anchor_idx=args.anchor,
        anchor_k=args.kanchor
    )
    
    # Read molecule to get anchor relative position
    mol_apos, mol_enames = read_xyz(mol_path)
    
    # Find anchor and opposite atom (for tilting diagnostic)
    ios = [i for i, e in enumerate(mol_enames) if e == 'O']
    if not ios:
        ios = list(range(len(mol_enames)))
    
    # Sort oxygens by lateral position to find "corners"
    xy_sum = mol_apos[ios, 0] + mol_apos[ios, 1]
    i_anchor = ios[np.argmin(xy_sum)]
    i_opposite = ios[np.argmax(xy_sum)]
    
    # Override anchor if specified
    if args.anchor is not None:
        i_anchor = args.anchor
        
    r_anchor = mol_apos[i_anchor]
    r_opposite = mol_apos[i_opposite]
    print(f"Anchor atom   ({i_anchor:2d}, {mol_enames[i_anchor]}): rel_pos={r_anchor}")
    print(f"Opposite atom ({i_opposite:2d}, {mol_enames[i_opposite]}): rel_pos={r_opposite}")
    
    # Calculate initial COM position to place the anchor atom exactly at the tip start
    tip_pos0 = np.array([args.x0, args.y0, args.z0])
    initial_com = tip_pos0 - r_anchor
    print(f"Initial COM pos: {initial_com}")
    
    # Place molecule at initial position
    afm.prepare(n_bodies=1, initial_positions=np.array([initial_com], dtype=np.float32))
    afm.anchor_idx = i_anchor # ensure correct anchor in AFM object if it was changed
    afm.set_anchor_positions(tip_pos0) # Set initial anchor

    
    xs = np.linspace(args.x0, args.x1, args.nscan)
    ys = np.full_like(xs, args.y0)
    zs = np.full_like(xs, args.z0)
    
    traj_path = "afm_scan.xyz"
    print(f"Running scan, saving to {traj_path}")
    
    com_traj = []
    anchor_traj = []
    opposite_traj = []
    tip_traj = []
    forces_norm = []
    energies = []
    
    with open(traj_path, 'w') as f:
        for i in range(args.nscan):
            tip_pos = [xs[i], ys[i], zs[i]]
            afm.set_anchor_positions(tip_pos)
            
            # Relax
            outputs, converged = afm.relax_to_constraint(
                nsteps=args.nsteps, dt=args.dt, fconv=args.fconv, tconv=args.tconv
            )
            
            # Save frame
            mol_apos = outputs['atom_positions'][0]
            all_apos = np.vstack([sub_apos, mol_apos[:,:3]])
            all_enames = list(sub_enames) + list(afm.rbd.enames)
            
            com_pos = outputs['pos'][0,:3]
            anchor_pos = mol_apos[i_anchor, :3]
            opp_pos = mol_apos[i_opposite, :3]
            f_body = outputs['body_force'][0,:3]
            f_norm = np.linalg.norm(f_body)
            # Energy is stored in the 4th component of atom_force or apos_world (fe.w)
            energy = np.sum(outputs['atom_force'][0, :, 3])
            
            com_traj.append(com_pos.copy())
            anchor_traj.append(anchor_pos.copy())
            opposite_traj.append(opp_pos.copy())
            tip_traj.append(np.array(tip_pos))
            forces_norm.append(f_norm)
            energies.append(energy)
            
            write_xyz_frame(f, all_enames, all_apos, comment=f"frame {i} tip {tip_pos} F_norm {f_norm:.4e} E {energy:.4f} converged {converged[0]}")
            print(f"Step {i:3d}/{args.nscan}: tip=({tip_pos[0]:.2f}, {tip_pos[1]:.2f}, {tip_pos[2]:.2f}) F={f_norm:.4e} E={energy:10.4f} conv={converged[0]} anc_z={anchor_pos[2]:.2f} opp_z={opp_pos[2]:.2f}")

    # --- Plotting Diagnostics ---
    com_traj = np.array(com_traj)
    anchor_traj = np.array(anchor_traj)
    opposite_traj = np.array(opposite_traj)
    tip_traj = np.array(tip_traj)
    
    fig, axs = plt.subplots(4, 1, figsize=(8, 14), sharex=True)
    
    # Plot 1: Z trajectories
    axs[0].plot(xs, com_traj[:, 2], label='COM Z')
    axs[0].plot(xs, anchor_traj[:, 2], label='Anchor Atom Z')
    axs[0].plot(xs, opposite_traj[:, 2], label='Opposite Atom Z')
    axs[0].plot(xs, tip_traj[:, 2], 'k--', label='Tip Z')
    axs[0].set_ylabel('Z Position (A)')
    axs[0].set_title('Trajectories and Tilt during X-scan')
    axs[0].legend()
    axs[0].grid(True)
    
    # Plot 2: Lateral Displacement
    axs[1].plot(xs, anchor_traj[:, 0] - tip_traj[:, 0], label='Anchor dx (from tip)')
    axs[1].plot(xs, anchor_traj[:, 1] - tip_traj[:, 1], label='Anchor dy (from tip)')
    axs[1].set_ylabel('Lateral Displacement (A)')
    axs[1].legend()
    axs[1].grid(True)
    
    # Plot 3: Force norm
    axs[2].semilogy(xs, forces_norm, label='Body Force Norm')
    axs[2].set_ylabel('Force (eV/A)')
    axs[2].legend()
    axs[2].grid(True)

    # Plot 4: Energy
    axs[3].plot(xs, energies, label='Total GridFF Energy', color='red')
    axs[3].set_xlabel('Tip X Position (A)')
    axs[3].set_ylabel('Energy (eV)')
    axs[3].legend()
    axs[3].grid(True)
    
    plt.tight_layout()
    plt.savefig('afm_scan_diagnostics.png')
    print("Saved diagnostics to afm_scan_diagnostics.png")

    # --- Plotting Projections (XY, XZ) ---
    fig2, axs2 = plt.subplots(2, 1, figsize=(10, 10))
    
    # Plot XY projection (Top View)
    axs2[0].plot(com_traj[:, 0], com_traj[:, 1], label='COM')
    axs2[0].plot(anchor_traj[:, 0], anchor_traj[:, 1], label='Anchor Atom')
    axs2[0].plot(opposite_traj[:, 0], opposite_traj[:, 1], label='Opposite Atom')
    axs2[0].plot(tip_traj[:, 0], tip_traj[:, 1], 'k--', alpha=0.5, label='Tip Path')
    axs2[0].set_xlabel('X Position (A)')
    axs2[0].set_ylabel('Y Position (A)')
    axs2[0].set_title('Top View (XY Projection)')
    axs2[0].legend()
    axs2[0].grid(True)
    axs2[0].set_aspect('equal')
    
    # Plot XZ projection (Side View)
    axs2[1].plot(com_traj[:, 0], com_traj[:, 2], label='COM')
    axs2[1].plot(anchor_traj[:, 0], anchor_traj[:, 2], label='Anchor Atom')
    axs2[1].plot(opposite_traj[:, 0], opposite_traj[:, 2], label='Opposite Atom')
    axs2[1].plot(tip_traj[:, 0], tip_traj[:, 2], 'k--', alpha=0.5, label='Tip Path')
    axs2[1].set_xlabel('X Position (A)')
    axs2[1].set_ylabel('Z Position (A)')
    axs2[1].set_title('Side View (XZ Projection)')
    axs2[1].legend()
    axs2[1].grid(True)
    # Don't force equal aspect for XZ if Z range is small
    
    plt.tight_layout()
    plt.savefig('afm_scan_projections.png')
    print("Saved projections to afm_scan_projections.png")

if __name__ == '__main__':
    main()
