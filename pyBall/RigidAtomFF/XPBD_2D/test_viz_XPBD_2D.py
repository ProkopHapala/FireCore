#!/usr/bin/env python3
"""
test_viz_XPBD_2D.py - Lightweight visualization test for 2D XPBD

Demonstrates real-time interactive visualization using FFutils.MolViz2D
with BLIT-based matplotlib animation and mouse picking.

Usage:
    python test_viz_XPBD_2D.py
    python test_viz_XPBD_2D.py --xyz /path/to/molecule.xyz
"""

import os
import sys
import argparse
import numpy as np

_THIS_DIR = os.path.abspath(os.path.dirname(__file__))
_REPO_ROOT = os.path.abspath(os.path.join(_THIS_DIR, '..'))
for _p in (_THIS_DIR, _REPO_ROOT):
    if _p not in sys.path:
        sys.path.insert(0, _p)

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from pyBall.FFutils import MolViz2D, attach_picker_2d, check_nans, print_run_header
from pyBall.AtomicSystem import AtomicSystem

# Import XPBD_2D simulation
from XPBD_2D import XPBD_2D, build_neighs_bk_from_bonds_2d, make_bk_slots_2d, make_stiffness_from_bonds_2d


def setup_simple_chain(n_atoms=5, bond_length=1.5):
    """Create a simple chain molecule for testing."""
    elems = ['C'] + ['H'] * (n_atoms - 1)
    pos = np.zeros((n_atoms, 2), dtype=np.float32)
    for i in range(n_atoms):
        pos[i, 0] = i * bond_length
    bonds = [(i, i+1) for i in range(n_atoms - 1)]
    return elems, pos, bonds


def setup_from_xyz(xyz_path):
    """Load molecule from XYZ file."""
    mol = AtomicSystem(fname=xyz_path)
    if mol.bonds is None or len(mol.bonds) == 0:
        mol.findBonds()
    elems = list(mol.enames)
    pos = np.asarray(mol.apos[:, :2], dtype=np.float32)  # Use only X,Y
    bonds = [(int(b[0]), int(b[1])) for b in mol.bonds]
    return elems, pos, bonds


def main():
    ap = argparse.ArgumentParser(description="XPBD_2D visualization test using FFutils")
    ap.add_argument('--xyz', type=str, default=None, help='Path to XYZ file')
    ap.add_argument('--n-atoms', type=int, default=5, help='Number of atoms for chain (if no xyz)')
    ap.add_argument('--iters', type=int, default=500, help='Number of simulation steps')
    ap.add_argument('--dt', type=float, default=0.02, help='Time step')
    ap.add_argument('--inner-iters', type=int, default=5, help='Inner PBD iterations')
    ap.add_argument('--view-scale', type=float, default=15.0, help='Viewport scale')
    ap.add_argument('--verbose', type=int, default=0, help='Verbosity level')
    args = ap.parse_args()
    
    # Setup molecule
    if args.xyz:
        elems, pos, bonds = setup_from_xyz(args.xyz)
    else:
        elems, pos, bonds = setup_simple_chain(args.n_atoms)
    
    natoms = len(elems)
    nnode = sum(1 for e in elems if e != 'H')  # Non-H atoms are nodes
    if nnode == 0:
        nnode = 1
    
    print_run_header(natoms=natoms, nnode=nnode, nbonds=len(bonds), title="XPBD_2D Visualization Test")
    
    # Build topology
    neighs, bks = build_neighs_bk_from_bonds_2d(natoms, bonds)
    bkSlots = make_bk_slots_2d(neighs, nnode=nnode, natoms=natoms)
    stiffness = make_stiffness_from_bonds_2d(natoms, neighs, k_bond=200.0)
    
    # Initialize simulator
    sim = XPBD_2D(natoms)
    sim.upload_topology(neighs, bkSlots, stiffness)
    
    # Initialize state
    rot = np.zeros((nnode, 2), dtype=np.float32)
    rot[:, 0] = 1.0  # cos(0)
    vel = np.zeros((natoms, 2), dtype=np.float32)
    omega = np.zeros((nnode,), dtype=np.float32)
    sim.upload_state(pos, rot, vel, omega)
    
    # Create visualization using FFutils
    print("\n[INFO] Creating visualization with FFutils.MolViz2D")
    viz = MolViz2D(elems=elems, show_labels=True, view_scale=args.view_scale)
    
    # Attach mouse picker
    print("[INFO] Attaching mouse picker (click to drag atoms)")
    
    def on_pick(pick):
        print(f"  Picked atom {pick['idx']}: {elems[pick['idx']]}")
    
    def on_drag(pick):
        # Update atom position during drag
        if pick['idx'] is not None and pick['active']:
            new_pos = pick['mouse'].copy()
            sim.set_atom_pos(pick['idx'], new_pos)
    
    pick = attach_picker_2d(viz, pick_radius=1.0, verbose=args.verbose,
                           on_pick=on_pick, on_drag=on_drag)
    
    # Simulation loop with visualization
    print(f"\n[INFO] Running {args.iters} simulation steps with real-time visualization")
    print("[INFO] Click and drag atoms to interact")
    print("[INFO] Close window to exit\n")
    
    for step in range(args.iters):
        # Run simulation step
        sim.step_pbd(nnode=nnode, iterations=args.inner_iters, relax=0.5, bmix=0.0, reset_hb=(step == 0))
        
        # Download state
        pos_cur, rot_cur, vel_cur, omega_cur = sim.download_state()
        
        # Check for NaNs
        nan_counts = check_nans(pos_cur, vel=vel_cur)
        if nan_counts['pos'] > 0:
            print(f"[WARN] Step {step}: NaN detected in positions!")
            break
        
        # Update visualization
        info = f"Step: {step}\n|v|_max: {np.max(np.linalg.norm(vel_cur, axis=1)):.4f}"
        viz.update(pos_cur, bonds=bonds, title=f"XPBD_2D Step {step}", info=info)
        
        # Small pause for animation
        viz.plt.pause(0.01)
        
        # Check if window closed
        if not viz.plt.fignum_exists(viz.fig.number):
            print("\n[INFO] Window closed by user")
            break
    
    print("\n[INFO] Simulation complete")
    viz.plt.show(block=True)


if __name__ == "__main__":
    main()
