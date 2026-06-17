#!/usr/bin/env python3
"""
test_viz_rigid_XPBD.py - Lightweight visualization test for 3D rigid XPBD

Demonstrates real-time interactive visualization using FFutils.MolViz3D
with raycasting-based mouse picking.

Usage:
    python test_viz_rigid_XPBD.py
    python test_viz_rigid_XPBD.py --xyz /path/to/molecule.xyz
"""

import os
import sys
import argparse
import numpy as np

_THIS_DIR = os.path.abspath(os.path.dirname(__file__))
_SHARED_DIR = os.path.abspath(os.path.join(_THIS_DIR, '..', 'shared'))
_REPO_ROOT = os.path.abspath(os.path.join(_THIS_DIR, os.pardir))
for _p in (_THIS_DIR, _SHARED_DIR, _REPO_ROOT):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from pyBall.FFutils import MolViz3D, attach_picker_3d, check_nans, print_run_header
from pyBall.AtomicSystem import AtomicSystem

# Import XPDB simulation
from XPDB_new import XPDB_new, build_neighs_bk_from_bonds


def setup_simple_molecule():
    """Create a simple water-like molecule for testing."""
    elems = ['O', 'H', 'H']
    pos = np.array([
        [0.0, 0.0, 0.0],
        [0.96, 0.0, 0.0],
        [-0.24, 0.93, 0.0],
    ], dtype=np.float32)
    bonds = [(0, 1), (0, 2)]
    return elems, pos, bonds, 1  # nnode=1


def setup_from_xyz(xyz_path, bAllNodes=False):
    """Load molecule from XYZ file."""
    mol = AtomicSystem(fname=xyz_path)
    if mol.bonds is None or len(mol.bonds) == 0:
        mol.findBonds()
    elems = list(mol.enames)
    pos = np.asarray(mol.apos[:, :3], dtype=np.float32)
    bonds = [(int(b[0]), int(b[1])) for b in mol.bonds]
    
    if bAllNodes:
        nnode = len(elems)
    else:
        nnode = sum(1 for e in elems if e != 'H')
    if nnode == 0:
        nnode = 1
    
    return elems, pos, bonds, nnode


def main():
    ap = argparse.ArgumentParser(description="Rigid XPBD 3D visualization test using FFutils")
    ap.add_argument('--xyz', type=str, default=None, help='Path to XYZ file')
    ap.add_argument('--iters', type=int, default=500, help='Number of simulation steps')
    ap.add_argument('--dt', type=float, default=0.1, help='Time step')
    ap.add_argument('--all-nodes', action='store_true', help='Treat all atoms as nodes')
    ap.add_argument('--verbose', type=int, default=0, help='Verbosity level')
    args = ap.parse_args()
    
    # Setup molecule
    if args.xyz:
        elems, pos, bonds, nnode = setup_from_xyz(args.xyz, bAllNodes=args.all_nodes)
    else:
        elems, pos, bonds, nnode = setup_simple_molecule()
    
    natoms = len(elems)
    
    print_run_header(natoms=natoms, nnode=nnode, nbonds=len(bonds), title="Rigid XPBD 3D Visualization Test")
    
    # Build topology
    neighs, bks = build_neighs_bk_from_bonds(natoms, bonds)
    
    # Create visualization using FFutils
    print("\n[INFO] Creating 3D visualization with FFutils.MolViz3D")
    viz = MolViz3D(elems=elems, show_labels=True)
    
    # Attach 3D mouse picker with raycasting
    print("[INFO] Attaching 3D mouse picker (raycasting-based)")
    
    picked_pos = {'original': None}
    
    def on_pick(pick):
        print(f"  Picked atom {pick['idx']}: {elems[pick['idx']]} at {pick['mouse3']}")
        picked_pos['original'] = pos[pick['idx']].copy()
    
    def on_drag(pick):
        if pick['idx'] is not None and pick['active']:
            # Update position during drag (for visualization)
            pos[pick['idx']] = pick['mouse3']
    
    def on_release(pick):
        if picked_pos['original'] is not None:
            # Restore original position on release (or keep new position)
            pass
        picked_pos['original'] = None
    
    pick = attach_picker_3d(viz, pick_radius_px=25, verbose=args.verbose,
                           on_pick=on_pick, on_drag=on_drag, on_release=on_release)
    
    # Simulation loop with visualization
    print(f"\n[INFO] Running {args.iters} visualization frames")
    print("[INFO] Click and drag atoms to interact (3D raycasting)")
    print("[INFO] Close window to exit\n")
    
    # Add small perturbation for visual interest
    vel = np.random.randn(natoms, 3).astype(np.float32) * 0.01
    
    for step in range(args.iters):
        # Simple damped dynamics for demo
        pos += vel * args.dt
        vel *= 0.99  # Damping
        
        # Check for NaNs
        nan_counts = check_nans(pos, vel=vel)
        if nan_counts['pos'] > 0:
            print(f"[WARN] Step {step}: NaN detected!")
            break
        
        # Update visualization
        viz.update(pos, bonds=bonds, title=f"Rigid XPBD 3D - Step {step}")
        
        # Check if window closed
        if not viz.plt.fignum_exists(viz.fig.number):
            print("\n[INFO] Window closed by user")
            break
    
    print("\n[INFO] Visualization complete")
    viz.plt.show(block=True)


if __name__ == "__main__":
    main()
