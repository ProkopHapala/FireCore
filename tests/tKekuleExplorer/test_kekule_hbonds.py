#!/usr/bin/env python3
"""
test_kekule_hbonds.py - Phase 2 headless test for Kekule Explorer
==================================================================
Tests: AtomQuery selection, H-bond detection between two fragments.
Outputs: Console report and PNG plot.

Run:  python test_kekule_hbonds.py
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))

import numpy as np
from pyBall.KekuleBackend import KekuleBackend
from pyBall.AtomicSystem import AtomicSystem
from pyBall import plotUtils as pu

OUT_DIR = os.path.join(os.path.dirname(__file__), 'out_hbonds')
os.makedirs(OUT_DIR, exist_ok=True)

def test_atom_query():
    """Test AtomQuery selection capabilities."""
    print("\n=== Test: AtomQuery Selection ===")
    backend = KekuleBackend()
    backend.toggle_ring(0, 0)
    # Pyridine node
    node_keys = sorted(backend.nodes.keys())
    backend.set_atom(node_keys[0], 'N')
    
    sys = backend.build_system(bPassivate=True)
    
    # Select Carbons with 2 neighbors (edge carbons in benzene)
    # Wait, in benzene all C have 2 heavy neighbors + 1 H. 
    # Total degree is 3.
    c_deg3 = sys.select_atoms("C deg{3}")
    print(f"  C with degree 3: {len(c_deg3)} (Expected 5 in pyridine)")
    
    # Select Nitrogen
    n_atoms = sys.select_atoms("N")
    print(f"  Nitrogen atoms: {len(n_atoms)} (Expected 1)")
    
    # Select Nitrogens with 1 H (pyrrolic by default)
    n_pyrrolic = sys.select_atoms("N n{H}{1}")
    print(f"  Pyrrolic N (1 H): {len(n_pyrrolic)} (Expected 1)")
    
    # Select H attached to N
    h_on_n = sys.select_atoms("H n{N}{1}")
    print(f"  H on N: {len(h_on_n)} (Expected 1)")
    
    # Toggle H on N (Remove H -> Pyridinic)
    sys.toggle_h_state(n_atoms[0])
    h_on_n_after = sys.select_atoms("H n{N}{1}")
    print(f"  H on N after toggle: {len(h_on_n_after)} (Expected 0)")
    
    # Select Nitrogens with 0 H (pyridinic)
    n_pyridinic_after = sys.select_atoms("N n{H}{0}")
    print(f"  Pyridinic N after toggle: {len(n_pyridinic_after)} (Expected 1)")
    
    assert len(n_pyrrolic) == 1
    assert len(h_on_n) == 1
    assert len(h_on_n_after) == 0
    assert len(n_pyridinic_after) == 1
    print("  PASSED")

def test_hbond_ribbons():
    """Test H-bond detection between two parallel fragments."""
    print("\n=== Test: H-bond Detection between ribbons ===")
    
    # Create two pyrrole-like rings separated by 3.5 A
    backend = KekuleBackend()
    backend.toggle_ring(0, 0)
    node_keys = sorted(backend.nodes.keys())
    # Top node is N (pyrrolic)
    top_key = max(node_keys, key=lambda k: k[1])
    backend.set_atom(top_key, 'N')
    sys1 = backend.build_system(bPassivate=True)
    
    # Create second ring shifted up
    backend2 = KekuleBackend()
    backend2.toggle_ring(0, 2) # shift in r-direction
    node_keys2 = sorted(backend2.nodes.keys())
    # Bottom node is O (ketone-like)
    bot_key = min(node_keys2, key=lambda k: k[1])
    backend2.set_atom(bot_key, 'O')
    sys2 = backend2.build_system(bPassivate=True)
    # Manually shift sys2 to get good H-bond distance
    # Ring 1 N is at y=1.42. H is at y=2.43
    # Ring 2 O is at y=2.84.
    # Shift sys2 up by 1.7 A -> O at y=4.54. H...O distance = 4.54 - 2.43 = 2.11 A
    sys2.apos[:,1] += 1.7
    sys2.apos[:,0] -= backend2.a_CC * np.sqrt(3) * 1.0 # Align x (q+r/2=0 for r=2 means q=-1)
    
    # Merge systems
    sys1.addSystems(sys2)
    sys1.neighs()
    
    print(f"  Total atoms: {len(sys1.apos)}")
    
    # Find H-bonds
    hbonds = sys1.find_hbonds(d_max=2.8, a_min=140.0)
    print(f"  H-bonds found: {len(hbonds)}")
    for d, h, a, dist, ang in hbonds:
        print(f"    {sys1.enames[d]}({d})-H({h}) ... {sys1.enames[a]}({a})  dist={dist:.3f} ang={ang:.1f}")
    
    # Plotting
    png_path = os.path.join(OUT_DIR, 'hbond_ribbons.png')
    # We'll use plotGeometry but we'd like to highlight H-bonds. 
    # For now just plot the geometry.
    pu.plotGeometry(sys1.apos, list(sys1.enames), bond_dist=1.8, title='H-bond Test', fname=png_path)
    print(f"  Saved PNG: {png_path}")
    
    assert len(hbonds) > 0, "No H-bonds found!"
    print("  PASSED")

if __name__ == '__main__':
    print("=" * 60)
    print("Kekule Explorer - Phase 2: Pattern Matching & H-Bonds")
    print("=" * 60)
    
    test_atom_query()
    test_hbond_ribbons()
    
    print("\n" + "=" * 60)
    print(f"All tests passed. Output files in: {OUT_DIR}")
    print("=" * 60)
