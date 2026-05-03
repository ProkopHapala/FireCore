#!/usr/bin/env python3
"""
test_kekule_topology.py - Phase 1 headless test for Kekule Explorer
====================================================================
Tests: hexagonal grid building, sp2 passivation, undercoordination detection.
Outputs: XYZ files and PNG plots for user review.

Run:  python test_kekule_topology.py
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))

import numpy as np
from pyBall.KekuleBackend import KekuleBackend
from pyBall.AtomicSystem import AtomicSystem
from pyBall import plotUtils as pu

OUT_DIR = os.path.join(os.path.dirname(__file__), 'out_topology')
os.makedirs(OUT_DIR, exist_ok=True)

def test_benzene():
    """Test 1: Single benzene ring - 6 C atoms, should get 6 H after passivation."""
    print("\n=== Test 1: Benzene ===")
    backend = KekuleBackend()
    backend.toggle_ring(0, 0)
    
    # Build without passivation first
    sys_bare = backend.build_system(bPassivate=False)
    print(f"  Bare: {len(sys_bare.apos)} atoms, {len(sys_bare.bonds)} bonds")
    assert len(sys_bare.apos) == 6, f"Expected 6 atoms, got {len(sys_bare.apos)}"
    
    # Check undercoordinated
    undercoord = sys_bare.find_undercoordinated()
    print(f"  Undercoordinated atoms: {len(undercoord)}")
    for ia, nm in undercoord:
        print(f"    atom {ia} ({sys_bare.enames[ia]}): missing {nm} bonds")
    
    # Build with passivation
    sys_h = backend.build_system(bPassivate=True)
    n_H = sum(1 for e in sys_h.enames if e == 'H')
    print(f"  Passivated: {len(sys_h.apos)} atoms ({n_H} H), {len(sys_h.bonds)} bonds")
    assert n_H == 6, f"Expected 6 H atoms, got {n_H}"
    
    # Save outputs
    xyz_path = os.path.join(OUT_DIR, 'benzene.xyz')
    sys_h.saveXYZ(xyz_path, blvec=False)
    print(f"  Saved XYZ: {xyz_path}")
    
    png_path = os.path.join(OUT_DIR, 'benzene.png')
    pu.plotGeometry(sys_h.apos, list(sys_h.enames), bond_dist=1.8, title='Benzene (C6H6)', fname=png_path)
    print(f"  Saved PNG: {png_path}")
    print("  PASSED")

def test_naphthalene():
    """Test 2: Two fused rings (naphthalene) - 10 C atoms, 8 H after passivation."""
    print("\n=== Test 2: Naphthalene ===")
    backend = KekuleBackend()
    backend.toggle_ring(0, 0)
    backend.toggle_ring(1, 0)
    
    sys_bare = backend.build_system(bPassivate=False)
    n_C = sum(1 for e in sys_bare.enames if e == 'C')
    print(f"  Bare: {n_C} C atoms, {len(sys_bare.bonds)} bonds")
    assert n_C == 10, f"Expected 10 C atoms, got {n_C}"
    
    sys_h = backend.build_system(bPassivate=True)
    n_H = sum(1 for e in sys_h.enames if e == 'H')
    n_C = sum(1 for e in sys_h.enames if e == 'C')
    print(f"  Passivated: {n_C} C + {n_H} H = {len(sys_h.apos)} atoms, {len(sys_h.bonds)} bonds")
    assert n_H == 8, f"Expected 8 H atoms, got {n_H}"
    
    xyz_path = os.path.join(OUT_DIR, 'naphthalene.xyz')
    sys_h.saveXYZ(xyz_path, blvec=False)
    print(f"  Saved XYZ: {xyz_path}")
    
    png_path = os.path.join(OUT_DIR, 'naphthalene.png')
    pu.plotGeometry(sys_h.apos, list(sys_h.enames), bond_dist=1.8, title='Naphthalene (C10H8)', fname=png_path)
    print(f"  Saved PNG: {png_path}")
    print("  PASSED")

def test_pyridine():
    """Test 3: Benzene with one C replaced by N (pyridine) - 5C + 1N, should get 5H."""
    print("\n=== Test 3: Pyridine ===")
    backend = KekuleBackend()
    backend.toggle_ring(0, 0)
    
    # Replace one node with N
    node_keys = sorted(backend.nodes.keys())
    backend.set_atom(node_keys[0], 'N')
    
    sys_bare = backend.build_system(bPassivate=False)
    n_C = sum(1 for e in sys_bare.enames if e == 'C')
    n_N = sum(1 for e in sys_bare.enames if e == 'N')
    print(f"  Bare: {n_C} C + {n_N} N = {len(sys_bare.apos)} atoms")
    
    # N in pyridine has 2 heavy neighbors (sp2, =N-), so target_valence N=2 (pyridinic)
    # By default target_valence N=3, so it would add 1H making it pyrrolic (-NH-).
    # For pyridinic N we need target_valence N=2 
    # Let's test both: default (pyrrolic) and pyridinic
    
    # Pyrrolic N (default: N needs 3 sigma bonds = 2 heavy + 1 H)
    sys_pyrrolic = backend.build_system(bPassivate=True)
    n_H_pyrrolic = sum(1 for e in sys_pyrrolic.enames if e == 'H')
    print(f"  Pyrrolic (N=3): {n_H_pyrrolic} H atoms")
    
    xyz_path = os.path.join(OUT_DIR, 'pyrrole_like.xyz')
    sys_pyrrolic.saveXYZ(xyz_path, blvec=False)
    png_path = os.path.join(OUT_DIR, 'pyrrole_like.png')
    pu.plotGeometry(sys_pyrrolic.apos, list(sys_pyrrolic.enames), bond_dist=1.8, title='Pyrrole-like (N gets H)', fname=png_path)
    print(f"  Saved PNG: {png_path}")
    
    # Pyridinic N (N=2, no H on nitrogen)
    sys_bare2 = backend.build_system(bPassivate=False)
    added = sys_bare2.add_capping_h_sp2(target_valence={'C': 3, 'N': 2, 'O': 2})
    n_H_pyridinic = sum(1 for e in sys_bare2.enames if e == 'H')
    print(f"  Pyridinic (N=2): {n_H_pyridinic} H atoms")
    
    xyz_path2 = os.path.join(OUT_DIR, 'pyridine.xyz')
    sys_bare2.saveXYZ(xyz_path2, blvec=False)
    png_path2 = os.path.join(OUT_DIR, 'pyridine.png')
    pu.plotGeometry(sys_bare2.apos, list(sys_bare2.enames), bond_dist=1.8, title='Pyridine (N=2, no H on N)', fname=png_path2)
    print(f"  Saved PNG: {png_path2}")
    print("  PASSED")

def test_anthracene_ribbon():
    """Test 4: Linear 3-ring ribbon (anthracene) - 14 C, 10 H."""
    print("\n=== Test 4: Anthracene ribbon ===")
    backend = KekuleBackend()
    backend.toggle_ring(0, 0)
    backend.toggle_ring(1, 0)
    backend.toggle_ring(2, 0)
    
    sys_bare = backend.build_system(bPassivate=False)
    n_C = sum(1 for e in sys_bare.enames if e == 'C')
    print(f"  Bare: {n_C} C atoms")
    assert n_C == 14, f"Expected 14 C atoms, got {n_C}"
    
    sys_h = backend.build_system(bPassivate=True)
    n_H = sum(1 for e in sys_h.enames if e == 'H')
    print(f"  Passivated: {n_C} C + {n_H} H = {len(sys_h.apos)} atoms")
    assert n_H == 10, f"Expected 10 H atoms, got {n_H}"
    
    xyz_path = os.path.join(OUT_DIR, 'anthracene.xyz')
    sys_h.saveXYZ(xyz_path, blvec=False)
    png_path = os.path.join(OUT_DIR, 'anthracene.png')
    pu.plotGeometry(sys_h.apos, list(sys_h.enames), bond_dist=1.8, title='Anthracene (C14H10)', fname=png_path)
    print(f"  Saved PNG: {png_path}")
    print("  PASSED")

def test_hbond_detection():
    """Test 5: H-bond detection using existing findHBonds on a pyrrolic system."""
    print("\n=== Test 5: H-bond detection ===")
    # Build two small pyrrole-like fragments stacked vertically 
    backend = KekuleBackend()
    backend.toggle_ring(0, 0)
    node_keys = sorted(backend.nodes.keys())
    # Make top node N (pyrrolic, will get H pointing up)
    top_key = max(node_keys, key=lambda k: k[1])
    backend.set_atom(top_key, 'N')
    sys1 = backend.build_system(bPassivate=True)
    
    # Use the existing findHBonds from AtomicSystem
    hbonds, hb_rs = sys1.findHBonds(Rb=1.5, Rh=3.0, angMax=60.0, bPrint=True)
    print(f"  H-bonds found: {len(hbonds)}")
    if len(hbonds) > 0:
        for hb in hbonds:
            print(f"    H-bond: atoms {hb}")
    
    xyz_path = os.path.join(OUT_DIR, 'hbond_test.xyz')
    sys1.saveXYZ(xyz_path, blvec=False)
    print(f"  Saved XYZ: {xyz_path}")
    print("  PASSED (H-bond detection runs, results depend on geometry)")

if __name__ == '__main__':
    print("=" * 60)
    print("Kekule Explorer - Phase 1: Topology Tests (Headless)")
    print("=" * 60)
    
    test_benzene()
    test_naphthalene()
    test_pyridine()
    test_anthracene_ribbon()
    test_hbond_detection()
    
    print("\n" + "=" * 60)
    print(f"All tests passed. Output files in: {OUT_DIR}")
    print("=" * 60)
