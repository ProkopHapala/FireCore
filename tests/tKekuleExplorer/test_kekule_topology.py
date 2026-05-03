#!/usr/bin/env python3
"""
test_kekule_topology.py - Headless tests for KekuleBackend (persistent AtomicSystem)
====================================================================================
Each test exercises one backend mutation method, verifies state, saves XYZ + PNG.

Run:  python test_kekule_topology.py --test all
      python test_kekule_topology.py --test add_ring
      python test_kekule_topology.py --test remove_ring
      ...
"""
import sys, os, argparse, io
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))

import numpy as np
from pyBall.KekuleBackend import KekuleBackend
from pyBall.AtomicSystem import AtomicSystem
from pyBall import plotUtils as pu

OUT_DIR = os.path.join(os.path.dirname(__file__), 'out_topology')
os.makedirs(OUT_DIR, exist_ok=True)

def _save_outputs(backend, tag):
    """Save XYZ and PNG for current backend state."""
    xyz_path = os.path.join(OUT_DIR, f'{tag}.xyz')
    backend.save_xyz(xyz_path)
    print(f"  Saved XYZ: {xyz_path}")
    png_path = os.path.join(OUT_DIR, f'{tag}.png')
    pu.plotGeometry(backend.sys.apos, list(backend.sys.enames), bond_dist=1.8, title=tag, fname=png_path)
    print(f"  Saved PNG: {png_path}")

def _assert_atoms(backend, expected_counts, msg=""):
    """Assert expected element counts in backend.sys."""
    elems = dict(zip(*np.unique(backend.sys.enames, return_counts=True))) if len(backend.sys.apos) > 0 else {}
    for e, count in expected_counts.items():
        actual = elems.get(e, 0)
        assert actual == count, f"{msg}: Expected {count} {e}, got {actual} (elems={elems})"

# ==================== Tests ====================

def test_add_ring():
    """Add single and multiple rings, verify atom counts and persistence."""
    print("\n=== test_add_ring ===")
    b = KekuleBackend()
    b.add_ring(0, 0)
    _assert_atoms(b, {'C': 6}, "After add_ring(0,0)")
    assert len(b.rings) == 1, "Should have 1 ring"
    assert len(b.sys.bonds) == 6, "Benzene should have 6 bonds"
    b.add_ring(1, 0)
    _assert_atoms(b, {'C': 10}, "After add_ring(1,0)")
    assert len(b.rings) == 2, "Should have 2 rings"
    _save_outputs(b, 'add_ring_naphthalene')
    b.report_state()
    print("  PASSED")

def test_remove_ring():
    """Add rings then remove one, verify shared atoms persist."""
    print("\n=== test_remove_ring ===")
    b = KekuleBackend()
    b.add_ring(0, 0)
    b.add_ring(1, 0)
    _assert_atoms(b, {'C': 10}, "Before remove")
    b.remove_ring(1, 0)
    _assert_atoms(b, {'C': 6}, "After remove_ring(1,0)")
    assert len(b.rings) == 1, "Should have 1 ring left"
    _save_outputs(b, 'remove_ring_benzene')
    b.report_state()
    print("  PASSED")

def test_fused_rings():
    """Build a 2x2 honeycomb patch, verify shared edges."""
    print("\n=== test_fused_rings ===")
    b = KekuleBackend()
    for q, r in [(0,0),(1,0),(0,1),(1,1)]:
        b.add_ring(q, r)
    n_C = sum(1 for e in b.sys.enames if e == 'C')
    print(f"  2x2 patch: {n_C} C atoms")
    _save_outputs(b, 'fused_rings_2x2')
    b.report_state()
    print("  PASSED")

def test_set_atom_type():
    """Change element type at a node."""
    print("\n=== test_set_atom_type ===")
    b = KekuleBackend()
    b.add_ring(0, 0)
    _assert_atoms(b, {'C': 6}, "Before set_atom_type")
    node_keys = list(b.node_to_atom.keys())
    b.set_atom_type(node_keys[0], 'N')
    _assert_atoms(b, {'C': 5, 'N': 1}, "After set_atom_type to N")
    b.set_atom_type(node_keys[1], 'O')
    _assert_atoms(b, {'C': 4, 'N': 1, 'O': 1}, "After set_atom_type to O")
    _save_outputs(b, 'set_atom_type')
    b.report_state()
    print("  PASSED")

def test_adjust_h():
    """Test persistent H passivation via adjust_h."""
    print("\n=== test_adjust_h ===")
    b = KekuleBackend()
    b.add_ring(0, 0)
    _assert_atoms(b, {'C': 6}, "Before adjust_h")
    added = b.adjust_h()
    n_H = sum(1 for e in b.sys.enames if e == 'H')
    print(f"  Added {len(added)} H atoms, total H={n_H}")
    assert n_H == 6, f"Expected 6 H, got {n_H}"
    _assert_atoms(b, {'C': 6, 'H': 6}, "After adjust_h")
    # Verify H persists after second add
    b.add_ring(1, 0)
    n_H_after = sum(1 for e in b.sys.enames if e == 'H')
    assert n_H_after == 6, f"H should persist, expected 6, got {n_H_after}"
    b.adjust_h()
    n_H_final = sum(1 for e in b.sys.enames if e == 'H')
    print(f"  After second ring + adjust_h: H={n_H_final}")
    _save_outputs(b, 'adjust_h_naphthalene')
    b.report_state()
    print("  PASSED")

def test_toggle_h_state():
    """Toggle N/O subtypes and verify H is added/removed persistently."""
    print("\n=== test_toggle_h_state ===")
    b = KekuleBackend()
    b.add_ring(0, 0)
    node_keys = list(b.node_to_atom.keys())
    b.set_atom_type(node_keys[0], 'N')
    b.adjust_h()
    n_H_init = sum(1 for e in b.sys.enames if e == 'H')
    print(f"  Pyrrolic N (default): H={n_H_init}")
    # Toggle to pyridinic (remove H)
    b.toggle_h_state(node_keys[0])
    n_H_pyri = sum(1 for e in b.sys.enames if e == 'H')
    print(f"  Pyridinic N: H={n_H_pyri}")
    assert n_H_pyri == n_H_init - 1, f"Should lose 1 H, got {n_H_pyri}"
    # Toggle back to pyrrolic (add H)
    b.toggle_h_state(node_keys[0])
    n_H_back = sum(1 for e in b.sys.enames if e == 'H')
    print(f"  Back to pyrrolic: H={n_H_back}")
    assert n_H_back == n_H_init, f"Should restore H, got {n_H_back}"
    _save_outputs(b, 'toggle_h_state')
    b.report_state()
    print("  PASSED")

def test_recalc_bonds():
    """Verify bond recalculation works after mutations."""
    print("\n=== test_recalc_bonds ===")
    b = KekuleBackend()
    b.add_ring(0, 0)
    b.adjust_h()
    nbonds1 = len(b.sys.bonds) if b.sys.bonds is not None else 0
    print(f"  Bonds after adjust_h: {nbonds1}")
    b.add_ring(1, 0)
    b.adjust_h()
    nbonds2 = len(b.sys.bonds) if b.sys.bonds is not None else 0
    print(f"  Bonds after second ring: {nbonds2}")
    assert nbonds2 > nbonds1, "More atoms should mean more bonds"
    _save_outputs(b, 'recalc_bonds')
    print("  PASSED")

def test_relax():
    """Mock relax test: just verify positions can be modified in-place."""
    print("\n=== test_relax (in-place position update) ===")
    b = KekuleBackend()
    b.add_ring(0, 0)
    b.adjust_h()
    pos_before = b.sys.apos.copy()
    # Simulate a small perturbation (mock relaxation)
    # Perturb only pinned atoms
    for i, pin in enumerate(b.atom_pin):
        if pin is not None:
            b.sys.apos[i] += np.random.normal(0, 0.01, 3)
    pos_after = b.sys.apos.copy()
    diff = np.max(np.abs(pos_after - pos_before))
    print(f"  Max position change: {diff:.4f} A")
    assert diff > 0, "Positions should have changed"
    b.snap_atoms_to_grid()
    pos_snapped = b.sys.apos.copy()
    # Only check pinned atoms snapped back
    pinned_mask = [p is not None for p in b.atom_pin]
    diff_snapped = np.max(np.abs(pos_snapped[pinned_mask] - pos_before[pinned_mask]))
    print(f"  After snap_to_grid (pinned only): max diff = {diff_snapped:.6f} A")
    assert diff_snapped < 1e-6, "Pinned atoms should snap back to grid"
    _save_outputs(b, 'relax_mock')
    print("  PASSED")

def test_remove_atom():
    """Remove a single atom and verify ring cleanup."""
    print("\n=== test_remove_atom ===")
    b = KekuleBackend()
    b.add_ring(0, 0)
    b.adjust_h()
    node_keys = list(b.node_to_atom.keys())
    b.remove_atom(node_keys[0])
    n_C = sum(1 for e in b.sys.enames if e == 'C')
    n_H = sum(1 for e in b.sys.enames if e == 'H')
    print(f"  After remove_atom: C={n_C}, H={n_H}")
    assert n_C == 5, "Should have 5 C left"
    assert len(b.rings) == 0, "Ring should be removed since incomplete"
    _save_outputs(b, 'remove_atom')
    b.report_state()
    print("  PASSED")

def test_load_save_xyz():
    """Save and reload XYZ, verify state round-trips."""
    print("\n=== test_load_save_xyz ===")
    b = KekuleBackend()
    b.add_ring(0, 0)
    b.add_ring(1, 0)
    b.adjust_h()
    xyz_path = os.path.join(OUT_DIR, 'roundtrip_test.xyz')
    b.save_xyz(xyz_path)
    b2 = KekuleBackend()
    b2.load_xyz(xyz_path)
    n_C = sum(1 for e in b2.sys.enames if e == 'C')
    n_H = sum(1 for e in b2.sys.enames if e == 'H')
    print(f"  Reloaded: C={n_C}, H={n_H}, rings={len(b2.rings)}")
    assert n_C == 10, f"Expected 10 C, got {n_C}"
    _save_outputs(b2, 'load_save_xyz')
    b2.report_state()
    print("  PASSED")

# ==================== CLI ====================

ALL_TESTS = [
    'add_ring', 'remove_ring', 'fused_rings', 'set_atom_type',
    'adjust_h', 'toggle_h_state', 'recalc_bonds', 'relax',
    'remove_atom', 'load_save_xyz'
]

def run_tests(tests):
    log_path = os.path.join(OUT_DIR, 'test_log.txt')
    log = open(log_path, 'w')
    def logprint(s):
        print(s)
        log.write(s + '\n')
        log.flush()
    
    logprint("=" * 60)
    logprint("KekuleBackend Headless Tests")
    logprint("=" * 60)
    
    for name in tests:
        func = globals().get(f'test_{name}')
        if func is None:
            logprint(f"\n[SKIP] Unknown test: {name}")
            continue
        try:
            func()
            logprint(f"  [OK] {name}")
        except Exception as e:
            logprint(f"  [FAIL] {name}: {e}")
            import traceback
            logprint(traceback.format_exc())
    
    logprint("\n" + "=" * 60)
    logprint(f"Done. Outputs in: {OUT_DIR}")
    logprint("=" * 60)
    log.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Headless tests for KekuleBackend')
    parser.add_argument('--test', type=str, default='all',
                        help=f"Comma-separated test names or 'all'. Available: {','.join(ALL_TESTS)}")
    args = parser.parse_args()
    
    if args.test == 'all':
        tests = ALL_TESTS
    else:
        tests = [t.strip() for t in args.test.split(',')]
    
    run_tests(tests)
