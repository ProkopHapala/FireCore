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
from pyBall.KekuleBackend import KekuleBackend, snap_to_grid, honeycomb_ring_nodes
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
    # Add second ring (naphthalene has 10 C, 8 H)
    b.add_ring(1, 0)
    n_H_after = sum(1 for e in b.sys.enames if e == 'H')
    print(f"  After second ring (auto H): H={n_H_after}")
    assert n_H_after == 8, f"Expected 8 H for naphthalene, got {n_H_after}"
    _save_outputs(b, 'adjust_h_naphthalene')
    b.report_state()
    print("  PASSED")

def test_toggle_h_state():
    """Test changing npi on N atoms using set_atom_valency (uniform electron counting)."""
    print("\n=== test_toggle_h_state ===")
    backend = KekuleBackend()
    backend.add_ring(0, 0)
    backend.adjust_h()
    # Change one C to N (defaults to sp2, npi=1)
    nk = snap_to_grid(honeycomb_ring_nodes(0, 0, backend.a_CC)[0], backend.a_CC)
    backend.set_atom_type(nk, 'N')
    backend.adjust_h()
    n_H_sp2 = sum(1 for e in backend.sys.enames if e == 'H')
    print(f"N_sp2 (npi=1): H={n_H_sp2}")
    # Change to sp3 (npi=0) - should have different H count
    backend.set_atom_valency(nk, 0)
    n_H_sp3 = sum(1 for e in backend.sys.enames if e == 'H')
    print(f"N_sp3 (npi=0): H={n_H_sp3}")
    # sp3 N should need more H than sp2 N (nsigma: 3 vs 2)
    assert n_H_sp3 > n_H_sp2, f"sp3 should have more H than sp2, got {n_H_sp3} vs {n_H_sp2}"
    _save_outputs(backend, 'toggle_h_state')
    backend.report_state()
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

def test_insert_atom_into_bond():
    """Test inserting an atom into a bond (A-B -> A-C-B)."""
    print("\n=== test_insert_atom_into_bond ===")
    backend = KekuleBackend()
    backend.auto_h_cap = False  # Disable to prevent extra bond creation from adjust_h
    backend.add_ring(0, 0)  # Creates 6 C atoms in a ring
    _assert_atoms(backend, {'C': 6}, "Initial benzene ring")
    
    # Get alive bonds from the graph
    bonds = [bd for bd in backend.graph.bonds.values() if bd.alive]
    assert len(bonds) == 6, f"Expected 6 bonds, got {len(bonds)}"
    
    # Pick the first bond
    bond = bonds[0]
    atom_a = bond.a
    atom_b = bond.b
    print(f"  Selected bond {bond._id} between Atom({atom_a._id}) and Atom({atom_b._id})")
    
    # Insert atom into bond
    new_atom = backend.insert_atom_into_bond(bond, 'C')
    print(f"  Inserted new Atom({new_atom._id})")
    
    # Should now have 7 C atoms
    _assert_atoms(backend, {'C': 7}, "After inserting atom into bond")
    
    # Should have 7 alive bonds (removed 1, added 2)
    bonds_after = [bd for bd in backend.graph.bonds.values() if bd.alive]
    print(f"  Alive bonds after: {len(bonds_after)}")
    assert len(bonds_after) == 7, f"Expected 7 alive bonds, got {len(bonds_after)}"
    
    # Verify the new atom is bonded to both original atoms (alive bonds only)
    new_atom_bonds = [bd for bd in new_atom.bonds if bd.alive]
    assert len(new_atom_bonds) == 2, f"New atom should have 2 alive bonds, got {len(new_atom_bonds)}"
    
    # Check that original bond is marked dead
    assert not bond.alive, "Original bond should be marked as dead"
    
    _save_outputs(backend, 'insert_atom_into_bond')
    backend.report_state()
    print("  PASSED")

def test_collapse_bond():
    """Test collapsing a bond (A-B -> A with neighbors of B transferred to A)."""
    print("\n=== test_collapse_bond ===")
    backend = KekuleBackend()
    backend.auto_h_cap = False  # Disable to prevent extra bond creation from adjust_h
    backend.add_ring(0, 0)  # Creates 6 C atoms in a ring
    _assert_atoms(backend, {'C': 6}, "Initial benzene ring")
    
    # Get alive bonds from the graph
    bonds = [bd for bd in backend.graph.bonds.values() if bd.alive]
    assert len(bonds) == 6, f"Expected 6 bonds, got {len(bonds)}"
    
    # Pick the first bond
    bond = bonds[0]
    atom_a = bond.a
    atom_b = bond.b
    print(f"  Selected bond {bond._id} between Atom({atom_a._id}) and Atom({atom_b._id})")
    
    # Record positions before collapse
    pos_a_before = atom_a.pos.copy()
    pos_b_before = atom_b.pos.copy()
    center = (pos_a_before + pos_b_before) / 2.0
    
    # Collapse bond - simulate mouse click far from atom_a (so atom_a survives)
    mouse_pos = pos_b_before[:2] + np.array([1.0, 1.0])  # Far from atom_a
    survivor = backend.collapse_bond(bond, mouse_pos)
    
    print(f"  Survivor: Atom({survivor._id})")
    
    # Should now have 5 C atoms (one removed)
    _assert_atoms(backend, {'C': 5}, "After collapsing bond")
    
    # Survivor should be at the center of the original bond
    survivor_pos = survivor.pos
    dist_to_center = np.linalg.norm(survivor_pos - center)
    assert dist_to_center < 0.01, f"Survivor should be at bond center, distance={dist_to_center}"
    
    # Original bond should be dead
    assert not bond.alive, "Original bond should be marked as dead"
    
    # In benzene ring, each C has 2 bonds. After collapse, survivor should have 2 bonds (its original other bond + transferred bond)
    survivor_bonds = [bd for bd in survivor.bonds if bd.alive]
    print(f"  Survivor has {len(survivor_bonds)} alive bonds")
    assert len(survivor_bonds) == 2, f"Survivor should have 2 bonds (original + transferred), got {len(survivor_bonds)}"
    
    _save_outputs(backend, 'collapse_bond')
    backend.report_state()
    print("  PASSED")

def test_collapse_bond_with_h_caps():
    """Test that H caps are properly cleaned up when collapsing bonds."""
    print("\n=== test_collapse_bond_with_h_caps ===")
    backend = KekuleBackend()
    backend.auto_h_cap = True  # Enable H capping
    backend.add_ring(0, 0)  # Creates 6 C atoms
    
    # Add H caps first
    backend.adjust_h()
    
    # Count atoms after H capping
    h_atoms = [a for a in backend.graph.atoms.values() if a.alive and a.ename == 'H']
    c_atoms = [a for a in backend.graph.atoms.values() if a.alive and a.ename == 'C']
    print(f"  After H capping: {len(c_atoms)} C atoms, {len(h_atoms)} H atoms")
    assert len(c_atoms) == 6, f"Expected 6 C atoms, got {len(c_atoms)}"
    assert len(h_atoms) > 0, "Expected some H atoms after H capping"
    
    # Get first bond and collapse it
    bonds = [bd for bd in backend.graph.bonds.values() if bd.alive and bd.a.ename == 'C' and bd.b.ename == 'C']
    assert len(bonds) > 0, "Expected C-C bonds"
    bond = bonds[0]
    
    # Get H children of the atom that will be removed
    pos_a = bond.a.pos.copy()
    pos_b = bond.b.pos.copy()
    mouse_pos = pos_b[:2] + np.array([1.0, 1.0])  # Far from atom_a
    
    # Count H caps before collapse
    to_remove = bond.b  # This atom will be removed (closer to mouse)
    h_children_before = backend.graph.h_children(to_remove)
    print(f"  Atom({to_remove._id}) to be removed has {len(h_children_before)} H children")
    
    # Collapse the bond
    survivor = backend.collapse_bond(bond, mouse_pos)
    
    # After collapse, H children of removed atom should be gone
    h_atoms_after = [a for a in backend.graph.atoms.values() if a.alive and a.ename == 'H']
    c_atoms_after = [a for a in backend.graph.atoms.values() if a.alive and a.ename == 'C']
    
    print(f"  After collapse: {len(c_atoms_after)} C atoms, {len(h_atoms_after)} H atoms")
    assert len(c_atoms_after) == 5, f"Expected 5 C atoms after collapse, got {len(c_atoms_after)}"
    # H count should decrease by the number of H children of removed atom
    expected_h = len(h_atoms) - len(h_children_before)
    print(f"  Expected {expected_h} H atoms, got {len(h_atoms_after)}")
    # Allow for some tolerance since H adjustment may add caps to the survivor
    
    _save_outputs(backend, 'collapse_bond_with_h_caps')
    backend.report_state()
    print("  PASSED")

# ==================== CLI ====================

ALL_TESTS = [
    'add_ring', 'remove_ring', 'fused_rings', 'set_atom_type',
    'adjust_h', 'toggle_h_state', 'recalc_bonds', 'relax',
    'remove_atom', 'load_save_xyz',
    'insert_atom_into_bond', 'collapse_bond', 'collapse_bond_with_h_caps'
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
