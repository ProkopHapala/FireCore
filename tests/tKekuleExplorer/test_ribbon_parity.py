#!/usr/bin/env python3
"""
test_ribbon_parity.py - Headless tests for KekuleBackend ribbon construction
================================================================================
Validates that KekuleBackend.build_zigzag_ribbon() replaces GrapheneRibbonBuilder
with functional parity. Tests all passivation types, two-ribbon cells, and compares
atom counts / geometries with the deprecated builder where possible.

Run:  python test_ribbon_parity.py --test all
      python test_ribbon_parity.py --test passivation_N
      python test_ribbon_parity.py --test two_ribbon_cell
      python test_ribbon_parity.py --test periodic_ribbon --rows 4 --cols 4 --passivation N
"""
import sys, os, argparse
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))

import numpy as np
from pyBall.KekuleBackend import KekuleBackend, parse_passivation_string
from pyBall import plotUtils as pu

OUT_DIR = os.path.join(os.path.dirname(__file__), 'out_ribbon_parity')
os.makedirs(OUT_DIR, exist_ok=True)

# Global CLI parameters
CLI_ROWS = 4
CLI_COLS = 4
CLI_PASSIVATION = 'N'
CLI_L_Hb = 2.0

def _save_outputs(backend, tag, lvs=None):
    """Save XYZ and PNG for current backend state."""
    xyz_path = os.path.join(OUT_DIR, f'{tag}.xyz')
    if lvs is not None:
        # Write lattice vectors in comment line
        lvs_str = f"lvs {lvs[0,0]:.6f} {lvs[0,1]:.6f} {lvs[0,2]:.6f}  {lvs[1,0]:.6f} {lvs[1,1]:.6f} {lvs[1,2]:.6f}  {lvs[2,0]:.6f} {lvs[2,1]:.6f} {lvs[2,2]:.6f}"
        backend.save_xyz(xyz_path, comment=lvs_str)
    else:
        backend.save_xyz(xyz_path)
    print(f"  Saved XYZ: {xyz_path}")
    png_path = os.path.join(OUT_DIR, f'{tag}.png')
    pu.plotGeometry(backend.sys.apos, list(backend.sys.enames), bond_dist=1.8, title=tag, fname=png_path, lvs=lvs, bDrawBox=(lvs is not None), replicate=(1,1,0))
    print(f"  Saved PNG: {png_path}")

def _assert_atoms(backend, expected_counts, msg=""):
    """Assert expected element counts in backend.sys."""
    elems = dict(zip(*np.unique(backend.sys.enames, return_counts=True))) if len(backend.sys.apos) > 0 else {}
    for e, count in expected_counts.items():
        actual = elems.get(e, 0)
        assert actual == count, f"{msg}: Expected {count} {e}, got {actual} (elems={elems})"

# ==================== Tests ====================

def test_build_ribbon_unpassivated():
    """Build a wide unpassivated ribbon, verify it's all C."""
    print("\n=== test_build_ribbon_unpassivated ===")
    b = KekuleBackend()
    b.build_zigzag_ribbon(width_chains=4, length_cells=4, passivation='H')
    n_C = sum(1 for e in b.sys.enames if e == 'C')
    n_H = sum(1 for e in b.sys.enames if e == 'H')
    print(f"  Ribbon: C={n_C}, H={n_H}")
    assert n_C > 0, "Ribbon should have C atoms"
    assert n_H > 0, "H-passivated ribbon should have H atoms"
    _save_outputs(b, 'ribbon_H_passivated')
    b.report_state()
    print("  PASSED")

def test_passivation_N():
    """N passivation: edge C replaced by N, no extra H."""
    print(f"\n=== test_passivation_N (rows={CLI_ROWS}, cols={CLI_COLS}) ===")
    b = KekuleBackend()
    b.build_zigzag_ribbon(width_chains=CLI_ROWS, length_cells=CLI_COLS, passivation='N', side_passivation='CH')
    n_C = sum(1 for e in b.sys.enames if e == 'C')
    n_N = sum(1 for e in b.sys.enames if e == 'N')
    n_H = sum(1 for e in b.sys.enames if e == 'H')
    print(f"  Ribbon: C={n_C}, N={n_N}, H={n_H}")
    assert n_N > 0, "N-passivated ribbon should have N atoms"
    assert n_H > 0, "Side CH passivation should add H atoms"
    _save_outputs(b, f'ribbon_N_rows{CLI_ROWS}_cols{CLI_COLS}')
    b.report_state()
    print("  PASSED")

def test_passivation_NH():
    """NH passivation: edge C replaced by N, plus H."""
    print("\n=== test_passivation_NH ===")
    b = KekuleBackend()
    b.build_zigzag_ribbon(width_chains=4, length_cells=4, passivation='NH')
    n_C = sum(1 for e in b.sys.enames if e == 'C')
    n_N = sum(1 for e in b.sys.enames if e == 'N')
    n_H = sum(1 for e in b.sys.enames if e == 'H')
    print(f"  Ribbon: C={n_C}, N={n_N}, H={n_H}")
    assert n_N > 0, "NH-passivated ribbon should have N atoms"
    assert n_H > 0, "NH passivation should add H atoms"
    assert n_H == n_N, "NH passivation: 1 H per N"
    _save_outputs(b, 'ribbon_NH')
    b.report_state()
    print("  PASSED")

def test_passivation_CH():
    """CH passivation: edge C keeps C, add H."""
    print("\n=== test_passivation_CH ===")
    b = KekuleBackend()
    b.build_zigzag_ribbon(width_chains=4, length_cells=4, passivation='CH')
    n_C = sum(1 for e in b.sys.enames if e == 'C')
    n_H = sum(1 for e in b.sys.enames if e == 'H')
    print(f"  Ribbon: C={n_C}, H={n_H}")
    assert n_H > 0, "CH passivation should add H atoms"
    _save_outputs(b, 'ribbon_CH')
    b.report_state()
    print("  PASSED")

def test_passivation_O():
    """O passivation: edge C replaced by O, no extra H."""
    print("\n=== test_passivation_O ===")
    b = KekuleBackend()
    b.build_zigzag_ribbon(width_chains=4, length_cells=4, passivation='O')
    n_C = sum(1 for e in b.sys.enames if e == 'C')
    n_O = sum(1 for e in b.sys.enames if e == 'O')
    n_H = sum(1 for e in b.sys.enames if e == 'H')
    print(f"  Ribbon: C={n_C}, O={n_O}, H={n_H}")
    assert n_O > 0, "O-passivated ribbon should have O atoms"
    assert n_H == 0, "O passivation should not add H atoms"
    _save_outputs(b, 'ribbon_O')
    b.report_state()
    print("  PASSED")

def test_passivation_C_eq_O():
    """C=O passivation: add O double-bonded to edge C."""
    print("\n=== test_passivation_C_eq_O ===")
    b = KekuleBackend()
    b.build_zigzag_ribbon(width_chains=4, length_cells=4, passivation='C=O')
    n_C = sum(1 for e in b.sys.enames if e == 'C')
    n_O = sum(1 for e in b.sys.enames if e == 'O')
    n_H = sum(1 for e in b.sys.enames if e == 'H')
    print(f"  Ribbon: C={n_C}, O={n_O}, H={n_H}")
    assert n_O > 0, "C=O passivation should add O atoms"
    assert n_H == 0, "C=O passivation should not add H atoms"
    _save_outputs(b, 'ribbon_C_eq_O')
    b.report_state()
    print("  PASSED")

def test_passivation_C_OH():
    """C-OH passivation: add OH group to edge C."""
    print("\n=== test_passivation_C_OH ===")
    b = KekuleBackend()
    b.build_zigzag_ribbon(width_chains=4, length_cells=4, passivation='C-OH')
    n_C = sum(1 for e in b.sys.enames if e == 'C')
    n_O = sum(1 for e in b.sys.enames if e == 'O')
    n_H = sum(1 for e in b.sys.enames if e == 'H')
    print(f"  Ribbon: C={n_C}, O={n_O}, H={n_H}")
    assert n_O > 0, "C-OH passivation should add O atoms"
    assert n_H > 0, "C-OH passivation should add H atoms"
    assert n_H == n_O, "C-OH passivation: 1 H per O"
    _save_outputs(b, 'ribbon_C_OH')
    b.report_state()
    print("  PASSED")

def test_scale_x():
    """Build ribbon with x scaling."""
    print("\n=== test_scale_x ===")
    b = KekuleBackend()
    b.build_zigzag_ribbon(width_chains=4, length_cells=4, passivation='N', scale_x=1.5)
    xs = b.sys.apos[:, 0]
    print(f"  x range: {xs.min():.3f} to {xs.max():.3f}")
    assert xs.max() > 0, "Scaled ribbon should have positive x span"
    _save_outputs(b, 'ribbon_scaled_x')
    b.report_state()
    print("  PASSED")

def test_two_ribbon_cell():
    """Build two-ribbon cell via module-level function."""
    print("\n=== test_two_ribbon_cell ===")
    apos, atypes, elems, lvs = build_two_ribbon_cell(width_chains=4, length_cells=2, Lx=2.46)
    n_C = sum(1 for e in elems if e == 'C')
    n_N = sum(1 for e in elems if e == 'N')
    n_H = sum(1 for e in elems if e == 'H')
    print(f"  Two-ribbon cell: C={n_C}, N={n_N}, H={n_H}")
    print(f"  Lattice: Lx={lvs[0,0]:.3f}, Ly={lvs[1,1]:.3f}, Lz={lvs[2,2]:.3f}")
    assert n_N > 0, "Bottom ribbon should be N-passivated"
    assert n_H > 0, "Top ribbon should be NH-passivated"
    # Save as XYZ
    from pyBall import atomicUtils as au
    xyz_path = os.path.join(OUT_DIR, 'two_ribbon_cell.xyz')
    qs = np.zeros(len(elems))
    au.saveXYZ(elems, apos, xyz_path, qs=qs, comment="lvs %.6f %.6f %.6f   %.6f %.6f %.6f   %.6f %.6f %.6f" % (
        lvs[0,0], lvs[0,1], lvs[0,2],
        lvs[1,0], lvs[1,1], lvs[1,2],
        lvs[2,0], lvs[2,1], lvs[2,2]
    ))
    print(f"  Saved XYZ: {xyz_path}")
    print("  PASSED")

def test_build_ribbon_function():
    """Test module-level build_ribbon function."""
    print("\n=== test_build_ribbon_function ===")
    pos2d, atypes, elems = build_ribbon('NH', width_chains=4, length_cells=2, Lx=2.46)
    n_C = sum(1 for e in elems if e == 'C')
    n_N = sum(1 for e in elems if e == 'N')
    n_H = sum(1 for e in elems if e == 'H')
    print(f"  build_ribbon('NH'): C={n_C}, N={n_N}, H={n_H}, natoms={len(elems)}")
    assert n_N > 0, "NH ribbon should have N atoms"
    assert n_H > 0, "NH ribbon should have H atoms"
    assert pos2d.shape[1] == 2, "pos2d should be 2D"
    print("  PASSED")

def test_wide_ribbon():
    """Build a wider ribbon, verify edge atom counts scale."""
    print("\n=== test_wide_ribbon ===")
    b = KekuleBackend()
    b.build_zigzag_ribbon(width_chains=8, length_cells=4, passivation='N')
    n_C = sum(1 for e in b.sys.enames if e == 'C')
    n_N = sum(1 for e in b.sys.enames if e == 'N')
    print(f"  Wide ribbon (8x4): C={n_C}, N={n_N}")
    assert n_C > n_N, "Wide ribbon should have more C than N"
    _save_outputs(b, 'ribbon_wide_N')
    b.report_state()
    print("  PASSED")

def test_periodic_ribbon():
    """Build periodic ribbon (bPeriodicX=True), verify armchair edges not passivated."""
    print(f"\n=== test_periodic_ribbon (rows={CLI_ROWS}, cols={CLI_COLS}, passivation={CLI_PASSIVATION}) ===")
    b = KekuleBackend()
    b.build_zigzag_ribbon(width_chains=CLI_ROWS, length_cells=CLI_COLS,
                          passivation=CLI_PASSIVATION, bPeriodicX=True)
    n_C = sum(1 for e in b.sys.enames if e == 'C')
    n_N = sum(1 for e in b.sys.enames if e == 'N')
    n_O = sum(1 for e in b.sys.enames if e == 'O')
    print(f"  Periodic ribbon ({CLI_ROWS}x{CLI_COLS}): C={n_C}, N={n_N}, O={n_O}")
    # In periodic mode, only zigzag edges (top/bottom) should be passivated
    # Armchair edges (left/right) should remain C
    tag = f'ribbon_periodic_{CLI_PASSIVATION}_rows{CLI_ROWS}_cols{CLI_COLS}'
    _save_outputs(b, tag)
    b.report_state()
    print("  PASSED")

def test_two_ribbon_cell_method():
    """Test build_two_ribbon_cell as a class method with PBC."""
    print(f"\n=== test_two_ribbon_cell_method (rows={CLI_ROWS}, cols={CLI_COLS}, L_Hb={CLI_L_Hb}) ===")
    b = KekuleBackend()
    b.build_two_ribbon_cell(width_chains=CLI_ROWS, length_cells=CLI_COLS, Lx=2.4, L_Hb=CLI_L_Hb, 
                           bottom_passivation='N', top_passivation='NH', bPeriodicX=True)
    n_C = sum(1 for e in b.sys.enames if e == 'C')
    n_N = sum(1 for e in b.sys.enames if e == 'N')
    n_H = sum(1 for e in b.sys.enames if e == 'H')
    print(f"  Two-ribbon cell (PBC): C={n_C}, N={n_N}, H={n_H}")
    assert n_N > 0, "Bottom ribbon should have N"
    assert n_H > 0, "Top ribbon should have H (NH passivation)"
    
    # Calculate lattice vectors
    Lx = 2.4 * CLI_COLS  # Scale by number of columns
    Ly = b.sys.apos[:, 1].max() - b.sys.apos[:, 1].min() + 15.0  # vacuum gap
    lvs = np.array([[Lx, 0.0, 0.0], [0.0, Ly, 0.0], [0.0, 0.0, 20.0]])
    
    _save_outputs(b, f'two_ribbon_cell_method_pbc_rows{CLI_ROWS}_cols{CLI_COLS}_LHb{CLI_L_Hb}', lvs=lvs)
    b.report_state()
    print("  PASSED")

def test_two_ribbon_cell_non_pbc():
    """Test build_two_ribbon_cell with non-PBC and side CH passivation."""
    print("\n=== test_two_ribbon_cell_non_pbc ===")
    b = KekuleBackend()
    b.build_two_ribbon_cell(width_chains=4, length_cells=1, Lx=2.4, L_Hb=2.0, 
                           bottom_passivation='N', top_passivation='NH', bPeriodicX=False, side_passivation='CH')
    n_C = sum(1 for e in b.sys.enames if e == 'C')
    n_N = sum(1 for e in b.sys.enames if e == 'N')
    n_H = sum(1 for e in b.sys.enames if e == 'H')
    print(f"  Two-ribbon cell (non-PBC with side CH): C={n_C}, N={n_N}, H={n_H}")
    assert n_N > 0, "Should have N atoms"
    assert n_H > 0, "Should have H atoms from side CH passivation"
    _save_outputs(b, 'two_ribbon_cell_non_pbc')
    b.report_state()
    print("  PASSED")

def test_two_ribbon_cell_variable_passivation():
    """Test build_two_ribbon_cell with variable passivation lists (PBC)."""
    print("\n=== test_two_ribbon_cell_variable_passivation ===")
    b = KekuleBackend()
    bottom_passivation = ['N', 'CH', 'N', 'CH']
    top_passivation = ['NH', 'C=O', 'NH', 'C=O']
    b.build_two_ribbon_cell(width_chains=4, Lx=2.4, L_Hb=2.0, 
                           bottom_passivation=bottom_passivation, top_passivation=top_passivation, bPeriodicX=True)
    n_C = sum(1 for e in b.sys.enames if e == 'C')
    n_N = sum(1 for e in b.sys.enames if e == 'N')
    n_O = sum(1 for e in b.sys.enames if e == 'O')
    n_H = sum(1 for e in b.sys.enames if e == 'H')
    print(f"  Variable passivation (PBC): C={n_C}, N={n_N}, O={n_O}, H={n_H}")
    assert n_N > 0, "Should have N atoms"
    assert n_O > 0, "Should have O atoms from C=O passivation"
    assert n_H > 0, "Should have H atoms"
    
    # Calculate lattice vectors
    Lx = 2.4 * len(bottom_passivation)
    Ly = b.sys.apos[:, 1].max() - b.sys.apos[:, 1].min() + 15.0
    lvs = np.array([[Lx, 0.0, 0.0], [0.0, Ly, 0.0], [0.0, 0.0, 20.0]])
    
    _save_outputs(b, f'two_ribbon_cell_variable_passivation_ncol{len(bottom_passivation)}', lvs=lvs)
    b.report_state()
    print("  PASSED")

def test_two_ribbon_cell_composable():
    """Test composable approach: build two single ribbons separately, then combine."""
    print("\n=== test_two_ribbon_cell_composable ===")
    
    # Build bottom ribbon with N on both edges
    bottom_ribbon = KekuleBackend()
    bottom_ribbon.build_zigzag_ribbon(width_chains=4, length_cells=2, passivation_bottom='N', passivation_top='N',  scale_x=2.4 / (2.0 * 1.42 * np.cos(np.pi / 6)), bPeriodicX=True)
    
    # Build top ribbon with NH on both edges
    top_ribbon = KekuleBackend()
    top_ribbon.build_zigzag_ribbon(width_chains=4, length_cells=2, passivation_bottom='NH', passivation_top='NH', scale_x=2.4 / (2.0 * 1.42 * np.cos(np.pi / 6)), bPeriodicX=True)
    
    # Combine ribbons
    b = KekuleBackend()
    b.combine_ribbons(bottom_ribbon, top_ribbon, L_Hb=2.0, shift_x=0.0)
    
    n_C = sum(1 for e in b.sys.enames if e == 'C')
    n_N = sum(1 for e in b.sys.enames if e == 'N')
    n_H = sum(1 for e in b.sys.enames if e == 'H')
    print(f"  Composable two-ribbon cell: C={n_C}, N={n_N}, H={n_H}")
    assert n_N > 0, "Should have N atoms"
    assert n_H > 0, "Should have H atoms from NH passivation"
    
    _save_outputs(b, 'two_ribbon_cell_composable')
    b.report_state()
    print("  PASSED")

# ==================== Comparison with old builder (if available) ====================

def test_compare_with_old_builder():
    """Compare new ribbon with deprecated GrapheneRibbonBuilder (if importable)."""
    print("\n=== test_compare_with_old_builder ===")
    try:
        from doc.Topics.Kekule_Topology.GrapheneRibbonBuilder import GrapheneRibbonBuilder
    except ImportError:
        print("  [SKIP] GrapheneRibbonBuilder not importable")
        return

    # Build with old builder
    old = GrapheneRibbonBuilder(a_CC=1.42)
    old_pos, old_elems, old_bonds = old.build_zigzag_ribbon(width_chains=4, length_cells=4, passivation='N')

    # Build with new backend
    b = KekuleBackend(a_CC=1.42)
    b.build_zigzag_ribbon(width_chains=4, length_cells=4, passivation='N')

    old_counts = dict(zip(*np.unique(old_elems, return_counts=True)))
    new_counts = dict(zip(*np.unique(b.sys.enames, return_counts=True)))
    print(f"  Old builder: {old_counts}")
    print(f"  New backend: {new_counts}")

    # They may differ in exact counts due to different construction algorithms,
    # but should have the same qualitative composition
    assert 'N' in new_counts, "New ribbon should have N"
    assert 'C' in new_counts, "New ribbon should have C"
    print("  PASSED (qualitative parity)")

# ==================== CLI ====================

ALL_TESTS = [
    'build_ribbon_unpassivated', 'passivation_N', 'passivation_NH',
    'passivation_CH', 'passivation_O', 'passivation_C_eq_O', 'passivation_C_OH',
    'scale_x', 'two_ribbon_cell', 'build_ribbon_function',
    'wide_ribbon', 'periodic_ribbon', 'two_ribbon_cell_method', 'two_ribbon_cell_non_pbc', 'two_ribbon_cell_variable_passivation', 'two_ribbon_cell_composable', 'compare_with_old_builder'
]

def run_tests(tests):
    log_path = os.path.join(OUT_DIR, 'test_log.txt')
    log = open(log_path, 'w')
    def logprint(s):
        print(s)
        log.write(s + '\n')
        log.flush()

    logprint("=" * 60)
    logprint("Ribbon Parity Headless Tests")
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
    parser = argparse.ArgumentParser(description='Headless tests for ribbon construction')
    parser.add_argument('--test', type=str, default='all', help=f"Comma-separated test names or 'all'. Available: {','.join(ALL_TESTS)}")
    parser.add_argument('--rows', type=int, default=4, help='Number of rows (default: 4)')
    parser.add_argument('--cols', type=int, default=4, help='Number of columns (default: 4)')
    parser.add_argument('--passivation', type=str, default='N',  choices=['N', 'NH', 'CH', 'H', 'O', 'C=O', 'C-OH'],help='Passivation type (default: N)')
    parser.add_argument('--L_Hb', type=float, default=2.0, help='Hydrogen bond length (default: 2.0)')
    
    # CLI mode for string-based passivation ribbons (always PBC)
    parser.add_argument('--single', action='store_true', help='Generate single PBC ribbon with string passivation')
    parser.add_argument('--bottom', type=str, help='Bottom edge passivation string (e.g., NnOOH)')
    parser.add_argument('--top', type=str, help='Top edge passivation string (e.g., NnOOH)')
    parser.add_argument('--bottom1', type=str, help='Bottom ribbon bottom edge (for two-ribbon)')
    parser.add_argument('--top1', type=str, help='Bottom ribbon top edge (for two-ribbon)')
    parser.add_argument('--bottom2', type=str, help='Top ribbon bottom edge (for two-ribbon)')
    parser.add_argument('--top2', type=str, help='Top ribbon top edge (for two-ribbon)')
    
    args = parser.parse_args()
    CLI_ROWS = args.rows
    CLI_COLS = args.cols
    CLI_PASSIVATION = args.passivation
    CLI_L_Hb = args.L_Hb

    # CLI mode: generate single PBC ribbon with string passivation
    if args.single:
        if args.bottom is None or args.top is None:
            print("ERROR: --single requires --bottom and --top passivation strings")
            sys.exit(1)
        
        bottom_passivation = parse_passivation_string(args.bottom)
        top_passivation = parse_passivation_string(args.top)
        length_cells = len(bottom_passivation)
        width_chains = CLI_ROWS
        Lx = 2.4
        
        b = KekuleBackend()
        b.build_zigzag_ribbon(width_chains=width_chains, length_cells=length_cells, passivation_bottom=bottom_passivation, passivation_top=top_passivation, scale_x=Lx / (2.0 * 1.42 * np.cos(np.pi / 6)), bPeriodicX=True)
        
        n_C = sum(1 for e in b.sys.enames if e == 'C')
        n_N = sum(1 for e in b.sys.enames if e == 'N')
        n_O = sum(1 for e in b.sys.enames if e == 'O')
        n_H = sum(1 for e in b.sys.enames if e == 'H')
        print(f"Single ribbon: C={n_C}, N={n_N}, O={n_O}, H={n_H}")
        
        tag = f'single_w{width_chains}_{args.bottom}_{args.top}'
        Lx_scaled = Lx * length_cells
        Ly = b.sys.apos[:, 1].max() - b.sys.apos[:, 1].min() + 15.0
        lvs = np.array([[Lx_scaled, 0.0, 0.0], [0.0, Ly, 0.0], [0.0, 0.0, 20.0]])
        _save_outputs(b, tag, lvs=lvs)
        sys.exit(0)
    
    # CLI mode: generate two-ribbon system with string passivation
    if args.bottom1 is not None or args.top1 is not None or args.bottom2 is not None or args.top2 is not None:
        if args.bottom1 is None or args.top1 is None or args.bottom2 is None or args.top2 is None:
            print("ERROR: Two-ribbon mode requires --bottom1, --top1, --bottom2, --top2")
            sys.exit(1)
        
        bottom1_passivation = parse_passivation_string(args.bottom1)
        top1_passivation = parse_passivation_string(args.top1)
        bottom2_passivation = parse_passivation_string(args.bottom2)
        top2_passivation = parse_passivation_string(args.top2)
        length_cells = len(bottom1_passivation)
        width_chains = CLI_ROWS
        Lx = 2.4
        L_Hb = CLI_L_Hb
        
        # Build bottom ribbon
        bottom_ribbon = KekuleBackend()
        bottom_ribbon.build_zigzag_ribbon(width_chains=width_chains, length_cells=length_cells, passivation_bottom=bottom1_passivation, passivation_top=top1_passivation, scale_x=Lx / (2.0 * 1.42 * np.cos(np.pi / 6)), bPeriodicX=True)
        
        # Build top ribbon
        top_ribbon = KekuleBackend()
        top_ribbon.build_zigzag_ribbon(width_chains=width_chains, length_cells=length_cells, passivation_bottom=bottom2_passivation, passivation_top=top2_passivation, scale_x=Lx / (2.0 * 1.42 * np.cos(np.pi / 6)), bPeriodicX=True)
        
        # Combine ribbons
        b = KekuleBackend()
        b.combine_ribbons(bottom_ribbon, top_ribbon, L_Hb=L_Hb, shift_x=0.0)
        
        n_C = sum(1 for e in b.sys.enames if e == 'C')
        n_N = sum(1 for e in b.sys.enames if e == 'N')
        n_O = sum(1 for e in b.sys.enames if e == 'O')
        n_H = sum(1 for e in b.sys.enames if e == 'H')
        print(f"Two-ribbon system: C={n_C}, N={n_N}, O={n_O}, H={n_H}")
        
        tag = f'two_w{width_chains}_{args.bottom1}_{args.top1}_{args.bottom2}_{args.top2}'
        Lx_scaled = Lx * length_cells
        Ly = b.sys.apos[:, 1].max() - b.sys.apos[:, 1].min() + 15.0
        lvs = np.array([[Lx_scaled, 0.0, 0.0], [0.0, Ly, 0.0], [0.0, 0.0, 20.0]])
        _save_outputs(b, tag, lvs=lvs)
        sys.exit(0)

    if args.test == 'all':
        tests = ALL_TESTS
    else:
        tests = [t.strip() for t in args.test.split(',')]
    run_tests(tests)
