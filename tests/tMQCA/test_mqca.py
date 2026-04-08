"""
MQCA test script
================
Produces three sets of images:

1. test1_groundstate_<input>.png
   Geometry + ground state of a small cluster for each input combination.

2. test2_logic_map.png
   2-D phase diagram: output logic type as function of (W1, W2).

3. test3_geometry_scan.png
   Which (W1,W2) values produce useful logic for multiple geometries;
   also prints a table of found logic functions.

Run from project root:
   python -u -m pyBall.OCL.test_mqca
or from this directory:
   PYTHONPATH=../.. python -u test_mqca.py
"""

import sys, os
import argparse
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from pyBall.OCL.MQCA import (
    MQCASolver,
    sq_lattice_sparse,
    apply_input_bias,
    eval_logic_table,
    scan_W1_W2,
    scan_W1_W2_top8,
    identify_logic,
    occ_mask_to_array,
    plot_ground_states,
    plot_logic_map,
    plot_logic_fraction_map,
    INPUT_COMBOS, LOGIC_NAMES, USEFUL_LOGIC, MAX_NEIGH, LOGIC_COLORS,
)

OUT_DIR = os.path.dirname(os.path.abspath(__file__))

# Parse CLI arguments
parser = argparse.ArgumentParser(description='MQCA test script')
parser.add_argument('--format', type=str, default='png', choices=['png', 'svg'],
                   help='Output format for plots (default: png)')
args = parser.parse_args()

PLOT_EXT = args.format

# ======================================================================
#  Cluster definitions
# ======================================================================

def make_cluster_T():
    """
    T-shaped 5-site cluster.

       *  *  *     (0,1) (1,1) (2,1)
          *           (1,0)
          *           (1,-1)

    input A: pad at (-1,1), neighbors site 0 (left tip)
    input B: pad at (3, 1), neighbors site 2 (right tip)
    output : site 3 (bottom stem)
    """
    positions = np.array([
        [0,1],[1,1],[2,1],
        [1,0],[1,-1],
    ], dtype=int)
    input_positions = [(-1, 1), (3, 1)]
    input_neighbors = [[0], [2]]
    output_site     = 4
    return positions, input_positions, input_neighbors, output_site


def make_cluster_T_extended_inputs():
    """
    T-shape with longer input lines (active cells on input arms).

    Layout: A and B are now active sites, not external pads.
    Input pads are further out, with 2-site input arms.

        A_X_B    y=3: Pads at (-1,3) and (3,3), sites at (0,3), (2,3)
        _X_X_    y=2: Sites (0,2), (2,2)
        _XXX_    y=1: Crossbar (0,1),(1,1),(2,1)
        __X__    y=0: Stem (1,0)
        __O__    y=-1: Output (1,-1)
    """
    positions = np.array([
        [0,3], [2,3],   # input arm sites (now part of cluster, sites 0,1)
        [0,2], [2,2],   # middle of input arms (sites 2,3)
        [0,1], [1,1], [2,1],  # crossbar (sites 4,5,6)
        [1,0],          # stem (site 7)
        [1,-1],         # output (site 8)
    ], dtype=int)
    input_positions = [(-1, 3), (3, 3)]  # external pads further out
    input_neighbors = [[0], [1]]  # site 0 neighbors A pad, site 1 neighbors B pad
    output_site = 8
    return positions, input_positions, input_neighbors, output_site


def make_cluster_T_extended_output():
    """
    T-shape with longer output line (2-site stem).

        A_X_B    y=3: Pads at (-1,3), (3,3), sites (0,3),(2,3)
        _XXX_    y=2: Crossbar (0,2),(1,2),(2,2)
        __X__    y=1: Upper stem (1,1)
        __X__    y=0: Lower stem (1,0)
        __O__    y=-1: Output (1,-1)
    """
    positions = np.array([
        [0,3], [2,3],   # under pads (sites 0,1)
        [0,2], [1,2], [2,2],  # crossbar (sites 2,3,4)
        [1,1],          # upper stem (site 5)
        [1,0],          # lower stem (site 6)
        [1,-1],         # output (site 7)
    ], dtype=int)
    input_positions = [(-1, 3), (3, 3)]
    input_neighbors = [[0], [1]]
    output_site = 7
    return positions, input_positions, input_neighbors, output_site


def make_cluster_T_no_center():
    """
    T-shape with central site removed (no site at cross-junction).

        A_X_B    y=2: Pads at (-1,2),(3,2), sites (0,2),(2,2)
        _X_X_    y=1: Arms only (0,1),(2,1) - NO CENTER
        __X__    y=0: Stem (1,0)
        __O__    y=-1: Output (1,-1)
    """
    positions = np.array([
        [0,2], [2,2],   # under pads (sites 0,1)
        [0,1], [2,1],   # arms only, no center (sites 2,3)
        [1,0],          # stem (site 4)
        [1,-1],         # output (site 5)
    ], dtype=int)
    input_positions = [(-1, 2), (3, 2)]
    input_neighbors = [[0], [1]]
    output_site = 5
    return positions, input_positions, input_neighbors, output_site


def make_cluster_T_fork():
    """
    Fork: Extended input arms + full crossbar + longer output.
    The ASCII art from the request:
        "A_B"    y=4: Pads at (0,4), (2,4)
        "X_X"    y=3: Sites (0,3), (2,3)
        "XXX"    y=2: Crossbar (0,2),(1,2),(2,2)
        "_X_"    y=1: Upper stem (1,1)
        "_O_"    y=0: Output (1,0)

    Actually removing center means crossbar has gap at (1,2).
    Correct fork: separate left/right arms that join at the stem.
    """
    positions = np.array([
        [0,3], [2,3],   # under pads (sites 0,1)
        [0,2], [2,2],   # crossbar arms, no center (sites 2,3)
        [1,1],          # upper stem connecting them (site 4)
        [1,0],          # output (site 5)
    ], dtype=int)
    input_positions = [(0, 4), (2, 4)]
    input_neighbors = [[0], [1]]
    output_site = 5
    return positions, input_positions, input_neighbors, output_site


def make_cluster_S():
    """
    S/Z-shaped 6-site cluster – asymmetric.

       .  *  *     (1,2) (2,2)
       *  *  .     (0,1) (1,1)
       *  *  .     (0,0) (1,0)

    input A: pad at (3,2), neighbors site 1 (top-right)
    input B: pad at (-1,0), neighbors site 4 (bottom-left)
    output : site 2 (mid-left)
    """
    positions = np.array([
        [1,2],[2,2],
        [0,1],[1,1],
        [0,0],[1,0],
    ], dtype=int)
    input_positions = [(3, 2), (-1, 0)]
    input_neighbors = [[1], [4]]
    output_site     = 2
    return positions, input_positions, input_neighbors, output_site


def make_cluster_zigzag():
    """
    Zigzag 5-site cluster.

    *  .  *  .  *
    .  *  .  *  .

    sites: (0,1)(1,0)(2,1)(3,0)(4,1)
    input A: pad at (-1,1), neighbors site 0
    input B: pad at (5, 1), neighbors site 4
    output : site 2
    """
    positions = np.array([
        [0,1],[1,0],[2,1],[3,0],[4,1]
    ], dtype=int)
    input_positions = [(-1, 1), (5, 1)]
    input_neighbors = [[0], [4]]
    output_site     = 2
    return positions, input_positions, input_neighbors, output_site


def make_cluster_L():
    """
    L-shaped 6-site cluster on a square grid.

      I_A  .  .
      *  *  *
      *  *  I_B  (output = site 4, upper-right active)

    Active site layout (grid coords):
      0:(0,0)  1:(1,0)  2:(2,0)
      3:(0,1)  4:(1,1)  5:(2,1)
         ↑  ↑
    input A: pad at (-1,0), neighbors site 0
    input B: pad at (2,2),  neighbors site 5
    output : site 4
    """
    positions = np.array([
        [0,0],[1,0],[2,0],
        [0,1],[1,1],[2,1],
    ], dtype=int)
    input_positions = [(-1, 0), (3, 1)]   # for plotting only
    input_neighbors = [[0], [5]]           # which active sites each input biases
    output_site     = 4
    return positions, input_positions, input_neighbors, output_site


def make_cluster_cross():
    """
    Cross / plus shaped 5-site cluster.

         *          (2,2)
      *  *  *    (1,1) (2,1) (3,1)
         *          (2,0)

    input A: pad at (2,3), neighbors site 4 (top)
    input B: pad at (0,1), neighbors site 1 (left)
    output : site 3 (right)
    """
    positions = np.array([
        [2,0],          # 0 bottom
        [1,1],[2,1],[3,1],  # 1 left, 2 centre, 3 right
        [2,2],          # 4 top
    ], dtype=int)
    input_positions = [(2, 3), (0, 1)]
    input_neighbors = [[4], [1]]
    output_site     = 3
    return positions, input_positions, input_neighbors, output_site


def make_cluster_chain():
    """
    Simple 4-site horizontal chain.

    I_A  *  *  *  *  I_B

    sites: 0 1 2 3
    input A biases site 0, input B biases site 3
    output: site 1 or 2
    """
    positions = np.array([[i,0] for i in range(4)], dtype=int)
    input_positions = [(-1, 0), (4, 0)]
    input_neighbors = [[0], [3]]
    output_site     = 2
    return positions, input_positions, input_neighbors, output_site

def make_cluster_straight_line_center():
    """
    5-site straight line with inputs on ends, output in center.

    I_A  *  *  O  *  *  I_B

    sites: 0 1 2 3 4
    input A biases site 0, input B biases site 4
    output: site 2 (center)
    """
    positions = np.array([[i,0] for i in range(5)], dtype=int)
    input_positions = [(-1, 0), (5, 0)]
    input_neighbors = [[0], [4]]
    output_site     = 2
    return positions, input_positions, input_neighbors, output_site


# ======================================================================
#  Shared solver
# ======================================================================
print("Initialising MQCASolver …")
solver = MQCASolver(preferred_vendor='nvidia', bPrint=True)

# ======================================================================
#  TEST 0 – sanity check: 2-site system, analytic energy known
# ======================================================================
print("\n=== TEST 0: 2-site sanity check ===")
# 2 sites: eps=[0.5, -0.5], W12=1.0
# States:
#  00: E=0
#  01: E=eps1=-0.5    (site 1 occupied)
#  10: E=eps0=0.5     (site 0 occupied)
#  11: E=eps0+eps1+W=0.5-0.5+1.0=1.0
# Ground state: E=-0.5, occ=01 (mask=0b10=2)
Esite2 = np.array([[0.5, -0.5]], dtype=np.float32)
W_val2 = np.zeros((2, MAX_NEIGH), dtype=np.float32)
W_idx2 = np.zeros((2, MAX_NEIGH), dtype=np.int32)
nNeigh2 = np.array([1, 1], dtype=np.int32)
W_val2[0,0] = 1.0;  W_idx2[0,0] = 1
W_val2[1,0] = 1.0;  W_idx2[1,0] = 0
E2, occ2 = solver.solve(Esite2, W_val2, W_idx2, nNeigh2, nSite=2)
print(f"  E_min={E2[0]:.4f} (expect -0.5)   occ_mask={occ2[0]:02b} (expect 10 = bit1 set = site1 occupied)")
assert abs(E2[0] - (-0.5)) < 1e-4, f"Sanity check failed: E={E2[0]}"
print("  PASSED")

# ======================================================================
#  TEST 1 – single cluster, plot all 4 ground states
# ======================================================================
print("\n=== TEST 1: single cluster ground states (Cross cluster at NAND) ===")

# Use Cross cluster which shows NAND with eps0=-1 physics
positions, input_positions, input_neighbors, output_site = make_cluster_cross()
nSite = len(positions)
Esite_base = -np.ones(nSite, dtype=np.float32)  # eps0 = -1 for all (prefer occupancy)

# First do a coarse scan to find NAND region for Cross cluster
print("  Scanning for NAND logic...")
nW_scan = 40
W1_scan = np.linspace(-3.0, 3.0, nW_scan)
W2_scan = np.linspace(-3.0, 3.0, nW_scan)
lmap, occ_all, E_all = scan_W1_W2(solver, positions, input_positions, Esite_base, input_neighbors, output_site,
                                   W1_scan, W2_scan, shift=-0.5)

# Find NAND (code 7) locations
nand_locs = np.argwhere(lmap == 7)
if len(nand_locs) > 0:
    # Pick middle of NAND region
    iy, ix = nand_locs[len(nand_locs)//2]
    W1, W2 = float(W1_scan[ix]), float(W2_scan[iy])
    print(f"  Found NAND at W1={W1:.3f}, W2={W2:.3f}")
else:
    # Fallback
    W1, W2 = 1.5, -0.5
    print(f"  NAND not found in scan, using W1={W1}, W2={W2}")

W_val, W_idx, nNeigh = sq_lattice_sparse(positions, W1, W2)

outputs_4, occ_4, E_4, Esite_4, code, logic_name = eval_logic_table(
    solver, Esite_base, W_val, W_idx, nNeigh, nSite,
    positions, input_positions, input_neighbors, output_site, W1, W2)

print(f"  W1={W1}  W2={W2}  → logic: {logic_name}  outputs={outputs_4.tolist()}")
print(f"  Energies: {[f'{e:.3f}' for e in E_4]}")
print(f"  On-site energies per input combo (eps = base + input_bias):")
for k in range(4):
    eps_str = '[' + ', '.join([f'{e:.2f}' for e in Esite_4[k]]) + ']'
    print(f"    In({INPUT_COMBOS[k]}) eps={eps_str}")
print(f"  Occupancy masks: {[bin(int(sum(o[i]<<i for i in range(len(o))))) for o in occ_4]}")

# Original visualization (backup)
fname1_orig = os.path.join(OUT_DIR, f'test1_groundstates_orig.{PLOT_EXT}')
plot_ground_states(positions, occ_4, E_4, outputs_4, input_positions,
                   output_site, logic_name,
                   fname=fname1_orig)
print(f"  Saved (original style): {fname1_orig}")

# New Tetris-style visualization
fname1 = os.path.join(OUT_DIR, f'test1_groundstates.{PLOT_EXT}')
plot_ground_states(positions, occ_4, E_4, outputs_4, input_positions,
                   output_site, logic_name,
                   input_neighbors=input_neighbors, Esite_4=Esite_4,
                   W1=W1, W2=W2, eps0=float(Esite_base[0]),
                   fname=fname1, tetris_style=True)
print(f"  Saved (Tetris style): {fname1}")


# ======================================================================
#  TEST 1b – Zigzag cluster at specific W1=0.5, W2=1
# ======================================================================
print("\n=== TEST 1b: Zigzag cluster at W1=0.5, W2=1 ===")

positions_z, input_positions_z, input_neighbors_z, output_site_z = make_cluster_zigzag()
nSite_z = len(positions_z)
Esite_base_z = -np.ones(nSite_z, dtype=np.float32)  # eps0 = -1

W1_z, W2_z = 0.5, 1.0
W_val_z, W_idx_z, nNeigh_z = sq_lattice_sparse(positions_z, W1_z, W2_z)

outputs_4_z, occ_4_z, E_4_z, Esite_4_z, code_z, logic_name_z = eval_logic_table(
    solver, Esite_base_z, W_val_z, W_idx_z, nNeigh_z, nSite_z,
    positions_z, input_positions_z, input_neighbors_z, output_site_z, W1_z, W2_z)

print(f"  W1={W1_z}  W2={W2_z}  → logic: {logic_name_z}  outputs={outputs_4_z.tolist()}")
print(f"  Energies: {[f'{e:.3f}' for e in E_4_z]}")
print(f"  On-site energies per input combo:")
for k in range(4):
    eps_str = '[' + ', '.join([f'{e:.2f}' for e in Esite_4_z[k]]) + ']'
    print(f"    In({INPUT_COMBOS[k]}) eps={eps_str}")

fname1b = os.path.join(OUT_DIR, f'test1b_zigzag_W1_0.5_W2_1.{PLOT_EXT}')
plot_ground_states(positions_z, occ_4_z, E_4_z, outputs_4_z, input_positions_z,
                   output_site_z, logic_name_z,
                   input_neighbors=input_neighbors_z, Esite_4=Esite_4_z,
                   W1=W1_z, W2=W2_z, eps0=float(Esite_base_z[0]),
                   fname=fname1b, tetris_style=True)
print(f"  Saved: {fname1b}")


# ======================================================================
#  DEBUG TEST – Check kernel stability for T_extended_inputs
# ======================================================================
print("\n=== DEBUG: T_extended_inputs stability check ===")

pos_d, inp_pos_d, inp_neigh_d, out_site_d = make_cluster_T_extended_inputs()
ns_d = len(pos_d)
Eb_d = -np.ones(ns_d, dtype=np.float32)

# Test at specific W1, W2 where NAND was found
W1_test, W2_test = 2.0, 0.5
Wv_d, Wi_d, Wn_d = sq_lattice_sparse(positions=pos_d, W1=W1_test, W2=W2_test, nSite=ns_d)

print(f"  Testing at W1={W1_test}, W2={W2_test}")
print(f"  nSite={ns_d}, positions shape={pos_d.shape}")
print(f"  W_val non-zero entries: {np.count_nonzero(Wv_d)}")

# Run solver 5 times and check consistency
results_debug = []
for run in range(5):
    outputs_4_d, occ_4_d, E_4_d, Esite_4_d, code_d, name_d = eval_logic_table(
        solver, Eb_d, Wv_d, Wi_d, Wn_d, ns_d,
        pos_d, inp_pos_d, inp_neigh_d, out_site_d, W1_test, W2_test)
    results_debug.append({
        'logic': name_d,
        'code': code_d,
        'outputs': outputs_4_d.copy(),
        'E': E_4_d.copy(),
        'occ': [o.copy() for o in occ_4_d]
    })
    print(f"  Run {run+1}: {name_d} outputs={outputs_4_d.tolist()} E={[f'{e:.3f}' for e in E_4_d]}")

# Check consistency
logic_codes = [r['code'] for r in results_debug]
if len(set(logic_codes)) > 1:
    print(f"  WARNING: Inconsistent results! Logic codes: {logic_codes}")
else:
    print(f"  OK: Consistent logic ({results_debug[0]['logic']}) across all runs")

# Print detailed energy/occupancy for first run
print(f"\n  Detailed energies (Run 1):")
for k in range(4):
    occ_mask = 0
    for i, occ in enumerate(results_debug[0]['occ'][k]):
        if occ:
            occ_mask |= (1 << i)
    print(f"    In({INPUT_COMBOS[k]}): E={results_debug[0]['E'][k]:.4f}  occ_mask=0b{occ_mask:09b}")

# Check scan stability
print(f"\n  Scan stability test (small grid, 3 runs):")
nW_test = 10
W1v_test = np.linspace(0.5, 2.5, nW_test)
W2v_test = np.linspace(0.1, 1.5, nW_test)
scan_results = []
for run in range(3):
    lmap, _, _ = scan_W1_W2(solver, pos_d, inp_pos_d, Eb_d, inp_neigh_d, out_site_d,
                            W1v_test, W2v_test, shift=-0.5)
    nand_count = np.sum(lmap == 7)
    true_count = np.sum(lmap == 15)
    scan_results.append({'nand': nand_count, 'true': true_count, 'lmap': lmap.copy()})
    print(f"    Run {run+1}: NAND={nand_count} cells, TRUE={true_count} cells")

# Check if all runs match
if all(np.array_equal(scan_results[0]['lmap'], r['lmap']) for r in scan_results[1:]):
    print(f"  OK: All scan runs identical")
else:
    print(f"  WARNING: Scan runs differ!")
    # Find first differing cell
    for i in range(nW_test):
        for j in range(nW_test):
            vals = [r['lmap'][i,j] for r in scan_results]
            if len(set(vals)) > 1:
                print(f"    First diff at ({i},{j}): {vals}")
                break
        else:
            continue
        break

# Permutation test - shuffle instance order and verify results identical
print(f"\n  Permutation test (detect race conditions/memory issues):")
nInst = 20  # Small batch
Esite_perm = np.random.randn(nInst, ns_d).astype(np.float32) * 0.5 - 1.0  # Random eps around -1
Wv_perm, Wi_perm, Wn_perm = sq_lattice_sparse(positions=pos_d, W1=2.0, W2=0.5, nSite=ns_d)

# Replicate W for batch
Wv_batch = np.tile(Wv_perm, (nInst, 1, 1))
Wi_batch = np.tile(Wi_perm, (nInst, 1, 1))
Wn_batch = np.tile(Wn_perm, (nInst, 1))

# Run 1: original order
E_orig, occ_orig = solver.solve_batch_W(Esite_perm, Wv_batch, Wi_batch, Wn_batch, ns_d)

# Run 2: shuffled order
perm = np.random.permutation(nInst)
E_shuf, occ_shuf = solver.solve_batch_W(Esite_perm[perm], Wv_batch[perm], Wi_batch[perm], Wn_batch[perm], ns_d)

# Unshuffle and compare
E_unshuf = E_shuf[np.argsort(perm)]
occ_unshuf = occ_shuf[np.argsort(perm)]

if np.allclose(E_orig, E_unshuf) and np.array_equal(occ_orig, occ_unshuf):
    print(f"  OK: Permutation test passed (results identical regardless of instance order)")
else:
    print(f"  WARNING: Permutation test FAILED!")
    diff_E = np.abs(E_orig - E_unshuf)
    diff_occ = occ_orig != occ_unshuf
    print(f"    Max energy diff: {np.max(diff_E):.6f}")
    print(f"    Different occupancy masks: {np.sum(diff_occ)} / {nInst}")
    for i in np.where(diff_occ)[0][:3]:
        print(f"    Instance {i}: orig=0b{occ_orig[i]:09b}, permuted=0b{occ_unshuf[i]:09b}")

# Python brute-force verification for single instance
print(f"\n  Python brute-force verification:")
def compute_energy_py(occ_mask, Esite, W_val, W_idx, nNeigh, nSite):
    """Compute energy for a given occupancy mask (Python version)."""
    E = 0.0
    for i in range(nSite):
        if not ((occ_mask >> i) & 1):
            continue
        E += Esite[i]
        for k in range(nNeigh[i]):
            j = W_idx[i, k]
            if j < i and ((occ_mask >> j) & 1):  # j<i to count each pair once
                E += W_val[i, k]
    return E

# Find ground state by brute force
print(f"    Brute-forcing all 2^{ns_d}={1<<ns_d} states...")
E_min_py = float('inf')
occ_min_py = 0
for occ in range(1 << ns_d):
    E = compute_energy_py(occ, Esite_perm[0], Wv_perm, Wi_perm, Wn_perm, ns_d)
    if E < E_min_py:
        E_min_py = E
        occ_min_py = occ

print(f"    Python: E_min={E_min_py:.4f}, occ=0b{occ_min_py:09b}")
print(f"    GPU:    E_min={E_orig[0]:.4f}, occ=0b{occ_orig[0]:09b}")
if abs(E_min_py - E_orig[0]) < 1e-4 and occ_min_py == occ_orig[0]:
    print(f"    OK: GPU matches Python brute-force")
else:
    print(f"    WARNING: GPU/CPU mismatch!")

# Compare two GPU kernels on same data
print(f"\n  Kernel comparison (mqca_groundstate vs mqca_groundstate_batch_W):")
E_single, occ_single = solver.solve(Esite_perm[0:1], Wv_perm, Wi_perm, Wn_perm, ns_d)
E_batch, occ_batch = solver.solve_batch_W(Esite_perm[0:1], Wv_batch[0:1], Wi_batch[0:1], Wn_batch[0:1], ns_d)
print(f"    Single kernel: E={E_single[0]:.4f}, occ=0b{occ_single[0]:09b}")
print(f"    Batch kernel:  E={E_batch[0]:.4f}, occ=0b{occ_batch[0]:09b}")
if abs(E_single[0] - E_batch[0]) < 1e-4 and occ_single[0] == occ_batch[0]:
    print(f"    OK: Both kernels agree")
else:
    print(f"    WARNING: Kernels disagree!")

# Test top8 kernel for degeneracy detection
print(f"\n  Top-8 kernel test (ground state uniqueness):")
Esite_top8_test = Esite_perm[0:1]  # Single instance
Wv_t8, Wi_t8, Wn_t8 = sq_lattice_sparse(positions=pos_d, W1=2.0, W2=0.5, nSite=ns_d)
Wv_batch_t8 = np.tile(Wv_t8, (1, 1, 1))
Wi_batch_t8 = np.tile(Wi_t8, (1, 1, 1))
Wn_batch_t8 = np.tile(Wn_t8, (1, 1))

E_top8, occ_top8 = solver.solve_batch_W_top8(Esite_top8_test, Wv_batch_t8, Wi_batch_t8, Wn_batch_t8, ns_d)
print(f"    Top 8 states:")
for k in range(8):
    print(f"      {k+1}. E={E_top8[0,k]:.4f}, occ=0b{occ_top8[0,k]:09b}")

# Check if ground state is well-defined (separated from first excited state)
ground_gap = E_top8[0, 1] - E_top8[0, 0]
print(f"    Ground state gap (E1 - E0): {ground_gap:.4f}")
if ground_gap < 0.01:
    print(f"    WARNING: Ground state is near-degenerate with first excited state!")
else:
    print(f"    OK: Ground state is well-defined (unique and separated)")
print(f"    Note: Excited state degeneracy among higher states is not a problem")

# Analyze energy gaps near phase boundary
print(f"\n  Energy gap analysis near phase boundary:")
# Pick a point where we saw the transition: around row 7, col 7
W1_test_gap = W1v_test[7]  # Use scan test values
W2_test_gap = W2v_test[7]
print(f"    Testing at W1={W1_test_gap:.2f}, W2={W2_test_gap:.2f} (transition region)")

# Build Esite for this (W1,W2)
Esite_gap = np.zeros((4, ns_d), dtype=np.float32)
for k, (A, B) in enumerate(INPUT_COMBOS):
    Esite_gap[k] = apply_input_bias(Eb_d, pos_d, inp_pos_d, inp_neigh_d, [A, B], W1_test_gap, W2_test_gap)

Wv_gap, Wi_gap, Wn_gap = sq_lattice_sparse(positions=pos_d, W1=W1_test_gap, W2=W2_test_gap, nSite=ns_d)

# Find top 5 states by brute force
states = []
for occ in range(1 << ns_d):
    E = compute_energy_py(occ, Esite_gap[0], Wv_gap, Wi_gap, Wn_gap, ns_d)  # Input (0,0)
    states.append((E, occ))
states.sort()

print(f"    Top 5 states for input (0,0):")
for i, (E, occ) in enumerate(states[:5]):
    print(f"      {i+1}. E={E:.4f}, occ=0b{occ:09b}")
print(f"    Energy gap: {states[1][0] - states[0][0]:.4f}")


# ======================================================================
#  TEST 1c – 4 T-shape variations with W1>0, W2>0, W1>W2
# ======================================================================
print("\n=== TEST 1c: 4 T-shape variations (W1>0, W2>0, W1>W2) ===")

T_variants = {
    'T_extended_inputs': make_cluster_T_extended_inputs,
    'T_extended_output': make_cluster_T_extended_output,
    'T_no_center': make_cluster_T_no_center,
    'T_fork': make_cluster_T_fork,
}

# Scan parameters
nW_f = 40
W1_vals_f = np.linspace(0.1, 3.0, nW_f)
W2_vals_f = np.linspace(0.1, 2.5, nW_f)

# Create 2x2 subplot for phase diagrams
fig, axes = plt.subplots(2, 2, figsize=(14, 12))
fig.suptitle('T-shape Variations | eps0=-1 | W1>0, W2>0, W1>W2 | Degenerate ground states masked (red hatched)', fontsize=14)

results = {}

for idx, (name, maker) in enumerate(T_variants.items()):
    ax = axes.flat[idx]
    row, col = idx // 2, idx % 2

    pos, inp_pos, inp_neigh, out_site = maker()
    ns = len(pos)
    Eb = -np.ones(ns, dtype=np.float32)

    # Scan with top8 kernel to detect degenerate ground states
    lmap, degenerate_mask = scan_W1_W2_top8(solver, pos, inp_pos, Eb, inp_neigh, out_site,
                                            W1_vals_f, W2_vals_f, nSite=ns, degeneracy_threshold=0.01, shift=-0.5)

    # Mask W1 <= W2 region
    W1_grid, W2_grid = np.meshgrid(W1_vals_f, W2_vals_f)
    mask_invalid = W1_grid <= W2_grid

    # Build image
    img = np.ones((nW_f, nW_f, 3)) * 0.9
    for code in range(16):
        c = np.array(LOGIC_COLORS[code][:3])
        mask = (lmap == code) & ~mask_invalid & ~degenerate_mask
        img[mask] = c
    img[mask_invalid] = [0.7, 0.7, 0.7]  # Darker gray for invalid (W1 <= W2)

    ax.imshow(img, origin='lower', extent=[W1_vals_f[0], W1_vals_f[-1], W2_vals_f[0], W2_vals_f[-1]])

    # Overlay hazard pattern for degenerate ground states
    ax.contour(W1_grid, W2_grid, degenerate_mask.astype(float), levels=[0.5], colors='red', linewidths=2)
    ax.contourf(W1_grid, W2_grid, degenerate_mask.astype(float), levels=[0.5, 1.5], colors='red', alpha=0.3)

    ax.plot([0, 3], [0, 3], 'k--', lw=1, alpha=0.5)
    ax.set_xlim(0, 3)
    ax.set_ylim(0, 2.5)
    ax.set_xlabel('W1 (Cartesian)', fontsize=9)
    ax.set_ylabel('W2 (Diagonal)', fontsize=9)
    ax.set_title(name.replace('_', ' '), fontsize=10)

    # Count logic functions (excluding degenerate areas)
    valid_mask = ~mask_invalid & ~degenerate_mask
    valid_codes = lmap[valid_mask]
    unique_codes = np.unique(valid_codes)
    logic_list = [LOGIC_NAMES[c] for c in unique_codes if c in LOGIC_NAMES]
    results[name] = logic_list
    useful = [l for l in logic_list if l in USEFUL_LOGIC]
    n_degenerate = np.sum(degenerate_mask & ~mask_invalid)
    print(f"  {name:20s}: {logic_list}")
    print(f"    Useful: {useful if useful else 'none'}")
    print(f"    Degenerate ground states: {n_degenerate} / {np.sum(~mask_invalid)} cells")

    # Print small ASCII version of logic map for T_extended_inputs
    if name == 'T_extended_inputs' and nW_f <= 40:
        print(f"    Raw logic map (codes: .=invalid, N=NAND(7), T=TRUE(15)):")
        symbols = {7: 'N', 15: 'T', -1: '.'}
        for iy in range(min(20, nW_f)):
            row = ""
            for ix in range(min(20, nW_f)):
                code = lmap[iy, ix]
                if mask_invalid[iy, ix]:
                    row += "."
                elif code == 7:
                    row += "N"
                elif code == 15:
                    row += "T"
                else:
                    row += str(code) if code < 10 else "X"
            print(f"      {row}")

# Overall legend
handles = [mpatches.Patch(color=LOGIC_COLORS[c], label=LOGIC_NAMES[c])
           for c in range(16) if c in LOGIC_NAMES]
fig.legend(handles=handles, loc='center right', fontsize=7, title='Logic', ncol=1, bbox_to_anchor=(1.02, 0.5))

plt.tight_layout(rect=[0, 0, 0.95, 1])
fname1c = os.path.join(OUT_DIR, f'test1c_T_variants_comparison.{PLOT_EXT}')
plt.savefig(fname1c, dpi=150)
plt.close()
print(f"  Saved: {fname1c}")


# ======================================================================
#  TEST 2 – 2-D (W1, W2) logic phase diagram for cross cluster
# ======================================================================
print("\n=== TEST 2: W1-W2 phase diagram ===")

nW = 60
W1_vals = np.linspace(0.0, 3.0, nW)  # Only positive W1
W2_vals = np.linspace(0.0, 3.0, nW)  # Only positive W2

logic_map, occ_map, E_map = scan_W1_W2(
    solver, positions, input_positions, Esite_base, input_neighbors, output_site,
    W1_vals, W2_vals, shift=-0.5)

# Print how many (W1,W2) points give each logic type
print("  Logic function counts:")
for code in np.unique(logic_map):
    cnt = np.sum(logic_map == code)
    print(f"    {LOGIC_NAMES[code]:8s}  {cnt:5d}  ({100.*cnt/nW**2:.1f}%)")

fname2a = os.path.join(OUT_DIR, f'test2_logic_map.{PLOT_EXT}')
plot_logic_map(W1_vals, W2_vals, logic_map,
               title=f'Cross cluster  |  eps0=-1  Logic map',
               fname=fname2a)
print(f"  Saved: {fname2a}")

fname2b = os.path.join(OUT_DIR, f'test2_useful_logic_map.{PLOT_EXT}')
plot_logic_fraction_map(W1_vals, W2_vals, logic_map,
                        title='Cross cluster  |  Useful logic regions (green)',
                        fname=fname2b)
print(f"  Saved: {fname2b}")


# ======================================================================
#  TEST 3 – scan multiple cluster geometries over (W1,W2), find which
#            ones implement useful logic functions
# ======================================================================
print("\n=== TEST 3: geometry scan for useful logic ===")

clusters = {
    'Cross'   : make_cluster_cross(),
    'T-shape' : make_cluster_T(),
    'L-shape' : make_cluster_L(),
    'S-shape' : make_cluster_S(),
    'Zigzag'  : make_cluster_zigzag(),
    'Chain'   : make_cluster_chain(),
    'Straight_center' : make_cluster_straight_line_center(),
}

nW3 = 40   # coarser grid for geometry scan
W1_scan = np.linspace(0.1, 3.0, nW3)  # Only positive W1
W2_scan = np.linspace(0.1, 3.0, nW3)  # Only positive W2

fig, axes = plt.subplots(len(clusters), 2, figsize=(12, 4*len(clusters)))
fig.suptitle('Geometry scan: logic phase diagrams', fontsize=13)

useful_codes = {c for c, n in LOGIC_NAMES.items() if n in USEFUL_LOGIC}

results_table = []

for row, (name, (pos, inp_pos, inp_neigh, out_site)) in enumerate(clusters.items()):
    ns = len(pos)
    Eb = -np.ones(ns, dtype=np.float32)  # eps0 = -1
    lmap, _, _ = scan_W1_W2(solver, pos, inp_pos, Eb, inp_neigh, out_site,
                              W1_scan, W2_scan, shift=-0.5)
    found_logic = [LOGIC_NAMES[c] for c in np.unique(lmap)]
    useful_found = [l for l in found_logic if l in USEFUL_LOGIC]
    print(f"  {name:12s}: logic found = {found_logic}")
    print(f"             useful = {useful_found}")
    results_table.append((name, useful_found))

    ax_map  = axes[row, 0]
    ax_use  = axes[row, 1]

    # phase diagram
    import matplotlib.patches as mpatches
    from pyBall.OCL.MQCA import LOGIC_COLORS
    nW2_, nW1_ = lmap.shape
    img = np.zeros((nW2_, nW1_, 3))
    for code in range(16):
        mask = lmap == code
        img[mask] = LOGIC_COLORS[code][:3]
    extent = [W1_scan[0], W1_scan[-1], W2_scan[0], W2_scan[-1]]
    ax_map.imshow(img, origin='lower', extent=extent, aspect='auto', interpolation='nearest')
    ax_map.set_title(f'{name}  –  logic map', fontsize=10)
    ax_map.set_xlabel('W1'); ax_map.set_ylabel('W2')
    present = np.unique(lmap)
    patches = [mpatches.Patch(color=LOGIC_COLORS[c], label=LOGIC_NAMES[c]) for c in present]
    ax_map.legend(handles=patches, fontsize=7, loc='upper right', framealpha=0.7)

    # useful fraction
    useful_mask = np.isin(lmap, list(useful_codes))
    ax_use.imshow(useful_mask.astype(float), origin='lower', extent=extent,
                  aspect='auto', cmap='RdYlGn', vmin=0, vmax=1, interpolation='nearest')
    frac = useful_mask.mean()*100
    ax_use.set_title(f'{name}  –  useful ({frac:.1f}%)', fontsize=10)
    ax_use.set_xlabel('W1'); ax_use.set_ylabel('W2')

plt.tight_layout()
fname3 = os.path.join(OUT_DIR, f'test3_geometry_scan.{PLOT_EXT}')
plt.savefig(fname3, dpi=150)
plt.close()
print(f"  Saved: {fname3}")

# ---- print summary table ----
print("\n=== SUMMARY ===")
print(f"{'Cluster':12s}  {'Useful logic functions found'}")
print('-'*50)
for name, useful in results_table:
    print(f"  {name:12s}  {useful if useful else 'none'}")

print("\nAll tests completed successfully.")
