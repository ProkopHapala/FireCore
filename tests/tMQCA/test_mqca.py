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
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

import numpy as np
import matplotlib.pyplot as plt

from pyBall.OCL.MQCA import (
    MQCASolver,
    sq_lattice_sparse,
    apply_input_bias,
    eval_logic_table,
    scan_W1_W2,
    identify_logic,
    occ_mask_to_array,
    plot_ground_states,
    plot_logic_map,
    plot_logic_fraction_map,
    INPUT_COMBOS, LOGIC_NAMES, USEFUL_LOGIC, MAX_NEIGH,
)

OUT_DIR = os.path.dirname(os.path.abspath(__file__))

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
print("\n=== TEST 1: single cluster ground states ===")

positions, input_positions, input_neighbors, output_site = make_cluster_T()
nSite = len(positions)
Esite_base = np.ones(nSite, dtype=np.float32)  # eps0 = 1 for all

W_in   = 1.0
# Try several (W1,W2) points and pick the most interesting one
test1_candidates = [(-1.5, 1.2), (-1.2, 0.3), (-1.0, -0.5), (-2.0, 0.5),
                    (-1.5, 0.5), (-1.3, 0.8), (1.0, 0.5), (-1.0, 1.5)]
best_code = 0;  best_params = test1_candidates[0]
for W1c, W2c in test1_candidates:
    Wv, Wi, Wn = sq_lattice_sparse(positions, W1c, W2c)
    o4, _, _, cc, cn = eval_logic_table(solver, Esite_base, Wv, Wi, Wn, nSite, input_neighbors, output_site, W_in=W_in)
    print(f"  Try W1={W1c}  W2={W2c}  → {cn}  outputs={o4.tolist()}")
    if cc not in (0, 15) and (best_code in (0,15) or cc == 7):   # prefer NAND, then any non-trivial
        best_code = cc;  best_params = (W1c, W2c)
W1, W2 = best_params
print(f"  Chosen W1={W1}  W2={W2}")
W_val, W_idx, nNeigh = sq_lattice_sparse(positions, W1, W2)

outputs_4, occ_4, E_4, code, logic_name = eval_logic_table(
    solver, Esite_base, W_val, W_idx, nNeigh, nSite,
    input_neighbors, output_site, W_in=W_in)

print(f"  W1={W1}  W2={W2}  → logic: {logic_name}  outputs={outputs_4.tolist()}")

fname1 = os.path.join(OUT_DIR, 'test1_groundstates.png')
plot_ground_states(positions, occ_4, E_4, outputs_4, input_positions,
                   output_site, logic_name,
                   W_val=W_val, W_idx=W_idx, nNeigh=nNeigh,
                   fname=fname1)
print(f"  Saved: {fname1}")


# ======================================================================
#  TEST 2 – 2-D (W1, W2) logic phase diagram for cross cluster
# ======================================================================
print("\n=== TEST 2: W1-W2 phase diagram ===")

nW = 60
W1_vals = np.linspace(-2.0, 2.0, nW)
W2_vals = np.linspace(-2.0, 2.0, nW)

logic_map, occ_map, E_map = scan_W1_W2(
    solver, positions, Esite_base, input_neighbors, output_site,
    W1_vals, W2_vals, W_in=W_in)

# Print how many (W1,W2) points give each logic type
print("  Logic function counts:")
for code in np.unique(logic_map):
    cnt = np.sum(logic_map == code)
    print(f"    {LOGIC_NAMES[code]:8s}  {cnt:5d}  ({100.*cnt/nW**2:.1f}%)")

fname2a = os.path.join(OUT_DIR, 'test2_logic_map.png')
plot_logic_map(W1_vals, W2_vals, logic_map,
               title=f'Cross cluster  |  eps0=1  W_in={W_in}  Logic map',
               fname=fname2a)
print(f"  Saved: {fname2a}")

fname2b = os.path.join(OUT_DIR, 'test2_useful_logic_map.png')
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
}

nW3 = 40   # coarser grid for geometry scan
W1_scan = np.linspace(-3.0, 3.0, nW3)
W2_scan = np.linspace(-3.0, 3.0, nW3)

fig, axes = plt.subplots(len(clusters), 2, figsize=(12, 4*len(clusters)))
fig.suptitle('Geometry scan: logic phase diagrams', fontsize=13)

useful_codes = {c for c, n in LOGIC_NAMES.items() if n in USEFUL_LOGIC}

results_table = []

for row, (name, (pos, inp_pos, inp_neigh, out_site)) in enumerate(clusters.items()):
    ns = len(pos)
    Eb = np.ones(ns, dtype=np.float32)
    lmap, _, _ = scan_W1_W2(solver, pos, Eb, inp_neigh, out_site,
                              W1_scan, W2_scan, W_in=1.0)
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
fname3 = os.path.join(OUT_DIR, 'test3_geometry_scan.png')
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
