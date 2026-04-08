#!/usr/bin/env python
"""Analyze degeneracy at specific (W1,W2) for T_extended_output cluster."""
import sys, os
sys.path.insert(0, '.')

from pyBall.OCL.MQCA import MQCASolver, sq_lattice_sparse, apply_input_bias, INPUT_COMBOS
import numpy as np

print("Creating solver...")
solver = MQCASolver(preferred_vendor='nvidia', bPrint=False)

# T_extended_output cluster
def make_cluster_T_extended_output():
    positions = np.array([
        [0,3], [2,3],   # under pads (sites 0,1)
        [0,2], [1,2], [2,2],  # crossbar (sites 2,3,4)
        [1,1],          # upper stem (site 5)
        [1,0],          # lower stem (site 6)
        [1,-1],         # output (site 7)
    ], dtype=int)
    input_positions = [(-1, 3), (3, 3)]
    input_neighbors = [[0], [1]]  # site 0 neighbors A pad, site 1 neighbors B pad
    output_site = 7
    return positions, input_positions, input_neighbors, output_site

pos, inp_pos, inp_neigh, out_site = make_cluster_T_extended_output()
nSite = len(pos)
Eb = -np.ones(nSite, dtype=np.float32)

W1, W2 = 2.0, 1.0
print(f"\n=== Analyzing T_extended_output at W1={W1}, W2={W2} ===\n")

# Build Esite for each input combo
from pyBall.OCL.MQCA import compute_input_bias
Esite_batch = np.zeros((4, nSite), dtype=np.float32)
for k, (A, B) in enumerate(INPUT_COMBOS):
    print(f"\nDEBUG Input ({A}, {B}):")
    print(f"  positions = {pos}")
    print(f"  input_positions = {inp_pos}")
    print(f"  input_neighbors = {inp_neigh}")
    print(f"  input_vals = [{A}, {B}]")
    print(f"  W1={W1}, W2={W2}")
    bias = compute_input_bias(pos, inp_pos, inp_neigh, [A, B], W1, W2)
    Esite_batch[k] = Eb + bias
    print(f"  bias = {bias}")
    print(f"  Esite = {Esite_batch[k]}")

# Build W
W_val, W_idx, nNeigh = sq_lattice_sparse(positions=pos, W1=W1, W2=W2, nSite=nSite)

# Reshape for batch
W_val_batch = np.tile(W_val, (4, 1, 1))
W_idx_batch = np.tile(W_idx, (4, 1, 1))
nNeigh_batch = np.tile(nNeigh, (4, 1))

# Get top8 states
E_top8, occ_top8 = solver.solve_batch_W_top8(Esite_batch, W_val_batch, W_idx_batch, nNeigh_batch, nSite)

print("Top 8 states for each input combination:")
print("=" * 80)

for k, (A, B) in enumerate(INPUT_COMBOS):
    print(f"\nInput ({A}, {B}):")
    print(f"  On-site energies: {Esite_batch[k]}")
    print(f"  Top 8 states:")
    for i in range(8):
        occ_mask = occ_top8[k, i]
        occ_array = [(occ_mask >> s) & 1 for s in range(nSite)]
        output_val = (occ_mask >> out_site) & 1
        print(f"    {i+1}. E={E_top8[k,i]:.6f}, output={output_val}, occ=0b{occ_mask:09b} = {occ_array}")
    
    # Check ground state gap
    ground_gap = E_top8[k, 1] - E_top8[k, 0]
    print(f"  Ground state gap: {ground_gap:.6f}")
    if ground_gap < 0.01:
        print(f"  >>> DEGENERATE: Ground state is not well-defined!")
    else:
        print(f"  >>> OK: Ground state is well-defined")

# Summary
print("\n" + "=" * 80)
print("SUMMARY:")
print("=" * 80)
degenerate_count = 0
for k in range(4):
    ground_gap = E_top8[k, 1] - E_top8[k, 0]
    if ground_gap < 0.01:
        degenerate_count += 1
        print(f"Input {INPUT_COMBOS[k]}: DEGENERATE (gap={ground_gap:.6f})")
    else:
        print(f"Input {INPUT_COMBOS[k]}: OK (gap={ground_gap:.6f})")

if degenerate_count > 0:
    print(f"\n>>> {degenerate_count}/4 input combinations have degenerate ground states")
    print(">>> This cluster is unreliable at this (W1,W2) point")
else:
    print(f"\n>>> All input combinations have well-defined ground states")
    print(">>> This cluster is reliable at this (W1,W2) point")
