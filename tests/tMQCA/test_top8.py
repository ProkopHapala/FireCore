#!/usr/bin/env python
"""Quick test of top8 kernel."""
import sys, os
sys.path.insert(0, '.')

print("Importing...")
from pyBall.OCL.MQCA import MQCASolver, sq_lattice_sparse
import numpy as np

print("Creating solver...")
solver = MQCASolver(preferred_vendor='nvidia', bPrint=False)
print("Solver created successfully")

# Test with T_extended_inputs cluster
positions = np.array([
    [0,3], [2,3], [0,2], [2,2], [0,1], [1,1], [2,1], [1,0], [1,-1]
], dtype=int)
nSite = len(positions)

Esite = -np.ones((1, nSite), dtype=np.float32)
W_val, W_idx, nNeigh = sq_lattice_sparse(positions, W1=2.0, W2=0.5, nSite=nSite)

# Reshape for batch
W_val_batch = W_val.reshape(1, nSite, 8)
W_idx_batch = W_idx.reshape(1, nSite, 8)
nNeigh_batch = nNeigh.reshape(1, nSite)

print(f"Testing top8 kernel with nSite={nSite}...")
E_top8, occ_top8 = solver.solve_batch_W_top8(Esite, W_val_batch, W_idx_batch, nNeigh_batch, nSite)

print("Top 8 states:")
for k in range(8):
    print(f"  {k+1}. E={E_top8[0,k]:.4f}, occ=0b{occ_top8[0,k]:09b}")

# Check if ground state is well-defined (separated from first excited state)
ground_gap = E_top8[0, 1] - E_top8[0, 0]
print(f"Ground state gap (E1 - E0): {ground_gap:.4f}")
if ground_gap < 0.01:
    print("WARNING: Ground state is near-degenerate with first excited state!")
else:
    print("OK: Ground state is well-defined (unique and separated)")
print(f"Note: Excited state degeneracy among higher states is not a problem")
