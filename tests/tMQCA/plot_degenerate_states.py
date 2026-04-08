#!/usr/bin/env python
"""Plot degenerate ground states for T_extended_output at W1=2.0, W2=1.0."""
import sys, os
sys.path.insert(0, '.')

from pyBall.OCL.MQCA import MQCASolver, sq_lattice_sparse, apply_input_bias, INPUT_COMBOS, plot_ground_states
import numpy as np
import matplotlib.pyplot as plt

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
    input_neighbors = [[0], [1]]
    output_site = 7
    return positions, input_positions, input_neighbors, output_site

pos, inp_pos, inp_neigh, out_site = make_cluster_T_extended_output()
nSite = len(pos)
Eb = -np.ones(nSite, dtype=np.float32)

W1, W2 = 2.0, 1.0

# Build Esite and W for each input combo
Esite_batch = np.zeros((4, nSite), dtype=np.float32)
for k, (A, B) in enumerate(INPUT_COMBOS):
    Esite_batch[k] = apply_input_bias(Eb, pos, inp_pos, inp_neigh, [A, B], W1, W2)

W_val, W_idx, nNeigh = sq_lattice_sparse(positions=pos, W1=W1, W2=W2, nSite=nSite)

# Get top8 states
W_val_batch = np.tile(W_val, (4, 1, 1))
W_idx_batch = np.tile(W_idx, (4, 1, 1))
nNeigh_batch = np.tile(nNeigh, (4, 1))

E_top8, occ_top8 = solver.solve_batch_W_top8(Esite_batch, W_val_batch, W_idx_batch, nNeigh_batch, nSite)

# Plot degenerate states for input (0,1) and (1,0)
fig, axes = plt.subplots(2, 6, figsize=(24, 8))
fig.suptitle(f'T_extended_output at W1={W1}, W2={W2} - Degenerate Ground States (E=-3.0)', fontsize=16)

# Input (0,1) - 6 degenerate states at E=-3.0
k = 1  # input (0,1)
for i in range(6):
    ax = axes[0, i]
    occ_mask = occ_top8[k, i]
    occ_array = [(occ_mask >> s) & 1 for s in range(nSite)]
    output_val = (occ_mask >> out_site) & 1
    ax.set_title(f'State {i+1}: out={output_val}\nocc={occ_array}', fontsize=9)
    
    # Plot sites as grid
    for s, (x, y) in enumerate(pos):
        color = 'red' if s == out_site else ('blue' if occ_array[s] else 'lightgray')
        marker = 's' if s in inp_neigh[0] + inp_neigh[1] else 'o'
        ax.scatter(x, y, s=200, c=color, marker=marker, edgecolors='black', linewidth=2)
        ax.text(x, y+0.2, str(s), ha='center', fontsize=8)
    
    # Plot input pads
    ax.scatter(inp_pos[0][0], inp_pos[0][1], s=150, c='green', marker='^', edgecolors='black', linewidth=2)
    ax.scatter(inp_pos[1][0], inp_pos[1][1], s=150, c='green', marker='^', edgecolors='black', linewidth=2)
    ax.text(inp_pos[0][0], inp_pos[0][1]+0.2, 'A', ha='center', fontsize=8)
    ax.text(inp_pos[1][0], inp_pos[1][1]+0.2, 'B', ha='center', fontsize=8)
    
    ax.set_aspect('equal')
    ax.set_xlim(-2, 4)
    ax.set_ylim(-2, 4)
    ax.grid(True, alpha=0.3)
    ax.set_xlabel('Input (0,1): A=0, B=1', fontsize=8)

# Input (1,0) - 6 degenerate states at E=-3.0
k = 2  # input (1,0)
for i in range(6):
    ax = axes[1, i]
    occ_mask = occ_top8[k, i]
    occ_array = [(occ_mask >> s) & 1 for s in range(nSite)]
    output_val = (occ_mask >> out_site) & 1
    ax.set_title(f'State {i+1}: out={output_val}\nocc={occ_array}', fontsize=9)
    
    # Plot sites as grid
    for s, (x, y) in enumerate(pos):
        color = 'red' if s == out_site else ('blue' if occ_array[s] else 'lightgray')
        marker = 's' if s in inp_neigh[0] + inp_neigh[1] else 'o'
        ax.scatter(x, y, s=200, c=color, marker=marker, edgecolors='black', linewidth=2)
        ax.text(x, y+0.2, str(s), ha='center', fontsize=8)
    
    # Plot input pads
    ax.scatter(inp_pos[0][0], inp_pos[0][1], s=150, c='green', marker='^', edgecolors='black', linewidth=2)
    ax.scatter(inp_pos[1][0], inp_pos[1][1], s=150, c='green', marker='^', edgecolors='black', linewidth=2)
    ax.text(inp_pos[0][0], inp_pos[0][1]+0.2, 'A', ha='center', fontsize=8)
    ax.text(inp_pos[1][0], inp_pos[1][1]+0.2, 'B', ha='center', fontsize=8)
    
    ax.set_aspect('equal')
    ax.set_xlim(-2, 4)
    ax.set_ylim(-2, 4)
    ax.grid(True, alpha=0.3)
    ax.set_xlabel('Input (1,0): A=1, B=0', fontsize=8)

plt.tight_layout()
plt.savefig('/home/prokophapala/git/FireCore/analyze_degeneracy_plot.png', dpi=150)
print("Saved: /home/prokophapala/git/FireCore/analyze_degeneracy_plot.png")
