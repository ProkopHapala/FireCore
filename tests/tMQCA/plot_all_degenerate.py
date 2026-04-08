#!/usr/bin/env python
"""Plot all 4 input combinations with degenerate states using Tetris view mode."""
import sys, os
sys.path.insert(0, '.')

from pyBall.OCL.MQCA import MQCASolver, sq_lattice_sparse, apply_input_bias, INPUT_COMBOS
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
    input_positions = [(0, 4), (2, 4)]  # MOVED: inputs above sites instead of from sides
    input_neighbors = [[0], [1]]
    output_site = 7
    return positions, input_positions, input_neighbors, output_site

pos, inp_pos, inp_neigh, out_site = make_cluster_T_extended_output()
nSite = len(pos)
Eb = -np.ones(nSite, dtype=np.float32)

W1, W2 = 0.0, 0.0  # Set W1=0, W2=0 (no coupling)

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

# Plot using Tetris-style visualization
fig, axes = plt.subplots(2, 4, figsize=(20, 10))
fig.suptitle(f'T_extended_output at W1={W1}, W2={W2} - All Input Combinations with Degenerate States (Tetris view)', fontsize=14)

# Tetris-style plotting function
def plot_tetris_state(ax, positions, occ_mask, output_site, input_positions, input_neighbors):
    """Plot a single state in Tetris style."""
    occ_array = [(occ_mask >> s) & 1 for s in range(len(positions))]
    
    # Draw grid
    ax.set_xlim(-2, 4)
    ax.set_ylim(-2, 4)
    ax.grid(True, alpha=0.3)
    
    # Draw sites as blocks
    for s, (x, y) in enumerate(positions):
        if occ_array[s]:
            # Occupied site - draw as filled block
            color = 'red' if s == output_site else 'blue'
            rect = plt.Rectangle((x-0.45, y-0.45), 0.9, 0.9, 
                                facecolor=color, edgecolor='black', linewidth=2)
            ax.add_patch(rect)
            ax.text(x, y, str(s), ha='center', va='center', 
                   color='white' if s == output_site else 'white', 
                   fontsize=10, fontweight='bold')
        else:
            # Empty site - draw as outline
            rect = plt.Rectangle((x-0.45, y-0.45), 0.9, 0.9, 
                                facecolor='lightgray', edgecolor='black', linewidth=1)
            ax.add_patch(rect)
            ax.text(x, y, str(s), ha='center', va='center', 
                   color='gray', fontsize=8)
    
    # Draw input pads
    ax.scatter(inp_pos[0][0], inp_pos[0][1], s=200, c='green', marker='^', 
              edgecolors='black', linewidth=2, zorder=10)
    ax.scatter(inp_pos[1][0], inp_pos[1][1], s=200, c='green', marker='^', 
              edgecolors='black', linewidth=2, zorder=10)
    ax.text(inp_pos[0][0], inp_pos[0][1]+0.3, 'A', ha='center', fontsize=10, fontweight='bold')
    ax.text(inp_pos[1][0], inp_pos[1][1]+0.3, 'B', ha='center', fontsize=10, fontweight='bold')
    
    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])

# Plot each input combination with its top 2 states (or more if degenerate)
for k, (A, B) in enumerate(INPUT_COMBOS):
    # Determine which states to plot (degenerate ones)
    ground_gap = E_top8[k, 1] - E_top8[k, 0]
    if ground_gap < 0.01:
        # Degenerate - plot top 2 states
        n_states_to_plot = 2
        degenerate_str = " (DEGENERATE)"
    else:
        # Non-degenerate - plot ground state only
        n_states_to_plot = 1
        degenerate_str = ""
    
    for i in range(n_states_to_plot):
        ax = axes[i, k]
        occ_mask = occ_top8[k, i]
        output_val = (occ_mask >> out_site) & 1
        ax.set_title(f'In({A},{B}) State {i+1}: E={E_top8[k,i]:.3f}, out={output_val}{degenerate_str}', 
                    fontsize=9)
        plot_tetris_state(ax, pos, occ_mask, out_site, inp_pos, inp_neigh)
    
    # If only 1 state plotted (non-degenerate), leave second row cell empty
    if n_states_to_plot == 1:
        axes[1, k].axis('off')

plt.tight_layout()
plt.savefig('/home/prokophapala/git/FireCore/plot_all_degenerate_tetris.png', dpi=150)
print("Saved: /home/prokophapala/git/FireCore/plot_all_degenerate_tetris.png")
