#!/usr/bin/env python
"""
CLI tool to analyze MQCA degeneracy for different cluster geometries and parameters.

Usage:
    python analyze_mqca_degeneracy.py --cluster T_extended_output --W1 2.0 --W2 0.0
    python analyze_mqca_degeneracy.py --cluster T_extended_inputs --W1 2.0 --W2 1.0 --input-pos side
"""
import sys, os
sys.path.insert(0, '.')

import argparse
from pyBall.OCL.MQCA import MQCASolver, sq_lattice_sparse, apply_input_bias, INPUT_COMBOS
import numpy as np
import matplotlib.pyplot as plt

# Cluster definitions
def make_cluster_T_extended_output(input_pos_type='above'):
    """T-shape with longer output line (2-site stem)."""
    positions = np.array([
        [0,3], [2,3],   # under pads (sites 0,1)
        [0,2], [1,2], [2,2],  # crossbar (sites 2,3,4)
        [1,1],          # upper stem (site 5)
        [1,0],          # lower stem (site 6)
        [1,-1],         # output (site 7)
    ], dtype=int)
    if input_pos_type == 'above':
        input_positions = [(0, 4), (2, 4)]
    elif input_pos_type == 'side':
        input_positions = [(-1, 3), (3, 3)]
    else:
        raise ValueError(f"Unknown input_pos_type: {input_pos_type}")
    input_neighbors = [[0], [1]]
    output_site = 7
    return positions, input_positions, input_neighbors, output_site

def make_cluster_T_extended_inputs(input_pos_type='side'):
    """T-shape with extended input line (2-site input stem)."""
    positions = np.array([
        [0,3], [2,3],   # under pads (sites 0,1)
        [0,2], [2,2],   # crossbar (sites 2,3)
        [0,1], [1,1], [2,1],  # upper stem (sites 4,5,6)
        [1,0],          # lower stem (site 7)
        [1,-1],         # output (site 8)
    ], dtype=int)
    if input_pos_type == 'above':
        input_positions = [(0, 4), (2, 4)]
    elif input_pos_type == 'side':
        input_positions = [(-1, 3), (3, 3)]
    else:
        raise ValueError(f"Unknown input_pos_type: {input_pos_type}")
    input_neighbors = [[0], [1]]
    output_site = 8
    return positions, input_positions, input_neighbors, output_site

def make_cluster_T_simple(input_pos_type='side'):
    """Simple T-shape (5 sites)."""
    positions = np.array([
        [0,1], [1,1], [2,1],
        [1,0], [1,-1],
    ], dtype=int)
    if input_pos_type == 'above':
        input_positions = [(0, 2), (2, 2)]
    elif input_pos_type == 'side':
        input_positions = [(-1, 1), (3, 1)]
    else:
        raise ValueError(f"Unknown input_pos_type: {input_pos_type}")
    input_neighbors = [[0], [2]]
    output_site = 4
    return positions, input_positions, input_neighbors, output_site

CLUSTER_MAKERS = {
    'T_extended_output': make_cluster_T_extended_output,
    'T_extended_inputs': make_cluster_T_extended_inputs,
    'T_simple': make_cluster_T_simple,
}

def plot_tetris_state(ax, positions, occ_mask, output_site, input_positions):
    """Plot a single state in Tetris style."""
    occ_array = [(occ_mask >> s) & 1 for s in range(len(positions))]
    
    ax.set_xlim(-2, 5)
    ax.set_ylim(-2, 5)
    ax.grid(True, alpha=0.3)
    
    for s, (x, y) in enumerate(positions):
        if occ_array[s]:
            color = 'red' if s == output_site else 'blue'
            rect = plt.Rectangle((x-0.45, y-0.45), 0.9, 0.9, 
                                facecolor=color, edgecolor='black', linewidth=2)
            ax.add_patch(rect)
            ax.text(x, y, str(s), ha='center', va='center', 
                   color='white', fontsize=10, fontweight='bold')
        else:
            rect = plt.Rectangle((x-0.45, y-0.45), 0.9, 0.9, 
                                facecolor='lightgray', edgecolor='black', linewidth=1)
            ax.add_patch(rect)
            ax.text(x, y, str(s), ha='center', va='center', 
                   color='gray', fontsize=8)
    
    ax.scatter(input_positions[0][0], input_positions[0][1], s=200, c='green', marker='^', 
              edgecolors='black', linewidth=2, zorder=10)
    ax.scatter(input_positions[1][0], input_positions[1][1], s=200, c='green', marker='^', 
              edgecolors='black', linewidth=2, zorder=10)
    ax.text(input_positions[0][0], input_positions[0][1]+0.3, 'A', ha='center', fontsize=10, fontweight='bold')
    ax.text(input_positions[1][0], input_positions[1][1]+0.3, 'B', ha='center', fontsize=10, fontweight='bold')
    
    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])

def main():
    parser = argparse.ArgumentParser(description='Analyze MQCA degeneracy')
    parser.add_argument('--cluster', type=str, default='T_extended_output',
                       choices=list(CLUSTER_MAKERS.keys()),
                       help='Cluster geometry')
    parser.add_argument('--W1', type=float, default=2.0, help='Cartesian coupling')
    parser.add_argument('--W2', type=float, default=0.0, help='Diagonal coupling')
    parser.add_argument('--input-pos', type=str, default='above',
                       choices=['above', 'side'],
                       help='Input position: above sites or from sides')
    parser.add_argument('--shift', type=float, default=0.0,
                       help='Shift applied to input values (0.0: 0→0,1→1; -0.5: 0→-0.5,1→+0.5; -1.0: 0→-1,1→0)')
    parser.add_argument('--output', type=str, default='plot_all_degenerate_tetris.png',
                       help='Output filename')
    args = parser.parse_args()
    
    print(f"Creating solver...")
    solver = MQCASolver(preferred_vendor='nvidia', bPrint=False)
    
    print(f"\n=== Configuration ===")
    print(f"Cluster: {args.cluster}")
    print(f"W1: {args.W1}, W2: {args.W2}")
    print(f"Input positions: {args.input_pos}")
    print(f"Shift: {args.shift}")
    
    # Get cluster geometry
    make_cluster = CLUSTER_MAKERS[args.cluster]
    pos, inp_pos, inp_neigh, out_site = make_cluster(args.input_pos)
    nSite = len(pos)
    Eb = -np.ones(nSite, dtype=np.float32)
    
    print(f"Number of sites: {nSite}")
    print(f"Positions:\n{pos}")
    print(f"Input positions: {inp_pos}")
    print(f"Input neighbors: {inp_neigh}")
    print(f"Output site: {out_site}")
    
    # Build Esite and W
    Esite_batch = np.zeros((4, nSite), dtype=np.float32)
    for k, (A, B) in enumerate(INPUT_COMBOS):
        Esite_batch[k] = apply_input_bias(Eb, pos, inp_pos, inp_neigh, [A, B], args.W1, args.W2, args.shift)
    
    W_val, W_idx, nNeigh = sq_lattice_sparse(positions=pos, W1=args.W1, W2=args.W2, nSite=nSite)
    
    # Get top8 states
    W_val_batch = np.tile(W_val, (4, 1, 1))
    W_idx_batch = np.tile(W_idx, (4, 1, 1))
    nNeigh_batch = np.tile(nNeigh, (4, 1))
    
    E_top8, occ_top8 = solver.solve_batch_W_top8(Esite_batch, W_val_batch, W_idx_batch, nNeigh_batch, nSite)
    
    # Analyze degeneracy
    print(f"\n=== Degeneracy Analysis ===")
    for k, (A, B) in enumerate(INPUT_COMBOS):
        ground_gap = E_top8[k, 1] - E_top8[k, 0]
        status = "DEGENERATE" if ground_gap < 0.01 else "OK"
        print(f"Input ({A},{B}): gap={ground_gap:.4f} [{status}]")
    
    # Plot
    fig, axes = plt.subplots(2, 4, figsize=(20, 10))
    fig.suptitle(f'{args.cluster} at W1={args.W1}, W2={args.W2} (inputs {args.input_pos})', fontsize=14)
    
    for k, (A, B) in enumerate(INPUT_COMBOS):
        ground_gap = E_top8[k, 1] - E_top8[k, 0]
        if ground_gap < 0.01:
            n_states_to_plot = 2
            degenerate_str = " (DEGENERATE)"
        else:
            n_states_to_plot = 1
            degenerate_str = ""
        
        for i in range(n_states_to_plot):
            ax = axes[i, k]
            occ_mask = occ_top8[k, i]
            output_val = (occ_mask >> out_site) & 1
            ax.set_title(f'In({A},{B}) State {i+1}: E={E_top8[k,i]:.3f}, out={output_val}{degenerate_str}', fontsize=9)
            plot_tetris_state(ax, pos, occ_mask, out_site, inp_pos)
        
        if n_states_to_plot == 1:
            axes[1, k].axis('off')
    
    plt.tight_layout()
    plt.savefig(args.output, dpi=150)
    print(f"\nSaved: {args.output}")

if __name__ == '__main__':
    main()
