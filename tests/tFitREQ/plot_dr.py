#!/usr/bin/env python3
import argparse
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

try:
    import pacmap
    HAS_PACMAP = True
except ImportError:
    HAS_PACMAP = False
    
try:
    import umap
    HAS_UMAP = True
except ImportError:
    HAS_UMAP = False
    
from sklearn.decomposition import PCA

def main():
    parser = argparse.ArgumentParser(description="Dimensionality Reduction on Trajectories")
    parser.add_argument("--traj_file", type=str, default="path_trajectory.npz",  help="Path to trajectory/ensemble NPZ (chain_initial, chain_relaxed, energies)")
    parser.add_argument("--trajectory_npz", type=str, default=None, help="Alias for --traj_file (used by run_barrier_pipeline.sh)")
    parser.add_argument("--out_prefix", type=str, default="dr")
    args = parser.parse_args()

    if args.trajectory_npz:
        args.traj_file = args.trajectory_npz
    
    if not os.path.exists(args.traj_file):
        print(f"File {args.traj_file} not found. Run string_scan.py first.")
        sys.exit(1)
        
    data = np.load(args.traj_file)

    if 'chain_initial' in data and 'chain_relaxed' in data:
        # Path/trajectory format (string_scan)
        chain_initial = data['chain_initial']
        chain_relaxed = data['chain_relaxed']
        energies_init = data['energies_init']
        energies_relax = data['energies_relax']
        X = np.vstack([chain_initial, chain_relaxed])
        E = np.concatenate([energies_init, energies_relax]) if len(energies_relax) > 0 else energies_init
        n_init = len(chain_initial)
    elif 'dofs' in data and 'error' in data:
        # Ensemble format (grid_search)
        X = np.asarray(data['dofs'])
        E = np.asarray(data['error'])
        chain_initial = X
        chain_relaxed = np.zeros((0, X.shape[1]))
        n_init = len(chain_initial)
    else:
        print(f"Unsupported NPZ format: missing chain_* or dofs/error in {args.traj_file}")
        sys.exit(1)
    
    # Align energy array length with X for scatter coloring
    if len(E) != len(X):
        if len(E) < len(X):
            pad_val = E[-1] if len(E) > 0 else 0.0
            E = np.concatenate([E, np.full(len(X) - len(E), pad_val)])
        else:
            E = E[:len(X)]
    # Avoid log of zero or negative
    E_log = np.log(E - E.min() + 1e-6)
    
    # 1. PCA
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X)
    
    plt.figure(figsize=(8,6))
    sc = plt.scatter(X_pca[:, 0], X_pca[:, 1], c=E_log, cmap='viridis', s=50)
    plt.plot(X_pca[:n_init, 0], X_pca[:n_init, 1], 'k--', alpha=0.5, label='Initial Path')
    if len(chain_relaxed) > 0:
        plt.plot(X_pca[n_init:, 0], X_pca[n_init:, 1], 'b-', alpha=0.5, label='Relaxed Path')
    plt.colorbar(sc, label="Log(Error)")
    plt.title(f"PCA Projection (Explained Var: {pca.explained_variance_ratio_.sum():.2f})")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.legend()
    plt.savefig(f"{args.out_prefix}_pca.png", dpi=150)
    plt.close()
    
    # 2. UMAP
    if HAS_UMAP and len(X) >= 5:
        # For small datasets, UMAP needs adjusted n_neighbors
        n_neighbors = min(15, len(X) - 1)
        reducer = umap.UMAP(n_components=2, n_neighbors=n_neighbors, min_dist=0.1)
        X_umap = reducer.fit_transform(X)
        
        plt.figure(figsize=(8,6))
        sc = plt.scatter(X_umap[:, 0], X_umap[:, 1], c=E_log, cmap='viridis', s=50)
        plt.plot(X_umap[:n_init, 0], X_umap[:n_init, 1], 'k--', alpha=0.5, label='Initial Path')
        if len(chain_relaxed) > 0:
            plt.plot(X_umap[n_init:, 0], X_umap[n_init:, 1], 'b-', alpha=0.5, label='Relaxed Path')
        plt.colorbar(sc, label="Log(Error)")
        plt.title(f"UMAP Projection")
        plt.xlabel("UMAP 1")
        plt.ylabel("UMAP 2")
        plt.legend()
        plt.savefig(f"{args.out_prefix}_umap.png", dpi=150)
        plt.close()
        
    # 3. PacMAP
    if HAS_PACMAP and len(X) >= 5:
        n_neighbors = min(10, len(X) - 1)
        reducer = pacmap.PaCMAP(n_components=2, n_neighbors=n_neighbors, MN_ratio=0.5, FP_ratio=2.0)
        X_pacmap = reducer.fit_transform(X, init="pca")
        
        plt.figure(figsize=(8,6))
        sc = plt.scatter(X_pacmap[:, 0], X_pacmap[:, 1], c=E_log, cmap='viridis', s=50)
        plt.plot(X_pacmap[:n_init, 0], X_pacmap[:n_init, 1], 'k--', alpha=0.5, label='Initial Path')
        if len(chain_relaxed) > 0:
            plt.plot(X_pacmap[n_init:, 0], X_pacmap[n_init:, 1], 'b-', alpha=0.5, label='Relaxed Path')
        plt.colorbar(sc, label="Log(Error)")
        plt.title(f"PacMAP Projection")
        plt.xlabel("PacMAP 1")
        plt.ylabel("PacMAP 2")
        plt.legend()
        plt.savefig(f"{args.out_prefix}_pacmap.png", dpi=150)
        plt.close()

    print(f"Saved DR plots to {args.out_prefix}_*.png")

if __name__ == "__main__":
    main()
