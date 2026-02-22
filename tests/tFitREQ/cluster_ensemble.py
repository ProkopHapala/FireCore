#!/usr/bin/env python3
import argparse
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

try:
    import umap
    HAS_UMAP = True
except ImportError:
    HAS_UMAP = False

def main():
    parser = argparse.ArgumentParser(description="Cluster ensemble of DOFs using UMAP")
    parser.add_argument("--data", type=str, default="grid_search_out/ensemble_data.npz")
    parser.add_argument("--out_prefix", type=str, default="grid_search_out/ensemble_umap")
    args = parser.parse_args()

    if not os.path.exists(args.data):
        print(f"Data file {args.data} not found. Run grid_search.py first.")
        sys.exit(1)

    data = np.load(args.data)
    dofs = data['dofs']
    kMorse = data['kMorse']
    Lepairs = data['Lepairs']
    weight_alpha = data['weight_alpha']
    error = data['error']
    dof_names = data['dof_names']

    # Log error for better color scaling
    log_err = np.log10(error)

    if not HAS_UMAP:
        print("UMAP not installed. Falling back to PCA.")
        from sklearn.decomposition import PCA
        reducer = PCA(n_components=2)
        proj = reducer.fit_transform(dofs)
    else:
        # Use small n_neighbors since we only have 48 points
        reducer = umap.UMAP(n_components=2, n_neighbors=5, min_dist=0.1, random_state=42)
        proj = reducer.fit_transform(dofs)

    # Plot colored by different properties
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()

    props = [
        (log_err, "Log10(Error)", "viridis"),
        (kMorse, "kMorse", "plasma"),
        (Lepairs, "Lepairs", "coolwarm"),
        (weight_alpha, "Weight Alpha", "magma")
    ]

    for ax, (cdata, title, cmap) in zip(axes, props):
        sc = ax.scatter(proj[:, 0], proj[:, 1], c=cdata, cmap=cmap, s=80, alpha=0.8, edgecolor='k')
        ax.set_title(f"Colored by {title}")
        fig.colorbar(sc, ax=ax, label=title)
        
        # Annotate a few points (e.g. min/max error)
        if title == "Log10(Error)":
            min_idx = np.argmin(cdata)
            max_idx = np.argmax(cdata)
            ax.annotate(f"Min Err ({min_idx})", (proj[min_idx, 0], proj[min_idx, 1]), xytext=(5,5), textcoords='offset points')
            ax.annotate(f"Max Err ({max_idx})", (proj[max_idx, 0], proj[max_idx, 1]), xytext=(5,5), textcoords='offset points')

    plt.suptitle("Ensemble Projection (UMAP/PCA) of Optimized Parameters", fontsize=16)
    plt.tight_layout()
    plt.savefig(f"{args.out_prefix}_all.png", dpi=150)
    print(f"Saved {args.out_prefix}_all.png")

if __name__ == "__main__":
    main()
