import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh

"""
Didactic Lecture: Load Balancing and Memory Reuse in GPU Molecular Dynamics

THE CHALLENGE:
GPU performance is limited by Memory Bandwidth. To maximize speed, we must load 
data into Shared Memory (LDS) and reuse it. For Density Matrix projection, we 
want "Dense Squares" in the adjacency matrix. 

THE Z-JUMP PROBLEM:
Z-curves map 3D -> 1D. However, Z-curves have "jumps" where two atoms adjacent 
in the 1D index are far apart in 3D. This creates "bi-nodal" clusters that 
waste LDS capacity.

THE SOLUTION:
We use PCA (Principal Component Analysis) to detect these stretched or jumped 
clusters and split them. We then re-home outliers and dissolve tiny clusters 
to ensure every GPU Work-group handles a dense, balanced workload.

ALGORITHM (didactic):
1) Z-order atoms to get a 1D sequence with spatial locality.
2) Chunk by target size; PCA-split clusters that are too stretched.
3) Prune atoms far from cluster COG (r_limit), dissolve tiny clusters, re-home ejected atoms.
4) Reorder clusters contiguously; compute adjacency.
5) Inter-cluster analysis:
   - If full block fill >= min_fill: mark green.
   - Else sample a pruned sub-block: keep atoms within rcut of midpoint between the two cluster COGs; if that sub-fill >= min_fill mark magenta.
Design goal: maximize dense blocks (better LDS/cache reuse) while keeping a fallback pruned view for borderline pairs.
"""

def get_z_indices(pts, res=1024):
    """Computes Morton (Z-curve) codes for 3D points."""
    def part_1d(n):
        n = int(n) & 0x3ff
        n = (n | (n << 16)) & 0xff0000ff
        n = (n | (n << 8))  & 0x0300f00f
        n = (n | (n << 4))  & 0x030c30c3
        n = (n | (n << 2))  & 0x09249249
        return n
    p_min, p_max = pts.min(0), pts.max(0)
    q = ((pts - p_min) / (p_max - p_min + 1e-9) * (res - 1)).astype(int)
    return np.array([(part_1d(p[0])<<2 | part_1d(p[1])<<1 | part_1d(p[2])) for p in q])

def pca_split_check(cluster_pts, spread_threshold):
    """
    Calculates the second moment tensor. 
    Returns (should_split, split_mask).
    """
    if len(cluster_pts) < 4:
        return False, None
    
    centered = cluster_pts - cluster_pts.mean(axis=0)
    cov = np.dot(centered.T, centered) / (len(cluster_pts) - 1)
    
    # Trace is the sum of eigenvalues (total variance)
    eigenvalues, eigenvectors = eigh(cov)
    trace = np.sum(eigenvalues)
    
    if trace > spread_threshold:
        # Split along the principal eigenvector (longest axis)
        v_main = eigenvectors[:, -1]
        projection = np.dot(centered, v_main)
        return True, projection > 0
    return False, None

def partition_system(pts, target_size=32, spread_limit=2.0, r_limit=2.5, min_size=8, rcut=3.5, min_fill=0.05):
    n = len(pts)
    # 1. Initial Z-Sort
    z_idx = np.argsort(get_z_indices(pts))
    pts_z = pts[z_idx]
    
    # 2. Initial Chunking
    n_initial = n // target_size
    labels = -1 * np.ones(n, dtype=int)
    cluster_list = []
    
    for i in range(n_initial):
        cluster_list.append(list(range(i*target_size, (i+1)*target_size)))
        
    # 3. PCA Splitting Pass
    refined_clusters = []
    for cluster in cluster_list:
        should_split, mask = pca_split_check(pts_z[cluster], spread_limit)
        if should_split:
            cluster = np.array(cluster)
            refined_clusters.append(cluster[mask].tolist())
            refined_clusters.append(cluster[~mask].tolist())
        else:
            refined_clusters.append(cluster)
            
    # Update labels based on refined clusters
    labels[:] = -1
    for idx, cluster in enumerate(refined_clusters):
        labels[cluster] = idx
        
    # 4. Pruning and Re-homing
    cogs = np.array([pts_z[labels == i].mean(0) if np.any(labels==i) else [0,0,0] for i in range(len(refined_clusters))])
    
    for i in range(n):
        l = labels[i]
        if l != -1:
            if np.linalg.norm(pts_z[i] - cogs[l]) > r_limit:
                labels[i] = -1 # Eject
                
    # 5. Dissolve small clusters (before re-homing so dissolved atoms can be re-homed)
    final_cluster_ids = np.unique(labels)
    for cid in final_cluster_ids:
        if cid == -1: continue
        members = np.where(labels == cid)[0]
        if len(members) < min_size:
            labels[members] = -1

    # Recompute COGs after dissolving
    cogs = np.array([pts_z[labels == i].mean(0) if np.any(labels==i) else [0,0,0] for i in range(len(refined_clusters))])

    # Re-home outsiders
    for i in np.where(labels == -1)[0]:
        dists = np.linalg.norm(pts_z[i] - cogs, axis=1)
        best = np.argmin(dists)
        if dists[best] < r_limit:
            labels[i] = best

    # 6. Cluster-pair dense check (between-cluster blocks)
    def cluster_pair_stats(pts_z, labels, rcut, r_limit, min_fill):
        """Identify inter-cluster pairs that are still dense enough to be considered together."""
        valid_labels = [l for l in np.unique(labels) if l != -1]
        if len(valid_labels) < 2:
            return []
        dist = np.sqrt(np.sum((pts_z[:, None, :] - pts_z[None, :, :])**2, -1))
        adj = dist < rcut
        dense_pairs = []
        for i_idx, li in enumerate(valid_labels):
            mask_i = labels == li
            if not np.any(mask_i): continue
            cog_i = pts_z[mask_i].mean(axis=0)
            for lj in valid_labels[i_idx+1:]:
                mask_j = labels == lj
                if not np.any(mask_j): continue
                cog_j = pts_z[mask_j].mean(axis=0)
                cog_dist = np.linalg.norm(cog_i - cog_j)
                # Skip if clusters are too far to interact
                if cog_dist > (2*r_limit + rcut):
                    continue
                contacts = adj[np.ix_(mask_i, mask_j)].sum()
                denom = mask_i.sum() * mask_j.sum()
                fill = contacts / denom if denom > 0 else 0.0
                if fill >= min_fill:
                    dense_pairs.append((li, lj, fill, contacts, mask_i.sum(), mask_j.sum()))
        return dense_pairs

    # Re-index to enforce contiguous ranges
    sort_key = labels.astype(float)
    sort_key[sort_key == -1] = np.inf # Sparse at the end
    final_order = np.lexsort((np.arange(n), sort_key))
    pts_final = pts_z[final_order]
    labels_final = labels[final_order]

    # Recompute dense pairs on the reordered system
    dense_pairs = cluster_pair_stats(pts_final, labels_final, rcut, r_limit, min_fill=min_fill)
    
    return pts_final, labels_final, dense_pairs

def plot_results(pts, labels, rcut, dense_pairs=None, min_fill=0.3):
    n = len(pts)
    dist = np.sqrt(np.sum((pts[:, None] - pts[None, :])**2, -1))
    adj = dist < rcut
    total_pairs = np.triu(adj, k=1).sum()
    adj_display = adj.astype(float)*0.6 + 0.2 # lighter gray fills for visibility
    
    fig, axes = plt.subplots(1, 2, figsize=(15, 7))
    
    # Adjacency Matrix (gray fills for contrast with boxes)
    axes[0].imshow(adj_display, cmap='Greys', interpolation='nearest', origin='lower', vmin=0.0, vmax=1.0)
    u_labels = np.unique(labels)
    cluster_ranges = {}
    sparse_start = n
    for l in u_labels:
        if l == -1: continue
        idx = np.where(labels == l)[0]
        start, end = idx[0], idx[-1]
        rect = plt.Rectangle((start, start), end-start, end-start, 
                             edgecolor='red', facecolor='none', lw=1, alpha=0.6)
        axes[0].add_patch(rect)
        sparse_start = min(sparse_start, start) if start > 0 else sparse_start
        cluster_ranges[l] = (start, end)
        
    axes[0].set_title("Optimized Partitioning\nRed Boxes = GPU Dense-Tile Tasks")
    pairs_drawn = 0
    dense_interactions = 0
    block_records = []
    valid_labels = [l for l in u_labels if l != -1]
    # Record red (diagonal) blocks
    for li in valid_labels:
        mask_i = labels == li
        sub_adj = adj[np.ix_(mask_i, mask_i)]
        sub_pairs = np.triu(sub_adj, k=1).sum()
        ni = mask_i.sum()
        denom_red = ni*(ni-1)/2 if ni>1 else 1
        fill_red = sub_pairs / denom_red
        dense_interactions += sub_pairs
        block_records.append(("red", li, li, fill_red, ni, ni))
    for idx_i, li in enumerate(valid_labels):
        mask_i = labels == li
        if not np.any(mask_i): continue
        for lj in valid_labels[idx_i+1:]:
            mask_j = labels == lj
            if not np.any(mask_j): continue
            i0, i1 = cluster_ranges[li]
            j0, j1 = cluster_ranges[lj]
            contacts = adj[np.ix_(mask_i, mask_j)].sum()
            denom = mask_i.sum() * mask_j.sum()
            fill = contacts / denom if denom > 0 else 0.0
            if fill >= min_fill:
                axes[0].add_patch( plt.Rectangle((j0, i0), (j1 - j0), (i1 - i0), edgecolor='green', facecolor='none', lw=1, alpha=0.6) )
                axes[0].add_patch( plt.Rectangle((i0, j0), (i1 - i0), (j1 - j0), edgecolor='green', facecolor='none', lw=1, alpha=0.6) )
                axes[0].text(j0, i0, f"{fill:.2f}", color='green', fontsize=7)
                pairs_drawn += 1
                dense_interactions += contacts
                block_records.append(("green", li, lj, fill, mask_i.sum(), mask_j.sum()))
            else:
                # magenta: pruned subset around midpoint of COGs
                idx_i = np.where(mask_i)[0]; idx_j = np.where(mask_j)[0]
                cog_i = pts[idx_i].mean(axis=0); cog_j = pts[idx_j].mean(axis=0)
                mid = 0.5*(cog_i + cog_j)
                sel_i = idx_i[np.linalg.norm(pts[idx_i]-mid, axis=1) <= rcut]
                sel_j = idx_j[np.linalg.norm(pts[idx_j]-mid, axis=1) <= rcut]
                if len(sel_i)==0 or len(sel_j)==0:
                    continue
                sub_contacts = adj[np.ix_(sel_i, sel_j)].sum()
                sub_denom = len(sel_i)*len(sel_j)
                sub_fill = sub_contacts / sub_denom if sub_denom>0 else 0.0
                if sub_fill >= min_fill:
                    si0, si1 = sel_i.min(), sel_i.max()
                    sj0, sj1 = sel_j.min(), sel_j.max()
                    axes[0].add_patch( plt.Rectangle((sj0, si0), (sj1 - sj0), (si1 - si0), edgecolor='magenta', facecolor='none', lw=1, alpha=0.6) )
                    axes[0].add_patch( plt.Rectangle((si0, sj0), (si1 - si0), (sj1 - sj0), edgecolor='magenta', facecolor='none', lw=1, alpha=0.6) )
                    axes[0].text(sj0, si0, f"{sub_fill:.2f}", color='magenta', fontsize=7)
                    pairs_drawn += 1
                    dense_interactions += sub_contacts
                    block_records.append(("magenta", li, lj, sub_fill, len(sel_i), len(sel_j)))
    if pairs_drawn:
        axes[0].set_xlabel(f"Dense/sampled inter-cluster pairs: {pairs_drawn}")
    
    # 3D Spatial View
    ax3d = fig.add_subplot(122, projection='3d')
    for l in u_labels:
        mask = labels == l
        color = 'gray' if l == -1 else None
        ax3d.scatter(pts[mask,0], pts[mask,1], pts[mask,2], s=15, alpha=0.6, label='Sparse' if l==-1 else None)
    ax3d.set_title("3D Clusters (PCA Refined)")
    
    plt.tight_layout()
    plt.show()

    # Print block summary with dims
    order_map = {"red":0, "green":1, "magenta":2}
    block_records_sorted = sorted(block_records, key=lambda x: (order_map.get(x[0], 3), x[3]), reverse=True)
    if block_records_sorted:
        print("Blocks:")
        for typ, li, lj, fill, si, sj in block_records_sorted:
            print(f"  {typ:7} ({li},{lj}) fill={fill:.3f} size=({si},{sj})")
    sparse_pairs = total_pairs - dense_interactions
    print(f"Dense interactions: {dense_interactions} | Sparse kernel pairs: {sparse_pairs}")

    # Histograms
    if block_records_sorted:
        fills = [r[3] for r in block_records_sorted]
        sizes = [r[4] for r in block_records_sorted] + [r[5] for r in block_records_sorted]
        fig2, axes2 = plt.subplots(1,2, figsize=(10,4))
        axes2[0].hist(fills, bins=20, color='gray', edgecolor='black')
        axes2[0].set_title("Fill factors")
        axes2[1].hist(sizes, bins=20, color='gray', edgecolor='black')
        axes2[1].set_title("Block side lengths")
        plt.tight_layout()
        plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser( description="PCA-Refined Molecular Partitioning for GPU")
    parser.add_argument("--n_atoms",  type=int, default=800)
    parser.add_argument("--target",   type=int, default=32, help="Target atoms per GPU work-group")
    parser.add_argument("--rcut",     type=float, default=3.5)
    parser.add_argument("--spread",   type=float, default=5.0, help="PCA Trace threshold for splitting")
    parser.add_argument("--r_limit",  type=float, default=3.0, help="Max distance from COG before ejection")
    parser.add_argument("--min_size", type=int,   default=10, help="Dissolve clusters smaller than this")
    parser.add_argument("--min_fill", type=float, default=0.5, help="Minimum fill fraction for dense pairs")
    parser.add_argument("--seed",     type=int, default=46)
    parser.add_argument("--walk_sigma", type=float, default=0.5)
    args = parser.parse_args()

    # Generate a polymer-like string to stress the Z-jump continuity
    pts = np.cumsum(np.random.normal(0, args.walk_sigma, (args.n_atoms, 3)), axis=0)
    
    pts_opt, labels, dense_pairs = partition_system(pts, args.target, args.spread, args.r_limit, args.min_size, args.rcut, args.min_fill)
    
    # Compute Metrics
    n_in_clusters = np.sum(labels != -1)
    print(f"Total Atoms: {args.n_atoms}")
    print(f"Atoms in Dense Clusters: {n_in_clusters} ({n_in_clusters/args.n_atoms:.1%})")
    print(f"Sparse Residual Atoms:  {args.n_atoms - n_in_clusters}")
    if dense_pairs:
        print("Dense inter-cluster blocks (li, lj, fill, contacts, size_i, size_j):")
        for li, lj, fill, contacts, si, sj in dense_pairs:
            print(f"  ({li},{lj}) fill={fill:.3f} contacts={contacts} sizes={si}/{sj}")
    
    plot_results(pts_opt, labels, args.rcut, dense_pairs=dense_pairs, min_fill=args.min_fill)