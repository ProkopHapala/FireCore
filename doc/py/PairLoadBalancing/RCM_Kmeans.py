import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import reverse_cuthill_mckee

# =============================================================================
# LECTURE: RCM vs. CLUSTERING
# =============================================================================

def get_z_indices(pts, res=1024):
    """ Standard Z-curve implementation for initial spatial grouping. """
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

def cluster_atoms_refined(pts, target_size=16, penalty_weight=2.0):
    """
    Refined Clustering:
    1. Initial Guess: Use Z-curve spans (inherently spatially local).
    2. Refinement: A K-means-like pass that penalizes cluster growth.
    """
    n = len(pts)
    z_idx = np.argsort(get_z_indices(pts))
    n_clusters = n // target_size
    
    # Initialize clusters based on Z-order sequence
    # (This is MUCH faster than random initialization)
    labels = np.zeros(n, dtype=int)
    for i in range(n_clusters):
        labels[z_idx[i*target_size : (i+1)*target_size]] = i
        
    # Refinement iterations
    for _ in range(2): # Only 1-2 passes needed because Z-init is so good
        pivots = np.array([pts[labels == i].mean(0) for i in range(n_clusters)])
        counts = np.bincount(labels, minlength=n_clusters)
        
        for i in range(n):
            dists = np.sum((pts[i] - pivots)**2, axis=1)
            # COST = Distance + Penalty for being a large cluster
            cost = dists + penalty_weight * (counts / target_size)
            new_l = np.argmin(cost)
            
            if new_l != labels[i]:
                counts[labels[i]] -= 1
                labels[i] = new_l
                counts[new_l] += 1
                
    return labels

# --- Setup System ---
N_ATOMS = 512
R_CUT = 2.0
CHUNK_SIZE = 32
pts = np.cumsum(np.random.normal(0, 0.7, (N_ATOMS, 3)), axis=0) # Polymer

# Create Adjacency
dist = np.sqrt(np.sum((pts[:, None] - pts[None, :])**2, -1))
adj = dist < R_CUT

# 1. RCM Ordering
rcm_idx = reverse_cuthill_mckee(csr_matrix(adj))

# 2. Cluster Ordering
cluster_labels = cluster_atoms_refined(pts, target_size=16)
cluster_idx = np.argsort(cluster_labels)

# --- Visualization ---
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

# A: Original Matrix
axes[0].imshow(adj, cmap='Greys')
axes[0].set_title("1. Original Order\n(Scattered Interactions)")

# B: RCM Matrix
adj_rcm = adj[np.ix_(rcm_idx, rcm_idx)]
axes[1].imshow(adj_rcm, cmap='Greys')
axes[1].set_title("2. RCM Ordering\n(The 'Snake' - Minimize Bandwidth)")

# C: Clustered Matrix
adj_clust = adj[np.ix_(cluster_idx, cluster_idx)]
axes[2].imshow(adj_clust, cmap='Greys')
# Draw grid lines for the dense blocks
for i in range(0, N_ATOMS, 16):
    axes[2].axhline(i, color='red', alpha=0.2, lw=0.5)
    axes[2].axvline(i, color='red', alpha=0.2, lw=0.5)
axes[2].set_title("3. Spatial Clustering\n(The 'Blocks' - Dense Squares)")

plt.tight_layout()
plt.show()

# --- Analysis of LDS Reuse ---
def calc_reuse(matrix, chunk):
    pi, pj = np.where(np.triu(matrix, k=1))
    uniques = [len(np.unique(np.concatenate([pi[i:i+chunk], pj[i:i+chunk]]))) 
               for i in range(0, len(pi)-chunk, chunk)]
    return np.mean(uniques)

print(f"Avg Unique Atoms (LDS load) for {CHUNK_SIZE} pairs:")
print(f"RCM Order:     {calc_reuse(adj_rcm, CHUNK_SIZE):.1f}")
print(f"Cluster Order: {calc_reuse(adj_clust, CHUNK_SIZE):.1f}")