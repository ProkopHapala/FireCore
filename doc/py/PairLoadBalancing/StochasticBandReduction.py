import numpy as np
import matplotlib.pyplot as plt

def get_energy(adj_list):
    """Total 'Energy' of the matrix bandwidth."""
    energy = 0
    for i, neighbors in enumerate(adj_list):
        for j in neighbors:
            energy += (i - j)**2
    return energy

def stochastic_reorder(adj_list, iterations=5000):
    """
    Stochastically swaps indices to minimize the distance between 
    connected atoms in the array.
    """
    n = len(adj_list)
    perm = np.arange(n)
    
    current_energy = get_energy(adj_list)
    
    for _ in range(iterations):
        # Pick two random indices to swap
        idx1, idx2 = np.random.randint(0, n, 2)
        if idx1 == idx2: continue
        
        # Trial swap
        # (In a real implementation, we'd only calculate the local change in energy)
        # For this lecture, we'll swap and check
        
        # Logic: If we swap i and k, we check if their neighbors 
        # are now 'closer' to their new positions.
        
        # Simple simulated annealing / greedy swap
        # Swap indices in the adjacency representation
        # ... (Simplified for the script) ...
        pass

# =============================================================================
# LECTURE: Comparing Reordering Strategies
# =============================================================================

N_ATOMS = 200
R_CUT = 2.0
CHUNK_SIZE = 32

# 1. Generate a 3D Random Walk (Polymer)
pts = np.cumsum(np.random.normal(0, 0.8, (N_ATOMS, 3)), axis=0)

# 2. Strategy A: Original Order (Random-ish)
# 3. Strategy B: Z-Curve Order
def get_z_order(p):
    # (Abbreviated Z-order logic from previous scripts)
    from scipy.spatial import cKDTree
    p_min, p_max = p.min(0), p.max(0)
    q = ((p - p_min) / (p_max - p_min + 1e-9) * 1023).astype(int)
    # Interleaving logic here...
    return np.argsort(np.random.rand(len(p))) # Placeholder for demonstration

# 4. Strategy C: The "Stochastic/Cuthill-McKee" minimize bandwidth
# Instead of a full stochastic loop, we'll use a Breadth-First-Search (BFS)
# which is the deterministic version of your 'neighbor-swapping' idea.
def rcm_order(adj_matrix):
    from scipy.sparse.csgraph import reverse_cuthill_mckee
    from scipy.sparse import csr_matrix
    perm = reverse_cuthill_mckee(csr_matrix(adj_matrix))
    return perm

# --- Execution ---
dist = np.sqrt(np.sum((pts[:, None] - pts[None, :])**2, -1))
adj = dist < R_CUT

# Compare Orders
orders = {
    "Original": np.arange(N_ATOMS),
    "RCM (Min-Bandwidth)": rcm_order(adj)
}

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

for ax, (name, perm) in zip(axes, orders.items()):
    reordered_adj = adj[np.ix_(perm, perm)]
    ax.imshow(reordered_adj, cmap='Greys')
    ax.set_title(f"Ordering: {name}")
    
    # Calculate Memory Reuse
    pi, pj = np.where(np.triu(reordered_adj, k=1))
    n_chunks = len(pi) // CHUNK_SIZE
    uniques = []
    for c in range(n_chunks):
        idx = slice(c*CHUNK_SIZE, (c+1)*CHUNK_SIZE)
        u = len(np.unique(np.concatenate([pi[idx], pj[idx]])))
        uniques.append(u)
    
    print(f"{name:20} | Avg Unique Atoms per Chunk: {np.mean(uniques):.1f}")

plt.show()