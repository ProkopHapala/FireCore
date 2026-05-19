import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

# =====================================================================
# 1. GEOMETRY: Build a 2D Hexagonal Flake (e.g., Polycyclic Aromatic)
# =====================================================================
def generate_graphene_flake(radius):
    """Generates coordinates for a honeycomb lattice within a given radius."""
    a = 1.42 # C-C bond length in Angstroms
    pts = []
    for i in range(-5, 6):
        for j in range(-5, 6):
            # Hexagonal lattice basis
            x_base = 1.5 * a * i
            y_base = np.sqrt(3) * a * (j + 0.5 * (i % 2))
            pts.append([x_base, y_base])           # Atom A
            pts.append([x_base + a, y_base])       # Atom B
            
    pts = np.array(pts)
    # Keep atoms inside the radius to form a finite flake
    distances = np.linalg.norm(pts, axis=1)
    pts_filtered = pts[distances < radius]
    # Adding a tiny random noise to coordinates to make exponential hopping dynamic
    noise = np.random.normal(0, 0.01, pts_filtered.shape) 
    return pts_filtered + noise

coords = generate_graphene_flake(radius=3.8) # Generates a flake of ~24-30 atoms
N = len(coords)

# =====================================================================
# 2. HAMILTONIAN: Exponential Hopping & Nitrogen Defect
# =====================================================================
t0 = -2.7        # Base hopping (eV)
r0 = 1.42        # Equilibrium distance (Angstrom)
beta = 3.0       # Exponential decay factor (1/Angstrom)
e_nitrogen = -1.5 # On-site energy for Nitrogen (more electronegative)

H = np.zeros((N, N))
bonds = []

# Populate hopping matrix
for i in range(N):
    for j in range(i + 1, N):
        dist = np.linalg.norm(coords[i] - coords[j])
        if dist < 1.8: # Cutoff for nearest neighbors
            t = t0 * np.exp(-beta * (dist - r0))
            H[i, j] = H[j, i] = t
            bonds.append((i, j))

# Introduce a Nitrogen Defect near the center
distances_to_center = np.linalg.norm(coords, axis=1)
defect_site = np.argmin(distances_to_center)
H[defect_site, defect_site] = e_nitrogen

# =====================================================================
# 3. EXACT SOLUTION: O(N^3) Diagonalization
# =====================================================================
evals, evecs = np.linalg.eigh(H)
rho_exact = np.zeros((N, N))
for lam in range(N):
    if evals[lam] < 0.0: # Half-filled (mu = 0)
        rho_exact += 2.0 * np.outer(evecs[:, lam], evecs[:, lam]) # x2 for spin

# =====================================================================
# 4. LINEAR SCALING BOP: Chebyshev Local Probe
# =====================================================================
def get_chebyshev_coeffs(M):
    """Jackson-damped Chebyshev coefficients for Step Function (E=0)"""
    c = np.zeros(M)
    for k in range(1, M, 2): 
        c_raw = -(2.0 / (np.pi * k)) * np.sin(k * np.pi / 2.0)
        ang = np.pi / (M + 1)
        damping = ((M - k + 1) * np.cos(k * ang) + np.sin(k * ang) / np.tan(ang)) / (M + 1)
        c[k] = c_raw * damping
    return c

# Scale H to [-1, 1] range. Spectral bounds for Graphene are roughly +/- 3*|t0|
E_max = 3.0 * np.abs(t0) + np.abs(e_nitrogen) + 0.5
H_scaled = H / E_max

def get_bond_density_convergence(H_scaled, start_site, target_site, M_max):
    """Drops probe at start_site, records density accumulation up to M_max hops."""
    c = get_chebyshev_coeffs(M_max)
    v_prev = np.zeros(N)
    v_curr = np.zeros(N)
    v_curr[start_site] = 1.0
    
    accumulated_density = np.zeros(M_max)
    current_density = 0.0
    
    # Hop 1
    v_next = H_scaled @ v_curr
    current_density += c[1] * v_next[target_site] * 2.0 # x2 for spin
    accumulated_density[1] = current_density
    
    # Hop 2 to M_max
    for k in range(2, M_max):
        v_next_next = 2.0 * (H_scaled @ v_next) - v_curr
        
        if k % 2 != 0: # Only odd paths contribute to bipartite bond density
            current_density += c[k] * v_next_next[target_site] * 2.0
            
        accumulated_density[k] = current_density
        
        v_prev, v_curr, v_next = v_curr, v_next, v_next_next
        
    return accumulated_density

# Run BOP for all bonds at M=15 to get the final map
M_target = 15
rho_bop = np.zeros(len(bonds))
for idx, (i, j) in enumerate(bonds):
    # We only need the final value at M_target
    conv = get_bond_density_convergence(H_scaled, i, j, M_target)
    rho_bop[idx] = conv[-1]

# =====================================================================
# 5. SELECT SPECIFIC BONDS FOR CONVERGENCE ANALYSIS
# =====================================================================
# We will track 3 different physical environments
edge_bond = None
core_bond = None
defect_bond = None

for idx, (i, j) in enumerate(bonds):
    if i == defect_site or j == defect_site:
        defect_bond = (idx, i, j)
    elif np.linalg.norm((coords[i] + coords[j])/2) > 2.5:
        edge_bond = (idx, i, j)
    elif np.linalg.norm((coords[i] + coords[j])/2) < 1.5:
        core_bond = (idx, i, j)

tracked_bonds = {
    "Edge Bond (Strong Kekulé)": edge_bond,
    "Core C-C Bond (Delocalized)": core_bond,
    "C-N Defect Bond": defect_bond
}

# =====================================================================
# 6. VISUALIZATION
# =====================================================================
fig = plt.figure(figsize=(15, 6))

# --- PLOT 1: The Molecular Graph (Exact vs BOP) ---
ax1 = plt.subplot(1, 2, 1)

# Prepare line segments and widths
segments = [(coords[i], coords[j]) for (i, j) in bonds]
# Density usually ranges from 0.4 (single) to ~0.9 (double)
# We normalize line width to make the Kekule structure obvious
widths_exact = (np.array([rho_exact[i, j] for i, j in bonds]) - 0.3) * 12
widths_bop = (rho_bop - 0.3) * 12

lc_bop = LineCollection(segments, linewidths=widths_bop, colors='blue', alpha=0.5, label='BOP Density')
lc_exact = LineCollection(segments, linewidths=widths_exact, colors='black', linestyle=':', alpha=0.8, label='Exact Density')

ax1.add_collection(lc_bop)
ax1.add_collection(lc_exact)

# Plot Atoms
ax1.scatter(coords[:, 0], coords[:, 1], c='black', s=50, zorder=3)
# Highlight Nitrogen Defect
ax1.scatter(coords[defect_site, 0], coords[defect_site, 1], c='red', s=120, zorder=4, label='Nitrogen Defect')

# Highlight Tracked Bonds
colors = ['green', 'purple', 'orange']
tracked_bonds_filtered = {k: v for k, v in tracked_bonds.items() if v is not None}
for c_idx, (name, (idx, i, j)) in enumerate(tracked_bonds_filtered.items()):
    ax1.plot([coords[i, 0], coords[j, 0]], [coords[i, 1], coords[j, 1]], 
             color=colors[c_idx], linewidth=5, zorder=2)

ax1.set_aspect('equal')
ax1.set_title(f"Aromatic Flake Bond Densities ($N={N}$)\nLine Thickness $\propto \\rho_{{ij}}$", fontsize=14)
ax1.axis('off')
ax1.legend(loc='lower left')

# --- PLOT 2: Convergence vs Iterations ---
ax2 = plt.subplot(1, 2, 2)
M_max = 20
M_range = np.arange(1, M_max, 2) # Only odd hops

for c_idx, (name, (idx, i, j)) in enumerate(tracked_bonds_filtered.items()):
    conv = get_bond_density_convergence(H_scaled, i, j, M_max)
    exact_val = rho_exact[i, j]
    
    # Plot BOP progression
    ax2.plot(M_range, conv[M_range], marker='o', color=colors[c_idx], linewidth=2, label=f'BOP: {name}')
    # Plot Exact horizontal line
    ax2.axhline(exact_val, color=colors[c_idx], linestyle='--', alpha=0.6)

ax2.set_title("Bond Density Convergence vs Max Scattering Path ($M$)", fontsize=14)
ax2.set_xlabel("Max Propagation Hops (Path Length $M$)", fontsize=12)
ax2.set_ylabel("Estimated Bond Density $\\rho_{ij}$", fontsize=12)
ax2.set_xticks(M_range)
ax2.grid(True, alpha=0.3)
ax2.legend(loc='lower right')

plt.tight_layout()
plt.show()