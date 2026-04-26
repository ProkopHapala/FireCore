"""
Topological grain boundary in Kekulé-distorted graphene
=======================================================
This script illustrates the formation of a conductive channel at a domain wall
between two different Kekulé bond patterns (aromatic vs quinoid) in graphene.

Physics background:
- Graphene π-electrons described by nearest-neighbour tight-binding.
- Kekulé bond alternation couples the two Dirac valleys and opens a gap.
- The effective low-energy Hamiltonian is a Dirac equation with a mass term m.
- The phase of the Kekulé pattern determines the sign/phase of m.
- A domain wall where m changes sign hosts a topologically protected zero-energy
  soliton state (Jackiw-Rebbi bound state) → 1D conducting channel.
- We compute the local density of states (LDOS) at E=0 to visualise the channel.

Author: illustrative script for learning purposes
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh

# ==============================================================================
# 1. Lattice generation: finite rectangular graphene flake
# ==============================================================================

# Carbon-carbon distance (in Angstrom, but we work in units where a_CC = 1)
a_CC = 1.0

# Lattice vectors of the honeycomb (A sublattice coordinates)
# Using orthogonal unit cell: a1 = (sqrt(3), 0), a2 = (sqrt(3)/2, 3/2)
# but we generate atoms by placing A and B sublattices directly.
# We'll build a flake that extends in x (width) and y (height).
# Domain wall is at x = 0.0, running along y-direction.
# We choose flake size: x from -Lx to Lx, y from -Ly to Ly.
Lx = 5.0   # half-width in units of a_CC
Ly = 3.0
dx = np.sqrt(3) * a_CC      # distance between neighbouring A atoms along x
dy = 1.0 * a_CC             # sub-unit in y (we'll use 3*dy for A spacing)
# Actually let's generate using integer cell indices.

# Build rectangular lattice of unit cells.
# Unit cell contains one A and one B atom.
# Lattice vectors: a1 = (sqrt(3)*a_CC, 0)
#                  a2 = (sqrt(3)/2 * a_CC, 1.5 * a_CC)
a1 = np.array([np.sqrt(3)*a_CC, 0.0])
a2 = np.array([np.sqrt(3)/2 * a_CC, 1.5 * a_CC])
# Basis: A at (0,0), B at (a1+a2)/3 ? Actually standard: B at (0, a_CC) in orthogonal?
# For graphene, A at (0,0), B at (a_CC, 0) if using the orthogonal cell?
# Let's use a simpler direct placement:
# Place A atoms at positions (sqrt(3)*i + sqrt(3)/2*j, 1.5*j) for integer i,j.
# Then B atoms at A + (0, a_CC). This yields honeycomb with vertical bonds.
# This orientation gives armchair edges in y-direction? We'll work with any.
# We'll define a function to generate atom positions.

def generate_flake(Lx, Ly):
    """
    Generate atom positions and sublattice indices for a graphene flake.
    Flake is a rectangle in Cartesian coordinates: |x| <= Lx, |y| <= Ly.
    Returns positions (Nx2) and sublattice (0 for A, 1 for B).
    """
    # Grid of A atom positions
    a1 = np.array([np.sqrt(3), 0.0])
    a2 = np.array([np.sqrt(3)/2, 1.5])
    # We'll iterate over integers n1, n2 covering the desired range.
    # Find min/max n1,n2 such that positions are within x,y bounds.
    n1_min = int(np.floor(-Lx / a1[0])) - 1
    n1_max = int(np.ceil( Lx / a1[0])) + 1
    n2_min = int(np.floor(-Ly / a2[1])) - 1
    n2_max = int(np.ceil( Ly / a2[1])) + 1
    
    pos_A = []
    for n1 in range(n1_min, n1_max+1):
        for n2 in range(n2_min, n2_max+1):
            r = n1*a1 + n2*a2
            if abs(r[0]) <= Lx and abs(r[1]) <= Ly:
                pos_A.append(r)
    pos_A = np.array(pos_A)
    # B atoms: shift by d = (0, a_CC)
    pos_B = pos_A + np.array([0.0, a_CC])
    # Combine
    positions = np.vstack([pos_A, pos_B])
    sublattice = np.array([0]*len(pos_A) + [1]*len(pos_B))
    return positions, sublattice

positions, sublattice = generate_flake(Lx, Ly)
N_atoms = len(positions)
print(f"Number of atoms: {N_atoms}")

# ==============================================================================
# 2. Define the domain wall and Kekulé hopping modulation
# ==============================================================================

# Tight-binding parameters
t0 = 2.8   # average hopping (eV)
Delta = 0.4  # bond alternation strength (modulation amplitude), must be < t0

# K-point vector for the Kekulé modulation.
# For graphene, the Dirac point K = (4π/(3√3 a_CC), 0) in Cartesian coordinates?
# Actually using a_CC=1, K = (4π/(3√3), 0) ≈ (2.418, 0).
# The modulation uses the dot product K·δ_j where δ_j are the three bond vectors.
K_vec = np.array([4*np.pi/(3*np.sqrt(3)), 0.0])

# Bond vectors from an A atom to its three B neighbours.
# With our geometry: B is at (0, a_CC) relative to A,
# plus shifts by lattice vectors.
# The three nearest-neighbour vectors:
delta1 = np.array([0.0, a_CC])           # straight up
delta2 = np.array([-np.sqrt(3)/2*a_CC, -0.5*a_CC])  # bottom-left
delta3 = np.array([ np.sqrt(3)/2*a_CC, -0.5*a_CC])  # bottom-right
bond_vectors = [delta1, delta2, delta3]

# Precompute K·δ for each bond direction
Kdot_delta = [np.dot(K_vec, d) for d in bond_vectors]

# Domain wall position (vertical line at x=0)
x_wall = 0.0

def phi_domain(r):
    """Domain phase: 0 for x < 0, pi for x > 0."""
    if r[0] < x_wall:
        return 0.0
    else:
        return np.pi

# We need to build a list of bonds (i, j) with hoppings.
# Bonds exist between A and B atoms that are nearest neighbours.
# Use a distance criterion: |r_i - r_j| ≈ a_CC (within tolerance).

tol = 0.1 * a_CC
bonds = []   # (i, j, t_ij)
for i in range(N_atoms):
    if sublattice[i] != 0:   # only A atoms (sublattice 0) as starting point to avoid duplicates
        continue
    pos_i = positions[i]
    for j in range(N_atoms):
        if sublattice[j] != 1:  # B atoms only
            continue
        pos_j = positions[j]
        vec = pos_j - pos_i
        dist = np.linalg.norm(vec)
        if abs(dist - a_CC) < tol:
            # Determine which bond vector this corresponds to
            # Compare angle/fidelity with the three bond_vectors
            # Use dot product to identify the closest direction.
            best_k = -1
            best_dot = -2.0
            for k, d in enumerate(bond_vectors):
                dot = np.dot(vec, d) / (dist * a_CC)
                if dot > best_dot:
                    best_dot = dot
                    best_k = k
            # If best_dot > 0.9 we trust it
            if best_dot > 0.9:
                # Midpoint of bond
                r_mid = 0.5*(pos_i + pos_j)
                phi = phi_domain(r_mid)
                # Hopping modulation
                t = t0 + Delta * np.cos(Kdot_delta[best_k] + phi)
                bonds.append((i, j, t))

print(f"Number of bonds: {len(bonds)}")

# ==============================================================================
# 3. Build the tight-binding Hamiltonian
# ==============================================================================
H = np.zeros((N_atoms, N_atoms), dtype=float)
for i, j, t in bonds:
    H[i, j] = t
    H[j, i] = t  # Hermitian

# ==============================================================================
# 4. Diagonalise Hamiltonian
# ==============================================================================
# Use dense eigh for N_atoms ~ a few hundred (fast enough)
E, psi = eigh(H)

# Find Fermi level (assume half-filling: number of electrons = N_atoms)
# In neutral graphene, each atom contributes one pz electron, so N_e = N_atoms.
# At half-filling, Fermi energy = (E[N/2-1] + E[N/2])/2 (if N even).
N_e = N_atoms
if N_e % 2 == 0:
    E_F = 0.5 * (E[N_e//2 - 1] + E[N_e//2])
else:
    E_F = E[(N_e-1)//2]
print(f"Fermi energy (approx): {E_F:.3f} eV")

# ==============================================================================
# 5. Local Density of States at E=0
# ==============================================================================
# Broadening parameter (eV)
sigma = 0.05
# Gaussian weight for each state
weights = np.exp(- (E - 0.0)**2 / (2*sigma**2)) / (np.sqrt(2*np.pi)*sigma)
# LDOS on each atom
ldos = np.zeros(N_atoms)
for n in range(N_atoms):
    # contribution of state n to LDOS at r: |psi_n(r)|^2 * weight
    ldos += np.abs(psi[:, n])**2 * weights[n]

# Normalise LDOS to max=1 for visualisation
ldos_norm = ldos / ldos.max()

# ==============================================================================
# 6. Density of States (DOS)
# ==============================================================================
E_range = np.linspace(-3.0, 3.0, 500)
dE = E_range[1] - E_range[0]
DOS = np.zeros_like(E_range)
sigma_dos = 0.02  # smaller broadening for smooth curve
for e_n in E:
    DOS += np.exp(- (E_range - e_n)**2 / (2*sigma_dos**2)) / (np.sqrt(2*np.pi)*sigma_dos)
DOS /= N_atoms  # per atom

# ==============================================================================
# 7. Visualisation
# ==============================================================================
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# --- (a) Lattice with LDOS as circles ---
ax = axes[0, 0]
ax.set_aspect('equal')
# Draw domain wall
ax.axvline(x=x_wall, color='red', linestyle='--', linewidth=2, alpha=0.7, label='Domain wall')
# Plot atoms with size proportional to LDOS
for i in range(N_atoms):
    x, y = positions[i]
    rad = 0.2 + 0.3 * ldos_norm[i]   # scale circle radius
    color = plt.cm.viridis(ldos_norm[i])
    ax.add_patch(plt.Circle((x, y), rad, color=color, ec='black', linewidth=0.2))
ax.set_xlim(-Lx-0.5, Lx+0.5)
ax.set_ylim(-Ly-0.5, Ly+0.5)
ax.set_title('LDOS at E=0 (soliton channel)')
ax.set_xlabel('x (a_CC units)')
ax.set_ylabel('y (a_CC units)')

# --- (b) Total density of states ---
ax1 = axes[0, 1]
ax1.plot(E_range, DOS, 'k-')
ax1.axvline(0.0, color='blue', linestyle=':', label='E=0')
ax1.set_xlabel('Energy (eV)')
ax1.set_ylabel('DOS (per atom)')
ax1.set_title('Total Density of States')
ax1.legend()
ax1.grid(alpha=0.3)

# --- (c) Cross-section of LDOS along x (averaged over y) ---
ax2 = axes[1, 0]
# Bin atoms by x coordinate
bins = np.linspace(-Lx, Lx, 50)
x_centers = 0.5*(bins[:-1]+bins[1:])
ldos_x_mean = np.zeros(len(x_centers))
for i in range(N_atoms):
    ix = np.digitize(positions[i,0], bins) - 1
    if 0 <= ix < len(x_centers):
        ldos_x_mean[ix] += ldos_norm[i]
# Normalise per bin (count)
counts, _ = np.histogram(positions[:,0], bins=bins)
ldos_x_mean = np.divide(ldos_x_mean, counts, out=np.zeros_like(ldos_x_mean), where=counts>0)
ax2.plot(x_centers, ldos_x_mean, 'b-o', markersize=4)
ax2.axvline(x_wall, color='red', linestyle='--')
ax2.set_xlabel('x (a_CC units)')
ax2.set_ylabel('Average LDOS (normalised)')
ax2.set_title('LDOS cross-section across domain wall')
ax2.grid(alpha=0.3)

# --- (d) Eigenvalue spectrum (energy levels close to zero) ---
ax3 = axes[1, 1]
ax3.plot(E, 'k.', markersize=2)
ax3.axhline(0.0, color='blue', linestyle=':')
ax3.set_ylabel('Energy (eV)')
ax3.set_xlabel('State index')
ax3.set_title('Eigenvalue spectrum (zoom around zero in inset)')
# Inset: zoom
inset_ax = ax3.inset_axes([0.2, 0.6, 0.35, 0.35])
inset_ax.plot(E, 'k.', markersize=1)
inset_ax.set_xlim(N_atoms//2 - 20, N_atoms//2 + 20)
inset_ax.set_ylim(-0.5, 0.5)
inset_ax.axhline(0.0, color='blue', linestyle=':')
inset_ax.set_title('Zoom around E=0')
ax3.grid(alpha=0.3)

plt.tight_layout()
plt.show()