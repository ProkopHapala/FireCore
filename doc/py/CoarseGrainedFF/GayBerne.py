import numpy as np
import matplotlib.pyplot as plt

# ==========================================
# 1. ATOMISTIC BENZENE MODEL
# ==========================================

class AtomisticBenzene:
    def __init__(self):
        # Carbon ring radius approx 1.4 Angstrom
        self.R = 1.4 
        # 6 carbons in a hexagon
        angles = np.radians(np.arange(0, 360, 60))
        self.local_pos = np.column_stack([
            self.R * np.cos(angles),
            self.R * np.sin(angles),
            np.zeros(6)
        ])
        
        # OPLS-AA parameters for Aromatic Carbon
        # epsilon in kcal/mol, sigma in Angstroms
        self.eps_C = 0.070 
        self.sig_C = 3.55  
        
    def get_global_pos(self, center, rotation_matrix):
        # Rotate then translate
        return (rotation_matrix @ self.local_pos.T).T + center

def lj_potential(r, eps, sig):
    # Standard 12-6 LJ
    if r < 0.1: return 0 # avoid singularity at 0
    rr = (sig / r) ** 6
    return 4 * eps * (rr**2 - rr)

def compute_atomistic_energy(mol1_pos, mol2_pos, eps, sig):
    # O(N^2) brute force sum
    energy = 0.0
    for p1 in mol1_pos:
        for p2 in mol2_pos:
            r = np.linalg.norm(p1 - p2)
            energy += lj_potential(r, eps, sig)
    return energy

# ==========================================
# 2. GAY-BERNE (GB) MATH
# ==========================================

def gay_berne_parameters(r_vec, u1, u2, 
                         sig_0, sig_face, sig_edge, 
                         eps_0, chi, chi_prime, mu, nu):
    """
    Computes the anisotropic sigma and epsilon for two ellipsoids.
    u1, u2: symmetry axes (Normals to the disk face)
    r_vec: vector from 1 to 2
    """
    r_len = np.linalg.norm(r_vec)
    r_hat = r_vec / r_len
    
    # Dot products
    c1 = np.dot(u1, r_hat)
    c2 = np.dot(u2, r_hat)
    c12 = np.dot(u1, u2)
    
    # --- 1. SHAPE TERM (Sigma) ---
    # Berne-Pechukas formula
    # For benzene (oblate), u is the normal. 
    # sig_face (thickness), sig_edge (diameter).
    
    # Chi parameter for shape
    # standard form uses u as the "symmetry axis". For a disk, 
    # the symmetry axis is the normal. Length along normal is L (thickness),
    # Width perp is D (diameter).
    # chi = (L^2 - D^2) / (L^2 + D^2)
    # For benzene, L < D, so chi is negative.
    
    sigma_term = 0.0
    denom_plus = 1.0 + chi * c12
    denom_minus = 1.0 - chi * c12
    
    if abs(denom_plus) > 1e-9:
        sigma_term += (c1 + c2)**2 / denom_plus
    if abs(denom_minus) > 1e-9:
        sigma_term += (c1 - c2)**2 / denom_minus
        
    sigma = sig_0 / np.sqrt(1.0 - 0.5 * chi * sigma_term)

    # --- 2. ENERGY TERM (Epsilon) ---
    # Eps_1: dependence on relative orientation of axes
    eps1 = 1.0 / np.sqrt(1.0 - (chi * c12)**2)
    
    # Eps_2: dependence on overlap direction
    eps2_term = 0.0
    d_plus = 1.0 + chi_prime * c12
    d_minus = 1.0 - chi_prime * c12
    
    if abs(d_plus) > 1e-9:
        eps2_term += (c1 + c2)**2 / d_plus
    if abs(d_minus) > 1e-9:
        eps2_term += (c1 - c2)**2 / d_minus
        
    eps2 = 1.0 - 0.5 * chi_prime * eps2_term
    
    epsilon = eps_0 * (eps1 ** nu) * (eps2 ** mu)
    
    return sigma, epsilon

def morse_potential(r, r_eq, depth, a=1.5):
    # Generalized Morse
    # r_eq comes from GB sigma
    # depth comes from GB epsilon
    # a is stiffness (width of well)
    return depth * ((1.0 - np.exp(-a * (r - r_eq)))**2 - 1.0)

# ==========================================
# 3. RUN SIMULATION
# ==========================================

# -- Setup Benzene --
benz = AtomisticBenzene()

# -- GB Parameters Tuned for Benzene --
# Benzene is a pancake: Thick ~3.4A, Width ~6.5A
sig_face = 3.4   # Thickness (Simulates stacking dist)
sig_edge = 6.0   # Diameter (Simulates edge-edge dist)
sig_0 = sig_face # Reference scale

# Anisotropy params
# chi = (l^2 - d^2)/(l^2 + d^2). l=face, d=edge.
chi = (sig_face**2 - sig_edge**2) / (sig_face**2 + sig_edge**2)

# Energy anisotropy
# Stacking is stronger than T-shape or Edge-Edge
eps_face = 3.0   # Strong stacking
eps_edge = 0.3   # Weak edge-touching
# chi_prime = (eps_face^(1/mu) - eps_edge^(1/mu)) / ...
# Empirically tuning chi_prime for this demo:
chi_prime = -0.8 # Negative favors parallel vectors (stacking)

mu = 2
nu = 1

# Scaling factor to match atomistic sum magnitude approximately
# (Atomistic sum has 36 pairs, GB has 1 pair)
GB_SCALE = 1.0 

# Define Path
distances = np.linspace(2.5, 9.0, 100)
u1 = np.array([0, 0, 1]) # Normal of Mol 1 (Fixed)

# --- SCENARIO 1: Face-to-Face (Stacking) ---
# Mol 2 moves along Z, Normal is [0,0,1] (Parallel)
u2_stack = np.array([0, 0, 1])
energies_atom_stack = []
energies_gb_stack = []

for d in distances:
    # Atomistic
    pos1 = benz.get_global_pos(np.array([0,0,0]), np.eye(3))
    pos2 = benz.get_global_pos(np.array([0,0,d]), np.eye(3))
    e_atom = compute_atomistic_energy(pos1, pos2, benz.eps_C, benz.sig_C)
    energies_atom_stack.append(e_atom)
    
    # Gay-Berne (Morse)
    r_vec = np.array([0, 0, d])
    sig_gb, eps_gb = gay_berne_parameters(r_vec, u1, u2_stack, 
                                          sig_0, sig_face, sig_edge, 
                                          eps_face, chi, chi_prime, mu, nu)
    
    # NOTE: In GB-LJ, we typically use shifted distance. 
    # For GB-Morse, we treat sig_gb as the equilibrium position r_eq.
    e_gb = morse_potential(d, sig_gb, eps_gb, a=1.8)
    energies_gb_stack.append(e_gb)

# --- SCENARIO 2: Edge-to-Edge (In-Plane) ---
# Mol 2 moves along X, Normal is [0,0,1] (Parallel normals, but side-by-side pos)
u2_edge = np.array([0, 0, 1])
energies_atom_edge = []
energies_gb_edge = []

distances_edge = np.linspace(5.0, 12.0, 100)

for d in distances_edge:
    # Atomistic
    pos1 = benz.get_global_pos(np.array([0,0,0]), np.eye(3))
    pos2 = benz.get_global_pos(np.array([d,0,0]), np.eye(3)) # Move in X
    e_atom = compute_atomistic_energy(pos1, pos2, benz.eps_C, benz.sig_C)
    energies_atom_edge.append(e_atom)
    
    # Gay-Berne (Morse)
    r_vec = np.array([d, 0, 0])
    # Note: Normal is still [0,0,1], but r_hat is [1,0,0]
    # This effectively samples the "Side" of the disk
    sig_gb, eps_gb = gay_berne_parameters(r_vec, u1, u2_edge, 
                                          sig_0, sig_face, sig_edge, 
                                          eps_face, chi, chi_prime, mu, nu)
    
    e_gb = morse_potential(d, sig_gb, eps_gb, a=1.8)
    energies_gb_edge.append(e_gb)

# ==========================================
# 4. PLOTTING
# ==========================================

plt.figure(figsize=(12, 5))

# Plot Stacking
plt.subplot(1, 2, 1)
plt.title("Face-to-Face Stacking (Z-axis)")
plt.plot(distances, energies_atom_stack, 'k-', lw=2, label="Atomistic Sum (LJ)")
plt.plot(distances, energies_gb_stack, 'r--', lw=2, label="Gay-Berne (Morse)")
plt.xlabel("Distance (A)")
plt.ylabel("Energy (kcal/mol)")
plt.legend()
plt.grid(True, alpha=0.3)
plt.ylim(-5, 2)

# Plot Edge-to-Edge
plt.subplot(1, 2, 2)
plt.title("Edge-to-Edge (X-axis)")
plt.plot(distances_edge, energies_atom_edge, 'k-', lw=2, label="Atomistic Sum (LJ)")
plt.plot(distances_edge, energies_gb_edge, 'r--', lw=2, label="Gay-Berne (Morse)")
plt.xlabel("Distance (A)")
plt.legend()
plt.grid(True, alpha=0.3)
plt.ylim(-5, 2)

plt.tight_layout()
plt.show()

# Print out check
print(f"Stacking Minimums -> Atomistic: {min(energies_atom_stack):.2f}, GB: {min(energies_gb_stack):.2f}")
print(f"Edge Minimums     -> Atomistic: {min(energies_atom_edge):.2f}, GB: {min(energies_gb_edge):.2f}")