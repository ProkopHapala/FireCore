import numpy as np
import matplotlib.pyplot as plt

# --- Physical Parameters for Polyacetylene (SSH Model) ---
N = 100          # Number of atoms (even)
t0 = -2.5        # Base hopping integral (eV)
alpha = 4.1      # Electron-phonon coupling (eV / Angstrom)
K = 21.0         # Sigma-bond spring constant (eV / Angstrom^2)

def build_hamiltonian(N, u):
    """
    Builds the 1D tight-binding Hamiltonian with Peierls distortion u.
    u > 0 implies dimerization (alternating short/long bonds).
    """
    H = np.zeros((N, N))
    for i in range(N - 1):
        # Alternating bond distortion: (-1)^i
        # If i is even, atoms i and i+1 move closer by 2u
        delta_r = 2 * u if i % 2 == 0 else -2 * u
        t = t0 - alpha * delta_r
        
        H[i, i+1] = t
        H[i+1, i] = t
    return H

# =====================================================================
# EXACT METHOD: O(N^3) Full Diagonalization
# =====================================================================
def get_energy_exact(H):
    evals, evecs = np.linalg.eigh(H)
    # Fill the bottom half of the states (mu = 0)
    occupied_evals = evals[evals < 0.0]
    # Factor of 2 for spin
    return 2.0 * np.sum(occupied_evals)

# =====================================================================
# METHOD 1: O(N) Chebyshev (Fermi Operator Expansion)
# =====================================================================
def get_energy_chebyshev(H, M=40):
    """
    Estimates total pi-energy using Chebyshev expansion of the density matrix.
    M is the polynomial order (determines interaction range).
    """
    N = H.shape[0]
    
    # 1. Scale H so its spectrum lies strictly inside [-1, 1]
    # For a 1D chain, the spectral radius is bounded by 2 * max(|t|)
    E_max = 2.0 * np.max(np.abs(H)) + 0.5 
    H_scaled = H / E_max
    
    # 2. Calculate Jackson-damped Chebyshev coefficients for Step function at E=0
    # Jackson damping prevents unphysical "Gibbs ringing" oscillations
    c = np.zeros(M)
    for k in range(1, M, 2): # Only odd terms are non-zero for E=0 step
        # Analytical integral for step function
        c_k_raw = -(2.0 / (np.pi * k)) * np.sin(k * np.pi / 2.0)
        
        # Jackson damping factor
        alpha_ang = np.pi / (M + 1)
        gk = ((M - k + 1) * np.cos(k * alpha_ang) + np.sin(k * alpha_ang) * (1 / np.tan(alpha_ang))) / (M + 1)
        
        c[k] = c_k_raw * gk
        
    # 3. Compute Energy efficiently using Matrix-Vector operations
    E_pi = 0.0
    
    # We only need rho_{i, i+1} to compute energy. We apply the polynomial to local basis vectors.
    for i in range(N - 1):
        # Start vector: |v0> = |i>
        v0 = np.zeros(N)
        v0[i] = 1.0
        
        # T_0(H)|v0> = |v0>
        v_prev = v0
        
        # T_1(H)|v0> = H|v0>
        v_curr = H_scaled @ v0
        
        rho_i_i1 = c[1] * v_curr[i+1] # Contribution from k=1
        
        for k in range(2, M):
            # T_{k}(H) = 2H T_{k-1}(H) - T_{k-2}(H)  -> Sparse Matrix-Vector multiply!
            v_next = 2.0 * (H_scaled @ v_curr) - v_prev
            
            if k % 2 != 0: # Only add odd terms
                rho_i_i1 += c[k] * v_next[i+1]
                
            v_prev = v_curr
            v_curr = v_next
            
        # Add to total energy: 2 (spin) * t_ij * rho_ji
        # Wait, step function coefficients give us half the density. 
        # Total density includes factor of 2 for spin.
        E_pi += 2.0 * 2.0 * H[i, i+1] * rho_i_i1
        
    return E_pi

# =====================================================================
# METHOD 2: O(N) Lanczos Recursion
# =====================================================================
def get_lanczos_coeffs(H, v0, M):
    """Standard Lanczos tridiagonalization: returns coefficients (a_n, b_n).
        These define the continued fraction representation of the local
        Green's function G(E) = 1 / (E - a_0 - b_1^2 / (E - a_1 - ...))."""
    N = H.shape[0]
    a = np.zeros(M)
    b = np.zeros(M)
    v_prev = np.zeros(N)
    v_curr = v0 / np.linalg.norm(v0)
    for k in range(M):
        w = H @ v_curr
        a[k] = np.dot(v_curr, w)
        w = w - a[k] * v_curr - (b[k-1] * v_prev if k > 0 else 0)
        if k < M - 1:
            b[k] = np.linalg.norm(w)
            if b[k] < 1e-10:
                break  # Krylov subspace exhausted
            v_prev = v_curr
            v_curr = w / b[k]
    return a, b

def eval_continued_fraction(E, a, b, M):
    """Evaluate continued fraction at complex energy E (bottom-up)."""
    G = 0.0 + 0.0j
    for k in range(M - 1, -1, -1):
        if k == M - 1:
            G = 1.0 / (E - a[k])
        else:
            G = 1.0 / (E - a[k] - (b[k] ** 2) * G)
    return G

def get_bond_density_CF(H, i, j, M, eta=0.2, E_range=(-10, 0.0), n_points=200):
    """Compute bond density rho_ij using Lanczos/Continued Fraction (Haydock).
        Uses symmetric |+> and antisymmetric |-> bond-centered states,
        evaluates local Green's function via continued fraction, and integrates
        the spectral density to get occupation."""
    N = H.shape[0]
    v_plus = np.zeros(N)
    v_plus[i] = 1.0 / np.sqrt(2)
    v_plus[j] = 1.0 / np.sqrt(2)
    a_plus, b_plus = get_lanczos_coeffs(H, v_plus, M)

    v_minus = np.zeros(N)
    v_minus[i] = 1.0 / np.sqrt(2)
    v_minus[j] = -1.0 / np.sqrt(2)
    a_minus, b_minus = get_lanczos_coeffs(H, v_minus, M)

    E_grid = np.linspace(E_range[0], E_range[1], n_points)
    dE = E_grid[1] - E_grid[0]
    rho_plus = 0.0
    rho_minus = 0.0
    for E in E_grid:
        z = E + 1j * eta
        G_plus = eval_continued_fraction(z, a_plus, b_plus, M)
        G_minus = eval_continued_fraction(z, a_minus, b_minus, M)
        rho_plus += -(1.0 / np.pi) * np.imag(G_plus) * dE
        rho_minus += -(1.0 / np.pi) * np.imag(G_minus) * dE
    # Algebraic identity: rho_ij = (rho_+ - rho_-)/2  (single-spin)
    return 0.5 * (rho_plus - rho_minus)

def get_energy_lanczos(H, M=15):
    """Total pi energy from Lanczos/Continued Fraction bond densities."""
    E_pi = 0.0
    N = H.shape[0]
    for i in range(N - 1):
        rho_i_i1 = get_bond_density_CF(H, i, i + 1, M)
        E_pi += 2.0 * 2.0 * H[i, i + 1] * rho_i_i1
    return E_pi

# =====================================================================
# RUN SIMULATION AND PLOT
# =====================================================================
u_vals = np.linspace(0.0, 0.1, 20)

E_exact = []
E_cheb = []
E_lanczos = []

for u in u_vals:
    H = build_hamiltonian(N, u)
    
    # 1. Pi Energy
    e_pi_ex = get_energy_exact(H)
    e_pi_ch = get_energy_chebyshev(H, M=40)
    e_pi_la = get_energy_lanczos(H, M=15)
    
    # 2. Sigma Spring Energy
    # K is spring constant. Spring energy per dimer is 2 * (1/2) K (2u)^2 ...
    # We sum over all bonds.
    E_sigma = np.sum(0.5 * K * (2 * u)**2 * np.ones(N-1))
    
    # We plot the change in energy relative to u=0
    E_exact.append(e_pi_ex + E_sigma)
    E_cheb.append(e_pi_ch + E_sigma)
    E_lanczos.append(e_pi_la + E_sigma)

# Normalize to Delta E
E_exact = np.array(E_exact) - E_exact[0]
E_cheb = np.array(E_cheb) - E_cheb[0]
E_lanczos = np.array(E_lanczos) - E_lanczos[0]

plt.figure(figsize=(8, 6))
plt.plot(u_vals, E_exact, 'k-', linewidth=2, label="Exact Diagonalization")
plt.plot(u_vals, E_cheb, 'r--', linewidth=2, label="Chebyshev FOE (O(N), M=40)")
plt.plot(u_vals, E_lanczos, 'b-.', linewidth=2, label="Lanczos Recursion (O(N), M=15)")
plt.axhline(0, color='gray', linestyle=':')
plt.xlabel("Peierls Distortion $u$ (Angstrom)", fontsize=12)
plt.ylabel("Change in Total Energy $\Delta E$ (eV)", fontsize=12)
plt.title("Peierls Dimerization in Polyacetylene ($N=100$)", fontsize=14)
plt.legend(fontsize=12)
plt.grid(True)
plt.show()