import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

# =====================================================================
# SYSTEM SETUP: Polyacetylene with Boundaries and a Defect
# =====================================================================
N = 40          # Finite chain length
t0 = -2.5       # Uniform hopping (eV)
u_peierls = 0.05 # Peierls distortion amplitude
alpha = 4.1      

# Construct Hamiltonian
H = np.zeros((N, N))
for i in range(N - 1):
    # Alternating distortion: short, long, short, long...
    t = t0 - alpha * (2 * u_peierls if i % 2 == 0 else -2 * u_peierls)
    H[i, i+1] = H[i+1, i] = t

# Introduce a Nitrogen Defect at atom 20
# Nitrogen is more electronegative, so its onsite energy is lower
H[20, 20] = -1.5 

# Exact Diagonalization (Baseline)
evals, evecs = np.linalg.eigh(H)
rho_exact_matrix = np.zeros((N, N))
for lam in range(N):
    if evals[lam] < 0.0: # Occupied states
        # Factor of 2 for spin
        rho_exact_matrix += 2.0 * np.outer(evecs[:, lam], evecs[:, lam])

rho_bonds_exact = np.array([rho_exact_matrix[i, i+1] for i in range(N-1)])

# =====================================================================
# METHOD A: Haydock Recursion (Algebraic Identity + Continued Fraction)
# =====================================================================
def get_lanczos_coeffs(H, v0, M):
    """Generates the a_n, b_n coefficients for the Continued Fraction."""
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
            if b[k] < 1e-10: break
            v_prev = v_curr
            v_curr = w / b[k]
    return a, b

def eval_continued_fraction(E, a, b, M):
    """Evaluates the Continued Fraction at a complex energy E."""
    # We evaluate from the bottom up
    G = 0.0 + 0.0j
    for k in range(M-1, -1, -1):
        if k == M-1:
            # Terminator (simulates infinite chain broadening)
            G = 1.0 / (E - a[k]) 
        else:
            G = 1.0 / (E - a[k] - (b[k]**2) * G)
    return G

def get_bond_density_CF(H, i, j, M, eta=0.2):
    """Computes rho_ij using the symmetric/antisymmetric trick."""
    N = H.shape[0]
    
    # 1. Symmetric state |+>
    v_plus = np.zeros(N)
    v_plus[i], v_plus[j] = 1.0/np.sqrt(2), 1.0/np.sqrt(2)
    a_plus, b_plus = get_lanczos_coeffs(H, v_plus, M)
    
    # 2. Antisymmetric state |->
    v_minus = np.zeros(N)
    v_minus[i], v_minus[j] = 1.0/np.sqrt(2), -1.0/np.sqrt(2)
    a_minus, b_minus = get_lanczos_coeffs(H, v_minus, M)
    
    # 3. Integrate Continued Fractions (Riemann sum for speed/robustness)
    # Using a small imaginary broadening eta
    E_grid = np.linspace(-10, 0.0, 200)
    dE = E_grid[1] - E_grid[0]
    
    rho_plus = 0.0
    rho_minus = 0.0
    for E in E_grid:
        z = E + 1j * eta
        G_plus = eval_continued_fraction(z, a_plus, b_plus, M)
        G_minus = eval_continued_fraction(z, a_minus, b_minus, M)
        
        rho_plus += -(1.0/np.pi) * np.imag(G_plus) * dE
        rho_minus += -(1.0/np.pi) * np.imag(G_minus) * dE
        
    # The Algebraic identity! (Multiply by 2 for spin)
    return 2.0 * 0.5 * (rho_plus - rho_minus)

# Compute Continued Fraction densities for all bonds at M=3
rho_bonds_CF_M3 = np.array([get_bond_density_CF(H, i, i+1, M=3) for i in range(N-1)])

# =====================================================================
# METHOD B: Analytical Convolution Cheat (O(1) per bond runtime)
# =====================================================================
def get_analytical_chi(max_m, t0, eta=0.05):
    """Precomputes the spatial decay of bond density chi(m) offline."""
    chi = np.zeros(max_m)
    E_grid = np.linspace(-4*np.abs(t0), 0.0, 1000)
    dE = E_grid[1] - E_grid[0]
    
    for m in range(max_m):
        integral = 0.0
        for E in E_grid:
            z = (E + 1j * eta) / (2.0 * t0)
            root = np.sqrt(z**2 - 1.0 + 0j)
            if np.imag(root) < 0: root = -root # Ensure positive imaginary part
            
            # Analytical Green's function for uniform chain
            def G(dist):
                return (-1.0 / (2.0 * t0 * root)) * (z - root)**dist
            
            G_m = G(m)
            G_mp1 = G(m+1)
            G_mm1 = G(np.abs(m-1))
            
            # Formula for bond-bond susceptibility
            integral += -(1.0/np.pi) * np.imag(G_m**2 + G_mp1*G_mm1) * dE
        chi[m] = integral * 2.0 # Factor of 2 for spin
    return chi

# Precompute convolution kernel
max_m_conv = 10
chi_kernel = get_analytical_chi(max_m_conv, t0)

# Runtime step: Apply Convolution (Super fast!)
# Base density for infinite uniform chain is 4/pi = 1.273 (for 2 spins)
rho_base = (4.0 / np.pi) 
rho_bonds_conv = np.ones(N-1) * rho_base

for i in range(N-1):
    delta_rho = 0.0
    # Slide the window
    for m in range(-max_m_conv+1, max_m_conv):
        k = i + m
        if 0 <= k < N-1:
            # How much does bond k deviate from uniform hopping t0?
            delta_t_k = H[k, k+1] - t0
            delta_rho += chi_kernel[np.abs(m)] * delta_t_k
            
    rho_bonds_conv[i] += delta_rho


# =====================================================================
# VISUALIZATION & CONVERGENCE PLOTS
# =====================================================================

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))

# --- PLOT 1: Spatial Bond Density (Defects and Boundaries) ---
bonds = np.arange(N-1)
ax1.plot(bonds, rho_bonds_exact, 'ko-', linewidth=2, label='Exact Diagonalization')
ax1.plot(bonds, rho_bonds_CF_M3, 'r^-', markersize=6, label='Continued Fraction (M=3)')
ax1.plot(bonds, rho_bonds_conv, 'bx--', markersize=6, label='Analytical Convolution')

ax1.axvline(20, color='gray', linestyle=':', label='Nitrogen Defect Location')
ax1.set_title("Bond Density Profiles along the Chain", fontsize=14)
ax1.set_xlabel("Bond Index", fontsize=12)
ax1.set_ylabel(r"Bond Density $\rho_{i, i+1}$", fontsize=12)
ax1.legend()
ax1.grid(True)

# Annotations for physical interpretation
ax1.text(2, 1.35, "Friedel Oscillations\nfrom Boundary", color='black', fontsize=10)
ax1.text(22, 1.4, "Defect breaks\nlocal symmetry", color='black', fontsize=10)

# --- PLOT 2: Convergence with Iterations (M) ---
M_vals = np.arange(1, 9)
error_middle_bond = []

# Check convergence at a specific bond far from the defect (e.g., bond 10)
exact_val = rho_bonds_exact[10]
for M in M_vals:
    val_CF = get_bond_density_CF(H, 10, 11, M=M, eta=0.2)
    error = np.abs(val_CF - exact_val)
    error_middle_bond.append(error)

ax2.plot(M_vals, error_middle_bond, 'ro-', linewidth=2, markersize=8)
ax2.set_yscale('log')
ax2.set_title("Convergence of Continued Fraction Method vs Recursion Depth (M)", fontsize=14)
ax2.set_xlabel("Number of Recursion Iterations / Vector Multiplies (M)", fontsize=12)
ax2.set_ylabel("Absolute Error in $\\rho_{10,11}$ (Log Scale)", fontsize=12)
ax2.set_xticks(M_vals)
ax2.grid(True, which="both", ls="--")

plt.tight_layout()
plt.show()