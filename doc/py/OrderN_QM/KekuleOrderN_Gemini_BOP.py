import numpy as np
import matplotlib.pyplot as plt

# =====================================================================
# SYSTEM SETUP: 1D Carbon Chain
# =====================================================================
N = 40          # Finite chain length
t0 = -2.5       # Uniform hopping (eV)
u_peierls = 0.05 # Peierls distortion (dimerization)
alpha = 4.1      

# Construct Hamiltonian with dimerization
H = np.zeros((N, N))
for i in range(N - 1):
    t = t0 - alpha * (2 * u_peierls if i % 2 == 0 else -2 * u_peierls)
    H[i, i+1] = H[i+1, i] = t

# We need the scaled Hamiltonian for Chebyshev propagation (Eigenvalues must be in [-1, 1])
E_max = 2.0 * np.max(np.abs(H)) + 0.1
H_scaled = H / E_max

# Calculate Exact Density Matrix for comparison
evals, evecs = np.linalg.eigh(H)
rho_exact = np.zeros((N, N))
for lam in range(N):
    if evals[lam] < 0.0:
        rho_exact += 2.0 * np.outer(evecs[:, lam], evecs[:, lam]) # x2 for spin

# =====================================================================
# ALGORITHM: Propagate Local Probe & Track Interference
# =====================================================================
def get_chebyshev_coefficients(M):
    """Calculates Jackson-damped coefficients for the Step Function at E=0"""
    c = np.zeros(M)
    for k in range(1, M, 2): # Only ODD coefficients matter for E=0 step
        # Base analytical coefficient for step function
        c_raw = -(2.0 / (np.pi * k)) * np.sin(k * np.pi / 2.0)
        
        # Jackson damping (removes Gibbs ringing for small M)
        ang = np.pi / (M + 1)
        damping = ((M - k + 1) * np.cos(k * ang) + np.sin(k * ang) / np.tan(ang)) / (M + 1)
        c[k] = c_raw * damping
    return c

def propagate_local_probe(H_scaled, start_site, M):
    """
    Drops a 'pebble' at start_site and watches the wave spread.
    Returns the history of the wave vector at each hop iteration.
    """
    N = H_scaled.shape[0]
    wave_history = np.zeros((M, N))
    
    # Hop 0: The initial localized pebble
    v_prev = np.zeros(N)
    v_curr = np.zeros(N)
    v_curr[start_site] = 1.0
    wave_history[0] = v_curr
    
    # Hop 1: First spread
    v_next = H_scaled @ v_curr
    wave_history[1] = v_next
    
    # Hop k: Chebyshev wave spreading (smart path counting)
    for k in range(2, M):
        v_next_next = 2.0 * (H_scaled @ v_next) - v_curr
        wave_history[k] = v_next_next
        
        # Shift vectors for next iteration
        v_curr = v_next
        v_next = v_next_next
        
    return wave_history

# =====================================================================
# EXPERIMENT: Watch the wave spread from Site 20 and measure Bond 20-21
# =====================================================================
M_max = 16
start_site = 20
target_site = 21 # We want the bond density between 20 and 21

# 1. Propagate the wave
wave_history = propagate_local_probe(H_scaled, start_site, M_max)

# 2. Get the coefficients
c = get_chebyshev_coefficients(M_max)

# 3. Sum the interferences at the target site to build the bond density
bond_density_contributions = np.zeros(M_max)
accumulated_bond_density = np.zeros(M_max)

current_density = 0.0
for k in range(M_max):
    # Contribution = Coefficient * Amplitude of wave at target site
    # Note: wave_history[k, target_site] is exactly <target | T_k(H) | start>
    contrib = c[k] * wave_history[k, target_site] * 2.0 # x2 for spin
    
    bond_density_contributions[k] = contrib
    current_density += contrib
    accumulated_bond_density[k] = current_density

# =====================================================================
# VISUALIZATIONS
# =====================================================================
fig = plt.figure(figsize=(14, 10))

# --- PLOT 1: The Wave Spreading (The "Splash") ---
ax1 = plt.subplot(2, 2, 1)
hops_to_plot = [1, 3, 5, 7]
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']

for i, k in enumerate(hops_to_plot):
    # Plot absolute wave amplitude
    ax1.plot(np.arange(N), np.abs(wave_history[k]), marker='o', 
             linestyle='-', color=colors[i], label=f'Hop $k={k}$', alpha=0.8)

ax1.axvline(start_site, color='k', linestyle='--', label='Start Site (20)')
ax1.set_xlim(10, 30)
ax1.set_title("1. Wave Propagation Outward from Local Probe")
ax1.set_xlabel("Atom Index")
ax1.set_ylabel("Wave Amplitude $|v_k|$")
ax1.legend()
ax1.grid(True, alpha=0.3)

# --- PLOT 2: How Paths Build the Bond Density ---
ax2 = plt.subplot(2, 2, 2)
# We only plot odd hops since even hops have c_k = 0 for bond density
odd_hops = np.arange(1, M_max, 2)
contributions = bond_density_contributions[odd_hops]

ax2.bar(odd_hops, contributions, color='purple', alpha=0.7)
ax2.axhline(0, color='k', linewidth=1)
ax2.set_title(f"2. Scattering Path Contributions to $\\rho_{{{start_site},{target_site}}}$")
ax2.set_xlabel("Scattering Path Length (Hops $k$)")
ax2.set_ylabel("Contribution to Bond Density")
ax2.set_xticks(odd_hops)
ax2.grid(axis='y', alpha=0.3)
ax2.text(1.5, contributions[0]*0.8, "Direct\nHop", ha='left', fontsize=10)
ax2.text(3.5, contributions[1]*0.8, "1st\nReflection", ha='left', fontsize=10)

# --- PLOT 3: Convergence of the Bond Density vs Iterations ---
ax3 = plt.subplot(2, 1, 2)
exact_val = rho_exact[start_site, target_site]

ax3.plot(np.arange(M_max), accumulated_bond_density, 'bo-', linewidth=2, label='Accumulated Density (Probe method)')
ax3.axhline(exact_val, color='r', linestyle='--', linewidth=2, label='Exact Diagonalization (Target)')
ax3.fill_between(np.arange(M_max), accumulated_bond_density, exact_val, color='blue', alpha=0.1)

ax3.set_title(f"3. Convergence of Bond Density $\\rho_{{{start_site},{target_site}}}$ vs Propagation Iterations (M)")
ax3.set_xlabel("Total Iterations / Max Hop Distance")
ax3.set_ylabel("Estimated Bond Density")
ax3.set_xticks(np.arange(M_max))
ax3.legend()
ax3.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()