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
# 3. EXACT SOLUTION: O(N^3) Diagonalization (reference)
# =====================================================================
evals, evecs = np.linalg.eigh(H)
rho_exact = np.zeros((N, N))
for lam in range(N):
    if evals[lam] < 0.0: # Half-filled (mu = 0)
        rho_exact += 2.0 * np.outer(evecs[:, lam], evecs[:, lam]) # x2 for spin

print(f"System size N={N}, bonds={len(bonds)}")
print(f"Spectral range: [{evals[0]:.2f}, {evals[-1]:.2f}] eV")

# =====================================================================
# 4. METHOD 1: Chebyshev / Fermi Operator Expansion (FOE)
#    (Original approach from v1)
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

# Scale H to [-1, 1] range
E_max_scaled = 3.0 * np.abs(t0) + np.abs(e_nitrogen) + 0.5
H_scaled = H / E_max_scaled

def get_bond_density_cheb(H_scaled, start_site, target_site, M):
    """Probe method: drop probe at start_site, measure amplitude at target_site
       using Chebyshev polynomial expansion of the density matrix.
       ρ_ij = Σ c_k * <j|T_k(H)|i> * 2 (for spin)

       Only odd k contribute because:
       - c_k = 0 for even k (step function at E=0)
       - For bipartite lattice, <j|T_k(H)|i> = 0 for even k when i,j neighbors"""
    N = H_scaled.shape[0]
    c = get_chebyshev_coeffs(M)

    # Need at least M >= 2 to have c[1] defined
    if M < 2:
        return 0.0

    v_prev = np.zeros(N)
    v_curr = np.zeros(N)
    v_curr[start_site] = 1.0
    current_density = 0.0

    # k=1: T_1(H) = H, first hop from start to neighbors
    v_next = H_scaled @ v_curr
    current_density += c[1] * v_next[target_site] * 2.0

    # k=2..M-1: T_k(H) = 2H·T_{k-1}(H) - T_{k-2}(H)
    for k in range(2, M):
        v_next_next = 2.0 * (H_scaled @ v_next) - v_curr
        if k % 2 != 0:  # Only odd terms contribute for bipartite lattice
            current_density += c[k] * v_next_next[target_site] * 2.0
        v_prev, v_curr, v_next = v_curr, v_next, v_next_next

    return current_density

# =====================================================================
# 5. METHOD 2: Lanczos / Continued Fraction (Haydock Recursion)
#    The proper Lanczos method from the reference:
#    - Uses bond-centered symmetric/antisymmetric states
#    - Builds tridiagonal Lanczos representation
#    - Evaluates local Green's function via continued fraction
#    - Integrates spectral density to get ρ_+ and ρ_-
#    - ρ_ij = (ρ_+ - ρ_-)/2
# =====================================================================
def get_lanczos_coeffs(H, v0, M):
    """Generates Lanczos tridiagonal coefficients (a_n, b_n) for starting vector v0.
       These define the continued fraction representation of the local Green's function:
       G(E) = 1 / (E - a_0 - b_1^2 / (E - a_1 - b_2^2 / (...)))"""
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
    """Evaluates the continued fraction at complex energy E (bottom-up).
       G(E) = 1 / (E - a_0 - b_1^2 / (E - a_1 - b_2^2 / (... - b_{M-1}^2/(E - a_{M-1})...)))"""
    G = 0.0 + 0.0j
    for k in range(M - 1, -1, -1):
        if k == M - 1:
            # Square-root terminator: G = 1 / (E - a - Σ(E))
            # Simple terminator: just truncate at level M-1
            G = 1.0 / (E - a[k])
        else:
            G = 1.0 / (E - a[k] - (b[k] ** 2) * G)
    return G

def get_bond_density_CF(H, i, j, M, eta=0.2, E_range=(-10, 0.0), n_points=200):
    """Computes ρ_ij using the Lanczos/Continued Fraction method:
       - Construct bond-centered symmetric |+> and antisymmetric |-> states
       - Get Lanczos coefficients for each
       - Evaluate G_+(E) and G_-(E) via continued fraction
       - Integrate spectral density: ρ = ∫ -(1/π) Im[G(E)] dE
       - ρ_ij = (ρ_+ - ρ_-)/2  (times 2 for spin)"""
    N = H.shape[0]

    # 1. Symmetric state |+> = (|i> + |j>)/√2
    v_plus = np.zeros(N)
    v_plus[i] = 1.0 / np.sqrt(2)
    v_plus[j] = 1.0 / np.sqrt(2)
    a_plus, b_plus = get_lanczos_coeffs(H, v_plus, M)

    # 2. Antisymmetric state |-> = (|i> - |j>)/√2
    v_minus = np.zeros(N)
    v_minus[i] = 1.0 / np.sqrt(2)
    v_minus[j] = -1.0 / np.sqrt(2)
    a_minus, b_minus = get_lanczos_coeffs(H, v_minus, M)

    # 3. Integrate LDOS via Riemann sum over occupied states (E < 0)
    # Note: E_range should cover all occupied eigenvalues. For safety, compute
    #       spectral bounds from Gershgorin: |E| ≤ max_i(Σ_j |H[i,j]|)
    #       Or use a generously wide range like [-10, 0] for typical C pi-systems
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

    # 4. Algebraic identity: ρ_ij = (ρ_+ - ρ_-)/2, times 2 for spin
    return 2.0 * 0.5 * (rho_plus - rho_minus)

# =====================================================================
# 6. COMPUTE BOND DENSITIES FOR ALL METHODS
# =====================================================================
M_target = 15

rho_cheb = np.zeros(len(bonds))
rho_cf = np.zeros(len(bonds))

print(f"\nComputing bond densities (M={M_target})...")
for idx, (i, j) in enumerate(bonds):
    rho_cheb[idx] = get_bond_density_cheb(H_scaled, i, j, M_target)
    rho_cf[idx] = get_bond_density_CF(H, i, j, M_target)
    if idx % 5 == 0:
        print(f"  Bond {idx}/{len(bonds)}: exact={rho_exact[i,j]:.4f}  cheb={rho_cheb[idx]:.4f}  cf={rho_cf[idx]:.4f}")

# Compute errors
rho_exact_vals = np.array([rho_exact[i, j] for i, j in bonds])
err_cheb = np.abs(rho_cheb - rho_exact_vals)
err_cf = np.abs(rho_cf - rho_exact_vals)
print(f"\nMean absolute errors:")
print(f"  Chebyshev:  {err_cheb.mean():.6f}  (max: {err_cheb.max():.6f})")
print(f"  Lanczos/CF: {err_cf.mean():.6f}  (max: {err_cf.max():.6f})")

# =====================================================================
# 7. SELECT SPECIFIC BONDS FOR DETAILED CONVERGENCE ANALYSIS
# =====================================================================
edge_bond = None
core_bond = None
defect_bond = None

for idx, (i, j) in enumerate(bonds):
    if i == defect_site or j == defect_site:
        defect_bond = (idx, i, j)
    elif np.linalg.norm((coords[i] + coords[j]) / 2) > 2.5:
        edge_bond = (idx, i, j)
    elif np.linalg.norm((coords[i] + coords[j]) / 2) < 1.5:
        core_bond = (idx, i, j)

tracked_bonds = {
    "Edge Bond (Strong Kekulé)": edge_bond,
    "Core C-C Bond (Delocalized)": core_bond,
    "C-N Defect Bond": defect_bond
}
tracked_bonds_filtered = {k: v for k, v in tracked_bonds.items() if v is not None}

# =====================================================================
# 8. CONVERGENCE STUDY BOTH METHODS
# =====================================================================
M_max_conv = 20
M_range = np.arange(1, M_max_conv + 1)

cheb_conv = {}
cf_conv = {}
exact_vals = {}

for name, (idx, i, j) in tracked_bonds_filtered.items():
    # Chebyshev convergence
    vals_cheb = []
    for M in M_range:
        vals_cheb.append(get_bond_density_cheb(H_scaled, i, j, M))
    cheb_conv[name] = np.array(vals_cheb)

    # Lanczos/CF convergence
    vals_cf = []
    for M in M_range:
        vals_cf.append(get_bond_density_CF(H, i, j, M))
    cf_conv[name] = np.array(vals_cf)

    exact_vals[name] = rho_exact[i, j]

# =====================================================================
# 9. VISUALIZATION
# =====================================================================
color_list = ['green', 'purple', 'orange']

fig = plt.figure(figsize=(18, 14))

# --- PLOT A: Molecular Graph - Overlay of all methods ---
axA = plt.subplot(2, 3, 1)
segments = [(coords[i], coords[j]) for (i, j) in bonds]

# Exact (black solid, reference)
widths_exact = (rho_exact_vals - 0.3) * 12
lc_exact = LineCollection(segments, linewidths=widths_exact, colors='black',
                          alpha=0.8, label='Exact (reference)', zorder=4)
axA.add_collection(lc_exact)

# Chebyshev (red dashed, offset slightly)
widths_cheb = (rho_cheb - 0.3) * 12
lc_cheb = LineCollection(segments, linewidths=widths_cheb, colors='red',
                         linestyle='--', alpha=0.4, label=f'Chebyshev (M={M_target})', zorder=3)
axA.add_collection(lc_cheb)

# Lanczos/CF (blue dotted)
widths_cf = (rho_cf - 0.3) * 12
lc_cf = LineCollection(segments, linewidths=widths_cf, colors='blue',
                        linestyle=':', alpha=0.6, label=f'Lanczos/CF (M={M_target})', zorder=5)
axA.add_collection(lc_cf)

# Atoms
axA.scatter(coords[:, 0], coords[:, 1], c='black', s=40, zorder=6)
axA.scatter(coords[defect_site, 0], coords[defect_site, 1], c='red', s=120, zorder=7, label='N Defect')

# Highlight tracked bonds
for c_idx, (name, (idx, i, j)) in enumerate(tracked_bonds_filtered.items()):
    axA.plot([coords[i, 0], coords[j, 0]], [coords[i, 1], coords[j, 1]],
             color=color_list[c_idx], linewidth=4, zorder=8)

axA.set_aspect('equal')
axA.set_title(f"Aromatic Flake Bond Densities ($N={N}$)\nThickness ∝ ρᵢⱼ", fontsize=12)
axA.axis('off')
axA.legend(loc='lower left', fontsize=7)

# --- PLOT B: Scatter - Method vs Exact ---
axB = plt.subplot(2, 3, 2)
axB.scatter(rho_exact_vals, rho_cheb, c='red', alpha=0.5, s=20, label='Chebyshev')
axB.scatter(rho_exact_vals, rho_cf, c='blue', alpha=0.5, s=20, label='Lanczos/CF')
axB.plot([0, 1.5], [0, 1.5], 'k--', lw=1, label='Perfect match')
axB.set_xlabel('Exact $\\rho_{ij}$')
axB.set_ylabel('Estimated $\\rho_{ij}$')
axB.set_title('Accuracy: Estimated vs Exact')
axB.legend(fontsize=8)
axB.grid(True, alpha=0.3)
axB.set_xlim(0.3, 1.2)
axB.set_ylim(0.3, 1.2)
axB.set_aspect('equal', adjustable='datalim')

# --- PLOT C: Error vs bond index ---
axC = plt.subplot(2, 3, 3)
axC.plot(err_cheb, 'r-', alpha=0.7, label=f'Chebyshev (mean={err_cheb.mean():.5f})')
axC.plot(err_cf, 'b-', alpha=0.7, label=f'Lanczos/CF (mean={err_cf.mean():.5f})')
axC.set_xlabel('Bond index')
axC.set_ylabel('Absolute error |ρ_est - ρ_exact|')
axC.set_title('Per-bond error comparison')
axC.legend(fontsize=8)
axC.grid(True, alpha=0.3)

# --- PLOT D: Convergence - Chebyshev ---
axD = plt.subplot(2, 3, 4)
for c_idx, (name, (idx, i, j)) in enumerate(tracked_bonds_filtered.items()):
    axD.plot(M_range, cheb_conv[name], marker='o', color=color_list[c_idx],
             linewidth=1.5, label=name, markersize=4)
    axD.axhline(exact_vals[name], color=color_list[c_idx], linestyle='--', alpha=0.4)
axD.set_title('Chebyshev: Convergence vs M')
axD.set_xlabel('Polynomial order M')
axD.set_ylabel('Estimated ρᵢⱼ')
axD.grid(True, alpha=0.3)
axD.legend(fontsize=7)
axD.set_xticks(range(0, M_max_conv + 1, 5))

# --- PLOT E: Convergence - Lanczos/CF ---
axE = plt.subplot(2, 3, 5)
for c_idx, (name, (idx, i, j)) in enumerate(tracked_bonds_filtered.items()):
    axE.plot(M_range, cf_conv[name], marker='s', color=color_list[c_idx],
             linewidth=1.5, label=name, markersize=4)
    axE.axhline(exact_vals[name], color=color_list[c_idx], linestyle='--', alpha=0.4)
axE.set_title('Lanczos/CF: Convergence vs M')
axE.set_xlabel('Recursion depth M')
axE.set_ylabel('Estimated ρᵢⱼ')
axE.grid(True, alpha=0.3)
axE.legend(fontsize=7)
axE.set_xticks(range(0, M_max_conv + 1, 5))

# --- PLOT F: Convergence rate (error vs M, log scale) ---
axF = plt.subplot(2, 3, 6)
for c_idx, (name, (idx, i, j)) in enumerate(tracked_bonds_filtered.items()):
    err_cheb_M = np.abs(cheb_conv[name] - exact_vals[name])
    err_cf_M = np.abs(cf_conv[name] - exact_vals[name])
    axF.semilogy(M_range, err_cheb_M, 'o-', color=color_list[c_idx],
                 linewidth=1.5, label=f'Cheb: {name}', alpha=0.6, markersize=3)
    axF.semilogy(M_range, err_cf_M, 's--', color=color_list[c_idx],
                 linewidth=1.5, label=f'CF: {name}', alpha=0.9, markersize=3)
axF.set_title('Convergence rate: Error vs M (log scale)')
axF.set_xlabel('M')
axF.set_ylabel('Absolute error (log scale)')
axF.grid(True, alpha=0.3, which='both')
axF.legend(fontsize=6, ncol=2)
axF.set_xticks(range(0, M_max_conv + 1, 5))

plt.tight_layout()
plt.savefig('KekuleOrderN_BOP_2D_comparison.png', dpi=150)
plt.show()

# =====================================================================
# 10. SUPPLEMENTARY: Quick 1D chain comparison (Peierls distortion)
#     for clean comparison between Chebyshev and Lanczos/CF
# =====================================================================
print("\n\n========== SUPPLEMENTARY: 1D Peierls chain comparison ==========")
N_1d = 40
t0_1d = -2.5
u_peierls = 0.05
alpha_1d = 4.1

H_1d = np.zeros((N_1d, N_1d))
for i in range(N_1d - 1):
    t = t0_1d - alpha_1d * (2 * u_peierls if i % 2 == 0 else -2 * u_peierls)
    H_1d[i, i + 1] = H_1d[i + 1, i] = t
H_1d[N_1d // 2, N_1d // 2] = e_nitrogen  # Nitrogen defect at center

# Exact
evals_1d, evecs_1d = np.linalg.eigh(H_1d)
rho_exact_1d = np.zeros((N_1d, N_1d))
for lam in range(N_1d):
    if evals_1d[lam] < 0.0:
        rho_exact_1d += 2.0 * np.outer(evecs_1d[:, lam], evecs_1d[:, lam])

# Scale for Chebyshev
E_max_1d = 2.0 * np.max(np.abs(H_1d)) + 0.5
H_1d_scaled = H_1d / E_max_1d

# Scan M from 1 to 15 and compute average error across all bonds
M_scan = np.arange(1, 16)
err_cheb_1d = []
err_cf_1d = []

for M in M_scan:
    err_cheb_bonds = []
    err_cf_bonds = []
    for i in range(N_1d - 1):
        rho_c = get_bond_density_cheb(H_1d_scaled, i, i + 1, M)
        rho_l = get_bond_density_CF(H_1d, i, i + 1, M, eta=0.15, n_points=300)
        err_cheb_bonds.append(abs(rho_c - rho_exact_1d[i, i + 1]))
        err_cf_bonds.append(abs(rho_l - rho_exact_1d[i, i + 1]))
    err_cheb_1d.append(np.mean(err_cheb_bonds))
    err_cf_1d.append(np.mean(err_cf_bonds))

fig2, (ax2a, ax2b) = plt.subplots(1, 2, figsize=(12, 5))

# 1D bond density profile
bonds_1d = np.arange(N_1d - 1)
exact_bonds_1d = np.array([rho_exact_1d[i, i + 1] for i in range(N_1d - 1)])
cheb_bonds_1d = np.array([get_bond_density_cheb(H_1d_scaled, i, i + 1, 15) for i in range(N_1d - 1)])
cf_bonds_1d = np.array([get_bond_density_CF(H_1d, i, i + 1, 15, eta=0.15, n_points=300) for i in range(N_1d - 1)])

ax2a.plot(bonds_1d, exact_bonds_1d, 'k-', lw=2, label='Exact')
ax2a.plot(bonds_1d, cheb_bonds_1d, 'r--', lw=1.5, label=f'Chebyshev (M=15)')
ax2a.plot(bonds_1d, cf_bonds_1d, 'b:', lw=1.5, label=f'Lanczos/CF (M=15)')
ax2a.axvline(N_1d // 2, color='gray', ls=':', label='N defect')
ax2a.set_xlabel('Bond index')
ax2a.set_ylabel('ρᵢ,ᵢ₊₁')
ax2a.set_title('1D Peierls chain: bond density profile')
ax2a.legend(fontsize=8)
ax2a.grid(True, alpha=0.3)

# 1D convergence
ax2b.semilogy(M_scan, err_cheb_1d, 'ro-', lw=1.5, label='Chebyshev')
ax2b.semilogy(M_scan, err_cf_1d, 'bs--', lw=1.5, label='Lanczos/CF')
ax2b.set_xlabel('M')
ax2b.set_ylabel('Mean absolute error (log scale)')
ax2b.set_title('1D chain: Mean error convergence vs M')
ax2b.grid(True, alpha=0.3, which='both')
ax2b.legend(fontsize=9)
ax2b.set_xticks(M_scan)

plt.tight_layout()
plt.savefig('KekuleOrderN_1D_comparison.png', dpi=150)
plt.show()
