import numpy as np
import matplotlib.pyplot as plt

# ==========================================
# 1. Math Helper Functions
# ==========================================

def get_smootherstep(r, r_min, r_max):
    """
    Computes blending factor S(r) and its derivative dS/dr.
    S goes from 0 (at r_min) to 1 (at r_max).
    """
    if r <= r_min:
        return 0.0, 0.0
    elif r >= r_max:
        return 1.0, 0.0
    
    # Normalized coordinate x in [0, 1]
    x = (r - r_min) / (r_max - r_min)
    
    # Smootherstep: 6x^5 - 15x^4 + 10x^3
    # This has 0 derivative at x=0 and x=1 (C2 continuous)
    S = x*x*x * (x * (x * 6 - 15) + 10)
    
    # Derivative dS/dx = 30x^4 - 60x^3 + 30x^2
    dSdx = 30*x*x*x*x - 60*x*x*x + 30*x*x
    
    # dS/dr = dS/dx * dx/dr
    dSdr = dSdx / (r_max - r_min)
    
    return S, dSdr

def vector_norm(v):
    return np.sqrt(np.sum(v**2))

# ==========================================
# 2. Physics & Logic
# ==========================================

class Cluster:
    def __init__(self, N, radius, center_pos, total_charge_zero=False):
        self.N = N
        # Random positions inside a sphere
        pos = np.random.rand(N, 3) * 2 - 1
        # Filter to sphere approx
        norms = np.linalg.norm(pos, axis=1)
        pos = pos / norms[:,None] * (np.random.rand(N,1)**(1/3)) * radius
        
        # Charges
        q = np.random.rand(N) * 2 - 1
        if total_charge_zero:
            q -= np.mean(q) # Make exactly neutral
        
        self.q_local = q
        self.pos_local = pos # Relative to cluster center
        self.center = np.array(center_pos, dtype=float)
        
        # Precompute Multipoles
        self.Q = np.sum(q)
        self.dipole = np.sum(q[:, None] * pos, axis=0) # Dipole relative to center

    def get_global_pos(self):
        return self.pos_local + self.center

def compute_exact(c1, c2):
    """
    Brute force N*M Coulomb interaction.
    Returns: Total Energy, Forces on C1 particles (Nx3)
    """
    pos1 = c1.get_global_pos()
    pos2 = c2.get_global_pos()
    
    E_total = 0.0
    forces1 = np.zeros_like(pos1)
    
    for i in range(c1.N):
        for j in range(c2.N):
            rij_vec = pos1[i] - pos2[j]
            r = np.linalg.norm(rij_vec)
            if r < 1e-9: continue
            
            # Coulomb Energy: q1*q2 / r
            E_total += c1.q_local[i] * c2.q_local[j] / r
            
            # Force on 1 due to 2: q1*q2 * vec / r^3
            # (Force pushes i away from j)
            f_vec = (c1.q_local[i] * c2.q_local[j] / (r**3)) * rij_vec
            forces1[i] += f_vec
            
    return E_total, forces1

def compute_multipole_expansion(c1, c2):
    """
    Computes:
    1. Approx Energy (Charge+Dipole)
    2. Force on Center of C1 (due to C2)
    3. Electric Field at Center of C1 (due to C2) acting on dipoles
    """
    R_vec = c1.center - c2.center # Vector from 2 to 1
    R = np.linalg.norm(R_vec)
    n = R_vec / R # Unit vector pointing 2 -> 1
    
    # ---------------------------
    # Energy Calculation
    # ---------------------------
    # T0: Charge-Charge
    E_qq = (c1.Q * c2.Q) / R
    
    # T1: Charge-Dipole
    # Pot from 2 at 1: V2 = (mu2 . n) / R^2 + Q2/R
    # Energy = Q1 * V2 + Q2 * V1 (careful with signs for V1)
    # Pot from 1 at 2: V1 = (mu1 . (-n)) / R^2
    term_c1_d2 = (c1.Q * np.dot(c2.dipole, n)) / (R**2) # Q1 seeing dipole 2
    term_c2_d1 = (c2.Q * np.dot(c1.dipole, -n)) / (R**2) # Q2 seeing dipole 1
    E_qd = term_c1_d2 + term_c2_d1
    
    # T2: Dipole-Dipole
    mu1_dot_mu2 = np.dot(c1.dipole, c2.dipole)
    mu1_dot_n   = np.dot(c1.dipole, n)
    mu2_dot_n   = np.dot(c2.dipole, n)
    
    E_dd = (mu1_dot_mu2 - 3 * mu1_dot_n * mu2_dot_n) / (R**3)
    
    E_approx = E_qq + E_qd + E_dd
    
    # ---------------------------
    # Force Calculation (Force on C1 CoM)
    # F = - Gradient_R (E)
    # ---------------------------
    
    # 1. Deriv of QQ: -(-1/R^2) * Q Q * n = Q Q / R^2 * n
    F_qq = (c1.Q * c2.Q / (R**2)) * n
    
    # 2. Deriv of QD
    # d/dr (1/r^2) = -2/r^3
    # d/dr (n) = (I - n*n)/r
    # This gets messy. Standard formula for Field from dipole 2 at dist R:
    # E_field_from_2 = (3(p.n)n - p) / R^3
    # Force on Q1 = Q1 * E_field_from_2
    field_from_2 = (3 * np.dot(c2.dipole, n) * n - c2.dipole) / (R**3)
    F_on_Q1 = c1.Q * field_from_2
    
    # Force on Dipole 1 due to Charge 2 (Field of monopole Q2)
    # Field from Q2 at R: (Q2 / R^2) * n
    # Gradient of (mu1 . E) ... Force on dipole = (mu . nabla) E
    # Result is: Q2 * (3(mu1.n)n - mu1) / R^3  (Symmetric to above)
    F_on_D1 = c2.Q * (3 * np.dot(c1.dipole, n) * n - c1.dipole) / (R**3)
    
    # 3. Deriv of DD
    # Force on dipole 1 due to dipole 2
    # F = 3/R^4 [ (m1.m2)n + (m1.n)m2 + (m2.n)m1 - 5(m1.n)(m2.n)n ]
    # Note: Signs depend on definition of n. Here n points 2->1 (repulsive direction)
    F_dd = (3 / R**4) * ( 
           mu1_dot_mu2 * n + 
           mu1_dot_n * c2.dipole + 
           mu2_dot_n * c1.dipole - 
           5 * mu1_dot_n * mu2_dot_n * n
           )

    F_com_on_1 = F_qq + F_on_Q1 + F_on_D1 + F_dd

    # ---------------------------
    # "Torque"/Field Terms (Internal distribution)
    # ---------------------------
    # We need the Field (E_loc) acting on C1 to distribute to atoms q_i
    # E_loc = - Gradient_mu1 (E)
    # Terms containing mu1:
    # 1. Q2 * (mu1 . -n) / R^2  -> Grad = -Q2/R^2 * n
    # 2. (mu1.mu2 - 3(mu1.n)(mu2.n)) / R^3
    
    term1 = - (c2.Q / (R**2)) * n
    term2 = (c2.dipole - 3 * mu2_dot_n * n) / (R**3)
    
    # E_field_at_1 is the external field felt by cluster 1
    # Note: The sign is flipped because E = -Grad V. 
    # The Energy term was - mu1 . E_ext.
    # So E_ext = - Gradient_mu1 ( Energy ) * -1 ? No.
    # U = - mu . E.  => dU/dmu = -E. => E = - dU/dmu.
    E_field_at_1 = -(term1 + term2)
    
    return E_approx, F_com_on_1, E_field_at_1

def compute_tiled_fmm(c1, c2, R_min, R_max):
    """
    Computes smoothed forces on C1 particles.
    """
    # 1. Calc Distance
    R_vec = c1.center - c2.center
    R = np.linalg.norm(R_vec)
    n = R_vec / R
    
    # 2. Get Blending Factors
    S, dSdr = get_smootherstep(R, R_min, R_max)
    
    # 3. Compute Exact (Always needed if S < 1 or for transition)
    # Optimization: if S=1 (far), we don't strictly need exact, 
    # but to plot the comparison we calculate it. 
    # In real code: if S==1: E_ex=0, F_ex=0.
    E_ex, F_ex = compute_exact(c1, c2)
    
    # 4. Compute Multipole
    E_app, F_com, E_field = compute_multipole_expansion(c1, c2)
    
    # 5. Blend
    # E_total = (1-S)E_ex + S*E_app
    E_total = (1-S)*E_ex + S*E_app
    
    forces_1 = np.zeros_like(F_ex)
    
    # Force Correction Term (Energy Barrier)
    # F_corr = - Gradient(S) * (E_app - E_ex)
    # Grad S = dS/dR * n (since n points radial outward from 2->1)
    # Note: Forces are on particles of 1. Moving 1 away increases R.
    # So dR/dr_i = n * (1/N) ? No, depends on mass. Assume equal mass -> 1/N.
    # Actually dR/dr_i = n * (1/N_1)
    
    dE_diff = (E_app - E_ex)
    F_switch_mag = - dSdr * dE_diff
    F_switch_vec = F_switch_mag * n * (1.0 / c1.N) # Distribute to all i
    
    for i in range(c1.N):
        # A. Exact Component
        f_exact_i = (1.0 - S) * F_ex[i]
        
        # B. Multipole Component
        # 1. CoM force distributed by mass (1/N)
        # 2. Local Field acting on charge q_i
        # 3. (Optional) Centrifugal terms for torque - ignored for point particles here
        f_multi_i = S * ( (F_com / c1.N) + (c1.q_local[i] * E_field) )
        
        # C. Switching Component
        f_switch_i = F_switch_vec
        
        forces_1[i] = f_exact_i + f_multi_i + f_switch_i
        
    return E_total, forces_1, E_ex, E_app, F_ex

# ==========================================
# 3. Simulation & Plotting
# ==========================================

def run_simulation():
    # Parameters
    N1, N2 = 5, 5
    Rad1, Rad2 = 2.0, 2.0
    R_start = 0.5  # Overlap!
    R_end = 25.0
    Steps = 100
    
    R_min_sw = 8.0
    R_max_sw = 12.0
    
    # Create Clusters
    c1 = Cluster(N1, Rad1, [0,0,0], total_charge_zero=False)
    c2 = Cluster(N2, Rad2, [R_end,0,0], total_charge_zero=False) # Start far
    
    distances = np.linspace(R_start, R_end, Steps)
    
    # History
    hist_R = []
    hist_E_ex = []
    hist_E_fmm = []
    hist_F_ex = []  # List of arrays [N, 3]
    hist_F_fmm = [] # List of arrays [N, 3]
    
    print(f"Cluster 1: Q={c1.Q:.2f}, mu={vector_norm(c1.dipole):.2f}")
    print(f"Cluster 2: Q={c2.Q:.2f}, mu={vector_norm(c2.dipole):.2f}")
    print(f"Switching Zone: {R_min_sw} -> {R_max_sw}")
    
    for r in distances:
        # Move Cluster 2
        c2.center = np.array([r, 0.0, 0.0])
        
        # Compute
        E_fmm, F_fmm, E_ex_val, E_app_val, F_ex_val = compute_tiled_fmm(c1, c2, R_min_sw, R_max_sw)
        
        # Store
        hist_R.append(r)
        hist_E_ex.append(E_ex_val)
        hist_E_fmm.append(E_fmm)
        hist_F_ex.append(F_ex_val) # Store force on all particles
        hist_F_fmm.append(F_fmm)

    # Convert to numpy for plotting
    hist_R = np.array(hist_R)
    hist_E_ex = np.array(hist_E_ex)
    hist_E_fmm = np.array(hist_E_fmm)
    hist_F_ex = np.array(hist_F_ex) # Shape [Steps, N1, 3]
    hist_F_fmm = np.array(hist_F_fmm)
    
    # ==========================
    # Plotting
    # ==========================
    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    
    # 1. Energy Plot
    ax = axes[0]
    ax.plot(hist_R, hist_E_ex, 'k:', lw=2, label='Exact Reference')
    ax.plot(hist_R, hist_E_fmm, 'r-', lw=1, label='Cluster FMM (Blended)')
    
    # Draw switching region
    ax.axvspan(R_min_sw, R_max_sw, color='yellow', alpha=0.2, label='Switching Zone')
    
    ax.set_ylabel('Total Potential Energy')
    ax.set_title('Energy Smoothness Check')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 2. Force Plot (Fx component)
    ax = axes[1]
    
    colors = plt.cm.tab10(np.linspace(0, 1, N1))
    
    method_lines = []
    # Proxy handles to explain ls/lw meaning (keep colors per particle)
    method_lines.append(plt.Line2D([0], [0], ls=':', lw=1.5, color='k', label='Exact'))
    method_lines.append(plt.Line2D([0], [0], ls='-', lw=0.5, color='k', label='Cluster FMM'))
    
    for i in range(N1):
        # Extract Fx for particle i across all steps
        fx_ex = hist_F_ex[:, i, 0]
        fx_fmm = hist_F_fmm[:, i, 0]
        
        c = colors[i]
        ax.plot(hist_R, fx_ex, ls=':', lw=2.0, color=c, alpha=0.9)
        ax.plot(hist_R, fx_fmm, ls='-', lw=1.5, color=c, label=f'Part {i}')
        
    ax.axvspan(R_min_sw, R_max_sw, color='yellow', alpha=0.2)
    ax.set_ylabel('Force X (on C1 particles)')
    ax.set_xlabel('Distance R')
    ax.set_title('Force Smoothness Check (Per Particle)')
    
    # Add legend: particles (colors) + methods (ls/lw)
    particle_legend = ax.legend(loc='upper right', title='Particles')
    ax.add_artist(particle_legend)
    ax.legend(handles=method_lines, loc='lower right', title='Methods')
    
    # Zoom Y axis to avoid divergence at R=0 ruining the scale
    # Find min/max in the valid range (e.g., R > 3)
    valid_idx = hist_R > 3.0
    if np.any(valid_idx):
        ymin = np.min(hist_F_ex[valid_idx, :, 0])
        ymax = np.max(hist_F_ex[valid_idx, :, 0])
        margin = (ymax - ymin) * 0.1
        ax.set_ylim(ymin - margin, ymax + margin)

    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    run_simulation()