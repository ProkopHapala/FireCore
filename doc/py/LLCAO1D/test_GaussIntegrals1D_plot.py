# /home/prokop/git/FireCore/doc/py/LLCAO1D/test_integrals_vs_distance.py

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf # For the analytical erf function in Coulomb interaction

from quantum_solver_1D import (
    gauss_to_grid,
    boys_function_F0,
    overlap_element_gauss,
    kinetic_element_gauss,
    coulomb_element_gauss
)

# ============== Gaussian Basis Functions and Integrals ================

def d2_gauss_dx2_grid(x_vals, Xc, a):
    """
    Computes -1/2 * d^2/dx^2 phi(x) for a normalized Gaussian phi(x) on a grid.
    """
    N = (2 * a / np.pi)**0.25
    term_par = a - 2 * a**2 * (x_vals - Xc)**2
    exp_val = np.exp(-a * (x_vals - Xc)**2)
    return N * term_par * exp_val

def num_d2_gauss_dx2_grid(x_vals, Xc, a):
    """
    Computes -1/2 * d^2/dx^2 phi(x) for a Gaussian phi(x) on a grid,
    where the Gaussian is first projected to the grid, and then derivatives are numerical.
    """
    phi_grid = gauss_to_grid(x_vals, Xc, a)
    d_phi_dx = np.gradient(phi_grid, x_vals, edge_order=2)
    d2_phi_dx2 = np.gradient(d_phi_dx, x_vals, edge_order=2)
    return -0.5 * d2_phi_dx2

def product_gauss_grid(x_vals, a1, X1, a2, X2):
    """
    Computes the product of two normalized Gaussians phi1*phi2 on a grid.
    """
    N1 = (2 * a1 / np.pi)**0.25
    N2 = (2 * a2 / np.pi)**0.25
    
    p = a1 + a2
    if abs(p) < 1e-12:
        Xp = (X1 + X2) / 2.0
        K_exp_term = 0
    else:
        Xp = (a1 * X1 + a2 * X2) / p
        K_exp_term = np.exp(-(a1 * a2 / p) * (X1 - X2)**2)

    K = N1 * N2 * K_exp_term
    return K * np.exp(-p * (x_vals - Xp)**2)

def num_product_gauss_grid(x_vals, a1, X1, a2, X2):
    """
    Computes the product of two normalized Gaussians phi1*phi2 on a grid
    by first projecting each Gaussian to the grid and then multiplying their grid representations.
    """
    phi1_grid = gauss_to_grid(x_vals, X1, a1)
    phi2_grid = gauss_to_grid(x_vals, X2, a2)
    return phi1_grid * phi2_grid

def smeared_coulomb_pot_grid(x_vals, X_nuc, Z_nuc, a_nuc_smear):
    """
    Computes potential V(x) from a Gaussian nuclear charge distribution.
    """
    if a_nuc_smear <= 1e-10:
        print(f"Warning: a_nuc_smear ({a_nuc_smear}) is very small. Potential may be unstable.")
        return -Z_nuc / (np.abs(x_vals - X_nuc) + 1e-7)

    t_arg = a_nuc_smear * (x_vals - X_nuc)**2
    f0_vals = np.array([boys_function_F0(t_val) for t_val in t_arg])
    return -Z_nuc * np.sqrt(2 * np.pi / a_nuc_smear) * f0_vals

def ana_gauss_charge_interact(Q1, a1, X1, Q2, a2, X2):
    """
    Analytical interaction energy between two Gaussian charge distributions.
    """
    R = np.abs(X1 - X2)
    if a1 <= 1e-10 or a2 <= 1e-10:
        print(f"Warning: Very small alpha ({a1}, {a2}) in Gaussian interaction.")
        return 0.0
    if R < 1e-9:
        return Q1 * Q2 * 2 * np.sqrt(a1 * a2 / (np.pi * (a1 + a2)))
    s2_sum = 0.5/a1 + 0.5/a2
    arg_erf = R / np.sqrt(s2_sum)
    return (Q1 * Q2 / R) * erf(arg_erf)

# --- Core numerical integration and plotting functions ---

def num_int_vs_dist(x_vals, f1_grid, f2_grid, max_sh_idx):
    """
    Computes Integral(func1(x) * func2(x-d))dx for various distances d.
    """
    dx = x_vals[1] - x_vals[0]
    nx = len(x_vals)
    
    if max_sh_idx >= nx // 2:
        print("Warning: max_shift_idx is large relative to grid size. Periodic boundary effects might be significant.")

    dists = np.arange(0, max_sh_idx + 1) * dx
    num_integrals = np.zeros_like(dists, dtype=float)

    for i, sh_idx in enumerate(range(max_sh_idx + 1)):
        f2_roll = np.roll(f2_grid, sh_idx)
        integrand = f1_grid * f2_roll
        num_integrals[i] = np.trapz(integrand, dx=dx)
        
    return dists, num_integrals

# ============== Plotting Functions ================

def plot1d(x_data, plots_data, title="", xlabel="x", ylabel="y", fig=None, ax=None):
    """
    Plots multiple functions on the same axes.

    Args:
        x_data (np.ndarray): The x-axis data.
        plots_data (list): A list of tuples, where each tuple contains
                           (y_data, label, format_string, linewidth).
                           Example: [(y1, 'Data 1', 'b-', 1.0), (y2, 'Data 2', 'r:', 1.5)]
        title (str): Plot title.
        xlabel (str): X-axis label.
        ylabel (str): Y-axis label.
        fig (matplotlib.figure.Figure, optional): Existing figure.
        ax (matplotlib.axes.Axes, optional): Existing axes.
    """
    show_plot_after = False
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 4)) # Default size
        show_plot_after = True

    for y_vals, lbl, fmt_str, lw in plots_data:
        ax.plot(x_data, y_vals, fmt_str, lw=lw, label=lbl)
    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend()
    ax.grid(True)
    
    if show_plot_after:
        plt.show()
    return fig, ax

def plot_integral_vs_distance(
        d_vals, ana_vals, num_vals,
        x_plot, f1_grid, f2_init_grid,
        title, int_lbl="Integral Value", #NOSONAR
        f1_lbl="Stationary Func1",
        f2_lbl="Initial Rolled Func2",
        show_diff=False, show_ratio=False, diff_scale=1.0, eps_ratio=1e-9
    ):
    """
    Plots analytical vs numerical integrals against distance, and the initial functions.
    Optionally plots difference and ratio.
    """
    num_subplots = 2
    if show_diff:
        num_subplots += 1
    if show_ratio:
        num_subplots += 1
    fig, axes = plt.subplots(num_subplots, 1, figsize=(12, 4 * num_subplots), sharex=False)
    fig.suptitle(title, fontsize=16)

    current_ax_idx = 0
    plot1d(d_vals, [
        (ana_vals, f'Analytical {int_lbl}', 'b-', 0.5),
        (num_vals, f'Numerical (roll) {int_lbl}', 'r:', 1.5)],
    title='Integral vs. Distance', xlabel='Distance d (Bohr)', ylabel=int_lbl, fig=fig, ax=axes[current_ax_idx])
    current_ax_idx+=1
    
    plot1d(x_plot,[
        (f1_grid,      f1_lbl, 'g-', 0.7),
        (f2_init_grid, f2_lbl, 'm:', 1.0)],
    title='Initial Functions on Grid', xlabel='x (Bohr)', ylabel='Amplitude', fig=fig, ax=axes[current_ax_idx])
    current_ax_idx+=1
    if show_diff:
        diff = (ana_vals - num_vals) * diff_scale
        plot1d(d_vals, [(diff, f'Difference (Ana - Num) * {diff_scale}', 'k-', 1.0)],title='Difference Plot', xlabel='Distance d (Bohr)', ylabel='Difference',fig=fig, ax=axes[current_ax_idx])
        current_ax_idx+=1
    if show_ratio:
        # Avoid division by zero or very small numbers for ratio
        safe_num_vals = num_vals + np.sign(num_vals)*eps_ratio + eps_ratio*(num_vals==0)
        ratio = ana_vals / safe_num_vals
        plot1d(d_vals, [(ratio, 'Ratio (Ana / Num)', 'purple', 1.0)], title='Ratio Plot', xlabel='Distance d (Bohr)', ylabel='Ratio',fig=fig, ax=axes[current_ax_idx])
        axes[current_ax_idx].axhline(1.0, color='gray', linestyle=':', linewidth=0.8)

    plt.tight_layout(rect=[0, 0, 1, 0.96]) # Adjust layout to make space for suptitle
    plt.show()


# ============== Test Functions ================

def test_overlap():
    # --- 1. Overlap Integral Test: <phi1(x) | phi2(x-d)> ---
    # phi1 is stationary (func1), phi2 is rolled (func2)
    print("\n--- Testing Overlap Integral S(d) ---")
    phi1_s_on_grid = gauss_to_grid(x_grid, X_center_stationary, alpha_1) # Stationary
    phi2_s_grid = gauss_to_grid(x_grid, X_center_stationary, alpha_2) # This will be rolled

    d_s, num_s_vals = num_int_vs_dist(x_grid, phi1_s_on_grid, phi2_s_grid, max_shift_idx)
    ana_s_vals = np.array([overlap_element_gauss(alpha_1, X_center_stationary, alpha_2, X_center_stationary + d) for d in d_s])
    
    plot_integral_vs_distance(
        d_s, ana_s_vals, num_s_vals,
        x_grid, phi1_s_on_grid, phi2_s_grid,
        title=f'Overlap Integral S(d) for $\\alpha_1=${alpha_1}, $\\alpha_2=${alpha_2}',
        int_lbl='S(d)',
        f1_lbl=f'Stationary $\\phi_1(x)$ ($\\alpha_1=${alpha_1})',
        f2_lbl=f'Initial $\\phi_2(x)$ ($\\alpha_2=${alpha_2})')

def test_kinetic_energy():
    # --- 2. Kinetic Energy Integral Test: <phi1(x-d) | -1/2 d^2/dx^2 | phi2(x)> ---
    print("\n--- Testing Kinetic Energy Integral T(d) ---")
    # func1: stationary, -1/2 d^2/dx^2 phi_stat(x)
    d2phi_stat_ana_grid = d2_gauss_dx2_grid(x_grid, X_center_stationary, alpha_1)
    d2phi_stat_num_grid = num_d2_gauss_dx2_grid(x_grid, X_center_stationary, alpha_1)

    # Plot comparison of analytical vs numerical second derivative
    plot1d(x_grid, [
        (d2phi_stat_ana_grid, r'Analytical $-1/2 d^2\phi/dx^2$ then grid', 'b-', 0.5),
        (d2phi_stat_num_grid, r'Grid then numerical $-1/2 d^2\phi/dx^2$', 'r:', 1.5)
        ],
        title=rf'Comparison of $-1/2 d^2\phi/dx^2$ for $\alpha={alpha_1}$',
        xlabel='x (Bohr)', ylabel='Amplitude',
    )

    # func2: rolled, phi_roll(x)
    phi_roll_t_grid = gauss_to_grid(x_grid, X_center_stationary, alpha_2)

    d_t, num_t_vals = num_int_vs_dist(x_grid, d2phi_stat_ana_grid, phi_roll_t_grid, max_shift_idx)

    ana_t_vals = np.zeros_like(d_t, dtype=float)
    for i, d_val in enumerate(d_t):
        s_overlap = overlap_element_gauss(alpha_2, X_center_stationary + d_val, alpha_1, X_center_stationary)
        ana_t_vals[i] = kinetic_element_gauss(alpha_2, X_center_stationary + d_val, alpha_1, X_center_stationary, s_overlap)
    
    plot_integral_vs_distance(
        d_t, ana_t_vals, num_t_vals,
        x_grid, d2phi_stat_ana_grid, phi_roll_t_grid,
        title=rf'Kinetic Integral T(d) for $\alpha_{{stat}}={alpha_1}$, $\alpha_{{roll}}={alpha_2}$',
        int_lbl='T(d)',
        f1_lbl=rf'Stationary $-1/2 d^2/dx^2 \phi_{{stat}}(x)$ ($\alpha_1={alpha_1}$)',
        f2_lbl=rf'Initial $\phi_{{roll}}(x)$ ($\alpha_2={alpha_2}$)')

def test_coulomb_integral():
    print("\n--- Testing Coulomb Integral V(d) ---")
    Z_nuc = 1.0
    X_nuc_stat = 0.0 
    a_nuc_smear  = 4.2 

    # func1: stationary, potential from a smeared nucleus
    V_smear_grid = smeared_coulomb_pot_grid(x_grid, X_nuc_stat, Z_nuc, a_nuc_smear)
    
    # func2: rolled, product Gaussian phi_P(x) formed from two Gaussians (alpha_p1, alpha_p2) centered at X_center_stationary
    ap1, ap2 = 0.6, 0.4
    p_el_prod = ap1 + ap2

    # Product Gaussian initially centered at X_center_stationary (0.0)
    phi_P_ana_grid = product_gauss_grid(x_grid, ap1, X_center_stationary, ap2, X_center_stationary)
    phi_P_num_grid = num_product_gauss_grid(x_grid, ap1, X_center_stationary, ap2, X_center_stationary)

    # Plot comparison of analytical vs numerical product Gaussian
    plot1d(x_grid, [
        (phi_P_ana_grid, r'Analytical Product Gaussian $\phi_P(x)$', 'b-', 0.5),
        (phi_P_num_grid, r'Numerical Product Gaussian (grid1 * grid2)', 'r:', 1.5)],
        title=rf'Comparison of Product Gaussian for $\alpha_{{p1}}={ap1}$, $\alpha_{{p2}}={ap2}$',
        xlabel='x (Bohr)', ylabel='Amplitude',
    )
    
    # For the integral test, func2 is the analytical product Gaussian
    phi_P_int_grid = phi_P_ana_grid 

    d_v, num_v_vals = num_int_vs_dist(x_grid, V_smear_grid, phi_P_int_grid, max_shift_idx)

    # Analytical calculation for interaction of smeared nucleus with electron product Gaussian
    ana_v_vals = np.zeros_like(d_v, dtype=float)
    Q_el_prod = overlap_element_gauss(ap1, X_center_stationary, ap2, X_center_stationary)
    Q_nuc_smear = -Z_nuc * np.sqrt(np.pi / a_nuc_smear)

    for i, d_val in enumerate(d_v):
        X_el_prod_c = X_center_stationary + d_val
        ana_v_vals[i] = ana_gauss_charge_interact(
            Q1=Q_nuc_smear, a1=a_nuc_smear, X1=X_nuc_stat,
            Q2=Q_el_prod,   a2=p_el_prod,   X2=X_el_prod_c
        )

    print( "V_smear_grid",   V_smear_grid  [:10] )   
    print( "phi_P_int_grid", phi_P_int_grid[:10] )                  

    plot_integral_vs_distance(
        d_v, ana_v_vals, num_v_vals,
        x_grid, V_smear_grid, phi_P_int_grid,
        title=rf'Coulomb Integral V(d) for Smeared Nucleus ($\alpha_{{nuc}}={a_nuc_smear}$) \& Product Gaussian ($\alpha_{{p1}}={ap1}$, $\alpha_{{p2}}={ap2}$)',
        int_lbl='V(d)', show_diff=True, show_ratio=True, diff_scale=1.0,
        f1_lbl=rf'Smeared Potential ($\alpha_{{nuc}}={a_nuc_smear}$, Z={Z_nuc})',
        f2_lbl=r'Initial Product Gaussian $\phi_P(x)$'
    )

if __name__ == '__main__':
    # --- Grid and General Parameters ---
    N_grid = 3001 
    hw_grid = 20.0 
    x_grid = np.linspace(-hw_grid, hw_grid, N_grid)
    dx_grid = x_grid[1] - x_grid[0]
    
    max_d_plot = 8.0 # Bohr
    max_shift_idx = int(max_d_plot / dx_grid)

    print(f"Grid: {N_grid} points from {-hw_grid} to {hw_grid} Bohr (step: {dx_grid:.4f})")
    print(f"Max shift index: {max_shift_idx} (corresponds to distance {max_shift_idx*dx_grid:.2f} Bohr)")

    # Basis function parameters (can be varied for different tests)
    alpha_1 = 0.8  # For Gaussian 1
    alpha_2 = 0.5  # For Gaussian 2
    X_center_stationary = 0.0 # Stationary function centered at 0

    test_overlap()
    test_kinetic_energy()
    test_coulomb_integral()
