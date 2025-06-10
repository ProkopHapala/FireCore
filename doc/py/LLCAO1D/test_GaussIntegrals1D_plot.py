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

# --- Plotting Helper ---
def plot_two_functions(x_data, y1_data, y2_data, title, xlabel, ylabel, y1_label, y2_label, fig=None, ax=None):
    """
    Plots two functions on the same axes.
    Uses thin solid blue line for y1_data and thicker dashed red line for y2_data.
    """
    show_plot_after = False
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 4)) # Default size if not part of a larger figure
        show_plot_after = True # Only call plt.show() if we created the figure here

    ax.plot(x_data, y1_data, 'b-', lw=0.5, label=y1_label)
    ax.plot(x_data, y2_data, 'r:', lw=1.5, label=y2_label)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend()
    ax.grid(True)
    if show_plot_after:
         plt.show()
    return fig, ax

# --- Helper functions specific to this test script ---
def d2_gaussian_dx2_on_grid(x_grid, X_center, alpha):
    """
    Computes -1/2 * d^2/dx^2 phi(x) for a normalized Gaussian phi(x) on a grid.
    phi(x) = N * exp(-alpha * (x-X_center)^2)
    d^2phi/dx^2 = N * (4*alpha^2*(x-X_center)^2 - 2*alpha) * exp(-alpha*(x-X_center)^2)
    So, -1/2 * d^2phi/dx^2 = N * (alpha - 2*alpha^2*(x-X_center)^2) * exp(-alpha*(x-X_center)^2)
    """
    N_norm = (2 * alpha / np.pi)**0.25
    term_in_parentheses = alpha - 2 * alpha**2 * (x_grid - X_center)**2
    exp_term = np.exp(-alpha * (x_grid - X_center)**2)
    return N_norm * term_in_parentheses * exp_term

def numerical_d2_gaussian_dx2_on_grid(x_grid, X_center, alpha):
    """
    Computes -1/2 * d^2/dx^2 phi(x) for a Gaussian phi(x) on a grid,
    where the Gaussian is first projected to the grid, and then derivatives are numerical.
    """
    phi_on_grid = gauss_to_grid(x_grid, X_center, alpha)
    # First derivative
    d_phi_dx = np.gradient(phi_on_grid, x_grid, edge_order=2)
    # Second derivative
    d2_phi_dx2 = np.gradient(d_phi_dx, x_grid, edge_order=2)
    return -0.5 * d2_phi_dx2

def get_product_gaussian_on_grid(x_grid, alpha1, X1, alpha2, X2):
    """
    Computes the product of two normalized Gaussians phi1*phi2 on a grid.
    phi1(x)*phi2(x) = K * exp(-p*(x-Xp)^2)
    where p = alpha1+alpha2, Xp = (a1X1+a2X2)/p
    K = N1*N2*exp(-(a1a2/p)*(X1-X2)^2)
    """
    N1 = (2 * alpha1 / np.pi)**0.25
    N2 = (2 * alpha2 / np.pi)**0.25
    
    p_prod = alpha1 + alpha2
    if abs(p_prod) < 1e-12: # Avoid division by zero if p is too small
        Xp_prod = (X1 + X2) / 2.0 # Or handle as an error/special case
        K_prod_exp_term = 0 # Effectively makes K_prod zero if exponents are tiny and opposite
    else:
        Xp_prod = (alpha1 * X1 + alpha2 * X2) / p_prod
        K_prod_exp_term = np.exp(-(alpha1 * alpha2 / p_prod) * (X1 - X2)**2)

    K_prod = N1 * N2 * K_prod_exp_term
    
    return K_prod * np.exp(-p_prod * (x_grid - Xp_prod)**2)

def numerical_product_gaussian_on_grid(x_grid, alpha1, X1, alpha2, X2):
    """
    Computes the product of two normalized Gaussians phi1*phi2 on a grid
    by first projecting each Gaussian to the grid and then multiplying their grid representations.
    """
    phi1_on_grid = gauss_to_grid(x_grid, X1, alpha1)
    phi2_on_grid = gauss_to_grid(x_grid, X2, alpha2)
    return phi1_on_grid * phi2_on_grid

def get_smeared_coulomb_potential_on_grid(x_grid, X_nuc, Z_nuc, alpha_nuc_smear):
    """
    Computes the potential V(x) = -Z_nuc * sqrt(2*pi/alpha_nuc_smear) * F0(alpha_nuc_smear * (x-X_nuc)^2)
    This is the potential at x due to a Gaussian nuclear charge distribution with exponent alpha_nuc_smear.
    """
    if alpha_nuc_smear <= 1e-10:
        print(f"Warning: alpha_nuc_smear ({alpha_nuc_smear}) is very small. Potential may be unstable or inaccurate.")
        # Fallback to a slightly regularized point charge potential if smearing is too weak
        epsilon_coulomb = 1e-7
        return -Z_nuc / (np.abs(x_grid - X_nuc) + epsilon_coulomb)

    t_boys_arg = alpha_nuc_smear * (x_grid - X_nuc)**2
    # Vectorize boys_function_F0 if it's not already (assuming it takes scalar t)
    f0_values = np.array([boys_function_F0(t_val) for t_val in t_boys_arg])
    
    return -Z_nuc * np.sqrt(2 * np.pi / alpha_nuc_smear) * f0_values

def analytical_gaussian_charges_interaction(Q1, alpha1, X1, Q2, alpha2, X2):
    """
    Analytical interaction energy between two Gaussian charge distributions.
    E = (Q1*Q2 / R) * erf(R / sqrt(1/(2*alpha1) + 1/(2*alpha2)))
    R = |X1 - X2|
    Assumes atomic units where 1/(4*pi*eps0) = 1.
    """
    R = np.abs(X1 - X2)

    if alpha1 <= 1e-10 or alpha2 <= 1e-10:
        # If one exponent is tiny, treat it as a point charge interacting with a Gaussian
        # This case needs careful handling or a different formula.
        # For simplicity, returning 0 or raising an error might be appropriate for extreme cases.
        print(f"Warning: Very small alpha ({alpha1}, {alpha2}) in Gaussian interaction.")
        return 0.0 # Or a more robust fallback

    if R < 1e-9: # Centers coincide
        # For R->0, (erf(x)/x) -> 2/sqrt(pi).
        # E = Q1*Q2 * (2/sqrt(pi)) / sqrt(1/(2*alpha1) + 1/(2*alpha2))
        return Q1 * Q2 * 2 * np.sqrt(alpha1 * alpha2 / (np.pi * (alpha1 + alpha2)))

    s_sq_sum = 0.5/alpha1 + 0.5/alpha2
    arg_erf = R / np.sqrt(s_sq_sum)
    return (Q1 * Q2 / R) * erf(arg_erf)

# --- Core numerical integration and plotting functions ---

def compute_numerical_integral_vs_distance(x_grid, func1_on_grid, func2_on_grid, max_shift_idx):
    """
    Computes Integral(func1(x) * func2(x-d))dx for various distances d.
    func1 is stationary. func2 is rolled.
    Distance d is achieved by rolling func2 by 'shift' indices.
    """
    dx_grid_step = x_grid[1] - x_grid[0]
    num_points = len(x_grid)
    
    if max_shift_idx >= num_points // 2:
        print("Warning: max_shift_idx is large relative to grid size. Periodic boundary effects might be significant.")

    distances = np.arange(0, max_shift_idx + 1) * dx_grid_step
    numerical_integrals = np.zeros_like(distances, dtype=float)

    for i, shift_idx in enumerate(range(max_shift_idx + 1)):
        func2_rolled = np.roll(func2_on_grid, shift_idx)
        integrand = func1_on_grid * func2_rolled
        numerical_integrals[i] = np.trapz(integrand, dx=dx_grid_step)
        
    return distances, numerical_integrals

def plot_integral_vs_distance(distances, analytical_values, numerical_values,
                              x_grid_plot, func1_on_grid, func2_initial_on_grid,
                              title, integral_label="Integral Value", #NOSONAR
                              func1_label="Stationary Func1",
                              func2_label="Initial Rolled Func2",
                              bDiff=False, bRatio=False, scDiff=1.0, epsilon_ratio=1e-9):
    """
    Plots analytical vs numerical integrals against distance, and the initial functions.
    Optionally plots difference and ratio.
    """
    num_subplots = 2
    if bDiff:
        num_subplots += 1
    if bRatio:
        num_subplots += 1
    fig, axes = plt.subplots(num_subplots, 1, figsize=(12, 4 * num_subplots), sharex=False)
    fig.suptitle(title, fontsize=16)


    # Plot 1: Integral values vs Distance
    current_ax_idx = 0
    axes[current_ax_idx].plot(distances, analytical_values, 'b-', lw=0.5, label=f'Analytical {integral_label}')
    axes[current_ax_idx].plot(distances, numerical_values,  'r:', lw=1.5, label=f'Numerical (roll) {integral_label}')
    axes[current_ax_idx].set_xlabel('Distance d (Bohr)')
    axes[current_ax_idx].set_ylabel(integral_label)
    axes[current_ax_idx].legend()
    axes[current_ax_idx].grid(True)
    axes[current_ax_idx].set_title('Integral vs. Distance')
    current_ax_idx+=1
    
    # Plot 2: Use the new plot_two_functions
    plot_two_functions(x_grid_plot, func1_on_grid, func2_initial_on_grid,
                       title='Initial Functions on Grid', xlabel='x (Bohr)', ylabel='Amplitude',
                       y1_label=func1_label, y2_label=func2_label, fig=fig, ax=axes[current_ax_idx])
    current_ax_idx+=1

    if bDiff:
        difference = (analytical_values - numerical_values) * scDiff
        axes[current_ax_idx].plot(distances, difference, 'k-', label=f'Difference (Ana - Num) * {scDiff}')
        axes[current_ax_idx].set_xlabel('Distance d (Bohr)')
        axes[current_ax_idx].set_ylabel('Difference')
        axes[current_ax_idx].legend()
        axes[current_ax_idx].grid(True)
        axes[current_ax_idx].set_title('Difference Plot')
        current_ax_idx+=1

    if bRatio:
        # Avoid division by zero or very small numbers for ratio
        ratio = analytical_values / (numerical_values + np.sign(numerical_values)*epsilon_ratio + epsilon_ratio*(numerical_values==0))
        axes[current_ax_idx].plot(distances, ratio, 'purple', ls='-.', label='Ratio (Ana / Num)')
        axes[current_ax_idx].set_xlabel('Distance d (Bohr)')
        axes[current_ax_idx].set_ylabel('Ratio')
        axes[current_ax_idx].axhline(1.0, color='gray', linestyle=':', linewidth=0.8)
        axes[current_ax_idx].legend()
        axes[current_ax_idx].grid(True)
        axes[current_ax_idx].set_title('Ratio Plot')

    plt.tight_layout(rect=[0, 0, 1, 0.96]) # Adjust layout to make space for suptitle
    plt.show()

# --- Main test execution ---
def run_all_tests():
    # --- Grid and General Parameters ---
    num_grid_points = 3001 # Odd number to have a point at zero
    grid_half_width = 20.0 # Make grid large enough: e.g., -20 to 20 Bohr
    x_grid = np.linspace(-grid_half_width, grid_half_width, num_grid_points)
    dx_grid_step = x_grid[1] - x_grid[0]
    
    # Max shift in terms of indices. Corresponds to max_distance = max_shift_idx * dx_grid_step
    max_distance_plot = 8.0 # Bohr
    max_shift_idx = int(max_distance_plot / dx_grid_step)

    print(f"Grid: {num_grid_points} points from {-grid_half_width} to {grid_half_width} Bohr (step: {dx_grid_step:.4f})")
    print(f"Max shift index: {max_shift_idx} (corresponds to distance {max_shift_idx*dx_grid_step:.2f} Bohr)")

    # Basis function parameters (can be varied for different tests)
    alpha_1 = 0.8  # For Gaussian 1
    alpha_2 = 0.5  # For Gaussian 2
    X_center_stationary = 0.0 # Stationary function centered at 0

    # --- 1. Overlap Integral Test: <phi1(x) | phi2(x-d)> ---
    # phi1 is stationary (func1), phi2 is rolled (func2)
    print("\n--- Testing Overlap Integral S(d) ---")
    phi1_s_on_grid = gauss_to_grid(x_grid, X_center_stationary, alpha_1) # Stationary
    phi2_s_on_grid = gauss_to_grid(x_grid, X_center_stationary, alpha_2) # This will be rolled

    distances_s, num_s = compute_numerical_integral_vs_distance(x_grid, phi1_s_on_grid, phi2_s_on_grid, max_shift_idx)
    ana_s = np.array([overlap_element_gauss(alpha_1, X_center_stationary, alpha_2, X_center_stationary + d) for d in distances_s])
    
    plot_integral_vs_distance(distances_s, ana_s, num_s,
                              x_grid, phi1_s_on_grid, phi2_s_on_grid,
                              title=f'Overlap Integral S(d) for $\\alpha_1=${alpha_1}, $\\alpha_2=${alpha_2}',
                              integral_label='S(d)',
                              func1_label=f'Stationary $\\phi_1(x)$ ($\\alpha_1=${alpha_1})',
                              func2_label=f'Initial $\\phi_2(x)$ ($\\alpha_2=${alpha_2})')

    # --- 2. Kinetic Energy Integral Test: <phi1(x-d) | -1/2 d^2/dx^2 | phi2(x)> ---
    # According to user: "-1/2 d^2/dx^2 phi_j" is stationary (func1). "phi_i" is rolled (func2).
    # So, Integral( (-1/2 d^2/dx^2 phi_stationary(x)) * phi_rolled(x-d) )
    print("\n--- Testing Kinetic Energy Integral T(d) ---")
    # func1: stationary, -1/2 d^2/dx^2 phi_stat(x)
    d2phi_stat_on_grid_analytical = d2_gaussian_dx2_on_grid(x_grid, X_center_stationary, alpha_1) # alpha_1 for stationary
    d2phi_stat_on_grid_numerical = numerical_d2_gaussian_dx2_on_grid(x_grid, X_center_stationary, alpha_1)

    # Plot comparison of analytical vs numerical second derivative
    plot_two_functions(x_grid, d2phi_stat_on_grid_analytical, d2phi_stat_on_grid_numerical,
                       title=rf'Comparison of $-1/2 d^2\phi/dx^2$ for $\alpha={alpha_1}$',
                       xlabel='x (Bohr)', ylabel='Amplitude',
                       y1_label=r'Analytical $-1/2 d^2\phi/dx^2$ then grid',
                       y2_label=r'Grid then numerical $-1/2 d^2\phi/dx^2$')

    # func2: rolled, phi_roll(x)
    phi_roll_t_on_grid = gauss_to_grid(x_grid, X_center_stationary, alpha_2) # alpha_2 for rolled

    distances_t, num_t = compute_numerical_integral_vs_distance(x_grid, d2phi_stat_on_grid_analytical, phi_roll_t_on_grid, max_shift_idx)
    

    ana_t = np.zeros_like(distances_t, dtype=float) 
    for i, d_val in enumerate(distances_t):
        # Analytical: <phi_rolled(x-d) | -1/2 d^2/dx^2 | phi_stationary(x)>
        # phi_rolled is at X_center_stationary + d_val
        # phi_stationary is at X_center_stationary
        s_for_t = overlap_element_gauss(alpha_2, X_center_stationary + d_val, alpha_1, X_center_stationary)
        ana_t[i] = kinetic_element_gauss(alpha_2, X_center_stationary + d_val, alpha_1, X_center_stationary, s_for_t)
    plot_integral_vs_distance(distances_t, ana_t, num_t,
                              x_grid, d2phi_stat_on_grid_analytical, phi_roll_t_on_grid,
                              title=rf'Kinetic Integral T(d) for $\alpha_{{stat}}={alpha_1}$, $\alpha_{{roll}}={alpha_2}$',
                              integral_label='T(d)',
                              func1_label=rf'Stationary $-1/2 d^2/dx^2 \phi_{{stat}}(x)$ ($\alpha_1={alpha_1}$)',
                              func2_label=rf'Initial $\phi_{{roll}}(x)$ ($\alpha_2={alpha_2}$)')

    # --- 3. Coulomb Integral Test: <phi_P(x-d) | -Z/|x-X_k| | phi_P(x-d)> ---
    # More precisely: Integral( (-Z/|x-X_k|) * phi_Product(x-d) )
    # User: "1/r potential" is stationary (func1). "gaussian" (product Gaussian) is rolled (func2).
    # So, Integral( (-Z/|x - X_nuc|) * phi_Product(x-d) )
    print("\n--- Testing Coulomb Integral V(d) ---")
    Z_nuc = 1.0
    X_nuc_stationary = 0.0 # Nucleus for potential is at origin
    alpha_nuc_smear  = 4.2 # Exponent for Gaussian smearing of the nucleus

    # func1: stationary, potential from a smeared nucleus
    smeared_potential_on_grid = get_smeared_coulomb_potential_on_grid(x_grid, X_nuc_stationary, Z_nuc, alpha_nuc_smear)
    
    # func2: rolled, product Gaussian phi_P(x) formed from two Gaussians (alpha_p1, alpha_p2) centered at X_center_stationary
    # These alpha_p1, alpha_p2 are components for the product Gaussian, not necessarily alpha_1, alpha_2 from above.
    alpha_p1 = 0.6 
    alpha_p2 = 0.4
    p_el_product = alpha_p1 + alpha_p2 # Exponent of the electron product Gaussian

    # Product Gaussian initially centered at X_center_stationary (0.0)
    phi_P_on_grid_analytical = get_product_gaussian_on_grid(x_grid, alpha_p1, X_center_stationary, alpha_p2, X_center_stationary)
    phi_P_on_grid_numerical = numerical_product_gaussian_on_grid(x_grid, alpha_p1, X_center_stationary, alpha_p2, X_center_stationary)

    # Plot comparison of analytical vs numerical product Gaussian
    plot_two_functions(x_grid, phi_P_on_grid_analytical, phi_P_on_grid_numerical,
                       title=rf'Comparison of Product Gaussian for $\alpha_{{p1}}={alpha_p1}$, $\alpha_{{p2}}={alpha_p2}$',
                       xlabel='x (Bohr)', ylabel='Amplitude',
                       y1_label=r'Analytical Product Gaussian $\phi_P(x)$', # No change needed here
                       y2_label=r'Numerical Product Gaussian (grid1 * grid2)') # No change needed here
    
    # For the integral test, func2 is the analytical product Gaussian
    phi_P_for_integral = phi_P_on_grid_analytical 

    distances_v, num_v = compute_numerical_integral_vs_distance(x_grid, smeared_potential_on_grid, phi_P_for_integral, max_shift_idx)

    # Analytical calculation for interaction of smeared nucleus with electron product Gaussian
    ana_v = np.zeros_like(distances_v, dtype=float)
    # The "charge" of the electron product Gaussian phi_P(x) = K_p1p2 * exp(-p_el * (x-X_P_center)^2)
    # is its integral: K_p1p2 * sqrt(pi/p_el).
    # This is equivalent to the overlap S_p1p2 if p1 and p2 were centered at the same point.
    Q_el_product = overlap_element_gauss(alpha_p1, X_center_stationary, alpha_p2, X_center_stationary)
    # Correct total charge for the smeared nucleus (source of V_smeared_nuc)
    Q_nuc_total_smeared = -Z_nuc * np.sqrt(np.pi / alpha_nuc_smear)

    for i, d_val in enumerate(distances_v):
        X_el_product_center = X_center_stationary + d_val # Center of the rolled electron product
        ana_v[i] = analytical_gaussian_charges_interaction(
            Q1=Q_nuc_total_smeared, alpha1=alpha_nuc_smear, X1=X_nuc_stationary,
            Q2=Q_el_product, alpha2=p_el_product, X2=X_el_product_center
        )
                                         
    plot_integral_vs_distance(distances_v, ana_v, num_v,
                              x_grid, smeared_potential_on_grid, phi_P_for_integral,
                              title=rf'Coulomb Integral V(d) for Smeared Nucleus ($\alpha_{{nuc}}={alpha_nuc_smear}$) \& Product Gaussian ($\alpha_{{p1}}={alpha_p1}$, $\alpha_{{p2}}={alpha_p2}$)',
                              integral_label='V(d)', bDiff=True, bRatio=True, scDiff=1.0, # Enable diff and ratio plots
                              func1_label=rf'Smeared Potential ($\alpha_{{nuc}}={alpha_nuc_smear}$, Z={Z_nuc})',
                              func2_label=r'Initial Product Gaussian $\phi_P(x)$')

if __name__ == '__main__':
    run_all_tests()
