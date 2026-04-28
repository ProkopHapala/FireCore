# /home/prokophapala/git/FireCore/doc/py/LLCAO1D/test_PolyIntegrals1D_plot.py

import numpy as np
import matplotlib.pyplot as plt

# trun of line-brakes when printing numpy arrays
np.set_printoptions(linewidth=1000)

# ##############################################################################
# ####### Generic Numerical Utilities (Normalization, Operators, Integrals) ######
# ##############################################################################

def numerical_normalize_on_grid(x_grid, phi_on_grid):
    """
    Normalizes a function phi(x) provided on a grid such that Int[phi(x)^2 dx] = 1.
    """
    if np.all(phi_on_grid == 0):
        return phi_on_grid # Avoid division by zero for null function
    norm_sq = np.trapz(phi_on_grid**2, x_grid)
    if norm_sq < 1e-12: # effectively zero
        print("Warning: Norm squared is very small, function might be effectively zero.")
        return phi_on_grid
    return phi_on_grid / np.sqrt(norm_sq)

def apply_kinetic_operator_on_grid(x_vals, phi_on_grid):
    """
    Computes -1/2 * d^2/dx^2 phi(x) for phi(x) given on a grid,
    using numerical differentiation.
    """
    d_phi_dx = np.gradient(phi_on_grid, x_vals, edge_order=2)
    d2_phi_dx2 = np.gradient(d_phi_dx, x_vals, edge_order=2)
    return -0.5 * d2_phi_dx2

def num_int_vs_dist(x_vals, f1_grid, f2_grid, max_sh_idx):
    """
    Computes Integral(func1(x) * func2(x-d))dx for various distances d
    by rolling f2_grid.
    f1_grid is stationary. f2_grid is shifted.
    """
    dx = x_vals[1] - x_vals[0]
    nx = len(x_vals)
    
    if max_sh_idx >= nx // 2:
        print("Warning: max_shift_idx is large relative to grid size. Periodic boundary effects might be significant.")

    dists = np.arange(0, max_sh_idx + 1) * dx
    num_integrals = np.zeros_like(dists, dtype=float)

    for i, sh_idx in enumerate(range(max_sh_idx + 1)):
        # Positive sh_idx rolls elements to the right, effectively f2(x - sh_idx*dx)
        f2_rolled = np.roll(f2_grid, sh_idx)
        integrand = f1_grid * f2_rolled
        num_integrals[i] = np.trapz(integrand, dx=dx)
        
    return dists, num_integrals

# ##############################################################################
# ####### Polynomial Basis Functions with Finite Support                 #######
# ##############################################################################

def phi_slater_poly(x_grid, Xc, w, n, normalize=True):
    """
    Symmetric polynomial basis function: (1 - |x-Xc|/w)^n for |x-Xc| < w.
    Approximates a Slater-like orbital (exponential decay) but with finite support.
    It is computationally cheaper than true Slater orbitals.
    Support is (Xc-w, Xc+w).
    w: half-width of the support.
    n: power, controls the sharpness of decay.
    """
    if w <= 0:
        raise ValueError("Width 'w' must be positive.")
    abs_rel_x = np.abs(x_grid - Xc) / w
    # Ensure base is non-negative before applying power
    core_val = np.maximum(0.0, 1.0 - abs_rel_x)
    phi_vals = np.where(abs_rel_x < 1, core_val**n, 0.0)
    
    if normalize:
        return numerical_normalize_on_grid(x_grid, phi_vals)
    return phi_vals

def phi_gauss_poly(x_grid, Xc, w, n, normalize=True):
    """
    Symmetric polynomial basis function: (1 - ((x-Xc)/w)^2)^n for |x-Xc| < w.
    Approximates a Gaussian-like orbital but with finite support.
    It is computationally cheaper than true Gaussian orbitals for some operations.
    Support is (Xc-w, Xc+w).
    w: half-width of the support.
    n: power, controls the sharpness of decay.
    """
    if w <= 0:
        raise ValueError("Width 'w' must be positive.")
    rel_x_sq = ((x_grid - Xc) / w)**2
    # Ensure base is non-negative before applying power (core_val > 0 implies rel_x_sq < 1)
    core_val = np.maximum(0.0, 1.0 - rel_x_sq)
    phi_vals = np.where(core_val > 0, core_val**n, 0.0) # Condition core_val > 0 is equivalent to rel_x_sq < 1

    if normalize:
        return numerical_normalize_on_grid(x_grid, phi_vals)
    return phi_vals

# ##############################################################################
# ####### Plotting Functions                                             #######
# ##############################################################################

# ##############################################################################
# ####### Plotting Functions                                             #######
# ##############################################################################

def plot1d(x_data, plots_data, title="", xlabel="x", ylabel="y", fig=None, ax=None, show_plot_after=False):
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 4))

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

def plot_poly_integrals_summary(x_grid, phi1_zero_grid, phi2_zero_grid,
                                d_vals, num_overlap_vals, num_kinetic_vals, title,
                                ana_overlap_vals=None,
                                ana_kinetic_vals=None
                                ):
    """
    Plots basis functions centered at zero, and their numerical overlap and kinetic integrals vs distance.
    Optionally plots analytical integrals if provided.
    """
    fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=False)
    fig.suptitle(title, fontsize=16)

    # Panel 1: Basis Functions centered at zero
    plot1d(x_grid, [(phi1_zero_grid, 'Basis 1 (centered)', 'b-', 1.5), (phi2_zero_grid, 'Basis 2 (centered)', 'r-', 1.5)],
           title='Basis Functions Centered at Zero', xlabel='x (Bohr)', ylabel='Amplitude', ax=axes[0], show_plot_after=False)

    # Panel 2: Overlap Integral vs. Distance
    overlap_plots_data = []
    if ana_overlap_vals is not None:
        overlap_plots_data.append((ana_overlap_vals, 'Analytical Overlap S(d)', 'b-', 0.7))
    num_overlap_label = 'Numerical Overlap S(d)' if ana_overlap_vals is None else 'Numerical S(d)'
    overlap_plots_data.append((num_overlap_vals, num_overlap_label, 'r:', 1.5))
    plot1d(d_vals, overlap_plots_data,
           title='Overlap Integral vs. Distance', xlabel='Distance d (Bohr)', ylabel='S(d)', ax=axes[1], show_plot_after=False)

    # Panel 3: Kinetic Energy Integral vs. Distance
    kinetic_plots_data = []
    if ana_kinetic_vals is not None:
        kinetic_plots_data.append((ana_kinetic_vals, 'Analytical Kinetic T(d)', 'b-', 0.7))
    num_kinetic_label = 'Numerical Kinetic T(d)' if ana_kinetic_vals is None else 'Numerical T(d)'
    kinetic_plots_data.append((num_kinetic_vals, num_kinetic_label, 'r:', 1.5))
    plot1d(d_vals, kinetic_plots_data,
           title='Kinetic Integral vs. Distance', xlabel='Distance d (Bohr)', ylabel='T(d)', ax=axes[2], show_plot_after=False)


    plt.tight_layout(rect=[0, 0, 1, 0.96]) # Adjust layout to make space for suptitle
    plt.show()

# ##############################################################################
# ####### Test Functions for Polynomial Basis                            #######
# ##############################################################################

def test_poly_overlap(x_grid_test, basis_func_type_1, params1, basis_func_type_2, params2, max_shift_idx_test):
    """
    Tests overlap integral <phi1(x) | phi2(x-d)> for polynomial basis functions.
    phi1 is stationary, phi2 is shifted.
    params = (Xc, w, n)
    """
    # This function is now integrated into run_poly_integral_test
    # Keeping it here for context if needed, but it won't be called directly
    pass

def test_poly_kinetic(x_grid_test, basis_func_type_1, params1, basis_func_type_2, params2, max_shift_idx_test):
    """
    Tests kinetic energy integral <phi1(x) | -1/2 d^2/dx^2 | phi2(x-d)> for polynomial basis functions.
    phi1 is stationary. Operator acts on phi2, which is then shifted.
    params = (Xc, w, n)
    """
    # This function is now integrated into run_poly_integral_test
    # Keeping it here for context if needed, but it won't be called directly
    pass

def run_poly_integral_test(x_grid_test, basis_func_type_1, params1, basis_func_type_2, params2, max_shift_idx_test):
    """
    Runs overlap and kinetic integral tests for a pair of polynomial basis functions
    and plots the results in a 3-panel summary figure.
    phi1 is stationary, phi2 is shifted.
    params = (Xc, w, n)
    """
    Xc1, w1, n1 = params1
    Xc2_initial, w2, n2 = params2 # Xc2_initial is the center before shifting

    title_str = f"{basis_func_type_1.__name__}(w={w1},n={n1}) vs {basis_func_type_2.__name__}(w={w2},n={n2})"
    print(f"\n--- Running Integral Tests for {title_str} ---")

    # --- Prepare functions for integral calculation ---
    # phi1 is stationary at Xc1
    phi1_grid = basis_func_type_1(x_grid_test, Xc1, w1, n1, normalize=True)
    # phi2 is initially at Xc2_initial, will be rolled
    phi2_grid_initial = basis_func_type_2(x_grid_test, Xc2_initial, w2, n2, normalize=True)

    # --- Calculate Overlap Integral ---
    # <phi1(x) | phi2(x-d)>
    # phi2 is initially at Xc2_initial. We roll it, so effectively its center becomes Xc2_initial + d
    d_vals, num_overlap_vals = num_int_vs_dist(x_grid_test, phi1_grid, phi2_grid_initial, max_shift_idx_test)

    # --- Calculate Kinetic Energy Integral ---
    # <phi1(x) | -1/2 d^2/dx^2 | phi2(x-d)>
    # Apply kinetic operator to the initial (unshifted) phi2
    kin_op_phi2_grid_initial = apply_kinetic_operator_on_grid(x_grid_test, phi2_grid_initial)
    # Op*phi2 is initially centered at Xc2_initial. We roll it.
    d_vals_kin, num_kinetic_vals = num_int_vs_dist(x_grid_test, phi1_grid, kin_op_phi2_grid_initial, max_shift_idx_test)
    # Note: d_vals_kin should be the same as d_vals, but we keep it separate for clarity

    # --- Prepare functions for visualization (centered at zero) ---
    phi1_zero_grid = basis_func_type_1(x_grid_test, 0.0, w1, n1, normalize=True)
    phi2_zero_grid = basis_func_type_2(x_grid_test, 0.0, w2, n2, normalize=True)

    # --- Plot Summary ---
    plot_poly_integrals_summary(
        x_grid_test,
        phi1_zero_grid,
        phi2_zero_grid,
        d_vals,
        num_overlap_vals,
        num_kinetic_vals,
        title=f'Integral Tests for {title_str}'
    )

# ##############################################################################
# ####### Main Execution Block                                           #######
# ##############################################################################

if __name__ == '__main__':
    # --- Grid and General Parameters ---
    N_grid = 2001  # Number of grid points
    hw_grid = 15.0 # Grid half-width (from -hw_grid to hw_grid)
    x_grid = np.linspace(-hw_grid, hw_grid, N_grid)
    dx_grid = x_grid[1] - x_grid[0]
    
    max_d_plot = 8.0 # Max distance for integral plots
    max_shift_idx = int(max_d_plot / dx_grid)

    print(f"Grid: {N_grid} points from {-hw_grid} to {hw_grid} Bohr (step: {dx_grid:.4f})")
    print(f"Max shift index: {max_shift_idx} (corresponds to distance {max_shift_idx*dx_grid:.2f} Bohr)")

    # --- Test Parameters for Polynomial Basis Functions ---
    # For integral tests, phi1 is stationary, phi2 is shifted.
    # Xc1_stat and Xc2_init define their initial positions before rolling phi2.
    # For symmetric functions, Xc is the center.
    Xc1_stat = 0.0 
    Xc2_init = 0.0 # Initial center before shifting for integral calculation

    # Ensure n >= 2 for smoother derivatives at the cutoff
    n_val_1 = 2.0
    n_val_2 = 2.0

    # --- Test phi_gauss_poly vs phi_gauss_poly ---
    print("\n" + "="*30 + " Testing phi_gauss_poly vs phi_gauss_poly " + "="*30)
    w_gp1 = 3.0
    w_gp2 = 2.5
    params_gp1 = (Xc1_stat, w_gp1, n_val_1)
    params_gp2 = (Xc2_init, w_gp2, n_val_2)

    run_poly_integral_test(x_grid, phi_gauss_poly, params_gp1, phi_gauss_poly, params_gp2, max_shift_idx)

    # --- Test phi_slater_poly vs phi_slater_poly ---
    print("\n" + "="*30 + " Testing phi_slater_poly vs phi_slater_poly " + "="*30)
    w_sp1 = 4.0
    w_sp2 = 3.0
    params_sp1 = (Xc1_stat, w_sp1, n_val_1)
    params_sp2 = (Xc2_init, w_sp2, n_val_2)

    run_poly_integral_test(x_grid, phi_slater_poly, params_sp1, phi_slater_poly, params_sp2, max_shift_idx)

    # --- Mixed Basis Functions Test (Optional, commented out as per request) ---
    # print("\n" + "="*30 + " Testing Mixed Basis Functions (GaussPoly vs SlaterPoly) " + "="*30)
    # # phi1: GaussPoly, stationary
    # # phi2: SlaterPoly, initial Xc=0.0, to be rolled
    # w_mixed_gp = 3.0
    # w_mixed_sp = 2.5
    # params_mixed_gp = (Xc1_stat, w_mixed_gp, n_val_1)
    # params_mixed_sp = (Xc2_init, w_mixed_sp, n_val_2)
    # run_poly_integral_test(x_grid, phi_gauss_poly, params_mixed_gp, phi_slater_poly, params_mixed_sp, max_shift_idx)

    print("\nAll tests completed.")