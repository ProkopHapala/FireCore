#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
# Import plotting utilities
from utils import plot_singular_values, plot_1d_profiles, plotFunctionApprox, morse_potential, gen_morse_prms, gen_morse_curves

def find_optimal_basis_svd(x, ys_list, n_basis):
    """Find optimal 1D basis functions using SVD.
    
    Args:
        x: Common x-coordinates (1D array)
        ys_list: List of sample functions (1D arrays)
        n_basis: Number of basis functions to extract
    
    Returns:
        basis: Optimal basis functions (columns)
        s: Singular values
        U: Full U matrix from SVD
    """
    if not ys_list:
        raise ValueError("Empty ys_list")
    nx = len(x)
    if any(len(y) != nx for y in ys_list):
        raise ValueError("All y must match x length")
    D = np.array(ys_list).T
    U, s, Vh = np.linalg.svd(D, full_matrices=False)
    if n_basis > U.shape[1]:
        print(f"Warning: Reducing n_basis from {n_basis} to {U.shape[1]}")
        n_basis = U.shape[1]
    
    return U[:, :n_basis], s, U

def reconstruct_functions(ys_list, basis):
    """Reconstruct functions using basis.
    
    Args:
        ys_list: Original functions
        basis: Basis functions (columns)
    
    Returns:
        rec_ys: Reconstructed functions
        coeffs: Coefficients for each function
    """
    rec_ys = []
    coeffs = []
    for y in ys_list:
        c = basis.T @ y
        rec_ys.append(basis @ c)
        coeffs.append(c)
    return rec_ys, coeffs

# --- Example Usage ---
if __name__ == "__main__":
    print("--- Finding Optimal Basis Functions for Morse Potentials using SVD ---")

    # 1. Define common x-coordinates
    z_values = np.linspace(1.0, 10.0, 200) 

    # 2. Generate a set of sample functions
    n_s = 16 
    np.random.seed(42) 

    D_val = 0.1  
    a_rng = (0.8, 2.5)
    r0_rng = (2.5, 4.0) # Vary r0 as well

    print(f"\nGenerating {n_s} Morse potential samples...")
    # Generate parameters first, then curves
    prms_list = gen_morse_prms(n_s, a_rng, r0_rng, D_val)
    sample_functions_ys, _ = gen_morse_curves(z_values, prms=prms_list)
    for i, p in enumerate(prms_list):
        print(f"  Generated Morse sample {i+1}: D={p['D']:.2f}, a={p['a']:.2f}, r0={p['r0']:.2f}")

    # Plot some of the sample functions
    plot_1d_profiles(z_values, np.array(sample_functions_ys), title='A Few Original Morse Potential Samples', max_plot=5, filename="sample_functions.png")

    # 3. Choose the number of basis functions to find
    num_optimal_basis = 4 
    print(f"\nAttempting to find {num_optimal_basis} optimal basis functions...")

    # 4. Find the optimal basis
    optimal_basis_cols, singular_values, U_full = find_optimal_basis_svd(  z_values, sample_functions_ys, num_optimal_basis )
    print(f"Shape of optimal_basis matrix (functions as columns): {optimal_basis_cols.shape}")
    print(f"Singular values (indicate importance): {singular_values[:num_optimal_basis*2]}") 
    plot_singular_values(singular_values, K_opt=num_optimal_basis, filename="singular_values.png")
    plot_1d_profiles(z_values, optimal_basis_cols.T, title=f'{num_optimal_basis} Optimal Basis Functions Derived by SVD', filename="optimal_basis_functions.png" )

    # 5. Reconstruct functions
    print(f"\nReconstructing sample functions using {num_optimal_basis} basis functions...")
    rec_ys, coeffs = reconstruct_functions(sample_functions_ys, optimal_basis_cols)

    # Prepare data for plotFunctionApprox
    nplt = min(5, len(sample_functions_ys))
    ys_approx = []
    total_rmse = 0
    for i in range(nplt):
        rmse = np.sqrt(np.mean((sample_functions_ys[i] - rec_ys[i])**2))
        total_rmse += rmse
        ys_approx.append((rec_ys[i], 0, f'Recon {i+1} (RMSE: {rmse:.1e})'))
    # Plot all comparisons in one figure
    ax1_recon, ax2_recon = plotFunctionApprox( z_values, sample_functions_ys[0], ys_approx, bError=True, errMax=0.1 )
    ax1_recon.set_title(f'Morse Potential Reconstruction (n_basis={num_optimal_basis})')
    plt.savefig("reconstruction_comparison.png")
    plt.show()

    print("\n--- Morse Potential Optimal Basis Demonstration Completed ---")
