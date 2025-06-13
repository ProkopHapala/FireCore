#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import lstsq, svd

def morse_potential(z, D, a, r0):
    """Defines the Morse potential function."""
    exp_term = np.exp(-a * (z - r0))
    return D * (exp_term**2 - 2 * exp_term)

def generate_polynomial_library_matrix(xs, degree):
    """
    Generates a matrix of polynomial basis functions up to a given degree.
    Phi_ij = xs_i^j

    Parameters:
    -----------
    xs : np.ndarray
        1D array of x-coordinates.
    degree : int
        The maximum degree of the polynomial (e.g., degree 8 means x^0 to x^8).

    Returns:
    --------
    phi_matrix : np.ndarray
        Matrix of shape (len(xs), degree + 1) where each column p is xs^p.
    """
    num_x_points = len(xs)
    num_library_funcs = degree + 1
    phi_matrix = np.zeros((num_x_points, num_library_funcs))
    for p in range(num_library_funcs):
        phi_matrix[:, p] = xs**p
    return phi_matrix

def find_optimal_analytical_basis(
    xs,
    list_of_sample_ys,
    library_degree,
    num_optimal_basis_K
):
    """
    Finds an optimal set of K analytical basis functions (linear combinations of
    a polynomial library) to represent a set of sample functions.

    Parameters:
    -----------
    xs : np.ndarray
        1D array of x-coordinates.
    list_of_sample_ys : list of np.ndarray
        List of sample functions (y-values).
    library_degree : int
        Maximum degree for the polynomial library (P-1).
        The library will have P = library_degree + 1 functions (x^0 to x^degree).
    num_optimal_basis_K : int
        The desired number of new optimal analytical basis functions (K).

    Returns:
    --------
    optimal_analytical_basis_evaluated : np.ndarray
        A (len(xs), K) matrix. Each column is one of the K new optimal
        analytical basis functions, evaluated on the xs grid.
    Uk_coeffs_for_new_basis : np.ndarray
        A (P, K) matrix. Each column k contains the P coefficients that define
        the k-th new basis function as a linear combination of the P library
        polynomials. (P = library_degree + 1).
    singular_values_coeffs : np.ndarray
        Singular values from the SVD of the coefficient matrix S_coeffs.
    """
    num_x_points = len(xs)
    num_samples_M = len(list_of_sample_ys)
    num_library_funcs_P = library_degree + 1

    # --- Step 0: Scale x-coordinates for numerical stability of polynomials ---
    # Scale xs to be roughly in [0, 1] or [-1, 1] if not already
    # For Morse, z_values are typically > 0. Let's scale to [0,1]
    xs_min = np.min(xs)
    xs_max = np.max(xs)
    xs_scaled = (xs - xs_min) / (xs_max - xs_min) if (xs_max - xs_min) > 1e-9 else xs

    # 1. Generate the library matrix Phi (Nx x P) using scaled xs
    # Each column p of phi_matrix is (xs_scaled)^p
    phi_matrix = generate_polynomial_library_matrix(xs_scaled, library_degree)
    print(f"Shape of polynomial library matrix Phi: {phi_matrix.shape}")

    # 2. Represent Sample Functions in the Full Analytical Library Basis
    # Construct sample_matrix_Y (Nx x M)
    sample_matrix_Y = np.array(list_of_sample_ys).T
    print(f"Shape of sample matrix Y: {sample_matrix_Y.shape}")

    # Solve Y approx Phi @ S_coeffs for S_coeffs (P x M)
    # S_coeffs[p, j] is the coefficient of the p-th library function for the j-th sample
    S_coeffs_matrix, residuals, rank, s_lstsq = lstsq(phi_matrix, sample_matrix_Y)
    print(f"Shape of S_coeffs_matrix (coefficients of samples in library basis): {S_coeffs_matrix.shape}")

    # 3. Find Principal Components in the Coefficient Space (SVD on S_coeffs)
    # S_coeffs = U_svd @ Sigma_svd @ Vh_svd
    # Columns of U_svd are principal directions in the P-dim coefficient space.
    U_svd_coeffs, singular_values_coeffs, Vh_svd_coeffs = svd(S_coeffs_matrix, full_matrices=False)
    # U_svd_coeffs has shape (P, min(P, M))

    if num_optimal_basis_K > U_svd_coeffs.shape[1]:
        print(f"Warning: num_optimal_basis_K ({num_optimal_basis_K}) is greater than available "
              f"components ({U_svd_coeffs.shape[1]}). Using {U_svd_coeffs.shape[1]}.")
        num_optimal_basis_K = U_svd_coeffs.shape[1]

    # Uk_coeffs_for_new_basis (P x K) contains coefficients defining new basis
    # Each column k of Uk_coeffs_for_new_basis defines the k-th new analytical basis function
    # b_k(x_scaled) = sum_{p=0}^{P-1} (Uk_coeffs_for_new_basis)_{pk} * (x_scaled)^p
    Uk_coeffs_for_new_basis = U_svd_coeffs[:, :num_optimal_basis_K]
    print(f"Shape of Uk_coeffs_for_new_basis (defines new analytical basis): {Uk_coeffs_for_new_basis.shape}")

    # 4. Construct the Optimal Analytical Basis Functions (evaluated on the grid)
    # B_new_evaluated (Nx x K) = phi_matrix @ Uk_coeffs_for_new_basis
    optimal_analytical_basis_evaluated = phi_matrix @ Uk_coeffs_for_new_basis
    print(f"Shape of optimal_analytical_basis_evaluated: {optimal_analytical_basis_evaluated.shape}")

    return optimal_analytical_basis_evaluated, Uk_coeffs_for_new_basis, singular_values_coeffs, (xs_min, xs_max)


def reconstruct_functions_analytical_basis(original_ys_list, optimal_analytical_basis_evaluated):
    """
    Reconstructs functions using the provided optimal analytical basis.
    This basis is NOT necessarily orthonormal, so we need to solve a small lstsq problem.
    """
    reconstructed_ys_list = []
    coefficients_list = []
    A_basis = optimal_analytical_basis_evaluated # Shape (Nx, K)

    for ys_original in original_ys_list: # ys_original is (Nx,)
        # Solve A_basis @ coeffs = ys_original for coeffs (K,)
        coeffs, _, _, _ = lstsq(A_basis, ys_original)
        reconstructed_ys = A_basis @ coeffs
        
        reconstructed_ys_list.append(reconstructed_ys)
        coefficients_list.append(coeffs)
    return reconstructed_ys_list, coefficients_list


# --- Example Usage ---
if __name__ == "__main__":
    print("--- Finding Optimal ANALYTICAL Basis for Morse Potentials ---")

    # 1. Define common z-coordinates
    z_values = np.linspace(1.0, 10.0, 200)

    # 2. Generate a set of sample Morse potentials
    num_samples = 16
    sample_morse_potentials_ys = []
    np.random.seed(42)
    D_morse = 0.1
    r0_morse = 3.0
    a_min_sample, a_max_sample = 0.8, 2.5
    a_values_sample = np.linspace(a_min_sample, a_max_sample, num_samples)

    print(f"\nGenerating {num_samples} Morse potential samples...")
    for i in range(num_samples):
        a_current = a_values_sample[i]
        ys = morse_potential(z_values, D_morse, a_current, r0_morse)
        sample_morse_potentials_ys.append(ys)

    # Plot some of the sample functions
    plt.figure(figsize=(10, 6))
    for i in range(min(5, num_samples)):
        plt.plot(z_values, sample_morse_potentials_ys[i], label=f'Morse Sample {i+1} (a={a_values_sample[i]:.2f})', alpha=0.7)
    plt.title('A Few Original Morse Potential Samples')
    plt.xlabel('z (Å)')
    plt.ylabel('Potential (eV)')
    plt.legend()
    plt.grid(True)
    plt.savefig("sample_morse_potentials_analytical.png")
    plt.show()

    # 3. Define parameters for the analytical basis
    polynomial_library_degree = 8  # P-1, so P = 9 library functions (z^0 to z^8)
    num_optimal_analytical_K = 4   # Desired number of new basis functions

    print(f"\nAttempting to find {num_optimal_analytical_K} optimal analytical basis functions "
          f"from a polynomial library of degree {polynomial_library_degree}...")

    # 4. Find the optimal analytical basis
    optimal_basis_eval, basis_coeffs_U_k, s_coeffs, (z_min_scale, z_max_scale) = \
        find_optimal_analytical_basis(
            z_values,
            sample_morse_potentials_ys,
            polynomial_library_degree,
            num_optimal_analytical_K
        )

    # Plot the singular values from SVD of S_coeffs
    plt.figure(figsize=(8, 5))
    plt.plot(range(1, len(s_coeffs) + 1), s_coeffs, 'o-')
    plt.title('Singular Values from SVD of Sample Coefficients Matrix (S_coeffs)')
    plt.xlabel('Component Number in Coefficient Space')
    plt.ylabel('Singular Value')
    plt.axvline(num_optimal_analytical_K, color='r', linestyle='--', label=f'Selected K={num_optimal_analytical_K}')
    plt.legend()
    plt.grid(True)
    plt.yscale('log') # Usually a wide range of values
    plt.savefig("singular_values_S_coeffs_analytical.png")
    plt.show()

    # Plot the derived optimal analytical basis functions
    plt.figure(figsize=(12, 7))
    for k in range(optimal_basis_eval.shape[1]):
        plt.plot(z_values, optimal_basis_eval[:, k], label=f'Optimal Analytical Basis {k+1}')
    plt.title(f'{num_optimal_analytical_K} Optimal Analytical Basis Functions (Polynomials in scaled z)')
    plt.xlabel('z (Å)')
    plt.ylabel('Basis Function Value')
    plt.legend()
    plt.grid(True)
    plt.savefig("optimal_analytical_basis_functions.png")
    plt.show()

    # Print the analytical expressions for the new basis functions
    print("\n--- Analytical Expressions for the Optimal Basis Functions ---")
    print(f"Note: z_scaled = (z - {z_min_scale:.3f}) / {(z_max_scale - z_min_scale):.3f}")
    for k in range(basis_coeffs_U_k.shape[1]): # Iterate through K columns
        expr_parts = []
        for p in range(basis_coeffs_U_k.shape[0]): # Iterate through P rows (polynomial degrees)
            coeff_val = basis_coeffs_U_k[p, k]
            if abs(coeff_val) > 1e-6: # Only include significant terms
                if p == 0:
                    expr_parts.append(f"{coeff_val:.3e}")
                elif p == 1:
                    expr_parts.append(f"({coeff_val:.3e} * z_scaled)")
                else:
                    expr_parts.append(f"({coeff_val:.3e} * z_scaled^{p})")
        print(f"b_{k+1}(z_scaled) = {' + '.join(expr_parts)}")


    # 5. Reconstruct the sample functions using the optimal analytical basis
    print(f"\nReconstructing sample Morse potentials using the {num_optimal_analytical_K} optimal analytical basis functions...")
    reconstructed_ys_list, _ = reconstruct_functions_analytical_basis(
        sample_morse_potentials_ys, optimal_basis_eval
    )

    # Plot a comparison for a few sample functions
    num_to_plot_comparison = 3
    plt.figure(figsize=(12, num_to_plot_comparison * 4))
    total_rmse = 0
    for i in range(num_to_plot_comparison):
        original = sample_morse_potentials_ys[i]
        reconstructed = reconstructed_ys_list[i]
        rmse = np.sqrt(np.mean((original - reconstructed)**2))
        total_rmse += rmse
        
        plt.subplot(num_to_plot_comparison, 1, i + 1)
        plt.plot(z_values, original, label=f'Original Morse (a={a_values_sample[i]:.2f})', color='blue', alpha=0.8)
        plt.plot(z_values, reconstructed, label=f'Reconstructed (RMSE: {rmse:.4e})', color='red', linestyle='--')
        plt.title(f'Morse Sample {i+1}: Original vs. Reconstructed (Analytical Basis)')
        plt.xlabel('z (Å)')
        plt.ylabel('Potential (eV)')
        plt.legend()
        plt.grid(True)
    
    avg_rmse = total_rmse / num_to_plot_comparison
    print(f"Average RMSE for the first {num_to_plot_comparison} plotted reconstructions: {avg_rmse:.4e}")

    plt.tight_layout()
    plt.savefig("reconstruction_comparison_analytical.png")
    plt.show()

    print("\n--- Optimal Analytical Basis Demonstration Completed ---")

ls