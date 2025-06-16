#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

def find_optimal_basis_svd(xs, list_of_ys, num_basis_functions):
    """
    Finds an optimal set of 1D basis functions for a given set of sample functions
    using Singular Value Decomposition (SVD).

    The basis functions are chosen to capture the maximum variance in the sample
    functions, thus minimizing the reconstruction error in a least-squares sense
    when approximating the original functions with this basis. The resulting
    basis functions are orthonormal.

    Parameters:
    -----------
    xs : np.ndarray
        1D array of x-coordinates, common to all sample functions.
    list_of_ys : list of np.ndarray
        A list where each element is a 1D numpy array representing a sample
        function (y-values). All y-arrays must have the same length as xs.
    num_basis_functions : int
        The desired number of optimal basis functions to extract.

    Returns:
    --------
    optimal_basis : np.ndarray
        A 2D numpy array of shape (len(xs), num_basis_functions) where each
        column is an optimal basis function.
    singular_values : np.ndarray
        The singular values corresponding to all components, ordered by importance.
        These indicate how much variance each component captures.
    U_full : np.ndarray
        The full U matrix from SVD (data_matrix = U_full @ np.diag(s) @ Vh_full).
        Its columns are all possible orthonormal basis functions for the provided samples.
    """
    if not list_of_ys:
        raise ValueError("list_of_ys cannot be empty.")
    
    num_x_points = len(xs)
    num_samples = len(list_of_ys)

    if any(len(ys) != num_x_points for ys in list_of_ys):
        raise ValueError("All sample functions in list_of_ys must have the same length as xs.")

    # 1. Construct the data matrix D
    # Each column in D is a sample function ys. Shape: (num_x_points, num_samples)
    data_matrix = np.array(list_of_ys).T  # Transpose to make functions columns

    # 2. Perform Singular Value Decomposition
    # data_matrix = U @ np.diag(s) @ Vh
    # Columns of U are the principal components (optimal basis functions for the y-space)
    # s contains singular values in descending order
    U_full, singular_values, Vh_full = np.linalg.svd(data_matrix, full_matrices=False)
    # full_matrices=False is more efficient and U_full will have shape (num_x_points, min(num_x_points, num_samples))

    if num_basis_functions > U_full.shape[1]:
        print(f"Warning: num_basis_functions ({num_basis_functions}) is greater than the number of "
              f"available non-zero singular value components ({U_full.shape[1]}). "
              f"Returning all {U_full.shape[1]} available components.")
        num_basis_functions = U_full.shape[1]

    # 3. Select the top num_basis_functions from U_full
    # These are the directions of greatest variance in the data.
    optimal_basis = U_full[:, :num_basis_functions]

    return optimal_basis, singular_values, U_full

def reconstruct_functions(original_ys_list, optimal_basis):
    """
    Reconstructs functions using the provided optimal basis.

    Parameters:
    -----------
    original_ys_list : list of np.ndarray
        The list of original y-values for each function.
    optimal_basis : np.ndarray
        The basis functions (columns are orthonormal basis vectors).

    Returns:
    --------
    reconstructed_ys_list : list of np.ndarray
        List of reconstructed y-values.
    coefficients_list : list of np.ndarray
        List of coefficients used for reconstruction for each function.
    """
    reconstructed_ys_list = []
    coefficients_list = []
    for ys_original in original_ys_list:
        # Project ys_original onto the basis: coeffs = optimal_basis.T @ ys_original
        # This works because columns of optimal_basis (from U of SVD) are orthonormal.
        coeffs = optimal_basis.T @ ys_original
        
        # Reconstruct: reconstructed_ys = optimal_basis @ coeffs
        reconstructed_ys = optimal_basis @ coeffs
        
        reconstructed_ys_list.append(reconstructed_ys)
        coefficients_list.append(coeffs)
    return reconstructed_ys_list, coefficients_list

# --- Morse Potential Function (from Fit1D.py) ---
def morse_potential(z, D, a, r0):
    """Defines the Morse potential function."""
    exp_term = np.exp(-a * (z - r0))
    return D * (exp_term**2 - 2 * exp_term)

# --- Example Usage ---
if __name__ == "__main__":
    print("--- Finding Optimal Basis Functions for Morse Potentials using SVD ---")

    # 1. Define common x-coordinates
    # Using z_values similar to Fit1D.py for consistency
    z_values = np.linspace(1.0, 10.0, 200) # Z range for functions

    # 2. Generate a set of sample functions
    num_samples = 16 # Number of Morse potentials to generate
    sample_functions_ys = []
    np.random.seed(42) # for reproducibility

    # Morse parameters (can be adjusted)
    D_morse = 0.1  # Morse depth (eV), kept constant for this example
    r0_morse = 3.0 # Equilibrium distance (Å), kept constant

    # Vary the 'a' parameter for different Morse potentials
    # Let's sample 'a' from a range, e.g., similar to Fit1D.py or a bit wider
    a_min_sample = 0.8
    a_max_sample = 2.5
    a_values_sample = np.linspace(a_min_sample, a_max_sample, num_samples)

    print(f"\nGenerating {num_samples} Morse potential samples...")
    for i in range(num_samples):
        a_current = a_values_sample[i]
        ys = morse_potential(z_values, D_morse, a_current, r0_morse)
        # No noise is added, functions will be smooth
        sample_functions_ys.append(ys)
        print(f"  Generated Morse sample {i+1}: D={D_morse:.2f}, a={a_current:.2f}, r0={r0_morse:.2f}")

    # Plot some of the sample functions
    plt.figure(figsize=(10, 6))
    for i in range(min(5, num_samples)): # Plot first 5 samples
        plt.plot(z_values, sample_functions_ys[i], label=f'Morse Sample {i+1} (a={a_values_sample[i]:.2f})', alpha=0.7)
    plt.title('A Few Original Morse Potential Samples')
    plt.xlabel('z (Å)')
    plt.ylabel('Potential (eV)')
    plt.legend()
    plt.grid(True)
    plt.savefig("sample_functions.png")
    plt.show()

    # 3. Choose the number of basis functions to find
    num_optimal_basis = 4 # As requested
    print(f"\nAttempting to find {num_optimal_basis} optimal basis functions...")

    # 4. Find the optimal basis
    optimal_basis, singular_values, U_full = find_optimal_basis_svd(
        z_values, sample_functions_ys, num_optimal_basis
    )
    print(f"Shape of optimal_basis matrix: {optimal_basis.shape}")
    print(f"Singular values (indicate importance): {singular_values[:num_optimal_basis*2]}") # Show a few more

    # Plot the singular values (scree plot)
    plt.figure(figsize=(8, 5))
    plt.plot(range(1, len(singular_values) + 1), singular_values, 'o-')
    plt.title('Singular Values (Scree Plot)')
    plt.xlabel('Component Number')
    plt.ylabel('Singular Value')
    plt.axvline(num_optimal_basis, color='r', linestyle='--', label=f'Selected: {num_optimal_basis}')
    plt.legend()
    plt.grid(True)
    plt.yscale('log')
    plt.savefig("singular_values.png")
    plt.show()
    
    # Plot the derived optimal basis functions
    plt.figure(figsize=(10, 6))
    for i in range(optimal_basis.shape[1]):
        plt.plot(z_values, optimal_basis[:, i], label=f'Optimal Basis {i+1}')
    plt.title(f'{num_optimal_basis} Optimal Basis Functions Derived by SVD')
    plt.xlabel('z (Å)')
    plt.ylabel('Basis Function Value (Arbitrary Units)')
    plt.legend()
    plt.grid(True)
    plt.savefig("optimal_basis_functions.png")
    plt.show()

    # 5. Reconstruct the sample functions using the optimal basis
    print(f"\nReconstructing sample functions using the {num_optimal_basis} optimal basis functions...")
    reconstructed_ys_list, coefficients_list = reconstruct_functions(sample_functions_ys, optimal_basis)

    # Plot a comparison for a few sample functions
    num_to_plot_comparison = 3
    plt.figure(figsize=(12, num_to_plot_comparison * 4))
    total_rmse = 0
    for i in range(num_to_plot_comparison):
        original = sample_functions_ys[i]
        reconstructed = reconstructed_ys_list[i]
        rmse = np.sqrt(np.mean((original - reconstructed)**2))
        total_rmse += rmse
        
        plt.subplot(num_to_plot_comparison, 1, i + 1)
        plt.plot(z_values, original, label=f'Original Morse (a={a_values_sample[i]:.2f})', color='blue', alpha=0.8)
        plt.plot(z_values, reconstructed, label=f'Reconstructed (RMSE: {rmse:.4e})', color='red', linestyle='--')
        plt.title(f'Morse Sample {i+1}: Original vs. Reconstructed with {num_optimal_basis} Basis Functions')
        plt.xlabel('z (Å)')
        plt.ylabel('Potential (eV)')
        plt.legend()
        plt.grid(True)
    
    avg_rmse = total_rmse / num_to_plot_comparison
    print(f"Average RMSE for the first {num_to_plot_comparison} plotted reconstructions: {avg_rmse:.4f}")

    plt.tight_layout()
    plt.savefig("reconstruction_comparison.png")
    plt.show()

    print("\n--- Morse Potential Optimal Basis Demonstration Completed ---")
