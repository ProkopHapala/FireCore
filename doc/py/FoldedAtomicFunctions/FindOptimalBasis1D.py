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

# --- Example Usage ---
if __name__ == "__main__":
    print("--- Finding Optimal Basis Functions using SVD ---")

    # 1. Define common x-coordinates
    xs = np.linspace(-5, 5, 200)

    # 2. Generate a set of sample functions
    # Let's create functions that are combinations of a few underlying shapes + noise
    def true_basis_1(x): return np.exp(-x**2 / 2)  # Gaussian
    def true_basis_2(x): return np.sin(x * 1.5)    # Sine wave
    def true_basis_3(x): return 0.2 * x            # Linear ramp

    num_samples = 20
    sample_functions_ys = []
    np.random.seed(42) # for reproducibility
    print(f"\nGenerating {num_samples} sample functions...")
    for i in range(num_samples):
        c1 = np.random.uniform(-2, 2)
        c2 = np.random.uniform(-1.5, 1.5)
        c3 = np.random.uniform(-1, 1)
        noise = np.random.normal(0, 0.1, len(xs)) # Add some noise
        ys = c1 * true_basis_1(xs) + \
             c2 * true_basis_2(xs) + \
             c3 * true_basis_3(xs) + \
             noise
        sample_functions_ys.append(ys)

    # Plot some of the sample functions
    plt.figure(figsize=(10, 6))
    for i in range(min(5, num_samples)): # Plot first 5 samples
        plt.plot(xs, sample_functions_ys[i], label=f'Sample {i+1}', alpha=0.7)
    plt.title('A Few Original Sample Functions')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.grid(True)
    #plt.savefig("sample_functions.png")
    #plt.show()

    # 3. Choose the number of basis functions to find
    num_optimal_basis = 3
    print(f"\nAttempting to find {num_optimal_basis} optimal basis functions...")

    # 4. Find the optimal basis
    optimal_basis, singular_values, U_full = find_optimal_basis_svd(
        xs, sample_functions_ys, num_optimal_basis
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
    #plt.savefig("singular_values.png")
    #plt.show()
    
    # Plot the derived optimal basis functions
    plt.figure(figsize=(10, 6))
    for i in range(optimal_basis.shape[1]):
        plt.plot(xs, optimal_basis[:, i], label=f'Optimal Basis {i+1}')
    plt.title(f'{num_optimal_basis} Optimal Basis Functions Derived by SVD')
    plt.xlabel('x')
    plt.ylabel('Basis Function Value')
    plt.legend()
    plt.grid(True)
    #plt.savefig("optimal_basis_functions.png")
    #plt.show()

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
        plt.plot(xs, original, label='Original Sample', color='blue', alpha=0.8)
        plt.plot(xs, reconstructed, label=f'Reconstructed (RMSE: {rmse:.4f})', color='red', linestyle='--')
        plt.title(f'Sample {i+1}: Original vs. Reconstructed')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.legend()
        plt.grid(True)
    
    avg_rmse = total_rmse / num_to_plot_comparison
    print(f"Average RMSE for the first {num_to_plot_comparison} plotted reconstructions: {avg_rmse:.4f}")

    plt.tight_layout()
    #plt.savefig("reconstruction_comparison.png")
    plt.show()

    print("\n--- Method Demonstration Completed ---")
