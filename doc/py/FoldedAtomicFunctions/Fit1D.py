import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import lstsq
from scipy.optimize import curve_fit

def gram_schmidt(vectors):
    """Orthogonalizes a set of vectors using the Gram-Schmidt process."""
    basis = []
    for v in vectors:
        w = v - np.sum(np.dot(v, b) * b for b in basis)
        if (w > 1.0e-10).any():  # Checks if the vector is not too small
            basis.append(w / np.linalg.norm(w))
    return np.array(basis)

def exponential_decay(z, alpha):
    """Defines an exponential decay function."""
    return np.exp(-alpha * z)

def morse_potential(z, D, a, r0):
    """Defines the Morse potential function."""
    exp_term = np.exp(-a * (z - r0))
    return D * (exp_term**2 - 2 * exp_term)

# --- Parameters ---
z_values = np.linspace(1.0, 10.0, 200) # Z range for functions
alpha_min = 1.5
alpha_max = 1.8

# --- 1) Orthogonalized Exponential Functions ---
print("\n--- 1) Orthogonalized Exponential Functions ---")
# Create two exponential functions with different alphas
f1 = exponential_decay(z_values, alpha_min)
f2 = exponential_decay(z_values, alpha_max)

# Orthogonalize using Gram-Schmidt
orthogonal_basis = gram_schmidt([f1, f2])

# Plot orthogonal basis functions
plt.figure(figsize=(8, 6))
plt.plot(z_values, orthogonal_basis[0], label=f"Orthogonal Function 1 (α={alpha_min:.1f})")
plt.plot(z_values, orthogonal_basis[1], label=f"Orthogonal Function 2 (α={alpha_max:.1f})")
plt.xlabel("z (Å)")
plt.ylabel("Function Value")
plt.title("Orthogonalized Exponential Functions")
plt.legend()
plt.grid(True)
#plt.savefig("orthogonal_exponentials.png")
#plt.show()

# --- 2) Test Functions ---
print("\n--- 2) Test Functions ---")
num_test_functions = 5
alphas_test = np.linspace(alpha_min, alpha_max, num_test_functions)
test_functions = [exponential_decay(z_values, alpha) for alpha in alphas_test]

# Plot test functions
plt.figure(figsize=(10, 6))
for i, func in enumerate(test_functions):
    plt.plot(z_values, func, label=f"Test Function {i+1} (α={alphas_test[i]:.2f})")
plt.xlabel("z (Å)")
plt.ylabel("Function Value")
plt.title("Test Exponential Functions")
plt.legend()
plt.grid(True)
#plt.savefig("test_exponentials.png")
#plt.show()

# --- 3) Fitting Morse Potentials ---
print("\n--- 3) Fitting Morse Potentials ---")

# Na and Cl parameters (example - replace with actual values)
Na_radius = 1.6
Cl_radius = 2.3
Na_energy = 0.01
Cl_energy = 0.01

# Combined parameters (example)
Rij = Na_radius + Cl_radius
Eij = np.sqrt(Na_energy * Cl_energy) # Geometric mean for energy

# Sample Morse parameters
num_morse_potentials = 3
morse_params_list = []
for i in range(num_morse_potentials):
    # Sample Morse parameters (you can modify the sampling logic)
    D = Eij  # Morse depth (eV)
    a = np.random.uniform(1.0, 2.0)  # Morse decay (1/Å)
    r0 = Rij  # Equilibrium distance (Å)
    morse_params_list.append((D, a, r0))
    print(f"  Morse potential {i+1}: D={D:.3f} eV, a={a:.2f} Å^-1, r0={r0:.2f} Å")

# Generate Morse potentials
morse_potentials = [morse_potential(z_values, D, a, r0) for D, a, r0 in morse_params_list]

# Basis functions for fitting
alphas_basis = [alpha_min, alpha_max, 2 * alpha_min, 2 * alpha_max]
basis_functions = [exponential_decay(z_values, alpha) for alpha in alphas_basis]
basis_matrix = np.array(basis_functions).T  # Shape: (n_z, n_basis)

# Fit and plot
plt.figure(figsize=(12, 8))
for i, V_morse in enumerate(morse_potentials):
    # Linear least squares fitting
    coefficients, residuals, rank, singular_values = lstsq(basis_matrix, V_morse)
    fitted_potential = basis_matrix @ coefficients
    
    # Print fitting results
    print(f"\n  Fitting Morse Potential {i+1}:")
    print(f"    Coefficients: {coefficients}")    

    # --- Weighted Fitting ---
    # Find the zero-crossing point of the Morse potential
    zero_crossing_idx = np.where(np.diff(np.sign(V_morse)) < 0)[0]
    if len(zero_crossing_idx) > 0:
        zero_crossing_z = z_values[zero_crossing_idx[0]]
    else:
        zero_crossing_z = None

    # Define weights: mask repulsive region
    weights = np.ones_like(V_morse)
    if zero_crossing_z is not None:
        repulsive_mask = (z_values < zero_crossing_z - 0.5)
        weights[repulsive_mask] = 0.1  # Reduce weight of repulsive region
        print(f"    Zero crossing at z ≈ {zero_crossing_z:.2f} Å, applying weights.")
    else:
        print("    No zero crossing found, using equal weights.")

    # Apply weights to the least squares problem
    weighted_V_morse = V_morse * weights
    weighted_basis_matrix = basis_matrix * weights[:, np.newaxis] # Apply weights to each basis function

    coefficients_weighted, residuals_weighted, rank_weighted, singular_values_weighted = lstsq(weighted_basis_matrix, weighted_V_morse)

    # Reconstruct fitted potential (using original basis matrix)
    fitted_potential_weighted = basis_matrix @ coefficients_weighted

    # Calculate RMSE with and without weights
    rmse_weighted = np.sqrt(np.mean(((fitted_potential_weighted - V_morse)*weights)**2))

    print(f"    Weighted Coefficients: {coefficients_weighted}")
    print(f"    Weighted RMSE: {rmse_weighted:.6f} eV")

    # --- Plotting ---
    plt.subplot(num_morse_potentials, 1, i + 1)  # (rows, cols, panel_number)
    min_V = np.min(V_morse)
    vmax = -1.5 * min_V
    vmin = 1.5 * min_V
    plt.ylim(vmin, vmax)

    rmse = np.sqrt(np.mean((fitted_potential - V_morse)**2))
    print(f"    RMSE: {rmse:.6f} eV")


    # Plotting
    plt.subplot(num_morse_potentials, 1, i + 1) # (rows, cols, panel_number)
    plt.plot(z_values, V_morse, label=f"Original Morse (D={morse_params_list[i][0]:.3f}, a={morse_params_list[i][1]:.2f}, r0={morse_params_list[i][2]:.2f})")
    plt.plot(z_values, fitted_potential_weighted, '-.', label=f"Fitted Potential (Weighted)")
    plt.plot(z_values, fitted_potential, '--', label="Fitted Potential")
    plt.xlabel("z (Å)")
    plt.ylabel("Potential (eV)")
    plt.title(f"Morse Potential Fitting {i+1}")
    plt.legend()
    plt.grid(True)

plt.tight_layout()  # Adjusts subplot params so that the subplot(s) fits in to the figure area.
#plt.savefig("morse_potential_fitting_weighted.png")
plt.show()

print("\n--- Script Completed ---")
