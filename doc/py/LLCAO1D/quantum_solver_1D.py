# quantum_solver.py

import numpy as np
from scipy.linalg import eigh
from scipy.special import erf

def boys_function_F0(t):
    """
    Computes the Boys function F0(t) for the nuclear attraction integral.
    F0(t) = 0.5 * sqrt(pi/t) * erf(sqrt(t)).
    A small epsilon is used to avoid division by zero for t=0.
    For t < 0, this formula is not directly applicable (argument of sqrt(t) would be complex).
    However, t = p*(X_p-X_k)^2 should always be >= 0 for positive p.
    """
    if t < 1e-10: # Threshold for t being effectively zero
        return 1.0
    else:
        # Ensure t is not negative to avoid issues with np.sqrt if t is extremely small negative due to precision
        return 0.5 * np.sqrt(np.pi / max(t, 0)) * erf(np.sqrt(max(t,0)))

def overlap_element_gauss(alpha_i, X_i, alpha_j, X_j):
    """ S_ij = <phi_i | phi_j> for normalized Gaussians. """
    p = alpha_i + alpha_j
    norm_i = (2 * alpha_i / np.pi)**0.25
    norm_j = (2 * alpha_j / np.pi)**0.25
    exp_term = np.exp(-(alpha_i * alpha_j / p) * (X_i - X_j)**2)
    integral_term = (np.pi / p)**0.5
    return norm_i * norm_j * exp_term * integral_term

def kinetic_element_gauss(alpha_i, X_i, alpha_j, X_j, s_ij_normalized):
    """ T_ij = <phi_i | -1/2 d^2/dx^2 | phi_j> for normalized Gaussians. """
    p = alpha_i + alpha_j
    # Using the symmetric form: T_ij = S_ij * (alpha_i*alpha_j/p) * (1 - 2*(alpha_i*alpha_j/p)*(X_i-X_j)^2)
    # This is equivalent to S_ij * [alpha_k - 2*alpha_k^2 * (1/(2p) + (X_p-X_k)^2)] for k=i or j
    # T_ij = S_ij * [ (alpha_i*alpha_j/p) - (2*alpha_i^2*alpha_j^2/p^2)*(X_i-X_j)^2 ]
    
    # If i == j, X_i - X_j = 0. T_ii = S_ii * alpha_i * alpha_j / p = 1.0 * alpha_i^2 / (2*alpha_i) = alpha_i / 2.0
    # This is correct: <phi_i | -1/2 d^2/dx^2 | phi_i> = alpha_i/2 for a normalized 1D Gaussian.
    if abs(p) < 1e-12: # Should not happen with positive exponents
        return 0.0 
        
    # The full expression derived from S_ij * [alpha_j - 2*alpha_j^2 * (1/(2p) + (X_p-X_j)^2)]
    # which simplifies to S_ij * [ (alpha_i*alpha_j/p) - (2*alpha_i^2*alpha_j^2/p**2)*(X_i-X_j)**2 ]
    # So, T_ij = s_ij_normalized * ( (alpha_i*alpha_j/p) - (2*(alpha_i**2)*(alpha_j**2)/p**2)*(X_i-X_j)**2 )
    term_const = (alpha_i * alpha_j) / p
    term_dist_sq = (2.0 * (alpha_i**2) * (alpha_j**2) / p**2) * (X_i - X_j)**2
    
    return s_ij_normalized * (term_const - term_dist_sq)

def coulomb_element_gauss(alpha_i, X_i, alpha_j, X_j, Z_nuc, X_nuc_pos):
    """ V_ij_k = <phi_i | -Z_k/|x-X_k| | phi_j> for one nucleus k. """
    p = alpha_i + alpha_j
    
    # Handle the case where p is zero (shouldn't happen with positive alphas, but for safety)
    if abs(p) < 1e-12:
        return 0.0
        
    X_p = (alpha_i * X_i + alpha_j * X_j) / p

    norm_i = (2 * alpha_i / np.pi)**0.25
    norm_j = (2 * alpha_j / np.pi)**0.25
    
    K_gaussian_product_factor = norm_i * norm_j * np.exp(-(alpha_i * alpha_j / p) * (X_i - X_j)**2)
    
    t_boys_arg = p * (X_p - X_nuc_pos)**2
    f0_value = boys_function_F0(t_boys_arg)
    
    # Handle the case where p is zero in the sqrt term (shouldn't happen)
    if abs(p) < 1e-12:
         sqrt_term = 0.0
    else:
         sqrt_term = np.sqrt(2 * np.pi / p)
    
    return -Z_nuc * K_gaussian_product_factor * sqrt_term * f0_value

# --- Standalone functions for numerical integration (for testing/validation) ---

def gauss_to_grid(x_grid, X_center, alpha):
    """ Evaluates a single normalized Gaussian basis function on a grid. """
    N_norm = (2 * alpha / np.pi)**0.25
    return N_norm * np.exp(-alpha * (x_grid - X_center)**2)

def numerical_overlap_gauss(x_grid, phi_i_on_grid, phi_j_on_grid):
    """ Numerically computes S_ij = <phi_i | phi_j> using trapezoidal rule. """
    integrand = phi_i_on_grid * phi_j_on_grid
    return np.trapz(integrand, x_grid)

def numerical_kinetic_gauss(x_grid, phi_i_on_grid, phi_j_on_grid):
    """ Numerically computes T_ij = <phi_i | -1/2 d^2/dx^2 | phi_j> using np.gradient and trapezoidal rule. """
    # Second derivative of phi_j using central differences (np.gradient)
    # Note: dx must be scalar for np.gradient if x_grid is not passed.
    # If x_grid is passed, dx is inferred.
    d_phi_j_dx = np.gradient(phi_j_on_grid, x_grid, edge_order=2)
    d2_phi_j_dx2 = np.gradient(d_phi_j_dx, x_grid, edge_order=2)
    
    integrand = phi_i_on_grid * (-0.5) * d2_phi_j_dx2
    return np.trapz(integrand, x_grid)

def numerical_coulomb_gauss(x_grid, phi_i_on_grid, phi_j_on_grid, nuclear_charges, nuclear_positions, epsilon_denom=1e-7):
    """
    Numerically computes V_ij = sum_k <phi_i | -Z_k/|x-X_k| | phi_j> using trapezoidal rule.
    A small epsilon_denom is added to |x-X_k| to prevent division by zero if a grid point
    coincides with a nuclear position. This makes it an approximation of the true 1/r potential.
    """
    V_ij_sum = 0.0
    phi_product_on_grid = phi_i_on_grid * phi_j_on_grid
    
    for Z_nuc, X_nuc_pos in zip(nuclear_charges, nuclear_positions):
        # Potential V_k(x) = -Z_k / (|x - X_k| + epsilon)
        # Using np.sqrt((x_grid - X_nuc_pos)**2 + epsilon_denom**2) is another way to soften
        # For 1/|x-X_k|, we add epsilon to the denominator directly.
        potential_values_k = -Z_nuc / (np.abs(x_grid - X_nuc_pos) + epsilon_denom)
        
        integrand_k = phi_product_on_grid * potential_values_k
        V_ij_sum += np.trapz(integrand_k, x_grid)
            
    return V_ij_sum


class QuantumSolver1D:
    """
    A simple 1D quantum solver based on the LCAO method with Gaussian basis functions.
    
    This solver calculates the molecular orbitals and their energies for a given set of
    1D atomic nuclei. It solves the generalized eigenvalue problem HC = SCE. The formulas
    used are based on standard integrals for s-type Gaussian orbitals, commonly
    used in 3D, which serve as a good illustrative model for our 1D case.
    """

    def __init__(self, atoms_config):
        """
        Initializes the solver with the system's configuration.

        Args:
            atoms_config (dict): A dictionary containing the system's properties.
                - 'positions': A NumPy array of nuclear positions (X_k).
                - 'charges': A NumPy array of nuclear charges (Z_k).
                - 'exponents': A NumPy array of Gaussian exponents for the basis functions (alpha_i).
        """
        self.positions = atoms_config['positions']
        self.charges = atoms_config['charges']
        self.exponents = atoms_config['exponents']
        self.num_basis_functions = len(self.positions)

        # Matrices to be computed by the solver
        self.s_matrix = np.zeros((self.num_basis_functions, self.num_basis_functions))
        self.t_matrix = np.zeros((self.num_basis_functions, self.num_basis_functions))
        self.v_matrix = np.zeros((self.num_basis_functions, self.num_basis_functions))
        self.h_matrix = None
        
        # Results
        self.energies = None
        self.coefficients = None
        
    def _calculate_matrices(self):
        """
        Calculates the Overlap (S), Kinetic (T), and Potential (V) matrices.
        """
        n = self.num_basis_functions

        # Iterate over all pairs of basis functions (i, j)
        for i in range(n):
            for j in range(n):
                # Parameters for basis functions i and j
                alpha_i, X_i = self.exponents[i], self.positions[i]
                alpha_j, X_j = self.exponents[j], self.positions[j]

                # --- Overlap Integral S_ij ---
                s_ij = overlap_element_gauss(alpha_i, X_i, alpha_j, X_j)
                self.s_matrix[i, j] = s_ij

                # --- Kinetic Energy Integral T_ij ---
                self.t_matrix[i, j] = kinetic_element_gauss(alpha_i, X_i, alpha_j, X_j, s_ij)

                # --- Nuclear Attraction Integral V_ij ---
                v_ij_sum = 0
                for k in range(len(self.charges)):
                    Z_nuc, X_nuc_pos = self.charges[k], self.positions[k]
                    v_ij_k = coulomb_element_gauss(
                        alpha_i, X_i, alpha_j, X_j, Z_nuc, X_nuc_pos
                    )
                    v_ij_sum += v_ij_k

                self.v_matrix[i, j] = v_ij_sum

        # Assemble the full Hamiltonian matrix
        self.h_matrix = self.t_matrix + self.v_matrix

    def solve(self):
        """
        Solves the generalized eigenvalue problem HC = SCE.
        
        The method calculates the matrices if they haven't been computed yet,
        then solves the equation to find orbital energies and coefficients.
        """
        self._calculate_matrices()

        # Solve the generalized eigenvalue problem using scipy.linalg.eigh,
        # which is optimized for Hermitian matrices (H and S are).
        self.energies, self.coefficients = eigh(self.h_matrix, self.s_matrix)

        return self.energies, self.coefficients

    def get_molecular_orbital(self, x_grid, mo_index):
        """
        Constructs the wavefunction of a specific molecular orbital on a given grid.

        Args:
            x_grid (np.ndarray): The 1D grid of points to evaluate the MO on.
            mo_index (int): The index of the molecular orbital to construct.

        Returns:
            np.ndarray: The values of the molecular orbital wavefunction on the grid.
        """
        if self.coefficients is None:
            raise RuntimeError("You must run the .solve() method before constructing an MO.")

        mo_coeffs = self.coefficients[:, mo_index]
        psi = np.zeros_like(x_grid, dtype=float)

        for i in range(self.num_basis_functions):
            c_i = mo_coeffs[i]
            alpha_i, X_i = self.exponents[i], self.positions[i]
            
            # Normalization constant for the i-th basis function
            N_i = (2 * alpha_i / np.pi)**0.25
            
            # Add the contribution of this basis function to the molecular orbital
            basis_function = N_i * np.exp(-alpha_i * (x_grid - X_i)**2)
            psi += c_i * basis_function
            
        return psi