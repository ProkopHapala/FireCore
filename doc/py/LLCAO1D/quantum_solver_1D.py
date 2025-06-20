# /home/prokophapala/git/FireCore/doc/py/LLCAO1D/quantum_solver_1D_new.py

import numpy as np
from scipy.linalg import eigh
from scipy.special import erf
# trun of line-brakes when printing numpy arrays
np.set_printoptions(linewidth=1000)

# --- Core Analytical Integral Functions (using alpha) ---

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
    if abs(p) < 1e-12: return 0.0 # Avoid division by zero if p is too small
    norm_i = (2 * alpha_i / np.pi)**0.25
    norm_j = (2 * alpha_j / np.pi)**0.25
    exp_term = np.exp(-(alpha_i * alpha_j / p) * (X_i - X_j)**2)
    integral_term = (np.pi / p)**0.5
    return norm_i * norm_j * exp_term * integral_term

def kinetic_element_gauss(alpha_i, X_i, alpha_j, X_j, s_ij_normalized):
    """ T_ij = <phi_i | -1/2 d^2/dx^2 | phi_j> for normalized Gaussians. """
    p = alpha_i + alpha_j
    if abs(p) < 1e-12: return 0.0
        
    term_const = (alpha_i * alpha_j) / p
    term_dist_sq = (2.0 * (alpha_i**2) * (alpha_j**2) / p**2) * (X_i - X_j)**2
    
    return s_ij_normalized * (term_const - term_dist_sq)

def coulomb_element_gauss(alpha_i, X_i, alpha_j, X_j, Z_nuc, X_nuc_pos):
    """ V_ij_k = <phi_i | -Z_k/|x-X_k| | phi_j> for one nucleus k. """
    p = alpha_i + alpha_j
    
    if abs(p) < 1e-12: return 0.0
        
    X_p = (alpha_i * X_i + alpha_j * X_j) / p

    norm_i = (2 * alpha_i / np.pi)**0.25
    norm_j = (2 * alpha_j / np.pi)**0.25
    
    K_gaussian_product_factor = norm_i * norm_j * np.exp(-(alpha_i * alpha_j / p) * (X_i - X_j)**2)
    
    t_boys_arg = p * (X_p - X_nuc_pos)**2
    f0_value = boys_function_F0(t_boys_arg)
    
    sqrt_term = np.sqrt(2 * np.pi / p) if abs(p) > 1e-12 else 0.0
    
    return -Z_nuc * K_gaussian_product_factor * sqrt_term * f0_value

def gauss_to_grid(x_grid, X_center, alpha):
    """ Evaluates a single normalized Gaussian basis function on a grid using alpha. """
    N_norm = (2 * alpha / np.pi)**0.25
    return N_norm * np.exp(-alpha * (x_grid - X_center)**2)

# --- Callback Wrappers for Gaussian Basis (using width 'w') ---

def gaussian_overlap_wrapper(w_i, X_i, w_j, X_j):
    return overlap_element_gauss(1/w_i**2, X_i, 1/w_j**2, X_j)

def gaussian_kinetic_wrapper(w_i, X_i, w_j, X_j, s_ij_normalized):
    return kinetic_element_gauss(1/w_i**2, X_i, 1/w_j**2, X_j, s_ij_normalized)

def gaussian_coulomb_wrapper(w_i, X_i, w_j, X_j, Z_nuc, X_nuc_pos):
    return coulomb_element_gauss(1/w_i**2, X_i, 1/w_j**2, X_j, Z_nuc, X_nuc_pos)

def gaussian_basis_eval_wrapper(w, X_center, x_grid):
    return gauss_to_grid(x_grid, X_center, 1/w**2)


# --- Standalone functions for numerical integration (for testing/validation) ---

def numerical_overlap_gauss(x_grid, phi_i_on_grid, phi_j_on_grid):
    """ Numerically computes S_ij = <phi_i | phi_j> using trapezoidal rule. """
    integrand = phi_i_on_grid * phi_j_on_grid
    return np.trapz(integrand, x_grid)

def numerical_kinetic_gauss(x_grid, phi_i_on_grid, phi_j_on_grid):
    """ Numerically computes T_ij = <phi_i | -1/2 d^2/dx^2 | phi_j> using np.gradient and trapezoidal rule. """
    d_phi_j_dx = np.gradient(phi_j_on_grid, x_grid, edge_order=2)
    d2_phi_j_dx2 = np.gradient(d_phi_j_dx, x_grid, edge_order=2)
    integrand = phi_i_on_grid * (-0.5) * d2_phi_j_dx2
    return np.trapz(integrand, x_grid)

def numerical_coulomb_gauss(x_grid, phi_i_on_grid, phi_j_on_grid, nuclear_charges, nuclear_positions, epsilon_denom=1e-7):
    """
    Numerically computes V_ij = sum_k <phi_i | -Z_k/|x-X_k| | phi_j> using trapezoidal rule.
    """
    V_ij_sum = 0.0
    phi_product_on_grid = phi_i_on_grid * phi_j_on_grid
    for Z_nuc, X_nuc_pos in zip(nuclear_charges, nuclear_positions):
        potential_values_k = -Z_nuc / (np.abs(x_grid - X_nuc_pos) + epsilon_denom)
        integrand_k = phi_product_on_grid * potential_values_k
        V_ij_sum += np.trapz(integrand_k, x_grid)
    return V_ij_sum


class QuantumSolver1D:
    """
    A simple 1D quantum solver based on the LCAO method, configurable via callbacks.
    """

    def __init__(self, atoms_config, overlap_fn, kinetic_fn, coulomb_fn, basis_eval_fn):
        """
        Initializes the solver with the system's configuration and integral callbacks.

        Args:
            atoms_config (dict): Contains system properties:
                - 'positions': NumPy array of nuclear positions (X_k).
                - 'charges': NumPy array of nuclear charges (Z_k).
                - 'widths': NumPy array of basis function widths (w_i).
                            (Basis functions are assumed to be centered on nuclear positions).
            overlap_fn (callable): fn(w_i, X_i, w_j, X_j) -> float
            kinetic_fn (callable): fn(w_i, X_i, w_j, X_j, s_ij_normalized) -> float
            coulomb_fn (callable): fn(w_i, X_i, w_j, X_j, Z_nuc, X_nuc_pos) -> float
            basis_eval_fn (callable): fn(w, X_center, x_grid) -> np.ndarray
        """
        # Nuclear properties
        self.nuclear_positions = atoms_config['positions']
        self.nuclear_charges = atoms_config['charges']
        
        # Basis function properties (assuming one basis function per nucleus, centered on it)
        self.basis_widths = atoms_config['widths']
        self.basis_centers = atoms_config['positions'] # Basis functions centered on nuclei
        self.num_basis_functions = len(self.basis_centers)

        if len(self.basis_widths) != self.num_basis_functions:
            raise ValueError("Number of basis widths must match number of basis centers/nuclear positions.")

        # Store callbacks
        self.overlap_fn = overlap_fn
        self.kinetic_fn = kinetic_fn
        self.coulomb_fn = coulomb_fn
        self.basis_eval_fn = basis_eval_fn

        # Matrices to be computed
        self.s_matrix = np.zeros((self.num_basis_functions, self.num_basis_functions))
        self.t_matrix = np.zeros((self.num_basis_functions, self.num_basis_functions))
        self.v_matrix = np.zeros((self.num_basis_functions, self.num_basis_functions))
        self.h_matrix = None
        
        # Results
        self.energies = None
        self.coefficients = None
        
    def _calculate_matrices(self):
        """
        Calculates the Overlap (S), Kinetic (T), and Potential (V) matrices using callbacks.
        """
        n = self.num_basis_functions

        for i in range(n):
            w_i, X_i = self.basis_widths[i], self.basis_centers[i]
            for j in range(n):
                w_j, X_j = self.basis_widths[j], self.basis_centers[j]

                # --- Overlap Integral S_ij ---
                s_ij = self.overlap_fn(w_i, X_i, w_j, X_j)
                self.s_matrix[i, j] = s_ij

                # --- Kinetic Energy Integral T_ij ---
                # The s_ij passed to kinetic_fn should be the normalized overlap,
                # which our overlap_fn already returns for normalized Gaussians.
                self.t_matrix[i, j] = self.kinetic_fn(w_i, X_i, w_j, X_j, s_ij)

                # --- Nuclear Attraction Integral V_ij ---
                v_ij_sum = 0
                for k_nuc in range(len(self.nuclear_charges)):
                    Z_nuc_k = self.nuclear_charges[k_nuc]
                    X_nuc_k_pos = self.nuclear_positions[k_nuc]
                    
                    v_ij_k = self.coulomb_fn( w_i, X_i, w_j, X_j, Z_nuc_k, X_nuc_k_pos )
                    v_ij_sum += v_ij_k
                self.v_matrix[i, j] = v_ij_sum

        self.h_matrix = self.t_matrix + self.v_matrix

    def solve(self):
        """ Solves the generalized eigenvalue problem HC = SCE. """
        self._calculate_matrices()
        self.energies, self.coefficients = eigh(self.h_matrix, self.s_matrix)
        return self.energies, self.coefficients

    def get_molecular_orbital(self, x_grid, mo_index):
        """
        Constructs the wavefunction of a specific MO on a grid using the basis_eval_fn callback.
        """
        if self.coefficients is None:
            raise RuntimeError("Run .solve() before constructing an MO.")

        mo_coeffs_vector = self.coefficients[:, mo_index]
        psi_on_grid = np.zeros_like(x_grid, dtype=float)

        for i in range(self.num_basis_functions):
            coeff_i = mo_coeffs_vector[i]
            w_i = self.basis_widths[i]
            X_center_i = self.basis_centers[i]
            
            basis_function_i_on_grid = self.basis_eval_fn(w_i, X_center_i, x_grid)
            psi_on_grid += coeff_i * basis_function_i_on_grid
            
        return psi_on_grid