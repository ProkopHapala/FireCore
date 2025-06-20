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
    Supports standard diagonalization and iterative localized orbital optimization.
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

        # Results for standard solver
        self.energies = None
        self.coefficients = None

        # Parameters for localized solver
        self.localization_centers = None
        self.localization_cutoff = None

        # Results for localized solver
        self.localized_energies = None
        self.localized_coefficients = None

    def _calculate_matrices(self, bDoOverlap=True, bDoKinetic=True, bDoPotential=True):
        """
        Calculates the Overlap (S), Kinetic (T), and Potential (V) matrices using callbacks.
        """
        n = self.num_basis_functions

        for i in range(n):
            w_i, X_i = self.basis_widths[i], self.basis_centers[i]
            for j in range(i+1):
                w_j, X_j = self.basis_widths[j], self.basis_centers[j]

                # --- Overlap Integral S_ij ---
                if bDoOverlap:
                    s_ij = self.overlap_fn(w_i, X_i, w_j, X_j)
                    self.s_matrix[i, j] = s_ij
                    self.s_matrix[j, i] = s_ij


                # --- Kinetic Energy Integral T_ij ---
                # The s_ij passed to kinetic_fn should be the normalized overlap,
                # which our overlap_fn already returns for normalized Gaussians.
                if bDoKinetic:
                    # Need to ensure s_ij is calculated even if bDoOverlap is False for kinetic_fn
                    if not bDoOverlap and self.s_matrix[i,j] == 0: # Recalculate if needed and not already done
                         s_ij = self.overlap_fn(w_i, X_i, w_j, X_j)
                    else:
                         s_ij = self.s_matrix[i,j]

                    t_ij = self.kinetic_fn(w_i, X_i, w_j, X_j, s_ij)
                    self.t_matrix[i, j] = t_ij
                    self.t_matrix[j, i] = t_ij

                if bDoPotential:
                    # --- Nuclear Attraction Integral V_ij ---
                    v_ij = 0
                    for k_nuc in range(len(self.nuclear_charges)):
                        Z_nuc_k = self.nuclear_charges[k_nuc]
                        X_nuc_k_pos = self.nuclear_positions[k_nuc]

                        v_ij_k = self.coulomb_fn( w_i, X_i, w_j, X_j, Z_nuc_k, X_nuc_k_pos )
                        v_ij += v_ij_k

                    self.v_matrix[i, j] = v_ij
                    self.v_matrix[j, i] = v_ij

        self.h_matrix = self.t_matrix + self.v_matrix

    def solve(self):
        """ Solves the generalized eigenvalue problem HC = SCE. """
        self._calculate_matrices()
        self.energies, self.coefficients = eigh(self.h_matrix, self.s_matrix)
        return self.energies, self.coefficients

    # --- Helper functions for matrix operations (designed for potential sparse adaptation) ---

    def _calculate_mo_overlap_matrix(self, C):
        """
        Calculates the molecular orbital overlap matrix S_MO = C^T S C.
        Can be adapted for sparse C and S later.
        """
        # Dense implementation:
        return C.T @ self.s_matrix @ C

    def _calculate_hamiltonian_action_on_mos(self, C):
        """
        Calculates the action of the Hamiltonian on the MO coefficients: G = H C.
        This is used for the energy gradient. Can be adapted for sparse H and C later.
        """
        # Dense implementation:
        return self.h_matrix @ C

    def _calculate_orbital_norm_sq(self, c_vec):
        """
        Calculates the norm squared of a single orbital vector: c_i^T S c_i.
        Can be adapted for sparse c_vec and S later.
        """
        # Dense implementation:
        return c_vec.T @ self.s_matrix @ c_vec

    def _calculate_orbital_energy_numerator(self, c_vec):
        """
        Calculates the numerator of the energy for a single orbital: c_i^T H c_i.
        Can be adapted for sparse c_vec and H later.
        """
        # Dense implementation:
        return c_vec.T @ self.h_matrix @ c_vec

    # --- Modular steps for the localized solver optimization loop ---

    def _initialize_localized_coefficients(self, num_mos, localization_centers, localization_cutoff):
        """
        Initializes localized coefficients, applies cutoff, and normalizes.
        Assumes one MO per basis function, centered at the basis function's position.
        """
        C = np.eye(self.num_basis_functions, num_mos)

        for i in range(num_mos):
            # Apply localization cutoff
            for mu in range(self.num_basis_functions):
                if np.linalg.norm(self.basis_centers[mu] - localization_centers[i]) > localization_cutoff:
                    C[mu, i] = 0.0
            # Normalize
            norm_sq = self._calculate_orbital_norm_sq(C[:, i])
            if norm_sq > 1e-12:
                 C[:, i] /= np.sqrt(norm_sq)
            else:
                 print(f"Warning: Initial orbital {i} zeroed out by localization cutoff.")
                 pass # Leave as zero
        return C

    def _apply_localization_cutoff(self, C, localization_centers, localization_cutoff):
        """ Applies the hard localization cutoff to the coefficient matrix C. """
        num_mos = C.shape[1]
        for i in range(num_mos):
            for mu in range(self.num_basis_functions):
                 if np.linalg.norm(self.basis_centers[mu] - localization_centers[i]) > localization_cutoff:
                     C[mu, i] = 0.0
        return C

    def _calculate_energy_gradient_update(self, C, loc_step_size, use_lagrange=False, alpha_loc=0.0, localization_centers=None):
        """
        Calculates the gradient-descent update for the coefficient matrix.

        Parameters
        ----------
        C : np.ndarray
            Current coefficient matrix  (num_basis × num_mos)
        loc_step_size : float
            Gradient-descent step size.
        use_lagrange : bool, optional
            If True include Lagrange multipliers so that the unconstrained
            gradient is  G = H C − S C Λ  with  Λ_ii = ⟨c_i|H|c_i⟩/⟨c_i|S|c_i⟩.
        alpha_loc : float, optional
            Strength of quadratic localization potential. 0 → off.
        localization_centers : np.ndarray, optional
            Centres (num_mos) used for the localization potential.
        """
        # --- electronic part
        G = self._calculate_hamiltonian_action_on_mos(C)   # H C

        if use_lagrange:
            # Λ_ii = E_i  with  E_i = (c_i^T H c_i)/(c_i^T S c_i)
            energies = np.zeros(C.shape[1])
            for i in range(C.shape[1]):
                denom = self._calculate_orbital_norm_sq(C[:, i])
                energies[i] = self._calculate_orbital_energy_numerator(C[:, i]) / max(denom, 1e-12)
            G -= self.s_matrix @ (C * energies)  # S C Λ where Λ stored along columns via broadcasting

        # --- localization potential  (quadratic in distance)
        if alpha_loc > 0.0 and localization_centers is not None:
            # For each basis function μ we pre-compute distance^2 to each centre i
            dist2 = (self.basis_centers[:, None] - localization_centers[None, :])**2  # shape (μ,i)
            G += 2.0 * alpha_loc * C * dist2   # derivative of α (r-R_i)^2 |c|^2 w.r.t c

        # gradient descent step
        return -loc_step_size * G

    def _perform_iterative_orthogonalization(self, C, ortho_damp, ortho_iter):
        """ Performs iterative Jacobi-like orthogonalization on C. """
        num_mos = C.shape[1]
        max_off_diag_S = 0.0 # Initialize max off-diagonal S_MO
        for o_iter in range(ortho_iter):
             S_MO = self._calculate_mo_overlap_matrix(C)
             max_off_diag_S = np.max(np.abs(S_MO - np.diag(np.diag(S_MO))))

             dC_ortho = np.zeros_like(C)
             for i in range(num_mos):
                 for j in range(num_mos):
                     if i == j: continue
                     Sij = S_MO[i, j]
                     dC_ortho[:, i] -= C[:, j] * Sij * ortho_damp
             C += dC_ortho
        return C, max_off_diag_S

    def _normalize_coefficients(self, C):
        """ Normalizes each column (orbital) in C. """
        num_mos = C.shape[1]
        for i in range(num_mos):
            norm_sq = self._calculate_orbital_norm_sq(C[:, i])
            if norm_sq > 1e-12:
                C[:, i] /= np.sqrt(norm_sq)
            else:
                 print(f"Warning: Orbital {i} collapsed during normalization.")
                 pass # Leave it potentially zero
        return C

    def _calculate_localized_orbital_energies(self, C):
        """ Calculates the energy for each localized orbital. """
        num_mos = C.shape[1]
        current_energies = np.zeros(num_mos)
        for i in range(num_mos):
             norm_sq = self._calculate_orbital_norm_sq(C[:, i])
             if norm_sq > 1e-12:
                 current_energies[i] = self._calculate_orbital_energy_numerator(C[:, i]) / norm_sq
             else:
                 current_energies[i] = np.nan # Indicate collapse
        return current_energies

    # --- Main localized solver method ---

    def solve_localized(self, localization_centers, localization_cutoff,
                        max_loc_iter=200, loc_energy_tol=1e-6, loc_coeff_tol=1e-6, loc_overlap_tol=1e-6,
                        loc_step_size=0.1, ortho_damp=0.5, ortho_iter=5,
                        use_lagrange=False, alpha_loc=0.0):
        """
        Solves for localized molecular orbitals using iterative energy minimization
        with normalization, orthogonalization, and localization constraints.

        Args:
            localization_centers (np.ndarray): Positions [num_mos] for localizing each MO.
                                               Assumed to match basis_centers for simplicity initially.
            localization_cutoff (float): Maximum distance from localization_center a basis function
                                         can be to have a non-zero coefficient.
            max_loc_iter (int): Maximum iterations for the outer optimization loop.
            loc_energy_tol (float): Convergence tolerance for energy change (not fully implemented yet).
            loc_coeff_tol (float): Convergence tolerance for max coefficient change per orbital.
            loc_overlap_tol (float): Convergence tolerance for max off-diagonal MO overlap.
            loc_step_size (float): Step size for the energy minimization gradient update.
            ortho_damp (float): Damping factor for the iterative orthogonalization step.
            ortho_iter (int): Number of internal iterations for orthogonalization per outer loop step.

        Returns:
            tuple: (localized_energies, localized_coefficients)
        """
        self._calculate_matrices() # Ensure H and S are computed

        num_mos = len(localization_centers) # Number of localized MOs
        if num_mos > self.num_basis_functions:
             raise ValueError("Number of localized orbitals cannot exceed number of basis functions.")
        # Enforce the assumption that localization centers match basis centers for this simple implementation
        if len(localization_centers) != self.num_basis_functions or not np.allclose(localization_centers, self.basis_centers):
             raise NotImplementedError("Current localized solver requires localization_centers to match basis_centers.")

        self.localization_centers = localization_centers
        self.localization_cutoff = localization_cutoff

        # --- Initialization ---
        C = self._initialize_localized_coefficients(num_mos, localization_centers, localization_cutoff)

        print(f"Starting localized solver with R_loc={localization_cutoff}, max_iter={max_loc_iter}, step_size={loc_step_size}, ortho_damp={ortho_damp}, ortho_iter={ortho_iter}")

        # Optimization loop
        for iter in range(max_loc_iter):
            C_old = C.copy()

            # --- Optimization Steps ---
            gradient_update_matrix = self._calculate_energy_gradient_update(C, loc_step_size, use_lagrange=use_lagrange, alpha_loc=alpha_loc, localization_centers=localization_centers)
            C += gradient_update_matrix
            C, max_off_diag_S = self._perform_iterative_orthogonalization(C, ortho_damp, ortho_iter)
            C = self._normalize_coefficients(C)
            C = self._apply_localization_cutoff(C, localization_centers, localization_cutoff)
            C = self._normalize_coefficients(C) # Re-normalize after cutoff

            # --- Calculate current energies and check for collapse ---
            current_energies = self._calculate_localized_orbital_energies(C)

            # --- Convergence Check ---
            coeff_change = np.max(np.abs(C - C_old))
            # Recalculate S_MO for final convergence check after all steps
            S_MO_final = self._calculate_mo_overlap_matrix(C)
            max_off_diag_S_final = np.max(np.abs(S_MO_final - np.diag(np.diag(S_MO_final))))
            
            # Calculate total energy (sum of occupied orbital energies)
            # Assuming all num_mos are occupied for simplicity here.
            # If only a subset were occupied, you'd sum only those.
            total_energy = np.nansum(current_energies) # Use nansum to handle potential NaN from collapsed orbitals

            # Calculate norm of optimization forces (gradient G = H C, before scaling by step_size)
            # The gradient_update_matrix is -loc_step_size * G, so G = -gradient_update_matrix / loc_step_size
            G_matrix = -gradient_update_matrix / loc_step_size if loc_step_size != 0 else np.zeros_like(gradient_update_matrix)
            norm_of_forces = np.linalg.norm(G_matrix)

            print(f"Iter {iter}: E_tot = {total_energy:.6f}, |Forces| = {norm_of_forces:.6f}, Max Coeff Chg = {coeff_change:.6f}, Max S_MO_offdiag = {max_off_diag_S_final:.6f}")

            if coeff_change < loc_coeff_tol and max_off_diag_S_final < loc_overlap_tol: # Add energy_change < loc_energy_tol if tracked
                 print("Localized solver converged.")
                 break

        # Store final results
        self.localized_coefficients = C
        self.localized_energies = current_energies # Store final energies

        # Note: Localized orbitals are not typically sorted by energy,
        # as their identity is tied to the localization center.

        return self.localized_energies, self.localized_coefficients

    def get_molecular_orbital(self, x_grid, mo_index):
        """
        Constructs the wavefunction of a specific MO on a grid using the basis_eval_fn callback.
        Uses coefficients from the standard diagonalization solver.
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

    def get_localized_molecular_orbital(self, x_grid, mo_index):
        """
        Constructs the wavefunction of a specific LOCALIZED MO on a grid
        using the basis_eval_fn callback and localized coefficients.
        """
        if self.localized_coefficients is None:
            raise RuntimeError("Run .solve_localized() before constructing a localized MO.")

        mo_coeffs_vector = self.localized_coefficients[:, mo_index]
        psi_on_grid = np.zeros_like(x_grid, dtype=float)

        for i in range(self.num_basis_functions):
            coeff_i = mo_coeffs_vector[i]
            w_i = self.basis_widths[i]
            X_center_i = self.basis_centers[i]
            basis_function_i_on_grid = self.basis_eval_fn(w_i, X_center_i, x_grid)
            psi_on_grid += coeff_i * basis_function_i_on_grid

        return psi_on_grid