#!/usr/bin/env python3
"""
Finds an optimal analytical basis set along the z-direction for a given set of
1D potential profiles using Singular Value Decomposition (SVD).

Method Outline:
1. Generate/Provide Sample Functions:
   - A set of M sample functions y_j(z), each defined on Nz points.
   - These can be generated from physical models (e.g., atomic potentials)
     or provided directly.
   - Represented as a matrix Y_T (M, Nz).

2. Generate/Provide Library Basis:
   - A library of P analytical basis functions phi_p(z).
   - Default: polynomials z^n. User can provide custom functions or a pre-computed matrix.
   - Represented as a matrix Phi_T (P, Nz).
   - Optionally, z-coordinates can be scaled for numerical stability (e.g., for polynomials).
   - Optionally, this basis can be orthogonalized (e.g., via Gram-Schmidt).

3. Represent Samples in Library Basis:
   - Express each sample function y_j(z) as a linear combination of phi_p(z):
     y_j(z) approx sum_p S_pj * phi_p(z)
   - Solve for the coefficient matrix S (P, M) using (weighted) least squares:
     Y_T.T approx Phi_T.T @ S  => S = lstsq(Phi_T.T, Y_T.T)[0]
   - Weights w(z) can be applied to emphasize certain regions of z.

4. Find Optimal Basis via SVD:
   - Perform SVD on the coefficient matrix S: S = U_svd @ Sigma @ Vh_svd.
   - The first K left singular vectors (columns of U_svd) define the K optimal
     linear combinations of the original library basis functions.
   - U_k = U_svd[:, :K]  (P, K) are the coefficients of the new optimal basis
     functions in terms of the old library basis.

5. Construct Optimal Basis:
   - The K new optimal basis functions B_opt_k(z) are:
     B_opt_k(z) = sum_p (U_k)_pk * phi_p(z)
   - Evaluated matrix: B_opt_T = U_k.T @ Phi_T  (K, Nz).
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import lstsq, svd
from scipy.interpolate import interp1d
import random

# Attempt to import from FoldedAtomicFunctions for sample generation
try:
    from FoldedAtomicFunctions import GridManager, PotentialCalculator
    FAF_AVAILABLE = True
except ImportError:
    FAF_AVAILABLE = False
    print("Warning: FoldedAtomicFunctions.py not found. Sample generation from atomic potentials will not be available.")
    GridManager, PotentialCalculator = None, None


# --- 1. Sample Function Generation ---

def define_system_configurations(
    num_systems: int,
    param_ranges: dict,
    fixed_params: dict = None
) -> list:
    """
    Generates a list of system configurations with randomized parameters.

    Parameters:
    -----------
    num_systems : int
        Number of system configurations to generate.
    param_ranges : dict
        Dictionary specifying parameter ranges. Example:
        {
            'L_x': (8.0, 12.0),  # Range for lattice constant
            'num_atom_types': 2,
            'atom_params': [  # List of dicts, one per atom type
                {'R': (1.5, 2.0), 'E': (0.005, 0.015), 'a': (1.4, 1.8), 'q': (-0.1, 0.1)}, # Type 1
                {'R': (2.0, 2.5), 'E': (0.005, 0.015), 'a': (1.4, 1.8), 'q': (-0.1, 0.1)}, # Type 2
            ],
            'atoms_per_cell_pattern': [[0], [1], [0,1], [1,0]] # Patterns of atom types in cell
                                                              # 0 refers to atom_params[0] etc.
                                                              # x-positions will be spread in cell
        }
    fixed_params : dict, optional
        Dictionary of parameters that are fixed for all systems.
        Example: {'atom_z_pos': 0.0}

    Returns:
    --------
    list_of_system_defs : list
        A list of dictionaries, each defining a system:
        {'L_x': float, 'atoms_specs': [{'x': float, 'z': float, 'q': float,
                                       'R': float, 'E': float, 'a': float}, ...]}
    """
    if fixed_params is None:
        fixed_params = {}
    
    list_of_system_defs = []
    default_atom_z = fixed_params.get('atom_z_pos', 0.0)

    for _ in range(num_systems):
        sys_def = {}
        
        # Lattice constant
        if 'L_x' in param_ranges:
            sys_def['L_x'] = random.uniform(*param_ranges['L_x'])
        elif 'L_x' in fixed_params:
            sys_def['L_x'] = fixed_params['L_x']
        else:
            raise ValueError("L_x must be in param_ranges or fixed_params")

        atoms_specs = []
        if 'atoms_per_cell_pattern' in param_ranges and 'atom_params' in param_ranges:
            pattern_choices = param_ranges['atoms_per_cell_pattern']
            chosen_pattern = random.choice(pattern_choices) # e.g., [0, 1] meaning atom type 0 then type 1
            
            num_atoms_in_cell = len(chosen_pattern)
            x_positions = [(i + 0.5) * sys_def['L_x'] / num_atoms_in_cell for i in range(num_atoms_in_cell)]

            for i, atom_type_idx in enumerate(chosen_pattern):
                atom_p_ranges = param_ranges['atom_params'][atom_type_idx]
                spec = {'x': x_positions[i], 'z': default_atom_z}
                for param_name, prange in atom_p_ranges.items():
                    spec[param_name] = random.uniform(*prange)
                atoms_specs.append(spec)
        
        sys_def['atoms_specs'] = atoms_specs
        list_of_system_defs.append(sys_def)
        
    return list_of_system_defs


def generate_1d_potential_profiles(
    system_definitions: list,
    common_z_coords: np.ndarray,
    x_slice_coords: list,
    potential_type: str = 'morse', # 'morse', 'coulomb', 'total'
    potential_calc_grid_step: float = 0.1,
    num_periodic_images: int = 5
) -> np.ndarray:
    """
    Generates 1D potential profiles along z for multiple systems and x-slices.

    Parameters:
    -----------
    system_definitions : list
        List of system definitions (output of `define_system_configurations` or user-provided).
        Each dict must contain 'L_x' (float) and 'atoms_specs' (list of dicts).
        Atom spec dict: {'x', 'z', 'q', 'R' (vdW radius), 'E' (vdW depth), 'a' (Morse exp)}.
    common_z_coords : np.ndarray
        1D array of z-coordinates (Nz) for the final profiles.
    x_slice_coords : list
        List of x-coordinates within the unit cell to take z-slices.
    potential_type : str
        'morse', 'coulomb', or 'total'.
    potential_calc_grid_step : float
        Grid step for internal 2D potential calculation.
    num_periodic_images : int
        Number of periodic images for PotentialCalculator.

    Returns:
    --------
    Y_T : np.ndarray
        Array of shape (M, Nz), where M = len(system_definitions) * len(x_slice_coords).
        Each row is a 1D potential profile.
    """
    if not FAF_AVAILABLE:
        raise RuntimeError("FoldedAtomicFunctions module is not available. Cannot generate potentials.")

    Nz = len(common_z_coords)
    min_z, max_z = np.min(common_z_coords), np.max(common_z_coords)
    
    all_profiles = []

    for i_sys, sys_def in enumerate(system_definitions):
        L_x = sys_def['L_x']
        atoms_specs = sys_def['atoms_specs']

        # Determine if x_slice_coords are fractional or absolute for the current L_x
        # Heuristic: if all original x_slice_coords are in [0,1), assume fractional.
        # A more robust method might involve an explicit flag.
        current_x_slices_absolute = []
        if x_slice_coords: # Ensure x_slice_coords is not None
            is_fractional_heuristic = all(0 <= x_val_orig < 1.0 for x_val_orig in x_slice_coords)
            if is_fractional_heuristic:
                current_x_slices_absolute = [x_val_orig * L_x for x_val_orig in x_slice_coords]
            else:
                current_x_slices_absolute = list(x_slice_coords) # Use as is (absolute)


        # Setup GridManager for 2D potential calculation
        # cell_z and z_offset are set to precisely cover the common_z_coords range
        gm = GridManager(cell_x=L_x,
                         cell_z=max_z - min_z, # This is the extent of the grid in z
                         grid_step=potential_calc_grid_step,
                         z_offset=min_z) # Grid starts at min_z

        pc = PotentialCalculator(gm)
        for atom_spec in atoms_specs:
            pc.add_atom(x=atom_spec['x'], z=atom_spec.get('z', 0.0), # Atoms typically at z=0 plane
                        charge=atom_spec.get('q', 0.0),
                        vdw_radius=atom_spec['R'],
                        vdw_depth=atom_spec['E'],
                        morse_a=atom_spec['a'])

        # Calculate 2D potential
        V_2D = np.zeros_like(gm.X)
        if potential_type == 'morse' or potential_type == 'total':
            V_2D += pc.calculate_morse_periodic(n_images=num_periodic_images)
        if potential_type == 'coulomb' or potential_type == 'total':
            V_2D += pc.calculate_coulomb_periodic(n_images=num_periodic_images)
        
        # Extract 1D slices and interpolate
        for x_val in current_x_slices_absolute:
            # Find closest x-index in the GridManager's x_coords
            # gm.x_coords are cell centers, X meshgrid is fine
            x_idx = np.argmin(np.abs(gm.x_coords - x_val))
            
            # Profile from 2D grid at this x_idx, along all z points of gm.z_coords
            profile_on_gm_z = V_2D[x_idx, :] 
            
            # Interpolate onto common_z_coords
            # gm.z_coords are the z-coordinates for profile_on_gm_z
            if len(gm.z_coords) > 1 :
                interp_func = interp1d(gm.z_coords, profile_on_gm_z,
                                       kind='cubic', fill_value="extrapolate")
                profile_on_common_z = interp_func(common_z_coords)
            elif len(gm.z_coords) == 1 and len(common_z_coords) == 1 and np.isclose(gm.z_coords[0], common_z_coords[0]):
                 profile_on_common_z = profile_on_gm_z # Single point match
            elif len(gm.z_coords) == 1: # Single point in grid, replicate if common_z_coords is also single point
                profile_on_common_z = np.full_like(common_z_coords, profile_on_gm_z[0])
            else: # Should not happen with reasonable common_z_coords and grid_step
                profile_on_common_z = np.zeros_like(common_z_coords)

            all_profiles.append(profile_on_common_z)
            # print(f"  System {i_sys+1}, x-slice {x_val:.2f}: profile generated.")

    if not all_profiles:
        return np.array([]).reshape(0, Nz)
    return np.array(all_profiles)


# --- 2. Weighting Function ---

def generate_z_weights(
    zs: np.ndarray,
    z0: float,
    z_decay_width: float = 1.0,
    min_weight: float = 0.0
) -> np.ndarray:
    """
    Generates weights for z-coordinates to damp repulsive regions.
    w = 1.0 for z > z0
    w = 1.0 - ((z - z0) / z_decay_width)^2 for z <= z0
    w = max(min_weight, w)

    Parameters:
    -----------
    zs : np.ndarray
        1D array of z-coordinates.
    z0 : float
        Transition point. Potential is typically attractive for z > z0.
    z_decay_width : float
        Controls the steepness of quadratic damping for z < z0.
        The weight becomes zero at z = z0 - z_decay_width.
    min_weight : float
        Minimum allowed weight.

    Returns:
    --------
    ws : np.ndarray
        1D array of weights, same shape as zs.
    """
    ws = np.ones_like(zs, dtype=float)
    mask_repulsive = zs <= z0
    
    # Normalized distance from z0 into the repulsive region
    # (z - z0) is negative or zero here.
    # We want weight to be 1 at z0, and decrease as z becomes smaller.
    # Let term = (zs[mask_repulsive] - z0) / z_decay_width. This is negative, range e.g. [-1, 0] if z_decay_width is appropriate.
    # term^2 is positive, range e.g. [0, 1].
    # 1 - term^2 gives damping.
    normalized_dist_sq = ((zs[mask_repulsive] - z0) / z_decay_width)**2
    ws[mask_repulsive] = 1.0 - normalized_dist_sq
    
    ws = np.maximum(ws, min_weight)
    return ws


# --- 3. Library Basis Set ---

def create_library_basis(
    zs: np.ndarray,
    config: dict
) -> tuple:
    """
    Creates the library basis matrix Phi_T (P, Nz) and scaling info.

    Parameters:
    -----------
    zs : np.ndarray
        1D array of z-coordinates (Nz).
    config : dict
        Configuration for basis generation. Examples:
        {'type': 'polynomial', 'degree': 8, 'scale_z': True}
        {'type': 'custom', 'generator_func': my_func, 'num_functions': P, 'params': {}}
          where my_func(z_array, index, num_functions, params) returns a 1D basis array.
        {'type': 'list_of_funcs', 'functions': [f1, f2, ...]}
          where fi(z_array) returns a 1D basis array.

    Returns:
    --------
    Phi_T : np.ndarray
        Library basis matrix (P, Nz). Each row is a basis function.
    z_scale_info : tuple or None
        (z_min_orig, z_range_orig) if scaling was applied, else None.
    basis_labels : list
        List of string labels for each basis function.
    """
    Nz = len(zs)
    basis_type = config.get('type', 'polynomial')
    z_coords_to_use = np.copy(zs)
    z_scale_info = None
    basis_labels = []

    if config.get('scale_z', False):
        z_min_orig = np.min(zs)
        z_max_orig = np.max(zs)
        z_range_orig = z_max_orig - z_min_orig
        if z_range_orig > 1e-9:
            z_coords_to_use = (zs - z_min_orig) / z_range_orig
            z_scale_info = (z_min_orig, z_range_orig)
        # else: z_coords_to_use remains zs (all points are same, range is zero)

    Phi_list = []

    if basis_type == 'polynomial':
        degree = config.get('degree', 8)
        P = degree + 1
        for p in range(P):
            Phi_list.append(z_coords_to_use**p)
            label = f"z^{p}"
            if z_scale_info: label += "_scaled"
            basis_labels.append(label)
    elif basis_type == 'custom':
        generator_func = config['generator_func']
        P = config['num_functions']
        params = config.get('params', {})
        for i in range(P):
            Phi_list.append(generator_func(z_coords_to_use, i, P, params))
            basis_labels.append(f"custom_func_{i}")
    elif basis_type == 'list_of_funcs':
        custom_functions = config['functions']
        P = len(custom_functions)
        for i, func in enumerate(custom_functions):
            Phi_list.append(func(z_coords_to_use)) # Pass scaled z if scale_z was true
            basis_labels.append(f"provided_func_{i}")
    else:
        raise ValueError(f"Unknown basis type: {basis_type}")

    if not Phi_list:
        return np.array([]).reshape(0, Nz), z_scale_info, []
        
    return np.array(Phi_list), z_scale_info, basis_labels


def orthogonalize_basis_gram_schmidt(
    Phi_T: np.ndarray,
    ws: np.ndarray = None
) -> np.ndarray:
    """
    Orthogonalizes the basis functions (rows of Phi_T) using Gram-Schmidt.
    If weights `ws` are provided, performs weighted Gram-Schmidt.
    The resulting basis Q_T will satisfy sum_i Q_ki Q_li w_i = delta_kl (approx).

    Parameters:
    -----------
    Phi_T : np.ndarray
        Basis matrix (P, Nz), where P is num_basis_funcs, Nz is num_z_points.
        Each row is a basis function.
    ws : np.ndarray, optional
        1D array of weights (Nz).

    Returns:
    --------
    Q_T : np.ndarray
        Orthogonalized basis matrix (P, Nz).
    """
    P, Nz = Phi_T.shape
    if P == 0: return np.array([]).reshape(0,Nz)
        
    Q_T = np.zeros_like(Phi_T)
    
    if ws is None:
        _ws = np.ones(Nz)
    else:
        _ws = np.copy(ws)
        if np.any(_ws < 0):
            print("Warning: Negative weights encountered in Gram-Schmidt. Clamping to 0.")
            _ws = np.maximum(_ws, 0)

    for k in range(P):
        vk = Phi_T[k, :]
        qk_prime = np.copy(vk)
        
        for j in range(k):
            qj = Q_T[j, :]
            # Inner product: <vk, qj>_w = sum(vk * qj * _ws)
            # Projection: (<vk, qj>_w / <qj, qj>_w) * qj
            # Since qj is already normalized (<qj,qj>_w = 1), it's just <vk, qj>_w * qj
            inner_prod_vk_qj = np.sum(vk * qj * _ws)
            qk_prime -= inner_prod_vk_qj * qj
            
        # Normalize qk_prime: ||qk_prime||_w = sqrt(sum(qk_prime^2 * _ws))
        norm_qk_prime = np.sqrt(np.sum(qk_prime**2 * _ws))
        
        if norm_qk_prime > 1e-10: # Avoid division by zero for linearly dependent vectors
            Q_T[k, :] = qk_prime / norm_qk_prime
        else:
            # This basis vector is linearly dependent on previous ones or zero.
            # Q_T[k,:] remains zeros, or handle as error/warning.
            print(f"Warning: Basis function {k} appears linearly dependent or zero during orthogonalization.")
            # Optionally, can try to fill with a random orthogonal vector, but usually indicates issues with input basis.
            
    return Q_T


# --- 4. Coefficient Computation & SVD ---

def compute_coefficients_in_library(
    Y_T: np.ndarray,
    Phi_T: np.ndarray,
    ws: np.ndarray = None,
    is_phi_orthogonal_weighted: bool = False
) -> np.ndarray:
    """
    Computes coefficients S (P, M) for Y_T approx S.T @ Phi_T.
    Equivalent to Y.T approx Phi_T.T @ S.

    Parameters:
    -----------
    Y_T : np.ndarray
        Sample functions (M, Nz). M samples, Nz z-points.
    Phi_T : np.ndarray
        Library basis functions (P, Nz). P basis functions.
    ws : np.ndarray, optional
        Weights (Nz).
    is_phi_orthogonal_weighted : bool
        If True, Phi_T rows are assumed orthonormal w.r.t. weights `ws`.

    Returns:
    --------
    S : np.ndarray
        Coefficient matrix (P, M).
    """
    M, Nz = Y_T.shape
    P = Phi_T.shape[0]

    if M == 0: return np.array([]).reshape(P,0)
    if P == 0: return np.array([]).reshape(0,M)

    if is_phi_orthogonal_weighted:
        # S_ps = sum_i Y_si * Phi_pi * w_i
        # S = (Phi_T * ws) @ Y_T.T
        if ws is None:
            _ws = np.ones(Nz)
        else:
            _ws = ws
        S = (Phi_T * _ws[np.newaxis, :]) @ Y_T.T
    else:
        # Solve Y_T.T approx Phi_T.T @ S using (weighted) least squares.
        # A x = b => Phi_T.T @ S_col = Y_T_col
        A = Phi_T.T # (Nz, P)
        b_matrix = Y_T.T # (Nz, M)
        
        if ws is not None:
            sqrt_w = np.sqrt(np.maximum(ws, 0)) # Ensure non-negative weights
            A_w = A * sqrt_w[:, np.newaxis] # Weighted A
            b_matrix_w = b_matrix * sqrt_w[:, np.newaxis] # Weighted b
            S, residuals, rank, s_vals_lstsq = lstsq(A_w, b_matrix_w)
        else:
            S, residuals, rank, s_vals_lstsq = lstsq(A, b_matrix)
    return S


def find_optimal_analytical_basis_svd(
    S: np.ndarray,
    Phi_T: np.ndarray,
    K: int
) -> tuple:
    """
    Performs SVD on coefficient matrix S to find K optimal basis functions.

    Parameters:
    -----------
    S : np.ndarray
        Coefficient matrix (P, M) from `compute_coefficients_in_library`.
    Phi_T : np.ndarray
        Original library basis matrix (P, Nz).
    K : int
        Desired number of new optimal basis functions.

    Returns:
    --------
    B_opt_T : np.ndarray
        Optimal basis functions (K, Nz). Each row is an optimal basis function.
    U_k_coeffs : np.ndarray
        Coefficients (P, K) defining new basis in terms of Phi_T.
        B_opt_T = U_k_coeffs.T @ Phi_T.
    s_vals_svd : np.ndarray
        Singular values from SVD of S.
    """
    P, M = S.shape
    if P == 0 or M == 0 :
        return np.array([]).reshape(min(K,P), Phi_T.shape[1]), np.array([]).reshape(P,min(K,P)), np.array([])

    U_svd, s_vals_svd, Vh_svd = svd(S, full_matrices=False)
    # U_svd is (P, min(P,M))

    num_available_components = U_svd.shape[1]
    if K > num_available_components:
        print(f"Warning: Requested K={K} optimal basis functions, but only "
              f"{num_available_components} available from SVD. Using K={num_available_components}.")
        K = num_available_components
    
    U_k_coeffs = U_svd[:, :K]  # (P, K)
    
    # B_opt_T = U_k_coeffs.T @ Phi_T
    # (K, P) @ (P, Nz) -> (K, Nz)
    B_opt_T = U_k_coeffs.T @ Phi_T
    
    return B_opt_T, U_k_coeffs, s_vals_svd


# --- 5. Plotting & Analysis Utilities ---

def plot_1d_profiles(
    zs: np.ndarray,
    profiles_T: np.ndarray, # (num_profiles, Nz)
    title: str,
    max_plot: int = 10,
    ws: np.ndarray = None,
    filename: str = None
):
    """Plots 1D profiles."""
    if profiles_T.size == 0:
        print(f"No profiles to plot for: {title}")
        return
        
    num_profiles = profiles_T.shape[0]
    plt.figure(figsize=(10, 6))
    for i in range(min(num_profiles, max_plot)):
        plt.plot(zs, profiles_T[i, :], label=f'Profile {i+1}', alpha=0.7)
    
    if ws is not None:
        # Scale weights to fit nicely on the plot
        ax2 = plt.gca().twinx()
        min_prof = np.min(profiles_T) if profiles_T.size > 0 else 0
        max_prof = np.max(profiles_T) if profiles_T.size > 0 else 1
        # Plot weights such that they are visible, e.g., in top 20% of y-axis range
        plot_ws = min_prof + 0.8 * (max_prof - min_prof) + 0.2 * (max_prof - min_prof) * (ws / np.max(ws) if np.max(ws) > 0 else ws)
        ax2.plot(zs, ws, 'k--', label='Weights (scaled)', alpha=0.5, linewidth=1)
        ax2.set_ylabel('Weights (scaled)')
        ax2.tick_params(axis='y')
        # plt.gca().legend(loc='upper left') # Combine legends if possible or adjust
        # ax2.legend(loc='upper right')


    plt.title(title)
    plt.xlabel('z (Å)')
    plt.ylabel('Potential / Value')
    if num_profiles > 0 : plt.legend()
    plt.grid(True)
    if filename:
        plt.savefig(filename)
        print(f"Plot saved to {filename}")
    plt.show()

def plot_singular_values(s_vals, K_opt, filename: str = None):
    """Plots singular values."""
    if s_vals.size == 0:
        print("No singular values to plot.")
        return
    plt.figure(figsize=(8, 5))
    plt.plot(range(1, len(s_vals) + 1), s_vals, 'o-')
    plt.title('Singular Values from SVD of Sample Coefficients Matrix (S)')
    plt.xlabel('Component Number')
    plt.ylabel('Singular Value')
    if K_opt > 0 and K_opt <= len(s_vals):
      plt.axvline(K_opt, color='r', linestyle='--', label=f'Selected K={K_opt}')
      plt.legend()
    plt.grid(True)
    plt.yscale('log')
    if filename:
        plt.savefig(filename)
        print(f"Plot saved to {filename}")
    plt.show()

def print_analytical_form_polynomial(
    U_k_coeffs: np.ndarray, # (P, K)
    z_scale_info: tuple, # (z_min, z_range)
    basis_labels: list = None, # Labels for the original Phi_T basis
    K_to_print: int = -1
):
    """Prints analytical form if original basis was polynomial and scaled."""
    if U_k_coeffs.size == 0: return

    P, K_actual = U_k_coeffs.shape
    if K_to_print < 0 or K_to_print > K_actual:
        K_to_print = K_actual

    print("\n--- Analytical Expressions for Optimal Basis Functions ---")
    if z_scale_info:
        z_min, z_range = z_scale_info
        if abs(z_range) < 1e-9:
             print(f"Note: z_scaled = (z - {z_min:.3f}) / (near_zero_range) => effectively z_scaled is const or z itself if z_min is also near zero")
        else:
             print(f"Note: z_scaled = (z - {z_min:.3f}) / {z_range:.3f}")
    else:
        print("Note: Original z-coordinates were used (no scaling).")

    for k_idx in range(K_to_print):
        expr_parts = []
        for p_idx in range(P):
            coeff_val = U_k_coeffs[p_idx, k_idx]
            if abs(coeff_val) > 1e-6: # Only include significant terms
                term_label = f"phi_{p_idx}"
                if basis_labels and p_idx < len(basis_labels):
                    term_label = basis_labels[p_idx]
                
                # Check if term_label implies polynomial for z_scaled
                is_poly_term = "z^" in term_label and ("_scaled" in term_label or not z_scale_info)

                if is_poly_term:
                    if "z^0" in term_label: # Constant term
                         expr_parts.append(f"{coeff_val:.3e}")
                    elif "z^1" in term_label and ("_scaled" in term_label or not z_scale_info): # Linear term
                         expr_parts.append(f"({coeff_val:.3e} * {'z_scaled' if z_scale_info else 'z'})")
                    else: # Higher order term
                        power = term_label.split('^')[-1].replace("_scaled","")
                        expr_parts.append(f"({coeff_val:.3e} * {'z_scaled' if z_scale_info else 'z'}^{power})")
                else: # General term from basis_labels
                    expr_parts.append(f"({coeff_val:.3e} * {term_label})")

        print(f"B_opt_{k_idx+1}({'z_scaled' if z_scale_info else 'z'}) = {' + '.join(expr_parts) if expr_parts else '0.0'}")


# --- Main Pipeline Function ---
def run_optimization_pipeline(
    # Z-coordinates
    common_z_coords: np.ndarray,
    # Sample functions: Provide Y_T directly OR generate them
    Y_T: np.ndarray = None,
    system_definitions: list = None, 
    num_systems_to_generate: int = 0,
    system_generation_config: dict = None,
    x_slice_coords_for_generation: list = None,
    # x_slice_coords_for_generation: list, e.g. [0.5] for fractional, or [5.0] for absolute for L_x=10
    potential_type_for_generation: str = 'morse', # 
    potential_calc_grid_step: float = 0.1,
    num_periodic_images_for_generation: int = 5,
    # Library basis: Provide Phi_T directly OR generate it
    Phi_T_library: np.ndarray = None,
    library_basis_labels: list = None, # If Phi_T_library is provided, labels can be too
    library_basis_config: dict = None, # e.g., {'type': 'polynomial', 'degree': 8, 'scale_z': True}
    z_scale_info_provided: tuple = None, # If Phi_T_library is provided and was scaled
    # Orthogonalization
    orthogonalize_library: bool = False,
    # Weighting
    z_weights: np.ndarray = None,
    z_weight_z0: float = None, 
    z_weight_decay_width: float = 1.0,
    z_weight_min_weight: float = 0.0,
    # SVD
    num_optimal_K: int = 4,
    # Plotting
    do_plots: bool = True,
    plot_filename_prefix: str = "opt_basis_"
):
    """
    Main pipeline to find optimal analytical basis.
    """
    print("--- Starting Optimal Basis Search Pipeline ---")

    zs = common_z_coords
    if zs is None or zs.ndim != 1 or zs.size == 0:
        raise ValueError("`common_z_coords` must be a non-empty 1D NumPy array.")
    Nz = len(zs)

    # 0. Define z-weights
    ws = z_weights
    if ws is None and z_weight_z0 is not None:
        print(f"Generating z-weights with z0={z_weight_z0}, decay_width={z_weight_decay_width}")
        ws = generate_z_weights(zs, z_weight_z0, z_weight_decay_width, z_weight_min_weight)
    
    # 1. Generate/Load Sample Functions Y_T (M, Nz)
    _Y_T = Y_T
    if _Y_T is None:
        if system_definitions is None and num_systems_to_generate > 0 and system_generation_config and x_slice_coords_for_generation:
            print(f"Generating {num_systems_to_generate} system configurations...")
            system_definitions = define_system_configurations(num_systems_to_generate, system_generation_config)
        
        if system_definitions and x_slice_coords_for_generation:
            print(f"Generating 1D potential profiles for {len(system_definitions)} systems and {len(x_slice_coords_for_generation)} x-slices...")
            _Y_T = generate_1d_potential_profiles(
                system_definitions, zs, x_slice_coords_for_generation,
                potential_type=potential_type_for_generation,
                potential_calc_grid_step=potential_calc_grid_step,
                num_periodic_images=num_periodic_images_for_generation
            )
            if _Y_T.size == 0:
                 raise ValueError("Sample function generation resulted in an empty set.")
            print(f"Generated Y_T with shape: {_Y_T.shape}")
        else:
            raise ValueError("No sample functions (Y_T) provided and cannot generate them with current parameters.")
    else:
        print(f"Using provided Y_T with shape: {_Y_T.shape}")

    if _Y_T.shape[1] != Nz:
        raise ValueError(f"Y_T shape {_Y_T.shape} inconsistent with common_z_coords (Nz={Nz})")

    if do_plots:
        plot_1d_profiles(zs, _Y_T, "Sample Functions (Y_T)", ws=ws, filename=f"{plot_filename_prefix}sample_functions.png" if plot_filename_prefix else None)

    # 2. Generate/Load Library Basis Phi_T (P, Nz)
    _Phi_T = Phi_T_library
    _z_scale_info = z_scale_info_provided
    _basis_labels = library_basis_labels if library_basis_labels else []

    if _Phi_T is None:
        if library_basis_config:
            print(f"Creating library basis with config: {library_basis_config}")
            _Phi_T, _z_scale_info, _basis_labels = create_library_basis(zs, library_basis_config)
            if _Phi_T.size == 0:
                raise ValueError("Library basis generation resulted in an empty set.")
            print(f"Generated Phi_T with shape: {_Phi_T.shape}")
        else:
            raise ValueError("No library basis (Phi_T) provided and no generation config given.")
    else:
        print(f"Using provided Phi_T with shape: {_Phi_T.shape}")
        if not _basis_labels: # Create generic labels if none provided for user's Phi_T
            _basis_labels = [f"provided_phi_{i}" for i in range(_Phi_T.shape[0])]


    if _Phi_T.shape[1] != Nz:
        raise ValueError(f"Phi_T shape {_Phi_T.shape} inconsistent with common_z_coords (Nz={Nz})")
    
    P_lib = _Phi_T.shape[0]
    if P_lib == 0:
        raise ValueError("Library basis Phi_T is empty.")

    is_phi_orthogonal_weighted = False
    if orthogonalize_library:
        print("Orthogonalizing library basis Phi_T using Gram-Schmidt...")
        _Phi_T = orthogonalize_basis_gram_schmidt(_Phi_T, ws)
        is_phi_orthogonal_weighted = True # After GS, it's ortho w.r.t weights used (or unweighted if ws=None)
        _basis_labels = [f"ortho_phi_{i}" for i in range(P_lib)] # Labels change after ortho
        _z_scale_info = None # Orthogonalization might change simple polynomial interpretation
        print(f"Orthogonalized Phi_T shape: {_Phi_T.shape}")


    if do_plots:
        plot_1d_profiles(zs, _Phi_T, "Library Basis Functions (Phi_T)", max_plot=P_lib, filename=f"{plot_filename_prefix}library_basis.png" if plot_filename_prefix else None)

    # 3. Compute Coefficients S (P, M)
    print("Computing coefficients S of samples in library basis...")
    S_coeffs = compute_coefficients_in_library(_Y_T, _Phi_T, ws, is_phi_orthogonal_weighted)
    print(f"Computed S_coeffs matrix with shape: {S_coeffs.shape}")

    # 4. Find Optimal Basis via SVD
    print(f"Finding {num_optimal_K} optimal basis functions via SVD...")
    B_opt_T, U_k_coeffs, s_vals_svd = find_optimal_analytical_basis_svd(
        S_coeffs, _Phi_T, num_optimal_K
    )
    print(f"Found optimal basis B_opt_T with shape: {B_opt_T.shape}")
    print(f"Coefficients U_k for new basis (in terms of Phi_T) shape: {U_k_coeffs.shape}")

    if do_plots:
        plot_singular_values(s_vals_svd, num_optimal_K, filename=f"{plot_filename_prefix}singular_values.png" if plot_filename_prefix else None)
        plot_1d_profiles(zs, B_opt_T, f"{num_optimal_K} Optimal Basis Functions (B_opt_T)", max_plot=num_optimal_K, filename=f"{plot_filename_prefix}optimal_basis.png" if plot_filename_prefix else None)

    # Print analytical form if applicable
    # This is most meaningful if Phi_T was polynomials and z_scale_info is available
    # and orthogonalize_library was False (or if labels correctly reflect the orthogonalized basis)
    if not orthogonalize_library and _z_scale_info and library_basis_config and library_basis_config.get('type') == 'polynomial':
         print_analytical_form_polynomial(U_k_coeffs, _z_scale_info, _basis_labels, K_to_print=num_optimal_K)
    elif _basis_labels: # Print general form using provided/generated labels
         print_analytical_form_polynomial(U_k_coeffs, _z_scale_info, _basis_labels, K_to_print=num_optimal_K)


    print("--- Pipeline Completed ---")
    return B_opt_T, U_k_coeffs, s_vals_svd, _Phi_T, _z_scale_info, _basis_labels


# --- Example Usage ---
if __name__ == "__main__":
    print("--- Example: Finding Optimal Basis for Morse-like Potentials ---")

    # 0. Define common z-coordinates
    zs_example = np.linspace(1.5, 8.0, 100) # Start from z_min=1.5 Angstrom

    # --- Option A: Generate sample functions ---
    if FAF_AVAILABLE:
        print("\nGenerating sample functions using FoldedAtomicFunctions...")
        # Define parameters for generating a few systems
        # For simplicity, let's vary Morse 'a' parameter for one atom type
        # and lattice constant L_x.
        
        # Define system generation parameters
        example_sys_gen_config = {
            'L_x': (9.0, 11.0),
            # To make generated potentials look like simple 1D Morse:
            # Use fixed L_x, one atom type, specific R and E, and slice at atom's x.
            # 'L_x': (10.0, 10.0), # Fixed L_x for predictable atom x-position
            'num_atom_types': 1, 
            'atom_params': [
                # Parameters to match the "simple analytical" example for comparison:
                # R (r0) = 2.5, E (D) = 0.1. Vary 'a'.
                {'R': (2.5, 2.5), 'E': (0.1, 0.1), 'a': (1.2, 2.0), 'q': (0.0,0.0)},
            ],
            'atoms_per_cell_pattern': [[0]] # Always a single atom in the cell (at L_x/2)
        }
        example_num_systems = 5 # Will create 5 system definitions
        
        # If L_x is e.g. 10.0, atom is at x=5.0. Slice there for dx=0.
        # Pass fractional 0.5; generate_1d_potential_profiles will convert to L_x/2.
        example_x_slices_fractional = [0.5] 
        # To get a "pure" Morse from a single atom, set num_periodic_images to 0.
        example_num_periodic_images = 0 

        # Define z-weights (optional, to damp repulsion)
        # For R=2.5, minimum is at z=2.5.
        example_z0_weights = 2.5 
        example_z_decay_width = 1.0 # weight becomes 0 at z0 - decay_width = 1.5 A

        # Define library basis (polynomials)
        example_lib_basis_config = {'type': 'polynomial', 'degree': 6, 'scale_z': True}
        example_K_optimal = 3

        B_opt_T_ex, U_k_ex, s_vals_ex, _, _, _ = run_optimization_pipeline(
            common_z_coords=zs_example,
            # Sample generation
            num_systems_to_generate=example_num_systems,
            system_generation_config=example_sys_gen_config,
            x_slice_coords_for_generation=example_x_slices_fractional,
            potential_type_for_generation='morse',
            num_periodic_images_for_generation=example_num_periodic_images,
            # Weights
            z_weight_z0=example_z0_weights, 
            z_weight_decay_width=example_z_decay_width,
            # Library basis
            library_basis_config=example_lib_basis_config,
            # Orthogonalization (False by default)
            # SVD
            num_optimal_K=example_K_optimal,
            # Plotting
            do_plots=True,
            plot_filename_prefix="example_faf_"
        )
        print(f"\nExample with FAF completed. Optimal basis shape: {B_opt_T_ex.shape}")

    # --- Option B: Provide simple analytical sample functions (if FAF not available or for quick test) ---
    else: # FAF_AVAILABLE is False
        print("\nFAF not available. Generating simple analytical sample functions for demonstration...")
        num_simple_samples = 10
        Y_T_simple = np.zeros((num_simple_samples, len(zs_example)))
        a_vals = np.linspace(0.8, 2.0, num_simple_samples)
        D_val, r0_val = 0.1, 2.5
        for i in range(num_simple_samples):
            exp_term = np.exp(-a_vals[i] * (zs_example - r0_val))
            Y_T_simple[i,:] = D_val * (exp_term**2 - 2 * exp_term)
        
        example_lib_basis_config_simple = {'type': 'polynomial', 'degree': 6, 'scale_z': True}
        example_K_optimal_simple = 3
        
        # Example z-weights
        example_z0_weights_simple = 2.5
        example_z_decay_width_simple = 1.0

        B_opt_T_simple, U_k_simple, s_vals_simple, _, _, _ = run_optimization_pipeline(
            common_z_coords=zs_example,
            Y_T=Y_T_simple,
            z_weight_z0=example_z0_weights_simple,
            z_weight_decay_width=example_z_decay_width_simple,
            library_basis_config=example_lib_basis_config_simple,
            num_optimal_K=example_K_optimal_simple,
            do_plots=True,
            plot_filename_prefix="example_simple_"
        )
        print(f"\nSimple example completed. Optimal basis shape: {B_opt_T_simple.shape}")

    # --- Example of providing everything manually ---
    print("\n--- Example: Manual specification of Y_T and Phi_T ---")
    # Create dummy Y_T (2 samples, 100 z-points)
    manual_zs = np.linspace(0, 5, 50)
    manual_Y_T = np.array([
        np.sin(manual_zs),
        np.cos(manual_zs) * np.exp(-0.5 * manual_zs)
    ])
    # Create dummy Phi_T (3 basis functions: 1, z, z^2)
    manual_Phi_T = np.array([
        np.ones_like(manual_zs),
        manual_zs,
        manual_zs**2
    ])
    manual_labels = ["1", "z", "z^2"]
    manual_K = 2

    try:
        B_opt_T_man, U_k_man, s_vals_man, _, _, _ = run_optimization_pipeline(
            common_z_coords=manual_zs,
            Y_T=manual_Y_T,
            Phi_T_library=manual_Phi_T,
            library_basis_labels=manual_labels, # Important for print_analytical_form
            num_optimal_K=manual_K,
            do_plots=True,
            plot_filename_prefix="example_manual_"
        )
        print(f"\nManual example completed. Optimal basis shape: {B_opt_T_man.shape}")
        # print_analytical_form_polynomial(U_k_man, None, manual_labels, K_to_print=manual_K) # Already called inside if conditions met
    except Exception as e:
        print(f"Error in manual example: {e}")
