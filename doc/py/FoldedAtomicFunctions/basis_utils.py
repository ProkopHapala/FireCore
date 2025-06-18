from __future__ import annotations
from scipy.linalg import qr
import warnings
# --- 5. Plotting & Analysis Utilities ---

import numpy as np # type: ignore
#import matplotlib.pyplot as plt
from numpy import newaxis
from optimize import fit_coefficients # Assuming optimize.py is in the same path or accessible


xyz_imgs = lambda xs, n, Lx: np.concatenate([xs + i * Lx for i in range(-n, n + 1)])


def morse_potential(z, D, a, r0):
    """Morse potential function V(z) = D[(1 − e^{−a(z−r0)})^2 − 1]."""
    e = np.exp(-a * (z - r0))
    return D * (e**2 - 2 * e)

def gen_morse_prms(n_s, a_rng = (0.8, 2.5), r0_rng = (2.5, 4.0), D_val=0.1):
    """Generates a list of Morse potential parameters."""
    prms_list = []
    a_vals = np.random.uniform(a_rng[0], a_rng[1], n_s)
    r0_vals = np.random.uniform(r0_rng[0], r0_rng[1], n_s)
    for i in range(n_s):
        prms_list.append({'D': D_val, 'a': a_vals[i], 'r0': r0_vals[i]})
    return prms_list

def gen_morse_curves(zs, prms=None, n_s_def=10, a_rng_def=(0.8, 2.5), r0_rng_def=(2.5, 3.5), D_def=0.1):
    """Generates Morse potential curves from parameters."""
    if prms is None: prms = gen_morse_prms(n_s_def, a_rng_def, r0_rng_def, D_def)
    ys_list = [morse_potential(zs, p['D'], p['a'], p['r0']) for p in prms]
    return ys_list, prms

def coulomb2D(X: np.ndarray, Z: np.ndarray, atoms: list[dict], n: int) -> np.ndarray:
    """Periodic 1/r summed over ±n images in x (vectorised)."""
    V = np.zeros_like(X)
    Lx = X.max()
    for at in atoms:
        if abs(at.get("q", 0.0)) < 1e-12:  continue
        xs_img = xyz_imgs(np.array([at["x"]]), n, Lx)
        R = np.hypot(X[..., newaxis] - xs_img, Z[..., newaxis] - at["z"])
        V += (at["q"] / np.where(R < 1e-8, 1e-8, R)).sum(-1)
    return V

def morse2D(X: np.ndarray, Z: np.ndarray, atoms: list[dict], n: int) -> np.ndarray:
    """Periodic Morse sum over images."""
    V = np.zeros_like(X)
    Lx = X.max()
    for at in atoms:
        D, a, r0 = at["D"], at["a"], at["r0"]
        xs_img = xyz_imgs(np.array([at["x"]]), n, Lx)
        R = np.hypot(X[..., newaxis] - xs_img, Z[..., newaxis] - at["z"])
        e = np.exp(-a * (R - r0))
        V += (D * (e ** 2 - 2 * e)).sum(-1)
    return V

def cos_exp_basis(X: np.ndarray, Z: np.ndarray, nx: int, nz: int, a0: float = 0.4) -> np.ndarray:
    """Return Φ with rows cos(2πk x/Lx)·exp(-a0 j z)."""
    Lx = X.max()
    rows = [
        np.cos(2 * np.pi * k * X / Lx) * np.exp(-a0 * j * Z)
        for k in range(nx + 1)
        for j in range(1, nz + 1)
    ]
    return np.vstack([r.ravel() for r in rows])  # (P, N)

################################################################################
# Polynomial basis
################################################################################

def poly_basis( z: np.ndarray, degree: int, *, scale: bool = True,) -> Tuple[np.ndarray, Tuple[float, float], List[str]]:
    """Return matrix ``Phi`` with rows z^n (n=0..degree) and optional scaling.

    Returns (Phi, (z_min, z_range), labels)
    """
    if scale:
        z_min, z_range = float(z.min()), float(z.max() - z.min())
        if z_range == 0.0:
            # fall back to unscaled
            z_scaled = z.copy()
            z_min, z_range = 0.0, 1.0
        else:
            z_scaled = (z - z_min) / z_range
    else:
        z_scaled = z
        z_min, z_range = 0.0, 0.0  # signifies *no* scaling

    rows = [np.ones_like(z_scaled)]
    labels = ["1"]
    for n in range(1, degree + 1):
        rows.append(z_scaled ** n)
        term = f"z^{n}_scaled" if scale else f"z^{n}"
        labels.append(term)

    Phi = np.vstack(rows)
    return Phi, (z_min, z_range), labels

################################################################################
# Cutoff Polynomial basis
################################################################################

def cutoff_poly_basis( z: np.ndarray,  z_cut: float, ns: list[int]|int,  z0: float=0.0, bLabels: bool = True) -> Tuple[np.ndarray, List[str]]:
    """
    Return matrix ``Phi`` with rows (z_cut - z)^(2*n_factor) for z < z_cut, 0 otherwise,
    for each n_factor in n.

    Parameters:
    -----------
    z : np.ndarray  The z-coordinates.
    z_cut : float   The cutoff distance.
    n : list[int]  A list of integers 'n' to generate powers (z_cut - z)^n.
    bLabels : bool  If True, return labels for each basis function. 

    Returns:
    --------
    Phi_rows : np.ndarray
        Matrix where each row is a basis function. Shape (len(valid_n), len(z)).
    labels : list[str]  Labels for each basis function.
    """
    rows   = []
    labels = [] if bLabels else None
    # Base term: (z_cut - z) for z < z_cut, 0 otherwise
    z_span = z_cut - z0
    base_term = np.maximum(0, 1-(z-z0)/z_span)
    if isinstance(ns, int): ns = range(1, ns + 1)
    for n in ns:
        rows  .append(base_term**n)
        if bLabels: labels.append(f"({z_cut:.2f}-z)^{n}")
    return np.array(np.vstack(rows)), labels

def construct_composite_cutoff_basis( z: np.ndarray,  basis_definitions: list[tuple[float, list[int]]], z0: float=0.0 ) -> Tuple[np.ndarray, List[str]]:
    """Constructs a composite basis by merging several cutoff_poly_basis sets."""
    all_phi_rows, all_labels = [], []
    for z_c, p_factors in basis_definitions:
        phi_curr, labels_curr = cutoff_poly_basis(z, z_c, p_factors, z0=z0)
        if phi_curr.size > 0:
            all_phi_rows.append(phi_curr)
            all_labels.extend(labels_curr)
    if not all_phi_rows:
        return np.array([]).reshape(0, z.shape[0] if z is not None and hasattr(z, "shape") else 0), []
    return np.vstack(all_phi_rows), all_labels

def get_basis_func_map(basis_definitions: list[tuple[float, list[int]]]) -> list[tuple[int, int]]:
    """
    Generates a map from a flat index of a composite basis function to its origin.

    Parameters:
    -----------
    basis_definitions : list[tuple[float, list[int]]]
        The definition of the composite basis, e.g., [(z_cut1, [pow1, pow2]), (z_cut2, [pow3])].

    Returns:
    --------
    func_map : list[tuple[int, int]]
        A list where func_map[k] = (izcut, ipow_idx_in_zcut_list).
        `izcut` is the index of the (z_cut, [powers]) tuple in basis_definitions.
        `ipow_idx_in_zcut_list` is the index of the specific power within that z_cut's power list.
    """
    func_map = []
    for izcut, (_, pows) in enumerate(basis_definitions):
        for ipow_idx_in_zcut_list in range(len(pows)):
            func_map.append((izcut, ipow_idx_in_zcut_list))
    return func_map

def remove_basis_func_by_flat_index(
    basis_definitions: list[tuple[float, list[int]]],
    flat_idx_to_remove: int,
    func_map: list[tuple[int, int]]
) -> list[tuple[float, list[int]]]:
    """
    Removes a basis function from a copy of basis_definitions based on its flat index.
    The original basis_definitions is not modified.
    """
    if flat_idx_to_remove < 0 or flat_idx_to_remove >= len(func_map):
        raise IndexError("flat_idx_to_remove is out of bounds for the provided func_map.")

    new_bd = [list(item) for item in basis_definitions] # Make mutable copy
    new_bd = [(item[0], list(item[1])) for item in new_bd] # Make power lists mutable

    izcut_to_modify, ipow_idx_to_remove = func_map[flat_idx_to_remove]
    new_bd[izcut_to_modify][1].pop(ipow_idx_to_remove)
    # The cleanup_bd function (called by the mutation handler) will remove empty z_cut entries.
    return new_bd
################################################################################
# Orthogonalisation
################################################################################

def gram_schmidt_weighted( Phi: np.ndarray, ws: np.ndarray | None = None) -> np.ndarray:
    """Return orthonormal basis Q with respect to the weights *ws*.

    If *ws* is ``None`` it reduces to QR with column pivoting.
    """
    if ws is None:
        # simple QR (rows = basis functions)
        Q, _ = qr(Phi.T, mode="economic")
        return Q.T

    # weighted Gram–Schmidt
    Q = []
    W = np.sqrt(ws)
    for v in Phi:
        v = v * W  # weight
        for q in Q:
            v = v - np.dot(q, v) * q
        norm = np.linalg.norm(v)
        if norm < 1e-12:
            continue
        Q.append(v / norm)
    return np.vstack(Q) / W  # bring back to unweighted representation

################################################################################
# Quick self-test
################################################################################

def print_analytical_form_polynomial(
    U_k_coeffs: np.ndarray,  # (P, K)
    z_scale_info: tuple | None,  # (z_min, z_range)
    basis_labels: list | None = None,
    K_to_print: int = -1,
):
    """Pretty-print analytical form of optimal basis if *Phi* was polynomial.

    The implementation is intentionally concise and skips negligible terms.
    """
    if U_k_coeffs.size == 0:
        return
    P, K_actual = U_k_coeffs.shape
    if K_to_print < 0 or K_to_print > K_actual:
        K_to_print = K_actual
    if z_scale_info:
        z_min, z_range = z_scale_info
        print( f"Note: z_scaled = (z − {z_min:.3f}) / {z_range:.3f} (used in original basis)" )
    else:
        print("Note: original z-coordinates used (no scaling).")
    for k in range(K_to_print):
        parts: list[str] = []
        for p in range(P):
            c = U_k_coeffs[p, k]
            if abs(c) < 1e-6:
                continue
            label = basis_labels[p] if basis_labels and p < len(basis_labels) else f"phi_{p}"
            parts.append(f"({c:+.3e} * {label})")
        print(f"B_opt_{k+1}(z) = " + " ".join(parts))

def get_base_zero_arg(k: float, n: int) -> float:
    """
    Calculates the value of z_arg where the base (1 - k*z_arg/(2n)) becomes zero.
    This is the "natural" cutoff for the (1 - k*z_arg/(2n))**(2n) form.
    """
    if k <= 0: raise ValueError("k must be positive.")
    if n <= 0: raise ValueError("n must be a positive integer.")
    return (2.0 * n) / k

def get_zcut(k: float, n: int, z0_approx_point: float) -> float:
    """Calculates z_cut for polynomial (z_cut - z)^(2n) to match exp(-kz) decay at z0_approx_point."""
    if k <= 0: raise ValueError("k must be positive.")
    if n <= 0: raise ValueError("n must be a positive integer.")
    return z0_approx_point + get_base_zero_arg(k, n)

def eval_poly_exp_native_form(z_arg: np.ndarray | float, k: float, n: int) -> tuple[np.ndarray | float, float]:
    """
    Calculates the polynomial approximation term (max(0, 1 - k_decay*z_arg/(2*n_factor)))**(2*n_factor).
    This approximates exp(-k_decay*z_arg).

    # Parameters:
    z_arg : np.ndarray    Argument of the exponential decay (e.g., 'z' or 'z-r0').
    k : float             The decay constant 'k'.
    n : int               The 'n' in the polynomial (1 - kz/(2n))**(2n).

    # Returns:
    approx_values : The polynomial approximation.
    z_cut_arg :  The value of z_arg where the base (1 - k*z_arg/(2n)) becomes zero.         z_cut_arg = (2 * n_factor) / k_decay.
    """
    if k <= 0: raise ValueError("k must be positive.")
    if n <= 0: raise ValueError("n must be a positive integer.")
    term_base = 1.0 - (k * z_arg) / (2.0 * n)
    approx_values = np.maximum(0, term_base)**(2 * n)
    z_cut_arg = get_base_zero_arg(k, n) # Natural cutoff for this form's argument
    return approx_values, z_cut_arg

def eval_poly_at_z0(zs: np.ndarray, k: float, n: int, z0_approx_point: float) -> tuple[np.ndarray | float, float]:
    """
    Evaluates polynomial C*(z_cut - z)^(2n) to approximate exp(-kz) optimally at z0_approx_point.

    Parameters:
    -----------
    zs : np.ndarray
        Array of z-coordinates for evaluation.
    k : float
        Decay constant in exp(-kz).
    n : int
        Power factor 'n' in (z_cut - z)^(2n).
    z0_approx_point : float
        The z-coordinate around which the approximation is optimized.

    Returns:
    --------
    poly_values : np.ndarray
        Values of the polynomial approximation.
    z_cut : float
        The calculated z_cut for the polynomial.
    """
    z_cut = get_zcut(k, n, z0_approx_point)
    
    # C = exp(-k*z0) / (z_cut - z0)^(2n) = exp(-k*z0) / ( (2n)/k )**(2n)
    C_coeff = np.exp(-k * z0_approx_point) / (get_base_zero_arg(k, n)**(2 * n))
    
    poly_values = C_coeff * (np.maximum(0, z_cut - zs)**(2 * n))
    return poly_values, z_cut

def get_taylor_coeffs(k: float, z0: float, order: int) -> np.ndarray:
    """
    Calculates Taylor expansion coefficients for exp(-kz) around z0 up to a given order.
    The polynomial is in terms of (z-z0).
    Coefficient c_i corresponds to the term (z-z0)^i / i!
    So, c_i = f^(i)(z0) = (-k)^i * exp(-kz0)
    The full term is (c_i / i!) * (z-z0)^i.
    This function returns coeffs_for_poly1d = [c_order/order!, c_{order-1}/(order-1)!, ..., c_1/1!, c_0/0!]
    for np.poly1d, which expects coefficients from highest power to lowest.
    """
    if order < 0: raise ValueError("Order must be non-negative.")
    exp_kz0 = np.exp(-k * z0)
    # Coefficients for P(x) = a_n x^n + ... + a_1 x + a_0
    # where x = (z-z0)
    coeffs_for_poly1d = [((-k)**i * exp_kz0) / np.math.factorial(i) for i in range(order)]
    return np.array(coeffs_for_poly1d)

def eval_taylor_exp_approx(zs: np.ndarray, k: float, z0: float, order: int) -> np.ndarray:
    """
    Evaluates Taylor expansion of exp(-kz) around z0 up to a given order.
    f(z) approx sum_{i=0 to order} [ f^(i)(z0) * (z-z0)^i / i! ]
    f^(i)(z0) = (-k)^i * exp(-kz0)
    """
    if order < 0: raise ValueError("Order must be non-negative.")
    
    # Get coefficients for polynomial in (z-z0)
    # coeffs are [c_order, c_{order-1}, ..., c_0] where c_i is for (z-z0)^i
    taylor_coeffs = get_taylor_coeffs(k, z0, order)
    print(f"Taylor O({order}) at z0={z0:.1f}: Coeffs (for (z-z0)^i terms): {taylor_coeffs}")
    poly = np.poly1d(taylor_coeffs[::-1])
    return poly(zs - z0)

def get_taylor_approx_data_for_plot(zs, k, z0, order_seq):
    """
    Helper to generate data for plotFunctionApprox for Taylor series.
    The 'z_cut' equivalent for Taylor isn't a sharp cutoff, so we pass z0
    as a reference point if the plotting function needs a coordinate.
    """
    ys_approx_taylor = []
    for order_val in order_seq:
        taylor_approx_values = eval_taylor_exp_approx(zs, k, z0, order_val)
        # Get and print coefficients for inspection
        coeffs = get_taylor_coeffs(k, z0, order_val)
        # print(f"Taylor O({order_val}) at z0={z0:.1f}:")
        # print(f"  Coeffs (for (z-z0)^i terms, highest power first): {coeffs}")
        # print(f"  Evaluated: {taylor_approx_values[:3]}...") # Print a few values
        # For plotFunctionApprox, the second element in the tuple is z_cut.
        # For Taylor, there isn't a direct equivalent. We can pass z0 or None.
        ys_approx_taylor.append((taylor_approx_values, z0, f'Taylor O({order_val})'))
    return ys_approx_taylor

################################################################################
# SVD Error Calculation
################################################################################

def calc_svd_reconstruction_errors(Y_samples: np.ndarray, weights: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Calculates SVD of Y.T and absolute SVD reconstruction error of Y_samples.

    Performs SVD on Y_samples.T to find optimal basis for rows of Y_samples.
    If weights are provided, SVD is performed on weighted Y_samples.
    Returns (eigenfunctions_U, singular_values_s, errors_vs_k).
    eigenfunctions_U columns are the principal components of Y_samples.
    errors_vs_k is sqrt(sum of squares of neglected singular values), an absolute error measure.
    """
    if Y_samples.ndim != 2 or Y_samples.shape[0] == 0 or Y_samples.shape[1] == 0:
        return np.array([]), np.array([]), np.array([])

    if weights is not None:
        if weights.ndim == 1 and weights.shape[0] == Y_samples.shape[1]: # (Nz)
            W_sqrt = np.sqrt(weights)
            Y_w = Y_samples * W_sqrt[np.newaxis, :]
        elif weights.ndim == 2 and weights.shape == Y_samples.shape: # (M, Nz)
            W_sqrt = np.sqrt(weights)
            Y_w = Y_samples * W_sqrt
        else:
            raise ValueError("Weights must be 1D (Nz) or 2D (M,Nz) matching Y_samples dimensions.")
        U_eigenfuncs, s_vals, _Vh_coeffs = np.linalg.svd(Y_w.T, full_matrices=False)
    else:
        # SVD on Y.T to get singular values for reconstructing rows of Y (the functions)
        # U_cols_are_eigenfunctions, s_vals, Vh_rows_are_coeffs_in_U
        U_eigenfuncs, s_vals, _Vh_coeffs = np.linalg.svd(Y_samples.T, full_matrices=False)

    #total_sq = np.sum(s_vals**2)
    #if total_sq < 1e-12: # Handle case where all singular values are zero
    #    return U_eigenfuncs, s_vals, np.zeros_like(s_vals)
    # Cumulative sum of squared singular values from the tail
    cumsum_sq_tail = np.cumsum(s_vals[::-1]**2)[::-1]
    # Relative error for keeping k components is sqrt(sum_sq_tail[k]) / sqrt(total_sq)
    return U_eigenfuncs, s_vals, np.sqrt(cumsum_sq_tail) # error is not relative yet, just sqrt of sum of squares of neglected s_vals

################################################################################
# Basis Set Analysis Utilities
################################################################################

def calc_fit_svd(Y: np.ndarray, phi: np.ndarray, weights: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """Calc fit coeffs S for Y ~ S.T @ phi, their SVD (U_S, s_vals_S), & initial fit error.

    Parameters:
    -----------
    Y : np.ndarray (M, Nz) Sample functions.
    phi : np.ndarray (P, Nz) Basis functions.
    weights : np.ndarray (Nz), optional. Weights for the z-points.

    Returns:
    --------
    S : Coefficient matrix (P, M).
    U_S : Left singular vectors of S (P, min(P,M)).
    s_vals_S : Singular values of S.
    initial_err : Absolute weighted error of fitting Y with full phi.
    """
    S                    = fit_coefficients(Y, phi, weights=weights)
    U_S, s_vals_S, _     = np.linalg.svd(S, full_matrices=False)
    Y_reconstructed_full = (S.T @ phi)

    if weights is not None:
        # W_sqrt is calculated based on weights' dimensions
        if weights.ndim == 1: # (Nz)
            W_sqrt = np.sqrt(weights)
            initial_err = np.linalg.norm((Y - Y_reconstructed_full) * W_sqrt[np.newaxis, :])
        elif weights.ndim == 2: # (M, Nz)
            W_sqrt = np.sqrt(weights)
            initial_err = np.linalg.norm((Y - Y_reconstructed_full) * W_sqrt)
        # else: error would have been raised by fit_coefficients
    else:
        initial_err = np.linalg.norm(Y - Y_reconstructed_full)
    return S, U_S, s_vals_S, initial_err

def eval_multi_basis_recon_err( Y: np.ndarray, phi_list: list[np.ndarray],  k_vals: list[int], weights: np.ndarray | None = None ) -> tuple[list[float], dict[int, list[float]]]:
    """Eval initial fit & k-component recon errors for a list of bases."""
    err_full = []
    err_k    = {k: [] for k in k_vals}

    for phi_current in phi_list:
        # S_loop and s_vals_loop are not directly used for error accumulation here,
        # but U_S_loop and initial_err_loop are.
        _S_loop, U_S_loop, _s_vals_loop, initial_err_loop = calc_fit_svd(Y, phi_current, weights=weights)
        err_full.append(initial_err_loop)

        for k in k_vals:
            recon_err_k = np.inf
            if k > 0 and k <= U_S_loop.shape[1]: # k must be positive and within available components
                U_k_loop        = U_S_loop[:, :k]
                B_opt_rows      = (U_k_loop.T @ phi_current)
                # The coefficients for Y in B_opt should ideally also be found via weighted lstsq
                # if B_opt_rows itself is not considered "final" but an intermediate.
                # However, B_opt_rows is derived from SVD of S, which came from weighted fit.
                # So, B_opt_rows is already "optimized" for the weighted data.
                # Using unweighted lstsq here to project Y onto B_opt_rows is a common simplification.
                coeffs_Y_in_B_opt = np.linalg.lstsq(B_opt_rows.T, Y.T, rcond=None)[0]
                Y_reconstructed_k = (coeffs_Y_in_B_opt.T @ B_opt_rows)
                if weights is not None:
                    if weights.ndim == 1: # (Nz)
                        W_sqrt = np.sqrt(weights)
                        recon_err_k = np.linalg.norm((Y - Y_reconstructed_k) * W_sqrt[np.newaxis, :])
                    elif weights.ndim == 2: # (M, Nz)
                        W_sqrt = np.sqrt(weights)
                        recon_err_k = np.linalg.norm((Y - Y_reconstructed_k) * W_sqrt)
                    # else: error would have been raised by calc_fit_svd
                else:
                    recon_err_k = np.linalg.norm(Y - Y_reconstructed_k)
            err_k[k].append(recon_err_k)

    return err_full, err_k

################################################################################
# Linear dependency / similarity utilities
################################################################################

def cosine_sim_sq_matrix(rows: np.ndarray) -> np.ndarray:
    """Return matrix of squared cosine similarities between *rows*.

    The input ``rows`` is expected to have shape (P, N) where each row is a
    basis function sampled on grid.  Using squared cosine similarity avoids an
    expensive square-root.
    """
    G = rows @ rows.T                          # Gram matrix (⟨f_i|f_j⟩)
    norms_sq = np.diag(G).astype(float)
    # Guard against divide-by-zero – zero-norm rows should not appear, but keep safe
    denom = np.outer(norms_sq, norms_sq) + 1e-30
    return (G**2) / denom

def calc_ld_metrics(rows: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (ld_sum, ld_max, csm) in one shot.
    - ld_sum: Σ_j cos²
    - ld_max: max_j cos²
    - csm   : full squared-cosine matrix (diagonal zeroed)
    """
    csm = cosine_sim_sq_matrix(rows)
    np.fill_diagonal(csm, 0.0)
    ld_sum = csm.sum(axis=1)
    ld_max = csm.max(axis=1)
    return ld_sum, ld_max, csm

# def linear_dependency_scores(rows: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
#     """
#     Calculates linear dependency scores for each row (basis function).

#     Returns:
#     --------
#     ld_sum_scores : np.ndarray
#         For each row i, LD_sum_i = Σ_{j≠i} cos²(f_i,f_j).
#         A higher value means row i is more expressible by others.
#     ld_max_pair_scores : np.ndarray
#         For each row i, LD_max_pair_i = max_{j≠i} cos²(f_i,f_j).
#         Indicates the strongest pairwise collinearity for row i.
#     """
#     ld_sum, ld_max, _ = calc_ld_metrics(rows)
#     return ld_sum, ld_max

################################################################################
# Generic utility functions (added during refactor)
################################################################################

def weighted_rmse(y_true: np.ndarray, y_pred: np.ndarray, weights: np.ndarray | None = None) -> float:
    """Compute Root Mean Square Error with optional weighting.

    Parameters
    ----------
    y_true, y_pred : np.ndarray
        Arrays of the same shape containing reference and predicted values.
    weights : np.ndarray | None
        If given, must be broadcastable to *y_true* and *y_pred*; each element
        scales the corresponding squared error. A value of zero excludes that
        point.

    Returns
    -------
    float
        The weighted RMSE (always positive, +1e-12 safety added).
    """
    diff_sq = (y_true - y_pred) ** 2
    if weights is not None:
        diff_sq = diff_sq * weights
        denom = float(np.sum(weights))
        if denom == 0.0:
            return float('inf')
    else:
        denom = diff_sq.size
    return float(np.sqrt(np.sum(diff_sq) / denom) + 1e-12)


def get_basis_metrics(basis_def: list[tuple[float, list[int]]]) -> tuple[int, int, int, int]:
    """Return nBasis, nZcut, nMax, nSum for a given *basis_def*.

    Notes
    -----
    * *basis_def* is a list of tuples ``(z_cut, [powers])``.
    """
    n_basis = sum(len(pows) for _, pows in basis_def)
    n_zcut = len(basis_def)
    all_pows = [p for _, pows in basis_def for p in pows]
    if not all_pows:
        return 0, n_zcut, 0, 0
    n_max = max(all_pows)
    n_sum = sum(all_pows)
    return n_basis, n_zcut, n_max, n_sum