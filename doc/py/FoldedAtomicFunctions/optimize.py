"""optimize.py
Fit sample functions to a library basis and find an optimal reduced basis via SVD.
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import lstsq, svd
from typing import Tuple, List

__all__ = [
    "fit_coefficients",
    "optimal_basis",
]

################################################################################
# Fitting coefficients
################################################################################

def fit_coefficients(
    Y: np.ndarray,  # (M, Nz)
    Phi: np.ndarray,  # (P, Nz)
    *,
    weights: np.ndarray | None = None,
) -> np.ndarray:
    """Return coefficient matrix S (P, M) solving Y ≈ S^T Φ.

    If *weights* are given, they are applied to both Y and Φ.
    """
    if weights is not None:
        W = np.sqrt(weights)
        Y = Y * W
        Phi = Phi * W

    # least-squares fit (over z-points)
    S, *_ = lstsq(Phi.T, Y.T)
    return S  # (P, M)

################################################################################
# Optimal basis via SVD
################################################################################

def optimal_basis(
    S: np.ndarray,  # (P, M)
    Phi: np.ndarray,  # (P, Nz)
    K: int,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (B_opt, U_k, singular_values).

    B_opt has shape (K, Nz).
    """
    U, s_vals, _ = svd(S, full_matrices=False)
    U_k = U[:, :K]  # (P, K)
    B_opt = U_k.T @ Phi  # (K, Nz)
    return B_opt, U_k, s_vals

################################################################################
# Quick self-test
################################################################################

if __name__ == "__main__":
    from basis import poly_basis
    z = np.linspace(0, 5, 50)
    Phi, _, _ = poly_basis(z, 5)
    # create toy samples
    rng = np.random.default_rng(0)
    coeffs = rng.normal(size=(Phi.shape[0], 3))
    Y = (coeffs.T @ Phi) + 0.05 * rng.normal(size=(3, z.size))

    S = fit_coefficients(Y, Phi)
    B_opt, U_k, s = optimal_basis(S, Phi, K=2)
    from utils import plot_1d_profiles, plot_singular_values

    plot_1d_profiles(z, B_opt, "Optimal basis")
    plot_singular_values(s, 2)
