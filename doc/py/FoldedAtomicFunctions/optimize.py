"""optimize.py
Fit sample functions to a library basis and find an optimal reduced basis via SVD.
"""
from __future__ import annotations

import numpy as np
from numpy.linalg import lstsq # Changed from scipy.linalg
from scipy.linalg import svd   # Keep svd from scipy if preferred, or change to numpy.linalg.svd
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
    weights: np.ndarray | None = None, # Can be (Nz,) or (M, Nz)
) -> np.ndarray:
    """Return coefficient matrix S (P, M) solving Y ≈ S^T Φ.

    If *weights* are given, they are applied to both Y and Φ.
    - If weights are (Nz,), they are applied uniformly to Y and Φ.
    - If weights are (M, Nz), they are applied per sample to Y and Φ.
    """
    M, Nz = Y.shape
    P = Phi.shape[0]
    S_out = np.zeros((P, M))

    if weights is not None:
        if weights.ndim == 1: # Uniform weights for all samples
            if weights.shape[0] != Nz:
                raise ValueError("1D weights must have length Nz.")
            W_sqrt = np.sqrt(weights)
            Y_w = Y * W_sqrt[np.newaxis, :] # Apply to all samples
            Phi_w = Phi * W_sqrt[np.newaxis, :] # Apply to basis
            # Solve Y_w.T approx Phi_w.T @ S using lstsq
            # (Phi_w.T is (Nz,P), Y_w.T is (Nz,M))
            S_out, *_ = lstsq(Phi_w.T, Y_w.T, rcond=None) # Added rcond=None
        elif weights.ndim == 2: # Per-sample weights
            if weights.shape != Y.shape:
                raise ValueError("2D weights must have the same shape as Y (M, Nz).")
            for m in range(M):
                W_sqrt_m = np.sqrt(weights[m, :])
                Y_m_w = Y[m, :] * W_sqrt_m
                Phi_m_w = Phi * W_sqrt_m[np.newaxis, :] # Weight basis for this specific sample's fit
                S_out[:, m], *_ = lstsq(Phi_m_w.T, Y_m_w.T, rcond=None)
        else:
            raise ValueError("Weights must be 1D (Nz,) or 2D (M,Nz).")
    else: # Unweighted
        S_out, *_ = lstsq(Phi.T, Y.T, rcond=None)
    return S_out  # (P, M)

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
    from plot_utils import plot1D, plot_SV
    import matplotlib.pyplot as plt
    plot1D(z, B_opt, "Optimal basis")
    plot_SV(s, 2)
    plt.show()
