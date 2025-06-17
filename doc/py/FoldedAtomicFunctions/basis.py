"""basis.py
Creation and orthogonalisation of analytical basis functions.
All functions are light-weight and purely functional.
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import qr
from typing import Tuple, List, Dict

__all__ = [
    "poly_basis",
    "cutoff_poly_basis",
    "gram_schmidt_weighted",
]

################################################################################
# Polynomial basis
################################################################################

def poly_basis(
    z: np.ndarray,
    degree: int,
    *,
    scale: bool = True,
) -> Tuple[np.ndarray, Tuple[float, float], List[str]]:
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

def cutoff_poly_basis(
    z: np.ndarray,
    z_cut: float,
    max_power_factor: int, # Corresponds to 'n' in (z_cut-z)^(2n), so 2n is the max power
    *,
    scale: bool = False, # Scaling is less common for this type but can be added
) -> Tuple[np.ndarray, Tuple[float, float], List[str]]:
    """Return matrix ``Phi`` with rows (z_cut - z_eff)^(2*n) for z_eff < z_cut_eff, 0 otherwise.

    n ranges from 1 to max_power_factor.
    z_eff and z_cut_eff are potentially scaled versions of z and z_cut.

    Returns (Phi, (z_min, z_range), labels)
    """
    z_eff = z
    z_cut_eff = z_cut
    scale_info = (0.0, 0.0) # Signifies no scaling by default
    scale_suffix = ""

    if scale:
        z_min, z_range = float(z.min()), float(z.max() - z.min())
        if z_range == 0.0:
            z_min, z_range = 0.0, 1.0 # Fallback, no effective scaling
        else:
            z_eff = (z - z_min) / z_range
            z_cut_eff = (z_cut - z_min) / z_range
        scale_info = (z_min, z_range)
        scale_suffix = "_scaled"

    rows   = []
    labels = []
    
    # Base term: (z_cut_eff - z_eff) for z_eff < z_cut_eff, 0 otherwise
    base_term = np.maximum(0, z_cut_eff - z_eff)

    for n_factor in range(1, max_power_factor + 1):
        rows.append(base_term**(2 * n_factor))
        labels.append(f"(z_cut{scale_suffix}-z{scale_suffix})^{2*n_factor}")

    return np.vstack(rows) if rows else np.array([]).reshape(0, len(z)), scale_info, labels

################################################################################
# Orthogonalisation
################################################################################

def gram_schmidt_weighted(
    Phi: np.ndarray, ws: np.ndarray | None = None
) -> np.ndarray:
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

if __name__ == "__main__":
    z = np.linspace(0, 5, 100)

    print("--- Polynomial Basis Test ---")
    Phi, scale_info, labels = poly_basis(z, 4)
    Q = gram_schmidt_weighted(Phi)
    print("Phi shape", Phi.shape, "Q shape", Q.shape)
    print("Labels:", labels)
    print("Scale Info:", scale_info)

    from utils import plot_1d_profiles
    plot_1d_profiles(z, Phi, "Raw polynomial basis")
    plot_1d_profiles(z, Q, "Orthogonalised basis")

    print("\n--- Cutoff Polynomial Basis Test ---")
    z_cut_test = 3.0
    max_n_test = 3
    Phi_cutoff, scale_info_cutoff, labels_cutoff = cutoff_poly_basis(z, z_cut_test, max_n_test)
    Q_cutoff = gram_schmidt_weighted(Phi_cutoff)
    print("Phi_cutoff shape", Phi_cutoff.shape)
    print("Scale Info:", scale_info_cutoff)
    print("Labels:", labels_cutoff)
    plot_1d_profiles(z, Phi_cutoff, f"Cutoff polynomial basis (z_cut={z_cut_test}, max_2n={2*max_n_test})")
    plot_1d_profiles(z, Q_cutoff, f"Orthogonalised cutoff polynomial basis (z_cut={z_cut_test}, max_2n={2*max_n_test})")
