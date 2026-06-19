"""Atomic form factors using Cromer-Mann parametrisation.

f(Q) = sum_i a_i * exp(-b_i * (Q / (4*pi))**2) + c

Coefficients from International Tables for Crystallography (Waasmaier & Kirfel 1995).
"""
import numpy as np

# Cromer-Mann coefficients: a1,b1,a2,b2,a3,b3,a4,b4,c
# Q convention: f(Q) = Σ a_i exp(-b_i * (Q/(4π))²) + c
_CROMER_MANN = {
    1:  np.array([0.489918, 20.6593, 0.262003, 7.74039, 0.196767, 49.5519, 0.049879, 2.20159, 0.001305]),  # H
    6:  np.array([2.31, 20.8439, 1.02, 10.2075, 1.5886, 0.5687, 0.865, 51.6512, 0.2156]),  # C
    14: np.array([6.2925, 2.4386, 3.0353, 32.3337, 1.9891, 0.6785, 1.541, 81.6937, 1.1407]),  # Si
}


def cromer_mann(Z: int, Q: np.ndarray) -> np.ndarray:
    """Return f(Q) for atomic number Z at scattering vectors Q (Å⁻¹)."""
    if Z not in _CROMER_MANN:
        raise KeyError(f"Cromer-Mann coefficients not available for Z={Z}")
    c = _CROMER_MANN[Z]
    s2 = (Q / (4.0 * np.pi)) ** 2
    f = c[8]  # c term
    for i in range(4):
        f += c[i * 2] * np.exp(-c[i * 2 + 1] * s2)
    return f.astype(np.float32)


def get_form_factor_table(unique_Z: np.ndarray, Q: np.ndarray) -> np.ndarray:
    """Build table ff_table[n_species, n_Q] of precomputed form factors.

    Args:
        unique_Z: sorted array of atomic numbers present in the structure.
        Q: (n_Q,) scattering vector magnitudes in Å⁻¹.

    Returns:
        ff_table: (n_species, n_Q) float32.
    """
    n_species = len(unique_Z)
    n_Q = len(Q)
    ff = np.empty((n_species, n_Q), dtype=np.float32)
    for i, Z in enumerate(unique_Z):
        ff[i, :] = cromer_mann(int(Z), Q)
    return ff
