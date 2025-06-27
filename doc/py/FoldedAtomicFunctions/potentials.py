"""potentials.py
Generation of sample potential functions in 1-D and 2-D.
All functions are kept minimalistic and functional.
"""

from __future__ import annotations

import numpy as np
from typing import Sequence, Dict, List

__all__ = [
    "morse_potential",
    "generate_morse_samples",
    "cut_2d_potential_profiles",
]

################################################################################
# 1-D Morse potential
################################################################################

def morse_potential(z: np.ndarray, D: float, a: float, r0: float) -> np.ndarray:
    """Return Morse potential V(z) = D[(1 − e^{−a(z−r0)})^2 − 1]."""
    exp_term = np.exp(-a * (z - r0))
    return D * (exp_term ** 2 - 2 * exp_term)


def generate_morse_samples(
    z: np.ndarray,
    params: Sequence[Dict[str, float]],
) -> np.ndarray:
    """Generate a matrix of Morse potentials.

    Parameters
    ----------
    z       : 1-D grid.
    params  : sequence of dicts with keys ``D``, ``a`` and ``r0``.

    Returns
    -------
    ndarray of shape (len(params), z.size)
    """
    return np.vstack([morse_potential(z, **p) for p in params])

################################################################################
# 2-D potential cut using FoldedAtomicFunctions (optional)
################################################################################

try:
    # heavy import, only if present
    from FoldedAtomicFunctions import GridManager, PotentialCalculator  # type: ignore

    _FAF_AVAILABLE = True
except ImportError:  # pragma: no cover – optional dependency
    GridManager = PotentialCalculator = None  # type: ignore
    _FAF_AVAILABLE = False


def _assert_faf():
    if not _FAF_AVAILABLE:
        raise ImportError(
            "FoldedAtomicFunctions not available – install / ensure it is on PYTHONPATH"
        )


def cut_2d_potential_profiles(
    system_definition: Dict,
    z: np.ndarray,
    x_slices: Sequence[float],
    num_periodic_images: int = 5,
) -> np.ndarray:
    """Return potential profiles V(x_i, z) for *x_slices* taken from a 2-D periodic
    system described by *system_definition*.

    The implementation is intentionally minimal.  It relies on the classes
    provided by ``FoldedAtomicFunctions.py`` and does **no** error handling.
    """
    _assert_faf()

    cell_x = float(system_definition.get("L_x", 10.0))
    cell_z = float(system_definition.get("L_z", 10.0))
    grid_step = float(system_definition.get("grid_step", 0.2))
    gm = GridManager(cell_x=cell_x, cell_z=cell_z, grid_step=grid_step)

    pc = PotentialCalculator(gm)
    for at in system_definition.get("atoms", []):
        pc.add_atom(**at)

    profiles: List[np.ndarray] = []
    for x in x_slices:
        vz = pc.potential_along_z(x, z)  # type: ignore[attr-defined]
        profiles.append(vz)

    return np.vstack(profiles)

################################################################################
# Quick self-test
################################################################################

if __name__ == "__main__":
    # Demonstration of 1-D Morse generation
    zs     = np.linspace(0.0, 10.0, 200)
    rng    = np.random.default_rng(0)
    params = [ {"D": 1.0,  "a": rng.uniform(1.5, 1.8),  "r0": rng.uniform(2.5, 4.0)} for _ in range(6) ]  # scatter in a and r0
    Y      = generate_morse_samples(zs, params)
    import matplotlib.pyplot as plt
    from plot_utils import plot1D
    plot1D(zs, Y, "Morse samples", scMin=1.2)
    plt.show()
