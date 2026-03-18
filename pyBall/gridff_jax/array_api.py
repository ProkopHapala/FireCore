"""Optional JAX backend helpers with NumPy/SciPy fallback."""

from __future__ import annotations

from dataclasses import dataclass
from importlib import import_module
from importlib.util import find_spec
from typing import Callable

import numpy as np
from scipy import ndimage as sp_ndimage

HAS_JAX = find_spec("jax") is not None
if HAS_JAX:
    jax = import_module("jax")
    jax.config.update("jax_enable_x64", True)
    jnp = import_module("jax.numpy")
    jax_map_coordinates = import_module("jax.scipy.ndimage").map_coordinates
else:
    jax = None
    jnp = None
    jax_map_coordinates = None

HAS_OPTAX = find_spec("optax") is not None
if HAS_OPTAX:
    optax = import_module("optax")
else:
    optax = None


@dataclass(frozen=True)
class ArrayBackend:
    name: str
    xp: object
    map_coordinates: Callable
    uses_jax: bool


NUMPY_BACKEND = ArrayBackend(
    name="numpy",
    xp=np,
    map_coordinates=sp_ndimage.map_coordinates,
    uses_jax=False,
)

JAX_BACKEND = (
    ArrayBackend(
        name="jax",
        xp=jnp,
        map_coordinates=jax_map_coordinates,
        uses_jax=True,
    )
    if HAS_JAX
    else None
)


def get_backend(prefer_jax: bool = True) -> ArrayBackend:
    if prefer_jax and HAS_JAX:
        return JAX_BACKEND
    return NUMPY_BACKEND


def backend_summary() -> dict:
    devices = []
    backend = get_backend().name
    default_backend = backend
    if HAS_JAX:
        devices = [f"{dev.platform}:{dev.id}" for dev in jax.devices()]
        default_backend = jax.default_backend()
    return {
        "jax": HAS_JAX,
        "optax": HAS_OPTAX,
        "default_backend": default_backend,
        "devices": devices,
    }


def as_numpy(array):
    if array is None:
        return None
    if HAS_JAX and hasattr(array, "__array__") is False and hasattr(array, "shape"):
        return np.asarray(array)
    return np.asarray(array)


def finite_difference_force(energy_fn, positions, eps: float = 1.0e-4):
    positions = np.asarray(positions, dtype=float)
    forces = np.zeros_like(positions)
    for ia in range(positions.shape[0]):
        for ic in range(3):
            dp = np.zeros_like(positions)
            dp[ia, ic] = eps
            ep = energy_fn(positions + dp)
            em = energy_fn(positions - dp)
            forces[ia, ic] = -(ep - em) / (2.0 * eps)
    return forces
