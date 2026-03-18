"""Geometry and serialization helpers used across the GridFF JAX workflow."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

import numpy as np


def ensure_dir(path: str | Path) -> Path:
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def save_json(payload, path: str | Path):
    path = Path(path)
    ensure_dir(path.parent)
    with path.open("w", encoding="utf8") as fout:
        json.dump(payload, fout, indent=2)


def normalize_quaternion(quaternion):
    quaternion = np.asarray(quaternion, dtype=float)
    norm = np.linalg.norm(quaternion)
    if norm < 1.0e-12:
        return np.array([1.0, 0.0, 0.0, 0.0], dtype=float)
    return quaternion / norm


def quaternion_to_matrix(quaternion):
    w, x, y, z = normalize_quaternion(quaternion)
    return np.array(
        [
            [1.0 - 2.0 * (y * y + z * z), 2.0 * (x * y - z * w), 2.0 * (x * z + y * w)],
            [2.0 * (x * y + z * w), 1.0 - 2.0 * (x * x + z * z), 2.0 * (y * z - x * w)],
            [2.0 * (x * z - y * w), 2.0 * (y * z + x * w), 1.0 - 2.0 * (x * x + y * y)],
        ],
        dtype=float,
    )


def random_quaternions(n: int, rng: np.random.Generator):
    u1 = rng.random(n)
    u2 = rng.random(n)
    u3 = rng.random(n)
    return np.stack(
        [
            np.sqrt(1.0 - u1) * np.sin(2.0 * np.pi * u2),
            np.sqrt(1.0 - u1) * np.cos(2.0 * np.pi * u2),
            np.sqrt(u1) * np.sin(2.0 * np.pi * u3),
            np.sqrt(u1) * np.cos(2.0 * np.pi * u3),
        ],
        axis=1,
    )[:, [3, 0, 1, 2]]


def fractional_to_cartesian(frac, cell):
    frac = np.asarray(frac, dtype=float)
    cell = np.asarray(cell, dtype=float)
    return frac @ cell


def cartesian_to_fractional(cart, cell):
    cart = np.asarray(cart, dtype=float)
    cell = np.asarray(cell, dtype=float)
    return np.linalg.solve(cell.T, cart.T).T


def wrap_fractional_uv(frac):
    frac = np.asarray(frac, dtype=float).copy()
    frac[..., 0:2] = frac[..., 0:2] - np.floor(frac[..., 0:2])
    return frac


def voxel_spacing_from_cell(cell, grid_shape_xyz: Iterable[int]):
    cell = np.asarray(cell, dtype=float)
    nx, ny, nz = [int(v) for v in grid_shape_xyz]
    return np.array(
        [
            np.linalg.norm(cell[0]) / max(nx, 1),
            np.linalg.norm(cell[1]) / max(ny, 1),
            np.linalg.norm(cell[2]) / max(nz, 1),
        ],
        dtype=float,
    )


def cluster_z_layers(zs, tol: float = 0.35):
    order = np.argsort(zs)
    clusters = []
    current = [int(order[0])]
    for idx in order[1:]:
        if abs(zs[idx] - zs[current[-1]]) <= tol:
            current.append(int(idx))
        else:
            clusters.append(current)
            current = [int(idx)]
    clusters.append(current)
    return clusters


def unique_rows_periodic(values, decimals: int = 5):
    values = np.asarray(values, dtype=float)
    if values.size == 0:
        return values.reshape(0, values.shape[-1] if values.ndim else 0)
    wrapped = values - np.floor(values)
    key = np.round(wrapped, decimals=decimals)
    _, indices = np.unique(key, axis=0, return_index=True)
    indices.sort()
    return wrapped[indices]
