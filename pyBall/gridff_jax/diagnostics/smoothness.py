"""Energy surface smoothness analysis for any MLIP or model.

Evaluates how smooth the potential energy surface is by computing
numerical derivatives at multiple z-heights over a lateral (u,v) grid
and checking for discontinuities, oscillations, and non-physical features.

Usage
-----
    from pyBall.gridff_jax.diagnostics.smoothness import smoothness_analysis
    result = smoothness_analysis(evaluate_fn, density, adsorbate, ...)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable

import numpy as np

from ..interfaces import AdsorbateDefinition, DensityData, PoseBatch
from ..pose_sampling.rigid import infer_reference_planes, transform_adsorbate
from ..utils import normalize_quaternion


@dataclass
class SmoothnessResult:
    """Results from energy surface smoothness analysis."""
    z_heights: list[float] = field(default_factory=list)
    n_grid: int = 0
    energy_maps: dict[float, np.ndarray] = field(default_factory=dict)
    grad_x_maps: dict[float, np.ndarray] = field(default_factory=dict)
    grad_y_maps: dict[float, np.ndarray] = field(default_factory=dict)
    laplacian_maps: dict[float, np.ndarray] = field(default_factory=dict)
    roughness_per_z: dict[float, float] = field(default_factory=dict)
    max_gradient_per_z: dict[float, float] = field(default_factory=dict)
    u_vals: np.ndarray | None = None
    v_vals: np.ndarray | None = None
    metadata: dict[str, Any] = field(default_factory=dict)


def _build_lateral_grid_poses(
    density: DensityData,
    adsorbate: AdsorbateDefinition,
    n_grid: int,
    z_height: float,
    quaternion: np.ndarray,
) -> tuple[PoseBatch, np.ndarray, np.ndarray]:
    """Build a uniform (u,v) grid of poses at fixed z and orientation."""
    planes = infer_reference_planes(density)
    z_ref = float(planes["z_ref"])
    quaternion = normalize_quaternion(quaternion)

    u_vals = np.linspace(0.0, 1.0, n_grid, endpoint=False)
    v_vals = np.linspace(0.0, 1.0, n_grid, endpoint=False)
    uu, vv = np.meshgrid(u_vals, v_vals, indexing="ij")
    uv_flat = np.column_stack([uu.ravel(), vv.ravel()])

    positions = []
    pose_params = []
    for uv in uv_flat:
        pos = transform_adsorbate(adsorbate, density, uv, z_height, quaternion, z_ref)
        positions.append(pos)
        pose_params.append(np.concatenate([uv, [z_height], quaternion]))

    poses = PoseBatch(
        adsorbate=adsorbate,
        positions=np.asarray(positions, dtype=float),
        pose_params=np.asarray(pose_params, dtype=float),
        site_labels=["grid"] * len(uv_flat),
        metadata={"z_ref": z_ref, "z_height": z_height, "scan_type": "smoothness"},
    )
    return poses, u_vals, v_vals


def smoothness_analysis(
    evaluate_fn: Callable[[DensityData, PoseBatch], Any],
    density: DensityData,
    adsorbate: AdsorbateDefinition,
    z_heights: list[float] | None = None,
    n_grid: int = 25,
    quaternion: np.ndarray | None = None,
    verbose: bool = True,
) -> SmoothnessResult:
    """Analyze energy surface smoothness over lateral (u,v) grid.

    Parameters
    ----------
    evaluate_fn : callable
        ``(density, poses) -> result`` where ``result.energies`` is an array.
    density : DensityData
        Substrate data.
    adsorbate : AdsorbateDefinition
        Molecule definition.
    z_heights : list[float], optional
        Heights above surface to analyze. Default: [2.0, 2.5, 3.0, 3.5, 4.0].
    n_grid : int
        Grid points per lateral dimension.
    quaternion : np.ndarray, optional
        Fixed orientation. Default: identity (flat on surface).
    verbose : bool
        Print progress.

    Returns
    -------
    SmoothnessResult
        Energy maps, gradient maps, Laplacian, and roughness metrics per z.
    """
    if z_heights is None:
        z_heights = [2.0, 2.5, 3.0, 3.5, 4.0]
    if quaternion is None:
        quaternion = np.array([1.0, 0.0, 0.0, 0.0], dtype=float)

    result = SmoothnessResult(
        z_heights=list(z_heights),
        n_grid=n_grid,
    )

    for z_h in z_heights:
        if verbose:
            print(f"[smoothness] z = {z_h:.2f} A: evaluating {n_grid}x{n_grid} grid ...",
                  flush=True)

        poses, u_vals, v_vals = _build_lateral_grid_poses(
            density, adsorbate, n_grid, z_h, quaternion,
        )
        result.u_vals = u_vals
        result.v_vals = v_vals

        eval_result = evaluate_fn(density, poses)
        energies = np.asarray(eval_result.energies, dtype=float)
        E_2d = energies.reshape(n_grid, n_grid)
        result.energy_maps[z_h] = E_2d

        du = 1.0 / n_grid
        grad_x = np.zeros_like(E_2d)
        grad_y = np.zeros_like(E_2d)
        for i in range(n_grid):
            ip = (i + 1) % n_grid
            im = (i - 1) % n_grid
            grad_x[i, :] = (E_2d[ip, :] - E_2d[im, :]) / (2.0 * du)
        for j in range(n_grid):
            jp = (j + 1) % n_grid
            jm = (j - 1) % n_grid
            grad_y[:, j] = (E_2d[:, jp] - E_2d[:, jm]) / (2.0 * du)

        result.grad_x_maps[z_h] = grad_x
        result.grad_y_maps[z_h] = grad_y

        lap = np.zeros_like(E_2d)
        for i in range(n_grid):
            ip = (i + 1) % n_grid
            im = (i - 1) % n_grid
            lap += (E_2d[ip, :] - 2.0 * E_2d + E_2d[im, :]) / (du ** 2)
        for j in range(n_grid):
            jp = (j + 1) % n_grid
            jm = (j - 1) % n_grid
            lap += (E_2d[:, jp] - 2.0 * E_2d + E_2d[:, jm]) / (du ** 2)
        result.laplacian_maps[z_h] = lap

        grad_magnitude = np.sqrt(grad_x ** 2 + grad_y ** 2)
        roughness = float(np.std(lap))
        max_grad = float(grad_magnitude.max())
        result.roughness_per_z[z_h] = roughness
        result.max_gradient_per_z[z_h] = max_grad

        if verbose:
            corrugation = float(E_2d.max() - E_2d.min()) * 1000.0
            print(f"  corrugation = {corrugation:.1f} meV, "
                  f"roughness(std(lap)) = {roughness:.4f}, "
                  f"max|grad| = {max_grad:.4f}")

    result.metadata = {
        "n_grid": n_grid,
        "z_heights": list(z_heights),
        "adsorbate": adsorbate.name,
    }
    return result
