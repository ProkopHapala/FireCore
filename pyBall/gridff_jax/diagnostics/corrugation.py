"""Lateral corrugation analysis for any MLIP or model.

Computes and compares the lateral energy corrugation (variation of
interaction energy across the surface unit cell at fixed z-heights)
for one or more evaluation methods.

Usage
-----
    from pyBall.gridff_jax.diagnostics.corrugation import corrugation_analysis
    result = corrugation_analysis(methods, density, adsorbate, ...)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable

import numpy as np

from ..interfaces import AdsorbateDefinition, DensityData, PoseBatch
from ..pose_sampling.rigid import infer_reference_planes, transform_adsorbate
from ..utils import normalize_quaternion


@dataclass
class CorrugationResult:
    """Results from lateral corrugation analysis."""
    z_heights: list[float] = field(default_factory=list)
    n_grid: int = 0
    u_vals: np.ndarray | None = None
    v_vals: np.ndarray | None = None
    energy_maps: dict[str, dict[float, np.ndarray]] = field(default_factory=dict)
    corrugation_amplitude: dict[str, dict[float, float]] = field(default_factory=dict)
    corrugation_rms: dict[str, dict[float, float]] = field(default_factory=dict)
    site_energies: dict[str, dict[float, dict[str, float]]] = field(default_factory=dict)
    metadata: dict[str, Any] = field(default_factory=dict)


def _build_grid_poses(
    density: DensityData,
    adsorbate: AdsorbateDefinition,
    n_grid: int,
    z_height: float,
    quaternion: np.ndarray,
) -> tuple[PoseBatch, np.ndarray, np.ndarray]:
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
        metadata={"z_ref": z_ref, "z_height": z_height},
    )
    return poses, u_vals, v_vals


def corrugation_analysis(
    methods: dict[str, Callable[[DensityData, PoseBatch], Any]],
    density: DensityData,
    adsorbate: AdsorbateDefinition,
    z_heights: list[float] | None = None,
    n_grid: int = 25,
    quaternion: np.ndarray | None = None,
    verbose: bool = True,
) -> CorrugationResult:
    """Compare lateral corrugation across multiple evaluation methods.

    Parameters
    ----------
    methods : dict[str, callable]
        Mapping of method name to evaluation function.
        Each function has signature ``(density, poses) -> result``
        where ``result.energies`` is an array of shape ``(n_poses,)``.
    density : DensityData
        Substrate data.
    adsorbate : AdsorbateDefinition
        Molecule definition.
    z_heights : list[float], optional
        Heights above surface. Default: [2.0, 2.5, 3.0, 4.0, 5.0].
    n_grid : int
        Grid points per lateral dimension.
    quaternion : np.ndarray, optional
        Fixed orientation. Default: identity.
    verbose : bool
        Print progress.

    Returns
    -------
    CorrugationResult
        Per-method, per-height energy maps and corrugation metrics.
    """
    if z_heights is None:
        z_heights = [2.0, 2.5, 3.0, 4.0, 5.0]
    if quaternion is None:
        quaternion = np.array([1.0, 0.0, 0.0, 0.0], dtype=float)

    result = CorrugationResult(
        z_heights=list(z_heights),
        n_grid=n_grid,
    )

    from ..pose_sampling.sites import infer_surface_sites
    sites = infer_surface_sites(density)
    site_uv_map = {s.label: np.asarray(s.uv, dtype=float) for s in sites}

    for method_name, eval_fn in methods.items():
        result.energy_maps[method_name] = {}
        result.corrugation_amplitude[method_name] = {}
        result.corrugation_rms[method_name] = {}
        result.site_energies[method_name] = {}

        for z_h in z_heights:
            if verbose:
                print(f"[corrugation] {method_name} @ z={z_h:.2f} A ...", flush=True)

            poses, u_vals, v_vals = _build_grid_poses(
                density, adsorbate, n_grid, z_h, quaternion,
            )
            result.u_vals = u_vals
            result.v_vals = v_vals

            eval_result = eval_fn(density, poses)
            energies = np.asarray(eval_result.energies, dtype=float)
            E_2d = energies.reshape(n_grid, n_grid)
            result.energy_maps[method_name][z_h] = E_2d

            amplitude = float(E_2d.max() - E_2d.min()) * 1000.0
            rms = float(np.std(E_2d)) * 1000.0
            result.corrugation_amplitude[method_name][z_h] = amplitude
            result.corrugation_rms[method_name][z_h] = rms

            site_e = {}
            for slabel, suv in site_uv_map.items():
                iu = int(np.argmin(np.abs(u_vals - suv[0])))
                iv = int(np.argmin(np.abs(v_vals - suv[1])))
                site_e[slabel] = float(E_2d[iu, iv]) * 1000.0
            result.site_energies[method_name][z_h] = site_e

            if verbose:
                print(f"  corrugation = {amplitude:.1f} meV (rms = {rms:.1f} meV)")

    result.metadata = {
        "n_grid": n_grid,
        "z_heights": list(z_heights),
        "methods": list(methods.keys()),
        "adsorbate": adsorbate.name,
        "sites": {k: v.tolist() for k, v in site_uv_map.items()},
    }
    return result
