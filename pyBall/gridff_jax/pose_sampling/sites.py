"""Surface-site inference for rigid adsorbate placement."""

from __future__ import annotations

import itertools

import numpy as np

from ..interfaces import DensityData, SurfaceSite
from ..utils import cartesian_to_fractional, cluster_z_layers, unique_rows_periodic


def _layer_indices(density: DensityData, tol: float = 0.35):
    zs = density.positions[:, 2]
    clusters = cluster_z_layers(zs, tol=tol)
    clusters.sort(key=lambda items: np.mean(zs[items]))
    return clusters


def _minimum_image_uv(delta_uv):
    delta = np.asarray(delta_uv, dtype=float)
    return delta - np.round(delta)


def infer_surface_sites(density: DensityData, layer_tol: float = 0.35, neighbor_scale: float = 1.18):
    layers = _layer_indices(density, tol=layer_tol)
    top_layer = layers[-1]
    top_positions = density.positions[top_layer]
    top_uv = cartesian_to_fractional(top_positions, density.cell)[:, :2]
    top_uv = top_uv - np.floor(top_uv)
    sites = [SurfaceSite(label="top", uv=uv) for uv in unique_rows_periodic(top_uv)]

    if len(top_uv) < 2:
        return sites

    cell2 = density.cell[:2, :2]
    pairs = []
    min_dist = None
    for i, j in itertools.combinations(range(len(top_uv)), 2):
        duv = _minimum_image_uv(top_uv[j] - top_uv[i])
        dxy = duv @ cell2
        dist = np.linalg.norm(dxy)
        if dist < 1.0e-8:
            continue
        if min_dist is None or dist < min_dist:
            min_dist = dist
        pairs.append((i, j, duv, dist))

    if min_dist is None:
        return sites

    bridges = []
    neighbor_pairs = []
    for i, j, duv, dist in pairs:
        if dist <= neighbor_scale * min_dist:
            neighbor_pairs.append((i, j))
            bridges.append((top_uv[i] + 0.5 * duv) % 1.0)
    for uv in unique_rows_periodic(np.asarray(bridges)):
        sites.append(SurfaceSite(label="bridge", uv=uv))

    hollows = []
    neighbor_set = {tuple(sorted(pair)) for pair in neighbor_pairs}
    for i, j, k in itertools.combinations(range(len(top_uv)), 3):
        if (
            tuple(sorted((i, j))) not in neighbor_set
            or tuple(sorted((j, k))) not in neighbor_set
            or tuple(sorted((i, k))) not in neighbor_set
        ):
            continue
        duv_ij = _minimum_image_uv(top_uv[j] - top_uv[i])
        duv_ik = _minimum_image_uv(top_uv[k] - top_uv[i])
        area = abs(np.cross(duv_ij, duv_ik))
        if area < 1.0e-6:
            continue
        hollows.append((top_uv[i] + (duv_ij + duv_ik) / 3.0) % 1.0)
    for uv in unique_rows_periodic(np.asarray(hollows)):
        sites.append(SurfaceSite(label="hollow", uv=uv))
    return sites
