"""Rigid pose generation and Cartesian placement."""

from __future__ import annotations

import math

import numpy as np

from ..interfaces import AdsorbateDefinition, DensityData, PoseBatch
from ..pbc import wrap_anchor_into_cell
from ..utils import (
    fractional_to_cartesian,
    normalize_quaternion,
    quaternion_to_matrix,
    random_quaternions,
)
from .sites import infer_surface_sites


def _axis_angle_quaternion(axis, angle):
    axis = np.asarray(axis, dtype=float)
    norm = np.linalg.norm(axis)
    if norm < 1.0e-12:
        return np.array([1.0, 0.0, 0.0, 0.0], dtype=float)
    axis = axis / norm
    half = 0.5 * float(angle)
    return np.array(
        [
            math.cos(half),
            axis[0] * math.sin(half),
            axis[1] * math.sin(half),
            axis[2] * math.sin(half),
        ],
        dtype=float,
    )


def _biased_linspace(n: int, lo: float, hi: float, power: float):
    if n <= 1:
        return np.array([lo], dtype=float)
    grid = np.linspace(0.0, 1.0, int(n), dtype=float)
    if abs(power - 1.0) > 1.0e-12:
        grid = np.power(grid, float(power))
    return lo + (hi - lo) * grid


def _sample_biased_uniform(rng: np.random.Generator, lo: float, hi: float, power: float):
    value = rng.random()
    if abs(power - 1.0) > 1.0e-12:
        value = value**float(power)
    return lo + (hi - lo) * value


def _systematic_quaternions(adsorbate: AdsorbateDefinition, count: int, tilt_deg: float):
    count = max(1, int(count))
    if adsorbate.natoms == 1:
        return np.array([[1.0, 0.0, 0.0, 0.0]], dtype=float)
    tilt = math.radians(min(max(float(tilt_deg), 0.0), 89.0))
    half_tilt = 0.5 * tilt if tilt > 1.0e-8 else 0.0
    base = [
        np.array([1.0, 0.0, 0.0, 0.0], dtype=float),
        _axis_angle_quaternion([0.0, 0.0, 1.0], 0.5 * math.pi),
        _axis_angle_quaternion([0.0, 0.0, 1.0], math.pi),
        _axis_angle_quaternion([0.0, 0.0, 1.0], 1.5 * math.pi),
        _axis_angle_quaternion([1.0, 0.0, 0.0], half_tilt),
        _axis_angle_quaternion([1.0, 0.0, 0.0], -half_tilt),
        _axis_angle_quaternion([0.0, 1.0, 0.0], half_tilt),
        _axis_angle_quaternion([0.0, 1.0, 0.0], -half_tilt),
    ]
    quats = np.vstack(base)
    if count <= len(quats):
        return quats[:count]
    repeats = int(np.ceil(count / len(quats)))
    tiled = np.tile(quats, (repeats, 1))
    return tiled[:count]


def infer_reference_planes(density: DensityData, layer_tol: float = 0.35):
    zs = np.asarray(density.positions[:, 2], dtype=float)
    z_ref = zs.max()
    counts = []
    for value in np.unique(np.round(zs, decimals=4)):
        count = np.count_nonzero(np.abs(zs - value) <= layer_tol)
        counts.append((count, value))
    if counts:
        counts.sort(key=lambda item: (item[1], item[0]))
        z_ref = counts[-1][1]
    z_image = z_ref + density.metadata.get("image_plane_offset", 1.0)
    return {"z_ref": float(z_ref), "z_image": float(z_image)}


def transform_adsorbate(
    adsorbate: AdsorbateDefinition,
    density: DensityData,
    uv: np.ndarray,
    z_height: float,
    quaternion: np.ndarray,
    z_ref: float,
):
    quaternion = normalize_quaternion(quaternion)
    rot = quaternion_to_matrix(quaternion)
    local = adsorbate.positions - adsorbate.positions[adsorbate.anchor_index]
    rotated = local @ rot.T
    # Wrap (u,v) into [0, 1) before converting to Cartesian
    uv_w = np.asarray(uv, dtype=float).copy()
    uv_w[:2] = uv_w[:2] - np.floor(uv_w[:2])
    anchor_frac = np.array([uv_w[0], uv_w[1], 0.0], dtype=float)
    anchor_cart = fractional_to_cartesian(anchor_frac, density.cell)
    anchor_cart[2] = z_ref + z_height
    positions = rotated + anchor_cart[None, :]
    # PBC: translate entire molecule so anchor is inside cell
    positions = wrap_anchor_into_cell(
        positions, adsorbate.anchor_index, density.cell, density.origin,
    )
    return positions


def _constrain_quaternions(adsorbate: AdsorbateDefinition, quats, tilt_deg: float):
    if adsorbate.natoms == 1:
        return np.tile(np.array([[1.0, 0.0, 0.0, 0.0]], dtype=float), (len(quats), 1))
    quats = np.asarray(quats, dtype=float).copy()
    if tilt_deg >= 179.0:
        return quats
    max_tilt = math.radians(tilt_deg)
    z_axis = np.array([0.0, 0.0, 1.0])
    for i, quat in enumerate(quats):
        rot = quaternion_to_matrix(quat)
        axis = rot @ z_axis
        angle = math.acos(np.clip(np.dot(axis, z_axis), -1.0, 1.0))
        if angle <= max_tilt:
            continue
        phi = math.atan2(axis[1], axis[0])
        half_tilt = 0.5 * max_tilt
        quats[i] = np.array(
            [
                math.cos(half_tilt),
                math.sin(half_tilt) * math.sin(phi),
                -math.sin(half_tilt) * math.cos(phi),
                0.0,
            ]
        )
    return quats


def build_pose_batch(
    density: DensityData,
    adsorbate: AdsorbateDefinition,
    samples_per_site: int = 48,
    systematic_orientations: int = 1,
    random_samples: int = 256,
    representative_sites_per_label: int = 1,
    seed: int = 12345,
    z_range=(1.5, 5.5),
    z_bias_power: float = 1.0,
    jitter_uv: float = 0.03,
    tilt_deg: float = 60.0,
):
    rng = np.random.default_rng(seed)
    raw_sites = infer_surface_sites(density)
    sites = []
    grouped = {}
    for site in raw_sites:
        grouped.setdefault(site.label, []).append(site)
    for label, label_sites in grouped.items():
        take = max(1, int(representative_sites_per_label))
        if len(label_sites) <= take:
            sites.extend(label_sites)
            continue
        indices = np.linspace(0, len(label_sites) - 1, take).round().astype(int)
        sites.extend([label_sites[i] for i in indices])
    planes = infer_reference_planes(density)
    systematic_z = _biased_linspace(max(2, int(samples_per_site)), z_range[0], z_range[1], max(float(z_bias_power), 1.0e-6))
    base_quats = _constrain_quaternions(
        adsorbate,
        _systematic_quaternions(adsorbate, systematic_orientations, tilt_deg=tilt_deg),
        tilt_deg=tilt_deg,
    )

    pose_rows = []
    cartesian = []
    labels = []
    for site in sites:
        for quat in base_quats:
            for z_height in systematic_z:
                uv = np.asarray(site.uv, dtype=float) + rng.normal(scale=jitter_uv, size=2)
                uv = uv - np.floor(uv)
                pose_rows.append(np.concatenate([uv, [z_height], quat]))
                cartesian.append(transform_adsorbate(adsorbate, density, uv, z_height, quat, planes["z_ref"]))
                labels.append(site.label)

    random_quats = _constrain_quaternions(
        adsorbate,
        random_quaternions(random_samples, rng),
        tilt_deg=tilt_deg,
    )
    for idx in range(random_samples):
        uv = rng.random(2)
        z_height = _sample_biased_uniform(rng, z_range[0], z_range[1], max(float(z_bias_power), 1.0e-6))
        quat = random_quats[idx]
        pose_rows.append(np.concatenate([uv, [z_height], quat]))
        cartesian.append(transform_adsorbate(adsorbate, density, uv, z_height, quat, planes["z_ref"]))
        labels.append("random")

    metadata = {
        "z_ref": planes["z_ref"],
        "z_image": planes["z_image"],
        "n_sites_total": len(raw_sites),
        "n_sites_used": len(sites),
        "surface_sites": [{"label": site.label, "uv": np.asarray(site.uv).tolist()} for site in sites],
    }
    return PoseBatch(
        adsorbate=adsorbate,
        positions=np.asarray(cartesian, dtype=float),
        pose_params=np.asarray(pose_rows, dtype=float),
        site_labels=labels,
        metadata=metadata,
    )
