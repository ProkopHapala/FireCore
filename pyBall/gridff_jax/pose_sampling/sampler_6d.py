"""Generic 6D rigid-body pose sampler for molecule-on-surface configurations.

Generates systematic configurations spanning the full 6D space:
    (u, v, z, alpha, beta, gamma)

where (u, v) are fractional lateral coordinates over one surface unit cell,
z is the height above the surface reference plane, and (alpha, beta, gamma)
are Euler angles (or equivalently quaternions) for molecular orientation.

The sampler is substrate-agnostic and molecule-agnostic. It reads molecule
definitions from built-in registries, XYZ files, or a batch molecule list file.

Usage
-----
    from pyBall.gridff_jax.pose_sampling.sampler_6d import (
        Sampler6DConfig, build_6d_pose_batch, load_molecule_list,
    )

    config = Sampler6DConfig(n_u=10, n_v=10, n_z=20, n_orient=8)
    molecules = load_molecule_list("molecules.txt")
    for mol in molecules:
        batch = build_6d_pose_batch(density, mol, config)
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from ..interfaces import AdsorbateDefinition, DensityData, PoseBatch
from ..pbc import (
    PBCInfo,
    check_molecule_cell_compatibility,
    compute_minimum_teacher_tiling,
    detect_supercell_multiplicity,
    primitive_uv_grid,
    wrap_anchor_into_cell,
)
from ..utils import fractional_to_cartesian, normalize_quaternion, quaternion_to_matrix
from .molecules import benchmark_adsorbates, load_adsorbate_from_xyz
from .rigid import infer_reference_planes
from .sites import infer_surface_sites


# ---------------------------------------------------------------------------
#  Configuration
# ---------------------------------------------------------------------------

@dataclass
class ExtraDimension:
    """Placeholder for future dimensions beyond 6D (e.g., bond stretch, dihedral).

    Not sampled in the current version, but the data structure is in place
    so that adding molecular deformations later does not require refactoring.
    """
    name: str
    range: tuple[float, float] = (0.0, 1.0)
    n_points: int = 1
    sampling_mode: str = "linspace"  # "linspace", "random", "grid"
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass
class Sampler6DConfig:
    """Configuration for the 6D rigid-body pose sampler.

    Parameters
    ----------
    n_u, n_v : int
        Number of lateral grid points along fractional a and b directions.
        The grid covers [0, 1) in each direction (one unit cell).
    n_z : int
        Number of height points per lateral position.
    z_range : tuple[float, float]
        (z_min, z_max) in Angstrom above the surface reference plane.
    z_bias_power : float
        Bias exponent for z-sampling.  >1 concentrates points near z_min
        (close to the surface).  1.0 = uniform.
    n_orient : int
        Number of molecular orientations. For monoatomic probes this is
        forced to 1.  For polyatomics, Fibonacci-sphere quaternions + yaw
        are generated.
    tilt_max_deg : float
        Maximum tilt angle from the surface normal (degrees).  Orientations
        with tilt > tilt_max_deg are projected back.
    n_yaw : int
        Number of yaw angles (rotation around z-axis) per tilt direction.
        Total orientations = n_tilt_samples * n_yaw (capped at n_orient).
    include_high_symmetry_sites : bool
        If True, ensure top/bridge/hollow sites are explicitly included
        in the (u,v) grid (in addition to the regular grid).
    random_fraction : float
        Fraction of additional random poses to add on top of the systematic
        grid.  0.0 = no random poses, 0.2 = 20% extra random.
    seed : int
        Random seed for reproducibility.
    extra_dimensions : list[ExtraDimension]
        Future-proof hooks for molecular deformation DOFs.
    """
    n_u: int = 10
    n_v: int = 10
    n_z: int = 20
    z_range: tuple[float, float] = (1.5, 5.5)
    z_bias_power: float = 1.5
    n_orient: int = 8
    tilt_max_deg: float = 60.0
    n_yaw: int = 4
    include_high_symmetry_sites: bool = True
    random_fraction: float = 0.1
    seed: int = 42
    extra_dimensions: list[ExtraDimension] = field(default_factory=list)


# ---------------------------------------------------------------------------
#  Molecule list loader
# ---------------------------------------------------------------------------

def load_molecule_list(
    path: str | Path,
    default_anchor_index: int = 0,
    use_input_charges: bool = False,
) -> list[AdsorbateDefinition]:
    """Load molecules from a text file.

    File format (one molecule per line):
        # Comments start with '#'
        H                          # built-in name
        CHONH2                     # built-in name
        /path/to/mol.xyz           # XYZ file path
        /path/to/mol.xyz  name=CO  anchor=0  charges=true

    Returns a list of AdsorbateDefinition objects.
    """
    path = Path(path)
    builtins = benchmark_adsorbates()
    molecules = []

    with path.open("r", encoding="utf-8") as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split()
            token = parts[0]

            kwargs: dict[str, Any] = {
                "anchor_index": default_anchor_index,
                "use_input_charges": use_input_charges,
            }
            name_override = None
            for part in parts[1:]:
                if part.startswith("name="):
                    name_override = part.split("=", 1)[1]
                elif part.startswith("anchor="):
                    kwargs["anchor_index"] = int(part.split("=", 1)[1])
                elif part.startswith("charges="):
                    kwargs["use_input_charges"] = part.split("=", 1)[1].lower() in ("true", "1", "yes")

            if token in builtins:
                molecules.append(builtins[token])
            elif Path(token).suffix.lower() in (".xyz", ".mol2"):
                mol = load_adsorbate_from_xyz(
                    xyz_path=token,
                    name=name_override or Path(token).stem,
                    anchor_index=kwargs["anchor_index"],
                    use_input_charges=kwargs["use_input_charges"],
                )
                molecules.append(mol)
            else:
                if token in builtins:
                    molecules.append(builtins[token])
                else:
                    raise ValueError(
                        f"Unknown molecule '{token}': not a built-in name "
                        f"and not a recognized file path. Built-in names: "
                        f"{sorted(builtins.keys())}"
                    )

    return molecules


# ---------------------------------------------------------------------------
#  Orientation generators
# ---------------------------------------------------------------------------

def _axis_angle_quaternion(axis: np.ndarray, angle: float) -> np.ndarray:
    axis = np.asarray(axis, dtype=float)
    norm = np.linalg.norm(axis)
    if norm < 1.0e-12:
        return np.array([1.0, 0.0, 0.0, 0.0], dtype=float)
    axis = axis / norm
    half = 0.5 * float(angle)
    return np.array(
        [math.cos(half), axis[0] * math.sin(half),
         axis[1] * math.sin(half), axis[2] * math.sin(half)],
        dtype=float,
    )


def _quaternion_multiply(q1: np.ndarray, q2: np.ndarray) -> np.ndarray:
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    return np.array([
        w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
        w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
        w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
        w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2,
    ], dtype=float)


def _fibonacci_sphere_quaternions(n_points: int, tilt_max_deg: float) -> np.ndarray:
    """Generate approximately uniform orientations using Fibonacci sphere sampling.

    Produces quaternions whose rotation axis is distributed on a spherical cap
    defined by tilt_max_deg from the z-axis.
    """
    if n_points <= 1:
        return np.array([[1.0, 0.0, 0.0, 0.0]], dtype=float)

    tilt_max = math.radians(min(max(float(tilt_max_deg), 0.0), 179.0))
    cos_max = math.cos(tilt_max)

    golden_ratio = (1.0 + math.sqrt(5.0)) / 2.0
    quats = []
    for i in range(n_points):
        cos_theta = 1.0 - (1.0 - cos_max) * i / max(n_points - 1, 1)
        cos_theta = np.clip(cos_theta, -1.0, 1.0)
        theta = math.acos(cos_theta)
        phi = 2.0 * math.pi * i / golden_ratio

        tilt_axis = np.array([math.cos(phi), math.sin(phi), 0.0])
        quats.append(_axis_angle_quaternion(tilt_axis, theta))

    return np.array(quats, dtype=float)


def _generate_orientations(
    adsorbate: AdsorbateDefinition,
    n_orient: int,
    tilt_max_deg: float,
    n_yaw: int,
) -> np.ndarray:
    """Generate molecular orientations as quaternions.

    For monoatomic probes: returns identity quaternion only.
    For polyatomic molecules: Fibonacci sphere tilts combined with yaw rotations.
    """
    if adsorbate.natoms == 1:
        return np.array([[1.0, 0.0, 0.0, 0.0]], dtype=float)

    n_tilt = max(1, n_orient // max(n_yaw, 1))
    tilt_quats = _fibonacci_sphere_quaternions(n_tilt, tilt_max_deg)

    if n_yaw <= 1:
        return tilt_quats[:n_orient]

    yaw_angles = np.linspace(0.0, 2.0 * math.pi, n_yaw, endpoint=False)
    combined = []
    for tq in tilt_quats:
        for yaw in yaw_angles:
            yaw_q = _axis_angle_quaternion(np.array([0.0, 0.0, 1.0]), yaw)
            combined.append(_quaternion_multiply(yaw_q, tq))
            if len(combined) >= n_orient:
                break
        if len(combined) >= n_orient:
            break

    return np.array(combined[:n_orient], dtype=float)


# ---------------------------------------------------------------------------
#  Z-sampling
# ---------------------------------------------------------------------------

def _biased_z_grid(n_z: int, z_min: float, z_max: float, power: float) -> np.ndarray:
    if n_z <= 1:
        return np.array([z_min], dtype=float)
    t = np.linspace(0.0, 1.0, n_z, dtype=float)
    if abs(power - 1.0) > 1.0e-12:
        t = np.power(t, power)
    return z_min + (z_max - z_min) * t


# ---------------------------------------------------------------------------
#  Lateral grid
# ---------------------------------------------------------------------------

def _lateral_grid(
    n_u: int,
    n_v: int,
    density: DensityData,
    include_high_symmetry: bool,
) -> np.ndarray:
    """Build a (u,v) grid over the surface unit cell.

    If include_high_symmetry is True, top/bridge/hollow sites are inserted
    into the grid, merging with nearby regular grid points (within 0.02 frac).
    """
    u_vals = np.linspace(0.0, 1.0, n_u, endpoint=False)
    v_vals = np.linspace(0.0, 1.0, n_v, endpoint=False)
    uu, vv = np.meshgrid(u_vals, v_vals, indexing="ij")
    uv_grid = np.column_stack([uu.ravel(), vv.ravel()])

    if include_high_symmetry:
        sites = infer_surface_sites(density)
        site_uvs = np.array([np.asarray(s.uv, dtype=float) for s in sites])
        merge_tol = 0.5 / max(n_u, n_v)
        extras = []
        for suv in site_uvs:
            dists = np.linalg.norm(uv_grid - suv[None, :], axis=1)
            if dists.min() > merge_tol:
                extras.append(suv)
        if extras:
            uv_grid = np.vstack([uv_grid, np.array(extras, dtype=float)])

    return uv_grid


# ---------------------------------------------------------------------------
#  Core: build the 6D pose batch
# ---------------------------------------------------------------------------

def _transform_adsorbate(
    adsorbate: AdsorbateDefinition,
    density: DensityData,
    uv: np.ndarray,
    z_height: float,
    quaternion: np.ndarray,
    z_ref: float,
) -> np.ndarray:
    """Place adsorbate at (u,v,z) with given orientation, PBC-wrapped.

    The anchor atom is guaranteed to be inside the cell via PBC wrapping.
    Other atoms may extend outside the cell (this is correct — the grid
    interpolation wraps them independently via modulo indexing).

    Returns (natoms, 3) Cartesian positions.
    """
    quaternion = normalize_quaternion(quaternion)
    rot = quaternion_to_matrix(quaternion)
    local = adsorbate.positions - adsorbate.positions[adsorbate.anchor_index]
    rotated = local @ rot.T
    # Wrap (u,v) into [0, 1) before converting to Cartesian
    uv_wrapped = np.asarray(uv, dtype=float).copy()
    uv_wrapped[:2] = uv_wrapped[:2] - np.floor(uv_wrapped[:2])
    anchor_frac = np.array([uv_wrapped[0], uv_wrapped[1], 0.0], dtype=float)
    anchor_cart = fractional_to_cartesian(anchor_frac, density.cell)
    anchor_cart[2] = z_ref + z_height
    positions = rotated + anchor_cart[None, :]
    # Wrap anchor atom into cell (translates entire molecule by lattice vector)
    positions = wrap_anchor_into_cell(
        positions, adsorbate.anchor_index, density.cell, density.origin,
    )
    return positions


def build_6d_pose_batch(
    density: DensityData,
    adsorbate: AdsorbateDefinition,
    config: Sampler6DConfig | None = None,
) -> PoseBatch:
    """Generate a systematic 6D pose batch for one molecule on a surface.

    Parameters
    ----------
    density : DensityData
        Substrate density / geometry information.
    adsorbate : AdsorbateDefinition
        Molecule to sample.
    config : Sampler6DConfig, optional
        Sampling configuration.  Uses defaults if None.

    Returns
    -------
    PoseBatch
        Poses with full 6D parametrization in ``pose_params`` columns:
        [u, v, z_height, qw, qx, qy, qz].
        Metadata includes dimension labels, grid sizes, and site info.
    """
    if config is None:
        config = Sampler6DConfig()

    rng = np.random.default_rng(config.seed)
    planes = infer_reference_planes(density)
    z_ref = float(planes["z_ref"])

    # --- PBC: validate cell size and detect supercell ---
    pbc_info = PBCInfo.from_density(density, adsorbate)
    if pbc_info.warnings:
        import warnings as _w
        for msg in pbc_info.warnings:
            _w.warn(f"[6D sampler PBC] {msg}")

    teacher_tile = compute_minimum_teacher_tiling(adsorbate, density.cell)

    # --- lateral grid ---
    uv_grid = _lateral_grid(
        config.n_u, config.n_v, density,
        include_high_symmetry=config.include_high_symmetry_sites,
    )

    # --- z grid ---
    z_grid = _biased_z_grid(
        config.n_z, config.z_range[0], config.z_range[1], config.z_bias_power,
    )

    # --- orientations ---
    orientations = _generate_orientations(
        adsorbate, config.n_orient, config.tilt_max_deg, config.n_yaw,
    )

    # --- systematic grid: (uv) x z x orient ---
    pose_rows = []
    cartesian = []
    labels = []

    n_systematic = len(uv_grid) * len(z_grid) * len(orientations)
    for uv in uv_grid:
        for z_h in z_grid:
            for quat in orientations:
                pose_rows.append(np.concatenate([uv, [z_h], quat]))
                cartesian.append(
                    _transform_adsorbate(adsorbate, density, uv, z_h, quat, z_ref)
                )
                labels.append("grid")

    # --- random poses ---
    n_random = max(0, int(config.random_fraction * n_systematic))
    if n_random > 0:
        random_uvs = rng.random((n_random, 2))
        random_zs = config.z_range[0] + (
            config.z_range[1] - config.z_range[0]
        ) * np.power(rng.random(n_random), max(config.z_bias_power, 1.0e-6))

        from ..utils import random_quaternions
        random_quats = random_quaternions(n_random, rng)

        for i in range(n_random):
            uv_r = random_uvs[i]
            z_r = random_zs[i]
            q_r = random_quats[i]
            pose_rows.append(np.concatenate([uv_r, [z_r], q_r]))
            cartesian.append(
                _transform_adsorbate(adsorbate, density, uv_r, z_r, q_r, z_ref)
            )
            labels.append("random")

    # --- high-symmetry site labels ---
    if config.include_high_symmetry_sites:
        sites = infer_surface_sites(density)
        site_uv_map = {s.label: np.asarray(s.uv, dtype=float) for s in sites}
        merge_tol = 0.5 / max(config.n_u, config.n_v)
        for idx in range(len(labels)):
            if labels[idx] != "grid":
                continue
            uv_pose = pose_rows[idx][:2]
            for slabel, suv in site_uv_map.items():
                if np.linalg.norm(uv_pose - suv) < merge_tol:
                    labels[idx] = slabel
                    break

    metadata = {
        "z_ref": z_ref,
        "z_image": float(planes["z_image"]),
        "sampler": "6d_systematic",
        "config": {
            "n_u": config.n_u,
            "n_v": config.n_v,
            "n_z": config.n_z,
            "n_orient": config.n_orient,
            "z_range": list(config.z_range),
            "z_bias_power": config.z_bias_power,
            "tilt_max_deg": config.tilt_max_deg,
            "n_yaw": config.n_yaw,
            "random_fraction": config.random_fraction,
            "seed": config.seed,
        },
        "n_systematic": n_systematic,
        "n_random": n_random,
        "n_total": len(pose_rows),
        "n_lateral": len(uv_grid),
        "n_z_per_lateral": len(z_grid),
        "n_orientations": len(orientations),
        "dimensions": ["u", "v", "z_height", "qw", "qx", "qy", "qz"],
        "extra_dimensions": [
            {"name": ed.name, "range": list(ed.range), "n_points": ed.n_points}
            for ed in config.extra_dimensions
        ],
        # PBC metadata
        "pbc": {
            "cell_lengths_A": pbc_info.cell_lengths.tolist(),
            "lateral_area_A2": pbc_info.lateral_area,
            "supercell_multiplicity": list(pbc_info.supercell_multiplicity),
            "min_teacher_tile": list(teacher_tile),
            "pbc_warnings": pbc_info.warnings,
        },
    }

    return PoseBatch(
        adsorbate=adsorbate,
        positions=np.asarray(cartesian, dtype=float),
        pose_params=np.asarray(pose_rows, dtype=float),
        site_labels=labels,
        metadata=metadata,
    )


# ---------------------------------------------------------------------------
#  Convenience: build batches for multiple molecules
# ---------------------------------------------------------------------------

def build_6d_batches_from_file(
    density: DensityData,
    molecule_list_path: str | Path,
    config: Sampler6DConfig | None = None,
) -> dict[str, PoseBatch]:
    """Load molecules from a file and generate 6D batches for each.

    Returns a dict mapping molecule name -> PoseBatch.
    """
    molecules = load_molecule_list(molecule_list_path)
    batches = {}
    for mol in molecules:
        batches[mol.name] = build_6d_pose_batch(density, mol, config)
    return batches
