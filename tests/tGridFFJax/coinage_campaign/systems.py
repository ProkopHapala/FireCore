from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from ase import Atoms
from ase.build import bulk, molecule

from .config import AdsorbateSpec, OrientationSeed


@dataclass(frozen=True)
class AdsorbatePlacement:
    atoms: Atoms
    anchor_index: int
    label: str
    mode: str
    tilt_deg: float
    azimuth_deg: float


def build_fcc_bulk(element: str, lattice_constant: float) -> Atoms:
    return bulk(element, "fcc", a=lattice_constant, cubic=True)


def _manual_formamide() -> Atoms:
    symbols = ["C", "O", "N", "H", "H", "H"]
    positions = np.array(
        [
            [0.00, 0.00, 0.00],
            [1.23, 0.00, 0.00],
            [-1.34, 0.00, 0.00],
            [0.00, 0.00, 1.10],
            [-1.87, 0.91, 0.00],
            [-1.87, -0.91, 0.00],
        ],
        dtype=float,
    )
    return Atoms(symbols=symbols, positions=positions)


def _manual_pyridine() -> Atoms:
    ring_radius = 1.39
    h_radius = 2.47
    ring_angles_deg = [90.0, 30.0, -30.0, -90.0, -150.0, 150.0]
    ring_symbols = ["N", "C", "C", "C", "C", "C"]
    positions = []
    symbols = []
    for angle_deg, symbol in zip(ring_angles_deg, ring_symbols):
        angle = math.radians(angle_deg)
        positions.append([ring_radius * math.cos(angle), ring_radius * math.sin(angle), 0.0])
        symbols.append(symbol)
    hydrogen_angles_deg = [30.0, -30.0, -90.0, -150.0, 150.0]
    for angle_deg in hydrogen_angles_deg:
        angle = math.radians(angle_deg)
        positions.append([h_radius * math.cos(angle), h_radius * math.sin(angle), 0.0])
        symbols.append("H")
    return Atoms(symbols=symbols, positions=np.array(positions, dtype=float))


def build_adsorbate(name: str) -> Atoms:
    if name == "H":
        return Atoms("H", positions=[[0.0, 0.0, 0.0]])
    if name == "CO":
        return molecule("CO")
    if name == "H2O":
        return molecule("H2O")
    if name == "NH3":
        return molecule("NH3")
    if name == "methanol":
        return molecule("CH3OH")
    if name in {"formamide", "HCONH2"}:
        return _manual_formamide()
    if name == "pyridine":
        return _manual_pyridine()
    raise ValueError(f"Unsupported adsorbate: {name}")


def default_adsorbate_specs() -> dict[str, AdsorbateSpec]:
    formamide_spec = AdsorbateSpec(
        name="HCONH2",
        formula="CH3NO",
        representative=False,
        orientation_seeds=(
            OrientationSeed("O_down", 1, "upright", 0.0, 6, 2.20),
            OrientationSeed("N_down", 2, "upright", 0.0, 6, 2.20),
            OrientationSeed("C_down", 0, "upright", 0.0, 6, 2.30),
            OrientationSeed("flat", 1, "flat", 0.0, 6, 2.60),
        ),
        scan_anchor_index=1,
    )
    return {
        "H": AdsorbateSpec(
            name="H",
            formula="H",
            representative=True,
            orientation_seeds=(
                OrientationSeed("H", 0, "upright", 0.0, 1, 1.10),
            ),
            scan_anchor_index=0,
        ),
        "CO": AdsorbateSpec(
            name="CO",
            formula="CO",
            representative=True,
            orientation_seeds=(
                OrientationSeed("C_down", 1, "upright", 0.0, 6, 1.80),
                OrientationSeed("O_down", 0, "upright", 0.0, 6, 1.90),
                OrientationSeed("flat", 1, "flat", 0.0, 6, 2.20),
            ),
            scan_anchor_index=1,
        ),
        "H2O": AdsorbateSpec(
            name="H2O",
            formula="H2O",
            representative=True,
            orientation_seeds=(
                OrientationSeed("O_down", 0, "upright", 0.0, 6, 2.10),
                OrientationSeed("H_down", 1, "upright", 0.0, 6, 2.10),
                OrientationSeed("flat", 0, "flat", 0.0, 6, 2.50),
            ),
            scan_anchor_index=0,
        ),
        "NH3": AdsorbateSpec(
            name="NH3",
            formula="H3N",
            representative=False,
            orientation_seeds=(
                OrientationSeed("N_down", 0, "upright", 0.0, 6, 2.10),
                OrientationSeed("H_down", 1, "upright", 0.0, 6, 2.20),
                OrientationSeed("flat", 0, "flat", 0.0, 6, 2.50),
            ),
            scan_anchor_index=0,
        ),
        "methanol": AdsorbateSpec(
            name="methanol",
            formula="CH4O",
            representative=False,
            orientation_seeds=(
                OrientationSeed("O_down", 1, "upright", 0.0, 6, 2.20),
                OrientationSeed("C_down", 0, "upright", 0.0, 6, 2.30),
                OrientationSeed("flat", 1, "flat", 0.0, 6, 2.60),
                OrientationSeed("OH_down", 3, "upright", 0.0, 6, 2.20),
            ),
            scan_anchor_index=1,
        ),
        "HCONH2": formamide_spec,
        "pyridine": AdsorbateSpec(
            name="pyridine",
            formula="C5H5N",
            representative=True,
            orientation_seeds=(
                OrientationSeed("N_down_upright", 0, "upright", 0.0, 6, 2.40),
                OrientationSeed("N_down_tilted", 0, "upright", 35.0, 6, 2.40),
                OrientationSeed("flat_ring", 0, "flat", 0.0, 6, 2.90),
                OrientationSeed("C_down", 1, "upright", 0.0, 6, 2.50),
            ),
            scan_anchor_index=0,
        ),
    }


def _rotation_from_vectors(source: np.ndarray, target: np.ndarray) -> np.ndarray:
    a = source / np.linalg.norm(source)
    b = target / np.linalg.norm(target)
    v = np.cross(a, b)
    c = float(np.dot(a, b))
    s = float(np.linalg.norm(v))
    if s < 1.0e-12 and c > 0.0:
        return np.eye(3)
    if s < 1.0e-12 and c < 0.0:
        trial = np.array([1.0, 0.0, 0.0])
        if abs(a[0]) > 0.8:
            trial = np.array([0.0, 1.0, 0.0])
        axis = np.cross(a, trial)
        axis = axis / np.linalg.norm(axis)
        cross = np.array(
            [
                [0.0, -axis[2], axis[1]],
                [axis[2], 0.0, -axis[0]],
                [-axis[1], axis[0], 0.0],
            ]
        )
        return np.eye(3) + 2.0 * cross @ cross
    vx = np.array(
        [
            [0.0, -v[2], v[1]],
            [v[2], 0.0, -v[0]],
            [-v[1], v[0], 0.0],
        ]
    )
    return np.eye(3) + vx + vx @ vx * ((1.0 - c) / (s * s))


def _rotation_about_axis(axis: np.ndarray, angle_deg: float) -> np.ndarray:
    angle = math.radians(angle_deg)
    axis = axis / np.linalg.norm(axis)
    cross = np.array(
        [
            [0.0, -axis[2], axis[1]],
            [axis[2], 0.0, -axis[0]],
            [-axis[1], axis[0], 0.0],
        ]
    )
    ident = np.eye(3)
    return ident + math.sin(angle) * cross + (1.0 - math.cos(angle)) * (cross @ cross)


def _plane_normal(positions: np.ndarray) -> np.ndarray:
    centered = positions - positions.mean(axis=0)
    _, _, vh = np.linalg.svd(centered, full_matrices=False)
    return vh[-1]


def orient_adsorbate(base_atoms: Atoms, anchor_index: int, mode: str, tilt_deg: float, azimuth_deg: float) -> Atoms:
    atoms = base_atoms.copy()
    positions = atoms.positions.copy()
    positions -= positions[anchor_index]
    rest = np.delete(positions, anchor_index, axis=0)
    if len(rest) == 0:
        atoms.positions = positions
        return atoms
    if mode == "flat":
        source_vector = _plane_normal(rest)
    else:
        source_vector = rest.mean(axis=0)
    if np.linalg.norm(source_vector) < 1.0e-10:
        source_vector = np.array([0.0, 0.0, 1.0])
    rotation = _rotation_from_vectors(source_vector, np.array([0.0, 0.0, 1.0]))
    positions = positions @ rotation.T
    if abs(tilt_deg) > 1.0e-12:
        positions = positions @ _rotation_about_axis(np.array([1.0, 0.0, 0.0]), tilt_deg).T
    if abs(azimuth_deg) > 1.0e-12:
        positions = positions @ _rotation_about_axis(np.array([0.0, 0.0, 1.0]), azimuth_deg).T
    atoms.positions = positions
    return atoms


def place_in_box(atoms: Atoms, box_length: float = 20.0) -> Atoms:
    boxed = atoms.copy()
    positions = boxed.positions.copy()
    positions -= positions.mean(axis=0)
    positions += box_length / 2.0
    boxed.positions = positions
    boxed.set_cell(np.diag([box_length, box_length, box_length]))
    boxed.pbc = True
    return boxed


def enumerate_seed_placements(spec: AdsorbateSpec) -> list[AdsorbatePlacement]:
    base = build_adsorbate(spec.name)
    placements: list[AdsorbatePlacement] = []
    for seed in spec.orientation_seeds:
        for idx in range(seed.azimuth_count):
            azimuth = 360.0 * idx / seed.azimuth_count
            oriented = orient_adsorbate(base, seed.anchor_index, seed.mode, seed.tilt_deg, azimuth)
            placements.append(
                AdsorbatePlacement(
                    atoms=oriented,
                    anchor_index=seed.anchor_index,
                    label=seed.label,
                    mode=seed.mode,
                    tilt_deg=seed.tilt_deg,
                    azimuth_deg=azimuth,
                )
            )
    return placements
