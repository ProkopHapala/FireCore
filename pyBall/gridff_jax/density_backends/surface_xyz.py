"""Surface XYZ backend for parity-style GridFF generation."""

from __future__ import annotations

import numpy as np

from pyBall import atomicUtils as au

from ..interfaces import DensityBackend, DensityData
from ..utils import voxel_spacing_from_cell


def _majority_top_layer_z(positions: np.ndarray) -> float:
    z_coords = np.asarray(positions[:, 2], dtype=float)
    unique_z, counts = np.unique(z_coords, return_counts=True)
    max_count = int(np.max(counts))
    return float(np.max(unique_z[counts == max_count]))


def _is_nice(n: int, allowed_factors=(2, 3, 5)) -> bool:
    value = int(n)
    for factor in allowed_factors:
        while value % factor == 0 and value > 1:
            value //= factor
    return value == 1


def _next_nice(n: int, allowed_factors=(2, 3, 5)) -> int:
    value = int(max(n, 1))
    while not _is_nice(value, allowed_factors):
        value += 1
    return value


def _adjust_dimensions(lengths, desired_voxel: float, allowed_factors=(2, 3, 5)):
    ns = []
    dg = []
    for length in lengths:
        target = int(np.ceil(float(length) / max(desired_voxel, 1.0e-6)))
        count = _next_nice(target, allowed_factors)
        ns.append(count)
        dg.append(float(length) / float(count))
    return ns, dg


class SurfaceXYZBackend(DensityBackend):
    def __init__(self, config, surface=None, grid=None):
        self.config = config
        self.surface = surface
        self.grid = grid

    def load(self) -> DensityData:
        if not self.config.xyz_path:
            raise ValueError("surface_xyz backend requires xyz_path")
        atoms = au.AtomicSystem(fname=self.config.xyz_path)
        positions = np.asarray(atoms.apos, dtype=float)
        charges = np.asarray(atoms.qs, dtype=float)
        symbols = [str(name) for name in atoms.enames]
        lvec = np.asarray(atoms.lvec, dtype=float)

        if self.grid is None:
            raise ValueError("surface_xyz backend requires grid configuration")

        if self.grid.surface_z0_mode == "top_layer_majority":
            z0 = _majority_top_layer_z(positions)
        elif self.grid.surface_z0_mode == "top_atom":
            z0 = float(np.max(positions[:, 2]))
        else:
            raise ValueError(f"Unsupported surface_z0_mode '{self.grid.surface_z0_mode}'")

        # Preserve the FULL lattice geometry (including off-diagonal entries).
        # Previously we collapsed to diag(lengths), which silently corrupted
        # skewed / stepped / rotated supercells. The voxel triplet is the
        # spacing-projected-along-each-lattice-vector, not the Cartesian
        # spacing — same convention as voxel_spacing_from_cell().
        lvec_full = np.asarray(lvec, dtype=float)
        lengths = tuple(float(np.linalg.norm(lvec_full[i])) for i in range(3))
        if self.config.grid_shape is None:
            ns, dg = _adjust_dimensions(lengths, float(self.grid.spacing))
            grid_shape_xyz = tuple(int(v) for v in ns)
        else:
            grid_shape_xyz = tuple(int(v) for v in self.config.grid_shape)
            dg = voxel_spacing_from_cell(lvec_full, grid_shape_xyz)
        cell = lvec_full
        # Origin: centre the cell laterally about the origin in the
        # diagonal/orthogonal case; for skewed cells the same centre-shift
        # along a and b vectors keeps the convention consistent.
        origin = (-0.5 * cell[0] - 0.5 * cell[1]).astype(float)
        origin[2] = z0
        voxel = np.asarray(dg, dtype=float)
        metadata = {
            "backend": "surface_xyz",
            "source_file": str(self.config.xyz_path),
            "grid_shape_xyz": grid_shape_xyz,
            "grid_shape_zyx": tuple(reversed(grid_shape_xyz)),
            "z0": float(z0),
            "grid_order": "zyx",
            "n_atoms": int(len(symbols)),
            "symbol_counts": {symbol: symbols.count(symbol) for symbol in sorted(set(symbols))},
        }
        return DensityData(
            cell=cell,
            origin=origin,
            voxel=voxel,
            symbols=symbols,
            positions=positions,
            charges=charges,
            rho_zyx=None,
            v_loc_zyx=None,
            grid_shape_xyz_hint=grid_shape_xyz,
            metadata=metadata,
        )
