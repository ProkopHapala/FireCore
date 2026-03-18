"""Analytic fallback density backend for tests and bring-up."""

from __future__ import annotations

import numpy as np
from ase.build import fcc111

from ..interfaces import DensityBackend, DensityData
from ..utils import voxel_spacing_from_cell
from .vasp_volumetric import read_vasp_volumetric


_VALENCE = {
    "H": 1.0,
    "C": 4.0,
    "N": 5.0,
    "O": 6.0,
    "Cu": 11.0,
    "Ag": 11.0,
    "Au": 11.0,
}


def _make_grid(shape_xyz):
    nx, ny, nz = [int(v) for v in shape_xyz]
    zs, ys, xs = np.meshgrid(
        np.arange(nz, dtype=float),
        np.arange(ny, dtype=float),
        np.arange(nx, dtype=float),
        indexing="ij",
    )
    return xs, ys, zs


class PseudoDensityBackend(DensityBackend):
    def __init__(self, config, surface=None, grid=None):
        self.config = config
        self.surface = surface
        self.grid = grid

    def _build_surface(self):
        if self.config.metadata.get("template_chgcar"):
            parsed = read_vasp_volumetric(self.config.metadata["template_chgcar"])
            return parsed.cell, parsed.symbols, parsed.positions, parsed.grid_xyz

        metal = getattr(self.surface, "metal", "Ag")
        size = getattr(self.surface, "size", (6, 6, 4))
        vacuum = float(getattr(self.surface, "vacuum", 18.0))
        slab = fcc111(metal, size=size, vacuum=vacuum, orthogonal=False)
        cell = slab.cell.array
        positions = slab.positions
        symbols = slab.get_chemical_symbols()
        if self.config.grid_shape is not None:
            grid_xyz = tuple(int(v) for v in self.config.grid_shape)
        else:
            spacing = getattr(self.grid, "spacing", 0.25)
            lengths = [np.linalg.norm(cell[i]) for i in range(3)]
            grid_xyz = tuple(max(8, int(np.ceil(length / spacing))) for length in lengths)
        return cell, symbols, positions, grid_xyz

    def load(self) -> DensityData:
        cell, symbols, positions, grid_xyz = self._build_surface()
        voxel = voxel_spacing_from_cell(cell, grid_xyz)
        xs, ys, zs = _make_grid(grid_xyz)
        nx, ny, nz = grid_xyz
        rho = np.zeros((nz, ny, nx), dtype=float)
        sigma = float(self.config.gaussian_sigma)
        sigma2 = sigma * sigma
        inv_cell = np.linalg.inv(cell.T)
        for symbol, position in zip(symbols, positions):
            frac = inv_cell @ position
            gx = frac[0] * nx
            gy = frac[1] * ny
            gz = frac[2] * nz
            dx = np.minimum(np.abs(xs - gx), nx - np.abs(xs - gx))
            dy = np.minimum(np.abs(ys - gy), ny - np.abs(ys - gy))
            dz = zs - gz
            rr = (dx * voxel[0]) ** 2 + (dy * voxel[1]) ** 2 + (dz * voxel[2]) ** 2
            rho += _VALENCE.get(symbol, 4.0) * np.exp(-0.5 * rr / sigma2)

        with np.errstate(divide="ignore", invalid="ignore"):
            height_axis = np.arange(nz, dtype=float) * voxel[2]
            z0 = positions[:, 2].max() + 0.5
            v_loc = -np.exp(-np.maximum(height_axis - z0, 0.0) / 1.5)[:, None, None] * rho.max()
        metadata = {
            "backend": "pseudo_density",
            "grid_shape_xyz": grid_xyz,
            "grid_shape_zyx": tuple(reversed(grid_xyz)),
            "grid_order": "zyx",
            "generated": True,
        }
        return DensityData(
            cell=np.asarray(cell, dtype=float),
            origin=np.zeros(3, dtype=float),
            voxel=voxel,
            symbols=list(symbols),
            positions=np.asarray(positions, dtype=float),
            charges=np.zeros(len(symbols), dtype=float),
            rho_zyx=rho,
            v_loc_zyx=v_loc,
            grid_shape_xyz_hint=tuple(int(v) for v in grid_xyz),
            metadata=metadata,
        )
