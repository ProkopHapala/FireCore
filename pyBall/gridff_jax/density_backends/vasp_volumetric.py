"""VASP volumetric readers for CHGCAR and LOCPOT."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy import ndimage as sp_ndimage

from ..interfaces import DensityBackend, DensityData
from ..utils import voxel_spacing_from_cell


@dataclass
class ParsedVaspVolumetric:
    comment: str
    cell: np.ndarray
    symbols: list[str]
    counts: list[int]
    positions: np.ndarray
    grid_xyz: tuple[int, int, int]
    values_zyx: np.ndarray


def _parse_symbol_and_count_lines(symbols_line: str, counts_line: str):
    left = symbols_line.split()
    right = counts_line.split()
    if left and all(token.lstrip("+-").isdigit() for token in left):
        counts = [int(token) for token in left]
        symbols = [f"X{i + 1}" for i in range(len(counts))]
        mode_line = counts_line
    else:
        symbols = left
        counts = [int(token) for token in right]
        mode_line = None
    return symbols, counts, mode_line


def read_vasp_volumetric(path: str | Path) -> ParsedVaspVolumetric:
    path = Path(path)
    with path.open("r", encoding="utf8", errors="ignore") as fin:
        lines = fin.readlines()

    cursor = 0
    comment = lines[cursor].rstrip("\n")
    cursor += 1
    scale = float(lines[cursor].split()[0])
    cursor += 1
    cell = np.array(
        [list(map(float, lines[cursor + i].split()[:3])) for i in range(3)],
        dtype=float,
    )
    cursor += 3
    cell *= scale

    symbols, counts, mode_line = _parse_symbol_and_count_lines(lines[cursor], lines[cursor + 1])
    if mode_line is None:
        cursor += 2
    else:
        cursor += 1
    if mode_line is not None:
        coord_mode = mode_line.strip()
    else:
        coord_mode = lines[cursor].strip()
        cursor += 1
    if coord_mode.lower().startswith("s"):
        coord_mode = lines[cursor].strip()
        cursor += 1

    natoms = int(sum(counts))
    raw_positions = np.array(
        [list(map(float, lines[cursor + i].split()[:3])) for i in range(natoms)],
        dtype=float,
    )
    cursor += natoms
    if coord_mode.lower().startswith("d"):
        positions = raw_positions @ cell
    else:
        positions = raw_positions * scale

    while cursor < len(lines) and not lines[cursor].strip():
        cursor += 1
    grid_xyz = tuple(int(token) for token in lines[cursor].split()[:3])
    cursor += 1
    nvalues = int(np.prod(grid_xyz))
    payload = np.fromstring(" ".join(lines[cursor:]), sep=" ", count=nvalues)
    if payload.size != nvalues:
        raise ValueError(f"{path} ended before the full volumetric payload was read")
    values_zyx = payload.reshape((grid_xyz[2], grid_xyz[1], grid_xyz[0]))
    expanded_symbols = []
    for symbol, count in zip(symbols, counts):
        expanded_symbols.extend([symbol] * count)
    return ParsedVaspVolumetric(
        comment=comment,
        cell=cell,
        symbols=expanded_symbols,
        counts=counts,
        positions=positions,
        grid_xyz=grid_xyz,
        values_zyx=values_zyx,
    )


def _compare_cells(left, right, tol: float = 1.0e-6):
    return np.allclose(np.asarray(left), np.asarray(right), atol=tol, rtol=0.0)


def _resample_zyx(values_zyx: np.ndarray, target_xyz: tuple[int, int, int]) -> np.ndarray:
    target_zyx = (int(target_xyz[2]), int(target_xyz[1]), int(target_xyz[0]))
    factors = tuple(float(target_zyx[i]) / float(values_zyx.shape[i]) for i in range(3))
    return sp_ndimage.zoom(values_zyx, zoom=factors, order=1, mode="nearest")


class VaspVolumetricBackend(DensityBackend):
    def __init__(self, config, surface=None, grid=None):
        self.config = config
        self.surface = surface
        self.grid = grid

    def load(self) -> DensityData:
        if not self.config.chgcar_path:
            raise ValueError("vasp_volumetric backend requires chgcar_path")
        chg = read_vasp_volumetric(self.config.chgcar_path)
        locpot = None
        if self.config.locpot_path:
            locpot = read_vasp_volumetric(self.config.locpot_path)
            if chg.grid_xyz != locpot.grid_xyz:
                raise ValueError(
                    f"CHGCAR grid {chg.grid_xyz} does not match LOCPOT grid {locpot.grid_xyz}"
                )
            if not _compare_cells(chg.cell, locpot.cell):
                raise ValueError("CHGCAR and LOCPOT cells do not match")

        grid_xyz = chg.grid_xyz
        rho_zyx = chg.values_zyx
        v_loc_zyx = locpot.values_zyx if locpot is not None else None
        if self.config.grid_shape is not None:
            grid_xyz = tuple(int(v) for v in self.config.grid_shape)
            rho_zyx = _resample_zyx(rho_zyx, grid_xyz)
            if v_loc_zyx is not None:
                v_loc_zyx = _resample_zyx(v_loc_zyx, grid_xyz)
        voxel = voxel_spacing_from_cell(chg.cell, grid_xyz)
        metadata = {
            "backend": "vasp_volumetric",
            "source_files": {
                "chgcar": str(self.config.chgcar_path),
                "locpot": str(self.config.locpot_path) if self.config.locpot_path else None,
            },
            "grid_shape_xyz": grid_xyz,
            "grid_shape_zyx": tuple(reversed(grid_xyz)),
            "comment": chg.comment,
            "grid_order": "zyx",
            "n_atoms": len(chg.symbols),
            "symbol_counts": {symbol: chg.symbols.count(symbol) for symbol in sorted(set(chg.symbols))},
            "source_grid_xyz": tuple(int(v) for v in chg.grid_xyz),
        }
        return DensityData(
            cell=chg.cell,
            origin=np.zeros(3, dtype=float),
            voxel=voxel,
            symbols=chg.symbols,
            positions=chg.positions,
            charges=None,
            rho_zyx=rho_zyx,
            v_loc_zyx=v_loc_zyx,
            grid_shape_xyz_hint=tuple(int(v) for v in grid_xyz),
            metadata=metadata,
        )
