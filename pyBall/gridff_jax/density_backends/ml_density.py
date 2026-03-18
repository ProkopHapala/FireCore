"""Adapter backend for externally generated ML density grids."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from ..interfaces import DensityBackend, DensityData


class MLDensityBackend(DensityBackend):
    def __init__(self, config, surface=None, grid=None):
        self.config = config
        self.surface = surface
        self.grid = grid

    def load(self) -> DensityData:
        if not self.config.volumetric_npz_path:
            raise ValueError("ml_density backend requires volumetric_npz_path")
        path = Path(self.config.volumetric_npz_path)
        data = np.load(path, allow_pickle=True)
        metadata = data.get("metadata", None)
        if metadata is None:
            metadata = {}
        else:
            metadata = metadata.item()
        metadata = dict(metadata)
        metadata.setdefault("backend", "ml_density")
        return DensityData(
            cell=np.asarray(data["cell"], dtype=float),
            origin=np.asarray(data["origin"], dtype=float),
            voxel=np.asarray(data["voxel"], dtype=float),
            symbols=[str(v) for v in data["symbols"].tolist()],
            positions=np.asarray(data["positions"], dtype=float),
            charges=np.asarray(data["charges"], dtype=float) if "charges" in data else None,
            rho_zyx=np.asarray(data["rho_zyx"], dtype=float) if "rho_zyx" in data else None,
            v_loc_zyx=np.asarray(data["v_loc_zyx"], dtype=float) if "v_loc_zyx" in data else None,
            grid_shape_xyz_hint=tuple(int(v) for v in data["grid_shape_xyz"].tolist()) if "grid_shape_xyz" in data else None,
            metadata=metadata,
        )
