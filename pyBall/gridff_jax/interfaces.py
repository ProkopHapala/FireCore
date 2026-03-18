"""Shared data structures and backend interfaces."""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Any

import numpy as np


@dataclass
class DensityData:
    cell: np.ndarray
    origin: np.ndarray
    voxel: np.ndarray
    symbols: list[str]
    positions: np.ndarray
    charges: np.ndarray | None = None
    rho_zyx: np.ndarray | None = None
    v_loc_zyx: np.ndarray | None = None
    grid_shape_xyz_hint: tuple[int, int, int] | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    @property
    def natoms(self) -> int:
        return len(self.symbols)

    @property
    def grid_shape_zyx(self) -> tuple[int, int, int]:
        array = self.v_loc_zyx if self.v_loc_zyx is not None else self.rho_zyx
        if array is not None:
            return tuple(int(v) for v in array.shape[:3])
        if self.grid_shape_xyz_hint is None:
            raise ValueError("DensityData does not contain rho_zyx, v_loc_zyx, or grid_shape_xyz_hint")
        nx, ny, nz = self.grid_shape_xyz_hint
        return (int(nz), int(ny), int(nx))

    @property
    def grid_shape_xyz(self) -> tuple[int, int, int]:
        nz, ny, nx = self.grid_shape_zyx
        return (nx, ny, nz)


@dataclass
class AdsorbateDefinition:
    name: str
    symbols: list[str]
    positions: np.ndarray
    charges: np.ndarray | None = None
    anchor_index: int = 0
    metadata: dict[str, Any] = field(default_factory=dict)

    @property
    def natoms(self) -> int:
        return len(self.symbols)


@dataclass
class SurfaceSite:
    label: str
    uv: np.ndarray
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass
class PoseBatch:
    adsorbate: AdsorbateDefinition
    positions: np.ndarray
    pose_params: np.ndarray
    site_labels: list[str]
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass
class TeacherResult:
    energies: np.ndarray
    forces: np.ndarray
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass
class ModelEvaluation:
    energies: np.ndarray
    forces: np.ndarray
    components: dict[str, np.ndarray] = field(default_factory=dict)
    charges: np.ndarray | None = None


class DensityBackend(ABC):
    @abstractmethod
    def load(self) -> DensityData:
        raise NotImplementedError


class TeacherBackend(ABC):
    @abstractmethod
    def evaluate_batch(
        self,
        density: DensityData,
        poses: PoseBatch,
    ) -> TeacherResult:
        raise NotImplementedError


class ReactiveModel(ABC):
    @abstractmethod
    def build_channels(
        self,
        density: DensityData,
        reactive_elements: tuple[str, ...],
    ) -> np.ndarray:
        raise NotImplementedError
