"""Adapter for external MAD-SURF-style teacher backends."""

from __future__ import annotations

import importlib
import math
import time

import numpy as np
from ase import Atoms

from ..interfaces import TeacherBackend, TeacherResult
from ..utils import quaternion_to_matrix


def _auto_tile(mol_positions: np.ndarray, cell: np.ndarray) -> tuple[int, int]:
    """Compute minimum tiling so the molecule fits without self-interaction.

    For each lateral axis (x, y), the required tile count is
    ``ceil(molecule_extent / cell_length) + 1`` where the +1 provides a
    safety margin so periodic images never overlap.

    Parameters
    ----------
    mol_positions : (N, 3) array
        Cartesian positions of the adsorbate atoms for a single pose.
    cell : (3, 3) array
        Unit cell vectors of the slab.

    Returns
    -------
    (n_tile_x, n_tile_y) : tuple[int, int]
    """
    cell = np.asarray(cell, dtype=float)
    mol_positions = np.asarray(mol_positions, dtype=float)
    tiles = []
    for axis in range(2):  # x, y only
        cell_length = float(np.linalg.norm(cell[axis]))
        mol_extent = float(mol_positions[:, axis].max() - mol_positions[:, axis].min())
        n = math.ceil(mol_extent / cell_length) + 1 if cell_length > 0.0 else 1
        tiles.append(max(n, 1))
    return (tiles[0], tiles[1])


class MADSurfTeacherBackend(TeacherBackend):
    def __init__(self, config, model_factory=None):
        self.config = config
        self.model_factory = model_factory
        self._callable = None
        self._calculator = None
        self._slab_energy = None
        self._slab_energy_tiled = {}  # keyed by (n_tile_x, n_tile_y)
        self._molecule_reference_cache = {}
        if model_factory is not None:
            self._callable = model_factory(config)
        elif config.model_module and config.model_callable:
            module = importlib.import_module(config.model_module)
            self._callable = getattr(module, config.model_callable)

    def _get_calculator(self):
        if self._calculator is not None:
            return self._calculator
        if self._callable is not None:
            return None
        if not self.config.model_path:
            raise RuntimeError("MAD-SURF backend requires teacher_backend.model_path when no callable is supplied")
        from mace.calculators import MACECalculator

        self._calculator = MACECalculator(
            model_paths=self.config.model_path,
            device=self.config.device,
            default_dtype=self.config.default_dtype,
        )
        return self._calculator

    @property
    def teacher_tile(self) -> tuple[int, int]:
        return tuple(self.config.teacher_tile) if hasattr(self.config, "teacher_tile") else (1, 1)

    def _slab_atoms(self, density):
        return Atoms(
            symbols=density.symbols,
            positions=np.asarray(density.positions, dtype=float),
            cell=np.asarray(density.cell, dtype=float),
            pbc=[True, True, True],
        )

    def _tiled_slab_atoms(self, density, n_tile_x: int, n_tile_y: int) -> Atoms:
        """Create a supercell by tiling the slab n_tile_x x n_tile_y times.

        The z-axis cell vector is unchanged; only the lateral (x, y)
        directions are replicated.
        """
        cell = np.asarray(density.cell, dtype=float)
        positions = np.asarray(density.positions, dtype=float)
        symbols = list(density.symbols)

        tiled_positions = []
        tiled_symbols = []
        for ix in range(n_tile_x):
            for iy in range(n_tile_y):
                shift = ix * cell[0] + iy * cell[1]
                tiled_positions.append(positions + shift)
                tiled_symbols.extend(symbols)

        big_cell = cell.copy()
        big_cell[0] = big_cell[0] * n_tile_x
        big_cell[1] = big_cell[1] * n_tile_y

        return Atoms(
            symbols=tiled_symbols,
            positions=np.vstack(tiled_positions),
            cell=big_cell,
            pbc=[True, True, True],
        )

    def _resolve_tile(self, density, mol_positions: np.ndarray) -> tuple[int, int]:
        """Determine the tiling to use for this evaluation.

        If teacher_tile is (0, 0) or "auto", compute the minimum tiling
        automatically from the molecule extent vs the cell size.
        Otherwise use the configured (n_tile_x, n_tile_y).
        """
        tx, ty = self.teacher_tile
        if tx <= 0 or ty <= 0:
            # auto mode
            return _auto_tile(mol_positions, np.asarray(density.cell, dtype=float))
        return (tx, ty)

    def _get_slab_energy(self, density, calc, n_tile_x: int, n_tile_y: int) -> float:
        """Return cached slab energy for the given tiling."""
        if n_tile_x == 1 and n_tile_y == 1:
            if self._slab_energy is None:
                slab = self._slab_atoms(density)
                slab.calc = calc
                self._slab_energy = float(slab.get_potential_energy())
            return self._slab_energy
        key = (n_tile_x, n_tile_y)
        if key not in self._slab_energy_tiled:
            slab = self._tiled_slab_atoms(density, n_tile_x, n_tile_y)
            slab.calc = calc
            self._slab_energy_tiled[key] = float(slab.get_potential_energy())
        return self._slab_energy_tiled[key]

    def _molecule_reference(self, poses):
        key = poses.adsorbate.name
        if key in self._molecule_reference_cache:
            return self._molecule_reference_cache[key]
        positions = np.asarray(poses.adsorbate.positions, dtype=float)
        positions = positions - positions.mean(axis=0, keepdims=True)
        cell = np.diag([40.0, 40.0, 40.0])
        atoms = Atoms(
            symbols=poses.adsorbate.symbols,
            positions=positions + np.array([20.0, 20.0, 20.0]),
            cell=cell,
            pbc=[True, True, True],
        )
        calc = self._get_calculator()
        atoms.calc = calc
        self._molecule_reference_cache[key] = {
            "energy": float(atoms.get_potential_energy()),
            "forces": np.asarray(atoms.get_forces(), dtype=float),
        }
        return self._molecule_reference_cache[key]

    def evaluate_batch(self, density, poses):
        if self._callable is not None:
            result = self._callable(density=density, poses=poses, config=self.config)
            if isinstance(result, TeacherResult):
                return result
            if isinstance(result, dict):
                return TeacherResult(
                    energies=result["energies"],
                    forces=result["forces"],
                    metadata=result.get("metadata", {}),
                )
            raise TypeError("MAD-SURF callable must return TeacherResult or a compatible dict")

        calc = self._get_calculator()
        molecule_reference = self._molecule_reference(poses) if self.config.interaction_energy else None

        # Determine tiling from the first pose (all poses share the same
        # adsorbate geometry so the extent is identical).
        n_tile_x, n_tile_y = self._resolve_tile(density, poses.positions[0])
        use_tiling = n_tile_x > 1 or n_tile_y > 1

        if use_tiling:
            tiled_slab = self._tiled_slab_atoms(density, n_tile_x, n_tile_y)
            slab_symbols = list(tiled_slab.get_chemical_symbols())
            slab_positions = np.asarray(tiled_slab.get_positions(), dtype=float)
            slab_cell = np.asarray(tiled_slab.get_cell(), dtype=float)
            # Shift to place molecules near the centre of the supercell
            # so they are far from all periodic boundaries.
            orig_cell = np.asarray(density.cell, dtype=float)
            centre_shift = 0.5 * (n_tile_x - 1) * orig_cell[0] + 0.5 * (n_tile_y - 1) * orig_cell[1]
        else:
            slab_symbols = list(density.symbols)
            slab_positions = np.asarray(density.positions, dtype=float)
            slab_cell = np.asarray(density.cell, dtype=float)
            centre_shift = np.zeros(3, dtype=float)

        slab_energy = self._get_slab_energy(density, calc, n_tile_x, n_tile_y)

        n_ads = poses.adsorbate.natoms
        energies = []
        forces = []
        t0 = time.perf_counter()
        for ipose, positions in enumerate(poses.positions):
            mol_positions = np.asarray(positions, dtype=float) + centre_shift
            combined = Atoms(
                symbols=slab_symbols + poses.adsorbate.symbols,
                positions=np.vstack([slab_positions, mol_positions]),
                cell=slab_cell,
                pbc=[True, True, True],
            )
            combined.calc = calc
            total_energy = float(combined.get_potential_energy())
            total_forces = np.asarray(combined.get_forces(), dtype=float)
            if self.config.interaction_energy:
                energies.append(total_energy - slab_energy - molecule_reference["energy"])
            else:
                energies.append(total_energy)
            ads_forces = total_forces[-n_ads:, :]
            if self.config.interaction_energy:
                quat = np.asarray(poses.pose_params[ipose, 3:7], dtype=float)
                rot = quaternion_to_matrix(quat)
                molecule_forces = molecule_reference["forces"] @ rot.T
                ads_forces = ads_forces - molecule_forces
            forces.append(ads_forces)
        elapsed = time.perf_counter() - t0
        return TeacherResult(
            energies=np.asarray(energies, dtype=float),
            forces=np.asarray(forces, dtype=float),
            metadata={
                "backend": "madsurf",
                "model_path": self.config.model_path,
                "device": self.config.device,
                "default_dtype": self.config.default_dtype,
                "interaction_energy": self.config.interaction_energy,
                "teacher_tile": [n_tile_x, n_tile_y],
                "teacher_eval_seconds": elapsed,
                "seconds_per_pose": elapsed / max(len(poses.positions), 1),
            },
        )
