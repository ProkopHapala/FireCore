"""Adapter for external MAD-SURF-style teacher backends."""

from __future__ import annotations

import importlib
import time

import numpy as np
from ase import Atoms

from ..interfaces import TeacherBackend, TeacherResult
from ..utils import quaternion_to_matrix


class MADSurfTeacherBackend(TeacherBackend):
    def __init__(self, config, model_factory=None):
        self.config = config
        self.model_factory = model_factory
        self._callable = None
        self._calculator = None
        self._slab_energy = None
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

    def _slab_atoms(self, density):
        return Atoms(
            symbols=density.symbols,
            positions=np.asarray(density.positions, dtype=float),
            cell=np.asarray(density.cell, dtype=float),
            pbc=[True, True, True],
        )

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
            pbc=[False, False, False],
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
        slab = self._slab_atoms(density)
        if self._slab_energy is None:
            slab.calc = calc
            self._slab_energy = float(slab.get_potential_energy())
        molecule_reference = self._molecule_reference(poses) if self.config.interaction_energy else None

        n_ads = poses.adsorbate.natoms
        energies = []
        forces = []
        t0 = time.perf_counter()
        for ipose, positions in enumerate(poses.positions):
            combined = Atoms(
                symbols=density.symbols + poses.adsorbate.symbols,
                positions=np.vstack([density.positions, positions]),
                cell=np.asarray(density.cell, dtype=float),
                pbc=[True, True, True],
            )
            combined.calc = calc
            total_energy = float(combined.get_potential_energy())
            total_forces = np.asarray(combined.get_forces(), dtype=float)
            if self.config.interaction_energy:
                energies.append(total_energy - self._slab_energy - molecule_reference["energy"])
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
                "teacher_eval_seconds": elapsed,
                "seconds_per_pose": elapsed / max(len(poses.positions), 1),
            },
        )
