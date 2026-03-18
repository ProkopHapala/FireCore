"""Synthetic teacher for end-to-end framework validation."""

from __future__ import annotations

import numpy as np

from ..hybrid_energy.model import HybridGridFFModel, default_hybrid_parameters
from ..interfaces import TeacherBackend, TeacherResult


class SyntheticTeacherBackend(TeacherBackend):
    def __init__(self, config, model_factory=None):
        self.config = config
        self.model_factory = model_factory

    def evaluate_batch(self, density, poses):
        elements = tuple(sorted(set(poses.adsorbate.symbols)))
        planes = poses.metadata
        params = default_hybrid_parameters(elements, z_image=planes["z_image"])
        for key in params.pauli:
            params.pauli[key] *= 1.08
            params.london[key] *= 0.92
            params.reactive[key] *= 1.10
        params.image_scale = 1.15
        params.reservoir_chi = -4.5
        params.reservoir_hardness = 7.0
        model = HybridGridFFModel(
            density=density,
            reactive_elements=elements,
            prefer_jax=False,
        )
        evaluation = model.evaluate_batch(poses, params=params)
        rng = np.random.default_rng(123)
        energies = evaluation.energies + rng.normal(scale=1.0e-4, size=evaluation.energies.shape)
        forces = evaluation.forces + rng.normal(scale=1.0e-4, size=evaluation.forces.shape)
        return TeacherResult(
            energies=energies,
            forces=forces,
            metadata={"backend": "synthetic", "elements": list(elements)},
        )
