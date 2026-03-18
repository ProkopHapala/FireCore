"""Reactive short-range grid channels derived from the slab density."""

from __future__ import annotations

import numpy as np

from ..interfaces import ReactiveModel


class DensityReactiveModel(ReactiveModel):
    def __init__(self, z_ref: float, voxel, power: float = 1.25, sigma_z: float = 1.2):
        self.z_ref = float(z_ref)
        self.voxel = np.asarray(voxel, dtype=float)
        self.power = float(power)
        self.sigma_z = float(sigma_z)

    def build_channels(self, density, reactive_elements):
        if density.rho_zyx is None:
            raise ValueError("Reactive channels require rho_zyx")
        rho = np.asarray(density.rho_zyx, dtype=float)
        rho_norm = rho / max(rho.max(), 1.0e-12)
        nz = rho.shape[0]
        z_axis = density.origin[2] + np.arange(nz, dtype=float) * self.voxel[2]
        height = np.exp(-0.5 * ((z_axis - self.z_ref) / max(self.sigma_z, 1.0e-6)) ** 2)
        base = np.power(np.clip(rho_norm, 0.0, None), self.power) * height[:, None, None]
        return np.stack([base.copy() for _ in reactive_elements], axis=-1)
