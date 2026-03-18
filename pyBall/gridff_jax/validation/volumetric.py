"""Validation checks for CHGCAR/LOCPOT consistency."""

from __future__ import annotations

import numpy as np

from ..substrate_fields import poisson_potential_from_density


def validate_vasp_volumetrics(density):
    metrics = {
        "has_rho": density.rho_zyx is not None,
        "has_locpot": density.v_loc_zyx is not None,
        "grid_shape_zyx": list(density.grid_shape_zyx),
        "natoms": density.natoms,
    }
    if density.rho_zyx is not None:
        metrics["rho_integral"] = float(np.sum(density.rho_zyx))
    if density.rho_zyx is not None and density.v_loc_zyx is not None:
        poiss = poisson_potential_from_density(density.rho_zyx, density.voxel)
        lhs = density.v_loc_zyx.reshape(-1) - float(np.mean(density.v_loc_zyx))
        rhs = poiss.reshape(-1) - float(np.mean(poiss))
        denom = np.linalg.norm(lhs) * np.linalg.norm(rhs)
        corr = float(np.dot(lhs, rhs) / denom) if denom > 0.0 else 0.0
        metrics["locpot_poisson_corr"] = corr
        metrics["locpot_poisson_shifted_rmse"] = float(np.sqrt(np.mean((lhs - rhs) ** 2)))
    return metrics
