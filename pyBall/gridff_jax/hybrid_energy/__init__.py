"""Hybrid student energy for metal-aware GridFF distillation."""

from __future__ import annotations

from .qeq import solve_qeq_with_reservoir
from .reactive import DensityReactiveModel

__all__ = ["HybridGridFFModel", "HybridParameters", "default_hybrid_parameters", "solve_qeq_with_reservoir", "DensityReactiveModel"]


def __getattr__(name: str):
    if name in {"HybridGridFFModel", "HybridParameters", "default_hybrid_parameters"}:
        from .model import HybridGridFFModel, HybridParameters, default_hybrid_parameters

        values = {
            "HybridGridFFModel": HybridGridFFModel,
            "HybridParameters": HybridParameters,
            "default_hybrid_parameters": default_hybrid_parameters,
        }
        return values[name]
    raise AttributeError(f"module '{__name__}' has no attribute '{name}'")
