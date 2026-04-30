"""Modular GridFF distillation utilities with optional JAX acceleration."""

from .array_api import HAS_JAX, HAS_OPTAX, backend_summary
from .config import (
    AdsorbateConfig,
    DensityBackendConfig,
    DiagnosticsConfig,
    ExportConfig,
    FeatureToggles,
    GridConfig,
    HybridModelConfig,
    RunConfig,
    SurfaceConfig,
    TeacherBackendConfig,
    TrainingConfig,
    load_config,
    save_config,
)
from .config import Sampler6DConfig as Sampler6DConfigJSON
from .pbc import PBCInfo
from .interfaces import (
    AdsorbateDefinition,
    DensityData,
    DensityBackend,
    ModelEvaluation,
    PoseBatch,
    SurfaceSite,
    TeacherBackend,
    TeacherResult,
)

__all__ = [
    "HAS_JAX",
    "HAS_OPTAX",
    "backend_summary",
    "AdsorbateConfig",
    "DensityBackendConfig",
    "ExportConfig",
    "FeatureToggles",
    "GridConfig",
    "HybridModelConfig",
    "RunConfig",
    "SurfaceConfig",
    "TeacherBackendConfig",
    "TrainingConfig",
    "load_config",
    "save_config",
    "AdsorbateDefinition",
    "DensityData",
    "DensityBackend",
    "ModelEvaluation",
    "PoseBatch",
    "SurfaceSite",
    "TeacherBackend",
    "TeacherResult",
]
