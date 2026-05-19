"""Dataclass-based configuration for the hybrid GridFF workflow."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field, is_dataclass
import json
from pathlib import Path
from typing import Any


@dataclass
class FeatureToggles:
    use_density_charge: bool = True
    use_locpot: bool = True
    use_ct_qeq: bool = False
    use_image_charge: bool = False
    use_image_charge_fixed: bool = False
    use_reactive_grid: bool = False
    use_teacher_residual: bool = True


@dataclass
class SurfaceConfig:
    metal: str = "Ag"
    facet: str = "111"
    size: tuple[int, int, int] = (12, 12, 4)
    layers: int = 4
    vacuum: float = 18.0
    rigid: bool = True


@dataclass
class AdsorbateConfig:
    name: str = "H"
    xyz_path: str | None = None
    samples_per_site: int = 48
    systematic_orientations: int = 1
    random_samples: int = 256
    representative_sites_per_label: int = 1
    z_range: tuple[float, float] = (1.2, 5.5)
    z_bias_power: float = 1.0
    jitter_uv: float = 0.03
    tilt_deg: float = 60.0


@dataclass
class DensityBackendConfig:
    kind: str = "vasp_volumetric"
    xyz_path: str | None = None
    chgcar_path: str | None = None
    locpot_path: str | None = None
    volumetric_npz_path: str | None = None
    grid_shape: tuple[int, int, int] | None = None
    gaussian_sigma: float = 0.65
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass
class TeacherBackendConfig:
    kind: str = "synthetic"
    model_module: str | None = None
    model_callable: str | None = None
    model_path: str | None = None
    batch_size: int = 64
    device: str = "cpu"
    default_dtype: str = "float64"
    interaction_energy: bool = True
    teacher_tile: tuple[int, int] = (1, 1)
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass
class GridConfig:
    spacing: float = 0.25
    padding_z: float = 6.0
    interpolation_order: int = 3
    export_shape: str = "xyzc"
    builder_mode: str = "metal_density_plq"
    alpha_morse: float = 1.5
    morse_rcut: float = 20.0
    coulomb_rcut: float = 80.0
    pairwise_point_chunk: int = 2048
    pairwise_atom_chunk: int = 256
    element_types_path: str = "cpp/common_resources/ElementTypes.dat"
    req_scale_radius: dict[str, float] = field(default_factory=dict)
    req_scale_energy: dict[str, float] = field(default_factory=dict)
    surface_z0_mode: str = "top_layer_majority"
    reactive_sigma_z: float = 1.2
    reactive_power: float = 1.25
    pauli_power: float = 1.5
    london_power: float = 0.75
    metal_density_pauli_power: float = 1.0
    metal_density_london_power: float = 0.5
    metal_density_rho_smoothing_sigma: float = 0.0
    metal_bulk_electron_density: float = 0.0530
    london_damping_d0: float = 0.0
    london_damping_width: float = 0.5
    use_pairwise_c6: bool = False
    c6_rcut: float = 15.0
    c6_s_r: float = 0.94
    # DISABLED(2026-03): density-derived auxiliary channels confirmed ineffective
    # in exhaustive 10-strategy ablation scan. All produce identical or worse RMSE
    # vs 37.1 meV PLQ baseline. Laplacian/gradient interfere with REQ optimization;
    # alt London is absorbed by existing REQ parameters. Code preserved for reference.
    use_density_laplacian: bool = False
    use_density_gradient: bool = False
    use_london_alt: bool = False
    london_alt_power: float = 1.5


@dataclass
class HybridModelConfig:
    """Hybrid-model configuration.

    .. note::
       The flags ``use_qeq``, ``use_image``, ``use_reactive_grid`` are
       **informational** — the runtime energy assembly in
       ``hybrid_energy/model.py`` reads ``FeatureToggles`` (see
       ``RunConfig.toggles``), NOT these fields. Only ``use_req_plq`` is read
       directly by the model (controls REQ-PLQ parameterisation). The flags
       are still set by ``apply_substrate_class()`` so that downstream
       exporters and provenance logs reflect the substrate-class choice
       consistently.
    """
    use_qeq: bool = False
    use_image: bool = False
    use_reactive_grid: bool = False
    use_req_plq: bool = False
    reactive_channels: tuple[str, ...] = ("H", "C", "O")
    reservoir_charge: float = 0.0
    image_plane_offset: float = 1.0
    image_damping: float = 0.5


@dataclass
class TrainingConfig:
    """Training-loop configuration.

    .. note::
       ``batch_size`` is currently **unused** — the JAX optimiser packs the
       full split into one JIT'd loss call (see ``fit/optimize.py``). For
       very large 6D datasets this can OOM the device; use
       ``--max-abs-teacher-eV`` to drop catastrophic outliers and/or reduce
       the pose-grid size. The field is retained for forward compatibility
       in case mini-batched training is added later.
    """
    seed: int = 12345
    batch_size: int = 64    # informational — see class docstring
    learning_rate: float = 1.0e-2
    max_steps: int = 400
    early_stopping_patience: int = 40
    force_weight: float = 10.0
    energy_weight: float = 1.0
    validation_split: float = 0.1
    test_split: float = 0.1
    fit_chi: bool = False
    fit_hardness: bool = False
    fit_image_plane: bool = False
    fit_reactive: bool = False
    fit_static_charge: bool = False
    fit_req_radius_offset: bool = False
    fit_req_energy_scale: bool = False
    fit_sample_shift_z: bool = False
    fit_coulomb_shift_z: bool = False
    fit_c6_coeff: bool = False
    # DISABLED(2026-03): see GridConfig comments — these channels were confirmed ineffective.
    fit_density_lap: bool = False
    fit_density_grad: bool = False
    fit_london_alt: bool = False
    req_radius_regularization: float = 5.0e-2
    req_energy_regularization: float = 5.0e-3
    c6_regularization: float = 1.0e-1
    # DISABLED(2026-03): regularization for disabled channels (kept for compatibility).
    density_lap_regularization: float = 0.0
    density_grad_regularization: float = 0.0
    london_alt_regularization: float = 0.0
    constraint_bound_fraction: float = 5.0e-2
    # Two-stage C₆ optimization: Stage 1 fits REQ (c6_coeff=0), Stage 2 freezes REQ and fits C₆
    two_stage_c6: bool = False
    stage2_max_steps: int = 400
    stage2_learning_rate: float = 5.0e-3
    stage2_force_weight: float = 10.0
    # Stage 3 (optional): gentle joint REQ+C₆ refinement with small learning rate
    stage3_refine: bool = False
    stage3_max_steps: int = 200
    stage3_learning_rate: float = 1.0e-4


@dataclass
class ExportConfig:
    out_dir: str = "tests/tGridFFJax/ag_hybrid_gridff"
    write_plq: bool = True
    write_sidecars: bool = True
    write_plots: bool = True


@dataclass
class Sampler6DConfig:
    """Configuration for the 6D rigid-body pose sampler.

    This mirrors ``pose_sampling.sampler_6d.Sampler6DConfig`` but lives here
    for JSON serialization alongside the rest of RunConfig.
    """
    n_u: int = 10
    n_v: int = 10
    n_z: int = 20
    z_range: tuple[float, float] = (1.5, 5.5)
    z_bias_power: float = 1.5
    n_orient: int = 8
    tilt_max_deg: float = 60.0
    n_yaw: int = 4
    include_high_symmetry_sites: bool = True
    random_fraction: float = 0.1
    seed: int = 42
    molecule_list_path: str | None = None


@dataclass
class DiagnosticsConfig:
    """Configuration for MLIP diagnostic tools."""
    fd_step_angstrom: float = 0.005
    fd_n_poses: int = 50
    smoothness_z_heights: list[float] = field(
        default_factory=lambda: [2.0, 2.5, 3.0, 3.5, 4.0]
    )
    smoothness_n_grid: int = 25
    corrugation_z_heights: list[float] = field(
        default_factory=lambda: [2.0, 2.5, 3.0, 4.0, 5.0]
    )
    corrugation_n_grid: int = 25
    report_format: str = "json"


@dataclass
class RunConfig:
    surface: SurfaceConfig = field(default_factory=SurfaceConfig)
    adsorbates: list[AdsorbateConfig] = field(
        default_factory=lambda: [
            AdsorbateConfig(name="H", z_range=(1.2, 5.0)),
            AdsorbateConfig(name="CO", z_range=(1.5, 5.5)),
            AdsorbateConfig(name="H2O", z_range=(1.5, 5.5)),
        ]
    )
    density_backend: DensityBackendConfig = field(default_factory=DensityBackendConfig)
    teacher_backend: TeacherBackendConfig = field(default_factory=TeacherBackendConfig)
    toggles: FeatureToggles = field(default_factory=FeatureToggles)
    grid: GridConfig = field(default_factory=GridConfig)
    hybrid_model: HybridModelConfig = field(default_factory=HybridModelConfig)
    training: TrainingConfig = field(default_factory=TrainingConfig)
    export: ExportConfig = field(default_factory=ExportConfig)
    sampler_6d: Sampler6DConfig = field(default_factory=Sampler6DConfig)
    diagnostics: DiagnosticsConfig = field(default_factory=DiagnosticsConfig)


def _coerce_dataclass(cls, payload: dict[str, Any]):
    values = {}
    for name, field_info in cls.__dataclass_fields__.items():
        if name not in payload:
            continue
        value = payload[name]
        ftype = field_info.type
        default = getattr(cls(), name) if callable(cls) else None
        if is_dataclass(default):
            values[name] = _coerce_dataclass(type(default), value)
        elif isinstance(default, list) and default and is_dataclass(default[0]):
            sub_cls = type(default[0])
            values[name] = [_coerce_dataclass(sub_cls, item) for item in value]
        else:
            values[name] = value
    return cls(**values)


def load_config(path: str | Path) -> RunConfig:
    path = Path(path)
    with path.open("r", encoding="utf8") as fin:
        payload = json.load(fin)
    return _coerce_dataclass(RunConfig, payload)


def save_config(config: RunConfig, path: str | Path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf8") as fout:
        json.dump(asdict(config), fout, indent=2)
