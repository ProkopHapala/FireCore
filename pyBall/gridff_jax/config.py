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


@dataclass
class HybridModelConfig:
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
    seed: int = 12345
    batch_size: int = 64
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
    req_radius_regularization: float = 5.0e-2
    req_energy_regularization: float = 5.0e-3
    constraint_bound_fraction: float = 5.0e-2


@dataclass
class ExportConfig:
    out_dir: str = "tests/tGridFFJax/ag_hybrid_gridff"
    write_plq: bool = True
    write_sidecars: bool = True
    write_plots: bool = True


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
