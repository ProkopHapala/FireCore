from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path


@dataclass(frozen=True)
class ProtocolConfig:
    name: str
    gga: str
    ivdw: int
    ltssurf: bool
    lasph: bool
    encut: int
    ediff: float
    ediffg_relax: float
    ismear_metal: int
    sigma_metal: float
    ismear_molecule: int
    sigma_molecule: float


@dataclass(frozen=True)
class HPCConfig:
    queue: str
    nodes: int
    ncpus: int
    mem: str
    walltime: str
    email_notifications: bool
    module_lines: tuple[str, ...]
    mpi_provider: str
    vasp_binary_env: str
    vasp_binary_default: str
    ncore_bulk: int
    ncore_slab: int
    ncore_molecule: int


@dataclass(frozen=True)
class SlabCandidate:
    label: str
    size_xy: tuple[int, int]
    n_layers: int
    vacuum: float
    n_fixed_layers: int


@dataclass(frozen=True)
class OrientationSeed:
    label: str
    anchor_index: int
    mode: str
    tilt_deg: float
    azimuth_count: int
    initial_height: float


@dataclass(frozen=True)
class AdsorbateSpec:
    name: str
    formula: str
    representative: bool
    orientation_seeds: tuple[OrientationSeed, ...]
    scan_anchor_index: int


@dataclass(frozen=True)
class ScanGrid:
    z_min: float = 2.0
    z_mid: float = 3.5
    z_upper_mid: float = 6.0
    z_max: float = 12.0
    step_near: float = 0.10
    step_mid: float = 0.25
    step_far: float = 0.50


@dataclass(frozen=True)
class CampaignConfig:
    repo_root: Path
    external_root: Path
    orr_root: Path
    potcar_root: Path
    local_vasp_bin: Path
    metals: tuple[str, ...]
    representative_adsorbates: tuple[str, ...]
    production_adsorbates: tuple[str, ...]
    candidate_slabs: tuple[SlabCandidate, ...]
    primary_protocol: ProtocolConfig
    audit_protocol: ProtocolConfig
    hpc: HPCConfig
    scan_grid: ScanGrid
    bulk_seed_lattice: dict[str, float] = field(default_factory=dict)
    bulk_relax_stages: tuple[str, ...] = ("relax_stage1", "relax_stage2", "final_static")
    slab_relax_stage1_name: str = "relax_stage1_nodipole"
    slab_relax_stage2_cycles: int = 3
    slab_final_stage_name: str = "final_static"
    slab_workfunc_stage_name: str = "workfunction"
    slab_bader_stage_name: str = "bader"
    gas_relax_stages: tuple[str, ...] = ("relax_stage1", "final_static")
    site_labels: tuple[str, ...] = ("ontop", "bridge", "fcc", "hcp")
    campaign_phases: tuple[str, ...] = (
        "00_bulk",
        "01_universal_slab_selection",
        "02_clean_slab_final",
        "03_gas_refs",
        "04_ads_seed_library",
        "05_ads_relax",
        "06_scan_rigid",
        "07_scan_relaxed",
        "08_volumetrics",
        "09_analysis",
    )

    def to_json(self) -> dict:
        payload = asdict(self)
        payload["repo_root"] = str(self.repo_root)
        payload["external_root"] = str(self.external_root)
        payload["orr_root"] = str(self.orr_root)
        payload["potcar_root"] = str(self.potcar_root)
        payload["local_vasp_bin"] = str(self.local_vasp_bin)
        return payload


def build_default_campaign_config(repo_root: Path, external_root: Path | None = None) -> CampaignConfig:
    campaign_root = external_root or Path("/home/niel/git/coinage_gridff_dft")
    primary = ProtocolConfig(
        name="PBE+vdWsurf",
        gga="PE",
        ivdw=2,
        ltssurf=True,
        lasph=True,
        encut=500,
        ediff=1.0e-6,
        ediffg_relax=-0.02,
        ismear_metal=1,
        sigma_metal=0.10,
        ismear_molecule=0,
        sigma_molecule=0.01,
    )
    audit = ProtocolConfig(
        name="RPBE-D3(BJ)",
        gga="RP",
        ivdw=12,
        ltssurf=False,
        lasph=True,
        encut=500,
        ediff=1.0e-6,
        ediffg_relax=-0.02,
        ismear_metal=1,
        sigma_metal=0.10,
        ismear_molecule=0,
        sigma_molecule=0.01,
    )
    hpc = HPCConfig(
        queue="luna",
        nodes=1,
        ncpus=96,
        mem="1280gb",
        walltime="160:00:00",
        email_notifications=True,
        module_lines=(
            "source /cvmfs/software.metacentrum.cz/modulefiles/5.1.0/loadmodules",
            "module purge",
            "module load intel-parallel-studio/cluster.2019.4-intel-19.0.4-iifs5gt",
            "unset I_MPI_FABRICS",
            "export I_MPI_OFI_PROVIDER=sockets",
            "export I_MPI_OFI_LIBRARY_INTERNAL=1",
        ),
        mpi_provider="sockets",
        vasp_binary_env="VASP_BIN_ON_HPC",
        vasp_binary_default="vasp_std",
        ncore_bulk=16,
        ncore_slab=16,
        ncore_molecule=16,
    )
    return CampaignConfig(
        repo_root=repo_root,
        external_root=campaign_root,
        orr_root=Path("/home/niel/git/ORR_HER_Ag_Colab"),
        potcar_root=Path("/home/niel/src/VASP/pseudoP/potpaw_pbe"),
        local_vasp_bin=Path("/home/niel/src/VASP/vasp.6.3.2_gpu/bin/vasp_std"),
        metals=("Ag", "Cu", "Au"),
        representative_adsorbates=("H", "CO", "H2O", "HCONH2", "pyridine"),
        production_adsorbates=("H", "CO", "H2O", "NH3", "methanol", "HCONH2", "pyridine"),
        candidate_slabs=(
            SlabCandidate("3x3x4", (3, 3), 4, 18.0, 2),
        ),
        primary_protocol=primary,
        audit_protocol=audit,
        hpc=hpc,
        scan_grid=ScanGrid(),
        bulk_seed_lattice={
            "Ag": 4.09,
            "Cu": 3.63,
            "Au": 4.08,
        },
    )
