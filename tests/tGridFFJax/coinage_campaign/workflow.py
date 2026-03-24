from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.build import add_adsorbate, fcc111
from ase.constraints import FixAtoms
from ase.io import read as ase_read
from ase.io import write as ase_write

from .config import AdsorbateSpec, CampaignConfig, SlabCandidate
from .systems import build_adsorbate, build_fcc_bulk, default_adsorbate_specs, enumerate_seed_placements, place_in_box
from .templates import (
    bader_incar,
    bulk_relax_incar,
    final_static_incar,
    gas_phase_incar,
    hpc_pbs_template,
    hpc_submit_script,
    kpoints_grid,
    local_run_script,
    relax_pipeline_script,
    relaxed_scan_incar,
    rigid_scan_incar,
    run_custodian_template,
    slab_relax_incar,
    workfunction_incar,
)


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def save_json(payload: dict, path: Path) -> None:
    ensure_dir(path.parent)
    path.write_text(json.dumps(payload, indent=2) + "\n")


def species_order(atoms: Atoms) -> list[str]:
    seen: list[str] = []
    for symbol in atoms.get_chemical_symbols():
        if symbol not in seen:
            seen.append(symbol)
    return seen


def assemble_potcar(outdir: Path, atoms: Atoms, potcar_root: Path) -> None:
    content = []
    for symbol in species_order(atoms):
        potcar_path = potcar_root / symbol / "POTCAR"
        content.append(potcar_path.read_text())
    (outdir / "POTCAR").write_text("".join(content))


def write_job_dir(
    job_dir: Path,
    atoms: Atoms,
    incar_text: str,
    kpoints_text: str,
    config: CampaignConfig,
    metadata: dict,
    local_job_name: str,
    hpc_ncore: int,
    hpc_ncpus: int | None = None,
    hpc_mem: str | None = None,
    hpc_scratch: str | None = None,
    hpc_binary_default: str | None = None,
    hpc_binary_env: str | None = None,
    use_custodian_hpc: bool = True,
) -> None:
    ensure_dir(job_dir)
    ase_write(job_dir / "POSCAR", atoms, format="vasp", direct=True, sort=False)
    (job_dir / "INCAR").write_text(incar_text)
    (job_dir / "KPOINTS").write_text(kpoints_text)
    assemble_potcar(job_dir, atoms, config.potcar_root)
    (job_dir / "run_local.sh").write_text(local_run_script(config.local_vasp_bin))
    (job_dir / "run.pbs").write_text(
        hpc_pbs_template(
            local_job_name,
            config.hpc,
            ncpus=hpc_ncpus,
            mem=hpc_mem,
            scratch_local=hpc_scratch,
            vasp_binary_default=hpc_binary_default,
            use_custodian=use_custodian_hpc,
        )
    )
    (job_dir / "run_custodian.py").write_text(
        run_custodian_template(
            config.hpc,
            hpc_ncore,
            binary_env_name=hpc_binary_env,
            binary_default=hpc_binary_default,
        )
    )
    (job_dir / "submit_hpc.sh").write_text(hpc_submit_script())
    (job_dir / "run_local.sh").chmod(0o755)
    (job_dir / "run.pbs").chmod(0o755)
    (job_dir / "run_custodian.py").chmod(0o755)
    (job_dir / "submit_hpc.sh").chmod(0o755)
    save_json(metadata, job_dir / "job_manifest.json")


def default_bulk_kmesh() -> tuple[int, int, int]:
    return (21, 21, 21)


def bulk_final_kmesh() -> tuple[int, int, int]:
    return (21, 21, 21)


def slab_kmesh(candidate: SlabCandidate) -> tuple[int, int, int]:
    if candidate.size_xy == (2, 2):
        return (6, 6, 1)
    if candidate.n_layers >= 5:
        return (4, 4, 1)
    return (4, 4, 1)


def slab_dense_kmesh(candidate: SlabCandidate) -> tuple[int, int, int]:
    relax_mesh = slab_kmesh(candidate)
    return (relax_mesh[0] * 3, relax_mesh[1] * 3, 1)


def gas_kmesh() -> tuple[int, int, int]:
    return (1, 1, 1)


def build_slab(element: str, lattice_constant: float, candidate: SlabCandidate) -> tuple[Atoms, list[int]]:
    slab = fcc111(
        element,
        size=(candidate.size_xy[0], candidate.size_xy[1], candidate.n_layers),
        a=lattice_constant,
        vacuum=candidate.vacuum,
        periodic=True,
    )
    z_positions = slab.positions[:, 2]
    unique_z = sorted(set(np.round(z_positions, 4)))
    z_cutoff = unique_z[candidate.n_fixed_layers - 1] + 1.0e-3
    fixed = [i for i, z in enumerate(z_positions) if z <= z_cutoff]
    slab.set_constraint(FixAtoms(indices=fixed))
    return slab, fixed


def lattice_from_concar(contcar_path: Path) -> float:
    atoms = ase_read(contcar_path)
    lengths = atoms.cell.lengths()
    if not np.allclose(lengths[0], lengths[1], atol=1.0e-6) or not np.allclose(lengths[1], lengths[2], atol=1.0e-6):
        raise ValueError(f"Expected cubic bulk cell in {contcar_path}, got lengths {lengths}")
    return float(lengths[0])


def resolve_bulk_lattice_constant(
    config: CampaignConfig,
    campaign_root: Path,
    metal: str,
    allow_reference_fallback: bool = False,
) -> tuple[float, str]:
    bulk_converged = campaign_root / "00_bulk" / metal / "final_static" / "CONTCAR"
    if bulk_converged.exists():
        return lattice_from_concar(bulk_converged), str(bulk_converged)
    if allow_reference_fallback:
        orr_converged = config.orr_root / f"{metal}_bulk" / "CONTCAR"
        if orr_converged.exists():
            return lattice_from_concar(orr_converged), str(orr_converged)
        return config.bulk_seed_lattice[metal], "config.bulk_seed_lattice"
    raise FileNotFoundError(f"Missing converged bulk CONTCAR for {metal}: {bulk_converged}")


def annotate_atoms(atoms: Atoms, adsorbate_name: str, anchor_index: int | None = None) -> Atoms:
    tagged = atoms.copy()
    symbols = tagged.get_chemical_symbols()
    ads_mask = np.zeros(len(symbols), dtype=int)
    metal_symbols = {"Ag", "Cu", "Au"}
    ads_mask[[i for i, sym in enumerate(symbols) if sym not in metal_symbols]] = 1
    tagged.set_array("gridff_adsorbate_mask", ads_mask)
    if anchor_index is not None:
        anchor = np.full(len(symbols), -1, dtype=int)
        anchor[len(symbols) - len([i for i, s in enumerate(symbols) if s not in metal_symbols]) + anchor_index] = 1
        tagged.set_array("gridff_anchor_marker", anchor)
    tagged.info["gridff_adsorbate"] = adsorbate_name
    return tagged


def place_adsorbate_on_slab(
    slab: Atoms,
    adsorbate_spec: AdsorbateSpec,
    placement,
    site_label: str,
    initial_height: float,
) -> tuple[Atoms, dict]:
    combined = slab.copy()
    add_adsorbate(
        combined,
        placement.atoms,
        height=initial_height,
        position=site_label,
        mol_index=placement.anchor_index,
    )
    slab_indices = list(range(len(slab)))
    combined.set_constraint(FixAtoms(indices=slab_indices))
    metadata = {
        "adsorbate": adsorbate_spec.name,
        "placement_label": placement.label,
        "orientation_mode": placement.mode,
        "tilt_deg": placement.tilt_deg,
        "azimuth_deg": placement.azimuth_deg,
        "anchor_index_local": placement.anchor_index,
        "site_label": site_label,
        "initial_height": initial_height,
    }
    return combined, metadata


def rigid_scan_values(config: CampaignConfig) -> list[float]:
    grid = config.scan_grid
    near = np.arange(grid.z_min, grid.z_mid + 0.0001, grid.step_near)
    mid = np.arange(grid.z_mid + grid.step_mid, grid.z_upper_mid + 0.0001, grid.step_mid)
    far = np.arange(grid.z_upper_mid + grid.step_far, grid.z_max + 0.0001, grid.step_far)
    return [float(x) for x in np.concatenate([near, mid, far])]


def rigid_scan_geometry(minimum_atoms: Atoms, anchor_global_index: int, target_z: float, metal_symbols: set[str]) -> Atoms:
    atoms = minimum_atoms.copy()
    positions = atoms.positions.copy()
    slab_mask = np.array([sym in metal_symbols for sym in atoms.get_chemical_symbols()])
    top_z = float(np.max(positions[slab_mask, 2]))
    current_z = float(positions[anchor_global_index, 2] - top_z)
    shift = target_z - current_z
    positions[~slab_mask, 2] += shift
    atoms.positions = positions
    return atoms


def write_phase_pipeline(root_dir: Path, stage_sequence: tuple[str, ...]) -> None:
    script = relax_pipeline_script(stage_sequence)
    path = root_dir / "run_pipeline_local.sh"
    path.write_text(script)
    path.chmod(0o755)


def slab_stage_sequence(config: CampaignConfig) -> tuple[str, ...]:
    dipole_cycles = tuple(
        f"relax_stage2_cycle{cycle}_dipole"
        for cycle in range(1, config.slab_relax_stage2_cycles + 1)
    )
    return (
        config.slab_relax_stage1_name,
        *dipole_cycles,
        config.slab_final_stage_name,
        config.slab_workfunc_stage_name,
        config.slab_bader_stage_name,
    )


def restart_source_for_stage(stage_sequence: tuple[str, ...], stage_name: str, config: CampaignConfig) -> str:
    if stage_name in {config.slab_workfunc_stage_name, config.slab_bader_stage_name}:
        return config.slab_final_stage_name
    return stage_sequence[stage_sequence.index(stage_name) - 1]


def candidate_from_label(config: CampaignConfig, label: str) -> SlabCandidate:
    for candidate in config.candidate_slabs:
        if candidate.label == label:
            return candidate
    raise ValueError(f"Unknown slab candidate: {label}")


def create_local_pilot_campaign(config: CampaignConfig, campaign_root: Path) -> dict:
    ensure_dir(campaign_root)
    ensure_dir(campaign_root / "launchers")
    manifest_jobs: list[dict] = []
    save_json(config.to_json(), campaign_root / "campaign_config.json")
    adsorbate_specs = default_adsorbate_specs()

    for phase in config.campaign_phases:
        ensure_dir(campaign_root / phase)

    for metal in config.metals:
        bulk_atoms = build_fcc_bulk(metal, config.bulk_seed_lattice[metal])
        bulk_root = campaign_root / "00_bulk" / metal
        bulk_metadata = {
            "phase": "00_bulk",
            "metal": metal,
            "job_kind": "bulk_reference",
            "protocol": config.primary_protocol.name,
            "local_execution_tier": "local_pilot",
        }
        bulk_stage_specs = (
            (
                "relax_stage1",
                bulk_relax_incar(config.primary_protocol, f"{metal}_bulk_stage1", "relax_stage1", config.hpc.ncore_bulk),
                kpoints_grid(default_bulk_kmesh(), gamma=True),
            ),
            (
                "relax_stage2",
                bulk_relax_incar(config.primary_protocol, f"{metal}_bulk_stage2", "relax_stage2", config.hpc.ncore_bulk),
                kpoints_grid(default_bulk_kmesh(), gamma=True),
            ),
            (
                "final_static",
                final_static_incar(
                    config.primary_protocol,
                    f"{metal}_bulk_final_static",
                    metallic=True,
                    write_volumetrics=False,
                    restart_mode="cold",
                    apply_dipole=False,
                    surface_screened=False,
                    ncore=config.hpc.ncore_bulk,
                    include_dispersion=False,
                ),
                kpoints_grid(bulk_final_kmesh(), gamma=True),
            ),
        )
        for stage_name, incar_text, kpoints_text in bulk_stage_specs:
            stage_dir = bulk_root / stage_name
            metadata = dict(bulk_metadata)
            metadata["stage"] = stage_name
            if stage_name != config.bulk_relax_stages[0]:
                previous_stage = config.bulk_relax_stages[config.bulk_relax_stages.index(stage_name) - 1]
                metadata["source_stage"] = previous_stage
                metadata["restart_files_from"] = str((bulk_root / previous_stage).relative_to(campaign_root))
            write_job_dir(
                stage_dir,
                bulk_atoms,
                incar_text,
                kpoints_text,
                config,
                metadata,
                f"{metal}_{stage_name}",
                config.hpc.ncore_bulk,
            )
            manifest_jobs.append({"path": str(stage_dir.relative_to(campaign_root)), **metadata})
        write_phase_pipeline(bulk_root, config.bulk_relax_stages)

    ag_lattice, ag_lattice_source = resolve_bulk_lattice_constant(config, campaign_root, "Ag", allow_reference_fallback=True)
    for candidate in config.candidate_slabs:
        slab, fixed_indices = build_slab("Ag", ag_lattice, candidate)
        slab_root = campaign_root / "01_universal_slab_selection" / "Ag" / candidate.label
        clean_root = slab_root / "clean_slab"
        stage_meta = {
            "phase": "01_universal_slab_selection",
            "metal": "Ag",
            "candidate_slab": candidate.label,
            "fixed_indices": fixed_indices,
            "bulk_lattice_constant": ag_lattice,
            "bulk_lattice_source": ag_lattice_source,
            "protocol": config.primary_protocol.name,
            "job_kind": "clean_slab",
        }
        stage_sequence = slab_stage_sequence(config)
        for stage_name in stage_sequence:
            stage_dir = clean_root / stage_name
            if stage_name == config.slab_final_stage_name:
                incar = final_static_incar(
                    config.primary_protocol,
                    f"Ag_{candidate.label}_slab_final",
                    metallic=True,
                    write_volumetrics=True,
                    restart_mode="cold",
                    apply_dipole=True,
                    surface_screened=True,
                    ncore=config.hpc.ncore_slab,
                )
                kpoints_text = kpoints_grid(slab_dense_kmesh(candidate))
            elif stage_name == config.slab_workfunc_stage_name:
                incar = workfunction_incar(config.primary_protocol, f"Ag_{candidate.label}_slab_workfunc", config.hpc.ncore_slab)
                kpoints_text = kpoints_grid(slab_dense_kmesh(candidate))
            elif stage_name == config.slab_bader_stage_name:
                incar = bader_incar(config.primary_protocol, f"Ag_{candidate.label}_slab_bader", config.hpc.ncore_slab)
                kpoints_text = kpoints_grid(slab_dense_kmesh(candidate))
            else:
                incar = slab_relax_incar(
                    config.primary_protocol,
                    f"Ag_{candidate.label}_clean_{stage_name}",
                    stage_name,
                    config.hpc.ncore_slab,
                    clean_slab=True,
                )
                kpoints_text = kpoints_grid(slab_kmesh(candidate))
            metadata = dict(stage_meta)
            metadata["stage"] = stage_name
            if stage_name != stage_sequence[0]:
                previous_stage = restart_source_for_stage(stage_sequence, stage_name, config)
                metadata["source_stage"] = previous_stage
                metadata["restart_files_from"] = str((clean_root / previous_stage).relative_to(campaign_root))
            write_job_dir(
                stage_dir,
                slab,
                incar,
                kpoints_text,
                config,
                metadata,
                f"Ag_{candidate.label}_{stage_name}",
                config.hpc.ncore_slab,
                use_custodian_hpc=True,
            )
            manifest_jobs.append({"path": str(stage_dir.relative_to(campaign_root)), **metadata})
        write_phase_pipeline(clean_root, stage_sequence)

        for ads_name in config.representative_adsorbates:
            spec = adsorbate_specs[ads_name]
            seed_root = slab_root / "representative_adsorbates" / ads_name
            for placement in enumerate_seed_placements(spec):
                for site_label in config.site_labels:
                    combined, placement_meta = place_adsorbate_on_slab(
                        slab,
                        spec,
                        placement,
                        site_label,
                        next(seed.initial_height for seed in spec.orientation_seeds if seed.label == placement.label),
                    )
                    initial_height = next(seed.initial_height for seed in spec.orientation_seeds if seed.label == placement.label)
                    az_label = f"{int(round(placement.azimuth_deg)):03d}"
                    variant_label = f"{placement.label}_{site_label}_{az_label}_{int(round(initial_height * 100)):03d}"
                    variant_root = seed_root / variant_label
                    base_meta = {
                        "phase": "01_universal_slab_selection",
                        "metal": "Ag",
                        "candidate_slab": candidate.label,
                        "adsorbate": ads_name,
                        "site_label": site_label,
                        "placement_label": placement.label,
                        "azimuth_deg": placement.azimuth_deg,
                        "bulk_lattice_constant": ag_lattice,
                        "bulk_lattice_source": ag_lattice_source,
                        "protocol": config.primary_protocol.name,
                        "job_kind": "representative_ads_seed",
                    }
                    base_meta.update(placement_meta)
                    for stage_name in stage_sequence:
                        stage_dir = variant_root / stage_name
                        if stage_name == config.slab_final_stage_name:
                            incar = final_static_incar(
                                config.primary_protocol,
                                f"Ag_{candidate.label}_{ads_name}_{variant_label}_final",
                                metallic=True,
                                write_volumetrics=True,
                                restart_mode="cold",
                                apply_dipole=True,
                                surface_screened=True,
                                ncore=config.hpc.ncore_slab,
                            )
                            kpoints_text = kpoints_grid(slab_dense_kmesh(candidate))
                        elif stage_name == config.slab_workfunc_stage_name:
                            incar = workfunction_incar(
                                config.primary_protocol,
                                f"Ag_{candidate.label}_{ads_name}_{variant_label}_workfunc",
                                config.hpc.ncore_slab,
                            )
                            kpoints_text = kpoints_grid(slab_dense_kmesh(candidate))
                        elif stage_name == config.slab_bader_stage_name:
                            incar = bader_incar(
                                config.primary_protocol,
                                f"Ag_{candidate.label}_{ads_name}_{variant_label}_bader",
                                config.hpc.ncore_slab,
                            )
                            kpoints_text = kpoints_grid(slab_dense_kmesh(candidate))
                        else:
                            incar = slab_relax_incar(
                                config.primary_protocol,
                                f"Ag_{candidate.label}_{ads_name}_{variant_label}_{stage_name}",
                                stage_name,
                                config.hpc.ncore_slab,
                            )
                            kpoints_text = kpoints_grid(slab_kmesh(candidate))
                        metadata = dict(base_meta)
                        metadata["stage"] = stage_name
                        if stage_name != stage_sequence[0]:
                            previous_stage = restart_source_for_stage(stage_sequence, stage_name, config)
                            metadata["source_stage"] = previous_stage
                            metadata["restart_files_from"] = str((variant_root / previous_stage).relative_to(campaign_root))
                        write_job_dir(
                            stage_dir,
                            combined,
                            incar,
                            kpoints_text,
                            config,
                            metadata,
                            f"Ag_{ads_name}_{stage_name}",
                            config.hpc.ncore_slab,
                            use_custodian_hpc=True,
                        )
                        manifest_jobs.append({"path": str(stage_dir.relative_to(campaign_root)), **metadata})
                    write_phase_pipeline(variant_root, stage_sequence)

    for ads_name in config.production_adsorbates:
        ads = place_in_box(build_adsorbate(ads_name))
        gas_root = campaign_root / "03_gas_refs" / ads_name
        is_atomic_h = ads_name == "H"
        gas_ncore = 1 if is_atomic_h else config.hpc.ncore_molecule
        gas_ncpus = 4 if is_atomic_h else config.hpc.ncpus_molecule
        gas_mem = "8gb" if is_atomic_h else config.hpc.mem_molecule
        gas_spin = is_atomic_h
        gas_magmom = "1" if is_atomic_h else None
        gas_nupdown = 1 if is_atomic_h else None
        gas_metadata = {
            "phase": "03_gas_refs",
            "adsorbate": ads_name,
            "job_kind": "gas_reference",
            "protocol": config.primary_protocol.name,
        }
        for stage_name in config.gas_relax_stages:
            stage_dir = gas_root / stage_name
            if stage_name == "final_static":
                incar_text = final_static_incar(
                    config.primary_protocol,
                    f"{ads_name}_gas_final_static",
                    metallic=False,
                    write_volumetrics=True,
                    restart_mode="cold",
                    apply_dipole=False,
                    surface_screened=False,
                    ncore=gas_ncore,
                    spin_polarized=gas_spin,
                    magmom=gas_magmom,
                    nupdown=gas_nupdown,
                )
            else:
                incar_text = gas_phase_incar(
                    config.primary_protocol,
                    f"{ads_name}_gas_relax",
                    stage_name,
                    gas_ncore,
                    spin_polarized=gas_spin,
                    magmom=gas_magmom,
                    nupdown=gas_nupdown,
                )
            metadata = dict(gas_metadata)
            metadata["stage"] = stage_name
            if stage_name != config.gas_relax_stages[0]:
                previous_stage = config.gas_relax_stages[config.gas_relax_stages.index(stage_name) - 1]
                metadata["source_stage"] = previous_stage
                metadata["restart_files_from"] = str((gas_root / previous_stage).relative_to(campaign_root))
            write_job_dir(
                stage_dir,
                ads,
                incar_text,
                kpoints_grid(gas_kmesh(), gamma=True),
                config,
                metadata,
                f"{ads_name}_{stage_name}",
                gas_ncore,
                hpc_ncpus=gas_ncpus,
                hpc_mem=gas_mem,
                hpc_scratch=config.hpc.scratch_local_molecule,
                hpc_binary_default=config.hpc.vasp_binary_default,
                hpc_binary_env=config.hpc.vasp_binary_env,
                use_custodian_hpc=True,
            )
            manifest_jobs.append({"path": str(stage_dir.relative_to(campaign_root)), **metadata})
        write_phase_pipeline(gas_root, config.gas_relax_stages)

    manifest = {
        "campaign_root": str(campaign_root),
        "workflow": "coinage_gridff_dft_local_first",
        "jobs": manifest_jobs,
        "notes": [
            "The reference slab is fixed to 3x3x4 for the current campaign.",
            "00_bulk and 01_universal_slab_selection are intended for immediate local validation or ORR-reference comparison.",
            "Bulk uses relax_stage1 -> relax_stage2 -> final_static with CONTCAR handoff and no surface-screened dispersion tags.",
            "Slabs and representative adsorbates use nodipole relax, repeated dipole-enabled cycles, then final_static, workfunction, and bader.",
            "Representative and production adsorption seeds are generated exhaustively across ontop/bridge/fcc/hcp.",
            "Rigid and relaxed scan jobs are generated later from converged minima using setup_scans_from_minima.py.",
        ],
    }
    save_json(manifest, campaign_root / "campaign_manifest.json")
    return manifest


def create_clean_slab_final_phase(
    config: CampaignConfig,
    campaign_root: Path,
    selected_slab_label: str,
) -> dict:
    candidate = candidate_from_label(config, selected_slab_label)
    manifest_jobs: list[dict] = []
    stage_sequence = slab_stage_sequence(config)
    for metal in config.metals:
        lattice_constant, lattice_source = resolve_bulk_lattice_constant(config, campaign_root, metal)
        slab, fixed_indices = build_slab(metal, lattice_constant, candidate)
        clean_root = campaign_root / "02_clean_slab_final" / metal / "clean_slab"
        base_meta = {
            "phase": "02_clean_slab_final",
            "metal": metal,
            "candidate_slab": candidate.label,
            "bulk_lattice_constant": lattice_constant,
            "bulk_lattice_source": lattice_source,
            "fixed_indices": fixed_indices,
            "protocol": config.primary_protocol.name,
            "job_kind": "clean_slab_final",
            "selected_slab": selected_slab_label,
        }
        for stage_name in stage_sequence:
            stage_dir = clean_root / stage_name
            if stage_name == config.slab_final_stage_name:
                incar = final_static_incar(
                    config.primary_protocol,
                    f"{metal}_{candidate.label}_clean_slab_final",
                    metallic=True,
                    write_volumetrics=True,
                    restart_mode="cold",
                    apply_dipole=True,
                    surface_screened=True,
                    ncore=config.hpc.ncore_slab,
                )
                kpoints_text = kpoints_grid(slab_dense_kmesh(candidate))
            elif stage_name == config.slab_workfunc_stage_name:
                incar = workfunction_incar(
                    config.primary_protocol,
                    f"{metal}_{candidate.label}_clean_slab_workfunc",
                    config.hpc.ncore_slab,
                )
                kpoints_text = kpoints_grid(slab_dense_kmesh(candidate))
            elif stage_name == config.slab_bader_stage_name:
                incar = bader_incar(
                    config.primary_protocol,
                    f"{metal}_{candidate.label}_clean_slab_bader",
                    config.hpc.ncore_slab,
                )
                kpoints_text = kpoints_grid(slab_dense_kmesh(candidate))
            else:
                incar = slab_relax_incar(
                    config.primary_protocol,
                    f"{metal}_{candidate.label}_clean_{stage_name}",
                    stage_name,
                    config.hpc.ncore_slab,
                    clean_slab=True,
                )
                kpoints_text = kpoints_grid(slab_kmesh(candidate))
            metadata = dict(base_meta)
            metadata["stage"] = stage_name
            if stage_name != stage_sequence[0]:
                previous_stage = restart_source_for_stage(stage_sequence, stage_name, config)
                metadata["source_stage"] = previous_stage
                metadata["restart_files_from"] = str((clean_root / previous_stage).relative_to(campaign_root))
            write_job_dir(
                stage_dir,
                slab,
                incar,
                kpoints_text,
                config,
                metadata,
                f"{metal}_{candidate.label}_{stage_name}",
                config.hpc.ncore_slab,
            )
            manifest_jobs.append({"path": str(stage_dir.relative_to(campaign_root)), **metadata})
        write_phase_pipeline(clean_root, stage_sequence)
    manifest = {
        "selected_slab": selected_slab_label,
        "jobs": manifest_jobs,
    }
    save_json(manifest, campaign_root / "02_clean_slab_final" / "phase_manifest.json")
    return manifest


def production_relax_stage_sequence(config: CampaignConfig) -> tuple[str, ...]:
    dipole_cycles = tuple(
        f"relax_stage2_cycle{cycle}_dipole"
        for cycle in range(1, config.slab_relax_stage2_cycles + 1)
    )
    return (
        config.slab_relax_stage1_name,
        *dipole_cycles,
        config.slab_final_stage_name,
    )


def create_production_adsorption_phases(
    config: CampaignConfig,
    campaign_root: Path,
    selected_slab_label: str,
) -> dict:
    candidate = candidate_from_label(config, selected_slab_label)
    adsorbate_specs = default_adsorbate_specs()
    seed_entries: list[dict] = []
    relax_jobs: list[dict] = []
    stage_sequence = production_relax_stage_sequence(config)

    for metal in config.metals:
        lattice_constant, lattice_source = resolve_bulk_lattice_constant(config, campaign_root, metal)
        slab, fixed_indices = build_slab(metal, lattice_constant, candidate)
        for ads_name in config.production_adsorbates:
            spec = adsorbate_specs[ads_name]
            seed_root = campaign_root / "04_ads_seed_library" / metal / ads_name
            relax_root = campaign_root / "05_ads_relax" / metal / ads_name
            for placement in enumerate_seed_placements(spec):
                initial_height = next(seed.initial_height for seed in spec.orientation_seeds if seed.label == placement.label)
                for site_label in config.site_labels:
                    combined, placement_meta = place_adsorbate_on_slab(
                        slab,
                        spec,
                        placement,
                        site_label,
                        initial_height,
                    )
                    az_label = f"{int(round(placement.azimuth_deg)):03d}"
                    variant_label = f"{placement.label}_{site_label}_{az_label}_{int(round(initial_height * 100)):03d}"
                    variant_seed_root = seed_root / variant_label
                    variant_relax_root = relax_root / variant_label
                    ensure_dir(variant_seed_root)
                    ase_write(variant_seed_root / "POSCAR", combined, format="vasp", direct=True, sort=False)
                    seed_meta = {
                        "phase": "04_ads_seed_library",
                        "metal": metal,
                        "adsorbate": ads_name,
                        "candidate_slab": candidate.label,
                        "fixed_indices": fixed_indices,
                        "bulk_lattice_constant": lattice_constant,
                        "bulk_lattice_source": lattice_source,
                        "site_label": site_label,
                        "placement_label": placement.label,
                        "azimuth_deg": placement.azimuth_deg,
                        "initial_height": initial_height,
                        "anchor_index_local": placement.anchor_index,
                        "protocol": config.primary_protocol.name,
                        "job_kind": "adsorption_seed",
                    }
                    seed_meta.update(placement_meta)
                    save_json(seed_meta, variant_seed_root / "seed_manifest.json")
                    seed_entries.append({"path": str(variant_seed_root.relative_to(campaign_root)), **seed_meta})

                    relax_meta = {
                        "phase": "05_ads_relax",
                        "metal": metal,
                        "adsorbate": ads_name,
                        "candidate_slab": candidate.label,
                        "fixed_indices": fixed_indices,
                        "bulk_lattice_constant": lattice_constant,
                        "bulk_lattice_source": lattice_source,
                        "site_label": site_label,
                        "placement_label": placement.label,
                        "azimuth_deg": placement.azimuth_deg,
                        "initial_height": initial_height,
                        "anchor_index_local": placement.anchor_index,
                        "protocol": config.primary_protocol.name,
                        "job_kind": "adsorption_relax",
                        "seed_source": str(variant_seed_root.relative_to(campaign_root)),
                    }
                    relax_meta.update(placement_meta)
                    for stage_name in stage_sequence:
                        stage_dir = variant_relax_root / stage_name
                        if stage_name == config.slab_final_stage_name:
                            incar = final_static_incar(
                                config.primary_protocol,
                                f"{metal}_{ads_name}_{variant_label}_final",
                                metallic=True,
                                write_volumetrics=False,
                                restart_mode="cold",
                                apply_dipole=True,
                                surface_screened=True,
                                ncore=config.hpc.ncore_slab,
                            )
                            kpoints_text = kpoints_grid(slab_dense_kmesh(candidate))
                        else:
                            incar = slab_relax_incar(
                                config.primary_protocol,
                                f"{metal}_{ads_name}_{variant_label}_{stage_name}",
                                stage_name,
                                config.hpc.ncore_slab,
                            )
                            kpoints_text = kpoints_grid(slab_kmesh(candidate))
                        metadata = dict(relax_meta)
                        metadata["stage"] = stage_name
                        if stage_name != stage_sequence[0]:
                            previous_stage = restart_source_for_stage(stage_sequence, stage_name, config)
                            metadata["source_stage"] = previous_stage
                            metadata["restart_files_from"] = str((variant_relax_root / previous_stage).relative_to(campaign_root))
                        write_job_dir(
                            stage_dir,
                            combined,
                            incar,
                            kpoints_text,
                            config,
                            metadata,
                            f"{metal}_{ads_name}_{stage_name}",
                            config.hpc.ncore_slab,
                        )
                        relax_jobs.append({"path": str(stage_dir.relative_to(campaign_root)), **metadata})
                    write_phase_pipeline(variant_relax_root, stage_sequence)

    seed_manifest = {
        "selected_slab": selected_slab_label,
        "entries": seed_entries,
    }
    relax_manifest = {
        "selected_slab": selected_slab_label,
        "jobs": relax_jobs,
    }
    save_json(seed_manifest, campaign_root / "04_ads_seed_library" / "phase_manifest.json")
    save_json(relax_manifest, campaign_root / "05_ads_relax" / "phase_manifest.json")
    return {
        "selected_slab": selected_slab_label,
        "seed_entries": seed_entries,
        "relax_jobs": relax_jobs,
    }


def bottom_fixed_metal_indices(full_atoms: Atoms, n_fixed_layers: int, metal_symbols: set[str]) -> list[int]:
    metal_indices = [i for i, sym in enumerate(full_atoms.get_chemical_symbols()) if sym in metal_symbols]
    z_positions = full_atoms.positions[:, 2]
    unique_z = sorted(set(np.round(z_positions[metal_indices], 4)))
    if n_fixed_layers <= 0:
        return []
    z_cutoff = unique_z[n_fixed_layers - 1] + 1.0e-3
    return [i for i in metal_indices if z_positions[i] <= z_cutoff]


def build_scan_constraints(
    full_atoms: Atoms,
    anchor_index_global: int,
    family: str,
    metal_symbols: set[str],
    n_fixed_layers: int,
) -> Atoms:
    atoms = full_atoms.copy()
    if family == "rigid":
        fixed_indices = [i for i, sym in enumerate(atoms.get_chemical_symbols()) if sym in metal_symbols]
        fixed_indices.extend([i for i, sym in enumerate(atoms.get_chemical_symbols()) if sym not in metal_symbols])
    elif family == "relaxed":
        fixed_indices = [i for i, sym in enumerate(atoms.get_chemical_symbols()) if sym in metal_symbols]
        fixed_indices.append(anchor_index_global)
    elif family == "relaxed_slab":
        fixed_indices = bottom_fixed_metal_indices(atoms, n_fixed_layers, metal_symbols)
        fixed_indices.append(anchor_index_global)
    else:
        raise ValueError(f"Unsupported scan family: {family}")
    atoms.set_constraint(FixAtoms(indices=sorted(set(fixed_indices))))
    return atoms


def detect_adsorbate_metadata(job_manifest_path: Path) -> dict:
    return json.loads(job_manifest_path.read_text())


def create_scan_jobs_from_minimum(
    minimum_dir: Path,
    out_root: Path,
    config: CampaignConfig,
    family: str,
) -> dict:
    job_manifest = detect_adsorbate_metadata(minimum_dir / "job_manifest.json")
    structure_path = minimum_dir / "CONTCAR"
    if not structure_path.exists():
        raise FileNotFoundError(f"Missing converged CONTCAR for scan source: {minimum_dir}")
    minimum_atoms = ase_read(structure_path)
    adsorbate_name = job_manifest["adsorbate"]
    ads_specs = default_adsorbate_specs()
    spec = ads_specs[adsorbate_name]
    slab_atom_count = sum(sym in {"Ag", "Cu", "Au"} for sym in minimum_atoms.get_chemical_symbols())
    anchor_index_global = slab_atom_count + job_manifest["anchor_index_local"]
    z_values = rigid_scan_values(config)
    family_root = (
        out_root
        / job_manifest["metal"]
        / job_manifest["candidate_slab"]
        / adsorbate_name
        / job_manifest["site_label"]
        / minimum_dir.parent.name
        / minimum_dir.name
        / family
    )
    ensure_dir(family_root)
    slab_candidate = candidate_from_label(config, job_manifest["candidate_slab"])
    jobs = []
    for z_value in z_values:
        label = f"z_{z_value:05.2f}".replace(".", "p")
        point_root = family_root / label
        geometry = rigid_scan_geometry(minimum_atoms, anchor_index_global, z_value, {"Ag", "Cu", "Au"})
        constrained = build_scan_constraints(
            geometry,
            anchor_index_global,
            family=family,
            metal_symbols={"Ag", "Cu", "Au"},
            n_fixed_layers=slab_candidate.n_fixed_layers,
        )
        point_meta = {
            "phase": (
                "06_scan_rigid"
                if family == "rigid"
                else ("07_scan_relaxed" if family == "relaxed" else config.scan_relaxed_slab_phase_name)
            ),
            "source_minimum": str(minimum_dir),
            "adsorbate": adsorbate_name,
            "site_label": job_manifest["site_label"],
            "placement_label": job_manifest["placement_label"],
            "scan_family": family,
            "z_value": z_value,
            "metal": job_manifest["metal"],
            "anchor_index_global": anchor_index_global,
            "anchor_index_local": job_manifest["anchor_index_local"],
            "protocol": config.primary_protocol.name,
            "candidate_slab": job_manifest["candidate_slab"],
        }
        if family == "rigid":
            static_dir = point_root / "full_system" / "static"
            write_job_dir(
                static_dir,
                constrained,
                rigid_scan_incar(config.primary_protocol, f"{job_manifest['metal']}_{adsorbate_name}_{label}_rigid", config.hpc.ncore_slab),
                kpoints_grid(slab_kmesh(slab_candidate)),
                config,
                dict(point_meta, stage="static", job_kind="rigid_scan"),
                f"{job_manifest['metal']}_{adsorbate_name}_rigid",
                config.hpc.ncore_slab,
            )
            jobs.append({"path": str(static_dir.relative_to(out_root)), **point_meta, "stage": "static"})
        else:
            relax_dir = point_root / "full_system" / "relax"
            final_dir = point_root / "full_system" / "final_static"
            write_job_dir(
                relax_dir,
                constrained,
                relaxed_scan_incar(
                    config.primary_protocol,
                    f"{job_manifest['metal']}_{adsorbate_name}_{label}_relax",
                    config.hpc.ncore_slab,
                ),
                kpoints_grid(slab_kmesh(slab_candidate)),
                config,
                dict(
                    point_meta,
                    stage="relax",
                    job_kind="relaxed_scan" if family == "relaxed" else "relaxed_slab_scan",
                ),
                f"{job_manifest['metal']}_{adsorbate_name}_relax",
                config.hpc.ncore_slab,
            )
            write_job_dir(
                final_dir,
                constrained,
                final_static_incar(
                    config.primary_protocol,
                    f"{job_manifest['metal']}_{adsorbate_name}_{label}_final",
                    metallic=True,
                    write_volumetrics=False,
                    restart_mode="restart",
                    apply_dipole=True,
                    surface_screened=True,
                    ncore=config.hpc.ncore_slab,
                ),
                kpoints_grid(slab_kmesh(slab_candidate)),
                config,
                dict(
                    point_meta,
                    stage="final_static",
                    job_kind="relaxed_scan" if family == "relaxed" else "relaxed_slab_scan",
                ),
                f"{job_manifest['metal']}_{adsorbate_name}_final",
                config.hpc.ncore_slab,
            )
            write_phase_pipeline(point_root / "full_system", ("relax", "final_static"))
            jobs.append({"path": str(relax_dir.relative_to(out_root)), **point_meta, "stage": "relax"})
            jobs.append({"path": str(final_dir.relative_to(out_root)), **point_meta, "stage": "final_static"})
    summary = {
        "minimum_dir": str(minimum_dir),
        "scan_family": family,
        "jobs": jobs,
    }
    save_json(summary, family_root / "scan_manifest.json")
    return summary


def create_interaction_reference_jobs(scan_root: Path, out_root: Path, config: CampaignConfig) -> dict:
    ensure_dir(out_root)
    jobs = []
    final_dirs = sorted(scan_root.rglob("full_system/final_static"))
    static_dirs = sorted(scan_root.rglob("full_system/static"))
    source_dirs = static_dirs + final_dirs
    for source_dir in source_dirs:
        structure_path = source_dir / "CONTCAR"
        if not structure_path.exists():
            raise FileNotFoundError(f"Missing converged CONTCAR for reference source: {source_dir}")
        atoms = ase_read(structure_path)
        metal_mask = np.array([sym in {"Ag", "Cu", "Au"} for sym in atoms.get_chemical_symbols()])
        slab_atoms = atoms[metal_mask]
        mol_atoms = atoms[~metal_mask]
        slab_atoms.cell = atoms.cell.copy()
        slab_atoms.pbc = atoms.pbc
        mol_atoms.cell = atoms.cell.copy()
        mol_atoms.pbc = atoms.pbc
        relative = source_dir.relative_to(scan_root)
        slab_dir = out_root / relative / "slab_only"
        mol_dir = out_root / relative / "molecule_only"
        source_meta = json.loads((source_dir / "job_manifest.json").read_text())
        slab_candidate = candidate_from_label(config, source_meta["candidate_slab"])
        write_job_dir(
            slab_dir,
            slab_atoms,
            final_static_incar(
                config.primary_protocol,
                f"{source_meta['metal']}_{source_meta['adsorbate']}_slab_ref",
                metallic=True,
                write_volumetrics=False,
                restart_mode="cold",
                apply_dipole=True,
                surface_screened=True,
                ncore=config.hpc.ncore_slab,
            ),
            kpoints_grid(slab_kmesh(slab_candidate)),
            config,
            dict(source_meta, phase="09_volumetrics", reference_kind="slab_only"),
            f"{source_meta['metal']}_{source_meta['adsorbate']}_slab_ref",
            config.hpc.ncore_slab,
        )
        write_job_dir(
            mol_dir,
            mol_atoms,
            final_static_incar(
                config.primary_protocol,
                f"{source_meta['metal']}_{source_meta['adsorbate']}_mol_ref",
                metallic=False,
                write_volumetrics=False,
                restart_mode="cold",
                apply_dipole=False,
                surface_screened=False,
                ncore=config.hpc.ncore_molecule,
            ),
            kpoints_grid(gas_kmesh(), gamma=True),
            config,
            dict(source_meta, phase="09_volumetrics", reference_kind="molecule_only"),
            f"{source_meta['metal']}_{source_meta['adsorbate']}_mol_ref",
            config.hpc.ncore_molecule,
            hpc_ncpus=config.hpc.ncpus_molecule,
            hpc_mem=config.hpc.mem_molecule,
            hpc_scratch=config.hpc.scratch_local_molecule,
            hpc_binary_default=config.hpc.vasp_binary_default,
            hpc_binary_env=config.hpc.vasp_binary_env,
        )
        jobs.append({"path": str(slab_dir.relative_to(out_root)), "phase": "09_volumetrics", "reference_kind": "slab_only"})
        jobs.append({"path": str(mol_dir.relative_to(out_root)), "phase": "09_volumetrics", "reference_kind": "molecule_only"})
    manifest = {
        "source_root": str(scan_root),
        "reference_root": str(out_root),
        "jobs": jobs,
    }
    save_json(manifest, out_root / "reference_manifest.json")
    save_json({"jobs": jobs}, out_root / "phase_manifest.json")
    return manifest
