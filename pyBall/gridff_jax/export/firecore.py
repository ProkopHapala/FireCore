"""Export hybrid GridFF artifacts to FireCore-compatible files."""

from __future__ import annotations

from dataclasses import asdict, is_dataclass

import numpy as np

from ..utils import ensure_dir, save_json


def _jsonable(value):
    if is_dataclass(value):
        return asdict(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, dict):
        return {key: _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    return value


def export_firecore_artifacts(out_dir, density, model, fit_result, toggles, teacher_backend_id: str):
    out_dir = ensure_dir(out_dir)
    pauli = np.asarray(model.grids["pauli_coeff_zyx"], dtype=float)
    london = np.asarray(model.grids["london_coeff_zyx"], dtype=float)
    coulomb = np.asarray(model.grids["coulomb_coeff_zyx"], dtype=float)
    reactive = np.asarray(model.grids["reactive_coeff_zyxc"], dtype=float)

    plq = np.stack([pauli, london, coulomb], axis=-1).transpose((2, 1, 0, 3))
    reactive_export = reactive.transpose((2, 1, 0, 3))
    np.save(out_dir / "Bspline_PLQd.npy", plq)
    np.save(out_dir / "ReactiveChannels.npy", reactive_export)

    # ---- Dispersion (C6) channel — only present if the model installed one.
    # The fitted student energy is  PLQ + Σ_el c6_el * dispersion_el(r) + offset.
    # FireCore C++ runtime currently consumes the PLQ grid only — so we ALSO
    # dump the dispersion grid, the per-element C6 coefficients, and the
    # energy offset to disk, and warn loudly when the fit includes them, so
    # callers know the exported PLQ grid alone DOES NOT reproduce the
    # fitted student.
    c6_dict = dict(getattr(fit_result.params, "c6_coeff", {}) or {})
    energy_offset = float(getattr(fit_result.params, "energy_offset", 0.0))
    disp_grid = model.grids.get("dispersion_zyxc")
    dispersion_present = (
        disp_grid is not None
        and bool(c6_dict)
        and (any(abs(v) > 1e-10 for v in c6_dict.values()) or abs(energy_offset) > 1e-10)
    )
    if dispersion_present:
        disp_array = np.asarray(disp_grid, dtype=float).transpose((2, 1, 0, 3))
        np.save(out_dir / "Bspline_Dispersion_zyxc.npy", disp_array)
        np.savez(
            out_dir / "DispersionFit.npz",
            c6_keys=np.asarray(sorted(c6_dict)),
            c6_values=np.asarray([c6_dict[k] for k in sorted(c6_dict)], dtype=float),
            energy_offset_eV=np.asarray([energy_offset], dtype=float),
            dispersion_channel_order=np.asarray(list(model.reactive_elements)),
        )
        import warnings as _w
        _w.warn(
            "[firecore export] fitted student includes a dispersion (C6) "
            "channel and/or non-zero energy_offset. Bspline_PLQd.npy alone "
            "does NOT reproduce the fitted student. Re-evaluating with the "
            "full Python pipeline requires Bspline_Dispersion_zyxc.npy + "
            "DispersionFit.npz alongside the PLQ grid.",
            RuntimeWarning, stacklevel=2,
        )
    np.savez_compressed(
        out_dir / "MetalResponse.npz",
        image_plane=np.array([fit_result.params.image_plane], dtype=float),
        image_scale=np.array([fit_result.params.image_scale], dtype=float),
        sample_shift_z=np.array([fit_result.params.sample_shift_z], dtype=float),
        coulomb_shift_z=np.array([fit_result.params.coulomb_shift_z], dtype=float),
        reservoir_chi=np.array([fit_result.params.reservoir_chi], dtype=float),
        reservoir_hardness=np.array([fit_result.params.reservoir_hardness], dtype=float),
        pauli_keys=np.asarray(sorted(fit_result.params.pauli)),
        pauli_values=np.asarray([fit_result.params.pauli[key] for key in sorted(fit_result.params.pauli)], dtype=float),
        london_keys=np.asarray(sorted(fit_result.params.london)),
        london_values=np.asarray([fit_result.params.london[key] for key in sorted(fit_result.params.london)], dtype=float),
        reactive_keys=np.asarray(sorted(fit_result.params.reactive)),
        reactive_values=np.asarray([fit_result.params.reactive[key] for key in sorted(fit_result.params.reactive)], dtype=float),
        chi_keys=np.asarray(sorted(fit_result.params.chi)),
        chi_values=np.asarray([fit_result.params.chi[key] for key in sorted(fit_result.params.chi)], dtype=float),
        hardness_keys=np.asarray(sorted(fit_result.params.hardness)),
        hardness_values=np.asarray([fit_result.params.hardness[key] for key in sorted(fit_result.params.hardness)], dtype=float),
    )

    grid_meta = {
        "cell": density.cell.tolist(),
        "grid_shape_xyz": list(density.grid_shape_xyz),
        "grid_shape_zyx": list(density.grid_shape_zyx),
        "voxel": density.voxel.tolist(),
        "origin": density.origin.tolist(),
        "channel_order": ["Pauli", "London", "Coulomb"],
        "reactive_channel_order": list(model.reactive_elements),
        "reference_plane": model.grids["metadata"]["z_ref"],
        "image_plane": fit_result.params.image_plane,
        "sample_shift_z": fit_result.params.sample_shift_z,
        "coulomb_shift_z": fit_result.params.coulomb_shift_z,
        "enabled_terms": _jsonable(toggles),
        "teacher_backend": teacher_backend_id,
        "interpolation_export": "xyzc_bspline_cubic_coefficients",
        "interpolation_runtime": model.interpolation_kind,
        "builder_metadata": _jsonable(model.grids.get("metadata", {})),
    }
    save_json(grid_meta, out_dir / "GridMeta.json")
    save_json(
        {
            "fit_metrics": _jsonable(fit_result.metrics),
            "optimizer": _jsonable(fit_result.optimizer_result),
            "history": _jsonable(fit_result.history),
        },
        out_dir / "metrics.json",
    )
    # Full HybridParameters checkpoint (complete reproducibility — chi,
    # hardness, image_plane, image_damping, sample_shifts, reservoir,
    # total_charge, etc. are not captured by fit_params.json on its own).
    p = fit_result.params
    np.savez(
        out_dir / "HybridParameters_full.npz",
        pauli_keys=np.asarray(sorted(p.pauli)),
        pauli_values=np.asarray([p.pauli[k] for k in sorted(p.pauli)], dtype=float),
        london_keys=np.asarray(sorted(p.london)),
        london_values=np.asarray([p.london[k] for k in sorted(p.london)], dtype=float),
        reactive_keys=np.asarray(sorted(p.reactive)),
        reactive_values=np.asarray([p.reactive[k] for k in sorted(p.reactive)], dtype=float),
        c6_keys=np.asarray(sorted(p.c6_coeff or {})),
        c6_values=np.asarray([(p.c6_coeff or {})[k] for k in sorted(p.c6_coeff or {})], dtype=float),
        static_charge_keys=np.asarray(sorted(p.static_charge or {})),
        static_charge_values=np.asarray(
            [(p.static_charge or {})[k] for k in sorted(p.static_charge or {})], dtype=float),
        chi_keys=np.asarray(sorted(p.chi or {})),
        chi_values=np.asarray([(p.chi or {})[k] for k in sorted(p.chi or {})], dtype=float),
        hardness_keys=np.asarray(sorted(p.hardness or {})),
        hardness_values=np.asarray([(p.hardness or {})[k] for k in sorted(p.hardness or {})], dtype=float),
        req_radius_offset_keys=np.asarray(sorted(p.req_radius_offset or {})),
        req_radius_offset_values=np.asarray(
            [(p.req_radius_offset or {})[k] for k in sorted(p.req_radius_offset or {})], dtype=float),
        req_energy_scale_keys=np.asarray(sorted(p.req_energy_scale or {})),
        req_energy_scale_values=np.asarray(
            [(p.req_energy_scale or {})[k] for k in sorted(p.req_energy_scale or {})], dtype=float),
        image_plane=np.asarray([p.image_plane], dtype=float),
        image_scale=np.asarray([p.image_scale], dtype=float),
        image_damping=np.asarray([p.image_damping], dtype=float),
        sample_shift_z=np.asarray([p.sample_shift_z], dtype=float),
        coulomb_shift_z=np.asarray([p.coulomb_shift_z], dtype=float),
        reservoir_chi=np.asarray([p.reservoir_chi], dtype=float),
        reservoir_hardness=np.asarray([p.reservoir_hardness], dtype=float),
        total_charge=np.asarray([p.total_charge], dtype=float),
        energy_offset_eV=np.asarray([p.energy_offset], dtype=float),
    )
    return {
        "plq_path": str(out_dir / "Bspline_PLQd.npy"),
        "reactive_path": str(out_dir / "ReactiveChannels.npy"),
        "metal_response_path": str(out_dir / "MetalResponse.npz"),
        "grid_meta_path": str(out_dir / "GridMeta.json"),
        "hybrid_full_path": str(out_dir / "HybridParameters_full.npz"),
        "dispersion_present": bool(dispersion_present),
        "dispersion_path": str(out_dir / "Bspline_Dispersion_zyxc.npy") if dispersion_present else None,
        "dispersion_fit_path": str(out_dir / "DispersionFit.npz") if dispersion_present else None,
    }
