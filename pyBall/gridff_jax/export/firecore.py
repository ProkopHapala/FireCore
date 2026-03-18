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
    return {
        "plq_path": str(out_dir / "Bspline_PLQd.npy"),
        "reactive_path": str(out_dir / "ReactiveChannels.npy"),
        "metal_response_path": str(out_dir / "MetalResponse.npz"),
        "grid_meta_path": str(out_dir / "GridMeta.json"),
    }
