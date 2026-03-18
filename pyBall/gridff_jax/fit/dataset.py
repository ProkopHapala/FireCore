"""Pose-dataset construction for teacher distillation."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from ..interfaces import AdsorbateDefinition, PoseBatch, TeacherResult
from ..pose_sampling import build_pose_batch, get_adsorbate
from ..utils import ensure_dir, save_json


def _dataset_payload(poses: PoseBatch, teacher: TeacherResult):
    return {
        "adsorbate_name": poses.adsorbate.name,
        "adsorbate_symbols": np.asarray(poses.adsorbate.symbols),
        "adsorbate_positions": np.asarray(poses.adsorbate.positions, dtype=float),
        "adsorbate_charges": np.asarray(poses.adsorbate.charges, dtype=float) if poses.adsorbate.charges is not None else np.zeros(0, dtype=float),
        "anchor_index": np.int32(poses.adsorbate.anchor_index),
        "positions": np.asarray(poses.positions, dtype=float),
        "pose_params": np.asarray(poses.pose_params, dtype=float),
        "site_labels": np.asarray(poses.site_labels),
        "energies": np.asarray(teacher.energies, dtype=float),
        "forces": np.asarray(teacher.forces, dtype=float),
    }


def build_teacher_dataset(density, adsorbate_config, teacher_backend, seed: int = 12345):
    adsorbate = get_adsorbate(
        adsorbate_config.name,
        xyz_path=adsorbate_config.xyz_path,
    )
    poses = build_pose_batch(
        density=density,
        adsorbate=adsorbate,
        samples_per_site=adsorbate_config.samples_per_site,
        systematic_orientations=adsorbate_config.systematic_orientations,
        random_samples=adsorbate_config.random_samples,
        representative_sites_per_label=adsorbate_config.representative_sites_per_label,
        seed=seed,
        z_range=adsorbate_config.z_range,
        z_bias_power=adsorbate_config.z_bias_power,
        jitter_uv=adsorbate_config.jitter_uv,
        tilt_deg=adsorbate_config.tilt_deg,
    )
    teacher = teacher_backend.evaluate_batch(density, poses)
    return poses, teacher


def save_pose_dataset(out_dir, poses: PoseBatch, teacher: TeacherResult, density_metadata=None):
    out_dir = ensure_dir(out_dir)
    payload = _dataset_payload(poses, teacher)
    np.savez_compressed(out_dir / f"{poses.adsorbate.name}.npz", **payload)
    metadata = {
        "adsorbate": poses.adsorbate.name,
        "pose_metadata": poses.metadata,
        "teacher_metadata": teacher.metadata,
        "density_metadata": density_metadata or {},
    }
    save_json(metadata, out_dir / f"{poses.adsorbate.name}.json")


def load_pose_dataset(path: str | Path):
    path = Path(path)
    data = np.load(path, allow_pickle=True)
    with path.with_suffix(".json").open("r", encoding="utf8") as fin:
        metadata = json.load(fin)
    adsorbate = AdsorbateDefinition(
        name=str(data["adsorbate_name"]),
        symbols=[str(v) for v in data["adsorbate_symbols"].tolist()],
        positions=np.asarray(data["adsorbate_positions"], dtype=float),
        charges=np.asarray(data["adsorbate_charges"], dtype=float) if data["adsorbate_charges"].size > 0 else None,
        anchor_index=int(data["anchor_index"]),
    )
    poses = PoseBatch(
        adsorbate=adsorbate,
        positions=np.asarray(data["positions"], dtype=float),
        pose_params=np.asarray(data["pose_params"], dtype=float),
        site_labels=[str(v) for v in data["site_labels"].tolist()],
        metadata=metadata["pose_metadata"],
    )
    teacher = TeacherResult(
        energies=np.asarray(data["energies"], dtype=float),
        forces=np.asarray(data["forces"], dtype=float),
        metadata=metadata["teacher_metadata"],
    )
    return poses, teacher, metadata
