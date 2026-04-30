"""Multi-dimensional fitting pipeline using 6D sampled data.

Replaces the z-scan-only fitting approach with a pipeline that uses
the full 6D (u,v,z,alpha,beta,gamma) pose space for parameter fitting.

Key improvements over the z-scan-only pipeline:
  - Training data includes lateral variation → parameters learn corrugation
  - Stratified train/val/test splits across all dimensions
  - Optional multi-molecule joint fitting for transferability
  - Per-dimension loss weighting

Usage
-----
    from pyBall.gridff_jax.fit.fit_6d import (
        generate_6d_reference_data,
        split_6d_dataset,
        fit_6d_parameters,
    )

    # Step 1: generate reference data
    dataset = generate_6d_reference_data(density, adsorbate, teacher, config)

    # Step 2: split into train/val/test
    train, val, test = split_6d_dataset(dataset)

    # Step 3: fit parameters
    result = fit_6d_parameters(model, train, val, training_config)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from ..interfaces import (
    AdsorbateDefinition,
    DensityData,
    PoseBatch,
    TeacherBackend,
    TeacherResult,
)
from ..pose_sampling.sampler_6d import Sampler6DConfig, build_6d_pose_batch
from ..utils import ensure_dir, save_json


# ---------------------------------------------------------------------------
#  Dataset container
# ---------------------------------------------------------------------------

@dataclass
class Dataset6D:
    """Container for 6D reference data from a teacher model."""
    adsorbate: AdsorbateDefinition
    poses: PoseBatch
    teacher: TeacherResult
    metadata: dict[str, Any] = field(default_factory=dict)

    @property
    def n_poses(self) -> int:
        return len(self.poses.positions)

    @property
    def energies(self) -> np.ndarray:
        return np.asarray(self.teacher.energies, dtype=float)

    @property
    def forces(self) -> np.ndarray:
        return np.asarray(self.teacher.forces, dtype=float)


@dataclass
class SplitDataset6D:
    """Train / validation / test split of a 6D dataset."""
    train: Dataset6D
    val: Dataset6D
    test: Dataset6D
    split_indices: dict[str, np.ndarray] = field(default_factory=dict)


# ---------------------------------------------------------------------------
#  Reference data generation
# ---------------------------------------------------------------------------

def generate_6d_reference_data(
    density: DensityData,
    adsorbate: AdsorbateDefinition,
    teacher_backend: TeacherBackend,
    sampler_config: Sampler6DConfig | None = None,
    chunk_size: int = 64,
    verbose: bool = True,
) -> Dataset6D:
    """Generate 6D reference data by sampling poses and evaluating with teacher.

    Parameters
    ----------
    density : DensityData
        Substrate information.
    adsorbate : AdsorbateDefinition
        Molecule to sample.
    teacher_backend : TeacherBackend
        Teacher model (MAD-SURF, DFT, or any MLIP).
    sampler_config : Sampler6DConfig, optional
        Sampling configuration.
    chunk_size : int
        Number of poses per teacher evaluation batch.
    verbose : bool
        Print progress.

    Returns
    -------
    Dataset6D
        Complete dataset with poses, energies, and forces.
    """
    if sampler_config is None:
        sampler_config = Sampler6DConfig()

    if verbose:
        print(f"[fit_6d] Generating 6D poses for {adsorbate.name} ...", flush=True)

    poses = build_6d_pose_batch(density, adsorbate, sampler_config)
    n = len(poses.positions)

    if verbose:
        print(f"[fit_6d] {n} poses generated. Evaluating teacher ...", flush=True)

    all_energies = []
    all_forces = []
    for start in range(0, n, chunk_size):
        stop = min(start + chunk_size, n)
        idx = np.arange(start, stop, dtype=int)
        partial_poses = PoseBatch(
            adsorbate=poses.adsorbate,
            positions=np.asarray(poses.positions[idx], dtype=float),
            pose_params=np.asarray(poses.pose_params[idx], dtype=float),
            site_labels=[poses.site_labels[int(i)] for i in idx],
            metadata=dict(poses.metadata),
        )
        result = teacher_backend.evaluate_batch(density, partial_poses)
        all_energies.append(np.asarray(result.energies, dtype=float))
        all_forces.append(np.asarray(result.forces, dtype=float))
        if verbose:
            print(f"  Teacher: {stop}/{n} ({100.0 * stop / n:.0f}%)", flush=True)

    teacher_result = TeacherResult(
        energies=np.concatenate(all_energies, axis=0),
        forces=np.concatenate(all_forces, axis=0),
        metadata={"teacher_type": type(teacher_backend).__name__},
    )

    if verbose:
        e = teacher_result.energies
        print(f"[fit_6d] Teacher done. E range: [{e.min():.4f}, {e.max():.4f}] eV")

    return Dataset6D(
        adsorbate=adsorbate,
        poses=poses,
        teacher=teacher_result,
        metadata={
            "sampler_config": sampler_config.__dict__,
            "n_poses": n,
        },
    )


# ---------------------------------------------------------------------------
#  Stratified splitting
# ---------------------------------------------------------------------------

def split_6d_dataset(
    dataset: Dataset6D,
    val_fraction: float = 0.1,
    test_fraction: float = 0.1,
    stratify_by_z: bool = True,
    stratify_by_site: bool = True,
    z_bins: int = 5,
    seed: int = 42,
) -> SplitDataset6D:
    """Split a 6D dataset into train/val/test with stratification.

    Stratification ensures each split contains representative samples from
    every z-window and every surface site type.

    Parameters
    ----------
    dataset : Dataset6D
        Full dataset to split.
    val_fraction : float
        Fraction for validation.
    test_fraction : float
        Fraction for test.
    stratify_by_z : bool
        Stratify by z-height bins.
    stratify_by_site : bool
        Stratify by site label.
    z_bins : int
        Number of z-bins for stratification.
    seed : int
        Random seed.

    Returns
    -------
    SplitDataset6D
    """
    rng = np.random.default_rng(seed)
    n = dataset.n_poses
    indices = np.arange(n)

    z_vals = dataset.poses.pose_params[:, 2]
    labels = dataset.poses.site_labels

    strata_keys = np.zeros(n, dtype=int)

    if stratify_by_z:
        z_quantiles = np.linspace(z_vals.min(), z_vals.max(), z_bins + 1)
        z_bin_idx = np.digitize(z_vals, z_quantiles[1:-1])
        strata_keys = strata_keys * (z_bins + 1) + z_bin_idx

    if stratify_by_site:
        unique_labels = sorted(set(labels))
        label_map = {lab: i for i, lab in enumerate(unique_labels)}
        label_idx = np.array([label_map[lab] for lab in labels])
        n_labels = len(unique_labels)
        strata_keys = strata_keys * (n_labels + 1) + label_idx

    train_idx, val_idx, test_idx = [], [], []
    unique_strata = np.unique(strata_keys)

    for stratum in unique_strata:
        mask = strata_keys == stratum
        stratum_indices = indices[mask]
        rng.shuffle(stratum_indices)

        n_s = len(stratum_indices)
        n_test = max(1, int(n_s * test_fraction))
        n_val = max(1, int(n_s * val_fraction))
        n_train = n_s - n_test - n_val

        if n_train < 1:
            n_train = max(1, n_s - 2)
            n_val = min(n_val, n_s - n_train)
            n_test = n_s - n_train - n_val

        test_idx.extend(stratum_indices[:n_test])
        val_idx.extend(stratum_indices[n_test:n_test + n_val])
        train_idx.extend(stratum_indices[n_test + n_val:])

    train_idx = np.array(train_idx, dtype=int)
    val_idx = np.array(val_idx, dtype=int)
    test_idx = np.array(test_idx, dtype=int)

    def _subset(idx_arr: np.ndarray) -> Dataset6D:
        return Dataset6D(
            adsorbate=dataset.adsorbate,
            poses=PoseBatch(
                adsorbate=dataset.adsorbate,
                positions=dataset.poses.positions[idx_arr],
                pose_params=dataset.poses.pose_params[idx_arr],
                site_labels=[dataset.poses.site_labels[int(i)] for i in idx_arr],
                metadata=dict(dataset.poses.metadata),
            ),
            teacher=TeacherResult(
                energies=dataset.teacher.energies[idx_arr],
                forces=dataset.teacher.forces[idx_arr],
                metadata=dict(dataset.teacher.metadata),
            ),
            metadata={"subset_size": len(idx_arr)},
        )

    return SplitDataset6D(
        train=_subset(train_idx),
        val=_subset(val_idx),
        test=_subset(test_idx),
        split_indices={
            "train": train_idx,
            "val": val_idx,
            "test": test_idx,
        },
    )


# ---------------------------------------------------------------------------
#  Multi-molecule dataset merging
# ---------------------------------------------------------------------------

def merge_datasets(datasets: list[Dataset6D]) -> Dataset6D:
    """Merge multiple 6D datasets (from different molecules) into one.

    This is used for multi-molecule joint fitting to produce transferable
    per-element parameters.
    """
    if len(datasets) == 1:
        return datasets[0]

    all_positions = []
    all_pose_params = []
    all_labels = []
    all_energies = []
    all_forces = []
    mol_names = []

    for ds in datasets:
        all_positions.append(ds.poses.positions)
        all_pose_params.append(ds.poses.pose_params)
        all_labels.extend(
            [f"{ds.adsorbate.name}:{lab}" for lab in ds.poses.site_labels]
        )
        all_energies.append(ds.teacher.energies)
        all_forces.append(ds.teacher.forces)
        mol_names.append(ds.adsorbate.name)

    first = datasets[0]

    max_atoms = max(p.shape[1] for p in all_positions)
    padded_positions = []
    for p in all_positions:
        if p.shape[1] < max_atoms:
            pad_width = ((0, 0), (0, max_atoms - p.shape[1]), (0, 0))
            padded_positions.append(np.pad(p, pad_width, constant_values=0.0))
        else:
            padded_positions.append(p)

    pad_param_dim = max(pp.shape[1] for pp in all_pose_params)
    padded_params = []
    for pp in all_pose_params:
        if pp.shape[1] < pad_param_dim:
            pad_width = ((0, 0), (0, pad_param_dim - pp.shape[1]))
            padded_params.append(np.pad(pp, pad_width, constant_values=0.0))
        else:
            padded_params.append(pp)

    merged_poses = PoseBatch(
        adsorbate=first.adsorbate,
        positions=np.concatenate(padded_positions, axis=0),
        pose_params=np.concatenate(padded_params, axis=0),
        site_labels=all_labels,
        metadata={
            "merged": True,
            "molecule_names": mol_names,
            "molecule_counts": [ds.n_poses for ds in datasets],
            "max_atoms": max_atoms,
        },
    )

    pad_dim = max(f.shape[1] for f in all_forces)
    padded_forces = []
    for f in all_forces:
        if f.shape[1] < pad_dim:
            pad_width = ((0, 0), (0, pad_dim - f.shape[1]), (0, 0))
            padded_forces.append(np.pad(f, pad_width, constant_values=0.0))
        else:
            padded_forces.append(f)

    merged_teacher = TeacherResult(
        energies=np.concatenate(all_energies, axis=0),
        forces=np.concatenate(padded_forces, axis=0),
        metadata={"merged_molecules": mol_names},
    )

    return Dataset6D(
        adsorbate=first.adsorbate,
        poses=merged_poses,
        teacher=merged_teacher,
        metadata={
            "merged": True,
            "molecule_names": mol_names,
            "total_poses": sum(ds.n_poses for ds in datasets),
        },
    )


# ---------------------------------------------------------------------------
#  Save / load 6D datasets
# ---------------------------------------------------------------------------

def save_6d_dataset(dataset: Dataset6D, out_dir: str | Path, prefix: str = ""):
    """Save a 6D dataset to disk as .npz + .json."""
    out_dir = ensure_dir(out_dir)
    name = prefix + dataset.adsorbate.name if prefix else dataset.adsorbate.name

    np.savez_compressed(
        out_dir / f"{name}_6d.npz",
        positions=dataset.poses.positions,
        pose_params=dataset.poses.pose_params,
        site_labels=np.asarray(dataset.poses.site_labels),
        energies=dataset.teacher.energies,
        forces=dataset.teacher.forces,
        adsorbate_symbols=np.asarray(dataset.adsorbate.symbols),
        adsorbate_positions=dataset.adsorbate.positions,
        adsorbate_charges=(
            dataset.adsorbate.charges
            if dataset.adsorbate.charges is not None
            else np.zeros(0, dtype=float)
        ),
        anchor_index=np.int32(dataset.adsorbate.anchor_index),
    )

    meta = {
        "adsorbate_name": dataset.adsorbate.name,
        "n_poses": dataset.n_poses,
        "pose_metadata": dataset.poses.metadata,
        "teacher_metadata": dataset.teacher.metadata,
        "dataset_metadata": dataset.metadata,
    }
    save_json(meta, out_dir / f"{name}_6d.json")


def load_6d_dataset(path: str | Path) -> Dataset6D:
    """Load a 6D dataset from .npz + .json files."""
    path = Path(path)
    data = np.load(path, allow_pickle=True)

    json_path = path.with_name(path.stem + ".json")
    import json
    with json_path.open("r", encoding="utf-8") as fh:
        meta = json.load(fh)

    charges_arr = data["adsorbate_charges"]
    adsorbate = AdsorbateDefinition(
        name=meta["adsorbate_name"],
        symbols=[str(s) for s in data["adsorbate_symbols"].tolist()],
        positions=np.asarray(data["adsorbate_positions"], dtype=float),
        charges=np.asarray(charges_arr, dtype=float) if charges_arr.size > 0 else None,
        anchor_index=int(data["anchor_index"]),
    )

    poses = PoseBatch(
        adsorbate=adsorbate,
        positions=np.asarray(data["positions"], dtype=float),
        pose_params=np.asarray(data["pose_params"], dtype=float),
        site_labels=[str(s) for s in data["site_labels"].tolist()],
        metadata=meta.get("pose_metadata", {}),
    )

    teacher = TeacherResult(
        energies=np.asarray(data["energies"], dtype=float),
        forces=np.asarray(data["forces"], dtype=float),
        metadata=meta.get("teacher_metadata", {}),
    )

    return Dataset6D(
        adsorbate=adsorbate,
        poses=poses,
        teacher=teacher,
        metadata=meta.get("dataset_metadata", {}),
    )
