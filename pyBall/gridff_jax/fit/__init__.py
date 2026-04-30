"""Dataset construction and fitting helpers."""

from .dataset import build_teacher_dataset, load_pose_dataset, save_pose_dataset
from .fit_6d import (
    Dataset6D,
    SplitDataset6D,
    generate_6d_reference_data,
    load_6d_dataset,
    merge_datasets,
    save_6d_dataset,
    split_6d_dataset,
)
from .optimize import FitResult, fit_hybrid_parameters, fit_two_stage_c6

__all__ = [
    "build_teacher_dataset",
    "load_pose_dataset",
    "save_pose_dataset",
    "FitResult",
    "fit_hybrid_parameters",
    "fit_two_stage_c6",
    "Dataset6D",
    "SplitDataset6D",
    "generate_6d_reference_data",
    "load_6d_dataset",
    "merge_datasets",
    "save_6d_dataset",
    "split_6d_dataset",
]
