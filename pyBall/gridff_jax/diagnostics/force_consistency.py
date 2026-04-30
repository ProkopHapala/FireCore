"""Force-energy derivative consistency check for any MLIP or model.

Compares analytical forces reported by the model with numerical forces
obtained via central finite differences of the energy.  This is the
"smoking gun" diagnostic for MLIP smoothness: if the energy surface is
not smooth, the analytical forces will deviate from dE/dr.

Usage
-----
    from pyBall.gridff_jax.diagnostics.force_consistency import (
        force_energy_consistency,
    )
    result = force_energy_consistency(
        evaluate_fn, density, poses, fd_step=0.005,
    )
    # result["max_inconsistency"]  -> worst-case |F_analytic - F_numerical|
    # result["rms_inconsistency"]  -> RMS over all poses and components
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable

import numpy as np

from ..interfaces import DensityData, PoseBatch, TeacherResult


@dataclass
class ConsistencyResult:
    """Results from force-energy consistency analysis."""
    n_poses: int = 0
    fd_step: float = 0.005
    max_inconsistency: float = 0.0
    rms_inconsistency: float = 0.0
    mae_inconsistency: float = 0.0
    per_atom_rms: np.ndarray | None = None
    per_component_rms: np.ndarray | None = None
    per_pose_max: np.ndarray | None = None
    worst_pose_index: int = 0
    worst_atom_index: int = 0
    worst_component: int = 0
    analytical_forces: np.ndarray | None = None
    numerical_forces: np.ndarray | None = None
    metadata: dict[str, Any] = field(default_factory=dict)


def _perturb_positions(positions: np.ndarray, pose_idx: int,
                       atom_idx: int, comp: int, delta: float) -> np.ndarray:
    """Create a copy of positions with one coordinate perturbed."""
    perturbed = positions.copy()
    perturbed[pose_idx, atom_idx, comp] += delta
    return perturbed


def force_energy_consistency(
    evaluate_fn: Callable[[DensityData, PoseBatch], TeacherResult],
    density: DensityData,
    poses: PoseBatch,
    fd_step: float = 0.005,
    max_poses: int = 50,
    verbose: bool = True,
) -> ConsistencyResult:
    """Check force-energy consistency via central finite differences.

    Parameters
    ----------
    evaluate_fn : callable
        Function with signature ``(density, poses) -> TeacherResult``.
        This can be a teacher backend's evaluate_batch, or any model.
    density : DensityData
        Substrate data.
    poses : PoseBatch
        Poses to check. If more than ``max_poses``, a uniform subset is used.
    fd_step : float
        Finite difference step in Angstrom.
    max_poses : int
        Maximum number of poses to check (for speed).
    verbose : bool
        Print progress.

    Returns
    -------
    ConsistencyResult
        Detailed consistency statistics.
    """
    n_total = len(poses.positions)
    if n_total > max_poses:
        indices = np.linspace(0, n_total - 1, max_poses, dtype=int)
    else:
        indices = np.arange(n_total)

    subset_positions = poses.positions[indices].copy()
    subset_params = poses.pose_params[indices].copy()
    subset_labels = [poses.site_labels[int(i)] for i in indices]
    n_check = len(indices)

    subset_poses = PoseBatch(
        adsorbate=poses.adsorbate,
        positions=subset_positions,
        pose_params=subset_params,
        site_labels=subset_labels,
        metadata=dict(poses.metadata),
    )

    if verbose:
        print(f"[force_consistency] Evaluating analytical forces for {n_check} poses ...",
              flush=True)
    ref_result = evaluate_fn(density, subset_poses)
    analytical_F = np.asarray(ref_result.forces, dtype=float)

    n_atoms = subset_positions.shape[1]
    numerical_F = np.zeros_like(analytical_F)

    total_evals = n_check * n_atoms * 3 * 2
    eval_count = 0

    for p in range(n_check):
        for a in range(n_atoms):
            for c in range(3):
                pos_plus = subset_positions.copy()
                pos_plus[p, a, c] += fd_step
                poses_plus = PoseBatch(
                    adsorbate=poses.adsorbate,
                    positions=pos_plus,
                    pose_params=subset_params.copy(),
                    site_labels=list(subset_labels),
                    metadata=dict(poses.metadata),
                )
                result_plus = evaluate_fn(density, poses_plus)
                e_plus = float(result_plus.energies[p])

                pos_minus = subset_positions.copy()
                pos_minus[p, a, c] -= fd_step
                poses_minus = PoseBatch(
                    adsorbate=poses.adsorbate,
                    positions=pos_minus,
                    pose_params=subset_params.copy(),
                    site_labels=list(subset_labels),
                    metadata=dict(poses.metadata),
                )
                result_minus = evaluate_fn(density, poses_minus)
                e_minus = float(result_minus.energies[p])

                numerical_F[p, a, c] = -(e_plus - e_minus) / (2.0 * fd_step)

                eval_count += 2
                if verbose and eval_count % (max(total_evals // 10, 1)) == 0:
                    print(f"  FD progress: {eval_count}/{total_evals} "
                          f"({100.0 * eval_count / total_evals:.0f}%)", flush=True)

    diff = analytical_F - numerical_F
    abs_diff = np.abs(diff)

    per_atom_rms = np.sqrt(np.mean(diff ** 2, axis=(0, 2)))
    per_component_rms = np.sqrt(np.mean(diff ** 2, axis=(0, 1)))
    per_pose_max = np.max(abs_diff.reshape(n_check, -1), axis=1)

    worst_flat = np.argmax(abs_diff)
    worst_pose = worst_flat // (n_atoms * 3)
    worst_atom = (worst_flat % (n_atoms * 3)) // 3
    worst_comp = worst_flat % 3

    result = ConsistencyResult(
        n_poses=n_check,
        fd_step=fd_step,
        max_inconsistency=float(abs_diff.max()),
        rms_inconsistency=float(np.sqrt(np.mean(diff ** 2))),
        mae_inconsistency=float(np.mean(abs_diff)),
        per_atom_rms=per_atom_rms,
        per_component_rms=per_component_rms,
        per_pose_max=per_pose_max,
        worst_pose_index=int(worst_pose),
        worst_atom_index=int(worst_atom),
        worst_component=int(worst_comp),
        analytical_forces=analytical_F,
        numerical_forces=numerical_F,
        metadata={
            "fd_step_angstrom": fd_step,
            "n_poses_checked": n_check,
            "n_atoms": n_atoms,
            "total_fd_evaluations": total_evals,
        },
    )

    if verbose:
        print(f"[force_consistency] Max |F_anal - F_num| = {result.max_inconsistency:.6f} eV/A")
        print(f"[force_consistency] RMS |F_anal - F_num| = {result.rms_inconsistency:.6f} eV/A")
        print(f"[force_consistency] MAE |F_anal - F_num| = {result.mae_inconsistency:.6f} eV/A")
        comp_names = ["x", "y", "z"]
        print(f"[force_consistency] Worst case: pose {result.worst_pose_index}, "
              f"atom {result.worst_atom_index}, "
              f"component {comp_names[result.worst_component]}")

    return result
