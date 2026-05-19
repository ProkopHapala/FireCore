#!/usr/bin/env python3
"""Minimal unit test for PrecomputedTeacherBackend.

Generates a 5-pose npz with hand-set energies/forces, instantiates the
backend, calls evaluate_batch on a PoseBatch whose pose_params match,
asserts the returned TeacherResult is bit-identical to the npz.

Also covers the missing-pose and adsorbate-mismatch error paths.

Run:
    JAX_PLATFORMS=cpu python3 tests/tGridFFJax/test_precomputed_teacher.py
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.append(str(ROOT))

from pyBall.gridff_jax.config import TeacherBackendConfig
from pyBall.gridff_jax.interfaces import AdsorbateDefinition, PoseBatch
from pyBall.gridff_jax.teacher_backends import make_teacher_backend


def make_test_npz(path: Path):
    rng = np.random.default_rng(42)
    n_pose, n_atom = 5, 3
    pose_params = np.array([
        [0.0, 0.0, 2.5, 1.0, 0.0, 0.0, 0.0],
        [0.5, 0.0, 3.0, 1.0, 0.0, 0.0, 0.0],
        [0.5, 0.5, 3.5, 0.92388, 0.0, 0.38268, 0.0],
        [0.0, 0.5, 4.0, 1.0, 0.0, 0.0, 0.0],
        [0.25, 0.25, 2.75, 1.0, 0.0, 0.0, 0.0],
    ], dtype=float)
    positions = rng.standard_normal((n_pose, n_atom, 3)).astype(float)
    energies = np.array([-0.3, -0.5, -0.4, -0.1, -0.6], dtype=float)
    forces = rng.standard_normal((n_pose, n_atom, 3)).astype(float)
    adsorbate_symbols = np.array(["O", "H", "H"], dtype=object)
    np.savez(path,
             positions=positions, pose_params=pose_params,
             energies=energies, forces=forces,
             adsorbate_symbols=adsorbate_symbols,
             metadata=np.array({"source": "test"}, dtype=object))
    return pose_params, energies, forces


def make_pose_batch(pose_params: np.ndarray, symbols: list[str]) -> PoseBatch:
    n_atom = len(symbols)
    n_pose = pose_params.shape[0]
    adsorbate = AdsorbateDefinition(
        name="water",
        symbols=symbols,
        positions=np.zeros((n_atom, 3)),
        anchor_index=0,
    )
    return PoseBatch(
        adsorbate=adsorbate,
        positions=np.zeros((n_pose, n_atom, 3)),
        pose_params=pose_params,
        site_labels=["test"] * n_pose,
        metadata={},
    )


def test_roundtrip():
    with tempfile.TemporaryDirectory() as td:
        npz_path = Path(td) / "teacher.npz"
        pose_params, energies, forces = make_test_npz(npz_path)

        cfg = TeacherBackendConfig(
            kind="precomputed",
            metadata={"npz_path": str(npz_path)},
        )
        backend = make_teacher_backend(cfg)

        poses = make_pose_batch(pose_params, ["O", "H", "H"])
        result = backend.evaluate_batch(density=None, poses=poses)

    np.testing.assert_array_equal(result.energies, energies)
    np.testing.assert_array_equal(result.forces, forces)
    assert result.metadata["backend"] == "precomputed"
    assert result.metadata["n_cached"] == 5
    print("test_roundtrip: PASSED")


def test_missing_pose():
    with tempfile.TemporaryDirectory() as td:
        npz_path = Path(td) / "teacher.npz"
        pose_params, _, _ = make_test_npz(npz_path)

        cfg = TeacherBackendConfig(
            kind="precomputed",
            metadata={"npz_path": str(npz_path)},
        )
        backend = make_teacher_backend(cfg)

        # Query a pose that's NOT in the cache
        bad_params = pose_params.copy()
        bad_params[0, 2] = 999.0
        poses = make_pose_batch(bad_params, ["O", "H", "H"])
        try:
            backend.evaluate_batch(density=None, poses=poses)
        except LookupError as e:
            assert "not found" in str(e)
            print("test_missing_pose: PASSED (raised LookupError as expected)")
            return
    raise AssertionError("expected LookupError for missing pose")


def test_adsorbate_mismatch():
    with tempfile.TemporaryDirectory() as td:
        npz_path = Path(td) / "teacher.npz"
        pose_params, _, _ = make_test_npz(npz_path)

        cfg = TeacherBackendConfig(
            kind="precomputed",
            metadata={"npz_path": str(npz_path)},
        )
        backend = make_teacher_backend(cfg)

        poses = make_pose_batch(pose_params, ["N", "H", "H"])  # wrong symbol
        try:
            backend.evaluate_batch(density=None, poses=poses)
        except ValueError as e:
            assert "adsorbate mismatch" in str(e)
            print("test_adsorbate_mismatch: PASSED (raised ValueError as expected)")
            return
    raise AssertionError("expected ValueError for adsorbate mismatch")


def test_energy_only():
    """NPZ with no 'forces' array → backend fills zeros, flags forces_available=False."""
    with tempfile.TemporaryDirectory() as td:
        npz_path = Path(td) / "energy_only.npz"
        pose_params = np.array([
            [0.0, 0.0, 2.5, 1.0, 0.0, 0.0, 0.0],
            [0.5, 0.5, 3.0, 1.0, 0.0, 0.0, 0.0],
        ], dtype=float)
        energies = np.array([-0.42, -0.17], dtype=float)
        np.savez(npz_path,
                 pose_params=pose_params, energies=energies,
                 adsorbate_symbols=np.array(["O", "H", "H"], dtype=object))

        cfg = TeacherBackendConfig(
            kind="precomputed",
            metadata={"npz_path": str(npz_path)},
        )
        backend = make_teacher_backend(cfg)

        poses = make_pose_batch(pose_params, ["O", "H", "H"])
        result = backend.evaluate_batch(density=None, poses=poses)

    np.testing.assert_array_equal(result.energies, energies)
    # Forces should be zeros with shape matching the pose batch's adsorbate
    assert result.forces.shape == (2, 3, 3)
    assert np.all(result.forces == 0.0)
    assert result.metadata["forces_available"] is False
    print("test_energy_only: PASSED (energies returned, forces zero-filled, "
          "metadata flag set)")


if __name__ == "__main__":
    test_roundtrip()
    test_missing_pose()
    test_adsorbate_mismatch()
    test_energy_only()
    print("\nAll tests PASSED")
